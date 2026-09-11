#!/usr/bin/env python3
"""分析比对到 reference 某个 region 上的 query base 的 dw / ar 差异分布。

对照 spec: troubleshooting/ref_region_ar_dw_analysis.md

流程:
  1. 两个输入 bam (对照组 / 实验组), 各自用 gsmm2 比对到 reference。
  2. 抽取比对到指定 region 上的 query base, 每个 base 取其 dw (dwell time) 与
     ar (arrival time) 值。dw/ar 是 BAM 中按 query 坐标存放的 per-base 数组
     (query 反向互补时 dw/ar 已随之反向, 见 spec 注意 #2)。
  3. 画差异分布直方图 (用 frequency 而非 count, 因为两个 bam 的 base 数差异大):
       - dw:  4 个子图 (A/C/G/T), 每图叠加 对照 + 实验 两组分布
       - ar: 16 个子图 (A->A, A->C, ... ), 每图叠加 对照 + 实验 两组分布

坐标系约定 (spec 注意 #1): 比对以 reference 方向为正, 分析以 query 方向为正。
  pysam 的 query_sequence 与 get_aligned_pairs 的 qpos 本身就处于 query 正向
  坐标, dw/ar 也按其存放, 因此全程在 query 坐标取 dw[qpos] / ar[qpos] 即可。

dw 语义: 每个 base 的 dw 表示该 base 测序时信号持续了多长时间。
ar 语义: 每个 base 的 ar 表示 前一个 base 信号结束 到 当前 base 信号开始 的时长,
  即 ar[qpos] 描述 base[qpos-1] -> base[qpos] 这一转移; 第一个 base (qpos=0)
  无前一 base, 其 ar 不计入转移分布。
"""

import matplotlib.pyplot as plt
import argparse
import logging
import os
import pathlib
import subprocess
from multiprocessing import cpu_count

import numpy as np
import pysam
import matplotlib

matplotlib.use("Agg")

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)

BASES = "ACGT"
AR_TRANSITIONS = [f"{p}->{c}" for p in BASES for c in BASES]
# 颜色: 对照组蓝色, 实验组橙色
COLOR_CONTROL = "#1f77b4"
COLOR_TREATMENT = "#ff7f0e"
GROUP_LABELS = ("control", "treatment")


def align_bam(bam_file: str, ref_fasta: str, prefix: str, threads: int,
              preset: str) -> str:
    """用 gsmm2 将 bam 比对到 reference, 返回坐标排序并建好索引的 aligned bam 路径。"""
    for fpath in [prefix + ext for ext in (".bam", ".bam.bai")]:
        if os.path.exists(fpath):
            os.remove(fpath)

    cmd = (f"gsmm2 --threads {threads} align -q {bam_file} -t {ref_fasta} "
           f"-p {prefix} --pt_tags dw,ar --noMar --query-forward")
    if preset:
        cmd += f" --preset {preset}"
    logging.info("cmd: %s", cmd)
    subprocess.check_call(cmd, shell=True)

    aligned_bam = f"{prefix}.bam"
    return aligned_bam


def _get_base_array(read, tag):
    """取某个 read 的 per-base dw/ar 数组; 缺失返回 None。"""
    if not read.has_tag(tag):
        return None
    arr = read.get_tag(tag)
    if isinstance(arr, (bytes, str)):
        arr = list(arr)
    return arr


def collect_core(bam, ref_name, ref_start, ref_end, dw_by_base, ar_by_trans):
    n_dw = 0
    n_ar = 0

    for read in bam.fetch(ref_name, ref_start, ref_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        seq = read.query_sequence
        dw = _get_base_array(read, "dw")
        ar = _get_base_array(read, "ar")
        if dw is None or ar is None:
            continue
        seq_len = len(seq)
        if len(dw) != seq_len or len(ar) != seq_len:
            logging.warning(
                "read %s: dw/ar length (%d/%d) != query length %d, skipped",
                read.qname, len(dw), len(ar), seq_len)
            continue

        prev_base = None
        for qpos, rpos in read.get_aligned_pairs():
            if rpos is None or qpos is None:
                continue
            if rpos < ref_start:
                continue
            if rpos >= ref_end:
                break
            if qpos is None:
                continue
            base = seq[qpos]
            if base not in dw_by_base:
                continue
            dw_by_base[base].append(float(dw[qpos]))
            n_dw += 1

            if prev_base is not None:
                trans = f"{prev_base}->{base}"
                if trans in ar_by_trans:
                    ar_by_trans[trans].append(float(ar[qpos]))
                    n_ar += 1
            prev_base = base

    return n_dw, n_ar


def collect_region_ar_dw(aligned_bam: str, ref_name: str,
                         ref_start: int, ref_end: int, fwd_only=False, rev_only=False):
    """遍历比对到 region [ref_start, ref_end) 上的 query base, 收集 (base,dw) 与 (prev->cur,ar)。

    Returns:
        dw_by_base:  {base: [dw_value, ...]}  (按 query 正序 base 分类)
        ar_by_trans: {trans: [ar_value, ...]} (trans 形如 'A->C', 前 base -> 当前 base)
        n_dw, n_ar:  实际收集到的 dw / ar 个数
    """
    dw_by_base = {b: [] for b in BASES}
    ar_by_trans = {t: [] for t in AR_TRANSITIONS}
    n_dw = 0
    n_ar = 0

    with pysam.AlignmentFile(aligned_bam, "rb", threads=cpu_count(),
                             check_sq=False) as bam:
        fwd_ref_name = f"{ref_name}___fwd"
        rev_ref_name = f"{ref_name}___rev"
        ref_len = bam.get_reference_length(fwd_ref_name)

        rev_ref_start = ref_len - ref_end
        rev_ref_end = ref_len - ref_start
        fwd_dw, fwd_ar = 0, 0
        if not rev_only:
            fwd_dw, fwd_ar = collect_core(bam=bam, ref_name=fwd_ref_name, ref_start=ref_start,
                                          ref_end=ref_end, dw_by_base=dw_by_base, ar_by_trans=ar_by_trans)

        rev_dw, rev_ar = 0, 0
        if not fwd_only:
            rev_dw, rev_ar = collect_core(bam=bam, ref_name=rev_ref_name, ref_start=rev_ref_start,
                                          ref_end=rev_ref_end, dw_by_base=dw_by_base, ar_by_trans=ar_by_trans)

        n_dw = fwd_dw + rev_dw
        n_ar = fwd_ar + rev_ar

    return dw_by_base, ar_by_trans, n_dw, n_ar


def _shared_bins(values, bins: int):
    """取覆盖两组数据 (合并后) 的分箱边界, 让对照/实验可比。"""
    combined = np.concatenate(values) if any(values) else np.array([0.0, 1.0])
    return np.histogram_bin_edges(combined, bins=bins)


def _plot_subplot(ax, groups, bins, title):
    for (label, vals), color in zip(groups, (COLOR_CONTROL, COLOR_TREATMENT)):
        if vals:
            ax.hist(vals, bins=bins, density=True, alpha=0.55, color=color,
                    label=label)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("value", fontsize=9)
    ax.set_ylabel("frequency", fontsize=9)
    ax.tick_params(labelsize=8)
    if any(vals for _, vals in groups):
        ax.legend(fontsize=7)


def plot_dw(dw_ctrl, dw_trt, out_png: str, bins: int):
    """2x2: A/C/G/T 的 dw 分布, 每图对照 vs 实验。"""
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for ax, base in zip(axes.ravel(), BASES):
        groups = [(GROUP_LABELS[0], dw_ctrl[base]),
                  (GROUP_LABELS[1], dw_trt[base])]
        _plot_subplot(ax, groups, _shared_bins([dw_ctrl[base], dw_trt[base]], bins),
                      f"dw on base {base}")
    fig.suptitle("Dwell time (dw) distribution per base in reference region",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_ar(ar_ctrl, ar_trt, out_png: str, bins: int):
    """4x4: 16 个 base 转移 (prev->cur) 的 ar 分布, 每图对照 vs 实验。"""
    fig, axes = plt.subplots(4, 4, figsize=(16, 14))
    for ax, trans in zip(axes.ravel(), AR_TRANSITIONS):
        groups = [(GROUP_LABELS[0], ar_ctrl[trans]),
                  (GROUP_LABELS[1], ar_trt[trans])]
        _plot_subplot(ax, groups, _shared_bins([ar_ctrl[trans], ar_trt[trans]], bins),
                      f"ar {trans}")
    fig.suptitle("Arrival time (ar) distribution per base transition in "
                 "reference region", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def analyze(control_bam: str, treatment_bam: str, ref_fasta: str, ref_name: str,
            ref_start: int, ref_end: int, outdir: str, threads: int,
            bins: int, preset: str = "", fwd_only=False, rev_only=False):
    threads = cpu_count() if threads is None else threads
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ctrl_prefix = os.path.join(outdir, "control.align")
    trt_prefix = os.path.join(outdir, "treatment.align")
    logging.info("aligning control bam ...")
    ctrl_aln = align_bam(control_bam, ref_fasta, ctrl_prefix, threads, preset)
    logging.info("aligning treatment bam ...")
    trt_aln = align_bam(treatment_bam, ref_fasta, trt_prefix, threads, preset)

    logging.info("collecting region %s:[%d,%d) dw/ar ...", ref_name, ref_start,
                 ref_end)
    dw_ctrl, ar_ctrl, n_dw_ctrl, n_ar_ctrl = collect_region_ar_dw(
        ctrl_aln, ref_name, ref_start, ref_end, fwd_only=fwd_only, rev_only=rev_only)
    dw_trt, ar_trt, n_dw_trt, n_ar_trt = collect_region_ar_dw(
        trt_aln, ref_name, ref_start, ref_end, fwd_only=fwd_only, rev_only=rev_only)

    logging.info("control   : %d dw, %d ar values", n_dw_ctrl, n_ar_ctrl)
    logging.info("treatment : %d dw, %d ar values", n_dw_trt, n_ar_trt)
    if n_dw_ctrl == 0 or n_dw_trt == 0:
        logging.warning(
            "one of the groups has no dw values in the region; check "
            "--ref-name / --ref-start / --ref-end and that the bam carries "
            "dw/ar tags")

    dw_png = os.path.join(outdir, "dw_by_base.png")
    ar_png = os.path.join(outdir, "ar_by_transition.png")
    plot_dw(dw_ctrl, dw_trt, dw_png, bins)
    plot_ar(ar_ctrl, ar_trt, ar_png, bins)
    logging.info("wrote %s", dw_png)
    logging.info("wrote %s", ar_png)
    return dw_png, ar_png


def main_cli():
    parser = argparse.ArgumentParser(
        description="Compare dw/ar distributions of query bases mapping into a "
                    "reference region, between a control and a treatment bam.")
    parser.add_argument("--control-bam", required=True,
                        help="对照组 bam (未比对, 将被 gsmm2 比对到 reference)")
    parser.add_argument("--treatment-bam", required=True,
                        help="实验组 bam (未比对, 将被 gsmm2 比对到 reference)")
    parser.add_argument("--ref", required=True, help="reference fasta")
    parser.add_argument("--ref-name", default=None,
                        help="region 所在的 contig 名 (reference 中的 rname)")
    parser.add_argument("--ref-start", type=int, required=True,
                        help="region 起点, 0-based (含)")
    parser.add_argument("--ref-end", type=int, required=True,
                        help="region 终点, 0-based 半开区间 [start, end)")

    parser.add_argument("--fwd-only", action="store_true")
    parser.add_argument("--rev-only", action="store_true")

    parser.add_argument("--outdir", default=None, help="输出目录")
    parser.add_argument("--threads", type=int, default=None)

    parser.add_argument("--bins", type=int, default=60, help="每个子图的分箱数")
    parser.add_argument("--preset", type=str, default="",
                        help="gsmm2 比对 preset, 如 map-ont (默认不指定)")

    args = parser.parse_args()

    ref_name = args.ref_name

    if args.ref_start < 0 or args.ref_end <= args.ref_start:
        raise ValueError("require 0 <= --ref-start < --ref-end")

    ref_fasta = pysam.FastaFile(args.ref)

    if ref_name is None:
        ref_name = ref_fasta.references[0]

    ref_seq = ref_fasta.fetch(ref_name)

    expanded_start = max(0, args.ref_start - 10)
    expanded_end = min(args.ref_end + 10, len(ref_seq))

    print(
        f"interested_region. {ref_seq[expanded_start:args.ref_start]}[{ref_seq[args.ref_start: args.ref_end]}]{ref_seq[args.ref_end:expanded_end]}")

    if args.outdir is None:
        args.outdir = os.path.join(
            os.path.dirname(args.ref) or ".",
            f"ref_region_ar_dw_{pathlib.Path(args.ref).stem}")

    dw_image, ar_image = analyze(args.control_bam, args.treatment_bam, args.ref, ref_name,
                                 args.ref_start, args.ref_end, args.outdir, args.threads, args.bins,
                                 args.preset, args.fwd_only, args.rev_only)
    images = [dw_image, ar_image]
    gsda_tmp_path = pathlib.Path("/root/projects/gsda/tmp-data-dir")
    import shutil
    if gsda_tmp_path.exists():
        for image_path in images:
            shutil.copy(image_path, gsda_tmp_path /
                        pathlib.Path(image_path).name)


if __name__ == "__main__":
    main_cli()
