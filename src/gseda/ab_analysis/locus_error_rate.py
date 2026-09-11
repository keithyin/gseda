#!/usr/bin/env python3
"""reference 不同位点准确率的 A/B 分析。

对照 spec: ab_analysis/locus_error_rate_design.md

流程:
  1. query 分实验组 / 对照组两组, 各自用 gsmm2 比对到 reference。比对时施加
     read-accuracy 阈值 (rq-range), 阈值可以是 1 个 (两组共用) 或 2 个 (各自一个)。
  2. 每组的 aligned bam 调 gsetl aligned-bam, 产出
     fact_aligned_bam_ref_locus_info.csv (每个 reference 位点的 eq/diff/ins/del/depth)。
  3. 在两张位点表上各自增加一列 “位点准确率” locus_accuracy = eq/(eq+diff+ins+del)。
  4. 把实验组与对照组两张表按 (refname, pos) join 成一张表输出, 供后续判断。

query 格式 bam (unmapped) / fastq 均可, gsmm2 align 的 -q 直接接受。
rq-range 只在 query 为带 rq 字段的 bam 时真正生效; fastq 输入时该参数被 gsmm2 忽略。
"""

import argparse
import logging
import os
import subprocess
import sys
from multiprocessing import cpu_count

import polars as pl

from gseda.fact_table_ana.polars_init import polars_env_init

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)
log = logging.getLogger(__name__)


def _parse_thr_list(s: str, name: str) -> list:
    """把 '0.99' 或 '0.99,0.97' 解析成 list[float] (1~2 个)。"""
    vals = [float(x) for x in s.replace(" ", "").split(",") if x != ""]
    if len(vals) not in (1, 2):
        raise ValueError(f"{name} 只能填 1 或 2 个值, 得到: {vals}")
    for v in vals:
        if v < 0:
            raise ValueError(f"{name} 每个值需 >= 0, 得到: {v}")
    return vals


def parse_rq_thr(s: str) -> list:
    """read-accuracy 阈值, 1~2 个值, 每个需在 [0,1]。"""
    vals = _parse_thr_list(s, "rq-thr")
    for v in vals:
        if v > 1.0:
            raise ValueError(f"rq-thr 每个值需在 [0,1], 得到: {v}")
    return vals


def parse_np_thr(s: str) -> list:
    """number-of-passes 阈值, 1~2 个整数值 (无上限校验)。"""
    vals = [int(float(x)) for x in s.replace(" ", "").split(",") if x != ""]
    if len(vals) not in (1, 2):
        raise ValueError(f"np-thr 只能填 1 或 2 个值, 得到: {vals}")
    for v in vals:
        if v < 0:
            raise ValueError(f"np-thr 每个值需 >= 0, 得到: {v}")
    return vals


def _resolve_range(thrs: list, group: str, upper: str) -> str:
    """1 个阈值 -> 两组共用; 2 个 -> [0]=对照组, [1]=实验组。返回 '{thr}:{upper}' 串。"""
    thr = thrs[0] if len(thrs) == 1 else thrs[0 if group == "control" else 1]
    return f"{thr}:{upper}"


def resolve_rq_range(rq_thrs: list, group: str) -> str:
    """rq-range, 上界 1.1 (accuracy 上限, 永不因上限丢 read)。"""
    return _resolve_range(rq_thrs, group, "1.1")


def resolve_np_range(np_thrs: list, group: str) -> str:
    """np-range, 上界取极大值 (passes 上限, 永不因上限丢 read)。"""
    return _resolve_range(np_thrs, group, "100000000")


def gsmm2_align(query: str, ref: str, prefix: str, rq_range: str,
                np_range, threads: int) -> str:
    """gsmm2 将 query (bam/unmapped 或 fastq) 比对到 ref, 返回 aligned bam 路径。

    rq_range / np_range 是 '{thr}:{upper}' 串, 仅在 query 为带对应字段的 bam 时生效;
    np_range 传 None 则不加 --np-range。
    """
    for fpath in (prefix + ".bam", prefix + ".bam.bai"):
        if os.path.exists(fpath):
            os.remove(fpath)
    cmd = ["gsmm2", "--threads", str(threads), "align", "-q", query,
           "-t", ref, "-p", prefix, "--noMar", "--rq-range", rq_range]
    if np_range is not None:
        cmd += ["--np-range", np_range]
    log.info("gsmm2 align (rq-range=%s np-range=%s): %s -> %s.bam",
             rq_range, np_range, query, prefix)
    log.info("$ %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    sys.stdout.write(proc.stdout or "")
    sys.stderr.write(proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError(f"gsmm2 align 失败 (rc={proc.returncode})")
    return prefix + ".bam"


def gsetl_locus(aligned_bam: str, ref: str, outdir: str) -> str:
    """gsetl aligned-bam, 只出 fact_ref_locus_info 表, 返回该 csv 路径。"""
    cmd = ["gsetl", "--outdir", outdir, "aligned-bam", "--bam", aligned_bam,
           "--ref-file", ref,
           "--factRecordStat", "0", "--factRefLocusInfo", "1",
           "--factBamBasic", "0", "--factErrorQueryLocusInfo", "0",
           "--factBaseQStat", "0", "--factPolyInfo", "0"]
    log.info("gsetl aligned-bam: %s -> %s", aligned_bam, outdir)
    log.info("$ %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    sys.stdout.write(proc.stdout or "")
    sys.stderr.write(proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError(f"gsetl aligned-bam 失败 (rc={proc.returncode})")
    return os.path.join(outdir, "fact_aligned_bam_ref_locus_info.csv")


def add_locus_accuracy(df: pl.DataFrame) -> pl.DataFrame:
    """在位点表上增加位点准确率列: eq/(eq+diff+ins+del)。

    小数指标统一保留 6 位小数。
    """
    total = pl.col("eq") + pl.col("diff") + pl.col("ins") + pl.col("del")
    return df.with_columns([
        (pl.col("eq") / total).round(6).alias("locus_accuracy"),
        (pl.col("eq") / pl.col("depth")).round(6).alias("locus_accuracy_by_depth"),
    ])


def process_group(query: str, ref: str, rq_range: str, np_range, threads: int,
                  group: str, outdir: str) -> pl.DataFrame:
    """对一个组跑完 比对 -> gsetl -> 加准确率列, 返回带 group 标记的位点 df。"""
    align_prefix = os.path.join(outdir, f"{group}.aligned")
    aligned_bam = gsmm2_align(query, ref, align_prefix, rq_range, np_range, threads)

    gsetl_dir = os.path.join(outdir, f"{group}-gsetl")
    locus_csv = gsetl_locus(aligned_bam, ref, gsetl_dir)

    df = pl.read_csv(locus_csv, separator="\t")
    df = add_locus_accuracy(df)
    # 只保留判断需要的列, 打上组标记
    keep = ["refname", "pos", "eq", "diff", "ins", "del", "depth",
            "aroundBases",
            "locus_accuracy", "locus_accuracy_by_depth"]
    df = df.select(keep).with_columns(pl.lit(group).alias("group"))
    log.info("group=%s 位点数=%d", group, df.height)
    return df


def join_groups(ctrl: pl.DataFrame, exp: pl.DataFrame) -> pl.DataFrame:
    """按 (refname, pos) join 对照组与实验组, 组内列加后缀区分。"""
    c = ctrl.drop("group")
    e = exp.drop("group")
    merged = c.join(e, on=["refname", "pos"], how="inner", suffix="_exp")
    # c 侧列补 _ctrl 后缀 (除 join key)
    c_suf = {n: f"{n}_ctrl" for n in c.columns if n not in ("refname", "pos")}
    merged = merged.rename(c_suf)
    return merged


def main_cli(argv=None):
    polars_env_init()
    parser = argparse.ArgumentParser(
        prog="locus_error_rate",
        description="reference 不同位点准确率的 A/B 分析: 两组 query 分别 gsmm2 比对 + "
                    "gsetl 出位点表, 加位点准确率后 join 成一张表。")
    parser.add_argument("--control", required=True,
                        help="对照组 query (unmapped bam 或 fastq)")
    parser.add_argument("--exp", required=True,
                        help="实验组 query (unmapped bam 或 fastq)")
    parser.add_argument("--ref", required=True, help="reference fasta")
    parser.add_argument("--rq-thr", required=True,
                        help="read-accuracy 阈值, 1 或 2 个值 (逗号分隔)。"
                             "1 个 -> 两组共用; 2 个 -> [对照组, 实验组]。"
                             "例: 0.99 或 0.99,0.97")
    parser.add_argument("--np-thr", default="5",
                        help="number-of-passes 阈值 (gsmm2 --np-range 下界), 1 或 2 个值 "
                             "(逗号分隔), 同 rq-thr 的分组规则。不给则默认 5。"
                             "例: 5 或 5,3")
    parser.add_argument("--outdir", default="locus_error_rate_out", help="输出目录")
    parser.add_argument("--threads", type=int, default=None,
                        help="gsmm2 线程数 (默认 CPU 核数)")
    args = parser.parse_args(argv)

    threads = cpu_count() if args.threads is None else args.threads
    rq_thrs = parse_rq_thr(args.rq_thr)
    np_thrs = parse_np_thr(args.np_thr) if args.np_thr else None
    os.makedirs(args.outdir, exist_ok=True)

    rq_ctrl = resolve_rq_range(rq_thrs, "control")
    rq_exp = resolve_rq_range(rq_thrs, "experiment")
    log.info("rq-range: control=%s experiment=%s", rq_ctrl, rq_exp)
    if np_thrs is not None:
        np_ctrl = resolve_np_range(np_thrs, "control")
        np_exp = resolve_np_range(np_thrs, "experiment")
        log.info("np-range: control=%s experiment=%s", np_ctrl, np_exp)
    else:
        np_ctrl = np_exp = None

    ctrl_df = process_group(args.control, args.ref, rq_ctrl, np_ctrl, threads,
                            "control", args.outdir)
    exp_df = process_group(args.exp, args.ref, rq_exp, np_exp, threads,
                           "experiment", args.outdir)

    joined = join_groups(ctrl_df, exp_df)

    # 每张单组表也落一份, 方便单独看
    ctrl_path = os.path.join(args.outdir, "control_locus_accuracy.csv")
    exp_path = os.path.join(args.outdir, "experiment_locus_accuracy.csv")
    ctrl_df.write_csv(ctrl_path, separator="\t")
    exp_df.write_csv(exp_path, separator="\t")

    joined_path = os.path.join(args.outdir, "locus_accuracy_joined.csv")
    joined.write_csv(joined_path, separator="\t")

    log.info("完成。输出目录: %s", args.outdir)
    log.info("  control   : %s", ctrl_path)
    log.info("  experiment: %s", exp_path)
    log.info("  joined    : %s  (%d loci)", joined_path, joined.height)
    return joined_path


if __name__ == "__main__":
    main_cli()
