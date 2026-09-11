"""
将 smc 共识 reads 按“在 reference 感兴趣区域是否有错误”分成 错误组 / 正确组 / (丢弃) 三组，
并据此把 adapter.bam 按 channel (ch tag) 切成 错误组 / 正确组 两个子 BAM。

流程 (对应 troubleshooting/error_channels.md):
  1. 用 gsmm2 将 smc 输出的 fastq 比对到 reference。
  2. 对每条 (primary) alignment，按 reference 的感兴趣 region 分类其 query：
        - error  : query 在 region 内有错误 (mismatch / 缺失即 deletion / insertion)
        - eq     : query 在 region 内无错误 (完全一致)
        - (丢弃) : query 的比对未覆盖 region (没有比对到该 region)
  3. 由 query 名字抽取 channel (query 名形如 `{run}/{channel}/{subread}`，channel = split("/")[1])，
     得到 错误组 channels、正确组 channels 两组；丢弃组不保留。
  4. 遍历 adapter.bam，按 record 的 `ch` tag 关联，把属于 错误组 / 正确组 channel 的
     record 分别写入 错误组 / 正确组 两个子 BAM。

注意：
  - region 使用 reference 坐标，0-based 半开区间 [start, end)，如 [1091, 1092) 即第 1091 个碱基。
  - gsmm2 输出的 CIGAR 只含 `= X I D S`（显式 eqx 操作，没有裸 `M`），因此可直接用 CIGAR 判断
    region 处是 mismatch (=error) 还是 match (=eq)，无需解析有损的 cs / md tag。
  - 每个 channel 对应一个 query (一条 primary alignment)，与需求中“把 query 分三组”一一对应。
  - adapter.bam 通常是整 run 的 subread 文件，可能很大 (数十 GB)，无 .bai，
    因此 step 4 是对整文件的一次流式扫描；step 1-3 很快，可先得到 channel 名单。

用法:
    python error_channels.py \
        --reference Group_0_Adaptor-barcode277-0.consensus.fasta \
        --smc-fastq Group_0_Adaptor-barcode277-0.fastq \
        --adapter-bam 20260831_240601Y0014_Run0005_called_demuxed_v4.bam \
        --region 1091:1092 \
        --outdir out

    # 只想要 channel 名单 (跳过对大 adapter.bam 的扫描): 去掉 --adapter-bam 即可。
    # 已有 gsmm2 比对结果时可用 --aligned-bam 跳过 step 1。
    # 可用 --phred-thr 30 在比对前按 read phreq 过滤 smc-fastq (低于阈值的 read 丢弃)。
"""

import argparse
import logging
import math
import os
import subprocess
import sys
from collections import defaultdict

import pysam
from tqdm import tqdm

from gseda.fastx_ana.fastx_basic_stat import phred33_to_rq

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)
log = logging.getLogger(__name__)

# pysam CIGAR op 码
_OP_I, _OP_D, _OP_M, _OP_EQ, _OP_X = 1, 2, 0, 7, 8


def parse_region(region_str: str):
    """把 '1091:1092' / '1091,1092' / '1091-1092' 解析为 (start, end)，0-based 半开。"""
    for sep in (":", ",", "-"):
        if sep in region_str:
            a, b = region_str.split(sep, 1)
            return int(a), int(b)
    raise ValueError(f"无法解析 region: {region_str!r}，请用 start:end (0-based 半开)")


def channel_of(qname: str):
    """从 smc query 名抽取 channel。形如 `{run}/{channel}/{subread}` -> split("/")[1]。"""
    parts = qname.split("/")
    if len(parts) < 2:
        return None
    try:
        return int(parts[1])
    except ValueError:
        return None


def classify_region(rec, start: int, end: int) -> str:
    """按 region [start, end) 对一条 primary alignment 分类: 'error' / 'eq' / 'noalign'。

    沿 CIGAR 从 rec.pos 起累加 reference 坐标：
      - 匹配列 (= / M) 覆盖 region  -> 'eq'  (其中 X 列覆盖 -> 'error')
      - 缺失 D 覆盖 region          -> 'error'
      - 插入 I 落在 region 边界上    -> 'error'
      - 比对未覆盖 region            -> 'noalign'
    """
    if rec.is_unmapped or rec.cigar is None:
        return "noalign"
    refpos = rec.pos  # 下一个 reference 碱基的 0-based 坐标
    for op, ln in rec.cigartuples:
        if op in (_OP_M, _OP_EQ, _OP_X):      # 消耗 reference
            lo, hi = refpos, refpos + ln - 1
            if hi >= start and lo < end:
                return "error" if op == _OP_X else "eq"
            refpos += ln
        elif op == _OP_D:                     # 缺失，消耗 reference
            lo, hi = refpos, refpos + ln - 1
            if hi >= start and lo < end:
                return "error"
            refpos += ln
        elif op == _OP_I:                     # 插入，不消耗 reference，坐标为 refpos
            if start <= refpos <= end:
                return "error"
        # S/H/N/P 不消耗 reference，忽略
    return "noalign"


def rq_to_phreq(rq: float) -> float:
    """rq (accuracy) -> 平均 read-level PHRED 分 (与 fastx_basic_stat 的 phreq 口径一致)。"""
    if rq >= 1 - 1e-10:
        rq = 1 - 1e-10
    if rq <= 0:
        return 0.0
    return -10.0 * math.log10(1.0 - rq)


def filter_fastq_by_phreq(smc_fastq: str, out_fastq: str, phreq_thr: float):
    """读取 smc fastq，按 read 的 phreq 过滤，phreq < phreq_thr 的 read 丢弃。

    返回 (total, kept, dropped)。phreq 口径与 fastx_basic_stat 一致：
    rq = 1 - mean(10^(q/-10))，phreq = -10*log10(1 - rq)。
    """
    total = kept = 0
    with pysam.FastxFile(smc_fastq) as reader, open(out_fastq, "w", encoding="utf-8") as out:
        for rec in tqdm(reader, desc=f"phreq filter (thr={phreq_thr})"):
            total += 1
            rq = phred33_to_rq(rec.quality) if rec.quality else 0.0
            if rq_to_phreq(rq) < phreq_thr:
                continue
            kept += 1
            out.write(f"@{rec.name}\n{rec.sequence}\n+\n{rec.quality}\n")
    return total, kept, total - kept


def gsmm2_align(smc_fastq: str, reference: str, prefix: str, threads: int):
    """gsmm2 将 fastq 比对到 reference，输出 `${prefix}.bam` (并建 .bai)。返回输出路径。"""
    # 注意 --threads / --preset 是 gsmm2 的顶层参数，需放在子命令 align 之前
    cmd = ["gsmm2", "--threads", str(threads), "align", "-q", smc_fastq,
           "--target", reference, "-p", prefix, "--noMar"]
    log.info("Step 1: gsmm2 比对  %s -> %s.bam", smc_fastq, prefix)
    log.info("$ %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    sys.stdout.write(proc.stdout or "")
    sys.stderr.write(proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError(f"gsmm2 align 失败 (rc={proc.returncode})")
    return prefix + ".bam"


def classify_channels(aligned_bam: str, start: int, end: int):
    """遍历 aligned.bam 的 primary alignment，返回 (error_set, eq_set, stats)。"""
    error_set, eq_set = set(), set()
    stats = defaultdict(int)
    seen_channels = set()
    with pysam.AlignmentFile(aligned_bam, "rb", check_sq=False) as af:
        for rec in af.fetch(until_eof=True):
            # 只用 primary alignment：每个 channel 一个 query，避免 secondary/supplementary 干扰
            if rec.is_secondary or rec.is_supplementary:
                continue
            stats["aligned"] += 1
            ch = channel_of(rec.qname)
            if ch is None:
                stats["no_channel"] += 1
                continue
            if ch in seen_channels:
                stats["dup_channel"] += 1
                continue
            seen_channels.add(ch)
            c = classify_region(rec, start, end)
            if c == "error":
                error_set.add(ch)
            elif c == "eq":
                eq_set.add(ch)
            else:
                stats["dropped"] += 1
    return error_set, eq_set, dict(stats)


def write_channel_lists(outdir: str, error_set, eq_set, dropped: int):
    paths = {}
    for name, s in (("error", error_set), ("eq", eq_set)):
        path = os.path.join(outdir, f"{name}_channels.txt")
        with open(path, "w", encoding="utf-8") as f:
            for ch in sorted(s):
                f.write(f"{ch}\n")
        paths[name] = path
    summary = os.path.join(outdir, "channel_groups_summary.txt")
    with open(summary, "w", encoding="utf-8") as f:
        f.write(f"error_channels\t{len(error_set)}\n")
        f.write(f"eq_channels\t{len(eq_set)}\n")
        f.write(f"dropped_channels\t{dropped}\n")
        f.write(f"both_error_and_eq\t{len(error_set & eq_set)}\n")
    paths["summary"] = summary
    return paths


def split_adapter_bam(adapter_bam: str, out_error_bam: str, out_eq_bam: str,
                      error_set, eq_set, threads: int = 40):
    """一次流式扫描 adapter.bam，按 ch tag 把 record 分流到 错误组 / 正确组 BAM。"""
    log.info("Step 4: 扫描 adapter.bam  按 ch 切分 (error=%d, eq=%d channels)",
             len(error_set), len(eq_set))
    counts = defaultdict(int)
    with pysam.AlignmentFile(adapter_bam, "rb", threads=threads, check_sq=False) as in_bam:
        header = in_bam.header
        with pysam.AlignmentFile(out_error_bam, "wb", header=header,
                                 threads=threads, check_sq=False) as eout, \
             pysam.AlignmentFile(out_eq_bam, "wb", header=header,
                                 threads=threads, check_sq=False) as qout:
            for rec in tqdm(in_bam.fetch(until_eof=True), desc="splitting adapter.bam"):
                counts["total"] += 1
                if not rec.has_tag("ch"):
                    counts["no_ch"] += 1
                    continue
                ch = int(rec.get_tag("ch"))
                if ch in error_set:
                    eout.write(rec)
                    counts["error"] += 1
                elif ch in eq_set:
                    qout.write(rec)
                    counts["eq"] += 1
                else:
                    counts["dropped"] += 1
    return dict(counts)


def main_cli(argv=None):
    parser = argparse.ArgumentParser(
        prog="error_channels",
        description="按 reference 感兴趣区域是否有错误，把 smc reads 分成 错误组/正确组，"
                    "并据此切分 adapter.bam。",
    )
    parser.add_argument("--reference", required=True, help="reference FASTA (consensus)")
    parser.add_argument("--smc-fastq", required=True,
                        help="smc 输出的 fastq (除非用 --aligned-bam)")
    parser.add_argument("--phred-thr", type=float, default=20., dest="phred_thr",
                        help="read phreq 阈值；phreq 低于该值的 read 在读取 smc-fastq 时直接丢弃 "
                             "(仅在走 gsmm2 比对、即未提供 --aligned-bam 时生效)")
    parser.add_argument("--region", required=True,
                        help="感兴趣区域，0-based 半开 start:end，如 1091:1092")
    parser.add_argument("--adapter-bam", default=None,
                        help="要切分的 adapter.bam；不给则只产出 channel 名单")
    parser.add_argument("--aligned-bam", default=None,
                        help="已存在的 gsmm2 比对 BAM，给了则跳过 step 1")
    parser.add_argument("--outdir", default="error_channels_out", help="输出目录")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1,
                        help="并行线程数 (默认 CPU 核数)")
    args = parser.parse_args(argv)

    start, end = parse_region(args.region)
    os.makedirs(args.outdir, exist_ok=True)
    base = os.path.join(args.outdir, "error_channels")

    # Step 1: 比对 (可跳过)
    if args.aligned_bam:
        aligned_bam = args.aligned_bam
        log.info("使用已有 aligned bam: %s", args.aligned_bam)
        if args.phred_thr is not None:
            log.warning("--phred-thr 已给定，但使用了 --aligned-bam，phreq 过滤被忽略")
    else:
        query_fastq = args.smc_fastq
        if args.phred_thr is not None:
            query_fastq = base + ".phred_filtered.fastq"
            log.info("Step 0: 按 phreq 过滤 smc fastq (phred-thr=%s)", args.phred_thr)
            total, kept, dropped = filter_fastq_by_phreq(
                args.smc_fastq, query_fastq, args.phred_thr)
            log.info("phreq 过滤: total=%d kept=%d dropped=%d -> %s",
                     total, kept, dropped, query_fastq)
            if kept == 0:
                raise RuntimeError("phreq 过滤后没有剩余 read，无法继续")
        aligned_bam = gsmm2_align(query_fastq, args.reference, base + ".aligned", args.threads)

    # Step 2-3: 分类 -> channel 名单
    log.info("Step 2: 按 region [%d, %d) 分类 primary alignments", start, end)
    error_set, eq_set, stats = classify_channels(aligned_bam, start, end)
    log.info("分类结果: aligned=%s dropped=%s (error=%d eq=%d)",
             stats.get("aligned"), stats.get("dropped"), len(error_set), len(eq_set))
    paths = write_channel_lists(args.outdir, error_set, eq_set, stats.get("dropped", 0))
    log.info("channel 名单已写出: %s", paths)

    # Step 4: 切分 adapter.bam (可选)
    if args.adapter_bam:
        sc = split_adapter_bam(
            args.adapter_bam,
            base + ".adapter_error.bam",
            base + ".adapter_eq.bam",
            error_set, eq_set, threads=args.threads,
        )
        log.info("adapter 切分完成: %s", sc)

    log.info("完成。输出目录: %s", args.outdir)
    return paths


if __name__ == "__main__":
    main_cli()
