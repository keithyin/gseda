#!/usr/bin/env python3
"""
STR / SSR / homopolymer 区域检测

流程:
1. 粗筛: 调用 trfrs --minscore <N> --str 对 fasta 做 tandem repeat 搜索
2. 粗筛结果后处理:
   - 去除 overlap (同一序列上重叠的区域, 按 score 保留最优者)
   - motif 归一化: 如 GCGC -> GC (还原成最小重复单元)
3. 精筛, 提供两种模式:
   - --str 模式: 直接使用粗筛(后处理)结果
   - --ssr 模式: 整段干净时直接保留; 否则在区间内提取所有足够长的干净 unit 重复子区间

用法:
    python str_ssr_detect.py <fa|dir> [more...] [--mode str|ssr] [--minscore 9] [--out result.tsv]
"""

import argparse
import csv
import pathlib
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Dict, List

FASTA_SUFFIXES = {".fa", ".fasta", ".fna"}


# ---------------------------------------------------------------- fasta


def collect_fasta_files(paths: List[str]) -> List[pathlib.Path]:
    """收集输入: 支持单个 fasta 文件或目录(递归)"""
    files: List[pathlib.Path] = []
    for p in paths:
        p = pathlib.Path(p)
        if p.is_dir():
            for suffix in FASTA_SUFFIXES:
                files.extend(sorted(p.rglob(f"*{suffix}")))
        elif p.is_file():
            files.append(p)
        else:
            raise FileNotFoundError(f"input not found: {p}")
    if not files:
        raise FileNotFoundError("no fasta file found in inputs")
    return files


def parse_fasta(path: pathlib.Path) -> Dict[str, str]:
    """返回 {序列名: 序列}"""
    records: Dict[str, str] = {}
    name = None
    chunks: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


# ---------------------------------------------------------------- trfrs 粗筛


@dataclass
class Repeat:
    seq_name: str
    start: int          # 1-based inclusive
    end: int            # 1-based inclusive
    period: int         # motif 长度
    copy_number: float
    score: int
    percent_matches: float
    motif: str          # trfrs 原始 motif
    repeat_seq: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def run_trfrs(fasta_files: List[pathlib.Path], minscore: int, workdir: pathlib.Path) -> pathlib.Path:
    """合并输入后运行 trfrs 粗筛, 返回 csv 结果文件路径"""
    combined = workdir / "combined_input.fa"
    with open(combined, "w") as out:
        for fa in fasta_files:
            for name, seq in parse_fasta(fa).items():
                out.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + "\n")

    out_dir = workdir / "trfrs_out"
    cmd = ["trfrs", "--minscore", str(minscore), "--str", "--csv", "-o", str(out_dir), str(combined)]
    subprocess.check_call(cmd)

    csv_files = sorted(out_dir.glob("*.csv"))
    if not csv_files:
        raise RuntimeError(f"trfrs produced no csv output in {out_dir}")
    return csv_files[0]


def parse_trfrs_csv(csv_path: pathlib.Path) -> List[Repeat]:
    repeats: List[Repeat] = []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            repeats.append(Repeat(
                seq_name=row["sequence_name"],
                start=int(row["start"]),
                end=int(row["end"]),
                period=int(row["period"]),
                copy_number=float(row["copy_number"]),
                score=int(row["score"]),
                percent_matches=float(row["percent_matches"]),
                motif=row["repeat_unit"].upper(),
                repeat_seq=row["repeat_sequence"].upper(),
            ))
    return repeats


# ---------------------------------------------------------------- 后处理


def normalize_motif(motif: str) -> str:
    """归一化 motif: 还原成最小重复单元, 例如 GCGC -> GC"""
    n = len(motif)
    for d in range(1, n // 2 + 1):
        if n % d == 0 and motif == motif[:d] * (n // d):
            return motif[:d]
    return motif


def cyclic_rotations(unit: str) -> List[str]:
    """返回 unit 的所有不同循环移位, 例如 GA -> [GA, AG]。

    TRF 的 motif 只确定到循环移位: GA* 与 AG* 描述同一个串联重复, 只是相位
    相差一个碱基。精筛必须同时考虑所有相位, 否则会按报告里的相位截断拷贝数。
    """
    if not unit:
        return []
    result: List[str] = []
    for k in range(len(unit)):
        rot = unit[k:] + unit[:k]
        if rot not in result:
            result.append(rot)
    return result


def resolve_overlaps(repeats: List[Repeat]) -> List[Repeat]:
    """去除 overlap: 按 score 从高到低贪心选择互不重叠的区域"""
    picked: List[Repeat] = []
    for r in sorted(repeats, key=lambda x: (-x.score, -x.length)):
        if any(not (r.end < p.start or r.start > p.end) for p in picked if p.seq_name == r.seq_name):
            continue
        picked.append(r)
    return sorted(picked, key=lambda x: (x.seq_name, x.start, x.end))


# ---------------------------------------------------------------- 精筛


def is_clean_ssr(repeat: Repeat, ref_seq: str) -> bool:
    """区域序列必须是某个 unit 循环移位的完全干净重复 (相位无关)"""
    region = ref_seq[repeat.start - 1:repeat.end]
    unit = normalize_motif(repeat.motif)
    if not unit or len(region) % len(unit) != 0:
        return False
    copies = len(region) // len(unit)
    return any(region == rot * copies for rot in cyclic_rotations(unit))


def find_clean_subregions(repeat: Repeat, ref_seq: str, min_copies: int = 3) -> List[Repeat]:
    """从不干净的 call 中提取所有 >= min_copies 个连续 unit 的干净子区间。

    TRF 的 motif 只确定到循环移位 (GA* 与 AG* 是同一个重复, 相位差一个碱基),
    这里对 motif 的每一种循环移位分别扫描。同一段物理重复在不同相位下会产生
    相互重叠的候选 (如 GA*8 与 AG*7), 最终只保留最长者, 避免少算一个 copy。
    """
    unit = normalize_motif(repeat.motif)
    if not unit:
        return []

    step = len(unit)
    region = ref_seq[repeat.start - 1:repeat.end]

    candidates: List[Repeat] = []
    for rot in cyclic_rotations(unit):
        i = 0
        while i + step <= len(region):
            if region[i:i + step] != rot:
                i += 1
                continue
            j = i
            while j + step <= len(region) and region[j:j + step] == rot:
                j += step
            copies = (j - i) // step
            if copies >= min_copies:
                candidates.append(Repeat(
                    seq_name=repeat.seq_name,
                    start=repeat.start + i,
                    end=repeat.start + j - 1,
                    period=step,
                    copy_number=float(copies),
                    score=repeat.score,
                    percent_matches=100.0,
                    motif=rot,
                    repeat_seq=rot * copies,
                ))
            i = j

    # 相位不同的候选可能相交但不完全相同 (如 GA[1,14] 与 AG[2,15]), 它们描述的是
    # 同一段物理重复而不是两个位点, 不能都保留。按拷贝数取最长者; 长度相同时取
    # 起点更靠前者, 不偏袒 TRF 报告的那个相位。
    best: List[Repeat] = []
    for cand in sorted(candidates, key=lambda r: (-r.length, r.start)):
        if any(not (cand.end < b.start or cand.start > b.end) for b in best):
            continue
        best.append(cand)
    return sorted(best, key=lambda r: (r.seq_name, r.start, r.end))


def fine_filter(repeats: List[Repeat], refs: Dict[str, str], mode: str,
                min_copies: int = 3) -> List[Repeat]:
    if mode != "ssr":
        return repeats
    result: List[Repeat] = []
    for repeat in repeats:
        result.extend(find_clean_subregions(repeat, refs[repeat.seq_name], min_copies))
    return sorted(result, key=lambda r: (r.seq_name, r.start, r.end))


# ---------------------------------------------------------------- 输出


def write_result(repeats: List[Repeat], out) -> None:
    header = ["seq_name", "start", "end", "length", "period", "motif",
              "copy_number", "score", "percent_matches", "repeat_seq"]
    out.write("\t".join(header) + "\n")
    for r in repeats:
        out.write("\t".join([
            r.seq_name, str(r.start), str(r.end), str(r.length),
            str(r.period), normalize_motif(r.motif),
            f"{r.length / r.period:.2f}", str(r.score),
            f"{r.percent_matches:.1f}", r.repeat_seq,
        ]) + "\n")


def write_patterns(repeats: List[Repeat], out) -> None:
    """按 pattern 汇总输出: 先按 pattern 字典序, 相同 pattern 按 repeats count 排序"""
    header = ["pattern", "repeats_count", "seq_name", "start", "end"]
    out.write("\t".join(header) + "\n")
    for r in sorted(repeats, key=lambda r: (normalize_motif(r.motif), r.copy_number,
                                            r.seq_name, r.start)):
        out.write("\t".join([
            normalize_motif(r.motif),
            f"{r.copy_number:g}",
            r.seq_name, str(r.start), str(r.end),
        ]) + "\n")


# ---------------------------------------------------------------- main


def main():
    parser = argparse.ArgumentParser(description="STR/SSR/homopolymer region detection based on trfrs")
    parser.add_argument("inputs", nargs="+", help="fasta 文件或目录")
    parser.add_argument("--mode", choices=["str", "ssr"], default="str",
                        help="精筛模式: str=直接用粗筛结果; ssr=整段干净则保留, 否则提取区间内干净的 "
                             "unit 重复子区间 (default: str)")
    parser.add_argument("--min-copies", type=int, default=3,
                        help="ssr 模式下子区间的最小 unit 拷贝数 (default: 3)")
    parser.add_argument("--minscore", type=int, default=9, help="trfrs 粗筛最小得分 (default: 9)")
    parser.add_argument("--patterns", action="store_true",
                        help="只输出 SSR pattern 汇总: 按 pattern 字典序排序, "
                             "相同 pattern 按 repeats count 排序")
    parser.add_argument("--out", default=None, help="输出 tsv 路径, 缺省打印到 stdout")
    args = parser.parse_args()

    fasta_files = collect_fasta_files(args.inputs)
    refs: Dict[str, str] = {}
    for fa in fasta_files:
        refs.update(parse_fasta(fa))

    with tempfile.TemporaryDirectory(prefix="str_ssr_") as tmp:
        csv_path = run_trfrs(fasta_files, args.minscore, pathlib.Path(tmp))
        repeats = parse_trfrs_csv(csv_path)

    print(f"[info] coarse repeats: {len(repeats)}", file=sys.stderr)

    repeats = resolve_overlaps(repeats)
    print(f"[info] after overlap resolve: {len(repeats)}", file=sys.stderr)

    repeats = fine_filter(repeats, refs, args.mode, args.min_copies)
    print(f"[info] after {args.mode} fine filter: {len(repeats)}", file=sys.stderr)

    writer = write_patterns if args.patterns else write_result
    if args.out:
        with open(args.out, "w") as f:
            writer(repeats, f)
        print(f"[info] result written to {args.out}", file=sys.stderr)
    else:
        writer(repeats, sys.stdout)


if __name__ == "__main__":
    main()
