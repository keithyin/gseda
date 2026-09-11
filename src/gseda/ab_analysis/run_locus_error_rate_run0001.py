#!/usr/bin/env python3
"""按 barcode 驱动 locus_error_rate.py 的 A/B 分析并汇总。

Run0001 的 15 个真 barcode 各自有一代参考 (STR<p>.fa)。对每个 barcode:
  control = ctrl/BarcodeNN.fastq
  exp     = exp/BarcodeNN.fastq
  ref     = STR<p>.fa
调用 locus_error_rate.main_cli 出单 barcode 的三张 csv, 最后按 (barcode, refname, pos)
拼成跨 barcode 的 control / experiment / joined 三张总表。
"""

import os
import sys

import polars as pl

from gseda.ab_analysis.locus_error_rate import main_cli as ler_main

BASE = "/data1/ccs_data/str-optimization/second-batch-of-data"
EXP_DEMUXED = f"{BASE}/20260805_250804Y0004_Run0001/" \
              "20260805_250804Y0004_Run0001_called-barcode-v4-2026Q2Model/demuxed"
CTRL_DEMUXED = f"{BASE}/20260805_250804Y0004_Run0001/" \
               "20260805_250804Y0004_Run0001_called-barcode-v4-baseline/demuxed"
REF_DIR = f"{BASE}/STR第二批一代测序/STR第二批一代测序/merged_output"
MAP_TSV = f"{BASE}/plasmid_2_barcode.tsv"
RUN = "20260805_250804Y0004_Run0001"
OUT_BASE = f"{BASE}/20260805_250804Y0004_Run0001/locus_error_rate"

# 24标签-N -> BarcodeNN
def label_to_barcode(label: str) -> str:
    n = int(label.rsplit("-", 1)[1])
    return f"Barcode{n:02d}"


def load_pairs():
    rows = []
    with open(MAP_TSV) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("plasmid\t"):
                continue
            plasmid, label, run = line.split("\t")
            if run != RUN:
                continue
            rows.append((plasmid, label_to_barcode(label)))
    rows.sort(key=lambda x: int(x[1].replace("Barcode", "")))
    return rows


def main():
    os.makedirs(OUT_BASE, exist_ok=True)
    pairs = load_pairs()
    print(f"[driver] {len(pairs)} barcodes for {RUN}")

    for plasmid, barcode in pairs:
        ctrl_fq = os.path.join(CTRL_DEMUXED, f"{barcode}.fastq")
        exp_fq = os.path.join(EXP_DEMUXED, f"{barcode}.fastq")
        ref_fa = os.path.join(REF_DIR, f"STR{plasmid}.fa")
        outdir = os.path.join(OUT_BASE, barcode)
        print(f"\n==== {barcode} -> {plasmid} ====")
        ler_main([
            "--control", ctrl_fq,
            "--exp", exp_fq,
            "--ref", ref_fa,
            "--rq-thr", "0",       # FASTQ 无 rq 字段, 该阈值被 gsmm2 忽略
            "--outdir", outdir,
        ])

    # ---- 汇总: 每个 barcode 一张 csv, 加 barcode 列后 vstack ----
    def collect(name):
        frames = []
        for _, barcode in pairs:
            p = os.path.join(OUT_BASE, barcode, name)
            if not os.path.exists(p):
                continue
            df = pl.read_csv(p, separator="\t")
            frames.append(df.with_columns(pl.lit(barcode).alias("barcode")))
        if not frames:
            return None
        out = pl.concat(frames)
        out = out.select(["barcode"] + [c for c in out.columns if c != "barcode"])
        return out

    for name, tag in [
        ("control_locus_accuracy.csv", "control_all"),
        ("experiment_locus_accuracy.csv", "experiment_all"),
        ("locus_accuracy_joined.csv", "joined_all"),
    ]:
        df = collect(name)
        if df is None:
            continue
        dest = os.path.join(OUT_BASE, f"{tag}.csv")
        df.write_csv(dest, separator="\t")
        print(f"\n[driver] 汇总 {tag}: {df.height} rows -> {dest}")

    print(f"\n[driver] 全部完成。输出目录: {OUT_BASE}")


if __name__ == "__main__":
    main()
