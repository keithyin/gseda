#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""按 motif join 两张 merged.gsmm2-hp-tr-aggr.csv（A/B 两侧），产出 motif 级 A/B 对比 CSV。

输入：同一 Barcode 下两侧 smc_all_reads-metric 目录里的 merged.gsmm2-hp-tr-aggr.csv
（列：motif, eq, diff, ins, del, aligned_span, eq_rate, ...，其中 eq_rate = eq / aligned_span）。
A/B 两侧统一称为 baseline 与 model。

输出列（固定头，见 OUT_COLS）：
  motif
  aligned_span_baseline / aligned_span_model   两侧比对跨度（覆盖量），整数
  Eq_rate_baseline      / Eq_rate_model        两侧碱基一致率
  eq_rate_rel_change_pct  = (Eq_model - Eq_base) / Eq_base * 100
  err_rate_rel_change_pct = ((1 - Eq_model) - (1 - Eq_base)) / (1 - Eq_base) * 100
                            = (Eq_base - Eq_model) / (1 - Eq_base) * 100
                            （错误率的相对变化，负值 = 错误率下降）
两个相对变化百分比最多保留 6 位小数（去掉多余的 0）；Eq_rate 两列直接沿用源文件数值；
Eq_rate_base == 1 时错误率相对变化无定义，留空。

join 口径：motif 上的 inner join，保持 baseline 原始行序；单边独有的 motif 打印告警（不静默丢）。

路径推导：不给 --baseline 或 --model 时，按
  <data-dir>/<run>.<barcode>.<side>.smc_all_reads-metric/merged.gsmm2-hp-tr-aggr.csv
推导。baseline 侧的 side 取 --base-tag（默认 baseline），model 侧取 --model-tag（默认 model）。

用法：
  motif_accuracy_ab.py --baseline a.csv --model b.csv -o out/x.csv
  motif_accuracy_ab.py --data-dir /data/str-ab --run RUN1 --barcode Barcode08
  motif_accuracy_ab.py --data-dir /data/str-ab --run RUN1 --barcode Barcode08 --model-tag new-model
"""

import argparse
import os
import re
import sys

import polars as pl

OUT_COLS = ["motif",
            "aligned_span_baseline", "aligned_span_model",
            "Eq_rate_baseline", "Eq_rate_model",
            "eq_rate_rel_change_pct", "err_rate_rel_change_pct"]

HERE = os.path.dirname(os.path.abspath(__file__))
METRIC_DIR_RE = re.compile(r".*?\.(Barcode\d+)\.(.+)\.smc_all_reads-metric$")


def metric_csv(data_dir, run, barcode, side):
    return os.path.join(
        data_dir, f"{run}.{barcode}.{side}.smc_all_reads-metric", "merged.gsmm2-hp-tr-aggr.csv")


def describe(path):
    """从 smc_all_reads-metric 目录名里取出 (barcode, side)；取不到则返回 (None, None)。"""
    parent = os.path.basename(os.path.dirname(os.path.abspath(path)))
    m = METRIC_DIR_RE.match(parent)
    if not m:
        return None, None
    return m.group(1), m.group(2)


def load(path, tag):
    if not os.path.exists(path):
        sys.exit(f"[error] {tag} 文件不存在: {path}")
    df = pl.read_csv(path, separator="\t")
    missing = {"motif", "aligned_span", "eq_rate"} - set(df.columns)
    if missing:
        sys.exit(f"[error] {tag} 缺少列 {sorted(missing)}: {path}")
    if df["motif"].n_unique() != df.height:
        sys.exit(f"[error] {tag} 的 motif 有重复: {path}")
    return df.with_row_index("_ord")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", default=None, help="baseline 侧 aggr csv")
    ap.add_argument("--model", default=None, help="model 侧 aggr csv")
    ap.add_argument("--data-dir", default=None,
                    help="推导路径时使用：smc-ab-test-dataset 根目录")
    ap.add_argument("--run", default=None, help="推导路径时使用：run 目录名前缀")
    ap.add_argument("--barcode", default=None, help="推导路径时使用：barcode 名，如 Barcode08")
    ap.add_argument("--base-tag", default="baseline",
                    help="baseline 侧 smc_all_reads-metric 目录名的 side 字段，默认 baseline")
    ap.add_argument("--model-tag", default="model",
                    help="model 侧 smc_all_reads-metric 目录名的 side 字段，默认 model")
    ap.add_argument("-o", "--out", default=None,
                    help="输出 csv（默认脚本目录下 out/<barcode>.motif-ab.csv）")
    args = ap.parse_args()

    def derive(side_tag, explicit, label):
        if explicit:
            return explicit
        missing = [f"--{name}" for name, value in
                   (("data-dir", args.data_dir), ("run", args.run), ("barcode", args.barcode))
                   if not value]
        if missing:
            sys.exit(f"[error] 未给 --{label}，且推导路径缺少 {', '.join(missing)}")
        return metric_csv(args.data_dir, args.run, args.barcode, side_tag)

    base_csv = derive(args.base_tag, args.baseline, "baseline")
    model_csv = derive(args.model_tag, args.model, "model")
    barcode = args.barcode or describe(base_csv)[0] or describe(model_csv)[0]
    out_csv = args.out or os.path.join(HERE, "out", f"{barcode or 'all'}.motif-ab.csv")

    print(f"baseline: {base_csv}")
    print(f"model   : {model_csv}")

    base = load(base_csv, "baseline")
    model = load(model_csv, "model")

    only_base = sorted(set(base["motif"]) - set(model["motif"]))
    only_model = sorted(set(model["motif"]) - set(base["motif"]))
    if only_base or only_model:
        print(f"[warn] 只在 baseline 里的 motif ({len(only_base)}): {only_base}")
        print(f"[warn] 只在 model 里的 motif ({len(only_model)}): {only_model}")
    else:
        print(f"motif 集合一致，共 {base.height} 个")

    joined = (base
              .join(model.select("motif", "aligned_span", "eq_rate"),
                    on="motif", how="inner", suffix="_model")
              .rename({"aligned_span": "aligned_span_baseline",
                       "eq_rate": "Eq_rate_baseline"})
              .rename({"eq_rate_model": "Eq_rate_model"})
              .sort("_ord").drop("_ord"))

    eq_b, eq_m = pl.col("Eq_rate_baseline"), pl.col("Eq_rate_model")
    joined = joined.with_columns([
        (((eq_m - eq_b) / eq_b) * 100).alias("eq_rate_rel_change_pct"),
        pl.when(1 - eq_b == 0)
          .then(None)
          .otherwise(((eq_b - eq_m) / (1 - eq_b)) * 100)
          .alias("err_rate_rel_change_pct"),
    ])

    # 百分比最多 6 位小数（去掉多余的 0），-0 -> 0，无定义的留空
    def fmt6(v):
        if v is None:
            return ""
        t = f"{round(v, 6):.6f}".rstrip("0").rstrip(".")
        return "0" if t in ("", "-0") else t

    for c in ("eq_rate_rel_change_pct", "err_rate_rel_change_pct"):
        joined = joined.with_columns(
            pl.col(c).map_elements(fmt6, return_dtype=pl.String).alias(c))

    os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
    joined.select(OUT_COLS).write_csv(out_csv)
    print(f"\n写出: {out_csv}  ({joined.height} 行 × {len(OUT_COLS)} 列)")


if __name__ == "__main__":
    main()
