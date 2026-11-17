#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""比较两个 merged.gsmm2-hp-tr-v3-aggr.csv 在指定 tag 下 ratio_within_motif_tag 的相对变化。

适用场景：同一 Barcode 的 A/B 两侧 smc_all_reads-metric 目录（baseline vs model），
想看某类 tag（常见是字面值 "pure+mixed"）在 correct call 行上
ratio_within_motif_tag 的变化百分比。

筛选条件（两条都必须满足）：
  tag == --tag（默认 "pure+mixed"，是一个字面取值，不是 pure 与 mixed 的并集）
  called == true_cnt（correct call）

join 口径：按 (motif, ref, true_base, true_cnt, called, tag) 做 inner join。
因为同 motif 同 tag 下只有 called == true_cnt 这一行满足条件，所以输出的
每一行对应一个 motif。
行顺序：沿用 base 文件中的原始行序；先校验两个文件筛选后的行序一致。
单边独有的键打印告警，不静默丢弃。

输出列：
  motif, ref, true_base, true_cnt, called, tag,
  num_base, num_new,
  ratio_within_motif_tag_base, ratio_within_motif_tag_new,
  delta_pct
delta_pct = (new - base) / base * 100，输出为固定 4 位小数的字符串，
形如 "1.2345%"；base 为 0 时输出 "NA"。

除 CSV 外，脚本直接打印「提升 Top N」与「下降 Top N」两张表。
两张表按绝对变化排序：abs_change_pp = (new - base) * 100，单位百分点。
N 由 --top 指定，默认 10。

用法：
  motif_repeat_count_ab.py BASE.csv NEW.csv
  motif_repeat_count_ab.py BASE.csv NEW.csv --tag pure -o out/x.csv
  motif_repeat_count_ab.py BASE.csv NEW.csv --top 20
"""

import argparse
import os
import re
import sys

import polars as pl

COL = "ratio_within_motif_tag"
KEYS = ["motif", "ref", "true_base", "true_cnt", "called", "tag"]
REQUIRED = set(KEYS) | {"num", COL}
OUT_COLS = ["motif", "ref", "true_base", "true_cnt", "called", "tag",
            "num_base", "num_new",
            f"{COL}_base", f"{COL}_new", "delta_pct"]

HERE = os.path.dirname(os.path.abspath(__file__))


def describe(path):
    """从 smc_all_reads-metric 目录名里取出 (barcode, side)；取不到则退回目录名。"""
    parent = os.path.basename(os.path.dirname(os.path.abspath(path)))
    m = re.match(r".*?\.(Barcode\d+)\.(.+)\.smc_all_reads-metric$", parent)
    if m:
        return m.group(1), m.group(2)
    return None, parent


def load(path, tag, side):
    if not os.path.exists(path):
        sys.exit(f"[error] {side} 文件不存在: {path}")
    df = pl.read_csv(path, separator="\t")
    missing = REQUIRED - set(df.columns)
    if missing:
        sys.exit(f"[error] {side} 缺少列 {sorted(missing)}: {path}")

    all_tag = df.filter(pl.col("tag") == tag)
    kept = all_tag.filter(pl.col("called") == pl.col("true_cnt"))
    if kept.height == 0:
        sys.exit(f"[error] {side} 在 tag={tag!r} 下没有 called==true_cnt 的行: {path}")

    dup = kept.select(KEYS).is_duplicated().sum()
    if dup:
        sys.exit(f"[error] {side} 在 tag={tag!r} 下筛选后键重复 {dup} 行: {path}")

    print(f"{side:<5}: tag={tag!r} {all_tag.height} 行 -> called==true_cnt {kept.height} 行")
    return kept.with_row_index("_ord"), [tuple(r) for r in kept.select(KEYS).iter_rows()]


def fmt4(v):
    """固定 4 位小数 + %；None -> NA；-0.0000 -> 0.0000。"""
    if v is None:
        return "NA"
    if round(v, 4) == 0:
        v = 0.0
    return f"{v:.4f}%"


def print_top(joined, n):
    """按绝对变化（百分点）打印提升与下降 Top N。"""
    rows = []
    for r in joined.iter_rows(named=True):
        rb, rn = r[f"{COL}_base"], r[f"{COL}_new"]
        rows.append((r["motif"], r["num_base"], r["num_new"], rb, rn,
                     (rn - rb) * 100, r["delta_pct"]))

    ups = sorted((x for x in rows if x[5] > 0), key=lambda x: x[5], reverse=True)[:n]
    downs = sorted((x for x in rows if x[5] < 0), key=lambda x: x[5])[:n]

    def show(title, sel):
        print(f"\n{title}（按绝对变化 pp 排序）")
        print(f"  {'#':>2}  {'motif':<8} {'num (base→new)':<16} "
              f"{'ratio (base→new)':<23} {'abs_pp':>10}  delta_pct")
        for i, (m, nb, nn, rb, rn, a, d) in enumerate(sel, 1):
            num_s = f"{nb}→{nn}"
            ratio_s = f"{rb:.6f}→{rn:.6f}"
            print(f"  {i:>2}  {m:<8} {num_s:<16} {ratio_s:<23} {a:>+10.4f}  {d}")

    show(f"提升 Top {n}", ups)
    show(f"下降 Top {n}", downs)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("base", help="baseline 的 merged.gsmm2-hp-tr-v3-aggr.csv")
    ap.add_argument("new", help="对比侧的 merged.gsmm2-hp-tr-v3-aggr.csv")
    ap.add_argument("--tag", default="pure+mixed",
                    help='tag 列的字面取值，默认 "pure+mixed"')
    ap.add_argument("--top", type=int, default=10,
                    help="打印的提升/下降 Top N 行数，默认 10")
    ap.add_argument("-o", "--out", default=None, help="输出 csv（默认放脚本目录的 out/）")
    args = ap.parse_args()

    base_bc, base_side = describe(args.base)
    new_bc, new_side = describe(args.new)
    print(f"base : {args.base}\nnew  : {args.new}")
    if base_bc and new_bc and base_bc != new_bc:
        print(f"[warn] barcode 不一致: base={base_bc} new={new_bc}")

    base, base_order = load(args.base, args.tag, "base")
    new, new_order = load(args.new, args.tag, "new")

    if base_order != new_order:
        print("[warn] 两个文件筛选后的行序不一致，输出按 base 的行序")
        same = len(set(base_order) & set(new_order))
        print(f"[warn] 行序一致的前缀长度不做保证；共同键 {same} 个")

    only_base = set(base_order) - set(new_order)
    only_new = set(new_order) - set(base_order)
    if only_base or only_new:
        print(f"[warn] 只在 base 里的键 ({len(only_base)}): {sorted(only_base)}")
        print(f"[warn] 只在 new  里的键 ({len(only_new)}): {sorted(only_new)}")
    else:
        print(f"键集合一致，共 {len(base_order)} 行")

    joined = (base
              .select(KEYS + ["_ord", "num", COL])
              .rename({"num": "num_base", COL: f"{COL}_base"})
              .join(new.select(KEYS + ["num", COL])
                        .rename({"num": "num_new", COL: f"{COL}_new"}),
                    on=KEYS, how="inner")
              .sort("_ord").drop("_ord"))

    b, n = pl.col(f"{COL}_base"), pl.col(f"{COL}_new")
    joined = joined.with_columns(
        pl.when(b == 0)
          .then(None)
          .otherwise(((n - b) / b) * 100)
          .map_elements(fmt4, return_dtype=pl.String)
          .alias("delta_pct"))

    deltas = [float(x.rstrip("%")) for x in joined["delta_pct"] if x != "NA"]
    up = sum(1 for x in deltas if x > 0)
    down = sum(1 for x in deltas if x < 0)
    print(f"delta: 上升 {up} 行，下降 {down} 行，不变 {len(deltas) - up - down} 行")

    print_top(joined, args.top)

    out_csv = args.out or os.path.join(
        HERE, "out",
        f"{base_bc or 'all'}.{new_side}_vs_{base_side}."
        f"{args.tag.replace('+', '-')}.{COL}-ab.csv")
    os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
    joined.select(OUT_COLS).write_csv(out_csv)
    print(f"\n写出: {out_csv}  ({joined.height} 行 × {len(OUT_COLS)} 列)")


if __name__ == "__main__":
    main()
