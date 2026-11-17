# str_motif_ab

STR 相关 A/B 评估脚本。目录里有两个脚本。每个脚本读取一类上游指标文件，输出一份 A/B 对比表。

## 1. 上游流程

两个脚本不直接读 BAM。它们读 `reads_quality_stats_hp_tr_v2.py` 或 `reads_quality_stats_hp_tr_v3.py` 的聚合结果。

```
BAM + ref FASTA
  |
  +-- reads_quality_stats_hp_tr_v2.py  --> merged.gsmm2-hp-tr-aggr.csv
  |                                          --> motif_accuracy_ab.py
  |
  +-- reads_quality_stats_hp_tr_v3.py  --> merged.gsmm2-hp-tr-v3-aggr.csv
                                             --> motif_repeat_count_ab.py
```

上游命令示例（v2）：

```bash
python /root/projects/gsda/third_party/gseda/src/gseda/ppl/reads_quality_stats_hp_tr_v2.py \
  --bams 20260723_250302Y0001_Run0003_called-demuxed.Barcode01.2026Q3StrModel.StrPoli-gapaffine.smc_all_reads.bam \
  --refs merged_output/Barcode01-ref.fasta \
  --np-range 5:100000
```

上游命令示例（v3）：

```bash
python /root/projects/gsda/third_party/gseda/src/gseda/ppl/reads_quality_stats_hp_tr_v3.py \
  --bams 20260723_250302Y0001_Run0003_called-demuxed.Barcode01.2026Q3StrModel.StrPoli-gapaffine.smc_all_reads.bam \
  --refs merged_output/Barcode01-ref.fasta \
  --np-range 5:100000
```

两条命令都把结果写到 BAM 所在目录下的 `<bam 文件名去掉扩展名>-metric/`：

| 上游脚本 | 聚合文件 | 其他输出 |
|---|---|---|
| `reads_quality_stats_hp_tr_v2.py` | `merged.gsmm2-hp-tr-aggr.csv` | `<stem>.gsmm2-hp-tr-fact.csv` |
| `reads_quality_stats_hp_tr_v3.py` | `merged.gsmm2-hp-tr-v3-aggr.csv` | `<stem>.gsmm2-hp-tr-v3-fact.csv`、`hp-tr-v3-pure-mixed-ratio.png` |

## 2. v2 与 v3 的关注点

| 上游脚本 | 关注点 | 用的指标 |
|---|---|---|
| v2 | STR 区域的碱基准确率 | `eq_rate = eq / aligned_span`，其中 `aligned_span = eq + diff + ins + del` |
| v3 | STR 预估的 repeats 个数是否正确 | `called`（模型估计的 repeats 数）相对 `true_cnt` / `ref`（参考 repeats 数）的分布。`ratio_within_motif_tag` 表示同一 motif、同一 tag 下，落在某个 `called` 取值上的 reads 占比 |

v3 在聚合前把 `called` 截断到 `[ref - 5, ref + 5]`。超出范围的 `called` 值都并到边界。

脚本名称按关注点区分：

- `motif_accuracy_ab.py` 处理 v2 的结果，比较 STR 区域准确率。
- `motif_repeat_count_ab.py` 处理 v3 的结果，比较 repeats 个数估计。

## 3. motif_accuracy_ab.py

按 motif 做 inner join，比较 A/B 两侧的 STR 区域准确率。

### 输入

两侧各一份 `merged.gsmm2-hp-tr-aggr.csv`。必需列：`motif`、`aligned_span`、`eq_rate`。同侧 `motif` 不得重复。

### 输出列

| 列 | 含义 |
|---|---|
| `motif` | motif 字符串 |
| `aligned_span_baseline` / `aligned_span_model` | 两侧比对跨度（覆盖量），整数 |
| `Eq_rate_baseline` / `Eq_rate_model` | 两侧碱基一致率，直接沿用源文件数值 |
| `eq_rate_rel_change_pct` | `(Eq_model - Eq_base) / Eq_base * 100` |
| `err_rate_rel_change_pct` | `(Eq_base - Eq_model) / (1 - Eq_base) * 100`，负值表示错误率下降 |

`eq_rate_rel_change_pct` 和 `err_rate_rel_change_pct` 最多保留 6 位小数。`Eq_rate_base == 1` 时错误率相对变化无定义，留空。

### 用法

```bash
# 直接给两份 csv
python motif_accuracy_ab.py --baseline a.csv --model b.csv -o out/x.csv

# 按目录名推导两侧路径
python motif_accuracy_ab.py --data-dir /data/str-ab --run RUN1 --barcode Barcode08

# 覆盖 model 侧的 side 字段（例如 2026Q3Model）
python motif_accuracy_ab.py --data-dir /data1/ccs_data/smc-ab-test-dataset/corn-str \
  --run 20260930_250302Y0006_Run0001_called-demuxed \
  --barcode Barcode08 --model-tag 2026Q3Model
```

推导路径的规则：

```
<data-dir>/<run>.<barcode>.<side>.smc_all_reads-metric/merged.gsmm2-hp-tr-aggr.csv
```

- baseline 侧的 `<side>` 取 `--base-tag`，默认 `baseline`。
- model 侧的 `<side>` 取 `--model-tag`，默认 `model`。
- 不给 `--baseline` 或 `--model` 时，`--data-dir`、`--run`、`--barcode` 都必须给出。

默认输出：脚本目录下 `out/<barcode>.motif-ab.csv`。`barcode` 从 `--barcode` 取；取不到时从输入路径的目录名解析。

## 4. motif_repeat_count_ab.py

比较两个 `merged.gsmm2-hp-tr-v3-aggr.csv` 在指定 tag 下 `ratio_within_motif_tag` 的相对变化。

### 输入

两个文件路径，作为位置参数：`BASE.csv NEW.csv`。

必需列：`motif`、`ref`、`true_base`、`true_cnt`、`called`、`tag`、`num`、`ratio_within_motif_tag`。

### 筛选条件

两条必须同时满足：

1. `tag == --tag`。默认 `"pure+mixed"`。v3 会把全部 tag 合并后的结果打上这个字面值。它不是 pure 与 mixed 两个 tag 的筛选并集。
2. `called == true_cnt`。即 correct call。

同一 motif 同一 tag 下只有 `called == true_cnt` 这一行满足条件，所以输出的每一行对应一个 motif。

### 输出列

| 列 | 含义 |
|---|---|
| `motif`、`ref`、`true_base`、`true_cnt`、`called`、`tag` | 连接键 |
| `num_base` / `num_new` | 两侧的 `num` |
| `ratio_within_motif_tag_base` / `ratio_within_motif_tag_new` | 两侧的 `ratio_within_motif_tag` |
| `delta_pct` | `(new - base) / base * 100`，固定 4 位小数并带 `%`；`base == 0` 时输出 `NA` |

脚本同时打印「提升 Top N」和「下降 Top N」两张表。排序按 `abs_change_pp = (new - base) * 100`，单位是百分点。`N` 由 `--top` 指定，默认 10。

### 用法

```bash
python motif_repeat_count_ab.py BASE.csv NEW.csv
python motif_repeat_count_ab.py BASE.csv NEW.csv --tag pure -o out/x.csv
python motif_repeat_count_ab.py BASE.csv NEW.csv --top 20
```

默认输出：脚本目录下

```
out/<barcode>.<new_side>_vs_<base_side>.<tag>.ratio_within_motif_tag-ab.csv
```

`barcode` 和 `side` 从输入路径的 `smc_all_reads-metric` 目录名解析。解析不到时，`barcode` 用 `all`，`side` 用目录名。

## 5. A/B 对比的前提

两侧必须使用相同的上游参数。至少包括：

- `--np-range`
- `--rq-range`（如果使用）
- `--ref-anchored`（如果使用）
- 参考 FASTA

参数不同时，`delta_pct` 会混入方法和参数的差异。

## 6. 依赖

两个脚本都用 `polars`。在项目的 `gseda` conda 环境中运行：

```bash
/root/miniconda3/envs/gseda/bin/python motif_accuracy_ab.py --help
/root/miniconda3/envs/gseda/bin/python motif_repeat_count_ab.py --help
```

## 7. 注意

- `out/` 下现有的 CSV 是旧版脚本生成的。表头仍含 `2026Q3Model`（例如 `aligned_span_2026Q3Model`）。用当前脚本重新生成才会得到 `_model` 列名。
- 推导路径要求 BAM stem 的格式为 `<run>.<barcode>.<side>.smc_all_reads`。格式不同时，请直接传 `--baseline` 和 `--model`。
