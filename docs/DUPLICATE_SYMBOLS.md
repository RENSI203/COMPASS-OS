# 重复 gene symbol 审计（release blocker 核查结论）

> 审计脚本：`tools/audit_duplicate_symbols.py`；数据表：`docs/duplicate_gene_symbols.tsv`；
> 契约测试：`tests/test_duplicate_genes.py`。

## 结论（一句话）

**冻结词表 ``feature_name`` 15,672 项全部唯一、零重复** ⇒ 官方 package 的
「列名 = gene symbol」输入契约**不存在不可逆映射歧义**；曾报告的「~55 个重复」
位于**源矩阵构建阶段**，不影响 package 输入。

## 1. 词表侧（模型实际使用的 15,672 个 token）

| 项 | 值 |
| --- | --- |
| 15,672 个 feature 中 unique gene symbols | **15,672** |
| duplicate rows | **0** |
| duplicate symbol 种类数 | **0** |
| ``gene_tokens_long.json`` 的基因项 | 15,672（+`CLS`/`CANCER`），与 `feature_name` **逐项相同**，同样零重复 |

⇒ 一个用户 symbol 列 → **恰好一个 frozen token**；既不需要 broadcast，也不需要 collapse。

## 2. 源矩阵侧（为什么会有「55 个重复」的说法）

原项目 ``01_prepare_expression.py`` 的 `gene_match_report.json` 记录
``n_duplicated_symbols = 55``。逐项复现结果：

| 项 | 值 |
| --- | --- |
| TCGA RSEM 源矩阵基因行 | 60,498（索引 = Ensembl ID） |
| 经 `probeMap_gencode.v23` 映射且符号在词表内 | 15,709 行 |
| 唯一符号 | **15,654** |
| **重复符号种类** | **55**（每种恰好 2 行 ⇒ 55 个多余行） |
| 词表中完全没有源行的基因 | 18（矩阵中填 `TPM = 0`，`coverage_c = 0.998851`） |
| 与 `01` 报告核对 | **55/55 完全一致，0 处不符** |

重复对的构成：多为 ``ENSG…`` + ``ENSGR…``（PAR_Y / 参考重复行），少数为旁系同源或通读对。例：

| symbol | 两个 Ensembl ID |
| --- | --- |
| `ACSL6` | `ENSG00000164398.12` ; `ENSG00000281938.1` |
| `AKAP17A` | `ENSG00000197976.10` ; `ENSGR0000197976.10` |
| `ASMT` | `ENSG00000196433.11` ; `ENSGR0000196433.11` |
| `C11orf71` | `ENSG00000180425.10` ; `ENSG00000282682.1` |

## 3. 原项目的消解规则（训练数据构建口径）

``01_prepare_expression.py:172-179``：

```python
if len(lst) > 1:
    best = int(np.argmax([float(r.mean()) for r in lst]))   # 保留行均值最大者
    chosen = lst[best]
```

即 **保留平均表达最大的那一行**，另一行丢弃。`docs/duplicate_gene_symbols.tsv`
记录了每个重复符号的 `original_source_ids`、`kept_source_id_max_mean`、
`dropped_source_id` 与 `max_over_min_row_mean`（数值审计：流式扫描原始 1.32 GB 矩阵，
逐 Ensembl 行统计；明细见 `docs/duplicate_source_row_stats.tsv`）。

**数值审计结果**（55 个符号，保留行均值 / 丢弃行均值）：

| 比值 | 符号数 | 含义 |
| --- | --- | --- |
| `> 2` | **52 / 55** | 两行表达差异显著；被丢弃的一行基本不表达（多数比值为 `inf`，即丢弃行均值 ≤ 0） |
| `1.01 – 1.2` | 1 / 55 | 两行近乎等价，选哪一行几乎无影响 |

⇒ 这 55 个符号占词表 **0.35%**（55/15,672）；`max-mean` 规则在 52/55 的情形下
是"选真正表达的那一条"，**选择本身不模糊**。

## 4. 对 package 的影响

| 场景 | 结论 |
| --- | --- |
| 用户按 symbol 提供**一列** | 与冻结 token **一一对应**，与生产口径完全一致（golden 测试逐位通过） |
| 用户按 Ensembl 自行 **sum/average** 那 55 个符号的多行 | 取值**可能**与该 token 训练时所见（最大均值行）不同 ⇒ README 已提示按 symbol 提供 |
| 用户提供**重复列名** | package **报 `InputError`**，不猜测（`test_duplicate_columns_are_rejected_not_guessed`） |

⇒ **不构成 release blocker**：symbol-only 契约无歧义，且重复列名被显式拒绝。
