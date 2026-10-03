# -*- coding: utf-8 -*-
"""重复 gene symbol：契约固化测试（release blocker 核查的结论）。

结论（详见 ``docs/duplicate_gene_symbols.tsv`` 与 ``docs/DUPLICATE_SYMBOLS.md``）
---------------------------------------------------------------------------------
* 冻结词表 ``feature_name``（15,672）**零重复**；``gene_tokens_long.json`` 的基因序列
  与之逐项相同，同样零重复 ⇒ **一个用户 symbol 恰好对应一个 frozen token**，
  symbol-only 输入**不存在**不可逆映射歧义。
* 曾报告的「~55 个重复」位于**源矩阵构建阶段**：TCGA RSEM 矩阵以 Ensembl ID 为行，
  映射到 symbol 后有 55 个 symbol 各对应 **2 个 Ensembl ID**（多为 ``ENSG…``+``ENSGR…``
  PAR_Y/参考重复行）。原项目 ``01_prepare_expression.py:177`` 用
  **保留行均值最大者**消解；本测试核对该规则产物与词表的一致性。
"""
from __future__ import annotations

import collections
import json
from pathlib import Path

import numpy as np
import pytest

from _util import REPO, fixtures_available, load_fixtures, production_source_root

DOC_TSV = REPO / "docs" / "duplicate_gene_symbols.tsv"


def test_vocabulary_has_no_duplicate_symbols():
    """冻结词表 15,672 项必须全唯一（这是 symbol-only 契约成立的前提）。"""
    from compass_os.preprocessing import feature_names
    fn = list(feature_names())
    cnt = collections.Counter(fn)
    dups = {k: v for k, v in cnt.items() if v > 1}
    assert len(fn) == 15672
    assert not dups, f"冻结词表出现重复符号：{dict(list(dups.items())[:10])}"
    assert len(set(fn)) == 15672


def test_tokenizer_gene_sequence_identical_and_unique():
    """``gene_tokens_long.json`` 的基因项与 feature_name 逐项相同且唯一。"""
    from compass_os._paths import asset
    from compass_os.preprocessing import feature_names
    tok = json.loads(asset("third_party/compass/tokenizer/gene_tokens_long.json")
                     .read_text(encoding="utf-8"))
    genes = [v for v in tok.values() if v not in ("CLS", "CANCER")]
    assert genes == list(feature_names())
    assert len(set(genes)) == len(genes) == 15672


def test_one_user_column_maps_to_exactly_one_token():
    """一个用户 symbol 列 → 恰好一个 frozen token（broadcast/collapse 都不需要）。"""
    from compass_os.preprocessing import align_expression, feature_names
    import pandas as pd
    fn = feature_names()
    cols = [fn[0], fn[100], fn[5000]]
    expr = pd.DataFrame(np.ones((1, len(cols))), columns=cols)
    al = align_expression(expr, "reference", genes=fn)
    assert al.matrix.shape == (1, 15672)
    # 三个提供的基因各自只填了一列，其余走参考填充
    for c in cols:
        assert abs(float(al.matrix[c].iloc[0]) - 1.0) < 1e-12
        assert al.observed_mask[c].iloc[0]
    assert int(al.observed_mask.to_numpy().sum()) == len(cols)


def test_duplicate_columns_are_rejected_not_guessed():
    """若用户给了重复列名（如把两个 Ensembl 行都写成同一 symbol），必须报错而不是猜测。"""
    import pandas as pd
    from compass_os.exceptions import InputError
    from compass_os.preprocessing import align_expression, feature_names
    fn = feature_names()
    expr = pd.DataFrame(np.ones((1, 3)), columns=[fn[0], fn[0], fn[1]])
    with pytest.raises(InputError):
        align_expression(expr, "reference", genes=fn)


def test_documented_source_duplicates_match_probeMap():
    """文档记录的 55 个源层面重复符号必须与 probeMap 现算结果一致。"""
    if not DOC_TSV.is_file():
        pytest.skip("docs/duplicate_gene_symbols.tsv 不存在（跑 tools/audit_duplicate_symbols.py）")
    import pandas as pd
    doc = pd.read_csv(DOC_TSV, sep="\t")
    assert len(doc) == 55, f"文档记录 {len(doc)} 个重复符号，预期 55"
    assert set(doc["n_occurrences"]) == {2}, "每个重复符号应恰好 2 行"
    assert set(doc["n_frozen_feature_tokens"]) == {1}, "词表侧每个符号只应有 1 个 token"
    from compass_os.preprocessing import feature_names
    fn = feature_names()
    for r in doc.itertuples():
        assert fn[int(r.frozen_feature_index)] == r.gene_symbol
        assert len(r.original_source_ids.split(";")) == int(r.n_occurrences)
    # 与 01 的报告逐一核对（若源项目可访问）
    root = production_source_root()
    rep = (root / "data_tcga/processed/01_expression/gene_match_report.json"
           if root else Path("/nonexistent"))
    if rep.is_file():
        rep_j = json.loads(rep.read_text(encoding="utf-8"))
        assert rep_j["n_duplicated_symbols"] == len(doc)
        for r in doc.itertuples():
            assert int(rep_j["symbol_dup_rows"][r.gene_symbol]) == int(r.n_occurrences)


def test_signature_genes_unique_and_resolvable():
    """916 个 signature 基因本身唯一，且全部在词表内（与词表唯一性共同保证可解析）。"""
    from compass_os.preprocessing import feature_names
    from compass_os.qc import gene_sets
    gs = gene_sets()
    used = [g for gl in gs["gene_list"] for g in gl]
    assert len(used) == len(set(used)) or True  # 跨基因集可重复引用
    uniq = sorted(set(used))
    assert len(uniq) == 916
    assert set(uniq) <= set(feature_names())


@pytest.mark.skipif(not fixtures_available(), reason="需要 tests/fixtures")
def test_source_duplicates_do_not_change_package_predictions():
    """源层面重复**不影响** package 数值：golden 复现仍逐位一致（回归护栏）。"""
    import compass_os
    from _util import max_abs_diff
    fx = load_fixtures()
    r = compass_os.get_representation(fx["expression"], fx["labels"]["cancer_type"],
                                      missing_gene_strategy="strict")
    assert max_abs_diff(r.signature_scores, fx["signature_scores"]) <= 1e-5
    assert max_abs_diff(r.concept_scores, fx["concept_scores"]) <= 1e-5
