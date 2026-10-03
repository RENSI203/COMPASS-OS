# 43 concept 分数是什么（语义审计摘要）

> 来源：原项目 `docs/COMPASS_concept_semantics_audit.md` 的结论摘要（该审计逐行核对上游
> 源码与冻结权重）。**本包不改模型**，只是把结论固化到文档与命名里。

## 结论（三句话）

1. **43 个「概念」不是表达量分数，而是模型内部的 learned latent 坐标。**
   基因集分数 = 编码器**嵌入**上的一个 `nn.Linear(32→1)`（**全 132 个基因集共用同一投影方向**）
   + 标量 bias，且 `ReLU` 被注释掉 ⇒ 符号完全不受约束；
   43 概念 = 其成员基因集分数的 **softmax 注意力凸组合**。
2. **训练目标里没有任何「概念分数 ↔ 基因集表达」的监督**：只有自监督对比/三元组损失
   与任务损失；也没有对概念分数做 Reference 差分/中心化 ⇒ **latent 轴的正负方向没有被生物学识别**。
3. ⇒ 概念名是 **semantic anchor（语义锚点）**，不是丰度/活性读数。
   后续机制解释**应以 `132 gene sets → pathway` 为主**，43 概念只能作为模型内部表征。

## 对本包的直接影响

| 项 | 要求 |
| --- | --- |
| 命名 | `signature_scores` / `concept_scores`；覆盖度只能叫 `concept_input_coverage`，禁止 confidence 系列用词 |
| 输出 | 只支持 association / prioritisation；**不输出** activation / suppression / causal mechanism |
| 文档 | README 与示例必须写明"hypothesis-generating features" |
| 能力 | 43 概念可作**模型内部表征**并列展示；机制优先级建议以 132 基因为主线 |

## 上游证据（源码位置）

| 事实 | 位置 |
| --- | --- |
| 43 由 132 经 softmax 注意力聚合（凸组合，权重非负且和为 1） | `compass/projector/cellpathwayaggregator.py:75`；`projector.py:197` |
| 基因集分数 = `nn.Linear(feature_dim→1)` 作用在**嵌入**上（不是表达量） | `compass/projector/genesetscorer.py:62` |
| ReLU 被注释 ⇒ 符号不受约束（关键证据） | `compass/projector/genesetscorer.py:76` |
| 无 Reference 差分/中心化施加在概念分数上 | `compass/projector/projector.py:183` |
| 审计实测：`softmax(w)·geneset_scores ≈ concept` Pearson 中位 1.0000、max\|Δ\| 中位 3.2e-08 | 本包对应测试：`test_concept_aggregation_identity` |
| 审计实测：**20/43 个概念在模型层面反向编码** | 原项目审计 PART C |
| 注意力集中度：最大权重中位 0.800（范围 0.288–1.000） | `concept_weights()` 可复核 |
