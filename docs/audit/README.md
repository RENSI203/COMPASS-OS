# 发布审查汇总包 ｜ release review bundle

**软件判定**：READY WITH LIMITATIONS ｜ **论文声明**：2 项已核实 / 5 项待核实（见 FINAL_RELEASE_READINESS.md §6）

| 文件 | 大小 | SHA256(前16) | 内容 |
|---|---|---|---|
| `00_CANDIDATE_BASELINE.md` | 28 KB | `f19e056cff9c9ff0` | 阶段0 候选冻结、版本裁定、release_audit 职责与裁定 |
| `01_DATA_LIMITATIONS.md` | 32 KB | `8049ac528deab5b5` | 阶段1 数据限制与适用范围（10 项） |
| `03_CORE_FUNCTIONS.md` | 17 KB | `a8247293dfbd7122` | 阶段3 主干功能与数值一致性（含容差术语纪律） |
| `04_DERIVED_SCRIPTS.md` | 25 KB | `1b2bb4c7959673c5` | 阶段4 衍生脚本与结果生产链 |
| `05_MAIN_RESULTS_CLAIMS.md` | 22 KB | `970857a274930d1b` | 阶段5 主要结果与论文声明 claim–evidence（含外部性能溯源核对） |
| `CANDIDATE_MANIFEST.tsv` | 22 KB | `1582c73bed163aa8` | 阶段0 候选 manifest |
| `FINAL_ARTIFACT_MANIFEST.tsv` | 0 KB | `36e21d129652d489` | 最终**分发产物** manifest（wheel/sdist SHA256，针对最终 wheel） |
| `FINAL_CANDIDATE_MANIFEST.tsv` | 24 KB | `f3790ccbcc77aac0` | 最终候选源码 manifest（213 文件） |
| `FINAL_RELEASE_READINESS.md` | 14 KB | `b237ecefaf2446f0` | 最终发布就绪评估：软件判定 + 论文声明核实状态 |

## 阶段 2 说明

**阶段 2（依赖项与安装审查）未作为独立报告执行**：其核心内容已由阶段 0
（wheel 构建、METADATA 依赖、可重复性、顺序要求）与阶段 1（输入尺度与依赖边界）覆盖；
最终分发与安装验收在 `FINAL_RELEASE_READINESS.md` §3.1 完成（针对最终 wheel、全新环境、仓库外）。
如需独立阶段 2 报告，应另行补做。
