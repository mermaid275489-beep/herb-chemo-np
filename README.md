# 中药–化疗网络药理预测

在浏览器里完成的网络药理学分析：输入中药或复方、化疗药或方案、癌种，计算中药靶点与以下基因集之间的重叠和网络邻近度：

- 癌种基因模块
- 化疗药作用靶点（PD）
- 化疗药代谢/转运基因（PK）
- 29 个不良反应模块

此外还做通路富集，对规则库收录的组合给出 PK/PD 机制推理。

网页是纯静态页面（`docs/`），由 GitHub Pages 托管，打开即用，不需要登录，也不需要后端。内置列表以外的癌种和药物，会从网页直接实时查询 Open Targets。

*Browser-based network pharmacology for traditional Chinese medicine (TCM) herbs combined with chemotherapy. Enter herbs or formulas, drugs or regimens, and a cancer type. The page computes target overlap and Guney network proximity against cancer, drug-target, drug-ADME and 29 adverse-event gene modules, and runs pathway enrichment.*

## 数据

| 内容 | 来源 | 处理 |
|---|---|---|
| 中药–成分–靶点 | TCMSP 2.3（old.tcmsp-e.com），502 味药 | 默认 OB ≥ 30%、DL ≥ 0.18；靶点名经 UniProt 映射为人类基因（已审阅条目，名称或词集完全一致），非人类靶点不计 |
| 蛋白互作网络 | STRING v12.0，人类，综合得分 ≥ 900 | 12,323 个基因、100,856 条边，同时作为统计背景 |
| 癌种基因 | Open Targets Platform | 关联得分前 300 个靶点；内置 50 个癌种，其余疾病实时查询 |
| 化疗药 | Open Targets（作用机制、FAERS）、final27 整理的 PD 核心基因、规则库 | 113 种常用抗肿瘤药，支持 FOLFOX 等方案缩写 |
| 不良反应模块 | final27 分析 | 29 个模块：人工整理的核心机制基因（标注促损伤/保护性），加 HPO 表型基因作为扩展 |
| 通路 | MSigDB v2026.1：Reactome、GO 生物过程 | 因许可限制未纳入 KEGG |

## 方法

与 final27 分析一致，`docs/app.js` 中有注释：

- **重叠**：超几何检验。
- **网络邻近度**：closest 距离（Guney 等，2016）。用度匹配随机抽样估计零分布：度窗口逐步放宽，直到候选节点 ≥ 20 个，重复 1,000 次，得到 Z 值和经验 P 值。
- **多重检验**：一次查询内的全部检验做 BH 校正，重叠和邻近度分开校正。
- **显著标准**：（重叠 > 0 且 q < 0.05）或（Z < −2 且 q < 0.05）。
- **富集**：只用位于模块内、或与模块基因直接相连的中药靶点。

`build/validate/compare_final27.js` 把 final27 某个试验的基因集输入网页的计算引擎，与论文结果逐项比对：

- 模块大小、重叠数、超几何 P、观测距离应完全一致。
- Z 值只因随机序列不同而有小幅差异。

```bash
python build/validate/make_input.py Sun_2024 COLORECTAL raw/val_input.json
node build/validate/compare_final27.js raw/val_input.json
```

## 重新构建数据

`build/` 中的脚本依赖本地项目里的 final27 输出、STRING 缓存和 MSigDB 文件。这些文件不在本仓库中。

```bash
python build/fetch_tcmsp.py --workers 5      # 抓取 TCMSP（可断点续传）
python build/map_targets.py                  # 靶点名 → 人类基因
python build/fetch_opentargets.py            # 癌种、化疗药
python build/build_static.py                 # 网络、通路、AE 模块
python build/build_site_data.py              # 中药、化疗药、癌种、规则库
python -m http.server 8770 --directory docs  # 本地预览
```

## 使用须知

本工具用于科研假设生成。

- 网络关联不等于临床作用，也不能区分抑制与诱导、加重与减轻。
- TCMSP 的预测靶点存在假阳性。
- 结果不能作为用药依据。

各数据源的许可与引用要求以其官网为准：

- TCMSP：Ru 等，J Cheminform 2014
- STRING：CC BY 4.0
- Open Targets：CC0
- MSigDB：CC BY 4.0
- HPO：使用时需注明
