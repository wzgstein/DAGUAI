# Round 23 独立复核：complex 级配对重分析

**日期：** 2026-08-28  
**性质：** 事后审计性重分析，不是预注册确认分析。  
**输入：** 两个成功 GitHub Actions artifact：zero-value decoupled 与 mask-value decoupled。

## 1. 为什么重分析

原工作流的 headline 主要按基因均值排列表示。实际统计结构不是 182 个独立观测：

- 每条细胞系评估 91 个相同 target genes；
- 这些 target 被 `assign_complex` 归入 38 个 complex；
- K562 与 RPE1 共享 target ontology；
- 同一 complex 内的基因显著相关。

因此，本复核把 **complex** 作为外层重采样单位，重点估计同一 target 上的 `contextual gain - SG_Static gain`，而不是只比较两个独立均值。

## 2. 支持与方法

- common covered pool：172 genes；
- covered CORUM-eligible：101 genes；
- 实际评估：每条细胞系 91 targets；
- assigned complexes：38；
- 主方法：`raw_knn`；
- 汇总：先在每个 complex 内平均 K562/RPE1 与成员基因的配对差，再对 38 个 complex 等权平均；
- 区间：50,000 次 complex-cluster bootstrap；
- 检验：50,000 次 complex-level sign flip；
- `complex_id` 来自原脚本的单一 complex 指派；CORUM 重叠结构未被完全表达。

## 3. 关键表示的 raw gain

| representation | equal-complex mean gain | 95% CI | positive-complex fraction |
|---|---:|---:|---:|
| GenePT_NCBI_UniProt | 0.0601 | [0.0348, 0.0893] | 0.7895 |
| GO_C | 0.0511 | [0.0296, 0.0754] | 0.7368 |
| SG_Static | 0.0353 | [0.0169, 0.0579] | 0.7895 |
| SG_H00 | 0.0295 | [0.0121, 0.0507] | 0.7105 |
| SG_D12 | 0.0278 | [0.0105, 0.0489] | 0.6842 |
| SG_H03 | 0.0238 | [0.0055, 0.0459] | 0.6316 |
| SG_H06 | 0.0225 | [0.0090, 0.0387] | 0.6579 |
| SG_H12 | 0.0174 | [0.0046, 0.0330] | 0.5789 |
| SG_CLSM06 | 0.0168 | [0.0065, 0.0280] | 0.6842 |
| SG_MEANM12 | 0.0157 | [0.0065, 0.0258] | 0.6842 |
| SG_DECM12 | 0.0094 | [0.0012, 0.0184] | 0.5263 |
| SG_TARGETM12 | 0.0057 | [-0.0017, 0.0132] | 0.6316 |
| SG_Random | 0.0055 | [-0.0002, 0.0109] | 0.7368 |

这支持一个窄边界结论：在当前 172-gene 交集上，`SG_Static` 的 relation-aligned gain 为正；但 GenePT 与 GO-C 的等 complex 均值更高。

## 4. mask-value contextual 相对 static 的配对增量

| representation | equal-complex difference | 95% CI | sign-flip p |
|---|---:|---:|---:|
| SG_H00 | -0.0058 | [-0.0137, 0.0025] | 0.1883 |
| SG_D12 | -0.0075 | [-0.0154, 0.0005] | 0.0765 |
| SG_D06 | -0.0089 | [-0.0191, 0.0008] | 0.0936 |
| SG_D09 | -0.0091 | [-0.0195, 0.0007] | 0.0913 |
| SG_D03 | -0.0104 | [-0.0205, -0.0009] | 0.0423 |
| SG_H03 | -0.0115 | [-0.0213, -0.0031] | 0.0101 |
| SG_H06 | -0.0127 | [-0.0243, -0.0026] | 0.0230 |
| SG_H12 | -0.0179 | [-0.0283, -0.0083] | 0.0008 |
| SG_H09 | -0.0182 | [-0.0295, -0.0083] | 0.0008 |
| SG_CLSM06 | -0.0185 | [-0.0392, -0.0013] | 0.0585 |
| SG_MEANM12 | -0.0196 | [-0.0395, -0.0034] | 0.0311 |
| SG_DECM12 | -0.0258 | [-0.0464, -0.0095] | 0.0011 |
| SG_TARGETM12 | -0.0296 | [-0.0500, -0.0133] | 0.0001 |

没有任何 contextual 读出的 95% bootstrap CI 完全位于零以上。若只看最接近 static 的两个：

- `SG_H00 - SG_Static`：**-0.0058**，95% CI **[-0.0137, 0.0025]**；
- `SG_D12 - SG_Static`：**-0.0075**，95% CI **[-0.0154, 0.0005]**。

早期与最终 hidden state、CLS/mean/decoder/target mask readout 的等 complex 均值全部为负。若把 0 视为“无增量”，现有数据不支持 contextual promotion。

## 5. zero-value 与 mask-value sensitivity

zero-value 版本同样没有 contextual representation 的 95% CI 完全位于零以上。所有对应 readout 的 `mask-minus-zero` 95% CI 均跨 0。因此 query value 会改变个别 target 的邻域，但没有证据表明它改变总体科学结论。

## 6. 语义表示相对 static

| representation | equal-complex difference versus SG_Static | 95% CI |
|---|---:|---:|
| GenePT_NCBI_UniProt | 0.0248 | [0.0058, 0.0477] |
| GO_C | 0.0158 | [-0.0031, 0.0374] |
| SG_GeneNorm | -0.0013 | [-0.0088, 0.0057] |
| SG_Value | -0.0259 | [-0.0491, -0.0072] |
| SG_Random | -0.0298 | [-0.0534, -0.0108] |

GenePT 相对 `SG_Static` 的等 complex 增量为正；GO-C 点估计为正但区间跨 0。`SG_GeneNorm` 与 static 基本相同，单独 baseline-expression value 不能解释 static 几何。

## 7. 本轮新增的有效性问题

### 7.1 pseudo-state 不是 fold-pure

cell-line consensus 使用全部 perturbation raw pseudobulk；外层测试 target 的原始表达也参与了 median consensus。

### 7.2 response feature selection 不是 fold-pure

`select_response` 在完整 K562/RPE1 response matrices 上按方差选择 response genes。测试 target 的响应参与了 feature selection。

### 7.3 significant-only 改变目标总体

`energy_test_p_value < 0.001` 在建模前筛除扰动。这不是普通缺失过滤，而是按响应显著性选择 target；主分析必须改成 all-finite，significant-only 只能作敏感性。

### 7.4 gene KFold 的正确解释

raw KNN 明确允许同 complex 的伙伴响应进入训练；这是 **within-relation interpolation** 的任务定义，并非无意泄漏。因而：

- 不能把它解释为 unseen-complex generalization；
- 统计推断必须以 complex/relation component 聚类；
- complex-held-out 应用于 learned adapter/generalization 轨道，不能直接替换 partner reference，否则测试 target 将没有关系伙伴可用。

### 7.5 “ceiling / utilization ratio” 命名不成立

CORUM partner mean 是一个 reference predictor，不是数学上界。GenePT 在部分条件可超过它，因此原 `aggregate_utilization_ratio` 不是有界利用率。

### 7.6 其他次要风险

- random embedding 只有 128 维，未与各表示逐一做维度匹配；
- 1199-gene sequence 只代表高表达、在词表内、在 CORUM 中、且被多个 embedding 覆盖的 target；
- layer/readout 多重搜索的 oracle 只能描述异质性；
- CORUM、GenePT 与 GO-C 的知识来源存在重叠，尚未测 semantic residual increment；
- checkpoint 暴露与数据污染未系统登记。

## 8. 复核裁决

**static accessibility：窄支持上通过；contextual-over-static：当前门失败。**

这不是证明 scGPT 的 contextual mechanism 永远无效。它证明的是：在当前 pseudo-state、target-value query、支持集与 within-relation retrieval estimand 下，固定 contextual readout 没有提供超过同一模型 static token table 的可复现增量。

下一轮必须同时改变**评价纯度**与**确认数据**，而不是继续在同一 91 个 target 上挑层。
