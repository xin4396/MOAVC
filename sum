截至 **2026年11月21日**，moavc 已完成“响应标签整理 → MoA environment embedding → 下游 latent DDPM 表达预测”的一整条实验链。**目前最明确的阶段结论是：模型可以拟合训练 cell 的响应关系，但尚未证明能从新 cell 的 control 表达中提取出超过简单基线的、可迁移的 MoA 表示；下游 DDPM 也尚未显示可靠的 MoA 条件收益。**

下面依据本地代码、实验配置、逐轮报告与保存的审计结果整理，没有新增训练，也没有重新登录服务器核验运行状态。现有报告记录的 R11、R12、R13 和 latent DDPM v1 均已完成。

这里把 **harness** 理解为整套实验约定：数据与输入、训练和测试的信息边界、split、候选集、基线、指标、聚合、选模、稳定性和复现条件。

最完整的现成交接入口是 [MOAVC_EXPERIMENT_HANDOFF_2026-11-20.md](F:/CodexWork/moavc/MOAVC_EXPERIMENT_HANDOFF_2026-11-20.md)。虽然文件名是11月20日，但正文已补充11月21日的 DDPM 结果。

---

**一、项目的 rationale：为什么要做 MoA environment embedding**

项目的核心问题是：

> 对一个没有给药数据的新 cell line，能否只根据它的未处理 control 表达，推断它在不同 drug×dose 下的响应环境，并将这种表示用于虚拟细胞的扰动预测？

目前的基本设计是：

\[
z_c=f_\theta(x_c^{control}),\qquad
h_{c,d,\delta}=P_\phi(z_c,e_d,\delta).
\]

其中：

- \(x_c^{control}\)：cell line 的未处理表达。
- \(z_c\)：drug-independent 的基础 cell 表示。
- \(e_d\)：512维 Uni-Mol 药物表示。
- \(\delta\)：剂量。
- \(h_{c,d,\delta}\)：drug×dose 条件下的 cell environment embedding。

这个设计把两个问题分开处理：

1. **control 表达反映 cell 的基础状态。**  
   新 cell 推理时能拿到 control，因此基础 encoder 必须只依赖 control。

2. **cell 之间的响应关系可以随 drug×dose 改变。**  
   两个 cell 对某种药相似，对另一种药可能不同。因此不能只学一份全局 cell similarity，也不能把不同药物的上下调沿 drug 维度平均掉。

训练时利用已知 cell 的处理响应构造监督：同一种 drug×dose 下，比较两个 cell 的转录响应向量，得到它们的条件相似度。模型则从 control 出发，生成能够重现这些关系的 \(h\)。

**这里的“MoA”是研究目标，实际监督仍是转录响应相似性，并不是直接测量的机制标签。** 因此，拟合 similarity 并不能自动证明模型学到了药物作用机制。

当前路线尤其强调：下游要拿到的是**可复用的 embedding \(h\)**。如果一个模型必须把 query 和 candidate 一起输入 pair MLP 才能得到分数，即使分数不错，也不能直接说明单独导出的 \(h\) 有用。这是 R11 之后转向固定几何监督的重要原因。

后续已经实际开展下游验证：

\[
\text{control expression}
\rightarrow \text{control latent}
\rightarrow \text{conditional DDPM}
\rightarrow \text{treated-state latent}
\rightarrow \text{treated expression}.
\]

两组 DDPM 分别使用冻结的 MoA \(h\)，或者直接使用 Uni-Mol+drug dose，检验 MoA 表示是否带来额外价值。

当前设计共识见 [CURRENT_DESIGN_DECISIONS.md](F:/CodexWork/moavc/cellmoa_v1/CURRENT_DESIGN_DECISIONS.md)。

---

**二、数据基础、预处理与标签版本**

**1. 当前建模使用的数据规模**

| 项目 | 当前设置 |
|---|---|
| 数据来源 | Tahoe；当前模型没有联合训练 L1000 |
| Cell lines | 50 |
| Cell split | 35 train / 7 validation / 8 test |
| 药物 | 377 |
| 剂量 | 0.05、0.5、5 μM |
| Canonical 条件 | 377×3＝1,131 |
| 基因 | 62,710 |
| DDPM 有效 cell×drug×dose 单位 | 54,180 |
| DDPM train / val / test 单位 | 37,287 / 7,887 / 9,006 |

需要区分三个覆盖数字：

- Tahoe cell-based 数据管线观察到 **56,827** 个药物条件。
- 其中 **54,470** 个满足细胞数要求、可检验。
- 模型限定到377个 canonical drugs 后，DDPM 实际使用 **54,180** 个有效单位。

另外，1,131种条件和大量 pair 并不等于有大量独立 cell 背景。control→environment 的跨 cell 外推，主要仍依赖 **35个训练 cell**。

有些 cell 的有效条件非常少：一个训练 cell 只有2个有效列表，另一个约210个，大多数 cell 则有约1,126–1,128个。采用 cell 等权宏平均时，这种覆盖差异需要一起解释。

固定划分见 [split.json](F:/CodexWork/moavc/cellmoa_v1/artifacts/split.json)。8个 test cell 为：

`CVCL_1125、CVCL_0179、CVCL_0366、CVCL_0371、CVCL_1055、CVCL_1717、CVCL_1724、CVCL_1098`。

**2. 三种表达处理方式必须分清**

| 用途 | 实际处理 |
|---|---|
| MoA encoder 的 control 输入 | 同 block 汇总原始 counts → 按全基因总量 CP10k → log1p → 合格 block 等权平均 |
| 响应 LFC / DEG | 先逐细胞 CP10k → 每板 treated/control 分别取细胞均值 → 匹配板等权 → 计算效应 |
| 下游 DDPM 输入与目标 | 读取上述配对 treated/control mean CP10k，再取自然 log1p |

因此，MoA encoder 与 DDPM 的 control 聚合方式存在已知差异。这个差异需要披露，但目前没有证据证明它就是失败的主要原因。

MoA control 训练还使用 block bootstrap 构造增强 views。全基因版本采用 train control 均值中心化，不逐基因 z-score、不额外 upper clip，也不重复 log。

注意：[preprocessing.json](F:/CodexWork/moavc/cellmoa_v1/artifacts/preprocessing.json) 中仍保留初版“4096输入基因”的信息，不能直接作为 R11–R13 的配置说明；后期全基因模型应以对应 checkpoint 和 `plan.json` 为准。

**3. 当前 Tahoe DEG 定义**

当前有效版本是 `cell_based_v1`：

- 原始 UMI，使用 `pass_filter=full` 的细胞。
- treated 与 DMSO 按 plate、cell line、time、vehicle 匹配。
- 每个有效板两组均至少30个细胞。
- 板内 Wilcoxon/Mann–Whitney U，跨板 van Elteren。
- 每条件对全部62,710基因做 BH。
- 主 DEG：`padj<0.05 & abs(effect)>0.5`。
- 不设独立 well 数门槛。
- 缺失时间保留 unknown，不推断为24H。

效应定义为：

\[
LFC_g=\log_2\frac{\mu^{treated}_g+10^{-9}}
{\mu^{DMSO}_g+10^{-9}}.
\]

旧的 Tahoe well-level edgeR DEG 已按历史要求删除；pseudobulk 表达保留。L1000 继续使用 RNA well 层面的 limma，不能套用 Tahoe 的单细胞检验。

方法及结果：

- [cell_based_v1/METHOD.md](F:/CodexWork/moavc/chat_version/scripts/cell_based_v1/METHOD.md)
- [cell_based_v1/RESULTS.md](F:/CodexWork/moavc/chat_version/tahoe/cell_based_v1/reports/RESULTS.md)
- [DEG_VERSION_INDEX.json](F:/CodexWork/moavc/chat_version/DEG_VERSION_INDEX.json)

**4. Similarity 标签经过四个主要版本**

定义后期使用的固定 kernel：

\[
K(a,b)=\frac{2a^\top b}{\|a\|^2+\|b\|^2}.
\]

它同时受方向和相对幅度影响。例如同方向、幅度差2倍时，分数为0.8。

| 标签 | 响应特征 | Similarity | 使用阶段 |
|---|---|---|---|
| A | DEG支持上的 LFC，其余置0 | signed min/max overlap | 初版、v2/v3、R4 |
| B | DEG支持上的 asinh(LFC) | 固定 \(K\)，另试 z-score/CDF | R5–R9 |
| C | 全部基因 raw LFC，不按显著性筛选 | 固定 \(K\) | R10 |
| D | 全部基因 asinh(LFC/1)，不筛显著性 | 固定 \(K\) | R11–R13 |

**不同标签版本的分数不能直接横向排名。** R10→R11 既改标签又改训练，不能把指标差异全部归因于模型结构。

当前 D 标签构建见 [round11/build_labels.py](F:/CodexWork/moavc/cellmoa_v1/round11/build_labels.py)，规范见 [label_spec.json](F:/CodexWork/moavc/cellmoa_v1/artifacts/round11_tokens_asinh/label_spec.json)。

---

**三、已经做过的实验：研究问题、结果与阶段含义**

**1. 数据审核与管线重建**

最早的工作是审查 Tahoe/L1000 的 DEG、跨 cell overlap 和跨系列复现性分析。

审核发现的主要问题包括：

- 跨 L1000 系列计算相关性时没有按 gene ID 对齐。
- 单处理孔的零方差被 floor 成极小值，造成虚假精度。
- 部分聚合遗漏 dose，或混合 time。
- 不同代码使用不同 overlap 定义，导致结果不可直接比较。

例如，原来全部3,929个共同条件跨系列相关性均值约0.0036；按 gene ID 对齐后约0.1249。修正后仍不算高复现，但原来的“几乎零相关”论据失效。

随后建立了 `chat_version` 的数据与统计管线，整理表达、DEG、匹配条件、不可检验状态和验证报告，再切换到当前 Tahoe cell-based 版本。

这部分的意义是：**先把数据定义和统计单位理清，避免模型学习管线错误。** 旧分析报告应作为历史记录阅读。

证据入口：[审核报告.md](F:/CodexWork/moavc/deg_audit/审核报告.md)。

**2. MoA embedding 历史实验总览**

以下数字保留各轮原有任务口径，不能跨标签作统一排行榜。

| 阶段 | 实验目的与主要改动 | 主要结果 | 阶段结论 |
|---|---|---|---|
| 初版 | 4096基因、256维 cell encoder、drug-conditioned mask | train-reference 重评 test Spearman约0.373；去mask后接近不变 | 条件模块作用弱 |
| optimization_v2 | 排序目标、mask/projection 对比 | test rho约0.586，但 MAE约0.580 | 排名与数值校准严重脱节 |
| optimization_v3 | query-blind prior＋embedding residual | test rho约0.680；MAE0.1452，高于 prior0.0191 | 收益大量来自 prior，不能归为纯 embedding |
| R4 | 全62,710基因、低秩输入、条件投影 | A标签 test rho0.530；列表宏R²−37.12 | 排序有信号，校准不足 |
| R5 | B标签、逐步降维全连接、多个读出与变换 | test rho0.526；标准宏R²−117.05 | 换标签和加大 encoder 未解决外推 |
| R5局部拟合 | one-hot/free code/旧低秩投影/完整矩阵对照 | 48条件完整矩阵版本 train R²约0.9992–0.9999 | 小面板可拟合，旧 projection 容量是局部限制之一 |
| R6 | 完整矩阵条件投影、冻结/慢学习率/anchor/算子和读出约束 | 保留分支 test rho0.532；宏R²约−105.99 | 通过当轮工程检查，未证明可用于VC |
| R7 | 48→96→192→384条件受控扩展 | 加宽后384条件 train R²约0.9883 | 训练暴露和容量会影响拟合；只有train诊断 |
| R8 | 全377×3面板，plain/robust/anchored | 某分支最高train R²0.977；选中anchored test rho0.578、宏R²−51.37 | 全量拟合与跨cell泛化之间仍有明显差距 |
| R9 | 不训练，检查极端负R²和低方差列表 | 确认稀疏DEG支持导致近零方差；增加epsilon R²诊断 | 极端负R²不一定是实现错误 |
| R10 | 全基因raw LFC；scalar与固定kernel比较 | scalar test rho0.471；发现零表达边界支配案例 | 去掉DEG阈值并不能自动获得可靠响应特征 |
| R11 | 同D标签重训旧单向量模型与多token attention | legacy0.570；tokens0.371；prior0.633 | 多token方案失败；旧完整头也未超过prior |
| R12 | 取消pair MLP/calibrator，固定 \(K\) 直接监督 raw \(h\)；joint/slow encoder | test rho0.370–0.411，均未过稳定性门槛 | drug信息进入几何，但新cell校准与泛化差 |
| R13 | 两个seed，普通encoder vs谱范数约束 | 谱约束两个seed均降低test rho | 本轮谱约束未改善泛化 |

**3. R11：分清“分数预测好”和“embedding 好”**

R11 使用相同 D 标签、split、loss、optimizer 和3000次全条件更新，对比原单向量架构与多 token attention。

| 方法 | Train rho | Val rho | Test rho | Test PCC | Test Top5 |
|---|---:|---:|---:|---:|---:|
| Query-blind prior | 0.6150 | 0.7078 | **0.6332** | **0.6895** | **0.4615** |
| Legacy 完整输出头 | 0.7785 | 0.5395 | 0.5697 | 0.6459 | 0.4038 |
| 多token attention | 0.5194 | 0.3849 | 0.3711 | 0.5311 | 0.2511 |

重要诊断：

- token 模型打乱 drug 后，test rho从0.37111变为0.37038，几乎不变。
- token 模型跨 cell 的条件 \(h\) 余弦中位数约0.9857，输出过度相似。
- token 模型 test MAE约0.90075，校准很差。
- legacy 完整头 test rho0.570，但直接读取其 \(h\) 几何只有约0.343。

因此，**legacy 的部分能力存在于输出头，不能全部归给 \(h\)**；token bank 有一定多样性，也不代表条件读取成功。

报告：[R11 RESULTS.md](F:/CodexWork/moavc/cellmoa_v1/artifacts/round11_tokens_asinh/RESULTS.md)。

**4. R12：直接监督可复用的 \(h\)**

R12 改成：

- 全基因低秩输入：`62710→128→4096→1024→256→64`。
- drug×dose 生成完整64×64条件变换。
- 输出 raw \(h\)，保留范数。
- 固定 \(K(h_q,h_j)\) 直接预测 similarity。
- 无 pair MLP、无 scalar calibrator、无 teacher。
- 比较 encoder 联合更新与后期慢更新。

| 方法 | Train rho | Val rho | Test rho | Test PCC | Test MAE |
|---|---:|---:|---:|---:|---:|
| Prior | 0.6150 | 0.7078 | **0.6332** | 0.6895 | 0.01324 |
| R12 joint | 0.7443 | 0.4324 | 0.4106 | 0.2215 | 0.10667 |
| R12 slow，主seed | 0.7456 | 0.4086 | 0.3882 | 0.2243 | 0.10168 |
| R12 slow，重复seed | 0.7335 | 0.3956 | 0.3703 | 0.1950 | 0.10064 |

主seed test 的误差不只是一个整体 offset：

- GT列表均值0.08165，pred均值0.17814。
- GT列表内SD均值0.01643，pred约0.13462。
- MSE中约38%来自列表均值偏移，约62%来自去均值后的幅度/形状误差。
- 正确drug→打乱drug，test rho从0.3882降到0.2757，说明条件确实进入了几何。
- 但减去 candidate prior 后的残差rho由train0.5406降至test0.1517，query-specific迁移仍弱。

三组全部未通过固定 validation bootstrap 稳定性门槛。后期降低 encoder LR，在这次受控比较里没有改善结果。

报告：[R12 RESULTS.md](F:/CodexWork/moavc/cellmoa_v1/artifacts/round12_direct_kernel/RESULTS.md)。

**5. R13：谱约束是否能改善 control encoder 外推**

这轮保留 R12 的固定几何路线，两版共同加入两个残差块，配对比较是否启用 encoder 线性层谱范数上限：

\[
cap=1.05\times\text{initial estimated spectral norm}.
\]

两个配对seed，各3000次全条件更新；同参数量、初始化、dropout随机性、增强、loss和预算。

| Seed | Plain test rho | Spectral test rho | Spectral−Plain |
|---|---:|---:|---:|
| 20261120 | 0.3833 | 0.3620 | −0.0213 |
| 20261121 | 0.3925 | 0.3782 | −0.0144 |

谱约束确实激活，但完整 encoder 的 control→\(z\) 敏感性没有一致下降，四组仍全部未过稳定性门槛。

该实验只支持：**本轮近似线性层谱约束无效。** 它不是完整 SNGP，也没有证明全局距离保持，不能据此否定所有谱约束或所有泛化方法。

报告：[R13 RESULTS.md](F:/CodexWork/moavc/cellmoa_v1/artifacts/round13_spectral_encoder/RESULTS.md)。

---

**四、下游 latent DDPM：已经实际做了什么**

这部分已经完成，早期文档中“尚未搭建下游VC”的描述已过时。

**1. 比较问题**

两组都接收同一个 control latent，额外条件分别为：

- **Drug组**：512维 Uni-Mol＋1维标准化 log-dose。
- **MoA组**：冻结的64维 \(h\)。

MoA \(h\) 来自 **R12 joint** 的既有 validation-selected checkpoint，**不是 R13 模型，也不是 `latest_model.json` 指向的 R6**。

本轮扩散的目标是 **treated-state latent 本身**，不是 LFC，也不是 latent delta。

**2. 共享重建模型**

先做了 attention AE：

- 训练4000步。
- validation 重建MSE为 no-change 的5.06倍。
- 判为失败 pilot，保留产物。

之后做 train-only 线性 AE/SVD 容量对照：

- 测试128到8192维。
- validation选择8192维。
- 全部62,710基因参与。
- 8192坐标分成64×128 tokens。
- 这只是计算分组，不是生物学模块。
- E/D由两组共享并冻结。

AE oracle，即 \(D(E(\text{真实treated}))\)，完整test delta PCC约0.6718，说明重建空间中仍有两组DDPM没有学到的响应信号。但它使用真实 treated，只能当诊断。

**3. DDPM 设置**

| 项目 | 实际设置 |
|---|---|
| 去噪主干 | 4层、4 head、width128 Transformer |
| Control输入 | 逐token融合＋cross-attention memory |
| Forward diffusion | 200步 cosine schedule |
| 目标 | v-prediction MSE |
| 训练 | 每组6000 updates，batch256 |
| 抽样 | 先均匀选训练cell，再选该cell可用条件 |
| Optimizer | AdamW，lr2e−4 |
| EMA | warm-start EMA |
| 参数量 | Drug 1,198,976；MoA 1,141,504 |
| 选模 | 最低固定噪声、cell-macro validation v-MSE |
| 最终step | Drug3500；MoA3000 |
| 推理 | 50步DDIM，eta0，4次采样，guidance1 |
| 多次采样聚合 | 每次decode到完整state，在线性表达上平均后再log1p |
| 冻结模块 | 共享E/D与MoA F均冻结 |

两组主干初始化一致，adapter因输入维度不同而参数量不同。MoA还使用了历史监督预训练，因此只能说**本轮DDPM训练预算匹配**，不能说总探索或总预训练成本相同。

此前还有一组逐PC whitening、仅cross-attention读取control的两臂pilot，各6000步。验证改善很小，未进入test；后续两组共同改为全局latent scale、逐token control融合与warm EMA。这些共同修改没有逐项隔离归因。

**4. 完整面板结果**

这里 delta 指 matched log-expression change。

| 方法 | Train delta PCC | Val delta PCC | Test delta PCC | Test表达MAE | Test delta R² |
|---|---:|---:|---:|---:|---:|
| Drug-dose平均响应 | **0.6198** | **0.6511** | **0.6172** | **0.01100** | **0.3834** |
| DDPM Drug | 0.4247 | 0.2140 | 0.1925 | 0.03238 | −7.8508 |
| DDPM MoA | 0.4040 | 0.2122 | 0.1945 | 0.03179 | −7.6194 |
| AE oracle，GT诊断 | 0.9810 | 0.7032 | 0.6718 | 0.01490 | −0.1660 |
| No-change | 未定义 | 未定义 | 未定义 | 0.01294 | −0.0017 |

MoA−Drug 的test delta PCC差值为：

\[
0.001945,\quad 95\%\text{ paired-cell CI }
[-0.000511,\ 0.004234].
\]

只有一个训练seed，区间包含0，**没有可靠MoA增益证据**。两组不仅低于平均响应基线，表达MAE也高于 no-change。

而且 train 本身就低于均值基线，因此失败不只是 unseen-cell generalization，也涉及生成转移拟合和条件读取。

增加采样的validation诊断中，4→64次采样仅把PCC从：

- Drug：0.2086→0.2141。
- MoA：0.2062→0.2110。

这不足以解释或填补基线差距。

报告：[latent DDPM RESULTS.md](F:/CodexWork/moavc/cellmoa_v1/artifacts/latent_ddpm_v1/RESULTS.md)。

---

**五、harness 条件：如何保证结果有明确含义**

**A. MoA similarity harness**

| 条件 | 当前约定 |
|---|---|
| 泛化对象 | 按cell身份留出的新cell；没有该cell的treated训练输入 |
| Split | 固定35/7/8，历史split |
| 一个评估unit | query cell×drug×dose |
| Candidates | 始终使用可用train cells |
| Train query | 排除自身 |
| 有效性 | 同步按GT可用mask筛选，至少4个candidate |
| Missing | 不可检验或缺失保持missing，不能填0响应 |
| 条件 | 保留drug×dose，不沿drug平均响应 |
| 当前标签 | 全基因asinh(LFC)固定kernel；R11–R13共用 |
| 主指标 | Spearman；Top1/Top5邻居一致性补充 |
| 辅助指标 | Pearson、MAE、MSE、标准R²、epsilon R² |
| 聚合 | 先对query内条件平均，再对query等权平均 |
| Robustness | query control换bootstrap views1/5/9，reference保持clean，测score漂移 |
| 负对照 | Drug shuffle保持dose；query-control shuffle；mean-control |
| 复现 | 输入/标签hash、plan、源码快照、history、selection lock、逐列表预测 |

Spearman按固定candidate ID对应关系计算，不能分别排序GT和pred后再求相关。

标准列表R²：

\[
R^2=1-\frac{MSE}{Var(GT)}.
\]

它没有负下界；稳定化版本为：

\[
R^2_\epsilon=1-\frac{MSE}{Var(GT)+10^{-4}}.
\]

两者是不同指标。epsilon不改变模型或标签，也不能作为能力改善。

Top5处理GT边界并列，避免任意打破ties。各指标的有效列表可能不同，所以必须同时看覆盖数。

实现主要在：

- [round4/metrics.py](F:/CodexWork/moavc/cellmoa_v1/round4/metrics.py)
- [round9/metrics.py](F:/CodexWork/moavc/cellmoa_v1/round9/metrics.py)
- [round11/metrics.py](F:/CodexWork/moavc/cellmoa_v1/round11/metrics.py)

**最重要的基线：query-blind candidate prior**

固定条件 \(p\)，candidate \(j\) 的预测分数是它与其他训练cell真实响应similarity的均值：

\[
\hat S_{q,j,p}^{prior}
=\operatorname{mean}_{i\in Train\setminus\{j,q\}} S_{i,j,p}.
\]

它不看新query的control，但使用训练响应参考库，并区分drug×dose。train评估时还排除query对prior的贡献。

这个基线强，说明 candidate 共性排序包含大量信息。模型需要展示相对于它的query-specific增量；在预测上加prior所得的收益不能全部归给embedding。

实现：[round11/audit.py](F:/CodexWork/moavc/cellmoa_v1/round11/audit.py)。

**R11–R13训练与选模条件**

- R11–R13正式run均为3000次全1,131条件更新；chunk48做梯度累积。
- AdamW：lr3e−4，cosine降到1e−5，weight decay1e−4，clip10。
- clean与control-bootstrap交替。
- 每250步保存检查点。
- R11 loss：0.70 soft listwise KL＋0.25 Pearson项＋0.05标准化Huber。
- R12/R13 loss：0.50连续pairwise BCE＋0.25 Pearson项＋0.25标准化Huber，query宏聚合。
- R12/R13按最大validation rho且bootstrap drift≤**0.002847219**选模。
- 无checkpoint通过时，按预声明fallback报告val rho最高者，标记`stable=false`。
- 在本轮test前锁定权重，不根据test挑seed或重新校准。

基准配置可查：

- [R11 benchmark_contract.json](F:/CodexWork/moavc/cellmoa_v1/artifacts/round11_tokens_asinh/benchmark_contract.json)
- [R12 benchmark_contract.json](F:/CodexWork/moavc/cellmoa_v1/artifacts/round12_direct_kernel/benchmark_contract.json)
- [R13 benchmark_contract.json](F:/CodexWork/moavc/cellmoa_v1/artifacts/round13_spectral_encoder/benchmark_contract.json)

**B. 下游VC harness**

完整自定义评估：

- 全62,710基因、全部54,180有效单位。
- 每unit计算表达误差、delta PCC/Spearman、未clip的标准R²和top-effect overlap。
- 条件先在cell内平均，再对cell等权平均。
- 未定义相关性保留并统计，不任意记0。
- 不用全表达的高相关性代替药物变化预测能力。

基线包括：

1. No-change。
2. 固定drug×dose，训练cell平均log-expression变化，转移到query control。
3. Control reconstruction。
4. AE oracle，仅GT诊断。

第二个基线在train评估时排除query。它与similarity任务的query-blind prior是两种不同对象。

**官方 Cell-Eval 部分**

- Cell-Eval **0.8.2**。
- 独立Python3.12环境。
- seed9127固定选32 drugs×3 doses。
- 全基因。
- 仅使用适用于pseudobulk均值表达的指标。
- 不伪造单细胞，不计算缺乏依据的DE p值或单细胞分布指标。
- 覆盖：train34cell/3178unit，val7/672，test8/767。
- 一个稀疏训练cell在官方面板无有效条件；完整自定义评估仍保留35个train cells。
- 这是适配的均值表达benchmark，不能称为完整STATE单细胞benchmark复现。

保留了两种control口径：

| 版本 | Control参考 | 解读 |
|---|---|---|
| Native官方 | 每cell统一control | No-change也可有非零delta PCC，混入参考差异 |
| Matched官方 | 每条件matched control，传入gene-effect向量 | 更贴合本项目配对任务 |

Native no-change 的test delta PCC达0.3814、PDS约0.8731，不能解释为预测到了药物响应。

Matched官方test结果：

| 方法 | Delta PCC | MAE | Cosine PDS |
|---|---:|---:|---:|
| Drug-dose平均响应 | **0.6307** | **0.01060** | **0.9964** |
| DDPM Drug | 0.1884 | 0.03224 | 0.7164 |
| DDPM MoA | 0.1906 | 0.03167 | 0.7136 |
| No-change | 未定义 | 0.01249 | 0.5052 |

Matched官方与同面板自定义MAE/PCC已数值核验。因此目前负结果并非仅由官方control参考差异造成。

相关配置和审计：

- [evaluation/protocol.json](F:/CodexWork/moavc/cellmoa_v1/artifacts/latent_ddpm_v1/evaluation/protocol.json)
- [comparison_audit.json](F:/CodexWork/moavc/cellmoa_v1/artifacts/latent_ddpm_v1/comparison_audit.json)
- [evaluation_audit.json](F:/CodexWork/moavc/cellmoa_v1/artifacts/latent_ddpm_v1/evaluation_audit.json)

**C. 结论适用范围**

当前所有近期结果都基于同一反复查看的历史split，属于**回顾性固定split比较**，不是全新外部盲测。8个test cell限制了统计结论；bootstrap重采样的是cell，不应把大量条件或两个seed当成更多独立cell。

Control bootstrap衡量采样稳定性，不能替代跨实验、跨平台或独立给药重复验证。

---

**六、实验脚本和结果位置索引**

本地根目录：[F:/CodexWork/moavc](F:/CodexWork/moavc)。

服务器对应根目录：

```text
/home/jiaxin/jiaxin/moavc/
├── chat_version/       数据、表达、DEG
└── cellmoa_v1/         建模、训练与实验产物
```

本地保存代码、报告、图、预测和部分checkpoint；原始大数据、完整历史checkpoint主要在服务器。

**数据与统计管线**

| 用途 | 脚本 | 报告/结果入口 |
|---|---|---|
| 旧管线审核 | [deg_audit](F:/CodexWork/moavc/deg_audit) 下三个 `probe_*.py` | [审核报告.md](F:/CodexWork/moavc/deg_audit/审核报告.md) |
| 表达/DEG主流水线 | [run_pipeline.py](F:/CodexWork/moavc/chat_version/scripts/run_pipeline.py)；`01`–`13`阶段脚本 | [reports/RESULTS.md](F:/CodexWork/moavc/chat_version/reports/RESULTS.md)，其中Tahoe DEG为历史版本 |
| 当前Tahoe cell-based | [run.py](F:/CodexWork/moavc/chat_version/scripts/cell_based_v1/run.py)、`core.py`、`validate_outputs.py` | [当前RESULTS.md](F:/CodexWork/moavc/chat_version/tahoe/cell_based_v1/reports/RESULTS.md) |
| Tahoe细胞数统计 | [14_summarize_tahoe_cell_numbers.py](F:/CodexWork/moavc/chat_version/scripts/14_summarize_tahoe_cell_numbers.py) | [tahoe_cell_numbers](F:/CodexWork/moavc/chat_version/reports/tahoe_cell_numbers) |
| L1000零DEG审计 | [15_audit_l1000_zero_deg.py](F:/CodexWork/moavc/chat_version/scripts/15_audit_l1000_zero_deg.py) | [REPORT.md](F:/CodexWork/moavc/chat_version/reports/l1000_zero_deg_audit/REPORT.md) |
| L1000连续响应导出 | [16_export_l1000_unfiltered_lfc.py](F:/CodexWork/moavc/chat_version/scripts/16_export_l1000_unfiltered_lfc.py) | [README.md](F:/CodexWork/moavc/chat_version/reports/l1000_lfc_unfiltered/README.md) |

服务器关键数据：

```text
chat_version/tahoe/pseudobulk_counts/plate*.h5
chat_version/tahoe/expression/<cell>/
chat_version/tahoe/cell_based_v1/deg/<cell>/
  all_results.h5
  condition_status.parquet
  significant.parquet
  genes.parquet

chat_version/l1000/<series>/expression/<cell>/
chat_version/l1000/<series>/deg/<cell>/
chat_version/l1000/<series>/lfc_unfiltered/<cell>/lfc.h5
```

`all_results.h5`中的核心字段包括：`gene_id、condition_id、effect、pvalue、padj、drug_mean_cp10k、dmso_mean_cp10k、pct_drug、pct_dmso`。读取必须按ID对齐。

**MoA模型实验**

| 实验 | 主要脚本 | 报告/产物 |
|---|---|---|
| 初版 | [prepare.py](F:/CodexWork/moavc/cellmoa_v1/prepare.py)、[train.py](F:/CodexWork/moavc/cellmoa_v1/train.py)、`model.py` | [MODEL_CARD.md](F:/CodexWork/moavc/cellmoa_v1/artifacts/MODEL_CARD.md) |
| 训练参考库重评 | [evaluate_anchor.py](F:/CodexWork/moavc/cellmoa_v1/evaluate_anchor.py) | [ANCHOR_RANKING.md](F:/CodexWork/moavc/cellmoa_v1/artifacts/ANCHOR_RANKING.md) |
| v2 / v3 | [optimize_v2.py](F:/CodexWork/moavc/cellmoa_v1/optimize_v2.py)、[optimize_v3.py](F:/CodexWork/moavc/cellmoa_v1/optimize_v3.py) | [optimization_v2](F:/CodexWork/moavc/cellmoa_v1/artifacts/optimization_v2)、[optimization_v3](F:/CodexWork/moavc/cellmoa_v1/artifacts/optimization_v3) |
| 数值列表审计 | [evaluate_similarity_lists.py](F:/CodexWork/moavc/cellmoa_v1/evaluate_similarity_lists.py) | [similarity_list_evaluation](F:/CodexWork/moavc/cellmoa_v1/artifacts/similarity_list_evaluation) |
| R4 | [round4](F:/CodexWork/moavc/cellmoa_v1/round4) | [round4_allgenes](F:/CodexWork/moavc/cellmoa_v1/artifacts/round4_allgenes) |
| R5 | [round5/train.py](F:/CodexWork/moavc/cellmoa_v1/round5/train.py) | [round5_asinh_progressive](F:/CodexWork/moavc/cellmoa_v1/artifacts/round5_asinh_progressive) |
| 小面板/投影容量 | R5的 `small_fit.py、grand_ablation.py、onehot_projection_test.py` | R5目录下 `small_fit_tests/、onehot_projection_tests/` |
| R6 | [round6](F:/CodexWork/moavc/cellmoa_v1/round6) | [round6_r2_robust](F:/CodexWork/moavc/cellmoa_v1/artifacts/round6_r2_robust) |
| R7 | [round7/train.py](F:/CodexWork/moavc/cellmoa_v1/round7/train.py)、`capacity.py` | [round7_controlled_scale](F:/CodexWork/moavc/cellmoa_v1/artifacts/round7_controlled_scale) |
| R8 | [round8](F:/CodexWork/moavc/cellmoa_v1/round8) | [round8_full_panel](F:/CodexWork/moavc/cellmoa_v1/artifacts/round8_full_panel) |
| R9 | [diagnose_variance.py](F:/CodexWork/moavc/cellmoa_v1/round9/diagnose_variance.py)、`report_epsilon.py` | [round9_variance_audit](F:/CodexWork/moavc/cellmoa_v1/artifacts/round9_variance_audit) |
| R10 | [round10](F:/CodexWork/moavc/cellmoa_v1/round10)，含标签、来源贡献与药物审计 | [round10_full_lfc](F:/CodexWork/moavc/cellmoa_v1/artifacts/round10_full_lfc) |
| R11 | [round11](F:/CodexWork/moavc/cellmoa_v1/round11)：`build_labels/train/model/loss/audit/infer.py` | [round11_tokens_asinh](F:/CodexWork/moavc/cellmoa_v1/artifacts/round11_tokens_asinh) |
| R12 | [round12](F:/CodexWork/moavc/cellmoa_v1/round12)：`train/model/loss/audit/diagnose/infer.py` | [round12_direct_kernel](F:/CodexWork/moavc/cellmoa_v1/artifacts/round12_direct_kernel) |
| R13 | [round13](F:/CodexWork/moavc/cellmoa_v1/round13)：`train/model/audit/infer/status_snapshot.py` | [round13_spectral_encoder](F:/CodexWork/moavc/cellmoa_v1/artifacts/round13_spectral_encoder) |

模型数据缓存主要位于服务器：

```text
cellmoa_v1/artifacts/
  dataset.npz
  controls_full.npz
  drug_embeddings.npz
  conditions.csv
  split.json
  label_asinh_inner_v1/labels.npz
  round10_full_lfc/full_lfc.npy
  round10_full_lfc/labels.npz
  round11_tokens_asinh/labels.npz
```

其中 `full_lfc.npy` 形状为 `[50,1131,62710]`，float32约14.185GB。

R11–R13产物中，复核时优先看：

- `plan.json`、`benchmark_contract.json`。
- `history.jsonl`、`result.json`。
- `audit/selection_lock.json`。
- `summary.csv`、`per_cell_summary.csv`。
- `audit/results.json`。
- 各split逐列表GT/pred。
- `cell_representation.npz` 和 `conditional_environment.npz`。
- `inference_validation.json`。
- 训练曲线、parity、benchmark图。

**Latent DDPM脚本**

代码目录：[latent_ddpm_v1](F:/CodexWork/moavc/cellmoa_v1/latent_ddpm_v1)。

| 功能 | 脚本 |
|---|---|
| 构造配对表达数据 | `prepare.py` |
| Attention AE pilot | `train_autoencoder.py` |
| 线性重建容量评估与最终冻结 | `probe_linear_reconstruction.py`、`finalize_linear.py` |
| 生成MoA/Drug条件 | `conditions.py` |
| DDPM模型与训练 | `model.py`、`train_diffusion.py` |
| 推理与完整评估 | `infer.py`、`evaluate.py`、`metrics.py` |
| 官方评估 | `official_smoke.py`、`official_eval.py` |
| 比较与数值审计 | `audit_comparison.py`、`audit_evaluation.py` |
| 条件使用与采样诊断 | `diagnose_conditioning.py`、`sampling_audit.py` |
| 运行编排与状态 | `worker.py`、`launch.py`、`status.py`、`finish_evaluation.py` |
| 汇总导出 | `report.py`、`export_reports.py`、`validate_delivery.py` |

结果根目录：[artifacts/latent_ddpm_v1](F:/CodexWork/moavc/cellmoa_v1/artifacts/latent_ddpm_v1)。

主要内容：

```text
data/manifest.json
autoencoder/
conditions/
drug/plan.json、selection_lock.json、history.jsonl
moa/plan.json、selection_lock.json、history.jsonl
attention_autoencoder_pilot/
whitened_diffusion_pilot/

comparison_audit.json
evaluation_audit.json
validation_sampling_audit.json
RESULTS.md

evaluation/
  summary.csv
  per_cell.csv
  drug_units.parquet
  moa_units.parquet
  baseline_units.parquet
  paired_comparison.json
  official_panel.json
  official_cell_eval/
  official_cell_eval_paired/
  learning_curves.png
  test_cells.png
  gt_pred_changes.png
```

直接查看数字最方便的是 [summary.csv](F:/CodexWork/moavc/cellmoa_v1/artifacts/latent_ddpm_v1/evaluation/summary.csv) 和 [per_cell.csv](F:/CodexWork/moavc/cellmoa_v1/artifacts/latent_ddpm_v1/evaluation/per_cell.csv)。

训练环境记录为PyTorch `2.11.0+cu128`；服务器Python为：

```text
/home/jiaxin/jiaxin/moavc/cellmoa_v1/.venv/bin/python
```

官方评估使用该实验目录下的独立Python3.12环境。

另外，本地 [latent_ddpm_mlp_v2](F:/CodexWork/moavc/cellmoa_v1/latent_ddpm_mlp_v2) 目前只看到 `GeneDrug_paper_text.txt`，没有实现和训练结果，不能计为已完成的MLP v2实验。

---

**七、阶段性的 conclusion，以及还没有解决的问题**

**1. 已经建立了可复核的实验链，但还没有建立有效的可迁移MoA表示。**

项目目前的成果主要是：

- 整理了数据定义和条件匹配。
- 明确了不同标签的统计含义。
- 建立了train-reference ranking harness。
- 引入了强prior与条件/query负对照。
- 区分了score能力与embedding几何。
- 完成了真实下游表达预测比较。

这些使失败可以被定位和复核，但不等于模型已经成功。

**2. 小面板可拟合，不能推出跨cell可预测。**

完整条件矩阵、加宽和增加训练暴露确实能改善train拟合。R5/R7说明标签在局部可以用几何表示，旧projection存在容量限制。

但R8、R12、R13说明，即使train很好，新cell依然明显偏移。因而“继续加容量”还没有成为得到证据支持的泛化解决方案。

**3. Query-blind prior是当前相似度任务最重要的参照。**

同一最新标签上prior test rho0.633，高于R11 legacy以及R12/R13。

这说明candidate共性很强，但其来源还没分解清楚，可能涉及共享响应程序、响应范数、检出率、噪声、共用control或批次。不能据此直接断言query没有信息，也不能断言标签全部无效。

当前模型尚未展示稳定的、超越这些共性的query-specific迁移能力。

**4. 标签有两种不同风险，连续压缩没有自动解决它们。**

- **DEG硬阈值风险**：稀疏支持导致大量exact-zero similarity和近零列表方差，使标准R²极端负。
- **Raw LFC风险**：`1e−9` pseudocount在单侧零表达处产生很大LFC。

两个已审计案例中，单侧零表达基因贡献了约99.6%的raw LFC平方幅度；asinh后仍约91%。这是两个案例的证据，尚不能外推为全数据比例。

Asinh压缩幅度，但不估计效应可靠性。R11–R13仍继承这一问题。

**5. 细胞层显著性不等于独立给药可重复性。**

当前cell-based管线在A549的14个板内DMSO孔对孔诊断中，主DEG数中位数117.5，范围32–741，说明well间背景差异也能被细胞层检验检测到。

这与旧管线“少量DMSO拆分诊断0个DEG”属于不同统计设计，不能拿旧结果为当前全部标签背书。

**6. 下游DDPM负结果同时涉及生成模型拟合和条件读取。**

- AE oracle test delta PCC约0.672。
- 两组DDPM约0.19。
- 两组train也低于均值响应基线。
- validation条件打乱诊断影响较小。
- 增加MC采样不能填补差距。

因此，这次比较不能单独区分“MoA表示不足”与“下游模型没有有效使用MoA”。但它已经足以说明：**当前组合没有达到可用VC效果，也没有可靠MoA增益。**

**7. 目前仍缺少的关键证据**

- 全面、独立well/plate层面的响应和similarity重现性审计。
- 全量pseudocount敏感性，以及低表达/检出率对kernel分子分母贡献的分层分析。
- 在可信标签上，简单query-control模型能提供多少超过candidate prior的信息。
- 下游确定性回归等对照，帮助分离条件可预测性与扩散优化问题。
- 多训练seed的下游验证。
- 全新未触碰cell或跨平台外部验证。

最后，读项目时最容易误解的是 [latest_model.json](F:/CodexWork/moavc/cellmoa_v1/artifacts/latest_model.json)：它仍指向 **R6 `bounded_readout_outer`**，只是历史保留指针；最新embedding实验是R13，最新完成的下游实验是latent DDPM v1。它们均未被推广为已验证可用的模型。
