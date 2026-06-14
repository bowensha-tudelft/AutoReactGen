# AutoReactGen 问题与方案

## 背景

AutoReactGen 从 GROMACS ITP/GRO 文件生成 LAMMPS `fix bond/react` 的 `.pre/.post/.map` 模板。当前仅实现了 `topo.py`（~230行）完成模板生成，但完整的 LAMMPS 模拟还需要 **data 文件**（体系拓扑+坐标+参数）和配套的 **in 文件**。

**现有 data 文件**（`test/primary_core.data`，1.5MB）由 InterMol 转换 → 外部脚本清洗得来。实际运行 `test/in.primary_core` 时报错：
```
ERROR: Invalid bond type OS-CG in Bonds section of molecule file
```
根因：模板用字符串类型标签，但 LAMMPS 输入缺少 `labelmap`。

**体系**：BCD（147原子，6种类型）×20 + MDI（29原子，7种类型）×100 → 5840原子，6080键，10340角，18720二面角，**0 improper**（已知缺陷）

**两种 ITP 格式并存**：

| 来源 | 键格式 | 角格式 | 二面角格式 |
|------|--------|--------|-----------|
| BCD (Q4MC-CD) | `ai aj funct c0(r0) c1(K)` | 同左格式 | functype=1: `phase K n` |
| MDI / BCD-MDI (Sobtop) | `ai aj functype r0 K` | `a0 K`（显式列名） | functype=9 propers, functype=4 impropers |

---

## 核心决策：三种方案对比

### 方案 A：修补 InterMol + 保留外部脚本

修改 InterMol 的 `lammps_parser.py` 3 处（改 2 行 class 检查 + 取消 dead branch 注释），继续用清洗脚本组装多分子体系。

- **优点**：改动量极小（~5 行）；ParmEd 自动单位换算；二面角转换已有成熟逻辑
- **缺点**：依赖外部项目（InterMol v0.1.0.dev0，10 年未出正式版）；碎片化流程（InterMol→清洗→topo）；类型编号需手动对齐

### 方案 B：完全自研，零外部依赖

从头实现 ITP 解析、参数提取、单位换算、data 写入、多分子组装。

- **优点**：完全掌控；类型编号天然统一；无版本兼容风险
- **缺点**：约 300-400 行新代码；需自行实现所有转换公式；二面角转换可能踩坑

### 方案 C：混合策略 — 借 InterMol 基础设施，自研输出层 ★推荐

```
┌─────────────────────────────────────────────────┐
│ InterMol (借用)           自研 (新增)            │
│ ┌──────────────┐    ┌──────────────────────┐    │
│ │ ITP 解析      │───▶│ TypeRegistry         │    │
│ │ (unit 系统)   │    │ 类型去重 & 编号映射    │    │
│ │              │    └──────┬───────────────┘    │
│ │ 参数提取      │         │                     │
│ │ (force kwds) │         ▼                     │
│ │              │    ┌──────────────────────┐    │
│ │ 单位换算      │    │ write_data.py        │    │
│ │ (ParmEd)     │───▶│ Data 文件各 section   │    │
│ │              │    │ (含 Impropers ✅)     │    │
│ └──────────────┘    └──────┬───────────────┘    │
│                            │                     │
│ topo.py (已有)             ▼                     │
│ ┌──────────────┐    ┌──────────────────────┐    │
│ │ parse_itp    │    │ build_system.py       │    │
│ │ parse_gro    │───▶│ 多分子坐标 & 拓扑组装  │    │
│ │ BFS 拓扑查询  │    └──────┬───────────────┘    │
│ │ gen_templates│          │                     │
│ │ write_map    │          ▼                     │
│ └──────────────┘    ┌──────────────────────┐    │
│                      │ write_in.py           │    │
│                      │ LAMMPS 输入文件生成    │    │
│                      └──────────────────────┘    │
└─────────────────────────────────────────────────┘
```

**具体分工**：

| 组件 | 来源 | 用途 |
|------|------|------|
| ITP 解析 | InterMol `gromacs_parser` | 解析 `[atoms]`, `[bonds]`, `[angles]`, `[dihedrals]` |
| 参数 + unit | InterMol `forcedata` | 参数名→列位置映射、ParmEd 单位对象 |
| 单位换算 | ParmEd `units` | `value_in_unit(units.kilocalorie_per_mole)` 等 |
| 二面角转换 | InterMol `convert_dihedrals` | RB→Trig 展开（如需 multi/harmonic 输出） |
| TypeRegistry | **自研** | 去重、分配数字类型 ID |
| Data 写入 | **自研** `write_data.py` | 正确区分 propers/impropers；按 LAMMPS 格式输出 |
| 多分子组装 | **自研** `build_system.py` | atom/bond/angle/dihedral index 偏移；GRO 坐标平移 |
| 模板生成 | `topo.py`（已有，微调） | 改用 TypeRegistry 数字类型 ID 输出 |
| in 文件 | **自研** `write_in.py` | pair_coeff, bond/angle/dihedral style, molecule 定义 |

**方案 C 优势**：
- 单位换算、二面角数学不重造轮子——用 InterMol + ParmEd 经过上千个体系验证的路径
- Improper bug 自然规避——我们只借 InterMol 的**解析和单位**，不借它的**写入**逻辑
- 类型编号全程掌控——TypeRegistry 在 data 和模板间共享
- 如果未来 InterMol 修了 bug 出新版，我们也可以切换回去

**需要给 InterMol 打的补丁**（可先 monkey-patch，不必改源码）：
- `gromacs_parser.py:181` 把 `'4': ProperPeriodicDihedral` 改为 `'4': ImproperHarmonicDihedral`——或者我们自己写 wrapper，在解析后把 functype=4 的条目从 dihedrals 分离出来

### 二面角简化策略

InterMol 将 GROMACS proper periodic 转换为 `multi/harmonic`（A1-A5 系数）或 `charmm`（K n d weight），取决于 phase 是否为 0/180。`test/primary_core.data` 中 30 种二面角类型中有 28 种是 multi/harmonic，仅 2 种是 charmm。

自研时可**统一使用 `dihedral_style fourier`**（GAFF 原生风格）：
- 每条 GROMACS proper periodic 项 → 一个 fourier 二面角类型
- `fourier` 风格是 LAMMPS 中对 GAFF 力场的标准选择，支持任意 phase 值
- 格式：`ID fourier K n d`（K: kcal/mol, n: integer, d: degrees）
- 同一四元组上多条不同 n 的项 → 多条 fourier 条目叠加（LAMMPS 自动求和）
- 转换：`K_lmp = K_gro / 4.184`，n 和 d 不变
- 避免了 RB→multi/harmonic 的 Chebyshev 多项式展开

### InterMol 源码级分析

InterMol 的转换管道是：GROMACS ITP → `ProperPeriodicDihedral` → `TrigDihedral`（8系数规范形式）→ `charmm`/`multi/harmonic`。

**Improper bug 精确追踪**：

| 位置 | 代码 | 问题 |
|------|------|------|
| `gromacs_parser.py:997` | `improper = numeric_dihedraltype in ['2','4']` | functype=4 正确检测为 improper ✓ |
| `gromacs_parser.py:181` | `'4': ProperPeriodicDihedral` | 映射为 ProperPeriodicDihedral（非 ImproperHarmonicDihedral） |
| `gromacs_parser.py:1007-1010` | `['1','3','4','5','9']` → `TrigDihedral` | functype=4 成为 TrigDihedral，improper 标志保留在 `.improper` 属性 |
| `lammps_parser.py:192` | `if False: #dihedral.improper` | **dead branch** — 预留的修复路径从未启用 |
| `lammps_parser.py:952-953` | `force.__class__ != ImproperHarmonicDihedral` → Dihedrals | TrigDihedral ≠ ImproperHarmonicDihedral，==> improper 错入 Dihedrals |
| `lammps_parser.py:962-963` | `force.__class__ == ImproperHarmonicDihedral` → Impropers | TrigDihedral ≠ ImproperHarmonicDihedral，improper 被排除 |

**修复 InterMol 只需改 2 行**（`lammps_parser.py:952-953,962-963`）：
```python
# write_dihedrals: 从 class 判断改为 improper 属性判断
dihedral_forces = {f for f in mol_type.dihedral_forces if not getattr(f, 'improper', False)}
# write_impropers: 同上
improper_forces = {f for f in mol_type.dihedral_forces if getattr(f, 'improper', False)}
```
同时取消 `lammps_parser.py:192` 的 dead branch 注释，使 functype=4 improper 走正确的写入路径。

**1/2 缩放因子**（`lammps_parser.py:58-59`）：
- `SCALE_INTO = 2.0`，`SCALE_FROM = 0.5`
- 仅用于 harmonic bond/angle/improper 的 `k` 参数（`canonical_bond:92-93`, `canonical_angle:126-127`, `canonical_dihedral:210`）
- 二面角 K / LJ ε 不经过此缩放

**GROMACS ↔ LAMMPS 符号约定差异**（`gromacs_parser.py:236-238`）：
- GROMACS RB 二面角转换时对 C1/C3/C5 取反（psi → phi 符号约定），LAMMPS 侧不做此处理
- 这仅影响 RB（functype=3）二面角，GAFF 使用 functype=1/9/4 不受影响

---

## 问题 1：Labelmap 缺失 → 模板类型标签无法解析 [P0]

**现象**：`BCD-MDI.pre` 用字符串类型标签（`OS-CG` 等），`test/in.primary_core` 无 `labelmap`，LAMMPS 报错。

**根因**：LAMMPS molecule 模板的 Types/Bonds/Angles 等 section 可使用字符串标签或数字 ID。字符串标签必须通过 `labelmap`（在 `molecule` 命令或 data 文件的 Type Labels section 中定义）映射到数字 ID。

**短期修复**：在 `molecule` 命令添加 labelmap：
```lammps
molecule pre_bcd_mdi ../BCD-MDI.pre
  labelmap bond OS-CG 1 CG-OS 2 CG-CG 3 ...
  labelmap angle OS-CG-CG 1 ...
```
**长期方案**：自研生成器中模板直接输出数字类型 ID，彻底消除 labelmap 需求。通过 `TypeRegistry` 统一 data 和模板的类型编号。

---

## 问题 2：Improper 二面角缺失 [P1]

**现象**：`test/primary_core.data` 显示 `0 impropers`，但 MDI.itp 的 `[ dihedrals ] ; impropers` 包含 12 条 functype=4 improper。

**根因**：
- InterMol bug：functype=4 被标记 `improper=True` 但映射为 `ProperPeriodicDihedral → TrigDihedral` 类，写入 LAMMPS 时只有 `ImproperHarmonicDihedral` 类才进入 Impropers section
- 上游清洗脚本额外硬编码 `n_impropers = 0`

**自研方案修复**：解析时检测 `'improper' in section_header.lower()` → 直接归类为 improper；生成 data 时独立写入 Impropers 和 Improper Coeffs section。转换：`K_lammps = K_gro / 8.368`（LAMMPS harmonic improper 内嵌 1/2 因子，与 angle 同理）。

---

## 问题 3：缺少完整的 data 文件生成器 [P0]

**现状**：`topo.py` 只生成 pre/post/map 模板。data 文件由外部流水线生成，不在仓库内。

**目标**：新增 `write_data()` 模块，输入 box.top + 各 ITP + GRO → 输出完整 data 文件。需实现：

1. **原子类型枚举**：从 box.top 的 `[atomtypes]` 合并去重各 ITP 的 atom types → 分配 1..N 编号
2. **参数提取 & 单位转换**：从 ITP 各 section 提取 bond/angle/dihedral/improper 参数，应用转换公式
3. **拓扑类型去重**：每种独特的 (原子类型组合, 参数) → 一个 bond/angle/dihedral/improper 类型 ID
4. **多分子体系组装**：按 box.top 中 `[molecules]` 的计数，偏移 atom index、bond/angle/dihedral index
5. **坐标生成**：从 GRO 读取单分子坐标，为每个分子实例平移（或使用预建的大 GRO）
6. **格式化输出**：按 LAMMPS data 格式写入各 section

同时生成配套的 in 文件（或 in 文件模板），包含 pair_coeff、bond/angle/dihedral style 定义。

---

## 问题 4：模板与 data 类型编号一致性 [P0]

**原子类型映射**（需在 data 和模板间共享）：

| ID | Name | Mass | σ (Å) | ε (kcal/mol) | 来源 |
|----|------|------|-------|-------------|------|
| 1 | CG | 12.010 | 3.3996695 | 0.1094000 | BCD |
| 2 | H2 | 1.008 | 2.2931733 | 0.0157000 | BCD |
| 3 | OS | 16.000 | 3.0000123 | 0.1700000 | BCD |
| 4 | H1 | 1.008 | 2.4713530 | 0.0157000 | BCD |
| 5 | OH | 16.000 | 3.0664734 | 0.2104000 | BCD |
| 6 | HO | 1.008 | 0.0000000 | 0.0000000 | BCD |
| 7 | c3 | 12.011 | 3.3996700 | 0.1094000 | MDI |
| 8 | hc | 1.008 | 2.6495330 | 0.0157000 | MDI |
| 9 | ca | 12.011 | 3.3996700 | 0.0860000 | MDI |
| 10 | ha | 1.008 | 2.5996420 | 0.0150000 | MDI |
| 11 | ne | 14.007 | 3.2499990 | 0.1700000 | MDI |
| 12 | c1 | 12.011 | 3.3996700 | 0.2100000 | MDI |
| 13 | o | 15.999 | 2.9599220 | 0.2100000 | MDI |
| 14 | hn | 1.008 | 1.0690780 | 0.0157000 | post |
| 15 | n | 14.007 | 3.2499990 | 0.1700000 | post |
| 16 | c | 12.011 | 3.3996700 | 0.0860000 | post |

（与 `test/in.primary_core` 中 `pair_coeff` 顺序一致。Post 类型 hn/n/c 仅在交联产物 BCD-MDI.itp 中出现。）

**注意大小写去重**：BCD-MDI.itp 中同时有 `CG`（来自 BCD）和 `c3`（来自 MDI）、`OS` 和 `os`、`H1` 和 `h1` 等近似类型。去重策略：**atom type 名称大小写敏感**，因此 `OS` 与 `os`、`H1` 与 `h1` 必须分开处理；仅当 **名称完全相同且参数在容差内一致** 时才复用类型 ID。

**当前进展**：`TypeRegistry` 第一版已在 `topo.py` 中实现，用于为 atom/bond/angle/dihedral/improper 分配稳定数字 ID。当前仅作为注册与 debug 层使用，尚未接入 `.pre/.post` 数字输出或 `box.data` writer。调试命令：
```bash
python topo.py -mol "BCD,MDI" -atom "146,28" --debug-types
```
当前统计：`atom=21, bond=19, angle=30, dihedral=71, improper=6`。

**解法**：后续 data 生成和模板生成共享同一 `TypeRegistry` 实例，保证 atom/bond/angle/dihedral/improper 类型 ID 完全一致。

---

## 附录：单位换算速查表（GROMACS → LAMMPS `real`）

| 参数 | GROMACS 单位 | LAMMPS 单位 | 换算 |
|------|-------------|------------|------|
| 坐标 / 键长 r₀ / LJ σ | nm | Å | **×10** |
| 键力常数 K | kJ/(mol·nm²) | kcal/(mol·Å²) | **÷ 836.8** ① |
| 角力常数 K | kJ/(mol·rad²) | kcal/(mol·rad²) | **÷ 8.368** ① |
| 角 θ₀ / χ₀ / 二面角相位 d | deg | deg | 不变 |
| 二面角力常数 K | kJ/mol | kcal/mol | **÷ 4.184** |
| 二面角多重度 n | — | — | 不变 |
| LJ ε | kJ/mol | kcal/mol | **÷ 4.184** |
| 电荷 / 质量 | e / amu | e / g/mol | 不变 |

> ① **关键**：LAMMPS harmonic 风格能量公式为 `E = K(x-x₀)²`，已内嵌 1/2 因子；GROMACS 公式为 `V = ½K(r-r₀)²`。因此 LAMMPS K 需额外 ÷2。

**转换推导示例**（键）：
K_lmp = K_gro × (1/2) × (1/4.184) × (1/100) = K_gro / 836.8

**验证**（来自 `test/primary_core.data`）：BCD bond `17-18`（CG-H1）：r₀=0.1092nm, K_gro=343088 → r₀_lmp=1.092Å ✓, K_lmp=343088/836.8=410.0 ✓

### 二面角风格选择

自研生成器建议统一使用 `dihedral_style charmm`：
- data 文件 Dihedral Coeffs 格式：`ID charmm K n d weight`
- 每条 GROMACS proper periodic（functype=1 或 9）→ 一条 charmm 项
- 同一四元组上多条（不同 n）→ 多条 charmm 项，LAMMPS 自动叠加

### GROMACS functype 对照

| GROMACS | 含义 | LAMMPS style | 写入 section |
|---------|------|-------------|-------------|
| bond functype 1 | harmonic bond | `harmonic` | Bond Coeffs |
| angle functype 1 | harmonic angle | `harmonic` | Angle Coeffs |
| dihedral functype 1 / 9 | proper periodic | `charmm` | Dihedral Coeffs |
| dihedral functype 4 | improper harmonic | `harmonic` | **Improper Coeffs** |

---

## 下一步行动计划

| 优先级 | 任务 | 预计代码量 | 依赖 |
|--------|------|----------|------|
| **P0** | 实现 `write_data.py`：ITP参数提取 + 单位转换 + data 文件输出 | ~150 行 | TypeRegistry |
| **P0** | 实现 `TypeRegistry`：原子/键/角/二面角类型全局注册与去重 | 已完成第一版 | 无 |
| **P1** | 扩展 `topo.py`：模板改用 TypeRegistry 的数字类型 ID 输出 | ~30 行修改 | TypeRegistry |
| **P1** | 多分子体系组装 + GRO 坐标读取与平移 | ~80 行 | data 生成器 |
| **P1** | in 文件生成（pair_coeff, bond/angle/dihedral style, molecule 定义等） | ~60 行 | TypeRegistry |
| **P2** | 在有 LAMMPS 的环境中验证：运行 `test/in.primary_core` | — | 以上全部 |
| **P2** | 短期验证用：手动在 `test/in.primary_core` 加 labelmap，确认模板格式正确 | 手动 | 无 |

### 当前决策

- **采用方案 C（混合策略）**：借 InterMol 解析 ITP + 单位换算，自研 data 写入和体系组装
- `parse_itp_rich()` / `parse_top()` 已完成第一版；`TypeRegistry` 已完成第一版，当前用于稳定类型注册和 debug 统计
- Atom type 名称大小写敏感：`OS`/`os`、`H1`/`h1` 等分开处理，不跨大小写合并
- `[ pairs ]` 暂不解析，先使用自动产生/`special_bonds` 路径
- InterMol 需小修补：将 functype=4 的 `ProperPeriodicDihedral` 映射改为正确路径（或事后从 dihedral_forces 中按 `.improper` 属性筛出）
- 二面角统一用 `fourier` 风格（GAFF 原生 proper periodic 形式），避免 `multi/harmonic` 的 5 项级数展开
- 键/角/improper 力常数注意 1/2 因子（÷836.8 / ÷8.368）
- 下一步优先：模板可选数字 ID 输出，随后实现 `box.data` writer
