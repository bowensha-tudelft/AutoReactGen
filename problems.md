# AutoReactGen 当前问题记录

## 背景

AutoReactGen 是一个自动化生成 LAMMPS `fix bond/react` 分子模板文件的工具，从 GROMACS 拓扑 (ITP) 和坐标 (GRO) 文件生成 `.pre`、`.post`、`.map` 模板文件。

当前生成的 `BCD-MDI.pre/post/map` 已接入 `test/in.primary_core`，但 LAMMPS 读取 `BCD-MDI.pre` 时失败。经过排查，存在以下三个核心问题。

---

## 问题 1：Labelmap 缺失导致 LAMMPS 无法解析模板

### 现象描述

`BCD-MDI.pre` / `.post` 模板文件中使用了字符串类型标签（如 `OS-CG`、`CG-OS-CG`），但 `test/in.primary_core` 中没有对应的 `labelmap` 指令。运行后 LAMMPS 报错：

```
ERROR: Invalid bond type OS-CG in Bonds section of molecule file
```

### 根因分析

LAMMPS `molecule` 文件支持两种类型标识方式：
1. **数字类型 ID**：直接引用 `data` 文件中的数值类型
2. **字符串标签**：通过 `labelmap` 指令将字符串映射为数字类型 ID

当前 `test/in.primary_core` 的 `molecule` 定义中没有 `labelmap` 指令，因此 LAMMPS 将模板文件中的字符串标签直接解析为数字，导致类型匹配失败。

### 影响范围

- **阻塞 smoke test**：任何包含字符串标签的模板文件均无法被 LAMMPS 读取
- **影响文件**：`test/in.primary_core` 中所有 `molecule` 定义（`BCD_Molecule`、`MDI_Molecule`、`BCD-MDI_Molecule`）

### 修复建议

在 `test/in.primary_core` 中每个 `molecule` 定义后添加 `labelmap` 指令，映射模板中的字符串标签到 `primary_core.data` 中的数字 type ID：

```lammps
labelmap bond OS-CG 7      # 对应 primary_core.data 中的 Bond Coeffs 索引
labelmap angle CG-OS-CG 11 # 对应 primary_core.data 中的 Angle Coeffs 索引
# ... 补充所有在模板中使用的 bond/angle/dihedral/improper 类型标签
```

---

## 问题 2：InterMol 将 Impropers 错误归类到 Dihedrals 中

### 现象描述

`primary_core.data` 中显示 `0 impropers / 18720 dihedrals`，但原始 GROMACS 拓扑 `MDI.itp` 中 `[ dihedrals ] ; impropers` 实际包含 12 条 functype=4 的 improper。这些 improper 被错误地写入了 Dihedrals section。

### 根因分析

**InterMol GROMACS Parser** (`gromacs_parser.py:997-1018`)：
- functype=4 的 dihedrals 被标记了 `improper=True`
- 但在 `['1','3','4','5','9']` 分支中被映射为 `TrigDihedral` 类，而非 `ImproperHarmonicDihedral` 类

**LAMMPS Writer** (`lammps_parser.py:950-965`)：
- `write_dihedrals()` 排除 `ImproperHarmonicDihedral` 类 → `TrigDihedral` 类的 impropers 被包含进 Dihedrals
- `write_impropers()` 只选择 `ImproperHarmonicDihedral` 类 → `TrigDihedral` 类的 impropers 被排除

**结果**：functype=4 的 impropers 被写入了 Dihedrals section，导致 `Impropers` 和 `Improper Coeffs` 均为 0。

### 影响范围

- **分子拓扑不完整**：MDI 残基中的 improper 约束丢失
- **下游脚本链断裂**：即使手动修正 data 文件，上游脚本（见问题 3）也会将其丢弃
- **模板生成**：当前 `topo.py` 已能生成含 `Impropers` section 的模板，但无法与 `0 impropers` 的 data 文件匹配

### 修复建议

**方案 A（推荐）：修改 InterMol 源码**

在 LAMMPS output 时，根据 `force.improper` 标志位（而非 `__class__`）来分类 dihedral/improper：

```python
# lammps_parser.py:write_dihedrals()
# 原逻辑：
# if not isinstance(force, ImproperHarmonicDihedral):

# 修复后：
# if not force.improper:
#     # 写入 Dihedrals

# lammps_parser.py:write_impropers()
# 原逻辑：
# if isinstance(force, ImproperHarmonicDihedral):

# 修复后：
# if force.improper:
#     # 写入 Impropers
```

**方案 B（绕过 InterMol）**

在生成 LAMMPS data 文件后，编写后处理脚本：
1. 扫描 Dihedrals section，识别 functype=4 的条目（impropers）
2. 将其移动到 Impropers section
3. 同步移动对应的 Improper Coeffs

**方案 C（临时规避）**

在 InterMol 修复前，让 `topo.py` 暂不输出 impropers，与当前 InterMol 产物保持一致。但这只是权宜之计， improper 的物理约束仍缺失。

---

## 问题 3：上游脚本硬编码 `0 impropers`

### 现象描述

在 `/home/bsha/md/crosslink/BCD-MDI/` 目录下，`clean_lammps_data.py` 和 `make_reaction_assets.py` 硬编码了 `0 impropers`。即使 InterMol 正确输出了 impropers，这些脚本也会将其丢弃。

### 根因分析

这些脚本在生成/清理 LAMMPS data 文件时，直接写入：

```python
# clean_lammps_data.py / make_reaction_assets.py（简化示意）
f.write(f"{n_impropers} impropers\n")  # n_impropers 硬编码为 0
```

没有读取或传递 InterMol 输出的 impropers 信息，也没有对应的 `Improper Coeffs` 和 `Impropers` section 写入逻辑。

### 影响范围

- **数据文件不完整**：即使 InterMol 修正了 improper 分类，这些脚本仍会将其过滤掉
- **模板与数据不匹配**：`BCD-MDI.pre/post` 中若包含 `Impropers`，则无法与 `0 impropers` 的 data 文件配合使用
- **需要重新生成数据**：修复脚本后，必须重新运行整个数据处理流水线

### 修复建议

1. **修改 `clean_lammps_data.py`**：从 InterMol 输出中读取 impropers 数量，保留 `Improper Coeffs` 和 `Impropers` section
2. **修改 `make_reaction_assets.py`**：传递 `n_impropers` 参数，不再硬编码为 0
3. **重新生成 `primary_core.data`**：修复后重新运行流水线，确保包含正确的 impropers

---

## 下一步行动计划

| 优先级 | 任务 | 目标文件/位置 | 预期结果 |
|--------|------|-------------|----------|
| P0 | 在 `in.primary_core` 中添加 labelmap 指令 | `test/in.primary_core` | LAMMPS 能解析模板中的字符串标签 |
| P1 | 修复 InterMol improper 分类逻辑 | `lammps_parser.py:950-965` | functype=4 条目正确归入 Impropers |
| P1 | 修复上游脚本硬编码问题 | `clean_lammps_data.py`, `make_reaction_assets.py` | data 文件正确保留 impropers |
| P2 | 重新生成 `primary_core.data` | `/home/bsha/md/crosslink/BCD-MDI/` | 包含 12 个 impropers 的完整数据文件 |
| P2 | 验证 smoke test | `test/in.primary_core` | LAMMPS 成功读取模板并完成反应模拟 |

### 当前决策

- **立即执行 P0**：添加 labelmap 指令，验证非 improper 类型（bond/angle/dihedral）的 labelmap 是否正确
- **P1 方案选择**：优先尝试方案 A（修改 InterMol 源码），若改动较大则采用方案 B（后处理脚本）
- **P1 同步修复**：上游脚本与 InterMol 修复同步进行，避免重复工作
