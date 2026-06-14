# AutoReactGen 方案与进度

## 目标

将当前 `topo.py` 从“反应邻域模板生成脚本”逐步升级为由 `box.top`/`box.gro` 驱动的静态生成工具：

```bash
python topo.py -mol "BCD,MDI" -atom "146,28" -box box
```

最终输出：

- `BCD-MDI.pre`
- `BCD-MDI.post`
- `BCD-MDI.map`
- `box.data`
- `box.in`

## 核心路线

1. **增强 parser**：保留现有 `parse_itp()` 行为，同时新增 rich parser，完整保留 atomtypes、atoms、bonds、angles、proper dihedrals、improper dihedrals 的参数。✅ 已完成第一版
2. **引入 TypeRegistry**：所有 atom/bond/angle/dihedral/improper 类型统一分配数字 ID。✅ 已完成第一版
3. **模板改用数字 ID**：`.pre/.post` 不再写 `OS-CG` 这类字符串 label，避免 LAMMPS `labelmap` 问题。⬜ 待做
4. **实现 box.data**：按 `box.top [molecules]` 展开 BCD×20 + MDI×100，使用 `box.gro` 坐标。⬜ 待做
5. **实现 box.in**：生成保守的 LAMMPS 输入文件，先覆盖 force-field style、coeff、read_data、molecule、bond/react 骨架。⬜ 待做
6. **无 LAMMPS 静态验证**：用 Python 检查计数、引用、类型 ID、坐标、单位转换和 proper/improper 分离。🔶 部分完成：parser/TypeRegistry debug 检查已加入

## 已完成：parser 增强第一版

已在 `topo.py` 中加入 rich parser，同时保持旧接口兼容。

### 新增数据结构

- `AtomType`：保留 atom type 的原子序数、质量、电荷、ptype、sigma、epsilon、注释。
- `Atom`：保留 atom id、type、residue/resname、atom name、cgnr、charge、mass、注释。
- `Bond`：保留 ai/aj、funct、参数列、注释。
- `Angle`：保留 ai/aj/ak、funct、参数列、注释。
- `Dihedral`：保留 ai/aj/ak/al、funct、参数列、proper/improper 标记、注释。
- `ItpTopology`：聚合单个 ITP 的 atomtypes、atoms、bonds、angles、propers、impropers、moleculetype、nrexcl。
- `TopologyFile`：聚合 `.top` 的 defaults、atomtypes、include、system、molecules。

### 新增 parser

- `parse_itp_rich(filepath)`：解析 ITP 并保留完整参数。
- `parse_top(filepath)`：解析 `box.top` 的系统层信息。
- `parse_itp(filepath)`：保留原来的 7 元组返回值，继续服务现有 BFS/template 逻辑。

### Parser 验证结果

rich parser 统计：

```text
BCD rich 147 6 154 287 546 0
MDI rich 29 7 30 46 66 12
BCD-MDI rich 176 21 185 337 622 14
```

legacy parser 兼容统计：

```text
BCD legacy 147 154 287 462 0
MDI legacy 29 30 46 66 12
BCD-MDI legacy 176 185 337 535 14
```

`rich` 与 `legacy` 的 proper 数量不同是预期行为：`rich` 保留重复/多项 dihedral，`legacy` 继续用 set 去重，避免破坏现有模板生成逻辑。

`box.top` 解析结果：

```text
includes ['BCD.itp', 'MDI.itp']
molecules [('BCD', 20), ('MDI', 100)]
atomtypes 13
defaults {'nbfunc': '1', 'comb-rule': '2', 'gen-pairs': 'yes', 'fudgeLJ': '0.5', 'fudgeQQ': '0.8333'}
system BCD+MDI
```

通过检查：

- `python -m py_compile topo.py`
- `python topo.py -mol "BCD,MDI" -atom "146,28" -within 4`

## 已完成：TypeRegistry 第一版

已在 `topo.py` 中加入 `TypeRegistry`，作为后续 `box.data`、`box.in`、`.pre/.post` 数字类型输出的共享类型源。

### 当前能力

- 独立 namespace：`atom`、`bond`、`angle`、`dihedral`、`improper`。
- 稳定递增 ID：同一 signature 重复注册返回同一 ID。
- atom type 按 **大小写敏感 name + 参数** 注册；`OS` 和 `os` 分开处理，不能因参数接近合并。
- bond type 端点对称：`A-B` 与 `B-A` 复用同一类型。
- angle type 保留中心原子，外侧两端对称。
- proper 与 improper 使用独立 namespace，互不合并。
- key 中使用 LAMMPS real 单位换算后的参数，方便后续生成 coeff section。
- 新增 `--debug-types` 调试入口。

### 当前单位换算

- atom σ：nm × 10 → Å
- atom ε：kJ/mol ÷ 4.184 → kcal/mol
- bond r0：nm × 10 → Å
- bond K：÷ 836.8
- angle K：÷ 8.368
- proper dihedral K：÷ 4.184
- improper K：÷ 8.368

### TypeRegistry 验证结果

命令：

```bash
python topo.py -mol "BCD,MDI" -atom "146,28" --debug-types
```

当前输出：

```text
TypeRegistry counts:
  atom: 21
  bond: 19
  angle: 30
  dihedral: 71
  improper: 6
```

单独 topology 注册统计：

```text
BCD.itp     {'atom': 6,  'bond': 6,  'angle': 12, 'dihedral': 40, 'improper': 0}
MDI.itp     {'atom': 7,  'bond': 7,  'angle': 9,  'dihedral': 13, 'improper': 3}
BCD-MDI.itp {'atom': 21, 'bond': 19, 'angle': 30, 'dihedral': 71, 'improper': 6}
```

行为检查：

- `OS` 与 `os` 分配不同 atom type ID。
- 同名同参数重复注册复用 ID。
- 同名但参数不同分配新 ID。
- `Bond(i,j)` 与 `Bond(j,i)` 复用 ID。
- `Angle(i,j,k)` 与 `Angle(k,j,i)` 复用 ID。

## 已知限制

- `parse_top()` 目前只记录 include，不递归展开 include；后续 `build_system()` 阶段处理。
- bonded 参数当前保留为 float tuple；如果未来遇到宏或非数字参数，应补充原始 token 保存。
- `[ pairs ]` 暂不解析；初版使用自动产生/`[defaults]`/`special_bonds` 路径，后续如需精确 1-4 pair 再补。
- include 路径当前按文件中写法保存；若支持子目录输入，需要 resolve 到 `.top` 所在目录。
- TypeRegistry 目前仅用于注册和 debug，不会默认改变 `.pre/.post` 的字符串类型输出。

## 下一步

推荐继续做 **模板数字 ID 输出** 或 **box.data 骨架**。更稳的顺序是：

1. 让模板生成路径可选使用 TypeRegistry 输出数字 ID；默认仍保持旧字符串输出。
2. 将模板生成从 legacy `parse_itp()` 逐步切到 `parse_itp_rich()`，这样 bonded type ID 能包含 functype 和参数。
3. 实现 `box.data` writer：使用同一 TypeRegistry 写 Masses、Pair Coeffs、Bond/Angle/Dihedral/Improper Coeffs 与拓扑 sections。
4. 实现 `box.in` writer：从 TypeRegistry 输出 styles、coeff/read_data/molecule/bond-react 骨架。
