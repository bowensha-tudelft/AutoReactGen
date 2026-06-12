# 当前问题调查记录

当前生成的 `BCD-MDI.pre/post/map` 已经接入 `test/in.primary_core`，但 LAMMPS 读取 `BCD-MDI.pre` 时失败。初步判断是：分子模板里的类型标签（例如 `OS-CG`）需要与 `primary_core.data` 中的数值 bond/angle/dihedral/improper type 对应的 labelmap 定义一致，否则模板无法正确解析。

另外，`primary_core.data` 目前是 `0 impropers / 0 improper types`，而新版本 `topo.py` 生成的 `BCD-MDI.pre/post` 里包含 `Impropers` 段，这很可能与现有数据文件结构不匹配，直接导致 smoke test 失败。

初步排查还发现，`/home/bsha/md/crosslink/BCD-MDI` 下的 `clean_lammps_data.py` 和 `make_reaction_assets.py` 要么硬编码、要么省略 impropers，因此这些脚本链路本身就可能把 impropers 丢掉。

InterMol 本地测试也显示了类似问题：使用 conda 环境 `apr`，把 `/home/bsha/AutoReactGen/MDI.gro + MDI.itp` 转成 `/tmp/mdi_intermol_test/mdi_test.lmp` 后，源 `MDI.itp` 里 `[ dihedrals ] ; impropers` 有 12 条，但输出的 LAMMPS data 里 `Improper Coeffs/Impropers` 都是 0。说明在这条本地转换路径上，InterMol 也没有保留这些 impropers。

结论上，要尽快跑通 smoke test，只有两条路：
1. 让 `topo.py` 暂时不输出 impropers，与当前 InterMol 产物保持一致；
2. 修复转换/清理流水线，保留 LAMMPS data 里的 improper coeffs 和 topology，然后继续在生成模板中保留 impropers。

## 下一步
- 先验证非 improper 类型的 labelmap 是否正确。
- 决定是否临时禁用模板中的 impropers。
- 继续排查/修补 InterMol 的 improper 分类，或者改用别的转换器、`parmed`、自定义 parser。
- 如果最终保留 impropers，则需要重新生成 `primary_core.data`，确保包含 `Improper Coeffs` 和 `Impropers`。
