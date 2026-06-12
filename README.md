# AutoReactGen

Automated generation of LAMMPS `fix bond/react` molecule template files from GROMACS topology (ITP) and coordinate (GRO) files.

## Overview

Setting up a reactive MD simulation in LAMMPS requires pre- and post-reaction molecule templates — files describing the topology and geometry of the reaction site before and after the reaction. Writing these by hand is tedious, error-prone, and scales poorly with molecular complexity.

**AutoReactGen** automates this: given multi-molecule GROMACS topologies, a reaction atom pair, and a bond-distance cutoff, it:

1. Performs BFS over the bond graph to identify all atoms within *N* bonds of the reaction site
2. Extracts the corresponding topology (bonds, angles, proper/improper dihedrals), atom types, and charges from the ITP files
3. Reads coordinates from GRO files
4. Generates a matching pair of LAMMPS-format template files (`.pre` / `.post`) with correct atom renumbering, type strings, and topology filtering

### Example

```bash
python topo.py -mol "BCD,MDI" -atom "146,28" -within 4
```

Reads `BCD.itp`, `MDI.itp`, `BCD-MDI.itp`, and corresponding `.gro` files. Generates `BCD-MDI.pre` and `BCD-MDI.post` — containing the 23 atoms within 4 bonds of reaction atoms 146 (BCD) and 28 (MDI), with all their topology and force field parameters.

### Data files

| File | Content |
|------|---------|
| `BCD.itp` / `MDI.itp` | Individual molecule topologies (147 / 29 atoms) |
| `BCD-MDI.itp` | Crosslinked product topology (176 atoms, BCD+MDI) |
| `BCD.gro` / `MDI.gro` | Individual molecule coordinates (nm) |
| `BCD-MDI.gro` | Product coordinates after crosslinking |
| `box.top` | System-wide atom type definitions |
| `template` | Reference for the LAMMPS molecule template format |

### Key design decisions

- **Min-distance BFS**: in multi-molecule mode, each atom's distance is `min(d_to_A, d_to_B)`. This ensures the template covers the full neighborhood of both reaction sites.
- **Pre/post atom correspondence**: the same atom set appears in both templates, ordered by combined ITP index (molecule 1 atoms first, then molecule 2). This satisfies LAMMPS's requirement that pre/post atoms are 1:1 mapped.
- **Atom types from ITP `[atoms]`**: charges read from column 7 (not `[atomtypes]`, which often contains zeros for GROMACS).
- **Coordinates**: pre-template uses individual GRO files; post-template uses the combined GRO for post-reaction geometry.

### LAMMPS compatibility

The generated templates follow the standard `fix bond/react` molecule file format. For multi-fragment systems, users should define `Fragments` sections and per-fragment RMSD constraints to handle molecules built in separate coordinate frames (see LAMMPS documentation for `fix bond/react` constraint syntax).

---

## Code structure

```
topo.py          — single-file tool (~230 lines)
```

### Parsing layer

| Function | Purpose |
|----------|---------|
| `parse_itp()` | Single-pass ITP parser. Extracts atom count, types, charges, bonds, angles, proper/improper dihedrals from `[atoms]`, `[bonds]`, `[angles]`, `[dihedrals]` sections. Detects improper sections via header comment keyword. |
| `parse_gro()` | GRO coordinate parser (fixed-width format, nm → Å conversion). |

### Graph layer

| Function | Purpose |
|----------|---------|
| `build_adj()` | Builds adjacency list from bond pairs — used by BFS. |
| `_bfs_distmat()` | Multi-source BFS core. Returns distances and per-distance atom groups. Supports optional `max_dist` cutoff. |
| `bfs_within()` | Single-atom wrapper around `_bfs_distmat`; prints distance matrix and returns within-set. |

### Display layer

| Function | Purpose |
|----------|---------|
| `show()` | Prints bonds/angles/dihedrals filtered by atom membership. When `within_set` is provided, keeps only entries where **all** atoms are in the set (induced subgraph view). |

### Template generation

| Function | Purpose |
|----------|---------|
| `_tp_type()` | Builds LAMMPS-style type strings (e.g. `CG-OS-CG`) from atom IDs and the type lookup table. |
| `write_template()` | Writes a `.pre` or `.post` file in the standard LAMMPS molecule template format: header counts, `Coords`, `Types`, `Charges`, `Molecules`, `Bonds`, `Angles`, `Dihedrals`, `Impropers`. |
| `gen_templates()` | Orchestrates the full pipeline: gathers pre/post data from individual and combined ITP/GRO files, calls `write_template()` for both templates. |

### Flow

```
main() → parse all ITPs → individual BFS (per molecule) → multi-source BFS (combined)
     → union atom set → gen_templates() → write_template() × 2 (.pre + .post)
     → show() filtered topology (optional)
```
