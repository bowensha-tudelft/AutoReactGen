#!/usr/bin/env python3
"""Query topology (bonds/angles/dihedrals) involving a given atom from a GROMACS ITP file.
Also generates LAMMPS reaction molecule templates (pre/post)."""

import argparse
from collections import deque
from pathlib import Path


def parse_itp(filepath):
    """Return (natoms, bonds, angles, propers, impropers, types, charges)."""
    natoms = 0
    bonds, angles = set(), set()
    propers, impropers = set(), set()
    types, charges = {}, {}
    section = None
    improper = False

    with open(filepath) as f:
        for raw in f:
            line = raw.split(';')[0].strip()
            if not line:
                continue
            if line.startswith('['):
                section = line.strip('[]').strip().split()[0]
                improper = 'improper' in raw.lower()
                continue

            tok = line.split()
            if section == 'atoms':
                aid = int(tok[0])
                types[aid] = tok[1]
                charges[aid] = float(tok[6])
                natoms += 1
            elif section == 'bonds':
                bonds.add((int(tok[0]), int(tok[1])))
            elif section == 'angles':
                angles.add((int(tok[0]), int(tok[1]), int(tok[2])))
            elif section == 'dihedrals':
                d = (int(tok[0]), int(tok[1]), int(tok[2]), int(tok[3]))
                (impropers if improper else propers).add(d)

    return natoms, bonds, angles, propers, impropers, types, charges


def parse_gro(filepath):
    """Return list of (x, y, z) tuples, 1-indexed."""
    with open(filepath) as f:
        lines = f.readlines()
    natoms = int(lines[1].strip())
    coords = [None] * (natoms + 1)
    for i in range(natoms):
        line = lines[2 + i]
        # GRO coords in nm → convert to Angstrom
        coords[i + 1] = (float(line[20:28]) * 10, float(line[28:36]) * 10, float(line[36:44]) * 10)
    return coords


def build_adj(bonds):
    adj = {}
    for a, b in bonds:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    return adj


def _bfs_distmat(starts, bonds, max_dist):
    """Multi-source BFS. Returns (dist dict, by_dist dict)."""
    adj = build_adj(bonds)
    dist = {}
    by_dist = {}
    q = deque()
    for s in starts:
        if s not in adj:
            continue
        dist[s] = 0
        by_dist.setdefault(0, []).append(s)
        q.append(s)

    while q:
        u = q.popleft()
        for v in adj.get(u, []):
            if v not in dist:
                d = dist[u] + 1
                if max_dist is not None and d > max_dist:
                    continue
                dist[v] = d
                by_dist.setdefault(d, []).append(v)
                q.append(v)

    return dist, by_dist


def bfs_within(atom, bonds, max_dist=None, silent=False):
    """BFS from atom; print distance matrix and return within-set."""
    dist, by_dist = _bfs_distmat([atom], bonds, max_dist)

    if not silent:
        print(f"Bond-distance matrix from atom {atom}:")
        for d in sorted(by_dist):
            print(f"  d={d} ({len(by_dist[d])}): {', '.join(map(str, sorted(by_dist[d])))}")

    return set(dist.keys())


def show(atom, bonds, angles, propers, impropers, flags, within_set=None):
    """Print topology entries. If within_set, keep only entries where all atoms are in it."""
    def keep(entry):
        return all(a in within_set for a in entry) if within_set else (atom in entry)

    if flags.bond:
        print("Bonds:")
        for a, b in sorted(bonds):
            if keep((a, b)):
                print(f"  {a}-{b}")
    if flags.angle:
        print("Angles:")
        for a, b, c in sorted(angles):
            if keep((a, b, c)):
                print(f"  {a}-{b}-{c}")
    if flags.dihedral:
        print("Proper dihedrals:")
        for a, b, c, d in sorted(propers):
            if keep((a, b, c, d)):
                print(f"  {a}-{b}-{c}-{d}")
        print("Improper dihedrals:")
        for a, b, c, d in sorted(impropers):
            if keep((a, b, c, d)):
                print(f"  {a}-{b}-{c}-{d}")


# ── template writer ──────────────────────────────────────────────

def _tp_type(entry, types):
    """Build type string from atom tuple and types lookup."""
    return '-'.join(types[a] for a in entry)


def write_template(fpath, atoms, idmap, types_lu, charges_lu, coords_lu,
                   mol_ids, bonds, angles, propers, impropers, title):
    """Write a LAMMPS molecule template file."""
    n = len(atoms)
    aset = set(atoms)

    def filt(entries, nreq):
        return sorted(e for e in entries if len({*e} & aset) == nreq)

    fbonds = filt(bonds, 2)
    fangles = filt(angles, 3)
    fpropers = filt(propers, 4)
    fimpropers = filt(impropers, 4)

    with open(fpath, 'w') as f:
        f.write(f"{title}\n\n")
        f.write(f"{n:>6} atoms\n")
        f.write(f"{len(fbonds):>6} bonds\n")
        f.write(f"{len(fangles):>6} angles\n")
        f.write(f"{len(fpropers):>6} dihedrals\n")
        f.write(f"{len(fimpropers):>6} impropers\n")

        # Coords
        f.write("\nCoords\n\n")
        for cid in atoms:
            tid = idmap[cid]
            x, y, z = coords_lu[cid]
            f.write(f"{tid:>6} {x:>15.6f} {y:>15.6f} {z:>15.6f}\n")

        # Types
        f.write("\nTypes\n\n")
        for cid in atoms:
            f.write(f"{idmap[cid]:>6} {types_lu[cid]}\n")

        # Charges
        f.write("\nCharges\n\n")
        for cid in atoms:
            f.write(f"{idmap[cid]:>6} {charges_lu[cid]:>12.6f}\n")

        # Molecules
        f.write("\nMolecules\n\n")
        for cid in atoms:
            f.write(f"{idmap[cid]:>6} {mol_ids[cid]}\n")

        # Bonds
        f.write("\nBonds\n\n")
        for i, (a, b) in enumerate(fbonds, 1):
            f.write(f"{i:>6} {_tp_type((a, b), types_lu):<12} {idmap[a]:>6} {idmap[b]:>6}\n")

        # Angles
        f.write("\nAngles\n\n")
        for i, (a, b, c) in enumerate(fangles, 1):
            f.write(f"{i:>6} {_tp_type((a, b, c), types_lu):<18} {idmap[a]:>6} {idmap[b]:>6} {idmap[c]:>6}\n")

        # Dihedrals
        f.write("\nDihedrals\n\n")
        for i, (a, b, c, d) in enumerate(fpropers, 1):
            f.write(f"{i:>6} {_tp_type((a, b, c, d), types_lu):<22} {idmap[a]:>6} {idmap[b]:>6} {idmap[c]:>6} {idmap[d]:>6}\n")

        # Impropers
        f.write("\nImpropers\n\n")
        for i, (a, b, c, d) in enumerate(fimpropers, 1):
            f.write(f"{i:>6} {_tp_type((a, b, c, d), types_lu):<22} {idmap[a]:>6} {idmap[b]:>6} {idmap[c]:>6} {idmap[d]:>6}\n")

        f.write("\n")


def gen_templates(mols, atom_list, individual, offsets, within_set, args):
    """Generate pre and post LAMMPS molecule template files."""
    # --- gather data ---
    atoms = sorted(within_set)
    idmap = {cid: i for i, cid in enumerate(atoms, 1)}  # combined -> template

    # Pre: types, charges, topology from individual ITPs
    pre_types, pre_charges, pre_coords, pre_mol = {}, {}, {}, {}
    pre_bonds, pre_angles, pre_propers, pre_impropers = set(), set(), set(), set()

    for i, mol in enumerate(mols):
        _, bonds, angles, propers, impropers, types, charges = individual[i]
        shift = offsets[i]
        gro = parse_gro(Path(f"{mol.strip()}.gro"))
        for aid in types:
            cid = aid + shift
            pre_types[cid] = types[aid]
            pre_charges[cid] = charges[aid]
            pre_coords[cid] = gro[aid]
            pre_mol[cid] = i + 1  # molecule 1 or 2
        for a, b in bonds:
            pre_bonds.add((a + shift, b + shift))
        for a, b, c in angles:
            pre_angles.add((a + shift, b + shift, c + shift))
        for a, b, c, d in propers:
            pre_propers.add((a + shift, b + shift, c + shift, d + shift))
        for a, b, c, d in impropers:
            pre_impropers.add((a + shift, b + shift, c + shift, d + shift))

    # Post: types, charges, topology from combined ITP
    combined_name = '-'.join(m.strip() for m in mols)
    _, cbonds, cangles, cpropers, cimpropers, ctypes, ccharges = parse_itp(
        Path(f"{combined_name}.itp"))

    # Coords for post: from combined GRO if available, else individual GROs
    cgro = Path(f"{combined_name}.gro")
    if cgro.exists():
        post_gro = parse_gro(cgro)
        post_coords = {cid: post_gro[cid] for cid in atoms if cid < len(post_gro)}
    else:
        post_coords = dict(pre_coords)

    post_mol = dict(pre_mol)  # same molecule assignment

    # --- write ---
    pre_title = f"Pre-reaction template: {'+'.join(mols)} atoms {'+'.join(map(str,atom_list))} within={args.within}"
    post_title = f"Post-reaction template: {combined_name} atoms {'+'.join(map(str,atom_list))} within={args.within}"

    write_template(f"{combined_name}.pre", atoms, idmap,
                   pre_types, pre_charges, pre_coords, pre_mol,
                   pre_bonds, pre_angles, pre_propers, pre_impropers, pre_title)

    write_template(f"{combined_name}.post", atoms, idmap,
                   ctypes, ccharges, post_coords, post_mol,
                   cbonds, cangles, cpropers, cimpropers, post_title)

    print(f"Generated {combined_name}.pre ({len(atoms)} atoms)")
    print(f"Generated {combined_name}.post ({len(atoms)} atoms)")


# ── main ─────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description="Query ITP topology for a given atom")
    p.add_argument('-mol', required=True, help='Molecule name(s), comma-separated for multi (reads <mol>.itp)')
    p.add_argument('-atom', required=True, help='Atom index(es), comma-separated, one per molecule')
    p.add_argument('-bond', action='store_true')
    p.add_argument('-angle', action='store_true')
    p.add_argument('-dihedral', action='store_true')
    p.add_argument('-within', type=int, nargs='?', const=None, default=-1,
                   help='Bond-distance cutoff; with multi-mol, uses combined topology')
    args = p.parse_args()

    mols = [m.strip() for m in args.mol.split(',')]
    atom_list = [int(x.strip()) for x in args.atom.split(',')]
    if len(atom_list) != len(mols):
        print(f"Error: got {len(mols)} mols but {len(atom_list)} atoms (must match)")
        return

    # --- Parse individual ITPs ---
    individual = []
    for mol in mols:
        fp = Path(f"{mol}.itp")
        if not fp.exists():
            print(f"Error: {fp} not found")
            return
        individual.append(parse_itp(fp))

    offsets = [0]
    for natoms, *_ in individual:
        offsets.append(offsets[-1] + natoms)

    # --- Single-molecule mode ---
    if len(mols) == 1:
        _, bonds, angles, propers, impropers, *_ = individual[0]
        within_set = None
        if args.within != -1:
            within_set = bfs_within(atom_list[0], bonds, args.within)
        show(atom_list[0], bonds, angles, propers, impropers, args, within_set)
        return

    # --- Multi-molecule mode ---
    combined_name = '-'.join(mols)
    cfp = Path(f"{combined_name}.itp")
    if not cfp.exists():
        print(f"Error: combined ITP {cfp} not found")
        return
    _, cbonds, cangles, cpropers, cimpropers, *_ = parse_itp(cfp)

    # Individual BFS, print with native numbering, then shift for union
    within_set = set()
    for i, (_, bonds, *_) in enumerate(individual):
        dist, by_dist = _bfs_distmat([atom_list[i]], bonds, args.within)
        s = set(dist.keys())
        print(f"Bond-distance matrix from {mols[i]} atom {atom_list[i]}:")
        for d in sorted(by_dist):
            print(f"  d={d} ({len(by_dist[d])}): {', '.join(map(str, sorted(by_dist[d])))}")
        within_set |= {x + offsets[i] for x in s}

    # Multi-source BFS on combined ITP
    starts = [atom_list[i] + offsets[i] for i in range(len(mols))]
    cdist, cby = _bfs_distmat(starts, cbonds, args.within)

    # Union
    final_set = within_set | set(cdist.keys())

    # Print combined distance matrix
    by_dist = {}
    for v in final_set:
        d = cdist.get(v, 1e9)
        by_dist.setdefault(d, []).append(v)

    label = ', '.join(f'{a}({m})' for a, m in zip(atom_list, mols))
    print(f"[Combined] bond-distance matrix from atoms {label}:")
    for d in sorted(by_dist):
        if args.within is not None and d > args.within:
            break
        atoms = sorted(by_dist[d])
        print(f"  d={d} ({len(atoms)}): {', '.join(map(str, atoms))}")

    # Show topology filtered by final_set
    show(atom_list[0], cbonds, cangles, cpropers, cimpropers, args, final_set)

    # Generate LAMMPS templates
    gen_templates(mols, atom_list, individual, offsets, final_set, args)


if __name__ == '__main__':
    main()
