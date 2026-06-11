#!/usr/bin/env python3
"""Query topology (bonds/angles/dihedrals) involving a given atom from a GROMACS ITP file."""

import argparse
from collections import deque
from pathlib import Path


def parse_itp(filepath):
    """Return (natoms, bonds, angles, proper_dihedrals, improper_dihedrals)."""
    natoms = 0
    bonds, angles = set(), set()
    propers, impropers = set(), set()
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
                natoms += 1
            elif section == 'bonds':
                bonds.add((int(tok[0]), int(tok[1])))
            elif section == 'angles':
                angles.add((int(tok[0]), int(tok[1]), int(tok[2])))
            elif section == 'dihedrals':
                d = (int(tok[0]), int(tok[1]), int(tok[2]), int(tok[3]))
                (impropers if improper else propers).add(d)

    return natoms, bonds, angles, propers, impropers


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

    mols = args.mol.split(',')
    atom_list = [int(x.strip()) for x in args.atom.split(',')]
    if len(atom_list) != len(mols):
        print(f"Error: got {len(mols)} mols but {len(atom_list)} atoms (must match)")
        return

    # --- Parse individual ITPs ---
    individual = []
    for mol in mols:
        fp = Path(f"{mol.strip()}.itp")
        if not fp.exists():
            print(f"Error: {fp} not found")
            return
        individual.append(parse_itp(fp))

    offsets = [0]
    for natoms, _, _, _, _ in individual:
        offsets.append(offsets[-1] + natoms)

    # --- Single-molecule mode ---
    if len(mols) == 1:
        _, bonds, angles, propers, impropers = individual[0]
        within_set = None
        if args.within != -1:
            within_set = bfs_within(atom_list[0], bonds, args.within)
        show(atom_list[0], bonds, angles, propers, impropers, args, within_set)
        return

    # --- Multi-molecule mode ---
    combined_name = '-'.join(m.strip() for m in mols)
    cfp = Path(f"{combined_name}.itp")
    if not cfp.exists():
        print(f"Error: combined ITP {cfp} not found")
        return
    _, cbonds, cangles, cpropers, cimpropers = parse_itp(cfp)

    # Individual BFS, print with native numbering, then shift for union
    within_set = set()
    for i, (_, bonds, _, _, _) in enumerate(individual):
        dist, by_dist = _bfs_distmat([atom_list[i]], bonds, args.within)
        s = set(dist.keys())
        print(f"Bond-distance matrix from {mols[i].strip()} atom {atom_list[i]}:")
        for d in sorted(by_dist):
            print(f"  d={d} ({len(by_dist[d])}): {', '.join(map(str, sorted(by_dist[d])))}")
        within_set |= {x + offsets[i] for x in s}

    # Multi-source BFS on combined ITP
    starts = [atom_list[i] + offsets[i] for i in range(len(mols))]
    cdist, cby = _bfs_distmat(starts, cbonds, args.within)

    # Union
    cmb_set = set(cdist.keys())
    final_set = within_set | cmb_set

    # Recompute by_dist over combined reachable, using combined distances
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


if __name__ == '__main__':
    main()
