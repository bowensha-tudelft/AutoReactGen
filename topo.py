#!/usr/bin/env python3
"""Query topology (bonds/angles/dihedrals) involving a given atom from a GROMACS ITP file."""

import argparse
from pathlib import Path

def parse_itp(filepath):
    """Return (bonds, angles, proper_dihedrals, improper_dihedrals) as sets of tuples."""
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
            if section == 'bonds':
                bonds.add((int(tok[0]), int(tok[1])))
            elif section == 'angles':
                angles.add((int(tok[0]), int(tok[1]), int(tok[2])))
            elif section == 'dihedrals':
                d = (int(tok[0]), int(tok[1]), int(tok[2]), int(tok[3]))
                (impropers if improper else propers).add(d)

    return bonds, angles, propers, impropers


def build_adj(bonds):
    """Build adjacency dict from bond pairs."""
    adj = {}
    for a, b in bonds:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    return adj


def show_distmat(atom, bonds, max_dist=None):
    """BFS from atom over bond graph, print atoms grouped by bond-distance."""
    adj = build_adj(bonds)
    if atom not in adj:
        print(f"Atom {atom} has no bonds.")
        return

    from collections import deque
    dist = {atom: 0}
    q = deque([atom])
    by_dist = {0: [atom]}

    while q:
        u = q.popleft()
        for v in adj.get(u, []):
            if v not in dist:
                dist[v] = dist[u] + 1
                if max_dist is not None and dist[v] > max_dist:
                    continue
                by_dist.setdefault(dist[v], []).append(v)
                q.append(v)

    print(f"Bond-distance matrix from atom {atom}:")
    for d in sorted(by_dist):
        if max_dist is not None and d > max_dist:
            break
        atoms = sorted(by_dist[d])
        print(f"  d={d} ({len(atoms)}): {', '.join(map(str, atoms))}")


def show(atom, bonds, angles, propers, impropers, flags):
    if flags.bond:
        print("Bonds:")
        for a, b in sorted(bonds):
            if atom in (a, b):
                print(f"  {a}-{b}")
    if flags.angle:
        print("Angles:")
        for a, b, c in sorted(angles):
            if atom in (a, b, c):
                print(f"  {a}-{b}-{c}")
    if flags.dihedral:
        print("Proper dihedrals:")
        for a, b, c, d in sorted(propers):
            if atom in (a, b, c, d):
                print(f"  {a}-{b}-{c}-{d}")
        print("Improper dihedrals:")
        for a, b, c, d in sorted(impropers):
            if atom in (a, b, c, d):
                print(f"  {a}-{b}-{c}-{d}")


def main():
    p = argparse.ArgumentParser(description="Query ITP topology for a given atom")
    p.add_argument('-mol', required=True, help='Molecule name (reads <mol>.itp)')
    p.add_argument('-atom2show', type=int, required=True, help='Atom index to query')
    p.add_argument('-bond', action='store_true')
    p.add_argument('-angle', action='store_true')
    p.add_argument('-dihedral', action='store_true')
    p.add_argument('-within', type=int, nargs='?', const=None, default=-1,
                   help='Show bond-distance matrix from atom, optional max depth')
    args = p.parse_args()

    itp = Path(f"{args.mol}.itp")
    if not itp.exists():
        print(f"Error: {itp} not found")
        return

    bonds, angles, propers, impropers = parse_itp(itp)

    if args.within != -1:
        show_distmat(args.atom2show, bonds, args.within)

    show(args.atom2show, bonds, angles, propers, impropers, args)


if __name__ == '__main__':
    main()
