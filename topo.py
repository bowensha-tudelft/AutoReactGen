#!/usr/bin/env python3
"""Query topology (bonds/angles/dihedrals) involving a given atom from a GROMACS ITP file.
Also generates LAMMPS reaction molecule templates (pre/post)."""

import argparse
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class AtomType:
    name: str
    atomic_number: int | None
    mass: float
    charge: float
    ptype: str
    sigma: float
    epsilon: float
    comment: str = ''


@dataclass(frozen=True)
class Atom:
    id: int
    type: str
    residue: int
    resname: str
    name: str
    cgnr: int
    charge: float
    mass: float | None = None
    comment: str = ''


@dataclass(frozen=True)
class Bond:
    ai: int
    aj: int
    funct: str
    params: tuple[float, ...] = field(default_factory=tuple)
    comment: str = ''

    @property
    def atoms(self):
        return (self.ai, self.aj)


@dataclass(frozen=True)
class Angle:
    ai: int
    aj: int
    ak: int
    funct: str
    params: tuple[float, ...] = field(default_factory=tuple)
    comment: str = ''

    @property
    def atoms(self):
        return (self.ai, self.aj, self.ak)


@dataclass(frozen=True)
class Dihedral:
    ai: int
    aj: int
    ak: int
    al: int
    funct: str
    params: tuple[float, ...] = field(default_factory=tuple)
    improper: bool = False
    comment: str = ''

    @property
    def atoms(self):
        return (self.ai, self.aj, self.ak, self.al)


@dataclass
class ItpTopology:
    filepath: Path
    atomtypes: dict[str, AtomType] = field(default_factory=dict)
    atoms: dict[int, Atom] = field(default_factory=dict)
    bonds: list[Bond] = field(default_factory=list)
    angles: list[Angle] = field(default_factory=list)
    propers: list[Dihedral] = field(default_factory=list)
    impropers: list[Dihedral] = field(default_factory=list)
    moleculetype: str | None = None
    nrexcl: int | None = None

    @property
    def natoms(self):
        return len(self.atoms)


@dataclass
class TopologyFile:
    filepath: Path
    defaults: dict[str, str] = field(default_factory=dict)
    atomtypes: dict[str, AtomType] = field(default_factory=dict)
    includes: list[Path] = field(default_factory=list)
    system: str | None = None
    molecules: list[tuple[str, int]] = field(default_factory=list)


def _split_data_comment(raw):
    data, sep, comment = raw.partition(';')
    return data.strip(), comment.strip() if sep else ''


def _float_tokens(tokens):
    vals = []
    for tok in tokens:
        try:
            vals.append(float(tok))
        except ValueError:
            break
    return tuple(vals)


def _int_or_none(tok):
    try:
        return int(tok)
    except ValueError:
        return None


def parse_itp_rich(filepath):
    """Parse a GROMACS ITP file and preserve force-field parameters.

    The older parse_itp() wrapper below intentionally keeps its compact return
    shape for the existing template/BFS code. New data/in writers should use
    this richer object so no parameter columns or duplicate dihedral terms are
    lost before TypeRegistry assignment.
    """
    topo = ItpTopology(Path(filepath))
    section = None
    section_improper = False

    with open(filepath) as f:
        for raw in f:
            data, comment = _split_data_comment(raw)
            if not data:
                continue
            if data.startswith('['):
                section = data.strip('[]').strip().split()[0].lower()
                section_improper = 'improper' in raw.lower()
                continue
            if data.startswith('#'):
                continue

            tok = data.split()
            _parse_common_topology_line(topo, section, section_improper, tok, comment)

    return topo


def parse_top(filepath):
    """Parse a GROMACS .top file system layer.

    This captures the parts needed for box-driven assembly: defaults,
    atomtypes, includes, system name, and [molecules] counts. Included ITP files
    are recorded but not recursively expanded here.
    """
    top = TopologyFile(Path(filepath))
    section = None

    with open(filepath) as f:
        for raw in f:
            data, comment = _split_data_comment(raw)
            if not data:
                continue
            if data.startswith('#include'):
                inc = data.split(maxsplit=1)[1].strip().strip('"<>')
                top.includes.append(Path(inc))
                continue
            if data.startswith('#'):
                continue
            if data.startswith('['):
                section = data.strip('[]').strip().split()[0].lower()
                continue

            tok = data.split()
            if section == 'defaults' and len(tok) >= 5:
                top.defaults = {
                    'nbfunc': tok[0],
                    'comb-rule': tok[1],
                    'gen-pairs': tok[2],
                    'fudgeLJ': tok[3],
                    'fudgeQQ': tok[4],
                }
            elif section == 'atomtypes' and len(tok) >= 7:
                name = tok[0]
                top.atomtypes[name] = AtomType(
                    name=name,
                    atomic_number=_int_or_none(tok[1]),
                    mass=float(tok[2]),
                    charge=float(tok[3]),
                    ptype=tok[4],
                    sigma=float(tok[5]),
                    epsilon=float(tok[6]),
                    comment=comment,
                )
            elif section == 'system' and top.system is None:
                top.system = data
            elif section == 'molecules' and len(tok) >= 2:
                top.molecules.append((tok[0], int(tok[1])))

    return top


def _parse_common_topology_line(topo, section, section_improper, tok, comment):
    if section == 'atomtypes' and len(tok) >= 7:
        name = tok[0]
        topo.atomtypes[name] = AtomType(
            name=name,
            atomic_number=_int_or_none(tok[1]),
            mass=float(tok[2]),
            charge=float(tok[3]),
            ptype=tok[4],
            sigma=float(tok[5]),
            epsilon=float(tok[6]),
            comment=comment,
        )
    elif section == 'moleculetype' and len(tok) >= 2:
        topo.moleculetype = tok[0]
        topo.nrexcl = int(tok[1])
    elif section == 'atoms' and len(tok) >= 7:
        aid = int(tok[0])
        topo.atoms[aid] = Atom(
            id=aid,
            type=tok[1],
            residue=int(tok[2]),
            resname=tok[3],
            name=tok[4],
            cgnr=int(tok[5]),
            charge=float(tok[6]),
            mass=float(tok[7]) if len(tok) > 7 else None,
            comment=comment,
        )
    elif section == 'bonds' and len(tok) >= 3:
        topo.bonds.append(Bond(
            ai=int(tok[0]),
            aj=int(tok[1]),
            funct=tok[2],
            params=_float_tokens(tok[3:]),
            comment=comment,
        ))
    elif section == 'angles' and len(tok) >= 4:
        topo.angles.append(Angle(
            ai=int(tok[0]),
            aj=int(tok[1]),
            ak=int(tok[2]),
            funct=tok[3],
            params=_float_tokens(tok[4:]),
            comment=comment,
        ))
    elif section == 'dihedrals' and len(tok) >= 5:
        funct = tok[4]
        improper = section_improper or funct == '4'
        d = Dihedral(
            ai=int(tok[0]),
            aj=int(tok[1]),
            ak=int(tok[2]),
            al=int(tok[3]),
            funct=funct,
            params=_float_tokens(tok[5:]),
            improper=improper,
            comment=comment,
        )
        (topo.impropers if improper else topo.propers).append(d)


def parse_itp(filepath):
    """Return (natoms, bonds, angles, propers, impropers, types, charges)."""
    topo = parse_itp_rich(filepath)
    bonds = {b.atoms for b in topo.bonds}
    angles = {a.atoms for a in topo.angles}
    propers = {d.atoms for d in topo.propers}
    impropers = {d.atoms for d in topo.impropers}
    types = {aid: atom.type for aid, atom in topo.atoms.items()}
    charges = {aid: atom.charge for aid, atom in topo.atoms.items()}
    return topo.natoms, bonds, angles, propers, impropers, types, charges


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


def write_map(fpath, title, initiator_ids, edge_ids, equivalences, delete_ids=()):
    """Write a LAMMPS fix bond/react map file."""
    with open(fpath, 'w') as f:
        f.write(f"{title}\n\n")
        f.write(f"{len(edge_ids):>6} edgeIDs\n")
        f.write(f"{len(equivalences):>6} equivalences\n")
        f.write(f"{len(delete_ids):>6} deleteIDs\n")

        f.write("\nInitiatorIDs\n\n")
        for tid in initiator_ids:
            f.write(f"{tid:>6}\n")

        f.write("\nEdgeIDs\n\n")
        for tid in edge_ids:
            f.write(f"{tid:>6}\n")

        f.write("\nEquivalences\n\n")
        for pre_id, post_id in equivalences:
            f.write(f"{pre_id:>6} {post_id:>6}\n")

        f.write("\nDeleteIDs\n\n")
        for tid in delete_ids:
            f.write(f"{tid:>6}\n")

        f.write("\n")


def gen_templates(mols, atom_list, individual, offsets, within_set, starts, cdist, args):
    """Generate pre/post LAMMPS molecule template files and reaction map."""
    # --- gather data ---
    atoms = sorted(within_set)
    idmap = {cid: i for i, cid in enumerate(atoms, 1)}  # combined -> template

    if len(starts) != 2:
        print("Error: .map generation currently requires exactly two initiator atoms")
        return

    initiator_ids = [idmap[cid] for cid in starts]
    selected_dist = {cid: d for cid, d in cdist.items() if cid in within_set}
    edge_depth = max(selected_dist.values()) if selected_dist else 0
    edge_ids = sorted(idmap[cid] for cid, d in selected_dist.items() if d == edge_depth)
    equivalences = [(i, i) for i in range(1, len(atoms) + 1)]
    delete_ids = []

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
    map_title = f"Reaction map: {combined_name} atoms {'+'.join(map(str,atom_list))} within={args.within}"

    write_template(f"{combined_name}.pre", atoms, idmap,
                   pre_types, pre_charges, pre_coords, pre_mol,
                   pre_bonds, pre_angles, pre_propers, pre_impropers, pre_title)

    write_template(f"{combined_name}.post", atoms, idmap,
                   ctypes, ccharges, post_coords, post_mol,
                   cbonds, cangles, cpropers, cimpropers, post_title)

    write_map(f"{combined_name}.map", map_title,
              initiator_ids, edge_ids, equivalences, delete_ids)

    print(f"Generated {combined_name}.pre ({len(atoms)} atoms)")
    print(f"Generated {combined_name}.post ({len(atoms)} atoms)")
    print(f"Generated {combined_name}.map ({len(equivalences)} equivalences, {len(edge_ids)} edgeIDs)")


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
    gen_templates(mols, atom_list, individual, offsets, final_set, starts, cdist, args)


if __name__ == '__main__':
    main()
