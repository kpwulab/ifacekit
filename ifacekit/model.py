"""Load a PDB or mmCIF file into flat numpy arrays and split it into units.

A *unit* is the object that interfaces are drawn between:
  * one polymer chain (protein, DNA or RNA)            -> "A"
  * one non-polymer residue (ligand, ion, cofactor)     -> "A:ATP501"
  * all branched residues of a chain (glycans)          -> "A:branched"
Waters are kept as atoms (for bridging analysis) but never form units.
"""
from dataclasses import dataclass, field

import gemmi
import numpy as np
from scipy.sparse import coo_matrix
from scipy.spatial import cKDTree

from . import chem


@dataclass
class Model:
    path: str
    fmt: str
    model_index: int
    # per atom arrays
    xyz: np.ndarray
    element: np.ndarray          # str
    atom_name: np.ndarray        # str
    res_idx: np.ndarray          # int index into residues
    unit_idx: np.ndarray         # int index into units, -1 for water
    radius: np.ndarray
    donor: np.ndarray
    acceptor: np.ndarray
    charge: np.ndarray
    ambiguous: np.ndarray
    apolar: np.ndarray           # apolar carbon / sulfur
    is_metal: np.ndarray
    is_water: np.ndarray
    bfactor: np.ndarray
    # tables
    residues: list = field(default_factory=list)   # dicts
    units: list = field(default_factory=list)      # dicts
    bonds: np.ndarray = None                       # (n,2) covalent bonds
    occupancy: np.ndarray = None
    flags: dict = field(default_factory=dict)
    alt_extra: np.ndarray = None                    # True for atoms of the loaded alternate conformer that the first-conformer model does not have

    @property
    def n_atoms(self):
        return len(self.xyz)


def residue_id(chain, resnum, icode, resname):
    ic = icode.strip() if icode and icode.strip() else ""
    return f"{chain}:{resnum}{ic}:{resname}"


def _has_backbone(res):
    """True if a non-standard residue carries a peptide or nucleotide backbone.

    PDB files carry no explicit entity table, so an ion or ligand that
    follows a chain can be mistaken for a polymer residue; the backbone atom
    names settle it.
    """
    names = {chem.norm_atom_name(a.name) for a in res}
    protein = {"N", "CA", "C"} <= names
    nucleic = "P" in names and ("O3'" in names or "C3'" in names or "C4'" in names)
    return protein or nucleic


def _classify(res, info):
    """Return one of polymer / branched / ligand / water."""
    if res.name in chem.WATER_NAMES or (info is not None and info.is_water()):
        return "water"
    et = res.entity_type
    if et == gemmi.EntityType.Water:
        return "water"
    if et == gemmi.EntityType.Branched:
        return "branched"
    if info is not None and (info.is_amino_acid() or info.is_nucleic_acid()):
        # Standard-looking residues in a ligand chain (free amino acid) stay ligands
        if et == gemmi.EntityType.NonPolymer:
            return "ligand"
        return "polymer"
    if et == gemmi.EntityType.Polymer and _has_backbone(res):
        return "polymer"          # modified residue inside a polymer
    return "ligand"


def build_bonds(xyz, element):
    """Covalent bonds from interatomic distances (metals are not bonded)."""
    n = len(xyz)
    cov = np.array([chem.covalent_radius(e) for e in element])
    metal = np.array([gemmi.Element(e).is_metal for e in element])
    tree = cKDTree(xyz)
    pairs = tree.query_pairs(r=2.8, output_type="ndarray")
    if len(pairs) == 0:
        return np.zeros((0, 2), dtype=int)
    i, j = pairs[:, 0], pairs[:, 1]
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    ok = (d < cov[i] + cov[j] + chem.BOND_TOLERANCE) & ~metal[i] & ~metal[j]
    return pairs[ok]


def load(path, model_index=0, select=None, merge_by_chain=False, unit_chain=None, conformer=None):
    """Read a structure.

    select          optional callable(chain_name, kind, residue_id) -> bool;
                    residues for which it returns False are skipped (kind is
                    one of polymer / branched / ligand / water)
    merge_by_chain  one unit per chain (polymer, ligands and glycans of that
                    chain together) instead of one unit per ligand residue
    unit_chain      optional callable(chain_name, kind, residue_id) -> chain
                    name of the unit a residue joins when merge_by_chain is
                    set (lets a ligand or glycan that carries its own chain
                    label join the polymer chain it belongs to)
    conformer       None: the first conformer of every atom (default). A letter: the atoms of that conformer plus
                    the atoms shared by all conformers; model.alt_extra marks the atoms that the default model lacks
    """
    st = gemmi.read_structure(str(path))
    fmt = "mmcif" if str(st.input_format).lower().endswith("mmcif") or str(path).lower().endswith((".cif", ".mmcif")) else "pdb"
    st.setup_entities()
    flags = {"n_models_in_file": len(st)}
    if model_index >= len(st):
        raise ValueError(f"model {model_index} not in file ({len(st)} models)")
    # drop other models, keep one
    model = st[model_index]
    flags["hydrogens_removed"] = sum(
        1 for ch in model for r in ch for a in r if a.element.is_hydrogen)
    n_heavy_all = sum(1 for ch in model for r in ch for a in r
                      if not a.element.is_hydrogen)
    # Keep the first conformer of every atom, drop hydrogens
    st2 = st.clone()
    st2.remove_hydrogens()
    first_alt, letters = {}, set()
    for ch in st2[model_index]:
        for res in ch:
            fa = next((a.altloc for a in res if a.altloc != "\0"), None)
            if fa is not None:
                first_alt[(ch.name, res.seqid.num, res.seqid.icode)] = fa
                letters.update(a.altloc for a in res if a.altloc != "\0")
    flags["altloc_letters"] = sorted(letters)
    if conformer is None:
        st2.remove_alternative_conformations()
    else:
        for ch in st2[model_index]:
            for res in ch:
                for k in range(len(res) - 1, -1, -1):
                    if res[k].altloc not in ("\0", conformer):
                        del res[k]
    model = st2[model_index]
    flags["altloc_atoms_dropped"] = n_heavy_all - sum(
        1 for ch in model for r in ch for a in r)

    xyz, element, aname, res_idx, bfac, occ = [], [], [], [], [], []
    residues = []
    units = []
    unit_lookup = {}
    atom_unit = []
    alt_extra = []
    for ch in model:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            kind = _classify(res, info)
            rid = residue_id(ch.name, res.seqid.num, res.seqid.icode, res.name)
            if select is not None and not select(ch.name, kind, rid):
                continue
            uc = unit_chain(ch.name, kind, rid) if (unit_chain and merge_by_chain) else ch.name
            if kind == "water":
                ukey = None
            elif merge_by_chain:
                ukey = (uc, "chain")
            elif kind == "polymer":
                ukey = (ch.name, "polymer")
            elif kind == "branched":
                ukey = (ch.name, "branched")
            else:
                ukey = (ch.name, f"{res.name}{res.seqid.num}{res.seqid.icode.strip()}")
            if ukey is not None and ukey not in unit_lookup:
                if merge_by_chain:
                    uid, ukind = uc, "chain"
                elif kind == "polymer":
                    uid, ukind = ch.name, "polymer"
                elif kind == "branched":
                    uid, ukind = f"{ch.name}:branched", "branched"
                else:
                    uid, ukind = f"{ch.name}:{ukey[1]}", "ligand"
                unit_lookup[ukey] = len(units)
                units.append({"id": uid, "kind": ukind, "chain": uc if merge_by_chain else ch.name,
                              "polymer_type": None, "n_residues": 0,
                              "n_atoms": 0, "nonstandard": [], "ligands": []})
            uidx = unit_lookup[ukey] if ukey is not None else -1
            is_std = bool(info is not None and info.is_standard())
            rec = {
                "id": rid, "chain": ch.name, "resnum": res.seqid.num,
                "icode": res.seqid.icode.strip(), "resname": res.name,
                "kind": kind, "unit": uidx,
                "protein": bool(kind == "polymer" and info is not None and info.is_amino_acid()),
                "nucleic": bool(kind == "polymer" and info is not None and info.is_nucleic_acid()),
                "standard": is_std,
            }
            residues.append(rec)
            ridx = len(residues) - 1
            if uidx >= 0:
                u = units[uidx]
                u["n_residues"] += 1
                if kind == "polymer" and not is_std:
                    u["nonstandard"].append(rid)
                if kind in ("ligand", "branched") and merge_by_chain:
                    u["ligands"].append(rid)
                if kind == "polymer" and u["polymer_type"] is None and info is not None:
                    if info.is_amino_acid():
                        u["polymer_type"] = "protein"
                    elif info.is_nucleic_acid():
                        u["polymer_type"] = "rna" if "U" in res.name else "dna" if res.name.startswith("D") else "nucleic"
                elif kind == "polymer" and info is not None and info.is_nucleic_acid():
                    # RNA if any residue is not deoxy
                    if not res.name.startswith("D") and u["polymer_type"] == "dna":
                        u["polymer_type"] = "rna"
            for a in res:
                xyz.append([a.pos.x, a.pos.y, a.pos.z])
                element.append(a.element.name.upper())
                aname.append(chem.norm_atom_name(a.name))
                res_idx.append(ridx)
                atom_unit.append(uidx)
                bfac.append(a.b_iso)
                occ.append(a.occ)
                alt_extra.append(conformer is not None and a.altloc == conformer
                                 and first_alt.get((ch.name, res.seqid.num, res.seqid.icode)) != conformer)
                if uidx >= 0:
                    units[uidx]["n_atoms"] += 1

    if merge_by_chain:
        for u in units:
            u["kind"] = "polymer" if u["polymer_type"] else "ligand"

    xyz = np.array(xyz, dtype=float).reshape(-1, 3)
    element = np.array(element)
    aname = np.array(aname)
    res_idx = np.array(res_idx, dtype=int)
    unit_idx = np.array(atom_unit, dtype=int)

    n = len(xyz)
    radius = np.array([chem.vdw_radius(e) for e in element])
    is_metal = np.array([gemmi.Element(e).is_metal for e in element])
    is_water = np.array([residues[r]["kind"] == "water" for r in res_idx]) if n else np.zeros(0, bool)

    donor = np.zeros(n, bool)
    acceptor = np.zeros(n, bool)
    charge = np.zeros(n, int)
    ambiguous = np.zeros(n, bool)
    for k in range(n):
        r = residues[res_idx[k]]
        d, a, c, amb = chem.atom_role(r["resname"], aname[k], element[k],
                                      r["protein"], r["nucleic"], bool(is_water[k]))
        donor[k], acceptor[k], charge[k], ambiguous[k] = d, a, c, amb
    # metals are positively charged but take no part in donor/acceptor logic
    charge[is_metal] = 1

    bonds = build_bonds(xyz, element)
    # apolar carbon/sulfur: not bonded to N or O (S and halogens treated as apolar)
    has_polar_nb = np.zeros(n, bool)
    polar_el = (element == "N") | (element == "O")
    if len(bonds):
        b0, b1 = bonds[:, 0], bonds[:, 1]
        has_polar_nb[b0[polar_el[b1]]] = True
        has_polar_nb[b1[polar_el[b0]]] = True
    apolar = ((element == "C") & ~has_polar_nb) | (element == "S") | np.isin(
        element, ["F", "CL", "BR", "I"])

    try:
        with open(path, "r", errors="ignore") as fh:
            text = fh.read(400000)
        low = text.lower()
        pred_words = [w for w in ("alphafold", "boltz", "chai-1", "chai_", "colabfold", "openfold", "esmfold", "rosettafold", "_ma_qa_metric") if w in low]
        exp_marks = [w for w in ("_exptl.method", "_refine.", "_reflns.", "_em_3d_reconstruction", "_pdbx_nmr", "expdta", "remark 210") if w in low]
        if pred_words:
            flags["predicted_model"], flags["model_source"] = True, "prediction keyword in file: " + pred_words[0]
        elif exp_marks:
            flags["predicted_model"], flags["model_source"] = False, "experimental record in file: " + exp_marks[0]
        else:
            flags["predicted_model"], flags["model_source"] = False, "no experimental record, assumed experimental"
    except OSError:
        flags["predicted_model"], flags["model_source"] = False, "file not readable"
    bf = np.array(bfac, dtype=float)
    flags["bfactor_scaled"] = False
    if bf.size and bf.max() <= 1.0 and bf.min() >= 0.0 and bf.max() > 0.0 and not flags["model_source"].startswith("experimental record"):
        # confidence stored as a fraction (Boltz, Chai and similar): bring it to the 0-100 pLDDT scale.
        # Not for a file with an experimental record, and not when every B-factor is 0 (NMR ensembles carry no B-factors)
        bf = bf * 100.0
        flags["bfactor_scaled"] = True
        flags["predicted_model"] = True
        flags["model_source"] = "confidence stored as a 0 to 1 fraction in the B-factor column"
    bfac = list(bf)
    flags["nonstandard_residues"] = sorted(
        {r["resname"] for r in residues if r["kind"] == "polymer" and not r["standard"]})
    flags["n_ligand_residues"] = sum(1 for r in residues if r["kind"] == "ligand")
    flags["n_water_residues"] = sum(1 for r in residues if r["kind"] == "water")
    flags["unknown_ligand_roles"] = bool(
        ambiguous[~is_water].any()) if n else False

    return Model(path=str(path), fmt=fmt, model_index=model_index, xyz=xyz,
                 element=element, atom_name=aname, res_idx=res_idx,
                 unit_idx=unit_idx, radius=radius, donor=donor,
                 acceptor=acceptor, charge=charge, ambiguous=ambiguous,
                 apolar=apolar, is_metal=is_metal, is_water=is_water,
                 bfactor=np.array(bfac), residues=residues, units=units,
                 bonds=bonds, flags=flags, occupancy=np.array(occ),
                 alt_extra=np.array(alt_extra, dtype=bool))


def within_n_bonds(bonds, n_atoms, depth=3):
    """Set of pair keys (i*n+j, i<j) that are 1-2, 1-3 or 1-4 neighbours."""
    if len(bonds) == 0:
        return np.zeros(0, dtype=np.int64)
    rows = np.concatenate([bonds[:, 0], bonds[:, 1]])
    cols = np.concatenate([bonds[:, 1], bonds[:, 0]])
    A = coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)),
                   shape=(n_atoms, n_atoms)).tocsr()
    A.data[:] = 1
    reach = A.copy()
    power = A.copy()
    for _ in range(depth - 1):
        power = power @ A
        power.data[:] = 1
        reach = reach + power
    reach = reach.tocoo()
    m = reach.row < reach.col
    keys = reach.row[m].astype(np.int64) * n_atoms + reach.col[m].astype(np.int64)
    return np.unique(keys)
