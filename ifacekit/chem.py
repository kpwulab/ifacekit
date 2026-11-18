"""Chemistry tables: radii and donor/acceptor/charge assignment.

Everything here is a plain lookup so that the engine does not depend on any
restraint library. Atoms that are not in the tables (ligands, modified
residues) fall back to element based rules and are marked ``ambiguous``.
"""

# Bondi van der Waals radii (heavy atoms only; hydrogens are not used).
VDW = {
    "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "SE": 1.90,
    "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98,
}
VDW_DEFAULT = 1.80

# Radii used only for building the covalent bond graph.
COVALENT = {
    "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07, "SE": 1.20,
    "F": 0.57, "CL": 1.02, "BR": 1.20, "I": 1.39, "H": 0.31,
}
COVALENT_DEFAULT = 1.10
BOND_TOLERANCE = 0.45

# Atom roles: D donor, A acceptor, B both (hydroxyl, imidazole), + / - charge.
_PROTEIN_SIDECHAIN = {
    "ARG": {"NE": "D+", "NH1": "D+", "NH2": "D+"},
    "ASN": {"OD1": "A", "ND2": "D"},
    "ASP": {"OD1": "A-", "OD2": "A-"},
    "GLN": {"OE1": "A", "NE2": "D"},
    "GLU": {"OE1": "A-", "OE2": "A-"},
    "HIS": {"ND1": "B", "NE2": "B"},
    "LYS": {"NZ": "D+"},
    "SER": {"OG": "B"},
    "THR": {"OG1": "B"},
    "TRP": {"NE1": "D"},
    "TYR": {"OH": "B"},
}
_PROTEIN_BACKBONE = {"N": "D", "O": "A", "OXT": "A-"}

_NUC_BACKBONE = {
    "OP1": "A-", "OP2": "A-", "OP3": "A-", "O1P": "A-", "O2P": "A-",
    "O5'": "A", "O3'": "A", "O4'": "A", "O2'": "B",
}
_NUC_BASE = {
    "A": {"N6": "D", "N1": "A", "N3": "A", "N7": "A"},
    "G": {"N1": "D", "N2": "D", "O6": "A", "N3": "A", "N7": "A"},
    "C": {"N4": "D", "O2": "A", "N3": "A"},
    "U": {"N3": "D", "O2": "A", "O4": "A"},
    "T": {"N3": "D", "O2": "A", "O4": "A"},
}
_NUC_PARENT = {
    "A": "A", "DA": "A", "G": "G", "DG": "G", "C": "C", "DC": "C",
    "U": "U", "DU": "U", "DT": "T", "T": "T",
}

WATER_NAMES = {"HOH", "WAT", "H2O", "DOD", "TIP", "TIP3", "SOL"}

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
STANDARD_NUC = set(_NUC_PARENT)


def norm_atom_name(name):
    return name.replace("*", "'")


def vdw_radius(element):
    return VDW.get(element.upper(), VDW_DEFAULT)


def covalent_radius(element):
    return COVALENT.get(element.upper(), COVALENT_DEFAULT)


def atom_role(resname, atomname, element, in_protein, in_nucleic, is_water):
    """Return (donor, acceptor, charge, ambiguous) for one atom.

    ``ambiguous`` means the role is a guess from the element only (no
    hydrogens and no restraint library), so mismatch rules are not applied.
    """
    el = element.upper()
    if is_water:
        return True, True, 0, True
    name = norm_atom_name(atomname)
    role = None
    if in_protein:
        role = _PROTEIN_SIDECHAIN.get(resname, {}).get(name)
        if role is None and name in _PROTEIN_BACKBONE:
            if not (name == "N" and resname == "PRO"):
                role = _PROTEIN_BACKBONE[name]
            else:
                return False, False, 0, False
    elif in_nucleic:
        role = _NUC_BACKBONE.get(name)
        if role is None:
            parent = _NUC_PARENT.get(resname)
            if parent:
                role = _NUC_BASE[parent].get(name)
    if role is not None:
        donor = role[0] in "DB"
        acceptor = role[0] in "AB"
        charge = 1 if role.endswith("+") else (-1 if role.endswith("-") else 0)
        return donor, acceptor, charge, False
    if el in ("N", "O"):
        # In a standard residue every N/O that can hydrogen bond is listed
        # above, so an unlisted one is inert (Pro N, glycosidic N9/N1).
        if (in_protein and resname in STANDARD_AA) or (in_nucleic and resname in _NUC_PARENT):
            return False, False, 0, False
        # Anything else (ligand, modified residue) is unknown chemistry.
        return True, True, 0, True
    return False, False, 0, False


# Ring atom names for the pi analysis (standard residues only).
RINGS = {
    "PHE": [("CG", "CD1", "CE1", "CZ", "CE2", "CD2")],
    "TYR": [("CG", "CD1", "CE1", "CZ", "CE2", "CD2")],
    "HIS": [("CG", "ND1", "CE1", "NE2", "CD2")],
    "TRP": [("CG", "CD1", "NE1", "CE2", "CD2"),
            ("CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3")],
}
_PURINE = [("N9", "C4", "C5", "N7", "C8"), ("C4", "C5", "C6", "N1", "C2", "N3")]
_PYRIMIDINE = [("N1", "C2", "N3", "C4", "C5", "C6")]
for _n in ("A", "DA", "G", "DG"):
    RINGS[_n] = _PURINE
for _n in ("C", "DC", "U", "DU", "T", "DT"):
    RINGS[_n] = _PYRIMIDINE

# Cation atoms for cation-pi.
CATIONS = {("LYS", "NZ"), ("ARG", "CZ")}
