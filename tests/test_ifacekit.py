"""Tests. Real structures are public PDB entries that tests/pdbdata.py downloads from RCSB on first use (cached in tests/.pdb_cache);
tests that need them are skipped when there is no internet."""
import json
import os
from pathlib import Path

import numpy as np
import pytest

import ifacekit
from ifacekit import model as mdl
from ifacekit.sasa import atom_sasa
from pdbdata import need, prediction_like        # public PDB entries, downloaded on first use (tests/pdbdata.py)


# ---------------------------------------------------------------- helpers
def pdb_text(atoms):
    """atoms: list of (serial-less) tuples (name, resname, chain, resnum, element, x, y, z, hetatm)."""
    lines = []
    for k, (name, resname, chain, resnum, el, x, y, z, het) in enumerate(atoms, 1):
        rec = "HETATM" if het else "ATOM  "
        lines.append(f"{rec}{k:5d} {name:<4s} {resname:>3s} {chain}{resnum:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {el:>2s}")
    lines.append("END")
    return "\n".join(lines) + "\n"


def pair(tmp_path, a, b, d, chain_a="A", chain_b="B"):
    """Two single-atom residues at distance d. a/b = (atom, resname, element)."""
    atoms = [(a[0], a[1], chain_a, 1, a[2], 0.0, 0.0, 0.0, False),
             (b[0], b[1], chain_b, 1, b[2], d, 0.0, 0.0, False)]
    p = tmp_path / "pair.pdb"
    p.write_text(pdb_text(atoms))
    return ifacekit.analyze(str(p), min_area=0.0, n_points=200)


def only(res):
    assert len(res["contacts"]) == 1, res["contacts"]
    return res["contacts"][0]


# ---------------------------------------------------------------- SASA
def test_single_atom_sasa_is_analytic():
    xyz = np.zeros((1, 3))
    r = np.array([1.7])
    got = atom_sasa(xyz, r, n_points=960)[0]
    assert got == pytest.approx(4 * np.pi * 3.1 ** 2, rel=1e-9)


def test_two_touching_atoms_lose_area_symmetrically():
    xyz = np.array([[0, 0, 0], [3.0, 0, 0]], float)
    r = np.array([1.7, 1.7])
    s = atom_sasa(xyz, r, n_points=960)
    # a finite point set is not rotationally symmetric, so allow ~2% sampling noise
    assert s[0] == pytest.approx(s[1], rel=0.02)
    assert s[0] < 4 * np.pi * 3.1 ** 2


def test_sasa_matches_freesasa_same_radii():
    fs = pytest.importorskip("freesasa")
    m = mdl.load(need("1A8O.pdb"))
    keep = np.where(~m.is_water)[0]
    ours = atom_sasa(m.xyz[keep], m.radius[keep], n_points=960).sum()
    s = fs.Structure()
    for k in keep:
        r = m.residues[m.res_idx[k]]
        s.addAtom(f"{m.atom_name[k]:<4s}", r["resname"], r["resnum"], r["chain"], *map(float, m.xyz[k]))
    for n, k in enumerate(keep):
        s.setRadius(n, float(m.radius[k]))
    ref = fs.calc(s, fs.Parameters({"algorithm": fs.ShrakeRupley, "probe-radius": 1.4,
                                    "n-points": 960})).totalArea()
    assert ours == pytest.approx(ref, rel=0.01)


# ---------------------------------------------------------------- contact rules
def test_salt_bridge_favorable(tmp_path):
    c = only(pair(tmp_path, ("OD1", "ASP", "O"), ("NZ", "LYS", "N"), 2.85))
    assert (c["type"], c["verdict"], c["hbond"]) == ("salt_bridge", "favorable", True)


def test_salt_bridge_long_range(tmp_path):
    c = only(pair(tmp_path, ("OD1", "ASP", "O"), ("NZ", "LYS", "N"), 3.8))
    assert (c["type"], c["verdict"], c["hbond"]) == ("salt_bridge", "favorable", False)


def test_impossibly_short_polar_pair(tmp_path):
    c = only(pair(tmp_path, ("OD1", "ASP", "O"), ("NZ", "LYS", "N"), 2.2))
    assert c["verdict"] == "unfavorable"


def test_hbond_ser_asn(tmp_path):
    c = only(pair(tmp_path, ("OG", "SER", "O"), ("ND2", "ASN", "N"), 2.9))
    assert (c["type"], c["verdict"]) == ("hbond", "favorable")


def test_long_hbond_is_ok_not_favorable(tmp_path):
    c = only(pair(tmp_path, ("OG", "SER", "O"), ("ND2", "ASN", "N"), 3.45))
    assert (c["type"], c["verdict"], c["reason"]) == ("hbond", "ok", "long_hbond")


@pytest.mark.parametrize("d,verdict", [(3.9, "ok"), (3.1, "borderline"), (2.9, "borderline"), (2.7, "unfavorable")])
def test_apolar_carbon_overlap(tmp_path, d, verdict):
    c = only(pair(tmp_path, ("CD1", "LEU", "C"), ("CD1", "ILE", "C"), d))
    assert c["verdict"] == verdict
    assert c["type"] in ("hydrophobic", "vdw")


@pytest.mark.parametrize("d,verdict", [(2.4, "unfavorable"), (2.7, "borderline"), (3.0, "borderline")])
def test_acceptor_acceptor_mismatch(tmp_path, d, verdict):
    c = only(pair(tmp_path, ("O", "GLY", "O"), ("O", "GLY", "O"), d))
    assert (c["type"], c["verdict"]) == ("polar_mismatch", verdict)


def test_backbone_n_to_o_is_hbond(tmp_path):
    c = only(pair(tmp_path, ("N", "GLY", "N"), ("O", "GLY", "O"), 2.95))
    assert (c["type"], c["verdict"]) == ("hbond", "favorable")


def test_pro_n_is_not_a_donor(tmp_path):
    c = pair(tmp_path, ("N", "PRO", "N"), ("O", "GLY", "O"), 2.95)
    assert not any(x["type"] == "hbond" for x in c["contacts"])


def test_like_charge_flagged(tmp_path):
    c = only(pair(tmp_path, ("NZ", "LYS", "N"), ("NH1", "ARG", "N"), 3.5))
    assert (c["type"], c["verdict"]) == ("polar_mismatch", "borderline") or c["type"] == "like_charge"


def test_bonded_atoms_are_not_compared(tmp_path):
    atoms = [("CA", "ALA", "A", 1, "C", 0.0, 0.0, 0.0, False),
             ("CB", "ALA", "A", 1, "C", 1.53, 0.0, 0.0, False)]
    p = tmp_path / "b.pdb"
    p.write_text(pdb_text(atoms))
    res = ifacekit.analyze(str(p), n_points=100)
    assert res["contacts"] == []


def test_metal_coordination(tmp_path):
    atoms = [("NE2", "HIS", "A", 1, "N", 0.0, 0.0, 0.0, False),
             ("ZN", "ZN", "A", 2, "ZN", 2.05, 0.0, 0.0, True)]
    p = tmp_path / "m.pdb"
    p.write_text(pdb_text(atoms))
    res = ifacekit.analyze(str(p), min_area=0.0, n_points=100)
    assert [c["type"] for c in res["contacts"]] == ["metal"]
    assert res["contacts"][0]["verdict"] == "ok"


def test_bridging_water(tmp_path):
    atoms = [("OD1", "ASP", "A", 1, "O", 0.0, 0.0, 0.0, False),
             ("O", "HOH", "A", 2, "O", 2.8, 0.0, 0.0, True),
             ("ND2", "ASN", "B", 1, "N", 5.6, 0.0, 0.0, False)]
    p = tmp_path / "w.pdb"
    p.write_text(pdb_text(atoms))
    res = ifacekit.analyze(str(p), min_area=0.0, n_points=100)
    # no direct interface area needed; the water must be recognised as bridging
    # when an interface exists, so add a direct contact too
    atoms.append(("CB", "ASP", "A", 1, "C", 0.0, 3.2, 0.0, False))
    atoms.append(("CB", "ASN", "B", 1, "C", 5.6, 3.2, 0.0, False))
    p.write_text(pdb_text(atoms))
    res = ifacekit.analyze(str(p), min_area=0.0, n_points=100)
    waters = [c for c in res["contacts"] if "HOH" in c["atom1"] + c["atom2"]]
    assert len(waters) == 2 and all(c["type"] == "hbond" for c in waters)


# ---------------------------------------------------------------- real files
def _key(res):
    ifs = sorted((tuple(sorted((i["unit_a"], i["unit_b"]))), i["area"]) for i in res["interfaces"])
    cts = sorted((tuple(sorted((c["atom1"], c["atom2"]))), c["type"], c["verdict"])
                 for c in res["contacts"])
    return ifs, cts, res["summary"]["verdicts"]


@pytest.mark.parametrize("stem", ["1A8O", "1LCD", "2XHE"])
def test_pdb_and_mmcif_agree(stem):
    a = ifacekit.analyze(need(stem + ".pdb"), n_points=200)
    b = ifacekit.analyze(need(stem + ".cif"), n_points=200)
    ka, kb = _key(a), _key(b)
    assert ka[0] == kb[0]
    assert ka[2] == kb[2]
    assert ka[1] == kb[1]


def test_interface_area_matches_freesasa_pairwise():
    fs = pytest.importorskip("freesasa")
    res = ifacekit.analyze(need("7CFN.cif"), n_points=960)
    it = next(i for i in res["interfaces"] if {i["unit_a"], i["unit_b"]} == {"B", "G"})
    m = mdl.load(need("7CFN.cif"))
    ia = np.where(m.unit_idx == next(k for k, u in enumerate(m.units) if u["id"] == "B"))[0]
    ib = np.where(m.unit_idx == next(k for k, u in enumerate(m.units) if u["id"] == "G"))[0]

    def total(idx):
        s = fs.Structure()
        for n, k in enumerate(idx):
            r = m.residues[m.res_idx[k]]
            s.addAtom(f"{m.atom_name[k]:<4s}", r["resname"], r["resnum"], r["chain"], *map(float, m.xyz[k]))
            s.setRadius(n, float(m.radius[k]))
        return fs.calc(s, fs.Parameters({"algorithm": fs.ShrakeRupley, "probe-radius": 1.4,
                                         "n-points": 960})).totalArea()

    ref = 0.5 * (total(ia) + total(ib) - total(np.concatenate([ia, ib])))
    assert it["area"] == pytest.approx(ref, rel=0.02)


def test_json_roundtrip(tmp_path):
    res = ifacekit.analyze(need("1LCD.pdb"), n_points=200)
    s = json.dumps(res)
    assert json.loads(s)["schema_version"] == res["schema_version"]


def test_summary_counts_cover_unrecorded_pairs():
    res = ifacekit.analyze(need("1A8O.pdb"), n_points=100)
    s = res["summary"]
    assert s["pairs_evaluated"] > s["contact_records"]
    assert sum(s["verdicts"].values()) == s["pairs_evaluated"]


# ---------------------------------------------------------------- chain selection
def _intra_problems(res, unit):
    return sorted((tuple(sorted((c["atom1"], c["atom2"]))), c["verdict"], c["reason"])
                  for c in res["contacts"]
                  if c["unit1"] == unit and c["unit2"] == unit
                  and c["verdict"] in ("borderline", "unfavorable"))


def test_quality_mode_matches_whole_file_for_a_chain():
    whole = ifacekit.analyze(need("1LCD.pdb"), n_points=100)
    only_a = ifacekit.analyze(need("1LCD.pdb"), chain="A", n_points=100)
    assert only_a["input"]["mode"] == "quality"
    assert len(only_a["units"]) == 1
    assert _intra_problems(only_a, "A") == _intra_problems(whole, "A")


def test_pair_mode_matches_whole_file_interface():
    whole = ifacekit.analyze(need("7CFN.cif"), n_points=200)
    pr = ifacekit.analyze(need("7CFN.cif"), pair=("B", "G"), n_points=200)
    w = next(i for i in whole["interfaces"] if {i["unit_a"], i["unit_b"]} == {"B", "G"})
    assert len(pr["interfaces"]) == 1
    assert pr["interfaces"][0]["area"] == pytest.approx(w["area"], abs=0.5)
    key = lambda r: sorted((tuple(sorted((c["atom1"], c["atom2"]))), c["type"], c["verdict"])
                           for c in r["contacts"] if {c["unit1"], c["unit2"]} == {"B", "G"})
    assert key(pr) == key(whole)


def test_pair_mode_reports_only_interactions_between_the_chains():
    pr = ifacekit.analyze(need("7CFN.cif"), pair=("A", "B"), n_points=100)
    assert all(it["inter_unit"] for x in pr["issues"] for it in x["items"])
    assert {u["id"] for u in pr["units"]} == {"A", "B"}


def test_chain_ligands_join_the_chain():
    q = ifacekit.analyze(need("7CFN.cif"), chain="R", n_points=100)
    assert [u["id"] for u in q["units"]] == ["R"]
    assert len(q["units"][0]["ligands"]) == 4


def test_unknown_chain_lists_available():
    with pytest.raises(ValueError, match="chains in file: A, B, G, N, R"):
        ifacekit.analyze(need("7CFN.cif"), chain="Z")


def test_bad_selections_rejected():
    with pytest.raises(ValueError):
        ifacekit.analyze(need("7CFN.cif"), pair=("A", "A"))
    with pytest.raises(ValueError):
        ifacekit.analyze(need("7CFN.cif"), chain="A", pair=("A", "B"))


def test_cli_without_selection_lists_chains(capsys):
    from ifacekit.cli import main
    rc = main([need("1LCD.pdb")])
    out = capsys.readouterr().out
    assert rc == 2 and "Chains in this file" in out and "--pair" in out


def test_cli_pair_and_chains_write_json(tmp_path, capsys):
    from ifacekit.cli import main
    f = tmp_path / "pair.json"
    assert main([need("1LCD.pdb"), "--pair", "A", "B", "--points", "100", "-o", str(f)]) == 0
    assert json.loads(f.read_text())["input"]["mode"] == "pair"
    g = tmp_path / "chains.json"
    assert main([need("1LCD.pdb"), "--chains", "A", "B", "--points", "100", "-o", str(g)]) == 0
    runs = json.loads(g.read_text())["runs"]
    assert [r["input"]["chains"] for r in runs] == [["A"], ["B"]]
    out = capsys.readouterr().out
    assert "QUALITY of chain A" in out and "INTERACTION between chain A and chain B" in out


def test_pair_without_contact_says_so(tmp_path):
    atoms = [("N", "GLY", "A", 1, "N", 0.0, 0.0, 0.0, False),
             ("N", "GLY", "B", 1, "N", 30.0, 0.0, 0.0, False)]
    p = tmp_path / "far.pdb"
    p.write_text(pdb_text(atoms))
    from ifacekit.report import text_report
    res = ifacekit.analyze(str(p), pair=("A", "B"), n_points=100)
    assert res["interfaces"] == []
    assert "No interface" in text_report(res)


def test_chains_all_checks_every_chain_with_overview(tmp_path, capsys):
    from ifacekit.cli import main
    f = tmp_path / "all.json"
    assert main([need("1LCD.pdb"), "--chains", "all", "--points", "100", "-o", str(f)]) == 0
    runs = json.loads(f.read_text())["runs"]
    assert sorted(r["input"]["chains"][0] for r in runs) == ["A", "B", "C"]
    out = capsys.readouterr().out
    assert out.startswith("Overview") and out.count("QUALITY of chain") == 3


def test_single_chain_result_unchanged_by_lazy_sasa():
    # a lone chain has no interface, so the SASA step is skipped; results must not depend on it
    a = ifacekit.analyze(need("1LCD.pdb"), chain="A", n_points=100)
    assert a["interfaces"] == []
    assert a["summary"]["pairs_evaluated"] > 0


# ---------------------------------------------------------------- attachment
def _relabel(src, dst, waters_to="W"):
    import gemmi
    st = gemmi.read_structure(str(src))
    st.setup_entities()
    model = st[0]
    wat = gemmi.Chain(waters_to)
    for ch in model:
        keep = []
        for r in ch:
            if r.name == "HOH":
                wat.add_residue(r)
        # drop waters from the original chain
        idx = [i for i, r in enumerate(ch) if r.name == "HOH"]
        for i in reversed(idx):
            del ch[i]
    model.add_chain(wat)
    st.write_pdb(str(dst))
    return dst


def test_waters_in_separate_chain_follow_contacts(tmp_path):
    src = need("5hpt.cif")
    new = _relabel(src, tmp_path / "w.pdb")
    ref = ifacekit.analyze(src, chain="A")
    got = ifacekit.analyze(str(new), chain="A")
    lab = ifacekit.analyze(str(new), chain="A", attach="label")
    assert got["flags"]["n_water_residues"] == ref["flags"]["n_water_residues"] > 0
    assert lab["flags"]["n_water_residues"] == 0
    assert got["summary"]["verdicts"] == ref["summary"]["verdicts"]
    with pytest.raises(ValueError, match="only ligands or waters"):
        ifacekit.analyze(str(new), chain="W")
    rows = {r["chain"]: r for r in ifacekit.list_chains(str(new))}
    assert rows["W"]["has_polymer"] is False and rows["A"]["has_polymer"] is True


def test_covalent_glycan_group_stays_with_its_chain(tmp_path):
    import gemmi
    import numpy as np
    src = need("5hpt.cif")
    st = gemmi.read_structure(src)
    st.setup_entities()
    model = st[0]
    # pick an Asn ND2 on chain A and grow a 2-residue chain away from the protein
    asn = next(r for r in model["A"] if r.name == "ASN")
    nd2 = asn["ND2"][0].pos
    cg = asn["CG"][0].pos
    u = np.array([nd2.x - cg.x, nd2.y - cg.y, nd2.z - cg.z])
    u /= np.linalg.norm(u)
    # walk outward until nothing is within 4.5 A, then bond
    ns = gemmi.NeighborSearch(model, st.cell, 6).populate()
    base = np.array([nd2.x, nd2.y, nd2.z])
    ch = gemmi.Chain("G")
    for k, resname in enumerate(["NAG", "NAG"]):
        r = gemmi.Residue()
        r.name, r.seqid, r.het_flag = resname, gemmi.SeqId(k + 1, " "), "H"
        for j, nm in enumerate(["C1", "O5"]):
            a = gemmi.Atom()
            a.name, a.element = nm, gemmi.Element(nm[0])
            p = base + u * (1.45 + 1.5 * (2 * k + j))
            a.pos = gemmi.Position(*p)
            a.occ = 1.0
            r.add_atom(a)
        ch.add_residue(r)
    model.add_chain(ch)
    out = tmp_path / "g.pdb"
    st.write_pdb(str(out))
    from ifacekit.attach import plan_attachments
    own = plan_attachments(str(out))
    owners = {k: v for k, v in own.items() if ":NAG" in k}
    assert len(owners) == 2
    assert all(v and v[0] == "A" for v in owners.values()), owners


# ---------------------------------------------------------------- waters
def _water_case(tmp_path):
    atoms = [("OG", "SER", "A", 1, "O", 0.0, 0.0, 0.0, False)]
    n = 2

    def wat(x, y, z):
        nonlocal n
        atoms.append(("O", "HOH", "A", n, "O", x, y, z, True))
        n += 1
    wat(2.9, 0, 0)          # one partner
    wat(0, 2.2, 0)          # too close to OG
    wat(30, 0, 0)           # isolated
    for cx, k in ((100.0, 6), (200.0, 5)):
        atoms.append(("MG", "MG", "A", n, "MG", cx, 0.0, 0.0, True))
        n += 1
        for dx, dy, dz in [(2.1, 0, 0), (-2.1, 0, 0), (0, 2.1, 0), (0, -2.1, 0), (0, 0, 2.1), (0, 0, -2.1)][:k]:
            wat(cx + dx, dy, dz)
    p = tmp_path / "w.pdb"
    p.write_text(pdb_text(atoms))
    return str(p)


def test_water_categories_and_metal_sites(tmp_path):
    r = ifacekit.analyze(_water_case(tmp_path), n_points=200)
    w = r["waters"]
    cat = {x["water"]: x for x in w["listed"]}
    assert w["counts"]["unsupported"] == 1
    assert any(x["category"] == "unsupported" for x in w["listed"])
    assert "too_close" in cat["A:3:HOH"]["flags"]
    assert "A:2:HOH" not in cat                      # one partner at 2.9 A is just weak, not flagged
    assert w["counts"]["weak"] >= 1
    assert w["counts"]["metal_bound"] == 11
    sites = {x["metal"].split(":")[1]: x for x in w["metals"]}
    assert len(sites) == 2
    ok = next(x for x in w["metals"] if x["n_ligands"] == 6)
    bad = next(x for x in w["metals"] if x["n_ligands"] == 5)
    assert ok["verdict"] == "ok" and ok["octahedral_rms_deg"] < 1
    assert bad["verdict"] == "borderline" and "5 ligands" in bad["notes"][0]


def test_water_judged_with_neighbours_outside_the_scope(tmp_path):
    # the water sits between two chains; inspecting chain A alone must still count chain B's partner
    atoms = [("OG", "SER", "A", 1, "O", 0.0, 0.0, 0.0, False),
             ("OG", "SER", "B", 1, "O", 5.6, 0.0, 0.0, False),
             ("O", "HOH", "W", 1, "O", 2.8, 0.0, 0.0, True)]
    p = tmp_path / "s.pdb"
    p.write_text(pdb_text(atoms))
    r = ifacekit.analyze(str(p), chain="A", n_points=200)
    assert r["waters"]["n"] == 1 and r["waters"]["counts"].get("supported") == 1
    assert ifacekit.analyze(str(p), chain="A", attach="label")["waters"]["n"] == 0


def test_metal_site_consistency_between_identical_chains(tmp_path):
    from ifacekit.consistency import metal_consistency
    atoms, n = [], 1
    for ch, cx, k in (("A", 0.0, 6), ("B", 50.0, 6), ("C", 100.0, 5)):
        atoms.append(("OG", "SER", ch, 1, "O", cx + 20, 20.0, 0.0, False))
        atoms.append(("MG", "MG", ch, 2, "MG", cx, 0.0, 0.0, True))
        for dx, dy, dz in [(2.1, 0, 0), (-2.1, 0, 0), (0, 2.1, 0), (0, -2.1, 0), (0, 0, 2.1), (0, 0, -2.1)][:k]:
            atoms.append(("O", "HOH", ch, 10 + n, "O", cx + dx, dy, dz, True))
            n += 1
    p = tmp_path / "c.pdb"
    p.write_text(pdb_text(atoms))
    c = metal_consistency(str(p))
    g = c["groups"][0]
    assert g["chains"] == ["A", "B", "C"]
    s = g["sites"][0]
    assert s["consensus"]["n_water"] == 6 and s["consensus"]["n_chains"] == 2
    assert [d["chain"] for d in s["deviants"]] == ["C"] and s["deviants"][0]["n_water"] == 5


# ---------------------------------------------------------------- geometry and protonation
def _backbone(tmp_path, n_ca, name="g.pdb"):
    # N at origin, CA n_ca along x, C placed with an ideal N-CA-C angle of 111 degrees
    import math
    ca = (n_ca, 0.0, 0.0)
    ang = math.radians(111.2)
    c = (ca[0] + 1.525 * math.cos(math.pi - ang), 1.525 * math.sin(math.pi - ang), 0.0)
    atoms = [("N", "ALA", "A", 1, "N", 0.0, 0.0, 0.0, False), ("CA", "ALA", "A", 1, "C", *ca, False),
             ("C", "ALA", "A", 1, "C", *c, False)]
    p = tmp_path / name
    p.write_text(pdb_text(atoms))
    return str(p)


def test_backbone_geometry_flags_stretched_bond(tmp_path):
    good = ifacekit.analyze(_backbone(tmp_path, 1.459, "a.pdb"), chain="A", n_points=200)["geometry"]
    bad = ifacekit.analyze(_backbone(tmp_path, 1.80, "b.pdb"), chain="A", n_points=200)["geometry"]
    assert good["n_outliers"] == 0
    assert any(o["kind"] == "bond" and o["atoms"] == "N-CA" for o in bad["outliers"])


def test_acid_oxygen_pair_is_possible_protonation_not_unfavorable(tmp_path):
    r = pair(tmp_path, ("OD1", "ASP", "O"), ("O", "GLY", "O"), 2.4)
    c = only(r)
    assert c["type"] == "polar_mismatch" and c["verdict"] == "borderline" and c["reason"] == "possible_protonation"
    r2 = pair(tmp_path, ("O", "ALA", "O"), ("O", "GLY", "O"), 2.4)
    assert only(r2)["verdict"] == "unfavorable"


# ---------------------------------------------------------------- AF3 scores and HTML report
def test_ipsae_limits():
    from ifacekit import af3
    n = 10
    full = {"pae": np.zeros((n, n)), "chain": np.array(list("AAAAABBBBB")), "res": np.arange(1, n + 1)}
    assert af3.ipsae(full, "A", "B") == pytest.approx(1.0)
    assert af3.ipsae({**full, "pae": np.full((n, n), 30.0)}, "A", "B") == 0.0
    ip = af3.interface_pae(full, "A", "B", {1, 2}, {6, 7})
    assert ip["mean"] == 0.0 and ip["fraction_below_5"] == 1.0


def test_html_report_is_written(tmp_path):
    from ifacekit import htmlreport
    atoms = [("N", "ALA", "A", 1, "N", 0.0, 0.0, 0.0, False), ("CA", "ALA", "A", 1, "C", 1.45, 0.0, 0.0, False),
             ("C", "ALA", "A", 1, "C", 2.0, 1.4, 0.0, False), ("CB", "ALA", "A", 1, "C", 2.0, -0.8, 1.2, False),
             ("N", "ALA", "B", 1, "N", 3.5, -1.0, 1.0, False), ("CA", "ALA", "B", 1, "C", 4.9, -1.0, 1.0, False),
             ("C", "ALA", "B", 1, "C", 5.4, 0.4, 1.0, False), ("CB", "ALA", "B", 1, "C", 5.4, -1.8, 2.2, False)]
    p = tmp_path / "m.pdb"
    p.write_text(pdb_text(atoms))
    out = tmp_path / "r.html"
    out3 = tmp_path / "r_3d.html"
    htmlreport.build(str(p), "A", "B", str(out), molstar_path=str(out3))
    text = out.read_text()
    assert "<h2>1. AI confidence metrics</h2>" in text and "<h2>3. Quality of the prediction</h2>" in text
    assert "molstar" not in text.lower() and "data-a=" not in text      # the plain page is unchanged
    oe = tmp_path / "e.html"
    htmlreport.build(str(p), "A", "B", str(oe), kind="experimental")
    te = oe.read_text()
    assert "AI confidence" not in te and "B-factor along the sequence" in te and "1. Interface interactions" in te
    t3 = out3.read_text()
    assert "molstar" in t3 and "data-a=" in t3 and "mvpanel" in t3


def test_ligand_diagram_and_report(tmp_path):
    pytest.importorskip("rdkit")
    from ifacekit import ligplot, ligreport
    from pdbdata import no_dictionary
    f = no_dictionary()                                                  # 8JUX without the bond table of its ligand V00
    res = ifacekit.analyze(f)
    d = ligplot.build(f, res, "A:V00401")
    from rdkit.Chem import rdMolDescriptors
    assert rdMolDescriptors.CalcMolFormula(d["mol"]) == "C26H29FN2O2"      # bond orders rebuilt from the geometry give the formula of the dictionary
    assert len(d["nodes"]) == 13
    assert sum(1 for n in d["nodes"] for c in n["contacts"] if c["type"] == "hbond") == 2
    svg = ligplot.svg(d)
    assert svg.startswith("<svg") and "Tyr268" in svg and "Gln269" in svg
    ligreport.build(f, str(tmp_path / "l.html"), "A:V00401", molstar_path=str(tmp_path / "l3.html"))
    assert "Ligand binding" in (tmp_path / "l.html").read_text() and "mvpanel" in (tmp_path / "l3.html").read_text()


def test_ligand_dictionary_and_views(tmp_path):
    pytest.importorskip("rdkit")
    from ifacekit import ligplot, ligreport
    f = need("4QZS.cif")                                   # file carries the chemical component dictionary
    d = ligplot.build(f, ifacekit.analyze(f), "A:JQ1201")
    assert d["bond_source"] == "file dictionary"
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    assert rdMolDescriptors.CalcMolFormula(d["mol"]) == "C23H25ClN4O2S"      # the real JQ1 formula, hydrogens filled in from the correct bond orders
    assert any(n["rn"] == "HOH" for n in d["nodes"])                           # waters are drawn
    html_path = tmp_path / "j.html"
    ligreport.build(f, str(html_path))
    t = html_path.read_text()
    for token in ("dg-modern", "dg-classic", "dgSave('png')", "dgSave('pdf')", "dgSave('svg')", "feDropShadow"):
        assert token in t
    assert t.count('aria-label="Interaction diagram of') == 4            # modern and classic, simple and all contacts
    for token in ("id='dg-modern-all'", "id='dg-classic-all'", "dgDetail('all')", "id='mediated'", 'H-bonds via water'):
        assert token in t


def test_model_type_detection_and_labels(tmp_path):
    pytest.importorskip("rdkit")
    from ifacekit import ligplot, ligreport
    exp, pred = need("4QZS.cif"), prediction_like()               # 4QZS changed to look like an AI prediction
    assert ifacekit.analyze(exp)["flags"]["predicted_model"] is False
    assert ifacekit.analyze(pred)["flags"]["predicted_model"] is True
    for f, word, other in ((exp, "B-factor", "pLDDT"), (pred, "pLDDT", "B-factor")):
        ligreport.build(f, str(tmp_path / "x.html"))
        t = (tmp_path / "x.html").read_text()
        assert f"see distance and {word}" in t and f"see distance and {other}" not in t
    ligreport.build(pred, str(tmp_path / "y.html"), predicted=False)       # the user can override the detection
    assert "see distance and B-factor" in (tmp_path / "y.html").read_text()


def test_figures_written_by_the_program(tmp_path):
    pytest.importorskip("rdkit")
    from ifacekit import cli
    f = need("4QZS.cif")
    rc = cli.main([f, "--ligand", "A:JQ1201", "--html", str(tmp_path / "j_report.html"), "-q"])
    assert rc == 0
    names = sorted(x.name for x in tmp_path.iterdir())
    assert "j_map_modern.svg" in names and "j_map_classic.svg" in names and "j_report.html" in names
    try:
        import cairosvg  # noqa: F401
    except ImportError:
        return
    assert "j_map_modern.png" in names and "j_map_classic.pdf" in names


def test_aromatic_rings_drawn_as_alternating_bonds_in_both_styles():
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from ifacekit import ligplot
    f = need("8JUX.cif")                                       # naphthalene, fluorophenyl, piperidine and one C=O
    d = ligplot.build(f, ifacekit.analyze(f), "A:V00401")
    n_bonds = d["mol"].GetNumBonds()
    km = ligplot.kekule(d["mol"])
    n_double = sum(1 for b in km.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE)
    assert n_double == 9 and not any(b.GetBondTypeAsDouble() == 1.5 for b in km.GetBonds())
    assert ligplot.svg(d).count('class="bond"') == n_bonds + n_double          # each double bond is two lines, nothing else


def test_ramachandran_reference():
    import os
    from ifacekit import rama
    files = [need(n) for n in ("1A7G.cif", "2XHE.cif", "4QZS.cif", "5hpt.cif")]
    summ, rows = rama.analyze(files[0])                # a fixed file, not "the first one the file system lists"
    assert summ["n"] == len(rows) > 0
    assert sum(summ["counts"].values()) == summ["n"]
    # The category follows the class cut-offs. The density in the rows is rounded to 5 decimals while the category comes from the
    # unrounded value, so a residue whose density lies within that rounding of a cut-off may differ: allow exactly that.
    for f in files:
        for r in rama.analyze(f)[1]:
            if r["category"] != rama.classify_density(r["class"], r["density"]):
                near = min(abs(r["density"] - rama.FAV_T[r["class"]]), abs(r["density"] - rama.ALW_T[r["class"]]))
                assert near <= 1e-5, (f, r)
    # a canonical alpha-helix residue is favored, a left-handed point for a general residue is not an outlier-free guarantee
    assert rama.classify_density("general", rama.density("general", -63.0, -43.0)) == "favored"
    assert rama.classify_density("general", rama.density("general", 60.0, -120.0)) in ("allowed", "outlier")


def test_ramachandran_in_pair_report(tmp_path):
    import os
    import subprocess
    import sys
    src = need("5hpt.cif")
    out = tmp_path / "r.html"
    subprocess.run([sys.executable, "-m", "ifacekit", src, "--pair", "D", "E", "--experimental", "--html", str(out)],
                   check=True, capture_output=True)
    t = out.read_text()
    assert "Ramachandran (chains D, E)" in t and "<svg" in t


def test_quality_html_report(tmp_path):
    import os
    import subprocess
    import sys
    src = need("5hpt.cif")
    out = tmp_path / "q.html"
    subprocess.run([sys.executable, "-m", "ifacekit", src, "--chains", "D", "E", "--html", str(out), "-q"],
                   check=True, capture_output=True)
    t = out.read_text()
    assert "chain quality" in t and 'id="c-D"' in t and 'id="c-E"' in t
    assert "Ramachandran (chains D)" in t and "Ramachandran (chains E)" in t
    assert "Mol*" not in t and "cdn.jsdelivr" not in t          # plain report works offline


def test_rotamer_reference():
    import os
    import gemmi
    from ifacekit import rotamer
    path = need("1A7G.cif")
    summ, rows = rotamer.analyze(path)
    assert summ["n"] == len(rows) > 0
    assert sum(summ["counts"].values()) == summ["n"]
    for r in rows:
        assert r["category"] == rotamer.classify_density(r["table"], r["density"])
        assert all(0 <= x < 360 for x in r["chi"])
    # number of chi angles follows the definitions; a gauche-minus Leu rotamer (mt) is favored, a strained one is not
    assert rotamer.classify_density("leu", rotamer.density("leu", [295.0, 175.0])) == "favored"
    assert rotamer.classify_density("leu", rotamer.density("leu", [20.0, 20.0])) == "outlier"
    st = gemmi.read_structure(path)
    for ch in st[0]:
        for res in ch:
            if res.name in rotamer.CHI:
                c = rotamer.chi_values(res)
                assert c is None or len(c) == len(rotamer.CHI[res.name])


def test_rotamer_in_quality_report(tmp_path):
    import os
    import subprocess
    import sys
    src = need("5hpt.cif")
    out = tmp_path / "q.html"
    subprocess.run([sys.executable, "-m", "ifacekit", src, "--chains", "D", "--html", str(out), "-q"], check=True, capture_output=True)
    assert "Rotamers (chains D)" in out.read_text()


def test_flip_check_and_flipped_rotamers():
    import os
    from ifacekit import flips, rotamer
    path = need("1A7G.cif")
    summ, rows = flips.analyze(path)
    assert summ["n"] == len(rows) and sum(summ["counts"].values()) == len(rows)
    for r in rows:
        assert abs(r["delta"] - (r["score_flipped"] - r["score_modelled"])) < 0.011
        u, f = flips.MARGIN[r["type"]]
        assert r["verdict"] == ("flip" if r["delta"] >= f else "uncertain" if r["delta"] >= u else "keep")
    _, rrows = rotamer.analyze(path)
    flipped = [r for r in rrows if r["table"] in ("asn", "gln", "his")]
    assert all("category_flipped" in r and r["flip"] in ("keep", "uncertain", "flip") for r in flipped)
    assert all("category_flipped" not in r for r in rrows if r["table"] not in ("asn", "gln", "his"))


def test_flip_changes_the_right_chi():
    # flipping an Asn amide turns chi2 by 180 degrees, Gln chi3, His chi2 (ring)
    from ifacekit import rotamer
    assert rotamer.FLIP_CHI == {"ASN": 1, "HIS": 1, "GLN": 2}


def test_rotamer_plots_in_quality_report(tmp_path):
    import os
    import subprocess
    import sys
    src = need("5hpt.cif")
    out, rot = tmp_path / "q.html", tmp_path / "rot.html"
    subprocess.run([sys.executable, "-m", "ifacekit", src, "--chains", "D", "--html", str(out), "--rotamer-html", str(rot), "-q"],
                   check=True, capture_output=True)
    for f in (out, rot):
        t = f.read_text()
        assert "Leu rotamer plot" in t and "Val rotamer plot" in t and "Asn, Gln and His flips" in t


def test_plots_pdf_and_link(tmp_path):
    import os
    import subprocess
    import sys
    import pytest
    pytest.importorskip("cairosvg")
    pypdf = pytest.importorskip("pypdf")
    src = need("5hpt.cif")
    pdf, html = tmp_path / "plots.pdf", tmp_path / "q.html"
    subprocess.run([sys.executable, "-m", "ifacekit", src, "--chains", "D", "E", "--html", str(html), "--plots-pdf", str(pdf), "-q"],
                   check=True, capture_output=True)
    r = pypdf.PdfReader(str(pdf))
    assert len(r.pages) >= 8                     # summary, then Ramachandran, rotamer and outlier pages for each chain
    text = " ".join(p.extract_text() or "" for p in r.pages)
    assert "Chain D: Ramachandran" in text and "Chain E: rotamers" in text
    t = html.read_text()
    assert "href='plots.pdf'" in t and "bold 15px" in t       # download link, and the chain bars of the sequence plot
    i = t.index("Ramachandran (chains D)")
    assert 0 < t.index("href='plots.pdf'") - i < 400           # the link sits directly under the Ramachandran header
    assert "a.pdfdl:visited" in t and "color:var(--acc)" in t   # explicit colour so that it is readable in dark mode


def test_no_pdf_link_without_pdf(tmp_path):
    import os
    import subprocess
    import sys
    src = need("5hpt.cif")
    html = tmp_path / "q.html"
    subprocess.run([sys.executable, "-m", "ifacekit", src, "--chains", "D", "--html", str(html), "-q"], check=True, capture_output=True)
    assert "Download all" not in html.read_text()


def test_portal_keeps_structure_name():
    import importlib.util
    import os
    import pytest
    pytest.importorskip("flask")
    spec = importlib.util.spec_from_file_location("ifk_portal", os.path.join("server", "app.py"))
    app = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(app)
    f = app.safe_structure_name
    assert f("8afb.cif", "cif") == "8afb.cif"
    assert f("my model (v2).cif", "cif") == "my_model__v2.cif"
    assert f("../../evil name.cif", "cif") == "evil_name.cif"          # no path parts
    assert f(".hidden.pdb", "pdb") == "hidden.pdb"                        # no leading dot
    assert f("<script>.cif", "cif") == "script.cif" or f("<script>.cif", "cif").startswith("_script")
    assert f("", "pdb") == "structure.pdb" and f(".cif", "cif") == "structure.cif"
    assert f("model.pdb", "cif") == "model.cif"                           # extension follows the checked content
    assert len(f("x" * 300 + ".cif", "cif")) <= 64


def test_simple_ligand_map_is_a_subset():
    pytest.importorskip("rdkit")
    from ifacekit import ligplot
    f = need("4QZS.cif")
    r = ifacekit.analyze(f)
    full = ligplot.build(f, r, "A:JQ1201", detail="full")
    simple = ligplot.build(f, r, "A:JQ1201", detail="simple")
    key = lambda n: (n["ch"], n["num"], n["rn"])
    assert {key(n) for n in simple["nodes"]} <= {key(n) for n in full["nodes"]}
    assert len(simple["nodes"]) < len(full["nodes"])
    for n in simple["nodes"]:
        if n["rn"] == "HOH":                                          # a water is drawn only when it H-bonds the ligand
            assert any(c["type"] == "hbond" for c in n["contacts"])
        else:                                                         # a residue: polar contact, or two contacts with the closest within 3.8 A
            assert n["polar"] or (len(n["contacts"]) >= 2 and min(c["d"] for c in n["contacts"]) <= 3.8)


def test_hbonds_direct_and_via_water():
    from ifacekit import hbonds
    def c(a1, a2, u1, u2, typ="hbond", d=2.9):
        return {"atom1": a1, "atom2": a2, "unit1": u1, "unit2": u2, "type": typ, "distance": d, "inter_unit": u1 != u2}
    contacts = [c("A:1:X:N", "B:2:Y:O", "A", "B"),                                   # direct
                c("A:3:X:O", "A:100:HOH:O", "A", "water"), c("B:4:Y:N", "A:100:HOH:O", "B", "water"),   # water 100 bridges A and B
                c("A:5:X:O", "A:101:HOH:O", "A", "water"),                                              # water 101 only touches A
                c("A:6:X:C", "B:7:Y:C", "A", "B", typ="vdw")]
    h = hbonds.split(contacts, lambda u: u == "A", lambda u: u == "B")
    assert (h["direct"], h["water_mediated"], h["total_connections"]) == (1, 1, 2)
    assert h["mediated_waters"][0]["water"] == "A:100:HOH:O" and h["side_a_to_water"] == 2
    assert hbonds.split([contacts[0]], lambda u: u == "A", lambda u: u == "B")["water_mediated"] == 0   # no waters: zero


def test_hbond_fields_in_results():
    f = need("5hpt.cif")
    r = ifacekit.analyze(f, pair=("D", "E"))
    it = r["interfaces"][0]
    assert it["hbonds"]["direct_polar_contacts"] == it["contacts"].get("hbond", 0)      # the distance-only count is the polar contacts of type hbond
    assert it["hbonds"]["direct_polar_contacts"] == 4 and it["hbonds"]["direct"] == 2   # two backbone bonds (2.97 and 3.00 A); the contacts of 3.37 and 3.48 A are not hydrogen bonds
    assert it["hbonds"]["total_connections"] == it["hbonds"]["direct"] + it["hbonds"]["water_mediated"]


def _portal_module():
    import importlib.util
    import os
    pytest.importorskip("flask")
    spec = importlib.util.spec_from_file_location("ifk_portal2", os.path.join("server", "app.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_reports_carry_licence_footer(tmp_path):
    import os
    import subprocess
    import sys
    src = need("5hpt.cif")
    pair, qual = tmp_path / "p.html", tmp_path / "q.html"
    subprocess.run([sys.executable, "-m", "ifacekit", src, "--pair", "D", "E", "--experimental", "--html", str(pair), "-q"], check=True, capture_output=True)
    subprocess.run([sys.executable, "-m", "ifacekit", src, "--chains", "D", "--html", str(qual), "-q"], check=True, capture_output=True)
    for f in (pair, qual):
        t = f.read_text()
        assert "MIT licence" in t and "Kuen-Phon Wu" in t and ifacekit.__version__ in t


def test_ligand_report_has_footer_and_water_split(tmp_path):
    pytest.importorskip("rdkit")
    from ifacekit import ligreport
    f = need("4QZS.cif")
    out = tmp_path / "l.html"
    d = ligreport.build(f, str(out), "A:JQ1201")
    t = out.read_text()
    assert "MIT licence" in t and "H-bonds to waters" in t and "Hydrogen bonds via water" in t
    assert d["hbonds"]["direct"] == 1 and d["hbonds"]["side_a_to_water"] == 2      # Asn140, and two waters (2.83 and 2.92 A)


def test_molstar_local_address(tmp_path):
    import os
    import subprocess
    import sys
    src = need("5hpt.cif")
    out = tmp_path / "p.html"
    subprocess.run([sys.executable, "-m", "ifacekit", src, "--pair", "D", "E", "--experimental", "--html", str(out), "--molstar",
                    "--molstar-url", "/static/molstar/", "-q"], check=True, capture_output=True)
    t3 = (tmp_path / "p_3d.html").read_text()
    assert "/static/molstar/" in t3 and "cdn.jsdelivr.net" not in t3
    subprocess.run([sys.executable, "-m", "ifacekit", src, "--pair", "D", "E", "--experimental", "--html", str(out), "--molstar", "-q"], check=True, capture_output=True)
    assert "cdn.jsdelivr.net" in (tmp_path / "p_3d.html").read_text()          # default: the CDN


def test_portal_cleanup_health_home_and_limits(tmp_path):
    import os
    import time
    p = _portal_module()
    p.DATA = str(tmp_path)
    old, new = "a" * 32, "b" * 32
    for jid in (old, new):
        os.makedirs(tmp_path / jid / "out")
    t = time.time() - (p.KEEP_DAYS + 1) * 86400
    os.utime(tmp_path / old, (t, t))
    p.purge_old()
    assert not (tmp_path / old).exists() and (tmp_path / new).exists()          # only the old job folder goes
    (tmp_path / "keepme.txt").write_text("x")
    p.purge_old()
    assert (tmp_path / "keepme.txt").exists()                                   # nothing else is touched
    c = p.app.test_client()
    page = c.get("/").get_data(as_text=True)
    assert "MIT licence" in page and ifacekit.__version__ in page and "__VERSION__" not in page
    h = c.get("/healthz").get_json()
    assert h["status"] == "ok" and h["version"] == ifacekit.__version__ and h["jobs_stored"] == 1 and h["max_jobs"] == p.MAX_JOBS
    p.MEM_MB = 1024
    cmd = p._limited(["python", "x"])
    assert cmd[0] == "/bin/sh" and "ulimit -v" in cmd[2] and cmd[3] == str(1024 * 1024) and cmd[-2:] == ["python", "x"]
    p.MEM_MB = 0
    assert p._limited(["python", "x"]) == ["python", "x"]
    assert p.out_of_memory(type("R", (), {"returncode": 1, "stderr": "OpenBLAS error: Memory allocation still failed"})())
    assert not p.out_of_memory(type("R", (), {"returncode": 0, "stderr": "Memory allocation"})())


def test_deploy_files_present():
    import os
    for f in ("server/ifacekit-server.sh", "server/deploy/ifacekit.service", "server/deploy/nginx-ifacekit.conf.example",
              "server/deploy/README_DEPLOY.txt", "server/static/molstar/molstar.js", "server/static/molstar/molstar.css",
              "server/static/molstar/LICENSE", "LICENSE", "CITATION.cff"):
        assert os.path.exists(f), f
    assert os.access("server/ifacekit-server.sh", os.X_OK)


def test_detailed_ligand_diagram(tmp_path):
    pytest.importorskip("rdkit")
    from ifacekit import ligdetail, ligplot
    f = need("4QZS.cif")                                        # JQ1: Asn140 side chain and two waters hydrogen-bond the ligand
    d = ligplot.build(f, ifacekit.analyze(f), None, detail="simple")
    s = ligdetail.svg_detail(d)
    assert s.startswith("<svg") and s.count("<svg") == 1
    assert "Asn140" in s and 'class="dd-b"' in s                # the residue is a structure with bonds, not a box
    assert 'class="dd-hb"' in s and "2.83" in s and "3.16" in s    # hydrogen bonds with their distances
    assert 'class="b"' not in s and 'class="al"' not in s       # class names do not clash with the other diagrams of the page
    # the Asn140 contacting atom is ND2 and a ring marks it
    assert s.count('stroke="#2f855a" stroke-width="2" opacity=".55"') >= 1


def test_detailed_option_in_report_and_figures(tmp_path):
    pytest.importorskip("rdkit")
    import subprocess
    import sys
    f = need("4QZS.cif")
    plain, det = tmp_path / "plain.html", tmp_path / "det.html"
    subprocess.run([sys.executable, "-m", "ifacekit", f, "--ligand", "--html", str(plain)], check=True, capture_output=True)
    subprocess.run([sys.executable, "-m", "ifacekit", f, "--ligand", "--detailed", "--html", str(det)], check=True, capture_output=True)
    assert "id='dg-detailed'" not in plain.read_text() and "Detailed</button>" not in plain.read_text()
    t = det.read_text()
    assert "id='dg-detailed'" in t and "id='dg-detailed-all'" in t and "id='v-detailed'" in t
    names = {p.name for p in tmp_path.iterdir()}
    assert "det_map_detailed.svg" in names and "det_map_all_detailed.svg" in names
    assert not any("detailed" in n for n in names if n.startswith("plain"))


def test_detailed_all_residues_option(tmp_path):
    pytest.importorskip("rdkit")
    from ifacekit import ligdetail, ligplot
    f = need("4QZS.cif")
    d = ligplot.build(f, ifacekit.analyze(f), None, detail="simple")
    few, many = ligdetail.svg_detail(d), ligdetail.svg_detail(d, all_residues=True)
    # Leu92 only has van der Waals and hydrophobic contacts: a box by default, a structure with the option
    assert 'rx="9"' in few and few.count("<rect") > many.count("<rect")
    assert many.count('class="dd-b"') > few.count('class="dd-b"')
    import subprocess
    import sys
    out = tmp_path / "a.html"
    subprocess.run([sys.executable, "-m", "ifacekit", f, "--ligand", "--detailed", "--detailed-all", "--html", str(out)],
                   check=True, capture_output=True)
    assert (tmp_path / "a_map_detailed.svg").read_text().count('class="dd-b"') == many.count('class="dd-b"')


def test_pair_report_lists_only_contacts_between_the_chains(tmp_path):
    """Regression: the polar-contact table and the interface clash/issue counts once included every protein-water contact of the model."""
    import re
    from ifacekit import htmlreport
    f = need("5hpt.cif")
    r = ifacekit.analyze(f, pair=("D", "E"))
    water_rows = [c for c in r["contacts"] if "water" in (c["unit1"], c["unit2"]) and c["type"] in ("hbond", "salt_bridge", "metal")]
    assert water_rows                                           # the example does contain protein-water contacts, so the test is meaningful
    out = tmp_path / "p.html"
    htmlreport.build(f, "D", "E", str(out), kind="experimental")
    h = out.read_text()
    seg = h[h.index("id='polar'"):]
    seg = seg[:seg.index("<h3", 20)]
    between = [c for c in r["contacts"] if {c["unit1"], c["unit2"]} == {"D", "E"} and c["type"] in ("hbond", "salt_bridge", "metal")]
    assert "HOH" not in seg and seg.count("<tr") - 1 == len(between)
    assert f"between D and E ({len(between)})" in h


def test_split_problems_units_ignores_waters():
    from ifacekit.htmlreport import split_problems
    n = iter(range(1, 100))
    def c(u1, u2, reason="steric_overlap"):
        k = next(n)                                              # a different atom pair every time (a pair is counted once)
        return {"atom1": f"{u1}:{k}:X:C", "atom2": f"{u2}:{k}:Y:O", "unit1": u1, "unit2": u2, "verdict": "unfavorable", "reason": reason,
                "inter_unit": u1 != u2, "overlap": 0.7, "distance": 2.5}
    cs = [c("A", "B"), c("A", "water"), c("B", "water"), c("A", "A"), c("A", "B", "short_hbond")]
    both = split_problems(cs, only_between=True)               # between any two units: a water counts as a unit (used by the ligand report)
    strict = split_problems(cs, units=("A", "B"))              # only between the two chains
    assert (len(both[0]), len(both[1])) == (3, 1) and (len(strict[0]), len(strict[1])) == (1, 1)


def test_bridge_flags_from_geometry():
    import math
    from ifacekit import hbonds
    W = (0.0, 0.0, 0.0)

    def flags(a_xyz, a_d, b_xyz, b_d):
        pts = {"A:1:X:O": a_xyz, "B:2:Y:O": b_xyz}
        return hbonds.bridge_flags([("A:1:X:O", a_d)], [("B:2:Y:O", b_d)], W, pts.get)
    assert flags((2.8, 0, 0), 2.8, (0, 2.8, 0), 2.8) == []                                   # 90 degrees, partners 4 A apart: a clear bridge
    t = flags((2.8, 0, 0), 2.8, (2.8 * math.cos(math.radians(60)), 2.8 * math.sin(math.radians(60)), 0), 2.8)
    assert t == ["triangle", "narrow_angle"]                                                 # 60 degrees: partners 2.8 A apart
    assert flags((2.0, 0, 0), 2.0, (-2.8, 0, 0), 2.8) == ["too_short"]                       # 2.0 A: too close for a water (the partners are 4.8 A apart)
    assert flags((3.4, 0, 0), 3.4, (0, 2.8, 0), 2.8) == ["long_bond"]                        # 3.4 A: weak
    assert hbonds.bridge_flags([("A:1:X:O", 2.0)], [("B:2:Y:O", 2.8)]) == ["too_short"]      # without coordinates only the distance tests run


def test_split_counts_clear_and_doubtful_bridges():
    from ifacekit import hbonds
    def c(a1, a2, u1, u2, d):
        return {"atom1": a1, "atom2": a2, "unit1": u1, "unit2": u2, "type": "hbond", "distance": d, "inter_unit": u1 != u2}
    xyz = {"A:100:HOH:O": (0, 0, 0), "A:1:X:O": (2.8, 0, 0), "B:2:Y:O": (0, 2.8, 0),                       # a clear bridge
           "A:101:HOH:O": (10, 0, 0), "A:3:X:O": (12.0, 0, 0), "B:4:Y:O": (10, 2.8, 0)}                    # a water 2.0 A from its partner
    cs = [c("A:1:X:O", "B:9:Y:N", "A", "B", 2.9),
          c("A:1:X:O", "A:100:HOH:O", "A", "water", 2.8), c("B:2:Y:O", "A:100:HOH:O", "B", "water", 2.8),
          c("A:3:X:O", "A:101:HOH:O", "A", "water", 2.0), c("B:4:Y:O", "A:101:HOH:O", "B", "water", 2.8)]
    h = hbonds.split(cs, lambda u: u == "A", lambda u: u == "B", xyz=xyz.get)
    assert (h["direct"], h["water_mediated"], h["water_mediated_clear"], h["water_mediated_doubtful"]) == (1, 2, 1, 1)
    assert (h["total_connections"], h["total_connections_clear"]) == (3, 2)
    q = {x["water"]: (x["quality"], x["flags"]) for x in h["mediated_waters"]}
    assert q["A:100:HOH:O"] == ("clear", []) and q["A:101:HOH:O"][0] == "doubtful" and "too_short" in q["A:101:HOH:O"][1]


def test_reports_show_the_check_of_the_bridges(tmp_path):
    pytest.importorskip("rdkit")
    from ifacekit import ligreport
    f = need("4QZS.cif")
    out = tmp_path / "l.html"
    d = ligreport.build(f, str(out), "A:JQ1201")
    h = d["hbonds"]
    assert h["water_mediated_clear"] + h["water_mediated_doubtful"] == h["water_mediated"] > 0
    assert all(x["quality"] in ("clear", "doubtful") and (x["quality"] == "clear") == (not x["flags"]) for x in h["mediated_waters"])
    t = out.read_text()
    assert ">check<" in t and "clear" in t and "not calibrated" in t


def test_rotamer_scores_every_alternate_conformer():
    import gemmi
    from ifacekit import rotamer
    f = need("4QZS.cif")                                        # Asn52 has two conformers (A and B)
    summ, rows = rotamer.analyze(f)
    alt = [r for r in rows if r["altloc"]]
    assert alt and {r["altloc"] for r in alt} <= {"A", "B", "C"}
    ids = {}
    for r in alt:
        ids.setdefault(r["residue"], set()).add(r["altloc"])
    assert any(len(v) >= 2 for v in ids.values())               # both conformers of a residue are rows of their own
    st = gemmi.read_structure(f)
    for ch in st[0]:
        for res in ch:
            for r in alt:
                if r["residue"] == f"{ch.name}:{res.seqid.num}{res.seqid.icode.strip()}:{res.name}":
                    c = rotamer.chi_values(res, alt=r["altloc"])
                    assert [round(x, 1) for x in c] == r["chi"]  # the chi angles of that conformer
    assert summ["n"] == len(rows) and rotamer.label(alt[0]).endswith("alt " + alt[0]["altloc"])


def test_alternate_conformer_contacts_are_added_and_can_be_switched_off():
    from ifacekit import model as mdl
    f = need("4QZS.cif")
    on = ifacekit.analyze(f, n_points=120)
    off = ifacekit.analyze(f, n_points=120, alt_conformers=False)
    assert on["flags"]["alt_conformers_scanned"] and off["flags"]["alt_conformers_scanned"] == []
    extra = [c for c in on["contacts"] if c.get("altloc")]
    assert extra and not any(c.get("altloc") for c in off["contacts"])
    assert len(on["contacts"]) == len(off["contacts"]) + len(extra) and on["summary"]["alt_conformer_contacts"] == len(extra)
    assert on["summary"]["pairs_evaluated"] > off["summary"]["pairs_evaluated"]
    m = mdl.load(f, conformer="B")
    assert m.alt_extra.any() and m.n_atoms >= mdl.load(f).n_atoms - 5
    g = need("5hpt.cif")                                        # no alternate conformers: nothing changes
    a, b = ifacekit.analyze(g, n_points=120), ifacekit.analyze(g, n_points=120, alt_conformers=False)
    assert a["flags"]["alt_conformers_scanned"] == [] and len(a["contacts"]) == len(b["contacts"])


def test_footnotes_state_the_comparison_with_molprobity(tmp_path):
    import os
    import subprocess
    import sys
    from ifacekit import about
    f = need("5hpt.cif")
    out = tmp_path / "p.html"
    subprocess.run([sys.executable, "-m", "ifacekit", f, "--pair", "D", "E", "--html", str(out), "-q"], check=True, capture_output=True)
    t = out.read_text()
    assert "293 high-resolution X-ray structures" in t and "five structures" not in t and "not yet validated on many" not in t
    import importlib.util
    pytest.importorskip("flask")
    spec = importlib.util.spec_from_file_location("ifk_portal3", os.path.join("server", "app.py"))
    app = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(app)
    page = app.app.test_client().get("/").get_data(as_text=True)
    assert about.COMPARISON_SHORT in page and "five structures" not in page and "__COMPARISON_SHORT__" not in page


def test_calibrated_rotamer_cutoffs_file():
    import json
    import os
    from ifacekit import rotamer
    enc = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools", "rotamer_reference", "enclosure_calibrated.json")))
    assert set(enc) == set(rotamer.TABLES)
    assert all(0.3 < v[0] < v[1] < 1.0 for v in enc.values())
    assert enc["leu"][1] > 0.98 and enc["ile"][1] > 0.98          # the types that were over-called before the calibration


def test_factsheet_page_content_and_command(tmp_path):
    import os
    import re
    import subprocess
    import sys
    from ifacekit import about
    path = about.factsheet_path()
    assert os.path.exists(path)
    h = open(path, encoding="utf-8").read()
    for sec in ("summary", "clashes", "outliers", "calibration", "limits", "cases"):
        assert f'id="{sec}"' in h
    ids = re.findall(r"rcsb\.org/structure/([0-9A-Z]{4})", h)
    assert len(ids) == 300 and len(set(ids)) == 300                      # the list of test cases
    assert "The 900 calibration structures" in h and "Download the table as CSV" in h
    assert "MolProbity failed" in h                                      # the cases where MolProbity itself failed are listed
    # the numbers of the footnote and of the fact sheet agree
    cards = dict(re.findall(r'<div class="k">([^<]*)</div><div class="v">([^<]*)</div>', h))
    assert cards["heavy-atom clashes flagged"] in about.COMPARISON and cards["between two different residues"] in about.COMPARISON
    assert cards["Ramachandran outliers found"] in about.COMPARISON and cards["rotamer outliers found"] in about.COMPARISON
    assert cards["clashes found, hydrogens added"] in about.COMPARISON and cards["clashes confirmed, hydrogens added"] in about.COMPARISON
    assert cards["clashes found, hydrogens added"] in about.COMPARISON_SHORT and 'id="hclashes"' in h and 'id="chimerax"' in h
    assert "Why the two lists differ" in h and "only MolProbity reports" in h and "are also found by" in h
    out = tmp_path / "fs.html"
    subprocess.run([sys.executable, "-m", "ifacekit", "--factsheet", str(out)], check=True, capture_output=True)
    assert out.read_bytes() == open(path, "rb").read()
    r = subprocess.run([sys.executable, "-m", "ifacekit"], capture_output=True, text=True)
    assert r.returncode == 2 and "--factsheet" in r.stderr


def test_reports_link_to_the_factsheet_when_its_address_is_known(tmp_path):
    import subprocess
    import sys
    f = need("5hpt.cif")
    a, b = tmp_path / "a.html", tmp_path / "b.html"
    subprocess.run([sys.executable, "-m", "ifacekit", f, "--pair", "D", "E", "--html", str(a), "-q"], check=True, capture_output=True)
    subprocess.run([sys.executable, "-m", "ifacekit", f, "--pair", "D", "E", "--html", str(b), "--factsheet-url", "/factsheet", "-q"], check=True, capture_output=True)
    assert "ifacekit --factsheet FILE.html" in a.read_text() and "href='/factsheet'" not in a.read_text()
    assert "href='/factsheet'" in b.read_text() and "Fact sheet: the comparison and the list of test cases" in b.read_text()


def test_portal_serves_the_factsheet_and_links_it():
    p = _portal_module()
    c = p.app.test_client()
    r = c.get("/factsheet")
    assert r.status_code == 200 and b"ifacekit fact sheet" in r.data and b"The test cases" in r.data
    assert "href='/factsheet'" in c.get("/").get_data(as_text=True)


# ---------------------------------------------------------------- hydrogens (clash check with hydrogens added)
def _loaded_with_hydrogens(name="1A8O.cif"):
    from ifacekit import hydrogens as hy
    from ifacekit import model
    m = model.load(need(name))
    H, env, groups = hy.add_hydrogens(m)
    return m, H, env, groups, hy


def test_hydrogen_geometry_and_scope():
    m, H, env, groups, hy = _loaded_with_hydrogens("5hpt.cif")
    assert H.n > 1000 and not m.is_water[H.parent[H.active]].any()
    kinds = {r["kind"] for r in (m.residues[k] for k in set(H.res.tolist()))}
    assert kinds == {"polymer"}                                           # no hydrogens on ligands, waters or anything else
    d = np.linalg.norm(H.xyz - m.xyz[H.parent], axis=1)
    ok = np.zeros(len(d), bool)
    for length in (0.97, 0.86, 0.84, 1.20):                               # C-H, N-H, O-H, S-H (a flipped amide or ring carries its H on the other atom)
        ok |= np.isclose(d, length, atol=0.01)
    assert ok.all()
    assert np.allclose(d[m.element[H.parent] == "S"], 1.20, atol=0.01)
    # methyl groups: three hydrogens at the tetrahedral angle, staggered against the neighbour
    by_parent = {}
    for h in range(H.n):
        by_parent.setdefault(int(H.parent[h]), []).append(h)
    meth = [p for p, hs in by_parent.items() if len(hs) == 3 and m.element[p] == "C"]
    assert meth
    for p in meth[:20]:
        hs = by_parent[p]
        for a, b in ((0, 1), (0, 2), (1, 2)):
            u, v = H.xyz[hs[a]] - m.xyz[p], H.xyz[hs[b]] - m.xyz[p]
            ang = np.degrees(np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v)))
            assert 105 < ang < 114
    # the heavy atoms are never moved
    from ifacekit import model
    assert np.array_equal(model.load(need("5hpt.cif")).xyz, m.xyz)


def test_hydrogen_clash_rules():
    m, H, env, groups, hy = _loaded_with_hydrogens("1A8O.cif")
    base = len(hy.clashes(m, H, env))
    h = int(np.where((~H.polar) & H.active)[0][100])                      # an aliphatic hydrogen far from the residue we borrow an atom from
    far = [k for k in range(m.n_atoms) if m.residues[m.res_idx[k]]["resnum"] > m.residues[H.res[h]]["resnum"] + 30 and m.element[k] == "C"][0]
    m.xyz[far] = H.xyz[h] + np.array([1.9, 0.0, 0.0])                      # a carbon 1.9 A from the hydrogen
    env = hy.Env(m, H)
    want = env.rad[env.n + h] + env.rad[far] - 1.9                          # radii of the Probe program: overlap of 0.85 to 1.0 A
    found = hy.clashes(m, H, env)
    assert want > 0.8 and any(c["i"] == env.n + h and c["j"] == far and abs(c["overlap"] - want) < 0.01 for c in found)
    assert len(found) > base
    # a polar hydrogen close to an acceptor at a good angle is a hydrogen bond, not a clash
    p = int(np.where(H.polar & H.active & (m.element[H.parent] == "N"))[0][0])
    acc = [k for k in range(m.n_atoms) if m.element[k] == "O" and m.residues[m.res_idx[k]]["resnum"] > m.residues[H.res[p]]["resnum"] + 30][0]
    u = H.xyz[p] - m.xyz[H.parent[p]]
    m.xyz[acc] = H.xyz[p] + u / np.linalg.norm(u) * 1.85                   # in line with N-H, overlap 1.05 + 1.40 - 1.85 = 0.6 A
    env = hy.Env(m, H)
    assert not any(c["i"] == env.n + p and c["j"] == acc for c in hy.clashes(m, H, env))


def test_hydrogens_change_only_the_clash_check():
    from ifacekit.core import analyze
    f = need("5hpt.cif")
    a = analyze(f, chain="D")
    b = analyze(f, chain="D", hydrogens=True)
    assert a["hydrogens"] is None and "all_atom_clashes" not in a["summary"]
    h = b["hydrogens"]
    assert h["n_hydrogens"] > 1000 and h["n_clashes"] == len(h["clashes"]) == b["summary"]["all_atom_clashes"]
    assert h["clashes"] == sorted(h["clashes"], key=lambda c: -c["overlap"]) and all(c["overlap"] >= 0.4 for c in h["clashes"])
    for k in ("contacts", "issues", "interfaces", "ramachandran", "rotamer", "flips", "geometry", "waters"):
        assert a[k] == b[k], k                                                # every other result comes from the file as given
    assert analyze(f, chain="D", hydrogens=True)["hydrogens"]["clashes"] == h["clashes"]     # deterministic
    pair = analyze(f, pair=("D", "E"), hydrogens=True)["hydrogens"]
    assert all(c["inter_unit"] for c in pair["clashes"])                   # a pair analysis keeps the clashes between units


def test_hydrogens_option_in_command_line_and_reports(tmp_path):
    import json
    import subprocess
    import sys
    f = need("5hpt.cif")
    out, js, q = tmp_path / "p.html", tmp_path / "r.json", tmp_path / "q.html"
    r = subprocess.run([sys.executable, "-m", "ifacekit", f, "--pair", "D", "E", "--hydrogens", "--html", str(out), "-o", str(js)], check=True, capture_output=True, text=True)
    assert "All-atom clashes with hydrogens added" in r.stdout and json.load(open(js))["hydrogens"]["n_hydrogens"] > 1000
    t = out.read_text()
    assert "All-atom clashes with hydrogens added" in t and "Hydrogens of ligands, waters" in t
    subprocess.run([sys.executable, "-m", "ifacekit", f, "--chains", "D", "--hydrogens", "--html", str(q), "-q"], check=True, capture_output=True)
    assert "All-atom clashes with hydrogens added" in q.read_text()
    r = subprocess.run([sys.executable, "-m", "ifacekit", f, "--pair", "D", "E", "--html", str(out)], check=True, capture_output=True, text=True)
    assert "All-atom clashes" not in r.stdout and "All-atom clashes" not in out.read_text()      # off by default


def test_portal_passes_the_hydrogens_option():
    p = _portal_module()
    src = open(os.path.join("server", "app.py")).read()
    assert src.count('"--hydrogens"') == 2                                 # interface and quality
    html = p.app.test_client().get("/").get_data(as_text=True)
    assert 'id="hy"' in html and "Also check clashes with hydrogens added" in html


def test_hydrogen_bond_geometry_rules():
    """The geometry test of a hydrogen bond: donor direction, acceptor angle, distance (no hydrogens)."""
    import types
    from ifacekit import contacts
    # donor X at the origin with two neighbours: the expected N-H direction is +x
    xyz = {"X": (0, 0, 0), "N1": (-1.0, 0.6, 0), "N2": (-1.0, -0.6, 0)}
    def check(acc, y):
        m = types.SimpleNamespace(xyz=np.array([xyz["X"], xyz["N1"], xyz["N2"], acc, y], float))
        nb = [[1, 2], [0], [0], [4], [3]]
        return contacts._geometry_ok(m, nb, 0, 3)
    assert check((2.9, 0, 0), (3.6, 1.1, 0))          # in line with N-H, acceptor angle 122 degrees
    assert not check((0, 2.9, 0), (0.7, 3.6, 0))      # 90 degrees from the N-H direction: the donor does not point to the acceptor
    assert not check((2.9, 0, 0), (2.0, 0, 0))        # the acceptor's own neighbour lies between: angle 0 degrees at the acceptor
    # a donor with one neighbour needs the angle N-X...A of at least 100 degrees
    m = types.SimpleNamespace(xyz=np.array([(0, 0, 0), (-1.0, 0, 0), (2.9, 0, 0)], float))
    assert contacts._geometry_ok(m, [[1], [0], []], 0, 2)
    m = types.SimpleNamespace(xyz=np.array([(0, 0, 0), (1.0, 0, 0), (2.9, 0, 0)], float))
    assert not contacts._geometry_ok(m, [[1], [0], []], 0, 2)
    # two atoms on the same position are not neighbours (duplicated chains)
    m = types.SimpleNamespace(xyz=np.array([(0, 0, 0), (0, 0, 0), (1.5, 0, 0)], float), bonds=np.array([[0, 1], [0, 2]]), n_atoms=3)
    assert contacts._neighbours(m) == [[2], [], [0]]


# ---------------------------------------------------------------- NMR ensembles
def test_nmr_without_bfactors_is_not_taken_for_a_prediction():
    """An NMR ensemble has all B-factors 0; that is not a confidence stored as a fraction (9M8X, 20 models)."""
    m = mdl.load(need("9M8X.cif"))
    assert m.flags["n_models_in_file"] == 20
    assert not m.flags["predicted_model"] and not m.flags["bfactor_scaled"]
    assert m.flags["model_source"].startswith("experimental record")


def test_ensemble_every_model_average_and_first_model(tmp_path):
    from ifacekit import ensemble
    ens = ensemble.analyze_ensemble(need("9M8X.cif"), chain="A", workers=2)
    assert ens["n_models_in_file"] == 20 and [r["model"] for r in ens["models"]] == list(range(1, 21))
    av = ens["average"]["rama_favored_pct"]
    vals = [r["rama_favored_pct"] for r in ens["models"]]
    assert av["n"] == 20 and av["min"] == min(vals) and av["max"] == max(vals)
    assert av["mean"] == pytest.approx(sum(vals) / 20, abs=1e-3) and av["sd"] > 0
    # the first model is the complete result of model 0, the same as a single-model run
    single = ifacekit.analyze(need("9M8X.cif"), chain="A", model_index=0)
    assert ens["first_model"]["summary"] == single["summary"]
    assert ens["models"][0]["heavy_borderline"] == single["summary"]["verdicts"].get("borderline", 0)
    # a different model gives different numbers
    assert len({r["heavy_borderline"] for r in ens["models"]}) > 1
    path = tmp_path / "ens.csv"
    ensemble.write_csv(ens, str(path))
    rows = path.read_text().strip().splitlines()
    assert len(rows) == 1 + 20 + 4 and rows[0].startswith("target,model") and rows[-4].split(",")[1] == "mean"
    assert "20 of 20 models" in ensemble.text_table(ens)


def test_ensemble_models_selection_and_pair_mode():
    from ifacekit import ensemble
    ens = ensemble.analyze_ensemble(need("1LCD.cif"), pair=("A", "B"), models=[0, 2, 7])      # model 7 does not exist (3 models)
    assert ens["target"] == "A+B" and ens["models_analysed"] == [1, 3]
    assert all("interface_area" in r and r["interface_area"] > 0 for r in ens["models"])
    with pytest.raises(ValueError):
        ensemble.analyze_ensemble(need("1LCD.cif"), chain="A", models=[9])


def test_ensemble_cli(tmp_path, capsys):
    from ifacekit.cli import main
    out = tmp_path / "e.json"
    csv_path = tmp_path / "e.csv"
    assert main([need("1LCD.cif"), "--chains", "A", "--ensemble", "--ensemble-csv", str(csv_path), "-o", str(out)]) == 0
    text = capsys.readouterr().out
    assert "Ensemble analysis of A: 3 of 3 models" in text and "Model 1\n" in text and "mean" in text
    d = json.loads(out.read_text())
    assert d["mode"] == "ensemble" and len(d["ensembles"]) == 1 and len(d["ensembles"][0]["models"]) == 3
    assert csv_path.exists()
    # without --ensemble a multi-model file is still read as the first model only
    assert main([need("1LCD.cif"), "--chains", "A", "-q", "-o", str(out)]) == 0
    assert json.loads(out.read_text())["runs"][0]["input"]["model_index"] == 0


@pytest.mark.parametrize("entry, n_models", [("1PJW", 15), ("9M8W", 20), ("9M8X", 20)])
def test_ensemble_of_three_nmr_entries(entry, n_models):
    from ifacekit import ensemble
    f = need(f"{entry}.cif")
    ens = ensemble.analyze_ensemble(f, chain="A", workers=2)
    assert ens["n_models_in_file"] == n_models and len(ens["models"]) == n_models
    assert ens["first_model"]["summary"] == ifacekit.analyze(f, chain="A", model_index=0)["summary"]
    for k in ("rama_favored_pct", "rotamer_outlier_pct", "rmsz_bonds"):
        assert ens["average"][k]["n"] == n_models and ens["average"][k]["min"] <= ens["average"][k]["mean"] <= ens["average"][k]["max"]


def test_residue_rmsd_of_an_nmr_ensemble_and_html_plots(tmp_path):
    from ifacekit import ensemble
    f = need("9M8X.cif")
    rows = ensemble.residue_rmsd(f)
    assert len(rows) == 83 and all(r["n_models"] == 20 and r["b"] == r["dev_mean"] and r["dev_min"] <= r["dev_mean"] <= r["dev_max"] and r["dev_sd"] >= 0 for r in rows)
    core = sorted(r["rmsd_backbone"] for r in rows)[len(rows) // 2]
    assert core < 1.5 and max(r["rmsd_backbone"] for r in rows) > 3 * core        # a defined core and mobile ends
    # two identical models: no spread at all
    two = tmp_path / "two.pdb"
    import gemmi
    st = gemmi.read_structure(f)
    while len(st) > 1:
        del st[len(st) - 1]
    st.add_model(st[0])
    st[1].num = 2
    st.write_pdb(str(two))
    z = ensemble.residue_rmsd(str(two))
    assert max(r["rmsd_all"] for r in z) < 1e-3
    from ifacekit.qualityreport import build
    out = tmp_path / "q.html"
    build(f, ["A"], str(out), ensemble=True)
    text = out.read_text(encoding="utf-8")
    assert "Backbone deviation along the sequence" in text and "B-factor along the sequence" not in text
    from ifacekit.htmlreport import build as pair_build
    out2 = tmp_path / "p.html"
    pair_build(need("1LCD.cif"), "A", "B", str(out2), kind="experimental", ensemble=True)
    assert "Backbone deviation along the sequence" in out2.read_text(encoding="utf-8")


def test_parse_models_numbers_start_at_one():
    from ifacekit.ensemble import parse_models
    assert parse_models(None, 20) == [0] and parse_models("first", 20) == [0] and parse_models("all", 3) == [0, 1, 2]
    assert parse_models("3", 20) == [2] and parse_models("1,4,7", 20) == [0, 3, 6] and parse_models("2-5", 20) == [1, 2, 3, 4] and parse_models("2,2,1", 5) == [1, 0]
    for bad in ("0", "21", "x", "3-"):
        with pytest.raises(ValueError):
            parse_models(bad, 20)


def test_model_reports_one_model_all_models_and_separate(tmp_path):
    import re
    from ifacekit.cli import main
    f = need("1LCD.cif")                                   # 3 models
    one = tmp_path / "one.html"
    assert main([f, "--pair", "A", "B", "--models", "2", "--html", str(one), "-q"]) == 0
    text = one.read_text(encoding="utf-8")
    assert "Model 2 of 3" in text and "Backbone deviation along the sequence" in text
    allp = tmp_path / "all.html"
    assert main([f, "--chains", "A", "--models", "all", "--html", str(allp), "--workers", "2", "-q"]) == 0
    page = allp.read_text(encoding="utf-8")
    assert "Ensemble summary" in page and all(f"data-v='{k}'" in page for k in (1, 2, 3)) and "mvp" not in page
    data = json.loads(re.search(r"id='ens-data'>(.*?)</script>", page, re.S).group(1).replace("<\\/", "</"))
    assert sorted(data["m"]) == ["1", "2", "3"] and "Ensemble analysis" not in data["sum"] and "<table" in data["sum"] and "<nav class='top'" in data["m"]["3"]["body"]
    assert page.count("<h1>") == 1 and "<!doctype" not in data["m"]["2"]["body"]
    bar = page.split("class='etabs'>")[1].split("id='ens-body'")[0]
    assert "light / dark" in bar and "class='theme'" in bar          # the switch is in the bar, so it is there on the summary and on every model
    assert "button[data-v]" in page                                    # the buttons of the models are the only ones the page rewires                   # one frame; each model only brings its tables and sections
    p3d = tmp_path / "all_3d.html"
    assert main([f, "--chains", "A", "--models", "all", "--html", str(allp), "--molstar", "--workers", "2", "-q"]) == 0
    t3 = p3d.read_text(encoding="utf-8")
    assert "id='mvp'" in t3 and "id='mv-all'" in t3 and "id='ms-open'" in t3 and t3.count("id='mvp'") == 1
    data3 = json.loads(re.search(r"id='ens-data'>(.*?)</script>", t3, re.S).group(1).replace("<\\/", "</"))
    assert "data-a=" in data3["m"]["2"]["body"]                                                    # rows can be clicked to show the atoms in Mol*
    sep = tmp_path / "s.html"
    assert main([f, "--pair", "A", "B", "--models", "1,3", "--separate", "--html", str(sep), "-q"]) == 0
    assert (tmp_path / "s_model01.html").exists() and (tmp_path / "s_model03.html").exists() and not sep.exists()
    assert main([f, "--chains", "A", "--models", "9", "--html", str(tmp_path / "x.html"), "-q"]) == 2      # no such model
    two = tmp_path / "two.html"
    assert main([f, "--pair", "A", "B", "--models", "2", "--html", str(two), "--molstar", "-q"]) == 0
    assert "Model 2 of 3" in two.read_text(encoding="utf-8") and "id='mvp'" in (tmp_path / "two_3d.html").read_text(encoding="utf-8")


def test_portal_ensemble_option(tmp_path, monkeypatch):
    """The portal reads the number of models, shows the ensemble option for a multi-model file and runs --ensemble --models all."""
    import importlib.util
    pytest.importorskip("flask")
    monkeypatch.setenv("IFACEKIT_DATA", str(tmp_path))
    spec = importlib.util.spec_from_file_location("ifk_portal_ens", os.path.join("server", "app.py"))
    app = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(app)
    html = (Path("server") / "index.html").read_text()
    assert 'id="ew"' in html and "Analyse every model of the file" in html and "n_models" in html
    client = app.app.test_client()
    f = need("1LCD.cif")                                            # 3 models
    with open(f, "rb") as fh:
        r = client.post("/api/inspect", data={"structure": (fh, "1LCD.cif")}, content_type="multipart/form-data").get_json()
    assert r["n_models"] == 3 and r["job"]
    a = client.post("/api/analyze", data={"job": r["job"], "type": "quality", "chains": "A", "ensemble": "1"}).get_json()
    assert "error" not in a, a
    labels = [x["label"] for x in a["files"]]
    assert any("Ensemble table" in x for x in labels) and any("RMSD of each residue" in x for x in labels)
    page = client.get(a["report"]).get_data(as_text=True)
    assert "Ensemble summary" in page and "data-v='3'" in page and "id='mvp'" not in page
    assert a["report_3d"] and "id='mvp'" in client.get(a["report_3d"]).get_data(as_text=True) and "id='mv-all'" in client.get(a["report_3d"]).get_data(as_text=True)
    one = client.post("/api/analyze", data={"job": r["job"], "type": "quality", "chains": "A", "ensemble": "0"}).get_json()
    assert not any("Ensemble table" in x["label"] for x in one["files"])


def test_ensemble_ramachandran_rotamer_and_pdf_show_every_model(tmp_path):
    from ifacekit import rama, rotamer
    from ifacekit.cli import main
    f = need("1LCD.cif")                                          # 3 models
    one, rows1 = rama.analyze(f, 0, ["A"])
    allm, rowsa = rama.analyze(f, 0, ["A"], models=[0, 1, 2])
    assert allm["n"] == 3 * one["n"] and allm["n_models"] == 3 and {r["model"] for r in rowsa} == {1, 2, 3} and "model" not in rows1[0]
    ts, trows = rotamer.analyze(f, 0, ["A"], models=[0, 1, 2])
    assert {r["model"] for r in trows} == {1, 2, 3} and "model 2" in " ".join(rotamer.label(r) for r in trows)
    rp, tp, pdfp = tmp_path / "r.html", tmp_path / "t.html", tmp_path / "p.pdf"
    assert main([f, "--chains", "A", "--ensemble", "--rama-html", str(rp), "--rotamer-html", str(tp), "-q"]) == 0
    assert "in 3 models" in rp.read_text() and "All 3 models of the ensemble" in tp.read_text()
    # without --ensemble the pages are for one model, as before
    rp1 = tmp_path / "r1.html"
    assert main([f, "--chains", "A", "--rama-html", str(rp1), "-q"]) == 0
    assert "models" not in rp1.read_text().split("<h1>")[1].split("</p>")[0]
    # the ensemble page holds the Ramachandran and rotamer plots of every model and the band plot of the deviations
    page = tmp_path / "e.html"
    assert main([f, "--chains", "A", "--models", "all", "--html", str(page), "-q"]) == 0
    text = page.read_text(encoding="utf-8")
    assert "All 3 models of the ensemble in one plot" in text and "Backbone deviation along the sequence" in text and '\\"lo\\":' in text and '\\"s1\\":' in text
