"""Compare ifacekit with a PISA interfaces XML for chain pairs inside the model.

usage: python tools/compare_pisa.py data/7d47.cif [data/7d47_pisa_interfaces.xml]
Only interfaces where both molecules are the identity copy (no symmetry mate)
are compared. PISA hydrogen bonds are shown at all distances and at <= 3.5 A.
"""
import sys
import xml.etree.ElementTree as ET

import ifacekit


def pisa_pairs(xml):
    root = ET.parse(xml).getroot()
    out = []
    for i in root.findall("interface"):
        ms = i.findall("molecule")
        def identity(m):
            val = lambda t: float(m.findtext(t) or 0)
            rot = [val(t) for t in ("rxx", "rxy", "rxz", "ryx", "ryy", "ryz", "rzx", "rzy", "rzz")]
            tr = [val(t) for t in ("tx", "ty", "tz")]
            cell = [m.findtext(t) for t in ("cell_i", "cell_j", "cell_k")]
            return rot == [1, 0, 0, 0, 1, 0, 0, 0, 1] and not any(tr) and cell == ["0", "0", "0"]
        if not all(identity(m) for m in ms):
            continue
        a, b = [m.findtext("chain_id") for m in ms]
        if a == b:
            continue
        d = [float(x.findtext("dist")) for x in i.findall("h-bonds/bond")]
        out.append(dict(a=a, b=b, area=float(i.findtext("int_area")), hb=len(d),
                        hb35=sum(x <= 3.5 for x in d),
                        sb=int(i.findtext("salt-bridges/n_bonds") or 0),
                        nres=sum(int(m.findtext("int_nres")) for m in ms),
                        n_total=len(root.findall("interface"))))
    return out


def main():
    cif = sys.argv[1]
    xml = sys.argv[2] if len(sys.argv) > 2 else cif.rsplit(".", 1)[0] + "_pisa_interfaces.xml"
    pairs = pisa_pairs(xml)
    print(f"{cif}: PISA interfaces inside the model: {len(pairs)}")
    if not pairs:
        print("  none (single chain, or only crystal-contact interfaces); nothing to compare")
        return
    print("pair   PISA_A2  ifk_A2   diff%  | hb PISA(all/<=3.5)  ifk | salt PISA  ifk | nres PISA  ifk")
    for p in pairs:
        try:
            r = ifacekit.analyze(cif, pair=(p["a"], p["b"]))
        except ValueError as e:                      # PISA names a molecule that is not a chain of the model (ligand, other label)
            print(f"{p['a']}-{p['b']:<4}{p['area']:8.1f}  not compared: {e}")
            continue
        it = r["interfaces"][0] if r["interfaces"] else None
        area = it["area"] if it else 0.0
        hb = it["contacts"].get("hbond", 0) if it else 0
        sb = it["contacts"].get("salt_bridge", 0) if it else 0
        nres = it["n_residues_a"] + it["n_residues_b"] if it else 0
        print(f"{p['a']}-{p['b']:<4}{p['area']:8.1f}{area:8.1f}{100 * (area - p['area']) / p['area']:7.1f}  "
              f"| {p['hb']:3d} /{p['hb35']:3d}          {hb:3d} | {p['sb']:3d}   {sb:3d} | {p['nres']:3d}  {nres:3d}")


if __name__ == "__main__":
    main()
