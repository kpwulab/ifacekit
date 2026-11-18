"""Hydrogen-bond counts with and without water.

direct          H-bonds between the two sides themselves (no water involved)
water-mediated  waters that make an H-bond to an atom of side A and an H-bond to an atom of side B (a bridge A - water - B)
total           direct H-bonds plus water-mediated bridges, counted as connections between the two sides

H-bond means a donor-acceptor pair with geometry (contact field hbond_geometry, see ifacekit.contacts): X...A of at most 3.3 A, the donor pointing to the
acceptor (angle between the expected N-H direction and X->A at most 50 degrees, or N-X...A of at least 100 degrees for a donor with one neighbour) and every angle
Y-A...X at the acceptor at least 100 degrees. No hydrogens are used. The older distance-only count (type "hbond", X...A up to 3.5 A, no angle) is kept as
direct_polar_contacts. A bridge through two waters is not counted.
A model without waters (an AI prediction) always has zero water-mediated bonds.

Every bridge is also checked with four simple geometric tests on its closest partner atom of each side (the limits below are
uncalibrated, chosen to catch the obvious cases):
  too_short     an H-bond shorter than BOND_MIN: too close for a water, the water may be misplaced or may be another atom
  long_bond     an H-bond longer than BOND_MAX: weak
  triangle      the two partner atoms are closer than PARTNER_MIN to each other: they touch directly, the water only sits beside that contact
  narrow_angle  the angle partner - water - partner is smaller than ANGLE_MIN
A bridge without flags is "clear", otherwise "doubtful". The counts of direct bonds and of all bridges are not changed by the checks.
"""
import math

BOND_MIN = 2.5        # A
BOND_MAX = 3.3        # A
PARTNER_MIN = 3.6     # A
ANGLE_MIN = 80.0      # degrees
FLAG_TEXT = {"too_short": f"H-bond shorter than {BOND_MIN} A", "long_bond": f"H-bond longer than {BOND_MAX} A",
             "triangle": f"partners within {PARTNER_MIN} A of each other", "narrow_angle": f"angle at the water under {ANGLE_MIN:.0f} deg"}


def xyz_lookup(model, atom_id_fn):
    """A function atom id -> (x, y, z) for an ifacekit model; the index is built on first use."""
    cache = {}

    def get(atom_id):
        if not cache:
            for k in range(model.n_atoms):
                cache[atom_id_fn(model, k)] = k
        k = cache.get(atom_id)
        return None if k is None else tuple(float(v) for v in model.xyz[k])
    return get


def _angle(w, a, b):
    v1 = [a[i] - w[i] for i in range(3)]
    v2 = [b[i] - w[i] for i in range(3)]
    n = math.dist((0, 0, 0), v1) * math.dist((0, 0, 0), v2)
    return math.degrees(math.acos(max(-1.0, min(1.0, sum(x * y for x, y in zip(v1, v2)) / n)))) if n > 1e-9 else 0.0


def bridge_flags(side_a, side_b, water_xyz=None, xyz=None):
    """side_a, side_b: lists of (atom id, distance). Returns the list of flags (see the module text)."""
    a = min(side_a, key=lambda t: t[1])
    b = min(side_b, key=lambda t: t[1])
    flags = []
    if min(a[1], b[1]) < BOND_MIN:
        flags.append("too_short")
    if max(a[1], b[1]) > BOND_MAX:
        flags.append("long_bond")
    if xyz is not None and water_xyz is not None:
        pa, pb = xyz(a[0]), xyz(b[0])
        if pa is not None and pb is not None:
            if math.dist(pa, pb) < PARTNER_MIN:
                flags.append("triangle")
            if _angle(water_xyz, pa, pb) < ANGLE_MIN:
                flags.append("narrow_angle")
    return flags


def split(contacts, is_a, is_b, xyz=None):
    """contacts: contact dicts of an analysis result; is_a / is_b: functions telling whether a unit id belongs to side A / B.
    xyz: optional function atom id -> (x, y, z); without it the checks that need coordinates (triangle, narrow_angle) are skipped."""
    direct, direct_polar, waters, to_water = 0, 0, {}, 0
    for c in contacts:
        geo = c.get("hbond_geometry", c["type"] == "hbond")
        polar = c["type"] == "hbond"
        if not (geo or polar):
            continue
        u1, u2 = c["unit1"], c["unit2"]
        if u1 == "water" and u2 == "water":
            continue
        if u1 == "water" or u2 == "water":
            if not geo:
                continue
            w, o, ou = (c["atom1"], c["atom2"], u2) if u1 == "water" else (c["atom2"], c["atom1"], u1)
            e = waters.setdefault(w, {"a": [], "b": []})
            if is_a(ou):
                e["a"].append((o, c["distance"]))
                to_water += 1
            if is_b(ou):
                e["b"].append((o, c["distance"]))
            continue
        if c["inter_unit"] and ((is_a(u1) and is_b(u2)) or (is_a(u2) and is_b(u1))):
            direct += 1 if geo else 0
            direct_polar += 1 if polar else 0
    bridges = []
    for w, e in sorted(waters.items()):
        if not (e["a"] and e["b"]):
            continue
        flags = bridge_flags(e["a"], e["b"], xyz(w) if xyz else None, xyz)
        bridges.append({"water": w, "side_a": [{"atom": a, "distance": d} for a, d in e["a"]],
                        "side_b": [{"atom": a, "distance": d} for a, d in e["b"]],
                        "flags": flags, "quality": "doubtful" if flags else "clear"})
    clear = sum(1 for x in bridges if x["quality"] == "clear")
    return {"direct": direct, "direct_polar_contacts": direct_polar, "water_mediated": len(bridges), "water_mediated_clear": clear,
            "water_mediated_doubtful": len(bridges) - clear, "total_connections": direct + len(bridges),
            "total_connections_clear": direct + clear, "mediated_waters": bridges, "side_a_to_water": to_water}
