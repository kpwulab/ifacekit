"""Shrake-Rupley solvent accessible surface area with a numpy inner loop.

Written here (instead of calling a library) so that the same Bondi radii are
used for SASA, burial and clash tests, and so nucleic acids, ligands and
modified residues need no special naming.
"""
import numpy as np
from scipy.spatial import cKDTree

PROBE = 1.4


def sphere_points(n=960):
    """Evenly spread unit vectors (golden spiral)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


def atom_sasa(xyz, radii, targets=None, env=None, probe=PROBE, n_points=960):
    """Per-atom SASA in square angstrom.

    xyz, radii : coordinates and vdW radii of the *environment* atoms
    targets    : indices (into xyz) of the atoms to evaluate; default all
    env        : unused placeholder for symmetry with older callers
    Returns an array aligned with ``targets``.
    """
    n = len(xyz)
    if targets is None:
        targets = np.arange(n)
    targets = np.asarray(targets, dtype=int)
    out = np.zeros(len(targets))
    if n == 0 or len(targets) == 0:
        return out
    pts = sphere_points(n_points)
    rp = radii + probe
    tree = cKDTree(xyz)
    rmax = rp.max()
    for k, i in enumerate(targets):
        ri = rp[i]
        nb = tree.query_ball_point(xyz[i], ri + rmax)
        nb = [j for j in nb if j != i]
        sp = xyz[i] + ri * pts
        if nb:
            nbx = xyz[nb]
            nbr = rp[nb]
            # squared distance from every surface point to every neighbour centre
            d2 = ((sp[:, None, :] - nbx[None, :, :]) ** 2).sum(axis=2)
            buried = (d2 < (nbr ** 2)[None, :]).any(axis=1)
            frac = 1.0 - buried.mean()
        else:
            frac = 1.0
        out[k] = 4.0 * np.pi * ri * ri * frac
    return out
