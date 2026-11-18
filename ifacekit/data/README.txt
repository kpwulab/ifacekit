rama_ifacekit.npz  - ifacekit's own Ramachandran reference (built 2026-10-27, tools/rama_reference/).
Source: PDB entries (PDB archive data are CC0). X-ray, resolution <= 1.8 A, one representative per 30% sequence-identity
group (13115 entries listed, RCSB search); entries with R-free <= 0.23 used. Residues: backbone B < 30, occupancy 1,
no alternate location, linked neighbours (3.74 M residues). Six classes (general, Gly, cis-Pro, trans-Pro, pre-Pro, Ile/Val).
Wrapped Gaussian kernel density, sigma 4 degrees, 2-degree grid. Favored / allowed cut-offs = density levels enclosing
98% / 99.95% of the reference residues. Not the MolProbity Top8000 data: those tables are not shipped.

factsheet.html - the fact sheet page (comparison with MolProbity, clashes with hydrogens added, list of test cases), made by tools/make_factsheet.py.

rotamer_ifacekit.npz - ifacekit's own side-chain rotamer reference (built 2026-10-27, tools/rotamer_reference/).
Same 13115 listed entries and R-free <= 0.23 as above. Residues: every atom used in the chi definitions has occupancy 1, no
alternate location and B < 30 (2.66 M residues, 17 tables: Arg, Asn, Asp, Cys, Gln, Glu, His, Ile, Leu, Lys, Met, Phe/Tyr, Pro,
Ser, Thr, Trp, Val; Phe and Tyr share a table; chi2 of Asp and Phe/Tyr and chi3 of Glu folded to 0-180 degrees). Periodic
Gaussian kernel density (smoothing 5 / 5 / 8 / 10 degrees for 1 / 2 / 3 / 4 chi angles, times 1.5), grid layout as in the
MolProbity tables. Favored / allowed cut-offs per type = density levels enclosing a share of the reference residues; the shares are calibrated
against MolProbity's rotalyze on 900 structures with a clean split (design document, section 5.12; tools/rotamer_reference/
enclosure_calibrated.json). The Top8000 data themselves are not shipped.
