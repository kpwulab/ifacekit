---
title: "ifacekit design document"
subtitle: "Interface, contact and clash analysis for PDB and mmCIF models"
author: "Version 0.8.0 (schema 0.12), updated 18 November 2026"
---

# 1. Purpose and scope

ifacekit analyses one model at a time (a multi-model file, such as an NMR ensemble, can be analysed model by model with `--ensemble`). It reports how chains touch each other, where atoms are too close or badly paired, and how good the model looks. It was started as a replacement for PISA that needs no licence, runs on a laptop or a server, and works for both experimental structures and AI predictions.

The tool looks at protein, nucleic acid, ligand, glycan, metal and water contacts inside the model it is given. It does not read or build crystal lattices, and it does not infer biological assemblies. This keeps one engine for an experimental model and for an AlphaFold-style prediction.

The first target use is the review of AI-designed protein binders. For each binder complex the report shows (1) AI confidence metrics, (2) interface interactions and (3) quality of the prediction. A second report type covers experimental structures without any AI metrics.

PISA and MolProbity (including phenix.ramalyze) are used as references for checking only. Neither is called by the tool.

Since version 0.4 the tool has a web portal for internal users (section 7), a Ramachandran analysis (section 3.10) and, from 0.4.7, a rotamer analysis (section 3.11), both with reference distributions that ifacekit built itself, from 0.4.8 rotamer plots and a flip check for Asn, Gln and His (section 3.12), and from 0.5.0 hydrogen-bond counts without and with water (section 3.14), a simpler ligand map (section 3.15), a production set-up of the portal with a limit on simultaneous analyses, a memory limit, clean-up of old jobs and Mol* served by the portal itself (section 7), and the licence and copyright in the footers. The licence is MIT and citation is suggested (section 10).

Version 0.5.5 (28 October 2026) puts the iFaceKit banner on the home page of the web portal and at the top of every HTML report. In the reports the banner is embedded in the page (ifacekit/banner.py, one small WebP image), so a report is still a single file that works offline. The release also has a step-by-step installation guide (INSTALL.txt, see section 8) and the archive is shipped with empty data and results folders. The model type (prediction or experimental) is now also detected from the file for pair reports, as it already was for ligand reports, and in a ligand report of an experimental structure the 3D panel no longer offers the pLDDT colouring; --experimental and --predicted still override the detection. A new option --detailed (and a box in the portal) adds a detailed ligand diagram in which the residues that hydrogen-bond or salt-bridge the ligand are drawn as chemical structures (section 3.16).

Versions 0.7.1 and 0.7.6 (28 October 2026) bring the comparison with MolProbity on 293 structures (section 5.11), the calibrated rotamer cut-offs and the alternate conformers (sections 3.11, 3.17, 5.12), the fact sheet page that the footnotes of the reports link to, and, new in 0.7.6, the optional clash check with hydrogens added (section 3.18, evaluated in section 5.13). ifacekit places the hydrogens itself, and only that check uses them.

Version 0.8.0 (18 November 2026) changes the hydrogen-bond count to a geometry test without hydrogens (section 3.14), adds the analysis of NMR ensembles and other multi-model files (section 3.19: every model analysed, mean and standard deviation, the deviation of each residue from the mean structure, the Ramachandran and rotamer plots with every model, one page whose buttons swap the tables of the models, with the Mol* panel in the 3D version, and the same in the portal), and corrects the reading of NMR files (all B-factors 0) as predictions. Every report shows the version in a box under the banner. The tests use public PDB entries that are downloaded on first use, and the evaluation tables and examples of this document use public entries only (unpublished complexes are named by codes).

# 2. Design criteria

| Criterion | What it means in practice |
|---|---|
| No licence barrier | Only gemmi, numpy and scipy are required. No Phenix, CCP4, CCP4 monomer library or network access is needed to run an analysis. The licence of ifacekit itself is still undecided and must be chosen before it is shared. |
| One engine for both model types | The same code analyses a crystal structure, a cryo-EM model and an AF3 prediction. AI metrics are an extra layer that is only read when the AI files are given. |
| The user chooses what is inspected | With no option the tool only lists the chains. The user then picks single chains, a chain pair, or the whole file. Nothing is analysed that was not asked for. |
| Every verdict has a reason | A contact is never just "bad". It carries a type (for example hbond or salt bridge), a verdict (favorable, ok, borderline, unfavorable) and a reason (for example steric overlap, short hydrogen bond, polar mismatch). The reader can see why. |
| Clashes and interface issues are kept apart | Steric clashes are shown separately from chemical issues such as a too short hydrogen bond or two like charges. Some of the second group may be explained by protonation, which a hydrogen-free tool cannot see. |
| Honest about what is not checked | Every report states that clashes that involve hydrogens are checked only when hydrogens are added (section 3.18), and the notes name what gets no hydrogens, that the flip check for Asn, Gln and His is a weak heavy-atom test (section 3.12), and that rotamer calls on these three are given as modelled and as flipped. Colour cut-offs are described as rules of thumb. |
| Calibrated against a reference, not copied from it | Thresholds were adjusted until heavy-atom problems found by MolProbity were also found by ifacekit, without flagging large numbers of contacts that MolProbity accepts. MolProbity is not part of the pipeline. |
| Self-contained output | A report is one HTML file with the data embedded. It works offline. The optional 3D page needs Mol*, which is loaded from a CDN today and is to be served from the internal server (section 7). |
| Mac and Linux, and a web service | Pure Python, no compiled extras beyond the three packages. The engine returns plain dictionaries that serialise to JSON, and the portal in server/ runs it as a subprocess. |
| Reproducible | Results are deterministic. The JSON output carries the tool version and a schema version. |

# 3. Design

## 3.1 Modules

| File | Role |
|---|---|
| model.py | Reads PDB or mmCIF with gemmi. Builds the atom table (coordinates, radii, element, occupancy, B-factor), residues and units. Flags a predicted model. |
| attach.py | Decides which waters, ligands and glycans belong to which chain. |
| contacts.py | Finds atom pairs within the candidate cutoff and gives each a type, verdict and reason. |
| chem.py | Atom chemistry tables: polar atoms, donors, acceptors, charges, radii. |
| rings.py | Aromatic rings and pi interactions (stacking, T-shaped, cation-pi). |
| sasa.py | Solvent accessible surface area with a Shrake-Rupley style point method. |
| waters.py | Water categories and flags, metal coordination, support by neighbours. |
| geometry.py | Backbone bond lengths and angles against approximate Engh and Huber values, chain breaks. |
| consistency.py | Compares metal sites between chains with identical sequence. |
| core.py | The analysis entry points. Combines the above and returns the result dictionary. |
| report.py | Text reports for the command line. |
| af3.py | Reads AlphaFold 3 summary and full-data JSON. Computes ipSAE, interface PAE and pDockQ. |
| htmlreport.py | The self-contained HTML report, the experimental variant and the Mol* 3D page. |
| banner.py | The iFaceKit banner as an inline image for the top of every HTML report (data/banner.webp). |
| rama.py | Ramachandran analysis: phi and psi, six residue classes, favored, allowed and outlier from the packaged reference grids. |
| ramaplot.py | Ramachandran plots as a stand-alone HTML page and as a block inside the pair report. |
| rotamer.py, rotamerhtml.py | Side-chain rotamer analysis: chi angles for 17 tables, density from the packaged reference grids, favored, allowed and outlier; the HTML block with plots by type, and the stand-alone page. |
| flips.py | Flip check for Asn, Gln and His side chains from the heavy-atom environment. |
| hydrogens.py | Places hydrogens on the protein and checks all-atom clashes (option `--hydrogens`, section 3.18). Used by the clash check only. |
| qualityreport.py | The HTML quality report of chosen chains. |
| data/rama_ifacekit.npz, data/rotamer_ifacekit.npz | The packaged Ramachandran and rotamer references (sections 3.10 and 3.11). |
| cli.py | The command line. |
| server/app.py | The web portal (Flask): upload checks, PDB download, analysis as a subprocess, result files. |
| server/index.html, server/inspect_file.py | The portal home page, and the helper that reads chains and ligands of an uploaded file in a separate process. |
| tools/ | Scripts that compare with PISA, MolProbity and phenix.ramalyze, and tools/rama_reference/ and tools/rotamer_reference/ that build the two references and run the held-out tests, and compare_ramalyze.py, compare_rotalyze.py and compare_reduce.py. Development only. |

## 3.2 Units and attachment

A unit is the thing that is compared with another unit. With `--chains` and `--pair` a unit is a whole chain, including its ligands, glycans and waters. With `--all` each polymer chain is a unit and each ligand residue is its own unit.

Waters, ligands and glycans are attached by structure and not by their chain label. A ligand or glycan covalently bonded to a chain belongs to it, and a glycan tree moves as one group. Anything else with an atom within 4 Å of the chain is attached to it. A ligand that touches two chains appears in both. The option `--attach label` restores the older rule that uses the label only.

For waters and metals the whole model is analysed first and the result is then cut to the chosen scope. This avoids calling a water unsupported only because its partners are in another chain.

## 3.3 Contact engine

Only heavy atoms are used. Atom pairs closer than 4.5 Å are candidates. Pairs that are covalently linked within four bonds are excluded.

| Quantity | Value |
|---|---|
| Hydrogen bond, maximum donor-acceptor distance | 3.5 Å |
| Hydrogen bond, short warning | below 2.55 Å |
| Hydrogen bond, short and unfavorable | below 2.30 Å (2.10 Å if a water is involved) |
| Salt bridge, maximum distance | 4.0 Å |
| Metal contact, maximum distance | 3.2 Å |
| Polar mismatch, unfavorable | below 2.5 Å (candidate up to 3.1 Å) |
| Like charges, unfavorable | below 3.0 Å (candidate up to 4.0 Å) |
| Steric overlap, unfavorable | 0.60 Å or more |
| Steric overlap, borderline | 0.25 to 0.60 Å |

Overlap is the sum of the two van der Waals radii minus the distance. A carboxylate close to another acceptor is reported as possible protonation, because the proton cannot be seen without hydrogens. These contacts are hints and not errors.

The number "unfavorable per 1000 atoms" is a convenient size-corrected count. It is not a MolProbity clashscore and the report says so.

## 3.4 Interface and surface

Surface area is calculated with a Shrake-Rupley point method, Bondi radii and a 1.4 Å probe. The same radii are used for burial and for the clash tests. The interface area is the average of the area buried on each side, which equals half of (area of A + area of B - area of A and B together). Interface residues, most buried residues, polar fraction, buried polar atoms without a partner, and bridging waters follow from the same calculation.

## 3.5 Quality of a model

* Contact problems in the whole model, split into clashes and other issues.
* Backbone bond and angle deviations as Z-scores against approximate Engh and Huber values, and chain break classification.
* Ramachandran categories (favored, allowed, outlier) per residue class, with counts, percentages and the outlier list (section 3.10).
* Rotamer categories (favored, allowed, outlier) per residue type, with counts, percentages and the outlier list (section 3.11).
* Water categories and metal coordination.
* Per-residue pLDDT for a prediction, or B-factor for an experimental model.

## 3.6 AI metrics (AlphaFold 3)

From the summary file: ipTM, pTM, ranking score, chain-pair ipTM, chain-pair minimum PAE, clash flag and fraction disordered. From the full-data file: the PAE matrix, from which ipSAE (PAE cut-off 10 Å, d0 = max(1.0, 1.24 (max(n,27) - 15)^(1/3) - 1.8)), the mean interface PAE and the fraction of interface residue pairs below 5 Å are computed. pDockQ is computed from CB contacts within 8 Å and pLDDT. In an AF3 CIF file the B-factor column holds pLDDT.

## 3.7 Reports

| Output | Content |
|---|---|
| Text | Command line reports for chains, pairs or the whole file. |
| JSON (`-o`) | The full result with a schema version. |
| HTML, predicted complex | Part 1 AI confidence metrics, part 2 interface interactions, part 3 quality of the prediction. |
| HTML, experimental | Part 1 interface interactions, part 2 quality of the model with B-factors. |
| HTML 3D page (`--molstar`) | The same report with a Mol* panel. Clicking a row of any table, or a residue on the plot, shows it in 3D, with a distance line for atom pairs. |

Common features of the HTML pages: a sticky navigation bar with links to each part, summary chips with the key numbers, a "Needs attention" box, tinted section backgrounds, light and dark themes, a PAE heat map, and a zoomable pLDDT or B-factor plot with hover information and interface residues marked.

Colour rules (rules of thumb, not validated limits):

| Metric | Green | Amber | Red |
|---|---|---|---|
| ipSAE | 0.6 or more | 0.3 to 0.6 | below 0.3 |
| ipTM | 0.8 or more | 0.6 to 0.8 | below 0.6 |
| pTM | 0.8 or more | 0.5 to 0.8 | below 0.5 |
| Interface PAE | 5 Å or less | up to 10 Å | above 10 Å |
| Interface pLDDT | 90 or more | 70 to 90 | below 70 |
| Contact problems | none | only borderline | any unfavorable |
| AF3 clash flag | no | | yes |

pDockQ is shown without a colour because its constants have not been verified.

## 3.8 Ligand interaction diagram and report

For a ligand bound to a protein, `--ligand` writes a report page and the figures. Everything is made by ifacekit itself in one command.

**Ligand chemistry.** Files from the wwPDB carry the chemical component dictionary. When it is present, ifacekit takes the bonds, bond orders and aromatic flags from the file, so the structure is exact (the JQ1 in 4QZS and the V00 inhibitor in 8JUX were read correctly this way). Only when the dictionary is missing, for example in a prediction, are the bonds and bond orders guessed from the coordinates (covalent radii, planarity, valence, checked by RDKit). The report states which source was used. Aromatic rings are always drawn as alternating single and double bonds.

**Two views.** The report holds both drawings and a switch between them. The LigPlot-style view has ellipses for residues, dashed hydrogen bonds with distances and dotted lines for the closest other contact. The modern view has coloured cards for residues, a dot for the residue type, curved hydrogen bonds with distance labels, small circles for waters, hover to isolate the contacts of one residue, and a ball-and-stick ligand with soft shading. In both views the ligand is turned to match its 3D orientation, and contacting residues are placed on the side where they are in the model and spread apart so they do not overlap. Waters in contact with the ligand are drawn (see 9.1).

**Figures and downloads.** The command writes SVG, 300 dpi PNG and PDF for each view next to the report. PNG and PDF need the optional package cairosvg, and without it (or without the system library libcairo that it loads) only SVG is written. The HTML page also has download buttons for SVG, PNG 300 dpi and PDF 300 dpi, made in the browser from the view on screen on a white background. The PDF from the page is a 300 dpi raster, the SVG is the vector form.

**Model type.** The confidence column is read as pLDDT for a predicted model and as B-factor for an experimental one. The type is detected from the file (prediction keywords such as AlphaFold, Boltz or Chai, a confidence stored as a 0 to 1 fraction, or the experimental records of the wwPDB), and a file with no sign of either is treated as experimental. The report says how it decided. `--experimental` and `--predicted` override the detection.

**Tables.** The page adds per-atom burial, pocket residues, buried polar atoms without a partner and ligand-protein clashes, all linked to Mol* in the 3D version. RDKit is an optional extra and is needed only for this report.

## 3.9 Command line

```
ifacekit model.cif                       list the chains
ifacekit model.cif --chains A C          quality of single chains
ifacekit model.cif --chains all          every chain with a comparison table
ifacekit model.cif --pair A B -o r.json  interaction between two chains
ifacekit model.cif --all                 everything
ifacekit model.cif --pair A B --html r.html [--af3-summary S --af3-full F]
ifacekit model.cif --pair A B --html r.html --experimental
ifacekit model.cif --chains A --html q.html                quality report of chain A (HTML, same look as the other reports)
ifacekit model.cif --chains A --rama-html rama.html           also a stand-alone Ramachandran page
ifacekit model.cif --chains A --rotamer-html rot.html         also a stand-alone rotamer page (plots, outliers, flip check)
ifacekit model.cif --chains A C --plots-pdf plots.pdf         one PDF with all Ramachandran and rotamer plots (also with --pair)
ifacekit model.cif --pair A B --html r.html --molstar --molstar-url /static/molstar/   3D page that loads Mol* from your own server
ifacekit model.cif --ligand [ID] --html lig_report.html    ligand report, plus lig_map_modern and lig_map_classic as SVG, PNG, PDF
  (ID is the ligand unit such as A:V00401, the first ligand is used when it is left out)
  --experimental or --predicted overrides the detected model type
add --molstar to also write r_3d.html
```

## 3.10 Ramachandran analysis

Phi and psi are calculated for every amino acid whose neighbours are linked by a peptide bond (C to N distance under 2.0 A). For alternate conformations the first conformer is used. A residue belongs to one of six classes, the same as in MolProbity so that numbers can be compared: general, glycine, cis-proline (omega within 30 degrees of 0), trans-proline, pre-proline (the next residue is proline) and isoleucine or valine. The density of the class is read from a 2 degree grid with bilinear interpolation on the periodic surface. A residue is favored above the density level that encloses 98% of the reference residues of its class, allowed above the level that encloses 99.95%, and an outlier below. The goals shown with the result are over 98% favored and under 0.2% outliers.

The reference grids were built by ifacekit from the PDB, not taken from another program. The steps (scripts in tools/rama_reference/, all inputs and intermediate files kept in a folder called rama_source, here /your_path/rama_source):

1. The RCSB search service gave X-ray entries with resolution 1.8 A or better, R-free 0.25 or lower and protein chains, one representative per 30% sequence identity group: 13,115 entries. The files were downloaded from RCSB and PDBj. The PDB archive is CC0, so the derived reference belongs to ifacekit.
2. Entries with R-free above 0.23 were dropped, read from each file. Residues were kept when the backbone atoms N, CA, C and O had occupancy 1, no alternate location and B below 30, and both neighbours were linked: 3.74 million residues (general 2,655,304; glycine 284,892; cis-proline 9,495; trans-proline 163,878; pre-proline 153,794; Ile/Val 476,295).
3. Per class a periodic Gaussian kernel density with sigma 4 degrees was calculated on the 2 degree grid, and the two cut-off levels were set from the 98% and 99.95% enclosure.

The grid file also stores these settings, and the result of every run names the reference used. Another grid file can be selected in code with `ifacekit.rama.use_reference(path)`.

Output: a block `ramachandran` in the JSON (counts, percentages, counts by class, list of outliers), a text block in the chain quality report, the option `--rama-html FILE` for a stand-alone page with one plot per class, and a section in the quality part of the pair report (chains of the pair only). The plot shades the favored and allowed regions and shows every residue as a point with its name on hover.

## 3.11 Rotamer analysis

Chi angles are calculated for every amino acid except Ala and Gly, for every alternate conformer separately (before 0.7.1 only the first conformer was used). The definitions and the 17 tables follow MolProbity so that numbers can be compared: Arg and Lys (four chi angles), Gln, Glu and Met (three), Asn, Asp, His, Ile, Leu, Trp and Phe/Tyr (two; Phe and Tyr share one table), Cys, Pro, Ser, Thr and Val (one). For the symmetric end groups chi2 of Asp and Phe/Tyr and chi3 of Glu are folded to 0 to 180 degrees. The density at the chi angles is read from the periodic grid by multilinear interpolation. A residue is favored at or above the favored level of its type, allowed at or above the allowed level, and an outlier below. The goal shown is under 0.3% outliers. No hydrogens are used. Asn, Gln and His are scored as modelled and also with the group flipped (the chi angle that changes by 180 degrees is chi2 for Asn and His, chi3 for Gln), and the summary also gives the result after the flips that the flip check suggests with confidence (section 3.12).

The reference was built the same way as the Ramachandran reference, from the same 13,115 entries with R-free 0.23 or lower (scripts in tools/rotamer_reference/). A residue is kept when every atom of its chi definitions has occupancy 1, no alternate location and B below 30: 2.66 million residues (for example Leu 337,449, Arg 123,773, Lys 79,761, Met 54,472, Trp 63,545). Per type a periodic Gaussian kernel density was calculated on a grid with the layout of the MolProbity tables (10 degree bins for four angles, 8 degrees for three, 5 degrees for two, 1 degree for one). The smoothing widths are 5, 5, 8 and 10 degrees for one to four angles, multiplied by 1.5 (see section 5.8).

The cut-offs are not the same fixed number for all types. MolProbity defines them as density relative to the peak (0.02 favored, 0.003 allowed), and the share of the reference residues that this encloses differs strongly by type (for example 78% favored for Lys, 99% for Pro). Fixed peak-relative cut-offs on our own grids gave 84% agreement with phenix.rotalyze, because the result depends on the smoothing and on how much data a type has. Therefore the cut-offs are density levels that enclose a share of the reference residues of the type. Up to version 0.6 the two shares per type were taken from MolProbity's Top8000 tables. From 0.7.1 they are calibrated against rotalyze with a clean split of the data (section 5.12): the smoothing factor of 1.5 was kept, and for 13 of the 17 types the two shares were set on 900 structures; types with fewer than 40 rotalyze outliers in the calibration set (Cys, His, Pro, Trp) keep the earlier shares. The shares are in `tools/rotamer_reference/enclosure_calibrated.json`. The grids themselves are ifacekit's own.

Alternate conformers: every conformer of a residue is scored on its own, as MolProbity does (the atoms of that conformer plus the atoms shared by all conformers). A row of the result has the field `altloc` (the letter, empty for a residue without alternate conformers), the counts count conformers, and the reports show "alt B" behind the residue.

Output: a block `rotamer` in the JSON (counts, percentages, counts by type, outlier list with chi angles), a text block in the chain quality report, and a section with a table by type and the outliers in the HTML quality report (per chain) and the pair report (chains of the pair). Rotamer plots (below) are in the HTML reports and in a stand-alone page (`--rotamer-html FILE`). Rotamer names (such as mt or tp) are not assigned.

The plots have one panel per residue type. For types with two or more chi angles the panel shows chi1 against chi2, shaded with the highest reference density over the remaining chi angles (favored dark, allowed light), so for three and four chi angles the shaded region is larger than the true region in the full space, and every residue is a point coloured by its category in the full space. For types with one chi angle the panel shows the reference density along chi1 with the two cut-off levels and the residues on the curve.

## 3.12 Flip check for Asn, Gln and His

The end atoms of an Asn or Gln amide (O and N) and the two sides of the His ring cannot be told apart in an electron density map, so a model can have the group the wrong way round. For every such group the environment is scored twice, as modelled and with the group turned by 180 degrees (O and N positions swapped; for His the four ring positions swapped). The score sums over the polar neighbours of the two end atoms: a good hydrogen bond adds up to 1 (distance within 3.2 A, or up to 3.5 A for half a point), with the hydrogens of an amide nitrogen placed from the geometry of the group (C-N-H 120 degrees, in the plane), the lone-pair side of an amide oxygen and the ring nitrogens treated with angle limits; a too close contact, two acceptors or two donors facing each other, and a non-polar atom pressing on a polar end atom each subtract 1. The environment is the whole model, first conformer, with waters and ligands, but without crystal contacts. The verdict is keep, uncertain (gain of 1 or more for the flipped group) or flip (gain of 2 or more). These thresholds were set on the training set of section 5.9.

This is much simpler than Reduce, which adds hydrogens, scores contacts with the full environment and optimises the whole hydrogen-bond network, and it is weaker (section 5.9). The tool therefore only suggests flips and says so in every place where they appear (JSON block `flips`, text report, HTML section). The flipped alternative is also given for the rotamer score.

## 3.13 Plots PDF and the chain labels of the sequence plot

`--plots-pdf FILE` writes one A4 PDF with all Ramachandran and rotamer plots of the chosen chains (the two chains of a pair with `--pair`). The pages are a summary page (counts for every chain, the references used, a short key), then for every chain its Ramachandran pages, its rotamer pages (up to six plots per page, two columns) and a page that lists the outliers and the suggested flips. The plots are the same SVG panels as in the HTML reports, with the page colours written out, converted with cairosvg and joined with pypdf, so the PDF needs the optional packages (`pip install ".[pdf]"`, which also needs the system library libcairo). When the PDF is written next to an HTML report, the report gets a download link to it, shown as a larger outlined button directly under the header of each Ramachandran section (with its own colours, so that it is readable in light and dark mode); when it cannot be written (a package is missing) the tool says so and the report has no link.

In the B-factor or pLDDT plot along the sequence each chain has its own colour for the line and, under the axis, a thick bar of the same colour with the chain name in a large bold label ("Chain D", shortened to the letter when the chain is narrow on the screen), and a thick divider between chains. This holds for the pair report and the quality report.

## 3.14 Hydrogen bonds without and with water

Hydrogen-bond numbers are given in three parts (module hbonds.py, JSON field `hbonds` of every interface). Since 0.8 a hydrogen bond is a donor-acceptor pair with geometry (contact field `hbond_geometry`): X···A of at most 3.3 Å, the donor pointing to the acceptor (angle between the expected N-H direction, minus the sum of the unit vectors to the donor's heavy neighbours, and X→A at most 50°; a donor with one neighbour needs N-X···A of at least 100°) and every angle Y-A···X at the acceptor at least 100°. Waters have no neighbours and pass; atoms of ligands with an uncertain role are judged by the distance only; atoms on the same position are not neighbours. The limits were set by hand after the distance-only count had turned out to be about 1.4 times the count of PISA, and compared with PISA, ChimeraX and PyMOL on 279 interfaces (manuscript Table 13: mean bonds per interface 8.3 for ifacekit, 7.1 for PISA, 6.1 to 6.7 for ChimeraX, 7.3 for the PyMOL assessment). The contact type hbond (complementary roles, up to 3.5 Å, no angle) and its verdicts are not changed; the pairs of that type are now called polar contacts and are counted as `direct_polar_contacts`.

* **direct**: H-bonds between the two sides themselves, with no water involved (a hydrogen bond with geometry between an atom of side A and an atom of side B); `direct_polar_contacts` is the older distance-only count;
* **via water**: the number of waters that make an H-bond to an atom of side A and an H-bond to an atom of side B (a bridge A - water - B); the waters are listed with their partner atoms and distances;
* **connections in total**: direct plus via water, counted as connections between the sides, not as bonds (a bridge has two bonds but counts once).

In the pair report the interface cards show the three numbers, the chips show direct and via water, and a table lists the bridging waters. In the ligand report the cards show hydrogen bonds from the ligand to the protein without water, H-bonds from the ligand to waters, and H-bonds via water (a water bonded to the ligand and to the protein), the table of polar contacts has a column that tells whether the partner is a protein or a water, and the bridging waters have their own table. The text report prints the same counts. **Check of every bridge.** A water that is H-bonded to both sides is not always a real bridge, so each bridge is tested on its closest partner atom of each side with four simple tests (limits in `ifacekit/hbonds.py`): `too_short` (an H-bond shorter than 2.5 A: too close for a water, the water may be misplaced or be another atom), `long_bond` (longer than 3.3 A: weak), `triangle` (the two partner atoms are closer than 3.6 A to each other: they touch directly and the water only sits beside that contact) and `narrow_angle` (the angle partner - water - partner is under 80 degrees). A bridge without flags is **clear**, otherwise **doubtful**. The JSON gives for every bridge `flags` and `quality`, and for the interface `water_mediated_clear`, `water_mediated_doubtful` and `total_connections_clear` (direct plus clear bridges); `water_mediated` and `total_connections` are not changed and still count all bridges. The table of bridges lists the bridging waters only, that is waters bonded to both sides; a water bonded to one side only is not an interface bond, is neither listed nor counted there, and belongs to the hydration of that chain (model quality part; in the ligand report a water bonded to the ligand only is counted under "H-bonds to waters" and drawn in the map). The pair and ligand reports show "7 clear, 4 doubtful" on the card and chip, the table of bridges has a column "check" with the reason, and the text report prints the same. The limits are uncalibrated judgement calls (section 5.10): a doubtful bridge is a reason to look at the density, not a verdict.

A bridge through two waters is not counted. A model without waters (an AI prediction) cannot have bonds via water, and every report says so instead of showing only a zero. The older list of waters that touch both sides (`bridging_waters`) is kept; it asks only for contact within range, not for an H-bond to each side. The quality report does not show H-bond counts and is unchanged.

## 3.15 Simple and full ligand map

The interaction map of the ligand report is drawn in two levels. The simple level (default) draws a water only when it makes an H-bond to the ligand, and a residue only when it has a polar contact (H-bond, salt bridge, metal) or two or more contacts of which the closest is within 3.8 A. The full level draws everything within 4 A, as before. The report has a switch between them next to the switch between the modern and the LigPlot-style view, the legend says how many residues and waters the simple map leaves out, and the tables always list every residue. The command line writes the figure files of both levels (`..._map_modern`, `..._map_classic` for the simple map and `..._map_all_modern`, `..._map_all_classic` for the full map, as SVG, PNG and PDF). On 4QZS (JQ1) the simple map has 12 residues and waters instead of 22 (8 waters and 2 residues left out), on protein-ligand complex B0001 13 instead of 19. The rule is a heuristic chosen to remove clutter; it was judged by eye on these two examples.

## 3.16 Detailed ligand diagram

With `--ligand --detailed`, or the box "Also draw the detailed diagram" in the portal, the ligand report gets a third view, Detailed, next to Modern and LigPlot-style, with the same switch between Simple and All contacts. The command also writes the figure files `*_map_detailed` and `*_map_all_detailed` (SVG, and PNG and PDF when libcairo is present). A residue that hydrogen-bonds or salt-bridges the ligand is drawn as its chemical structure: the side chain from the C-alpha (grey dot), plus the backbone atoms that take part in a contact, built from fixed templates of the 20 standard amino acids. The contacting atom is ringed and joined to the ligand atom by a dashed line with its distance. Other residues stay boxes, waters stay small circles, and residues without a template (modified residues such as the 3-fluorotyrosine of 4QZS) stay grey boxes. The report without the option is unchanged. With `--detailed-all` (a second box in the portal, off by default) every standard residue of the diagram is drawn as a structure, also those with only van der Waals or hydrophobic contacts; these get a dotted line to their closest contact atom and no ring. The result is more crowded.

Layout. The ligand is drawn flat, as in the other views, and a residue is placed on the side where it is in 3D. If the line to its ligand atom would cross a bond of the ligand or pass close to another ligand atom (a contact seen from behind the ligand), the residue is moved round the ligand by the smallest angle that gives a clear line. The pictures are then pushed apart by their real size and kept off the ligand, and each structure is turned, or mirrored, so that the contacting atom faces the ligand and its atoms keep clear of the ligand and of the pictures already placed. Not done: the depth of a residue (near or far from the viewer) is not shown, lines of different residues are not checked against each other, and the ligand is drawn flat and not as the ball-and-stick of the Modern view. Checked by eye on 4QZS (JQ1) and 8JUX (V00) and by two automatic tests; the layout has not been tried on many ligands.

## 3.17 Alternate conformers in the contact analysis

The contact analysis works on one conformer of every atom (the first, as before), which gives the interfaces, areas and bonds. Since 0.7.1 it also checks the other conformers: for every further conformer letter the model is loaded with the atoms of that conformer plus the atoms shared by all conformers, and the contacts that involve at least one atom of that conformer are added to the result (the pairs between shared atoms were already evaluated). These contacts carry the field `altloc` in the JSON and the contact and issue lists, are counted in the verdict totals (`summary.alt_conformer_contacts` gives their number), and appear in the clash tables of the reports with a column "conformer". `flags.alt_conformers_scanned` lists the letters that were scanned. The option `alt_conformers=False` of `ifacekit.analyze` switches it off. Interface areas, hydrogen-bond counts, water analysis, Ramachandran and the flip check still use the first conformer only. A clash between two atoms of different conformers (A of one residue with B of its neighbour) is not an event that can occur and is not evaluated.

## 3.18 Clashes with hydrogens added

Option `--hydrogens` (box "Also check clashes with hydrogens added" in the portal, interface and chain reports; module hydrogens.py, JSON block `hydrogens`, summary fields `all_atom_clashes` and `all_atom_clashscore`, schema 0.11). It is off by default. **Only this check uses the hydrogens; every other result (contacts, interfaces, rotamers, flips, waters, geometry) comes from the file as given**, and the heavy atoms are never moved. Any hydrogens in the file are ignored (the loader drops them), so all models are judged in the same way. ifacekit places the hydrogens itself and calls no other program.

**Placement.** The 20 standard amino acids, the backbone amide and the N terminus (NH3+, or NH2+ for Pro) get hydrogens at ideal positions: aromatic and sp2 groups in the plane, CH2 and CH groups tetrahedral, NH2 groups of Asn, Gln and Arg planar. The X-H distances are the shortened X-ray values that Reduce uses for X-ray models (C-H 0.97, N-H 0.86, O-H 0.84, S-H 1.20 A). A cysteine in a disulfide gets no thiol hydrogen. Waters, ligands, nucleic acids and modified residues get no hydrogens, so their clashes with hydrogens are not seen (a limit that the reports state).

**Alternate conformers** are not scanned by this check: only the first conformer of every atom is used (the contact tables of the reports do scan the other conformers, section 3.17).

**Rotating groups.** Methyl groups are staggered and not turned (Reduce does the same: in a check on 129 methyls every one was exactly staggered). The hydroxyl of Ser and Thr and the thiol of Cys are turned in steps of 15 degrees, the NH3 of Lys and the N terminus in steps of 10 degrees, the hydroxyl of Tyr has two positions (0 and 180 degrees); each group takes the position with the lowest contact energy against everything around it, in three passes over the model.

**Flips.** An Asn or Gln amide has two states (the NH2 on the atom labelled N, or on the one labelled O, the coordinates unchanged); a His ring has six (two orientations times HID, HIE and HIP). The state is chosen by the same contact energy: overlaps are penalised, a hydrogen bond (polar H, acceptor, distance and an angle at the H of 110 degrees or more) is rewarded by 2 units, HIP costs 1 unit, and a flip must gain 4 units over the given orientation (the value that agrees best with Reduce's choices, see section 5.13). The flips chosen are listed in the result (`hydrogens.flips`), but the rest of the report does not use them.

**Clash rule.** An overlap (sum of the radii minus the distance) of 0.4 A or more between atoms more than four bonds apart is a clash; the bond rule counts the hydrogen as one bond, and two hydrogens are also left out when their heavy atoms are three bonds apart or fewer. Radii are those of the Probe program that MolProbity uses (aliphatic H 1.22, aromatic and polar H 1.05, C 1.70, aromatic C 1.75, N 1.55, O 1.40, S 1.80 A). A polar hydrogen next to an acceptor at a hydrogen-bond angle is not a clash up to an overlap of 1.0 A (1.2 A if one of the two atoms is charged); two N or O atoms are not a clash below 0.55 A. The clashscore is the number of clash pairs per 1000 atoms including the hydrogens. The limits were fitted against MolProbity's clash lists (section 5.13) and are constants at the top of hydrogens.py.


## 3.19 Multi-model files and NMR ensembles

Default: one model of the file is read (model 0, or `--model N`). With `--ensemble` the module `ifacekit/ensemble.py` runs the normal analysis (`core.analyze`) on every model, for one chain, one pair of chains or the whole file; the models are independent, so they can be analysed in parallel (`--workers`).
For each model a flat row of numbers is taken from the result: Ramachandran favored, allowed and outlier percent, rotamer favored and outlier percent, bond and angle RMSZ, geometry outliers, heavy-atom unfavorable and borderline contacts, the all-atom clashes and clashscore when the hydrogen check is on,
the suggested flips, and for interfaces the summed area, the direct hydrogen bonds and the number of interfaces. The average block holds, for each number, the mean, the sample standard deviation (n - 1), minimum, maximum and n. Reports of the models: `--models LIST` (first, all, or model numbers counted from 1 as the file does) chooses the models that get an HTML report (`ifacekit/ensemblereport.py`). The report of one model is the normal chain-quality or interface report of that model with the plots of the ensemble and a note "Model i of N". For several models the page keeps the frame of the normal report (banner, title, notes) and has a bar with a button for the ensemble summary and one for each model. The part of a report that belongs to one model (the navigation, the overview, the tables and sections, up to the notes at the end) is marked in the HTML (`BODY_START`, `BODY_END`), taken from the report of each model, stored as JSON in the page and put in place of the visible part when a button is pressed (its scripts are run again). The summary holds the table of all models with mean and sd, the band plot, and the Ramachandran and rotamer plots with every model in one plot. With `--molstar` the page is also written with one Mol* panel (`FILE_3d.html`): the structure file with all models is embedded once, the panel follows the model chosen on the page (`IFK.setModel`), a click on a table row shows those atoms, and the button "all models" overlays every model with a transparent cartoon. `--separate` writes one file per model instead. The models are rendered in parallel (`--workers`).
The Ramachandran and rotamer analyses take a list of models (`models=`): the rows of all the models are analysed together and each row carries its model number, so that the Ramachandran page, the rotamer page and the plots PDF of an ensemble show every model in one plot (the counts are residue-model pairs). The spread of each residue (the plot that replaces the B-factor / pLDDT plot in the HTML reports of an ensemble: a line with the mean over the models of the distance of the residue from the mean structure, a band of one standard deviation and a band of the range) is computed by `residue_rmsd`: every chain is superposed on its own on the CA atoms of the least mobile 70% of its residues (iterated three times, Kabsch), the mean position of each heavy atom over the models is taken, and the RMSD of a residue is the root of the mean squared distance of its atoms to their mean positions, for the backbone atoms (N, CA, C, O; the plotted value), for CA alone and for all heavy atoms. Chains are fitted separately, so the movement of one chain against another is not part of the number.
The result of the first model is kept in full (`first_model`) and is what the text report, the HTML and the
other pages show. A mean is the mean of per-model results; no averaged structure is built, and contacts between models are not analysed. Choosing the number for a figure or a table (mean, one model) is left to the user; the ensemble table makes the spread visible.
NMR models have no B-factors (all 0); the rule that reads a B-factor column in the range 0 to 1 as a fraction confidence (Boltz, Chai) therefore needs a non-zero B-factor and no experimental record in the file, so that an NMR ensemble is not taken for a prediction.
Tests: 1LCD (3 models, protein and DNA), 1PJW (15 models), 9M8W and 9M8X (20 models each).

# 4. Functions at a glance

1. Chain overview and per-chain quality, with waters, ligands and glycans attached by structure.
2. Pair analysis with interface area, hydrogen bonds, salt bridges, hydrophobic and van der Waals contacts, pi interactions, buried polar atoms and bridging waters.
3. Clash and issue lists with a reason for each entry.
4. Water and metal analysis, including a metal consistency check between identical chains.
5. Backbone geometry check.
6. AlphaFold 3 metrics, ipSAE, interface PAE, pDockQ, PAE map and pLDDT plot.
7. HTML reports for predicted and experimental complexes, with an optional Mol* 3D page.
8. A Python interface (`ifacekit.analyze`, `ifacekit.list_chains`) and JSON output.
9. Ramachandran analysis with own reference, text, JSON and plots.
9b. Rotamer analysis with own reference, text, JSON and HTML tables.
9c. Alternate conformers: rotamers and contacts are scored for every conformer.
10. HTML quality report for chosen chains (`--chains ... --html`), with one section per chain.
11. One PDF with all Ramachandran and rotamer plots (`--plots-pdf`), linked from the reports and the portal.
12. Hydrogen-bond counts without and with water, in the pair and ligand reports and the JSON.
13. A simple and a full ligand map.
14. A web portal with four analysis types, file upload with checks, and download of an entry from RCSB, PDBe or PDBj.

# 5. Evaluation

## 5.1 What was tested

The test suite has more than 100 tests. They use public PDB entries that are downloaded from RCSB on first use (and skipped without internet). They cover SASA against freesasa with the same radii, pairwise interface area against freesasa, equal results for PDB and mmCIF input, chain and pair modes, attachment, waters, metal consistency, backbone geometry, the protonation rule, the ipSAE limits, the ligand diagram on protein-ligand complex B0001 (rebuilt formula and contact counts), the exact ligand chemistry read from the file dictionary for JQ1 in 4QZS, aromatic rings drawn as alternating bonds in both views, detection of prediction versus experiment with the matching labels and the user override, the figure files written by the command, and the writing of the HTML reports (including that the plain report contains no Mol* content), the Ramachandran categories against the cut-off rule, the Ramachandran block in the pair report, the rotamer reference (chi definitions, category rule, a known favored and a strained Leu rotamer), the rotamer section, plots and flip table in the HTML quality report, the flip check (verdict rule, flipped chi), and the HTML quality report (a section per chain, no Mol* content). The tests that use real structures (for example 1A8O and 7CFN from the PDB) read them from a folder given by the environment variable IFACEKIT_TESTDATA and skip themselves when the file is missing. The folder is included in the archive as testdata.

## 5.2 Interface area against PISA

PISA was run from CCP4 on the same files. Only interfaces between two identity copies inside the model were compared (no symmetry mates).

| Model | Pair | PISA area (Å²) | ifacekit area (Å²) | Difference |
|---|---|---|---|---|
| 7d47 | B-A | 515.4 | 525.9 | +2.0 % |
| 9lqk | B-A | 540.7 | 547.5 | +1.3 % |
| 5hpt | B-A | 927.9 | 917.5 | -1.1 % |
| 5hpt | C-A | 945.4 | 943.8 | -0.2 % |
| 5hpt | D-A | 978.8 | 936.8 | -4.3 % |
| 5hpt | C-G | 226.2 | 212.8 | -5.9 % |
| 5hpt | B-D | 56.2 | 50.5 | -10.2 % |
| Protein complex A0001 (AF3) | A-B | 1150.1 | 1158.3 | +0.7 % |
| Protein complex A0002 (AF3) | A-B | 1134.9 | 1132.8 | -0.2 % |

Large interfaces agree within about 1 to 2 %, with one larger deviation of 4.3 %. Small interfaces differ by 6 to 10 %. Hydrogen bond and salt bridge counts are close but not identical, because PISA and ifacekit use different distance limits and different rules for which atom pairs count. 7yqm has no interface inside the model, so there was nothing to compare. The comparison script stops with an error on interfaces that involve a ligand molecule (seen with 7d47 and 9lqk beyond the first interface), which is a limit of the script.

## 5.3 Clashes against MolProbity

MolProbity (Phenix 2.0) was run on five structures. Its severe clashes were separated into heavy-atom pairs and pairs that involve a hydrogen. Matching is done by residue pair, because symmetric atoms such as CD1 and CD2 can be named differently.

| Model | MolProbity severe clashes (heavy atom / with H) | Heavy-atom clashes whose residue pair ifacekit also flags | Clashes with H where ifacekit flags a residue anyway | ifacekit unfavorable pairs not in the MolProbity list |
|---|---|---|---|---|
| 5hpt | 97 / 334 | 77 of 97 | 206 of 334 | 5 |
| Protein complex A0003 (with Mg, ADP, water) | 12 / 65 | 12 of 12 | 2 of 65 | 4 |
| Protein complex A0001 | 9 / 43 | 9 of 9 | 16 of 43 | 0 |
| Protein complex A0002 | 12 / 54 | 12 of 12 | 30 of 54 | 0 |
| 9lqk | 4 / 2 | 3 of 4 | 2 of 2 | 0 |

The heavy-atom agreement is good. The clashes that involve hydrogens are mostly not seen, which is expected for a heavy-atom method. The weakest case is protein complex A0003, where only 2 of 65 hydrogen clashes are caught. The 5hpt model has 20 heavy-atom clashes that ifacekit does not flag at the residue pair level.

An earlier version flagged 122 unfavorable pairs on the A0003 model. After the thresholds were recalibrated it flags 16. The thresholds were therefore tuned on the same five structures that are used here for comparison. The agreement above is not an independent test.

## 5.4 Backbone geometry

The residues with the largest backbone deviations were compared with the worst residues listed by MolProbity and matched. This was checked on the same small set and was not turned into a numerical score.

## 5.5 Reports

The HTML reports were opened in Chrome and Safari on a Mac and checked in headless Chrome on Linux, including the interactive plot, navigation, light and dark themes, and the Mol* panel (click on a hydrogen bond row, interface view and pLDDT colouring). Windows browsers were not tested.

## 5.6 Ligand reports

Four ligand cases were run with the command line and inspected by eye.

* Protein-ligand complex B0001 (AI prediction, no dictionary). The bond orders were guessed and gave the right formula. Three hydrogen bonds were found (2.94, 2.59 and 3.25 A) among 19 pocket residues. Detected as a prediction because the confidence is stored as a fraction.
* 4QZS, BRD4 with JQ1 (X-ray). The structure was read from the file dictionary and matched the known JQ1 skeleton and formula (C23H25ClN4O2S). One hydrogen bond to Asn140 (3.16 A) and two to waters (2.83 and 2.92 A). Detected as experimental. The modified residue 3-fluorotyrosine (YOF) is listed with the residue code and is not given a residue type, so it is grey in the picture.
* 8JUX, SARS-CoV-2 papain-like protease with the inhibitor V00 (X-ray). Both copies (A and B) were run, and a file cut to chain A alone. Hydrogen bonds to Asp164 (3.26 A) and Tyr268 (3.04 A). No clashes in the ligand.
* The ligand charts were not compared with LigPlot+ or PLIP output yet. A PDBsum LigPlot for 4QZS could not be fetched, so the comparison is still to do.

## 5.10 Hydrogen bonds without and with water

Run on the example structures (direct: no water; via water: waters bonded to both sides, split into clear and doubtful by the check of section 3.14):

| Structure, interface | direct | via water | clear | doubtful | note |
|---|---|---|---|---|---|
| 5hpt, chains D and E | 4 | 0 | 0 | 0 | no bridging waters |
| Protein complex A0002 (AI prediction), chains A and B | 8 | 0 | 0 | 0 | no waters in the model |
| 9lqk, chains A and B | 6 | 11 | 7 | 4 | doubtful: A:201 (2.06 A bond), B:102 (2.33 A), A:221 and A:215 (triangles) |
| 7d47, chains A and B | 0 | 6 | 1 | 5 | no direct H-bond and five of the six bridges flagged, as expected: this interface is crystal packing, not a biological assembly |
| Protein complex A0003, chains F and I | 24 | 18 | 14 | 4 | |
| 4QZS, BRD4 and ligand JQ1 | 1 | 2 | 2 | 0 | JQ1 also bonds two waters directly, which are the two bridging waters |
| 8AFB, ligand LXD203 | 8 | 2 | 2 | 0 | 3 H-bonds to waters |

The direct count equals the number of H-bonds that the earlier version of the pair report gave, in every case tested. A test checks the rule on constructed contacts (a direct bond, a bridging water, a water that touches one side only). The earlier ligand report counted waters among its hydrogen bonds (for 4QZS it gave 3, of which 2 were to waters); it now gives 1 plus the two water numbers.

Correction after version 0.5.5 and accuracy check of the protein-protein counts. In the HTML pair report the table of polar contacts and the counts of interface clashes and issues had wrongly included all contacts with waters (a water is a different unit from the chains, so it passed the test "between units"). With the fix the table of 9lqk A-B has 7 rows instead of 204 and that of 7d47 A-B 2 instead of 405; the interface clashes and issues of 7d47 A-B drop from 4 and 46 to 0 and 2. The cards, the JSON and the text report were not affected. The protein-protein counts were checked against PISA with tools/compare_pisa.py (identity-copy interfaces only; PISA counts salt bridges separately from H-bonds): 5hpt D-E 4 H-bonds (PISA 4), E-D 4 (4), H-G 8 (8), F-D 6 H-bonds and 5 salt bridges (8, 7), C-A 4 and 5 (7, 5); 7d47 A-B 0 H-bonds and 1 salt bridge (0, 1); 9lqk A-B 6 H-bonds and 1 salt bridge (3, 3); A0002 8 and 5 (9, 7); A0001 12 and 14 (13, 15). Counting every N or O pair within 3.5 A between the chains gives 4 (5hpt D-E, all counted), 2 (7d47, all counted), 13 of 15 (protein complex A0002) and 12 of 15 (9lqk, the rest are O to O or N to N pairs, called polarity mismatches). Differences of one to three bonds are what is expected between two sets of rules; ifacekit does not use symmetry mates and PISA does.

How often the check matters (screening of 160 random structures from the 13,115 entries of section 3.10, with 2 to 8 chains and files up to 400 kB; every interface of the whole file, chains and ligands):

| | interfaces | with a bridging water | bridges | doubtful | interfaces with a doubtful bridge | interfaces without any clear bridge |
|---|---|---|---|---|---|---|
| chain - chain | 351 | 329 | 2,703 | 797 (29%) | 268 (81% of those with bridges) | 13 |
| ligand - protein or other | 1,334 | 562 | 1,021 | 365 (36%) | 267 (48%) | 132 |

The flags that occur (a bridge can carry several): chain - chain: long_bond 462, narrow_angle 322, triangle 238, too_short 167; ligand: narrow_angle 189, long_bond 160, triangle 151, too_short 103. The share of doubtful bridges depends strongly on the limits, which were chosen to catch the obvious cases (the 9lqk waters above) and were not calibrated against density or expert judgement: with the longest allowed bond set to 3.5 A instead of 3.3 A the share is 20% (758 of 3,744 bridges), and with also 3.2 A for the partners and 70 degrees for the angle it is 14% (532). The check is therefore a screening that points at waters worth looking at, not a validated test. The counts of direct bonds and of all bridges were not changed by it.

## 5.11 Comparison with MolProbity on 293 structures

Until this test the contact thresholds had been set and checked on five structures (section 5.3). The comparison was repeated on a larger and independent set. Method: 300 random entries of the 13,115 entries of section 3.10 (X-ray, resolution 1.8 A or better, R-free 0.25 or lower, one per 30% sequence-identity cluster; files up to 450 kB, 1 to 8 chains with more than 20 non-water residues; the structures used earlier excluded), 150 released from 2024-01-01 and 150 before. phenix.molprobity (Phenix 2.0) ran on 293 of them; on seven (9sqk, 9gir, 5bmt, 5nak, 5lom, 4l9b, 5ko4) Phenix stopped with an atom selection syntax error on TLS selections written in braces in the file header (for example {A|* E|*}; other files with TLS records ran normally), which is not related to ifacekit. ifacekit analysed all 293. The script is tools/compare_molprobity_batch.py and the result file of the run is kept with the other release files. The contact thresholds were not changed or fitted in this test. The test was run twice: first with version 0.6.2 and then with 0.7.1, which scores alternate conformers (sections 3.11 and 3.17) and has the calibrated rotamer cut-offs (section 5.12). Both are given.

**Clashes.** MolProbity's severe clashes (overlap 0.4 A or more, after Reduce adds hydrogens) were compared with ifacekit by residue pair. MolProbity reported 7,159 severe clashes: 1,693 between heavy atoms and 5,466 involving a hydrogen, which ifacekit cannot see.

| Heavy-atom clashes of MolProbity, by kind | flagged by 0.6.2 | flagged by 0.7.1 |
|---|---|---|
| between two different residues, not neighbours | 305 of 313 (97.4%) | 305 of 313 (97.4%) |
| with a water | 678 of 830 (81.7%) | 796 of 830 (95.9%) |
| involving an alternate conformer | 87 of 210 (41.4%) | 202 of 210 (96.2%) |
| between neighbouring residues (i and i+1) | 57 of 111 (51.4%) | 72 of 111 (64.9%) |
| inside one residue | 54 of 229 (23.6%) | 73 of 229 (31.9%) |
| all | 1,181 of 1,693 (69.8%) | 1,448 of 1,693 (85.5%) |

With the alternate conformers scored, ifacekit flags 85.5% of MolProbity's heavy-atom clashes and 93.9% of those between two different residues (1,464 clashes; 77.0% before). What is left is mostly outside what ifacekit evaluates: it skips atom pairs up to five bonds apart (inside a residue and between neighbours; the rise from 23.6% and 51.4% comes from alternate-conformer clashes in these groups). The share found falls with the overlap (85.5% at 0.4 A or more, 83.7% at 0.6 A, 73.5% at 0.8 A) because the largest overlaps are dominated by the excluded kinds. By element pair, N-O, O-O and C-O clashes are found in 96%, 94% and 87%, C-C and C-N in 62% and 46%. Of the 5,466 clashes that involve a hydrogen, the residue pair is flagged in 35.4%, and one of the two residues is flagged somewhere in 59.0%.

The other way round, 79.1% of ifacekit's 628 unfavorable residue pairs (497) are also in MolProbity's severe list (72.1% in its heavy-atom list; 79.0% of the 404 steric-clash pairs; 78.1% of 502 pairs before). The 12,563 borderline pairs agree in only 14.6%: borderline is a "look at it" category (overlap from 0.25 A, other radii) and not a clash list. Results were similar in the two halves of the set (heavy-atom recall 86.9% for entries released from 2024 and 84.6% for earlier ones; unfavorable pairs confirmed 82.3% and 77.5%).

Per structure, the number of ifacekit unfavorable pairs per 1000 atoms follows MolProbity's clashscore only moderately (Spearman 0.56; with borderline pairs added 0.53) and the count of unfavorable pairs follows MolProbity's heavy-atom clash count better (0.76). This is expected because three quarters of MolProbity's clashes involve hydrogens. In 112 of the 293 structures ifacekit has no unfavorable pair although MolProbity reports severe clashes in all but 5 (median clashscore 3.8; median ifacekit unfavorable pairs per 1000 atoms 0.42). So ifacekit does not replace a clashscore. It finds most of the heavy-atom clashes that MolProbity finds, and flags hardly anything that MolProbity does not.

**Ramachandran** (not changed between the two runs). The 148 entries released from 2024 were scored with the reference built from entries released before 2024 only, so this part is independent of the reference (51,036 residues). MolProbity (ramalyze) reports 37 outliers, ifacekit 29, 27 are shared: 73.0% of MolProbity's outliers are found and 93.1% of ifacekit's are confirmed. The favored percentage per structure is 98.1% (MolProbity) and 98.6% (ifacekit), a mean absolute difference of 0.45 points (Spearman 0.87), so ifacekit's favored regions are somewhat more generous, as noted in section 5.7. With the packaged reference, which contains these entries, ifacekit reports 22 outliers (21 shared), so a reference that contains the entries is less sensitive for them. The 145 earlier entries (in the reference set): MolProbity 60 outliers, ifacekit 61, 42 shared. Alternate conformers of Ramachandran are not scored separately. The numbers are small (37 outliers in 148 structures).

**Rotamers.** Same design (independent subset, reference built from entries before 2024). MolProbity (rotalyze) reports 369 outliers on the 148 entries. With version 0.6.2 ifacekit reported 340, 220 shared: 59.6% found, 64.7% confirmed (per structure, Spearman 0.59). With 0.7.1 (calibrated cut-offs and alternate conformers): ifacekit reports 344, 290 shared: 78.6% of MolProbity's outliers found, 84.3% of ifacekit's confirmed, Spearman 0.85 (earlier entries: 626 outliers, ifacekit 526, 491 shared, 78.4% found, 93.3% confirmed, Spearman 0.93). Of MolProbity's 369 outliers, 99 are in residues with alternate conformers; with 0.6.2 only 21% of these were found, against 78% of the 270 others, which is why scoring every conformer matters. By type (shared / only ifacekit / only MolProbity): Leu 47 / 3 / 15, Lys 33 / 6 / 13, Val 31 / 6 / 6, Arg 24 / 1 / 9, Ser 28 / 3 / 2, Asn 11 / 17 / 3, Ile 12 / 1 / 6. Leu and Ile, which had been over-called (101 of 120 outliers that only ifacekit reported), are now close; Asn is still over-called (17 outliers only ifacekit reports, 39% confirmed) and is the weakest type.

**Fact sheet.** The results of this section and the list of the 300 cases (PDB ID, release date, resolution, size, and the results per structure; the seven cases where MolProbity failed are marked; the 900 calibration structures are listed at the end) are collected on one page, `ifacekit/data/factsheet.html`. The notes at the end of the reports and the footer of the portal home page link to it: the portal serves it at `/factsheet` and passes this address to the reports (`--factsheet-url`); without an address the note says how to write the page to a file (`ifacekit --factsheet FILE`). The page has a filter and a CSV download of the case table. It is made by `tools/make_factsheet.py` from the result files, so no number is typed by hand, and it must be made again whenever a new comparison is run; it states the version of ifacekit that produced the numbers.

**What the test does not cover.** All structures are high-resolution X-ray entries of the kind used to build the references. Differences are likely to be larger for low-resolution models, cryo-EM models, AI predictions and nucleic acids.

## 5.12 Calibration of the rotamer cut-offs

Setting. rotalyze output (phenix.rotalyze, every residue with its score, chi angles and category) was collected for 900 calibration entries (600 released before 2024, 300 from 2024; none of the 293 comparison cases) and for the comparison cases. The chi angles of ifacekit agree with rotalyze: of 511,537 angles compared on 245,200 residues, 22 differ by more than 1.5 degrees (all are 0 against 180 degrees edge cases of symmetric groups), and for the alternate conformers 6 of the angles of 4,615 conformers differ. So the definitions are not the cause of disagreements. Residues with alternate conformers were left out of the calibration.

Method (tools/rotamer_reference/05_calibrate_cutoffs.py). For calibration the reference grids were built from entries released before 2024 without the calibration entries, so that no residue is scored with a grid that contains it. For each type and each smoothing factor (1.0, 1.5, 2.0, 3.0) the share of reference residues enclosed by the allowed cut-off was searched for the best F1 of the outlier class against rotalyze, and then the favored share for the best three-class agreement. Older and recent calibration entries have the same weight (the score of a setting is the mean of the two), because older entries have more rotamer outliers (in a first run on older entries only, Asn had 136 rotalyze outliers in 600 entries against 11 in 150 recent test entries, and a calibration on them over-called Asn on recent structures). Types with fewer than 40 rotalyze outliers keep the earlier shares. The factor 1.5 gave the best sum of outliers times F1 (1,739 / 1,774 / 1,637 / 1,207 for 1.0 / 1.5 / 2.0 / 3.0).

| Type | rotalyze outliers (calibration) | F1 before | F1 after | favored / allowed share |
|---|---|---|---|---|
| Arg | 228 | 0.69 | 0.79 | 0.85 / 0.97 |
| Asn | 167 | 0.61 | 0.74 | 0.95 / 0.985 |
| Asp | 123 | 0.69 | 0.86 | 0.97 / 0.992 |
| Gln | 142 | 0.81 | 0.84 | 0.90 / 0.985 |
| Glu | 245 | 0.82 | 0.80 | 0.90 / 0.985 |
| Ile | 125 | 0.55 | 0.87 | 0.95 / 0.992 |
| Leu | 321 | 0.54 | 0.80 | 0.95 / 0.99 |
| Lys | 348 | 0.76 | 0.76 | 0.85 / 0.975 |
| Met | 69 | 0.82 | 0.82 | 0.90 / 0.98 |
| Phe/Tyr | 66 | 0.75 | 0.82 | 0.95 / 0.992 |
| Ser | 115 | 0.83 | 0.93 | 0.95 / 0.994 |
| Thr | 79 | 0.54 | 0.77 | 0.95 / 0.99 |
| Val | 119 | 0.80 | 0.87 | 0.95 / 0.985 |
| Cys, His, Pro, Trp | 19, 36, 25, 13 | not changed | | 0.981 / 0.997, 0.970 / 0.995, 0.995 / 0.9997, 0.965 / 0.995 |

Tests. Test A: 150 entries released from 2024 (40,595 residues without alternate conformers, 269 rotalyze outliers), scored with grids from entries before 2024. Test B: 150 older entries (42,395 residues, 477 outliers) that were left out of the grids scoring them. Neither was used for the calibration.

| | rotalyze outliers | ifacekit outliers before | after | shared before | after | found before | after | confirmed before | after | three-class agreement before | after |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Test A | 269 | 309 | 256 | 198 | 208 | 73.6% | 77.3% | 64.1% | 81.2% | 97.0% | 98.4% |
| Test B | 477 | 468 | 412 | 357 | 378 | 74.8% | 79.2% | 76.3% | 91.7% | 97.4% | 98.6% |

A calibration that gave every residue the same weight (older entries dominate) gave 78.8% found and 76.5% confirmed on test A, and 81.3% and 89.8% on test B; the equal weight of the two groups was chosen because it keeps more of the precision on recent entries (81.2% against 76.5%). Asn stays the weak type on recent entries (test A: 23 outliers by ifacekit, 10 shared, rotalyze 11). The calibration used rotalyze output on the raw files, without Reduce flips. The cut-offs were therefore fitted to MolProbity's calls, and the packaged reference states that.

## 5.13 Clashes with hydrogens added, against MolProbity

**Placement against Reduce.** Reduce (used only here, as a reference) was run on three structures (9cdy, 8zmo, 8y7v). Of 4,562 hydrogens that both programs place on the same atom, 96.9% lie within 0.15 A of the Reduce position (median 0.02 A); the rest are turned groups (hydroxyl, thiol, Tyr OH) and flipped groups. Asn, Gln and His differ in the choices of flip and protonation, which are compared below.

**Flips.** On the first 40 cases Reduce reports its choice for every amide and ring (USER MOD lines). With the gain a flip must show set to 0.05 units, ifacekit flipped 15% of the amides, Reduce about 7%. The gain was raised until the two agreed best: 2 units gave 1,029 of 1,105 agreeing amides, 4 units 1,060 (95.9%; 6 flips only ifacekit makes, 39 only Reduce makes) and 6 or 10 units 1,061 or 1,062, so 4 was kept. His rings agree in 305 of 313.

**Clash lists.** tools/compare_hydrogen_clashes.py compares, for every case, the clash pairs of ifacekit with hydrogens added with the pairs of MolProbity's probe list (residue-pair level, clashes of 0.4 A or more). The limits (overlap allowed for hydrogen bonds, the bond count that excludes close hydrogens, the cap for N/O pairs, the gain of a flip) were fitted on the first 40 of the 293 cases; the other 253 were not used for any setting. Reduce's own hydrogens, scored with the same rules, reach an F1 of 0.79 against the list, which is the ceiling of this scoring. On the first 20 cases (first conformer only) Reduce's hydrogens give 84.3% found and 74.2% confirmed. A rule closer to Probe (dots on each atom surface, each dot given to the neighbour it overlaps most, dots inside bonded atoms dropped, hydrogens five bonds apart compared) was tried with Reduce's hydrogens and reached 85.7% and 75.6% on the same cases: too little to adopt. Agreement of 90% or more would need Reduce and Probe to be copied in detail, including the hydrogens they put on ligands.

| Structures | cases | MolProbity clash pairs | of which with a hydrogen | ifacekit pairs | shared | found | found, with a hydrogen | found, heavy atoms only | confirmed | correlation of the counts |
|---|---|---|---|---|---|---|---|---|---|---|
| all | 293 | 4,781 | 3,891 | 4,914 | 3,779 | 79.0% | 82.8% | 62.8% | 76.9% | 0.93 |
| tuned on (first 40) | 40 | 559 | 453 | 578 | 446 | 79.8% | 84.1% | 61.3% | 77.2% | 0.96 |
| independent | 253 | 4,222 | 3,438 | 4,336 | 3,333 | 78.9% | 82.6% | 63.0% | 76.9% | 0.93 |

Alternate conformers are left out of the comparison: the hydrogen check works on the first conformer of every atom, and the 735 further clash pairs of MolProbity that involve an atom of another conformer (13% of its list) are not counted. (Counting them, the same ifacekit pairs give 68.7% found.) In plain words: of the 4,781 clash pairs MolProbity reports, ifacekit finds 3,779 (79%); of the 4,914 pairs ifacekit reports, 3,779 (77%) are also in MolProbity's list. The fact sheet shows this as two overlapping circles.

The agreement is the same on the cases that were not used for tuning, so the limits are not overfitted to the 40. The fact sheet shows the table and the numbers per case.

**Why the two lists differ by 20 to 30%.** The two programs agree on 3,779 pairs. 1,002 pairs (21% of MolProbity's 4,781) are only in MolProbity's list, and 1,135 (23% of ifacekit's 4,914) only in ifacekit's. tools/explain_hydrogen_disagreement.py gives every such pair the most likely cause by fixed rules (the causes are attributions, not proofs), over all 293 cases:

| Only MolProbity reports | pairs | share | Most likely cause |
|---|---|---|---|
| hydrogens of waters, ligands, other non-protein residues | 300 | 30% | MolProbity places hydrogens there; ifacekit places none |
| hydrogen-bond rule | 217 | 22% | ifacekit counts a polar hydrogen next to an acceptor as a hydrogen bond, MolProbity does not for these pairs |
| bond-count rule | 135 | 13% | hydrogens of one residue or neighbours whose heavy atoms are three bonds apart or fewer: ifacekit leaves them out, MolProbity counts some |
| overlap just under the limit | 111 | 11% | 0.30 to 0.40 A in ifacekit, 0.40 A or more in MolProbity |
| two N or O atoms, overlap under 0.55 A | 103 | 10% | ifacekit calls it a short hydrogen bond |
| hydrogen placed elsewhere | 74 | 7% | turned hydroxyl, thiol or NH3, flip, protonation |
| two heavy atoms (protein / water or ligand) | 37 / 25 | 4% / 2% | no overlap of 0.3 A in ifacekit |

| Only ifacekit reports | pairs | share | Most likely cause |
|---|---|---|---|
| overlap 0.40 to 0.45 A | 505 | 44% | MolProbity's value for the same atoms is just under 0.40 A |
| two heavy atoms, one a water | 307 | 27% | in MolProbity the water carries hydrogens, so the pair is a hydrogen bond there; ifacekit has no water hydrogens and exempts only overlaps under 0.55 A |
| hydrogen with fixed geometry (C-H, backbone N-H) | 249 | 22% | the neighbouring group is turned or flipped differently, or Probe gives the overlap to another atom |
| hydrogen of a turned or flipped group | 32 | 3% | state chosen differently from Reduce |
| two heavy atoms (ligand or ion / protein) | 19 / 23 | 2% / 2% | different radii; protein heavy atoms |

So the disagreement is mostly a matter of definitions and of the hydrogens ifacekit does not place, not of wrong geometry. Of the pairs only MolProbity reports, 30% are hydrogens on waters and ligands (a limit of ifacekit that could be lifted by placing hydrogens there), 45% are choices between two reasonable rules (hydrogen bond, bond count, short N/O pair), 11% are a hair under the limit, and only 7% come from a hydrogen placed in another position. Of the pairs only ifacekit reports, 44% are within 0.05 A of the limit: the two programs do not measure overlap in the same way (MolProbity's Probe measures it on dots on the atom surfaces and gives each dot to the neighbour it overlaps most; ifacekit uses the distance between the atom centres), and a pair that is 0.42 A in one is 0.38 A in the other. 27% are pairs with waters, to which MolProbity gives hydrogens. The remaining third has no single cause. A limit of 0.42 A instead of 0.40 A would remove many of the first group but lose more of MolProbity's pairs (F1 0.72 against 0.74 for the pairs with a hydrogen, first 40 cases), so 0.40 A was kept. With Reduce's own hydrogens the same rules still disagree with MolProbity on about a quarter of the pairs (see below), which shows how much of the difference belongs to the rules and not to the hydrogens.

MolProbity also puts hydrogens on ligands (it has a dictionary for them), which ifacekit does not. Clashes between two heavy atoms are found less well here (62.8%) than by the heavy-atom contact tables (section 5.11, 85.5%, which also scan alternate conformers), so for heavy-atom problems the contact tables remain the better guide. The counts per structure follow MolProbity's closely (correlation 0.93) but the number is not MolProbity's clashscore. Time: a median of 0.4 s per structure for placement and the clash check (2.6 s for the largest of the cases, 11,800 atoms; about 0.2 ms per atom, mostly Python loops). Not tested: low-resolution and cryo-EM models, predictions, structures with nucleic acids.

**ChimeraX as a hydrogen source (one-off check, not part of the package).** On all 293 cases (first conformer, no alternate conformers) ChimeraX 1.12 added hydrogens (`addh hbond true`, which also protonates waters and ligands: 1.01 million hydrogens against 0.80 million placed by ifacekit) and the clashes were scored in three ways against the same MolProbity pairs. (1) With ifacekit's clash rules on the ChimeraX positions as they are (C-H 1.09 A, the nuclear distances): 79.8% of MolProbity's pairs found, but only 16.7% of the 22,819 pairs confirmed, because the radii of the rules belong to the shortened X-ray positions. (2) The same after moving every hydrogen to the shortened X-ray distance: 74.4% found, 43.7% of 8,144 pairs confirmed. (3) With the `clashes` command of ChimeraX itself (overlap 0.4 A, hydrogen-bond allowance 0.6 A, bond separation 4): 42.8% found, 17.6% of 11,631 pairs confirmed. ifacekit's own hydrogens give 79.0% found and 76.9% of 4,914 pairs confirmed. ChimeraX hydrogens are therefore no better for this purpose: they do not flip Asn, Gln and His, put hydrogens on waters and ligands that MolProbity treats differently, and its clash rule is another one. This is why the hydrogens are placed by ifacekit. The data are kept (the cleaned models, the ChimeraX models with hydrogens and the ChimeraX clash lists, 266 MB, in the comparison folder next to the MolProbity results), the scripts are tools/compare_chimerax_hydrogens.py and tools/chimerax_addh_clashes.py, and the fact sheet shows the table, the overlap diagrams and a column per structure.

## 5.7 Ramachandran

Comparison with phenix.ramalyze (Phenix 2.0), residue by residue, on 8,815 residues of seven structures (4QZS, 5hpt, 7yqm, 8afb, 9lqk, protein complexes A0003 and A0002). Residue classes agreed for every residue. The category agreed for 8,666 residues (98.3%). All differences lie at the boundary between favored and allowed: 135 residues that ramalyze calls allowed are favored here and 14 go the other way. The three outliers (all in 5hpt) are the same in both programs. ifacekit also reports 11 residues next to alternate conformations that ramalyze leaves out (first conformer used).

The first version of the module used the Top8000 tables from the MolProbity data folder with MolProbity's fixed cut-offs and gave 100% agreement on the same residues. That version was replaced because the redistribution terms of those tables were not clear.

Sensitivity of the own reference, on 8,826 residues of the same structures (agreement refers to the Top8000 categories; Top8000 itself gives 96.67% favored and 0.03% outliers):

| R-free limit | sigma | favored | outliers | agreement |
|---|---|---|---|---|
| 0.25 | 3 | 97.68% | 0.06% | 98.58% |
| 0.25 | 4 | 98.03% | 0.03% | 98.32% |
| 0.25 | 5 | 98.27% | 0.03% | 98.13% |
| 0.23 | 3 | 97.64% | 0.07% | 98.61% |
| 0.23 | 4 | 98.05% | 0.03% | 98.30% |
| 0.23 | 5 | 98.29% | 0.03% | 98.11% |
| 0.20 | 3 | 97.81% | 0.11% | 98.39% |
| 0.20 | 4 | 98.24% | 0.05% | 98.37% |
| 0.20 | 5 | 98.43% | 0.03% | 98.22% |

The R-free limit hardly changes the regions (most listed entries are already well refined), and a limit of 0.20 removes 40% of the residues and leaves only 6,080 cis-proline residues. R-free 0.23 and sigma 4 degrees were chosen. The own regions are somewhat larger than Top8000: for general residues 5,742 favored grid cells against 5,273. The reason is not known. Candidates are the homology filter (70% in Top8000, 30% clusters here), a clash filter in Top8000 and the smoothing method.

Held-out test (tools/rama_reference/04_heldout_test.py). Entries released from 2024-01-01 (1,634 entries, 450,428 residues) were held out and the grids rebuilt from the 11,481 earlier entries with the same method. The held-out residues fall into the categories at the same rates as the reference residues, so the regions generalise to newer structures:

| Set | residues | favored | allowed | outlier |
|---|---|---|---|---|
| reference (before 2024) | 3,293,230 | 98.45% | 1.52% | 0.029% |
| held out (2024 and later) | 450,428 | 98.49% | 1.48% | 0.033% |
| Top8000 tables on the held-out residues | 450,428 | 98.09% | 1.86% | 0.048% |

Per class the held-out favored rates are 98.3% to 99.0%. The cis-proline class has only 1,173 held-out residues (three outliers). The held-out residues passed the same quality filters as the reference (B below 30, occupancy 1, R-free 0.23 or lower), so the test shows that the regions generalise to new well-refined structures, not how they behave on poorly refined models.

## 5.8 Rotamers

Comparison with phenix.rotalyze (Phenix 2.0) on 7,546 residues of the same seven structures, residues with a single conformer. The chi angles agree for every residue (after folding the symmetric end groups, which rotalyze prints unfolded). The category agrees for 7,110 residues (94.2%). Of the 154 rotalyze outliers, 119 are outliers here, 32 are allowed and 3 favored; ifacekit calls 56 more residues outliers that rotalyze calls allowed. Agreement by type is 88.8% to 99.3%, lowest for Leu (88.8%), then Gln, Lys, Ile and Val (about 90% to 91%), and 98.7% or more for Asn, Phe/Tyr, Pro, Trp and Met. The disagreements are mostly between favored and allowed.

The comparison was also used to choose the smoothing: fixed peak-relative cut-offs gave 83.9%, enclosure cut-offs 91.9%, and widening the smoothing by 1.5, 2 and 3 gave 94.2%, 92.0% and 89.9%. The factor 1.5 is therefore one parameter set on the same seven structures that are used for the comparison, and the agreement figure is not an independent test.

Note (0.7.1): the comparison with rotalyze in the paragraphs of this section was made on seven structures with the first cut-offs. The cut-offs were later calibrated on 900 entries and tested on 300 others, see sections 5.11 and 5.12 for the current numbers.

Held-out test (tools/rotamer_reference/04_heldout_test.py), same split as for Ramachandran (entries from 2024-01-01 held out, grids rebuilt from the earlier entries): 312,154 held-out residues are 97.60% favored, 2.11% allowed and 0.290% outliers, against 97.44%, 2.15% and 0.409% for the reference residues. By type the held-out outlier rates are 0.05% to 1.1%. Leu stands out with 1.1% outliers (1.5% in the reference residues) because the MolProbity Leu table removes some rotamers (tt and mp) and the cut-offs follow it. The same limits as for Ramachandran apply: well-refined residues only.

## 5.9 Flips

The flip check was compared with Reduce 3.24 (the version distributed with MolProbity, run with -BUILD, which optimises flips), residue by residue, on 70 random entries released before 2024 (training set: 2,380 groups, 165 flipped by Reduce, 6.9%) and 70 random entries released from 2024 (test set: 2,948 groups, 91 flipped by Reduce, 3.1%), all from the reference entries of section 3.10. The two thresholds of section 3.12 were set on the training set and then fixed. A group counts as a Reduce flip when Reduce prints FLIP for it.

| Score gain needed | Set | suggested right | suggested wrong | missed | precision | recall |
|---|---|---|---|---|---|---|
| 2 (flip) | training | 70 | 31 | 95 | 0.69 | 0.42 |
| 2 (flip) | test | 20 | 49 | 71 | 0.29 | 0.22 |
| 1 (flip or uncertain) | training | 92 | 97 | 73 | 0.49 | 0.56 |
| 1 (flip or uncertain) | test | 36 | 122 | 55 | 0.23 | 0.40 |

An earlier version that used only distances (no hydrogen geometry) gave precision 0.18 and recall 0.15 on the test set at a gain of 2. The check is useful to point at groups worth a look (about 9 times more often a Reduce flip than average at a gain of 2 on the test set) but it is a weak predictor of Reduce's decision, and much weaker on the newer entries, in which most groups had probably already been placed by a Reduce-style step during refinement, so that the remaining flips are the hard cases (this explanation was not tested). It must not be used to flip groups without looking at the density. A better check would need hydrogens and the hydrogen-bond network optimisation of Reduce, which ifacekit does not do by design.

# 6. Known limitations

* The contact thresholds were set on five structures and compared with MolProbity on 293 independent high-resolution X-ray entries (section 5.11): 85.5% of MolProbity's heavy-atom clashes (93.9% of those between different residues) are flagged, but the heavy-atom tables see no hydrogens (three quarters of MolProbity's clashes) and skip atom pairs up to five bonds apart; the optional hydrogen check (section 3.18, 5.13) finds 79% of MolProbity's clash pairs (first conformer only), and ifacekit is not a clashscore. Low-resolution models, cryo-EM models, predictions and nucleic acids were not part of the comparison.
* The rotamer cut-offs were calibrated against MolProbity's rotalyze on 900 entries and tested on 300 others (section 5.12): about 79% of MolProbity's rotamer outliers are found and 84 to 92% of ifacekit's are confirmed. Asn is still over-called on recent entries, Cys, His, Pro and Trp keep uncalibrated cut-offs (too few outliers), and alternate conformers are scored for rotamers and contacts but not for Ramachandran, interface areas, bonds and the flip check.
* Hydrogens are used only in the optional clash check (section 3.18); they are placed on proteins only (no waters, ligands, nucleic acids or modified residues), so clashes of those with hydrogens are not seen, and placement is validated on X-ray models only. The Asn, Gln and His flip check of the rotamer page is still a heavy-atom test that agrees weakly with Reduce (section 5.9); it only suggests flips. The flips chosen by the hydrogen check agree with Reduce's on 96% of the amides but are not used by the rest of the report.
* The Ramachandran reference was checked against phenix.ramalyze on seven structures, and the test structures were not checked for membership in the reference set. A held-out test on 2024 and later entries was done (section 5.7) but only for well-refined residues. The favored regions are somewhat larger than in MolProbity, so favored percentages can come out higher. Residues are classified with the first conformer.
* The ipSAE code is a re-implementation of the published definition and has not been compared against the reference code. The pDockQ constants have not been verified.
* The colour cut-offs for AI metrics are common rules of thumb.
* A report covers one chain pair. There is no overview of all pairs and no table to compare many designs.
* The "whole model" clash count includes all chains of the model, also those that are not in the chosen pair.
* PISA's free energy of dissociation and assembly prediction are not implemented, by design.
* The 3D page links from tables to Mol* but not back from Mol* to the tables. The residue labels overlap when the panel is small.
* Without a dictionary in the file the ligand bond orders are guessed from geometry and can be wrong for unusual chemistry. Four ligand cases were checked so far. Modified residues are not mapped to their parent amino acid, so they appear as an unclassified residue.
* In the modern ligand view a few residue cards can still overlap each other or sit close to a ligand atom, and a long label can touch a hydrogen bond label.
* The model type detection assumes an experimental structure when the file shows no sign of either. A prediction from a program with no keyword in the file and no fractional confidence needs `--predicted`.
* The browser downloads of SVG, PNG and PDF were tested in headless Chromium only, not yet in Safari.
* The check of the bridging waters uses four simple limits that were not calibrated (section 5.10); 14 to 31% of the bridges are flagged depending on the limits, and a flag is only a reason to look at the density.
* Bonds via water need waters in the model; water positions in X-ray models are less certain than protein atoms, and a bridge through two waters is not counted. The rule for the simple ligand map was judged by eye on two ligands only.
* The plots PDF was only checked here by turning its pages into images, not in other PDF viewers, and the plots on a page are small (six per page). The shaded backgrounds are smoothed when scaled.
* The portal has no log-in. It is meant for a trusted internal network (a private address range of the research group), which is the decision for the lab server; add nginx authentication if that changes. The limit on simultaneous analyses is counted per server process, so the launcher runs one process with many threads.

# 7. Web portal

The portal (server/app.py, home page server/index.html) is for internal users. It offers four kinds of analysis:

1. AI-predicted complex: a structure plus the optional AlphaFold summary_confidences and full_data JSON files. Chain pair chosen by the user or, when left empty, the pair with the largest interface.
2. Experimental interface: one PDB or mmCIF file, two chains.
3. Ligand binding site: one ligand chosen from the list read from the file. The model type is detected from the file. Report, 3D page and interaction maps (SVG, PNG, PDF).
4. Chain quality: chosen chains or all. The main result is an HTML report like the pair and ligand reports (overview and comparison table, B-factor or pLDDT plot, then one section per chain with contact problems, backbone geometry, Ramachandran plots, waters and metal sites). The text report, the Ramachandran page, the rotamer page (plots, outliers, flip check), the PDF with all Ramachandran and rotamer plots and the JSON are offered as extra downloads. The two interface types also give the PDF.

The work goes in two steps. The file is uploaded or fetched first and read in a separate process, and the page then shows the chains and ligands that are really in the file. The analysis is a second request that names the job.

Archives are not accepted, by design: an unpacked archive can carry anything. Every upload is checked for size (100 MB for a structure and the full_data file, 5 MB for the summary), for plain text without binary content (zip, gzip, executables and files with null bytes are refused, also under a changed extension), for the expected content (atom records for a structure, a parsable JSON object with the expected AlphaFold fields) and is stored under its own name reduced to letters, digits, dot, dash and underscore (no path parts, no leading dot, at most 60 characters, extension set from the checked content; an entry fetched by PDB ID is stored as for example 8AFB.cif). The reports, the PDF and the text therefore show what was analysed. The AlphaFold JSON files keep fixed names. The name is only used as a file name inside the job folder and is escaped wherever it is shown. Files are only read, never executed. The job folders are deleted after 7 days when the server starts.

A structure can also be fetched by PDB ID from RCSB, PDBe or PDBj, or from any of them in turn. Only the validated four-character ID is placed into one of three fixed HTTPS addresses, redirects to other hosts are refused, the size is limited and the PDBj answer (gzip) is unpacked with a size limit. The server needs outgoing HTTPS to files.rcsb.org, www.ebi.ac.uk and pdbj.org.

Server operation (0.5.0). `server/ifacekit-server.sh` starts the portal with gunicorn (one worker process with 16 threads) when gunicorn is installed, and with Flask's development server otherwise. Settings come from environment variables: `IFACEKIT_DATA` (job folder), `IFACEKIT_KEEP_DAYS` (7) and `IFACEKIT_PURGE_MINUTES` (60), `IFACEKIT_MAX_JOBS` (4) and `IFACEKIT_QUEUE_WAIT` (600 s), `IFACEKIT_MEM_MB` (8192), `IFACEKIT_THREADS` (16).

* Clean-up: a background thread of the server process deletes job folders older than the retention time every hour, so the 7 days named in the footer hold on a server that runs for months (before 0.5.0 the clean-up ran only at start-up). Only folders named like a job are touched.
* Simultaneous analyses: at most `IFACEKIT_MAX_JOBS` analyses run at once; further requests wait up to the queue time and then get the message that the server is busy.
* Memory: every analysis process runs under `ulimit -v` with `IFACEKIT_MEM_MB` (Linux); a structure that does not fit gets the message that it is too large for the memory limit. Time limit as before (10 minutes).
* Mol*: the 3D pages use `server/static/molstar` (molstar.js and molstar.css of Mol* 5.0.0, MIT licence file included) through the option `--molstar-url`, so they work without internet; a report made from the command line still uses the CDN unless the option is given.
* Footer: the home page and every report show the version, the MIT licence and the copyright holder, and the citation suggestion.
* `/healthz` returns the version, the analyses running and waiting, the number of stored jobs, the free disk space and whether local Mol* is present.
* `server/deploy/` has a systemd unit, an nginx example (TLS, upload size, long timeouts) and a README for installing it as a service. The portal needs the root of its own host name because its pages use absolute addresses. It still needs outgoing HTTPS to files.rcsb.org, www.ebi.ac.uk and pdbj.org for the PDB download.

Tested here: a gunicorn run with one analysis slot (one analysis running and one waiting at the same moment, both finished), the memory limit (a limit of 250 MB gave the message above), the clean-up thread (a job older than the retention time was deleted while the server was running), the local Mol* page (no address on the CDN in the 3D page, and molstar.js served), the home page and `/healthz`. Not tested: systemd, nginx, TLS, and a longer run under real load.

Tested with curl: the four analysis types, AF files, the PDB download from the three sources (entry 8AFB, and a missing entry 0000), and the refusal of archives, renamed archives, text files with a wrong extension and wrong JSON files. The home page was opened in the browser pane and showed no console errors, but the forms were not clicked through in a browser. The chain quality type was also run through the portal after a PDBj download, and its Ramachandran page was served.

Open items, in the order I would do them:

1. Readers for other predictors (Boltz, Chai, local AF3), because their file names and JSON fields differ.
2. Optional overview of all chain pairs and a table to compare many designs.
3. Comparison of MolProbity and ifacekit on the 1000 computed de novo proteins. This would be the first test that is independent of the calibration set. The same set can serve as the held-out test of the Ramachandran reference.
4. Rotamer names, and a test of the rotamer and Ramachandran calls on poorly refined structures, for example the 1000 computed de novo proteins. Ramachandran and rotamer sections in the ligand report, if wanted.
5. Ligand map follow-ups. Optional SMILES or chemical component input for ligands whose file has no dictionary, with a check of atom and element counts and a warning that a wrong SMILES gives a wrong drawing. A better card layout in the modern view. Compare the maps with LigPlot+ and PLIP on the same files. Then the protein to protein, DNA and RNA maps in the plan.

# 8. Installing on the server

The step-by-step guide for Linux and macOS (environment, optional parts, libcairo, check, portal, update, common problems) is INSTALL.txt in the archive. The summary below is kept for the design record. In the archive the `data` and `results` folders are empty and `testdata` holds the PDB models used by the tests.

Requirements: Python 3.9 or newer, and the packages gemmi (0.6.5 or newer), numpy and scipy for the core analysis. The portal also needs Flask (and gunicorn for the production launcher). The ligand report (`--ligand`) also needs rdkit. The PNG and PDF figure files also need cairosvg and the system library libcairo; without libcairo, `pip install cairosvg` succeeds but ifacekit writes SVG only and prints one warning. The 3D page needs a browser with internet access (Mol* from jsDelivr), unless it is made with `--molstar-url` or by the portal, which serves its own copy of Mol*. Phenix, CCP4, PISA and MolProbity are not needed. They were used only on the development machines to produce the files that the comparison scripts read. For the tests, pytest, freesasa, rdkit and cairosvg are also needed. Without freesasa two tests skip.

```
tar xzf ifacekit.tgz
cd ifacekit
python3 -m venv .venv
. .venv/bin/activate
pip install .
pip install ".[ligand]"            # rdkit + cairosvg, only for the ligand report (cairosvg only for PNG and PDF files)
# cairosvg also needs the system library libcairo (pip does not install it):
#   macOS: brew install cairo | Debian/Ubuntu: sudo apt install libcairo2 | Fedora/RHEL: sudo dnf install cairo
pip install ".[server]"            # flask + gunicorn, only for the web portal (see server/deploy/README_DEPLOY.txt)
pip install ".[pdf]"               # cairosvg + pypdf, only for --plots-pdf (also needs libcairo, as above)
pip install pytest freesasa        # only for the tests
python -m pytest -q tests
curl -O https://files.rcsb.org/download/5HPT.cif
ifacekit 5HPT.cif --pair A B --experimental --html /tmp/r.html --molstar
python server/app.py --port 8000    # the portal, then open http://127.0.0.1:8000
```

On a lab server the program folder can be, for example, /your_path/ifacekit (replace /your_path by your own folder), with a small launcher script for the command line and one for the portal (`ifacekit-server --port 8000`) on the PATH of the users. The portal listens on 127.0.0.1 only unless `--host 0.0.0.0` is given. Job files go to the folder given by IFACEKIT_DATA (for example /your_path/jobs).

Mol* is loaded from jsDelivr (version 5.0.0) in the current 3D pages, so the person opening the page needs internet until item 1 above is done.


# 9. Design notes to keep for the paper

## 9.1 Water in the interaction map

ifacekit treats water as part of the binding site. Waters that contact the ligand are drawn in the interaction map as small numbered circles, and a water that makes a hydrogen bond to the ligand is joined by a labelled line like any residue. They are also listed in the interaction tables with the same distance and verdict rules as protein atoms. In the 4QZS (BRD4 with JQ1) test, two waters make hydrogen bonds to the ligand (2.83 and 2.92 A), and the other contacting waters are shown without a bond.

Reason: ordered water is often part of the real recognition pattern, and a map without it can look as if a polar group of the ligand has no partner. The author's observation is that PLIP leaves water out of this picture. Before this statement goes into a paper, it must be checked against the PLIP version used, because PLIP does have a water-bridge interaction type in some versions. The safe wording is to say what ifacekit does (waters are drawn and counted), and to compare with PLIP and LigPlot+ by running them on the same files.

Open point: by default only waters with a hydrogen bond to the ligand should be drawn, and the others switched on by an option.

## 9.2 Ligand chemistry

When the input file carries the chemical component dictionary (files from the wwPDB do), ifacekit uses the bonds and aromatic flags written in the file. The 4QZS and 8JUX files give the correct ligand structures this way. Only when the dictionary is missing are bond orders guessed from the coordinates, and the report states which source was used ("file dictionary" or "geometry guess"). A user-supplied SMILES is a planned optional input and must be checked against the atom and element counts of the model.

## 9.3 Ramachandran reference

The reference is built from the PDB by a documented, repeatable procedure (section 3.10), so it can be described in the methods without relying on the licence of another program's tables. Statements for the paper: the criteria for entries and residues, the class definitions (as in MolProbity), the density method and the 98% and 99.95% enclosure rule, and the comparison with phenix.ramalyze (section 5.7) with its limits. Do not state that the results equal MolProbity: they agree on 98.3% of residues in the small test, with all differences at the favored and allowed boundary. The held-out test on entries from 2024 and later (section 5.7) supports that the regions generalise.

# 10. Licence and citation

ifacekit is released under the MIT licence (file LICENSE), the same type as ProteinMPNN and RFdiffusion. The version, the licence, the copyright holder and the citation suggestion are shown in the footer of the portal home page and of every HTML report. The Mol* files that the portal serves are under the MIT licence of the Mol* contributors, and their licence file is kept next to them. The copyright holder named in the file is Kuen-Phon Wu, Institute of Biological Chemistry, Academia Sinica; check the institution's policy on software ownership before the code is shared. Citation is a suggestion, not a condition. CITATION.cff and the portal footer carry a placeholder for the bioRxiv preprint, to be completed when it is posted. The PDB data used for the Ramachandran reference are CC0. AlphaFold Server output has its own terms of use (non-commercial), which apply to the users of the portal.

## 9.4 Rotamer reference

Statements for the paper: the entry and residue criteria (as for Ramachandran, with all chi atoms required to be well ordered), the 17 tables and chi definitions (as in MolProbity), the density method, and how the cut-offs were set: per type, the enclosure that MolProbity's peak-relative cut-offs give on its Top8000 tables, applied to ifacekit's own grids. This must be said openly, because two numbers per type come from MolProbity's data. Report the comparison with rotalyze (94.2%, section 5.8) with its limits: one smoothing factor was set on the same seven structures, Asn, Gln and His are not flipped, and ifacekit calls more outliers than rotalyze. The held-out test (97.6% favored, 0.29% outliers on 2024 and later entries) shows that the reference generalises, not that it is correct for poorly refined structures.
