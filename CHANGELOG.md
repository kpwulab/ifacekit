# Changelog

All notable changes to ifacekit. Versions follow `MAJOR.MINOR.PATCH`; the version of the JSON result layout (`schema_version`)
is noted separately because it changes less often. Newest first. Changes that are not released yet are listed under *Unreleased*.

## [Unreleased]

## [0.8.0] - 2026-11-18

Result layout: `schema_version` 0.12 (new contact field `hbond_geometry`, `hbonds.direct_polar_contacts`, the `ensemble` result mode, the flag `ensemble`).

### Added
- NMR ensembles and other multi-model files: option `--ensemble` analyses every model of the file (with `--chains`, `--pair` or `--all`) and gives the numbers of each model (Ramachandran, rotamers, bond and angle RMSZ,
  heavy-atom contact verdicts, clashes with hydrogens when `--hydrogens` is on, interface area and direct hydrogen bonds for interfaces), their mean, standard deviation, minimum and maximum, and the complete analysis of the first model.
  `--ensemble-csv FILE` writes one row per model plus the mean, sd, min and max; `--workers N` runs the models in parallel. Python: `ifacekit.ensemble.analyze_ensemble(path, chain=, pair=, models=, workers=)`;
  result mode `ensemble` with `models`, `average` and `first_model`. The HTML reports (`--html` with `--chains` or `--pair`) show the tables of the first model, and the B-factor / pLDDT plot is replaced by the deviation of each residue from the mean structure over the models: a line (mean over the models of the root-mean-square distance of the backbone atoms N, CA, C, O from their mean positions, each chain superposed on its least mobile 70% of residues) with a band of one standard deviation and a lighter band from the lowest to the highest model; the RMSD is in the tooltip (`ifacekit.ensemble.residue_rmsd`, `residue_rmsd` in the JSON, `FILE_residues.csv` next to `--ensemble-csv`). With `--ensemble` the Ramachandran page (`--rama-html`), the rotamer page (`--rotamer-html`) and the plots PDF (`--plots-pdf`) show every model of the file in one plot (hover for the model; the counts are residue-model pairs), and the summary of the ensemble page holds the same plots. The report of a single model keeps the plots of that model. `--models LIST` chooses the models that get a report (first, all, or the model numbers of the PDB file, e.g. `3`, `1,4,7`, `2-5`; implies `--ensemble`): one model gives the report of that model; several models with `--html` give one page in the frame of the normal report (banner, title, notes) with a bar of buttons: the ensemble summary (table of all models with mean and standard deviation, the band plot of the deviations, and the Ramachandran and rotamer plots with every model in one plot) and one button per model, which swaps the tables and sections below the bar for those of that model (`--separate` writes one file per model instead, `FILE_model03.html`). With `--molstar` the page is also written with the Mol* 3D panel as `FILE_3d.html`: the panel shows the model chosen on the page, a click on a table row shows those atoms, and the button "all models" overlays all models. Module `ifacekit.ensemblereport`. The portal offers it for chain quality and experimental interfaces when the uploaded file has more than one model (check box "Analyse every model of the file"; the 3D page and CSV tables are offered for download). Without `--ensemble` only model 0 (or `--model N`) is read, as before.

### Fixed
- A file with an experimental record and all B-factors 0 (NMR models) was taken for a prediction whose confidence is stored as a fraction (pLDDT 0). The fraction rule now needs a non-zero B-factor and no experimental record.

### Changed
- Hydrogen bonds are counted with a geometry test, still without hydrogens: X...A of at most 3.3 A, the donor points to the acceptor (the expected N-H direction within 50 degrees of X->A, or N-X...A of at least 100 degrees for a
  donor with one neighbour) and every angle at the acceptor at least 100 degrees (`ifacekit.contacts`, constants `HBGEO_*`). The new contact field `hbond_geometry` marks them; the direct count and the water bridges of
  `hbonds` use it, and the older distance-only count is kept as `hbonds.direct_polar_contacts`. The contact type `hbond` (polar contact up to 3.5 A) and its verdicts are unchanged, so the clash and verdict results do not change.
  On 279 interfaces the mean number of direct bonds falls from 10.2 to 8.3 per interface (PISA 7.1, ChimeraX 6.1 to 6.7, PyMOL assessment 7.3).
- Atoms on the same position (duplicated chains) are no covalent neighbours in the geometry test.
- Every HTML report and the portal show a small box with the ifacekit version under the banner.
- The tests use public PDB entries that are downloaded from RCSB on first use (cache `tests/.pdb_cache` or `IFACEKIT_TEST_CACHE`); the repository carries no structure files. Documentation examples and the evaluation tables of the design document use public entries only (unpublished cases are named by codes such as protein complex A0001).
- Double-click launchers for Windows (`launchers/windows/*.bat`) and macOS (`launchers/mac/*.command`): start the web portal on the user's own computer, or make the HTML report of every chain of a structure file (every model of an NMR file); the first start installs ifacekit into a private Python environment. Written for Windows 10/11 and macOS 12 or newer, not yet tested on those systems.
- `tools/md2docs.py` makes the plain-text and Word versions of the design document from its Markdown file (the Word file needs python-docx).
- New tools for the benchmarks: `tools/compare_pisa_batch.py` (interface area, residues and hydrogen bonds against PISA) and `tools/make_benchmark_figures.py`.

## [0.7.6] - 2026-10-28

### Added
- Clash check with hydrogens added: option `--hydrogens` (check box "Also check clashes with hydrogens added" in the portal, interface and chain reports). ifacekit places the hydrogens of the protein itself
  (`ifacekit/hydrogens.py`; no Reduce or other program): ideal geometry with the shortened X-ray X-H distances, methyl groups staggered, hydroxyl, thiol and NH3 groups turned to the best contacts, Asn, Gln and His
  flipped and (His) protonated to the best state, in three passes over the model. Only this check uses the hydrogens; every other result comes from the file as given. It reports all-atom clashes (overlap of 0.4 A or
  more, Probe radii, atoms more than four bonds apart; a polar hydrogen next to an acceptor is a hydrogen bond), a clashscore per 1000 atoms, and the flips it chose. Waters, ligands, nucleic acids and modified residues
  get no hydrogens. New result block `hydrogens` and summary fields `all_atom_clashes` and `all_atom_clashscore` (schema 0.11); tables in the pair and quality reports and lines in the text reports.
  Against MolProbity on 293 structures (alternate conformers left out) it finds 3,779 of 4,781 clash pairs = 79% (83% of those with a hydrogen) and 77% of its own 4,914 pairs are confirmed; the count per structure follows MolProbity's (r = 0.93). The settings were
  fitted on 40 of the cases and the result is the same on the other 253 (78.9% and 76.9%). Placement agrees with Reduce (97% of 4,562 hydrogens on three structures within 0.15 A) and its flip choices on 96% of 1,105 Asn and Gln amides (40 cases). Design document section 5.13.
- `tools/compare_hydrogen_clashes.py`: compares the clashes with hydrogens added with the clash lists of MolProbity, case by case, and writes the result file the fact sheet reads.
- `tools/compare_chimerax_hydrogens.py` (with `tools/chimerax_addh_clashes.py`, run inside ChimeraX): a check with hydrogens added by ChimeraX instead of ifacekit's own, scored against MolProbity (79% found and 77% confirmed with ifacekit's hydrogens; ChimeraX hydrogens find 74 to 80% but confirm only 17 to 44% of their pairs, ChimeraX's own clash finder finds 43%). The fact sheet shows it (section "Hydrogens added by ChimeraX instead", table, overlap diagrams, a column in the case table).
- `tools/explain_hydrogen_disagreement.py`: gives every clash pair on which ifacekit and MolProbity disagree its most likely cause (hydrogens on waters and ligands, hydrogen-bond rule, bond-count rule, overlap just under the limit, hydrogen placed elsewhere ...). The fact sheet shows the counts and the causes, and the design document (section 5.13) explains why 20 to 30% of the pairs differ.
- The fact sheet has a section "Clashes with hydrogens added" and a column for it in the case table (`tools/make_factsheet.py --hydrogens FILE`).

## [0.7.1] - 2026-10-28

### Added
- Alternate conformers: rotamers are scored for every conformer (rows have `altloc`, the reports show "alt B"), and the contact analysis checks the contacts of
  every further conformer (`altloc` in contacts and issues, `flags.alt_conformers_scanned`, `summary.alt_conformer_contacts`, a column "conformer" in the clash tables;
  `alt_conformers=False` switches it off). Against MolProbity, the share of its heavy-atom clashes that are flagged rises from 69.8% to 85.5% (alternate-conformer
  clashes from 41% to 96%) and the share of its rotamer outliers that are found from 59.6% to 78.6%.
- `tools/rotamer_reference/05_calibrate_cutoffs.py` and `enclosure_calibrated.json`; design document section 5.12.
- Fact sheet page (`ifacekit/data/factsheet.html`, made by `tools/make_factsheet.py`): the comparison with MolProbity and the list of the 300 test cases, with a filter and a CSV download. The notes of the reports and the portal home page link to it (portal address `/factsheet`, option `--factsheet-url`); `ifacekit --factsheet FILE` writes it to a file.

### Changed
- Rotamer cut-offs are calibrated against MolProbity's rotalyze (900 calibration entries, tests on 300 other entries): confirmed outliers 64% to 81% (recent entries) and
  76% to 92% (older entries); Leu and Ile are no longer over-called. Types with too few outliers (Cys, His, Pro, Trp) keep the earlier cut-offs. The rotamer reference is rebuilt.
- The footnotes of the reports and of the portal home page now state the comparison with MolProbity on 293 structures (single source: `ifacekit/about.py`) instead of "five structures".
- The result layout has new fields (see above), schema 0.10.

### Added
- `tools/compare_molprobity_batch.py`: compares ifacekit with `phenix.molprobity` on many structures (clashes, Ramachandran and rotamer outliers), with an
  independent test for entries released after the references were built; `tools/rama_reference/03_build_grids.py` can build the Ramachandran grids from entries
  before a cut-off date.
- Design document, section 5.11: comparison with MolProbity on 293 high-resolution X-ray structures (clashes between residues that are not bonded: 97% found;
  Ramachandran outliers: 73% found, 93% confirmed; rotamer outliers: about 60% agreement, mainly Leu and Ile over-called), and the limits that follow.

## [0.6.2] - 2026-10-28

### Changed
- The table of hydrogen bonds via water (pair and ligand reports, text report) now says in its heading and note that it lists the bridging waters only,
  waters bonded to both sides, and where waters bonded to one side only are found (not an interface bond; hydration of that chain, or "H-bonds to waters" for a ligand).

## [0.6.1] - 2026-10-28

### Added
- Check of every bridging water (hydrogen bonds via water): too short, too long, partner atoms touching each other (triangle), narrow angle.
  A bridge is *clear* or *doubtful*; the pair and ligand reports show "N clear, M doubtful" on the card and chip, a column "check" in the table of
  bridges and in the text report, and the JSON has `flags` and `quality` per bridge and `water_mediated_clear`, `water_mediated_doubtful`,
  `total_connections_clear` per interface (the existing counts are unchanged). The limits are in `ifacekit/hbonds.py` and are not calibrated:
  in a screening of 160 structures 14 to 31% of the bridges are flagged, depending on the limits (design document, section 5.10).

## [0.6.0] - 2026-10-28

### Added
- `docs/USAGE.md`, a guide to the commands with examples.
- Repository tooling for GitHub: `tools/bump_version.py` (keeps the version in all files in step and checks the changelog),
  `tools/update_install.sh` (update an installed copy from a tag or branch), tests and release workflows in `.github/`,
  `CONTRIBUTING.md` (branches, pull requests, releases).

### Changed
- The tests now find their example structures in `testdata/` of the repository, with no path from the original development machine, so all of them run in a fresh clone and in the automatic checks.

### Fixed
- Pair report (experimental and predicted): the table "Hydrogen bonds, salt bridges and metal contacts" and the counts of "interface clashes" and
  "interface issues" included every protein-water contact of the whole model, so an interface with few protein-protein contacts looked as if it had
  hundreds (for example 9lqk A-B: 204 rows, 7 of them between the chains; 7d47 A-B: 405 rows, 2 between the chains; 4 "interface clashes" and
  46 "issues" that were all, or almost all, with waters). They now list only contacts between the two chains. The numbers on the cards, the
  JSON and the text report were already correct. Waters that bridge the chains stay in their own table ("via water"), and other protein-water
  problems stay in the model quality part. The protein-protein counts agree with PISA (see the design document, section 5.10).
- `--ligand` with an ID that is not in the file printed a Python traceback. It now says so and lists the ligand IDs of the file.

## [0.5.5] - 2026-10-28

### Added
- Banner on the portal home page and in the README.
- Detailed ligand diagram: residues that hydrogen-bond or salt-bridge the ligand are drawn as chemical structures
  (`--detailed`, `--detailed-all`; two check boxes in the portal; a third view in the ligand report and extra figure files).
- Ubuntu installation guide in `server/deploy/README_DEPLOY.txt`.

### Changed
- Release layout: the folders `data/` and `results/` are empty in the release (`.gitkeep`); example structures are in `testdata/`.

## [0.5.0] - 2026-10-27   (schema 0.9)

### Added
- Hydrogen-bond counts without and with water: direct, via water (waters bonded to both sides, listed by name) and connections in total,
  in the pair report, the ligand report, the text report and the JSON (`interfaces[].hbonds`).
- Simple ligand map (default; a water only when it bonds the ligand, a residue only with a polar contact or two close contacts) and the full map
  with a switch in the report; figure files for both.
- Portal for a lab server: gunicorn launcher, periodic deletion of old jobs, a limit on simultaneous analyses, a memory limit per analysis,
  `/healthz`, systemd and nginx examples, Mol* served by the portal (`--molstar-url`), the structure keeps its own file name.
- Version, MIT licence and copyright in the footer of the portal and of every HTML report.

## [0.4.8] - 2026-10-27   (schema 0.8)

### Added
- Rotamer plots per residue type; flip check for Asn, Gln and His (a suggestion, weak against Reduce); one PDF with all Ramachandran and rotamer plots
  (`--plots-pdf`); coloured chain bars and large labels in the B-factor / pLDDT plot.

## [0.4.7] - 2026-10-27   (schema 0.7)

### Added
- Rotamer analysis with ifacekit's own reference (17 residue-type tables); held-out test of the Ramachandran reference.

## [0.4.5] - 2026-10-27   (schema 0.6)

### Added
- HTML quality report for chosen chains (`--chains ... --html`).
- Ramachandran analysis with a reference that ifacekit built from 13,115 PDB entries (replacing the first version, which used external tables).

## [0.4.0] - 2026-10-27

### Added
- Web portal (four analysis types, upload checks, fetch by PDB ID), Ramachandran analysis (first version), MIT licence.

## [0.3.0] - 2026-10-27   (schema 0.5)

First packaged version: chain overview and quality, pair analysis with interface area, contacts, clashes, waters and metal sites,
backbone geometry, AlphaFold 3 metrics (ipSAE, interface PAE, pDockQ), HTML reports with an optional Mol* 3D page, ligand report and map.
