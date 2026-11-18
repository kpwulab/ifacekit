# Using ifacekit: a guide to the commands

ifacekit looks at **one model at a time** (a PDB or mmCIF file) and tells you how its chains touch each other, where atoms are too close
or badly paired, how a ligand sits in its pocket, and how good the backbone and the side chains look. It works the same for an
experimental structure and for an AI prediction (AlphaFold 3, Boltz and similar). It uses heavy atoms only (no hydrogens) and needs no
licence key and no network.

This page is a tour of the commands. The details of how everything is calculated are in [ifacekit_design.md](ifacekit_design.md).

---

## 1. Install

```bash
git clone https://github.com/<user>/ifacekit.git
cd ifacekit
python3 -m venv .venv && . .venv/bin/activate
pip install ".[ligand,pdf]"            # core + ligand report + the PDF of plots
ifacekit --version
```

* The core needs only `gemmi`, `numpy` and `scipy`. `.[ligand]` adds `rdkit` and `cairosvg` (ligand report, PNG/PDF figures), `.[pdf]` adds
  `pypdf` (one PDF with all plots), `.[server]` adds `flask` and `gunicorn` (the web portal).
* PNG and PDF figures need the system library **libcairo**: `sudo apt install libcairo2` (Ubuntu), `brew install cairo` (macOS). Without it
  everything still works, you just get SVG only.
* Step-by-step instructions, the portal and common problems: [INSTALL.txt](../INSTALL.txt).

Try it on a public structure from the PDB:

```bash
curl -O https://files.rcsb.org/download/5HPT.cif
```
```bash
ifacekit 5HPT.cif
```

---

## 2. The idea: say what you want to look at

If you give only a file, ifacekit **lists the chains** and tells you what to ask for. It never analyses anything you did not choose.

```
$ ifacekit 5HPT.cif

Chains in this file:
  chain  type     residues  atoms  ligands  waters
  A      protein  373       3111   -        3
  B      protein  76        607    -        0
  ...
```

Then choose one of four kinds of analysis:

| You want to know | Command | Result |
|---|---|---|
| **1. How a predicted complex looks** (AI confidence, interface, quality) | `ifacekit model.cif --pair A B --html report.html --af3-summary S.json --af3-full F.json` | HTML report with ipTM, ipSAE, PAE, pLDDT, interface contacts, model quality |
| **2. How two chains of an experimental structure touch** | `ifacekit model.cif --pair A B --html report.html --experimental` | HTML report: interface area, hydrogen bonds, clashes, waters, B-factors, model quality |
| **3. How a ligand binds** | `ifacekit model.cif --ligand A:JQ1201 --html lig.html` | HTML report with interaction map, contact tables, contact quality; figure files |
| **4. How good a chain (or all chains) is** | `ifacekit model.cif --chains A --html quality.html` | HTML quality report: clashes, backbone geometry, Ramachandran, rotamers, waters, metals |

Add `-q` to switch off the text summary on screen, and `-o result.json` to write everything to a JSON file.

---

## 3. The four analyses

### 3.1 Chain quality: `--chains`

```bash
ifacekit model.cif --chains A                 # text report for chain A
ifacekit model.cif --chains A C F             # each chain on its own
ifacekit model.cif --chains all               # every chain, with a comparison table first
ifacekit model.cif --chains D E --html quality.html     # HTML report, one section per chain
```

Every chain is analysed **together with its own ligands, glycans and waters** (assigned by bond or contact, not by chain label;
`--attach label` uses the chain label instead). Contacts with *other* chains are not part of this analysis: use `--pair` for those.

You get: atom-pair contacts and clashes inside the chain, backbone bond and angle deviations, **Ramachandran** categories, **rotamer**
categories, a check of Asn/Gln/His flips (a suggestion only), waters that look unsupported and metal sites.

Extra files for the chosen chains:

```bash
ifacekit model.cif --chains A C --rama-html rama.html        # Ramachandran plots, one panel per residue class
ifacekit model.cif --chains A C --rotamer-html rotamer.html  # rotamer plots, outliers, flip check
ifacekit model.cif --chains A C --plots-pdf plots.pdf        # all Ramachandran and rotamer plots in one PDF
```

### 3.2 Interface between two chains: `--pair A B`

```bash
ifacekit model.cif --pair A B                          # text report
ifacekit model.cif --pair A B -o ab.json               # also a JSON file
ifacekit model.cif --pair A B --html report.html       # HTML report
ifacekit model.cif --pair A B --html report.html --molstar   # plus report_3d.html with a 3D viewer
```

The report gives the buried **interface area**, polar fraction, **hydrogen bonds counted three ways** (direct; via water, meaning waters that
bond to both chains; connections in total), salt bridges, hydrophobic and van der Waals contacts, ring interactions, buried polar atoms
without a partner, the most buried residues, and contacts that look wrong (clashes, too short hydrogen bonds, polarity mismatches).

*Experimental or predicted?* ifacekit **detects** this from the file. To force it: `--experimental` (B-factors, no AI metrics) or
`--predicted` (the confidence column is read as pLDDT).

### 3.3 AI-predicted complexes: AlphaFold 3 files

AlphaFold 3 gives, next to the model, two JSON files. Give them to get the AI metrics (ipTM, ipSAE, interface PAE, pDockQ, PAE map):

```bash
ifacekit fold_x_model_0.cif --pair A B --html report.html \
    --af3-summary fold_x_summary_confidences_0.json --af3-full fold_x_full_data_0.json
```

If the two files sit next to the model and have the same stem in their names, ifacekit finds them by itself. Without them you still get the
interface and quality sections, with pLDDT from the model file, but no PAE plot or ipSAE.

### 3.4 Ligands: `--ligand`

```bash
ifacekit model.cif --ligand --html lig.html            # the first ligand in the file
ifacekit model.cif --ligand A:JQ1201 --html lig.html   # a given ligand
ifacekit model.cif --ligand A:JQ1201 --html lig.html --detailed   # also a diagram with residues drawn as chemical structures
```

The **ligand ID** is `chain:name+residue number`, for example `A:JQ1201`. If you give an ID that is not in the file, ifacekit lists the
valid ones. The report has:

* an **interaction map** (a modern view and a LigPlot-style view). By default it is *simple*: a water appears only if it hydrogen-bonds the
  ligand, and a residue only if it has a polar contact or at least two close contacts. A switch shows **all contacts**. The tables always list
  every residue;
* hydrogen bonds from the ligand to the protein, to waters, and **via water**, salt bridges, hydrophobic contacts, buried polar atoms
  without a partner, clashes;
* figure files next to the report (SVG, PNG 300 dpi, PDF): `lig_map_modern.*`, `lig_map_classic.*`, `lig_map_all_*` (all contacts) and,
  with `--detailed`, `lig_map_detailed.*`.

If the file carries the chemical dictionary of the ligand (files from the wwPDB do) its exact bonds are used; otherwise bond orders are
guessed from the geometry and the report says so. Check the drawing against the real compound.

---

## 4. All the options

| Option | Meaning |
|---|---|
| `--chains ID [ID ...]` / `--chains all` | quality of chains, each on its own |
| `--pair A B` | interaction between two chains |
| `--all` | every chain and every ligand as units (large output) |
| `--ligand [ID]` | ligand report (with `--html`) |
| `--html FILE` | write a self-contained HTML report |
| `--molstar` | with `--html`: also `FILE_3d.html` with a Mol* 3D panel (needs internet unless `--molstar-url` is given) |
| `--molstar-url URL` | where Mol* is served from, for example `/static/molstar/` on the portal |
| `--experimental` / `--predicted` | say what kind of model it is (otherwise detected) |
| `--af3-summary JSON` / `--af3-full JSON` | AlphaFold 3 confidence files |
| `--detailed` / `--detailed-all` | ligand: also the diagram with structures (only H-bonding and salt-bridging residues, or all) |
| `--hydrogens` | also check clashes with hydrogens added (see below); every other result uses the file as given |
| `--rama-html FILE` / `--rotamer-html FILE` | with `--chains`: stand-alone Ramachandran or rotamer page |
| `--plots-pdf FILE` | with `--chains` or `--pair`: one PDF with all Ramachandran and rotamer plots |
| `-o FILE` | write the full result as JSON |
| `-q` | no text summary on screen |
| `--attach {contact,label}` | how ligands, glycans and waters are given to chains (default `contact`) |
| `--model N` | model number in a multi-model file (default 0; only one model is analysed) |
| `--ensemble` | multi-model file (NMR ensemble): analyse every model; prints a table with one row per model, the mean and standard deviation, and then the complete analysis of the first model. Works with `--chains`, `--pair` and `--all`; `--html`, `--rama-html`, `--rotamer-html` and `--plots-pdf` are for the first model |
| `--models LIST` | which models get a report: `first` (default), `all`, or model numbers as the file counts them (`MODEL 1`, `2`, ...), e.g. `3`, `1,4,7` or `2-5`; implies `--ensemble`. One model: the report of that model. Several models with `--html`: one page in the frame of the normal report with a bar of buttons, the ensemble summary and one button per model that swaps the tables and sections for those of that model; with `--molstar` also `FILE_3d.html` with the Mol* panel that follows the model (and shows all models together on request) |
| `--separate` | with `--models` and `--html`: one HTML file per model (`FILE_model03.html`) instead of one page |
| `--ensemble-csv FILE` | with `--ensemble`: one row per model plus mean, sd, min and max as CSV (one file per chain when several chains are given: `FILE` gets the chain name before the extension) |
| `--workers N` | with `--ensemble`: number of models analysed in parallel (default 1) |
| `--min-area A2` | smallest interface to report (default 10) |
| `--points N` | sphere points for the surface area (default 960) |
| `--all-intra` | also record ordinary contacts inside a chain (large JSON) |
| `--version`, `-h` | version, help |

---

## 5. Reading the results

**Contact verdicts.** Every pair of atoms close together gets a *type* (hydrogen bond, salt bridge, hydrophobic, van der Waals, metal,
polarity mismatch, like charge) and a *verdict*: `favorable`, `ok`, `borderline`, `unfavorable`, with a short reason. Green is acceptable, amber
means look at it, red likely a problem. The cut-offs are rules of thumb. "Unfavorable per 1000 atoms" is **not** a MolProbity clashscore.

**Hydrogen bonds.** A hydrogen bond needs geometry (no hydrogens are used): donor-acceptor distance of at most 3.3 Å, the donor pointing to the acceptor (the expected N-H direction, opposite to the sum of the directions to the donor's heavy neighbours, within 50° of the line to the acceptor; for a donor with one neighbour the angle N-X···A is at least 100°) and every angle neighbour-A···X at the acceptor at least 100°. The older distance-only count (up to 3.5 Å, no angle) is still given as *polar contacts*. *Direct* = between the two sides, no water. *Via water* = number of waters that bond to both sides (a bridge counts once). *Connections in total*
= direct + via water. A bridge through two waters is not counted. A model without waters (an AI prediction) says that bonds via water cannot be counted.
The table of *via water* bonds lists **bridging waters only**, those bonded to both sides; a water bonded to one side only is not an interface bond and is not listed (it is part of the hydration of that chain, see the quality part of the report; for a ligand it is counted under "H-bonds to waters").
Every bridge is also checked with four simple tests: a bond shorter than 2.5 A (too close for a water), longer than 3.3 A (weak), the two partner atoms touching each other
(a triangle, the water only sits beside a direct contact), or a narrow angle at the water (under 80 degrees). Bridges without a flag are **clear**, the others **doubtful**, shown as
"7 clear, 4 doubtful" and in a column "check" of the table. The limits are simple and not calibrated: a doubtful bridge is a reason to look at the density.

**Ramachandran.** Each residue is *favored*, *allowed* or *outlier* in its class (general, Gly, cis-Pro, trans-Pro, pre-Pro, Ile/Val), judged against a
reference that ifacekit built from 13,115 high-resolution PDB entries. Goals: more than 98 % favored, under 0.2 % outliers.

**Rotamers.** Side chains of 17 residue types are *favored*, *allowed* or *outlier* by their chi angles, every alternate conformer on its own (goal: under 0.3 % outliers).
Asn, Gln and His are also scored with the group flipped.

**Asn/Gln/His flips.** The ends of an Asn or Gln amide, and the sides of a His ring, cannot be told apart in a map. ifacekit scores the
hydrogen bonds around the group as modelled and flipped and **suggests** a flip when the flipped group fits clearly better. This test is weaker than
Reduce (see the design document, section 5.9). Always look at the density.

**pLDDT / B-factor plot.** Each chain has its own colour and a labelled bar under the axis; red ticks mark interface residues. Scroll to zoom,
drag to pan, hover to read a residue.

**The JSON file** (`-o`) has `schema_version`, `tool`, `input`, `units`, `interfaces` (each with `area`, `contacts`, `verdicts`, `hbonds`,
`residues_a/b`, `bridging_waters`), `contacts`, `issues`, `geometry`, `waters`, `ramachandran`, `rotamer`, `flips`, `summary`, and `hydrogens` (null unless `--hydrogens` was given).

**Clashes with hydrogens added** (`--hydrogens`, or the box in the portal). The contact tables above see heavy atoms only, and three quarters of MolProbity's clashes involve a hydrogen.
With this option ifacekit places hydrogens on the protein itself (no other program is called) and adds a table of all-atom clashes (an overlap of 0.4 A or more) to the pair and
quality reports, lines to the text report, and the block `hydrogens` to the JSON (`clashes`, `clashscore` per 1000 atoms, the `flips` it chose). Only this check uses the hydrogens;
all other results come from the file as given. Waters, ligands, nucleic acids and modified residues get no hydrogens, so their clashes with hydrogens are not seen.
Example: `ifacekit model.cif --chains A --hydrogens --html quality.html`. On 293 structures it found 79% of MolProbity's clash pairs (alternate conformers left out) and 77% of its own pairs were confirmed
(design document, sections 3.18 and 5.13, and the fact sheet). Alternate conformers are not scanned by this check: only the first conformer of every atom is used.

---

**How it compares with MolProbity.** On 293 high-resolution X-ray structures ifacekit flagged 86% of the heavy-atom clashes that MolProbity reports
(94% of those between two different residues; clashes that involve hydrogens, three quarters of MolProbity's, are checked only with `--hydrogens`, which found 79% of its clash pairs and confirmed 77% of its own), found 73% of its
Ramachandran outliers and 79% of its rotamer outliers, and 93% and 84% of its own outliers were confirmed. Low-resolution and cryo-EM models,
predictions and nucleic acids were not part of the comparison. Details: [ifacekit_design.md](ifacekit_design.md), sections 5.11 and 5.12.
The **fact sheet** page shows the comparison and the list of all test cases: the portal serves it at `/factsheet` (the notes of the reports link to it), and
`ifacekit --factsheet facts.html` writes it to a file.

## 6. From Python

```python
import ifacekit

ifacekit.list_chains("model.cif")                       # one row per chain
res = ifacekit.analyze("model.cif", chain="A")          # chain quality
res = ifacekit.analyze("model.cif", pair=("A", "B"))    # interface
res["interfaces"][0]["hbonds"]                          # {'direct': ..., 'water_mediated': ..., 'total_connections': ...}
res["ramachandran"]["percent"]                          # {'favored': ..., 'allowed': ..., 'outlier': ...}
```

The result is a plain dictionary (it can be written as JSON).

---

## 7. The web portal

For people who prefer a browser: upload a file or give a PDB ID, pick one of the four analyses, get the report.

```bash
pip install ".[server]"
server/ifacekit-server.sh --host 127.0.0.1 --port 8000        # then open http://127.0.0.1:8000
```

The page has the same four analyses. After the file is read, it lists the chains and ligands that are really in it, so you choose from a list.
Archives (zip, gz) are refused, every upload is checked, and old jobs are deleted after 7 days. To run it as a service for a group, with limits, a
systemd unit and an nginx example, see [server/deploy/README_DEPLOY.txt](../server/deploy/README_DEPLOY.txt). The portal has no login: put it on a
private network.

---

## 8. Good to know, and when something goes wrong

* **One model unless `--ensemble`, no lattice.** Only the first model of a multi-model file is analysed (`--model N` to choose) unless `--ensemble` is given, which analyses every model. Crystal contacts are not considered.
* **Alternate conformers**: the contacts and the rotamers are checked for every conformer (rows and contacts carry the conformer letter, shown as "alt B" or in a column "conformer"); the interface areas, hydrogen bonds, waters, Ramachandran and the flip check use the first conformer. **Hydrogens** in the file are ignored.
* **`chain X not found`**: use the chain letters shown by `ifacekit model.cif`. For files from the AlphaFold server these are usually A, B, ...
* **`ligand ... not found`**: the message lists the ligand IDs of the file.
* **No PNG or PDF files next to the ligand report**: install libcairo (section 1).
* **The 3D page stays empty**: it loads Mol* from the internet unless you made it with `--molstar-url`.
* **Results for poorly refined or predicted models**: the Ramachandran and rotamer references come from well-refined, low-B residues;
  low-resolution models will show more outliers by design.
* More about every method, its tests and its limits: [ifacekit_design.md](ifacekit_design.md).

---

## 9. Licence and citation

MIT licence (see [LICENSE](../LICENSE)). If ifacekit is useful in your work, a citation is suggested: see [CITATION.cff](../CITATION.cff).


## NMR ensembles

Model numbers in `--models` start at 1, as in the PDB file (`--model N`, the older option, counts from 0).

```
ifacekit 9M8X.cif --chains A --ensemble --hydrogens --workers 4 --ensemble-csv ensemble_A.csv -o ensemble_A.json
```

The text output starts with a table (one row per model: Ramachandran favored and outlier %, rotamer outlier %, bond and angle RMSZ, heavy-atom unfavorable and borderline contacts, clashes and clashscore with hydrogens added,
and for interfaces the area and the direct hydrogen bonds), followed by the mean and the standard deviation, and then the full report of the first model. The JSON file has `mode: "ensemble"` and one entry per chain or pair in
`ensembles`, each with `models` (the table), `average` (mean, sd, min, max, n of every number) and `first_model` (the complete single-model result). The standard deviation is the sample standard deviation (n - 1).
With `--html` the report shows the tables of the first model and, instead of the B-factor or pLDDT plot, the deviation of every residue from the mean structure: a line with the mean over the models of the root-mean-square distance of the backbone atoms (N, CA, C, O) from their mean positions after superposing each chain on the CA atoms of its least mobile 70% of residues (iterated three times), a dark band of one standard deviation, and a light band from the lowest to the highest model; the RMSD is in the tooltip. `--ensemble-csv FILE` also writes these numbers and the RMSD of each residue (backbone, CA only and all heavy atoms) to `FILE_residues.csv`. With `--ensemble` the Ramachandran page (`--rama-html`), the rotamer page (`--rotamer-html`) and the plots PDF (`--plots-pdf`) show every model of the file in one plot (smaller, see-through points; hover for the model; the numbers count residue-model pairs); the page made by `--models` holds the same plots in its summary, while the report of one model shows the plots of that model.
The numbers are those of the single-model analysis of each model; nothing is averaged before the analysis, so a mean is the mean of per-model results, not the analysis of an averaged structure.
Models are analysed independently: the contacts between models (as in an ensemble of alternative states) are not considered.

```
ifacekit 9M8X.cif --chains A --models all --html 9M8X_ensemble.html --workers 8     # one page: summary, a button per model that swaps the tables; add --molstar for the 3D page
ifacekit 9M8X.cif --chains A --models 3 --html model3.html                          # report of model 3
ifacekit 9M8X.cif --chains A --models 2-4 --separate --html m.html                  # m_model02.html, m_model03.html, m_model04.html
```

In the portal the check box "Analyse every model of the file" appears for chain quality and experimental interface when the uploaded file holds several models; it runs `--ensemble --models all` and offers the ensemble table and the RMSD of every residue as CSV files. The 3D (Mol*) page is not made for ensembles.
