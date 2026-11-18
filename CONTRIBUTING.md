# Working on ifacekit: branches, updates and releases

## The idea in one paragraph
`main` always works and always has a version number. You change things on a short-lived branch, open a pull request, the automatic tests run,
and you merge. When a set of changes is ready to be a release you set the new version, tag it, and GitHub builds the release files.
Users and servers update by switching to a newer tag.

## Branches
| Branch | What it is | Rule |
|---|---|---|
| `main` | the stable line; every release is a tag on it | never commit directly; changes come by pull request; protect it in the GitHub settings (require the tests to pass) |
| `feature/<topic>` | new work, for example `feature/usage-statistics` | branch from `main`, merge back by pull request, delete afterwards |
| `fix/<topic>` | a bug fix, for example `fix/ligand-id-error` | same |
| `docs/<topic>` | documentation only | same |
| `release/X.Y.x` | only when an old release needs a bug fix after `main` has moved on | branch from the tag, cherry-pick the fix, tag `vX.Y.Z` on it |

```bash
git switch main && git pull
git switch -c feature/usage-statistics          # new work
# ... edit, test ...
git add -A && git commit -m "Count portal uses per analysis type"
git push -u origin feature/usage-statistics     # then open a pull request on GitHub
```

## Every change
1. Run the tests: `python -m pytest -q tests` (about two minutes; the first run downloads a few public PDB entries, about 6 MB, into `tests/.pdb_cache`).
2. Add a line to **CHANGELOG.md** under `## [Unreleased]` (Added, Changed, Fixed). Say what a user will notice.
3. If the structure of the JSON result changes, change `SCHEMA_VERSION` in `ifacekit/core.py` and say so in the changelog.
4. Update `docs/USAGE.md` and the design document `docs/ifacekit_design.md` when a command or a behaviour changes.
5. Open the pull request. The `tests` workflow runs on Python 3.9 and 3.12; merge when it is green.

## Version numbers
`MAJOR.MINOR.PATCH`. Patch: bug fixes only. Minor: new features that do not break existing commands or result fields. Major: something that
breaks them. The version is written in five files; never edit them by hand:

```bash
python tools/bump_version.py --check            # do all places agree, does CHANGELOG.md have the section?
python tools/bump_version.py X.Y.Z              # set the new version everywhere, turn [Unreleased] into [X.Y.Z]
python tools/bump_version.py X.Y.Z --schema A.B # also the JSON schema version
```

## Making a release
```bash
git switch main && git pull
python tools/bump_version.py X.Y.Z
$EDITOR CHANGELOG.md                            # read the new section, add the date and wording
# if you keep docs/ifacekit_design.txt and .docx, regenerate them from the .md
python -m pytest -q tests
git add -A && git commit -m "Release X.Y.Z"
git push                                        # (through a pull request if main is protected)
git tag -a vX.Y.Z -m "ifacekit X.Y.Z" && git push origin vX.Y.Z
```
Pushing the tag starts the `release` workflow: it checks that the tag matches the version, builds `ifacekit_X.Y.Z.tgz`, a wheel and a source
distribution, and creates the GitHub release with them. Optionally connect the repository to Zenodo (zenodo.org, GitHub integration): every
release then gets a DOI that you can cite and put in `CITATION.cff`.

## How users and servers update
| Who | How |
|---|---|
| Someone who cloned the repository | `git pull` for the development version, or `tools/update_install.sh` for the newest release tag, or `tools/update_install.sh vX.Y.5` for a given one |
| Someone who installs without cloning | `pip install --upgrade "git+https://github.com/<user>/ifacekit.git@vX.Y.Z"` (or download the archive from the release page) |
| A lab server running the portal | `tools/update_install.sh`, then restart the service (`sudo systemctl restart ifacekit`) or reload gunicorn without downtime (`kill -HUP <master process id>`) |

`tools/update_install.sh` refuses to run when you have uncommitted changes, and leaves the checkout on the tag; go back with `git switch main`.
Keep site-specific files (addresses, launcher scripts) outside the repository: the code is public.

## The fact sheet
`ifacekit/data/factsheet.html` shows the comparison with MolProbity and the list of test cases; the footnotes of the reports link to it. It is built from the result files of a
comparison and is not rebuilt by a release. After a new comparison (`tools/compare_molprobity_batch.py`, `tools/compare_hydrogen_clashes.py` for the clash check with hydrogens added, `tools/explain_hydrogen_disagreement.py` for the causes of the differences, and `tools/compare_chimerax_hydrogens.py` for the ChimeraX check, which needs ChimeraX) and a calibration, run `tools/make_factsheet.py` (its first lines give
the arguments), check the page in a browser and commit it together with the changed numbers in `ifacekit/about.py` (`COMPARISON`) and in the design document. A test compares
the numbers of the footnote with those of the page.

## Things that must never be committed
Local data and results (`data/`, `results/`, they are ignored), portal job folders, passwords or tokens, addresses of your own servers,
and structure files that are not allowed to be public. The tests use only public PDB entries, which they download when needed (`tests/pdbdata.py`); no structure file is part of the repository.
