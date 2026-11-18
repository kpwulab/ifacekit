#!/usr/bin/env python3
"""Show, check or set the version of ifacekit in every file that names it.

  python tools/bump_version.py --check          all places agree and CHANGELOG.md has a section for that version (used by CI)
  python tools/bump_version.py 0.6.0            set a new version everywhere; turns [Unreleased] in CHANGELOG.md into [0.6.0]
  python tools/bump_version.py 0.6.0 --schema 1.0   also set the JSON schema version (ifacekit/core.py and the design document)

Places: ifacekit/__init__.py, pyproject.toml, CITATION.cff, docs/ifacekit_design.md (first line of the header). The README-type files carry no version number.
After a bump regenerate docs/ifacekit_design.txt and .docx from the .md if you keep them: python tools/md2docs.py docs/ifacekit_design.md docs/ifacekit_design.txt docs/ifacekit_design.docx (the Word file needs python-docx).
"""
import datetime
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLACES = [
    ("ifacekit/__init__.py", r'(__version__ = ")([^"]+)(")'),
    ("pyproject.toml", r'(?m)^(version = ")([^"]+)(")'),
    ("CITATION.cff", r"(?m)^(version: )(\S+)()"),
    ("docs/ifacekit_design.md", r"(Version )(\d[\w.]*)( \(schema)"),
]
VER = re.compile(r"^\d+\.\d+\.\d+$")


def read_versions(root=ROOT):
    out = {}
    for rel, pat in PLACES:
        m = re.search(pat, open(os.path.join(root, rel), encoding="utf-8").read())
        out[rel] = m.group(2) if m else None
    return out


def check(root=ROOT):
    """Return a list of problems (empty when everything agrees)."""
    v = read_versions(root)
    problems = [f"{rel}: no version found" for rel, x in v.items() if x is None]
    found = {x for x in v.values() if x}
    if len(found) > 1:
        problems.append("versions differ: " + ", ".join(f"{rel}={x}" for rel, x in v.items()))
    if len(found) == 1:
        ver = found.pop()
        if not VER.match(ver):
            problems.append(f"'{ver}' is not of the form MAJOR.MINOR.PATCH")
        log = open(os.path.join(root, "CHANGELOG.md"), encoding="utf-8").read()
        if not re.search(r"(?m)^## \[" + re.escape(ver) + r"\]", log):
            problems.append(f"CHANGELOG.md has no section '## [{ver}]'")
    return problems


def bump(new, root=ROOT, schema=None, today=None):
    if not VER.match(new):
        raise ValueError("the version must look like 0.6.0")
    old = {x for x in read_versions(root).values() if x}
    if len(old) != 1:
        raise ValueError("the files do not agree on the current version; fix that first (python tools/bump_version.py --check)")
    old = old.pop()
    for rel, pat in PLACES:
        p = os.path.join(root, rel)
        s = open(p, encoding="utf-8").read()
        s = re.sub(pat, lambda m: m.group(1) + new + m.group(3), s, count=1)
        open(p, "w", encoding="utf-8").write(s)
    if schema:
        p = os.path.join(root, "ifacekit/core.py")
        s = open(p, encoding="utf-8").read()
        s, n = re.subn(r'(SCHEMA_VERSION = ")([^"]+)(")', lambda m: m.group(1) + schema + m.group(3), s, count=1)
        if n != 1:
            raise ValueError("SCHEMA_VERSION not found in ifacekit/core.py")
        open(p, "w", encoding="utf-8").write(s)
        p = os.path.join(root, "docs/ifacekit_design.md")
        s = open(p, encoding="utf-8").read()
        s = re.sub(r"(\(schema )([^)]+)(\))", lambda m: m.group(1) + schema + m.group(3), s, count=1)
        open(p, "w", encoding="utf-8").write(s)
    p = os.path.join(root, "CHANGELOG.md")
    s = open(p, encoding="utf-8").read()
    day = (today or datetime.date.today()).isoformat()
    if "## [Unreleased]" in s:
        s = s.replace("## [Unreleased]", f"## [Unreleased]\n\n## [{new}] - {day}", 1)
    else:
        s = s.replace("\n## [", f"\n## [{new}] - {day}\n\n(write the changes here)\n\n## [", 1)
    open(p, "w", encoding="utf-8").write(s)
    return old


def main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return 0
    if argv[0] == "--check":
        probs = check()
        for k, v in read_versions().items():
            print(f"  {k:28s} {v}")
        print("OK" if not probs else "PROBLEMS:\n  " + "\n  ".join(probs))
        return 1 if probs else 0
    schema = argv[argv.index("--schema") + 1] if "--schema" in argv else None
    old = bump(argv[0], schema=schema)
    print(f"{old} -> {argv[0]}" + (f", schema {schema}" if schema else ""))
    print("Next: check CHANGELOG.md, run the tests, commit, tag v" + argv[0] + " (see CONTRIBUTING.md).")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
