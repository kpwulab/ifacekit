"""Tests of the repository tooling: version bump and check, and the message for a wrong ligand ID."""
import importlib.util
import os
import shutil
import subprocess
import sys

import pytest

from pdbdata import need

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _tool():
    spec = importlib.util.spec_from_file_location("bump_version", os.path.join(ROOT, "tools", "bump_version.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_versions_agree_in_this_checkout():
    t = _tool()
    assert t.check(ROOT) == [], t.check(ROOT)
    import ifacekit
    assert set(t.read_versions(ROOT).values()) == {ifacekit.__version__}


def _copy_tree(tmp_path):
    t = _tool()
    for rel, _ in t.PLACES + [("CHANGELOG.md", None), ("ifacekit/core.py", None)]:
        dst = tmp_path / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(os.path.join(ROOT, rel), dst)
    return t


def test_bump_changes_every_place_and_the_changelog(tmp_path):
    t = _copy_tree(tmp_path)
    old = next(iter(set(t.read_versions(str(tmp_path)).values())))
    t.bump("9.8.7", root=str(tmp_path), schema="9.9")
    assert set(t.read_versions(str(tmp_path)).values()) == {"9.8.7"}
    assert t.check(str(tmp_path)) == []                                    # the changelog got its section
    log = (tmp_path / "CHANGELOG.md").read_text()
    assert "## [Unreleased]\n\n## [9.8.7] - " in log
    assert 'SCHEMA_VERSION = "9.9"' in (tmp_path / "ifacekit/core.py").read_text()


def test_check_finds_a_mismatch_and_a_missing_changelog_section(tmp_path):
    t = _copy_tree(tmp_path)
    p = tmp_path / "pyproject.toml"
    p.write_text(p.read_text().replace('version = "', 'version = "1.', 1))
    assert any("differ" in x for x in t.check(str(tmp_path)))
    t2 = _copy_tree(tmp_path)
    log = tmp_path / "CHANGELOG.md"
    log.write_text(log.read_text().replace("## [", "## [x", 3))
    assert any("CHANGELOG" in x for x in t2.check(str(tmp_path)))
    with pytest.raises(ValueError):
        t.bump("not-a-version", root=str(tmp_path))


def test_wrong_ligand_id_gives_a_clean_message():
    pytest.importorskip("rdkit")
    f = need("4QZS.cif")
    r = subprocess.run([sys.executable, "-m", "ifacekit", f, "--ligand", "A:NOPE", "--html", os.devnull + ".html"],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 2
    assert "Traceback" not in r.stderr and "A:JQ1201" in r.stderr and "not found" in r.stderr


def test_workflow_files_are_valid_and_scripts_executable():
    yaml = pytest.importorskip("yaml")
    for f in ("tests.yml", "release.yml"):
        d = yaml.safe_load(open(os.path.join(ROOT, ".github", "workflows", f)))
        assert d["jobs"]
    for f in ("tools/bump_version.py", "tools/update_install.sh", "server/ifacekit-server.sh"):
        assert os.access(os.path.join(ROOT, f), os.X_OK), f
    for f in ("CHANGELOG.md", "CONTRIBUTING.md", "docs/USAGE.md", ".gitignore", "LICENSE"):
        assert os.path.exists(os.path.join(ROOT, f)), f


def test_readme_type_files_carry_no_version_or_lab_path():
    """The README-type files say just `ifacekit` and use /your_path in examples, so a release never has to edit them."""
    import ifacekit
    files = ["README.md", "INSTALL.txt", "docs/USAGE.md", "server/deploy/README_DEPLOY.txt", "server/deploy/ifacekit.service",
             "server/deploy/nginx-ifacekit.conf.example", "ifacekit/data/README.txt"]
    for rel in files:
        text = open(os.path.join(ROOT, rel), encoding="utf-8").read()
        assert ifacekit.__version__ not in text, rel
        for bad in ("/opt/software", "/home/", "/data/rama_source", "ifacekit_project"):
            assert bad not in text, (rel, bad)
    assert open(os.path.join(ROOT, "INSTALL.txt")).readline().strip() == "ifacekit - installation guide"
    assert "/your_path/ifacekit" in open(os.path.join(ROOT, "server/deploy/ifacekit.service")).read()


def test_double_click_launchers():
    """The launchers exist, keep the line endings and the executable bit that Windows and macOS need, and carry no path of the author's computers."""
    import os
    import shutil
    import subprocess
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    win = ["launchers/windows/ifacekit-portal.bat", "launchers/windows/ifacekit-check-file.bat"]
    mac = ["launchers/mac/ifacekit-portal.command", "launchers/mac/ifacekit-check-file.command"]
    for f in win:
        raw = open(os.path.join(root, f), "rb").read()
        assert raw.count(b"\r\n") == raw.count(b"\n") > 10, f          # every line ends with CR LF
    for f in mac:
        p = os.path.join(root, f)
        raw = open(p, "rb").read()
        assert b"\r" not in raw and raw.startswith(b"#!/bin/bash"), f
        assert os.access(p, os.X_OK), f"{f} is not executable"
        if shutil.which("bash"):
            assert subprocess.run(["bash", "-n", p]).returncode == 0, f
    for f in win + mac + ["launchers/README.txt"]:
        t = open(os.path.join(root, f), encoding="utf-8").read()
        for bad in ("/opt/software", "/home/", "/data/", "ifacekit_project", "kpwu", "192.168"):
            assert bad not in t, (f, bad)
    gitattr = open(os.path.join(root, ".gitattributes")).read()
    assert "*.bat text eol=crlf" in gitattr and "*.command text eol=lf" in gitattr
