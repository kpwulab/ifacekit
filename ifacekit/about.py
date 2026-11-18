"""Version, licence and citation text shown in the footers of the reports and of the portal."""
import html
import os

from . import __version__

FACTSHEET_URL = os.environ.get("IFACEKIT_FACTSHEET_URL") or None   # where the fact sheet page is served (set by --factsheet-url; the portal uses /factsheet)

LICENCE = "MIT licence"
COPYRIGHT = "© 2026 Kuen-Phon Wu, Institute of Biological Chemistry, Academia Sinica"
CITE = "If ifacekit is useful in your work, we suggest citing: [authors, title, bioRxiv preprint, DOI: to be added]."


# What the reports say about how ifacekit compares with MolProbity (numbers: design document, section 5.11, 293 structures)
COMPARISON = ("How it compares with MolProbity: on 293 high-resolution X-ray structures ifacekit flagged 86% of the heavy-atom clashes that MolProbity reports "
              "(94% of those between two different residues), and it found 73% of MolProbity's Ramachandran outliers and 79% of its rotamer outliers, while 93% and 84% "
              "of its own outliers were confirmed (design document, section 5.11, and the fact sheet). Clashes that involve hydrogens are three quarters of MolProbity's: they are checked only with the option that adds hydrogens, "
              "which found 79% of MolProbity's clash pairs (alternate conformers left out) and confirmed 77% of its own. Waters, ligands and nucleic acids get no hydrogens. "
              "Low-resolution and cryo-EM models, predictions and nucleic acids were not part of the comparison.")
COMPARISON_SHORT = ("Heavy-atom contacts, Ramachandran and rotamer outliers were compared with MolProbity on 293 high-resolution X-ray structures "
                    "(design document, section 5.11); clashes that involve hydrogens are checked only when hydrogens are added (option --hydrogens: 79% of MolProbity's clash pairs found).")


def factsheet_path():
    """The packaged fact sheet page: the comparison with MolProbity and the list of test cases."""
    return os.path.join(os.path.dirname(__file__), "data", "factsheet.html")


def _link(url=None):
    url = url or FACTSHEET_URL
    if url:
        return f" <a href='{html.escape(url)}' target='_blank' rel='noopener'>Fact sheet: the comparison and the list of test cases</a>."
    return " The fact sheet with the comparison and the list of test cases can be written to a file with: ifacekit --factsheet FILE.html."


def comparison_html(url=None):
    """The note on the comparison with MolProbity for the end of a report, with the link to the fact sheet when its address is known."""
    return COMPARISON + _link(url)


def comparison_short_html(url=None):
    return COMPARISON_SHORT + _link(url)


def footer_html():
    return (f"<footer style='margin-top:32px;padding-top:10px;border-top:1px solid var(--line);font-size:12.5px;color:var(--mut)'>"
            f"ifacekit {__version__} &middot; {LICENCE} &middot; {COPYRIGHT}.<br>{CITE}</footer>")
