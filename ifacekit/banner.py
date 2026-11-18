"""The iFaceKit banner as an inline image, so every HTML report stays one self-contained file that works offline."""
import base64
import os

from . import __version__

_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "banner.webp")
_CACHE = []


def banner_html():
    """<img> tag with the banner embedded as a data URI, and under it a small box with the version of ifacekit."""
    if not _CACHE:
        try:
            with open(_PATH, "rb") as fh:
                _CACHE.append("data:image/webp;base64," + base64.b64encode(fh.read()).decode("ascii"))
        except OSError:
            _CACHE.append("")
    box = (f'<div class="ifk-version" style="display:inline-block;margin:0 0 12px;padding:1px 9px;border:1px solid currentColor;border-radius:6px;'
           f'font:12px/1.6 sans-serif;opacity:.7">version {__version__}</div>')
    if not _CACHE[0]:
        return box
    return ('<img class="ifk-banner" alt="iFaceKit" width="800" height="267" src="' + _CACHE[0] +
            '" style="display:block;width:100%;max-width:520px;height:auto;margin:0 0 6px">' + box)
