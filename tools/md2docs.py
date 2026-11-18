#!/usr/bin/env python3
"""Make the plain-text and Word versions of a Markdown document (the design document of ifacekit).

  python tools/md2docs.py docs/ifacekit_design.md docs/ifacekit_design.txt docs/ifacekit_design.docx

Handles what the document uses: a YAML header (title, subtitle, author), headings (#, ##, ###), paragraphs, bullet lists (* or -), numbered lists,
pipe tables, fenced code blocks, **bold**, *italic* and `code`. The Word file needs python-docx (pip install python-docx); without it only the text file is written.
"""
import re
import sys
import textwrap

WIDTH = 92
AUTHOR = "Kuen-Phon Wu"


def parse(md):
    """Blocks: ("meta", dict), ("h", level, text), ("p", text), ("ul", [text]), ("ol", [(n, text)]), ("table", [rows]), ("code", [lines])."""
    lines = md.replace("\r\n", "\n").split("\n")
    meta = {}
    if lines and lines[0].strip() == "---":
        end = lines.index("---", 1)
        for ln in lines[1:end]:
            if ":" in ln:
                k, v = ln.split(":", 1)
                meta[k.strip()] = v.strip().strip('"')
        lines = lines[end + 1:]
    blocks = [("meta", meta)]
    i = 0
    while i < len(lines):
        ln = lines[i]
        if not ln.strip():
            i += 1
        elif ln.startswith("```"):
            j = i + 1
            code = []
            while j < len(lines) and not lines[j].startswith("```"):
                code.append(lines[j])
                j += 1
            blocks.append(("code", code))
            i = j + 1
        elif re.match(r"#{1,4} ", ln):
            m = re.match(r"(#{1,4}) (.*)", ln)
            blocks.append(("h", len(m.group(1)), m.group(2).strip()))
            i += 1
        elif ln.lstrip().startswith("|"):
            rows = []
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                cells = [c.strip() for c in re.split(r"(?<!\\)\|", lines[i].strip().strip("|"))]
                cells = [c.replace("\\|", "|") for c in cells]
                if not all(re.fullmatch(r":?-{2,}:?", c) for c in cells if c):
                    rows.append(cells)
                i += 1
            blocks.append(("table", rows))
        elif re.match(r"[*-] ", ln):
            items = []
            while i < len(lines) and (re.match(r"[*-] ", lines[i]) or (lines[i].startswith("  ") and lines[i].strip() and items)):
                if re.match(r"[*-] ", lines[i]):
                    items.append(lines[i][2:].strip())
                else:
                    items[-1] += " " + lines[i].strip()
                i += 1
            blocks.append(("ul", items))
        elif re.match(r"\d+\. ", ln):
            items = []
            while i < len(lines) and (re.match(r"\d+\. ", lines[i]) or (lines[i].startswith("   ") and lines[i].strip() and items)):
                m = re.match(r"(\d+)\. (.*)", lines[i])
                if m:
                    items.append((m.group(1), m.group(2).strip()))
                else:
                    items[-1] = (items[-1][0], items[-1][1] + " " + lines[i].strip())
                i += 1
            blocks.append(("ol", items))
        else:
            para = [ln.strip()]
            i += 1
            while i < len(lines) and lines[i].strip() and not re.match(r"(#{1,4} |```|[*-] |\d+\. |\|)", lines[i]):
                para.append(lines[i].strip())
                i += 1
            blocks.append(("p", " ".join(para)))
    return blocks


def plain(t):
    t = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1 (\2)", t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"\1", t)
    t = re.sub(r"(?<![\w*])\*([^*\n]+)\*(?![\w*])", r"\1", t)
    return t.replace("`", "")


def to_txt(blocks):
    out = []
    for b in blocks:
        k = b[0]
        if k == "meta":
            m = b[1]
            out += [m.get("title", ""), m.get("author", ""), "", ""] if m else []
        elif k == "h":
            out += ["", plain(b[2]), ""] if b[1] <= 2 else ["", plain(b[2])]
            if b[1] > 2:
                out.append("")
        elif k == "p":
            out += textwrap.wrap(plain(b[1]), WIDTH) + [""]
        elif k == "ul":
            for t in b[1]:
                out += textwrap.wrap(plain(t), WIDTH, initial_indent="  - ", subsequent_indent="    ")
            out.append("")
        elif k == "ol":
            for n, t in b[1]:
                out += textwrap.wrap(plain(t), WIDTH, initial_indent=f"  {n}. ", subsequent_indent=" " * (len(n) + 5))
            out.append("")
        elif k == "code":
            out += ["    " + c for c in b[1]] + [""]
        elif k == "table":
            rows = [[plain(c) for c in r] for r in b[1]]
            n = max(len(r) for r in rows)
            rows = [r + [""] * (n - len(r)) for r in rows]
            w = [max(len(r[i]) for r in rows) for i in range(n)]
            fmt = lambda r: "  " + "  ".join(c.ljust(w[i]) for i, c in enumerate(r)).rstrip()
            out.append(fmt(rows[0]))
            out.append("  " + "  ".join("-" * x for x in w))
            out += [fmt(r) for r in rows[1:]]
            out.append("")
    return re.sub(r"\n{4,}", "\n\n\n", "\n".join(out)).rstrip() + "\n"


def add_runs(par, text, size=None, code_font="Courier New"):
    from docx.shared import Pt
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1 (\2)", text)
    for part in re.split(r"(\*\*[^*]+\*\*|`[^`]+`|(?<![\w*])\*[^*\n]+\*(?![\w*]))", text):
        if not part:
            continue
        if part.startswith("**"):
            r = par.add_run(part[2:-2])
            r.bold = True
        elif part.startswith("`"):
            r = par.add_run(part[1:-1])
            r.font.name = code_font
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            r = par.add_run(part[1:-1])
            r.italic = True
        else:
            r = par.add_run(part)
        if size:
            r.font.size = Pt(size)


def to_docx(blocks, path):
    from docx import Document
    from docx.shared import Cm, Pt
    doc = Document()
    for sec in doc.sections:
        sec.left_margin = sec.right_margin = Cm(2.0)
        sec.top_margin = sec.bottom_margin = Cm(2.0)
    st = doc.styles["Normal"]
    st.font.name = "Arial"
    st.font.size = Pt(10)
    for b in blocks:
        k = b[0]
        if k == "meta" and b[1]:
            doc.add_heading(b[1].get("title", ""), level=0)
            for key in ("subtitle", "author"):
                if b[1].get(key):
                    p = doc.add_paragraph()
                    add_runs(p, b[1][key])
        elif k == "h":
            doc.add_heading(plain(b[2]), level=min(b[1], 3))
        elif k == "p":
            add_runs(doc.add_paragraph(), b[1])
        elif k == "ul":
            for t in b[1]:
                add_runs(doc.add_paragraph(style="List Bullet"), t)
        elif k == "ol":
            for n, t in b[1]:
                p = doc.add_paragraph()
                p.paragraph_format.left_indent = Cm(0.9)
                p.paragraph_format.first_line_indent = Cm(-0.6)
                add_runs(p, f"{n}. {t}")
        elif k == "code":
            for c in b[1]:
                p = doc.add_paragraph()
                p.paragraph_format.space_after = Pt(0)
                p.paragraph_format.left_indent = Cm(0.6)
                r = p.add_run(c)
                r.font.name = "Courier New"
                r.font.size = Pt(8.5)
            doc.add_paragraph()
        elif k == "table":
            rows = b[1]
            n = max(len(r) for r in rows)
            t = doc.add_table(rows=len(rows), cols=n)
            t.style = "Table Grid"
            for i, r in enumerate(rows):
                for j in range(n):
                    cell = t.cell(i, j)
                    cell.text = ""
                    par = cell.paragraphs[0]
                    add_runs(par, r[j] if j < len(r) else "", size=8.5)
                    if i == 0:
                        for run in par.runs:
                            run.bold = True
            doc.add_paragraph()
    cp = doc.core_properties
    cp.author = AUTHOR
    cp.last_modified_by = AUTHOR
    for b in blocks:
        if b[0] == "meta" and b[1].get("title"):
            cp.title = b[1]["title"]
    cp.comments = ""
    doc.save(path)


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    blocks = parse(open(argv[0], encoding="utf-8").read())
    open(argv[1], "w", encoding="utf-8").write(to_txt(blocks))
    print("written", argv[1])
    if len(argv) > 2:
        try:
            to_docx(blocks, argv[2])
            print("written", argv[2])
        except ImportError:
            print("python-docx is not installed: the Word file was not written (pip install python-docx)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
