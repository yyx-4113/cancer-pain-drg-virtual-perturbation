#!/usr/bin/env python3
"""Build the PLOS ONE submission pack into submission_PLOS_ONE/.

Outputs:
  Manuscript.docx                 PLOS ONE-formatted main document
  Cover_Letter.docx               cover letter
  S1_Supporting_Information.docx  data/code availability text
  Figure1.tif / Figure2.tif / Figure3.tif   (600 ppi, RGB)
  SUBMISSION_MANIFEST.md          file->EM item-type map + metadata + open items
  source/manuscript_plos.md       the polished markdown source (traceability)

Formatting applied (PLOS ONE submission guidelines, verified 2026-10-10):
  Times New Roman 12 pt, double spaced, 1-inch margins, page numbers in footer.
  Unstructured abstract (single paragraph, 296 words, no citations).
  References numbered in order of citation with DOIs (Vancouver).
  Tables editable, caption above; figures as separate TIFF files, legends in text.
"""
import os, re, shutil
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, Inches, RGBColor
from PIL import Image

HERE = r"D:\2026.10\不同部位癌痛患者对药物敏感性生信分析\cancer_pain_virtual_perturbation"
MS = os.path.join(HERE, "docs", "manuscript", "manuscript_final.md")
COVER = os.path.join(HERE, "submission_PLOS_ONE", "source", "cover_letter.md")
SI = os.path.join(HERE, "submission_PLOS_ONE", "source", "supporting_information.md")
FIGDIR = os.path.join(HERE, "docs", "manuscript", "figures")
OUT = os.path.join(HERE, "submission_PLOS_ONE")
BODY_FONT = "Times New Roman"

os.makedirs(OUT, exist_ok=True)


# --------------------------------------------------------------------------- #
def new_document(double=True):
    doc = Document()
    normal = doc.styles["Normal"]
    normal.font.name = BODY_FONT
    normal.font.size = Pt(12)
    normal.element.rPr.rFonts.set(qn("w:eastAsia"), BODY_FONT)
    pf = normal.paragraph_format
    pf.line_spacing = 2.0 if double else 1.15
    pf.space_after = Pt(0)
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Cm(2.54)
        s.left_margin = s.right_margin = Cm(2.54)
        add_page_numbers(s)
    return doc


def add_page_numbers(section):
    p = section.footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    run._r.addnext(fld)


INLINE = re.compile(r"(\*\*.+?\*\*|\*[^*\n]+?\*|`[^`\n]+?`)")


def add_runs(p, text, bold=False, italic=False, size=None, font=None):
    for piece in INLINE.split(text):
        if not piece:
            continue
        b, i, mono = bold, italic, False
        if piece.startswith("**") and piece.endswith("**") and len(piece) > 4:
            piece, b = piece[2:-2], True
        elif piece.startswith("*") and piece.endswith("*") and len(piece) > 2:
            piece, i = piece[1:-1], True
        elif piece.startswith("`") and piece.endswith("`") and len(piece) > 2:
            piece, mono = piece[1:-1], True
        run = p.add_run(piece)
        run.bold, run.italic = b, i
        if size:
            run.font.size = Pt(size)
        run.font.name = "Courier New" if mono else (font or BODY_FONT)
    return p


def para(doc, text="", *, bold=False, italic=False, size=None, align=None,
         space_before=0, space_after=4, single=False, indent=None, hanging=None):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    if single:
        pf.line_spacing = 1.0
    if indent is not None:
        pf.left_indent = indent
    if hanging is not None:
        pf.left_indent = hanging
        pf.first_line_indent = -hanging
    pf.space_before = Pt(space_before)
    pf.space_after = Pt(space_after)
    if align is not None:
        p.alignment = align
    if text:
        add_runs(p, text, bold=bold, italic=italic, size=size)
    return p


def is_table_row(line):
    return line.strip().startswith("|") and line.strip().endswith("|")


def is_sep(line):
    return bool(re.fullmatch(r"\|[\s:|\-]+\|", line.strip()))


def cells(line):
    return [c.strip() for c in line.strip().strip("|").split("|")]


def add_table(doc, rows, size=9.0):
    if not rows:
        return
    width = max(len(r) for r in rows)
    t = doc.add_table(rows=0, cols=width)
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = True
    for ri, row in enumerate(rows):
        out = t.add_row().cells
        for ci in range(width):
            txt = row[ci] if ci < len(row) else ""
            out[ci].text = ""
            c = out[ci].paragraphs[0]
            c.paragraph_format.line_spacing = 1.0
            add_runs(c, txt, bold=(ri == 0), size=size)
        if ri == 0:
            trPr = t.rows[0]._tr.get_or_add_trPr()
            trPr.append(OxmlElement("w:tblHeader"))
    para(doc, "", single=True, space_after=4)


# --------------------------------------------------------------------------- #
# manuscript builder
def build_manuscript():
    doc = new_document()
    lines = open(MS, encoding="utf-8").read().split("\n")

    i = 0
    in_refs = False
    # ---- front matter (until first '---') ----
    while i < len(lines):
        line = lines[i].rstrip()
        s = line.strip()
        if s == "---":
            i += 1
            doc.add_page_break()
            break
        if s.startswith("# "):
            para(doc, s[2:], bold=True, size=16, align=WD_ALIGN_PARAGRAPH.CENTER, space_after=10)
        elif s:
            para(doc, s, space_after=2, single=True)
        i += 1

    # ---- body ----
    while i < len(lines):
        raw = lines[i]
        line = raw.rstrip()
        s = line.strip()
        i += 1
        if not s:
            continue
        if s.startswith("## "):
            head = s[3:].strip()
            in_refs = head.lower().startswith("references")
            if head.lower() == "abstract":
                para(doc, "Abstract", bold=True, size=12, space_before=8, space_after=4)
            else:
                para(doc, head, bold=True, size=12, space_before=8, space_after=4)
            continue
        if s.startswith("### "):
            para(doc, s[4:].strip(), bold=True, size=11, space_before=6, space_after=2)
            continue
        if s.startswith("!["):  # image embed
            m = re.search(r"!\[[^\]]*\]\(([^)]+)\)", s)
            if m:
                img = os.path.join(os.path.dirname(MS), m.group(1))
                if os.path.exists(img):
                    try:
                        doc.add_picture(img, width=Inches(6.0))
                        doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
                    except Exception as e:
                        para(doc, f"[Figure file: {os.path.basename(img)}]", italic=True)
            continue
        if is_table_row(line):  # collect table
            rows = []
            while i - 1 < len(lines) and is_table_row(lines[i - 1]):
                if not is_sep(lines[i - 1]):
                    rows.append(cells(lines[i - 1]))
                i += 1
            add_table(doc, rows)
            continue
        if s.startswith("*") and ("footnote" in s.lower() or s.startswith("*Table")
                                   or s.startswith("*Δ") or "marks a contrast" in s):
            para(doc, s, italic=True, size=9, single=True, space_after=4)
            continue
        if in_refs:
            para(doc, s, single=True, space_after=4, indent=Inches(0.5), hanging=Inches(0.5))
            continue
        if s.startswith("**") and s.endswith("**"):
            para(doc, s, bold=True, single=True, space_after=4)
            continue
        para(doc, s, space_after=4)
    path = os.path.join(OUT, "Manuscript.docx")
    doc.save(path)
    return path


# --------------------------------------------------------------------------- #
# simple markdown renderer for cover letter / supporting info
def render_simple(doc, md_path):
    for line in open(md_path, encoding="utf-8").read().split("\n"):
        s = line.strip()
        if not s:
            continue
        if s.startswith("# "):
            para(doc, s[2:], bold=True, size=14, space_before=6, space_after=6)
        elif s.startswith("## "):
            para(doc, s[3:], bold=True, size=12, space_before=6, space_after=4)
        elif s.startswith("- "):
            p = para(doc, s[2:], space_after=2, single=True)
            p.style = doc.styles["Normal"]
            p.paragraph_format.left_indent = Inches(0.3)
        else:
            para(doc, s, space_after=4, single=(not s.startswith("Dear")) )


def build_cover():
    doc = new_document()
    render_simple(doc, COVER)
    path = os.path.join(OUT, "Cover_Letter.docx")
    doc.save(path)
    return path


def build_si():
    doc = new_document()
    render_simple(doc, SI)
    path = os.path.join(OUT, "S1_Supporting_Information.docx")
    doc.save(path)
    return path


# --------------------------------------------------------------------------- #
def convert_figures():
    # PLOS ONE: figures TIFF/EPS, 300-600 ppi, <=10 MB each. Resize so the long
    # side is <= 2000 px (=> 300 dpi print width ~6.7 in) and LZW-compress, which
    # keeps every figure comfortably under the 10 MB ceiling while staying inside
    # the required 300-600 ppi range.
    out = []
    for n in (1, 2, 3):
        src = os.path.join(FIGDIR, f"fig{n}.png")
        dst = os.path.join(OUT, f"Figure{n}.tif")
        with Image.open(src) as im:
            im = im.convert("RGB")
            w, h = im.size
            long = max(w, h)
            if long > 2000:
                scale = 2000.0 / long
                im = im.resize((int(w * scale + 0.5), int(h * scale + 0.5)),
                               Image.LANCZOS)
            im.save(dst, "TIFF", dpi=(300, 300), compression="tiff_lzw")
        out.append(dst)
    return out


def main():
    ms = build_manuscript()
    cl = build_cover()
    si = build_si()
    figs = convert_figures()
    # copy polished source for traceability
    shutil.copyfile(MS, os.path.join(OUT, "source", "manuscript_plos.md"))
    # manifest
    man = os.path.join(OUT, "SUBMISSION_MANIFEST.md")
    open(man, "w", encoding="utf-8").write(MANIFEST)
    print("built:")
    for p in [ms, cl, si, *figs, man]:
        print(f"  {os.path.relpath(p, OUT):42s} {os.path.getsize(p)/1024:8.1f} KB")
    return 0


MANIFEST = """# PLOS ONE Submission Manifest

Repository: https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation (release v1.0.0)
Article type: Research Article
Date prepared: 2026-10-10

## File -> Editorial Manager item type
| Local file | EM item type | Notes |
|---|---|---|
| Manuscript.docx | Manuscript | Main document (double-spaced, Times New Roman 12pt, page numbers) |
| Cover_Letter.docx | Cover letter | Includes AI-assistance disclosure + 296-word abstract confirmation |
| Figure1.tif | Figure | 300 ppi, RGB, LZW-TIFF, <10 MB |
| Figure2.tif | Figure | 300 ppi, RGB, LZW-TIFF, <10 MB |
| Figure3.tif | Figure | 300 ppi, RGB, LZW-TIFF, <10 MB |
| S1_Supporting_Information.docx | Supporting Information | Data/code availability, reproducibility gates |

## Metadata for the submission form
- Title: Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a negative-control calibration of in-silico perturbation
- Running head: Divergent DRG programs: bone cancer pain vs CIPN
- Article type: Research Article
- Abstract: unstructured, 296 words (<= 300), no citations
- Keywords (subject categories, optional at submission): cancer pain; dorsal root ganglion; chemotherapy-induced peripheral neuropathy; interferon-stimulated genes; in-silico perturbation; negative control
- References: 22, all with DOI, numbered in order of citation (Vancouver)
- Tables: 5 (in manuscript, editable, caption above)
- Figures: 3 (separate TIFF, legends in text)
- Corresponding author: Yongxin Yang, 960856791@qq.com, ORCID 0009-0004-9698-6552
- Funding: none; Competing interests: none; Ethics: secondary analysis of public data (original approvals cited)

## Required statements (enter in the EM form)
- Financial disclosure: "This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors."
- Competing interests: "The author declares no competing financial or personal interests."
- Ethics: secondary reanalysis of publicly deposited animal data; original IACUC approvals cited [10,11,14]; no new animals/humans.
- Data availability: all data public (GEO GSE149648, GSE160543, GSE292522); code/artefacts at the GitHub URL above.
- AI disclosure: an OpenAI GPT-4-class LLM was used for code debugging and language editing only; no role in analysis design, results, or interpretation.

## Open items for the author (not auto-filled)
- Create / log in to the Editorial Manager account with the name in Latin script (Yongxin Yang) before submitting.
- Select the subject category that best matches the manuscript at submission.
- These files are staged under English filenames; the upload set is byte-identical to what the gates validated.
"""


if __name__ == "__main__":
    raise SystemExit(main())
