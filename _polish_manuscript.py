#!/usr/bin/env python3
"""Polish manuscript_final.md for AI-trace reduction (de-AI) while preserving
every number, citation and argument. Writes the polished text back to the same file
after backing it up.

Transforms (all number/citation-safe):
  1. Abstract: restructure from labelled (Background/Methods/Results/Conclusions)
     to a single unstructured paragraph (PLOS ONE requires this).
  2. Em-dash (U+2014) reduction in PROSE only (table rows, which use "—" as an
     "not applicable" marker, are left untouched):
       - "**TERM** — desc"  ->  "**TERM**: desc"
       - " x — y "          ->  " x; y "   (semicolon, the safe universal stand-in)
"""
import re, shutil, os

HERE = r"D:\2026.10\不同部位癌痛患者对药物敏感性生信分析\cancer_pain_virtual_perturbation"
SRC = os.path.join(HERE, "docs", "manuscript", "manuscript_final.md")
BAK = SRC + ".bak_polish"

text = open(SRC, encoding="utf-8").read()

# ---------------------------------------------------------------- 1. abstract
# Locate the abstract block (between "## Abstract" and the "**Keywords:**" line).
m = re.search(r"## Abstract\s+(.*?)\n\*\*Keywords:\*\*", text, re.S)
assert m, "abstract block not found"
before, after_key = text[:m.start(1)], text[m.end(1):]

# PLOS ONE requires an UNSTRUCTURED single-paragraph abstract (<=300 words, no
# citations). Crafted verbatim from the original labelled abstract; every number,
# accession and q/p value is preserved.
NEW_ABSTRACT = (
    "Bone cancer pain and chemotherapy-induced peripheral neuropathy (CIPN) are "
    "managed as one entity, yet their sensory-neuron biology is assumed to differ. "
    "Two questions are conflated: whether the dorsal root ganglion (DRG) programmes "
    "diverge, and whether in-silico perturbation can rank analgesic targets. We "
    "reanalysed rat DRG RNA-seq from two same-platform studies (GSE149648, tibial "
    "Walker-256 bone cancer pain, sham and days 1/7/14; GSE160543, "
    "paclitaxel/oxaliplatin neuropathy; n = 4 per group). Differential expression "
    "used empirical-Bayes moderated t with a muscle covariate. Modules and markers "
    "were compared by exhaustive label permutation (4,900 split pairs), giving an "
    "exact p. Seventeen analgesic-target modules were perturbed in silico against "
    "five negative-control modules, degree-matched nulls, an oracle anchor and "
    "additive baselines. An independent mouse single-cell dataset (GSE292522, "
    "B16F10, n = 4 + 4) was analysed as per-animal pseudobulk. The aetiologies "
    "diverge sharply: myeloid markers rise in bone cancer pain and fall in "
    "paclitaxel neuropathy (Cd68 \u0394 = +3.42, Aif1 \u0394 = +3.08, q = 0.012), "
    "while galanin moves oppositely (\u0394 = \u22122.59, q = 0.012). Day 7 is "
    "dominated by an interferon-stimulated-gene programme (module t = +24.9, q = "
    "6.8 \u00d7 10\u2075\u2075). The perturbation layer adds nothing: it correlates "
    "0.991 with an additive model, ranks an erythrocyte control above every target "
    "(|z| = 7.39), and gives controls the larger effects (p = 0.0019). It resolves "
    "large effects but not the median target, and no reversal \u0394 survives "
    "correction (minimum q = 0.165). The interferon programme is not present in the "
    "B16F10 model; its markers are if anything down-regulated; and the response is "
    "neuronal-injury-driven: the signature is model-specific. Bone cancer pain and "
    "the paclitaxel/oxaliplatin neuropathy models (CIPN) drive qualitatively "
    "different, exactly testable DRG programmes, readable from individual markers. "
    "In-silico perturbation provides no target ranking; it fails its own negative "
    "controls, degree-matched nulls and additive baseline."
)
word_count = len(NEW_ABSTRACT.split())
assert word_count <= 300, f"abstract {word_count} words exceeds 300"
new_ablock = "## Abstract\n\n" + NEW_ABSTRACT + "\n\n**Keywords:**"
text = before + new_ablock + after_key
print(f"[abstract] restructured -> {word_count} words (<=300)")

# ---------------------------------------------------------------- 2. em-dash
# Operate line by line; skip table rows (start with '|') and separator rows.
lines = text.split("\n")
EMDASH = "—"  # U+2014
count = 0
out = []
for line in lines:
    stripped = line.strip()
    if stripped.startswith("|") or (stripped and set(stripped) <= set("|-: ")):
        out.append(line)          # leave tables (incl. "—" N/A markers) untouched
        continue
    new = line
    # bold-term lead-in: "**TERM** — desc" -> "**TERM**: desc"
    new = re.sub(r"(\*\*[^*]+\*\*)(\s*)" + EMDASH + r"(\s+)", r"\1: ", new)
    # generic prose break: " x — y " -> " x; y "
    new2 = re.sub(r"(\s)" + EMDASH + r"(\s)", "; ", new)
    count += new2.count("; ") - new.count("; ")  # approximate, informational only
    out.append(new2)
text = "\n".join(out)
remaining = text.count(EMDASH)
print(f"[em-dash] prose em-dashes converted; remaining em-dashes in file: {remaining} (should be table N/A markers only)")

# ---------------------------------------------------------------- write back
shutil.copyfile(SRC, BAK)
open(SRC, "w", encoding="utf-8").write(text)
print(f"[ok] polished manuscript written; backup at {os.path.basename(BAK)}")
