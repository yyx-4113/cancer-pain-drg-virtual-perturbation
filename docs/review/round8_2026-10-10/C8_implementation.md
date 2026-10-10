# Reviewer C8 — Implementation / Provenance Audit (Round 8, 2026-10-10)

**Manuscript:** `docs/manuscript/manuscript_final.md`
**Repo:** `D:\2026.10\不同部位癌痛患者对药物敏感性生信分析\cancer_pain_virtual_perturbation`
**Remote (cited):** `https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation`, release tag `v1.0.0`

I treated this as a first submission and verified every judgement independently from the text, the deposited artefacts (`results/tables/`), and the git/GitHub release state. I did **not** read any prior review, the author-verification statement, the SOP, or any `manuscript_v0*.md`.

---

## Findings

### Finding 1 — Release contains the submitted manuscript (verified), but local checkout ≠ released commit
【Problem】 The local repository's `HEAD` and local `v1.0.0` tag both resolve to commit `e0941b7`, whereas the **remote** `v1.0.0` tag (and remote `main`) resolve to a *different* commit `88ddb527`. The checkout I audited is not the released commit.

【Evidence】
- `git log --oneline -1` → `e0941b7 (HEAD -> main, tag: v1.0.0) Round 7 revisions...`
- `git rev-list -n 1 v1.0.0` → `e0941b7230abc2717bbe6b7fe6708c1e56de24ae` (= local HEAD)
- GitHub API `git/refs/tags/v1.0.0` → annotated tag object `c254d7a2…` → commit `88ddb5278adf684b3192ebed0c403b999cd6e480`
- GitHub API `commits/main` → `88ddb527…` (so **remote tag == remote main**)
- Decisive check: manuscript blob SHA-1 at the remote tag (`contents/.../manuscript_final.md?ref=v1.0.0`) = `011a0df3596a8ddf270565fb5d94ae1e33ee69d6`; local `git hash-object docs/manuscript/manuscript_final.md` = `011a0df3596a8ddf270565fb5d94ae1e33ee69d6`. **Identical.** (The raw-byte sha256 differs only because the working copy is CRLF, 66,250 B, while the repo/release stores LF, 65,898 B — a line-ending normalization Git applies on hash; content is identical.)

【Why it matters】 The core mandate — *"the cited release must contain the submitted manuscript"* — is **satisfied**: the manuscript in release `v1.0.0` (`88ddb527`) is byte-for-byte identical (after normalization) to the submitted file. However, the local analysis scripts/tables I verified the numbers against belong to commit `e0941b7`, which is **not** what a reader reaches via the citation. I could not enumerate the remote tree (GitHub API rate-limited me at 403 after the earlier calls), so I cannot independently confirm that the released commit's `results/tables/` reproduce the manuscript numbers. The manuscript text itself matches, so this is a reproducibility gap rather than a content defect.

【Specific fix】 (a) Author should run `scripts/98_check_manuscript_numbers.py` and `scripts/99_unit_tests.py` against a fresh checkout of the released `v1.0.0` commit (`88ddb527`) and confirm both pass; (b) align the local `v1.0.0` tag to the released commit, or explicitly document that the review snapshot is a working/earlier commit. **Tier 2 — not a manuscript-content blocker** because the manuscript blob is confirmed present in the cited release.

---

### Finding 2 — "Python 3.13" asserted in prose but not pinned anywhere
【Problem】 Manuscript §2.9 and the README claim **Python 3.13**, but `requirements.txt` pins only the five numeric libraries, and `CITATION.cff`'s `version: 1.0.0` is the *release* version, not the interpreter. The interpreter version is therefore asserted but not enforced.

【Evidence】
- `manuscript_final.md` §2.9: *"Python 3.13 with numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0 and matplotlib 3.11.1 (`requirements.txt`)."*
- `README.md` line 183: *"Python 3.13 (see `requirements.txt` for the exact numeric stack)."*
- `requirements.txt` lists only `numpy==2.4.6 / scipy==1.18.1 / pandas==2.3.3 / statsmodels==0.15.0 / matplotlib==3.11.1` (no `python` pin, no `requires-python`).
- `CITATION.cff`: `version: 1.0.0` (release), no interpreter field.

【Why it matters】 Minor reproducibility gap: a future environment with a different Python minor version could behave differently. Not a number-traceability defect.

【Specific fix】 Add a `requires-python = ">=3.13,<3.14"` (or pin `python==3.13.x`) to `pyproject`/`requirements.txt`, or soften the prose to "Python 3.13 used". **Tier 3 — trivial.**

---

### Finding 3 — No defect in the numeric/DOI/cross-file layers (positive)
【Problem】 None. (Recorded to document coverage; see "Stands up" and "What I actually checked".)

【Evidence】 See below: 291/291 number-gate checks pass; unit tests pass; ≥10 independent table-cell recomputations reproduce the printed values; the three added DOIs resolve to the exact cited papers; library versions are consistent across manuscript/README/requirements; accessions and n-values agree between README and manuscript.

【Why it matters】 The implementation/provenance surface is largely sound; the only open item is the local-vs-released commit gap in Finding 1.

【Specific fix】 None required for acceptance.

---

## Stands up (verified independently)

1. **Number-gate integrity.** `scripts/98_check_manuscript_numbers.py` reports **`291 checks, 0 failed`** and *"ALL MANUSCRIPT NUMBERS VERIFIED AGAINST ARTEFACTS"* (exit 0). It re-derives key quantities at run time (T1 median |z|, T2 binomial p, T5 MDE, T1/T2/T5 regenerated vs JSON) rather than comparing to hardcoded literals, and asserts all 19 references carry a DOI.

2. **Unit-test integrity.** `scripts/99_unit_tests.py` reports **`ALL TESTS PASSED`** (exit 0), covering exactly the bug classes the manuscript claims: no self-loops, knockout/overexpression as exact mirrors, rejection of identifier matrices, D14≠D1, degree-matched draws staying in-bin, exhaustive splits containing the true assignment, per-column batched-cosine normalization.

3. **Independent recomputation of ≥10 cells across Tables 1–5 matches the printed values** (see below — all within reported rounding):
   - **T1 ISG (D7):** mean t = 24.923 → printed `+24.92`; q = 6.759e-05 → `6.8×10⁻⁵`. Glial/immune mean t = 11.412 → `+11.41`, q 0.00267 → `0.0027`. Housekeeping mean t = 0.653 → `+0.65`, q 0.5237 → `0.52`.
   - **T3 Cd68 (D7 vs PTX):** logFC 2.839/−0.579, Δ 3.418 → `+3.42`, q 0.0125 → `0.012`. **Aif1:** 2.192/−0.883, Δ 3.076 → `+3.08`, q 0.012. **Gal:** −0.672/+1.920, Δ −2.592 → `−2.59`, q 0.012. **Ccl2 (vs OXA):** 0.897/−0.493, Δ 1.391 → `+1.39`, q 0.0173 → `0.017`. **Il1b:** 0.800/−0.431, Δ 1.230 → `+1.23`, q 0.0182 → `0.018`. **Gfap:** 1.498 → `+1.50` (Bonferroni ns, q 0.0907).
   - **T4 Glial/immune (D7 vs PTX):** 2.053/−0.643, Δ 2.695 → `+2.70`, q 0.0073. **Cytokine/chemokine:** 1.396/−0.306, Δ 1.702 → `+1.70`, q 0.0073. **TLR:** 0.854/−0.421, Δ 1.275 → `+1.28`, q 0.0073. **Neuropeptides:** −0.515/+0.830, Δ −1.345 → `−1.35`, q 0.0145. **CTRL erythrocyte (vs OXA):** −1.054/+1.429, Δ −2.483 → `−2.48`, q 0.0073.
   - **T5 myeloid fraction:** Con mean 0.1170% (printed `0.12%`), B16F10 0.0517% (`0.05%`); Cohen d (n−2) = 1.267 → `1.27`; MWU against degree-matched null: two-sided p = 0.114286 → `0.114`, one-sided (correct "control > tumour" tail) p = 0.05714 → `0.057`; U = 2/14 (the two U statistics are 2 and 14). All match.

4. **DOI authenticity — all three added references resolve and match exactly:**
   - **[17]** `10.1007/s00011-022-01663-2` → *Inflammation Research* 2023;72(1):117–132, "Activation of the STING pathway induces peripheral sensitization via neuroinflammation in a rat model of bone cancer pain", authors Zhang Y, Wang W, Gong Z, Peng Y, Li X, Zhang Z, Zhang X, You X, Wu J (Crossref) — matches manuscript line 236.
   - **[18]** `10.1073/pnas.2503779122` → *PNAS* 2025;122(34):e2503779122, "Unraveling the neuroimmune mechanisms in cancer-induced bone pain…", authors Wang C, Ji H, Wang H, Chen Z, Zhou L, Yang Y, Jiang Y, Yu G, Jiang L, Tang Z (Crossref) — matches manuscript line 237.
   - **[19]** `10.1016/j.brainres.2007.06.066` → *Brain Research* 2007;1168:46–59, "An evolving cellular pathology occurs in dorsal root ganglia… paclitaxel in the rat", authors Peters CM, Jimenez-Andrade JM, Kuskowski MA, Ghilardi JR, Mantyh PW (Crossref) — matches manuscript line 238.

5. **Cited release contains the submitted manuscript.** Remote `v1.0.0` (`88ddb527` = remote `main`) manuscript blob SHA-1 = `011a0df3…` = local `git hash-object` of `manuscript_final.md`. The release therefore holds the exact submitted file.

6. **Cross-file version/accession consistency.** Library versions `numpy 2.4.6 / scipy 1.18.1 / pandas 2.3.3 / statsmodels 0.15.0 / matplotlib 3.11.1` are identical across manuscript §2.9, `requirements.txt`, and `README.md` (line 114). `CITATION.cff` `version: 1.0.0` / `date-released: 2026-10-10` matches the tag. Accessions (GSE149648, GSE160543, GSE292522; GSE249443 as PoC) and n-values (n=4 per group; n=4+4 sc) agree between README and manuscript.

---

## Questions for the author (non-blocking)

1. Between local commit `e0941b7` and the released `88ddb527`, what changed? The manuscript blob is identical, but please confirm the released `results/tables/` and `scripts/` reproduce the reported numbers (i.e., `98`/`99` pass on a checkout of `88ddb527`).
2. Would you pin the Python interpreter in `requirements.txt`/`CITATION.cff` to remove the only version assertion that is currently unenforced?

---

## What I actually checked (commands, outputs, discrepancies)

**Environment:** PowerShell (per mandate); managed interpreter `C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe`; `$env:PYTHONIOENCODING="utf-8"`.

1. **Number gate** — `python scripts/98_check_manuscript_numbers.py manuscript_final.md` (note: the script joins its own `docs/manuscript/` prefix, so the bare filename must be passed; passing the full path doubled it and raised `FileNotFoundError` — a process error on my side, not a manuscript defect). Correct run tail:
   ```
   [PASS] abstract <= 300 words: True
   [PASS] all references carry a DOI: 19
   [PASS] figure file present (fig1.png): True  (fig2/fig3 same)
   ========================================================================
   291 checks, 0 failed
   ALL MANUSCRIPT NUMBERS VERIFIED AGAINST ARTEFACTS
   ```
   Exit 0.

2. **Unit tests** — `python scripts/99_unit_tests.py` → tail: `T14e …`, `T14f (2,4) shape`, `T14g …`, `T13c …`; then `====== ALL TESTS PASSED`. Exit 0.

3. **Table recomputation** — wrote a standalone script reading `results/tables/module_level_signal.csv`, `etiology_marker_cross.csv`, `etiology_module_logfc_cross.csv`, `sc_GSE292522_composition.csv`, `sc_GSE292522_summary.json`; recomputed every cell listed in Stands-up #3 from raw artefact columns. All matched the manuscript to the printed precision. (Script deleted after use; not left in repo.)

4. **DOIs** — WebFetch on `doi.org` for [17] returned a Springer bot-challenge page and [19] returned empty; the PNAS page did not yield clean metadata. Used the **Crossref REST API** (`api.crossref.org/works/<doi>`) via `curl`, which returned authoritative title/author/journal/year/volume/issue/pages for all three — all matching. [18] additionally confirmed against the live PNAS page title.

5. **Git / release provenance:**
   - `git status -s` → only **untracked** helper files (`_*.py`, `tmp_*`, `deposited_neg_ref.json`, `err.txt`, `docs/review/round8_2026-10-10/`, plus my own temp artefacts); **no modified tracked file** — `manuscript_final.md` is clean. Working tree is clean for tracked content.
   - `git log --oneline -1` and `git rev-list -n 1 v1.0.0` → `e0941b7` (local tag == local HEAD).
   - `git ls-remote …` against github.com **failed** (`Recv failure: Connection was reset` — Git's HTTPS transport to GitHub is blocked in this sandbox), so I used the GitHub REST API instead.
   - API: `git/refs/tags/v1.0.0` → annotated tag `c254d7a2…` → commit `88ddb527…`; `commits/main` → `88ddb527…` (remote tag == remote main).
   - Manuscript blob at remote tag (`contents/.../manuscript_final.md?ref=v1.0.0`) SHA-1 `011a0df3…` == local `git hash-object` `011a0df3…`. **Confirmed identical.**
   - Discrepancy investigated: raw sha256 of the working file (`1f6b10b3…`, 66,250 B, CRLF) ≠ decoded remote content sha256 (`009f57e9…`, 65,898 B, LF). Reconciled as **line-ending normalization only** — the git blob SHA-1 (which hashes `blob <size>\0<bytes>` after Git's clean filter) is identical, so content is identical. Not a defect.
   - Attempted a recursive remote-tree comparison (one `git/trees` call) to scope the `e0941b7`/`88ddb527` divergence, but the GitHub API returned `403 rate limit exceeded`; the manuscript-blob match already settles the decisive question (release contains the manuscript).

6. **Cross-file drift** — `requirements.txt`, `README.md` (line 114), and `manuscript_final.md` §2.9 all carry identical library versions. `CITATION.cff` carries only the release version `1.0.0` (no lib pins). README and manuscript agree on accessions and n-values. Only gap: "Python 3.13" is in prose but unpinned (Finding 2).

---

## VERDICT

**ACCEPTABLE — no remaining implementation/provenance blocker.** The number gate (291/0), the unit-test suite (ALL PASSED), an independent recomputation of ≥10 Table 1–5 cells, all three added DOIs, and the presence of the submitted manuscript in the cited release `v1.0.0` (remote commit `88ddb527`) all hold. One Tier-2 provenance item remains (local checkout `e0941b7` ≠ released commit `88ddb527`; recommend confirming `98`/`99` pass on the released checkout) and one Tier-3 item (pin Python 3.13) — neither blocks acceptance.
