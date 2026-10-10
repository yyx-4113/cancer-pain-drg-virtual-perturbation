import re, json
import numpy as np, pandas as pd
ROOT = r"D:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation"

print("="*70)
print("A) ERYTHROCYTE COMPARTMENT group counts (script37 / celltype artefacts)")
ct = pd.read_csv(ROOT + "/results/tables/sc_GSE292522_celltype_DE.csv")
print("  celltype_DE cols:", list(ct.columns))
print(ct.to_string())
try:
    loo = json.load(open(ROOT + "/results/tables/sc_GSE292522_compartment_loo.json"))
    print("  compartment_loo keys:", json.dumps(loo, indent=1)[:1500])
except Exception as e:
    print("  loo err", e)

print("="*70)
print("B) myeloid fraction range from composition")
comp = pd.read_csv(ROOT + "/results/tables/sc_GSE292522_composition.csv")
my = comp[["sample","Myeloid","group"]]
print(my.to_string())
print("  min %.5f max %.5f  (as pct: min %.3f%% max %.3f%%)" %
      (my["Myeloid"].min(), my["Myeloid"].max(), 100*my["Myeloid"].min(), 100*my["Myeloid"].max()))

print("="*70)
print("C) mean percentile of targets within all 132")
b = pd.read_csv(ROOT + "/results/tables/etiology_bootstrap_differences.csv")
uni = b[(b["network_scope"]=="per_contrast") & (b["perturbation"]=="knockout")].copy()
from scipy.stats import percentileofscore
z = uni["z_vs_randset"].abs().values
tgt_pct = [percentileofscore(z, abs(v), kind='rank') for v in uni[uni["kind"]=="target"]["z_vs_randset"]]
print("  mean target percentile (rank) = %.2f ; median = %.2f" % (np.mean(tgt_pct), np.median(tgt_pct)))

print("="*70)
print("D) TEXT CHECKS in manuscript")
txt = open(ROOT + "/docs/manuscript/manuscript_final.md", encoding="utf-8").read()
for needle in ["Table 4", "Table 5", "Tables 4", "Tables 5", "(Table 4", "(Table 5",
               "55.7", "44.7", "two tumour and four control", "two control and four tumour",
               "CPM", "per-10k", "per 10k", "10,000", "1e4", "×10"]:
    n = txt.count(needle)
    print("  %-32s -> %d" % (needle, n))

print("="*70)
print("E) reference citation order")
# find in-text [n] order of first appearance after References? just count refs
refs = re.findall(r"^(\d+)\. ", txt, re.M)
print("  ref numbers present:", refs)
# first appearance order of bracket citations in body (before References)
body = txt.split("## References")[0]
order = []
for m in re.finditer(r"\[(\d+(?:,\d+)*(?:–\d+)?|(?:\d+)[,–]\d+)\]", body):
    for tok in re.split(r"[,–]", m.group(1)):
        tok = tok.strip()
        if tok.isdigit() and tok not in order:
            order.append(tok)
print("  first-appearance order:", order)
