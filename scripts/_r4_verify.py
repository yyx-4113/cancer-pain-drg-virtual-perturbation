import re, json, math
import numpy as np, pandas as pd
ROOT = r"D:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation"
MS = ROOT + "/docs/manuscript/manuscript_final.md"
txt = open(MS, encoding="utf-8").read()

print("="*70)
print("1) ABSTRACT WORD COUNT")
# abstract between '## Abstract' and '**Keywords'
m = re.search(r"## Abstract\s*(.*?)\*\*Keywords", txt, re.S)
abs_body = m.group(1)
# strip headings labels
body = abs_body
label_tokens = re.findall(r"\*\*(Background|Methods|Results|Conclusions)\.\*\*", body)
clean = re.sub(r"\*\*(Background|Methods|Results|Conclusions)\.\*\*", " ", body)
clean = re.sub(r"[*\u00b9\u2082\u2083\u2084\u2085\u207b\u2070-\u2079]+", " ", clean)
words_body = len([w for w in clean.split() if any(c.isalnum() for c in w)])
print("  body words (labels stripped):", words_body, " labels:", len(label_tokens))
with_labels = words_body + len(label_tokens)
print("  with labels:", with_labels)
# alt: count all alnum tokens including labels
clean2 = re.sub(r"[*\u00b9\u2082\u2083\u2084\u2085\u207b\u2070-\u2079]+", " ", re.sub(r"\*\*(Background|Methods|Results|Conclusions)\.\*\*", r"\1.", body))
words_all = len([w for w in clean2.split() if any(c.isalnum() for c in w)])
print("  all tokens incl label words:", words_all)

print("="*70)
print("2) TARGET MEDIAN |z| PERCENTILE (§3.6 '55.7th percentile')")
b = pd.read_csv(ROOT + "/results/tables/etiology_bootstrap_differences.csv")
print("  cols:", [c for c in b.columns])
uni = b[(b["network_scope"]=="per_contrast") & (b["perturbation"]=="knockout")].copy()
print("  dedup KO universe n =", len(uni), " kinds:", uni["kind"].value_counts().to_dict())
z = uni["z_vs_randset"].abs().values
tgt = uni[uni["kind"]=="target"]["z_vs_randset"].abs().values
allz = np.sort(z)
med_t = np.median(tgt)
# percentile of target median within all
below = (allz < med_t).sum(); eq = (allz == med_t).sum()
print("  target median |z| = %.4f ; all median = %.4f" % (med_t, np.median(z)))
print("  rank-below pct = %.1f ; pct_gt (scipy-style) = %.1f" %
      (100*below/len(allz), 100*(below+0.5*eq)/len(allz)))
# standard percentileofscore 'rank' = (below + 0.5*eq)/n
from scipy.stats import percentileofscore
print("  scipy percentileofscore rank =", round(percentileofscore(z, med_t, kind='rank'),1))
print("  scipy percentileofscore mean =", round(percentileofscore(z, med_t, kind='mean'),1))

print("="*70)
print("3) ERYTHROCYTE COMPARTMENT 2v4 GROUP COUNTS")
comp = pd.read_csv(ROOT + "/results/tables/sc_GSE292522_composition.csv")
print("  composition cols:", list(comp.columns))
print(comp.to_string())

print("="*70)
print("4) MYELOID FRACTION RANGE")
if "cell_type" in comp.columns:
    my = comp[comp["cell_type"].str.lower().str.contains("myeloid")]
    print(my.to_string())
