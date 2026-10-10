"""Look up DOIs for the manuscript references lacking them, via Crossref."""
import json, urllib.request, urllib.parse, time

TITLES = {
    "1": "Management of cancer pain in adult patients ESMO Clinical Practice Guidelines",
    "2": "Adult cancer pain Version 3.2019 NCCN",
    "3": "Bone cancer pain causes consequences and therapeutic opportunities Mantyh",
    "5": "Chemotherapy-induced peripheral neuropathy a current review Staff Grisold",
    "6": "Prevention and management of chemotherapy-induced peripheral neuropathy in survivors of adult cancers ASCO guideline update",
    "7": "scGen predicts single-cell perturbation responses",
    "8": "Predicting transcriptional outcomes of novel multigene perturbations with GEARS",
    "12": "Linear models and empirical Bayes methods for assessing differential expression in microarray experiments",
    "13": "Controlling the false discovery rate a practical and powerful approach to multiple testing",
    "14": "Bone cancer pain from model to mechanism to therapy",
}

for k, t in TITLES.items():
    try:
        u = ("https://api.crossref.org/works?rows=1&select=DOI,title,container-title,issued"
             "&query.bibliographic=" + urllib.parse.quote(t))
        req = urllib.request.Request(u, headers={"User-Agent": "ref-check/1.0 (mailto:960856791@qq.com)"})
        with urllib.request.urlopen(req, timeout=20) as r:
            d = json.load(r)
        it = d["message"]["items"][0]
        print(k, "|", (it.get("title") or [""])[0][:60], "|",
              it.get("DOI"), "|", (it.get("container-title") or [""])[0][:30])
    except Exception as e:
        print(k, "ERR", e)
    time.sleep(0.6)
