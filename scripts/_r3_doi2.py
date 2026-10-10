"""Targeted Crossref DOI lookups with journal+year filters."""
import json, urllib.request, urllib.parse, time

Q = {
    "1_ESMO": ("Annals of Oncology", 2018, "Management of cancer pain in adult patients ESMO Clinical Practice Guidelines"),
    "2_NCCN": ("Journal of the National Comprehensive Cancer Network", 2019, "Adult Cancer Pain Version 3.2019"),
    "3_Mantyh": ("Pain", 2013, "Bone cancer pain causes consequences therapeutic opportunities"),
    "6_Loprinzi": ("Journal of Clinical Oncology", 2020, "Prevention and management of chemotherapy-induced peripheral neuropathy survivors adult cancers ASCO"),
    "7_scGen": ("Nature Methods", 2019, "scGen predicts single-cell perturbation responses"),
    "14_Luger": ("Journal of Pain and Symptom Management", 2005, "Bone cancer pain from model to mechanism to therapy"),
}

for k, (jrn, yr, t) in Q.items():
    try:
        u = ("https://api.crossref.org/works?rows=3&select=DOI,title,container-title,issued"
             "&query.bibliographic=" + urllib.parse.quote(t) +
             "&query.container-title=" + urllib.parse.quote(jrn))
        req = urllib.request.Request(u, headers={"User-Agent": "ref-check/1.0 (mailto:960856791@qq.com)"})
        with urllib.request.urlopen(req, timeout=20) as r:
            d = json.load(r)
        for it in d["message"]["items"]:
            y = (it.get("issued", {}).get("date-parts", [[None]])[0][0])
            ct = (it.get("container-title") or [""])[0]
            print(k, "|", y, "|", ct[:35], "|", it.get("DOI"), "|", (it.get("title") or [""])[0][:55])
        print("-" * 80)
    except Exception as e:
        print(k, "ERR", e)
    time.sleep(0.6)
