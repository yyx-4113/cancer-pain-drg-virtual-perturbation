import re
P = r"D:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation/docs/manuscript/manuscript_final.md"
t = open(P, encoding="utf-8").read()

head, rest = t.split("## Tables", 1)
tables_sec, figures_sec = rest.split("## Figures", 1)

# --- reorder table blocks: old [1,2,3,4,5] -> [1,4,2,3,5] ---
parts = tables_sec.split("**Table ")
leading = parts[0]
blocks = ["**Table " + p for p in parts[1:]]
assert len(blocks) == 5, len(blocks)
order = [0, 3, 1, 2, 4]            # indices into old blocks
new_blocks = []
for k, idx in enumerate(order, start=1):
    b = re.sub(r"^\*\*Table \d+\.", "**Table %d." % k, blocks[idx])
    new_blocks.append(b)
tables_new = leading + "".join(new_blocks)

# --- remap in-text citations in the body ---
mapping = {2: 3, 3: 4, 4: 2}       # old -> new

def remap_single(m):
    n = int(m.group(1))
    return "Table %d" % mapping.get(n, n)

head = re.sub(r"\bTables (\d)[–-](\d)\b",
              lambda m: "Tables %d–%d" % (mapping.get(int(m.group(1)), int(m.group(1))),
                                          mapping.get(int(m.group(2)), int(m.group(2)))), head)
head = re.sub(r"\bTable (\d)\b", remap_single, head)

out = head + "## Tables" + tables_new + "## Figures" + figures_sec
open(P, "w", encoding="utf-8").write(out)

# report
print("in-text Table refs now:")
for m in re.finditer(r"\bTables? \d[^\n]{0,25}", head):
    print("  ...", m.group(0)[:40])
print("table headers in new order:")
for m in re.finditer(r"^\*\*Table \d+\..*$", open(P, encoding="utf-8").read(), re.M):
    print("  ", m.group(0)[:70])
