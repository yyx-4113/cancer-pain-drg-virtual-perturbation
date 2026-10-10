import re
P = r"D:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation/docs/manuscript/manuscript_final.md"
t = open(P, encoding="utf-8").read()

# revert the intro insertion of [14] (Luger will be renumbered to last)
t = t.replace("inside the DRG [3,4,14].", "inside the DRG [3,4].")

body, refs = t.split("## References")
tables = refs.split("## Tables")
refblock = tables[0]
rest = "## Tables" + tables[1]

# parse reference entries
entries = {}
for m in re.finditer(r"^(\d+)\. (.+)$", refblock, re.M):
    entries[int(m.group(1))] = m.group(2)

# desired new order: 1..13 unchanged, then old15, old16, old14
new_order = list(range(1, 14)) + [15, 16, 14]
mapping = {old: new for new, old in enumerate(new_order, start=1)}
# -> {14:16, 15:14, 16:15}

# remap in-text citations in body
def remap(m):
    nums = [int(x) for x in m.group(1).split(",")]
    return "[" + ",".join(str(mapping[n]) for n in nums) + "]"
body = re.sub(r"\[(\d+(?:,\d+)*)\]", remap, body)

# rebuild reference list
lines = []
for new_n, old_n in enumerate(new_order, start=1):
    lines.append(f"{new_n}. {entries[old_n]}")
newrefblock = "\n" + "\n".join(lines) + "\n\n"

out = body + "## References" + newrefblock + rest
open(P, "w", encoding="utf-8").write(out)
print("mapping:", mapping)
print("rewritten")
