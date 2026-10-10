#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Verify the remote v1.0.0 tag points to the Round 9 commit whose tree (minus
.github/workflows, which the API cannot write) equals the local ca44643 tree."""
import os, subprocess, json, urllib.request, urllib.error

CWD = r"D:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation"
OWNER, REPO = "yyx-4113", "cancer-pain-drg-virtual-perturbation"

p = subprocess.run(["git", "credential", "fill"], cwd=CWD,
                   input="protocol=https\nhost=github.com\n", capture_output=True, text=True)
token = [l.split("=", 1)[1].strip() for l in p.stdout.splitlines()
         if l.lower().startswith("password=")][0]


def api(path):
    req = urllib.request.Request("https://api.github.com" + path)
    req.add_header("Authorization", "Bearer " + token)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", "yyx-verify")
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


ref = api(f"/repos/{OWNER}/{REPO}/git/refs/tags/v1.0.0")
tag_sha = ref["object"]["sha"]
print("remote tag object:", tag_sha[:12], "| type:", ref["object"]["type"])
tagobj = api(f"/repos/{OWNER}/{REPO}/git/tags/{tag_sha}")
commit_sha = tagobj["object"]["sha"]
print("tagged commit    :", commit_sha[:12])
print("tag message      :", tagobj["message"].strip())
commit = api(f"/repos/{OWNER}/{REPO}/git/commits/{commit_sha}")
remote_tree = commit["tree"]["sha"]
print("remote tree      :", remote_tree[:12])

# local ca44643 tree minus .github/workflows
idx = os.path.join(CWD, ".git", "verify_idx")
env = dict(os.environ, GIT_INDEX_FILE=idx)
subprocess.run(["git", "read-tree", "HEAD^{tree}"], cwd=CWD, env=env, capture_output=True)
subprocess.run(["git", "rm", "--cached", "-r", "-q", ".github/workflows"], cwd=CWD, env=env, capture_output=True)
local_norm = subprocess.run(["git", "write-tree"], cwd=CWD, env=env, capture_output=True).stdout.decode().strip()
try:
    os.remove(idx)
except OSError:
    pass
print("local norm tree  :", local_norm[:12])
print("RESULT:", "MATCH (release == submitted manuscript)" if local_norm == remote_tree
      else "MISMATCH (re-run sync)")
