"""Static reference data for the site: STRING network, pathway libraries, AE modules.

  docs/data/network.json   STRING v12 human, combined score >= 900 (as in final27):
                           {"genes": [...], "edges": [i0, j0, i1, j1, ...]}
  docs/data/pathways.json  MSigDB v2026.1 Reactome + GO BP, restricted to network genes
                           (KEGG is left out: its licence does not allow public redistribution)
  docs/data/ae.json        the 29 AE modules (core genes with direction, extended genes)
"""
from __future__ import annotations

import gzip
import sys

import pandas as pd

from common import FINAL27, GMT_DIR, STRING_CACHE, STRING_MIN_SCORE, clean_gene, string_symbols, write_json


def network():
    sym = string_symbols()
    edges = set()
    with gzip.open(STRING_CACHE / "9606.protein.links.v12.0.txt.gz", "rt", encoding="utf-8") as fh:
        next(fh, None)
        for line in fh:
            a, b, score = line.split()
            if int(score) < STRING_MIN_SCORE:
                continue
            ga, gb = sym.get(a), sym.get(b)
            if ga and gb and ga != gb:
                edges.add((ga, gb) if ga < gb else (gb, ga))
    genes = sorted({g for e in edges for g in e})
    idx = {g: i for i, g in enumerate(genes)}
    flat = [x for a, b in sorted(edges) for x in (idx[a], idx[b])]
    print(f"[network] {len(genes)} genes, {len(edges)} edges")
    write_json("network.json", {"source": f"STRING v12.0 Homo sapiens, combined score >= {STRING_MIN_SCORE}",
                                "genes": genes, "edges": flat})
    return idx


def pathways(idx: dict[str, int]):
    libs = {"Reactome": "c2.cp.reactome.v2026.1.Hs.symbols.gmt", "GO_BP": "c5.go.bp.v2026.1.Hs.symbols.gmt"}
    out = {"source": "MSigDB v2026.1.Hs (CC BY 4.0): Reactome (C2:CP:REACTOME), GO Biological Process (C5:GO:BP)",
           "libraries": {}}
    for name, fn in libs.items():
        terms = []
        for line in (GMT_DIR / fn).read_text(encoding="utf-8").splitlines():
            parts = line.split("\t")
            members = sorted({idx[g] for g in map(clean_gene, parts[2:]) if g in idx})
            if len(members) >= 3:
                pretty = parts[0].split("_", 1)[1].replace("_", " ").capitalize() if "_" in parts[0] else parts[0]
                terms.append([pretty, members])
        out["libraries"][name] = terms
        print(f"[pathways] {name}: {len(terms)} terms")
    write_json("pathways.json", out)


def ae_modules():
    m = pd.read_csv(FINAL27 / "ae_modules.csv")
    g = pd.read_csv(FINAL27 / "ae_module_genes.csv")
    mods = []
    for _, r in m.iterrows():
        sub = g[g["AE_Module_ID"] == r["AE_Module_ID"]]
        core = sub[sub["模块层级"] == "核心"]
        ext = sub[sub["模块层级"] != "核心"]
        mods.append({
            "id": r["AE_Module_ID"], "cn": r["中文名"], "en": r["英文名"], "short": r["图中简称"], "group": r["系统分组"],
            "synonyms": [s.strip() for s in str(r["同义词"]).split(";") if s.strip()],
            "core": [[clean_gene(x), d] for x, d in zip(core["Gene"], core["基因方向"]) if clean_gene(x)],
            "extended": sorted({clean_gene(x) for x in ext["Gene"] if clean_gene(x)}),
        })
    print(f"[ae] {len(mods)} modules")
    write_json("ae.json", {"source": "final27 AE modules: curated mechanism genes (core) + HPO phenotype genes (extended)",
                           "modules": mods})


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    idx = network()
    pathways(idx)
    ae_modules()
