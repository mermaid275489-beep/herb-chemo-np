"""Map TCMSP target names to human gene symbols.

TCMSP gives only protein names. Order of resolution:
  1. the target-name -> gene table already used in the final27 analysis (tcmspswiss_targets.csv)
  2. UniProt REST, reviewed human entries whose recommended or alternative name matches exactly
     (case-insensitive). Anything else, including non-human targets such as bacterial enzymes,
     stays unmapped, as in final27.
Results are cached in raw/target_map.json, so reruns only query new names.

  python build/map_targets.py
"""
from __future__ import annotations

import json
import re
import sys
import time

import pandas as pd
import requests

from common import LEGACY_TARGET_MAP, RAW, clean_gene

CACHE = RAW / "target_map.json"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"


def tcmsp_target_names() -> set[str]:
    names = set()
    for f in (RAW / "tcmsp" / "herbs").glob("*.json"):
        for e in json.loads(f.read_text(encoding="utf-8"))["edges"]:
            if e.get("target_name"):
                names.add(e["target_name"].strip())
    return names


def names_of(entry: dict) -> list[str]:
    d = entry.get("proteinDescription") or {}
    out = []
    for block in [d.get("recommendedName")] + (d.get("alternativeNames") or []):
        if block and block.get("fullName"):
            out.append(block["fullName"]["value"])
    for sub in (d.get("contains") or []) + (d.get("includes") or []):
        if sub.get("recommendedName", {}).get("fullName"):
            out.append(sub["recommendedName"]["fullName"]["value"])
    return out


def tokens(s: str) -> tuple[str, ...]:
    """Name as a sorted multiset of alphanumeric tokens, so that punctuation ([NAD+] vs [NAD(+)],
    'epsilon 1' vs 'epsilon-1') and word order ('A2a receptor' vs 'receptor A2a') do not matter."""
    return tuple(sorted(re.findall(r"[a-z0-9]+", s.casefold())))


def query_uniprot(q: str) -> list[dict] | None:
    for attempt in range(4):
        try:
            r = requests.get(UNIPROT, params={"query": f"({q}) AND (organism_id:9606) AND (reviewed:true)",
                                              "fields": "accession,gene_primary,protein_name", "format": "json",
                                              "size": 25}, timeout=60)
            r.raise_for_status()
            return r.json().get("results", [])
        except requests.RequestException as e:
            print(f"  retry {attempt + 1}: {e}", flush=True)
            time.sleep(3 * (attempt + 1))
    return None


def pick(name: str, entries: list[dict], match) -> dict | None:
    """Entries whose names match; prefer a match on the recommended name when several genes qualify."""
    hits = []
    for entry in entries:
        genes = entry.get("genes") or []
        gene = clean_gene(genes[0].get("geneName", {}).get("value")) if genes else ""
        names = names_of(entry)
        if gene and any(match(n) for n in names):
            hits.append((gene, entry["primaryAccession"], bool(names) and match(names[0])))
    genes = sorted({h[0] for h in hits})
    if len(genes) > 1:
        rec = sorted({h[0] for h in hits if h[2]})
        if len(rec) == 1:
            genes = rec
        else:
            return {"gene": "", "status": "uniprot_ambiguous", "candidates": genes}
    if len(genes) == 1:
        return {"gene": genes[0], "accession": next(h[1] for h in hits if h[0] == genes[0])}
    return None


def uniprot(name: str) -> dict:
    entries = query_uniprot(f'protein_name:"{name}"')
    if entries is None:
        return {"gene": "", "status": "query_failed"}
    want = name.casefold()
    hit = pick(name, entries, lambda n: n.casefold() == want)
    if hit and hit["gene"]:
        return {**hit, "status": "uniprot_exact"}
    # second pass: token-set match, query by all tokens
    tk = tokens(name)
    if tk:
        entries = query_uniprot(" AND ".join(f"protein_name:{t}" for t in dict.fromkeys(re.findall(r"[a-z0-9]+", name.casefold()))))
        hit2 = pick(name, entries or [], lambda n: tokens(n) == tk)
        if hit2 and hit2["gene"]:
            return {**hit2, "status": "uniprot_token_match"}
        if hit2:
            return hit2
    return hit or {"gene": "", "status": "unmapped_or_nonhuman"}


def main():
    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    legacy = pd.read_csv(LEGACY_TARGET_MAP).dropna(subset=["Target_Name", "Gene"])
    legacy_map = {}
    for t, g in zip(legacy["Target_Name"], legacy["Gene"]):
        if clean_gene(g):
            legacy_map.setdefault(t.strip(), set()).add(clean_gene(g))
    names = sorted(tcmsp_target_names())
    retry = {"unmapped_or_nonhuman", "uniprot_ambiguous", "query_failed"}
    todo = [n for n in names if n not in cache or (cache[n]["status"] in retry and not cache[n].get("v2"))]
    print(f"[targets] {len(names)} TCMSP target names, {len(todo)} new", flush=True)
    for i, n in enumerate(todo, 1):
        if n in legacy_map and len(legacy_map[n]) == 1:
            cache[n] = {"gene": next(iter(legacy_map[n])), "status": "final27_table"}
        else:
            cache[n] = {**uniprot(n), "v2": True}
            time.sleep(0.15)
        if i % 100 == 0:
            CACHE.write_text(json.dumps(cache, ensure_ascii=False, indent=0), encoding="utf-8")
            print(f"  {i}/{len(todo)}", flush=True)
    CACHE.write_text(json.dumps(cache, ensure_ascii=False, indent=0), encoding="utf-8")
    status = pd.Series([cache[n]["status"] for n in names]).value_counts()
    print(status.to_string())


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
