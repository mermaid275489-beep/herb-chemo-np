"""Assemble the per-entity data files the page loads.

  docs/data/herbs.json    TCMSP herbs -> ingredients (OB, DL) -> human target genes
  docs/data/drugs.json    anticancer drugs: PD genes, ADME genes, FAERS signals mapped to AE modules
  docs/data/cancers.json  Open Targets top-300 gene modules for the curated cancer list
  docs/data/kb.json       rule-based PK/PD knowledge base from hdi_evidence_kg (for the mechanism panel)

Run after fetch_tcmsp.py, map_targets.py and fetch_opentargets.py.
"""
from __future__ import annotations

import json
import re
import sys

import pandas as pd

from common import FINAL27, KB, LEXICON, RAW, clean_gene, write_json

# PD core genes curated in the final27 analysis (CHEMO_CORE); these take precedence over
# Open Targets so that results match the paper.
FINAL27_PD = {
    "doxorubicin": ["TOP2A", "TOP2B"], "cyclophosphamide": ["CYP2B6", "CYP2C19", "CYP3A4", "ALDH1A1"],
    "fluorouracil": ["TYMS", "DPYD", "UMPS"], "docetaxel": ["TUBB", "TUBB3"], "paclitaxel": ["TUBB", "TUBB3"],
    "oxaliplatin": ["ERCC1", "ERCC2", "XPA", "SLC31A1", "ATP7A", "ATP7B"],
    "cisplatin": ["ERCC1", "ERCC2", "XPA", "SLC31A1", "ATP7A", "ATP7B"],
    "carboplatin": ["ERCC1", "ERCC2", "XPA", "SLC31A1", "ATP7A", "ATP7B"],
    "leucovorin": ["SLC19A1", "FOLR1", "MTHFR", "MTHFD1"], "capecitabine": ["TYMS", "DPYD", "CES2", "CDA"],
    "gemcitabine": ["RRM1", "RRM2", "DCK", "CDA", "SLC29A1"], "vinorelbine": ["TUBB", "TUBB3"],
    "etoposide": ["TOP2A", "TOP2B"], "pemetrexed": ["TYMS", "DHFR", "GART"],
    "irinotecan": ["TOP1", "UGT1A1", "ABCB1", "ABCG2"], "bevacizumab": ["VEGFA"],
}
# Same-class stand-ins for drugs whose Open Targets record has no protein target.
CLASS_PROXY = {
    "nedaplatin": "cisplatin", "lobaplatin": "cisplatin", "epirubicin": "doxorubicin", "pirarubicin": "doxorubicin",
    "vindesine": "vinorelbine", "ifosfamide": "cyclophosphamide",
}
CURATED_PD = {
    "temozolomide": (["MGMT", "MLH1", "MSH2", "MSH6", "PMS2"], "MGMT 修复与错配修复决定敏感性（文献共识）"),
    "dacarbazine": (["MGMT", "MLH1", "MSH2", "MSH6", "PMS2"], "与替莫唑胺同为 O6-甲基鸟嘌呤烷化剂"),
    "lomustine": (["MGMT"], "O6-烷基鸟嘌呤由 MGMT 修复"),
    "fludarabine": (["POLA1", "RRM1", "RRM2", "DCK"], "抑制 DNA 聚合酶 α 与核糖核苷酸还原酶，经 dCK 活化"),
    "trifluridine": (["TYMS", "TYMP"], "抑制胸苷酸合成酶；替匹嘧啶抑制 TYMP"),
}

# MedDRA preferred-term keywords (FAERS) -> AE module; substring match on lower-cased PT names.
FAERS_TO_AE = {
    "AE01": ["nausea", "vomiting"], "AE02": ["decreased appetite", "anorexia"], "AE03": ["diarrhoea", "diarrhea"],
    "AE04": ["constipation"], "AE05": ["abdominal pain"], "AE06": ["stomatitis", "mucosal inflammation"],
    "AE07": ["oral pain"], "AE08": ["oesophagitis", "gastritis", "mucositis"], "AE09": ["dyspepsia", "abdominal discomfort"],
    "AE10": ["anaemia", "anemia"], "AE11": ["leukopenia", "white blood cell count decreased"],
    "AE12": ["neutropenia", "neutrophil count decreased"], "AE13": ["thrombocytopenia", "platelet count decreased"],
    "AE14": ["hypercalcaemia"], "AE15": ["hypermagnesaemia"],
    "AE16": ["neuropathy", "neurotoxicity", "paraesthesia", "hypoaesthesia"], "AE17": ["dizziness", "vertigo"],
    "AE18": ["headache"], "AE19": ["pain"], "AE20": ["fatigue", "asthenia"], "AE21": ["pyrexia", "fever"],
    "AE22": ["palmar-plantar", "hand-foot"], "AE23": ["hypersensitivity", "anaphyla", "infusion related reaction"],
    "AE24": ["injection site", "infusion site", "extravasation"], "AE25": ["pruritus"], "AE26": ["rash"],
    "AE27": ["hepat", "transaminase", "liver"], "AE28": ["hypophosphataemia"], "AE29": ["dry mouth"],
}
# Rule-KB AE keys -> AE modules
KB_AE_TO_MODULE = {
    "neutropenia": ["AE12", "AE11"], "thrombocytopenia": ["AE13"], "diarrhea": ["AE03"], "nausea_vomiting": ["AE01"],
    "mucositis": ["AE06", "AE08"], "hepatotoxicity": ["AE27"], "peripheral_neuropathy": ["AE16"],
    "hypersensitivity": ["AE23"], "hand_foot": ["AE22"], "rash": ["AE26"], "constipation_ileus": ["AE04"],
}


def herb_aliases() -> dict[str, list[str]]:
    """Chinese name -> English/Latin synonyms from the hdi_evidence_kg lexicon, so that e.g.
    'Astragalus' finds 黄芪 (TCMSP lists it as 'Hedysarum Multijugum Maxim.')."""
    lex = json.loads(LEXICON.read_text(encoding="utf-8"))
    out = {}
    for key, h in lex["herbs"].items():
        ascii_names = [key.replace("_", " ")] + [s for s in h["synonyms"] if s.isascii()]
        for s in h["synonyms"]:
            if re.search(r"[一-鿿]", s):
                out.setdefault(s, []).extend(ascii_names)
    return out


def herbs(target_map: dict) -> dict:
    aliases = herb_aliases()
    gene_idx, genes, comp_idx, compounds, out = {}, [], {}, [], []
    skipped = []
    for f in sorted((RAW / "tcmsp" / "herbs").glob("*.json")):
        d = json.loads(f.read_text(encoding="utf-8"))
        h, ingredients, edges = d["herb"], d["ingredients"], d["edges"]
        targets_by_mol: dict[str, dict[int, int]] = {}
        for e in edges:
            g = target_map.get((e.get("target_name") or "").strip(), {}).get("gene", "")
            if not g:
                continue
            if g not in gene_idx:
                gene_idx[g] = len(genes)
                genes.append(g)
            validated = 1 if str(e.get("validated") or "").strip() not in {"", "0", "None"} else 0
            t = targets_by_mol.setdefault(e["MOL_ID"], {})
            t[gene_idx[g]] = max(t.get(gene_idx[g], 0), validated)
        comps = []
        for ing in ingredients:
            mol = ing["MOL_ID"]
            if mol not in comp_idx:
                ob = float(ing.get("ob") or 0)
                dl = float(ing.get("dl") or 0)
                tg = targets_by_mol.get(mol, {})
                comp_idx[mol] = len(compounds)
                compounds.append([mol, ing["molecule_name"], round(ob, 2), round(dl, 3),
                                  sorted(tg), sorted(g for g, v in tg.items() if v)])
            else:
                # targets can differ by herb page only if TCMSP is inconsistent; union them
                c = compounds[comp_idx[mol]]
                tg = targets_by_mol.get(mol, {})
                c[4] = sorted(set(c[4]) | set(tg))
                c[5] = sorted(set(c[5]) | {g for g, v in tg.items() if v})
            comps.append(comp_idx[mol])
        if not ingredients:
            skipped.append(h["herb_cn_name"])
        out.append({"cn": h["herb_cn_name"], "py": h["herb_pinyin"], "la": h["herb_en_name"],
                    "alias": sorted(set(aliases.get(h["herb_cn_name"], []))),
                    "cat": h.get("child_cn_name") or "", "cat_en": h.get("child_en_name") or "", "c": comps})
    print(f"[herbs] {len(out)} herbs, {len(compounds)} compounds, {len(genes)} mapped genes; "
          f"{len(skipped)} herbs without ingredients")
    return {"source": "TCMSP 2.3 (old.tcmsp-e.com), retrieved 2026-10; target names mapped to human genes via UniProt",
            "ob_min": 30, "dl_min": 0.18, "genes": genes, "compounds": compounds, "herbs": out}


def adme_table() -> dict[str, list[list[str]]]:
    """drug -> [[gene, role, source]] from final27 chemotherapy_adme.csv and the rule KB."""
    out: dict[str, dict[str, list[str]]] = {}
    a = pd.read_csv(FINAL27 / "chemotherapy_adme.csv")
    for drug, gene, role, system in zip(a["Chemotherapy_drug"], a["Gene"], a["Role"], a["System"]):
        g = clean_gene(gene)
        if g:
            out.setdefault(drug.lower(), {}).setdefault(g, [g, f"{role}（{system}）", "final27"])
    kb = json.loads(KB.read_text(encoding="utf-8"))
    alias = {"P-GP": "ABCB1", "DPD": "DPYD", "CES1": "CES1", "OAT1": "SLC22A6", "OAT3": "SLC22A8", "OCT2": "SLC22A2",
             "MATE1": "SLC47A1", "OATP1B1": "SLCO1B1", "OATP1B3": "SLCO1B3", "HENT1": "SLC29A1", "DCK": "DCK", "CBR1": "CBR1"}
    role_cn = {"elimination": "清除", "efflux": "外排", "hepatic_uptake": "肝摄取", "activation": "活化",
               "toxic_metabolite": "生成毒性代谢物", "toxicity_uptake": "器官摄取（毒性）", "toxicity_efflux": "器官外排",
               "gut_reactivation": "肠道再活化", "cellular_uptake": "细胞摄取", "mixed": "活化/失活并存"}
    for drug, d in kb["drugs"].items():
        for p in d["pathways"]:
            t = p["target"].upper()
            g = alias.get(t, t)
            if clean_gene(g) and not g.startswith("GUT_"):
                out.setdefault(drug, {}).setdefault(g, [g, f"{role_cn.get(p['role'], p['role'])}（{'主要' if p['importance'] == 'major' else '次要'}）", "rule KB"])
    return {k: sorted(v.values()) for k, v in out.items()}


def drugs() -> dict:
    res = json.loads((RAW / "opentargets" / "resolved.json").read_text(encoding="utf-8"))
    kb = json.loads(KB.read_text(encoding="utf-8"))
    adme = adme_table()
    out = []
    by_name = {d["name"]: d for d in res["drugs"]}
    for d in res["drugs"]:
        name = d["name"]
        if name in FINAL27_PD:
            pd_genes, pd_src = FINAL27_PD[name], "final27 文中整理的化疗药药效核心基因"
        elif d["targets"]:
            pd_genes, pd_src = d["targets"], "Open Targets 作用机制靶点"
        elif name in CLASS_PROXY:
            ref = CLASS_PROXY[name]
            pd_genes = FINAL27_PD.get(ref) or by_name[ref]["targets"]
            pd_src = f"同类药物代理（参照 {ref}）"
        elif name in CURATED_PD:
            pd_genes, pd_src = CURATED_PD[name]
        else:
            pd_genes, pd_src = [], "无蛋白靶点（直接作用于 DNA 等），药效学网络分析不可用"
        ae_hits = {}
        for pt, llr in d["faers"]:
            for mod, kws in FAERS_TO_AE.items():
                if any(k in pt.lower() for k in kws) and (mod not in ae_hits or llr > ae_hits[mod][1]):
                    ae_hits[mod] = [pt, llr]
        kb_drug = kb["drugs"].get(name)
        for k in (kb_drug or {}).get("ae", {}):
            for mod in KB_AE_TO_MODULE.get(k, []):
                ae_hits.setdefault(mod, [f"规则库：{kb['ae'].get(k, {}).get('label', k)}", None])
        out.append({"cn": d["cn"], "aliases": d["aliases"], "en": name, "id": d["id"], "type": d["type"],
                    "moa": d["moa"], "pd": sorted(set(pd_genes)), "pd_src": pd_src, "adme": adme.get(name, []),
                    "ae": ae_hits, "faers": d["faers"][:15], "kb": name if kb_drug else None})
    regimens = {k: v for k, v in kb["regimens"].items()}
    print(f"[drugs] {len(out)} drugs; {sum(1 for d in out if d['pd'])} with PD genes, "
          f"{sum(1 for d in out if d['adme'])} with ADME genes")
    return {"source": "Open Targets Platform (mechanisms of action, FAERS), final27 curated PD core, hdi_evidence_kg rule KB",
            "drugs": out, "regimens": regimens}


def cancers() -> dict:
    res = json.loads((RAW / "opentargets" / "resolved.json").read_text(encoding="utf-8"))
    out = [{"cn": c["cn"], "aliases": c["aliases"], "en": c["name"], "id": c["id"], "genes": c["genes"]}
           for c in res["cancers"]]
    print(f"[cancers] {len(out)}")
    return {"source": "Open Targets Platform associatedTargets, top 300 by overall score", "cancers": out}


def kb() -> dict:
    """Rule KB plus the herb Chinese-name index so TCMSP herbs can be matched to KB entries."""
    pharm = json.loads(KB.read_text(encoding="utf-8"))
    lex = json.loads(LEXICON.read_text(encoding="utf-8"))
    herb_cn = {}
    for key, h in lex["herbs"].items():
        for s in h["synonyms"]:
            if re.search(r"[一-鿿]", s):
                herb_cn[s] = key
    return {"pharm": pharm, "herb_cn": herb_cn}


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    target_map = json.loads((RAW / "target_map.json").read_text(encoding="utf-8"))
    write_json("herbs.json", herbs(target_map))
    write_json("drugs.json", drugs())
    write_json("cancers.json", cancers())
    write_json("kb.json", kb())
