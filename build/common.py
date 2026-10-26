"""Shared paths and helpers for the site build."""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]           # herb_chemo_np/
PROJECT = ROOT.parent                                  # Clinical_HDS_Cancer_Search_2026/
RAW = ROOT / "raw"
SITE_DATA = ROOT / "docs" / "data"

FINAL27 = PROJECT / "outputs" / "network_pharmacology_final27_2026-09-19" / "data"
STRING_CACHE = PROJECT / "outputs" / "network_pharmacology_17_studies" / "rebuilt" / "data" / "string_v12_cache"
GMT_DIR = PROJECT / "outputs" / "network_pharmacology_new_formulas_2026-08-25" / "raw" / "enrichr_libraries"
LEGACY_TARGET_MAP = PROJECT / "outputs" / "network_pharmacology_17_studies" / "rebuilt" / "data" / "tcmspswiss_targets.csv"
KB = PROJECT / "hdi_evidence_kg" / "data" / "pharmacology.json"
LEXICON = PROJECT / "hdi_evidence_kg" / "data" / "lexicon.json"

STRING_MIN_SCORE = 900   # same as the final27 analysis
OB_MIN, DL_MIN = 30.0, 0.18
TOP_N_CANCER = 300
OT_API = "https://api.platform.opentargets.org/api/v4/graphql"


def clean_gene(value) -> str:
    """Identical to clean_gene() in the final27 pipeline, so the network matches it node for node."""
    s = str(value or "").strip().upper()
    return s if re.fullmatch(r"[A-Z0-9][A-Z0-9._-]*", s) else ""


def string_symbols() -> dict[str, str]:
    out = {}
    with gzip.open(STRING_CACHE / "9606.protein.info.v12.0.txt.gz", "rt", encoding="utf-8") as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                out[parts[0]] = clean_gene(parts[1])
    return out


def write_json(name: str, obj, compact: bool = True) -> Path:
    SITE_DATA.mkdir(parents=True, exist_ok=True)
    path = SITE_DATA / name
    path.write_text(json.dumps(obj, ensure_ascii=False, separators=(",", ":") if compact else None,
                               indent=None if compact else 1), encoding="utf-8")
    print(f"  wrote {path.relative_to(ROOT)} ({path.stat().st_size / 1e6:.2f} MB)")
    return path
