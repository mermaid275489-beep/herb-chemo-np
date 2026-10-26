"""Write the input for compare_final27.js: one final27 study's gene set, its PD and cancer
modules, and the final27 network_results rows to compare against.

  python build/validate/make_input.py [Study_ID] [Cancer_Module_ID] [out.json]
  node build/validate/compare_final27.js raw/val_input.json
"""
import json
import sys
from pathlib import Path

import pandas as pd

D = Path(__file__).resolve().parents[3] / "outputs" / "network_pharmacology_final27_2026-09-19" / "data"
study = sys.argv[1] if len(sys.argv) > 1 else "Sun_2024"
cancer_id = sys.argv[2] if len(sys.argv) > 2 else "COLORECTAL"
out = Path(sys.argv[3] if len(sys.argv) > 3 else Path(__file__).resolve().parents[2] / "raw" / "val_input.json")
g = pd.read_csv(D / "intervention_gene_relations.csv")
cm = pd.read_csv(D / "cancer_modules.csv")
pdg = pd.read_csv(D / "chemo_pd.csv")
n = pd.read_csv(D / "network_results.csv")
ref = n[(n.Study_ID == study) & (n.Context_type.ne("CANCER") | n.Context_ID.eq(cancer_id))]
ref = ref.drop_duplicates(["Context_type", "Context_ID"])[["Context_type", "Context_ID", "Set_A_n", "Set_B_n", "Overlap_n",
                                                           "Hypergeometric_P", "Distance", "Random_mean", "Random_SD", "Z"]]
out.write_text(json.dumps({"genes": sorted(set(g[g.Study_ID == study].Gene.dropna().str.upper())),
           "cancer": sorted(set(cm[cm.Cancer_Module_ID == cancer_id].Gene)),
           "pd": sorted(set(pdg[pdg.Study_ID == study].Gene)), "ref": ref.to_dict("records")},
          ensure_ascii=False), encoding="utf-8")
print(f"wrote {out}")
