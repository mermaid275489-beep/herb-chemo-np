"""Resolve the curated cancer and drug lists against Open Targets and cache the results.

Cancers: top 300 associated targets by overall score (as in final27); the nine final27
cancer modules keep their exact MONDO IDs. Drugs: mechanism-of-action targets and
significant FAERS adverse events. Anything not listed here can still be queried live
from the page.

  python build/fetch_opentargets.py            # resolve + fetch (cached under raw/opentargets)
"""
from __future__ import annotations

import json
import sys
import time

import requests

from common import OT_API, RAW, TOP_N_CANCER, clean_gene

CACHE = RAW / "opentargets"

# (Chinese name, aliases, English query or fixed ontology ID)
CANCERS = [
    ("肺癌", ["肺恶性肿瘤"], "MONDO_0005138"),
    ("非小细胞肺癌", ["NSCLC"], "MONDO_0005233"),
    ("小细胞肺癌", ["SCLC"], "small cell lung carcinoma"),
    ("肺腺癌", [], "lung adenocarcinoma"),
    ("肺鳞癌", ["肺鳞状细胞癌"], "squamous cell lung carcinoma"),
    ("乳腺癌", [], "MONDO_0004989"),
    ("三阴性乳腺癌", ["TNBC"], "triple-negative breast cancer"),
    ("结直肠癌", ["大肠癌", "CRC"], "MONDO_0005575"),
    ("结肠癌", [], "colon carcinoma"),
    ("直肠癌", [], "rectum carcinoma"),
    ("胃癌", [], "gastric carcinoma"),
    ("食管癌", ["食道癌"], "esophageal carcinoma"),
    ("食管鳞癌", ["食管鳞状细胞癌"], "esophageal squamous cell carcinoma"),
    ("肝癌", ["肝细胞癌", "原发性肝癌", "HCC"], "hepatocellular carcinoma"),
    ("胆管癌", [], "cholangiocarcinoma"),
    ("胆囊癌", [], "gallbladder carcinoma"),
    ("胰腺癌", ["胰腺导管腺癌"], "MONDO_0006047"),
    ("鼻咽癌", [], "MONDO_0015459"),
    ("头颈部鳞癌", ["头颈癌", "头颈部肿瘤"], "head and neck squamous cell carcinoma"),
    ("口腔癌", [], "oral cavity carcinoma"),
    ("喉癌", [], "laryngeal carcinoma"),
    ("甲状腺癌", [], "thyroid carcinoma"),
    ("卵巢癌", [], "MONDO_0005140"),
    ("宫颈癌", ["子宫颈癌"], "MONDO_0005131"),
    ("子宫内膜癌", [], "MONDO_0002447"),
    ("前列腺癌", [], "prostate carcinoma"),
    ("膀胱癌", ["尿路上皮癌"], "urinary bladder carcinoma"),
    ("肾癌", ["肾细胞癌"], "renal cell carcinoma"),
    ("胶质瘤", [], "glioma"),
    ("胶质母细胞瘤", [], "glioblastoma"),
    ("黑色素瘤", [], "melanoma"),
    ("骨肉瘤", [], "osteosarcoma"),
    ("软组织肉瘤", ["肉瘤"], "soft tissue sarcoma"),
    ("急性髓系白血病", ["AML", "急性髓细胞白血病"], "acute myeloid leukemia"),
    ("急性淋巴细胞白血病", ["ALL"], "acute lymphoblastic leukemia"),
    ("慢性髓系白血病", ["CML", "慢性粒细胞白血病"], "chronic myelogenous leukemia"),
    ("慢性淋巴细胞白血病", ["CLL"], "chronic lymphocytic leukemia"),
    ("白血病", [], "leukemia"),
    ("多发性骨髓瘤", ["骨髓瘤"], "plasma cell myeloma"),
    ("弥漫大B细胞淋巴瘤", ["DLBCL"], "diffuse large B-cell lymphoma"),
    ("霍奇金淋巴瘤", [], "Hodgkins lymphoma"),
    ("非霍奇金淋巴瘤", [], "non-Hodgkin lymphoma"),
    ("淋巴瘤", [], "lymphoma"),
    ("骨髓增生异常综合征", ["MDS"], "myelodysplastic syndrome"),
    ("间皮瘤", [], "mesothelioma"),
    ("神经母细胞瘤", [], "neuroblastoma"),
    ("胃肠间质瘤", ["GIST"], "gastrointestinal stromal tumor"),
    ("基底细胞癌", [], "basal cell carcinoma"),
    ("皮肤鳞癌", ["皮肤鳞状细胞癌"], "cutaneous squamous cell carcinoma"),
    ("睾丸生殖细胞肿瘤", ["睾丸癌"], "testicular germ cell tumor"),
]

# (Chinese name, aliases, English name for the Open Targets search)
DRUGS = [
    ("顺铂", ["DDP", "CDDP"], "cisplatin"), ("卡铂", [], "carboplatin"), ("奥沙利铂", ["L-OHP"], "oxaliplatin"),
    ("奈达铂", [], "nedaplatin"), ("洛铂", [], "lobaplatin"),
    ("紫杉醇", ["白蛋白紫杉醇", "紫杉醇脂质体", "nab-paclitaxel"], "paclitaxel"), ("多西他赛", ["多西紫杉醇"], "docetaxel"),
    ("卡巴他赛", [], "cabazitaxel"), ("艾立布林", ["艾日布林"], "eribulin"),
    ("长春新碱", [], "vincristine"), ("长春瑞滨", ["去甲长春花碱"], "vinorelbine"), ("长春碱", [], "vinblastine"),
    ("长春地辛", [], "vindesine"),
    ("氟尿嘧啶", ["5-FU", "5-氟尿嘧啶"], "fluorouracil"), ("卡培他滨", [], "capecitabine"), ("替加氟", ["替吉奥", "S-1"], "tegafur"),
    ("吉西他滨", [], "gemcitabine"), ("培美曲塞", [], "pemetrexed"), ("甲氨蝶呤", ["MTX"], "methotrexate"),
    ("阿糖胞苷", [], "cytarabine"), ("雷替曲塞", [], "raltitrexed"), ("氟达拉滨", [], "fludarabine"),
    ("地西他滨", [], "decitabine"), ("阿扎胞苷", [], "azacitidine"), ("巯嘌呤", ["6-巯基嘌呤"], "mercaptopurine"),
    ("羟基脲", [], "hydroxyurea"), ("曲氟尿苷", ["曲氟尿苷替匹嘧啶", "TAS-102"], "trifluridine"),
    ("亚叶酸钙", ["亚叶酸", "四氢叶酸", "LV"], "leucovorin"),
    ("伊立替康", ["CPT-11"], "irinotecan"), ("拓扑替康", [], "topotecan"), ("依托泊苷", ["VP-16"], "etoposide"),
    ("多柔比星", ["阿霉素", "脂质体多柔比星"], "doxorubicin"), ("表柔比星", ["表阿霉素"], "epirubicin"),
    ("吡柔比星", [], "pirarubicin"), ("柔红霉素", [], "daunorubicin"), ("伊达比星", [], "idarubicin"),
    ("米托蒽醌", [], "mitoxantrone"),
    ("环磷酰胺", [], "cyclophosphamide"), ("异环磷酰胺", [], "ifosfamide"), ("替莫唑胺", [], "temozolomide"),
    ("达卡巴嗪", [], "dacarbazine"), ("苯达莫司汀", [], "bendamustine"), ("美法仑", ["马法兰"], "melphalan"),
    ("白消安", [], "busulfan"), ("卡莫司汀", [], "carmustine"), ("洛莫司汀", [], "lomustine"),
    ("博来霉素", [], "bleomycin"), ("丝裂霉素", [], "mitomycin"), ("曲贝替定", [], "trabectedin"),
    ("三氧化二砷", ["亚砷酸"], "arsenic trioxide"), ("维A酸", ["全反式维甲酸", "ATRA"], "tretinoin"),
    ("吉非替尼", [], "gefitinib"), ("厄洛替尼", [], "erlotinib"), ("埃克替尼", [], "icotinib"), ("阿法替尼", [], "afatinib"),
    ("奥希替尼", [], "osimertinib"), ("克唑替尼", [], "crizotinib"), ("阿来替尼", [], "alectinib"), ("劳拉替尼", [], "lorlatinib"),
    ("索拉非尼", [], "sorafenib"), ("仑伐替尼", [], "lenvatinib"), ("瑞戈非尼", [], "regorafenib"), ("舒尼替尼", [], "sunitinib"),
    ("阿帕替尼", [], "rivoceranib"), ("安罗替尼", [], "anlotinib"), ("培唑帕尼", [], "pazopanib"), ("阿昔替尼", [], "axitinib"),
    ("卡博替尼", [], "cabozantinib"), ("呋喹替尼", [], "fruquintinib"),
    ("伊马替尼", [], "imatinib"), ("达沙替尼", [], "dasatinib"), ("尼洛替尼", [], "nilotinib"),
    ("拉帕替尼", [], "lapatinib"), ("吡咯替尼", [], "pyrotinib"), ("奈拉替尼", [], "neratinib"),
    ("哌柏西利", [], "palbociclib"), ("阿贝西利", [], "abemaciclib"), ("瑞波西利", [], "ribociclib"),
    ("奥拉帕利", [], "olaparib"), ("尼拉帕利", [], "niraparib"), ("伊布替尼", [], "ibrutinib"), ("依维莫司", [], "everolimus"),
    ("硼替佐米", [], "bortezomib"), ("来那度胺", [], "lenalidomide"), ("沙利度胺", [], "thalidomide"),
    ("维莫非尼", [], "vemurafenib"), ("达拉非尼", [], "dabrafenib"), ("曲美替尼", [], "trametinib"),
    ("贝伐珠单抗", [], "bevacizumab"), ("西妥昔单抗", [], "cetuximab"), ("曲妥珠单抗", [], "trastuzumab"),
    ("帕妥珠单抗", [], "pertuzumab"), ("利妥昔单抗", [], "rituximab"), ("雷莫西尤单抗", [], "ramucirumab"),
    ("帕博利珠单抗", ["K药"], "pembrolizumab"), ("纳武利尤单抗", ["O药"], "nivolumab"), ("阿替利珠单抗", [], "atezolizumab"),
    ("度伐利尤单抗", [], "durvalumab"), ("信迪利单抗", [], "sintilimab"), ("卡瑞利珠单抗", [], "camrelizumab"),
    ("替雷利珠单抗", [], "tislelizumab"), ("特瑞普利单抗", [], "toripalimab"), ("伊匹木单抗", [], "ipilimumab"),
    ("他莫昔芬", ["三苯氧胺"], "tamoxifen"), ("来曲唑", [], "letrozole"), ("阿那曲唑", [], "anastrozole"),
    ("依西美坦", [], "exemestane"), ("氟维司群", [], "fulvestrant"), ("比卡鲁胺", [], "bicalutamide"),
    ("恩扎卢胺", [], "enzalutamide"), ("阿比特龙", [], "abiraterone"), ("戈舍瑞林", [], "goserelin"), ("亮丙瑞林", [], "leuprolide"),
]


def gql(query: str, variables: dict) -> dict:
    for attempt in range(4):
        try:
            r = requests.post(OT_API, json={"query": query, "variables": variables}, timeout=90)
            r.raise_for_status()
            out = r.json()
            if out.get("errors"):
                raise RuntimeError(out["errors"])
            return out["data"]
        except (requests.RequestException, RuntimeError) as e:
            print(f"  retry {attempt + 1}: {e}", flush=True)
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"Open Targets query failed: {variables}")


def cached(name: str, fn):
    path = CACHE / f"{name}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    data = fn()
    path.write_text(json.dumps(data, ensure_ascii=False, indent=1), encoding="utf-8")
    time.sleep(0.3)
    return data


SEARCH = """query($q: String!, $e: [String!]!) { search(queryString: $q, entityNames: $e, page: {index: 0, size: 5}) {
  hits { id name entity } } }"""
DISEASE = """query($id: String!) { disease(efoId: $id) { id name
  associatedTargets(page: {index: 0, size: %d}) { count rows { score target { approvedSymbol } } } } }""" % TOP_N_CANCER
DRUG = """query($id: String!) { drug(chemblId: $id) { id name drugType maximumClinicalStage
  mechanismsOfAction { rows { mechanismOfAction actionType targets { approvedSymbol } } }
  adverseEvents(page: {index: 0, size: 30}) { count criticalValue rows { name meddraCode count logLR } } } }"""


def resolve(query: str, entity: str) -> dict | None:
    if "_" in query and query.split("_")[0] in {"MONDO", "EFO", "HP", "Orphanet"}:
        return {"id": query, "name": None}
    hits = cached(f"search_{entity}_{query}", lambda: gql(SEARCH, {"q": query, "e": [entity]}))["search"]["hits"]
    exact = [h for h in hits if (h["name"] or "").lower() == query.lower()]
    return (exact or hits or [None])[0]


def cancers() -> list[dict]:
    out = []
    for cn, aliases, q in CANCERS:
        hit = resolve(q, "disease")
        if not hit:
            print(f"  [cancer] not found: {cn} ({q})")
            continue
        d = cached(f"disease_{hit['id']}", lambda: gql(DISEASE, {"id": hit["id"]}))["disease"]
        genes = [clean_gene(r["target"]["approvedSymbol"]) for r in d["associatedTargets"]["rows"]]
        out.append({"cn": cn, "aliases": aliases, "query": q, "id": d["id"], "name": d["name"],
                    "n_assoc": d["associatedTargets"]["count"], "genes": [g for g in genes if g]})
        print(f"  [cancer] {cn:<10} -> {d['id']:<16} {d['name']}  ({len(genes)} genes)")
    return out


def drugs() -> list[dict]:
    out = []
    for cn, aliases, q in DRUGS:
        hit = resolve(q, "drug")
        if not hit:
            print(f"  [drug] not found: {cn} ({q})")
            continue
        d = cached(f"drug_{hit['id']}", lambda: gql(DRUG, {"id": hit["id"]}))["drug"]
        moa = d.get("mechanismsOfAction") or {"rows": []}
        targets = sorted({clean_gene(t["approvedSymbol"]) for r in moa["rows"] for t in (r["targets"] or [])} - {""})
        ae = d.get("adverseEvents") or {}
        out.append({"cn": cn, "aliases": aliases, "query": q, "id": d["id"], "name": d["name"].lower(),
                    "type": d["drugType"], "phase": d.get("maximumClinicalStage"),
                    "moa": sorted({r["mechanismOfAction"] for r in moa["rows"]}), "targets": targets,
                    "faers": [[r["name"], round(r["logLR"], 1)] for r in (ae.get("rows") or [])]})
        print(f"  [drug] {cn:<8} -> {d['id']:<14} {d['name']:<22} targets={','.join(targets)[:60]}")
    return out


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    c, d = cancers(), drugs()
    (CACHE / "resolved.json").write_text(json.dumps({"cancers": c, "drugs": d}, ensure_ascii=False, indent=1),
                                         encoding="utf-8")
    print(f"[opentargets] {len(c)}/{len(CANCERS)} cancers, {len(d)}/{len(DRUGS)} drugs")
