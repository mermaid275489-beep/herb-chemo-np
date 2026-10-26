"use strict";
/* 中药–化疗网络药理预测 — everything runs in the browser.
   Methods mirror the final27 analysis: STRING v12 (combined score >= 900) as the network and
   universe; hypergeometric overlap; Guney "closest" proximity with degree-matched randomization
   (pool widened until >= 20 nodes); BH across all tests of one query, separately for overlap and
   proximity; significant = (overlap > 0 and q < 0.05) or (Z < -2 and q < 0.05); enrichment of
   module-adjacent herb genes against Reactome and GO BP with BH per library. */

const OT_API = "https://api.platform.opentargets.org/api/v4/graphql";
const SEED_BASE = 20260829;
const S = { data: null, G: null, herbs: [], drugs: [], cancer: null, pathways: null, result: null };
const $ = id => document.getElementById(id);
const esc = s => String(s ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const fmt = (x, d = 2) => x == null || !Number.isFinite(x) ? "—" : x.toFixed(d);
const fmtP = p => p == null || !Number.isFinite(p) ? "—" : p < 1e-3 ? p.toExponential(1) : p.toFixed(3);
/* yield so the progress text can paint; in a hidden tab timers are throttled to ~1 s, so use a
   message-channel task there instead */
const sleep0 = () => new Promise(r => {
  if (document.visibilityState === "visible") return setTimeout(r, 0);
  const ch = new MessageChannel(); ch.port1.onmessage = () => r(); ch.port2.postMessage(0);
});
const norm = s => String(s || "").toLowerCase().replace(/[\s\-_'’.,()]/g, "");
const isAscii = s => /^[\x20-\x7e]+$/.test(s);

const REGIMENS = {
  FOLFOX: ["fluorouracil", "oxaliplatin", "leucovorin"], mFOLFOX6: ["fluorouracil", "oxaliplatin", "leucovorin"],
  FOLFIRI: ["fluorouracil", "irinotecan", "leucovorin"], FOLFIRINOX: ["fluorouracil", "irinotecan", "oxaliplatin", "leucovorin"],
  CAPOX: ["capecitabine", "oxaliplatin"], XELOX: ["capecitabine", "oxaliplatin"], SOX: ["tegafur", "oxaliplatin"],
  FLOT: ["fluorouracil", "leucovorin", "oxaliplatin", "docetaxel"], GEMOX: ["gemcitabine", "oxaliplatin"],
  AC: ["doxorubicin", "cyclophosphamide"], "AC-T": ["doxorubicin", "cyclophosphamide", "paclitaxel"],
  "TC（紫杉醇+卡铂）": ["paclitaxel", "carboplatin"], "TC（多西他赛+环磷酰胺）": ["docetaxel", "cyclophosphamide"],
  TP: ["paclitaxel", "cisplatin"], GP: ["gemcitabine", "cisplatin"], EP: ["etoposide", "cisplatin"], PF: ["cisplatin", "fluorouracil"],
  CHOP: ["cyclophosphamide", "doxorubicin", "vincristine"], "R-CHOP": ["rituximab", "cyclophosphamide", "doxorubicin", "vincristine"],
  ABVD: ["doxorubicin", "bleomycin", "vinblastine", "dacarbazine"], BEP: ["bleomycin", "etoposide", "cisplatin"],
};
const EXAMPLES = [
  { herbs: ["黄芪", "党参", "白术", "茯苓", "甘草"], drugs: ["FOLFOX"], cancer: "结直肠癌", label: "四君子汤加黄芪 × FOLFOX × 结直肠癌" },
  { herbs: ["人参"], drugs: ["伊马替尼"], cancer: "慢性髓系白血病", label: "人参 × 伊马替尼 × CML" },
  { herbs: ["黄连"], drugs: ["他莫昔芬"], cancer: "乳腺癌", label: "黄连 × 他莫昔芬 × 乳腺癌" },
  { herbs: ["丹参"], drugs: ["多西他赛"], cancer: "胃癌", label: "丹参 × 多西他赛 × 胃癌" },
  { herbs: ["黄芩", "甘草", "白芍", "大枣"], drugs: ["伊立替康"], cancer: "结直肠癌", label: "黄芩汤 × 伊立替康 × 结直肠癌" },
];

/* ------------------------------------------------------------------ data */
async function loadData() {
  const names = ["herbs", "drugs", "cancers", "network", "ae", "kb"];
  const parts = await Promise.all(names.map(async n => {
    const r = await fetch(`data/${n}.json`);
    if (!r.ok) throw new Error(`data/${n}.json: HTTP ${r.status}`);
    return r.json();
  }));
  S.data = Object.fromEntries(names.map((n, i) => [n, parts[i]]));
  S.G = buildGraph(S.data.network);
  const h = S.data.herbs;
  S.herbIndex = h.herbs.map((x, i) => ({ i, cn: x.cn, py: x.py.toLowerCase(), la: x.la.toLowerCase(),
    keys: [x.cn, x.py, x.la].map(norm), alias: (x.alias || []).map(norm) }));
  S.drugByEn = new Map(S.data.drugs.drugs.map(d => [d.en, d]));
}

function buildGraph(net) {
  const n = net.genes.length, E = net.edges, deg = new Int32Array(n);
  for (let i = 0; i < E.length; i += 2) { deg[E[i]]++; deg[E[i + 1]]++; }
  const off = new Int32Array(n + 1);
  for (let i = 0; i < n; i++) off[i + 1] = off[i] + deg[i];
  const nb = new Int32Array(off[n]), pos = off.slice(0, n);
  for (let i = 0; i < E.length; i += 2) { const a = E[i], b = E[i + 1]; nb[pos[a]++] = b; nb[pos[b]++] = a; }
  const byDeg = new Map();
  for (let i = 0; i < n; i++) { if (!byDeg.has(deg[i])) byDeg.set(deg[i], []); byDeg.get(deg[i]).push(i); }
  return { n, genes: net.genes, idx: new Map(net.genes.map((g, i) => [g, i])), deg, off, nb, byDeg, pools: new Map() };
}

/* ------------------------------------------------------------------ statistics */
let LF = null;
function lchoose(n, k) {
  if (!LF || LF.length <= n) { LF = new Float64Array(n + 2); for (let i = 1; i < LF.length; i++) LF[i] = LF[i - 1] + Math.log(i); }
  return k < 0 || k > n ? -Infinity : LF[n] - LF[k] - LF[n - k];
}
function hyperSF(k, N, K, n) {          // P(X >= k), X ~ Hypergeometric(N, K, n)
  if (k <= 0) return 1;
  const hi = Math.min(K, n);
  if (k > hi) return 0;
  const base = lchoose(N, n);
  let s = 0;
  for (let i = k; i <= hi; i++) s += Math.exp(lchoose(K, i) + lchoose(N - K, n - i) - base);
  return Math.min(1, s);
}
function bh(ps) {
  const p = ps.map(x => x == null || !Number.isFinite(x) ? 1 : Math.min(Math.max(x, 0), 1));
  const order = p.map((_, i) => i).sort((a, b) => p[a] - p[b]);
  const out = new Array(p.length);
  let run = 1;
  for (let r = p.length - 1; r >= 0; r--) { run = Math.min(run, p[order[r]] * p.length / (r + 1)); out[order[r]] = run; }
  return out;
}
function mulberry32(a) {
  return () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
}
const seedOf = key => SEED_BASE + [...key].reduce((s, c) => s + c.codePointAt(0), 0);

function pool(G, d) {
  if (G.pools.has(d)) return G.pools.get(d);
  let out = [];
  for (let w = 0; w <= 100; w++) {
    out = [];
    for (let k = Math.max(1, d - w); k <= d + w; k++) { const b = G.byDeg.get(k); if (b) for (const x of b) out.push(x); }
    if (out.length >= 20) break;
  }
  const arr = Int32Array.from(out);
  G.pools.set(d, arr);
  return arr;
}
function bfs(G, sources) {
  const dist = new Int32Array(G.n).fill(-1), q = new Int32Array(G.n);
  let h = 0, t = 0;
  for (const s of sources) if (dist[s] < 0) { dist[s] = 0; q[t++] = s; }
  while (h < t) {
    const u = q[h++];
    for (let k = G.off[u]; k < G.off[u + 1]; k++) { const v = G.nb[k]; if (dist[v] < 0) { dist[v] = dist[u] + 1; q[t++] = v; } }
  }
  return dist;
}
function meanDist(dist, idx) {
  let s = 0, c = 0;
  for (const a of idx) if (dist[a] >= 0) { s += dist[a]; c++; }
  return c ? s / c : null;
}
/* Herb-background specificity: share of TCMSP herbs (same filters) whose closest distance to the
   module is <= the query's. All TCMSP target sets lean towards hub genes, so this separates
   "close to this module" from "close to everything". */
function herbPercentile(dist, obs, bg, exclude) {
  let n = 0, le = 0;
  for (const h of bg) {
    if (exclude.has(h.i)) continue;
    const d = meanDist(dist, h.idx);
    if (d === null) continue;
    n++; if (d <= obs) le++;
  }
  return n ? { pct: 100 * le / n, nBg: n } : {};
}
function proximity(G, A, B, key, nRand, bg, exclude) {
  if (!A.length || !B.length) return { status: "无可映射基因" };
  const dist = bfs(G, B);
  const obs = meanDist(dist, A);
  if (obs === null) return { status: "网络不连通" };
  const rng = mulberry32(seedOf(key)), pools = A.map(a => pool(G, G.deg[a])), vals = [];
  for (let r = 0; r < nRand; r++) {
    let ss = 0, cc = 0;
    for (const p of pools) { const d = dist[p[Math.floor(rng() * p.length)]]; if (d >= 0) { ss += d; cc++; } }
    if (cc) vals.push(ss / cc);
  }
  const mu = vals.reduce((x, y) => x + y, 0) / vals.length;
  const sd = Math.sqrt(vals.reduce((x, y) => x + (y - mu) ** 2, 0) / vals.length) || 1e-9;
  const inB = new Uint8Array(G.n); for (const b of B) inB[b] = 1;
  let oneHop = 0;
  for (const a of A) for (let k = G.off[a]; k < G.off[a + 1]; k++) if (inB[G.nb[k]]) oneHop++;
  return { d: obs, mu, sd, z: (obs - mu) / sd, p: (vals.filter(v => v <= obs).length + 1) / (vals.length + 1), oneHop, status: "ok",
    ...(bg ? herbPercentile(dist, obs, bg, exclude) : {}) };
}

/* ------------------------------------------------------------------ query sets */
function herbTargets(herbIdxs, opt) {
  const H = S.data.herbs, geneComp = new Map(), per = [];
  for (const hi of herbIdxs) {
    const herb = H.herbs[hi];
    let pass = 0;
    const genes = new Set();
    for (const ci of herb.c) {
      const c = H.compounds[ci];
      if (opt.obdl && !(c[2] >= H.ob_min && c[3] >= H.dl_min)) continue;
      const tg = opt.validated ? c[5] : c[4];
      if (!tg.length) continue;
      pass++;
      for (const g of tg) {
        const sym = H.genes[g];
        genes.add(sym);
        if (!geneComp.has(sym)) geneComp.set(sym, new Map());
        geneComp.get(sym).set(ci, herb.cn);
      }
    }
    per.push({ cn: herb.cn, nIngredients: herb.c.length, nComp: pass, nGenes: genes.size });
  }
  return { genes: [...geneComp.keys()].sort(), geneComp, per };
}
function drugSets(drugs) {
  const pd = new Map(), adme = [];
  for (const d of drugs) {
    for (const g of d.pd || []) { if (!pd.has(g)) pd.set(g, []); pd.get(g).push(d.cn || d.en); }
    for (const [g, role, src] of d.adme || []) adme.push({ gene: g, role, src, drug: d.cn || d.en });
  }
  return { pd, adme };
}
const toIdx = genes => [...new Set(genes)].map(g => S.G.idx.get(g)).filter(i => i !== undefined).sort((a, b) => a - b);
/* every single TCMSP herb with >= 5 network genes under the current filters (cached per filter setting) */
function herbBackground(opt) {
  const key = `${opt.obdl}|${opt.validated}`;
  if (S.bg?.key === key) return S.bg.sets;
  const sets = S.data.herbs.herbs.map((_, i) => ({ i, idx: Int32Array.from(toIdx(herbTargets([i], opt).genes)) })).filter(h => h.idx.length >= 5);
  S.bg = { key, sets };
  return sets;
}

/* ------------------------------------------------------------------ analysis */
async function analyse(opt, progress) {
  const G = S.G, N = G.n;
  const herb = herbTargets(S.herbs, opt);
  const { pd, adme } = drugSets(S.drugs);
  const A = toIdx(herb.genes), Aset = new Set(herb.genes);
  const contexts = [];
  if (S.cancer) contexts.push({ type: "CANCER", id: S.cancer.id, label: `癌种：${S.cancer.cn || S.cancer.en}`, genes: S.cancer.genes });
  if (pd.size) contexts.push({ type: "CHEMO_PD", id: "PD", label: "化疗药作用靶点", genes: [...pd.keys()] });
  for (const m of S.data.ae.modules) {
    const core = m.core.map(x => x[0]);
    contexts.push({ type: "AE_CORE", id: m.id, label: `${m.cn}（核心）`, ae: m, genes: core });
    contexts.push({ type: "AE_EXTENDED", id: m.id, label: `${m.cn}（扩展）`, ae: m, genes: [...new Set([...core, ...m.extended])] });
  }
  const herbKey = S.herbs.map(i => S.data.herbs.herbs[i].cn).join("+");
  const nRand = opt.fast ? 200 : 1000;
  progress("准备中药背景分布…"); await sleep0();
  const bg = herbBackground(opt), exclude = new Set(S.herbs);
  for (let i = 0; i < contexts.length; i++) {
    const c = contexts[i];
    const B = toIdx(c.genes), Bset = new Set(c.genes.filter(g => G.idx.has(g)));
    const overlap = herb.genes.filter(g => G.idx.has(g) && Bset.has(g));
    c.nB = B.length;
    c.overlap = overlap;
    c.pOverlap = overlap.length ? hyperSF(overlap.length, N, B.length, A.length) : 1;
    c.prox = proximity(G, A, B, `${herbKey}|${c.type}|${c.id}`, nRand, bg, exclude);
    if (i % 4 === 3) { progress(`网络邻近度计算 ${i + 1}/${contexts.length}`); await sleep0(); }
  }
  const qo = bh(contexts.map(c => c.pOverlap)), qp = bh(contexts.map(c => c.prox.p ?? 1));
  contexts.forEach((c, i) => {
    c.qOverlap = qo[i]; c.qProx = qp[i];
    c.sigOverlap = c.overlap.length > 0 && c.qOverlap < 0.05;
    c.sigProx = c.prox.status === "ok" && c.prox.z < -2 && c.qProx < 0.05;
    c.sig = c.sigOverlap || c.sigProx;
  });
  const pk = adme.filter(x => Aset.has(x.gene));
  const known = {};
  for (const d of S.drugs) for (const [mod, v] of Object.entries(d.ae || {})) (known[mod] ||= []).push([d.cn || d.en, ...v]);
  return { opt, herb, pd, adme, pk, A, contexts, known, herbKey, nRand, kb: kbPredict() };
}

/* near genes: herb genes inside the module or one STRING step from it (as in final27) */
function nearGenes(A, moduleGenes) {
  const G = S.G, inB = new Uint8Array(G.n);
  for (const g of moduleGenes) { const i = G.idx.get(g); if (i !== undefined) inB[i] = 1; }
  return A.filter(a => { if (inB[a]) return true; for (let k = G.off[a]; k < G.off[a + 1]; k++) if (inB[G.nb[k]]) return true; return false; });
}
function enrich(genesIdx) {
  if (genesIdx.length < 3 || !S.pathways) return null;
  const N = S.G.n, inSet = new Uint8Array(N);
  for (const g of genesIdx) inSet[g] = 1;
  const out = {};
  for (const [lib, terms] of Object.entries(S.pathways.libraries)) {
    const rows = [];
    for (const [name, members] of terms) {
      let k = 0;
      for (const m of members) if (inSet[m]) k++;
      if (k < 3) continue;
      rows.push({ name, k, K: members.length, p: hyperSF(k, N, members.length, genesIdx.length), members });
    }
    const q = bh(rows.map(r => r.p));
    rows.forEach((r, i) => r.q = q[i]);
    rows.sort((a, b) => a.q - b.q || b.k - a.k || a.p - b.p);
    out[lib] = { n: rows.length, top: rows.slice(0, 12).map(r => ({ ...r, genes: r.members.filter(m => inSet[m]).map(m => S.G.genes[m]) })) };
  }
  return out;
}

/* ------------------------------------------------------------------ rule-based KB (hdi_evidence_kg) */
const KB_CANCER = [[/肺/, "lung_cancer"], [/乳腺/, "breast_cancer"], [/结直肠|结肠|直肠|大肠/, "colorectal_cancer"], [/胃癌/, "gastric_cancer"],
  [/肝/, "liver_cancer"], [/胰腺/, "pancreatic_cancer"], [/卵巢/, "ovarian_cancer"], [/白血病/, "leukemia"], [/淋巴瘤/, "lymphoma"],
  [/鼻咽|头颈|口腔|喉/, "head_neck_cancer"], [/食管/, "esophageal_cancer"]];
const KB_EV = { clinical: 1, in_vivo: .5, in_vitro: .3 }, KB_EVL = { clinical: "人体", in_vivo: "动物", in_vitro: "体外" };
const KB_IMP = { major: 1, minor: .5 }, KB_SEV = { DLT: 1, serious: 1, common: .8, rare: .5 }, KB_OBS = { L1: 1, L2a: 1, L2b: 1, L3: .8, L4: .5 };
const KB_ROLE = { elimination: "清除", efflux: "外排", hepatic_uptake: "肝摄取", activation: "活化", toxic_metabolite: "生成毒性代谢物",
  toxicity_uptake: "器官摄取（毒性）", toxicity_efflux: "器官外排（解毒）", gut_reactivation: "肠道再活化", cellular_uptake: "肿瘤细胞摄取", mixed: "活化/失活并存" };

function kbPredict() {
  const P = S.data.kb.pharm, cnMap = S.data.kb.herb_cn;
  const herbs = [...new Set(S.herbs.map(i => cnMap[S.data.herbs.herbs[i].cn]).filter(k => k && P.herbs[k]))];
  const drugs = Object.fromEntries(S.drugs.map(d => d.kb || d.en).filter(k => P.drugs[k]).map(k => [k, P.drugs[k]]));
  if (!herbs.length || !Object.keys(drugs).length) return null;
  const cancerName = S.cancer ? (S.cancer.cn || "") : "";
  const ck = (KB_CANCER.find(([re]) => re.test(cancerName)) || [])[1];
  const cancer = (ck && P.cancers[ck]) || {};
  const label = a => (P.ae[a] || {}).label || a;
  const profiles = h => { const p = P.herbs[h]; const out = [[h, p]]; for (const c of p.components || []) if (P.herbs[c]) out.push([`${c}（${h} 组分）`, P.herbs[c]]); return out; };
  const effect = (role, eff) => { const s = eff === "inhibits" ? 1 : -1;
    if (["elimination", "efflux", "hepatic_uptake"].includes(role)) return ["exposure", s]; if (role === "activation") return ["exposure", -s];
    if (role === "cellular_uptake") return ["efficacy", -s]; if (["toxic_metabolite", "toxicity_uptake", "gut_reactivation"].includes(role)) return ["organ", -s];
    if (role === "toxicity_efflux") return ["organ", s]; return ["mixed", 0]; };
  const pk = [];
  for (const h of herbs) for (const [dk, d] of Object.entries(drugs)) {
    const paths = Object.fromEntries(d.pathways.map(p => [p.target, p]));
    for (const [src, prof] of profiles(h)) for (const m of prof.modulates) {
      const p = paths[m.target]; if (!p) continue;
      const [kind, sign] = effect(p.role, m.effect);
      let score = KB_EV[m.evidence] * KB_IMP[p.importance] * (src !== h ? .7 : 1);
      let cons;
      if (kind === "exposure") { score *= cancer.pk_amplify || 1; if (p.role === "efflux" && d.route === "iv") score *= .6;
        cons = sign > 0 ? `暴露↑ → ${d.exposure_ae.map(label).join("、")}风险↑` : "暴露↓ → 疗效可能下降"; }
      else if (kind === "organ") cons = `${label(p.ae)}${sign > 0 ? "↑" : "↓"}`;
      else if (kind === "efficacy") cons = sign < 0 ? "肿瘤细胞内药物↓ → 疗效可能下降" : "肿瘤细胞内药物↑";
      else cons = "方向不确定";
      pk.push({ herb: h, src, m, p, kind, sign, drug: dk, cons, score: Math.round(score * 1000) / 1000 });
    }
  }
  const flags = new Set(cancer.flags || []), pdHits = new Map();
  for (const h of herbs) for (const [src, prof] of profiles(h)) for (const e of prof.pd) {
    const ev = KB_EV[e.evidence] * (src !== h ? .7 : 1);
    for (const r of P.pd_rules) {
      if (r.herb_effect !== e.effect) continue;
      const ch = flags.has(r.cancer_flag); let matched = false;
      for (const [dk, d] of Object.entries(drugs)) {
        const aes = (r.drug_ae || []).filter(a => a in d.ae);
        if (!aes.length && !d.tags.includes(r.drug_tag)) continue;
        matched = true;
        const sev = aes.length ? Math.max(...aes.map(a => KB_SEV[d.ae[a]])) : .8;
        const hit = { herb: h, src, e, drug: dk, r, score: Math.min(ev * r.weight * sev * (ch ? 1.2 : 1), 1.2) };
        const k = [h, e.effect, dk, r.result, r.direction].join("|");
        if (!pdHits.has(k) || hit.score > pdHits.get(k).score) pdHits.set(k, hit);
      }
      if (ch && !matched) { const k = [h, e.effect, "", r.result, r.direction].join("|");
        const hit = { herb: h, src, e, drug: null, r, score: ev * r.weight * .8 };
        if (!pdHits.has(k) || hit.score > pdHits.get(k).score) pdHits.set(k, hit); }
    }
  }
  const risks = new Map();
  const add = (ae, dir, score, text) => {
    if (!risks.has(ae)) risks.set(ae, { ae, label: label(ae), up: 0, down: 0, why: [] });
    const r = risks.get(ae);
    if (dir === "increase") { score *= (cancer.ae_boost || {})[ae] || 1; r.up = Math.max(r.up, score) + (r.up ? .1 : 0); } else r.down = Math.max(r.down, score);
    r.why.push([dir, text]);
  };
  const pairs = new Map();
  for (const s of pk) { const k = s.herb + "|" + s.drug; if (!pairs.has(k)) pairs.set(k, []); pairs.get(k).push(s); }
  for (const steps of pairs.values()) {
    const { herb, drug } = steps[0], d = drugs[drug], exp = steps.filter(s => s.kind === "exposure");
    if (exp.length) {
      const top = Math.max(...exp.map(s => KB_EV[s.m.evidence])), lead = exp.filter(s => KB_EV[s.m.evidence] === top);
      const v = lead.reduce((a, s) => a + s.sign * s.score, 0) + .3 * exp.filter(s => !lead.includes(s)).reduce((a, s) => a + s.sign * s.score, 0);
      const dir = v > .05 ? "↑" : v < -.05 ? "↓" : "≈";
      const text = `${herb} 经 ${[...new Set(exp.map(s => s.m.target))].join("/")} 使 ${drug} 暴露${dir}（以${KB_EVL[lead[0].m.evidence]}证据为主）`;
      if (v > .05) for (const ae of d.exposure_ae) add(ae, "increase", Math.min(v, 1) * (KB_SEV[d.ae[ae] || "common"] || .8), text);
      else if (v < -.05) add("efficacy_loss", "increase", Math.min(-v, 1), text);
    }
    for (const s of steps) if (s.kind === "organ") add(s.p.ae, s.sign > 0 ? "increase" : "decrease", s.score, `${s.src}：${s.cons}`);
  }
  for (const h of pdHits.values()) add(h.r.result, h.r.direction, h.score, `${h.src}（${h.e.effect}）× ${h.drug || "癌种背景"}：${h.r.explain}`);
  const observed = P.observed_pairs.filter(o => herbs.includes(o.herb) && o.drug in drugs);
  for (const o of observed) if (o.ae) add(o.ae, o.direction, KB_OBS[o.level] ?? .5, `${o.finding}（${o.level}，${o.ref}）`);
  const ae = [...risks.values()].map(r => ({ ...r, up: Math.min(r.up, 1.5) })).sort((a, b) => b.up - a.up || b.down - a.down);
  return { herbs, drugs: Object.keys(drugs), ae, pk, observed };
}

/* ------------------------------------------------------------------ live Open Targets */
async function ot(query, variables) {
  const r = await fetch(OT_API, { method: "POST", headers: { "content-type": "application/json" }, body: JSON.stringify({ query, variables }) });
  const j = await r.json();
  if (j.errors) throw new Error(j.errors[0].message);
  return j.data;
}
const OT_SEARCH = `query($q:String!,$e:[String!]!){search(queryString:$q,entityNames:$e,page:{index:0,size:8}){hits{id name entity description}}}`;
const OT_DISEASE = `query($id:String!){disease(efoId:$id){id name associatedTargets(page:{index:0,size:300}){rows{target{approvedSymbol}}}}}`;
const OT_DRUG = `query($id:String!){drug(chemblId:$id){id name drugType mechanismsOfAction{rows{mechanismOfAction targets{approvedSymbol}}}
  adverseEvents(page:{index:0,size:30}){rows{name logLR}}}}`;
async function otDisease(id) {
  const d = (await ot(OT_DISEASE, { id })).disease;
  return { id: d.id, en: d.name, cn: "", live: true, genes: d.associatedTargets.rows.map(r => (r.target.approvedSymbol || "").toUpperCase()).filter(Boolean) };
}
async function otDrug(id) {
  const d = (await ot(OT_DRUG, { id })).drug;
  const genes = [...new Set(d.mechanismsOfAction.rows.flatMap(r => r.targets.map(t => (t.approvedSymbol || "").toUpperCase())))].filter(Boolean);
  return { cn: "", en: d.name.toLowerCase(), id: d.id, type: d.drugType, moa: d.mechanismsOfAction.rows.map(r => r.mechanismOfAction),
    pd: genes, pd_src: genes.length ? "Open Targets 作用机制靶点（实时检索）" : "Open Targets 未列出蛋白靶点", adme: [], ae: {}, live: true,
    faers: (d.adverseEvents?.rows || []).map(r => [r.name, Math.round(r.logLR * 10) / 10]) };
}

/* ------------------------------------------------------------------ pickers */
function searchHerbs(q) {
  const n = norm(q);
  if (!n) return [];
  const scored = [];
  for (const h of S.herbIndex) {
    let s = 0;
    if (h.keys[0] === n || h.keys[1] === n) s = 100;
    else if (h.cn.startsWith(q.trim())) s = 80;
    else if (h.keys[1].startsWith(n)) s = 70;
    else if (h.cn.includes(q.trim())) s = 60;
    else if (h.keys[2].includes(n)) s = 50;
    else if (n.length > 2 && h.alias.some(a => a.includes(n))) s = 45;
    if (s) scored.push([s - h.cn.length * .1, h.i]);
  }
  return scored.sort((a, b) => b[0] - a[0]).slice(0, 12).map(x => x[1]);
}
function searchDrugs(q) {
  const n = norm(q), out = [];
  if (!n) return out;
  for (const [name, list] of Object.entries(REGIMENS)) if (norm(name).startsWith(n)) out.push({ kind: "reg", name, list });
  for (const d of S.data.drugs.drugs) {
    const keys = [d.cn, d.en, ...(d.aliases || [])].map(norm);
    if (keys.some(k => k === n)) out.unshift({ kind: "drug", d });
    else if (keys.some(k => k.startsWith(n) || (n.length > 1 && k.includes(n)))) out.push({ kind: "drug", d });
  }
  return out.slice(0, 12);
}
function searchCancers(q) {
  const n = norm(q);
  if (!n) return [];
  return S.data.cancers.cancers.filter(c => [c.cn, c.en, ...(c.aliases || [])].map(norm).some(k => k.includes(n))).slice(0, 10);
}

function setupPicker(kind) {
  const input = $(`${kind}-input`), menu = $(`${kind}-menu`);
  let items = [], active = -1, timer = null, liveSeq = 0;
  const close = () => { menu.hidden = true; active = -1; };
  const render = (extraNote = "") => {
    menu.innerHTML = items.map((it, i) => `<div class="opt${i === active ? " active" : ""}" role="option" data-i="${i}">${it.html}</div>`).join("")
      + (extraNote ? `<div class="note">${extraNote}</div>` : "");
    menu.hidden = !items.length && !extraNote;
  };
  const update = () => {
    const q = input.value.trim();
    items = [];
    if (!q) { close(); return; }
    if (kind === "herb") {
      items = searchHerbs(q).map(i => { const h = S.data.herbs.herbs[i];
        return { html: `<b>${esc(h.cn)}</b><span class="sub">${esc(h.py)} · ${esc(h.la)}</span><span class="tag">${esc(h.cat)}</span>`, pick: () => addHerb(i) }; });
      render(items.length ? "" : "TCMSP 中未找到该药材（可试拼音或拉丁名）");
    } else if (kind === "drug") {
      items = searchDrugs(q).map(x => x.kind === "reg"
        ? { html: `<b>${esc(x.name)}</b><span class="sub">${esc(x.list.map(e => (S.drugByEn.get(e) || {}).cn || e).join(" + "))}</span><span class="tag">方案</span>`, pick: () => x.list.forEach(e => S.drugByEn.get(e) && addDrug(S.drugByEn.get(e))) }
        : { html: `<b>${esc(x.d.cn)}</b><span class="sub">${esc(x.d.en)}${x.d.aliases.length ? " · " + esc(x.d.aliases.join("、")) : ""}</span><span class="tag">${esc(x.d.type || "")}</span>`, pick: () => addDrug(x.d) });
      if (isAscii(q) && q.length >= 3) items.push({ html: `<span class="sub">在 Open Targets 中实时搜索药物 “${esc(q)}”</span><span class="tag">联网</span>`, pick: () => liveSearch("drug", q) });
      render(items.length ? "" : "未找到；其他药物请输入英文通用名以联网检索");
    } else {
      items = searchCancers(q).map(c => ({ html: `<b>${esc(c.cn)}</b><span class="sub">${esc(c.en)}</span><span class="tag">${esc(c.id)}</span>`, pick: () => setCancer(c) }));
      if (isAscii(q) && q.length >= 3) items.push({ html: `<span class="sub">在 Open Targets 中实时搜索疾病 “${esc(q)}”</span><span class="tag">联网</span>`, pick: () => liveSearch("disease", q) });
      render(items.length ? "" : "内置列表中没有；请输入英文病名以联网检索任意疾病");
    }
  };
  const liveSearch = async (entity, q) => {
    const seq = ++liveSeq;
    items = []; render("正在检索 Open Targets…");
    try {
      const hits = (await ot(OT_SEARCH, { q, e: [entity] })).search.hits;
      if (seq !== liveSeq) return;
      items = hits.map(h => ({ html: `<b>${esc(h.name)}</b><span class="sub">${esc(h.description || "")}</span><span class="tag">${esc(h.id)}</span>`,
        pick: async () => { status(`正在从 Open Targets 获取 ${h.name}…`);
          try { entity === "disease" ? setCancer(await otDisease(h.id)) : addDrug(await otDrug(h.id)); status(""); }
          catch (e) { status(`Open Targets 获取失败：${e.message}`, true); } } }));
      active = 0; render(items.length ? "" : "Open Targets 没有找到匹配项");
      input.focus();
    } catch (e) { items = []; render(`Open Targets 检索失败：${esc(e.message)}`); }
  };
  input.addEventListener("input", () => { clearTimeout(timer); timer = setTimeout(update, 120); });
  input.addEventListener("focus", update);
  input.addEventListener("keydown", e => {
    if (e.key === "ArrowDown" && items.length) { active = (active + 1) % items.length; render(); e.preventDefault(); }
    else if (e.key === "ArrowUp" && items.length) { active = (active - 1 + items.length) % items.length; render(); e.preventDefault(); }
    else if (e.key === "Enter") { const it = items[active >= 0 ? active : 0]; if (it) { e.preventDefault(); it.pick(); input.value = ""; close(); } }
    else if (e.key === "Escape") close();
    else if (e.key === "Backspace" && !input.value) { if (kind === "herb") S.herbs.pop(); else if (kind === "drug") S.drugs.pop(); else S.cancer = null; renderChips(); }
  });
  menu.addEventListener("mousedown", e => {
    const el = e.target.closest(".opt"); if (!el) return;
    e.preventDefault();
    const it = items[+el.dataset.i];
    if (it) { const isLive = /实时搜索/.test(el.textContent); it.pick(); if (!isLive) { input.value = ""; close(); } }
  });
  input.addEventListener("blur", () => setTimeout(close, 150));
}
function addHerb(i) { if (!S.herbs.includes(i)) S.herbs.push(i); renderChips(); }
function addDrug(d) { if (!S.drugs.some(x => x.en === d.en)) S.drugs.push(d); renderChips(); }
function setCancer(c) { S.cancer = c; renderChips(); }
function renderChips() {
  const chip = (label, sub, onx) => { const el = document.createElement("span"); el.className = "chip";
    el.innerHTML = `${esc(label)}${sub ? ` <small>${esc(sub)}</small>` : ""}<button type="button" aria-label="移除 ${esc(label)}">×</button>`;
    el.querySelector("button").onclick = onx; return el; };
  const hc = $("herb-chips"); hc.replaceChildren(...S.herbs.map(i => chip(S.data.herbs.herbs[i].cn, "", () => { S.herbs = S.herbs.filter(x => x !== i); renderChips(); })));
  const dc = $("drug-chips"); dc.replaceChildren(...S.drugs.map(d => chip(d.cn || d.en, d.cn ? d.en : d.id, () => { S.drugs = S.drugs.filter(x => x !== d); renderChips(); })));
  const cc = $("cancer-chips"); cc.replaceChildren(...(S.cancer ? [chip(S.cancer.cn || S.cancer.en, S.cancer.cn ? S.cancer.en : S.cancer.id, () => { S.cancer = null; renderChips(); })] : []));
}
function status(msg, err = false) { const el = $("status"); el.textContent = msg; el.classList.toggle("err", err); }

/* ------------------------------------------------------------------ rendering */
/* significance from the final27 rule, graded by herb-background percentile */
function verdict(c) {
  if (c.prox.status !== "ok" && !c.overlap.length) return `<span class="pill na">${esc(c.prox.status === "ok" ? "无重叠" : c.prox.status)}</span>`;
  const pct = c.prox.pct;
  if (c.sig && pct != null && pct <= 10) return `<span class="pill sig">显著 · 高特异</span>`;
  if (c.sig && pct != null && pct <= 25) return `<span class="pill near">显著 · 较特异</span>`;
  if (c.sig) return `<span class="pill none">显著 · 非特异</span>`;
  if (c.overlap.length || (c.prox.z ?? 0) < -1) return `<span class="pill none">弱关联</span>`;
  return `<span class="pill none">无明显关联</span>`;
}
const pctTxt = c => c.prox.pct == null ? "—" : `${c.prox.pct < 1 ? "<1" : Math.round(c.prox.pct)}%`;
const geneList = (genes, max = 40) => genes.length ? genes.slice(0, max).map(g => `<span class="gene">${esc(g)}</span>`).join(", ") + (genes.length > max ? ` 等 ${genes.length} 个` : "") : "—";

function render(res) {
  const R = $("results"), H = S.data.herbs, G = S.G;
  const herbNames = S.herbs.map(i => H.herbs[i].cn).join(" + ");
  const drugNames = S.drugs.map(d => d.cn || d.en).join(" + ");
  const inNet = res.A.length;
  const cancerCtx = res.contexts.find(c => c.type === "CANCER"), pdCtx = res.contexts.find(c => c.type === "CHEMO_PD");
  const out = [];
  out.push(`<section class="panel"><div class="sectionhead"><h2>${esc(herbNames)} × ${esc(drugNames || "—")}${S.cancer ? " × " + esc(S.cancer.cn || S.cancer.en) : ""}</h2>
    <p class="hint">随机化 ${res.nRand} 次 · ${res.opt.obdl ? "OB/DL 筛选" : "未做 OB/DL 筛选"} · ${res.opt.validated ? "仅验证靶点" : "验证 + 预测靶点"}</p></div>
    <div class="strip">
      <div class="stat"><div class="k">中药活性成分</div><div class="v">${res.herb.per.reduce((s, x) => s + x.nComp, 0)}</div><div class="d">${res.herb.per.map(x => `${esc(x.cn)} ${x.nComp}/${x.nIngredients}`).join("；")}</div></div>
      <div class="stat"><div class="k">中药靶点（在 STRING 网络中）</div><div class="v">${res.herb.genes.length}<small class="hint"> / ${inNet}</small></div><div class="d">人类基因，经 UniProt 映射</div></div>
      <div class="stat"><div class="k">化疗药作用靶点</div><div class="v">${res.pd.size}</div><div class="d">${S.drugs.length ? S.drugs.map(d => `<span title="${esc(d.pd_src)}">${esc(d.cn || d.en)} ${d.pd.length}</span>`).join("；") : "未选择化疗药"}</div></div>
      <div class="stat"><div class="k">化疗药代谢/转运基因</div><div class="v">${new Set(res.adme.map(x => x.gene)).size}</div><div class="d">${res.adme.length ? "来自 final27 表与规则库" : "该药暂无收录"}</div></div>
      <div class="stat"><div class="k">癌种基因模块</div><div class="v">${S.cancer ? S.cancer.genes.length : "—"}</div><div class="d">${S.cancer ? `Open Targets top 300 · ${esc(S.cancer.id)}` : "未选择癌种"}</div></div>
    </div>
    ${inNet < 3 ? `<p class="note" style="margin-top:10px">可用于网络分析的中药靶点太少（${inNet} 个）。可取消 OB/DL 筛选或“只用验证靶点”后重试。</p>` : ""}
    ${res.contexts.filter(c => c.sig).length > 0.6 * res.contexts.length ? `<p class="note" style="margin-top:10px">${res.contexts.filter(c => c.sig).length}/${res.contexts.length} 项检验达到显著。TCMSP 靶点普遍偏向网络中的枢纽基因，几乎任何中药都会与多数模块“显著接近”（final27 敏感性分析也观察到这一点）。请主要看<b>中药间百分位</b>：它把本输入与全部 TCMSP 中药比较，百分位越低，说明越特异地靠近该模块。</p>` : ""}
  </section>`);

  const card = (title, c, extra = "") => !c ? "" : `<div class="card"><div class="sectionhead"><h3>${esc(title)}</h3>${verdict(c)}</div>
    <div class="nums"><span>重叠 ${c.overlap.length}/${c.nB}</span><span>重叠 q ${fmtP(c.qOverlap)}</span><span>Z ${fmt(c.prox.z)}</span><span>邻近度 q ${fmtP(c.qProx)}</span><span>d ${fmt(c.prox.d)} / 随机 ${fmt(c.prox.mu)}</span><span title="在 ${c.prox.nBg ?? "—"} 味中药中，与该模块距离不大于本输入的比例">中药间百分位 ${pctTxt(c)}</span></div>
    <div class="genes">重叠基因：${geneList(c.overlap)}</div>${extra}</div>`;
  const pkCard = `<div class="card"><div class="sectionhead"><h3>中药靶点 × 化疗药代谢/转运基因（PK）</h3>${res.pk.length ? `<span class="pill sig">${res.pk.length} 处重叠</span>` : `<span class="pill none">无直接重叠</span>`}</div>
    ${res.pk.length ? `<div class="tbl"><table><tr><th>基因</th><th>化疗药</th><th>该基因在药物处置中的角色</th><th>贡献成分</th></tr>${res.pk.map(x => `<tr><td class="gene">${esc(x.gene)}</td><td>${esc(x.drug)}</td><td>${esc(x.role)}</td><td>${esc(topCompounds(res, x.gene, 3))}</td></tr>`).join("")}</table></div>
      <p class="hint">TCMSP 只给出“成分–靶点”关联，不区分抑制或诱导；方向请看下方规则库推理或查原始文献。</p>`
      : `<p class="hint">${res.adme.length ? "中药靶点中没有该化疗药的代谢酶或转运体。" : "所选化疗药暂未收录代谢/转运基因。"}</p>`}</div>`;
  out.push(`<section class="panel"><h2>关键网络关系</h2><div class="cards">
    ${card(cancerCtx ? `中药靶点 × ${cancerCtx.label.replace("癌种：", "")}基因模块` : "", cancerCtx, `<p class="hint">与癌种模块显著邻近提示中药可能作用于该肿瘤相关网络（潜在协同或拮抗，需实验验证）。</p>`)}
    ${card("中药靶点 × 化疗药作用靶点（PD）", pdCtx, `<p class="hint">邻近提示药效学层面可能存在相互影响。</p>`)}
    ${pkCard}
    ${!cancerCtx ? `<div class="card"><h3>癌种</h3><p class="hint">未选择癌种，跳过癌种模块分析。</p></div>` : ""}
  </div></section>`);

  const aeRows = [];
  const mods = S.data.ae.modules.map(m => ({ m, core: res.contexts.find(c => c.type === "AE_CORE" && c.id === m.id), ext: res.contexts.find(c => c.type === "AE_EXTENDED" && c.id === m.id), known: res.known[m.id] }));
  const best = x => [x.core, x.ext].filter(c => c.sig && c.prox.pct != null).reduce((m, c) => Math.min(m, c.prox.pct), 101);
  mods.sort((a, b) => best(a) - best(b) || Math.min(a.core.prox.z ?? 9, a.ext.prox.z ?? 9) - Math.min(b.core.prox.z ?? 9, b.ext.prox.z ?? 9));
  const pickC = x => [x.core, x.ext].filter(c => c.sig).sort((a, b) => (a.prox.pct ?? 101) - (b.prox.pct ?? 101))[0] || (x.core.overlap.length ? x.core : x.ext);
  for (const x of mods) {
    const dirOf = new Map(x.m.core);
    const coreGenes = x.core.overlap.map(g => `<span class="gene ${dirOf.get(g) === "保护性" ? "g-pro" : "g-inj"}">${esc(g)}</span>`).join(", ") || "—";
    aeRows.push(`<tr class="ae-row" data-id="${x.m.id}" tabindex="0"><td><b>${esc(x.m.cn)}</b><div class="hint">${esc(x.m.group)}</div></td>
      <td>${x.known ? `<span class="pill known" title="${esc(x.known.map(k => `${k[0]}：${k[1]}${k[2] != null ? `（logLR ${k[2]}）` : ""}`).join("；"))}">已知</span>` : ""}</td>
      <td class="num">${x.core.overlap.length}/${x.core.nB}</td><td class="num">${fmt(x.core.prox.z)} <span class="hint">q ${fmtP(x.core.qProx)}</span></td><td class="num">${pctTxt(x.core)}</td>
      <td class="num">${x.ext.overlap.length}/${x.ext.nB}</td><td class="num">${fmt(x.ext.prox.z)} <span class="hint">q ${fmtP(x.ext.qProx)}</span></td><td class="num">${pctTxt(x.ext)}</td>
      <td>${verdict(pickC(x))}</td></tr>
      <tr class="detail" data-for="${x.m.id}" hidden><td colspan="9">
        <div>核心重叠基因（<span class="g-inj">红=促损伤</span>，<span class="g-pro">绿=保护性</span>）：${coreGenes}</div>
        <div>扩展模块重叠基因：${geneList(x.ext.overlap, 60)}</div>
        ${x.known ? `<div>所选化疗药的已知相关不良事件：${esc(x.known.map(k => `${k[0]} — ${k[1]}${k[2] != null ? `（FAERS logLR ${k[2]}）` : ""}`).join("；"))}</div>` : ""}
      </td></tr>`);
  }
  out.push(`<section class="panel"><div class="sectionhead"><h2>不良反应模块（29 个）</h2><p class="hint">点击行展开重叠基因；“已知”表示所选化疗药在 FAERS 或规则库中有该类不良事件</p></div>
    <div class="tbl"><table><thead><tr><th>不良反应</th><th>化疗药</th><th>核心 重叠</th><th>核心 Z</th><th>核心 百分位</th><th>扩展 重叠</th><th>扩展 Z</th><th>扩展 百分位</th><th>网络信号</th></tr></thead><tbody>${aeRows.join("")}</tbody></table></div>
    <p class="hint">按特异性排序。百分位指与该模块网络距离不大于本输入的中药占全部中药的比例，越低越特异（≤ 10% 记为高特异，≤ 25% 较特异）。网络信号只说明中药靶点与该不良反应相关基因在网络上接近或重叠，不代表加重或减轻；方向需要结合基因功能（核心基因的促损伤/保护性标注）与临床证据判断。</p></section>`);

  out.push(`<section class="panel"><div class="sectionhead"><h2>机制网络</h2><p class="hint">展示与癌种、化疗靶点、PK 基因或特异（百分位 ≤ 25%）不良反应模块重叠的靶点；可拖动缩放</p></div>
    <div class="legend"><span><i class="dot" style="background:var(--herb)"></i>中药</span><span><i class="dot" style="background:#9bc8a8"></i>成分</span><span><i class="dot" style="background:var(--target)"></i>靶点</span><span><i class="dot" style="background:var(--drug)"></i>化疗药</span><span><i class="dot" style="background:#6b7b8c"></i>癌种</span><span><i class="dot" style="background:var(--module)"></i>不良反应模块</span></div>
    <div id="graph"></div></section>`);
  out.push(`<section class="panel" id="enrich-panel"><div class="sectionhead"><h2>通路富集</h2><p class="hint">Reactome 与 GO 生物过程；超几何检验，BH 校正，背景为 STRING 网络基因</p></div><div id="enrich-body"><p class="hint">正在载入通路库…</p></div></section>`);
  out.push(renderKB(res.kb));
  out.push(`<section class="panel"><h2>下载</h2><div class="downloads">
    <button class="btn ghost" data-dl="network">网络检验结果 CSV</button><button class="btn ghost" data-dl="targets">中药靶点–成分 CSV</button>
    <button class="btn ghost" data-dl="enrich">富集结果 CSV</button><button class="btn ghost" id="copy-link">复制本次分析链接</button></div></section>`);
  R.innerHTML = out.join("");
  R.querySelectorAll("tr.ae-row").forEach(tr => {
    const toggle = () => { const d = R.querySelector(`tr.detail[data-for="${tr.dataset.id}"]`); d.hidden = !d.hidden; };
    tr.onclick = toggle; tr.onkeydown = e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); toggle(); } };
  });
  R.querySelectorAll("[data-dl]").forEach(b => b.onclick = () => download(b.dataset.dl));
  $("copy-link").onclick = copyLink;
  drawGraph(res);
  renderEnrichment(res);
}

function topCompounds(res, gene, n) {
  const m = res.herb.geneComp.get(gene);
  if (!m) return "";
  return [...m.keys()].slice(0, n).map(ci => S.data.herbs.compounds[ci][1]).join("、") + (m.size > n ? ` 等 ${m.size} 个` : "");
}

function renderKB(kb) {
  if (!kb) return `<section class="panel"><h2>规则库机制推理</h2><p class="hint">所选中药或化疗药不在人工整理的药理规则库（约 20 种草药/制剂 × 20 种化疗药）中，本节跳过。</p></section>`;
  const lv = x => x >= .8 ? "高" : x >= .45 ? "中" : x > 0 ? "低" : "";
  return `<section class="panel"><div class="sectionhead"><h2>规则库机制推理（PK/PD）</h2><p class="hint">${esc(kb.herbs.join("、"))} × ${esc(kb.drugs.join("、"))}；来自人工整理的药理知识库，等级代表证据可信度而非效应大小；文献引用已于 2026-10-26 逐条核对 PubMed 并附 PMID</p></div>
    <div class="cols"><div class="tbl"><table><tr><th>可能的不良反应</th><th>等级</th><th>机制</th></tr>${kb.ae.map(r => `<tr><td><b>${esc(r.label)}</b></td>
      <td>${r.up ? `<span class="pill ${r.up >= .8 ? "sig" : "near"}">${lv(r.up)}</span>` : `<span class="pill known">可能减轻</span>`}</td>
      <td>${r.why.map(w => `<div>${w[0] === "increase" ? "↑" : "↓"} ${esc(w[1])}</div>`).join("")}</td></tr>`).join("") || `<tr><td colspan="3" class="hint">未推导出相互作用</td></tr>`}</table></div>
    <div class="tbl"><table><tr><th>中药/成分</th><th>靶点</th><th>证据</th><th>化疗药通路</th><th>推导</th></tr>${kb.pk.sort((a, b) => b.score - a.score).map(s => `<tr><td>${esc(s.src)}<div class="hint">${esc(s.m.constituent || "")}</div></td>
      <td class="gene">${s.m.effect === "inhibits" ? "抑制" : "诱导"} ${esc(s.m.target)}</td><td>${KB_EVL[s.m.evidence]}${s.m.ref ? `<div class="hint">${esc(s.m.ref)}</div>` : ""}</td>
      <td>${esc(s.drug)}：${KB_ROLE[s.p.role]}</td><td>${esc(s.cons)}</td></tr>`).join("") || `<tr><td colspan="5" class="hint">无 PK 通路交集</td></tr>`}</table></div></div>
    ${kb.observed.length ? `<p style="margin-top:8px"><b>已知人体研究：</b>${kb.observed.map(o => `${esc(o.herb)} + ${esc(o.drug)}（${esc(o.level)}）${esc(o.finding)} — ${esc(o.ref)}`).join("；")}</p>` : ""}</section>`;
}

function drawGraph(res) {
  const el = $("graph");
  if (!window.vis) { el.outerHTML = `<p class="hint">网络图组件未能加载，表格结果不受影响。</p>`; return; }
  const fg = getComputedStyle(document.documentElement).getPropertyValue("--fg").trim();
  const nodes = new Map(), edges = [];
  const node = (id, label, color, shape = "dot", size = 10) => { if (!nodes.has(id)) nodes.set(id, { id, label, color, shape, size, font: { color: shape === "box" ? "#fff" : fg, size: 12 } }); };
  const gene2ctx = new Map();
  const link = (g, ctxNode) => { if (!gene2ctx.has(g)) gene2ctx.set(g, new Set()); gene2ctx.get(g).add(ctxNode); };
  const cancerCtx = res.contexts.find(c => c.type === "CANCER");
  if (cancerCtx) cancerCtx.overlap.forEach(g => link(g, "c:" + S.cancer.id));
  for (const [g, ds] of res.pd) if (res.herb.geneComp.has(g)) ds.forEach(d => link(g, "d:" + d));
  res.pk.forEach(x => link(x.gene, "d:" + x.drug));
  const specific = c => c.sig && (c.prox.pct ?? 100) <= 25;
  res.contexts.filter(c => c.type === "AE_CORE" && (specific(c) || specific(res.contexts.find(e => e.type === "AE_EXTENDED" && e.id === c.id))))
    .forEach(c => c.overlap.forEach(g => link(g, "a:" + c.id)));
  const genes = [...gene2ctx.keys()].sort((a, b) => gene2ctx.get(b).size - gene2ctx.get(a).size || (S.G.deg[S.G.idx.get(b)] || 0) - (S.G.deg[S.G.idx.get(a)] || 0)).slice(0, 45);
  if (!genes.length) { el.outerHTML = `<p class="hint">没有与癌种、化疗药或显著不良反应模块重叠的靶点可画。</p>`; return; }
  S.herbs.forEach(i => node("h:" + S.data.herbs.herbs[i].cn, S.data.herbs.herbs[i].cn, "#2f8f5b", "box"));
  S.drugs.forEach(d => node("d:" + (d.cn || d.en), d.cn || d.en, "#c0392b", "box"));
  if (cancerCtx) node("c:" + S.cancer.id, S.cancer.cn || S.cancer.en, "#6b7b8c", "ellipse");
  const compCount = new Map();
  for (const g of genes) for (const ci of (res.herb.geneComp.get(g) || new Map()).keys()) compCount.set(ci, (compCount.get(ci) || 0) + 1);
  const comps = new Set([...compCount.entries()].sort((a, b) => b[1] - a[1]).slice(0, 25).map(x => x[0]));
  for (const g of genes) {
    node("g:" + g, g, "#7d6bc2", "dot", 9);
    for (const ctx of gene2ctx.get(g)) {
      if (ctx.startsWith("a:")) { const m = S.data.ae.modules.find(x => x.id === ctx.slice(2)); node(ctx, m.short || m.cn, "#e08e1b", "diamond", 12); }
      edges.push({ from: "g:" + g, to: ctx, color: { color: "#bbb" } });
    }
    for (const [ci, herbCn] of res.herb.geneComp.get(g) || []) {
      if (!comps.has(ci)) continue;
      const c = S.data.herbs.compounds[ci];
      node("m:" + ci, c[1].length > 22 ? c[1].slice(0, 20) + "…" : c[1], "#9bc8a8", "dot", 7);
      edges.push({ from: "m:" + ci, to: "g:" + g, color: { color: "#9bc8a8" } });
      edges.push({ from: "h:" + herbCn, to: "m:" + ci, color: { color: "#2f8f5b" } });
    }
  }
  const kept = [...new Map(edges.map(e => [e.from + ">" + e.to, e])).values()].filter(e => nodes.has(e.from) && nodes.has(e.to));
  const linked = new Set(kept.flatMap(e => [e.from, e.to]));
  const net = new vis.Network(el, { nodes: new vis.DataSet([...nodes.values()].filter(n => linked.has(n.id))), edges: new vis.DataSet(kept) },
    { edges: { smooth: false, width: 1 }, nodes: { borderWidth: 0 }, physics: { solver: "forceAtlas2Based", stabilization: { iterations: 250 } }, interaction: { hover: true } });
  net.once("stabilizationIterationsDone", () => { net.setOptions({ physics: false }); net.fit(); });
}

async function renderEnrichment(res) {
  const body = $("enrich-body");
  try {
    if (!S.pathways) { const r = await fetch("data/pathways.json"); if (!r.ok) throw new Error(`HTTP ${r.status}`); S.pathways = await r.json(); }
  } catch (e) { body.innerHTML = `<p class="hint">通路库载入失败：${esc(e.message)}</p>`; return; }
  const sets = [["中药全部靶点", res.A]];
  const cancerCtx = res.contexts.find(c => c.type === "CANCER"), pdCtx = res.contexts.find(c => c.type === "CHEMO_PD");
  if (cancerCtx) sets.push([`癌种模块邻域（${cancerCtx.label.replace("癌种：", "")}）`, nearGenes(res.A, cancerCtx.genes)]);
  if (pdCtx) sets.push(["化疗药靶点邻域", nearGenes(res.A, pdCtx.genes)]);
  res.enrich = sets.map(([label, idx]) => ({ label, n: idx.length, r: enrich(idx) }));
  body.innerHTML = `<div class="cols">${res.enrich.map(e => `<div><h3>${esc(e.label)}<span class="hint">（输入 ${e.n} 个基因）</span></h3>${!e.r ? `<p class="hint">基因少于 3 个，跳过。</p>` :
    Object.entries(e.r).map(([lib, v]) => `<div class="tbl"><table><tr><th>${lib === "GO_BP" ? "GO 生物过程" : "Reactome"}</th><th>重叠</th><th>q</th></tr>${v.top.slice(0, 8).map(t => `<tr title="${esc(t.genes.join(", "))}"><td>${esc(t.name)}</td><td class="num">${t.k}/${t.K}</td><td class="num">${fmtP(t.q)}</td></tr>`).join("") || `<tr><td colspan="3" class="hint">无可检验条目</td></tr>`}</table></div>`).join("")}</div>`).join("")}</div>
    <p class="hint">“邻域”指位于该模块内或与模块基因在 STRING 中直接相连的中药靶点，与 final27 分析一致。鼠标悬停可看重叠基因。</p>`;
}

/* ------------------------------------------------------------------ downloads & links */
function csv(rows) { return "﻿" + rows.map(r => r.map(x => { const s = String(x ?? ""); return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; }).join(",")).join("\n"); }
function download(kind) {
  const res = S.result; if (!res) return;
  let rows, name;
  if (kind === "network") {
    rows = [["context_type", "context_id", "label", "module_genes_in_network", "herb_genes_in_network", "overlap_n", "overlap_genes", "hypergeom_p", "overlap_q", "distance", "random_mean", "random_sd", "z", "empirical_p", "proximity_q", "one_hop_edges", "significant", "herb_background_percentile", "herb_background_n"]];
    for (const c of res.contexts) rows.push([c.type, c.id, c.label, c.nB, res.A.length, c.overlap.length, c.overlap.join(";"), c.pOverlap, c.qOverlap, c.prox.d, c.prox.mu, c.prox.sd, c.prox.z, c.prox.p, c.qProx, c.prox.oneHop, c.sig, c.prox.pct, c.prox.nBg]);
    name = "network_tests.csv";
  } else if (kind === "targets") {
    rows = [["gene", "in_string_network", "compounds", "herbs"]];
    for (const [g, m] of res.herb.geneComp) rows.push([g, S.G.idx.has(g), [...m.keys()].map(ci => S.data.herbs.compounds[ci][1]).join(";"), [...new Set(m.values())].join(";")]);
    name = "herb_targets.csv";
  } else {
    rows = [["gene_set", "library", "term", "overlap", "term_size", "p", "q", "genes"]];
    for (const e of res.enrich || []) for (const [lib, v] of Object.entries(e.r || {})) for (const t of v.top) rows.push([e.label, lib, t.name, t.k, t.K, t.p, t.q, t.genes.join(";")]);
    name = "enrichment_top.csv";
  }
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([csv(rows)], { type: "text/csv;charset=utf-8" }));
  a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
function stateHash() {
  const p = new URLSearchParams();
  p.set("h", S.herbs.map(i => S.data.herbs.herbs[i].cn).join(","));
  p.set("d", S.drugs.map(d => d.live ? d.id : d.en).join(","));
  if (S.cancer) p.set("c", S.cancer.id);
  if (!$("opt-obdl").checked) p.set("obdl", "0");
  if ($("opt-validated").checked) p.set("v", "1");
  return "#" + p.toString();
}
async function copyLink() {
  const url = location.href.split("#")[0] + stateHash();
  try { await navigator.clipboard.writeText(url); status("已复制分析链接"); } catch { prompt("复制下面的链接：", url); }
}
async function restoreFromHash() {
  if (location.hash.length < 3) return false;
  const p = new URLSearchParams(location.hash.slice(1));
  const herbs = (p.get("h") || "").split(",").filter(Boolean);
  S.herbs = herbs.map(cn => S.data.herbs.herbs.findIndex(h => h.cn === cn)).filter(i => i >= 0);
  S.drugs = [];
  for (const k of (p.get("d") || "").split(",").filter(Boolean)) {
    if (S.drugByEn.has(k)) S.drugs.push(S.drugByEn.get(k));
    else if (/^CHEMBL/.test(k)) { try { S.drugs.push(await otDrug(k)); } catch { /* skip */ } }
  }
  const cid = p.get("c");
  S.cancer = cid ? (S.data.cancers.cancers.find(c => c.id === cid) || await otDisease(cid).catch(() => null)) : null;
  $("opt-obdl").checked = p.get("obdl") !== "0";
  $("opt-validated").checked = p.get("v") === "1";
  renderChips();
  return S.herbs.length > 0;
}

/* ------------------------------------------------------------------ run */
async function run() {
  if (!S.herbs.length) { status("请至少选择一味中药。", true); return; }
  if (!S.drugs.length && !S.cancer) { status("请至少选择化疗药或癌种之一。", true); return; }
  const btn = $("run"); btn.disabled = true;
  try {
    const opt = { obdl: $("opt-obdl").checked, validated: $("opt-validated").checked, fast: $("opt-fast").checked };
    status("计算中…"); await sleep0();
    const t0 = performance.now();
    S.result = await analyse(opt, status);
    render(S.result);
    history.replaceState(null, "", stateHash());
    status(`完成，用时 ${((performance.now() - t0) / 1000).toFixed(1)} 秒。`);
  } catch (e) { console.error(e); status(`分析出错：${e.message}`, true); }
  finally { btn.disabled = false; }
}
function useExample(ex) {
  S.herbs = ex.herbs.map(cn => S.data.herbs.herbs.findIndex(h => h.cn === cn)).filter(i => i >= 0);
  S.drugs = [];
  for (const name of ex.drugs) { if (REGIMENS[name]) REGIMENS[name].forEach(e => S.drugByEn.get(e) && S.drugs.push(S.drugByEn.get(e)));
    else { const d = S.data.drugs.drugs.find(x => x.cn === name || x.en === name); if (d) S.drugs.push(d); } }
  S.cancer = S.data.cancers.cancers.find(c => c.cn === ex.cancer) || null;
  renderChips(); run();
}
function renderAbout() {
  const H = S.data.herbs, nComp = H.compounds.length;
  $("about-body").innerHTML = `
  <p><b>数据</b></p><ul>
    <li>中药–成分–靶点：TCMSP 2.3（old.tcmsp-e.com）全部 ${H.herbs.length} 味药、${nComp.toLocaleString()} 个成分；默认按 OB ≥ 30%、DL ≥ 0.18 筛选口服活性成分；靶点含 TCMSP 实验验证与预测两类（可切换为只用验证靶点）。TCMSP 只给蛋白名，已通过 UniProt（人类、已审阅条目、名称完全匹配）映射为基因符号，非人类靶点不计。</li>
    <li>蛋白互作网络：STRING v12.0 人类，综合得分 ≥ 900（${S.G.n.toLocaleString()} 个基因、${(S.data.network.edges.length / 2).toLocaleString()} 条边），同时作为统计背景。</li>
    <li>癌种基因模块：Open Targets Platform 疾病关联靶点前 300 名（内置 ${S.data.cancers.cancers.length} 个常见癌种；其他疾病联网实时查询）。</li>
    <li>化疗药：${S.data.drugs.drugs.length} 种常用抗肿瘤药。药效靶点优先采用 final27 分析中整理的核心基因，其次为 Open Targets 作用机制靶点；铂类、多数烷化剂等直接作用于 DNA 的药物没有蛋白靶点，按同类药物代理或标注为不可用。代谢/转运基因来自 final27 表和人工规则库。不良事件来自 Open Targets 收录的 FAERS 显著信号。</li>
    <li>不良反应模块：final27 分析的 29 个模块，核心基因为人工整理的机制基因（标注促损伤/保护性），扩展基因来自 HPO 表型–基因关联。</li>
    <li>通路：MSigDB v2026.1 中的 Reactome 与 GO 生物过程（KEGG 因许可限制未纳入）。</li></ul>
  <p><b>方法</b>（与 final27 分析一致）：中药靶点与各模块的重叠用超几何检验；网络邻近度用 Guney 等的 closest 距离（每个中药靶点到最近模块基因的最短路径均值），以度匹配随机抽样（度窗口逐步放宽直到候选 ≥ 20 个节点）重复 1,000 次得到 Z 值和经验 P 值；一次查询内的全部检验分别对重叠 P 和邻近度 P 做 BH 校正；重叠 > 0 且 q < 0.05，或 Z < −2 且 q < 0.05 视为显著。富集只使用位于模块内或与模块基因直接相连的中药靶点。随机数种子固定，同一输入结果可复现，但与论文中 Python 实现的随机序列不同，Z 值会有小幅差异。</p>
  <p><b>特异性（本网页新增，不在 final27 中）：</b>TCMSP 靶点普遍偏向网络枢纽基因，度匹配随机化无法完全消除这种偏倚，因此大多数中药与多数模块都会“显著接近”。为此，每个模块再计算“中药间百分位”：在同样筛选条件下、网络中靶点 ≥ 5 个的全部单味中药里，与该模块的 closest 距离不大于本输入的比例（不含输入本身）。百分位 ≤ 10% 记为高特异，≤ 25% 记为较特异。多味复方是与单味药比较，仅作参考。</p>
  <p class="note"><b>使用须知：</b>本工具用于科研假设生成。网络关联不等于临床作用，不能区分抑制或诱导、加重或减轻；TCMSP 预测靶点存在假阳性。结果不能作为用药依据，临床决策请咨询医生和临床药师。</p>
  <p class="hint">主要参考：Ru J 等，J Cheminform 2014（TCMSP）；Szklarczyk D 等，Nucleic Acids Res 2023（STRING v12）；Ochoa D 等，Nucleic Acids Res 2023（Open Targets）；Guney E 等，Nat Commun 2016（网络邻近度）；Liberzon A 等，Cell Syst 2015（MSigDB）；Köhler S 等，Nucleic Acids Res 2021（HPO）。</p>`;
}

async function main() {
  status("正在载入数据…");
  try { await loadData(); } catch (e) { status(`数据载入失败：${e.message}`, true); return; }
  ["herb", "drug", "cancer"].forEach(setupPicker);
  const ex = $("examples");
  for (const e of EXAMPLES) { const b = document.createElement("button"); b.type = "button"; b.textContent = e.label; b.onclick = () => useExample(e); ex.appendChild(b); }
  $("run").onclick = run;
  renderAbout();
  status(`已载入 ${S.data.herbs.herbs.length} 味中药、${S.data.drugs.drugs.length} 种化疗药、${S.data.cancers.cancers.length} 个癌种。`);
  if (await restoreFromHash()) run(); else useExample(EXAMPLES[0]);
}
main();
