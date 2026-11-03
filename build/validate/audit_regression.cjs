/* Regression tests for the defects found in the 2026-10-26 site audit.
 * Each case is a counter-example from that report; the test runs the page's own JavaScript
 * engine (sliced out of docs/app.js) so a fix cannot pass here while the site still misbehaves.
 *
 *   node build/validate/audit_regression.cjs
 */
const fs = require('fs');
const path = require('path');

const DOCS = path.resolve(__dirname, '../../docs');
const app = fs.readFileSync(path.join(DOCS, 'app.js'), 'utf8');
const slice = (a, b) => {
  const i = app.indexOf(a), j = app.indexOf(b);
  if (i < 0 || j < 0 || j <= i) throw new Error(`cannot slice app.js between ${a} and ${b}`);
  return app.slice(i, j);
};
const read = f => JSON.parse(fs.readFileSync(path.join(DOCS, 'data', f), 'utf8'));

/* globals the sliced code expects */
var S = { data: {}, herbs: [], drugs: [], cancer: null };
const document = { visibilityState: 'visible' };

eval(slice('const SEED_BASE', 'const S =')
  + slice('function buildGraph', '/* ------------------------------------------------------------------ statistics */')
  + slice('/* ------------------------------------------------------------------ statistics */',
          '/* ------------------------------------------------------------------ query sets */').replace(/^let LF/m, 'var LF')
  + slice('const toIdx =', '/* ------------------------------------------------------------------ analysis */')
  + slice('/* near genes:', '/* ------------------------------------------------------------------ rule-based KB')
  + slice('/* ------------------------------------------------------------------ rule-based KB',
          '/* ------------------------------------------------------------------ live Open Targets'));

S.data = { herbs: read('herbs.json'), drugs: read('drugs.json'), kb: read('kb.json'), network: read('network.json') };
S.pathways = read('pathways.json');
S.G = buildGraph(S.data.network);
const drugByEn = new Map(S.data.drugs.drugs.map(d => [d.en, d]));
const herbIdx = cn => S.data.herbs.herbs.findIndex(h => h.cn === cn);

let failed = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${name}${ok ? '' : `\n        ${detail}`}`);
  if (!ok) failed++;
};
const run = (herbNames, drugEn) => {
  S.herbs = herbNames.map(herbIdx).filter(i => i >= 0);
  S.drugs = [drugByEn.get(drugEn)].filter(Boolean);
  S.cancer = null;
  return kbPredict();
};

/* ---- issue 1: BH must span the whole pathway library, not only terms with >= 3 overlap.
   The audit's counter-example is the 5-gene chemo-target neighbourhood of the default query. */
{
  const genes = ['DDB1', 'F2', 'PCNA', 'RASSF1', 'TOP2A'];
  const idx = genes.map(g => S.G.idx.get(g)).filter(i => i !== undefined).sort((a, b) => a - b);
  check('issue 1: counter-example gene set resolves', idx.length === 5, `got ${idx.length}`);
  const r = enrich(idx);
  const go = r.GO_BP, re = r.Reactome;
  check('issue 1: BH family is the full library',
    go.nTests === S.pathways.libraries.GO_BP.length && re.nTests === S.pathways.libraries.Reactome.length,
    `GO ${go.nTests}/${S.pathways.libraries.GO_BP.length}, Reactome ${re.nTests}/${S.pathways.libraries.Reactome.length}`);
  const ddr = go.top.find(t => /dna damage response/i.test(t.name));
  check('issue 1: "DNA damage response" q matches the independent check (~0.1624)',
    ddr && Math.abs(ddr.q - 0.1624325129616696) < 1e-6, ddr ? `q=${ddr.q}` : 'term not shown');
  const ptm = re.top.find(t => /post translational protein modification/i.test(t.name));
  check('issue 1: "Post translational protein modification" q matches (~0.1023)',
    ptm && Math.abs(ptm.q - 0.10228859744898414) < 1e-6, ptm ? `q=${ptm.q}` : 'term not shown');
  check('issue 1: no term in this 5-gene example survives q < 0.05',
    go.nSig === 0 && re.nSig === 0, `GO ${go.nSig}, Reactome ${re.nSig}`);
}

/* ---- issue 2: curated evidence must be reachable from the names the page actually offers */
for (const [cn, expect] of [['贯叶金丝桃', 'st_johns_wort'], ['南五味子', 'schisandra_sphenanthera'],
                            ['北五味子', 'schisandra_chinensis']]) {
  const kb = run([cn], 'irinotecan');
  check(`issue 2: ${cn} reaches the rule KB as ${expect}`,
    !kb.empty && kb.herbs.includes(expect), kb.empty ? 'kbPredict returned empty' : kb.herbs.join(','));
}
{
  /* the two Schisandra species carry different human evidence and must not be merged */
  const nan = run(['南五味子'], 'irinotecan'), bei = run(['北五味子'], 'irinotecan');
  const top = k => (k.ae[0] || {}).up || 0;
  check('issue 2: S. sphenanthera (CYP3A4/tacrolimus) ranks above S. chinensis (P-gp only)',
    top(nan) > top(bei) + 0.5, `sphenanthera ${top(nan)}, chinensis ${top(bei)}`);
}

/* ---- issue 3: a decoction piece must not inherit an injection's profile */
{
  const kb = run(['薏苡仁'], 'paclitaxel');
  check('issue 3: 薏苡仁 does not inherit Kanglaite injection',
    kb.empty && (kb.uncoveredHerbs || []).includes('薏苡仁'),
    kb.empty ? 'ok but name missing from uncovered list' : `matched ${kb.herbs.join(',')}`);
}

/* ---- issue 4: opposite directions are condition-dependent, not cancelling */
{
  const kb = run(['丹参'], 'docetaxel');
  const c = (kb.conflicts || []).find(x => x.targets.includes('CYP3A4'));
  check('issue 4: danshen CYP3A4 conflict is reported', !!c, JSON.stringify(kb.conflicts));
  check('issue 4: both branches are kept', c && c.branches.length >= 2,
    c ? `${c.branches.length} branch(es)` : 'no conflict');
  check('issue 4: no net-exposure AE is emitted for the conflicting pair',
    !kb.ae.some(a => a.why.some(w => /暴露[↑↓]/.test(w[1]))),
    kb.ae.map(a => a.label).join(','));
}

/* ---- issue 5: a null human study of the same pair constrains the in vitro extrapolation */
{
  const kb = run(['水飞蓟'], 'irinotecan');
  const worst = Math.max(0, ...kb.ae.map(a => a.up));
  check('issue 5: milk thistle x irinotecan no longer ranks as moderate', worst <= 0.35,
    `highest up = ${worst}`);
  check('issue 5: the human-null caveat is surfaced',
    kb.pk.some(s => (s.caveats || []).some(c => /未检出/.test(c))),
    JSON.stringify(kb.pk.map(s => s.caveats)));
}

/* ---- issue 6: inhibition and induction are not symmetrically consequential */
{
  const kb = run(['黄连'], 'osimertinib');
  const worst = Math.max(0, ...kb.ae.map(a => a.up));
  check('issue 6: berberine x osimertinib is no longer ranked high', worst < 0.45, `highest up = ${worst}`);
  check('issue 6: the asymmetry is explained',
    kb.pk.some(s => (s.caveats || []).some(c => /不对称/.test(c))),
    JSON.stringify(kb.pk.map(s => s.caveats)));
  check('issue 6: probe-drug extrapolation is labelled',
    kb.pk.some(s => (s.caveats || []).some(c => /探针药/.test(c))), 'no probe caveat');
}

console.log(failed ? `\n${failed} check(s) FAILED` : '\nall checks passed');
process.exit(failed ? 1 : 0);
