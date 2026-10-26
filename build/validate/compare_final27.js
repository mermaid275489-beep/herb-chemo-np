const fs=require('fs'); const site=require('path').resolve(__dirname,'../../docs')+'/';
const app=fs.readFileSync(site+'app.js','utf8');
const pick=(a,b)=>app.slice(app.indexOf(a),app.indexOf(b));
eval(pick('const SEED_BASE','const S =')+pick('function buildGraph','/* ------------------------------------------------------------------ statistics */')+pick('/* ------------------------------------------------------------------ statistics */','/* ------------------------------------------------------------------ query sets */').replace(/^let LF/m,'var LF'));
const G=buildGraph(JSON.parse(fs.readFileSync(site+'data/network.json','utf8')));
const ae=JSON.parse(fs.readFileSync(site+'data/ae.json','utf8')).modules;
const V=JSON.parse(fs.readFileSync(process.argv[2],'utf8'));
const idx=gs=>[...new Set(gs)].map(g=>G.idx.get(g)).filter(i=>i!==undefined).sort((a,b)=>a-b);
const A=idx(V.genes); console.log('A in network', A.length);
const sets={CHEMO_PD:{'实际化疗PD核心':V.pd}, CANCER:new Proxy({},{get:()=>V.cancer}), AE_CORE:{}, AE_EXTENDED:{}};
for(const m of ae){const core=m.core.map(x=>x[0]); sets.AE_CORE[m.id]=core; sets.AE_EXTENDED[m.id]=[...new Set([...core,...m.extended])];}
const rows=[];
for(const r of V.ref){const B=idx(sets[r.Context_type][r.Context_ID]); const Bs=new Set(B);
  const ov=A.filter(a=>Bs.has(a)).length; const p=ov?hyperSF(ov,G.n,B.length,A.length):1;
  const pr=proximity(G,A,B,'val|'+r.Context_type+'|'+r.Context_ID,1000);
  rows.push([r.Context_type,r.Context_ID,`B ${B.length}/${r.Set_B_n}`,`ov ${ov}/${r.Overlap_n}`,`p ${p.toExponential(2)}/${r.Hypergeometric_P.toExponential(2)}`,`d ${pr.d.toFixed(4)}/${r.Distance.toFixed(4)}`,`mu ${pr.mu.toFixed(3)}/${r.Random_mean.toFixed(3)}`,`sd ${pr.sd.toFixed(3)}/${r.Random_SD.toFixed(3)}`,`Z ${pr.z.toFixed(2)}/${r.Z.toFixed(2)}`].join('  '));}
console.log(rows.join('\n'));
