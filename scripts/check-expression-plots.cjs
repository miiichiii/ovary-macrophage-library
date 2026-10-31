// Mathematical fixtures + actual gzip gene vectors and VM interaction contracts. Not a visual test.
const fs=require('node:fs'),vm=require('node:vm'),path=require('node:path'),assert=require('node:assert/strict'),zlib=require('node:zlib');
const root=path.resolve(__dirname,'..'),viewer=path.join(root,'studies/morris-2022/viewer'),atlas=JSON.parse(fs.readFileSync(path.join(viewer,'manifest.json')));
const api=require(path.join(viewer,'expression-plots.js')),checks=[];
function check(name,fn){fn();checks.push(name)}
function near(a,b,tol=1e-12){assert(Math.abs(a-b)<=tol,`${a} != ${b}`)}
check('means include zero counts; odd/even medians and detection fraction',()=>{
 const a=api.math.summarize([0,0,2,6]);assert.equal(a.n,4);assert.equal(a.mean,2);assert.equal(a.median,1);assert.equal(a.positive,2);assert.equal(api.math.summarize([0,2,6]).median,2);
});
check('empty groups are null, not zero; invalid finite/nonnegative values rejected',()=>{
 const a=api.math.summarize([]);assert.equal(a.n,0);assert.equal(a.mean,null);assert.equal(a.median,null);assert.equal(api.math.summarize([0,0]).mean,0);assert.throws(()=>api.math.summarize([NaN]));assert.throws(()=>api.math.summarize([-1]));
});
check('normalization uses full-gene totals and target sum, preserving zero',()=>{
 const a=api.math.normalize([0,3,5],[2,0,1],[20,30,100],10000);assert.equal(a[0],Math.log1p(500));assert.equal(a[1],0);assert.equal(a[2],Math.log1p(1000));assert.throws(()=>api.math.normalize([0],[0],[0],10000));
});
check('row Z uses population SD and preserves nulls and constant rows',()=>{
 const z=api.math.rowZscores([[1,2,3,null],[0,0,null],[null,null]]);near(z[0][0],-Math.sqrt(1.5));near(z[0][1],0);near(z[0][2],Math.sqrt(1.5));assert.equal(z[0][3],null);assert.deepEqual(z[1],[0,0,null]);assert.deepEqual(z[2],[null,null]);
});
check('bounded distribution retains zero observations and honest one-cell/constant groups',()=>{
 assert.equal(api.math.density([],5),null);assert.equal(api.math.density([0],5).constant,0);assert.equal(api.math.density([2,2,2],5).constant,2);
 const d=api.math.density([0,0,0,1],1),one=api.math.density([0,1,1,1],1);assert.equal(d.x.length,80);assert.equal(d.width.length,80);assert(d.width[0]>d.width.at(-1));assert(one.width[0]<one.width.at(-1));assert.equal(d.x[0],0);assert.equal(d.x.at(-1),1);assert(d.width.every(v=>Number.isFinite(v)&&v>=0&&v<=.38));
});
class Element{constructor(tag='div'){this.tagName=tag;this.children=[];this.dataset={};this.value='';this.style={};this.textContent='';this.classList={toggle(){}}}setAttribute(k,v){this[k]=v;if(k.startsWith('data-'))this.dataset[k.slice(5)]=v}append(...xs){this.children.push(...xs);if(this.tagName==='select'&&this.children.length)this.value=this.children[0].value}replaceChildren(...xs){this.children=[];this.append(...xs)}addEventListener(k,fn){this['on'+k]=fn}scrollIntoView(){}querySelectorAll(){return[]}on(k,fn){this[k]=fn}remove(){}click(){}}
const els={},purged=[];const document={documentElement:{},head:new Element(),createElement:t=>new Element(t),querySelector:s=>els[s]??=new Element(),querySelectorAll:()=>[],getElementById:()=>new Element()};
const Plotly={newPlot:(el,traces,layout)=>{el.last={traces,layout};return Promise.resolve()},react:(el,traces,layout)=>{el.last={traces,layout};return Promise.resolve()},restyle(){},Plots:{resize(){}},purge(plot){purged.push(plot)}};
const context={window:{OVARY:atlas},document,location:{search:'?lang=en'},localStorage:{getItem(){return null},setItem(){}},URLSearchParams,Plotly,matchMedia:()=>({matches:true}),setTimeout,clearTimeout,atob,Blob,Response,DecompressionStream};vm.createContext(context);
vm.runInContext(fs.readFileSync(path.join(viewer,'cluster-labels.js'),'utf8'),context);
vm.runInContext(fs.readFileSync(path.join(viewer,'expression-plots.js'),'utf8'),context);
vm.runInContext(fs.readFileSync(path.join(viewer,'viewer.js'),'utf8'),context);
const qc=vm.runInContext('state.get(atlas.qc_run)',context),p=qc.expressionPlots;
function rawGene(symbol){const gene=atlas.features.findIndex(f=>f.symbol===symbol);let entries;vm.runInNewContext(fs.readFileSync(path.join(viewer,'genes',Math.floor(gene/atlas.gene_chunk_size)+'.js'),'utf8'),{window:{registerGenes:(c,e)=>entries=e},registerGenes:(c,e)=>entries=e});const record=entries[gene],buf=zlib.gunzipSync(Buffer.from(record.data,'base64')),counts=new Float64Array(atlas.n_cells);for(let i=0;i<record.n;i++)counts[buf.readUInt32LE(i*4)]=buf.readUInt32LE((record.n+i)*4);return{gene,counts}}
const mafb=rawGene('Mafb'),c1qa=rawGene('C1qa');
function filter(stage='',source='analysis',direct=false){qc.stageSelect.value=stage;qc.stageSource.value=source;qc.evidenceSelect.value=direct?'direct':'all';qc.stageSelect.onchange()}
(async()=>{
 // Feed actual gzip chunks through the production browser decoder, not just synthetic arrays.
 for(const {gene}of[mafb,c1qa])vm.runInContext(fs.readFileSync(path.join(viewer,'genes',Math.floor(gene/atlas.gene_chunk_size)+'.js'),'utf8'),context);
 p.input.value='Mafb, C1qa';await context.window.ExpressionPlots.load(p);
 check('real on-demand gzip decoding + two gene vectors render seven QC clusters',()=>{
  assert.equal(p.error.textContent,'');assert.equal(p.models.length,2);assert.equal(p.groups.length,7);assert.equal(p.groups.reduce((n,g)=>n+g.indices.length,0),1377);assert.equal(p.plot.hidden,false);assert(p.rendered.traces.length>7);
  p.models[0].values.forEach((v,i)=>near(v,Math.log1p(mafb.counts[qc.run.rows[i]]/atlas.totals[qc.run.rows[i]]*qc.run.recipe.target_sum)));
 });
 check('real summaries equal independent zero-inclusive count normalization by cluster',()=>{
  for(const g of p.models[0].groups){const indices=qc.run.clusters.flatMap((c,i)=>String(c)===g.key?[i]:[]),vals=indices.map(i=>Math.log1p(mafb.counts[qc.run.rows[i]]/atlas.totals[qc.run.rows[i]]*qc.run.recipe.target_sum));assert.equal(g.summary.n,indices.length);near(g.summary.mean,vals.reduce((a,b)=>a+b,0)/vals.length);assert.equal(g.summary.positive,vals.filter(v=>v>0).length)}
 });
 check('stage grouping uses current author/inferred labels and direct evidence only',()=>{
  p.group.value='stage';context.window.ExpressionPlots.render(p);assert.deepEqual(Array.from(p.groups,g=>g.indices.length),[317,486,193,381,0]);filter('','author');assert.deepEqual(Array.from(p.groups,g=>g.indices.length),[217,266,161,283,450]);filter('estrus','analysis',true);assert.deepEqual(Array.from(p.groups,g=>g.indices.length),[0,266,0,0,0]);
 });
 check('single-cluster scope intersects cycle and exact lasso row IDs',()=>{
  filter();p.scope.value=String(qc.run.clusters[0]);p.selected.checked=true;qc.selected=qc.run.rows.slice(0,30);vm.runInContext('selectionText(state.get(atlas.qc_run));window.ExpressionPlots.refresh(state.get(atlas.qc_run))',context);const expected=qc.run.rows.slice(0,30).filter((row,i)=>String(qc.run.clusters[i])===p.scope.value).length;assert.equal(p.groups.reduce((n,g)=>n+g.indices.length,0),expected);assert(expected>0);
 });
 check('empty lasso scope does not silently fall back to all cells',()=>{
  qc.selected=[];vm.runInContext('selectionText(state.get(atlas.qc_run));window.ExpressionPlots.refresh(state.get(atlas.qc_run))',context);assert.equal(p.plot.hidden,true);assert.equal(p.export.disabled,true);assert(p.status.textContent.includes('No lasso'));p.selected.checked=false;p.scope.value='';
 });
 check('heatmap means, Z-score, null cells and readable annotations',()=>{
  filter('estrus','analysis',true);p.kind.value='heatmap';context.window.ExpressionPlots.render(p);const tr=p.rendered.traces[0];assert.equal(tr.type,'heatmap');assert.equal(tr.z.length,2);assert.equal(tr.z[0][0],null);assert.equal(tr.z[0][1],p.models[0].groups[1].summary.mean);assert(p.rendered.layout.annotations.length===8);p.scale.value='zscore';context.window.ExpressionPlots.render(p);assert.equal(p.rendered.traces[0].z[0][1],0);assert.equal(p.rendered.traces[0].z[0][0],null);assert.equal(p.rendered.traces[0].customdata[0][1][1],266);assert.equal(p.rendered.traces[0].customdata[0][1][4],'Mafb');
 });
 check('CSV records scientific values, empty groups, and filter provenance',()=>{
  const old=p.a.downloadCSV;let output;p.a.downloadCSV=(name,rows)=>output=rows;p.export.onclick();p.a.downloadCSV=old;assert.equal(output.length,11);assert.equal(output[0][8],'mean_lognorm');assert.equal(output[1][6],0);assert.equal(output[1][8],'');assert.equal(output[2][6],266);assert.equal(output[2][11],'estrus');assert.equal(output[2][12],'analysis');assert.equal(output[2][13],true);
 });
 check('language refresh updates controls, stage names, plots, and cell-type labels',()=>{
  filter();p.group.value='cluster';vm.runInContext('language="ja";refreshLanguage()',context);assert.equal(p.submit.textContent,'図を表示');assert(p.scope.children[1].textContent.includes('クラスター'));assert(p.rendered.traces[0].customdata[0][0][0].includes('クラスター'));vm.runInContext('language="en";refreshLanguage()',context);assert.equal(p.submit.textContent,'Draw plots');
 });
 filter();
 for(const run of atlas.runs){
  context.testRunId=run.id;vm.runInContext('openRun(testRunId)',context);const other=vm.runInContext('state.get(testRunId)',context),op=other.expressionPlots;op.input.value='Mafb';await context.window.ExpressionPlots.load(op);
  assert.equal(op.error.textContent,'');assert.equal(op.groups.reduce((n,g)=>n+g.indices.length,0),run.n);
  const weighted=op.models[0].groups.reduce((sum,g)=>sum+(g.summary.mean??0)*g.summary.n,0)/run.n;
  const independent=run.rows.reduce((sum,row)=>sum+Math.log1p(mafb.counts[row]/atlas.totals[row]*run.recipe.target_sum),0)/run.n;near(weighted,independent,1e-10);
 }
 checks.push('all 14 saved analyses render real Mafb and weighted means match independent full-row normalization');
 p.input.value='not-a-gene';await context.window.ExpressionPlots.load(p);
 check('invalid gene errors are visible with no stale export or plot',()=>{assert(p.error.textContent.includes('not-a-gene'));assert.equal(p.export.disabled,true);assert.equal(p.plot.hidden,true)});
 p.input.value=Array.from({length:13},(_,i)=>'x'+i).join(',');await context.window.ExpressionPlots.load(p);
 check('request bounded to 12 gene queries',()=>assert(p.error.textContent.includes('12')));
 const original=p.a.geneCounts;let finish;
 p.a.geneCounts=()=>new Promise(resolve=>finish=resolve);p.input.value='Mafb';const pending=context.window.ExpressionPlots.load(p);
 p.a.geneCounts=original;p.input.value='C1qa';await context.window.ExpressionPlots.load(p);finish(mafb.counts);await pending;
 check('late result from older query cannot overwrite newer query',()=>{assert.equal(p.models.length,1);assert.equal(p.models[0].feature.symbol,'C1qa');assert.equal(p.loading,false)});
 let finishFilter;p.a.geneCounts=()=>new Promise(resolve=>finishFilter=resolve);p.input.value='Mafb';const pendingFilter=context.window.ExpressionPlots.load(p);filter('proestrus','analysis',true);finishFilter(mafb.counts);await pendingFilter;p.a.geneCounts=original;
 check('filter changed during load is applied to newly loaded plots',()=>assert.equal(p.groups.reduce((n,g)=>n+g.indices.length,0),217));
 let fail;p.a.geneCounts=()=>new Promise((_,reject)=>fail=reject);p.input.value='Mafb';const failed=context.window.ExpressionPlots.load(p);fail(Error('Network fixture'));await failed;p.a.geneCounts=original;
 check('network failure is visible and hides prior results/export',()=>{assert.equal(p.error.textContent,'Network fixture');assert.equal(p.plot.hidden,true);assert.equal(p.export.disabled,true)});
 p.input.value='Mafb';await context.window.ExpressionPlots.load(p);
 const doomed=p.plot;els['#reset'].onclick();
 check('reset disposes plot state and purges charts',()=>{assert(p.disposed);assert(purged.includes(doomed));assert.equal(p.models,null)});
 const result={type:'math + production gzip decoder + non-browser interaction contracts',passed:checks.length,checks};fs.writeFileSync(path.join(root,'docs/EXPRESSION_PLOTS_CHECKS.json'),JSON.stringify(result,null,2));console.log(JSON.stringify(result));
})().catch(e=>{console.error(e);process.exitCode=1});
