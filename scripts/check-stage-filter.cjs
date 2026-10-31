// Data and interaction-contract verification using DOM stubs; not a browser-render test.
const fs=require('node:fs'),vm=require('node:vm'),path=require('node:path'),assert=require('node:assert/strict'),zlib=require('node:zlib');
const root=path.resolve(__dirname,'..'),viewer=path.join(root,'studies/morris-2022/viewer'),atlas=JSON.parse(fs.readFileSync(path.join(viewer,'manifest.json')));
class Element{constructor(tag='div'){this.tagName=tag;this.children=[];this.dataset={};this.value='';this.textContent='';this.classList={toggle(){}}}setAttribute(k,v){this[k]=v;if(k.startsWith('data-'))this.dataset[k.slice(5)]=v}append(...xs){this.children.push(...xs);if(this.tagName==='select'&&this.children.length)this.value=this.children[0].value}replaceChildren(...xs){this.children=xs}addEventListener(k,fn){this['on'+k]=fn}scrollIntoView(){}querySelectorAll(){return[]}on(){}remove(){}click(){}}
const els={},plots=[],purged=[];const document={documentElement:{},head:new Element(),createElement:t=>new Element(t),querySelector:s=>els[s]??=new Element(),querySelectorAll:()=>[],getElementById:()=>new Element()};const Plotly={newPlot:(el,traces,layout)=>{el.last={traces,layout};plots.push(el);return Promise.resolve()},react:(el,traces,layout)=>{el.last={traces,layout};return Promise.resolve()},restyle(){},Plots:{resize(){}},purge(plot){purged.push(plot)}};
const context={window:{OVARY:atlas},document,location:{search:'?lang=en'},localStorage:{getItem(){return null},setItem(){}},URLSearchParams,Plotly,matchMedia:()=>({matches:true}),setTimeout,clearTimeout};vm.createContext(context);vm.runInContext(fs.readFileSync(path.join(viewer,'viewer.js'),'utf8'),context);
const run=atlas.runs.find(r=>r.id===atlas.qc_run),qc=vm.runInContext('state.get(atlas.qc_run)',context);const checks=[];
function check(name,fn){fn();checks.push(name)}
function filter(stage='',source='analysis',direct=false){qc.stageSelect.value=stage;qc.stageSource.value=source;qc.evidenceSelect.value=direct?'direct':'all';qc.stageSelect.onchange();return qc.plot.last.traces.flatMap(t=>t.customdata.map((d,i)=>({row:d[2],x:t.x[i],y:t.y[i]})))}
check('QC default candidate count / evidence',()=>{assert.equal(filter().length,1377);assert(qc.filterCount.textContent.includes('927 direct'));assert(qc.filterCount.textContent.includes('450 inferred'))});
for(const [stage,n,directN]of[['proestrus',317,217],['estrus',486,266],['metestrus',193,161],['diestrus',381,283]]){
 check(stage+' analysis stage count + exact saved coordinates',()=>{const points=filter(stage);assert.equal(points.length,n);for(const p of points){const i=run.rows.indexOf(p.row);assert.equal(run.metadata.analysis_stage[i],stage);assert.equal(p.x,run.coordinates[i][0]);assert.equal(p.y,run.coordinates[i][1])}});
 check(stage+' direct-only count',()=>{assert.equal(filter(stage,'analysis',true).length,directN)});
}
check('author stage retains450 unknown; direct-only excludes inferred',()=>{assert.equal(filter('unknown','author').length,450);assert.equal(filter('','analysis',true).length,927);assert.equal(filter('unknown','author',true).length,0);assert(qc.empty.textContent.includes('No cells'))});
// Load one real gene chunk and reconstruct raw counts independently using Node gzip.
const gene=atlas.features.findIndex(f=>f.symbol==='Mafb');assert(gene>=0);let entries;const chunk=Math.floor(gene/atlas.gene_chunk_size);vm.runInNewContext(fs.readFileSync(path.join(viewer,'genes',chunk+'.js'),'utf8'),{window:{registerGenes:(c,e)=>entries=e},registerGenes:(c,e)=>entries=e});const record=entries[gene],buf=zlib.gunzipSync(Buffer.from(record.data,'base64')),counts=new Float64Array(atlas.n_cells);for(let i=0;i<record.n;i++)counts[buf.readUInt32LE(i*4)]=buf.readUInt32LE((record.n+i)*4);
context.testGene=gene;context.testCounts=counts;filter('estrus');vm.runInContext('addExpression(state.get(atlas.qc_run),testGene,testCounts)',context);const expression=qc.expressions.get(gene);expression.colors.value='Magma';expression.range=[0,3];
check('Mafb expression aligns486 estrus cells with exact raw-count normalization',()=>{const tr=expression.plot.last.traces[0];assert.equal(tr.x.length,486);tr.customdata.forEach((d,i)=>{const expected=Math.log1p(counts[d[2]]/atlas.totals[d[2]]*run.recipe.target_sum);assert.equal(d[3],expected);assert.equal(tr.marker.color[i],expected)})});
check('changing stage updates existing expression and preserves palette/manual range',()=>{filter('proestrus','analysis',true);const tr=expression.plot.last.traces[0];assert.equal(tr.x.length,217);assert.equal(tr.marker.colorscale,'Magma');assert.equal(tr.marker.cmin,0);assert.equal(tr.marker.cmax,3)});
check('filtered expression CSV includes only displayed cells plus evidence columns',()=>{vm.runInContext('downloadCSV=(name,rows)=>{window.csv=rows}',context);expression.card.children.at(-1).onclick();const csv=context.window.csv;assert.equal(csv.length,218);assert(csv[0].includes('stage_evidence'));for(const row of csv.slice(1)){assert.equal(row.at(-1),'DIRECT_AUTHOR_CELL_MATCH');assert.equal(row.at(-2),'proestrus')}});
check('child run inherits stage/evidence and parent remains visible',()=>{vm.runInContext('openRun("55528ce54cd769a955a2d978",state.get(atlas.qc_run).filter)',context);const child=vm.runInContext('state.get("55528ce54cd769a955a2d978")',context);assert.equal(child.filter.stage,'proestrus');assert(child.filter.direct);const n=child.run.metadata.analysis_stage.reduce((acc,v,i)=>acc+(v==='proestrus'&&child.run.metadata.stage_evidence[i]==='DIRECT_AUTHOR_CELL_MATCH'),0);assert.equal(child.plot.last.traces.reduce((a,t)=>a+t.x.length,0),n);assert.equal(vm.runInContext('state.has(atlas.qc_run)',context),true)});
check('direct evidence is not invented for whole ovary',()=>{vm.runInContext('openRun(atlas.default_run)',context);const full=vm.runInContext('state.get(atlas.default_run)',context);assert.equal(full.evidenceSelect.children[1].disabled,true);assert(full.filterCount.textContent.includes('54,748 evidence not recorded'))});
check('subset without own stage metadata inherits by source row identity',()=>{const child=atlas.runs.find(r=>r.n===3503),parent=atlas.runs.find(r=>r.id===child.parent);context.testId=child.id;const info=vm.runInContext('stageRecords(byId.get(testId))',context);assert.equal(info.length,child.n);child.rows.forEach((row,i)=>assert.equal(info[i].author,parent.metadata.stage[parent.rows.indexOf(row)]))});
check('Japanese filter captions update',()=>{vm.runInContext('language="ja";refreshLanguage()',context);assert(qc.filterCount.textContent.startsWith('表示'));assert.equal(qc.stageSelect.children[1].textContent,'発情前期 (proestrus)')});
check('four-stage comparison keeps all stages while single-stage expression stays filtered',()=>{
 filter('estrus');expression.range=null;expression.compareButton.onclick();
 assert.equal(expression.comparing,true);assert.equal(expression.compareButton['aria-expanded'],'true');assert.equal(expression.comparisonPanels.length,4);
 assert.deepEqual(Array.from(expression.comparisonPanels,p=>p.plot.last.traces[0].x.length),[317,486,193,381]);assert.equal(expression.plot.last.traces[0].x.length,486);
});
check('comparison plots share full-run axes, full-run automatic color scale and exact saved coordinates',()=>{
 const panels=expression.comparisonPanels,low=Math.min(...expression.values),high=Math.max(...expression.values);
 for(const panel of panels){const {traces:[tr],layout}=panel.plot.last;assert.equal(layout.xaxis.fixedrange,true);assert.equal(layout.yaxis.fixedrange,true);assert.equal(tr.marker.cauto,false);assert.equal(tr.marker.cmin,low);assert.equal(tr.marker.cmax,high);assert.deepEqual(layout.xaxis.range,panels[0].plot.last.layout.xaxis.range);assert.deepEqual(layout.yaxis.range,panels[0].plot.last.layout.yaxis.range);tr.customdata.forEach((d,i)=>{const row=run.rows.indexOf(d[2]);assert.equal(tr.x[i],run.coordinates[row][0]);assert.equal(tr.y[i],run.coordinates[row][1]);assert.equal(tr.marker.color[i],expression.values[row])})}
 assert(panels[0].plot.last.layout.xaxis.range[0]<Math.min(...run.coordinates.map(p=>p[0])));assert(panels[0].plot.last.layout.yaxis.range[1]>Math.max(...run.coordinates.map(p=>p[1])));
});
check('comparison follows author label source and direct evidence without changing full-run auto scale',()=>{
 const max=expression.comparisonPanels[0].plot.last.traces[0].marker.cmax;
 filter('unknown','author');assert.equal(expression.plot.last.traces[0].x.length,450);assert.deepEqual(Array.from(expression.comparisonPanels,p=>p.plot.last.traces[0].x.length),[217,266,161,283]);
 filter('estrus','analysis',true);assert.deepEqual(Array.from(expression.comparisonPanels,p=>p.plot.last.traces[0].x.length),[217,266,161,283]);assert.equal(expression.comparisonPanels[0].plot.last.traces[0].marker.cmax,max);
});
check('palette, manual range and automatic reset refresh every comparison panel',()=>{
 expression.colors.value='Plasma';expression.colors.onchange();const form=expression.card.children[2];form.children[1].children[1].value='0';form.children[2].children[1].value='2.75';form.onsubmit({preventDefault(){}});
 for(const panel of expression.comparisonPanels){const marker=panel.plot.last.traces[0].marker;assert.equal(marker.colorscale,'Plasma');assert.equal(marker.cmin,0);assert.equal(marker.cmax,2.75)}
 form.children[4].onclick();assert.equal(expression.range,null);for(const panel of expression.comparisonPanels)assert.equal(panel.plot.last.traces[0].marker.cmax,Math.max(...expression.values));
});
check('comparison captions and toggle refresh in both languages',()=>{
 vm.runInContext('language="en";refreshLanguage()',context);assert(expression.compareButton.textContent.includes('Hide'));assert(expression.comparisonPanels[0].title.textContent.startsWith('Proestrus'));
 vm.runInContext('language="ja";refreshLanguage()',context);assert.equal(expression.compareButton.textContent,'4周期の比較を閉じる');assert(expression.comparisonPanels[0].title.textContent.startsWith('発情前期'));
});
check('zero-cell comparison panels are explicit and retain common scale',()=>{
 const synthetic={...run,id:'empty-stage-fixture',rows:run.rows.slice(0,1),coordinates:run.coordinates.slice(0,1),clusters:run.clusters.slice(0,1),n:1,metadata:{analysis_stage:['proestrus'],stage:['proestrus'],stage_evidence:['DIRECT_AUTHOR_CELL_MATCH']}};
 context.emptyFixture=synthetic;vm.runInContext('addPanel(emptyFixture);addExpression(state.get(emptyFixture.id),testGene,testCounts)',context);const panel=vm.runInContext('state.get(emptyFixture.id)',context),g=panel.expressions.get(gene);g.compareButton.onclick();assert.equal(g.comparisonPanels[0].plot.last.traces[0].x.length,1);for(const empty of g.comparisonPanels.slice(1)){assert.equal(empty.plot.last.traces[0].x.length,0);assert(empty.empty.textContent.includes('該当'));assert.equal(empty.plot.last.traces[0].marker.cmax,g.comparisonPanels[0].plot.last.traces[0].marker.cmax)}
 panel.article.children[0].children[1].children[2].onclick();assert(g.comparisonPanels.length===0);
});
check('hide and expression close purge comparison plots; reopening does not duplicate panels',()=>{
 const old=expression.comparisonPanels.map(p=>p.plot);expression.compareButton.onclick();assert.equal(expression.comparisonHost.hidden,true);assert.equal(expression.comparisonPanels.length,0);old.forEach(p=>assert(purged.includes(p)));
 expression.compareButton.onclick();assert.equal(expression.comparisonPanels.length,4);const reopened=expression.comparisonPanels.map(p=>p.plot);expression.card.children[0].children[1].onclick();reopened.forEach(p=>assert(purged.includes(p)));assert(purged.includes(expression.plot));assert(!qc.expressions.has(gene));
});
check('reset purges active comparison plots and restores initial view',()=>{
 vm.runInContext('addExpression(state.get(atlas.qc_run),testGene,testCounts)',context);const g=qc.expressions.get(gene);g.compareButton.onclick();const panels=g.comparisonPanels.map(p=>p.plot);els['#reset'].onclick();panels.forEach(p=>assert(purged.includes(p)));assert.equal(g.comparisonPanels.length,0);assert.equal(vm.runInContext('state.get(atlas.qc_run).expressions.size',context),0);
});
const result={test_type:'non-browser VM data/interaction contracts',passed:checks.length,checks};fs.writeFileSync(path.join(root,'docs/STAGE_FILTER_CHECKS.json'),JSON.stringify(result,null,2));console.log(JSON.stringify(result));
