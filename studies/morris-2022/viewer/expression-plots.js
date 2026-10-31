/* Browser-only expression summaries. No sampling: all in-scope cells, including zero counts. */
(function (root) {
 'use strict';
 const MAX_GENES = 12, BINS = 80;
 function summarize(values) {
  const sorted=Array.from(values).sort((a,b)=>a-b),n=sorted.length;
  if(sorted.some(v=>!Number.isFinite(v)||v<0))throw Error('Invalid normalized expression');
  if(!n)return {n:0,mean:null,median:null,min:null,max:null,positive:0};
  const mid=Math.floor(n/2);
  return {n,mean:sorted.reduce((a,b)=>a+b,0)/n,median:n%2?sorted[mid]:(sorted[mid-1]+sorted[mid])/2,min:sorted[0],max:sorted[n-1],positive:sorted.filter(v=>v>0).length};
 }
 function normalize(counts,rows,totals,target) {
  if(!(target>0&&Number.isFinite(target)))throw Error('Invalid normalization target');
  return Float64Array.from(rows,row=>{
   if(!(totals[row]>0&&Number.isFinite(totals[row]))||!(counts[row]>=0&&Number.isFinite(counts[row])))throw Error('Invalid expression counts or total');
   return Math.log1p(counts[row]/totals[row]*target);
  });
 }
 function rowZscores(matrix) {
  return matrix.map(row=>{
   const valid=row.filter(v=>v!==null&&Number.isFinite(v));
   if(!valid.length)return row.map(()=>null);
   const mean=valid.reduce((a,b)=>a+b,0)/valid.length;
   const sd=Math.sqrt(valid.reduce((a,b)=>a+(b-mean)**2,0)/valid.length);
   return row.map(v=>v===null?null:sd>1e-12?(v-mean)/sd:0);
  });
 }
 // Fixed-grid, Gaussian-smoothed histogram density: O(n + BINS²), independent of cell count.
 // The curve is clipped to the nonnegative support; constant groups are not smoothed.
 function density(values,maximum,bins=BINS) {
  if(!values.length)return null;
  let low=Infinity,high=-Infinity;for(const v of values){low=Math.min(low,v);high=Math.max(high,v)}
  if(low===high)return {constant:low,x:[low],width:[0]};
  const end=maximum>0?maximum:high,step=end/(bins-1),hist=new Float64Array(bins);
  for(const value of values)hist[Math.min(bins-1,Math.max(0,Math.round(value/step)))]+=1;
  // Bounded 2-bin bandwidth; all cells contribute, without passing tens of thousands of points to Plotly.
  const smooth=Array.from({length:bins},(_,i)=>{
   let sum=0;for(let j=Math.max(0,i-8);j<=Math.min(bins-1,i+8);j++)sum+=hist[j]*Math.exp(-.5*((i-j)/2)**2);return sum;
  });
  const peak=Math.max(...smooth);
  return {constant:null,x:Array.from({length:bins},(_,i)=>i*step),width:smooth.map(v=>v/peak*.38)};
 }
 const words={en:{title:'Violin plots / heatmap',genes:'Genes (up to 12)',draw:'Draw plots',kind:'Plot',violin:'Violin plot',heatmap:'Heatmap',group:'Group by',cluster:'Cluster',stage:'Cycle stage',scope:'Cell group',all:'All clusters',selected:'Only lasso / box-selected cells',scale:'Heatmap scale',mean:'Mean log-normalized expression',zscore:'Row Z-score of group means',gene:'Violin gene',empty:'No cells match these filters.',noGenes:'Enter one or more gene symbols or feature IDs.',limit:'Choose at most 12 genes.',loading:'Loading gene data…',ready:'cells · zeros included',noData:'No cells',csv:'Summary CSV',note:'Uses the stage and evidence filters above. Violin width shows relative density within each group, not cell number; dots mark medians. Constant groups are shown as lines.',zNote:'Row Z-scores compare group means within each gene (population SD); constant rows are 0. Grey cells have no data.',meanNote:'Arithmetic mean of per-cell ln(1 + count / full-gene total × target sum), including zeros. Grey cells have no data.',value:'Log-normalized expression',n:'Cells',positive:'Expressing (%)',median:'Median',pending:'Loading new genes; current plots hidden.',error:'Could not draw plots. Retry or choose fewer genes.',selectionEmpty:'No lasso / box selection. Select cells above or turn off the selection filter.',allLabel:'all',file:'expression-summary'},ja:{title:'バイオリンプロット / ヒートマップ',genes:'遺伝子（12個まで）',draw:'図を表示',kind:'表示',violin:'バイオリンプロット',heatmap:'ヒートマップ',group:'比較する群',cluster:'クラスター',stage:'性周期',scope:'対象クラスター',all:'全クラスター',selected:'投げ縄・矩形で選択した細胞のみ',scale:'ヒートマップの値',mean:'平均発現（log正規化）',zscore:'遺伝子ごとの群平均 Zスコア',gene:'バイオリンの遺伝子',empty:'この条件に該当する細胞はありません。',noGenes:'遺伝子名またはfeature IDを入力してください。',limit:'遺伝子は12個以内で指定してください。',loading:'遺伝子データを読み込み中…',ready:'細胞 · 発現ゼロも集計',noData:'細胞なし',csv:'集計CSV',note:'上の性周期・対応根拠の絞り込みを反映します。バイオリンの幅は各群内の相対密度で、細胞数ではありません。点は中央値、一定値の群は線で示します。',zNote:'遺伝子ごとに群平均をZスコア化（母標準偏差）。一定の行は0、細胞がないマスは灰色です。',meanNote:'各細胞の ln(1 + count / 全遺伝子合計 × target sum) の算術平均。発現ゼロを含み、細胞がないマスは灰色です。',value:'発現（log正規化）',n:'細胞数',positive:'発現細胞 (%)',median:'中央値',pending:'新しい遺伝子を読み込み中です。',error:'図を表示できませんでした。再試行するか遺伝子数を減らしてください。',selectionEmpty:'投げ縄・矩形の選択がありません。上の図で細胞を選ぶか、選択細胞のみの条件を外してください。',allLabel:'全て',file:'expression-summary'}};
 function mount(s,a) {
  const p={s,a,request:0,data:null,disposed:false,models:null};s.expressionPlots=p;
  const tx=k=>words[a.language()][k]||k;p.tx=tx;
  const node=(tag,props={},children=[])=>a.el(tag,props,children);
  const local=(tag,key,props={})=>{const e=node(tag,{...props,text:tx(key)});e.plotKey=key;p.labels.push(e);return e};
  p.labels=[];
  const select=(key,options)=>node('select',{'aria-label':tx(key)},options.map(([value,k])=>local('option',k,{value})));
  const wrap=(key,input)=>node('label',{},[local('span',key),input]);
  p.input=node('input',{type:'text',placeholder:'Mafb, Alox15, C1qa','aria-label':tx('genes'),autocomplete:'off'});
  p.kind=select('kind',[['violin','violin'],['heatmap','heatmap']]);
  p.group=select('group',[['cluster','cluster'],['stage','stage']]);
  p.scope=node('select',{'aria-label':tx('scope')});
  p.selected=node('input',{type:'checkbox'});
  p.scale=select('scale',[['mean','mean'],['zscore','zscore']]);
  p.gene=node('select',{'aria-label':tx('gene')});
  p.scaleWrap=wrap('scale',p.scale);p.geneWrap=wrap('gene',p.gene);
  p.submit=local('button','draw',{type:'submit'});
  p.status=node('p',{class:'expression-status',role:'status','aria-live':'polite'});
  p.error=node('p',{class:'error',role:'alert'});
  p.plot=node('div',{class:'plot summary-plot','data-kind':'expression-summary'});p.plot.hidden=true;
  p.note=local('p','note',{class:'hint summary-note'});p.scaleNote=node('p',{class:'hint summary-note'});
  p.export=local('button','csv',{class:'secondary',type:'button',onclick:()=>exportSummary(p)});p.export.disabled=true;
  p.form=node('form',{class:'controls summary-form'},[wrap('genes',p.input),p.submit]);
  const options=node('div',{class:'controls summary-options'},[wrap('kind',p.kind),wrap('group',p.group),wrap('scope',p.scope),p.geneWrap,p.scaleWrap]);
  p.host=node('section',{class:'expression-summaries'},[local('h3','title'),p.form,options,node('label',{class:'summary-selection'},[p.selected,local('span','selected')]),p.status,p.error,p.plot,p.note,p.scaleNote,p.export]);
  p.form.onsubmit=e=>{e.preventDefault();load(p)};
  for(const control of[p.kind,p.group,p.scope,p.selected,p.scale,p.gene])control.onchange=()=>render(p);
  p.gene.value='';p.kind.value='violin';p.group.value='cluster';p.scale.value='mean';
  p.onResize=()=>{if(!p.disposed)render(p)};root.addEventListener?.('resize',p.onResize);
  refresh(s);return p.host;
 }
 async function load(p) {
  const {a,s}=p,tx=p.tx,request=++p.request;
  p.error.textContent='';p.status.textContent='';p.loading=false;p.data=null;p.models=null;p.plot.hidden=true;p.export.disabled=true;
  let genes;
  try {
   const queries=[...new Set(p.input.value.split(/[\s,;、]+/u).filter(Boolean))];
   if(!queries.length)throw Error(tx('noGenes'));
   if(queries.length>MAX_GENES)throw Error(tx('limit'));
   genes=[...new Set(queries.map(a.resolve))];
  }catch(e){p.error.textContent=e.message;return}
  p.loading=true;p.plot.hidden=true;p.export.disabled=true;p.status.textContent=tx('loading');
  try {
   const data=[];
   // Sequential bounded loading shares the existing chunk and decoded-gene caches.
   for(const gene of genes){
    const counts=await a.geneCounts(gene);
    if(p.disposed||request!==p.request||!a.alive())return;
    data.push({gene,feature:a.atlas.features[gene],values:normalize(counts,s.run.rows,a.atlas.totals,s.run.recipe.target_sum)});
   }
   if(p.disposed||request!==p.request||!a.alive())return;
   p.data=data;p.loading=false;
   const previous=p.gene.value;
   p.gene.replaceChildren(...data.map(d=>a.el('option',{value:String(d.gene),text:d.feature.symbol+' · '+d.feature.feature_id})));
   p.gene.value=data.some(d=>String(d.gene)===previous)?previous:String(data[0].gene);
   render(p);
  }catch(e){if(!p.disposed&&request===p.request&&a.alive()){p.loading=false;p.status.textContent='';p.error.textContent=e.message||tx('error');p.plot.hidden=true;p.export.disabled=true}}
 }
 function makeGroups(p) {
  const {s,a}=p,byStage=p.group.value==='stage',selection=new Set(s.selected),scope=p.scope.value;
  let keys=byStage?['proestrus','estrus','metestrus','diestrus','unknown']:[...new Set(s.run.clusters.map(String))];
  if(!byStage&&scope!=='')keys=keys.filter(k=>k===scope);
  const groups=keys.map(key=>({key,label:byStage?a.stageName(key):a.clusterCaption(s.run,key),indices:[]})),lookup=new Map(groups.map(g=>[g.key,g]));
  for(const i of a.visibleIndices(s.run,s.filter)){
   if(scope!==''&&String(s.run.clusters[i])!==scope)continue;
   if(p.selected.checked&&!selection.has(s.run.rows[i]))continue;
   const key=byStage?a.stageInfo(s.run,i)[s.filter.source]:String(s.run.clusters[i]);
   if(!lookup.has(key)){const g={key,label:a.stageName(key),indices:[]};groups.push(g);lookup.set(key,g)}
   lookup.get(key).indices.push(i);
  }
  return groups;
 }
 function render(p) {
  const {a,s,tx}=p;
  p.scaleWrap.hidden=p.kind.value!=='heatmap';p.geneWrap.hidden=p.kind.value!=='violin';
  p.scaleNote.textContent=p.kind.value==='heatmap'?tx(p.scale.value==='zscore'?'zNote':'meanNote'):'';
  if(p.disposed||p.loading||!p.data)return;
  const groups=makeGroups(p),models=p.data.map(d=>({...d,groups:groups.map(g=>{const values=g.indices.map(i=>d.values[i]);return {...g,values,summary:summarize(values)}})}));
  p.models=models;p.groups=groups;
  const total=groups.reduce((n,g)=>n+g.indices.length,0);
  p.status.textContent=total?`${total.toLocaleString()} ${tx('ready')}`:p.selected.checked&&!s.selected.length?tx('selectionEmpty'):tx('empty');
  p.export.disabled=!total;p.plot.hidden=!total;
  if(!total)return;
  p.error.textContent='';
  const compact=(p.host.clientWidth||root.innerWidth||1000)<600;
  const clusterTick=g=>/^\d+$/.test(g.key)?'C'+g.key:g.key;
  const axisLabel=g=>{const label=g.label.replace(/^(?:クラスター|Cluster)\s+/u,'');return label.length>22?label.slice(0,20)+'…':label};
  const base={autosize:true,margin:{l:compact?110:280,r:compact?20:35,t:20,b:70},paper_bgcolor:'#fff',plot_bgcolor:'#eef1ee',showlegend:false,hovermode:'closest',uirevision:s.run.id+':summary:'+p.group.value+':'+p.kind.value};
  let traces;
  if(p.kind.value==='heatmap') {
   const means=models.map(m=>m.groups.map(g=>g.summary.mean)),isZ=p.scale.value==='zscore',z=isZ?rowZscores(means):means;
   const labels=groups.map(g=>g.label),custom=models.map(m=>m.groups.map(g=>[g.label,g.summary.n,g.summary.mean,g.summary.positive&&g.summary.n?g.summary.positive/g.summary.n*100:0,m.feature.symbol]));
   traces=[{type:'heatmap',x:groups.map((_,i)=>i),y:models.map((m,i)=>i),z,customdata:custom,colorscale:isZ?'RdBu':'Viridis',reversescale:isZ,zmid:isZ?0:undefined,zmin:isZ?-Math.max(1,...z.flat().filter(Number.isFinite).map(Math.abs)):0,zmax:isZ?Math.max(1,...z.flat().filter(Number.isFinite).map(Math.abs)):Math.max(1e-9,...means.flat().filter(Number.isFinite)),hoverongaps:false,colorbar:{title:{text:isZ?'Z':tx('value')},thickness:14},hovertemplate:'%{customdata[0]}<br>%{customdata[4]}<br>'+tx('n')+': %{customdata[1]}<br>'+tx('mean')+': %{customdata[2]:.4f}<br>'+tx('positive')+': %{customdata[3]:.2f}<br>'+(isZ?'Z: %{z:.3f}':'')+'<extra></extra>'}];
   base.xaxis={tickmode:'array',tickvals:groups.map((_,i)=>i),ticktext:groups.map(g=>p.group.value==='cluster'?clusterTick(g):g.label),tickangle:-35,range:[-.5,groups.length-.5]};
   base.yaxis={tickmode:'array',tickvals:models.map((_,i)=>i),ticktext:models.map(m=>m.feature.symbol),autorange:'reversed'};
   base.margin.l=compact?65:110;base.margin.r=compact?65:100;base.xaxis.tickfont={size:compact?10:12};base.height=Math.max(330,models.length*38+160);
   base.annotations=[];models.forEach((m,i)=>m.groups.forEach((g,j)=>{if(!g.summary.n)base.annotations.push({x:j,y:i,text:'—',showarrow:false,font:{color:'#687870',size:12}})}));
  } else {
   const model=models.find(m=>String(m.gene)===p.gene.value)||models[0];let maximum=0;
   for(const g of model.groups)if(g.summary.max!==null)maximum=Math.max(maximum,g.summary.max);
   maximum=maximum||1;traces=[];
   model.groups.forEach((g,i)=>{
    if(!g.summary.n)return;
    const d=density(g.values,maximum),color=['#257c74','#ce724b','#7875aa','#bf9a3e','#6699ae','#a35c7c','#7e9661'][i%7];
    if(d.constant===null){traces.push({type:'scatter',mode:'lines',x:d.x.concat(d.x.slice().reverse()),y:d.width.map(w=>i+w).concat(d.width.slice().reverse().map(w=>i-w)),fill:'toself',fillcolor:color,line:{color,width:1},opacity:.65,hoverinfo:'skip'})}
    else traces.push({type:'scatter',mode:'lines',x:[d.constant,d.constant],y:[i-.22,i+.22],line:{color,width:3},hoverinfo:'skip'});
    traces.push({type:'scatter',mode:'markers',x:[g.summary.median],y:[i],marker:{color:'#163c37',size:8},customdata:[[g.label,g.summary.n,g.summary.mean,g.summary.positive/g.summary.n*100]],hovertemplate:g.label+'<br>'+model.feature.symbol+'<br>'+tx('n')+': %{customdata[1]}<br>'+tx('median')+': %{x:.4f}<br>'+tx('mean')+': %{customdata[2]:.4f}<br>'+tx('positive')+': %{customdata[3]:.2f}<extra></extra>'});
   });
   base.xaxis={title:{text:model.feature.symbol+' · '+(compact?(a.language()==='ja'?'log正規化':'log-normalized'):tx('value'))},range:[-maximum*.025,maximum*1.05],zeroline:true};
   base.yaxis={tickmode:'array',tickvals:groups.map((_,i)=>i),ticktext:groups.map(g=>(compact?(p.group.value==='cluster'?clusterTick(g):g.label.split(' (')[0]):axisLabel(g))+' (n='+g.indices.length+')'),tickfont:{size:compact?10:12},range:[groups.length-.4,-.6]};
   base.height=Math.max(360,groups.length*38+100);base.annotations=groups.flatMap((g,i)=>g.indices.length?[]:[{x:maximum*.5,y:i,text:tx('noData'),showarrow:false,font:{color:'#687870'}}]);
  }
  p.plot.style.height=base.height+'px';
  p.rendered={traces,layout:base};
  try {const result=Plotly.react(p.plot,traces,base,{...a.config,toImageButtonOptions:{format:'png',filename:'expression-'+p.kind.value,scale:2}});if(result?.catch)result.catch(()=>{if(!p.disposed)p.error.textContent=tx('error')})}catch(e){p.error.textContent=tx('error')}
 }
 function exportSummary(p) {
  if(!p.models||p.export.disabled)return;
  const {s,a}=p,z=rowZscores(p.models.map(m=>m.groups.map(g=>g.summary.mean)));
  const rows=[['run_id','gene_symbol','feature_id','group_by','group','cell_type_or_stage','n_cells','n_expressing','mean_lognorm','median_lognorm','group_mean_row_zscore','stage_filter','stage_label_source','direct_author_only','cluster_scope','selection_only','normalization_target_sum']];
  p.models.forEach((m,i)=>m.groups.forEach((g,j)=>rows.push([s.run.id,m.feature.symbol,m.feature.feature_id,p.group.value,g.key,g.label,g.summary.n,g.summary.positive,g.summary.mean??'',g.summary.median??'',z[i][j]??'',s.filter.stage||'all',s.filter.source,s.filter.direct,p.scope.value||'all',!!p.selected.checked,s.run.recipe.target_sum])));
  a.downloadCSV('expression-summary-'+s.run.id+'.csv',rows);
 }
 function refresh(s) {
  const p=s.expressionPlots;if(!p||p.disposed)return;
  for(const label of p.labels)label.textContent=p.tx(label.plotKey);
  for(const [control,key]of[[p.input,'genes'],[p.kind,'kind'],[p.group,'group'],[p.scope,'scope'],[p.scale,'scale'],[p.gene,'gene']])control.setAttribute('aria-label',p.tx(key));
  const scope=p.scope.value;
  p.scope.replaceChildren(p.a.el('option',{value:'',text:p.tx('all')}),...[...new Set(s.run.clusters)].map(c=>p.a.el('option',{value:String(c),text:p.a.clusterCaption(s.run,c)})));
  p.scope.value=scope||'';
  if(p.loading)p.status.textContent=p.tx('loading');
  render(p);
 }
 function dispose(s) {const p=s.expressionPlots;if(p){p.disposed=true;++p.request;root.removeEventListener?.('resize',p.onResize);Plotly.purge(p.plot);p.data=null;p.models=null}}
 const api={mount,refresh,dispose,load,render,makeGroups,math:{summarize,normalize,rowZscores,density},MAX_GENES};
 if(typeof module!=='undefined'&&module.exports)module.exports=api;
 root.ExpressionPlots=api;
})(typeof window!=='undefined'?window:globalThis);
