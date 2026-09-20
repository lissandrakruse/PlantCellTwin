const $ = (s) => document.querySelector(s);
const $$ = (s) => document.querySelectorAll(s);
let DATA = null, running = true, start = Date.now(), phase = 0, condition = 'flight', surveys = [];
const tissueNames = { root: 'Root', hypocotyl: 'Hypocotyl', shoot: 'Shoot / leaves', 'whole plant': 'Whole plant' };
const tissueNotes = {
  root: 'Roots showed a tissue-specific response that is diluted by whole-plant averaging. Candidate effects include signaling, oxygen-response and cell-wall-associated genes.',
  hypocotyl: 'The hypocotyl contrast produces a distinct candidate ranking, supporting an organ-aware twin rather than one universal plant-cell state.',
  shoot: 'Shoot tissue has its own flight–ground effect profile; the twin treats photosynthetic tissue separately from roots and hypocotyls.',
  'whole plant': 'Whole-plant averaging has fewer replicates and can conceal organ-specific responses, matching the main biological lesson of the NASA study.'
};
const tissue = $('#tissue'), exposure = $('#exposure'), light = $('#light');
const currentResult = () => DATA.tissues[tissue.value];
function projection() {
  const e = condition === 'flight' ? +exposure.value : 0;
  return Math.min(1, currentResult().median_abs_log2fc * e + Math.abs(+light.value) / 100 * .22);
}
function update() {
  if (!DATA) return;
  const d = currentResult(), p = projection(), idx = Math.round(p * 100);
  $('#exposureOut').value = Math.round(+exposure.value * 100) + '%';
  $('#lightOut').value = (+light.value > 0 ? '+' : '') + light.value + '%';
  $('#replicateLabel').textContent = `${d.n_flight} flight + ${d.n_ground} ground samples`;
  $('#sampleCount').textContent = d.n_flight + d.n_ground;
  $('#fdrHits').textContent = d.significant_total;
  $('#candidateCount').textContent = d.exploratory_candidates;
  $('#candidateLabel').textContent = d.exploratory_candidates;
  $('#effectLabel').textContent = (d.median_abs_log2fc * (condition === 'flight' ? +exposure.value : 0)).toFixed(3);
  $('#responseIndex').textContent = idx;
  $('#responseBar').style.width = idx + '%';
  $('#cellState').textContent = condition === 'ground' ? 'GROUND REFERENCE' : (p > .35 ? 'AMPLIFIED FLIGHT HYPOTHESIS' : 'FLIGHT RESPONSE');
  const extra = +light.value === 0 ? '' : ` The ${+light.value > 0 ? 'increased' : 'reduced'}-light term is an explicit model extrapolation, not a measured GLDS-7 variable.`;
  $('#explanation').textContent = tissueNotes[tissue.value] + ` Independent Welch tests found ${d.significant_total} probes meeting FDR < 0.05 and |log₂FC| ≥ 1; ${d.exploratory_candidates} meet the unadjusted candidate screen.` + extra;
  drawGeneChart();
}
tissue.addEventListener('change', update);
[exposure, light].forEach(i => i.addEventListener('input', update));
$$('[data-condition]').forEach(b => b.onclick = () => {
  condition = b.dataset.condition;
  $$('[data-condition]').forEach(x => x.classList.toggle('active', x === b));
  exposure.disabled = condition === 'ground'; update();
});
$('#resetBtn').onclick = () => {
  tissue.value = 'root'; condition = 'flight'; exposure.value = 1; light.value = 0; exposure.disabled = false;
  $$('[data-condition]').forEach(x => x.classList.toggle('active', x.dataset.condition === 'flight')); update();
};
function renderSurveys() {
  $('#surveyStatus').textContent = `${surveys.length} registro${surveys.length === 1 ? '' : 's'}`;
  $('#surveyRows').innerHTML = surveys.length ? surveys.map((r, i) => `<tr><td>${i + 1}</td><td>${r.tissue}</td><td>${r.condition === 'flight' ? 'ISS flight' : '1 g ground'}</td><td>${r.exposure}%</td><td>${r.light > 0 ? '+' : ''}${r.light}%</td><td>${r.index}</td><td>${r.gene} (${r.log2fc > 0 ? '+' : ''}${r.log2fc})</td></tr>`).join('') : '<tr class="empty-row"><td colspan="7">Registre um cenário para iniciar o levantamento comparativo.</td></tr>';
}
$('#recordBtn').onclick = () => {
  const d = currentResult(), gene = d.top_genes[0];
  surveys.push({ tissue: tissueNames[tissue.value], condition, exposure: Math.round(+exposure.value * 100), light: +light.value, index: Math.round(projection() * 100), gene: gene.symbol, log2fc: gene.log2fc });
  renderSurveys();
};
$('#clearBtn').onclick = () => { surveys = []; renderSurveys(); };
$('#exportBtn').onclick = () => {
  if (!surveys.length) return;
  const head = ['id','tissue','condition','exposure_percent','light_perturbation_percent','twin_index','top_probe_gene','observed_log2fc'];
  const rows = surveys.map((r,i) => [i+1,r.tissue,r.condition,r.exposure,r.light,r.index,r.gene,r.log2fc]);
  const csv = [head, ...rows].map(row => row.map(v => `"${String(v).replaceAll('"','""')}"`).join(',')).join('\n');
  const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([csv], {type:'text/csv'})); a.download = 'glds7_digital_twin_survey.csv'; a.click(); URL.revokeObjectURL(a.href);
};
setInterval(() => {
  if (!running) return;
  const s = Math.floor((Date.now() - start) / 1000);
  $('#clock').textContent = `T+ ${String(Math.floor(s/3600)).padStart(2,'0')}:${String(Math.floor(s%3600/60)).padStart(2,'0')}:${String(s%60).padStart(2,'0')}`;
}, 1000);
const cell = $('#cellCanvas'), cc = cell.getContext('2d'), chart = $('#chartCanvas'), gc = chart.getContext('2d');
const particles = Array.from({length:45}, () => ({a:Math.random()*6.28,r:.35+Math.random()*.6,z:Math.random()*2-1,s:.002+Math.random()*.005}));
function sizeCanvas(c) { const d=devicePixelRatio||1,r=c.getBoundingClientRect(); if(c.width!==r.width*d||c.height!==r.height*d){c.width=r.width*d;c.height=r.height*d;c.getContext('2d').setTransform(d,0,0,d,0,0)} return r; }
function drawCell(){
  const r=sizeCanvas(cell),w=r.width,h=r.height,x=w/2,y=h/2+15,R=Math.min(w,h)*.27,rr=DATA?projection():0;
  cc.clearRect(0,0,w,h); phase+=running?.008:0; cc.save();cc.translate(x,y);cc.rotate(Math.sin(phase*.25)*.025);
  cc.beginPath();cc.roundRect(-R*1.25,-R*.9,R*2.5,R*1.8,R*.22);cc.fillStyle='rgba(25,82,51,.2)';cc.fill();cc.lineWidth=5;cc.strokeStyle=rr>.35?'#ff9a51':'#8dda68';cc.shadowColor=cc.strokeStyle;cc.shadowBlur=18;cc.stroke();cc.shadowBlur=0;cc.lineWidth=1.5;cc.strokeStyle='rgba(186,252,74,.55)';cc.strokeRect(-R*1.17,-R*.82,R*2.34,R*1.64);
  const vac=cc.createLinearGradient(-R,0,R,0);vac.addColorStop(0,'rgba(43,183,162,.12)');vac.addColorStop(1,'rgba(118,224,195,.28)');cc.beginPath();cc.ellipse(R*.18,0,R*.67,R*.57,0,0,Math.PI*2);cc.fillStyle=vac;cc.fill();cc.strokeStyle='rgba(96,231,204,.55)';cc.stroke();
  cc.beginPath();cc.arc(-R*.48,R*.08,R*.22,0,Math.PI*2);cc.fillStyle=rr>.35?'rgba(255,154,81,.58)':'rgba(186,252,74,.48)';cc.fill();cc.strokeStyle='#d4ff75';cc.stroke();
  for(let i=0;i<8;i++){const a=i/8*Math.PI*2+phase*.2,px=Math.cos(a)*R*.82,py=Math.sin(a)*R*.55;cc.save();cc.translate(px,py);cc.rotate(a);cc.beginPath();cc.ellipse(0,0,R*.16,R*.075,0,0,Math.PI*2);cc.fillStyle=rr>.35?'rgba(255,154,81,.7)':'rgba(113,221,72,.78)';cc.fill();cc.strokeStyle=rr>.35?'#ff9a51':'#bafc4a';cc.stroke();cc.beginPath();cc.moveTo(-R*.1,0);cc.lineTo(R*.1,0);cc.strokeStyle='rgba(5,40,20,.75)';cc.stroke();cc.restore();}
  particles.forEach(p=>{if(running)p.a+=p.s*(1+rr*4);const px=Math.cos(p.a)*R*p.r,py=Math.sin(p.a)*R*p.r*.7+p.z*8;cc.beginPath();cc.arc(px,py,1.2+(p.z+1),0,7);cc.fillStyle=rr>.35&&p.r>.7?'#ff9a51':'rgba(190,255,170,.72)';cc.fill();});cc.restore();requestAnimationFrame(drawCell);
}
function drawGeneChart(){
  if(!DATA)return; const r=sizeCanvas(chart),w=r.width,h=r.height,d=currentResult(),genes=d.top_genes.slice(0,8);gc.clearRect(0,0,w,h);
  const mid=w/2,max=Math.max(1,...genes.map(g=>Math.abs(g.log2fc))),row=h/(genes.length+1);gc.strokeStyle='rgba(135,164,156,.35)';gc.beginPath();gc.moveTo(mid,8);gc.lineTo(mid,h-5);gc.stroke();
  genes.forEach((g,i)=>{const y=(i+1)*row,obs=g.log2fc,proj=obs*(condition==='flight'?+exposure.value:0)+(Math.sign(obs)||1)*(+light.value/100*.25);gc.font='10px monospace';gc.textBaseline='middle';gc.fillStyle='#a7beb8';gc.textAlign=obs<0?'left':'right';gc.fillText(g.symbol,obs<0?mid+5:mid-5,y-5);gc.fillStyle='#3ce5c2';gc.fillRect(mid,y,obs/max*(w*.34),4);gc.fillStyle='#bafc4a';gc.fillRect(mid,y+6,proj/max*(w*.34),3);});
}
function renderValidation(v) {
  $('#replicatedCount').textContent = v.summary.replicated;
  $('#partialCount').textContent = v.summary.partial;
  const label = {replicated:'REPLICADO', partial:'PARCIAL', 'not replicated':'NÃO REPLICADO'};
  const shown = v.comparisons.slice(0, 10);
  $('#validationRows').innerHTML = shown.map(r => `<tr><td>${tissueNames[r.tissue]}</td><td>${r.study}</td><td>${r.overlap}</td><td>${Math.round(r.direction_agreement*100)}%</td><td>${r.overlap_fdr < .001 ? '&lt;0.001' : r.overlap_fdr.toFixed(3)}</td><td>${r.direction_fdr < .001 ? '&lt;0.001' : r.direction_fdr.toFixed(3)}</td><td><span class="status ${r.status.replace(' ','-')}">${label[r.status]}</span></td></tr>`).join('');
}
function renderAdvanced(v) {
  $('#goCount').textContent = v.go_summary.fdr_005;
  $('#predictiveCount').textContent = v.transfer_summary.predictive;
  $('#rankCount').textContent = v.transfer_summary.rank_consistent;
  const go = v.go_summary.top.slice(0, 6);
  $('#goRows').innerHTML = go.map(x => `<article><b>${tissueNames[x.tissue]} · ${x.direction === 'up' ? 'aumentado' : 'reduzido'}</b><span>${x.term}</span><small>FDR ${x.fdr < .001 ? '&lt;0,001' : x.fdr.toFixed(3).replace('.',',')} · ${x.overlap} genes</small></article>`).join('');
  const tr = v.transfer_summary.top.filter(x => x.transfer_status !== 'not supported').slice(0, 6);
  const labels = {predictive:'PREDITIVO', 'rank-consistent':'RANK CONSISTENTE'};
  $('#transferRows').innerHTML = tr.map(x => `<article><b>${tissueNames[x.tissue]} → ${x.study}</b><span>${labels[x.transfer_status]}</span><small>ρ ${x.spearman_rho.toFixed(3)} · AUC ${x.direction_auc == null ? '—' : x.direction_auc.toFixed(3)}</small></article>`).join('');
}
function renderDefense(v) {
  $('#defensePairCount').textContent = v.total_replicated_pairs;
  const supported = v.summary.filter(x => x.replicated_pathways > 0);
  $('#defenseStudies').innerHTML = supported.map(s => `<article><header><b>${s.study}</b><strong>${s.replicated_pathways} via${s.replicated_pathways === 1 ? '' : 's'}</strong></header>${s.pathways.slice(0,7).map(p => `<div><span>${p.term}</span><small>${p.overlap} genes · FDR ${p.fdr < .001 ? '&lt;0,001' : p.fdr.toFixed(3).replace('.',',')}</small></div>`).join('')}</article>`).join('');
}
function renderRaw(raw, plates) {
  const root = raw.tissue_summary.find(x => x.tissue === 'root');
  const sensitivity = raw.pipeline_sensitivity.find(x => x.tissue === 'root');
  $('#rawRootHits').textContent = root.fdr_005;
  $('#rawPathways').textContent = raw.defense_confirmation.confirmed_terms;
  $('#stablePathways').textContent = plates.stable_terms;
  $('#rootPipelineRho').textContent = sensitivity.effect_spearman.toFixed(3).replace('.', ',');
  $('#stableTermRows').innerHTML = plates.terms.map(x => `<span>${x.term}<small>${x.positive_plates}/3 placas</small></span>`).join('');
}
function renderExpansion(v){
  const role={1:'VALIDAÇÃO DIRETA',2:'CONTEXTO',3:'DOSE DE GRAVIDADE',4:'REGULAÇÃO',5:'MECANISMO',6:'GRAVISSENSIBILIDADE',7:'CÉLULA',8:'CONTROLE RADIAÇÃO',9:'CONTROLE HIPOBARIA'};
  const completed={'GLDS-208':'ANÁLISE CONCLUÍDA','GLDS-251':'ANÁLISE CONCLUÍDA'};
  $('#expansionRoadmap').innerHTML=v.expansion_studies.sort((a,b)=>a.priority-b.priority).map(x=>`<article><header><b>${x.study}</b><span>${completed[x.study]||role[x.priority]}</span></header><strong>${x.material}</strong><p>${x.rationale}</p></article>`).join('');
}
Promise.all([fetch('glds7-results.json').then(r=>r.json()), fetch('nasa-cross-validation.json').then(r=>r.json()), fetch('advanced-validation.json').then(r=>r.json()), fetch('defense-validation.json').then(r=>r.json()), fetch('raw-rma-confirmation.json').then(r=>r.json()), fetch('plate-stability.json').then(r=>r.json()), fetch('bioconductor-summary.json').then(r=>r.json()), fetch('expanded-evidence.json').then(r=>r.json())]).then(([data, validation, advanced, defense, raw, plates, bioc, expansion])=>{DATA=data;update();renderValidation(validation);renderAdvanced(advanced);renderDefense(defense);renderRaw(raw,plates);renderExpansion(expansion);const root=bioc.concordance.find(x=>x.tissue==='root');$('#biocRootHits').textContent=root.bioc_fdr_005;$('#biocRootRho').textContent=Number(root.effect_spearman).toFixed(3);$('#biocTopOverlap').textContent=root.top300_overlap+'/300';$('#biocDirection').textContent=Math.round(root.top300_same_direction*100)+'%';drawCell();}).catch(()=>{$('#explanation').textContent='The NASA analysis files could not be loaded. Reload the page to retry.';});
