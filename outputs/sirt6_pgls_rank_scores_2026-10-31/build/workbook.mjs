import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const outArg=process.argv.find(x=>x.startsWith('--out='));
const out=outArg?path.resolve(outArg.slice(6)):path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const d=JSON.parse(await fs.readFile(path.join(out,'analysis.json'),'utf8'));
const shared=d.lambda_mode==='shared_by_assay';
const lambdaLabel=shared?'Shared lambda':'ML lambda';
const wb=Workbook.create();
const names=['Summary','Gene tests','Rank scores','Species groups'];
const sheets=Object.fromEntries(names.map(n=>[n,wb.worksheets.add(n)]));
const col=n=>{let s='';for(;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
const value=x=>x??'n.a.';
const assay=a=>a==='rna'?'Liver RNA':'ATAC TSS ±5 kb';
function base(sh,end,widths,title){
  sh.showGridLines=false;
  sh.getRange(`A1:${col(widths.length)}${end}`).format={font:{name:'Arial',size:10,color:'#243746'},rowHeight:21,verticalAlignment:'center'};
  widths.forEach((w,i)=>sh.getRange(`${col(i+1)}1:${col(i+1)}${end}`).format.columnWidth=w);
  sh.getRange('A2').values=[[title]];
  sh.getRange('A2').format.font={name:'Arial',size:15,bold:true,color:'#183F62'};
}
function table(sh,start,headers,rows){
  sh.getRange(`A${start}:${col(headers.length)}${start+rows.length}`).values=[headers,...rows.map(r=>r.map(value))];
  sh.getRange(`A${start}:${col(headers.length)}${start}`).format={fill:'#183F62',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center',rowHeight:45};
  for(let i=start+1;i<=start+rows.length;i+=2)sh.getRange(`A${i}:${col(headers.length)}${i}`).format.fill='#F2F6FA';
}

const results=[...d.results].sort((a,b)=>a.assay.localeCompare(b.assay)||(a.p??2)-(b.p??2)||a.gene.localeCompare(b.gene));
const tests=sheets['Gene tests'];
const fields=['gene','assay','status','n','n_rearranged','n_unrearranged','reference_genes','unique_scores','beta','se','ci_low','ci_high','lambda_ml','p','q','df','t_stat','intercept','gls_r2','score_median_rearranged','score_median_unrearranged','ols_beta','ols_p','bm_beta','bm_p','wilcoxon_p','wilcoxon_q','log_likelihood'];
base(tests,results.length+9,[17,20,59,12,16,18,16,16,17,17,17,17,16,17,17,12,17,17,16,21,21,17,17,17,17,17,17,20,110],'PGLS rank score by rearrangement group');
tests.getRange('A3').values=[[shared?'Primary model: rank score ~ rearranged, with one pooled Pagel lambda shared by all genes within each assay.':'Primary model: rank score ~ rearranged, with a separate maximum-likelihood Pagel lambda for each gene and assay.']];
tests.getRange('A4').values=[['Positive coefficient = higher rank score in rearranged species. Confidence intervals and p-values use the GLS t approximation.']];
tests.getRange('A5').values=[['Saved statistical results. BH correction is separate by assay; constant scores have no estimable test. Fixed-lambda columns are sensitivity checks.']];
table(tests,7,['Gene','Assay','Fit status','Species n','Rearranged n','Unrearranged n','Shared genes','Distinct scores','Group coefficient','Standard error','95% CI lower','95% CI upper',lambdaLabel,'PGLS p','PGLS BH q','Residual df','t statistic','Intercept','GLS R²','Rearranged median score','Unrearranged median score','OLS coefficient','OLS p','Lambda 1 coefficient','Lambda 1 p','Wilcoxon p','Wilcoxon BH q','Log likelihood','Shared ranking genes'],results.map(r=>[...fields.map(k=>k==='assay'?assay(r[k]):r[k]),d.panels[r.assay+'|'+r.gene].join(', ')||'n.a.']));
tests.getRange(`D8:AB${results.length+7}`).format.horizontalAlignment='right';
tests.getRange(`I8:AB${results.length+7}`).setNumberFormat('0.000000');
tests.getRange(`P8:P${results.length+7}`).setNumberFormat('0');
tests.getRange(`M8:M${results.length+7}`).setNumberFormat('0.000');
tests.getRange(`O8:O${results.length+7}`).conditionalFormats.add('cellIs',{operator:'lessThan',formula:.05,format:{fill:'#D8EADD',font:{bold:true}}});
tests.getRange(`AC8:AC${results.length+7}`).format.wrapText=true;
for(let i=0;i<results.length;i++){
  let n=d.panels[results[i].assay+'|'+results[i].gene].join(', ').length;
  tests.getRange(`A${i+8}:AC${i+8}`).format.rowHeight=Math.max(21,Math.ceil(n/100)*14+6);
  if(results[i].gene==='SIRT6')tests.getRange(`A${i+8}:AB${i+8}`).format.font.bold=true;
}
tests.freezePanes.freezeRows(7);tests.freezePanes.freezeColumns(2);
if(shared){
  tests.getRange(`AD7:AF${results.length+7}`).values=[['Prior gene-specific lambda','Prior gene-specific p','Prior gene-specific BH q'],...results.map(r=>[r.individual_lambda,r.individual_p,r.individual_q].map(value))];
  tests.getRange(`AD7:AF${results.length+7}`).format={font:{name:'Arial',size:10,color:'#243746'},verticalAlignment:'center',columnWidth:22};
  tests.getRange('AD7:AF7').format={fill:'#183F62',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center'};
  tests.getRange(`AD8:AF${results.length+7}`).setNumberFormat('0.000000');
}

const summary=sheets.Summary;
base(summary,110,[24,18,18,18,18,18,18,18,18,18],shared?'Rank scores with a shared phylogenetic parameter':'Rank scores after phylogenetic adjustment');
summary.tabColor='#183F62';
summary.getRange('A3').values=[['10 rearranged and 20 unrearranged species; muskrat, human, macaque and bowhead whale excluded. Missing values remain unavailable.']];
summary.getRange('A4').values=[['Rank score: 1 = highest, 0 = lowest. PGLS tests the group effect with estimated Pagel lambda; no body-mass or MLS predictor.']];
if(shared)summary.getRange('A5').values=[[`Shared lambda: RNA = ${d.shared_lambda.rna.lambda_value.toFixed(4)}; ATAC = ${d.shared_lambda.atac.lambda_value.toFixed(4)}. Lambda = 0 gives ordinary regression with no phylogenetic correction.`]];
table(summary,6,['Assay','PGLS fits','Constant scores','Nominal p < 0.05','PGLS BH q < 0.05','Prior Wilcoxon BH q < 0.05'],[]);
const last=results.length+7;
for(let i=0;i<2;i++){
  const row=7+i,a=['rna','atac'][i];summary.getRange(`A${row}`).values=[[assay(a)]];
  summary.getRange(`B${row}:F${row}`).formulas=[[
    `=COUNTIFS('Gene tests'!$B$8:$B$${last},A${row},'Gene tests'!$C$8:$C$${last},"Tested")`,
    `=COUNTIFS('Gene tests'!$B$8:$B$${last},A${row},'Gene tests'!$C$8:$C$${last},"Not estimable: constant rank score")`,
    `=COUNTIFS('Gene tests'!$B$8:$B$${last},A${row},'Gene tests'!$N$8:$N$${last},"<0.05",'Gene tests'!$C$8:$C$${last},"Tested")`,
    `=COUNTIFS('Gene tests'!$B$8:$B$${last},A${row},'Gene tests'!$O$8:$O$${last},"<0.05",'Gene tests'!$C$8:$C$${last},"Tested")`,
    `=COUNTIFS('Gene tests'!$B$8:$B$${last},A${row},'Gene tests'!$AA$8:$AA$${last},"<0.05")`
  ]];
}
summary.getRange('A10').values=[['SIRT6']];summary.getRange('A10').format.font.bold=true;
table(summary,11,['Assay','Species n','Shared genes','Group coefficient','95% CI lower','95% CI upper',lambdaLabel,'PGLS p','PGLS BH q','Prior Wilcoxon p'],[]);
for(let i=0;i<2;i++){
  const a=['rna','atac'][i],sourceRow=results.findIndex(r=>r.assay===a&&r.gene==='SIRT6')+8,row=12+i;
  summary.getRange(`A${row}:J${row}`).formulas=[['B','D','G','I','K','L','M','N','O','Z'].map(c=>`='Gene tests'!${c}${sourceRow}`)];
}
summary.getRange('D12:J13').setNumberFormat('0.0000');summary.getRange('B12:J13').format.horizontalAlignment='right';
const sirt6Atac=results.find(r=>r.assay==='atac'&&r.gene==='SIRT6');
summary.getRange('A15').values=[[`SIRT6 ATAC group effect = ${sirt6Atac.beta.toFixed(4)}; p = ${sirt6Atac.p.toFixed(4)}; BH q = ${sirt6Atac.q.toFixed(4)}. ${shared?'Shared lambda is pooled across genes.':'Lambda is estimated for SIRT6.'}`]];
summary.getRange('A16').values=[['RNA uses only three shared genes and two observed score levels. Gaussian PGLS inference on these tied scores is exploratory.']];
summary.getRange('A17').values=[['The ten rearranged species form one clade. The group effect cannot establish causation by an inversion.']];
summary.getRange('A19').values=[['Genes with PGLS BH q < 0.05']];summary.getRange('A19').format.font.bold=true;
const significant=results.filter(r=>r.q!==null&&r.q<.05);
table(summary,20,['Gene','Assay','Direction in rearranged','Group coefficient','95% CI lower','95% CI upper',lambdaLabel,'PGLS p','PGLS BH q'],[]);
for(let i=0;i<significant.length;i++){
  let r=significant[i],row=21+i,sourceRow=results.indexOf(r)+8;
  summary.getRange(`A${row}:B${row}`).formulas=[[`='Gene tests'!A${sourceRow}`,`='Gene tests'!B${sourceRow}`]];
  summary.getRange(`C${row}`).formulas=[[`=IF(D${row}>0,"Higher","Lower")`]];
  summary.getRange(`D${row}:I${row}`).formulas=[['I','K','L','M','N','O'].map(c=>`='Gene tests'!${c}${sourceRow}`)];
}
summary.getRange(`D21:I${20+significant.length}`).setNumberFormat('0.0000');
summary.getRange(`H21:I${20+significant.length}`).setNumberFormat('0.000000');
const methodStart=Math.max(32,significant.length+23);
summary.getRange(`A${methodStart}`).values=[['Analysis definitions and limitations']];summary.getRange(`A${methodStart}`).format.font.bold=true;
let nr=methodStart+2;
for(const note of d.notes){
  let line='';for(const w of note.split(' ')){if((line+' '+w).length>155){summary.getRange(`A${nr++}`).values=[[line]];line=w;}else line+=(line?' ':'')+w;}
  if(line)summary.getRange(`A${nr++}`).values=[[line]];nr++;
}

const raw=sheets['Rank scores'];
base(raw,d.measurements.length+12,[20,17,32,30,19,20,20,18,18,18,59],'Source medians, shared ranks and rank scores');
raw.getRange('A3').values=[['Source: '+d.source]];
raw.getRange('A4').values=[['Original medians: RNA source liver log.tpm; ATAC integrated TSS ±5 kb CPM·bp. Blank score inputs stay unavailable.']];
raw.getRange('A5').values=[['Rank score is calculated here from shared-panel rank. PGLS results are saved external fits; changing a cell does not refit PGLS.']];
const resultByKey=new Map(d.results.map(r=>[r.assay+'|'+r.gene,r]));
table(raw,7,['Assay','Gene','Species','Scientific name','Group','Source median','Shared-panel rank','Shared genes','Rearranged indicator','Rank score','Gene fit status'],d.measurements.map(r=>[assay(r.assay),r.gene,r.species,r.latin,r.group,r.value,r.rank,r.reference_genes,r.rearranged,null,resultByKey.get(r.assay+'|'+r.gene).status]));
const scoreFormulas=d.measurements.map((r,i)=>[`=IF(AND(ISNUMBER(G${i+8}),H${i+8}>1),1-(G${i+8}-1)/(H${i+8}-1),"n.a.")`]);
raw.getRange(`J8:J${d.measurements.length+7}`).formulas=scoreFormulas;
raw.getRange(`F8:F${d.measurements.length+7}`).setNumberFormat('0.0000');
raw.getRange(`G8:G${d.measurements.length+7}`).setNumberFormat('0.0');
raw.getRange(`J8:J${d.measurements.length+7}`).setNumberFormat('0.0000');
raw.getRange(`F8:J${d.measurements.length+7}`).format.horizontalAlignment='right';
raw.freezePanes.freezeRows(7);raw.freezePanes.freezeColumns(3);

const species=sheets['Species groups'];
const order={Rearranged:0,Unrearranged:1,Excluded:2};
const spp=[...d.species].sort((a,b)=>order[a.group]-order[b.group]||a.common.localeCompare(b.common));
base(species,55,[33,31,20,56,90,23,80],'Species groups and assignment provenance');
species.getRange('A3').values=[['Sources: local annotation-based gene order and user assignments on 31 October 2026. Muskrat is retained only for documentation.']];
species.getRange('A4').values=[['No species substitutions in TimeTree. Rearranged species form one clade among the included taxa.']];
table(species,6,['Species','Scientific name','Analysis group','Assignment source','Assignment reason','Annotation-only group','Annotation-only reason'],spp.map(s=>[s.common,s.latin,s.group,s.evidence,s.reason,s.annotation_group,s.annotation_reason]));
species.getRange('A40').values=[['Covariance source: '+d.references[0]]];
species.getRange('A41').values=[['Local tree: '+d.tree_file]];
species.getRange('A42').values=[['Model reference: '+d.references[1]]];
species.getRange('A43').values=[['Lambda reference: '+d.references[2]]];
species.getRange('A44').values=[['Measurements source: '+d.measurement_source]];
species.getRange('A45').values=[['Source SHA-256: '+d.measurement_sha256]];
species.freezePanes.freezeRows(6);species.freezePanes.freezeColumns(3);

wb.recalculate();
const expected=['rna','atac'].map(a=>[d.summary[a].tested,d.summary[a].constant,d.summary[a].nominal.length,d.summary[a].fdr.length,d.summary[a].prior_wilcoxon_fdr.length]);
if(JSON.stringify(summary.getRange('B7:F8').values)!==JSON.stringify(expected))throw Error('Summary mismatch');
const scores=raw.getRange(`J8:J${d.measurements.length+7}`).values;
d.measurements.forEach((r,i)=>{if(r.rank_score===null?scores[i][0]!=='n.a.':Math.abs(scores[i][0]-r.rank_score)>1e-12)throw Error('Rank score mismatch at row '+(i+8));});
for(let i=0;i<2;i++){
  const r=results.find(r=>r.assay===['rna','atac'][i]&&r.gene==='SIRT6');
  if(Math.abs(summary.getRange(`H${i+12}`).values[0][0]-r.p)>1e-12)throw Error('SIRT6 mismatch');
}
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:8},maxChars:800})).ndjson);
console.log((await wb.inspect({kind:'table',range:'Summary!A11:J13',include:'values,formulas',tableMaxRows:3,tableMaxCols:10,maxChars:2500})).ndjson);
for(const [name,range] of [['Summary',`A1:J${Math.max(29,20+significant.length)}`],['Gene tests','A1:O15'],['Rank scores','A1:K14'],['Species groups','A1:E16']]){
  if(process.argv.includes('--summary-preview')&&name!=='Summary')continue;
  if(shared&&['Rank scores','Species groups'].includes(name))continue;
  const p=await wb.render({sheetName:name,range,scale:1.25,format:'png'});
  await fs.writeFile(path.join(out,'build',name.replaceAll(' ','_')+'.png'),new Uint8Array(await p.arrayBuffer()));
}
await(await SpreadsheetFile.exportXlsx(wb)).save(path.join(out,shared?'PGLS_shared_lambda_rank_scores.xlsx':'PGLS_rearrangement_rank_scores.xlsx'));
console.log('Workbook saved; all rank-score formulas and summary counts reconciled.');
