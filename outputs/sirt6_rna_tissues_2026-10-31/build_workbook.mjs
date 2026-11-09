import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const out=path.dirname(fileURLToPath(import.meta.url));
const d=JSON.parse(await fs.readFile(path.join(out,'analysis.json'),'utf8'));
const wb=Workbook.create();
const sheets=Object.fromEntries(['Summary','Gene tests','Rank scores','Species groups','Sources'].map(n=>[n,wb.worksheets.add(n)]));
const col=n=>{let s='';for(;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
const value=x=>x??'n.a.';
function base(sh,end,widths,title){
  sh.showGridLines=false;
  sh.getRange(`A1:${col(widths.length)}${end}`).format={font:{name:'Arial',size:10,color:'#243746'},rowHeight:21,verticalAlignment:'center'};
  widths.forEach((w,i)=>sh.getRange(`${col(i+1)}1:${col(i+1)}${end}`).format.columnWidth=w);
  sh.getRange('A2').values=[[title]];
  sh.getRange('A2').format.font={name:'Arial',size:15,bold:true,color:'#183F62'};
}
function table(sh,start,headers,rows,{stripes=false}={}){
  sh.getRange(`A${start}:${col(headers.length)}${start+rows.length}`).values=[headers,...rows.map(r=>r.map(value))];
  sh.getRange(`A${start}:${col(headers.length)}${start}`).format={fill:'#183F62',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center',rowHeight:42};
  if(stripes)for(let i=start+1;i<=start+rows.length;i+=2)sh.getRange(`A${i}:${col(headers.length)}${i}`).format.fill='#F2F6FA';
}

const results=[...d.results].sort((a,b)=>d.tissues.indexOf(a.tissue)-d.tissues.indexOf(b.tissue)||(a.p??2)-(b.p??2)||a.gene.localeCompare(b.gene));
const fields=['tissue','gene','status','n','n_rearranged','n_unrearranged','reference_genes','unique_scores','beta','se','ci_low','ci_high','lambda_ml','p','q_tissue','q_global','q_sirt6_tissues','df','intercept','median_rearranged','median_unrearranged','bm_beta','bm_p'];
const gt=sheets['Gene tests'];
const headers=['Tissue','Gene','Fit status','Species n','Rearranged n','Unrearranged n','Shared ranking genes n','Distinct scores','Group coefficient','Standard error','95% CI lower','95% CI upper','Shared lambda','p','BH q within tissue','BH q across tissues','SIRT6-only q across tissues','Residual df','Intercept','Rearranged median score','Unrearranged median score','Lambda 1 coefficient','Lambda 1 p','Shared ranking gene names'];
base(gt,results.length+9,[17,17,56,13,16,18,19,17,18,18,18,18,16,18,19,19,22,14,17,23,23,20,20,95],'RNA rank-score tests by tissue');
gt.getRange('A3').values=[['Model: rank score ~ rearranged. One lambda is shared across genes within each tissue. Positive coefficient = higher rank score in rearranged species.']];
gt.getRange('A4').values=[['Saved external PGLS results. BH across tissues covers all 120 valid gene–tissue tests; constant and unavailable outcomes are not tested.']];
gt.getRange('A5').values=[['The same TimeTree and species assignments are used throughout. Tissue and target coverage can differ. Lambda 1 is a sensitivity check.']];
table(gt,7,headers,results.map(r=>[...fields.map(k=>r[k]),d.panels[r.tissue+'|'+r.gene].join(', ')||'n.a.']));
const end=results.length+7;
gt.getRange(`D8:W${end}`).format.horizontalAlignment='right';
gt.getRange(`I8:W${end}`).setNumberFormat('0.0000');
gt.getRange(`N8:Q${end}`).setNumberFormat('0.000000');
gt.getRange(`R8:R${end}`).setNumberFormat('0');
gt.getRange(`X8:X${end}`).format.wrapText=true;
for(let i=0;i<results.length;i++){
  const r=results[i];
  if(r.gene==='SIRT6')gt.getRange(`A${i+8}:W${i+8}`).format.font.bold=true;
  const length=d.panels[r.tissue+'|'+r.gene].join(', ').length;
  if(length>85)gt.getRange(`A${i+8}:X${i+8}`).format.rowHeight=Math.ceil(length/85)*14+6;
}
gt.getRange(`O8:P${end}`).conditionalFormats.add('cellIs',{operator:'lessThan',formula:.05,format:{fill:'#D8EADD',font:{bold:true}}});
gt.freezePanes.freezeRows(7);gt.freezePanes.freezeColumns(2);

const s=sheets.Summary;
base(s,67,[21,17,18,17,18,18,19,19,17,49],'RNA expression ranks across tissues');
s.tabColor='#183F62';
s.getRange('A3').values=[['SIRT6-locus rearranged versus unrearranged species. Human, macaque, bowhead whale and muskrat remain excluded.']];
s.getRange('A4').values=[['No gene passes FDR < 0.05 within a tissue or across tissues. Shared lambda = 0 in every estimable tissue, so these fits apply no phylogenetic correction.']];
s.getRange('A6').values=[['SIRT6']];s.getRange('A6').format.font.bold=true;
table(s,7,['Tissue','Rearranged n','Unrearranged n','Shared genes n','Group coefficient','p','BH q within tissue','BH q across tissues','Shared lambda','Fit status'],[]);
for(let i=0;i<d.tissues.length;i++){
  const tissue=d.tissues[i],rr=results.findIndex(r=>r.tissue===tissue&&r.gene==='SIRT6')+8,row=i+8;
  s.getRange(`A${row}:J${row}`).formulas=[['A','E','F','G','I','N','O','P','M','C'].map(c=>`='Gene tests'!${c}${rr}`)];
}
s.getRange('E8:I14').setNumberFormat('0.0000');
s.getRange('B8:I14').format.horizontalAlignment='right';
s.getRange('A16').values=[['Brain and heart: SIRT6 ranks 2 of 3 in every measured species. Liver, kidney, lung and skin: no significant SIRT6 rank difference.']];
s.getRange('A17').values=[['The three shared SIRT6 reference genes are SIRT6, CHAF1A and HDGFL2. Missing data substantially limit ranking resolution.']];
s.getRange('A18').values=[['Spleen has only beaver data and cannot support a group test. A constant rank is not evidence that absolute RNA expression is identical.']];
s.getRange('A20').values=[['Tissue coverage and locus-wide tests']];s.getRange('A20').format.font.bold=true;
table(s,21,['Tissue','Species with locus data','Samples with locus data','Valid gene tests','Constant scores','Nominal p < 0.05','BH q within tissue < 0.05','BH q across tissues < 0.05','Shared lambda','Coverage note'],[]);
for(let i=0;i<d.summary.length;i++){
  const z=d.summary[i],row=i+22;
  s.getRange(`A${row}:C${row}`).values=[[z.tissue,z.species,z.samples]];
  s.getRange(`D${row}:H${row}`).formulas=[[
    `=COUNTIFS('Gene tests'!$A$8:$A$${end},A${row},'Gene tests'!$C$8:$C$${end},"Tested")`,
    `=COUNTIFS('Gene tests'!$A$8:$A$${end},A${row},'Gene tests'!$C$8:$C$${end},"Not estimable: constant rank score")`,
    ...['N','O','P'].map(c=>`=COUNTIFS('Gene tests'!$A$8:$A$${end},A${row},'Gene tests'!$${c}$8:$${c}$${end},"<0.05",'Gene tests'!$C$8:$C$${end},"Tested")`)
  ]];
  s.getRange(`I${row}:J${row}`).values=[[value(z.lambda_value),z.tested?'Coverage varies by target gene':'Only one unrearranged species']];
}
s.getRange('B22:I28').format.horizontalAlignment='right';s.getRange('I22:I28').setNumberFormat('0.0000');
s.getRange('A30').values=[['Nominal associations; none survives FDR correction']];s.getRange('A30').format.font.bold=true;
table(s,31,['Tissue','Gene','Direction in rearranged','Group coefficient','p','BH q within tissue','BH q across tissues'],[]);
const nominal=results.filter(r=>r.p!==null&&r.p<.05);
for(let i=0;i<nominal.length;i++){
  const rr=results.indexOf(nominal[i])+8,row=i+32;
  s.getRange(`A${row}:B${row}`).formulas=[[`='Gene tests'!A${rr}`,`='Gene tests'!B${rr}`]];
  s.getRange(`C${row}`).formulas=[[`=IF(D${row}>0,"Higher","Lower")`]];
  s.getRange(`D${row}:G${row}`).formulas=[['I','N','O','P'].map(c=>`='Gene tests'!${c}${rr}`)];
}
s.getRange(`D32:G${31+nominal.length}`).setNumberFormat('0.0000');
s.getRange('A38').values=[['Interpretation']];s.getRange('A38').format.font.bold=true;
[
  'Rank score ranges from 0 (lowest) to 1 (highest). Positive coefficients indicate higher relative expression ranks in rearranged species.',
  'Each gene uses a reference set observed in all species included in that test. Sparse measurements can leave very few genes to rank.',
  'Shared lambda is estimated separately by tissue from summed gene-level profile likelihoods. Lambda = 0 gives ordinary regression in this tree.',
  'All rearranged species form one clade. These tests cannot separate inversion effects from other changes inherited along that lineage.',
  'The bounded and often tied scores make Gaussian PGLS inference exploratory. An absent rank difference does not rule out an expression difference.',
  'Full methods, source-file provenance and the preserved species assignments are provided in the supporting sheets.'
].forEach((txt,i)=>s.getRange(`A${40+i}`).values=[[txt]]);

const measurements=d.measurements.filter(r=>r.median_log_tpm!==null);
const rs=sheets['Rank scores'];
base(rs,measurements.length+9,[17,17,34,31,21,22,18,17,20,18,18,17],'Species medians and within-species ranks');
rs.getRange('A3').values=[['Source: CrossSpecies_Count.xlsx with missing mouse records from rnaLongTPM.xlsx. Species medians are source log.tpm without another logarithm.']];
rs.getRange('A4').values=[['Only measured gene–species–tissue records are listed. Missing values are never zero. Samples are aggregated within species and tissue.']];
rs.getRange('A5').values=[['The score formula is shown here; editing cells does not rerun the external PGLS fits. Shared ranking panels are listed in Gene tests.']];
table(rs,7,['Tissue','Gene','Species','Scientific name','Group','Median source log.tpm','Samples n','Shared rank','Shared genes n','Rank score','Rearranged indicator','Used in test'],measurements.map(r=>[r.tissue,r.gene,r.species,r.latin,r.group,r.median_log_tpm,r.n_samples,r.rank,r.reference_genes,null,r.rearranged,r.used_in_test]));
const mEnd=measurements.length+7;
rs.getRange(`J8:J${mEnd}`).formulas=measurements.map((r,i)=>[`=IF(AND(ISNUMBER(H${i+8}),I${i+8}>1),1-(H${i+8}-1)/(I${i+8}-1),"n.a.")`]);
rs.getRange(`F8:F${mEnd}`).setNumberFormat('0.0000');rs.getRange(`H8:H${mEnd}`).setNumberFormat('0.0');rs.getRange(`J8:J${mEnd}`).setNumberFormat('0.0000');
rs.getRange(`F8:K${mEnd}`).format.horizontalAlignment='right';rs.freezePanes.freezeRows(7);rs.freezePanes.freezeColumns(3);

const sg=sheets['Species groups'];
base(sg,48,[34,32,21,57,86],'Preserved species groups');
sg.getRange('A3').values=[['Ten rearranged and twenty unrearranged species. Muskrat is retained below for documentation but excluded from all tests.']];
sg.getRange('A4').values=[['Groups use the previously inferred marker order plus the explicit user assignments. Missing tissue data remain unavailable.']];
table(sg,6,['Species','Scientific name','Group','Assignment source','Assignment reason'],d.species.map(r=>[r.common,r.latin,r.group,r.evidence,r.reason]),{stripes:true});
sg.freezePanes.freezeRows(6);sg.freezePanes.freezeColumns(3);

const src=sheets.Sources;
base(src,d.samples.length+95,[27,23,18,19,34,31,21,31,23,23,18],'Source workbooks and sample coverage');
table(src,4,['Source workbook','Data rows','Size in bytes','SHA-256'],d.sources.map(z=>[z.workbook,z.total_rows,z.bytes,z.sha256]));
src.getRange('D5:D6').format.columnWidth=72;
src.getRange('A8').values=[['Source overlap: 9,061 locus records matched across workbooks; all duplicate TPM and log.tpm values agree. Each measurement is counted once.']];
src.getRange('A9').values=[['BMR_Turk_Lu is analyzed as lung, following the earlier confirmed tissue correction. The original tissue label is retained below.']];
let noteRow=12;
src.getRange(`A${noteRow++}`).values=[['Methods and limitations']];src.getRange('A12').format.font.bold=true;
for(const note of d.notes){
  let line='';
  for(const word of note.split(' ')){
    if((line+' '+word).length>150){src.getRange(`A${noteRow++}`).values=[[line]];line=word;}
    else line+=(line?' ':'')+word;
  }
  if(line)src.getRange(`A${noteRow++}`).values=[[line]];
  noteRow++;
}
for(const url of [...new Set(d.references)])src.getRange(`A${noteRow++}`).values=[['Reference: '+url]];
noteRow+=2;
const sampleHeader=noteRow;
table(src,sampleHeader,['Workbook','Sample','Analysis tissue','Source tissue','Species','Scientific name','Group','Source study','Source rows','Selected locus genes','Used in analysis'],d.samples.map(r=>[r.workbook,r.sample,r.tissue,r.source_tissue,r.species,r.latin,r.group,r.source,r.total_source_rows,r.selected_locus_genes,r.used]));
src.getRange(`D${sampleHeader}:D${sampleHeader+d.samples.length}`).format.columnWidth=20;
src.getRange('D5:D6').format.wrapText=true;
src.getRange('A5:D6').format.rowHeight=60;
src.getRange('B5:C6').setNumberFormat('#,##0');

wb.recalculate();
for(let i=0;i<d.summary.length;i++){
  const z=d.summary[i],expected=[z.tested,z.constant,z.nominal,z.fdr_tissue,z.fdr_global];
  if(JSON.stringify(s.getRange(`D${i+22}:H${i+22}`).values[0])!==JSON.stringify(expected))throw Error('Tissue summary mismatch: '+z.tissue);
  const r=d.results.find(r=>r.tissue===z.tissue&&r.gene==='SIRT6');
  const got=s.getRange(`F${i+8}`).values[0][0];
  if(r.p===null?got!=='n.a.':Math.abs(got-r.p)>1e-12)throw Error('SIRT6 p mismatch');
}
const calculated=rs.getRange(`J8:J${mEnd}`).values;
measurements.forEach((r,i)=>{if(r.rank_score===null?calculated[i][0]!=='n.a.':Math.abs(calculated[i][0]-r.rank_score)>1e-12)throw Error('Rank formula mismatch '+i);});
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:10},maxChars:1000})).ndjson);
console.log((await wb.inspect({kind:'table',range:'Summary!A7:J14',include:'values,formulas',tableMaxRows:8,tableMaxCols:10,maxChars:2400})).ndjson);
await fs.mkdir(path.join(out,'build'),{recursive:true});
for(const [name,range,file] of [['Summary','A1:J45','Summary'],['Gene tests','A1:P13','Gene_tests'],['Rank scores','A1:L13','Rank_scores'],['Species groups','A1:E15','Species_groups'],['Sources','A1:D14','Sources'],['Sources',`A${sampleHeader}:K${sampleHeader+7}`,'Source_samples']]){
  if(process.argv.includes('--sources-preview')&&file!=='Sources')continue;
  const p=await wb.render({sheetName:name,range:file==='Sources'?'A1:J14':range,scale:1.15,format:'png'});
  await fs.writeFile(path.join(out,'build',file+'.png'),new Uint8Array(await p.arrayBuffer()));
}
await(await SpreadsheetFile.exportXlsx(wb)).save(path.join(out,'RNA_tissue_rank_comparison.xlsx'));
console.log('Workbook saved. All tissue counts, SIRT6 p-values and rank-score formulas reconciled.');
