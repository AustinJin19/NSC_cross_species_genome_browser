import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const out=path.dirname(fileURLToPath(import.meta.url));
const d=JSON.parse(await fs.readFile(path.join(out,'analysis.json'),'utf8'));
const wb=Workbook.create();
const names=['Summary','Test statistics','Species ranks','Pairwise correlations','Species coverage','Methods'];
const sh=Object.fromEntries(names.map(n=>[n,wb.worksheets.add(n)]));
const col=n=>{let s='';for(;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
const val=x=>x??'Not tested';
function base(s,rows,widths,title){
  s.showGridLines=false;
  s.getRange(`A1:${col(widths.length)}${rows}`).format={font:{name:'Arial',size:10,color:'#243746'},rowHeight:22,verticalAlignment:'center'};
  widths.forEach((w,i)=>s.getRange(`${col(i+1)}1:${col(i+1)}${rows}`).format.columnWidth=w);
  s.getRange('A2').values=[[title]];
  s.getRange('A2').format.font={name:'Arial',size:16,bold:true,color:'#183F62'};
}
function table(s,row,headers,rows){
  s.getRange(`A${row}:${col(headers.length)}${row+rows.length}`).values=[headers,...rows.map(r=>r.map(val))];
  s.getRange(`A${row}:${col(headers.length)}${row}`).format={fill:'#183F62',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},rowHeight:43,wrapText:true,horizontalAlignment:'center'};
  if(rows.length<50)for(let r=row+1;r<=row+rows.length;r+=2)s.getRange(`A${r}:${col(headers.length)}${r}`).format.fill='#F0F5F9';
}
function text(s,row,content,end='I',height=32){
  const r=s.getRange(`A${row}:${end}${row}`);r.merge();r.values=[[content]];r.format.wrapText=true;r.format.rowHeight=height;
}
function freeze(s,row,columns=2){s.freezePanes.freezeRows(row);s.freezePanes.freezeColumns(columns);}

const t=sh['Test statistics'];
base(t,27,[23,20,12,12,14,18,12,18,18,18,18,41,20,19,20,17,17,15,39,95],'Within-group concordance tests');
text(t,3,'Species are blocks; genes are ranked outcomes. Both groups use the same complete gene panel within each dataset.','L');
text(t,4,'Primary p-values use exact rank-label permutations for RNA and 200,000 Monte Carlo permutations for ATAC. BH correction covers 14 valid tests.','L');
text(t,5,'ATAC: zero of 200,000 permutations exceeded the observed statistic. The reported Monte Carlo p-value is its simulation resolution, not an exact p-value.','L');
const fields=['dataset','group','n_species','n_genes','kendall_w','friedman_chi2','df','p_chi2','q_chi2','p_permutation','q_permutation','permutation_method','permutations','exceedances','mc_upper95','seed','tie_correction','all_measured_species','status'];
table(t,7,['Dataset','Group','Species n','Genes n','Kendall W','Friedman chi-square','df','Asymptotic p','Asymptotic BH q','Permutation p','Permutation BH q','Permutation method','Permutations','Exceedances','MC p 95% upper bound','Random seed','Tie correction','Species with any locus data','Status','Primary panel genes'],d.summary.map(z=>[...fields.map(k=>z[k]),d.panels[z.dataset].join('; ')]));
t.getRange('E8:F23').setNumberFormat('0.0000');t.getRange('H8:K23').setNumberFormat('0.000E+00');t.getRange('O8:O23').setNumberFormat('0.000E+00');
t.getRange('Q8:Q23').setNumberFormat('0.0000');t.getRange('C8:R23').format.horizontalAlignment='right';
t.getRange('K8:K23').conditionalFormats.add('cellIs',{operator:'lessThan',formula:.05,format:{fill:'#DDEEE3'}});
freeze(t,7);

const s=sh.Summary;
base(s,65,[25,15,15,13,16,16,19,19,42],'Agreement in SIRT6-locus gene rankings');s.tabColor='#183F62';
text(s,3,'Within-group agreement: rearranged versus unrearranged species, using the previous SIRT6 comparison subsets.');
text(s,4,'All 14 testable comparisons have permutation BH q < 0.05. Kendall W describes the strength of agreement; 1 means identical rankings.', 'I',36);
text(s,5,'RNA uses only SIRT6, CHAF1A and HDGFL2; ATAC uses 30 shared genes. These are conventional tests without phylogenetic correction.', 'I',36);
table(s,7,['Dataset','Rearranged n','Unrearranged n','Shared genes n','Rearranged W','Unrearranged W','Rearranged BH q','Unrearranged BH q','Interpretation'],[]);
for(let i=0;i<d.datasets.length;i++){
  const ds=d.datasets[i],a=d.summary.findIndex(z=>z.dataset===ds&&z.group==='Rearranged')+8,b=d.summary.findIndex(z=>z.dataset===ds&&z.group==='Unrearranged')+8,r=i+8;
  s.getRange(`A${r}:H${r}`).formulas=[[`='Test statistics'!A${a}`,`='Test statistics'!C${a}`,`='Test statistics'!C${b}`,`='Test statistics'!D${a}`,`='Test statistics'!E${a}`,`='Test statistics'!E${b}`,`='Test statistics'!K${a}`,`='Test statistics'!K${b}`]];
  s.getRange(`I${r}`).values=[[ds==='RNA Spleen'?'Not tested: only one species':'Within-group concordance detected']];
}
s.getRange('E8:F15').setNumberFormat('0.000');s.getRange('G8:H15').setNumberFormat('0.000E+00');s.getRange('B8:H15').format.horizontalAlignment='right';
text(s,17,'What the result means', 'I',24);s.getRange('A17').format.font.bold=true;
text(s,18,'Species within each group tend to preserve the gene hierarchy. Significance does not mean every species has exactly the same ranks.');
text(s,19,'Separate significant within-group tests do not establish that the two groups have the same ranking pattern or test an effect of the rearrangement.');
text(s,20,'RNA liver consensus in both groups: HDGFL2 > SIRT6 > CHAF1A. This describes only three genes and cannot characterize the full 63-gene locus.');
text(s,21,'All rearranged species form one clade among the included taxa. Shared ancestry can contribute to agreement; these p-values assume independent species.');
text(s,23,'Broader pairwise similarity using all available shared locus genes', 'I',26);s.getRange('A23').format.font.bold=true;
text(s,24,'Descriptive medians of pairwise Spearman correlations. Gene sets differ between pairs; pairs are not independent and no pooled significance test is applied.');
table(s,26,['Dataset','Rearranged median rho','Unrearranged median rho','Between-group median rho','Rearranged pairs','Unrearranged pairs','Between-group pairs','Shared genes range','Use'],d.datasets.map(ds=>{
  const zz=['Rearranged','Unrearranged','Between groups'].map(g=>d.pair_summary.find(z=>z.dataset===ds&&z.comparison===g));
  const valid=zz.filter(z=>z.min_shared_genes!=null);
  return [ds,...zz.map(z=>z.median_rho),...zz.map(z=>z.pairs),valid.length?`${Math.min(...valid.map(z=>z.min_shared_genes))}–${Math.max(...valid.map(z=>z.max_shared_genes))}`:'Not tested','Descriptive only'];
}));
s.getRange('B27:D34').setNumberFormat('0.000');s.getRange('B27:H34').format.horizontalAlignment='right';
text(s,36,'Verification and scope', 'I',24);s.getRange('A36').format.font.bold=true;
s.getRange('A37').values=[['Valid concordance tests']];s.getRange('B37').formulas=[['=COUNTIF(\'Test statistics\'!S8:S23,"Tested")']];
s.getRange('A38').values=[['Permutation BH q < 0.05']];s.getRange('B38').formulas=[['=COUNTIFS(\'Test statistics\'!S8:S23,"Tested",\'Test statistics\'!K8:K23,"<0.05")']];
text(s,40,'Independent verification: R stats::friedman.test reproduces all 14 chi-square statistics and p-values; maximum Kendall W difference is 1.12 × 10^-16.');
text(s,41,'Human, macaque, bowhead whale and muskrat are excluded. Missing measurements are not imputed. Coverage and group assignment provenance are included.');
text(s,42,'Sheets: test statistics; individual species ranks; pairwise correlations; species coverage; methods, source hashes and species assignments.');
freeze(s,7,1);

const p=sh['Species ranks'];
base(p,d.profiles.length+9,[24,21,35,27,19,23,18,13,18,18],'Ranks in the complete shared panels');
text(p,3,'Rank 1 is the highest measurement within that species. Average ranks are used for ties. No values are imputed.','J');
text(p,4,'RNA uses species-median source log.tpm; ATAC uses species-median integrated CPM·bp in TSS ±5 kb. Scores are calculated below; statistics are saved external results.','J',36);
table(p,6,['Dataset','Group','Species','Species ID','Gene','Source median','Source units','Panel genes n','Descending rank','Rank score'],d.profiles.map(z=>[z.dataset,z.group,z.species,z.species_id,z.gene,z.source_median,z.dataset.startsWith('RNA')?'source log.tpm':'CPM·bp',d.panels[z.dataset].length,z.rank,null]));
const pEnd=d.profiles.length+6;
p.getRange(`J7:J${pEnd}`).formulas=d.profiles.map((z,i)=>[`=IF(H${i+7}>1,1-(I${i+7}-1)/(H${i+7}-1),"Not tested")`]);
p.getRange(`F7:F${pEnd}`).setNumberFormat('0.0000');p.getRange(`I7:I${pEnd}`).setNumberFormat('0.0');p.getRange(`J7:J${pEnd}`).setNumberFormat('0.0000');freeze(p,6,3);

const pc=sh['Pairwise correlations'];
base(pc,d.pairs.length+9,[24,35,35,21,21,22,18,18,32,110],'Pairwise species ranking similarity');
text(pc,3,'Spearman rho uses every available shared gene from the original 63-gene locus, with at least three shared genes per pair.','I');
text(pc,4,'These are descriptive correlations using different gene sets across pairs. A species pair is not an independent biological replicate.','I');
table(pc,6,['Dataset','Species A','Species B','Group A','Group B','Comparison','Shared genes n','Spearman rho','Status','Shared gene names'],d.pairs.map(z=>[z.dataset,z.species_a,z.species_b,z.group_a,z.group_b,z.comparison,z.shared_genes,z.spearman_rho,z.status,z.gene_names]));
pc.getRange(`H7:H${d.pairs.length+6}`).setNumberFormat('0.0000');freeze(pc,6,3);

const cv=sh['Species coverage'];
base(cv,d.coverage.length+9,[24,35,31,21,19,18,21,65],'Species coverage and panel eligibility');
text(cv,3,'The primary concordance analysis uses species with SIRT6 measured and the complete panel shared across both groups within that dataset.','H');
text(cv,4,'Keeping all RNA species with any locus data leaves zero genes shared across all unrearranged species; the SIRT6 subset preserves the earlier comparison.','H',36);
table(cv,6,['Dataset','Species','Species ID','Group','Measured locus genes n','SIRT6 measured','Primary panel included','Inclusion or exclusion reason'],d.coverage.map(z=>[z.dataset,z.species,z.species_id,z.group,z.measured_genes,z.sirt6_measured?'Yes':'No',z.primary_profile_included?'Yes':'No',z.reason]));freeze(cv,6,3);

const m=sh.Methods;
base(m,150,[28,38,22,32,57,85],'Methods, provenance and species assignments');
let row=4;
for(const note of d.notes){text(m,row,note,'F',note.length>240?52:38);row+=2;}
text(m,row,'Independent validation: validate_friedman.R and R_friedman_validation.tsv reproduce all 14 valid tests using base R.','F');row+=2;
text(m,row,'The workbook contains saved results from analyze.py. Editing workbook cells does not rerun the permutations or change the source data.','F');row+=3;
table(m,row,['Source analysis file','SHA-256'],d.sources.map(z=>[z.file,z.sha256]));
for(let i=1;i<=d.sources.length;i++){m.getRange(`A${row+i}:B${row+i}`).format.wrapText=true;m.getRange(`A${row+i}:B${row+i}`).format.rowHeight=110;}
row+=d.sources.length+3;
for(const url of d.references){text(m,row,url,'F',30);row++;}
row+=2;
const assignmentHeader=row;
table(m,row,['Species ID','Common name','Group','Scientific name','Assignment source','Assignment reason'],d.species.map(z=>[z.id,z.common,z.group,z.latin,z.evidence,z.reason]));
m.getRange(`E${row+1}:F${row+d.species.length}`).format.wrapText=true;
m.getRange(`A${row+1}:F${row+d.species.length}`).format.rowHeight=43;

wb.recalculate();
if(s.getRange('B37').values[0][0]!==14||s.getRange('B38').values[0][0]!==14)throw Error('Summary counts failed');
for(let i=0;i<d.datasets.length;i++)for(const [group,c] of [['Rearranged','E'],['Unrearranged','F']]){
  const z=d.summary.find(z=>z.dataset===d.datasets[i]&&z.group===group),got=s.getRange(`${c}${i+8}`).values[0][0];
  if(z.kendall_w==null?got!=='Not tested':Math.abs(got-z.kendall_w)>1e-12)throw Error('W mismatch');
}
const scores=p.getRange(`J7:J${pEnd}`).values;
d.profiles.forEach((z,i)=>{if(Math.abs(scores[i][0]-z.rank_score)>1e-12)throw Error('Rank score mismatch');});
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:10},maxChars:1200})).ndjson);
await fs.mkdir(path.join(out,'build'),{recursive:true});
for(const [name,range,file] of [['Summary','A1:I21','Summary'],['Summary','A23:I42','Summary_pairwise'],['Test statistics','A1:L12','Tests'],['Species ranks','A1:J13','Ranks'],['Pairwise correlations','A1:I13','Pairs'],['Species coverage','A1:H13','Coverage'],['Methods','A1:F16','Methods'],['Methods',`A${assignmentHeader}:F${assignmentHeader+7}`,'Assignments']]){
  const img=await wb.render({sheetName:name,range,scale:1.1,format:'png'});
  await fs.writeFile(path.join(out,'build',file+'.png'),new Uint8Array(await img.arrayBuffer()));
}
await(await SpreadsheetFile.exportXlsx(wb)).save(path.join(out,'Rank_concordance_within_groups.xlsx'));
console.log('Saved workbook. Counts, W values and all rank-score formulas verified.');
