import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const out=process.argv[2]?path.resolve(process.argv[2]):path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const d=JSON.parse(await fs.readFile(path.join(out,'analysis.json'),'utf8'));
const wb=Workbook.create();
const names=['SIRT6','All gene tests','Species groups','Measurements','Rank panels'];
const sheets=Object.fromEntries(names.map(n=>[n,wb.worksheets.add(n)]));
const col=n=>{let s='';for(;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
const val=x=>x??'n.a.';
const assay=x=>x==='rna'?'Liver RNA':'ATAC TSS ±5 kb';
function base(sh,end,widths,title){
 sh.showGridLines=false;
 sh.getRange(`A1:${col(widths.length)}${end}`).format={font:{name:'Arial',size:10,color:'#243746'},rowHeight:21,verticalAlignment:'center'};
 widths.forEach((w,i)=>sh.getRange(`${col(i+1)}1:${col(i+1)}${end}`).format.columnWidth=w);
 sh.getRange('A2').values=[[title]];
 sh.getRange('A2').format.font={name:'Arial',size:15,bold:true,color:'#183F62'};
}
function table(sh,start,headers,rows){
 sh.getRange(`A${start}:${col(headers.length)}${start+rows.length}`).values=[headers,...rows.map(r=>r.map(val))];
 sh.getRange(`A${start}:${col(headers.length)}${start}`).format={fill:'#183F62',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center',rowHeight:45};
 for(let row=start+1;row<=start+rows.length;row+=2)sh.getRange(`A${row}:${col(headers.length)}${row}`).format.fill='#F2F6FA';
}

const all=sheets['All gene tests'];
const tests=[...d.results].sort((a,b)=>a.assay.localeCompare(b.assay)||(a.p??2)-(b.p??2)||a.gene.localeCompare(b.gene));
const fields=['gene','assay','status','reference_genes','n_rearranged','n_unrearranged','median_rearranged','median_unrearranged','delta_higher','p','q','u','permutations','original_rank_p','original_rank_q','percentile_p','percentile_q','raw_value_p','raw_value_q'];
base(all,tests.length+8,[16,20,57,15,15,15,18,18,17,15,15,13,18,16,16,16,16,16,16],'Gene-rank tests by gene-order group');
all.tabColor='#183F62';
all.getRange('A3').values=[['Primary test uses the same reference genes within every species entering that target-gene comparison. Rank 1 is highest.']];
all.getRange('A4').values=[['Positive effect = nearer the top in rearranged species. BH correction is separate by assay. P-values assume independent species.']];
all.getRange('A5').values=[['Sensitivity columns use original ranks, relative ranks, or raw values. Saved statistical results; rerun analysis after input changes.']];
table(all,7,['Gene','Assay','Test status','Shared genes','Rearranged n','Unrearranged n','Rearranged median rank','Unrearranged median rank','Rank-biserial effect','Exact p','BH q','U (rearranged)','Label permutations','Original rank p','Original rank q','Relative rank p','Relative rank q','Raw value p','Raw value q'],tests.map(r=>fields.map(k=>k==='assay'?assay(r[k]):r[k])));
all.getRange(`G8:I${tests.length+7}`).setNumberFormat('0.000');
for(const c of ['J','K','N','O','P','Q','R','S'])all.getRange(`${c}8:${c}${tests.length+7}`).setNumberFormat('0.000000');
all.getRange(`L8:L${tests.length+7}`).setNumberFormat('0.0');
all.getRange(`M8:M${tests.length+7}`).setNumberFormat('#,##0');
all.getRange(`D8:S${tests.length+7}`).format.horizontalAlignment='right';
all.freezePanes.freezeRows(7);all.freezePanes.freezeColumns(2);
for(let i=0;i<tests.length;i++)if(tests[i].gene==='SIRT6')all.getRange(`A${i+8}:S${i+8}`).format.font.bold=true;
all.getRange(`K8:K${tests.length+7}`).conditionalFormats.add('cellIs',{operator:'lessThan',formula:.05,format:{fill:'#D8EADD',font:{bold:true}}});

const groupOrder={'Rearranged':0,'Unrearranged':1,'Unclassified':2,'Excluded':3};
const species=[...d.species].sort((a,b)=>groupOrder[a.group]-groupOrder[b.group]||a.common.localeCompare(b.common));
const focus=sheets.SIRT6;
base(focus,92,[30,19,17,16,16,16,17,16,16,16,18,18,18],'SIRT6: gene-order group and local rank');
focus.tabColor='#183F62';
focus.getRange('A3').values=[[`${d.species.length} species documented: ${d.group_counts.Rearranged} rearranged, ${d.group_counts.Unrearranged} unrearranged, ${d.group_counts.Excluded??0} excluded, ${d.group_counts.Unclassified??0} unclassified. Human, macaque and bowhead whale also excluded.`]];
focus.getRange('A4').values=[['Groups use marker order with explicit user assignments where recorded. These tests are not phylogenetically adjusted.']];
focus.getRange('A5').values=[[`Rank 1 is highest. The SIRT6 RNA comparison uses ${d.panels['rna|SIRT6'].length} shared genes because of missing measurements; ATAC uses ${d.panels['atac|SIRT6'].length}.`]];
table(focus,7,['Assay','Shared genes','Rearranged n','Unrearranged n','Rearranged median rank','Unrearranged median rank','Exact p','BH q','Rank-biserial effect','Original rank p','Original rank q','Raw value p','Raw value q'],[]);
for(let i=0;i<2;i++){
 const a=['rna','atac'][i],target=tests.findIndex(r=>r.gene==='SIRT6'&&r.assay===a)+8,row=8+i;
 const refs=['B','D','E','F','G','H','J','K','I','N','O','R','S'];
 focus.getRange(`A${row}:M${row}`).formulas=[refs.map(c=>`='All gene tests'!${c}${target}`)];
}
focus.getRange('B8:D9').setNumberFormat('0');focus.getRange('E8:F9').setNumberFormat('0.0');focus.getRange('G8:M9').setNumberFormat('0.000000');
focus.getRange('A11').values=[['RNA units: source liver log.tpm. ATAC units: integrated CPM·bp. Shared rank = primary comparison; original rank = previous workbook.']];
const lookup=Object.fromEntries(d.measurements.filter(r=>r.gene==='SIRT6').map(r=>[r.assay+'|'+r.species_id,r]));
const focusRows=species.map(s=>{const r=lookup['rna|'+s.id],a=lookup['atac|'+s.id];const status=(x)=>s.group==='Excluded'?'Excluded by user':x.used_in_test?'Included':x.value===null?'Unavailable':s.group==='Unclassified'?'Group unresolved':'Not tested';return [s.common,s.group,r.value,r.rank,r.original_rank,r.original_genes_measured,a.value,a.rank,a.original_rank,a.original_genes_measured,status(r),status(a)];});
table(focus,13,['Species','Group','RNA log.tpm','RNA shared rank','RNA original rank','RNA genes measured','ATAC CPM·bp','ATAC shared rank','ATAC original rank','ATAC genes measured','RNA test status','ATAC test status'],focusRows);
focus.getRange('C14:C44').setNumberFormat('0.000');focus.getRange('G14:G44').setNumberFormat('#,##0.0');
focus.getRange('C14:J44').format.horizontalAlignment='right';
for(const c of ['D','E','H','I'])focus.getRange(`${c}14:${c}44`).setNumberFormat('0.0');
focus.getRange('A47').values=[['Shared RNA panel: '+d.panels['rna|SIRT6'].join(', ')]];
focus.getRange('A49').values=[['Across the locus']];focus.getRange('A49').format.font.bold=true;
table(focus,50,['Assay','Tested genes','Nominal p < 0.05','BH q < 0.05'],[]);
for(let i=0;i<2;i++){
 const row=51+i,a=assay(['rna','atac'][i]),last=tests.length+7;
 focus.getRange(`A${row}`).values=[[a]];
 focus.getRange(`B${row}:D${row}`).formulas=[[
  `=COUNTIFS('All gene tests'!$B$8:$B$${last},A${row},'All gene tests'!$C$8:$C$${last},"Tested")`,
  `=COUNTIFS('All gene tests'!$B$8:$B$${last},A${row},'All gene tests'!$J$8:$J$${last},"<0.05",'All gene tests'!$C$8:$C$${last},"Tested")`,
  `=COUNTIFS('All gene tests'!$B$8:$B$${last},A${row},'All gene tests'!$K$8:$K$${last},"<0.05",'All gene tests'!$C$8:$C$${last},"Tested")`
 ]];
}
focus.getRange('A55').values=[['Analysis definitions and limitations']];focus.getRange('A55').format.font.bold=true;
let noteRow=57;
for(const note of d.notes){
 let line='';
 for(const word of note.split(' ')){if((line+' '+word).length>190){focus.getRange(`A${noteRow++}`).values=[[line]];line=word;}else line+=(line?' ':'')+word;}
 if(line)focus.getRange(`A${noteRow++}`).values=[[line]];
 noteRow++;
}
focus.getRange(`A57:M${noteRow}`).format={font:{name:'Arial',size:10,color:'#243746'},rowHeight:21};

const groups=sheets['Species groups'];
base(groups,44,[30,29,19,87,29,17,13,29,17,13,29,17,13,20,20,55,20,80],'Gene-order assignments used for comparison');
groups.getRange('A3').values=[['Sources: local browser transcript indexes and user assignments, 31 October 2026. TSS positions are 1-based inclusive.']];
groups.getRange('A4').values=[['User assignments override unresolved marker order. Original annotation groups and reasons are retained in the last two columns.']];
const grows=species.map(s=>[s.common,s.latin,s.group,s.reason,...['SIRT6','CREB3L3','TLE6'].flatMap(g=>{let l=s[g+'_loci'];return [l.map(x=>x.chrom).join('; ')||null,l.length===1?l[0].tss:null,l.map(x=>x.strand).join('; ')||null];}),s.CREB3L3_relative_bp,s.TLE6_relative_bp,s.evidence,s.annotation_group??s.group,s.annotation_reason??s.reason]);
table(groups,6,['Species','Scientific name','Group','Assignment basis','SIRT6 chromosome','SIRT6 TSS','SIRT6 strand','CREB3L3 chromosome','CREB3L3 TSS','CREB3L3 strand','TLE6 chromosome','TLE6 TSS','TLE6 strand','CREB3L3 from SIRT6 (bp)','TLE6 from SIRT6 (bp)','Assignment source','Annotation-only group','Annotation-only reason'],grows);
for(const c of ['F','I','L','N','O'])groups.getRange(`${c}7:${c}37`).setNumberFormat('#,##0');
groups.freezePanes.freezeRows(6);groups.freezePanes.freezeColumns(3);
groups.getRange('C7:C37').conditionalFormats.add('containsText',{text:'Unclassified',format:{fill:'#FFF0CF'}});
groups.getRange('C7:C37').conditionalFormats.add('containsText',{text:'Excluded',format:{fill:'#FFF0CF'}});

const raw=sheets.Measurements;
const mfields=['assay','gene','species','group','value','original_rank','original_genes_measured','rank','reference_genes','used_in_test'];
base(raw,d.measurements.length+10,[20,17,30,19,19,17,18,17,18,18],'Species-median measurements and ranks');
raw.getRange('A3').values=[['Source: 29 October 2026 locus matrices used for the 30 October ordered-gene workbook. Original measurements are preserved.']];
raw.getRange('A4').values=[['RNA is source liver log.tpm; ATAC is integrated CPM·bp. n.a. means unavailable. No zero imputation.']];
raw.getRange('A5').values=[['Original ranks include all measured candidate labels; primary ranks use a shared, target-specific reference panel.']];
table(raw,7,['Assay','Gene','Species','Gene-order group','Species median value','Original rank','Original genes measured','Shared-panel rank','Shared-panel genes','Used in primary test'],d.measurements.map(r=>mfields.map(k=>k==='assay'?assay(r[k]):k==='used_in_test'?(r[k]?'Yes':'No'):r[k])));
raw.getRange(`E8:E${d.measurements.length+7}`).setNumberFormat('0.0000');
raw.getRange(`E8:I${d.measurements.length+7}`).format.horizontalAlignment='right';
for(const c of ['F','H'])raw.getRange(`${c}8:${c}${d.measurements.length+7}`).setNumberFormat('0.0');
raw.freezePanes.freezeRows(7);raw.freezePanes.freezeColumns(3);

const refs=sheets['Rank panels'];
base(refs,tests.length+12,[20,17,16,115],'Shared genes defining each rank test');
refs.getRange('A3').values=[['For each target, these genes are measured in all classified species with a value for that target. All species in that test use this same list.']];
table(refs,5,['Assay','Target gene','Reference genes','Shared reference panel'],tests.map(r=>[assay(r.assay),r.gene,r.reference_genes,d.panels[r.assay+'|'+r.gene].join(', ')||'n.a.']));
refs.getRange(`D6:D${tests.length+5}`).format.wrapText=true;
for(let i=0;i<tests.length;i++)refs.getRange(`A${i+6}:D${i+6}`).format.rowHeight=Math.max(24,Math.ceil(d.panels[tests[i].assay+'|'+tests[i].gene].join(', ').length/104)*15+8);
refs.freezePanes.freezeRows(5);refs.freezePanes.freezeColumns(2);
refs.getRange(`A${tests.length+8}`).values=[['Test reference: '+d.references[0]]];
refs.getRange(`A${tests.length+9}`).values=[['Phylogenetic dependence: '+d.references[1]]];

wb.recalculate();
if(JSON.stringify(focus.getRange('B51:D52').values)!==JSON.stringify(['rna','atac'].map(a=>[d.summary[a].tested,d.summary[a].nominal_p05.length,d.summary[a].fdr05.length])))throw Error('Summary does not match independent analysis');
for(let i=0;i<2;i++){let r=d.results.find(r=>r.assay===['rna','atac'][i]&&r.gene==='SIRT6');if(Math.abs(focus.getRange(`G${8+i}`).values[0][0]-r.p)>1e-12)throw Error('SIRT6 p-value mismatch');}
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:8},maxChars:800})).ndjson);
console.log((await wb.inspect({kind:'table',range:'SIRT6!A7:I9',include:'values,formulas',tableMaxRows:3,tableMaxCols:9,maxChars:3000})).ndjson);
for(const [name,range] of [['SIRT6','A1:I19'],['All gene tests','A1:K18'],['Species groups','A1:F16'],['Measurements','A1:J16'],['Rank panels','A1:D12']]){
 const p=await wb.render({sheetName:name,range,scale:1.25,format:'png'});
 await fs.writeFile(path.join(out,'build',name.replaceAll(' ','_')+'.png'),new Uint8Array(await p.arrayBuffer()));
}
await(await SpreadsheetFile.exportXlsx(wb)).save(path.join(out,'SIRT6_rearrangement_gene_rank_tests.xlsx'));
console.log('Workbook saved; all source and test records retained; summary reconciled.');
