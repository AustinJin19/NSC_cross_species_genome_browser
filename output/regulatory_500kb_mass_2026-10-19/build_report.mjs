import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook,SpreadsheetFile} from '@oai/artifact-tool';
const here=path.dirname(fileURLToPath(import.meta.url));
const out=path.resolve(here,'../../outputs/regulatory-500kb-mass-2026-10-19');
const data=JSON.parse(await fs.readFile(path.join(here,'report_data.json'),'utf8'));
const go=JSON.parse(await fs.readFile(path.join(here,'go/report_data.json'),'utf8'));
const comparison=JSON.parse(await fs.readFile(path.join(here,'comparison.json'),'utf8'));
const wb=Workbook.create();
const summary=wb.worksheets.add('Summary'), genes=wb.worksheets.add('All genes'), terms=wb.worksheets.add('GO terms'), members=wb.worksheets.add('GO genes'), mapping=wb.worksheets.add('Gene mapping'), species=wb.worksheets.add('Species');
const navy='#123B60',ink='#203446';
const labels={FDR05_positive:'FDR positive',FDR05_negative:'FDR negative',nominal_p05_positive:'Nominal positive',nominal_p05_negative:'Nominal negative'};
function base(s,range){s.showGridLines=false;s.getRange(range).format.font={name:'Arial',size:10,color:ink};s.getRange(range).format.verticalAlignment='center';s.getRange(range).format.rowHeight=20;}
function title(s,text){s.getRange('A2').values=[[text]];s.getRange('A2').format.font={name:'Arial',size:15,color:navy,bold:true};s.getRange('A2').format.rowHeight=30;}
function header(s,range){s.getRange(range).format={fill:navy,font:{name:'Arial',size:10,color:'#FFFFFF',bold:true},wrapText:true,horizontalAlignment:'center',verticalAlignment:'center',rowHeight:34};s.getRange(range).format.borders={insideVertical:{style:'thin',color:'#FFFFFF'}};}
function widths(s,values,n){values.forEach((w,i)=>s.getRangeByIndexes(0,i,n,1).format.columnWidth=w);}
function table(s,records,headers,name,titleText,note,ws){
 const end=6+Math.max(1,records.length), last=String.fromCharCode(64+headers.length);base(s,`A1:${last}${end}`);title(s,titleText);s.getRange('A4').values=[[note]];s.getRange(`A6:${last}6`).values=[headers];
 if(records.length){s.getRange(`A7:${last}${end}`).values=records;s.tables.add(`A6:${last}${end}`,true,name);}else s.getRange('A7').values=[['No records']];
 header(s,`A6:${last}6`);widths(s,ws,end);s.freezePanes.freezeRows(6);s.freezePanes.freezeColumns(1);return end;
}
const geneRows=data.results.map(r=>[r.gene,r.n_species,r.partial_r,r.p_value,r.BH_FDR,r.MLS_coefficient,r.mass_coefficient_per_log10g,r.df,null,null]);
const ge=table(genes,geneRows,['Gene','Species (n)','Partial r','p-value','Gene BH FDR','MLS coefficient','Mass coefficient','df','FDR < 0.05','Nominal p < 0.05'],'GeneResults','Regulatory activity correlations','TSS ±500 kb peak count; log10 count and MLS; log10 body mass covariate. Positive/negative refer to the MLS association.',[22,14,14,17,17,18,18,9,18,23]);
genes.getRange('I7').formulas=[['=IF(E7<0.05,IF(C7>0,"Positive","Negative"),"No")']];genes.getRange(`I7:I${ge}`).fillDown();
genes.getRange('J7').formulas=[['=IF(D7<0.05,IF(C7>0,"Positive","Negative"),"No")']];genes.getRange(`J7:J${ge}`).fillDown();
genes.getRange(`B7:B${ge}`).setNumberFormat('0');genes.getRange(`C7:C${ge}`).setNumberFormat('0.000');genes.getRange(`D7:E${ge}`).setNumberFormat('0.00E+00');genes.getRange(`F7:G${ge}`).setNumberFormat('0.0000');
genes.getRange(`E7:E${ge}`).conditionalFormats.add('cellIs',{operator:'lessThan',formula:.05,format:{fill:'#E4F1EB'}});
const goRows=Object.values(go.groups).flatMap(g=>g.terms);
const te=table(terms,goRows.map(r=>[labels[r.list],r.go_id,r.term,r.go_fdr,null,r.overlap,r.query_size,r.term_in_background,r.background_size,r.raw_p]),['Input gene list','GO ID','Biological process','GO BH FDR','Fold enrichment','Overlap (k)','Mapped query (n)','Background term (K)','Background (N)','Raw p-value'],'GoResults','GO Biological Process enrichment','FDR lists use gene-level FDR < 0.05. Nominal lists are exploratory. GO-term BH FDR is corrected separately per input list.',[23,16,80,17,18,14,18,18,18,17]);
if(goRows.length){
 terms.getRange('E7').formulas=[['=(F7/G7)/(H7/I7)']];terms.getRange(`E7:E${te}`).fillDown();terms.getRange(`D7:D${te}`).setNumberFormat('0.00E+00');terms.getRange(`E7:E${te}`).setNumberFormat('0.00');terms.getRange(`F7:I${te}`).setNumberFormat('#,##0');terms.getRange(`J7:J${te}`).setNumberFormat('0.00E+00');terms.getRange(`C7:C${te}`).format.wrapText=true;terms.getRange(`A7:J${te}`).format.autofitRows();
 terms.getRange(`D7:D${te}`).conditionalFormats.add('cellIs',{operator:'lessThan',formula:.05,format:{fill:'#E4F1EB'}});
}
const me=table(members,go.memberships.map(r=>[labels[r.list],r.go_id,r.term,r.gene,r.ensembl_id]),['Input gene list','GO ID','Biological process','Input gene','Human Ensembl ID'],'GoMembers','Genes contributing to enriched terms','Only terms with GO-term FDR < 0.05. Filter by list and GO ID to see the contributing genes.',[23,16,88,24,27]);
if(go.memberships.length){members.getRange(`C7:C${me}`).format.wrapText=true;members.getRange(`A7:E${me}`).format.autofitRows();}
const mapEnd=table(mapping,go.mappings.map(r=>[labels[r.list],r.gene,r.status,r.human_ensembl_ids]),['Input gene list','Input gene','Mapping status','Human Ensembl ID / ambiguous candidates'],'Mapping','GO input gene mapping','Ambiguous and unmapped genes excluded. Duplicate human mappings are counted once.',[23,26,35,85]);
const se=table(species,data.species.map(s=>[s.species,s.scientific_name,s.MLS,s.body_mass_g,s.loci,s.selected_genes,s.excluded_region,s.excluded_zero]),['Species ID','Scientific name','MLS (years)','Body mass (g)','Indexed loci','Selected genes','Fail region rule','Zero median count'],'Species','Species and available regions','30 eligible species. Region and zero-count exclusions vary by gene and can overlap.',[31,30,15,19,17,19,19,20]);
species.getRange(`C7:D${se}`).setNumberFormat('#,##0.0');species.getRange(`E7:H${se}`).setNumberFormat('#,##0');

base(summary,'A1:F70');widths(summary,[65,17,18,19,18,18],70);title(summary,'Regulatory activity, lifespan and GO enrichment');summary.tabColor=navy;
summary.getRange('A4').values=[['Peak count within TSS ±500 kb (1 Mb total); analysis dated 19 October 2026.']];
summary.getRange('A6:C6').values=[['Gene association threshold','Positive','Negative']];header(summary,'A6:C6');
summary.getRange('A7:A8').values=[['Gene BH FDR < 0.05'],['Nominal p < 0.05']];
summary.getRange('B7:C8').formulas=[
 [`=COUNTIFS('All genes'!I7:I${ge},"Positive")`,`=COUNTIFS('All genes'!I7:I${ge},"Negative")`],
 [`=COUNTIFS('All genes'!J7:J${ge},"Positive")`,`=COUNTIFS('All genes'!J7:J${ge},"Negative")`]
];
summary.getRange('A10:A12').values=[['Genes tested'],['Maximum eligible species'],['Minimum gene BH FDR']];
summary.getRange('B10:B12').formulas=[[`=COUNTA('All genes'!A7:A${ge})`],[`=COUNTA('Species'!A7:A${se})`],[`=MIN('All genes'!E7:E${ge})`]];summary.getRange('B12').setNumberFormat('0.0000');
summary.getRange('A15:F15').values=[['GO input list','Input genes','Mapped genes','GO FDR < 0.05','Lowest GO FDR','Status']];header(summary,'A15:F15');
Object.entries(go.groups).forEach(([name,g],i)=>{
 const row=16+i, s=g.summary;
 summary.getRange(`A${row}:F${row}`).values=[[labels[name],s.input_genes,s.mapped_genes??null,s.run?null:null,s.minimum_go_fdr??null,s.run?'Completed':'Empty list']];
 if(s.run)summary.getRange(`D${row}`).formulas=[[`=COUNTIFS('GO terms'!A7:A${te},A${row},'GO terms'!D7:D${te},"<0.05")`]];
});summary.getRange('E16:E19').setNumberFormat('0.00E+00');
let leading=data.results.filter(r=>r.BH_FDR<.05);const strict=leading.length>0;
if(!strict)leading=data.results.filter(r=>r.p_value<.05);
leading=leading.slice(0,10);
summary.getRange('A22').values=[[strict?'Leading FDR-significant genes':'Leading nominal candidates (not FDR-significant)']];summary.getRange('A22').format.font.bold=true;
summary.getRange('A24:E24').values=[['Gene','Species','Partial r','p-value','Gene BH FDR']];header(summary,'A24:E24');
leading.forEach((r,i)=>{const source=data.results.indexOf(r)+7;summary.getRange(`A${25+i}:E${25+i}`).formulas=[[`='All genes'!A${source}`,`='All genes'!B${source}`,`='All genes'!C${source}`,`='All genes'!D${source}`,`='All genes'!E${source}`]];});
summary.getRange('C25:C34').setNumberFormat('0.000');summary.getRange('D25:E34').setNumberFormat('0.00E+00');
summary.getRange('A37').values=[['Methods and interpretation']];summary.getRange('A37').format.font.bold=true;
const notes=[
 'Response: log10(species-median peak count); predictor: log10(MLS in years).',
 'Covariate: log10(body mass in grams). Gene length is not included.',
 'Human, macaque, bowhead and species without a verified MLS are excluded.',
 'At least 50% of the requested genomic span must be available. Exactly 50% is retained.',
 'Highest-median locus is selected before the region filter; no lower-locus substitution.',
 'Zero counts excluded in log mode; no pseudocount. Missing peak files are not zero.',
 'At least 5 species per gene. Two-sided partial-correlation t test; df = n − 3.',
 'Gene-level BH FDR across all tested genes for this selected model.',
 'GO: human Biological Process annotations; custom background is all tested genes.',
 'GO-term BH FDR separately within each input list. Nominal-list enrichment is exploratory.',
 'GO terms overlap. Peak sharing can also make neighboring gene scores dependent.',
 'Associations do not establish causality. Phylogeny and sequencing depth are not adjusted.',
 'Model and GO statistics are static exports; rerun analysis to refresh the results.',
 'Source: local ATAC peak calls, indexed loci and supplied species traits with AnAge fallbacks.',
 'g:Profiler hsapiens GO:BP; version '+Object.values(go.groups).find(g=>g.summary.version)?.summary.version,
 'API: https://biit.cs.ut.ee/gprofiler/page/apis',
 'GO guidance: https://geneontology.org/docs/go-enrichment-analysis/',
 'The ±500 kb window spans 1 Mb. Peaks can be assigned to several neighboring genes.',
 'FDR is calculated for this window only, not jointly across all windows previously examined.',
 'Human gene-symbol annotations are a functional proxy, not verified one-to-one orthology.'
];notes.forEach((v,i)=>summary.getRange(`A${38+i}`).values=[[v]]);
summary.getRange('A60').values=[['Follow-up of the six FDR-significant genes from TSS ±50 kb']];
summary.getRange('A60').format.font.bold=true;
summary.getRange('A62:F62').values=[['Gene','Species (500 kb)','Partial r (500 kb)','p-value (500 kb)','FDR (500 kb)','Prior FDR (50 kb)']];header(summary,'A62:F62');
comparison.forEach((r,i)=>{
 const row=63+i,source=data.results.findIndex(g=>g.gene===r.gene)+7;
 summary.getRange(`A${row}`).values=[[r.gene]];
 if(source>=7)summary.getRange(`B${row}:E${row}`).formulas=[[`='All genes'!B${source}`,`='All genes'!C${source}`,`='All genes'!D${source}`,`='All genes'!E${source}`]];
 summary.getRange(`F${row}`).values=[[r.FDR_50kb]];
});
summary.getRange('C63:C68').setNumberFormat('0.000');summary.getRange('D63:F68').setNumberFormat('0.00E+00');
summary.getRange('A70').values=[['Prior FDR source: 19 October 2026 analysis of the same model at TSS ±50 kb.']];
wb.recalculate();
if(JSON.stringify(summary.getRange('B7:C8').values)!==JSON.stringify([[data.summary.FDR_positive,data.summary.FDR_negative],[data.summary.nominal_positive,data.summary.nominal_negative]]))throw Error('Gene summary mismatch');
for(const [i,g] of Object.values(go.groups).entries())if(g.summary.run&&summary.getRange(`D${16+i}`).values[0][0]!==g.summary.significant_terms)throw Error('GO summary mismatch');
console.log((await wb.inspect({kind:'table',range:'Summary!A6:F19',include:'values,formulas',tableMaxRows:14,tableMaxCols:6,maxChars:3000})).ndjson);
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},summary:'Formula errors',maxChars:1000})).ndjson);
for(const [sheet,range,file] of [['Summary','A1:F34','summary'],['Summary','A37:F57','methods'],['Summary','A60:F70','comparison'],['All genes','A1:J13','genes'],['GO terms','A1:J13','go'],['GO genes','A1:E13','members'],['Gene mapping','A1:D13','mapping'],['Species','A1:H13','species']]){
 if(process.argv.length>2&&!process.argv.slice(2).includes(file))continue;
 const image=await wb.render({sheetName:sheet,range,scale:1.3,format:'png'});await fs.writeFile(path.join(here,`${file}_preview.png`),new Uint8Array(await image.arrayBuffer()));
}
const book=await SpreadsheetFile.exportXlsx(wb);await book.save(path.join(out,'Regulatory_500kb_MLS_and_GO_2026-10-19.xlsx'));
console.log('Workbook saved.');
