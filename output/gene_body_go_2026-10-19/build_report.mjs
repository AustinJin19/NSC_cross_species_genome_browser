import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const here=path.dirname(fileURLToPath(import.meta.url));
const out=path.resolve(here,'../../outputs/go-gene-body-2026-10-19');
const data=JSON.parse(await fs.readFile(path.join(here,'report_data.json'),'utf8'));
const wb=Workbook.create();
const summary=wb.worksheets.add('Summary');
const pos=wb.worksheets.add('Positive GO');
const neg=wb.worksheets.add('Negative GO');
const members=wb.worksheets.add('Term genes');
const mapping=wb.worksheets.add('Gene mapping');
const navy='#123B60', ink='#203446';
function base(sheet,range){
 sheet.showGridLines=false;
 sheet.getRange(range).format.font={name:'Arial',size:10,color:ink};
 sheet.getRange(range).format.verticalAlignment='center';
 sheet.getRange(range).format.rowHeight=20;
}
function heading(sheet,cell,title){
 sheet.getRange(cell).values=[[title]];
 sheet.getRange(cell).format.font={name:'Arial',size:15,bold:true,color:navy};
 sheet.getRange(cell).format.rowHeight=29;
}
function header(sheet,range){
 sheet.getRange(range).format={fill:navy,font:{name:'Arial',size:10,color:'#FFFFFF',bold:true},wrapText:true,horizontalAlignment:'center',verticalAlignment:'center',rowHeight:32};
 sheet.getRange(range).format.borders={insideVertical:{style:'thin',color:'#FFFFFF'}};
}
function widths(sheet,ws,n){ws.forEach((width,i)=>sheet.getRangeByIndexes(0,i,n,1).format.columnWidth=width);}
function terms(sheet,group){
 const rows=data.groups[group].terms, end=6+rows.length;
 base(sheet,`A1:J${end}`);
 heading(sheet,'A2',`${group==='positive'?'Positive':'Negative'} correlation: GO Biological Process`);
 sheet.getRange('A4').values=[[group==='positive'?'No terms pass GO-term BH FDR < 0.05.':'212 terms pass GO-term BH FDR < 0.05. Contributing genes are listed in Term genes.']];
 sheet.getRange('A6:J6').values=[['GO ID','Biological process','BH FDR','Fold enrichment','Overlap (k)','Mapped query (n)','Background term (K)','Background (N)','Raw p-value','FDR < 0.05']];
 sheet.getRange(`A7:J${end}`).values=rows.map(r=>[r.go_id,r.term,r.fdr,null,r.overlap,r.query_size,r.term_in_background,r.background_size,r.raw_p,null]);
 sheet.getRange('D7').formulas=[['=(E7/F7)/(G7/H7)']];sheet.getRange(`D7:D${end}`).fillDown();
 sheet.getRange('J7').formulas=[['=IF(C7<0.05,"Yes","No")']];sheet.getRange(`J7:J${end}`).fillDown();
 sheet.tables.add(`A6:J${end}`,true,`${group}GoTerms`);header(sheet,'A6:J6');
 widths(sheet,[16,78,16,17,13,17,18,18,16,15],end);
 sheet.getRange(`B7:B${end}`).format.wrapText=true;
 sheet.getRange(`A7:J${end}`).format.autofitRows();
 sheet.getRange(`C7:C${end}`).setNumberFormat('0.00E+00');
 sheet.getRange(`D7:D${end}`).setNumberFormat('0.00');
 sheet.getRange(`E7:H${end}`).setNumberFormat('#,##0');
 sheet.getRange(`I7:I${end}`).setNumberFormat('0.00E+00');
 sheet.getRange(`C7:I${end}`).format.horizontalAlignment='right';
 sheet.getRange(`C7:C${end}`).conditionalFormats.add('cellIs',{operator:'lessThan',formula:0.05,format:{fill:'#E4F1EB'}});
 sheet.freezePanes.freezeRows(6);sheet.freezePanes.freezeColumns(2);
 return end;
}
const posEnd=terms(pos,'positive'), negEnd=terms(neg,'negative');
const termRows=[];
for(const r of data.groups.negative.terms.filter(r=>r.significant)){
 const genes=r.genes.split(';'), ids=r.ensembl_ids.split(';');
 if(genes.length!==r.overlap||ids.length!==r.overlap)throw Error('Term membership count mismatch');
 genes.forEach((g,i)=>termRows.push([r.go_id,r.term,g,ids[i]]));
}
const te=6+termRows.length;
base(members,`A1:D${te}`);heading(members,'A2','Genes contributing to enriched terms');
members.getRange('A4').values=[['Negative list; one row per human gene and significant GO term. Genes can belong to multiple terms.']];
members.getRange('A6:D6').values=[['GO ID','Biological process','Input gene symbol','Human Ensembl gene']];
members.getRange(`A7:D${te}`).values=termRows;
members.tables.add(`A6:D${te}`,true,'TermMembership');header(members,'A6:D6');
widths(members,[16,94,27,27],te);
members.getRange(`B7:B${te}`).format.wrapText=true;
members.getRange(`A7:D${te}`).format.autofitRows();
members.freezePanes.freezeRows(6);members.freezePanes.freezeColumns(1);
const me=6+data.mapping.length;
base(mapping,`A1:D${me}`);heading(mapping,'A2','Input gene mapping');
mapping.getRange('A4').values=[['Ambiguous IDs are candidates only; they were excluded. Duplicate human gene mappings were counted once.']];
mapping.getRange('A6:D6').values=[['Correlation list','Input symbol','Mapping status','Human Ensembl ID / ambiguous candidates']];
mapping.getRange(`A7:D${me}`).values=data.mapping.map(r=>[r.group,r.symbol,r.status,r.human_ensembl_ids]);
mapping.tables.add(`A6:D${me}`,true,'InputGeneMapping');header(mapping,'A6:D6');
widths(mapping,[21,25,48,75],me);mapping.freezePanes.freezeRows(6);mapping.freezePanes.freezeColumns(2);
summary.tabColor=navy;pos.tabColor='#007F68';neg.tabColor='#A14E39';

base(summary,'A1:D45');widths(summary,[77,17,18,18],45);
heading(summary,'A2','GO enrichment of lifespan-associated genes');
summary.getRange('A4').values=[['Gene body + 2 kb promoter; nominal association p < 0.05; analysis dated 19 October 2026.']];
summary.getRange('A6:C6').values=[['List metrics','Positive','Negative']];header(summary,'A6:C6');
summary.getRange('A7:A15').values=[['Input symbols'],['Mapped human genes'],['Unmapped symbols'],['Ambiguous symbols'],['Duplicate mapping'],['GO terms returned'],['GO terms at BH FDR < 0.05'],['Lowest GO-term BH FDR'],['Effective mapped background']];
const groups=`'Gene mapping'!$A$7:$A$${me}`, status=`'Gene mapping'!$C$7:$C$${me}`;
for(const [group,col,sheet,end] of [['positive','B','Positive GO',posEnd],['negative','C','Negative GO',negEnd]]){
 const rr=`'${sheet}'!$C$7:$C$${end}`;
 summary.getRange(`${col}7:${col}15`).formulas=[
  [`=COUNTIFS(${groups},"${group}")`],
  [`=COUNTIFS(${groups},"${group}",${status},"Mapped")`],
  [`=COUNTIFS(${groups},"${group}",${status},"Unmapped; excluded")`],
  [`=COUNTIFS(${groups},"${group}",${status},"Ambiguous; excluded")`],
  [`=COUNTIFS(${groups},"${group}",${status},"Duplicate mapping; not counted again")`],
  [`=COUNTA('${sheet}'!A7:A${end})`],
  [`=COUNTIFS(${rr},"<0.05")`],
  [`=MIN(${rr})`],
  [`='${sheet}'!H7`]
 ];
}
summary.getRange('B7:C15').setNumberFormat('#,##0');summary.getRange('B14:C14').setNumberFormat('0.00E+00');
summary.getRange('A13:C13').format.fill='#EAF1F8';summary.getRange('A13:C13').format.font.bold=true;
heading(summary,'A18','Selected enriched processes in the negative list');
summary.getRange('A20:D20').values=[['Biological process','Genes','Fold enrichment','BH FDR']];header(summary,'A20:D20');
data.representative_ids.forEach((gid,i)=>{
 const idx=data.groups.negative.terms.findIndex(r=>r.go_id===gid)+7;
 summary.getRange(`A${21+i}:D${21+i}`).formulas=[[`='Negative GO'!B${idx}`,`='Negative GO'!E${idx}`,`='Negative GO'!D${idx}`,`='Negative GO'!C${idx}`]];
});
summary.getRange('B21:B28').setNumberFormat('#,##0');summary.getRange('C21:C28').setNumberFormat('0.00');summary.getRange('D21:D28').setNumberFormat('0.00E+00');
summary.getRange('A30').values=[['Interpretation and methods']];summary.getRange('A30').format.font.bold=true;
const notes=[
 'Exploratory: no individual input gene passes association-level FDR < 0.05.',
 'The 212 enriched terms overlap; they are not 212 independent pathways.',
 'Negative means lower adjusted ATAC activity with longer MLS, not pathway inhibition.',
 'Model: log10(activity) versus log10(MLS), with log10(length) and log10(mass) covariates.',
 'Current species exclusions and the 50% region-availability rule were retained.',
 'GO test: hypergeometric overrepresentation; BH FDR separately within each list.',
 'Background: all 20,200 tested symbols; 19,350 mapped human genes.',
 'Human symbols are a functional proxy; cross-species orthology has not been verified.',
 'Statistics are static exports. Rerun the analysis to refresh GO tests and mappings.',
 'Source: g:Profiler, hsapiens, GO:BP; version e114_eg62_p19_9dc92ac8.',
 'https://biit.cs.ut.ee/gprofiler/page/apis',
 'GO background guidance: https://geneontology.org/docs/go-enrichment-analysis/',
 'Full requests, responses and numerical verification were retained with the analysis.'
];
notes.forEach((note,i)=>summary.getRange(`A${31+i}`).values=[[note]]);
summary.getRange('A40:A43').format.font.color='#526477';
wb.recalculate();
const metrics=summary.getRange('B7:C15').values;
if(JSON.stringify(metrics.slice(0,7))!==JSON.stringify([[836,1936],[801,1894],[26,28],[9,13],[0,1],[5986,8326],[0,212]]))throw Error('Summary mismatch '+JSON.stringify(metrics));
for(const [sheet,rows] of [[pos,data.groups.positive.terms],[neg,data.groups.negative.terms]]){
 for(const i of [0,Math.floor(rows.length/2),rows.length-1]){
  if(Math.abs(sheet.getRange(`D${7+i}`).values[0][0]-rows[i].fold_enrichment)>1e-12)throw Error('Fold enrichment mismatch');
 }
}
console.log((await wb.inspect({kind:'table',range:'Summary!A6:C15',include:'values,formulas',tableMaxRows:10,tableMaxCols:3,maxChars:2500})).ndjson);
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},summary:'Final formula error scan',maxChars:1000})).ndjson);
for(const [sheet,range,name] of [['Summary','A1:D43','summary'],['Positive GO','A1:J12','positive'],['Negative GO','A1:J12','negative'],['Term genes','A1:D12','members'],['Gene mapping','A1:D12','mapping']]){
 const image=await wb.render({sheetName:sheet,range,scale:1.4,format:'png'});
 await fs.writeFile(path.join(here,`${name}_preview.png`),new Uint8Array(await image.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(path.join(out,'GO_gene_body_positive_negative_2026-10-19.xlsx'));
console.log(JSON.stringify({file:path.join(out,'GO_gene_body_positive_negative_2026-10-19.xlsx'),positive_terms:5986,negative_terms:8326,significant_negative:212,gene_term_memberships:termRows.length}));
