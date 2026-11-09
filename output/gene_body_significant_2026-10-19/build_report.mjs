import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const here=path.dirname(fileURLToPath(import.meta.url));
const root=path.resolve(here,'../..');
const output=path.join(root,'outputs/gene-body-mls-2026-10-19');
const data=JSON.parse(await fs.readFile(path.join(here,'report_data.json'),'utf8'));
const rows=data.rows;
const positive=rows.filter(r=>r.p_value<.05&&r.partial_r>0);
const negative=rows.filter(r=>r.p_value<.05&&r.partial_r<0);
const wb=Workbook.create();
const summary=wb.worksheets.add('Summary');
const pos=wb.worksheets.add('Nominal positive');
const neg=wb.worksheets.add('Nominal negative');
const all=wb.worksheets.add('All tested');
const species=wb.worksheets.add('Species');
const navy='#003F7D', pale='#EAF1F8', ink='#172B3A';

function base(sheet,range){
 sheet.showGridLines=false;
 sheet.getRange(range).format.font={name:'Arial',size:10,color:ink};
 sheet.getRange(range).format.verticalAlignment='center';
 sheet.getRange(range).format.rowHeight=19;
}
function heading(sheet,cell,title){
 sheet.getRange(cell).values=[[title]];
 sheet.getRange(cell).format.font={name:'Arial',size:15,bold:true,color:navy};
 sheet.getRange(cell).format.rowHeight=27;
}
function header(sheet,range){
 sheet.getRange(range).format={fill:navy,font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},horizontalAlignment:'center',verticalAlignment:'center',wrapText:true,rowHeight:30};
 sheet.getRange(range).format.borders={insideVertical:{style:'thin',color:'#FFFFFF'}};
}
const fields=['gene','n_species','partial_r','p_value','BH_FDR','MLS_coefficient','length_coefficient_per_log10kb','mass_coefficient_per_log10g','df'];
const labels=['Gene','Species (n)','Partial r','p-value','BH FDR','MLS coefficient','Length coefficient','Mass coefficient','df'];
const widths=[20,12,13,15,15,16,18,18,9];
function geneTable(sheet,records,title,tableName){
 const end=6+records.length;
 base(sheet,`A1:I${end}`);
 heading(sheet,'A2',title);
 sheet.getRange('A4').values=[['Nominal p < 0.05 is exploratory. No gene passes BH FDR < 0.05.']];
 sheet.getRange('A4').format.font={italic:true,color:'#526477'};
 sheet.getRange('A6:I6').values=[labels];
 sheet.getRange(`A7:I${end}`).values=records.map(r=>fields.map(k=>r[k]));
 sheet.tables.add(`A6:I${end}`,true,tableName);
 header(sheet,'A6:I6');
 for(let i=0;i<widths.length;i++)sheet.getRangeByIndexes(0,i,end,1).format.columnWidth=widths[i];
 sheet.getRange(`B7:I${end}`).format.horizontalAlignment='right';
 sheet.getRange(`B7:B${end}`).setNumberFormat('0');
 sheet.getRange(`C7:C${end}`).setNumberFormat('0.000');
 sheet.getRange(`D7:E${end}`).setNumberFormat('0.000E+00');
 sheet.getRange(`F7:H${end}`).setNumberFormat('0.0000');
 sheet.getRange(`I7:I${end}`).setNumberFormat('0');
 sheet.freezePanes.freezeRows(6);
 sheet.freezePanes.freezeColumns(1);
 return end;
}
geneTable(pos,positive,'Positive nominal associations','PositiveGenes');
geneTable(neg,negative,'Negative nominal associations','NegativeGenes');
const last=geneTable(all,rows,'All tested genes','AllTestedGenes');
all.getRange('A4').values=[['Imported model estimates; one test per gene. Covariates: log10(length) and log10(body mass).']];
pos.tabColor='#007A65';neg.tabColor='#9E4B34';summary.tabColor=navy;

const speciesRows=data.species.map(s=>[s.species,s.scientific_name,s.MLS,s.MLS_source.body_mass_g,s.loci,s.selected_genes_excluded_by_region,s.MLS_source.source||'',s.MLS_source.note||'']);
const se=6+speciesRows.length;
base(species,`A1:H${se}`);heading(species,'A2','Species and trait provenance');
species.getRange('A4').values=[['Eligible species before gene-specific availability, positive-activity and annotation filters.']];
species.getRange('A6:H6').values=[['Species ID','Scientific name','MLS (years)','Body mass (g)','Scored loci','Genes failing region rule','Trait source','Source note']];
species.getRange(`A7:H${se}`).values=speciesRows;
species.tables.add(`A6:H${se}`,true,'SpeciesTraits');header(species,'A6:H6');
[32,29,14,18,15,24,65,65].forEach((width,i)=>species.getRangeByIndexes(0,i,se,1).format.columnWidth=width);
species.getRange(`C7:D${se}`).setNumberFormat('#,##0.0');
species.getRange(`E7:F${se}`).setNumberFormat('#,##0');
species.getRange(`G7:H${se}`).format.wrapText=true;
species.getRange(`A7:H${se}`).format.autofitRows();
species.freezePanes.freezeRows(6);species.freezePanes.freezeColumns(2);

base(summary,'A1:F26');
for(const [column,width] of [['A',33],['B',13],['C',13],['D',4],['E',25],['F',52]])summary.getRange(`${column}1:${column}26`).format.columnWidth=width;
heading(summary,'A2','Gene activity & lifespan');
summary.getRange('A4').values=[['Gene body + 2 kb upstream; analysis dated 19 Oct 2026']];
summary.getRange('A6:C8').values=[['Threshold','Positive','Negative'],['BH FDR < 0.05',null,null],['Nominal p < 0.05',null,null]];
header(summary,'A6:C6');
const rr=`'All tested'!$C$7:$C$${last}`, pp=`'All tested'!$D$7:$D$${last}`, qq=`'All tested'!$E$7:$E$${last}`;
summary.getRange('B7:C8').formulas=[
 [`=COUNTIFS(${rr},">0",${qq},"<0.05")`,`=COUNTIFS(${rr},"<0",${qq},"<0.05")`],
 [`=COUNTIFS(${rr},">0",${pp},"<0.05")`,`=COUNTIFS(${rr},"<0",${pp},"<0.05")`]
];
summary.getRange('A11:A14').values=[['Genes tested'],['Eligible species'],['Minimum species per gene'],['Minimum BH FDR']];
summary.getRange('B11:B14').formulas=[
 [`=COUNTA('All tested'!A7:A${last})`],[`=COUNTA('Species'!A7:A${se})`],
 [`=MIN('All tested'!B7:B${last})`],[`=MIN(${qq})`]
];
summary.getRange('B7:C13').setNumberFormat('#,##0');summary.getRange('B14').setNumberFormat('0.0000');
summary.getRange('A7:C7').format.fill=pale;
summary.getRange('A7:C7').format.font.bold=true;
summary.getRange('E6:F14').values=[
 ['Model','Current browser settings'],['Response','log10(median activity)'],['Predictor','log10(MLS in years)'],
 ['Covariates','log10(gene length) + log10(body mass)'],['Region completeness','At least 50% of the requested genomic span'],
 ['Locus selection','Highest median, then availability filter'],['Excluded species','Human, macaque, bowhead whale'],
 ['Test','Two-sided partial-correlation t test'],['Multiple testing','BH across all 20,200 log10-MLS tests']
];
header(summary,'E6:F6');
summary.getRange('A18').values=[['Interpretation']];summary.getRange('A18').format.font.bold=true;
summary.getRange('A19').values=[['The nominal lists do not meet FDR < 0.05. Positive and negative refer to the adjusted association with log10(MLS).']];
summary.getRange('A20').values=[['Gene symbols are not verified one-to-one orthologs. Phylogeny is not controlled.']];
summary.getRange('A22').values=[['Data and refresh']];summary.getRange('A22').format.font.bold=true;
summary.getRange('A23').values=[['Source: local ATAC bigWigs and gene annotations; supplied species traits with AnAge fallbacks (see Species).']];
summary.getRange('A24').values=[['Model estimates are a static analysis export. Rerun the analysis to regenerate the gene lists and estimates.']];
summary.getRange('A26').values=[['Source tables: log10_MLS/all_tested.tsv and species_provenance.json, generated 19 Oct 2026.']];
summary.getRange('A26').format.font={size:10,color:'#526477',italic:true};

wb.recalculate();
const counts=summary.getRange('B7:C8').values;
const metrics=summary.getRange('B11:B14').values;
if(JSON.stringify(counts)!==JSON.stringify([[0,0],[836,1936]]))throw Error('Summary counts disagree: '+JSON.stringify(counts));
if(metrics[0][0]!==20200||metrics[1][0]!==30||Math.abs(metrics[3][0]-data.summary.minimum_BH_FDR)>1e-12)throw Error('Summary verification failed');
console.log((await wb.inspect({kind:'table',range:'Summary!A6:F14',include:'values,formulas',tableMaxRows:9,tableMaxCols:6,maxChars:4000})).ndjson);
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30},summary:'Formula error scan',maxChars:2000})).ndjson);
await fs.mkdir(output,{recursive:true});
for(const [sheet,range,file] of [['Summary','A1:F26','summary'],['Nominal positive','A1:I15','positive'],['Nominal negative','A1:I15','negative'],['All tested','A1:I15','all_tested'],['Species','A1:F15','species'],['Species','G6:H12','species_sources']]){
 if(process.argv.length>2&&!process.argv.slice(2).includes(sheet))continue;
 const preview=await wb.render({sheetName:sheet,range,scale:1.5,format:'png'});
 await fs.writeFile(path.join(here,file+'_preview.png'),new Uint8Array(await preview.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(path.join(output,'Gene_body_MLS_correlations_2026-10-19.xlsx'));
console.log(JSON.stringify({workbook:path.join(output,'Gene_body_MLS_correlations_2026-10-19.xlsx'),tested:rows.length,nominalPositive:positive.length,nominalNegative:negative.length,FDR_significant:0}));
