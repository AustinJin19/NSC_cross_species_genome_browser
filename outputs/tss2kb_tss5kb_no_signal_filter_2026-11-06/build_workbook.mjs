import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const dir=path.dirname(fileURLToPath(import.meta.url));
const data=JSON.parse(await fs.readFile(path.join(dir,'results.json'),'utf8'));
const wb=Workbook.create();
const blue='#173E54', pale='#EDF4F7', ink='#203B4A', green='#E5F3EC', red='#FBEAEA';
const summary=wb.worksheets.add('Summary');
const listInfo=[];

function base(name,title,note,headers,rows,widths){
  const sh=wb.worksheets.add(name), cols=headers.length, last=String.fromCharCode(64+cols);
  sh.showGridLines=false;
  sh.getRange(`A1:${last}${Math.max(8,rows.length+5)}`).format.font={name:'Arial',size:10,color:ink};
  sh.getRange('A1').values=[[title]];
  sh.getRange('A1').format.font={name:'Arial',size:16,bold:true,color:blue};
  sh.getRange('A2').values=[[note]];
  sh.getRange('A3').values=[['← Summary: counts, analysis settings and interpretation']];
  sh.getRange(`A5:${last}5`).values=[headers];
  sh.getRange(`A5:${last}5`).format={fill:blue,font:{bold:true,color:'#FFFFFF'},wrapText:true,rowHeight:32};
  if(rows.length){
    sh.getRange(`A6:${last}${rows.length+5}`).values=rows;
    const table=sh.tables.add(`A5:${last}${rows.length+5}`,true,name.replace(/[^A-Za-z]/g,'')+'Table');
    table.showFilterButton=true;
  } else sh.getRange('A6').values=[['No genes meet this threshold.']];
  widths.forEach((w,i)=>sh.getRange(`${String.fromCharCode(65+i)}:${String.fromCharCode(65+i)}`).format.columnWidth=w);
  sh.getRange('1:1').format.rowHeight=27;
  sh.freezePanes.freezeRows(5);
  sh.freezePanes.freezeColumns(1);
  return sh;
}

for(const mass of ['No','Yes']){
  for(const nominal of [false,true]){
    const name=nominal?(mass==='Yes'?'Nominal with mass':'Nominal without mass'):(mass==='Yes'?'With body mass':'Without body mass');
    const rows=data.results.filter(r=>r.body_mass_correction===mass&&(nominal?r.p<.05&&r.BH_FDR>=.05:r.BH_FDR<.05))
      .sort((a,b)=>a.window_kb-b.window_kb||a.p-b.p||a.gene.localeCompare(b.gene));
    const sh=base(name,
      `${nominal?'Nominal-only hits':'FDR-significant genes'} · ${mass==='Yes'?'with':'without'} body-mass correction`,
      nominal?'P < 0.05 but BH FDR ≥ 0.05. Exploratory only; FDR-significant genes are on the primary list.':
      'BH FDR < 0.05. Each row is one gene–window result; a gene may appear at both windows.',
      ['Gene','TSS ± (kb)','Direction','Correlation r','P value','BH FDR','Species n','MLS coefficient','Log-mass coefficient'],
      rows.map(r=>[r.gene,r.window_kb,r.direction,r.r,r.p,r.BH_FDR,r.n_species,r.MLS_beta,r.mass_beta]),
      [21,13,13,17,16,16,13,20,23]);
    const end=Math.max(6,rows.length+5);
    sh.getRange(`D6:D${end}`).setNumberFormat('0.000');
    sh.getRange(`E6:F${end}`).setNumberFormat('0.00E+00');
    sh.getRange(`H6:I${end}`).setNumberFormat('0.000');
    sh.getRange(`B6:B${end}`).setNumberFormat('0');
    sh.getRange(`G6:G${end}`).setNumberFormat('0');
    if(rows.length){
      sh.getRange(`C6:C${end}`).conditionalFormats.add('containsText',{text:'Positive',format:{fill:green}});
      sh.getRange(`C6:C${end}`).conditionalFormats.add('containsText',{text:'Negative',format:{fill:red}});
    }
    listInfo.push({name,mass,nominal,rows:rows.length,end,genes:new Set(rows.map(r=>r.gene)).size});
  }
}

const species=base('Species','Species and traits',
  '31 species have usable MLS values before gene-specific coverage filtering. Common names follow the browser catalog.',
  ['Species ID','Common name','Scientific name','MLS (years)','Body mass (g)','Analysis status'],
  data.inventory.map(r=>[r.species,r.common_name,r.scientific_name,r.MLS_years,r.body_mass_g,r.status]),[30,29,30,16,19,31]);
species.getRange('D6:E39').setNumberFormat('#,##0.00');

summary.showGridLines=false;
summary.getRange('A1:I37').format.font={name:'Arial',size:10,color:ink};
summary.getRange('A1').values=[['TSS gene activity and maximum lifespan']];
summary.getRange('A1').format.font={size:18,bold:true,color:blue};
summary.getRange('1:1').format.rowHeight=30;
summary.getRange('A2').values=[['TSS ±2 kb and ±5 kb · positive-signal coverage filter removed · 6 November 2026']];
summary.getRange('A3').values=[['Primary significance: BH FDR < 0.05. Nominal P < 0.05 counts include both FDR hits and nominal-only hits.']];
summary.getRange('A6:I6').values=[['TSS ± (kb)','Body-mass correction','Genes tested','FDR positive','FDR negative','FDR total','Nominal positive','Nominal negative','Nominal total']];
summary.getRange('A6:I6').format={fill:blue,font:{bold:true,color:'#FFFFFF'},wrapText:true,rowHeight:36};
for(let i=0;i<data.summary.length;i++){
  const r=data.summary[i],rr=7+i;
  const fdr=listInfo.find(s=>s.mass===r.body_mass_correction&&!s.nominal),nom=listInfo.find(s=>s.mass===r.body_mass_correction&&s.nominal);
  const count=(s,d)=>`COUNTIFS('${s.name}'!$B$6:$B$${s.end},A${rr},'${s.name}'!$C$6:$C$${s.end},"${d}")`;
  summary.getRange(`A${rr}:C${rr}`).values=[[r.window_kb,r.body_mass_correction==='Yes'?'With log10(mass)':'Without',r.genes_tested]];
  summary.getRange(`D${rr}:I${rr}`).formulas=[[`=${count(fdr,'Positive')}`,`=${count(fdr,'Negative')}`,`=SUM(D${rr}:E${rr})`,
    `=D${rr}+${count(nom,'Positive')}`,`=E${rr}+${count(nom,'Negative')}`,`=SUM(G${rr}:H${rr})`]];
}
summary.getRange('C7:I10').setNumberFormat('#,##0');
summary.getRange('D7:F10').format.fill=pale;
summary.getRange('A12').values=[['Two primary gene lists']];
summary.getRange('A12').format.font={bold:true,size:12,color:blue};
summary.getRange('A13:C15').values=[['Worksheet','Unique genes','Gene–window hits'],
  ...listInfo.filter(x=>!x.nominal).map(s=>[s.name,s.genes,s.rows])];
summary.getRange('A13:C13').format={fill:pale,font:{bold:true,color:blue}};
const notes=[
  ['Model','Ordinary regression: log10(activity) ~ log10(MLS), with or without log10(body mass in g). No phylogenetic correction.'],
  ['Activity','Per-sample sum of CPM × bp across the full TSS window, followed by the species median; log10 is applied afterward.'],
  ['Direction','Positive / negative refers to the MLS coefficient. With body mass, r is the partial correlation controlling for log10(mass).'],
  ['Species','Human and bowhead whale excluded; macaque included. Cave nectar bat excluded because MLS is reported only as a lower bound.'],
  ['Coverage','Retain a species if at least 50% of the full requested window is inside chromosome/scaffold bounds. Exactly 50% passes.'],
  ['Signal filter','No positive-base coverage threshold is applied. Species with nonpositive median activity remain excluded because log10 is undefined.'],
  ['Multiple loci','Select the highest species-median activity locus before the genomic-span filter; do not substitute another locus after exclusion.'],
  ['Sample policy','Current browser sample catalog, including low-TSS samples; no additional sample-QC filtering was applied.'],
  ['Testable genes','At least 5 eligible species and nondegenerate regression. No gene-length correction is used for these fixed TSS windows.'],
  ['Multiple testing','Benjamini–Hochberg correction is performed separately for each of the four window/model combinations.'],
  ['Gene matching','Current browser gene-name mapping (HGNC aliases, supported HCOP orthologs and projected annotations); ambiguous matches are not pooled.'],
  ['Interpretation','Associations are exploratory. Phylogeny, technical differences and overlapping gene windows can affect these results.'],
  ['Source','Current project bigWigs, gene-activity annotation index and browser species MLS/body-mass table; full source fingerprints accompany the workbook.'],
  ['Reproducibility','Source-validated measurements reused; correlations and BH FDR recomputed. Full tested-gene tables and validation accompany this workbook.']
];
summary.getRange('A17').values=[['Methods and reading the results']];
summary.getRange('A17').format.font={bold:true,size:12,color:blue};
notes.forEach(([label,value],i)=>{
  const r=18+i;summary.getRange(`A${r}`).values=[[label]];summary.getRange(`A${r}`).format.font={bold:true,color:blue};
  summary.getRange(`B${r}:I${r}`).merge();summary.getRange(`B${r}`).values=[[value]];
  summary.getRange(`A${r}:I${r}`).format={rowHeight:32,wrapText:true,verticalAlignment:'center'};
});
[21,25,17,17,17,16,21,21,19].forEach((w,i)=>summary.getRange(`${String.fromCharCode(65+i)}:${String.fromCharCode(65+i)}`).format.columnWidth=w);
summary.freezePanes.freezeRows(6);

wb.recalculate();
const calculated=summary.getRange('D7:I10').values;
data.summary.forEach((r,i)=>{
  const expected=[r.FDR_positive,r.FDR_negative,r.FDR_positive+r.FDR_negative,r.nominal_positive,r.nominal_negative,r.nominal_positive+r.nominal_negative];
  if(JSON.stringify(calculated[i])!==JSON.stringify(expected))throw Error(`Count mismatch ${i}: ${JSON.stringify(calculated[i])} vs ${JSON.stringify(expected)}`);
});
await fs.writeFile(path.join(dir,'summary_inspection.json'),JSON.stringify(await wb.inspect({kind:'table',range:'Summary!A6:I15',include:'values,formulas',tableMaxRows:10,tableMaxCols:9}),null,2));
await fs.writeFile(path.join(dir,'workbook_inspection.json'),JSON.stringify(await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#N/A',options:{useRegex:true,maxResults:100}}),null,2));
await fs.mkdir(path.join(dir,'previews'),{recursive:true});
for(const name of ['Summary',...listInfo.map(s=>s.name),'Species']){
  const range=name==='Summary'?'A1:I31':name==='Species'?'A1:F16':'A1:I16';
  const img=await wb.render({sheetName:name,range,scale:1.3,format:'png'});
  await fs.writeFile(path.join(dir,'previews',name.replaceAll(' ','_')+'.png'),new Uint8Array(await img.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);
const dest=path.join(dir,'TSS_2kb_5kb_MLS_no_positive_signal_filter.xlsx');
await xlsx.save(dest);
console.log(JSON.stringify({output:dest,calculated,lists:listInfo},null,2));
