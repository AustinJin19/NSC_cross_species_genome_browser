import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const dir=path.dirname(fileURLToPath(import.meta.url));
const d=JSON.parse(await fs.readFile(path.join(dir,'report_data.json'),'utf8'));
const wb=Workbook.create();
const overview=wb.worksheets.add('Comparison');
const sig=wb.worksheets.add('Significant terms');
const plot=wb.worksheets.add('Plot terms');
const inputs=wb.worksheets.add('Input genes');
const methods=wb.worksheets.add('Methods');
const font={name:'Arial',size:11,color:'#23394B'};
const navy='#173B55', pale='#F0F5F8';
const col=n=>{let s='';for(n++;n>0;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
function table(sheet,rows,widths,start,name){
  const keys=Object.keys(rows[0]),last=col(keys.length-1),end=start+rows.length;
  sheet.getRange(`A${start}:${last}${end}`).values=[keys,...rows.map(r=>keys.map(k=>r[k]??null))];
  sheet.getRange(`A${start}:${last}${end}`).format={font,wrapText:true,verticalAlignment:'center'};
  const t=sheet.tables.add(`A${start}:${last}${end}`,true,name);
  t.style='TableStyleMedium2';t.showFilterButton=true;
  widths.forEach((w,i)=>sheet.getRange(`${col(i)}1:${col(i)}${end}`).format.columnWidth=w);
  sheet.getRange(`A${start}:${last}${start}`).format={fill:navy,font:{...font,bold:true,color:'#FFFFFF'},horizontalAlignment:'center',rowHeight:58};
  rows.forEach((r,i)=>{
    const lines=Math.max(...keys.map((k,j)=>String(r[k]??'').split('\n').reduce((n,t)=>n+Math.max(1,Math.ceil(t.length/(widths[j]*.82))),0)));
    sheet.getRange(`A${i+start+1}:${last}${i+start+1}`).format={fill:i%2?pale:'#FFFFFF',rowHeight:Math.max(29,lines*16+12)};
  });
  sheet.showGridLines=false;
  sheet.freezePanes.freezeRows(start);sheet.freezePanes.freezeColumns(2);
  return {end,last};
}
function terms(sheet,rows,name){
  const data=rows.map(r=>({Group:r.Group,Library:r.Library,'Term ID':r.Term_ID,Term:r.Term,
    'BH FDR':r.BH_FDR,'Fold enrichment':r.Fold_enrichment,'Overlap genes':r.Overlap,
    'Query genes':r.Query_size,'Gene ratio':r.Gene_ratio,'Term background':r.Term_background,
    'Background genes':r.Background_size,'Raw P':r.Raw_P,'Joint BH FDR':r.Joint_BH_FDR_six_queries,
    'FDR < 0.05':r.Significant_FDR05?'Yes':'No',Genes:r.Genes,'Term URL':r.Source_URL}));
  const t=table(sheet,data,[23,12,19,58,16,18,16,14,16,18,18,16,17,15,78,66],1,name);
  for(const c of ['E','L','M'])sheet.getRange(`${c}2:${c}${t.end}`).setNumberFormat('0.000E+00');
  sheet.getRange(`F2:F${t.end}`).setNumberFormat('0.00');
  sheet.getRange(`I2:I${t.end}`).setNumberFormat('0.0%');
  sheet.getRange(`E2:E${t.end}`).conditionalFormats.addCustom('AND(ISNUMBER(E2),E2<0.05)',{fill:'#DCEEDB',font:{bold:true,color:'#285C30'}});
  sheet.getRange(`P2:P${t.end}`).format.font={...font,color:'#176A99'};
  return t;
}
terms(sig,d.significant,'SignificantEnrichment');
terms(plot,d.comparison,'ComparedEnrichmentTerms');
const genes=d.input_genes.map(r=>({Group:r.Group,Gene:r.Gene,'Human symbol':r.Human_symbol,
  'Submitted ID':r.Submitted_ID,'Corrected significant windows (±kb)':r.Corrected_windows_kb,
  'Uncorrected significant windows (±kb)':r.Uncorrected_windows_kb,
  'Corrected regulatory direction':r.Corrected_direction,'Uncorrected regulatory direction':r.Uncorrected_direction,
  'Mapping status':r.Mapping_status,'Mapped Ensembl IDs':r.Mapped_Ensembl_IDs}));
table(inputs,genes,[23,18,19,28,30,30,28,28,22,36],1,'EnrichmentInputGenes');
inputs.getRange('I2:I331').conditionalFormats.addCustom('I2<>"Mapped"',{fill:'#FFF0D9',font:{bold:true,color:'#85580C'}});

const themes={Shared:'Chemokine signaling, acetylcholine signaling and angiogenesis',
  'Corrected only':'Neurotransmitter receptor internalization and GABA transport',
  'Uncorrected only':'No terms pass the primary FDR threshold'};
const groups=['Shared','Corrected only','Uncorrected only'];
const summary=groups.map(group=>({Group:group,'Input genes':0,'Mapped genes':0,
  'GO BP terms':0,'KEGG terms':0,'Main enriched themes':themes[group]}));
table(overview,summary,[24,16,17,17,17,64],6,'GroupEnrichmentSummary');
overview.tabColor=navy;
overview.getRange('A2').values=[['Pathway enrichment by body-mass correction status']];
overview.getRange('A2').format={font:{...font,bold:true,size:17},wrapText:false,rowHeight:33};
overview.getRange('A3').values=[['Regulatory activity and MLS. Significant genes pooled across TSS ±20, ±50, ±100, ±200, ±500 and ±1000 kb.']];
overview.getRange('A4').values=[['Shared = significant in both models in at least one window; the significant windows need not be the same.']];
overview.getRange('A3:F4').format={font:{...font,size:10,color:'#526A7A'},wrapText:false,rowHeight:24};
for(let i=0;i<3;i++){
  const r=7+i;
  overview.getRange(`B${r}`).formulas=[[`=COUNTIF('Input genes'!A2:A331,A${r})`]];
  overview.getRange(`C${r}`).formulas=[[`=COUNTIFS('Input genes'!A2:A331,A${r},'Input genes'!I2:I331,"Mapped")`]];
  overview.getRange(`D${r}`).formulas=[[`=COUNTIFS('Significant terms'!A2:A28,A${r},'Significant terms'!B2:B28,"GO:BP")`]];
  overview.getRange(`E${r}`).formulas=[[`=COUNTIFS('Significant terms'!A2:A28,A${r},'Significant terms'!B2:B28,"KEGG")`]];
}
overview.getRange('B7:E9').format.horizontalAlignment='center';
overview.getRange('B7:E9').setNumberFormat('#,##0');
const notes=[
  'Significant = BH FDR < 0.05, calculated separately within each group and library (six correction families).',
  'Same effective background for every query: 19,196 mapped genes tested in at least one regulatory window.',
  'Corrected = log₁₀(body mass) covariate. Human, macaque and bowhead whale are excluded from the source correlations.',
  'MS4A20 and ZFP865 lack usable human annotation mapping; 71 of 73 corrected-only genes enter enrichment.',
  'Different significance across groups does not establish a statistically significant difference between the groups.',
  'GO terms share genes. Overlapping regulatory windows can give nearby genes correlated signals.',
  'Use Significant terms for all 27 hits, Plot terms for cross-group values and Methods for analysis details.',
];
notes.forEach((text,i)=>{
  overview.getRange(`A${12+i}`).values=[[text]];
  overview.getRange(`A${12+i}:F${12+i}`).format={font:{...font,size:10,color:'#526A7A'},wrapText:false,rowHeight:25};
});
const methodRows=[...d.notes,
  {Item:'Files',Details:'Significant terms and comparison plot accompany this workbook. all_terms.tsv.gz retains every tested term and both FDR definitions. Requests, API responses and scripts are saved beside the report.'},
  {Item:'Direction labels',Details:'Input genes retain source regulatory directions. ** = gene-level FDR < 0.05; * = nominal P < 0.05 only. The enrichment groups combine both positive and negative directions.'},
];
table(methods,methodRows,[34,142],1,'EnrichmentMethods');
methods.freezePanes.freezeColumns(1);
wb.recalculate();
const expected=[[143,143,18,2],[73,71,7,0],[114,114,0,0]];
if(JSON.stringify(overview.getRange('B7:E9').values)!==JSON.stringify(expected))throw new Error('Summary formulas differ from validated enrichment counts');
const checks=[];
checks.push((await wb.inspect({kind:'region',sheetId:'Comparison',range:'A6:F9',maxChars:3500})).ndjson);
checks.push((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30}})).ndjson);
await fs.writeFile(path.join(dir,'workbook_checks.txt'),checks.join('\n'));
const previews=[['Comparison','A1:F19','summary_preview.png'],['Significant terms','A1:I9','significant_preview.png'],
  ['Significant terms','L1:P6','term_genes_preview.png'],['Plot terms','A1:I10','comparison_preview.png'],
  ['Input genes','A1:J9','genes_preview.png'],['Methods','A1:B9','methods_preview.png'],
  ['Methods','A10:B16','methods_end_preview.png']];
for(const [sheetName,range,file] of previews){
  const blob=await wb.render({sheetName,range,scale:1,format:'png'});
  await fs.writeFile(path.join(dir,file),new Uint8Array(await blob.arrayBuffer()));
}
const out=await SpreadsheetFile.exportXlsx(wb);
await out.save(path.join(dir,'Regulatory_pathway_enrichment_comparison.xlsx'));
console.log('Exported 5 sheets: comparison, 27 significant terms, 36 comparison rows, 330 input genes, methods.');
