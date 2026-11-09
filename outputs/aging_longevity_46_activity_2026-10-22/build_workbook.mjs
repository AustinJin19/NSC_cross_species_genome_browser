import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const dir=path.dirname(fileURLToPath(import.meta.url));
const d=JSON.parse(await fs.readFile(path.join(dir,'results.json'),'utf8'));
const wb=Workbook.create();
const summary=wb.worksheets.add('Summary');
const genes=wb.worksheets.add('46 genes');
const values=wb.worksheets.add('Species values');
const methods=wb.worksheets.add('Methods');
const navy='#173B55',pale='#EEF4F8';
const font={name:'Arial',size:11,color:'#243B4B'};
const col=i=>{let s='';for(i++;i>0;i=Math.floor((i-1)/26))s=String.fromCharCode(65+(i-1)%26)+s;return s;};

function table(sheet, rows, keys, widths, start=1, name) {
  const last=col(keys.length-1),end=start+rows.length;
  sheet.getRange(`A${start}:${last}${end}`).values=[keys,...rows.map(r=>keys.map(k=>r[k]??null))];
  sheet.getRange(`A${start}:${last}${end}`).format={font,wrapText:true,verticalAlignment:'center'};
  if(name)sheet.tables.add(`A${start}:${last}${end}`,true,name);
  widths.forEach((w,i)=>sheet.getRange(`${col(i)}1:${col(i)}${end}`).format.columnWidth=w);
  sheet.getRange(`A${start}:${last}${start}`).format={fill:navy,font:{...font,bold:true,color:'#FFFFFF'},horizontalAlignment:'center',rowHeight:48};
  rows.forEach((r,i)=>{
    const lines=Math.max(...keys.map((k,j)=>String(r[k]??'').split('\n').reduce((n,t)=>n+Math.max(1,Math.ceil(t.length/(widths[j]*.82))),0)));
    sheet.getRange(`A${start+i+1}:${last}${start+i+1}`).format={fill:i%2?pale:'#FFFFFF',rowHeight:Math.max(29,lines*15+10)};
  });
  sheet.showGridLines=false;
  return end;
}
const rows=d.wide.map(r=>({
  Gene:r.Gene,
  'TSS ±2 kb direction':r.tss2kb_Direction,'TSS ±2 kb partial r':r.tss2kb_Partial_r,
  'TSS ±2 kb P':r.tss2kb_P,'TSS ±2 kb BH FDR':r.tss2kb_Genomewide_BH_FDR,'TSS ±2 kb N':r.tss2kb_N_species,
  'Body + 2 kb direction':r.body_Direction,'Body + 2 kb partial r':r.body_Partial_r,
  'Body + 2 kb P':r.body_P,'Body + 2 kb BH FDR':r.body_Genomewide_BH_FDR,'Body + 2 kb N':r.body_N_species,
  'Regulatory positive windows (±kb)':r.Regulatory_significant_windows_kb,
  'Aging / longevity sources':r.Resources,
  'Gene name':r.Gene_name,
}));
const ge=table(genes,rows,Object.keys(rows[0]),[16,18,15,14,16,12,19,16,14,17,12,27,33,60],1,'GeneActivityResults');
genes.freezePanes.freezeRows(1);genes.freezePanes.freezeColumns(1);genes.tabColor=navy;
for(const c of ['C','H']) {
  genes.getRange(`${c}2:${c}${ge}`).setNumberFormat('+0.000;-0.000;0.000');
  genes.getRange(`${c}2:${c}${ge}`).conditionalFormats.add('colorScale',{
    colors:['#67A9CF','#FFFFFF','#EF8A62'],thresholds:[-.65,0,.65]});
}
for(const c of ['D','E','I','J'])genes.getRange(`${c}2:${c}${ge}`).setNumberFormat('0.0000');
for(const c of ['D','I'])genes.getRange(`${c}2:${c}${ge}`).conditionalFormats.add('cellIs',{
  operator:'lessThan',formula:.05,format:{fill:'#FFF0CF',font:{bold:true}}});
for(const c of ['F','K']){
  genes.getRange(`${c}2:${c}${ge}`).setNumberFormat('0');
  genes.getRange(`${c}2:${c}${ge}`).format.horizontalAlignment='center';
}
genes.getRange(`L2:L${ge}`).format.horizontalAlignment='center';

summary.showGridLines=false;summary.tabColor=navy;
summary.getRange('A1:G28').format={font,rowHeight:26,verticalAlignment:'center'};
summary.getRange('A2').values=[['Gene activity in the 46 aging/longevity overlaps']];
summary.getRange('A2').format={font:{...font,bold:true,size:16},rowHeight:30,wrapText:false};
summary.getRange('A3').values=[['log₁₀(activity) versus log₁₀(MLS); body-mass adjusted; gene body also length adjusted']];
summary.getRange('A3').format={font:{...font,size:10,color:'#536775'},wrapText:false};
const counts=d.summary.map(r=>({Method:r.Method,'Positive':r.Positive,'Negative':r.Negative,
  'Nominal positive (P < .05)':r.Nominal_positive,'Nominal negative (P < .05)':r.Nominal_negative,
  'BH FDR < .05':r.FDR_positive+r.FDR_negative,'Tested genome-wide':r.Genomewide_tests}));
table(summary,counts,Object.keys(counts[0]),[34,15,15,20,20,17,23],5);
for(const [row,rc,pc,qc]of [[6,'C','D','E'],[7,'H','I','J']]) {
  summary.getRange(`B${row}:F${row}`).formulas=[[
    `=COUNTIF('46 genes'!${rc}2:${rc}47,">0")`,
    `=COUNTIF('46 genes'!${rc}2:${rc}47,"<0")`,
    `=COUNTIFS('46 genes'!${rc}2:${rc}47,">0",'46 genes'!${pc}2:${pc}47,"<0.05")`,
    `=COUNTIFS('46 genes'!${rc}2:${rc}47,"<0",'46 genes'!${pc}2:${pc}47,"<0.05")`,
    `=COUNTIF('46 genes'!${qc}2:${qc}47,"<0.05")`
  ]];
}
summary.getRange('B6:G7').setNumberFormat('#,##0');
summary.getRange('B6:G7').format.horizontalAlignment='center';
summary.getRange('A9').values=[['Nominal associations — none pass genome-wide FDR < 0.05']];
summary.getRange('A9').format={font:{...font,bold:true},wrapText:false};
const nominal=d.nominal.map(r=>({Gene:r.Gene,Region:r.Method==='TSS ±2 kb'?'TSS ±2 kb':'Body + 2 kb',
  Direction:r.Direction,'Partial r':r.Partial_r,'P':r.P,'BH FDR':r.Genomewide_BH_FDR,'Species N':r.N_species}));
table(summary,nominal,Object.keys(nominal[0]),[34,15,15,20,20,17,23],11);
summary.getRange('D12:D17').setNumberFormat('+0.000;-0.000;0.000');
summary.getRange('E12:F17').setNumberFormat('0.0000');
summary.getRange('G12:G17').setNumberFormat('0');
summary.getRange('A19').values=[['Direction agreement across methods']];
summary.getRange('A19').format={font:{...font,bold:true},wrapText:false};
const agreement=Object.entries(d.direction_counts).map(([direction,n])=>[direction,n]);
summary.getRange('A20:B23').values=agreement;
summary.getRange('A20:A23').format.wrapText=false;
summary.getRange('B20:B23').setNumberFormat('0');
summary.getRange('D20').values=[['Order: TSS ±2 kb / gene body + 2 kb']];
summary.getRange('D21').values=[['All 46 were positive in the regulatory peak-count search.']];
summary.getRange('D22').values=[['Direction alone does not establish statistical significance.']];
summary.getRange('D20:D22').format={wrapText:false,font:{...font,size:10,color:'#536775'}};
summary.getRange('A25').values=[['Human, macaque and bowhead whale excluded; ≥50% available region required.']];
summary.getRange('A26').values=[[`${d.same_species_genes}/46 genes use identical species sets across methods; per-gene N is retained.`]];
summary.getRange('A27').values=[['Activity is integrated CPM·bp within the complete region, not peak counts or RNA expression.']];
summary.getRange('A28').values=[['Results are exploratory cross-species associations; see Methods for sources and assumptions.']];
summary.getRange('A25:A28').format={wrapText:false,font:{...font,size:10,color:'#536775'}};

const sr=d.species_values.map(r=>({Gene:r.Gene,Method:r.Method,Species:r.Species,
  'Median activity (CPM·bp)':r.Median_activity_CPM_bp,'log₁₀(activity)':r.Log10_median_activity,
  'MLS (years)':r.MLS_years,'Body mass (g)':r.Body_mass_g,'Gene length (kb)':r.Gene_length_kb,
  'Available region':r.Available_fraction,'Samples':r.Sample_count,'Fit status':r.Status,
  Locus:r.Locus,'Scientific name':r.Scientific_name,'Species ID':r.Species_ID}));
const se=table(values,sr,Object.keys(sr[0]),[16,31,29,23,19,16,20,21,19,14,30,89,32,32],1,'SpeciesActivityData');
values.freezePanes.freezeRows(1);values.freezePanes.freezeColumns(2);
values.getRange(`D2:D${se}`).setNumberFormat('#,##0.00');
values.getRange(`E2:E${se}`).setNumberFormat('0.0000');
values.getRange(`F2:H${se}`).setNumberFormat('#,##0.00');
values.getRange(`I2:I${se}`).setNumberFormat('0.0%');
values.getRange(`J2:J${se}`).setNumberFormat('0');
values.getRange(`J2:J${se}`).format.horizontalAlignment='center';
values.getRange(`K2:K${se}`).conditionalFormats.add('containsText',{text:'below',format:{fill:'#FFF0CF'}});
table(methods,d.notes,['Item','Details'],[32,128],1,'AnalysisMethods');
methods.freezePanes.freezeRows(1);

wb.recalculate();
const checks=[];
checks.push((await wb.inspect({kind:'region',sheetId:'Summary',range:'A5:G7',maxChars:2500})).ndjson);
checks.push((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30}})).ndjson);
await fs.writeFile(path.join(dir,'workbook_checks.txt'),checks.join('\n'));
const previews=process.argv.includes('--changed-only')
  ? [['46 genes','A1:M12','genes_preview.png'],['Species values','A1:K9','species_preview.png']]
  : [['Summary','A1:G28','summary_preview.png'],['46 genes','A1:M12','genes_preview.png'],
    ['46 genes','A34:N47','genes_bottom_preview.png'],['Species values','A1:K9','species_preview.png'],['Methods','A1:B15','methods_preview.png']];
for(const [sheet,range,file]of previews) {
  const preview=await wb.render({sheetName:sheet,range,scale:1.3,format:'png'});
  await fs.writeFile(path.join(dir,file),new Uint8Array(await preview.arrayBuffer()));
}
const output=await SpreadsheetFile.exportXlsx(wb);
await output.save(path.join(dir,'46_genes_activity_and_MLS.xlsx'));
console.log(checks[0]);
console.log('Saved workbook and previews.');
