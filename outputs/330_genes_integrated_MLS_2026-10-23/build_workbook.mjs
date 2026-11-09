import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const dir=path.dirname(fileURLToPath(import.meta.url));
const d=JSON.parse(await fs.readFile(path.join(dir,'results.json'),'utf8'));
const wb=Workbook.create();
const overview=wb.worksheets.add('330 genes');
const regulatory=wb.worksheets.add('Regulatory windows');
const details=wb.worksheets.add('Correlation details');
const evidence=wb.worksheets.add('Aging evidence');
const roles=wb.worksheets.add('Role details');
const methods=wb.worksheets.add('Methods');
const font={name:'Arial',size:11,color:'#23394B'};
const navy='#173B55', pale='#F0F5F8';
const col=n=>{let s='';for(n++;n>0;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};

function table(sheet,rows,widths,start,name){
  const keys=Object.keys(rows[0]),last=col(keys.length-1),end=start+rows.length;
  sheet.getRange(`A${start}:${last}${end}`).values=[keys,...rows.map(r=>keys.map(k=>r[k]??null))];
  sheet.getRange(`A${start}:${last}${end}`).format={font,wrapText:true,verticalAlignment:'center'};
  const t=sheet.tables.add(`A${start}:${last}${end}`,true,name);
  t.style='TableStyleMedium2';
  t.showFilterButton=true;
  widths.forEach((w,i)=>sheet.getRange(`${col(i)}1:${col(i)}${end}`).format.columnWidth=w);
  sheet.getRange(`A${start}:${last}${start}`).format={fill:navy,font:{...font,bold:true,color:'#FFFFFF'},horizontalAlignment:'center',rowHeight:62,
    borders:{insideVertical:{style:'thin',color:'#FFFFFF'}}};
  rows.forEach((r,i)=>{
    const lines=Math.max(...keys.map((k,j)=>String(r[k]??'').split('\n').reduce((n,t)=>n+Math.max(1,Math.ceil(t.length/(widths[j]*.83))),0)));
    sheet.getRange(`A${i+start+1}:${last}${i+start+1}`).format={fill:i%2?pale:'#FFFFFF',rowHeight:Math.max(28,lines*15+10)};
  });
  sheet.showGridLines=false;
  sheet.freezePanes.freezeRows(start);
  sheet.freezePanes.freezeColumns(2);
  return {end,last};
}

function directionFormat(sheet,range){
  sheet.getRange(range).conditionalFormats.add('beginsWith',{text:'Positive',format:{fill:'#FBE5DE',font:{color:'#8A3024'}}});
  sheet.getRange(range).conditionalFormats.add('beginsWith',{text:'Negative',format:{fill:'#DEECF6',font:{color:'#235C85'}}});
  sheet.getRange(range).conditionalFormats.add('beginsWith',{text:'Not ',format:{fill:'#ECEFF1',font:{color:'#65757F'}}});
}

const rows=d.overview.map(r=>({Gene:r.Gene,'Human symbol':r.Human_symbol,
  'Regulatory significance':r.Regulatory_significance,
  'Significant windows: corrected (±kb)':r.Corrected_windows_kb,
  'Significant windows: uncorrected (±kb)':r.Uncorrected_windows_kb,
  Experiment:r.Experiment,GenAge:r.GenAge,CellAge:r.CellAge,GWAS:r.GWAS,
  'Reported aging / longevity role':r.Reported_aging_longevity_role,
  'Regulatory MLS: corrected':r.Regulatory_corrected,
  'Regulatory MLS: uncorrected':r.Regulatory_uncorrected,
  'Gene activity TSS ±2 kb':r.tss2kb,'Gene activity body + 2 kb':r.body,
  'Transcriptome MLS':r.rna,'Proteome MLS':r.protein}));
const ov=table(overview,rows,[16,19,23,25,25,13,12,12,12,47,23,23,22,23,21,21],6,'CombinedGeneOverview');
overview.tabColor=navy;
overview.getRange('A2').values=[['330 regulatory hits: body-mass correction, aging evidence and MLS correlations']];
overview.getRange('A2').format={font:{...font,bold:true,size:16},wrapText:false,rowHeight:31};
overview.getRange('A3').formulas=[[`="Significant in both: "&COUNTIFS(C7:C${ov.end},"Both")&"     Corrected only: "&COUNTIFS(C7:C${ov.end},"Corrected only")&"     Uncorrected only: "&COUNTIFS(C7:C${ov.end},"Uncorrected only")`]];
overview.getRange('A3').format={font:{...font,bold:true},wrapText:false,rowHeight:25};
overview.getRange('A4').values=[['** FDR < 0.05   |   * nominal P < 0.05 only   |   Unmarked directions are nonsignificant   |   Evidence: Yes / not found in compiled reference']];
overview.getRange('A5').values=[['Corrected = log₁₀(body mass) covariate. Gene activity, transcriptome and proteome retain body-mass correction; gene-body activity also adjusts for log₁₀(gene length).']];
overview.getRange('A4:P5').format={font:{...font,size:10,color:'#526A7A'},wrapText:false,rowHeight:23};
overview.getRange(`C7:P${ov.end}`).format.horizontalAlignment='center';
overview.getRange(`J7:J${ov.end}`).format.horizontalAlignment='left';
for(const [label,color,ink] of [['Both','#DCEEDB','#285C30'],['Corrected only','#E2EAF6','#24537A'],['Uncorrected only','#FFF0D9','#85580C']]){
  overview.getRange(`C7:C${ov.end}`).conditionalFormats.addCustom(`C7="${label}"`,{fill:color,font:{bold:true,color:ink}});
}
overview.getRange(`F7:I${ov.end}`).conditionalFormats.add('containsText',{text:'Yes',format:{fill:'#DCEEDB',font:{bold:true,color:'#285C30'}}});
overview.getRange(`F7:I${ov.end}`).conditionalFormats.add('containsText',{text:'No',format:{font:{color:'#86939D'}}});
directionFormat(overview,`K7:P${ov.end}`);

const regRows=d.regulatory.map(r=>({Gene:r.Gene,'Human symbol':r.Human_symbol,Model:r.Model,
  'Body mass correction':r.Body_mass_correction,'TSS flank (±kb)':r.Window_flank_kb,
  'Full width (kb)':r.Full_width_kb,Direction:r.Direction,'Correlation r':r.Correlation_r,P:r.P,
  'BH FDR within window':r.BH_FDR_window,'Significant FDR < 0.05':r.Significant_FDR05,
  'Species N':r.N_species,DF:r.DF,'MLS coefficient':r.MLS_coefficient,
  'FDR family':r.FDR_family,'Species IDs':r.Species_IDs}));
const re=table(regulatory,regRows,[16,19,20,20,18,18,18,17,17,21,22,14,12,22,68,100],1,'RegulatoryWindows');
regulatory.getRange(`H2:H${re.end}`).setNumberFormat('+0.000;-0.000;0.000');
regulatory.getRange(`I2:J${re.end}`).setNumberFormat('0.000E+00');
regulatory.getRange(`N2:N${re.end}`).setNumberFormat('+0.0000;-0.0000;0.0000');
for(const c of ['E','F','L','M'])regulatory.getRange(`${c}2:${c}${re.end}`).setNumberFormat('0');
regulatory.getRange(`J2:J${re.end}`).conditionalFormats.addCustom('AND(ISNUMBER(J2),J2<0.05)',{fill:'#DCEEDB',font:{bold:true}});
directionFormat(regulatory,`G2:G${re.end}`);

const detailRows=d.details.map(r=>({Gene:r.Gene,'Human symbol':r.Human_symbol,Method:r.Method,
  'Body mass correction':r.Body_mass_correction,Covariates:r.Covariates,Direction:r.Direction,
  'Partial r':r.Partial_r,P:r.P,'BH FDR':r.BH_FDR,'Species N':r.N_species,DF:r.DF,
  'MLS coefficient':r.MLS_coefficient,Status:r.Status,'Source gene symbol':r.Source_gene_symbol,
  'Symbol matching':r.Symbol_matching,'FDR family':r.FDR_family,'Species IDs':r.Species_IDs}));
const de=table(details,detailRows,[16,19,35,20,48,18,16,17,17,14,12,22,44,22,40,54,100],1,'AssayCorrelations');
details.getRange(`G2:G${de.end}`).setNumberFormat('+0.000;-0.000;0.000');
details.getRange(`H2:I${de.end}`).setNumberFormat('0.000E+00');
details.getRange(`L2:L${de.end}`).setNumberFormat('+0.0000;-0.0000;0.0000');
details.getRange(`J2:K${de.end}`).setNumberFormat('0');
details.getRange(`I2:I${de.end}`).conditionalFormats.addCustom('AND(ISNUMBER(I2),I2<0.05)',{fill:'#DCEEDB',font:{bold:true}});
directionFormat(details,`F2:F${de.end}`);

const evidenceRows=d.evidence.map(r=>({...r,'P value':r['P value']!==''&&Number.isFinite(Number(r['P value']))?Number(r['P value']):null}));
const ev=table(evidence,evidenceRows,[18,17,33,21,28,19,52,34,46,100,19,24,29,26,100,23,100],1,'AgingEvidenceRecords');
evidence.getRange(`K2:K${ev.end}`).setNumberFormat('0.000E+00');
evidence.getRange(`O2:O${ev.end}`).format.font={...font,color:'#176A99'};
const roleRows=d.roles.map(r=>({Gene:r.Gene,'Human symbol':r.Human_symbol,'Reported role':r.Reported_role,
  'Longevity role':r.Longevity_role,'Cellular senescence role':r.Cellular_senescence_role,
  'Directional evidence':r.Directional_evidence?'Yes':'No','Evidence basis':r.Evidence_basis,
  Interpretation:r.Interpretation,'Source URLs':r.Source_URLs}));
table(roles,roleRows,[16,19,47,31,31,20,100,95,110],1,'AgingRoleDetails');
table(methods,d.notes,[36,135],1,'AnalysisMethodNotes');
methods.freezePanes.freezeColumns(1);

wb.recalculate();
const checks=[];
checks.push((await wb.inspect({kind:'region',sheetId:'330 genes',range:'A2:C12',maxChars:2500})).ndjson);
checks.push((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30}})).ndjson);
await fs.writeFile(path.join(dir,'workbook_checks.txt'),checks.join('\n'));
if(!String(overview.getRange('A3').values[0][0]).includes('143'))throw new Error('Summary formula did not recalculate');
const fullPreviews=[['330 genes','A1:P16','overview_preview.png'],['330 genes',`A${ov.end-7}:P${ov.end}`,'overview_end_preview.png'],
  ['Regulatory windows','A1:K10','regulatory_preview.png'],['Correlation details','A1:M9','correlations_preview.png'],
  ['Aging evidence','A1:F8','evidence_preview.png'],['Aging evidence','J1:Q5','evidence_links_preview.png'],
  ['Role details','A1:G8','roles_preview.png'],['Methods','A1:B12','methods_preview.png'],['Methods','A13:B22','methods_end_preview.png']];
const previews=process.argv.includes('--overview-only')?fullPreviews.filter(p=>p[0]==='330 genes'):fullPreviews;
for(const gene of ['EIF4EBP1','IL11','SIRT6']){
  const row=7+d.overview.findIndex(r=>r.Gene===gene);
  previews.push(['330 genes',`A${row}:P${row}`,`${gene}_preview.png`]);
}
for(const [sheet,range,file] of previews){
  const blob=await wb.render({sheetName:sheet,range,scale:1,format:'png'});
  await fs.writeFile(path.join(dir,file),new Uint8Array(await blob.arrayBuffer()));
}
const output=await SpreadsheetFile.exportXlsx(wb);
await output.save(path.join(dir,'330_genes_corrected_and_uncorrected_MLS.xlsx'));
console.log('Exported workbook: 330 genes, 3,960 regulatory rows, 1,320 assay rows, '+d.evidence.length+' evidence rows.');
