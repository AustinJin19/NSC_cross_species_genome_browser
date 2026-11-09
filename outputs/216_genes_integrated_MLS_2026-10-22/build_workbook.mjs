import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const dir=path.dirname(fileURLToPath(import.meta.url));
const d=JSON.parse(await fs.readFile(path.join(dir,'results.json'),'utf8'));
const wb=Workbook.create();
const overview=wb.worksheets.add('216 genes');
const details=wb.worksheets.add('Correlation details');
const evidence=wb.worksheets.add('Aging evidence');
const roles=wb.worksheets.add('Role details');
const methods=wb.worksheets.add('Methods');
const font={name:'Arial',size:11,color:'#23394B'};
const navy='#173B55',pale='#F0F5F8';
const col=n=>{let s='';for(n++;n>0;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};

function table(sheet,rows,widths,start,name){
  const keys=Object.keys(rows[0]),last=col(keys.length-1),end=start+rows.length;
  sheet.getRange(`A${start}:${last}${end}`).values=[keys,...rows.map(r=>keys.map(k=>r[k]??null))];
  sheet.getRange(`A${start}:${last}${end}`).format={font,wrapText:true,verticalAlignment:'center'};
  sheet.tables.add(`A${start}:${last}${end}`,true,name);
  widths.forEach((w,i)=>sheet.getRange(`${col(i)}1:${col(i)}${end}`).format.columnWidth=w);
  sheet.getRange(`A${start}:${last}${start}`).format={fill:navy,font:{...font,bold:true,color:'#FFFFFF'},horizontalAlignment:'center',rowHeight:58};
  rows.forEach((r,i)=>{
    const lines=Math.max(...keys.map((k,j)=>String(r[k]??'').split('\n').reduce((n,t)=>n+Math.max(1,Math.ceil(t.length/(widths[j]*.83))),0)));
    sheet.getRange(`A${i+start+1}:${last}${i+start+1}`).format={fill:i%2?pale:'#FFFFFF',rowHeight:Math.max(28,lines*15+10)};
  });
  sheet.showGridLines=false;sheet.freezePanes.freezeRows(start);sheet.freezePanes.freezeColumns(1);
  return {end,last};
}

const rows=d.overview.map(r=>({
  Gene:r.Gene,'Human symbol':r.Human_symbol,'Significant regulatory windows (±kb)':r.Significant_windows_kb,
  Experiment:r.Experiment,GenAge:r.GenAge,CellAge:r.CellAge,GWAS:r.GWAS,
  'Reported aging / longevity role':r.Reported_aging_longevity_role,
  'Regulatory MLS':r.Regulatory,'Gene activity TSS ±2 kb':r.tss2kb,
  'Gene activity body + 2 kb':r.body,'Transcriptome MLS':r.rna,'Proteome MLS':r.protein,
}));
const ov=table(overview,rows,[16,19,27,13,12,12,12,47,19,21,23,20,20],5,'IntegratedGeneOverview');
overview.tabColor=navy;
overview.getRange('A2').values=[['216 regulatory hits: aging evidence and MLS correlation directions']];
overview.getRange('A2').format={font:{...font,bold:true,size:16},wrapText:false,rowHeight:31};
overview.getRange('A3').values=[['** FDR < 0.05  |  * nominal P < 0.05 only  |  Unmarked directions are nonsignificant  |  Evidence: Yes / not found in compiled reference']];
overview.getRange('A3').format={font:{...font,size:10,color:'#526A7A'},wrapText:false,rowHeight:25};
overview.getRange(`C6:M${ov.end}`).format.horizontalAlignment='center';
overview.getRange(`H6:H${ov.end}`).format.horizontalAlignment='left';
overview.getRange(`D6:G${ov.end}`).conditionalFormats.add('containsText',{text:'Yes',format:{fill:'#DCEEDB',font:{bold:true,color:'#285C30'}}});
overview.getRange(`D6:G${ov.end}`).conditionalFormats.add('containsText',{text:'No',format:{font:{color:'#86939D'}}});
overview.getRange(`I6:M${ov.end}`).conditionalFormats.add('beginsWith',{text:'Positive',format:{fill:'#FBE5DE',font:{color:'#8A3024'}}});
overview.getRange(`I6:M${ov.end}`).conditionalFormats.add('beginsWith',{text:'Negative',format:{fill:'#DEECF6',font:{color:'#235C85'}}});
overview.getRange(`I6:M${ov.end}`).conditionalFormats.add('containsText',{text:'Not tested',format:{fill:'#ECEFF1',font:{color:'#65757F'}}});

const dr=d.details.map(r=>({Gene:r.Gene,'Human symbol':r.Human_symbol,Method:r.Method,Direction:r.Direction,
  'Partial r':r.Partial_r,P:r.P,'BH FDR':r.BH_FDR,'Species N':r.N_species,Status:r.Status,
  'Representative regulatory window (±kb)':r.Representative_window_kb,'Source gene symbol':r.Matched_source_symbol,
  'Symbol matching':r.Symbol_matching,'MLS coefficient':r.MLS_beta,DF:r.DF,
  'FDR family':r.FDR_family,'Species IDs':r.Species_IDs}));
const de=table(details,dr,[16,19,35,18,14,16,16,14,41,26,22,40,22,12,54,92],1,'CorrelationStatistics');
details.getRange(`E2:E${de.end}`).setNumberFormat('+0.000;-0.000;0.000');
details.getRange(`F2:G${de.end}`).setNumberFormat('0.000E+00');
details.getRange(`M2:M${de.end}`).setNumberFormat('+0.0000;-0.0000;0.0000');
for(const c of ['H','J','N']){
  details.getRange(`${c}2:${c}${de.end}`).setNumberFormat('0');
  details.getRange(`${c}2:${c}${de.end}`).format.horizontalAlignment='center';
}
details.getRange(`G2:G${de.end}`).conditionalFormats.addCustom('AND(ISNUMBER(G2),G2<0.05)',
  {fill:'#DCEEDB',font:{bold:true}});

const ev=table(evidence,d.evidence,[18,17,33,21,28,19,52,34,46,100,19,24,29,26,95,23,85],1,'AgingEvidenceRecords');
// GWAS P values are text in the source TSV; preserve numeric ones for sorting.
d.evidence.forEach((r,i)=>{if(r['P value']!==''&&Number.isFinite(Number(r['P value'])))evidence.getRange(`K${i+2}`).values=[[Number(r['P value'])]];});
evidence.getRange(`K2:K${ev.end}`).setNumberFormat('0.000E+00');
evidence.getRange(`O2:O${ev.end}`).format.font={...font,color:'#176A99'};
const roleRows=d.roles.map(r=>({Gene:r.Gene,'Human symbol':r.Human_symbol,'Reported role':r.Reported_role,
  'Longevity role':r.Longevity_role,'Cellular senescence role':r.Cellular_senescence_role,
  'Directional evidence':r.Directional_evidence?'Yes':'No','Evidence basis':r.Evidence_basis,
  Interpretation:r.Interpretation,'Source URLs':r.Source_URLs}));
table(roles,roleRows,[16,19,47,31,31,20,95,95,105],1,'AgingRoleDetails');
table(methods,d.notes,[36,128],1,'AnalysisMethodNotes');

wb.recalculate();
const checks=[];
checks.push((await wb.inspect({kind:'region',sheetId:'216 genes',range:'A5:M10',maxChars:2500})).ndjson);
checks.push((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30}})).ndjson);
await fs.writeFile(path.join(dir,'workbook_checks.txt'),checks.join('\n'));
const previews=process.argv.includes('--changed-only')
  ? [['216 genes','A5:M14','overview_preview.png'],['Correlation details','A1:J11','correlations_preview.png']]
  : [['216 genes','A1:M14','overview_preview.png'],['216 genes',`A${ov.end-9}:M${ov.end}`,'overview_end_preview.png'],
    ['Correlation details','A1:J11','correlations_preview.png'],['Aging evidence','A1:F8','evidence_preview.png'],
    ['Aging evidence','J1:O6','evidence_links_preview.png'],['Methods','A1:B11','methods_preview.png'],['Methods','A12:B20','methods_end_preview.png']];
if(!process.argv.includes('--changed-only'))previews.push(['Role details','A1:G8','roles_preview.png'],['Methods','A21:B24','role_methods_preview.png']);
for(const [sheet,range,file] of previews){
  const blob=await wb.render({sheetName:sheet,range,scale:1.25,format:'png'});
  await fs.writeFile(path.join(dir,file),new Uint8Array(await blob.arrayBuffer()));
}
const output=await SpreadsheetFile.exportXlsx(wb);
await output.save(path.join(dir,'216_genes_aging_evidence_and_MLS.xlsx'));
console.log('Exported 216-gene workbook, 1,080 correlation rows and '+d.evidence.length+' evidence rows.');
