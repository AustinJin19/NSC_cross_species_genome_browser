import fs from 'node:fs/promises';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const dir = new URL('.', import.meta.url).pathname;
const data = JSON.parse(await fs.readFile(`${dir}comparison.json`, 'utf8'));
const wb = Workbook.create();
const main = wb.worksheets.add('Assembly comparison');
const evidence = wb.worksheets.add('Evidence');
const navy = '#193E58', ink = '#203444';
for (const sheet of [main, evidence]) {
  sheet.showGridLines = false;
  sheet.getRange(sheet === main ? 'A1:M53' : 'A1:J39').format.font = {name:'Arial',size:11,color:ink};
  sheet.getRange(sheet === main ? 'A1:M53' : 'A1:J39').format.verticalAlignment='center';
}
main.tabColor=navy;
main.getRange('A2').values=[['Genome assemblies across four resources']];
main.getRange('A2').format.font={name:'Arial',size:16,bold:true,color:navy};
main.getRange('A2:M2').format.rowHeight=29;
main.getRange('A3:M3').format.borders={bottom:{style:'thin',color:'#A8BAC8'}};
main.getRange('A4').values=[['34 browser species, including human, macaque and bowhead whale. Checked 7 November 2026.']];
main.getRange('A4').format.font.italic=true;
main.getRange('A6:E6').values=[['Compared with your assembly','Species present','Matching build','Different build/version','Species absent']];
main.getRange('A7:A9').values=[['TOGA2 catalogue'],['Multiz-470way MSA'],['Cactus-447way MSA']];
const cols=['G','J','M'];
for (let i=0;i<3;i++) {
  const r=i+7,c=cols[i];
  main.getRange(`B${r}:E${r}`).formulas=[[
    `=COUNTA(A19:A52)-E${r}`,
    `=COUNTIFS(${c}19:${c}52,"Matching build")`,
    `=SUM(COUNTIFS(${c}19:${c}52,"Different build/version"),COUNTIFS(${c}19:${c}52,"Different version; chromosomes match"))`,
    `=COUNTIFS(${c}19:${c}52,"Absent")`
  ]];
}
main.getRange('B7:E9').setNumberFormat('0');
main.getRange('B7:E9').format.horizontalAlignment='right';
main.getRange('A7:E9').format.rowHeight=25;
main.getRange('A11').values=[['Matching build means accession, paired RefSeq/GenBank, named reference build or complete scaffold-size agreement. Sequence bases were not compared.']];
main.getRange('A12').values=[['TOGA2 may list several builds. A matching build is selected when available; otherwise a related build is shown. Multiz alternatives are retained.']];
main.getRange('A13').values=[['Eastern gray squirrel is counted as a different version in TOGA2 and Multiz, although all 21 assembled chromosome accessions and lengths match.']];
main.getRange('A14').values=[['Matching across all four']];
main.getRange('B14').formulas=[['=COUNTIFS(G19:G52,"Matching build",J19:J52,"Matching build",M19:M52,"Matching build")']];
main.getRange('A14:B14').format.font.bold=true;
main.getRange('B14').setNumberFormat('0');
main.getRange('A15').values=[['GCA and GCF versions are shown exactly as recorded. Chromosome/scaffold aliases still need checking before transferring genomic intervals.']];

const headers=['Common name','Scientific name','Your assembly','Your accession','TOGA2 assembly selected','TOGA2 accession','TOGA2 versus yours','Multiz-470way assembly','Multiz accession(s)','Multiz versus yours','Cactus-447way assembly','Cactus accession','Cactus versus yours'];
main.getRange('A18:M18').values=[headers];
main.getRange('A19:M52').values=data.rows.map(r=>[
  r.common_name,r.scientific_name,r.local_assembly,r.local_accession,
  r.toga2_assembly_id==='Absent'?'Absent':`${r.toga2_assembly_id}: ${r.toga2_assembly}`,
  r.toga2_accession,r.toga2_status,r.multiz470_assembly,r.multiz470_accession,r.multiz470_status,
  r.cactus447_assembly,r.cactus447_accession,r.cactus447_status
]);
main.getRange('A18:M52').format.wrapText=true;
main.getRange('A18:M18').format.rowHeight=40;
const widths={A:30,B:31,C:39,D:23,E:43,F:23,G:25,H:50,I:24,J:25,K:36,L:23,M:25};
for(const [c,width] of Object.entries(widths)) main.getRange(`${c}1:${c}53`).format.columnWidth=width;
for(let i=0;i<34;i++){
  const values=main.getRange(`A${i+19}:M${i+19}`).values[0];
  const lines=Math.max(...values.map((v,j)=>Math.ceil(String(v||'').length/(Object.values(widths)[j]*1.05))));
  main.getRange(`A${i+19}:M${i+19}`).format.rowHeight=30+15*(lines-1);
}
main.getRange('B19:B52').format.font.italic=true;
const t=main.tables.add('A18:M52',true,'AssemblyComparison');
t.style='TableStyleMedium2';t.showFilterButton=true;
main.freezePanes.freezeRows(18);
main.freezePanes.freezeColumns(2);

evidence.getRange('A2').values=[['Assembly matching evidence and resources']];
evidence.getRange('A2').format.font={name:'Arial',size:16,bold:true,color:navy};
evidence.getRange('A2:J2').format.rowHeight=29;
evidence.getRange('A3:J3').format.borders={bottom:{style:'thin',color:'#A8BAC8'}};
evidence.getRange('A4').values=[['Local source: NSC_genome_assemblies_and_annotations.xlsx, reconciled to the current 34-species inventory.']];
evidence.getRange('F2').formulas=[['=HYPERLINK("https://genome.senckenberg.de/download/TOGA2/assemblies_and_species.tsv","TOGA2 assembly catalogue")']];
evidence.getRange('H2').formulas=[['=HYPERLINK("https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cons470way","Multiz-470way assembly table")']];
evidence.getRange('I2').formulas=[['=HYPERLINK("https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/makeDb/trackDb/human/hg38/cactus447way.html","Cactus-447way assembly table")']];
evidence.getRange('A5:D5').values=[['Species','Local GenBank pair','TOGA2 match evidence','Coordinate interpretation']];
evidence.getRange('F5:J5').values=[['Your assembly resource','TOGA2 genome resource','Multiz assembly resource(s)','Cactus assembly resource(s)','Additional sequence-size evidence']];
evidence.getRange('A6:D39').values=data.rows.map(r=>[r.common_name,r.local_paired_genbank||'Not recorded / not applicable',r.toga2_match_basis,r.coordinate_note]);
evidence.getRange('F6:J39').values=data.rows.map(r=>[r.local_assembly_source,r.toga2_genome_url,r.multiz470_assembly_links,r.cactus447_assembly_links,r.support_url]);
for (const [c,width] of Object.entries({A:30,B:26,C:52,D:70,E:3,F:52,G:60,H:52,I:52,J:60})) evidence.getRange(`${c}1:${c}39`).format.columnWidth=width;
evidence.getRange('A5:D39').format.wrapText=true;
evidence.getRange('F5:J39').format.wrapText=true;
for(let i=0;i<34;i++){
  const values=evidence.getRange(`A${i+6}:J${i+6}`).values[0];
  const widths=[30,26,52,70,3,52,60,52,52,60];
  const lines=Math.max(...values.map((v,j)=>Math.ceil(String(v||'').length/(widths[j]*1.05))));
  evidence.getRange(`A${i+6}:J${i+6}`).format.rowHeight=30+15*(lines-1);
}
evidence.getRange('A5:J5').format.rowHeight=40;
evidence.getRange('F6:J39').format.font.color='#1765A1';
const et=evidence.tables.add('A5:D39',true,'MatchEvidence');et.style='TableStyleMedium2';et.showFilterButton=true;
evidence.freezePanes.freezeRows(5);
evidence.freezePanes.freezeColumns(1);

for(const [sheet,range] of [[main,'A6:E6'],[main,'A18:M18'],[evidence,'A5:D5'],[evidence,'F5:J5']]){
  const rr=sheet.getRange(range);
  rr.format.fill=navy;rr.format.font={name:'Arial',size:11,bold:true,color:'#FFFFFF'};
  rr.format.horizontalAlignment='center';rr.format.verticalAlignment='center';rr.format.wrapText=true;
  rr.format.borders={insideVertical:{style:'thin',color:'#FFFFFF'}};
}
main.getRange('A6:E6').format.rowHeight=36;
for(const c of cols){
  const range=main.getRange(`${c}19:${c}52`);
  range.conditionalFormats.add('containsText',{text:'Different',format:{fill:'#FFF0DC',font:{color:'#81450B'}}});
  range.conditionalFormats.add('containsText',{text:'Absent',format:{fill:'#ECEFF2',font:{color:'#5B6470'}}});
}

wb.recalculate();
const totals=main.getRange('B7:E9').values;
const expected=[[33,28,5,1],[29,10,19,5],[23,8,15,11]];
if(JSON.stringify(totals)!==JSON.stringify(expected)) throw new Error(`Summary mismatch: ${JSON.stringify(totals)}`);
if(main.getRange('B14').values[0][0]!==8) throw new Error('All-four match total mismatch');
console.log((await wb.inspect({kind:'table',range:'Assembly comparison!A6:E9',include:'values,formulas',tableMaxRows:4,tableMaxCols:5,maxChars:2200})).ndjson);
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30},summary:'Formula error scan'})).ndjson);
for(const [sheet,range,name] of [
  ['Assembly comparison','A1:G9','preview_summary'],
  ['Assembly comparison','A18:G24','preview_toga'],
  ['Assembly comparison','H18:M24','preview_msa'],
  ['Evidence','A5:D10','preview_evidence'],
  ['Evidence','F5:J8','preview_resources']
]){
  const image=await wb.render({sheetName:sheet,range,scale:1,format:'png'});
  await fs.writeFile(`${dir}${name}.png`,new Uint8Array(await image.arrayBuffer()));
}
const output=await SpreadsheetFile.exportXlsx(wb);
await output.save(`${dir}Genome_assembly_comparison_Local_TOGA2_Multiz470_Cactus447.xlsx`);
console.log('Exported workbook with 34 species and verified summary counts.');
