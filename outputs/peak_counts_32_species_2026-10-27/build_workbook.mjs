import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const out = path.dirname(fileURLToPath(import.meta.url));
const data = JSON.parse(await fs.readFile(path.join(out, 'report_data.json'), 'utf8'));
const wb = Workbook.create();
const summary = wb.worksheets.add('Species summary');
const samples = wb.worksheets.add('Sample counts');
const navy = '#173B50', light = '#F0F5F8';
function base(sheet, range) {
  sheet.showGridLines = false;
  sheet.getRange(range).format.font = {name:'Arial',size:11,color:'#203744'};
  sheet.getRange(range).format.verticalAlignment = 'center';
  sheet.getRange(range).format.rowHeight = 23;
}
function headers(sheet, range) {
  sheet.getRange(range).format = {fill:navy,font:{name:'Arial',size:11,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center',verticalAlignment:'center',rowHeight:42};
}
base(summary,'A1:J43'); base(samples,'A1:L86');
summary.tabColor = navy;
summary.getRange('A2').values = [['ATAC peak counts by species']];
summary.getRange('A2').format.font = {name:'Arial',size:16,bold:true,color:navy};
summary.getRange('A3').values = [['32 species; 80 available samples. Human, rhesus macaque and bowhead whale excluded. No TSS/QC exclusion.']];
summary.getRange('A4').values = [['Peak entries count narrowPeak records, including multiple summits in the same interval.']];
summary.getRange('A5').values = [['Distinct regions count unique chromosome/start/end coordinates per sample. Medians and ranges are across samples.']];
summary.getRange('A3:J5').format.font = {name:'Arial',size:11,color:'#5F6C76'};
const sh = ['Common name','Scientific name','Samples counted','Median peak entries','Min peak entries','Max peak entries','Median distinct regions','Min distinct regions','Max distinct regions','Missing sample'];
summary.getRange('A7:J7').values = [sh];
headers(summary,'A7:J7');
summary.getRange('A8:J39').values = data.summary.map(r => [r.common_name,r.scientific_name,r.samples_counted,r.median_peaks,r.min_peaks,r.max_peaks,r.median_unique_regions,r.min_unique_regions,r.max_unique_regions,data.samples.filter(s=>s.species===r.species&&s.peak_count===null).map(s=>s.sample).join(', ')]);
summary.getRange('C8:C39').setNumberFormat('#,##0');
summary.getRange('D8:I39').setNumberFormat('#,##0.0');
summary.getRange('E8:F39').setNumberFormat('#,##0');
summary.getRange('H8:I39').setNumberFormat('#,##0');
summary.getRange('C8:I39').format.horizontalAlignment='right';
summary.getRange('B8:B39').format.font.italic=true;
for (let r=8;r<=39;r+=2) summary.getRange(`A${r}:J${r}`).format.fill=light;
const widths=[245,240,90,115,105,105,130,115,115,120];
widths.forEach((v,i)=>summary.getRangeByIndexes(0,i,43,1).format.columnWidthPx=v);
summary.freezePanes.freezeRows(7);

samples.getRange('A2').values=[['Per-sample peak counts']];
samples.getRange('A2').format.font={name:'Arial',size:16,bold:true,color:navy};
samples.getRange('A3').values=[['Sources: canonical narrowPeak files listed in manifest.local.tsv; QC metadata from bundle_sample_qc_summary.tsv. Counted 27 October 2026.']];
samples.getRange('A3:L3').format.font={name:'Arial',size:11,color:'#5F6C76'};
const ph=['Common name','Scientific name','Sample','Peak entries','Distinct regions','Extra entries at same interval','Peak-file status','TSS enrichment','Final read pairs','Assembly label','','Source peak file'];
samples.getRange('A5:L5').values=[ph];
samples.getRange('A6:L86').values=data.samples.map(r=>[r.common_name,r.scientific_name,r.sample,r.peak_count,r.unique_peak_regions,r.peak_count===null?null:r.duplicate_intervals,r.status,r.tss_enrichment,r.final_read_pairs,r.assembly,null,r.peak_file==='-'?'Not supplied':path.basename(r.peak_file)]);
headers(samples,'A5:J5'); headers(samples,'L5:L5');
samples.getRange('D6:F86').setNumberFormat('#,##0');
samples.getRange('H6:H86').setNumberFormat('0.00');
samples.getRange('I6:I86').setNumberFormat('#,##0');
samples.getRange('D6:F86').format.horizontalAlignment='right';
samples.getRange('H6:I86').format.horizontalAlignment='right';
for(let r=6;r<=86;r+=2) samples.getRange(`A${r}:J${r}`).format.fill=light;
const pw=[245,240,155,110,120,135,145,105,130,135,22,360];
pw.forEach((v,i)=>samples.getRangeByIndexes(0,i,86,1).format.columnWidthPx=v);
samples.getRange('B6:B86').format.font.italic=true;
samples.getRange('G6:G86').conditionalFormats.add('containsText',{text:'Missing',format:{fill:'#FFF1D6',font:{color:'#8C4B00',bold:true}}});
samples.freezePanes.freezeRows(5);

wb.recalculate();
const inspected = await wb.inspect({kind:'table',range:"'Species summary'!A7:J12",include:'values',tableMaxRows:6,tableMaxCols:10,maxChars:2800});
await fs.writeFile(path.join(out,'workbook_inspection.ndjson'),inspected.ndjson);
const err=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},summary:'Error scan'});
await fs.writeFile(path.join(out,'workbook_errors.ndjson'),err.ndjson);
for(const [sheetName,range,file] of [
  ['Species summary','A1:J39','species_summary_preview.png'],
  ['Sample counts','A1:L15','sample_counts_preview.png']
]) {
  const preview=await wb.render({sheetName,range,scale:1.3,format:'png'});
  await fs.writeFile(path.join(out,file),new Uint8Array(await preview.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(path.join(out,'ATAC_peak_counts_32_species.xlsx'));
console.log(JSON.stringify({output:'ATAC_peak_counts_32_species.xlsx',species:data.summary.length,samples:data.samples.length,available_samples:data.provenance.samples_counted}));
