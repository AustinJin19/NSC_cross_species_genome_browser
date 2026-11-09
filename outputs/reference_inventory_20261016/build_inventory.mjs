import fs from 'node:fs/promises';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const dir = '/Users/zjin27/Downloads/CS_Genome_Browser/outputs/reference_inventory_20261016';
const data = JSON.parse(await fs.readFile(`${dir}/inventory.json`, 'utf8'));
const wb = Workbook.create();
const refs = wb.worksheets.add('Species references');
const files = wb.worksheets.add('Local file evidence');
const navy='#003F7D', ink='#233247', light='#F1F5FA';
const col = n => {let s=''; for(n++;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s; return s;};
function styleSheet(sheet, headers, values, widths, tableName) {
  const bottom=6+values.length, last=col(headers.length-1);
  sheet.showGridLines=false;
  sheet.getRange(`A1:${last}${bottom}`).format.font={name:'Arial',size:11,color:ink};
  sheet.getRange(`A6:${last}${bottom}`).values=[headers,...values];
  const table=sheet.tables.add(`A6:${last}${bottom}`,true,tableName);
  table.style='TableStyleLight9';
  table.showFilterButton=true;
  sheet.getRange(`A6:${last}6`).format={fill:navy,font:{name:'Arial',size:11,bold:true,color:'#FFFFFF'},horizontalAlignment:'center',verticalAlignment:'center',wrapText:true,rowHeight:34};
  sheet.getRange(`A7:${last}${bottom}`).format.wrapText=true;
  sheet.getRange(`A7:${last}${bottom}`).format.verticalAlignment='top';
  widths.forEach((w,i)=>sheet.getRange(`${col(i)}1:${col(i)}${bottom}`).format.columnWidth=w);
  for(let i=0;i<values.length;i++){
    const lines=Math.max(...values[i].map((v,c)=>String(v??'').split('\n').reduce((n,x)=>n+Math.ceil(Math.max(1,x.length)/(widths[c]-4)),0)));
    sheet.getRange(`A${i+7}:${last}${i+7}`).format.rowHeight=Math.max(58,lines*15+12);
  }
  sheet.getRange('A2').format.font={name:'Arial',size:16,bold:true,color:navy};
  sheet.getRange('A2').format.rowHeight=25;
  sheet.getRange(`A3:${last}4`).format.font={name:'Arial',size:11,color:'#59687B'};
  sheet.getRange(`A3:${last}4`).format.wrapText=false;
  sheet.getRange(`A4:${last}4`).format.borders={bottom:{style:'thin',color:'#C9D4E3'}};
  sheet.freezePanes.freezeRows(6);
  sheet.freezePanes.freezeColumns(2);
  return bottom;
}
const head=['Common name','Scientific name','Browser status','Genome assembly used','Assembly accession','Gene annotation used','Version recorded locally','Assembly resource','Annotation resource','What the links provide','Evidence and limitations'];
const values=data.species.map(r=>[r.common,r.scientific,r.status,r.assembly,r.accession||'Unresolved',r.annotation,r.release,r.assembly_url?'Open assembly resource':'No verified URL',r.annotation_url?'Open annotation resource':'No verified URL',r.link_scope,r.notes]);
const end=styleSheet(refs,head,values,[25,29,23,35,24,46,43,26,28,59,77],'SpeciesReferences');
refs.tabColor=navy;
refs.getRange('A2').values=[['Genome assemblies and gene annotations']];
refs.getRange('A3').values=[['NSC cross-species ATAC project: 31 active browser species and 4 additional species without bigWig tracks.']];
refs.getRange('A4').values=[['Inventory: October 16, 2026. Resource links do not establish an exact file match unless stated. Local file evidence includes retained alternatives.']];
// The renderer does not implement HYPERLINK formulas. Keep labels as plain
// values here; add native OOXML hyperlink relationships after export.
refs.getRange(`H7:I${end}`).format.font={name:'Arial',size:11,color:'#0563C1',underline:'single'};
refs.getRange(`E7:E${end}`).conditionalFormats.add('containsText',{text:'Unresolved',format:{fill:'#FFF1CE',font:{color:'#805300',bold:true}}});
refs.getRange(`C7:C${end}`).conditionalFormats.add('containsText',{text:'Not active',format:{fill:'#F1F3F6',font:{color:'#667085'}}});

const fhead=['Scientific name','Use in this project','FASTA file','Gene annotation file','Active browser model index','FASTA sequences','Indexed FASTA length (bp)','First FASTA header','GTF header / provenance','Feature source example','Transcript example','GTF SHA-256'];
const fvals=data.files.map(r=>[r.scientific,r.status,r.fasta,r.gtf,r.index||'Not active',r.sequence_count,r.total_bp,r.header,r.annotation_header,r.feature_source,r.example_transcript,r.gtf_sha256]);
const fend=styleSheet(files,fhead,fvals,[29,43,68,72,65,16,22,90,90,22,52,72],'ReferenceFileEvidence');
files.tabColor='#8294AA';
files.getRange('A2').values=[['Local reference files and provenance']];
files.getRange('A3').values=[['44 FASTA–annotation pairs. Paths are relative to the CS_Genome_Browser project folder.']];
files.getRange('A4').values=[['Sequence counts and lengths come from .fai indexes. SHA-256 identifies the local compressed GTF. Alternatives are not assumed to match the active annotation.']];
files.getRange(`F7:G${fend}`).setNumberFormat('#,##0');
files.getRange(`F7:G${fend}`).format.horizontalAlignment='right';
files.getRange(`B7:B${fend}`).conditionalFormats.add('containsText',{text:'not active',format:{fill:'#F1F3F6',font:{color:'#667085'}}});
wb.recalculate();
for(const [sheet,range] of [['Species references','A6:G10'],['Species references','H7:I9'],['Local file evidence','A6:G9']]){
  console.log((await wb.inspect({kind:'table',range:`'${sheet}'!${range}`,include:'values,formulas',tableMaxRows:5,tableMaxCols:7,maxChars:3200})).ndjson);
}
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:50},summary:'Final formula error scan',maxChars:2000})).ndjson);
for(const [sheet,range,file] of [
 ['Species references','A1:F12','preview_references.png'],
 ['Species references','F6:K10','preview_resources.png'],
 ['Species references','A36:F41','preview_pending.png'],
 ['Local file evidence','A1:D12','preview_files.png'],
 ['Local file evidence','H6:L10','preview_provenance.png'],
]){
  const image=await wb.render({sheetName:sheet,range,scale:1.4,format:'png'});
  await fs.writeFile(`${dir}/${file}`,new Uint8Array(await image.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(`${dir}/NSC_genome_assemblies_and_annotations.xlsx`);
console.log('Saved NSC_genome_assemblies_and_annotations.xlsx');
