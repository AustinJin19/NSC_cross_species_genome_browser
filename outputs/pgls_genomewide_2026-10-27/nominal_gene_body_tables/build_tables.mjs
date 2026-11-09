import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import assert from 'node:assert/strict';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const out = path.dirname(fileURLToPath(import.meta.url));
const base = path.dirname(out);
const text = await fs.readFile(path.join(base, 'results/body/nominal_p05.tsv'), 'utf8');
const lines = text.trimEnd().split(/\r?\n/);
const keys = lines.shift().split('\t');
const records = lines.map(line => Object.fromEntries(line.split('\t').map((v,i) => [keys[i],v])));
assert.equal(records.length, 1949);
assert.equal(new Set(records.map(r => r.gene)).size, 1949);
for (const r of records) {
  assert.equal(r.method, 'body');
  assert.ok(Number(r.p_value) < 0.05);
  assert.ok(Number(r.BH_FDR) >= 0.05);
  assert.equal(Number(r.MLS_beta) > 0 ? 'Positive' : 'Negative', r.direction);
}

const columns = [
  ['Gene','gene','text',18],
  ['Direction','direction','text',14],
  ['Species (n)','n_species','0',13],
  ['MLS coefficient','MLS_beta','0.0000',17],
  ['Coefficient SE','MLS_SE','0.0000',17],
  ['Nominal p-value','p_value','0.00E+00',18],
  ['BH FDR','BH_FDR','0.0000',14],
  ['Estimated lambda','lambda','0.000',18],
  ['Residual df','df','0',13],
  ['GLS model R²','GLS_R2','0.0000',16],
  ['Body-mass coefficient','mass_beta','0.0000',20],
  ['Gene-length coefficient','length_beta','0.0000',21],
  ['OLS MLS coefficient','OLS_beta','0.0000',20],
  ['OLS p-value','OLS_p_value','0.00E+00',16],
  ['OLS BH FDR','OLS_BH_FDR','0.0000',16],
];
const wb = Workbook.create();
const verification = {};
for (const [direction, expected, accent] of [['Positive',569,'#164E73'],['Negative',1380,'#863D32']]) {
  const rows = records.filter(r => r.direction === direction).sort((a,b) => Number(a.p_value)-Number(b.p_value) || a.gene.localeCompare(b.gene));
  assert.equal(rows.length, expected);
  const sh = wb.worksheets.add(direction);
  sh.showGridLines = false;
  sh.tabColor = accent;
  const last = rows.length + 9;
  sh.getRange(`A1:O${last}`).format.font = {name:'Arial',size:10,color:'#203040'};
  sh.getRange(`A1:O${last}`).format.verticalAlignment = 'center';
  sh.getRange(`A1:O${last}`).format.rowHeight = 21;
  sh.getRange('A2').values = [[`Gene-body PGLS: ${direction.toLowerCase()} MLS associations`]];
  sh.getRange('A2').format.font = {name:'Arial',size:14,bold:true,color:accent};
  sh.getRange('A2:H2').format.rowHeight = 28;
  sh.getRange('A2:H2').format.borders = {bottom:{style:'thin',color:accent}};
  const notes = [
    `${expected.toLocaleString('en-US')} genes at nominal p < 0.05. None pass BH FDR < 0.05.`,
    'log10(activity) ~ log10(MLS) + log10(body mass) + log10(gene length)',
    "Pagel's lambda estimated per gene. Human, macaque and bowhead whale excluded.",
    'Source: genome-wide PGLS scan, 27 October 2026. Gene body + 2 kb promoter.',
    'Activity: integrated CPM·bp. BH FDR calculated across 20,200 tested genes.'
  ];
  sh.getRange('A3:A7').values = notes.map(v=>[v]);
  sh.getRange('A3:A7').format.font = {name:'Arial',size:10,color:'#4B5868'};
  sh.getRange('A3').format.font = {name:'Arial',size:10,bold:true,color:accent};
  sh.getRange('A9:O9').values = [columns.map(c=>c[0])];
  const data = rows.map(r=>columns.map(c=>c[2]==='text' ? r[c[1]] : Number(r[c[1]])));
  sh.getRange(`A10:O${last}`).values = data;
  const table = sh.tables.add(`A9:O${last}`,true,`${direction}GeneResults`);
  table.showFilterButton = true;
  sh.getRange('A9:O9').format = {
    fill:accent,font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},
    horizontalAlignment:'center',verticalAlignment:'center',wrapText:true,
    rowHeight:34,
    borders:{insideVertical:{style:'thin',color:'#FFFFFF'}}
  };
  for (let i=0;i<columns.length;i++) {
    const c=String.fromCharCode(65+i);
    sh.getRange(`${c}1:${c}${last}`).format.columnWidth=columns[i][3];
    sh.getRange(`${c}10:${c}${last}`).setNumberFormat(columns[i][2]==='text'?'@':columns[i][2]);
    sh.getRange(`${c}10:${c}${last}`).format.horizontalAlignment=columns[i][2]==='text'?'left':'right';
  }
  sh.freezePanes.freezeRows(9);
  sh.freezePanes.freezeColumns(2);
  assert.deepEqual(sh.getRange(`A10:O${last}`).values,data);
  verification[direction]={count:rows.length,first:rows[0].gene,last:rows.at(-1).gene,minP:Number(rows[0].p_value),maxP:Number(rows.at(-1).p_value)};
}
wb.recalculate();
for (const name of ['Positive','Negative']) {
  const check=await wb.inspect({kind:'table',range:`${name}!A9:H12`,include:'values',tableMaxRows:4,tableMaxCols:8,maxChars:1800});
  console.log(check.ndjson);
  const preview=await wb.render({sheetName:name,range:'A1:H19',scale:1.5,format:'png'});
  await fs.writeFile(path.join(out,`${name.toLowerCase()}_preview.png`),new Uint8Array(await preview.arrayBuffer()));
}
const err=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},maxChars:1000});
console.log(err.ndjson);
const output=await SpreadsheetFile.exportXlsx(wb);
const destination=path.join(out,'Gene_body_PGLS_nominal_569_positive_1380_negative.xlsx');
await output.save(destination);
await fs.writeFile(path.join(out,'verification.json'),JSON.stringify(verification,null,2));
console.log(JSON.stringify({destination,verification}));
