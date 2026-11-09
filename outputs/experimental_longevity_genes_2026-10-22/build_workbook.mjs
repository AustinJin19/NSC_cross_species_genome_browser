import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const root = path.dirname(fileURLToPath(import.meta.url));
const data = JSON.parse(await fs.readFile(path.join(root, 'collection.json'), 'utf8'));
const orthology = JSON.parse(await fs.readFile(path.join(root, 'orthology.json'), 'utf8'));
const wb = Workbook.create();
const overview = wb.worksheets.add('Overview');
const mappings = wb.worksheets.add('Ortholog mappings');
const mappingEvidence = wb.worksheets.add('Ortholog evidence');
const evidence = wb.worksheets.add('Experiments');
const screen = wb.worksheets.add('GenAge screen');
const navy = '#173B55', blue = '#EAF2F8', teal = '#E3F3EB', amber = '#FFF2D8';
const baseFont = {name: 'Aptos', size: 11, color: '#203547'};

function col(n) {
  let s = '';
  for (n++; n > 0; n = Math.floor((n - 1) / 26)) s = String.fromCharCode(65 + (n - 1) % 26) + s;
  return s;
}
function tableSheet(sheet, records, tableName, widths, height) {
  const keys = Object.keys(records[0]);
  const last = col(keys.length - 1);
  const end = records.length + 1;
  const values = [keys, ...records.map(r => keys.map(k => r[k]))];
  sheet.getRange(`A1:${last}${end}`).values = values;
  sheet.getRange(`A1:${last}${end}`).format = {
    font: baseFont, wrapText: true, verticalAlignment: 'top', rowHeight: height,
  };
  sheet.tables.add(`A1:${last}${end}`, true, tableName);
  sheet.getRange(`A1:${last}1`).format = {
    fill: navy, font: {...baseFont, bold: true, color: '#FFFFFF'}, rowHeight: 48,
    verticalAlignment: 'center', wrapText: true,
  };
  widths.forEach((width, i) => sheet.getRange(`${col(i)}1:${col(i)}${end}`).format.columnWidth = width);
  for (let i = 2; i <= end; i++) {
    sheet.getRange(`A${i}:${last}${i}`).format.fill = i % 2 === 0 ? '#F3F7FA' : '#FFFFFF';
  }
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(2);
  sheet.showGridLines = false;
  return {keys, last, end};
}

const mp = tableSheet(mappings, orthology.summary, 'OrthologSummary',
  [15, 13, 40, 43, 53, 58, 24, 20, 79, 68, 58, 18, 55, 65, 65], 98);
mappings.tabColor = '#216FA0';
for (let i = 0; i < orthology.summary.length; i++) {
  const record = orthology.summary[i];
  mappings.getRange(`C${i+2}:D${i+2}`).format.fill = teal;
  mappings.getRange(`A${i+2}`).format.font = {...baseFont, bold: true};
  if (record['Source gene'] === 'lipl-4') {
    mappings.getRange(`A${i+2}:${mp.last}${i+2}`).format.rowHeight = 156;
  }
  if (record['Source gene'] === 'pat-10') {
    mappings.getRange(`D${i+2}`).format.fill = amber;
    mappings.getRange(`I${i+2}`).format.fill = amber;
  }
}
const me = tableSheet(mappingEvidence, orthology.pairs, 'OrthologSupport',
  [15, 29, 19, 25, 17, 18, 21, 19, 25, 15, 18, 20, 18, 18, 43, 72, 27, 20, 18, 18, 65, 55], 58);
mappingEvidence.tabColor = '#216FA0';
for (let i = 0; i < orthology.pairs.length; i++) {
  mappingEvidence.getRange(`O${i+2}`).format.fill =
    orthology.pairs[i]['Mapping decision'] === 'Retained' ? teal : amber;
}

const ev = tableSheet(evidence, data.experiments, 'PrimaryExperiments',
  [12, 19, 14, 25, 24, 46, 42, 66, 54, 34, 72, 14, 78, 60, 60], 68);
evidence.tabColor = '#14856C';
evidence.getRange(`L2:L${ev.end}`).setNumberFormat('0');
for (let i = 0; i < data.experiments.length; i++) {
  const d = data.experiments[i]['Beneficial perturbation direction'];
  evidence.getRange(`E${i+2}`).format.fill = d === 'Decrease' ? blue : d === 'Increase' ? teal : amber;
  evidence.getRange(`E${i+2}`).format.font = {...baseFont, bold: true};
}

const sc = tableSheet(screen, data.screen, 'GenAgeScreening',
  [13, 22, 29, 60, 19, 29, 30, 66, 33, 60, 53], 42);
screen.tabColor = '#BF8A30';
for (let i = 0; i < data.screen.length; i++) {
  if (data.screen[i]['Lifespan effect (GenAge)'] === 'Increase and Decrease') {
    screen.getRange(`F${i+2}`).format.fill = amber;
  }
}

// Plain, copyable URL values: the runtime does not calculate HYPERLINK.
for (const [sheet, records, mapping] of [
  [evidence, data.experiments, [['N', 'Primary source URL'], ['O', 'Additional primary source URL']]],
  [screen, data.screen, [['J', 'Source URL'], ['K', 'NCBI Gene URL']]],
  [mappings, orthology.summary, [['M', 'Source NCBI Gene URL'], ['N', 'Human mapping URL'], ['O', 'Mouse mapping URL']]],
  [mappingEvidence, orthology.pairs, [['U', 'DIOPT API URL'], ['V', 'Target NCBI Gene URL']]],
]) {
  for (const [letter, field] of mapping) {
    sheet.getRange(`${letter}2:${letter}${records.length+1}`).format.font = {...baseFont, color: '#176A99'};
  }
}

overview.showGridLines = false;
overview.tabColor = navy;
overview.getRange('A1:C39').format = {font: baseFont, wrapText: true, verticalAlignment: 'center', rowHeight: 28};
overview.getRange('A1:A39').format.columnWidth = 33;
overview.getRange('B1:B39').format.columnWidth = 107;
overview.getRange('C1:C39').format.columnWidth = 22;
overview.getRange('A1:C1').values = [['ANIMAL LONGEVITY', 'Experimental gene perturbations', '22 October 2026']];
overview.getRange('A1:C1').format = {fill: navy, font: {...baseFont, size: 17, bold: true, color: '#FFFFFF'}, rowHeight: 52, wrapText: true};
overview.getRange('A2:C2').values = [['Scope', 'Selected primary studies and a separately labelled GenAge screening set. Not an exhaustive or systematic replication review.', '']];
overview.getRange('A2:C2').format = {fill: blue, rowHeight: 46, wrapText: true};
overview.getRange('A4:C4').values = [['COLLECTION', 'What is counted', 'Number']];
overview.getRange('A4:C4').format = {fill: navy, font: {...baseFont, color: '#FFFFFF', bold: true}};
overview.getRange('A5:B8').values = [
  ['Primary studies', 'Intervention records'],
  ['Primary studies', 'Distinct gene-species pairs'],
  ['GenAge screen', 'Animal gene records with a reported increase'],
  ['GenAge screen', 'Of these: mixed increase and decrease outcomes'],
];
overview.getRange('C5').formulas = [[`=COUNTA(Experiments!A2:A${ev.end})`]];
overview.getRange('C6').values = [[data.stats.unique_gene_species_pairs]];
overview.getRange('C7').formulas = [[`=COUNTA('GenAge screen'!A2:A${sc.end})`]];
overview.getRange('C8').formulas = [[`=COUNTIF('GenAge screen'!F2:F${sc.end},"Increase and Decrease")`]];
overview.getRange('A5:C6').format.fill = teal;
overview.getRange('A7:C8').format.fill = amber;
overview.getRange('C5:C8').format.font = {...baseFont, bold: true, size: 15};
overview.getRange('A10:C10').values = [['Animal', 'Distinct genes within this species', 'Interventions']];
overview.getRange('A10:C10').format = {fill: navy, font: {...baseFont, bold: true, color: '#FFFFFF'}};
const order = ['Mus musculus', 'Caenorhabditis elegans', 'Drosophila melanogaster', 'Anopheles stephensi'];
const names = ['Mouse', 'Worm', 'Fruit fly', 'Mosquito'];
order.forEach((s, i) => {
  const row = 11 + i;
  overview.getRange(`A${row}:C${row}`).values = [[names[i], `${data.stats.unique_genes_by_species[s]} | ${s}`, null]];
  overview.getRange(`A${row}:C${row}`).format.fill = i % 2 === 0 ? '#F3F7FA' : '#FFFFFF';
  overview.getRange(`C${row}`).formulas = [[`=COUNTIF(Experiments!D2:D${ev.end},"${s}")`]];
});
overview.getRange('A16:C16').values = [['INTERPRETATION', 'How to use the evidence', '']];
overview.getRange('A16:C16').format = {fill: navy, font: {...baseFont, bold: true, color: '#FFFFFF'}};
const notes = [
  ['Direct lifespan evidence', 'Experiments lists knockout, partial loss, RNAi, activating alleles, overexpression, gene delivery or direct targeted inhibition with a reported animal lifespan increase. GWAS, expression correlations, healthspan alone and progeria rescue alone are insufficient.'],
  ['Direction matters', 'Decrease / Increase describes the target gene or protein, not its downstream pathway. Activate / alter function retains allele or isoform details. Full knockout is not interchangeable with partial loss.'],
  ['Context matters', 'Keep tissue, sex, age, genetic background and diet with each gene. Some benefits occur only in particular conditions. Effect sizes are study-specific and are not pooled estimates or a ranking of interventions.'],
  ['Scope of verification', 'Primary reports were checked for the intervention and lifespan outcome. Numerical effect sizes, sample sizes, sex-stratified results and replication history were not exhaustively extracted. Qualitative entries do not mean the paper lacks numerical data.'],
  ['ATAC interpretation', 'Accessibility is not identical to RNA abundance, protein activity or experimental gene perturbation. These experiments do not by themselves predict the sign of your cross-species liver ATAC-MLS correlation.'],
  ['Cross-species matching', 'Ortholog mappings provides human and mouse candidates for the 17 curated worm/fly genes. Ortholog evidence retains support scores and exclusions. Original animal lifespan evidence stays in Experiments; homology does not establish lifespan extension in mammals. Has2 was a naked-mole-rat transgene tested in mice.'],
  ['GenAge screening', 'GenAge Build 21, released 28 August 2023, downloaded 22 October 2026. Retained animal records whose lifespan-effect field contains Increase, including mixed outcomes. Fungal records excluded.'],
  ['Screening limitations', 'The 852 database records have not each been checked against their underlying experiments. Original longevity-influence labels are retained; they are not exact KO/OE instructions. The 54 mixed-outcome records need particular care.'],
  ['Source and provenance', 'GenAge: https://genomics.senescence.info/genes/models.html . Primary paper links appear on each Experiments row. The original source archive, release notes, TSV exports and checksum manifest are preserved beside this workbook.'],
];
notes.forEach(([label, text], i) => {
  const row = 17 + i;
  overview.getRange(`A${row}:C${row}`).values = [[label, text, '']];
  overview.getRange(`A${row}:C${row}`).format.rowHeight = 74;
  overview.getRange(`A${row}`).format = {fill: blue, font: {...baseFont, bold: true}, wrapText: true};
});
overview.getRange('A27:C27').values = [['ORTHOLOG MAPPING', 'Curated worm and fly genes only', 'Number']];
overview.getRange('A27:C27').format = {fill: navy, font: {...baseFont, bold: true, color: '#FFFFFF'}};
overview.getRange('A28:C32').values = [
  ['Source genes', '11 worm genes and 6 fruit-fly genes', orthology.stats.source_genes],
  ['Human candidates', 'Distinct retained human genes; all 17 source genes mapped', orthology.stats.distinct_targets.Human],
  ['Mouse candidates', 'Distinct retained mouse genes; 16 source genes mapped', orthology.stats.distinct_targets.Mouse],
  ['Retained pairs', 'Source gene to target gene links across both target species', orthology.stats.retained_pairs],
  ['Unresolved mouse mapping', 'pat-10: strongest mouse predictions score only 2/13', 1],
];
overview.getRange('A28:C31').format.fill = teal;
overview.getRange('A32:C32').format.fill = amber;
overview.getRange('C28:C32').format.font = {...baseFont, bold: true, size: 15};
const mappingNotes = [
  ['Mapping rule', 'DIOPT score >=4, high/moderate confidence, and a current protein-coding NCBI record. This is our screening cutoff, not a probability or proof of function. Multiple supported candidates are retained rather than forcing one gene.'],
  ['Dataset used', 'DIOPT API v9 returned dataset 9.1 (13 methods maximum), retrieved 22 October 2026. The advertised v10 endpoint was unavailable. NCBI gene type and current symbols were checked; this is not a claim to use the latest DIOPT release.'],
  ['Ambiguous families', 'daf-2 maps to IGF1R/INSR/INSRR; lipl-4 maps to several lipases. pat-10 maps to human CALML6, but its mouse mapping is unresolved at our cutoff. FOXO3B is retained under its current protein-coding NCBI/HGNC annotation.'],
  ['Evidence boundary', 'Lifespan effects and beneficial perturbation directions were measured in the source worm or fly. They must not be assigned automatically to every human/mouse candidate. Gene family membership is not evidence of a conserved intervention effect.'],
  ['Mapping sources', 'DIOPT documentation: https://fgr.hms.harvard.edu/diopt-documentation . Per-gene DIOPT API and NCBI links are in the mapping sheets. Raw responses and 550 predictions are preserved in the companion audit export.'],
];
mappingNotes.forEach(([label, text], i) => {
  const row = 34 + i;
  overview.getRange(`A${row}:C${row}`).values = [[label, text, '']];
  overview.getRange(`A${row}:C${row}`).format.rowHeight = 74;
  overview.getRange(`A${row}`).format = {fill: blue, font: {...baseFont, bold: true}, wrapText: true};
});
overview.freezePanes.freezeRows(2);
wb.recalculate();

const report = await wb.inspect({kind: 'region', sheetId: 'Overview', range: 'A4:C14', maxChars: 3500, tableMaxRows: 11, tableMaxCols: 3});
await fs.writeFile(path.join(root, 'ortholog_workbook_inspection.txt'), report.ndjson ?? JSON.stringify(report));
for (const [sheetName, range, filename] of [
  ['Overview', 'A27:C38', 'preview_ortholog_overview.png'],
  ['Ortholog mappings', 'A1:F9', 'preview_ortholog_mappings.png'],
  ['Ortholog mappings', 'G1:K5', 'preview_ortholog_notes.png'],
  ['Ortholog evidence', 'E1:O7', 'preview_ortholog_evidence.png'],
]) {
  const p = await wb.render({sheetName, range, scale: 1.25, format: 'png'});
  await fs.writeFile(path.join(root, filename), new Uint8Array(await p.arrayBuffer()));
}
const output = await SpreadsheetFile.exportXlsx(wb);
const outputPath = path.join(root, 'Experimental_animal_lifespan_genes_with_orthologs.xlsx');
await output.save(outputPath);
console.log(JSON.stringify({output: outputPath, records: data.experiments.length, screening: data.screen.length, orthologSources: orthology.summary.length, orthologPairs: orthology.stats.retained_pairs}));
