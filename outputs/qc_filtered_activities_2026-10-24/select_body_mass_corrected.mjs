import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { fileURLToPath } from 'node:url';

const root = path.dirname(fileURLToPath(import.meta.url));
const output = path.join(root, 'body_mass_corrected');
const web = path.resolve(root, '../../cross_species_visual_tool/web_app/data');
const sha = filename => crypto.createHash('sha256').update(fs.readFileSync(filename)).digest('hex');
const parse = filename => {
  const lines = fs.readFileSync(filename, 'utf8').trimEnd().split(/\r?\n/);
  const fields = lines[0].split('\t');
  return { fields, rows: lines.slice(1).map(line => Object.fromEntries(line.split('\t').map((value, i) => [fields[i], value]))) };
};
const write = (filename, fields, rows) => fs.writeFileSync(filename,
  [fields.join('\t'), ...rows.map(row => fields.map(field => row[field]).join('\t'))].join('\n') + '\n');

fs.mkdirSync(output, { recursive: true });
let summaries = [], summaryFields;
for (const [family, manifestName] of [['gene_activity', 'activity'], ['regulatory_activity', 'regulatory']]) {
  const manifest = JSON.parse(fs.readFileSync(path.join(root, `${manifestName}_manifest.json`), 'utf8'));
  if (manifest.qc_sha256 !== sha(path.join(web, 'sample_qc_exclusions.json')) ||
      manifest.traits_sha256 !== sha(path.join(web, 'species_mls.json'))) throw Error('QC/traits changed since recalculation');
  const summary = parse(path.join(root, `${family}_summary.tsv`));
  const significant = parse(path.join(root, `${family}_significant_genes.tsv`));
  const selectedSummary = summary.rows.filter(row => row.body_mass_correction === 'True');
  const selectedGenes = significant.rows.filter(row => row.body_mass_correction === 'True');
  for (const row of selectedSummary) {
    const model = selectedGenes.filter(gene => gene.method === row.method && gene.gene_length_correction === row.gene_length_correction);
    for (const direction of ['Positive', 'Negative']) {
      if (model.filter(gene => gene.direction === direction).length !== Number(row[`FDR_${direction.toLowerCase()}`])) throw Error('Count mismatch');
    }
    if (model.some(gene => Number(gene.BH_FDR) >= 0.05 || Number(gene.n_species) < 5)) throw Error('Invalid significance filter');
  }
  write(path.join(output, `${family}_significant_genes.tsv`), significant.fields, selectedGenes);
  summaries.push(...selectedSummary); summaryFields = summary.fields;
  console.log(JSON.stringify({family, models: selectedSummary.length, significant_gene_window_rows: selectedGenes.length,
    unique_genes: new Set(selectedGenes.map(row => row.gene)).size,
    positive_unique: new Set(selectedGenes.filter(row => row.direction === 'Positive').map(row => row.gene)).size,
    negative_unique: new Set(selectedGenes.filter(row => row.direction === 'Negative').map(row => row.gene)).size}));
}
write(path.join(output, 'summary.tsv'), summaryFields, summaries);
fs.writeFileSync(path.join(output, 'METHODS.txt'), [
  'Body-mass-corrected results from QC-filtered ATAC samples, 2026-10-24.',
  'Selected from the completed recalculation; no additional model or threshold changes.',
  'Response: log10(species-median activity). Predictor: log10(MLS in years). Covariate: log10(adult body mass in grams).',
  'Gene-body models with and without log10(annotated gene length) are both listed and explicitly flagged.',
  'Regulatory activity here means overlapping peak count, not CPM signal within peaks.',
  'Exclude samples with TSS enrichment <4; retain missing TSS. Human and macaque are excluded from fits; no bowhead samples pass QC.',
  'Highest-median locus selected after sample exclusion. Require >=50% genomic region availability and >=5 species.',
  'Significance: BH FDR <0.05, calculated separately per method and covariate model. Direction is the partial correlation with MLS.',
  'Each gene list has one row per significant gene/method/model; the same gene can appear in multiple windows.',
  'Source: ../activity_manifest.json, ../regulatory_manifest.json and the corresponding *_significant_genes.tsv / *_summary.tsv files.'
].join('\n') + '\n');
