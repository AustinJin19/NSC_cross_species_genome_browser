import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const root=path.dirname(fileURLToPath(import.meta.url));
const d=JSON.parse(await fs.readFile(path.join(root,'reference.json'),'utf8'));
const wb=Workbook.create();
const overview=wb.worksheets.add('Overview');
const genes=wb.worksheets.add('Genes');
const evidence=wb.worksheets.add('Evidence');
const experiments=wb.worksheets.add('Experiments');
const models=wb.worksheets.add('GenAge models');
const gwas=wb.worksheets.add('GWAS associations');
const traits=wb.worksheets.add('GWAS trait selection');
const methods=wb.worksheets.add('Methods');
const navy='#173B55', pale='#EEF4F8', green='#E7F3ED', amber='#FFF2D9';
const font={name:'Aptos',size:11,color:'#23394B'};
const col=n=>{let s='';for(n++;n>0;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};

function table(sheet,records,name,widths) {
  const keys=Object.keys(records[0]), end=records.length+1,last=col(keys.length-1);
  sheet.getRange(`A1:${last}${end}`).values=[keys,...records.map(r=>keys.map(k=>r[k]??null))];
  sheet.getRange(`A1:${last}${end}`).format={font,wrapText:true,verticalAlignment:'center'};
  sheet.tables.add(`A1:${last}${end}`,true,name);
  widths.forEach((w,i)=>sheet.getRange(`${col(i)}1:${col(i)}${end}`).format.columnWidth=w);
  sheet.getRange(`A1:${last}1`).format={fill:navy,font:{...font,bold:true,color:'#FFFFFF'},rowHeight:52,wrapText:true};
  records.forEach((r,i)=>{
    const lines=Math.max(...keys.map((k,j)=>String(r[k]??'').split('\n').reduce((n,t)=>n+Math.max(1,Math.ceil(t.length/(widths[j]*0.84))),0)));
    sheet.getRange(`A${i+2}:${last}${i+2}`).format={fill:i%2?pale:'#FFFFFF',rowHeight:Math.max(28,lines*16+12)};
  });
  sheet.showGridLines=false; sheet.freezePanes.freezeRows(1); sheet.freezePanes.freezeColumns(1);
  sheet.tabColor=navy;
  return {keys,end,last};
}

const gr=table(genes,d.master,'HumanGeneReference',[17,39,15,14,14,14,15,61,19,17,28,17,21,21,21,21,19,19,20,20,21,34,27,39,16]);
genes.getRange(`C2:F${gr.end}`).format.horizontalAlignment='center';
genes.getRange(`G2:G${gr.end}`).setNumberFormat('0');
genes.getRange(`G2:G${gr.end}`).format.horizontalAlignment='center';
genes.getRange(`U2:U${gr.end}`).setNumberFormat('0.00E+00');
genes.getRange(`Y2:Y${gr.end}`).setNumberFormat('0');
d.master.forEach((r,i)=>{
  ['Experiment','GenAge','CellAge','GWAS'].forEach((k,j)=>{
    if(r[k]==='Yes') genes.getRange(`${col(j+2)}${i+2}`).format={fill:green,font:{...font,bold:true}};
  });
});
const ev=table(evidence,d.evidence,'EvidenceRecords',[17,16,30,17,28,20,54,42,43,104,17,21,27,20,73]);
evidence.getRange(`K2:K${ev.end}`).setNumberFormat('0.00E+00');
evidence.getRange(`O2:O${ev.end}`).format.font={...font,color:'#176A99'};
const ex=table(experiments,d.experiments,'AnimalExperiments',[54,31,16,23,16,29,27,49,49,78,68,37,80,16,88,73,73]);
const mo=table(models,d.models,'GenAgeOriginalModels',[16,22,30,18,65,29,28,72,31,68]);
const gw=table(gwas,d.gwas,'GWASAssociations',[16,33,65,24,17,70,54,54,25,30,24,65,39,22,19,95,22,95,90,72]);
gwas.getRange(`E2:E${gw.end}`).setNumberFormat('0.00E+00');
gwas.getRange(`D2:E${gw.end}`).format.horizontalAlignment='center';
d.gwas.forEach((r,i)=>{if(r['Included at P <= 5e-8']==='Yes')gwas.getRange(`D${i+2}`).format.fill=green;});
table(traits,d.traits,'TraitSelection',[69,68,34,106,21]);

const notes=[
  {Topic:'Scope',Method:'A broad source-labelled human gene reference for aging and longevity. The human list is a union of curated animal lifespan experiments, GenAge human candidates and model homologs, CellAge senescence manipulation/expression genes, and selected human GWAS associations. It is not a list of proven lifespan-extending interventions in humans.',Resource:'Compilation',Version:'22 October 2026',URL:''},
  {Topic:'Source labels',Method:'Experiment means the previously curated primary animal lifespan-extension studies, not every experimental record in a database. GenAge and CellAge retain their own source labels even when based on experiments. A gene can have several labels; counts overlap, and shared papers are not independent replication.',Resource:'Compilation',Version:'This collection',URL:''},
  {Topic:'Human names',Method:'Current HGNC-approved symbols are resolved preferentially by human Entrez ID, then approved symbol or Ensembl ID, then an unambiguous HGNC previous/alias symbol. Noncoding genes and pseudogenes are retained and labelled by gene type. Unresolved identifiers are preserved in the TSV audit, not guessed.',Resource:'HGNC',Version:'Downloaded 22 October 2026',URL:'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'},
  {Topic:'Experimental collection',Method:`${d.stats.experiment_records} intervention records involving ${d.stats.experiment_gene_species_pairs} original gene–species pairs. Human candidates derive from 22 mouse genes plus 17 worm/fly genes. The mosquito PTEN record is preserved without a forced human mapping. Intervention, tissue, direction and primary paper remain on Experiments. This curated subset is not exhaustive.`,Resource:'Experiment',Version:'Prior collection, 22 October 2026',URL:'Per-record primary paper URLs on Experiments'},
  {Topic:'Experimental orthologs',Method:'Mouse genes use HCOP mappings supported by at least 3 databases. The prior worm/fly DIOPT 9.1 candidates are retained at score >=4, high/moderate confidence and current protein-coding status. One-to-many candidates remain explicit. The intervention benefit was tested in the source animal; it is not transferred to every human homolog.',Resource:'HCOP / DIOPT',Version:'HCOP files updated 1 October 2026; DIOPT dataset 9.1',URL:'https://www.genenames.org/help/hcop/ ; https://fgr.hms.harvard.edu/diopt-documentation'},
  {Topic:'GenAge human',Method:'All 307 human candidate records are included. The source why field is retained in Evidence. These categories range from human/mammal/model evidence to functional, upstream/downstream and putative evidence; they do not all establish a human lifespan effect.',Resource:'GenAge',Version:'Build 21, 28 August 2023',URL:'https://genomics.senescence.info/genes/human_genes.zip'},
  {Topic:'GenAge models',Method:`All 2,205 original model-organism records are retained on GenAge models. ${d.stats.genage_model_mapped} have at least one retained human mapping; ${2205-d.stats.genage_model_mapped} remain without one. Animal and fungal homologs have separate flags. Source lifespan effect and pro-/anti-longevity or fitness annotations are retained; reduced lifespan or fitness alone is not proof of accelerated aging.`,Resource:'GenAge',Version:'Build 21, 28 August 2023',URL:'https://genomics.senescence.info/genes/models_genes.zip'},
  {Topic:'GenAge mapping',Method:'Human candidates are the union of the official GenAge human-homolog export and HCOP predictions supported by >=3 databases for mouse, worm, fly, budding yeast, fission yeast and zebrafish. This is a screening cutoff, not a probability. GenAge-supplied homologs have no comparable support score. Mapping sources and support are retained in the companion ortholog TSV.',Resource:'GenAge / HCOP',Version:'Retrieved 22 October 2026; HCOP files updated 1 October 2026',URL:'https://genomics.senescence.info/genes/search.php?export=true&homologs=true&organism=Homo+sapiens&show=5'},
  {Topic:'HCOP downloads',Method:'Bulk fifteen-column files for human_mouse, human_c.elegans, human_fruitfly, human_s.cerevisiae, human_s.pombe and human_zebrafish. Source gene IDs are used where available; species-specific source symbols are a fallback. Human names are harmonized after mapping.',Resource:'HCOP',Version:'File timestamps 1 October 2026',URL:'https://www.genenames.org/download/hcop/tsv/'},
  {Topic:'CellAge manipulation',Method:'866 unique genes across 949 senescence-type records. Gene-level Induces/Inhibits/Unclear labels and cancer-cell context are preserved with PMID links. These are senescence experiments in cells, not animal lifespan tests. A source gene-effect label must not be interpreted as a KO/OE recipe.',Resource:'CellAge',Version:'Build 3, 23 April 2023',URL:'https://genomics.senescence.info/cells/cellAge.zip'},
  {Topic:'CellAge expression',Method:'The 1,259-gene replicative-senescence expression signature is included but separately flagged. It is observational expression evidence, not perturbation evidence. Source p values are retained without relabelling them as FDR. The downloadable signature archive carries Build 2 notes; the Build 3 release says this signature was unchanged.',Resource:'CellAge',Version:'Signature archive Build 2, 23 December 2020; unchanged in Build 3',URL:'https://genomics.senescence.info/cells/cellSignatures.zip'},
  {Topic:'GWAS source',Method:'Full ontology-annotated association download from the NHGRI-EBI GWAS Catalog, with reported and mapped genes, variants, study accessions, PMID, sample descriptions and source statistics. The complete download is preserved locally. These are catalog associations, not newly calculated association tests.',Resource:'GWAS Catalog',Version:`API metadata release ${d.stats.gwas_release}; downloaded 22 October 2026`,URL:'https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip'},
  {Topic:'GWAS phenotype scope',Method:'Direct and parental longevity/lifespan/age at death; healthspan and global healthy-aging traits; biological/epigenetic age or age acceleration; combined aging traits. Excludes reproductive lifespan, cosmetic skin traits, disease-only, cognition-only and component-biomarker traits. Trait selection records the exact inclusions and exclusions. This is a scoped collection, not every possible aging-associated phenotype.',Resource:'GWAS Catalog',Version:'Trait rules in compile_reference.py',URL:'https://www.ebi.ac.uk/gwas/docs/methods/curation/'},
  {Topic:'GWAS significance',Method:`Only associations with catalog P <= 5e-8 contribute genes to the master list: ${d.stats.gwas_significant_associations} association rows from ${d.stats.gwas_significant_studies} study accessions and ${d.stats.gwas_significant_publications} publications. Catalog p values are rounded, so threshold-boundary values are approximate. No new FDR or replication claim is imposed. All ${d.stats.gwas_in_scope_associations} in-scope rows remain on GWAS associations, with inclusion flagged.`,Resource:'GWAS Catalog',Version:'P threshold 5e-8',URL:'https://www.ebi.ac.uk/gwas/docs/methods/curation/'},
  {Topic:'GWAS gene assignment',Method:'Union of author-reported genes and Catalog overlapping or flanking genes, with the assignment method retained per gene–association row. Entrez/Ensembl IDs recover resolvable symbols. Explicit paired-gene shorthand is expanded (CHRNA3/5, NOL4L/-DT, slash-separated symbols). Gene proximity does not identify a causal target. Multiple genes can share one locus, and expression or ATAC direction is not inferred from a variant effect.',Resource:'GWAS Catalog / HGNC',Version:'GRCh38.p14; Ensembl 116 from API metadata',URL:'https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/'},
  {Topic:'Unresolved GWAS records',Method:`${d.stats.gwas_unmapped_significant_associations} significant association rows have no resolvable current human gene and remain in the association table. Unresolved tokens, including source spelling errors and nonspecific RNA-family labels, are listed in the audit TSV. Missing gene assignment is not a negative association.`,Resource:'Audit',Version:'This collection',URL:''},
  {Topic:'Interpretation of direction',Method:'Animal perturbation direction, GenAge longevity influence, CellAge senescence effect and a GWAS allele effect are different measurements. Do not collapse them into one universal pro-/anti-aging label. None alone establishes the expected sign of cross-species liver ATAC versus MLS.',Resource:'Compilation',Version:'This collection',URL:''},
  {Topic:'Evidence subsets',Method:'Filter CellAge manipulation=Yes for cell-based perturbation evidence; GWAS longevity=Yes for direct/parental longevity rather than all aging traits; GenAge animal homolog=Yes to distinguish animal from fungal evidence. Expression-only candidates remain identifiable through the flags. Resource count is descriptive, not an evidence-quality score.',Resource:'Compilation',Version:'This collection',URL:''},
  {Topic:'Existing paper expression lists',Method:'The earlier Tyshkovskiy expression signatures and literature shortlist are separate project outputs and are not silently relabelled as animal lifespan experiments or added to this four-resource union. This workbook focuses on the four sources requested here.',Resource:'Scope',Version:'This collection',URL:''},
  {Topic:'Reproducibility',Method:'The companion ZIP contains the gene matrix, per-source symbol lists, evidence, original GenAge model records, experimental records, GWAS associations, trait decisions, mapping audit, unresolved identifiers and source checksums. Original downloaded archives and compilation scripts are retained in this output folder.',Resource:'Files',Version:'22 October 2026',URL:''},
];
table(methods,notes,'MethodsAndSources',[29,119,27,48,95]);

overview.showGridLines=false; overview.tabColor=navy;
overview.getRange('A1:D32').format={font,wrapText:true,verticalAlignment:'center',rowHeight:30};
for(const [i,w] of [31,17,66,62].entries())overview.getRange(`${col(i)}1:${col(i)}32`).format.columnWidth=w;
overview.getRange('A1:D1').values=[['Aging and longevity','Gene reference','Human symbols with original organism evidence','22 October 2026']];
overview.getRange('A1:D1').format={fill:navy,font:{...font,size:16,bold:true,color:'#FFFFFF'},rowHeight:64};
overview.getRange('A3:D3').values=[['Unique human genes',null,'Union across four labelled resources','Each gene appears once on Genes.']];
overview.getRange('B3').formulas=[[`=COUNTA(Genes!A2:A${gr.end})`]];
overview.getRange('A4:D4').values=[['Two or more resources',null,'Source memberships overlap','Source count is not a quality or replication score.']];
overview.getRange('B4').formulas=[[`=COUNTIF(Genes!G2:G${gr.end},">=2")`]];
overview.getRange('A3:D4').format={fill:green,rowHeight:45};
overview.getRange('B3:B4').format.font={...font,bold:true,size:18};
overview.getRange('A6:D6').values=[['Resource','Human genes','Included evidence','Meaning of the count']];
overview.getRange('A6:D6').format={fill:navy,font:{...font,bold:true,color:'#FFFFFF'},rowHeight:38};
const sourceRows=[
  ['Experiment','C','Previously curated animal lifespan-extension studies','62 human candidates from mouse genes and worm/fly homolog mappings.'],
  ['GenAge','D','Human aging candidates and human homologs of model genes','307 human entries plus mapped animal/fungal candidates; overlapping sets.'],
  ['CellAge','E','Senescence manipulations and a separate expression signature','866 manipulation genes and 1,259 expression genes; 1,955 in their union.'],
  ['GWAS','F','Genome-wide significant longevity, healthspan and biological-age associations','Reported and mapped candidate genes, not proven causal genes.'],
];
sourceRows.forEach(([s,c,inc,note],i)=>{
  const row=7+i;overview.getRange(`A${row}:D${row}`).values=[[s,null,inc,note]];
  overview.getRange(`B${row}`).formulas=[[`=COUNTIF(Genes!${c}2:${c}${gr.end},"Yes")`]];
  overview.getRange(`A${row}:D${row}`).format={fill:i%2?pale:'#FFFFFF',rowHeight:68};
});
overview.getRange('A12:D12').values=[['GWAS subset','Human genes','Trait class','Subsets overlap']];
overview.getRange('A12:D12').format={fill:navy,font:{...font,bold:true,color:'#FFFFFF'},rowHeight:38};
[['Longevity / lifespan','Q','Direct or parental longevity, lifespan or age at death'],['Healthspan','R','Healthspan or global healthy aging'],['Biological age','S','Age clocks and age acceleration'],['Combined aging traits','T','General or multivariate aging phenotypes']].forEach(([s,c,t],i)=>{
  const row=13+i;overview.getRange(`A${row}:D${row}`).values=[[s,null,t,'P <= 5e-8 in the Catalog']];
  overview.getRange(`B${row}`).formulas=[[`=COUNTIF(Genes!${c}2:${c}${gr.end},"Yes")`]];
  overview.getRange(`A${row}:D${row}`).format={fill:i%2?pale:'#FFFFFF',rowHeight:44};
});
overview.getRange('A18:D18').values=[['Use the gene list','', 'Filter the four resource columns and the more specific evidence flags.','Evidence links each gene to source records, organisms, phenotype and reference.']];
overview.getRange('A18:D18').format={fill:green,rowHeight:74};
overview.getRange('A20:D20').values=[['Interpretation','', 'Aging association, cellular senescence and lifespan extension are distinct endpoints.','Homolog mapping and GWAS proximity identify candidates; they do not prove a human lifespan effect.']];
overview.getRange('A20:D20').format={fill:amber,rowHeight:88};
overview.getRange('A22:D22').values=[['Mapping coverage',d.stats.genage_model_mapped,'Of 2,205 GenAge model records have a retained human mapping.','All original model records remain available, including unmapped genes.']];
overview.getRange('A22:D22').format={fill:pale,rowHeight:70};
overview.getRange('A24:D24').values=[['GWAS coverage',d.stats.gwas_significant_associations,'Significant association rows from 68 study accessions and 34 publications.','Ten significant rows have no resolvable human gene; all are retained on GWAS associations.']];
overview.getRange('A24:D24').format={fill:pale,rowHeight:70};
overview.getRange('A26:D26').values=[['Source versions','','GenAge Build 21; CellAge Build 3; HGNC and HCOP downloaded 22 October 2026.','GWAS Catalog API release metadata: 22 October 2026. Methods provides scope and download URLs.']];
overview.getRange('A26:D26').format={fill:pale,rowHeight:70};
overview.getRange('B3:B24').setNumberFormat('#,##0');
overview.getRange('B3:B24').format.horizontalAlignment='center';
overview.freezePanes.freezeRows(1);

wb.recalculate();
const inspections=[];
for(const range of ['A3:D10','A12:D16'])inspections.push(await wb.inspect({kind:'region',sheetId:'Overview',range,maxChars:5000,tableMaxRows:10,tableMaxCols:4}));
inspections.push(await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30},summary:'Formula error scan'}));
await fs.writeFile(path.join(root,'workbook_inspection.txt'),inspections.map(r=>r.ndjson??JSON.stringify(r)).join('\n'));
for(const [sheetName,range,filename] of [
  ['Overview','A1:D16','overview'],['Genes','A1:H9','genes'],['Evidence','A1:F6','evidence'],
  ['Experiments','A1:G5','experiments'],['GenAge models','A1:E6','models'],
  ['GWAS associations','A1:F6','gwas'],['GWAS trait selection','A1:C7','traits'],['Methods','A1:C4','methods'],
]){
  const preview=await wb.render({sheetName,range,scale:1.25,format:'png'});
  await fs.writeFile(path.join(root,`preview_${filename}.png`),new Uint8Array(await preview.arrayBuffer()));
}
const out=await SpreadsheetFile.exportXlsx(wb);
const output=path.join(root,'Aging_and_longevity_gene_reference.xlsx');
await out.save(output);
await fs.writeFile(path.join(root,'README.txt'),
  'AGING AND LONGEVITY GENE REFERENCE\nPrepared 22 October 2026\n\n'+
  JSON.stringify(d.stats,null,2)+'\n\n'+notes.map(n=>n.Topic+'\n'+n.Method+'\nSource: '+n.Resource+'; '+n.Version+'\n'+n.URL).join('\n\n')+'\n');
console.log(JSON.stringify({output,genes:d.master.length,evidence:d.evidence.length,sheets:8}));
