import fs from 'node:fs/promises';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';
const out=new URL('.',import.meta.url).pathname;
const d=JSON.parse(await fs.readFile(out+'locus_analysis.json','utf8'));
const wb=Workbook.create();
function col(n){let s='';for(n++;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;}
function sheet(name,title,subtitle,headers,rows,widths){
  const s=wb.worksheets.add(name); s.showGridLines=false;
  const end=col(headers.length-1), last=rows.length+6;
  s.getRange(`A1:${end}${last}`).format.font={name:'Arial',size:10,color:'#243746'};
  s.getRange(`A1:${end}${last}`).format.rowHeight=22;
  s.getRange(`A1:${end}${last}`).format.verticalAlignment='center';
  s.getRange('A2').values=[[title]];s.getRange('A2').format.font={size:16,bold:true,color:'#183D57'};
  s.getRange('A2').format.rowHeight=28;
  s.getRange('A4').values=[[subtitle]];s.getRange('A4').format.font={size:10,italic:true,color:'#526779'};
  s.getRange(`A6:${end}6`).values=[headers];
  if(rows.length)s.getRange(`A7:${end}${last}`).values=rows;
  s.tables.add(`A6:${end}${last}`,true,name.replaceAll(' ','')+'Table');
  s.getRange(`A6:${end}6`).format={fill:'#183D57',font:{name:'Arial',size:10,color:'#FFFFFF',bold:true},wrapText:true,horizontalAlignment:'center',verticalAlignment:'center',rowHeight:38};
  for(let i=0;i<widths.length;i++)s.getRange(`${col(i)}1:${col(i)}${last}`).format.columnWidth=widths[i];
  for(let i=7;i<=last;i+=2)s.getRange(`A${i}:${end}${i}`).format.fill='#F1F5F8';
  s.freezePanes.freezeRows(6);s.freezePanes.freezeColumns(name==='Species loci'?2:1);
  return s;
}
const cr=[...d.clusters].sort((a,b)=>b.n_genes-a.n_genes||a.cluster.localeCompare(b.cluster));
const clusters=sheet('Locus groups','Locations of the 71 regulatory hits',
  '69 mapped human genes: 18 multi-gene groups containing 63 genes, plus 6 singletons. Reference: hg38.',
  ['Group','Chromosome','First TSS (1-based)','Last TSS (1-based)','TSS span (kb)','Genes (n)','Common window (kb)','Genes in genomic order'],
  cr.map(r=>[r.cluster,r.chrom,r.first_tss_1based,r.last_tss_1based,null,r.n_genes,null,r.genes]),
  [9,13,20,20,16,11,21,83]);
clusters.tabColor='#183D57';
clusters.getRange(`C7:D30`).setNumberFormat('#,##0');
clusters.getRange('E7:E30').setNumberFormat('0.0');clusters.getRange('G7:G30').setNumberFormat('0.0');
clusters.getRange('H7:H30').format.wrapText=true;clusters.getRange('A7:H30').format.rowHeight=32;
for(let i=0;i<cr.length;i++){
 const r=i+7;clusters.getRange(`E${r}`).formulas=[[`=(D${r}-C${r})/1000`]];
 clusters.getRange(`G${r}`).formulas=[[`=MAX(0,1000-E${r})`]];
}
const gs=sheet('Human genes','Human reference coordinates and original statistics',
  'Coordinates are 1-based inclusive. TSS is the browser-selected human transcript anchor; it is not an original fitted-species window.',
 ['Original gene','Human symbol','Group','Genes in group','Chromosome','Cytoband','Gene start (bp)','Gene end (bp)','TSS (bp)','Strand','Species in fit','Partial r','P value','BH FDR','Mapping','Transcript','TSS selection','UCSC link'],
 d.genes.map(r=>[r.gene,r.human_symbol,r.cluster||'Unmapped',r.cluster_genes??null,r.chrom||'Unmapped',r.cytoband,r.start_1based??null,r.end_1based??null,r.tss_1based??null,r.strand||'',r.n_species,r.partial_r,r.p_value,r.BH_FDR,r.mapping,r.transcript||'',r.tss_method||'',r.ucsc_url||'']),
 [16,16,12,15,14,23,18,18,18,9,15,12,13,13,40,28,28,80]);
gs.getRange('G7:I77').setNumberFormat('#,##0');gs.getRange('L7:L77').setNumberFormat('0.000');gs.getRange('M7:M77').setNumberFormat('0.000E+00');gs.getRange('N7:N77').setNumberFormat('0.0000');
gs.getRange('O7:O77').format.wrapText=true;
gs.getRange('A76:R77').format.rowHeight=34;
const pairs=sheet('Window overlaps','Pairs with overlapping windows in human or fitted species',
 'Species fractions use species included for both genes. Overlap fraction uses the shorter clipped window, among overlapping species.',
 ['Gene A','Gene B','Human TSS distance (bp)','Human overlap (bp)','Species both fitted','Same contig (n)','Overlapping (n)','Overlapping species (%)','Median overlap (%)','Overlapping species'],
 d.pairs.map(r=>[r.gene_a,r.gene_b,r.human_tss_distance_bp,r.human_window_overlap_bp,r.species_both_in_fit,r.species_same_contig,r.species_overlapping,r.overlap_species_fraction,r.median_overlap_fraction,r.overlapping_species]),
 [16,16,23,21,20,18,18,24,22,105]);
pairs.getRange(`C7:D${d.pairs.length+6}`).setNumberFormat('#,##0');pairs.getRange(`H7:I${d.pairs.length+6}`).setNumberFormat('0.0%');
const sr=[...d.species_loci].sort((a,b)=>a.gene.localeCompare(b.gene)||a.species.localeCompare(b.species));
const native=sheet('Species loci','Selected loci from the original 500 kb analysis',
 'Native assembly coordinates. Window start is 0-based; end is exclusive. Fit inclusion matches the original 2026-10-19 results.',
 ['Gene','Species ID','Scientific name','Contig','Gene start (1-based)','Gene end (1-based)','Strand','TSS (1-based)','Clipped window start0','Clipped window end','Median peak count','Available fraction','Included in fit','Requested start0','Requested end'],
 sr.map(r=>[r.gene,r.species,r.scientific_name,r.chrom,r.start_1based,r.end_1based,r.strand,r.tss_1based,r.window_start0,r.window_end,r.median_peak_count,r.available_fraction,r.included_in_original_fit?'Yes':'No',r.requested_window_start0,r.requested_window_end]),
 [16,29,29,27,23,23,10,22,24,24,21,20,18,23,23]);
for(const c of ['E','F','H','I','J','N','O'])native.getRange(`${c}7:${c}${sr.length+6}`).setNumberFormat('#,##0');
native.getRange(`K7:K${sr.length+6}`).setNumberFormat('0.0');native.getRange(`L7:L${sr.length+6}`).setNumberFormat('0.0%');
native.getRange(`M7:M${sr.length+6}`).conditionalFormats.add('containsText',{text:'No',format:{font:{color:'#9C3B2E'}}});
const methods=[
 ['Hit list','FDR05_positive.tsv from Regulatory_500kb_MLS_and_GO_2026-10-19.xlsx analysis: 71 positive genes, no negative genes at BH FDR < 0.05.'],
 ['Original model','log10(species median peak count) versus log10(MLS), with log10(body mass) as a covariate. Ordinary partial correlation, not PGLS.'],
 ['Original filters','Human, macaque and bowhead whale excluded. Require at least 50% of requested region, positive peak count, at least 5 species. Statistics are preserved, not re-estimated.'],
 ['Human annotation','GENCODE v50 primary assembly hg38. Gene bodies combine transcript spans on the same chromosome and strand. TSS uses the locally selected MANE/canonical anchor.'],
 ['GENCODE source','https://www.gencodegenes.org/human/release_50.html'],
 ['Human symbols','HGNC complete-set snapshot in the project. Exact symbol matching except B3GNTL1 -> QTGAL (previous approved symbol). Same symbols are not proof of orthology.'],
 ['Unmapped human genes','MS4A20 and ZFP865 have no unique human symbol/anchor match in these resources. No speculative human ortholog was substituted.'],
 ['Mouse coordinates for MS4A20','mm10 chr19:11,096,814–11,130,889 (-), TSS 11,130,889. Original selected locus.'],
 ['Mouse coordinates for ZFP865','mm10 chr7:5,020,376–5,033,226 (+), TSS 5,020,376. Original selected locus.'],
 ['Group definition','Connect genes on the same human chromosome when adjacent TSSs are less than 1,000,000 bp apart. This yields 24 groups, including 6 singletons.'],
 ['Common interval','The common overlap of all theoretical ±500 kb windows in a group. One chain (L21) spans >1 Mb and has no interval shared by every member.'],
 ['Human coordinate scope','Human is a common positional reference only and remains excluded from the original correlation. Human grouping is not an independent-locus count across all species.'],
 ['Native locus selection','Original cached species_scores retain all loci. Select the largest median per symbol/species, using lexicographic first for ties, before availability filtering; no fallback.'],
 ['Native clipping','Cached requested start/end can extend beyond a contig. Actual start = max(0, requested start); actual end = actual start + cached available bp.'],
 ['Pair table coverage','Includes a pair if its human ±500 kb windows overlap or its original windows overlap in at least one jointly fitted species. Missing human distance is blank.'],
 ['Pair overlap fraction','Intersection bp / smaller available window length. Median is calculated only among species where the pair overlaps; the number of such species is shown separately.'],
 ['Correlation verification',`All 71 original partial r values reproduced from the recovered loci and archived traits; maximum absolute difference ${d.summary.maximum_original_r_reproduction_error.toExponential(3)}.`],
 ['Interpretation','Nearby genes can count many of the same peaks. These are correlated regional signals, and proximity alone does not identify which gene a peak regulates.'],
 ['Clustering inference','Observed location and overlap summary only. No enrichment test against genomic gene density, and no claim that 71 genes represent 71 independent loci.'],
 ['Sources','Original result: FDR05_positive.tsv; archived species_scores and species_provenance.json; local Human_models.tsv.gz and Human_anchors.json; hgnc_complete_set.tsv.'],
];
const notes=sheet('Methods','Definitions and sources','Analysis date: 2026-11-03. Original significant-gene result: 2026-10-19.',
 ['Item','Definition or source'],methods,[32,130]);
notes.getRange(`B7:B${methods.length+6}`).format.wrapText=true;
notes.getRange(`A7:B${methods.length+6}`).format.rowHeight=40;
notes.tabColor='#8D9AA5';
wb.recalculate();
console.log((await wb.inspect({kind:'table',range:"'Locus groups'!A6:H11",include:'values,formulas',tableMaxRows:6,tableMaxCols:8,maxChars:2200})).ndjson);
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30},summary:'Final formula error scan',maxChars:1000})).ndjson);
for(const [name,range,suffix] of [['Locus groups','A1:H13','clusters'],['Human genes','A1:N12','human'],['Human genes','O70:R77','mapping'],['Window overlaps','A1:I13','pairs'],['Species loci','A1:J12','species'],['Species loci','K6:O12','species_values'],['Methods','A1:B15','methods']]){
 const blob=await wb.render({sheetName:name,range,scale:1,format:'png'});
 await fs.writeFile(out+'preview_'+suffix+'.png',new Uint8Array(await blob.arrayBuffer()));
}
const file=await SpreadsheetFile.exportXlsx(wb);await file.save(out+'Regulatory_71_gene_locations.xlsx');
console.log('Workbook saved');
