import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const root=path.dirname(fileURLToPath(import.meta.url));
const data=JSON.parse(await fs.readFile(path.join(root,'report_data.json'),'utf8'));
const windows=data.summary.map(r=>r.window_kb);
const wb=Workbook.create();
const summary=wb.worksheets.add('Summary');
const sig=wb.worksheets.add('Significant genes');
const joint=wb.worksheets.add('Joint FDR genes');
const comparison=wb.worksheets.add('Across windows');
const species=wb.worksheets.add('Species');
const methods=wb.worksheets.add('Methods');
const font={name:'Arial',size:11,color:'#243746'};
const navy='#183B55',blue='#EAF1F7',green='#E2F0E9',red='#FBE8E5';
const col=i=>{let s='';for(i++;i;i=Math.floor((i-1)/26))s=String.fromCharCode(65+(i-1)%26)+s;return s;};
function header(s,range){s.getRange(range).format={fill:navy,font:{...font,bold:true,color:'#FFFFFF'},wrapText:true,rowHeight:38,verticalAlignment:'center',horizontalAlignment:'center'};}
function table(s,name,keys,rows,widths){
  const end=rows.length+1,last=col(keys.length-1);
  s.getRange(`A1:${last}${end}`).values=[keys,...rows];
  s.getRange(`A1:${last}${end}`).format={font,rowHeight:23,verticalAlignment:'center'};
  s.tables.add(`A1:${last}${end}`,true,name);
  s.getRange(`A2:${last}${end}`).format.fill='#FFFFFF';
  for(let i=2;i<=end;i+=2)s.getRange(`A${i}:${last}${i}`).format.fill='#F3F7FA';
  header(s,`A1:${last}1`);
  widths.forEach((w,i)=>s.getRange(`${col(i)}1:${col(i)}${end}`).format.columnWidth=w);
  s.freezePanes.freezeRows(1);s.freezePanes.freezeColumns(2);s.showGridLines=false;
  return end;
}

const sigKeys=['TSS flank (kb)','Gene','Direction','Species (n)','Partial r','p-value','BH FDR within window','BH FDR across all windows','MLS coefficient','Mass coefficient','df'];
const sigRows=data.significant.map(r=>[r.window_kb,r.gene,r.direction,r.n_species,r.partial_r,r.p_value,r.BH_FDR_window,r.BH_FDR_all_windows,r.MLS_coefficient,r.mass_coefficient_per_log10g,r.df]);
const se=table(sig,'SignificantGenes',sigKeys,sigRows,[16,22,15,15,15,17,21,23,20,20,12]);
sig.getRange(`E2:E${se}`).setNumberFormat('0.000');
sig.getRange(`F2:H${se}`).setNumberFormat('0.000E+00');
sig.getRange(`I2:J${se}`).setNumberFormat('0.0000');
sig.getRange(`C2:C${se}`).conditionalFormats.add('containsText',{text:'Positive',format:{fill:green}});
sig.getRange(`C2:C${se}`).conditionalFormats.add('containsText',{text:'Negative',format:{fill:red}});
const jointRows=data.joint_significant.map(r=>[r.window_kb,r.gene,r.direction,r.n_species,r.partial_r,r.p_value,r.BH_FDR_window,r.BH_FDR_all_windows,r.MLS_coefficient,r.mass_coefficient_per_log10g,r.df]);
const je=table(joint,'JointSignificantGenes',sigKeys,jointRows,[16,22,15,15,15,17,21,23,20,20,12]);
joint.getRange(`E2:E${je}`).setNumberFormat('0.000');
joint.getRange(`F2:H${je}`).setNumberFormat('0.000E+00');
joint.getRange(`I2:J${je}`).setNumberFormat('0.0000');
joint.getRange(`C2:C${je}`).conditionalFormats.add('containsText',{text:'Positive',format:{fill:green}});
joint.getRange(`C2:C${je}`).conditionalFormats.add('containsText',{text:'Negative',format:{fill:red}});

const compKeys=['Gene','Significant windows (n)','Significant TSS flanks (kb)',...windows.flatMap(k=>[`Partial r ±${k} kb`,`BH FDR ±${k} kb`,`Species n ±${k} kb`])];
const compRows=data.comparison.map(r=>[r.gene,r.significant_windows,r.windows_kb,...windows.flatMap(k=>[r[`r_${k}kb`]??'n.a.',r[`FDR_${k}kb`]??'n.a.',r[`n_${k}kb`]??'n.a.'])]);
const ce=table(comparison,'AcrossWindows',compKeys,compRows,[22,22,34,...windows.flatMap(()=>[19,19,19])]);
windows.forEach((k,i)=>{
  const c=3+i*3;
  comparison.getRange(`${col(c)}2:${col(c)}${ce}`).setNumberFormat('0.000');
  comparison.getRange(`${col(c+1)}2:${col(c+1)}${ce}`).setNumberFormat('0.000E+00');
  comparison.getRange(`${col(c+1)}2:${col(c+1)}${ce}`).conditionalFormats.add('cellIs',{operator:'lessThan',formula:0.05,format:{fill:green,font:{bold:true}}});
});

const speKeys=['Species ID','Common name','Scientific name','MLS (years)','Body mass (g)','Analysis status','Trait source','Trait note'];
const speRows=data.species.map(r=>[r.species,r.common_name,r.scientific_name,r.MLS_years??'n.a.',r.body_mass_g??'n.a.',r.analysis_status,r.trait_source,r.trait_note]);
const pe=table(species,'SpeciesTraits',speKeys,speRows,[30,28,31,17,20,37,80,74]);
species.getRange(`D2:E${pe}`).setNumberFormat('#,##0.0');
species.getRange(`F2:H${pe}`).format={wrapText:true,rowHeight:56};

summary.showGridLines=false;summary.tabColor=navy;
summary.getRange('A1:J35').format={font,rowHeight:25,verticalAlignment:'center'};
[22,23,20,20,22,22,22,22,22,22].forEach((w,i)=>summary.getRange(`${col(i)}1:${col(i)}35`).format.columnWidth=w);
summary.getRange('A2').values=[['Regulatory activity and MLS across TSS windows']];
summary.getRange('A2').format.font={...font,size:16,bold:true};
summary.getRange('A3').values=[['Log₁₀ median peak count versus log₁₀ MLS, adjusted for log₁₀ body mass.']];
summary.getRange('A3').format.font={...font,italic:true};
summary.getRange('A5:J5').values=[['TSS flank (kb)','Full width (kb)','Genes tested','FDR positive','FDR negative','Joint FDR positive','Joint FDR negative','Nominal positive','Nominal negative','Minimum window FDR']];
header(summary,'A5:J5');
summary.getRange('A6:J11').values=data.summary.map(r=>[r.window_kb,r.full_width_kb,r.tested,null,null,r.joint_FDR_positive,r.joint_FDR_negative,r.nominal_positive,r.nominal_negative,r.minimum_FDR]);
data.summary.forEach((r,i)=>{
 const row=6+i;
 summary.getRange(`D${row}`).formulas=[[`=COUNTIFS('Significant genes'!$A$2:$A$${se},A${row},'Significant genes'!$C$2:$C$${se},"Positive")`]];
 summary.getRange(`E${row}`).formulas=[[`=COUNTIFS('Significant genes'!$A$2:$A$${se},A${row},'Significant genes'!$C$2:$C$${se},"Negative")`]];
 summary.getRange(`F${row}`).formulas=[[`=COUNTIFS('Joint FDR genes'!$A$2:$A$${je},A${row},'Joint FDR genes'!$C$2:$C$${je},"Positive")`]];
 summary.getRange(`G${row}`).formulas=[[`=COUNTIFS('Joint FDR genes'!$A$2:$A$${je},A${row},'Joint FDR genes'!$C$2:$C$${je},"Negative")`]];
});
summary.getRange('A6:I11').setNumberFormat('#,##0');
summary.getRange('J6:J11').setNumberFormat('0.000E+00');
summary.getRange('D6:E11').format.fill=blue;
summary.getRange('A13').values=[['Distinct genes significant in at least one window']];
summary.getRange('D13').formulas=[[`=COUNTA('Across windows'!A2:A${ce})`]];
summary.getRange('A14').values=[['Genes significant in all six windows']];
summary.getRange('D14').formulas=[[`=COUNTIFS('Across windows'!B2:B${ce},6)`]];
summary.getRange('A16').values=[['Shared significant genes between windows']];
summary.getRange('A16').format.font={...font,bold:true};
summary.getRange('A17:G17').values=[['TSS flank (kb)',...windows.map(k=>`±${k} kb`)]];
header(summary,'A17:G17');
summary.getRange('A18:G23').values=data.overlap.map(r=>[r.window_kb,...windows.map(k=>r[`overlap_${k}kb`])]);
summary.getRange('B18:G23').conditionalFormats.add('colorScale',{colors:['#FFFFFF','#4F8EB9'],thresholds:['min','max']});
summary.getRange('A25').values=[['SIRT6 across windows']];
summary.getRange('A25').format.font={...font,bold:true};
summary.getRange('A26:F26').values=[['TSS flank (kb)','Species (n)','Partial r','p-value','Window FDR','Joint FDR']];
header(summary,'A26:F26');
summary.getRange('A27:F32').values=data.sirt6.map(r=>r?[r.window_kb,r.n_species,r.partial_r,r.p_value,r.BH_FDR_window,r.BH_FDR_all_windows]:[null,null,null,null,null,null]);
summary.getRange('C27:C32').setNumberFormat('0.000');
summary.getRange('D27:F32').setNumberFormat('0.000E+00');
summary.getRange('A34').values=[['FDR < 0.05. Window FDR is corrected within each window; joint FDR covers all six windows.']];
summary.getRange('A35').values=[['Species and selected loci can differ across windows because the available-region filter is applied separately.']];

methods.showGridLines=false;
const notes=[
 ['Scope','Genome-wide regulatory activity at TSS ±20, ±50, ±100, ±200, ±500 and ±1,000 kb. These are flanks; total widths are 40 kb to 2 Mb.'],
 ['Measurement','Count each sample’s ATAC peaks overlapping the half-open TSS window. A peak crossing a window edge is counted once. No CPM integration or distance weighting.'],
 ['Species summary','Take the median across available samples at each locus. For multiple loci, choose the highest median first, with lexical locus tie-breaking. Apply eligibility afterward, without falling back to another locus.'],
 ['Availability','At least 50% of the requested window must fall within the matching chromosome/scaffold bounds for every measured sample. This checks genomic span, not ATAC coverage or internal assembly gaps.'],
 ['Regression','OLS: log₁₀(median peak count) = intercept + β₁ log₁₀(MLS in years) + β₂ log₁₀(body mass in g). Partial r residualizes both activity and MLS against log₁₀ body mass. No gene-length covariate.'],
 ['Filters','Exclude human, macaque, and bowhead whale. Exclude missing/lower-bound MLS, missing body mass, nonpositive median count, insufficient available span and degenerate fits. At least 5 species per gene; df = n − 3.'],
 ['Multiplicity','Primary results use Benjamini–Hochberg FDR < 0.05 separately within each window, matching the earlier analysis. A second BH FDR pools all tested gene–window associations. Nominal p < 0.05 is reported separately.'],
 ['Cross-window comparison','Species sets and selected loci can change across windows; counts are not a controlled comparison using identical species at every window. Nested windows and neighboring genes share peaks.'],
 ['Biological interpretation','Gene symbols are matched across species; this is not a verified one-to-one orthology set. The model does not adjust for phylogeny, batch, sequencing depth or global peak-call differences. Association is exploratory, not causal.'],
 ['Source data','Local gene_activity.sqlite, species_mls.json, configured ATAC peak files and matching bigWig chromosome bounds, as of 22 October 2026. Original source files were not modified. Input metadata and source paths are recorded in manifest.json.'],
 ['Validation',`${data.validation.browser_scorer_locus_windows.toLocaleString()} locus-window checks matched the browser scorer. The 500 kb analysis reproduced every tested gene and the original 71 significant positive genes. Maximum r/p/FDR difference: ${data.validation.maximum_500kb_difference.toExponential(2)}.`],
 ['Full results','The companion archive includes all tested associations and positive/negative gene lists for each window, nominal and FDR thresholds, coverage summaries, and the source manifest. The workbook includes window-FDR and joint-FDR lists, with cross-window comparisons for genes significant within any window.'],
];
table(methods,'AnalysisMethods',['Item','Definition'],notes,[28,132]);
methods.getRange(`A2:B${notes.length+1}`).format={wrapText:true,rowHeight:66};
methods.freezePanes.unfreeze();
wb.recalculate();
const inspected=await wb.inspect({kind:'region',sheetId:'Summary',range:'A5:J14',maxChars:3500,tableMaxRows:10,tableMaxCols:10});
await fs.writeFile(path.join(root,'workbook_inspection.txt'),inspected.ndjson??JSON.stringify(inspected));
for(const [sheetName,range,file] of [
 ['Summary','A2:J14','preview_summary.png'],
 ['Significant genes','A1:H8','preview_significant.png'],
 ['Joint FDR genes','A1:H8','preview_joint.png'],
 ['Across windows','A1:L8','preview_comparison.png'],
 ['Species','A1:F8','preview_species.png'],
 ['Methods','A1:B8','preview_methods.png']
]){
 const b=await wb.render({sheetName,range,scale:1.25,format:'png'});
 await fs.writeFile(path.join(root,file),new Uint8Array(await b.arrayBuffer()));
}
const output=await SpreadsheetFile.exportXlsx(wb);
const filename=path.join(root,'Regulatory_activity_MLS_six_windows.xlsx');
await output.save(filename);
console.log(JSON.stringify({output:filename,significant_gene_window_pairs:sigRows.length,distinct_significant_genes:compRows.length}));
