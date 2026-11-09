import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {createRequire} from 'node:module';
const here=path.dirname(fileURLToPath(import.meta.url));
const require=createRequire(import.meta.url);
const stats=require('../../cross_species_visual_tool/web_app/www/statistics.js');
const validation=JSON.parse(await fs.readFile(path.join(here,'validation_points.json'),'utf8'));
for(const [gene,{points,result}] of Object.entries(validation)){
 const fit=stats.massAdjustedFit(points);
 for(const [a,b] of [[fit.r,result.partial_r],[fit.p,result.p_value],[fit.slope,result.MLS_coefficient],[fit.massSlope,result.mass_coefficient_per_log10g]]){
  if(Math.abs(a-b)>1e-9)throw Error(`${gene} does not match browser statistics: ${a} vs ${b}`);
 }
 if(fit.df!==result.df)throw Error('Degrees of freedom mismatch');
 console.log(gene,JSON.stringify({n:points.length,r:fit.r,p:fit.p,MLS_slope:fit.slope,mass_slope:fit.massSlope}));
}
const api=JSON.parse(await fs.readFile(path.join(here,'USP9X_browser.json'),'utf8'));
const groups=new Map();
for(const r of api.records){if(!groups.has(r.species))groups.set(r.species,[]);groups.get(r.species).push(r);}
for(const p of validation.USP9X.points){
 const selected=groups.get(p.species).sort((a,b)=>Number(Number.isFinite(b.median))-Number(Number.isFinite(a.median))||(b.median-a.median)||a.locus.localeCompare(b.locus))[0];
 if(selected.locus!==p.locus||selected.median!==p.median||Math.abs(selected.available_fraction-p.available_fraction)>1e-12||selected.region_eligible!==true)throw Error('API score mismatch for '+p.species);
}
console.log('USP9X species scores, locus selection and availability match the live API.');
