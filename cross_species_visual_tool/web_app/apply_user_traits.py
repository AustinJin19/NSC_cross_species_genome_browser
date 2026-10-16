"""Apply the user's 2026-10-15 MLS/body-mass table, retaining censored values."""
import argparse,csv,json,math,zipfile,io
from pathlib import Path
from server import SPECS
ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'data/metadata/species_traits_user_2026-10-15.tsv'
DEST=Path(__file__).parent/'data/species_mls.json'
def apply():
 rows={r['source_name']:r for r in csv.DictReader(SOURCE.open(),delimiter='\t')}
 aliases={'Nannospalax galili':'Spalax galili','Bramus lutescens':'Ellobius lutescens'}
 with zipfile.ZipFile(ROOT/'data/omics/lifespan/Anage.zip') as z:
  anage={r['Genus']+' '+r['Species']:r for r in csv.DictReader(io.StringIO(z.read('anage_data.txt').decode()),delimiter='\t')}
 names=set(rows)|{s[2] for s in SPECS};result={}
 for name in sorted(names):
  source=aliases.get(name,name);r=rows.get(source)
  if r is None:
   result[name]={'years':None,'body_mass_g':None,'source':'Not matched to user-provided table','note':'Excluded from correlation; no substitution from older references'};continue
  text=r['MLS_display'];bound=text.startswith('>');uncertain='*' in text;value=float(text.lstrip('>').rstrip('*'));mass=float(r['adult_body_mass_g'])
  result[name]={'years':None if bound or uncertain else value,'log10_years':None if bound or uncertain else math.log10(value),'display_years':text,'body_mass_g':mass,'log10_body_mass_g':math.log10(mass),'lower_bound_years':value if bound else None,'source':'User-provided figure, 2026-10-15; species_traits_user_2026-10-15.tsv','source_species':source,'note':('Lower bound, excluded from ordinary regression' if bound else 'Asterisk qualifier not defined in supplied image' if uncertain else 'Exact value transcribed from supplied figure')+('; mapped taxonomic name' if source!=name else '')}
 # User requested original AnAge values for ASM and ambiguous/missing MLS.
 for name in sorted(names):
  entry=result[name]
  if name!='Acomys cahirinus' and entry.get('years') is not None:continue
  lookup='Spalax ehrenbergi' if name=='Nannospalax galili' else name
  a=anage.get(lookup,{})
  if a.get('Maximum longevity (yrs)'):
   value=float(a['Maximum longevity (yrs)'])
   entry.update(years=value,log10_years=math.log10(value),display_years=str(value),lower_bound_years=None,source='AnAge build 15, supplied July 2023 snapshot',source_species=lookup,note='AnAge fallback requested by user')
   if name=='Nannospalax galili':entry['note']='Original project BMR value; Spalax ehrenbergi complex proxy, not a species-specific galili record. AnAge attributes 20.2 years to golani.'
  if (name=='Acomys cahirinus' or entry.get('body_mass_g') is None) and a.get('Adult weight (g)'):
   mass=float(a['Adult weight (g)']);entry.update(body_mass_g=mass,log10_body_mass_g=math.log10(mass))
 DEST.write_text(json.dumps(result,indent=2)+'\n')
 print('Applied supplied exact values with requested AnAge fallbacks')
if __name__=='__main__':apply()
