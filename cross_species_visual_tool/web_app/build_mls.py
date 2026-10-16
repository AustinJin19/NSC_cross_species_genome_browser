"""Build local lifespan metadata; preserve source provenance and exclude unavailable values."""
import csv,sqlite3,zipfile,io,json,math
from pathlib import Path
from server import SPECS
root=Path(__file__).resolve().parent
species=[r[0] for r in sqlite3.connect(root/'data/liver_omics.sqlite').execute('select distinct species from measurements order by species')]
species=sorted(set(species) | {s[2] for s in SPECS})
lab_path=root.parents[2]/'Open4Gene-main'/'species_MLS.csv'
lab={r['species']:r for r in csv.DictReader(lab_path.open())} if lab_path.exists() else {}
existing=json.loads((root/'data/species_mls.json').read_text()) if (root/'data/species_mls.json').exists() else {}
with zipfile.ZipFile(root.parents[1]/'data/omics/lifespan/Anage.zip') as z:
 anage={r['Genus']+' '+r['Species']:r for r in csv.DictReader(io.StringIO(z.read('anage_data.txt').decode()),delimiter='\t')}
result={}
for sp in species:
 if not lab_path.exists() and existing.get(sp,{}).get('years') is not None:
  result[sp]=existing[sp];continue
 row=lab.get(sp)
 if row and 'unverified' not in row['status']:
  value=float(row['MLS_corrected']);source='Open4Gene-main/species_MLS.csv · MLS_corrected';note=row['source_note']
 elif sp in anage and anage[sp]['Maximum longevity (yrs)']:
  value=float(anage[sp]['Maximum longevity (yrs)']);source='Anage.zip · anage_data.txt · Maximum longevity (yrs)';note='Supplied July 2023 snapshot; exact scientific-name match'
 else:
  result[sp]={'years':None,'source':'No verified local reference','note':'Excluded from correlation'};continue
 result[sp]={'years':value,'log10_years':math.log10(value),'source':source,'note':note}
(root/'data/species_mls.json').write_text(json.dumps(result,indent=2)+'\n')
print('Matched',sum(r['years'] is not None for r in result.values()),'of',len(result))

# User-supplied exact traits and requested AnAge fallbacks take precedence.
if (root.parents[1]/"data/metadata/species_traits_user_2026-10-15.tsv").exists():
 from apply_user_traits import apply
 apply()
