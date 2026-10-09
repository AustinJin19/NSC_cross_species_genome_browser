"""Build local lifespan metadata; preserve source provenance and exclude unavailable values."""
import csv,sqlite3,zipfile,io,json,math
from pathlib import Path
root=Path(__file__).resolve().parent
species=[r[0] for r in sqlite3.connect(root/'data/liver_omics.sqlite').execute('select distinct species from measurements order by species')]
lab={r['species']:r for r in csv.DictReader(open(root.parents[2]/'Open4Gene-main'/'species_MLS.csv'))}
with zipfile.ZipFile(root.parents[1]/'Anage.zip') as z:
 anage={r['Genus']+' '+r['Species']:r for r in csv.DictReader(io.StringIO(z.read('anage_data.txt').decode()),delimiter='\t')}
result={}
for sp in species:
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
