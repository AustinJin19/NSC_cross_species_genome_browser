"""Read-only liver RNA/protein queries from the workbook-derived SQLite index."""
import csv
import io
import json
import sqlite3
from contextlib import closing
from functools import lru_cache
from pathlib import Path

DB = Path(__file__).resolve().parent / 'data' / 'liver_omics.sqlite'


def connect():
    if not DB.exists():
        raise ValueError('Liver data index is not ready. Run build_omics.py first.')
    con=sqlite3.connect(f'file:{DB}?mode=ro',uri=True)
    con.row_factory=sqlite3.Row
    return con


@lru_cache(maxsize=2)
def catalog_for_version(version):
    with closing(connect()) as con:
        catalog=[dict(r) for r in con.execute('SELECT assay,species,count(DISTINCT sample) AS samples FROM measurements WHERE NOT (assay="rna" AND sample="BMR_Turk_Lu") GROUP BY assay,species ORDER BY species')]
        sources=json.loads(con.execute("SELECT value FROM metadata WHERE key='sources'").fetchone()[0])
    return catalog,sources


def query(gene):
    gene=str(gene).strip().upper()
    if not gene or len(gene)>100: raise ValueError('Enter a gene symbol')
    with closing(connect()) as con:
        records=[dict(r) for r in con.execute('SELECT * FROM measurements WHERE symbol=? AND NOT (assay="rna" AND sample="BMR_Turk_Lu") ORDER BY assay,species,sample,source_row',(gene,))]
    catalog,sources=catalog_for_version(DB.stat().st_mtime_ns)
    stale=[]
    for source in sources.values():
        path=Path(__file__).resolve().parent.parent.parent/source['file']
        if not path.exists() or path.stat().st_size!=source['size'] or path.stat().st_mtime_ns!=source['mtime_ns']:
            stale.append(source['file'])
    lifespans=json.loads((DB.parent/'species_mls.json').read_text())
    return dict(gene=gene,tissue='Liver',records=records,catalog=catalog,sources=sources,stale=stale,lifespans=lifespans)


def export_csv(gene):
    data=query(gene)
    out=io.StringIO();fields=['assay','symbol','sample','species','raw','log','normalized','norm_on','source','source_row','workbook','MLS_years','log10_MLS','MLS_source']
    writer=csv.DictWriter(out,fieldnames=['tissue']+fields,extrasaction='ignore');writer.writeheader()
    for row in data['records']:
        lifespan=data['lifespans'].get(row['species'],{})
        writer.writerow(dict(tissue='Liver',**row,MLS_years=lifespan.get('years'),log10_MLS=lifespan.get('log10_years'),MLS_source=lifespan.get('source')))
    return out.getvalue()
