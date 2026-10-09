"""Stream liver records from the supplied XLSX files into a queryable SQLite index.
Run with Python and lxml. Source workbooks are never modified.
"""
import json
import math
import sqlite3
import time
from pathlib import Path
from zipfile import ZipFile
from lxml import etree
from mouse_supplement import add_mouse_records

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent.parent
DB = HERE / 'data' / 'liver_omics.sqlite'
NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'


def rows(path):
    # Source workbook dimensions incorrectly say A1; stream the actual rows.
    with ZipFile(path) as z:
        strings = [''.join(si.itertext()) for si in etree.fromstring(z.read('xl/sharedStrings.xml'))]
        with z.open('xl/worksheets/sheet1.xml') as f:
            for _, row in etree.iterparse(f, events=('end',), tag=NS+'row'):
                values = {}
                for cell in row:
                    col = ''.join(c for c in cell.get('r','') if c.isalpha())
                    value = cell.findtext(NS+'v')
                    if cell.get('t') == 's' and value is not None:
                        value = strings[int(value)]
                    values[col] = value
                yield int(row.get('r')), values
                row.clear()
                while row.getprevious() is not None:
                    del row.getparent()[0]


def number(value):
    try:
        n = float(value)
        return n if math.isfinite(n) else None
    except (ValueError, TypeError):
        return None


def build():
    DB.parent.mkdir(exist_ok=True)
    temporary = DB.with_suffix('.building')
    temporary.unlink(missing_ok=True)
    con = sqlite3.connect(temporary)
    con.executescript('CREATE TABLE measurements(assay TEXT, symbol TEXT, sample TEXT, species TEXT, common TEXT, raw REAL, log REAL, normalized REAL, source TEXT, source_row INTEGER, norm_on TEXT); CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);')
    metadata = {}
    for assay, name in [('rna','CrossSpecies_Count.xlsx'),('protein','protLongDF.xlsx')]:
        path=SOURCE/name; started=time.time(); batch=[]; liver=0; total=0
        iterator=rows(path); _, header=next(iterator); columns={v:k for k,v in header.items()}
        for row_number, row in iterator:
            total+=1
            value=lambda key: row.get(columns.get(key))
            if str(value('tissue')).strip().lower() == 'liver':
                symbol=value('symbol'); species=value('species'); sample=value('sample' if assay=='rna' else 'sample.y')
                # User-confirmed tissue correction: this RNA sample is Lung.
                if assay=='rna' and sample=='BMR_Turk_Lu': continue
                if not symbol or not species or not sample: raise ValueError(f'{name} row {row_number}: missing identifiers')
                batch.append((assay,symbol.upper(),sample,species,value('common') or species,
                              number(value('tpm' if assay=='rna' else 'quant')),number(value('log.tpm' if assay=='rna' else 'log.quant')),
                              number(value('quant.norm')),value('source') or name,row_number,value('norm.on')))
                liver+=1
            if len(batch)>=10000:
                con.executemany('INSERT INTO measurements VALUES (?,?,?,?,?,?,?,?,?,?,?)',batch);con.commit();batch=[]
            if total%250000==0: print(f'{name}: {total:,} rows scanned, {liver:,} liver records',flush=True)
        if batch: con.executemany('INSERT INTO measurements VALUES (?,?,?,?,?,?,?,?,?,?,?)',batch)
        metadata[assay]=dict(file=name,rows=total,liver_rows=liver,size=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns)
        con.commit();print(f'{name}: complete in {time.time()-started:.1f}s; {liver:,} liver records',flush=True)
    mouse_path=SOURCE/'rnaLongTPM.xlsx'
    def mouse_records():
        iterator=rows(mouse_path); _,header=next(iterator); columns={v:k for k,v in header.items()}
        for n,row in iterator:
            value=lambda key:row.get(columns.get(key))
            if value('species')=='Mus musculus' and str(value('tissue')).lower()=='liver':
                yield value('symbol'),value('species'),value('sample'),number(value('tpm')),number(value('log.tpm')),n
    metadata['rna_mouse']=add_mouse_records(con,mouse_records(),mouse_path)
    con.execute('CREATE INDEX gene_lookup ON measurements(symbol,assay,species)')
    con.execute('CREATE INDEX assay_species ON measurements(assay,species)')
    con.execute('INSERT INTO metadata VALUES (?,?)',('sources',json.dumps(metadata)))
    con.commit();con.close();temporary.replace(DB)
    print(DB,flush=True)

if __name__=='__main__': build()
