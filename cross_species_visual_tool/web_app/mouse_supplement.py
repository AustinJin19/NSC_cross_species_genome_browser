"""Add missing mouse liver records while keeping primary RNA measurements."""
import json

def add_mouse_records(con, records, path):
    if 'workbook' not in {r[1] for r in con.execute('PRAGMA table_info(measurements)')}:
        con.execute('ALTER TABLE measurements ADD COLUMN workbook TEXT')
        con.execute("UPDATE measurements SET workbook=CASE assay WHEN 'rna' THEN 'CrossSpecies_Count.xlsx' ELSE 'protLongDF.xlsx' END")
    keys=set(con.execute("SELECT symbol,species,sample FROM measurements WHERE assay='rna'"))
    batch=[]
    for symbol,species,sample,raw,log,source_row in records:
        key=(symbol.upper(),species,sample)
        if species!='Mus musculus' or key in keys: continue
        keys.add(key)
        batch.append(('rna',symbol.upper(),sample,species,species,raw,log,None,path.name,source_row,None,path.name))
    con.executemany('INSERT INTO measurements VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',batch)
    return dict(file=path.name,liver_rows=len(batch),scope='Mouse liver supplement only',size=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns)
