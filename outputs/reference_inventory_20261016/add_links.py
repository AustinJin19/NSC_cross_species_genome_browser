"""Add native Excel links, a capability missing from the authoring API.

All workbook contents and formatting are authored by artifact-tool. This
bounded OOXML step only supplies hyperlink relationships for resource cells.
"""
import json, pathlib, zipfile, xml.etree.ElementTree as ET

root=pathlib.Path(__file__).resolve().parent
path=root/'NSC_genome_assemblies_and_annotations.xlsx'
data=json.loads((root/'inventory.json').read_text())
ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
rel='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
pkg='http://schemas.openxmlformats.org/package/2006/relationships'
ET.register_namespace('',ns)
ET.register_namespace('r',rel)
with zipfile.ZipFile(path) as z:
    items={x.filename:z.read(x.filename) for x in z.infolist()}
sheet=ET.fromstring(items['xl/worksheets/sheet1.xml'])
linkpath='xl/worksheets/_rels/sheet1.xml.rels'
rels=ET.fromstring(items[linkpath]) if linkpath in items else ET.Element(f'{{{pkg}}}Relationships')
links=ET.Element(f'{{{ns}}}hyperlinks')
count=0
for i,s in enumerate(data['species'],7):
    for col,key in [('H','assembly_url'),('I','annotation_url')]:
        url=s[key]
        if not url: continue
        count+=1
        rid=f'rIdReferenceLink{count}'
        ET.SubElement(rels,f'{{{pkg}}}Relationship',Id=rid,Type=rel+'/hyperlink',Target=url,TargetMode='External')
        ET.SubElement(links,f'{{{ns}}}hyperlink',{'ref':f'{col}{i}',f'{{{rel}}}id':rid})
# Place hyperlinks after data validations and before print/drawing/table nodes.
later={'printOptions','pageMargins','pageSetup','headerFooter','rowBreaks','colBreaks','customProperties','cellWatches','ignoredErrors','smartTags','drawing','legacyDrawing','legacyDrawingHF','picture','oleObjects','controls','webPublishItems','tableParts','extLst'}
pos=next((i for i,e in enumerate(sheet) if e.tag.split('}')[-1] in later),len(sheet))
sheet.insert(pos,links)
items['xl/worksheets/sheet1.xml']=ET.tostring(sheet,encoding='utf-8',xml_declaration=True)
items[linkpath]=ET.tostring(rels,encoding='utf-8',xml_declaration=True)
tmp=path.with_suffix('.xlsx.tmp')
with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as z:
    for name,value in items.items(): z.writestr(name,value)
tmp.replace(path)
with zipfile.ZipFile(path) as z:
    assert z.testzip() is None
    for name in z.namelist():
        if name.endswith('.xml') or name.endswith('.rels'): ET.fromstring(z.read(name))
    sh=ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
    assert len(sh.find(f'{{{ns}}}hyperlinks'))==66
    assert not sh.findall(f'.//{{{ns}}}f')
    assert not sh.findall(f'.//{{{ns}}}c[@t="e"]')
    print('Native resource hyperlinks:',count)
    print('Frozen panes:',[x.attrib for x in sh.findall(f'.//{{{ns}}}pane')])
    print('Workbook tables:',[n for n in z.namelist() if n.startswith('xl/tables/') and n.endswith('.xml')])
    print('Workbook bytes:',path.stat().st_size)
