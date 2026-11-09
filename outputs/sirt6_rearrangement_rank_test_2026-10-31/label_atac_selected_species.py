"""Label the top two rearranged and bottom four unrearranged SIRT6 ranks."""
from io import BytesIO
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

SOURCE = Path(__file__).resolve().parent
ROOT = SOURCE.parents[1]
OUT = ROOT / 'output/pdf/SIRT6_rank_comparison_selected_species.pdf'
OUT.parent.mkdir(parents=True, exist_ok=True)
data = json.loads((SOURCE / 'analysis.json').read_text())
svg = ET.parse(SOURCE / 'SIRT6_rank_comparison.svg').getroot()
ns = {'s': 'http://www.w3.org/2000/svg'}
reader = PdfReader(SOURCE / 'SIRT6_rank_comparison.pdf')
assert len(reader.pages) == 1
page = reader.pages[0]
width, height = float(page.mediabox.width), float(page.mediabox.height)
assert abs(width - 720) < 1e-6 and abs(height - 374.4) < 1e-6
pdfmetrics.registerFont(TTFont('ArialEmbedded', '/System/Library/Fonts/Supplemental/Arial.ttf'))
overlay = BytesIO()
c = canvas.Canvas(overlay, pagesize=(width, height), initialFontName='ArialEmbedded')
labels = []

for group, collection, color in [('Rearranged', 'PathCollection_3', '#0072B2'),
                                  ('Unrearranged', 'PathCollection_4', '#9A6A00')]:
    rows = sorted([r for r in data['measurements']
                   if r['gene'] == 'SIRT6' and r['assay'] == 'atac'
                   and r['used_in_test'] and r['group'] == group], key=lambda r: r['rank'])
    points = svg.find(f'.//s:g[@id="{collection}"]', ns).findall('.//s:use', ns)
    assert len(rows) == len(points)
    # Smaller rank means higher relative activity, matching the inverted plot axis.
    chosen = rows[:2] if group == 'Rearranged' else rows[-4:]
    selected_ids = {r['species_id'] for r in chosen}
    for r, point in zip(rows, points):
        if r['species_id'] not in selected_ids:
            continue
        x, svg_y = float(point.attrib['x']), float(point.attrib['y'])
        expected_y = 61.729912 + (r['rank'] - 1) * (332.357495 - 61.729912) / 34
        assert abs(svg_y - expected_y) < 1e-5
        y = height - svg_y
        right_x, baseline = x - 13, y + 8
        if r['species'] == 'African spiny mouse':
            right_x, baseline = 500, y + 10
        elif r['species'] == 'Jamaican fruit bat':
            right_x = 600
        text_width = pdfmetrics.stringWidth(r['species'], 'ArialEmbedded', 8.5)
        assert right_x - text_width > 403 and right_x < 717
        c.setStrokeColor(HexColor(color))
        c.setFillColor(HexColor(color))
        c.setLineWidth(0.65)
        line_end_x = right_x + 3
        if r['species'] == 'African spiny mouse':
            line_end_x = x + 4
        c.line(x - 3.5, y + 1.5, line_end_x, baseline - 2)
        c.setFont('ArialEmbedded', 8.5)
        c.drawRightString(right_x, baseline, r['species'])
        labels.append(dict(species=r['species'], group=group, rank=r['rank']))

assert [r['species'] for r in labels] == ['Meadow vole', 'African spiny mouse',
    'Jamaican fruit bat', 'Nutria', 'North American beaver', 'Cave nectar bat']
c.save()
page.merge_page(PdfReader(overlay).pages[0])
writer = PdfWriter()
writer.add_page(page)
writer.add_metadata({'/Title': 'SIRT6 rank comparison - selected ATAC species',
                     '/Subject': 'Top two rearranged and bottom four unrearranged species by relative SIRT6 activity; original data and statistics retained.'})
with OUT.open('wb') as f:
    writer.write(f)
check = PdfReader(OUT)
text = check.pages[0].extract_text()
assert all(r['species'] in text for r in labels)
assert len(check.pages) == 1
print(json.dumps({'pdf': str(OUT), 'labeled': labels}, indent=2))
