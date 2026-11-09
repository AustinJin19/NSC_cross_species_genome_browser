"""Label Tukey outliers on the original PDF without changing its plotted data."""
from io import BytesIO
import json
from pathlib import Path
from statistics import quantiles
import xml.etree.ElementTree as ET

from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

SOURCE = Path(__file__).resolve().parent
ROOT = SOURCE.parents[1]
OUT = ROOT / 'output/pdf/SIRT6_rank_comparison_outliers_labeled.pdf'
OUT.parent.mkdir(parents=True, exist_ok=True)
data = json.loads((SOURCE / 'analysis.json').read_text())
svg = ET.parse(SOURCE / 'SIRT6_rank_comparison.svg').getroot()
ns = {'s': 'http://www.w3.org/2000/svg'}
reader = PdfReader(SOURCE / 'SIRT6_rank_comparison.pdf')
assert len(reader.pages) == 1
page = reader.pages[0]
width, height = float(page.mediabox.width), float(page.mediabox.height)
assert abs(width - 720) < 1e-6 and abs(height - 374.4) < 1e-6
overlay = BytesIO()
pdfmetrics.registerFont(TTFont('ArialEmbedded', '/System/Library/Fonts/Supplemental/Arial.ttf'))
c = canvas.Canvas(overlay, pagesize=(width, height), initialFontName='ArialEmbedded')
audit = []
labels = []

for group, collection, color in [('Rearranged', 'PathCollection_3', '#0072B2'),
                                  ('Unrearranged', 'PathCollection_4', '#D59B00')]:
    rows = sorted([r for r in data['measurements']
                   if r['gene'] == 'SIRT6' and r['assay'] == 'atac'
                   and r['used_in_test'] and r['group'] == group], key=lambda r: r['rank'])
    points = svg.find(f'.//s:g[@id="{collection}"]', ns).findall('.//s:use', ns)
    assert len(rows) == len(points)
    q1, _, q3 = quantiles([r['rank'] for r in rows], n=4, method='inclusive')
    lower, upper = q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1)
    for r, point in zip(rows, points):
        x, svg_y = float(point.attrib['x']), float(point.attrib['y'])
        # Original SVG and PDF were exported together. Confirm plotted rank mapping.
        expected_y = 61.729912 + (r['rank'] - 1) * (332.357495 - 61.729912) / 34
        assert abs(svg_y - expected_y) < 1e-5
        is_outlier = r['rank'] < lower or r['rank'] > upper
        audit.append(dict(species=r['species'], group=group, rank=r['rank'],
                          q1=q1, q3=q3, lower_fence=lower, upper_fence=upper,
                          outlier=is_outlier))
        if not is_outlier:
            continue
        y = height - svg_y
        label_x, label_y = x + 15, y + 8
        c.setStrokeColor(HexColor(color))
        c.setLineWidth(0.65)
        c.line(x + 3.5, y + 1.5, label_x - 3, label_y - 1)
        c.setFillColor(HexColor(color))
        c.setFont('ArialEmbedded', 9.5)
        c.drawString(label_x, label_y, r['species'])
        labels.append(r['species'])

assert labels == ['Meadow vole'], labels
c.save()
page.merge_page(PdfReader(overlay).pages[0])
writer = PdfWriter()
writer.add_page(page)
writer.add_metadata({'/Title': 'SIRT6 rank comparison - ATAC outliers labeled',
                     '/Subject': 'Original provisional groups and statistics preserved; ATAC outliers use 1.5 IQR fences.'})
with OUT.open('wb') as f:
    writer.write(f)
check = PdfReader(OUT)
assert len(check.pages) == 1
assert 'Meadow vole' in check.pages[0].extract_text()
assert float(check.pages[0].mediabox.width) == width
assert float(check.pages[0].mediabox.height) == height
audit_path = ROOT / 'tmp/pdfs/sirt6_outlier_labels/outlier_audit.json'
audit_path.parent.mkdir(parents=True, exist_ok=True)
audit_path.write_text(json.dumps(audit, indent=2) + '\n')
print(json.dumps({'pdf': str(OUT), 'labeled': labels,
                  'unrearranged_outliers': sum(r['outlier'] for r in audit if r['group'] == 'Unrearranged')}, indent=2))
