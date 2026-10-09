import io
import unittest
import server  # Enables optional local PDF dependencies.
from pdf_export import export_pdf, safe_drawing
from pypdf import PdfReader

SVG = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 900 179"><line x1="0" x2="900" y1="67" y2="67" stroke="#003f7f"/><text x="700" y="10" font-size="8">0.00-0.500 CPM</text></svg>'

class PDFExport(unittest.TestCase):
    def test_pagination_retains_samples_and_limits(self):
        panels = [dict(name='Mouse', sample=f'sample-{i}', coordinates='chr1:1-10000 (-)', svg=SVG) for i in range(8)]
        blob = export_pdf(dict(gene='HAS2', summary='8 panels', ruler=SVG, panels=panels))
        reader = PdfReader(io.BytesIO(blob))
        self.assertEqual(len(reader.pages),2)
        for page_index,page in enumerate(reader.pages):
            page_text=page.extract_text()
            self.assertEqual(sum(f'sample-{i}' in page_text for i in range(8)),4)
            for i in range(page_index*4,page_index*4+4):self.assertIn(f'sample-{i}',page_text)
        text = '\n'.join(p.extract_text() for p in reader.pages)
        for i in range(8): self.assertIn(f'sample-{i}', text)
        self.assertEqual(text.count('0.00-0.500 CPM'),8+len(reader.pages))
        self.assertIn('UPSTATE NATHAN SHOCK CENTER',text)

    def test_tall_panels_and_last_page(self):
        panels=[dict(name=f'Species {i}',sample=f'replicate-{i}',svg=SVG.replace('900 179','900 310')) for i in range(9)]
        reader=PdfReader(io.BytesIO(export_pdf(dict(gene='SOX2',ruler=SVG,panels=panels))))
        self.assertEqual(len(reader.pages),3)
        self.assertEqual([p.extract_text().count('replicate-') for p in reader.pages],[4,4,1])

    def test_rejects_external_svg_content(self):
        for svg in ['<!DOCTYPE svg><svg/>', '<svg><image href="https://example.com"/></svg>', '<svg><script/></svg>']:
            with self.assertRaises(ValueError): safe_drawing(svg)

    def test_empty_view_rejected(self):
        with self.assertRaises(ValueError): export_pdf({'panels':[]})
