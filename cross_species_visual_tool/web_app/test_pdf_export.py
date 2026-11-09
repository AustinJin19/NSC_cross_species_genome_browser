import io
import unittest
import server  # Enables optional local PDF dependencies.
from pdf_export import export_pdf, safe_drawing, slide_layout, compact_slide_text
from reportlab.graphics.shapes import Drawing, String
from pypdf import PdfReader

SVG = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 900 179"><line x1="0" x2="900" y1="67" y2="67" stroke="#003f7f"/><text x="700" y="10" font-size="8">0.00-0.500 CPM</text></svg>'

class PDFExport(unittest.TestCase):
    def test_single_page_retains_all_samples_and_limits(self):
        panels = [dict(name='Mouse', sample=f'sample-{i}', coordinates='chr1:1-10000 (-)', svg=SVG) for i in range(8)]
        blob = export_pdf(dict(gene='HAS2', summary='8 panels', ruler=SVG, panels=panels))
        reader = PdfReader(io.BytesIO(blob))
        self.assertEqual(len(reader.pages),1)
        text = '\n'.join(p.extract_text() for p in reader.pages)
        for i in range(8): self.assertIn(f'sample-{i}', text)
        self.assertEqual(text.count('0.00-0.500 CPM'),8+len(reader.pages))
        self.assertIn('UPSTATE NATHAN SHOCK CENTER',text)

    def test_tall_panels_expand_page_without_squeezing(self):
        panels=[dict(name=f'Species {i}',sample=f'replicate-{i}',svg=SVG.replace('900 179','900 1000')) for i in range(9)]
        reader=PdfReader(io.BytesIO(export_pdf(dict(gene='SOX2',ruler=SVG,panels=panels))))
        self.assertEqual(len(reader.pages),1)
        self.assertEqual(reader.pages[0].extract_text().count('replicate-'),9)
        self.assertGreater(float(reader.pages[0].mediabox.height), 6500)
        from pypdf.generic import ContentStream
        transforms = [operands for operands, op in ContentStream(reader.pages[0].get_contents(), reader).operations if op == b'cm']
        for a, b, c, d, e, f in transforms:
            self.assertAlmostEqual(abs(float(a)), abs(float(d)))

    def test_loop_metadata_exports_but_unsafe_attributes_still_rejected(self):
        svg = SVG.replace('</svg>', '<path data-loop-track="human-liver" d="M 10 120 Q 100 20 200 120" fill="none" stroke="#1864ab" stroke-dasharray="4 3"/></svg>')
        reader = PdfReader(io.BytesIO(export_pdf({'ruler': SVG, 'panels': [{'name': 'Human', 'svg': svg}]})))
        self.assertEqual(len(reader.pages), 1)
        self.assertIn('Human', reader.pages[0].extract_text())
        for bad in [svg.replace('data-loop-track', 'onclick'), svg.replace('data-loop-track', 'href'), svg.replace('#1864ab', 'url(https://example.com)')]:
            with self.assertRaises(ValueError): safe_drawing(bad)

    def test_rejects_external_svg_content(self):
        for svg in ['<!DOCTYPE svg><svg/>', '<svg><image href="https://example.com"/></svg>', '<svg><script/></svg>']:
            with self.assertRaises(ValueError): safe_drawing(svg)

    def test_spaced_gene_labels_export(self):
        svg = SVG.replace('font-size="8"', 'font-size="8" data-gene="APOE" textLength="40" lengthAdjust="spacingAndGlyphs"')
        reader = PdfReader(io.BytesIO(export_pdf({'ruler': SVG, 'panels': [{'name': 'Human', 'svg': svg}]})))
        self.assertIn('0.00-0.500 CPM', reader.pages[0].extract_text())

    def test_empty_view_rejected(self):
        with self.assertRaises(ValueError): export_pdf({'panels':[]})

    def test_compact_track_translation_exports_safely(self):
        svg = SVG.replace('font-size="8"', 'font-size="8" transform="translate(0,-2.500)"')
        reader = PdfReader(io.BytesIO(export_pdf({'ruler': SVG, 'panels': [{'name': 'Human', 'svg': svg}]})))
        self.assertIn('0.00-0.500 CPM', reader.pages[0].extract_text())
        for transform in ['rotate(30)', 'scale(2)', 'translate(10,-2)', 'translate(0,NaN)', 'translate(0,-1000001)']:
            with self.assertRaises(ValueError): safe_drawing(svg.replace('translate(0,-2.500)', transform))

    def test_slide_formats_retain_five_species_on_one_page(self):
        heights = [448, 330, 230, 252, 418]
        panels = [dict(name=f'Species-{i}', sample=f'sample-{i}', coordinates='chr19:2951598-4551598 (-)',
                       svg=SVG.replace('900 179', f'900 {h}')) for i, h in enumerate(heights)]
        drawings = [safe_drawing(p['svg']) for p in panels]
        for page_layout, width in [('16:9', 960), ('4:3', 720)]:
            layout = slide_layout(panels, drawings, (width, 540))
            self.assertGreater(layout['plot_width'], width * .75)
            self.assertAlmostEqual(sum(layout['row_heights']), layout['body_height'])
            self.assertGreater(layout['y_scale'], 0)
            for drawing, row_height in zip(drawings, layout['row_heights']):
                self.assertLessEqual(drawing.height * layout['y_scale'] + layout['row_padding'] * 2, row_height + .001)
            reader = PdfReader(io.BytesIO(export_pdf(dict(gene='SIRT6', ruler=SVG, panels=panels, pageLayout=page_layout))))
            self.assertEqual(len(reader.pages), 1)
            self.assertEqual(tuple(float(v) for v in reader.pages[0].mediabox[2:]), (width, 540))
            text = reader.pages[0].extract_text()
            for i in range(5):
                self.assertEqual(text.count(f'Species-{i}'), 1)
                self.assertEqual(text.count(f'sample-{i}'), 1)
            self.assertEqual([text.index(f'Species-{i}') for i in range(5)], sorted(text.index(f'Species-{i}') for i in range(5)))
            self.assertNotIn('NSC Genome Browser |', text)
            # One ruler, and all five plot origins share the same horizontal position.
            from pypdf.generic import ContentStream
            origins = [operands for operands, op in ContentStream(reader.pages[0].get_contents(), reader).operations
                       if op == b'cm' and abs(float(operands[4]) - layout['plot_x']) < .001]
            self.assertEqual(len(origins), 6)
            self.assertEqual([float(v[5]) for v in origins], sorted((float(v[5]) for v in origins), reverse=True))

    def test_compact_text_preserves_anchor_and_glyph_proportions(self):
        drawing = Drawing(900, 400)
        label = String(200, 120, 'SIRT6', fontSize=8)
        drawing.add(label)
        compact_slide_text(drawing, .8, .25)
        wrapper = drawing.contents[0]
        a, _, _, d, e, f = wrapper.transform
        self.assertAlmostEqual(a * label.x + e, label.x)
        self.assertAlmostEqual(d * label.y + f, label.y)
        self.assertAlmostEqual(a * .8, d * .25)
        self.assertGreater(d * .25 * label.fontSize, 2)
        self.assertIs(wrapper.contents[0], label)

    def test_slide_missing_track_and_invalid_layout(self):
        payload = dict(gene='SIRT6', ruler=SVG, panels=[dict(name='Missing', error='No gene annotation')], pageLayout='4:3')
        reader = PdfReader(io.BytesIO(export_pdf(payload)))
        self.assertIn('No gene annotation', reader.pages[0].extract_text())
        payload['pageLayout'] = 'invalid'
        with self.assertRaises(ValueError): export_pdf(payload)
