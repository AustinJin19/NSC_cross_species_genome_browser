"""Vector PDF export of the exact browser-rendered tracks, without remote assets."""
import io
import math
import xml.etree.ElementTree as ET
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import landscape, A4
from reportlab.lib.colors import HexColor
from reportlab.graphics import renderPDF
from svglib.svglib import svg2rlg


def clean_text(value, limit=300):
    return str(value)[:limit].replace('\u2013', '-').replace('\u2212', '-').replace('\u2192', '>').replace('\u2190', '<').replace('\u2032', "'")


def safe_drawing(svg):
    # Only the small SVG vocabulary emitted by browser.js is accepted. No URLs,
    # images, stylesheets, scripts, external entities, or arbitrary SVG documents.
    if not isinstance(svg, str) or len(svg) > 2_000_000 or '<!' in svg or '<?' in svg:
        raise ValueError('Invalid track SVG')
    root = ET.fromstring(svg)
    tags = {'svg', 'line', 'text', 'polygon', 'rect', 'path'}
    attrs = {'viewBox', 'preserveAspectRatio', 'width', 'height', 'x', 'y', 'x1', 'x2', 'y1', 'y2',
             'd', 'stroke-opacity', 'fill', 'fill-opacity', 'stroke', 'stroke-width', 'stroke-dasharray', 'font-size',
             'font-family', 'text-anchor', 'points', 'rx', 'role', 'aria-label', 'style'}
    for element in root.iter():
        if element.tag.split('}')[-1] not in tags:
            raise ValueError('Unsupported SVG element')
        for attr, value in list(element.attrib.items()):
            if attr not in attrs or 'url(' in value.lower() or 'href' in attr.lower():
                raise ValueError('Unsupported SVG attribute')
            if attr in {'style', 'role', 'aria-label', 'preserveAspectRatio'}:
                del element.attrib[attr]
        if element.text:
            element.text = clean_text(element.text)
    view = [float(v) for v in root.attrib.get('viewBox', '').split()]
    if len(view) != 4 or not all(math.isfinite(v) for v in view) or not 1 <= view[2] <= 2000 or not 1 <= view[3] <= 2000:
        raise ValueError('Invalid SVG dimensions')
    root.set('width', str(view[2])); root.set('height', str(view[3]))
    return svg2rlg(io.BytesIO(ET.tostring(root)))


def export_pdf(payload):
    panels = payload.get('panels', [])
    if not isinstance(panels, list) or not 1 <= len(panels) <= 40:
        raise ValueError('Choose between 1 and 40 panels')
    ruler = safe_drawing(payload['ruler'])
    stream = io.BytesIO()
    page_width, page_height = A4
    c = canvas.Canvas(stream, pagesize=(page_width, page_height), pageCompression=1)
    c.setTitle('NSC Genome Browser - ' + clean_text(payload.get('gene', '')))
    c.setAuthor('Upstate Nathan Shock Center')
    margin, label_width = 28, 120
    plot_width = page_width - 2 * margin - label_width
    page_number = 0

    def text(x, y, value, size=9, color='#212529', bold=False):
        c.setFillColor(HexColor(color)); c.setFont('Helvetica-Bold' if bold else 'Helvetica', size)
        c.drawString(x, y, clean_text(value))

    def footer():
        text(margin, 29, 'TSS / symbol matching; equal TSS distance does not establish sequence homology.', 7, '#6c757d')
        text(margin, 18, 'Coordinates: 1-based inclusive. Source signal units; displayed limits may clip peaks.', 7, '#6c757d')
        c.setFont('Helvetica', 8); c.drawRightString(page_width-margin, 20, str(page_number))

    def start_page():
        nonlocal page_number
        page_number += 1
        text(margin, page_height-35, 'UPSTATE NATHAN SHOCK CENTER', 10, '#003f7f', True)
        text(margin, page_height-60, 'NSC Genome Browser | ' + clean_text(payload.get('gene', '')), 18, '#003f7f', True)
        text(margin, page_height-79, payload.get('summary', ''), 9, '#6c757d')
        c.setStrokeColor(HexColor('#ffd100')); c.setLineWidth(2)
        c.line(margin, page_height-90, page_width-margin, page_height-90)
        scale = plot_width / ruler.width
        c.saveState(); c.translate(margin+label_width, page_height-120)
        c.scale(scale, scale); renderPDF.draw(ruler,c,0,0); c.restoreState()
        return page_height-129

    cursor = start_page()
    slot_height = (page_height - 129 - 55) / 4
    for panel_index, panel in enumerate(panels):
        drawing = safe_drawing(panel['svg']) if panel.get('svg') else None
        height = slot_height - 10
        if panel_index and panel_index % 4 == 0:
            footer(); c.showPage(); cursor=start_page()
        text(margin, cursor-14, panel.get('name',''), 9, '#003f7f', True)
        text(margin, cursor-30, panel.get('sample',''), 8)
        # Wrap long assembly coordinate labels without cutting accession strings.
        coords = clean_text(panel.get('coordinates', ''))
        if ':' in coords:
            chrom, interval = coords.split(':', 1)
            text(margin, cursor-46, chrom, 7, '#6c757d')
            text(margin, cursor-58, interval, 7, '#6c757d')
        else:
            text(margin, cursor-46, coords, 7, '#6c757d')
        if drawing:
            c.saveState(); c.translate(margin+label_width, cursor-height+10)
            scale=plot_width/drawing.width
            c.scale(scale,min(scale,(height-18)/drawing.height))
            renderPDF.draw(drawing,c,0,0); c.restoreState()
        else:
            text(margin+label_width,cursor-30,panel.get('error','No data'),9,'#6c757d')
        c.setStrokeColor(HexColor('#dee2e6')); c.setLineWidth(.5)
        c.line(margin,cursor-height,page_width-margin,cursor-height)
        cursor -= height+10
    footer(); c.save()
    return stream.getvalue()
