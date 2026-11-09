"""Vector PDF export of the exact browser-rendered tracks, without remote assets."""
import io
import math
import re
import xml.etree.ElementTree as ET
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import landscape, A4
from reportlab.lib.colors import HexColor
from reportlab.graphics import renderPDF
from reportlab.graphics.shapes import Group, String
from reportlab.lib.utils import simpleSplit
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
             'font-family', 'text-anchor', 'textLength', 'lengthAdjust', 'data-gene', 'data-loop-track', 'points', 'rx', 'role', 'aria-label', 'style', 'transform'}
    for element in root.iter():
        if element.tag.split('}')[-1] not in tags:
            raise ValueError('Unsupported SVG element')
        for attr, value in list(element.attrib.items()):
            if attr not in attrs or 'url(' in value.lower() or 'href' in attr.lower():
                raise ValueError('Unsupported SVG attribute')
            if attr == 'transform':
                match = re.fullmatch(r'translate\(0,(-?\d+(?:\.\d+)?)\)', value)
                if not match or abs(float(match.group(1))) > 1_000_000:
                    raise ValueError('Unsupported SVG transform')
            if attr in {'style', 'role', 'aria-label', 'preserveAspectRatio', 'data-gene', 'data-loop-track'}:
                del element.attrib[attr]
        if element.text:
            element.text = clean_text(element.text)
    view = [float(v) for v in root.attrib.get('viewBox', '').split()]
    if len(view) != 4 or not all(math.isfinite(v) for v in view) or not 1 <= view[2] <= 2000 or not 1 <= view[3] <= 1000000:
        raise ValueError('Invalid SVG dimensions')
    root.set('width', str(view[2])); root.set('height', str(view[3]))
    return svg2rlg(io.BytesIO(ET.tostring(root)))


def slide_layout(panels, drawings, page_size):
    """Fit one vertical stack with a shared plot width and left-hand labels."""
    page_width, page_height = page_size
    margin, ruler_height, row_padding = 14, 16, 1.5
    label_width = page_width * .18
    plot_x = margin + label_width + 10
    plot_width = page_width - margin - plot_x
    body_top = page_height - margin - ruler_height
    body_height = body_top - margin
    # Reduce label size only when many panels require it. Plot geometry shares
    # one vertical compression factor; labels stay beside their own tracks.
    label_scale = 1
    while True:
        labels = []
        for panel in panels:
            lines = []
            fields = [(panel.get('name', ''), 8.5, True), (panel.get('sample', ''), 7, False)]
            coords = clean_text(panel.get('coordinates', ''))
            fields.extend((part, 6.5, False) for part in coords.split(':', 1))
            for value, size, bold in fields:
                size *= label_scale
                for line in simpleSplit(clean_text(value), 'Helvetica-Bold' if bold else 'Helvetica', size, label_width):
                    lines.append((line, size, bold))
            labels.append(lines)
        minimum_heights = [sum(size * 1.15 for _, size, _ in lines) + row_padding * 2 for lines in labels]
        if sum(minimum_heights) <= body_height or label_scale < .1:
            break
        label_scale *= .9

    def heights(scale):
        return [max(minimum, drawing.height * scale + row_padding * 2) if drawing else max(minimum, 12 * label_scale)
                for minimum, drawing in zip(minimum_heights, drawings)]

    total_drawing_height = sum(drawing.height for drawing in drawings if drawing)
    low, high = 0, body_height / total_drawing_height if total_drawing_height else 0
    for _ in range(60):
        mid = (low + high) / 2
        if sum(heights(mid)) <= body_height:
            low = mid
        else:
            high = mid
    return dict(plot_x=plot_x, plot_width=plot_width, body_top=body_top,
                body_height=body_height, ruler_height=ruler_height, margin=margin,
                row_padding=row_padding, labels=labels, row_heights=heights(low), y_scale=low)


def compact_slide_text(drawing, x_scale, y_scale, minimum_size=4.25):
    """Keep SVG text glyphs proportional when track geometry is compressed.

    Text anchors stay at their original genomic positions. A height cap tied to
    the 12-pixel gene-label offset keeps glyphs below the exon bars.
    """
    def visit(group, sx, sy):
        a, _, _, d, _, _ = group.transform
        sx, sy = sx * abs(a), sy * abs(d)
        for index, child in enumerate(group.contents):
            if isinstance(child, String) and child.fontSize and sx and sy:
                size = min(child.fontSize * sx, max(child.fontSize * sy, min(minimum_size, 12 * sy)))
                ax, ay = size / (child.fontSize * sx), size / (child.fontSize * sy)
                wrapper = Group(child)
                wrapper.transform = (ax, 0, 0, ay, child.x * (1-ax), child.y * (1-ay))
                group.contents[index] = wrapper
            elif isinstance(child, Group):
                visit(child, sx, sy)
    visit(drawing, x_scale, y_scale)


def export_slide(payload, panels, drawings, ruler, page_size):
    layout = slide_layout(panels, drawings, page_size)
    stream = io.BytesIO()
    width, height = page_size
    c = canvas.Canvas(stream, pagesize=page_size, pageCompression=1)
    c.setTitle('NSC Genome Browser - ' + clean_text(payload.get('gene', '')))
    c.setAuthor('Upstate Nathan Shock Center')

    def text(x, y, value, size=9, color='#212529', bold=False):
        c.setFillColor(HexColor(color))
        c.setFont('Helvetica-Bold' if bold else 'Helvetica', size)
        c.drawString(x, y, clean_text(value))

    # A white figure with one ruler and no title block maximizes slide area.
    c.setFillColor(HexColor('#ffffff')); c.rect(0, 0, width, height, fill=1, stroke=0)
    x_scale, y_scale = layout['plot_width'] / ruler.width, layout['ruler_height'] / ruler.height
    compact_slide_text(ruler, x_scale, y_scale, minimum_size=6)
    c.saveState(); c.translate(layout['plot_x'], layout['body_top']); c.scale(x_scale, y_scale)
    renderPDF.draw(ruler, c, 0, 0); c.restoreState()
    cursor = layout['body_top']
    for panel, drawing, lines, row_height in zip(panels, drawings, layout['labels'], layout['row_heights']):
        label_y = cursor - layout['row_padding']
        for line, size, bold in lines:
            label_y -= size
            text(layout['margin'], label_y, line, size, '#003f7f' if bold else '#6c757d', bold)
            label_y -= size * .15
        if drawing:
            x_scale, y_scale = layout['plot_width'] / drawing.width, layout['y_scale']
            compact_slide_text(drawing, x_scale, y_scale)
            c.saveState(); c.translate(layout['plot_x'], cursor-layout['row_padding']-drawing.height*y_scale)
            c.scale(x_scale, y_scale); renderPDF.draw(drawing, c, 0, 0); c.restoreState()
        else:
            text(layout['plot_x'], cursor-layout['row_padding']-7, panel.get('error', 'No data'), 7, '#6c757d')
        cursor -= row_height
        c.setStrokeColor(HexColor('#dee2e6')); c.setLineWidth(.4)
        c.line(layout['margin'], cursor, width-layout['margin'], cursor)
    c.save()
    return stream.getvalue()


def export_pdf(payload):
    panels = payload.get('panels', [])
    if not isinstance(panels, list) or not 1 <= len(panels) <= 40:
        raise ValueError('Choose between 1 and 40 panels')
    ruler = safe_drawing(payload['ruler'])
    drawings = [safe_drawing(panel['svg']) if panel.get('svg') else None for panel in panels]
    page_layout = payload.get('pageLayout', 'full-height')
    if page_layout in {'16:9', '4:3'}:
        return export_slide(payload, panels, drawings, ruler, (960 if page_layout == '16:9' else 720, 540))
    if page_layout != 'full-height':
        raise ValueError('Choose a 16:9, 4:3 or full-height PDF layout')
    stream = io.BytesIO()
    # One continuous vector figure: preserve each SVG's aspect ratio and grow
    # the page vertically instead of clipping tall tracks or paginating species.
    page_width, minimum_height = landscape(A4)
    margin, label_width = 28, 132
    plot_width = page_width - 2 * margin - label_width
    row_heights = [max(56, drawing.height * plot_width / drawing.width + 6) if drawing else 56 for drawing in drawings]
    header_height, footer_height = 129, 55
    page_height = max(minimum_height, header_height + sum(row_heights) + footer_height)
    c = canvas.Canvas(stream, pagesize=(page_width, page_height), pageCompression=1)
    c.setTitle('NSC Genome Browser - ' + clean_text(payload.get('gene', '')))
    c.setAuthor('Upstate Nathan Shock Center')

    def text(x, y, value, size=9, color='#212529', bold=False):
        c.setFillColor(HexColor(color)); c.setFont('Helvetica-Bold' if bold else 'Helvetica', size)
        c.drawString(x, y, clean_text(value))

    def footer():
        text(margin, 29, 'TSS / symbol matching; equal TSS distance does not establish sequence homology.', 7, '#6c757d')
        text(margin, 18, 'Coordinates: 1-based inclusive. Source signal units; displayed limits may clip peaks.', 7, '#6c757d')

    def start_page():
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
    for panel, drawing, height in zip(panels, drawings, row_heights):
        text(margin, cursor-10, panel.get('name',''), 9, '#003f7f', True)
        text(margin, cursor-22, panel.get('sample',''), 8)
        # Wrap long assembly coordinate labels without cutting accession strings.
        coords = clean_text(panel.get('coordinates', ''))
        if ':' in coords:
            chrom, interval = coords.split(':', 1)
            text(margin, cursor-34, chrom, 7, '#6c757d')
            text(margin, cursor-44, interval, 7, '#6c757d')
        else:
            text(margin, cursor-34, coords, 7, '#6c757d')
        if drawing:
            scale=plot_width/drawing.width
            c.saveState(); c.translate(margin+label_width, cursor-3-drawing.height*scale)
            c.scale(scale,scale)
            renderPDF.draw(drawing,c,0,0); c.restoreState()
        else:
            text(margin+label_width,cursor-30,panel.get('error','No data'),9,'#6c757d')
        c.setStrokeColor(HexColor('#dee2e6')); c.setLineWidth(.5)
        c.line(margin,cursor-height,page_width-margin,cursor-height)
        cursor -= height
    footer(); c.save()
    return stream.getvalue()
