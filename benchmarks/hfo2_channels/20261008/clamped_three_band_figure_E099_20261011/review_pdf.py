"""Render the actual reading build with existing Python PDFium; no DFT."""
import argparse
import hashlib
import json
from pathlib import Path

import pypdfium2 as pdfium
from PIL import Image, ImageDraw


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdf', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh PDF review directory required')
    args.output.mkdir(parents=True)
    pdf = pdfium.PdfDocument(args.pdf)
    texts, previews, detail_pages = [], [], []
    for i in range(len(pdf)):
        page = pdf[i]
        text_page = page.get_textpage()
        text = text_page.get_text_range()
        texts.append(text)
        preview = page.render(scale=70/72).to_pil().convert('RGB')
        preview.save(args.output/f'page-{i+1:02d}.png')
        previews.append(preview)
        if 'Three ordinary-converged' in text or 'Matched-substrate' in text or '(d--f)' in text:
            page.render(scale=130/72).to_pil().save(args.output/f'detail-{i+1:02d}.png')
            detail_pages.append(i+1)
        text_page.close(); page.close()
    for start in range(0, len(previews), 9):
        sheet = Image.new('RGB', (900, 1314), 'white')
        draw = ImageDraw.Draw(sheet)
        for j, original in enumerate(previews[start:start+9]):
            thumb = original.copy(); thumb.thumbnail((288, 405))
            x, y = (j%3)*300+6, (j//3)*438+24
            sheet.paste(thumb, (x, y)); draw.text((x, y-18), f'Page {start+j+1}', fill='black')
        sheet.save(args.output/f'contact-{start//9+1}.png')
    result = dict(pdf_sha256=hashlib.sha256(args.pdf.read_bytes()).hexdigest(),
        pdf_bytes=args.pdf.stat().st_size, pages=len(previews),
        rendered_pages=list(range(1, len(previews)+1)), detail_pages=detail_pages,
        backend='Python / existing pypdfium2', figure5_caption_found=bool(detail_pages),
        visual_review='pending actual image inspection')
    print(json.dumps(result))
    pdf.close()


if __name__ == '__main__':
    main()
