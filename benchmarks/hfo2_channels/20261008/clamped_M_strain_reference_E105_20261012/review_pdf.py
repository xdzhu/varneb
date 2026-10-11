"""Render the actual E105 reading build; visual acceptance remains manual."""
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
        raise FileExistsError('fresh review directory required')
    args.output.mkdir(parents=True)
    pdf = pdfium.PdfDocument(args.pdf)
    previews, texts, details = [], [], []
    for index in range(len(pdf)):
        page = pdf[index]
        textpage = page.get_textpage()
        text = textpage.get_text_range()
        texts.append(text)
        preview = page.render(scale=70/72).to_pil().convert('RGB')
        previews.append(preview)
        preview.save(args.output/f'page-{index+1:02d}.png')
        if index >= len(pdf)-4 or any(marker in text for marker in
                ('Table 5.', 'Table 6.', 'Table 7.', '3.11.')):
            page.render(scale=130/72).to_pil().save(args.output/f'detail-{index+1:02d}.png')
            details.append(index+1)
        textpage.close()
        page.close()
    for start in range(0, len(previews), 9):
        sheet = Image.new('RGB', (900, 1314), 'white')
        draw = ImageDraw.Draw(sheet)
        for j, original in enumerate(previews[start:start+9]):
            thumb = original.copy()
            thumb.thumbnail((288, 405))
            x, y = j%3*300+6, j//3*438+24
            sheet.paste(thumb, (x, y))
            draw.text((x, y-18), f'Page {start+j+1}', fill='black')
        sheet.save(args.output/f'contact-{start//9+1}.png')
    full = '\n'.join(texts)
    numbers = ('97.227719', '211.974692', '-114.746973', '54.9251', '3.9008',
               '91.5873', '75.922579', '21.305139', '2.658882', '-1.230403')
    result = dict(pdf_sha256=hashlib.sha256(args.pdf.read_bytes()).hexdigest(),
        pdf_bytes=args.pdf.stat().st_size, pages=len(previews),
        rendered_pages=list(range(1, len(previews)+1)), detail_pages=details,
        backend='Python / existing pypdfium2',
        required_numbers_in_actual_pdf={n:n in full for n in numbers},
        figure9_caption_pages=[i+1 for i,t in enumerate(texts) if 'Figure 9.' in t],
        visual_review='pending actual image inspection')
    pdf.close()
    print(json.dumps(result))


if __name__ == '__main__':
    main()
