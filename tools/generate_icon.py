#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "figure"
FONT = "/usr/share/fonts/liberation-mono/LiberationMono-Bold.ttf"

ASCII_LOGO = [
    r"  ___   ____    ___    ____  _   _  ____  _  __  ___  _____ ",
    r" / _ \ | __ )  / _ \  / ___|| | | |/ ___|| |/ / |_ _||_   _|",
    r"| /_\ ||  _ \ | /_\ || |    | | | |\___ \| ' /   | |   | |  ",
    r"|  _  || |_) ||  _  || |___ | |_| | ___) | . \   | |   | |  ",
    r"|_| |_||____/ |_| |_| \____| \___/ |____/|_|\_\ |___|  |_|  ",
]


def load_font(size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(FONT, size=size)


def draw_centered_multiline(
    image: Image.Image,
    lines: list[str],
    font: ImageFont.FreeTypeFont,
    fill: tuple[int, int, int, int],
    y: int,
    x_offset: int = 0,
) -> None:
    draw = ImageDraw.Draw(image)
    line_height = int(font.size * 1.07)
    for row, line in enumerate(lines):
        bbox = draw.textbbox((0, 0), line, font=font)
        width = bbox[2] - bbox[0]
        x = (image.width - width) // 2 + x_offset
        draw.text((x, y + row * line_height), line, font=font, fill=fill)


def add_scanlines(image: Image.Image, alpha: int = 26) -> None:
    overlay = Image.new("RGBA", image.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    for y in range(0, image.height, 6):
        draw.line((0, y, image.width, y), fill=(255, 255, 255, alpha), width=1)
    image.alpha_composite(overlay)


def make_banner() -> None:
    width, height = 1360, 238
    bg = (8, 10, 12, 255)
    image = Image.new("RGBA", (width, height), bg)
    font = load_font(31)
    subtitle_font = load_font(18)

    logo_y = 20
    glow = Image.new("RGBA", image.size, (0, 0, 0, 0))
    draw_centered_multiline(glow, ASCII_LOGO, font, (38, 218, 255, 160), logo_y, x_offset=2)
    glow = glow.filter(ImageFilter.GaussianBlur(3.2))
    image.alpha_composite(glow)

    draw_centered_multiline(image, ASCII_LOGO, font, (217, 222, 226, 255), logo_y)
    draw_centered_multiline(image, ASCII_LOGO, font, (94, 240, 255, 92), logo_y, x_offset=2)

    draw = ImageDraw.Draw(image)
    subtitle = "ABACUS + DeepMD workflow toolkit  |  v1.0  |  Han Enci, Xi'an University of Technology"
    bbox = draw.textbbox((0, 0), subtitle, font=subtitle_font)
    draw.text(
        ((width - (bbox[2] - bbox[0])) // 2, 200),
        subtitle,
        font=subtitle_font,
        fill=(118, 205, 224, 230),
    )
    add_scanlines(image)
    image.save(FIGURE_DIR / "abacuskit_logo.png")


def make_square_icon() -> None:
    size = 512
    image = Image.new("RGBA", (size, size), (8, 10, 12, 255))
    draw = ImageDraw.Draw(image)
    font_big = load_font(138)
    font_small = load_font(25)

    # Soft terminal glow.
    glow = Image.new("RGBA", image.size, (0, 0, 0, 0))
    gd = ImageDraw.Draw(glow)
    gd.rounded_rectangle((50, 64, 462, 420), radius=34, outline=(38, 218, 255, 180), width=6)
    gd.text((126, 142), "AK", font=font_big, fill=(38, 218, 255, 180))
    glow = glow.filter(ImageFilter.GaussianBlur(8))
    image.alpha_composite(glow)

    draw.rounded_rectangle((50, 64, 462, 420), radius=34, outline=(217, 222, 226, 255), width=5)
    draw.text((126, 142), "AK", font=font_big, fill=(217, 222, 226, 255))
    draw.text((128, 142), "AK", font=font_big, fill=(94, 240, 255, 120))
    label = "abacuskit"
    bbox = draw.textbbox((0, 0), label, font=font_small)
    draw.text(((size - (bbox[2] - bbox[0])) // 2, 438), label, font=font_small, fill=(118, 205, 224, 235))
    add_scanlines(image, alpha=22)
    image.save(FIGURE_DIR / "abacuskit_icon.png")
    image.save(FIGURE_DIR / "abacuskit_icon.ico", sizes=[(256, 256), (128, 128), (64, 64), (32, 32), (16, 16)])


def main() -> int:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    make_banner()
    make_square_icon()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
