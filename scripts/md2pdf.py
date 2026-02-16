#!/usr/bin/env python
"""Convert a Markdown file to PDF using pypandoc + LaTeX.

Preserves LaTeX math formatting and embeds images.

Usage:
    python scripts/md2pdf.py docs/RESULTS.md
    python scripts/md2pdf.py docs/RESULTS.md -o results/RESULTS.pdf
"""
import argparse
import os
import sys
from pathlib import Path

import pypandoc


def convert(input_path: str, output_path: str | None = None) -> None:
    input_file = Path(input_path).resolve()
    if not input_file.exists():
        print(f"Error: {input_file} not found", file=sys.stderr)
        sys.exit(1)

    if output_path is None:
        output_file = input_file.with_suffix(".pdf")
    else:
        output_file = Path(output_path).resolve()

    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Run pandoc from the markdown file's directory so relative image
    # paths (e.g. ../results/figures/foo.png) resolve correctly.
    original_dir = os.getcwd()
    os.chdir(input_file.parent)

    extra_args = [
        "--pdf-engine=xelatex",
        "-V", "geometry:margin=1in",
        "-V", "colorlinks=true",
        "-V", "linkcolor=blue",
        "-V", "urlcolor=blue",
        "-V", "mainfont=Latin Modern Roman",
        "-V", "mathfont=Latin Modern Math",
        "--standalone",
    ]

    try:
        pypandoc.convert_file(
            str(input_file),
            "pdf",
            outputfile=str(output_file),
            extra_args=extra_args,
        )
        os.chdir(original_dir)
        print(f"PDF written to: {output_file}")
    except Exception as e:
        os.chdir(original_dir)
        print(f"Error during conversion: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Markdown to PDF")
    parser.add_argument("input", help="Path to the Markdown file")
    parser.add_argument("-o", "--output", help="Output PDF path (default: same name, .pdf)")
    args = parser.parse_args()
    convert(args.input, args.output)
