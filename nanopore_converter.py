#!/usr/bin/env python3
"""
Nanopore Format Converter
=========================
Converts single_read_fast5 and multi_read_fast5 to pod5 format.

Accepts a folder, compressed archive (.tar.gz, .tar.bz2, .zip), or a text
file listing paths (one per line, mixed folders and archives).

Dependencies:
    pip install ont_fast5_api pod5

Usage:
    python nanopore_converter.py /path/to/input --output-dir /path/to/output
    python nanopore_converter.py input_list.txt --output-dir /path/to/output
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path

logger = logging.getLogger(__name__)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert nanopore fast5 files to pod5 format."
    )
    parser.add_argument(
        "input",
        help="Folder, compressed archive, or text file listing paths (one per line).",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Base output directory for converted pod5 files.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Thread count for conversion tools (default: 4).",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable detailed logging.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    logger.info("Nanopore converter starting.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
