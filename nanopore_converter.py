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


def check_tools():
    """Verify that required CLI tools are on PATH."""
    missing = []
    for tool in ("single_to_multi_fast5", "pod5"):
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        logger.error("Required tools not found on PATH: %s", ", ".join(missing))
        logger.error("Install with: pip install ont_fast5_api pod5")
        sys.exit(2)


ARCHIVE_EXTENSIONS = {".tar.gz", ".tar.bz2", ".tgz", ".zip"}
LIST_EXTENSIONS = {".txt", ".csv", ".list"}


def detect_input_type(path_str: str) -> str | None:
    """Determine if input is a directory, archive, or list file.

    Returns 'directory', 'archive', 'listfile', or None.
    """
    p = Path(path_str)
    if not p.exists():
        return None
    if p.is_dir():
        return "directory"
    if p.is_file():
        # Check for double extensions like .tar.gz
        name = p.name.lower()
        for ext in ARCHIVE_EXTENSIONS:
            if name.endswith(ext):
                return "archive"
        if p.suffix.lower() in LIST_EXTENSIONS:
            return "listfile"
    return None


def extract_archive(archive_path: str) -> str:
    """Extract a compressed archive to a temporary directory.

    Returns the path to the temporary directory. Caller is responsible for
    cleanup.
    """
    tmp_dir = tempfile.mkdtemp(prefix="nanopore_conv_")
    p = Path(archive_path)
    name = p.name.lower()

    logger.info("Extracting %s to %s", archive_path, tmp_dir)

    if name.endswith((".tar.gz", ".tgz", ".tar.bz2")):
        mode = "r:gz" if name.endswith((".tar.gz", ".tgz")) else "r:bz2"
        with tarfile.open(archive_path, mode) as tar:
            tar.extractall(tmp_dir)
    elif name.endswith(".zip"):
        with zipfile.ZipFile(archive_path, "r") as zf:
            zf.extractall(tmp_dir)
    else:
        shutil.rmtree(tmp_dir)
        raise ValueError(f"Unsupported archive format: {archive_path}")

    return tmp_dir


def parse_list_file(path: str) -> list[str]:
    """Read a text file of paths, one per line.

    Blank lines and lines starting with # are skipped.
    """
    paths = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                paths.append(line)
    return paths


def resolve_inputs(input_path: str) -> list[str]:
    """Resolve the input argument to a list of paths to process.

    Handles directories, archives (returned as-is), and list files.
    """
    input_type = detect_input_type(input_path)
    if input_type == "directory" or input_type == "archive":
        return [input_path]
    elif input_type == "listfile":
        return parse_list_file(input_path)
    else:
        logger.error("Input not found or unrecognized: %s", input_path)
        return []


def classify_fast5_dir(dir_path: str) -> str | None:
    """Classify a directory of fast5 files as single_read or multi_read.

    Samples one .fast5 file and uses the 1 MB size threshold:
    - < 1 MB: single_read_fast5
    - >= 1 MB: multi_read_fast5

    Returns None if no .fast5 files are found.
    """
    root = Path(dir_path)
    for dirpath, _, filenames in os.walk(root):
        for fname in filenames:
            if fname.lower().endswith(".fast5"):
                fpath = Path(dirpath) / fname
                try:
                    size = fpath.stat().st_size
                except OSError:
                    continue
                if size < 1_000_000:
                    return "single_read_fast5"
                else:
                    return "multi_read_fast5"
    return None


def convert_run(input_dir: str, output_base: str, threads: int = 4) -> None:
    """Convert a single run directory from fast5 to pod5.

    For single_read_fast5: two-step process via intermediate multi_read_fast5.
    For multi_read_fast5: direct pod5 conversion.

    Raises subprocess.CalledProcessError on conversion failure.
    """
    fmt = classify_fast5_dir(input_dir)
    if fmt is None:
        logger.warning("No .fast5 files found in %s, skipping.", input_dir)
        return

    run_name = Path(input_dir).name
    pod5_out = os.path.join(output_base, run_name, "pod5")
    os.makedirs(pod5_out, exist_ok=True)

    if fmt == "single_read_fast5":
        # Step 1: single -> multi
        multi_tmp = os.path.join(output_base, run_name, "_multi_fast5_tmp")
        os.makedirs(multi_tmp, exist_ok=True)
        logger.info("[%s] Step 1/2: single_to_multi_fast5", run_name)
        result = subprocess.run(
            [
                "single_to_multi_fast5",
                "-i", input_dir,
                "-s", multi_tmp,
                "-t", str(threads),
                "--recursive",
            ],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            logger.error("single_to_multi_fast5 failed:\n%s", result.stderr)
            raise subprocess.CalledProcessError(
                result.returncode, "single_to_multi_fast5"
            )

        # Step 2: multi -> pod5
        logger.info("[%s] Step 2/2: pod5 convert fast5", run_name)
        result = subprocess.run(
            [
                "pod5", "convert", "fast5",
                multi_tmp + "/",
                "--output", pod5_out + "/",
                "--threads", str(threads),
                "--recursive",
            ],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            logger.error("pod5 convert failed:\n%s", result.stderr)
            raise subprocess.CalledProcessError(
                result.returncode, "pod5 convert fast5"
            )

        # Cleanup intermediate files
        logger.info("[%s] Cleaning up intermediate multi_fast5.", run_name)
        shutil.rmtree(multi_tmp, ignore_errors=True)

    elif fmt == "multi_read_fast5":
        logger.info("[%s] Converting multi_read_fast5 -> pod5", run_name)
        result = subprocess.run(
            [
                "pod5", "convert", "fast5",
                input_dir + "/",
                "--output", pod5_out + "/",
                "--threads", str(threads),
                "--recursive",
            ],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            logger.error("pod5 convert failed:\n%s", result.stderr)
            raise subprocess.CalledProcessError(
                result.returncode, "pod5 convert fast5"
            )

    logger.info("[%s] Conversion complete -> %s", run_name, pod5_out)


def get_output_name(input_path: str) -> str:
    """Derive an output folder name from the input path.

    Strips archive extensions to get a clean directory name.
    """
    name = Path(input_path).name
    # Strip known archive suffixes
    for ext in (".tar.gz", ".tar.bz2", ".tgz", ".zip"):
        if name.lower().endswith(ext):
            name = name[: -len(ext)]
            break
    return name


def process_item(item_path: str, output_dir: str, threads: int = 4) -> bool:
    """Process a single input item (directory or archive).

    Returns True on success, False on failure.
    """
    item_type = detect_input_type(item_path)
    tmp_dir = None

    try:
        if item_type == "archive":
            tmp_dir = extract_archive(item_path)
            # The extracted archive may contain a single top-level folder
            # or files directly. Find the directory with fast5 files.
            work_dir = tmp_dir
            # If there is exactly one subdirectory, use it
            entries = list(Path(tmp_dir).iterdir())
            if len(entries) == 1 and entries[0].is_dir():
                work_dir = str(entries[0])
        elif item_type == "directory":
            work_dir = item_path
        else:
            logger.error("Cannot process: %s (not a directory or archive)", item_path)
            return False

        convert_run(work_dir, output_dir, threads=threads)
        return True

    except subprocess.CalledProcessError as e:
        logger.error("Conversion failed for %s: %s", item_path, e)
        return False
    except Exception as e:
        logger.error("Unexpected error processing %s: %s", item_path, e)
        return False
    finally:
        if tmp_dir and os.path.isdir(tmp_dir):
            logger.info("Cleaning up temp directory: %s", tmp_dir)
            shutil.rmtree(tmp_dir, ignore_errors=True)


def main(argv=None):
    args = parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    check_tools()

    inputs = resolve_inputs(args.input)
    if not inputs:
        logger.error("No valid inputs found.")
        sys.exit(2)

    os.makedirs(args.output_dir, exist_ok=True)

    succeeded = 0
    failed = 0
    for item in inputs:
        logger.info("Processing: %s", item)
        if process_item(item, args.output_dir, threads=args.threads):
            succeeded += 1
        else:
            failed += 1

    logger.info("Done. %d succeeded, %d failed.", succeeded, failed)
    if failed > 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
