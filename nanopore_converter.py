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

try:
    import h5py
    HAS_H5PY = True
except ImportError:  # ont_fast5_api depends on h5py, so this is rare
    HAS_H5PY = False

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
        "--force",
        action="store_true",
        help="Re-convert inputs whose output pod5/ folder already exists "
             "(default: skip them). The old output is replaced only after "
             "the new conversion succeeds.",
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


# Working folders inside <output-dir>/<run>/ (removed after each conversion)
PARTIAL_DIRNAME = ".pod5_partial"
MULTI_TMP_DIRNAME = "_multi_fast5_tmp"
OLD_DIRNAME = ".pod5_old"

ARCHIVE_EXTENSIONS = (".tar.gz", ".tar.bz2", ".tgz", ".tar", ".zip")
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

    try:
        if name.endswith((".tar.gz", ".tgz", ".tar.bz2", ".tar")):
            # "r:*" auto-detects gzip/bzip2/uncompressed
            with tarfile.open(archive_path, "r:*") as tar:
                if hasattr(tarfile, "data_filter"):
                    # Reject absolute paths, ".." traversal and device files
                    tar.extractall(tmp_dir, filter="data")
                else:  # Python < 3.10.12 / 3.11.4 without extraction filters
                    tar.extractall(tmp_dir)
        elif name.endswith(".zip"):
            # zipfile sanitises absolute paths and ".." components itself
            with zipfile.ZipFile(archive_path, "r") as zf:
                zf.extractall(tmp_dir)
        else:
            raise ValueError(f"Unsupported archive format: {archive_path}")
    except BaseException:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise

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


def _classify_fast5_structure(fpath: Path) -> str | None:
    """Classify a fast5 file by its HDF5 layout; None if undeterminable.

    Multi-read files have top-level read_<uuid> groups (and usually a
    file_type=multi-read attribute); single-read files have Raw/ and
    UniqueGlobalKey/ at the root.
    """
    if not HAS_H5PY:
        return None
    try:
        with h5py.File(fpath, "r") as f:
            file_type = f.attrs.get("file_type")
            if file_type is not None:
                if isinstance(file_type, bytes):
                    file_type = file_type.decode("utf-8", errors="replace")
                file_type = str(file_type).lower()
                if file_type.startswith("multi"):
                    return "multi_read_fast5"
                if file_type.startswith("single"):
                    return "single_read_fast5"
            for key in f:
                if key.startswith("read_"):
                    return "multi_read_fast5"
            if "Raw" in f or "UniqueGlobalKey" in f:
                return "single_read_fast5"
    except Exception:
        return None
    return None


def classify_fast5_dir(dir_path: str) -> str | None:
    """Classify a directory of fast5 files as single_read or multi_read.

    Samples one .fast5 file.  When h5py is available the HDF5 layout is
    inspected (reliable for small multi-read files, e.g. fast5_fail/).
    Otherwise, or if the file cannot be parsed, the 1 MB size threshold
    is used:
    - < 1 MB: single_read_fast5
    - >= 1 MB: multi_read_fast5

    Returns None if no .fast5 files are found.
    """
    root = Path(dir_path)
    for dirpath, _, filenames in os.walk(root):
        for fname in filenames:
            if fname.lower().endswith(".fast5"):
                fpath = Path(dirpath) / fname
                by_structure = _classify_fast5_structure(fpath)
                if by_structure is not None:
                    return by_structure
                try:
                    size = fpath.stat().st_size
                except OSError:
                    continue
                if size < 1_000_000:
                    return "single_read_fast5"
                else:
                    return "multi_read_fast5"
    return None


def _has_fast5(directory: str) -> bool:
    for _, _, filenames in os.walk(directory):
        if any(f.lower().endswith(".fast5") for f in filenames):
            return True
    return False


def convert_run(
    input_dir: str, output_base: str, threads: int = 4, run_name: str | None = None
) -> bool:
    """Convert a single run directory from fast5 to pod5.

    For single_read_fast5: two-step process via intermediate multi_read_fast5.
    For multi_read_fast5: direct pod5 conversion.

    run_name sets the output folder name (default: the input directory name).
    Returns False if the input contains no .fast5 files, True otherwise.
    Raises subprocess.CalledProcessError on conversion failure.

    Output is written to a hidden <run>/.pod5_partial/ folder and moved to
    <run>/pod5/ only once every step has succeeded, so an interrupted or
    failed run never leaves a half-written pod5/ behind, and an existing
    pod5/ is only replaced by a complete new conversion.
    """
    fmt = classify_fast5_dir(input_dir)
    if fmt is None:
        return False

    if not run_name:
        run_name = Path(input_dir).resolve().name
    run_out = os.path.join(output_base, run_name)
    pod5_out = os.path.join(run_out, "pod5")
    pod5_partial = os.path.join(run_out, PARTIAL_DIRNAME)
    multi_tmp = os.path.join(run_out, MULTI_TMP_DIRNAME)
    # Leftovers from an interrupted earlier run would otherwise be mixed in
    for stale in (pod5_partial, multi_tmp):
        if os.path.isdir(stale):
            logger.info("[%s] Removing leftover %s", run_name, stale)
            shutil.rmtree(stale)
    os.makedirs(pod5_partial)

    try:
        _convert_to(fmt, input_dir, pod5_partial, multi_tmp, run_name, threads)
    except BaseException:
        shutil.rmtree(pod5_partial, ignore_errors=True)
        raise
    finally:
        shutil.rmtree(multi_tmp, ignore_errors=True)

    if os.path.isdir(pod5_out):
        # Move the old output aside (fast rename) before swapping the new one
        # in, then delete it, so pod5/ never disappears during a slow rmtree.
        logger.info("[%s] Replacing existing output %s", run_name, pod5_out)
        old_out = os.path.join(run_out, OLD_DIRNAME)
        shutil.rmtree(old_out, ignore_errors=True)
        os.replace(pod5_out, old_out)
        os.replace(pod5_partial, pod5_out)
        shutil.rmtree(old_out, ignore_errors=True)
    else:
        os.replace(pod5_partial, pod5_out)
    logger.info("[%s] Conversion complete -> %s", run_name, pod5_out)
    return True


def _convert_to(
    fmt: str, input_dir: str, pod5_out: str, multi_tmp: str, run_name: str, threads: int
) -> None:
    """Run the conversion tool(s) for one input, writing pod5 into pod5_out."""
    if fmt == "single_read_fast5":
        # Step 1: single -> multi
        os.makedirs(multi_tmp)
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
        # single_to_multi_fast5 logs per-file errors (e.g. when given
        # multi-read files) but still exits 0; make sure it produced output.
        if not _has_fast5(multi_tmp):
            logger.error(
                "single_to_multi_fast5 produced no output for %s:\n%s",
                input_dir, result.stderr,
            )
            raise subprocess.CalledProcessError(
                1, "single_to_multi_fast5", stderr="no multi-read fast5 written"
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


def get_output_name(input_path: str) -> str:
    """Derive an output folder name from the input path.

    Strips archive extensions to get a clean directory name.  The path is
    resolved first so inputs such as "." or "run/" get a real name.
    """
    name = Path(input_path).resolve().name
    # Strip known archive suffixes
    for ext in ARCHIVE_EXTENSIONS:
        if name.lower().endswith(ext):
            name = name[: -len(ext)]
            break
    return name


def is_converted(output_dir: str, run_name: str) -> bool:
    """True if <output_dir>/<run_name>/pod5/ already holds .pod5 files.

    pod5/ is only ever created by a fully successful conversion (see
    convert_run), so its presence marks a finished earlier run.
    """
    pod5_out = Path(output_dir) / run_name / "pod5"
    try:
        return any(p.suffix == ".pod5" for p in pod5_out.rglob("*"))
    except OSError:
        return False


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
            # Name the output after the archive, not the temp directory
            run_name = get_output_name(item_path)
        elif item_type == "directory":
            work_dir = item_path
            run_name = get_output_name(item_path)
        else:
            logger.error("Cannot process: %s (not a directory or archive)", item_path)
            return False

        if not convert_run(work_dir, output_dir, threads=threads, run_name=run_name):
            logger.error("No .fast5 files found in %s", item_path)
            return False
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
    skipped = 0
    failed = 0
    seen_names: dict[str, str] = {}
    for item in inputs:
        logger.info("Processing: %s", item)
        run_name = get_output_name(item)
        if run_name in seen_names:
            # Two inputs in this invocation map to the same output folder;
            # converting both would make the second overwrite (or skip) the first.
            logger.error(
                "Output name collision: %s and %s both map to %s; rename one of them.",
                seen_names[run_name], item, os.path.join(args.output_dir, run_name),
            )
            failed += 1
            continue
        seen_names[run_name] = item
        if not args.force and is_converted(args.output_dir, run_name):
            logger.info(
                "Skipping %s: already converted in %s (use --force to re-convert).",
                item, os.path.join(args.output_dir, run_name, "pod5"),
            )
            skipped += 1
            continue
        if process_item(item, args.output_dir, threads=args.threads):
            succeeded += 1
        else:
            failed += 1

    logger.info("Done. %d converted, %d skipped, %d failed.", succeeded, skipped, failed)
    if failed > 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
