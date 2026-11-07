#!/usr/bin/env python3
"""
Nanopore Format Converter
=========================
Converts single_read_fast5 and multi_read_fast5 to pod5 format.

Accepts a folder, compressed archive (.tar.gz, .tar.bz2, .zip), or a text
file listing paths (one per line, mixed folders and archives).

Dependencies (a plain `pip install ont_fast5_api` fails on current Python):
    pip install h5py pod5 progressbar2 packaging "setuptools<81"
    pip install --no-deps ont_fast5_api

Usage:
    python nanopore_converter.py /path/to/input --output-dir /path/to/output
    python nanopore_converter.py input_list.txt --output-dir /path/to/output
"""

import argparse
import logging
import os
import re
import shutil
import signal
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

if sys.version_info < (3, 12):
    sys.exit(f"{Path(__file__).name} requires Python 3.12 or newer "
             f"(running {sys.version.split()[0]})")

try:
    import h5py
    HAS_H5PY = True
except ImportError:  # ont_fast5_api depends on h5py, so this is rare
    HAS_H5PY = False

logger = logging.getLogger(__name__)

# Unreadable fast5 files skipped during this invocation (for the summary)
SKIPPED_UNREADABLE: list[Path] = []


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
        logger.error('Install with: pip install h5py pod5 progressbar2 packaging "setuptools<81" '
                     '&& pip install --no-deps ont_fast5_api')
        sys.exit(2)


# Working folders inside <output-dir>/<run>/ (removed after each conversion)
PARTIAL_DIRNAME = ".pod5_partial"
MULTI_TMP_DIRNAME = "_multi_fast5_tmp"
STAGE_DIRNAME = "_fast5_stage"
OLD_DIRNAME = ".pod5_old"
# Inside pod5/: the resolved input path the output was converted from
SOURCE_MARKER = ".source"

ARCHIVE_EXTENSIONS = (".tar.gz", ".tar.bz2", ".tar.xz", ".tgz", ".tbz2", ".txz", ".tar", ".zip")
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
                # Reject absolute paths, ".." traversal and device files
                tar.extractall(tmp_dir, filter="data")
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

    Blank lines and lines starting with # are skipped.  "~" is expanded and
    relative entries are resolved against the list file's own directory (not
    the current directory), so a list keeps working wherever it is run from.
    """
    base = Path(path).resolve().parent
    paths = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                entry = Path(line).expanduser()
                if not entry.is_absolute():
                    entry = base / entry
                paths.append(str(entry))
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


def _inspect_fast5(fpath: Path) -> tuple[str, int] | None:
    """Return ("single", 1) / ("multi", n_reads) for a readable fast5 file,
    or None if it cannot be opened or has neither layout."""
    try:
        with h5py.File(fpath, "r") as f:
            reads = sum(1 for key in f if key.startswith("read_"))
            if reads:
                return "multi", reads
            if "Raw" in f or "UniqueGlobalKey" in f:
                return "single", 1
    except Exception:
        return None
    return None


@dataclass
class Fast5Scan:
    """Every .fast5 under an input directory, grouped by kind."""
    single: list[Path] = field(default_factory=list)
    multi: list[Path] = field(default_factory=list)
    unreadable: list[Path] = field(default_factory=list)
    expected_reads: int | None = 0  # None when reads could not be counted

    @property
    def empty(self) -> bool:
        return not (self.single or self.multi)

    def first_readable(self) -> Path | None:
        return (self.single or self.multi or [None])[0]


def scan_fast5_dir(dir_path: str, exclude: list[str] | None = None) -> Fast5Scan:
    """Classify every .fast5 file under dir_path and count its reads.

    Each file is opened with h5py (about 0.2 ms per file), so runs that mix
    single- and multi-read files (e.g. single-read fast5_pass/ with a
    multi-read fast5_fail/) are converted completely, and truncated or
    corrupt files are reported instead of derailing the whole run.
    Directories in `exclude` (e.g. the output folder) are not descended into.

    Without h5py each file falls back to the 1 MB size heuristic and reads
    cannot be counted (expected_reads is None).
    """
    excluded = {os.path.realpath(e) for e in (exclude or [])}
    scan = Fast5Scan()
    if not HAS_H5PY:
        scan.expected_reads = None
    for dirpath, dirnames, filenames in os.walk(dir_path):
        dirnames[:] = sorted(
            d for d in dirnames
            if os.path.realpath(os.path.join(dirpath, d)) not in excluded
        )
        for fname in sorted(filenames):
            if not fname.lower().endswith(".fast5"):
                continue
            fpath = Path(dirpath) / fname
            if not HAS_H5PY:
                try:
                    size = fpath.stat().st_size
                except OSError:
                    scan.unreadable.append(fpath)
                    continue
                (scan.single if size < 1_000_000 else scan.multi).append(fpath)
                continue
            info = _inspect_fast5(fpath)
            if info is None:
                scan.unreadable.append(fpath)
            elif info[0] == "single":
                scan.single.append(fpath)
                scan.expected_reads += 1
            else:
                scan.multi.append(fpath)
                scan.expected_reads += info[1]
    return scan


def classify_fast5_dir(dir_path: str) -> str | None:
    """Summarise a directory's fast5 content: "single_read_fast5",
    "multi_read_fast5", "mixed_fast5", or None if no readable .fast5."""
    scan = scan_fast5_dir(dir_path)
    if scan.single and scan.multi:
        return "mixed_fast5"
    if scan.single:
        return "single_read_fast5"
    if scan.multi:
        return "multi_read_fast5"
    return None


def _stage(files: list[Path], directory: str) -> None:
    """Symlink files into a flat directory (unique names, original files
    untouched) so each conversion tool sees exactly the files meant for it."""
    os.makedirs(directory)
    for i, src in enumerate(files):
        os.symlink(os.path.abspath(src), os.path.join(directory, f"{i:07d}_{src.name}"))


_TOOL_PROBLEM = re.compile(r"error|warning|failed|not a|cannot|could not|exception", re.I)


def _run(cmd: list[str], label: str) -> str:
    """Run a conversion tool in its own process group and return its stderr.

    Lines that look like problems are logged at WARNING even when the tool
    exits 0 (both tools report per-file failures that way).  If this process
    is interrupted (Ctrl-C, or SIGTERM via the handler set in main()), the
    tool and all its worker processes are terminated before re-raising.
    """
    logger.debug("Running: %s", " ".join(cmd))
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        start_new_session=True,
    )
    try:
        stdout, stderr = proc.communicate()
    except BaseException:
        _kill_group(proc)
        raise
    problems = [ln.strip() for ln in (stderr or "").splitlines() if _TOOL_PROBLEM.search(ln)]
    for line in problems[:20]:
        logger.warning("[%s] %s", label, line)
    if len(problems) > 20:
        logger.warning("[%s] ... %d more problem lines", label, len(problems) - 20)
    if proc.returncode != 0:
        logger.error("%s failed (exit %d):\n%s", label, proc.returncode, stderr)
        raise subprocess.CalledProcessError(proc.returncode, label, stdout, stderr)
    return stderr or ""


def _kill_group(proc: subprocess.Popen) -> None:
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(proc.pid, sig)
        except (ProcessLookupError, PermissionError):
            return
        try:
            proc.wait(timeout=5)
            return
        except subprocess.TimeoutExpired:
            continue


def count_pod5_reads(directory: str) -> int:
    """Total number of reads in all .pod5 files under directory."""
    import pod5
    total = 0
    for p5 in Path(directory).rglob("*.pod5"):
        with pod5.Reader(p5) as reader:
            total += reader.num_reads
    return total


# RunInfo fields dorado uses for automatic model selection, and the key each
# one has in the chemistry dict produced by the format checker.
_RUNINFO_FIELDS = (
    ("flow_cell_product_code", "flowcell"),
    ("sequencing_kit", "kit"),
    ("sample_rate", "sample_rate"),
)


def fill_missing_pod5_metadata(pod5_dir: str, chemistry: dict) -> int:
    """Fill empty flowcell/kit/sample-rate RunInfo fields in converted pod5.

    `pod5 convert fast5` only reads tracking_id, so pod5 files converted
    from older fast5 (flowcell stored only in context_tags/flowcell_type,
    or kit only as experiment_kit) lack the fields dorado needs for
    automatic model selection.  Values already present are never changed.
    Best-effort: a file that cannot be patched is left as written and a
    warning is logged.  Returns the number of files rewritten.
    """
    try:
        import dataclasses
        import pod5
    except ImportError:
        logger.warning("pod5 Python package not available; skipping metadata fill")
        return 0
    patched_files = 0
    for p5 in sorted(Path(pod5_dir).rglob("*.pod5")):
        tmp = p5.with_name(p5.name + ".tmp")
        try:
            with pod5.Reader(p5) as reader:
                needs = any(
                    not getattr(rec.run_info, field) and chemistry.get(key)
                    for rec in reader.reads()
                    for field, key in _RUNINFO_FIELDS
                )
            if not needs:
                continue
            patched = {}
            with pod5.Reader(p5) as reader, pod5.Writer(tmp) as writer:
                for rec in reader.reads():
                    read = rec.to_read()
                    acq = read.run_info.acquisition_id
                    if acq not in patched:
                        updates = {
                            field: chemistry[key]
                            for field, key in _RUNINFO_FIELDS
                            if not getattr(read.run_info, field) and chemistry.get(key)
                        }
                        patched[acq] = dataclasses.replace(read.run_info, **updates)
                    read.run_info = patched[acq]
                    writer.add_read(read)
            os.replace(tmp, p5)
            patched_files += 1
        except Exception as e:
            logger.warning("Could not fill metadata in %s: %s", p5, e)
            tmp.unlink(missing_ok=True)
    return patched_files


def _detect_chemistry(sample: Path | None) -> dict | None:
    """Read flowcell/kit/sample rate from one readable source fast5
    (best-effort).

    Reuses the format checker's multi-fallback extraction when it is
    importable (it lives next to this script).
    """
    if sample is None:
        return None
    try:
        from nanopore_format_checker import extract_chemistry_fast5
    except ImportError:
        logger.debug("nanopore_format_checker not importable; skipping metadata fill")
        return None
    try:
        return extract_chemistry_fast5(sample)
    except Exception:
        return None


class ConversionError(RuntimeError):
    """A conversion finished but its output cannot be trusted."""


def convert_run(
    input_dir: str, output_base: str, threads: int = 4, run_name: str | None = None,
    source: str | None = None,
) -> bool:
    """Convert every fast5 file under input_dir to pod5.

    Single-read files go through single_to_multi_fast5 first; multi-read
    files go to pod5 convert directly; a run may contain both.  Unreadable
    files are skipped with a warning.  The pod5 output must contain exactly
    as many reads as the readable input files, otherwise ConversionError is
    raised and nothing is published.

    run_name sets the output folder name (default: the input directory
    name); source is recorded in pod5/.source (default: input_dir).
    Returns False if the input contains no readable .fast5 files.
    Raises subprocess.CalledProcessError or ConversionError on failure.

    Output is written to a hidden <run>/.pod5_partial/ folder and moved to
    <run>/pod5/ only once every step has succeeded, so an interrupted or
    failed run never leaves a half-written pod5/ behind, and an existing
    pod5/ is only replaced by a complete new conversion.
    """
    if not run_name:
        run_name = Path(input_dir).resolve().name
    run_out = os.path.join(output_base, run_name)
    scan = scan_fast5_dir(input_dir, exclude=[output_base])
    if scan.unreadable:
        SKIPPED_UNREADABLE.extend(scan.unreadable)
        shown = ", ".join(str(p) for p in scan.unreadable[:5])
        more = f" (+{len(scan.unreadable) - 5} more)" if len(scan.unreadable) > 5 else ""
        logger.warning("[%s] Skipping %d unreadable/corrupt fast5 file(s): %s%s",
                       run_name, len(scan.unreadable), shown, more)
    if scan.empty:
        return False
    logger.info("[%s] Found %d single-read and %d multi-read fast5 file(s)%s",
                run_name, len(scan.single), len(scan.multi),
                f", {scan.expected_reads} reads" if scan.expected_reads is not None else "")

    pod5_out = os.path.join(run_out, "pod5")
    pod5_partial = os.path.join(run_out, PARTIAL_DIRNAME)
    multi_tmp = os.path.join(run_out, MULTI_TMP_DIRNAME)
    stage = os.path.join(run_out, STAGE_DIRNAME)
    # Leftovers from an interrupted earlier run would otherwise be mixed in
    for stale in (pod5_partial, multi_tmp, stage):
        if os.path.lexists(stale):
            logger.info("[%s] Removing leftover %s", run_name, stale)
            shutil.rmtree(stale)
    os.makedirs(pod5_partial)

    try:
        pod5_inputs = []
        if scan.single:
            _stage(scan.single, os.path.join(stage, "single"))
            os.makedirs(multi_tmp)
            logger.info("[%s] Merging %d single-read files (single_to_multi_fast5)",
                        run_name, len(scan.single))
            _run(["single_to_multi_fast5", "-i", os.path.join(stage, "single"),
                  "-s", multi_tmp, "-t", str(threads), "--recursive"],
                 "single_to_multi_fast5")
            pod5_inputs.append(multi_tmp + "/")
        if scan.multi:
            _stage(scan.multi, os.path.join(stage, "multi"))
            pod5_inputs.append(os.path.join(stage, "multi") + "/")
        logger.info("[%s] Converting to pod5 (pod5 convert fast5)", run_name)
        _run(["pod5", "convert", "fast5", *pod5_inputs, "--output", pod5_partial + "/",
              "--threads", str(threads), "--recursive"], "pod5 convert fast5")

        _verify_read_count(pod5_partial, scan, run_name)

        chemistry = _detect_chemistry(scan.first_readable())
        if chemistry:
            n = fill_missing_pod5_metadata(pod5_partial, chemistry)
            if n:
                logger.info(
                    "[%s] Filled missing pod5 metadata in %d file(s) "
                    "(flowcell=%s, kit=%s, sample_rate=%s)", run_name, n,
                    chemistry.get("flowcell"), chemistry.get("kit"),
                    chemistry.get("sample_rate"),
                )
        Path(pod5_partial, SOURCE_MARKER).write_text(
            str(Path(source or input_dir).resolve()) + "\n")
    except BaseException:
        shutil.rmtree(pod5_partial, ignore_errors=True)
        raise
    finally:
        shutil.rmtree(multi_tmp, ignore_errors=True)
        shutil.rmtree(stage, ignore_errors=True)

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


def _verify_read_count(pod5_dir: str, scan: Fast5Scan, run_name: str) -> None:
    """Fail unless the pod5 output holds exactly the reads found in the input.

    Both conversion tools exit 0 after skipping files they cannot handle,
    so this is the only reliable signal that no reads were lost.
    """
    if scan.expected_reads is None:
        logger.warning("[%s] h5py unavailable: cannot verify the read count", run_name)
        return
    got = count_pod5_reads(pod5_dir)
    if got != scan.expected_reads:
        raise ConversionError(
            f"pod5 output has {got} reads but the input fast5 files contain "
            f"{scan.expected_reads}; see the tool warnings above"
        )
    logger.info("[%s] Verified %d reads in pod5 output", run_name, got)


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


def converted_source(output_dir: str, run_name: str) -> str | None:
    """The input path recorded in pod5/.source, or None (absent/older output)."""
    try:
        return (Path(output_dir) / run_name / "pod5" / SOURCE_MARKER).read_text().strip()
    except OSError:
        return None


def _overlap_problem(item: str, run_out: str) -> str | None:
    """Explain why converting directory `item` into `run_out` is unsafe.

    An output folder inside the input is fine (the scan skips it); the
    input itself, or a folder containing it, is not: pod5/ there may be the
    run's own original data, which --force would replace.
    """
    src = Path(item).resolve()
    out = Path(run_out).resolve()
    if out == src or out in src.parents:
        return (f"output folder {out} is the input itself or contains it; "
                f"converting would replace data inside the input")
    return None


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

        if not convert_run(work_dir, output_dir, threads=threads, run_name=run_name,
                           source=item_path):
            logger.error("No readable .fast5 files found in %s", item_path)
            return False
        return True

    except (subprocess.CalledProcessError, ConversionError) as e:
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

    # Turn SIGTERM (kill, batch-scheduler time limits) into an exception so
    # temp folders are cleaned up and child tools are terminated, as on Ctrl-C.
    def _on_sigterm(signum, frame):
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, _on_sigterm)

    inputs = resolve_inputs(args.input)
    if not inputs:
        logger.error("No valid inputs found.")
        sys.exit(2)

    os.makedirs(args.output_dir, exist_ok=True)

    succeeded = 0
    skipped = 0
    failed = 0
    seen_names: dict[str, str] = {}
    seen_paths: set[Path] = set()
    for item in inputs:
        resolved = Path(item).resolve()
        if resolved in seen_paths:
            logger.warning("Skipping duplicate entry: %s", item)
            continue
        seen_paths.add(resolved)
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
        run_out = os.path.join(args.output_dir, run_name)
        if os.path.isdir(item):
            problem = _overlap_problem(item, run_out)
            if problem:
                logger.error("Refusing to convert %s: %s. Use a different --output-dir.",
                             item, problem)
                failed += 1
                continue
        previous = converted_source(args.output_dir, run_name)
        if previous is not None and previous != str(resolved):
            logger.error(
                "Refusing to convert %s: %s holds output converted from a "
                "different input (%s). Use a different --output-dir or rename.",
                item, os.path.join(run_out, "pod5"), previous,
            )
            failed += 1
            continue
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
    if SKIPPED_UNREADABLE:
        logger.warning(
            "%d unreadable/corrupt fast5 file(s) were not converted (listed above); "
            "their reads are missing from the output.", len(SKIPPED_UNREADABLE))
    if failed > 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
