#!/usr/bin/env python3
"""
Nanopore Run Format Checker
============================
Scans a target folder containing nanopore run directories and classifies each
run by its read format: multi_read_fast5, single_read_fast5, or pod5.

Usage:
    python nanopore_format_checker.py /path/to/runs_folder [--convert-to pod5|single_fast5]

Dependencies:
    - h5py (for inspecting fast5 files to distinguish single vs multi)
    - pod5 (optional, for pod5 validation)

Install:
    pip install h5py ont_fast5_api pod5
"""

import argparse
import os
import re
import sys
from pathlib import Path

if sys.version_info < (3, 12):
    sys.exit(f"{Path(__file__).name} requires Python 3.12 or newer "
             f"(running {sys.version.split()[0]})")

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False

try:
    import pod5
    HAS_POD5 = True
except ImportError:
    HAS_POD5 = False

FLOWCELL_PORE = {
    # R10.4.1
    "FLO-MIN114": "R10.4.1",
    "FLO-PRO114": "R10.4.1",
    "FLO-FLG114": "R10.4.1",
    "FLO-PRO114M": "R10.4.1",
    # R10.4.1 HD (high-duplex flowcells)
    "FLO-MIN114HD": "R10.4.1",
    "FLO-PRO114HD": "R10.4.1",
    "FLO-FLG114HD": "R10.4.1",
    # R10.4 (intermediate, pre-R10.4.1)
    "FLO-MIN112": "R10.4",
    "FLO-PRO112": "R10.4",
    # R10.3 (short-lived, early 2020)
    "FLO-MIN111": "R10.3",
    # R9.4.1
    "FLO-MIN106": "R9.4.1",
    "FLO-MIN106D": "R9.4.1",
    "FLO-FLG001": "R9.4.1",
    "FLO-PRO002": "R9.4.1",
    "FLO-PRO002M": "R9.4.1",
    # R9.4.1 special
    "FLO-MINSP6": "R9.4.1",
    "FLO-PRO001": "R9.4.1",
    "FLO-PRO002-ECO": "R9.4.1",
    # R9.5.1 (1D^2 flowcell; dorado only supports it for RNA002, not DNA)
    "FLO-MIN107": "R9.5.1",
    # RNA (direct RNA flowcells)
    "FLO-MIN004RA": "RNA004",
    "FLO-PRO004RA": "RNA004",
}

# Fallback: infer pore type from sequencing kit when flowcell product code is
# absent (common with older MinKNOW versions that did not write
# flow_cell_product_code).
KIT_PORE = {
    # R9.4.1 ligation / rapid / barcoding kits
    "SQK-LSK109": "R9.4.1",
    "SQK-LSK109-XL": "R9.4.1",
    "SQK-RAD004": "R9.4.1",
    "SQK-RBK004": "R9.4.1",
    "SQK-RPB004": "R9.4.1",
    "SQK-PBK004": "R9.4.1",
    "SQK-PCS109": "R9.4.1",
    "SQK-PCB109": "R9.4.1",
    "SQK-DCS109": "R9.4.1",
    "SQK-CAS109": "R9.4.1",
    "SQK-CS9109": "R9.4.1",
    # Kit 10/11 (R9.4.1 e8 in dorado's chemistry table)
    "SQK-LSK110": "R9.4.1",
    "SQK-LSK110-XL": "R9.4.1",
    "SQK-LSK111": "R9.4.1",
    "SQK-LSK111-XL": "R9.4.1",
    "SQK-NBD111-24": "R9.4.1",
    "SQK-NBD111-96": "R9.4.1",
    "SQK-RBK110-96": "R9.4.1",
    "SQK-RBK111-24": "R9.4.1",
    "SQK-RBK111-96": "R9.4.1",
    "SQK-PCB110": "R9.4.1",
    "SQK-PCB111-24": "R9.4.1",
    "SQK-PCS111": "R9.4.1",
    "SQK-MLK111-96-XL": "R9.4.1",
    "SQK-ULK001": "R9.4.1",
    "SQK-16S024": "R9.4.1",
    # Older R9.4.1 kits
    "SQK-RAD002": "R9.4.1",
    "SQK-RAD003": "R9.4.1",
    "SQK-LSK108": "R9.4.1",
    "SQK-LWP001": "R9.4.1",
    "SQK-LWB001": "R9.4.1",
    "SQK-RAS201": "R9.4.1",
    "SQK-RLB001": "R9.4.1",
    "VSK-VSK002": "R9.4.1",
    # Kit 12 (Q20+ chemistry for R10.4 / FLO-MIN112; dorado has no model)
    "SQK-LSK112": "R10.4",
    "SQK-LSK112-XL": "R10.4",
    "SQK-NBD112-24": "R10.4",
    "SQK-NBD112-96": "R10.4",
    # R10.4.1 kits (Kit 14)
    "SQK-LSK114": "R10.4.1",
    "SQK-LSK114-XL": "R10.4.1",
    "SQK-RAD114": "R10.4.1",
    "SQK-RBK114-24": "R10.4.1",
    "SQK-RBK114-96": "R10.4.1",
    "SQK-NBD114-24": "R10.4.1",
    "SQK-NBD114-96": "R10.4.1",
    "SQK-PCS114": "R10.4.1",
    "SQK-PCB114-24": "R10.4.1",
    "SQK-ULK114": "R10.4.1",
    "SQK-RPB114-24": "R10.4.1",
    "SQK-MLK114-96-XL": "R10.4.1",
    "SQK-16S114-24": "R10.4.1",
    "SQK-MAB114-24": "R10.4.1",
    "SQK-APK114": "R10.4.1",
    # RNA kits
    "SQK-RNA002": "RNA002",
    "SQK-RNA004": "RNA004",
}


# Gzipped files that are clearly not raw-signal archives (MinKNOW summaries,
# basecalled reads, reports).  These must not be mistaken for compressed
# single-read fast5 archives.
NON_RAW_GZ_SUFFIXES = (
    ".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz", ".bam.gz", ".sam.gz",
    ".txt.gz", ".tsv.gz", ".csv.gz", ".json.gz", ".html.gz", ".md.gz",
    ".log.gz", ".pdf.gz",
)


# Archive suffixes the converter can extract (longest first, for stripping)
RAW_ARCHIVE_SUFFIXES = (".tar.gz", ".tar.bz2", ".tar.xz", ".tgz", ".tbz2", ".txz",
                        ".tar", ".zip")


def is_raw_archive(filename: str) -> bool:
    """Return True if a filename looks like a compressed raw-data archive.

    Accepts tar archives (plain, gzip, bzip2, xz), .zip and plain .gz files,
    but rejects gzipped text/sequence outputs such as sequencing_summary.txt.gz
    or reads.fastq.gz.
    """
    name = filename.lower()
    if name.endswith(RAW_ARCHIVE_SUFFIXES):
        return True
    if name.endswith(".gz"):
        return not name.endswith(NON_RAW_GZ_SUFFIXES)
    return False


def _sh_quote(value: str) -> str:
    """Always single-quote a string for bash, escaping embedded single quotes."""
    return "'" + str(value).replace("'", "'\"'\"'") + "'"


def is_nanopore_run_dir(dirname: str) -> bool:
    """Check if directory name matches nanopore run naming convention (starts with date)."""
    return bool(re.match(r"^\d{8}_", dirname))


def find_named_subdirs(
    run_path: Path,
    name: str | None = None,
    prefix: str | None = None,
    max_depth: int = 2,
) -> list[Path]:
    """
    Find subdirectories by exact name or prefix using os.scandir.

    Only descends into directories, skipping all files at every level.
    This avoids enumerating large numbers of files (100k+ in single-read
    fast5 directories) that glob-based approaches would touch.
    """
    matches = []

    def _walk(directory: str, depth: int) -> None:
        if depth > max_depth:
            return
        try:
            with os.scandir(directory) as it:
                for entry in it:
                    if not entry.is_dir(follow_symlinks=False):
                        continue
                    matched = False
                    if name is not None and entry.name == name:
                        matches.append(Path(entry.path))
                        matched = True
                    elif prefix is not None and entry.name.startswith(prefix):
                        matches.append(Path(entry.path))
                        matched = True
                    if not matched:
                        _walk(entry.path, depth + 1)
        except PermissionError:
            pass

    _walk(str(run_path), 1)
    return matches


def discover_run_structure(run_path: Path, max_depth: int = 5) -> dict[str, list[Path]]:
    """
    Walk the run directory tree once and categorize all format-related subdirectories.

    Returns a dict with keys:
        "pod5"         - directories named exactly "pod5"
        "pod5_prefix"  - directories starting with "pod5_" (pod5_pass, pod5_fail, etc.)
        "fast5"        - directories named exactly "fast5"
        "fast5_prefix" - directories starting with "fast5_" (fast5_pass, fast5_fail, etc.)
        "fastq_prefix" - directories starting with "fastq" (fastq_pass, fastq_fail, etc.)
    """
    result = {
        "pod5": [],
        "pod5_prefix": [],
        "fast5": [],
        "fast5_prefix": [],
        "fastq_prefix": [],
    }

    # Symlinked folders are followed (data is often linked in from other
    # storage); real paths already visited are skipped to avoid loops.
    visited = {os.path.realpath(run_path)}

    def _walk(directory: str, depth: int) -> None:
        if depth > max_depth:
            return
        try:
            with os.scandir(directory) as it:
                for entry in it:
                    try:
                        if not entry.is_dir():
                            continue
                    except OSError:
                        continue
                    name = entry.name
                    if name.startswith("."):
                        # Hidden work folders (.convert_tmp, .pod5_partial) hold
                        # incomplete conversions, not run data
                        continue
                    real = os.path.realpath(entry.path)
                    if real in visited:
                        continue
                    visited.add(real)
                    path = Path(entry.path)
                    matched = False
                    if name == "pod5":
                        result["pod5"].append(path)
                        matched = True
                    elif name.startswith("pod5_"):
                        result["pod5_prefix"].append(path)
                        matched = True
                    elif name == "fast5":
                        result["fast5"].append(path)
                        matched = True
                    elif name.startswith("fast5_"):
                        result["fast5_prefix"].append(path)
                        matched = True
                    elif name.startswith("fastq"):
                        result["fastq_prefix"].append(path)
                        matched = True

                    if not matched:
                        _walk(entry.path, depth + 1)
        except PermissionError:
            pass

    _walk(str(run_path), 1)
    return result


def _is_dir_readable(directory: Path) -> bool:
    """Check if a directory is readable without enumerating its contents."""
    try:
        with os.scandir(directory) as it:
            next(it, None)
        return True
    except PermissionError:
        return False


def find_files_with_ext(directory: Path, ext: str, limit: int = 5) -> list[Path]:
    """Find files with given extension in directory (non-recursive), up to limit."""
    found = []
    try:
        for f in directory.iterdir():
            if f.is_file() and f.suffix == ext:
                found.append(f)
                if len(found) >= limit:
                    break
    except PermissionError:
        pass
    return found


def fast_count_files(
    directory: Path,
    ext: str | None = None,
    recursive: bool = False,
    extensions: tuple[str, ...] | None = None,
) -> tuple[int, int]:
    """
    Count files matching extension(s) using os.scandir for speed.

    Accepts either a single ext (".fast5") or a tuple of extensions
    (".fastq", ".fastq.gz", ".fq", ".fq.gz") for single-pass counting.

    Returns (file_count, total_size_bytes). The size is accumulated from
    stat info that is already cached by the OS after the is_file() check.
    """
    count = 0
    size = 0
    try:
        with os.scandir(directory) as it:
            for entry in it:
                if entry.is_file(follow_symlinks=False):
                    matched = False
                    if ext is not None and entry.name.endswith(ext):
                        matched = True
                    elif extensions is not None and any(
                        entry.name.endswith(e) for e in extensions
                    ):
                        matched = True
                    if matched:
                        count += 1
                        try:
                            size += entry.stat(follow_symlinks=False).st_size
                        except OSError:
                            pass
                elif recursive and entry.is_dir(follow_symlinks=False):
                    sub_count, sub_size = fast_count_files(
                        Path(entry.path), ext=ext, recursive=True, extensions=extensions
                    )
                    count += sub_count
                    size += sub_size
    except PermissionError:
        pass
    return count, size


def estimate_dir_size(
    directory: Path,
    ext: str,
    sample_size: int = 50,
    recursive: bool = False,
) -> tuple[int, int, bool]:
    """
    Count matching files and estimate total size via sampling.

    Iterates all directory entries (cheap readdir) but limits stat() calls
    to sample_size files. Remaining file sizes are extrapolated from the
    sampled average. For directories with fewer matching files than
    sample_size, the returned size is exact.

    Returns (file_count, estimated_size_bytes, is_estimated).
    """
    count = 0
    sampled_sizes = []

    def _scan(path: str) -> None:
        nonlocal count
        try:
            with os.scandir(path) as it:
                for entry in it:
                    if entry.is_file(follow_symlinks=False) and entry.name.endswith(ext):
                        count += 1
                        if len(sampled_sizes) < sample_size:
                            try:
                                sampled_sizes.append(entry.stat(follow_symlinks=False).st_size)
                            except OSError:
                                pass
                    elif recursive and entry.is_dir(follow_symlinks=False):
                        _scan(entry.path)
        except PermissionError:
            pass

    _scan(str(directory))

    if not sampled_sizes:
        return count, 0, False
    if count <= len(sampled_sizes):
        return count, sum(sampled_sizes), False

    avg_size = sum(sampled_sizes) / len(sampled_sizes)
    exact_portion = sum(sampled_sizes)
    estimated_portion = int(avg_size * (count - len(sampled_sizes)))
    return count, exact_portion + estimated_portion, True


def classify_fast5_by_size(fast5_file: Path) -> str:
    """
    Determine if a fast5 file is single-read or multi-read based on file size.
    
    Single-read fast5: typically 1-50 KB (one read per file)
    Multi-read fast5: typically 1+ MB (hundreds/thousands of reads per file)
    
    This is much faster than opening the HDF5 file, especially for large directories.
    """
    try:
        size_bytes = fast5_file.stat().st_size
        # Threshold: files under 1 MB are likely single-read
        if size_bytes < 1_000_000:  # 1 MB
            return "single_read_fast5"
        else:
            return "multi_read_fast5"
    except (OSError, PermissionError):
        return "fast5_unknown"


def classify_fast5_by_structure(fast5_file: Path) -> str | None:
    """
    Determine single vs multi-read by inspecting the HDF5 layout (needs h5py).

    Multi-read files have top-level ``read_<uuid>`` groups (and usually a
    ``file_type = multi-read`` root attribute); single-read files have
    ``Raw/`` and ``UniqueGlobalKey/`` at the top level.  Only the root group
    is read, so this costs a single file open.

    Returns None when h5py is unavailable or the file cannot be parsed.
    """
    if not HAS_H5PY:
        return None
    try:
        with h5py.File(fast5_file, "r") as f:
            file_type = f.attrs.get("file_type")
            if file_type is not None:
                file_type = _decode_attr(file_type).lower()
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


def find_fast5_samples(directories: list[Path], limit: int = 5, max_depth: int = 4,
                       max_entries: int = 50_000) -> list[Path]:
    """Find up to `limit` .fast5 files under the given directories.

    Files in a directory are taken before descending, so a 100k-file
    single-read folder costs only a few entries.  Sub-folders are searched
    up to max_depth levels (barcode01/0/, pass/0/, ...), symlinked folders
    are followed with a loop guard, and at most max_entries directory
    entries are examined in total.
    """
    samples: list[Path] = []
    visited: set[str] = set()
    budget = max_entries

    def _search(path: str, depth: int) -> None:
        nonlocal budget
        real = os.path.realpath(path)
        if depth > max_depth or real in visited or budget <= 0:
            return
        visited.add(real)
        subdirs = []
        try:
            with os.scandir(path) as it:
                for entry in it:
                    budget -= 1
                    if budget <= 0 or len(samples) >= limit:
                        break
                    try:
                        if entry.is_file() and entry.name.lower().endswith(".fast5"):
                            samples.append(Path(entry.path))
                        elif entry.is_dir():
                            subdirs.append(entry.path)
                    except OSError:
                        continue
        except OSError:
            return
        for sub in sorted(subdirs):
            if len(samples) >= limit:
                return
            _search(sub, depth + 1)

    for d in directories:
        if len(samples) >= limit:
            break
        _search(str(d), 1)
    return samples


def classify_fast5_samples(samples: list[Path]) -> tuple[str, Path | None]:
    """Classify a run from several sample files.

    Uses the first sample whose HDF5 layout can be read, so one truncated or
    corrupt file (common after a MinKNOW crash) does not decide the whole
    run.  Only if no sample can be parsed (or h5py is missing) does the size
    heuristic apply.  Returns (classification, sample used).
    """
    for sample in samples:
        by_structure = classify_fast5_by_structure(sample)
        if by_structure:
            return by_structure, sample
    for sample in samples:
        by_size = classify_fast5_by_size(sample)
        if by_size != "fast5_unknown":
            return by_size, sample
    return "fast5_unknown", (samples[0] if samples else None)


def classify_fast5(fast5_file: Path) -> str:
    """
    Determine if a fast5 file is single-read or multi-read.

    Inspects the HDF5 root group when h5py is available (one file open),
    which correctly handles small multi-read files such as those in
    fast5_fail/ or short runs.  Falls back to the file size heuristic
    (single-read <1MB, multi-read >=1MB) when h5py is missing or the file
    cannot be parsed.
    """
    return classify_fast5_by_structure(fast5_file) or classify_fast5_by_size(fast5_file)


def _decode_attr(val):
    """Decode an HDF5 attribute value to a string."""
    if isinstance(val, bytes):
        return val.decode("utf-8", errors="replace")
    return str(val)


def extract_chemistry_fast5(file_path: Path) -> dict | None:
    """Extract flowcell chemistry metadata from a fast5 file.

    Reads context_tags and tracking_id attributes from the HDF5 file.
    Handles both single-read (UniqueGlobalKey/) and multi-read
    (<read_id>/) layouts. Falls back to tracking_id/flow_cell_product_code
    when context_tags/flowcell_type is empty or missing.

    Returns {"flowcell": "FLO-MIN114", "kit": "SQK-LSK114",
    "sample_rate": 5000} or None on failure.
    """
    if not HAS_H5PY:
        return None
    try:
        with h5py.File(file_path, "r") as f:
            ctx = None
            trk = None

            if "UniqueGlobalKey" in f:
                if "context_tags" in f["UniqueGlobalKey"]:
                    ctx = f["UniqueGlobalKey/context_tags"]
                if "tracking_id" in f["UniqueGlobalKey"]:
                    trk = f["UniqueGlobalKey/tracking_id"]
            else:
                # Multi-read layout: check first few reads for metadata groups
                for key in list(f.keys())[:5]:
                    try:
                        grp = f[key]
                        if ctx is None and "context_tags" in grp:
                            ctx = grp["context_tags"]
                        if trk is None and "tracking_id" in grp:
                            trk = grp["tracking_id"]
                    except Exception:
                        continue
                    if ctx is not None and trk is not None:
                        break

            # Get flowcell from context_tags first, fall back to tracking_id
            flowcell = ""
            kit = ""
            sample_rate = 0

            if ctx is not None:
                flowcell = _decode_attr(ctx.attrs.get("flowcell_type", "")).upper()
                kit = _decode_attr(ctx.attrs.get("sequencing_kit", "")).upper()
                # Older MinKNOW versions used experiment_kit instead of sequencing_kit
                if not kit:
                    kit = _decode_attr(ctx.attrs.get("experiment_kit", "")).upper()
                rate_str = _decode_attr(ctx.attrs.get("sample_frequency", "0"))
                sample_rate = int(rate_str) if rate_str.isdigit() else 0

            # Fallback to tracking_id for flowcell product code
            if not flowcell and trk is not None:
                flowcell = _decode_attr(trk.attrs.get("flow_cell_product_code", "")).upper()
            # Also try tracking_id for sample_frequency if not found
            if sample_rate == 0 and trk is not None:
                rate_str = _decode_attr(trk.attrs.get("sample_frequency", "0"))
                sample_rate = int(rate_str) if rate_str.isdigit() else 0

            # Last resort: get sampling_rate from channel_id group (older
            # multi-read fast5 files that lack context_tags/tracking_id
            # entirely).  channel_id/sampling_rate is a float (e.g. 4000.0).
            if sample_rate == 0:
                for key in list(f.keys())[:5]:
                    try:
                        grp = f[key]
                        if "channel_id" in grp:
                            raw_rate = grp["channel_id"].attrs.get("sampling_rate", 0)
                            sample_rate = int(float(raw_rate))
                            if sample_rate > 0:
                                break
                    except Exception:
                        continue

            if not flowcell and not kit and sample_rate == 0:
                return None
            return {"flowcell": flowcell, "kit": kit, "sample_rate": sample_rate}
    except Exception:
        return None


def extract_chemistry_pod5(file_path: Path) -> dict | None:
    """Extract flowcell chemistry metadata from a pod5 file.

    Reads the RunInfo from the first read.  Checks direct fields first,
    then falls back to context_tags and tracking_id dicts (which may
    retain metadata from the original fast5 conversion).

    Returns {"flowcell": "FLO-MIN114", "kit": "SQK-LSK114",
    "sample_rate": 5000} or None.
    """
    if not HAS_POD5:
        return None
    try:
        with pod5.Reader(file_path) as reader:
            for read in reader.reads():
                info = read.run_info
                flowcell = (info.flow_cell_product_code or "").upper()
                kit = (info.sequencing_kit or "").upper()
                sample_rate = int(info.sample_rate or 0)

                # Fallback: check context_tags dict (may hold fast5 metadata)
                ctx = info.context_tags or {}
                if not flowcell:
                    flowcell = ctx.get("flowcell_type", "").upper()
                if not kit:
                    kit = ctx.get("sequencing_kit", "").upper()
                if not kit:
                    kit = ctx.get("experiment_kit", "").upper()
                if sample_rate == 0:
                    rate_str = ctx.get("sample_frequency", "0")
                    sample_rate = int(rate_str) if rate_str.isdigit() else 0

                # Fallback: check tracking_id dict
                trk = info.tracking_id or {}
                if not flowcell:
                    flowcell = trk.get("flow_cell_product_code", "").upper()

                if not flowcell and not kit and sample_rate == 0:
                    return None
                return {
                    "flowcell": flowcell,
                    "kit": kit,
                    "sample_rate": sample_rate,
                }
        return None
    except Exception:
        return None


def extract_chemistry(file_path: Path, file_format: str) -> dict | None:
    """Extract flowcell chemistry from a data file.

    Dispatches to the pod5 or fast5 extractor based on file_format.
    Returns None on any failure -- chemistry extraction is best-effort.
    """
    try:
        if file_format == "pod5":
            return extract_chemistry_pod5(file_path)
        elif file_format in ("fast5", "single_read_fast5", "multi_read_fast5"):
            return extract_chemistry_fast5(file_path)
    except Exception:
        pass
    return None


# Explicit dorado model names, for runs where automatic model selection
# ("dorado basecaller sup ...") cannot work (fast5 input, or flowcell/kit/
# sample-rate combinations missing from dorado's chemistry table).
# Source: dorado 0.9.6 model list (last release with R9.4.1/4 kHz/RNA002).
DORADO_LEGACY_MODELS = {
    "R9.4.1": "dna_r9.4.1_e8_sup@v3.6",
    "R10.4.1_4khz": "dna_r10.4.1_e8.2_400bps_sup@v4.1.0",
    "R10.4.1_260bps": "dna_r10.4.1_e8.2_260bps_sup@v4.1.0",
    "RNA002": "rna002_70bps_hac@v3",  # no sup model exists for RNA002
}


def classify_chemistry(chemistry: dict) -> dict:
    """Map raw chemistry metadata to pore type and dorado version recommendation.

    Input: {"flowcell": "FLO-MIN114", "kit": "SQK-LSK114", "sample_rate": 5000}
    Output: {"pore", "analyte", "dorado_version", "model_hint",
             "explicit_model", "note"}

    Support matrix (verified against the dorado changelog and the chemistry
    table in dorado 0.9.6's dorado/models/kits.cpp):
      - dorado >=1.0: DNA R10.4.1 E8.2 at 5 kHz and RNA004 only; POD5 input
        only (FAST5 support was removed in 1.0.0).
      - dorado 0.9.6 (last 0.x): additionally R9.4.1 E8 (4 kHz), R10.4.1 at
        4 kHz (400bps and 260bps), and RNA002.  Reads FAST5, but automatic
        model selection needs POD5 with a recognised flowcell/kit pair.
      - No dorado release has models for R10.3 (FLO-MIN111), R10.4
        (FLO-MIN112/PRO112, Kit 12) or DNA on R9.5.1 (FLO-MIN107); those
        need the discontinued Guppy basecaller.

    dorado_version is None when no dorado release can basecall the data.
    """
    flowcell = chemistry.get("flowcell", "")
    kit = chemistry.get("kit", "")
    sample_rate = chemistry.get("sample_rate", 0)

    pore = FLOWCELL_PORE.get(flowcell, "unknown")
    # Fallback: infer pore type from kit when flowcell code is absent
    if pore == "unknown" and kit:
        pore = KIT_PORE.get(kit, "unknown")
    # Last resort: infer from sample rate when both flowcell and kit are truly
    # absent (empty strings, not unrecognized codes).  Only runs with zero
    # metadata reach this path -- overwhelmingly pre-2020 R9.4.1 data.
    # PromethION R9.4.1 used 6 kHz sampling in 2018-2019, so any rate other
    # than exactly 5000 Hz (exclusive to R10.4.1) is assumed R9.4.1.
    if pore == "unknown" and not flowcell and not kit and sample_rate > 0:
        if sample_rate == 5000:
            pore = "R10.4.1"
        else:
            pore = "R9.4.1"
    analyte = "rna" if kit.startswith("SQK-RNA") else "dna"

    dorado_version = None
    model_hint = "sup"
    explicit_model = None
    note = None
    # Verified with dorado 0.9.6: it reads multi-read fast5 only, and only with
    # an explicit model path; single-read fast5 fails ("Unable to open the
    # group channel_id") and automatic model selection needs pod5.
    legacy_input = ("automatic model selection needs pod5 input; dorado 0.9.6 reads "
                    "multi-read fast5 only with an explicit model, and cannot read "
                    "single-read fast5")

    if kit.startswith("SQK-RNA004") or pore == "RNA004":
        analyte = "rna"
        dorado_version = ">=1.0"
    elif kit.startswith("SQK-RNA002"):
        analyte = "rna"
        dorado_version = "0.9.6"
        model_hint = "hac"
        explicit_model = DORADO_LEGACY_MODELS["RNA002"]
        note = ("RNA002 support was dropped in dorado 1.0; use dorado 0.9.6 "
                "(no sup model exists for RNA002)")
    elif pore == "R9.4.1":
        dorado_version = "0.9.6"
        explicit_model = DORADO_LEGACY_MODELS["R9.4.1"]
        note = "R9.4.1 support was dropped in dorado 1.0; use dorado 0.9.6"
        if sample_rate and sample_rate != 4000:
            note += (f"; dorado's only R9.4.1 model expects 4 kHz but this run "
                     f"was sampled at {sample_rate} Hz -- check results carefully")
    elif pore == "R10.4.1":
        is_260bps = kit.endswith("-260")
        if sample_rate == 5000 and not is_260bps:
            dorado_version = ">=1.0"
        elif sample_rate == 4000 or is_260bps:
            dorado_version = "0.9.6"
            explicit_model = DORADO_LEGACY_MODELS[
                "R10.4.1_260bps" if is_260bps else "R10.4.1_4khz"]
            note = ("R10.4.1 at 4 kHz requires dorado 0.9.6; only 5 kHz data "
                    "is supported in >=1.0")
        else:
            dorado_version = ">=1.0"
            note = ("Sample rate unknown: dorado >=1.0 supports R10.4.1 at "
                    "5 kHz only; use dorado 0.9.6 if the run was sampled at 4 kHz")
    elif pore in ("R10.4", "R10.3", "R9.5.1"):
        names = {"R10.4": "R10.4 (FLO-MIN112/FLO-PRO112, Kit 12)",
                 "R10.3": "R10.3 (FLO-MIN111)",
                 "R9.5.1": "DNA on R9.5.1 (FLO-MIN107)"}
        note = (f"No dorado release supports {names[pore]}; basecall with the "
                f"legacy Guppy basecaller (discontinued)")
    # pore == "unknown" -> dorado_version stays None

    if dorado_version == "0.9.6" and note:
        note += f". {legacy_input}"

    return {
        "pore": pore,
        "analyte": analyte,
        "dorado_version": dorado_version,
        "model_hint": model_hint if dorado_version else None,
        "explicit_model": explicit_model,
        "note": note,
    }


def dorado_advice_lines(chem_class: dict, formats: list[str]) -> list[str]:
    """Human-readable dorado recommendation for one run (verbose output).

    Checked against real dorado 0.9.6 runs: a model complex (fast/hac/sup)
    is resolved from pod5 RunInfo only; an explicit model name is treated as
    a folder path, so it must be downloaded first; single-read fast5 cannot
    be read by any dorado version.
    """
    lines = []
    ver = chem_class.get("dorado_version")
    note = chem_class.get("note")
    if not ver:
        if note:
            lines.append(f"dorado: {note}")
        return lines
    has_pod5 = "pod5" in formats
    lines.append(f"dorado: use version {ver}")
    if note:
        lines.append(f"        {note}")
    if not has_pod5 and "single_read_fast5" in formats:
        lines.append("        single-read fast5 cannot be read by dorado: convert to pod5 "
                     "first (--convert-to pod5 or nanopore_converter.py)")
    elif not has_pod5 and ver == ">=1.0" and any(
            f in formats for f in ("multi_read_fast5", "fast5_unknown")):
        lines.append("        dorado >=1.0 cannot read fast5: convert to pod5 first "
                     "(--convert-to pod5 or nanopore_converter.py)")
    lines.append(f"model: dorado basecaller {chem_class['model_hint']} <pod5_dir>/ > calls.bam")
    model = chem_class.get("explicit_model")
    if model:
        # dorado treats a model *name* as a folder path, so download it first
        # (into the current folder) and run basecaller from the same folder.
        lines.append("       if auto-selection fails (missing flowcell/kit metadata), "
                     "use the explicit model:")
        lines.append(f"         dorado download --model {model}")
        lines.append(f"         dorado basecaller {model} <pod5_dir>/ > calls.bam")
    return lines


def count_files_recursive(directory: Path, ext: str) -> tuple[int, int]:
    """Count files with given extension recursively using fast os.scandir.

    Returns (file_count, total_size_bytes).
    """
    return fast_count_files(directory, ext, recursive=True)


def compute_dir_size(directory: Path) -> int:
    """Recursively sum the size of all files in a directory tree via os.scandir."""
    total = 0

    def _walk(path: str) -> None:
        nonlocal total
        try:
            with os.scandir(path) as it:
                for entry in it:
                    if entry.is_file(follow_symlinks=False):
                        try:
                            total += entry.stat(follow_symlinks=False).st_size
                        except OSError:
                            pass
                    elif entry.is_dir(follow_symlinks=False):
                        _walk(entry.path)
        except PermissionError:
            pass

    _walk(str(directory))
    return total


def format_size(size_bytes: int) -> str:
    """Format a byte count as a human-readable string (e.g. '1.2 GB')."""
    if size_bytes < 1024:
        return f"{size_bytes} B"
    for unit in ("KB", "MB", "GB", "TB"):
        size_bytes /= 1024
        if size_bytes < 1024 or unit == "TB":
            return f"{size_bytes:.1f} {unit}"
    return f"{size_bytes:.1f} TB"


def _find_first_pod5(directory: Path) -> Path | None:
    """Find the first .pod5 file in a directory, descending into one level of subdirs.

    Barcoded runs store pod5 files inside barcode subdirectories
    (pod5_pass/barcode01/*.pod5), so a flat scan of the pod5 directory
    may yield only subdirectories. This helper checks the top level first,
    then peeks into immediate subdirectories.
    """
    try:
        with os.scandir(directory) as it:
            subdirs = []
            for entry in it:
                if entry.is_file(follow_symlinks=False) and entry.name.endswith(".pod5"):
                    return Path(entry.path)
                if entry.is_dir(follow_symlinks=False):
                    subdirs.append(entry.path)
        for subdir in subdirs:
            try:
                with os.scandir(subdir) as it:
                    for entry in it:
                        if entry.is_file(follow_symlinks=False) and entry.name.endswith(".pod5"):
                            return Path(entry.path)
            except PermissionError:
                pass
    except PermissionError:
        pass
    return None


def analyze_run(run_path: Path, quick: bool = False) -> dict:
    """
    Analyze a nanopore run directory and determine its read format(s).

    When quick=True, skip file counting and size calculation -- only report
    the detected format(s).

    Detection logic:
        - pod5:             pod5/ subfolder 1-2 levels down
        - multi_read_fast5: fast5/ subfolder 1-2 levels down, files contain multiple reads
        - single_read_fast5: .fast5 files in run root or in numeric subfolders (0/, 1/, 2/, ...)
                             each file contains a single read

    Returns dict with:
        - formats: list of detected formats
        - details: dict with per-format info (paths, file counts)
    """
    result = {"formats": [], "details": {}, "chemistry": None, "chemistry_classification": None, "run_path": str(run_path)}

    # --- Check permissions on run directory and its subdirectories ---
    # Use os.scandir to extract only subdirectories, skipping all files.
    # This avoids materializing 100k+ Path objects in single-read fast5 dirs.
    permission_errors = []
    subdirs = []
    try:
        with os.scandir(run_path) as it:
            for entry in it:
                if entry.is_dir(follow_symlinks=False):
                    subdirs.append(Path(entry.path))
    except PermissionError:
        result["formats"].append("permission_denied")
        result["details"]["permission_denied"] = {
            "reasons": [f"Cannot read run directory: {run_path}"],
        }
        return result

    for sd in subdirs:
        if not _is_dir_readable(sd):
            permission_errors.append(sd.name)

    if permission_errors and not any(
        sd.name not in permission_errors for sd in subdirs
    ):
        # ALL subdirectories are unreadable
        result["formats"].append("permission_denied")
        result["details"]["permission_denied"] = {
            "reasons": [f"All subdirectories are unreadable ({len(permission_errors)} dir(s))"],
            "inaccessible_dirs": permission_errors,
        }
        return result

    # --- Single-pass directory discovery ---
    structure = discover_run_structure(run_path)
    pod5_dirs = structure["pod5"]
    pod5_variant_dirs = structure["pod5_prefix"]
    fast5_dirs = structure["fast5"]
    fast5_variant_dirs = structure["fast5_prefix"]
    fastq_dirs = structure["fastq_prefix"]

    # --- Check for pod5 (pod5/ or pod5_pass/pod5_fail/pod5_skip, up to 5 levels down) ---
    all_pod5_dirs = pod5_dirs + pod5_variant_dirs
    if all_pod5_dirs:
        variant_names = sorted(set(d.name for d in all_pod5_dirs))
        unreadable_pod5 = [str(d) for d in all_pod5_dirs if not _is_dir_readable(d)]
        pod5_detail = {
            "directories": [str(d) for d in all_pod5_dirs],
            "folder_variants": variant_names,
        }
        if not quick:
            counts = [count_files_recursive(d, ".pod5") for d in all_pod5_dirs]
            total_pod5 = sum(c for c, _ in counts)
            total_pod5_size = sum(s for _, s in counts)
            pod5_detail["file_count"] = total_pod5
            pod5_detail["data_size_bytes"] = total_pod5_size
        if unreadable_pod5:
            pod5_detail["note"] = f"Permission denied on {len(unreadable_pod5)} pod5 dir(s)"
            pod5_detail["inaccessible_dirs"] = unreadable_pod5
        elif not quick and total_pod5 == 0:
            pod5_detail["note"] = "Empty pod5 folder(s) found"
        result["formats"].append("pod5")
        result["details"]["pod5"] = pod5_detail

        # Extract chemistry from first pod5 file
        if result["chemistry"] is None:
            for d in all_pod5_dirs:
                if str(d) in unreadable_pod5:
                    continue
                pod5_file = _find_first_pod5(d)
                if pod5_file:
                    chem = extract_chemistry(pod5_file, "pod5")
                    if chem:
                        result["chemistry"] = chem
                        result["chemistry_classification"] = classify_chemistry(chem)
                        break

    # --- Check for fast5 (fast5/ or fast5_skip/fast5_pass/fast5_fail, up to 5 levels down) ---
    all_fast5_dirs = fast5_dirs + fast5_variant_dirs
    if all_fast5_dirs:
        unreadable_fast5 = [str(d) for d in all_fast5_dirs if not _is_dir_readable(d)]

        # Sample a few fast5 files (descending into barcode/numeric subfolders)
        readable_fast5_dirs = [d for d in all_fast5_dirs if str(d) not in unreadable_fast5]
        samples = find_fast5_samples(readable_fast5_dirs)
        classification, sample_file = classify_fast5_samples(samples)
        fast5_dir_set = {str(d) for d in all_fast5_dirs}
        has_subdirs = any(str(p.parent) not in fast5_dir_set for p in samples)

        variant_names = sorted(set(d.name for d in all_fast5_dirs))
        fast5_info = {
            "directories": [str(d) for d in all_fast5_dirs],
            "folder_variants": variant_names,
        }
        if unreadable_fast5:
            fast5_info["inaccessible_dirs"] = unreadable_fast5

        if sample_file:
            fast5_info["sampled_file"] = sample_file.name

            if classification == "single_read_fast5":
                if not quick:
                    counts = [estimate_dir_size(d, ".fast5", recursive=True) for d in all_fast5_dirs]
                    total_count = sum(c for c, _, _ in counts)
                    total_size = sum(s for _, s, _ in counts)
                    is_estimated = any(e for _, _, e in counts)
                    fast5_info["file_count"] = total_count
                    fast5_info["data_size_bytes"] = total_size
                    if is_estimated:
                        fast5_info["size_estimated"] = True
                if has_subdirs:
                    fast5_info["note"] = "fast5/ contains subdirs with single-read files"
                else:
                    fast5_info["note"] = "fast5/ contains single-read files"
                result["formats"].append("single_read_fast5")
                result["details"]["single_read_fast5"] = fast5_info
            elif classification == "multi_read_fast5":
                # Multi-read -- count is fast (few large files)
                if not quick:
                    counts = [count_files_recursive(d, ".fast5") for d in all_fast5_dirs]
                    total_fast5 = sum(c for c, _ in counts)
                    total_fast5_size = sum(s for _, s in counts)
                    fast5_info["file_count"] = total_fast5
                    fast5_info["data_size_bytes"] = total_fast5_size
                result["formats"].append("multi_read_fast5")
                result["details"]["multi_read_fast5"] = fast5_info
            else:
                result["formats"].append("fast5_unknown")
                result["details"]["fast5_unknown"] = fast5_info
                if not HAS_H5PY:
                    fast5_info["note"] = "Install h5py to distinguish single/multi read fast5"

            if sample_file and result["chemistry"] is None:
                chem = extract_chemistry(sample_file, "fast5")
                if chem:
                    result["chemistry"] = chem
                    result["chemistry_classification"] = classify_chemistry(chem)
        elif unreadable_fast5:
            # Found fast5 dirs but can't read into them — can't determine format without access
            fast5_info["note"] = f"Permission denied on {len(unreadable_fast5)} fast5 dir(s) — format unknown"
            result["formats"].append("fast5_unknown")
            result["details"]["fast5_unknown"] = fast5_info
        else:
            # Found fast5 folders but they're empty — can't determine single vs multi without sample files
            fast5_info["note"] = f"Empty fast5 folder(s) found — format unknown"
            result["formats"].append("fast5_unknown")
            result["details"]["fast5_unknown"] = fast5_info

    # --- Check for fastq (fastq_pass/fastq_fail, up to 5 levels down) ---
    if fastq_dirs:
        unreadable_fastq = [str(d) for d in fastq_dirs if not _is_dir_readable(d)]
        variant_names = sorted(set(d.name for d in fastq_dirs))
        fastq_detail = {
            "directories": [str(d) for d in fastq_dirs],
            "folder_variants": variant_names,
        }
        if not quick:
            _fastq_exts = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
            counts = [fast_count_files(d, extensions=_fastq_exts, recursive=True) for d in fastq_dirs]
            total_fastq = sum(c for c, _ in counts)
            total_fastq_size = sum(s for _, s in counts)
            fastq_detail["file_count"] = total_fastq
            fastq_detail["data_size_bytes"] = total_fastq_size
        if unreadable_fastq:
            fastq_detail["note"] = f"Permission denied on {len(unreadable_fastq)} fastq dir(s)"
            fastq_detail["inaccessible_dirs"] = unreadable_fastq
        elif not quick and total_fastq == 0:
            fastq_detail["note"] = "Empty fastq folder(s) found"
        result["formats"].append("fastq")
        result["details"]["fastq"] = fastq_detail

    # --- Check for single_read_fast5 (.fast5 in run root or numeric subdirs 0/, 1/, 2/, ...) ---
    # Only check if we didn't already find fast5 via a fast5/ subfolder
    # Uses os.scandir for speed — single_fast5 dirs can have 100k+ files
    if not all_fast5_dirs:
        root_files = []
        numeric_dirs = []
        try:
            with os.scandir(run_path) as it:
                for entry in it:
                    try:
                        if entry.is_file() and entry.name.lower().endswith(".fast5"):
                            if len(root_files) < 5:
                                root_files.append(Path(entry.path))
                        elif entry.is_dir() and re.match(r"^\d+$", entry.name):
                            numeric_dirs.append(Path(entry.path))
                    except OSError:
                        continue
        except PermissionError:
            pass
        # Every numeric folder is a candidate, not just the first one listed
        samples = root_files or find_fast5_samples(sorted(numeric_dirs), max_depth=2)
        classification, sample_file = classify_fast5_samples(samples)

        if sample_file:
            fast5_info = {
                "directories": [str(run_path)],
                "sampled_file": sample_file.name,
                "layout": "root" if root_files else "numeric_subdirs",
            }
            if not quick:
                if classification == "multi_read_fast5":
                    count, size = fast_count_files(run_path, ".fast5", recursive=True)
                    is_estimated = False
                else:
                    count, size, is_estimated = estimate_dir_size(
                        run_path, ".fast5", recursive=True
                    )
                fast5_info["file_count"] = count
                fast5_info["data_size_bytes"] = size
                if is_estimated:
                    fast5_info["size_estimated"] = True

            if result["chemistry"] is None:
                chem = extract_chemistry(sample_file, "fast5")
                if chem:
                    result["chemistry"] = chem
                    result["chemistry_classification"] = classify_chemistry(chem)

            if classification in ("single_read_fast5", "multi_read_fast5"):
                result["formats"].append(classification)
                result["details"][classification] = fast5_info
            else:
                result["formats"].append("fast5_unknown")
                result["details"]["fast5_unknown"] = fast5_info

    # --- Check for compressed archives (tar, gz, tar.gz) -- treat as single_read_fast5 ---
    # Only when no raw data (pod5/fast5) was found: a pod5 or fast5 run that
    # merely contains e.g. sequencing_summary.txt.gz is not an archive run.
    raw_formats = ("pod5", "single_read_fast5", "multi_read_fast5", "fast5_unknown")
    if not all_fast5_dirs and not any(f in result["formats"] for f in raw_formats):
        archives = []
        try:
            with os.scandir(run_path) as it:
                for entry in it:
                    if entry.is_file(follow_symlinks=False) and is_raw_archive(entry.name):
                        archives.append(Path(entry.path))
        except PermissionError:
            pass

        if archives:
            result["formats"].append("single_read_fast5")
            result["details"]["single_read_fast5"] = {
                "directories": [str(run_path)],
                "file_count": 0,
                "archive_files": [f.name for f in archives],
                "archive_count": len(archives),
                "note": "Compressed archives found — assumed single-read fast5",
            }

    if not result["formats"]:
        # Gather diagnostic info about what IS in the directory
        diag = diagnose_unknown(run_path)
        result["formats"].append("unknown")
        result["details"]["unknown"] = diag

    return result


def _has_file_with_ext(directory: Path, ext: str, max_depth: int = 5,
                       max_entries: int = 10_000) -> bool:
    """Check if any file with the given extension exists, with depth and iteration limits."""
    budget = max_entries

    def _search(path: str, depth: int) -> bool:
        nonlocal budget
        if depth > max_depth or budget <= 0:
            return False
        try:
            with os.scandir(path) as it:
                for entry in it:
                    budget -= 1
                    if budget <= 0:
                        return False
                    if entry.is_file(follow_symlinks=False) and entry.name.endswith(ext):
                        return True
                    if entry.is_dir(follow_symlinks=False) and _search(entry.path, depth + 1):
                        return True
        except PermissionError:
            pass
        return False
    return _search(str(directory), 1)


def diagnose_unknown(run_path: Path) -> dict:
    """Gather diagnostic info for a run directory that matched no known format."""
    diag = {"reasons": []}

    # Separate files and directories in a single os.scandir pass
    MAX_DIAG_ENTRIES = 5_000
    subdirs = []
    file_count = 0
    extensions = {}
    truncated = False
    try:
        entry_count = 0
        with os.scandir(run_path) as it:
            for entry in it:
                entry_count += 1
                if entry_count >= MAX_DIAG_ENTRIES:
                    truncated = True
                    break
                try:
                    is_dir = entry.is_dir()
                    is_file = not is_dir and entry.is_file()
                except OSError:
                    is_dir = is_file = False
                if is_dir:
                    subdirs.append(entry.name)
                elif is_file:
                    file_count += 1
                    dot_pos = entry.name.rfind(".")
                    ext = entry.name[dot_pos:].lower() if dot_pos >= 0 else "(no extension)"
                    extensions[ext] = extensions.get(ext, 0) + 1
    except PermissionError:
        diag["reasons"].append("Permission denied reading directory")
        return diag

    if truncated:
        diag["reasons"].append(f"Sampling stopped at {MAX_DIAG_ENTRIES} entries (directory may contain more)")

    if not subdirs and file_count == 0:
        diag["reasons"].append("Directory is empty")
        return diag

    diag["subdirectory_names"] = subdirs
    diag["file_count"] = file_count
    diag["file_extensions"] = extensions

    # Build human-readable reasons
    if not subdirs and file_count > 0:
        if extensions:
            ext_summary = ", ".join(f"{count}x {ext}" for ext, count in sorted(extensions.items()))
            diag["reasons"].append(f"No subdirectories; {file_count} file(s) in root: {ext_summary}")
        else:
            diag["reasons"].append(f"No subdirectories; {file_count} file(s) with no recognized extensions")

        if not any(ext in extensions for ext in (".fast5", ".pod5")):
            diag["reasons"].append("No .fast5 or .pod5 files found")
    else:
        dir_names = ", ".join(subdirs[:10])
        diag["reasons"].append(f"Subdirectories found: {dir_names}")

        # Check which subdirs are readable
        unreadable = []
        readable_subdirs = []
        for sd_name in subdirs:
            sd_path = run_path / sd_name
            if _is_dir_readable(sd_path):
                readable_subdirs.append(sd_path)
            else:
                unreadable.append(sd_name)

        if unreadable:
            diag["inaccessible_dirs"] = unreadable
            diag["reasons"].append(
                f"Permission denied on {len(unreadable)} subdir(s): {', '.join(unreadable[:5])}"
                + (f" (+{len(unreadable)-5} more)" if len(unreadable) > 5 else "")
            )

        if not readable_subdirs and unreadable:
            diag["reasons"].append("All subdirectories are unreadable -- cannot determine format")
            return diag

        # Check if readable subdirs contain fast5/pod5 deeper than we searched.
        # Uses bounded os.scandir walk instead of unbounded rglob.
        has_deep_fast5 = any(_has_file_with_ext(sd, ".fast5") for sd in readable_subdirs[:10])
        has_deep_pod5 = any(_has_file_with_ext(sd, ".pod5") for sd in readable_subdirs[:10])

        if has_deep_fast5:
            diag["reasons"].append("Found .fast5 files deeper in tree but not in expected locations (fast5/, root, or numeric subdirs)")
        if has_deep_pod5:
            diag["reasons"].append("Found .pod5 files deeper in tree but not in expected pod5/ subdirectory")
        if not has_deep_fast5 and not has_deep_pod5:
            diag["reasons"].append("No .fast5 or .pod5 files found anywhere in directory tree")

    return diag


def print_conversion_help(fmt: str):
    """Print conversion instructions for a given format."""
    if fmt == "multi_read_fast5":
        print("  Conversion options:")
        print("    -> To pod5:         pod5 convert fast5 ./run_dir/ --output ./pod5_dir/ --threads 20 --recursive")
        print("                        # https://github.com/nanoporetech/pod5-file-format")
        print("    -> To single_fast5: multi_to_single_fast5 --input_path multi/ --save_path single/ --recursive")
        print("                        # https://github.com/nanoporetech/ont_fast5_api")
    elif fmt == "single_read_fast5":
        print("  Conversion options:")
        print("    -> To pod5 (two steps required):")
        print("       1. single_to_multi_fast5 -i run_dir/ -s multi_dir/ -t 4 --recursive")
        print("       2. pod5 convert fast5 multi_dir/ --output pod5_dir/ --threads 20 --recursive")
        print("       # single-read fast5 must first be merged to multi-read format")
        print("    -> To multi_fast5:  single_to_multi_fast5 -i single_dir/ -s multi_dir/ -t 4 --recursive")
    elif fmt == "pod5":
        print("  Conversion options:")
        print("    -> To fast5:        pod5 convert to_fast5 input.pod5 --output ./fast5_dir/")
        print("    -> Inspect:         pod5 inspect reads input.pod5")
        print("    -> Merge:           pod5 merge input1.pod5 input2.pod5 --output merged.pod5")


def _pod5_metadata_fix_lines(pod5_dir: str, chemistry: dict) -> list[str]:
    """Generate bash lines that patch pod5 metadata using the pod5 Python API.

    Rewrites each pod5 file in pod5_dir, injecting the detected flowcell,
    kit, and sample_rate into RunInfo.  This is necessary for old fast5 runs
    where the converter produces pod5 files with empty or incorrect metadata.

    pod5.RunInfo is a frozen dataclass, so a patched copy is built with
    dataclasses.replace() (once per acquisition) and assigned to each Read.
    """
    flowcell = chemistry.get("flowcell", "")
    kit = chemistry.get("kit", "")
    sample_rate = chemistry.get("sample_rate", 0)
    if not flowcell and not kit and not sample_rate:
        return []
    updates = {}
    if flowcell:
        updates["flow_cell_product_code"] = flowcell
    if kit:
        updates["sequencing_kit"] = kit
    if sample_rate:
        updates["sample_rate"] = int(sample_rate)
    summary = f"flowcell={flowcell}, kit={kit}, sample_rate={sample_rate}"
    lines = []
    lines.append(f"echo {_sh_quote(f'  Patching pod5 metadata ({summary})...')}")
    # Quoted here-document: bash performs no expansion inside it.  Values are
    # embedded with repr() so arbitrary paths/codes are valid Python literals.
    lines.append("python3 << 'PATCH_EOF'")
    lines.append("import dataclasses, os, sys")
    lines.append("from pathlib import Path")
    lines.append("try:")
    lines.append("    import pod5")
    lines.append("except ImportError:")
    lines.append("    print('Warning: pod5 Python package not installed, skipping metadata patch')")
    lines.append("    sys.exit(0)")
    lines.append(f"pod5_dir = Path({pod5_dir!r})")
    lines.append(f"updates = {updates!r}")
    lines.append("files = sorted(pod5_dir.rglob('*.pod5'))")
    lines.append("if not files:")
    lines.append("    print('  No pod5 files found to patch')")
    lines.append("    sys.exit(0)")
    lines.append("def matches(info):")
    lines.append("    for k, v in updates.items():")
    lines.append("        cur = getattr(info, k)")
    lines.append("        if (str(cur).upper() != str(v).upper()) if isinstance(v, str) else cur != v:")
    lines.append("            return False")
    lines.append("    return True")
    lines.append("failed = 0")
    lines.append("for p5 in files:")
    lines.append("    tmp = p5.with_name(p5.name + '.tmp')")
    lines.append("    patched = {}")
    lines.append("    try:")
    lines.append("        # Skip files that already carry these values (case-insensitive,")
    lines.append("        # as dorado compares them): rewriting a large pod5 is expensive.")
    lines.append("        with pod5.Reader(p5) as reader:")
    lines.append("            if all(matches(r.run_info) for r in reader.reads()):")
    lines.append("                print(f'  Already correct: {p5.name}')")
    lines.append("                continue")
    lines.append("        if tmp.exists():")
    lines.append("            tmp.unlink()")
    lines.append("        with pod5.Reader(p5) as reader, pod5.Writer(tmp) as writer:")
    lines.append("            for record in reader.reads():")
    lines.append("                read = record.to_read()")
    lines.append("                key = read.run_info.acquisition_id")
    lines.append("                if key not in patched:")
    lines.append("                    patched[key] = dataclasses.replace(read.run_info, **updates)")
    lines.append("                read.run_info = patched[key]")
    lines.append("                writer.add_read(read)")
    lines.append("        os.replace(tmp, p5)")
    lines.append("        print(f'  Patched: {p5.name}')")
    lines.append("    except Exception as e:")
    lines.append("        failed += 1")
    lines.append("        print(f'  Warning: failed to patch {p5.name}: {e}')")
    lines.append("        if tmp.exists():")
    lines.append("            tmp.unlink()")
    lines.append("if failed:")
    lines.append("    print(f'  Warning: {failed} pod5 file(s) could not be patched')")
    lines.append("PATCH_EOF")
    return lines


def _verify_read_count_lines(inputs: list[str], pod5_dir: str, exclude: str) -> list[str]:
    """Bash lines that fail the script unless pod5_dir holds every read found
    in the fast5 files under `inputs` (ignoring the `exclude` work folder).

    pod5 convert and single_to_multi_fast5 exit 0 after skipping files they
    cannot read (mixed single/multi-read runs, truncated files), so without
    this check reads would be lost silently.
    """
    return [
        "python3 - << 'VERIFY_EOF'",
        "import os, sys",
        "from pathlib import Path",
        "try:",
        "    import h5py, pod5",
        "except ImportError:",
        "    print('  Warning: h5py/pod5 Python packages missing; read count not verified')",
        "    sys.exit(0)",
        f"inputs = {inputs!r}",
        f"exclude = os.path.realpath({exclude!r})",
        "expected = unreadable = 0",
        "for top in inputs:",
        "    for root, dirs, files in os.walk(top):",
        "        dirs[:] = [d for d in dirs if os.path.realpath(os.path.join(root, d)) != exclude]",
        "        for name in files:",
        "            if not name.lower().endswith('.fast5'):",
        "                continue",
        "            try:",
        "                with h5py.File(os.path.join(root, name), 'r') as f:",
        "                    n = sum(1 for k in f if k.startswith('read_'))",
        "                    expected += n if n else (1 if ('Raw' in f or 'UniqueGlobalKey' in f) else 0)",
        "            except Exception:",
        "                unreadable += 1",
        f"got = sum(pod5.Reader(p).num_reads for p in Path({pod5_dir!r}).rglob('*.pod5'))",
        "if unreadable:",
        "    print(f'  Warning: {unreadable} unreadable fast5 file(s) skipped')",
        "if got != expected:",
        "    print(f'  ERROR: pod5 output has {got} reads but the fast5 input has {expected}.')",
        "    print('  Reads were lost (mixed single/multi-read files or corrupt files?).')",
        "    print('  Convert this run with nanopore_converter.py, which handles both.')",
        "    sys.exit(1)",
        "print(f'  Verified {got} reads')",
        "VERIFY_EOF",
    ]


def _guarded_block(run_name: str, body: list[str], cleanup: str | None = None) -> list[str]:
    """Run body in a fail-fast subshell; a failure is recorded and the script
    continues with the next run (exit status 1 at the end).

    The subshell's status is checked on a separate line on purpose: inside
    `( ... ) || ...` bash would silently disable `set -e`.
    """
    q = _sh_quote
    lines = ["(", "  set -e"]
    in_heredoc = False
    for line in body:  # heredoc bodies must stay unindented
        if in_heredoc:
            lines.append(line)
            if line in ("VERIFY_EOF", "PATCH_EOF"):
                in_heredoc = False
        else:
            lines.append("  " + line)
            if line.endswith(("<< 'VERIFY_EOF'", "<< 'PATCH_EOF'")):
                in_heredoc = True
    lines.append(")")
    on_fail = f"rm -rf {q(cleanup)}; " if cleanup else ""
    lines.append(f"if [ $? -ne 0 ]; then {on_fail}FAILED+=({q(run_name)}); "
                 f"echo {q(f'FAILED: {run_name}')} >&2; fi")
    return lines


def _pod5_conversion_lines(run_name: str, info: dict, fmt: str, detail: dict,
                           output_dir: str | None) -> list[str]:
    """Re-runnable bash block converting one run's fast5 to pod5.

    - Output: <output_dir>/<run>/pod5, or <run>/pod5_from_fast5 in place.  The
      in-place name is fixed (not chosen by what exists) so re-checking a
      converted run never schedules the same conversion again.
    - Skipped when the output already holds .pod5 files, so a failed run can
      be fixed and the script simply re-run.
    - Work happens in <base>/.convert_tmp and is moved into place only after
      the read count is verified; the input is only the detected fast5
      folders (or the run root for root/numeric-folder layouts).
    """
    q = _sh_quote
    run_path = info.get("run_path", "")
    inputs = detail.get("directories") or [run_path]
    if output_dir:
        base = os.path.join(output_dir, run_name)
        out = os.path.join(base, "pod5")
    else:
        base = run_path
        out = os.path.join(run_path, "pod5_from_fast5")
    work = os.path.join(base, ".convert_tmp")
    partial = os.path.join(work, "pod5")
    multi_tmp = os.path.join(work, "multi_fast5")

    lines = [f"# --- {run_name} ({fmt} -> pod5) ---"]
    lines.append(f"if compgen -G {q(out)}'/*.pod5' > /dev/null; then")
    lines.append(f"  echo {q(f'Skipping {run_name}: already converted in {out}')}")
    lines.append("else")
    steps = ", two steps" if fmt == "single_read_fast5" else ""
    body = [f"echo {q(f'Converting {run_name} ({fmt} -> pod5{steps})...')}",
            f"rm -rf {q(work)}",
            f"mkdir -p {q(partial)}"]
    if fmt == "single_read_fast5":
        body.append(f"mkdir -p {q(multi_tmp)}")
        for i, d in enumerate(inputs):
            # Distinct filename prefix per input folder so batches don't collide
            body.append(f"single_to_multi_fast5 -i {q(d)} -s {q(multi_tmp)} -f batch{i} -t 4 --recursive")
        pod5_inputs = [multi_tmp + "/"]
    else:
        pod5_inputs = [d.rstrip("/") + "/" for d in inputs]
    body.append("pod5 convert fast5 " + " ".join(q(i) for i in pod5_inputs)
                + f" --output {q(partial + '/')} --threads 20 --recursive")
    body.extend(_verify_read_count_lines(inputs, partial, work))
    chem = info.get("chemistry")
    if chem:
        body.extend(_pod5_metadata_fix_lines(partial, chem))
    body.append(f"mkdir -p {q(os.path.dirname(out))}")
    body.append(f"if [ -d {q(out)} ]; then rmdir {q(out)}; fi  # empty leftover (non-empty: fail)")
    body.append(f"mv {q(partial)} {q(out)}")
    body.append(f"rm -rf {q(work)}")
    lines.extend(_guarded_block(run_name, body, cleanup=work))
    lines.append("fi")
    lines.append("")
    return lines


def generate_conversion_script(runs: dict, target_format: str, output_dir: str | None = None):
    """Generate a bash conversion script.

    When output_dir is provided, all converted files are written under
    <output_dir>/<run_name>/<format>/ instead of alongside the originals.

    For runs where chemistry was detected from the source fast5, a
    post-conversion metadata patching step is included to inject the
    correct flowcell, kit, and sample_rate into the pod5 RunInfo.

    All paths are shell-quoted.  Runs whose data exists only as compressed
    archives are listed as comments (they must be extracted first).
    """
    q = _sh_quote
    lines = [
        "#!/usr/bin/env bash",
        "# Auto-generated nanopore format conversion script",
        f"# Target format: {target_format}",
        "# Safe to re-run: runs whose output already holds .pod5 files are skipped.",
        "# Each run executes fail-fast in its own subshell; a failed run is",
        "# reported and the script continues, exiting 1 at the end.",
        "set -uo pipefail",
        "FAILED=()",
        "",
    ]

    for run_name, info in runs.items():
        for fmt in info["formats"]:
            if fmt == target_format:
                continue  # Already in target format

            detail = info["details"].get(fmt, {})
            run_path = info.get("run_path", "")
            if fmt in ("single_read_fast5", "multi_read_fast5") and (
                detail.get("archive_files") or not run_path
            ):
                archives = detail.get("archive_files") or [run_name]
                lines.append(f"# SKIPPED {run_name}: data is in compressed archive(s) "
                             f"({', '.join(archives)}).")
                lines.append("#   Extract first, then re-run the format checker, or use "
                             "nanopore_converter.py which handles archives directly.")
                lines.append("")
                continue

            if target_format == "pod5" and fmt in ("multi_read_fast5", "single_read_fast5"):
                lines.extend(_pod5_conversion_lines(run_name, info, fmt, detail, output_dir))

            elif target_format == "single_fast5" and fmt == "multi_read_fast5":
                dirs = detail.get("directories") or [run_path]
                base = (os.path.join(output_dir, run_name, "single_fast5") if output_dir
                        else os.path.join(run_path, "single_fast5"))
                for d in dirs:
                    # One output folder per input folder: fast5_pass and
                    # fast5_fail are not mixed together.
                    sub = "root" if os.path.normpath(d) == os.path.normpath(run_path) \
                        else os.path.basename(os.path.normpath(d))
                    out = os.path.join(base, sub)
                    lines.extend(_guarded_block(f"{run_name}/{sub}", [
                        f"echo {q(f'Converting {run_name} ({fmt} -> single_fast5): {sub}')}",
                        f"mkdir -p {q(out)}",
                        f"multi_to_single_fast5 --input_path {q(d)} --save_path {q(out)} --recursive",
                    ]))
                    lines.append("")

    lines.append('if [ ${#FAILED[@]} -gt 0 ]; then')
    lines.append('  echo "${#FAILED[@]} run(s) failed: ${FAILED[*]}" >&2')
    lines.append('  exit 1')
    lines.append('fi')
    lines.append('echo "All conversions finished."')
    return "\n".join(lines)


def _tsv_field(value) -> str:
    """Make a value safe for a TSV cell (tabs/newlines would split rows)."""
    return re.sub(r"[\t\r\n]+", " ", str(value))


def write_stats_tsv(all_runs: dict, output_path: str) -> None:
    """Write per-run statistics to a TSV file.

    Produces one row per format per run, so a run containing both pod5 and
    fastq output generates two rows.

    Columns: run_name, format, file_count, data_size_bytes, size_estimated,
    directories, notes, flowcell_code, sequencing_kit, sample_rate, pore_type,
    dorado_version, dorado_model (explicit legacy model name, if any),
    dorado_note
    """
    header = ["run_name", "format", "file_count", "data_size_bytes",
              "size_estimated", "directories", "notes",
              "flowcell_code", "sequencing_kit", "sample_rate",
              "pore_type", "dorado_version", "dorado_model", "dorado_note"]
    with open(output_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for run_name, info in all_runs.items():
            for fmt in info["formats"]:
                detail = info["details"].get(fmt, {})
                file_count = detail.get("file_count", "")
                data_size = detail.get("data_size_bytes", "")
                estimated = detail.get("size_estimated", "")
                if estimated is True:
                    estimated = "True"
                elif estimated is False or estimated == "":
                    estimated = ""
                dirs = ";".join(detail.get("directories", []))
                # Aggregate notes from multiple fields
                notes_parts = []
                if detail.get("note"):
                    notes_parts.append(detail["note"])
                for reason in detail.get("reasons", []):
                    notes_parts.append(reason)
                for af in detail.get("archive_files", []):
                    notes_parts.append(f"archive: {af}")
                notes = "; ".join(notes_parts)
                chem = info.get("chemistry") or {}
                chem_class = info.get("chemistry_classification") or {}
                row = [
                    run_name,
                    fmt,
                    str(file_count),
                    str(data_size),
                    str(estimated),
                    dirs,
                    notes,
                    chem.get("flowcell", ""),
                    chem.get("kit", ""),
                    str(chem.get("sample_rate", "")) if chem.get("sample_rate") else "",
                    chem_class.get("pore", ""),
                    chem_class.get("dorado_version", "") or "",
                    chem_class.get("explicit_model", "") or "",
                    chem_class.get("note", "") or "",
                ]
                fh.write("\t".join(_tsv_field(v) for v in row) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Scan nanopore run folders and classify read formats (pod5, single/multi fast5)."
    )
    parser.add_argument(
        "target_folder",
        help="Path to folder containing nanopore run directories",
    )
    parser.add_argument(
        "--convert-to",
        choices=["pod5", "single_fast5"],
        help="Generate a bash conversion script for the specified target format",
    )
    parser.add_argument(
        "--script-output",
        default="convert_runs.sh",
        help="Output path for generated conversion script (default: convert_runs.sh)",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        metavar="DIR",
        help="Base output directory for converted files (default: alongside originals)",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed information per run",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Skip file counting and size calculation; only report format classification",
    )
    parser.add_argument(
        "--output-stats", "-o",
        default=None,
        metavar="FILE",
        help="Write per-run statistics to a TSV file",
    )
    args = parser.parse_args()

    # Absolute paths, so the generated script works from any directory
    target = Path(args.target_folder).resolve()
    if args.output_dir:
        args.output_dir = os.path.abspath(args.output_dir)
    if not target.is_dir():
        print(f"Error: '{target}' is not a directory", file=sys.stderr)
        sys.exit(1)

    if not HAS_H5PY:
        print("Warning: h5py not installed. Cannot distinguish single vs multi read fast5.",
              file=sys.stderr)
        print("  Install with: pip install h5py\n", file=sys.stderr)

    # Find nanopore run directories
    run_dirs = sorted(
        [d for d in target.iterdir() if d.is_dir() and is_nanopore_run_dir(d.name)]
    )

    # Find compressed archives in target folder — treat as single_read_fast5 runs
    archive_files = sorted(
        [f for f in target.iterdir()
         if f.is_file() and is_raw_archive(f.name)
         and is_nanopore_run_dir(f.name)]
    )

    if not run_dirs and not archive_files:
        print(f"No nanopore run directories found in '{target}'")
        print("  (Expected directories starting with YYYYMMDD_ pattern)")
        sys.exit(0)

    print(f"Found {len(run_dirs)} nanopore run(s) and {len(archive_files)} archive(s) in '{target}'\n")
    print(f"{'Run / Archive':<55} {'Format(s)':<25} {'Chemistry':<18} {'Data size':<12} {'Files'}")
    print("-" * 120)

    all_runs = {}
    format_summary = {"pod5": 0, "multi_read_fast5": 0, "single_read_fast5": 0,
                       "fastq": 0, "fast5_unknown": 0, "permission_denied": 0, "unknown": 0}

    for run_dir in run_dirs:
        result = analyze_run(run_dir, quick=args.quick)

        # Sum recognized format data sizes (unless --quick)
        if not args.quick:
            result["total_size_bytes"] = sum(
                detail.get("data_size_bytes", 0)
                for detail in result["details"].values()
            )

        all_runs[run_dir.name] = result

        formats_str = ", ".join(result["formats"])
        file_counts = []
        for fmt in result["formats"]:
            detail = result["details"].get(fmt, {})
            count = detail.get("file_count")
            if count:
                file_counts.append(f"{count}")
            format_summary[fmt] = format_summary.get(fmt, 0) + 1

        count_str = " / ".join(file_counts) if file_counts else "-"
        size_str = format_size(result["total_size_bytes"]) if "total_size_bytes" in result else "-"
        if any(d.get("size_estimated") for d in result["details"].values()):
            size_str = "~" + size_str

        chem_class = result.get("chemistry_classification")
        if chem_class and chem_class["pore"] != "unknown":
            chem_str = chem_class["pore"]
            rate = result.get("chemistry", {}).get("sample_rate", 0)
            if rate:
                chem_str += f" {rate // 1000}kHz"
        else:
            chem_str = "-"

        # For unknown/permission_denied runs, show a brief reason even in non-verbose mode
        if "unknown" in result["formats"] or "permission_denied" in result["formats"]:
            fmt_key = "permission_denied" if "permission_denied" in result["formats"] else "unknown"
            reasons = result["details"].get(fmt_key, {}).get("reasons", [])
            short_reason = reasons[0] if reasons else "?"
            print(f"{run_dir.name:<55} {fmt_key:<25} {'-':<18} {size_str:<12} {short_reason}")
        else:
            print(f"{run_dir.name:<55} {formats_str:<25} {chem_str:<18} {size_str:<12} {count_str}")

        if args.verbose:
            for fmt in result["formats"]:
                detail = result["details"].get(fmt, {})
                if fmt in ("unknown", "permission_denied"):
                    for reason in detail.get("reasons", []):
                        print(f"  ⚠  {reason}")
                    if detail.get("inaccessible_dirs"):
                        for d in detail["inaccessible_dirs"][:10]:
                            print(f"  🔒 {d}")
                    if detail.get("subdirectory_names"):
                        print(f"  └── subdirs: {', '.join(detail['subdirectory_names'])}")
                    if detail.get("file_extensions"):
                        ext_str = ", ".join(f"{c}x {e}" for e, c in detail["file_extensions"].items())
                        print(f"  └── files: {ext_str}")
                else:
                    if "directories" in detail:
                        for d in detail["directories"]:
                            print(f"  +-- {d}")
                    if "archive_files" in detail:
                        for a in detail["archive_files"]:
                            print(f"  +-- {a}")
                    if "data_size_bytes" in detail:
                        prefix = "~" if detail.get("size_estimated") else ""
                        print(f"  data size: {prefix}{format_size(detail['data_size_bytes'])}")
                    if "note" in detail:
                        print(f"  note: {detail['note']}")
                    print_conversion_help(fmt)
            # Show chemistry info in verbose mode
            chem = result.get("chemistry")
            chem_class = result.get("chemistry_classification")
            if chem and chem_class:
                rate = chem.get("sample_rate", 0)
                rate_str = f", {rate} Hz" if rate else ""
                print(f"  chemistry: {chem_class['pore']} ({chem['flowcell']}, {chem['kit']}{rate_str})")
                for line in dorado_advice_lines(chem_class, result["formats"]):
                    print(f"  {line}")
            print()

    # Process compressed archives in target folder as single_read_fast5
    for archive in archive_files:
        # Strip archive extensions to get the run name
        name = archive.name
        for ext in RAW_ARCHIVE_SUFFIXES + (".gz",):
            if name.lower().endswith(ext):
                name = name[: -len(ext)]
                break
        if name in all_runs:
            # A run folder of the same name exists: keep both entries
            name = archive.name

        result = {
            "formats": ["single_read_fast5"],
            "details": {
                "single_read_fast5": {
                    "directories": [str(archive.parent)],
                    "file_count": 0,
                    "archive_files": [archive.name],
                    "archive_count": 1,
                    "note": "Compressed archive — assumed single-read fast5",
                }
            },
        }
        # Compute archive file size (unless --quick)
        if not args.quick:
            try:
                archive_size = archive.stat().st_size
                result["total_size_bytes"] = archive_size
            except OSError:
                pass

        all_runs[name] = result
        format_summary["single_read_fast5"] = format_summary.get("single_read_fast5", 0) + 1
        arc_size_str = format_size(result["total_size_bytes"]) if "total_size_bytes" in result else "-"
        print(f"{archive.name:<55} {'single_read_fast5 (archive)':<25} {'-':<18} {arc_size_str:<12} -")

        if args.verbose:
            print(f"  ℹ  Compressed archive — assumed single-read fast5")
            print(f"  Extraction: tar xzf '{archive}' or gunzip '{archive}'")
            print()

    # Summary
    print(f"\n{'='*120}")
    print("Summary:")
    for fmt, count in format_summary.items():
        if count > 0:
            print(f"  {fmt:<25} {count} run(s)")

    # Generate conversion script if requested
    if args.convert_to:
        script = generate_conversion_script(all_runs, args.convert_to, args.output_dir)
        script_path = Path(args.script_output)
        script_path.write_text(script)
        script_path.chmod(0o755)
        print(f"\nConversion script written to: {script_path}")
        print(f"  Review and run: bash {script_path}")

    if args.output_stats:
        write_stats_tsv(all_runs, args.output_stats)
        print(f"\nStatistics written to: {args.output_stats}")


if __name__ == "__main__":
    main()
