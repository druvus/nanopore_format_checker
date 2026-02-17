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

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False


def is_nanopore_run_dir(dirname: str) -> bool:
    """Check if directory name matches nanopore run naming convention (starts with date)."""
    return bool(re.match(r"^\d{8}_", dirname))


def find_subdirs(run_path: Path, target_name: str, max_depth: int = 2) -> list[Path]:
    """Find subdirectories matching target_name up to max_depth levels down."""
    matches = []
    for depth in range(1, max_depth + 1):
        pattern = "/".join(["*"] * depth)
        for p in run_path.glob(pattern):
            if p.is_dir() and p.name == target_name:
                matches.append(p)
    return matches


def find_subdirs_prefix(run_path: Path, prefix: str, max_depth: int = 3) -> list[Path]:
    """Find subdirectories whose name starts with prefix, up to max_depth levels down."""
    matches = []
    for depth in range(1, max_depth + 1):
        pattern = "/".join(["*"] * depth)
        for p in run_path.glob(pattern):
            if p.is_dir() and p.name.startswith(prefix):
                matches.append(p)
    return matches


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


def fast_count_files(directory: Path, ext: str, recursive: bool = False) -> int:
    """
    Count files with given extension using os.scandir for speed.
    For large directories (100k+ single-read fast5), this is much faster
    than Path.rglob or Path.iterdir with stat calls.
    Supports compound extensions like .fastq.gz.
    """
    count = 0
    try:
        with os.scandir(directory) as it:
            for entry in it:
                if entry.is_file(follow_symlinks=False) and entry.name.endswith(ext):
                    count += 1
                elif recursive and entry.is_dir(follow_symlinks=False):
                    count += fast_count_files(Path(entry.path), ext, recursive=True)
    except PermissionError:
        pass
    return count


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


def classify_fast5(fast5_file: Path) -> str:
    """
    Determine if a fast5 file is single-read or multi-read.
    
    Uses file size heuristic for speed - single-read files are typically <1MB,
    multi-read files are typically >1MB.
    """
    return classify_fast5_by_size(fast5_file)


def count_files_recursive(directory: Path, ext: str) -> int:
    """Count files with given extension recursively using fast os.scandir."""
    return fast_count_files(directory, ext, recursive=True)


def analyze_run(run_path: Path) -> dict:
    """
    Analyze a nanopore run directory and determine its read format(s).

    Detection logic:
        - pod5:             pod5/ subfolder 1-2 levels down
        - multi_read_fast5: fast5/ subfolder 1-2 levels down, files contain multiple reads
        - single_read_fast5: .fast5 files in run root or in numeric subfolders (0/, 1/, 2/, ...)
                             each file contains a single read

    Returns dict with:
        - formats: list of detected formats
        - details: dict with per-format info (paths, file counts)
    """
    result = {"formats": [], "details": {}}

    # --- Check permissions on run directory and its subdirectories ---
    permission_errors = []
    try:
        contents = list(run_path.iterdir())
    except PermissionError:
        result["formats"].append("permission_denied")
        result["details"]["permission_denied"] = {
            "reasons": [f"Cannot read run directory: {run_path}"],
        }
        return result

    subdirs = [c for c in contents if c.is_dir()]
    for sd in subdirs:
        try:
            list(sd.iterdir())
        except PermissionError:
            permission_errors.append(sd.name)

    if permission_errors and not any(
        # Check if at least some subdirs are readable
        sd.name not in permission_errors for sd in subdirs if sd.is_dir()
    ):
        # ALL subdirectories are unreadable
        result["formats"].append("permission_denied")
        result["details"]["permission_denied"] = {
            "reasons": [f"All subdirectories are unreadable ({len(permission_errors)} dir(s))"],
            "inaccessible_dirs": permission_errors,
        }
        return result

    # --- Check for pod5 (pod5/ or pod5_pass/pod5_fail/pod5_skip, up to 3 levels down) ---
    pod5_dirs = find_subdirs(run_path, "pod5", max_depth=5)
    pod5_variant_dirs = find_subdirs_prefix(run_path, "pod5_", max_depth=5)
    all_pod5_dirs = pod5_dirs + pod5_variant_dirs
    if all_pod5_dirs:
        total_pod5 = sum(count_files_recursive(d, ".pod5") for d in all_pod5_dirs)
        variant_names = sorted(set(d.name for d in all_pod5_dirs))
        # Check if directories are readable
        unreadable_pod5 = []
        for d in all_pod5_dirs:
            try:
                list(d.iterdir())
            except PermissionError:
                unreadable_pod5.append(str(d))
        pod5_detail = {
            "directories": [str(d) for d in all_pod5_dirs],
            "folder_variants": variant_names,
            "file_count": total_pod5,
        }
        if unreadable_pod5:
            pod5_detail["note"] = f"Permission denied on {len(unreadable_pod5)} pod5 dir(s)"
            pod5_detail["inaccessible_dirs"] = unreadable_pod5
        elif total_pod5 == 0:
            pod5_detail["note"] = "Empty pod5 folder(s) found"
        result["formats"].append("pod5")
        result["details"]["pod5"] = pod5_detail

    # --- Check for fast5 (fast5/ or fast5_skip/fast5_pass/fast5_fail, up to 5 levels down) ---
    fast5_dirs = find_subdirs(run_path, "fast5", max_depth=5)
    fast5_variant_dirs = find_subdirs_prefix(run_path, "fast5_", max_depth=5)
    all_fast5_dirs = fast5_dirs + fast5_variant_dirs
    if all_fast5_dirs:
        # Check which directories are readable
        unreadable_fast5 = []
        for d in all_fast5_dirs:
            try:
                list(d.iterdir())
            except PermissionError:
                unreadable_fast5.append(str(d))

        # Sample a fast5 file using os.scandir — also peek into numeric subdirs (0/, 1/, ...)
        sample_file = None
        has_numeric_subdirs = False
        for d in all_fast5_dirs:
            if str(d) in unreadable_fast5:
                continue
            try:
                with os.scandir(d) as it:
                    for entry in it:
                        if sample_file is None and entry.is_file(follow_symlinks=False) and entry.name.endswith(".fast5"):
                            sample_file = Path(entry.path)
                            break
                        elif entry.is_dir(follow_symlinks=False) and re.match(r"^\d+$", entry.name):
                            has_numeric_subdirs = True
                            if sample_file is None:
                                try:
                                    with os.scandir(entry.path) as sub_it:
                                        for sub_entry in sub_it:
                                            if sub_entry.is_file(follow_symlinks=False) and sub_entry.name.endswith(".fast5"):
                                                sample_file = Path(sub_entry.path)
                                                break
                                except PermissionError:
                                    pass
                        if sample_file is not None:
                            break
            except PermissionError:
                pass
            if sample_file is not None:
                break

        variant_names = sorted(set(d.name for d in all_fast5_dirs))
        fast5_info = {
            "directories": [str(d) for d in all_fast5_dirs],
            "folder_variants": variant_names,
        }
        if unreadable_fast5:
            fast5_info["inaccessible_dirs"] = unreadable_fast5

        if sample_file:
            classification = classify_fast5(sample_file)
            fast5_info["sampled_file"] = sample_file.name

            if classification == "single_read_fast5":
                # Single-read fast5 — skip file counting (too many files)
                if has_numeric_subdirs:
                    fast5_info["note"] = "fast5/ contains numeric subdirs with single-read files"
                else:
                    fast5_info["note"] = "fast5/ contains single-read files"
                result["formats"].append("single_read_fast5")
                result["details"]["single_read_fast5"] = fast5_info
            elif classification == "multi_read_fast5":
                # Multi-read — count is fast (few large files)
                total_fast5 = sum(count_files_recursive(d, ".fast5") for d in all_fast5_dirs)
                fast5_info["file_count"] = total_fast5
                result["formats"].append("multi_read_fast5")
                result["details"]["multi_read_fast5"] = fast5_info
            else:
                result["formats"].append("fast5_unknown")
                result["details"]["fast5_unknown"] = fast5_info
                if not HAS_H5PY:
                    fast5_info["note"] = "Install h5py to distinguish single/multi read fast5"
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

    # --- Check for fastq (fastq_pass/fastq_fail, up to 3 levels down) ---
    fastq_dirs = find_subdirs_prefix(run_path, "fastq", max_depth=5)
    if fastq_dirs:
        total_fastq = 0
        unreadable_fastq = []
        for d in fastq_dirs:
            try:
                list(d.iterdir())
            except PermissionError:
                unreadable_fastq.append(str(d))
            total_fastq += count_files_recursive(d, ".fastq")
            total_fastq += count_files_recursive(d, ".fastq.gz")
            total_fastq += count_files_recursive(d, ".fq")
            total_fastq += count_files_recursive(d, ".fq.gz")
        variant_names = sorted(set(d.name for d in fastq_dirs))
        fastq_detail = {
            "directories": [str(d) for d in fastq_dirs],
            "folder_variants": variant_names,
            "file_count": total_fastq,
        }
        if unreadable_fastq:
            fastq_detail["note"] = f"Permission denied on {len(unreadable_fastq)} fastq dir(s)"
            fastq_detail["inaccessible_dirs"] = unreadable_fastq
        elif total_fastq == 0:
            fastq_detail["note"] = "Empty fastq folder(s) found"
        result["formats"].append("fastq")
        result["details"]["fastq"] = fastq_detail

    # --- Check for single_read_fast5 (.fast5 in run root or numeric subdirs 0/, 1/, 2/, ...) ---
    # Only check if we didn't already find fast5 via a fast5/ subfolder
    # Uses os.scandir for speed — single_fast5 dirs can have 100k+ files
    if not all_fast5_dirs:
        sample_file = None
        classification = None
        has_numeric_dirs = False
        fast5_in_root = False

        try:
            with os.scandir(run_path) as it:
                for entry in it:
                    # Check for .fast5 file in root (just need one for classification)
                    if not sample_file and entry.is_file(follow_symlinks=False) and entry.name.endswith(".fast5"):
                        # Quick size check without creating Path object
                        size_bytes = entry.stat(follow_symlinks=False).st_size
                        classification = "single_read_fast5" if size_bytes < 1_000_000 else "multi_read_fast5"
                        sample_file = entry.name
                        fast5_in_root = True
                    # Check for numeric subdirectory
                    elif not has_numeric_dirs and entry.is_dir(follow_symlinks=False) and re.match(r"^\d+$", entry.name):
                        has_numeric_dirs = True
                        # Try to find and classify a .fast5 inside this numeric dir
                        if not sample_file:
                            try:
                                with os.scandir(entry.path) as sub_it:
                                    for sub_entry in sub_it:
                                        if sub_entry.is_file(follow_symlinks=False) and sub_entry.name.endswith(".fast5"):
                                            size_bytes = sub_entry.stat(follow_symlinks=False).st_size
                                            classification = "single_read_fast5" if size_bytes < 1_000_000 else "multi_read_fast5"
                                            sample_file = sub_entry.name
                                            break
                            except PermissionError:
                                pass
                    # Stop early once we have what we need
                    if sample_file:
                        break
        except PermissionError:
            pass

        if sample_file and classification:
            fast5_info = {
                "sampled_file": sample_file,
                "layout": "root" if fast5_in_root else "numeric_subdirs",
                "size_based_classification": True,
            }

            if classification in ("single_read_fast5", "multi_read_fast5"):
                result["formats"].append(classification)
                result["details"][classification] = fast5_info
            else:
                result["formats"].append("fast5_unknown")
                result["details"]["fast5_unknown"] = fast5_info

    # --- Check for compressed archives (tar, gz, tar.gz) → treat as single_read_fast5 ---
    if not fast5_dirs and "single_read_fast5" not in result["formats"]:
        archive_exts = (".tar", ".gz", ".tar.gz", ".tgz")
        archives = []
        try:
            for f in run_path.iterdir():
                if f.is_file() and any(f.name.endswith(ext) for ext in archive_exts):
                    archives.append(f)
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


def diagnose_unknown(run_path: Path) -> dict:
    """Gather diagnostic info for a run directory that matched no known format."""
    diag = {"reasons": []}

    try:
        contents = list(run_path.iterdir())
    except PermissionError:
        diag["reasons"].append("Permission denied reading directory")
        return diag

    if not contents:
        diag["reasons"].append("Directory is empty")
        return diag

    # Categorize contents
    subdirs = [c for c in contents if c.is_dir()]
    files = [c for c in contents if c.is_file()]

    diag["subdirectory_names"] = [d.name for d in subdirs]
    diag["file_count"] = len(files)

    # Check file extensions present
    extensions = {}
    for f in files:
        ext = f.suffix.lower() if f.suffix else "(no extension)"
        extensions[ext] = extensions.get(ext, 0) + 1
    diag["file_extensions"] = extensions

    # Build human-readable reasons
    if not subdirs and not files:
        diag["reasons"].append("Directory is empty")
    elif not subdirs and files:
        if extensions:
            ext_summary = ", ".join(f"{count}x {ext}" for ext, count in sorted(extensions.items()))
            diag["reasons"].append(f"No subdirectories; {len(files)} file(s) in root: {ext_summary}")
        else:
            diag["reasons"].append(f"No subdirectories; {len(files)} file(s) with no recognized extensions")

        # Specific hints
        if not any(ext in extensions for ext in (".fast5", ".pod5")):
            diag["reasons"].append("No .fast5 or .pod5 files found")
    else:
        dir_names = ", ".join(d.name for d in subdirs[:10])
        diag["reasons"].append(f"Subdirectories found: {dir_names}")

        # Check which subdirs are readable
        unreadable = []
        readable_subdirs = []
        for sd in subdirs:
            try:
                list(sd.iterdir())
                readable_subdirs.append(sd)
            except PermissionError:
                unreadable.append(sd.name)

        if unreadable:
            diag["inaccessible_dirs"] = unreadable
            diag["reasons"].append(
                f"Permission denied on {len(unreadable)} subdir(s): {', '.join(unreadable[:5])}"
                + (f" (+{len(unreadable)-5} more)" if len(unreadable) > 5 else "")
            )

        if not readable_subdirs and unreadable:
            diag["reasons"].append("All subdirectories are unreadable — cannot determine format")
            return diag

        # Check if readable subdirs contain fast5/pod5 deeper than we searched
        has_deep_fast5 = False
        has_deep_pod5 = False
        for sd in readable_subdirs[:10]:
            try:
                for f in sd.rglob("*.fast5"):
                    has_deep_fast5 = True
                    break
                for f in sd.rglob("*.pod5"):
                    has_deep_pod5 = True
                    break
            except PermissionError:
                continue
            if has_deep_fast5 or has_deep_pod5:
                break

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
        print("    → To pod5:          pod5 convert fast5 input.fast5 --output output.pod5")
        print("                        pod5 convert fast5 ./fast5_dir/ --output ./pod5_dir/")
        print("                        # https://github.com/nanoporetech/pod5-file-format")
        print("    → To single_fast5:  multi_to_single_fast5 --input_path multi/ --save_path single/")
        print("                        # https://github.com/nanoporetech/ont_fast5_api")
    elif fmt == "single_read_fast5":
        print("  Conversion options:")
        print("    → To pod5:          pod5 convert fast5 ./fast5_dir/ --output ./pod5_dir/")
        print("    → To multi_fast5:   single_to_multi_fast5 --input_path single/ --save_path multi/")
    elif fmt == "pod5":
        print("  Conversion options:")
        print("    → To fast5:         pod5 convert to_fast5 input.pod5 --output ./fast5_dir/")
        print("    → Inspect:          pod5 inspect reads input.pod5")
        print("    → Merge:            pod5 merge output.pod5 input1.pod5 input2.pod5")


def generate_conversion_script(runs: dict, target_format: str, output_dir: str = "converted"):
    """Generate a bash conversion script."""
    lines = [
        "#!/usr/bin/env bash",
        "# Auto-generated nanopore format conversion script",
        f"# Target format: {target_format}",
        "set -euo pipefail",
        "",
    ]

    for run_name, info in runs.items():
        for fmt in info["formats"]:
            if fmt == target_format:
                continue  # Already in target format

            if target_format == "pod5" and fmt in ("multi_read_fast5", "single_read_fast5"):
                dirs = info["details"].get(fmt, {}).get("directories", [])
                for d in dirs:
                    out = os.path.join(os.path.dirname(d), "pod5")
                    lines.append(f"echo 'Converting {run_name} ({fmt} → pod5)...'")
                    lines.append(f"mkdir -p '{out}'")
                    lines.append(f"pod5 convert fast5 '{d}/' --output '{out}/' --threads 4")
                    lines.append("")

            elif target_format == "single_fast5" and fmt == "multi_read_fast5":
                dirs = info["details"].get(fmt, {}).get("directories", [])
                for d in dirs:
                    out = os.path.join(os.path.dirname(d), "single_fast5")
                    lines.append(f"echo 'Converting {run_name} ({fmt} → single_fast5)...'")
                    lines.append(f"mkdir -p '{out}'")
                    lines.append(f"multi_to_single_fast5 --input_path '{d}' --save_path '{out}'")
                    lines.append("")

    return "\n".join(lines)


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
        "--verbose", "-v",
        action="store_true",
        help="Show detailed information per run",
    )
    args = parser.parse_args()

    target = Path(args.target_folder)
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
    archive_exts = (".tar", ".gz", ".tar.gz", ".tgz")
    archive_files = sorted(
        [f for f in target.iterdir()
         if f.is_file() and any(f.name.endswith(ext) for ext in archive_exts)
         and is_nanopore_run_dir(f.name)]
    )

    if not run_dirs and not archive_files:
        print(f"No nanopore run directories found in '{target}'")
        print("  (Expected directories starting with YYYYMMDD_ pattern)")
        sys.exit(0)

    print(f"Found {len(run_dirs)} nanopore run(s) and {len(archive_files)} archive(s) in '{target}'\n")
    print(f"{'Run / Archive':<55} {'Format(s)':<30} {'Files'}")
    print("-" * 100)

    all_runs = {}
    format_summary = {"pod5": 0, "multi_read_fast5": 0, "single_read_fast5": 0,
                       "fastq": 0, "fast5_unknown": 0, "permission_denied": 0, "unknown": 0}

    for run_dir in run_dirs:
        result = analyze_run(run_dir)
        all_runs[run_dir.name] = result

        formats_str = ", ".join(result["formats"])
        file_counts = []
        for fmt in result["formats"]:
            detail = result["details"].get(fmt, {})
            count = detail.get("file_count", 0)
            if count:
                file_counts.append(f"{count}")
            format_summary[fmt] = format_summary.get(fmt, 0) + 1

        count_str = " / ".join(file_counts) if file_counts else "-"
        # For unknown/permission_denied runs, show a brief reason even in non-verbose mode
        if "unknown" in result["formats"] or "permission_denied" in result["formats"]:
            fmt_key = "permission_denied" if "permission_denied" in result["formats"] else "unknown"
            reasons = result["details"].get(fmt_key, {}).get("reasons", [])
            short_reason = reasons[0] if reasons else "?"
            print(f"{run_dir.name:<55} {fmt_key:<30} {short_reason}")
        else:
            print(f"{run_dir.name:<55} {formats_str:<30} {count_str}")

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
                            print(f"  └── {d}")
                    if "archive_files" in detail:
                        for a in detail["archive_files"]:
                            print(f"  └── {a}")
                    if "note" in detail:
                        print(f"  ℹ  {detail['note']}")
                    print_conversion_help(fmt)
            print()

    # Process compressed archives in target folder as single_read_fast5
    for archive in archive_files:
        # Strip archive extensions to get the run name
        name = archive.name
        for ext in (".tar.gz", ".tgz", ".tar", ".gz"):
            if name.endswith(ext):
                name = name[: -len(ext)]
                break

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
        all_runs[name] = result
        format_summary["single_read_fast5"] = format_summary.get("single_read_fast5", 0) + 1
        print(f"{archive.name:<55} {'single_read_fast5 (archive)':<30} -")

        if args.verbose:
            print(f"  ℹ  Compressed archive — assumed single-read fast5")
            print(f"  Extraction: tar xzf '{archive}' or gunzip '{archive}'")
            print()

    # Summary
    print(f"\n{'='*100}")
    print("Summary:")
    for fmt, count in format_summary.items():
        if count > 0:
            print(f"  {fmt:<25} {count} run(s)")

    # Generate conversion script if requested
    if args.convert_to:
        script = generate_conversion_script(all_runs, args.convert_to)
        script_path = Path(args.script_output)
        script_path.write_text(script)
        script_path.chmod(0o755)
        print(f"\nConversion script written to: {script_path}")
        print(f"  Review and run: bash {script_path}")


if __name__ == "__main__":
    main()
