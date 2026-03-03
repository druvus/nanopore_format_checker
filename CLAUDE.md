# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Single-file Python CLI tool (`nanopore_format_checker.py`) that scans a folder of Oxford Nanopore sequencing run directories and classifies each run by its read format: `pod5`, `multi_read_fast5`, or `single_read_fast5`. It also detects `fastq` output directories and compressed archives. Detects flowcell chemistry (R9.4.1, R10.3, R10.4, R10.4.1, RNA004) from data file metadata and recommends the appropriate dorado version for basecalling. Optionally generates bash conversion scripts with post-conversion metadata patching.

## Running

```bash
python nanopore_format_checker.py /path/to/runs_folder
python nanopore_format_checker.py /path/to/runs_folder --verbose
python nanopore_format_checker.py /path/to/runs_folder --quick
python nanopore_format_checker.py /path/to/runs_folder --convert-to pod5 --script-output convert.sh
python nanopore_format_checker.py /path/to/runs_folder --convert-to pod5 --output-dir /scratch/converted
python nanopore_format_checker.py /path/to/runs_folder -o stats.tsv
```

### Tests

```bash
python test_optimizations.py
```

### CLI flags

| Flag | Short | Description |
|------|-------|-------------|
| `--verbose` | `-v` | Show detailed per-run info and conversion hints |
| `--quick` | | Skip file counting and size calculation |
| `--convert-to {pod5,single_fast5}` | | Generate a bash conversion script |
| `--script-output FILE` | | Path for generated script (default: `convert_runs.sh`) |
| `--output-dir DIR` | | Base output directory for converted files (default: alongside originals) |
| `--output-stats FILE` | `-o` | Write per-run statistics to a TSV file |

## Dependencies

```bash
pip install h5py ont_fast5_api pod5
```

`h5py` is the only hard dependency (for fast5 inspection). `pod5` is optional (needed for pod5 chemistry extraction and metadata patching).

## Architecture

The tool is a single script with no module structure. Key functions:

- `is_nanopore_run_dir()` - Identifies run directories by the `YYYYMMDD_` naming convention
- `discover_run_structure()` - Single-pass directory walk that categorizes all format-related subdirectories (pod5, fast5, fastq variants)
- `analyze_run()` - Core function that classifies a single run directory. Detection priority: pod5 dirs -> fast5 dirs -> fast5 in root/numeric subdirs -> compressed archives -> unknown
- `classify_fast5()` / `classify_fast5_by_size()` - Distinguishes single vs multi-read fast5 using a 1 MB file size threshold (avoids slow HDF5 parsing)
- `fast_count_files()` / `estimate_dir_size()` - Uses `os.scandir` for performance with large directories (100k+ single-read fast5 files); sampling-based size estimation for very large directories
- `diagnose_unknown()` - Gathers diagnostic info when a run matches no known format
- `generate_conversion_script()` - Emits bash scripts using `pod5 convert` or `ont_fast5_api` tools; includes post-conversion metadata patching when chemistry was detected from source files
- `write_stats_tsv()` - Writes per-run statistics to a TSV file (one row per format per run)
- `extract_chemistry()` / `extract_chemistry_fast5()` / `extract_chemistry_pod5()` - Multi-fallback chemistry extraction from data files (one file per run)
- `classify_chemistry()` - Maps flowcell/kit/sample-rate to pore type and recommends dorado version

### Chemistry detection strategy

Chemistry extraction uses a multi-level fallback chain:

1. `context_tags/flowcell_type` and `sequencing_kit` (standard MinKNOW metadata)
2. `tracking_id/flow_cell_product_code` (when context_tags is empty)
3. `context_tags/experiment_kit` (older MinKNOW attribute name)
4. `channel_id/sampling_rate` (very old fast5 files with no context_tags or tracking_id)

For pod5 files, the same fallback applies to RunInfo fields, then `context_tags` and `tracking_id` dicts.

Pore type inference: `FLOWCELL_PORE` table (23 entries) -> `KIT_PORE` table (35 entries) -> sample rate (exact 5000 Hz = R10.4.1, all others = R9.4.1). The sample-rate fallback only triggers when both flowcell and kit are empty.

### Detection strategy

The tool searches for format-specific subdirectories (`pod5/`, `pod5_pass/`, `fast5/`, `fast5_pass/`, `fastq_pass/`, etc.) up to 5 levels deep via `discover_run_structure()`. For fast5, it samples one file and classifies by size rather than opening the HDF5 structure. Barcoded layouts (e.g., `fast5_pass/barcode13/`) are handled by descending into subdirectories during sampling.

## Conventions

- Permission errors are handled gracefully at every level; inaccessible directories are reported rather than causing failures
- File counting uses `os.scandir` instead of `Path.rglob` for performance on large directories
- The 1 MB threshold for single vs multi-read fast5 classification is a heuristic (single-read: typically 1-50 KB; multi-read: typically 1+ MB)
- `pod5 convert fast5` does not support single-read fast5 directly; conversion requires a two-step workflow: `single_to_multi_fast5` first, then `pod5 convert fast5`
- Conversion scripts include pod5 metadata patching (via the pod5 Python API) when chemistry was detected from source fast5 files -- this corrects empty or incorrect metadata in converted pod5 files
- Chemistry extraction is best-effort; failures never block format detection
- Flowcell product codes are normalized to uppercase (fast5 stores them as lowercase bytes)
- R9.4.1/R10.3 support was dropped in dorado 1.0; the tool recommends dorado 0.9.6 for these pore types
- `FLOWCELL_PORE` dict maps flowcell product codes to pore generations; update when new flowcells are released
- `KIT_PORE` dict maps sequencing kit codes to pore types; used as fallback when flowcell code is absent
