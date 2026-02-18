# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Single-file Python CLI tool (`nanopore_format_checker.py`) that scans a folder of Oxford Nanopore sequencing run directories and classifies each run by its read format: `pod5`, `multi_read_fast5`, or `single_read_fast5`. It can also detect `fastq` output directories and compressed archives. Detects flowcell chemistry (R9.4.1, R10.4.1) from data file metadata and recommends the appropriate dorado version for basecalling. Optionally generates bash conversion scripts to convert between formats.

## Running

```bash
python nanopore_format_checker.py /path/to/runs_folder
python nanopore_format_checker.py /path/to/runs_folder --verbose
python nanopore_format_checker.py /path/to/runs_folder --quick
python nanopore_format_checker.py /path/to/runs_folder --convert-to pod5 --script-output convert.sh
python nanopore_format_checker.py /path/to/runs_folder --convert-to pod5 --output-dir /scratch/converted
python nanopore_format_checker.py /path/to/runs_folder -o stats.tsv
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

`h5py` is the only hard dependency (for fast5 inspection). `pod5` is optional.

## Architecture

The tool is a single script with no module structure. Key functions:

- `is_nanopore_run_dir()` - Identifies run directories by the `YYYYMMDD_` naming convention
- `analyze_run()` - Core function that classifies a single run directory. Detection priority: pod5 dirs -> fast5 dirs -> fast5 in root/numeric subdirs -> compressed archives -> unknown
- `classify_fast5()` / `classify_fast5_by_size()` - Distinguishes single vs multi-read fast5 using a 1 MB file size threshold (avoids slow HDF5 parsing)
- `fast_count_files()` - Uses `os.scandir` for performance with large directories (100k+ single-read fast5 files)
- `diagnose_unknown()` - Gathers diagnostic info when a run matches no known format
- `generate_conversion_script()` - Emits bash scripts using `pod5 convert` or `ont_fast5_api` tools; supports `--output-dir` for writing to a separate location
- `write_stats_tsv()` - Writes per-run statistics to a TSV file (one row per format per run)
- `extract_chemistry()` / `extract_chemistry_fast5()` / `extract_chemistry_pod5()` - Reads flowcell product code, sequencing kit, and sample rate from data files (one file per run)
- `classify_chemistry()` - Maps flowcell code to pore type (R9.4.1/R10.4.1) and recommends dorado version

### Detection strategy

The tool searches for format-specific subdirectories (`pod5/`, `pod5_pass/`, `fast5/`, `fast5_pass/`, `fastq_pass/`, etc.) up to 5 levels deep. For fast5, it samples one file and classifies by size rather than opening the HDF5 structure. Single-read fast5 directories skip file counting due to the large number of files.

## Conventions

- Permission errors are handled gracefully at every level; inaccessible directories are reported rather than causing failures
- File counting uses `os.scandir` instead of `Path.rglob` for performance on large directories
- The 1 MB threshold for single vs multi-read fast5 classification is a heuristic (single-read: typically 1-50 KB; multi-read: typically 1+ MB)
- `pod5 convert fast5` does not support single-read fast5 directly; conversion requires a two-step workflow: `single_to_multi_fast5` first, then `pod5 convert fast5`
- Chemistry extraction is best-effort; failures never block format detection
- Flowcell product codes are normalized to uppercase (fast5 stores them as lowercase bytes)
- R9.4.1 support was dropped in dorado 1.0; the tool recommends dorado 0.9.6 for R9.4.1 data
- `FLOWCELL_PORE` dict maps flowcell product codes to pore generations; update when new flowcells are released
