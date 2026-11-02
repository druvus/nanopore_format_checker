# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Two standalone Python CLI tools for Oxford Nanopore sequencing data:

- **`nanopore_format_checker.py`** -- Scans a folder of nanopore run directories and classifies each run by its read format (`pod5`, `multi_read_fast5`, `single_read_fast5`). Detects `fastq` output directories and compressed archives. Extracts flowcell chemistry (R9.4.1, R9.5.1, R10.3, R10.4, R10.4.1, RNA002, RNA004) from data file metadata and recommends the appropriate dorado version and model for basecalling. Optionally generates bash conversion scripts with post-conversion metadata patching.

- **`nanopore_converter.py`** -- Converts fast5 files to pod5 format. Accepts a folder, compressed archive (`.tar.gz`, `.tgz`, `.tar.bz2`, `.tar`, `.zip`), or a text file listing paths (one per line, mixed folders and archives). Handles the two-step single_read_fast5 conversion automatically (`single_to_multi_fast5` then `pod5 convert fast5`). Calls conversion tools as subprocesses.

## Running

### Format checker

```bash
python nanopore_format_checker.py /path/to/runs_folder
python nanopore_format_checker.py /path/to/runs_folder --verbose
python nanopore_format_checker.py /path/to/runs_folder --quick
python nanopore_format_checker.py /path/to/runs_folder --convert-to pod5 --script-output convert.sh
python nanopore_format_checker.py /path/to/runs_folder --convert-to pod5 --output-dir /scratch/converted
python nanopore_format_checker.py /path/to/runs_folder -o stats.tsv
```

### Converter

```bash
python nanopore_converter.py /path/to/run_folder --output-dir /scratch/converted
python nanopore_converter.py /path/to/archive.tar.gz --output-dir /scratch/converted
python nanopore_converter.py input_list.txt --output-dir /scratch/converted --threads 8 --verbose
```

### Tests

```bash
python test_optimizations.py          # format checker tests (81 tests)
python -m pytest test_nanopore_converter.py -v  # converter tests (36 tests)
```

### Format checker CLI flags

| Flag | Short | Description |
|------|-------|-------------|
| `--verbose` | `-v` | Show detailed per-run info and conversion hints |
| `--quick` | | Skip file counting and size calculation |
| `--convert-to {pod5,single_fast5}` | | Generate a bash conversion script |
| `--script-output FILE` | | Path for generated script (default: `convert_runs.sh`) |
| `--output-dir DIR` | | Base output directory for converted files (default: alongside originals) |
| `--output-stats FILE` | `-o` | Write per-run statistics to a TSV file |

### Converter CLI flags

| Flag | Short | Description |
|------|-------|-------------|
| `input` | | Positional: folder, compressed archive, or text file listing paths |
| `--output-dir DIR` | | Required. Base output directory for converted pod5 files |
| `--threads N` | | Thread count for conversion tools (default: 4) |
| `--force` | | Re-convert inputs that already have output (default: skip them) |
| `--verbose` | `-v` | Enable detailed logging |

## Dependencies

```bash
pip install h5py ont_fast5_api pod5
```

`h5py` is the only hard dependency for the format checker (fast5 inspection). `pod5` is optional for the checker (needed for pod5 chemistry extraction and metadata patching). The converter requires `ont_fast5_api` (provides `single_to_multi_fast5` CLI) and `pod5` (provides `pod5 convert fast5` CLI) to be installed and on PATH.

## Architecture

### Format checker (`nanopore_format_checker.py`)

Single script with no module structure. Key functions:

- `is_nanopore_run_dir()` - Identifies run directories by the `YYYYMMDD_` naming convention
- `discover_run_structure()` - Single-pass directory walk that categorizes all format-related subdirectories (pod5, fast5, fastq variants)
- `analyze_run()` - Core function that classifies a single run directory. Detection priority: pod5 dirs -> fast5 dirs -> fast5 in root/numeric subdirs -> compressed archives -> unknown
- `classify_fast5()` / `classify_fast5_by_structure()` / `classify_fast5_by_size()` - Distinguishes single vs multi-read fast5 by reading the HDF5 root group of one file (`read_*` groups = multi-read); falls back to a 1 MB size threshold when h5py is missing or the file is unparseable
- `is_raw_archive()` - Decides whether a file is a raw-data archive (`.tar`, `.tar.gz`, `.tgz`, `.gz`), excluding gzipped text/sequence outputs like `sequencing_summary.txt.gz`
- `fast_count_files()` / `estimate_dir_size()` - Uses `os.scandir` for performance with large directories (100k+ single-read fast5 files); sampling-based size estimation for very large directories
- `diagnose_unknown()` - Gathers diagnostic info when a run matches no known format
- `generate_conversion_script()` - Emits bash scripts using `pod5 convert` or `ont_fast5_api` tools; includes post-conversion metadata patching when chemistry was detected from source files
- `write_stats_tsv()` - Writes per-run statistics to a TSV file (one row per format per run), including `dorado_version`, `dorado_model` and `dorado_note`
- `extract_chemistry()` / `extract_chemistry_fast5()` / `extract_chemistry_pod5()` - Multi-fallback chemistry extraction from data files (one file per run)
- `classify_chemistry()` - Maps flowcell/kit/sample-rate to pore type, dorado version (`None` when no dorado release supports it), model hint and explicit legacy model name (`DORADO_LEGACY_MODELS`)
- `dorado_advice_lines()` - Verbose per-run dorado advice (version, fast5 vs pod5 input, auto-selection vs explicit model)

### Chemistry detection strategy

Chemistry extraction uses a multi-level fallback chain:

1. `context_tags/flowcell_type` and `sequencing_kit` (standard MinKNOW metadata)
2. `tracking_id/flow_cell_product_code` (when context_tags is empty)
3. `context_tags/experiment_kit` (older MinKNOW attribute name)
4. `channel_id/sampling_rate` (very old fast5 files with no context_tags or tracking_id)

For pod5 files, the same fallback applies to RunInfo fields, then `context_tags` and `tracking_id` dicts.

Pore type inference: `FLOWCELL_PORE` table (21 entries) -> `KIT_PORE` table (55 entries) -> sample rate (exact 5000 Hz = R10.4.1, all others = R9.4.1). The sample-rate fallback only triggers when both flowcell and kit are empty.

### Detection strategy

The tool searches for format-specific subdirectories (`pod5/`, `pod5_pass/`, `fast5/`, `fast5_pass/`, `fastq_pass/`, etc.) up to 5 levels deep via `discover_run_structure()`. For fast5, it samples one file and classifies it by its HDF5 root layout (size fallback). Archive detection only runs when no pod5/fast5 data was found. Barcoded layouts (e.g., `fast5_pass/barcode13/`) are handled by descending into subdirectories during sampling.

### Converter (`nanopore_converter.py`)

Single script that orchestrates format conversion via subprocess calls. Key functions:

- `parse_args()` - CLI argument parsing (input, --output-dir, --threads, --force, --verbose)
- `check_tools()` - Verifies `single_to_multi_fast5` and `pod5` are on PATH; exits with code 2 if missing
- `detect_input_type()` - Auto-detects whether input is a directory, archive (`.tar.gz`, `.tgz`, `.tar.bz2`, `.tar`, `.zip`), or list file (`.txt`, `.csv`, `.list`)
- `extract_archive()` - Extracts compressed archives to a temp directory (tar uses the `data` extraction filter to block path traversal); caller handles cleanup
- `parse_list_file()` - Reads a text file of paths (one per line, skips blanks and `#` comments)
- `resolve_inputs()` - Resolves the input argument to a flat list of paths to process
- `classify_fast5_dir()` - Walks a directory, samples one `.fast5` file, classifies by HDF5 layout with the 1 MB size threshold as fallback (same approach as the format checker)
- `convert_run()` - Core conversion: for single_read_fast5, runs `single_to_multi_fast5` then `pod5 convert fast5`; for multi_read_fast5, runs `pod5 convert fast5` directly. Raises if `single_to_multi_fast5` writes no output (it exits 0 even when every file fails). Writes into `<run>/.pod5_partial/` and swaps it into `<run>/pod5/` only on success (old output is renamed aside, then deleted); clears leftovers from interrupted runs; always removes the intermediate `_multi_fast5_tmp/`
- `fill_missing_pod5_metadata()` - After conversion, fills empty `flow_cell_product_code` / `sequencing_kit` / `sample_rate` in pod5 RunInfo from the source fast5 (via the checker's `extract_chemistry_fast5`); never overwrites, never fails the conversion
- `get_output_name()` - Derives output folder name from the resolved input path, stripping archive extensions
- `is_converted()` - True if `<output-dir>/<name>/pod5/` already holds `.pod5` files; `main()` skips such inputs unless `--force`, and fails duplicate output names within one invocation
- `process_item()` - Orchestrates one input item (directory or archive) with try/except/finally cleanup of temp directories
- `main()` - Entry point: parses args, checks tools, resolves inputs, iterates over items, tracks success/failure counts

### Conversion pipeline

```
single_read_fast5 -> single_to_multi_fast5 -> pod5 convert fast5 -> cleanup temp
multi_read_fast5  -> pod5 convert fast5
```

Output structure: `<output-dir>/<input_basename>/pod5/*.pod5`. A `pod5/` folder only exists after a fully successful conversion, which is what makes skip-on-rerun safe.

### Exit codes (converter)

- 0: all conversions succeeded
- 1: some conversions failed (summary printed)
- 2: fatal error (bad arguments, missing tools)

## Conventions

- Permission errors are handled gracefully at every level; inaccessible directories are reported rather than causing failures
- File counting uses `os.scandir` instead of `Path.rglob` for performance on large directories
- The 1 MB threshold for single vs multi-read fast5 classification is only a fallback heuristic (single-read: typically 1-50 KB; multi-read: typically 1+ MB, but small multi-read files exist, e.g. in `fast5_fail/`)
- `pod5.RunInfo` is a frozen dataclass: metadata patching must use `dataclasses.replace()`, not attribute assignment
- Generated scripts single-quote every path via `_sh_quote()`; archive-only runs are emitted as `# SKIPPED` comments
- `pod5 convert fast5` does not support single-read fast5 directly; conversion requires a two-step workflow: `single_to_multi_fast5` first, then `pod5 convert fast5`
- Conversion scripts include pod5 metadata patching (via the pod5 Python API) when chemistry was detected from source fast5 files -- this corrects empty or incorrect metadata in converted pod5 files
- Chemistry extraction is best-effort; failures never block format detection
- Flowcell product codes are normalized to uppercase (fast5 stores them as lowercase bytes)
- Dorado support matrix (verified against the dorado changelog and `dorado/models/kits.cpp`): >=1.0 supports R10.4.1 5 kHz and RNA004 only, with POD5 input only (FAST5 removed in 1.0.0); 0.9.6 is the last release with R9.4.1 (4 kHz), R10.4.1 4 kHz (400 and 260 bps) and RNA002; no dorado release supports R10.3, R10.4 (Kit 12) or DNA on R9.5.1 (FLO-MIN107) -- those get `dorado_version = None` and a Guppy note
- dorado's automatic model selection looks up (flowcell, kit, sample_rate) from pod5 RunInfo (case-insensitive) in its chemistry table; that is why conversion patches/fills these fields. Never invent a flowcell code that the source data does not contain
- RNA002 has no sup model (`model_hint` is `hac`)
- `FLOWCELL_PORE` dict maps flowcell product codes to pore generations; update when new flowcells are released
- `KIT_PORE` dict maps sequencing kit codes to pore types; used as fallback when flowcell code is absent
