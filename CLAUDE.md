# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Standalone Python CLI tools for Oxford Nanopore sequencing data:

- **`nanopore_format_checker.py`** -- Scans a folder of nanopore run directories and classifies each run by its read format (`pod5`, `multi_read_fast5`, `single_read_fast5`). Detects `fastq` output directories and compressed archives. Extracts flowcell chemistry (R9.4.1, R9.5.1, R10.3, R10.4, R10.4.1, RNA002, RNA004) from data file metadata and recommends the appropriate dorado version and model for basecalling. Optionally generates bash conversion scripts with post-conversion metadata patching.

- **`nanopore_converter.py`** -- Converts fast5 files to pod5 format. Accepts a folder, compressed archive (`.tar.gz`, `.tgz`, `.tar.bz2`, `.tar`, `.zip`), or a text file listing paths (one per line, mixed folders and archives). Handles the two-step single_read_fast5 conversion automatically (`single_to_multi_fast5` then `pod5 convert fast5`). Calls conversion tools as subprocesses.

- **`compress_folders.py`** -- Archives each subfolder of an input directory into its own `.tar.gz` (copies existing archives, recompresses plain `.tar`), in parallel, with atomic writes and skip-if-exists re-runs. `_add_tree()` walks the tree itself instead of `tarfile.add()` so an unreadable file or subfolder is skipped and counted rather than failing the whole archive; an unlistable or untraversable top-level folder is an error (an empty archive would be skipped forever on re-runs). A top-level symlink is archived as its target's content; items sharing a destination archive are errors; the output dir is never archived into itself; Ctrl-C cancels queued items (`_STOP`); a hard link to a skipped file is stored as a regular file.

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
python test_optimizations.py          # format checker tests (94 tests)
python -m pytest test_nanopore_converter.py -v  # converter tests (46 tests)
python -m pytest test_compress_folders.py -v    # archiver tests (12; two are skipped when run as root)
python -m pytest test_integration.py -v         # real pod5/single_to_multi_fast5 end-to-end (5; skipped without the tools)
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

Python 3.12 or newer is required; each script exits with an error on older versions (version guard right after the stdlib imports). CI tests 3.12, 3.13 and 3.14.

```bash
pip install h5py pod5 progressbar2 packaging "setuptools<81"
pip install --no-deps ont_fast5_api
```

A plain `pip install ont_fast5_api` fails on current Python (pinned `progressbar33` fails to build; it imports `pkg_resources`, gone from setuptools 81+); the two commands above are what CI runs.

`h5py` is the only hard dependency for the format checker (fast5 inspection). `pod5` is optional for the checker (needed for pod5 chemistry extraction and metadata patching). The converter requires `ont_fast5_api` (provides `single_to_multi_fast5` CLI) and `pod5` (provides `pod5 convert fast5` CLI) to be installed and on PATH.

## Architecture

### Format checker (`nanopore_format_checker.py`)

Single script with no module structure. Key functions:

- `is_nanopore_run_dir()` - Identifies run directories by the `YYYYMMDD_` naming convention
- `discover_run_structure()` - Single-pass directory walk that categorizes all format-related subdirectories (pod5, fast5, fastq variants)
- `analyze_run()` - Core function that classifies a single run directory. Detection priority: pod5 dirs -> fast5 dirs -> fast5 in root/numeric subdirs -> compressed archives -> unknown
- `find_fast5_samples()` / `classify_fast5_samples()` - Collect up to 5 sample files (files before subfolders, up to 4 levels, symlinks followed with a loop guard, entry budget) and classify from the first one h5py can parse; size fallback only if none parses
- `classify_fast5()` / `classify_fast5_by_structure()` / `classify_fast5_by_size()` - Single-file classification by HDF5 root group (`read_*` groups = multi-read), 1 MB size threshold as fallback
- `is_raw_archive()` - Decides whether a file is a raw-data archive (tar in any compression, `.zip`, `.gz`), excluding gzipped text/sequence outputs like `sequencing_summary.txt.gz`
- `fast_count_files()` / `estimate_dir_size()` - Uses `os.scandir` for performance with large directories (100k+ single-read fast5 files); sampling-based size estimation for very large directories
- `diagnose_unknown()` - Gathers diagnostic info when a run matches no known format
- `generate_conversion_script()` / `_pod5_conversion_lines()` / `_guarded_block()` - Emit re-runnable bash: per run, skip if output holds .pod5; otherwise work in `<base>/.convert_tmp`, convert only the detected fast5 folders, verify read counts (`_verify_read_count_lines`), patch metadata, move into place. Each run is a fail-fast subshell whose status is checked on its own line (`( ... ) ||` would disable `set -e` inside); failures are collected and the script exits 1 at the end. In-place output is always `<run>/pod5_from_fast5`
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

The tool searches for format-specific subdirectories (`pod5/`, `pod5_pass/`, `fast5/`, `fast5_pass/`, `fastq_pass/`, etc.) up to 5 levels deep via `discover_run_structure()`. For fast5, it samples up to five files and classifies by the first readable HDF5 root layout (size fallback). Symlinked folders are followed; hidden folders are skipped. Archive detection only runs when no pod5/fast5 data was found. Barcoded layouts (e.g., `fast5_pass/barcode13/`) are handled by descending into subdirectories during sampling.

### Converter (`nanopore_converter.py`)

Single script that orchestrates format conversion via subprocess calls. Key functions:

- `parse_args()` - CLI argument parsing (input, --output-dir, --threads, --force, --verbose)
- `check_tools()` - Verifies `single_to_multi_fast5` and `pod5` are on PATH; exits with code 2 if missing
- `detect_input_type()` - Auto-detects whether input is a directory, archive (`.tar.gz`, `.tgz`, `.tar.bz2`, `.tbz2`, `.tar.xz`, `.txz`, `.tar`, `.zip`), or list file (`.txt`, `.csv`, `.list`)
- `extract_archive()` - Extracts compressed archives to a temp directory (tar uses the `data` extraction filter to block path traversal); caller handles cleanup
- `parse_list_file()` - Reads a text file of paths (one per line, skips blanks and `#` comments)
- `resolve_inputs()` - Resolves the input argument to a flat list of paths to process
- `scan_fast5_dir()` - Opens every `.fast5` with h5py and groups it as single-read, multi-read (with read count) or unreadable; `expected_reads` drives verification. `classify_fast5_dir()` summarises it (`mixed_fast5` possible)
- `_run()` - Runs a tool in its own process group, logs problem lines from stderr at WARNING even on exit 0, kills the whole group on interrupt
- `convert_run()` - Core conversion: stages single-read files (symlinks) for `single_to_multi_fast5` and multi-read files for `pod5 convert fast5` (a run may contain both), then `_verify_read_count()` requires the pod5 output to hold exactly `expected_reads` (raises `ConversionError`). Unreadable files are skipped and reported (`SKIPPED_UNREADABLE`, repeated in the summary). Writes `pod5/.source`. Writes into `<run>/.pod5_partial/` and swaps it into `<run>/pod5/` only on success (old output is renamed aside, then deleted); clears leftovers from interrupted runs; always removes the intermediate `_multi_fast5_tmp/`
- `fill_missing_pod5_metadata()` - After conversion, fills empty `flow_cell_product_code` / `sequencing_kit` / `sample_rate` in pod5 RunInfo from the source fast5 (via the checker's `extract_chemistry_fast5`); never overwrites, never fails the conversion
- `get_output_name()` - Derives output folder name from the resolved input path, stripping archive extensions
- `is_converted()` / `converted_source()` / `_overlap_problem()` - `main()` skips converted inputs unless `--force`, refuses output converted from a different source (`pod5/.source`), refuses output on top of the input, dedupes resolved inputs, and fails duplicate output names within one invocation. SIGTERM is turned into `SystemExit` so cleanup runs
- `parse_list_file()` resolves relative entries against the list file's folder and expands `~`
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
- Never trust a conversion tool's exit code alone: `pod5 convert` and `single_to_multi_fast5` exit 0 after skipping files they cannot handle. Verify read counts (converter: `_verify_read_count`; scripts: embedded check)
- `pod5 convert --output DIR/` writes a *file* named DIR if the folder does not exist; always create the output folder first
- Hidden folders (starting with `.`) are work folders and never run data
- Tests must not depend on `os.scandir` order; `test_optimizations.py`'s runner gives each test its own directory
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
- Verified with real dorado 0.9.6: a model *name* passed to `dorado basecaller` is treated as a folder path (download it first with `dorado download --model NAME`, run from the same folder); no dorado reads single-read fast5; 0.9.6 reads multi-read fast5 only with an explicit model; auto-selection on pod5 without a flowcell code aborts with `__UNKNOWN_FLOWCELL__`
- `FLOWCELL_PORE` dict maps flowcell product codes to pore generations; update when new flowcells are released
- `KIT_PORE` dict maps sequencing kit codes to pore types; used as fallback when flowcell code is absent
