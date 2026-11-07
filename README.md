# Nanopore Format Checker & Converter

Tools for inspecting and converting Oxford Nanopore sequencing data.

## Features

- **Format detection**: Classifies runs as `pod5`, `multi_read_fast5`, `single_read_fast5`, `fastq`, or compressed archives
- **Chemistry detection**: Identifies pore type (R9.4.1, R9.5.1, R10.3, R10.4, R10.4.1, RNA002, RNA004) from flowcell codes, kit codes, or sample rate
- **Dorado version guidance**: Recommends dorado >=1.0 for R10.4.1 5 kHz and RNA004, dorado 0.9.6 for R9.4.1, R10.4.1 4 kHz and RNA002, and flags chemistries no dorado release supports (R10.3, R10.4)
- **Conversion scripts**: Generates bash scripts for fast5-to-pod5 conversion with correct metadata
- **Direct conversion**: Converts fast5 to pod5 directly, handling single-read and multi-read formats automatically
- **Flexible input**: Accepts folders, compressed archives, or text files listing multiple inputs
- **TSV export**: Per-run statistics including format, file counts, data sizes, and chemistry

## Installation

```bash
pip install h5py ont_fast5_api pod5
```

`h5py` is the only hard dependency for the format checker. The converter requires `ont_fast5_api` (provides `single_to_multi_fast5`) and `pod5` (provides `pod5 convert fast5`) to be installed and on PATH.

Requires Python 3.12 or newer; the scripts exit with an error on older versions. On current Python versions `pip install ont_fast5_api` fails: its pinned `progressbar33` dependency no longer builds, and it imports `pkg_resources`, which setuptools 81+ no longer ships. Install it without dependencies, with the maintained `progressbar2` and an older setuptools (this is what CI does):

```bash
pip install h5py pod5 progressbar2 packaging "setuptools<81"
pip install --no-deps ont_fast5_api
```

To basecall the results you need [dorado](https://github.com/nanoporetech/dorado): the version depends on the chemistry, see [Dorado version recommendations](#dorado-version-recommendations).

## Usage

### Format checker

```bash
# Scan and classify all runs
python nanopore_format_checker.py /path/to/runs_folder

# Detailed output with chemistry info and conversion hints
python nanopore_format_checker.py /path/to/runs_folder --verbose

# Fast scan (skip file counting and size calculation)
python nanopore_format_checker.py /path/to/runs_folder --quick

# Generate a conversion script to pod5
python nanopore_format_checker.py /path/to/runs_folder --convert-to pod5

# Convert to a separate output directory
python nanopore_format_checker.py /path/to/runs_folder --convert-to pod5 --output-dir /scratch/converted

# Export per-run statistics to TSV
python nanopore_format_checker.py /path/to/runs_folder -o stats.tsv
```

#### Checker options

| Flag | Short | Description |
|------|-------|-------------|
| `--verbose` | `-v` | Show detailed per-run info and conversion hints |
| `--quick` | | Skip file counting and size calculation |
| `--convert-to {pod5,single_fast5}` | | Generate a bash conversion script |
| `--script-output FILE` | | Path for generated script (default: `convert_runs.sh`) |
| `--output-dir DIR` | | Base output directory for converted files |
| `--output-stats FILE` | `-o` | Write per-run statistics to a TSV file |

#### Example output

```
Found 12 nanopore run(s) and 0 archive(s) in '/data/sequencing'

Run / Archive                                           Format(s)                 Chemistry          Data size    Files
------------------------------------------------------------------------------------------------------------------------
20240115_PAM12345_sample1                               pod5                      R10.4.1 5kHz       42.3 GB      120
20231020_FAP92655_sample2                               multi_read_fast5          R10.4.1 5kHz       28.1 GB      85
20200210_MN19414_old_run                                single_read_fast5         R9.4.1 4kHz        ~12.5 GB     156203
20180501_MN17089_very_old                               multi_read_fast5          R9.4.1 4kHz        8.2 GB       24
```

### Converter

```bash
# Convert a single run folder
python nanopore_converter.py /path/to/run_folder --output-dir /scratch/converted

# Convert from a compressed archive
python nanopore_converter.py /path/to/run.tar.gz --output-dir /scratch/converted

# Convert multiple runs from a list file
python nanopore_converter.py input_list.txt --output-dir /scratch/converted

# With more threads and verbose output
python nanopore_converter.py /path/to/input --output-dir /scratch/converted --threads 8 --verbose

# Re-convert inputs that were already converted
python nanopore_converter.py input_list.txt --output-dir /scratch/converted --force
```

The converter auto-detects the input type and opens **every** fast5 file (about 0.2 ms each) to sort it into single-read, multi-read or unreadable, so runs that mix both kinds (e.g. single-read `fast5_pass/` with multi-read `fast5_fail/`) are converted completely:
- **single-read files**: merged with `single_to_multi_fast5`, then converted with `pod5 convert fast5`
- **multi-read files**: converted with `pod5 convert fast5` directly
- **Unreadable/corrupt files** (e.g. truncated after a MinKNOW crash): skipped with a warning naming them, and counted again in the final summary
- **Read-count verification**: both tools exit 0 after silently skipping files they cannot handle, so the pod5 output must contain exactly as many reads as the readable input files, otherwise the item fails and nothing is published. Tool warnings and errors are shown even when the tool succeeds.
- **Metadata fill**: after conversion, empty flowcell / kit / sample-rate fields in the pod5 RunInfo are filled from the source fast5 (values pod5 wrote are never changed). `pod5 convert` only reads `tracking_id`, so without this older runs that store the flowcell in `context_tags` lose it, and dorado's automatic model selection fails. Needs `nanopore_format_checker.py` next to the converter; skipped with a warning otherwise.
- **Compressed archives** (`.tar.gz`, `.tgz`, `.tar.bz2`, `.tbz2`, `.tar.xz`, `.txz`, `.tar`, `.zip`): safely extracted to a temp directory, converted, then cleaned up; output is named after the archive

Output is written to `<output-dir>/<input_name>/pod5/`.

Re-running is safe: inputs whose `pod5/` folder already holds `.pod5` files are skipped (archives are not even extracted), so an interrupted batch can simply be started again. Use `--force` to re-convert them. Each conversion is written to a hidden `.pod5_partial/` folder and only moved to `pod5/` once it succeeds, so a failed or interrupted conversion never leaves a half-written `pod5/` behind, and with `--force` the previous output is kept if the new conversion fails. Two inputs that map to the same output name in one invocation are reported as a failure. Each `pod5/` records its source in `pod5/.source`, so a *different* input with the same name (e.g. `run2/fast5_pass` after `run1/fast5_pass`) is refused instead of being skipped or overwritten. An `--output-dir` that would place the output on top of the input (e.g. the input's parent folder, which could replace the run's own `pod5/`) is refused; an output folder inside the input is fine and excluded from the scan.

On Ctrl-C or SIGTERM (e.g. a batch-scheduler time limit) the conversion tools and their worker processes are stopped and temporary folders, including extracted archives, are removed.

The input list file supports one path per line (folders and archives can be mixed), with blank lines and `#` comments skipped. Relative entries are resolved against the list file's own folder, `~` is expanded, and duplicate entries are converted once.

#### Converter options

| Flag | Short | Description |
|------|-------|-------------|
| `input` | | Positional: folder, compressed archive, or text file listing paths |
| `--output-dir DIR` | | Required. Base output directory for converted pod5 files |
| `--threads N` | | Thread count for conversion tools (default: 4) |
| `--force` | | Re-convert inputs that already have output (default: skip them) |
| `--verbose` | `-v` | Enable detailed logging |

#### Exit codes

| Code | Meaning |
|------|---------|
| 0 | All conversions succeeded |
| 1 | Some conversions failed (summary printed) |
| 2 | Fatal error (bad arguments, missing tools) |

### Archiving run folders

`compress_folders.py` packs every subfolder of an input directory into its own `.tar.gz` in an output directory (copying existing archives and recompressing plain `.tar` files), in parallel:

```bash
python compress_folders.py /data/raw /data/compressed --dry-run   # preview
python compress_folders.py /data/raw /data/compressed -w 4        # 4 parallel workers
```

- **Safe re-runs:** each archive is written to a temp file and renamed into place when complete, and items whose archive already exists are skipped.
- **Symlinked folders** at the top level are archived with their content (not as a bare link).
- **Name clashes:** two items that would produce the same archive (e.g. folder `run1/` and file `run1.tar.gz`) are both reported as errors and neither is written. An output directory inside the input directory is never archived into itself.
- **Ctrl-C** cancels queued items and stops archives in progress (their temp files are removed).
- **Unreadable files and subfolders** are skipped with a warning and counted in the summary; the rest of the folder is still archived, and an unreadable subfolder is kept as an empty directory entry. A hard link to a skipped file is stored as a regular file, so the archive still extracts. A top-level folder that cannot be listed or entered is reported as an error and no archive is written. Delete an archive to rebuild it after fixing permissions.
- **Exit code** is 1 if any item failed.

## Format detection

The tool identifies run directories by the `YYYYMMDD_` naming convention and searches for format-specific subdirectories up to 5 levels deep:

- **pod5**: `pod5/`, `pod5_pass/`, `pod5_fail/` directories
- **fast5**: `fast5/`, `fast5_pass/`, `fast5_fail/` directories (with barcode/numeric subfolders up to 4 levels deep), or `.fast5` files in the run root or numeric subdirectories
- **fastq**: `fastq_pass/`, `fastq_fail/` directories
- **Archives**: `.tar.gz`, `.tgz`, `.tar.bz2`, `.tar.xz`, `.tar`, `.zip`, `.gz` files (treated as single-read fast5) -- only when no pod5/fast5 data was found; gzipped text/sequence outputs such as `sequencing_summary.txt.gz` or `*.fastq.gz` are ignored

Fast5 runs are classified from up to five sample files: the first one whose HDF5 root group h5py can read decides (`read_*` groups = multi-read; `Raw/`/`UniqueGlobalKey/` = single-read), so one truncated file does not misclassify a run. Only if no sample can be parsed (or h5py is unavailable) is a size heuristic used (<1 MB single-read, >=1 MB multi-read). Symlinked data folders are followed (with a loop guard), and hidden folders such as the converters' `.convert_tmp/` and `.pod5_partial/` work folders are ignored. A run folder and a top-level archive with the same name are reported separately. Barcoded layouts (e.g., `fast5_pass/barcode01/`) are handled automatically.

## Chemistry detection

Chemistry is extracted from one data file per run using a multi-level fallback:

1. **Flowcell product code** from `context_tags` or `tracking_id` -- mapped via `FLOWCELL_PORE` table (21 flowcell codes)
2. **Sequencing kit code** -- mapped via `KIT_PORE` table (55 kit codes) when flowcell is absent
3. **Sample rate** -- last resort for very old runs with no metadata (5000 Hz = R10.4.1, all others = R9.4.1)

### Dorado version recommendations

Dorado 1.0 (May 2025) removed R9.4.1, R10.4.1 4 kHz and RNA002 models, **and removed FAST5 input**. Dorado 0.9.6 is the last release that supports them. Dorado never shipped models for R10.3 or R10.4. Sources: the [dorado changelog](https://github.com/nanoporetech/dorado/blob/master/CHANGELOG.md) and the chemistry table in dorado 0.9.6 (`dorado/models/kits.cpp`).

| Pore / condition | Dorado version | Example flowcells / kits | Explicit model (if auto-selection fails) |
|------------------|----------------|--------------------------|------------------------------------------|
| R10.4.1, 5 kHz | >=1.0 | FLO-MIN114, FLO-PRO114, FLO-MIN114HD | -- |
| R10.4.1, 4 kHz (400 bps) | 0.9.6 | FLO-MIN114 + SQK-LSK114 at 4 kHz | `dna_r10.4.1_e8.2_400bps_sup@v4.1.0` |
| R10.4.1, 260 bps | 0.9.6 | kits ending `-260` | `dna_r10.4.1_e8.2_260bps_sup@v4.1.0` |
| R9.4.1 (4 kHz) | 0.9.6 | FLO-MIN106, FLO-PRO002, FLO-FLG001; Kit 9/10/11 | `dna_r9.4.1_e8_sup@v3.6` |
| RNA002 | 0.9.6 | SQK-RNA002 | `rna002_70bps_hac@v3` (no sup model exists) |
| RNA004 | >=1.0 | FLO-MIN004RA, FLO-PRO004RA, SQK-RNA004 | -- |
| R10.4 | none (legacy Guppy) | FLO-MIN112, FLO-PRO112; Kit 12 | -- |
| R10.3 | none (legacy Guppy) | FLO-MIN111 | -- |
| R9.5.1 (DNA) | none (legacy Guppy) | FLO-MIN107 | -- |

How to basecall once you know the version:

- **dorado >=1.0** needs POD5 input. Convert fast5 runs first (`--convert-to pod5` or `nanopore_converter.py`), then run `dorado basecaller sup pod5_dir/ > calls.bam`.
- **dorado 0.9.6**: automatic model selection (`dorado basecaller sup ...`) only works on POD5 whose RunInfo has a flowcell + kit pair and sample rate that dorado recognises; otherwise dorado stops with `No supported chemistry found for flowcell_code: '__UNKNOWN_FLOWCELL__'`. The conversion scripts' metadata patching and the converter's metadata fill exist to make auto-selection work for old runs whose fast5 lacks these fields. If a run has no flowcell code at all, use the explicit model from the table. dorado treats a model *name* as a folder path, so download it first and run basecaller from the same folder:

  ```bash
  dorado download --model dna_r9.4.1_e8_sup@v3.6
  dorado basecaller dna_r9.4.1_e8_sup@v3.6 pod5_dir/ > calls.bam
  ```
- **fast5 input**: no dorado version reads single-read fast5, so convert it to pod5. dorado 0.9.6 reads multi-read fast5 only with an explicit model (as above); converting to pod5 is still recommended.

These points were checked by running dorado 0.9.6 on converted test data: automatic selection worked for pod5 with full or filled-in metadata, failed for pod5 without a flowcell code (and for any fast5), and the two-step explicit-model command worked. dorado 2.1.0 auto-selected `dna_r10.4.1_e8.2_400bps_fast@v5.2.0` for converted R10.4.1 5 kHz data and refused R9.4.1 data ("chemistry has been deprecated since Dorado v1.0.0").
- An R10.4.1 run whose sample rate could not be read is reported as `>=1.0` with a note to use 0.9.6 if it was sampled at 4 kHz. R9.4.1 runs at a rate other than 4 kHz (e.g. early 6 kHz PromethION) are flagged, as dorado's R9.4.1 model expects 4 kHz.

`--verbose` prints this advice per run, and the TSV includes `dorado_version`, `dorado_model` and `dorado_note` columns.

## Conversion scripts

Generated scripts handle:

- **multi_read_fast5 to pod5**: Direct conversion of the detected fast5 folders with `pod5 convert fast5 --threads 20 --recursive`
- **single_read_fast5 to pod5**: Two-step workflow (`single_to_multi_fast5` per detected folder, then `pod5 convert fast5`)
- **Output**: `<output-dir>/<run>/pod5/`, or `<run>/pod5_from_fast5/` without `--output-dir`. Paths are absolute, so the script can be run from any directory.
- **Safe re-runs**: runs whose output already holds `.pod5` files are skipped; work happens in a hidden `.convert_tmp/` folder that is moved into place only on success.
- **Failures are isolated**: each run executes fail-fast in its own subshell; a failed run is reported, the script continues with the next one and exits 1 at the end, listing the failed runs.
- **Read-count verification**: before publishing, an embedded check compares the reads in the fast5 input with the pod5 output. The script detects a run's format from a few sample files, so a run that mixes single- and multi-read files (or has corrupt files) fails this check rather than losing reads silently; convert such runs with `nanopore_converter.py`, which classifies every file.
- **Metadata patching**: When chemistry was detected from the source fast5 files, the script includes a post-conversion step that writes the detected flowcell, kit, and sample rate into the pod5 RunInfo using the pod5 Python API. dorado's automatic model selection (`dorado basecaller sup ...`) looks these three fields up in its chemistry table, so without them older runs can only be basecalled with an explicit model name. A run with no flowcell code anywhere in its fast5 still needs the explicit model; the tool does not guess one. Files that already carry the right values (compared case-insensitively, as dorado does) are not rewritten.
- **Unsupported conversions are not emitted**: runs whose data is only in compressed archives become `# SKIPPED` comments; `--convert-to single_fast5` only converts multi-read fast5 (one output folder per input folder, so pass and fail stay separate).

## Tests

```bash
# Format checker tests (94 tests)
python test_optimizations.py

# Converter tests (46 tests)
python -m pytest test_nanopore_converter.py -v

# Archiver tests (12 tests; 2 need a non-root user)
python -m pytest test_compress_folders.py -v

# End-to-end tests with the real pod5 / single_to_multi_fast5 tools (5 tests;
# skipped if the tools are not installed)
python -m pytest test_integration.py -v
```

CI (`.github/workflows/tests.yml`) installs the real conversion tools and runs all four suites with pytest on Python 3.12, 3.13 and 3.14 for every pull request and push to `main`.

The format checker tests cover format detection, chemistry extraction, the dorado support matrix, conversion script generation, and metadata patching. The converter tests cover input detection, archive extraction, fast5 classification, conversion logic, and end-to-end pipeline integration.
