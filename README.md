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

Requires Python 3.12 or newer; the scripts exit with an error on older versions. On Python 3.13, `pip install ont_fast5_api` fails because its `progressbar33` dependency no longer builds; install it without dependencies and use the maintained `progressbar2` instead:

```bash
pip install h5py pod5 progressbar2
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

The converter auto-detects the input type and fast5 format:
- **single_read_fast5**: runs `single_to_multi_fast5` then `pod5 convert fast5`, cleans up intermediate files (fails the item if the merge step produces no output)
- **multi_read_fast5**: runs `pod5 convert fast5` directly
- **Metadata fill**: after conversion, empty flowcell / kit / sample-rate fields in the pod5 RunInfo are filled from the source fast5 (values pod5 wrote are never changed). `pod5 convert` only reads `tracking_id`, so without this older runs that store the flowcell in `context_tags` lose it, and dorado's automatic model selection fails. Needs `nanopore_format_checker.py` next to the converter; skipped with a warning otherwise.
- **Compressed archives** (`.tar.gz`, `.tgz`, `.tar.bz2`, `.tar`, `.zip`): safely extracted to a temp directory, converted, then cleaned up; output is named after the archive

Output is written to `<output-dir>/<input_name>/pod5/`.

Re-running is safe: inputs whose `pod5/` folder already holds `.pod5` files are skipped (archives are not even extracted), so an interrupted batch can simply be started again. Use `--force` to re-convert them. Each conversion is written to a hidden `.pod5_partial/` folder and only moved to `pod5/` once it succeeds, so a failed or interrupted conversion never leaves a half-written `pod5/` behind, and with `--force` the previous output is kept if the new conversion fails. Two inputs that map to the same output name in one invocation are reported as a failure.

The input list file supports one path per line (folders and archives can be mixed), with blank lines and `#` comments skipped.

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

## Format detection

The tool identifies run directories by the `YYYYMMDD_` naming convention and searches for format-specific subdirectories up to 5 levels deep:

- **pod5**: `pod5/`, `pod5_pass/`, `pod5_fail/` directories
- **fast5**: `fast5/`, `fast5_pass/`, `fast5_fail/` directories, or `.fast5` files in the run root or numeric subdirectories
- **fastq**: `fastq_pass/`, `fastq_fail/` directories
- **Archives**: `.tar.gz`, `.tgz`, `.tar`, `.gz` files (treated as single-read fast5) -- only when no pod5/fast5 data was found; gzipped text/sequence outputs such as `sequencing_summary.txt.gz` or `*.fastq.gz` are ignored

Fast5 files are classified by inspecting the HDF5 root group of one sampled file (`read_*` groups = multi-read; `Raw/`/`UniqueGlobalKey/` = single-read). If h5py is unavailable or the file cannot be parsed, a size heuristic is used (<1 MB single-read, >=1 MB multi-read). Barcoded layouts (e.g., `fast5_pass/barcode01/`) are handled automatically.

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
- **dorado 0.9.6** still reads fast5, but automatic model selection (`dorado basecaller sup ...`) only works on POD5 whose RunInfo has a flowcell + kit pair and sample rate that dorado recognises. Otherwise pass the explicit model name from the table, e.g. `dorado basecaller dna_r9.4.1_e8_sup@v3.6 reads/ > calls.bam`. The conversion scripts' metadata patching exists to make auto-selection work for old runs whose fast5 lacks these fields.
- An R10.4.1 run whose sample rate could not be read is reported as `>=1.0` with a note to use 0.9.6 if it was sampled at 4 kHz. R9.4.1 runs at a rate other than 4 kHz (e.g. early 6 kHz PromethION) are flagged, as dorado's R9.4.1 model expects 4 kHz.

`--verbose` prints this advice per run, and the TSV includes `dorado_version`, `dorado_model` and `dorado_note` columns.

## Conversion scripts

Generated scripts handle:

- **multi_read_fast5 to pod5**: Direct conversion with `pod5 convert fast5 --threads 20 --recursive`
- **single_read_fast5 to pod5**: Two-step workflow (`single_to_multi_fast5` then `pod5 convert fast5`)
- **Metadata patching**: When chemistry was detected from the source fast5 files, the script includes a post-conversion step that writes the detected flowcell, kit, and sample rate into the pod5 RunInfo using the pod5 Python API. dorado's automatic model selection (`dorado basecaller sup ...`) looks these three fields up in its chemistry table, so without them older runs can only be basecalled with an explicit model name. A run with no flowcell code anywhere in its fast5 still needs the explicit model; the tool does not guess one.
- **Unsupported conversions are not emitted**: runs whose data is only in compressed archives become `# SKIPPED` comments; `--convert-to single_fast5` only converts multi-read fast5.

## Tests

```bash
# Format checker tests (81 tests)
python test_optimizations.py

# Converter tests (36 tests)
python -m pytest test_nanopore_converter.py -v
```

CI (`.github/workflows/tests.yml`) runs both suites with pytest on Python 3.12, 3.13 and 3.14 for every pull request and push to `main`.

The format checker tests cover format detection, chemistry extraction, the dorado support matrix, conversion script generation, and metadata patching. The converter tests cover input detection, archive extraction, fast5 classification, conversion logic, and end-to-end pipeline integration.
