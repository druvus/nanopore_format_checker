# Nanopore Run Format Checker

Scans a folder of Oxford Nanopore sequencing run directories and classifies each run by its read format. Detects flowcell chemistry and recommends the appropriate dorado basecaller version. Optionally generates conversion scripts with metadata patching.

## Features

- **Format detection**: Classifies runs as `pod5`, `multi_read_fast5`, `single_read_fast5`, `fastq`, or compressed archives
- **Chemistry detection**: Identifies pore type (R9.4.1, R10.3, R10.4, R10.4.1, RNA004) from flowcell codes, kit codes, or sample rate
- **Dorado version guidance**: Recommends dorado 0.9.6 for R9.4.1/R10.3 data, dorado >=1.0 for R10.4.1
- **Conversion scripts**: Generates bash scripts for fast5-to-pod5 conversion with correct metadata
- **TSV export**: Per-run statistics including format, file counts, data sizes, and chemistry

## Installation

```bash
pip install h5py ont_fast5_api pod5
```

`h5py` is the only hard dependency. `pod5` is optional (enables pod5 chemistry extraction and metadata patching in conversion scripts).

## Usage

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

### Options

| Flag | Short | Description |
|------|-------|-------------|
| `--verbose` | `-v` | Show detailed per-run info and conversion hints |
| `--quick` | | Skip file counting and size calculation |
| `--convert-to {pod5,single_fast5}` | | Generate a bash conversion script |
| `--script-output FILE` | | Path for generated script (default: `convert_runs.sh`) |
| `--output-dir DIR` | | Base output directory for converted files |
| `--output-stats FILE` | `-o` | Write per-run statistics to a TSV file |

### Example output

```
Found 12 nanopore run(s) and 0 archive(s) in '/data/sequencing'

Run / Archive                                           Format(s)                 Chemistry          Data size    Files
------------------------------------------------------------------------------------------------------------------------
20240115_PAM12345_sample1                               pod5                      R10.4.1 5kHz       42.3 GB      120
20231020_FAP92655_sample2                               multi_read_fast5          R10.4.1 5kHz       28.1 GB      85
20200210_MN19414_old_run                                single_read_fast5         R9.4.1 4kHz        ~12.5 GB     156203
20180501_MN17089_very_old                               multi_read_fast5          R9.4.1 4kHz        8.2 GB       24
```

## Format detection

The tool identifies run directories by the `YYYYMMDD_` naming convention and searches for format-specific subdirectories up to 5 levels deep:

- **pod5**: `pod5/`, `pod5_pass/`, `pod5_fail/` directories
- **fast5**: `fast5/`, `fast5_pass/`, `fast5_fail/` directories, or `.fast5` files in the run root or numeric subdirectories
- **fastq**: `fastq_pass/`, `fastq_fail/` directories
- **Archives**: `.tar.gz`, `.tgz`, `.tar`, `.gz` files (treated as single-read fast5)

Fast5 files are classified as single-read (<1 MB) or multi-read (>=1 MB) based on file size. Barcoded layouts (e.g., `fast5_pass/barcode01/`) are handled automatically.

## Chemistry detection

Chemistry is extracted from one data file per run using a multi-level fallback:

1. **Flowcell product code** from `context_tags` or `tracking_id` -- mapped via `FLOWCELL_PORE` table (23 flowcell codes)
2. **Sequencing kit code** -- mapped via `KIT_PORE` table (35 kit codes) when flowcell is absent
3. **Sample rate** -- last resort for very old runs with no metadata (5000 Hz = R10.4.1, all others = R9.4.1)

### Supported pore types

| Pore | Dorado version | Example flowcells |
|------|---------------|-------------------|
| R10.4.1 | >=1.0 | FLO-MIN114, FLO-PRO114, FLO-MIN114HD |
| R10.4 | >=1.0 (5kHz) / 0.9.6 (4kHz) | FLO-MIN112, FLO-PRO112 |
| R10.3 | 0.9.6 | FLO-MIN111 |
| R9.4.1 | 0.9.6 | FLO-MIN106, FLO-MIN107, FLO-PRO002 |
| RNA004 | >=1.0 | FLO-MIN004RA, FLO-PRO004RA |

## Conversion scripts

Generated scripts handle:

- **multi_read_fast5 to pod5**: Direct conversion with `pod5 convert fast5 --threads 20 --recursive`
- **single_read_fast5 to pod5**: Two-step workflow (`single_to_multi_fast5` then `pod5 convert fast5`)
- **Metadata patching**: When chemistry was detected from the source fast5 files, the script includes a post-conversion step that injects the correct flowcell, kit, and sample rate into the pod5 RunInfo using the pod5 Python API. This is necessary for old runs where the converter produces pod5 files with empty or incorrect metadata.

## Tests

```bash
python test_optimizations.py
```

75 tests covering format detection, chemistry extraction, conversion script generation, and metadata patching.
