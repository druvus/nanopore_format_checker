# Nanopore Converter Design

**Date:** 2026-03-03
**Status:** Approved

## Overview

Standalone Python CLI script (`nanopore_converter.py`) that converts Oxford Nanopore
single_read_fast5 and multi_read_fast5 files to pod5 format. Handles folders,
compressed archives, and input list files.

## Requirements

- Accept folder, compressed file (.tar.gz, .tar.bz2, .zip), or text file with
  paths (one per line, mixed folders and archives allowed)
- Convert single_read_fast5 via two-step process: single_to_multi_fast5, then
  pod5 convert fast5
- Convert multi_read_fast5 directly via pod5 convert fast5
- Output to a user-specified directory using input basename as subdirectory name
- Extract compressed archives to temp directory, convert, then clean up
- Call conversion tools as subprocesses (not Python API)
- No metadata patching or basecalling
- Sequential processing (conversion tools are already multi-threaded)

## CLI Interface

```
python nanopore_converter.py INPUT --output-dir /scratch/converted [--threads 4] [--verbose]
```

- `INPUT`: folder path, compressed archive, or text file listing paths
- `--output-dir`: required, base directory for converted pod5 output
- `--threads`: thread count passed to conversion tools (default: 4)
- `--verbose`: detailed logging

## Input Detection

1. **Directory** -- walk to find .fast5 files, classify as single or multi-read
   using the 1MB file size threshold
2. **Compressed archive** -- extract to tempfile.mkdtemp(), treat as directory
3. **Text file** (.txt, .csv, .list) -- read lines, each line is a folder or
   archive path

## Output Structure

```
<output-dir>/
  <input_basename>/
    pod5/
      *.pod5
```

For example:
- Input: `/data/20230101_run1/` -> `<output-dir>/20230101_run1/pod5/`
- Input: `/data/20230101_run1.tar.gz` -> `<output-dir>/20230101_run1/pod5/`

## Conversion Pipeline

```
detect_format(path)
  |
  +-- single_read_fast5 -> single_to_multi_fast5 -> pod5 convert fast5 -> cleanup temp
  +-- multi_read_fast5  -> pod5 convert fast5
  +-- unknown/empty     -> skip with warning
```

Subprocess commands:
- `single_to_multi_fast5 -i <input> -s <temp_multi> -t <threads> --recursive`
- `pod5 convert fast5 <input>/ --output <output>/pod5/ --threads <threads> --recursive`

## Format Detection

Classify fast5 files using file size heuristic (same as nanopore_format_checker.py):
- Files < 1 MB: single-read fast5 (typically 1-50 KB)
- Files >= 1 MB: multi-read fast5 (typically 1+ MB)

Sample one .fast5 file from the directory to determine format.

## Error Handling

- Verify `single_to_multi_fast5` and `pod5` are on PATH at startup
- Each run wrapped in try/except; failures log and continue
- Subprocess return codes checked; non-zero raises per-run error
- Temp directories cleaned in finally blocks
- Summary of successes/failures printed at end

## Exit Codes

- 0: all conversions succeeded
- 1: some conversions failed
- 2: fatal error (bad arguments, missing tools)

## Dependencies

- `ont_fast5_api` (provides `single_to_multi_fast5` CLI)
- `pod5` (provides `pod5 convert fast5` CLI)
- Python standard library only for the script itself (subprocess, tempfile,
  tarfile, zipfile, shutil, argparse, logging, pathlib, os)
