# Chemistry Detection Design

Detect flowcell chemistry (R9.4.1, R10.4.1, RNA) from pod5/fast5 file
metadata and recommend the appropriate dorado version and model for
basecalling.

## Motivation

Users need to know which dorado version to use for basecalling. R9.4.1
support was dropped in dorado 1.0, so R9.4.1 data requires dorado 0.9.6.
R10.4.1 at 4kHz (older MinKNOW) also requires 0.9.6. The tool already
opens run directories and samples files -- extracting chemistry from those
files adds minimal overhead.

## Data sources

Chemistry information is embedded in the sequencing data files:

- **Pod5 `RunInfo`**: `flow_cell_product_code`, `sequencing_kit`,
  `sample_rate` via the pod5 Python API. This is the same source dorado
  uses for automatic model selection.
- **Fast5 HDF5 attributes**: `context_tags/flowcell_type`,
  `context_tags/sequencing_kit`, `context_tags/sample_frequency` via h5py.
  Values are stored as bytes in lowercase and need decoding/normalization.
- **Archives**: no metadata accessible without extraction.

### Fast5 HDF5 structure differences

- Single-read: attributes at `/UniqueGlobalKey/context_tags/` and
  `/UniqueGlobalKey/tracking_id/`
- Multi-read: attributes at `/<read_id>/context_tags/` and
  `/<read_id>/tracking_id/`

## New functions

### `extract_chemistry_pod5(file_path) -> dict | None`

Opens one pod5 file, reads first read's `run_info`. Returns
`{"flowcell": "FLO-MIN114", "kit": "SQK-LSK114", "sample_rate": 5000}`
or `None`. Requires `HAS_POD5`.

### `extract_chemistry_fast5(file_path) -> dict | None`

Opens one fast5 with h5py. Tries `UniqueGlobalKey/context_tags`
(single-read), then `<first_key>/context_tags` (multi-read). Decodes
bytes, normalizes to uppercase. Returns same dict shape. Requires
`HAS_H5PY`.

### `extract_chemistry(file_path, file_format) -> dict | None`

Dispatcher: calls pod5 or fast5 variant based on format. Wrapped in
try/except -- never blocks format detection on failure.

### `classify_chemistry(chemistry) -> dict`

Pure lookup, no I/O. Maps raw chemistry dict to:

```python
{
    "pore": "R10.4.1",
    "analyte": "dna",
    "dorado_version": ">=1.0",
    "model_hint": "sup",
    "note": None,
}
```

## Flowcell lookup table

```python
FLOWCELL_PORE = {
    # R10.4.1
    "FLO-MIN114": "R10.4.1",
    "FLO-PRO114": "R10.4.1",
    "FLO-FLG114": "R10.4.1",
    "FLO-PRO114M": "R10.4.1",
    # R9.4.1
    "FLO-MIN106": "R9.4.1",
    "FLO-MIN106D": "R9.4.1",
    "FLO-MIN107": "R9.4.1",
    "FLO-FLG001": "R9.4.1",
    "FLO-PRO002": "R9.4.1",
    "FLO-PRO002M": "R9.4.1",
}
```

## Dorado version logic

- R9.4.1 -> dorado 0.9.6
- R10.4.1 + sample_rate < 5000 -> dorado 0.9.6 (4kHz legacy)
- R10.4.1 + sample_rate >= 5000 -> dorado >=1.0
- RNA002 kit (SQK-RNA002) -> dorado 0.9.6
- RNA004 kit (SQK-RNA004) -> dorado >=1.0
- Unknown flowcell -> report raw values, no recommendation

## Integration into analyze_run()

- After pod5 detection: find first `.pod5` file via `os.scandir`, call
  `extract_chemistry_pod5()`.
- After fast5 detection: use the already-sampled `sample_file` Path, call
  `extract_chemistry_fast5()`. Fix root/numeric path to store full Paths
  (currently stores bare filename).
- Store at `result["chemistry"]` (run-level, not per-format).
- Prefer pod5 chemistry when both pod5 and fast5 exist.

## Output changes

### Main table

Add "Chemistry" column: compact string like `R10.4.1 5kHz`, `R9.4.1`,
`RNA004`, or `-` when unavailable.

### Verbose mode

```
  chemistry: R9.4.1 (FLO-MIN106, SQK-LSK109, 4000 Hz)
  dorado: use version 0.9.6 -- R9.4.1 support was dropped in 1.0
  model: dorado basecaller sup <input>/
```

### TSV output

Add columns: `flowcell_code`, `sequencing_kit`, `sample_rate`,
`pore_type`, `dorado_version`.

## Behavior

- Always on (no opt-in flag). Opens one file per run.
- Best-effort: failures in chemistry extraction are logged in verbose mode
  but never prevent format detection.
- `HAS_POD5` / `HAS_H5PY` guards: skip extraction when libraries are
  missing.

## Tests

- `test_extract_chemistry_fast5_multi_read` -- mock multi-read fast5,
  verify extraction
- `test_extract_chemistry_fast5_single_read` -- mock single-read fast5,
  verify extraction
- `test_classify_chemistry_r10` -- R10.4.1 5kHz -> dorado >=1.0
- `test_classify_chemistry_r9` -- R9.4.1 -> dorado 0.9.6
- `test_classify_chemistry_r10_4khz` -- R10.4.1 4kHz -> dorado 0.9.6
- `test_classify_chemistry_rna` -- RNA004 -> dorado >=1.0, RNA002 -> 0.9.6
- `test_classify_chemistry_unknown` -- unknown flowcell -> no recommendation
- `test_analyze_run_includes_chemistry` -- full integration with mock
  fast5 file
- `test_tsv_includes_chemistry_columns` -- verify new TSV columns
