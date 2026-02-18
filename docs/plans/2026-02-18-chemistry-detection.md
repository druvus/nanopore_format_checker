# Chemistry Detection Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Detect flowcell chemistry from pod5/fast5 metadata and recommend the correct dorado version and basecalling model.

**Architecture:** Sample one data file per run (already done for fast5 classification, new for pod5). Extract flowcell product code, sequencing kit, and sample rate. Map to pore type via lookup table. Recommend dorado version based on pore + sample rate + kit.

**Tech Stack:** Python, h5py (existing dep), pod5 (optional dep), os.scandir

**Design doc:** `docs/plans/2026-02-18-chemistry-detection-design.md`

---

### Task 1: Add HAS_POD5 guard and FLOWCELL_PORE lookup table

**Files:**
- Modify: `nanopore_format_checker.py:25-29` (imports) and after line 29 (new constants)

**Step 1: Add pod5 import guard after the h5py guard**

After line 29, add:

```python
try:
    import pod5
    HAS_POD5 = True
except ImportError:
    HAS_POD5 = False
```

**Step 2: Add FLOWCELL_PORE lookup table**

After the new `HAS_POD5` block (before `is_nanopore_run_dir`), add:

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

**Step 3: Verify syntax**

Run: `python -c "import nanopore_format_checker"`
Expected: No errors

---

### Task 2: Implement extract_chemistry_fast5()

**Files:**
- Modify: `nanopore_format_checker.py` (new function after `classify_fast5`, ~line 277)
- Test: `test_optimizations.py`

**Step 1: Write failing tests**

Add to `test_optimizations.py` (imports and test functions before `main()`):

```python
def test_extract_chemistry_fast5_single_read(tmp_path: Path) -> None:
    """Extract chemistry from a single-read fast5 file."""
    f5 = tmp_path / "20240101_chem_sr" / "read.fast5"
    f5.parent.mkdir(parents=True)
    import h5py
    with h5py.File(f5, "w") as f:
        ctx = f.create_group("UniqueGlobalKey/context_tags")
        ctx.attrs["flowcell_type"] = "flo-min106"
        ctx.attrs["sequencing_kit"] = "sqk-lsk109"
        ctx.attrs["sample_frequency"] = "4000"
    result = extract_chemistry_fast5(f5)
    assert result is not None
    assert result["flowcell"] == "FLO-MIN106"
    assert result["kit"] == "SQK-LSK109"
    assert result["sample_rate"] == 4000
    print("  PASS: extract_chemistry_fast5 single-read")


def test_extract_chemistry_fast5_multi_read(tmp_path: Path) -> None:
    """Extract chemistry from a multi-read fast5 file."""
    f5 = tmp_path / "20240101_chem_mr" / "read.fast5"
    f5.parent.mkdir(parents=True)
    import h5py
    with h5py.File(f5, "w") as f:
        ctx = f.create_group("read_001/context_tags")
        ctx.attrs["flowcell_type"] = "flo-min114"
        ctx.attrs["sequencing_kit"] = "sqk-lsk114"
        ctx.attrs["sample_frequency"] = "5000"
    result = extract_chemistry_fast5(f5)
    assert result is not None
    assert result["flowcell"] == "FLO-MIN114"
    assert result["kit"] == "SQK-LSK114"
    assert result["sample_rate"] == 5000
    print("  PASS: extract_chemistry_fast5 multi-read")
```

Add `extract_chemistry_fast5` to the import list and the test list in `main()`.

**Step 2: Run tests to verify they fail**

Run: `python test_optimizations.py`
Expected: FAIL -- `extract_chemistry_fast5` not defined

**Step 3: Implement extract_chemistry_fast5()**

Add after `classify_fast5()` (~line 277) in `nanopore_format_checker.py`:

```python
def _decode_attr(val):
    """Decode an HDF5 attribute value to a string."""
    if isinstance(val, bytes):
        return val.decode("utf-8", errors="replace")
    return str(val)


def extract_chemistry_fast5(file_path: Path) -> dict | None:
    """Extract flowcell chemistry metadata from a fast5 file.

    Reads context_tags attributes from the HDF5 file. Handles both
    single-read (UniqueGlobalKey/context_tags) and multi-read
    (<read_id>/context_tags) layouts.

    Returns {"flowcell": "FLO-MIN114", "kit": "SQK-LSK114",
    "sample_rate": 5000} or None on failure.
    """
    if not HAS_H5PY:
        return None
    try:
        with h5py.File(file_path, "r") as f:
            ctx = None
            if "UniqueGlobalKey" in f and "context_tags" in f["UniqueGlobalKey"]:
                ctx = f["UniqueGlobalKey/context_tags"]
            else:
                for key in f.keys():
                    if "context_tags" in f[key]:
                        ctx = f[f"{key}/context_tags"]
                        break
            if ctx is None:
                return None
            flowcell = _decode_attr(ctx.attrs.get("flowcell_type", "")).upper()
            kit = _decode_attr(ctx.attrs.get("sequencing_kit", "")).upper()
            rate_str = _decode_attr(ctx.attrs.get("sample_frequency", "0"))
            sample_rate = int(rate_str) if rate_str.isdigit() else 0
            if not flowcell:
                return None
            return {"flowcell": flowcell, "kit": kit, "sample_rate": sample_rate}
    except Exception:
        return None
```

**Step 4: Run tests to verify they pass**

Run: `python test_optimizations.py`
Expected: All pass including the two new chemistry fast5 tests

**Step 5: Commit**

```bash
git add nanopore_format_checker.py test_optimizations.py
git commit -m "Add extract_chemistry_fast5() for reading flowcell metadata from fast5 files"
```

---

### Task 3: Implement extract_chemistry_pod5()

**Files:**
- Modify: `nanopore_format_checker.py` (new function after `extract_chemistry_fast5`)
- Test: `test_optimizations.py`

**Step 1: Write failing test**

Note: pod5 may not be installed in the test environment. The test should
skip if pod5 is unavailable. However, we can test the `None` return path
and mock-based tests.

```python
def test_extract_chemistry_pod5_missing_lib(tmp_path: Path) -> None:
    """extract_chemistry_pod5 returns None when pod5 is not installed or file is invalid."""
    fake = tmp_path / "20240101_pod5_chem" / "fake.pod5"
    fake.parent.mkdir(parents=True)
    fake.write_bytes(b"\x00" * 100)
    result = extract_chemistry_pod5(fake)
    # Either None (pod5 not installed) or None (invalid file)
    assert result is None
    print("  PASS: extract_chemistry_pod5 missing lib or invalid file")
```

Add `extract_chemistry_pod5` to imports and test list.

**Step 2: Run test to verify it fails**

Run: `python test_optimizations.py`
Expected: FAIL -- `extract_chemistry_pod5` not defined

**Step 3: Implement extract_chemistry_pod5()**

Add after `extract_chemistry_fast5()`:

```python
def extract_chemistry_pod5(file_path: Path) -> dict | None:
    """Extract flowcell chemistry metadata from a pod5 file.

    Reads the RunInfo from the first read. Returns {"flowcell":
    "FLO-MIN114", "kit": "SQK-LSK114", "sample_rate": 5000} or None.
    """
    if not HAS_POD5:
        return None
    try:
        with pod5.Reader(file_path) as reader:
            for read in reader.reads():
                info = read.run_info
                flowcell = info.flow_cell_product_code or ""
                kit = info.sequencing_kit or ""
                sample_rate = info.sample_rate or 0
                if not flowcell:
                    return None
                return {
                    "flowcell": flowcell.upper(),
                    "kit": kit.upper(),
                    "sample_rate": int(sample_rate),
                }
        return None
    except Exception:
        return None
```

**Step 4: Run tests to verify they pass**

Run: `python test_optimizations.py`
Expected: All pass

**Step 5: Commit**

```bash
git add nanopore_format_checker.py test_optimizations.py
git commit -m "Add extract_chemistry_pod5() for reading flowcell metadata from pod5 files"
```

---

### Task 4: Implement extract_chemistry() dispatcher and classify_chemistry()

**Files:**
- Modify: `nanopore_format_checker.py` (two new functions after extract_chemistry_pod5)
- Test: `test_optimizations.py`

**Step 1: Write failing tests**

```python
def test_classify_chemistry_r10_5khz(tmp_path: Path) -> None:
    """R10.4.1 at 5kHz should recommend dorado >=1.0."""
    chem = {"flowcell": "FLO-MIN114", "kit": "SQK-LSK114", "sample_rate": 5000}
    result = classify_chemistry(chem)
    assert result["pore"] == "R10.4.1"
    assert result["dorado_version"] == ">=1.0"
    assert result["analyte"] == "dna"
    print("  PASS: classify_chemistry R10.4.1 5kHz")


def test_classify_chemistry_r9(tmp_path: Path) -> None:
    """R9.4.1 should recommend dorado 0.9.6."""
    chem = {"flowcell": "FLO-MIN106", "kit": "SQK-LSK109", "sample_rate": 4000}
    result = classify_chemistry(chem)
    assert result["pore"] == "R9.4.1"
    assert result["dorado_version"] == "0.9.6"
    assert result["note"] is not None and "0.9.6" in result["note"]
    print("  PASS: classify_chemistry R9.4.1")


def test_classify_chemistry_r10_4khz(tmp_path: Path) -> None:
    """R10.4.1 at 4kHz (legacy) should recommend dorado 0.9.6."""
    chem = {"flowcell": "FLO-MIN114", "kit": "SQK-LSK114", "sample_rate": 4000}
    result = classify_chemistry(chem)
    assert result["pore"] == "R10.4.1"
    assert result["dorado_version"] == "0.9.6"
    print("  PASS: classify_chemistry R10.4.1 4kHz")


def test_classify_chemistry_rna004(tmp_path: Path) -> None:
    """RNA004 kit should recommend dorado >=1.0."""
    chem = {"flowcell": "FLO-MIN114", "kit": "SQK-RNA004", "sample_rate": 4000}
    result = classify_chemistry(chem)
    assert result["analyte"] == "rna"
    assert result["dorado_version"] == ">=1.0"
    print("  PASS: classify_chemistry RNA004")


def test_classify_chemistry_rna002(tmp_path: Path) -> None:
    """RNA002 kit should recommend dorado 0.9.6."""
    chem = {"flowcell": "FLO-MIN106", "kit": "SQK-RNA002", "sample_rate": 3000}
    result = classify_chemistry(chem)
    assert result["analyte"] == "rna"
    assert result["dorado_version"] == "0.9.6"
    print("  PASS: classify_chemistry RNA002")


def test_classify_chemistry_unknown(tmp_path: Path) -> None:
    """Unknown flowcell code should report unknown pore with no version recommendation."""
    chem = {"flowcell": "FLO-PROTOTYPE", "kit": "SQK-BETA", "sample_rate": 5000}
    result = classify_chemistry(chem)
    assert result["pore"] == "unknown"
    assert result["dorado_version"] is None
    print("  PASS: classify_chemistry unknown flowcell")
```

Add `classify_chemistry`, `extract_chemistry` to imports and all tests to list.

**Step 2: Run tests to verify they fail**

Run: `python test_optimizations.py`
Expected: FAIL -- `classify_chemistry` not defined

**Step 3: Implement extract_chemistry() and classify_chemistry()**

Add after `extract_chemistry_pod5()`:

```python
def extract_chemistry(file_path: Path, file_format: str) -> dict | None:
    """Extract flowcell chemistry from a data file.

    Dispatches to the pod5 or fast5 extractor based on file_format.
    Returns None on any failure -- chemistry extraction is best-effort.
    """
    try:
        if file_format == "pod5":
            return extract_chemistry_pod5(file_path)
        elif file_format in ("fast5", "single_read_fast5", "multi_read_fast5"):
            return extract_chemistry_fast5(file_path)
    except Exception:
        pass
    return None


def classify_chemistry(chemistry: dict) -> dict:
    """Map raw chemistry metadata to pore type and dorado version recommendation.

    Input: {"flowcell": "FLO-MIN114", "kit": "SQK-LSK114", "sample_rate": 5000}
    Output: {"pore", "analyte", "dorado_version", "model_hint", "note"}
    """
    flowcell = chemistry.get("flowcell", "")
    kit = chemistry.get("kit", "")
    sample_rate = chemistry.get("sample_rate", 0)

    pore = FLOWCELL_PORE.get(flowcell, "unknown")
    analyte = "rna" if kit.startswith("SQK-RNA") else "dna"

    dorado_version = None
    model_hint = "sup"
    note = None

    # RNA kit overrides
    if kit == "SQK-RNA004":
        dorado_version = ">=1.0"
        note = None
    elif kit == "SQK-RNA002":
        dorado_version = "0.9.6"
        note = "RNA002 support was dropped in dorado 1.0; use dorado 0.9.6"
    elif pore == "R9.4.1":
        dorado_version = "0.9.6"
        note = "R9.4.1 support was dropped in dorado 1.0; use dorado 0.9.6"
    elif pore == "R10.4.1":
        if sample_rate >= 5000:
            dorado_version = ">=1.0"
        else:
            dorado_version = "0.9.6"
            note = "R10.4.1 at 4kHz requires dorado 0.9.6; 5kHz data supported in >=1.0"
    # pore == "unknown" -> dorado_version stays None

    return {
        "pore": pore,
        "analyte": analyte,
        "dorado_version": dorado_version,
        "model_hint": model_hint,
        "note": note,
    }
```

**Step 4: Run tests to verify they pass**

Run: `python test_optimizations.py`
Expected: All pass

**Step 5: Commit**

```bash
git add nanopore_format_checker.py test_optimizations.py
git commit -m "Add extract_chemistry() dispatcher and classify_chemistry() lookup"
```

---

### Task 5: Integrate chemistry extraction into analyze_run()

**Files:**
- Modify: `nanopore_format_checker.py:321-607` (analyze_run function)
- Test: `test_optimizations.py`

**Step 1: Write failing test**

```python
def test_analyze_run_includes_chemistry(tmp_path: Path) -> None:
    """analyze_run() should include chemistry from fast5 metadata."""
    import h5py
    run = tmp_path / "20240101_run_chem"
    fast5_dir = run / "fast5"
    fast5_dir.mkdir(parents=True)
    f5 = fast5_dir / "batch_001.fast5"
    # Create a multi-read fast5 (>1MB) with chemistry metadata
    with h5py.File(f5, "w") as f:
        ctx = f.create_group("read_001/context_tags")
        ctx.attrs["flowcell_type"] = "flo-min106"
        ctx.attrs["sequencing_kit"] = "sqk-lsk109"
        ctx.attrs["sample_frequency"] = "4000"
        # Pad to >1MB for multi-read classification
        f.create_dataset("read_001/Raw/Signal", data=b"\x00" * 1_500_000)

    result = analyze_run(run)
    assert "chemistry" in result, f"Missing chemistry key: {result.keys()}"
    assert result["chemistry"]["flowcell"] == "FLO-MIN106"
    assert result["chemistry"]["kit"] == "SQK-LSK109"
    clas = result.get("chemistry_classification")
    assert clas is not None
    assert clas["pore"] == "R9.4.1"
    assert clas["dorado_version"] == "0.9.6"
    print("  PASS: analyze_run includes chemistry")
```

Add to test list.

**Step 2: Run test to verify it fails**

Run: `python test_optimizations.py`
Expected: FAIL -- "Missing chemistry key"

**Step 3: Integrate into analyze_run()**

Three changes in `analyze_run()`:

**3a.** Initialize chemistry in result dict. Change line 338 from:
```python
result = {"formats": [], "details": {}}
```
to:
```python
result = {"formats": [], "details": {}, "chemistry": None, "chemistry_classification": None}
```

**3b.** After pod5 detection block (after line 401), add pod5 file sampling and chemistry extraction:
```python
    # Extract chemistry from first pod5 file
    if all_pod5_dirs and result["chemistry"] is None:
        for d in all_pod5_dirs:
            if str(d) in unreadable_pod5:
                continue
            try:
                with os.scandir(d) as it:
                    for entry in it:
                        if entry.is_file(follow_symlinks=False) and entry.name.endswith(".pod5"):
                            chem = extract_chemistry(Path(entry.path), "pod5")
                            if chem:
                                result["chemistry"] = chem
                                result["chemistry_classification"] = classify_chemistry(chem)
                            break
            except PermissionError:
                pass
            if result["chemistry"]:
                break
```

**3c.** After fast5 classification (where `sample_file` is a full Path, around line 447), add:
```python
        if sample_file and result["chemistry"] is None:
            chem = extract_chemistry(sample_file, "fast5")
            if chem:
                result["chemistry"] = chem
                result["chemistry_classification"] = classify_chemistry(chem)
```

**3d.** For fast5 in root/numeric subdirs (line 532 area), fix `sample_file` to store full paths. Change:
```python
sample_file = entry.name
```
to:
```python
sample_file = Path(entry.path)
```
And same for the numeric subdir case (line 545 area):
```python
sample_file = sub_entry.name
```
to:
```python
sample_file = Path(sub_entry.path)
```

Then after line 555 (`if sample_file and classification:`), before building `fast5_info`, add:
```python
        if sample_file and classification and result["chemistry"] is None:
            chem = extract_chemistry(sample_file if isinstance(sample_file, Path) else run_path / sample_file, "fast5")
            if chem:
                result["chemistry"] = chem
                result["chemistry_classification"] = classify_chemistry(chem)
```

Also update the `fast5_info["sampled_file"]` to handle Path objects:
```python
"sampled_file": sample_file.name if isinstance(sample_file, Path) else sample_file,
```

**Step 4: Run tests to verify they pass**

Run: `python test_optimizations.py`
Expected: All pass including the new integration test

**Step 5: Commit**

```bash
git add nanopore_format_checker.py test_optimizations.py
git commit -m "Integrate chemistry extraction into analyze_run()"
```

---

### Task 6: Update main() output -- table, verbose, and TSV

**Files:**
- Modify: `nanopore_format_checker.py:925-994` (main output loop)
- Modify: `nanopore_format_checker.py:810-852` (write_stats_tsv)
- Test: `test_optimizations.py`

**Step 1: Write failing test for TSV**

```python
def test_tsv_includes_chemistry_columns(tmp_path: Path) -> None:
    """TSV output should include chemistry columns."""
    all_runs = {
        "20240101_run_chem": {
            "formats": ["pod5"],
            "details": {
                "pod5": {
                    "file_count": 5,
                    "data_size_bytes": 1000,
                    "directories": ["/data/run/pod5"],
                }
            },
            "chemistry": {"flowcell": "FLO-MIN114", "kit": "SQK-LSK114", "sample_rate": 5000},
            "chemistry_classification": {
                "pore": "R10.4.1", "analyte": "dna",
                "dorado_version": ">=1.0", "model_hint": "sup", "note": None,
            },
        }
    }
    tsv_path = str(tmp_path / "stats_chem.tsv")
    write_stats_tsv(all_runs, tsv_path)

    with open(tsv_path) as fh:
        lines = fh.readlines()
    header = lines[0].strip().split("\t")
    assert "flowcell_code" in header
    assert "sequencing_kit" in header
    assert "sample_rate" in header
    assert "pore_type" in header
    assert "dorado_version" in header
    row = lines[1].strip().split("\t")
    fc_idx = header.index("flowcell_code")
    assert row[fc_idx] == "FLO-MIN114"
    pore_idx = header.index("pore_type")
    assert row[pore_idx] == "R10.4.1"
    print("  PASS: TSV includes chemistry columns")
```

Add to test list.

**Step 2: Run test to verify it fails**

Run: `python test_optimizations.py`
Expected: FAIL -- "flowcell_code" not in header

**Step 3: Update write_stats_tsv()**

In `write_stats_tsv()`, update the header and row construction:

```python
header = ["run_name", "format", "file_count", "data_size_bytes",
          "size_estimated", "directories", "notes",
          "flowcell_code", "sequencing_kit", "sample_rate",
          "pore_type", "dorado_version"]
```

After the `notes` variable and before building `row`, add:

```python
chem = info.get("chemistry") or {}
chem_class = info.get("chemistry_classification") or {}
```

Extend the row list with:

```python
row = [
    run_name,
    fmt,
    str(file_count),
    str(data_size),
    str(estimated),
    dirs,
    notes,
    chem.get("flowcell", ""),
    chem.get("kit", ""),
    str(chem.get("sample_rate", "")) if chem.get("sample_rate") else "",
    chem_class.get("pore", ""),
    chem_class.get("dorado_version", ""),
]
```

**Step 4: Update main table header and row format**

In `main()`, update the table header (line ~926):

```python
print(f"{'Run / Archive':<55} {'Format(s)':<25} {'Chemistry':<18} {'Data size':<12} {'Files'}")
print("-" * 120)
```

After building `size_str` and before printing, build chemistry string:

```python
chem_class = result.get("chemistry_classification")
if chem_class and chem_class["pore"] != "unknown":
    chem_str = chem_class["pore"]
    rate = result.get("chemistry", {}).get("sample_rate", 0)
    if rate:
        chem_str += f" {rate // 1000}kHz"
else:
    chem_str = "-"
```

Update the print lines to include chemistry:

For unknown/permission_denied:
```python
print(f"{run_dir.name:<55} {fmt_key:<25} {'-':<18} {size_str:<12} {short_reason}")
```

For normal runs:
```python
print(f"{run_dir.name:<55} {formats_str:<25} {chem_str:<18} {size_str:<12} {count_str}")
```

**Step 5: Update verbose output**

After the existing verbose format/detail block (around line 994), before `print()`, add chemistry verbose output:

```python
            # Show chemistry info in verbose mode
            chem = result.get("chemistry")
            chem_class = result.get("chemistry_classification")
            if chem and chem_class:
                rate = chem.get("sample_rate", 0)
                rate_str = f", {rate} Hz" if rate else ""
                print(f"  chemistry: {chem_class['pore']} ({chem['flowcell']}, {chem['kit']}{rate_str})")
                if chem_class.get("dorado_version"):
                    ver = chem_class["dorado_version"]
                    print(f"  dorado: use version {ver}")
                    if chem_class.get("note"):
                        print(f"          {chem_class['note']}")
                    print(f"  model: dorado basecaller {chem_class['model_hint']} <input>/")
```

Also update the archive section table line (~line 943 area) and the summary
separator to use 120 width.

**Step 6: Run tests to verify they pass**

Run: `python test_optimizations.py`
Expected: All pass

**Step 7: Commit**

```bash
git add nanopore_format_checker.py test_optimizations.py
git commit -m "Add chemistry column to table output, verbose dorado hints, and TSV columns"
```

---

### Task 7: Update CLAUDE.md and final verification

**Files:**
- Modify: `CLAUDE.md`

**Step 1: Update CLAUDE.md**

Add to the key functions list:
- `extract_chemistry()` / `extract_chemistry_fast5()` / `extract_chemistry_pod5()` -- reads flowcell/kit/sample_rate from data files
- `classify_chemistry()` -- maps flowcell code to pore type and dorado version recommendation

Add to conventions:
- Chemistry extraction is best-effort; failures never block format detection
- Flowcell product codes are normalized to uppercase (fast5 stores them lowercase)

**Step 2: Run all tests**

Run: `python test_optimizations.py`
Expected: All pass (should be ~50+ tests)

**Step 3: Test with real data (manual)**

```bash
python nanopore_format_checker.py /path/to/runs --verbose
python nanopore_format_checker.py /path/to/runs -o stats.tsv
```

Verify: chemistry column populated, verbose dorado hints shown, TSV has new columns.

**Step 4: Commit**

```bash
git add CLAUDE.md
git commit -m "Update CLAUDE.md with chemistry detection documentation"
```

---

### Task 8: Push all commits

**Step 1: Push**

```bash
git push
```
