# Nanopore Converter Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a standalone Python CLI script that converts single_read_fast5 and multi_read_fast5 files to pod5 format, accepting folders, compressed archives, or input list files.

**Architecture:** Single-file sequential pipeline. Auto-detects input type (directory, archive, list file), classifies fast5 format using a 1 MB size heuristic, and shells out to `single_to_multi_fast5` and `pod5 convert fast5` as subprocesses. Archives are extracted to temp directories and cleaned up after conversion.

**Tech Stack:** Python 3.10+ standard library (subprocess, tempfile, tarfile, zipfile, shutil, argparse, logging, pathlib). External CLI tools: `single_to_multi_fast5` (from ont_fast5_api), `pod5` (from pod5 package).

---

### Task 1: Script skeleton and argument parsing

**Files:**
- Create: `nanopore_converter.py`
- Test: `test_nanopore_converter.py`

**Step 1: Write the failing test for argument parsing**

```python
#!/usr/bin/env python3
"""Tests for nanopore_converter.py"""

import subprocess
import sys

def test_help_flag():
    """Script should print help and exit 0."""
    result = subprocess.run(
        [sys.executable, "nanopore_converter.py", "--help"],
        capture_output=True, text=True
    )
    assert result.returncode == 0
    assert "--output-dir" in result.stdout
    assert "--threads" in result.stdout
    assert "--verbose" in result.stdout

def test_missing_output_dir():
    """Script should fail with exit 2 when --output-dir is missing."""
    result = subprocess.run(
        [sys.executable, "nanopore_converter.py", "/tmp/fake"],
        capture_output=True, text=True
    )
    assert result.returncode == 2
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest test_nanopore_converter.py::test_help_flag -v`
Expected: FAIL (file does not exist)

**Step 3: Write minimal implementation**

```python
#!/usr/bin/env python3
"""
Nanopore Format Converter
=========================
Converts single_read_fast5 and multi_read_fast5 to pod5 format.

Accepts a folder, compressed archive (.tar.gz, .tar.bz2, .zip), or a text
file listing paths (one per line, mixed folders and archives).

Dependencies:
    pip install ont_fast5_api pod5

Usage:
    python nanopore_converter.py /path/to/input --output-dir /path/to/output
    python nanopore_converter.py input_list.txt --output-dir /path/to/output
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path

logger = logging.getLogger(__name__)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert nanopore fast5 files to pod5 format."
    )
    parser.add_argument(
        "input",
        help="Folder, compressed archive, or text file listing paths (one per line).",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Base output directory for converted pod5 files.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Thread count for conversion tools (default: 4).",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable detailed logging.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    logger.info("Nanopore converter starting.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest test_nanopore_converter.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add nanopore_converter.py test_nanopore_converter.py
git commit -m "feat: add converter script skeleton with argument parsing"
```

---

### Task 2: Tool availability check

**Files:**
- Modify: `nanopore_converter.py`
- Test: `test_nanopore_converter.py`

**Step 1: Write the failing test**

```python
import unittest.mock as mock

def test_check_tools_available(tmp_path):
    """Should succeed when both tools are on PATH."""
    from nanopore_converter import check_tools
    # Mock shutil.which to return a path for both tools
    with mock.patch("nanopore_converter.shutil.which") as mock_which:
        mock_which.side_effect = lambda cmd: f"/usr/bin/{cmd}"
        check_tools()  # Should not raise

def test_check_tools_missing():
    """Should raise SystemExit when tools are missing."""
    from nanopore_converter import check_tools
    import pytest
    with mock.patch("nanopore_converter.shutil.which", return_value=None):
        with pytest.raises(SystemExit):
            check_tools()
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest test_nanopore_converter.py::test_check_tools_available -v`
Expected: FAIL (check_tools not defined)

**Step 3: Write minimal implementation**

Add to `nanopore_converter.py`:

```python
def check_tools():
    """Verify that required CLI tools are on PATH."""
    missing = []
    for tool in ("single_to_multi_fast5", "pod5"):
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        logger.error("Required tools not found on PATH: %s", ", ".join(missing))
        logger.error("Install with: pip install ont_fast5_api pod5")
        sys.exit(2)
```

Call `check_tools()` at the start of `main()`.

**Step 4: Run tests to verify they pass**

Run: `python -m pytest test_nanopore_converter.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add nanopore_converter.py test_nanopore_converter.py
git commit -m "feat: add tool availability check for single_to_multi_fast5 and pod5"
```

---

### Task 3: Input type detection

**Files:**
- Modify: `nanopore_converter.py`
- Test: `test_nanopore_converter.py`

**Step 1: Write the failing tests**

```python
def test_detect_input_directory(tmp_path):
    """Directory input should be detected."""
    from nanopore_converter import detect_input_type
    d = tmp_path / "run1"
    d.mkdir()
    assert detect_input_type(str(d)) == "directory"

def test_detect_input_archive(tmp_path):
    """Archive files should be detected."""
    from nanopore_converter import detect_input_type
    for ext in (".tar.gz", ".tar.bz2", ".zip"):
        f = tmp_path / f"run1{ext}"
        f.touch()
        assert detect_input_type(str(f)) == "archive", f"Failed for {ext}"

def test_detect_input_listfile(tmp_path):
    """Text files should be detected as list files."""
    from nanopore_converter import detect_input_type
    for ext in (".txt", ".csv", ".list"):
        f = tmp_path / f"inputs{ext}"
        f.write_text("/some/path\n")
        assert detect_input_type(str(f)) == "listfile", f"Failed for {ext}"

def test_detect_input_nonexistent():
    """Non-existent path should return None."""
    from nanopore_converter import detect_input_type
    assert detect_input_type("/nonexistent/path") is None
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest test_nanopore_converter.py::test_detect_input_directory -v`
Expected: FAIL

**Step 3: Write minimal implementation**

```python
ARCHIVE_EXTENSIONS = {".tar.gz", ".tar.bz2", ".tgz", ".zip"}
LIST_EXTENSIONS = {".txt", ".csv", ".list"}


def detect_input_type(path_str: str) -> str | None:
    """Determine if input is a directory, archive, or list file.

    Returns 'directory', 'archive', 'listfile', or None.
    """
    p = Path(path_str)
    if not p.exists():
        return None
    if p.is_dir():
        return "directory"
    if p.is_file():
        # Check for double extensions like .tar.gz
        name = p.name.lower()
        for ext in ARCHIVE_EXTENSIONS:
            if name.endswith(ext):
                return "archive"
        if p.suffix.lower() in LIST_EXTENSIONS:
            return "listfile"
    return None
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest test_nanopore_converter.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add nanopore_converter.py test_nanopore_converter.py
git commit -m "feat: add input type detection (directory, archive, listfile)"
```

---

### Task 4: Archive extraction

**Files:**
- Modify: `nanopore_converter.py`
- Test: `test_nanopore_converter.py`

**Step 1: Write the failing tests**

```python
import tarfile as _tarfile
import zipfile as _zipfile

def test_extract_tar_gz(tmp_path):
    """Should extract .tar.gz to temp dir and return the path."""
    from nanopore_converter import extract_archive
    # Create a tar.gz with a file inside
    src = tmp_path / "data"
    src.mkdir()
    (src / "test.fast5").write_bytes(b"fake")
    archive = tmp_path / "data.tar.gz"
    with _tarfile.open(archive, "w:gz") as tar:
        tar.add(src / "test.fast5", arcname="data/test.fast5")
    result = extract_archive(str(archive))
    assert Path(result).is_dir()
    # Find the fast5 inside
    fast5_files = list(Path(result).rglob("*.fast5"))
    assert len(fast5_files) == 1
    shutil.rmtree(result)

def test_extract_zip(tmp_path):
    """Should extract .zip to temp dir and return the path."""
    from nanopore_converter import extract_archive
    archive = tmp_path / "data.zip"
    with _zipfile.ZipFile(archive, "w") as zf:
        zf.writestr("data/test.fast5", b"fake")
    result = extract_archive(str(archive))
    assert Path(result).is_dir()
    fast5_files = list(Path(result).rglob("*.fast5"))
    assert len(fast5_files) == 1
    shutil.rmtree(result)
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest test_nanopore_converter.py::test_extract_tar_gz -v`
Expected: FAIL

**Step 3: Write minimal implementation**

```python
def extract_archive(archive_path: str) -> str:
    """Extract a compressed archive to a temporary directory.

    Returns the path to the temporary directory. Caller is responsible for
    cleanup.
    """
    tmp_dir = tempfile.mkdtemp(prefix="nanopore_conv_")
    p = Path(archive_path)
    name = p.name.lower()

    logger.info("Extracting %s to %s", archive_path, tmp_dir)

    if name.endswith((".tar.gz", ".tgz", ".tar.bz2")):
        mode = "r:gz" if name.endswith((".tar.gz", ".tgz")) else "r:bz2"
        with tarfile.open(archive_path, mode) as tar:
            tar.extractall(tmp_dir)
    elif name.endswith(".zip"):
        with zipfile.ZipFile(archive_path, "r") as zf:
            zf.extractall(tmp_dir)
    else:
        shutil.rmtree(tmp_dir)
        raise ValueError(f"Unsupported archive format: {archive_path}")

    return tmp_dir
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest test_nanopore_converter.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add nanopore_converter.py test_nanopore_converter.py
git commit -m "feat: add archive extraction for tar.gz, tar.bz2, and zip"
```

---

### Task 5: Fast5 format classification

**Files:**
- Modify: `nanopore_converter.py`
- Test: `test_nanopore_converter.py`

**Step 1: Write the failing tests**

```python
def test_classify_fast5_single(tmp_path):
    """Small fast5 files should be classified as single_read."""
    from nanopore_converter import classify_fast5_dir
    d = tmp_path / "fast5"
    d.mkdir()
    # Single-read fast5: < 1MB
    (d / "read1.fast5").write_bytes(b"x" * 500)
    assert classify_fast5_dir(str(d)) == "single_read_fast5"

def test_classify_fast5_multi(tmp_path):
    """Large fast5 files should be classified as multi_read."""
    from nanopore_converter import classify_fast5_dir
    d = tmp_path / "fast5"
    d.mkdir()
    # Multi-read fast5: >= 1MB
    (d / "batch0.fast5").write_bytes(b"x" * 2_000_000)
    assert classify_fast5_dir(str(d)) == "multi_read_fast5"

def test_classify_fast5_empty(tmp_path):
    """Empty directory should return None."""
    from nanopore_converter import classify_fast5_dir
    d = tmp_path / "empty"
    d.mkdir()
    assert classify_fast5_dir(str(d)) is None

def test_classify_fast5_nested(tmp_path):
    """Should find fast5 in nested subdirectories."""
    from nanopore_converter import classify_fast5_dir
    d = tmp_path / "run" / "fast5_pass" / "barcode01"
    d.mkdir(parents=True)
    (d / "read1.fast5").write_bytes(b"x" * 500)
    assert classify_fast5_dir(str(tmp_path / "run")) == "single_read_fast5"
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest test_nanopore_converter.py::test_classify_fast5_single -v`
Expected: FAIL

**Step 3: Write minimal implementation**

```python
def classify_fast5_dir(dir_path: str) -> str | None:
    """Classify a directory of fast5 files as single_read or multi_read.

    Samples one .fast5 file and uses the 1 MB size threshold:
    - < 1 MB: single_read_fast5
    - >= 1 MB: multi_read_fast5

    Returns None if no .fast5 files are found.
    """
    root = Path(dir_path)
    for dirpath, _, filenames in os.walk(root):
        for fname in filenames:
            if fname.lower().endswith(".fast5"):
                fpath = Path(dirpath) / fname
                try:
                    size = fpath.stat().st_size
                except OSError:
                    continue
                if size < 1_000_000:
                    return "single_read_fast5"
                else:
                    return "multi_read_fast5"
    return None
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest test_nanopore_converter.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add nanopore_converter.py test_nanopore_converter.py
git commit -m "feat: add fast5 format classification using size heuristic"
```

---

### Task 6: Single-run conversion logic

**Files:**
- Modify: `nanopore_converter.py`
- Test: `test_nanopore_converter.py`

**Step 1: Write the failing tests**

```python
def test_convert_run_multi_read(tmp_path, monkeypatch):
    """Multi-read fast5 should call pod5 convert directly."""
    from nanopore_converter import convert_run
    calls = []
    def mock_run(cmd, **kwargs):
        calls.append(cmd)
        result = mock.Mock()
        result.returncode = 0
        return result
    monkeypatch.setattr(subprocess, "run", mock_run)

    input_dir = tmp_path / "run1"
    input_dir.mkdir()
    (input_dir / "batch.fast5").write_bytes(b"x" * 2_000_000)

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    convert_run(str(input_dir), str(output_dir), threads=4)

    # Should have called pod5 convert fast5 once
    assert len(calls) == 1
    assert "pod5" in calls[0][0]

def test_convert_run_single_read(tmp_path, monkeypatch):
    """Single-read fast5 should call single_to_multi_fast5 then pod5 convert."""
    from nanopore_converter import convert_run
    calls = []
    def mock_run(cmd, **kwargs):
        calls.append(cmd)
        result = mock.Mock()
        result.returncode = 0
        return result
    monkeypatch.setattr(subprocess, "run", mock_run)

    input_dir = tmp_path / "run1"
    input_dir.mkdir()
    (input_dir / "read.fast5").write_bytes(b"x" * 500)

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    convert_run(str(input_dir), str(output_dir), threads=4)

    # Should have called single_to_multi_fast5, then pod5 convert
    assert len(calls) == 2
    assert "single_to_multi_fast5" in calls[0][0]
    assert "pod5" in calls[1][0]
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest test_nanopore_converter.py::test_convert_run_multi_read -v`
Expected: FAIL

**Step 3: Write minimal implementation**

```python
def convert_run(input_dir: str, output_base: str, threads: int = 4) -> None:
    """Convert a single run directory from fast5 to pod5.

    For single_read_fast5: two-step process via intermediate multi_read_fast5.
    For multi_read_fast5: direct pod5 conversion.

    Raises subprocess.CalledProcessError on conversion failure.
    """
    fmt = classify_fast5_dir(input_dir)
    if fmt is None:
        logger.warning("No .fast5 files found in %s, skipping.", input_dir)
        return

    run_name = Path(input_dir).name
    pod5_out = os.path.join(output_base, run_name, "pod5")
    os.makedirs(pod5_out, exist_ok=True)

    if fmt == "single_read_fast5":
        # Step 1: single -> multi
        multi_tmp = os.path.join(output_base, run_name, "_multi_fast5_tmp")
        os.makedirs(multi_tmp, exist_ok=True)
        logger.info("[%s] Step 1/2: single_to_multi_fast5", run_name)
        result = subprocess.run(
            [
                "single_to_multi_fast5",
                "-i", input_dir,
                "-s", multi_tmp,
                "-t", str(threads),
                "--recursive",
            ],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            logger.error("single_to_multi_fast5 failed:\n%s", result.stderr)
            raise subprocess.CalledProcessError(
                result.returncode, "single_to_multi_fast5"
            )

        # Step 2: multi -> pod5
        logger.info("[%s] Step 2/2: pod5 convert fast5", run_name)
        result = subprocess.run(
            [
                "pod5", "convert", "fast5",
                multi_tmp + "/",
                "--output", pod5_out + "/",
                "--threads", str(threads),
                "--recursive",
            ],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            logger.error("pod5 convert failed:\n%s", result.stderr)
            raise subprocess.CalledProcessError(
                result.returncode, "pod5 convert fast5"
            )

        # Cleanup intermediate files
        logger.info("[%s] Cleaning up intermediate multi_fast5.", run_name)
        shutil.rmtree(multi_tmp, ignore_errors=True)

    elif fmt == "multi_read_fast5":
        logger.info("[%s] Converting multi_read_fast5 -> pod5", run_name)
        result = subprocess.run(
            [
                "pod5", "convert", "fast5",
                input_dir + "/",
                "--output", pod5_out + "/",
                "--threads", str(threads),
                "--recursive",
            ],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            logger.error("pod5 convert failed:\n%s", result.stderr)
            raise subprocess.CalledProcessError(
                result.returncode, "pod5 convert fast5"
            )

    logger.info("[%s] Conversion complete -> %s", run_name, pod5_out)
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest test_nanopore_converter.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add nanopore_converter.py test_nanopore_converter.py
git commit -m "feat: add single-run conversion logic (single/multi fast5 to pod5)"
```

---

### Task 7: List file parsing and batch processing

**Files:**
- Modify: `nanopore_converter.py`
- Test: `test_nanopore_converter.py`

**Step 1: Write the failing tests**

```python
def test_parse_list_file(tmp_path):
    """Should read paths from a list file, skipping blanks and comments."""
    from nanopore_converter import parse_list_file
    listf = tmp_path / "inputs.txt"
    listf.write_text("/path/to/run1\n\n# comment\n/path/to/run2.tar.gz\n")
    paths = parse_list_file(str(listf))
    assert paths == ["/path/to/run1", "/path/to/run2.tar.gz"]

def test_resolve_inputs_directory(tmp_path):
    """A directory input should resolve to a single-item list."""
    from nanopore_converter import resolve_inputs
    d = tmp_path / "run1"
    d.mkdir()
    inputs = resolve_inputs(str(d))
    assert inputs == [str(d)]

def test_resolve_inputs_listfile(tmp_path):
    """A list file should resolve to multiple paths."""
    from nanopore_converter import resolve_inputs
    d1 = tmp_path / "run1"
    d1.mkdir()
    d2 = tmp_path / "run2"
    d2.mkdir()
    listf = tmp_path / "inputs.txt"
    listf.write_text(f"{d1}\n{d2}\n")
    inputs = resolve_inputs(str(listf))
    assert len(inputs) == 2
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest test_nanopore_converter.py::test_parse_list_file -v`
Expected: FAIL

**Step 3: Write minimal implementation**

```python
def parse_list_file(path: str) -> list[str]:
    """Read a text file of paths, one per line.

    Blank lines and lines starting with # are skipped.
    """
    paths = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                paths.append(line)
    return paths


def resolve_inputs(input_path: str) -> list[str]:
    """Resolve the input argument to a list of paths to process.

    Handles directories, archives (returned as-is), and list files.
    """
    input_type = detect_input_type(input_path)
    if input_type == "directory" or input_type == "archive":
        return [input_path]
    elif input_type == "listfile":
        return parse_list_file(input_path)
    else:
        logger.error("Input not found or unrecognized: %s", input_path)
        return []
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest test_nanopore_converter.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add nanopore_converter.py test_nanopore_converter.py
git commit -m "feat: add list file parsing and input resolution"
```

---

### Task 8: Main pipeline orchestration

**Files:**
- Modify: `nanopore_converter.py`
- Test: `test_nanopore_converter.py`

**Step 1: Write the failing tests**

```python
def test_process_item_directory(tmp_path, monkeypatch):
    """process_item should call convert_run for a directory."""
    from nanopore_converter import process_item
    calls = []
    monkeypatch.setattr(
        "nanopore_converter.convert_run",
        lambda *a, **kw: calls.append(("convert_run", a, kw))
    )
    d = tmp_path / "run1"
    d.mkdir()
    (d / "test.fast5").write_bytes(b"x" * 500)
    result = process_item(str(d), str(tmp_path / "out"), threads=4)
    assert result is True
    assert len(calls) == 1

def test_process_item_archive(tmp_path, monkeypatch):
    """process_item should extract archive, convert, then cleanup."""
    from nanopore_converter import process_item
    calls = []
    monkeypatch.setattr(
        "nanopore_converter.convert_run",
        lambda *a, **kw: calls.append(("convert_run", a, kw))
    )
    # Create a tar.gz with a fake fast5
    src = tmp_path / "data"
    src.mkdir()
    (src / "test.fast5").write_bytes(b"x" * 500)
    archive = tmp_path / "data.tar.gz"
    with _tarfile.open(archive, "w:gz") as tar:
        tar.add(src / "test.fast5", arcname="data/test.fast5")
    result = process_item(str(archive), str(tmp_path / "out"), threads=4)
    assert result is True
    assert len(calls) == 1
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest test_nanopore_converter.py::test_process_item_directory -v`
Expected: FAIL

**Step 3: Write minimal implementation**

```python
def get_output_name(input_path: str) -> str:
    """Derive an output folder name from the input path.

    Strips archive extensions to get a clean directory name.
    """
    name = Path(input_path).name
    # Strip known archive suffixes
    for ext in (".tar.gz", ".tar.bz2", ".tgz", ".zip"):
        if name.lower().endswith(ext):
            name = name[: -len(ext)]
            break
    return name


def process_item(item_path: str, output_dir: str, threads: int = 4) -> bool:
    """Process a single input item (directory or archive).

    Returns True on success, False on failure.
    """
    item_type = detect_input_type(item_path)
    tmp_dir = None

    try:
        if item_type == "archive":
            tmp_dir = extract_archive(item_path)
            # The extracted archive may contain a single top-level folder
            # or files directly. Find the directory with fast5 files.
            work_dir = tmp_dir
            # If there is exactly one subdirectory, use it
            entries = list(Path(tmp_dir).iterdir())
            if len(entries) == 1 and entries[0].is_dir():
                work_dir = str(entries[0])
        elif item_type == "directory":
            work_dir = item_path
        else:
            logger.error("Cannot process: %s (not a directory or archive)", item_path)
            return False

        convert_run(work_dir, output_dir, threads=threads)
        return True

    except subprocess.CalledProcessError as e:
        logger.error("Conversion failed for %s: %s", item_path, e)
        return False
    except Exception as e:
        logger.error("Unexpected error processing %s: %s", item_path, e)
        return False
    finally:
        if tmp_dir and os.path.isdir(tmp_dir):
            logger.info("Cleaning up temp directory: %s", tmp_dir)
            shutil.rmtree(tmp_dir, ignore_errors=True)


def main(argv=None):
    args = parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    check_tools()

    inputs = resolve_inputs(args.input)
    if not inputs:
        logger.error("No valid inputs found.")
        sys.exit(2)

    os.makedirs(args.output_dir, exist_ok=True)

    succeeded = 0
    failed = 0
    for item in inputs:
        logger.info("Processing: %s", item)
        if process_item(item, args.output_dir, threads=args.threads):
            succeeded += 1
        else:
            failed += 1

    logger.info("Done. %d succeeded, %d failed.", succeeded, failed)
    if failed > 0:
        return 1
    return 0
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest test_nanopore_converter.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add nanopore_converter.py test_nanopore_converter.py
git commit -m "feat: add main pipeline orchestration with archive handling"
```

---

### Task 9: Integration test and final polish

**Files:**
- Modify: `test_nanopore_converter.py`

**Step 1: Write integration test**

```python
def test_full_pipeline_dry_run(tmp_path, monkeypatch):
    """End-to-end test with mocked subprocess calls."""
    # Create a list file with a directory and an archive
    run_dir = tmp_path / "run1"
    run_dir.mkdir()
    (run_dir / "reads.fast5").write_bytes(b"x" * 500)

    src2 = tmp_path / "run2"
    src2.mkdir()
    (src2 / "batch.fast5").write_bytes(b"x" * 2_000_000)
    archive = tmp_path / "run2.tar.gz"
    with _tarfile.open(archive, "w:gz") as tar:
        tar.add(src2, arcname="run2")

    listf = tmp_path / "inputs.txt"
    listf.write_text(f"{run_dir}\n{archive}\n")

    output = tmp_path / "output"

    # Mock subprocess.run and shutil.which
    def mock_run(cmd, **kwargs):
        result = mock.Mock()
        result.returncode = 0
        return result
    monkeypatch.setattr(subprocess, "run", mock_run)
    monkeypatch.setattr(shutil, "which", lambda cmd: f"/usr/bin/{cmd}")

    from nanopore_converter import main
    rc = main([str(listf), "--output-dir", str(output)])
    assert rc == 0
```

**Step 2: Run tests**

Run: `python -m pytest test_nanopore_converter.py -v`
Expected: PASS

**Step 3: Commit**

```bash
git add test_nanopore_converter.py
git commit -m "test: add end-to-end integration test"
```

---

## Team Assignment

The 4 team members should work on these tasks:

| Member | Role | Tasks |
|--------|------|-------|
| **architect** | Script skeleton, arg parsing, main orchestration | Task 1, Task 8 |
| **detector** | Input detection, fast5 classification | Task 3, Task 5 |
| **converter** | Tool check, archive extraction, conversion logic | Task 2, Task 4, Task 6 |
| **tester** | List file parsing, integration tests, code review | Task 7, Task 9 |

Dependencies:
- Task 2 (tool check) blocks Task 6 (conversion logic)
- Tasks 3, 5 block Task 8 (main pipeline needs detect + classify)
- Tasks 6, 7 block Task 8 (main pipeline needs convert_run + resolve_inputs)
- Tasks 1-8 all block Task 9 (integration test)
