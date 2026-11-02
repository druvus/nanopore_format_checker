#!/usr/bin/env python3
"""Tests for nanopore_converter.py"""

import shutil
import subprocess
import sys
import tarfile as _tarfile
import unittest.mock as mock
import zipfile as _zipfile
from pathlib import Path


def _emulate_tool(cmd):
    """Mimic the side effect of single_to_multi_fast5 (writes into -s DIR)."""
    if cmd and cmd[0] == "single_to_multi_fast5":
        save = Path(cmd[cmd.index("-s") + 1])
        save.mkdir(parents=True, exist_ok=True)
        (save / "batch_0.fast5").write_bytes(b"x" * 10)

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

def test_convert_run_multi_read(tmp_path, monkeypatch):
    """Multi-read fast5 should call pod5 convert directly."""
    from nanopore_converter import convert_run
    calls = []
    def mock_run(cmd, **kwargs):
        calls.append(cmd)
        _emulate_tool(cmd)
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
        _emulate_tool(cmd)
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

def test_process_item_directory(tmp_path, monkeypatch):
    """process_item should call convert_run for a directory."""
    from nanopore_converter import process_item
    calls = []
    monkeypatch.setattr(
        "nanopore_converter.convert_run",
        lambda *a, **kw: calls.append(("convert_run", a, kw)) or True
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
        lambda *a, **kw: calls.append(("convert_run", a, kw)) or True
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
        _emulate_tool(cmd)
        result = mock.Mock()
        result.returncode = 0
        return result
    monkeypatch.setattr(subprocess, "run", mock_run)
    monkeypatch.setattr(shutil, "which", lambda cmd: f"/usr/bin/{cmd}")

    from nanopore_converter import main
    rc = main([str(listf), "--output-dir", str(output)])
    assert rc == 0


# ---------------------------------------------------------------------------
# Regression tests from the code audit
# ---------------------------------------------------------------------------

def test_classify_small_multi_read_by_structure(tmp_path):
    """A <1 MB multi-read fast5 must be detected as multi_read via HDF5."""
    import pytest
    h5py = pytest.importorskip("h5py")
    from nanopore_converter import classify_fast5_dir
    d = tmp_path / "fast5_fail"
    d.mkdir()
    with h5py.File(d / "batch_0.fast5", "w") as f:
        f.attrs["file_type"] = "multi-read"
        f.create_group("read_00000000-0000-0000-0000-000000000000")
    assert (d / "batch_0.fast5").stat().st_size < 1_000_000
    assert classify_fast5_dir(str(d)) == "multi_read_fast5"


def test_single_to_multi_without_output_fails(tmp_path, monkeypatch):
    """If single_to_multi_fast5 writes nothing, the run must fail."""
    import pytest
    from nanopore_converter import convert_run

    def mock_run(cmd, **kwargs):  # exits 0 but writes nothing
        result = mock.Mock()
        result.returncode = 0
        result.stderr = "Cannot input MultiRead files to single_to_multi"
        return result
    monkeypatch.setattr(subprocess, "run", mock_run)

    input_dir = tmp_path / "run1"
    input_dir.mkdir()
    (input_dir / "read.fast5").write_bytes(b"x" * 500)
    with pytest.raises(subprocess.CalledProcessError):
        convert_run(str(input_dir), str(tmp_path / "out"))


def test_archive_output_named_after_archive(tmp_path, monkeypatch):
    """Archives with files at the top level must not produce temp-dir names."""
    from nanopore_converter import process_item
    calls = []
    monkeypatch.setattr(
        "nanopore_converter.convert_run",
        lambda *a, **kw: calls.append(kw) or True,
    )
    archive = tmp_path / "20240101_run.tar.gz"
    src = tmp_path / "a.fast5"
    src.write_bytes(b"x" * 500)
    src2 = tmp_path / "b.fast5"
    src2.write_bytes(b"x" * 500)
    with _tarfile.open(archive, "w:gz") as tar:
        tar.add(src, arcname="a.fast5")
        tar.add(src2, arcname="b.fast5")
    assert process_item(str(archive), str(tmp_path / "out")) is True
    assert calls[0]["run_name"] == "20240101_run"


def test_process_item_no_fast5_is_failure(tmp_path):
    """A directory without fast5 files is reported as a failure."""
    from nanopore_converter import process_item
    d = tmp_path / "empty_run"
    d.mkdir()
    assert process_item(str(d), str(tmp_path / "out")) is False


def test_extract_rejects_path_traversal(tmp_path):
    """Tar members escaping the extraction dir must be rejected."""
    import io
    import pytest
    from nanopore_converter import extract_archive
    if not hasattr(_tarfile, "data_filter"):
        pytest.skip("tarfile extraction filters unavailable")
    archive = tmp_path / "evil.tar.gz"
    with _tarfile.open(archive, "w:gz") as tar:
        info = _tarfile.TarInfo("../escaped.fast5")
        info.size = 3
        tar.addfile(info, io.BytesIO(b"abc"))
    with pytest.raises(_tarfile.TarError):
        extract_archive(str(archive))
    assert not (Path(_tarfile.__file__).parent / "escaped.fast5").exists()
