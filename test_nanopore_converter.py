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
    archive = tmp_path / "evil.tar.gz"
    with _tarfile.open(archive, "w:gz") as tar:
        info = _tarfile.TarInfo("../escaped.fast5")
        info.size = 3
        tar.addfile(info, io.BytesIO(b"abc"))
    with pytest.raises(_tarfile.TarError):
        extract_archive(str(archive))
    assert not (Path(_tarfile.__file__).parent / "escaped.fast5").exists()


# ---------------------------------------------------------------------------
# Re-running into an existing output directory
# ---------------------------------------------------------------------------

def _fake_tools(monkeypatch, calls, fail_pod5=False):
    """Mock subprocess.run: single_to_multi writes a batch, pod5 writes output.pod5."""
    def mock_run(cmd, **kwargs):
        calls.append(cmd)
        _emulate_tool(cmd)
        result = mock.Mock()
        result.returncode = 0
        result.stderr = ""
        if cmd[0] == "pod5":
            if fail_pod5:
                result.returncode = 1
                return result
            out = Path(cmd[cmd.index("--output") + 1])
            # Mimic pod5's refusal to overwrite an existing file
            assert not (out / "output.pod5").exists(), "pod5 would refuse to overwrite"
            (out / "output.pod5").write_bytes(b"new")
        return result
    monkeypatch.setattr(subprocess, "run", mock_run)
    monkeypatch.setattr(shutil, "which", lambda cmd: f"/usr/bin/{cmd}")


def _make_run(tmp_path, name="run1"):
    d = tmp_path / name
    d.mkdir()
    (d / "reads.fast5").write_bytes(b"x" * 2_000_000)  # multi-read by size
    return d


def test_rerun_skips_converted_input(tmp_path, monkeypatch):
    from nanopore_converter import main
    calls = []
    _fake_tools(monkeypatch, calls)
    run = _make_run(tmp_path)
    out = tmp_path / "out"
    assert main([str(run), "--output-dir", str(out)]) == 0
    assert len(calls) == 1
    assert (out / "run1" / "pod5" / "output.pod5").read_bytes() == b"new"
    # Second run: nothing is converted, exit code still 0
    assert main([str(run), "--output-dir", str(out)]) == 0
    assert len(calls) == 1


def test_rerun_skips_archive_without_extracting(tmp_path, monkeypatch):
    from nanopore_converter import main
    calls = []
    _fake_tools(monkeypatch, calls)
    src = _make_run(tmp_path, "20240101_run")
    archive = tmp_path / "20240101_run.tar.gz"
    with _tarfile.open(archive, "w:gz") as tar:
        tar.add(src, arcname="20240101_run")
    out = tmp_path / "out"
    (out / "20240101_run" / "pod5").mkdir(parents=True)
    (out / "20240101_run" / "pod5" / "output.pod5").write_bytes(b"old")
    extracted = []
    monkeypatch.setattr("nanopore_converter.extract_archive",
                        lambda p: extracted.append(p) or str(tmp_path))
    assert main([str(archive), "--output-dir", str(out)]) == 0
    assert extracted == [] and calls == []


def test_force_replaces_existing_output(tmp_path, monkeypatch):
    from nanopore_converter import main
    calls = []
    _fake_tools(monkeypatch, calls)
    run = _make_run(tmp_path)
    pod5_dir = tmp_path / "out" / "run1" / "pod5"
    pod5_dir.mkdir(parents=True)
    (pod5_dir / "output.pod5").write_bytes(b"old")
    assert main([str(run), "--output-dir", str(tmp_path / "out"), "--force"]) == 0
    assert (pod5_dir / "output.pod5").read_bytes() == b"new"
    assert sorted(p.name for p in pod5_dir.parent.iterdir()) == ["pod5"]


def test_force_keeps_old_output_when_conversion_fails(tmp_path, monkeypatch):
    from nanopore_converter import main
    calls = []
    _fake_tools(monkeypatch, calls, fail_pod5=True)
    run = _make_run(tmp_path)
    pod5_dir = tmp_path / "out" / "run1" / "pod5"
    pod5_dir.mkdir(parents=True)
    (pod5_dir / "output.pod5").write_bytes(b"old")
    assert main([str(run), "--output-dir", str(tmp_path / "out"), "--force"]) == 1
    assert (pod5_dir / "output.pod5").read_bytes() == b"old"
    assert sorted(p.name for p in pod5_dir.parent.iterdir()) == ["pod5"]


def test_failed_conversion_leaves_no_pod5_dir(tmp_path, monkeypatch):
    """A failure must not create pod5/, or the next run would skip it as done."""
    from nanopore_converter import is_converted, main
    calls = []
    _fake_tools(monkeypatch, calls, fail_pod5=True)
    run = _make_run(tmp_path)
    out = tmp_path / "out"
    assert main([str(run), "--output-dir", str(out)]) == 1
    assert not (out / "run1" / "pod5").exists()
    assert not is_converted(str(out), "run1")


def test_leftovers_from_interrupted_run_are_cleared(tmp_path, monkeypatch):
    from nanopore_converter import convert_run
    calls = []
    _fake_tools(monkeypatch, calls)
    run = tmp_path / "run1"
    run.mkdir()
    (run / "read.fast5").write_bytes(b"x" * 500)  # single-read by size
    out = tmp_path / "out"
    (out / "run1" / ".pod5_partial").mkdir(parents=True)
    (out / "run1" / ".pod5_partial" / "output.pod5").write_bytes(b"partial")
    (out / "run1" / "_multi_fast5_tmp").mkdir()
    (out / "run1" / "_multi_fast5_tmp" / "stale.fast5").write_bytes(b"stale")
    assert convert_run(str(run), str(out)) is True
    assert (out / "run1" / "pod5" / "output.pod5").read_bytes() == b"new"
    assert sorted(p.name for p in (out / "run1").iterdir()) == ["pod5"]


def test_duplicate_output_names_in_one_invocation_fail(tmp_path, monkeypatch):
    from nanopore_converter import main
    calls = []
    _fake_tools(monkeypatch, calls)
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    a = _make_run(tmp_path / "a")  # a/run1
    b = _make_run(tmp_path / "b")  # b/run1 -> same output name
    listf = tmp_path / "inputs.txt"
    listf.write_text(f"{a}\n{b}\n")
    assert main([str(listf), "--output-dir", str(tmp_path / "out")]) == 1
    assert len(calls) == 1  # only the first one was converted


# ---------------------------------------------------------------------------
# Filling missing pod5 metadata (dorado automatic model selection needs it)
# ---------------------------------------------------------------------------

def _write_pod5(path, flowcell="", kit="", sample_rate=4000, n_reads=3):
    import uuid
    from datetime import datetime, timezone
    import numpy as np
    import pod5
    now = datetime.now(timezone.utc)
    run_info = pod5.RunInfo(
        acquisition_id="acq1", acquisition_start_time=now, adc_max=4095, adc_min=-4096,
        context_tags={}, experiment_name="", flow_cell_id="FAK1",
        flow_cell_product_code=flowcell, protocol_name="", protocol_run_id="",
        protocol_start_time=now, sample_id="", sample_rate=sample_rate,
        sequencing_kit=kit, sequencer_position="", sequencer_position_type="",
        software="", system_name="", system_type="", tracking_id={},
    )
    with pod5.Writer(path) as writer:
        for _ in range(n_reads):
            writer.add_read(pod5.Read(
                read_id=uuid.uuid4(), pore=pod5.Pore(1, 1, "pore"),
                calibration=pod5.Calibration(0.0, 1.0), read_number=1,
                start_sample=0, median_before=200.0,
                end_reason=pod5.EndReason.from_reason_with_default_forced(
                    pod5.EndReasonEnum.SIGNAL_POSITIVE),
                run_info=run_info, signal=np.zeros(100, dtype=np.int16),
            ))


def test_fill_missing_pod5_metadata_only_fills_empty_fields(tmp_path):
    import pytest
    pod5 = pytest.importorskip("pod5")
    from nanopore_converter import fill_missing_pod5_metadata
    d = tmp_path / "pod5"
    d.mkdir()
    _write_pod5(d / "output.pod5", flowcell="", kit="sqk-lsk109")
    chem = {"flowcell": "FLO-MIN106", "kit": "SQK-OTHER", "sample_rate": 4000}
    assert fill_missing_pod5_metadata(str(d), chem) == 1
    with pod5.Reader(d / "output.pod5") as reader:
        infos = [(r.run_info.flow_cell_product_code, r.run_info.sequencing_kit)
                 for r in reader.reads()]
    # flowcell filled; kit kept exactly as pod5 wrote it
    assert infos == [("FLO-MIN106", "sqk-lsk109")] * 3
    assert not list(d.glob("*.tmp"))
    # Nothing missing any more -> file is not rewritten again
    assert fill_missing_pod5_metadata(str(d), chem) == 0


def test_fill_missing_pod5_metadata_never_raises(tmp_path):
    from nanopore_converter import fill_missing_pod5_metadata
    d = tmp_path / "pod5"
    d.mkdir()
    (d / "output.pod5").write_bytes(b"not a pod5 file")
    assert fill_missing_pod5_metadata(str(d), {"flowcell": "FLO-MIN106"}) == 0
    assert (d / "output.pod5").read_bytes() == b"not a pod5 file"
    assert not list(d.glob("*.tmp"))
