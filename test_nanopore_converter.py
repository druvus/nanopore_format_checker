#!/usr/bin/env python3
"""Tests for nanopore_converter.py"""

import shutil
import signal
import subprocess
import sys
import tarfile as _tarfile
import unittest.mock as mock
import zipfile as _zipfile
from pathlib import Path

import pytest

import nanopore_converter


# ---------------------------------------------------------------------------
# Minimal but structurally real fast5 files, and fake conversion tools that
# keep honest read counts (so read-count verification is exercised).
# ---------------------------------------------------------------------------

def _single_fast5(path, read_id=None, run_id="run_a"):
    import uuid
    import h5py
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        read = f.create_group("Raw/Reads/Read_1")
        read.attrs["read_id"] = str(read_id or uuid.uuid4())
        read.create_dataset("Signal", data=[1, 2, 3])
        f.create_group("UniqueGlobalKey/tracking_id").attrs["run_id"] = run_id
        f.create_group("UniqueGlobalKey/channel_id")
    return path


def _multi_fast5(path, n_reads=3, read_ids=None, run_id="run_a"):
    import uuid
    import h5py
    path.parent.mkdir(parents=True, exist_ok=True)
    ids = read_ids or [str(uuid.uuid4()) for _ in range(n_reads)]
    with h5py.File(path, "w") as f:
        f.attrs["file_type"] = "multi-read"
        for rid in ids:
            f.create_group(f"read_{rid}/tracking_id").attrs["run_id"] = run_id
            f.create_dataset(f"read_{rid}/Raw/Signal", data=[1, 2, 3])
    return path


def _read_ids_in(paths):
    """Read IDs in all fast5 files under the given dirs (follows symlinks)."""
    import os
    import h5py
    ids = []
    for p in paths:
        for root, _dirs, files in os.walk(p, followlinks=True):
            for name in sorted(files):
                if not name.endswith(".fast5"):
                    continue
                with h5py.File(Path(root) / name, "r") as h:
                    keys = [k[5:] for k in h if k.startswith("read_")]
                    ids.extend(keys or [h["Raw/Reads/Read_1"].attrs["read_id"]])
    return [i.decode() if isinstance(i, bytes) else str(i) for i in ids]


def _reads_in(paths):
    return len(_read_ids_in(paths))


class FakeTools:
    """Stand-in for nanopore_converter._run.

    single_to_multi_fast5 writes one multi-read file holding the input
    files' reads; pod5 convert writes a real pod5 holding exactly the
    input reads (minus `drop`, or with `duplicate` extra copies), refusing
    to overwrite like the real tool.
    """

    def __init__(self, monkeypatch, fail_pod5=False, drop=0, duplicate=0):
        self.calls = []
        self.fail_pod5 = fail_pod5
        self.drop = drop
        self.duplicate = duplicate
        monkeypatch.setattr(nanopore_converter, "_run", self)
        monkeypatch.setattr(shutil, "which", lambda cmd: f"/usr/bin/{cmd}")

    def __call__(self, cmd, label):
        self.calls.append(cmd)
        if cmd[0] == "single_to_multi_fast5":
            src = cmd[cmd.index("-i") + 1]
            save = Path(cmd[cmd.index("-s") + 1])
            _multi_fast5(save / "batch_0.fast5", read_ids=_read_ids_in([src]))
        elif cmd[0] == "pod5":
            if self.fail_pod5:
                raise subprocess.CalledProcessError(1, label)
            out_idx = cmd.index("--output")
            inputs = cmd[3:out_idx]
            out = Path(cmd[out_idx + 1])
            assert not (out / "output.pod5").exists(), "pod5 would refuse to overwrite"
            ids = _read_ids_in(inputs)
            ids = ids[: len(ids) - self.drop] + ids[: self.duplicate]
            _write_pod5(out / "output.pod5", read_ids=ids)
        return ""


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
    from nanopore_converter import classify_fast5_dir
    _single_fast5(tmp_path / "fast5" / "read1.fast5")
    assert classify_fast5_dir(str(tmp_path / "fast5")) == "single_read_fast5"


def test_classify_fast5_multi(tmp_path):
    from nanopore_converter import classify_fast5_dir
    _multi_fast5(tmp_path / "fast5" / "batch0.fast5")
    assert classify_fast5_dir(str(tmp_path / "fast5")) == "multi_read_fast5"


def test_classify_fast5_mixed(tmp_path):
    from nanopore_converter import classify_fast5_dir
    _single_fast5(tmp_path / "run" / "fast5_pass" / "r.fast5")
    _multi_fast5(tmp_path / "run" / "fast5_fail" / "b.fast5")
    assert classify_fast5_dir(str(tmp_path / "run")) == "mixed_fast5"


def test_classify_fast5_empty(tmp_path):
    from nanopore_converter import classify_fast5_dir
    (tmp_path / "empty").mkdir()
    assert classify_fast5_dir(str(tmp_path / "empty")) is None


def test_classify_fast5_nested(tmp_path):
    from nanopore_converter import classify_fast5_dir
    _single_fast5(tmp_path / "run" / "fast5_pass" / "barcode01" / "0" / "read1.fast5")
    assert classify_fast5_dir(str(tmp_path / "run")) == "single_read_fast5"


def test_scan_counts_reads_and_flags_unreadable(tmp_path):
    from nanopore_converter import scan_fast5_dir
    run = tmp_path / "run"
    _single_fast5(run / "a" / "r1.fast5")
    _single_fast5(run / "a" / "r2.fast5")
    _multi_fast5(run / "b" / "m.fast5", 7)
    (run / "b" / "truncated.fast5").write_bytes(b"\x89HDF\r\n\x1a\n" + b"x" * 100)
    scan = scan_fast5_dir(str(run))
    assert len(scan.single) == 2 and len(scan.multi) == 1
    assert [p.name for p in scan.unreadable] == ["truncated.fast5"]
    assert scan.expected_reads == 9


def test_convert_run_multi_read(tmp_path, monkeypatch):
    """Multi-read fast5 should call pod5 convert directly."""
    from nanopore_converter import convert_run
    tools = FakeTools(monkeypatch)
    _multi_fast5(tmp_path / "run1" / "batch.fast5")
    assert convert_run(str(tmp_path / "run1"), str(tmp_path / "output")) is True
    assert [c[0] for c in tools.calls] == ["pod5"]


def test_convert_run_single_read(tmp_path, monkeypatch):
    """Single-read fast5 should call single_to_multi_fast5 then pod5 convert."""
    from nanopore_converter import convert_run
    tools = FakeTools(monkeypatch)
    _single_fast5(tmp_path / "run1" / "read.fast5")
    assert convert_run(str(tmp_path / "run1"), str(tmp_path / "output")) is True
    assert [c[0] for c in tools.calls] == ["single_to_multi_fast5", "pod5"]


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
    assert process_item(str(d), str(tmp_path / "out"), threads=4) is True
    assert len(calls) == 1


def test_process_item_archive(tmp_path, monkeypatch):
    """process_item should extract archive, convert, then cleanup."""
    from nanopore_converter import process_item
    calls = []
    monkeypatch.setattr(
        "nanopore_converter.convert_run",
        lambda *a, **kw: calls.append(("convert_run", a, kw)) or True
    )
    src = _single_fast5(tmp_path / "data" / "test.fast5")
    archive = tmp_path / "data.tar.gz"
    with _tarfile.open(archive, "w:gz") as tar:
        tar.add(src, arcname="data/test.fast5")
    assert process_item(str(archive), str(tmp_path / "out"), threads=4) is True
    assert len(calls) == 1


def test_full_pipeline_dry_run(tmp_path, monkeypatch):
    """End-to-end with fake tools: a directory and an archive from a list."""
    FakeTools(monkeypatch)
    run_dir = tmp_path / "run1"
    _single_fast5(run_dir / "reads.fast5")
    src2 = tmp_path / "run2"
    _multi_fast5(src2 / "batch.fast5")
    archive = tmp_path / "run2.tar.gz"
    with _tarfile.open(archive, "w:gz") as tar:
        tar.add(src2, arcname="run2")
    listf = tmp_path / "inputs.txt"
    listf.write_text(f"{run_dir}\n{archive}\n")
    from nanopore_converter import main, count_pod5_reads
    assert main([str(listf), "--output-dir", str(tmp_path / "output")]) == 0
    assert count_pod5_reads(str(tmp_path / "output" / "run1" / "pod5")) == 1
    assert count_pod5_reads(str(tmp_path / "output" / "run2" / "pod5")) == 3


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
        f.create_dataset("read_00000000-0000-0000-0000-000000000000/Raw/Signal", data=[1])
    assert (d / "batch_0.fast5").stat().st_size < 1_000_000
    assert classify_fast5_dir(str(d)) == "multi_read_fast5"


def test_lost_reads_fail_the_run(tmp_path, monkeypatch):
    """Both tools exit 0 after dropping files; the read count must catch it."""
    from nanopore_converter import ConversionError, convert_run
    FakeTools(monkeypatch, drop=1)
    _single_fast5(tmp_path / "run1" / "a.fast5")
    _single_fast5(tmp_path / "run1" / "b.fast5")
    with pytest.raises(ConversionError, match="has 1 reads .* contain 2"):
        convert_run(str(tmp_path / "run1"), str(tmp_path / "out"))
    assert not (tmp_path / "out" / "run1" / "pod5").exists()


def test_mixed_single_and_multi_converts_all_reads(tmp_path, monkeypatch):
    """Regression: only the sampled kind used to be converted."""
    from nanopore_converter import convert_run, count_pod5_reads
    tools = FakeTools(monkeypatch)
    run = tmp_path / "runA"
    for i in range(4):
        _single_fast5(run / "fast5_pass" / f"r{i}.fast5")
    _multi_fast5(run / "fast5_fail" / "batch.fast5", 5)
    assert convert_run(str(run), str(tmp_path / "out")) is True
    assert [c[0] for c in tools.calls] == ["single_to_multi_fast5", "pod5"]
    assert count_pod5_reads(str(tmp_path / "out" / "runA" / "pod5")) == 9


def test_corrupt_file_is_skipped_and_reported(tmp_path, monkeypatch):
    """Regression: a truncated first file used to turn a run into 0 reads."""
    from nanopore_converter import convert_run, count_pod5_reads
    nanopore_converter.SKIPPED_UNREADABLE.clear()
    FakeTools(monkeypatch)
    run = tmp_path / "runM"
    (run / "fast5_pass").mkdir(parents=True)
    (run / "fast5_pass" / "aaa_truncated.fast5").write_bytes(b"\x89HDF" + b"x" * 50)
    _multi_fast5(run / "fast5_pass" / "batch_1.fast5", 6)
    assert convert_run(str(run), str(tmp_path / "out")) is True
    assert count_pod5_reads(str(tmp_path / "out" / "runM" / "pod5")) == 6
    assert [p.name for p in nanopore_converter.SKIPPED_UNREADABLE] == ["aaa_truncated.fast5"]
    nanopore_converter.SKIPPED_UNREADABLE.clear()


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

def _make_run(tmp_path, name="run1"):
    d = tmp_path / name
    _multi_fast5(d / "reads.fast5")
    return d


def _pod5_reads(path):
    import pod5
    with pod5.Reader(path) as r:
        return r.num_reads


def test_rerun_skips_converted_input(tmp_path, monkeypatch):
    from nanopore_converter import main
    tools = FakeTools(monkeypatch)
    run = _make_run(tmp_path)
    out = tmp_path / "out"
    assert main([str(run), "--output-dir", str(out)]) == 0
    assert len(tools.calls) == 1
    assert _pod5_reads(out / "run1" / "pod5" / "output.pod5") == 3
    # Second run: nothing is converted, exit code still 0
    assert main([str(run), "--output-dir", str(out)]) == 0
    assert len(tools.calls) == 1


def test_rerun_skips_archive_without_extracting(tmp_path, monkeypatch):
    from nanopore_converter import main
    tools = FakeTools(monkeypatch)
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
    assert extracted == [] and tools.calls == []


def test_force_replaces_existing_output(tmp_path, monkeypatch):
    from nanopore_converter import main
    FakeTools(monkeypatch)
    run = _make_run(tmp_path)
    pod5_dir = tmp_path / "out" / "run1" / "pod5"
    pod5_dir.mkdir(parents=True)
    (pod5_dir / "output.pod5").write_bytes(b"old")
    assert main([str(run), "--output-dir", str(tmp_path / "out"), "--force"]) == 0
    assert _pod5_reads(pod5_dir / "output.pod5") == 3
    assert sorted(p.name for p in pod5_dir.parent.iterdir()) == ["pod5"]


def test_force_keeps_old_output_when_conversion_fails(tmp_path, monkeypatch):
    from nanopore_converter import main
    FakeTools(monkeypatch, fail_pod5=True)
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
    FakeTools(monkeypatch, fail_pod5=True)
    run = _make_run(tmp_path)
    out = tmp_path / "out"
    assert main([str(run), "--output-dir", str(out)]) == 1
    assert not (out / "run1" / "pod5").exists()
    assert not is_converted(str(out), "run1")


def test_leftovers_from_interrupted_run_are_cleared(tmp_path, monkeypatch):
    from nanopore_converter import convert_run
    FakeTools(monkeypatch)
    run = tmp_path / "run1"
    _single_fast5(run / "read.fast5")
    out = tmp_path / "out"
    (out / "run1" / ".pod5_partial").mkdir(parents=True)
    (out / "run1" / ".pod5_partial" / "output.pod5").write_bytes(b"partial")
    (out / "run1" / "_multi_fast5_tmp").mkdir()
    _multi_fast5(out / "run1" / "_multi_fast5_tmp" / "stale.fast5", 9)
    (out / "run1" / "_fast5_stage").mkdir()
    assert convert_run(str(run), str(out)) is True
    assert _pod5_reads(out / "run1" / "pod5" / "output.pod5") == 1
    assert sorted(p.name for p in (out / "run1").iterdir()) == ["pod5"]


def test_duplicate_output_names_in_one_invocation_fail(tmp_path, monkeypatch):
    from nanopore_converter import main
    tools = FakeTools(monkeypatch)
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    a = _make_run(tmp_path / "a")  # a/run1
    b = _make_run(tmp_path / "b")  # b/run1 -> same output name
    listf = tmp_path / "inputs.txt"
    listf.write_text(f"{a}\n{b}\n")
    assert main([str(listf), "--output-dir", str(tmp_path / "out")]) == 1
    assert len(tools.calls) == 1  # only the first one was converted


def test_same_name_from_different_source_is_refused(tmp_path, monkeypatch):
    """Regression: run2/fast5_pass was skipped as run1/fast5_pass's output."""
    from nanopore_converter import main
    tools = FakeTools(monkeypatch)
    one = _make_run(tmp_path / "run1", "fast5_pass")
    two = _make_run(tmp_path / "run2", "fast5_pass")
    out = str(tmp_path / "out")
    assert main([str(one), "--output-dir", out]) == 0
    for extra in ([], ["--force"]):
        assert main([str(two), "--output-dir", out, *extra]) == 1
    assert len(tools.calls) == 1
    assert (tmp_path / "out" / "fast5_pass" / "pod5" / ".source").read_text().strip() \
        == str(one.resolve())


def test_output_overlapping_input_is_refused(tmp_path, monkeypatch):
    """Regression: --output-dir <parent> --force deleted the run's own pod5/."""
    from nanopore_converter import main
    tools = FakeTools(monkeypatch)
    run = _make_run(tmp_path / "runs", "runP")
    (run / "pod5").mkdir()
    (run / "pod5" / "original.pod5").write_bytes(b"original")
    assert main([str(run), "--output-dir", str(tmp_path / "runs"), "--force"]) == 1
    assert (run / "pod5" / "original.pod5").read_bytes() == b"original"
    assert tools.calls == []


def test_output_inside_input_is_not_rescanned(tmp_path, monkeypatch):
    from nanopore_converter import main, count_pod5_reads
    FakeTools(monkeypatch)
    run = _make_run(tmp_path, "run1")
    out = run / "converted"
    assert main([str(run), "--output-dir", str(out)]) == 0
    assert main([str(run), "--output-dir", str(out), "--force"]) == 0
    assert count_pod5_reads(str(out / "run1" / "pod5")) == 3


def test_list_file_paths_relative_to_list_and_deduplicated(tmp_path, monkeypatch):
    from nanopore_converter import main, parse_list_file
    tools = FakeTools(monkeypatch)
    data = tmp_path / "data"
    _make_run(data, "run1")
    _make_run(data, "run2")
    listf = data / "list.txt"
    listf.write_bytes(b"run1\r\nrun2/\r\n./run2\r\n")
    assert parse_list_file(str(listf)) == [str(data / "run1"), str(data / "run2/"),
                                           str(data / "./run2")]
    monkeypatch.chdir(tmp_path)  # run from elsewhere
    assert main([str(listf), "--output-dir", str(tmp_path / "out")]) == 0
    assert len(tools.calls) == 2


def test_run_surfaces_tool_problems_and_failures(caplog):
    from nanopore_converter import _run
    with caplog.at_level("WARNING"):
        _run(["sh", "-c", "echo 'progress 50%' >&2; echo 'ERROR: bad file x' >&2"], "tool")
    assert "ERROR: bad file x" in caplog.text and "progress" not in caplog.text
    with pytest.raises(subprocess.CalledProcessError):
        _run(["sh", "-c", "exit 3"], "tool")


def test_run_kills_tool_process_group_when_interrupted(tmp_path):
    """Regression: SIGTERM left the conversion tools (and workers) running."""
    import os
    import time
    from nanopore_converter import _run
    pidfile = tmp_path / "pids"

    class Interrupted(Exception):
        pass

    def on_alarm(signum, frame):
        raise Interrupted
    old = signal.signal(signal.SIGALRM, on_alarm)
    signal.setitimer(signal.ITIMER_REAL, 1.0)
    try:
        with pytest.raises(Interrupted):
            # parent tool plus a background "worker" in the same group
            _run(["sh", "-c", f"sleep 60 & echo $! > {pidfile}; sleep 60"], "tool")
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)
    worker = int(pidfile.read_text())
    time.sleep(0.2)
    # Dead means gone, or a zombie (killed but not yet reaped: in containers
    # PID 1 may never reap orphans).
    try:
        state = Path(f"/proc/{worker}/stat").read_text().split()[2]
    except FileNotFoundError:
        state = "gone"
    assert state in ("gone", "Z"), f"worker still running (state {state})"


# ---------------------------------------------------------------------------
# Filling missing pod5 metadata (dorado automatic model selection needs it)
# ---------------------------------------------------------------------------

def _write_pod5(path, flowcell="", kit="", sample_rate=4000, n_reads=3, read_ids=None,
                acquisition_id="run_a"):
    import uuid
    from datetime import datetime, timezone
    import numpy as np
    import pod5
    now = datetime.now(timezone.utc)
    run_info = pod5.RunInfo(
        acquisition_id=acquisition_id, acquisition_start_time=now, adc_max=4095, adc_min=-4096,
        context_tags={}, experiment_name="", flow_cell_id="FAK1",
        flow_cell_product_code=flowcell, protocol_name="", protocol_run_id="",
        protocol_start_time=now, sample_id="", sample_rate=sample_rate,
        sequencing_kit=kit, sequencer_position="", sequencer_position_type="",
        software="", system_name="", system_type="", tracking_id={},
    )
    with pod5.Writer(path) as writer:
        for rid in (read_ids if read_ids is not None else [uuid.uuid4() for _ in range(n_reads)]):
            writer.add_read(pod5.Read(
                read_id=uuid.UUID(str(rid)), pore=pod5.Pore(1, 1, "pore"),
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


# ---------------------------------------------------------------------------
# Audit round 3
# ---------------------------------------------------------------------------

def test_symlinked_folders_are_followed_once(tmp_path, monkeypatch):
    """Regression: a symlinked fast5_fail/ was silently left out."""
    import os
    from nanopore_converter import convert_run, count_pod5_reads
    FakeTools(monkeypatch)
    run = tmp_path / "run"
    for i in range(3):
        _single_fast5(run / "fast5_pass" / f"r{i}.fast5")
    _multi_fast5(tmp_path / "otherdisk" / "fast5_fail" / "b.fast5", 4)
    os.symlink(tmp_path / "otherdisk" / "fast5_fail", run / "fast5_fail")
    os.symlink(run, run / "loop")  # must not recurse forever or double count
    assert convert_run(str(run), str(tmp_path / "out")) is True
    assert count_pod5_reads(str(tmp_path / "out" / "run" / "pod5")) == 7


def test_unreadable_folder_fails_the_run(tmp_path, monkeypatch):
    """Regression: os.walk silently skipped folders it could not list."""
    import os
    from nanopore_converter import ConversionError, convert_run
    FakeTools(monkeypatch)
    run = tmp_path / "run"
    _single_fast5(run / "fast5_pass" / "r.fast5")
    _single_fast5(run / "fast5_fail" / "r.fast5")
    locked = str(run / "fast5_fail")
    real_scandir = os.scandir

    def scandir(path="."):
        if str(path) == locked:
            raise PermissionError(13, "Permission denied", locked)
        return real_scandir(path)

    monkeypatch.setattr(os, "scandir", scandir)
    with pytest.raises(ConversionError, match="cannot read 1 folder"):
        convert_run(str(run), str(tmp_path / "out"))
    assert not (tmp_path / "out" / "run" / "pod5").exists()


def test_duplicate_copies_are_converted_once(tmp_path, monkeypatch):
    """Regression: a second copy of the reads passed verification twice."""
    from nanopore_converter import convert_run
    FakeTools(monkeypatch)
    run = tmp_path / "run"
    ids = []
    for i in range(4):
        f = _single_fast5(run / "fast5_pass" / f"r{i}.fast5")
        ids.append(_read_ids_in([f.parent])[-1])
    _multi_fast5(run / "backup" / "batch.fast5", read_ids=ids)           # full copy
    for i in range(2):                                                    # single copies
        _single_fast5(run / "zz_copy" / f"r{i}.fast5", read_id=ids[i])
    _multi_fast5(run / ".convert_tmp" / "multi" / "b.fast5", read_ids=ids)  # hidden: ignored
    assert convert_run(str(run), str(tmp_path / "out")) is True
    out_ids = _pod5_ids(tmp_path / "out" / "run" / "pod5" / "output.pod5")
    assert sorted(out_ids) == sorted(ids)


def test_partial_duplicate_multi_file_fails(tmp_path, monkeypatch):
    from nanopore_converter import ConversionError, convert_run
    FakeTools(monkeypatch)
    run = tmp_path / "run"
    f = _single_fast5(run / "a" / "r0.fast5")
    shared = _read_ids_in([f.parent])[0]
    _multi_fast5(run / "b" / "batch.fast5", read_ids=[shared, "00000000-0000-0000-0000-000000000001"])
    with pytest.raises(ConversionError, match="share only some reads"):
        convert_run(str(run), str(tmp_path / "out"))


def test_duplicated_reads_in_output_fail_verification(tmp_path, monkeypatch):
    from nanopore_converter import ConversionError, convert_run
    FakeTools(monkeypatch, drop=1, duplicate=1)  # same count, one read doubled
    _multi_fast5(tmp_path / "run" / "b.fast5", 3)
    with pytest.raises(ConversionError, match="duplicated"):
        convert_run(str(tmp_path / "run"), str(tmp_path / "out"))


def test_macos_metadata_files_are_ignored(tmp_path, monkeypatch):
    from nanopore_converter import scan_fast5_dir
    run = tmp_path / "run"
    _multi_fast5(run / "fast5_pass" / "batch_0.fast5")
    (run / "fast5_pass" / "._batch_0.fast5").write_bytes(b"\x00\x05\x16\x07AppleDouble")
    (run / "__MACOSX" / "fast5_pass").mkdir(parents=True)
    (run / "__MACOSX" / "fast5_pass" / "._x.fast5").write_bytes(b"junk")
    scan = scan_fast5_dir(str(run))
    assert scan.unreadable == [] and len(scan.multi) == 1


@pytest.mark.parametrize("suffix,mode", [(".tar.xz", "w:xz"), (".txz", "w:xz"),
                                         (".tbz2", "w:bz2"), (".tar", "w")])
def test_every_accepted_archive_suffix_extracts(tmp_path, suffix, mode):
    """Regression: .tar.xz/.txz/.tbz2 were accepted but 'Unsupported'."""
    from nanopore_converter import detect_input_type, extract_archive
    src = _single_fast5(tmp_path / "data" / "r.fast5")
    archive = tmp_path / f"run{suffix}"
    with _tarfile.open(archive, mode) as tar:
        tar.add(src, arcname="run/r.fast5")
    assert detect_input_type(str(archive)) == "archive"
    out = extract_archive(str(archive))
    assert list(Path(out).rglob("*.fast5"))
    shutil.rmtree(out)


def test_metadata_fill_is_per_acquisition(tmp_path):
    """Regression: one run's flowcell was written into every run's reads."""
    import pod5
    from nanopore_converter import fill_missing_pod5_metadata
    d = tmp_path / "pod5"
    d.mkdir()
    _write_pod5(d / "a.pod5", acquisition_id="run_a", n_reads=2)
    _write_pod5(d / "b.pod5", acquisition_id="run_b", n_reads=2)
    per_acq = {"run_a": {"flowcell": "FLO-MIN106"}}  # run_b chemistry unknown
    assert fill_missing_pod5_metadata(str(d), per_acquisition=per_acq) == 1
    with pod5.Reader(d / "b.pod5") as r:
        assert {x.run_info.flow_cell_product_code for x in r.reads()} == {""}
    with pod5.Reader(d / "a.pod5") as r:
        assert {x.run_info.flow_cell_product_code for x in r.reads()} == {"FLO-MIN106"}
    # A single chemistry is not applied across several acquisitions
    assert fill_missing_pod5_metadata(str(d), chemistry={"flowcell": "FLO-MIN111"}) == 0


def test_sighup_stops_tools_and_cleans_up(tmp_path):
    """Regression: closing the terminal left the tools running as orphans."""
    import os
    import time
    stub = tmp_path / "bin"
    stub.mkdir()
    for tool in ("single_to_multi_fast5", "pod5"):
        f = stub / tool
        f.write_text(f"#!/bin/sh\necho $$ > {tmp_path}/{tool}.pid\nexec sleep 60\n")
        f.chmod(0o755)
    _single_fast5(tmp_path / "run" / "r.fast5")
    env = dict(os.environ, PATH=f"{stub}:{os.environ['PATH']}")
    proc = subprocess.Popen([sys.executable, str(Path(__file__).parent / "nanopore_converter.py"),
                             str(tmp_path / "run"), "--output-dir", str(tmp_path / "out")],
                            env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    pidfile = tmp_path / "single_to_multi_fast5.pid"
    for _ in range(100):
        if pidfile.exists() and pidfile.read_text().strip():
            break
        time.sleep(0.1)
    proc.send_signal(signal.SIGHUP)
    proc.wait(timeout=20)
    time.sleep(0.2)
    tool_pid = int(pidfile.read_text())
    try:
        state = Path(f"/proc/{tool_pid}/stat").read_text().split()[2]
    except FileNotFoundError:
        state = "gone"
    assert state in ("gone", "Z"), f"tool still running (state {state})"
    assert sorted(p.name for p in (tmp_path / "out" / "run").iterdir()) == []


def _pod5_ids(path):
    import pod5
    with pod5.Reader(path) as r:
        return [str(i) for i in r.read_ids]


# ---------------------------------------------------------------------------
# Audit round 4
# ---------------------------------------------------------------------------

def test_unconvertible_files_are_skipped_not_fatal(tmp_path, monkeypatch):
    """Regression: one file without read_id/Signal made the tool abort a batch
    and failed the whole run ('1 reads but ... 5')."""
    import h5py
    from nanopore_converter import convert_run, scan_fast5_dir
    run = tmp_path / "run"
    for i in range(3):
        _single_fast5(run / f"good{i}.fast5")
    no_id = _single_fast5(run / "no_id.fast5")
    with h5py.File(no_id, "a") as f:
        del f["Raw/Reads/Read_1"].attrs["read_id"]
    no_signal = _single_fast5(run / "no_signal.fast5")
    with h5py.File(no_signal, "a") as f:
        del f["Raw/Reads/Read_1/Signal"]
    no_raw = _multi_fast5(run / "multi_no_signal.fast5", 2)
    with h5py.File(no_raw, "a") as f:
        for k in list(f):
            del f[k]["Raw"]
    empty = run / "empty_multi.fast5"
    with h5py.File(empty, "w") as f:
        f.attrs["file_type"] = "multi-read"
    scan = scan_fast5_dir(str(run))
    assert sorted(p.name for p in scan.unreadable) == [
        "multi_no_signal.fast5", "no_id.fast5", "no_signal.fast5"]  # empty_multi: just no reads
    assert scan.expected_reads == 3 and scan.read_ids is not None
    FakeTools(monkeypatch)
    assert convert_run(str(run), str(tmp_path / "out")) is True


def test_stage_names_are_neutral(tmp_path):
    """Regression: .FAST5 was ignored by the tools; long names hit NAME_MAX."""
    import os
    from nanopore_converter import _stage
    src = [_single_fast5(tmp_path / "READ.FAST5"), _single_fast5(tmp_path / ("x" * 240 + ".fast5"))]
    _stage(src, str(tmp_path / "stage"))
    assert sorted(os.listdir(tmp_path / "stage")) == ["0000000.fast5", "0000001.fast5"]


def test_multi_read_file_from_two_runs_records_both(tmp_path):
    import uuid
    import h5py
    from nanopore_converter import scan_fast5_dir
    f = tmp_path / "run" / "merged.fast5"
    f.parent.mkdir()
    with h5py.File(f, "w") as h:
        for rid, run_id in ((uuid.uuid4(), "run_a"), (uuid.uuid4(), "run_b")):
            h.create_group(f"read_{rid}/tracking_id").attrs["run_id"] = run_id
            h.create_dataset(f"read_{rid}/Raw/Signal", data=[1])
    assert sorted(scan_fast5_dir(str(tmp_path / "run")).run_samples) == ["run_a", "run_b"]


def test_nohup_is_respected(tmp_path):
    """Regression: the SIGHUP handler overrode nohup, so logging out killed it."""
    import os
    import time
    stub = tmp_path / "bin"
    stub.mkdir()
    for tool in ("single_to_multi_fast5", "pod5"):
        (stub / tool).write_text("#!/bin/sh\nexec sleep 30\n")
        (stub / tool).chmod(0o755)
    _single_fast5(tmp_path / "run" / "r.fast5")
    env = dict(os.environ, PATH=f"{stub}:{os.environ['PATH']}")
    proc = subprocess.Popen(
        ["nohup", sys.executable, str(Path(__file__).parent / "nanopore_converter.py"),
         str(tmp_path / "run"), "--output-dir", str(tmp_path / "out")],
        env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=tmp_path)
    try:
        time.sleep(2)
        proc.send_signal(signal.SIGHUP)
        time.sleep(1)
        assert proc.poll() is None, "converter died on SIGHUP under nohup"
    finally:
        proc.terminate()
        proc.wait(timeout=20)
