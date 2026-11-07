#!/usr/bin/env python3
"""Tests for compress_folders.py"""

import os
import tarfile
from pathlib import Path

import pytest

import compress_folders as cf


def _tree(root: Path) -> Path:
    """run/ with nested dirs, files and symlinks to a file and a directory."""
    src = root / "run"
    (src / "a" / "b").mkdir(parents=True)
    (src / "empty").mkdir()
    (src / "top.txt").write_text("top")
    (src / "a" / "one.txt").write_text("one")
    (src / "a" / "b" / "two.txt").write_text("two")
    os.symlink("top.txt", src / "link_to_file")
    os.symlink("a", src / "link_to_dir")
    return src


def _members(archive: Path) -> list[tuple[str, str]]:
    kinds = {tarfile.REGTYPE: "f", tarfile.DIRTYPE: "d", tarfile.SYMTYPE: "l"}
    with tarfile.open(archive) as tar:
        return [(m.name, kinds.get(m.type, "?")) for m in tar.getmembers()]


def _compress(src: Path, out: Path) -> cf.ProgressTracker:
    out.mkdir(exist_ok=True)
    progress = cf.ProgressTracker(total=1)
    cf.compress_directory(src, out / f"{src.name}.tar.gz", dry_run=False, progress=progress)
    return progress


def test_add_tree_matches_tar_add(tmp_path):
    """On a fully readable tree, the archive matches tarfile.add()."""
    src = _tree(tmp_path)
    ref = tmp_path / "ref.tar.gz"
    with tarfile.open(ref, "w:gz") as tar:
        tar.add(src, arcname="run")
    progress = _compress(src, tmp_path / "out")
    ours = tmp_path / "out" / "run.tar.gz"
    assert _members(ours) == _members(ref)
    assert progress.errors == 0 and progress._permission_warnings == 0
    with tarfile.open(ours) as tar:
        assert tar.extractfile("run/a/b/two.txt").read() == b"two"
        assert tar.getmember("run/link_to_dir").linkname == "a"  # not followed


def test_unreadable_subdirectory_is_skipped_not_fatal(tmp_path, monkeypatch):
    """Regression: one unlistable subfolder used to fail the whole archive."""
    src = _tree(tmp_path)
    locked = str(src / "a" / "b")
    real_listdir = os.listdir

    def listdir(path):
        if str(path) == locked:
            raise PermissionError(13, "Permission denied", path)
        return real_listdir(path)

    monkeypatch.setattr(cf.os, "listdir", listdir)
    progress = _compress(src, tmp_path / "out")
    names = dict(_members(tmp_path / "out" / "run.tar.gz"))
    assert progress.errors == 0
    assert progress._permission_warnings == 1
    assert names["run/a/b"] == "d"            # kept as an empty directory
    assert "run/a/b/two.txt" not in names     # its contents are skipped
    assert names["run/a/one.txt"] == "f"      # everything else is archived
    assert names["run/top.txt"] == "f"


def test_unreadable_file_is_skipped(tmp_path, monkeypatch):
    src = _tree(tmp_path)
    secret = str(src / "a" / "one.txt")

    def fake_open(path, *args, **kwargs):
        if str(path) == secret:
            raise PermissionError(13, "Permission denied", path)
        return open(path, *args, **kwargs)

    monkeypatch.setattr(cf, "open", fake_open, raising=False)
    progress = _compress(src, tmp_path / "out")
    names = dict(_members(tmp_path / "out" / "run.tar.gz"))
    assert progress.errors == 0 and progress._permission_warnings == 1
    assert "run/a/one.txt" not in names
    assert names["run/a/b/two.txt"] == "f"


def test_unreadable_top_level_directory_is_an_error(tmp_path, monkeypatch):
    """No (empty) archive for a folder that cannot be listed at all,
    otherwise every re-run would skip it as 'output exists'."""
    src = _tree(tmp_path)
    real_listdir = os.listdir

    def listdir(path):
        if str(path) == str(src):
            raise PermissionError(13, "Permission denied", path)
        return real_listdir(path)

    monkeypatch.setattr(cf.os, "listdir", listdir)
    out = tmp_path / "out"
    progress = _compress(src, out)
    assert progress.errors == 1
    assert list(out.iterdir()) == []  # no archive and no leftover .tmp


@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0,
                    reason="root bypasses file permissions")
def test_real_permissions(tmp_path):
    """End to end with real chmod 000 entries (runs on non-root CI)."""
    src = _tree(tmp_path)
    locked_dir = src / "a" / "b"
    secret = src / "top.txt"
    locked_dir.chmod(0)
    secret.chmod(0)
    try:
        progress = _compress(src, tmp_path / "out")
    finally:
        locked_dir.chmod(0o755)
        secret.chmod(0o644)
    names = dict(_members(tmp_path / "out" / "run.tar.gz"))
    assert progress.errors == 0 and progress._permission_warnings == 2
    assert names["run/a/b"] == "d" and "run/a/b/two.txt" not in names
    assert "run/top.txt" not in names
    assert names["run/a/one.txt"] == "f"


# ---------------------------------------------------------------------------
# Audit round 2
# ---------------------------------------------------------------------------

def test_top_level_symlink_archives_its_content(tmp_path):
    """Regression: a symlinked run folder was archived as a bare link."""
    real = _tree(tmp_path / "store")
    (tmp_path / "in").mkdir()
    os.symlink(real, tmp_path / "in" / "linked")
    _compress(tmp_path / "in" / "linked", tmp_path / "out")
    names = dict(_members(tmp_path / "out" / "linked.tar.gz"))
    assert names["linked"] == "d"
    assert names["linked/a/b/two.txt"] == "f"


def test_items_sharing_a_destination_are_errors(tmp_path):
    """Regression: run1/ and run1.tar.gz both wrote out/run1.tar.gz."""
    src = tmp_path / "in"
    (src / "run1").mkdir(parents=True)
    (src / "run1" / "a.txt").write_text("a")
    (src / "run1.tar.gz").write_bytes(b"existing archive")
    (src / "run2").mkdir()
    errors = cf.process(src, tmp_path / "out", dry_run=False, workers=2)
    assert errors == 2
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == ["run2.tar.gz"]


def test_output_inside_input_is_not_archived(tmp_path):
    src = tmp_path / "in"
    (src / "runA").mkdir(parents=True)
    (src / "runA" / "x.txt").write_text("x")
    out = src / "archives"
    assert cf.process(src, out, dry_run=False, workers=2) == 0
    assert sorted(p.name for p in out.iterdir()) == ["runA.tar.gz"]


def test_hardlink_to_skipped_file_is_stored_as_file(tmp_path, monkeypatch):
    """Regression: the second link pointed at a member missing from the tar."""
    src = tmp_path / "runH"
    src.mkdir()
    (src / "a.fast5").write_text("data")
    os.link(src / "a.fast5", src / "b.fast5")
    locked = str(src / "a.fast5")
    calls = {"n": 0}

    def fake_open(path, *args, **kwargs):
        if str(path) == locked:
            raise PermissionError(13, "Permission denied", path)
        return open(path, *args, **kwargs)

    monkeypatch.setattr(cf, "open", fake_open, raising=False)
    _compress(src, tmp_path / "out")
    archive = tmp_path / "out" / "runH.tar.gz"
    assert dict(_members(archive)) == {"runH": "d", "runH/b.fast5": "f"}
    with tarfile.open(archive) as tar:  # extracts without dangling hard links
        tar.extractall(tmp_path / "x", filter="data")
    assert (tmp_path / "x" / "runH" / "b.fast5").read_text() == "data"


def test_untraversable_top_level_directory_is_an_error(tmp_path, monkeypatch):
    """Regression: mode r-- (listable, not traversable) gave an empty archive."""
    src = _tree(tmp_path)
    real_lstat = os.lstat

    def lstat(path, *a, **kw):
        if os.path.dirname(str(path)) == str(src):
            raise PermissionError(13, "Permission denied", path)
        return real_lstat(path, *a, **kw)

    monkeypatch.setattr(cf.os, "lstat", lstat)
    out = tmp_path / "out"
    progress = _compress(src, out)
    assert progress.errors == 1 and list(out.iterdir()) == []


def test_interrupt_cancels_queued_items(tmp_path, monkeypatch):
    """Regression: Ctrl-C still ran every queued item to completion."""
    import time
    src = tmp_path / "in"
    for i in range(6):
        (src / f"run{i}").mkdir(parents=True)
    started = []

    def slow(srcdir, dest, dry_run, progress):
        started.append(srcdir.name)
        time.sleep(0.3)

    monkeypatch.setattr(cf, "compress_directory", slow)

    def interrupted(futures):
        raise KeyboardInterrupt
        yield  # pragma: no cover

    monkeypatch.setattr(cf, "as_completed", interrupted)
    cf._STOP.clear()
    try:
        with pytest.raises(KeyboardInterrupt):
            cf.process(src, tmp_path / "out", dry_run=False, workers=1)
        assert cf._STOP.is_set()
        assert len(started) < 6
    finally:
        cf._STOP.clear()


@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0,
                    reason="root bypasses file permissions")
def test_real_untraversable_top_level_directory(tmp_path):
    src = _tree(tmp_path)
    src.chmod(0o644)  # r-- for owner: listable, not traversable
    try:
        progress = _compress(src, tmp_path / "out")
    finally:
        src.chmod(0o755)
    assert progress.errors == 1
    assert not (tmp_path / "out" / "run.tar.gz").exists()


# ---------------------------------------------------------------------------
# Audit round 3
# ---------------------------------------------------------------------------

def test_output_dir_equal_to_input_dir_reruns_cleanly(tmp_path):
    """Regression: the second run saw its own archives as clashing inputs."""
    d = tmp_path / "d"
    (d / "run1").mkdir(parents=True)
    (d / "run1" / "a.txt").write_text("a")
    assert cf.process(d, d, dry_run=False, workers=2) == 0
    assert cf.process(d, d, dry_run=False, workers=2) == 0
    assert (d / "run1.tar.gz").exists()


@pytest.mark.parametrize("worker", ["recompress_tar", "copy_archive"])
def test_cancel_stops_recompress_and_copy_quietly(tmp_path, worker, caplog):
    """Regression: Ctrl-C did not stop .tar recompression or archive copies,
    and cancelled items were logged as errors with tracebacks."""
    src_dir = tmp_path / "data"
    src_dir.mkdir()
    (src_dir / "f.bin").write_bytes(b"x" * 1000)
    tar_path = tmp_path / "in.tar"
    with tarfile.open(tar_path, "w") as tar:
        tar.add(src_dir / "f.bin", arcname="f.bin")
    out = tmp_path / "out"
    out.mkdir()
    dest = out / ("in.tar.gz" if worker == "recompress_tar" else "in.tar")
    progress = cf.ProgressTracker(total=1)
    cf._STOP.set()
    try:
        with caplog.at_level("INFO"):
            getattr(cf, worker)(tar_path, dest, False, progress)
    finally:
        cf._STOP.clear()
    assert not dest.exists() and list(out.iterdir()) == []
    assert progress.errors == 0
    assert "cancelled: in.tar" in caplog.text and "Traceback" not in caplog.text


# ---------------------------------------------------------------------------
# Audit round 4
# ---------------------------------------------------------------------------

def test_many_unreadable_files_stay_linear(tmp_path, monkeypatch):
    """Regression: each unreadable file scanned all of tar.inodes (quadratic)."""
    import time
    src = tmp_path / "run"
    (src / "a_ok").mkdir(parents=True)
    (src / "z_locked").mkdir()
    # tarfile only tracks inodes of hard-linked files (st_nlink > 1), so give
    # every file a second link outside the archived tree.
    (tmp_path / "links").mkdir()
    # unreadable files come after the readable ones (sorted order)
    for sub, n in (("a_ok", 20000), ("z_locked", 10000)):
        for i in range(n):
            f = src / sub / f"f{i}"
            f.write_bytes(b"x")
            os.link(f, tmp_path / "links" / f"{sub}{i}")

    def fake_open(path, *a, **k):
        if "/z_locked/" in str(path):
            raise PermissionError(13, "denied", path)
        return open(path, *a, **k)

    monkeypatch.setattr(cf, "open", fake_open, raising=False)
    start = time.perf_counter()
    progress = _compress(src, tmp_path / "out")
    assert progress._permission_warnings == 10000
    assert time.perf_counter() - start < 15  # quadratic version: >20 s


def test_symlink_loop_at_top_level_is_skipped(tmp_path):
    src = tmp_path / "in"
    (src / "run1").mkdir(parents=True)
    os.symlink("loop", src / "loop")  # self-referencing
    assert cf.process(src, tmp_path / "out", dry_run=False, workers=1) == 0
    assert (tmp_path / "out" / "run1.tar.gz").exists()


def test_sigterm_cleans_up_like_ctrl_c(tmp_path):
    """Regression: SIGTERM left .tmp files in the output."""
    import signal as _signal
    import subprocess
    import sys
    import time
    src = tmp_path / "in" / "big"
    src.mkdir(parents=True)
    for i in range(40):
        (src / f"f{i}").write_bytes(os.urandom(2_000_000))
    out = tmp_path / "out"
    proc = subprocess.Popen([sys.executable, str(Path(cf.__file__)), str(tmp_path / "in"), str(out)],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(100):
        if out.exists() and any(p.name.endswith(".tmp") for p in out.iterdir()):
            break
        time.sleep(0.05)
    proc.send_signal(_signal.SIGTERM)
    assert proc.wait(timeout=30) == 130
    assert not [p for p in out.iterdir() if p.name.endswith(".tmp")]


def test_unrelated_archive_with_same_name_is_not_skipped(tmp_path):
    """Regression: with output = input, an unrelated run1.tar.gz counted as done."""
    d = tmp_path / "d"
    (d / "run1").mkdir(parents=True)
    (d / "run1" / "a.txt").write_text("a")
    other = tmp_path / "other"
    other.mkdir()
    (other / "x.txt").write_text("x")
    with tarfile.open(d / "run1.tar.gz", "w:gz") as tar:
        tar.add(other, arcname="other")
    before = (d / "run1.tar.gz").read_bytes()
    progress = _compress(d / "run1", d)
    assert progress._errors == 1 and progress._skipped == 0
    assert (d / "run1.tar.gz").read_bytes() == before  # never overwritten
    # An archive made from the folder itself is still skipped on re-runs
    (d / "run1.tar.gz").unlink()
    assert _compress(d / "run1", d)._errors == 0
    progress = _compress(d / "run1", d)
    assert progress._skipped == 1 and progress._errors == 0
