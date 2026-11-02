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
