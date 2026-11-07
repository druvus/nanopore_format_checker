#!/usr/bin/env python3
"""
Compress subfolders (and copy existing archives) from an input directory to an output directory.

For each item directly inside the input folder:
  - Subdirectory  → compress it into a .tar.gz archive in the output folder
  - .tar file     → recompress to .tar.gz in the output folder
  - Compressed file (.tar.gz, .zip, etc.) → copy as-is to the output folder

Features:
  - Parallel processing with configurable worker threads
  - Thread-safe progress reporting
  - Atomic writes (temp file + rename) to avoid partial outputs
  - Skips items already present in the output directory (safe re-runs)
  - Graceful handling of permission errors (skips unreadable files and
    directories, keeping the rest of the folder)
  - Fast --dry-run that does not scan file trees
"""

import argparse
import logging
import os
import shutil
import signal
import sys
import tarfile
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

if sys.version_info < (3, 12):
    sys.exit(f"{Path(__file__).name} requires Python 3.12 or newer "
             f"(running {sys.version.split()[0]})")

COMPRESSED_EXTENSIONS = {
    ".tar.gz", ".tgz", ".tar.bz2", ".tbz2", ".tar.xz", ".txz",
    ".zip", ".gz", ".bz2", ".xz", ".7z", ".rar",
}

log = logging.getLogger("compress_folders")

# Set on Ctrl-C so archives in progress stop at the next entry
_STOP = threading.Event()


class Cancelled(Exception):
    """Raised inside a worker when the run is being interrupted."""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_archive_suffix(path: Path) -> str | None:
    """Return the compound or simple suffix if the file is a recognised compressed archive."""
    name_lower = path.name.lower()
    for ext in COMPRESSED_EXTENSIONS:
        if name_lower.endswith(ext):
            return ext
    return None


def _is_uncompressed_tar(path: Path) -> bool:
    """Return True if the file is an uncompressed .tar (not .tar.gz etc.)."""
    name_lower = path.name.lower()
    if not name_lower.endswith(".tar"):
        return False
    for ext in COMPRESSED_EXTENSIONS:
        if name_lower.endswith(ext):
            return False
    return True


def _format_size(nbytes: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(nbytes) < 1024:
            return f"{nbytes:.1f} {unit}"
        nbytes /= 1024
    return f"{nbytes:.1f} PiB"


def _dir_size(path: Path) -> tuple[int, int]:
    """Return (total_bytes, file_count) for a directory tree.
    Silently skips files/dirs that cannot be stat'd (permission errors)."""
    total = 0
    count = 0
    for root, _dirs, files in os.walk(path):
        if _STOP.is_set():
            raise Cancelled("interrupted")
        for fname in files:
            try:
                total += os.path.getsize(os.path.join(root, fname))
                count += 1
            except OSError:
                pass  # permission denied, broken symlink, etc.
    return total, count


def _dest_name_for_dir(src: Path) -> str:
    """Compute the output archive name for a directory.
    Handles dirs whose name already ends in .tar to avoid .tar.tar.gz."""
    name = src.name
    if name.lower().endswith(".tar"):
        return name + ".gz"  # foo.tar → foo.tar.gz
    return name + ".tar.gz"


def _dest_name_for_tar(src: Path) -> str:
    """Compute the output archive name for an uncompressed .tar file.
    foo.tar → foo.tar.gz"""
    return src.name + ".gz"


def _atomic_tmp(dest: Path) -> Path:
    """Return a unique temp path next to dest to avoid collisions between workers."""
    return dest.parent / f".{dest.name}.{uuid.uuid4().hex[:8]}.tmp"


# ---------------------------------------------------------------------------
# Thread-safe progress tracker
# ---------------------------------------------------------------------------

@dataclass
class ProgressTracker:
    total: int
    _done: int = 0
    _errors: int = 0
    _skipped: int = 0
    _bytes_in: int = 0
    _bytes_out: int = 0
    _permission_warnings: int = 0
    _lock: threading.Lock = field(default_factory=threading.Lock)
    _start: float = field(default_factory=time.monotonic)

    def complete(self, name: str, bytes_in: int = 0, bytes_out: int = 0) -> None:
        with self._lock:
            self._done += 1
            self._bytes_in += bytes_in
            self._bytes_out += bytes_out
            elapsed = time.monotonic() - self._start
            pct = (self._done / self.total) * 100
            log.info(
                "  ✓ %-55s  [%d/%d  %3.0f%%]  (%.1fs elapsed)",
                name, self._done, self.total, pct, elapsed,
            )

    def error(self, name: str) -> None:
        with self._lock:
            self._done += 1
            self._errors += 1
            log.error("  ✗ %-55s  [%d/%d]", name, self._done, self.total)

    def skip(self, name: str, reason: str = "already exists") -> None:
        with self._lock:
            self._done += 1
            self._skipped += 1
            log.info(
                "  ⏭ %-55s  [%d/%d] (%s)", name, self._done, self.total, reason
            )

    def add_permission_warning(self) -> None:
        with self._lock:
            self._permission_warnings += 1

    def summary(self, dry_run: bool) -> None:
        elapsed = time.monotonic() - self._start
        processed = self._done - self._skipped - self._errors
        log.info("─" * 72)
        log.info(
            "Finished: %d processed, %d skipped, %d errors  (%.1fs total)",
            processed, self._skipped, self._errors, elapsed,
        )
        if self._bytes_in:
            log.info(
                "Data: %s in → %s out",
                _format_size(self._bytes_in), _format_size(self._bytes_out),
            )
        if self._permission_warnings:
            log.warning(
                "Permission warnings: %d unreadable files/directories were skipped "
                "inside archives. Check source permissions if this is unexpected.",
                self._permission_warnings,
            )
        if dry_run:
            log.info("[DRY-RUN — no files were written]")

    @property
    def errors(self) -> int:
        return self._errors


# ---------------------------------------------------------------------------
# Task functions (each runs in its own thread)
# ---------------------------------------------------------------------------

def _add_tree(tar: tarfile.TarFile, src: Path, arcname: str, on_skip) -> None:
    """Add a directory tree to an open tar, skipping anything unreadable.

    Equivalent to tar.add(src, arcname) (sorted entries, symlinks stored as
    links and never followed), except that a file that cannot be opened or a
    directory that cannot be listed is reported via on_skip(arcname, error)
    instead of aborting the whole archive.  An unreadable directory is still
    stored as an (empty) directory entry so the tree structure is kept.
    """
    # A top-level symlink to a folder is archived as the folder's content,
    # not as a bare link (which would hold no data at all).
    stack = [(os.path.realpath(src) if os.path.islink(src) else str(src), arcname)]
    while stack:
        if _STOP.is_set():
            raise Cancelled("interrupted")
        path, arc = stack.pop()
        try:
            tarinfo = tar.gettarinfo(path, arc)
        except OSError as e:  # vanished or cannot be stat'ed
            on_skip(arc, e)
            continue
        if tarinfo is None:  # sockets etc., which tar.add also ignores
            continue
        if tarinfo.isreg():
            try:
                fh = open(path, "rb")
            except OSError as e:
                # gettarinfo() already registered this inode, so a later hard
                # link to it would be stored as a link to a missing member.
                # Look the inode up directly: scanning tar.inodes for every
                # unreadable file made large trees quadratic.
                try:
                    st = os.lstat(path)  # gettarinfo keys tar.inodes by (ino, dev)
                    key = (st.st_ino, st.st_dev)
                    if tar.inodes.get(key) == arc:
                        del tar.inodes[key]
                except OSError:
                    pass
                on_skip(arc, e)
                continue
            with fh:
                tar.addfile(tarinfo, fh)
            continue
        tar.addfile(tarinfo)
        if tarinfo.isdir():
            try:
                children = sorted(os.listdir(path))
            except OSError as e:
                on_skip(arc + "/", e)
                continue
            # Reverse so the stack pops entries in sorted order
            for child in reversed(children):
                stack.append((os.path.join(path, child), f"{arc}/{child}"))


def _copy_with_stop(src: Path, dest: Path, chunk: int = 16 * 1024 * 1024) -> None:
    """shutil.copy2, but in chunks so Ctrl-C (_STOP) can interrupt a huge copy."""
    with open(src, "rb") as fin, open(dest, "wb") as fout:
        while True:
            if _STOP.is_set():
                raise Cancelled("interrupted")
            block = fin.read(chunk)
            if not block:
                break
            fout.write(block)
    shutil.copystat(src, dest)


def _archive_top_folder(archive: Path) -> str | None:
    """Name of the first entry's top-level folder in a tar archive (archives
    made by compress_directory start with '<folder>/'); None if unreadable.
    Only the first header is read, so this is fast even for huge archives."""
    try:
        with tarfile.open(archive, "r:*") as tar:
            first = tar.next()
    except (tarfile.TarError, OSError, EOFError):
        return None
    return first.name.split("/", 1)[0] if first is not None else None


def compress_directory(
    src: Path, dest_archive: Path, dry_run: bool, progress: ProgressTracker
) -> None:
    """Compress a directory into a .tar.gz archive, skipping unreadable files and directories."""
    name = src.name
    try:
        # Skip for re-runs, but only if the existing archive was made from
        # this folder: an unrelated run1.tar.gz next to run1/ (common when the
        # output folder is the input folder) must not count as done.
        if dest_archive.exists():
            top = _archive_top_folder(dest_archive)
            if top == name:
                progress.skip(name, "output exists")
            else:
                found = f"its top-level entry is {top!r}" if top else "it cannot be read"
                log.error("  %s already exists but was not made from folder %s (%s); "
                          "rename one of them", dest_archive.name, name, found)
                progress.error(name)
            return

        if dry_run:
            log.info("  [DRY-RUN] Would compress dir: %s → %s", name, dest_archive.name)
            progress.complete(name)
            return

        # An unlistable (or listable but untraversable, mode r--) top-level
        # folder is an error, not an empty archive: an empty archive would be
        # skipped as "output exists" on every re-run.
        children = os.listdir(src)
        if children:
            os.lstat(os.path.join(src, children[0]))

        # Scan size (tolerant of permission errors)
        size, file_count = _dir_size(src)
        log.info("Compressing: %s (%s, %d files) → %s",
                 name, _format_size(size), file_count, dest_archive.name)

        skipped: list[str] = []

        def _on_skip(arcname: str, err: OSError) -> None:
            skipped.append(arcname)
            log.warning("    Unreadable, skipping: %s (%s)", arcname, err.strerror or err)

        tmp = _atomic_tmp(dest_archive)
        try:
            with tarfile.open(tmp, "w:gz") as tar:
                _add_tree(tar, src, name, _on_skip)
            tmp.rename(dest_archive)
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise

        if skipped:
            for _ in skipped:
                progress.add_permission_warning()
            log.warning("    %d unreadable files/directories skipped in %s", len(skipped), name)

        result_size = dest_archive.stat().st_size
        progress.complete(name, bytes_in=size, bytes_out=result_size)

    except Cancelled:
        log.info("  cancelled: %s", name)
    except Exception:
        log.exception("  ERROR compressing %s", name)
        progress.error(name)


def copy_archive(
    src: Path, dest: Path, dry_run: bool, progress: ProgressTracker
) -> None:
    """Copy an existing compressed file to the output directory."""
    name = src.name
    try:
        if dest.exists():
            progress.skip(name, "output exists")
            return

        if dry_run:
            size = src.stat().st_size
            log.info("  [DRY-RUN] Would copy: %s (%s) → %s",
                     name, _format_size(size), dest.name)
            progress.complete(name, bytes_in=size, bytes_out=size)
            return

        size = src.stat().st_size
        log.info("Copying: %s (%s)", name, _format_size(size))

        tmp = _atomic_tmp(dest)
        try:
            _copy_with_stop(src, tmp)
            tmp.rename(dest)
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise

        progress.complete(name, bytes_in=size, bytes_out=size)

    except Cancelled:
        log.info("  cancelled: %s", name)
    except Exception:
        log.exception("  ERROR copying %s", name)
        progress.error(name)


def recompress_tar(
    src: Path, dest: Path, dry_run: bool, progress: ProgressTracker
) -> None:
    """Recompress an uncompressed .tar file into a .tar.gz archive."""
    name = src.name
    try:
        if dest.exists():
            progress.skip(name, "output exists")
            return

        size = src.stat().st_size

        if dry_run:
            log.info("  [DRY-RUN] Would recompress: %s (%s) → %s",
                     name, _format_size(size), dest.name)
            progress.complete(name, bytes_in=size)
            return

        log.info("Recompressing: %s (%s) → %s", name, _format_size(size), dest.name)

        tmp = _atomic_tmp(dest)
        try:
            with tarfile.open(src, "r:") as src_tar:
                with tarfile.open(tmp, "w:gz") as dest_tar:
                    for member in src_tar:
                        if _STOP.is_set():
                            raise Cancelled("interrupted")
                        try:
                            fileobj = src_tar.extractfile(member) if member.isreg() else None
                            dest_tar.addfile(member, fileobj)
                        except (PermissionError, OSError) as e:
                            progress.add_permission_warning()
                            log.warning("    Skipping member %s: %s", member.name, e)
            tmp.rename(dest)
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise

        result_size = dest.stat().st_size
        progress.complete(name, bytes_in=size, bytes_out=result_size)

    except Cancelled:
        log.info("  cancelled: %s", name)
    except Exception:
        log.exception("  ERROR recompressing %s", name)
        progress.error(name)


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def _build_task_list(
    input_dir: Path, output_dir: Path
) -> tuple[list[tuple[str, Path, Path]], list[Path]]:
    """Classify items and compute (action, src, dest) tuples.
    Returns (tasks, skipped_files)."""
    items = sorted(input_dir.iterdir())
    tasks: list[tuple[str, Path, Path]] = []
    skipped: list[Path] = []
    out_real = output_dir.resolve()

    for p in items:
        try:
            real = p.resolve()
        except (OSError, RuntimeError) as e:  # symlink loop (RuntimeError on 3.12)
            log.warning("Skipping %s: %s", p.name, e)
            skipped.append(p)
            continue
        if real == out_real or real in out_real.parents:
            # The output folder (or a folder containing it) would be archived
            # into itself, including other workers' partial .tmp files.
            log.warning("Skipping %s: it contains the output directory", p.name)
            skipped.append(p)
            continue
        if p.is_dir():
            dest = output_dir / _dest_name_for_dir(p)
            tasks.append(("compress", p, dest))
        elif p.is_file() and _is_uncompressed_tar(p):
            dest = output_dir / _dest_name_for_tar(p)
            tasks.append(("recompress", p, dest))
        elif p.is_file() and _get_archive_suffix(p):
            dest = output_dir / p.name
            if real == dest.resolve():
                # Output dir == input dir: this archive is already in place
                # (possibly made by an earlier run of this tool).
                skipped.append(p)
                continue
            tasks.append(("copy", p, dest))
        else:
            skipped.append(p)

    return tasks, skipped


def _split_conflicts(
    tasks: list[tuple[str, Path, Path]]
) -> tuple[list[tuple[str, Path, Path]], list[tuple[str, Path, Path]]]:
    """Separate tasks whose destination is shared with another task (e.g.
    folder run1/ and file run1.tar.gz both -> run1.tar.gz); running both
    would let the last writer silently replace the other's archive."""
    by_dest: dict[Path, list] = {}
    for task in tasks:
        by_dest.setdefault(task[2], []).append(task)
    ok = [t for t in tasks if len(by_dest[t[2]]) == 1]
    conflicts = [t for t in tasks if len(by_dest[t[2]]) > 1]
    return ok, conflicts


def process(
    input_dir: Path, output_dir: Path, dry_run: bool, workers: int
) -> int:
    """Scan input_dir and process items in parallel. Returns error count."""
    if not input_dir.is_dir():
        log.error("Input path is not a directory: %s", input_dir)
        sys.exit(1)

    if not dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)

    log.info("Input:   %s", input_dir)
    log.info("Output:  %s", output_dir)
    log.info("Workers: %d", workers)

    tasks, skipped = _build_task_list(input_dir, output_dir)
    tasks, conflicts = _split_conflicts(tasks)

    # Pre-check what already exists in the output
    existing_outputs = set()
    if output_dir.is_dir():
        existing_outputs = {p.name for p in output_dir.iterdir()}

    # Also check for in-progress .tmp files that might indicate a crashed run
    tmp_files = [n for n in existing_outputs if n.endswith(".tmp")]
    if tmp_files:
        log.warning(
            "Found %d .tmp files in output (possibly from a crashed run). "
            "Consider removing them: %s",
            len(tmp_files), output_dir,
        )

    n_already_done = sum(1 for _, _, dest in tasks if dest.name in existing_outputs)

    # Count by type
    counts = {"compress": 0, "recompress": 0, "copy": 0}
    for action, _, _ in tasks:
        counts[action] += 1

    log.info(
        "Found %d directories, %d compressed archives, %d uncompressed .tar files, "
        "%d other files (ignored)",
        counts["compress"], counts["copy"], counts["recompress"], len(skipped),
    )
    if n_already_done:
        log.info(
            "Already in output: %d/%d items (will be skipped)",
            n_already_done, len(tasks),
        )

    total = len(tasks) + len(conflicts)
    if total == 0:
        log.warning("Nothing to process.")
        return 0

    progress = ProgressTracker(total=total)
    log.info("─" * 72)
    for _action, src, dest in conflicts:
        log.error("  %s and another item would both be written to %s; "
                  "rename one of them", src.name, dest.name)
        progress.error(src.name)

    dispatch = {
        "compress": compress_directory,
        "recompress": recompress_tar,
        "copy": copy_archive,
    }

    pool = ThreadPoolExecutor(max_workers=workers, thread_name_prefix="worker")
    try:
        futures = []
        for action, src, dest in tasks:
            fn = dispatch[action]
            futures.append(pool.submit(fn, src, dest, dry_run, progress))

        for fut in as_completed(futures):
            exc = fut.exception()
            if exc:
                log.error("Unexpected worker error: %s", exc)
    except KeyboardInterrupt:
        # Without this, the executor would still run every queued item.
        log.warning("Interrupted: cancelling queued items, stopping running ones...")
        _STOP.set()
        pool.shutdown(wait=True, cancel_futures=True)
        raise
    finally:
        pool.shutdown(wait=True)

    if skipped:
        log.debug(
            "Ignored non-archive files: %s", ", ".join(s.name for s in skipped)
        )

    progress.summary(dry_run)
    return progress.errors


def _install_signal_handlers() -> None:
    """SIGTERM (scheduler time limits, kill) and SIGHUP get the same handling
    as Ctrl-C: queued items are cancelled, archives in progress are stopped
    and their temp files removed. Signals already ignored (nohup) stay so."""
    def handler(signum, frame):
        raise KeyboardInterrupt
    for sig in (signal.SIGTERM, signal.SIGHUP):
        if signal.getsignal(sig) is not signal.SIG_IGN:
            signal.signal(sig, handler)


def main() -> None:
    _install_signal_handlers()
    cpu_count = os.cpu_count() or 4
    default_workers = min(cpu_count, 8)

    parser = argparse.ArgumentParser(
        description="Compress subfolders and collect archives from an input directory.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  %(prog)s /data/raw /data/compressed --dry-run     # preview
  %(prog)s /data/raw /data/compressed                # run (skips existing)
  %(prog)s /data/raw /data/compressed -w 4 -v        # 4 workers, verbose
""",
    )
    parser.add_argument("input_dir", type=Path, help="Source directory to scan")
    parser.add_argument("output_dir", type=Path, help="Destination for archives")
    parser.add_argument(
        "-n", "--dry-run", action="store_true",
        help="Preview actions without creating or copying any files",
    )
    parser.add_argument(
        "-w", "--workers", type=int, default=default_workers,
        help=f"Number of parallel workers (default: {default_workers})",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable debug-level logging",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)-5s] %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.dry_run:
        log.info("═" * 72)
        log.info("  DRY-RUN MODE — no files will be written")
        log.info("═" * 72)

    try:
        errors = process(
            args.input_dir.resolve(),
            args.output_dir.resolve(),
            args.dry_run,
            args.workers,
        )
    except KeyboardInterrupt:
        log.warning("Interrupted; archives in progress were discarded (no partial files left). "
                    "Re-run to continue: finished archives are skipped.")
        sys.exit(130)
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()