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
  - Graceful handling of permission errors (skips unreadable files)
  - Fast --dry-run that does not scan file trees
"""

import argparse
import logging
import os
import shutil
import sys
import tarfile
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

COMPRESSED_EXTENSIONS = {
    ".tar.gz", ".tgz", ".tar.bz2", ".tbz2", ".tar.xz", ".txz",
    ".zip", ".gz", ".bz2", ".xz", ".7z", ".rar",
}

log = logging.getLogger("compress_folders")


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
                "Permission warnings: %d files were skipped inside archives "
                "(unreadable files). Check source permissions if this is unexpected.",
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

def compress_directory(
    src: Path, dest_archive: Path, dry_run: bool, progress: ProgressTracker
) -> None:
    """Compress a directory into a .tar.gz archive, skipping unreadable files."""
    name = src.name
    try:
        # Check if output already exists (fast skip for re-runs)
        if dest_archive.exists():
            progress.skip(name, "output exists")
            return

        if dry_run:
            log.info("  [DRY-RUN] Would compress dir: %s → %s", name, dest_archive.name)
            progress.complete(name)
            return

        # Scan size (tolerant of permission errors)
        size, file_count = _dir_size(src)
        log.info("Compressing: %s (%s, %d files) → %s",
                 name, _format_size(size), file_count, dest_archive.name)

        # Custom error handler for tar.add — logs and skips unreadable files
        perm_errors = 0

        def _onerror(tarinfo: tarfile.TarInfo) -> tarfile.TarInfo | None:
            """Filter callback: skip files we can't read."""
            nonlocal perm_errors
            fpath = os.path.join(src, tarinfo.name) if not os.path.isabs(tarinfo.name) else tarinfo.name
            if tarinfo.isfile():
                try:
                    with open(os.path.join(str(src.parent), tarinfo.name), "rb") as f:
                        f.read(1)
                except PermissionError:
                    perm_errors += 1
                    log.warning("    Permission denied, skipping: %s", tarinfo.name)
                    return None
                except OSError as e:
                    perm_errors += 1
                    log.warning("    OS error, skipping: %s (%s)", tarinfo.name, e)
                    return None
            return tarinfo

        tmp = _atomic_tmp(dest_archive)
        try:
            with tarfile.open(tmp, "w:gz") as tar:
                tar.add(src, arcname=name, filter=_onerror)
            tmp.rename(dest_archive)
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise

        if perm_errors:
            for _ in range(perm_errors):
                progress.add_permission_warning()
            log.warning("    %d files skipped due to permissions in %s", perm_errors, name)

        result_size = dest_archive.stat().st_size
        progress.complete(name, bytes_in=size, bytes_out=result_size)

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
            shutil.copy2(src, tmp)
            tmp.rename(dest)
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise

        progress.complete(name, bytes_in=size, bytes_out=size)

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

    for p in items:
        if p.is_dir():
            dest = output_dir / _dest_name_for_dir(p)
            tasks.append(("compress", p, dest))
        elif p.is_file() and _is_uncompressed_tar(p):
            dest = output_dir / _dest_name_for_tar(p)
            tasks.append(("recompress", p, dest))
        elif p.is_file() and _get_archive_suffix(p):
            dest = output_dir / p.name
            tasks.append(("copy", p, dest))
        else:
            skipped.append(p)

    return tasks, skipped


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

    total = len(tasks)
    if total == 0:
        log.warning("Nothing to process.")
        return 0

    progress = ProgressTracker(total=total)
    log.info("─" * 72)

    dispatch = {
        "compress": compress_directory,
        "recompress": recompress_tar,
        "copy": copy_archive,
    }

    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="worker") as pool:
        futures = []
        for action, src, dest in tasks:
            fn = dispatch[action]
            futures.append(pool.submit(fn, src, dest, dry_run, progress))

        for fut in as_completed(futures):
            exc = fut.exception()
            if exc:
                log.error("Unexpected worker error: %s", exc)

    if skipped:
        log.debug(
            "Ignored non-archive files: %s", ", ".join(s.name for s in skipped)
        )

    progress.summary(dry_run)
    return progress.errors


def main() -> None:
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

    errors = process(
        args.input_dir.resolve(),
        args.output_dir.resolve(),
        args.dry_run,
        args.workers,
    )
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()