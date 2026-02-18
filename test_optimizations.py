#!/usr/bin/env python3
"""
Tests for the optimized nanopore_format_checker.

Creates mock directory trees to verify all format detection paths work
correctly after the os.scandir optimization changes.
"""

import os
import stat
import sys
import tempfile
from pathlib import Path

# Import the checker module from the same directory
sys.path.insert(0, str(Path(__file__).parent))
import io

from nanopore_format_checker import (
    _has_file_with_ext,
    analyze_run,
    compute_dir_size,
    diagnose_unknown,
    discover_run_structure,
    estimate_dir_size,
    fast_count_files,
    find_named_subdirs,
    format_size,
    generate_conversion_script,
    print_conversion_help,
    write_stats_tsv,
)


def make_file(path: Path, size: int = 100) -> None:
    """Create a file with the given size in bytes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x00" * size)


def test_pod5_subdirectory(tmp_path: Path) -> None:
    """Run with pod5/ subdirectory containing .pod5 files."""
    run = tmp_path / "20240101_run_pod5"
    (run / "pod5").mkdir(parents=True)
    make_file(run / "pod5" / "reads_001.pod5", size=5_000_000)
    make_file(run / "pod5" / "reads_002.pod5", size=5_000_000)

    result = analyze_run(run)
    assert "pod5" in result["formats"], f"Expected pod5, got {result['formats']}"
    assert result["details"]["pod5"]["file_count"] == 2
    print("  PASS: pod5 subdirectory")


def test_pod5_variant_dirs(tmp_path: Path) -> None:
    """Run with pod5_pass/ and pod5_fail/ variant directories."""
    run = tmp_path / "20240101_run_pod5var"
    (run / "pod5_pass").mkdir(parents=True)
    (run / "pod5_fail").mkdir(parents=True)
    make_file(run / "pod5_pass" / "reads.pod5", size=5_000_000)

    result = analyze_run(run)
    assert "pod5" in result["formats"], f"Expected pod5, got {result['formats']}"
    assert "pod5_fail" in result["details"]["pod5"]["folder_variants"]
    assert "pod5_pass" in result["details"]["pod5"]["folder_variants"]
    print("  PASS: pod5 variant directories")


def test_multi_read_fast5(tmp_path: Path) -> None:
    """Run with fast5/ containing multi-read fast5 files (>1MB each)."""
    run = tmp_path / "20240101_run_multi"
    (run / "fast5").mkdir(parents=True)
    # Multi-read fast5 files are >1MB
    make_file(run / "fast5" / "batch_001.fast5", size=2_000_000)
    make_file(run / "fast5" / "batch_002.fast5", size=2_000_000)

    result = analyze_run(run)
    assert "multi_read_fast5" in result["formats"], f"Expected multi_read_fast5, got {result['formats']}"
    assert result["details"]["multi_read_fast5"]["file_count"] == 2
    print("  PASS: multi-read fast5 in fast5/ subdirectory")


def test_single_read_fast5_in_fast5_dir(tmp_path: Path) -> None:
    """Run with fast5/ containing single-read fast5 files (<1MB each)."""
    run = tmp_path / "20240101_run_single_f5dir"
    (run / "fast5").mkdir(parents=True)
    # Single-read fast5 files are <1MB (typically 1-50KB)
    for i in range(10):
        make_file(run / "fast5" / f"read_{i:04d}.fast5", size=30_000)

    result = analyze_run(run)
    assert "single_read_fast5" in result["formats"], f"Expected single_read_fast5, got {result['formats']}"
    print("  PASS: single-read fast5 in fast5/ subdirectory")


def test_single_read_fast5_in_root(tmp_path: Path) -> None:
    """Run with .fast5 files directly in root (no fast5/ subdirectory)."""
    run = tmp_path / "20240101_run_single_root"
    run.mkdir(parents=True)
    for i in range(10):
        make_file(run / f"read_{i:04d}.fast5", size=30_000)

    result = analyze_run(run)
    assert "single_read_fast5" in result["formats"], f"Expected single_read_fast5, got {result['formats']}"
    assert result["details"]["single_read_fast5"]["layout"] == "root"
    print("  PASS: single-read fast5 in root")


def test_single_read_fast5_numeric_subdirs(tmp_path: Path) -> None:
    """Run with numeric subdirs (0/, 1/, 2/) containing .fast5 files."""
    run = tmp_path / "20240101_run_single_numeric"
    run.mkdir(parents=True)
    for subdir in range(3):
        d = run / str(subdir)
        d.mkdir()
        for i in range(5):
            make_file(d / f"read_{i:04d}.fast5", size=30_000)

    result = analyze_run(run)
    assert "single_read_fast5" in result["formats"], f"Expected single_read_fast5, got {result['formats']}"
    assert result["details"]["single_read_fast5"]["layout"] == "numeric_subdirs"
    print("  PASS: single-read fast5 in numeric subdirs")


def test_fast5_variant_dirs(tmp_path: Path) -> None:
    """Run with fast5_pass/ and fast5_fail/ variant directories."""
    run = tmp_path / "20240101_run_f5var"
    (run / "fast5_pass").mkdir(parents=True)
    (run / "fast5_fail").mkdir(parents=True)
    make_file(run / "fast5_pass" / "batch.fast5", size=2_000_000)

    result = analyze_run(run)
    assert "multi_read_fast5" in result["formats"], f"Expected multi_read_fast5, got {result['formats']}"
    assert "fast5_fail" in result["details"]["multi_read_fast5"]["folder_variants"]
    assert "fast5_pass" in result["details"]["multi_read_fast5"]["folder_variants"]
    print("  PASS: fast5 variant directories")


def test_compressed_archives(tmp_path: Path) -> None:
    """Run with compressed archives in run directory."""
    run = tmp_path / "20240101_run_archives"
    run.mkdir(parents=True)
    make_file(run / "reads.tar.gz", size=1000)
    make_file(run / "more_reads.tgz", size=1000)

    result = analyze_run(run)
    assert "single_read_fast5" in result["formats"], f"Expected single_read_fast5, got {result['formats']}"
    assert result["details"]["single_read_fast5"]["archive_count"] == 2
    print("  PASS: compressed archives")


def test_empty_directory(tmp_path: Path) -> None:
    """Run with empty directory."""
    run = tmp_path / "20240101_run_empty"
    run.mkdir(parents=True)

    result = analyze_run(run)
    assert "unknown" in result["formats"], f"Expected unknown, got {result['formats']}"
    print("  PASS: empty directory")


def test_fastq_directory(tmp_path: Path) -> None:
    """Run with fastq_pass/ directory."""
    run = tmp_path / "20240101_run_fastq"
    (run / "fastq_pass").mkdir(parents=True)
    make_file(run / "fastq_pass" / "reads_001.fastq.gz", size=1000)
    make_file(run / "fastq_pass" / "reads_002.fastq.gz", size=1000)
    make_file(run / "fastq_pass" / "reads_003.fq", size=1000)

    result = analyze_run(run)
    assert "fastq" in result["formats"], f"Expected fastq, got {result['formats']}"
    assert result["details"]["fastq"]["file_count"] == 3
    print("  PASS: fastq directory")


def test_mixed_formats(tmp_path: Path) -> None:
    """Run with both pod5 and fast5 directories (both formats detected)."""
    run = tmp_path / "20240101_run_mixed"
    (run / "pod5").mkdir(parents=True)
    (run / "fast5").mkdir(parents=True)
    make_file(run / "pod5" / "reads.pod5", size=5_000_000)
    make_file(run / "fast5" / "reads.fast5", size=2_000_000)

    result = analyze_run(run)
    assert "pod5" in result["formats"], f"Expected pod5 in {result['formats']}"
    assert "multi_read_fast5" in result["formats"], f"Expected multi_read_fast5 in {result['formats']}"
    print("  PASS: mixed formats (pod5 + multi_read_fast5)")


def test_nested_pod5(tmp_path: Path) -> None:
    """Run with pod5 directory nested several levels deep."""
    run = tmp_path / "20240101_run_nested"
    nested = run / "output" / "basecalling" / "pod5"
    nested.mkdir(parents=True)
    make_file(nested / "reads.pod5", size=5_000_000)

    result = analyze_run(run)
    assert "pod5" in result["formats"], f"Expected pod5, got {result['formats']}"
    print("  PASS: nested pod5 directory")


def test_find_named_subdirs_skips_files(tmp_path: Path) -> None:
    """Verify find_named_subdirs only matches directories, not files."""
    run = tmp_path / "20240101_run_findtest"
    run.mkdir(parents=True)
    # Create a file named "fast5" (not a directory)
    make_file(run / "fast5", size=100)
    # Create a directory named "fast5_pass"
    (run / "fast5_pass").mkdir()

    exact = find_named_subdirs(run, name="fast5")
    prefix = find_named_subdirs(run, prefix="fast5_")
    assert len(exact) == 0, f"Should not match file named fast5, got {exact}"
    assert len(prefix) == 1, f"Should match fast5_pass directory, got {prefix}"
    print("  PASS: find_named_subdirs skips files")


def test_empty_fast5_dir(tmp_path: Path) -> None:
    """Run with empty fast5/ directory -- should report fast5_unknown."""
    run = tmp_path / "20240101_run_empty_f5"
    (run / "fast5").mkdir(parents=True)

    result = analyze_run(run)
    assert "fast5_unknown" in result["formats"], f"Expected fast5_unknown, got {result['formats']}"
    print("  PASS: empty fast5 directory")


def test_fast5_with_numeric_subdirs(tmp_path: Path) -> None:
    """Run with fast5/ containing numeric subdirs with single-read files."""
    run = tmp_path / "20240101_run_f5_numeric"
    fast5_dir = run / "fast5"
    for subdir in range(3):
        d = fast5_dir / str(subdir)
        d.mkdir(parents=True)
        for i in range(5):
            make_file(d / f"read_{i:04d}.fast5", size=30_000)

    result = analyze_run(run)
    assert "single_read_fast5" in result["formats"], f"Expected single_read_fast5, got {result['formats']}"
    assert "numeric subdirs" in result["details"]["single_read_fast5"].get("note", "")
    print("  PASS: fast5/ with numeric subdirs")


def test_many_files_in_root_performance(tmp_path: Path) -> None:
    """Simulate a directory with many small files to verify the scan is not slow."""
    import time
    run = tmp_path / "20240101_run_perf"
    run.mkdir(parents=True)
    # Create 1000 small files (smaller than real 100k but enough to test the path)
    for i in range(1000):
        make_file(run / f"read_{i:06d}.fast5", size=100)

    start = time.monotonic()
    result = analyze_run(run)
    elapsed = time.monotonic() - start

    assert "single_read_fast5" in result["formats"], f"Expected single_read_fast5, got {result['formats']}"
    # Should complete in well under 5 seconds even with 1000 files
    assert elapsed < 5.0, f"Took too long: {elapsed:.2f}s"
    print(f"  PASS: 1000 files in root ({elapsed:.3f}s)")


def test_diagnose_unknown_with_deep_fast5(tmp_path: Path) -> None:
    """Unknown run with .fast5 files buried deep in non-standard subdirs."""
    run = tmp_path / "20240101_run_deep"
    # Create a non-standard directory structure
    deep = run / "custom_output" / "data" / "reads"
    deep.mkdir(parents=True)
    make_file(deep / "read_001.fast5", size=30_000)

    result = analyze_run(run)
    # This should be classified as unknown because the fast5 isn't in expected locations
    assert "unknown" in result["formats"], f"Expected unknown, got {result['formats']}"
    diag = result["details"]["unknown"]
    assert any("deeper in tree" in r for r in diag["reasons"]), f"Should mention deep files: {diag['reasons']}"
    print("  PASS: diagnose_unknown with deep fast5 files")


def test_discover_run_structure(tmp_path: Path) -> None:
    """Verify single-walk discovery finds all format categories."""
    run = tmp_path / "20240101_run_discover"
    (run / "pod5").mkdir(parents=True)
    (run / "pod5_pass").mkdir(parents=True)
    (run / "fast5").mkdir(parents=True)
    (run / "fast5_fail").mkdir(parents=True)
    (run / "fastq_pass").mkdir(parents=True)
    # Nested pod5 inside a subdir
    (run / "output" / "pod5").mkdir(parents=True)
    # Non-matching directory
    (run / "logs").mkdir(parents=True)

    structure = discover_run_structure(run)
    assert len(structure["pod5"]) == 2, f"Expected 2 pod5 dirs, got {structure['pod5']}"
    assert len(structure["pod5_prefix"]) == 1, f"Expected 1 pod5_prefix, got {structure['pod5_prefix']}"
    assert len(structure["fast5"]) == 1, f"Expected 1 fast5 dir, got {structure['fast5']}"
    assert len(structure["fast5_prefix"]) == 1, f"Expected 1 fast5_prefix, got {structure['fast5_prefix']}"
    assert len(structure["fastq_prefix"]) == 1, f"Expected 1 fastq_prefix, got {structure['fastq_prefix']}"
    print("  PASS: discover_run_structure")


def test_fast_count_files_returns_size(tmp_path: Path) -> None:
    """Verify fast_count_files returns (count, size) tuple."""
    d = tmp_path / "20240101_run_count_size" / "data"
    d.mkdir(parents=True)
    make_file(d / "a.pod5", size=1000)
    make_file(d / "b.pod5", size=2000)
    make_file(d / "c.txt", size=500)  # should not match

    count, size = fast_count_files(d, ext=".pod5")
    assert count == 2, f"Expected count=2, got {count}"
    assert size == 3000, f"Expected size=3000, got {size}"
    print("  PASS: fast_count_files returns (count, size)")


def test_fast_count_files_recursive_size(tmp_path: Path) -> None:
    """Verify recursive counting also accumulates size."""
    base = tmp_path / "20240101_run_recsize" / "data"
    (base / "sub1").mkdir(parents=True)
    (base / "sub2").mkdir(parents=True)
    make_file(base / "a.pod5", size=100)
    make_file(base / "sub1" / "b.pod5", size=200)
    make_file(base / "sub2" / "c.pod5", size=300)

    count, size = fast_count_files(base, ext=".pod5", recursive=True)
    assert count == 3, f"Expected count=3, got {count}"
    assert size == 600, f"Expected size=600, got {size}"
    print("  PASS: fast_count_files recursive with size")


def test_format_size(tmp_path: Path) -> None:
    """Verify human-readable size formatting."""
    assert format_size(0) == "0 B"
    assert format_size(512) == "512 B"
    assert format_size(1024) == "1.0 KB"
    assert format_size(1536) == "1.5 KB"
    assert format_size(1_048_576) == "1.0 MB"
    assert format_size(1_073_741_824) == "1.0 GB"
    assert format_size(1_099_511_627_776) == "1.0 TB"
    print("  PASS: format_size")


def test_compute_dir_size(tmp_path: Path) -> None:
    """Verify total directory size calculation."""
    run = tmp_path / "20240101_run_dirsize"
    (run / "sub").mkdir(parents=True)
    make_file(run / "a.txt", size=100)
    make_file(run / "sub" / "b.txt", size=200)

    total = compute_dir_size(run)
    assert total == 300, f"Expected 300, got {total}"
    print("  PASS: compute_dir_size")


def test_quick_mode(tmp_path: Path) -> None:
    """Verify quick mode returns formats without counts or sizes."""
    run = tmp_path / "20240101_run_quick"
    (run / "pod5").mkdir(parents=True)
    make_file(run / "pod5" / "reads.pod5", size=5_000_000)

    result = analyze_run(run, quick=True)
    assert "pod5" in result["formats"], f"Expected pod5, got {result['formats']}"
    detail = result["details"]["pod5"]
    assert "file_count" not in detail, f"Quick mode should not have file_count, got {detail}"
    assert "data_size_bytes" not in detail, f"Quick mode should not have data_size_bytes, got {detail}"
    print("  PASS: quick mode")


def test_size_in_pod5_result(tmp_path: Path) -> None:
    """Verify data_size_bytes appears in pod5 detail when not in quick mode."""
    run = tmp_path / "20240101_run_pod5size"
    (run / "pod5").mkdir(parents=True)
    make_file(run / "pod5" / "reads_001.pod5", size=1000)
    make_file(run / "pod5" / "reads_002.pod5", size=2000)

    result = analyze_run(run)
    assert "pod5" in result["formats"]
    detail = result["details"]["pod5"]
    assert detail["file_count"] == 2
    assert detail["data_size_bytes"] == 3000, f"Expected 3000, got {detail.get('data_size_bytes')}"
    print("  PASS: data_size_bytes in pod5 result")


def test_size_in_fastq_result(tmp_path: Path) -> None:
    """Verify data_size_bytes appears in fastq detail."""
    run = tmp_path / "20240101_run_fqsize"
    (run / "fastq_pass").mkdir(parents=True)
    make_file(run / "fastq_pass" / "reads.fastq.gz", size=500)

    result = analyze_run(run)
    assert "fastq" in result["formats"]
    detail = result["details"]["fastq"]
    assert detail["file_count"] == 1
    assert detail["data_size_bytes"] == 500, f"Expected 500, got {detail.get('data_size_bytes')}"
    print("  PASS: data_size_bytes in fastq result")


def test_discover_skips_data_dirs(tmp_path: Path) -> None:
    """Verify discover_run_structure does not recurse into matched data directories."""
    run = tmp_path / "20240101_run_skipdata"
    fast5_pass = run / "fast5_pass"
    fast5_pass.mkdir(parents=True)
    # Place a nested pod5/ inside fast5_pass -- should NOT be discovered
    (fast5_pass / "pod5").mkdir()
    # Also place a nested fast5/ inside pod5_pass -- should NOT be discovered
    pod5_pass = run / "pod5_pass"
    pod5_pass.mkdir(parents=True)
    (pod5_pass / "fast5").mkdir()

    structure = discover_run_structure(run)
    assert len(structure["pod5"]) == 0, f"Should not find pod5 inside fast5_pass: {structure['pod5']}"
    assert len(structure["fast5"]) == 0, f"Should not find fast5 inside pod5_pass: {structure['fast5']}"
    assert len(structure["fast5_prefix"]) == 1, f"Should find fast5_pass: {structure['fast5_prefix']}"
    assert len(structure["pod5_prefix"]) == 1, f"Should find pod5_pass: {structure['pod5_prefix']}"
    print("  PASS: discover_run_structure skips data dirs")


def test_find_named_subdirs_skips_matched(tmp_path: Path) -> None:
    """Verify find_named_subdirs does not recurse into matched directories."""
    run = tmp_path / "20240101_run_findskip"
    fast5_dir = run / "fast5"
    fast5_dir.mkdir(parents=True)
    # Nest another fast5 inside -- should NOT be found if we skip matched dirs
    (fast5_dir / "fast5").mkdir()

    results = find_named_subdirs(run, name="fast5")
    assert len(results) == 1, f"Should find only top-level fast5, got {results}"
    print("  PASS: find_named_subdirs skips matched dirs")


def test_has_file_with_ext_budget(tmp_path: Path) -> None:
    """Verify _has_file_with_ext stops searching after budget is exhausted."""
    d = tmp_path / "20240101_run_budget"
    d.mkdir(parents=True)
    # Create 20 non-matching files
    for i in range(20):
        make_file(d / f"file_{i:04d}.txt", size=10)
    # Place a matching file that would only be found after budget runs out
    make_file(d / "zzz_last.fast5", size=10)

    # With budget of 5, should not find the .fast5 file (it comes after many .txt files)
    found = _has_file_with_ext(d, ".fast5", max_entries=5)
    assert not found, "Should not find .fast5 with budget=5"

    # With default budget, should find it
    found = _has_file_with_ext(d, ".fast5")
    assert found, "Should find .fast5 with default budget"
    print("  PASS: _has_file_with_ext budget")


def test_diagnose_unknown_capped(tmp_path: Path) -> None:
    """Verify diagnose_unknown caps its scan and reports truncation."""
    run = tmp_path / "20240101_run_diagcap"
    run.mkdir(parents=True)
    # Create more files than the diagnostic cap
    for i in range(5_100):
        make_file(run / f"data_{i:06d}.bin", size=10)

    diag = diagnose_unknown(run)
    assert any("Sampling stopped" in r for r in diag["reasons"]), \
        f"Should mention sampling cap: {diag['reasons']}"
    print("  PASS: diagnose_unknown capped")


def test_data_size_summing(tmp_path: Path) -> None:
    """Verify that total_size_bytes is the sum of data_size_bytes from details."""
    run = tmp_path / "20240101_run_datasum"
    (run / "pod5").mkdir(parents=True)
    (run / "fastq_pass").mkdir(parents=True)
    make_file(run / "pod5" / "reads.pod5", size=1000)
    make_file(run / "fastq_pass" / "reads.fastq.gz", size=500)

    result = analyze_run(run)
    expected = sum(
        detail.get("data_size_bytes", 0)
        for detail in result["details"].values()
    )
    assert expected == 1500, f"Expected 1500, got {expected}"
    print("  PASS: data_size_bytes summing")


def test_estimate_dir_size_exact(tmp_path: Path) -> None:
    """Few files -- count and size should be exact (no estimation)."""
    d = tmp_path / "20240101_run_est_exact" / "data"
    d.mkdir(parents=True)
    make_file(d / "a.fast5", size=100)
    make_file(d / "b.fast5", size=200)
    make_file(d / "c.txt", size=999)  # should not match

    count, size, is_estimated = estimate_dir_size(d, ".fast5")
    assert count == 2, f"Expected count=2, got {count}"
    assert size == 300, f"Expected size=300, got {size}"
    assert not is_estimated, "Should be exact, not estimated"
    print("  PASS: estimate_dir_size exact")


def test_estimate_dir_size_estimated(tmp_path: Path) -> None:
    """More files than sample_size -- size should be estimated."""
    d = tmp_path / "20240101_run_est_sampled" / "data"
    d.mkdir(parents=True)
    file_size = 500
    num_files = 20
    for i in range(num_files):
        make_file(d / f"read_{i:04d}.fast5", size=file_size)

    count, size, is_estimated = estimate_dir_size(d, ".fast5", sample_size=5)
    assert count == num_files, f"Expected count={num_files}, got {count}"
    assert is_estimated, "Should be estimated with sample_size=5"
    # All files are the same size, so the estimate should be exact
    assert size == num_files * file_size, f"Expected {num_files * file_size}, got {size}"
    print("  PASS: estimate_dir_size estimated")


def test_estimate_dir_size_recursive(tmp_path: Path) -> None:
    """Files in subdirectories -- recursive counting."""
    base = tmp_path / "20240101_run_est_rec" / "data"
    (base / "sub1").mkdir(parents=True)
    (base / "sub2").mkdir(parents=True)
    make_file(base / "a.fast5", size=100)
    make_file(base / "sub1" / "b.fast5", size=200)
    make_file(base / "sub2" / "c.fast5", size=300)

    count, size, is_estimated = estimate_dir_size(base, ".fast5", recursive=True)
    assert count == 3, f"Expected count=3, got {count}"
    assert size == 600, f"Expected size=600, got {size}"
    assert not is_estimated, "Should be exact with only 3 files"
    print("  PASS: estimate_dir_size recursive")


def test_single_read_fast5_has_size(tmp_path: Path) -> None:
    """Single-read fast5 in fast5/ should report data_size_bytes."""
    run = tmp_path / "20240101_run_sf5_size"
    (run / "fast5").mkdir(parents=True)
    for i in range(10):
        make_file(run / "fast5" / f"read_{i:04d}.fast5", size=30_000)

    result = analyze_run(run)
    assert "single_read_fast5" in result["formats"]
    detail = result["details"]["single_read_fast5"]
    assert "data_size_bytes" in detail, f"Missing data_size_bytes: {detail}"
    assert detail["file_count"] == 10, f"Expected 10, got {detail.get('file_count')}"
    assert detail["data_size_bytes"] == 300_000, f"Expected 300000, got {detail['data_size_bytes']}"
    print("  PASS: single_read_fast5 in fast5/ has size")


def test_single_read_fast5_root_has_size(tmp_path: Path) -> None:
    """Single-read fast5 in root should report data_size_bytes."""
    run = tmp_path / "20240101_run_sf5_root_size"
    run.mkdir(parents=True)
    for i in range(10):
        make_file(run / f"read_{i:04d}.fast5", size=30_000)

    result = analyze_run(run)
    assert "single_read_fast5" in result["formats"]
    detail = result["details"]["single_read_fast5"]
    assert "data_size_bytes" in detail, f"Missing data_size_bytes: {detail}"
    assert detail["file_count"] == 10, f"Expected 10, got {detail.get('file_count')}"
    assert detail["data_size_bytes"] == 300_000, f"Expected 300000, got {detail['data_size_bytes']}"
    print("  PASS: single_read_fast5 in root has size")


def test_write_stats_tsv_basic(tmp_path: Path) -> None:
    """Verify TSV output contains correct header and row values."""
    all_runs = {
        "20240101_run_pod5": {
            "formats": ["pod5"],
            "details": {
                "pod5": {
                    "file_count": 10,
                    "data_size_bytes": 50000000,
                    "directories": ["/data/run/pod5"],
                    "note": "All good",
                }
            },
        }
    }
    tsv_path = str(tmp_path / "stats.tsv")
    write_stats_tsv(all_runs, tsv_path)

    with open(tsv_path) as fh:
        lines = fh.readlines()
    assert len(lines) == 2, f"Expected 2 lines (header + 1 row), got {len(lines)}"
    header = lines[0].strip().split("\t")
    assert header == ["run_name", "format", "file_count", "data_size_bytes",
                      "size_estimated", "directories", "notes"]
    row = lines[1].strip().split("\t")
    assert row[0] == "20240101_run_pod5"
    assert row[1] == "pod5"
    assert row[2] == "10"
    assert row[3] == "50000000"
    assert row[5] == "/data/run/pod5"
    assert "All good" in row[6]
    print("  PASS: write_stats_tsv basic")


def test_write_stats_tsv_multi_format_run(tmp_path: Path) -> None:
    """Run with pod5 + fastq should produce two rows."""
    all_runs = {
        "20240101_run_mixed": {
            "formats": ["pod5", "fastq"],
            "details": {
                "pod5": {
                    "file_count": 5,
                    "data_size_bytes": 1000,
                    "directories": ["/data/run/pod5"],
                },
                "fastq": {
                    "file_count": 3,
                    "data_size_bytes": 500,
                    "directories": ["/data/run/fastq_pass"],
                },
            },
        }
    }
    tsv_path = str(tmp_path / "stats_multi.tsv")
    write_stats_tsv(all_runs, tsv_path)

    with open(tsv_path) as fh:
        lines = fh.readlines()
    assert len(lines) == 3, f"Expected 3 lines (header + 2 rows), got {len(lines)}"
    row1 = lines[1].strip().split("\t")
    row2 = lines[2].strip().split("\t")
    assert row1[1] == "pod5"
    assert row2[1] == "fastq"
    print("  PASS: write_stats_tsv multi-format run")


def test_generate_conversion_single_to_pod5(tmp_path: Path) -> None:
    """single_read_fast5 -> pod5 should produce a two-step script."""
    runs = {
        "20240101_run_single": {
            "formats": ["single_read_fast5"],
            "details": {
                "single_read_fast5": {
                    "directories": ["/data/run/fast5"],
                }
            },
        }
    }
    script = generate_conversion_script(runs, "pod5")
    assert "single_to_multi_fast5" in script, "Should contain single_to_multi_fast5 step"
    assert "pod5 convert fast5" in script, "Should contain pod5 convert step"
    assert "multi_fast5_tmp" in script, "Should use multi_fast5_tmp intermediate dir"
    assert "--threads 1" in script, "Should use --threads 1"
    assert "two steps" in script, "Should mention two steps"
    print("  PASS: generate_conversion_script single_read_fast5 -> pod5")


def test_generate_conversion_multi_to_pod5(tmp_path: Path) -> None:
    """multi_read_fast5 -> pod5 should produce a direct conversion."""
    runs = {
        "20240101_run_multi": {
            "formats": ["multi_read_fast5"],
            "details": {
                "multi_read_fast5": {
                    "directories": ["/data/run/fast5"],
                }
            },
        }
    }
    script = generate_conversion_script(runs, "pod5")
    assert "pod5 convert fast5" in script, "Should contain pod5 convert"
    assert "--threads 1" in script, "Should use --threads 1"
    assert "--recursive" in script, "Should use --recursive"
    assert "single_to_multi_fast5" not in script, "Should NOT contain single_to_multi_fast5"
    print("  PASS: generate_conversion_script multi_read_fast5 -> pod5")


def test_print_conversion_help_single_fast5(tmp_path: Path) -> None:
    """Conversion help for single_read_fast5 should mention two steps."""
    captured = io.StringIO()
    old_stdout = sys.stdout
    sys.stdout = captured
    try:
        print_conversion_help("single_read_fast5")
    finally:
        sys.stdout = old_stdout
    output = captured.getvalue()
    assert "two steps" in output.lower(), f"Should mention 'two steps', got: {output}"
    assert "single_to_multi_fast5" in output, f"Should mention single_to_multi_fast5, got: {output}"
    assert "pod5 convert fast5" in output, f"Should mention pod5 convert fast5, got: {output}"
    print("  PASS: print_conversion_help single_read_fast5")


def main():
    print("Running format detection tests...\n")
    passed = 0
    failed = 0

    tests = [
        test_pod5_subdirectory,
        test_pod5_variant_dirs,
        test_multi_read_fast5,
        test_single_read_fast5_in_fast5_dir,
        test_single_read_fast5_in_root,
        test_single_read_fast5_numeric_subdirs,
        test_fast5_variant_dirs,
        test_compressed_archives,
        test_empty_directory,
        test_fastq_directory,
        test_mixed_formats,
        test_nested_pod5,
        test_find_named_subdirs_skips_files,
        test_empty_fast5_dir,
        test_fast5_with_numeric_subdirs,
        test_many_files_in_root_performance,
        test_diagnose_unknown_with_deep_fast5,
        test_discover_run_structure,
        test_fast_count_files_returns_size,
        test_fast_count_files_recursive_size,
        test_format_size,
        test_compute_dir_size,
        test_quick_mode,
        test_size_in_pod5_result,
        test_size_in_fastq_result,
        test_discover_skips_data_dirs,
        test_find_named_subdirs_skips_matched,
        test_has_file_with_ext_budget,
        test_diagnose_unknown_capped,
        test_data_size_summing,
        test_estimate_dir_size_exact,
        test_estimate_dir_size_estimated,
        test_estimate_dir_size_recursive,
        test_single_read_fast5_has_size,
        test_single_read_fast5_root_has_size,
        test_write_stats_tsv_basic,
        test_write_stats_tsv_multi_format_run,
        test_generate_conversion_single_to_pod5,
        test_generate_conversion_multi_to_pod5,
        test_print_conversion_help_single_fast5,
    ]

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        for test_fn in tests:
            try:
                test_fn(tmp_path)
                passed += 1
            except Exception as e:
                print(f"  FAIL: {test_fn.__name__}: {e}")
                failed += 1

    print(f"\n{passed} passed, {failed} failed out of {passed + failed} tests")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
