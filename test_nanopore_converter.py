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
