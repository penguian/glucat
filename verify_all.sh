#!/bin/bash
set -e

if [ "$1" = "--coverage" ]; then
  echo "=== Running C++ Block Coverage ==="
  if ! command -v llvm-profdata &> /dev/null || ! command -v llvm-cov &> /dev/null; then
    echo "Error: llvm-profdata or llvm-cov not found. Please install LLVM tools." >&2
    exit 1
  fi
  shift
  bash test_coverage/src/run_clang_doctest_coverage.sh "$@"
  exit 0
fi

if [ ! -f Makefile ]; then
  echo "Error: Makefile not found. Please bootstrap and configure the project first:" >&2
  echo "  make -f admin/Makefile.common bootstrap" >&2
  echo "  ./configure [options]" >&2
  exit 1
fi

# C++ linting & static analysis
make lint-check
make cppcheck

# Primary C++ regression suite
make check-local

# Run PyClical, demos, and benchmarks linting
ruff check pyclical/ pyclical/demos/ benchmarks/
pylint pyclical/ pyclical/demos/ benchmarks/

# Notebook validation
python3 pyclical/demos/validate_notebooks.py
