#!/bin/bash
set -e

# Run PyClical, demos, and benchmarks linting
ruff check pyclical/ pyclical/demos/ benchmarks/
pylint pyclical/ pyclical/demos/ benchmarks/

# Notebook validation (to be enabled in Stage 4)
python3 pyclical/demos/validate_notebooks.py

# C++ linting & static analysis
make lint-check
make cppcheck

# Primary C++ regression suite
make check-local
