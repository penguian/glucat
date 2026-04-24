#!/bin/bash
set -e

# Move to the project root directory
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

DOCTEST_DIR="/home/leopardi/src-downloaded/doctest/doctest/doctest"

echo "=== Bootstrapping ==="
make -f admin/Makefile.common bootstrap || true

echo "=== Configuring for Clang Doctest Coverage ==="
export CXX="clang++"
export CXXFLAGS="-O0 -g -fprofile-instr-generate -fcoverage-mapping"
export LDFLAGS="-fprofile-instr-generate"

# Configure with doctest and without pyclical
./configure --disable-pyclical --with-doctest=$DOCTEST_DIR

echo "=== Cleaning previous builds ==="
make clean || true
rm -f *.profraw test_doctest/*.profraw test_doctest/*.profdata

echo "=== Building and running doctest ==="
export LLVM_PROFILE_FILE="test_doctest.profraw"
make -C test_doctest -j$(nproc)
./test_doctest/test_doctest

echo "=== Merging profiling data ==="
llvm-profdata merge -sparse test_doctest.profraw -o glucat_doctest.profdata

echo "=== Generating llvm-cov HTML report ==="
mkdir -p coverage_html_doctest
llvm-cov show -format=html -output-dir=coverage_html_doctest \
    -instr-profile=glucat_doctest.profdata \
    test_doctest/test_doctest \
    -- glucat/

llvm-cov report -instr-profile=glucat_doctest.profdata test_doctest/test_doctest -- glucat/

echo "Coverage report generated in coverage_html_doctest/index.html"
