#!/bin/bash
set -e

# Move to the project root directory
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"


echo "=== Bootstrapping ==="
make -f admin/Makefile.common bootstrap || true

echo "=== Configuring for Clang Doctest Coverage ==="
export CXX="clang++"
export CXXFLAGS="-O0 -g -fprofile-instr-generate -fcoverage-mapping"
export LDFLAGS="-fprofile-instr-generate"

# Configure with doctest and without pyclical
./configure --disable-pyclical --with-doctest

echo "=== Cleaning previous builds ==="
make clean || true
find . -name "*.profraw" -delete
find . -name "*.profdata" -delete

echo "=== Building and running doctest ==="
# We set LLVM_PROFILE_FILE as an absolute path to avoid directory confusion
export LLVM_PROFILE_FILE="$(pwd)/test_doctest.profraw"
make -C test_doctest check -j$(nproc)

echo "=== Merging profiling data ==="
llvm-profdata merge -sparse test_doctest.profraw -o glucat_doctest.profdata

echo "=== Generating llvm-cov HTML report ==="
OUTPUT_DIR="test_coverage/results/coverage_html_doctest"
mkdir -p "$OUTPUT_DIR"
llvm-cov show -format=html -output-dir="$OUTPUT_DIR" \
    -instr-profile=glucat_doctest.profdata \
    test_doctest/test_doctest \
    -- glucat/

llvm-cov report -instr-profile=glucat_doctest.profdata test_doctest/test_doctest -- glucat/

echo "Coverage report generated in $OUTPUT_DIR/index.html"
