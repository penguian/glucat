#!/bin/bash
set -e

# Move to the project root directory
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

echo "=== Bootstrapping ==="
make -f admin/Makefile.common bootstrap || true

echo "=== Configuring for Clang Coverage ==="
export CXX="clang++"
export CXXFLAGS="-O0 -g -fprofile-instr-generate -fcoverage-mapping"
export LDFLAGS="-fprofile-instr-generate"

# Configure without pyclical to focus on C++ library coverage
./configure --disable-pyclical

echo "=== Cleaning previous builds ==="
make clean || true
# Explicitly remove all old profiling files to prevent merging stale data
find . -name "*.profraw" -delete
find . -name "*.profdata" -delete

echo "=== Building and running tests in parallel ==="
export LLVM_PROFILE_FILE="%p.profraw"
# This builds the library and runs all tests (test00-test17) using all cores
make check -j$(( $(nproc) / 2 ))

echo "=== Merging profiling data ==="
# Find all profraw files created by the tests
profraw_files=$(find . -name "*.profraw")
if [ -n "$profraw_files" ]; then
    llvm-profdata merge -sparse $profraw_files -o glucat.profdata
else
    echo "No .profraw files found. Exiting."
    exit 1
fi

echo "=== Generating llvm-cov HTML report ==="
OUTPUT_DIR="test_coverage/results/coverage_html_clang"
mkdir -p "$OUTPUT_DIR"

# We must pass every compiled test binary to llvm-cov so it can read the coverage map
# Even if they are all merged into one profdata, llvm-cov needs the binary for the debug symbols/mapping
TESTS="test00 test01 test02 test03 test04 test05 test06 test07 test08 test09 test10 test11 test12 test13 test14 test15 test16 test17"
OBJECT_ARGS=""
for t in $TESTS; do
    if [ -f "${t}/${t}" ]; then
        OBJECT_ARGS="$OBJECT_ARGS -object ${t}/${t}"
    else
        echo "Warning: binary ${t}/${t} not found."
    fi
done

# First binary is positional; rest use -object; source directory follows --
FIRST_OBJ=""
REST_ARGS=""
for arg in $OBJECT_ARGS; do
    if [ -z "$FIRST_OBJ" ] && [ "$arg" != "-object" ]; then
        FIRST_OBJ="$arg"
    elif [ "$arg" != "-object" ]; then
        REST_ARGS="$REST_ARGS -object $arg"
    fi
done

llvm-cov show -format=html -output-dir="$OUTPUT_DIR" \
    -instr-profile=glucat.profdata \
    $FIRST_OBJ $REST_ARGS \
    -- glucat/

llvm-cov report -instr-profile=glucat.profdata $FIRST_OBJ $REST_ARGS -- glucat/

echo "Coverage report generated in $OUTPUT_DIR/index.html"
