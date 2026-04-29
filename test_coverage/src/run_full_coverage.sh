#!/bin/bash
set -e

# Move to the project root directory
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

echo "=== Bootstrapping ==="
make -f admin/Makefile.common bootstrap || true

echo "=== Configuring for Full Coverage (Clang) ==="
export CXX="clang++"
export CXXFLAGS="-O0 -g -fprofile-instr-generate -fcoverage-mapping"

# Cleanup previous data
rm -f *.profraw *.profdata
rm -rf test_coverage/results/coverage_html_full
mkdir -p test_coverage/results/coverage_html_full

# Directory to preserve QD-enabled legacy binaries before Armadillo make clean
LEGACY_BIN_DIR="$(pwd)/test_coverage/results/legacy_bins"
rm -rf "$LEGACY_BIN_DIR"
mkdir -p "$LEGACY_BIN_DIR"

TESTS="test00 test01 test02 test03 test04 test05 test06 test07 test08 test09 test10 test11 test12 test13 test14 test15 test16 test17"

# 1. Run Legacy Test Suite (test00-test17)
echo "=== 1. Running Legacy Test Suite ==="
./configure --disable-pyclical --with-doctest --with-qd
make clean || true
# We use a unique name for legacy profile
export LLVM_PROFILE_FILE="$(pwd)/legacy_%p.profraw"
make check -j$(( $(nproc) / 2 ))

# Back up all legacy binaries NOW, before any subsequent make clean wipes them.
# This preserves the QD-enabled instrumented binaries for llvm-cov.
# Without this, the Armadillo pass's make clean destroys the QD-instantiated
# test11 (and other) binaries, causing qd.h to show 0% coverage even though
# the profraw data contains QD execution records.
echo "=== Backing up legacy binaries ==="
for t in $TESTS; do
    if [ -f "${t}/${t}" ]; then
        cp "${t}/${t}" "$LEGACY_BIN_DIR/${t}"
        echo "  Backed up ${t}/${t}"
    fi
done

# 2. Run Doctest Suite (Eigen Backend)
echo "=== 2. Running Doctest (Eigen) ==="
export LLVM_PROFILE_FILE="$(pwd)/doctest_eigen.profraw"
make -C test_doctest check
# Backup binary
cp test_doctest/test_doctest test_doctest/test_doctest_eigen

# 3. Run Doctest Suite (Armadillo Backend)
if [ -d "/usr/include/armadillo" ] || [ -f "/usr/include/armadillo" ] || [ -d "/usr/local/include/armadillo" ]; then
    echo "=== 3. Running Doctest (Armadillo) ==="
    ./configure --disable-pyclical --with-doctest --with-armadillo --with-qd
    # Only clean test_doctest, NOT the entire tree.
    # A top-level 'make clean' would wipe the QD-enabled legacy binaries above;
    # we only need to rebuild test_doctest with the Armadillo backend.
    make -C test_doctest clean || true
    export LLVM_PROFILE_FILE="$(pwd)/doctest_arma.profraw"
    make -C test_doctest check
    # Backup binary
    cp test_doctest/test_doctest test_doctest/test_doctest_arma
else
    echo "Armadillo not found. Skipping Armadillo pass."
fi

echo "=== Merging ALL profiling data ==="
profraw_files=$(find . -name "*.profraw")
llvm-profdata merge -sparse $profraw_files -o glucat_full.profdata

echo "=== Generating Comprehensive Coverage Report ==="
# Collect backed-up legacy binaries (QD-enabled, matching the legacy profraw data).
# Using the backed-up copies ensures llvm-cov can decode QD coverage even though
# the working-tree binaries may have been cleaned/rebuilt since.
BINARIES=""
for t in $TESTS; do
    if [ -f "$LEGACY_BIN_DIR/${t}" ]; then
        BINARIES="$BINARIES -object $LEGACY_BIN_DIR/${t}"
    else
        echo "Warning: backed-up legacy binary for ${t} not found."
    fi
done

# Add doctest binaries (we backed them up)
if [ -f "test_doctest/test_doctest_eigen" ]; then
    BINARIES="$BINARIES -object test_doctest/test_doctest_eigen"
fi
if [ -f "test_doctest/test_doctest_arma" ]; then
    BINARIES="$BINARIES -object test_doctest/test_doctest_arma"
fi

# The first binary is positional, rest use -object
FIRST_OBJ=$(echo $BINARIES | awk '{print $2}')
REST_OBJS=$(echo $BINARIES | cut -d' ' -f3-)

llvm-cov show -format=html -output-dir="test_coverage/results/coverage_html_full" \
    -instr-profile=glucat_full.profdata \
    $FIRST_OBJ $REST_OBJS \
    -- glucat/

llvm-cov report -instr-profile=glucat_full.profdata $FIRST_OBJ $REST_OBJS -- glucat/

echo "Full coverage report generated in test_coverage/results/coverage_html_full/index.html"
