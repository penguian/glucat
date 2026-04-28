#!/bin/bash
set -e

# Move to the project root directory
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

echo "=== Bootstrapping ==="
./admin/Makefile.common bootstrap || true

echo "=== Configuring for Full Coverage (Clang) ==="
export CXX="clang++"
export CXXFLAGS="-O0 -g -fprofile-instr-generate -fcoverage-mapping"
export LDFLAGS="-fprofile-instr-generate"

# Cleanup previous data
rm -f *.profraw *.profdata
rm -rf test_coverage/results/coverage_html_full
mkdir -p test_coverage/results/coverage_html_full

# 1. Run Legacy Test Suite (test00-test17)
echo "=== 1. Running Legacy Test Suite ==="
./configure --disable-pyclical --with-doctest --with-qd
make clean || true
# We use a unique name for legacy profile
export LLVM_PROFILE_FILE="$(pwd)/legacy_%p.profraw"
make check -j$(( $(nproc) / 2 ))

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
    make clean || true
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
# Collect all binaries
TESTS="test00 test01 test02 test03 test04 test05 test06 test07 test08 test09 test10 test11 test12 test13 test14 test15 test16 test17"
BINARIES=""

# Add legacy binaries
for t in $TESTS; do
    if [ -f "${t}/${t}" ]; then
        BINARIES="$BINARIES -object ${t}/${t}"
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
