#!/bin/bash
set -e

# Move to the project root directory
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

# Parse arguments
BACKEND="both"
for i in "$@"; do
  case $i in
    --backend=*)
      BACKEND="${i#*=}"
      shift
      ;;
    *)
      ;;
  esac
done

echo "=== Backend selected: $BACKEND ==="

# Check for tools
if ! command -v llvm-profdata &> /dev/null || ! command -v llvm-cov &> /dev/null; then
    echo "Error: llvm-profdata or llvm-cov not found. Please install LLVM tools."
    exit 1
fi

echo "=== Bootstrapping ==="
make -f admin/Makefile.common bootstrap || true

# Prepare for coverage
export CXX="clang++"
export CXXFLAGS="-O0 -g -fprofile-instr-generate -fcoverage-mapping"

# Clean previous
find . -name "*.profraw" -delete
find . -name "*.profdata" -delete
rm -rf test_doctest/test_doctest_*

run_coverage() {
    local name=$1
    local flags=$2

    echo "=== Running coverage for $name ==="
    ./configure --disable-pyclical --with-doctest $flags
    make clean || true

    # We set LLVM_PROFILE_FILE as an absolute path to avoid directory confusion
    export LLVM_PROFILE_FILE="$(pwd)/${name}.profraw"
    make -C test_doctest check -j$(( $(nproc) / 2 ))

    # Backup binary for llvm-cov
    cp test_doctest/test_doctest test_doctest/test_doctest_${name}
}

if [[ "$BACKEND" == "eigen" || "$BACKEND" == "both" ]]; then
    run_coverage "eigen" ""
fi

if [[ "$BACKEND" == "arma" || "$BACKEND" == "both" ]]; then
    # Check if armadillo is available
    if [ -d "/usr/include/armadillo" ] || [ -f "/usr/include/armadillo" ] || [ -d "/usr/local/include/armadillo" ]; then
        run_coverage "arma" "--with-armadillo"
    else
        echo "Warning: Armadillo not found. Skipping arma backend."
        if [[ "$BACKEND" == "arma" ]]; then exit 1; fi
    fi
fi

echo "=== Merging profiling data ==="
profraw_files=$(find . -name "*.profraw")
if [ -n "$profraw_files" ]; then
    llvm-profdata merge -sparse $profraw_files -o glucat_doctest.profdata
else
    echo "No .profraw files found. Exiting."
    exit 1
fi

echo "=== Generating llvm-cov HTML report ==="
OUTPUT_DIR="test_coverage/results/coverage_html_doctest"
mkdir -p "$OUTPUT_DIR"

# Collect all binaries
OBJECT_ARGS=""
for f in test_doctest/test_doctest_*; do
    if [ -f "$f" ]; then
        if [ -z "$OBJECT_ARGS" ]; then
            OBJECT_ARGS="$f"
        else
            OBJECT_ARGS="$OBJECT_ARGS -object $f"
        fi
    fi
done

llvm-cov show -format=html -output-dir="$OUTPUT_DIR" \
    -instr-profile=glucat_doctest.profdata \
    $OBJECT_ARGS \
    -- glucat/

llvm-cov report -instr-profile=glucat_doctest.profdata $OBJECT_ARGS -- glucat/

echo "Coverage report generated in $OUTPUT_DIR/index.html"
