#!/bin/bash
set -e

# Move to the project root directory
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "=== Configuring for PyClical Coverage ==="
export CXX="clang++"
export CXXFLAGS="-O0 -g -fprofile-instr-generate -fcoverage-mapping"
export LDFLAGS="-fprofile-instr-generate"
export USER_LDFLAGS="-fprofile-instr-generate"

# Make sure we have Cython and Python headers
./configure --enable-pyclical

echo "=== Cleaning and Building ==="
make clean || true
rm -f $(find . -name "*.profraw") glucat.profdata
make -j$(nproc)

echo "=== Running PyClical tests ==="
export LLVM_PROFILE_FILE="$(pwd)/%p.profraw"
# Run PyClical tests
cd pyclical
python3 test.py
cd ..

echo "=== Merging profiling data ==="
profraw_files=$(find . -name "*.profraw")
if [ -n "$profraw_files" ]; then
    llvm-profdata merge -sparse $profraw_files -o pyclical.profdata
else
    echo "No .profraw files found. Exiting."
    exit 1
fi

echo "=== Generating llvm-cov report for PyClical ==="
SO_FILE=$(find pyclical -name "PyClical.*.so" | head -n 1)

if [ -z "$SO_FILE" ]; then
    echo "PyClical extension not found."
    exit 1
fi

llvm-cov report -instr-profile=pyclical.profdata "$SO_FILE" -- glucat/
