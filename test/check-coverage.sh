#!/usr/bin/env bash
#
#    check-coverage.sh: Clang source-based coverage gap report (test00-test19 vs test_doctest)
#    copyright            : (C) 2026 by Paul C. Leopardi
#
#    This library is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library.  If not, see <http://www.gnu.org/licenses/>.

set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
package_dir=$(cd "$here/.." && pwd)
output_dir="${1:-${package_dir}/coverage_report}"

COVERAGE_FLAGS="-fprofile-instr-generate -fcoverage-mapping"
CLANG="${CXX:-clang++}"

# Ensure MAKEFLAGS respects the machine limit of -j8
export MAKEFLAGS="-j8"

echo "=== 1. Preparing coverage build in out-of-tree directory ==="
build_dir="${package_dir}/../glucat.coverage"
rm -rf "$build_dir"
cp -a "$package_dir" "$build_dir"

pushd "$build_dir" > /dev/null
  echo "Running ./configure with Clang coverage instrumentation..."
  ./configure CXX="$CLANG" CXXFLAGS="-O1 -g $COVERAGE_FLAGS" --disable-pyclical > /dev/null
  echo "Building test binaries (test00..test19, test_doctest)..."
  make -j8 $(seq -f 'test%02g/test%02g' 0 19 | tr '\n' ' ') test_doctest/test_doctest > /dev/null
popd > /dev/null

mkdir -p "${output_dir}/raw"

echo "=== 2. Running test00..test19 and gathering raw profiles ==="
for i in $(seq -f '%02g' 0 19); do
  echo -n "Running test${i}... "
  LLVM_PROFILE_FILE="${output_dir}/raw/test${i}.profraw" \
    "$build_dir/test${i}/test${i}" > /dev/null 2>&1 || true
  echo "done."
done

echo "=== 3. Running test_doctest and gathering raw profile ==="
LLVM_PROFILE_FILE="${output_dir}/raw/test_doctest.profraw" \
  "$build_dir/test_doctest/test_doctest" > /dev/null 2>&1 || true

echo "=== 4. Merging profile data ==="
llvm-profdata merge -sparse \
  "${output_dir}/raw"/test*.profraw \
  -o "${output_dir}/all_tests.profdata"

llvm-profdata merge -sparse \
  "${output_dir}/raw/test_doctest.profraw" \
  -o "${output_dir}/test_doctest_only.profdata"

echo "=== 5. Generating LCOV coverage exports for glucat/ headers ==="
all_bins=$(seq -f "$build_dir/test%02g/test%02g" 0 19)
doctest_bin="$build_dir/test_doctest/test_doctest"

# Build --object flags for llvm-cov
object_args=()
for b in $all_bins; do
  object_args+=(--object "$b")
done

llvm-cov export \
  --format=lcov \
  --instr-profile="${output_dir}/all_tests.profdata" \
  --object "$doctest_bin" \
  "${object_args[@]}" \
  "$build_dir/glucat/" \
  > "${output_dir}/all_tests.lcov"

llvm-cov export \
  --format=lcov \
  --instr-profile="${output_dir}/test_doctest_only.profdata" \
  --object "$doctest_bin" \
  "$build_dir/glucat/" \
  > "${output_dir}/test_doctest_only.lcov"

echo "=== 6. Computing coverage gap report ==="
python3 "${here}/coverage_gap_report.py" \
  "${output_dir}/all_tests.lcov" \
  "${output_dir}/test_doctest_only.lcov" \
  > "${output_dir}/coverage_gap_report.txt"

echo
cat "${output_dir}/coverage_gap_report.txt"
echo
echo "Coverage gap report saved to ${output_dir}/coverage_gap_report.txt"

# Cleanup build directory
rm -rf "$build_dir"
