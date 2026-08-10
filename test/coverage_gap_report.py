#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
coverage_gap_report.py:
Compare two LCOV coverage files and report lines covered in baseline (test00-test19)
that are NOT covered in comparison target (test_doctest).
"""
#
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

import sys


def parse_lcov(filepath: str) -> dict[str, dict[int, int]]:
    """Parse LCOV file into a dict mapping filename -> {line_number: hit_count}."""
    file_coverage: dict[str, dict[int, int]] = {}
    current_file = None

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith("SF:"):
                current_file = line[3:].strip()
                if current_file not in file_coverage:
                    file_coverage[current_file] = {}
            elif line.startswith("DA:") and current_file:
                parts = line[3:].split(",")
                if len(parts) >= 2:
                    lineno = int(parts[0])
                    hits = int(parts[1])
                    file_coverage[current_file][lineno] = hits
            elif line == "end_of_record":
                current_file = None

    return file_coverage


def main() -> int:
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <all_tests.lcov> <test_doctest_only.lcov>", file=sys.stderr)
        return 1

    all_tests_cov = parse_lcov(sys.argv[1])
    doctest_cov = parse_lcov(sys.argv[2])

    uncovered_count = 0

    for fname in sorted(all_tests_cov.keys()):
        lines_all = all_tests_cov[fname]
        lines_doc = doctest_cov.get(fname, {})

        file_uncovered = []
        for lineno in sorted(lines_all.keys()):
            hits_all = lines_all[lineno]
            hits_doc = lines_doc.get(lineno, 0)
            if hits_all > 0 and hits_doc == 0:
                file_uncovered.append(lineno)

        if file_uncovered:
            print(f"File: {fname}")
            for lineno in file_uncovered:
                print(f"  {fname}:{lineno}  [UNCOVERED BY DOCTEST]")
                uncovered_count += 1
            print()

    print(f"Total lines covered by runtime tests but uncovered by doctest: {uncovered_count}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
