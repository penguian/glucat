#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_license_headers.py:
Verify that source files contain the GluCat LGPL header.
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

import os
import sys

IGNORED_DIRS = {
    ".git",
    ".pytest_cache",
    ".deps",
    ".libs",
    ".venvs",
    "venv",
    "venvs",
    "build",
    "dist",
    "__pycache__",
    "node_modules",
    "artifacts",
}
EXCLUDED_FILES = {"config.h", "glucat_config.h", "PyClical.cpp"}
VALID_EXTENSIONS = {".py", ".h", ".cpp", ".hpp", ".pyx", ".pxd"}
LICENSE_SUBSTRINGS = ["GNU Lesser General Public License", "CC BY-SA 3.0"]


def collect_files(root_dir="."):
    """
    Recursively collect source files with valid extensions, excluding hidden/venv dirs.
    """
    collected = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        dirnames[:] = [
            d
            for d in dirnames
            if d not in IGNORED_DIRS and not d.startswith(".")
        ]
        for filename in filenames:
            if filename in EXCLUDED_FILES:
                continue
            ext = os.path.splitext(filename)[1]
            if ext in VALID_EXTENSIONS:
                collected.append(os.path.join(dirpath, filename))
    return collected


def check_file(filepath):
    """
    Check if the given file contains valid license header substrings.
    """
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            lines = [f.readline() for _ in range(50)]
        content = "".join(lines)

        normalized_path = os.path.normpath(filepath).replace("\\", "/")
        exempt_suffix = "/pyclical/demos/plotting/plotting_demo_dialog.py"
        is_exempt = (
            normalized_path == "pyclical/demos/plotting/plotting_demo_dialog.py"
            or normalized_path.endswith(exempt_suffix)
        )

        if is_exempt:
            if "BSD" + " Style" in content:
                return True
            print(
                f"Error: Missing BSD license header in {filepath}",
                file=sys.stderr,
            )
            return False

        if "BSD" + " Style" in content:
            msg = (
                "Error: BSD"
                f" Style license header found in non-exempt file {filepath}"
            )
            print(msg, file=sys.stderr)
            return False

        if not any(sub in content for sub in LICENSE_SUBSTRINGS):
            msg = (
                "Error: Missing LGPL or CC BY-SA license header "
                f"in {filepath}"
            )
            print(msg, file=sys.stderr)
            return False
        return True
    except (OSError, UnicodeDecodeError) as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return False


def main():
    """
    Main entry point to iterate over arguments and verify license headers.
    """
    if len(sys.argv) >= 2:
        files_to_check = sys.argv[1:]
    else:
        files_to_check = collect_files(".")

    failed = False
    for filepath in files_to_check:
        if not os.path.isfile(filepath):
            continue
        ext = os.path.splitext(filepath)[1]
        if ext not in VALID_EXTENSIONS:
            continue
        if not check_file(filepath):
            failed = True

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
