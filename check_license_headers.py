#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# check_license_headers.py:
# Verify that source files contain the GluCat LGPL header.
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
import os

LICENSE_SUBSTRINGS = ["GNU Lesser General Public License", "CC BY-SA 3.0"]

def check_file(filepath):
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            lines = [f.readline() for _ in range(50)]
        content = "".join(lines)

        normalized_path = os.path.normpath(filepath).replace("\\", "/")
        exempt_suffix = "/pyclical/demos/plotting_demo_dialog.py"
        is_exempt = (
            normalized_path == "pyclical/demos/plotting_demo_dialog.py"
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
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return False

def main():
    if len(sys.argv) < 2:
        sys.exit(0)

    failed = False
    for filepath in sys.argv[1:]:
        if not os.path.isfile(filepath):
            continue
        ext = os.path.splitext(filepath)[1]
        if ext not in [".py", ".h", ".cpp", ".hpp", ".pyx", ".pxd"]:
            continue
        if not check_file(filepath):
            failed = True

    if failed:
        sys.exit(1)

if __name__ == "__main__":
    main()
