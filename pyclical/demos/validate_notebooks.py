#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# validate_notebooks.py: Build and validate Jupyter notebook files.
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
#

import glob
import json
import os
import subprocess
import sys
import nbformat
from nbformat.validator import validate, NotebookValidationError


def main():
    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 1. Run the build script
    print("Building notebooks...")
    try:
        subprocess.run(
            [sys.executable, "build_pyclical_notebooks.py"],
            cwd=script_dir,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Failed to build notebooks: {e}", file=sys.stderr)
        print(f"Stdout:\n{e.stdout}", file=sys.stderr)
        print(f"Stderr:\n{e.stderr}", file=sys.stderr)
        return 1

    # 2. Find all generated .ipynb files
    notebook_files = glob.glob(os.path.join(script_dir, "*.ipynb"))
    if not notebook_files:
        print("No notebook files found to validate.", file=sys.stderr)
        return 1

    # 3. Validate JSON syntax and notebook format
    has_errors = False
    for nb_path in sorted(notebook_files):
        nb_name = os.path.basename(nb_path)
        print(f"Validating {nb_name}...")

        # Check JSON syntax
        try:
            with open(nb_path, "r", encoding="utf-8") as f:
                json.load(f)
        except json.JSONDecodeError as e:
            print(
                f"  [ERROR] Invalid JSON syntax in {nb_name}: {e}",
                file=sys.stderr,
            )
            has_errors = True
            continue

        # Check notebook schema
        try:
            # Re-read with nbformat to validate schema
            with open(nb_path, "r", encoding="utf-8") as f:
                nb = nbformat.read(f, as_version=4)
            validate(nb)
            print(f"  {nb_name} is valid.")
        except NotebookValidationError as e:
            print(
                f"  [ERROR] Schema validation failed in {nb_name}: {e}",
                file=sys.stderr,
            )
            has_errors = True
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(
                f"  [ERROR] Unexpected error validating {nb_name}: {e}",
                file=sys.stderr,
            )
            has_errors = True

    return 1 if has_errors else 0


if __name__ == "__main__":
    sys.exit(main())
