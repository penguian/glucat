#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# GluCat: Generic library of universal Clifford algebra templates
#
# format_helper.py : Helper script to run clang-format with custom preprocessing
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
import re
import subprocess

def preprocess(content):
    # Transforms:
    # auto blah()
    # { return foo; }
    # into:
    # auto blah()
    # {
    #   return foo;
    # }
    pattern = re.compile(r"\n(\s*)\{\s*(return\s+[^;]+;)\s*\}")
    content = pattern.sub(r"\n\1{\n\1  \2\n\1}", content)
    # Transforms:
    # { } or {}
    # into:
    # {
    # }
    empty_pattern = re.compile(r"\n(\s*)\{\s*\}")
    return empty_pattern.sub(r"\n\1{\n\1}", content)

def postprocess(content):
    # Transforms:
    # auto blah()
    # {
    #   return foo;
    # }
    # back into:
    # auto blah()
    # { return foo; }
    pattern = re.compile(r"\n(\s*)\{\n\s*(return\s+[^;]+;)\n\s*\}")
    content = pattern.sub(r"\n\1{ \2 }", content)
    # Transforms:
    # {
    # }
    # back into:
    # { }
    empty_pattern = re.compile(r"\n(\s*)\{\n\1\}")
    return empty_pattern.sub(r"\n\1{ }", content)


def check_file(filepath):
    if "glucat_config.h" in filepath:
        return True
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read()
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return False

    preprocessed = preprocess(content)

    process = subprocess.Popen(
        [
            "clang-format",
            "-dry-run",
            "-Werror",
            f"--assume-filename={filepath}",
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    stdout, stderr = process.communicate(input=preprocessed)

    if process.returncode != 0:
        # Print compiler-like diagnostics for the file
        print(f"Formatting violations in {filepath}:", file=sys.stderr)
        sys.stderr.write(stderr)
        return False
    return True

def fix_file(filepath):
    if "glucat_config.h" in filepath:
        return True
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read()
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return False


    preprocessed = preprocess(content)

    process = subprocess.Popen(
        ["clang-format", f"--assume-filename={filepath}"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    formatted, stderr = process.communicate(input=preprocessed)

    if process.returncode != 0:
        print(f"clang-format error on {filepath}:", file=sys.stderr)
        sys.stderr.write(stderr)
        return False

    postprocessed = postprocess(formatted)

    if postprocessed != content:
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(postprocessed)
        print(f"Formatted {filepath}")
    return True

def main():
    if len(sys.argv) < 3:
        print("Usage: format_helper.py [check|fix] [files...]", file=sys.stderr)
        sys.exit(1)

    mode = sys.argv[1]
    files = sys.argv[2:]

    success = True
    if mode == "check":
        for filepath in files:
            if not check_file(filepath):
                success = False
    elif mode == "fix":
        for filepath in files:
            if not fix_file(filepath):
                success = False
    else:
        print(f"Unknown mode: {mode}", file=sys.stderr)
        sys.exit(1)

    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
