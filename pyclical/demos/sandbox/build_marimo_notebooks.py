#!/usr/bin/env python3
"""
Convert generated Jupyter notebooks (.ipynb) into Marimo reactive Python notebooks.
"""
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# build_marimo_notebooks.py: Build Marimo reactive notebooks from Jupyter notebooks.
#
#    copyright            : (C) 2026 by Paul C. Leopardi
#
# Licensed under CC BY-SA 3.0 http://creativecommons.org/licenses/by-sa/3.0/

import glob
import json
import os
import shutil
import subprocess
import sys


def convert_ipynb_to_marimo_fallback(ipynb_path, out_py_path):
    """
    Programmatically convert a Jupyter notebook JSON structure into a Marimo .py app
    when the marimo CLI tool is not installed.
    """
    with open(ipynb_path, "r", encoding="utf-8") as f:
        nb = json.load(f)

    lines = []
    lines.append("import marimo")
    lines.append("")
    lines.append('app = marimo.App(title="PyClical Marimo Sandbox")')
    lines.append("")

    cell_index = 0
    for cell in nb.get("cells", []):
        cell_type = cell.get("cell_type")
        source = cell.get("source", [])
        if isinstance(source, list):
            source_text = "".join(source)
        else:
            source_text = str(source)

        if not source_text.strip():
            continue

        cell_index += 1
        lines.append("@app.cell")
        if cell_type == "markdown":
            # Wrap markdown in marimo.md(...)
            escaped_md = repr(source_text)
            lines.append(f"def cell_{cell_index}(mo):")
            lines.append(f"    return mo.md({escaped_md}),")
        else:
            # Code cell
            indented_code = "\n    ".join(source_text.splitlines())
            lines.append(f"def cell_{cell_index}():")
            if indented_code.strip():
                lines.append(f"    {indented_code}")
            else:
                lines.append("    pass")
            lines.append("    return")
        lines.append("")

    lines.append('if __name__ == "__main__":')
    lines.append("    app.run()")
    lines.append("")

    with open(out_py_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def main():
    """
    Main conversion loop for all .ipynb files in the sandbox directory.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    marimo_cli = shutil.which("marimo")

    notebook_files = glob.glob(os.path.join(script_dir, "*.ipynb"))
    if not notebook_files:
        print("No Jupyter notebooks found to convert.", file=sys.stderr)
        return 1

    print(f"Found {len(notebook_files)} notebooks for Marimo conversion.")
    for nb_path in sorted(notebook_files):
        base_name = os.path.splitext(os.path.basename(nb_path))[0]
        out_py = os.path.join(script_dir, f"{base_name}_marimo.py")
        print(f"Generating {os.path.basename(out_py)}...")

        if marimo_cli:
            try:
                subprocess.run(
                    [marimo_cli, "convert", nb_path, "-o", out_py],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                print("  Used marimo CLI converter.")
                continue
            except subprocess.CalledProcessError:
                print("  CLI conversion failed; falling back to python generator.")

        convert_ipynb_to_marimo_fallback(nb_path, out_py)
        print("  Generated using internal AST/template fallback generator.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
