"""
Utilities for generating and formatting PyClical IPython notebook cells.
"""
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_notebook_utils.py:
#   This file contains utilities for use with PyClical notebooks.
#
#    copyright            : (C) 2013-2014 by Paul C. Leopardi
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

import json as js
import textwrap
from pyclical_tutorial_utils import *


def print_cell_markdown(source):
    """
Print a JSON-formatted markdown cell.
"""
    print(
        js.dumps(
            {"cell_type": "markdown", "metadata": {}, "source": [source]},
            indent=1,
        ),
        end="",
    )


def print_cell_code(code_input, execution_count):
    """
Print a JSON-formatted code cell.
"""
    print(
        js.dumps(
            {
                "cell_type": "code",
                "source": [code_input],
                "metadata": {},
                "outputs": [],
                "execution_count": execution_count,
            },
            indent=1,
        ),
        end="",
    )


def print_metadata_name(name):
    """
Print metadata section name for notebook.
"""
    print(js.dumps({"name": name}, indent=1), end="")


class NotebookContext(InteractionContext):
    """
Execution context for compiling tutorials into Jupyter/IPython notebooks.
"""

    def __init__(self, dictionary):
        super().__init__(dictionary)
        self.execution_count = 1

    # pylint: disable=unused-argument
    def print_notebook_header(self, notebook_title):
        """
Print notebook JSON header format.
"""
        print("{")
        print(' "nbformat":', 4, ",")
        print(' "nbformat_minor":', 2, ",")
        print(' "metadata": {')
        print('   "kernelspec": {')
        print('    "display_name": "Python 3",')
        print('    "language": "python",')
        print('    "name": "python3"')
        print("   },")
        print('   "language_info": {')
        print('    "codemirror_mode": {')
        print('     "name": "ipython",')
        print('     "version": 3')
        print("    },")
        print('    "file_extension": ".py",')
        print('    "mimetype": "text/x-python",')
        print('    "name": "python"')
        print("   }")
        print("  },")
        print(' "cells": [')

    def print_notebook_footer(self):
        """
Print notebook JSON footer format.
"""
        print_cell_markdown(" ")
        print("")
        print(" ]")
        print("}")

    def pause(self):
        """
No-op for notebooks.
"""

    def print_head(self, output_str, indent=""):
        """
Print notebook header cell.
"""
        print_cell_markdown("# " + output_str)
        print(",")

    def print_fill(self, output_str, indent="    "):
        """
Print notebook markdown content cell.
"""
        print_cell_markdown(output_str)
        print(",")

    def print_line(self):
        """
No-op for notebook line breaks.
"""

    def print_exec(self, command_str):
        """
Print notebook code cell.
"""
        print_cell_code(command_str, self.execution_count)
        print(",")
        self.execution_count += 1

    def display_image(self, url):
        """
        Print notebook code cells for importing Image and rendering path or URL.
        """
        self.print_exec("from IPython.display import Image")
        if url.startswith("http://") or url.startswith("https://"):
            self.print_exec(f"Image(url='{url}')")
        else:
            self.print_exec(f"Image(filename='{url}')")

    def input_exec(self, prompt, sandbox):
        """
No-op for notebook interactive input.
"""

    def check_exec(self, prompt, var_name, value_str):
        """
Print exercise solution block in notebook format.
"""
        self.print_fill("Exercise: Enter a Python statement to " + prompt)
        self.print_exec("")
        self.print_fill(
            "Here is one way to do this, and then print the result:"
        )
        self.print_exec(
            var_name + " = " + value_str + "; print(" + var_name + ")"
        )

    def input_eval(self, prompt):
        """
No-op for notebook expression evaluation.
"""

    def check_eval(self, prompt, value_str, command_str):
        """
Print exercise evaluation block in notebook format.
"""
        self.print_fill("Exercise: Enter a Python expression to " + prompt)
        self.print_exec("")
        self.print_fill("Here is one way to use such an expression:")
        command_str = command_str.format(value_str)
        self.print_exec(command_str)


# Backward-compatibility alias
notebook_context = NotebookContext  # pylint: disable=invalid-name


def _get_marimo_version():
    """Return the installed marimo version string or a safe fallback."""
    try:
        import marimo  # pylint: disable=import-outside-toplevel

        return getattr(marimo, "__version__", "0.23.5")
    except ImportError:
        return "0.23.5"


def _format_markdown_cell(text):
    """Format Markdown text as a Marimo mo.md(...) call matching marimo check."""
    cleaned = textwrap.dedent(
        "\n".join(line.rstrip() for line in text.splitlines())
    ).strip()
    if (
        '"""' not in cleaned
        and '""' not in cleaned
        and not cleaned.endswith(("\\", '"'))
        and "\r" not in cleaned
    ):
        indented = "\n".join(
            f"    {line}" if line else "" for line in cleaned.splitlines()
        )
        return f'    mo.md(r"""\n{indented}\n    """)'
    return f"    mo.md({repr(cleaned)})"


def _emit_marimo_item(lines, item):
    """Append Marimo @app.cell definitions for a single recorded notebook item."""
    kind = item[0]
    if kind == "markdown":
        lines.extend(
            [
                "",
                "",
                "@app.cell(hide_code=True)",
                "def _(mo):",
                _format_markdown_cell(item[1]),
            ]
        )
        return

    if kind == "image":
        lines.extend(
            [
                "",
                "",
                "@app.cell(hide_code=True)",
                "def _(mo, pmu):",
                f"    pmu.render_image(mo, {repr(item[1])})",
            ]
        )
        return

    if kind == "code":
        _, btn_id, step_index, cmd_str = item
        lines.extend(
            [
                "",
                "",
                "@app.cell(hide_code=True)",
                "def _(mo, pmu):",
                f"    btn_{btn_id} = pmu.make_code_button(mo, {repr(cmd_str)})",
                f"    return (btn_{btn_id},)",
                "",
                "",
                "@app.cell(hide_code=True)",
                f"def _(btn_{btn_id}, mo, pmu, session):",
                f"    pmu.render_code_step(mo, btn_{btn_id}, session, {step_index})",
            ]
        )
        return

    if kind == "check_exec":
        _, ex_id, step_index, prompt, var_name, value_str = item
        placeholder = f"{var_name} = ..."
        lines.extend(
            [
                "",
                "",
                "@app.cell(hide_code=True)",
                "def _(mo, pmu):",
                (
                    f"    ex_{ex_id} = pmu.make_exercise_form("
                    f"mo, placeholder={repr(placeholder)})"
                ),
                f"    return (ex_{ex_id},)",
                "",
                "",
                "@app.cell(hide_code=True)",
                f"def _(ex_{ex_id}, mo, pmu, session):",
                "    pmu.render_check_exec(",
                "        mo,",
                "        session,",
                f"        {step_index},",
                f"        ex_{ex_id},",
                f"        {repr(prompt)},",
                f"        {repr(var_name)},",
                f"        {repr(value_str)},",
                "    )",
            ]
        )
        return

    if kind == "check_eval":
        _, ex_id, step_index, prompt, value_str, command_str = item
        lines.extend(
            [
                "",
                "",
                "@app.cell(hide_code=True)",
                "def _(mo, pmu):",
                (
                    f"    ex_{ex_id} = pmu.make_exercise_form("
                    'mo, placeholder="Enter a Python expression...")'
                ),
                f"    return (ex_{ex_id},)",
                "",
                "",
                "@app.cell(hide_code=True)",
                f"def _(ex_{ex_id}, mo, pmu, session):",
                "    pmu.render_check_eval(",
                "        mo,",
                "        session,",
                f"        {step_index},",
                f"        ex_{ex_id},",
                f"        {repr(prompt)},",
                f"        {repr(value_str)},",
                f"        {repr(command_str)},",
                "    )",
            ]
        )


class MarimoContext(InteractionContext):
    """Execution context for compiling tutorials into Marimo (.py) notebooks."""

    def __init__(self, dictionary):
        """Initialize buffers for single-pass Marimo notebook generation."""
        super().__init__(dictionary)
        self._title = ""
        self._module_name = ""
        self._items = []
        self._steps = []
        self._pending_md = []
        self._counts = {"button": 0, "exercise": 0}

    def _flush_markdown(self):
        """Combine adjacent Markdown fragments into a single Marimo cell."""
        if self._pending_md:
            joined = "\n\n".join(self._pending_md).strip()
            if joined:
                self._items.append(("markdown", joined))
            self._pending_md = []

    def print_notebook_header(self, notebook_title, module_name=""):
        """Reset per-notebook state and record title/module metadata."""
        self._title = notebook_title
        self._module_name = module_name or notebook_title
        self._items = []
        self._steps = []
        self._pending_md = []
        self._counts = {"button": 0, "exercise": 0}

    def print_notebook_footer(self):
        """Emit the complete Marimo reactive Python notebook to stdout."""
        self._flush_markdown()
        safe_title = self._title.replace('"""', '"')
        lines = [
            "# -*- coding: utf-8 -*-",
            "#",
            "# PyClical: Python interface to GluCat:",
            "#           Generic library of universal Clifford algebra templates",
            "#",
            f"# {self._module_name}_marimo.py: {safe_title}",
            "#",
            "#    copyright            : (C) 2026 by Paul C. Leopardi",
            "#",
            (
                "# Licensed under CC BY-SA 3.0 "
                "http://creativecommons.org/licenses/by-sa/3.0/"
            ),
            '"""Generated Marimo reactive notebook."""',
            (
                "# pylint: disable=invalid-name,import-outside-toplevel,"
                "reimported,line-too-long,too-many-lines,useless-return"
            ),
            "",
            "import marimo",
            "",
            f"__generated_with = {repr(_get_marimo_version())}",
            f"app = marimo.App(app_title={repr(self._title)})",
            "",
            "",
            "@app.cell(hide_code=True)",
            "def _():",
            "    import marimo as mo",
            "",
            "    return (mo,)",
            "",
            "",
            "@app.cell(hide_code=True)",
            "def _(mo):",
            "    import pathlib as _pathlib",
            "    import sys as _sys",
            "",
            (
                "    _nb_dir = _pathlib.Path("
                "mo.notebook_dir() or _pathlib.Path.cwd()).resolve()"
            ),
            "    for _entry in (_nb_dir, _nb_dir.parent, _nb_dir.parent.parent):",
            "        _entry_str = str(_entry)",
            "        if _entry_str not in _sys.path:",
            "            _sys.path.insert(0, _entry_str)",
            "",
            "    import pyclical_marimo_utils as pmu",
            "",
            "    _steps = [",
        ]
        for cmd_str, save_snap in self._steps:
            lines.append(f"        ({repr(cmd_str)}, {repr(save_snap)}),")
        lines.extend(
            [
                "    ]",
                "    session = pmu.create_session(_steps)",
                "    return pmu, session",
            ]
        )

        for item in self._items:
            _emit_marimo_item(lines, item)

        lines.extend(
            [
                "",
                "",
                'if __name__ == "__main__":',
                "    app.run()",
            ]
        )
        print("\n".join(lines))

    def pause(self):
        """No-op for Marimo notebooks."""

    def print_head(self, output_str, indent=""):
        """Accumulate heading Markdown text."""
        self._pending_md.append("# " + output_str)

    def print_fill(self, output_str, indent="    "):
        """Accumulate body Markdown text."""
        self._pending_md.append(output_str)

    def print_line(self):
        """No-op for Marimo notebooks."""

    def print_exec(self, command_str):
        """Record an immutable executable code step."""
        if not command_str.strip():
            return
        self._flush_markdown()
        step_index = len(self._steps)
        self._steps.append((command_str, False))
        self._counts["button"] += 1
        self._items.append(
            ("code", self._counts["button"], step_index, command_str)
        )

    def display_image(self, url):
        """Record an image cell for local or remote media."""
        self._flush_markdown()
        self._items.append(("image", url))

    # pylint: disable=unused-argument
    def input_exec(self, prompt, sandbox):
        """No-op for Marimo interactive statement input."""

    def check_exec(self, prompt, var_name, value_str):
        """Record an interactive statement exercise with solution step."""
        self._flush_markdown()
        step_index = len(self._steps)
        solution_cmd = f"{var_name} = {value_str}; print({var_name})"
        self._steps.append((solution_cmd, True))
        self._counts["exercise"] += 1
        self._items.append(
            (
                "check_exec",
                self._counts["exercise"],
                step_index,
                prompt,
                var_name,
                value_str,
            )
        )

    # pylint: disable=unused-argument
    def input_eval(self, prompt):
        """No-op for Marimo interactive expression evaluation."""

    def check_eval(self, prompt, value_str, command_str):
        """Record an interactive expression exercise with solution step."""
        self._flush_markdown()
        step_index = len(self._steps)
        solution_cmd = command_str.format(value_str)
        self._steps.append((solution_cmd, True))
        self._counts["exercise"] += 1
        self._items.append(
            (
                "check_eval",
                self._counts["exercise"],
                step_index,
                prompt,
                value_str,
                command_str,
            )
        )


marimo_context = MarimoContext  # pylint: disable=invalid-name
