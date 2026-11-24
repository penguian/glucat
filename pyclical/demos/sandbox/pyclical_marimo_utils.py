"""
Utilities for running and rendering PyClical interactive Marimo notebooks.
"""
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_marimo_utils.py:
#   Runtime helpers for PyClical Marimo reactive notebooks.
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

import contextlib
import io
from pathlib import Path
import re


def _max_backtick_run(text):
    """Return the length of the longest contiguous run of backticks in text."""
    return max((len(match) for match in re.findall(r"`+", text)), default=0)


def _code_fence(text, lang="text"):
    """Wrap text in a fenced Markdown code block safe against embedded backticks."""
    fence = "`" * max(3, _max_backtick_run(text) + 1)
    return f"{fence}{lang}\n{text}\n{fence}"


def _clean_prompt(prompt):
    """Normalize internal whitespace in an exercise prompt string."""
    return " ".join(prompt.split())


def _clone_scope(scope, pyclical_types):
    """Return an isolated shallow-plus-PyClical copy of a tutorial scope dict."""
    cloned = {}
    for key, val in scope.items():
        if key == "__builtins__":
            cloned[key] = val
        elif pyclical_types and isinstance(val, pyclical_types):
            cloned[key] = type(val)(val)
        elif isinstance(val, (list, dict, set)):
            cloned[key] = val.copy()
        else:
            cloned[key] = val
    return cloned


class MarimoSession:
    """Container holding pre-executed outputs, errors, and exercise snapshots."""

    def __init__(self, steps, pyclical_types=(), init_error=None):
        """Initialize an empty session for a sequence of tutorial steps."""
        self.steps = list(steps)
        self.pyclical_types = pyclical_types
        self.init_error = init_error
        self.outputs = [""] * len(self.steps)
        self.errors = [None] * len(self.steps)
        self.snapshots = {}

    def get_output(self, step_index):
        """Return captured stdout text for the given step index."""
        return self.outputs[step_index]

    def get_error(self, step_index):
        """Return captured error message for the given step index, or None."""
        return self.errors[step_index]

    def get_snapshot(self, step_index):
        """Return a fresh isolated scope snapshot captured before step_index."""
        return _clone_scope(self.snapshots[step_index], self.pyclical_types)


def create_session(steps):
    """Execute tutorial steps sequentially and capture outputs and snapshots."""
    step_list = list(steps)
    try:
        # Lazy imports keep Markdown cells decoupled from compiled PyClical.
        import PyClical  # pylint: disable=import-outside-toplevel
        import pyclical_tutorial_utils as ptu  # pylint: disable=import-outside-toplevel
    except Exception as exc:
        return MarimoSession(
            step_list,
            init_error=f"{type(exc).__name__}: {exc}",
        )

    pyclical_types = (PyClical.clifford, PyClical.index_set)
    session = MarimoSession(step_list, pyclical_types=pyclical_types)
    scope = ptu.get_allowed_scope()

    for idx, (code_str, save_snapshot) in enumerate(step_list):
        if save_snapshot:
            session.snapshots[idx] = _clone_scope(scope, pyclical_types)
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                ptu.allowed_exec(code_str, scope)
        except Exception as exc:
            session.errors[idx] = f"{type(exc).__name__}: {exc}"
        session.outputs[idx] = buf.getvalue()

    return session


def make_code_button(mo, code_str):
    """Create a compact toggleable play button for a Python code step."""
    return mo.ui.button(
        value=False,
        on_click=lambda shown: not shown,
        label="▶",
        kind="neutral",
        tooltip=f"Click to show or hide output: {code_str}",
    )


def render_code_step(mo, btn, session, step_index):
    """Render a play button beside an immutable code block and toggled output."""
    code_str = session.steps[step_index][0]
    play_btn = mo.style(btn, {"flex-shrink": "0"})
    code_block = mo.style(
        mo.md(_code_fence(code_str, lang="python")),
        {"flex": "1", "min-width": "0"},
    )
    row = mo.hstack(
        [play_btn, code_block],
        justify="start",
        align="center",
        gap=0.5,
    )
    if not btn.value:
        return row
    if session.init_error:
        msg = f"**PyClical unavailable:** `{session.init_error}`"
        return mo.vstack([row, mo.callout(mo.md(msg), kind="danger")])
    err = session.get_error(step_index)
    if err:
        return mo.vstack([row, mo.callout(mo.md(_code_fence(err)), kind="warn")])
    out = session.get_output(step_index).rstrip("\n")
    if out:
        return mo.vstack([row, mo.md(_code_fence(out))])
    return mo.vstack([row, mo.md("_Executed (no output)._")])


def make_exercise_form(mo, placeholder):
    """Create a code editor form with a Check button for tutorial exercises."""
    return mo.ui.code_editor(
        value="",
        language="python",
        placeholder=placeholder,
        min_height=60,
        show_copy_button=False,
    ).form(
        submit_button_label="Check",
        bordered=True,
        clear_on_submit=False,
    )


def _build_solution_accordion(mo, intro_text, solution_code, session, step_index):
    """Build a collapsible accordion showing the reference solution and output."""
    parts = [intro_text, _code_fence(solution_code, lang="python")]
    if session.init_error:
        parts.append(f"**PyClical unavailable:** `{session.init_error}`")
    else:
        err = session.get_error(step_index)
        if err:
            parts.append("**Output error:**\n\n" + _code_fence(err))
        else:
            out = session.get_output(step_index).rstrip("\n")
            parts.append(_code_fence(out) if out else "_Executed (no output)._")
    return mo.accordion({"Show solution": mo.md("\n\n".join(parts))})


def _eval_reference(ptu, expr_str, scope):
    """Evaluate a reference expression, returning (value, error_str)."""
    try:
        return ptu.allowed_eval(expr_str, scope), None
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"


# pylint: disable=too-many-arguments,too-many-positional-arguments
def _format_exec_result(mo, ptu, sandbox, expected, var_name, captured):
    """Format grading feedback after successful student statement execution."""
    if var_name not in sandbox:
        msg = f"**Not quite.** Your statement did not assign `{var_name}`."
        if captured:
            msg += "\n\n" + _code_fence(captured)
        return mo.callout(mo.md(msg), kind="warn")

    is_correct = ptu.is_near(sandbox[var_name], expected)
    header = (
        "**That's right.**" if is_correct else "**Not quite.** Your statement gave:"
    )
    parts = [header, _code_fence(f"{var_name} = {sandbox[var_name]}")]
    if captured:
        parts.append(_code_fence(captured))
    return mo.callout(
        mo.md("\n\n".join(parts)),
        kind="success" if is_correct else "warn",
    )


# pylint: disable=too-many-arguments,too-many-positional-arguments
def _grade_exec(mo, session, step_index, user_code, var_name, value_str):
    """Execute and grade a student statement for a check_exec exercise."""
    if not user_code:
        return mo.callout(
            mo.md("Enter a Python statement and click **Check**."),
            kind="info",
        )
    if session.init_error:
        return mo.callout(
            mo.md(f"**PyClical unavailable:** `{session.init_error}`"),
            kind="danger",
        )

    import pyclical_tutorial_utils as ptu  # pylint: disable=import-outside-toplevel

    expected, expected_err = _eval_reference(
        ptu, value_str, session.get_snapshot(step_index)
    )
    sandbox = session.get_snapshot(step_index)
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            ptu.allowed_exec(user_code, sandbox)
    except Exception as exc:
        return mo.callout(
            mo.md(
                "**Not quite.** Error running your statement:\n\n"
                + _code_fence(f"{type(exc).__name__}: {exc}")
            ),
            kind="warn",
        )

    captured = buf.getvalue().rstrip("\n")
    if expected_err:
        parts = [f"Reference expression raised `{expected_err}`."]
        if var_name in sandbox:
            parts.append(_code_fence(f"{var_name} = {sandbox[var_name]}"))
        if captured:
            parts.append(_code_fence(captured))
        return mo.callout(mo.md("\n\n".join(parts)), kind="warn")

    return _format_exec_result(mo, ptu, sandbox, expected, var_name, captured)


def _grade_eval(mo, session, step_index, user_expr, value_str):
    """Evaluate and grade a student expression for a check_eval exercise."""
    if not user_expr:
        return mo.callout(
            mo.md("Enter a Python expression and click **Check**."),
            kind="info",
        )
    if session.init_error:
        return mo.callout(
            mo.md(f"**PyClical unavailable:** `{session.init_error}`"),
            kind="danger",
        )

    import pyclical_tutorial_utils as ptu  # pylint: disable=import-outside-toplevel

    expected, expected_err = _eval_reference(
        ptu, value_str, session.get_snapshot(step_index)
    )
    sandbox = session.get_snapshot(step_index)
    try:
        user_val = ptu.allowed_eval(user_expr, sandbox)
    except Exception as exc:
        err_msg = f"{type(exc).__name__}: {exc}"
        return mo.callout(
            mo.md(
                "**Not quite.** Error evaluating your expression:\n\n"
                + _code_fence(err_msg)
            ),
            kind="warn",
        )

    if expected_err:
        parts = [
            f"Reference expression raised `{expected_err}`.",
            _code_fence(str(user_val)),
        ]
        return mo.callout(mo.md("\n\n".join(parts)), kind="warn")

    is_correct = ptu.is_near(user_val, expected)
    header = (
        "**Your expression gives the right value.**"
        if is_correct
        else "**Not quite.** Your expression evaluated to:"
    )
    return mo.callout(
        mo.md(header + "\n\n" + _code_fence(str(user_val))),
        kind="success" if is_correct else "warn",
    )


# pylint: disable=too-many-arguments,too-many-positional-arguments
def render_check_exec(
    mo, session, step_index, exercise_form, prompt, var_name, value_str
):
    """Render an interactive statement exercise with grading and solution."""
    prompt_text = _clean_prompt(prompt)
    elements = [
        mo.md(f"**Exercise:** Enter a Python statement to {prompt_text}"),
        exercise_form,
    ]
    if exercise_form.value is not None:
        elements.append(
            _grade_exec(
                mo,
                session,
                step_index,
                exercise_form.value.strip(),
                var_name,
                value_str,
            )
        )
    solution_code = f"{var_name} = {value_str}; print({var_name})"
    elements.append(
        _build_solution_accordion(
            mo,
            "Here is one way to do this, and then print the result:",
            solution_code,
            session,
            step_index,
        )
    )
    return mo.vstack(elements)


# pylint: disable=too-many-arguments,too-many-positional-arguments
def render_check_eval(
    mo, session, step_index, exercise_form, prompt, value_str, command_str
):
    """Render an interactive expression exercise with grading and solution."""
    prompt_text = _clean_prompt(prompt)
    elements = [
        mo.md(f"**Exercise:** Enter a Python expression to {prompt_text}"),
        exercise_form,
    ]
    if exercise_form.value is not None:
        elements.append(
            _grade_eval(
                mo,
                session,
                step_index,
                exercise_form.value.strip(),
                value_str,
            )
        )
    solution_code = command_str.format(value_str)
    elements.append(
        _build_solution_accordion(
            mo,
            "Here is one way to use such an expression:",
            solution_code,
            session,
            step_index,
        )
    )
    return mo.vstack(elements)


def render_image(mo, url):
    """Render an online URL or sandbox-relative image using marimo.image."""
    if url.startswith(("http://", "https://")):
        src = url
    else:
        nb_dir = Path(mo.notebook_dir() or Path.cwd()).resolve()
        src = str((nb_dir / url).resolve())
    return mo.image(src=src, alt=Path(url).name)
