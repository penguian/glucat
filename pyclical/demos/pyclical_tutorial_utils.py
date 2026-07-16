"""
Utilities for PyClical tutorials and interactive context management.
"""
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_utils.py:
#   This file contains utilities for use with PyClical tutorials.
#
#    copyright            : (C) 2012-2026 by Paul C. Leopardi
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

import importlib
import math
import numbers
import shutil
from textwrap import TextWrapper

from builtins import input, range
import PyClical
from PyClical import *


# Stubs for tutorial context functions populated dynamically at runtime
# pylint: disable=unused-argument
def pause(*args, **kwargs):
    """
Stub for pause tutorial function.
"""


def print_line(*args, **kwargs):
    """
Stub for print_line tutorial function.
"""


def print_head(*args, **kwargs):
    """
Stub for print_head tutorial function.
"""


def print_fill(*args, **kwargs):
    """
Stub for print_fill tutorial function.
"""


def print_exec(*args, **kwargs):
    """
Stub for print_exec tutorial function.
"""


def input_exec(*args, **kwargs):
    """
Stub for input_exec tutorial function.
"""


def input_eval(*args, **kwargs):
    """
Stub for input_eval tutorial function.
"""


def check_exec(*args, **kwargs):
    """
Stub for check_exec tutorial function.
"""


def check_eval(*args, **kwargs):
    """
Stub for check_eval tutorial function.
"""


# Allowed builtins

allowed_builtins = {
    "False": False,
    "None": None,
    "True": True,
    "abs": abs,
    "bool": bool,
    "complex": complex,
    "dict": dict,
    "divmod": divmod,
    "enumerate": enumerate,
    "filter": filter,
    "float": float,
    "format": format,
    "int": int,
    "len": len,
    "list": list,
    "map": map,
    "math": math,
    "max": max,
    "min": min,
    "pow": pow,
    "print": print,
    "range": range,
    "repr": repr,
    "reversed": reversed,
    "round": round,
    "set": set,
    "sorted": sorted,
    "str": str,
    "sum": sum,
    "tuple": tuple,
    "type": type,
    "zip": zip,
}


def allowed_import(name, globals_dict=None, locals_dict=None, fromlist=(), level=0):
    """
Restricted __import__ function allowing only PyClical and IPython.
"""
    if name == "PyClical":
        return PyClical
    if name.startswith("IPython"):
        return importlib.import_module(name)
    raise ImportError(f"Import of '{name}' is not allowed in tutorial sandbox.")


allowed_builtins["__import__"] = allowed_import


def get_allowed_scope():
    """
Return a dictionary suitable for use as a restricted execution namespace.
"""
    scope = {"__builtins__": allowed_builtins}
    for name in dir(PyClical):
        if not name.startswith("_"):
            scope[name] = getattr(PyClical, name)
    return scope


def allowed_exec(command_str, scope=None):
    """
Restricted exec function that executes within an allowed scope.
"""
    if scope is None:
        scope = get_allowed_scope()
    exec(command_str, scope)  # pylint: disable=exec-used


def allowed_eval(expression_str, scope=None):
    """
Restricted eval function that evaluates within an allowed scope.
"""
    if scope is None:
        scope = get_allowed_scope()
    return eval(expression_str, scope)  # pylint: disable=eval-used


def get_console_width():
    """
    Determine terminal console width or default to 80.
    """
    default_console_width = 80
    try:
        return shutil.get_terminal_size((default_console_width, 24)).columns
    except Exception:
        return default_console_width


def fill(output_str, indent="    "):
    """
Format and wrap text according to console width.
"""
    console_width = get_console_width()
    wrapper = TextWrapper(
        initial_indent=indent, subsequent_indent=indent, width=console_width
    )
    return wrapper.fill(output_str)


def is_near(x, y):  # pylint: disable=too-many-return-statements
    """
Determine whether two elements or lists of elements are near equal.
"""
    try:
        if isinstance(x, list):
            if not isinstance(y, list) or len(x) != len(y):
                return False
            for elem_x, elem_y in zip(x, y):
                if not is_near(elem_x, elem_y):
                    return False
            return True
        if isinstance(x, bool):
            return x == y
        if isinstance(x, (numbers.Real, PyClical.clifford)):
            tol = 4.0 * PyClical.scalar_epsilon
            if abs(x) > tol:
                return abs(x - y) / abs(x) < tol
            return abs(x - y) < tol
        return x == y
    except Exception:
        return False


def get_object_methods(obj):
    """
Retrieve all callable public methods of an object as a dictionary.
"""
    return {
        method: getattr(obj, method)
        for method in dir(obj)
        if callable(getattr(obj, method)) and not method.startswith("_")
    }


class InteractionContext:
    """
Base class for managing tutorial interaction execution contexts.
"""

    def __init__(self, dictionary):
        self.object_names = dictionary

    def pause(self):
        """
Pause interaction.
"""
        print("pause")

    def print_line(self):
        """
Print a line separation.
"""
        print("print_line")

    def print_head(self, output_str, indent="    "):
        """
Print a header string.
"""
        print("print_head: ", output_str)

    def print_fill(self, output_str, indent="    "):
        """
Print wrapped text string.
"""
        print("print_fill: ", output_str)

    def print_exec(self, command_str):
        """
Execute and print a command.
"""
        print("print_exec: ", command_str)

    def input_exec(self, prompt, sandbox):
        """
Execute interactive code input.
"""
        print("input_exec: ", prompt)

    def input_eval(self, prompt):
        """
Evaluate interactive code input.
"""
        print("input_eval: ", prompt)

    def check_exec(self, prompt, var_name, value_str):
        """
Validate execution result against target variable.
"""
        print("check_exec: ", prompt, var_name, value_str)

    def check_eval(self, prompt, value_str, command_str):
        """
Validate evaluation result against expected value.
"""
        print("check_exec: ", prompt, value_str, command_str)


class TutorialContext(InteractionContext):
    """
Interactive execution context used for running PyClical tutorials.
"""

    def __init__(self, dictionary=None):
        super().__init__(get_allowed_scope())

    def pause(self):
        """
Prompt user to press ENTER if running in terminal.
"""
        if sys.stdin.isatty() and sys.stdout.isatty():
            input("Press ENTER to continue:")

    def print_line(self):
        """
Print an empty line.
"""
        print("")

    def print_head(self, output_str, indent=""):
        """
Print filled header string.
"""
        print(fill(output_str, indent))

    def print_fill(self, output_str, indent="    "):
        """
Print filled body string.
"""
        print(fill(output_str, indent))

    def print_exec(self, command_str):
        """
Display command and execute in allowed scope.
"""
        print(">>>", command_str)
        allowed_exec(command_str, self.object_names)

    def input_exec(self, prompt, sandbox):
        """
Prompt user for input and execute in sandbox scope.
"""
        input_str = input(prompt + "\n>>> ")
        allowed_exec(input_str, sandbox)

    def input_eval(self, prompt):
        """
Prompt user for expression input and evaluate in allowed scope.
"""
        input_str = input(prompt + "\n>>> ")
        return allowed_eval(input_str, self.object_names)

    def check_exec(self, prompt, var_name, value_str):
        """
Interactively test user command against reference value.
"""
        try:
            sandbox = self.object_names.copy()
            filled_prompt = fill(
                "Exercise: Enter a Python statement to " + prompt, ""
            )
            print("")
            self.input_exec(filled_prompt, sandbox)
            value = allowed_eval(value_str, self.object_names)
        except Exception:
            value = None
        if var_name in sandbox and is_near(sandbox[var_name], value):
            print("\nThat's right.\n")
        else:
            print("\nNot quite.\n")
        self.print_fill(
            "Here is one way to do this, and then print the result:"
        )
        self.print_exec(
            var_name + " = " + value_str + "; print(" + var_name + ")"
        )

    def check_eval(self, prompt, value_str, command_str):
        """
Interactively test user expression against reference value.
"""
        try:
            value = allowed_eval(value_str, self.object_names)
        except Exception:
            value = None
        try:
            filled_prompt = fill(
                "Exercise: Enter a Python expression to " + prompt, ""
            )
            print("")
            input_value = self.input_eval(filled_prompt)
            if is_near(input_value, value):
                print("\nYour expression gives the right value.\n")
            else:
                print("\nNot quite.\n")
        except Exception:
            print("\nNot quite.\n")
        self.print_fill("Here is one way to use such an expression:")
        command_str = command_str.format(value_str)
        self.print_exec(command_str)


# Aliases for backwards compatibility with legacy tutorials
interaction_context = InteractionContext  # pylint: disable=invalid-name
tutorial_context = TutorialContext  # pylint: disable=invalid-name
