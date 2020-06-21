# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_utils.py: This file contains utilities for use with PyClical tutorials.
#
#    copyright            : (C) 2012-2014 by Paul C. Leopardi
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

from builtins import input
from builtins import range
from builtins import object
import sys
import numbers
import numpy as np
from PyClical import *

def get_console_width():
    import os
    default_console_width = 80
    try:
        height_str, width_str = os.popen('stty size', 'r').read().split()
        console_width = int(width_str)
        if console_width < 1:
            console_width = default_console_width
    except:
        console_width = default_console_width
    return console_width

def fill(output_str, indent = "    "):
    from textwrap import TextWrapper
    console_width = get_console_width()
    wrapper = TextWrapper(initial_indent    = indent,
                          subsequent_indent = indent,
                          width = console_width)
    return wrapper.fill(output_str)

def is_near(x, y):
    try:
        if isinstance(x, list):
            result = isinstance(y, list) and (len(x) == len(y))
            if result:
                for k in range(len(x)):
                    if not is_near(x[k], y[k]):
                        result = False
                        break
            return result
        elif isinstance(x, bool):
            return x == y
        elif isinstance(x, numbers.Real) or isinstance(x, clifford):
            tol = 4.0 * np.finfo(np.double).eps
            if abs(x) > tol:
                return abs(x-y) / abs(x) < tol
            else:
                return abs(x-y) < tol
        else:
            return x == y
    except:
         return False

def get_object_methods(obj):
    return dict([
        (method, getattr(obj, method))
        for method in dir(obj)
        if callable(getattr(obj, method)) and not method.startswith('_')])

class interaction_context(object):
    def __init__(self, dictionary):
        self.object_names = dictionary

    def pause(self):
        print("pause")

    def print_line(self):
        print("print_line")

    def print_head(self, output_str, indent = "    "):
        print("print_head: ", output_str)

    def print_fill(self, output_str, indent = "    "):
        print("print_fill: ", output_str)

    def print_exec(self, command_str):
        print("print_exec: ", command_str)

    def input_exec(self, prompt, sandbox):
        print("input_exec: ", prompt)

    def input_eval(self, prompt):
        print("input_eval: ", prompt)

    def check_exec(self, prompt, var_name, value_str):
        print("check_exec: ", prompt, var_name, value_str)

    def check_eval(self, prompt, value_str, command_str):
        print("check_exec: ", prompt, value_str, command_str)

class tutorial_context(interaction_context):
    def __init__(self, dictionary):
        self.object_names = dictionary

    def pause(self):
        if sys.stdin.isatty() and sys.stdout.isatty():
            input("Press ENTER to continue:")

    def print_line(self):
        print("")

    def print_head(self, output_str, indent = ""):
        print(fill(output_str, indent))

    def print_fill(self, output_str, indent = "    "):
        print(fill(output_str, indent))

    def print_exec(self, command_str):
        print(">>>", command_str)
        exec(command_str, self.object_names)

    def input_exec(self, prompt, sandbox):
        input_str = input(prompt + "\n>>> ")
        exec(input_str, sandbox)

    def input_eval(self, prompt):
        input_str = input(prompt + "\n>>> ")
        return eval(input_str, self.object_names)

    def check_exec(self, prompt, var_name, value_str):
        try:
            sandbox = self.object_names.copy()
            filled_prompt = fill("Exercise: Enter a Python statement to " + prompt, "")
            print("")
            self.input_exec(filled_prompt, sandbox)
            value = eval(value_str, self.object_names)
        except KeyboardInterrupt:
            raise
        except:
            value = None
        if var_name in sandbox and is_near(sandbox[var_name], value):
            print("\nThat's right.\n")
        else:
            print("\nNot quite.\n")
        self.print_fill("Here is one way to do this, and then print(the result:)")
        self.print_exec(var_name+" = "+value_str+"; print("+var_name+")")

    def check_eval(self, prompt, value_str, command_str):
        try:
            value = eval(value_str, self.object_names)
        except:
            value = None
        try:
            filled_prompt = fill("Exercise: Enter a Python expression to " + prompt, "")
            print("")
            input_value = self.input_eval(filled_prompt)
            if is_near(input_value, value):
                print("\nYour expression gives the right value.\n")
            else:
                print("\nNot quite.\n")
        except KeyboardInterrupt:
            raise
        except:
            print("\nNot quite.\n")
        self.print_fill("Here is one way to use such an expression:")
        command_str = command_str.format(value_str)
        self.print_exec(command_str)
