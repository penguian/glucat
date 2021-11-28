# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_notebook_utils.py: This file contains utilities for use with PyClical notebooks.
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
from pyclical_tutorial_utils import *

def print_cell_markdown(source):
    print(js.dumps({
           "cell_type": "markdown",
           "metadata": {},
           "source": [
            source
           ]
          }, indent=1), end='')

def print_cell_code(code_input, execution_count):
    print(js.dumps({
           "cell_type": "code",
           "source": [
            code_input
           ],
           "metadata": {},
           "outputs": [],
           "execution_count": execution_count
          }, indent=1), end='')

def print_metadata_name(name):
    print(js.dumps({
           "name": name
          }, indent=1), end='')

class notebook_context(interaction_context):

    def __init__(self, dictionary):
        self.object_names = dictionary
        self.execution_count = 1

    def print_notebook_header(self, notebook_title):
        print('{')
        print(' "metadata": {},')
        print(' "nbformat":', 4, ',')
        print(' "nbformat_minor":', 2, ',')
        print(' "metadata": {')
        print('   "kernelspec": {')
        print('    "display_name": "Python 3",')
        print('    "language": "python",')
        print('    "name": "python3"')
        print('   },')
        print('   "language_info": {')
        print('    "codemirror_mode": {')
        print('     "name": "ipython",')
        print('     "version": 3')
        print('    },')
        print('    "file_extension": ".py",')
        print('    "mimetype": "text/x-python",')
        print('    "name": "python"')
        print('   }')
        print('  },')
        print(' "cells": [')

    def print_notebook_footer(self):
        print_cell_markdown(" ")
        print('')
        print(' ]')
        print('}')

    def pause(self):
        pass

    def print_head(self, output_str, indent = ""):
        print_cell_markdown("# " + output_str)
        print(',')

    def print_fill(self, output_str, indent = "    "):
        print_cell_markdown(output_str)
        print(',')

    def print_line(self):
        pass

    def print_exec(self, command_str):
        print_cell_code(command_str, self.execution_count)
        print(',')
        self.execution_count += 1

    def input_exec(self, prompt, sandbox):
        pass

    def check_exec(self, prompt, var_name, value_str):
        self.print_fill("Exercise: Enter a Python statement to " + prompt)
        self.print_exec("")
        self.print_fill("Here is one way to do this, and then print the result:")
        self.print_exec(var_name + " = " + value_str + "; print(" + var_name + ")")

    def input_eval(self, prompt):
        pass

    def check_eval(self, prompt, value_str, command_str):
        self.print_fill("Exercise: Enter a Python expression to " + prompt)
        self.print_exec("")
        self.print_fill("Here is one way to use such an expression:")
        command_str = command_str.format(value_str)
        self.print_exec(command_str)
