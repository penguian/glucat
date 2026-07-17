#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# test_pytest_doctests.py : pytest discovery wrapper for PyClical doctests
#
#    copyright            : (C) 2008-2026 by Paul C. Leopardi
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

"""
pytest discovery wrapper for PyClical doctests.
"""

import doctest
import PyClical
import pytest


def test_pyclical_doctests():
    """
Run all doctests defined within PyClical module using doctest.testmod.
"""
    results = doctest.testmod(PyClical, verbose=False)
    if results.failed > 0:
        pytest.fail(f"PyClical doctests failed: {results}")
