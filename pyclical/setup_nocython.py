# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# setup_nocython.py: Use Distutils to build PyClical from C++ source..
#
#    copyright            : (C) 2008-2022 by Paul C. Leopardi
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

from setuptools import setup
from setup_ext import setup_ext, cxx_build_ext
import os
ext_name = os.environ['ext_name']
source   = os.environ['source_cpp']
ext = setup_ext(ext_name, source)
setup(
    name = ext_name,
    cmdclass = {'build_ext': cxx_build_ext},
    ext_modules = [ext]
)
