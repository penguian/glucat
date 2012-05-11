# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# setup_nocython.py: Use Distutils to build PyClical from C++ source..
#
#    copyright            : (C) 2008-2012 by Paul C. Leopardi
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

from distutils.core import setup
from distutils.extension import Extension
import os
ext_name = os.environ['ext_name']
source   = os.environ['source_cpp']
cxxflags = os.environ['CXXFLAGS']
includes = os.environ['INCLUDES']
ldflags  = os.environ['LDFLAGS']
ext = Extension(
    ext_name,               # name of extension
    language="c++",
    sources=[source], # filename of our C++ source
    include_dirs=[".",".."],
    extra_compile_args=includes.split()+cxxflags.split(),
    extra_link_args=ldflags.split(),
)
setup(
    name = ext_name,
    ext_modules = [ext]
)
