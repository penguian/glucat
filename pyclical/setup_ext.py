# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# setup_ext.py: Use Distutils to set up an extension to use to build PyClical.
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

from distutils.extension import Extension
import os
# From stackoveflow user subdir 2012-03-16
from distutils.sysconfig import get_config_vars
(opt,) = get_config_vars('OPT')
os.environ['OPT'] = " ".join(
    flag for flag in opt.split() if flag != '-Wstrict-prototypes'
)
#
cxxflags = os.environ['CXXFLAGS']
includes = os.environ['INCLUDES']
ldflags  = os.environ['LDFLAGS']
#
def setup_ext(ext_name, source):
    ext = Extension(
        ext_name,         # name of extension
        sources=[source], # filename of our Cython source
        include_dirs=[".",".."],
        extra_compile_args=includes.split() + cxxflags.split(),
        extra_link_args=ldflags.split(),
    )
    return ext

