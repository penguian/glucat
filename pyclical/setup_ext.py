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
am_cppflags = os.environ['AM_CPPFLAGS']
ldflags  = os.environ['LDFLAGS']
#
def setup_ext(ext_name, source):
    ext = Extension(
        ext_name,         # name of extension
        sources=[source], # filename of our Cython source
        include_dirs=[".",".."],
        extra_compile_args=am_cppflags.split() + cxxflags.split(),
        extra_link_args=ldflags.split(),
    )
    return ext

# From https://github.com/SublimeCodeIntel/silvercity/blob/master/setup.py
# as at 2015-08-31
from distutils.command.build_ext import build_ext
class cxx_build_ext(build_ext):
    def build_extensions(self):
        # Allow a custom C++ compiler through the environment variables.
        compiler = os.environ.get('CXX')
        if compiler is not None:
            import sysconfig
            (ccshared, cflags) = sysconfig.get_config_vars(
                'CCSHARED', 'CFLAGS')
            args = {}
            args['compiler_so'] = compiler + ' ' + ccshared + ' ' + cflags
            self.compiler.set_executables(**args)
        build_ext.build_extensions(self)
