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
        new_compiler = os.environ.get('CXX')
        cxxversion = os.environ.get('CXXVERSION')
        cxxversion_split = cxxversion.split(".")
        try:
            cxxmajor = int(cxxversion_split[0])
        except ValueError:
            cxxmajor = 0
        if new_compiler is not None:
            # From stackoveflow user subdir 2012-03-16
            # See also https://docs.python.org/2/distutils/apiref.html
            ignore_flags = (
                {
                '-Wstrict-prototypes'
                }
                if new_compiler.startswith("g++") else
                {
                '-fstack-protector-strong',
                '-Wdate-time',
                '-Wno-unused-result',
                '-Wstrict-prototypes'
                }
            )
            map_flag = "-ffile-prefix-map"
            ignore_map_flag = (
                (new_compiler.startswith("g++") and cxxmajor < 8) or
                (new_compiler.startswith("clang++") and cxxmajor < 10)
            )
            new_compiler_flags = [
                word for word in self.compiler.compiler_so[1:]
                if not word in ignore_flags]
            if ignore_map_flag:
                new_compiler_flags = [
                    word for word in new_compiler_flags
                    if not word.startswith(map_flag)]
            new_compiler_flags.append('-fstack-protector')

            new_linker_flags = [
                word for word in self.compiler.linker_so[1:]
                if not word in ignore_flags]
            if ignore_map_flag:
                new_linker_flags = [
                    word for word in new_linker_flags
                    if not word.startswith(map_flag)]

            new_compiler_so = [new_compiler] + new_compiler_flags
            new_linker_so =   [new_compiler] + new_linker_flags
            args = {}
            args['compiler_so'] = " ".join(new_compiler_so)
            args['linker_so'] = " ".join(new_linker_so)
            self.compiler.set_executables(**args)
        build_ext.build_extensions(self)
