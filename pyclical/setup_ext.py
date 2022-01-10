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

from setuptools.extension import Extension
import os
#
cxxflags = os.environ["CXXFLAGS"]
am_cppflags = os.environ["AM_CPPFLAGS"]
libraries = os.environ["LIBRARIES"].replace("-l", "")
libraries_split = libraries.split()
#
def setup_ext(ext_name, source):
    ext = Extension(
        ext_name,         # name of extension
        sources=[source], # filename of our Cython source
        include_dirs=[".",".."],
        extra_compile_args=am_cppflags.split() + cxxflags.split(),
        libraries=libraries_split
    )
    return ext

# From https://github.com/SublimeCodeIntel/silvercity/blob/master/setup.py
# as at 2015-08-31
from setuptools.command.build_ext import build_ext
class cxx_build_ext(build_ext):

    def build_extensions(self):

        def executable_args(self, compiler):

            args = {}

            ignore_flags = (
                {
                "-Wstrict-prototypes"
                }
                if compiler.startswith("g++") else
                {
                "-fstack-protector-strong",
                "-Wdate-time",
                "-Wno-unused-result",
                "-Wstrict-prototypes"
                }
            )

            # Define compiler_so

            # From stackoveflow user subdir 2012-03-16
            # See also https://docs.python.org/2/distutils/apiref.html
            compiler_flags = [
                word for word in self.compiler.compiler_so[1:]
                if not word in ignore_flags]

            map_flag = "-ffile-prefix-map"
            ignore_map_flag = (
                (compiler.startswith("g++") and cxxmajor < 8) or
                (compiler.startswith("clang++") and cxxmajor < 10)
            )
            if ignore_map_flag:
                compiler_flags = [
                    word for word in compiler_flags
                    if not word.startswith(map_flag)]

            compiler_flags.append("-fstack-protector")
            if (compiler.startswith("icpx") and cxxmajor >= 2020):
                compiler_flags.append("-Wno-unused-command-line-argument")
            compiler_so = [compiler] + compiler_flags
            args["compiler_so"] = " ".join(compiler_so)

            # Define linker_so

            linker_flags = [
                word for word in ["-Wl,-z,notext"] + self.compiler.linker_so[1:]
                if not word in ignore_flags]
            if ignore_map_flag:
                linker_flags = [
                    word for word in linker_flags
                    if not word.startswith(map_flag)]
            linker_so =   [compiler] + linker_flags
            args["linker_so"] = " ".join(linker_so)

            return args


        def add_libraries(self):

            global libraries_split

            mkl_libraries = [
                "mkl_gf",
                "mkl_gf_lp64",
                "mkl_intel",
                "mkl_intel_lp64",
                "mkl_core"]
            using_mkl = False
            substituted = False
            for word in mkl_libraries:
                try:
                    word_index = libraries_split.index(word)
                    using_mkl = True
                    if not substituted:
                        libraries_split[word_index] = "mkl_rt"
                        substituted = True
                except ValueError:
                    pass
            if using_mkl:
                libraries_split = [
                    word for word in libraries_split
                    if not word in mkl_libraries]
            for lib in libraries_split:
                self.compiler.add_library(lib)


        # Allow a custom C++ compiler through the environment variables.
        compiler = os.environ.get("CXX")
        cxxversion = os.environ.get("CXXVERSION")
        cxxversion_split = cxxversion.split(".")
        try:
            cxxmajor = int(cxxversion_split[0])
        except ValueError:
            cxxmajor = 0
        if compiler is not None:

            args = executable_args(self, compiler)
            self.compiler.set_executables(**args)
            add_libraries(self)

        build_ext.build_extensions(self)
