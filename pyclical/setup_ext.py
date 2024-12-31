# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# setup_ext.py: Use Distutils to set up an extension to use to build PyClical.
#
#    copyright            : (C) 2008-2024 by Paul C. Leopardi
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


def filtered_libraries():
    libraries_list = os.environ["LIBRARIES"].replace("-l", "").split()
    filtered_libraries_list = []
    mkl_libraries = [
        "mkl_gf",
        "mkl_gf_lp64",
        "mkl_intel",
        "mkl_intel_lp64"]
    substituted = False
    for lib in mkl_libraries:
        try:
            lib_index = libraries_list.index(lib)
            if not substituted:
                libraries_list[lib_index] = "mkl_rt"
                substituted = True
        except ValueError:
            pass
    for lib in libraries_list:
        if lib not in mkl_libraries:
            filtered_libraries_list.append(lib)
    return filtered_libraries_list


def setup_ext(ext_name, source):
    all_includes_list = os.environ["all_includes"].replace("-I", "").split()
    ext = Extension(
        ext_name,
        sources=[source],
        include_dirs=all_includes_list,
        libraries=filtered_libraries())
    return ext

