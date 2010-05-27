# -*- coding: utf-8 -*-
import Cython.Compiler.Options
Cython.Compiler.Options.annotate=True
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os
cxxflags=os.environ['CXXFLAGS']
includes=os.environ['INCLUDES']
ldflags=os.environ['LDFLAGS']
ext = Extension(
"PyCliCal", # name of extension
    ["PyCliCal.pyx"], # filename of our Cython source
    language="c++", # this causes Cython to create C++ source
    include_dirs=[".",".."],
    extra_compile_args=includes.split()+cxxflags.split(),
    extra_link_args=ldflags.split(),
)
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [ext]
)
