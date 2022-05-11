INSTALL for GluCat 0.12.0 with PyClical
========================================

Prerequisites: Before You Begin
===============================

GluCat uses the C++ Standard Library and the Boost Library. The PyClical Python
extension module is built using Cython and Python. Make sure that you have all
of these installed and working before attempting to use GluCat with PyClical.
Use the instructions at http://www.boost.org/more/download.html to obtain
the Boost Library. Make sure that you are able to build the Boost library with
the same C++ compiler version as you will be using to build programs that will
use the GluCat library.

The default configuration of GluCat, as well as most combinations of
configuration options, require a compiler that supports the C++ 2011 standard.

Scroll to the end of these instructions to see a list of successful builds,
including version numbers of various software components, and notes on software
versions.


Installation from Git clone
===========================

You can install GluCat in one of two ways:
1. By cloning one of the Git repositories for GluCat;
2. By downloading a tarball.

To install the first way, from (e.g.) GitHub, run the following commands on a
Linux machine or equivalent Posix environment connected to the Internet:
```
> git clone git@github.com:penguian/glucat.git glucat-0.12.0
> cd glucat-0.12.0
> make -f admin/Makefile.common cvs
```
This results in a directory structure that includes glucat-0.12.0/configure,
allowing you to make and install GluCat in the same way as if you had downloaded
and unzipped the tarball glucat-0.12.0.tar.gz.


Directory Structure
===================

Once you have downloaded, unzipped and untarred the source code, or followed
the instructions above to install from Git clone, you should have a directory,
glucat-0.12.0. Under this directory you should see a number of subdirectories,
including `./admin`, `./doc`, `./glucat`, `./gfft_test`, `./products`,
`./pyclical`, `./squaring`, `./test`, `./test_runtime`, `./testxx`, and
`./transforms`.

The following instructions are meant to be used with `glucat-0.12.0` as the
current directory.


Installation Summary
====================

There are four different types of program that you can make using the GluCat
library with PyClical. These are:

 1. The PyClical Python extension module.
 2. The timing test programs.
 3. The regression test programs.
 4. Your own programs, written in C++.

There are four different things to install. These are:

 1. The PyClical Python extension module.
 2. The GluCat C++ header files.
 3. The basic documentation, including the PyClical demonstration scripts and
    notebooks.
 4. The GluCat API documentation.

The simplest way to install GluCat with PyClical is to run the following
commands, in order:
```
 ./configure
 make
 sudo make install
 make clean
```
This makes and installs the PyClical Python extension module and also
installs the GluCat header files, as well as the basic doumentation.

You can also make and install the GluCat API documentation, after running
`./configure`, by running the following commands, in order:
```
 make doc
 sudo make install-doc
 make clean
```
More details of the installation process are given below.


The Four Different Types of Programs
=====================================

The PyClical Python extension module
------------------------------------

GluCat includes a Python extension module called PyClical.

The subdirectory `./pyclical` contains source code for the PyClical extension
module, and the subdirectory `./pyclical/demos` contains demo and tutorial
scripts written in Python.

PyClical is written in C++ and Cython, and is defined in the files
`pyclical/glucat.pxd`, `pyclical/PyClical.h`, `pyclical/PyClical.pxd`, and
`pyclical/PyClical.pyx`.

The two types of test programs
------------------------------

GluCat includes two types of test programs, timing (benchmark) tests, and
functionality (regression) tests. The regression tests are also used as
programming examples. The source code for the timing tests is in `./gfft_test`,
`./products`, `./squaring` and `./transforms`. The source code for the regression
tests is in the subdirectories `./test00` to `./test17`. Compilation of these test
programs uses C++ headers defined in `./test` as well as `./glucat`.

Your own programs, written in C++
---------------------------------

Once you have installed  the GluCat C++ library, you use this library to build
your own C++ programs. To do this, you will need to follow some or all of the
same steps you took to build the test programs and the PyClical extension
module. This includes setting compiler flags, and including relevant headers
and libraries.

Pay special attention to the `-D` flags described in the configuration section
below, as these control optional parts of the compilation of the GluCat library.
In particular, if you compile your own programs without setting any of these
`-D` flags, please be aware that the default random number generator used by
GluCat requires the use of the C++ 2011 standard. See details below.

To automate the build process, you should use GNU make with your own Makefile.
For very elaborate software, you may also want to consider using GNU Autotools.


To Configure
============

The build process for PyClical and the test programs uses a single `./configure`
shell script and a single overall Makefile, plus Makefiles in the relevant
subdirectories.

As briefly described above, the simplest way to install this package is:

 1. `cd` to the `glucat-0.12.0` directory containing the source code and type
    `./configure` to configure GluCat with PyClical for your system.
    If you are using `csh` on an old version of System V, you might need to type
    `sh ./configure` instead to prevent `csh` from trying to execute
    `./configure` itself.

    Running `./configure` takes a while.  While running, it prints some messages
    telling which features it is checking for.

 2. Type `make` to compile PyClical and the timing tests.

 3. Type `make install` to install PyClical, the PyClical demos, the `glucat`
    header files, and any remaining data files and documentation. As with
    `make install-doc`, you may need to use `sudo` or you may need to log
    in as `root` to perform this operation.

 4. You can remove the program binaries and object files from the source code
    directory by typing `make clean`.

This section describes the first step in detail.

The `./configure` shell script attempts to guess correct values for various
system-dependent variables used during compilation.  It uses those values to
create a `Makefile` in each relevant subdirectory.  It may also create one or
more `.h` files containing system-dependent definitions.  Finally, it creates a
shell script `./config.status` that you can run in the future to recreate the
current configuration, a file `./config.cache` that saves the results of its
tests to speed up reconfiguring, and a file `./config.log` containing compiler
output (useful mainly for debugging `./configure`).

If you need to do unusual things to compile GluCat with PyClical, please try to
figure out how `./configure` could check whether to do them, and email diffs or
instructions to the address given in the `./README` file so they can be considered
for the next release.  If at some point `./config.cache` contains results you don't
want to keep, you may remove or edit it.

The file `./configure.ac` is used to create `./configure` via a program called
`autoconf`.  You only need `./configure.ac` if you want to change or regenerate
`./configure` using a newer version of `autoconf`.


Compilers and Options
---------------------

Some systems require unusual options for compilation or linking that the
`./configure` script does not know about.  You can give `./configure` initial
values for variables by setting them in the environment.  Using a
Bourne-compatible shell, you can do that on the command line like this:
```
    CC=c89 CFLAGS=-O2 LIBS=-lposix ./configure
```
Or on systems that have the `env` program, you can do it like this:
```
    env CPPFLAGS=-I/usr/local/include LDFLAGS=-s ./configure
```

Optional Features
-----------------

Defaults for the options are specified in brackets.

Installation directories:

```
  --prefix=PREFIX         install architecture-independent files in PREFIX
                          [/usr/local]
```
By default, `make install` will install the GluCat C++ headers in
`/usr/local/include/glucat`, and will install the `PyClical.so` shared object file
in a directory under `/usr/local` whose full pathname is determined by Python, for
example `/usr/local/python3.7/site-packages`. You can specify an installation
prefix other than /usr/local by using the option `--prefix=PREFIX`, for example
`--prefix=$HOME`.

```
  --exec-prefix=EPREFIX   install architecture-dependent files in EPREFIX
                          [PREFIX]
```
The option `--exec-prefix=EPREFIX` changes the installation prefix for PyClical,
without affecting the installation prefix for the GluCat C++ headers.

```
  --enable-debug[=ARG]    enables debug symbols (yes|no|full) [default=no]
```
This option controls debugging and optimization options. ARG can be `yes`, `no`
or `full`. The default is `no`.

The option `--enable-debug=no` adds the compiler flag `-DNDEBUG` to `CXXFLAGS` in
the Makefiles, and turns on optimizations, such as adding the compiler flag
`-O3` to `CXXFLAGS` in the Makefiles, as well as a number of compiler and version
dependent optimization flags. The optimized compilation may take a long time to
complete.

The option `--enable-debug=yes` adds the compiler flag `-DNDEBUG` to `CXXFLAGS` in
the Makefiles, and turns on debugging, by adding the compiler flags `-O1 -g` to
`CXXFLAGS` in the Makefiles.

The option `--enable-debug=full` turns on full debugging, by adding the compiler
flags `-O0 -g3` to `CXXFLAGS` in the Makefiles, and does not add the compiler flag
`-DNDEBUG` to `CXXFLAGS` in the Makefiles.

The preprocessor symbol `NDEBUG` is used by Boost uBLAS to control debugging.
If NDEBUG is defined, then uBLAS compiles in release mode, including the use of
expression templates. If `NDEBUG` is not defined, then uBLAS compiles in debug
mode.


If you are compiling your own programs using the GluCat library, to control
debugging and performance, your Makefile needs to pass the following flags to
the C++ compiler.

 For full production and timing tests:
  `-DNDEBUG -O3`

 For testing with some debugging capability:
  `-DNDEBUG -O1 -g`

 For full debugging capability, including the use of debug code in uBLAS:
  `-O0 -g3`

```
  --disable-debug         disables debug output and debug symbols [default=yes]
```
Since the default for `--enable-debug` is `no`, the option `--disable-debug`
does nothing.

```
  --enable-pyclical       uses Cython to build PyClical Python extension module
                          [default=yes]
```
This option determines whether the PyClical Python extension module is built.
The option `--enable-pyclical` sets the conditional flag `make_pyclical` in the
Makefile for GluCat, and adds the directory `./pyclical` to the `SUBDIRS` list,
which is the list of subdirectories to be built. This is the default behaviour.
If you do not want to build PyClical, use the option `--disable-pyclical`.

For PyClical to build successfully using `pyclical/Makefile`, you will need a
recent version of Python and (optionally) a recent version of Cython. You will
also need to ensure that the compiler can find the header file <Python.h>, by
using the `--with-extra-includes` option, if necessary. If your computer has
Python installed, but not Cython, you can build PyClical from the provided
`PyClical.cpp` file. See details in `To Build` below.


Optional Packages
-----------------

```
  --with-extra-includes=DIR
```
This option adds non standard include paths.

```
  --with-extra-libs=DIR
```
This option adds non standard library paths.

```
  --with-demo-dir=DIR     [default=$DATAROOTDIR/pyclical/demos]
```
This option defines the installation directory for the PyClical demos.

```
  --with-qd               uses dd_real and qd_real [default=no]
```
This option controls the use of the QD high precision floating point library.

The option `--with-qd` adds `-D_GLUCAT_USE_QD` to `CXXFLAGS` and adds the flag
`-lqd` to the list of libraries, `LIBS` in the Makefiles, if the header
file `<qd/qd_real.h>` and the library `libqd` are usable.


If `_GLUCAT_USE_QD` is defined, `glucat/qd.h` includes `<qd/qd_real.h>` and supports
the use of QD in GluCat by defining specializations for `numeric_traits<dd_real>`
and `numeric_traits<qd_real>`.


To compile your own programs using the GluCat library with QD, your Makefile
needs to pass the flags `-D_GLUCAT_USE_QD` and `-lqd` to the C++ compiler.
You will also need to ensure that the include path used by the compiler sees
`<qd/qd_real.h>` and the library path sees `libqd.*`.

```
  --with-eig[=ARG]        library to use for eigenvalues
                          (no|bindings|blaze) [default=no]
```
This option is used to control `_GLUCAT_USE_EIGENVALUES` and determine which
libraries to use. ARG can be `no` or `bindings`. The default is `no`.

The option `--with-eig=bindings` adds
`-D_GLUCAT_USE_EIGENVALUES -D_GLUCAT_USE_BINDINGS` to `CXXFLAGS` and adds the
flags `-llapack -lblas` to the list of libraries, `LIBS` in the
Makefiles, if the header file `<boost/numeric/bindings/driver/lapack/gees.hpp>`
and the libraries `liblapack` and `libblas` are usable. To accomplish this, the
configure script uses the `AX_LAPACK` and `AX_BLAS macros`, as documented at
at https://www.gnu.org/software/autoconf-archive/The-Macros.html
This, in turn means that you will need to have a Fortran compiler installed,
and preferably have the F77 environment variable set to refer to this compiler.

The option `--with-eig=blaze` adds `-D_GLUCAT_USE_EIGENVALUES -D_GLUCAT_USE_BLAZE`
to `CXXFLAGS` and adds the flags `-llapack -lblas` to the list of libraries,
`LIBS` in the Makefiles, if the header file `<blaze/Math.h>` and the libraries
`liblapack` and `libblas` are usable. To accomplish this, the configure script
uses the `AX_LAPACK` and `AX_BLAS macros`, as mentioned above.

The preprocessor symbol `_GLUCAT_USE_EIGENVALUES` controls whether the `sqrt()`
and `log()` functions in `glucat/matrix_multi_imp.h`  detect and handle negative
real eigenvalues and imaginary eigenvalues correctly. If `_GLUCAT_USE_EIGENVALUES`
is defined, then `sqrt()` and `log()` call the function `classify_eigenvalues()`,
defined in `glucat/matrix_imp.h`, to detect negative real eigenvalues and imaginary
eigenvalues, and handle negative real eigenvalues by expanding the algebra.
Otherwise, `sqrt()` and `log()` operate as per GluCat 0.5.0 and earlier, which gives
incorrect results in the case of negative real eigenvalues.

The function `eigenvalues()` in `glucat/matrix_imp.h` calls an external function
to obtain the eigenvalues of a matrix. Which function is used depends on one of
a number of preprocessor symbols:

If `_GLUCAT_USE_BINDINGS` is defined, `glucat/matrix_imp.h` includes
`<boost/numeric/bindings/lapack/driver/gees.hpp>` and uses the Boost Numeric
Bindings library. To use this library, you will need to download and install
it yourself, preferably from https://github.com/uBLAS/numeric_bindings

If `_GLUCAT_USE_BLAZE` is defined, `glucat/matrix_imp.h` includes
`<blaze/Math.h>` and related Blaze include files, as per the Blaze template
library. To use this library, you will need to install it yourself, preferably
from https://bitbucket.org/blaze-lib/blaze/src/master/

To compile your own programs using the GluCat library, to detect and correctly
handle negative real eigenvalues in the `sqrt()` and `log()` functions, your
Makefile needs to pass the flag `-D_GLUCAT_USE_EIGENVALUES` to the C++ compiler,
as well as one of the following choices of flags, and the corresponding header
files and libraries must be usable.

* For Boost Numeric Bindings:
  `-D_GLUCAT_USE_BINDINGS -llapack -lblas`
  You will also need to ensure that the include path used by the compiler sees
  `<boost/numeric/bindings/lapack/driver/gees.hpp>` and the library path sees
  `liblapack.*` and `libblas.*`.
* For Blaze:
  `-D_GLUCAT_USE_BLAZE -llapack -lblas`
  You will also need to ensure that the include path used by the compiler sees
  `<blaze/Math.h>` etc. and the library path sees `liblapack.*` and `libblas.*`.
  Blaze also requires C++14, so your Makefile needs to use `-std=c++14` or the
  equivalent for your C++ compiler.


Operation Controls
------------------

 `./configure` recognizes the following options to control how it operates.
```
 --cache-file=FILE
      Use and save the results of the tests in FILE instead of
      ./config.cache.  Set FILE to `/dev/null` to disable caching, to debug the
      ./configure script.

 --help
      Print a summary of the options to ./configure, and exit.

 --quiet
 --silent
 -q
      Do not print messages saying which checks are being made.

 --srcdir=DIR
      Look for source code in directory DIR.  Usually ./configure can determine
      that directory automatically.

 --version
      Print the version of Autoconf used to generate the ./configure script,
      and exit.
```
 `./configure` also accepts some other, not widely useful, options.


To Build
========

Building PyClical and the timing test programs
----------------------------------------------

To build PyClical and the timing test programs, set the environment variable
CXX to indicate your C++ compiler, eg. `g++` for GNU C++, `icpc` or `icpx`
for Intel C++, then run `./configure` as above, and then run `make`.

Make uses the headers in `./glucat` and `./test` and the source in `./pyclical`,
`./gfft_test`, `./products`, `./squaring` and `./transforms` to build PyClical
and the timing test programs `./gfft_test/gfft_test`, `./products/products`,
`./squaring/squaring` and `./transforms/transforms`.

The following build steps will be performed if you have not selected the
configuration option `--disable-pyclical` (equivalently, `--enable-pyclical=no`).

If Cython is installed then `make` builds PyClical by running the command
```
  ext_name=`PyClical` source_pyx=`PyClical.pyx` \
  CXX=`$(CXX)` CXXVERSION=`$(CXXVERSION)` CFLAGS=`` \
  CXXFLAGS=`$(EXTCXXFLAGS)` AM_CPPFLAGS=`$(EXTAM_CPPFLAGS)` \
  LDFLAGS=`$(USER_LDFLAGS)` LIBRARIES=`$(LIBS)` \
 $(PYTHON) setup.py build_ext --inplace
```
with `EXTCXXFLAGS = $(glucat_extra_cxxflags_pyclical) $(CXXFLAGS)`,
`EXTAM_CPPFLAGS=$(all_includes)`, and the values of the other environment variables
set by `./configure`.

You can run `pyclical/setup.py` yourself, but you must set the environment variables
to appropriate values. See `To Configure` above to determine these values.

Alternatively, if you have Python installed but do not have Cython, then
`./configure` will recognize this, and make will build PyClical via the command
```
  ext_name=`PyClical` source_cpp=`PyClical_nocython.cpp` \
  CXX=`$(CXX)` CXXVERSION=`$(CXXVERSION)` CFLAGS=`` \
  CXXFLAGS=`$(EXTCXXFLAGS)` AM_CPPFLAGS=`$(EXTAM_CPPFLAGS)` \
  LDFLAGS=`$(USER_LDFLAGS)` LIBRARIES=`$(LIBS)` \
  $(PYTHON) setup_nocython.py build_ext --inplace
```
with all variables set as above. Again, you can run `pyclical/setup_nocython.py`
yourself, but you must set all of the relevant environment variables to
appropriate values.

In the `pyclical/demos` directory, there is a Python script
`build_pyclical_notebooks.py`. This script builds a set of Jupyter notebooks
that correspond to the PyClical tutorials and some of the demos. Running `make`
will use this Python script to build these notebooks.


Building and running the regression test programs
-------------------------------------------------

To build and run the regression test programs, set the environment variable `CXX`
to indicate your C++ compiler, eg. `g++` for GNU C++, `icpc` or `icpx` for Intel
C++, then run `./configure` as above, and then run `make check`. This builds and
runs the executable files `./test00/test00` to `./test17/test17`. This produces
the intermediate output files `./test00/test00.out` to `./test17/test17.out`,
and the final test output file ./test_runtime/test.out. You can use a parallel
make for `make check` , e.g. `make check -j 6`. This is especially useful on
modern multicore machines.

Warning:If you use too many jobs with parallel make, the compiler will have
problems obtaining enough memory to run efficiently.


To Test
=======

The test_runtime directory
--------------------------

The test runtime directory `./test_runtime` contains sample test output files.

The sample test output files include `eg3.res`, `gfft_test-11.out`, `products-8.out`,
`squaring-11.out` and `transforms-8.out`. There are also 20 versions of the output
of the regression tests. These are described below.

`./test_runtime` also contains the test input file `eg8.txt`. This file is needed
by programming example 8 (reading multivectors from input).


Re-running the regression tests
-------------------------------

The main regression test script is `./test/test.sh`. Once you have built and run
the regression tests via `make check`, you can use this script to re-run the
regression tests. The script `./test/test.sh` reruns tests `./test00/test00` to
`./test17/test17`, using relative pathnames, so it is best to leave `test.sh`
where it is and invoke it using its full path name. This allows it to find
`test00` to `test17`.

The test script `./test/test.sh` takes any number (including zero) of numeric
parameters. Parameters in the range `00` to `17` correspond to coding examples
`./test00/test00` to `./test17/test17`. These examples are run in numerical
order. With zero parameters, all examples from `00` to `17` are run in order.
Many of the examples are run twice - once with `framed_multi<Scalar_T>` and once
with `matrix_multi<Scalar_T>`.

The `./test_runtime` directory contains 20 sample versions of the regression test
results, corresponding to 10 different combinations of configuration parameters,
for two different sets of tests, the complete set of 18 tests, and a subset of 3
tests. The tests were all run on an Intel(R) Core(TM) i7 CPU 870  @ 2.93GHz+ with
```
    Linux 5.13.0-25-generic #26-Ubuntu SMP x86_64
    Kubuntu 21.10
    gcc version 11.2.0 (Ubuntu 11.2.0-7ubuntu2)
    Boost 1.74.0
    Boost Numeric Bindings
    GSL 2.6
    QD 2.3.22
    Cython 0.29.21
    Python 3.9.7
```
The test output file names and corresponding configuration commands are defined
in `./test/config-options.txt` and are:

  1. `test.configure.default.out`:
```
./configure
```
  2. `test.configure.debug-full.out`:
```
./configure --enable-debug=full
```
  3. `test.configure.debug-yes.out`:
```
./configure --enable-debug=yes
```
  4. `test.configure.disable-dependency.out`:
```
./configure --disable-dependency-tracking
```
  5. `test.configure.disable-pyclical.out`:
```
./configure --disable-pyclical
```
  6. `test.configure.prefix-home-opt.out`:
```
./configure --prefix=$HOME/opt
```
  7. `test.configure.eig-bindings.out`:
```
./configure --with-eig=bindings --with-extra-includes=$PATHTO/numeric_bindings
```
  8. `test.configure.eig-bindings-qd.out`:
```
./configure --with-eig=bindings --with-extra-includes=$PATHTO/numeric_bindings \
            --with-qd
```
  9. `test.configure.eig-blaze.out`:
```
./configure --with-eig=blaze
```
 10. `test.configure.eig-blaze-debug-full.out`:
```
./configure --with-eig=blaze --enable-debug=full
```
 11. `test.configure.eig-blaze-debug-yes.out`:
```
./configure --with-eig=blaze --enable-debug=yes
```
 12. `test.configure.eig-blaze-qd.out`:
```
./configure --with-eig=blaze --with-qd
```
For each of the 12 `test.configure.*.out` files, there is a corresponding
`fast-test.configure.*.out` file, making a total of 24 files.

When you run your own test using `./test/test.sh`, you should compare its output
to the output file corresponding to the closest match to the configuration
options you used to build your copy of the GluCat library.

The reason why sample test results corresponding to 12 different combinations
of configuration parameters are included in `test_runtime` is that the test output
strongly depends on the configuration options chosen. In particular:

* If `--with-qd` is chosen, extra tests in `./test00/test00` and `./test11/test11`
  are done using the `dd_real` and `qd_real` scalar types.

* If either `--with-eig=bindings` or `--with-eig=blaze` is chosen the algorithms
  used for the square root, logarithm and inverse trig functions will become much
  more accurate, and most tests in `./test11/test11` will succeed. Even if this
  option is chosen, some tests in `./test11/test11` fail due to insufficient
  accuracy. This is most likely caused by a combination of excessive round off
  and truncation error with respect to the condition numbers of the matrices used
  in calculating these functions.

The tests typically use floating point arithmetic, and `./test00/test00` and
`./test11/test11` in particular also use random number generators. Therefore if
you run the tests using different architecture, compilers or random number
generators, you should expect to have different floating point arithmetic
results, but generally, still within acceptable error tolerances, except as
noted for `./test11/test11` above.

The regression tests `./test00/test00` to `./test17/test17` recognize the program
arguments `--help`, `--no-catch`, and `--verbose`. The `--no-catch` argument
disables the default exception catching behaviour of a regression test, to
allow prgram crashes to be more easily debugged. For `./test00/test00` and
`./test11/test11` the argument `--verbose` produces verbose output essentially by
setting the error tolerance to zero. Verbose output can become quite large.

The test script `./test/test_optional.sh` runs all examples 00 to 17 in order,
with the given parameters as program arguments.


Systematic regression testing
-----------------------------

The `./test` directory contains the scripts `./test/test-all-config-options.sh`,
`./test/fast-test-all-config-options.sh` and
`./test/pyclical-test-all-config-options.sh`. Each of these scripts use the
configuration given by the file `./test/config-options.txt` to build and run the
regression tests and the PyClical doctests once for each set of configuration
options specified in the file. The syntax of each line of the file
`./test/config-options.txt` is simply
```
abbreviation: options
```
This corresponds to running the regression tests and doctests using the
configure command
```
./configure $options
```
then copying the output to `./test_runtime/test.configure.$abbreviation.out`
or `./test_runtime/fast-test.configure.$abbreviation.out` respectively.

For example, for `./test/test-all-config-options.sh`

The line
```
default:
```
specifies the abbreviation `default` and the options `""`, i.e. no options.

This corresponds to running the regression tests and doctests using the
configure command
```
./configure
```
and copying the output to `./test_runtime/test.configure.default.out`

The line
```
disable-dependency:          --disable-dependency-tracking
```
specifies the abbreviation `disable-dependency` and the option
`--disable-dependency-tracking`. This corresponds to running the
regression tests and doctests using the configure command
```
./configure --disable-dependency-tracking
```
and copying the output to `./test_runtime/test.configure.disable-dependency.out`

If you want to test with options different from those in the file
`./test/config-options.txt`, you can write your own file using the same syntax,
for example, `some-other-file.txt`, and invoke your tests with the command
```
config_options_file=some-other-file.txt ./test/test-all-config-options.sh
```
You can also give parameters to `./test/test-all-config-options.sh` and these
are passed to the make command. In particular, invoking (e.g.)
```
./test/test-all-config-options.sh -j 6
```
performs a parallel `make check` for each configuration option, potentially
speeding up the entire testing process.

Rather than running the regression tests in-place and copying the output
directly into `./test_runtime`, the script `./test/test-all-config-options.sh`
produces as many copies of the whole direcory `glucat-0.12.0` as there are lines
in `./test/config-options.txt`, naming them `glucat-0.12.0.1` to `glucat-0.12.0.12`,
in the parent directory of `glucat-0.12.0`. This allows the effect of each set
of configuration options to be directly compared, and also ensures that any
side-effect of a configuration does not affect the test results of another
configuration.

The script `./test/diff-all-config-outputs.sh` compares each relevant test output
file with the corresponding file in `./test_runtime` or `./pyclical`. For example,
line 4 of `./test/config-options.txt`
```
disable-dependency:          --disable-dependency-tracking
```
causes `./test/diff-all-config-outputs.sh` to use diff to compare
`glucat-0.12.0.4/test_runtime/test.configure.disable-dependency.out` to
`glucat-0.12.0/test_runtime/test.configure.disable-dependency.out`, and compare
`glucat-0.12.0.4/pyclical/test.out` to `glucat-0.12.0/pyclical/test.out`.

Each comparison should only produce a line containing the line number of
the configuration being compared: 1 to 12.

The exceptional cases are:

1. If the configuration options have caused any sort of error.
2. Differences in compilers and libraries causing different floating point
   results. This currently occurs with the Intel C++ compiler, which produces
   output different from either the GNU C++ compiler or the Clang compiler.
   For example, QD version 2.3.16 contains an update that fixes a problem with
   `tanh`.

If the output of your systematic tests differs due to a difference in compilers
or libraries, you may want to copy this output to `./test_runtime` to ease future
comparisons. To do so, run the script `./test/copy-all-config-outputs.sh`.

The difference between `./test/test-all-config-options.sh` and
`./test/fast-test-all-config-options.sh` is that the former runs all of the tests
`./test00` to `./test17`, whereas the latter runs only `./test00`, `./test10` and
./test11. Both scripts build and check `./pyclical`. More specifically, the
`./test/test-one-config-option.sh` script runs `make check`, and the
`./test/fast-test-one-config-option.sh` script runs `make fast-check`. For the
definitions of these arguments to `make`, see the file `./Makefile.am.in`. To
compare or copy the output of `./test/fast-test-all-config-options.sh`, use either
`./test/fast-diff-all-config-outputs.sh` or `./test/fast-copy-all-config-outputs.sh`
respectively.

The script `./test/pyclical-test-all-config-options.sh` just builds and checks
`./pyclical` without running any of the other regression tests. To examine the
output of `./test/python-test-all-config-options.sh` just examine all of the
`pyclical/pyclicat-test.check.out` files for errors.


Running the timing (benchmark) tests
------------------------------------

The test program `./gfft_test/gfft_test` takes a parameter `n`, and transforms
larger and larger multivectors within the subalgebra defined by the frame of
the index set `{-n, ..., -1, 1, ..., n}`.

The test program `./products/products` takes a parameter `n`, and runs a timing
test which uses the products `*`, `^`, `%` and `&` to multply larger and larger
multivectors within subalgebras defined by frames limited by the value of `n`.

The test program `./squaring/squaring` takes a parameter `n`, and runs a timing
test which squares larger and larger multivectors within subalgebras defined
by frames limited by the value of `n`.

The test program `./transforms/transforms` takes a parameter `n`, and transforms
larger and larger multivectors within the subalgebras defined by  frames
limited by the value of `n`.

The test script `./test/timing_tests.sh` takes up to 4 numeric parameters.
The command `./test/timing_tests.sh $a $b $c $d` runs
```
 ./products/products $a
 ./squaring/squaring $b
 ./gfft_test/gfft_test $c
 ./transforms/transforms $d
```
The default is:
```
 ./products/products 8
 ./squaring/squaring 11
 ./gfft_test/gfft_test 11
 ./transforms/transforms 8
```
The sample timing test results in `./test_runtime` are from programs
built and run using the configure command:
```
./configure --with-eig=bindings --with-extra-includes=$PATHTO/numeric_bindings \
            --with-qd
```
on `Intel(R) Core(TM) i7 CPU 870  @ 2.93GHz+` with
```
    Linux 5.15.0-27-generic #28-Ubuntu SMP x86_64
    Kubuntu 22.04 LTS
    gcc version 11.2.0 (Ubuntu 11.2.0-19ubuntu1)
    Boost 1.74.0
    GSL 2.7.1
    QD 2.3.23
```

Testing PyClical
----------------
Once you have built PyClical, run the doctests. In `python3` or `ipython3`, etc.:
```
 >>> import PyClical
 >>> PyClical._test()
 TestResults(failed=0, attempted=647)
 >>> quit()
```
Alternatively, in the directory `pyclical`, run the script `test.py` using:
```
ipython3 --classic --no-banner < test.py
```
and compare the output with `test.out`.

In the directory `pyclical/demos`, the Python script files `pyclical_demo.py` and
`sqrt_log_demo.py` have corresponding output files `pyclical_demo.out` and
`sqrt_log_demo.out`. Run these two Python script files and compare their output to
the contents of the two output files. See the file `README.md` under "Using the
PyClical extension module with Python" for instructions on how to run these files.


To Install
==========

Once you are satisfied that GluCat works, you can run `make install`.
This install the headers from `./glucat` into the header installation directory,
`$PREFIX/include`, which defaults to `/usr/local/include`. If you have chosen to
build PyClical, `make install` also installs the file `PyClical.so` into
a directory under `$EPREFIX` that is determined by the version of Python you are
running, for example `/usr/local/python3.7/site-packages/`.

Running `make install` will also install some documentation. Specifically, the
GNU top level documentation files `AUTHORS.md`, `ChangeLog`, `COPYING`,
`glucat.lsm`, `INSTALL.md`, `NEWS`, `README.md` and `TODO.md` are installed into
`$DATAROOTDIR/doc/glucat`. This defaults to `$PREFIX/share/doc/glucat`, i.e.
`/usr/local/share/doc/glucat`. The PyClical demos and notebooks are installed by
default into the directory `$DATAROOTDIR/pyclical/demos`. This defaults to
`$PREFIX/share/pyclical/demos`, i.e. `/usr/local/share/pyclical/demos`. You can
change this directory by using the configuration option `--with-demo-dir=DIR`.

You will need permission to update the installation directories, so you may need
to use `sudo`, login as `root`, or `su` to `root` before you run `make install`.


List of Successful Builds
=========================

GluCat 0.12.0 with PyClical has so far been built and tested using:

 1) Pensieri:
    4 core `Intel(R) Core(TM) i7 CPU 870  @ 2.93GHz` with
    ```
    Linux 5.15.0-27-generic #28-Ubuntu SMP x86_64
    Kubuntu 22.04 LTS
    Boost 1.74.0
    Boost Numeric Bindings
    GSL 2.7.1
    QD 2.3.23
    Cython 0.29.28
    Python 3.10.4
    Numpy 1.21.5
    Matplotlib 3.5.1
@   Mayavi2 4.7.4
    VTK 9.1.0
    Doxygen 1.9.1
    pdfTeX 3.141592653-2.6-1.40.22 (TeX Live 2022/dev/Debian)
    ```

    `./test/test-all-config-options.sh`:
    All 12 configuration commands corresponding to each of the 12
    `test.configure*.out` files in `./test_runtime`
    tested with the following compiler versions:
    1) `gcc version 11.2.0 (Ubuntu 11.2.0-19ubuntu1)`
    2) `clang version 14.0.0 (14.0.0-1ubuntu1)`
@   3) `icpx version 2022.0.0 (2022.0.0.20211123)`

    `./test/fast-test-all-config-options.sh`:
    All 12 configuration commands corresponding to each of the 12
    `fast-test.configure*.out` files in `./test_runtime`
    tested with the following compiler versions:
@   1) `gcc version 7.5.0 (Ubuntu 7.5.0-6ubuntu4)`
@   2) `gcc version 8.5.0 (Ubuntu 8.5.0-0ubuntu4)`
@   3) `gcc version 9.4.0 (Ubuntu 9.4.0-3ubuntu1)`
@   4) `gcc version 10.3.0 (Ubuntu 10.3.0-11ubuntu1)`
    5) `clang version 9.0.1 (9.0.1-16.1ubuntu1)`
@   6) `clang version 11.0.1 (11.0.1-2ubuntu5)`
@   7) `clang version 12.0.1 (12.0.1-8build1)`

 2) Pensieri (VirtualBox):
    Virtual 1 core `Intel(R) Core(TM) i7 CPU 870 @ 2.93GHz` with
    ```
    Linux 5.16.8-1-default #1 SMP 2022
    openSUSE Tumbleweed Snapshot TBD 20211215
    g++ (SUSE Linux) 11.2.1 20220103
    Boost 1.78.0
    Boost Numeric Bindings
    GSL 2.6-6.4
    QD 2.3.22-1.13
    Python 3.8.12
    Numpy 1.21.4
    Matplotlib 3.5.1
    Mayavi2 4.7.4
    VTK 9.1.0-1.8
    Doxygen 1.9.2-1.2
    pdfTeX 3.141592653-2.6-1.40.22 (TeX Live 2021/TeX Live for SUSE Linux)
    ```
    `./test/fast-test-all-config-options.sh`
    All 12 configuration commands corresponding to each of the 12
    `fast-test.configure*.out` files in `./test_runtime`

 3) NCI Gadi:
    48 core `Intel(R) Xeon(R) Platinum 8274 CPU @ 3.20GHz` with
    ```
    Linux  4.18.0-348.2.1.el8.nci.x86_64 SMP x86_64
    Rocky Linux release 8.4 (Green Obsidian)
    ```
    1) `icpc version 2021.5.0`
    2) `icpx version 2021.5.0`
    ```
    Boost 1.77.0
    Boost Numeric Bindings
    Cython 0.29.24
    Python 3.8.8
    Numpy 1.20.2
    ```
    `./test/fast-test-all-config-options.sh`
    All 12 configuration commands corresponding to each of the 12
    `fast-test.configure*.out` files in `./test_runtime`

 4) CoCalc:
    Virtual 8 core `Intel(R) Xeon(R) CPU @ 2.80GHz` with
    ```
    Linux 5.11.0-1020-gcp #22~20.04.1-Ubuntu SMP x86_64
    Ubuntu 18.04.5 LTS
    gcc version 7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04)
    Boost 1.65.1
    Boost Numeric Bindings
    Cython 0.29.21
    Python 3.6.9
    Numpy 1.18.5
    CXXFLAGS=`-I/home/user/usr/local/include`
    USER_LDFLAGS=`-L/home/user/usr/local/lib`
    LD_LIBRARY_PATH=`/home/user/usr/local/lib`
    ```
    `./test/fast-test-all-config-options.sh`
    All 12 configuration commands corresponding to each of the 12
    `fast-test.configure*.out` files in `./test_runtime`

 5) AWS Graviton:
    Virtual 4 core `ARM Cortex-A72 Model 3 (AWS Graviton A1 image)` with
    ```
    Linux 5.13.0-1019-aws #21~20.04.1-Ubuntu SMP aarch64
    Ubuntu 20.04.4 LTS
    g++ (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
    Boost 1.71.1
    Boost Numeric Bindings
    GSL Version: 2.5+dfsg-6build1
    QD Version: 2.3.22+dfsg.1-3build1
    Cython 0.29.14
    Python 3.8.10
    ```
    `./test/fast-test-all-config-options.sh`
    All 12 configuration commands corresponding to each of the 12
    `fast-test.configure*.out` files in `./test_runtime`

                
Notes on software versions
==========================

Building the documentation requires recent versions of both `doxygen` and
`latex`, (e.g. `texlive-2021`).


Using `gcc` version 11.X with Cython from any Conda environment that uses
`libstdc++` from `gcc` version 10.X and earlier (e.g. version 4.8.3) results
in incompatible `libstdc++` versions.
see https://github.com/conda/conda/issues/10757
and https://github.com/cython/cython/issues/4218
and https://github.com/stan-dev/pystan/issues/294


PyClical is now compatible with Python 3 and is backwards incompatible
with Python 2.


The PyClical plotting demos use Numpy and either Matplotlib or Mayavi2.
The versions used need to be compatible with each other and with Python.
For Mayavi2 the versions of VTK and TVTK used also need to be compatible.


The use of Mayavi2 4.7.4 on Kubuntu 22.04 depends on the following packages:

Ubuntu packages
---------------
```
envisage 6.0.1
python3-attr 21.2.0-1
python3-configobj 5.0.6-5
python3-idna 3.3-1
python3-importlib-metadata 4.6.4-1
python3-matplotlib 3.5.1-2
python3-numpy 1.21.5-1build2
python3-pygments 2.11.2
python3-setuptools 59.6.0-1.2
python3-zipp 1.0.0-3
```
pip packages
------------
```
aiohttp 3.8.1
aiosignal 1.2.0
apptools 5.1.0
async-timeout 4.0.2
charset-normalizer 2.0.12
frozenlist 1.3.0
multidict 6.0.2
pyface 7.4.1
traits 6.3.2
traitsui 7.3.1
wslink 1.6.4
yarl 1.7.2
```

As of May 2022, installation of `vtk 9.1.0` as a pip package on Kubuntu 22.04 does not work
because a Python 3.10 compatible wheel is not available.
The `vtk 9.1.0` package can instead be installed using
```
pip install --find-links https://wheels.pyvista.org/ pyvista
```
as per [Pyvista instructions](https://github.com/pyvista/pyvista/discussions/2064).

After installing `vtk 9.1.0`, the `mayavi 4.7.4` package can be installed from source
using the following instructions
```
git clone https://github.com/enthought/mayavi.git
cd mayavi
python setup.py install --prefix=~/.local
```
See also https://pypi.org/project/mayavi/


The use of Mayavi2 4.7.4 with Python 3.8.12 on openSUSE Tumbleweed
requires the following RPM packages:
```
mayavi 4.7.4-1.2
python3-tvtk 4.7.4-1.2
python3-vtk 9.1.0-1.3
python38-apptools 4.5.0-1.10
python38-configobj 5.0.6-3.12
python38-envisage 6.0.1-1.2
python38-importlib-metadata-4.8.2-1.1
python38-importlib-resources 5.4.0-1.1
python38-numpy 1.21.2-4.1
python38-pyface 7.3.0-3.4
python38-Pygments 2.9.0-2.1
python38-setuptools 57.4.0-1.2
python38-six 1.16.0-2.2
python38-traits 6.3.1-1.1
python38-traitsui 7.1.1-1.6
python38-zipp 3.6.0-1.1
```
The use of `jupyter-notebook-6.2.0` on Kubuntu 21.10 depends on the following
packages:

Ubuntu packages
---------------
```
ipython3 7.20.0-1
jupyter-nbformat 5.1.2-1
jupyter-notebook 6.2.0-1
python3-ipykernel 5.4.3-1
python3-ipython-genutils 0.2.0-4
python3-nbconvert 5.6.1-3
python3-nbformat 5.1.2-1
python3-pexpect 4.8.0-2ubuntu1
python3-send2trash 1.6.0~b1+git20210122.2eb3242-1
python3-traitlets 5.0.5-1
```
pip packages
------------
```
ipywidgets 7.6.5
jupyterlab-widgets 1.0.2
widgetsnbextension 3.5.2
```

The use of `jupyter-notebook-6.4.6` with Python 3.8.12 on openSUSE Tumbleweed
requires the following packages:

RPM packages:
-------------
```
jupyter-widgetsnbextension-3.5.2-1.1.noarch
python38-ipywidgets-7.6.5-1.1.noarch
python38-jupyterlab-widgets-1.0.2-1.1.noarch
python38-widgetsnbextension-3.5.2-1.1.noarch
```
pip packages
------------
```
jupyter-nbextensions-configurator-0.4.1
nbconvert-5.6.1
send2trash 1.8.0
```
Version incompatibilities discovered in testing `pyclical/demos`:

 1. Using Mayavi2 4.7.4 with VTK 9.1.0 on Kubuntu 21.10 results in the following
    warning message when running `pyclical/demos/plotting_demo_dialog.py` and
    `pyclical/demos/plotting_demo_mayavi.py`:
```
WARNING: Imported VTK version (9.1) does not match the one used
         to build the TVTK classes (9.0). This may cause problems.
         Please rebuild TVTK.
```
    The warning can be ignored. The plotting demos run normally.

 2. Using `jupyter-notebook-6.4.6` with Python 3.8.12 on openSUSE Tumbleweed
    with `nbconvert-6.0.7` results in warning messages such as the following:
```
Config option `template_path` not recognized by `LenvsHTMLExporter`.
Did you mean one of: `extra_template_paths, template_name, template_paths`?
```
    This is due to an incompatibility between `nbconvert-6.0`+ and `jupyter_latex_envs`.
    The workaround is to use `nbconvert-5.6.1`.
    See https://github.com/ipython-contrib/jupyter_contrib_nbextensions/issues/1529
    and https://github.com/jfbercher/jupyter_latex_envs/pull/58

 3. Using `jupyter-notebook-6.4.6` with Python 3.8.12 on openSUSE Tumbleweed
    results in the following warning message:
```
404 GET /notebooks/biblio.bib (127.0.0.1): No such file or directory: biblio.bib
```
    The warning can be ignored. The notebooks work normally.


The following bugs and workarounds apply to earlier versions of GluCat,
and may still be applicable to GluCat 0.12.0, but have not been checked
for this version.


 1. Using Mayavi2 4.7.1 with VTK 7.1.1 as per Kubuntu 21.04 yields two bugs
    likely caused by version mismatch:
    1. Running `pyclical/demos/plotting_demo_mayavi.py` results in:
```
Warning: In ./Common/ExecutionModel/vtkAlgorithm.cxx, line 1419
vtkGlyph3D (): Attempt to get connection index 0 for input port 0,
which has 0 connections.
```
       A similar issue:
       https://vtk.org/pipermail/vtk-developers/2014-May/014965.html
    2. The 3D plots produced by pyclical/demos/plotting_demo_mayavi.py have
       an incorrect z-order. A similar issue:
       https://github.com/enthought/mayavi/issues/656

 2. Using Mayavi2 4.7.2 with VTK 9.0.1 and Python 3.8 on openSUSE
    Tumbleweed results in the following warning message when running
    `pyclical/demos/plotting_demo_mayavi.py`
```
/usr/lib64/python3.8/site-packages/vtkmodules/numpy_interface/algorithms.py:
209: SyntaxWarning: `is` with a literal. Did you mean `==`?
  if max_dims is 1:
```
 3. GluCat needs an include library which either defines the macro `isnan` or
    defines `std::isnan`. The C++ standard library included with `gcc` 4.5.2 and
    above defines `std::isnan`.

 4. Cython versions earlier than 0.15 do not build PyClical correctly,
    because PyClical uses generators, which were only implemented with
    Cython 0.15.

 5. Cython versions to and including 0.16 do not build PyClical correctly
    for C++11. If you try to use `g++` with `-std=c++11` you will see
    an error message like:
```
In function `void __Pyx_RaiseArgtupleInvalid(...)`:
error: unable to find string literal operator ‘operator PY_FORMAT_SIZE_T’
```
    See https://github.com/cython/cython/pull/109

    The workaround is to edit `PyClical.cpp` and put a space before and after each
    occurrence of `PY_FORMAT_SIZE_T`. This was fixed some time after Cython 0.16.

 6. GluCat will not work with QD versions earlier than 2.3.10, because older
    versions of QD do not zero-initialize `dd_real` and `qd_real` as required by
    `ublas::clear()`.

 7. With clang++ 3.2, building PyClical results in the warning
```
clang: warning: argument unused during compilation: '--param ssp-buffer-size=4'
```
    This is harmless, and was fixed after Clang version 3.2.
    See http://llvm.org/bugs/show_bug.cgi?id=15327

 8. The following version incompatibility was observed during testing with
    GluCat 0.8.2:

    With `g++` 5.3.1 and Boost 1.53.0 or Boost 1.55.0, the header file
    `<boost/smart_ptr/shared_ptr.hpp>` generates multiple warnings of the form:
```
warning: ‘template<class> class std::auto_ptr’ is deprecated [-Wdeprecated-declarations]
```
    This does not occur with `g++` 4.8.5 and Boost 1.53.0, 1.55.0 or 1.61.0,
    nor with `g++` 5.3.1 and Boost 1.58.0 or Boost 1.61.0.

    See https://svn.boost.org/trac/boost/ticket/11411
    and https://svn.boost.org/trac/boost/ticket/11622
