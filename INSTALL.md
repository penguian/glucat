INSTALL for GluCat 0.13.0 with PyClical
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
> git clone git@github.com:penguian/glucat.git glucat-0.13.0
> cd glucat-0.13.0
> make -f admin/Makefile.common cvs
```
This results in a directory structure that includes glucat-0.13.0/configure,
allowing you to make and install GluCat in the same way as if you had downloaded
and unzipped the tarball glucat-0.13.0.tar.gz.


Directory Structure
===================

Once you have downloaded, unzipped and untarred the source code, or followed
the instructions above to install from Git clone, you should have a directory,
glucat-0.13.0. Under this directory you should see a number of subdirectories,
including `./admin`, `./doc`, `./glucat`, `./gfft_test`, `./products`,
`./pyclical`, `./squaring`, `./test`, `./test_move`, `./test_runtime`,
`./testxx`, and `./transforms`.

The following instructions are meant to be used with `glucat-0.13.0` as the
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
installs the GluCat header files, as well as the basic documentation.

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

 1. `cd` to the `glucat-0.13.0` directory containing the source code and type
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
                          (no|blaze) [default=no]
```
This option is used to control `_GLUCAT_USE_EIGENVALUES` and determine which
libraries to use. ARG can be `no` or `blaze`. The default is `no`.

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

If `_GLUCAT_USE_BLAZE` is defined, `glucat/matrix_imp.h` includes
`<blaze/Math.h>` and related Blaze include files, as per the Blaze template
library. To use this library, you will need to install it yourself, preferably
from https://bitbucket.org/blaze-lib/blaze/src/master/

To compile your own programs using the GluCat library, to detect and correctly
handle negative real eigenvalues in the `sqrt()` and `log()` functions, your
Makefile needs to pass the flag `-D_GLUCAT_USE_EIGENVALUES` to the C++ compiler,
as well as one of the following choices of flags, and the corresponding header
files and libraries must be usable.

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
  all_includes="$(all_includes)" LDFLAGS=`$(USER_LDFLAGS)` LIBRARIES=`$(LIBS)` \
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
  all_includes="$(all_includes)" LDFLAGS=`$(USER_LDFLAGS)` LIBRARIES=`$(LIBS)` \
  $(PYTHON) setup_nocython.py build_ext --inplace
```
with all variables set as above. Again, you can run `pyclical/setup_nocython.py`
yourself, but you must set all of the relevant environment variables to
appropriate values.

Please note that the ability to build PyClical without Cython is deprecated and
will be removed in future versions.

In the `pyclical/demos` directory, there is a Python script
`build_pyclical_notebooks.py`. This script builds a set of Jupyter notebooks
that correspond to the PyClical tutorials and some of the demos. Running `make`
will use this Python script to build these notebooks.


Building and running the regression test programs
-------------------------------------------------

To build and run the regression test programs, set the environment variable `CXX`
to indicate your C++ compiler, eg. `g++` for GNU C++, `icpc` or `icpx` for Intel
C++, then run `./configure` as above, and then run `make check`. This builds and
runs the executable files `./test_move/test_move` and `./test00/test00` to
`./test17/test17`. This produces the intermediate output files
`./test_move/test_move.out` and `./test00/test00.out` to `./test17/test17.out`,
and the final test output file `./test_runtime/test.out`. You can use a parallel
make for `make check` , e.g. `make check -j 4`. This is especially useful on
modern multicore machines.

Warning: If you use too many jobs with parallel make, the compiler will have
problems obtaining enough memory to run efficiently.


To Test
=======

The test_runtime directory
--------------------------

The test runtime directory `./test_runtime` contains sample test output files.

The sample test output files include `eg3.res`, `gfft_test-11.out`, `products-8.out`,
`squaring-11.out` and `transforms-8.out`. There are also 22 versions of the output
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

The `./test_runtime` directory contains 22 sample versions of the regression test
results, corresponding to 11 different combinations of configuration parameters,
for two different sets of tests: the complete set of 18 tests, and a subset of 3
tests. The tests were all run on an 8 core 
`AMD Ryzen 7 8840HS w/ Radeon 780M Graphics` @ 3.3 GHz with
```
    Linux 6.11.0-14-generic #15-Ubuntu SMP 2025
    Kubuntu 24.10
    g++ 14.2.0 (Ubuntu 14.2.0-4ubuntu2)
    Blaze 3.9.0
    Boost 1.83.0
    GSL 2.8
    QD 2.3.23
    Cython 3.0.11
    Python 3.12.7
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
  7. `test.configure.qd.out`:

```
./configure --with-qd
```
  8. `test.configure.eig-blaze.out`:

```
./configure --with-eig=blaze
```
  9. `test.configure.eig-blaze-debug-full.out`:

```
./configure --with-eig=blaze --enable-debug=full
```
 10. `test.configure.eig-blaze-debug-yes.out`:

```
./configure --with-eig=blaze --enable-debug=yes
```
 11. `test.configure.eig-blaze-qd.out`:

```
./configure --with-eig=blaze --with-qd
```
For each of the 11 `test.configure.*.out` files, there is a corresponding
`fast-test.configure.*.out` file, making a total of 22 files.

When you run your own test using `./test/test.sh`, you should compare its output
to the output file corresponding to the closest match to the configuration
options you used to build your copy of the GluCat library.

The reason why sample test results corresponding to 11 different combinations
of configuration parameters are included in `test_runtime` is that the test output
strongly depends on the configuration options chosen. In particular:

* If `--with-qd` is chosen, extra tests in `./test00/test00` and `./test11/test11`
  are done using the `dd_real` and `qd_real` scalar types.

* If `--with-eig=blaze` is chosen the algorithms used for the square root,
  logarithm and inverse trig functions will become much more accurate, and most
  tests in `./test11/test11` will succeed. Even if this option is chosen, some
  tests in `./test11/test11` fail due to insufficient accuracy. This is most
  likely caused by a combination of excessive round off and truncation error
  with respect to the condition numbers of the matrices used in calculating
  these functions.

The tests typically use floating point arithmetic, and `./test00/test00` and
`./test11/test11` in particular also use random number generators. Therefore if
you run the tests using different architecture, compilers or random number
generators, you should expect to have different floating point arithmetic
results, but generally, still within acceptable error tolerances, except as
noted for `./test11/test11` above.

The regression tests `./test00/test00` to `./test17/test17` recognize the program
arguments `--help`, `--no-catch`, and `--verbose`. The `--no-catch` argument
disables the default exception catching behaviour of a regression test, to
allow program crashes to be more easily debugged. For `./test00/test00` and
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
./test/test-all-config-options.sh -j 4
```
performs a parallel `make check` for each configuration option, potentially
speeding up the entire testing process.

Rather than running the regression tests in-place and copying the output
directly into `./test_runtime`, the script `./test/test-all-config-options.sh`
produces as many copies of the whole directory `glucat-0.13.0` as there are lines
in `./test/config-options.txt`, naming them `glucat-0.13.0.1` to `glucat-0.13.0.11`,
in the parent directory of `glucat-0.13.0`. This allows the effect of each set
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
`glucat-0.13.0.4/test_runtime/test.configure.disable-dependency.out` to
`glucat-0.13.0/test_runtime/test.configure.disable-dependency.out`, and compare
`glucat-0.13.0.4/pyclical/test.out` to `glucat-0.13.0/pyclical/test.out`.

Each comparison should only produce a line containing the line number of
the configuration being compared: 1 to 11.

The exceptional cases are:

1. If the configuration options have caused any sort of error.
2. Differences in compilers and libraries causing different floating point
   results. A compiler difference previously occurred with the Intel C++ compiler,
   which produced output different from either the GNU C++ compiler or the Clang
   compiler. As an example of a library difference, QD version 2.3.16 contains an
   update that fixes a problem with `tanh`.

If the output of your systematic tests differs due to a difference in compilers
or libraries, you may want to copy this output to `./test_runtime` to ease future
comparisons. To do so, run the script `./test/copy-all-config-outputs.sh`.

The difference between `./test/test-all-config-options.sh` and
`./test/fast-test-all-config-options.sh` is that the former runs all of the tests
`./test00` to `./test17`, whereas the latter runs only `./test00`, `./test10` and
`./test11`. Both scripts build and check `./pyclical`. More specifically, the
`./test/test-one-config-option.sh` script runs `make check`, and the
`./test/fast-test-one-config-option.sh` script runs `make fast-check`. For the
definitions of these arguments to `make`, see the file `./Makefile.am.in`. To
compare or copy the output of `./test/fast-test-all-config-options.sh`, use either
`./test/fast-diff-all-config-outputs.sh` or `./test/fast-copy-all-config-outputs.sh`
respectively.

The script `./test/pyclical-test-all-config-options.sh` just builds and checks
`./pyclical` without running any of the other regression tests. To examine the
output of `./test/pyclical-test-all-config-options.sh` just examine all of the
`pyclical/pyclical-test.check.out` files for errors.


Running the timing (benchmark) tests
------------------------------------

The test program `./gfft_test/gfft_test` takes a parameter `n`, and transforms
larger and larger multivectors within the subalgebra defined by the frame of
the index set `{-n, ..., -1, 1, ..., n}`.

The test program `./products/products` takes a parameter `n`, and runs a timing
test which uses the products `*`, `^`, `%` and `&` to multiply larger and larger
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
./configure --with-eig=blaze --with-qd
```
on an 8 core `AMD Ryzen 7 8840HS w/ Radeon 780M Graphics` @ 3.3 GHz with
```
    Linux 6.11.0-14-generic #15-Ubuntu SMP 2025
    Kubuntu 24.10
    g++ 14.2.0 (Ubuntu 14.2.0-4ubuntu2)
    Blaze 3.9.0
    Boost 1.83.0
    GSL 2.8
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

GluCat 0.13.0 with PyClical has so far been built and tested using:

 1) Tempesta:
    8 core `AMD Ryzen 7 8840HS w/ Radeon 780M Graphics` @ 3.3 GHz with

    ```
    Linux 6.11.0-14-generic #15-Ubuntu SMP 2025
    Kubuntu 24.10
    Boost 1.83.0
    GSL 2.8
    QD 2.3.23
    Cython 3.0.11
    Python 3.12.7

    Numpy 1.21.5
    Matplotlib 3.5.1
    Mayavi2 4.8.1
    VTK 9.1.0
    Doxygen 1.9.8
    pdfTeX 3.141592653-2.6-1.40.25 (TeX Live 2023/Debian)
    ```

    `./test/test-all-config-options.sh`:
    All 11 configuration commands corresponding to each of the 11
    `test.configure*.out` files in `./test_runtime`
    tested with the following compiler versions:

    1) `g++ 14.2.0 (Ubuntu 14.2.0-4ubuntu2)`
    2) `Ubuntu clang version 19.1.1 (1ubuntu1)`

 2) Pensieri:
    4 core `Intel(R) Core(TM) i7 CPU 870  @ 2.93GHz` with

    ```
    Linux 6.11.0-13-generic #14-Ubuntu SMP 2024
    Kubuntu 24.10
    Blaze 3.9.0
    Boost 1.83.0
    GSL 2.8
    QD 2.3.23
    Cython 3.0.11
    Python 3.12.7

    Numpy 1.21.5
    Matplotlib 3.5.1
    Mayavi2 4.8.1
    VTK 9.1.0
    Doxygen 1.9.8
    pdfTeX 3.141592653-2.6-1.40.25 (TeX Live 2023/Debian)
    ```

    `./test/fast-test-all-config-options.sh`:
    All 11 configuration commands corresponding to each of the 11
    `fast-test.configure*.out` files in `./test_runtime`
    tested with the following compiler versions:

    1) `g++ 14.2.0 (Ubuntu 14.2.0-4ubuntu2)`
    2) `Ubuntu clang version 19.1.1 (1ubuntu1)`
    3) `Intel(R) oneAPI DPC++/C++ Compiler 2025.0.4 (2025.0.4.20241205)`

    Note: One test in test_runtime/fast-test.configure.eig-blaze-qd.out
    fails due to a difference between Intel and AMD x86-64 floating
    point arithmetic.

 3) Vincitor (Pensieri running VirtualBox):
    Virtual 1 core `Intel(R) Core(TM) i7 CPU 870 @ 2.93GHz` with

    ```
    Linux 6.12.6-1-default #1 SMP 2024
    openSUSE Tumbleweed Release 20241224
    g++ (SUSE Linux) 14.2.1 20241007
    Blaze 3.8.2 (blaze-devel cblas-devel libcblas3)
    Boost 1.86.0
    Cython version 3.0.12 (python313-Cython)
    GSL 2.8
    QD 2.3.24 (qd-devel libqd0)
    Python 3.11.11
    Numpy 2.1.3
    Matplotlib 3.10.0
    Mayavi2 4.8.2
    VTK 9.4.1
    Doxygen 1.12.0 (doxygen)
    pdfTeX 3.141592653-2.6-1.40.26 (TeX Live 2024/TeX Live for SUSE Linux)
    ```
    `./test/fast-test-all-config-options.sh`
    All 11 configuration commands corresponding to each of the 11
    `fast-test.configure*.out` files in `./test_runtime`

    Note: One test in test_runtime/fast-test.configure.eig-blaze-qd.out
    fails due to a difference between Intel and AMD x86-64 floating
    point arithmetic.

 4) AWS Graviton:
    Virtual 4 core `ARM Cortex-A72 Model 3 (AWS Graviton A1 image)` with

    ```
    Linux 5.15.0-1077-aws #84~20.04.1-Ubuntu SMP aarch64
    Ubuntu 20.04.6 LTS
    g++ (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0
    Blaze 3.9.0
    Boost 1.71.0
    GSL Version: 2.5+dfsg-6build1
    QD Version: 2.3.22+dfsg.1-3build1
    Cython 0.29.14
    Python 3.8.10
    ```
    `./test/fast-test-all-config-options.sh`
    All 11 configuration commands corresponding to each of the 11
    `fast-test.configure*.out` files in `./test_runtime`.
    Note: All tests involving `long double` give different answers
    from the same tests on x86-64 hardware, because `long double` is
    128 bits on ARM 64 hardware vs 80 bits on x86-64.
    https://github.com/ARM-software/abi-aa/blob/main/aapcs64/aapcs64.rst#arithmetic-types

 5) CoCalc:
    Virtual 2 core `Intel(R) Xeon(R) CPU @ 2.80GHz` with

    ```
    Linux 5.15.0-1046-gcp #54~20.04.1-Ubuntu SMP x86 64
    Ubuntu 22.04.4 LTS
    g++ (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
    Blaze 3.9.0
    Boost 1.74.0
    Cython 0.29.30
    Python 3.10.12
    Numpy 1.23.5
    CXXFLAGS=`-I/home/user/include`
    USER_LDFLAGS=`-L/home/user/lib`
    LD_LIBRARY_PATH=`/home/user/lib`
    ```
    `./test/fast-test-all-config-options.sh`
    All 11 configuration commands corresponding to each of the 11
    `fast-test.configure*.out` files in `./test_runtime`

    Note: One test in test_runtime/fast-test.configure.eig-blaze-qd.out
    fails due to a difference between Intel and AMD x86-64 floating
    point arithmetic.

 6) GitHub codespaces:
    AMD EPYC 7763 64-Core Processor with

    ```
    Linux codespaces-50eecb 6.5.0-1025-azure #26~22.04.1-Ubuntu SMP x86_64
    Ubuntu 20.04.6 LTS
    g++ (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0
    Blaze 3.8.2
    Boost 1.71.1
    GSL Version: 2.5+dfsg-6build1
    QD Version: 2.3.22+dfsg.1-3build1
    Cython 0.29.14
    Python 3.12.1
    ```
    `./test/fast-test-all-config-options.sh`
    All 11 configuration commands corresponding to each of the 11
    `fast-test.configure*.out` files in `./test_runtime`.

     Note: all configuration commands other than
     ```
     ./configure --disable-pyclical
     ```
     result in
     ```
     configure: WARNING: Cannot build using Cython.
     configure: WARNING: Cannot build PyClical.
     ```

Notes on software versions
==========================

Building the documentation requires recent versions of both `doxygen` and
`latex`, (e.g. `texlive-2024`).


PyClical is now compatible with Python 3 and is backwards incompatible
with Python 2.


The PyClical plotting demos use Numpy and either Matplotlib or Mayavi2.
The versions used need to be compatible with each other and with Python.
For Mayavi2 the versions of VTK and TVTK used also need to be compatible.

The use of Mayavi2-based plotting demos on Kubuntu 24.10 is achieved via the
following procedure:

1. Install Conda from the Anaconda distribution.
   https://www.anaconda.com/download

2. Run `pyclical/demos/kubuntu-24-conda-install-mayavi.sh`

3. Run `./configure` with your preferred options.

4. Run make clean.

5. Run make.

6. Change directory to `pyclical/demos`.

7. Run `./kubuntu-mayavi-env.sh` before running either
   `./plotting_demo_dialog.py` or `./plotting_demo_mayavi.py`.

For Kubuntu 25.10, Conda is not needed, but you still need to run
`./kubuntu-mayavi-env.sh` before running either
`./plotting_demo_dialog.py` or `./plotting_demo_mayavi.py`.

The use of Mayavi2 4.8.2 with Python 3.11.11 on openSUSE Tumbleweed involves
the installation following and other related RPM packages:

    ```
    mayavi 4.8.2
    python3-vtk 9.4.1
    python311-apptools 5.3.0
    python311-configobj 5.0.9
    python311-envisage 6.1.1
    python311-importlib-metadata-8.6.1
    python311-importlib-resources 6.1.1
    python311-numpy 2.1.3
    python311-pyface 8.0.0
    python311-Pygments 2.19.1
    python311-setuptools 75.8.0
    python311-six 1.17.0
    python311-traits 6.4.3
    python311-traitsui 8.0.0
    python311-zipp 3.21.0
    ```
The exact versions needed will change as openSUSE Tumbleweed is updated.
