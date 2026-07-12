INSTALL for GluCat 0.98a2 with PyClical
========================================

Prerequisites: Before You Begin
===============================

GluCat uses the C++ Standard Library, the Boost Library, and either Eigen or
Armadillo for linear algebra. The PyClical Python extension module is built
using Cython and Python. Make sure that you have all of these installed and
working before attempting to use GluCat with PyClical.

The PyClical plotting demos (`plotting_demo_mayavi.py`,
`plotting_demo_dialog.py`) additionally require a specially configured
Mayavi/VTK environment. On x86-64, this is provided via a Conda environment;
on ARM aarch64 (Asahi Linux), via a system venv. See "Setting Up the Mayavi
Plotting Environment" below for full instructions.

Use the instructions at http://www.boost.org/more/download.html to obtain
the Boost Library. Make sure that you are able to build the Boost library with
the same C++ compiler version as you will be using to build programs that will
use the GluCat library.

GluCat uses the Eigen C++ library by default for its linear algebra backend.
You can obtain Eigen from http://eigen.tuxfamily.org/. Alternatively, you can
configure GluCat to use the Armadillo C++ library, available from
http://arma.sourceforge.net/.

The default configuration of GluCat, as well as most combinations of
configuration options, require a compiler that supports the C++ 2023 standard.

Note: The legacy Intel Classic C++ compilers (`icc`/`icpc`) are deprecated and will be removed in a future release. Users in Intel hardware environments should compile using the modern oneAPI LLVM-based C++ compilers (`icx`/`icpx`) instead.

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
> git clone git@github.com:penguian/glucat.git glucat-0.98a2
> cd glucat-0.98a2
> make -f admin/Makefile.common bootstrap
```
This results in a directory structure that includes glucat-0.98a2/configure,
allowing you to make and install GluCat in the same way as if you had downloaded
and unzipped the tarball glucat-0.98a2.tar.gz.


Directory Structure
===================

Once you have downloaded, unzipped and untarred the source code, or followed
the instructions above to install from Git clone, you should have a directory,
glucat-0.98a2. Under this directory you should see a number of subdirectories,
including `./admin`, `./doc`, `./glucat`, `./gfft_test`, `./products`,
`./pyclical`, `./squaring`, `./test`, `./test_coverage`, `./test_doctest`,
`./test_move`, `./test_runtime.*`, `./test00` to `./test19`, and `./transforms`.

The following instructions are meant to be used with `glucat-0.98a2` as the
current directory.


Installation Summary
====================

There are several different types of program that you can make using the GluCat
library with PyClical. These are:

 1. The PyClical Python extension module.
 2. The PyClical demos.
 3. The PyClical plotting demos (requires a special Mayavi/VTK environment).
 4. The two types of test programs (timing and regression).
 5. Code coverage and unit tests.
 6. Your own programs, written in C++.

The following items can be installed:

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


The Different Types of Programs
===============================

The PyClical Python extension module
------------------------------------

GluCat includes a Python extension module called PyClical.

The subdirectory `./pyclical` contains source code for the PyClical extension
module.

PyClical is written in C++ and Cython, and is defined in the files
`pyclical/glucat.pxd`, `pyclical/PyClical.h`, `pyclical/PyClical.pxd`, and
`pyclical/PyClical.pyx`.

The PyClical demos
------------------

The subdirectory `./pyclical/demos` contains demo and tutorial scripts written
in Python. These scripts rely on the built PyClical extension module.

The PyClical plotting demos
---------------------------

Two of the demos in `./pyclical/demos/` (`plotting_demo_mayavi.py` and
`plotting_demo_dialog.py`) use Mayavi2 and VTK for 3D visualization.
Because VTK and Mayavi are highly sensitive to library version mismatches,
they require a special execution environment (see "Setting Up the Mayavi
Plotting Environment" below).

The two types of test programs
------------------------------

GluCat includes two types of test programs, timing (benchmark) tests, and
functionality (regression) tests. The regression tests are also used as
programming examples. The source code for the timing tests is in `./gfft_test`,
`./products`, `./squaring` and `./transforms`. The source code for the regression
tests is in the subdirectories `./test00` to `./test19`. Compilation of these test
programs uses C++ headers defined in `./test` as well as `./glucat`.

Code Coverage and Unit Tests
----------------------------

GluCat now includes a modern unit testing suite based on `doctest`, and
infrastructure for generating code coverage reports.

The `./test_doctest` directory contains the build system for the doctest-based
unit tests. These tests are integrated directly into the header files and can
be enabled by defining `GLUCAT_DOCTEST`.

The `./test_coverage` directory contains shell scripts for generating code
coverage reports using Clang/LLVM. See "Running Coverage Tests" below.

Your own programs, written in C++
---------------------------------

Once you have installed the GluCat C++ library, you use this library to build
your own C++ programs. See "Building your own C++ programs" in the `To Build`
section below for the required compiler flags, include paths, and linker
flags.


To Configure
============

The build process for PyClical and the test programs uses a single `./configure`
shell script and a single overall Makefile, plus Makefiles in the relevant
subdirectories.

As briefly described above, the simplest way to install this package is:

 1. `cd` to the `glucat-0.98a2` directory containing the source code and type
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

The preprocessor symbol `NDEBUG` is used to control debugging in the C++
Standard Library and the linear algebra backends (Eigen or Armadillo).

If `NDEBUG` is defined, the library compiles in release mode, which typically
disables runtime assertions and enables optimizations. If `NDEBUG` is not
defined, the library compiles in debug mode, providing more extensive runtime
checks at the cost of performance.


If you are compiling your own programs using the GluCat library, to control
debugging and performance, your Makefile needs to pass the following flags to
the C++ compiler.

 For full production and timing tests:
  `-DNDEBUG -O3`

 For testing with some debugging capability:
  `-DNDEBUG -O1 -g`

 For full debugging capability:
  `-O0 -g3`

```
  --disable-debug         disables debug output and debug symbols [default=yes]
```
Since the default for `--enable-debug` is `no`, the option `--disable-debug`
does nothing.

```
  --enable-strict         compile with strict compiler options (may not work!)
```
This option adds strict compiler flags to `CXXFLAGS`, such as `-pedantic`. Use this
with caution as it may cause the build to fail on warnings.

```
  --disable-warnings      disable compilation with -Wall and similar
```
By default, the compiler flag `-Wall` is added to `CXXFLAGS`. This option prevents
`-Wall` from being added.

```
  --enable-profile        create profiling infos [default=no]
```
This option adds profiling flags (e.g. `-pg`) to `CXXFLAGS` and `LDFLAGS` to
enable performance profiling with tools like `gprof`.

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
  --with-armadillo        use Armadillo library [default=no]
```
This option controls the use of the Armadillo C++ linear algebra library. The only valid values are `yes` or `no`.

The value `yes` (or `--with-armadillo` without an explicit value) adds `-D_GLUCAT_USE_ARMADILLO` to `CXXFLAGS` and the flag `-larmadillo` to the list of libraries, `LIBS` in the Makefiles.

To compile your own programs using GluCat headers with Armadillo, your Makefile needs to pass the flags `-D_GLUCAT_USE_ARMADILLO` and `-larmadillo` (along with `-lflexiblas` or `-lopenblas` if those backends are used) to the C++ compiler. You will also need to ensure that the include path used by the compiler sees `<armadillo>` and the library path sees `libarmadillo.*`.

```
  --with-blas=yes|no|flexiblas|openblas
                          use specified BLAS/LAPACK library backend [default=no]
```
This option controls the BLAS/LAPACK library backend used in conjunction with Armadillo. This option requires that the Armadillo library is also used.

The value `yes` automatically probes the host system's optimized backends:
1. It checks for the `flexiblas` library. If found, it enables FlexiBLAS, which appends `-lflexiblas` to the libraries.
2. If `flexiblas` is not found, it checks for the `openblas` library. If found, it enables OpenBLAS, which appends `-lopenblas` to the libraries.
3. If both are absent, it defaults to standard dynamic linkage with `-larmadillo`.

Specifying `flexiblas` or `openblas` directly skips the auto-probing and forces the configuration to use the respective library.

To compile your own programs using the GluCat library with OpenBLAS or FlexiBLAS, your Makefile needs to pass the flag `-lopenblas` or `-lflexiblas` respectively to the C++ compiler.

```
  --with-openmp           use OpenMP (requires Armadillo) [default=no]
```
This option controls the use of OpenMP for parallel processing. This option
requires that the Armadillo library is also used.

The option `--with-openmp` adds `-D_GLUCAT_USE_OPENMP` and OpenMP compiler flags
(e.g. `-fopenmp` or `-qopenmp`) to `CXXFLAGS` in the Makefiles.


To compile your own programs using the GluCat library with OpenMP, your Makefile
needs to pass the flags `-D_GLUCAT_USE_OPENMP` and the OpenMP compiler flags
to the C++ compiler.

Note: If you are using the clang++ compiler you will need to ensure that libomp
is installed.

```
  --with-doctest          enable the doctest unit testing suite [default=no]
```
This option enables the doctest unit testing suite. The option `--with-doctest`
defines `GLUCAT_DOCTEST` in the header files and enables compiling the doctest-based
unit tests in `test_doctest/`.

```
  --with-unordered-map=boost|std
                          use specified map implementation [default=boost]
```
This option chooses the implementation of the hash map used in `framed_multi`.
The default is `boost`, which uses `boost::unordered_flat_map` (requires
`Boost 1.83.0` or later). The value `std` uses `std::unordered_map`.
Use `boost` for potentially better performance.


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

Building the PyClical extension module
--------------------------------------

To build the PyClical extension module, set the environment variable CXX
to indicate your C++ compiler (e.g. `g++` for GNU C++, `icpx` for Intel C++),
run `./configure`, and then run `make`.

Make uses the headers in `./glucat` and the source in `./pyclical` to build
the PyClical extension module.

The following build steps will be performed if you have not selected the
configuration option `--disable-pyclical` (equivalently, `--enable-pyclical=no`).

If Cython is installed then `make` builds PyClical by running the command:

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

You can run `pyclical/setup.py` yourself, but you must set the environment
variables to appropriate values. See `To Configure` above to determine these
values.

Alternatively, if you have Python installed but do not have Cython, then
`./configure` will recognize this, and make will build PyClical via the command:

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


Building and running the PyClical demos
---------------------------------------

In the `pyclical/demos` directory, there is a Python script
`build_pyclical_notebooks.py`. This script builds a set of Jupyter notebooks
that correspond to the PyClical tutorials and some of the demos. Running `make`
will use this Python script to build these notebooks.

Once built, you can run standard PyClical demos in Python directly:

```bash
cd pyclical/demos
python3 clifford_demo.py
```

To launch the Jupyter notebooks for the PyClical tutorials and demos:

```bash
jupyter notebook pyclical/demos/
```

If you are using the system Python (not a Conda or venv environment), you may want to
register a dedicated Jupyter kernel so that the notebooks use the correct PyClical build:

```bash
make -C pyclical install-pyclical-kernel
```


Building the PyClical plotting demos
------------------------------------

The plotting demos (`plotting_demo_mayavi.py` and `plotting_demo_dialog.py`)
require a specially configured environment to satisfy their Mayavi/VTK dependencies.
Building these demos requires creating/activating the appropriate environment (via Conda
on x86-64 or a venv on ARM), configuring/building PyClical against the correct Python interpreter, and exporting the required runtime plotting environment variables.

Please refer to the detailed section "Setting Up the Mayavi Plotting Environment"
below for the full build and execution procedures for these demos.


Building and running the test programs
--------------------------------------

To build and run the timing and functionality (regression) test programs, set the
environment variable CXX to indicate your C++ compiler (e.g. `g++` for GNU C++,
`icpx` for Intel C++), run `./configure` as above, and then run `make` (for timing tests)
or `make check` (for regression tests).

Make uses the C++ headers in `./glucat` and `./test` to compile the timing test programs
`./gfft_test/gfft_test`, `./products/products`, `./squaring/squaring`, and
`./transforms/transforms` using the source code in those respective directories.

Running `make check` builds and runs the functionality regression tests `./test_move/test_move`
and `./test00/test00` to `./test19/test19`. This produces the intermediate output files
`./test_move/test_move.out` and `./test00/test00.out` to `./test19/test19.out`,
and the final test output file `./test_runtime/test.out`.

Note the difference between `make check` and `make check-local`:
* `make check` is a standard recursive Automake target. It first enters each
  subdirectory in the `SUBDIRS` list (including timing tests like `gfft_test`)
  and runs `make check` there, before running the legacy regression tests.
* `make check-local` is a faster alternative that skips the recursive pass
  through `SUBDIRS` and only runs the legacy functionality tests.

When you want to generate `test_runtime/test.out` to compare with the existing
sample files in `test_runtime.*`, you should use parallel make with the
`check-local` target:

```bash
make -j${NPROCS} check-local
```
where `${NPROCS}` is the number of processes you want to use (e.g. `make -j4 check-local`).

Warning: If you use too many jobs with parallel make, the compiler will have
problems obtaining enough memory to run efficiently.


Code Coverage and Unit Tests
----------------------------

To build and run the modern unit tests, you must configure with the
`--with-doctest` option. Then run:

```
 make -C test_doctest check
```

This will compile and run the unit tests integrated into the library headers.
The output will be written to `test_doctest/test_doctest.out`.

If you have Clang and LLVM installed, you can generate code coverage reports
using the scripts in `test_coverage/src/`.

To run the coverage tests for the doctest suite:

```bash
bash test_coverage/src/run_clang_doctest_coverage.sh [--backend=eigen|arma|both]
```

By default, the script tests both backends and combines their results into a
unified report. You can use the `--backend` argument to select a specific
backend.

The resulting HTML report will be available in
`test_coverage/results/coverage_html_doctest/index.html`.


Building your own C++ programs
------------------------------

To compile C++ applications against the installed GluCat library:

**Include paths**: Add the GluCat installation prefix to your compiler's
include search path, e.g. `-I/usr/local/include`.

**Linker flags**: Link against the same linear algebra backend you configured
GluCat with:
- Eigen (default): no extra `-l` flag needed; headers-only.
- Armadillo: add `-larmadillo` (and `-lflexiblas` or `-lopenblas` if used).
- QD: add `-lqd`.
- OpenMP: add the appropriate OpenMP flag (e.g. `-fopenmp` for GCC/Clang).

**Preprocessor macros**: Pass the same `-D` flags you used when building
GluCat, or omit them for defaults. Key macros are:
- `_GLUCAT_USE_ARMADILLO` — use Armadillo instead of Eigen.
- `_GLUCAT_USE_QD` — enable QD high-precision types.
- `_GLUCAT_USE_OPENMP` — enable OpenMP parallelism (requires Armadillo).
- No macros set — defaults to the Eigen backend and the C++ 2023 standard
  random number generator.

**Build system**: Use GNU make with your own Makefile. For larger projects,
consider GNU Autotools.


Setting Up the Mayavi Plotting Environment
==========================================

The PyClical plotting demos (`plotting_demo_mayavi.py` and
`plotting_demo_dialog.py`) use Mayavi2 and VTK for 3D visualization.
Because VTK and Mayavi are highly sensitive to library version mismatches,
they require a special execution environment.

The setup procedure depends on your hardware architecture because Conda's
`linux-aarch64` VTK/Mayavi binaries are compiled with 4 KB memory page
alignment, which is incompatible with the 16 KB page size required by
Apple Silicon (Asahi Linux). On ARM aarch64, use the system RPM packages
(which are rebuilt for 16 KB alignment by the Asahi team) via a Python
virtual environment instead.

Confirmed working versions:

| Machine | Architecture | OS | Python | Mayavi | VTK |
|:---|:---|:---|:---|:---|:---|
| Tempesta (AMD Ryzen) | x86-64 | Kubuntu 26.04 | 3.12.13 (Conda) | 4.8.3 | 9.4.2 |
| Pensieri (Intel Core) | x86-64 | Kubuntu 25.04 | 3.12.x (Conda) | 4.8.x | 9.4.x |
| Ginestra (Apple M2 Pro) | aarch64 | Fedora Asahi Remix 43 | 3.14.x (system) | 4.8.x | 9.4.x |


When a special environment is not needed
----------------------------------------

- **C++ only** (GluCat headers, test programs `test00`–`test19`, doctest,
  benchmarks): No Python environment is required at all. Configure with
  `--disable-pyclical` to skip PyClical entirely:
  ```bash
  ./configure --disable-pyclical
  make check-local
  ```

- **PyClical without plotting demos** (doctests, tutorials, all demos except
  the Mayavi-based ones): The system Python with NumPy and Cython is
  sufficient. No Conda or venv is required.


x86-64 (Ubuntu, Kubuntu) — Conda path
-------------------------------------

Install [Miniforge](https://github.com/conda-forge/miniforge) or
[Anaconda](https://www.anaconda.com/download) if you do not already have
Conda or Mamba available.

1.  Set up the Conda environment. Run this command from the repository root:
    ```bash
    source pyclical/demos/plotting/setup-plotting-env.sh
    ```
    This script creates or updates the `pyclical-plotting` Conda environment from
    `pyclical/demos/plotting/plotting-env.yml`, activates it, and removes the
    conda-forge `mesalib` if a native GPU is detected (via `/dev/dri/card0`).

2.  Bootstrap the build system (git clone only, skip for release tarballs):
    ```bash
    make -f admin/Makefile.common bootstrap
    ```

3.  Configure and build PyClical (ensure the environment is active first):
    ```bash
    ./configure
    make -C pyclical -j$(($(nproc)/2))
    ```

4.  Export the runtime plotting environment variables:
    ```bash
    source pyclical/demos/plotting/export-plotting-vars.sh
    ```


openSUSE Tumbleweed — RPM path
------------------------------

The use of Mayavi2 on openSUSE Tumbleweed involves installing the
following RPM packages (or their current equivalents):

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

After installing the RPM packages, follow the same venv and build steps as
the ARM aarch64 section below, using `zypper` in place of `dnf`
where appropriate. The `export-plotting-vars.sh` script and the "Running the
Mayavi demos" steps below apply to openSUSE as well.


ARM aarch64 (Fedora Asahi Remix) — system venv path
---------------------------------------------------

Conda's `linux-aarch64` VTK/Mayavi binaries are 4 KB page-aligned and
will segfault immediately on import on Asahi Linux, which requires 16 KB
page alignment. Use the system package manager to install VTK and Mayavi,
then create a virtual environment that inherits those packages.

1.  Install the 16 KB-aligned VTK/Mayavi packages from Fedora Asahi RPMs:
    ```bash
    sudo dnf install python3-mayavi python3-vtk python3-qt5
    ```

2.  Create a virtual environment that inherits those system packages.
    Run from the glucat repository root (`.venvs/` is listed in `.gitignore`):
    ```bash
    python3 -m venv --system-site-packages .venvs/pyclical-mayavi
    source .venvs/pyclical-mayavi/bin/activate
    ```

3.  Configure and build PyClical against the venv's Python interpreter.
    The venv must be active before running `./configure`, because
    `$(which python3)` resolves to whichever interpreter is on `$PATH` at the
    time `./configure` runs. If the venv is not active, PyClical will be built
    against the system Python ABI and will not have access to the venv's packages.

    Run from the repository root. If working from a git clone (not a release
    tarball), run the Autotools bootstrap first:
    ```bash
    make -f admin/Makefile.common bootstrap
    ```

    Then configure and build:
    ```bash
    ./configure PYTHON=$(which python3)
    make -C pyclical -j$(($(nproc)/2))
    ```

After completing these steps, proceed to "Running the Mayavi demos (all platforms)" below, beginning with:
```bash
source pyclical/demos/plotting/export-plotting-vars.sh
```


Running the Mayavi demos (all platforms)
----------------------------------------

Once the environment is set up and PyClical is built, the remaining steps
are identical on all platforms.

1.  From the repository root, source the runtime environment script.
    This sets `PYTHONPATH`, `QT_API`, and `QT_QPA_PLATFORM`:
    ```bash
    source pyclical/demos/plotting/export-plotting-vars.sh
    ```

2.  Run either plotting demo from the `pyclical/demos` directory:
    ```bash
    cd pyclical/demos
    python3 plotting_demo_mayavi.py
    python3 plotting_demo_dialog.py
    ```

3.  To run in non-interactive mode (e.g. for automated testing):
    ```bash
    cd pyclical/demos
    GLUCAT_NON_INTERACTIVE=1 python3 plotting_demo_mayavi.py
    ```




To Test
=======

The test_runtime.* directories
------------------------------

The static test runtime baseline directories `./test_runtime.x86-64` and
`./test_runtime.aarch64` contain sample test output files.

The sample test output files include `expressions-8.out`, `gfft_test-11.out`,
`products-8.out`, `squaring-11.out`, and `transforms-8.out`.

The `./test_runtime.x86-64` directory contains 24 sample versions of the regression
test results (corresponding to 12 different combinations of configuration
parameters, for two different sets of tests: the complete set of 20 tests, and a
subset of 3 tests), and also includes the sample output file `eg3.res` and the
test input file `eg8.txt` (needed by programming example 8, reading multivectors from input).

The `./test_runtime.aarch64` directory contains 12 sample versions of the full
regression test results (corresponding to 12 different combinations of configuration
parameters for the complete set of 20 tests).


Re-running the regression tests
-------------------------------

The main regression test script is `./test/test.sh`. Once you have built and run
the regression tests via `make check`, you can use this script to re-run the
regression tests. The script `./test/test.sh` reruns tests `./test00/test00` to
`./test19/test19`, using relative pathnames, so it is best to leave `test.sh`
where it is and invoke it using its full path name. This allows it to find
`test00` to `test19`.

The test script `./test/test.sh` takes any number (including zero) of numeric
parameters. Parameters in the range `00` to `19` correspond to coding examples
`./test00/test00` to `./test19/test19`. These examples are run in numerical
order. With zero parameters, all examples from `00` to `19` are run in order.
Many of the examples are run twice - once with `framed_multi<Scalar_T>` and once
with `matrix_multi<Scalar_T>`.

The baseline test results are structured as follows:

The tests results in `./test_runtime.aarch64` come from runs on an Apple M2 Pro
with 6 Avalanche performance cores and 4 Icestorm efficiency cores with
    ```
    Linux 6.19.14-400.asahi.fc43.aarch64+16k SMP aarch64
    Fedora Linux Asahi Remix 43 (KDE Plasma Desktop Edition)
    g++ 15.2.1 20260123 (Red Hat 15.2.1-7)
    Armadillo 12.8.1
    Doxygen 1.14.0
    Eigen3 3.4.0
    GSL 2.8
    QD 2.3.24
    Cython 3.1.3
    Numpy 2.3.5
    Python 3.14.5
    ```
The tests results in `./test_runtime.x86-64` come from runs on an 8 core
AMD Ryzen 7 8840HS w/ Radeon 780M Graphics` @ 3.3 GHz with
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

  7. `test.configure.armadillo.out`:

```
./configure --with-armadillo
```
  7. `test.configure.armadillo-debug-yes.out`:

```
./configure --with-armadillo --enable-debug=yes
```
  9. `test.configure.armadillo-openmp.out`:

```
./configure --with-armadillo --with-openmp
```
 10. `test.configure.armadillo-qd-std.out`:

```
./configure --with-armadillo --with-qd --with-unordered-map=std
```
 11. `test.configure.qd.out`:

```
./configure --with-qd
```
 12. `test.configure.std.out`:

```
./configure --with-unordered-map=std
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

The tests typically use floating point arithmetic, and `./test00/test00` and
`./test11/test11` in particular also use random number generators. Therefore if
you run the tests using different architecture, compilers or random number
generators, you should expect to have different floating point arithmetic
results, but generally, still within acceptable error tolerances.

The regression tests `./test00/test00` to `./test19/test19` recognize the program
arguments `--help`, `--no-catch`, and `--verbose`. The `--no-catch` argument
disables the default exception catching behaviour of a regression test, to
allow program crashes to be more easily debugged. For `./test00/test00` and
`./test11/test11` the argument `--verbose` produces verbose output essentially by
setting the error tolerance to zero. Verbose output can become quite large.

The test script `./test/test_optional.sh` runs all examples 00 to 19 in order,
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
or `./test_runtime/fast-test.configure.$abbreviation.out` respectively where the
directory `./test_runtime` is created if it does not already exist.

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
You can also give parameters to `./test/test-all-config-options.sh` and these are
passed to the make command. In particular, invoking (e.g.)

```
./test/test-all-config-options.sh -j 4
```
performs a parallel `make check` for each configuration option, potentially
speeding up the entire testing process.

Rather than running the regression tests in-place and copying the output
directly into `./test_runtime`, the script `./test/test-all-config-options.sh`
produces as many copies of the whole directory `glucat-0.98a2` as there are lines
in `./test/config-options.txt`, naming them `glucat-0.98a2.1` to `glucat-0.98a2.12`,
in the parent directory of `glucat-0.98a2`. This allows the effect of each set
of configuration options to be directly compared, and also ensures that any
side-effect of a configuration does not affect the test results of another
configuration.

The script `./test/diff-all-config-outputs.sh ${arch}` compares each relevant test
output file with the corresponding file in `./test_runtime.${arch}` or `./pyclical`.
For example, line 4 of `./test/config-options.txt`

```
disable-dependency:          --disable-dependency-tracking
```
causes `./test/diff-all-config-outputs.sh` to use diff to compare
`glucat-0.98a2.4/test_runtime/test.configure.disable-dependency.out` to
`glucat-0.98a2/test_runtime.${arch}/test.configure.disable-dependency.out`, and
compare `glucat-0.98a2.4/pyclical/test.out` to `glucat-0.98a2/pyclical/test.out`.

Each comparison should only produce a line containing the line number of
the configuration being compared: 1 to 12.

The exceptional cases are:

1. If the configuration options have caused any sort of error.
2. Differences in compilers and libraries causing different floating point
   results. A compiler difference previously occurred with the Intel C++ compiler,
   which produced output different from either the GNU C++ compiler or the Clang
   compiler. As an example of a library difference, QD version 2.3.16 contains an
   update that fixes a problem with `tanh`.

If the output of your systematic tests differs due to a difference in compilers
or libraries, you may want to copy this output to the staging `./test_runtime`
directory to ease future comparisons without clobbering the baseline
`./test_runtime.${arch}` directory directly. To do so, run the script
`./test/copy-all-config-outputs.sh`.

The difference between `./test/test-all-config-options.sh` and
`./test/fast-test-all-config-options.sh` is that the former runs all of the tests
`./test00` to `./test19`, whereas the latter runs only `./test00`, `./test10` and
`./test11`. Both scripts build and check `./pyclical`. More specifically, the
`./test/test-one-config-option.sh` script runs `make check`, and the
`./test/fast-test-one-config-option.sh` script runs `make fast-check`. For the
definitions of these arguments to `make`, see the file `./Makefile.am.in`. To
compare or copy the output of `./test/fast-test-all-config-options.sh`, use either
`./test/fast-diff-all-config-outputs.sh ${arch}` or `./test/fast-copy-all-config-outputs.sh`
respectively.

The script `./test/pyclical-test-all-config-options.sh` just builds and checks
`./pyclical` without running any of the other regression tests. To examine the
output of `./test/pyclical-test-all-config-options.sh` just examine all of the
`pyclical/pyclical-test.check.out` files for errors.


Running the timing (benchmark) tests
------------------------------------

All timing test outputs report execution times in milliseconds (ms).

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

The test script `./test/timing_tests.sh` takes up to 5 numeric parameters.
The command `./test/timing_tests.sh $a $b $c $d $e` runs

```
 ./products/products $a
 ./squaring/squaring $b
 ./gfft_test/gfft_test $c
 ./transforms/transforms $d
 ./expressions/expressions $e
```
The default is:

```
 ./products/products 8
 ./squaring/squaring 11
 ./gfft_test/gfft_test 11
 ./transforms/transforms 8
 ./expressions/expressions 8
```

The sample timing test results are organized in target-specific baseline directories:
* `./test_runtime.x86-64` (for AMD/Intel x86-64 platforms)
* `./test_runtime.aarch64` (for ARM64/Apple Silicon platforms)

The sample outputs were generated using:

```
./configure
```
on an 8 core `AMD Ryzen 7 8840HS w/ Radeon 780M Graphics` @ 3.3 GHz with
```
    Linux 6.17.0-12-generic #12-Ubuntu SMP UTC
    Kubuntu 25.10
    g++ 15.2.0 (Ubuntu 15.2.0-4ubuntu4)
    Boost 1.88.0
    Eigen 3.4.0
    GSL 2.8
    QD 2.3.23
```

Systematic Benchmarking
-----------------------

In addition to individual timing tests, the test suite includes scripts for systematic performance benchmarking across 16 different library configurations (varying backends such as Eigen/Armadillo, parallel and serial BLAS, and OpenMP settings). Note that these benchmarking directories and scripts are included in the Git repository but are deliberately excluded from the installation tarball to prevent bloat.

These benchmarks are driven by the following scripts in the `./test` directory:
* `./test/benchmark-all-config-options.sh`: Builds and runs benchmarks for all 16 configurations specified in `./test/benchmark-config-options.txt`.
* `./test/benchmark-one-config-option.sh`: Runs a single configuration by line number.
* `./test/copy-all-benchmark-outputs.sh` / `./test/copy-one-benchmark-output.sh`: Copies benchmark results to the `./benchmarks` directory.

To ensure consistent processor affinity, multithreading, and cache behavior, the benchmarks source configuration-specific environment scripts from the `./benchmarks` directory, which route environment variables (such as `OMP_NUM_THREADS` and `OPENBLAS_NUM_THREADS`) depending on the profile:
* `env-*.sh`: A set of 16 wrappers corresponding to each configuration abbreviation.
* `env_setup_common.sh`: A common script performing platform-independent CPU architecture probing (isolating Performance cores on hybrid Apple Silicon Asahi Linux platforms, and targeting physical cores on homogeneous AMD/Intel x86-64 Linux).

Results and comparative reports from these systematic benchmarks are stored in the `./doc/benchmarks/` directory. In particular, the file `./doc/benchmarks/compiler_architecture_comparison_report.md` contains a comprehensive performance analysis across three hardware architectures (Intel Core i7-870, AMD Ryzen 7 8840HS, Apple Avalanche M2 Pro) and three compilers (GCC, Clang, Intel oneAPI), accompanied by performance scaling plots.

Testing PyClical
----------------
Once you have built PyClical, run the doctests. In `python3` or `ipython3`, etc.:

```
 >>> import PyClical
 >>> PyClical._test()
 TestResults(failed=0, attempted=661)
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

GluCat 0.98a2 with PyClical has so far been built and tested using:

 1) Tempesta:
    8 core `AMD Ryzen 7 8840HS w/ Radeon 780M Graphics` @ 3.3 GHz with

    ```
    Linux 7.0.0-15-generic #15-Ubuntu SMP 2026
    Kubuntu 26.04 LTS
    Armadillo 15.2.1
    Doxygen 1.15.0
    Eigen3 3.4.0
    GSL 2.8
    QD 2.3.23
    TeXLive 2025.20260124

    Conda 26.3.2 env containing:
    Boost 1.88.0
    Cython 3.1.6
    Jupyter Server 2.18.2
    Matplotlib 3.10.9
    Mayavi2 4.8.3
    Notebook 7.5.6
    Numpy 1.26.4
    PyQT 5.15.11
    Python 3.12.13
    VTK 9.4.2
    ```

    `./test/test-all-config-options.sh`:
    All 12 configuration commands corresponding to each of the 12
    `test.configure*.out` files in `./test_runtime`
    tested with the following compiler versions:

    1) `g++ (Ubuntu 15.2.0-16ubuntu1) 15.2.0`
    2) `Ubuntu clang version 21.1.8 (6ubuntu1)`

 2) Pensieri:
    4 core `Intel(R) Core(TM) i7 CPU 870  @ 2.93GHz` with

    ```
    Linux 6.14.0-35-generic #35-Ubuntu SMP 2025
    Kubuntu 25.04 LTS
    Armadillo 15.2.1
    Boost 1.83.0
    Cython 3.0.11
    Doxygen 1.9.8
    Eigen3 3.4.0
    GSL 2.8
    Numpy 2.2.3
    QD 2.3.23
    Python 3.13.3
    TeXLive 2024.20250309

    ```
    `./test/fast-test-all-config-options.sh`:
    All 12 configuration commands corresponding to each of the 12
    `fast-test.configure*.out` files in `./test_runtime`
    tested with the following compiler versions:

    1) `g++ 14.2.0 (Ubuntu 14.2.0-19ubuntu2)`
    2) `Ubuntu clang version 20.1.2 (0ubuntu1)`
    3) `Intel(R) oneAPI DPC++/C++ Compiler 2025.0.4 (2025.0.4.20241205)`

 3) Vincitor (Pensieri running VirtualBox):
    Virtual 1 core `Intel(R) Core(TM) i7 CPU 870 @ 2.93GHz` with

    ```
    Linux 6.12.6-1-default #1 SMP 2024
    openSUSE Tumbleweed Release 20241224
    g++ (SUSE Linux) 14.2.1 20241007
    Armadillo 15.2.4
    Boost 1.86.0
    Cython 3.0.12
    Doxygen 1.12.0
    Eigen 5.0.1
    GSL 2.8
    Numpy 2.1.3
    Python 3.11.11
    QD 2.3.24
    TeX Live 2026.20260301
    ```
    `./test/fast-test-all-config-options.sh`
    All 12 configuration commands corresponding to each of the 12
    `fast-test.configure*.out` files in `./test_runtime.x86-64`

 4) Ginestra
    Apple M2 Pro with 6 Avalanche performance cores and 4 Icestorm efficiency cores with
    ```
    Linux 6.19.14-400.asahi.fc43.aarch64+16k SMP aarch64
    Fedora Linux Asahi Remix 43 (KDE Plasma Desktop Edition)
    g++ 15.2.1 20260123 (Red Hat 15.2.1-7)
    Armadillo 12.8.1
    Cython 3.1.3
    Doxygen 1.14.0
    Eigen3 3.4.0
    GSL 2.8
    Numpy 2.3.5
    Python 3.14.5
    QD 2.3.24
    ```
    `./test/test-all-config-options.sh`
    All 12 configuration commands corresponding to each of the 12
    `test.configure*.out` files in `./test_runtime.aarch64`


Notes on software versions
==========================

This section documents currently encountered building bugs and workarounds:

*   **LaTeX/Doxygen version constraints**:
    Building the documentation requires recent versions of both `doxygen` and
    `latex` (e.g. `texlive-2024`). Older versions may cause build errors.

*   **GCC 15 false-positive warnings**:
    Under GCC 15 on Fedora Asahi Remix, false-positive `-Wmaybe-uninitialized`
    warnings can be encountered during optimized compilation of Eigen/Armadillo headers.
    To compile with strict options enabled, either configure with `--disable-strict` or
    suppress warnings in your compiler flags.
