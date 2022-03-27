README for GluCat 0.11.0 with PyClical
======================================

GluCat is a library of C++ template classes for calculations with the universal
Clifford algebras over the real field. The dimension of the algebras you can
use is limited only by computer word size. Unlike some software packages for
Clifford algebras, the maximum signature is fixed by default or by user defined
template parameter values, and calculations are performed in the appropriate
subalgebra.

GluCat classes are meant to be used as numeric classes with other template
libraries. To do this, GluCat classes need to look and behave like floating
point  numeric types such as float and double. The GluCat classes do this by
providing a definition for std::numeric_limits<>, most or all the expected
operators like +, -, *, /, and functions like sqrt(). In use with other template
libraries, you will need to look out for the differences between GluCat classes
and floating point numeric types. Two key differences are that Clifford algebras
do not have a total ordering, and that multiplication in Clifford algebras is
not necessarily commutative. A third difference is that not all non-zero
elements of a Clifford algebra have a multiplicative inverse.

PyClical is a Python extension module that is designed to make it easy for you
to use Python to perform computations in Clifford algebras, and to implement
algorithms and programs that use these computations. The syntax and semantics
of the main index_set and clifford classes in PyClical is designed to be as
similar as possible to that of the GluCat template classes index_set<> and
matrix_multi<double>, so that you should be able to use PyClical to teach
yourself how to use calculations in Clifford algebras, and then transfer that
knowledge to your ability to use GluCat template classes in C++.


Before you go any further
=========================

To make the best use of the GluCat template library and the PyClical extension
module, you will need to become familiar with Geometric Algebra and Clifford
algebras in general. The AUTHORS file includes lists of recommended web pages,
books, and software that can help you to learn about Clifford algebras. Here is
a short list to get you started.

* Wikipedia: [Geometric algebra
  ](http://en.wikipedia.org/wiki/Geometric_algebra)

* Alan Macdonald, [A Survey of Geometric Algebra and Geometric Calculus
  ](http://faculty.luther.edu/~macdonal/GA&GC.pdf)

* Alan Macdonald, [Linear and Geometric Algebra, 2011.
  ](http://books.google.com.au/books?q=isbn:9781453854938)

You should also install PyClical and run the tutorials. Once you have run the
PyClical tutorials, run the PyClical demos and look at the sample demo output.
Now examine the C++ code for the GluCat tests, as well as the sample output.
Finally, look at Glucat API documentation. Instructions for each of these steps
are given below.


Getting started
===============

To use the software to its full potential, you will want to:
1. Download and update the software,
2. Understand the source code directory structure,
3. Resolve dependencies, install and test the software,
4. Use the PyClical extension module with Python,
5. Write C++ programs that use GluCat classes,
6. Compile and run C++ programs that use GluCat classes, and finally,
7. Help to improve the software by requesting features, filing bug reports and
   writing your own source code patches.

These tasks are treated in some detail below.


1 Downloading and updating the software
---------------------------------------

You can download new versions of GluCat from
http://sourceforge.net/projects/glucat/files/
or alternatively, use the Git version:
From SourceForge: http://sourceforge.net/p/glucat/git/ci/master/tree/
From GitHub: https://github.com/penguian/glucat
For more on installing from Git, see `INSTALL.md`.


2 Understanding the source code directory structure
---------------------------------------------------

GluCat is a C++ template library, similar to the C++ Standard Library or the
Boost uBLAS Library (uBLAS). It consists of source code header files, a
suite of test routines, and the PyClical Python extension module and associated
files.

Once you have downloaded, unzipped and untarred the GluCat source code,
you should have a directory, `glucat-xxx`, where `xxx` is the version number.
Under `glucat-xxx` you should see a number of directories, including `./admin`,
`./doc`, `./gfft_test`, `./glucat`, `./m4`, `./products`, `./pyclical`,
`./squaring`, `./test`, `./test_runtime`, `./testxx`,  and `./transforms`.

The `./glucat` directory contains all the header files that define the GluCat C++
template library.

The `./pyclical` directory contains the C++ and Cython source code for the
PyClical Python extension module, as well as a subdirectory `./pyclical/demos`
containing Python source code for the PyClical demos and tutorials, and sample
demo output.

The `./admin` and `./m4` directories are part of the `autotools` infrastructure for
building GluCat, and should normally be left unchanged.

The `./doc` directory contains documentation. Currently only the GluCat API
Reference Manual can be found here, under `./doc/api`.

The `./gfft_test`, `./products`, `./squaring` and `./transform` directories
contain the C++ source code for timing tests for GluCat.

The `./test` and `./testxx` directories contain the C++ source code for
programming examples and regression tests for GluCat.

The `./test_runtime` directory contains regression test input and sample output
for the GluCat timing and regression tests.


3 Resolving dependencies, installing and testing the software
-------------------------------------------------------------

Detailed instructions for these tasks are included in the `./INSTALL.md` file in the
same directory that contains this `README.md` file.


4 Using the PyClical extension module with Python
-------------------------------------------------

The PyClical Python extension module is written in C++ and Cython, and is defined
in the files `pyclical/glucat.pxd`, `pyclical/PyClical.h`, `pyclical/PyClical.pxd`,
and `pyclical/PyClical.pyx`. PyClical is designed to be installed using `make`. For
details on building PyClical, see the `./INSTALL.md` file.

The following instructions assume that you have already installed PyClical. If
you have only run `make` within the PyClical directory, but have not yet
installed PyClical, then, assuming you are using the `bash` interpreter on Linux,
you will need to set the `PYTHONPATH` environment variable so that Python can
find your newly built copy of PyClical. If the `make` has succeeded, you should
have the file `./pyclical/PyClical.so`. Set `PYTHONPATH` to include the full
`./pyclical` directory path name before any other path names. For example:
```
  export PYTHONPATH=/home/leopardi/src/glucat/pyclical:$PYTHONPATH
```
or you can change the `PYTHONPATH` variable for just one command, e.g.
```
  PYTHONPATH=/home/leopardi/src/glucat/pyclical:$PYTHONPATH python3
```
PyClical is designed to be used within a Python environment. You will usually
need to run a Python IDE or interpreter, such as IDLE, `ipython3` or `python3`. The
following instructions use the standard `python3` interpreter.

To use the capabilities of PyClical from within Python, you must either import
the PyClical extension module or import objects from this module.  The simplest
way to do this is to use the following Python statement:
```
>>> from PyClical import *
```
Probably the easiest way to get familiar with PyClical is to make a copy of the
`pyclical/demos` directory and run the tutorials and demos. By default, the
`pyclical/demos` directory installs into `/usr/local/share/pyclical/demos`.

For example, assuming you are using the Bash shell on Linux, and have installed
PyClical, use the following commands:
```
% cp /usr/local/share/pyclical/demos /path/to/my/demos
% cd /path/to/my/demos
% python3 pyclical_tutorials.py
```
where you must replace `/path/to/my/demos` with the real pathname you want to use.

The `pyclical_tutorials` program starts by displaying
```
    Currently available PyClical tutorials:

    0.0 Notation.
    0.1 Index sets.
    0.2 Operations.
    0.3 Algebraic functions.
    0.4 Square root and transcendental functions.
    1.0 Plane geometry.
    1.1 Complex numbers.
    1.2 Space geometry and vector algebra.
    1.3 Electromagnetism and Lorentz transformations.
    1.4 The fourth dimension.
    1.5 Conformal Geometric Algebra.
    2.0 Exterior product.

    Enter the number for the tutorial you want to run, or Q to quit:
```
Just enter one of the numbers, e.g. `0.0`, or `Q` to quit. At any point in an
individual tutorial, you can interrupt by entering `CTRL-c`.

If you are running Linux or a Unix equivalent, the following should also work:
```
% cd /path/to/my/demos
% chmod +x pyclical_tutorials.py
% ./pyclical_tutorials.py
```

For more usage examples, see the example Python files `clifford_demo.py`,
`pyclical_demo.py`, `plotting_demo.py`, `plotting_demo_dialog.py`,
`plotting_demo_mayavi.py`, and `sqrt_log_demo.py`, and the example output files
`pyclical_demo.out` and `sqrt_log_demo.out`.

To run `clifford_demo.py`, `pyclical_demo.py`, or `sqrt_log_demo.py`, use the
following commands:
```
% cd /path/to/my/demos
% python3 $DEMO.py
```
where `$DEMO` is one of `clifford_demo`, `pyclical_demo` or `sqrt_log_demo`.

If you are running Linux or a Unix equivalent, the following should also work:
```
% cd /path/to/my/demos
% chmod +x $DEMO.py
% ./$DEMO.py
```
To run `plotting_demo.py`, run `ipython3 --pylab` and enter the following command
at the IPython prompt:
```
In [1]: %run plotting_demo
```
This demo uses Matplotlib to produce a number of plots.

To run `plotting_demo_mayavi.py`, first ensure that you have Mayavi2 and wxPython
installed and working.
(See http://code.enthought.com/projects/mayavi/ and http://www.wxpython.org/ )
Then run `ipython3` and enter the following command at the `ipython3` prompt:
```
In [1]: %run plotting_demo_mayavi
```
This demo uses Matplotlib to produce a number of plots. With Mayavi and wxPython,
these plots are displayed in interactive windows, you can rotate, zoom and pan
them. See (e.g.) http://mayavi.sourceforge.net/docs/guide/ch03s04.html

You can also run the Mayavi plotting demo from a graphical user interface.
To do this, run `./plotting_demo_dialog.py`.

The tutorials and demos are also accompanied by a corresponding set of Jupyter
notebooks. To build the notebooks, see `INSTALL.md`. To run the notebooks, assuming
that you have a copy of the notebooks in the directory `/path/to/my/demos`, you have
`ipython3` installed, you have installed PyClical, or set `PYTHONPATH` appropriately,
and you are able to use Jupyter notebooks via a web browser, use the following
commands:
```
% cd /path/to/my/demos
% jupyter notebook
```
Your web browser should open a new window or tab, displaying a page that lists
all of the available tutorials and demos as notebooks. To select a notebook,
click on the corresponding name in the list.


5 Writing C++ programs that use GluCat classes
----------------------------------------------

Once you have familiarized yourself with Clifford algebras and have tried using
PyClical, take a good look at the test C++ code in `./test00` to `./test17` and the
test output in `./test_runtime`.

A good way to begin writing your own C++ code using GluCat is to start with the
programming example code in `./test01`. The file `test01/peg01.cpp` includes
`test/driver.h` and `test01/peg01.h`.

The key lines of code in `test/driver.h` are:
```
 #include "glucat/glucat.h"
 #include "glucat/glucat_imp.h"
 #include "test/tuning.h"
```
The first line includes `"glucat/glucat.h"`, a convenience header that includes
all the GluCat declarations.

The second line includes `"glucat/glucat_imp.h"`, a convenience header that
includes all the GluCat definitions.

The third line includes `"test/tuning.h"`, a convenience header that defines
specific instances of `tuning<>` for use as `Tune_P`, the tuning policy class.
This class is used as a template parameter in in the template classes
`framed_multi<>` and `matrix_multi<>`. The `tuning<>` template class is defined in
`"glucat/tuning.h"`. See this file and `"test/tuning.h"` for examples of values
that you can use for the template parameters of `tuning<>` for use as `Tune_P`.

The key line of code in `test01/peg01.h` is:
```
 using namespace glucat;
```
Glucat defines and uses the ``glucat::` namespace. You can use names from this
namespace by using the `glucat::` prefix or by the `using` declaration above.


To obtain detailed information on the GluCat namespaces, classes and functions,
see the Doxygen documentation in `doc/api/GluCat-API-reference-manual.pdf` (PDF)
and `doc/api/html/` (HTML). By default, this documentation is installed in
the directories `/usr/local/doc/glucat/api/pdf` and `/usr/local/doc/glucat/api/html`
respectively.

5.1 Truncation
--------------

Because the transform from `framed_multi<>` to `matrix_multi<>` may involve many
multiplications and additions, rounding errors can accumulate and result in one
or more non-zero matrix entries where exact arithmetic would have yielded zero.
Similarly, the inverse transform from `matrix_multi<>` to `framed_multi<>` can
result in small non-zero terms which would be zero in exact arithmetic.

The `clifford_algebra::truncated()` member function returns the multivector
`*this` with all relatively small terms set to zero. The function compares the
absolute value of each term to the largest absolute value of any term, and sets
to zero any term whose relative absolute value is smaller than the given limit.
The default limit is `std::numeric_limits<Scalar_T>::epsilon()`. The `truncated()`
matrix function operates similarly to set relatively small matrix entries to zero.

5.1.1 Printing uses truncation
------------------------------

Printing using `operator<<()` uses truncation to suppress the printing of
relatively small terms. Depending on the floating point format used, the
truncation used is calculated using the floating point precision as follows.

  scientific:
  ```
    truncation = 10 ** -(precision+1)
  ```

  fixed:
  ```
    truncation = (10 ** -precision) / max_abs
  ```
  This has the effect of deleting all terms whose absolute value is less than
  `10 ** -precision`.

  hexfloat:
  ```
    truncation = default_truncation
  ```

  default:
  ```
    truncation = 10 ** -precision
  ```

See (e.g.) https://www.cplusplus.com/reference/ios/ios_base/precision/
and see `./test17` for examples.


6 Compiling and running C++ programs that use GluCat classes
------------------------------------------------------------

You can use the `Makefile` for `./test01` as the starting point for your own
`Makefile`. Your `Makefile` needs to pass the appropriate flags to your compiler.

 * Make sure that the header installation directory is in your include path.
   The default header installation directory is `/usr/local/include/`.

 * GluCat uses preprocessor defines to control header files.
   See "To Configure" in the `./INSTALL.md` file for details of these preprocessor
   defines.


7 Helping to improve the software
---------------------------------

To request features, file bug reports or submit source code patches:

The fastest way to ask for help with GluCat or to submit patches is to send an
email to the project manager ( Paul Leopardi <leopardi@users.sourceforge.net> ).

The SourceForge page http://sourceforge.net/projects/glucat/support also
provides project forums, project mailing lists and project trackers.

The GitHub project page at https://github.com/penguian/glucat provides a
convenient interface that allows you to fork the code and create pull requests.

If you are thinking of writing patches, please try to match the programming
style used in the relevant source files, and be aware of the coding standards
used, as listed below.


Coding standards
================

The headers are split into declarations and definitions. The software was
developed using the GNU g++ compiler. Separate compilation is possible, if
you include both the declarations and definitions in each compilation unit,
but compilation is slow and the resulting binary is large.

The code follows many, but not all of the guidelines in the
GNU C++ Standard Library Style Guidelines at
https://gcc.gnu.org/onlinedocs/libstdc++/manual/source_code_style.html

The code also follows much of Scott Meyers' advice from "Effective C++",
Second Edition.

Some code conventions are:
* `name_t`: The name of a type, either global or local to a class.
* `Sentence_Case_T`: The name of a type used as a template parameter.
* `Other_Sentence_Case`: Other template parameters, including template template
parameters.
* `ALL_CAPS_WITH_UNDERSCORES`: A global constant defined in `<glucat/global.h>`
