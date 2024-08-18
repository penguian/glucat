TODO for GluCat 0.12.1 with PyClical
====================================

Documentation:
* Write a programmer's guide with descriptions of usage via use cases.
* Provide better user documentation for PyClical.

Packaging:
* Improve the packaging of the example and test programs.

Portability:
* Build and test GluCat and PyClical using Cygwin on Windows
  (requested by Alan Bromborsky).
* Port to other architectures and compilers which support template template
  parameters.

Interfaces:

*  Downwards:
   * Remove the Boost bindings interface.
   * Try using Blaze, Eigen or Armadillo as a replacement for uBLAS.

*  Upwards:
   * Expand the Cython-based Python extension module PyClical into a Sage interface.
   * Try defining Boost concepts and more numeric traits so that GluCat can
     eventually become a Boost library.

Transcendental functions:
* Devise better algorithms and better implementations of existing algorithms for
  the transcendental functions. In particular, pay more attention to radius of
  convergence, condition number of matrices and poking out of the subalgebra.
* Investigate the use of matrix decompositions in the evaluation of
  transcendental functions. See N. J. Higham, Functions of Matrices: Theory and
  Computation, 2008.

Experimental:
* Expand the use of numeric type promotion and demotion from transcendental
  functions to operations such as division.
* Port to C++17 to take maximum advantage of C++17 semantics and idioms.
* Investigate the use of expression templates.
* Try refactoring the relationship between `matrix_multi`, `framed_multi` and
  `clifford_algebra` to allow more flexibility with template parameters.
  Possibly use enable_if and SFINAE to do this.
* Try adding a `Matrix_Tag` template parameter to `framed_multi` and `matrix_multi`,
  to determine if `matrix_t` is `compressed`, `dense`, etc.
* Try removing the template parameters `LO` and `HI` from `framed_multi` and
  `matrix_multi`, and using `DEFAULT_LO` and `DEFAULT_HI` where these are needed in
  `framed_multi.h`, `matrix_multi.h`, etc.
* Add convenience constructors to `index_set<>`: `index_set<>(int, int)`,
  `index_set<>(int, int, int)`, ... etc.
* Try replacing multiplication by +/-1 within inner products by addition and
  subtraction.
* Try creating a class, `vector_multi<>`, which uses `std::vector` rather than
  `std::map`. This should be faster than `framed_multi<>`, if tuned properly, for
  full multivectors. For sparse multivectors, it may be slower.
