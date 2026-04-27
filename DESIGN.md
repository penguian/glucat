GluCat design notes 2016-07-10, updated 2026-04-25
==================================================

This document describes some of the decisions that underly the design of GluCat,
concentrating mainly on the use of linear algebra.


Early design history of GluCat
------------------------------

Here I describe some of the early history of GluCat, and the reasons why it
initially used MTL and then used uBLAS for linear algebra.

I started GluCat as a coursework Masters project in mathematics at UNSW. Chapter
5 of my coursework Masters thesis contains a description of GluCat in its early
form, as at March 2002, before GluCat used uBLAS. Section 10 lists design
decisions, including the criteria that led to the choice of MTL as the initial
linear algebra library for GluCat. See also:
* [[MTL] "Evaluation of MTL for GluCat" Paul C. Leopardi (2001-12-02 20:12:48) - archived
  ](https://web.archive.org/web/20020121044809/http://www.osl.iu.edu/MailArchives/mtl-devel/msg00244.php),
* [[MTL] "GluCat 0.0.4: Generic library of universal Clifford algebra templates" Paul C. Leopardi (2002-01-24 23:58:04) - archived
  ](https://web.archive.org/web/20020306221057/http://www.osl.iu.edu/MailArchives/mtl-devel/msg00280.php).


Around May 2002, I decided to migrate GluCat from MTL. For my initial questions
about uBLAS, see
* [[boost] "uBLAS vs MTL for sparse matrices" Paul C. Leopardi (2002-05-03 14:38:58) - archived
  ](http://lists.boost.org/Archives/boost/2002/05/29119.php).


In June 2002, Joerg Walter, one of the key developers of uBLAS at the time,
contacted me with a port of GluCat from MTL to uBLAS, based on GluCat 0.0.6.
In the file ublas.h, this contained the functions `lu_factorize` and `lu_solve`, to
replace the use of `mtl::lu_factorize` and `mtl::lu_solve` in `operator/=()` in
`glucat/matrix_multi_imp.h`.

Joerg Walter continued to help the porting effort until and after I released
the first version of GluCat to use uBLAS:
[GluCat 0.1.0 on 30 December 2002.
](https://sourceforge.net/p/glucat/news/2002/12/glucat-010-uses-ublas-and-gcc-32-or-intel-c-70/)


So, in fact, the two main reasons that I moved GluCat from MTL to uBLAS were:

1. The operator syntax and deep copy semantics of uBLAS are closer to Matlab and
mathematical notation than the function calls and shallow copying used by MTL.

2. Joerg Walter did the vast bulk of the porting work, including writing
`lu_factorize` and `lu_solve`, as well as using GluCat as a test case for uBLAS
compressed matrices.

Some other reasons are listed in the SourceForge news item linked above.


Up until GluCat 0.2.3 in 2007, GluCat used compressed matrices as `matrix_t` in
`matrix_multi_t`. The reason why GluCat switched to using dense matrices was that
uBLAS operations on compressed matrices, especially addition, was (and probably
still is) asymptotically slow. For details, see:
* [[ublas] "Worst case time complexity of ublas::compressed_matrix operations?" Paul C. Leopardi (2007-01-21 08:04:56)
  ](http://lists.boost.org/ublas/2007/01/1687.php),
* [[ublas] "matrix addition speed" Per Abrahamsen (2008-01-21 12:55:40)
  ](http://lists.boost.org/ublas/2008/01/2585.php)


Modernization phase (2024-2026)
-------------------------------

Starting in 2024, a major refactoring effort was undertaken to modernize the
GluCat codebase and address the long-standing goal of moving away from a hard
dependency on uBLAS.

The key design decisions during this phase were:

1. Adoption of the C++23 standard, enabling the use of modern features such as
`std::numbers::pi`, concepts (via `requires` expressions), and improved move
semantics.

2. Introduction of a CRTP-based (Curiously Recurring Template Pattern) matrix
abstraction layer. The new base class `matrix_base<>` (in `glucat/matrix_base.h`)
provides a unified interface that delegates to specific implementation wrappers.
To ensure template robustness and avoid issues with circular dependencies or 
complex deduction, these wrappers and the base class use explicit or trailing 
return types for all public methods (e.g., `trace()`, `norm_inf()`, `nnz()`, 
`classify_eigenvalues()`).

3. Implementation of wrappers for Eigen (`eigen_matrix_wrapper`) and Armadillo
(`arma_matrix_wrapper`), allowing GluCat to use these high-performance libraries
as backends.

4. Migration of the core `matrix_multi` and `framed_multi` classes to use these
new wrappers, effectively decoupling the Clifford algebra logic from the
underlying linear algebra library.

5. Integration of a modern unit testing framework (`doctest`) and code coverage
analysis infrastructure.

6. Hardening of the C++ API by marking cross-representation and string-based 
constructors as `explicit`. This prevents unintended implicit conversions 
that could lead to subtle performance regressions or logical errors. 
To maintain mathematical expressiveness, templated assignment and compound 
assignment operators (e.g., `=`, `+=`) are used to allow mixed-representation 
operations without requiring explicit casts, provided the scalar types match. 
This ensures that `A += B` remains interchangeable with `A = A + B`. 
Global operator templates continue to be used for mixed-type binary 
algebraic expressions.

7. Reinforcement of representation-independence. To ensure that the public API
only exposes essential Clifford algebra properties, matrix-specific
metadata functions such as nbr_rows() and nbr_cols() in matrix_multi
are defined as protected. This design choice, consistent with the previous
uBLAS-based implementation which did not expose matrix dimensions,
prevents user code from relying on internal representation details and
ensures a consistent interface across all multivector representations
(e.g., framed_multi).

8. Optimization of sparse diagonal operations. To address performance issues 
with manual diagonal write loops (which are particularly slow for compressed 
sparse formats), the `unit(rows, cols)` method was added to the backend 
wrappers, and `unit_helper` was specialized for all wrappers to ensure the 
generic `unit(dim)` free function uses backend-optimized methods (e.g., 
`setIdentity()` for Eigen, `eye()` for Armadillo). The `trace()` function 
for sparse wrappers was also optimized to use backend-specific efficient 
implementations.

9. Modernization of the public API by replacing trailing return types 
(`auto func() -> T`) with standard ones (`T func()`) across all core 
classes and matrix wrappers. This improves API reasoning and consistency 
while maintaining the stability of template instantiations. Where 
necessary for sparse matrix proxy objects (e.g., in Armadillo), 
`decltype(auto)` is used to preserve reference semantics without trailing 
syntax.


Split of code between glucat/matrix_imp.h and glucat/matrix_multi_imp.h
-----------------------------------------------------------------------

The source code of GluCat that uses linear algebra is concentrated in two files,
`glucat/matrix_imp.h` and `glucat/matrix_multi_imp.h`. The reason for the split
is to hide some of the linear algebra implementation details from the
`matrix_multi<>` class by placing them within template functions in the
namespace `matrix`. In general, functions in the namespace `matrix` operate on
operands of type `LHS_T`, `RHS_T` or `Matrix_T`, where each of these is treated
("duck typed") as a generic matrix type.

Functions such as `classify_eigenvalues()` and the underlying eigenvalue
solvers are provided by the backend-specific wrappers (Eigen or Armadillo).
This allows GluCat to accurately compute the matrix square root and
logarithm functions, which in turn are used to compute the corresponding
functions for elements of a Clifford algebra.
See [2] for further details.


The header file `glucat/matrix_multi_imp.h` is used to implement the Clifford
algebra interface of `matrix_multi<>`, including algebraic operations and
transcendental functions. Besides the scaling and squaring type algorithms used
for the square root and logarithm, Pade' approximations are used, which in turn
depend on being able to quickly calculate and divide matrix polynomials. See the
comment header of `glucat/matrix_multi_imp.h` for references.


Historically, Joerg Walter initially conducted the port of GluCat to
uBLAS. He coded `ublas::lu_factorize` and `lu_substitute` as part of the port, then
moved them into uBLAS. In `operator/()`, GluCat performs division in a Clifford
algebra by using this LU decomposition and back substitution to compute X=B/A by
solving AT*XT=BT, with a right hand side consisting of the whole square matrix
BT. The function uses iterative refinement to obtain a more accurate solution.
See the reference given in the source code.


GluCat Linear Algebra Architecture
--------------------------------

Right now, the design of GluCat depends on being able to use sparse matrices
easily. The reason why sparse matrices are used at all comes down to two things:

1. When representing an element of a Clifford algebra as a linear combination of
basis matrices, these basis matrices are usually *monomial*: they have only one
non-zero per row or per column, and are thus sparse in the algebraic sense.

2. The number of basis elements in a Clifford algebra grows exponentially with
dimension. Thus, for GluCat to deal with large and arbitrary dimensions, if
basis elements are to be stored in e.g. a cache, they must be stored
efficiently.


To support multiple backends (Eigen, Armadillo), GluCat uses a CRTP-based
abstraction layer. The `matrix_multi<>` template class no longer directly
references specific library types. Instead, it uses type selectors to choose
between the available wrappers (e.g., `eigen_matrix_wrapper` or
`arma_matrix_wrapper`).


Essentially, GluCat uses a matrix for one of two things:

1. As a matrix which the value of a generator or a basis element of a Clifford
algebra (relative to a frame).

This is the class `matrix_multi<>::basis_matrix_t`.

This is the class `matrix_multi<>::matrix_t`, `framed_multi<>::matrix_t`, etc.
To store its terms efficiently, `framed_multi` uses a high-performance hash map.
Specifically, GluCat uses `boost::unordered_flat_map` when available (requires
Boost 1.83.0+), falling back to `std::unordered_map` otherwise.

Thus, in `glucat/matrix_multi.h`, the definition of the template class `matrix_multi<>`
contains:
```
  private:
    using basis_matrix_t = matrix::sparse_matrix_t<int>;
    using matrix_t = matrix::matrix_t<Scalar_T>;
    using matrix_index_t = matrix::matrix_index_t;
```
where `matrix_t` and `sparse_matrix_t` are selected based on the scalar type and
via `--with-armadillo`). To ensure template robustness and avoid issues with 
circular dependencies or complex deduction, these wrappers use explicit return 
types for all public methods (e.g., `trace()`, `norm_inf()`, `nnz()`).
and in `glucat/framed_multi.h`, the definition of the template class `framed_multi<>`
contains:
```
    using matrix_t = typename matrix_multi_t::matrix_t;
```
In `glucat/framed_multi_imp.h`, the function `fast()` uses the class
`matrix_multi_t::basis_matrix_t`, and this class is also explicitly used in
`glucat/matrix_multi_imp.h`.

In `glucat/matrix_imp.h` and `glucat/matrix_multi_imp.h`, there are loops that
use iterators over whole matrices or over one dimension of a matrix. In the
library backends (Eigen or Armadillo), these iterate over the
structural non-zeros of a compressed matrix or over all elements of a dense
matrix.
 Thus in `glucat/matrix_imp.h` and
`glucat/matrix_multi_imp.h`, the same functions are often used for both sparse
and dense matrices. One exception is the function `matrix::mono_prod()` in
`glucat/matrix_imp.h`, which assumes that the matrix is sparse, in fact
monomial.


The reason why the generator and basis matrices are treated separately to
general matrices in GluCat comes down to a fundamental design decision, driven
by the need to support efficient multiplication of real Clifford algebra
elements with arbitrary signatures, the need to support multiplicative inverses
and division wherever possible, and the need to support geometric algebra
operations, such as the contraction, Hestenes inner product and wedge
(Grassmann) product. These requirements meant that *both* the `framed_multi_t`
class (best for geometric algebra operations) and the `matrix_multi_t` class (best
for multiplication, inverse and division) were needed in GluCat, and the
`matrix_multi_t` class had to use matrices of various different sizes. Not only
that, but the conversion between the two classes needed to be efficient for both
small signatures (e.g. Cl(1,1)) and large signatures (e.g. Cl(8,8)).


The reason why the template class `basis_matrix_t` is defined using
`matrix::sparse_matrix_t<int>` rather than using the multivector's scalar type
(e.g., `matrix::sparse_matrix_t<Scalar_T>`)
comes down to space and speed.

The type Scalar_T can be `float` (32 bits), `double` (64 bits), `long double`
(80 bits on x86 architectures), or (if the QD library is used) `dd_real` (128
bits) or `qd_real` (256 bits). Thus for any `Scalar_T` other than float, a
`basis_matrix_t` based on `int` should take up less space than one based on
`Scalar_T`. The use of `int` for `basis_matrix_t` was introduced in GluCat
0.8.2. Extensive testing showed that the new `basis_matrix_t` produces code that
is not significantly slower than the old `basis_matrix_t` that used `Scalar_T`.


Naive and fast representation and inverse representation algorithms
-------------------------------------------------------------------

As described in the paper, "A generalized FFT for Clifford algebras" [1], there
are essentially two types of algorithm that convert between the two classes,
naive and fast, with the fast algorithm being faster asymptotically as the
signature gets larger.

In GluCat, the naive algorithm makes extensive use of the class `basis_matrix_t`
and the functions `matrix::mono_prod()` and `matrix::mono_kron()`, since conversion
from `framed_multi_t` to `matrix_multi_t` relies on explicitly expanding a Clifford
algebra element as a linear combination of basis matrices, where each basis
matrix is obtained by multiplying together a a specific sequence of generator
matrices, and the generators themselves are obtained as Kronecker products of
smaller generators.

Rather than just always creating each required basis matrix on the fly, GluCat
instead uses two caches, a generator table of template class `generator_table<>`,
declared in `glucat/generation.h`, and a basis cache of template class `basis_table<>`,
declared and defined in `glucat/matrix_multi_imp.h`. Each of these two caches is
implemented as a singleton, in generator_table::generator() in `glucat/generation_imp.h`,
and in `basis_table::basis()` in `glucat/matrix_multi_imp.h`, respectively.

The key function used in the naive conversion from `framed_multi_t` to
`matrix_multi_t` is `matrix_multi_t::operator+=(term_t)`, defined in
`glucat/matrix_multi_imp.h`. This function uses the underlying matrix's
`plus_assign()` member function, along with a call to
`matrix_multi_t::basis_element()`.

Conversion in the other direction is given by the `framed_multi_t` constructor
from `matrix_multi_t`, defined in `glucat/framed_multi_imp.h`. This uses the function
`matrix::inner()`, defined in `glucat/matrix_imp.h` as well as
`matrix_multi_t::basis_element()`.


The fast algorithm to convert from `framed_multi_t` to `matrix_multi_t`, as
implemented in the function `framed_multi_t::fast()`, in
`glucat/framed_multi_imp.h`, also uses `basis_matrix_t`, but in a different way.
The key function used in `framed_multi_t::fast()` is `matrix::kron()`, defined
in `glucat/matrix_imp.h`. This is a Kronecker product with a dense result. Fast
conversion in the other direction is handled by `matrix_multi_t::fast()` in
`glucat/matrix_multi_imp.h`. This function also uses `basis_matrix_t`. It relies
on the function `matrix::signed_perm_nork()` defined in `glucat/matrix_imp.h`,
which performs a left Kronecker quotient, similar to that described in
[1, p. 679].


Choice between naive and fast representation and inverse representation algorithms
----------------------------------------------------------------------------------

In the definition of `template struct tuning<>` in `glucat/tuning.h` you will find
the following lines
```
  // Tuning for FFT
    // Minimum map size needed to invoke generalized FFT
    enum { fast_size_threshold = Fast_Size_Threshold };
    // Minimum matrix dimension needed to invoke inverse generalized FFT
    enum { inv_fast_dim_threshold = Inv_Fast_Dim_Threshold };
```
These parameters are used to tune the choice between the naive and the fast algorithms
in each direction.

The values used by a particular instance of `framed_multi_t` or `matrix_multi_t`
are set in the policy class template parameter `Tune_P`, which is assumed to be an
explicit specialization of `template struct tuning<>`. See `test/tuning.h` for examples.


In `glucat/matrix_multi_imp.h`, in the constructor of `matrix_multi_t` from
`framed_multi_t`, there appears an `if` statement starting with:
```
if (val.size() >= Tuning_Values_P::fast_size_threshold)
```
The effect of the if statement is that if the number of terms in `val` is at least
`fast_size_threshold`, then the fast algorithm is used to convert `framed_multi_t`
to `matrix_multi_t`, otherwise the naive algorithm is used.


In `glucat/framed_multi_imp.h`, in the constructor of `framed_multi_t` from
`matrix_multi_t`, appears the line:
```
const auto dim = val.m_matrix.nbr_rows();
```
In other words, `dim` is set to the first dimension of the `matrix_t` representing
`val` (and the matrix is assumed to be square).

If `dim >= Tune_P::inv_fast_dim_threshold`, then the fast transformation from
`matrix_multi_t` to `framed_multi_t` is used, otherwise the naive algorithm is used.


The definition of `template struct tuning` also contains the tuning parameter
`basis_max_count`:
```
  // Tuning for basis cache
    // Maximum index count of folded frames in basis cache
    enum { basis_max_count = Basis_Max_Count };
```
This parameter is used in the function `matrix_multi_t::basis_element()` to
determine if the basis cache will be used to construct a basis element. If the
size of the frame of the `matrix_multi_t` is at most `basis_max_count`, then the
basis cache is used.


GluCat uses two tests to assess the speed of the naive and the fast algorithms.
The test results are used to set the default values of `fast_size_threshold` and
`inv_fast_dim_threshold`.

The first test is `gfft_test`, which tests the fast transform for speed and
accuracy. The test calls the fast transform functions directly, rather than
going through the constructors. The test output looks like:
```
Clifford algebra transform test:
Fill: 0.5
R_{-11,11} in R_{-11,11}:
 CPU = mm:        0.012 fm:        0.009 trials mm:     1 fm:     1  diff: 0.00e+00
R_{-10,10} in R_{-10,10}:
 CPU = mm:        0.002 fm:        0.001 trials mm:     1 fm:     1  diff: 0.00e+00
[...]
R_{-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10} in R_{-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10}:
 CPU = mm:      931.432 fm:     1509.958 trials mm:     1 fm:     1  diff: 1.98e-16
R_{-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11} in R_{-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11}:
 CPU = mm:     5133.609 fm:     7969.989 trials mm:     1 fm:     1  diff: 2.09e-16
```
The `mm` CPU time is time to convert from `framed_multi_t` to `matrix_multi_t`.

The `fm` CPU time is time to convert from `matrix_multi_t` to `framed_multi_t`.

The `diff` is the relative difference of `abs(new_a-a)` with respect to `a`, where `a`
is the initial random `framed_multi_t` value, and `new_a` is the result of
converting to `matrix_multi_t` and back again.

The second test is `transforms`, which compares the naive and fast transforms.
The test output looks like:
```
framed_multi<double,DEFAULT_LO,DEFAULT_HI,tuning_naive_p>
Clifford algebra transform test:
Fill: 0.5
Cl( 1, 0) in Cl(16, 0) CPU = mm:       1.90 (old)       4.08 (new) fm:    3887.53 (old)       1.29 (new)  diff: old: 6.22e-15 new: 0.00e+00 fm: 0.00e+00 mm: 6.22e-15
Cl( 2, 0) in Cl(16, 0) CPU = mm:       1.40 (old)       3.93 (new) fm:    3919.84 (old)       1.19 (new)  diff: old: 5.10e-15 new: 3.36e-17 fm: 0.00e+00 mm: 5.10e-15
```
Here, `mm` is again conversion from `framed_multi_t` to `matrix_multi_t`, and `fm` is
conversion from `matrix_multi_t` to `framed_multi_t`. The naive algorithm is called
`(old)` and the fast algorithm is called `(new)`.











In future
---------

Most of the goals identified in the 2016-2022 period have been addressed by the
move to the CRTP-based architecture and the support for Eigen and Armadillo.

Remaining areas for development:

* [Eigen 3 does not depend on other linear algebra libraries,
  ](http://eigen.tuxfamily.org/index.php?title=Main_Page)
[but there is some support for CUDA.
](http://eigen.tuxfamily.org/dox-devel/TopicCUDA.html)
* [Armadillo can support OpenBLAS, NVBLAS and ACML as back ends.
  ](http://arma.sourceforge.net/faq.html#speed)

Here are some comparisons of the two:
[1](http://nghiaho.com/?p=936),
[2](http://nghiaho.com/?p=954),
[3](http://nghiaho.com/?p=1726),
[4](http://scicomp.stackexchange.com/questions/351/recommendations-for-a-usable-fast-c-matrix-library),
[5](https://web.archive.org/web/20160620155944/http://fenicsproject.org/pipermail/fenics/2013-September/000591.html).

There are now other options, but I have not had time to look at them in depth:
* [Blaze](https://bitbucket.org/blaze-lib/blaze)
* [MTL4](https://www.simunova.com/en/mtl4/)
* [ViennaCL](http://viennacl.sourceforge.net/)


References
----------

[1] Paul Leopardi, ["A generalized FFT for Clifford algebras", Bulletin of the
Belgian Mathematical Society - Simon Stevin, Volume 11, Number 5, 2005, pp.
663-688. MR 2130632.
](http://projecteuclid.org/euclid.bbms/1110205626)

[2] Paul Leopardi, ["Approximating the square root and logarithm functions in
Clifford algebras: what to do in the case of negative eigenvalues?", (extended
abstract) AGACSE 2010, June 2010.
](https://sites.google.com/site/paulleopardi/Leopardi-AGACSE10-Abstract.pdf)

Describes how the Clifford algebras over the real numbers can be treated as real
matrices, except in the case of negative real eigenvalues, when the square root
and logarithm functions may take values in a larger Clifford algebra.
