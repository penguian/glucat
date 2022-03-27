AUTHORS for GluCat 0.11.0 with PyClical
=======================================

Paul C. Leopardi <leopardi@users.sourceforge.net>

Primary references
------------------

1. Paul Leopardi, ["Practical computation with Clifford algebras"
   ](https://drive.google.com/open?id=0BynZQtPnxGf-Zm5NR1B1dXBTdTA),
   (unpublished incomplete coursework Masters major project report) UNSW 2002.

2. Paul Leopardi, ["Practical computation with Clifford algebras"
   ](https://drive.google.com/open?id=0BynZQtPnxGf-eVJJTkIyUEowOXM),
   (unpublished presentation) UNSW 2002.

3. Paul Leopardi, ["A generic library of universal Clifford algebra templates"
   ](https://drive.google.com/open?id=0BynZQtPnxGf-dkJMZDFFLTlYSTQ),
   (poster) ICCA 6, Cookville, 2002.

4. Paul Leopardi, ["A generalized FFT for Clifford algebras"
   ](https://projecteuclid.org/journalArticle/Download?urlId=10.36045%2Fbbms%2F1110205626),
   Bulletin of the Belgian Mathematical Society - Simon Stevin,
   Volume 11, Number 5, 2005, pp. 663-688.

5. Paul Leopardi, ["Approximating functions in Clifford algebras"
   ](https://drive.google.com/open?id=0BynZQtPnxGf-aXRoU1owWktyLTA),
   (unpublished presentation) ANZMC, Christchurch, 2008.

6. Paul Leopardi, ["Approximating the square root and logarithm functions in
   Clifford algebras: what to do in the case of negative eigenvalues?"
   ](https://drive.google.com/open?id=0BynZQtPnxGf-MmRxLThDOE9TdmM),
   (extended abstract) AGACSE 2010, Amsterdam, 2010.

Project supervisor
------------------

* Dr. William McLean, School of Mathematics, UNSW.

Prototype
---------

* Arvind Raja

Arvind Raja's original header comments and references follow.
```
clifford algebra package,  Arvind.Raja@hut.fi
ref: Press et.al. "Numerical Recipes in C", 2nd ed., C.U.P., 1992.
ref: LEDA, v 3.0, Stefan N\"aher, Max-Planck-Institut f\"ur Informatik
ref: Stroustrup B., "The C++ Programming Language", 2nd ed.,
     Addison-Wesley, 1991.
ref: R. Sedgewick, "Algorithms in C++", Addison-Wesley, 1992.
ref: S. Meyers, "Effective C++ ", Addison-Wesley, 1992.
```
Arvind Raja's prototype was modelled on Clical, by Pertti Lounesto, R. Mikkola,
V. Vierros.

Definitions
-----------

* Milton Abramowicz and Irene A. Stegun:
  Transcendental functions.
* Mark Ashdown, G. P. Wene:
  Negative integer indices.
* Chris Doran, David Hestenes, Anthony Lasenby, Pertti Lounesto, Ian R. Porteous, Garret Sobczyk:
  Clifford algebra sum, difference, inner product, outer product,
  geometric product, involutions and anti-involutions, norms.
* Leo Dorst:
  Left contraction.
* Leo Dorst, Arvind Raja:
  Index sets.
* Yorick Hardy, Nikos P. Pitsianis, Charles F. van Loan:
  Kronecker quotient.
* David Hestenes, Garret Sobczyk,  G. P. Wene:
  Countably infinite dimensional Clifford algebra.
* Pertti Lounesto, Ian R. Porteous:
  Matrix basis.

Algorithms
----------

* Milton Abramowicz and Irene A. Stegun, Gene H. Golub, Charles F. van Loan,
C.F. Gerald,  P.O. Wheatley:
  Transcendental functions.
* Joerg Arndt:
  Bit functions on index sets.
* Sheung Hun Cheng, Nicholas J. Higham, Charles S. Kenney, Alan J. Laub :
  Product form of Denman-Beavers square root iteration, incomplete square root cascade logarithm.
* Leo Dorst:
  Left contraction.
* P. Fleckenstein:
  `framed_multi<>` geometric product
* Gene H. Golub, Charles F. van Loan, Nicholas J. Higham:
  Matrix division, iterative refinement.
* Gene H. Golub, Charles F. van Loan:
  Matrix exponential.
* Gene H. Golub, Charles F. van Loan, C.F. Gerald, P.O. Wheatley:
  Pade approximation.
* Pertti Lounesto, Ian R. Porteous:
  Matrix basis.
* A. Lumsdaine, J. Siek, MTL project:
  Matrix algorithms and linear algebra.
* Arvind Raja:
  `framed_multi<>` sum, difference, inner product, outer product,
  geometric product, involutions and anti-involutions, norms.
* Joerg Walter, uBLAS project:
  LU factor and solve, matrix algorithms and linear algebra.
* Jan Cnops:
  Recursive expressions for matrix representations of Clifford algebras.
* Michael Clausen, Ulrich Baum, David K. Maslen, Daniel N. Rockmore:
  Generalized FFTs for finite groups.

uBLAS Interface
---------------

* Joerg Walter: uBLAS interface.

QD Interface
------------

* David H. Bailey and QD team: Changes to constructors in QD 2.3.10.

Configuration and Building
--------------------------

* Stephan Kulow, Walter Tasin, KDevelop team:
  Original code used in `admin/*` and `configure.ac.in`
* John Calcote:
  Autools examples and ideas for `make check`, `make dist`, `make doc`, and `make install`,
  used in `configure.ac.in`, `Makefile.am.in` and various `Makefile.am` files.
* Chris Miceli:
  Integrating Doxygen with Autotools.
* [Clang-Tidy:](https://clang.llvm.org/extra/clang-tidy/)
  C++ code modernization

PyClical Python extension module
--------------------------------

* Paul Zimmermann:
  Hosting Sage Days 10 in Nancy, 2008
* William Stein:
  Sage Days 10 in Nancy, 2008
* Attendees of Sage Days 10 in Nancy, 2008:
  Advice on use of Cython with C++.
* Robert Bradshaw, Mark Florisson, cython-users:
  More advice on use of Cython with C++.
* Brian Quinlan:
  Code in `pyclical/setup_ext.py` adapted from
  [setup.py](https://github.com/SublimeCodeIntel/CodeIntel/blob/master/silvercity/setup.py).

PyClical Tutorials
------------------

* Perttu Puska:
  Permission to re-implement the CLICAL tutorials in PyClical.

Testing and patches
-------------------

* Alan Bromborsky, Johannes Brunen, John Fletcher, Henk Jansen, Joerg Walter,
Patrick Welche.

openSUSE Build Service
----------------------

* Atri Bhattacharya, Adrian Schröter:
  [openSUSE Build Service](https://build.opensuse.org/package/show/science/glucat).

Workarounds
-----------

* Carlos O'Ryan: workaround used in `generator.h` for:
"... only defines a private destructor and has no friends"

Advice
------

* A. Alexandrescu, N. Josuttis, A. Koenig, Scott Meyers, B. Moo,
Bjarne Stroustrup, T. Veldhuizen:
  C++ usage.
* Joerg Arndt:
  Bit wizardry.
* Paul Ivanov:
  IPython notebooks.
* Duraid Madina, Russell Standish:
  C++, debugging.
* Christian Perwass:
  Operators, user interface.

Related Projects
----------------

* John Fletcher: [BoostCliffordDiscussion
  ](http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?BoostCliffordDiscussion)

Specific references
-------------------

### References for identities

1. L. Dorst, ["Honing geometric algebra for its use in the computer sciences"
](https://staff.fnwi.uva.nl/l.dorst/clifford/sommer.pdf),
in Geometric Computing with Clifford Algebras, (G. Sommer, ed.)
Springer 2001, Chapter 6, pp. 127-152.

2. L. Dorst, ["The inner products of geometric algebra"
](https://staff.fnwi.uva.nl/l.dorst/clifford/inner.ps),
in Applications of Geometric Algebra in Computer Science and Engineering
(Dorst, Doran, Lasenby, eds), Birkhauser, 2002.

3. D. Hestenes, G. Sobczyk, Clifford Algebra to Geometric Calculus, D. Reidel, 1984.

### References for generalized FFTs on finite groups

1. L. Babai, L. Ronyai, "Computing irreducible representations of finite groups",
Math. Comp. 55, 1990, pp. 705--722.

2. U. Baum, M. Clausen, "Computing irreducible representations of supersolvable
groups", Math. Comp. 63, 1994, pp. 351--359.

3. T. Beth, "On the computational complexity of the general discrete Fourier
transform", Theoretical Computer Science, 51 (3) 1987, pp. 331--339.

4. G. S. Chirikjian, A. B. Kyatkin, Engineering applications of noncommutative
harmonic analysis. With emphasis on rotation and motion groups. CRC Press, 2001.

5. M. Clausen, "Fast generalized Fourier transforms", Theoretical Computer
Science 67, 1989, pp. 55--63.

6. M. Clausen, U. Baum, "Fast Fourier Transforms," Bibliographisches Institut &
F. A. Brockhaus AG, Mannheim, 1993.

7. P. Diaconis, D. Rockmore, "Efficient computation of the Fourier transform
of finite groups," J. Amer. Math. Soc., 3 (2), 1990, pp. 297--332.

8. D. Maslen, D. Rockmore, "Generalized FFTs: A survey of some recent results",
in Proc. 1995 DIMACS Workshop in Groups and Computation, L. Finkelstein
and W. Kantor (eds.).

9. D. Maslen, D. Rockmore, "Separation of variables and the computation of
Fourier transforms on finite groups", I. J. of the Amer. Math. Soc., 10,
No. 1, 1997, pp. 169--214.

10. D. Rockmore, "Some applications of generalized FFTs" (Appendix with D.
Healy), in Groups and Computation II, DIMACS Series in Discrete Math and
Computer Science, Vol. 28, L. Finkelstein and W. Kantor (eds.), 1997, pp.
329--369.

### References for Kronecker quotients

1. C. F. van Loan and N. Pitsianis, "Approximation with Kronecker products",
in Linear Algebra for Large Scale and Real-Time Applications, Marc S. Moonen,
Gene H. Golub, and Bart L. R. Moor (eds.), 1993, pp. 293--314.

2. Y. Hardy, "On Kronecker quotients", Electronic Journal of Linear Algebra,
Vol. 27, 2014, pp. 172--189.

Recommended general references
------------------------------

### Books

1. Rafal Ablamowicz, Pertti Lounesto, Josep Parra (eds.),
"Clifford Algebras with Numeric and Symbolic Computations", Birkhauser, 1996.

2. Chris Doran and Anthony Lasenby,
["Geometric Algebra for Physicists",  Cambridge University Press, 2003, 2007.
](http://geometry.mrao.cam.ac.uk/2007/01/geometric-algebra-for-physicists/)

3. David Hestenes, Garret Sobczyk, "Clifford Algebra to Geometric Calculus",
Kluwer, first published in 1984; reprinted with corrections in 1992.

4. Pertti Lounesto, "Clifford Algebras and Spinors", Cambridge University Press,
Second enhanced edition, 2001.

5. Alan Macdonald, "Linear and Geometric Algebra", 2011.

6. Ian R. Porteous, "Clifford Algebras and the Classical groups", Cambridge
University Press, 1995.

### Web sites

* [Rafal Ablamowitz - archived
  ](https://web.archive.org/web/20210211165810/https://math.tntech.edu/rafal/)
* [William Baylis
  ](http://cronus.uwindsor.ca/baylis-research)
* [The Clifford Research Group at Ghent University - archived
  ](https://web.archive.org/web/20190917051859/http://cage.ugent.be/crg/)
* [Leo Dorst
  ](https://staff.fnwi.uva.nl/l.dorst/clifford/index.html)
* [The Geometric Algebra Research Group at Cambridge
  ](http://geometry.mrao.cam.ac.uk/)
* [David Hestenes
  ](http://geocalc.clas.asu.edu/)
* [Soeren Krausshar
  ](https://www.uni-erfurt.de/erziehungswissenschaftliche-fakultaet/fakultaet/profil/fachgebiete-und-professuren/mathematik-und-mathematikdidaktik/team/soeren-krausshar)
* [Pertti Lounesto - archived
  ](https://users.aalto.fi/~ppuska/mirror/Lounesto/)
* [Alan Macdonald
  ](http://faculty.luther.edu/~macdonal/laga/)
* [Christian Perwass
  ](http://www.researchgate.net/profile/Christian_Perwass/publications)
* [Wikipedia: Geometric algebra
  ](http://en.wikipedia.org/wiki/Geometric_algebra)


Other Clifford algebra software
-------------------------------

### Current
 * [GABLE, a Matlab package with a geometric algebra tutorial
   ](https://staff.fnwi.uva.nl/l.dorst/GABLE/index.html)
 * [Gaigen 2.5, a C++ code optimized geometric algebra implementation generator
   ](http://sourceforge.net/projects/g25/)
 * [galgebra: geometric algebra/calculus modules for sympy
   ](https://github.com/pygae/galgebra)

### Historical
 * [CLICAL, a Clifford algebra calculator - archived
   ](https://users.aalto.fi/~ppuska/mirror/Lounesto/CLICAL.htm)
 * [CLIFFORD 2017 for Maple 2017
   ](https://web.archive.org/web/20171004075723/http://math.tntech.edu/rafal/cliff2017/index.html)
 * [CLUCalc, 3D calculation and visualization software
   ](https://github.com/foobarz/CLUCalcSource-4.3.3)
 * [GA Package for Maple
   ](http://www.mrao.cam.ac.uk/~maja1/software/GA/)

### swMATH
 * https://www.swmath.org/?term=clifford+algebras
 
