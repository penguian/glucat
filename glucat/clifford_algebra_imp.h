#ifndef _GLUCAT_CLIFFORD_ALGEBRA_IMP_H
#define _GLUCAT_CLIFFORD_ALGEBRA_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    clifford_algebra_imp.h : Implement common Clifford algebra functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2007 by Paul C. Leopardi
 ***************************************************************************

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library.  If not, see <http://www.gnu.org/licenses/>.

 ***************************************************************************
 This library is based on a prototype written by Arvind Raja and was
 licensed under the LGPL with permission of the author. See Arvind Raja,
 "Object-oriented implementations of Clifford algebras in C++: a prototype",
 in Ablamowicz, Lounesto and Parra (eds.)
 "Clifford algebras with numeric and symbolic computations", Birkhauser, 1996.
 ***************************************************************************
     See also Arvind Raja's original header comments in glucat.h
 ***************************************************************************/

// References for algorithms:
// [AS]:
// Milton Abramowicz and Irene A. Stegun, "Handbook of mathematical functions",
// Dover 1972, first published 1965.
// [CHKL]:
// Sheung Hun Cheng, Nicholas J. Higham, Charles S. Kenney and Alan J. Laub,
// "Approximating the Logarithm of a Matrix to Specified Accuracy", 1999.
// ftp://ftp.ma.man.ac.uk/pub/narep/narep353.ps.gz
// [GL]:
// Gene H. Golub and Charles F. van Loan,
// "Matrix Computations", 3rd ed., Johns Hopkins UP, 1996.
// [GW]:
// C.F. Gerald and P.O. Wheatley, "Applied Numerical Analysis",
// 6th Edition, Addison-Wesley, 1999.
// [H]:
// Nicholas J. Higham 
// "The Scaling and Squaring Method for the Matrix Exponential Revisited",
// SIAM Journal on Matrix Analysis and Applications,
// Vol. 26,  Issue 4 (2005), pp. 1179-1193.
// [Z]:
// Doron Zeilberger, "PADE" (Maple code), 2002.
// http://www.math.rutgers.edu/~zeilberg/tokhniot/PADE

namespace glucat
{
  template< typename Scalar_T, typename Index_Set_T, typename Multivector_T>
  const std::string
  clifford_algebra<Scalar_T,Index_Set_T,Multivector_T>::
  classname()
  { return "clifford_algebra"; }

  /// Test for inequality of multivectors
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  inline
  bool
  operator!= (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs)
  { return !(lhs == rhs); }


  /// Test for inequality of multivector and scalar
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  bool
  operator!= (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr)
  { return lhs.frame().count() != 0 || scalar(lhs) != scr; }

  /// Test for inequality of scalar and multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  bool
  operator!= (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs)
  { return rhs.frame().count() != 0 || scalar(rhs) != scr; }

  /// Geometric sum of multivector and scalar
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator+ (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr)
  {
    Multivector<Scalar_T,LO,HI> result = lhs;
    return result += scr;
  }

  /// Geometric sum of scalar and multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator+ (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs)
  {
    return rhs + scr;
  }

  /// Geometric sum
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator+ (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs)
  {
    Multivector<Scalar_T,LO,HI> result = lhs;
    return result += rhs;
  }

  /// Geometric difference of multivector and scalar
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator- (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr)
  {
    Multivector<Scalar_T,LO,HI> result = lhs;
    return result += -scr;
  }

  /// Geometric difference of scalar and multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator- (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs)
  { return -rhs + scr; }

  /// Geometric difference
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator- (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs)
  {
    Multivector<Scalar_T,LO,HI> result = lhs;
    return result -= rhs;
  }

  /// Product of multivector and scalar
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator* (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr)
  {
    Multivector<Scalar_T,LO,HI> result = lhs;
    return result *= scr;
  }

  /// Product of scalar and multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator* (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs)
  { // Note: this assumes that scalar commutes with multivector.
    // This excludes Clifford algebras over non-commuting rings.
    return rhs * scr;
  }

  /// Geometric product
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator* (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs)
  {
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    return lhs * multivector_t(rhs);
  }

  /// Outer product
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator^ (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs)
  {
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    return lhs ^ multivector_t(rhs);
  }

  /// Inner product
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator& (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs)
  {
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    return lhs & multivector_t(rhs);
  }

  /// Left contraction
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator% (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs)
  {
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    return lhs % multivector_t(rhs);
  }

  /// Hestenes scalar product
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  inline
  Scalar_T
  star (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs)
  {
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    return star(lhs, multivector_t(rhs));
  }

  /// Quotient of multivector and scalar
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator/ (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr)
  {
    Multivector<Scalar_T,LO,HI> result = lhs;
    return result /= scr;
  }

  /// Quotient of scalar and multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator/ (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs)
  {
    Multivector<Scalar_T,LO,HI> result = scr;
    return result /= rhs;
  }

  /// Geometric quotient
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  inline
  const Multivector<Scalar_T,LO,HI>
  operator/ (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs)
  {
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    return lhs / multivector_t(rhs);
  }

  /// Geometric multiplicative inverse
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  inv(const Multivector<Scalar_T,LO,HI>& val)
  { return val.inv(); }

  /// Main involution, each {i} is replaced by -{i} in each term, eg. {1}*{2} -> (-{2})*(-{1})
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  involute(const Multivector<Scalar_T,LO,HI>& val)
  { return val.involute(); }

  /// Reversion, eg. {1}*{2} -> {2}*{1}
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  reverse(const Multivector<Scalar_T,LO,HI>& val)
  { return val.reverse(); }

  /// Conjugation, rev o invo == invo o rev
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  conj(const Multivector<Scalar_T,LO,HI>& val)
  { return val.conj(); }

  /// Scalar_T quadratic form == (rev(x)*x)(0)
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  quad(const Multivector<Scalar_T,LO,HI>& val)
  { return val.quad(); }

  /// Scalar part
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  scalar(const Multivector<Scalar_T,LO,HI>& val)
  { return val.scalar(); }

  /// Pure part
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  pure(const Multivector<Scalar_T,LO,HI>& val)
  { return val - val.scalar(); }

  /// Even part
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  even(const Multivector<Scalar_T,LO,HI>& val)
  { return val.even(); }

  /// Odd part
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  odd(const Multivector<Scalar_T,LO,HI>& val)
  { return val.odd(); }

  /// Vector part of multivector, as a vector_t with respect to frame()
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const std::vector<Scalar_T>
  vector_part(const Multivector<Scalar_T,LO,HI>& val)
  { return val.vector_part(); }

  /// Absolute value == sqrt(norm)
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  abs(const Multivector<Scalar_T,LO,HI>& val)
  { return numeric_traits<Scalar_T>::sqrt(val.norm()); }

  /// Scalar_T norm == sum of norm of coordinates
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  norm(const Multivector<Scalar_T,LO,HI>& val)
  { return val.norm(); }

  /// Real part of scalar part
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  real(const Multivector<Scalar_T,LO,HI>& val)
  { return std::real(scalar(val)); }

  /// Imaginary part of scalar part
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  imag(const Multivector<Scalar_T,LO,HI>& val)
  { return std::imag(scalar(val)); }

  /// Integer power of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  pow(const Multivector<Scalar_T,LO,HI>& lhs, int rhs)
  { 
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    multivector_t a;
    if (rhs < 0)
    {
      if (lhs == Scalar_T(0))
        return numeric_traits<Scalar_T>::NaN();
      rhs = -rhs;
      a = lhs.inv();
    }
    else
      a = lhs;
    multivector_t result = Scalar_T(1);
    for (;
        rhs != 0;
        rhs >>= 1, a *= a)
      if (rhs & 1)
        result *= a;
    return result;
  } 

  /// Multivector power of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  inline
  const Multivector<Scalar_T,LO,HI>
  pow(const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs)
  { 
    if (lhs == Scalar_T(0))
    {
      const Scalar_T m = rhs.scalar();
      if (rhs == m)
        return (m < 0)
               ? numeric_traits<Scalar_T>::NaN()
               : (m == 0)
                 ? Scalar_T(1)
                 : Scalar_T(0);
      else
        return Scalar_T(0);
    }
    return exp(log(lhs) * rhs); 
  }

  /// Square root of -1 which commutes with all members of the frame of the given multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  elliptic(const Multivector<Scalar_T,LO,HI>& val)
  {
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    typedef typename multivector_t::index_set_t index_set_t;
    index_set_t frm = val.frame();
    index_t incp[] = {0, 2, 1, 0};
    index_t incq[] = {1, 0, 0, 0};
    index_t bott = pos_mod((frm.count_pos() - frm.count_neg()), 4);
    for (index_t 
        k = 0; 
        k != incp[bott]; 
        k++)
      for (index_t 
          idx = 1; 
          idx != HI+1; 
          ++idx)
        if (!frm[idx])
        {
          frm.set(idx);
          break;
        }
    for (index_t 
        k = 0; 
        k != incq[bott]; 
        k++)
      for (index_t 
          idx = -1; 
          idx != LO-1; 
          --idx)
        if (!frm[idx])
        {
          frm.set(idx);
          break;
        }
    index_t new_bott = pos_mod(frm.count_pos() - frm.count_neg(), 4);

    if ((incp[new_bott] == 0) && (incq[new_bott] == 0))
      return multivector_t(frm, Scalar_T(1));
    else
      // Return IEEE NaN or -Inf
      return numeric_traits<Scalar_T>::NaN();
  }

  /// Check that i is a valid complexifier for val
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  static
  void
  check_complex(const Multivector<Scalar_T,LO,HI>& val, const Multivector<Scalar_T,LO,HI>& i)
  {
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    typedef typename multivector_t::error_t error_t;
    typedef typename multivector_t::index_set_t index_set_t;
    const index_set_t i_frame = i.frame();
    // We need i to be a complexifier whose frame is large enough to represent val
    if (elliptic(i) != i ||
       (val.frame() | i_frame) != i_frame ||
        elliptic(val).frame().count() > i_frame.count())
      throw error_t("check_complex(val, i): i is not a valid complexifier for val");
  }

  /// Pade' approximation
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  static
  const Multivector<Scalar_T,LO,HI>
  pade_approx(const Scalar_T a[], const Scalar_T b[], const Multivector<Scalar_T,LO,HI>& X)
  {
    // Pade' approximation
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [GL], Section 11.3, p572-576.

    typedef Multivector<Scalar_T,LO,HI> multivector_t;

    if (X.isnan())
      return numeric_traits<Scalar_T>::NaN();

    const multivector_t& X2 = X*X;
    const multivector_t& X4 = X2*X2;
    const multivector_t& X6 = X4*X2;
    const multivector_t& N = a[0]+X2*a[2]+X4*a[4]+X6*(a[6]+X2*a[8]+X4*a[10]) +
                          X*(a[1]+X2*a[3]+X4*a[5]+X6*(a[7]+X2*a[9]));
    const multivector_t& D = b[0]+X2*b[2]+X4*b[4]+X6*(b[6]+X2*b[8]+X4*b[10]) +
                          X*(b[1]+X2*b[3]+X4*b[5]+X6*(b[7]+X2*b[9]));
    return N / D;
  }

  /// Single step of product form of Denman-Beavers square root iteration
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  static
  void
  db_step(Multivector<Scalar_T,LO,HI>& M, Multivector<Scalar_T,LO,HI>& Y)
  {
    // Reference: [CHKL]
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    const multivector_t& iM = inv(M);
    M = (Scalar_T(1) + (M + iM)/Scalar_T(2)) / Scalar_T(2);
    Y *= (Scalar_T(1) + iM) / Scalar_T(2);
  }

  /// Product form of Denman-Beavers square root iteration
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  static
  const Multivector<Scalar_T,LO,HI>
  db_sqrt(const Multivector<Scalar_T,LO,HI>& val)
  {
    // Reference: [CHKL]
    typedef Multivector<Scalar_T,LO,HI> multivector_t;

    if (val == 0)
      return val;

    multivector_t M = val;
    multivector_t Y = val;
    const Scalar_T tol = std::numeric_limits<Scalar_T>::epsilon() / Scalar_T(2);
    const Scalar_T tol2 = tol * tol;
    int step = 0;
    for (; 
        step != Tune_P::sqrt_max_steps && !(Y.isnan()) && norm(Scalar_T(1) - M) > tol2; 
        ++step)
      db_step(M, Y);
    if (Y.isnan() || (step == Tune_P::sqrt_max_steps && norm(Scalar_T(1) - M) > tol2))
      return numeric_traits<Scalar_T>::NaN();

    return Y;
  }

  /// Square root of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  static
  const Multivector<Scalar_T,LO,HI>
  sqrt(const Multivector<Scalar_T,LO,HI>& val, const Multivector<Scalar_T,LO,HI>& i, bool prechecked)
  {
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [GL], Section 11.3, p572-576
    // Reference: [Z], Pade1

    // Pade approximation produced by Pade1(sqrt(1+x),x,10,10):
    // numer := 1+21/4*x+189/16*x^2+119/8*x^3
    //          +735/64*x^4+5733/1024*x^5+7007/4096*x^6+1287/4096*x^7
    //          +2079/65536*x^8+385/262144*x^9+21/1048576*x^10;
    // denom := 1+19/4*x+153/16*x^2+85/8*x^3
    //          +455/64*x^4+3003/1024*x^5+3003/4096*x^6+429/4096*x^7
    //          +495/65536*x^8+55/262144*x^9+1/1048576*x^10;

    static const Scalar_T a[] =
    {
         1.0,          21.0/4.0,      189.0/16.0,     119.0/8.0,
       735.0/64.0,   5733.0/1024.0,  7007.0/4096.0,  1287.0/4096.0,
      2079.0/65536.0, 385.0/262144.0,  21.0/1048576.0
    };
    static const Scalar_T b[] =
    {
         1.0,          19.0/4.0,      153.0/16.0,      85.0/8.0,
       455.0/64.0,   3003.0/1024.0,  3003.0/4096.0,   429.0/4096.0,
       495.0/65536.0,  55.0/262144.0,   1.0/1048576.0
    };

    if (val == 0)
      return val;

    if (val.isnan())
      return numeric_traits<Scalar_T>::NaN();

    const Scalar_T& realval = real(val);
    if (val == realval && realval < Scalar_T(0))
    {
      if (!prechecked)
        check_complex(val, i);
      return i * numeric_traits<Scalar_T>::sqrt(-realval);
    }
    // Scale val towards circle abs(val) == 1
    const Scalar_T& scale = abs(val);
    const Scalar_T& sqrt_scale = numeric_traits<Scalar_T>::sqrt(scale);
    if (numeric_traits<Scalar_T>::isNaN_or_isInf(sqrt_scale))
      return numeric_traits<Scalar_T>::NaN();

    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    const multivector_t& unitval = val / scale;
    if (unitval == 0)
      return numeric_traits<Scalar_T>::NaN();

    if (norm(unitval - Scalar_T(1)) < Scalar_T(1.0/4.0))
      // Pade' approximation of square root
      return pade_approx(a, b, unitval - Scalar_T(1)) * sqrt_scale;
    else
    {
      // Product form of Denman-Beavers square root iteration
      multivector_t scaled_result = db_sqrt(unitval);
      if (!scaled_result.isnan())
        return scaled_result * sqrt_scale;
      else
      {
        if (!prechecked)
          check_complex(val, i);
        scaled_result = db_sqrt(-unitval);
        if (!scaled_result.isnan())
          return scaled_result * sqrt_scale * i;
        else
          return numeric_traits<Scalar_T>::NaN();
      }
    }
  }

  /// Square root of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  sqrt(const Multivector<Scalar_T,LO,HI>& val)
  { return sqrt(val, elliptic(val), true); }

  /// Exponential of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  exp(const Multivector<Scalar_T,LO,HI>& val)
  {
    // Scaling and squaring Pade' approximation of matrix exponential
    // Reference: [GL], Section 11.3, p572-576
    // Reference: [H]

    typedef Multivector<Scalar_T,LO,HI> multivector_t;

    const Scalar_T scalar_val = scalar(val);
    const Scalar_T scalar_exp = Scalar_T(std::exp(scalar_val));
    if (numeric_traits<Scalar_T>::isNaN_or_isInf(scalar_exp))
      return numeric_traits<Scalar_T>::NaN();

    if (val.isnan())
      return numeric_traits<Scalar_T>::NaN();

    multivector_t A = val - scalar_val;
    const Scalar_T pure_scale = A.norm();
    if (numeric_traits<Scalar_T>::isNaN_or_isInf(pure_scale))
      return numeric_traits<Scalar_T>::NaN();

    const int ilog2_scale = 
      std::max(0, int(ceil((log2(pure_scale) + Scalar_T(A.frame().count()))/Scalar_T(2)) - 3));
    const Scalar_T i_scale = Scalar_T(std::pow(Scalar_T(2), ilog2_scale));
    if (numeric_traits<Scalar_T>::isNaN_or_isInf(i_scale))
      return numeric_traits<Scalar_T>::NaN();

    A /= i_scale;
    multivector_t pure_exp;
    {
      const int q = 13;
      static Scalar_T c[q+1];
      if (c[0] != Scalar_T(1))
      {
        c[0] = Scalar_T(1);
        for (int 
            k = 0; 
            k != q; 
            ++k)
          c[k+1] = c[k]*(q-k) / ((2*q-k)*(k+1));
      }
      const multivector_t& A2 = A*A;
      const multivector_t& A4 = A2*A2;
      const multivector_t& A6 = A4*A2;
      const multivector_t& U =     c[0]+A2*c[2]+A4*c[4]+A6*(c[6]+A2*c[8]+A4*c[10]+A6*c[12]);
      const multivector_t& AV = A*(c[1]+A2*c[3]+A4*c[5]+A6*(c[7]+A2*c[9]+A4*c[11]+A6*c[13]));
      pure_exp = (U+AV) / (U-AV);
    }
    for (int 
        k = 0; 
        k != ilog2_scale; 
        ++k)
      pure_exp *= pure_exp;
    return pure_exp * scalar_exp;
  }

  /// Scaled Pade' approximation of log
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  static
  const Multivector<Scalar_T,LO,HI>
  pade_log(const Multivector<Scalar_T,LO,HI>& val)
  {
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [CHKL]
    // Reference: [GL], Section 11.3, p572-576
    // Reference: [Z], Pade1

    // Pade approximation produced by Pade1(log(1-x),x,10,10):
    // numer := -x+9/2*x^2-484/57*x^3
    //          +497/57*x^4-34167/6460*x^5+14917/7752*x^6-2761/6783*x^7
    //          +831/18088*x^8-4861/2116296*x^9+671/21162960*x^10;
    // denom := 1-5*x+405/38*x^2-240/19*x^3
    //          +2940/323*x^4-1323/323*x^5+735/646*x^6-60/323*x^7
    //          +135/8398*x^8-5/8398*x^9+1/184756*x^10;


    static const Scalar_T a[] =
    {
        0.0,            -1.0,              9.0/2.0,       -484.0/57.0,
      497.0/57.0,   -34167.0/6460.0,   14917.0/7752.0,   -2761.0/6783.0,
      831.0/18088.0, -4861.0/2116296.0,  671.0/21162960.0
    };
    static const Scalar_T b[] =
    {
         1.0,           -5.0,            405.0/38.0,      -240.0/19.0,
      2940.0/323.0,  -1323.0/323.0,      735.0/646.0,      -60.0/323.0,
       135.0/8398.0,    -5.0/8398.0,       1.0/184756.0
    };
    typedef Multivector<Scalar_T,LO,HI> multivector_t;

    if (val == 0 || val.isnan())
      return numeric_traits<Scalar_T>::NaN();

    // Scale val towards circle abs(A) == 1
    Scalar_T scale = abs(val);

    return pade_approx(a, b, Scalar_T(1) - val / scale) + Scalar_T(std::log(scale)) ;
  }

  /// Incomplete square root cascade and scaled Pade' approximation of log
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  static
  const Multivector<Scalar_T,LO,HI>
  cascade_log(const Multivector<Scalar_T,LO,HI>& val)
  {
    // Reference: [CHKL]
    typedef Multivector<Scalar_T,LO,HI> multivector_t;

    if (val == 0 || val.isnan())
      return numeric_traits<Scalar_T>::NaN();

    multivector_t Y = val;
    multivector_t E = Scalar_T(0);
    int outer_step;
    for (outer_step = 0;
        outer_step != Tune_P::log_max_outer_steps &&
        norm(Scalar_T(1) - Y) > Scalar_T(1.0/8.0);
        ++outer_step)
    {
      if (Y == 0 || Y.isnan())
        return numeric_traits<Scalar_T>::NaN();

      // Incomplete product form of Denman-Beavers square root iteration
      multivector_t M = Y;
      for (int 
          inner_step = 0;
          inner_step != Tune_P::log_max_inner_steps && 
          norm(Scalar_T(1) - M) > 0;
          ++inner_step)
        db_step(M, Y);
      E += E + Scalar_T(1) - M;
    }
    if (outer_step == Tune_P::log_max_outer_steps && norm(Scalar_T(1) - Y) > Scalar_T(1.0/8.0))
      std::cerr << "Warning: log iteration did not converge. norm = "
                << norm(Scalar_T(1) - Y) << std::endl;
    return pade_log(Y) * Scalar_T(std::pow(Scalar_T(2), outer_step)) + E;
  }

  /// Natural logarithm of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  log(const Multivector<Scalar_T,LO,HI>& val, const Multivector<Scalar_T,LO,HI>& i, bool prechecked)
  {
    // Scaled incomplete square root cascade and scaled Pade' approximation of log
    // Reference: [CHKL]

    if (val == 0 || val.isnan())
      return numeric_traits<Scalar_T>::NaN();

    const Scalar_T& realval = real(val);
    if (val == realval && realval < Scalar_T(0))
    {
      if (!prechecked)
        check_complex(val, i);
      return Scalar_T(std::log(-realval)) + i * Scalar_T(l_pi);
    }
    // Scale val towards abs(A) == 1
    const Scalar_T& scale = abs(val);
    const Scalar_T& log_scale = std::log(scale);
    if (val == scale)
      return log_scale;
    else
    {
      typedef Multivector<Scalar_T,LO,HI> multivector_t;
      const multivector_t unitval = val/scale;
      if (inv(unitval).isnan())
        return numeric_traits<Scalar_T>::NaN();
      multivector_t scaled_result = cascade_log(unitval);
      if (!scaled_result.isnan())
        return scaled_result + log_scale;
      else
      {
        if (!prechecked)
          check_complex(val, i);
        scaled_result = cascade_log(-unitval);
        if (!scaled_result.isnan())
          return scaled_result + log_scale + i * Scalar_T(l_pi);
        else
          return numeric_traits<Scalar_T>::NaN();
      }
    }
  }

  /// Natural logarithm of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  log(const Multivector<Scalar_T,LO,HI>& val)
  { return log(val, elliptic(val), true); }

  /// Hyperbolic cosine of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  cosh(const Multivector<Scalar_T,LO,HI>& val)
  { 
    return val.isnan()
         ? numeric_traits<Scalar_T>::NaN()
         : (exp(val)+exp(-val)) / Scalar_T(2); 
  }

  /// Inverse hyperbolic cosine of multivector
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  acosh(const Multivector<Scalar_T,LO,HI>& val, const Multivector<Scalar_T,LO,HI>& i, bool prechecked)
  { 
    if (!prechecked)
      check_complex(val, i);

    return val.isnan()
         ? numeric_traits<Scalar_T>::NaN()
         : log(val + sqrt(val*val - Scalar_T(1), i, true), i, true); 
  }

  /// Inverse hyperbolic cosine of multivector
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  acosh(const Multivector<Scalar_T,LO,HI>& val)
  { return acosh(val, elliptic(val), true); }

  /// Cosine of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  cos(const Multivector<Scalar_T,LO,HI>& val, const Multivector<Scalar_T,LO,HI>& i, bool prechecked)
  {
    typedef Multivector<Scalar_T,LO,HI> multivector_t;

    if (val.isnan())
      return numeric_traits<Scalar_T>::NaN();

    if (!prechecked)
      check_complex(val, i);

    const Scalar_T& s = scalar(val);
    const Scalar_T& twopi = Scalar_T(2) * Scalar_T(l_pi);
    const multivector_t& z = i *
      (val - s + Scalar_T(std::fmod(Scalar_T(std::real(s)), twopi)) +
                 Scalar_T(std::imag(s)));
    return (exp(z)+exp(-z)) / Scalar_T(2);
  }

  /// Cosine of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  cos(const Multivector<Scalar_T,LO,HI>& val)
  { return cos(val, elliptic(val), true); }

  /// Inverse cosine of multivector
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  acos(const Multivector<Scalar_T,LO,HI>& val, const Multivector<Scalar_T,LO,HI>& i, bool prechecked)
  { 
    if (!prechecked)
      check_complex(val, i);

    return val.isnan()
         ? numeric_traits<Scalar_T>::NaN()
         : i * acosh(val, i, true); 
  }

  /// Inverse cosine of multivector
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  acos(const Multivector<Scalar_T,LO,HI>& val)
  { return acos(val, elliptic(val), true); }

  /// Hyperbolic sine of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  sinh(const Multivector<Scalar_T,LO,HI>& val)
  { 
    return val.isnan()
           ? numeric_traits<Scalar_T>::NaN()
           : (exp(val)-exp(-val)) / Scalar_T(2); 
  }

  /// Inverse hyperbolic sine of multivector
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  asinh(const Multivector<Scalar_T,LO,HI>& val, const Multivector<Scalar_T,LO,HI>& i, bool prechecked)
  { 
    if (!prechecked)
      check_complex(val, i);

    return val.isnan()
         ? numeric_traits<Scalar_T>::NaN()
         : log(val + sqrt(val*val + Scalar_T(1), i, true), i, true); 
  }

  /// Inverse hyperbolic sine of multivector
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  asinh(const Multivector<Scalar_T,LO,HI>& val)
  { return asinh(val, elliptic(val), true); }

  /// Sine of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  sin(const Multivector<Scalar_T,LO,HI>& val, const Multivector<Scalar_T,LO,HI>& i, bool prechecked)
  {
    if (val.isnan())
      return numeric_traits<Scalar_T>::NaN();
    if (!prechecked)
      check_complex(val, i);
    return i * (cos(val) - exp(i * val));
  }

  /// Sine of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  sin(const Multivector<Scalar_T,LO,HI>& val)
  { return sin(val, elliptic(val), true); }

  /// Inverse sine of multivector
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  asin(const Multivector<Scalar_T,LO,HI>& val, const Multivector<Scalar_T,LO,HI>& i, bool prechecked)
  {
    if (!prechecked)
      check_complex(val, i);

    return val.isnan()
         ? numeric_traits<Scalar_T>::NaN()
         : -i * asinh(i * val, i, true);
  }

  /// Inverse sine of multivector
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  asin(const Multivector<Scalar_T,LO,HI>& val)
  { return asin(val, elliptic(val), true); }

  /// Hyperbolic tangent of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  tanh(const Multivector<Scalar_T,LO,HI>& val)
  { 
    return val.isnan()
         ? numeric_traits<Scalar_T>::NaN()
         : sinh(val) / cosh(val); 
  }

  /// Inverse hyperbolic tangent of multivector
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  atanh(const Multivector<Scalar_T,LO,HI>& val, const Multivector<Scalar_T,LO,HI>& i, bool prechecked)
  { 
    if (!prechecked)
      check_complex(val, i);

   return val.isnan()
        ? numeric_traits<Scalar_T>::NaN()
        : log((Scalar_T(1) + val) / (Scalar_T(1) - val), i, true) / Scalar_T(2); 
  }

  /// Inverse hyperbolic tangent of multivector
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  atanh(const Multivector<Scalar_T,LO,HI>& val)
  { return atanh(val, elliptic(val), true); }

  /// Tangent of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  tan(const Multivector<Scalar_T,LO,HI>& val, const Multivector<Scalar_T,LO,HI>& i, bool prechecked)
  { 
    if (!prechecked)
      check_complex(val, i);

    return val.isnan()
         ? numeric_traits<Scalar_T>::NaN()
         : sin(val, i, true) / cos(val, i, true); 
  }

  /// Tangent of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  tan(const Multivector<Scalar_T,LO,HI>& val)
  { return tan(val, elliptic(val), true); }

  /// Inverse tangent of multivector
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  atan(const Multivector<Scalar_T,LO,HI>& val, const Multivector<Scalar_T,LO,HI>& i, bool prechecked)
  {
    if (!prechecked)
      check_complex(val, i);

    return val.isnan()
          ? numeric_traits<Scalar_T>::NaN()
          : -i * atanh(i * val, i, true);
  }

  /// Inverse tangent of multivector
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  atan(const Multivector<Scalar_T,LO,HI>& val)
  { return atan(val, elliptic(val), true); }

}
#endif  // _GLUCAT_CLIFFORD_ALGEBRA_IMP_H
