#ifndef _GLUCAT_CLIFFORD_ALGEBRA_IMP_H
#define _GLUCAT_CLIFFORD_ALGEBRA_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    clifford_algebra_imp.h : Implement common Clifford algebra functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2009 by Paul C. Leopardi
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
  { return !(lhs == scr); }

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
  { return numeric_traits<Scalar_T>::real(scalar(val)); }

  /// Imaginary part of scalar part
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  imag(const Multivector<Scalar_T,LO,HI>& val)
  { return numeric_traits<Scalar_T>::imag(scalar(val)); }

  /// Integer power of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  pow(const Multivector<Scalar_T,LO,HI>& lhs, int rhs)
  { 
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    typedef numeric_traits<Scalar_T> traits_t;

    multivector_t a;
    if (rhs < 0)
    {
      if (lhs == Scalar_T(0))
        return traits_t::NaN();
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
    typedef numeric_traits<Scalar_T> traits_t;

    if (lhs == Scalar_T(0))
    {
      const Scalar_T m = rhs.scalar();
      if (rhs == m)
        return (m < 0)
               ? traits_t::NaN()
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
    typedef numeric_traits<Scalar_T> traits_t;
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
      return traits_t::NaN();
  }

  /// Check that i is a valid complexifier for val
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  static
  void
  check_complex(const Multivector<Scalar_T,LO,HI>& val, 
                const Multivector<Scalar_T,LO,HI>& i, const bool prechecked = false)
  {
    if (!prechecked)
    {
      typedef Multivector<Scalar_T,LO,HI> multivector_t;
      typedef typename multivector_t::index_set_t index_set_t;
      typedef typename multivector_t::error_t error_t;

      const index_set_t i_frame = i.frame();
      // We need i to be a complexifier whose frame is large enough to represent val
      if (elliptic(i) != i ||
         (val.frame() | i_frame) != i_frame ||
          elliptic(val).frame().count() > i_frame.count())
        throw error_t("check_complex(val, i): i is not a valid complexifier for val");
    }
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
    typedef numeric_traits<Scalar_T> traits_t;

    if (X.isnan())
      return traits_t::NaN();

    const multivector_t& X2 = X*X;
    const multivector_t& X4 = X2*X2;
    const multivector_t& X6 = X4*X2;
    const multivector_t& N = a[0]+X2*a[2]+X4*a[4]+X6*(a[6]+X2*a[8]+X4*a[10]+X6*a[12]) +
                          X*(a[1]+X2*a[3]+X4*a[5]+X6*(a[7]+X2*a[9]+X4*a[11]+X6*a[13]));
    const multivector_t& D = b[0]+X2*b[2]+X4*b[4]+X6*(b[6]+X2*b[8]+X4*b[10]+X6*b[12]) +
                          X*(b[1]+X2*b[3]+X4*b[5]+X6*(b[7]+X2*b[9]+X4*b[11]+X6*b[13]));
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
    M = ((M + iM)/Scalar_T(2) + Scalar_T(1)) / Scalar_T(2);
    Y *= (iM + Scalar_T(1)) / Scalar_T(2);
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
    typedef numeric_traits<Scalar_T> traits_t;

    if (val == Scalar_T(0))
      return val;

    multivector_t M = val;
    multivector_t Y = val;
    static const Scalar_T tol = std::numeric_limits<Scalar_T>::epsilon() / Scalar_T(1024);
    static const Scalar_T tol2 = tol * tol;
    static const int sqrt_max_steps = Tune_P::sqrt_max_steps;
    Scalar_T norm_M_1 = norm(M - Scalar_T(1));

    for (int step = 0;
        step != sqrt_max_steps && norm_M_1 > tol2; 
        ++step)
    {
      if (Y.isnan())
        return traits_t::NaN();
      db_step(M, Y);
      norm_M_1 = norm(M - Scalar_T(1));
    }
    if (norm_M_1 > tol2)
      return traits_t::NaN();
    else
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

    // Pade approximation produced by Pade1(sqrt(1+x),x,13,13):
    // numer := 1+27/4*x+81/4*x^2+2277/64*x^3
    //          +10395/256*x^4+32319/1024*x^5+8721/512*x^6+26163/4096*x^7
    //          +53703/32768*x^8+36465/131072*x^9+3861/131072*x^10
    //          +7371/4194304*x^11+819/16777216*x^12+27/67108864*x^13;
    // denom := 1+25/4*x+69/4*x^2+1771/64*x^3
    //          +7315/256*x^4+20349/1024*x^5+4845/512*x^6+12597/4096*x^7
    //          +21879/32768*x^8+12155/131072*x^9+1001/131072*x^10
    //          +1365/4194304*x^11+91/16777216*x^12+1/67108864*x^13;

    static const Scalar_T a[] =
    {
          1.0,             27.0/4.0,         81.0/4.0,      2277.0/64.0,
      10395.0/256.0,    32319.0/1024.0,    8721.0/512.0,   26163.0/4096.0,
      53703.0/32768.0,  36465.0/131072.0,  3861.0/131072.0,
       7371.0/4194304.0,  819.0/16777216.0,  27.0/67108864.0
    };
    static const Scalar_T b[] =
    {
          1.0,             25.0/4.0,         69.0/4.0,      1771.0/64.0,
       7315.0/256.0,    20349.0/1024.0,    4845.0/512.0,   12597.0/4096.0, 
      21879.0/32768.0,  12155.0/131072.0,  1001.0/131072.0,
       1365.0/4194304.0,   91.0/16777216.0,   1.0/67108864.0
    };

    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    typedef numeric_traits<Scalar_T> traits_t;

    if (val == Scalar_T(0))
      return val;


    if (val.isnan())
      return traits_t::NaN();

    check_complex(val, i, prechecked);

    const Scalar_T realval = real(val);
    if (val == realval)
    {
      if (realval < Scalar_T(0))
        return i * traits_t::sqrt(-realval);
      else
        return traits_t::sqrt(realval);
    }
    static const Scalar_T sqrt_2 = traits_t::sqrt(Scalar_T(2));
    const multivector_t val2 = val*val;
    const Scalar_T realval2 = real(val2);
    if (val2 == realval2 && realval2 > 0)
      return sqrt(-i * val) * (i + Scalar_T(1)) / sqrt_2;

    // Scale val towards abs(A) == 1 or towards A == 1 as appropriate
    const Scalar_T max_norm = Scalar_T(1.0/9.0);
    const Scalar_T scale = 
      (realval != Scalar_T(0) && norm(val/realval - Scalar_T(1)) < Scalar_T(1))
        ? realval 
        : (realval < Scalar_T(0)) 
          ? -abs(val) 
          : abs(val);
    const Scalar_T sqrt_scale = traits_t::sqrt(traits_t::abs(scale));
    if (traits_t::isNaN_or_isInf(sqrt_scale))
      return traits_t::NaN();

    multivector_t rescale = sqrt_scale;
    if (scale < Scalar_T(0))
      rescale = i * sqrt_scale;

    const multivector_t& unitval = val / scale;
    const multivector_t& scaled_result = 
      (norm(unitval - Scalar_T(1)) < max_norm)
        // Pade' approximation of square root
        ? pade_approx(a, b, unitval - Scalar_T(1))
        // Product form of Denman-Beavers square root iteration
        : db_sqrt(unitval);
    if (scaled_result.isnan())
    {
      const multivector_t& mi_unitval = -i * unitval;
      const multivector_t& scaled_mi_result = 
        (norm(mi_unitval - Scalar_T(1)) < max_norm) 
          ? pade_approx(a, b, mi_unitval - Scalar_T(1))
          : db_sqrt(mi_unitval);
      if (scaled_mi_result.isnan())
        return traits_t::NaN();
      else
        return scaled_mi_result * rescale * (i + Scalar_T(1)) / sqrt_2;
    }
    else 
      return scaled_result * rescale;
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
    typedef numeric_traits<Scalar_T> traits_t;

    if (val == Scalar_T(0))
      return Scalar_T(1);
    if (val.isnan())
      return traits_t::NaN();

    const Scalar_T scalar_val = scalar(val);
    const Scalar_T scalar_exp = traits_t::exp(scalar_val);
    if (traits_t::isNaN_or_isInf(scalar_exp))
      return traits_t::NaN();
    if (val == scalar_val)
      return scalar_exp;

    multivector_t A = val - scalar_val;
    const Scalar_T pure_scale = A.norm();

    if (traits_t::isNaN_or_isInf(pure_scale))
      return traits_t::NaN();
    if (pure_scale == Scalar_T(0))
      return scalar_exp;

    const int ilog2_scale = 
      std::max(0, traits_t::to_int(ceil((log2(pure_scale) + Scalar_T(A.frame().count()))/Scalar_T(2)) - 3));
    const Scalar_T i_scale = traits_t::pow(Scalar_T(2), ilog2_scale);
    if (traits_t::isNaN_or_isInf(i_scale))
      return traits_t::NaN();

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

  /// Pade' approximation of log
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

    // Pade approximation produced by Pade1(log(1-x),x,13,13):
    // numer := -x+6*x^2-4741/300*x^3
    //          +1441/60*x^4-107091/4600*x^5+8638/575*x^6-263111/40250*x^7
    //          +153081/80500*x^8-395243/1101240*x^9+28549/688275*x^10
    //          -605453/228813200*x^11+785633/10296594000*x^12
    //          -1145993/1873980108000*x^13;
    // denom := 1-13/2*x+468/25*x^2-1573/50*x^3
    //          +1573/46*x^4-11583/460*x^5+10296/805*x^6-2574/575*x^7
    //          +11583/10925*x^8-143/874*x^9+572/37145*x^10
    //          -117/148580*x^11+13/742900*x^12-1/10400600*x^13;

    static const Scalar_T a[] =
    {
            0.0,                  -1.0,                      6.0,            -4741.0/300.0,
         1441.0/60.0,        -107091.0/4600.0,            8638.0/575.0,    -263111.0/40250.0, 
       153081.0/80500.0,     -395243.0/1101240.0,        28549.0/688275.0, 
      -605453.0/228813200.0,  785633.0/10296594000.0, -1145993.0/1873980108000.0
    };
    static const Scalar_T b[] =
    {
         1.0,                    -13.0/2.0,                468.0/25.0,       -1573.0/50.0, 
      1573.0/46.0,            -11583.0/460.0,            10296.0/805.0,      -2574.0/575.0,
     11583.0/10925.0,           -143.0/874.0,              572.0/37145.0,
      -117.0/148580.0,            13.0/742900.0,            -1.0/10400600.0
    };
    typedef numeric_traits<Scalar_T> traits_t;

    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();
    else
      return pade_approx(a, b, -val + Scalar_T(1));
  }

  /// Incomplete square root cascade and Pade' approximation of log
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  static
  const Multivector<Scalar_T,LO,HI>
  cascade_log(const Multivector<Scalar_T,LO,HI>& val)
  {
    // Reference: [CHKL]
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    typedef numeric_traits<Scalar_T> traits_t;

    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();

    multivector_t Y = val;
    multivector_t E = Scalar_T(0);
    static const Scalar_T epsilon = std::numeric_limits<Scalar_T>::epsilon();
    static const Scalar_T max_inner_norm = traits_t::pow(epsilon, 4);
    static const Scalar_T max_outer_norm = Scalar_T(1.0/9.0);
    Scalar_T norm_Y_1;
    int outer_step;
    Scalar_T pow_2_outer_step = Scalar_T(1);
    Scalar_T pow_4_outer_step = Scalar_T(1);
    for (outer_step = 0, norm_Y_1 = norm(Y - Scalar_T(1));
        outer_step != Tune_P::log_max_outer_steps && norm_Y_1 * pow_2_outer_step > max_outer_norm;
        ++outer_step,    norm_Y_1 = norm(Y - Scalar_T(1)))
    {
      if (Y == Scalar_T(0) || Y.isnan())
        return traits_t::NaN();

      // Incomplete product form of Denman-Beavers square root iteration
      multivector_t M = Y;
      for (int 
          inner_step = 0;
          inner_step != Tune_P::log_max_inner_steps && 
            norm(M - Scalar_T(1)) * pow_4_outer_step > max_inner_norm;
          ++inner_step)
        db_step(M, Y);

      E += (M - Scalar_T(1)) * pow_2_outer_step;
      pow_2_outer_step *= Scalar_T(2);
      pow_4_outer_step *= Scalar_T(4);
    }
    if (outer_step == Tune_P::log_max_outer_steps && norm_Y_1 * pow_2_outer_step > max_outer_norm)
      return traits_t::NaN();
    else
      return pade_log(Y) * pow_2_outer_step - E;
  }

  /// Natural logarithm of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  log(const Multivector<Scalar_T,LO,HI>& val, const Multivector<Scalar_T,LO,HI>& i, bool prechecked)
  {
    // Scaled incomplete square root cascade and scaled Pade' approximation of log
    // Reference: [CHKL]

    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    typedef numeric_traits<Scalar_T> traits_t;

    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();

    const Scalar_T realval = real(val);
    if (val == realval)
    {
      if (realval < Scalar_T(0))
      {
        check_complex(val, i, prechecked);
        return i * traits_t::pi() + traits_t::log(-realval);
      }
      else
        return traits_t::log(realval);
    } 
    const multivector_t val2 = val*val;
    const Scalar_T realval2 = real(val2);
    if (val2 == realval2 && realval2 > 0)
    {
      check_complex(val, i, prechecked);
      return log(-i * val) + i * traits_t::pi()/Scalar_T(2);
    } 
    // Scale val towards abs(A) == 1 or towards A == 1 as appropriate
    const Scalar_T max_norm = Scalar_T(1.0/9.0);
    const Scalar_T scale = 
      (realval != Scalar_T(0) && norm(val/realval - Scalar_T(1)) < max_norm)
        ? realval 
        : (realval < Scalar_T(0)) 
          ? -abs(val) 
          : abs(val);
    if (scale == Scalar_T(0))
      return traits_t::NaN();

    const Scalar_T log_scale = traits_t::log(traits_t::abs(scale));
    multivector_t rescale = log_scale;
    if (scale < Scalar_T(0))
    {
      check_complex(val, i, prechecked);
      rescale = i * traits_t::pi() + log_scale;
    }
    const multivector_t unitval = val/scale;
    if (inv(unitval).isnan())
      return traits_t::NaN();
    multivector_t scaled_result = cascade_log(unitval);
    if (scaled_result.isnan())
      return traits_t::NaN();
    else 
      return scaled_result + rescale;
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
    typedef numeric_traits<Scalar_T> traits_t;
    return val.isnan()
         ? traits_t::NaN()
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
    typedef numeric_traits<Scalar_T> traits_t;
    check_complex(val, i, prechecked);
    return val.isnan()
         ? traits_t::NaN()
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
    typedef numeric_traits<Scalar_T> traits_t;

    if (val.isnan())
      return traits_t::NaN();

    check_complex(val, i, prechecked);
    const Scalar_T& s = scalar(val);
    static const Scalar_T& twopi = Scalar_T(2) * traits_t::pi();
    const multivector_t& z = i *
      (val - s + traits_t::fmod(traits_t::real(s), twopi) + traits_t::imag(s));
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
    typedef numeric_traits<Scalar_T> traits_t;
    check_complex(val, i, prechecked);
    return val.isnan()
         ? traits_t::NaN()
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
    typedef numeric_traits<Scalar_T> traits_t;
    return val.isnan()
           ? traits_t::NaN()
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
    typedef numeric_traits<Scalar_T> traits_t;
    check_complex(val, i, prechecked);
    return val.isnan()
         ? traits_t::NaN()
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
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    typedef numeric_traits<Scalar_T> traits_t;

    if (val.isnan())
      return traits_t::NaN();

    check_complex(val, i, prechecked);
    const Scalar_T& s = scalar(val);
    static const Scalar_T& twopi = Scalar_T(2) * traits_t::pi();
    const multivector_t& z = i *
      (val - s + traits_t::fmod(traits_t::real(s), twopi) + traits_t::imag(s));
    return i * (exp(-z)-exp(z)) / Scalar_T(2);
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
    typedef numeric_traits<Scalar_T> traits_t;
    check_complex(val, i, prechecked);
    return val.isnan()
         ? traits_t::NaN()
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
    typedef numeric_traits<Scalar_T> traits_t;
    return val.isnan()
         ? traits_t::NaN()
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
    typedef numeric_traits<Scalar_T> traits_t;
    check_complex(val, i, prechecked);
    return val.isnan()
        ? traits_t::NaN()
        : (norm(val + Scalar_T(1)) > norm(val - Scalar_T(1)))
          ? (log(val + Scalar_T(1), i, true) - log(-val + Scalar_T(1), i, true)) / Scalar_T(2)
          : log((val + Scalar_T(1)) / (-val + Scalar_T(1)), i, true) / Scalar_T(2); 
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
    typedef numeric_traits<Scalar_T> traits_t;
    check_complex(val, i, prechecked);
    return val.isnan()
         ? traits_t::NaN()
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
    typedef numeric_traits<Scalar_T> traits_t;
    check_complex(val, i, prechecked);
    return val.isnan()
          ? traits_t::NaN()
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
