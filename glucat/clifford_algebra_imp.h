#ifndef _GLUCAT_CLIFFORD_ALGEBRA_IMP_H
#define _GLUCAT_CLIFFORD_ALGEBRA_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    clifford_algebra_imp.h : Implement common Clifford algebra functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001 by Paul C. Leopardi
    email                : leopardi@bigpond.net.au
 ***************************************************************************
 *   This library is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2.1 of the  *
 *   License, or (at your option) any later version.                       *
 *   See http://www.fsf.org/copyleft/lesser.html for details               *
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

namespace glucat
{
  template< typename Scalar_T, typename Index_Set_T, typename Multivector_T>
  const std::string
  clifford_algebra<Scalar_T,Index_Set_T,Multivector_T>::
  classname()
  { return "clifford_algebra"; }

  // Test for inequality of multivectors
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


  // Test for inequality of multivector and scalar
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  bool
  operator!= (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr)
  { return lhs.frame().count() != 0 || scalar(lhs) != scr; }

  // Test for inequality of scalar and multivector
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
    Multivector<Scalar_T,LO,HI> result = lhs;
    return result *= rhs;
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
    Multivector<Scalar_T,LO,HI> result = lhs;
    return result %= rhs;
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
    Multivector<Scalar_T,LO,HI> result = lhs;
    return result &= rhs;
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
    Multivector<Scalar_T,LO,HI> result = lhs;
    return result ^= rhs;
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
  star(const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs)
  { return scalar(lhs * rhs); }

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
    Multivector<Scalar_T,LO,HI> result = lhs;
    return result /= rhs;
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
  { return val[index_set<LO,HI>()]; }

  /// Pure part
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  pure(const Multivector<Scalar_T,LO,HI>& val)
  { return val - scalar(val); }

  /// Even part
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  even(const Multivector<Scalar_T,LO,HI>& val)
  { return val.even(); }

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
  { return std::sqrt(val.norm()); }

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
  inline
  const Multivector<Scalar_T,LO,HI>
  pow(const Multivector<Scalar_T,LO,HI>& lhs, int rhs)
  { return lhs.pow(rhs); }

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
  { return exp(log(lhs) * rhs); }

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
    for (index_t k = 0; k != incp[bott]; k++)
      for (index_t idx = 1; idx != HI+1; ++idx)
        if (!frm[idx])
        {
          frm.set(idx);
          break;
        }
    for (index_t k = 0; k != incq[bott]; k++)
      for (index_t idx = -1; idx != LO-1; --idx)
        if (!frm[idx])
        {
          frm.set(idx);
          break;
        }
    index_t new_bott = pos_mod(frm.count_pos() - frm.count_neg(), 4);

    if ((incp[new_bott] == 0) && (incq[new_bott] == 0))
      return multivector_t(frm, Scalar_T(1));
    else
      return Scalar_T(std::log(0.0)); // This actually returns -Inf;
  }

  /// Pade' approximation
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  pade_approx(const Scalar_T a[], const Scalar_T b[], const Multivector<Scalar_T,LO,HI>& A)
  {
    // Pade' approximation
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [GL], Section 11.3, p572-576.

    typedef Multivector<Scalar_T,LO,HI> multivector_t;

    const multivector_t& A2 = A*A;
    const multivector_t& A4 = A2*A2;
    const multivector_t& N = a[0]+A2*a[2]+(a[4]+A2*a[6]+A4*a[8])*A4 +
                          A*(a[1]+A2*a[3]+(a[5]+A2*a[7])*A4);
    const multivector_t& D = b[0]+A2*b[2]+(b[4]+A2*b[6]+A4*b[8])*A4 +
                          A*(b[1]+A2*b[3]+(b[5]+A2*b[7])*A4);
    return N / D;
  }

  /// Single step of product form of Denman-Beavers square root iteration
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
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
  const Multivector<Scalar_T,LO,HI>
  db_sqrt(const Multivector<Scalar_T,LO,HI>& val)
  {
    // Reference: [CHKL]
    typedef Multivector<Scalar_T,LO,HI> multivector_t;

    if (val == 0)
      return val;
    multivector_t M = val;
    multivector_t Y = val;
    int step = 0;
    for (; step != Tune_P::sqrt_max_steps && norm(Scalar_T(1) - M) > 0; ++step)
      db_step(M, Y);
    if (step == Tune_P::sqrt_max_steps && norm(Scalar_T(1) - M) > 0)
      std::cerr << "Warning: sqrt iteration did not converge. norm = "
           << norm(Scalar_T(1) - M) << std::endl;
    return Y;
  }

  /// Square root of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  sqrt(const Multivector<Scalar_T,LO,HI>& val)
  {
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [GL], Section 11.3, p572-576.
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    static const Scalar_T a[] =
    {
       1.0,          17.0/4.0,      119.0/16.0,
     221.0/32.0,    935.0/256.0,    561.0/512.0,
     357.0/2048.0,   51.0/4096.0,    17.0/65536.0
    };
    static const Scalar_T b[] =
    {
       1.0,         15.0/4.0,       91.0/16.0,
     143.0/32.0,   495.0/256.0,    231.0/512.0,
     105.0/2048.0,   9.0/4096.0,     1.0/65536.0
    };

    if (val == 0)
      return val;
    const Scalar_T realval = real(val);
    if (val == realval && realval < Scalar_T(0))
      return Scalar_T(std::sqrt(-realval)) * elliptic(val);

    // Scale val towards circle abs(val) == 1
    const Scalar_T scale = abs(val);
    const multivector_t& unitval = val / scale;
    if (norm(unitval - Scalar_T(1)) < Scalar_T(1))
      // Pade' approximation of square root
      return pade_approx(a, b, unitval - Scalar_T(1)) * Scalar_T(std::sqrt(scale));
    else
      // Product form of Denman-Beavers square root iteration
      return db_sqrt(unitval) * Scalar_T(std::sqrt(scale));
  }

  /// Exponential of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  exp(const Multivector<Scalar_T,LO,HI>& val)
  {
    // Scaling and squaring Pade' approximation of matrix exponential
    // Reference: [GL], Section 11.3, p572-576.
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    multivector_t A = pure(val);
    const int j = std::max(0, 1 + int(floor(log2(A.max_abs()))));
    A /= Scalar_T(std::pow(Scalar_T(2),j));
    multivector_t result;
    {
      const int q = 8;
      static Scalar_T c[q+1];
      if (c[0] != Scalar_T(1))
      {
        c[0] = Scalar_T(1);
        for (int k = 0; k != q; ++k)
          c[k+1] = c[k]*(q-k) / ((2*q-k)*(k+1));
      }
      const multivector_t& A2 = A*A;
      const multivector_t& A4 = A2*A2;
      const multivector_t& U =     c[0]+A2*c[2]+(c[4]+A2*c[6]+A4*c[8])*A4;
      const multivector_t& AV = A*(c[1]+A2*c[3]+(c[5]+A2*c[7])*A4);
      result = (U+AV) / (U-AV);
    }
    for (int k = 0; k != j && !result.isnan(); ++k)
      result *= result;
    return Scalar_T(std::exp(scalar(val))) * result;
  }

  /// Scaled Pade' approximation of log
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  pade_log(const Multivector<Scalar_T,LO,HI>& val)
  {
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [CHKL]
    // Reference: [GL], Section 11.3, p572-576.
    static const Scalar_T a[] =
    {
        0.0,        -1.0,          7.0/2.0,
      -73.0/15.0,   41.0/12.0,  -743.0/585.0,
       31.0/130.0, -37.0/1925.0,  43.0/101810.0
    };
    static const Scalar_T b[] =
    {
        1.0,        -4.0,        98.0/15.0,
      -28.0/5.0,    35.0/13.0,  -28.0/39.0,
       14.0/143.0, -4.0/715.0,    1.0/12870.0
    };
    typedef Multivector<Scalar_T,LO,HI> multivector_t;

    if (val == 0)
      return Scalar_T(std::log(0.0));
    // Scale val towards circle abs(A) == 1
    Scalar_T scale = abs(val);
    return pade_approx(a, b, Scalar_T(1) - val / scale) + Scalar_T(std::log(scale)) ;
  }

  /// Incomplete square root cascade and scaled Pade' approximation of log
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  cascade_log(const Multivector<Scalar_T,LO,HI>& val)
  {
    // Reference: [CHKL]
    typedef Multivector<Scalar_T,LO,HI> multivector_t;

    if (val == 0)
      return Scalar_T(std::log(0.0));
    multivector_t Y = val;
    multivector_t E = Scalar_T(0);
    int outer_step;
    for (outer_step = 0;
        outer_step != Tune_P::log_max_outer_steps &&
        norm(Scalar_T(1) - Y) > Scalar_T(0.5);
        ++outer_step)
    {
      if (Y == 0)
        return Scalar_T(std::log(0.0));
      // Incomplete product form of Denman-Beavers square root iteration
      multivector_t M = Y;
      for (int inner_step = 0;
          inner_step != Tune_P::log_max_inner_steps && norm(Scalar_T(1) - M) > 0;
          ++inner_step)
        db_step(M, Y);
      E += E + Scalar_T(1) - M;
    }
    if (outer_step == Tune_P::log_max_outer_steps && norm(Scalar_T(1) - Y) > Scalar_T(0.5))
      std::cerr << "Warning: log iteration did not converge. norm = "
                << norm(Scalar_T(1) - Y) << std::endl;
    return pade_log(Y) * Scalar_T(std::pow(Scalar_T(2), outer_step)) + E;
  }

  /// Natural logarithm of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  log(const Multivector<Scalar_T,LO,HI>& val)
  {
    // Scaled incomplete square root cascade and scaled Pade' approximation of log
    // Reference: [CHKL]
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    if (val == 0)
      return Scalar_T(std::log(0.0));
    const Scalar_T realval = real(val);
    if (val == realval && realval < Scalar_T(0))
      return Scalar_T(std::log(-realval)) + elliptic(val) * Scalar_T(l_pi);

    // Scale val towards abs(A) == 1
    const Scalar_T scale = abs(val);
    const Scalar_T logscale = std::log(scale);
    if (val == scale)
      return logscale;
    else
      return cascade_log(val / scale) + logscale;
  }

  /// Cosine of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  cos(const Multivector<Scalar_T,LO,HI>& val)
  {
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    const Scalar_T s = scalar(val);
    const Scalar_T twopi = Scalar_T(2) * Scalar_T(l_pi);
    const multivector_t& i = elliptic(val);
    const multivector_t& z = i *
      (val - s + Scalar_T(std::fmod(Scalar_T(std::real(s)), twopi)) +
                 Scalar_T(std::imag(s)));
    return (exp(z)+exp(-z)) / Scalar_T(2);
  }

  /// Inverse cosine of multivector
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  acos(const Multivector<Scalar_T,LO,HI>& val)
  { return elliptic(val)*acosh(val); }

  /// Hyperbolic cosine of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  cosh(const Multivector<Scalar_T,LO,HI>& val)
  { return (exp(val)+exp(-val)) / Scalar_T(2); }

  /// Inverse hyperbolic cosine of multivector
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  acosh(const Multivector<Scalar_T,LO,HI>& val)
  { return log(val + sqrt(val*val - Scalar_T(1))); }

  /// Sine of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  sin(const Multivector<Scalar_T,LO,HI>& val)
  {
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    const Scalar_T s = scalar(val);
    const Scalar_T twopi = Scalar_T(2) * Scalar_T(l_pi);
    const multivector_t& i = elliptic(val);
    const multivector_t& z = i *
      (val - s + Scalar_T(std::fmod(Scalar_T(std::real(s)), twopi)) +
                 Scalar_T(std::imag(s)));
    return (exp(z)-exp(-z)) / (i*Scalar_T(2));
  }

  /// Inverse sine of multivector
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  asin(const Multivector<Scalar_T,LO,HI>& val)
  {
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    const multivector_t& i = elliptic(val);
    return -i*asinh(i*val);
  }

  /// Hyperbolic sine of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  sinh(const Multivector<Scalar_T,LO,HI>& val)
  { return (exp(val)-exp(-val)) / Scalar_T(2); }

  /// Inverse hyperbolic sine of multivector
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  asinh(const Multivector<Scalar_T,LO,HI>& val)
  { return log(val + sqrt(val*val + Scalar_T(1))); }

  /// Tangent of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  tan(const Multivector<Scalar_T,LO,HI>& val)
  { return sin(val) / cos(val); }

  /// Inverse tangent of multivector
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  atan(const Multivector<Scalar_T,LO,HI>& val)
  {
    typedef Multivector<Scalar_T,LO,HI> multivector_t;
    const multivector_t& i = elliptic(val);
    return -i*atanh(i*val);
  }

  /// Hyperbolic tangent of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  tanh(const Multivector<Scalar_T,LO,HI>& val)
  { return sinh(val) / cosh(val); }

  /// Inverse hyperbolic tangent of multivector
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const Multivector<Scalar_T,LO,HI>
  atanh(const Multivector<Scalar_T,LO,HI>& val)
  { return log((Scalar_T(1) + val) / (Scalar_T(1) - val)) / Scalar_T(2); }
}
#endif  // _GLUCAT_CLIFFORD_ALGEBRA_IMP_H
