#ifndef _GLUCAT_CLIFFORD_ALGEBRA_H
#define _GLUCAT_CLIFFORD_ALGEBRA_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    clifford_algebra.h : Declare the operations of a Clifford algebra
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

namespace glucat
{
  /// clifford_algebra<> declares the operations of a Clifford algebra
  template< typename Scalar_T, typename Index_Set_T, typename Multivector_T>
  class clifford_algebra
  {
  public:
    typedef Scalar_T                      scalar_t;
    typedef Index_Set_T                   index_set_t;
    typedef Multivector_T                 multivector_t;
    typedef std::pair< const index_set_t, Scalar_T > pair_t;
    typedef std::vector<Scalar_T>         vector_t;
    static  char *      classname();

    virtual ~clifford_algebra() {};

  // clifford_algebra operations
    /// Test for equality of multivectors
    virtual bool                operator==  (const multivector_t& val) const = 0;
    /// Test for equality of multivector and scalar
    virtual bool                operator==  (const Scalar_T& scr) const =0;
    /// Geometric sum
    virtual multivector_t&      operator+=  (const multivector_t& rhs) =0;
    /// Geometric sum of multivector and scalar
    virtual multivector_t&      operator+=  (const Scalar_T& scr) =0;
    /// Geometric difference
    virtual multivector_t&      operator-=  (const multivector_t& rhs) =0;
    /// Unary -
    virtual const multivector_t operator-   () const =0;
    /// Product of multivector and scalar
    virtual multivector_t&      operator*=  (const Scalar_T& scr) =0;
    /// Geometric product
    virtual multivector_t&      operator*=  (const multivector_t& rhs) =0;
    /// Contraction
    virtual multivector_t&      operator%=  (const multivector_t& rhs) =0;
    /// Inner product
    virtual multivector_t&      operator&=  (const multivector_t& rhs) =0;
    /// Outer product
    virtual multivector_t&      operator^=  (const multivector_t& rhs) =0;
    /// Quotient of multivector and scalar
    virtual multivector_t&      operator/=  (const Scalar_T& scr) =0;
    /// Geometric quotient
    virtual multivector_t&      operator/=  (const multivector_t& rhs) =0;
    /// Geometric multiplicative inverse
    virtual const multivector_t inv         () const =0;
    /// *this to the m
    virtual const multivector_t pow         (int m) const =0;
    /// Outer product power
    virtual const multivector_t outer_pow   (int m) const =0;
    /// Subscripting: map from index set to scalar coordinate
    virtual Scalar_T            operator[]  (const index_set_t& ist) const =0;
    /// Pure grade-vector part
    virtual const multivector_t operator()  (index_t grade) const =0;
    /// Even part of multivector, sum of even grade parts
    virtual const multivector_t even()      const =0;
    /// Vector part of multivector, as a vector_t with respect to frame()
    virtual const vector_t      vector_part () const = 0;
		/// Main involution, each {i} is replaced by -{i} in each term, eg. {1}*{2} -> (-{2})*(-{1})
    virtual const multivector_t involute()  const =0;
    /// Reversion, eg. {1}*{2} -> {2}*{1}
    virtual const multivector_t reverse()   const =0;
    /// Conjugation, reverse o involute == involute o reverse
    virtual const multivector_t conj()      const =0;
    /// Scalar_T norm == sum of norm of coordinates
    virtual Scalar_T            norm()      const =0;
    /// Scalar_T quadratic form == (rev(x)*x)(0)
    virtual Scalar_T            quad()      const =0;
    virtual Scalar_T            max_abs()   const =0;
    virtual const multivector_t truncated
                  (const Scalar_T& limit = Scalar_T(DEFAULT_TRUNCATION)) const =0;
    /// Check if a multivector contains any IEEE NaN values
    virtual bool                isnan()     const =0;
    /// Subalgebra generated by all generators of terms of given multivector
    virtual const index_set_t   frame()     const =0;
    /// Write formatted multivector to output
    virtual void                write       (char* msg="") const =0;
    /// Write formatted multivector to file
    virtual void                write       (std::ofstream& ofile, char* msg="") const =0;
  };

#ifndef _GLUCAT_CLIFFORD_ALGEBRA_OPERATIONS
#define _GLUCAT_CLIFFORD_ALGEBRA_OPERATIONS \
    bool                operator==  (const multivector_t& val) const; \
    bool                operator==  (const Scalar_T& scr) const;      \
    multivector_t&      operator+=  (const multivector_t& rhs);       \
    multivector_t&      operator+=  (const Scalar_T& scr);            \
    multivector_t&      operator-=  (const multivector_t& rhs);       \
    const multivector_t operator-   () const;                         \
    multivector_t&      operator*=  (const Scalar_T& scr);            \
    multivector_t&      operator*=  (const multivector_t& rhs);       \
    multivector_t&      operator%=  (const multivector_t& rhs);       \
    multivector_t&      operator&=  (const multivector_t& rhs);       \
    multivector_t&      operator^=  (const multivector_t& rhs);       \
    multivector_t&      operator/=  (const Scalar_T& scr);            \
    multivector_t&      operator/=  (const multivector_t& rhs);       \
    const multivector_t inv         () const;                         \
    const multivector_t pow         (int m) const;                    \
    const multivector_t outer_pow   (int m) const;                    \
    Scalar_T            operator[]  (const index_set_t& ist) const;   \
    const multivector_t operator()  (index_t grade) const;            \
    const multivector_t even()      const;                            \
    const vector_t      vector_part() const;                          \
    const multivector_t involute()  const;                            \
    const multivector_t reverse()   const;                            \
    const multivector_t conj()      const;                            \
    Scalar_T            norm()      const;                            \
    Scalar_T            quad()      const;                            \
    const index_set_t   frame()     const;                            \
    Scalar_T            max_abs()   const;                            \
    const multivector_t truncated                                     \
          (const Scalar_T& limit = Scalar_T(DEFAULT_TRUNCATION)) const; \
    bool                isnan       () const;                         \
    void                write       (char* msg="") const;             \
    void                write       (std::ofstream& ofile, char* msg="") const;
#endif // _GLUCAT_CLIFFORD_ALGEBRA_OPERATIONS

  // Test for inequality of multivectors
  template
	<
		template<typename, const index_t, const index_t> class Multivector,
		template<typename, const index_t, const index_t> class RHS,
		typename Scalar_T, const index_t LO, const index_t HI
	>
  bool
  operator!= (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs);

  // Test for inequality of multivector and scalar
  template< template<typename, const index_t, const index_t> class Multivector, typename Scalar_T, const index_t LO, const index_t HI >
  bool
  operator!= (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr);

  // Test for inequality of scalar and multivector
  template< template<typename, const index_t, const index_t> class Multivector, typename Scalar_T, const index_t LO, const index_t HI >
  bool
  operator!= (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs);

  /// Geometric sum of multivector and scalar
  template< template<typename, const index_t, const index_t> class Multivector, typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  operator+ (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr);

  /// Geometric sum of scalar and multivector
  template< template<typename, const index_t, const index_t> class Multivector, typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  operator+ (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs);

  /// Geometric sum
  template
	<
		template<typename, const index_t, const index_t> class Multivector,
		template<typename, const index_t, const index_t> class RHS,
		typename Scalar_T, const index_t LO, const index_t HI
	>
  const Multivector<Scalar_T,LO,HI>
  operator+ (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs);

  /// Geometric difference of multivector and scalar
  template< template<typename, const index_t, const index_t> class Multivector, typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  operator- (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr);

  /// Geometric difference of scalar and multivector
  template< template<typename, const index_t, const index_t> class Multivector, typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  operator- (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs);

  /// Geometric difference
  template
	<
		template<typename, const index_t, const index_t> class Multivector,
		template<typename, const index_t, const index_t> class RHS,
		typename Scalar_T, const index_t LO, const index_t HI
	>
  const Multivector<Scalar_T,LO,HI>
  operator- (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs);

  /// Product of multivector and scalar
  template< template<typename, const index_t, const index_t> class Multivector, typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  operator* (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr);

  /// Product of scalar and multivector
  template< template<typename, const index_t, const index_t> class Multivector, typename Scalar_T, const index_t LO, const index_t HI >
  const Multivector<Scalar_T,LO,HI>
  operator* (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs);

  /// Geometric product
  template
	<
		template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
		typename Scalar_T, const index_t LO, const index_t HI
	>
  const Multivector<Scalar_T,LO,HI>
  operator* (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs);

  /// Left contraction
  template
	<
		template<typename, const index_t, const index_t> class Multivector,
		template<typename, const index_t, const index_t> class RHS,
		typename Scalar_T, const index_t LO, const index_t HI
	>
  const Multivector<Scalar_T,LO,HI>
  operator% (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs);

  /// Inner product
  template
	<
		template<typename, const index_t, const index_t> class Multivector,
		template<typename, const index_t, const index_t> class RHS,
		typename Scalar_T, const index_t LO, const index_t HI
	>
  const Multivector<Scalar_T,LO,HI>
  operator& (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs);

  /// Outer product
  template
	<
		template<typename, const index_t, const index_t> class Multivector,
		template<typename, const index_t, const index_t> class RHS,
		typename Scalar_T, const index_t LO, const index_t HI
	>
  const Multivector<Scalar_T,LO,HI>
  operator^ (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs);

  /// Hestenes scalar product
  template
	<
		template<typename, const index_t, const index_t> class Multivector,
		template<typename, const index_t, const index_t> class RHS,
		typename Scalar_T, const index_t LO, const index_t HI
	>
  Scalar_T
  star(const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs);

  /// Quotient of multivector and scalar
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  operator/ (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr);

  /// Quotient of scalar and multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  operator/ (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs);

  /// Geometric quotient
  template
	<
		template<typename, const index_t, const index_t> class Multivector,
		template<typename, const index_t, const index_t> class RHS,
		typename Scalar_T, const index_t LO, const index_t HI
	>
  const Multivector<Scalar_T,LO,HI>
  operator/ (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs);

  /// Geometric multiplicative inverse
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  inv(const Multivector<Scalar_T,LO,HI>& val);

  /// Main involution, each {i} is replaced by -{i} in each term, eg. {1}*{2} -> (-{2})*(-{1})
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  involute(const Multivector<Scalar_T,LO,HI>& val);

  /// Reversion, eg. {1}*{2} -> {2}*{1}
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  reverse(const Multivector<Scalar_T,LO,HI>& val);

  /// Conjugation, rev o invo == invo o rev
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  conj(const Multivector<Scalar_T,LO,HI>& val);

  /// Scalar_T quadratic form == (rev(x)*x)(0)
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  Scalar_T
  quad(const Multivector<Scalar_T,LO,HI>& val);

  /// Scalar part
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  Scalar_T
  scalar(const Multivector<Scalar_T,LO,HI>& val);

  /// Pure part
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  pure(const Multivector<Scalar_T,LO,HI>& val);

  /// Even part
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  even(const Multivector<Scalar_T,LO,HI>& val);

  /// Vector part of multivector, as a vector_t with respect to frame()
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const std::vector<Scalar_T>
  vector_part(const Multivector<Scalar_T,LO,HI>& val);

  /// Absolute value == sqrt(norm)
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  Scalar_T
  abs(const Multivector<Scalar_T,LO,HI>& val);

  /// Scalar_T norm == sum of norm of coordinates
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  Scalar_T
  norm(const Multivector<Scalar_T,LO,HI>& val);

  /// Real part of scalar part
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  Scalar_T
  real(const Multivector<Scalar_T,LO,HI>& val);

  /// Imaginary part of scalar part
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  Scalar_T
  imag(const Multivector<Scalar_T,LO,HI>& val);

  /// Integer power of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  pow(const Multivector<Scalar_T,LO,HI>& lhs, int rhs);

  /// Multivector power of multivector
  template
	<
		template<typename, const index_t, const index_t> class Multivector,
		template<typename, const index_t, const index_t> class RHS,
		typename Scalar_T, const index_t LO, const index_t HI
	>
  const Multivector<Scalar_T,LO,HI>
  pow(const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs);

  /// Square root of -1 which commutes with all members of the frame of the given multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  elliptic(const Multivector<Scalar_T,LO,HI>& val);

  /// Square root of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  sqrt(const Multivector<Scalar_T,LO,HI>& val);

  // Transcendental functions
  /// Exponential of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  exp(const Multivector<Scalar_T,LO,HI>& val);

  /// Natural logarithm of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  log(const Multivector<Scalar_T,LO,HI>& val);

  /// Cosine of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  cos(const Multivector<Scalar_T,LO,HI>& val);

  /// Inverse cosine of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  acos(const Multivector<Scalar_T,LO,HI>& val);

  /// Hyperbolic cosine of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  cosh(const Multivector<Scalar_T,LO,HI>& val);

  /// Inverse hyperbolic cosine of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  acosh(const Multivector<Scalar_T,LO,HI>& val);

  /// Sine of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  sin(const Multivector<Scalar_T,LO,HI>& val);

  /// Inverse sine of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  asin(const Multivector<Scalar_T,LO,HI>& val);

  /// Hyperbolic sine of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  sinh(const Multivector<Scalar_T,LO,HI>& val);

  /// Inverse hyperbolic sine of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  asinh(const Multivector<Scalar_T,LO,HI>& val);

  /// Tangent of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  tan(const Multivector<Scalar_T,LO,HI>& val);

  /// Inverse tangent of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  atan(const Multivector<Scalar_T,LO,HI>& val);

  /// Hyperbolic tangent of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  tanh(const Multivector<Scalar_T,LO,HI>& val);

  /// Inverse hyperbolic tangent of multivector
  template<	template<typename, const index_t, const index_t> class Multivector,
						typename Scalar_T, const index_t LO, const index_t HI	>
  const Multivector<Scalar_T,LO,HI>
  atanh(const Multivector<Scalar_T,LO,HI>& val);
}
#endif  // _GLUCAT_CLIFFORD_ALGEBRA_H
