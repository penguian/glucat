#ifndef _GLUCAT_SCALAR_H
#define _GLUCAT_SCALAR_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    scalar.h : Define functions for scalar_t
                             -------------------
    begin                : 2001-12-20
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
 "Clifford algebras with numeric and symbolic computations, Birkhauser, 1996."
 ***************************************************************************
 See also Arvind Raja's original header comments and references in glucat.h
 ***************************************************************************/

namespace glucat
{
  /// Extra traits which extend numeric limits
  // Reference: [AA], 2.4, p. 30-31
  template< typename Scalar_T >
  class numeric_traits
  {
  private:
    /// Smart isinf specialised for Scalar_T without infinity
    inline
    static 
    bool
    isInf(Scalar_T val, bool_to_type<false>)
    { return false; }

    /// Smart isnan specialised for Scalar_T with quiet NaN
    inline
    static 
    bool
    isInf(Scalar_T val, bool_to_type<true>)
    { return std::isinf(val); }

    /// Smart isnan specialised for Scalar_T without quiet NaN
    inline
    static 
    bool
    isNaN(Scalar_T val, bool_to_type<false>)
    { return false; }

    /// Smart isnan specialised for Scalar_T with quiet NaN
    inline
    static 
    bool
    isNaN(Scalar_T val, bool_to_type<true>)
    { return std::isnan(val); }

  public:
    /// Smart isinf
    inline
    static 
    bool
    isInf(Scalar_T val)
    {
      return isInf(val,
             bool_to_type< std::numeric_limits<Scalar_T>::has_infinity >() );
    }

    /// Smart isnan
    inline
    static 
    bool
    isNaN(Scalar_T val)
    {
      return isNaN(val,
             bool_to_type< std::numeric_limits<Scalar_T>::has_quiet_NaN >() );
    }

    /// Smart isnan or isinf
    inline
    static 
    bool
    isNaN_or_isInf(Scalar_T val)
    {
      return isNaN(val,
             bool_to_type< std::numeric_limits<Scalar_T>::has_quiet_NaN >() )
          || isInf(val,
             bool_to_type< std::numeric_limits<Scalar_T>::has_infinity >() );
    }

    /// Smart NaN
    inline
    static 
    const Scalar_T
    NaN()
    {
      return std::numeric_limits<Scalar_T>::has_quiet_NaN 
           ? std::numeric_limits<Scalar_T>::quiet_NaN() 
           : Scalar_T(std::log(0.0));
    }

    /// Absolute value of scalar
    inline
    static 
    const Scalar_T
    abs(Scalar_T val)
    { return boost::numeric::ublas::type_traits<Scalar_T>::type_abs(val); }

    /// Square root of scalar
    inline
    static 
    const Scalar_T
    sqrt(Scalar_T val)
    { return boost::numeric::ublas::type_traits<Scalar_T>::type_sqrt(val); }

    /// Logarithm of scalar
    inline
    static 
    const Scalar_T
    log(Scalar_T val)
    { return std::log(val); }
  };

  /// Log base 2 for Scalar_T
  template< typename Scalar_T >
  inline
  Scalar_T
  log2(const Scalar_T x)
  { return std::log(x)/Scalar_T(l_ln2); }
}
#endif // _GLUCAT_SCALAR_H
