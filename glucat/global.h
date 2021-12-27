#ifndef _GLUCAT_GLOBAL_H
#define _GLUCAT_GLOBAL_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    global.h : Global declarations
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2021 by Paul C. Leopardi
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

#include "glucat/portability.h"

#include <limits>
#include <climits>

namespace glucat
{
  // References:
  // [AA]: A. Alexandrescu, "Modern C++ Design", Addison-Wesley, 2001.

  /// Compile time assertion
  // Reference: [AA], p. 25
  template<bool> struct CTAssertion;
  template<> struct CTAssertion<true> { };
  #define _GLUCAT_CTAssert(expr, msg) \
      namespace { struct msg { glucat::CTAssertion<(expr)> ERROR_##msg; }; }

  /// Type comparison
  // Reference: [AA], pp. 34--37
  template < typename LHS_T, typename RHS_T >
  class compare_types
  {
  public:
    enum { are_same = false };
  };
  template < typename T >
  class compare_types<T, T>
  {
  public:
    enum { are_same = true };
  };

  /// Bool to type
  // Reference: [AA], 2.4, p. 29
  template< bool truth_value >
  class bool_to_type
  {
  private:
    enum { value = truth_value };
  };

  // Global types which determine sizes
  /// Size of index_t should be enough to represent LO, HI
  using index_t = int;
  /// Size of set_value_t should be enough to contain index_set<LO,HI>
  using set_value_t = unsigned long;

  // Global constants
  /// Timing constant: deprecated here - moved to test/timing.h
  const double MS_PER_S = 1000.0;

  // Constants which determine sizes

  // Bits per unsigned long
  #if   (ULONG_MAX == (4294967295UL))
  #define _GLUCAT_BITS_PER_ULONG 32
  #elif (ULONG_MAX == (18446744073709551615UL))
  #define _GLUCAT_BITS_PER_ULONG 64
  #elif defined(__WORDSIZE)
  #define _GLUCAT_BITS_PER_ULONG __WORDSIZE
  #endif

  /// If radix of unsigned char is not 2, we can't easily determine number of bits from sizeof
  _GLUCAT_CTAssert(std::numeric_limits<unsigned char>::radix == 2, CannotDetermineBitsPerChar)

  /// Number of bits per char is used to determine number of bits in set_value_t
  const index_t BITS_PER_CHAR = std::numeric_limits<unsigned char>::digits;

  /// Number of bits in set_value_t
  const index_t BITS_PER_SET_VALUE = std::numeric_limits<set_value_t>::digits;

  _GLUCAT_CTAssert(_GLUCAT_BITS_PER_ULONG == BITS_PER_SET_VALUE, BitsPerULongDoesNotMatchSetValueT)

  // Constants which are determined by size
  /// Default lowest index in an index set
  const index_t DEFAULT_LO = -index_t(BITS_PER_SET_VALUE / 2);
  /// Default highest index in an index set
  const index_t DEFAULT_HI =  index_t(BITS_PER_SET_VALUE / 2);

  /// Modulo function which works reliably for lhs < 0
  template< typename LHS_T, typename RHS_T >
  inline
  auto
  pos_mod(LHS_T lhs, RHS_T rhs) -> LHS_T
  { return lhs > 0? lhs % rhs : (-lhs) % rhs == 0 ? 0 : rhs - (-lhs) % rhs; }

}
#endif // _GLUCAT_GLOBAL_H
