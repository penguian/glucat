#ifndef _GLUCAT_GLOBAL_H
#define _GLUCAT_GLOBAL_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    global.h : Global declarations
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2016 by Paul C. Leopardi
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

  /// Default for truncation
  const double DEFAULT_TRUNCATION = std::numeric_limits<float>::epsilon();

  /// Precision policy
  enum precision_t
  {
    precision_demoted,
    precision_same,
    precision_promoted
  };

  // Tuning policy default constants
  const unsigned int DEFAULT_Mult_Matrix_Threshold   =       8;
  const unsigned int DEFAULT_Div_Max_Steps           =       4;
  const unsigned int DEFAULT_Sqrt_Max_Steps          =     256;
  const unsigned int DEFAULT_Log_Max_Outer_Steps     =     256;
  const unsigned int DEFAULT_Log_Max_Inner_Steps     =      32;
  const unsigned int DEFAULT_Basis_Max_Count         =      12;
  const unsigned int DEFAULT_Fast_Size_Threshold     = 1 <<  6;
  const unsigned int DEFAULT_Inv_Fast_Dim_Threshold  = 1 <<  3;
  const unsigned int DEFAULT_Products_Size_Threshold = 1 << 22;
  const precision_t  DEFAULT_Function_Precision      = precision_same;


  /// Tuning policy
  template
  <
    unsigned int Mult_Matrix_Threshold   = DEFAULT_Mult_Matrix_Threshold,
    unsigned int Div_Max_Steps           = DEFAULT_Div_Max_Steps,
    unsigned int Sqrt_Max_Steps          = DEFAULT_Sqrt_Max_Steps,
    unsigned int Log_Max_Outer_Steps     = DEFAULT_Log_Max_Outer_Steps,
    unsigned int Log_Max_Inner_Steps     = DEFAULT_Log_Max_Inner_Steps,
    unsigned int Basis_Max_Count         = DEFAULT_Basis_Max_Count,
    unsigned int Fast_Size_Threshold     = DEFAULT_Fast_Size_Threshold,
    unsigned int Inv_Fast_Dim_Threshold  = DEFAULT_Inv_Fast_Dim_Threshold,
    unsigned int Products_Size_Threshold = DEFAULT_Products_Size_Threshold,
    precision_t  Function_Precision      = DEFAULT_Function_Precision
  >
  struct tuning
  {
  // Tuning for multiplication
    /// Minimum index count needed to invoke matrix multiplication algorithm
    enum { mult_matrix_threshold = Mult_Matrix_Threshold };
  // Tuning for division
    /// Maximum steps of iterative refinement in division algorithm
    enum { div_max_steps = Div_Max_Steps };
  // Tuning for sqrt
    /// Maximum number of steps in square root iteration
    enum { sqrt_max_steps = Sqrt_Max_Steps };
  // Tuning for log
    /// Maximum number of incomplete square roots in cascade log algorithm
    enum { log_max_outer_steps = Log_Max_Outer_Steps };
    /// Maximum number of steps in incomplete square root within cascade log algorithm
    enum { log_max_inner_steps = Log_Max_Inner_Steps };
  // Tuning for basis cache
    /// Maximum index count of folded frames in basis cache
    enum { basis_max_count = Basis_Max_Count };
  // Tuning for FFT
    /// Minimum map size needed to invoke generalized FFT
    enum { fast_size_threshold = Fast_Size_Threshold };
    /// Minimum matrix dimension needed to invoke inverse generalized FFT
    enum { inv_fast_dim_threshold = Inv_Fast_Dim_Threshold };
  // Tuning for products (other than geometric product)
    /// Minimum size needed for to invoke faster products algorithms
    enum { products_size_threshold = Products_Size_Threshold };
  // Tuning for precision of exp, log and sqrt functions
    /// Precision used for exp, log and sqrt functions
    static const precision_t function_precision = Function_Precision;
  };

  /// Modulo function which works reliably for lhs < 0
  template< typename LHS_T, typename RHS_T >
  inline
  LHS_T
  pos_mod(LHS_T lhs, RHS_T rhs)
  { return lhs > 0? lhs % rhs : (-lhs) % rhs == 0 ? 0 : rhs - (-lhs) % rhs; }

}
#endif // _GLUCAT_GLOBAL_H
