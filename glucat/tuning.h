#ifndef GLUCAT_TUNING_H
#define GLUCAT_TUNING_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    tuning.h : Policy classes to control tuning
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
 "Clifford algebras with numeric and symbolic computations", Birkhauser, 1996.
 ***************************************************************************
     See also Arvind Raja's original header comments in glucat.h
 ***************************************************************************/

// If radix of int is not 2, we can't easily set thresholds
_GLUCAT_CTAssert(std::numeric_limits<unsigned int>::radix == 2, CannotSetThresholds)

namespace glucat
{
  /// Base class for policies
  struct policy{};

  /// Precision policy
  enum precision_t
  {
    precision_demoted,
    precision_same,
    precision_promoted
  };
  const precision_t  Tuning_Default_Function_Precision      = precision_same;

  // Tuning policy default constants
  const unsigned int Tuning_Default_Mult_Matrix_Threshold   =       8;
  const unsigned int Tuning_Default_Div_Max_Steps           =       4;
  const unsigned int Tuning_Default_CR_Sqrt_Max_Steps       =     256;
  const unsigned int Tuning_Default_DB_Sqrt_Max_Steps       =     256;
  const unsigned int Tuning_Default_Log_Max_Outer_Steps     =     256;
  const unsigned int Tuning_Default_Log_Max_Inner_Steps     =      32;
  const unsigned int Tuning_Default_Basis_Max_Count         =      12;
  const unsigned int Tuning_Default_Fast_Size_Threshold     = 1 <<  6;
  const unsigned int Tuning_Default_Inv_Fast_Dim_Threshold  = 1 <<  3;
  const unsigned int Tuning_Default_Products_Size_Threshold = 1 << 22;
  const unsigned int Tuning_Default_Denom_Different_Bits    =       8;
  const unsigned int Tuning_Default_Extra_Different_Bits    =       8;

  /// Tuning policy values
  template
  <
  unsigned int Mult_Matrix_Threshold   = Tuning_Default_Mult_Matrix_Threshold,
  unsigned int Div_Max_Steps           = Tuning_Default_Div_Max_Steps,
  unsigned int CR_Sqrt_Max_Steps       = Tuning_Default_CR_Sqrt_Max_Steps,
  unsigned int DB_Sqrt_Max_Steps       = Tuning_Default_DB_Sqrt_Max_Steps,
  unsigned int Log_Max_Outer_Steps     = Tuning_Default_Log_Max_Outer_Steps,
  unsigned int Log_Max_Inner_Steps     = Tuning_Default_Log_Max_Inner_Steps,
  unsigned int Basis_Max_Count         = Tuning_Default_Basis_Max_Count,
  unsigned int Fast_Size_Threshold     = Tuning_Default_Fast_Size_Threshold,
  unsigned int Inv_Fast_Dim_Threshold  = Tuning_Default_Inv_Fast_Dim_Threshold,
  unsigned int Products_Size_Threshold = Tuning_Default_Products_Size_Threshold,
  unsigned int Denom_Different_Bits    = Tuning_Default_Denom_Different_Bits,
  unsigned int Extra_Different_Bits    = Tuning_Default_Extra_Different_Bits
  >
  struct tuning_values : policy
  {
    using tuning_values_p = tuning_values
    <
    Mult_Matrix_Threshold,
    Div_Max_Steps,
    CR_Sqrt_Max_Steps,
    DB_Sqrt_Max_Steps,
    Log_Max_Outer_Steps,
    Log_Max_Inner_Steps,
    Basis_Max_Count,
    Fast_Size_Threshold,
    Inv_Fast_Dim_Threshold,
    Products_Size_Threshold,
    Denom_Different_Bits,
    Extra_Different_Bits
    >;
    // Tuning for multiplication
    /// Minimum index count needed to invoke matrix multiplication algorithm
    enum { mult_matrix_threshold = Mult_Matrix_Threshold };
    // Tuning for division
    /// Maximum steps of iterative refinement in division algorithm
    enum { div_max_steps = Div_Max_Steps };
    // Tuning for sqrt
    /// Maximum number of steps in cyclic reduction square root iteration
    enum { cr_sqrt_max_steps = CR_Sqrt_Max_Steps };
    /// Maximum number of steps in Denman-Beavers square root iteration
    enum { db_sqrt_max_steps = DB_Sqrt_Max_Steps };
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
    /// Denominator of proportion of different bits allowed in approximate equality
    enum { denom_different_bits = Denom_Different_Bits };
    /// Extra number of different bits allowed in approximate equality
    enum { extra_different_bits = Extra_Different_Bits };
  };

  using default_tuning_values_p = tuning_values<>;

  /// Tuning policy constants
  template
  <
    typename Tuning_Values_P = default_tuning_values_p,
    precision_t  Function_Precision = Tuning_Default_Function_Precision
  >
  struct tuning : policy
  {
    using tune_p = tuning
    <
      Tuning_Values_P,
      Function_Precision
    >;
    using tuning_values_p = Tuning_Values_P;
    /// Precision used for exp, log and sqrt functions
    static const precision_t function_precision = Function_Precision;
    /// Tuning used for return values of exp, log and sqrt functions
    using tuning_same_p = tuning
    <
      Tuning_Values_P,
      precision_same
    >;
  };

  using default_tuning_demoted_p = tuning
  <
    default_tuning_values_p,
    precision_demoted
  >;

  using default_tuning_same_p = tuning
  <
    default_tuning_values_p,
    precision_same
  >;

  using default_tuning_promoted_p = tuning
  <
    default_tuning_values_p,
    precision_promoted
  >;
}

#endif // GLUCAT_TUNING_H
