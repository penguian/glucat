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
  // Base class for policies
  struct policy{};

  // Precision policy
  enum precision_t
  {
    precision_demoted,
    precision_same,
    precision_promoted
  };
  const precision_t  Tuning_Default_Function_Precision      = precision_same;

  // Tuning policy default constants
  const unsigned int Tuning_Default_Mult_Matrix_Threshold   =       4;
  const unsigned int Tuning_Default_CR_Sqrt_Max_Steps       =     256;
  const unsigned int Tuning_Default_DB_Sqrt_Max_Steps       =     256;
  const unsigned int Tuning_Default_Log_Max_Outer_Steps     =     256;
  const unsigned int Tuning_Default_Log_Max_Inner_Steps     =      32;
  const unsigned int Tuning_Default_Basis_Max_Count         =      12;
  const unsigned int Tuning_Default_Fast_Size_Threshold     =      16;
  const unsigned int Tuning_Default_Inv_Fast_Dim_Threshold  =       3;
  const unsigned int Tuning_Default_Denom_Different_Bits    =       8;
  const unsigned int Tuning_Default_Extra_Different_Bits    =       8;
  const unsigned int Tuning_Default_Products_Different_Bits   =     2;
  const unsigned int Tuning_Default_Products_Matrix_Threshold =    14;

  /// Tuning policy values
  template
  <
  unsigned int Mult_Matrix_Threshold   = Tuning_Default_Mult_Matrix_Threshold,
  unsigned int CR_Sqrt_Max_Steps       = Tuning_Default_CR_Sqrt_Max_Steps,
  unsigned int DB_Sqrt_Max_Steps       = Tuning_Default_DB_Sqrt_Max_Steps,
  unsigned int Log_Max_Outer_Steps     = Tuning_Default_Log_Max_Outer_Steps,
  unsigned int Log_Max_Inner_Steps     = Tuning_Default_Log_Max_Inner_Steps,
  unsigned int Basis_Max_Count         = Tuning_Default_Basis_Max_Count,
  unsigned int Fast_Size_Threshold     = Tuning_Default_Fast_Size_Threshold,
  unsigned int Inv_Fast_Dim_Threshold  = Tuning_Default_Inv_Fast_Dim_Threshold,
  unsigned int Denom_Different_Bits    = Tuning_Default_Denom_Different_Bits,
  unsigned int Extra_Different_Bits    = Tuning_Default_Extra_Different_Bits,
  unsigned int Products_Different_Bits = Tuning_Default_Products_Different_Bits,
  unsigned int Products_Matrix_Threshold = Tuning_Default_Products_Matrix_Threshold
  >
  struct tuning_values : policy
  {
    using tuning_values_p = tuning_values
    <
    Mult_Matrix_Threshold,
    CR_Sqrt_Max_Steps,
    DB_Sqrt_Max_Steps,
    Log_Max_Outer_Steps,
    Log_Max_Inner_Steps,
    Basis_Max_Count,
    Fast_Size_Threshold,
    Inv_Fast_Dim_Threshold,
    Denom_Different_Bits,
    Extra_Different_Bits,
    Products_Different_Bits,
    Products_Matrix_Threshold
    >;
    // Tuning for multiplication
    // Minimum index count needed to invoke matrix multiplication algorithm
    static constexpr unsigned int mult_matrix_threshold = Mult_Matrix_Threshold;
    // Tuning for sqrt
    // Maximum number of steps in cyclic reduction square root iteration
    static constexpr unsigned int cr_sqrt_max_steps = CR_Sqrt_Max_Steps;
    // Maximum number of steps in Denman-Beavers square root iteration
    static constexpr unsigned int db_sqrt_max_steps = DB_Sqrt_Max_Steps;
    // Tuning for log
    // Maximum number of incomplete square roots in cascade log algorithm
    static constexpr unsigned int log_max_outer_steps = Log_Max_Outer_Steps;
    // Maximum number of steps in incomplete square root within cascade log algorithm
    static constexpr unsigned int log_max_inner_steps = Log_Max_Inner_Steps;
    // Tuning for basis cache
    // Maximum index count of folded frames in basis cache
    static constexpr unsigned int basis_max_count = Basis_Max_Count;
    // Tuning for FFT
    // Minimum map size needed to invoke generalized FFT
    static constexpr unsigned int fast_size_threshold = Fast_Size_Threshold;
    // Minimum matrix dimension needed to invoke inverse generalized FFT
    static constexpr unsigned int inv_fast_dim_threshold = Inv_Fast_Dim_Threshold;
    // Tuning for precision of exp, log and sqrt functions
    // Denominator of proportion of different bits allowed in approximate equality
    static constexpr unsigned int denom_different_bits = Denom_Different_Bits;
    // Extra number of different bits allowed in approximate equality
    static constexpr unsigned int extra_different_bits = Extra_Different_Bits;
    // Tuning for matrix operators pruning threshold: number of different bits allowed
    static constexpr unsigned int products_different_bits = Products_Different_Bits;
    // Minimum index count needed to invoke matrix non-geometric multiplication algorithm
    static constexpr unsigned int products_matrix_threshold = Products_Matrix_Threshold;
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
    // Precision used for exp, log and sqrt functions
    static const precision_t function_precision = Function_Precision;
    // Tuning used for return values of exp, log and sqrt functions
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
