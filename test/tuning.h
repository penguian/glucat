#ifndef GLUCAT_TEST_TUNING_H
#define GLUCAT_TEST_TUNING_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    tuning.h : Class definitions to control test tuning
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

namespace glucat
{
  const unsigned int Tuning_Int_Digits = std::numeric_limits<int>::digits;
  const unsigned int Tuning_Max_Threshold = 1 << Tuning_Int_Digits;

  // Specific tuning policy constants and tuning policies

  const unsigned int Tuning_Slow_Mult_Matrix_Threshold   = Tuning_Max_Threshold;
  const unsigned int Tuning_Slow_Basis_Max_Count         =       0;
  const unsigned int Tuning_Slow_Fast_Size_Threshold     = Tuning_Max_Threshold;
  const unsigned int Tuning_Slow_Inv_Fast_Dim_Threshold  = Tuning_Max_Threshold;
  const unsigned int Tuning_Slow_Products_Size_Threshold = Tuning_Max_Threshold;

  using tuning_values_slow_p = tuning_values
    <
      Tuning_Slow_Mult_Matrix_Threshold,
      Tuning_Default_Div_Max_Steps,
      Tuning_Default_CR_Sqrt_Max_Steps,
      Tuning_Default_DB_Sqrt_Max_Steps,
      Tuning_Default_Log_Max_Outer_Steps,
      Tuning_Default_Log_Max_Inner_Steps,
      Tuning_Slow_Basis_Max_Count,
      Tuning_Slow_Fast_Size_Threshold,
      Tuning_Slow_Inv_Fast_Dim_Threshold,
      Tuning_Slow_Products_Size_Threshold,
      Tuning_Default_Denom_Different_Bits,
      Tuning_Default_Extra_Different_Bits
    >;

  using tuning_slow_p = tuning<tuning_values_slow_p, precision_same>;

  const unsigned int Tuning_Naive_Mult_Matrix_Threshold   =       0;
  const unsigned int Tuning_Naive_Basis_Max_Count         = Tuning_Max_Threshold;
  const unsigned int Tuning_Naive_Fast_Size_Threshold     = Tuning_Max_Threshold;
  const unsigned int Tuning_Naive_Inv_Fast_Dim_Threshold  = Tuning_Max_Threshold;

  using tuning_values_naive_p = tuning_values
    <
      Tuning_Naive_Mult_Matrix_Threshold,
      Tuning_Default_Div_Max_Steps,
      Tuning_Default_CR_Sqrt_Max_Steps,
      Tuning_Default_DB_Sqrt_Max_Steps,
      Tuning_Default_Log_Max_Outer_Steps,
      Tuning_Default_Log_Max_Inner_Steps,
      Tuning_Naive_Basis_Max_Count,
      Tuning_Naive_Fast_Size_Threshold,
      Tuning_Naive_Inv_Fast_Dim_Threshold,
      Tuning_Default_Products_Size_Threshold,
      Tuning_Default_Denom_Different_Bits,
      Tuning_Default_Extra_Different_Bits
    >;

  using tuning_naive_p = tuning<tuning_values_naive_p, precision_same>;

  const unsigned int Tuning_Fast_Mult_Matrix_Threshold   =       0;
  const unsigned int Tuning_Fast_Div_Max_Steps           =       0;
  const unsigned int Tuning_Fast_CR_Sqrt_Max_Steps       =     256;
  const unsigned int Tuning_Fast_DB_Sqrt_Max_Steps       =     256;
  const unsigned int Tuning_Fast_Log_Max_Outer_Steps     =      16;
  const unsigned int Tuning_Fast_Log_Max_Inner_Steps     =       8;
  const unsigned int Tuning_Fast_Basis_Max_Count         =       1;
  const unsigned int Tuning_Fast_Fast_Size_Threshold     =       0;
  const unsigned int Tuning_Fast_Inv_Fast_Dim_Threshold  =       0;
  const unsigned int Tuning_Fast_Products_Size_Threshold =       0;

  using tuning_values_fast_p = tuning_values
    <
      Tuning_Fast_Mult_Matrix_Threshold,
      Tuning_Fast_Div_Max_Steps,
      Tuning_Fast_CR_Sqrt_Max_Steps,
      Tuning_Fast_DB_Sqrt_Max_Steps,
      Tuning_Fast_Log_Max_Outer_Steps,
      Tuning_Fast_Log_Max_Inner_Steps,
      Tuning_Fast_Basis_Max_Count,
      Tuning_Fast_Fast_Size_Threshold,
      Tuning_Fast_Inv_Fast_Dim_Threshold,
      Tuning_Fast_Products_Size_Threshold,
      Tuning_Default_Denom_Different_Bits,
      Tuning_Default_Extra_Different_Bits
    >;

  using tuning_fast_p = tuning<tuning_values_fast_p, precision_same>;

}
#endif // GLUCAT_TEST_TUNING_H
