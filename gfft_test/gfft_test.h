#ifndef GLUCAT_GFFT_TEST_H
#define GLUCAT_GFFT_TEST_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    gfft_test.h: GFFT test
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2012 by Paul C. Leopardi
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

#define _GLUCAT_TEST_REPEAT
#include "glucat/glucat.h"

// If radix of int is not 2, we can't easily set thresholds
_GLUCAT_CTAssert(std::numeric_limits<int>::radix == 2, CannotSetThresholds)

const int DRIVER_BASIS_MAX_COUNT        = 1;
const int DRIVER_FAST_SIZE_THRESHOLD    = 1 << (std::numeric_limits<int>::digits - 1);
const int DRIVER_INV_FAST_DIM_THRESHOLD = 1 << (std::numeric_limits<int>::digits - 1);
typedef glucat::tuning
  <
    glucat::DEFAULT_Mult_Matrix_Threshold,
    glucat::DEFAULT_Div_Max_Steps,
    glucat::DEFAULT_Sqrt_Max_Steps,
    glucat::DEFAULT_Log_Max_Outer_Steps,
    glucat::DEFAULT_Log_Max_Inner_Steps,
    DRIVER_BASIS_MAX_COUNT,
    DRIVER_FAST_SIZE_THRESHOLD,
    DRIVER_INV_FAST_DIM_THRESHOLD
  >
  Tune_P;

#include "glucat/glucat_imp.h"
#include "test/timing.h"
#include "test/try_catch.h"
#include <iomanip>

namespace glucat_gfft_test
{
  using namespace glucat;
  using namespace timing;
  using namespace std;

  const index_t max_n = DEFAULT_HI;

  template< typename Index_Set_T >
  inline
  static
  void
  print_times(const Index_Set_T& frame1, const Index_Set_T& frame2,
              const double mm_cpu_time,  const double fm_cpu_time,
              const int mm_trials,       const int fm_trials)
  {
    cout << "R_" << frame1 << " in "
         << "R_" << frame2 << ":" << endl;
    cout << " CPU = ";
    const ios::fmtflags& old_flags = cout.flags();
    const streamsize width = 12;
    const streamsize old_prec = cout.precision();
    const streamsize prec = 3;
    const int trial_width = 5;
    cout.setf(ios_base::fixed);
    cout.setf(ios_base::showpoint);
    cout << setprecision(prec)
         << "mm: "
         << setw(width) << mm_cpu_time << " "
         << "fm: "
         << setw(width) << fm_cpu_time << " "
         << "trials "
         << "mm: "
         << setw(trial_width) << mm_trials << " "
         << "fm: "
         << setw(trial_width) << fm_trials << " "
         << setprecision(old_prec);
    cout.flags(old_flags);
  }

  template< typename Scalar_T >
  inline
  static
  Scalar_T
  relative(const Scalar_T abs_lhs_minus_rhs, const Scalar_T abs_rhs)
  {
    return
      (abs_rhs == Scalar_T(0))
      ? abs_lhs_minus_rhs
      : abs_lhs_minus_rhs / abs_rhs;
  }

  template< typename Multivector_T >
  static
  void
  time_fast(const typename Multivector_T::index_set_t inner_frame,
            const typename Multivector_T::index_set_t outer_frame,
            const typename Multivector_T::scalar_t fill)
  {
    typedef typename Multivector_T::framed_multi_t framed_multi_t;
    typedef typename Multivector_T::matrix_multi_t matrix_multi_t;
    typedef typename Multivector_T::scalar_t scalar_t;

    Multivector_T a = Multivector_T::random(inner_frame, fill);
    clock_t cpu_time = clock();
      matrix_multi_t A = a.template fast_matrix_multi<scalar_t>(outer_frame);
    double mm_cpu_time = elapsed(cpu_time);
    int mm_trials = 1;
#ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; mm_cpu_time == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
        for (int trials = 0; trials != nbr_trials; ++trials)
          A = a.template fast_matrix_multi<scalar_t>(outer_frame);
      mm_cpu_time = elapsed(cpu_time) / nbr_trials;
      mm_trials = nbr_trials;
    }
#endif
    cpu_time = clock();
      framed_multi_t new_a = A.template fast_framed_multi<scalar_t>();
    double fm_cpu_time = elapsed(cpu_time);
    int fm_trials = 1;
#ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; fm_cpu_time == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
        for (int trials = 0; trials != nbr_trials; ++trials)
          new_a = A.template fast_framed_multi<scalar_t>();
      fm_cpu_time = elapsed(cpu_time) / nbr_trials;
      fm_trials = nbr_trials;
    }
#endif
    print_times(inner_frame, outer_frame, mm_cpu_time, fm_cpu_time, mm_trials, fm_trials);
    const ios::fmtflags& old_flags = cout.flags();
    const streamsize width = 8;
    const streamsize old_prec = cout.precision();
    const streamsize prec = 2;
    cout.setf(ios_base::scientific);
    cout.setf(ios_base::showpoint);
    cout << setprecision(prec)
         << " diff: " << setw(width) << relative(abs(new_a - a), abs(a))
         << setprecision(old_prec)  << endl;
    cout.flags(old_flags);
  }

  template< typename Multivector_T >
  static
  void
  fast_test_with_fill(const index_t n, const index_t max_n,
                      const typename Multivector_T::scalar_t fill)
  {
    typedef typename Multivector_T::index_set_t index_set_t;

    const index_t max_index = min(n, max_n);
    cout << "Fill: " << fill << endl;
    for (index_t k = max_index; k != 0; --k)
    {
      index_set_t inner_frame = index_set_t();
      for (index_t i = k; i < max_index+1; i += k)
      {
        inner_frame |= index_set_t(-i) | index_set_t(i);
        time_fast<Multivector_T>(inner_frame, inner_frame, fill);
      }
    }
  }

  template< typename Multivector_T >
  static
  void
  fast_test(const index_t n, const index_t max_n)
  {
    cout << "Clifford algebra transform test:" << endl;

    fast_test_with_fill<Multivector_T>(n, max_n, 0.5);
    fast_test_with_fill<Multivector_T>(n, max_n, 1.0);
  }
}

int gfft_test(const int n);

#endif // GLUCAT_GFFT_TEST_H
