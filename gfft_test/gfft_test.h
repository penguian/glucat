#ifndef GLUCAT_GFFT_TEST_H
#define GLUCAT_GFFT_TEST_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    gfft_test.h: GFFT test
                             -------------------
    begin                : Sun 2001-12-09
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
 "Clifford algebras with numeric and symbolic computations", Birkhauser, 1996.
 ***************************************************************************
     See also Arvind Raja's original header comments in glucat.h
 ***************************************************************************/

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
#include "test/try_catch.h"
#include <iomanip>

namespace glucat_gfft_test
{
  using namespace glucat;
  using namespace std;

  const index_t max_n = DEFAULT_HI;

  inline
  static
  double
  elapsed( clock_t cpu_time )
  { return ((clock() - cpu_time)*MS_PER_S) / CLOCKS_PER_SEC; }

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
    const streamsize width = 10;
    const streamsize old_prec = cout.precision();
    const streamsize prec = 2;
    const int trial_width = 4;
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

  template< typename Multivector_T >
  static
  void
  time_fast(Multivector_T& a, const Multivector_T& b,
            const typename Multivector_T::index_set_t& inner_frame,
            const typename Multivector_T::index_set_t& outer_frame)
  {
    typedef typename Multivector_T::matrix_multi_t matrix_multi_t;
    typedef typename Multivector_T::framed_multi_t framed_multi_t;
    typedef typename Multivector_T::index_set_t index_set_t;
    typedef typename Multivector_T::scalar_t    scalar_t;

    int trials;
    const double min_cpu_time = 20.0;
    const int min_trials_if_small = 1000;
    a = a*b*(scalar_t(1.0*rand())/RAND_MAX) + a*(scalar_t(1.0*rand())/RAND_MAX);
    int mm_trials = 1;
    clock_t cpu_time = clock();
     matrix_multi_t A = a.fast_matrix_multi(outer_frame);
    double mm_cpu_time = elapsed(cpu_time);
    if (mm_cpu_time < MS_PER_S)
    {
      const int max_trials =
       mm_cpu_time < min_cpu_time ?
         min_trials_if_small :
         int(ceil(MS_PER_S/mm_cpu_time));
      cpu_time = clock();
      for (trials = 0; trials != max_trials; ++trials)
        A = a.fast_matrix_multi(outer_frame);
      double trial_time = elapsed(cpu_time);
      matrix_multi_t new_A;
      cpu_time = clock();
      for (mm_trials = 0; mm_trials != trials; ++mm_trials)
        new_A = A;
      double for_time = elapsed(cpu_time);
      mm_cpu_time = max(0.0,trial_time - for_time)/trials;
    }
    int fm_trials = 1;
    cpu_time = clock();
     framed_multi_t new_a = A.fast_framed_multi();
    double fm_cpu_time = elapsed(cpu_time);
    if (fm_cpu_time < MS_PER_S)
    {
      const int max_trials =
       fm_cpu_time < min_cpu_time ?
         min_trials_if_small :
         int(ceil(MS_PER_S/fm_cpu_time));
      cpu_time = clock();
      for (trials = 0; trials != max_trials; ++trials)
        new_a = A.fast_framed_multi();
      double trial_time = elapsed(cpu_time);
      cpu_time = clock();
      for (fm_trials = 0; fm_trials != trials; ++fm_trials)
        new_a = a;
      double for_time = elapsed(cpu_time);
      fm_cpu_time = max(0.0,trial_time - for_time)/trials;
    }
    print_times(inner_frame, outer_frame, mm_cpu_time, fm_cpu_time, mm_trials, fm_trials);
    const ios::fmtflags& old_flags = cout.flags();
    const streamsize width = 8;
    const streamsize old_prec = cout.precision();
    const streamsize prec = 2;
    cout.setf(ios_base::scientific);
    cout.setf(ios_base::showpoint);
    cout << setprecision(prec)
         << " diff: " << setw(width) << abs(A.fast_framed_multi() - a)/abs(a)
         << setprecision(old_prec)  << endl;
    cout.flags(old_flags);
  }

  template< class Multivector_T >
  static
  void
  fast_test(const index_t n, const index_t max_n)
  {
    cout << "Clifford algebra transform test:" << endl;

    typedef Multivector_T m_;
    typedef typename m_::matrix_multi_t matrix_multi_t;
    typedef typename m_::framed_multi_t framed_multi_t;
    typedef typename m_::index_set_t e_;
    typedef typename m_::scalar_t    scalar_t;

    const index_t max_index = min(n, max_n);
    for (index_t k = max_index; k != 0; --k)
    {
      m_ a = 1;
      e_ inner_frame = e_();
      for (index_t i = k; i < max_index+1; i += k)
      {
        inner_frame |= e_(i);
        a = a*m_(e_(i) , 1.0)*(scalar_t(1.0*rand())/RAND_MAX) + a*(scalar_t(1.0*rand())/RAND_MAX);
        inner_frame |= e_(-i);
        time_fast(a, m_(e_(-i), 1.0), inner_frame, inner_frame);
      }
    }
  }
}

int gfft_test(const int n);

#endif // GLUCAT_GFFT_TEST_H
