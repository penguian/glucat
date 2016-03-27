#ifndef GLUCAT_TRANSFORMS_H
#define GLUCAT_TRANSFORMS_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    transforms.h : Transforms
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
#include "test/undefine.h"
#define _GLUCAT_TEST_TUNING_NAIVE
#include "test/tuning.h"
#include "glucat/glucat_imp.h"
#include "test/timing.h"
#include "test/try_catch.h"
#include <stdio.h>
#include <iomanip>

namespace glucat_fast_test
{
  using namespace glucat;
  using namespace glucat::timing;
  using namespace std;

  const index_t max_n = DEFAULT_HI;
  typedef index_set< -max_n, max_n > index_set_t;

  inline
  static
  void
  print_times(const index_set_t& frame1, const index_set_t& frame2,
              const double old_mm_cpu_time,
              const double new_mm_cpu_time,
              const double old_fm_cpu_time,
              const double new_fm_cpu_time)
  {
    const int index_width = 2;
    cout << "Cl(" << setw(index_width) <<  max_pos(frame1) << ","
                  << setw(index_width) << -min_neg(frame1) << ") in "
         << "Cl(" << setw(index_width) <<  max_pos(frame2) << ","
                  << setw(index_width) << -min_neg(frame2) << ")"
         << " CPU = ";
    const ios::fmtflags& old_flags = cout.flags();
    const streamsize width = 10;
    cout.setf(ios_base::fixed);
    cout.setf(ios_base::showpoint);
    const streamsize old_prec = cout.precision(2);
    cout << "mm: "
         << setw(width) << old_mm_cpu_time << " (old) "
         << setw(width) << new_mm_cpu_time << " (new) "
         << "fm: "
         << setw(width) << old_fm_cpu_time << " (old) "
         << setw(width) << new_fm_cpu_time << " (new) ";
    cout.precision(old_prec);
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
  time_fast(const index_set_t inner_frame,
            const index_set_t outer_frame,
            const typename Multivector_T::scalar_t fill)
  {
    typedef typename Multivector_T::matrix_multi_t matrix_multi_t;
    typedef typename Multivector_T::framed_multi_t framed_multi_t;
    typedef typename Multivector_T::scalar_t scalar_t;

    Multivector_T a = Multivector_T::random(inner_frame, fill);
    clock_t cpu_time = clock();
      matrix_multi_t old_A = matrix_multi_t(a, outer_frame);
    double old_mm_cpu_time = elapsed(cpu_time);
#ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; old_mm_cpu_time == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
        for (int trials = 0; trials != nbr_trials; ++trials)
          old_A = matrix_multi_t(a, outer_frame);
      old_mm_cpu_time = elapsed(cpu_time) / nbr_trials;
    }
#endif
    cpu_time = clock();
      matrix_multi_t new_A = a.template fast_matrix_multi<scalar_t>(outer_frame);
    double new_mm_cpu_time = elapsed(cpu_time);
#ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; new_mm_cpu_time == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
        for (int trials = 0; trials != nbr_trials; ++trials)
          new_A = a.template fast_matrix_multi<scalar_t>(outer_frame);
      new_mm_cpu_time = elapsed(cpu_time) / nbr_trials;
    }
#endif
    cpu_time = clock();
      framed_multi_t old_a = old_A;
    double old_fm_cpu_time = elapsed(cpu_time);
#ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; old_fm_cpu_time == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
        for (int trials = 0; trials != nbr_trials; ++trials)
          old_a = old_A;
      old_fm_cpu_time = elapsed(cpu_time) / nbr_trials;
    }
#endif
    cpu_time = clock();
      framed_multi_t new_a = new_A.template fast_framed_multi<scalar_t>();
    double new_fm_cpu_time = elapsed(cpu_time);
#ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; new_fm_cpu_time == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
        for (int trials = 0; trials != nbr_trials; ++trials)
          new_a = new_A.template fast_framed_multi<scalar_t>();
      new_fm_cpu_time = elapsed(cpu_time) / nbr_trials;
    }
#endif
    print_times(inner_frame, outer_frame,
                old_mm_cpu_time, new_mm_cpu_time,
                old_fm_cpu_time, new_fm_cpu_time);
    const ios::fmtflags& old_flags = cout.flags();
    const streamsize width = 8;
    cout.setf(ios_base::scientific);
    cout.setf(ios_base::showpoint);
    const streamsize old_prec = cout.precision(2);
    cout << " diff: old: " << setw(width) << relative(abs(old_A - old_a), abs(old_a))
               << " new: " << setw(width) << relative(abs(new_A - new_a), abs(new_a))
               << " fm: "  << setw(width) << relative(abs(new_A - old_A), abs(old_A))
               << " mm: "  << setw(width) << relative(abs(new_a - old_a), abs(old_a))
               << endl;
    cout.precision(old_prec);
    cout.flags(old_flags);
  }

  template< typename Multivector_T >
  static
  void
  fast_test_with_fill(const index_t max_pos, const index_t max_index,
                      const typename Multivector_T::scalar_t fill)
  {
    cout << "Fill: " << fill << endl;

    typedef typename Multivector_T::index_set_t index_set_t;
    const index_set_t pos_frame   = index_set_t(make_pair(1, max_pos));
    const index_set_t outer_frame = index_set_t(make_pair(-max_index, max_index));

    index_set_t inner_frame = index_set_t();
    for (index_t i = 1; i != max_pos+1; ++i)
    {
      inner_frame |= index_set_t(i);
      time_fast<Multivector_T>(inner_frame, pos_frame, fill);
    }
    inner_frame = index_set_t();
    for (index_t i = 1; i != max_index+1; ++i)
    {
      inner_frame |= index_set_t(i);
      time_fast<Multivector_T>(inner_frame, outer_frame, fill);
      inner_frame |= index_set_t(-i);
      time_fast<Multivector_T>(inner_frame, outer_frame, fill);
    }
    inner_frame = index_set_t();
    for (index_t i = 1; i != max_index+1; ++i)
    {
      inner_frame |= index_set_t(i);
      time_fast<Multivector_T>(inner_frame, inner_frame, fill);
      inner_frame |= index_set_t(-i);
      time_fast<Multivector_T>(inner_frame, inner_frame, fill);
    }
  }

  template< typename Multivector_T >
  static
  void
  fast_test(const index_t n, const index_t max_n)
  {
    cout << "Clifford algebra transform test:" << endl;

    const index_t max_pos = min(index_t(2*n), max_n);
    const index_t max_index = min(n, max_n);

    fast_test_with_fill<Multivector_T>(max_pos, max_index, 0.5);
    fast_test_with_fill<Multivector_T>(max_pos, max_index, 1.0);
  }
}

int transforms(const int n);

#endif // GLUCAT_TRANSFORMS_H
