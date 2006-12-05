#ifndef GLUCAT_TRANSFORMS_H
#define GLUCAT_TRANSFORMS_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    transforms.h : Transforms
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

#include "glucat/glucat.h"
const int DRIVER_FAST_SIZE_THRESHOLD    = 1 << 30;
const int DRIVER_INV_FAST_DIM_THRESHOLD = 1 << 15;
typedef glucat::tuning
  <
    glucat::DEFAULT_Mult_Matrix_Threshold,
    glucat::DEFAULT_Div_Max_Steps,
    glucat::DEFAULT_Sqrt_Max_Steps,
    glucat::DEFAULT_Log_Max_Outer_Steps,
    glucat::DEFAULT_Log_Max_Inner_Steps,
    glucat::DEFAULT_Basis_Max_Count,
    DRIVER_FAST_SIZE_THRESHOLD,
    DRIVER_INV_FAST_DIM_THRESHOLD
  >
  Tune_P;
#include "glucat/glucat_imp.h"
#include "test/try_catch.h"
#include <stdio.h>
#include <iomanip>

namespace glucat_fast_test
{
  using namespace glucat;
  using namespace std;
  const index_t max_n = DEFAULT_HI;
  typedef double Scalar_T;
  typedef index_set< -max_n, max_n > index_set_t;

  inline
  static
  double
  elapsed( clock_t cpu_time )
  { return ((clock() - cpu_time)*MS_PER_S) / CLOCKS_PER_SEC; }

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

  template< typename Multivector_T >
  static
  void
  time_fast(Multivector_T& a, const Multivector_T& b,
            const index_set_t& inner_frame,
            const index_set_t& outer_frame)
  {
    typedef typename Multivector_T::matrix_multi_t matrix_multi_t;
    typedef typename Multivector_T::framed_multi_t framed_multi_t;
    int trials;
    const int min_trials_if_zero = 1000;
    a = a*b*(Scalar_T(1.0*rand())/RAND_MAX) + a*(Scalar_T(1.0*rand())/RAND_MAX);
    clock_t cpu_time = clock();
     matrix_multi_t old_A = matrix_multi_t(a, outer_frame);
    double old_mm_cpu_time = elapsed(cpu_time);
    if (old_mm_cpu_time < MS_PER_S)
    {
      const int max_trials =
       old_mm_cpu_time == 0.0 ?
         min_trials_if_zero :
         int(floor(MS_PER_S/old_mm_cpu_time));
      cpu_time = clock();
      for (trials = 0; trials != max_trials; ++trials)
        old_A = matrix_multi_t(a, outer_frame);
      old_mm_cpu_time = elapsed(cpu_time)/trials;
    }
    cpu_time = clock();
     matrix_multi_t new_A = a.fast_matrix_multi(outer_frame);
    double new_mm_cpu_time = elapsed(cpu_time);
    if (new_mm_cpu_time < MS_PER_S)
    {
      const int max_trials =
       new_mm_cpu_time == 0.0 ?
         min_trials_if_zero :
         int(floor(MS_PER_S/new_mm_cpu_time));
      cpu_time = clock();
      for (trials = 0; trials != max_trials; ++trials)
        new_A = a.fast_matrix_multi(outer_frame);
      new_mm_cpu_time = elapsed(cpu_time)/trials;
    }
    cpu_time = clock();
     framed_multi_t old_a = old_A;
    double old_fm_cpu_time = elapsed(cpu_time);
    if (old_fm_cpu_time < MS_PER_S)
    {
      const int max_trials =
       old_fm_cpu_time == 0.0 ?
         min_trials_if_zero :
         int(floor(MS_PER_S/old_fm_cpu_time));
      cpu_time = clock();
      for (trials = 0; trials != max_trials; ++trials)
        old_a = old_A;
      old_fm_cpu_time = elapsed(cpu_time)/trials;
    }
    cpu_time = clock();
     framed_multi_t new_a = new_A.fast_framed_multi();
    double new_fm_cpu_time = elapsed(cpu_time);
    if (new_fm_cpu_time < MS_PER_S)
    {
      const int max_trials =
       new_fm_cpu_time == 0.0 ?
         min_trials_if_zero :
         int(floor(MS_PER_S/new_fm_cpu_time));
      cpu_time = clock();
      for (trials = 0; trials != max_trials; ++trials)
        new_a = new_A.fast_framed_multi();
      new_fm_cpu_time = elapsed(cpu_time)/trials;
    }
    print_times(inner_frame, outer_frame,
      old_mm_cpu_time, new_mm_cpu_time, old_fm_cpu_time, new_fm_cpu_time);
    const ios::fmtflags& old_flags = cout.flags();
    const streamsize width = 8;
    cout.setf(ios_base::scientific);
    cout.setf(ios_base::showpoint);
    const streamsize old_prec = cout.precision(2);
    cout << " diff: old: " << setw(width) << abs(old_A - old_a)/abs(old_a)
               << " new: " << setw(width) << abs(new_A - new_a)/abs(new_a)
               << " fm: "  << setw(width) << abs(new_A - old_A)/abs(old_A)
               << " mm: "  << setw(width) << abs(new_a - old_a)/abs(old_a)
               << endl;
    cout.precision(old_prec);
    cout.flags(old_flags);
  }

  template< class Multivector_T >
  static
  void
  fast_test(const index_t n, const index_t max_n)
  {
    cout << "Clifford algebra transform test:" << endl;

    typedef Multivector_T m_;
    typedef index_set_t e_;
    typedef typename m_::matrix_multi_t matrix_multi_t;
    typedef typename m_::framed_multi_t framed_multi_t;
    static const index_t v_lo = e_::v_lo;
    static const index_t v_hi = e_::v_hi;

    m_ a;
    const index_t max_pos = min(index_t(2*n), max_n);
    const e_ pos_frame = index_range<v_lo,v_hi>(1, max_pos);

    const index_t max_index = min(n, max_n);
    const e_ outer_frame = index_range<v_lo,v_hi>(-max_index, max_index);
    a = m_(1, pos_frame);
    e_ inner_frame;
    for (index_t i = 1; i != max_pos+1; i++)
    {
      inner_frame |= e_(i);
      time_fast(a, m_(e_(i) , 1.0, pos_frame), inner_frame, pos_frame);
    }
    inner_frame = e_();
    a = m_(1, outer_frame);
    for (index_t i = 1; i != max_index+1; i++)
    {
      inner_frame |= e_(i);
      time_fast(a, m_(e_(i) , 1.0, outer_frame), inner_frame, outer_frame);
      inner_frame |= e_(-i);
      time_fast(a, m_(e_(-i), 1.0, outer_frame), inner_frame, outer_frame);
    }
    a = 1;
    inner_frame = e_();
    for (index_t i = 1; i != max_index+1; i++)
    {
      inner_frame |= e_(i);
      time_fast(a, m_(e_(i) , 1.0), inner_frame, inner_frame);
      inner_frame |= e_(-i);
      time_fast(a, m_(e_(-i), 1.0), inner_frame, inner_frame);
    }
  }
}

int transforms(const int n);

#endif // GLUCAT_TRANSFORMS_H
