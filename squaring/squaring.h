#ifndef GLUCAT_TEST_SQUARING_H
#define GLUCAT_TEST_SQUARING_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    squaring.h : multiplication speed test
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
#include "squaring/tuning.h"
#include "glucat/glucat_imp.h"
#include <stdio.h>
#include <iomanip>

namespace glucat_mult_test
{
  using namespace glucat;
  using namespace std;
  const index_t max_n = DEFAULT_HI;
  typedef double Scalar_T;
  typedef index_set< -max_n, max_n > index_set_t;

  inline
  double
  elapsed( clock_t cpu_time )
  { return ((clock() - cpu_time)*MS_PER_S) / CLOCKS_PER_SEC; }

  inline
  void
  print_times(const index_set_t& frame1, const index_set_t& frame2,
              const double setup_cpu_time, const double mult_cpu_time)
  {
    const int index_width = 2;
    cout << "Cl(" << setw(index_width) <<  max_pos(frame1) << ","
                  << setw(index_width) << -min_neg(frame1) << ") in "
         << "Cl(" << setw(index_width) <<  max_pos(frame2) << ","
                  << setw(index_width) << -min_neg(frame2) << ")"
         << " CPU = ";
    const ios::fmtflags& old_flags = cout.flags();
    const int width = 12;
    const int old_prec = cout.precision();
    const int new_prec = 2;
    cout.setf(ios_base::fixed);
    cout.setf(ios_base::showpoint);
    cout << setprecision(new_prec)
         << setw(width) << setup_cpu_time + mult_cpu_time << " ms: "
         << setw(width) << setup_cpu_time << " (setup) + "
         << setw(width) << mult_cpu_time  << " (mult). "
         << setprecision(old_prec)
         << endl;
    cout.flags(old_flags);
  }

  template< typename Multivector_T >
  void
  time_mult(Multivector_T& a, const Multivector_T& b,
            const index_set_t& inner_frame,
            const index_set_t& outer_frame)
  {
    const int min_trials_if_zero = 1000;
    Multivector_T c;
    clock_t cpu_time = clock();
     c = a * b*(Scalar_T(1.0*rand())/RAND_MAX) + a*(Scalar_T(1.0*rand())/RAND_MAX);
    double setup_cpu_time = elapsed(cpu_time);
    if (setup_cpu_time < MS_PER_S)
    {
      const int max_trials =
       setup_cpu_time == 0.0 ?
         min_trials_if_zero :
         int(floor(MS_PER_S/setup_cpu_time));
      clock_t cpu_time = clock();
      int trials;
      for (trials = 0; trials != max_trials; ++trials)
        c = a * b*(Scalar_T(1.0*rand())/RAND_MAX) + a*(Scalar_T(1.0*rand())/RAND_MAX);
      setup_cpu_time = elapsed(cpu_time)/trials;
    }
    a = c;
    cpu_time = clock();
     c = a * a;
    double mult_cpu_time = elapsed(cpu_time);
    if (mult_cpu_time < MS_PER_S)
    {
      const int max_trials =
       mult_cpu_time == 0.0 ?
         min_trials_if_zero :
         int(floor(MS_PER_S/mult_cpu_time));
      clock_t cpu_time = clock();
      int trials;
      for (trials = 0; trials != max_trials; ++trials)
        c = a * a;
      mult_cpu_time = elapsed(cpu_time)/trials;
    }
    print_times(inner_frame, outer_frame, setup_cpu_time, mult_cpu_time);
  }

  template< class Multivector_T >
  void
  mult_test(const index_t n, const index_t max_n)
  {
    cout << "Clifford algebra squaring test:" << endl;

    typedef Multivector_T multivector_t;
    typedef index_set_t e_;
    static const index_t v_lo = e_::v_lo;
    static const index_t v_hi = e_::v_hi;
    const index_t max_index = min(n, max_n);
    const e_ outer_frame = index_range<v_lo,v_hi>(-max_index, max_index);
    
    clock_t cpu_time = clock();
    multivector_t a = multivector_t(outer_frame, 1);
    const double setup_cpu_time = elapsed(cpu_time);
    
    cout << "Square of unit pseudoscalar (" << a << ")" << endl;
    
    cpu_time = clock();
     a *= a;
    const double mult_cpu_time = elapsed(cpu_time);
    
    print_times(outer_frame, outer_frame, setup_cpu_time, mult_cpu_time);
    cout << "Square is " << a << endl;
    e_ inner_frame;
    a = 1;
    inner_frame = e_();
    for (index_t i = 1; i != max_index+1; i++)
    {
      inner_frame |= e_(i);
      time_mult(a, multivector_t(e_(i) , 1.0), inner_frame, inner_frame);
      inner_frame |= e_(-i);
      time_mult(a, multivector_t(e_(-i), 1.0), inner_frame, inner_frame);
    }
  }
}

int squaring(const int n)
{
  using namespace glucat_mult_test;
  if (n > max_n)
  {
    cout << "Value " << n << " is too big." << endl;
    cout << "Maximum value allowed is " << max_n << "." << endl;
    return 1;
  }
  cout << "framed_multi<double>" << endl;
  mult_test< framed_multi<double> >(n, max_n);
  cout << "matrix_multi<double>" << endl;
  mult_test< matrix_multi<double> >(n, max_n);

  return 0;
}
#endif // GLUCAT_TEST_SQUARING_H
