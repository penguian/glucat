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

#include <iomanip>

namespace glucat_mult_test
{
  using namespace glucat;
  const int width = 12;
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
    cout << "Cl(" <<  max_pos(frame1) << ","
                  << -min_neg(frame1) << ") in "
         << "Cl(" <<  max_pos(frame2) << ","
                  << -min_neg(frame2) << ")"
         << " CPU = "
         << setw(width) << setup_cpu_time + mult_cpu_time << " ms: "
         << setw(width) << setup_cpu_time   << " (setup) + "
         << setw(width) << mult_cpu_time    << " (mult). "
         << endl;
  }

  template< typename Multivector_T >
  void
  time_mult(Multivector_T& a, const Multivector_T& b,
            const index_t i, const index_t j,
            const index_set_t& inner_frame,
            const index_set_t& outer_frame)
  {
    Multivector_T c;
    clock_t cpu_time = clock();
     a = a * b*(Scalar_T(1.0*rand())/RAND_MAX) + a*(Scalar_T(1.0*rand())/RAND_MAX);
    const double setup_cpu_time = elapsed(cpu_time);
    cpu_time = clock();
     c = a * a;
    const double mult_cpu_time = elapsed(cpu_time);
    print_times(inner_frame, outer_frame, setup_cpu_time, mult_cpu_time);
  }

  template< class Multivector_T >
  void
  mult_test(const index_t n, const index_t max_n)
  {
    cout << "Clifford algebra squaring test:" << endl;

    typedef Multivector_T m_;
    typedef index_set_t e_;
    static const index_t v_lo = e_::v_lo;
    static const index_t v_hi = e_::v_hi;

    m_ a;
    const index_t max_pos = min(index_t(2*n), max_n);
    const e_ pos_frame = index_range<v_lo,v_hi>(1, max_pos);

    const index_t max_index = min(n, max_n);
    const e_ outer_frame = index_range<v_lo,v_hi>(-max_index, max_index);
    clock_t cpu_time = clock();
     a = m_(outer_frame, 1);
    const double setup_cpu_time = elapsed(cpu_time);
    cout << "Square of unit pseudoscalar (" << a << ")" << endl;
    cpu_time = clock();
     a *= a;
    const double mult_cpu_time = elapsed(cpu_time);
    print_times(outer_frame, outer_frame, setup_cpu_time, mult_cpu_time);
    cout << "Square is " << a << endl;
    a = m_(1, pos_frame);
    e_ inner_frame;
    for (index_t i = 1; i != max_pos+1; i++)
    {
      inner_frame |= e_(i);
      time_mult(a, m_(e_(i) , 1.0, pos_frame), i, 0,   inner_frame, pos_frame);
    }
    inner_frame = e_();
    a = m_(1, outer_frame);
    for (index_t i = 1; i != max_index+1; i++)
    {
      inner_frame |= e_(i);
      time_mult(a, m_(e_(i) , 1.0, outer_frame), i, i-1, inner_frame, outer_frame);
      inner_frame |= e_(-i);
      time_mult(a, m_(e_(-i), 1.0, outer_frame), i, i,   inner_frame, outer_frame);
    }
    a = 1;
    inner_frame = e_();
    for (index_t i = 1; i != max_index+1; i++)
    {
      inner_frame |= e_(i);
      time_mult(a, m_(e_(i) , 1.0), i, i-1, inner_frame, inner_frame);
      inner_frame |= e_(-i);
      time_mult(a, m_(e_(-i), 1.0), i, i,   inner_frame, inner_frame);
    }
  }
}

int squaring(const long int n)
{
  using namespace glucat_mult_test;

  cout << "framed_multi<double>" << endl;
  mult_test< framed_multi<double> >(n, max_n);
  cout << "matrix_multi<double>" << endl;
  mult_test< matrix_multi<double> >(n, max_n);

  return 0;
}
#endif // GLUCAT_TEST_SQUARING_H
