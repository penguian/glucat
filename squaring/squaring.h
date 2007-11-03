#ifndef GLUCAT_TEST_SQUARING_H
#define GLUCAT_TEST_SQUARING_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    squaring.h : multiplication speed test
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
#include "squaring/tuning.h"
#include "glucat/glucat_imp.h"
#include "test/try_catch.h"
#include <stdio.h>
#include <iomanip>

namespace glucat_mult_test
{
  using namespace glucat;
  using namespace std;

  const index_t max_n = DEFAULT_HI;

  inline
  static
  double
  elapsed( clock_t cpu_time )
  { return ((clock() - cpu_time)*MS_PER_S) / CLOCKS_PER_SEC; }

  template < typename Index_Set_T >
  inline
  void
  print_times(const Index_Set_T& frame1, const Index_Set_T& frame2,
              const double setup_cpu_time, 
              const double prod_cpu_time)
  {
    const int index_width = 2;
    cout << "Cl(" << setw(index_width) <<  max_pos(frame1) << ","
                  << setw(index_width) << -min_neg(frame1) << ") in "
         << "Cl(" << setw(index_width) <<  max_pos(frame2) << ","
                  << setw(index_width) << -min_neg(frame2) << ")"
         << " CPU = ";
    const ios::fmtflags& old_flags = cout.flags();
    const streamsize width = 12;
    const streamsize old_prec = cout.precision();
    const streamsize new_prec = 3;
    cout.setf(ios_base::fixed);
    cout.setf(ios_base::showpoint);
    cout << setprecision(new_prec)
         << setw(width) << setup_cpu_time << " ms (setup) "
         << setw(width) << prod_cpu_time  << " (*) "
         << setprecision(old_prec)
         << endl;
    cout.flags(old_flags);
  }

#ifdef GLUCAT_TEST_REPEAT
  static
  int
  enough_trials(double cpu_time)
  {
    const int min_trials_if_zero = 10000;
    return (cpu_time == 0.0)
        ? min_trials_if_zero 
        : int(floor(MS_PER_S/cpu_time));
  }
#endif

  template< typename Multivector_T >
  static
  void
  time_mult(Multivector_T& a, const Multivector_T& b,
            const typename Multivector_T::index_set_t& inner_frame,
            const typename Multivector_T::index_set_t& outer_frame)
  {
    typedef Multivector_T multivector_t;
    typedef typename multivector_t::index_set_t index_set_t;
    typedef typename multivector_t::scalar_t scalar_t;

    const scalar_t RAND_SCALE = 1.0/RAND_MAX;
    const scalar_t rand1 = rand()*RAND_SCALE;

    Multivector_T c;
    clock_t cpu_time = clock();
     c = a * (b * rand1 + scalar_t(1.0));
    double setup_cpu_time = elapsed(cpu_time);
#ifdef GLUCAT_TEST_REPEAT
    if (setup_cpu_time < MS_PER_S)
    {
      const int max_trials = enough_trials(setup_cpu_time);
      cpu_time = clock();
      int trials;
      for (trials = 0; trials != max_trials; ++trials)
        c = a * (b * rand1 + scalar_t(1.0));
      setup_cpu_time = elapsed(cpu_time)/trials;
    }
#endif
    a = c;
    cpu_time = clock();
     c = a * a;
    double prod_cpu_time = elapsed(cpu_time);
#ifdef GLUCAT_TEST_REPEAT
    if (prod_cpu_time < MS_PER_S)
    {
      const int max_trials = enough_trials(prod_cpu_time);
      cpu_time = clock();
      int trials;
      for (trials = 0; trials != max_trials; ++trials)
        c = a * a;
      prod_cpu_time = elapsed(cpu_time)/trials;
    }
#endif
    print_times(inner_frame, outer_frame, setup_cpu_time, prod_cpu_time);
  }

  template< class Multivector_T >
  static
  void
  mult_test(const index_t n, const index_t max_n)
  {
    typedef Multivector_T multivector_t;
    typedef typename multivector_t::index_set_t e_;
    static const index_t v_lo = e_::v_lo;
    static const index_t v_hi = e_::v_hi;
    static const index_t index_lim = std::min(-v_lo, v_hi);
    if (n > index_lim)
    {
      cout << "Value " << n << " is too big." << endl;
      cout << "Maximum value possible is " << index_lim << "." << endl;
      return;
    }
    cout << "Clifford algebra squaring test:" << endl;
    const index_t max_index = min(n, max_n);
    const e_ outer_frame = e_(make_pair(-max_index, max_index));

    multivector_t a = multivector_t(outer_frame, 1);

    cout << a << " * " << a << " == " << (a * a) << endl;
    cout << endl;

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

int squaring(const int n);

#endif // GLUCAT_TEST_SQUARING_H
