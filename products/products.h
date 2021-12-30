#ifndef GLUCAT_TEST_PRODUCTS_H
#define GLUCAT_TEST_PRODUCTS_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    products.h : products speed test
                             -------------------
    begin                : Sun 2007-09-01
    copyright            : (C) 2012-2021 by Paul C. Leopardi
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
#include "test/timing.h"
#include <stdio.h>
#include <iomanip>

namespace glucat_prod_test
{
  using namespace glucat;
  using namespace glucat::timing;
  using namespace std;

  const index_t max_n = DEFAULT_HI;

  template < typename Index_Set_T >
  inline
  void
  print_times(const Index_Set_T& frame1, const Index_Set_T& frame2,
              const double setup_cpu_time,
              const double prod_cpu_time,
              const double outer_cpu_time,
              const double inner_cpu_time,
              const double contract_cpu_time)
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
         << setw(width) << outer_cpu_time << " (^) "
         << setw(width) << inner_cpu_time << " (&) "
         << setw(width) << contract_cpu_time  << " (%) "
         << setprecision(old_prec)
         << endl;
    cout.flags(old_flags);
  }

  template< typename Multivector_T >
  static
  void
  time_mult(const typename Multivector_T::index_set_t inner_frame,
            const typename Multivector_T::index_set_t outer_frame,
            const typename Multivector_T::scalar_t fill)
  {
    typedef Multivector_T multivector_t;

    clock_t cpu_time = clock();
      multivector_t a = multivector_t(multivector_t::random(inner_frame, fill), outer_frame, true);
    double setup_cpu_time = elapsed(cpu_time);
#ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; setup_cpu_time == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
        for (int trials = 0; trials != nbr_trials; ++trials)
          a = multivector_t(multivector_t::random(inner_frame, fill), outer_frame, true);
      setup_cpu_time = elapsed(cpu_time) / nbr_trials;
    }
#endif
    multivector_t c;
    cpu_time = clock();
      c = a * a;
    double prod_cpu_time = elapsed(cpu_time);
#ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; prod_cpu_time == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
        for (int trials = 0; trials != nbr_trials; ++trials)
          c = a * a;
      prod_cpu_time = elapsed(cpu_time) / nbr_trials;
    }
#endif
    cpu_time = clock();
      c = a ^ a;
    double outer_cpu_time = elapsed(cpu_time);
#ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; outer_cpu_time == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
      for (int trials = 0; trials != nbr_trials; ++trials)
        c = a ^ a;
      outer_cpu_time = elapsed(cpu_time) / nbr_trials;
    }
#endif
    cpu_time = clock();
      c = a & a;
    double inner_cpu_time = elapsed(cpu_time);
#ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; inner_cpu_time == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
        for (int trials = 0; trials != nbr_trials; ++trials)
          c = a & a;
      inner_cpu_time = elapsed(cpu_time) / nbr_trials;
    }
#endif
    cpu_time = clock();
      c = a % a;
    double contract_cpu_time = elapsed(cpu_time);
#ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; contract_cpu_time == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
        for (int trials = 0; trials != nbr_trials; ++trials)
          c = a % a;
      contract_cpu_time = elapsed(cpu_time) / nbr_trials;
    }
#endif
    print_times(inner_frame, outer_frame,
                setup_cpu_time, prod_cpu_time,
                outer_cpu_time, inner_cpu_time, contract_cpu_time);
  }

  template< class Multivector_T >
  static
  void
  prod_test_with_fill(const index_t max_index, const typename Multivector_T::scalar_t fill)
  {
    cout << "Fill: " << fill << endl;

    typedef Multivector_T multivector_t;
    typedef typename multivector_t::index_set_t index_set_t;

    const index_set_t outer_frame = index_set_t(make_pair(-max_index, max_index));
    index_set_t inner_frame = index_set_t();
    for (index_t i = 1; i != max_index+1; i++)
    {
      inner_frame |= index_set_t(i);
      time_mult<multivector_t>(inner_frame, inner_frame, fill);
      inner_frame |= index_set_t(-i);
      time_mult<multivector_t>(inner_frame, inner_frame, fill);
    }
    inner_frame = index_set_t();
    for (index_t i = 1; i != max_index+1; i++)
    {
      inner_frame |= index_set_t(i);
      time_mult<multivector_t>(inner_frame, outer_frame, fill);
      inner_frame |= index_set_t(-i);
      time_mult<multivector_t>(inner_frame, outer_frame, fill);
    }
  }

  template< class Multivector_T >
  static
  void
  prod_test(const index_t n, const index_t max_n)
  {
    static const index_t v_lo = Multivector_T::v_lo;
    static const index_t v_hi = Multivector_T::v_hi;
    static const index_t index_lim = std::min(-v_lo, v_hi);
    if (n > index_lim)
    {
      cout << "Value " << n << " is too big." << endl;
      cout << "Maximum value possible is " << index_lim << "." << endl;
      return;
    }
    cout << "Clifford algebra products test:" << endl;
    const index_t max_index = min(n, max_n);

    prod_test_with_fill<Multivector_T>(max_index, 0.5);
    prod_test_with_fill<Multivector_T>(max_index, 1.0);
  }
}

int products(const int n);

#endif // GLUCAT_TEST_PRODUCTS_H
