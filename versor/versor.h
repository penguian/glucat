#ifndef GLUCAT_VERSOR_H
#define GLUCAT_VERSOR_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    versor.h : Versor product and timing benchmark
                             -------------------
    begin                : Thu Jun 11 2026
    copyright            : (C) 2001-2026 by Paul C. Leopardi
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

 ***************************************************************************/

#define _GLUCAT_TEST_REPEAT
#include "test/timing.h"
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <vector>

namespace glucat_versor_test
{
  using namespace glucat;
  using namespace glucat::timing;
  using namespace std;

  const index_t max_n = DEFAULT_HI;
  typedef index_set< -max_n, max_n > index_set_t;

  template<typename Multivector_T>
  void time_versor(const index_set_t& frame)
  {
    // Create random multivector
    Multivector_T X = Multivector_T::random(frame);

    // Create a random bivector generator A and its exponentiated rotor R
    Multivector_T A = Multivector_T::random(frame)(2);
    Multivector_T R = exp(A);

    // 0. Time naive sandwich product (slow solver path)
    clock_t cpu_time = clock();
    Multivector_T result_sand = R * X / R.involute();
    double time_sand = elapsed(cpu_time);
    #ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; time_sand == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
      for (int trials = 0; trials != nbr_trials; ++trials)
        result_sand = R * X / R.involute();
      time_sand = elapsed(cpu_time) / nbr_trials;
    }
    #endif

    // 1. Time operator| (slow solver path)
    cpu_time = clock();
    Multivector_T result_vert = X | R;
    double time_vert = elapsed(cpu_time);
    #ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; time_vert == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
      for (int trials = 0; trials != nbr_trials; ++trials)
        result_vert = X | R;
      time_vert = elapsed(cpu_time) / nbr_trials;
    }
    #endif

    // 2. Time versor (fast reversion path)
    cpu_time = clock();
    Multivector_T result_versor = X.versor(R);
    double time_trans = elapsed(cpu_time);
#ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; time_trans == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
      for (int trials = 0; trials != nbr_trials; ++trials)
        result_versor = X.versor(R);
      time_trans = elapsed(cpu_time) / nbr_trials;
    }
#endif

    // 3. Time versor_exp (fast reversion path directly with bivector generator)
    cpu_time = clock();
    Multivector_T result_versor_exp = X.versor_exp(A);
    double time_trans_exp = elapsed(cpu_time);
#ifdef _GLUCAT_TEST_REPEAT
    for (int nbr_trials = EXTRA_TRIALS; time_trans_exp == 0.0; nbr_trials *= EXTRA_TRIALS)
    {
      cpu_time = clock();
      for (int trials = 0; trials != nbr_trials; ++trials)
        result_versor_exp = X.versor_exp(A);
      time_trans_exp = elapsed(cpu_time) / nbr_trials;
    }
#endif

    // Print results in a nice format
    const int width = 10;
    cout << "Cl(" << setw(2) << max_pos(frame) << "," << setw(2) << -min_neg(frame) << ") | "
         << "naive : " << setw(width) << fixed << setprecision(3) << time_sand << " ms | "
         << "operator| : " << setw(width) << time_vert << " ms | "
         << "versor : " << setw(width) << time_trans << " ms | "
         << "versor_exp : " << setw(width) << time_trans_exp << " ms" << endl;
  }

  template<typename Multivector_T>
  void run_typed_benchmark(const int n)
  {
    for (int d = 1; d <= n; ++d)
    {
      index_set_t frame;
      for (int i = 1; i <= d; ++i)
        frame.set(i);
      time_versor<Multivector_T>(frame);
    }
    for (int d = 1; d <= n; ++d)
    {
      index_set_t frame;
      for (int i = -1; i >= -d; --i)
        frame.set(i);
      time_versor<Multivector_T>(frame);
    }
  }

  template<typename Scalar_T>
  void run_benchmark(const int n)
  {
    cout << "\n--- Benchmarking framed_multi ---" << endl;
    run_typed_benchmark< framed_multi<Scalar_T, -max_n, max_n> >(n);

    cout << "\n--- Benchmarking matrix_multi ---" << endl;
    run_typed_benchmark< matrix_multi<Scalar_T, -max_n, max_n> >(n);
  }
}

int versor_test(const int n);

#endif // GLUCAT_VERSOR_H
