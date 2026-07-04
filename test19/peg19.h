#ifndef GLUCAT_TEST_PEG19_H
#define GLUCAT_TEST_PEG19_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg19.h : equivalence test for sandwich products and transforms
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

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "glucat/glucat.h"

using namespace std;
using namespace glucat;

template<typename Multivector_T>
bool check_mv(const string& name, const Multivector_T& A, const Multivector_T& B, typename Multivector_T::scalar_t tol) {
    auto diff = A - B;
    auto n = diff.norm();
    bool passed = (n < tol);
    if (passed) {
        cout << "  " << name << ": OK" << endl;
    } else {
        cout << "  " << name << ": FAILED (norm of diff: " << n << ")" << endl;
    }
    return passed;
}

template<typename Scalar_T>
void versor_equivalence_test() {
    using mm_t = matrix_multi<Scalar_T, -16, 16>;
    using index_set_t = typename mm_t::index_set_t;
    static const Scalar_T eps = std::numeric_limits<Scalar_T>::epsilon();
    const Scalar_T tol = eps * 4096.0;

    index_set_t frm;
    frm.set(1);
    frm.set(2);
    frm.set(3);

    mm_t X = mm_t::random(frm);

    // Bivector generator
    mm_t A = mm_t::random(frm)(2);
    mm_t R = exp(A);

    cout << "Testing bivector generator equivalence on Cl(3,0):" << endl;

    // Explicit sandwich: R * X * involute(exp(-A))
    mm_t explicit_sand = R * X * involute(exp(-A));
    // Solver-based operator|: X | R
    mm_t operator_or = X | R;
    // versor(): X.versor(R)
    mm_t versor_val = X.versor(R);
    // versor_exp(): X.versor_exp(A)
    mm_t versor_exp_val = X.versor_exp(A);

    bool all_passed = true;
    all_passed &= check_mv("operator|  == explicit", operator_or, explicit_sand, tol);
    all_passed &= check_mv("versor  == explicit", versor_val, explicit_sand, tol);
    all_passed &= check_mv("versor_exp == explicit", versor_exp_val, explicit_sand, tol);
    all_passed &= check_mv("versor_exp(prechecked) == explicit", X.versor_exp(A, true), explicit_sand, tol);

    if (all_passed) {
        cout << "All equivalence tests PASSED" << endl;
    } else {
        cout << "Some equivalence tests FAILED" << endl;
    }
}

int test19();

#endif
