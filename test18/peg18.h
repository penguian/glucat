#ifndef GLUCAT_TEST_PEG18_H
#define GLUCAT_TEST_PEG18_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg18.h : consistency test between framed_multi and matrix_multi
                             -------------------
    begin                : Wed May 06 2026
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

template<typename fm_t, typename mm_t>
bool check_scalar(const string& name, typename mm_t::scalar_t A, typename fm_t::scalar_t a, bool verbose) {
    auto diff = std::abs(A - a);
    bool passed = (diff < 1e-12);
    if (verbose && passed) {
        cout << "  " << name << ": OK" << endl;
    } else if (!passed) {
        if (!verbose) cout << endl;
        cout << "  " << name << ": FAILED (val1: " << A << ", val2: " << a << ", diff: " << diff << ")" << endl;
    }
    return passed;
}

template<typename fm_t, typename mm_t>
bool check_mv(const string& name, const mm_t& A, const fm_t& a, bool verbose) {
    mm_t B(a, A.frame());
    auto diff = A - B;
    auto n = diff.norm();
    bool passed = (n < 1e-12);
    if (verbose && passed) {
        cout << "  " << name << ": OK" << endl;
    } else if (!passed) {
        if (!verbose) cout << endl;
        cout << "  " << name << ": FAILED (norm of diff: " << n << ")" << endl;
    }
    return passed;
}

template<typename Scalar_T>
void consistency_test() {
    using fm_t = framed_multi<Scalar_T, -16, 16>;
    using mm_t = matrix_multi<Scalar_T, -16, 16>;
    using index_set_t = typename fm_t::index_set_t;

    bool verbose = control_t::verbose();
    std::vector<std::pair<int, int>> signatures;
    // (0,2), (1,1), (2,0)
    for (int p=0; p<=2; ++p) signatures.push_back({p, 2-p});
    // (4,0) to (0,4)
    for (int p=4; p>=0; --p) signatures.push_back({p, 4-p});
    // (8,0) to (0,8)
    for (int p=8; p>=0; --p) signatures.push_back({p, 8-p});

    for (const auto& sig : signatures) {
        if (verbose) {
          cout << "Testing signature (" << sig.first << "," << sig.second << ")" << endl;
        } else {
          cout << "Testing signature (" << sig.first << "," << sig.second << ") : " << flush;
        }
        index_set_t frm;
        for (int i=1; i<=sig.first; ++i) frm.set(i);
        for (int i=1; i<=sig.second; ++i) frm.set(-i);

        fm_t a = fm_t::random(frm);
        fm_t b = fm_t::random(frm);
        mm_t A(a, frm);
        mm_t B(b, frm);

        bool all_passed = true;
        // Scalar-yielding operators
        all_passed &= check_scalar<fm_t, mm_t>("norm()", A.norm(), a.norm(), verbose);
        all_passed &= check_scalar<fm_t, mm_t>("quad()", A.quad(), a.quad(), verbose);
        all_passed &= check_scalar<fm_t, mm_t>("max_abs()",  A.max_abs(),  a.max_abs(), verbose);
        all_passed &= check_scalar<fm_t, mm_t>("scalar()", A.scalar(), a.scalar(), verbose);
        all_passed &= check_scalar<fm_t, mm_t>("star()", star(A, B), star(a, b), verbose);
        all_passed &= check_scalar<fm_t, mm_t>("hstar()", hstar(A, B), hstar(a, b), verbose);

        // Multivector-yielding operators (Unary)
        all_passed &= check_mv<fm_t, mm_t>("operator-", -A, -a, verbose);
        all_passed &= check_mv<fm_t, mm_t>("inv()",      A.inv(), a.inv(), verbose);
        all_passed &= check_mv<fm_t, mm_t>("involute()", A.involute(), a.involute(), verbose);
        all_passed &= check_mv<fm_t, mm_t>("reverse()",  A.reverse(), a.reverse(), verbose);
        all_passed &= check_mv<fm_t, mm_t>("conj()",     A.conj(), a.conj(), verbose);

        // Multivector-yielding operators (Binary)
        all_passed &= check_mv<fm_t, mm_t>("operator+", A + B, a + b, verbose);
        all_passed &= check_mv<fm_t, mm_t>("operator-", A - B, a - b, verbose);
        all_passed &= check_mv<fm_t, mm_t>("operator*", A * B, a * b, verbose);
        all_passed &= check_mv<fm_t, mm_t>("operator/", A / B, a / b, verbose);
        all_passed &= check_mv<fm_t, mm_t>("operator^", A ^ B, a ^ b, verbose);
        all_passed &= check_mv<fm_t, mm_t>("operator&", A & B, a & b, verbose);
        all_passed &= check_mv<fm_t, mm_t>("operator%", A % B, a % b, verbose);
        all_passed &= check_mv<fm_t, mm_t>("operator|", A | B, a | b, verbose);

        if (!verbose) {
          if (all_passed) cout << "PASSED" << endl;
          else cout << "FAILED" << endl;
        }
    }
}

int test18();

#endif
