#ifndef GLUCAT_TEST_PEG18_H
#define GLUCAT_TEST_PEG18_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg18.h : consistency test between framed_multi and matrix_multi
                             -------------------
    begin                : Wed May 06 2026
    copyright            : (C) 2001-2026 by Paul C. Leopardi
 ***************************************************************************/

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "glucat/glucat.h"

using namespace std;
using namespace glucat;

template<typename fm_t, typename mm_t>
bool check_scalar(const string& name, typename mm_t::scalar_t A, typename fm_t::scalar_t a) {
    auto diff = std::abs(A - a);
    if (diff < 1e-12) {
        cout << "  " << name << ": OK" << endl;
        return true;
    } else {
        cout << "  " << name << ": FAILED (val1: " << A << ", val2: " << a << ", diff: " << diff << ")" << endl;
        return false;
    }
}

template<typename fm_t, typename mm_t>
bool check_mv(const string& name, const mm_t& A, const fm_t& a) {
    mm_t B(a, A.frame());
    auto diff = A - B;
    auto n = diff.norm();
    if (n < 1e-12) {
        cout << "  " << name << ": OK" << endl;
        return true;
    } else {
        cout << "  " << name << ": FAILED (norm of diff: " << n << ")" << endl;
        return false;
    }
}

template<typename Scalar_T>
void consistency_test() {
    using fm_t = framed_multi<Scalar_T, -16, 16>;
    using mm_t = matrix_multi<Scalar_T, -16, 16>;
    using index_set_t = typename fm_t::index_set_t;

    std::vector<std::pair<int, int>> signatures;
    // (0,2), (1,1), (2,0)
    for (int p=0; p<=2; ++p) signatures.push_back({p, 2-p});
    // (4,0) to (0,4)
    for (int p=4; p>=0; --p) signatures.push_back({p, 4-p});
    // (8,0) to (0,8)
    for (int p=8; p>=0; --p) signatures.push_back({p, 8-p});

    for (const auto& sig : signatures) {
        cout << "Testing signature (" << sig.first << "," << sig.second << ")" << endl;
        index_set_t frm;
        for (int i=1; i<=sig.first; ++i) frm.set(i);
        for (int i=1; i<=sig.second; ++i) frm.set(-i);

        fm_t a = fm_t::random(frm);
        fm_t b = fm_t::random(frm);
        mm_t A(a, frm);
        mm_t B(b, frm);

        // Scalar-yielding operators
        check_scalar<fm_t, mm_t>("norm()", A.norm(), a.norm());
        check_scalar<fm_t, mm_t>("quad()", A.quad(), a.quad());
        check_scalar<fm_t, mm_t>("max_abs()",  A.max_abs(),  a.max_abs());
        check_scalar<fm_t, mm_t>("scalar()", A.scalar(), a.scalar());
        check_scalar<fm_t, mm_t>("star()", star(A, B), star(a, b));
        check_scalar<fm_t, mm_t>("hstar()", hstar(A, B), hstar(a, b));

        // Multivector-yielding operators (Unary)
        check_mv<fm_t, mm_t>("operator-", -A, -a);
        check_mv<fm_t, mm_t>("inv()",      A.inv(), a.inv());
        check_mv<fm_t, mm_t>("involute()", A.involute(), a.involute());
        check_mv<fm_t, mm_t>("reverse()",  A.reverse(), a.reverse());
        check_mv<fm_t, mm_t>("conj()",     A.conj(), a.conj());

        // Multivector-yielding operators (Binary)
        check_mv<fm_t, mm_t>("operator+", A + B, a + b);
        check_mv<fm_t, mm_t>("operator-", A - B, a - b);
        check_mv<fm_t, mm_t>("operator*", A * B, a * b);
        check_mv<fm_t, mm_t>("operator/", A / B, a / b);
        check_mv<fm_t, mm_t>("operator^", A ^ B, a ^ b);
        check_mv<fm_t, mm_t>("operator&", A & B, a & b);
        check_mv<fm_t, mm_t>("operator%", A % B, a % b);
        check_mv<fm_t, mm_t>("operator|", A | B, a | b);
    }
}

int test18();

#endif
