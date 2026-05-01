#ifndef _GLUCAT_TEST_DOCTEST_TEST_DOCTEST_H
#define _GLUCAT_TEST_DOCTEST_TEST_DOCTEST_H

#include <doctest.h>
#include "glucat/glucat.h"
#include <iostream>
#include <fstream>

namespace glucat {

  template< typename Multivector_T >
  bool is_error(const Multivector_T& lhs, const Multivector_T& rhs, typename Multivector_T::scalar_t tol)
  {
    return ( (abs(lhs) < tol) || (abs(rhs) < tol) )
           ? (abs(lhs - rhs) > tol)
           : (abs(lhs - rhs) > abs(rhs) * tol);
  }

  template< typename Multivector_T >
  void test_idents_templated(Multivector_T& a, Multivector_T& b, Multivector_T& c)
  {
    typedef typename Multivector_T::scalar_t scalar_t;
    static const scalar_t eps = std::numeric_limits<scalar_t>::epsilon();
    const scalar_t tol = eps * 512.0;

    Multivector_T lhs;
    Multivector_T rhs;

    const index_t a_grade = a.frame().count();
    const index_t b_grade = b.frame().count();

    SUBCASE("Identity [HS] (1.21a)") {
      for (index_t r = 1; r <= a_grade; ++r) {
        const Multivector_T a_r = a(r);
        for (index_t s = 1; s <= b_grade; ++s) {
          const Multivector_T b_s = b(s);
          lhs = a_r & b_s;
          rhs = (a_r * b_s)(index_t(std::abs(int(r-s))));
          CHECK_FALSE(is_error(lhs, rhs, tol));
        }
      }
    }

    SUBCASE("Identity [HS] (1.22a)") {
      for (index_t r = 0; r <= a_grade; ++r) {
        const Multivector_T a_r = a(r);
        for (index_t s = 0; s <= b_grade; ++s) {
          const Multivector_T b_s = b(s);
          lhs = a_r ^ b_s;
          rhs = (a_r * b_s)(r+s);
          CHECK_FALSE(is_error(lhs, rhs, tol));
        }
      }
    }

    SUBCASE("Identity [HS] (1.25a)") {
      lhs = (a ^ b) ^ c;
      rhs = a ^ (b ^ c);
      CHECK_FALSE(is_error(lhs, rhs, tol));
    }

    SUBCASE("Identity [HS] (1.31)") {
      const Multivector_T a_1 = a(1);
      lhs = a_1 * b;
      rhs = (a_1 & b) + (a_1 ^ b);
      CHECK_FALSE(is_error(lhs, rhs, tol));
    }
  }

  template< typename Multivector_T >
  void test_transcendental_templated(const Multivector_T& A)
  {
    using scalar_t = typename Multivector_T::scalar_t;
    static const scalar_t eps = std::numeric_limits<scalar_t>::epsilon();
    const scalar_t tol = eps * 1024.0;

    const Multivector_T exp_A = exp(A);
    const Multivector_T exp_mA = exp(-A);
    const Multivector_T cos_A = cos(A);
    const Multivector_T sin_A = sin(A);
    const Multivector_T cosh_A = cosh(A);
    const Multivector_T sinh_A = sinh(A);
    const Multivector_T sqrt_A = sqrt(A);

    SUBCASE("exp(A)*exp(-A) == 1") {
      CHECK_FALSE(is_error(exp_A * exp_mA, Multivector_T(1), tol));
    }

    SUBCASE("cosh(A)+sinh(A) == exp(A)") {
      CHECK_FALSE(is_error(cosh_A + sinh_A, exp_A, tol));
    }

    SUBCASE("sqrt(A)*sqrt(A) == A") {
      if (!A.isnan() && !A.isinf()) {
        CHECK_FALSE(is_error(sqrt_A * sqrt_A, A, tol));
      }
    }
  }

  template< typename Multivector_T >
  void run_peg11_test()
  {
    test_transcendental_templated(Multivector_T(0));
    test_transcendental_templated(Multivector_T(1));
    test_transcendental_templated(Multivector_T::random(typename Multivector_T::index_set_t(1)));
  }

  template< typename Multivector_T >
  void run_peg00_test(index_t max_index)
  {
    using index_set_t = typename Multivector_T::index_set_t;
    using scalar_t = typename Multivector_T::scalar_t;
    
    scalar_t fill = 0.5;
    index_set_t frm = index_set_t();
    for (index_t i = 1; i != max_index+1; i++)
    {
      frm |= index_set_t(i);
      Multivector_T a = Multivector_T::random(frm, fill);
      Multivector_T b = Multivector_T::random(frm, fill);
      Multivector_T c = Multivector_T::random(frm, fill);
      test_idents_templated(a, b, c);
      
      frm |= index_set_t(-i);
      a = Multivector_T::random(frm, fill);
      b = Multivector_T::random(frm, fill);
      c = Multivector_T::random(frm, fill);
      test_idents_templated(a, b, c);
    }
  }

} // namespace glucat

#endif // _GLUCAT_TEST_DOCTEST_TEST_DOCTEST_H
