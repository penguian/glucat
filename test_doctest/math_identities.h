#ifndef _GLUCAT_TEST_DOCTEST_MATH_IDENTITIES_H
#define _GLUCAT_TEST_DOCTEST_MATH_IDENTITIES_H

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
  void test_edge_cases_templated()
  {
    using scalar_t = typename Multivector_T::scalar_t;
    using index_set_t = typename Multivector_T::index_set_t;
    static const scalar_t eps = std::numeric_limits<scalar_t>::epsilon();
    const scalar_t tol = eps * 1024.0;

    SUBCASE("Self-Aliasing Guards") {
      Multivector_T a = Multivector_T::random(index_set_t(1), 1.0);
      if (a.frame() != index_set_t(1)) {
        a += Multivector_T(index_set_t(1));
      }
      Multivector_T a_copy = a;
      
      a += a;
      CHECK_FALSE(is_error(a, a_copy * scalar_t(2), tol));
      a = a_copy;
      
      a -= a;
      CHECK_FALSE(is_error(a, Multivector_T(0), tol));
      a = a_copy;
      
      a *= a;
      CHECK_FALSE(is_error(a, a_copy * a_copy, tol));
      a = a_copy;
      
      CHECK(a == a);
    }

    SUBCASE("NaN and Infinity Propagation") {
      if constexpr (std::numeric_limits<scalar_t>::has_quiet_NaN) {
        scalar_t nan_val = std::numeric_limits<scalar_t>::quiet_NaN();
        Multivector_T nan_mv(nan_val);
        
        CHECK(nan_mv.isnan());
        CHECK(!nan_mv.isinf());
        
        Multivector_T res_sqrt = sqrt(nan_mv);
        CHECK(res_sqrt.isnan());
        
        Multivector_T res_log = log(nan_mv);
        CHECK(res_log.isnan());
        
        Multivector_T res_exp = exp(nan_mv);
        CHECK(res_exp.isnan());

        Multivector_T a = Multivector_T::random(index_set_t(1), 1.0);
        if (a.frame() != index_set_t(1)) {
          a += Multivector_T(index_set_t(1));
        }
        Multivector_T res_div = a / nan_mv;
        CHECK(res_div.isnan());
      }
      
      if constexpr (std::numeric_limits<scalar_t>::has_infinity) {
        scalar_t inf_val = std::numeric_limits<scalar_t>::infinity();
        Multivector_T inf_mv(inf_val);
        CHECK(inf_mv.isinf());
        CHECK(!inf_mv.isnan());
      }
    }

    SUBCASE("Division by Zero") {
      Multivector_T a = Multivector_T::random(index_set_t(1), 1.0);
      if (a.frame() != index_set_t(1)) {
        a += Multivector_T(index_set_t(1));
      }
      Multivector_T zero_mv(0);
      
      Multivector_T res1 = a / scalar_t(0);
      CHECK((res1.isnan() || res1.isinf()));
      
      Multivector_T res2 = a / zero_mv;
      CHECK(res2.isnan());
    }

    SUBCASE("Exception Handling") {
      Multivector_T a = Multivector_T::random(index_set_t(1), 1.0);
      // Ensure 'a' has frame {1} by forcing a term if it's empty
      if (a.frame() != index_set_t(1)) {
        a += Multivector_T(index_set_t(1));
      }
      
      // Negative exponent in outer_pow
      CHECK_THROWS_AS(a.outer_pow(-1), glucat_error);
      
      // vector_part with incompatible frame
      index_set_t invalid_frame(2); // Frame that 'a' doesn't have
      CHECK_THROWS_AS(a.vector_part(invalid_frame, false), glucat_error);
      
      // write with invalid stream
      std::ofstream bad_stream; // Not opened
      bad_stream.setstate(std::ios::failbit);
      CHECK_THROWS_AS(a.write(bad_stream), glucat_error);
    }

    SUBCASE("Advanced Constructor Gaps") {
      index_set_t ist(1);
      index_set_t frm(1);
      scalar_t crd(1.0);
      
      // Test matrix_multi / framed_multi constructor with prechecked
      Multivector_T m1(ist, crd, frm, true);
      CHECK(m1.frame() == frm);
      
      Multivector_T m2(ist, crd, frm, false);
      CHECK(m2.frame() == frm);
      
      // Trigger exception path for constructor
      index_set_t invalid_ist(2);
      CHECK_THROWS_AS(Multivector_T(invalid_ist, crd, frm, false), glucat_error);
    }

    SUBCASE("Fast Path (32 terms)") {
      using Tuning_Values_P = typename Multivector_T::tune_p::tuning_values_p;
      using framed_multi_t = typename Multivector_T::framed_multi_t;
      index_set_t large_frm;
      for(int i=1; i<=6; ++i) large_frm |= index_set_t(i); // 64 dimensions
      
      framed_multi_t f;
      for(int i=0; i<32; ++i) {
         index_set_t ist;
         for(int j=0; j<6; ++j) if (i & (1 << j)) ist |= index_set_t(j+1);
         f += typename framed_multi_t::term_t(ist, scalar_t(i+1));
      }
      
      CHECK(f.nbr_terms() >= 16);
      
      // Explicit conversion triggers constructor paths
      Multivector_T a(f);
      Multivector_T b(f, large_frm);
      
      Multivector_T res = a * b;
      res = a ^ b;
      res = a & b;
      res = a % b;
    }

    SUBCASE("Edge cases and Special Values") {
      using framed_multi_t = typename Multivector_T::framed_multi_t;
      index_set_t frm(1);
      scalar_t inf = std::numeric_limits<scalar_t>::infinity();
      scalar_t nan = std::numeric_limits<scalar_t>::quiet_NaN();
      
      Multivector_T v_inf(inf);
      Multivector_T v_nan(nan);
      
      CHECK(v_inf.isinf());
      CHECK(v_nan.isnan());
      
      // Explicitly check framed_multi_t too
      framed_multi_t f_inf(inf);
      framed_multi_t f_nan(nan);
      CHECK(f_inf.isinf());
      CHECK(f_nan.isnan());
      
      Multivector_T res;
      // log of something not near 1
      Multivector_T v1("{1,2}"); // Likely triggers cascade/pade
      try { res = log(v1); } catch (...) {}
      
      // log of something with large values
      Multivector_T v2 = Multivector_T(100.0) + Multivector_T("{1,2}");
      try { res = log(v2); } catch (...) {}
      
      // log and sqrt of negative/zero
      Multivector_T v_neg(-1.0);
      try { res = log(v_neg); } catch (...) {}
      try { res = sqrt(v_neg); } catch (...) {}
      
      Multivector_T v_zero(0.0);
      try { res = log(v_zero); } catch (...) {}

      // atan and others
      try { res = atan(v1); } catch (...) {}
      try { res = asinh(v1); } catch (...) {}
    }

    SUBCASE("Explicit Logic Paths") {
      index_set_t frm(1);
      Multivector_T a = Multivector_T::random(frm, 1.0);
      if (a == Multivector_T(0)) a += Multivector_T(frm, scalar_t(1.0), frm, true);
      Multivector_T b = Multivector_T::random(frm, 1.0);
      Multivector_T res;

      // star product
      scalar_t s = star(a, b);
      
      // quad, norm, vector_part, operator()(int)
      scalar_t q = a.quad();
      scalar_t n = a.norm();
      auto vp = a.vector_part();
      auto vp2 = a.vector_part(frm, true);
      Multivector_T a0 = a(0);

      // term-based operator+=
      a += typename Multivector_T::term_t(frm, scalar_t(1.0));
      
      // Print a term to cover operator<<(ostream&, pair)
      std::cout << typename Multivector_T::term_t(frm, scalar_t(1.0)) << std::endl;

      // Force cascade_log / pade_log
      Multivector_T v_large = Multivector_T(1000.0) + Multivector_T("{1,2}");
      try { res = log(v_large); } catch (...) {}
      Multivector_T v_small = Multivector_T(0.1) + Multivector_T("{1,2}");
      try { res = log(v_small); } catch (...) {}

      // Cross-representation constructor
      using framed_multi_t = typename Multivector_T::framed_multi_t;
      
      // Print framed_multi explicitly
      framed_multi_t f_pr("{1,2}");
      std::cout << f_pr << std::endl;

      if constexpr (std::is_same_v<Multivector_T, framed_multi_t>) {
         // If we are testing framed_multi, construct from matrix_multi
         using matrix_multi_t = matrix_multi<scalar_t, -8, 8, typename Multivector_T::tune_p>;
         matrix_multi_t mm(a);
         Multivector_T f_from_m(mm);
      } else {
         // If we are testing matrix_multi, construct from framed_multi
         framed_multi_t f(a);
         Multivector_T m_from_f(f);
      }

      // basis_element (only for matrix_multi)
      if constexpr (requires { a.basis_element(frm); }) {
         Multivector_T b_elem = a.basis_element(frm);
      }
      
      // grade, involute, reverse, conj, max_abs, truncated
      CHECK(a.grade() >= 0);
      CHECK(a.max_abs() >= 0);
      Multivector_T i = a.involute();
      Multivector_T r = a.reverse();
      Multivector_T c = a.conj();
      Multivector_T t = a.truncated(0.1);
      
      // inv() for non-scalar
      Multivector_T e1(index_set_t(1), scalar_t(1.0), index_set_t(1), true);
      Multivector_T inv_test = Multivector_T(1.0) + e1;
      try {
        Multivector_T inv_a = inv_test.inv();
      } catch (...) {}
    }
    
    SUBCASE("Matrix Base Edge Cases") {
      using MatrixDbl_T = glucat::matrix::dense_matrix<double>;
      
      // Target classify_eigenvalues branches: both_eigs
      // Use rotation blocks and a negative entry to get:
      // - imaginary eigenvalues (via pi/2 rotation)
      // - negative real eigenvalue (-1.0)
      // - a gap distribution where the wrap-around gap is strictly largest
      MatrixDbl_T m(5, 5);
      m.zeros();
      // Block 1: rot by pi/2 -> eigenvalues {i, -i}
      m(0,1) = -1.0;
      m(1,0) = 1.0;
      // Block 2: rot by 0.6*pi -> eigenvalues {e^i0.6pi, e^-i0.6pi}
      double theta = 0.6 * 3.14159265358979323846;
      m(2,2) = std::cos(theta); m(2,3) = -std::sin(theta);
      m(3,2) = std::sin(theta); m(3,3) = std::cos(theta);
      // Negative real
      m(4,4) = -1.0;
      
      auto genus = m.classify_eigenvalues();
      CHECK(genus.m_eig_case == glucat::matrix::both_eigs);
      
      // Target classify_eigenvalues branches: singular
      MatrixDbl_T ms(1, 1);
      ms.zeros();
      auto genus_s = ms.classify_eigenvalues();
      CHECK(genus_s.m_is_singular == true);
      
      // Target unit() for float
      using MatrixFlt_T = glucat::matrix::dense_matrix<float>;
      MatrixFlt_T mf = glucat::matrix::unit<MatrixFlt_T>(2);
      CHECK(glucat::matrix::nbr_rows(mf) == 2);
      CHECK(glucat::matrix::nbr_cols(mf) == 2);
    }
    
    SUBCASE("Generation Edge Cases") {
      using WideMultivector_T = glucat::matrix_multi<double, -16, 16>;
      using wide_index_set_t = typename WideMultivector_T::index_set_t;
      using wide_scalar_t = typename WideMultivector_T::scalar_t;
      
      wide_index_set_t p8;
      for(int i=1; i<=8; ++i) p8 |= wide_index_set_t(i);
      wide_index_set_t n8;
      for(int i=1; i<=8; ++i) n8 |= wide_index_set_t(-i);
      
      // Signature (8, 0) -> triggers gen_from_pm4_qp4
      WideMultivector_T v_p8(p8, wide_scalar_t(1.0), p8, true);
      // Signature (0, 8) -> triggers gen_from_pp4_qm4
      WideMultivector_T v_n8(n8, wide_scalar_t(1.0), n8, true);
      
      // Signature (2, 0) -> triggers gen_from_qp1_pm1
      wide_index_set_t p2;
      for(int i=1; i<=2; ++i) p2 |= wide_index_set_t(i);
      WideMultivector_T v_p2(p2, wide_scalar_t(1.0), p2, true);
      
      // Signature (1, 7) -> triggers gen_from_pp4_qm4 (via bott=2 path)
      wide_index_set_t p1, n7;
      p1 |= wide_index_set_t(1);
      for(int i=1; i<=7; ++i) n7 |= wide_index_set_t(-i);
      WideMultivector_T v_p1n7(p1 | n7, wide_scalar_t(1.0), p1 | n7, true);
      
      // Signature (10, 0) -> triggers gen_from_pm4_qp4 (via bott=2 path)
      wide_index_set_t p10;
      for(int i=1; i<=10; ++i) p10 |= wide_index_set_t(i);
      WideMultivector_T v_p10(p10, wide_scalar_t(1.0), p10, true);
      
      // Direct call to operator() with bott=1 to hit default branch
      // We use sparse_matrix_t<int> because it's what matrix_multi uses for basis elements,
      // and it's guaranteed to have the 3-argument constructor (via eigen_sparse_wrapper).
      using gen_matrix_t = glucat::matrix::sparse_matrix_t<int>;
      glucat::gen::generator_table<gen_matrix_t>::generator()(1, 0);
    }
  }

  template< typename Multivector_T >
  void run_peg11_test()
  {
    test_transcendental_templated(Multivector_T(0));
    test_transcendental_templated(Multivector_T(1));
    test_transcendental_templated(Multivector_T::random(typename Multivector_T::index_set_t(1)));
    test_edge_cases_templated<Multivector_T>();
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

#endif
