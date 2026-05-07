#ifndef _GLUCAT_TEST_DOCTEST_TEST_MATRIX_METHODS_H
#define _GLUCAT_TEST_DOCTEST_TEST_MATRIX_METHODS_H

#include "test_doctest.h"

namespace glucat {

  template< typename Multivector_T >
  void test_matrix_methods_templated(const Multivector_T& mm)
  {
    typedef typename Multivector_T::scalar_t scalar_t;
    typedef typename Multivector_T::framed_multi_t framed_multi_t;
    static const scalar_t eps = std::numeric_limits<scalar_t>::epsilon();
    const scalar_t tol = eps * 512.0;

    framed_multi_t fm(mm);

    SUBCASE("involute()") {
      CHECK_FALSE(is_error(Multivector_T(fm.involute()), mm.involute(), tol));
    }

    SUBCASE("even() and odd()") {
      CHECK_FALSE(is_error(Multivector_T(fm.even()), mm.even(), tol));
      CHECK_FALSE(is_error(Multivector_T(fm.odd()), mm.odd(), tol));
      CHECK_FALSE(is_error(mm.even() + mm.odd(), mm, tol));
    }

    SUBCASE("reverse() and conj()") {
      CHECK_FALSE(is_error(Multivector_T(fm.reverse()), mm.reverse(), tol));
      CHECK_FALSE(is_error(Multivector_T(fm.conj()), mm.conj(), tol));
    }
    
    SUBCASE("decompose() and project()") {
      auto grades = mm.decompose();
      Multivector_T sum(scalar_t(0), mm.frame());
      for (size_t i = 0; i < grades.size(); ++i) {
        sum += grades[i];
        // Verify individual projection matches decomposition
        CHECK_FALSE(is_error(mm.project(i), grades[i], tol));
      }
      // Verify sum matches original
      CHECK_FALSE(is_error(sum, mm, tol));

      // Verify against reference framed_multi projection
      for (size_t i = 0; i < grades.size(); ++i) {
        Multivector_T ref_i(fm(i), mm.frame());
        CHECK_FALSE(is_error(grades[i], ref_i, tol));
      }
    }
  }

  template< typename Multivector_T >
  void run_checkerboard_test()
  {
    using index_set_t = typename Multivector_T::index_set_t;
    
    // Representative signatures for the 8 types (p-q mod 8)
    struct Signature { int p; int q; };
    Signature signatures[] = {
      {0, 0}, // Type 0: R
      {1, 0}, // Type 1: C (or R+R?)
      {2, 0}, // Type 2: H
      {3, 0}, // Type 3: H+H
      {0, 4}, // Type 4: H
      {0, 3}, // Type 5: C
      {0, 2}, // Type 6: R
      {0, 1}  // Type 7: R+R
    };

    for (const auto& sig : signatures) {
      // Construct frame for Cl(p,q)
      index_set_t frame;
      for (int i = 1; i <= sig.p; ++i) frame |= index_set_t(i);
      for (int i = 1; i <= sig.q; ++i) frame |= index_set_t(-i);

      DOCTEST_SUBCASE((std::string("Signature Cl(") + std::to_string(sig.p) + "," + std::to_string(sig.q) + ")").c_str()) {
        Multivector_T mm = Multivector_T::random(frame);
        test_matrix_methods_templated(mm);
      }
    }
  }

} // namespace glucat

#endif // _GLUCAT_TEST_DOCTEST_TEST_MATRIX_METHODS_H
