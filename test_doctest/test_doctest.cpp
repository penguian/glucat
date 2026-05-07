#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

// Define this so the headers compile the TEST_CASE blocks
#define GLUCAT_DOCTEST

#include <glucat/glucat_imp.h>
#include "test_doctest.h"
#include "test_matrix_methods.h"
TEST_CASE("math::identities") {
  using namespace glucat;

  SUBCASE("matrix_multi<double>") {
    run_peg00_test<matrix_multi<double, -8, 8>>(8);
    run_peg11_test<matrix_multi<double, -8, 8>>();
    run_checkerboard_test<matrix_multi<double, -8, 8>>();
  }

  SUBCASE("framed_multi<double>") {
    run_peg00_test<framed_multi<double, -8, 8>>(2);
    run_peg11_test<framed_multi<double, -8, 8>>();
  }
}
