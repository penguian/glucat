#ifndef GLUCAT_TEST_PEG11_H
#define GLUCAT_TEST_PEG11_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg11.cpp : programming example 11 : Square root and transcendental functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2021 by Paul C. Leopardi
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

#include <iomanip>
#include <cstdio>

namespace peg11
{
  using namespace std;
  using namespace glucat;

  template< class Multivector_T >
  static
  void
  check(const Multivector_T& lhs, const Multivector_T& rhs, const string& msg, const bool need_inv = false)
  {
    const auto lhs_isinf = lhs.isinf();
    const auto lhs_isnan = lhs.isnan();
    const auto rhs_isinf = rhs.isinf();
    const auto rhs_isnan = rhs.isnan();
    const auto exceptional =
      (rhs_isinf && !lhs_isinf) ||
      (rhs_isnan && !lhs_isnan) ||
      (!rhs_isinf && !need_inv && lhs_isinf) ||
      (!rhs_isnan && !need_inv && lhs_isnan) ||
      (!rhs_isinf && need_inv && !inv(rhs).isinf() && lhs_isinf) ||
      (!rhs_isnan && need_inv && !inv(rhs).isnan() && lhs_isnan);

    if (exceptional)
    {
      cout << "Test failed: " << msg << endl;
      cout << "lhs == " << lhs << endl;
      return;
    }
    else
    {
      const auto different = control_t::verbose()
        ? (lhs != rhs)
        : !approx_equal(lhs, rhs);
      if (different)
      {
        using scalar_t = typename Multivector_T::scalar_t;
        const auto norm_diff = norm_of_diff(rhs, lhs);
        const auto rhs_tol = norm_tol(rhs);
        const auto relative = (norm(rhs) > rhs_tol) && (norm(lhs) > rhs_tol);
        const auto prec = cout.precision(5);
        cout << "Test failed: " << msg << endl;
        cout << ((relative) ? "Relative" : "Absolute");
        cout << " norm of difference == "
             << numeric_traits<scalar_t>::sqrt(norm_diff) << endl;
        if (!control_t::verbose())
        {
          cout.precision(numeric_limits<scalar_t>::digits10);
          cout << "lhs == " << lhs << endl;
          cout << "rhs == " << rhs << endl;
        }
        cout.precision(prec);
      }
    }
  }

  template< class Multivector_T >
  static
  void
  transcendtest(const Multivector_T& A, const bool random=false)
  {
    using m_ = Multivector_T;
    using scalar_t = typename Multivector_T::scalar_t;
    const auto prec = cout.precision(numeric_limits<scalar_t>::digits10);
    if (random)
      cout << "Random A in " << A.frame() << endl;
    else
      A.write("A");

    check(exp(A)*exp(-A), m_(1),    "exp(A)*exp(-A) != 1");
    check(exp(-A)*exp(A), m_(1),    "exp(-A)*exp(A) != 1");
    check(exp(scalar(A))*exp(pure(A)), exp(A),
                                    "exp(scalar(A))*exp(pure(A)) != exp(A)");
    check(cos(A)+complexifier(A)*sin(A), exp(complexifier(A)*A),
                                    "cos(A)+complexifier(A)*sin(A) != exp(complexifier(A)*A)");
    check(cosh(A)+sinh(A), exp(A),  "cosh(A)+sinh(A) != exp(A)");
    check(cos(A)*tan(A), sin(A),    "cos(A)*tan(A) != sin(A)");
    check(cosh(A)*tanh(A), sinh(A), "cosh(A)*tanh(A) != sinh(A)");

    // if ((A == scalar(A)) || !((inv(A)).isnan()))
    check(sqrt(A)*sqrt(A), A,       "sqrt(A)*sqrt(A) != A");
    if (!((inv(A)).isnan()))
      check(exp(log(A)), A,         "exp(log(A)) != A", true);
    check(cos(acos(A)), A,          "cos(acos(A)) != A", true);
    check(cosh(acosh(A)), A,        "cosh(acosh(A)) != A", true);
    check(sin(asin(A)), A,          "sin(asin(A)) != A", true);
    check(sinh(asinh(A)), A,        "sinh(asinh(A)) != A", true);
    check(tan(atan(A)), A,          "tan(atan(A)) != A", true);
    if (!(log(m_(1)+A).isnan() || log(m_(1)-A).isnan()))
      check(tanh(atanh(A)), A,      "tanh(atanh(A)) != A", true);
    cout << endl;
    cout.precision(prec);
  }

  template< class Multivector_T >
  static
  void
  rand_transcendtest(int n)
  {
    using index_set_t = typename Multivector_T::index_set_t;
    using index_pair_t = typename index_set_t::index_pair_t;

    static const auto v_lo = index_set_t::v_lo;
    static const auto v_hi = index_set_t::v_hi;
    static const auto max_n = min(-v_lo, v_hi);
    const auto max_index = index_t(n);
    if (max_index > max_n)
    {
      cout << "Value " << max_index << " is too big." << endl;
      cout << "Maximum value possible is " << max_n << "." << endl;
      return;
    }
    auto frm = index_set_t();

    for (index_t p = 0; p != max_index+1; ++p)
    {
      frm = index_set_t(0);
      if (p != 0)
        frm = index_set_t(index_pair_t(1,p),true);

      for (index_t q = max(0,p-2); q != p+1; ++q)
      {
        if (q != 0)
        {
          if (q == p-2)
            frm |= index_set_t(index_pair_t(-q,-1),true);
          else
            frm |= index_set_t(-q);
        }
        for (int k=0; k!=2; ++k)
          transcendtest(Multivector_T::random(frm), true);
      }
    }
  }

  template< class Multivector_T >
  static
  void
  do_test11()
  {
    using m_ = Multivector_T;
    using scalar_t = typename Multivector_T::scalar_t;

    cout << "Epsilon ==" << numeric_limits<scalar_t>::epsilon() << endl;

    const auto pi = numeric_traits<scalar_t>::pi();
    const auto i = m_("{-1}");
    transcendtest(m_(0));
    transcendtest(m_(1));
    transcendtest(-m_(1));
    transcendtest(m_(2));
    transcendtest(-m_(2));
    transcendtest(m_(scalar_t(0.5)));
    transcendtest(-m_(scalar_t(0.5)));
    transcendtest(m_(pi));
    // exp(i*pi) should be -1
    transcendtest(pi*m_("{1,2,3}"));
    transcendtest(-pi*m_("{1,2,3}"));
    transcendtest(pi*m_("{-1,1,2,3,4}"));
    transcendtest(-pi*m_("{-1,1,2,3,4}"));
    transcendtest(pi*m_("{-2,-1,1,2,3,4,5}"));
    transcendtest(-pi*m_("{-2,-1,1,2,3,4,5}"));
    transcendtest(pi*m_("{-2,-1,1,2,3,4,5}")*scalar_t(10.0));
    transcendtest(-pi*m_("{-2,-1,1,2,3,4,5}")*scalar_t(10.0));
    transcendtest(pi*m_("{1,2}"));
    transcendtest(pi*i/scalar_t(2.0));
    transcendtest(-pi*i/scalar_t(2.0));
    transcendtest(pi*i*scalar_t(100.0));
    transcendtest(-pi*i*scalar_t(100.0));
    transcendtest(m_("{1}"));
    transcendtest(-m_("{1}"));
    transcendtest(m_("{-1,1}"));
    transcendtest(-m_("{-1,1}"));
    transcendtest(m_("1+{-1,1}"));
    transcendtest(m_("1-{-1,1}"));
    transcendtest(m_("{-1}+{1,2}"));
    transcendtest(m_("{-1}-{1,2}"));
    transcendtest(m_("{1}+{1,2}"));
    transcendtest(m_("{1}-{1,2}"));
    transcendtest(m_("{-2,-1,1,2}"));
    transcendtest(-m_("{-2,-1,1,2}"));
    transcendtest(m_("{-1}+{1}"));

    rand_transcendtest<m_>(3);
  }

  template <typename Scalar_T, const index_t LO = DEFAULT_LO, const index_t HI = DEFAULT_HI, typename Tune_P = tuning<> >
  void
  do_test11_tuned(const std::string tuning_name, const std::string scalar_typename)
  {
    std::cout << tuning_name << std::endl;
    std::cout << std::endl << "framed_multi<" << scalar_typename << "," << LO << "," << HI << ">" << std::endl;
    do_test11< framed_multi<Scalar_T,LO,HI,Tune_P> >();
    std::cout << std::endl << "matrix_multi<" << scalar_typename << "," << LO << "," << HI << ">" << std::endl;
    do_test11< matrix_multi<Scalar_T,LO,HI,Tune_P> >();
  }
}

int test11();

#endif
