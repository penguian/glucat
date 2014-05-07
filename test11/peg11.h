#ifndef GLUCAT_TEST_PEG11_H
#define GLUCAT_TEST_PEG11_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg11.cpp : programming example 11 : Square root and transcendental functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2012 by Paul C. Leopardi
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

#include "glucat/glucat.h"

typedef glucat::tuning
  <
    glucat::DEFAULT_Mult_Matrix_Threshold,
    glucat::DEFAULT_Div_Max_Steps,
    glucat::DEFAULT_Sqrt_Max_Steps,
    glucat::DEFAULT_Log_Max_Outer_Steps,
    glucat::DEFAULT_Log_Max_Inner_Steps,
    glucat::DEFAULT_Basis_Max_Count,
    glucat::DEFAULT_Fast_Size_Threshold,
    glucat::DEFAULT_Inv_Fast_Dim_Threshold,
    glucat::DEFAULT_Products_Size_Threshold,
    glucat::precision_same,
    glucat::precision_promoted
  >
  Tune_P;

#include "glucat/glucat_imp.h"
#include <iomanip>
#include <cstdio>
#include "test/try_catch.h"
#include "test/control.h"

namespace peg11
{
  using namespace std;
  using namespace glucat;

  template< class Multivector_T >
  static
  void
  check(const Multivector_T& A, const Multivector_T& B, const string& msg, const bool need_inv = false)
  {
    typedef typename Multivector_T::scalar_t scalar_t;

    static const scalar_t scalar_eps  = numeric_limits<scalar_t>::epsilon();
    scalar_t tol2;
    if (test_control.m_verbose_output)
      tol2 = 0.0;
    else
    {
     typedef typename Multivector_T::framed_multi_t framed_multi_t;
      const double nbr_terms = framed_multi_t(A).truncated(scalar_eps).nbr_terms();
      scalar_t tol = scalar_eps *
                     numeric_traits<scalar_t>::pow(scalar_t(2), numeric_limits<scalar_t>::digits / 16 + 4);
      tol2 = tol * tol * scalar_t(std::max(nbr_terms, 1.0));
    }
    const bool relative = (norm(A) > tol2) && (norm(B) > tol2);
    const scalar_t abs_norm_diff = norm(A-B);
    const scalar_t norm_diff = (relative) ? abs_norm_diff/norm(A) : abs_norm_diff;
    const bool A_isnan = A.isnan();
    const bool B_isnan = B.isnan();
    if ((A_isnan && !B_isnan)
     ||(!A_isnan && !need_inv && B_isnan)
     ||(!A_isnan && need_inv && !inv(A).isnan() && B_isnan)
     || norm_diff > tol2)
    {
      cout << "Test failed: " << msg << endl;
      if (norm_diff > tol2)
      {
        const streamsize prec = cout.precision(5);
        cout << ((relative) ? "Relative" : "Absolute");
        cout << " norm of difference == "
             << numeric_traits<scalar_t>::sqrt(norm_diff) << endl;
         if (!test_control.m_verbose_output)
         {
           cout.precision(numeric_limits<scalar_t>::digits10);
           cout << "lhs==" << B << endl;
           cout << "rhs==" << A << endl;
         }
        cout.precision(prec);
      }
      else
        cout << "lhs == " << B << endl;
    }
  }

  template< class Multivector_T >
  static
  void
  transcendtest(const Multivector_T& A, const bool random=false)
  {
    typedef Multivector_T m_;
    typedef typename m_::scalar_t scalar_t;
    const streamsize prec = cout.precision(numeric_limits<scalar_t>::digits10);
    if (random)
      cout << "Random A in " << A.frame() << endl;
    else
      A.write("A");

    check(m_(1), exp(A)*exp(-A),    "exp(A)*exp(-A) != 1");
    check(m_(1), exp(-A)*exp(A),    "exp(-A)*exp(A) != 1");
    check(exp(A), exp(scalar(A))*exp(pure(A)),
                                    "exp(scalar(A))*exp(pure(A)) != exp(A)");
    check(exp(complexifier(A)*A), cos(A)+complexifier(A)*sin(A),
                                    "cos(A)+complexifier(A)*sin(A) != exp(complexifier(A)*A)");
    check(exp(A), cosh(A)+sinh(A),  "cosh(A)+sinh(A) != exp(A)");
    check(sin(A), cos(A)*tan(A),    "cos(A)*tan(A) != sin(A)");
    check(sinh(A), cosh(A)*tanh(A), "cosh(A)*tanh(A) != sinh(A)");

    if ((A == scalar(A)) || !((inv(A)).isnan()))
      check(A, sqrt(A)*sqrt(A), "sqrt(A)*sqrt(A) != A");
    if (!((inv(A)).isnan()))
      check(A, exp(log(A)),   "exp(log(A)) != A", true);
    check(A, cos(acos(A)),    "cos(acos(A)) != A", true);
    check(A, cosh(acosh(A)),  "cosh(acosh(A)) != A", true);
    check(A, sin(asin(A)),    "sin(asin(A)) != A", true);
    check(A, sinh(asinh(A)),  "sinh(asinh(A)) != A", true);
    check(A, tan(atan(A)),    "tan(atan(A)) != A", true);
    if (!(log(m_(1)+A).isnan() || log(m_(1)-A).isnan()))
      check(A, tanh(atanh(A)),"tanh(atanh(A)) != A", true);
    cout << endl;
    cout.precision(prec);
  }

  template< class Multivector_T >
  static
  void
  rand_transcendtest(int n)
  {
    typedef Multivector_T m_;
    typedef typename m_::index_set_t index_set_t;
    typedef typename index_set_t::index_pair_t index_pair_t;

    static const index_t v_lo = index_set_t::v_lo;
    static const index_t v_hi = index_set_t::v_hi;
    static const index_t max_n = min(-v_lo, v_hi);
    const index_t max_index = index_t(n);
    if (max_index > max_n)
    {
      cout << "Value " << max_index << " is too big." << endl;
      cout << "Maximum value possible is " << max_n << "." << endl;
      return;
    }
    index_set_t frm;

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
          transcendtest(m_::random(frm), true);
      }
    }
  }

  template< class Multivector_T >
  static
  void
  do_test11()
  {
    typedef Multivector_T m_;
    typedef typename m_::scalar_t scalar_t;

    cout << "Epsilon ==" << numeric_limits<scalar_t>::epsilon() << endl;

    const scalar_t pi = numeric_traits<scalar_t>::pi();
    const m_ i("{-1}");
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
    transcendtest(m_("{-2,-1,1,2}"));
    transcendtest(-m_("{-2,-1,1,2}"));
    transcendtest(m_("{-1}+{1}"));

    rand_transcendtest<m_>(3);
  }
}

int test11();

#endif
