#ifndef GLUCAT_TEST_PEG11_H
#define GLUCAT_TEST_PEG11_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg11.cpp : programming example 11 : Square root and transcendental functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2009 by Paul C. Leopardi
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

namespace peg11
{
  using namespace std;
  using namespace glucat;

  template< class Multivector_T >
  static
  void 
  check(const Multivector_T& A, const Multivector_T& B, const std::string& msg, const bool need_inv = false)
  {
    typedef Multivector_T m_;
    typedef typename m_::scalar_t s_;
    static const s_ tol = numeric_traits<s_>::sqrt(numeric_limits<s_>::epsilon());
    const s_ tol2 = tol*tol;
    const bool relative = norm(A) > tol2 && norm(B) > tol2;

    const s_ abs_norm_diff = norm(A-B);
    const s_ norm_diff = (relative) ? abs_norm_diff/norm(A) : abs_norm_diff;
    const bool A_isnan = A.isnan();
    const bool B_isnan = B.isnan();
    if ((A_isnan && !B_isnan) 
     ||(!A_isnan && !need_inv && B_isnan) 
     ||(!A_isnan && need_inv && !inv(A).isnan() && B_isnan) 
     || norm_diff > tol2)
      std::cout << "Test failed: " << msg << std::endl;
    if (norm_diff > tol2)
    {
      const streamsize prec = cout.precision(10);
      std::cout <<   "Absolute norm of difference == " << sqrt(abs_norm_diff) << std::endl;
      if (relative)
        std::cout << "Relative norm of difference == " << sqrt(norm_diff) << std::endl;
      cout.precision(prec);
    }
  }

  template< class Multivector_T >
  static
  void 
  transcendtest(const Multivector_T& A, const bool random=false)
  {
    typedef Multivector_T m_;
    typedef typename m_::scalar_t s_;
    const streamsize prec = cout.precision(numeric_limits<s_>::digits10);
    if (random)
      cout << "Random A in " << A.frame() << endl;
    else
      A.write("A");

    check(m_(1), exp(A)*exp(-A),    "exp(A)*exp(-A) != 1");
    check(m_(1), exp(-A)*exp(A),    "exp(-A)*exp(A) != 1");
    check(exp(A), exp(scalar(A))*exp(pure(A)),
                                    "exp(scalar(A))*exp(pure(A)) != exp(A)");
    check(exp(elliptic(A)*A), cos(A)+elliptic(A)*sin(A), 
                                    "cos(A)+elliptic(A)*sin(A) != exp(elliptic(A)*A)");
    check(exp(A), cosh(A)+sinh(A),  "cosh(A)+sinh(A) != exp(A)");
    check(sin(A), cos(A)*tan(A),    "cos(A)*tan(A) != sin(A)");
    check(sinh(A), cosh(A)*tanh(A), "cosh(A)*tanh(A) != sinh(A)");

    check(A, sqrt(A)*sqrt(A), "sqrt(A)*sqrt(A) != A");
    check(A, exp(log(A)),     "exp(log(A)) != A", true);
    check(A, cos(acos(A)),    "cos(acos(A)) != A", true);
    check(A, cosh(acosh(A)),  "cosh(acosh(A)) != A", true);
    check(A, sin(asin(A)),    "sin(asin(A)) != A", true);
    check(A, sinh(asinh(A)),  "sinh(asinh(A)) != A", true);
    check(A, tan(atan(A)),    "tan(atan(A)) != A", true);
    check(A, tanh(atanh(A)),  "tanh(atanh(A)) != A", true);
    cout << endl;
    cout.precision(prec);
  }

  template< class Multivector_T >
  static
  void
  rand_transcendtest(int n)
  {
    typedef Multivector_T m_;
    typedef typename m_::scalar_t s_;
    typedef typename m_::index_set_t e_;
    static const index_t v_lo = e_::v_lo;
    static const index_t v_hi = e_::v_hi;
    static const index_t max_n = std::min(-v_lo, v_hi);
    const index_t max_index = index_t(n);
    if (max_index > max_n)
    {
      cout << "Value " << max_index << " is too big." << endl;
      cout << "Maximum value possible is " << max_n << "." << endl;
      return;
    }
    m_ a = s_(1);
    const s_ RAND_SCALE = 1.0/RAND_MAX;
    for (index_t i = 1; i != max_index+1; ++i)
      for (int k = 0; k != 2; ++k)
      {
        a *= (m_(e_((2*k-1)*i),  s_(rand()*RAND_SCALE)) + s_(rand()*RAND_SCALE*2.0));
        transcendtest(a, true);
      }
  }

  template< class Multivector_T >
  static
  void 
  do_test11()
  {
    typedef Multivector_T m_;
    typedef typename m_::scalar_t s_;
    const s_ pi = numeric_traits<s_>::pi();
    const m_ i("{-1}");
    transcendtest(m_(0));
    transcendtest(m_(1));
    transcendtest(-m_(1));
    transcendtest(m_(2));
    transcendtest(-m_(2));
    transcendtest(m_(s_(0.5)));
    transcendtest(-m_(s_(0.5)));
    transcendtest(m_(pi));
    // exp(i*pi) should be -1
    transcendtest(pi*m_("{1,2,3}"));
    transcendtest(-pi*m_("{1,2,3}"));
    transcendtest(pi*m_("{-1,1,2,3,4}"));
    transcendtest(-pi*m_("{-1,1,2,3,4}"));
    transcendtest(pi*m_("{-2,-1,1,2,3,4,5}"));
    transcendtest(-pi*m_("{-2,-1,1,2,3,4,5}"));
    transcendtest(pi*m_("{-2,-1,1,2,3,4,5}")*s_(10.0));
    transcendtest(-pi*m_("{-2,-1,1,2,3,4,5}")*s_(10.0));
    transcendtest(pi*m_("{1,2}"));
    transcendtest(pi*i/s_(2.0));
    transcendtest(-pi*i/s_(2.0));
    transcendtest(pi*i*s_(100.0));
    transcendtest(-pi*i*s_(100.0));
    transcendtest(m_("{1}"));
    transcendtest(-m_("{1}"));
    transcendtest(m_("{-1,1}"));
    transcendtest(-m_("{-1,1}"));
    transcendtest(m_("{-2,-1,1,2}"));
    transcendtest(-m_("{-2,-1,1,2}"));
    transcendtest(m_("{-1}+{1}"));
    // the following should produce NaN values
    transcendtest(m_("{-1}+{1}")*s_(numeric_limits<s_>::max()));
    rand_transcendtest<m_>(4);
  }
}

int test11();

#endif
