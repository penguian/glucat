#ifndef GLUCAT_TEST_PEG11_H
#define GLUCAT_TEST_PEG11_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg11.cpp : programming example 11 : Square root and transcendental functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2007 by Paul C. Leopardi
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
    const s_ tol = s_(numeric_limits<s_>::epsilon()) * s_(1 << 7);
    const s_ tol2 = tol*tol;
    const s_ norm_diff = (norm(A) < tol2 || norm(B) < tol2) ? norm(A-B) : norm(A-B)/norm(A);
    const bool A_isnan = A.isnan();
    const bool B_isnan = B.isnan();
    if ((A_isnan && !B_isnan) 
     ||(!A_isnan && !need_inv && B_isnan) 
     ||(!A_isnan && need_inv && !inv(A).isnan() && B_isnan) 
     || norm_diff > tol2)
    {
      std::cout << "Test failed: ";
      B.write(msg);
      if (norm_diff > tol2)
      {
        if (norm(A) < tol2 || norm(B) < tol2)
          std::cout << "Norm";
        else
          std::cout << "Relative norm";
        std::cout << " of difference == " << sqrt(norm_diff) << std::endl;
      }
    } 
  }

  template< class Multivector_T >
  static
  void 
  transcendtest(const Multivector_T& A)
  {
    typedef Multivector_T m_;
    typedef typename m_::scalar_t s_;
    const streamsize prec = cout.precision(numeric_limits<s_>::digits10);
    A.write("A");
    check(m_(1), exp(A)*exp(-A),    "exp(A)*exp(-A) != 1");
    check(m_(1), exp(-A)*exp(A),    "exp(-A)*exp(A) != 1");
    check(exp(A), std::exp(scalar(A))*exp(pure(A)),
                                    "std::exp(scalar(A))*exp(pure(A)) != exp(A)");
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
  do_test11()
  {
    typedef Multivector_T m_;
    typedef typename m_::scalar_t s_;
    const s_ pi = s_(l_pi);
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
  }
}

int test11();

#endif
