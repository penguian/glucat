#ifndef GLUCAT_TEST_PEG10_H
#define GLUCAT_TEST_PEG10_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg10.cpp : programming example 10 : Norms and involutions
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

#include <iomanip>
#include <time.h>

namespace peg10
{
  using namespace glucat;
  using namespace std;

  template< class Multivector_T >
  static
  void 
  test_norms_and_involutions(const Multivector_T& a, const Multivector_T& b)
  {
    typedef Multivector_T number;
    typedef typename number::scalar_t Scalar_T;

    const int cwidth = 24;
    const int width = 14;
    const int prec = 12;

    number ab = a*b;
    number ba = b*a;

    Scalar_T norm_a = a.norm();
    Scalar_T norm_b = b.norm();
    Scalar_T norm_ab = ab.norm();

    cout << setw(cwidth) << "a= "
    << setw(width) << setprecision(prec)
    << a.norm() << " "
    << a << endl;
    cout << setw(cwidth) << "involute(a)= "
    << setw(width) << setprecision(prec)
    << involute(a).norm() << " "
    << involute(a) << endl;
    cout << setw(cwidth) << "reverse(a)=  "
    << setw(width) << setprecision(prec)
    << reverse(a).norm() << " "
    << reverse(a) << endl;
    cout << setw(cwidth) << "conj(a)= "
    << setw(width) << setprecision(prec)
    << conj(a).norm() << " "
    << conj(a) << endl;
    cout << setw(cwidth) << "b= "
    << setw(width) << setprecision(prec)
    << b.norm() << " "
    << b << endl;
    cout << setw(cwidth) << "involute(b)= "
    << setw(width) << setprecision(prec)
    << involute(b).norm() << " "
    << involute(b) << endl;
    cout << setw(cwidth) << "reverse(b)=  "
    << setw(width) << setprecision(prec)
    << reverse(b).norm() << " "
    << reverse(b) << endl;
    cout << setw(cwidth) << "conj(b)= "
    << setw(width) << setprecision(prec)
    << conj(b).norm() << " "
    << conj(b) << endl;

    cout << setw(cwidth) << "norm(a)norm(b)= "
    << setw(width) << setprecision(prec)
    << norm_a*norm_b << endl;
    cout << setw(cwidth) << "(pn(ab)+pn(reverse(a)b))/2= "
    << setw(width) << setprecision(prec)
    << (norm_ab+(reverse(a)*b).norm())/2.0 << endl;

    cout << setw(cwidth) << "ab= "
    << setw(width) << setprecision(prec)
    << ab.norm() << " "
    << ab << endl;
    cout << setw(cwidth) << "involute(a)involute(b)= "
    << setw(width) << setprecision(prec)
    << (involute(a)*involute(b)).norm() << " "
    << involute(a)*involute(b) << endl;
    cout << setw(cwidth) << "reverse(b)reverse(a)= "
    << setw(width) << setprecision(prec)
    << (reverse(b)*reverse(a)).norm() << " "
    << reverse(b)*reverse(a) << endl;
    cout << setw(cwidth) << "conj(b)conj(a)= "
    << setw(width) << setprecision(prec)
    << (conj(b)*conj(a)).norm() << " "
    << conj(b)*conj(a) << endl;

    cout << setw(cwidth) << "ba= "
    << setw(width) << setprecision(prec)
    << ba.norm() << " "
    << ba << endl;
    cout << setw(cwidth) << "involute(b)involute(a)= "
    << setw(width) << setprecision(prec)
    << (involute(b)*involute(a)).norm() << " "
    << involute(b)*involute(a) << endl;
    cout << setw(cwidth) << "reverse(a)reverse(b)= "
    << setw(width) << setprecision(prec)
    << (reverse(a)*reverse(b)).norm() << " "
    << reverse(a)*reverse(b) << endl;
    cout << setw(cwidth) << "conj(a)conj(b)= "
    << setw(width) << setprecision(prec)
    << (conj(a)*conj(b)).norm() << " "
    << conj(a)*conj(b) << endl;

    cout << setw(cwidth) << "conj(a)b= "
    << setw(width) << setprecision(prec)
    << (conj(a)*b).norm() << " "
    << conj(a)*b << endl;
    cout << setw(cwidth) << "reverse(a)involute(b)= "
    << setw(width) << setprecision(prec)
    << (reverse(a)*involute(b)).norm() << " "
    << reverse(a)*involute(b) << endl;
    cout << setw(cwidth) << "conj(b)a= "
    << setw(width) << setprecision(prec)
    << (conj(b)*a).norm() << " "
    << conj(b)*a << endl;
    cout << setw(cwidth) << "reverse(b)involute(a)= "
    << setw(width) << setprecision(prec)
    << (reverse(b)*involute(a)).norm() << " "
    << reverse(b)*involute(a) << endl;

    cout << setw(cwidth) << "involute(a)reverse(b)= "
    << setw(width) << setprecision(prec)
    << (involute(a)*reverse(b)).norm() << " "
    << involute(a)*reverse(b) << endl;
    cout << setw(cwidth) << "a*conj(b)= "
    << setw(width) << setprecision(prec)
    << (a*conj(b)).norm() << " "
    << a*conj(b) << endl;
    cout << setw(cwidth) << "involute(b)reverse(a)= "
    << setw(width) << setprecision(prec)
    << (involute(b)*reverse(a)).norm() << " "
    << involute(b)*reverse(a) << endl;
    cout << setw(cwidth) << "b*conj(a)= "
    << setw(width) << setprecision(prec)
    << (b*conj(a)).norm() << " "
    << b*conj(a) << endl;

    cout << setw(cwidth) << "involute(a)b= "
    << setw(width) << setprecision(prec)
    << (involute(a)*b).norm() << " "
    << involute(a)*b << endl;
    cout << setw(cwidth) << "a*involute(b)= "
    << setw(width) << setprecision(prec)
    << (a*involute(b)).norm() << " "
    << a*involute(b) << endl;
    cout << setw(cwidth) << "conj(b)reverse(a)= "
    << setw(width) << setprecision(prec)
    << (conj(b)*reverse(a)).norm() << " "
    << conj(b)*reverse(a) << endl;
    cout << setw(cwidth) << "reverse(b)conj(a)= "
    << setw(width) << setprecision(prec)
    << (reverse(b)*conj(a)).norm() << " "
    << reverse(b)*conj(a) << endl;

    cout << setw(cwidth) << "involute(b)a= "
    << setw(width) << setprecision(prec)
    << (involute(b)*a).norm() << " "
    << involute(b)*a << endl;
    cout << setw(cwidth) << "b*involute(a)= "
    << setw(width) << setprecision(prec)
    << (b*involute(a)).norm() << " "
    << b*involute(a) << endl;
    cout << setw(cwidth) << "conj(a)reverse(b)= "
    << setw(width) << setprecision(prec)
    << (conj(a)*reverse(b)).norm() << " "
    << conj(a)*reverse(b) << endl;
    cout << setw(cwidth) << "reverse(a)conj(b)= "
    << setw(width) << setprecision(prec)
    << (reverse(a)*conj(b)).norm() << " "
    << reverse(a)*conj(b) << endl;

    cout << setw(cwidth) << "reverse(a)b=  "
    << setw(width) << setprecision(prec)
    << (reverse(a)*b).norm() << " "
    << reverse(a)*b << endl;
    cout << setw(cwidth) << "conj(a)involute(b)= "
    << setw(width) << setprecision(prec)
    << (conj(a)*involute(b)).norm() << " "
    << conj(a)*involute(b) << endl;
    cout << setw(cwidth) << "reverse(b)a=  "
    << setw(width) << setprecision(prec)
    << (reverse(b)*a).norm() << " "
    << reverse(b)*a << endl;
    cout << setw(cwidth) << "conj(b)involute(a)= "
    << setw(width) << setprecision(prec)
    << (conj(b)*involute(a)).norm() << " "
    << conj(b)*involute(a) << endl;

    cout << setw(cwidth) << "a*reverse(b)= "
    << setw(width) << setprecision(prec)
    << (a*reverse(b)).norm() << " "
    << a*reverse(b) << endl;
    cout << setw(cwidth) << "involute(a)conj(b)= "
    << setw(width) << setprecision(prec)
    << (involute(a)*conj(b)).norm() << " "
    << involute(a)*conj(b) << endl;
    cout << setw(cwidth) << "b*reverse(a)= "
    << setw(width) << setprecision(prec)
    << (b*reverse(a)).norm() << " "
    << b*reverse(a) << endl;
    cout << setw(cwidth) << "involute(b)conj(a)= "
    << setw(width) << setprecision(prec)
    << (involute(b)*conj(a)).norm() << " "
    << involute(b)*conj(a) << endl;
  }
}

int test10();

#endif
