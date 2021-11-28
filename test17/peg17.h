#ifndef GLUCAT_TEST_PEG17_H
#define GLUCAT_TEST_PEG17_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg17.cpp : programming example 17 : Truncation and printing
                             -------------------
    begin                : Sun 2021-11-24
    copyright            : (C) 2001-2025 by Paul C. Leopardi
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

#include <string>

using namespace std;
using namespace glucat;

template <typename T>
void print(const string& str, const T& val)
{
  cout << str << ": " << endl
       << val << endl;
}

template <typename T>
void print_truncated(
  const string& str, const T& val,
  const typename T::scalar_t truncation = T::default_truncation)
{
  cout << "truncation: "
       <<  truncation << endl;
  cout << str << ".truncated(truncation): " << endl
       << val.truncated(truncation) << endl;
}

template <typename Scalar_T>
void printing_test(const streamsize precision=9)
{
  using cm = matrix_multi<Scalar_T>;
  using cf = framed_multi<Scalar_T>;

  const auto floatfield = ios_base::floatfield;
  const auto fixed = ios_base::fixed;
  const auto scientific = ios_base::scientific;
  const auto hexfloat = fixed | scientific;

  const auto a = cf("{-3}+{-2}+{-1}");
  const auto b = cf("{-2}+1.e-4{1}+{2}+1.e4{3}");

  const auto frame = a.frame() | b.frame();
  const auto A = cm(a, frame);
  const auto B = cm(b, frame);

  cout << precision << " digit precision" << endl;
  const streamsize old_precision = cout.precision(precision);

  print("a", a);
  print("A = cm(a)", A);
  print("b", b);
  print("B = cm(b)", B);

  cout << "Fixed format" << endl;
  cout.setf(fixed, floatfield);
  print("a", a);
  print("A", A);
  print("b", b);
  print("B", B);

  cout << "Scientific format" << endl;
  cout.setf(scientific, floatfield);
  print("a", a);
  print("A", A);
  print("b", b);
  print("B", B);

  cout << "Hexadecimal floating point format" << endl;
  cout.setf(hexfloat, floatfield);
  print("a", a);
  print("A", A);
  print("b", b);
  print("B", B);

  cout << "Default floating point format" << endl;
  cout.unsetf(floatfield);
  print("a", a);
  print("A", A);
  print("b", b);
  print("B", B);

  cout << "Default truncation" << endl;

  cout << "Fixed format" << endl;
  cout.setf(fixed, floatfield);
  print_truncated("a", a);
  print_truncated("A", A);
  print_truncated("b", b);
  print_truncated("B", B);

  cout << "Scientific format" << endl;
  cout.setf(scientific, floatfield);
  print_truncated("a", a);
  print_truncated("A", A);
  print_truncated("b", b);
  print_truncated("B", B);

  cout << "Hexadecimal floating point format" << endl;
  cout.setf(hexfloat, floatfield);
  print_truncated("a", a);
  print_truncated("A", A);
  print_truncated("b", b);
  print_truncated("B", B);

  cout << "Default floating point format" << endl;
  cout.unsetf(floatfield);
  print_truncated("a", a);
  print_truncated("A", A);
  print_truncated("b", b);
  print_truncated("B", B);

  cout.precision(old_precision);
}

int test17();

#endif
