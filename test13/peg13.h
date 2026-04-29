#ifndef GLUCAT_TEST_PEG13_H
#define GLUCAT_TEST_PEG13_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg13.cpp : programming example 13 : Multiplication and division
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

template <typename Scalar_T>
void
do_test13(const std::string scalar_typename)
{
    using namespace glucat;
    using namespace std;

    std::cout << std::endl << scalar_typename << ":" << std::endl;

    using cm = matrix_multi<Scalar_T>;
    using cf = framed_multi<Scalar_T>;
    const cf a("{-3}+{-2}+{-1}");
    const cf b("{-2}+1.e-8{1}+{2}+1.e8{3}");

    const typename cf::index_set_t sub = a.frame() | b.frame();
    const cm A( a, sub );
    const cm B( b, sub );
    cout << "a: "           << endl << a << endl;
    cout << "A = cm(a): "   << endl << A << endl;
    cout << "b: "           << endl << b << endl;
    cout << "B = cm(b): "   << endl << B << endl;

    // Multiplication
    const cf c = a*b;
    cout << "a*b: "         << endl << c << endl;
    const cm C = A * B;
    cout << "A*B: "         << endl << C << endl;

    // Division
    cout << "B/B: "         << endl << (B / B) << endl;
    cout << "(B/B)*B: "     << endl << (B / B) * B << endl;
    cout << "(A*B)/(A*B): " << endl << (C / C) << endl;
    cout << "inv(A): "      << endl << inv(A) << endl;

    // Truncate
    const streamsize prec = cout.precision(9);
    cout << "(b/b).truncated(): " << endl
    <<  (b/b).truncated()    << endl;
    cout << "((b/b)*b).truncated(): " << endl
    <<  ((b/b)*b).truncated()    << endl;
    cout << "((a*b)/(a*b)).truncated(): " << endl
    << (c/c).truncated() << endl;
    cout << "inv(a).truncated(): " << endl
    <<  inv(a).truncated()    << endl;
    cout.precision(prec);
}

int test13();

#endif
