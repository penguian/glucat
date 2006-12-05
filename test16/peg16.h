#ifndef GLUCAT_TEST_PEG16_H
#define GLUCAT_TEST_PEG16_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg16.cpp : programming example 16 : Matrices of multivectors
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001 by Paul C. Leopardi
    email                : leopardi@bigpond.net.au
 ***************************************************************************
 *   This library is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2.1 of the  *
 *   License, or (at your option) any later version.                       *
 *   See http://www.fsf.org/copyleft/lesser.html for details               *
 ***************************************************************************
 This library is based on a prototype written by Arvind Raja and was
 licensed under the LGPL with permission of the author. See Arvind Raja,
 "Object-oriented implementations of Clifford algebras in C++: a prototype",
 in Ablamowicz, Lounesto and Parra (eds.)
 "Clifford algebras with numeric and symbolic computations", Birkhauser, 1996.
 ***************************************************************************
     See also Arvind Raja's original header comments in glucat.h
 ***************************************************************************/

namespace peg16
{
  using namespace glucat;
  using namespace std;

  template< class Multivector_T >
  static
  void 
  do_test16()
  {
    typedef Multivector_T number;
    typedef typename number::index_set_t index_set_t;
    const index_t hi = index_set_t::v_hi;
    typedef ublas::matrix<number> matrix_of_mv;
    const int dim = 5;
    matrix_of_mv a(dim,dim);
    for (int i = 0; i != dim; ++i)
      for (int j = 0; j != dim; ++j)
        a(i,j) = number("0");
    a(1,1) = number("1");
    std::cout << a << std::endl;
    matrix_of_mv b(dim,dim);
    for (int i = 0; i != dim; ++i)
    {
      for (int j = 0; j != dim; ++j)
        b(i,j) = number("0");
      b(i,i) = number(index_set_t((i % hi) +1),1.0);
    }
    std::cout << b << std::endl;
    matrix_of_mv c(dim,dim);
    c = ublas::prod(a, b);
    std::cout << c << std::endl;
    c += b;
    std::cout << c << std::endl;
  }
}

int test16();

#endif
