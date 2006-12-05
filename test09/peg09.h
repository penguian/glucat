#ifndef GLUCAT_TEST_PEG09_H
#define GLUCAT_TEST_PEG09_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg09.cpp : programming example 9 : vector_part
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

#include <iomanip>

namespace peg09
{
  using namespace glucat;
  using namespace std;
  template< class Multivector_T >
  static
  void 
  do_test9()
  {
    typedef Multivector_T number;
    typedef typename number::vector_t vector_t;
    typedef typename number::index_set_t index_set_t;

    const number a = number("3+{-1}+3{4}+4{5}+2{1,3}+{-1}+9{1,3,6}");
    a.write("a =");
    const index_set_t a_frm = a.frame();
    cout << "a_frm = a.frame() = " << a_frm << endl;
    cout << "a_frm.count() = " << a_frm.count() << endl;
    const vector_t a_vec = a.vector_part();
    cout << "a_vec = a.vector_part() =" << endl;
    index_t idx = a_frm.min();
    const index_t a_frm_end = a_frm.max()+1;
    for (typename vector_t::const_iterator scvec = a_vec.begin();
         scvec != a_vec.end(); ++scvec)
      {
        while(idx != a_frm_end && !a_frm[idx])
          ++idx;
        cout << "[" << idx++ << "] " << *scvec << endl;
      }
    cout << "a_vec.size() = " << a_vec.size() << endl;
    const number A = number(a_vec, a_frm);
    A.write("A = number(a_vec, a_frm) =");
    (a(1)).write("a(1)");
    cout << endl;
  }
}

int test09();

#endif
