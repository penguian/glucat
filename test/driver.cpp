/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    driver.cpp : Example and timing test driver
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

#include "test/driver.h"
using namespace glucat;

int main(int argc, char ** argv)
{
  using namespace std;
  typedef int (*intfn)();
  intfn testfn[] =
  {
    test1, test2, test3, test4, test5, test6, test7, test8,
    test9, test10,test11,test12,test13,test14,test15,test16
  };
  const long int NBRTESTS = sizeof(testfn)/sizeof(intfn *);
  bitset<NBRTESTS> tests;
  if (argc == 1)
    tests.set();
  for (argc--, argv++; argc != 0; argc--, argv++)
  {
    long int test_nbr = 0;
    sscanf(*argv, "%ld", &test_nbr);
    if (test_nbr > 0 && test_nbr <= NBRTESTS)
      tests.set(test_nbr-1);
    else
      try
        { squaring(-test_nbr); }
      catch (const glucat_error& e)
        { e.print_error_msg(); }
  }
  for (index_t test_nbr = 0; test_nbr != NBRTESTS; ++test_nbr)
    if (tests[test_nbr])
    {
      cout << endl << "Test " << test_nbr+1 << ":" << endl;
      try
        { testfn[test_nbr](); }
      catch (const glucat_error& e)
        { e.print_error_msg(); }
      catch (bad_alloc)
        { cerr << "bad_alloc" << endl; }
      catch (...)
        { cerr << "unexpected exception" << endl; }
    }
  return 0;
}
