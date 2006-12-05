/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    squaring.cpp : Timing test driver
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

#include "squaring/squaring.h"

int squaring(const int n)
{
  using namespace glucat_mult_test;
  if (n > max_n)
  {
    cout << "Value " << n << " is too big." << endl;
    cout << "Maximum value allowed is " << max_n << "." << endl;
    return 1;
  }
  cout << "framed_multi<double>" << endl;
  mult_test< framed_multi<double> >(n, max_n);
  cout << "matrix_multi<double>" << endl;
  mult_test< matrix_multi<double> >(n, max_n);

  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{
  using namespace std;
  for (argc--, argv++; argc != 0; argc--, argv++)
  {
    int n = 0;
    sscanf(*argv, "%d", &n);
    try_catch(squaring, n);
  }
  return 0;
}
