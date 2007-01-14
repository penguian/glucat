/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    transforms.cpp : Example and timing test driver
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2007 by Paul C. Leopardi
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

#include "transforms/transforms.h"

int transforms(const int n)
{
  using namespace glucat_fast_test;
  if (n > max_n)
  {
    cout << "Value " << n << " is too big." << endl;
    cout << "Maximum value allowed is " << max_n << "." << endl;
    return 1;
  }
  cout << "framed_multi<double>" << endl;
  fast_test< framed_multi<double> >(n, max_n);
  cout << "matrix_multi<double>" << endl;
  fast_test< matrix_multi<double> >(n, max_n);
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{
  using namespace std;
  for (argc--, argv++; argc != 0; argc--, argv++)
  {
    int max_index = 0;
    sscanf(*argv, "%d", &max_index);
    if (max_index > 0)
      try_catch( transforms, max_index );
  }
  return 0;
}
