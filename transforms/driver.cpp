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

#include "transforms/driver.h"
using namespace glucat;

int main(int argc, char ** argv)
{
  using namespace std;
  for (argc--, argv++; argc != 0; argc--, argv++)
  {
    long int max_index = 0;
    sscanf(*argv, "%ld", &max_index);
    if (max_index > 0)
      try
        { transforms(max_index); }
      catch (const glucat_error& e)
        { e.print_error_msg(); }
      catch (bad_alloc)
        { cerr << "bad_alloc" << endl; }
      catch (...)
        { cerr << "unexpected exception" << endl; }
  }
  return 0;
}
