#ifndef _GLUCAT_TRY_CATCH_H
#define _GLUCAT_TRY_CATCH_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    errors_imp.h : Define error functions
                             -------------------
    begin                : Sun 2001-12-20
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

namespace glucat
{
  /// For exception catching: pointer to function returning int
  typedef int (*intfn)();
  
  /// For exception catching: pointer to function of int returning int
  typedef int (*intintfn)(int);
  
  /// Exception catching for functions returning int
  int try_catch(intfn f)
  {
    int result = 0;
    try
      { result = (*f)(); }
    catch (const glucat_error& e)
      { e.print_error_msg(); }
    catch (std::bad_alloc)
      { std::cerr << "bad_alloc" << std::endl; }
    catch (...)
      { std::cerr << "unexpected exception" << std::endl; }
    return result;
  }
  
  /// Exception catching for functions of int returning int
  int try_catch(intintfn f, int arg)
  {
    int result = 0;
    try
      { result = (*f)(arg); }
    catch (const glucat_error& e)
      { e.print_error_msg(); }
    catch (std::bad_alloc)
      { std::cerr << "bad_alloc" << std::endl; }
    catch (...)
      { std::cerr << "unexpected exception" << std::endl; }
    return result;
  }
}
#endif // _GLUCAT_TRY_CATCH_H
