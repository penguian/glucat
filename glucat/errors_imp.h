#ifndef _GLUCAT_ERRORS_IMP_H
#define _GLUCAT_ERRORS_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    errors_imp.h : Define error functions
                             -------------------
    begin                : Sun 2001-12-20
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

namespace glucat
{
  /// Specific exception class
  template< class Class_T >
  error<Class_T>::
  error(const std::string& msg)
  : glucat_error(Class_T::classname(), msg)
  { }

  template< class Class_T >
  error<Class_T>::
  error(const std::string& context, const std::string& msg)
  : glucat_error(context, msg)
  { }

  template< class Class_T >
  const std::string
  error<Class_T>::
  heading() const throw()
  { return "Error in glucat::"; }

  template< class Class_T >
  const std::string
  error<Class_T>::
  classname() const throw()
  { return name; }

  template< class Class_T >
  void
  error<Class_T>::
  print_error_msg() const
  { std::cerr << heading() << classname() << std::endl << what() << std::endl; }
}
#endif // _GLUCAT_ERRORS_IMP_H
