#ifndef _GLUCAT_ERRORS_IMP_H
#define _GLUCAT_ERRORS_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    errors_imp.h : Define error functions
                             -------------------
    begin                : Sun 2001-12-20
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

#include "glucat/errors.h"

#include <string>
#include <iostream>
#include <ostream>

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
  auto
  error<Class_T>::
  heading() const noexcept -> const std::string
  { return "Error in glucat::"; }

  template< class Class_T >
  auto
  error<Class_T>::
  classname() const noexcept -> const std::string
  { return name; }

  template< class Class_T >
  void
  error<Class_T>::
  print_error_msg() const
  { std::cerr << heading() << classname() << std::endl << what() << std::endl; }
}
#endif // _GLUCAT_ERRORS_IMP_H
