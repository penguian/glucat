#ifndef _GLUCAT_PORTABILITY_H
#define _GLUCAT_PORTABILITY_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    portability.h : Work around non-standard compilers and libraries
                             -------------------
    begin                : Sun 2001-08-18
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

//***************************** fix to `ios_base':
#if defined (__GNUG__) && (__GNUC__ <= 2) && (__GNUC_MINOR__ <= 96)
#define ios_base ios // scope is different for standard C++
#endif

//***************************** workaround for ICC and G++ 3.3+
#if defined (__ICL) || defined (__ICC) || defined (__GNUG__) && (__GNUC__ >= 3) && (__GNUC_MINOR__ >= 3)
#define _GLUCAT_PRIVATE public
#else
#define _GLUCAT_PRIVATE private
#endif

//***************************** workaround for G++ 3.2 typename bug
#if defined (__GNUG__) && (__GNUC__ >= 3) && (__GNUC_MINOR__ <= 2)
#define _GLUCAT_USE_STRUCT_NAME(T)
#else
#define _GLUCAT_USE_STRUCT_NAME(T) T::
#endif

#endif // _GLUCAT_PORTABILITY_H
