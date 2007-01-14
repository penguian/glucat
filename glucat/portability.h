#ifndef _GLUCAT_PORTABILITY_H
#define _GLUCAT_PORTABILITY_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    portability.h : Work around non-standard compilers and libraries
                             -------------------
    begin                : Sun 2001-08-18
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

//***************************** workarounds for ICC
#if defined (BOOST_INTEL) || defined (__ICL) || defined (__ICC)
# pragma warning( disable: 177 ) // variable was declared but never referenced
# pragma warning( disable: 279 ) // controlling expression is constant
# pragma warning( disable: 383 ) // value copied to temporary, reference to temporary ...
# pragma warning( disable: 444 ) // destructor for base is not virtual
# pragma warning( disable: 593 ) // variable was set but never used
# pragma warning( disable: 810 ) // conversion from "double" to "int" may lose significant bits
# pragma warning( disable: 869 ) // parameter was never referenced
# pragma warning( disable: 981 ) // operands are evaluated in unspecified order
# pragma warning( disable: 1572 ) // floating-point equality and inequality comparisons ...
#endif

#endif // _GLUCAT_PORTABILITY_H
