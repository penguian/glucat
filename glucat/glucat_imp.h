#ifndef _GLUCAT_GLUCAT_IMP_H
#define _GLUCAT_GLUCAT_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    glucat_imp.h : Organize GluCat template definitions which cannot be compiled separately
                             -------------------
    begin                : Sun 2001-12-25
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
     Arvind Raja's original header comments and references follow.
 ***************************************************************************
// clifford algebra package,  Arvind.Raja@hut.fi
// ref: Press et.al. "Numerical Recipes in C", 2nd ed., C.U.P., 1992.
// ref: LEDA, v 3.0, Stefan N\"aher, Max-Planck-Institut f\"ur Informatik
// ref: Stroustrup B., "The C++ Programming Language", 2nd ed.,
//      Addison-Wesley, 1991.
// ref: R. Sedgewick, "Algorithms in C++", Addison-Wesley, 1992.
// ref: S. Meyers, "Effective C++ ", Addison-Wesley, 1992.
 ***************************************************************************/

// Template definitions which cannot be compiled separately
#include "glucat/errors_imp.h"

#include <sstream>
#include "glucat/index_set_imp.h"
#include "glucat/clifford_algebra_imp.h"

#include <mtl/norm.h>
#include "glucat/framed_multi_imp.h"

#include <mtl/lu.h>
#include "glucat/matrix_multi_imp.h"

#include "glucat/matrix_imp.h"

#endif  // _GLUCAT_GLUCAT_IMP_H
