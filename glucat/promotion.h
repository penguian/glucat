#ifndef _GLUCAT_PROMOTION_H
#define _GLUCAT_PROMOTION_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    promotion.h : Define promotion and demotion for specific scalar types
                             -------------------
    begin                : 2021-11-13
    copyright            : (C) 2021 by Paul C. Leopardi
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
 "Clifford algebras with numeric and symbolic computations, Birkhauser, 1996."
 ***************************************************************************
 See also Arvind Raja's original header comments and references in glucat.h
 ***************************************************************************/

#include "glucat/global.h"
#include "glucat/scalar.h"
#include "glucat/qd.h"

#include <cfloat>
#include <limits>

#if defined(_GLUCAT_USE_QD)
# include <qd/qd_real.h>
#endif

namespace glucat
{
  /// Extra traits which extend numeric limits
  // Reference: [AA], 2.4, p. 30-31

#if !defined(_GLUCAT_USE_QD) || !defined(QD_API)

# if DBL_MANT_DIG < LDBL_MANT_DIG

  /// Promoted type for double
  template<>
  struct
  numeric_traits<double>::
  promoted {using type = long double;};

  /// Demoted type for long double
  template<>
  struct
  numeric_traits<long double>::
  demoted {using type = double;};

# else

  /// Promoted type for double
  template<>
  struct
  numeric_traits<double>::
  promoted {using type = double;};

  /// Demoted type for long double
  template<>
  struct
  numeric_traits<long double>::
  demoted {using type = float;};

# endif // DBL_MANT_DIG < LDBL_MANT_DIG

  /// Promoted type for long double
  template<>
  struct
  numeric_traits<long double>::
  promoted {using type = long double;};

#else

# if (DBL_MANT_DIG < LDBL_MANT_DIG) && (LDBL_MANT_DIG < DBL_MANT_DIG*2)

  /// Promoted type for double
  template<>
  struct
  numeric_traits<double>::
  promoted {using type = long double;};

  /// Demoted type for long double
  template<>
  struct
  numeric_traits<long double>::
  demoted {using type = double;};

  /// Promoted type for long double
  template<>
  struct
  numeric_traits<long double>::
  promoted {using type = dd_real;};

  /// Demoted type for dd_real
  template<>
  struct
  numeric_traits<dd_real>::
  demoted {using type = long double;};

  /// Promoted type for dd_real
  template<>
  struct
  numeric_traits<dd_real>::
  promoted {using type = qd_real;};

  /// Demoted type for qd_real
  template<>
  struct
  numeric_traits<qd_real>::
  demoted {using type = dd_real;};

# elif (LDBL_MANT_DIG < DBL_MANT_DIG*2)

  /// Promoted type for double
  template<>
  struct
  numeric_traits<double>::
  promoted {using type = dd_real;};

  /// Demoted type for long double
  template<>
  struct
  numeric_traits<long double>::
  demoted {using type = float;};

  /// Promoted type for long double
  template<>
  struct
  numeric_traits<long double>::
  promoted {using type = dd_real;};

  /// Demoted type for dd_real
  template<>
  struct
  numeric_traits<dd_real>::
  demoted {using type = double;};

  /// Promoted type for dd_real
  template<>
  struct
  numeric_traits<dd_real>::
  promoted {using type = qd_real;};

  /// Demoted type for qd_real
  template<>
  struct
  numeric_traits<qd_real>::
  demoted {using type = dd_real;};

# else

  /// Promoted type for double
  template<>
  struct
  numeric_traits<double>::
  promoted {using type = dd_real;};

  /// Demoted type for dd_real
  template<>
  struct
  numeric_traits<dd_real>::
  demoted {using type = double;};

  /// Promoted type for dd_real
  template<>
  struct
  numeric_traits<dd_real>::
  promoted {using type = long double;};

  /// Demoted type for long double
  template<>
  struct
  numeric_traits<long double>::
  demoted {using type = dd_real;};

  /// Promoted type for long double
  template<>
  struct
  numeric_traits<long double>::
  promoted {using type = qd_real;};

  /// Demoted type for qd_real
  template<>
  struct
  numeric_traits<qd_real>::
  demoted {using type = long double;};

# endif // (DBL_MANT_DIG < LDBL_MANT_DIG) && (LDBL_MANT_DIG < DBL_MANT_DIG*2)

  /// Promoted type for qd_real
  template<>
  struct
  numeric_traits<qd_real>::
  promoted {using type = qd_real;};

#endif // !defined(_GLUCAT_USE_QD) || !defined(QD_API)

} // namespace glucat

#endif // _GLUCAT_PROMOTION_H
