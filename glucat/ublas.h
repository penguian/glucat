#ifndef _GLUCAT_UBLAS_H
#define _GLUCAT_UBLAS_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    ublas.h : Define special uBLAS functions
                             -------------------
    begin                : Sun 2002-12-21
    copyright            : (C) 2002 by Paul C. Leopardi
                         : uBLAS interface contributed by Joerg Walter
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
     For Arvind Raja's original header comments, see glucat.h
 ***************************************************************************/

namespace boost { namespace numeric { namespace ublas {

    template<class T, class A = unbounded_array<T> >
    class permutation_matrix:
        public vector<T, A> {
    public:
        typedef typename vector<T, A>::size_type size_type;

        BOOST_UBLAS_EXPLICIT
        BOOST_UBLAS_INLINE
        permutation_matrix (size_type size):
            vector<T, A> (size) {}
    };

    template<class PM, class MV>
    void 
    swap_rows (const PM &pm, MV &mv, vector_tag);

    template<class PM, class MV>
    void 
    swap_rows (const PM &pm, MV &mv, matrix_tag);

    // Dispatcher
    template<class PM, class MV>
    void
    swap_rows (const PM &pm, MV &mv);

    template<class M, class PM>
    typename M::size_type
    lu_factorize (M &m, PM &pm);

    template<class M, class PM, class MV>
    void
    lu_solve (M &m, const PM &pm, MV &mv);

} } }

namespace ublas = boost::numeric::ublas;

#endif  // _GLUCAT_UBLAS_H
