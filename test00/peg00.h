#ifndef GLUCAT_TEST_PEG00_H
#define GLUCAT_TEST_PEG00_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg00.cpp : programming example 00 : Geometric algebra identities
                             -------------------
    begin                : Sat 2007-09-01
    copyright            : (C) 2007-2012 by Paul C. Leopardi
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

#include <iomanip>

namespace peg00
{
  // References for identities:
  // [D01]: Leo Dorst, "Honing geometric algebra for its use in the computer sciences",
  // in Geometric Computing with Clifford Algebras, (G. Sommer, ed.)
  // Springer 2001, Chapter 6, pp. 127-152.
  // http://staff.science.uva.nl/~leo/clifford/index.html

  // [D02]: Leo Dorst, "The inner products of geometric algebra", in
  // Applications of Geometric Algebra in Computer Science and Engineering
  // (Dorst, Doran, Lasenby, eds), Birkhauser, 2002.
  // http://staff.science.uva.nl/~leo/clifford/index.html

  // [HS]: David Hestenes, Garret Sobczyk, Clifford Algebra to Geometric Calculus,
  // D. Reidel, 1984.

  using namespace glucat;

  template< typename Multivector_T >
  static
  bool
  is_error(const Multivector_T& lhs, const Multivector_T& rhs, typename Multivector_T::scalar_t tol)
  {
    return ( (abs(lhs) < tol) || (abs(rhs) < tol) )
           ? (abs(lhs - rhs) > tol)
           : (abs(lhs - rhs) > abs(rhs) * tol);
  }

  template< typename Multivector_T >
  static
  void
  print_error_lhs_rhs(const Multivector_T& lhs, const Multivector_T& rhs, typename Multivector_T::scalar_t tol)
  {
    std::cout
        << "  LHS == " << lhs << std::endl;
    std::cout
        << "  RHS == " << rhs << std::endl;
    std::cout
        << "  Relative norm of difference == " << abs(lhs - rhs)/abs(rhs) << std::endl;
    if (tol > 0.0)
      std::cout
        << "  Tolerance                   == " << tol << std::endl;
  }

  template< typename Multivector_T >
  static
  bool
  test_idents(Multivector_T& a, Multivector_T& b,
              Multivector_T& c, typename Multivector_T::index_set_t e)
  {
    typedef Multivector_T multivector_t;
    typedef typename multivector_t::scalar_t scalar_t;

    static const scalar_t eps = std::numeric_limits<scalar_t>::epsilon();
    const scalar_t tol = (control_t::verbose())
      ? 0.0
      : eps * 512.0;

    multivector_t lhs;
    multivector_t rhs;

    bool success = true;

    const index_t a_grade = a.frame().count();
    const index_t b_grade = b.frame().count();

    { // Identity [HS] (1.21a)
      for (index_t r = 1;
           r <= a_grade;
           ++r)
      {
        const multivector_t a_r = a(r);
        for (index_t s = 1;
             s <= b_grade;
             ++s)
        {
          const multivector_t b_s = b(s);
          lhs = a_r & b_s;
          rhs = (a_r * b_s)(index_t(std::abs(r-s)));

          if ( is_error(lhs, rhs, tol) )
          {
            std::cout << "Identity [HS] (1.21a) failed in " << e
                      << ": r == " << r
                      << ", s == " << s
                      << std::endl;
            print_error_lhs_rhs(lhs, rhs, tol);
            success = false;
          }
        }
      }
    }

    { // Identity [HS] (1.21b)
      for (index_t r = 0;
           r <= b_grade;
           ++r)
      {
        lhs = a(0) & b(r);
        rhs = a(r) & b(0);

        if ( (abs(lhs) > tol) || (abs(rhs) > tol) )
        {
          std::cout << "Identity [HS] (1.21b) failed in " << e
                    << ": grade == " << r
                    << std::endl;
          print_error_lhs_rhs(lhs, rhs, tol);
          success = false;
        }
      }
    }

    { // Identity [HS] (1.22a)
      for (index_t r = 0;
           r <= a_grade;
           ++r)
      {
        const multivector_t a_r = a(r);
        for (index_t s = 0;
             s <= b_grade;
             ++s)
        {
          const multivector_t b_s = b(s);
          lhs = a_r ^ b_s;
          rhs = (a_r * b_s)(r+s);

          if ( is_error(lhs, rhs, tol) )
          {
            std::cout << "Identity [HS] (1.22a) failed in " << e
                      << ": r == " << r
                      << ", s == " << s
                      << std::endl;
            print_error_lhs_rhs(lhs, rhs, tol);
            success = false;
          }
        }
      }
    }

    { // Identity [HS] (1.25a)
      lhs = (a ^ b) ^ c;
      rhs = a ^ (b ^ c);

      if ( is_error(lhs, rhs, tol) )
      {
        std::cout << "Identity [HS] (1.25a) failed in " << e
                  << std::endl;
        print_error_lhs_rhs(lhs, rhs, tol);
        success = false;
      }
    }

    { // Identity [HS] (1.31)
      const multivector_t a_1 = a(1);
      lhs = a_1 * b;
      rhs = (a_1 & b) + (a_1 ^ b);

      if ( is_error(lhs, rhs, tol) )
      {
        std::cout << "Identity [HS] (1.31) failed in " << e
                  << std::endl;
        print_error_lhs_rhs(lhs, rhs, tol);
        success = false;
      }
    }

    { // Identity [D01] (Section 2.3 Example 2, vector)
      const multivector_t a_1 = a(1);
      lhs = a_1 * b;
      rhs = (a_1 % b) + (a_1 ^ b);

      if ( is_error(lhs, rhs, tol) )
      {
        std::cout << "Identity [D01] (Section 2.3 Example 2, vector) failed in " << e
                  << std::endl;
            print_error_lhs_rhs(lhs, rhs, tol);
        success = false;
      }
    }

    { // Identity [HS] (1.63)
      const multivector_t a_2 = a(2);
      const multivector_t b_m_b_1 = b - b(1);
      lhs =  a_2 * b_m_b_1;
      rhs = (a_2 & b_m_b_1) + (a_2 * b_m_b_1 - b_m_b_1 * a_2)/scalar_t(2) + (a_2 ^ b_m_b_1);

      if ( is_error(lhs, rhs, tol) )
      {
        std::cout << "Identity [HS] (1.63) failed in " << e
                  << std::endl;
        print_error_lhs_rhs(lhs, rhs, tol);
        success = false;
      }
    }

    { // Identity [D01] (Section 2.3 Example 2, bivector)
      const multivector_t a_2 = a(2);
      lhs =  a_2 * b;
      rhs = (a_2 % b) + (a_2 * b - b * a_2)/scalar_t(2) + (a_2 ^ b);

      if ( is_error(lhs, rhs, tol) )
      {
        std::cout << "Identity [D01] (Section 2.3 Example 2, bivector) failed in " << e
                  << std::endl;
        print_error_lhs_rhs(lhs, rhs, tol);
        success = false;
      }
    }

    { // Identity [HS] (1.44)
      const scalar_t scalar_lhs = star(a, b);
      const scalar_t scalar_rhs = scalar(a * b);

      const scalar_t scalar_diff = scalar_lhs - scalar_rhs;
      if (scalar_diff*scalar_diff > tol*tol)
      {
        std::cout << "Identity [HS] (1.44) failed in " << e
                  << std::endl;
        print_error_lhs_rhs(multivector_t(scalar_lhs), multivector_t(scalar_rhs), tol);
        success = false;
      }
    }

    { // Identity [HS] (1.48)
      const scalar_t scalar_lhs = star(a, b);
      const scalar_t scalar_rhs = star(reverse(a), reverse(b) );

      const scalar_t scalar_diff = scalar_lhs - scalar_rhs;
      if (scalar_diff*scalar_diff > tol*tol)
      {
        std::cout << "Identity [HS] (1.48) failed in " << e
                  << std::endl;
        print_error_lhs_rhs(multivector_t(scalar_lhs), multivector_t(scalar_rhs), tol);
        success = false;
      }
    }

    { // Identity [D01] (2.5) (a.2)
      const multivector_t a_0 = a(0);
      lhs = a_0 % b;
      rhs = a_0 * b;

      if ( is_error(lhs, rhs, tol) )
      {
        std::cout << "Identity [D01] (2.5) (a.2) failed in " << e
                  << std::endl;
        print_error_lhs_rhs(lhs, rhs, tol);
        success = false;
      }
    }

    { // Identity [D01] (2.5) (a.3)
      const multivector_t a_1 = a(1);
      const multivector_t b_0 = b(0);
      lhs = a_1 % b_0;
      rhs = 0;

      if ( is_error(lhs, rhs, tol) )
      {
        std::cout << "Identity [D01] (2.5) (a.3) failed in " << e
                  << std::endl;
        print_error_lhs_rhs(lhs, rhs, tol);
        success = false;
      }
    }

    { // Identity [D01] (2.5) (c)
      const multivector_t a_1 = a(1);
      lhs = a_1 % (b ^ c);
      rhs = ( (a_1 % b) ^ c) + (involute(b) ^ (a_1 % c) );

      if ( is_error(lhs, rhs, tol) )
      {
        std::cout << "Identity [D01] (2.5) (c) failed in " << e
                  << std::endl;
        print_error_lhs_rhs(lhs, rhs, tol);
        success = false;
      }
    }

    { // Identity [D01] (2.5) (d), [D02] (2.12),
      lhs = (a ^ b) % c;
      rhs = a % (b % c);

      if ( is_error(lhs, rhs, tol) )
      {
        std::cout << "Identity [D01] (2.5) (d) failed in " << e
                  << std::endl;
        print_error_lhs_rhs(lhs, rhs, tol);
        success = false;
      }
    }
    return success;
  }

  template< class Multivector_T >
  static
  void
  do_test00(const index_t max_index)
  {
    typedef Multivector_T multivector_t;
    typedef typename multivector_t::index_set_t index_set_t;
    typedef typename multivector_t::scalar_t scalar_t;

    bool success = true;

    scalar_t fill = 0.5;
    index_set_t frm = index_set_t();
    for (index_t i = 1; i != max_index+1; i++)
    {
      frm |= index_set_t(i);
      multivector_t a = multivector_t::random(frm, fill);
      multivector_t b = multivector_t::random(frm, fill);
      multivector_t c = multivector_t::random(frm, fill);
      success &= test_idents(a, b, c, frm);
      frm |= index_set_t(-i);
      a = multivector_t::random(frm, fill);
      b = multivector_t::random(frm, fill);
      c = multivector_t::random(frm, fill);
      success &= test_idents(a, b, c, frm);
    }
    if (success)
      std::cout << "All tests passed." << std::endl;
  }
}

int test00();

#endif
