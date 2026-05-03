#ifndef _GLUCAT_ERRORS_IMP_H
#define _GLUCAT_ERRORS_IMP_H
/**************************************************************************
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
  /*
   * @brief Specific exception with an error message
   * @details
   * @tparam Class_T Base exception class
   * @param msg Error message
   */
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
  heading() const noexcept -> std::string_view
  { return "Error in glucat::"; }

  template< class Class_T >
  auto
  error<Class_T>::
  classname() const noexcept -> std::string_view
  { return name; }

  template< class Class_T >
  void
  error<Class_T>::
  print_error_msg() const
  { std::cerr << heading() << classname() << std::endl << what() << std::endl; }
}

#ifdef GLUCAT_DOCTEST
#include <sstream>
#include <string_view>
struct dummy_class {
  static constexpr std::string_view name = "dummy_class";
  static std::string_view classname() { return name; }
};

TEST_CASE("errors") {
  using namespace glucat;

  SUBCASE("Two-argument constructor") {
    error<dummy_class> e("test_context", "test_message");
    CHECK(e.heading() == "Error in glucat::");
    CHECK(e.classname() == "test_context");
    CHECK(std::string(e.what()) == "test_message");
  }

  SUBCASE("One-argument constructor") {
    error<dummy_class> e("test_message_only");
    CHECK(e.classname() == "dummy_class");
    CHECK(std::string(e.what()) == "test_message_only");
  }

  SUBCASE("Output verification (std::cerr redirection)") {
    error<dummy_class> e("context", "message");

    // Redirect cerr to a stringstream using RAII
    std::stringstream buffer;
    struct CerrRedirect {
      std::streambuf* old;
      CerrRedirect(std::streambuf* new_buf) : old(std::cerr.rdbuf(new_buf)) {}
      ~CerrRedirect() { std::cerr.rdbuf(old); }
    } redirect(buffer.rdbuf());

    e.print_error_msg();

    std::string output = buffer.str();
    CHECK(output.find("Error in glucat::") != std::string::npos);
    CHECK(output.find("context") != std::string::npos);
    CHECK(output.find("message") != std::string::npos);
  }
}
#endif

#endif // _GLUCAT_ERRORS_IMP_H
