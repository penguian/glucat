#ifndef _GLUCAT_CONTROL_H
#define _GLUCAT_CONTROL_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    control.h : Define and set parameters to control tests
                             -------------------
    begin                : 2010-04-21
    copyright            : (C) 2010-2016 by Paul C. Leopardi
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
#include "glucat/glucat_config.h"
#include "test/try_catch.h"

namespace glucat
{
  /// Parameters to control tests
  class control_t
  {
  private:
    /// Test parameters are valid
    bool m_valid;
    bool valid() const
    { return m_valid; }

    /// Catch exceptions
    bool m_catch_exceptions;
    bool catch_exceptions() const
    { return m_catch_exceptions; }

    /// Produce more detailed output from tests
    static bool m_verbose_output;
    
    /// Constructor from program arguments
    control_t(int argc, char ** argv);
    // Enforce singleton
    // Reference: A. Alexandrescu, "Modern C++ Design", Chapter 6
    control_t() = default;
    ~control_t() = default;
    control_t(const control_t&);
    control_t& operator= (const control_t&);

    /// Friend declaration to avoid compiler warning:
    /// "... only defines a private destructor and has no friends"
    /// Ref: Carlos O'Ryan, ACE http://doc.ece.uci.edu
    friend class friend_for_private_destructor;
  public:
    /// Single instance
    /// Ref: Scott Meyers, "Effective C++" Second Edition, Addison-Wesley, 1998.
    static const control_t& control(int argc, char ** argv)
    { static const control_t c(argc, argv); return c; }
    
    /// Call a function that returns int
    int call(intfn f) const;
    /// Call a function of int that returns int
    int call(intintfn f, int arg) const;

    /// Produce more detailed output from tests
    static bool verbose()
    { return m_verbose_output; }
  };

  /// Produce more detailed output from tests
  bool control_t::m_verbose_output = false;
  
  /// Test control constructor from program arguments
  control_t::
  control_t(int argc, char ** argv)
  : m_valid(true), m_catch_exceptions(true)
  {
    bool print_help = false;
    const std::string& arg_0_str = argv[0];
    const std::string program_name = arg_0_str.substr(arg_0_str.find_last_of('/')+1);
    for (int arg_ndx = 1; arg_ndx < argc; ++arg_ndx)
    {
      const std::string& arg_str = argv[arg_ndx];
      bool valid = false;
      if (arg_str.substr(0,2) == "--")
      {
        valid = true;
        const std::string& arg_name = arg_str.substr(2);
        if (arg_name == "help")
        {
          this->m_valid = false;
          print_help = true;
        }
        else if (arg_name == "verbose")
          this->m_verbose_output = true;
        else if (arg_name == "no-catch")
          this->m_catch_exceptions = false;
        else
          valid = false;
      }
      if (!valid)
      {
        std::cout << "Invalid argument: " << arg_str << std::endl;
        this->m_valid = false;
        print_help = true;
      }
    }
    if (print_help)
    {
      std::cout << program_name << " for " << GLUCAT_PACKAGE_NAME << " version " << GLUCAT_VERSION << ":" << std::endl;
      std::cout << "Usage:  " << program_name << " [option ...]" << std::endl;
      std::cout << "Options:" << std::endl;
      std::cout << "  --help      : Print this summary." << std::endl;
      std::cout << "  --no-catch  : Do not catch exceptions." << std::endl;
      std::cout << "  --verbose   : Produce more detailed test output." << std::endl;
    }
  }

  /// Call a function that returns int
  inline
  int 
  control_t::
  call(intfn f) const
  {
    if (valid())
      return (catch_exceptions())
        ? try_catch(f)
        : (*f)();
    else
      return 1;
  }
      
  /// Call a function of int that returns int
  inline
  int 
  control_t::
  call(intintfn f, int arg) const
  {
    if (valid())
      return (catch_exceptions())
        ? try_catch(f, arg)
        : (*f)(arg);
    else
      return 1;
  }
}
#endif // _GLUCAT_CONTROL_H
