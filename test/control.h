#ifndef _GLUCAT_CONTROL_H
#define _GLUCAT_CONTROL_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    control.h : Define and set parameters to control tests
                             -------------------
    begin                : 2010-04-21
    copyright            : (C) 2010 by Paul C. Leopardi
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
#include "config.h"
namespace glucat
{
  /// Parameters to control tests
  struct control_t
  {
    /// Test parameters are valid
    bool m_valid;
    /// Produce more detailed output from tests
    bool m_verbose_output;
    /// Catch exceptions
    bool m_catch_exceptions;

    /// Default constructor
    control_t()
    : m_valid(false), m_verbose_output(false), m_catch_exceptions(true)
    {};
    /// Constructor from program arguments
    control_t(int argc, char ** argv);
  };

  /// Test control constructor from program arguments
  control_t::
  control_t(int argc, char ** argv)
  : m_valid(true), m_verbose_output(false), m_catch_exceptions(true)
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
      std::cout << program_name << " for GluCat version " << VERSION << ":" << std::endl;
      std::cout << "Usage:  " << program_name << " [option ...]" << std::endl;
      std::cout << "Options:" << std::endl;
      std::cout << "  --verbose   : Produce more detailed test output." << std::endl;
      std::cout << "  --help      : Print this summary." << std::endl;
    }
  }

  /// Global variable to control this test_control
  control_t test_control = control_t();
}
#endif // _GLUCAT_CONTROL_H
