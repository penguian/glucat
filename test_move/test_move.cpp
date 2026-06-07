/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    test_move.cpp : Test driver for move constructor and move assignment
                             -------------------
    begin                : Wed May 06 2026
    copyright            : (C) 2026 by Paul C. Leopardi
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

 ***************************************************************************/

#include <iostream>
#include <type_traits>
#include "test/driver.h"

template <typename T>
void check_move_ability(const std::string& class_name) {
    std::cout << "Checking " << class_name << ":" << std::endl;
    std::cout << "  is_move_constructible: " << std::is_move_constructible<T>::value << "" << std::endl;
    std::cout << "  is_move_assignable:    " << std::is_move_assignable<T>::value << "" << std::endl;

    if (std::is_move_constructible<T>::value && std::is_move_assignable<T>::value) {
        std::cout << "  [PASS] Move semantics enabled." << std::endl;
    } else {
        std::cout << "  [FAIL] Move semantics NOT enabled." << std::endl;
    }
    std::cout << "----------------------------------------" << std::endl;
}

int main() {
    using namespace glucat;

    std::cout << "Move semantics test:" << std::endl;
    std::cout << std::endl;

    check_move_ability< index_set<DEFAULT_LO,DEFAULT_HI> >("index_set<DEFAULT_LO,DEFAULT_HI>");
    check_move_ability< framed_multi<double> >("framed_multi<double>");
    check_move_ability< matrix_multi<double> >("matrix_multi<double>");

    return 0;
}
