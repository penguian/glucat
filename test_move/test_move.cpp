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
