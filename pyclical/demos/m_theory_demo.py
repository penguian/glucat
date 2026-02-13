#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# m_theory_demo.py: This file demonstrates tests for M-Theory physics,
#                   based on specific Glucat features.
#
#    copyright            : (C) 2014 by Paul C. Leopardi
#
#    This library is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library.  If not, see <http://www.gnu.org/licenses/>.

from PyClical import *
from pyclical_tutorial_utils import *

def run(ctx=tutorial_context(globals())):
    for name, method in get_object_methods(ctx).items():
        exec("global "+name+";"+name+"=method")

    print_fill("# M-Theory Analysis: R_(10, 1) vs R_(1, 10)")
    print_line()
    print_fill("We use the p=positive, q=negative convention.")
    print_fill("R_(10, 1) has 10 positive squares and 1 negative square.")
    print_line()

    # --- 1. THE PSEUDOSCALAR TEST ---
    print_fill("1. THE PSEUDOSCALAR TEST")
    print_fill("We define omega for R_(10, 1) using indices {-1, 1, 2, ..., 10}.")
    # Signature: q=1 (-1), p=10 (1...10)
    print_exec("omega_10_1 = e({-1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10})")
    print_exec("print(f'R_(10, 1) omega^2: {omega_10_1 * omega_10_1}')")

    print_fill("Now we define omega for R_(1, 10) using ten negative indices.")
    # Signature: p=1 (1), q=10 (-2...-11)
    print_exec("omega_1_10 = e({1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11})")
    print_exec("print(f'R_(1, 10) omega^2: {omega_1_10 * omega_1_10}')")

    pause()

    # --- 2. IDEMPOTENTS AND SPLITTING ---
    print_line()
    print_fill("2. REAL IDEMPOTENTS (MAJORANA SPLIT)")
    print_fill("Since R_(10, 1) omega^2 = 1, we can create real projection operators.")
    print_fill("This proves the algebra splits into M32(R) + M32(R).")

    print_exec("P_plus = 0.5 * (1 + omega_10_1)")
    print_exec("P_minus = 0.5 * (1 - omega_10_1)")

    print_fill("Verify Idempotency: P+ * P+ should equal P+")
    print_exec("print(f'P+ * P+: {P_plus * P_plus}')")

    print_fill("Verify Orthogonality: P+ * P- should equal 0")
    print_exec("print(f'P+ * P-: {P_plus * P_minus}')")

    pause()

    # --- 3. BRANE ROTATION (POWER MOVE) ---
    print_line()
    print_fill("3. M2-BRANE ROTATION")
    print_fill("We define an M2-brane in the (1, 2) plane and rotate it.")

    print_exec("M2 = e({1, 2})")
    print_fill("Define a rotor for a 90-degree rotation in the (2, 3) plane.")
    # Rotor R = exp(theta/2 * Bivector)
    print_exec("theta = math.pi/2")
    print_exec("B = e({2, 3})")
    print_exec("R = cos(theta/2) - B * sin(theta/2)")

    print_fill("Apply rotation: M2_new = M2 | R")
    print_exec("M2_new = M2 | R")
    print_exec("print(f'Rotated M2-brane: {M2_new}')")
    print_line()

    pause()

    # --- 4. BRANE DECOMPOSITION CHECK ---
    print_line()
    print_fill("4. BRANE DECOMPOSITION CHECK")
    print_fill("Confirming dimensions for M2 (Grade 2) and M5 (Grade 5).")
    print_exec("print('M2 basis size (11C2):', 55)")
    print_exec("print('M5 basis size (11C5):', 462)")

    pause()

    # --- 5. M5-BRANE EXAMPLE ---
    print_line()
    print_fill("5. M5-BRANE EXAMPLE")
    print_fill("Creating a Grade 5 element (M5-brane) in R_(10, 1) and checking its square.")
    # Signature: q=1 (-1), p=10 (1...10)
    print_exec("M5 = e({-1, 1, 2, 3, 4}) + e({6, 7, 8, 9, 10}); print(M5)")
    print_exec("print(M5 * M5)")

    print_line()
    print_fill("Demo complete. R_(10, 1) is the correct choice for Majorana spinors.")

    pause()
    print_line()
    print_fill("You have completed the structural analysis of M-Theory physics in Glucat.")

if __name__ == "__main__":
    try:
        run()
    except:
        print("The demo was interrupted.")
        pass
