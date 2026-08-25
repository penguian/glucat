#!/usr/bin/env python3
"""
Multiparticle Spacetime Algebra (MSTA) 1, 2, and 3-Particle Simulation Demo.
"""
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# msta_demo.py: This file demonstrates Multiparticle Spacetime Algebra (MSTA)
#               for 1, 2, and 3-particle quantum systems, showing frame setups,
#               idempotents, entangled states, and higher-grade symmetries.
#
#    copyright            : (C) 2026 by Paul C. Leopardi
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

import math
from PyClical import *
from pyclical_tutorial_utils import *


def run(ctx=tutorial_context(globals())):
    """
    Run the Multiparticle Spacetime Algebra (MSTA) 1, 2, and 3-particle demo.
    """
    for name, method in get_object_methods(ctx).items():
        globals()[name] = method

    print_head("Multiparticle Spacetime Algebra (MSTA) Demo")
    print_line()
    print_fill(
        "This demo illustrates Multiparticle Spacetime Algebra (MSTA) for 1, 2,"
        + " and 3 particles in PyClical / GluCat."
    )
    print_fill(
        "MSTA assigns a 4D Minkowski spacetime frame G_(1,3) to each particle."
    )
    print_line()

    # =========================================================================
    # SECTION 1: 1-PARTICLE MSTA (n = 1, Signature: istpq(1, 3))
    # =========================================================================
    print_head("1. Single-Particle MSTA Space (n = 1, Frame: istpq(1, 3))")
    print_line()
    print_fill(
        "For 1 particle (1 qubit), the spacetime frame has 1 timelike generator"
        + " e(1) (e^2 = +1) and 3 spacelike generators e(-1), e(-2), e(-3) (e^2 = -1)."
    )
    print_line()

    print_exec("frame_1 = istpq(1, 3); print('1-Particle Frame:', frame_1)")
    print_exec("e1 = e(1); e_m1 = e(-1); e_m2 = e(-2); e_m3 = e(-3)")
    print_exec("print('e(1)^2  =', e1 * e1)")
    print_exec("print('e(-1)^2 =', e_m1 * e_m1)")

    print_line()
    print_fill(
        "We construct the 1-particle primitive idempotent P1 = 0.5 * (1 + e(1))."
    )
    print_exec("P1 = 0.5 * (1 + e1); print('P1 =', P1)")
    print_exec("print('P1 * P1 =', P1 * P1)")

    print_line()
    print_fill(
        "We define a spatial rotation rotor R in the (-1, -2) plane by angle theta = pi/3."
    )
    print_exec("theta = math.pi / 3")
    print_exec("B1_12 = e({-1, -2})")
    print_exec("R1 = math.cos(theta / 2) - B1_12 * math.sin(theta / 2)")
    print_exec("print('Rotor R1 =', R1)")

    print_line()
    print_fill(
        "State spinor transformation: Left product R1 * P1 vs Versor transformation P1 | R1."
    )
    print_exec("psi1 = R1 * P1; print('psi1 = R1 * P1 =', psi1)")
    print_exec("print('P1 | R1 =', P1 | R1)")

    pause()

    # =========================================================================
    # SECTION 2: 2-PARTICLE MSTA (n = 2, Signature: istpq(2, 6))
    # =========================================================================
    print_line()
    print_head("2. Two-Particle MSTA Space (n = 2, Frame: istpq(2, 6))")
    print_line()
    print_fill(
        "For 2 particles (2 qubits), the aggregate configuration space is 8D G_(2,6)."
    )
    print_fill(
        "Particle 1 uses generators { -3, -2, -1, 1 }. Particle 2 uses generators { -6, -5, -4, 2 }."
    )
    print_line()

    print_exec("frame_2 = istpq(2, 6); print('2-Particle Frame:', frame_2)")
    print_exec("P1 = 0.5 * (1 + e(1))")
    print_exec("P2 = 0.5 * (1 + e(2))")
    print_exec("P12 = P1 * P2; print('Joint Idempotent P12 =', P12)")
    print_exec("print('P12 * P12 =', P12 * P12)")

    print_line()
    print_fill(
        "We construct an inter-particle entangling rotor R_ent coupling Particle 1 and Particle 2."
    )
    print_fill(
        "Using the cross-particle bivector B_12 = e({-1, -4}) in plane (-1, -4)."
    )
    print_exec("B_cross = e({-1, -4})")
    print_exec("phi = math.pi / 4")
    print_exec("R_ent = math.cos(phi) - B_cross * math.sin(phi)")
    print_exec("print('Entangling Rotor R_ent =', R_ent)")

    print_line()
    print_fill("Applying the entangling rotor to prepare 2-particle entangled state psi12:")
    print_exec("psi12 = R_ent * P12; print('psi12 =', psi12)")

    print_line()
    print_fill("Grade decomposition of 2-particle state psi12:")
    print_exec("print('Grade 0 (Scalar):', psi12(0))")
    print_exec("print('Grade 2 (Bivector):', psi12(2))")
    print_exec("print('Grade 4 (4-Vector):', psi12(4))")

    pause()

    # =========================================================================
    # SECTION 3: 3-PARTICLE MSTA (n = 3, Signature: istpq(3, 9))
    # =========================================================================
    print_line()
    print_head("3. Three-Particle MSTA Space (n = 3, Frame: istpq(3, 9))")
    print_line()
    print_fill(
        "For 3 particles (3 qubits), the 12D aggregate Minkowski frame is G_(3,9)."
    )
    print_fill(
        "Particle 1: {-3..1}, Particle 2: {-6..2}, Particle 3: {-9..3}."
    )
    print_fill(
        "Multivector space dimension: 2^12 = 4096 real terms."
    )
    print_line()

    print_exec("frame_3 = istpq(3, 9); print('3-Particle Frame size:', frame_3.count())")
    print_exec("P3 = 0.5 * (1 + e(3))")
    print_exec("P123 = P1 * P2 * P3; print('3-Particle Idempotent P123 =', P123)")

    print_line()
    print_fill("Constructing 3-particle interaction rotor (Heisenberg spin chain coupling):")
    print_exec("B12 = e({-1, -4}); B23 = e({-4, -7})")
    print_exec("R_chain = (math.cos(math.pi/6) - B12 * math.sin(math.pi/6)) * (math.cos(math.pi/6) - B23 * math.sin(math.pi/6))")
    print_exec("psi123 = R_chain * P123; print('3-Particle State psi123 =', psi123)")

    print_line()
    print_fill("=========================================================================")
    print_fill("Joan Lasenby's 'Extra Room' Symmetries: Higher-Grade Inspection")
    print_fill("=========================================================================")
    print_fill("In QRA (Vašík), 3 qubits require 6 generators (G_(3,3), 64 terms).")
    print_fill("In MSTA (Lasenby), 3 particles occupy 12 generators (G_(3,9), 4096 terms).")
    print_fill("We extract higher-grade blades to isolate extra physical symmetries:")
    print_line()

    print_exec("print('Grade 0 Component  (Scalar):    ', psi123(0))")
    print_exec("print('Grade 2 Component  (Bivector):  ', psi123(2))")
    print_exec("print('Grade 4 Component  (4-Vector):  ', psi123(4))")
    print_exec("print('Grade 6 Component  (6-Vector):  ', psi123(6))")

    print_line()
    print_fill(
        "MSTA demonstration complete. MSTA successfully captures 1, 2, and 3-particle"
        + " quantum states and higher-grade geometric symmetries in PyClical!"
    )

    pause()
    print_line()
    print_fill(
        "You have completed the MSTA 1, 2, and 3-Particle PyClical demo (msta_demo.py)."
    )


if __name__ == "__main__":
    try:
        run()
    except (KeyboardInterrupt, Exception):  # pylint: disable=broad-exception-caught
        print("The demo was interrupted.")
