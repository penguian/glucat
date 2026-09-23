#!/usr/bin/env python3
"""
Multiparticle Spacetime Algebra (MSTA) 1 to 5-Particle Simulation Demo.
Includes MSTA Hamiltonians, Energy Eigenstates, Unitary Time Evolution,
and Lasenby Higher-Grade Symmetries.
"""
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# msta_demo.py: This file demonstrates Multiparticle Spacetime Algebra (MSTA)
#               for 1, 2, 3, 4, and 5-particle quantum systems, showing frame setups,
#               idempotents, Hamiltonians, energy eigenstates, time evolution rotors,
#               and higher-grade symmetries up to 20D.
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

from pyclical import *
from pyclical_tutorial_utils import *


def run(ctx=tutorial_context(globals())):
    """
    Run the Multiparticle Spacetime Algebra (MSTA) 1 to 5-particle demo.
    """
    for name, method in get_object_methods(ctx).items():
        globals()[name] = method

    print_head("Multiparticle Spacetime Algebra (MSTA) Demo (1 to 5 Particles)")
    print_line()
    print_fill(
        "This demo illustrates Multiparticle Spacetime Algebra (MSTA) for 1, 2, 3, 4,"
        + " and 5 particles in PyClical / GluCat."
    )
    print_fill(
        "MSTA assigns a 4D Minkowski spacetime frame G_(1,3) to each particle,"
        + " representing quantum observables as multivectors, Hamiltonians as bivector"
        + " generators, and states as idempotent-projected spinors."
    )
    print_line()
    print_exec("from pyclical import *")
    print_exec("from pyclical_tutorial_utils import *")

    pause()

    # =========================================================================
    # SECTION 1: 1-PARTICLE MSTA (n = 1, Signature: istpq(1, 3))
    # =========================================================================
    print_head("1. Single-Particle MSTA Space (n = 1, Frame: istpq(1, 3))")
    print_line()
    print_fill(
        "For 1 particle (1 qubit), the Minkowski spacetime frame has 1 timelike generator"
        + " e(1) (e^2 = +1) and 3 spacelike generators e(-1), e(-2), e(-3) (e^2 = -1)."
    )
    print_line()

    print_exec("frame_1 = istpq(1, 3); print('1-Particle Frame:', frame_1)")
    print_exec("e1 = e(1); e_m1 = e(-1); e_m2 = e(-2); e_m3 = e(-3)")
    print_exec("P1 = 0.5 * (1 + e1); print('Primitive Idempotent P1 =', P1)")

    print_line()
    print_fill("--- Pauli Spin Observables & Zeeman Spin Hamiltonian ---")
    print_fill(
        "In MSTA, particle 1 spatial Pauli spin observables are represented by relative bivectors:"
    )
    print_fill(
        "  sigma_1x = e({1, -1}),  sigma_1y = e({1, -2}),  sigma_1z = e({1, -3})"
    )
    print_exec("sigma_1x = e({1, -1}); sigma_1y = e({1, -2}); sigma_1z = e({1, -3})")

    print_line()
    print_fill(
        "The 1-particle Zeeman spin Hamiltonian for magnetic field along z-axis is:"
    )
    print_fill("  H1 = 0.5 * w0 * sigma_1z")
    print_exec("w0 = 1.0")
    print_exec("H1 = 0.5 * w0 * sigma_1z; print('1-Particle Hamiltonian H1 =', H1)")

    print_line()
    print_fill("--- Energy Eigenstates ---")
    print_fill("Spin UP eigenstate  psi_up   = P1  (Energy E_up = +w0/2):")
    print_exec("psi_up = P1; print('psi_up =', psi_up)")
    print_fill("Spin DOWN eigenstate psi_down = e(-1) * P1 (Energy E_down = -w0/2):")
    print_exec("psi_down = e_m1 * P1; print('psi_down =', psi_down)")

    print_line()
    print_fill("--- Unitary Time Evolution Rotor ---")
    print_fill("Time evolution rotor R1(t) = exp(-H1 * t / 2):")
    print_exec("t = 1.5707963267948966; # t = pi / 2")
    print_exec("R1_t = exp(-H1 * t / 2.0); print('Time Evolution Rotor R1(t) =', R1_t)")
    print_exec("psi1_t = R1_t * psi_up; print('Evolved State psi1(t) =', psi1_t)")

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
    print_fill("Particle 1: { -3, -2, -1, 1 }. Particle 2: { -6, -5, -4, 2 }.")
    print_line()

    print_exec("frame_2 = istpq(2, 6); print('2-Particle Frame:', frame_2)")
    print_exec("P1 = 0.5 * (1 + e(1)); P2 = 0.5 * (1 + e(2))")
    print_exec("P12 = P1 * P2; print('Joint Idempotent P12 =', P12)")

    print_line()
    print_fill("--- Particle 2 Spin Observables & Heisenberg Exchange Hamiltonian ---")
    print_exec("sigma_2x = e({2, -4}); sigma_2y = e({2, -5}); sigma_2z = e({2, -6})")
    print_fill("The 2-particle Heisenberg exchange Hamiltonian is:")
    print_fill(
        "  H12 = J * (sigma_1x * sigma_2x + sigma_1y * sigma_2y + sigma_1z * sigma_2z)"
    )
    print_exec("J = 1.0")
    print_exec(
        "H12 = J * (sigma_1x * sigma_2x + sigma_1y * sigma_2y + sigma_1z * sigma_2z)"
    )
    print_exec("print('Heisenberg Hamiltonian H12 =', H12)")

    print_line()
    print_fill("--- 2-Particle Energy Eigenstates ---")
    print_fill(
        "Singlet Eigenstate  psi_singlet = (1/sqrt(2)) * (e(-2) - e(-1)) * P12  (Energy E_S = -3J):"
    )
    print_exec("psi_singlet = (1.0 / math.sqrt(2.0)) * (e(-2) - e(-1)) * P12")
    print_exec("print('psi_singlet =', psi_singlet)")
    print_fill(
        "Triplet m=0 Eigenstate psi_triplet = (1/sqrt(2)) * (e(-2) + e(-1)) * P12  (Energy E_T = +J):"
    )
    print_exec("psi_triplet = (1.0 / math.sqrt(2.0)) * (e(-2) + e(-1)) * P12")
    print_exec("print('psi_triplet =', psi_triplet)")

    print_line()
    print_fill("Grade decomposition of 2-particle singlet eigenstate:")
    print_exec("print('Grade 0 (Scalar):', psi_singlet(0))")
    print_exec("print('Grade 2 (Bivector):', psi_singlet(2))")
    print_exec("print('Grade 4 (4-Vector):', psi_singlet(4))")

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
    print_fill("Particle 1: {-3..1}, Particle 2: {-6..2}, Particle 3: {-9..3}.")
    print_fill("Multivector space dimension: 2^12 = 4096 real terms.")
    print_line()

    print_exec(
        "frame_3 = istpq(3, 9); print('3-Particle Frame size:', frame_3.count())"
    )
    print_exec("P3 = 0.5 * (1 + e(3))")
    print_exec("P123 = P12 * P3; print('3-Particle Idempotent P123 =', P123)")

    print_line()
    print_fill("--- 3-Particle Heisenberg Spin Chain Hamiltonian ---")
    print_exec("sigma_3x = e({3, -7}); sigma_3y = e({3, -8}); sigma_3z = e({3, -9})")
    print_exec(
        "H23 = J * (sigma_2x * sigma_3x + sigma_2y * sigma_3y + sigma_2z * sigma_3z)"
    )
    print_exec(
        "H31 = J * (sigma_3x * sigma_1x + sigma_3y * sigma_1y + sigma_3z * sigma_1z)"
    )
    print_exec(
        "H123 = H12 + H23 + H31; print('3-Particle Chain Hamiltonian H123 =', H123)"
    )

    print_line()
    print_fill("--- 3-Particle Unitary Time Evolution & GHZ State ---")
    print_exec(
        "R123_t = exp(-H123 * 0.5); print('Time Evolution Rotor R123(t) =', R123_t)"
    )
    print_exec(
        "psi123 = R123_t * P123; print('3-Particle Evolved State psi123 =', psi123)"
    )

    print_line()
    print_fill("Joan Lasenby's 'Extra Room' Symmetries: 3-Particle Grade Inspection")
    print_exec("print('Grade 0 Component  (Scalar):    ', psi123(0))")
    print_exec("print('Grade 2 Component  (Bivector):  ', psi123(2))")
    print_exec("print('Grade 4 Component  (4-Vector):  ', psi123(4))")
    print_exec("print('Grade 6 Component  (6-Vector):  ', psi123(6))")

    pause()

    # =========================================================================
    # SECTION 4: 4-PARTICLE MSTA (n = 4, Signature: istpq(4, 12))
    # =========================================================================
    print_line()
    print_head("4. Four-Particle MSTA Space (n = 4, Frame: istpq(4, 12))")
    print_line()
    print_fill(
        "For 4 particles (4 qubits), the 16D aggregate Minkowski frame is G_(4,12)."
    )
    print_fill(
        "Particle 1: {-3..1}, Particle 2: {-6..2}, Particle 3: {-9..3}, Particle 4: {-12..4}."
    )
    print_fill("Multivector space dimension: 2^16 = 65,536 real terms.")
    print_line()

    print_exec(
        "frame_4 = istpq(4, 12); print('4-Particle Frame size:', frame_4.count())"
    )
    print_exec("P4 = 0.5 * (1 + e(4))")
    print_exec("P1234 = P123 * P4; print('4-Particle Idempotent P1234 =', P1234)")

    print_line()
    print_fill("--- 4-Particle Ring Hamiltonian & Eigenstate Evolution ---")
    print_exec("sigma_4x = e({4, -10}); sigma_4y = e({4, -11}); sigma_4z = e({4, -12})")
    print_exec(
        "H34 = J * (sigma_3x * sigma_4x + sigma_3y * sigma_4y + sigma_3z * sigma_4z)"
    )
    print_exec(
        "H41 = J * (sigma_4x * sigma_1x + sigma_4y * sigma_1y + sigma_4z * sigma_1z)"
    )
    print_exec(
        "H1234 = H12 + H23 + H34 + H41; print('4-Particle Ring Hamiltonian H1234 =', H1234)"
    )

    print_exec(
        "R1234_t = exp(-H1234 * 0.25); print('4-Particle Time Evolution Rotor =', R1234_t)"
    )
    print_exec(
        "psi1234 = R1234_t * P1234; print('4-Particle State psi1234 =', psi1234)"
    )

    print_line()
    print_fill("4-Particle Grade Component Extraction:")
    print_exec("print('Grade 0 Component  (Scalar):    ', psi1234(0))")
    print_exec("print('Grade 2 Component  (Bivector):  ', psi1234(2))")
    print_exec("print('Grade 4 Component  (4-Vector):  ', psi1234(4))")
    print_exec("print('Grade 6 Component  (6-Vector):  ', psi1234(6))")
    print_exec("print('Grade 8 Component  (8-Vector):  ', psi1234(8))")

    pause()

    # =========================================================================
    # SECTION 5: 5-PARTICLE MSTA (n = 5, Signature: istpq(5, 15))
    # =========================================================================
    print_line()
    print_head("5. Five-Particle MSTA Space (n = 5, Frame: istpq(5, 15))")
    print_line()
    print_fill(
        "For 5 particles (5 qubits), the 20D aggregate Minkowski frame is G_(5,15)."
    )
    print_fill(
        "Particle 1: {-3..1}, Particle 2: {-6..2}, Particle 3: {-9..3}, Particle 4: {-12..4}, Particle 5: {-15..5}."
    )
    print_fill("Multivector space dimension: 2^20 = 1,048,576 real terms.")
    print_line()

    print_exec(
        "frame_5 = istpq(5, 15); print('5-Particle Frame size:', frame_5.count())"
    )
    print_exec("P5 = 0.5 * (1 + e(5))")
    print_exec("P12345 = P1234 * P5; print('5-Particle Idempotent P12345 =', P12345)")

    print_line()
    print_fill("--- 5-Particle Pentagonal Ring Hamiltonian & Time Evolution ---")
    print_exec("sigma_5x = e({5, -13}); sigma_5y = e({5, -14}); sigma_5z = e({5, -15})")
    print_exec(
        "H45 = J * (sigma_4x * sigma_5x + sigma_4y * sigma_5y + sigma_4z * sigma_5z)"
    )
    print_exec(
        "H51 = J * (sigma_5x * sigma_1x + sigma_5y * sigma_1y + sigma_5z * sigma_1z)"
    )
    print_exec(
        "H12345 = H12 + H23 + H34 + H45 + H51; print('5-Particle Ring Hamiltonian H12345 =', H12345)"
    )

    print_exec(
        "R12345_t = exp(-H12345 * 0.1); print('5-Particle Time Evolution Rotor =', R12345_t)"
    )
    print_exec(
        "psi12345 = R12345_t * P12345; print('5-Particle State psi12345 =', psi12345)"
    )

    print_line()
    print_fill(
        "5-Particle Grade Component Extraction (Higher-Grade Physical Correlations):"
    )
    print_exec("print('Grade 0 Component  (Scalar):    ', psi12345(0))")
    print_exec("print('Grade 2 Component  (Bivector):  ', psi12345(2))")
    print_exec("print('Grade 4 Component  (4-Vector):  ', psi12345(4))")
    print_exec("print('Grade 6 Component  (6-Vector):  ', psi12345(6))")
    print_exec("print('Grade 8 Component  (8-Vector):  ', psi12345(8))")
    print_exec("print('Grade 10 Component (10-Vector): ', psi12345(10))")

    print_line()
    print_fill(
        "MSTA demonstration complete. Hamiltonians, energy eigenstates, time evolution rotors,"
        + " and higher-grade symmetries are successfully formulated across 1 to 5 particles!"
    )

    pause()
    print_line()
    print_fill(
        "You have completed the MSTA 1 to 5-Particle PyClical demo (msta_demo.py)."
    )


if __name__ == "__main__":
    try:
        run()
    except (KeyboardInterrupt, Exception):  # pylint: disable=broad-exception-caught
        print("The demo was interrupted.")
