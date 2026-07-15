#!/usr/bin/env python3
"""Time and relative directions in space demonstration script."""
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# clifford_demo.py: This is the demo file for the presentation
#                 Time and Relative Directions in Space.
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


def run(ctx):
    """Run the Time and Relative Directions in Space presentation demo."""
    for name, method in get_object_methods(ctx).items():
        globals()[name] = method

    print_fill("# Time and Relative Directions in Space.")

    pause()
    print_fill(
        "There is more to space and time than you might have learnt at school."
    )
    print_fill(
        "In this presentation, we explore space and time with the help of"
        + " Grassmann's and Clifford's geometric algebras."
    )
    print_fill(
        "This presentation is based on one given at Aston University in 2014,"
        + " hosted by John Fletcher."
    )
    print_exec("from PyClical import *")
    print_exec("from IPython.display import Image")

    pause()
    print_fill("A *vector* is a quantity with a *direction* and a *magnitude*.")
    print_line()
    print_fill("We know how to add vectors in space ...")
    print_exec("a = 3*e(1); print(a)")
    print_exec("b = 3*e(1)+4*e(2); print(b)")
    print_exec("c = a + b; print(c)")
    print_line()
    print_fill("... we know how to subtract vectors ...")
    print_exec("d = b - a; print(d)")
    print_line()
    print_fill(
        "... but how do we multiply vectors, divide vectors, rotate vectors !?"
    )

    pause()
    print_fill("There is more than one way to multiply vectors.")
    print_line()
    print_fill(
        r"The *scalar product* of $a$ and $b$ is denoted by $a \cdot b$"
        + r" (or $(a,b)$ or $\langle a,b \rangle$)."
    )
    print_fill(
        r"If we express $b$ as $b=b_{||}+b_{\perp}$, where $b_{||}$ is parallel"
        + r" to $a$, and $b_{\perp}$ is orthogonal to $a$, then we have"
    )
    print_fill(
        r"$|a \cdot b| = \| a \| \| b_{||} \|$, and $a \cdot b = b \cdot a$."
    )
    print_exec("print(a & b); print(b & a)")

    pause()
    print_fill(
        "Grassmann's *exterior (wedge) product* obeys different rules ..."
    )
    print_fill(
        r"$|a \wedge b| = \|a\| \|b_{\perp}\|$, and $a \wedge b = - b \wedge a$."
    )
    print_exec("print(a ^ b); print(b ^ a)")

    pause()
    print_exec(
        "Image(url="
        + "'http://upload.wikimedia.org/wikipedia/commons/6/6c/"
        + "Hermann_Gra%C3%9Fmann.jpg')"
    )
    print_fill("Hermann Grassmann (1809-1877) [Wikimedia Commons]")

    pause()
    print_fill(
        r"Given a basis $\{v_1, v_2, \ldots, v_n\} \subset \mathbb{R}^n$:"
    )
    print_fill(
        r"$v_1 \wedge v_2$ is a *bivector* - a directed area in a 2-plane."
    )
    print_fill(
        r"$v_1 \wedge v_2 \wedge v_3$ is a *trivector* - a directed volume"
        + r" of a 3-space."
    )
    print_fill(
        r"$v_1 \wedge v_2 \wedge \ldots \wedge v_n$ is an $n$-*vector* - a"
        + r" directed volume of an $n$-space."
    )
    print_line()
    print_fill(r"There is no $n+1$-vector in $n$-space. Why?")

    pause()
    print_fill("The wedge product detects linear dependence.")
    print_fill(
        r"The set of vectors $\{v_1, v_2, \ldots, v_m\}$ is linearly dependent"
        + r" if and only if the product $v_1 \wedge v_2 \wedge \ldots \wedge v_m = 0$."
    )
    print_exec("print(a)")
    print_exec("print(b)")
    print_exec("c = 3 * a - 2 * b; print(c)")
    print_exec("print(a ^ b)")
    print_exec("print(a ^ b ^ c)")

    pause()
    print_fill(
        r"Given a basis $\{v_1, v_2, \ldots, v_n\} \subset \mathbb{R}^n$,"
        + r" the *Grassmann algebra* of $\mathbb{R}^n$ is the $2^n$-dimensional"
        + r" space of *multivectors* of the form"
    )
    print_fill(
        r"$x_{\emptyset} + x_{\{1\}} v_1 + \ldots$"
        + r"$+ x_{\{n\}} v_n + x_{\{1,2\}} v_1 \wedge v_2 + \ldots + x_{\{n-1,n\}}$"
        + r" $v_{n-1} \wedge v_n + \ldots$"
        + r"$+ x_{\{1,\ldots,n\}} v_1 \wedge \ldots \wedge v_n$."
    )
    print_line()
    print_fill(
        "The wedge product and Grassmann algebras are used in differential"
        + " geometry and in physics and engineering."
    )

    pause()
    print_fill("William Kingdon Clifford had other ideas ...")
    print_exec(
        "Image(url="
        + "'http://upload.wikimedia.org/wikipedia/commons/e/eb/"
        + "Clifford_William_Kingdon.jpg')"
    )

    pause()
    print_fill(
        r"Clifford asked: what if $a a = (a,a) = \langle a,a \rangle = a \cdot a$"
        + r" for every vector $a \in \mathbb{R}^n$?"
    )
    print_fill(
        "Notice that if we expand $(a+b)(a+b)$ using the distributive law,"
        + " we get ..."
    )
    print_fill(
        r"$(a+b)(a+b) = a a + a b + b a + b b = a \cdot a + a b + b a + b \cdot b$."
    )
    print_fill(
        "On the other hand, by linearity of the scalar product,"
        + " we also get ..."
    )
    print_fill(
        r"$(a+b)(a+b) = (a+b) \cdot (a+b) = a \cdot a + 2 a \cdot b + b \cdot b$."
    )
    print_fill(r"So ... $a b + b a = 2 a \cdot b$.")
    print_fill(r"In particular, if $a \cdot b = 0$, then $a b = - b a$.")
    print_line()
    print_fill(
        "Vectors $a$ and $b$ commute if and only if they are parallel,"
        + " and anticommute if and only if they are orthogonal."
    )

    pause()
    print_fill(
        r"If $\{e_1, e_2, \ldots, e_n\}$ is an orthogonal basis of $\mathbb{R}^n$,"
        + " then the wedge and Clifford products of *disjoint* subsets of"
        + " this basis coincide ..."
    )
    print_fill(
        r"$e_1 e_2 = e_1 \wedge e_2, \ldots, e_1 e_2 \ldots e_n ="
        + r" e_1 \wedge e_2 \wedge \ldots \wedge e_n$"
    )
    print_exec("print(e(1)*e(2)); print(e(1) ^ e(2))")
    print_fill(
        r"... so that every multivector $x$ in the $2^n$-dimensional"
        + r" Grassmann algebra on $\mathbb{R}^n$ can be expressed as"
    )
    print_fill(
        r"$x_{\emptyset} + x_{\{1\}} e_1 + \ldots$"
        + r"$+ x_{\{n\}} e_n + x_{\{1,2\}} e_1 e_2 + \ldots + x_{\{n-1,n\}}$"
        + r" $e_{n-1} e_n + \ldots + x_{\{1,\ldots,n\}} e_1 e_2 \ldots e_n$."
    )
    print_fill(
        r"This $2^n$-dimensional vector space with the Clifford product is called"
        + r" the *Clifford algebra* on $\mathbb{R}^n$."
    )

    pause()
    print_fill("So ... what about *relative directions in space*?")
    print_fill(
        r"Well, now that we know that the Clifford square of every vector"
        + r" is the square of its length,"
        + r" we know how to calculate the Clifford inverse of every non-zero vector:"
        + r" $a^{-1} = \frac{a}{a a}$."
    )
    print_exec("print(a); print(inv(a)); print(inv(a)*a); print(a*inv(a))")

    pause()
    print_fill(
        r"With the convention that $a/b := a b^{-1}$, we can now divide vectors."
    )
    print_exec("print(a); print(a/b); print((a/b)*b)")

    pause()
    print_fill("Fine, but what use is that?")
    print_fill(
        "It turns out that the ability to divide vectors is very useful in expressing"
        + " orthogonal transformations."
    )
    print_fill(
        r"The expression $-a b / a = -a b a^{-1}$ is the reflection of $b$"
        + r" in the direction of $a$."
    )
    print_fill(
        r"$-a b a^{-1} = -a (b_{\perp} + b_{||}) a^{-1} = -a b_{\perp} a^{-1}$"
        + r" $-a b_{||} a^{-1} = a a^{-1} b_{\perp} - a a^{-1} b_{||}$"
        + r" $= b_{\perp} - b_{||}$."
    )
    print_exec("print(a)")
    print_exec("print(b); print(abs(b))")
    print_exec("c = b | a; print(c); print(abs(c))")

    pause()
    print_fill(
        r"The *Cartan-Dieudonne' Theorem* says that every orthogonal"
        + r" transformation on $\mathbb{R}^n$"
        + r" can be expressed as the product of at most $n$ vectors."
    )
    print_fill(
        r"In particular, $x \mapsto b a x (b a)^{-1}$ is a reflection in the"
        + r" direction of $a$ followed by a reflection in the direction of $b$."
    )
    print_fill(
        "This is also a rotation through *twice* the angle between $a$ and $b$."
    )
    print_exec("x = e(1)+2*e(2); print(x); print(x*x)")
    print_exec("y = x | (b * a); print(y); print(y*y)")
    print_exec("print(acos((a/abs(a)) & (b/abs(b))) * 180/pi)")
    print_exec("print(acos((x/abs(x)) & (y/abs(y))) * 180/pi)")

    pause()
    print_fill(
        "You can take the exponential of a multivector, using the usual Taylor"
        + " expansion ..."
    )
    print_fill(r"$\operatorname{e}^{x} = 1 + x + x^2/2 + x^3/(3!) + \ldots$")
    print_fill("The exponential of a vector is always a scalar plus a vector.")
    print_exec("print(a); print(exp(a))")
    print_exec("print(b); print(exp(b))")

    pause()
    print_fill(
        r"It turns out that every special orthogonal transformation can be expressed as"
        + r" the exponential of a bivector (when $n \geq 2$) ..."
    )
    print_fill(
        r"$x \mapsto \operatorname{e}^{B} x \operatorname{e}^{-B}"
        + r" = (-\operatorname{e}^{B}) x (-\operatorname{e}^{-B})$."
    )
    print_exec("print(x); print(x*x)")
    print_exec("y = exp(a ^ b) * x * exp(-a ^ b); print(y); print(y*y)")
    print_fill(
        r"... and the set of exponentials of bivectors forms a group called the"
        + r" *Spin group* $\operatorname{Spin}(n)$ over $\mathbb{R}^n$."
    )

    pause()
    print_fill("But whatever happened to *time* ...?")
    print_fill(
        r"Einstein and Minkowski showed that *space-time* can be described using"
        + r" the *Minkowski space* $\mathbb{R}^{3,1}$."
    )
    print_fill(
        r"This is a space with a basis $\{e_{\{-1\}},e_{\{1\}},e_{\{2\}},e_{\{3\}}\}$,"
        + r" where $e_{\{-1\}}^2 = -1$."
    )
    print_fill(
        "In this space, the special orthogonal transformations connected to the"
        + " identity are called the *restricted Lorentz transformations*."
    )
    print_fill(
        "Just like rotations in space, these can also be expressed as the"
        + " exponential of a bivector."
    )

    pause()
    print_exec("A = -3*e({-1,1})-e({-1,2})-2*e({-1,3}); print(A)")
    print_exec("s = exp(A/2); print(s)")
    print_exec("x = 2*e(-1)+2*e(1)+3*e(2)+5*e(3); print(x); print(x*x)")
    print_exec("y = s*x/s; print(y); print(y*y)")

    pause()
    print_fill(
        "There is even more to time and space than this, but *our* time has run out."
    )


if __name__ == "__main__":
    try:
        run(tutorial_context(globals()))
    except (KeyboardInterrupt, Exception):  # pylint: disable=broad-exception-caught
        print("The demo was interrupted.")
