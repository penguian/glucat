#!/usr/bin/env python
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

def run(ctx=tutorial_context(globals())):
    ctx_methods = get_object_methods(ctx)
    for method in ctx_methods:
        exec(method + " = ctx_methods['" + method +"']")

    print_fill("# Time and Relative Directions in Space.")

    pause()
    print_fill("There is more to space and time than you might have learnt at school.")
    print_fill("In this presentation, we explore space and time with the help of Grassmann's and Clifford's geometric algebras.")
    print_fill("This presentation is based on one given at Aston University in 2014, hosted by John Fletcher.")
    print_exec("from PyClical import *")
    print_exec("from IPython.display import Image")

    pause()
    print_fill("A *vector* is a quantity with a *direction* and a *magnitude*.")
    print_line()
    print_fill("We know how to add vectors in space ...")
    print_exec("a = 3*e(1); print a")
    print_exec("b = 3*e(1)+4*e(2); print b")
    print_exec("c = a + b; print c")
    print_line()
    print_fill("... we know how to subtract vectors ...")
    print_exec("d = b - a; print d")
    print_line()
    print_fill("... but how do we multiply vectors, divide vectors, rotate vectors !?")

    pause()
    print_fill("There is more than one way to multiply vectors.")
    print_line()
    print_fill("The *scalar product* of $a$ and $b$ is denoted by $a \cdot b$ (or $(a,b)$ or $\langle a,b \\rangle$).")
    print_fill("If we express $b$ as $b=b_{||}+b_{\perp}$, where $b_{||}$ is parallel to $a$, and"
              +" $b_{\perp}$ is orthogonal to $a$, then we have")
    print_fill("$|a \cdot b| = \| a \| \| b_{||} \|$, and $a \cdot b = b \cdot a$.")
    print_exec("print a & b; print b & a")

    pause()
    print_fill("Grassmann's *exterior (wedge) product* obeys different rules ...")
    print_fill("$|a \wedge b| = \|a\| \|b_{\perp}\|$, and $a \wedge b = - b \wedge a$.")
    print_exec("print a ^ b; print b ^ a")

    pause()
    print_exec("Image(url='http://upload.wikimedia.org/wikipedia/commons/6/6c/Hermann_Gra%C3%9Fmann.jpg')")
    print_fill("Hermann Grassmann (1809-1877) [Wikimedia Commons]")

    pause()
    print_fill("Given a basis $\\\{v_1, v_2, \ldots, v_n\\\} \subset \mathbb{R}^n$:")
    print_fill("$v_1 \wedge v_2$ is a *bivector* - a directed area in a 2-plane.")
    print_fill("$v_1 \wedge v_2 \wedge v_3$ is a *trivector* - a directed volume of a 3-space.")
    print_fill("$v_1 \wedge v_2 \wedge \ldots \wedge v_n$ is an $n$-*vector* - a directed volume of an $n$-space.")
    print_line()
    print_fill("There is no $n+1$-vector in $n$-space. Why?")

    pause()
    print_fill("The wedge product detects linear dependence.")
    print_fill("The set of vectors $\\\{v_1, v_2, \ldots, v_m\\\}$ is linearly dependent if and only if"
              +" the product $v_1 \wedge v_2 \wedge \ldots \wedge v_m = 0$.")
    print_exec("print a")
    print_exec("print b")
    print_exec("c = 3 * a - 2 * b; print c")
    print_exec("print a ^ b")
    print_exec("print a ^ b ^ c")

    pause()
    print_fill("Given a basis $\\\{v_1, v_2, \ldots, v_n\\\} \subset \mathbb{R}^n$,"
              +" the *Grassmann algebra* of $\mathbb{R}^n$ is the $2^n$-dimensional space of"
              +" *multivectors* of the form")
    print_fill("$x_{\\emptyset} + x_{\\\{1\\\}} v_1 + \ldots$"
              +"$+ x_{\\\{n\\\}} v_n + x_{\\\{1,2\\\}} v_1 \wedge v_2 + \ldots + x_{\\\{n-1,n\\\}} v_{n-1} \wedge v_n + \ldots$"
              +"$+ x_{\\\{1,\ldots,n\\\}} v_1 \wedge \ldots \wedge v_n$.")
    print_line()
    print_fill("The wedge product and Grassmann algebras are used in differential geometry and in physics and engineering.")

    pause()
    print_fill("William Kingdon Clifford had other ideas ...")
    print_exec("Image(url='http://upload.wikimedia.org/wikipedia/commons/e/eb/Clifford_William_Kingdon.jpg')")
    print_fill("William Kingdon Clifford (1845-1879) [Wikimedia Commons]")

    pause()
    print_fill("What if we *add* the inner product to the wedge product ... ?")
    print_exec("c = (a ^ b) + (a & b); print c")
    print_fill("This is called the *Clifford* (geometric) product of vectors.")
    print_fill("$a b := a \wedge b + a \cdot b$.")
    print_exec("c = a * b; print c")

    pause()
    print_fill("The Clifford product has remarkable properties.")
    print_fill("It detects when two vectors $v_1, v_2$ are orthogonal."
              +" In this case, they *anticommute*: $v_1 v_2 = - v_2 v_1$.")
    print_exec("v_1 = e(1) + 2 * e(2); v_2 = e(2) - 2*e(1)")
    print_exec("print v_1 * v_2; print v_2 * v_1; print v_1 & v_2")

    pause()
    print_fill("The Clifford square of a vector is the square of its length.")
    print_fill("$b b = b \wedge b + b \cdot b = b \cdot b = \|b\|^2$.")
    print_exec("b = 3*e(1)+4*e(2); print b; print b*b; print abs(b)**2")

    pause()
    print_fill("If $\\\{e_1, e_2, \ldots, e_n\\\}$ is an orthogonal basis of $\mathbb{R}^n$,"
              +" then the wedge and Clifford products of *disjoint* subsets of this basis coincide ...")
    print_fill("$e_1 e_2 = e_1 \wedge e_2, \ldots, e_1 e_2 \ldots e_n = e_1 \wedge e_2 \wedge \ldots \wedge e_n$")
    print_exec("print e(1)*e(2); print e(1) ^ e(2)")
    print_fill("... so that every multivector $x$ in the $2^n$-dimensional Grassmann algebra on $\mathbb{R}^n$ can be expressed as")
    print_fill("$x_{\\emptyset} + x_{\\\{1\\\}} e_1 + \ldots$"
              +"$+ x_{\\\{n\\\}} e_n + x_{\\\{1,2\\\}} e_1 e_2 + \ldots + x_{\\\{n-1,n\\\}} e_{n-1} e_n + \ldots$"
              +"$+ x_{\\\{1,\ldots,n\\\}} e_1 e_2 \ldots e_n$.")
    print_fill("This $2^n$-dimensional vector space with the Clifford product is called"
              +" the *Clifford algebra* on $\mathbb{R}^n$.")

    pause()
    print_fill("So ... what about *relative directions in space*?")
    print_fill("Well, now that we know that the Clifford square of every vector is the square of its length,"
              +" we know how to calculate the Clifford inverse of every non-zero vector: $a^{-1} = \\frac{a}{a a}$.")
    print_exec("print a; print inv(a); print inv(a)*a; print a*inv(a)")

    pause()
    print_fill("With the convention that $a/b := a b^{-1}$, we can now divide vectors.")
    print_exec("print a; print a/b; print (a/b)*b")

    pause()
    print_fill("Fine, but what use is that?")
    print_fill("It turns out that the ability to divide vectors is very useful in expressing orthogonal transformations.")
    print_fill("The expression $-a b / a = -a b a^{-1}$ is the reflection of $b$ in the direction of $a$.")
    print_fill("$-a b a^{-1} = -a (b_{\perp} + b_{||}) a^{-1} = -a b_{\perp} a^{-1} -a b_{||} a^{-1}$"
              +"$= a a^{-1} b_{\perp} - a a^{-1} b_{||} = b_{\perp} - b_{||}$.")
    print_exec("print a")
    print_exec("print b; print abs(b)")
    print_exec("c = -a * b / a; print c; print abs(c)")

    pause()
    print_fill("The *Cartan-Dieudonne' Theorem* says that every orthogonal transformation on $\mathbb{R}^n$"
              +" can be expressed as the product of at most $n$ vectors.")
    print_fill("In particular, $x \mapsto b a x (b a)^{-1}$ is a reflection in the direction of $a$"
              +" followed by a reflection in the direction of $b$.")
    print_fill("This is also a rotation throught *twice* the angle between $a$ and $b$.")
    print_exec("x = e(1)+2*e(2); print x; print x*x")
    print_exec("y = (b * a) * x * inv(b * a); print y; print y*y")
    print_exec("print acos((a/abs(a)) & (b/abs(b))) * 180/pi")
    print_exec("print acos((x/abs(x)) & (y/abs(y))) * 180/pi")

    pause()
    print_fill("You can take the exponential of a multivector, using the usual Taylor expansion ...")
    print_fill("$\operatorname{e}^{x} = 1 + x + x^2/2 + x^3/(3!) + \ldots$")
    print_fill("The exponential of a vector is always a scalar plus a vector.")
    print_exec("print a; print exp(a)")
    print_exec("print b; print exp(b)")

    pause()
    print_fill("It turns out that every special orthogonal transformation can be expressed as"
              +" the exponential of a bivector (when $n \geq 2$) ...")
    print_fill("$x \mapsto \operatorname{e}^{B} x \operatorname{e}^{-B} = (-\operatorname{e}^{B}) x (-\operatorname{e}^{-B})$.")
    print_exec("print x; print x*x")
    print_exec("y = exp(a ^ b) * x * exp(-a ^ b); print y; print y*y")
    print_fill("... and the set of exponentials of bivectors forms a group called the"
              +" *Spin group* $\operatorname{Spin}(n)$ over $\mathbb{R}^n$.")

    pause()
    print_fill("But whatever happened to *time* ...?")
    print_fill("Einstein and Minkowski showed that *space-time* can be described using"
              +" the *Minkowski space* $\mathbb{R}^{3,1}$.")
    print_fill("This is a space with a basis $\\\{e_{\\\{-1\\\}},e_{\\\{1\\\}},e_{\\\{2\\\}},e_{\\\{3\\\}}\\\}$,"
              +" where $e_{\\\{-1\\\}}^2 = -1$.")
    print_fill("In this space, the special orthogonal transformations connected to the identity"
              +" are called the *restricted Lorentz transformations*.")
    print_fill("Just like rotations in space, these can also be expressed as the exponential of a bivector.")

    pause()
    print_exec("A = -3*e({-1,1})-e({-1,2})-2*e({-1,3}); print A")
    print_exec("s = exp(A/2); print s")
    print_exec("x = 2*e(-1)+2*e(1)+3*e(2)+5*e(3); print x; print x*x")
    print_exec("y = s*x/s; print y; print y*y")

    pause()
    print_fill("There is even more to time and space than this, but *our* time has run out.")

if __name__ == "__main__":
    try:
        run()
    except:
        print("The demo was interrupted.")
        pass
