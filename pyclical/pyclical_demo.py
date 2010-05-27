# -*- coding: utf-8 -*-
import sys
from PyCliCal import *
pi = 3.14159265358979

def pause():
    if hasattr(sys, 'ps1') or hasattr(sys, 'ipcompleter'):
        if sys.stdin.isatty():
            raw_input("Press ENTER to continue")

print ""
print "This demonstration file contains examples on the following topics: "
print "  -Rotations in three dimensions  R^3"
print "  -Lorentz transformations of the Minkowski space-time  R^3,1 "
print "  -Electromagnetism in the Clifford algebras  Cl_3  and  Cl_3,1"
print "  -Complex Clifford algebra  CxCl_3,1"
print "  -Cayley algebra of octonions"
print "  -Conformal transformations and Vahlen matrices"
print "  -Elementary functions in the exterior algebra"
print "  -Pure spinors"
print "  -Conformal covariance of the Maxwell and Dirac equations"
pause()
print ""
print "ROTATIONS IN THREE DIMENSIONS  R^3"

print "Cl 3:"
e1 = clifford("{1}")
e2 = clifford("{2}")
e3 = clifford("{3}")

print "> e1 = ",e1
print "> e2 = ",e2
print "> e3 = ",e3

x          = 2*e1+3*e2+5*e3
print "> x = 2*e1+3*e2+5*e3"
print "   ==",x
a          = (3*e1+2*e2+e3)/10
print "> a = (3*e1+2*e2+e3)/10"
print "   ==",a
s          = exp(e1*e2*e3*a/2)
print "> s = exp(e1*e2*e3*a/2)"
print "   ==",s
print ""
pause()
print "Rotate the vector  x  about the axis  a  by the angle  ^a^"
y          = s*x/s
print "> y = s*x/s"
print "   ==",y
print ""
print "Check that the length is preserved by computing the squares of the vectors"
print      "> y*y"
print "   ==",y*y
print      "> x*x"
print "   ==",x*x
print ""
pause()
print "LORENTZ TRANSFORMATIONS IN THE MINKOWSKI SPACE-TIME  R^3,1"

print "Cl 3,1:"
e1 = clifford("{1}")
e2 = clifford("{2}")
e3 = clifford("{3}")
e4 = clifford("{-1}")

print "> e1 = ",e1
print "> e2 = ",e2
print "> e3 = ",e3
print "> e4 = ",e4

A          = 4*e1*e2+3*e1*e3+3*e1*e4+e2*e3+e2*e4+2*e3*e4
print "> A = 4*e1*e2+3*e1*e3+3*e1*e4+e2*e3+e2*e4+2*e3*e4"
print "   ==",A
s          = exp(A/2)
print "> s = exp(A/2)"
print "   ==",s
print ""
print "Lorentz transformation of a space-time point (or event)"
x          = 2*e1+3*e2+5*e3+2*e4
print "> x = 2*e1+3*e2+5*e3+2*e4"
print "   ==",x
y          = s*x/involute(s)
print "> y = s*x/involute(s)"
print "   ==",y
pause()
print "Check that the square-norm is preserved"
print      "> y*y"
print "   ==",y*y
print      "> x*x"
print "   ==",x*x
pause()
print "Lorentz transformation of an electromagnetic bivector  F"
F          = 5*e1*e2+e1*e3+2*e1*e4+7*e2*e3+e2*e4+2*e3*e4
print "> F = 5*e1*e2+e1*e3+2*e1*e4+7*e2*e3+e2*e4+2*e3*e4"
print "   ==",F
G          = s*F/s
print "> G = s*F/s"
print "   ==",G
print ""
print "Check that the Lorentz invariants are preserved (j = e1*e2*e3*e4)"
print      "> G*G/2"
print "   ==",G*G/2
print      "> F*F/2"
print "   ==",F*F/2
pause()
print "Compute the Poynting vector and the energy density"
print      "> -F*e4*F/2"
print "   ==",-F*e4*F/2
pause()
print "ELECTROMAGNETISM IN THE CLIFFORD ALGEBRAS  Cl_3  AND  Cl_3,1"
print ""
print "First, begin in the three-dimensional Euclidean space  R^3"
E          = e1+2*e2+4*e3
print "> E = e1+2*e2+4*e3"
print "   ==",E
B          = 3*e1+5*e2+7*e3
print "> B = 3*e1+5*e2+7*e3"
print "   ==",B
print ""
print "The role of the imaginary unit is played by the unit volume element "
i          = e1*e2*e3
print "> i = e1*e2*e3"
print "   ==",i
F          = E-i*B
print "> F = E-i*B"
print "   ==",F
pause()
print "Compute the Lorentz invariants"
print      "> (E*E-B*B)/2"
print "   ==",(E*E-B*B)/2
print      "> E&B"
print "   ==",E&B
print ""
print "The above computations can be combined in one formula"
print      "> F*F/2"
print "   ==",F*F/2
pause()
print "Compute the energy density and the Poynting vector"
print      "> (E*E+B*B)/2"
print "   ==",(E*E+B*B)/2
print      "> (E^B)/i"
print "   ==",(E^B)/i
print ""
print "The above computations can be combined in one formula"
print      "> -involute(F)*F/2"
print "   ==",-involute(F)*F/2
pause()
print "Consider a boost at half the velocity of light"
print "in the direction of the positive  x-axis"
v          = 0.5*e1
print "> v = 0.5*e1"
print "   ==",v
a          = atanh(v)
print "> a = atanh(v)"
print "   ==",a
s          = exp(a/2)
print "> s = exp(a/2)"
print "   ==",s
pause()
print "Lorentz transformation of the electromagnetic field  F"
G          = s*F/s
print "> G = s*F/s"
print "   ==",G
pause()
print "Extract out of  G  its magnetic induction"
print      "> i*G(2)"
print "   ==",i*G(2)
pause()
print "Lorentz transformation of a space-time point"
x          = 10+e1+e2
print "> x = 10+e1+e2"
print "   ==",x
y          = s*x/involute(s)
print "> y = s*x/involute(s)"
print "   ==",y
pause()
print "Check that the quadratic form is preserved"
print      "> y*conj(y)"
print "   ==",y*conj(y)
print      "> x*conj(x)"
print "   ==",x*conj(x)
pause()
print "Next, perform the same computations in the space-time  R^3,1"
pause()
E          = E*e4
print "> E = E*e4"
print "   ==",E
B          = B*e4
print "> B = B*e4"
print "   ==",B
print ""
print "The role of the imaginary unit is played by the unit volume element"
i          = e1*e2*e3*e4
print "> i = e1*e2*e3*e4"
print "   == ",i
F          = E-i*B
print "> F = E-i*B"
print "   ==",F
pause()
print "Compute the Lorentz invariants"
print      "> F*F/2"
print "   ==",F*F/2
pause()
print "Compute the Poynting vector and the energy density"
print      "> -F*e4*F/2"
print "   ==",-F*e4*F/2
pause()

print ""
print "COMPLEX CLIFFORD ALGEBRA  CxCl_3,1"

print "Cl 3,3:"
e1 = clifford("{1}")
e2 = clifford("{2}")
e3 = clifford("{3}")
e4 = clifford("{-3}")
e5 = clifford("{-2}")
e6 = clifford("{-1}")

print "> e1 = ",e1
print "> e2 = ",e2
print "> e3 = ",e3
print "> e4 = ",e4
print "> e5 = ",e5
print "> e6 = ",e6

print ""
print "The role of the imaginary unit is played by a unit bivector of  Cl_3,3"
i          = e5*e6
print "> i = e1*e2*e3*e4"
print "   == ",i
pause()
x          = 3*e1+2*e2+e4
print "> x = 3*e1+2*e2+e4"
print "   ==",x
y          = e3+5*e4
print "> y = e3+5*e4"
print "   ==",y
print ""
z          = x+i*y
print "> z = x+i*y"
print "   ==",z
pause()
A          = 2*e1*e2+2*e2*e4+3*e3*e4
print "> A = 2*e1*e2+2*e2*e4+3*e3*e4"
print "   ==",A
B          = e2*e3+2*e3*e1+3*e1*e4
print "> B = e2*e3+2*e3*e1+3*e1*e4"
print "   ==",B
print ""
C          = A+i*B
print "> C = A+i*B"
print "   ==",C
pause()
print "This computation might take a while:"
s          = exp(C/2)
print "> s = exp(C/2)"
print "   ==",s
pause()
print "Check that taking logarithm returns  C/2"
print      "> log(s)"
print "   ==",log(s)
pause()
print "Perform a complex rotation"
w          = s*z/s
print "> w = s*z/s"
print "   ==",w
pause()
print "Extract the real and imaginary parts"
u          = (w^i)/i
print "> u = (w^i)/i"
print "   ==",u
v          = (w-u)/i
print "> v = (w-u)/i"
print "   ==",v
pause()
print "Check that the invariants of the complex rotation are preserved"
print      "> x*x-y*y"
print "   ==",x*x-y*y
print      "> u*u-v*v"
print "   ==",u*u-v*v
print ""
print      "> x&y"
print "   ==",x&y
print      "> u&v"
print "   ==",u&v
pause()

print "CAYLEY ALGEBRA OF OCTONIONS"

print "Cl 0,7:"
e1 = clifford("{-7}")
e2 = clifford("{-6}")
e3 = clifford("{-5}")
e4 = clifford("{-4}")
e5 = clifford("{-3}")
e6 = clifford("{-2}")
e7 = clifford("{-1}")

print "> e1 = ",e1
print "> e2 = ",e2
print "> e3 = ",e3
print "> e4 = ",e4
print "> e5 = ",e5
print "> e6 = ",e6
print "> e7 = ",e7

print ""
print "Define a Cayley product of two paravectors in  R+R7 (= R+R^0,7)"
f          = (1-e1*e2*e4)*(1-e2*e3*e5)*(1-e3*e4*e6)*(1-e4*e5*e7)
print "> f = (1-e1*e2*e4)*(1-e2*e3*e5)*(1-e3*e4*e6)*(1-e4*e5*e7)"
print "   ==",f
o = lambda x,y:   real(x*y*f)+(x*y*f)(1)
print "> o(x,y) = real(x*y*f)+(x*y*f)(1)"
pause()
print      "> o(e1,2+3*e2)"
print "   ==",o(e1,2+3*e2)
pause()
a          = 1+2*e1+4*e2+2*e3
print "> a = 1+2*e1+4*e2+2*e3"
print "   ==",a
b          = 2+e1+2*e5
print "> b = 2+e1+2*e5"
print "   ==",b
pause()
c          = o(a,b)
print "> c = o(a,b)"
print "   ==",c
pause()
print "Check that the absolute value is preserved"
print      "> abs(c)"
print "   ==",abs(c)
print      "> abs(a)*abs(b)"
print "   ==",abs(a)*abs(b)
pause()
print "Consider the associator of three vectors in  R7"
a          = 3*e1+2*e4+e5
print "> a = 3*e1+2*e4+e5"
print "   ==",a
b          = 4*e2+3*e4+e5
print "> b = 4*e2+3*e4+e5"
print "   ==",b
c          = 3*e3+e4+2*e7
print "> c = 3*e3+e4+2*e7"
print "   ==",c
print ""
print "This computation might take a while:"
u          = o(o(a,b),c)-o(a,o(b,c))
print "> u = o(o(a,b),c)-o(a,o(b,c))"
print "   ==",u
pause()
print "Show that the associator is perpendicular to the factors"
V          = a^b^c
print "> V = a^b^c"
print "   ==",V
print      "> (u&V)/V"
print "   ==",(u&V)/V
pause()
print "Check the Moufang identity  (ab)*(ca) = a(bc)*a"
print      "> o(o(a,b),o(c,a))"
print "   ==",o(o(a,b),o(c,a))
print      "> o(o(a,o(b,c)),a)"
print "   ==",o(o(a,o(b,c)),a)
pause()

print "CONFORMAL TRANSFORMATIONS AND VAHLEN MATRICES"
print "in the Minkowski space-time  R^3,1"

print "Cl 4,2:"
e1 = clifford("{1}")
e2 = clifford("{2}")
e3 = clifford("{3}")
e4 = clifford("{4}")
e5 = clifford("{-2}")
e6 = clifford("{-1}")

print "> e1 = ",e1
print "> e2 = ",e2
print "> e3 = ",e3
print "> e4 = ",e4
print "> e5 = ",e5
print "> e6 = ",e6

print ""
print "We choose the two extra dimensions  e4, e6  so that our Minkowski space-time"
print "has a basis  e1, e2, e3, e5  where  e5*e5 = -1 "
epos          = e4
print "> epos = e4"
print "      == ",epos
eneg          = e6
print "> eneg = e6"
print "      == ",eneg
pause()
print "We give entries of a  2x2-matrix in the Lie algebra of the Vahlen group"
A          = 3+4*e1*e2+e2*e3+2*e3*e5
print "> A = 3+4*e1*e2+e2*e3+2*e3*e5"
print "   ==",A
B          = 7*e1+e2+3*e5
print "> B = 7*e1+e2+3*e5"
print "   ==",B
C          = 5*e2+e5
print "> C = 5*e2+e5"
print "   ==",C
D          = -reverse(A)
print "> D = -reverse(A)"
print "   ==",D
pause()
print "We give a matrix basis of the Lie algebra of the Vahlen group"
mat11          = (1+(epos^eneg))/2
print "> mat11 = (1+(epos^eneg))/2"
print "       ==",mat11
mat12          = (epos-eneg)/2
print "> mat12 = (epos-eneg)/2"
print "       ==",mat12
mat21          = (epos+eneg)/2
print "> mat21 = (epos+eneg)/2"
print "       ==",mat21
mat22          = (1-(epos^eneg))/2
print "> mat22 = (1-(epos^eneg))/2"
print "       ==",mat22
pause()
print "We form the following VAHLEN matrix"
print "in the Lie algebra of the Vahlen group"
VAHLEN          = A*mat11+B*mat12+involute(C)*mat21+involute(D)*mat22
print "> VAHLEN = A*mat11+B*mat12+involute(C)*mat21+involute(D)*mat22"
print "        ==",VAHLEN
pause()
print "By exponentiation we get the vahlen matrix"
vahlen          = exp(VAHLEN/2)
print "> vahlen = exp(VAHLEN/2)"
print "        ==",vahlen
pause()
Vah11          = (vahlen^epos^eneg)/(epos^eneg)
print "> Vah11 = (vahlen^epos^eneg)/(epos^eneg)"
print "       ==",Vah11
Vah12          = ((vahlen^eneg)/eneg-Vah11)/epos
print "> Vah12 = ((vahlen^eneg)/eneg-Vah11)/epos"
print "       ==",Vah12
Vah21          = ((vahlen^epos)/epos-Vah11)/eneg
print "> Vah21 = ((vahlen^epos)/epos-Vah11)/eneg"
print "       ==",Vah21
Vah22          = (vahlen-Vah11-Vah12*epos-Vah21*eneg)/(epos^eneg)
print "> Vah22 = (vahlen-Vah11-Vah12*epos-Vah21*eneg)/(epos^eneg)"
print "       ==",Vah22
pause()
print "Here are the entries ot the vahlen matrix"
a =          Vah11+Vah22
print "> a = Vah11+Vah22"
print "   ==",a
b =          Vah12-Vah21
print "> b = Vah12-Vah21"
print "   ==",b
c =          involute(Vah12+Vah21)
print "> c = involute(Vah12+Vah21)"
print "   ==",c
d =          involute(Vah11-Vah22)
print "> d = involute(Vah11-Vah22)"
print "   ==",d
pause()
print "Create a Mobius transformation"
g = lambda x:  (a*x+b)/(c*x+d)
print "> g(x)= (a*x+b)/(c*x+d)"
pause()
x          = 6*e1+e2+3*e3+4*e5
print "> x = 6*e1+e2+3*e3+4*e5"
print "   ==",x
y          = 3*e1+4*e2+e3+5*e5
print "> y = 3*e1+4*e2+e3+5*e5"
print "   ==",y
pause()
print      "> g(x)"
print "   ==",g(x)
print      "> g(y)"
print "   ==",g(y)
pause()
print "Check the formula  g(x)-g(y) = 1/reverse(c*x+d)*(x-y)/(c*y+d)"
print      "> g(x)-g(y)"
print "   ==",g(x)-g(y)
print      "> 1/reverse(c*x+d)*(x-y)/(c*y+d)"
print "   ==",1/reverse(c*x+d)*(x-y)/(c*y+d)
pause()

print "ELEMENTARY FUNCTIONS IN THE EXTERIOR ALGEBRA"
print "of the Minkowski space-time  R^3,1"
print ""

print "Cl 3,1:"
e1 = clifford("{1}")
e2 = clifford("{2}")
e3 = clifford("{3}")
e4 = clifford("{-1}")

print "> e1 = ",e1
print "> e2 = ",e2
print "> e3 = ",e3
print "> e4 = ",e4

x          = 3*e1+4*e2+2*e3+3*e4
print "> x = 3*e1+4*e2+2*e3+3*e4"
print "   ==",x
B          = pi/2*e1*e2+0.5*e3*e4
print "> B = pi/2*e1*e2+0.5*e3*e4"
print "   ==",B
pause()
print "Rotate  x  by the angle  pi/2  and perform a boost in the direction  e3"
u          = exp(B/2)
print "> u = exp(B/2)"
print "   ==",u
y          = u*x/u
print "> y = u*x/u"
print "   ==",y
pause()
print      "> y*y"
print "   ==",y*y
print      "> x*x"
print "   ==",x*x
pause()
print "Compare the above ordinary exponential  exp(A) = 1+A+AA/2+AAA/6+..."
print "to the following exterior exponential  expo(B) = 1+B+B^B/2+B^B^B/6:"
expo = lambda B: 1+B+(B^B)/2+(B^B^B)/6
v          = expo(B)
print "> v = expo(B)"
print "   ==",v
pause()
z          = v*x/v
print "> z = v*x/v"
print "   ==",z
print      "> z*z"
print "   ==",z*z
print      "> x*x"
print "   ==",x*x
pause()
print "Check that ordinary logarithm and exterior logarithm return  B/2  and  B"
print        "log(u)"
print "   ==",log(u)
pause()
logi = lambda x:  x-(x^x)/2+(x^x^x)/3
print "> logi(x)= x-(x^x)/2+(x^x^x)/3"
logo = lambda x:  log(real(x))+logi(x/real(x)-1)
print "> logo(x)= log(real(x))+logi(x/real(x)-1)"
pause()
print      "> logo(v)"
print "   ==",logo(v)
pause()
print      "> u*reverse(u)"
print "   ==",u*reverse(u)
print      "> v*reverse(v)"
print "   ==",v*reverse(v)
pause()

print "PURE SPINORS"

print "Cl 3,4:"
e1 = clifford("{1}")
e2 = clifford("{2}")
e3 = clifford("{3}")
e4 = clifford("{-4}")
e5 = clifford("{-3}")
e6 = clifford("{-2}")
e7 = clifford("{-1}")

print "> e1 = ",e1
print "> e2 = ",e2
print "> e3 = ",e3
print "> e4 = ",e4
print "> e5 = ",e5
print "> e6 = ",e6
print "> e7 = ",e7

v          = (e1+e4)*(e2+e5)*(e3+e6)/8
print "> v = (e1+e4)*(e2+e5)*(e3+e6)/8"
print "   ==",v
f          = (1+e1*e4)*(1+e2*e5)*(1+e3*e6)/8
print "> f = (1+e1*e4)*(1+e2*e5)*(1+e3*e6)/8"
print "   ==",f
print ""
print "Here  v = e1*e2*e3 f  as the following computation shows:"
print      "> e1*e2*e3*f"
print "   ==",e1*e2*e3*f
pause()
print "The following bivectors annul  v:"
print "e1*e4-e2*e5, e2*e5-e3*e6"
print "e1*e2+e1*e5, e1*e5+e2*e4, e1*e7+e4*e7"
print "e1*e3+e1*e6, e1*e6+e3*e4, e2*e7+e5*e7"
print "e2*e3+e2*e6, e2*e6+e3*e5, e3*e7+e6*e7"
print "and their exponential stabilizes  v."
pause()
print "As an example, take a linear combination of the above bivectors:"
B          = 4*(e1*e4-e2*e5)+7*(e1*e5+e2*e4)+3*(e3*e7+e6*e7)
print "> B = 4*(e1*e4-e2*e5)+7*(e1*e5+e2*e4)+3*(e3*e7+e6*e7)"
print "   ==",B
print ""
print "and exponentiate:        ... this computation might take a while ..."
s          = exp(B/10)
print "> s = exp(B/10)"
print "   ==",s
pause()
print "Check that  s v = v:"
print      "> s*v"
print "   ==",s*v
print "> v =",v
pause()
print "Take an arbitrary bivector:"
F          = 4*e1*e2+7*e2*e5+5*e3*e6+3*e4*e7+2*e5*e6
print "> F = 4*e1*e2+7*e2*e5+5*e3*e6+3*e4*e7+2*e5*e6"
print "   ==",F
print ""
print "and exponentiate:         ... this computation might take a while ..."
u          = exp(F/10)
print "> u = exp(F/10)"
print "   ==",u
pause()
print "The product  u*v  is a pure spinor"
w          = u*v
print "> w = u*v"
print "   ==",w
pause()
print "Check that  reverse(w)*Ak*w  vanishes for a pure spinor  w  and a  k-vector  Ak"
print "when  k = 0,1,2  but not when  k = 3  (v  and  reverse(w)*A3*w  are parallel)"
A1          = 2*e1+3*e2+5*e3+7*e4+e5+4*e6
print "> A1 = 2*e1+3*e2+5*e3+7*e4+e5+4*e6"
print "    ==",A1
A2          = 3*e1*e2+4*e1*e5+2*e2*e3+3*e2*e4+5*e3*e4+e3*e6+e4*e5
print "> A2 = 3*e1*e2+4*e1*e5+2*e2*e3+3*e2*e4+5*e3*e4+e3*e6+e4*e5"
print "    ==",A2
A3          = e1*e2*e3+4*e1*e2*e6+3*e2*e3*e4+2*e2*e4*e5+e3*e4*e5+5*e3*e5*e6+e4*e5*e6
print "> A3 = e1*e2*e3+4*e1*e2*e6+3*e2*e3*e4+2*e2*e4*e5+e3*e4*e5+5*e3*e5*e6+e4*e5*e6"
print "    ==",A3
pause()
print "These computations might take a while:"
print      "> reverse(w)*w"
print "   ==",reverse(w)*w
print      "> 10**18*reverse(w)*w"
print "   ==",1.0e18*reverse(w)*w
pause()
print      "> reverse(w)*A1*w"
print "   ==",reverse(w)*A1*w
print      "> 10**18*reverse(w)*A1*w"
print "   ==",1.0e18*reverse(w)*A1*w
pause()
print      "> reverse(w)*A2*w"
print "   ==",reverse(w)*A2*w
print      "> 10**18*reverse(w)*A2*w"
print "   ==",1.0e18*reverse(w)*A2*w
pause()
print      "> reverse(w)*A3*w"
print "   ==",reverse(w)*A3*w
print "> v =",v
pause()
print "You have completed the demonstration file pyclical_demo.py."
