# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# sqrt_log_demo.py: Demonstrate various sqrt and log calculations with PyClical.
#
#    copyright            : (C) 2008-2012 by Paul C. Leopardi
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

print '> e1=clifford("{1}")'
e1         =clifford("{1}")

print '> sqrt_e1=sqrt(e1)'
sqrt_e1         =sqrt(e1)
print '> print sqrt_e1'
print          sqrt_e1
print '> print sqrt_e1*sqrt_e1'
print          sqrt_e1*sqrt_e1
print '> print abs(sqrt_e1*sqrt_e1-e1)'
print          abs(sqrt_e1*sqrt_e1-e1)

print '> log_e1=log(e1)'
log_e1         =log(e1)
print '> print log_e1'
print          log_e1
print '> print exp(log_e1)'
print          exp(log_e1)
print '> print abs(exp(log_e1)-e1)'
print          abs(exp(log_e1)-e1)

print '> v=clifford("-2{1}+2{2}-3{3}")'
v         =clifford("-2{1}+2{2}-3{3}")

print '> sqrt_v=sqrt(v)'
sqrt_v         =sqrt(v)
print '> print sqrt_v'
print          sqrt_v
print '> print sqrt_v*sqrt_v'
print          sqrt_v*sqrt_v
print '> print abs(sqrt_v*sqrt_v-v)'
print          abs(sqrt_v*sqrt_v-v)

print '> log_v=log(v)'
log_v         =log(v)
print '> print log_v'
print          log_v
print '> print exp(log_v)'
print          exp(log_v)
print '> print abs(exp(log_v)-v)'
print          abs(exp(log_v)-v)

print '> x=clifford("-2{1}+2{2}-3{3}+4{-1,1,2}")'
x         =clifford("-2{1}+2{2}-3{3}+4{-1,1,2}")

print '> sqrt_x=sqrt(x)'
sqrt_x         =sqrt(x)
print '> print sqrt_x'
print          sqrt_x
print '> print sqrt_x*sqrt_x'
print          sqrt_x*sqrt_x
print '> print abs(sqrt_x*sqrt_x-x)'
print          abs(sqrt_x*sqrt_x-x)

print '> log_x=log(x)'
log_x         =log(x)
print '> print log_x'
print          log_x
print '> print exp(log_x)'
print          exp(log_x)
print '> print abs(exp(log_x)-x)'
print          abs(exp(log_x)-x)

print '> x=random_clifford(istpq(2,1))'
x         =random_clifford(istpq(2,1))
print '> print x'
print          x

print '> sqrt_x=sqrt(x)'
sqrt_x         =sqrt(x)
print '> print sqrt_x'
print          sqrt_x
print '> print sqrt_x*sqrt_x'
print          sqrt_x*sqrt_x
print '> print abs(sqrt_x*sqrt_x-x)'
print          abs(sqrt_x*sqrt_x-x)

print '> log_x=log(x)'
log_x         =log(x)
print '> print log_x'
print          log_x
print '> print exp(log_x)'
print          exp(log_x)
print '> print abs(exp(log_x)-x)'
print          abs(exp(log_x)-x)
