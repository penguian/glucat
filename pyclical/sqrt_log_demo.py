# -*- coding: utf-8 -*-
import sys
from PyCliCal import *
pi = 3.14159265358979

print '> e1=clifford("{1}")'
e1         =clifford("{1}")

print '> sqrt_e1=sqrt(e1)'
sqrt_e1         =sqrt(e1)
print '> print sqrt_e1'
print          sqrt_e1
print '> print sqrt_e1*sqrt_e1'
print          sqrt_e1*sqrt_e1
print '> print sqrt_e1*sqrt_e1-e1'
print          sqrt_e1*sqrt_e1-e1

print '> log_e1=log(e1)'
log_e1         =log(e1)
print '> print log_e1'
print          log_e1
print '> print exp(log_e1)'
print          exp(log_e1)
print '> print exp(log_e1)-e1'
print          exp(log_e1)-e1

print '> v=clifford("-2{1}+2{2}-3{3}")'
v         =clifford("-2{1}+2{2}-3{3}")

print '> sqrt_v=sqrt(v)'
sqrt_v         =sqrt(v)
print '> print sqrt_v'
print          sqrt_v
print '> print sqrt_v*sqrt_v'
print          sqrt_v*sqrt_v
print '> print sqrt_v*sqrt_v-v'
print          sqrt_v*sqrt_v-v

print '> log_v=log(v)'
log_v         =log(v)
print '> print log_v'
print          log_v
print '> print exp(log_v)'
print          exp(log_v)
print '> print exp(log_v)-v'
print          exp(log_v)-v

print '> x=clifford("-2{1}+2{2}-3{3}+4{-1,1,2}")'
x         =clifford("-2{1}+2{2}-3{3}+4{-1,1,2}")

print '> sqrt_x=sqrt(x)'
sqrt_x         =sqrt(x)
print '> print sqrt_x'
print          sqrt_x
print '> print sqrt_x*sqrt_x'
print          sqrt_x*sqrt_x
print '> print sqrt_x*sqrt_x-x'
print          sqrt_x*sqrt_x-x

print '> log_x=log(x)'
log_x         =log(x)
print '> print log_x'
print          log_x
print '> print exp(log_x)'
print          exp(log_x)
print '> print exp(log_x)-x'
print          exp(log_x)-x

print "> pq=lambda p,q:index_set(str(range(-q,p+1)).replace('[','{').replace(']','}'))"
pq=lambda p,q:index_set(str(range(-q,p+1)).replace('[','{').replace(']','}'))

print '> x=random_clifford(pq(2,1))'
x         =random_clifford(pq(2,1))
print '> print x'
print          x

print '> sqrt_x=sqrt(x)'
sqrt_x         =sqrt(x)
print '> print sqrt_x'
print          sqrt_x
print '> print sqrt_x*sqrt_x'
print          sqrt_x*sqrt_x
print '> print sqrt_x*sqrt_x-x'
print          sqrt_x*sqrt_x-x

print '> log_x=log(x)'
log_x         =log(x)
print '> print log_x'
print          log_x
print '> print exp(log_x)'
print          exp(log_x)
print '> print exp(log_x)-x'
print          exp(log_x)-x
