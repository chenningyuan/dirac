#!/usr/bin/env python
#
# utility to calculate various cube file derivatives (radovan bast - 31 aug 2007)
#
# you have a collection of cube files which plot something that depends on q
# in my case i have 11 cube files named cube+q and q ranges from -0.25 to +0.25
#
#                                       cube-0.25
#                                       cube-0.10
#                                       cube-0.05
#                                       cube-0.02
#                                       cube-0.01
#                                       cube+0.00
#                                       cube+0.01
#                                       cube+0.02
#                                       cube+0.05
#                                       cube+0.10
#                                       cube+0.25
#
# for another list q_l, edit q_l (below)
#
# example: ./cubicle_derivatives.py -p 3 -d 1
# will do: put a 3th order polynomial p(q) = a + b*q + c*q**2 + d*q**3
#          through the cubes and dump the first derivative to stdout (as cube)
#
# the script does not check for meanigless things like for instance: ./cubicle.py -p 1 -d 137
#
# WATCH OUT: the cube files should use the same grid (same points, not only same dimensions)!

from pylab import *
import sys, string, math, Numeric, re
from optparse import OptionParser

#-------------------------------------------------------------------------------

parser = OptionParser()

parser.add_option('-p',
                  '--polynomial_order',
                  dest = 'polynomial_order',
                  help = 'order of the polynomial')

parser.add_option('-d',
                  '--derivative',
                  dest = 'derivative',
                  help = 'derivative')

(options, args) = parser.parse_args()

#-------------------------------------------------------------------------------

def from_file(file):

    f = open(file, 'r')
    s = f.read()
    f.close()

    return s

#-------------------------------------------------------------------------------

def file_to_array(s1):

    s1 = string.split(s1, '\n')

    s2 = ''
    for i in range(len(s1)):
        if i > 1:
           s2 = s2 + string.replace(s1[i], '\n', ' ')

    for x in range(10):
        s2 = string.replace(s2, '  ', ' ')

    return string.split(s2, ' ')

#-------------------------------------------------------------------------------

def fact(x):
    return (1 if x==0 else x * fact(x-1))

#-------------------------------------------------------------------------------

def get_derivative(x_l, y_l, p_order, d_order):

    c_l = polyfit(x_l, y_l, p_order)

    d_l = []
    for i in range(d_order + 1):
        d_l.append(fact(i)*c_l[len(c_l) - i - 1])

    return d_l[d_order]

#-------------------------------------------------------------------------------

q_l = []

q_l.append(-0.25)
q_l.append(-0.10)
q_l.append(-0.05)
q_l.append(-0.02)
q_l.append(-0.01)
q_l.append( 0.00)
q_l.append( 0.01)
q_l.append( 0.02)
q_l.append( 0.05)
q_l.append( 0.10)
q_l.append( 0.25)

#-------------------------------------------------------------------------------

cube_l = []

for q in q_l:

    if q < 0:
       file = 'cube%4.2f' % q
    else:
       file = 'cube+%4.2f' % q

    cube_l.append(file_to_array(from_file(file)))

#-------------------------------------------------------------------------------

nr_atoms = int(cube_l[0][1])

dim_x = int(cube_l[0][ 5])
dim_y = int(cube_l[0][ 9])
dim_z = int(cube_l[0][13])

first = 17 + nr_atoms*5
last  = len(cube_l[0]) - 1

#-------------------------------------------------------------------------------

polynomial_order = int(options.polynomial_order)
derivative       = int(options.derivative)

new_cube = []

for i in range(len(cube_l[0])):

    if i < first:
       new_cube.append(cube_l[0][i])
    else:
       if i < last + 1:

          p_l = []
          for cube in cube_l:
              p_l.append(float(cube[i]))

          new_cube.append(get_derivative(q_l, p_l, polynomial_order, derivative))

#-------------------------------------------------------------------------------

header = string.split(from_file('cube+0.00'), '\n')
for x in range(6 + nr_atoms):
    print header[x]

index = first
for i in range(dim_x):
    for j in range(dim_y):

        s = ''
        count = 0
        for k in range(dim_z):

            if count < 6:
               s = s + '%13.5e' % new_cube[index]
               count = count + 1
            else:
               print s
               s = '%13.5e' % new_cube[index]
               count = 1

            if k == dim_z - 1:
               print s

            index = index + 1
