#!/usr/bin/env python
#
# utility to add cube files (radovan bast - 31 aug 2007)
#
# examples: ./cubicle.py -1 cube1 -f 0.5 -2 cube2 > new_cube
#           here new_cube = cube1 + 0.5*cube2
#
#           ./cubicle.py -1 cube1 -f "-1.0" -2 cube2 > new_cube
#           with a negative factor f it is possible to subtract cube files
#
# it assumes that the cube files have same dimensions, nr. of atoms, and so on

import sys, string, os, re
from optparse import OptionParser

#-------------------------------------------------------------------------------

parser = OptionParser()

parser.add_option('-f',
                  '--factor',
                  dest = 'factor',
                  help = 'factor')

parser.add_option('-1',
                  '--cube1',
                  dest = 'cube1',
                  help = 'filename of cube 1')

parser.add_option('-2',
                  '--cube2',
                  dest = 'cube2',
                  help = 'filename of cube 2')

(options, args) = parser.parse_args()

#-------------------------------------------------------------------------------

def to_file(s, file):

    f = open(file.lower(), 'w')
    f.write(s)
    f.close()

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

    return string.split(s2, ' '), s1

#-------------------------------------------------------------------------------

cube1 = from_file(options.cube1)
cube2 = from_file(options.cube2)

cube1_array, array = file_to_array(cube1)
cube2_array, array = file_to_array(cube2)

nr_atoms = int(cube1_array[1])

dim_x = int(cube1_array[ 5])
dim_y = int(cube1_array[ 9])
dim_z = int(cube1_array[13])

first = 17 + nr_atoms*5
last  = len(cube1_array) - 1

#-------------------------------------------------------------------------------

new_cube_array = []
f = float(options.factor)

for i in range(len(cube1_array)):

    if i < first:
       new_cube_array.append(cube1_array[i])
    else:
       if i < last + 1:
          sum = float(cube1_array[i]) + f*float(cube2_array[i])
          new_cube_array.append(sum)

#-------------------------------------------------------------------------------

for x in range(6 + nr_atoms):
    print array[x]

index = first
for i in range(dim_x):
    for j in range(dim_y):

        s = ''
        count = 0
        for k in range(dim_z):

            if count < 6:
               s = s + '  %8.5e' % new_cube_array[index]
               count = count + 1
            else:
               print s
               s = '  %8.5e' % new_cube_array[index]
               count = 1

            if k == dim_z - 1:
               print s

            index = index + 1
