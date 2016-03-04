#!/usr/bin/env python

# this script needs 'numerical_grid' with dft grid,
# and creates 'numerical_grid_cont' in a format needed for fde
#

import os

f = open('numerical_grid', 'r')
lines = f.readlines()
f.close()

g = open('numerical_grid_cont', 'w')

nr_points = 0
for line in lines:
    list = line.strip().split()
    if (len(list) == 1):
        if (list[0] != '-1'):
            nr_points = nr_points + int(list[0])

g.write(str(nr_points)+'\n')
for line in lines:
    list = line.strip().split()
    if (len(list) > 1):
        g.write(line)

g.write('-1')

g.close()

# now rename, to avoid confusion:
os.rename('numerical_grid', 'numerical_grid.ORIG')
os.rename('numerical_grid_cont', 'numerical_grid')

