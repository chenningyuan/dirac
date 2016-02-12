#!/usr/bin/env python

# this script reads in two or three files:
# (1) "FROZEN"
#    rx  ry  rz  w  |elpot   nucpot  n   gn_x gn_y gn_z
#     1  2   3   4  | 5       6      7    8    9    10
# and two files with perturbed density, with the following order of columns:
# (2) "pertden_direct_lao.FINAL"
#      rx  ry  rz w | s_b(s=x,y,z; B=x,y,z)   gs_b(r=x,y,z; s=x,y,z; B=x,y,z)
#       1  2   3  4 |         5-13                14-40
# (3) "pertden_reorth_lao.FINAL"
#      rx  ry  rz w | n_b(B=x,y,z)  gn_b(r=x,y,z; B=x,y,z)  s_b(s=x,y,z; B=x,y,z)   gs_b(r=x,y,z; s=x,y,z; B=x,y,z)
#       1  2   3  4 |   5-7             8-16                       17-25                    26-52
#   pertden_direct_lao.FINAL and pertden_reorth_lao.FINAL files: only one may be present
#
# resulting file ('result') has:
#  rx ry rz w| all densities from FROZEN file | sum of (the same) perturbed densities from both pertden files*
#  1  2  3  4 |   5 - ...                     |     11-58
#                                               * in order: n_b, gn_b, s_b, gs_b

import os
import sys
from optparse import OptionParser
import string

usage = \
         """%prog --pertden_direct=* --pertden_reorth=* --frozen=* --result=*"""

parser = OptionParser(usage)

parser.add_option('--pertden_direct', 
                  type='string',
                  action='store',
                  dest='pertden_direct_filename',
                  help='file with direct-lao part of perturbed density',
                  metavar='FILE')

parser.add_option('--pertden_reorth', 
                  type='string',
                  action='store',
                  dest='pertden_reorth_filename',
                  help='file with reorth-lao part of perturbed density',
                  metavar='FILE')

parser.add_option('--frozen', 
                  type='string',
                  action='store',
                  dest='frozen_filename',
                  help='file with frozen density from fde',
                  metavar='FILE')

parser.add_option('--result', 
                  type='string',
                  action='store',
                  dest='result_filename',
                  help='final file',
                  metavar='FILE')

(options, args) = parser.parse_args()


def read_file(lines):
    g = []
    p = []
    q = []
    for i, line in enumerate(lines):
        l = line.strip().split()
        t = []
        for j in range(len(l)):
            t.append(l[j])
        q.append(t)
#...'transpoze'
    for i in range(len(lines[1].strip().split())):
        v = []
        for j in range(len(lines)):
            v.append(float(q[j][i].replace('D','E')))
        p.append(v)
#...first 3 elements are always grid points, 4th is the weight on a grid point
    g.append(p[0])
    g.append(p[1])
    g.append(p[2])
    g.append(p[3])
    del(p[0:4])
    return g,p


# read-in file with FROZEN density
# --------------------------------
f_frozen = open(options.frozen_filename, 'r')
lines_frozen = f_frozen.readlines()
lines_frozen = lines_frozen[1:]
f_frozen.close()
gridp_frozen, prop_frozen = read_file(lines_frozen)
nrp_frozen = len(lines_frozen)

#print 'nrp_frozen = ', nrp_frozen
#for i in range(len(gridp_frozen)):
#    for j in range(nrp_frozen):
#        print 'j, i, g = ', j, i, gridp_frozen[i][j]

# read-in two files with perturbed density
# ----------------------------------------
f_direct = open(options.pertden_direct_filename, 'r')
lines_direct = f_direct.readlines()
f_direct.close()
nrp_direct = len(lines_direct)
gridp_direct, prop_direct = read_file(lines_direct)

f_reorth = open(options.pertden_reorth_filename, 'r')
lines_reorth = f_reorth.readlines()
f_reorth.close()
nrp_reorth = len(lines_reorth)
gridp_reorth, prop_reorth = read_file(lines_reorth)

#print 'nrp_reorth = ', nrp_reorth
#for i in range(len(prop_reorth)):
#    for j in range(nrp_reorth):
#        print 'j, i, g = ', j, i, prop_reorth[i][j]

# before writing, check whether all the files have the same number of lines (grid points)
# and whether grid points are in the same order
# check
if (nrp_direct != nrp_reorth):
    print 'nrp_direct, nrp_reorth = ',  nrp_direct, nrp_reorth
    sys.exit("different number of points on pertden and reorth files")

if (nrp_frozen != nrp_direct):
    print 'nrp_frozen, nrp_direct = ',  nrp_frozen, nrp_direct
    sys.exit("different number of points on pertden and frozen files")

for i in range(nrp_direct):
    all_ok = ((gridp_direct[0]==gridp_reorth[0]) and (gridp_direct[1]==gridp_reorth[1]) and (gridp_direct[2]==gridp_reorth[2]) and (gridp_direct[3]==gridp_reorth[3]))
    if not all_ok:
        sys.exit("different grid points or points in different order on pertden_direct and pertden_reorth files")

for i in range(nrp_direct):
    all_ok = ((gridp_direct[0]==gridp_frozen[0]) and (gridp_direct[1]==gridp_frozen[1]) and (gridp_direct[2]==gridp_frozen[2]) and (gridp_direct[3]==gridp_frozen[3]))
    if not all_ok:
        sys.exit("different grid points or points in different order on pertden and frozen files")


# write result:
# grid points + frozen property + perturbed property
# --------------------------------------------------
print 'len(gridp_frozen)=', len(gridp_frozen)
print 'len(prop_frozen)=', len(prop_frozen)
print 'len(prop_direct)=',len(prop_direct)
g = open(options.result_filename, 'w')
# ncol we need is not the number of all columns, but ncol - 4 (minus 3 columns with grid points and 1 column with weight)
ncol = len(gridp_frozen) + len(prop_frozen) + len(prop_direct)
ncol = ncol - 4
g.write("%20d %20d \n" %(nrp_frozen, ncol))
for i in range(nrp_frozen):
    gridp  = ''.join("%20.12e" %(gridp_frozen[j][i]) for j in range(len(gridp_frozen)))
    frozen = ''.join("%20.12e" %(prop_frozen[j][i]) for j in range(len(prop_frozen)))
    prop   = ''.join("%20.12e" %(prop_direct[j][i]+prop_reorth[j][i]) for j in range(len(prop_direct)))
    g.write(str(gridp)+str(frozen)+str(prop)+'\n')
g.close()

