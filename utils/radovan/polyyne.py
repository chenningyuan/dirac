#!/usr/bin/env python

# radovan (2009-11-23):
#
# create DIRAC molecule file for indealized polyyne
# with nr_unit_cells unit cells
# parameters from JCP 98, 8324
#
# example (8 unit cells):
# ./polyyne.py 8

import sys

#-------------------------------------------------------------------------------

nr_unit_cells = int(sys.argv[1])

c_basis = 'aug-cc-pVDZ'
h_basis = 'cc-pVDZ'

c_c        = 1.40
c_c_triple = 1.18
c_h        = 1.07

#-------------------------------------------------------------------------------

z_h = []
z_c = []

z = 0.0
z_h.append(z)
z += c_h
z_c.append(z)
z += c_c_triple
z_c.append(z)

for i in range(1, nr_unit_cells):
    z += c_c
    z_c.append(z)
    z += c_c_triple
    z_c.append(z)

z += c_h
z_h.append(z)

print 'INTGRL\n\n\nC%4i%15s' % (2, 'A')

print '       1.  %4i'   % 2
for z in z_h:
    print 'H    0.0  0.0  %f ' % z
print 'LARGE BASIS %s'   % h_basis

print '       6.  %4i'   % (2*nr_unit_cells)
for z in z_c:
    print 'C    0.0  0.0  %f ' % z
print 'LARGE BASIS %s'   % c_basis

print 'FINISH'
