#!/usr/bin/env python

import sys, string, os, re

#-------------------------------------------------------------------------------

d_atom = {}

d_atom['He'] = ' 2'
d_atom['Ne'] = '10'

#-------------------------------------------------------------------------------

d_symtxt = {}

d_symtxt['auto'] = 'C   1              A'
d_symtxt['c1'  ] = 'C   1    0         A'
d_symtxt['c2'  ] = 'C   1    1XY       A'
d_symtxt['c2h' ] = 'C   1    2  ZXYZ   A'
d_symtxt['c2v' ] = 'C   1    2X  Y     A'
d_symtxt['ci'  ] = 'C   1    1XYZ      A'
d_symtxt['cs'  ] = 'C   1    1  Z      A'
d_symtxt['d2'  ] = 'C   1    2XY  YZ   A'
d_symtxt['d2h' ] = 'C   1    3  Z  Y  XA'

#-------------------------------------------------------------------------------

d_operator = {}

ex = ' XDIPLEN'
ey = ' YDIPLEN'
ez = ' ZDIPLEN'

bx = ''' 'b field'
 XAVECTOR
 ZDIPLEN
 YDIPLEN
 COMFACTOR
 68.517999904721'''
by = ''' 'b field'
 YAVECTOR
 XDIPLEN
 ZDIPLEN
 COMFACTOR
 68.517999904721'''
bz = ''' 'b field'
 ZAVECTOR
 YDIPLEN
 XDIPLEN
 COMFACTOR
 68.517999904721'''

d_operator['e'] = [ex, ey, ez]
d_operator['b'] = [bx, by, bz]

#-------------------------------------------------------------------------------

d_quantity = {}

d_quantity['e'] = ['N']
d_quantity['b'] = ['J', 'M', 'DIVJ', 'DIVM', 'ROTJ', 'ROTM']

#-------------------------------------------------------------------------------

d_iq = {}

d_iq['N'       ] = ' 1'
d_iq['M'       ] = ' 2'
d_iq['DIVM'    ] = ' 3'
d_iq['ROTM'    ] = ' 4'

d_iq['J'       ] = ' 6'
d_iq['DIVJ'    ] = ' 7'
d_iq['ROTJ'    ] = ' 8'
d_iq['JPARA'   ] = ' 9'
d_iq['DIVJPARA'] = '10'
d_iq['ROTJPARA'] = '11'
d_iq['JDIA'    ] = '12'
d_iq['VPARA'   ] = '13'
d_iq['DIVVPARA'] = '14'
d_iq['ROTVPARA'] = '15'
d_iq['ELF'     ] = '16'
d_iq['K'       ] = '17'
d_iq['H'       ] = '18'
d_iq['GAMMA5'  ] = '19'

#-------------------------------------------------------------------------------

d_orientation = {}

d_orientation[0] = 'x'
d_orientation[1] = 'y'
d_orientation[2] = 'z'

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
