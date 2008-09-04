#!/usr/bin/env python

import sys, string, os, re

#-------------------------------------------------------------------------------

d_atom = {}

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

d_hamiltonian = {}

d_hamiltonian['dc'       ] = 'LVCORR'
d_hamiltonian['dc_urkbal'] = 'LVCORR\n.URKBAL'
d_hamiltonian['dc_sof'   ] = 'LVCORR\n.SPINFREE'
#d_hamiltonian['iotc2'    ] = 'IOTC'
#d_hamiltonian['iotc4'    ] = 'IOTC4'
d_hamiltonian['nr2'      ] = 'NONREL'

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
