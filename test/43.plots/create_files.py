#!/usr/bin/env python

import sys, string, os, re

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

operator_d = {}

operator_d['1x'] = """ 'operator'
 XAVECTOR
 'NEF 03'
 'NEF 02'
 COMFACTOR
 -137.03599980944"""
operator_d['1y'] = """ 'operator'
 YAVECTOR
 'NEF 01'
 'NEF 03'
 COMFACTOR
 -137.03599980944"""
operator_d['1z'] = """ 'operator'
 ZAVECTOR
 'NEF 02'
 'NEF 01'
 COMFACTOR
 -137.03599980944"""

operator_d['2x'] = """ 'operator'
 XAVECTOR
 'NEF 06'
 'NEF 05'
 COMFACTOR
 -137.03599980944"""
operator_d['2y'] = """ 'operator'
 YAVECTOR
 'NEF 04'
 'NEF 06'
 COMFACTOR
 -137.03599980944"""
operator_d['2z'] = """ 'operator'
 ZAVECTOR
 'NEF 05'
 'NEF 04'
 COMFACTOR
 -137.03599980944"""

operator_d['ex'] = " 'XDIPLEN'"
operator_d['ey'] = " 'YDIPLEN'"
operator_d['ez'] = " 'ZDIPLEN'"

operator_d['bx'] = """ 'operator'
 XAVECTOR
 'ZDIPLEN'
 'YDIPLEN'
 COMFACTOR
 -68.517999904721"""
operator_d['by'] = """ 'operator'
 YAVECTOR
 'XDIPLEN'
 'ZDIPLEN'
 COMFACTOR
 -68.517999904721"""
operator_d['bz'] = """ 'operator'
 ZAVECTOR
 'YDIPLEN'
 'XDIPLEN'
 COMFACTOR
 -68.517999904721"""

#-------------------------------------------------------------------------------

template_response = from_file('template_response')

for operator in operator_d:

    s = string.replace(template_response, 'operator', operator_d[operator])

    file = 'response_' + operator + '.inp'
    to_file(s, file)

#-------------------------------------------------------------------------------




sys.exit()




#for hamiltonian in hamiltonian_d:
#    for kappa in ['', '_para', '_dia']:
#
#        endung = '_integration_12' + kappa
#        template = from_file('template' + endung)
#
#        s = string.replace(template, 'hamiltonian', hamiltonian_d[hamiltonian][0])
#        file = hamiltonian + endung + '.inp'
#        to_file(s, file)
#
#        endung = '_visual_12' + kappa
#        template = from_file('template' + endung)
#
#        s = string.replace(template, 'hamiltonian', hamiltonian_d[hamiltonian][0])
#        file = hamiltonian + endung + '.inp'
#        to_file(s, file)
#
#        endung = '_visual_21' + kappa
#        template = from_file('template' + endung)
#
#        s = string.replace(template, 'hamiltonian', hamiltonian_d[hamiltonian][0])
#        file = hamiltonian + endung + '.inp'
#        to_file(s, file)
#
#        endung = '_visual_j_1x' + kappa
#        template = from_file('template' + endung)
#
#        s = string.replace(template, 'hamiltonian', hamiltonian_d[hamiltonian][0])
#        file = hamiltonian + endung + '.inp'
#        to_file(s, file)
#
##-------------------------------------------------------------------------------
#
#s = '''#!/bin/bash\n\n'''
#
#for hamiltonian in hamiltonian_d:
#    for operator in operator_d:
#
#        s = s + 'cp %s.DFCOEF DFCOEF\n' % hamiltonian
#        text = hamiltonian + '_' + operator
#        s = s + '/home/bast/bin/dirac/snapshot/pam -mw 60 -incmo -outcmo -get PAMXVC co %s\n' % text
#        s = s + 'mv DFCOEF %s.DFCOEF\n' % hamiltonian
#        s = s + 'mv PAMXVC %s_%s.PAMXVC\n\n' % (hamiltonian, operator)
#
#to_file(s, 'response.run')
#os.system('chmod 700 response.run')
#
##-------------------------------------------------------------------------------
#
#s = '''#!/bin/bash\n\n'''
#
#for hamiltonian in hamiltonian_d:
#    s = s + 'cp %s.DFCOEF DFCOEF\n' % hamiltonian
#
#   #for operator in operator_d:
#    for operator in ['2x', '2y', '2z']:
#        s = s + 'cp %s_%s.PAMXVC %s____\n' % (hamiltonian, operator, operator)
#
#    for kappa in ['', '_para', '_dia']:
#        endung = '_integration_12' + kappa
#        s = s + '/home/bast/bin/dirac/snapshot/pam -mw 60 -incmo -copy '
#        s = s + '"2x____ 2y____ 2z____" co %s%s\n' % (hamiltonian, endung)
#
#    for operator in operator_d:
#        s = s + 'rm -f  %s____\n' % operator
#
#to_file(s, 'integration.run')
#os.system('chmod 700 integration.run')
#
##-------------------------------------------------------------------------------
#
#s = '''#!/bin/bash\n\n'''
#
#for hamiltonian in hamiltonian_d:
#    s = s + 'cp %s.DFCOEF DFCOEF\n' % hamiltonian
#
#    for operator in operator_d:
#        s = s + 'cp %s_%s.PAMXVC %s____\n' % (hamiltonian, operator, operator)
#
#    for kappa in ['', '_para', '_dia']:
#      # endung = '_visual_12' + kappa
#      # s = s + '/home/bast/bin/dirac/snapshot/pam -mw 60 -incmo -copy '
#      # s = s + '"2x____ 2y____ 2z____" -get PLOT_2D co %s%s\n' % (hamiltonian, endung)
#      # s = s + 'mv PLOT_2D %s%s_mol.PLOT_2D\n'                 % (hamiltonian, endung)
#      # endung = '_visual_21' + kappa
#      # s = s + '/home/bast/bin/dirac/snapshot/pam -mw 60 -incmo -copy '
#      # s = s + '"1x____ 1y____ 1z____" -get PLOT_2D co %s%s\n' % (hamiltonian, endung)
#      # s = s + 'mv PLOT_2D %s%s_mol.PLOT_2D\n'                 % (hamiltonian, endung)
#        endung = '_visual_j_1x' + kappa
#        s = s + '/home/bast/bin/dirac/snapshot/pam -mw 60 -incmo -copy '
#        s = s + '"1x____ 1y____ 1z____" -get PLOT_2D co %s%s\n' % (hamiltonian, endung)
#        s = s + 'mv PLOT_2D %s%s_mol.PLOT_2D\n'                 % (hamiltonian, endung)
#
#    for operator in operator_d:
#        s = s + 'rm -f  %s____\n' % operator
#
#to_file(s, 'visual.run')
#os.system('chmod 700 visual.run')
