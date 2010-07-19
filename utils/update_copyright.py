#! /usr/bin/python

import string
import os
import re

#-------------------------------------------------------------------------------

dirac_copyright_fortran_start = '!dirac_copyright_start'
dirac_copyright_fortran_end   = '!dirac_copyright_end'

dirac_copyright_fortran = '''
!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.sdu.dk
'''

#-------------------------------------------------------------------------------

def from_file(file):

    f = open(file, 'r')
    s = f.read()
    f.close()
    return s

#-------------------------------------------------------------------------------

def to_file(s, file):

    f = open(file, 'w')
    f.write(s)
    f.close()

#-------------------------------------------------------------------------------

fortran_file_l = []
c_file_l       = []

pwd = os.getcwd()

for path, dir, files in os.walk(pwd):
    for file in files:
        if file[-4:] == '.F90':
           fortran_file_l.append(path + '/' + file)

# fortran files
p = re.compile('(?<=' + dirac_copyright_fortran_start + ').*(?=' + dirac_copyright_fortran_end + ')', re.DOTALL)
for file in fortran_file_l:
    s = from_file(file)
    if dirac_copyright_fortran_start in s:
       s = re.sub(p, dirac_copyright_fortran, s)
       to_file(s, file)
