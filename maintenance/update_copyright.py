#! /usr/bin/python

import string
import os
import re

#-------------------------------------------------------------------------------

dirac_copyright_fortran_start = '!dirac_copyright_start'
dirac_copyright_fortran_end   = '!dirac_copyright_end'

dirac_copyright_fortran = '''
!      Copyright (c) 2011 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.vu.nl
'''

dirac_copyright_c_start = '/\* dirac_copyright_start \*/'
dirac_copyright_c_end   = '/\* dirac_copyright_end \*/'

dirac_copyright_c = '''
/*
 *     Copyright (c) 2011 by the authors of DIRAC.
 *     All Rights Reserved.
 *
 *     This source code is part of the DIRAC program package.
 *     It is provided under a written license and may be used,
 *     copied, transmitted, or stored only in accordance to the
 *     conditions of that written license.
 *
 *     In particular, no part of the source code or compiled modules may
 *     be distributed outside the research group of the license holder.
 *     This means also that persons (e.g. post-docs) leaving the research
 *     group of the license holder may not take any part of Dirac,
 *     including modified files, with him/her, unless that person has
 *     obtained his/her own license.
 *
 *     For information on how to get a license, as well as the
 *     author list and the complete list of contributors to the
 *     DIRAC program, see: http://dirac.chem.vu.nl
 */
'''

dalton_copyright_fortran_start = '!dalton_copyright_start'
dalton_copyright_fortran_end   = '!dalton_copyright_end'

dalton_copyright_fortran = '''
!      Copyright (c) 2011 by the authors of Dalton (see below).
!      All Rights Reserved.
!   
!      The source code in this file is part of
!      "Dalton, a molecular electronic structure program,
!       Release DALTON2011 (2011), see http://daltonprogram.org"
!   
!      This source code is provided under a written licence and may be
!      used, copied, transmitted, or stored only in accord with that
!      written licence.
!   
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the licence holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the licence holder may not take any part of Dalton,
!      including modified files, with him/her, unless that person has
!      obtained his/her own licence.
!   
!      For further information, including how to get a licence, see:
!         http://daltonprogram.org
'''

dalton_copyright_c_start = '/\* dalton_copyright_start \*/'
dalton_copyright_c_end   = '/\* dalton_copyright_end \*/'

dalton_copyright_c = '''
/*
 *     Copyright (c) 2011 by the authors of Dalton (see below).
 *     All Rights Reserved.
 *  
 *     The source code in this file is part of
 *     "Dalton, a molecular electronic structure program,
 *      Release DALTON2011 (2011), see http://daltonprogram.org"
 *  
 *     This source code is provided under a written licence and may be
 *     used, copied, transmitted, or stored only in accord with that
 *     written licence.
 *  
 *     In particular, no part of the source code or compiled modules may
 *     be distributed outside the research group of the licence holder.
 *     This means also that persons (e.g. post-docs) leaving the research
 *     group of the licence holder may not take any part of Dalton,
 *     including modified files, with him/her, unless that person has
 *     obtained his/her own licence.
 *  
 *     For further information, including how to get a licence, see:
 *        http://daltonprogram.org
 */
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
        if file[-4:] == '.f90':
           fortran_file_l.append(path + '/' + file)
        if file[-2:] == '.F':
           fortran_file_l.append(path + '/' + file)
        if file[-2:] == '.c':
           c_file_l.append(path + '/' + file)

# fortran files dirac copyright
p = re.compile('(?<=' + dirac_copyright_fortran_start + ').*(?=' + dirac_copyright_fortran_end + ')', re.DOTALL)
for file in fortran_file_l:
    s = from_file(file)
    if dirac_copyright_fortran_start in s:
       s = re.sub(p, dirac_copyright_fortran, s)
       to_file(s, file)

# fortran files dalton copyright
p = re.compile('(?<=' + dalton_copyright_fortran_start + ').*(?=' + dalton_copyright_fortran_end + ')', re.DOTALL)
for file in fortran_file_l:
    s = from_file(file)
    if dalton_copyright_fortran_start in s:
       s = re.sub(p, dalton_copyright_fortran, s)
       to_file(s, file)

# c files dirac copyright
p = re.compile('(?<=' + dirac_copyright_c_start + ').*(?=' + dirac_copyright_c_end + ')', re.DOTALL)
for file in c_file_l:
    s = from_file(file)
    if string.replace(dirac_copyright_c_start, '\\', '') in s:
       s = re.sub(p, dirac_copyright_c, s)
       to_file(s, file)

# c files dalton copyright
p = re.compile('(?<=' + dalton_copyright_c_start + ').*(?=' + dalton_copyright_c_end + ')', re.DOTALL)
for file in c_file_l:
    s = from_file(file)
    if string.replace(dalton_copyright_c_start, '\\', '') in s:
       s = re.sub(p, dalton_copyright_c, s)
       to_file(s, file)
