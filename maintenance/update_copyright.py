#!/usr/bin/env python

"""
Script that updates copyright headers. (c) Radovan Bast
"""

import os
import re
import io


def convert_to_fortran_comments(s):
    s2 = []
    for line in s.split('\n'):
        if line == '/*':
            pass
        elif line == ' */':
            pass
        elif line == ' *':
            s2.append('!')
        elif line == '':
            pass
        else:
            s2.append('! ' + line[2:])
    return '\n' + '\n'.join(s2) + '\n'


copyright = {}

copyright[('dirac', 'c')] = '''
/*
 *     Copyright (c) 2019 by the authors of DIRAC.
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
 *     DIRAC program, see: http://www.diracprogram.org
 */
'''
copyright[('dirac', 'fortran')] = convert_to_fortran_comments(copyright[('dirac', 'c')])

copyright[('dalton', 'c')] = '''
/*
 *     Copyright (c) 2017 by the authors of Dalton (see below).
 *     All Rights Reserved.
 *
 *     The source code in this file is part of
 *     "Dalton, a molecular electronic structure program,
 *      Release DALTON2017 (2017), see http://daltonprogram.org"
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

copyright[('dalton', 'fortran')] = convert_to_fortran_comments(copyright[('dalton', 'c')])

files = {}
for lang in ['fortran', 'c']:
    files[lang] = []

for path, _, _files in os.walk(os.getcwd()):
    for f in _files:
        if f.endswith('.F90') or f.endswith('.f90') or f.endswith('.F'):
            files['fortran'].append(os.path.join(path, f))
        if f.endswith('.c'):
            files['c'].append(os.path.join(path, f))

for lang in ['fortran', 'c']:
    for code in ['dirac', 'dalton']:
        if lang == 'fortran':
            p = re.compile('(?<=!{}_copyright_start).*(?=!{}_copyright_end)'.format(code, code), re.DOTALL)
        else:
            p = re.compile('(?<=/\* {}_copyright_start \*/).*(?=/\* {}_copyright_end \*/)'.format(code, code), re.DOTALL)
        for file_name in files[lang]:
            with io.open(file_name, 'r', encoding="ascii", errors="ignore") as f:
                s = f.read()
            if '_copyright_start' in s:
                s = re.sub(p, copyright[(code, lang)], s)
                with io.open(file_name, 'w', encoding="ascii", errors="ignore") as f:
                    f.write(s)
