#!/usr/bin/env python

import sys, string, os, re
from common import *

#-------------------------------------------------------------------------------

os.system('rm -f *.inp *.mol *out* DFCOEF')

#-------------------------------------------------------------------------------

template_mol = from_file('template_mol')
template_inp = from_file('template_inp')

for atom in d_atom:
    for sym in d_symtxt:

        mol = string.replace(template_mol, 'zz', d_atom[atom])
        mol = string.replace(mol, 'el', atom)
        mol = string.replace(mol, 'symtxt', d_symtxt[sym])

        mol_file = atom + '_' + sym + '.mol'
        to_file(mol, mol_file)

for hamiltonian in d_hamiltonian:

        inp = string.replace(template_inp, 'hamiltonian', d_hamiltonian[hamiltonian])

        inp_file = hamiltonian + '.inp'
        to_file(inp, inp_file)

#-------------------------------------------------------------------------------

run_script = '#!/bin/sh\n\n'

for atom in d_atom:
    for sym in d_symtxt:
        for hamiltonian in d_hamiltonian:

            mol_file = atom + '_' + sym
            inp_file = hamiltonian

            run_script = run_script + 'MOL=' + mol_file.lower() + '; export MOL\n'
            run_script = run_script + 'MENU=' + inp_file.lower() + '; export MENU\n'
            run_script = run_script + '../run.sh\n\n'

run_script = run_script + 'exit 0'

to_file(run_script, 'run.sh')

os.system('chmod 700 run.sh')
