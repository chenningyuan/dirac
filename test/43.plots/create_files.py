#!/usr/bin/env python

import sys, string, os, re
from common import *

#-------------------------------------------------------------------------------

os.system('rm -f *.inp *.mol *out* *PLOT DFCOEF* PAMXVC* *tgz*')

#-------------------------------------------------------------------------------

template_mol      = from_file('template_mol')
template_response = from_file('template_response')
template_plot     = from_file('template_plot')

for atom in d_atom:
    for sym in d_symtxt:

        mol = string.replace(template_mol, 'zz', d_atom[atom])
        mol = string.replace(mol, 'el', atom)
        mol = string.replace(mol, 'symtxt', d_symtxt[sym])

        mol_file = atom + '_' + sym + '.mol'
        to_file(mol, mol_file)

for field in d_operator:
    for orientation in d_orientation:

        response = string.replace(template_response, 'operator', d_operator[field][orientation])

        response_file = 'response_' + field + d_orientation[orientation] + '.inp'
        to_file(response, response_file)

        for quantity in d_quantity[field]:

            plot = string.replace(template_plot, 'quantity', d_iq[quantity])

            plot_file = 'plot_' + field + d_orientation[orientation] + '_' + quantity + '.inp'
            to_file(plot, plot_file)

#-------------------------------------------------------------------------------

run_script = '#!/bin/sh\n\n'

for atom in d_atom:
    for sym in d_symtxt:
        for field in d_operator:
            for orientation in d_orientation:

                mol_file = atom + '_' + sym
                inp_file = 'response_' + field + d_orientation[orientation]

                run_script = run_script + 'MOL=' + mol_file.lower() + '; export MOL\n'
                run_script = run_script + 'MENU=' + inp_file.lower() + '; export MENU\n'
                run_script = run_script + '../run.sh -get "PAMXVC" -outcmo\n\n'

                for quantity in d_quantity[field]:

                    inp_file = 'plot_' + field + d_orientation[orientation] + '_' + quantity

                    run_script = run_script + 'MENU=' + inp_file.lower() + '; export MENU\n'
                    run_script = run_script + '../run.sh -copy "PAMXVC" -get "PLOT*" -incmo; for file in PLOT*; do echo $file; mv $file $MENU\_$MOL.PLOT; done\n\n'

run_script = run_script + 'exit 0'

to_file(run_script, 'run.sh')

os.system('chmod 700 run.sh')
