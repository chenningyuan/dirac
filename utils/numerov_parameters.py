#!/usr/bin/python

reduced_mass     = 9.7031
degree_potential = 8
degree_property  = 6
nr_steps         = 3000
q_min            = -1.0
q_max            = -q_min
energy_precision = 1e-10
nr_solutions     = 3

# units: reduced mass is in amu, as it comes out in GAUSSIAN
#        q is in angstrom (makes sense if you use GAUSSIAN)
#        potential and property are in hartree

# btw: this example is the HF PV shift (C-F stretch)
#      in CHFClBr (table II in Schwerdtfeger et al. PRA 71, 012103 (2005))

potential_l = []
#                     q     potential(q)
potential_l.append([-0.25,-0.61093334302e+03])
potential_l.append([-0.10,-0.61102237795e+03])
potential_l.append([-0.05,-0.61103061231e+03])
potential_l.append([-0.02,-0.61103259217e+03])
potential_l.append([-0.01,-0.61103285334e+03])
potential_l.append([ 0.00,-0.61103293760e+03])
potential_l.append([ 0.01,-0.61103285720e+03])
potential_l.append([ 0.02,-0.61103262365e+03])
potential_l.append([ 0.05,-0.61103110894e+03])
potential_l.append([ 0.10,-0.61102641160e+03])
potential_l.append([ 0.25,-0.61100326033e+03])

potential_harmonic_l = []
#                              q     potential(q)
potential_harmonic_l.append([-0.02,-0.61103259217e+03])
potential_harmonic_l.append([-0.01,-0.61103285334e+03])
potential_harmonic_l.append([ 0.00,-0.61103293760e+03])
potential_harmonic_l.append([ 0.01,-0.61103285720e+03])
potential_harmonic_l.append([ 0.02,-0.61103262365e+03])

property_l = []
#                    q     property(q)
property_l.append([-0.25, 1.525781348757E-17])
property_l.append([-0.10, 8.728600913297E-18])
property_l.append([-0.05, 7.002011374476E-18])
property_l.append([-0.02, 6.088383502529E-18])
property_l.append([-0.01, 5.804338628868E-18])
property_l.append([ 0.00, 5.529990630427E-18])
property_l.append([ 0.01, 5.266120321598E-18])
property_l.append([ 0.02, 5.012527317735E-18])
property_l.append([ 0.05, 4.309240003629E-18])
property_l.append([ 0.10, 3.325129801526E-18])
property_l.append([ 0.25, 1.610579307504E-18])
