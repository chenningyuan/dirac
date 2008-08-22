#!/usr/bin/python

#-------------------------------------------------------------------------------
#
# this script will calculate the vibrational levels (and wavefunctions)
# corresponding to a normal mode
# numerically using the Numerov-Cooley algorithm
#
# last modification: 2008-08-21 (radovan.bast@chem.uit.no)
#
# you give:        potential at several points along the normal mode
#                  property values at several points along the normal mode
#                  reduced mass
#                  (edit numerov_parameters.py)
#
# parameters:      degree_potential = 3
#                    will fit a third order polynomial to the potential points
#
#                  degree_property  = 4
#                    will fit a fourth order polynomial to the property points
#
#                  nr_steps         = 500       means 501 grid points
#                  q_min            = -1.0      defines the range (in a.u.)
#                  q_max            = -q_min    defines the range (in a.u.)
#                  energy_precision = 1e-12
#                    script will accept eigenenergy if energy step below this threshold
#
#                  nr_solutions     = 2         you want two solutions (zero-point plus one)
#
#                  (edit numerov_parameters.py)
#
# script will do:  increase the energy and count the nodes of the wave function
#                  if the number of nodes changes and stepsize is below energy_precision
#                  it will accept the solution, integrate the property along q
#                  and move on to the next solution until nr_solutions is reached
#                  it will also calculate the transition frequencies 0 -> n
#                  this is useful to check against the harmonic frequencies
#
# advice:        * energy_precision is often more important than number of grid points
#                * be careful with the range (q_min and q_max)
#                * if the script enters an endless loop probably the reduced mass is wrong
#                  or the range (q_min and q_max) is too large
#                * practice first with the harmonic oscillator
#                * it is a good idea to play with parameters to check
#                  convergence and numerical stability
#                * potential and property are approximated by polynomials
#                  that contain coefficients FROM ZEROTH to nth order (that's what polyfit gives)
#                  and you might not want that
#                  (for instance you might insist that the gradient of the potential is zero at equilibrium)
#                  in this case you can provide your own
#                  expansion coefficients or program an alternative interpolation scheme
#
# integration grid and optimization algorithm are extremely primitive
# but they do the job ...
#
#-------------------------------------------------------------------------------

import sys
import string
import time

from pylab              import *
from Numeric            import *
from numerov_parameters import *

ts = time.clock()

#-------------------------------------------------------------------------------

def run_numerov(potential_coef,
                property_coef,
                reduced_mass,
                nr_steps,
                q_min,
                q_max,
                energy_precision,
                nr_solutions):

    s = '%10s' % 'energy'                   \
      + '%29s' % 'transition frequency'     \
      + '%14s' % 'property'                 \
      + '%20s' % 'shift'                    \
      + '%21s' % 'shift'                    \
      + '%16s' % 'nodes'                    \
      + '%23s' % '99.9% probability'
    print s

    s = '%10s' % '(a.u.)'  \
      + '%22s' % '(cm^-1)' \
      + '%20s' % '(a.u.)'  \
      + '%22s' % '(a.u.)'  \
      + '%20s' % '(Hz)'    \
      + '%33s' % 'within'
    print s

    n     = nr_steps + 1
    step  = (q_max - q_min)/nr_steps
    step2 = step**2.0

    q                      = zeros(n, Float)
    potential              = zeros(n, Float)
    property               = zeros(n, Float)
    g                      = zeros(n, Float)
    psi                    = zeros(n, Float)
    psi_squared_save       = zeros([nr_solutions, n], Float)
    energy_save            = zeros(nr_solutions, Float)

    for i in range(n):
        q[i]         = q_min + i*step
        potential[i] = polyval(potential_coef, q[i])
        property[i]  = polyval(property_coef,  q[i])

    energy_guess  = 1e-4
    energy_step   = 1e-4
    nr_nodes_last = 0

    expectation_value_reference = 0.0
    zero_point_energy           = 0.0

    while nr_nodes_last < nr_solutions:

          g = reduced_mass*1822.888479031408*2.0*(potential - energy_guess)

          psi[0]   = 0.0
          psi[1]   = 1.0e-6

          for i in range(n):
              if i > 1:
                 psi[i] = (2.0*psi[i-1] - psi[i-2] + g[i-1]*psi[i-1]*step2*5.0/6.0 + g[i-2]*psi[i-2]*step2/12.0)/(1.0 - g[i]*step2/12.0)

          nr_nodes = 0
          i_save   = n
          for i in range(n):
              if psi[i-1] != 0.0:
                 if (psi[i]/psi[i-1]) < 0.0:
                    nr_nodes = nr_nodes + 1
                    i_save = i

          psi[i_save+1:] = 0.0

          psi = psi/sqrt(dot(psi, psi))

          if (abs(energy_step) < energy_precision) and (nr_nodes > nr_nodes_last):

             norm = 0.0
             done = 0
             for i in range(len(psi)):
                 norm = norm + psi[i]*psi[i]
                 if (norm > 0.001) and not done:
                    i_left = i
                    done   = 1

             norm = 0.0
             done = 0
             for i in range(n-1, -1, -1):
                 norm = norm + psi[i]*psi[i]
                 if (norm > 0.001) and not done:
                    i_right = i
                    done    = 1

             expectation_value = dot(psi, property*psi)
             psi_squared_save[nr_nodes-1] = 4.0*psi**2.0

             energy_save[nr_nodes-1] = energy_guess

             if nr_nodes == 1:
                expectation_value_reference = expectation_value
                zero_point_energy = energy_guess

             transition_frequency = (energy_guess - zero_point_energy)*2.194746313705e5

             shift_au = (expectation_value - expectation_value_reference)
             shift_hz = shift_au*6.579683920721e15

             s = '%12.8f'          % energy_guess               \
               + '%22.8f'          % transition_frequency       \
               + '%22.8e'          % expectation_value          \
               + '%22.8e'          % shift_au                   \
               + '%22.8e'          % shift_hz                   \
               + '%s'              % str(nr_nodes - 1).rjust(8) \
               + '%15.2f .. %1.2f' % (q[i_left], q[i_right])
             print s

             energy_step   = 1.0e-4
             nr_nodes_last = nr_nodes_last + 1

          if nr_nodes > nr_nodes_last:
             if energy_step > 0.0:
                energy_step = energy_step/(-10.0)
          else:
             if energy_step < 0.0:
                energy_step = energy_step/(-10.0)

          energy_guess = energy_guess + energy_step

    return q, psi_squared_save, energy_save

#-------------------------------------------------------------------------------

potential_x   = array(potential_l)[:, 0]
potential_y   = array(potential_l)[:, 1]

potential_h_x = array(potential_h_l)[:, 0]
potential_h_y = array(potential_h_l)[:, 1]

property_x    = array(property_l)[:, 0]
property_y    = array(property_l)[:, 1]

bohr2angstrom = 0.529177249

potential_x   = potential_x/bohr2angstrom
potential_h_x = potential_h_x/bohr2angstrom
property_x    = property_x/bohr2angstrom

potential_y   -= min(potential_y)
potential_h_y -= min(potential_h_y)

potential_coef   = polyfit(potential_x,   potential_y,   degree_potential)
potential_h_coef = polyfit(potential_h_x, potential_h_y, 2)
property_coef    = polyfit(property_x,    property_y,    degree_property)

#-------------------------------------------------------------------------------

print
print 'reduced_mass     = ', reduced_mass
print 'degree_potential = ', degree_potential
print 'degree_property  = ', degree_property
print 'nr_steps         = ', nr_steps
print 'q_min            = ', q_min
print 'q_max            = ', q_max
print 'energy_precision = ', energy_precision
print 'nr_solutions     = ', nr_solutions
print

print 'anharmonic'
print '=========='
print

q, psi_squared_save, energy_save = run_numerov(potential_coef,
                                               property_coef,
                                               reduced_mass,
                                               nr_steps,
                                               q_min,
                                               q_max,
                                               energy_precision,
                                               nr_solutions)

#plot_numerov(potential_x_l,
#             potential_y_l,
#             potential_coef,
#             property_x_l,
#             property_y_l,
#             property_coef,
#             q,
#             psi_squared_save,
#             energy_save,
#             'q / a.u.',
#             'E(q) / E~B~h~N~',
#             'E~B~PV~N~(q) / E~B~h~N~',
#             'plot')

print
print 'harmonic'
print '========'
print

run_numerov(potential_h_coef,
            property_coef,
            reduced_mass,
            nr_steps,
            q_min,
            q_max,
            energy_precision,
            nr_solutions)

te = time.clock()
print
print 'time = %8.2f' % (te -ts)
