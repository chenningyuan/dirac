#!/usr/bin/python

#-------------------------------------------------------------------------------
#
# this script will calculate the vibrational levels (and wavefunctions)
# corresponding to a normal mode
# numericaly using the Numerov-Cooley algorithm
#
# last modification: 24/01/2008 (radovan bast)
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

import sys, string, math, Numeric, Ngl
from pylab import *
from numerov_parameters import *
#sys.path.append('/home/bast/schreibtisch/python')
#from gp import *

#-------------------------------------------------------------------------------

def shift_minimum_to_zero(l):

    minimum = l[0]
    for i in range(len(l)):
        if l[i] < minimum:
           minimum = l[i]

    for i in range(len(l)):
        l[i] = l[i] - minimum

    return l

#-------------------------------------------------------------------------------

potential_x_l = []
potential_y_l = []
for potential in potential_l:
    potential_x_l.append(potential[0]/0.529177249)
    potential_y_l.append(potential[1])

shift_minimum_to_zero(potential_y_l)
potential_coefficients = polyfit(potential_x_l, potential_y_l, degree_potential)

property_x_l = []
property_y_l = []
for property in property_l:
    property_x_l.append(property[0]/0.529177249)
    property_y_l.append(property[1])

property_coefficients = polyfit(property_x_l, property_y_l, degree_property)

#plot_one_x_two_y(potential_x_l, potential_y_l,
#                 property_x_l, property_y_l,
#                 'q / a.u.', 'E(q) / E~B~h~N~', 'E~B~PNC~N~(q) / E~B~h~N~', 'plot')

#-------------------------------------------------------------------------------

nr_points = nr_steps + 1
step      = (q_max - q_min)/nr_steps

q                      = zeros(nr_points, Float)
potential              = zeros(nr_points, Float)
property               = zeros(nr_points, Float)
g                      = zeros(nr_points, Float)
psi                    = zeros(nr_points, Float)

for i in range(nr_points):
    q[i] = q_min + i*step
    for j in range(len(potential_coefficients)):
        k = len(potential_coefficients) - j - 1
        potential[i] = potential[i] + (q[i]**j)*potential_coefficients[k]
    for j in range(len(property_coefficients)):
        k = len(property_coefficients) - j - 1
        property[i] = property[i] + (q[i]**j)*property_coefficients[k]

#-------------------------------------------------------------------------------

energy_guess  = 1e-4
energy_step   = 1e-4
nr_nodes_last = 0

expectation_value_reference = 0.0
zero_point_energy           = 0.0

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

print

#-------------------------------------------------------------------------------

def normalize(psi):

    norm = 0.0
    for i in range(len(psi)):
        norm = norm + psi[i]*psi[i]
    for i in range(len(psi)):
        psi[i] = psi[i]/sqrt(norm)

    return psi

#-------------------------------------------------------------------------------

while nr_nodes_last < nr_solutions:

      psi[0]   = 0.0
      psi[1]   = 1.0e-6
      for i in range(nr_points):
          g[i] = reduced_mass*1836.152701*(2.0*potential[i] - 2.0*energy_guess)
          if i > 1:
             psi[i] = 2.0*psi[i-1] - psi[i-2] + g[i-1]*psi[i-1]*step*step*5.0/6.0 + g[i-2]*psi[i-2]*step*step/12.0
             psi[i] = psi[i]/(1.0 - g[i]*step*step/12.0)

      psi = normalize(psi)

      nr_nodes = 0
      i_save   = nr_points
      for i in range(nr_points):
          if (psi[i]/psi[i-1]) < 0.0:
             nr_nodes = nr_nodes + 1
             i_save = i

      for i in range(nr_points):
          if i > i_save:
             psi[i] = 0.0

      psi = normalize(psi)

      if (abs(energy_step) < energy_precision) and (nr_nodes > nr_nodes_last):

         norm = 0.0
         done = 0
         for i in range(len(psi)):
             norm = norm + psi[i]*psi[i]
             if (norm > 0.999) and not done:
                i_right = i
                done    = 1

         norm = 0.0
         done = 0
         for i in range(nr_points-1, -1, -1):
             norm = norm + psi[i]*psi[i]
             if (norm > 0.999) and not done:
                i_left = i
                done   = 1

         expectation_value = 0.0
         for i in range(nr_points):
             expectation_value = expectation_value + property[i]*psi[i]*psi[i]

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

print
