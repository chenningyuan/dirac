#!/usr/bin/python

#-------------------------------------------------------------------------------
#
# this script will calculate the vibrational levels (and wavefunctions)
# corresponding to a normal mode
# numericaly using the Numerov-Cooley algorithm
#
# last modification: 06/10/2007 (radovan bast)
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
#
# integration grid and optimization algorithm are extremely primitive
# but they do the job ...
#
#-------------------------------------------------------------------------------

import sys, string, math, Numeric, Ngl
from pylab import *
from numerov_parameters import *

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

x_l = []
y_l = []
for potential in potential_l:
    x_l.append(potential[0]/0.529177249)
    y_l.append(potential[1])

shift_minimum_to_zero(y_l)
potential_coefficients = polyfit(x_l, y_l, degree_potential)

x_l = []
y_l = []
for property in property_l:
    x_l.append(property[0]/0.529177249)
    y_l.append(property[1])

property_coefficients = polyfit(x_l, y_l, degree_property)

#-------------------------------------------------------------------------------

nr_points = nr_steps + 1
step      = (q_max - q_min)/nr_steps

q                      = zeros(nr_points, Float)
potential              = zeros(nr_points, Float)
property               = zeros(nr_points, Float)
g                      = zeros(nr_points, Float)
psi                    = zeros(nr_points, Float)
psi_normalized         = zeros(nr_points, Float)
psi_normalized_squared = zeros([nr_solutions, nr_points], Float)

for i in range(nr_points):
    q[i] = q_min + i*step
    for j in range(len(potential_coefficients)):
        k = len(potential_coefficients) - j - 1
        potential[i] = potential[i] + (q[i]**j)*potential_coefficients[k]
    for j in range(len(property_coefficients)):
        k = len(property_coefficients) - j - 1
        property[i] = property[i] + (q[i]**j)*property_coefficients[k]

#-------------------------------------------------------------------------------

energy_guess  = energy_precision
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

s = '%10s' % 'energy'               \
  + '%29s' % 'transition frequency' \
  + '%14s' % 'property'             \
  + '%20s' % 'shift'                \
  + '%21s' % 'shift'                \
  + '%16s' % 'nodes'
print s
s = '%10s' % '(a.u.)'  \
  + '%22s' % '(cm^-1)' \
  + '%20s' % '(a.u.)'  \
  + '%22s' % '(a.u.)'  \
  + '%20s' % '(Hz)'
print s

print

while nr_nodes_last < nr_solutions:

      nr_nodes = 0
      psi[0]   = 0
      psi[1]   = 1e-6

      for i in range(nr_points):
          g[i] = reduced_mass*1836.152701*(2*potential[i] - 2*energy_guess)
          if i > 1:
             psi[i] = 2*psi[i-1] - psi[i-2] + g[i-1]*psi[i-1]*step*step*5/6 + g[i-2]*psi[i-2]*step*step/12
             psi[i] = psi[i]/(1 - g[i]*step*step/12)
             if (psi[i]/psi[i-1]) < 0:
                nr_nodes = nr_nodes + 1
                i_save = i

      if nr_nodes > nr_nodes_last:
         if energy_step > 0:
            energy_step = energy_step/(-10)
      else:
         if energy_step < 0:
            energy_step = energy_step/(-10)

      if abs(energy_step) < energy_precision and nr_nodes > nr_nodes_last:

         for i in range(nr_points):
             if i > i_save:
                psi[i] = 0

         norm = 0
         for i in range(nr_points):
             norm = norm + psi[i]*psi[i]

         expectation_value = 0
         for i in range(nr_points):
             psi_normalized[i] = psi[i]/sqrt(norm)
             psi_normalized_squared[nr_nodes_last][i] = psi_normalized[i]**2
             expectation_value = expectation_value + property[i]*psi_normalized_squared[nr_nodes_last][i]

         if nr_nodes == 1:
            expectation_value_reference = expectation_value
            zero_point_energy = energy_guess

         transition_frequency = (energy_guess - zero_point_energy)*2.194746313705e5

         shift_au = (expectation_value - expectation_value_reference)
         shift_hz = shift_au*6.579683920721e15

         s = '%12.8f' % energy_guess               \
           + '%22.8f' % transition_frequency       \
           + '%22.8e' % expectation_value          \
           + '%22.8e' % shift_au                   \
           + '%22.8e' % shift_hz                   \
           + '%s'     % str(nr_nodes - 1).rjust(8)
         print s

         energy_step   = 1e-4
         nr_nodes_last = nr_nodes_last + 1

      energy_guess = energy_guess + energy_step

print
