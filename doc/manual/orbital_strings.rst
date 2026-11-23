:orphan:
 

.. _orbital_strings:

Specification of orbital strings
================================

A subset of orbitals can be specified in two ways. Usually it suffices
to specify the keyword energy, followed by three real numbers to specify
a lower limit, upper limit and quasi-degeneracy threshold. All orbitals
that have an energy between the two thresholds are then selected. The
quasi-degeneracy threshold is used to extend the range if the lower or
upper threshold cuts through a set of closely energy-spaced orbitals. It
will move the upper or lower limit until it encounters a large enough
gap between orbital energies. The default value used in the AO-MO
integral transformation program to specify the active space for
correlated calculations is e.g. given by the string::

  energy -10.0 20.0 1.0

This picks out all orbitals that have an energy between -10.0 and +20.0
hartree (a.u.), with a minimum gap of 1.0 hartree.

For the newer implementation of orbitals strings used in ExaCorr it is also
possible to replace the upper or lower limit value by the word Fermi, for
instance ::

 energy -10.0 Fermi

to pick up all occupied orbitals above -10.0 Hartree. 

More detailed control is possible by using an orbital string. Here a
character string is given where individual orbitals are listed separated
by a comma or as a range of orbitals. The range limit is then separated
by "..". Positive energy and negative energy orbitals are indicated
using positive and negative numbers, respectively. For instance the
orbital string::

  -5..5,7

This picks out the five upper negative energy ("positronic") and five lower
positive energy ("electronic") orbitals plus positive energy orbital number
seven.

In the newer implementation (so ExaCorr only) one may substitute one or more
numbers by a string H-m or L+n to indicate HOMO-m or LUMO+n. For instance::

 H-2..L+5

to pick up all orbitals from HOMO-2 up to (including) LUMO+5.

If all orbitals (both positive and negative energy) are wanted, then one may
write::

  all

The concept of orbitals strings gives great flexibility in the specification of
orbitals and serves as an alternative to editing of vector files.

The later ExaCorr functionality is at the moment intended and tested for closed
shell cases. In this case the Fermi level is defined as (eps_HOMO+eps_LUMO)/2 where
eps indicates the orbital energy and HOMO is the highest (doubly) occupied orbital
and LUMO is the lowest unoccupied orbital. These definitions are retained also for
open shell cases, but be aware that in that case the concept of a Fermi level is less
clear and the selection may yield unexpected results. Explicit specification of energies
or orbital indices may then give more consistent results.
