Troubleshooting of the SCF convergence
=======================================


Obtaining converged molecular orbitals (or better, molecular spinors-MS)
is an essential step in each *ab-intio* quantum chemical calculation. We
provide you with some tricks on how to 'convince' the system to get
converged MS for a chosen state.

Basis set requirements
~~~~~~~~~~~~~~~~~~~~~~

First and important step: make sure that you have appropriate basis set
for your chemical system. Basis set should be large enough to describe
properly the electronic structure.

Spinors occupation
~~~~~~~~~~~~~~~~~~

Next parameters you can vary is the occupation of shells. If you are not
getting converged molecular spinors, for instance, you have high-spin
open-shell system, remove one or more electrons from your system.

Repeat the removing if necessary in order to obtain converged MS for the
charged system.

Use them as the starting set (-incmo flag form pam) for your system.
Sometimes this could help. If not, continue playing with the occupation
of spinors (for open-shells prefer the average occupation of several
degenerate spinors instead of one) and/or increase your basis set.

Excited states
~~~~~~~~~~~~~~

Generally it is very hard to get converged excited state using the very
basic single-determinant Hartree-Fock method. You might try it either by
dynamical overlap selection method (see SCF keywords) or directly by
placing electrons into lowest unoccupied spinors.

Oscillating SCF energies
~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes it happens that the SCF energy is oscillating between two
values for many iterations. A purely pragmatic solution to this problem
is to restart the calculation from the oscillating state and experiment
with either dynamic overlap selection (.OVLSEL) or the number of error
vectors used in the DIIS procedure (.MXDIIS).

Hf\ :sup:`2 +`\ , Ta\ :sup:`3 +`\ , W\ :sup:`4 +`\  with broken orbital degeneracies
------------------------------------------------------------------------------------

A typical symmetry breaking problem was encountered in Tel Aviv for
Lanthanide ions.

Hf\ :sup:`2 +`\ , Ta\ :sup:`3 +`\ , W\ :sup:`4 +`\  are converging
closed shell systems with the electronic structure of Yb (Z=70). The
first lanthanide in the line with the same electronic structure, Lu+, is
not converging.

One observes broken degeneracies of their atomic p,d,f... shells (using
four- and two-component Hamiltonians), even for spin-free and
Levy-Leblond Hamiltonian. However, the Yb closed shell atom, whose
electronic structure is emulated, comes out with proper degeneracies.

The reason for this problem is that the start guess for orbitals is
based on a bare nucleus or other simple model Hamiltonian. In the first
diagonalization this Hamiltonian gives orbitals that do not have the
expected order 4f, 6p but rather 6p, 4f. The program then proceeds by
occupying the three 6p orbitals, complemented by a remaining four 4f
orbitals. This density does not have the correct spherical symmetry and
offsets the SCF procedure yielding a symmetry-broken solution.

*Work-around*

One should provide a better starting guess, e.g. by first calculating
the highly-charged ion with a Xenon configuration. This gives a
symmetric density that can be used as starting guess (using the DFCOEF
file) for the real calculation. Some care is needed with this restart as
also the orbitals for the Xe-like ion may have an undesired order. This
may require use of the MO selection options of the program (see the
manual).

*Solution for Hf\ :sup:`2 +`\ *

The idea behind (though not exactly following *Work-around* is to fix
the wrong order of orbitals which is causing the degeneracy.

The 'default' wrong order of forming molecular orbitals is for s-shell
(LUMO) beeing inserted after the first d-shell (HOMO).

One has to apply the MJ-selection for the linear symmetry, or the
orbital reordering scheme to have s- separated from d-orbitals.

It is possible to get converged MOs not only for linear, but also for
D2h and C2v symmetries. However, for the C2v it takes many iterations.

