:orphan:

.. _case_CmF:

===========================
Electronic structure of CmF
===========================

We want to investigate the structure of the diatomic molecule CmF (Curium Fluoride).
Following :cite:`Mochizuki_JCP2003` we assume that curium has valence configuration
in :math:`5f^77s^2` in the molecule, in particular since curium is a late actinide
and thus quite similar to lanthanides. The molecule will therefore have an atomic-like
open :math:`5f^7` shell. This may cause convergence problems since DIRAC by default
assumes open-shell orbitals to have energies above close-shell ones, which is likely to
not be the case here. In order to guide DIRAC we shall therefore employ the following
strategy:

1. We calculate the curium atom.
2. We import the :math:`5f^7` orbitals into the molecular calculation, place them in the open shell
   and keep them frozen in a first SCF calculation.
3. We then freeze all closed-shell orbitals, allowing the :math:`5f^7` orbitals to relax in the molecule.

The curium atom
---------------

The ground state configuration of curium (Cm, Z=96) is :math:`[Rn]5f^76d^17s^2` (see `here <https://www.webelements.com/curium/>`_ for general information 
and `also here <http://physics.nist.gov/PhysRefData/Handbook/Tables/curiumtable5.htm>`_), 
so with two open shells. We follow a two-step procedure to assure convergence:

1. We calculate the trivalent cation :math:`Cm^{3+}` with electron configuration :math:`[Rn]5f^7`.
2. Once we have identified the :math:`6d` and :math:`7s` orbitals amongst the virtuals we can easily calculate the neutral atom using reordering of orbitals and overlap selection.

For the molecular calcualtion we do not neccessarily need the fully converged :math:`5f` orbitals of the neutral atom, but we follow this strategy in part for educational purposes.
For ease of analysis we shall impose linear symmetry also in the atomic 
calculation. We achieve this through the introduction of a ghost center

.. literalinclude:: Cm.xyz

The menu file for the trivalent cation then reads

.. literalinclude:: CmIII.inp

Note that 

1. we perform Mulliken population on all positive-energy orbitals.
2. we use the keyword SHELL to only indicate atomic shell type since we do not need more detailed information
3. we specify that the ghost center has no basis.

We run the calculation using::

   pam --inp=CmIII --mol=Cm.xyz --put "cf.CmIII=DFCOEF"

(keywords for parallel run omitted). Here is a snippet of the Mulliken population analysis

.. literalinclude:: CmIII.pop

One clearly sees that the fractionally occupied :math:`5f` orbitals are followed by first virtual :math:`7s`, then virtual :math:`6d` orbitals.
We want to reorder them as :math:`\left(5f,7s,6d\right)\rightarrow\left(5f,6d,7s\right)` such that both open shells follow the closed :math:`7s` 
shell in the neutral atom. We therefore set up the input file

.. literalinclude:: Cm.inp

This assures that the orbitals are correctly ordered when we start the SCF iterations. In order to keep them in the desired order, we invoke 
overlap selection. We impose non-dynamic overlap selection, which means that orbitals in a given iteration are selected based on their overlap 
with the *initial* orbitals.

This calculation converges smoothly in 13 iterations.



