

Polarizable Continuum Model
=======================================================

The possibility to perform solvent calculations is currently under development in the 
``pcm-response`` branch. We will exploit a *Continuum Solvation Model* (CSM) namely the
*Polarizable Continuum Model* (PCM).
PCM is a *focused model*: the solute (a single molecule or a cluster containing the solute 
and some relevant solvent molecules) is described quantum mechanically, while the solvent is 
approximated as a structureless continuum whose interaction with the solute is mediated by 
its permittivity, :math:`\varepsilon`.
  
The solute is accomodated inside a molecular cavity, built as a set of interlocking spheres 
centered on the atoms constituting the molecule under investigation.
Although in principle we could exploit whatever level of QM theory, we will limit ourselves 
to the implementation of an SCF description of the solute.

Basic Theory
------------------------

In CSMs we write the Schrödinger equation for the solute as:

.. math::

   \left[ \mathcal{H}_0 + \mathcal{V}_{\mathrm{int}}(\rho_{\mathrm{M}})\right] | \psi \rangle = E | \phi \rangle

where :math:`\mathcal{H}_0` is the solute Hamiltonian *in vacuo* and :math:`\mathcal{V}_\mathrm{int}(\rho_\mathrm{M})`
is the solute-solvent interaction potential, **which depends on the solute wavefunction** through the first-order
density matrix. This introduces a *nonlinearity* which must be dealt with appropriately.
In standard molecular electronic-structure theory we write the energy as the expectation value of the Hamiltonian
(:math: `| \phi \rangle` being a trial wave function):

.. math::
   \langle \phi | \mathcal{H} | \phi \rangle

An upper-bound estimate to the exact energy of the system can then be obtained by a variational procedure 
from this functional.
When introducing a nonlinearity of the type above, the standard functional **does not** lead to an upper-bound 
estimate to the exact energy upon minimization. The appropriate functional is instead given by:

.. math::
   \langle \phi | \mathcal{H}_0 + \frac{1}{2}\mathcal{V}_\mathrm{int} | \phi \rangle

which has the status of a *free energy* (an extensive justification of this fact is given in [Tomasi1994]_).

 
Thus in the PCM framework, the basic energetic quantity is the *free energy of solvation* 
which is conveniently partitioned as follows ([Mennucci2007]_, [Tomasi2005]_, [Amovilli1998]_):

.. math:: 
   \Delta G_\mathrm{sol} = \Delta G_\mathrm{el} + G_\mathrm{cav} + G_\mathrm{dis} + G_\mathrm{rep} + \Delta G_\mathrm{Mm}

where:

- :math:`\Delta G_\mathrm{el}` accounts for the *electrostatic* solute-solvent interaction, arising from mutual polarization in the charge distributions;
- :math:`G_\mathrm{cav}` is the *cavitation* energy, needed to form the molecular cavity inside the continuum representing the solvent;
- :math:`G_\mathrm{dis}` is the *dispersion* energy, due to the solute-solvent dispesion interactions;
- :math:`G_\mathrm{rep}` is the *repulsion* energy, which accounts for Pauli repulsioni;
- :math:`\Delta G_\mathrm{Mm}` is due to molecular motion and accounts for *entropic* contributions to the free energy.

The electrostatic term is, usually, the largest contribution to the solvation energy. Thus, for now, we will
**exclusively** be concerned with its calculation (see also [Amovilli1998]_ for a discussion of the other terms).
In CSMs the calculation of this term requires the solution of the classical Poisson problem nested within the
QM calculation.

.. [Tomasi1994] J. Tomasi and M. Persico; 
   *Molecular Interactions in Solution: An Overview of Methods Based on Continuous Distributions of the Solvent* 
   Chem. Rev. **94**, 2027 (1994)

.. [Mennucci2007] B. Mennucci, R. Cammi (eds.);
   *Continuum Solvation Models in Chemical Physics*
   Wiley, New York (2007)

.. [Tomasi2005] J. Tomasi, B. Mennucci, R. Cammi;
   *Quantum Mechanical Continuum Solvation Models*
   Chem. Rev. **105**, 2999 (2005)

.. [Amovilli1998] C. Amovilli, V. Barone, R. Cammi, E. Cancès, M. Cossi, B. Mennucci, C. S. Pomelli and J .Tomasi;
   *Recent Advances in the Description of Solvent Effects with the Polarizable Continuum Model*
   Adv. Quantum Chem. **32**, 227 (1998)


