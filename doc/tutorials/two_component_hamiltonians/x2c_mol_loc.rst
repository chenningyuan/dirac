

X2C and local X2C
=================

It is beyond the scope of this tutorial to review the ideas behind the exact two-component Hamiltonian (X2C) approach.
For a recent and comprehensive review the reader is referred to [Saue2011]_. The central idea is to generate a two-component Hamiltonian 
that reproduces the positive-energy spectrum of the parent four-component Hamiltonian. This goal may be achieved by 
a unitary transformation :math:`\mathbf{U}` that will formally decouple the large and small components by block-diagonalization:

.. math::

    {\hat{h}_{bd}} = 
    \mathbf{U}^{\dagger}
    \left[\begin{array}{cc}
    \hat{h}_{LL} & \hat{h}_{LS} \\
    \hat{h}_{SL} & \hat{h}_{SS} \end{array}\right] 
    \mathbf{U} = 
    \left[\begin{array}{cc}
    \hat{h}_{++} & 0 \\
         0           & \hat{h}_{--} \end{array}\right],


Popular choices to obtain a unitary transformation matrix :math:`\mathbf{U}` are nicely discussed in Section 3 of Reference
[Saue2011]_ while a detailed abstract of the actual implementation of the X2C approach in Dirac is given in the last part of Section
3.3 in the same review paper. Property operators are transformed automatically on the fly to two-components within 
the X2C approach in Dirac in order to avoid the so-called picture-change errors. For a nice discussion on the latter
topic with the example of charge densities at a heavy nucleus see for example [Knecht2011]_. Transforming the 
the electronic Hamiltonian based on the above prescription 

.. math::

    {\hat{h}^{2c}} = 
    \left[\mathbf{U}^{\dagger}
    \hat{h}^{4c}\mathbf{U}\right]_{++}

it is mandatory to also transform the parent four-component property-operator :math:`\Omega^{4c}` 
using the same unitary transformation matrix :math:`\mathbf{U}`  


.. math::

    {\Omega^{2c}} = 
    \left[\mathbf{U}^{\dagger}
    \Omega^{4c}\mathbf{U}\right]_{++}

since otherwise we would introduce (significant) errors ('picture-change errors') by simply using the approximate relation 

.. math::

    {\Omega^{2c}} \approx \left[\Omega^{4c}\right]_{LL} .


X2C
---

To invoke the X2C Hamiltonian in Dirac which is based on the diagonalization of the 1-electron bare-nucleus Dirac
Hamiltonian, use the keyword .X2C under the input deck \*\*HAMILTONIAN::

  **HAMILTONIAN
  .X2C


This includes by default atomic-mean-field two-electron spin-same-orbit corrections from the
AMFI module. To add spin-other-orbit corrections (Gaunt term of the Breit interaction) use::

  **HAMILTONIAN
  .X2C
  .GAUNT


In the spinfree case one should use::

  **HAMILTONIAN
  .X2C
  .SPINFREE


local X2C
---------

For large systems the diagonalisation of the Dirac Hamiltonian --- which is an essential step in the derivation of the
unitary transformation matrix :math:`\mathbf{U}_{mol}` --- may become a bottleneck for the full molecule (hence to stress
this fact we have introduced the explicit label *mol*).  
This bottleneck can be avoided by extracting the coupling at the level of individual atoms as proposed in a working protocol by Peng
and Reiher [Peng2012]_. The central idea behind this atomic or local scheme is to approximate :math:`\mathbf{U}_{mol}` 
by a superposition of :math:`i=1,...,N_A` atomic (or in principle also fragment) matrices :math:`\mathbf{U}_{i}`

.. math::

   \mathbf{U}_{mol} \approx \sum\limits_{i=1}^{N_{A}} \mathbf{U}_{i} \oplus  \mathbf{U}_{i+1} \oplus \mathbf{U}_{i+2} \oplus \ldots, 


In order to benefit from the above approach in Dirac which at present works *only for C1 symmetry* we need to first carry
out the atomic calculations and save the atomic :math:`\mathbf{U}_{i}` matrices. This can be done as follows. 
We run the atomic X2C transformation (using a cc-pVDZ basis set) using the input::

  **DIRAC
  .TITLE
   atomic X2C - U matrix construction
  **HAMILTONIAN
  .X2C
  **WAVE FUNCTION
  *SCF
  **MOLECULE
  *BASIS
  .DEFAULT 
  cc-pVDZ
  *SYMMETRY
  .NOSYM
  **END OF

Note the use of .NOSYM to enforce C1 symmetry in the atomic run even for an atom. A possible xyz input for an atom (here
C) would look like::

   1
   carbon atom XYZ file
   C 0.0 0.0 0.0

We save the file X2CMAT containing the atomic :math:`\mathbf{U}_{i}` matrix and rename it according to the atomic number
of the respective atom, here X2CMAT.006. For uranium (Z=92) the ending would be X2CMAT.092::

  $ ./pam --inp=atom.inp --mol=atom.xyz --get=X2CMAT
  $ mv X2CMAT X2CMAT.006


After we have created all atomic :math:`\mathbf{U}_{i}` matrices for each individual atom type in our molecule, e.g.
W(CH:math:`_{3}`):math:`_6` we are ready to setup the molecular X2C calculation based on the local/atomic approach of the
decoupling unitary  :math:`\mathbf{U}_{i}` matrices. The relevant part of the input would then be::

  **DIRAC
  .TITLE
   local X2C based on atomic U matrices
  .WAVE FUNCTION
  **HAMILTONIAN
  .X2C
  *X2C
  .fragX2C
  **WAVE FUNCTION
  .SCF
  *SCF
  ...
  ...
  **MOLECULE
  *BASIS
  .DEFAULT 
  cc-pVDZ
  *SYMMETRY
  .NOSYM
  **END OF


This input should be used with the following pam line::

  $ ./pam ... --put=X2CMAT.* 

In the above approach the molecular :math:`\mathbf{U}_{mol}` matrix will be approximated as described above and
corresponds exactly to the DLU scheme proposed by Peng and Reiher [Peng2012]_. It works with both spinfree and
spin-dependent (spin-orbit) Hamiltonian schemes and includes the picture-change transformation of the molecular
four-component property operators :math:`\Omega^{4c}_{mol}`. 
