1s core-ionization of N2 by TD-DFT
==================================

Introduction
------------

We want to study the excitation of a 1s electron of the :math:`N_2` molecule to an empty orbital.
More precisely we shall look at the excitation of an electron from the bonding 
:math:`1s\sigma_g` or anti-bonding :math:`1s\sigma_u` -orbitals to the vacant :math:`2p\pi_g` or :math:`2p\sigma_u` orbitals (see MO-diagram below). 

.. image:: N2_MOdiagram.jpg
   :scale: 40
   :align: center
   :alt: N2_MOdiagram

Note that this diagram does not take spin-orbit into account, but we shall consider
this interaction later on. Let us first consider the possible final states. 
One electron leaves from one of four spin-orbitals and enters one of six spin-orbitals.
This gives 24 determinants which translates into the following states:

+-----------------------------------+------------------------+
| Configuration                     | States                 |
+===================================+========================+
| :math:`1s\sigma_g^{-1}2p\pi_g`    | :math:`^{1,3}\Pi_g`    |
+-----------------------------------+------------------------+
| :math:`1s\sigma_g^{-1}2p\sigma_u` | :math:`^{1,3}\Sigma_u` |
+-----------------------------------+------------------------+
| :math:`1s\sigma_u^{-1}2p\pi_g`    | :math:`^{1,3}\Pi_u`    |
+-----------------------------------+------------------------+
| :math:`1s\sigma_u^{-1}2p\sigma_u` | :math:`^{1,3}\Sigma_g` |
+-----------------------------------+------------------------+

Spin-orbit free calculation
---------------------------

Preparing the input files
~~~~~~~~~~~~~~~~~~~~~~~~~

We employ the following molecular input file `N2.mol` ::

  DIRAC
  N2  
  cc-pVTZ basis
  C   1              A
          7.    2
  N1         0.0        0.0       0.00000
  N2         0.0        0.0       1.09800
  LARGE BASIS cc-pVTZ
  FINISH                                                                                   

Here we do not provide any symmetry information, meaning that we ask DIRAC to detect it.
DIRAC will find that the full group is :math:`D_{\infty h}`. With spin-orbit coupling DIRAC
will then activate linear supersymmetry, but in the spin-orbit free case it will simply use
the highest Abelian single point group, that is :math:`D_{2h}`. For the final states DIRAC
will  employ the *total* symmetry, that is the combined spin and spatial symmetry.
Here we shall keep in mind that the singlet spin function is totally symmetric (:math:`A_g`), whereas the
triplet spin functions transform as rotations. We know the triplet functions as:

.. math::

   T_{-1} = \alpha_1\alpha_2;\quad T_0=\frac{1}{\sqrt{2}}\left(\alpha_1\beta_2-\beta_1\alpha_2\right);\quad T_{+1}=\beta_1\beta_2

but for our purposes it will be more convenient to form the combinations:

.. math::

   T_x=\frac{1}{\sqrt{2}}\left(T_{-1}-T_{+1}\right);\quad T_y=\frac{i}{\sqrt{2}}\left(T_{-1}+T_{+1}\right);\quad T_z=T_0

which transform as rotations :math:`R_x\left(B_{3g}\right)` , :math:`R_y\left(B_{2g}\right)` and :math:`R_z\left(B_{1g}\right)` , respectively. 

We can now set up the following correlation of states:

+------------------------+--------------------------------+------------------------+----------------------------------------------------------------------------+  
| State                  | Spin                           | Spatial                | Spin :math:`\otimes` Spatial                                               |
+========================+================================+========================+============================================================================+  
|  :math:`^1\Sigma_g`    | :math:`A_g`                    | :math:`A_g`            | :math:`A_g`                                                                |
+------------------------+--------------------------------+------------------------+----------------------------------------------------------------------------+  
|  :math:`^1\Sigma_u`    | :math:`A_g`                    | :math:`B_{1u}`         | :math:`B_{1u}`                                                             |
+------------------------+--------------------------------+------------------------+----------------------------------------------------------------------------+  
|  :math:`^1\Pi_{x,y;g}` | :math:`A_g`                    | :math:`B_{3g}, B_{2g}` | :math:`B_{3g}, B_{2g}`                                                     |
+------------------------+--------------------------------+------------------------+----------------------------------------------------------------------------+  
|  :math:`^1\Pi_{x,y;u}` | :math:`A_g`                    | :math:`B_{3u}, B_{2u}` | :math:`B_{3u}, B_{2u}`                                                     |
+------------------------+--------------------------------+------------------------+----------------------------------------------------------------------------+  
|  :math:`^3\Sigma_g`    | :math:`B_{3g}, B_{2g}, B_{1g}` | :math:`A_g`            | :math:`B_{3g}, B_{2g}, B_{1g}`                                             |
+------------------------+--------------------------------+------------------------+----------------------------------------------------------------------------+  
|  :math:`^3\Sigma_u`    | :math:`B_{3g}, B_{2g}, B_{1g}` | :math:`B_{1u}`         | :math:`B_{2u}, B_{3u}, A_{u}`                                              |
+------------------------+--------------------------------+------------------------+----------------------------------------------------------------------------+  
|  :math:`^3\Pi_{x,y;g}` | :math:`B_{3g}, B_{2g}, B_{1g}` | :math:`B_{3g}, B_{2g}` | :math:`\left(A_g, B_{1g}, B_{2g}\right), \left(B_{1g}, A_g, B_{3g}\right)` |
+------------------------+--------------------------------+------------------------+----------------------------------------------------------------------------+  
|  :math:`^3\Pi_{x,y;u}` | :math:`B_{3g}, B_{2g}, B_{1g}` | :math:`B_{3u}, B_{2u}` | :math:`\left(A_u, B_{1u}, B_{2u}\right), \left(B_{1u}, A_u, B_{3u}\right)` |
+------------------------+--------------------------------+------------------------+----------------------------------------------------------------------------+  

Counting total symmetries we find the 24 microstates are evenly distributed amongst the eight irreps of :math:`D_{2h}`:


+----------------+-------------------------------------------------------------+-----------------------+
| Irrep          | Core-ionized state                                          |                       |
+================+=============================================================+=======================+
| :math:`A_g`    | :math:`^1\Sigma_g, ^{3(x)}\Pi_{x;g}, ^{3(y)}\Pi_{y;g}`      | :math:`x^2, y^2, z^2` |
+----------------+-------------------------------------------------------------+-----------------------+
| :math:`B_{3u}` | :math:`^{1}\Pi_{x;u}, ^{3(y)}\Sigma_u, ^{3(z)}\Pi_{y;u}`    | :math:`x`             |
+----------------+-------------------------------------------------------------+-----------------------+
| :math:`B_{2u}` | :math:`^{1}\Pi_{y;u}, ^{3(x)}\Sigma_u, ^{3(z)}\Pi_{x;u}`    | :math:`y`             |
+----------------+-------------------------------------------------------------+-----------------------+
| :math:`B_{1g}` | :math:`^{3(z)}\Sigma_g, ^{3(y)}\Pi_{x;g}, ^{3(x)}\Pi_{y;g}` | :math:`xy`            |
+----------------+-------------------------------------------------------------+-----------------------+
| :math:`B_{1u}` | :math:`^1\Sigma_u, ^{3(y)}\Pi_{x;u}, ^{3(x)}\Pi_{y;u}`      | :math:`z`             |
+----------------+-------------------------------------------------------------+-----------------------+
| :math:`B_{2g}` | :math:`^{1}\Pi_{y;g}, ^{3(y)}\Sigma_g, ^{3(z)}\Pi_{x;g}`    | :math:`xz`            |
+----------------+-------------------------------------------------------------+-----------------------+
| :math:`B_{3g}` | :math:`^{1}\Pi_{x;g}, ^{3(x)}\Sigma_g, ^{3(z)}\Pi_{y;g}`    | :math:`yz`            |
+----------------+-------------------------------------------------------------+-----------------------+   
| :math:`A_{u}`  | :math:`^{3(z)}\Sigma_u, ^{3(x)}\Pi_{x;u}, ^{3(y)}\Pi_{y;u}` | :math:`xyz`           |
+----------------+-------------------------------------------------------------+-----------------------+

From these considerations we now set up the following menu file for our calculation::

  **DIRAC
  .WAVE FUNCTION
  .ANALYZE
  .PROPERTIES
  **INTEGRALS
  *READIN
  .UNCONTRACT
  **HAMILTONIAN
  .SPINFREE
  .DFT
  PBE0
  **WAVE FUNCTION
  .SCF
  *SCF
  .CLOSED SHELL
  6 8
  **ANALYZE
  .MULPOP
  *MULPOP
  .VECPOP
  1..5
  1..5
  **PROPERTIES
  .PRINT
  1
  *EXCITATION ENERGIES
  .OCCUP
  1
  1
  .VIRTUA
  4,5
  5
  .EXCITATIONS  !Ag
  1 3
  .EXCITATIONS  !B3u
  2 3
  !.EXCITATIONS  !B2u
  !3 3
  .EXCITATIONS  !B1g
  4 3
  .EXCITATIONS  !B1u
  5 3
  !.EXCITATIONS  !B2g
  !6 3
  .EXCITATIONS  !B3g
  7 3
  .EXCITATIONS  !Au
  8 3
  .OPERATOR
   XDIPLEN
  .OPERATOR
   YDIPLEN
  .OPERATOR
   ZDIPLEN
  .ANALYZE
  .SKIPEP
  *END OF INPUT

In the :ref:`*SCF` section we give the electron occupation of :math:`N_2`: 6 and 8 electrons in *gerade* and *ungerade* orbitals, respectively.
We also ask for a Mulliken population analysis (:ref:`.MULPOP`) for the occupied orbitals and the orbitals involved in the core excitation.

Let us  now look at how we set up the calculation of excitation energies under `*EXCITATION ENERGIES`. We have seen that there are three
excitations per boson irrep. Note that the numbering of irreps follow what you for instance find in the :math:`D_{2h}` direct product
table in the output::

  |   | Ag   B3u  B2u  B1g  B1u  B2g  B3g  Au 
  -----+----------------------------------------
  Ag  | Ag   B3u  B2u  B1g  B1u  B2g  B3g  Au 
  B3u | B3u  Ag   B1g  B2u  B2g  B1u  Au   B3g
  B2u | B2u  B1g  Ag   B3u  B3g  Au   B1u  B2g
  B1g | B1g  B2u  B3u  Ag   Au   B3g  B2g  B1u
  B1u | B1u  B2g  B3g  Au   Ag   B3u  B2u  B1g
  B2g | B2g  B1u  Au   B3g  B3u  Ag   B1g  B2u
  B3g | B3g  Au   B1u  B2g  B2u  B1g  Ag   B3u
  Au  | Au   B3g  B2g  B1u  B1g  B2u  B3u  Ag 

Note also that we skip excitations in :math:`B_{2u}` and :math:`B_{2g}`, since they are related by symmetry
to the excitations of :math:`B_{3u}` and :math:`B_{3g}`, respectively. 

If nothing further is specified the excitation energies are calculated by a "bottoms-up" approach and so we will
get valence excitations only, since the core-excitations are much higher in energy. We therefore restrict the
excitations to the occupied :math:`1s\sigma_g` and :math:`1s\sigma_u` orbitals. 
It is also possible to limit the virtual orbitals to the empty orbitals of the *2p* manifold::

  .VIRTUA                                                                                            
  4,5                                                                                                
  5                                                                                                  

but this is in fact an unnecessary restriction.

We furthermore ask for transition moments to be calculated with respect to the component of the dipole moment operator.
These will be non-zero only for excitations in irreps :math:`B_{3u}, B_{2u}` and :math:`B_{1u}`.
Finally we ask for analysis of what orbitals contribute to the various excitations. For this the Mulliken population
analysis may come in handy as reference.

Looking at the output
~~~~~~~~~~~~~~~~~~~~~

After running the calculation, let us now look at the output. We take the dipole-allowed transitions of :math:`B_{3u}` as
an example::

  *** Excitations of boson symmetry 2 : B3u
  
  --------------------------------------------------------------------------------
   Excitation no.   1 excitation energy    14.28739014 a.u.    4.06E-07 (converged)
   Nonrel. sym.:   TzB2u                      388.7797 eV
   *** Transition moments <0|A|n> for this excitation:
  
   A - XDIPLEN              B3u  T+         0.00000000D+00 a.u.
  --------------------------------------------------------------------------------
   Excitation no.   2 excitation energy    14.32916125 a.u.    3.86E-07 (converged)
   Nonrel. sym.:   S0B3u                      389.9163 eV
   *** Transition moments <0|A|n> for this excitation:
  
   A - XDIPLEN              B3u  T+        -0.108039773420 a.u.
  
       -> dipole oscillator strength (length)   :  0.11150564
       -> dipole oscillator strength (velocity) :  0.00000000
  
       Total dipole radiation rate (length)   7.35613E+11 s-1
        - corresponding radiation life time   1.35941E-12 s
  --------------------------------------------------------------------------------
   Excitation no.   3 excitation energy    14.69518720 a.u.    1.65E-06 (converged)
   Nonrel. sym.:   TyB1u                      399.8764 eV
   *** Transition moments <0|A|n> for this excitation:
  
   A - XDIPLEN              B3u  T+         0.00000000D+00 a.u.
  

There is sufficient symmetry in the calculation (symmetry distinct rotation) to allow DIRAC to pinpoint
the symmetry of the core-ionized state (see the entry `Nonrel. sym.`). The two triplet states have strictly
zero transition moment (and thus intensity), whereas the :math:`^1\Pi_u` state has a transition moment of 
0.1080 a.u.

Looking further down in the output we find dominant inactive and virtual orbitals. The first excited state :math:`^3\Pi_u`,
at 388.78 eV, is dominated by the excitation `1(i:E1u) --->      4(v:E1g)`, which, as can be inferred from the Mulliken 
population analysis, corresponds to :math:`1s\sigma_u \rightarrow 2p\pi_{y;g}`. The second excited state :math:`^1\Pi_u`, 
at 389.91 eV, corresponds to `1(i:E1u) --->      5(v:E1g)` (:math:`1s\sigma_u \rightarrow 2p\pi_{x;g}`), 
whereas the third excited state :math:`^3\Sigma_u`, at 399.88 eV, is dominated by `1(i:E1g) --->      5(v:E1u)` 
(:math:`1s\sigma_g \rightarrow 2p\sigma_u`).

Including spin-orbit
--------------------

Spin-orbit is included by simply commenting out the keyword :ref:`.SPINFREE` in the input above::

  **HAMILTONIAN
  !.SPINFREE

The output concerning excited states of total symmetry :math:`B_{2u}` in :math:`D_{2h}` now reads::

   *** Excitations of boson symmetry 2 : B3u
  
   --------------------------------------------------------------------------------
   Excitation no.   1 excitation energy    14.28738705 a.u.    4.44E-07 (converged)
   Linear symmetry:  1u                       388.7796 eV
   *** Transition moments <0|A|n> for this excitation:
  
   A - XDIPLEN              B3u  T+        -0.80477959D-03 a.u.
  
       -> dipole oscillator strength (length)   :  6.16900972E-06
       -> dipole oscillator strength (velocity) :  0.00000000E+00

       Total dipole radiation rate (length)   4.04606E+07 s-1
        - corresponding radiation life time   2.47154E-08 s
   --------------------------------------------------------------------------------
   Excitation no.   2 excitation energy    14.32916239 a.u.    3.86E-07 (converged)
   Linear symmetry:  1u                       389.9163 eV
   *** Transition moments <0|A|n> for this excitation:
  
   A - XDIPLEN              B3u  T+         0.108036713021 a.u.
  
       -> dipole oscillator strength (length)   :  0.11149933
       -> dipole oscillator strength (velocity) :  0.00000000
  
       Total dipole radiation rate (length)   7.35571E+11 s-1
        - corresponding radiation life time   1.35949E-12 s
   --------------------------------------------------------------------------------
   Excitation no.   3 excitation energy    14.69518602 a.u.    1.74E-06 (converged)
   Linear symmetry:  1u                       399.8763 eV
   *** Transition moments <0|A|n> for this excitation:
  
   A - XDIPLEN              B3u  T+         0.49111365D-05 a.u.
  
       -> dipole oscillator strength (length)   :  2.36291361E-10
       -> dipole oscillator strength (velocity) :  0.00000000E+00
  
       Total dipole radiation rate (length)    1639.5     s-1
        - corresponding radiation life time   6.09946E-04 s

To understand the effect of spin-orbit coupling we note that in this irrep we pick out the :math:`1_u` of the
original three :math:`^{2S+1}\Lambda_u` states, thus enabling all three states to mix. In practice the two original
triplet states steal intensity from the singlet state, such that the transition moment of the latter state is
reduced from 0.1080398 a.u. to 0.108037 a.u.. This change is not very spectacular since the nitrogen molecule is composed of
light atoms for which relativistic effects are not very strong. We can mimic a more strongly relativistic system by reducing
the speed of light to e.g. 20 a.u.::

  **GENERAL                                                                                           
  .CVALUE                                                                                             
  20.0D0

The transition moment of the singlet :math:`^1\Pi_u` state is then reduced from 0.102809 a.u. to 0.974999 a.u. 
upon inclusion of spin-orbit coupling.

Simulating the core-excitation spectrum using complex response
--------------------------------------------------------------

`Gas Phase Core Excitation Database <http://unicorn.mcmaster.ca/corex/cedb-title.html>`_
