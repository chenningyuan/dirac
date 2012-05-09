SCF niobium atom
================

DIRAC presently does not do open-shell Hartree-Fock in a strict sense, rather
average-of-configuration, which amounts to the optimization of the
average energy of a set of configurations (or determinants) generated
from the specification of a given number of open shells and their
electron occupation. As an example of such a calculation, showing in
particular tricks to help convergence, we consider the niobium atom.

The ground state configuration of niobium is
[Kr].4d\ :sup:`4`\ .5s\ :sup:`1`\  and so a reasonable input for an
average Hartree-Fock calculation would be

::

    **DIRAC
    .TITLE
     Nb atom
    .WAVE FUNCTION
    **WAVE FUNCTION
    .SCF
    *SCF
    .CLOSED SHELL
     18 18
    .OPEN S
     2
     4/10,0
     1/2,0
    *END OF

Here we create separate open shells for the four *4d* electrons and the
single *5s* electron to reduce the number of states in the average.
Please note that by defaul the (SS\|SS) type Coulomb integrals are not
calculated, but the energy corrected (see keyword
`.DOSSSS </dirac/index.php/Manual:HAMILTONIAN#.DOSSSS>`_), which reduces
the computational cost significantly!

We combine this menu file with the `molecular
file </dirac/index.php/Manual:Basis:MOL>`_:

::

    INTGRL
    Nb atom
    test
    C   1
           41.    1
    Nb      0.0000000000        0.0000000000        0.0000000000
    LARGE BASIS dyall.v3z
    FINISH

When using the basis set library one should normally add

::

    **INTEGRALS
    *READIN
    .UNCONTRACT

to the menu file to `decontract the
basis </dirac/index.php/Manual:INTEGRALS:READIN#.UNCONTRACT>`_, but in
this particular case the basis is already decontracted.

Unfortunately this calculation does not converge. The reason is that the
open-shell *4d* and *5s* orbitals mix. This is easily seen by performing
a `Mulliken population
analysis </dirac/index.php/Manual:ANALYZE:MULPOP>`_. Since we only want
to distinguish atomic *4d* and *5s* orbitals we use the
`.LABDEF </dirac/index.php/Manual:ANALYZE:MULPOP#.LABDEF>`_ option to
simplify the output from the Mulliken population analysis. In the SCF
output we find the labels

::

                                   GETLAB: SO-labels
                                   -----------------
      * Large components:   35
        1  L Ag Nb s        2  L Ag Nb dxx      3  L Ag Nb dyy      4  L Ag Nb dzz      5  L Ag Nb g400     6  L Ag Nb g220
        7  L Ag Nb g202     8  L Ag Nb g040     9  L Ag Nb g022    10  L Ag Nb g004    11  L B3uNb px      12  L B3uNb fxxx
       13  L B3uNb fxyy    14  L B3uNb fxzz    15  L B2uNb py      16  L B2uNb fxxy    17  L B2uNb fyyy    18  L B2uNb fyzz
       19  L B1gNb dxy     20  L B1gNb g310    21  L B1gNb g130    22  L B1gNb g112    23  L B1uNb pz      24  L B1uNb fxxz
       25  L B1uNb fyyz    26  L B1uNb fzzz    27  L B2gNb dxz     28  L B2gNb g301    29  L B2gNb g121    30  L B2gNb g103
       31  L B3gNb dyz     32  L B3gNb g211    33  L B3gNb g031    34  L B3gNb g013    35  L Au Nb fxyz

from which we generate the input:

::

    *MULPOP
    .LABDEF
    5
    s           1
    p           11,15,23
    d           2,3,4,19,27,31
    f           12,13,14,16,17,18,24,25,26,35
    g           5,6,7,8,9,10,20,21,22,28,29,30,32,33,34

The Mulliken population analysis for the open-shell orbitals show
(slightly edited):

::

    * Electronic eigenvalue no. 10: -0.2345545198473       (Occupation : f = 0.4000)  m_j=  1/2
    ==========================================================================================

     Gross     Total   |    s              d
     -----------------------------------------------------
     alpha    0.9952  |      0.5066         0.4886
     beta     0.0048  |      0.0000         0.0048

    * Electronic eigenvalue no. 11: -0.2341282163131       (Occupation : f = 0.4000)  m_j= -3/2
    ==========================================================================================

    Gross     Total   |    d
    --------------------------------------
     alpha    0.0668  |      0.0668
     beta     0.9332  |      0.9331

    * Electronic eigenvalue no. 12: -0.2196186152479       (Occupation : f = 0.4000)  m_j= -3/2
    ==========================================================================================

    Gross     Total   |    d
    --------------------------------------
     alpha    0.9331  |      0.9331
     beta     0.0669  |      0.0668

    * Electronic eigenvalue no. 13: -0.2189244235127       (Occupation : f = 0.4000)  m_j=  5/2
    ==========================================================================================

    Gross     Total   |    d
    --------------------------------------
     alpha    0.9999  |      0.9999
     beta     0.0001  |      0.0000

    * Electronic eigenvalue no. 14: -0.2182274767504       (Occupation : f = 0.4000)  m_j=  1/2
    ==========================================================================================

    Gross     Total   |    s              d
    -----------------------------------------------------
     alpha    0.9919  |      0.4932         0.4987
     beta     0.0081  |      0.0000         0.0080

    * Electronic eigenvalue no. 15: -0.2173465762720       (Occupation : f = 0.5000)  m_j=  1/2
    ==========================================================================================

    Gross     Total   |    s              d
    -----------------------------------------------------
     alpha    0.0129  |      0.0002         0.0126
     beta     0.9871  |      0.0000         0.9871

and we see significant mixing of *4d* and *5s* orbitals for
m\ :sub:`j`\  = 1/2. This mixing is forbidden from atomic symmetry. The
problem is that we run the atom only in linear supersymmetry and so the
mixing is introduced by numerical noise due to the energetic proximity
of these orbitals.

An efficient way to handle such a case is to strip off the open-shell
electrons and simply calculate first Nb\ :sup:`5+`\  with the electron
configuration of Kr:

::

    .CLOSED SHELL
     18 18

This calculation converges smoothly and by extending the Mulliken
population analysis to include the virtual orbitals we find in fermion
ircop E\ :sub:`1g`\  that orbitals 10 - 14 are purely *d*, whereas
orbital 15 is pure *s*. Restarting with the coefficients obtained from
this calculation we now converge to the ground state configuration of
the neutral atom by imposing vector selection by overlap:

::

    *SCF
    .CLOSED SHELL
     18 18 
    .OPEN S
     2
     4/10,0
     1/2,0
    .OVLSEL
    .NODYNSEL

This assures that the orbitals obtained in the calculation on
Nb\ :sup:`5+`\  stay in place. Note that DIRAC will always put
open-shell orbitals after closed-shell ones. If orbitals are not in the
right place, overlap selection can be combined with
`reordering </dirac/index.php/Manual:WAVE_FUNCTION#.REORDER_MO.27S>`_.

