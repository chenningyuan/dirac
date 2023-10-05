:orphan:
 

===================
Projection analysis
===================

------
Theory
------

A serious flaw of Mulliken population analysis is its basis set sensitivity. This is well known from
the non-relativistic domain. In the relativistic domain a further disadvantage is that the 
expansion of large and small components in separate basis sets means large and small component 
populations come separately and it is not easy to see how they are connected. Furthermore, if
the large and small component are expanded in scalar basis functions, as is the case of the 
DIRAC code, then the population analysis can not distinguish distributions from different 
spin-orbit components of atomic orbitals, for instance, :math:`p_{1/2}` and :math:`p_{3/2}`.

Projection analysis :cite:`Dubillard2006` eliminates all these inconveniences. It is based on the simple concept of
molecular orbitals as linear combination of atomic orbitals (LCAO) and allows analysis of 
electronic structure in terms of well-defined atomic (or fragment) orbitals. 

Projection analysis requires as input a set of pre-calculated orbitals {:math:`\psi^A_p`} 
for the constituent atoms (or fragment) of the molecule. Here :math:`A` refers to 
the atomic center (or fragment) and :math:`p` is on orbital index for that center. 
The molecular orbitals are then expanded as

.. math::

   \left|\psi^{MO}_i\right> = \sum_A \sum_{p\in A} \psi^A_p c^A_{pi} + \left|\psi^{pol}_i\right>

Typically one will select the occupied orbitals of the constituent atoms. They are not 
guaranteed to fully span the molecular orbitals and therefore the expansion includes
their orthogonal complement :math:`\left|\psi^{pol}_i\right>`, denoted the *polarization contribution*.
The expansion coefficients :math:`c^A_{pi}` are found by solving the equation

.. math::

   \left<\psi^A_p|\psi^{MO}_i\right> = \sum_B\sum_{q\in B}\left<\psi^A_p|\psi^B_q\right>c^B_{qi}

obtained by projection.

For projection analysis, arbitrary fragments can in principle be used. There is, however, an important restriction 
in that DIRAC can only do analysis in terms of fragments complying to the molecular point group. For instance, 
the two hydrogens in water are related by symmetry in the full :math:`C_{2v}` molecular point group. They then have 
to be treated as a single fragment, or they can be separated by lowering the symmetry. Starting from DIRAC14 
projection analysis can be carried out in terms of individual atoms using the :ref:`PROJECTION_.ATOMS` keyword, as shown in 
the example below.

------------------
Example 1: Methane
------------------

Introduction
------------

As a first simple example we shall consider the methane molecule.
The molecular input file ``CH4.xyz`` (:download:`download <../../../../test/tutorial_projection_analysis/CH4.xyz>`) is

.. literalinclude:: ../../../../test/tutorial_projection_analysis/CH4.xyz

Suppose that we have run a DFT(PBE) calculation using::

   pam --inp=PBE  --mol=CH4.xyz 

and the menu file ``PBE.inp`` (:download:`download <../../../../test/tutorial_projection_analysis/PBE.inp>`) is

.. literalinclude:: ../../../../test/tutorial_projection_analysis/PBE.inp

Preparing for analysis
----------------------

We would now like to analyze the electronic structure of the methane molecule in terms of 
the orbitals of the constituent atoms. We first generate the atomic fragments. 
We recommend calculating the atoms in the highest possible 
symmetry and the converting the coefficients to no (:math:`C_1`) symmetry. 
This allows a very precise identification
of the individual atomic orbitals (including :math:`m_j` quantum number). 
We accordingly run calculations where we save MO files, each as the unique file::

   pam --inp=H  --mol=H 
   pam --inp=C  --mol=C 

with corresponding input files ``H.inp`` (:download:`download <../../../../test/tutorial_projection_analysis/H.inp>`),
``H.mol`` (:download:`download <../../../../test/tutorial_projection_analysis/H.mol>`), 
``C.inp`` (:download:`download <../../../../test/tutorial_projection_analysis/C.inp>`),
and ``C.mol`` (:download:`download <../../../../test/tutorial_projection_analysis/C.mol>`).


Running the analysis
--------------------

We are now ready for the projection analysis itself. We prepare the input file ``prj.inp`` 
(:download:`download <../../../../test/tutorial_projection_analysis/prj.inp>`)

.. literalinclude:: ../../../../test/tutorial_projection_analysis/prj.inp

In the ``*PROJEC`` section we specify the two atomic types identified by the name of 
their checkpoint files as well as an orbital string (see :ref:`orbital_strings`) listing 
what orbitals to include for each atom. DIRAC will assume that each atom has been calculated 
in its proper basis.

We now run the projection analysis using::

   pam --inp=prj --mol=CH4.xyz --put "C_C.h5=C.h5 H_H.h5=H.h5 PBE_CH4.h5=CHECKPOINT.h5" --noh5

Note that we set the pam flag ``--noh5`` to avoid generating a new HDF5 checkpoint file, since nothing new is added to the existing one.   


Interpreting the output
-----------------------

The first thing to look for in the output 
(:download:`download <../../../../test/tutorial_projection_analysis/result/prj_CH4.out>`) is the section::
   
   * Total gross contributions:
      AFCX 1               6.4301   E -     6.4301   P -     0.0000
      AFHX 1               0.8716   E -     0.8716   P -     0.0000
      AFHX 2               0.8716   E -     0.8716   P -     0.0000
      AFHX 3               0.8716   E -     0.8716   P -     0.0000
      AFHX 4               0.8716   E -     0.8716   P -     0.0000
      Polarization:     0.0834

It shows the gross population (in Mulliken sense, but in terms of the fragment orbitals) of each atom 
and then, most importantly, the gross population of the polarization contribution. In this example it
is small and perfectly acceptable. If it becomes too large the analysis looses meaning because it 
impies that significant electron density has not been assigned to any fragment. In the case of
significant polarization contribution one can either increase the number of atomic orbitals used
in the analysis or turn on *repolarization*. For the latter it suffices to add the :ref:`PROJECTION_.POLREF` keyword
(see the corresponding :download:`input file <../../../../test/tutorial_projection_analysis/prj_polref.inp>` )
which activates *Intrinsic Atomic Orbitals (IAOs)*, :cite:`Knizia2013`, 
thus eliminating completely the polarization contribution. 
We then get (see the corresponding :download:`output <../../../../test/tutorial_projection_analysis/result/prj_polref_CH4.out>` file )::

  * Total gross contributions:
     AFCX 1               6.4782   E -     6.4782   P -     0.0000
     AFHX 1               0.8804   E -     0.8804   P -     0.0000
     AFHX 2               0.8804   E -     0.8804   P -     0.0000
     AFHX 3               0.8804   E -     0.8804   P -     0.0000
     AFHX 4               0.8804   E -     0.8804   P -     0.0000
     Polarization:     0.0000

from which we can conclude that at this level of theory the carbon and hydrogen charges in methane 
are -0.48e and +0.12e, respectively.

Having concluded that the polarization contribution is acceptable, we can proceed with the analysis.
We turn to this section::

   * Total reference orbital contributions:
       AFCX 1  E     1       1.99999     -0.10660951E+02
       AFCX 1  E     2       1.27874     -0.97828469E+00
       AFCX 1  E     3       1.05077     -0.65580701E+00
       AFCX 1  E     4       1.05030     -0.65541446E+00
       AFCX 1  E     5       1.05030     -0.65541446E+00
       AFHX 1  E     1       0.87162     -0.23118508E+00
       AFHX 2  E     1       0.87162     -0.23118508E+00
       AFHX 3  E     1       0.87162     -0.23118508E+00
       AFHX 4  E     1       0.87162     -0.23118508E+00
   
which gives the accumulated gross population of each fragment orbital. We can identify 
the fragments from the name of the checkpoint file and its orbital energy. For instance
the first orbital on the list is clearly carbon *1s* followed by carbon *2s*. Next comes three
orbitals that are almost degenerate. These are cleary the carbon *2p* orbitals, the first one
being :math:`2p_{1/2}` followed by the two :math:`2p_{3/2}`. By furthermore looking into the
output of the carbon atom calculation (note that we asked for Mulliken population analysis in the
input) we find that the first :math:`2p_{3/2}` orbital has :math:`m_j=1/2` and the 
second :math:`m_j=-3/2`. From these gross population we can deduce the atomic configurations of 
the atoms in the molecule. For carbon we find :math:`[He]2s^{1.3}2p^{3.2}` and for hydrogen 
:math:`1s^{0.9}`.

There is also detailed output for each molecular orbital in the output. For the third molecular
orbital we find for instance::

   * Electronic eigenvalue nr.  3: -0.3411330794649       (Occupation : f = 1.0000)
   ================================================
       Orbital            Total           Eigenvalue       Kramers partner 1   Kramers partner 2
       AFCX 1  E     3       0.57945     -0.65580701E+00   (-1.0000,-0.0000)   ( 0.0000, 0.0000)
       AFHX 1  E     1       0.31559     -0.23118508E+00   ( 0.0000,-0.5774)   (-0.5774,-0.5774)
       AFHX 2  E     1       0.31559     -0.23118508E+00   ( 0.0000, 0.5774)   ( 0.5774,-0.5774)
       AFHX 3  E     1       0.31559     -0.23118508E+00   ( 0.0000, 0.5774)   (-0.5774, 0.5774)
       AFHX 4  E     1       0.31559     -0.23118508E+00   ( 0.0000,-0.5774)   ( 0.5774, 0.5774)
   * Gross contributions:
      AFCX 1               0.5254   E -     0.5254   P -     0.0000
      AFHX 1               0.1165   E -     0.1165   P -     0.0000
      AFHX 2               0.1165   E -     0.1165   P -     0.0000
      AFHX 3               0.1165   E -     0.1165   P -     0.0000  
      AFHX 4               0.1165   E -     0.1165   P -     0.0000
      Polarization:     0.0085

To fully understand the output we have to keep in mind that all calculations are Kramers-restricted
such that all orbitals come in pairs. The expansion of the molecular orbitals should therefore
rather be written as

.. math::

   \left|\psi^{MO}_i\right> = \sum_A \sum_{p\in A} \left(\psi^A_p c^A_{pi} + \psi^A_\overline{p} c^A_{\overline{p}i}\right) + \left|\psi^{pol}_i\right>

where :math:`c^A_{pi}` and :math:`c^A_{\overline{p}i}` is the expansion coefficient of Kramers partner 
1 and 2, respectively. Let us now define the absolute value

.. math:: \left|c^A_{Pi}\right| = \sqrt{\left|c^A_{pi}\right|^2+\left|c^A_{\overline{p}i}\right|^2}

This is the real number reported under the heading ``Total``. Under the heading 
``Kramers partner 1`` is reported the complex number

.. math:: z_1 = \left(Re(z_1),Im(z_1)\right) = \frac{c^A_{pi}}{\left|c^A_{Pi}\right|},

and analogously for Kramers partner 2, from which we get the detailed decomposition of the molecular
orbital. We find that this molecular orbital is composed of the carbon :math:`2p_{1/2}` orbital and
the hydrogen :math:`1s` orbitals.


Exercises
---------

1. What is the effect of increasing the C-H bond distance on the atomic configurations of the C,H atoms in the methane molecule and 
   on the atomic charges of the constituent atoms ? Collect your calculated data in a table.

2. How different basis sets affect the atomic configurations of the C,H atoms in the methane molecule and 
   the atomic charges of the same constituent atoms ?

3. Investigate the effect of various DFT functionals (see :ref:`*DFT`) on the atomic configurations of the C and H atoms 
   in the methane molecule and on atomic charges of these constituent atoms.


-----------------
Example 2: Uranyl
-----------------

Introduction
------------

We next turn to a molecule with a heavy atom and a more complicated electronic structure, namely uranyl 
:math:`UO_2^{2+}`. Here the two oxygens are related by symmetry, and so one would at first consider 
running the projection analysis with lower (or no) symmetry or using the :ref:`PROJECTION_.ATOMS` keyword, as 
in the methane example above. However, we can reduce the symmetry from :math:`D_{\infty h}` to :math:`C_{\infty v}` 
by introducing a ghost center. Our molecular input file then looks like

.. literalinclude:: UO2.mol

Now the oxygen atoms are no longer connected by symmetry and we can carry out the projection 
analysis in :math:`C_{\infty v}` symmetry. 

Preparing the analysis
----------------------

Using the menu file ``HF.inp``

.. literalinclude:: UO2_HF.inp

we run 4-component Hartree-Fock::

   pam --mw=150 --inp=HF --mol=UO2 

Coefficients are saved on the h5 checkpoint file. Next we generate atomic reference orbitals. For oxygen this is quite straightforward::

   pam --inp=O --mol=O 

using molecule input ``O.mol``

.. literalinclude:: O.mol

and the menu file ``O.inp``

.. literalinclude:: O.inp

For the uranium we must be a bit more careful since the ground state configuration is an open-shell one: :math:`[Rn]5f^36d^17s^2`.
Using the molecular input ``U.mol`` 

.. literalinclude:: U.mol

and the inp file ``U.inp``

.. literalinclude:: U_KPSELE.inp

Note that :ref:`SCF_.KPSELE` is needed for getting the convergence in this case.
Convergence is now straightforward ::

   pam --inp=U --mol=U 

and the correct configuration is confirmed by Mulliken population analysis. 

Running the analysis
--------------------

We are now ready to do the projection analysis. We prepare the input file ``prj.inp``

.. literalinclude:: prj_UO2.inp

and prepare the checkpoint files::

   $ ln -s U_U.h5 U.h5
   $ ln -s O_O.h5 O1.h5
   $ ln -s O_O.h5 O2.h5
   $ ln -s HF_UO2.h5 CHECKPOINT.h5

We run the projection analysis::

   pam --incmo --inp=prj --mol=UO2 --copy="U.h5 O1.h5 O2.h5"

Interpreting the output
-----------------------

As in the methane case above we first look at the polarization contribution::

  * Total gross contributions:
     U.h5                89.0601   E -    89.0601   P -     0.0000
     O1.h5                8.2929   E -     8.2929   P -     0.0000
     O2.h5                8.2929   E -     8.2929   P -     0.0000
     Polarization:     0.3540

We see that it is more sizable than before. If we eliminate the polarization contribution using the :ref:`PROJECTION_.POLREF` keyword, we obtain::

  * Total gross contributions:
     U.h5                89.1616   E -    89.1616   P -     0.0000
     O1.h5                8.4192   E -     8.4192   P -     0.0000
     O2.h5                8.4192   E -     8.4192   P -     0.0000
     Polarization:     0.0000

The atomic gross populations tell us that the oxygen atoms has a charge -0.42, whereas the uranium atom has charge +2.84, considerably far from the formal oxidation state +VI.

We next turn to the electron configuration of the atoms in the molecule by looking at the section::

  * Total reference orbital contributions:
    U.h5    E     1       2.00000     -0.42540265E+04   1s
    U.h5    E     2       2.00000     -0.79454094E+03   2s
    U.h5    E     3       2.00000     -0.76701589E+03   2p-
    U.h5    E     4       2.00000     -0.62584540E+03   2p+
    U.h5    E     5       2.00000     -0.62584540E+03   2p+
    U.h5    E     6       2.00000     -0.20055934E+03   3s
    U.h5    E     7       2.00000     -0.18784209E+03   3p-
    U.h5    E     8       2.00000     -0.15532828E+03   3p+
    U.h5    E     9       2.00000     -0.15532828E+03   3p+
    U.h5    E    10       2.00000     -0.13488736E+03   3d-
    U.h5    E    11       2.00000     -0.13488736E+03   3d-
    U.h5    E    12       2.00000     -0.12826632E+03   3d+
    U.h5    E    13       2.00000     -0.12826632E+03   3d+
    U.h5    E    14       2.00000     -0.12826632E+03   3d+
    U.h5    E    15       2.00000     -0.51094916E+02   4s
    U.h5    E    16       2.00000     -0.45328009E+02   4p-
    U.h5    E    17       2.00000     -0.36852950E+02   4p+
    U.h5    E    18       2.00002     -0.36852950E+02   4p+
    U.h5    E    19       2.00000     -0.27589819E+02   4d-
    U.h5    E    20       2.00000     -0.27589819E+02   4d-
    U.h5    E    21       2.00000     -0.26025387E+02   4d+
    U.h5    E    22       2.00000     -0.26025387E+02   4d+
    U.h5    E    23       2.00000     -0.26025387E+02   4d+
    U.h5    E    24       1.99999     -0.13857922E+02   4f-
    U.h5    E    25       1.99999     -0.13857922E+02   4f-
    U.h5    E    26       1.99999     -0.13857922E+02   4f-
    U.h5    E    27       1.99999     -0.13449030E+02   4f+
    U.h5    E    28       1.99999     -0.13449030E+02   4f+
    U.h5    E    29       1.99999     -0.13449030E+02   4f+
    U.h5    E    30       1.99998     -0.13449030E+02   4f+
    U.h5    E    31       2.00009     -0.11340683E+02   5s 
    U.h5    E    32       2.00003     -0.90978058E+01   5p-
    U.h5    E    33       1.99993     -0.70954196E+01   5p+
    U.h5    E    34       1.99949     -0.70954196E+01   5p+
    U.h5    E    35       1.99993     -0.37973229E+01   5d-
    U.h5    E    36       1.99940     -0.37973229E+01   5d-
    U.h5    E    37       1.99996     -0.34992701E+01   5d+
    U.h5    E    38       1.99970     -0.34992701E+01   5d+
    U.h5    E    39       1.99898     -0.34992701E+01   5d+
    U.h5    E    40       1.98173     -0.17475069E+01   6s
    U.h5    E    41       1.93657     -0.10998688E+01   6p-
    U.h5    E    42       1.98095     -0.77331755E+00   6p+ (1/2)
    U.h5    E    43       1.75378     -0.77331755E+00   6p+ (3/2)
    U.h5    E    44       0.03763     -0.15556420E+00   7s
    U.h5    E    45       0.00001     -0.13518385E+00   5f- (1/2)
    U.h5    E    46       0.16298     -0.13518385E+00   5f- (3/2)
    U.h5    E    47       0.83970     -0.13518385E+00   5f- (5/2)
    U.h5    E    48       0.00001     -0.10440860E+00   5f+ (7/2)
    U.h5    E    49       0.00001     -0.10440860E+00   5f+ (5/2)
    U.h5    E    50       0.35514     -0.10440860E+00   5f+ (3/2)
    U.h5    E    51       0.95400     -0.10440860E+00   5f+ (1/2)
    U.h5    E    52       0.09263     -0.96703023E-01   6d- (1/2)
    U.h5    E    53       0.42441     -0.96703023E-01   6d- (3/2)
    U.h5    E    54       0.00004     -0.78515378E-01   6d+ (1/2)
    U.h5    E    55       0.31038     -0.78515378E-01   6d+ (3/2)
    U.h5    E    56       0.33421     -0.78515378E-01   6d+ (5/2)
    O1.h5   E     1       1.99960     -0.18913024E+02   1s
    O1.h5   E     2       1.89755     -0.87951960E+00   2s
    O1.h5   E     3       1.42472     -0.33082479E+00   2p-
    O1.h5   E     4       1.54866     -0.33082479E+00   2p+
    O1.h5   E     5       1.54866     -0.33082479E+00   2p+
    O2.h5   E     1       1.99960     -0.18913024E+02   1s
    O2.h5   E     2       1.89755     -0.87951960E+00   2s
    O2.h5   E     3       1.42472     -0.33082479E+00   2p-
    O2.h5   E     4       1.54866     -0.33082479E+00   2p+
    O2.h5   E     5       1.54866     -0.33082479E+00   2p+

We can to a large extent deduce the identity of the various atomic orbitals by just looking at orbital energies 
and I have added it by hand to the above output in a hopefully obvious notation. Even more precise information, 
such as :math:`m_j` value (indicated in parenthesis) can be obtained from looking at the Mulliken population analysis of each atom. 
From this list of gross population we extract the atomic valence configuration :math:`5f^{2.31}6d^{1.16}7s^{0.04}` for 
the uranium atom and :math:`2s^{1.91}2p^{4.28}` for oxygen. An interesting feature is that the uranium *6p* population 
adds up to *5.67* and not six electrons. This is a manifestation of the so-called *6p-hole*, due to overlap with the oxygen ligands, and mostly located to the :math:`6p_{3/2.1/2}` orbital.

It is also possible to carry out the projection analysis using the :ref:`PROJECTION_.ATOMS` keyword. We prepare the input file ``aproj.inp``

.. literalinclude:: aproj_UO2.inp

and prepare the checkpoint files::

   $ ln -s U_U.h5 U.h5
   $ ln -s O_O.h5 O.h5
   $ ln -s HF_UO2.h5 CHECKPOINT.h5

We run the projection analysis::

   pam --incmo --inp=aproj --mol=UO2 --copy="U.h5 O.h5"

The reader may convince himself that the output agrees with the projection analysis reported above.   

Projection analysis in localized orbitals
-----------------------------------------

It is instructive to carry out the projection analysis in *localized* orbitals; this allows us to "see" chemical bonds.
