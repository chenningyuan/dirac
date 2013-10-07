===================
Projection analysis
===================

------
Theory
------

A serious flaw of Mulliken population analysis is its basis set sensitivity. This is well known from
the non-relativistic domain. In the relativistic domain a further disadvantage is that the 
expansion of large and small components in separate basis setz means large and small component 
populations come separately and it is not easy to see how they are connected. Furthermore, if
the large and small component are expanded in scalar basis functions, as is the case of the 
DIRAC code, then the poulation analysis can not distinguish distributions from different 
spin-orbit components of atomic orbitals, for instance, :math:`p_{1/2}` and :math:`p_{3/2}`.

Projection analysis eliminates all these inconveniences. It is based on the simple concept of
molecular orbitals as linear combination of atomic orbitals and allows analysis of 
electronic structure in terms of well-defined atomic (or fragment) orbitals. 

Projection analysis requires as input a set of pre-calculated orbitals {:math:`\psi^A_p`} 
for the constituent atoms (or fragment) of the molecule. Here :math:`A` refers to 
the atomic center (or fragment) and :math:`p` is on orbital index for that center. 
The molecular orbitals are then expanded as

.. math::

   \left|\psi^{MO}_i\right> = \sum_A \sum_{p\in A} \psi^A_p c^A_{pi} + \left|\psi^{pol}_i\right>

Typically one will select the occupied orbitals of the constituent atoms. They are not 
guaranteed to fully span the molecular orbitals and therefore the expansion includes
their orthogonal complement :math:`\left|\psi^{pol}_i\right>`, denoted the polarization contribution.
The expansion coefficients :math:`c^A_{pi}` are found by solving the equation

.. math::

   \left<\psi^A_p|\psi^{MO}_i\right> = \sum_B\sum_{q\in B}\left<\psi^A_p|\psi^B_q\right>c^B_{qi}

obtained by projection.

------------------
Example 1: Methane
------------------

Introduction
------------

As a first simple example (albeit showing some subtle points) we shall consider the methane molecule.
The molecular input file ``CH4.mol`` is

.. literalinclude:: CH4.mol

Suppose that we have run a DFT(PBE) calculation using::

   pam --inp=PBE  --mol=CH4 --outcmo

where we keep the coefficients and the menu file ``PBE.inp`` is

.. literalinclude:: PBE.inp

Preparing for analysis
----------------------

We would now like to analyze the electronic structure of the methane molecule in terms of 
the orbitals of the constituent atoms. Now, looking in the output we find::

  Full group is:  T(d)           
  Represented as: D2 

This means that DIRAC correctly detected :math:`T_d` as the molecular point group symmetry. 
However, since the spatial symmetry of DIRAC is limited to :math:`D_{2h}` and subgroups the 
program will use the :math:`D_2` subgroup for the calculation. In this subgroup, as in the parent 
:math:`T_d` opint group, the four hydrogen  atoms are related by symmetry. Since DIRAC presently 
can only do projection analysis in terms of fragments complying to the molecular point group, 
we have to lower the symmetry if we want to carry out the projection analysis in terms of 
individual atoms. A simple solution is simply to have DIRAC dump the coeffients in :math:`C_1` 
symmetry. We then simply add::

  **GENERAL
  .ACMOUT

to ``PBE.inp`` and run the calculation as::

   pam --inp=PBE  --mol=CH4 --get "DFACMO=ac.CH4"

where the latter argument means that the file ``DFACMO`` of coefficients in :math:`C_1` (no) symmetry
is saved under the name ``ac.CH4``. It is also possible to convert an input ``DFCOEF`` file to corresponding ``DFACMO`` by using

.. literalinclude:: convert.inp

and the command::

  pam --incmo --inp=convert --mol=CH4 --get "DFACMO=ac.CH4"

However, before turning to the analysis itself there is one subtle point. Looking in the DIRAC
output again we find the lines::

    Original Coordinates                    
    --------------------                    
       6         0.00000000     0.00000000     0.00000000       1
       1         1.19031798     1.19031798     1.19031798       1
       1        -1.19031798    -1.19031798     1.19031798       1
       1         1.19031798    -1.19031798    -1.19031798       1
       1        -1.19031798     1.19031798    -1.19031798       1

    Symmetry class found: T(d)   

    Centered and Rotated                    
    --------------------                    
       1         1.19031798    -1.19031798     1.19031798       1
       1        -1.19031798     1.19031798     1.19031798       1
       1        -1.19031798    -1.19031798    -1.19031798       1
       1         1.19031798     1.19031798    -1.19031798       1
       6         0.00000000     0.00000000     0.00000000       1

showing that DIRAC will (translate and) rotate molecules during the process of symmetry detection.
For the projection analysis we can therefore not use the original coefficient file since we 
now will turn off symmetry detection and insist on running without symmetry. The simplest solution 
to this complication is to use the xyz-file that we also find in the output (we shall call it ``CH4.xyz``)

.. literalinclude:: CH4.xyz

(we have replaced ``H4`` by ``H``)
We now generate the atomic fragments. I recommend calculating the atoms in the highest possible 
symmetry and the converting the coefficients to no symmetry. This allows a very precise identification
of the individual atomic orbitals (including :math:`m_j` quantum number). We accordingly run::

   pam --inp=H  --mol=H --get "DFACMO=ac.H"
   pam --inp=C  --mol=C --get "DFACMO=ac.C"

with corresponding input files ``H.inp``, 

.. literalinclude:: H.inp 

``H.mol``, 

.. literalinclude:: H.mol

``C.inp``

.. literalinclude:: H.inp 

and ``C.mol``

.. literalinclude:: C.mol

Running the analysis
--------------------

We are now ready for the projection analysis itself. We prepare the input file ``prj.inp``

.. literalinclude:: prj.inp

In the ``*PROJEC`` section we specify five fragments identified by the name of their coefficient files.
 We also specify ``.OWNBAS``, which means that the fragments have been calculated in their own basis and not the full molecular one. Note also that since we will use a *xyz* - file for this calculation
we have added a ``**MOLECULE`` section where we specify basis *and* we explicitly turn off symmetry.

As a final preparation we make symbolic links to the needed coefficient files::

   $ ln -s ac.C AFCXXX
   $ ln -s ac.H AFH1XX
   $ ln -s ac.H AFH2XX
   $ ln -s ac.H AFH3XX
   $ ln -s ac.H AFH4XX
   $ ln -sf ac.CH4 DFCOEF

We now run the projection analysis using::

   pam --inp=prj --mol=CH4.xyz --incmo --copy="AF*"

Interpreting the output
-----------------------

The first thing to look for in the output is the section::
   
   * Total gross contributions:
      AFCXXX               6.4301   E -     6.4301   P -     0.0000
      AFH1XX               0.8716   E -     0.8716   P -     0.0000
      AFH2XX               0.8716   E -     0.8716   P -     0.0000
      AFH3XX               0.8716   E -     0.8716   P -     0.0000
      AFH4XX               0.8716   E -     0.8716   P -     0.0000
      Polarization:     0.0834

It shows the gross population (in Mulliken sense, but in terms of the fragment orbitals) of each atom 
and then, most importantly, the gross population of the polarization contribution. In this example it
is small and perfectly acceptable. If it becomes too large the analysis looses meaning because it 
impies that significant electron density has not been assigned to any fragment. In the case of
significant polarization contribution one can either increase the number of atomic orbitals used
in the analysis or turn on repolarization. 

Having concluded that the polarization contribution is acceptable, we can proceed with the analysis.
From the atomic gross population we can already extract atomic charges (keeping in mind that some 
density is not assigned to any atom). We can see that the carbon atom has a small negative charge of
-0.43 atomic units, whereas the hydrogens have positive charges  +0.13. 

We next turn to this section::

      * Total reference orbital contributions:
          AFCXXX  E     1       1.99999     -0.10660951E+02
          AFCXXX  E     2       1.27874     -0.97828469E+00
          AFCXXX  E     3       1.05077     -0.65580701E+00
          AFCXXX  E     4       1.05030     -0.65541446E+00
          AFCXXX  E     5       1.05030     -0.65541446E+00
          AFH1XX  E     1       0.87162     -0.23118508E+00
          AFH2XX  E     1       0.87162     -0.23118508E+00
          AFH3XX  E     1       0.87162     -0.23118508E+00
          AFH4XX  E     1       0.87162     -0.23118508E+00
   
which gives the accumulated gross population of each fragment orbital. We can identify 
the fragments from the name of the coefficient file and its orbital energy. For instance
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

   * Electronic eigenvalue nr.  3: -0.3411330794652       (Occupation : f = 1.0000)
   ================================================
       Orbital            Total           Eigenvalue       Kramers partner 1   Kramers partner 2
       AFCXXX  E     3       0.57945     -0.65580701E+00   (-1.0000,-0.0000)   ( 0.0000,-0.0000)
       AFH1XX  E     1       0.31559     -0.23118508E+00   ( 0.0000, 0.5774)   (-0.5774,-0.5774)
       AFH2XX  E     1       0.31559     -0.23118508E+00   ( 0.0000,-0.5774)   ( 0.5774,-0.5774)
       AFH3XX  E     1       0.31559     -0.23118508E+00   ( 0.0000,-0.5774)   (-0.5774, 0.5774)
       AFH4XX  E     1       0.31559     -0.23118508E+00   ( 0.0000, 0.5774)   ( 0.5774, 0.5774)
   * Gross contributions:
      AFCXXX               0.5254   E -     0.5254   P -     0.0000
      AFH1XX               0.1165   E -     0.1165   P -     0.0000
      AFH2XX               0.1165   E -     0.1165   P -     0.0000
      AFH3XX               0.1165   E -     0.1165   P -     0.0000
      AFH4XX               0.1165   E -     0.1165   P -     0.0000
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

-----------------
Example 2: Uranyl
-----------------

We next turn to a molecule with a heavy atom and a more complicated electronic structure, namely uranyl 
:math:`UO_2^{2+}`. Here the two oxygens are related by symmetry, and so one would at first consider 
running the projection analysis with no symmetry, following the procedure outlined in the methane 
example above. However, we can reduce the symmetry from :math:`D_{\infty h}` to :math:`C_{\infty v}` 
by introducing a ghost center. Our molecular input file then looks like

.. literalinclude:: UO2.mol

Now the oxygen atoms are no longer connected by symmetry and we can carry out the projection 
analysis in :math:`C_{\infty v}` symmetry. 
