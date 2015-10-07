:orphan:

.. _GRASP_computer_exercises:

************************
GRASP computer exercises
************************

The following exercises are to be done in groups of 2-3 students. Each group shall cover different groups/atoms of the periodic table of elements.
At the end of the day, all groups should give a short presentation on their findings to the other groups.

In order to run a GRASP calculation you need an input file of arbitrary name, e.g. :file:`grasp.inp` which will be setup in the
following. The run-command inside your cluster run script looks as::

  $ grasp.x < grasp.inp > grasp.out

where the GRASP output is re-directed to the file :file:`grasp.out`. Make sure that you give a proper name to your output
file if you want to take a look at what you did back home.

Unless otherwise stated we will always use a Gaussian nuclear model in GRASP.

Exercises - Day 1
=================

Nuclear models in GRASP
-----------------------

The electron-nucleus interaction in the electronic Hamiltonian is :math:`V_{en}=-e\phi_{nuc}` 
where :math:`\phi_{nuc}`  is the nuclear scalar potential

.. math::

   \phi_{nuc}(\mathbf{r}_1) = \frac{1}{4\pi\varepsilon_0}\int\frac{\rho_{nuc}(\mathbf{r}_2)}{\left|\mathbf{r}_1-\mathbf{r}_2\right|}d^3\mathbf{r}_2

GRASP offers four different models for the nuclear charge distribution :math:`\rho_{nuc}` : the point charge, the homogeneously charged sphere, 
the Fermi 2-parameter charge distribution and the Gaussian charge distribution.
You find more information about these models `at this URL <http://www.few.vu.nl/~visscher/FiniteNuclei/FiniteNuclei.htm>`_ .
(Note that there is a systematic error in this source in that the potentials come with a negative sign whereas this comes from multiplication of 
the electron charge to give the electron-nucleus interaction).
We will start with calculations on a one-electron system for each of the four nuclear models.

::

  One electron atom
  1 1 1
  1s 1
  ANG 7
  -1
  MCP 8 9
  MCDF
  9 0 10
  86 222.0
  NUC 0
  AL
  END

* point nucleus::

   NUC 0

* finite nucleus (Uniform Sphere)::

   NUC 1

* finite nucleus (Fermi)::

   NUC 2

* finite nucleus (Gaussian)::

   NUC 3

**Questions**:

1. How does the use of a finite-nucleus model affect the total energy of an atom compared to the point-nucleus model?

2. What is the effect on the spin-orbit splitting?

3. Do differences between finite-nuclear models increase or decrease with increasing nuclear charge *Z*?

Exercises - Day 2
=================

.. note:: A database for atomic data can be found `here <http://physics.nist.gov/PhysRefData/Handbook/periodictable.htm>`_ .

Atomic numerical calculations
-----------------------------

In order to get started on "real" atoms we calculate the average level structure of two possible ground state configurations of coinage-metals (Cu, Ag, Au, Rg), using GRASP.


::

  Cu average calculation
  2 7 1
  1s
  2s
  2p
  3s
  3p
  3d 10 9
  4s  1 2
  ANG 7
  -1
  -1
  MCP 8 9
  MCDF
  9 0 10 
  29 63.0
  NUC 3
  SCF NITIT=10
  AL
  END

Run all atoms of the coinage-metal series::

   Cu: Z = 29,  A =  63
   Ag: Z = 47,  A = 107
   Au: Z = 79,  A = 197
   Rg: Z = 111, A = 283

In order to generate a less biased average set of orbital solutions (degeneracy of the d :math:`^9` s :math:`^2`: 10; d :math:`^{10}` s :math:`^1`: 2), compare your results with an additional calculation using:

::

  AL 10 2

.. note:: We explicitly set the maximum number of SCF iterations ("SCF NITIT=10") in order to allow Cu to converge (default SCF iteration = 6)


**Questions**:

1. The two candidates for the ground state configuration are :math:`(n-1)d^{10}ns^{1}` and :math:`(n-1)d^{9}ns^{2}`. For each configuration what are the possible 
values of total angular momentum :math:`J` ? Your answer will help you identify the ground state configuration by looking at the list of energies of the individual 
atomic levels. 

2. Does the ground state configuration vary down the coinage-metal series?

3. How does the d-s excitation energy (:math:`(n-1)d^{10}ns^{1}\rightarrow (n-1)d^{9}ns^{2}`) change down the coinage-metal series?

Spin-orbit mixing within a group
--------------------------------

In this exercise we shall investigate how spin-orbit mixing evolves as we descend a column of the periodic table. 
We will consider group 14, that is, the carbon group. A sample input, for silicium is given below

::

  jj-LS transformation Si atom
  1 5 2
  1s
  2s
  2p
  3s
  3p 2
  ANG
  -1
  MCP 8 9
  MCDF
  9 0 10
  14 28.0
  NUC 3
  AL
  END

The :math:`3p^2` electron configuration of silicium gives rise the three LS-states: :math:`^3P`, :math:`^1D` and :math:`^1S`. 
With the introduction of spin-orbit coupling they are split into :math:`^3P_2`, :math:`^3P_1`, :math:`^3P_0`, :math:`^1D_2` and :math:`^1S_0` (verify this !) 
and states with the same :math:`J` may interact and mix. This means for instance that :math:`^3P_2` and :math:`^1D_2` interact, leading to singlet-triplet mixing, 
forbidden in a non-relativistic framework. There is also possibility of mixing between :math:`^3P_0` and :math:`^1S_0`.

In order to understand how to proceed, let us first look at what an average level (AL) calculation means: 

In the GRASP input we specify one or more electron configurations. Assume that we generate all possible Slater 
determinants :math:`\left\{\Phi_I\right\}` corresponding to the configurations, and then build and diagonalize the Configuration 
Interaction (CI) matrix in the space of the Slater determinants

.. math::

   H\mathbf{c}_K=E_K\mathbf{c}_K;\quad H_{IJ} = \left<\Phi_I\left|H\right|\Phi_J\right>

This gives us a set of electronic states expressed in terms of atomic state functions 

.. math::

   \Psi_K = \sum_I^{N}\Phi_I c_{IK}

whose number is determined by the degeneracy of the state and where :math:`N` refers to the total number of Slater determinants 
and therefore number of microstates/atomic state functions.

In order to build the CI matrix we need orbitals. Rather than optimizing them for an individual electronic state we shall 
determine them by optimizing the *average energy*

.. math::

   E^{av} = \frac{1}{N}\sum_K^N \left<\Psi_K\left|H\right|\Psi_K\right>

This looks rather impossible because the average energy contains the atomic state function that we do not know until we have diagonalized the CI 
matrix, which we can not build before we have orbitals available. However, if we expand the atomic state functions in terms of determinants 
we obtain

.. math::

   E^{av} = \frac{1}{N}\sum_K^N \sum_I^N\sum_J^N \left<\Phi_I\left|H\right|\Phi_J\right>c^*_{IK}c_{JK} = \frac{1}{N}\sum_K^N \left<\Phi_K\left|H\right|\Phi_K\right>

We arrive at the final result bwcause the atomic state functions are orthonormal such that the expansion coefficients :math:`\left\{c_{IK}\right\}` form 
a unitary matrix. Now the average energy is expressed in terms of our set of Slater determinants and we may proceed. In fact, GRASP will first optimize the
orbitals by an SCF procedure based on the above expression and then solve the CI problem to get the individual atomic state functions and atomic level energies.

The atomic state functions have well-defined total angular momentum :math:`J`. In the "ANG" input section of GRASP we specify the desired :math:`J` for 
each electron configuration. If we want all possible values, we simply set "-1" as in the present case. A single Slater determinant is in general not an 
eigenfunction of :math:`J^2`. GRASP therefore constructs configuration functions (CSFs), that is, linear combinations of Slater determinants with well-defined 
:math:`J` and parity. These are by default listed in the output, but in the exercises above we had set "ANG 7" to suppress this output. 
With the above input we find for instance that the first CSF is generated from the valence configuration :math:`3p_{1/2}^13p_{2/2}^1` and has :math:`J=2` and even parity.
At the end of the GRASP output the composition of atomic state functions in terms of CSFs is listed. Since the CSFs are orthonormal the weight CSF of index :math:`I` 
to an atomic state function of index :math:`K` is simply :math:`\left|c_{IK}\right|^2`.

In order to investigate singlet-triplet mixing we also need to know the composition of atomic state functions in terms of LS-coupled CSFs. 
The third number "2" in line 2 in the above input makes GRASP also set up CSFs in LS-coupling. The non-relativistic CSFs are listed after the relativistic ones. 
With the present input we find that the first non-relativistic CSF couples to :math:`^3P_2`. The weight of non-relativistic CSFs to atomic state functions is listed 
after the relativistic ones.

**Questions**:

1. Does the singlet-triplet mixing (*state interaction*) increase or decrease with increasing *Z*? You can investigate this by monitoring the mixing of :math:`^3P_2` and :math:`^1D_2` 
in the lowest atomic level of :math:`J=2`, making a graph showing how the weight of :math:`^1D_2` evolves with increasing *Z*.

2. How much of the ground state is described by the :math:`np_{1/2}^{2}` configuration with increasing *Z*?  You can investigate this by monitoring the composition of the ground-state 
in terms of relativistic CSFs. You can calculate the effective occupation number of the :math:`np_{1/2}` orbital by, for each CSF, multiplying its occupation by the weight of the CSF, 
then summing up these numbers. You can then plot the effective occupation number of the :math:`np_{1/2}` orbital as a function of increasing *Z*. How do you interpret your result ?

Speed of light and the non-relativistic limit
---------------------------------------------

* Run an average-of-configuration SCF calculation the gold atom:

::

  Au average calculation
  2 14 1
  1s
  2s
  2p
  3s
  3p
  3d
  4s
  4p
  4d
  5s
  5p
  4f
  5d 10 9
  6s 1  2
  ANG 7
  -1
  -1
  MCP 8 9
  MCDF
  9 0 10
  79 197.0
  NUC 3
  SCF CON=1.0D0
  AL
  END


* GRASP allows you to vary the speed of light *c* between *Z* < *c* < :math:`\infty`::

   SCF CON=x.xDx

* perform a non-relativistic calculation using a very large *c* value, *c* :math:`>> 137`::

   SCF CON=1.0D6

**Questions**:

1. Is it possible to set *c < Z*?

2. Can you argue from a comparison of the relativistic ("regular" speed-of-light)  and non-relativistic calculation why gold is *golden*? What about Ag?

3. Which spinors are stabilized by relativistic effects and which are destabilized?


Radial densities for the d-/p block
-----------------------------------

* Plot radial density for heavy elements, compare *s* and *d* functions (*d* block) resp. p :math:`_{1/2}` and p :math:`_{3/2}` (p block). 

A simple FORTRAN code :download:`orbpri.F` for doing it is available. Here we shall demonstrate how to use it for the copper atom. 
First we remind you that the form of relativistic atomic orbitals is

.. math::

   \left[\begin{array}{c}R^L\chi_{\kappa,m}\\iR^S\chi_{-\kappa,m}\\\end{array}\right] = \frac{1}{r}\left[\begin{array}{c}P\chi_{\kappa,m}\\iQ\chi_{-\kappa,m}\\\end{array}\right]

so the radial density is

.. math::

   n(r)=\left(\left\{R^L\right)^2+\left(R^S\right)^2\right\}r^2=\left(P^2+Q^2\right).

The input file is

::

  Cu average calculation
  2 7 1
  1s
  2s
  2p
  3s
  3p
  3d 10 9
  4s  1 2
  ANG 7
  -1
  -1
  MCP 8 9
  MCDF 28
  9 0 10 0 11
  29 63.0
  NUC 3
  ORBOUT 1 - 7
  AL
  END

The line ORBOUT specifies that the radial functions :math:`P` and :math:`Q` (and other information) is printed to file. This file will have unit number 11, as specified by the final number following the MCDF line.

We start up orbpri.x::

   $orbpri.x
   *** OUTPUT from ORBPRI ***
   Give unit number of orbital file to read

and we enter unit number 11::

   11
   *** OUTPUT from ORBPRI ***
   * Atomic number:               29.0
   * Number of grid points:       2080
   * First grid point:       0.345E-07
   * Grid step size :        0.100E-01
      
   * Number of orbitals read:       10
      
   Do you want info on orbitals (y/n) ?

Well, why not ?::

   y
   Orb.        E                   P0                     Q0
    1  1s       -332.845231246527    353.946570713740      0.000000000000
    2  2s        -41.883350988707    110.925718644588      0.000000000000
    3  2p-       -36.527466295898      0.975314932490      9.217465299856
    4  2p        -35.763812849690    769.447328175348      0.000000000000
    5  3s         -5.353239297458     41.725471888184      0.000000000000
    6  3p-        -3.614445494458      0.362775735823      3.428505650010
    7  3p         -3.511783710417    286.145123083999      0.000000000000
    8  3d-        -0.674391837066      0.265542615160      5.019157934192
    9  3d         -0.661421437065    264.801381125017      0.000000000000
   10  4s         -0.289259556935      9.151716417086      0.000000000000
   Select one of the following:
     1. Print radial functions P and Q.
     2. Print radial functions R(large) and R(small).
     3. Print radial distribution of individual orbitals.
     4. Exit.

We choose to look at the density of individual orbitals::

   3
   Name atom (A3)

The rest is then rather straightforward::

   Cu
   Outer radial point = :   36.862
   Give index of orbital to print(-1 to quit):
   6
   Give index of orbital to print(-1 to quit):
   7
   Give index of orbital to print(-1 to quit):
   8
   Give index of orbital to print(-1 to quit):
   9
   Give index of orbital to print(-1 to quit):
   10
   Give index of orbital to print(-1 to quit):
   -1
   *** EXITING ORBPRI ***

Looking in your working directory you now see a series of files with self-explanatory titles.
With your favorite plotting program you may now make a plot like this one

.. image:: Cu_raddens.png
   :scale: 80
   :alt: Copper radial density


Each group should pick one of their favorite block/group.

::

   calculate the 3d, 4d, 5d, 6d block (coinage metals)
   calculate the 3d, 4d, 5d, 6d block (group: "Zn, Cd, Hg, Cn")
   calculate the 2p, 3p, 4p, 5p, 6p, 7p block (group 14)
   calculate the 2p, 3p, 4p, 5p, 6p, 7p block (group 16)

**Questions**:

1. How does the radial overlap between *s* and *d* resp. p :math:`_{1/2}` and p :math:`_{3/2}` change down a group?

2. Where would you expect the most open-shell character (at the top or bottom of a group)?


Using a non-relativistic or relativistically optimized basis set for heavy atoms?
---------------------------------------------------------------------------------

.. note:: all basis set inputs are available online: http://wiki.chem.vu.nl/dirac/index.php/GRASP_and_DIRAC_inputs

* run the Thallium atom with a non-relativistic double-:math:`\zeta`  basis set:

::

  Thallium atom
   1 15 1
  1S
  2S
  2P
  3S
  3P
  3D
  4S
  4P
  4D
  4F
  5S
  5P
  5D
  6S
  6P 1
  ANG 1 7
  -1
  MCP 1 2
  MCDF 24
  2 0 3 0 14
  81 204.0
  NUC 3
  GAUSS 4
  1S  22
  2P- 17
  2P  17 2P- 1 17
  3D- 13
  3D  13 3D- 1 13
  4F-  8
  4F   8 4F- 1  8
  1.8623969E+07
  3.6867955E+06
  9.6378590E+05
  2.9791376E+05
  1.0338698E+05
  3.9115255E+04
  1.5831221E+04
  6.7653929E+03
  3.0236823E+03
  1.4021816E+03
  6.6607727E+02
  3.0047853E+02
  1.5516535E+02
  7.9776475E+01
  3.7180980E+01
  1.9582277E+01
  8.6490359E+00
  4.4732684E+00
  1.5727172E+00
  7.2742584E-01
  1.4401137E-01
  5.4696933E-02
  8.4903516E+04
  2.0108565E+04
  6.5288090E+03
  2.4952142E+03
  1.0578546E+03
  4.8183325E+02
  2.3099156E+02
  1.1460187E+02
  5.6345780E+01
  2.8858154E+01
  1.4257420E+01
  7.1249109E+00
  3.4904261E+00
  1.5777429E+00
  6.6465131E-01
  1.5479093E-01
  4.8201891E-02
  4.4351057E+03
  1.3418293E+03
  5.2396706E+02
  2.3190630E+02
  1.1060211E+02
  5.4985007E+01
  2.8038853E+01
  1.4122437E+01
  7.0047966E+00
  3.4026123E+00
  1.5589903E+00
  6.6349200E-01
  2.5380659E-01
  4.6883901E+02
  1.6331607E+02
  6.9830068E+01
  3.2496980E+01
  1.5735686E+01
  7.6001533E+00
  3.5619156E+00
  1.4998623E+00
  ORBOUT 1 - 15
  AL
  END

* run the Thallium atom with a relativistic double-:math:`\zeta`  basis set:

::

  Thallium atom
   1 15 1
  1S
  2S
  2P
  3S
  3P
  3D
  4S
  4P
  4D
  4F
  5S
  5P
  5D
  6S
  6P 1
  ANG 1 7
  -1
  MCP 1 2
  MCDF 24
  2 0 3 0 14
  81 204.0
  NUC 3
  GAUSS 4
  1S  24
  2P- 20
  2P  20 2P- 1 20
  3D- 13
  3D  13 3D- 1 13
  4F-  8
  4F   8 4F- 1  8
  5.8033828E+07
  1.4702562E+07
  4.6355575E+06
  1.5852740E+06
  5.7924626E+05
  2.2187944E+05
  8.8730780E+04
  3.6848078E+04
  1.5843298E+04
  7.0310514E+03
  3.2121888E+03
  1.5070783E+03
  7.2301062E+02
  3.3980515E+02
  1.7093901E+02
  8.8017323E+01
  4.3016225E+01
  2.2609863E+01
  1.0266569E+01
  5.2648225E+00
  1.8988944E+00
  8.8461372E-01
  2.0106574E-01
  7.3237091E-02
  1.3728265E+07
  2.3572128E+06
  5.2356678E+05
  1.3629076E+05
  4.0148089E+04
  1.3182021E+04
  4.7879862E+03
  1.9028441E+03
  8.1502288E+02
  3.6967602E+02
  1.7504377E+02
  8.4456267E+01
  4.2043229E+01
  2.1172075E+01
  1.0107599E+01
  4.9031493E+00
  1.9632888E+00
  8.0549699E-01
  2.0049629E-01
  5.5927884E-02
  7.8021591E+03
  2.0396979E+03
  7.3049097E+02
  3.0516736E+02
  1.3964975E+02
  6.7292643E+01
  3.3520599E+01
  1.6568385E+01
  8.0815299E+00
  3.8806549E+00
  1.6956350E+00
  6.9524060E-01
  2.5357650E-01
  4.8999559E+02
  1.6633197E+02
  6.9980085E+01
  3.2175090E+01
  1.5417577E+01
  7.3699578E+00
  3.4086613E+00
  1.4138469E+00
  ORBOUT 1 - 15
  AL
  END

**Questions**:

1. How does the spin-orbit splitting from both basis set calculation compare with the experimental reference?

2. Can you recommend to use a non-relativistically optimized basis set in relativistic (four-component) calculation?  Why/why not?

3. What are the major qualitative differences between both basis sets?


Basis set quality
-----------------

.. note:: all basis set inputs are available online: http://wiki.chem.vu.nl/dirac/index.php/GRASP_and_DIRAC_inputs

* run the Bi atom numerically (if not already done in `Radial densities for the d-/p block`_):

::

  Bismuth atom
  1 15 1
  1S
  2S
  2P
  3S
  3P
  3D
  4S
  4P
  4D
  4F
  5S
  5P
  5D
  6S
  6P 3
  ANG 7
  -1
  MCP 8 9
  MCDF 28
  9 0 10
  83 209.0
  NUC 3
  AL
  END

* repeat the calculation using Ken Dyall's

 i) double-:math:`\zeta` basis set

::

  Bi atom - DZ basis
  1 15 1
  1s
  2s
  2p
  3s
  3p
  3d
  4s
  4p
  4d
  4f
  5s
  5p
  5d
  6s
  6p 3
  ANG 1 7
  -1
  MCP 1 2
  MCDF 24
  2 0 3 0 14
  83 209.0
  NUC 3
  GAUSS 4
  1S  24
  2P- 20
  2P  20 2P- 1 20
  3D- 14
  3D  14 3D- 1 14
  4F-  8
  4F   8 4F- 1  8
   5.8022684E+07
   1.4819799E+07
   4.7223874E+06
   1.6321833E+06
   6.0232740E+05
   2.3262600E+05
   9.3657019E+04
   3.9098920E+04
   1.6878840E+04
   7.5127674E+03
   3.4393803E+03
   1.6158149E+03
   7.7585154E+02
   3.6612190E+02
   1.8440048E+02
   9.5065117E+01
   4.6466268E+01
   2.4477290E+01
   1.1190161E+01
   5.7978488E+00
   2.1587077E+00
   1.0422614E+00
   2.6590866E-01
   1.0073440E-01
   1.5267851E+07
   2.6863694E+06
   6.0456060E+05
   1.5848335E+05
   4.6785870E+04
   1.5325077E+04
   5.5345524E+03
   2.1842431E+03
   9.2942393E+02
   4.1932007E+02
   1.9785121E+02
   9.5696786E+01
   4.7435630E+01
   2.3929631E+01
   1.1351604E+01
   5.5346123E+00
   2.2410872E+00
   9.5317500E-01
   2.4829310E-01
   7.7915822E-02
   1.1456342E+04
   2.9347912E+03
   1.0364305E+03
   4.3032355E+02
   1.9654041E+02
   9.5235753E+01
   4.7867448E+01
   2.4364758E+01
   1.2223922E+01
   6.0607197E+00
   2.8412869E+00
   1.2973356E+00
   5.5200808E-01
   2.0223207E-01
   5.5102305E+02
   1.8692436E+02
   7.8801541E+01
   3.6351986E+01
   1.7527887E+01
   8.4583618E+00
   3.9675228E+00
   1.6865442E+00
  AL
  END


ii) triple-:math:`\zeta` basis set


::

 Bi atom - TZ basis
  1 15 1
  1s
  2s
  2p
  3s
  3p
  3d
  4s
  4p
  4d
  4f
  5s
  5p
  5d
  6s
  6p 3
  ANG 1 7
  -1
  MCP 1 2
  MCDF 24
  2 0 3 0 14
  83 209.0
  NUC 3
  GAUSS 4
  1S  30
  2P- 25
  2P  25 2P- 1 25
  3D- 17
  3D  17 3D- 1 17
  4F- 11
  4F  11 4F- 1 11
   6.0311989E+07
   1.6051490E+07
   5.4946022E+06
   2.0921783E+06
   8.7228608E+05
   3.8455448E+05
   1.7769146E+05
   8.4842022E+04
   4.1639260E+04
   2.0888520E+04
   1.0685288E+04
   5.5616491E+03
   2.9418744E+03
   1.5792819E+03
   8.5880021E+02
   4.7160384E+02
   2.6367206E+02
   1.4960897E+02
   8.5864339E+01
   4.9330236E+01
   2.9084112E+01
   1.7004392E+01
   9.4518823E+00
   5.3907940E+00
   2.9254412E+00
   1.5767178E+00
   8.2262394E-01
   3.3654907E-01
   1.5552971E-01
   6.7994583E-02
   4.5716312E+07
   1.2020013E+07
   3.5338508E+06
   1.1328857E+06
   3.8880908E+05
   1.4139323E+05
   5.4175877E+04
   2.1825254E+04
   9.2485101E+03
   4.1227683E+03
   1.9285313E+03
   9.4167499E+02
   4.7644513E+02
   2.4829591E+02
   1.3204380E+02
   7.0732928E+01
   3.8843957E+01
   2.1268489E+01
   1.1356211E+01
   6.0893913E+00
   3.1150062E+00
   1.5809224E+00
   7.6275107E-01
   2.9950363E-01
   1.2171302E-01
   3.6056471E+04
   8.7942690E+03
   2.9773896E+03
   1.2026340E+03
   5.4320948E+02
   2.6392698E+02
   1.3504708E+02
   7.1319062E+01
   3.8578165E+01
   2.0919821E+01
   1.1151610E+01
   5.9122803E+00
   3.0430319E+00
   1.5135542E+00
   7.2243078E-01
   3.2503953E-01
   1.2657989E-01
   1.0703029E+03
   3.6323880E+02
   1.5446734E+02
   7.3368355E+01
   3.6869724E+01
   1.9107106E+01
   9.8894189E+00
   5.0068176E+00
   2.3615820E+00
   8.9551709E-01
   2.7862580E-01
  AL
  END


iii) quadruple-:math:`\zeta` basis set


::

  Bi atom - QZ basis
  1 15 1
  1s
  2s
  2p
  3s
  3p
  3d
  4s
  4p
  4d
  4f
  5s
  5p
  5d
  6s
  6p 3
  ANG 1 7
  -1
  MCP 1 2
  MCDF 24
  2 0 3 0 14
  83 209.0
  NUC 3
  GAUSS 4
  1S  34
  2P- 31
  2P  31 2P- 1 31
  3D- 21
  3D  21 3D- 1 21
  4F- 14
  4F  14 4F- 1 14
   6.0327973E+07
   1.6059400E+07
   5.5002311E+06
   2.0965944E+06
   8.7614365E+05
   3.8775710E+05
   1.8029442E+05
   8.6844149E+04
   4.3076212E+04
   2.1842909E+04
   1.1293409E+04
   5.9507204E+03
   3.2014852E+03
   1.7641591E+03
   1.0016015E+03
   5.9564241E+02
   3.9672091E+02
   2.5845574E+02
   1.6068287E+02
   9.8444637E+01
   6.0146010E+01
   3.7934344E+01
   2.3912645E+01
   1.5276004E+01
   9.2087298E+00
   5.5574031E+00
   3.2023544E+00
   1.9180785E+00
   1.1731977E+00
   7.4463019E-01
   4.7115520E-01
   2.3902743E-01
   1.1886460E-01
   5.6025919E-02
   7.2548229E+07
   2.5127961E+07
   8.9943363E+06
   3.3978472E+06
   1.3406152E+06
   5.4867118E+05
   2.3165676E+05
   1.0069143E+05
   4.5050536E+04
   2.0784822E+04
   9.9000264E+03
   4.8664210E+03
   2.4655886E+03
   1.2855834E+03
   6.8790131E+02
   3.7620602E+02
   2.0950837E+02
   1.1794969E+02
   6.6689215E+01
   3.8460927E+01
   2.2109515E+01
   1.2482731E+01
   7.1083213E+00
   4.0003720E+00
   2.2284663E+00
   1.2125620E+00
   6.4365302E-01
   3.1488676E-01
   1.5158717E-01
   7.1607861E-02
   3.2913977E-02
   1.6463954E+05
   3.8667080E+04
   1.2458903E+04
   4.7981120E+03
   2.0858623E+03
   9.8974360E+02
   5.0061517E+02
   2.6520654E+02
   1.4561957E+02
   8.1829425E+01
   4.6861624E+01
   2.7006142E+01
   1.5460270E+01
   8.8034943E+00
   4.9821735E+00
   2.7644317E+00
   1.4934857E+00
   8.0288750E-01
   4.4714894E-01
   2.4637901E-01
   1.0017191E-01
   2.5971422E+03
   8.6492152E+02
   3.6744541E+02
   1.7603323E+02
   9.0982272E+01
   4.8997622E+01
   2.7207540E+01
   1.5262097E+01
   8.5485499E+00
   4.7315946E+00
   2.5336398E+00
   1.2638986E+00
   4.4762345E-01
   1.7508119E-01
  AL
  END

**Questions**:

1. Look at the ground and lowest excited states - does a double-:math:`\zeta`  basis set suffice to reach a good agreement with the numerical solution?

2. Compare with experiment!

3. R.95 values as a measure for the density at the nucleus: would a quadruple-:math:`\zeta`  basis set be large enough to reach agreement with the numerical reference value? What may be missing?

Using *j*- or *l*-optimized basis sets?
---------------------------------------

.. note:: all basis set inputs are available online: http://wiki.chem.vu.nl/dirac/index.php/GRASP_and_DIRAC_inputs

* use a *j*-optimized basis set to calculate the spin-orbit spliiting of the ground state of element 113.

::

  element 113 - j-optimized DZ basis
  1 19 1
  1s
  2s
  2p
  3s
  3p
  3d
  4s
  4p
  4d
  4f
  5s
  5p
  5d
  5f
  6s
  6p
  6d
  7s 2
  7p 1
  ANG 7
  -1
  MCP 8 9
  MCDF 11
  9 0 0 0 10
  113 287.0
  NUC 3
  SCF NITIT=6 CON=1.0E0
  GAUSS 4 CONTROL
  CONTROL CNVVEC 6.3 CNVENG 12
  1s  26
  2p- 26 1s 1 26
  2p  22
  3d- 22 2p 1 22
  3d  15
  4f- 15 3d 1 15
  4f  10
  5.410555451E+07   1.477674220E+07   5.046760516E+06   1.881036246E+06
  7.519670234E+05   3.121717712E+05   1.337537787E+05   5.861701238E+04
  2.623376918E+04   1.196778403E+04   5.565613895E+03   2.641110608E+03
  1.279308333E+03   6.350605841E+02   3.229251397E+02   1.668148069E+02
  8.855127743E+01   4.626843234E+01   2.491302232E+01   1.332236224E+01
  7.050429402E+00   3.566653213E+00   1.746480412E+00   8.012456416E-01
  2.363147173E-01   8.211626183E-02
  1.832973111E+06   3.677410004E+05   1.042558311E+05   3.551995083E+04
  1.365628428E+04   5.729961829E+03   2.571176326E+03   1.216215930E+03
  6.003294614E+02   3.061461712E+02   1.601715602E+02   8.574943752E+01
  4.600764505E+01   2.487057706E+01   1.355724873E+01   7.189512031E+00
  3.709687491E+00   1.753488637E+00   8.174464605E-01   3.431291818E-01
  1.196202297E-01   3.605428667E-02
  1.438935404E+04   4.111888007E+03   1.542160218E+03   6.637471705E+02
  3.111800230E+02   1.538758849E+02   7.883962728E+01   4.095259646E+01
  2.148494432E+01   1.119186568E+01   5.692583697E+00   2.791037755E+00
  1.265293623E+00   5.247952224E-01   1.904499576E-01
  1.250294088E+03   4.300054627E+02   1.834372720E+02   8.629363698E+01
  4.258654175E+01   2.151568823E+01   1.081918963E+01   5.287539266E+00
  2.459405736E+00   1.015550465E+00
  AL
  END

* use an *l*-optimized double-:math:`\zeta`  basis set to calculate the spin-orbit spliiting of the ground state of element 113.

::

  element 113 - l-optimized DZ basis
  1 19 1
  1s
  2s
  2p
  3s
  3p
  3d
  4s
  4p
  4d
  4f
  5s
  5p
  5d
  5f
  6s
  6p
  6d
  7s 2
  7p 1
  ANG 7
  -1
  MCP 8 9
  MCDF 11
  9 0 0 0 10
  113 287.0
  NUC 3
  SCF NITIT=6 CON=1.0E0
  GAUSS 4 CONTROL FILE
  FILE ORBOUT 15
  CONTROL CNVVEC 6.2
  1s  26
  2p- 23 INPUT 23
  2p  23 2p- 1 23
  3d- 16 INPUT 16
  3d  16 3d- 1 16
  4f- 10 INPUT 10
  4f  10 4f- 1 10
  5.265642072E+07   1.398668173E+07   4.744225153E+06   1.773651781E+06
  7.211584750E+05   3.072954606E+05   1.361714165E+05   6.186190016E+04
  2.870767672E+04   1.354593172E+04   6.493475213E+03   3.160915774E+03
  1.566402942E+03   7.828926638E+02   4.032058295E+02   2.129754819E+02
  1.130317381E+02   6.157486221E+01   3.081682545E+01   1.745962073E+01
  8.406496322E+00   4.569965531E+00   1.838611962E+00   8.939622677E-01
  2.577221550E-01   9.512584750E-02
  4.229507716E+07   1.141955059E+07   3.385481418E+06   1.077807478E+06
  3.626176904E+05   1.278139940E+05   4.700831502E+04   1.803351286E+04
  7.236424880E+03   3.047220410E+03   1.346097801E+03   6.212385772E+02
  2.974547576E+02   1.444719950E+02   7.281587369E+01   3.559708607E+01
  1.847008238E+01   8.782215331E+00   4.364365927E+00   1.807522394E+00
  7.635699789E-01   2.235716213E-01   5.494900450E-02
  5.139188102E+04   1.215627481E+04   3.972026569E+03   1.553402513E+03
  6.804362135E+02   3.219877107E+02   1.600788432E+02   8.229169547E+01
  4.263535587E+01   2.234970892E+01   1.166226200E+01   5.953952373E+00
  2.940992133E+00   1.282899982E+00   5.318130871E-01   1.946401187E-01
  1.424804508E+03   4.735984924E+02   1.985441055E+02   9.228016323E+01
  4.523210419E+01   2.271190558E+01   1.138109376E+01   5.546411754E+00
  2.573192362E+00   1.061568660E+00
  AL
  END

**Questions**:

1. Can you come up with an explanation for the differences in the spin-orbit splitting calculated with an *j*- and *l*-optimized basis set?

2. Which calculation yields a value closer to the numerical reference and why?

