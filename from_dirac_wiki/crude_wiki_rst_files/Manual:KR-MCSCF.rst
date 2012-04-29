Manual:KR-MCSCF
===============

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

*2- and 4-component relativistic KR-MCSCF module*
\ :sup:``[1] <#_note-jensen:mcscf>`_`\ 
\ :sup:``[2] <#_note-thyssen:thesis>`_`\ 
\ :sup:``[3] <#_note-thyssen:mcscf>`_`\ 
\ :sup:``[4] <#_note-knecht:thesis>`_`\  original implementation written
by Joern Thyssen, Timo Fleig and Hans Joergen Aagaard Jensen,
parallelization and continuous improvement by Stefan Knecht and Hans
Joergen Aagaard Jensen. Linear symmetry adaptation by Stefan Knecht and
Hans Joergen Aagaard Jensen.

**stefan: FIXME: a more detailed introduction would be most welcome.**

Contents
--------

-  `1 Mandatory keywords <#Mandatory_keywords>`_

   -  `1.1 .MK2REF <#.MK2REF>`_
   -  `1.2 .MK2DEL <#.MK2DEL>`_
   -  `1.3 .GASSH <#.GASSH>`_
   -  `1.4 .GASSPC <#.GASSPC>`_

-  `2 optional keywords under
   \*KRMCSCF <#optional_keywords_under_.2AKRMCSCF>`_

   -  `2.1 .THRESH <#.THRESH>`_
   -  `2.2 .SYMMETRY <#.SYMMETRY>`_
   -  `2.3 .CI PROGRAM <#.CI_PROGRAM>`_
   -  `2.4 .INACTIVE <#.INACTIVE>`_
   -  `2.5 .MAX MACRO <#.MAX_MACRO>`_
   -  `2.6 .MAX MICRO <#.MAX_MICRO>`_
   -  `2.7 .DELETE <#.DELETE>`_
   -  `2.8 .SKIPEE <#.SKIPEE>`_
   -  `2.9 .SKIPEP <#.SKIPEP>`_
   -  `2.10 .PRINT <#.PRINT>`_

-  `3 optional keywords under
   \*OPTIMI <#optional_keywords_under_.2AOPTIMI>`_

   -  `3.1 .MAX CI <#.MAX_CI>`_
   -  `3.2 .MXCIVE <#.MXCIVE>`_
   -  `3.3 .RSTRCI <#.RSTRCI>`_
   -  `3.4 .NOOCCN <#.NOOCCN>`_
   -  `3.5 .ANALYZ <#.ANALYZ>`_

-  `4 References <#References>`_

Mandatory keywords
==================

**NOTE: the following keywords must be placed under the input deck**

::

    *KRMCSCF

.MK2REF
-------

Reference value |M\_K^{ref}|, given as |2\*M\_K^{ref}| (integer value);
where *M*\ :sub:`*K*`\  is the number of Kramers orbitals
φ\ :sub:`*p*`\  minus the number of Kramers-paired orbitals
|\\phi\_{\\overline{p}}| in a given determinant.

*Default:* None.

.MK2DEL
-------

value of 2 \* Δ*M*\ :sub:`*K*`\  (integer value). The range of allowed
*M*\ :sub:`*K*`\  values is |2\*\\left(M\_K^{ref}\\pm\\Delta
M\_K\\right)|

*Default:* None.

.GASSH
------

Number and specification of GAS orbitals. Line with the number of GA
spaces used (1-7), followed by one line per GAS with number of orbitals
per fermion correp.

*Default:* None.

.GASSPC
-------

specification of CI calculation. One line per GAS with 2 numbers each:
The first gives the minimal number of accumulated electrons after this
GAS, the second the corresponding maximum number, separated by blanks
(defining occupation constraints of each GAS).

*Default:* None.

optional keywords under \*KRMCSCF
=================================

.THRESH
-------

convergence threshold in the MCSCF gradient given as double precision
value. The present default value is quite tight and might be adapted if
one is only interested in MCSCF energies.

*Default:*

::

    0.5D-05

.SYMMETRY
---------

integer value giving the (boson/fermion) symmetry of the state to
optimize on.

*Default:*

::

    1

If you run the calculation in **linear symmetry** (**only available in
the development version**) you have to specify the **2 x OMEGA** value
of the state to optimize on. The doubling (*2 x*) stems from the fact
that we want to avoid non-integer input, e.g. in case of an odd number
of electrons we might have OMEGA=1/2, 3/2, 5/2, etc. values and the
corresponding input for the OMEGA=1/2 would then read as

::

    1

If we have a system with an even number of electrons and inversion
symmetry the input for an OMEGA=2g state would read as

::

    4g

.CI PROGRAM
-----------

character string reading either *LUCIAREL* or *GASCIP*. This keyword
enables to choose between the two CI programs that are implemented for
the MCSCF. Note that for large-scale MCSCF (| > 1.0 \\times 10^7|
determinants) only LUCIAREL can be used.

*Default:*

::

    GASCIP

.INACTIVE
---------

Inactive orbitals per fermion correp given as integer values separated
by blanks.

*Default:* All orbitals active.

.MAX MACRO
----------

integer value giving the maximum number of MACRO iterations in the MCSCF
optimization.

*Default:*

::

    25

.MAX MICRO
----------

integer value giving the maximum number of MICRO iterations in the MCSCF
optimization.

*Default:*

::

    50

.DELETE
-------

delete particular active-secondary e-e rotations (rotations between
electronic-electronic spinors) in the gradient and Hessian calculation
specified by an orbital string of virtual (secondary) orbitals for each
fermion correp.

*Default:*

::

    include all active-secondary e-e rotations.

*example:* (deleting all e-e rotations for secondary orbitals 20,21 and
22 in fermion correp 1 (*gerade*) and 24,25,26,27 and 28 in fermion
correp 2 (*ungerade*)):

::

    .DELETE
    20,21,22
    24..28

.SKIPEE
-------

skip e-e rotations (rotations between electronic-electronic spinors) in
the gradient and Hessian calculation.

*Default:*

::

    include e-e rotations.

.SKIPEP
-------

skip e-p rotations (rotations between electronic-positronic spinors) in
the gradient and Hessian calculation.

*Default:*

::

    include e-p rotations.

.PRINT
------

raise default print level to the given integer value. Please use with
care as you may get millions of output lines if you choose a too high
value.

*Default:*

::

    print level is 0.

optional keywords under \*OPTIMI
================================

**NOTE: the following keywords must be placed under the input deck**

::

    *OPTIMI

.MAX CI
-------

Number of initial CI iterations.

*Default:*

::

    .MAX CI
     5

.MXCIVE
-------

size of Davidson subspace.

example:

::

    .MXCIVE
     24

*Default:* 3 times # of states to optimize on (remember this is a
state-selective MCSCF), therefore:

::

    3 

.RSTRCI
-------

Optional restart from CI vector on file KRCI\_CVECS.x where **x** is
determined by the symmetry of the wave function (symmetry == 1 --> **x =
1**; symmetry == 2 --> **x = 2**; etc.). **NOTE:** in the DIRAC10
release version the **x** <==> symmetry relation is: symmetry == 1 -->
**x = a**; symmetry == 2 --> **x = b**; etc.

::

    .RSTRCI
     1

*Default:* No restart.

.NOOCCN
-------

compute natural orbital occupation numbers for the electronic state.

*Default:*

::

    do not compute natural orbital occupation numbers.

.ANALYZ
-------

analyze the final CI wave function printing the coefficients for each
determinant above a given threshold.

*Default:*

::

    do not analyze the final CI wave function.

References
==========

#. `↑ <#_ref-jensen:mcscf_0>`_ H. J. Aa. Jensen and K. G. Dyall and T.
   Saue and K. Fægri jr., *`Relativistic 4-component
   Multi-Configurational Self-Consistent Field Theory for Molecules:
   Formalism <http://link.aip.org/link/doi/10.1063/1.471644>`_*, J.
   Chem. Phys. **104**, *4083* (1996) .
#. `↑ <#_ref-thyssen:thesis_0>`_ J. Thyssen, *`PhD
   thesis <http://dirac.chem.sdu.dk/thesis/thesis-jth2001.pdf>`_*,
   University of Southern Denmark (Denmark) '***, '** (2004)* .
#. `↑ <#_ref-thyssen:mcscf_0>`_ J. Thyssen and T. Fleig and H. J. Aa.
   Jensen, *`A Direct Relativistic Four-Component Multi-Configuration
   Self-Consistent-Field Method for
   Molecules <http://link.aip.org/link/doi/10.1063/1.2943670>`_*, J.
   Chem. Phys. **129**, *034109* (2008) .
#. `↑ <#_ref-knecht:thesis_0>`_ S. R. Knecht, *`PhD
   thesis <http://dirac.chem.sdu.dk/thesis/09.PhD.StefanKnecht.pdf>`_*,
   University of Duesseldorf (Germany) '***, '** (2009)* .

