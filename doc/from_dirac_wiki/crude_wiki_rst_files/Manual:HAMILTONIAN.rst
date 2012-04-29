Manual:HAMILTONIAN
==================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*\*HAMILTONIAN -- Specification of the
   Hamiltonian <#.2A.2AHAMILTONIAN_--_Specification_of_the_Hamiltonian>`_

   -  `1.1 .LVCORR <#.LVCORR>`_
   -  `1.2 .GAUNT <#.GAUNT>`_
   -  `1.3 .DOSSSS <#.DOSSSS>`_
   -  `1.4 .X2C <#.X2C>`_
   -  `1.5 .X2Cmmf <#.X2Cmmf>`_
   -  `1.6 .SPINFREE <#.SPINFREE>`_
   -  `1.7 .LEVY-LEBLOND <#.LEVY-LEBLOND>`_
   -  `1.8 .NONREL <#.NONREL>`_
   -  `1.9 .DFT <#.DFT>`_
   -  `1.10 .DFTAUTO <#.DFTAUTO>`_
   -  `1.11 .LVNEW <#.LVNEW>`_
   -  `1.12 .X2COLD <#.X2COLD>`_
   -  `1.13 .X2C4 <#.X2C4>`_
   -  `1.14 .NOSPIN <#.NOSPIN>`_
   -  `1.15 .HFXFAC <#.HFXFAC>`_
   -  `1.16 .ALDA <#.ALDA>`_
   -  `1.17 .OPERATOR <#.OPERATOR>`_
   -  `1.18 .INTFLG <#.INTFLG>`_
   -  `1.19 .BSS <#.BSS>`_
   -  `1.20 .ONESTEP <#.ONESTEP>`_
   -  `1.21 .NOAMFI <#.NOAMFI>`_
   -  `1.22 .DO2C4C <#.DO2C4C>`_
   -  `1.23 .DO4C2C <#.DO4C2C>`_
   -  `1.24 .USE\_DF <#.USE_DF>`_
   -  `1.25 .CONT2C <#.CONT2C>`_
   -  `1.26 .ZORA <#.ZORA>`_
   -  `1.27 .SOLVENT <#.SOLVENT>`_
   -  `1.28 .FDE <#.FDE>`_
   -  `1.29 .NOSFMU <#.NOSFMU>`_
   -  `1.30 .PRINT <#.PRINT>`_
   -  `1.31 .ONESYS <#.ONESYS>`_
   -  `1.32 .URKBAL <#.URKBAL>`_
   -  `1.33 .FREEPJ <#.FREEPJ>`_
   -  `1.34 .VEXTPJ <#.VEXTPJ>`_
   -  `1.35 .NOSMLV <#.NOSMLV>`_
   -  `1.36 .SMLV1C <#.SMLV1C>`_
   -  `1.37 .ONECAP <#.ONECAP>`_
   -  `1.38 .ONECNV <#.ONECNV>`_
   -  `1.39 .CMPEIG <#.CMPEIG>`_

-  `2 \*DFT <#.2ADFT>`_
-  `3 \*AMFI <#.2AAMFI>`_
-  `4 \*FDE <#.2AFDE>`_

   -  `4.1 Sample inputs <#Sample_inputs>`_
   -  `4.2 References <#References>`_

\*\*HAMILTONIAN -- Specification of the Hamiltonian
===================================================

This section defines the electronic Hamiltonian that is to be used.
Internally the program will always work with 4-component operators that
are expanded using distinct large and small component basis sets. In the
transformation to an orthogonal basis set one may, however, combine
large and small component functions and/or functions of different
symmetry in order to obtain a matrix expansion of e.g. the spinfree
modified Dirac equation or the Lévy-Leblond equation
\ :sup:``[1] <#_note-visscher-saue-00>`_`\ .

In addition one can also modify the Hamiltonian by introduction of an
additional operator, e.g. describing an external field. Any additional
operator defined in this section must be totally symmetric under both
the molecular point group and time reversal symmetry. The latter
requirement precludes the introduction of external magnetic fields.

.LVCORR
-------

This keyword gives Dirac-Coulomb Hamiltonian in which (SS\|SS) integrals
are neglected and replaced by an interatomic SS correction (calculated
as a classical repulsion term of (tabulated) small component atomic
charges) \ :sup:``[2] <#_note-visscher-97>`_`\ .

This is currently the most economical and accurate approximation to the
full Dirac-Coulomb Hamiltonian and can certainly be used for the
calculation of spectroscopic constants and valence properties; for core
properties testing is recommended (see also `.LVNEW <#.LVNEW>`_). This
is the default Hamiltonian choice in DIRAC11.

.GAUNT
------

Add the Gaunt interaction to the Hamiltonian. This will increase the
computational time significantly but is important when studying
spin-orbit splittings and/or performing accurate studies of light
molecules. The current implementation is limited to including the Gaunt
interaction in the construction of the Fock matrix and works for
Hartree-Fock and DFT. For using Gaunt in combination with DFT see the
`\*DFT </dirac/index.php/Manual:HAMILTONIAN:DFT>`_ section of the
manual. Transformation of the Gaunt part of the two electron operator to
the MO basis is not yet implemented, for this purpose we recommend the
use of the molecular mean field approximation. This can be used for
example in MP2/CC/FSCC/IHFSCC calculations with RELCCSD by means of the
`X2Cmmf <#.X2Cmmf>`_ Hamiltonian (see also our FAQ page:
`[1] <http://tinyurl.com/cu8gomm>`_).

For the relativistic two-component mode (see `.X2C <#.X2C>`_ keyword)
with AMFI contributions it uses both spin-same and spin-other orbit
mean-field parts.

.DOSSSS
-------

This keyword gives unmodified Dirac-Coulomb Hamiltonian which was
default till the DIRAC11 release. Explicitly including (SS\|SS) type
Coulomb integrals does give the most accurate description of the system
but does increase computational cost significantly. Use this option for
high-accuracy calculations, preferably in conjunction with the .GAUNT
keyword to also include the Gaunt correction to the two-electron
interaction.

.X2C
----

Keyword referring to the new module for the *E*xact *Two*-*C*omponent
relativistic Hamiltonian. \ :sup:``[3] <#_note-ilias:iotc>`_`\ 
\ :sup:``[4] <#_note-knecht:x2c_mod>`_`\ .

DIRAC runs in 2-component mode. One should combine this option with an
AMFI correction to the unscreened spin-orbit operator. This is the
default in the development version unless `.NOAMFI <#.NOAMFI>`_
specified. See also the warning there.

.X2Cmmf
-------

Keyword referring to the 2-component molecular-mean-field Hamiltonian
apporach within the new module for the *E*xact *Two*-*C*omponent
relativistic Hamiltonian. \ :sup:``[5] <#_note-sikkema_mmf>`_`\ 
\ :sup:``[4] <#_note-knecht:x2c_mod>`_`\ .

DIRAC starts with a 4c-SCF run and performs a transformation to
2-component mode (based on the converged Fock operator) prior to a
post-HF correlation step. after the SCF . One can combine this option
with the .GAUNT keyword which activates the inclusion of
spin-other-orbit contributions in the Hamiltonian. The
X2Cmmf-Hamiltonian can at present (DIRAC11 release) only be used for
post-HF calculation within the RELCCSD module. Patches for other
correlation modules in DIRAC are planned to be released soon. For
further information and a detailed guideline please consult our `FAQ
page </dirac/index.php/FAQ>`_.

.SPINFREE
---------

Use Dyall's spinfree Hamiltonian \ :sup:``[6] <#_note-dyall-94>`_`\  to
obtain results without spin-orbit coupling for the four-component
Hamiltonian in the default restricted kinetic balance scheme. This
keyword works also for two-component relativistic Hamiltonians where one
can choose between two spin-free schemes - see the `.BSS <#.BSS>`_
keyword.

Note that this option should not be used for response calculations with
time-antisymmetric (magnetic) operators as it will eliminate important
contributions. Add more...

.LEVY-LEBLOND
-------------

Use the nonrelativistic Lévy-Leblond Hamiltonian
\ :sup:``[7] <#_note-levy-leblond-67>`_`\ .

Use this option before any additional one-electron operators are
specified, because it redefines the metric used in the calculation.
Lévy-Leblond provides nonrelativistic mode in DIRAC.

.NONREL
-------

Standard nonrelativistic calculation based on the Schrodinger equation.
Should give identical energy results as with the
`.LEVY-LEBLOND <#.LEVY-LEBLOND>`_ keyword.

DIRAC runs in the 2-component mode.

.DFT
----

Perform a Kohn-Sham calculation. In the following line you must specify
the desired DFT functional.

The functional can either be selected from a set of predefined
combinations of exchange and correlation functionals, e.g.:

::

    .DFT
     B3LYP

or composed by specifying GGAKEY followed by a list of the desired
functionals together with their weights:

::

    .DFT
     GGAKEY PW86X=1.0 P86C=1.0

*Exchange functionals:*

-  LDA exchange \ :sup:``[8] <#_note-lda:exchange>`_`\ : SLATER
-  Becke 1988 \ :sup:``[9] <#_note-gga:B88>`_`\ : BECKE
-  The correction term to LDA proposed by Becke
   \ :sup:``[9] <#_note-gga:B88>`_`\ : B88CORR
-  Perdew-Wang 1986 \ :sup:``[10] <#_note-gga:pw86x>`_`\ : PW86X

*Correlation functionals:*

-  The Vosko-Wilk-Nusair LDA correlation energy (VWN5)
   \ :sup:``[11] <#_note-lda:vwn>`_`\ : VWN
-  The Lee-Yang-Parr functional \ :sup:``[12] <#_note-gga:lyp>`_`\ : LYP
-  Perdew 1986 \ :sup:``[13] <#_note-gga:p86c>`_`\ : P86C

*Some of the predefined combinations (see the `Dalton
manual <http://www.daltonprogram.org>`_ for a more complete list):*

-  LDA
-  CAMB3LYP
-  BLYP
-  B3LYP
-  BP86
-  BPW91
-  PBE
-  PBE0

.DFTAUTO
--------

Perform a Kohn-Sham calculation using functionals provided by the XCFun
library. In the following line you must specify the desired DFT
functional.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.LVNEW
------

Modification of `.LVCORR <#.LVCORR>`_ that obtains the atomic small
component charge via a Mulliken analysis instead of the original table
look-up. (The problem with the table look-up is that the electrostatics
in the molecule will be wrong if you have specified a basis set which
does not give the correct small-electron charge because of deficiencies
in the core region.)

.X2COLD
-------

One-step Exact (infinite order) 2-Component relativistic Hamiltonian
\ :sup:``[3] <#_note-ilias:iotc>`_`\ .

**note: in the previous release version Dirac2008 this keyword is called
.X2C.**: DIRAC runs in the (memory saving) 2-component mode. Note that
one should combine this option with the SPINFREE option as X2C will only
provide an unscreened (bare nucleus) spin-orbit operator that gives
unphysically large spin-orbit contributions. To get a realistic screened
spin-orbit operator AMFI is added in the development version unless
`.NOAMFI <#.NOAMFI>`_ specified. See also the warning there.

.X2C4
-----

One-step Exact (infinite order) 2-Component relativistic Hamiltonian
\ :sup:``[3] <#_note-ilias:iotc>`_`\ .

DIRAC runs in the 4-component mode. This mode is useful if you wish to
restart from a previous 4-component calculation and vice versa. AMFI is
added in the development version unless `.NOAMFI <#.NOAMFI>`_ specified.
See also the warning there.

.NOSPIN
-------

Implies .SPINFREE, but also remove all spin-symmetry-breaking
(quaternion "imaginary" or "triplet" terms) from property gradients in
response calculations. Used for analyzing magnetic properties similarly
to how it is done with non-relativistic methods.

.HFXFAC
-------

Weight of exchange in Fock matrix construction.

*Default:*

::

    .HFXFAC
     1.0

.ALDA
-----

Invoke the adiabatic LDA approximation: anything beyond the Kohn-Sham
functional (exchange-correlation kernel etc.) is approximated by LDA
(exchange switched off in hybrid functionals). In the development
version this keyword has been moved after
`\*DFT </dirac/index.php/Manual:HAMILTONIAN:DFT>`_.

.OPERATOR
---------

Specification of an additional one-electron operator in the Hamiltonian.
The operator must be totally symmetric both under the molecular point
group and time reversal symmetry. The field strength of the operator is
specified with COMFACTOR. The keyword can be repeated for addition of
more than one operator.

See the `specification of one-electron
operators </dirac/index.php/Manual:One-electron_operators>`_ section for
more information and explicit examples.

.INTFLG
-------

Specify what two-electron integrals to include. All other modules use
this as the default value.

*Default:* Include LL, SL, and SS integrals (1 = include; 0 = do not
include):

::

    .INTFLG
     1 1 0

.BSS
----

Use the 2-component relativistic Hamiltonian obtained after the
Barysz-Sadlej-Snijders transformation of the Dirac Hamiltonian in the
finite basis set.\ :sup:``[14] <#_note-ilias:pbo>`_`\ .

Calculations using the 2-component BSS Hamiltonian are running only with
large component basis functions.

The user can select order of scalar and spin-orbit terms in the .BSS
Hamiltonian.

*Both scalar and spin-orbit relativistic effects up to infinite order:*

::

    .BSS
     099

*Scalar relativistic effects of type "from the beginning" up to the
infinite order, no spin-orbit interaction:*

::

    .SPINFREE
    .BSS
     109

*Scalar relativistic effects "from the end" up to the infinite order:*

::

    .SPINFREE
    .BSS
     009

*Traditional scalar relativistic "from the beginning" second-order
Douglas-Kroll-Hess Hamiltonian:*

::

    .SPINFREE
    .BSS
     102

*Scalar relativistic "from the end" second-order Douglas-Kroll-Hess
Hamiltonian:*

::

    .SPINFREE
    .BSS
     002

*Douglas-Kroll-Hess Hamiltonian with first-order spin-orbit and
second-order scalar relativistic effects:*

::

    .BSS
     012

*Douglas-Kroll-Hess Hamiltonian with second-order spin-orbit and
second-order scalar relativistic effects:*

::

    .BSS
     022

`|image1| </dirac/index.php/Image:Dirac.png>`_

**Only in the development version of DIRAC**

In this two-component variational scheme is possible to combine external
AMFI spin-orbit terms\ :sup:``[15] <#_note-ilias:amfi>`_`\  together
with BSS integrals - see the
`\*AMFI </dirac/index.php/Manual:HAMILTONIAN:AMFI>`_ section. Note that
AMFI provides only one-center integrals.

*Infinite order scalar and spin-orbit terms, (one-electron) mean-field
spin-same-orbit (MFSSO2) from AMFI*

::

    .BSS
    2999

This is in fact the 'maximum' of two-component relativity, resembling
Dirac-Coulomb Hamiltonian.

*Second-order Douglas-Kroll-Hess spin-free from 'the beginning' with
first order spin-orbit (SO1) term plus mean-field spin-same-orbit
(MFSSO) from AMFI*

::

    .BSS
    2112

SO1 may come either from BSS-transformation or from AMFI. In the latter
case it is only for one-center and for point nucleus.

*Infinite order scalar and spin-orbit terms, (one-electron) mean-field
spin-same and spin-other-orbit (MFSO2) from AMFI*

::

    .BSS
    3999

This mimics the Dirac-Coulomb-Gaunt Hamiltonian, as the spin-other-orbit
comes from the Gaunt interaction term.

DIRAC enables switching-off AMFI spin-orbit contributions from various
centers. See keyword
`.NOAMFIC </dirac/index.php/Manual:HAMILTONIAN:AMFI#.NOAMFC>`_.

*Infinite order scalar terms "from the beginning", (one-electron) spin
orbit (SO1) and the mean-field spin-same-orbit (MFSSO2) terms
exclusively from AMFI*

::

    .BSS
    4109

*Second order (DKH) scalar terms "from the beginning", (one-electron)
spin orbit (SO1) and the mean-field spin-same-orbit (MFSSO2) terms
exclusively from AMFI*

::

    .BSS
    4102

*Infinite order scalar terms "from the end", (one-electron) spin orbit
(SO1) and the mean-field spin-same-orbit and spin-other-orbit (MFSO2)
terms from AMFI*

::

    .BSS
    5009

*Infinite order scalar terms "from the end", and the (one-electron) spin
orbit term (SO1) from AMFI*

::

    .BSS
    6009

AMFI one-electron spin-orbit terms - SO1 - are currently for the point
nucleus.

For the BSS value 'axyz' of y=0, DIRAC employs spin-free picture change
transformation of property operators, although the system is not in the
boson (spin-free) symmetry for a>1.

*Rough first order, DKH1 (not recommended for practical calculations)*

::

    .BSS
    111

For comparison only between BSS-SO1 and AMFI-SO1 integrals (point
nucleus only) use

::

    .BSS
    001

.ONESTEP
--------

Together with the `.BSS <#.BSS>`_ keyword - for the infinite order only
- invokes the one-step infinite order method (which is otherwise called
by `.X2C <#.X2C>`_, `.X2C4 <#.X2C4>`_ keywords).

.NOAMFI
-------

Do not include the AMFI contribution in the development version where
AMFI is the default. This holds also for keywords `.X2C <#.X2C>`_ and
`.X2C4 <#.X2C4>`_). In the DIRAC08 distribution version .NOAMFI is the
default.

Warning: missing AMFI contribution (for `.X2C <#.X2C>`_,
`.X2C4 <#.X2C4>`_ and `.BSS <#.BSS>`_ keywords) may lead to
overestimation of spin-orbit effects since these would be represented by
one-electron terms only, without two-electron shielding.

.DO2C4C
-------

`|image2| </dirac/index.php/Image:Dirac.png>`_

**Only in the development version of DIRAC**

After iterations at the two-component level *ascend* to the
four-component level.

.DO4C2C
-------

`|image3| </dirac/index.php/Image:Dirac.png>`_

**Only in the development version of DIRAC**

After iterations at the four-component level do the relativistic
transformation to the two-component level.

.USE\_DF
--------

`|image4| </dirac/index.php/Image:Dirac.png>`_

**Only in the development version of DIRAC**

After the four-component DC-SCF method do the infinite order
transformation (either one- or two-step) upon the Fock-Dirac matrix.
Otherwise it is transforming Dirac bare nucleus.

Used only with keyword `.DO4C2C <#.DO4C2C>`_.

.CONT2C
-------

`|image5| </dirac/index.php/Image:Dirac.png>`_

**Only in the development version of DIRAC**

After four-component DC-SCF continue with two-component iterations.

Used only with keyword `.DO4C2C <#.DO4C2C>`_.

Integer (3,4,5) should be specified in free format on the line following
`.CONT2C <#.CONT2C>`_:

::

    .CONT2C
     3

.ZORA
-----

Use the zeroth-order regular approximation
\ :sup:``[16] <#_note-lenthe:zora1>`_`\ 
\ :sup:``[17] <#_note-lenthe:zora2>`_`\ 
\ :sup:``[1] <#_note-visscher-saue-00>`_`\  of the Dirac Hamiltonian in
the Hartree-Fock procedure. Works only for closed-shell systems. The
implementation offers only little computational advantages and is
intended chiefly for comparisons of methodology. Note that the
combination `.SPINFREE <#.SPINFREE>`_ and `.ZORA <#.ZORA>`_ gives a
spinfree formalism that differs from the conventional spinfree ZORA
formulation. Two integers should be specified in free format on the line
following `.ZORA <#.ZORA>`_:

::

    .ZORA
     1 1

The first number indicates whether the density is to be normalized over
the 2-component (0; ZORA) or 4-component metric (1; ZORA4).

The second number specifies whether the orbital energies should be
unmodified (0; normal ZORA) or scaled (1; scaled ZORA).

.SOLVENT
--------

Model solvent effects by placing the molecule in a spherical cavity in a
dielectric continuum.

.FDE
----

Activates the frozen density embedding (FDE) functionality. Options can
be specified under the \*FDE menu.

In order to use FDE, the user must generated an embedding potential with
the ADF code (see the developer's website [`[2] <http://www.scm.com>`_]
for further information).

.NOSFMU
-------

`|image6| </dirac/index.php/Image:Dirac.png>`_

**Only in the development version of DIRAC**

In spinfree correlated calculations group multiplication tables are by
default set up as direct products of spatial and spin symmetries. This
flag turns off this, and so the spinfree case is treated similar to the
spin-orbit case.

`|image7| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.ONESYS
-------

Ignore two-particle interactions.

The keyword ensures the diagonalization of the bare nuclei Hamiltonian
matrix without proceeding to the iterative SCF method.

.URKBAL
-------

Unrestricted kinetic balance.

The default is restricted kinetic balance. This is imposed by deleting
unphysical solutions from the free particle positronic spectrum. This
leads to a 1:1 ratio of electronic and positronic solutions. This
preprojection is sensitive to linear dependencies and should therefore
preferably be used in conjunction with the spherical transformation of
both large and small components.

.FREEPJ
-------

Project out all free positronic solutions from the MO space.

.VEXTPJ
-------

Project out all external field positronic solutions from the MO space.

.NOSMLV
-------

Delete SS nuclear attraction integrals. This will take out contributions
to the one-electron spin-orbit interaction and the Darwin interaction.

.SMLV1C
-------

Neglect potential for multi-center SS blocks, i.e. multi-center SS
nuclear attraction integrals and multi-center SS two-electron integrals.

.ONECAP
-------

Consider taking the `.SMLV1C <#.SMLV1C>`_ model one step further. Only
one-center contributions to the LS and SS two-electron integrals and SS
nuclear attraction integrals are calculated explicitly. The
electrostatic effects of the terms neglected this way are included by
calculating the classical repulsion from small component charges based
on a Mulliken population analysis. Note that we therefore only need to
calculate the derivative of the LL integrals when calculating the
molecular gradient.

.ONECNV
-------

::

    .ONECNV
     THRESHOLD

Employ the one-center model as given by `.ONECAP <#.ONECAP>`_ until
convergence to a specified THRESHOLD, whereafter the full set of
two-electron integrals will be used. This threshold applies to whatever
convergence criteria has been selected
(`.EVCCNV </dirac/index.php/Manual:WAVE_FUNCTION:DHFCAL#.EVCCNV>`_,
`.ERGCNV </dirac/index.php/Manual:WAVE_FUNCTION:DHFCAL#.ERGCNV>`_, or
`.FCKCNV </dirac/index.php/Manual:WAVE_FUNCTION:DHFCAL#.FCKCNV>`_).

.CMPEIG
-------

When some two-component relativistic Hamiltonian is chosen, compare
eigenvalues between the 'parent' four-component Dirac and derived
two-component one-electron Hamiltonians.

For the infinite order (one- and two-step) two-component Hamiltonians
eigenvalues are identical with four-component Dirac counterparts. For
the second-order (and lower order) Douglas-Krol-Hess Hamiltonian they
slightly differ.

`\*DFT </dirac/index.php/Manual:HAMILTONIAN:DFT>`_
==================================================

In this link the directives for modifying the DFT calculation are
described (e.g. modifying grid, using ALDA and much more)

`\*AMFI </dirac/index.php/Manual:HAMILTONIAN:AMFI>`_
====================================================

AMFI directives

`\*FDE </dirac/index.php/Manual:HAMILTONIAN:FDE>`_
==================================================

Frozen density embedding (FDE) directives.

`Sample inputs </dirac/index.php/Manual:HAMILTONIAN:Sample_inputs>`_
--------------------------------------------------------------------

Some examples of the Hamiltonian input.

References
----------

#. ↑ \ :sup:``1.0 <#_ref-visscher-saue-00_0>`_`\ 
   \ :sup:``1.1 <#_ref-visscher-saue-00_1>`_`\  L. Visscher and T. Saue,
   *Approximate relativistic electronic structure methods based on the
   quaternion modified Dirac equation*, J. Chem. Phys. **113**, *3996*
   (2000) `electronic version <http://dx.doi.org/10.1063/1.1288371>`_.
#. `↑ <#_ref-visscher-97_0>`_ L. Visscher, *Approximate molecular
   Dirac-Coulomb calculations using a simple Coulombic correction*,
   Theor. Chem. Acc. **98**, *68* (1997) `electronic
   version <http://dx.doi.org/10.1007/s002140050280>`_.
#. ↑ \ :sup:``3.0 <#_ref-ilias:iotc_0>`_`\ 
   \ :sup:``3.1 <#_ref-ilias:iotc_1>`_`\ 
   \ :sup:``3.2 <#_ref-ilias:iotc_2>`_`\  M. Iliaš and T. Saue, *An
   infinite-order two-component relativistic Hamiltonian by a simple
   one-step transformation*, J. Chem. Phys. **126**, *064102* (2007)
   `electronic version <http://dx.doi.org/10.1063/1.2436882>`_.
#. ↑ \ :sup:``4.0 <#_ref-knecht:x2c_mod_0>`_`\ 
   \ :sup:``4.1 <#_ref-knecht:x2c_mod_1>`_`\  S. Knecht and T. Saue,
   *tba*, xxx **xxx**, *Universities of Strasbourg and Odense*
   (2010-2011) .
#. `↑ <#_ref-sikkema_mmf_0>`_ J. Sikkema, L. Visscher, T. Saue, and M.
   Iliaš, *The molecular mean-field approach for correlated relativistic
   calculations*, J. Chem. Phys. **131**, *124116* (2009) .
#. `↑ <#_ref-dyall-94_0>`_ K. G. Dyall, *An exact separation of the
   spin-free and spin-dependent terms of the Dirac-Coulomb-Breit
   hamiltonian*, J. Chem. Phys. **100**, *2118* (1994) `electronic
   version <http://dx.doi.org/10.1063/1.466508>`_.
#. `↑ <#_ref-levy-leblond-67_0>`_ J.-M. Lévy-Leblond, *Nonrelativistic
   particles and wave equations*, Commun. Math. Phys. **6**, *286*
   (1967) `electronic version <http://dx.doi.org/10.1007/BF01646020>`_.
#. `↑ <#_ref-lda:exchange_0>`_ P. A. M. Dirac, *Note on Exchange
   Phenomena in the Thomas Atom*, Proc. Roy. Soc. London **26**, *376*
   (1930) .
#. ↑ \ :sup:``9.0 <#_ref-gga:B88_0>`_`\ 
   \ :sup:``9.1 <#_ref-gga:B88_1>`_`\  A. D. Becke, *Density-Functional
   Exchange-Energy Approximation With Correct Asymptotic Behavior*,
   Phys. Rev. A **38**, *3098* (1988) `electronic
   version <http://dx.doi.org/10.1103/PhysRevA.38.3098>`_.
#. `↑ <#_ref-gga:pw86x_0>`_ J. P. Perdew and Y. Wang, *Accurate and
   Simple Density Functional for the Electronic Exchange Energy:
   Generalized Gradient Approximation*, Phys. Rev. B **33**, *8800*
   (1986) `electronic
   version <http://dx.doi.org/10.1103/PhysRevB.33.8800>`_.
#. `↑ <#_ref-lda:vwn_0>`_ S. J. Vosko, L. Wilk, and M. Nusair, *Accurate
   Spin-Dependent Electron Liquid Correlation Energies for Local Spin
   Density Calculations: A Critical Analysis*, Can. J. Phys. **58**,
   *1200* (1980) .
#. `↑ <#_ref-gga:lyp_0>`_ C. Lee, W. Yang, and R. G. Parr, *Development
   of the Colle-Salvetti Correlation-Energy Formula into a Functional of
   the Electron Density*, Phys. Rev. B **37**, *785* (1988) `electronic
   version <http://dx.doi.org/10.1103/PhysRevB.37.785>`_.
#. `↑ <#_ref-gga:p86c_0>`_ J. P. Perdew, *Density-Functional
   Approximation for the Correlation Energy of the Inhomogenous Electron
   Gas*, Phys. Rev. B **33**, *8822* (1986) `electronic
   version <http://dx.doi.org/10.1103/PhysRevB.33.8822>`_.
#. `↑ <#_ref-ilias:pbo_0>`_ M. Iliaš and H. J. Aa. Jensen and V. Kello
   and B. O. Roos and M. Urban, *Theoretical study of PbO and the PbO
   anion*, Chem.Phys.Lett. **408**, *210* (2005) `electronic
   version <http://dx.doi.org/doi:10.1016/j.cplett.2005.04.027>`_.
#. `↑ <#_ref-ilias:amfi_0>`_ M. Iliaš and V. Kello and L. Visscher and
   B. Schimmelpfennig, *Inclusion of mean-field spin–orbit effects based
   on all-electron two-component spinors: Pilot calculations on atomic
   and molecular properties*, J.Chem.Phys. **115**, *9667* (2001)
   `electronic version <http://dx.doi.org/DOI:10.1063/1.1413510>`_.
#. `↑ <#_ref-lenthe:zora1_0>`_ E. van Lenthe, E. J. Baerends, and J. G.
   Snijders, *Relativistic total energies using regular approximations*,
   J. Chem. Phys. **101**, *9783* (1994) `electronic
   version <http://dx.doi.org/10.1063/1.467943>`_.
#. `↑ <#_ref-lenthe:zora2_0>`_ E. van Lenthe, J. G. Snijders, and E. J.
   Baerends, *The zero-order regular approximation for relativistic
   effects: The effect of spin-orbit coupling in closed shell
   molecules*, J. Chem. Phys. **105**, *6505* (1996) `electronic
   version <http://dx.doi.org/10.1063/1.472460>`_.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:HAMILTONIAN <http://wiki.chem.vu.nl/dirac/index.php/Manual:HAMILTONIAN>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:HAMILTONIAN>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:HAMILTONIAN&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:HAMILTONIAN&action=edit>`_
-  `History </dirac/index.php?title=Manual:HAMILTONIAN&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:HAMILTONIAN>`_

` </dirac/index.php/Dirac_Program>`_

Navigation
''''''''''

-  `Home </dirac/index.php/Dirac_Program>`_
-  `Documentation </dirac/index.php/Documentation>`_
-  `Reference manual </dirac/index.php/Manual>`_
-  `Discussion forum </dirac/index.php/DiscussionForum>`_
-  `News Archive </dirac/index.php/News_archive>`_
-  `Basis sets </dirac/index.php/Basis_Sets>`_

Search
''''''

Toolbox
'''''''

-  `What links
   here </dirac/index.php/Special:Whatlinkshere/Manual:HAMILTONIAN>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:HAMILTONIAN>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:HAMILTONIAN&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:HAMILTONIAN&oldid=34641>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 10:40, 27 April 2012.
-  This page has been accessed 5,104 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image1| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image2| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image3| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image4| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image5| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image6| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image7| image:: /dirac/images/thumb/6/65/HJAaJ.png/41px-HJAaJ.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
