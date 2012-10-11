.. _hamiltonian:

.. _\*\*HAMILTONIAN:

\*\*HAMILTONIAN
===============

Specification of the Hamiltonian.

This section defines the electronic Hamiltonian that is to be used.
Internally the program will always work with 4-component operators that
are expanded using distinct large and small component basis sets. In the
transformation to an orthogonal basis set one may, however, combine
large and small component functions and/or functions of different
symmetry in order to obtain a matrix expansion of e.g. the spinfree
modified Dirac equation or the Lévy-Leblond equation [Visscher2000]_.

In addition one can also modify the Hamiltonian by introduction of an
additional operator, e.g. describing an external field. Any additional
operator defined in this section must be totally symmetric under both
the molecular point group and time reversal symmetry. The latter
requirement precludes the introduction of external magnetic fields.

.. _.LVCORR:

.LVCORR
-------

This keyword gives Dirac-Coulomb Hamiltonian in which (SS|SS) integrals
are neglected and replaced by an interatomic SS correction (calculated
as a classical repulsion term of (tabulated) small component atomic
charges) [Visscher1997a]_ .

This is currently the most economical and accurate approximation to the
full Dirac-Coulomb Hamiltonian and can certainly be used for the
calculation of spectroscopic constants and valence properties; for core
properties testing is recommended (see also :ref:`.LVNEW`). This
is the default Hamiltonian choice in DIRAC11.

.. _.GAUNT:

.GAUNT
------

Add the Gaunt interaction to the Hamiltonian. This will increase the
computational time significantly but is important when studying
spin-orbit splittings and/or performing accurate studies of light
molecules. The current implementation is limited to including the Gaunt
interaction in the construction of the Fock matrix and works for
Hartree-Fock and DFT. For using Gaunt in combination with DFT see the
:ref:`\*DFT` section of the
manual. Transformation of the Gaunt part of the two electron operator to
the MO basis is not yet implemented, for this purpose we recommend the
use of the molecular mean field approximation. This can be used for
example in MP2/CC/FSCC/IHFSCC calculations with RELCCSD by means of the
:ref:`.X2Cmmf` Hamiltonian (see also our FAQ/Tutorial pages)

For the relativistic two-component mode (see :ref:`.X2C` keyword)
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

.. _.X2C:

.X2C
----

Keyword referring to the new module for the *E* xact *Two*- *C* omponent
relativistic Hamiltonian, Refs. [Ilias2007]_, [Knecht2010]_.

DIRAC runs in 2-component mode. One should combine this option with an
AMFI correction to the unscreened spin-orbit operator. This is the
default in the development version unless :ref:`.NOAMFI` specified. 
See also the warning there.

.. _.X2Cmmf:

.X2Cmmf
-------

Keyword referring to the 2-component molecular-mean-field Hamiltonian
apporach within the new module for the *E* xact *Two*- *C* omponent
relativistic Hamiltonian, Refs. [Sikkema2009]_, [Knecht2010]_.

DIRAC starts with a 4c-SCF run and performs a transformation to
2-component mode (based on the converged Fock operator) prior to a
post-HF correlation step. after the SCF . One can combine this option
with the .GAUNT keyword which activates the inclusion of
spin-other-orbit contributions in the Hamiltonian. The
X2Cmmf-Hamiltonian can at present (DIRAC11 release) only be used for
post-HF calculation within the RELCCSD module. Patches for other
correlation modules in DIRAC are planned to be released soon. 

For further information and a detailed guideline please consult our FAQ
page.

.. _.SPINFREE:

.SPINFREE
---------

Use Dyall's spinfree Hamiltonian, Ref. [Dyall1994]_ ,
to obtain results without spin-orbit coupling for the four-component
Hamiltonian in the default restricted kinetic balance scheme. This
keyword works also for two-component relativistic Hamiltonians where one
can choose between two spin-free schemes - see the :ref:`.BSS` 
keyword.

Note that this option should not be used for response calculations with
time-antisymmetric (magnetic) operators as it will eliminate important
contributions. Add more...

.. _.LEVY-LEBLOND:

.LEVY-LEBLOND
-------------

Use the nonrelativistic Lévy-Leblond Hamiltonian [Levy1967]_.

Use this option before any additional one-electron operators are
specified, because it redefines the metric used in the calculation.
Lévy-Leblond provides nonrelativistic mode in DIRAC.

.. _.NONREL:

.NONREL
-------

Standard nonrelativistic calculation based on the Schrodinger equation.
Should give identical energy results as with the :ref:`.LEVY-LEBLOND` keyword.

DIRAC runs in the 2-component spinfree mode, what in fact represents traditional one-component mode.

.. _.DFT:

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

-  LDA exchange  [Dirac1930]_
-  Becke 1988 [Becke1988]_
-  The correction term to LDA proposed by Becke [Becke1988]_
-  Perdew-Wang 1986 [Perdew1986]_

*Correlation functionals:*

-  The Vosko-Wilk-Nusair LDA correlation energy (VWN5) [Vosko1980]_
-  The Lee-Yang-Parr functional  [Lee1988]_
-  Perdew 1986 [Perdew1986]_

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

**Advanced options**

.. _.LVNEW:

.LVNEW
------

Modification of :ref:`.LVCORR` that obtains the atomic small
component charge via a Mulliken analysis instead of the original table
look-up. (The problem with the table look-up is that the electrostatics
in the molecule will be wrong if you have specified a basis set which
does not give the correct small-electron charge because of deficiencies
in the core region.)

.X2COLD
-------

One-step Exact (infinite order) 2-Component relativistic Hamiltonian [Ilias2007]_.

**note:** in the previous release version Dirac2008 this keyword is called .X2C.

DIRAC runs in the (memory saving) 2-component mode. Note that
one should combine this option with the SPINFREE option as X2C will only
provide an unscreened (bare nucleus) spin-orbit operator that gives
unphysically large spin-orbit contributions. To get a realistic screened
spin-orbit operator AMFI is added in the development version unless
:ref:`.NOAMFI` specified. See also the warning there.

.. _.X2C4:

.X2C4
-----

One-step Exact (infinite order) 2-Component relativistic Hamiltonian [Ilias2007]_ .

DIRAC runs in the 4-component mode. This mode is useful if you wish to
restart from a previous 4-component calculation and vice versa. AMFI is
added in the development version unless :ref:`.NOAMFI` specified.
See also the warning there.

.NOSPIN
-------

Implies :ref:`.SPINFREE`, but also remove all spin-symmetry-breaking
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
version this keyword has been moved after :ref:`.DFT`.

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

.. _.INTFLG:

.INTFLG
-------

Specify what two-electron integrals to include. All other modules use
this as the default value.

*Default:* Include LL, SL, and SS integrals (1 = include; 0 = do not
include):

::

    .INTFLG
     1 1 0

.. _.BSS:

.BSS
-----

Use the 2-component relativistic Hamiltonian obtained after the
Barysz-Sadlej-Snijders transformation of the Dirac Hamiltonian in the
finite basis set, see Ref. [Ilias2005]_. Calculations using the 2-component BSS Hamiltonian are running only with
large component basis functions.

For more see tutorial.

.ONESTEP
--------

Together with the :ref:`.BSS` keyword - for the infinite order only
- invokes the one-step infinite order method (which is otherwise called
by :ref:`.X2C`, :ref:`.X2C4` keywords).

.. _.NOAMFI:

.NOAMFI
-------

Do not include the AMFI contribution in the development version where
AMFI is the default. This holds also for keywords :ref:`.X2C` and
:ref:`.X2C4`). In the DIRAC08 distribution version .NOAMFI is the
default.

Warning: missing AMFI contribution (for :ref:`.X2C`,
:ref:`.X2C4` and :ref:`.BSS` keywords) may lead to
overestimation of spin-orbit effects since these would be represented by
one-electron terms only, without two-electron shielding.

.DO2C4C
-------

**Only in the development version of DIRAC**

After iterations at the two-component level *ascend* to the
four-component level.

.. _.DO4C2C:

.DO4C2C
-------

**Only in the development version of DIRAC**

After iterations at the four-component level do the relativistic
transformation to the two-component level.

.USE\_DF
--------

**Only in the development version of DIRAC**

After the four-component DC-SCF method do the infinite order
transformation (either one- or two-step) upon the Fock-Dirac matrix.
Otherwise it is transforming Dirac bare nucleus.

Used only with keyword :ref:`.DO4C2C`.

.. _.CONT2C:

.CONT2C
-------

**Only in the development version of DIRAC**

After four-component DC-SCF continue with two-component iterations.

Used only with keyword :ref:`.DO4C2C`.

Integer (3,4,5) should be specified in free format on the line following
:ref:`.CONT2C`:

::

    .CONT2C
     3

.. _.ZORA:	 

.ZORA
-----

Use the zeroth-order regular approximation [vanLenthe1994]_,[vanLenthe1996]_,[Visscher2000]_,
of the Dirac Hamiltonian in
the Hartree-Fock procedure. Works only for closed-shell systems. The
implementation offers only little computational advantages and is
intended chiefly for comparisons of methodology. Note that the
combination :ref:`.SPINFREE` and :ref:`.ZORA` gives a
spinfree formalism that differs from the conventional spinfree ZORA
formulation. Two integers should be specified in free format on the line
following :ref:`.ZORA`:

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
the ADF code (see the developer's website http://www.scm.com
for further information).

.NOSFMU
-------

**Only in the development version of DIRAC**

In spinfree correlated calculations group multiplication tables are by
default set up as direct products of spatial and spin symmetries. This
flag turns off this, and so the spinfree case is treated similar to the
spin-orbit case.

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

Employ the one-center model as given by :ref:`.ONECAP` until
convergence to a specified THRESHOLD, whereafter the full set of
two-electron integrals will be used. This threshold applies to whatever
convergence criteria has been selected (:ref:`.EVCCNV`, :ref:`.ERGCNV` or :ref:`.FCKCNV`).

.CMPEIG
-------

When some two-component relativistic Hamiltonian is chosen, compare
eigenvalues between the 'parent' four-component Dirac and derived
two-component one-electron Hamiltonians.

For the infinite order (one- and two-step) two-component Hamiltonians
eigenvalues are identical with four-component Dirac counterparts. For
the second-order (and lower order) Douglas-Krol-Hess Hamiltonian they
slightly differ.

.. _\*DFT:

\*DFT
=====

In this link the directives for modifying the DFT calculation are
described (e.g. modifying grid, using ALDA and much more)

.SAOP The more general version of .SAOP! - followed by a single line of
input.

*typical input (here correcting LDA):*

::

    **HAMILTONIAN
    .DFT
    LDA
    *DFT
    .SAOP
    LBalpha

Single line following .SAOP : asymptotic potential (currently LB94 or
LBalpha).


The specification of the functional is described in connection with the
.DFT option in the :ref:`.DFT` section

The DFT module is described in references [Schipper2000]_, [Saue2002]_, [Fossgaard2003]_.

It uses standard non-relativistic functionals. This is justified in many
cases since various studies, like refs. [Varga1999]_, [Varga2000]_, [Mayer1996]_,
indicate that
relativistic corrections to the exchange-correlation functionals have a
negligible effect on spectroscopic constants; however, for example for
core properties further studies are necessary.

Integration algorithm and thresholds
====================================

.SCREENING
----------

Specify the screening threshold used in the BLAS-3 integration scheme.

*Default:*

::

    .SCREENING
     1.0D-18

.POINTWISE
----------

Switch to old (unscreened) BLAS-2 integration scheme (default in
DIRAC10) instead of the new BLAS-3 scheme (default in DIRAC11).

.TINYDENS
---------

Specify the tiny density threshold. Grid points with a density smaller
than this threshold do not contribute to XC matrix elements.

*Default:*

::

    .TINYDENS
     1.0D-14

Asymptotic correction
=====================

.GRAC
-----

Activate the gradient regulated asymptotic correction scheme (GRAC), [Gruning2001]_,
followed by two lines of input.

First line: Asymptotic potential (currently LB94 or LBalpha).

Second line (free format): Parameters α and β (see reference [Gruning2001]_)
the ionization
potential, and the threshold for difference in HOMO eigenvalue below
which asymptotic correction is switched on. A sufficiently converged
density is needed before the correction is activated (before the bulk
potential is shifted in order to reproduce the desired IP (see reference [Gruning2001]_)

*Typical input (here correcting PBE0):*

::

    **HAMILTONIAN
    .DFT
    PBE0
    *DFT
    .GRAC
    LB94
    0.5 40.0 0.79248 1.0D-6

.. _.SAOP!:
	
.SAOP!
------

Activate the statistical averaging of (model) orbital potentials (SAOP)
as defined in reference [Schipper2000]_.
This implies (and
activates) the ALDA kernel. The functional under :ref:`.DFT` is expected to be
GLLBhole. The asymptotic potential is LBalpha.

*Recommended input:*

::

    **HAMILTONIAN
    .DFT
    GLLBhole
    *DFT
    .SAOP!

.SAOP
-----

The more general version of :ref:`.SAOP!`

Followed by a single line of input: asymptotic potential (currently LB94
or LBalpha).

*Typical input (here correcting LDA):*

::

    **HAMILTONIAN
    .DFT
    LDA
    *DFT
    .SAOP
    LBalpha

Spin magnetization TDDFT
========================

.NOSDFT
-------

Turn off spin density contribution to XC response.

.COLLINEAR
----------

Use the collinear approximation as a definition of the spin density.

.. math::

    s = m_z = \sum_i \phi_i^\dagger \Sigma_z \phi_i

instead of the default noncollinear definition


.. math::

    s = \vert \mathbf{m} \vert

.BETASIGMA
----------

use

.. math::

    \mathbf{m} = \sum_{i} \phi_{i}^{\dagger} \beta \mathbf{\Sigma} \phi_{i}

instead of the default form

.. math::

   \mathbf{m} = \sum_i \phi_i^\dagger \mathbf{\Sigma} \phi_i

Adiabatic local density approximation (ALDA)
============================================

.. _.ALDA:

.ALDA
-----

Approximate all functional derivatives beyond the xc potential by SVWN
derivatives. For hybrid functionals exact exchange is switched off in
the solution of the response equation.

.. _.XALDA:

.XALDA
------

Same as :ref:`.ALDA` but keep the fraction of exact exchange of the
xc functional under :ref:`.DFT`.

.ALDA+
------

Use :ref:`.ALDA` only for the Hermitian contribution (density
contribution) and use the proper xc kernel for the anti-Hermitian part
(spin density contribution).

.ALDA-
------

Use :ref:`.ALDA` only for the anti-Hermitian contribution (spin
density contribution) and use the proper xc kernel for the Hermitian
part (density contribution).

.XALDA+
-------

Use :ref:`.XALDA` only for the Hermitian contribution (density
contribution) and use the proper xc kernel for the anti-Hermitian part
(spin density contribution).

.XALDA-
-------

Use :ref:`.XALDA` only for the anti-Hermitian contribution (spin
density contribution) and use the proper xc kernel for the Hermitian
part (density contribution).

Other functionality
===================

**Programmers options**

.GAUNTSCALE
-----------

Scale Gaunt integrals (if included) with the same factor as for
Hartree-Fock exchange. This means that hybrid functionals include
fractional HF Gaunt interaction, and pure functionals no HF Gaunt
interaction at all. If this option is not given the Gaunt integrals will
be included to 100%, meaning that even an LDA calculation will include
full Hartree-Fock Gaunt interaction when the .GAUNT keyword is given in
:ref:`\*\*HAMILTONIAN`.

.. _\*AMFI:

\*AMFI 
======

Specification of the AMFI, RELSCF directives.

**This page describes unreleased functionality. The keywords may not be
available in your version of DIRAC.**

The AMFI program of B.Schimmelpfennig serves for calculating spin-orbit
contributions to various quasirelativistic two-component Hamiltonians.
Its important part is the independent scalar atomic code - RELSCF (named
as AT34), which is providing atomic orbital coefficients for the AMFI
mean-field summation.

The user may control not only the extent of the output, but also provide
here some parameters for the RELSCF code. Note that closed/open-shell
HF-SCF code gives both the nonrelativistic and the spin-free
second-order Douglass-Krol-Hess SCF energies, which can be compared with
DIRAC results.

RELSCF is using only decontrated basis sets restricted to occupied
shells only; for example, if you set some spdf-basis for a DIRAC
BSS+AMFI run of Magnesium (12 electrons), RELSCF will use, for the sake
of computational effectivity, this basis set with only s- and p-
exponents.

**Only in the development version of DIRAC**

.PRNT\_A
--------

Print level for AMFI calculations.

*Default:*

::

    .PRNT_A
     0

.PRNT\_S
--------

Print level for the RELSCF part.

*Default:*

::

    .PRNT_S
     0

.MXITER
-------

Maximum number of iterations in RELSCF.

*Default:*

::

    .MXITER
     50

.AMFICH
-------

Set up an artificial charge of the calculated system which has effect
only on the AMFI mean-field summation(s). This does not change the
electronic occupation for the DIRAC SCF procedure.

The reason is that AMFI-attached scalar relativistic Hartree-Fock SCF
module might not be converging for certain high-spin open-shell atoms
(like Pt,Yb), what is anyway expected for a single-determinant code. A
small positive charge (+1, +2) usually helps to get converged molecular
orbitals of given atom with a negligible change in results due to AMFI
mean-field contribution. The AMFI mean-field contributions are large for
core shells, but small for valence orbitals.

In the case of polyatomic systems the overall AMFI-charge is
proportionally distributed over individual atoms and rounded to the
closest integer value.

Example:

::

    .AMFICH
    +1

The user can check the mean-field occupation/charge of each handled atom
in the AMFI/RELSCF output, together with the RELSCF convergence status.

.. _.NOAMFC:

.NOAMFC
--------

Discard selected centers from AMFI calculations. Related the
:ref:`.BSS` keyword, whose
numerical value *axyz* must have y=0 (i.e. 'spin-free' from the DIRAC
BSS side, while AMFI provides SO terms of first order in potential).

First line specifies how many atoms are neglected, the second line gives
centers to be omitted. Number for each center of calculated system is
determined in the DIRAC output file.

Infinite order scalar terms "from the end", and the (one-electron) spin
orbit term (SO1) from AMFI. Centers 1 and 3 are neglected.

::

    .BSS
    5009
    .NOAMFC
    2
    1  3

.ONLSCF
-------

Stop Dirac calculation after RELSCF program has finished. The
(converged) AO coefficients are written to the file *RELSCF\_COEF* and
can be used to restart from.

This may be useful in combination with the keyword **.AMFICH**, *i.e.*
the user may first run the RELSCF step for an atom with an artificial
charge.

The obtained coefficients can then be used to start the RELSCF step for
the neutral atom for which otherwise one would not obtain convergence.

Input example:

::

    .ONLSCF


Frozen density embedding (FDE) directives.

FDE
===

Frozen density embedding (FDE) directives.

Data import options
-------------------

.EMBPOT
-------

Specifies the name of the file containing the embedding potential.

Extend !
