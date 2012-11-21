

.. _hamiltonian:

===============
\*\*HAMILTONIAN
===============
.. index:: **HAMILTONIAN
Specification of the Hamiltonian.

This section defines the electronic Hamiltonian that is to be used.  Internally
the program will always work with 4-component operators that are expanded using
distinct large and small component basis sets. In the transformation to an
orthogonal basis set one may, however, combine large and small component
functions and/or functions of different symmetry in order to obtain a matrix
expansion of e.g. the spin-free modified Dirac equation or the Lévy-Leblond
equation [Visscher2000]_.

In addition one can also modify the Hamiltonian by introduction of an
additional operator, e.g. describing an external field. Any additional operator
defined in this section must be totally symmetric under both the molecular
point group and time reversal symmetry. The latter requirement precludes the
introduction of external magnetic fields.




.. _.LVCORR:

.LVCORR
-------
.. index:: .LVCORR

This keyword activates the Dirac--Coulomb Hamiltonian in which (SS|SS) integrals
are neglected and replaced by an interatomic SS correction (calculated
as a classical repulsion term of (tabulated) small component atomic
charges) [Visscher1997a]_ .

This is currently the most economical and accurate approximation to the
full Dirac--Coulomb Hamiltonian and can certainly be used for the
calculation of spectroscopic constants and valence properties; for core
properties, testing is recommended (see also :ref:`.LVNEW`). This
is the default Hamiltonian choice since DIRAC11.


.. _.GAUNT:

.GAUNT
------
.. index:: .GAUNT
Add the Gaunt interaction to the Hamiltonian. This will increase the
computational time significantly but is important when studying
spin-orbit splittings and/or performing accurate studies of light
molecules. The current implementation is limited to including the Gaunt
interaction in the construction of the Fock matrix and works for
Hartree--Fock and DFT. For using Gaunt in combination with DFT see the
:ref:`\*DFT` section of the
manual. Transformation of the Gaunt part of the two electron operator to
the MO basis is not yet implemented, for this purpose we recommend the
use of the molecular mean field approximation. This can be used for
example in MP2/CC/FSCC/IHFSCC calculations with RELCCSD by means of the
:ref:`.X2Cmmf` Hamiltonian (see also our FAQ/Tutorial pages)

For the relativistic two-component mode (see :ref:`.X2C` keyword)
with AMFI contributions it uses both spin-same and spin-other orbit
mean-field parts.



.. _.DOSSSS:
.DOSSSS
-------
.. index:: .DOSSSS
This keyword gives unmodified Dirac-Coulomb Hamiltonian which was
the default untill the DIRAC11 release. Explicitly including (SS\|SS) type
Coulomb integrals does give the most accurate description of the system
but does increase computational cost significantly. Use this option for
high-accuracy calculations, preferably in conjunction with :ref:`.GAUNT`
to also include the Gaunt correction to the two-electron
interaction.


.. _.X2C:

.X2C
----
.. index:: .X2C

Keyword referring to the new module for the Exact 2-Component
relativistic Hamiltonian, Refs. [Ilias2007]_, [Knecht2010]_.

One should combine this option with an
AMFI correction to the unscreened spin-orbit operator. This is the
default in the development version unless :ref:`.NOAMFI` is specified.


.. _.X2Cmmf:

.X2Cmmf
-------
.. index:: .X2Cmmf

Keyword referring to the 2-component molecular-mean-field Hamiltonian approach
within the new module for the Exact 2-Component relativistic Hamiltonian, Refs.
[Sikkema2009]_, [Knecht2010]_.

DIRAC starts with a 4c-SCF run and performs a transformation to 2-component
mode (based on the converged Fock operator) prior to a post-HF correlation
step. after the SCF. One can combine this option with :ref:`\.GAUNT` which
activates the inclusion of spin-other-orbit contributions in the Hamiltonian.
The X2Cmmf-Hamiltonian can at present only be used for post-HF calculation
within the RELCCSD module. Patches for other correlation modules in DIRAC are
planned to be released soon.

.. _.ECP:

.ECP
----
.. index:: .ECP
.. warning:: documentation missing


.. index:: .PCM
.PCM
----

.. warning:: documentation missing


.. _.SPINFREE:

.SPINFREE
---------
.. index:: .SPINFREE
Use Dyall's spin-free Hamiltonian, Ref. [Dyall1994]_ ,
to obtain results without spin-orbit coupling for the four-component
Hamiltonian in the default restricted kinetic balance scheme. This
keyword works also for two-component relativistic Hamiltonians where one
can choose between two spin-free schemes - see the :ref:`.BSS`
keyword.

Note that this option should not be used for response calculations with
time-antisymmetric (magnetic) operators as it will eliminate important
contributions.


.. _.LEVY-LEBLOND:

.LEVY-LEBLOND
-------------
.. index:: .LEVY-LEBLOND

Use the nonrelativistic Lévy-Leblond Hamiltonian [Levy1967]_.

Use this option before any additional one-electron operators are specified,
because it redefines the metric used in the calculation.


.. _.NONREL:

.NONREL
-------
.. index:: .NONREL
Standard nonrelativistic calculation based on the Schrodinger equation.  Should
give identical energy results as with the :ref:`.LEVY-LEBLOND` keyword.

DIRAC runs in the 2-component spin-free mode, which in fact represents the
traditional one-component mode.


.. _.DFT:
.DFT
----
.. index:: .DFT

Perform a Kohn--Sham density functional theory calculation. In the following line you must specify
the desired DFT functional.

The functional can either be selected from a set of predefined
combinations of exchange and correlation functionals, e.g.::

  .DFT
   B3LYP

or composed by specifying GGAKEY followed by a list of the desired
functionals together with their weights::

  .DFT
   GGAKEY PW86X=1.0 P86C=1.0

*Exchange functionals:*

-  LDA exchange  [Dirac1930]_
-  Becke 1988 [Becke1988]_
-  The correction term to LDA proposed by Becke [Becke1988]_
-  Perdew--Wang 1986 [Perdew1986]_

*Correlation functionals:*

-  The Vosko--Wilk--Nusair LDA correlation energy (VWN5) [Vosko1980]_
-  The Lee--Yang--Parr functional  [Lee1988]_
-  Perdew 1986 [Perdew1986]_

*Predefined combinations:*

-  LDA
-  CAMB3LYP
-  BLYP
-  B3LYP
-  BP86
-  BPW91
-  PBE
-  PBE0

.. _.DFTAUTO:

.DFTAUTO
--------
.. index:: .DFTAUTO
Perform a Kohn--Sham calculation using functionals provided by the XCFun
library. In the following line you must specify the desired DFT
functional.


.. _.LVNEW:

.LVNEW
------
.. index:: .LVNEW

Modification of :ref:`.LVCORR` that obtains the atomic small
component charge via a Mulliken analysis instead of the original table
look-up. (The problem with the table look-up is that the electrostatics
in the molecule will be wrong if you have specified a basis set which
does not give the correct small-electron charge because of deficiencies
in the core region.)

.. _.X2COLD:

.X2COLD
-------
.. index:: .X2COLD
One-step Exact (infinite order) 2-Component relativistic Hamiltonian
[Ilias2007]_.  This keyword was called .X2C prior to DIRAC10.  DIRAC runs in
the (memory saving) 2-component mode. Note that one should combine this option
with the spin-free option as X2C will only provide an unscreened (bare nucleus)
spin-orbit operator that gives unphysically large spin-orbit contributions. To
get a realistic screened spin-orbit operator AMFI is added in the development
version unless :ref:`.NOAMFI` specified.


.. _.X2C4:

.X2C4
-----
.. index:: .X2C4
One-step Exact (infinite order) 2-Component relativistic Hamiltonian
[Ilias2007]_ .

DIRAC runs in the 4-component mode. This mode is useful if you wish to restart
from a previous 4-component calculation and vice versa. AMFI is added in the
development version unless :ref:`.NOAMFI` specified.

.. _.NOSPIN:

.NOSPIN
-------
.. index:: .NOSPIN
Implies :ref:`.SPINFREE`, but also remove all spin-symmetry-breaking
(quaternion "imaginary" or "triplet" terms) from property gradients in response
calculations. Used for analyzing magnetic properties similarly to how it is
done with non-relativistic methods.

.. _.HFXFAC:

.HFXFAC
-------
.. index:: .HFXFAC
Weight of exchange in Fock matrix construction. Default::

  .HFXFAC
   1.0

.. _.OPERATOR:

.OPERATOR
---------
.. index:: .OPERATOR
Specification of an additional one-electron operator in the Hamiltonian.
The operator must be totally symmetric both under the molecular point
group and time reversal symmetry. The field strength of the operator is
specified with COMFACTOR. The keyword can be repeated for addition of
more than one operator.

See the `specification of one-electron
operators </dirac/index.php/Manual:One-electron_operators>`_ section for
more information and explicit examples.

.. warning:: broken link


.. _.INTFLG:

.INTFLG
-------
.. index:: .INTFLG
Specify what two-electron integrals to include. All other modules use
this as the default value. By default
include LL and SL, exclude SS integrals (1 = include; 0 = do not
include)::

  .INTFLG
   1 1 0


.. _.BSS:

.BSS
----
.. index:: .BSS

Use the 2-component relativistic Hamiltonian obtained after the
Barysz--Sadlej--Snijders transformation of the Dirac Hamiltonian in the finite
basis set, see Ref. [Ilias2005]_. Calculations using the 2-component BSS
Hamiltonian are running only with large component basis functions.

.. _.ONESTEP:

.ONESTEP
--------
.. index:: .ONESTEP
Together with the :ref:`.BSS` keyword - for the infinite order only
- invokes the one-step infinite order method (which is otherwise called
by :ref:`.X2C`, :ref:`.X2C4` keywords).


.. _.NOAMFI:

.NOAMFI
-------
.. index:: .NOAMFI
Do not include the AMFI contribution where
AMFI is the default. This holds also for keywords :ref:`.X2C` and
:ref:`.X2C4`). In the DIRAC08 distribution version .NOAMFI was the
default.

.. warning:: missing AMFI contribution (for :ref:`.X2C`,
             :ref:`.X2C4` and :ref:`.BSS` keywords) may lead to
             overestimation of spin-orbit effects since these would be represented by
             one-electron terms only, without two-electron shielding.

.. _.DO2C4C:

.DO2C4C
-------
.. index:: .DO2C4C
After iterations at the two-component level ascend to the
four-component level.

.. warning :: only in development version


.. _.DO4C2C:

.DO4C2C
-------
.. index:: .DO4C2C
After iterations at the four-component level do the relativistic
transformation to the two-component level.

.. warning :: only in development version

.. _.USE\_DF:

.USE\_DF
--------
.. index:: .USE_DF
.. warning :: only in development version

After the four-component DC-SCF method do the infinite order transformation
(either one- or two-step) upon the Fock-Dirac matrix.  Otherwise it is
transforming Dirac bare nucleus.

Used only with keyword :ref:`.DO4C2C`.


.. _.CONT2C:

.CONT2C
-------
.. index:: .CONT2C
.. warning :: only in development version

After four-component DC-SCF continue with two-component iterations.

Used only with keyword :ref:`.DO4C2C`.

Integer (3,4,5) should be specified in free format on the line following
:ref:`.CONT2C`::

  .CONT2C
   3

.. warning:: and what does this integer mean?


.. _.ZORA:

.ZORA
-----
.. index:: .ZORA
Use the zeroth-order regular approximation [vanLenthe1994]_,[vanLenthe1996]_,[Visscher2000]_,
of the Dirac Hamiltonian in
the Hartree-Fock procedure. Works only for closed-shell systems. The
implementation offers only little computational advantages and is
intended chiefly for comparisons of methodology. Note that the
combination :ref:`.SPINFREE` and :ref:`.ZORA` gives a
spin-free formalism that differs from the conventional spin-free ZORA
formulation. Two integers should be specified in free format on the line
following :ref:`.ZORA`::

  .ZORA
   1 1

The first number indicates whether the density is to be normalized over
the 2-component (0; ZORA) or 4-component metric (1; ZORA4).

The second number specifies whether the orbital energies should be
unmodified (0; normal ZORA) or scaled (1; scaled ZORA).

.. _.SOLVENT:

.SOLVENT
--------
.. index:: .SOLVENT
Model solvent effects by placing the molecule in a spherical cavity in a
dielectric continuum.

.. _.FDE:

.FDE
----
.. index:: .FDE
Activates the frozen density embedding (FDE) functionality. Options can
be specified under the \*FDE menu.

In order to use FDE, the user must generated an embedding potential with
the ADF code (see the developer's website http://www.scm.com
for further information).

.. _.NOSFMU:

.NOSFMU
-------
.. index:: .NOSFMU
In spin-free correlated calculations group multiplication tables are by
default set up as direct products of spatial and spin symmetries. This
flag turns off this, and so the spin-free case is treated similar to the
spin-orbit case.

.. warning :: only in development version

.. _.PRINT:

.PRINT
------
.. index:: .PRINT
Print level. Default::

  .PRINT
   0

.. _.ONESYS:

.ONESYS
-------
.. index:: .ONESYS
Ignore two-particle interactions.

The keyword ensures the diagonalization of the bare nuclei Hamiltonian
matrix without proceeding to the iterative SCF method.

.. _.URKBAL:

.URKBAL
-------
.. index:: .URKBAL
Unrestricted kinetic balance.

The default is restricted kinetic balance. This is imposed by deleting
unphysical solutions from the free particle positronic spectrum. This
leads to a 1:1 ratio of electronic and positronic solutions. This
preprojection is sensitive to linear dependencies and should therefore
preferably be used in conjunction with the spherical transformation of
both large and small components.

.. _.FREEPJ:

.FREEPJ
-------
.. index:: .FREEPJ
Project out all free positronic solutions from the MO space.

.. _.VEXTPJ:

.VEXTPJ
-------
.. index:: .VEXTPJ
Project out all external field positronic solutions from the MO space.

.. _.NOSMLV:

.NOSMLV
-------
.. index:: .NOSMLV
Delete SS nuclear attraction integrals. This will take out contributions
to the one-electron spin-orbit interaction and the Darwin interaction.

.. _.SMLV1C:

.SMLV1C
-------
.. index:: .SMLV1C
Neglect potential for multi-center SS blocks, i.e. multi-center SS
nuclear attraction integrals and multi-center SS two-electron integrals.

.. _.ONECAP:

.ONECAP
-------
.. index:: .ONECAP
Consider taking the `.SMLV1C <#.SMLV1C>`_ model one step further. Only
one-center contributions to the LS and SS two-electron integrals and SS
nuclear attraction integrals are calculated explicitly. The
electrostatic effects of the terms neglected this way are included by
calculating the classical repulsion from small component charges based
on a Mulliken population analysis. Note that we therefore only need to
calculate the derivative of the LL integrals when calculating the
molecular gradient.

.. warning:: broken link

.. _.ONECNV:

.ONECNV
-------
.. index:: .ONECNV
Employ the one-center model as given by :ref:`.ONECAP` until
convergence to a specified THRESHOLD, whereafter the full set of
two-electron integrals will be used::

  .ONECNV
   THRESHOLD

This threshold applies to whatever
convergence criteria has been selected (:ref:`.EVCCNV`, :ref:`.ERGCNV` or :ref:`.FCKCNV`).

.. _.CMPEIG:

.CMPEIG
-------
.. index:: .CMPEIG
When some two-component relativistic Hamiltonian is chosen, compare
eigenvalues between the 'parent' four-component Dirac and derived
two-component one-electron Hamiltonians.

For the infinite order (one- and two-step) two-component Hamiltonians
eigenvalues are identical with four-component Dirac counterparts. For
the second-order (and lower order) Douglas--Kroll--Hess Hamiltonian they
slightly differ.

.. _.SCQSET:

.SCQSET
-------
.. index:: .SCQSET
.. warning:: documentation missing

.. _.SPINF2:

.SPINF2
-------
.. index:: .SPINF2
.. warning:: documentation missing

.. _.HFXMU:

.HFXMU
------
.. index:: .HFXMU
.. warning:: documentation missing

.. _.YREQ1:

.YREQ1
------
.. index:: .YREQ1
.. warning:: documentation missing

.. _.BLOCKD:

.BLOCKD
-------
.. index:: .BLOCKD
.. warning:: documentation missing

.. _.BEG_2C:

.BEG_2C
-------
.. index:: .BEG_2C
.. warning:: documentation missing

.. _.VACPOL:

.VACPOL
-------
.. index:: .VACPOL
.. warning:: documentation missing

.. _.MIXVAC:

.MIXVAC
-------
.. index:: .MIXVAC
.. warning:: documentation missing

.. _.CAP:

.CAP
----
.. index:: .CAP
.. warning:: documentation missing

.. _.NOVREF:

.NOVREF
-------
.. index:: .NOVREF
.. warning:: documentation missing

.. _.POSVAC:

.POSVAC
-------
.. index:: .POSVAC
.. warning:: documentation missing

.. _.IOTC4:

.IOTC4
------
.. index:: .IOTC4
.. warning:: documentation missing

.._.BSS4:

.BSS4
-----
.. index:: .BSS4
.. warning:: documentation missing

.. _.DKH2:

.DKH2
-----
.. index:: .DKH2
.. warning:: documentation missing

.. _.DKH1:

.DKH1
-----
.. index:: .DKH1
.. warning:: documentation missing

.. _.QDOTS:

.QDOTS
------
.. index:: .QDOTS
.. warning:: documentation missing

.. _.ECP:

.ECP
----
.. index:: .IOTC

.. warning:: documentation missing

.. _.MMF:

.MMF
----
.. index:: .MMF
.. warning:: documentation missing

.. _.MO4C2C:

.MO4C2C
-------
.. index:: .MO4C2C
.. warning:: documentation missing

.. _.HFXATT:

.HFXATT
-------
.. index:: .HFXATT
.. warning:: documentation missing

.. _.xB:

.xB
---
.. index:: .xB

.. warning:: documentation missing

.. _.xC:

.xC
---
.. index:: .xC

.. warning:: documentation missing
