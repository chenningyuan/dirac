:orphan:
 

star(FDE)

Frozen density embedding (FDE) directives. 

General information on the FDE method can be found in recent 
reviews (e.g. :cite:`Jacob2013`, :cite:`Gomes2012a`), whereas details 
concerning the Dirac implementation are found in :cite:`Gomes2008`, 
:cite:`Hofener2012` and :cite:`Hofener2013`. 

**Data import**
===============

If Dirac is used to calculate the subsystem of interest, data from
the subsystem(s) making up the environment must be fed into the program.

keyword(EMBPOT)

Specifies the name of the file containing the embedding potential over an integration grid.

*Default:*

::

    .EMBPOT
     EMBPOT


keyword(FRDENS)

Specifies the name of the file containing the "frozen" quantities (integration grid, electrostatic potential, 
electrostatic potential due to external fields (e.g. point charges), density and (x,y,z) components of 
the density gradient). 

*Default:*

::

    .FRDENS
     FRZDNS


**Embedding potential**
=======================

If the :ref:`FDE_.EMBPOT` option was selected, the program will simple compute matrix elements of the
the AO basis used by the program, and add it to the one-electron Fock matrix. Because of that,
this option can be used with any wavefunction in Dirac. The ground-state calculation in this 
case corresponds to the "fixed potential" scheme of :cite:`Gomes2008`.

However, if :ref:`FDE_.EMBPOT` is not specified, the program implicitly selects the 
:ref:`FDE_.UPDATE` option, and assumes a frozen density has been provided (see :ref:`FDE_.FRDENS`).

keyword(UPDATE)

This option requires that the embedding potential for the ground-state be calculated using the 
"frozen" quantities and those for the subsystem treated with Dirac. The non-additive 
contributions are calculted with the orbital-free kinetic energy and exchange-correlation 
density functionals specified under :ref:`FDE_.NAKEF` and :ref:`FDE_.NAXCF`, respectively. 

Currently this option is only supported during the SCF procedure, so it does not have any 
effect for correlated wavefunctions (e.g. MP2, CI, CC). 

The :ref:`FDE_.EMBPOT` and this option cannot be specified at the same time.


keyword(NAKEF)

Specifies the orbital-free kinetic energy density functional used to calculate the non-additive
kinetic energy contributions to the ground-state embedding potential (see :ref:`FDE_.UPDATE`) and/or
response contributions (see :ref:`FDE_.RSP` or :ref:`FDE_.RSPLDA`)

*Default:*

::

    .NAKEF
     kin_tf


Which corresponds to the Thomas-Fermi kinetic energy functional.
Other available functionals are: PW91k (kin_pw91).


keyword(NAXCF)

Specifies the exchange-correlation density functional used to calculate the non-additive
kinetic energy contributions to the ground-state embedding potential (see :ref:`FDE_.UPDATE`)
and/or response contributions (see :ref:`FDE_.RSP` or :ref:`FDE_.RSPLDA`)

*Default:*

::

    .NAXCF
     lda


Where lda is equivalent to specifying "slaterx=1.0 vwnc=1.0" (without quotes). 
Other available functionals are: pbe (equivalent to "pbex=1.0 pbec=1.0"), 
blyp (equivalent to "beckex=1.0 lypc=1.0"), pp86 (equivalent to "pw86x=1.0 p86c=1.0") 

**Response contributions**
==========================

For response-based approaches (TDHF, TDDFT), contributions from FDE to the active
subsystem can be included through the keywords below (the default is to not include
any). For correlated wavefunctions (e.g. MP2, CI, CC) these do not yet have any effect. 

These options can be used together with either :ref:`FDE_.UPDATE` or :ref:`FDE_.EMBPOT`, but 
require that "frozen" quantities are present (see :ref:`FDE_.FRDENS`) 


keyword(RSP)

Specifies that FDE response contributions should be calculated employing
the density functionals specified under :ref:`FDE_.NAKEF` and :ref:`FDE_.NAXCF`, if any. 

keyword(RSPLDA)

Specifies that the non-additive exchange-correlation and kinetic energy 
FDE response contributions should be calculated with the LDA and Thomas-Fermi 
functionals, respectively.


**Data export**
===============

Dirac can also provide ground-state data (density and density gradient, electrostatic potential etc) 
from the subsystem in question over a grid, so that these can be used by other 
programs.

keyword(GRIDOU)

.. warning:: Development version only. 

Specifies the grid over which the quantities will be calculated and exported, and enables
the export. The input is a XYZW file, and the output is in XML format. The original file
is overwritten.

*Default:*

::

    .GRIDOU
     GRIDOUT 


keyword(LEVEL)

.. warning:: Development version only. 

Specifies which kind of wavefunction will be used to obtain the exported quantities

*Default:*

::

    .LEVEL
     DHF 


which corresponds to SCF (Hartree-Fock, DFT) ones. Also supported: MP2.

keyword(EXONLY)

.. warning:: Development version only. 

This option forces the program to skip the calculation of any FDE contributions, and proceed to
exporting the 

**Debug/expert options**
========================

keyword(SKIPX)

Specifies that the non-additive exchange-correlation contributions are not to be calculated. 

keyword(SKIPK)

Specifies that the non-additive kinetic energy contributions are not to be calculated. 
