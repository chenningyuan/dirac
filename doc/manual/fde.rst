star(FDE)

Frozen density embedding (FDE) directives. 

information on the FDE formalism can be found in [Gomes2012a]_ and references therein,
whereas details concerning the Dirac implementation are found in [Gomes2008]_, 
[Hofener2012]_ and [Hofener2013]_. 

keyword(EMBPOT)

Specifies the name of the file containing the embedding potential over an integration grid.

*Default:*

::

    .EMBPOT
     EMBPOT


By specifying this option, the program will compute matrix elements of this embedding potential on
the AO basis used by the program, and add it to the one-electron Fock matrix. Because of that, 
this option can be used with any wavefunction in Dirac. 

If only this option is specified, the calculation corresponds to the "fixed potential" scheme of [Gomes2008]_. 

keyword(FRZDNS)

Specifies the name of the file containing the "frozen" quantities (integration grid, electrostatic potential, 
electrostatic potential due to external fields (e.g. point charges), density and (x,y,z) components of 
the density gradient). 

*Default:*

::

    .FRZDNS
     FRZDNS

If :ref:`.EMBPOT` is not specified, the program implicitly assums this option has been selected and
proceeds to the calculation of the embedding potential during the SCF procedure. 

In doing so, the non-additive parts will be calculted with the orbital-free kinetic energy and 
exchange-correlation density functionals specified under :ref:`.NAKEF` and :ref:`.NAXCF`, 
respectively. 

This option does not yet have any effect for correlated wavefunctions (e.g. MP2, CI, CC).

keyword(NAKEF)

Specifies the orbital-free kinetic energy density functional used to calculate the non-additive
kinetic energy contributions to the ground-state embedding potential (see ref:`.FRZDNS`) and/or
response contributions (see :ref:`.RSP` or :ref:`.RSPLDA`)

*Default:*

::

    .NAKEF
     kin_tf


Which corresponds to the Thomas-Fermi kinetic energy functional.
Other available functionals are: PW91k (kin_pw91).


keyword(NAXCF)

Specifies the exchange-correlation density functional used to calculate the non-additive
kinetic energy contributions to the ground-state embedding potential (see ref:`.FRZDNS`)
and/or response contributions (see :ref:`.RSP` or :ref:`.RSPLDA`)

*Default:*

::

    .NAXCF
     lda


Where lda is equivalent to specifying "slaterx=1.0 vwnc=1.0" (without quotes). 
Other available functionals are: pbe (equivalent to "pbex=1.0 pbec=1.0"), 
blyp (equivalent to "beckex=1.0 lypc=1.0"), pp86 (equivalent to "pw86x=1.0 p86c=1.0") 

keyword(RSP)

Specifies that FDE response contributions should be calculated, as this is not done by default. 
The density functionals used will be those specified under :ref:`.NAKEF` and :ref:`.NAXCF`, if any. 

One should not that this option can be used together with either :ref:`.FRZDNS` or :ref:`.EMBPOT`, 
but it requires that "frozen" quantities are given (as in :ref:`.FRZDNS`) 

This option does not yet have any effect for correlated wavefunctions (e.g. MP2, CI, CC).


keyword(RSPLDA)

Specifies that the non-additive exchange-correlation and kinetic energy 
FDE response contributions should be calculated with the LDA and Thomas-Fermi 
functionals, respectively.

This option does not yet have any effect for correlated wavefunctions (e.g. MP2, CI, CC).
