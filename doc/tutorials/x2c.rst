Run and restart an X2C calculation
==================================

These examples shows how one can run a Hartree-Fock (HF) calculation in
2-component mode with the exact two component (X2C) Hamiltonian and
restart a property calculation at the HF level in a second run without
the need to re-do the X2C transformation.

The Br2 molecule (`Br2.xyz <#Br2.xyz>`_) has linear symmetry (D\_oo,h),
and the program will detect this automatically. In most parts of the
program the highest abelian subgroup (D2h) will be used instead of the
full linear symmetry. In this example we use Dyall's augmented TZ basis
dyall.av3z in uncontracted form.

The file (`scf\_x2c.inp <#scf_x2c.inp>`_) provides an input for the
initial HF run. Prior to the HF the Hamiltonian transformation from 4-
to 2-components will be performed and for restart purposes one should
therefore save the files **X2CMAT** and **AOMOMAT** containing the
Hamiltonian information and the file **DFCOEF** containing the
(converged) MO-coefficients (keyword **--outcmo** below).

Run the first step in the example using:

::

    pam --mol=Br2.xyz --inp=scf_x2c.inp --outcmo --get "X2CMAT AOMOMAT" 

Output will be written to the file scf\_x2c\_Br2.out, which can be
monitored during the calculation.

scf\_x2c.inp
''''''''''''

::

    **DIRAC
    .WAVE FUNCTION
    **MOLECULE
    *BASIS
    ! Use the same basis for all atoms
    .DEFAULT
    dyall.av3z
    **HAMILTONIAN
    ! Use the exact two-component Hamiltonian
    ! with atomic mean field spin-orbit
    .X2C
    **INTEGRALS
    *READIN
    ! A limitation in the AMFI code requires an uncontracted basis
    .UNCONTRACT
    **WAVE FUNCTION
    .SCF
    *END OF

Br2.xyz
'''''''

::

    2
    Br 0.0 0.0  1.244
    Br 0.0 0.0 -1.244

In order to proceed in a second step with a property calculation of,
e.g. the contact density at the Br nucleus, at the HF level, we consider
the file (`property\_x2c.inp <#property_x2c.inp>`_). Restarting from the
(converged) MO-coefficients and reading the transformed X2C-Hamiltonian
from file is done automatically once the file are present. We therefore
restart the calculation using:

::

    pam --mol=Br2.xyz --inp=property_x2c.inp --incmo --copy "X2CMAT AOMOMAT" 

Output will be written to the file property\_x2c\_Br2.out, which can be
monitored during the calculation.

property\_x2c.inp
'''''''''''''''''

::

    **DIRAC
    .PROPERTIES
    **MOLECULE
    *BASIS
    ! Use the same basis for all atoms
    .DEFAULT
    dyall.av3z
    **HAMILTONIAN
    ! Use the exact two-component Hamiltonian
    ! with atomic mean field spin-orbit
    .X2C
    **INTEGRALS
    *READIN
    ! A limitation in the AMFI code requires an uncontracted basis
    .UNCONTRACT
    **PROPERTIES
    .RHONUC
    *END OF

