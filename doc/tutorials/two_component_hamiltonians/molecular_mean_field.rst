

Molecular mean-field X2C
========================

The molecular-mean-field X2C Hamiltonian (X2Cmmf) is **available** for the following methods:

-  Coupled Cluster (CC): module RELCCSD
-  multi-reference (MR) CC: Fock-space CC (including the intermediate Hamiltonian IHFSCC version): module RELCCSD

The remaining post-SCF modules are **not supported** yet:

-  KR-MCSCF
-  KR-CI
-  TD-DFT

The X2Cmmf Hamiltonian scheme can be used for the Dirac-Coulomb as well
as the Dirac-Coulomb-Gaunt Hamiltonian (add .GAUNT to \*\*HAMILTONIAN).

One-step X2Cmmf
---------------

In order to perform a two-component (MR)CC after a four-component SCF calculation add the keyword .X2Cmmf 
to the \*\*HAMILTONIAN keyword section:

::

    **HAMILTONIAN
    .X2Cmmf

After the SCF calculation Dirac will carry out the transformation to two-component prior to the 4-index transformation of
the two-electron integrals and the CC run.

Two-step X2Cmmf
---------------

Besides the one-step appraoch to an X2Cmmf calculation as described above one may also proceed in seperate steps. 
This may be useful if one would like to use less/more MPI processes in the SCF step than in the 4-index/CC calculation. In order to
proceed run the SCF first while adding .X2Cmmf to the \*\*HAMILTONIAN keyword section:

::

    **HAMILTONIAN
    .X2Cmmf

Save the files DFCOEF, AOMOMAT and X2CMAT (use the keywords *--outcmo --get="AOMOMAT X2CMAT"* with our default submit
script *pam*) . Before restarting any calculation to proceed with the 4-index transformation (MOLTRA) and/or the (MR)CC run,
set the keywords combination as indicated below under the \*\*HAMILTONIAN section and add the
following keyword to the namelist RELCCSSD input::

    **HAMILTONIAN
    .X2C
    *X2C
    .mmf-restart
    ...
    ...
    &CCSORT USEOE=T, NORECMP=T &END


or to the new RELCCSD deck input::

    **HAMILTONIAN
    .X2C
    *X2C
    .mmf-restart
    ...
    ...
    *CCSORT
    .NORCMP

Dirac will automatically read the two-component coefficients from the file DFCOEF (use *--incmo --put="AOMOMAT X2CMAT"* on
the *pam* command line) and proceed in two-component mode with the 4-index integral transformation and/or the (MR)CC
calculation, respectively.
