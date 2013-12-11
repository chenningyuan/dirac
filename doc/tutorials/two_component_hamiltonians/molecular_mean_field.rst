

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

Besides the simple, recommended, one-step appraoch to a X2Cmmf calculation as described above one may also proceed in seperate steps. 
This may be useful if one would like to use less/more MPI processes in the SCF step than in the 4-index/CC calculation. In order to
proceed run the SCF first while adding .X2Cmmf to the \*\*HAMILTONIAN keyword section:

::

    **HAMILTONIAN
    .X2Cmmf

Save the file DFCOEF (keyword *--outcmo* with our default submit
script *pam*). Before restarting calculation to proceed with the 4-index transformation (MOLTRA) and the (MR)CC run,
remove the keyword .X2Cmmf from the \*\*HAMILTONIAN section (**seems non-intuitive (!)**) and add the
following keywords to the \*\*MOLTRA and either the namelist RELCCSSD input::

    **HAMILTONIAN
    **MOLTRA
    .INTFL2
    1 0 0 0
    .INTFL4
    1 0 0 0
    ...
    ...
    &CCSORT USEOE=T, NORECMP=T &END


or the new RELCCSD deck input::

    **MOLTRA
    .INTFL2
    1 0 0 0
    .INTFL4
    1 0 0 0
    **RELCC
    ...
    ...
    *CCSORT
    .NORCMP

Dirac will automatically read the two-component coefficients from the file DFCOEF (to add it to the run use *--incmo* on
the *pam* command line) and proceed in two-component mode with the 4-index integral transformation and (MR)CC
calculation, respectively.
