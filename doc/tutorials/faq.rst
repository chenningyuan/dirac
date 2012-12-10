I would like to use the 2c-molecular-mean-field Hamiltonian (X2Cmmf) for my post-HF (Coupled-Cluster [CC], Fock-Space CC, KR-MCSCF or KR-CI) calculations - what do i need to consider?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The 2c-molecular-mean-field Hamiltonian (X2Cmmf) is at present (Dirac11
release Nov 2011) **available** for the following methods:

-  Coupled-Cluster
-  Fock-space Coupled Cluster (including the intermediate Hamiltonian
   version)

The remaining modules are **not supported** at the Dirac11 release date
but patches are planned in the course of 2012:

-  KR-MCSCF
-  KR-CI
-  TD-DFT

The X2Cmmf Hamiltonian scheme can be used for the Dirac-Coulomb as well
as the Dirac-Coulomb-Gaunt Hamiltonian (add .GAUNT to \*\*HAMILTONIAN).

As follows we will consider the necessary inputs for an X2Cmmf
calculation where we perform the SCF + post-SCF calculation in **one**
shot or **two separate** runs:

-  **one-step** procedure:

   -  add .X2Cmmf to the \*\*HAMILTONIAN keyword section:

::

    **HAMILTONIAN
    .X2Cmmf

that's it. After the SCF calculation the program will proceed with the
transformation to 2-components prior to the post-HF steps.

-  **two-step** procedure (**EXPERT mode**!):

1. ***SCF step***. Add .X2Cmmf to the \*\*HAMILTONIAN keyword section:

::

    **HAMILTONIAN
    .X2Cmmf

and save the file DFCOEF (keyword *--outcmo* with our default submit
script *pam*)

2. ***post-HF step***. Remove .X2Cmmf from the \*\*HAMILTONIAN section
(**non-intuitive (!)** step which will be solved in Dirac12) and add the
following keywords to the \*\*MOLTRA and namelist RELCCSSD input
section.

::

    **MOLTRA
    .INTFL2
    1 0 0 0 
    .INTFL4
    1 0 0 0   
    ...
    ...
    &CCSORT USEOE=T, NORECMP=T &END

You are now ready to restart your post-HF calculation but do not forget
to copy your MO-coefficient file DFCOEF to your scratch directory, i.e.
add *--incmo* to your *pam* command line.

