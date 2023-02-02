:orphan:

xamfX2C: Two-electon picture change corrections for the X2C Hamiltonian
=======================================================================

This tutorial demonstrates the usage and usefulness of the (extended) atomic mean-field 
two-electron picture-change corrections for the X2C Hamiltonian. 
For more information on the actual details of the theory and implementation 
see Ref. :cite:`Knecht2022`.

Methane (:math:`{\rm CH}_4`) - The ultrarelativistic case  
---------------------------------------------------------

In the following example, we will calculate total energies and spinor energies of 
CH :math:`_4` with several different Hamiltonians within the framework of DFT/PBE/dyall.v2z. 

For convenience, in the Table below (see also Table IX in Ref. :cite:`Knecht2022`), 
we summarise the SCF total energy (E) and spinor energies (:math:`\epsilon`) of the doubly 
degenerate occupied spinors for :math:`{\rm CH}_4` as obtained from DFT/PBE/v2z calculations within 
a four-component Dirac–Coulomb (:math:`^4{\rm DC}`) as well as a two-component Hamiltonian framework, 
including the :math:`{\rm (e)amfX2C}_{\rm DC}` models. All energies are given in Hartree. 
The speed of light :math:`c` was reduced by a factor 10.

+-----------------------+-----------------------------+------------------------------+------------------------------+-------------------------------+-----------------------+
|                       | :math:`{\rm 1eX2C}_{\rm D}` |:math:`{\rm AMFIX2C}_{\rm D}` | :math:`{\rm amfX2C}_{\rm DC}`| :math:`{\rm eamfX2C}_{\rm DC}`| :math:`{}^{4}{\rm DC}`|
+=======================+=============================+==============================+==============================+===============================+=======================+
| E                     |     -42.14220               |    -42.14039                 |     -42.26469                |    -42.25774                  |    -42.25850          |
+-----------------------+-----------------------------+------------------------------+------------------------------+-------------------------------+-----------------------+ 
|:math:`\epsilon_1`     |     -10.22206               |    -10.22401                 |     -10.36154                |    -10.36028                  |    -10.35794          |
+-----------------------+-----------------------------+------------------------------+------------------------------+-------------------------------+-----------------------+ 
|:math:`\epsilon_2`     |     -0.65939                |     -0.65966                 |     -0.66288                 |    -0.66269                   |    -0.66274           |
+-----------------------+-----------------------------+------------------------------+------------------------------+-------------------------------+-----------------------+
|:math:`\epsilon_3`     |     -0.35705                |     -0.35288                 |     -0.35649                 |    -0.35291                   |    -0.35290           |
+-----------------------+-----------------------------+------------------------------+------------------------------+-------------------------------+-----------------------+
|:math:`\epsilon_{4-5}` |     -0.33313                |     -0.33503                 |     -0.33282                 |    -0.33527                   |    -0.33544           |
+-----------------------+-----------------------------+------------------------------+------------------------------+-------------------------------+-----------------------+

In order to make use of the (e)amf two-electron picture-change corrections, we will need to run atomic 
calculations for each element / nuclei type in the molecule. Hence, in the present case, 
we will run atomic calculations for C and H, respectively. The latter seems odd as H is by definition a 
one-electron system but as we will see, the MO coefficients will be needed for the :math:`{\rm eamfX2C}_{\rm DC}` 
model.

Discussion
++++++++++

For a detailed discussion of the performance of the various picture-change correction models for the X2C Hamiltonian 
in comparison with the four-component reference data, we refer the interested user to Ref. :cite:`Knecht2022`. 
Let us highlight, though, that for the particular case of :math:`{\rm eamfX2C}_{\rm DC}` all 
picture-change correction terms for DFT (and likewise HF) 
are derived in molecular basis on the basis of molecular densities which are built 
from a superposition of atomic input densities. 
Consequently, in contrast to the simpler :math:`{\rm amfX2C}_{\rm DC}` model, the former comprises 
atomic contributions regardless of the actual atom type, that is, both “light” (one-electron) and
“heavy” (many-electron) atoms contribute equally. 
Bearing the latter in mind, the total SCF energies as well as
spinor energies compiled in the above Table confirm
the unique numerical performance of the :math:`{\rm eamfX2C}_{\rm DC}` 
Hamiltonian model in comparison to the :math:`{}^4{\rm DC}` reference data.

Atomic calculations
-------------------

Note that, besides scaling the speed of light by a factor of :math:`c`/10, we also 
use a point nucleus as well as a different set of constants as default in our calculations. 
This was done for Ref. :cite:`Knecht2022` to match calculations with `ReSpect <http://www.respectprogram.org/>`__. 
Since we would like to reproduce in this tutorial the results from Ref. :cite:`Knecht2022`, we will keep 
those details. Needless to say, the :math:`{\rm eamf}` model works for finite nuclei just as for point nuclei.

1. H atom

For the atomic run of H, we use the following input (``h.inp``)

:: 
   
   **DIRAC
   .WAVE FUNCTION
   **INTEGRAL
   .NUCMOD
    1
   *READIN
   .UNCONTRACTED
   **GENERAL
   .CVALUE
   13.703599907400d0
   .CODATA
    CODATA10
   **HAMILTONIAN
   .X2Cmmf
   .DOSSSS
   .DFT
    PBE
   .ONESYS
   *X-AMFI
   .amf
   .eamf
   **WAVE FUNCTION
   .SCF
   *SCF
   **MOLECULE
   *BASIS
   .DEFAULT
    dyall.v2z
   **END OF

and the xyz file (``h.xyz``)::

   1
   atom xyz file
   H          0.0000000000   0.0000000000   0.0000000000

and run the calculation with::

   $ ./pam --inp=h.inp --mol=h.xyz --get="amfPCEC.h5" && mv amfPCEC.h5 amfPCEC.001.h5

In the latter step, we rename the file ``amfPCEC.h5`` for later convenience.

2. C atom

For the atomic run of C, we use the following input (``c.inp``)

:: 
   
   **DIRAC
   .WAVE FUNCTION
   **INTEGRAL
   .NUCMOD
    1
   *READIN
   .UNCONTRACTED
   **GENERAL
   .CVALUE
   13.703599907400d0
   .CODATA
    CODATA10
   **HAMILTONIAN
   .X2Cmmf
   .DOSSSS
   .DFT
    PBE
   *X-AMFI
   .amf
   .eamf
   **WAVE FUNCTION
   .SCF
   *SCF
   .CLOSED SHELL
   4 0
   .OPEN SHELL
   1
   2/0,6
   **MOLECULE
   *BASIS
   .DEFAULT
    dyall.v2z
   **END OF

and the xyz file (``c.xyz``)::

   1
   atom xyz file
   C          0.0000000000   0.0000000000   0.0000000000

and run the calculation with::

   $ ./pam --inp=c.inp --mol=c.xyz --get="amfPCEC.h5" && mv amfPCEC.h5 amfPCEC.006.h5

In the latter step, we rename the file ``amfPCEC.h5`` for later convenience.

Molecular calculations
----------------------

Before embarking on the molecular calculations, we need to make a single preparatory step, that is::

   $ $BUILD_DIR/merge_amf.py amfPCEC.006.h5 amfPCEC.001.h5

which will combine the **atomic** ``amfPCEC.xxx.h5`` files into a single ``amfPCEC.h5`` file 
ready for use in the molecular :math:`{\rm (e)amfX2C}` calculations below. The python 
script ``merge_amf.py`` is located in the build directory (named in the above code line 
``$BUILD_DIR``) of DIRAC.

The xyz coordinates for :math:`{\rm CH}_4` are given below (``ch4.xyz``)::

   5

   C      0.0000000000   0.0000000000   0.0000000000
   H      0.6298891440   0.6298891440  -0.6298891440
   H      0.6298891440  -0.6298891440   0.6298891440
   H     -0.6298891440   0.6298891440   0.6298891440
   H     -0.6298891440  -0.6298891440  -0.6298891440 


1. :math:`{\rm 1eX2C}_{\rm D}`

To run the molecular calculation for :math:`{\rm CH}_4` with the :math:`{\rm 1eX2C}_{\rm D}` Hamiltonian, 
that is, neglecting any two-electron picture-change corrections, we use the following input (``1eX2C.inp``)::

   **DIRAC
   .WAVE FUNCTION
   **GENERAL
   .CVALUE
   13.703599907400d0
   .CODATA
    CODATA10
   **INTEGRAL
   .NUCMOD
    1
   *READIN
   .UNCONTRACTED
   **HAMILTONIAN
   .X2C
   .DFT
    PBE
   .NOAMFI
   **WAVE FUNCTION
   .SCF
   *SCF
   .CLOSED SHELL
    10
   **MOLECULE
   *BASIS
   .DEFAULT
    dyall.v2z
   **END OF

The calculation can be run with::

   $ ./pam --inp=1eX2C.inp --mol=ch4.xyz

2. :math:`{\rm AMFIX2C}_{\rm D}`

To run the molecular calculation for :math:`{\rm CH}_4` with the :math:`{\rm AMFIX2C}_{\rm D}` Hamiltonian, 
that is, with two-electron picture-change corrections from the \*AMFI module, 
we use the following input (``AMFIX2C.inp``)::

   **DIRAC
   .WAVE FUNCTION
   **GENERAL
   .CVALUE
   13.703599907400d0
   .CODATA
    CODATA10
   **INTEGRAL
   .NUCMOD
    1
   *READIN
   .UNCONTRACTED
   **HAMILTONIAN
   .X2C
   .DFT
    PBE
   **WAVE FUNCTION
   .SCF
   *SCF
   .CLOSED SHELL
    10
   **MOLECULE
   *BASIS
   .DEFAULT
    dyall.v2z
   **END OF

The calculation can be run with::

   $ ./pam --inp=AMFIX2C.inp --mol=ch4.xyz

3. :math:`{\rm amfX2C}_{\rm DC}`

To run the molecular calculation for :math:`{\rm CH}_4` with the :math:`{\rm amfX2C}_{\rm DC}` Hamiltonian, 
that is, with atomic-mean-field two-electron picture-change corrections, 
we use the following input (``amfX2C.inp``)::

   **DIRAC
   .WAVE FUNCTION
   **GENERAL
   .CVALUE
   13.703599907400d0
   .CODATA
    CODATA10
   **INTEGRAL
   .NUCMOD
    1
   *READIN
   .UNCONTRACTED
   **HAMILTONIAN
   .X2C
   .DFT
    PBE
   *X-AMFI
   .amf
   **WAVE FUNCTION
   .SCF
   *SCF
   .CLOSED SHELL
    10
   **MOLECULE
   *BASIS
   .DEFAULT
    dyall.v2z
   **END OF

The calculation can be run with::

   $ ./pam --inp=amfX2C.inp --mol=ch4.xyz --put="amfPCEC.h5"

4. :math:`{\rm eamfX2C}_{\rm DC}`

To run the molecular calculation for :math:`{\rm CH}_4` with the :math:`{\rm eamfX2C}_{\rm DC}` Hamiltonian, 
that is, with full extended atomic-mean-field two-electron picture-change corrections, 
we use the following input (``eamfX2C.inp``)::

   **DIRAC
   .WAVE FUNCTION
   **GENERAL
   .CVALUE
   13.703599907400d0
   .CODATA
    CODATA10
   **INTEGRAL
   .NUCMOD
    1
   *READIN
   .UNCONTRACTED
   **HAMILTONIAN
   .X2C
   .DFT
    PBE
   *X-AMFI
   .eamf
   **WAVE FUNCTION
   .SCF
   *SCF
   .CLOSED SHELL
    10
   **MOLECULE
   *BASIS
   .DEFAULT
    dyall.v2z
   **END OF

The calculation can be run with::

   $ ./pam --inp=eamfX2C.inp --mol=ch4.xyz --put="amfPCEC.h5"

5. :math:`^{4}{\rm DC}`

To run the molecular calculation for :math:`{\rm CH}_4` with the :math:`{}^4{\rm DC}` Hamiltonian, 
that is, in full 4-component mode,  we use the following input (``4DC.inp``)::

   **DIRAC
   .WAVE FUNCTION
   **GENERAL
   .CVALUE
   13.703599907400d0
   .CODATA
    CODATA10
   **INTEGRAL
   .NUCMOD
    1
   *READIN
   .UNCONTRACTED
   **HAMILTONIAN
   .DOSSSS
   .DFT
    PBE
   **WAVE FUNCTION
   .SCF
   *SCF
   .CLOSED SHELL
    10
   **MOLECULE
   *BASIS
   .DEFAULT
    dyall.v2z
   **END OF

The calculation can be run with::

   $ ./pam --inp=4DC.inp --mol=ch4.xyz