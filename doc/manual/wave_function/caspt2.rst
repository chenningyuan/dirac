:orphan:
 

star(CASPT2)

Relativistic IVO/CASCI/CASPT2 and IVO/RASCI/RASPT2 module, Ref. :cite:`Masuda2025`

| Initial implementation by Minori Abe,
| parallelization, adaptation to DIRAC, and implement RASCI/RASPT2 functionality by Kohei Noda,
| implement :ref:`**MOLTRA` interface and IVO functionality by Sumika Iwamuro,
| and optimization by Kohei Noda and Yasuto Masuda.

**required keywords**
=====================

keyword(NINACT)

Specify the nubmer of inactive spinors.

keyword(NACT)

Specify the number of active (CAS or RAS) spinors.

keyword(NELEC)

Specify the number of electrons in the active space.

keyword(NSEC)

Specify the number of secondary spinors.

keyword(SUBPROGRAMS)

| Specify the subprograms to be used in the CASPT2 calculation.
| Currently the following subprograms are available
| - IVO
| - CASCI (perform RASCI or CASCI calculation based on your input)
| - CASPT2 (perform RASPT2 or CASPT2 calculation based on your input)
|
| for example, if you want to perform IVO/CASCI/CASPT2 calculation, you need to specify the subprograms as follows:

::

    .SUBPROGRAMS
    IVO
    CASCI
    CASPT2

**required keywords for CASCI/CASPT2 and RASCI/RASPT2**
=======================================================

keyword(CASPT2_CIROOTS)

| Specify the total symmetry number and the number of roots to be calculated in the CASCI/CASPT2 calculation.
| for example, if you want to calculate Ag symmetry in C2h molecule (total symmetry = 5) with 1,2,4,5 and 6 lowest root,
| you need to specify the keyword as follows:

::

    .CASPT2_CIROOTS
    5 1 2 4..6

| first number is the total symmetry number, and the rest are the roots to be calculated.
| you can also additional lines to specify the roots for other symmetries.
| for example:

::

    .CASPT2_CIROOTS
    5 1 2 4..6
    6 1 2

**required keywords for IVO**
=============================

keyword(NOCCG)

| Specify the number of gerade occupied orbitals (MOs).
| This keyword is required when the calculation is performed as a molecule with inversion center symmetry.
| Conversely, this keyword must not be used when the molecule does not have inversion center symmetry.


keyword(NOCCU)

| Specify the number of ungerade occupied orbitals (MOs).
| This keyword is required when the calculation is performed as a molecule with inversion center symmetry.
| Conversely, this keyword must not be used when the molecule does not have inversion center symmetry.


keyword(NOCC)

| Specify the number of occupied molecular orbitals (MOs).
| This keyword is required when the calculation is performed as a molecule without inversion center symmetry.
| Conversely, this keyword must not be used when the molecule has inversion center symmetry.

**required keyword for RASCI/RASPT2**
======================================

keyword(RAS2)

| Specify RAS2 spinor numbers.
| if you want to set the RAS2 orbitals 7, 8, 11, 12, 13, and 14 lowest spinors,
| you need to specify the keyword as follows:

::

    .RAS2
    7,8,11..14

| With this keyword, you can also include spinols that are not energetically continuous in the CAS of the CASCI/CASPT2, RASCI/RASPT2 calculation.
| for example, if you want to include 7, 8, 11, 12, 13, and 14 lowest spinors,
| you can specify the keyword as follows:

::

    .NINACT
    6
    .NACT
    6
    .RAS2
    7,8,11..14
    .SECONDARY
    6

| In this case, the first 6 spinors are inactive,
| 7, 8, 11, 12, 13, and 14 are active,
| and 9, 10, 15 and later are secondary spinors.

**optional keywords for RASCI/RASPT2**
======================================

keyword(RAS1)

| Specify RAS1 spinor numbers (row 1), and the maximum number of holes in RAS1 (row 2).
| If you want to set the RAS1 orbitals 3, 4, 5 and 6 lowest spinors, and the maximum number of holes in RAS1 is 2,
| you need to specify the keyword as follows:

::

    .RAS1
    3..6
    2

keyword(RAS3)

| Specify RAS3 spinor numbers (row 1), and the maximum number of electrons in RAS3 (row 2).
| If you want to set the RAS3 orbitals 15, 16, 17 and 18 lowest spinors, and the maximum number of electrons in RAS3 is 2,
| you need to specify the keyword as follows:

::

    .RAS3
    15..18
    2

**optional keywords for IVO**
=============================

keyword(NHOMO)

Specify the number of HOMO-like spinors

*Default:*

::

    .NHOMO
    0

keyword(NVCUTG)

| Specify the number of virtual cut (not used in IVO) gerade orbitals (MOs).
| This keyword can be used when the calculation is performed as a molecule with inversion center symmetry.
| Conversely, this keyword must not be used when the molecule does not have inversion center symmetry.


*Default:*

::

    .NVCUTG
    0

keyword(NVCUTU)

| Specify the number of virtual cut (not used in IVO) ungerade orbitals (MOs).
| This keyword can be used when the calculation is performed as a molecule with inversion center symmetry.
| Conversely, this keyword must not be used when the molecule does not have inversion center symmetry.


*Default:*

::

    .NVCUTU
    0

keyword(NVCUT)

| Specify the number of virtual cut (not used in IVO) orbitals (MOs).
| This keyword can be used when the calculation is performed as a molecule without inversion center symmetry.
| Conversely, this keyword must not be used when the molecule havs inversion center symmetry.


*Default:*

::

    .NVCUT
    0

**general optional keywords for this module**
=============================================

keyword(ESHIFT)

Specify the energy shift for the CASPT2 calculation.

*Default:*
::

    .ESHIFT
    0.0

keyword(MINHOLERAS1)
    
| Specify the minimum number of holes in RAS1.
| This keyword is only used in RASCI/RASPT2 calculation.
| NOTE: It is useful when you want to calculate excitation state energy calculation (including core excitation) with RASCI/RASPT2.

*Default:*
::

    .MINHOLERAS1
    0


keyword(SCHEME)

| Specify the MOLTRA scheme to be used in the CASPT2 calculation.
| Usually, you can use the same scheme as :ref:`MOLTRA_.SCHEME` under :ref:`**MOLTRA`.

*Default:*

::

    .SCHEME
    4

keyword(COUNTNDET)

Count and print the number of determinants for each total symmetry with current CASPT2 input.
and skip the any other calculations in this module.
(default: false, don't count and print the number of determinants).

keyword(DEBUGPRINT)

Print debugging information.
(default: false, don't print debugging information).

keyword(RESTART)

| Restart the CASPT2 calculation from the previous checkpoint file. (default: false, don't restart).
| NOTE: currently, restart functionality is only available for CASPT2/RASPT2 subprograms.
| In addition, `gen_caspt2_restart <https://raw.githubusercontent.com/RQC-HU/dirac_caspt2/refs/heads/main/tools/gen_dcaspt2_restart>`_ script is required to generate the restart file for the CASPT2/RASPT2 calculation.
| The following commands are required to generate the restart file and restart the CASPT2/RASPT2 calculation:

::

    $ gen_caspt2_restart <prev_output_file>
    $ pam --inp <input_file> --mol <molfile> --put "caspt2_restart_*"