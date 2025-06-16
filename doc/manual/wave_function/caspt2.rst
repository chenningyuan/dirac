:orphan:
 

star(CASPT2)

Relativistic IVO/CASCI/CASPT2 and IVO/RASCI/RASPT2 module, Ref. :cite:`Masuda2025`,Ref. :cite:`Abe2006`,Ref. :cite:`Anderson1990`,Ref. :cite:`Huzinaga1970`

| Initial implementation: Minori Abe
| Parallelisation, DIRAC integration, and RASCI/RASPT2 extension: Kohei Noda
| :ref:`**MOLTRA` interface and IVO functionality: Sumika Iwamuro
| Code optimisation: Kohei Noda and Yasuto Masuda
|
| The program is based on the CASPT2D formalism reported in :cite:`Anderson1990`.

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
| Currently the following subprograms are available:
| - IVO
| - CASCI (perform RASCI or CASCI calculation based on your input)
| - CASPT2 (perform RASPT2 or CASPT2 calculation based on your input)
|
| For example, if you want to carry out an IVO/CASCI/CASPT2 calculation sequentially, you must specify the sub-programs as follows:

::

    .SUBPROGRAMS
    IVO
    CASCI
    CASPT2

**required keywords for CASCI/CASPT2 and RASCI/RASPT2**
=======================================================

keyword(CASPT2_CIROOTS)

| Specify the total symmetry index and the roots to be included in the CASCI/CASPT2 calculation.
| For example, to compute the Ag symmetry of a C₂h molecule (symmetry index 5) and request the 1st, 2nd, 4th, 5th, and 6th lowest-energy roots, 
| you need to set the keyword as follows:

::

    .CASPT2_CIROOTS
    5 1 2 4..6

| The first number is the total-symmetry index; the remaining numbers specify the roots to be calculated.
| Additional lines can be added to define roots for other symmetries.
| For example:

::

    .CASPT2_CIROOTS
    5 1 2 4..6
    6 1 2

**required keywords for IVO**
=============================

keyword(NOCCG)

| Specify the number of occupied gerade Kramers‐pair orbitals (count pairs, not individual spinors).
| This keyword must be specified for molecules that possess an inversion center. 
| Conversely, it should not be used for molecules without inversion symmetry.


keyword(NOCCU)

| Specify the number of occupied ungerade Kramers‐pair orbitals (count pairs, not individual spinors).
| This keyword must be specified for molecules that possess an inversion center. 
| Conversely, it should not be used for molecules without inversion symmetry.


keyword(NOCC)

| Specify the number of occupied Kramers‐pair orbitals (count pairs, not individual spinors).
| This keyword must be specified for molecules that lack an inversion center; 
| it should not be used for molecules with inversion symmetry.

**required keyword for RASCI/RASPT2**
======================================

keyword(RAS2)

| Specify the RAS2 spinor indices.
| For example, to assign the spinors ranked 7th, 8th, 11th, 12th, 13th, and 14th lowest in energy to the RAS2 space,
| set the keyword as follows:
::

    .RAS2
    7,8,11..14

| By using this keyword, you can select spinors that are not contiguous in energy for inclusion in the CAS of a CASCI/CASPT2 or RASCI/RASPT2 calculation.
| For example, to include the spinors ranked 7th, 8th, 11th, 12th, 13th, and 14th lowest in energy, 
| set the keyword as follows:

::

    .NINACT
    6
    .NACT
    6
    .RAS2
    7,8,11..14
    .SECONDARY
    6

| In this case, spinors 1–6 are inactive; 
| spinors 7, 8, 11, 12, 13, and 14 are active;
| and spinors 9, 10, 15-18 are secondary.

**optional keywords for RASCI/RASPT2**
======================================

keyword(RAS1)

| Specify the RAS1 spinor indices on row 1 and the maximum number of holes allowed in RAS1 on row 2.
| For example, to assign spinors 3, 4, 5, and 6 (counted in ascending energy) to RAS1 and to allow up to two holes,
| set the keyword as follows:

::

    .RAS1
    3..6
    2

keyword(RAS3)

| Specify the RAS3 spinor indices on row 1 and the maximum number of electrons allowed in RAS3 on row 2.
| For example, to assign spinors 15, 16, 17, and 18 (counted in ascending energy) to RAS3 and to allow up to two electrons, 
| set the keyword as follows:

::

    .RAS3
    15..18
    2

**optional keywords for IVO**
=============================

keyword(NHOMO)

| Specify the number of HOMO-like spinors. 
| By default (NHOMO = 0), this value is determined automatically by counting all occupied spinors whose orbital energies lie within 0.1 a.u. of the HOMO energy.
| If desired, you can override the default by specifying NHOMO explicitly in the input file (e.g., 2, 4, or 6).

*Default:*

::

    .NHOMO
    0

keyword(NVCUTG)

| Specify the number of gerade virtual Kramers-pair orbitals to be excluded in the IVO scheme.
| This keyword is applicable only to molecules with an inversion center; do not use it for molecules without inversion symmetry.


*Default:*

::

    .NVCUTG
    0

keyword(NVCUTU)

| Specify the number of ungerade virtual Kramers-pair orbitals to be excluded in the IVO scheme.
| This keyword is applicable only to molecules with an inversion center; do not use it for molecules without inversion symmetry.

*Default:*

::

    .NVCUTU
    0

keyword(NVCUT)

| Specify the number of virtual Kramers-pair orbitals to be excluded in the IVO scheme.
| This keyword should be used only for molecules without an inversion center; do not use it for molecules with inversion symmetry.


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
    
| Specify the minimum number of holes allowed in RAS1.
| This keyword is used only in RASCI/RASPT2 calculations.
| Note: Setting this parameter may be helpful for excited-state calculations (including core excitations) performed with RASCI/RASPT2.

*Default:*
::

    .MINHOLERAS1
    0


keyword(SCHEME)

| Specify the MOLTRA scheme to be used for the CASPT2 calculation.
| In most cases, you can simply reuse the scheme defined under :ref:`**MOLTRA` (See :ref:`MOLTRA_.SCHEME`).

*Default:*

::

    .SCHEME
    4

keyword(COUNTNDET)

| Count and print the number of determinants within each irreducible representation as defined by the current CASPT2 input,
| skipping all other calculations in this module.
| (Default: false — no determinant counting is performed.)

keyword(DEBUGPRINT)

| Print debugging information.
| (Default: false — no debugging information is printed.)

keyword(RESTART)

| Restart the CASPT2 calculation from an existing checkpoint file.
| (Default: false — the calculation starts from scratch.)
| Note: Restart capability is currently implemented only for the CASPT2 or RASPT2 sub-programs.
| To create a valid restart file, run the helper script
|  `gen_caspt2_restart <https://raw.githubusercontent.com/RQC-HU/dirac_caspt2/refs/heads/main/tools/gen_dcaspt2_restart>`_ .
| The typical workflow is:

::

    $ gen_caspt2_restart <prev_output_file>
    $ pam --inp <input_file> --mol <molfile> --put "caspt2_restart_*"