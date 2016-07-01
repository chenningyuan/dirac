:orphan:

UF6 molecule
============

The UF6 (146 electrons, Oh) is known compound connected with the uranium enrichment.
We shall investigate UF6 as the closed-shell (singlet) state 
in the highest D2h* computational symmetry 
(UF6 mol file :download:`UF6.v2z.D2h.mol <../../../../../test/tutorial_UF6/UF6.v2z.D2h.mol>`).

How do we obtain converged molecular orbitals with the all-electron method ?
The default SCF automatic occupation scheme does not work.
Iteartions are interrupted with message about too many changes in the SCF occupation, meaning
that no convergence was achieved.

Clearly, another route must be taken.

Atomic start and the autooccupation 
-----------------------------------

We demostrate the atomic start combined with the automatic occupation scheme,
which successfully give converged state with the occupation of (E1g,E1u)=(72,74),
see the input file, 
:download:`x2c_a4p.atstart_scf_autoocc_2fs.inp <../../../../../test/tutorial_UF6/x2c_a4p.atstart_scf_autoocc_2fs.inp>`, and 
the :download:`output file <../../../../../test/tutorial_UF6/result/x2c_a4p.atstart_scf_autoocc_2fs_UF6.v2z.D2h.out>`.

Input files for corresponding F,U atomic runs are
:download:`F.v2z.D2h.mol <../../../../../test/tutorial_UF6/F.v2z.D2h.mol>`,
:download:`F.x2c.scf_2fs.inp <../../../../../test/tutorial_UF6/F.x2c.scf_2fs.inp>`,
:download:`U.v2z.D2h.mol <../../../../../test/tutorial_UF6/U.v2z.D2h.mol>`,
:download:`x2c.scf_autoocc.inp <../../../../../test/tutorial_UF6/U.v2z.D2h.mol>`.
Output files from these (F,U) atomic runs are
:download:`F.x2c.scf_2fs_F.v2z.D2h.out <../../../../../test/tutorial_UF6/result/F.x2c.scf_2fs_F.v2z.D2h.out>` and
:download:`x2c.scf_autoocc_U.v2z.D2h.out <../../../../../test/tutorial_UF6/result/x2c.scf_autoocc_U.v2z.D2h.out>`.

Note that it is not necessary to specify proper atomic configuration of the U atom as the auto-occupation 
scheme does the job and produces suitable starting MOs.

Verification of the occupation
------------------------------

Next, we have to verify that this occupation of spinors gives the lowest energy.
For that, we perform several atomic-start based SCF calculations with different
occupation schemes. 
Results - occupations and SCF energy - are collected in the following table:

==========  ===================
occupation  SCF energy
==========  ===================
68   78     not conv
70   76     -28636.011248088082
72   74     -28636.624880692940 
74   72     -28635.715232162675
76   70     -28634.793937853457
==========  ===================

We see that manually assigned occupation of 72,74 gives the lowest energy, what is the
proof of correct autooccupation settings (see one input file, :download:`UF6.x2c_a4p.atstart_scf_autoocc_72_74 <../../../../../test/tutorial_UF6/UF6.x2c_a4p.atstart_scf_autoocc_72_74.inp>` )

Geometry optimization
---------------------

The UF6 molecule spans high Oh symmetry. Therefore, the geometry optimization
requires change of one parameter only - the U-F bond distance.

To perform the UF6 geometry optimization with the BP86 functional effectively, let us first obtain 
starting set of (DFT) molecular orbitals. For that, we first run single point DFT step with molecular
orbitals file from the previous SCF step and reuse them for the DFT-geometry optimization. 
See the BP86 input file, :download:`UF6.x2c_a4p.bp86_72_74.inp <../../../../../test/tutorial_UF6/UF6.x2c_a4p.bp86_72_74.inp>`, and the geometry optimization input, 
:download:`UF6.geomopt.x2c_a4p.bp86_72_74.inp <../../../../../test/tutorial_UF6/UF6.geomopt.x2c_a4p.bp86_72_74.inp>`.

Thanks to this MO start, number SCF iterations during the geometry optimization cycle is less than 10,
see the corresponding :download:`output file <../../../../../test/tutorial_UF6/result/UF6.geomopt.x2c_a4p.bp86_72_74_UF6.v2z.D2h.out>` . The final U-F bond distance obtained with the X2C-BP86-v2z method is 2.020A (experimental distance is 1.999A) 

Testing Python script
---------------------

The comprehensive testing script of this tutorial is :download:`test <../../../../../test/tutorial_UF6/test>`.
Pay attention to comments and commented files.
