:orphan:

.. _relccsd_developer:

RELCCSD: theory references and code map
========================================

This page links the equations of the papers describing the RELCCSD module
(``src/relccsd``) to the routines that implement them. It is meant for developers;
the input is described in the :ref:`manual <relccsd>`.

Throughout, the papers are referred to by the short keys below. The same keys are
used in comments in the source code, for instance ``VLD96 Eq. (5), term 4``.

====== =========================== ===================================================
Key    Reference                   Covers
====== =========================== ===================================================
VLD96  :cite:`Visscher1996`        CCSD and CCSD(T) energy (closed and open shell),
                                   double group symmetry, integral classes and memory
VEK01  :cite:`Visscher2001`        Fock-space CC: sectors, effective Hamiltonian,
                                   folded diagrams
SVS16  :cite:`Shee2016`            Λ equations, one-particle density matrix, tensor
                                   contraction routines (``contraction_xyz``)
SSVG18 :cite:`Shee2018`            EOM-CC (IP, EA, EE): σ-vectors, Davidson solver,
                                   odd-rank and non-totally symmetric contractions
PV03   :cite:`Pernpointner2003`    MPI parallelization of the CCSD code: symmetry-driven
                                   distribution of the integrals and of the work
====== =========================== ===================================================

The CCSD equations of VLD96, Eqs. (3)--(16), are identical to Eqs. (3)--(16) of
the earlier paper :cite:`Visscher1995`, whose remaining part describes a
Kramers-restricted closed-shell formalism that is *not* used in RELCCSD.
The intermediate Hamiltonian variants of Fock-space CC are described in
:cite:`Landau2004` and :cite:`Eliav2005`.


Ground-state CCSD (VLD96)
-------------------------

The iterations are driven by ``cceqn_driver_amplitudes``.
In each iteration the residuals S1 and S2 are built, and the new amplitudes are obtained
by dividing them by the orbital energy denominators. For this reason the orbital energies
are subtracted from the diagonals of ``HO`` and ``HV`` (``HOINTM``, ``HVINTM``).

===================================== ========================= ===============================
Quantity                              VLD96                     Routine / array
===================================== ========================= ===============================
CCSD energy                           Eq. (3)                   ``CCENG``
:math:`\tau` amplitudes               Eq. (8)                   ``GETTAU``
:math:`H_C^A`                         Eq. (9)                   ``HVINTM`` → ``HV``
:math:`H_I^K`                         Eq. (10)                  ``HOINTM`` → ``HO``
:math:`H_C^K`                         Eq. (11)                  ``HOVINTM`` → ``HOV``
:math:`G_C^A`                         Eq. (12)                  ``GVINTM`` → ``GV``
:math:`G_I^K`                         Eq. (13)                  ``GOINTM`` → ``GO``
:math:`A_{IJ}^{KL}`                   Eq. (14)                  ``AINTM`` → ``AZ``
:math:`B_{CD}^{AB}`                   Eq. (15)                  ``BINTM`` (contracted directly)
:math:`H_{IC}^{AK}`                   Eq. (16)                  ``HINTM`` → ``H``
T1 equations                          Eq. (4)                   ``amplitude_equation_T1``
T2 equations                          Eq. (5)                   ``amplitude_equation_T2``
first-order terms                     Eqs. (4)/(5), term 1      ``S1S2_ADDFO``
(T), +T and -T corrections            Eqs. (17)--(22)           ``cceqn_t3_energy_correction.F``
===================================== ========================= ===============================

Term numbering
~~~~~~~~~~~~~~

The papers do not number the individual terms. In RELCCSD they are counted in order of
appearance, splitting bracketed expressions into their parts. Implicit summation over
repeated indices is assumed.

T1 equations, VLD96 Eq. (4):

==== ==================================================== =====================================
Term Expression                                           Where
==== ==================================================== =====================================
1    :math:`F_I^A`                                        ``S1S2_ADDFO``
2    :math:`-2F_C^K T_K^A T_I^C`                          ``amplitude_equation_T1``, with 5b
3    :math:`+H_C^A T_I^C`                                 ``amplitude_equation_T1``
4    :math:`-H_I^K T_K^A`                                 ``amplitude_equation_T1``
5a   :math:`+H_C^K T_{IK}^{AC}`                           ``amplitude_equation_T1``
5b   :math:`+H_C^K T_K^A T_I^C`                           ``amplitude_equation_T1``, with 2
6    :math:`+V_{IC}^{AK} T_K^C`                           ``amplitude_equation_T1``
7    :math:`+\sum_{C<D} V_{CD}^{AK}\tau_{IK}^{CD}`        ``T1EQN2`` (VOVV, batched)
8    :math:`-\sum_{K<L} V_{IC}^{KL}\tau_{KL}^{AC}`        ``T1EQN1``
==== ==================================================== =====================================

T2 equations, VLD96 Eq. (5), with :math:`P_{PQ}f(P,Q) = f(P,Q) - f(Q,P)`, Eq. (6):

==== ==================================================================== ===============================
Term Expression                                                           Where
==== ==================================================================== ===============================
1    :math:`V_{IJ}^{AB}`                                                  ``S1S2_ADDFO``
2    :math:`+P_{AB}\,G_C^A T_{IJ}^{CB}`                                   ``amplitude_equation_T2``
3    :math:`-P_{AB}\,V_{IJ}^{AK} T_K^B`                                   ``amplitude_equation_T2``, with 9
4    :math:`-P_{IJ}\,G_I^K T_{KJ}^{AB}`                                   ``amplitude_equation_T2``
5    :math:`+P_{IJ}\,V_{IC}^{AB} T_J^C`                                   ``amplitude_equation_T2`` (VOVV, batched)
6    :math:`+\sum_{K<L} A_{IJ}^{KL}\tau_{KL}^{AB}`                        ``amplitude_equation_T2``
7    :math:`+\sum_{C<D} B_{CD}^{AB}\tau_{IJ}^{CD}`                        ``BINTM`` (VVVV, batched)
8    :math:`+P_{IJ}P_{AB}\,H_{IC}^{AK} T_{JK}^{BC}`                       ``amplitude_equation_T2``
9    :math:`-P_{IJ}P_{AB}\,V_{IC}^{AK} T_J^C T_K^B`                       ``amplitude_equation_T2``, with 3
==== ==================================================================== ===============================

Several terms are coded in a different but equivalent form, for instance
:math:`-P_{IJ}\,T_{IK}^{AB} G_J^K` for term 4. Under the permutation operator these are
identical to the expressions above.

Integral classes and memory
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two-electron integrals are sorted into six classes (OOOO, VOOO, VVOO, VOVO, VOVV, VVVV;
VLD96 Sec. IV). The VVVV and VOVV classes are in general too large to be held in memory;
the terms using them (T1 term 7, T2 terms 5 and 7) read them in fixed-size batches, so that the
memory requirement stays at a few times the size of the T2 vector (VLD96 Sec. IV).


MPI parallelization (PV03)
--------------------------

The parallelization of the CCSD code is described in PV03. It is symmetry driven: all
integral and amplitude arrays are blocked by the irrep of a compound (pair) index, and the
work is divided by distributing, within each irrep, the columns of these blocks over the
MPI processes. This needs no communication during the contractions, and the large integral
classes are divided over the processes, so that both the memory and the disk space per
process decrease with the number of processes.

Distribution lists
   ``DISLIST`` (``ccdriv.F``) sets up, for every irrep, the range of columns handled by each
   process and stores it in ``IDIST`` (common block in ``ccpar.inc``, together with
   ``MASTER``, ``MYPROC`` and ``NMPROC``). The columns are divided evenly; the highest
   processes get one column more when the number does not divide. ``IDIST(:,1,:)`` and
   ``IDIST(:,2,:)`` describe the VVVV and VOVV integrals, which are distributed over their
   second (VV pair) index; they hold the local column range and the offsets that translate
   the absolute integral index into the index in the local array. Other entries give the
   column ranges used to split contractions that involve only replicated arrays:
   ``IDIST(:,3,1)`` for T1 terms 5a and 6, ``IDIST(:,4,:)`` for the VOVO contraction in
   ``HINTM``, ``IDIST(:,5,:)`` for T2 term 6 (:math:`A\cdot\tau`) and ``IDIST(:,10,:)`` for
   T2 term 8 (:math:`H\cdot T_2`). Entries 6--9 are set up but not used at present; the
   corresponding terms are computed by the master. ``NV5PAR`` and ``NV6PAR`` are the numbers
   of VOVV and VVVV integrals held by the process.

Integral sorting
   The transformed integrals (``MDCINT``, or the ``MDCINX*`` files of a parallel 4-index
   transformation, which are read by each process in turn) are broadcast by the reading
   process (``GET2INTS``, ``ccints.F``). Every process keeps all OOOO, VOOO, VVOO and VOVO
   integrals but only its own slices of the VOVV and VVVV classes, and writes these to its
   own sorted integral files. ``GETVOVV`` and ``GETVVVV`` (``ccgetv.F``) return only the
   local integrals.

CCSD iterations
   The amplitudes T1 and T2 and the small intermediates are replicated on all processes.
   Each process computes the contributions of its local VOVV and VVVV integrals (T1 term 7,
   T2 terms 5 and 7) and its share of the distributed contractions. Terms that are cheap or
   cannot be split are computed by the master alone (blocks marked *TASK SPLITTING*). The
   partial residuals S1 and S2 are summed on the master and broadcast to all processes in
   ``PARTS`` (``cceqn_utils.F``).

   Because ``PARTS`` adds up the residuals of all processes, each contribution must be
   computed exactly once over all processes. Two examples show what can go wrong:

   * T2 term 6 (:math:`A\cdot\tau`) is computed on every process, but each process
     contracts only its own share of the summation index: ``IDIST(1,5,IREP)`` columns,
     starting at the offsets ``IDIST(2,5,IREP)`` in :math:`\tau` and ``IDIST(3,5,IREP)`` in
     :math:`A`. If a process looped over all columns instead, every process would add the
     complete term, and after ``PARTS`` the term would appear ``NMPROC`` times in S2.
   * T1 term 7 (``T1EQN2``) uses the VOVV integrals, which ``GETVOVV`` returns only for the
     local slice. It is therefore computed on every process, and only the sum over the
     processes is complete. If this call were moved into a *TASK SPLITTING* block, only the
     master's slice would be used and the term would be incomplete in parallel runs, while
     serial runs would still be correct.

   So a term that is computed on every process must be restricted to the local integrals or
   the local column range, and a term in a *TASK SPLITTING* block may only use arrays that
   are complete on the master. Errors of either kind do not show up in serial runs, so new
   code should always be tested with more than one process.

(T) correction
   The parallelization of the triples correction is not part of PV03. It uses the same
   distribution of the VOVV integrals: each process computes its slice of the VVVO
   contribution for every triple, the triples are assigned round-robin to the processes,
   and the energy contributions are summed at the end (see the header of
   ``cceqn_t3_energy_correction.F``).

The Λ equations, density matrices and EOM-CC use a different strategy: the VOVV integrals
are gathered on every process (see the next sections).


Fock-space CC (VEK01)
---------------------

VEK01 does not restate the CCSD equations; it refers to them as Eq. (VLD-n) and lists the
changes for the higher sectors:

* The T1 and T2 arrays are enlarged to hold the amplitudes and the effective Hamiltonian
  elements of all sectors up to (2h,2p), so that the same kernels
  (``amplitude_equation_T1/T2``) are used in every sector.
* ``DELFCK`` and ``DELINT`` (``cceqn_utils.F``) remove the matrix elements that do not belong
  to the sector being solved.
* The effective Hamiltonian elements are the corresponding elements of the residual,
  :math:`(H_\mathrm{eff})_{EF} = S_{EF}`.
* The folded diagrams (item 2 of VEK01 Sec. II B) are added in ``FOLDED`` (``ccfspc.F``), giving
  the amplitude update of VEK01 Eq. (6). ``DENOMF`` applies the energy denominators.
* After convergence of a sector, :math:`H_\mathrm{eff}` is diagonalized (``ccfspc.F``).


Λ equations and density matrix (SVS16)
--------------------------------------

The steps of SVS16 Sec. III B correspond to:

======================================== ======================= ================================================
Quantity                                 SVS16                   Routine
======================================== ======================= ================================================
:math:`\bar F_m^i`                       Eq. (A4)                ``f_bar_mi`` (``intermediates.F90``)
:math:`\bar F_a^e`                       Eq. (A5)                ``f_bar_ae``
:math:`\bar F_e^m`                       Eq. (A6)                ``f_bar_me``
:math:`W_{mn}^{ij}`                      Eq. (A7)                ``w_ijmn``
:math:`W_{ej}^{mb}`                      Eq. (A8)                ``w_ejmb``
:math:`W_{mn}^{ie}`                      Eq. (A9)                ``w_iemn``
:math:`W_{am}^{ef}`                      Eq. (A10)               ``w_efam``
:math:`W_{ab}^{ef}`                      Eq. (A11)               ``w_efab_out_of_core`` (VVVV, batched)
:math:`\bar W_{ej}^{mb}`                 Eq. (A12)               ``wbar_mbej``
:math:`W_{ie}^{mn}`                      Eq. (A13)               ``w_mnie``
:math:`W_{ef}^{am}`                      Eq. (A14)               ``w_amef``
:math:`G_a^e`, :math:`G_m^i`             Eqs. (A15), (A16)       ``G_intm``
:math:`\Lambda_1` equation               Eq. (A1)                ``lambda_equation_L1`` (``lambda_equation.F90``)
:math:`\Lambda_2` equation               Eq. (A2)                ``lambda_equation_L2``
one-particle density matrix              Eqs. (A17)--(A20)       ``cc_density`` (``ccgrad.F``)
======================================== ======================= ================================================

The iterations are driven by ``cceqn_driver_lambda``.
Routines with the suffix ``_MBPT2`` are the corresponding MP2-level variants.

As stated in SVS16 Sec. III B, the VOVV integrals, distributed as in PV03, are *gathered on every
process* for these steps (``GETVOVV_incore``), while the VVVV integrals stay distributed and
only local contributions of :math:`W_{ab}^{ef}` are used. The in-core routines therefore must
not be used in terms of the energy code that are computed on every process and summed in
``PARTS``.


EOM-CC (SSVG18)
---------------

``eom_driver.F90`` sets up the calculation, ``IM_container`` (``sigma_eom.F90``) builds the
intermediates of SSVG18 Eqs. (B13)--(B25), which are those of SVS16, and the Davidson solver for
non-Hermitian matrices (Hirao--Nakatsuji, SSVG18 Sec. III B) is in ``davidson/davidson.F90``.
Left and right eigenvectors are obtained in separate runs.

============ ================== ==================== =========================
EOM type     right σ-vector     left σ-vector        SSVG18
============ ================== ==================== =========================
IP           ``sigma_IP``       ``Lambda_IP``        Eqs. (B1)--(B2) / (B3)--(B4)
EA           ``sigma_EA``       *not implemented*    Eqs. (B5)--(B6) / (B7)--(B8)
EE           ``sigma_ee_1h1p``  ``Lambda_EE``        Eqs. (B9)--(B10) / (B11)--(B12)
============ ================== ==================== =========================

The left EOM-EA σ-vector is not computed (the call of ``Lambda_EA`` in
``create_sigma_vector_left`` is commented out), so left EOM-EA eigenvectors, and hence
EOM-EA properties (``.EOMPROP``), are not available.

For IP and EA the odd-rank tensors are made even-rank by adding a "continuum" orbital
(index ``c1``), which carries the symmetry of the removed or added electron
(SSVG18 Sec. III A, ``set_continuous_orbital``). Trial vectors are generated in
``trial_vectors_eom.F90``.


Tensor contraction routines (SVS16 App. B, SSVG18 App. A)
---------------------------------------------------------

The routines ``contraction_xyz`` in ``contraction.F90`` contract a rank-x tensor with a
rank-y tensor into a rank-z tensor (xyz = 222, 242, 422, 442, 444, 244, 424). They map the
contraction onto a matrix multiplication by sorting the input tensors, call BLAS, sort the
result back, and, where needed, antisymmetrize it and pack it into triangular storage.
Example from SVS16 Eq. (B1):

.. code-block:: fortran

   call contraction_444((/"o1","o3","p1","p3"/), (/"o2","p3","p2","o3"/), &
                        (/"o1","o2","p1","p2"/), ResL2, 1.0d0, 1.0d0, nrep, &
                        RightTensor=w_voov, LeftTensor=L2)

Arguments, in order: the index strings of the left, right and output tensors, in the order
in which they appear in the diagram (``o`` occupied, ``p`` virtual, ``c`` continuum orbital,
numbered); the output tensor; the diagram factor; the scaling factor of the product; ``nrep``;
optionally the input tensors (if absent, the integrals are fetched from file); and optionally
the irreps of the left and right input tensors when these are not totally symmetric
(SSVG18 App. A).

Points to keep in mind when using these routines:

* The output tensor is an assumed-shape array, and the routines use ``size(ImTensor)``
  when they add the sorted product to it. The actual argument must therefore be declared
  with exactly the size of the tensor, not as a larger work buffer.
* The source comments call the second scaling factor "beta in DGEMM". This is only true
  when no sorting of the product is needed. Otherwise the product is added to the output
  tensor multiplied by this factor, and the old content of the output tensor is not scaled.
  Use 1.0 and zero the output tensor beforehand if needed.
