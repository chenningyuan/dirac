=================================
Problems with the diagonalization 
=================================

Diagonalizations are important part of (relativistic) quantum chemistry methods.

There is the diagonalization of the Lowdin matrix, then of the Jz MO matrix and, finally,
of the Fock MO matrix.

In DIRAC, there are few routines to accomplish these tasks:

* lapack' routine DSYEVR (for real symmetric matrixes only)
 
* eispack's routine RS  (for real symmetric matrixes only)

* own DIRAC routine JACOBI (for real symmetric matrices only)

* own DIRAC routine QJACOBI (for quaternion Hermitian matrices)

* and own DIRAC routine based on Householder diagonalization (for quaternion Hermitian matrices)

Sometimes one has to carefully choose the diagonalization method, depending either on the quality (accuracy)
of obtained results - eigenvalues and eigenvectors - or on the time performance.


Problem case - DSYEVR
=====================

The lapack' routine for diagonalization, DSYEVR, which is fast enough, is giving - for certain input matrices - 
numerically unstable eigenvectors.

This happens in the LINSYM routine, when the linear symmetry is turned on.

O2 test system
--------------
The real test system, where numerical discrepancies with DSYEVR are observed, is the O2 molecule of the *fscc_highspin* test: ::
  
 pam --noarch --inp=hf.inp --mol=O2.mol

(Input files for download are from the DIRAC test directory, namely :download:`hf.inp <../../../test/fscc_highspin/hf.inp>` and
:download:`O2.mol <../../../test/fscc_highspin/O2.mol>` )

The proper Total SCF energy of this system is, ::

  -149.68666145184687

while the DSYEVR routine employed in the LINSYM procedure gives different value by about 10^-7, like (obtained by Miro) ::

  -149.68666198782267

together with some positron states intruding.

Molecular symmetries lower than linear symmetry (D2h, C2v ...) give proper and identical SCF energies.


The smallest test system
------------------------
By "descending" to the smallest system where error is clearly showing up we came to the artificial X2.1s1p1d system: ::

 pam --noarch --inp=dc.1el.2fs.inp  --mol=X2.1s1p1d.asd.mol --get "Jz_SS_matrix.fermirp=2-2"
 pam --noarch --inp=dc.1el.2fs.inp  --mol=X2.1s1p1d.D2h.mol

(Input files for download are :download:`dc.1el.2fs.inp <../../../test/fscc_highspin/dc.1el.2fs.inp>`,
:download:`X2.1s1p1d.asd.mol <../../../test/fscc_highspin/X2.1s1p1d.asd.mol>` and
:download:`X2.1s1p1d.D2h.mol <../../../test/fscc_highspin/X2.1s1p1d.D2h.mol>`.)

The affected value is the "Inactive energy (Output from ONEERG)" in the output: ::

  -15622.826231066265 (wrong value,    linear symmetry)
  -15622.826368346343 (correct value,  linear symmetry)
  -15622.826368346305 (correct value,  D2h symmetry)

The demonstration program 
-------------------------
After tracking the source of the problem we get to the 9x9 matrix of the X2.1s1p1d test system, *Jz_SS_matrix.fermirp=2-2*, 
for the diagonalization testing, which is saved in the LINSYM routine.

In the linear symmetry eigenvectors obtained with the DSYEVR routine happen to be numerically inaccurate. 
This can be verified with the abovementioned matrix.

With the DIRAC standalone diagonalization testing program, *diag.x* (its source code is in the utils directory, 
together with the text input file :download:`DIAGONALIZE_TESTS.INP <../../../utils/DIAGONALIZE_TESTS.INP>`), we found discrepancies 
in numerical accuracy of eigenvectors (U)
obtained with the DSYEVR routine in comparison to other diagonalization methods. Eigenvalues (eps) are unaffected. 
A is the original matrix for diagonalization, I is the unit matrix. The norm is defined as sum of absolute values of elements divided by total number of elements.
We split matrix elements into diagonal and off-diagonal.
Values around 10^(-6) clearly indicate deviation from zero threshold, where one would expect values around 10^(-14).

Intel+MKL-i8:
~~~~~~~~~~~~~
::

 U^{+}*A*U - eps =  0 = C: 
 C_norm_diag:0.1684D-14    C_norm_offdiag:0.6623D-06
 U^{+}*U -I =  0 = AXr: 
 AXr_norm_diag:0.2097D-15  AXr_norm_offdiag:0.1325D-05
 U*U^{+} -I =  0 = AY: 
 AYr_norm_diag:0.6293D-05  AYr_norm_offdiag:0.4174D-05

GNU+ownmath-i8:
~~~~~~~~~~~~~~~
::

 U^{+}*A*U - eps =  0 = C: 
 C_norm_diag:0.1739D-14  C_norm_offdiag:0.5091D-06  
 U^{+}*U -I =  0 = AXr: 
 AXr_norm_diag:0.2097D-15  AXr_norm_offdiag:0.1018D-05
 U*U^{+} -I =  0 = AY: 
 AYr_norm_diag:0.5239D-05  AYr_norm_offdiag:0.3407D-05

The solution is to keep non-DSYEVR diagonalization routine in the LINSYM part to avoid possibility of highly degenerate values with numerically unstable eigenvectors.

