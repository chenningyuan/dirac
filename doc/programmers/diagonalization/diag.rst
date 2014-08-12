=============================
Problems with diagonalization 
=============================

Diagonalizations are important part of (relativistic) quantum chemistry methods.

There is the diagonalization of the Lowdin matrix, then of the Jz MO matrix and, finally,
of the Fock MO matrix.

In DIRAC, there are few routines to accomplish these tasks:

 lapack' routine DSYEVR (real symmetric matrixes only)
 
 eispack's routine RS  (real symmetric matrixes only)

 own DIRAC routine JACOBI (real symmetric matrices only)

 own DIRAC routine QJACOBI (quaternion Hermitian matrices)

 own DIRAC routine based on Householder diagonalization (quaternion Hermitian matrices)

Sometimes one has to carefully choose the diagonalization method, depending either on the quality (accuracy)
of obtained results - eigenvalues and eigenvectors - or on the performance.


Problem case - DSYEVR
=====================

The lapack' routine DSYEVR, which is fast enough, is giving - for certain input matrices - 
numerically unstable eigenvectors.

This happens in the LINSYM routine, when the linear symmetry is on.

O2 test system
--------------
The test system, where discrepancies with DSYEVR are observed, is the O2 molecule of the fscc_highspin" test: ::
  
 pam --noarch --inp=hf.inp --mol=O2.mol

(Input files for download are from the DIRAC test directory, namely :download:`hf.inp<../../../../test/fscc_highspin/hf.inp>` and
:download:`O2.mol <../../../../test/fscc_highspin/O2.mol>` )

The proper Total SCF energy of this system is, ::

  -149.68666145184687

while the DSYEVR routine employed in the LINSYM part gives different value by about 10^-7,like ::

  -149.68666198782267

together with positron states intruding.

Molecular symmetries lower than linear symmetry (D2h, C2v ...) give proper SCF energies.


The smallest test system
------------------------
By "descending" to the smallest system where error is showing up we came to the artificial X2.1s1p1d system: ::

 pam --noarch --inp=dc.1el.2fs.inp  --mol=X2.1s1p1d.asd.mol --get "Jz_SS_matrix.fermirp=2-2"
 pam --noarch --inp=dc.1el.2fs.inp  --mol=X2.1s1p1d.D2h.mol

(Input files for download are :download:`dc.1el.2fs.inp <../../../../test/fscc_highspin/dc.1el.2fs.inp>`,
:download:`X2.1s1p1d.asd.mol <../../../../test/fscc_highspin/X2.1s1p1d.asd.mol>` and
:download:`X2.1s1p1d.D2h.mol <../../../../test/fscc_highspin/X2.1s1p1d.D2h.mol>`.)

In the linear symmetry eigenvectors obtained with the DSYEVR routine happen to be numerically inaccurate. 
This can be verified with the saved *Jz_SS_matrix.fermirp=2-2* matrix.

The demonstration matrix
~~~~~~~~~~~~~~~~~~~~~~~~
After tracking the source of the problem we converged to the 9x9 matrix, *Jz_SS_matrix.fermirp=2-2*, 
for the diagonalization testing, which is saved from the LINSYM routine.

With the DIRAC standalone diagonalization testing program, *diag.x*, we found discrepancies 
in numerical accuracy of eigenvectors
obtained with the DSYEVR routine in comparison to other diagonalization methods.

