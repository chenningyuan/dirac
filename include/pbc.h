C
C Variables required for periodic boundary conditions.
C
C
C We need a maximum multipole order and a maximum number of
C multipole moments. Note that qsiz = (qlmax+1)^2.
C (N.B. max_qlmax_pbc should ideally be equal to MXQNM in the
C       common block include/maxmom.h)
C
      integer      max_qlmax_pbc, max_qsiz_pbc
      parameter    ( max_qlmax_pbc = 21 )
      parameter    ( max_qsiz_pbc = (1+max_qlmax_pbc)**2 )
C
C We need a place to store the multipole moments of the central
C cell.
C The variable mm_q_loss_thres specifies how much charge may be left
C unaccounted for by the multipole moments before it is considered
C a serious error.
C
      double precision cell_mulmom(max_qsiz_pbc)
      double precision ovl_mulmom(max_qsiz_pbc)
      double precision mm_q_loss_thres
C
C
C The variable max_ovlfcomp stores the maximum element of an
C overlap matrix. This helps determine the the necessary size
C of the near-field. max_kinfcomp and max_hamfcomp are the corresponding
C variables for kinetic energy and (part of the) one-electron hamiltonian.
C
      double precision max_ovlfcomp, max_kinfcomp, max_hamfcomp
C
C The variables arg_latindex (=i,j,k) and
C arg_latcoord = ia+jb+kc, with a,b,c being lattice vectors,
C are passed as arguments to the integral computations. These
C variables specify an offset for the B-center of one-electron
C integrals.
C
      double precision rspc_vec_std(3), nf_vec_std(3),
     &                 nf2_vec_std(3)
      integer          rspc_vec_index(3), nf_vec_index(3),
     &                 nf2_vec_index(3)
C
C threshold for basis function diffuseness
C
      double precision  pbc_diffthres
C
C For the far-field multipole expansion we need a parameter specifying
C the actual order of the expansion
C
      integer      qlmax_pbc
C
C We need an lmax for W and T tensors.
C
      integer      twlmax_pbc
C
C The T tensor representing the cell-far field interaction should
C be iterated a specific number of times. We store that number in
C
      integer      num_tlat_its
C
C we need a file handle when writing files specific to the PBC
C calculations
C
      integer         LUPBC
C
C the variable rspc_cutoff_1op is set automatically to the number
C of real-space layers of 1-electron operators needed to come
C within a requested accuracy;
C the variable rspc_cutoff_dmat is the maximum number of real-space
C layers of density matrix components that will be computed and is
C set when the input files are read.
C
      integer         rspc_cutoff_mm_dmat
C
C dopbc is set when the .dal file is read and indicates whether
C periodic boundary conditions should be used.
C
      logical         dopbc, pbcactive
C
C the flag atom_coord_lat indicates if the user gives the atom
C positions in lattice coordinates (true) or standard coordinates
C (false).
C
      logical         atom_coord_lat
C
C flag indicating whether diffuse functions should be removed from basis
C
      logical        pbc_no_diffuse
C
C flag indicating Hartree calculation (i.e. no DFT)
      logical        pbcdohartree
C
C flag indicating whether each R tensor in the near-field Coulomb calc.
C should be handled separately. (Doing this is slow; it is better to
C add up R tensors for different near-field translations and only contract
C with the sum.)
C
      logical        pbc_sep_Rtensor
C
      COMMON /pbcl/ dopbc, pbcactive, atom_coord_lat,
     &              pbc_no_diffuse, pbcdohartree,
     &              pbc_sep_Rtensor
      COMMON /pbci/
     &             LUPBC, rspc_cutoff_mm_dmat,
     &             rspc_vec_index, nf_vec_index,
     &             nf2_vec_index,
     &             qlmax_pbc, twlmax_pbc,
     &             num_tlat_its
      COMMON /pbcd/
     &             rspc_vec_std, nf_vec_std,
     &             nf2_vec_std,
     &             max_ovlfcomp, max_kinfcomp, max_hamfcomp,
     &             cell_mulmom, ovl_mulmom,
     &             mm_q_loss_thres, pbc_diffthres
C
C


C
C     this is kept to exemplify syntax:
C     defaults:
C
C      PARAMETER ( MMSKIPDF = .FALSE., MMCUTDF = .TRUE. )
C      PARAMETER ( MMPRDF = 0, LMAXDF = 8 )


