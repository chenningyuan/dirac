module foldy_wouthuysen

! Small numerical helpers for Foldy-Wouthuysen/X2C coefficient transforms.
! This module deliberately avoids dependencies on ROSE data types and file I/O.

  use rose_global, only: luout
  use linear_algebra, only: hermitian_matrix_sqrt_pair, matrix_mul, &
                            transpose_conjg, diff_norm2, hermitian_matrix_power, diag_matrix

  implicit none

  private

  public :: construct_fw_coefficients

contains

   subroutine construct_fw_coefficients(CL,CS,SLL,SSS,CFW)

      ! Foldy-Wouthuysen transformed coefficients, formulated in the Gram matrices of the spinors:
      !
      !    GL = CL^+ SLL CL ,   GS = CS^+ SSS CS            (both nspinor x nspinor)
      !    CFW = CL GL^(-1/2) (1 + GL^(-1/2) GS GL^(-1/2))^(1/2) GL^(1/2)
      !
      ! This is identical in exact arithmetic to the construction via X = CS (CL^+ SLL CL)^(-1) CL^+ SLL and
      ! R = (SLL + X^+ SSS X)^(1/2), but SLL and SSS enter only through GL and GS. No square root or inverse of
      ! an AO overlap matrix is needed, so (nearly) linearly dependent AO sets, such as Cartesian basis
      ! functions with their lower-l contaminants, do not matter. Only matrices of the dimension of the number
      ! of spinors are diagonalized. CFW^+ SLL CFW = GL + GS (the four-component norm of the input spinors).

      implicit none

      real(8), intent(in)              :: CL(:,:,:), CS(:,:,:)
      real(8), intent(in)              :: SLL(:,:,:), SSS(:,:,:)
      real(8), allocatable, intent(out):: CFW(:,:,:)

      ! Guards. The input must be orthonormal four-component spinors, for which GL + GS = 1 and 0 < GL <= 1.
      ! The deviation from orthonormality is a complete test of the input: a wrong basis, overlap or coefficient
      ! set (e.g. a checkpoint in which the spinors do not live in the space spanned by the overlap matrices
      ! that are passed in) shows up in it. Electronic spinors have GL eigenvalues close to 1 (never below ~0.5,
      ! even for superheavy elements); eigenvalues far below that mean negative-energy (positronic) spinors.
      real(8), parameter               :: NormWarn = 1.D-6, NormStop = 1.D-3
      real(8), parameter               :: GLWarn = 1.D-1, GLStop = 1.D-8
      real(8), allocatable             :: GLew(:), GLcopy(:,:,:)
      real(8), allocatable             :: CLbar(:,:,:), CSbar(:,:,:)
      real(8), allocatable             :: GL(:,:,:), GS(:,:,:), GLsqrt(:,:,:), GLinvsqrt(:,:,:)
      real(8), allocatable             :: W(:,:,:), Wsqrt(:,:,:), T(:,:,:), Id(:,:,:), GLS(:,:,:)
      real(8)                          :: dev
      integer                          :: nL, nS, npos, rcq, i

      nL   = size(CL,1)
      npos = size(CL,2)
      rcq  = size(CL,3)
      nS   = size(CS,1)

      if (size(CS,2) .ne. npos) stop 'construct_fw_coefficients: CL and CS have different number of spinors.'
      if (size(CS,3) .ne. rcq) stop 'construct_fw_coefficients: CL and CS work in different algebra.'
      if (size(SLL,1) .ne. nL .or. size(SLL,2) .ne. nL) &
         stop 'construct_fw_coefficients: SLL dimensions are inconsistent with CL.'
      if (size(SSS,1) .ne. nS .or. size(SSS,2) .ne. nS) &
         stop 'construct_fw_coefficients: SSS dimensions are inconsistent with CS.'
      if (size(SLL,3) .ne. rcq .or. size(SSS,3) .ne. rcq) &
         stop 'construct_fw_coefficients: overlaps and coefficients work in different algebra.'

      if (npos > nL) then
         write(luout,'(A,I0,A,I0,A)') ' FW transformation: more spinors (', npos, &
            ') than large-component basis functions (', nL, ').'
         stop 'construct_fw_coefficients: the spinors cannot be linearly independent.'
      end if

      allocate(CLbar(npos,nL,rcq),CSbar(npos,nS,rcq))
      allocate(GL(npos,npos,rcq),GS(npos,npos,rcq),GLS(npos,npos,rcq))
      allocate(GLsqrt(npos,npos,rcq),GLinvsqrt(npos,npos,rcq))
      allocate(W(npos,npos,rcq),Wsqrt(npos,npos,rcq),T(npos,npos,rcq),Id(npos,npos,rcq))
      allocate(GLew(npos),GLcopy(npos,npos,rcq))
      allocate(CFW(nL,npos,rcq))

      call transpose_conjg(CL,CLbar)
      call transpose_conjg(CS,CSbar)
      call matrix_mul(CLbar,SLL,CL,GL)
      call matrix_mul(CSbar,SSS,CS,GS)

      ! consistency of the input: for orthonormal four-component spinors GL + GS is the unit matrix
      Id = 0.D0
      forall(i=1:npos) Id(i,i,1) = 1.D0
      GLS = GL + GS
      call diff_norm2(GLS,Id,dev)
      if (dev > NormWarn) then
         write(luout,*)
         write(luout,'(A,E12.4)') ' **WARNING** FW transformation: the input spinors are not orthonormal, ||GL+GS-1|| =', dev
      end if
      if (dev > NormStop) then
         write(luout,'(A,E12.4)') ' FW transformation: ||GL+GS-1|| =', dev
         write(luout,*) 'The spinors are not orthonormal in the metric of the given large- and small-component overlap matrices.'
         write(luout,*) 'Check that coefficients and overlaps belong to the same basis (ordering, normalisation, Cartesian versus'
         write(luout,*) 'pure spherical functions) and that the spinors are converged four-component spinors.'
         stop 'construct_fw_coefficients: inconsistent input spinors.'
      end if

      ! conditioning of the large-component Gram matrix (physical eigenvalues, independent of the AO basis)
      GLcopy = GL
      call diag_matrix(GLcopy,GLew)
      if (minval(GLew) < GLStop) then
         write(luout,'(A,E12.4)') ' FW transformation: smallest eigenvalue of the large-component Gram matrix GL =', minval(GLew)
         stop 'construct_fw_coefficients: a spinor has no large component.'
      end if
      if (minval(GLew) < GLWarn) then
         write(luout,*)
         write(luout,'(A,E12.4)') ' **WARNING** FW transformation: small eigenvalue of the large-component Gram matrix GL =', &
              minval(GLew)
         write(luout,*) 'This indicates negative-energy (positronic) spinors among the input spinors.'
      end if

      call hermitian_matrix_sqrt_pair(GL,GLsqrt,GLinvsqrt,'construct_fw_coefficients (large-component Gram matrix)')

      ! W = GL^(-1/2) GS GL^(-1/2),  T = GL^(-1/2) (1+W)^(1/2) GL^(1/2),  CFW = CL T
      call matrix_mul(GLinvsqrt,GS,GLinvsqrt,W)
      Wsqrt = W + Id
      call hermitian_matrix_power(Wsqrt,0.5D0,'construct_fw_coefficients (1+W)')
      call matrix_mul(GLinvsqrt,Wsqrt,GLsqrt,T)
      call matrix_mul(CL,T,CFW)

      if (any(CFW /= CFW)) stop 'construct_fw_coefficients: NaN in the FW coefficients.'

      deallocate(CLbar,CSbar,GL,GS,GLS,GLsqrt,GLinvsqrt,W,Wsqrt,T,Id,GLew,GLcopy)

   end subroutine construct_fw_coefficients

end module foldy_wouthuysen
