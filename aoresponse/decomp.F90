module decompMod

   use matrix_module
   use precision

type DecompItem
      !> Logical unit number for DALTON.OUT
      integer     :: lupri
      !> Logical unit number for DALTON.ERR
      integer     :: luerr
   !CONFIGURATIONS:
   !===============
      logical     :: cholesky_decomp
      logical     :: lowdin_diagonalize
      logical     :: lowdin_iterative
      logical     :: lowdin_qiterative
      logical     :: basis123
      logical     :: cfg_do_in_oao
      logical     :: cfg_lcv
      logical     :: cfg_lcm
      logical     :: cfg_lcvbf
      logical     :: lcv_basis
      logical     :: cfg_gcbasis
      logical     :: cfg_unres
      ! Molecule info:
      !> Number of occupied orbitals (if restricted)
      integer     :: nocc
      !> Number of occupied alpha orbitals (if unrestricted)
      integer     :: nocca
      !> Number of occupied beta orbitals (if unrestricted)
      integer     :: noccb
      !> Number of active electrons
      integer     :: nactive
   !SETTINGS:
   !=========
   !Matrices:
   !=========
      type(Matrix) :: lcv_CMO
      type(Matrix) :: U
      type(Matrix) :: U_inv
      type(Matrix) :: U_sav
      type(Matrix) :: U_inv_sav
      type(Matrix) :: DU
      type(Matrix) :: FU
      type(Matrix) :: FUP
      type(Matrix) :: FUQ
      type(Matrix) :: QU
end type DecompItem

contains

subroutine decomp_set_default_config(decomp)
implicit none
   type(decompItem) :: decomp

   !CONFIGURATIONS:
   !===============
      decomp%cholesky_decomp      = .false.
      decomp%lowdin_diagonalize   = .true.
      decomp%lowdin_iterative     = .false.
      decomp%lowdin_qiterative    = .false.
      decomp%basis123             = .false.
      decomp%cfg_do_in_oao        = .false.
      decomp%cfg_lcv              = .false.
      decomp%cfg_lcm              = .false.
      decomp%cfg_lcvbf            = .false.
      decomp%lcv_basis            = .false.
      decomp%cfg_gcbasis          = .false.
      decomp%cfg_unres            = .false.

end subroutine decomp_set_default_config

subroutine chol_save_decomposition(decomp)
use matrix_operations
implicit none
     type(DecompItem)    :: decomp
 
     call mat_init(decomp%U_sav,decomp%U%nrow,decomp%U%ncol)
     call mat_init(decomp%U_inv_sav,decomp%U_inv%nrow,decomp%U_inv%ncol)
                    
     call mat_assign(decomp%U_sav,decomp%U) 
     call mat_assign(decomp%U_inv_sav,decomp%U_inv)
  
end subroutine chol_save_decomposition


subroutine chol_restore_decomposition(decomp)
use matrix_operations
implicit none
     type(DecompItem)    :: decomp

     call mat_assign(decomp%U,decomp%U_sav) 
     call mat_assign(decomp%U_inv, decomp%U_inv_sav)

     call mat_free(decomp%U_sav)
     call mat_free(decomp%U_inv_sav)
  
end subroutine chol_restore_decomposition

   subroutine decomp_init(ndim,decomp)
   use matrix_operations
     implicit none
     integer, intent(in) :: ndim
     type(DecompItem)    :: decomp

     call mat_init(decomp%U,ndim,ndim)
     call mat_init(decomp%U_inv,ndim,ndim)
     call mat_init(decomp%DU,ndim,ndim)
     call mat_init(decomp%FU,ndim,ndim)
     call mat_init(decomp%FUP,ndim,ndim)
     call mat_init(decomp%FUQ,ndim,ndim)
     call mat_init(decomp%QU,ndim,ndim)

   end subroutine decomp_init

   subroutine decomp_shutdown(decomp)
   use matrix_operations
   implicit none
     type(DecompItem)    :: decomp

     call mat_free(decomp%U)
     call mat_free(decomp%U_inv)
     call mat_free(decomp%DU)
     call mat_free(decomp%FU)
     call mat_free(decomp%FUP)
     call mat_free(decomp%FUQ)
     call mat_free(decomp%QU)

   end subroutine decomp_shutdown

   !Get cholesky decomposition of S = Ut*U
   !   or Lowdin decomposition of S = S^1/2 * S^1/2
   subroutine decomposition(S,decomp)
   use matrix_operations
#ifndef PRG_DIRAC
   use lowdin_module
   use lstiming
#endif /* PRG_DIRAC */
   implicit none
         type(Matrix)                              :: S, S_sav
         real(realk),dimension(:), allocatable     :: work1
!Cholesky variables
         real(realk),dimension(:,:), allocatable   :: U_full, U_inv_full !, Ut_inv_full, Ut_full
         integer, dimension(:),allocatable         :: IPVT
         real(realk)                               :: RCOND, dummy(2), tmstart, tmend
         !real(realk)                               :: cpu1, cpu2, wall1, wall2
         integer                                   :: IERR, i, j, fulldim, ludecomp, ndim
         type(DecompItem)                          :: decomp

!  call lstimer('START ',tmstart,tmend,decomp%lupri)
   call qenter('DECOMP')

   ndim = S%nrow
   if (decomp%cfg_unres) then 
      fulldim = 2*ndim
   else
      fulldim = ndim 
   endif

   if (decomp%cfg_unres) then
      if (.not. decomp%cholesky_decomp) then
         write(decomp%lupri,*) 'FALLBACK: open shell currently works only with cholesky decomposition'
         write(decomp%lupri,*) '- Ask Dr. Jansik to fix it for Lowdin!'
         decomp%cholesky_decomp = .true.
      endif
   endif

   if (decomp%lcv_basis) then
      !compute overlap S in lcv basis, save S in zeta_A basis
      call mat_init(S_sav,S%ncol,S%nrow)
      call mat_assign(S_sav,S)
      call mat_mul(decomp%lcv_CMO,S,'t','n',1d0,0d0,decomp%U)
      call mat_mul(decomp%U,decomp%lcv_CMO,'n','n',1d0,0d0,S)
   endif

   if (decomp%cholesky_decomp .or. ndim < 3) then
      !print *, "Preparing to do cholesky decomposition..."
      write (decomp%lupri,*) "Preparing to do cholesky decomposition..."
      allocate(U_full(fulldim,fulldim)) 
      allocate(U_inv_full(fulldim,fulldim))
      allocate(work1(fulldim))
      allocate(IPVT(fulldim))

        !write(lupri,*) 'Overlap S:'
        !call mat_print(S,1,S%nrow,1,S%nrow,lupri)

      U_full = 0.0d0
      call mat_to_full(S,1.0d0,U_full)

        !write (lupri,*) 'Overlap S, full:'
        !call OUTPUT(U_full, 1, fulldim, 1, fulldim, fulldim, fulldim, 1, lupri)

      !Set lower half of U = 0:
      do i = 1, fulldim
         do j = 1, i-1
            U_full(i,j) = 0.0d0
         enddo
      enddo

!     call dchdc(U_full,fulldim,fulldim,work1,0,0,IERR)

      !Invert U and Ut:
      IPVT = 0 ; RCOND = 0.0d0  

      U_inv_full = U_full
      call DGECO(U_inv_full,fulldim,fulldim,IPVT,RCOND,work1)
      call DGEDI(U_inv_full,fulldim,fulldim,IPVT,dummy,work1,01)

      IPVT = 0 ; RCOND = 0.0d0  

      !Convert framework:
      call mat_set_from_full(U_full,1.0d0,decomp%U)
      call mat_set_from_full(U_inv_full,1.0d0,decomp%U_inv)
      
      deallocate(U_full) 
      deallocate(U_inv_full) 
      deallocate(work1)
      deallocate(IPVT)
!     call lstimer('CHOLDC',tmstart,tmend,decomp%lupri)

   else if (decomp%lowdin_diagonalize) then
      !print *, "Preparing to do S^1/2 decomposition..."
      write (decomp%lupri,*) "Preparing to do S^1/2 decomposition..."

      if(decomp%cfg_unres) call quit('lowdin it not tested for unres')

      call chol_lowdin_diag(decomp%lupri,decomp%cfg_unres,S,decomp%U,decomp%U_inv)
!     call lstimer('LWDIAG',tmstart,tmend,decomp%lupri)
   else if (decomp%lowdin_iterative) then
      !print *, "Preparing to do iterative S^1/2 decomposition..."
      write (decomp%lupri,*) "Preparing to do iterative S^1/2 decomposition..."
!     call lowdin_schulz(S,decomp%U,decomp%U_inv,decomp%lupri)

!     call lstimer('LWITER',tmstart,tmend,decomp%lupri)

   else if (decomp%lowdin_qiterative) then
      !print *, "Preparing to do iterative S^1/2 decomposition..."
      write (decomp%lupri,*) "Preparing to do iterative S^1/2 decomposition in quadruple precision..."

! Quadruple precision not supported by gfortran-4.1.1
#ifndef GFORTRAN
      call lowdin_qschulz(S,decomp%U,decomp%U_inv)
#else
      call quit('Qschultz does not support quadruple precision')
#endif

!     call lstimer('LWQITER',tmstart,tmend,decomp%lupri)
   else if (decomp%basis123) then
      write (decomp%lupri,*) ".123BASIS: Decomposition matrices read from disk..."

      ludecomp = -1
      CALL GPOPEN(ludecomp,'U.in','unknown','SEQUENTIAL','UNFORMATTED',i,j)
!     call mat_read_from_disk(ludecomp,decomp%U)
      CALL GPCLOSE(ludecomp,'KEEP')

      ludecomp = -1
      CALL GPOPEN(ludecomp,'U_inv.in','unknown','SEQUENTIAL','UNFORMATTED',i,j)
!     call mat_read_from_disk(ludecomp,decomp%U_inv)
      CALL GPCLOSE(ludecomp,'KEEP')

!     call lstimer('123BAS',tmstart,tmend,decomp%lupri)
   else 
     WRITE(decomp%LUPRI,'(/A)') &
     &     'No type of decomposition specified '
          CALL QUIT(' No type of decomposition specified ')
   endif

    if (decomp%lcv_basis) then
      !compute  U = lc_CMO*(lc_CMO'*S*lc_CMO)^-0.5
      !        iU = U'*S
     
      call mat_assign(S,S_sav)
      call mat_assign(S_sav,decomp%U_inv) 

      call mat_mul(decomp%lcv_CMO,S_sav,'n','n',1d0,0d0,decomp%U_inv)
      call mat_mul(decomp%U_inv,S,'t','n',1d0,0d0,decomp%U)

      call mat_free(S_sav)
      call mat_free(decomp%lcv_CMO)
    endif

   call qexit('DECOMP')
   end subroutine decomposition 

   subroutine res_to_oao_basis(decomp, res, resU)
   use matrix_operations
   implicit none
        type(decompItem),intent(in)  :: decomp
        type(Matrix), intent(in)     :: res
        type(Matrix)                 :: wrk
        type(Matrix)                 :: resU !Intent(out)

   call MAT_INIT(wrk,res%nrow,res%ncol)

   call mat_mul(decomp%U_inv,res,'t','n',1.d0,0.d0,wrk)
   call mat_mul(wrk,decomp%U_inv,'n','n',1.d0,0.d0,resU)

   call mat_free(wrk)
   end subroutine res_to_oao_basis

   subroutine res_from_oao_basis(decomp, resU, res)
   use matrix_operations
   implicit none
        type(decompItem),intent(in)  :: decomp
        type(Matrix), intent(in)     :: resU
        type(Matrix)                 :: wrk
        type(Matrix)                 :: res !Intent(out)

   call MAT_INIT(wrk,res%nrow,res%ncol)

   call mat_mul(decomp%U,resU,'t','n',1.d0,0.d0,wrk)
   call mat_mul(wrk,decomp%U,'n','n',1.d0,0.d0,res)

   call mat_free(wrk)
   end subroutine res_from_oao_basis

   subroutine x_to_oao_basis(decomp, x, xU)  !Used for arh
   use matrix_operations
   implicit none
        type(decompItem),intent(in)  :: decomp
        type(Matrix), intent(in)     :: x
        type(Matrix)                 :: wrk
        type(Matrix)                 :: xU !Intent(out)
                                                                                                                     
   call MAT_INIT(wrk,xU%nrow,xU%ncol)
   call mat_mul(decomp%U,x,'n','n',1.d0,0.d0,wrk)
   call mat_mul(wrk,decomp%U,'n','t',1.d0,0.d0,xU)
   call mat_free(wrk)
   end subroutine x_to_oao_basis

   subroutine x_from_oao_basis(decomp, xU, x)
   use matrix_operations
   implicit none
        type(decompItem),intent(in)  :: decomp
        type(Matrix), intent(in)     :: xU
        type(Matrix)                 :: wrk
        type(Matrix)                 :: x !Intent(out)

   call MAT_INIT(wrk,xU%nrow,xU%ncol)
   call mat_mul(decomp%U_inv,xU,'n','n',1.d0,0.d0,wrk)
   call mat_mul(wrk,decomp%U_inv,'n','t',1.d0,0.d0,x)
   call mat_free(wrk)
   end subroutine x_from_oao_basis

!   subroutine x_to_oao_basis(decomp, x, xU) !Used only for testing
!   implicit none
!        type(decompItem),intent(in)  :: decomp
!        type(Matrix), intent(in)     :: x
!        type(Matrix)                 :: wrk
!        type(Matrix)                 :: xU !Intent(out)
!
!   call MAT_INIT(wrk,xU%nrow,xU%ncol)
!   call mat_mul(decomp%U,x,'n','n',1.d0,0.d0,wrk)
!   call mat_mul(wrk,decomp%U,'n','t',1.d0,0.d0,xU)
!   call mat_free(wrk)
!   end subroutine x_to_oao_basis

   subroutine get_oao_transformed_matrices(decomp,F,D,S)
   use matrix_operations
   implicit none
        type(decompItem)                        :: decomp
        type(Matrix), intent(in), target        :: F, D, S
        type(Matrix)                            :: wrk
        integer                                 :: ndim

      ndim = S%nrow
      !write (lupri,*) "Getting transformed matrices..."
      call mat_init(wrk,ndim,ndim)

      if (decomp%cfg_do_in_oao) then
         decomp%DU = D
         decomp%FU = F
      else
         call mat_mul(decomp%U,D,'n','n',1.d0,0.d0,wrk)
         call mat_mul(wrk,decomp%U,'n','t',1.d0,0.d0,decomp%DU)
         !write(lupri,*) 'DU:'
         !call mat_print(DU,1,FU%nrow,1,FU%ncol,lupri)

         call mat_mul(decomp%U_inv,F,'t','n',1.d0,0.d0,wrk)
         call mat_mul(wrk,decomp%U_inv,'n','n',1.d0,0.d0,decomp%FU)
      endif

!      call mat_add_identity(1.0d0, -1.0d0, decomp%DU, decomp%QU)
      !write(lupri,*) 'FU:'
      !call mat_print(FU,1,FU%nrow,1,FU%ncol,lupri)

       call mat_mul(decomp%DU,decomp%FU,'n','n',1.d0,0.d0,wrk)
       call mat_mul(wrk,decomp%DU,'n','n',1.d0,0.d0,decomp%FUP)
      !write(lupri,*) 'FUP:'
      !call mat_print(FUP,1,FUP%nrow,1,FUP%ncol,lupri)

       call mat_mul(decomp%QU,decomp%FU,'n','n',1.d0,0.d0,wrk)
       call mat_mul(wrk,decomp%QU,'n','n',1.d0,0.d0,decomp%FUQ)
      !write(lupri,*) 'FUQ:'
      !call mat_print(FUQ,1,FUQ%nrow,1,FUQ%ncol,lupri)

   call mat_free(wrk)
   end subroutine get_oao_transformed_matrices

   subroutine project_oao_basis(decomp, X, symm, X_proj)
   use matrix_operations
   implicit none
       type(decompItem),intent(in) :: decomp       
       TYPE(matrix), intent(in) :: X
       integer, intent(in)      :: symm
       TYPE(matrix)             :: scr, PXQ, QXP
       TYPE(matrix)             :: X_proj   !Intent(out)

   call MAT_INIT(PXQ,X%nrow,X%ncol)
   call MAT_INIT(scr,X%nrow,X%ncol)

   !Project in Cholesky basis P(X) = Du X Qu + Qu X Du
   if (symm == 1 .or. symm == 2) then
      call MAT_MUL(decomp%DU,X,'n','n',1.0d0,0.0d0,scr)
      !write(lupri,*) 'DU*x_in, chol_PROJECT_cholbasis:'
      !call mat_print(scr,1,DU%nrow,1,DU%ncol,lupri)

      call MAT_MUL(scr,decomp%QU,'n','n',1.0d0,0.0d0,PXQ)
      !write(lupri,*) 'PXQ, chol_PROJECT_cholbasis:'
      !call mat_print(PXQ,1,DU%nrow,1,DU%ncol,lupri)

      X_proj = PXQ

      !call mat_trans(PXQ, QXP)
      call mat_trans(PXQ, scr)
       if (symm == 2) then 
         !call mat_scal(-1.0d0,QXP)
         call mat_scal(-1.0d0,scr)
      endif
      !write(lupri,*) 'QXP, chol_PROJECT_cholbasis:'
      !call mat_print(QXP,1,DU%nrow,1,DU%ncol,lupri)

      call mat_daxpy(1.0d0,scr,X_proj)
   else
      call MAT_INIT(QXP,X%nrow,X%ncol)

      call MAT_MUL(decomp%DU,X,'n','n',1.0d0,0.0d0,scr)
      call MAT_MUL(scr,decomp%QU,'n','n',1.0d0,0.0d0,PXQ)
      call MAT_MUL(decomp%QU,X,'n','n',1.0d0,0.0d0,scr)
      call MAT_MUL(scr,decomp%DU,'n','n',1.0d0,0.0d0,QXP)

      call MAT_ADD(1.0d0,PXQ,1.0d0,QXP,X_proj)

      call MAT_FREE(QXP)
   endif

   call MAT_FREE(PXQ)
   call MAT_FREE(scr)
   end subroutine project_oao_basis

end module decompMod

subroutine chol_lowdin_diag(lupri,unres,S,U,U_inv)
  use matrix_operations
  use precision
#ifndef PRG_DIRAC
  use lowdin_module
#endif /* PRG_DIRAC */
  implicit none
  integer, intent(in) :: lupri
  logical, intent(in) :: unres
  type(matrix) :: S, U, U_inv
  real(realk), allocatable :: S_full(:,:), Ut_full(:,:), Ut_inv_full(:,:)
  integer dim

  dim = S%ncol

  ! We decompose the alpha part of the matrix which we expect to be spin
  ! independent.
  if(unres) then
     allocate(S_full(dim,2*dim))
!    call mat_to_full2(S, 1d0, S_full)
  else
     allocate(S_full(dim,dim))
     call mat_to_full(S,1.0d0,S_full)
  end if

  allocate(Ut_inv_full(dim,dim), Ut_full(dim,dim))
 !call lowdin_diag(S%nrow,S_full,Ut_full,Ut_inv_full,lupri)
  call mat_set_from_full(Ut_inv_full, 1d0, U_inv)
  call mat_set_from_full(Ut_full,     1d0, U)
  deallocate(S_full, Ut_inv_full, Ut_full)
end subroutine chol_lowdin_diag


