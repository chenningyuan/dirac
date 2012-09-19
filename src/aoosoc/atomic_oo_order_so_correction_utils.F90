!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end
!
!
! this module contains all the utility routines required to add two-electron
! spin-orbit corrections to the hamiltonian in DIRAC-sorted SA-AO basis. 
! the spin-orbit corrections will be obtained as the spin-dependent parts of the 
! (converged) 4c-Fock operator in an atomic Dirac run.
!
! written by sknecht august 2012
!
module atomic_oo_order_so_correction_utils

  use picture_change_operations
  use module_interest_interface

  implicit none

  public get_scso_2e_fock
  public make_so_2e_fock
  public get_normalized_2e_soc
  public dump_normalized_2e_soc
  public read_normalized_2e_soc

  private

contains

!----------------------------------------------------------------------
  subroutine get_scso_2e_fock(           &
                             dmat,       &
                             fmat,       &
                             nrows,      &
                             ncols,      &
                             nr_2e_fock, &
                             intflg,     &
                             mdirac,     &
                             nz,         &
                             print_lvl   &
                             )
!**********************************************************************
     real(8), intent(in)                :: dmat(nrows,ncols,nz,nr_2e_fock)
     real(8), intent(out)               :: fmat(nrows,ncols,nz,nr_2e_fock)
     integer, intent(in)                :: nr_2e_fock
     integer, intent(in)                :: nrows
     integer, intent(in)                :: ncols
     integer, intent(in)                :: nz
     integer, intent(in)                :: intflg
     integer, intent(in)                :: print_lvl
     logical, intent(in)                :: mdirac
!----------------------------------------------------------------------
     real(8), parameter                 :: hfxfac = 1.0d0
!----------------------------------------------------------------------

!     take care of present "limitations" (being pragmatic for now)
      if(.not.mdirac   ) stop 'mdirac required for aoosoc module' 
      if(nz /= 4       ) stop 'c1 symmetry required for aoosoc module' 
      if(nr_2e_fock > 1) stop 'only 1 2e-fock matrix allowed for aoosoc module at present'

!     2e-fock matrix construction
#ifdef MOD_INTEREST
      call process_two_electron_integrals_4comp(            &
                                                fmat,       &
                                                dmat,       &
                                                intflg,     &
                                                hfxfac      &
                                               ) 
#endif

  end subroutine get_scso_2e_fock
!----------------------------------------------------------------------

  subroutine make_so_2e_fock(                            &
                             Gmat,                       & 
                             nrow,                       &
                             ncol,                       &
                             nz                          &
                            )
!----------------------------------------------------------------------
!
!    purpose: retain only two-electron spin-orbit elements in 2-el Fock matrix
!
!----------------------------------------------------------------------
     real(8), intent(inout)   :: Gmat(2*nrow,2*ncol,4)
     integer, intent(in)      :: nrow
     integer, intent(in)      :: ncol
     integer, intent(in)      :: nz
!----------------------------------------------------------------------
     integer                  :: i, j
!----------------------------------------------------------------------

!      check dimensionalities
       if(nz /= 4) stop '2-el Fock matrix is not in quaternion format'
 
       if(nrow /= ncol) stop '2-el Fock matrix is not in MDIRAC format'

!      take out scalar parts
       do j=1,nrow
         do i=1,nrow
           !> quaternion(0) 
           Gmat(i     ,j     ,1) = 0.0d0                 ! ll
           Gmat(i+nrow,j     ,1) = 0.0d0                 ! sl
           Gmat(i     ,j+nrow,1) = 0.0d0                 ! sl
           Gmat(i+nrow,j+nrow,1) = 0.0d0                 ! ss
           !> quaternion(i) 
           Gmat(i     ,j     ,2) = Gmat(i     ,j     ,2) ! ll
           Gmat(i+nrow,j     ,2) = Gmat(i+nrow,j     ,2) ! sl
           Gmat(i     ,j+nrow,2) = Gmat(i     ,j+nrow,2) ! sl
           Gmat(i+nrow,j+nrow,2) = Gmat(i+nrow,j+nrow,2) ! ss
           !> quaternion(j) 
           Gmat(i     ,j     ,3) = Gmat(i     ,j     ,3) ! ll
           Gmat(i+nrow,j     ,3) = Gmat(i+nrow,j     ,3) ! sl
           Gmat(i     ,j+nrow,3) = Gmat(i     ,j+nrow,3) ! sl
           Gmat(i+nrow,j+nrow,3) = Gmat(i+nrow,j+nrow,3) ! ss
           !> quaternion(k) 
           Gmat(i     ,j     ,4) = Gmat(i     ,j     ,4) ! ll
           Gmat(i+nrow,j     ,4) = Gmat(i+nrow,j     ,4) ! sl
           Gmat(i     ,j+nrow,4) = Gmat(i     ,j+nrow,4) ! sl
           Gmat(i+nrow,j+nrow,4) = Gmat(i+nrow,j+nrow,4) ! ss
         enddo
       enddo

  end subroutine make_so_2e_fock
!----------------------------------------------------------------------

  subroutine get_normalized_2e_soc(                      &
                                   fmat,                 &
                                   pct_mat,              &
                                   dmat,                 &
                                   nrows,                &
                                   ncols,                &
                                   nr_2e_fock,           &
                                   nr_ao_total_aoo,      &
                                   nr_ao_large_aoo,      &
                                   nr_ao_all,            &
                                   nr_ao_l,              &
                                   nr_fsym,              &
                                   nr_quat,              &
                                   aoo_cb_pq_to_uq,      &
                                   ioff_aomat_x,         &
                                   aoo_bs_to_fs,         &
                                   print_lvl)
!**********************************************************************
     real(8), intent(inout)             :: dmat(*)
     real(8), intent(in)                :: pct_mat(*)
     real(8), intent(inout)             :: fmat(*)
     integer, intent(in)                :: nrows
     integer, intent(in)                :: ncols
     integer, intent(in)                :: nr_2e_fock
     integer, intent(in)                :: nr_ao_total_aoo
     integer, intent(in)                :: nr_ao_large_aoo
     integer, intent(in)                :: nr_ao_all(nr_fsym)
     integer, intent(in)                :: nr_ao_l(nr_fsym)
     integer, intent(in)                :: nr_fsym
     integer, intent(in)                :: nr_quat
     integer, intent(in)                :: ioff_aomat_x(2,2)
     integer, intent(in)                :: aoo_cb_pq_to_uq(4, 0:7)
     integer, intent(in)                :: aoo_bs_to_fs(0:7,2)
     integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
     integer                            :: fock_mat_offset
     integer                            :: i
!----------------------------------------------------------------------

!     normalize (picture-change transform) 2e-SO-fock matrix
      fock_mat_offset = 1

      do i = 1, nr_2e_fock

        dmat(1:nr_ao_total_aoo**2 * nr_quat) = 0.0d0
        call perform_pct_saao_bas(fmat(fock_mat_offset),        &
                                  pct_mat,                      &
                                  dmat,                         &
                                  nr_ao_total_aoo,              &
                                  nr_ao_large_aoo,              &
                                  nr_ao_all,                    &
                                  nr_ao_l,                      &
                                  nr_fsym,                      &
                                  nr_quat,                      &
                                  aoo_cb_pq_to_uq,              &
                                  ioff_aomat_x,                 &
                                  0,                            &
                                  aoo_bs_to_fs,                 &
                                  .true.,                       &
                                  print_lvl)

        fock_mat_offset = fock_mat_offset + nrows*ncols*nr_quat

      end do
  
!     todo: combine fock matrices in case nr_2e_fock > 1

  end subroutine get_normalized_2e_soc
!----------------------------------------------------------------------

  subroutine dump_normalized_2e_soc(                            &
                                    Gmat,                       & 
                                    nrow,                       &
                                    ncol,                       &
                                    nz,                         &
                                    file_base,                  &
                                    element,                    &
                                    fh                          &
                                   )
!----------------------------------------------------------------------
!
!    purpose: dump two-electron spin-orbit corrections (normalized)
!
!----------------------------------------------------------------------
     real(8), intent(inout)        :: Gmat(nrow,ncol,4)
     integer, intent(in)           :: nrow
     integer, intent(in)           :: ncol
     integer, intent(in)           :: nz
     integer, intent(in)           :: fh
     integer, intent(in)           :: element
     character (len=8), intent(in) :: file_base
!----------------------------------------------------------------------
     integer                       :: i
     character (len=3)             :: file_extension
     character (len=12)            :: file_name
!----------------------------------------------------------------------

!    determine full file name
     write(file_extension,'(i3)') element
     do i = 1, 3
        if(file_extension(i:i) .eq. ' ') file_extension(i:i) = '0'
     end do
     write(file_name,'(a8,a1,a3)') file_base,'.',file_extension

!    dump 2e-soc on file
     open(fh,file=file_name(1:12),status='replace',form='unformatted',      &
          access='sequential',action="readwrite",position='rewind')

     write(fh) nrow, ncol, nz
     write(fh) Gmat
     close(fh,status='keep')

  end subroutine dump_normalized_2e_soc
!----------------------------------------------------------------------

  subroutine read_normalized_2e_soc(                            &
                                    Gmat,                       & 
                                    nrow,                       &
                                    ncol,                       &
                                    nz,                         &
                                    file_base,                  &
                                    element,                    &
                                    fh                          &
                                   )
!----------------------------------------------------------------------
!
!    purpose: read two-electron spin-orbit corrections (normalized)
!
!----------------------------------------------------------------------
     real(8), intent(inout)        :: Gmat(nrow,ncol,4)
     integer, intent(in)           :: nrow
     integer, intent(in)           :: ncol
     integer, intent(in)           :: nz
     integer, intent(in)           :: fh
     integer, intent(in)           :: element
     character (len=8), intent(in) :: file_base
!----------------------------------------------------------------------
     integer                       :: i
     integer                       :: ncol_tmp
     integer                       :: nrow_tmp
     integer                       :: nz_tmp
     character (len=3)             :: file_extension
     character (len=12)            :: file_name
!----------------------------------------------------------------------

!    determine full file name
     write(file_extension,'(i3)') element
     do i = 1, 3
        if(file_extension(i:i) .eq. ' ') file_extension(i:i) = '0'
     end do
     write(file_name,'(a8,a1,a3)') file_base,'.',file_extension

!    dump 2e-soc on file
     open(fh,file=file_name(1:12),status='old',form='unformatted',      &
          access='sequential',action="readwrite",position='rewind')

     read(fh) nrow_tmp, ncol_tmp, nz_tmp
     if(nrow_tmp >  nrow) stop 'read_normalized_2e_soc: row dimension does not fit'
     if(ncol_tmp >  ncol) stop 'read_normalized_2e_soc: col dimension does not fit'
     if(nz_tmp   /= 4   ) stop 'read_normalized_2e_soc: nz  dimension does not fit'

     read(fh) Gmat

     close(fh,status='keep')

  end subroutine read_normalized_2e_soc
!----------------------------------------------------------------------

 end module
