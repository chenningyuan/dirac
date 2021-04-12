!#define DEBUG_SOC
!dirac_copyright_start
!      Copyright (c) 2019 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
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
module xamfi_utils

  use picture_change_operations
  use x_fock
  use xamfi_internal_parameters, only: aoo_dfopen, aoo_nopen, aoo_n2tmt,        &
                                       nr_fsym, nr_ao_all, nr_tmo, ioff_aomat_x,&
                                       lutmat_aoo, N2BBASXQ_dim_aoo,            &
                                       ioff_tmot, ioff_tmt, nzt_aoo
  
#ifdef DEBUG_SOC
  use x2c_utils, only: print_x2cmat
#endif


  implicit none

  public get_normalized_2e_soc
  public dump_normalized_2e_soc
  public read_normalized_2e_soc
  public xamfi_get_energy_correction
  public hello_xamfi
  public goodbye_xamfi

  private

contains

!----------------------------------------------------------------------

  subroutine get_normalized_2e_soc(                      &
                                   fmat,                 &
                                   pct_mat,              &
                                   dmat,                 &
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
     integer, allocatable               :: lvalue_array(:)
!----------------------------------------------------------------------

!     normalize (picture-change transform) 2e-SO-fock matrix
#ifdef DEBUG_SOC
      open(99,file='soc-contributions',status='replace',form='formatted',  &
      access='sequential',action="readwrite",position='rewind')

      call print_x2cmat(fmat,                   &
                        nr_ao_total_aoo,        &
                        nr_ao_total_aoo,        &
                        nr_quat,                &
                        aoo_cb_pq_to_uq(1,0),   &
                        'aoosoc - 2e-SOC-unnrm',&
                        99                      &
                       )
#endif

      dmat(1:nr_ao_total_aoo**2 * nr_quat) = 0.0d0
      call perform_pct_saao_bas(fmat,                         &
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
#ifdef DEBUG_SOC
!     debug print
      if(print_lvl > 2)then
        call print_x2cmat(fmat,                   &
                          nr_ao_large_aoo,        &
                          nr_ao_large_aoo,        &
                          nr_quat,                &
                          aoo_cb_pq_to_uq(1,0),   &
                          'aoosoc - 2e-SOC',      &
                          6                       &
                         )
      end if
#endif

#ifdef DEBUG_SOC
      call print_x2cmat(fmat,                   &
                        nr_ao_large_aoo,        &
                        nr_ao_large_aoo,        &
                        nr_quat,                &
                        aoo_cb_pq_to_uq(1,0),   &
                        'aoosoc - 2e-SOC',      &
                        99                      &
                       )
      close(99,status='keep')
#endif



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

!    read 2e-soc from file
     open(fh,file=file_name(1:12),status='old',form='unformatted',      &
          access='sequential',action="readwrite",position='rewind')

     read(fh) nrow_tmp, ncol_tmp, nz_tmp
     if(nrow_tmp >  nrow)then
       print *, 'read_normalized_2e_soc: row dimension does not fit',nrow_tmp,nrow
       call quit( 'read_normalized_2e_soc: row dimension does not fit' )
     end if
     if(ncol_tmp >  ncol)then
       print *, 'read_normalized_2e_soc: col dimension does not fit',ncol_tmp,ncol
       call quit( 'read_normalized_2e_soc: col dimension does not fit' )
     end if
     if(nz_tmp   /= nz   )then
       print *, 'read_normalized_2e_soc: nz  dimension does not fit',nz_tmp
       call quit( 'read_normalized_2e_soc: nz  dimension does not fit' )
     end if

     read(fh) Gmat

     close(fh,status='keep')

  end subroutine read_normalized_2e_soc
!----------------------------------------------------------------------
  function xamfi_get_energy_correction(DMAT,XAMFI_FMAT,      &
                                       DFRAC2,JSHELL,ISHELL, &
                                       NZ,N2BBASX,NFMAT) result(e2tmp)
  
    real*8,  intent(in) :: DMAT(N2BBASX,NZ,NFMAT),XAMFI_FMAT(N2BBASX,NZ,1)
    real*8,  intent(in) :: DFRAC2
    integer, intent(in) :: JSHELL, ISHELL, NZ, N2BBASX, NFMAT                                  
    real*8              :: e2tmp
    real*8,  external   :: ddot
    integer             :: iz

    e2tmp = 0.0d0
    DO iz = 1,NZ
      e2tmp = e2tmp + DFRAC2*DDOT(N2BBASX,DMAT(1,IZ,JSHELL+1),1,XAMFI_FMAT(1,IZ,ISHELL+1),1)
    END DO

  end function xamfi_get_energy_correction
!----------------------------------------------------------------------

  subroutine hello_xamfi()
!**********************************************************************

    print '(/18x,a)', ' *********************************************************************'
    print '( 18x,a)', ' ***                          X-AMFI                               ***'
    print '(18x,a )', ' ***                                                               ***'
    print '( 18x,a)', ' *** a module for the generation of atomic two-electron            ***'
    print '(18x,a )', ' *** scalar- and spin-orbit corrections                            ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** library version: X-AMFI v1.3 (July 2020)                      ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** authors:         - Stefan Knecht (GSI Darmstadt, DE)          ***'
    print '(18x,a )', ' ***                  - Trond Saue ( U Toulouse, F)                ***'
    print '(18x,a )', ' ***                  - Hans Joergen Aa. Jensen (SDU Odense, DK)   ***'
    print '(18x,a )', ' ***                  - Michal Repisky (U Tromsoe, NO)             ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** features:        - atomic two-electron scalar corrections     ***'
    print '(18x,a )', ' ***                  - atomic two-electron spin-orbit corrections ***'
    print '(18x,a )', ' ***                    @ spin-same-orbit                          ***'
    print '(18x,a )', ' ***                    @ spin-same + spin-other-orbit             ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** contact: s.knecht@gsi.de                                      ***'
    print '(18x,a/)', ' *********************************************************************'

  end subroutine hello_xamfi
!**********************************************************************

  subroutine goodbye_xamfi(Z)
!**********************************************************************
    real*8, INTENT(IN) :: Z
    print '(/18x,a)', ' *********************************************************************'
    print '( 18x,a)', ' ***                          X-AMFI                               ***'
    print '( 18x,a)', ' *** Atomic 2e-scalar and spin-orbit corrections are ready.        ***'
    print '( 18x,a)', ' *** The corrections have been stored on the file aoo2esoc.Z       ***'
    print '( 18x,a)', ' *** where Z is the nuclear charge in 3-digit format.              ***'
    print '( 18x,a,i3,a,i0.3,a)', ' *** Example for the present case with Z = ',NINT(Z),': aoo2esoc.',NINT(Z),'       ***'
    print '(18x,a/)', ' *********************************************************************'

  end subroutine goodbye_xamfi
!**********************************************************************

#if !defined(PRG_DIRAC) && !defined(PRG_DALTON)
  subroutine quit(text)
!**********************************************************************

#ifdef VAR_MPI
    use mpi
#endif

    character*(*), inten(in) :: text

#ifdef VAR_MPI
    print *, text
    call MPI_ABORT(-100)
#else
    stop text
#endif

  end subroutine quit
#endif
!**********************************************************************

 end module
