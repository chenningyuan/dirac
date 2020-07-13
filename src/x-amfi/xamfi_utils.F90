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
  use xamfi_internal_parameters, only: aoo_dfopen, aoo_nopen, aoo_n2tmt,        &
                                       nr_fsym, nr_ao_all, nr_tmo, ioff_aomat_x,&
                                       lutmat_aoo, N2BBASXQ_dim_aoo,            &
                                       ioff_tmot, ioff_tmt, nzt_aoo
#ifdef MOD_INTEREST
  use module_interest_interface
#endif
#ifdef DEBUG_SOC
  use x2c_utils, only: print_x2cmat
#endif


  implicit none

  public get_scso_2e_fock
  public get_normalized_2e_soc
  public dump_normalized_2e_soc
  public read_normalized_2e_soc
  public hello_xamfi
  public goodbye_xamfi

  private

contains

!----------------------------------------------------------------------
  subroutine get_scso_2e_fock(                       &
                             dmat,                   &
                             cmo,                    &
                             fmat,                   &
                             nrows,                  &
                             ncols,                  &
                             nr_2e_fock_matrices,    &
                             intflg,                 &
                             mdirac,                 &
                             nz,                     &
                             ipqtoq,                 &
                             ihqmat,                 &
                             lwork_ext,              &
                             print_lvl               &
                             )
!**********************************************************************
#ifdef VAR_MPI
#include "mpif.h"
#endif
     real(8), intent(inout)             :: dmat(nrows,ncols,nz,nr_2e_fock_matrices)
     real(8), intent(inout)             :: cmo(*)
     real(8), intent(out)               :: fmat(nrows,ncols,nz,nr_2e_fock_matrices)
     integer, intent(in)                :: nr_2e_fock_matrices
     integer, intent(in)                :: nrows
     integer, intent(in)                :: ncols
     integer, intent(in)                :: nz
     integer, intent(in)                :: intflg
     integer, intent(in)                :: ipqtoq(4,0:7)
     integer, intent(in)                :: ihqmat(4,-1:1)
     integer, intent(in)                :: lwork_ext
     integer, intent(in)                :: print_lvl
     logical, intent(in)                :: mdirac
!----------------------------------------------------------------------
     integer                            :: i, npos, local_lwork, j, kwork, kfree, lfree
     integer                            :: kdtemp, kdevec, kfmo, kfao, ktmat, ksaomo
     real(8), allocatable               :: work(:)
     integer, allocatable               :: isymop(:), ifckop(:), ihrmop(:), pos(:)
     logical, allocatable               :: saveflags(:)
     logical                            :: parcal
!----------------------------------------------------------------------

!    take care of present "limitations" (being pragmatic for now)
     if(nz /= 4       ) call quit( 'c1 symmetry required for aoosoc module' )

#ifdef BUGGY_AOC_XAMFI
!    !> in order for this routine to work with AOC as well we would need to transform 
     !> the AO Fock matrices to MO basis via MKMOFK and backtransform to AO basis afterwards
     !> FIXME!
     !> For the time being we circumvent it by reusing the Fock matrices from the preceeding 
     !> X2C transformation.
     if(nr_2e_fock_matrices > 1)then

       call get_scso_2e_fock_aoc(                       &
                                 dmat,                  &
                                 fmat,                  &
                                 nrows,                 &
                                 ncols,                 &
                                 1,                     &
                                 nz                     &
                                )


     else
#endif

       local_lwork = lwork_ext - nrows*ncols*nz*nr_2e_fock_matrices - 4 - npos - 3*nr_2e_fock_matrices - 10000
       allocate(work(local_lwork))

       !> 2e-fock matrix construction
       allocate(isymop(nr_2e_fock_matrices))
       allocate(ifckop(nr_2e_fock_matrices))
       allocate(ihrmop(nr_2e_fock_matrices))

       ihrmop = 0
       ifckop = 0
       isymop = 0

       do i = 1,nr_2e_fock_matrices
!        totally symmetric operator
         isymop(i) = 1
!        fock matrix type
         ifckop(i) = 1
!        hermitian operator
         ihrmop(i) = 1
       end do

       allocate(saveflags(4))

       call SaveTaskDistribFlags(saveflags)
       call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. ,.TRUE. /))
       parcal = .false.
#ifdef VAR_MPI
       call mpi_initialized(parcal,i)
#endif
       call SetIntTaskArrayDimension(npos,parcal)
       allocate(pos(npos))


!      2e-fock matrices (Dirac driver routine)
       call twofck(isymop,ihrmop,ifckop,fmat,dmat,nr_2e_fock_matrices,   &
                   pos,intflg,print_lvl,work,local_lwork)

       if(parcal) call SetTaskDistribFlags(saveflags)
       deallocate(pos)
       deallocate(saveflags)
       deallocate(ihrmop)
       deallocate(ifckop)
       deallocate(isymop)

       if(nr_2e_fock_matrices == 1) goto 9999

       !> Add active two-electron Fock matrices to closed-shell Fock matrix
       do i = 1,aoo_nopen
         !> The factor is DF and not D2*DF because of the multiplication with D2 in TWOFCK.
         call daxpy(nrows*ncols*nz,aoo_dfopen(i),fmat(1,1,1,i+1),1,fmat(1,1,1,1),1)
       end do

       !> Transform AO Fock to MO Fock (taking care of coupling of open shells)
       !> ---------------------------------------------------------------------

       KWORK = 1; KFREE = KWORK; LFREE = local_lwork
       CALL MEMGET2('REAL','TMAT' ,KTMAT ,aoo_N2TMT,WORK,KFREE,LFREE)
       CALL MEMGET2('REAL','SAOMO',KSAOMO,aoo_N2TMT,WORK,KFREE,LFREE)

       CALL OPNFIL(LUTMAT_aoo,'AOMOMATaoo','OLD','XAMFI ')

       print *, 'AOMOMAT is open...',lutmat_aoo
       write(6,*) ' dimensions for AOMOMAT reading ... ',aoo_n2tmt

       !> read AO2MO transformation matrix TMAT and LUTMAT rec#2 with SAOMO
       CALL READT(LUTMAT_aoo,aoo_N2TMT,work(ktmat)); CALL READT(LUTMAT_aoo,aoo_N2TMT,work(KSAOMO))
       CLOSE(LUTMAT_aoo,STATUS='KEEP')

       CALL MEMGET2('REAL','DTEMP',KDTEMP,N2BBASXQ_dim_aoo*nr_2e_fock_matrices,WORK,KFREE,LFREE)
       CALL MEMGET2('REAL','DEVEC',KDEVEC,N2BBASXQ_dim_aoo,WORK,KFREE,LFREE)
       CALL MEMGET2('REAL','FMO  ',KFMO  ,N2BBASXQ_dim_aoo,WORK,KFREE,LFREE)
       CALL MKMOFK(work(kFMO),fmat,dmat,work(kTMAT),cmo,work(KSAOMO),WORK(KDTEMP),WORK(KDEVEC),WORK(KFREE),LFREE)

       !> Backtransform MO Fock to AO Fock
       !> --------------------------------

       CALL OPNFIL(LUTMAT_aoo,'AOMOMATaoo','OLD','XAMFI ')
       CALL READT(LUTMAT_aoo,aoo_N2TMT,WORK(KTMAT)); CALL READT(LUTMAT_aoo,aoo_N2TMT,work(KSAOMO))
       CLOSE(LUTMAT_aoo,STATUS='KEEP')

       CALL MEMGET2('REAL','FAO  ',KFAO  ,nrows*ncols*nz,WORK,KFREE,LFREE)

       do i = 1,nr_fsym
         CALL QTRANS('MOAO','S',0.0d0,nr_ao_all(I),nr_ao_all(I),                              &
                     nr_tmo(I),nr_tmo(I),                                                     &
                     work(kfAO+ioff_aomat_x(I,I)),nrows,ncols,NZ,IPQTOQ(1,0),                 &
                     WORK(kfmo  +ioff_tmot(I)),nr_tmo(I)   ,nr_tmo(I),NZ,IPQTOQ(1,0),         &
                     WORK(ksaomo+ioff_tmt(I) ),nr_ao_all(I),nr_tmo(I),nzt_aoo,IPQTOQ(1,0),    &
                     WORK(ksaomo+ioff_tmt(I) ),nr_ao_all(I),nr_tmo(I),nzt_aoo,IPQTOQ(1,0),    &
                     WORK(KFREE),LFREE,30)
       end do

       fmat = 0; call dcopy(nrows*ncols*nz,work(kfao),1,fmat(1,1,1,1),1)

9999   deallocate(work)

#ifdef BUGGY_AOC_XAMFI
     end if
#endif

  end subroutine get_scso_2e_fock
!----------------------------------------------------------------------

  subroutine get_scso_2e_fock_aoc(                        &
                                  dmat,                   &
                                  fmat,                   &
                                  nrows,                  &
                                  ncols,                  &
                                  nr_2e_fock_matrices,    &
                                  nz                      &
                                 )
!**********************************************************************
     real(8), intent(inout)             :: dmat(nrows,ncols,nz,nr_2e_fock_matrices)
     real(8), intent(out)               :: fmat(nrows,ncols,nz,nr_2e_fock_matrices)
     integer, intent(in)                :: nr_2e_fock_matrices
     integer, intent(in)                :: nrows
     integer, intent(in)                :: ncols
     integer, intent(in)                :: nz
!----------------------------------------------------------------------
     integer                            :: i
!----------------------------------------------------------------------

     dmat = 0
     fmat = 0

     !> read the F1 Fock matrix in AO basis from file
     open(81,file='4cf1AO',status='old',form='unformatted',         &
          access='sequential',action='read',position='rewind')
     read(81) dmat
     close(81,status="delete")

     !> read the F1+F2 Fock matrix in AO basis from file
     open(81,file='4cf1f2AO',status='old',form='unformatted',       &
          access='sequential',action='read',position='rewind')
     read(81) fmat
     close(81,status="delete")
 
     !> subtract F1 from full Fock matrix to get F2
     call daxpy(nrows*ncols*nz*nr_2e_fock_matrices,-1.0d0,dmat,1,fmat,1)

  end subroutine get_scso_2e_fock_aoc
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

  subroutine hello_xamfi()
!**********************************************************************

    print '(/18x,a)', ' *********************************************************************'
    print '( 18x,a)', ' ***                          X-AMFI                               ***'
    print '(18x,a )', ' ***                                                               ***'
    print '( 18x,a)', ' *** a module for the generation of atomic two-electron            ***'
    print '(18x,a )', ' *** scalar- and spin-orbit corrections                            ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** library version: X-AMFI v1.2 (June 2016)                      ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** authors:         - Stefan Knecht (ETH Zuerich, CH)            ***'
    print '(18x,a )', ' ***                  - Trond Saue ( U Toulouse, F)                ***'
    print '(18x,a )', ' ***                  - Hans Joergen Aa. Jensen (SDU Odense, DK)   ***'
    print '(18x,a )', ' ***                  - Michal Repisky (U Tromsoe, NO)             ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** features:        - atomic two-electron scalar corrections     ***'
    print '(18x,a )', ' ***                  - atomic two-electron spin-orbit corrections ***'
    print '(18x,a )', ' ***                    @ spin-same-orbit                          ***'
    print '(18x,a )', ' ***                    @ spin-same + spin-other-orbit             ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** contact: stefan.knecht@phys.chem.ethz.ch                      ***'
    print '(18x,a/)', ' *********************************************************************'

  end subroutine hello_xamfi
!**********************************************************************

  subroutine goodbye_xamfi()
!**********************************************************************

    print '(/18x,a)', ' *********************************************************************'
    print '( 18x,a)', ' ***                          X-AMFI                               ***'
    print '( 18x,a)', ' *** Atomic 2e-scalar and spin-orbit corrections are ready.        ***'
    print '( 18x,a)', ' *** The corrections have been stored on the file aoo2esoc.Z where ***'
    print '( 18x,a)', ' *** Z is the nuclear charge in 3-digit format.                    ***'
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
