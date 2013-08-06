!#define DEBUG_SOC
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
  use x2c_utils, only: print_x2cmat


  implicit none

  public get_scso_2e_fock
  public make_so_2e_fock
  public pure_so_2e
  public get_normalized_2e_soc
  public dump_normalized_2e_soc
  public read_normalized_2e_soc
  public hello_aoosoc
  public goodbye_aoosoc

  private

contains

!----------------------------------------------------------------------
  subroutine get_scso_2e_fock(           &
                             dmat,       &
                             fmat,       &
                             nrows,      &
                             ncols,      &
                             nopen,      &
                             aoo_dfopen, &
                             intflg,     &
                             mdirac,     &
                             nz,         &
                             ipqtoq,     &
                             ihqmat,     &
                             naosh_ls,   &
                             print_lvl   &
                             )
!**********************************************************************
     real(8), intent(inout)             :: dmat(nrows,ncols,nz,1+nopen)
     real(8), intent(out)               :: fmat(nrows,ncols,nz,1+nopen)
     real(8), intent(in)                :: aoo_dfopen(nopen)
     integer, intent(in)                :: nopen
     integer, intent(in)                :: nrows
     integer, intent(in)                :: ncols
     integer, intent(in)                :: nz
     integer, intent(in)                :: intflg
     integer, intent(in)                :: ipqtoq(4,0:7)
     integer, intent(in)                :: ihqmat(4,-1:1)
     integer, intent(in)                :: naosh_ls
     integer, intent(in)                :: print_lvl
     logical, intent(in)                :: mdirac
!----------------------------------------------------------------------
     real(8), parameter                 :: hfxfac = 1.0d0
     integer                            :: irep, iz, iq, ih
     integer                            :: i, intflg_local
!----------------------------------------------------------------------

!     take care of present "limitations" (being pragmatic for now)
      if(.not.mdirac   ) call quit( 'mdirac required for aoosoc module' )
      if(nz /= 4       ) call quit( 'c1 symmetry required for aoosoc module' )

      do i = 1, nopen
        call daxpy(nrows*ncols*nz,aoo_dfopen(i),dmat(1,1,1,i+1),1,dmat(1,1,1,1),1)
      end do

!     2e-fock matrix construction

#ifdef MOD_INTEREST
      intflg_local = intflg - mod(intflg,2) ! remove the bit for (LL|LL) integrals

      call process_two_electron_integrals_4comp(               &
                                                fmat(1,1,1,1), &
                                                dmat(1,1,1,1), &
                                                intflg_local,  &
                                                hfxfac         &
                                               )
#else
      call quit( 'InteRest required for aoosoc module' )
#endif

!     Symmetrize FMAT matrix
!     ======================
!     IH = 0   Non-symmetric
!     IH = 1   Symmetric
!     IH = 2   Anti-symmetric

      irep = 0
      do iz = 1,nz
        iq = ipqtoq(iz,irep)
        ih = ihqmat(iq,1)
        if(ih == 1)then
          call fulmat('S',naosh_ls,naosh_ls,fmat(1,1,iz,1))
        else if(ih == 2)then
          call fulmat('A',naosh_ls,naosh_ls,fmat(1,1,iz,1))
        end if
      end do

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
       if(nz /= 4) call quit( '2-el Fock matrix is not in quaternion format' )

       if(nrow /= ncol) call quit( '2-el Fock matrix is not in MDIRAC format' )

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

  subroutine pure_so_2e(                            &
                        Gmat,                       &
                        nrow,                       &
                        ncol,                       &
                        nz                          &
                       )
!----------------------------------------------------------------------
!
!    purpose: retain only two-electron spin-orbit elements in imaginary parts
!
!----------------------------------------------------------------------
     real(8), intent(inout)   :: Gmat(nrow,ncol,4)
     integer, intent(in)      :: nrow
     integer, intent(in)      :: ncol
     integer, intent(in)      :: nz
!----------------------------------------------------------------------
     integer                  :: i, j
!----------------------------------------------------------------------

!      check dimensionalities
       if(nz /= 4) call quit( '2-el soc matrix is not in quaternion format' )

       if(nrow /= ncol) call quit( '2-el soc matrix is not in MDIRAC format' )

!      take out real parts
       do j=1,nrow
         do i=1,nrow
           !> quaternion(0)
           Gmat(i     ,j     ,1) = 0.0d0                 ! ll
         enddo
       enddo

  end subroutine pure_so_2e
!----------------------------------------------------------------------

  subroutine lval_so_2e(                            &
                        Gmat,                       &
                        cart2lval,                  &
                        nrow,                       &
                        ncol,                       &
                        nz                          &
                       )
!----------------------------------------------------------------------
!
!    purpose: retain only two-electron spin-orbit elements in imaginary parts between elements of the same lvalue
!
!----------------------------------------------------------------------
     real(8), intent(inout)   :: Gmat(nrow,ncol,4)
     integer, intent(in)      :: cart2lval(nrow)
     integer, intent(in)      :: nrow
     integer, intent(in)      :: ncol
     integer, intent(in)      :: nz
!----------------------------------------------------------------------
     integer                  :: i, j, iz
!----------------------------------------------------------------------

!      take out offdiagonal l-value parts
     do iz = 2, nz
       do j=1,nrow
         do i=1,ncol
           !> lz(i) /= lz(j)
           if(cart2lval(i) /= cart2lval(j))then
             Gmat(i     ,j     ,iz) = 0.0d0
             Gmat(j     ,i     ,iz) = 0.0d0
           end if
         enddo
       enddo
     enddo

  end subroutine lval_so_2e
!----------------------------------------------------------------------

  subroutine get_normalized_2e_soc(                      &
                                   fmat,                 &
                                   pct_mat,              &
                                   dmat,                 &
                                   nrows,                &
                                   ncols,                &
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

      call get_ss_so_2e_fock(                   &
                             fmat,              &
                             nr_ao_large_aoo,   &
                             nr_ao_large_aoo,   &
                             nr_quat            &
                            )

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
!     pure 2e-SO matrix
      call pure_so_2e(                            &
                      fmat,                       &
                      nr_ao_large_aoo,            &
                      nr_ao_large_aoo,            &
                      nr_quat                     &
                      )

#ifdef DEBUG_SOC
!     get lvalues for cartesian functions
      allocate(lvalue_array(nr_ao_large_aoo))
      lvalue_array(1:nr_ao_large_aoo) = 0.0d0
      call transfer_gto_info(                     &
                             lvalue_array         &
                          )

!     check for l-values in 2e-SO matrix
      call lval_so_2e(                            &
                      fmat,                       &
                      lvalue_array,               &
                      nr_ao_large_aoo,            &
                      nr_ao_large_aoo,            &
                      nr_quat                     &
                      )

      deallocate(lvalue_array)
#endif

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

  subroutine get_ss_so_2e_fock(                            &
                               Gmat,                       &
                               nrow,                       &
                               ncol,                       &
                               nz                          &
                              )
!----------------------------------------------------------------------
!
!    purpose: retain only ss-block of two-electron spin-orbit elements in 2-el Fock matrix
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
       if(nz /= 4) call quit( '2-el Fock matrix is not in quaternion format' )

       if(nrow /= ncol) call quit( '2-el Fock matrix is not in MDIRAC format' )

!      take out ss/sl parts
       do j=1,nrow
         do i=1,nrow
           !> quaternion(i)
!          Gmat(i     ,j     ,2) = Gmat(i     ,j     ,2) ! ll
           Gmat(i     ,j     ,2) = 0.0d0
!          Gmat(i+nrow,j     ,2) = 0.0d0
!          Gmat(i     ,j+nrow,2) = 0.0d0
!          Gmat(i+nrow,j+nrow,2) = 0.0d0
           !> quaternion(j)
!          Gmat(i     ,j     ,3) = Gmat(i     ,j     ,3) ! ll
           Gmat(i     ,j     ,3) = 0.0d0
!          Gmat(i+nrow,j     ,3) = 0.0d0
!          Gmat(i     ,j+nrow,3) = 0.0d0
!          Gmat(i+nrow,j+nrow,3) = 0.0d0
           !> quaternion(k)
!          Gmat(i     ,j     ,4) = Gmat(i     ,j     ,4) ! ll
           Gmat(i     ,j     ,4) = 0.0d0
!          Gmat(i+nrow,j     ,4) = 0.0d0
!          Gmat(i     ,j+nrow,4) = 0.0d0
!          Gmat(i+nrow,j+nrow,4) = 0.0d0
         enddo
       enddo

  end subroutine get_ss_so_2e_fock
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
     if(nz_tmp   /= 4   )then
       print *, 'read_normalized_2e_soc: nz  dimension does not fit',nz_tmp
       call quit( 'read_normalized_2e_soc: nz  dimension does not fit' )
     end if

     read(fh) Gmat

     close(fh,status='keep')

  end subroutine read_normalized_2e_soc
!----------------------------------------------------------------------

  subroutine hello_aoosoc()
!**********************************************************************

    print '(/18x,a)', ' *********************************************************************'
    print '( 18x,a)', ' ***   Preparation of oo-order atomic 2e-spin-orbit corrections    ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** library version:  aoo2eSOC v1.0 (January 2013)                ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** authors:          - Stefan Knecht                             ***'
    print '(18x,a )', ' ***                   - Trond Saue                                ***'
    print '(18x,a )', ' ***                   - Hans Joergen Aagaard Jensen               ***'
    print '(18x,a )', ' *** contributor:      - Michal Repisky                            ***'
    print '(18x,a )', ' *** features:         - atomic 2e-SOC corrections:                ***'
    print '(18x,a )', ' ***                      --> spin-same-orbit                      ***'
    print '(18x,a )', ' ***                     (--> spin-other-orbit will be available   ***'
    print '(18x,a )', ' ***                          in aoo2eSOC v1.1)                    ***'
    print '(18x,a )', ' ***                                                               ***'
    print '( 18x,a)', ' ***        Universities of Odense, Toulouse, and Tromsoe          ***'
    print '(18x,a/)', ' *********************************************************************'

  end subroutine hello_aoosoc
!**********************************************************************

  subroutine goodbye_aoosoc()
!**********************************************************************

    print '(/18x,a)', ' *********************************************************************'
    print '( 18x,a)', ' *** oo-order atomic 2e-spin-orbit corrections have been prepared. ***'
    print '( 18x,a)', ' *** corrections are stored on file aoo2esoc.Z (Z=nuclear charge)  ***'
    print '(18x,a/)', ' *********************************************************************'

  end subroutine goodbye_aoosoc
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
