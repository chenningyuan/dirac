!#define DEBUG_SOC
!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
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
                                       ioff_tmot, ioff_tmt, nzt_aoo, amfh5Data
  use amfH5Interface
  
#ifdef DEBUG_SOC
  use x2c_utils, only: print_x2cmat
#endif


  implicit none

  public get_normalized_2e_soc
  public dump_normalized_2e_soc
  public dump_Xc_c1_dens
  public dump_delta_xc_energy
  public read_xamfi_general_matrix
  public read_xamfi_general_scalar
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
                                    fh,                         &
                                    final_dump                  &
                                   )
!----------------------------------------------------------------------
!
!    purpose: dump two-electron spin-orbit corrections (normalized)
!
!----------------------------------------------------------------------
     use xamfi_internal_parameters, only: aoo_cb_pq_to_uq,  aoo_bs_to_fs,&
                                          aoo_bs_irrep_mat, aoo_pointer_quat,&
                                          aoo_iqmult_trip_q, aoo_cb_uq_to_pq
#ifdef DEBUG_SOC
     use x2c_utils, only: print_x2cmat
#endif
     real(8), intent(inout)        :: Gmat(nrow,ncol,nz)
     integer, intent(in)           :: nrow
     integer, intent(in)           :: ncol
     integer, intent(in)           :: nz
     integer, intent(in)           :: fh
     integer, intent(in)           :: element
     character (len=8), intent(in) :: file_base
     logical, intent(in)           :: final_dump
!----------------------------------------------------------------------
     integer                       :: i, irep, iz, iq, ipar, irepd, ipq, j, k
     character (len=3)             :: file_extension
     character (len=12)            :: file_name
     real*8, ALLOCATABLE           :: Gmat_c1(:,:,:)
     integer, PARAMETER            :: nzc1 = 4
!----------------------------------------------------------------------

    if(final_dump)then
        allocate(Gmat_c1(nrow,ncol,nzc1)); Gmat_c1 = 0
#ifdef DEBUG_SOC
        call print_x2cmat(Gmat, &
        nrow,                   &
        ncol,                   &
        nz,                     &
        aoo_cb_pq_to_uq(1,0),   &
        'aoosoc - 2e-SOC 0',    &
        6                       &
       )
#endif
        ! a. insert phases
        IREP=0
        IF (NZ < 4) THEN
            DO  IZ = 1, NZ
                IQ = aoo_cb_pq_to_uq(IZ,IREP)
                CALL Q2BPHASE('F',IQ,1,Gmat(1,1,iz))
            END DO
        END IF
#ifdef DEBUG_SOC
       call print_x2cmat(Gmat,                  &
                        nrow,                   &
                        ncol,                   &
                        nz,                     &
                        aoo_cb_pq_to_uq(1,0),   &
                        'aoosoc - 2e-SOC a',    &
                        6                       &
                       )
#endif

        ! b. transform to unsorted SO-AO basis
        call bstobu_no_work(Gmat,nz)

#ifdef DEBUG_SOC
        call print_x2cmat(Gmat, &
        nrow,                   &
        ncol,                   &
        nz,                     &
        aoo_cb_pq_to_uq(1,0),   &
        'aoosoc - 2e-SOC b',    &
        6                       &
       )
#endif
        ! c. transform from unsorted SO-AO basis to AO (C1* symmetry) basis
        IPAR = aoo_bs_to_fs(IREP,1)
        DO IZ = 1, nzc1
           IREPD = aoo_bs_irrep_mat(IZ,IREP)
           IQ    = aoo_iqmult_trip_q(1,aoo_pointer_quat(IREPD,IPAR),IZ)
           IPQ   = aoo_cb_uq_to_pq(IQ,IREP)
           CALL MTSOAO(Gmat(1,1,IPQ),     &
                       Gmat_c1(1,1,IZ),   &
                       nrow,IREPD,0)
        END DO


        ! get rid of numerical noise
        do iz = 1, nzc1
          do k = 1, ncol
            do j= 1, nrow
              if(abs(Gmat_c1(j,k,iz)) < 1.0d-12) Gmat_c1(j,k,iz) = 0.0d0
            end do
          end do
        end do

#ifdef DEBUG_SOC
        call print_x2cmat(Gmat_c1, &
        nrow,                      &
        ncol,                      &
        nzc1,                      &
        aoo_cb_pq_to_uq(1,0),      &
        'aoosoc - 2e-SOC c',       &
        6                          &
       )
#endif
    endif !> final_dump

!> HDF5 data dump
if(final_dump)then
        call amf_write_hdf5(amfh5Data(1),Gmat_c1,nrow,ncol,nzc1,file_base)
        DEALLOCATE(Gmat_c1)
else
        call amf_write_hdf5(amfh5Data(1),Gmat,nrow,ncol,nz,file_base)
end if


  end subroutine dump_normalized_2e_soc
!----------------------------------------------------------------------
!#define DEBUG_SOC
  subroutine dump_Xc_c1_dens(   &
    Dmat,                       &
    nrow,                       &
    ncol,                       &
    nz,                         &
    nr_2e_fock_matrices,        &
    x_nopen,                    &
    x_dfopen,                   &
    file_base,                  &
    element,                    &
    fh,                         &
    final_dump                  &
   )
!----------------------------------------------------------------------
!
!    purpose: dump (atomic Xc [X=2,4]) density matrix in C1 symmetry
!
!----------------------------------------------------------------------
use xamfi_internal_parameters, only: aoo_cb_pq_to_uq,  aoo_bs_to_fs,&
          aoo_bs_irrep_mat, aoo_pointer_quat,&
          aoo_iqmult_trip_q, aoo_cb_uq_to_pq
#ifdef DEBUG_SOC
use x2c_utils, only: print_x2cmat
#endif
real*8 , intent(inout)        :: Dmat(nrow,ncol,nz,nr_2e_fock_matrices)
integer, intent(in)           :: nrow
integer, intent(in)           :: ncol
integer, intent(in)           :: nz
integer, intent(in)           :: nr_2e_fock_matrices
integer, INTENT(IN)           :: x_nopen
real*8 , INTENT(IN)           :: x_dfopen(0:x_nopen)
integer, intent(in)           :: fh
integer, intent(in)           :: element
character (len=8), intent(in) :: file_base
logical, intent(in)           :: final_dump
!----------------------------------------------------------------------
integer                       :: i, irep, iz, iq, ipar, irepd, ipq
character (len=3)             :: file_extension
character (len=12)            :: file_name
real*8, ALLOCATABLE           :: Dmat_c1(:,:,:)
integer, PARAMETER            :: nzc1 = 4
!----------------------------------------------------------------------


!> Start by constructing a fractional occupation density (eqivalent to what is used in an atomic start guess)
!> that is, add active density matrices to closed-shell density matrix

if(nr_2e_fock_matrices > 1)then
do i = 1,x_nopen
  call daxpy(nrow*ncol*nz,x_dfopen(i),Dmat(1,1,1,i+1),1,Dmat(1,1,1,1),1)
end do
endif

if(final_dump)then
allocate(Dmat_c1(nrow,ncol,nzc1)); Dmat_c1 = 0
#ifdef DEBUG_SOC
call print_x2cmat(Dmat, &
nrow,                   &
ncol,                   &
nz,                     &
aoo_cb_pq_to_uq(1,0),   &
'aoosoc - DMAT 0',      &
6                       &
)
#endif

! a. insert phases
IREP=0
IF (NZ < 4) THEN
DO  IZ = 1, NZ
IQ = aoo_cb_pq_to_uq(IZ,IREP)
CALL Q2BPHASE('D',IQ,1,Dmat(1,1,iz,1))
END DO
END IF
#ifdef DEBUG_SOC
call print_x2cmat(Dmat,                  &
nrow,                   &
ncol,                   &
nz,                     &
aoo_cb_pq_to_uq(1,0),   &
'aoosoc - DMAT a',      &
6                       &
)
#endif

! b. transform to unsorted SO-AO basis
call bstobu_no_work(Dmat,nz)

#ifdef DEBUG_SOC
call print_x2cmat(Dmat, &
nrow,                   &
ncol,                   &
nz,                     &
aoo_cb_pq_to_uq(1,0),   &
'aoosoc - DMAT b',      &
6                       &
)
#endif

! c. transform from unsorted SO-AO basis to AO (C1* symmetry) basis
IPAR = aoo_bs_to_fs(IREP,1)
DO IZ = 1, nzc1
IREPD = aoo_bs_irrep_mat(IZ,IREP)
IQ    = aoo_iqmult_trip_q(1,aoo_pointer_quat(IREPD,IPAR),IZ)
IPQ   = aoo_cb_uq_to_pq(IQ,IREP)
CALL DTSOAO(Dmat(1,1,IPQ,1),     &
Dmat_c1(1,1,IZ),   &
nrow,IREPD,0)
END DO

#ifdef DEBUG_SOC
call print_x2cmat(Dmat_c1, &
nrow,                      &
ncol,                      &
nzc1,                      &
aoo_cb_pq_to_uq(1,0),      &
'aoosoc - DMAT c',         &
6                          &
)
#undef DEBUG_SOC
#endif

end if

!> HDF5 data dump
if(final_dump)then
call amf_write_hdf5(amfh5Data(1),Dmat_c1,nrow,ncol,nzc1,file_base)
DEALLOCATE(Dmat_c1)
else
call amf_write_hdf5(amfh5Data(1),Dmat(1,1,1,1),nrow,ncol,nz,file_base)
end if


end subroutine dump_Xc_c1_dens
!----------------------------------------------------------------------

  subroutine dump_delta_xc_energy(                           &
                                  Gscalar,                   &
                                  file_base,                 &
                                  element,                   &
                                  fh                         &
                                  )
!----------------------------------------------------------------------
!
!    purpose: dump DELTA of XC energy
!
!----------------------------------------------------------------------
     real(8), intent(in)           :: Gscalar
     integer, intent(in)           :: fh
     integer, intent(in)           :: element
     character (len=8), intent(in) :: file_base
!----------------------------------------------------------------------
     integer                       :: i
     character (len=3)             :: file_extension
     character (len=12)            :: file_name
!----------------------------------------------------------------------

    !> HDF5 data dump
    call amf_write_hdf5(amfh5Data(1),Gscalar,1,1,1,file_base)

  end subroutine dump_delta_xc_energy
!----------------------------------------------------------------------

  subroutine read_xamfi_general_matrix(                            &
                                       Gmat,                       &
                                       nrow,                       &
                                       ncol,                       &
                                       nzc1,                       &
                                       file_base,                  &
                                       element,                    &
                                       fh,                         &
                                       atom_id                     &
                                      )
!----------------------------------------------------------------------
!
!    purpose: read two-electron spin-orbit corrections (normalized)
!
!----------------------------------------------------------------------
     real(8), intent(inout)        :: Gmat(nrow,ncol,nzc1)
     integer, intent(in)           :: nrow
     integer, intent(in)           :: ncol
     integer, intent(in)           :: nzc1
     integer, intent(in)           :: fh
     integer, intent(in)           :: element, atom_id
     character (len=8), intent(in) :: file_base
!----------------------------------------------------------------------
     integer                       :: i
     integer                       :: ncol_tmp
     integer                       :: nrow_tmp
     integer                       :: nz_tmp
     character (len=3)             :: file_extension
     character (len=12)            :: file_name
!----------------------------------------------------------------------

     call amf_read_hdf5(amfH5Data(atom_id),Gmat,file_base)

  end subroutine read_xamfi_general_matrix
!----------------------------------------------------------------------

  function read_xamfi_general_scalar(                            &
                                     file_base,                  &
                                     element,                    &
                                     fh,                         &
                                     atom_id                     &
                                    )result(Gscalar)
!----------------------------------------------------------------------
!
!    purpose: read some scalar value from file
!
!----------------------------------------------------------------------
     integer, intent(in)           :: fh
     integer, intent(in)           :: element, atom_id
     character (len=8), intent(in) :: file_base
!----------------------------------------------------------------------
     real*8                        :: Gscalar
     integer                       :: i
     character (len=3)             :: file_extension
     character (len=12)            :: file_name
!----------------------------------------------------------------------

     call amf_read_hdf5(amfH5Data(atom_id),Gscalar,file_base)

  end function read_xamfi_general_scalar
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
    print '(18x,a )', ' *** library version: X-AMFI v1.5 (February 2023)                  ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** authors:         - Stefan Knecht (Algorithmiq, FI)            ***'
    print '(18x,a )', ' ***                  - Trond Saue ( U Toulouse, F)                ***'
    print '(18x,a )', ' ***                  - Hans Joergen Aa. Jensen (SDU Odense, DK)   ***'
    print '(18x,a )', ' ***                  - Michal Repisky (U Tromsoe, NO)             ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** features:        - atomic two-electron scalar corrections     ***'
    print '(18x,a )', ' ***                  - atomic two-electron spin-orbit corrections ***'
    print '(18x,a )', ' ***                    @ spin-same-orbit                          ***'
    print '(18x,a )', ' ***                    @ spin-same + spin-other-orbit             ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** contact: stefan.knecht@gmail.com                              ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** If results obtained with this code are published, please cite ***'
    print '(18x,a )', ' ***                                                               ***'
    print '(18x,a )', ' *** S. Knecht, M. Repisky, H. J. Aa. Jensen, T.Saue,              ***'
    print '(18x,a )', ' *** J. Chem. Phys., 157, 114106 (2022);                           ***'
    print '(18x,a )', ' *** https://doi.org/10.1063/5.0095112                             ***'
    print '(18x,a )', ' ***                                                               ***'
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

    character*(*), intent(in) :: text

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
