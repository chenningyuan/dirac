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
! this module contains all the functionality required to add one- and two-electron
! spin-orbit corrections to the hamiltonian in DIRAC-sorted SA-AO basis.
! the spin-orbit corrections will be obtained as the spin-dependent parts of the
! (converged) 4c-Fock operator in an atomic Dirac run.
!
! written by sknecht august 2012
!
module xamfi_driver

  use x_fock
  use xamfi_global_parameters
  use xamfi_utils
  use xamfi_internal_parameters

!#define DEBUG_SOC
#ifdef  DEBUG_SOC
  use x2c_utils, only: print_x2cmat
#endif

  implicit none

  public get_xamfi_correction
  public put_xamfi_correction
  public put_atomic_2e_mean_field_correction
  public put_atomic_density_correction
  public collect_delta_xc_energy

  private

contains

!----------------------------------------------------------------------
  subroutine get_xamfi_correction(            &
                                  dmat,       &
                                  cmo,        &
                                  pct_mat,    &
                                  calc_2esoc, &
                                  isDFT,      &
                                  isLoA,      &
                                  lwork_ext,  &
                                  print_lvl   &
                                  )
!**********************************************************************
     real(8), intent(inout)             :: dmat(nrows,ncols,nr_quat,nr_2e_fock_matrices)
     real(8), intent(inout)             :: cmo(*)
     real(8), intent(in)                :: pct_mat(*)
     logical, intent(in)                :: calc_2esoc
     logical, intent(in)                :: isDFT
     logical, intent(in)                :: isLoA
     integer, intent(in)                :: lwork_ext
     integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
     real(8), allocatable               :: fmat(:), tmp(:)
!----------------------------------------------------------------------

      if(aooeamf)then
        !> save atomic density matrix (of atom(1) ... == atomic fragment)
        allocate(tmp(nrows*ncols*nr_quat*nr_2e_fock_matrices)); tmp = 0
        call dcopy(nrows*ncols*nr_quat*nr_2e_fock_matrices,dmat,1,tmp,1)
        call dump_Xc_c1_dens(                            &
                             tmp,                        &
                             nrows,                      &
                             ncols,                      &
                             nr_quat,                    &
                             nr_2e_fock_matrices,        &
                             aoo_nopen,                  &
                             aoo_dfopen(0:aoo_nopen),    &
                             'aooDM_4c',                 &
                             nint(atom(1)%charge),       &
                             102,                        &
                             .true.                      &
                             )
        deallocate(tmp);
      endif

      if(.not. calc_2esoc) return

      call hello_xamfi()

      allocate(fmat(nrows*ncols*nr_quat*nr_2e_fock_matrices)); fmat = 0

      !> get full (scalar+spin-orbit) 2e-fock matrix
      call get_scso_2e_fock(                                              &
                            dmat                  = dmat,                 &
                            cmo                   = cmo,                  &
                            fmat                  = fmat,                 &
                            nrows                 = nrows,                &
                            ncols                 = ncols,                &
                            nr_2e_fock_matrices   = nr_2e_fock_matrices,  &
                            intflg                = aoo_intflg,           &
                            nz                    = nr_quat,              &
                            ipqtoq                = aoo_cb_pq_to_uq(1,0), &
                            x_nopen               = aoo_nopen,            &
                            x_dfopen              = aoo_dfopen(0:aoo_nopen),&
                            x_n2tmt               = aoo_n2tmt,            &
                            nzt_x                 = nzt_aoo,              &
                            nr_fsym               = nr_fsym,              &
                            N2BBASXQ_dim_x        = N2BBASXQ_dim_aoo,     &
                            nr_tmo                = nr_tmo,               &
                            nr_ao_all             = nr_ao_all,            &
                            nr_mo_all             = dim_oshell,           &
                            ioff_tmot             = ioff_tmot,            &
                            ioff_tmt              = ioff_tmt,             &
                            ioff_aomat_x          = ioff_aomat_x,         &
                            myAOMOMAT             = 'AOMOMATaoo',         &
                            isDFT                 = isDFT,                &
                            isLoA                 = isLoA,                &
                            xcmode                = 4,                    &
                            lwork_ext             = lwork_ext,            &
                            print_lvl             = print_lvl             &
                           )
#ifdef DEBUG_SOC
      !> debug print
      if(print_lvl > 0)then
        call print_x2cmat(                             &
                          fmat,                        &
                          nr_ao_total_aoo,             &
                          nr_ao_total_aoo,             &
                          nr_quat,                     &
                          aoo_cb_pq_to_uq(1,0),        &
                          'X-AMFI - 2e-complete',      &
                          6                            &
                         )
      end if
#endif

      if(aoo2ediag)then
        !> save atomic 2e-Fock on file (of atom(1) ... == atomic fragment)
        !> never mind the confusing naming of the routine that we call...
        allocate(tmp(nrows*ncols*nr_quat*nr_2e_fock_matrices)); tmp = 0
        call dcopy(nrows*ncols*nr_quat*nr_2e_fock_matrices,fmat,1,tmp,1)
        call dump_normalized_2e_soc(                            &
                                    tmp,                        &
                                    nrows,                      &
                                    ncols,                      &
                                    nr_quat,                    &
                                    'aoo2e_4c',                 &
                                    nint(atom(1)%charge),       &
                                    101,                        &
                                    .true.                      &
                                    )
        aoo2ediag = .false.; deallocate(tmp)
      end if

      !> normalize (picture-change transform) 2e-SO-fock matrix
      call get_normalized_2e_soc(                      &
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
                                 print_lvl             &
                                )

#ifdef DEBUG_SOC
        !> debug print
        call print_x2cmat(                             &
                          fmat,                        &
                          nr_ao_large_aoo,             &
                          nr_ao_large_aoo,             &
                          nr_quat,                     &
                          aoo_cb_pq_to_uq(1,0),        &
                          'X-AMFI - 2e-normalized',    &
                          6                            &
                         )
#endif

      !> save spin-orbit parts on file (of atom(1) ... == atomic fragment)
      call dump_normalized_2e_soc(                            &
                                  fmat,                       &
                                  nr_ao_large_aoo,            &
                                  nr_ao_large_aoo,            &
                                  nr_quat,                    &
                                  file_name_aoo_contributions,&
                                  nint(atom(1)%charge),       &
                                  aoo_fsoc,                   &
                                  .false.                     &
                                 )

      deallocate(fmat)

      !call goodbye_xamfi()

  end subroutine get_xamfi_correction
!----------------------------------------------------------------------

  subroutine put_xamfi_correction(                             &
                                  Gmat,                        &
                                  naosh_l,                     &
                                  nrows_ext,                   &
                                  ncols_ext,                   &
                                  nr_quat_ext,                 &
                                  center_id,                   &
                                  center_total,                &
                                  print_lvl                    &
                                 )
!**********************************************************************
     real(8), intent(inout)             :: Gmat(*)
     integer, intent(in)                :: naosh_l
     integer, intent(in)                :: nrows_ext
     integer, intent(in)                :: ncols_ext
     integer, intent(in)                :: nr_quat_ext
     integer, intent(in)                :: center_id
     integer, intent(inout)             :: center_total
     integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
     real(8), allocatable               :: Gmat_frag(:)
     real(8), allocatable               :: rbuf(:)
     real(8), allocatable               :: cbuf(:)
     real(8), allocatable               :: ibuf(:)
     real(8), allocatable               :: jbuf(:)
     INTEGER, parameter                 :: nzc1 = 4
!----------------------------------------------------------------------


      if(nint(atom(center_total)%charge) > 1)then

        allocate(Gmat_frag(nrows_ext*ncols_ext*nzc1)); Gmat_frag = 0
!       !> get spin-orbit parts from file (of atom(i) ... == atomic fragment)
        call read_xamfi_general_matrix(                                 &
                                       Gmat_frag,                       &
                                       nrows_ext,                       &
                                       ncols_ext,                       &
                                       nzc1,                            &
                                       file_name_aoo_contributions,     &
                                       nint(atom(center_total)%charge), &
                                       aoo_fsoc                         &
                                      )

#ifdef DEBUG_SOC
      !> debug print
      if(print_lvl > 5)then
        call print_x2cmat(                                       &
                          Gmat_frag,                             &
                          nrows_ext,                             &
                          ncols_ext,                             &
                          nzc1,                                  &
                          aoo_cb_pq_to_uq(1,0),                  &
                          'X-AMFI - corrections from file',      &
                          6                                      &
                         )
      end if
#endif


        allocate(rbuf(naosh_l))
        allocate(cbuf(naosh_l))
        allocate(ibuf(nrows_ext))
        allocate(jbuf(ncols_ext))

        rbuf = 0
        cbuf = 0
        ibuf = 0
        jbuf = 0

        call put_aoblock(Gmat,                            &
                         naosh_l,                         &
                         naosh_l,                         &
                         center_total,                    &
                         aoo_nont(center_id),             &
                         aoo_nucdeg(center_total),        &
                         nzc1,                            &
                         Gmat_frag,                       &
                         nrows_ext,                       &
                         nrows_ext,                       &
                         -1,                              &
                         ncols_ext,                       &
                         ncols_ext,                       &
                         -1,                              &
                         rbuf,                            &
                         cbuf,                            &
                         ibuf,                            &
                         jbuf)

#ifdef DEBUG_SOC
      !> debug print
      if(print_lvl > 5)then
        call print_x2cmat(                                       &
                          Gmat,                                  &
                          naosh_l,                               &
                          naosh_l,                               &
                          nr_quat_ext,                           &
                          aoo_cb_pq_to_uq(1,0),                  &
                          'X-AMFI - corrections accumulated',    &
                          6                                      &
                         )
      end if
#endif
        deallocate(jbuf)
        deallocate(ibuf)
        deallocate(cbuf)
        deallocate(rbuf)
        deallocate(Gmat_frag)
        !write(6,*) 'done correction for center ',center_total,aoo_nont(center_id),center_id,aoo_nucdeg(center_total)
        !call flshfo(6)

      end if

!     update counter
      center_total = center_total + aoo_nont(center_id)

  end subroutine put_xamfi_correction
!----------------------------------------------------------------------
  subroutine put_atomic_2e_mean_field_correction(&
    Gmat,                        &
    naosh_ls,                    &
    nrows_ext,                   &
    ncols_ext,                   &
    nr_quat_ext,                 &
    center_id,                   &
    center_total,                &
    print_lvl                    &
   )
!**********************************************************************
real(8), intent(inout)             :: Gmat(*)
integer, intent(in)                :: naosh_ls
integer, intent(in)                :: nrows_ext
integer, intent(in)                :: ncols_ext
integer, intent(in)                :: nr_quat_ext
integer, intent(in)                :: center_id
integer, intent(inout)             :: center_total
integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
real(8), allocatable               :: Gmat_frag(:)
real(8), allocatable               :: rbuf(:)
real(8), allocatable               :: cbuf(:)
real(8), allocatable               :: ibuf(:)
real(8), allocatable               :: jbuf(:)
INTEGER, parameter                 :: nzc1 = 4
!----------------------------------------------------------------------


if(nint(atom(center_total)%charge) > 1)then

allocate(Gmat_frag(nrows_ext*ncols_ext*nzc1)); Gmat_frag = 0
!       !> get spin-orbit parts from file (of atom(i) ... == atomic fragment)

call read_xamfi_general_matrix(               &
             Gmat_frag,                       &
             nrows_ext,                       &
             ncols_ext,                       &
             nzc1,                            &
             'aoo2e_4c',                      &
             nint(atom(center_total)%charge), &
             101                              &
                             )

#ifdef DEBUG_SOC
!> debug print
if(print_lvl > 5)then
call print_x2cmat(                     &
Gmat_frag,                             &
nrows_ext,                             &
ncols_ext,                             &
nzc1,                                  &
aoo_cb_pq_to_uq(1,0),                  &
'2e-mean-field corrections from file', &
6                                      &
)
end if
#endif


allocate(rbuf(naosh_ls))
allocate(cbuf(naosh_ls))
allocate(ibuf(nrows_ext))
allocate(jbuf(ncols_ext))

rbuf = 0
cbuf = 0
ibuf = 0
jbuf = 0

call put_aoblock(Gmat,           &
naosh_ls,                        &
naosh_ls,                        &
center_total,                    &
aoo_nont(center_id),             &
aoo_nucdeg(center_total),        &
nzc1,                            &
Gmat_frag,                       &
nrows_ext,                       &
nrows_ext,                       &
-1,                              &
ncols_ext,                       &
ncols_ext,                       &
-1,                              &
rbuf,                            &
cbuf,                            &
ibuf,                            &
jbuf)

#ifdef DEBUG_SOC
!> debug print
if(print_lvl > 5)then
call print_x2cmat(                     &
Gmat,                                  &
naosh_ls,                              &
naosh_ls,                              &
nr_quat_ext,                           &
aoo_cb_pq_to_uq(1,0),                  &
'2e-mean-field corrections accumulated',    &
6                                      &
)
end if
#endif
deallocate(jbuf)
deallocate(ibuf)
deallocate(cbuf)
deallocate(rbuf)
deallocate(Gmat_frag)
!write(6,*) 'done correction for center ',center_total,aoo_nont(center_id),center_id,aoo_nucdeg(center_total)
!call flshfo(6)

end if

!     update counter
center_total = center_total + aoo_nont(center_id)

end subroutine put_atomic_2e_mean_field_correction
!----------------------------------------------------------------------

subroutine put_atomic_density_correction(&
  Gmat,                        &
  naosh_ls,                    &
  nrows_ext,                   &
  ncols_ext,                   &
  nr_quat_ext,                 &
  center_id,                   &
  center_total,                &
  xcmode,                      &
  print_lvl                    &
 )
!**********************************************************************
real(8), intent(inout)             :: Gmat(*)
integer, intent(in)                :: naosh_ls
integer, intent(in)                :: nrows_ext
integer, intent(in)                :: ncols_ext
integer, intent(in)                :: nr_quat_ext
integer, intent(in)                :: center_id
integer, intent(inout)             :: center_total
integer, intent(in)                :: xcmode
integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
real(8), allocatable               :: Gmat_frag(:)
real(8), allocatable               :: rbuf(:)
real(8), allocatable               :: cbuf(:)
real(8), allocatable               :: ibuf(:)
real(8), allocatable               :: jbuf(:)
INTEGER, parameter                 :: nzc1 = 4
character(len=8)                   :: prefix
!----------------------------------------------------------------------

allocate(Gmat_frag(nrows_ext*ncols_ext*nzc1)); Gmat_frag = 0
!       !> get atomic density from file (of atom(i) ... == atomic fragment)
                prefix = 'aooDM_2c'       
if(xcmode == 4) prefix = 'aooDM_4c'
call read_xamfi_general_matrix(                  &
                Gmat_frag,                       &
                nrows_ext,                       &
                ncols_ext,                       &
                nzc1,                            &
                prefix(1:8),                     &
                nint(atom(center_total)%charge), &
                101                              &
                              )

#ifdef DEBUG_SOC
!> debug print
if(print_lvl > 5)then
call print_x2cmat(                     &
Gmat_frag,                             &
nrows_ext,                             &
ncols_ext,                             &
nzc1,                                  &
aoo_cb_pq_to_uq(1,0),                  &
'atomic density from file',            &
6                                      &
)
end if
#endif


allocate(rbuf(naosh_ls))
allocate(cbuf(naosh_ls))
allocate(ibuf(nrows_ext))
allocate(jbuf(ncols_ext))

rbuf = 0
cbuf = 0
ibuf = 0
jbuf = 0

call put_aoblock(Gmat,           &
naosh_ls,                        &
naosh_ls,                        &
center_total,                    &
aoo_nont(center_id),             &
aoo_nucdeg(center_total),        &
nzc1,                            &
Gmat_frag,                       &
nrows_ext,                       &
nrows_ext,                       &
-1,                              &
ncols_ext,                       &
ncols_ext,                       &
-1,                              &
rbuf,                            &
cbuf,                            &
ibuf,                            &
jbuf)

#ifdef DEBUG_SOC
!> debug print
if(print_lvl > 5)then
call print_x2cmat(                           &
Gmat,                                        &
naosh_ls,                                    &
naosh_ls,                                    &
nr_quat_ext,                                 &
aoo_cb_pq_to_uq(1,0),                        &
'atomic density corrections accumulated',    &
6                                            &
)
end if
#endif
deallocate(jbuf)
deallocate(ibuf)
deallocate(cbuf)
deallocate(rbuf)
deallocate(Gmat_frag)

!     update counter
center_total = center_total + aoo_nont(center_id)

end subroutine put_atomic_density_correction
!----------------------------------------------------------------------

  function collect_delta_xc_energy(                             &
                                   center_id,                   &
                                   center_total,                &
                                   print_lvl                    &
   ) result(delta_xc_e)
!**********************************************************************
integer, intent(in)                :: center_id
integer, intent(inout)             :: center_total
integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
real*8                             :: delta_xc_e
!----------------------------------------------------------------------


delta_xc_e = 0.0d0

if(nint(atom(center_total)%charge) > 1)then

delta_xc_e =  read_xamfi_general_scalar(                                 &
                                        'aooDLXCE',                      &
                                        nint(atom(center_total)%charge), &
                                        101                              &
                                       )
!> debug print
if(print_lvl > 1)then
write(6,'(/a,f12.9)')&
'DELTA XC energy correction from file', delta_xc_e
end if


!write(6,*) 'done correction for center ',center_total,aoo_nont(center_id),center_id,aoo_nucdeg(center_total)
!call flshfo(6)

!> multiply the correction term with the degeneracy factor for this center
delta_xc_e = delta_xc_e * aoo_nucdeg(center_total)

!> debug print
if(print_lvl > 1)then
write(6,'(/a,f12.9)')&
'DELTA XC energy correction multiplied with degeneracy of nuclear center', delta_xc_e
end if

end if !> for nuclear charge of center > 1

!     update counter
center_total = center_total + aoo_nont(center_id)

end function collect_delta_xc_energy
!----------------------------------------------------------------------

#ifdef DEBUG_SOC
#undef DEBUG_SOC
#endif
  end module xamfi_driver
