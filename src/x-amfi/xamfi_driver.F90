!dirac_copyright_start
!      Copyright (c) 2018 by the authors of DIRAC.
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
#undef DEBUG_SOC
#endif

  implicit none

  public get_xamfi_correction
  public put_xamfi_correction

  private

contains

!----------------------------------------------------------------------
  subroutine get_xamfi_correction(            &
                                  dmat,       &
                                  cmo,        &
                                  pct_mat,    &
                                  calc_2esoc, &
                                  lwork_ext,  &
                                  print_lvl   &
                                  )
!**********************************************************************
     real(8), intent(inout)             :: dmat(nrows,ncols,nr_quat,nr_2e_fock_matrices)
     real(8), intent(inout)             :: cmo(*)
     real(8), intent(in)                :: pct_mat(*)
     logical, intent(in)                :: calc_2esoc
     integer, intent(in)                :: lwork_ext
     integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
     integer                            :: i
     real(8), allocatable               :: fmat(:)
!----------------------------------------------------------------------

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
                            ioff_tmot             = ioff_tmot,            &
                            ioff_tmt              = ioff_tmt,             &
                            ioff_aomat_x          = ioff_aomat_x,         &
                            myAOMOMAT             = 'AOMOMATaoo',         &
                            lwork_ext             = lwork_ext,            &
                            print_lvl             = print_lvl             &
                           )

#ifdef DEBUG_SOC
      !> debug print
      if(print_lvl > 5)then
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

!#define DEBUG_SOC
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
                                  aoo_fsoc                    &
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
     real(8), optional, intent(inout)   :: Gmat(*)
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
!----------------------------------------------------------------------


      if(nint(atom(center_total)%charge) > 1)then

        allocate(Gmat_frag(nrows_ext*ncols_ext*4)); Gmat_frag = 0
!       !> get spin-orbit parts from file (of atom(i) ... == atomic fragment)
        call read_normalized_2e_soc(                                 &
                                    Gmat_frag,                       &
                                    nrows_ext,                       &
                                    ncols_ext,                       &
                                    nr_quat_ext,                     &
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
                          nr_quat_ext,                           &
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
                         4,                               &
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
#ifdef DEBUG_SOC  
#undef DEBUG_SOC
#endif
  end module xamfi_driver