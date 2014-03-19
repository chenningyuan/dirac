!dirac_copyright_start
!      Copyright (c) 2013 by the authors of DIRAC.
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
module atomic_oo_order_so_correction

  use atomic_oo_order_so_correction_cfg
  use atomic_oo_order_so_correction_utils
  use module_aoosoc_pfg
  use x2c_utils, only: print_x2cmat

  implicit none

  public get_atomic_oo_order_so_correction
  public put_atomic_oo_order_so_correction

  private

contains

!----------------------------------------------------------------------
  subroutine get_atomic_oo_order_so_correction(            &
                                               dmat,       &
                                               pct_mat,    &
                                               calc_2esoc, &
                                               lwork_ext,  &
                                               print_lvl   &
                                               )
!**********************************************************************
     real(8), intent(inout)             :: dmat(nrows,ncols,nr_quat,nr_2e_fock_matrices)
     real(8), intent(in)                :: pct_mat(*)
     logical, intent(in)                :: calc_2esoc
     integer, intent(in)                :: lwork_ext
     integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
     integer                            :: i
     real(8), allocatable               :: fmat(:)
!----------------------------------------------------------------------
      call hello_aoosoc()

      if(nr_2e_fock_matrices > 1)call quit('aoosoc module does not work yet for atomic AOC-HF; use .FOCC in *SCF')

!     allocate(fmat(nrows*ncols*nr_quat*nr_2e_fock_matrices))
      allocate(fmat(nrows*ncols*nr_quat*1))
      fmat = 0.0d0

!     get full (scalar+spin-orbit) 2e-fock matrix
      if(calc_2esoc)then
        call get_scso_2e_fock(                           &
                              dmat,                      &
                              fmat,                      &
                              nrows,                     &
                              ncols,                     &
                              nr_2e_fock_matrices-1,     &
                              aoo_dfopen,                &
                              aoo_intflg,                &
                              aoo_mdirac,                &
                              nr_quat,                   &
                              aoo_cb_pq_to_uq(1,0),      &
                              aoo_ihqmat(1,-1),          &
                              nr_ao_total_aoo,           &
                              lwork_ext,                 &
                              print_lvl                  &
                             )


!       debug print
        if(print_lvl > 2)then
          call print_x2cmat(                             &
                            fmat,                        &
                            nr_ao_total_aoo,             &
                            nr_ao_total_aoo,             &
                            nr_quat,                     &
                            aoo_cb_pq_to_uq(1,0),        &
                            'aoosoc - 2e-complete',      &
                            6                            &
                           )
        end if

#ifdef OLD_VERSION_AOOSOC
!       get 2e-SO-fock matrix
        call make_so_2e_fock(                            &
                             fmat,                       &
                             nrows/2,                    &
                             ncols/2,                    &
                             nr_quat                     &
                            )
!       debug print
        if(print_lvl > 2)then
          call print_x2cmat(                             &
                            fmat,                        &
                            nr_ao_total_aoo,             &
                            nr_ao_total_aoo,             &
                            nr_quat,                     &
                            aoo_cb_pq_to_uq(1,0),        &
                            'aoosoc - 2e-soc-only',      &
                            6                            &
                           )
        end if
#endif

!       normalize (picture-change transform) 2e-SO-fock matrix
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

      end if

!     save spin-orbit parts on file (of atom(1) ... == atomic fragment)
      call dump_normalized_2e_soc(                            &
                                  fmat,                       &
!only true for mdirac             nrows/2,                    &
!only true for mdirac             ncols/2,                    &
                                  nr_ao_large_aoo,            &
                                  nr_ao_large_aoo,            &
                                  nr_quat,                    &
                                  file_name_aoo_contributions,&
                                  nint(atom(1)%charge),       &
                                  aoo_fsoc                    &
                                 )

      call goodbye_aoosoc()

  end subroutine get_atomic_oo_order_so_correction
!----------------------------------------------------------------------

  subroutine put_atomic_oo_order_so_correction(                             &
                                               Gmat,                        &
                                               naosh_l,                     &
                                               nrows_ext,                   &
                                               ncols_ext,                   &
                                               nr_quat_ext,                 &
                                               center_id,                   &
                                               center_total,                &
                                               max_center,                  &
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
     integer, intent(in)                :: max_center
     integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
     real(8), allocatable               :: Gmat_frag(:)
     real(8), allocatable               :: rbuf(:)
     real(8), allocatable               :: cbuf(:)
     real(8), allocatable               :: ibuf(:)
     real(8), allocatable               :: jbuf(:)
!----------------------------------------------------------------------


      if(nint(atom(center_total)%charge) > 1)then

        allocate(Gmat_frag(nrows_ext*ncols_ext*4))
!       print *, 'charge, nrows_ext,total,id',nint(atom(center_total)%charge),nrows_ext,center_total,center_id
!       get spin-orbit parts from file (of atom(i) ... == atomic fragment)
        call read_normalized_2e_soc(                                 &
                                    Gmat_frag,                       &
                                    nrows_ext,                       &
                                    ncols_ext,                       &
                                    nr_quat_ext,                     &
                                    file_name_aoo_contributions,     &
                                    nint(atom(center_total)%charge), &
                                    aoo_fsoc                         &
                                   )
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

        deallocate(jbuf)
        deallocate(ibuf)
        deallocate(cbuf)
        deallocate(rbuf)
        deallocate(Gmat_frag)

      end if

!     update counter
      center_total = center_total + aoo_nont(center_id)
!     print *, 'center_total,  aoo_nont(center_id)',center_total,  aoo_nont(center_id)

  end subroutine put_atomic_oo_order_so_correction
!----------------------------------------------------------------------

  end module atomic_oo_order_so_correction
