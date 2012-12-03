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

  implicit none

  public get_atomic_oo_order_so_correction
  public put_atomic_oo_order_so_correction

  private

contains

!----------------------------------------------------------------------
  subroutine get_atomic_oo_order_so_correction(            &
                                               dmat,       &
                                               pct_mat,    &
                                               print_lvl   &
                                               )
!**********************************************************************
     real(8), intent(inout)             :: dmat(nrows,ncols,nr_quat,nr_2e_fock_matrices)
     real(8), intent(in)                :: pct_mat(*)
     integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
     integer                            :: i
     real(8), allocatable               :: fmat(:)
!----------------------------------------------------------------------

      allocate(fmat(nrows*ncols*nr_quat*nr_2e_fock_matrices))
      fmat = 0.0d0

!     get full (scalar+spin-orbit) 2e-fock matrix
      call get_scso_2e_fock(                          &
                            dmat,                     &
                            fmat,                     &
                            nrows,                    &
                            ncols,                    &
                            nr_2e_fock_matrices,      &
                            aoo_intflg,               &
                            aoo_mdirac,               &
                            nr_quat,                  &
                            aoo_cb_pq_to_uq(1,0),     &
                            aoo_ihqmat(1,-1),         &
                            nr_ao_total_aoo,          &
                            print_lvl                 &
                           )


      if(print_lvl > 2)then
        CALL HEADER( &
        'blubb aoosoc twofock 2-electron Fock matrix for you',-1)
        CALL PRQMAT(fmat,nr_ao_total_aoo,nr_ao_total_aoo,nr_ao_total_aoo,nr_ao_total_aoo, &
                    nr_quat,aoo_cb_pq_to_uq(1,0),6)
      end if

!     get 2e-SO-fock matrix
      call make_so_2e_fock(                            &
                           fmat,                       & 
                           nrows/2,                    &
                           ncols/2,                    &
                           nr_quat                     &
                          )

!     normalize (picture-change transform) 2e-SO-fock matrix
      call get_normalized_2e_soc(                      &
                                 fmat,                 &
                                 pct_mat,              &
                                 dmat,                 &
                                 nrows,                &
                                 ncols,                &
                                 nr_2e_fock_matrices,  &
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

!     pure 2e-SO matrix
      call pure_so_2e(                            &
                      fmat,                       & 
                      nrows,                      &
                      ncols,                      &
                      nr_quat                     &
                      )

!     save spin-orbit parts on file (of atom(1) ... == atomic fragment)
      call dump_normalized_2e_soc(                            &
                                  fmat,                       & 
                                  nrows,                      &
                                  ncols,                      &
                                  nr_quat,                    &
                                  file_name_aoo_contributions,&
                                  nint(atom(1)%charge),       &
                                  aoo_fsoc                    &
                                 )
     
!     release memory for atom type
      if(allocated(atom)) deallocate(atom)                             

  end subroutine get_atomic_oo_order_so_correction
!----------------------------------------------------------------------

  subroutine put_atomic_oo_order_so_correction(                             &
                                               Gmat,                        &
                                               naosh_ls,                    &
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
     integer, intent(in)                :: naosh_ls
     integer, intent(in)                :: nrows_ext
     integer, intent(in)                :: ncols_ext
     integer, intent(in)                :: nr_quat_ext
     integer, intent(in)                :: center_id
     integer, intent(inout)             :: center_total
     integer, intent(in)                :: max_center
     integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
     real(8), allocatable               :: Gmat_frag(:)
     real(8), allocatable               :: work(:)
     integer                            :: lwork
     integer                            :: lfree
     integer                            :: kfree
!----------------------------------------------------------------------

      allocate(Gmat_frag(nrows_ext*ncols_ext*4))

!     get spin-orbit parts from file (of atom(i) ... == atomic fragment)
      call read_normalized_2e_soc(                                 &
                                  Gmat_frag,                       & 
                                  nrows_ext,                       &
                                  ncols_ext,                       &
                                  nr_quat_ext,                     &
                                  file_name_aoo_contributions,     &
                                  nint(atom(center_total)%charge), &
                                  aoo_fsoc                         &
                                 )
!     distribute  atomic data in molecular block 
      lwork = 2*naosh_ls + nrows_ext + ncols_ext + 24
      allocate(work(lwork))
      kfree = 1
      lfree = lwork
      call put_aoblock(Gmat,naosh_ls,center_total,aoo_nont(center_id),aoo_nucdeg(center_total),4, &
                       Gmat_frag,nrows_ext,ncols_ext,-1,nrows_ext,ncols_ext,-1,work,kfree,lfree   &
                      )

      deallocate(work)
      deallocate(Gmat_frag)

!     update counter
      center_total = center_total + aoo_nont(center_id)

!     release memory when last type has been added
      if(center_total == max_center)then
        deallocate(atom)
        deallocate(aoo_nont)
        deallocate(aoo_nucdeg)
      end if

  end subroutine put_atomic_oo_order_so_correction
!----------------------------------------------------------------------

  end module atomic_oo_order_so_correction
