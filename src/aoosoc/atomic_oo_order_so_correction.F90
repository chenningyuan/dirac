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
                            print_lvl                 &
                           )


      if(print_lvl > -1)then
        CALL HEADER( &
        'blubb aoosoc twofock 2-electron Fock matrix for you',-1)
        CALL PRQMAT(fmat,nr_ao_total_aoo,nr_ao_total_aoo,nr_ao_total_aoo,nr_ao_total_aoo, &
                    nr_quat,aoo_cb_pq_to_uq(1,0),6)
      end if

!     get pure 2e-SO-fock matrix
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
                                               nrows_ext,                   &
                                               ncols_ext,                   &
                                               nr_quat_ext,                 &
                                               center_id,                   &
                                               center_total,                &
                                               print_lvl                    &
                                              )    
!**********************************************************************
     real(8), optional, intent(inout)   :: Gmat(*)
     integer, intent(in)                :: nrows_ext
     integer, intent(in)                :: ncols_ext
     integer, intent(in)                :: nr_quat_ext
     integer, intent(in)                :: center_id
     integer, intent(in)                :: center_total
     integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
!----------------------------------------------------------------------

!     get spin-orbit parts from file (of atom(i) ... == atomic fragment)
      call read_normalized_2e_soc(                              &
                                  Gmat,                         & 
                                  nrows_ext,                    &
                                  ncols_ext,                    &
                                  nr_quat_ext,                  &
                                  file_name_aoo_contributions,  &
                                  nint(atom(center_id)%charge), &
                                  aoo_fsoc                      &
                                 )
!     distribute  atomic data in molecular hamiltonian block 
!     use tronds put_atomic_block...
      stop 'define distribution of atomic data in molecular hamiltonian block'

!     release memory for atom type
      if(center_id == center_total)then 
        deallocate(atom)                             
      end if

  end subroutine put_atomic_oo_order_so_correction
!----------------------------------------------------------------------

  end module atomic_oo_order_so_correction
