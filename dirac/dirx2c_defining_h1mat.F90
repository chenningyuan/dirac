!      Copyright (c) 2010 by the authors of DIRAC.
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
!
!
! module containing the selection of the defining matrix "h1" from which 
! the decoupling matrix shall be extracted.
!
! written by sknecht june 2010
!
module x2c_def_h1_4c

  use memory_allocator
  use x2c_fio
  use x2c_cb_interface, only:         &
      sel_wrk_x2c,                    &
      desel_wrk_x2c,                  &
      reset_x2c_cb_onefck

  implicit none

  public x2c_get_h1_defining_mat

  private

  real(8), parameter :: val_d1      =  1.0d0
  real(8), parameter :: val_d0      =  0.0d0
  real(8), parameter :: val_dm1     = -1.0d0

contains

!**********************************************************************
  subroutine x2c_get_h1_defining_mat(h1_onmo_bas,          &
                                     h1_saao_bas,          &
                                     tra_mat_4cAOMO,       &
                                     naosh_ls,             &
                                     naosh_all,            &
                                     norb_dim,             &
                                     ioff_aomx,            &
                                     nfsym,                &
                                     nz,                   &
                                     ipq_off_in,           &
                                     defining_h1mat,       &
                                     is_len_wrk,           &
                                     is_fragment,          &
                                     x2c_file_scr,         &
                                     linear_sym,           &
                                     print_lvl)
!**********************************************************************
!
!    purpose: driver routine to obtain a 'defining' h1 matrix in orthonormal MO basis
!             currently implemented: 
!               - defining_h1mat == 0: case a. bare-nucleus 1-el Dirac-Hamiltonian
!         TODO:(- defining_h1mat == 1: case b. (extended) Huckel-start guess)
!         TODO:(- defining_h1mat == 2: case c. 4c-Fock-Dirac operator matrix)
!
!----------------------------------------------------------------------
     real(8), intent(inout) :: h1_onmo_bas(*)
     real(8), intent(inout) :: h1_saao_bas(*)
     real(8), intent(inout) :: tra_mat_4cAOMO(*)
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: norb_dim(nfsym)
     integer, intent(in)    :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(in)    :: defining_h1mat
     integer, intent(in)    :: is_fragment
     integer, intent(in)    :: is_len_wrk
     integer, intent(in)    :: x2c_file_scr
     integer, intent(in)    :: print_lvl
     logical, intent(in)    :: linear_sym
!----------------------------------------------------------------------
     integer                :: ioff_osh1       
     integer                :: ioff_ao_o_sh1   
     integer                :: ioff_h1_saao_mat
     integer                :: i       
     integer                :: norb1_f       
     integer                :: nbas1_f       
     character (len=12)     :: flabel
!----------------------------------------------------------------------

!      step 1: get h1 ('defining matrix') in SA-AO basis according to value for defining_h1mat
!
!              in : defining_h1mat: integer value: 0 == case a. ; 1 == case b. 2 == case c.
!              out: 4c h1 matrix in SA-AO basis --> h1_saao_bas
!      ------------------------------------------------------------------------------
       call get_h1_4c_saao_basis(h1_saao_bas,            &
                                 defining_h1mat,         &
                                 linear_sym,             &
                                 is_len_wrk,             &
                                 print_lvl)

!      step 2: transform h1 to orthonormal MO basis if defining h1 matrix is 
!              case a. 1-el. bare nucleus Dirac-Hamiltonian 
!              case b. extended Huckel start guess
!
!              in : 4c h1 matrix in SA-AO basis --> h1_saao_bas; transformation matrix --> tra_mat_4cAOMO
!              out: h1 matrix in 4c-ON MO basis --> h1_onmo_bas
!      ------------------------------------------------------------------------------
       if(defining_h1mat < 2)then


!        initialize matrix offsets
         ioff_osh1     = 1
         ioff_ao_o_sh1 = 1

         do i = 1, nfsym

!          set dimensions and special offsets for matrices accessed in the subroutine
           nbas1_f          = naosh_all(i)
           norb1_f          = norb_dim(i)
           ioff_h1_saao_mat = ioff_aomx(i,i)

           if( norb_dim(i) > 0 )then

             call get_h1_4c_on_basis(h1_onmo_bas(ioff_osh1),            &
                                     h1_saao_bas(1+ioff_h1_saao_mat),   &
                                     tra_mat_4cAOMO(ioff_ao_o_sh1),     &
                                     naosh_ls,                          &
                                     nbas1_f,                           &
                                     norb1_f,                           &
                                     nz,                                &
                                     ipq_off_in,                        &
                                     print_lvl)

           end if

!          store defining h1_4cON matrix on file
           write(flabel,'(a7,i4,i1)') 'h1_4cON',is_fragment,i
           call x2c_write(flabel,h1_onmo_bas(ioff_osh1),norb1_f**2 * nz,x2c_file_scr)

!          update offsets for matrices
           ioff_ao_o_sh1 = ioff_ao_o_sh1 + nbas1_f * norb1_f * nz
           ioff_osh1     = ioff_osh1     + norb1_f**2 * nz

         end do
       else

!        step 2: case c.
         call quit('*** error in x2c_get_h1: def. h1 matrix == 4c-FD not implemented yet in the new X2C module. ***')
         
       end if

  end subroutine x2c_get_h1_defining_mat

!**********************************************************************
  subroutine get_h1_4c_saao_basis(h1_saao_bas,            &
                                  defining_h1mat,         &
                                  linear_sym,             &
                                  is_len_wrk,             &
                                  print_lvl)
!**********************************************************************
!
!    purpose: get 'defining' matrix h1 in symmetry-adapted (SA)-AO basis
!             currently implemented: 
!               - defining_h1mat == 0: bare-nucleus 1-el Dirac-Hamiltonian
!         TODO:(- defining_h1mat == 1: (extended) Huckel-start guess)
!         TODO:(- defining_h1mat == 2: 4c-Fock-Dirac operator matrix)
!
!----------------------------------------------------------------------
     real(8), intent(inout) :: h1_saao_bas(*)
     integer, intent(in)    :: defining_h1mat
     logical, intent(in)    :: linear_sym
     integer, intent(in)    :: is_len_wrk
     integer, intent(in)    :: print_lvl
!----------------------------------------------------------------------
     real(8), pointer       :: wrk(:)
     integer                :: lwrk
!**********************************************************************

!      give full length for WORK array allocations
       lwrk = is_len_wrk
       call sel_wrk_x2c(wrk,lwrk)

       if(defining_h1mat == 0)then

         call reset_x2c_cb_onefck(.false.)
         call onefck(h1_saao_bas,print_lvl,wrk,lwrk)
         call reset_x2c_cb_onefck(.true.)

       else if(defining_h1mat == 1)then
!        HUCKEL start guess...
       else if(defining_h1mat == 2)then
!        4c-Fock-Dirac Hamiltonian...
       end if

!      reset pointer to WORK array
       call desel_wrk_x2c(wrk,lwrk)

  end subroutine get_h1_4c_saao_basis

!**********************************************************************
  subroutine get_h1_4c_on_basis(h1_onmo_bas,            &
                                h1_saao_bas,            &
                                tra_mat_4cAOMO,         &
                                naosh_ls,               &
                                nbas1_f,                &
                                norb1_f,                &
                                nz,                     &
                                ipq_off_in,             &
                                print_lvl)
!**********************************************************************
!
!    purpose: transform the defining matrix h1 from 4c-SA-AO basis to 
!             4c-orthonormal (ON) MO-basis according to:
!             h1_4cON = V+ * h1_4cao * V
!
!----------------------------------------------------------------------
     real(8), intent(inout) :: h1_onmo_bas(*)
     real(8), intent(inout) :: h1_saao_bas(*)
     real(8), intent(inout) :: tra_mat_4cAOMO(*)
     integer, intent(in)    :: nbas1_f
     integer, intent(in)    :: norb1_f
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: nz
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(in)    :: print_lvl
!----------------------------------------------------------------------
!----------------------------------------------------------------------

!      h1_4cON = V+ * h1_4cao * V
       call qtrans90('AOMO','S',val_d0,          &
                     nbas1_f,                    &
                     nbas1_f,                    &
                     norb1_f,                    &
                     norb1_f,                    &
                     h1_saao_bas,                &
                     naosh_ls,                   &
                     naosh_ls,                   &
                     nz,                         &
                     ipq_off_in(1,0),            &
                     h1_onmo_bas,                &
                     norb1_f,                    &
                     norb1_f,                    &
                     nz,                         &
                     ipq_off_in(1,0),            &
                     tra_mat_4cAOMO,             &
                     nbas1_f,                    &
                     norb1_f,                    &
                     nz,                         &
                     ipq_off_in(1,0),            &
                     tra_mat_4cAOMO,             &
                     nbas1_f,                    &
                     norb1_f,                    &
                     nz,                         &
                     ipq_off_in(1,0),            &
                     print_lvl)

  end subroutine get_h1_4c_on_basis
!**********************************************************************

end module x2c_def_h1_4c
