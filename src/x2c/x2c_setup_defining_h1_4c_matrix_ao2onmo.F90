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
! module containing the selection of the defining matrix "h1" from which 
! the decoupling matrix shall be extracted.
!
! written by sknecht   june 2010 - original implementation
!
!            sknecht august 2012 - added InteRest interface for 1-e Hamiltonian 
!                                - added test support for X2C-spinfree (before the decoupling)
!                                  simply activate "x2c_spinfree_" in the next line
!#define x2c_spinfree_before_4c_to_2c
!                                  default for X2C-spinfree is "after the decoupling"
!
module x2c_def_h1_4c

  use x2c_import_exportC1
  use x2c_fio
  use x2c_cb_interface, only:         &
      reset_x2c_cb_onefck
  use x2c_utils,        only:         &
      print_x2cmat,                   &
!#ifdef x2c_spinfree_before_4c_to_2c
      make_spinfree_h1_onmo,          &
      get_boson_irrep_info,           &
!#endif
      read_4cfock_operator_x2c,       &
      get_saomo_x2c,                  &
      generic_interface_ao2mo_mo2ao
  use x2cmod_cfg
  use x2c_2e_soc_interface, only:     &
      x2c_add_2e_diagonal_atomic,     &
      x2c_add_2e_from_atomic_dens

  implicit none

  public x2c_get_h1_defining_mat
  public get_h1_4c_saao_basis
  public get_h1_4c_on_basis

  private

  real*8, parameter :: val_d1      =  1.0d0
  real*8, parameter :: val_d0      =  0.0d0
  real*8, parameter :: val_dm1     = -1.0d0

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
                                     nzt,                  &
                                     ipq_off_in,           &
                                     defining_h1mat,       &
                                     is_len_wrk,           &
                                     x2c_file_scr,         &
                                     linear_sym,           &
                                     mdirac,               &
                                     spinfree,             &
                                     isDFT,                &
                                     print_lvl)
!**********************************************************************
!
!    purpose: driver routine to obtain a 'defining' h1 matrix in orthonormal MO basis
!             currently implemented: 
!               - defining_h1mat == 0: case a. bare-nucleus 1-el Dirac-Hamiltonian
!               - defining_h1mat == 1: case b. bare-nucleus 1-el Dirac-Hamiltonian + diagonal (atomic) 2e-mean-field contributions
!               - defining_h1mat == 2: case c. 4c-Fock-Dirac operator matrix)
!               - defining_h1mat == 3: case d. free-particle matrix
!               - defining_h1mat == 4: case e. bare-nucleus 1-el Dirac-Hamiltonian + first-order (molecular) 2e-mean-field contributions
!
!----------------------------------------------------------------------
!
     real*8 , intent(inout) :: h1_onmo_bas(*)
     real*8 , intent(inout) :: h1_saao_bas(*)
     real*8 , intent(inout) :: tra_mat_4cAOMO(*)
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: norb_dim(nfsym)
     integer, intent(in)    :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: nzt
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(in)    :: defining_h1mat
     integer, intent(in)    :: is_len_wrk
     integer, intent(in)    :: x2c_file_scr
     integer, intent(in)    :: print_lvl
     logical, intent(in)    :: linear_sym
     logical, intent(in)    :: mdirac
     logical, intent(in)    :: spinfree
     logical, intent(in)    :: isDFT
!----------------------------------------------------------------------
     integer                :: ioff_osh1
     integer                :: ioff_osh2
     integer                :: ioff_ao_o_sh1   
     integer                :: ioff_ao_ao_sh1
     integer                :: i       
     integer                :: norb1_f       
     integer                :: nbas1_f
     integer                :: ndim_ao
     integer                :: mz
     integer, allocatable   :: boson_info(:)
!**********************************************************************

!     step 1: get h1 ('defining matrix') in SA-AO basis according to value for defining_h1mat
!
!             in : defining_h1mat as integer value: 
!                  0 == case a.; bare-nucleus 4c-h1 
!                  1 == case b.; bare-nucleus 4c-h1 + diagonal (atomic) 2e-mean-field contributions
!                  2 == case c.; 4c-Dirac-Fock matrix
!                  3 == case d.; free-particle matrix
!                  4 == case e.; bare nucleus Dirac-Hamiltonian + first-order (molecular) 2e-mean-field contributions
!             out: 4c h1 matrix in SA-AO basis --> h1_saao_bas
!     ------------------------------------------------------------------------------
      call get_h1_4c_saao_basis(h1_saao_bas,            &
                                h1_onmo_bas,            &
                                tra_mat_4cAOMO,         &
                                naosh_ls,               &
                                naosh_all,              &
                                norb_dim,               &
                                nfsym,                  &
                                nz,                     &
                                nzt,                    &
                                ipq_off_in,             &
                                defining_h1mat,         &
                                linear_sym,             &
                                mdirac,                 &
                                isDFT,                  &
                                is_len_wrk,             &
                                x2c_file_scr,           &
                                print_lvl)

!     step 2: transform h1 to orthonormal MO basis if defining h1 matrix is 
!             case 0. 1-el. bare nucleus Dirac-Hamiltonian 
!             case 1. 1-el. bare nucleus Dirac-Hamiltonian + diagonal (atomic) 2e-mean-field contributions (amfX2C)
!             case 3. free-particle matrix
!             case 4. 1-el. bare nucleus Dirac-Hamiltonian + first-order (molecular) 2e-mean-field contributions (eamfX2C)
!
!             in : 4c h1 matrix in SA-AO basis --> h1_saao_bas; transformation matrix --> tra_mat_4cAOMO
!             out: h1 matrix in 4c-ON MO basis --> h1_onmo_bas
!     ------------------------------------------------------------------------------
      mz = nz
      if (nfsym == 2) then
        allocate(boson_info(norb_dim(1)+norb_dim(2)))
      else
        allocate(boson_info(norb_dim(1)))
      endif
      if(spinfree)then
        call get_boson_irrep_info(boson_info,  &
                                  nfsym,       &
                                  norb_dim,    &
                                  print_lvl)
        mz = 1
      end if

!     initialize matrix offsets
      ioff_osh1      = 1
      ioff_osh2      = 1
      ioff_ao_o_sh1  = 1
      ioff_ao_ao_sh1 = 1

      do i = 1, nfsym

!       set dimensions and special offsets for matrices accessed in the subroutine
        nbas1_f = naosh_all(i)
        norb1_f = norb_dim(i)
        ndim_ao = naosh_all(i)
!       the h1 matrices (1el-dirac or free-particle matrix) have a different storage mode 
!       in memory in saao basis than the 4c-fock operator
        if(defining_h1mat /= 2)then 
          ioff_ao_ao_sh1 = ioff_aomx(i,i) + 1
          ndim_ao        = naosh_ls
        end if

!       call print_x2cmat( tra_mat_4cAOMO(ioff_ao_o_sh1),nbas1_f,norb1_f,nz,ipq_off_in,'blubb sf aomo',6)
!       call print_x2cmat( h1_saao_bas(ioff_ao_ao_sh1),nbas1_f,ndim_ao,nz,ipq_off_in,'blubb sf h1aoba',6)

        if( norb_dim(i) > 0 )then

          call get_h1_4c_on_basis(h1_onmo_bas(ioff_osh1),            &
                                  h1_saao_bas(ioff_ao_ao_sh1),       &
                                  tra_mat_4cAOMO(ioff_ao_o_sh1),     &
                                  ndim_ao,                           &
                                  nbas1_f,                           &
                                  norb1_f,                           &
                                  nz,                                &
                                  mz,                                &
                                  nzt,                               &
                                  ipq_off_in,                        &
                                  spinfree,                          &
                                  boson_info(ioff_osh2),             &
                                  print_lvl)
        end if
!       call print_x2cmat(h1_onmo_bas(ioff_osh1),norb1_f,norb1_f,nz,ipq_off_in,'blubb sf h1onmo',6)

!       update offsets for matrices
        ioff_ao_o_sh1  = ioff_ao_o_sh1  + nbas1_f * norb1_f * nzt
        ioff_ao_ao_sh1 = ioff_ao_ao_sh1 + nbas1_f * nbas1_f * nz 
        ioff_osh1      = ioff_osh1      + norb1_f**2 * nz
        ioff_osh2      = ioff_osh2      + norb1_f

      end do

      deallocate(boson_info)

  end subroutine x2c_get_h1_defining_mat

!**********************************************************************
  subroutine get_h1_4c_saao_basis(h1_saao_bas,            &
                                  scrmat1,                &
                                  scrmat2,                &
                                  naosh_ls,               &
                                  naosh_all,              &
                                  norb_dim,               &
                                  nfsym,                  &
                                  nz,                     &
                                  nzt,                    &
                                  ipq_off_in,             &
                                  defining_h1mat,         &
                                  linear_sym,             &
                                  mdirac,                 &
                                  isDFT,                  &
                                  is_len_wrk,             &
                                  x2c_file_scr,           &
                                  print_lvl)
!**********************************************************************
!
!    purpose: get 'defining' matrix h1 in symmetry-adapted (SA)-AO basis
!             currently implemented: 
!               - defining_h1mat == 0: bare-nucleus 1-el Dirac-Hamiltonian
!               - defining_h1mat == 1: bare-nucleus 1-el Dirac-Hamiltonian + diagonal (atomic) 2e-mean-field contributions
!               - defining_h1mat == 2: 4c-Fock-Dirac operator matrix
!               - defining_h1mat == 3: free-particle matrix
!               - defining_h1mat == 4: bare-nucleus 1-el Dirac-Hamiltonian + first-order (molecular) 2e-mean-field contributions
!
!             note: scrmat1 and scrmat2 are in use exclusively for
!             defining_h1mat == 2: 4c-fock operator as h1 (defining and/or 'basic' h1)
!
!----------------------------------------------------------------------
     real*8 , intent(inout) :: h1_saao_bas(*)
     real*8 , intent(inout) :: scrmat1(*)
     real*8 , intent(inout) :: scrmat2(*)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: norb_dim(nfsym)
     integer, intent(in)    :: nz
     integer, intent(in)    :: nzt
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(in)    :: defining_h1mat
     logical, intent(in)    :: linear_sym
     logical, intent(in)    :: mdirac
     logical, intent(in)    :: isDFT
     integer, intent(in)    :: is_len_wrk
     integer, intent(in)    :: x2c_file_scr
     integer, intent(in)    :: print_lvl
!----------------------------------------------------------------------
     real*8 , allocatable   :: wrk(:)
     integer                :: lwrk
     integer                :: integral_types
!**********************************************************************

!      give full length for WORK array allocations
       lwrk = is_len_wrk

       select case(defining_h1mat)
       case(0)
!        bare-nucleus h1 Dirac-Hamiltonian
#ifdef PRG_DIRAC
         allocate(wrk(lwrk))
         call reset_x2c_cb_onefck(.false.)
         call onefck(h1_saao_bas,print_lvl,wrk,lwrk)
         call reset_x2c_cb_onefck(.true.)
         deallocate(wrk)
#endif

       case(1)
!        bare-nucleus h1 Dirac-Hamiltonian...
#ifdef PRG_DIRAC
        allocate(wrk(lwrk))
        call reset_x2c_cb_onefck(.false.)
        call onefck(h1_saao_bas,print_lvl,wrk,lwrk)
        call reset_x2c_cb_onefck(.true.)
        deallocate(wrk)
#endif

#ifdef MOD_XAMFI
!       ... + diagonal (atomic) 2e-mean-field contributions
        call x2c_add_2e_diagonal_atomic(h1_saao_bas,                       &
                                        naosh_ls,                          &
                                        nfsym,                             &
                                        nz,                                &
                                        x2c_bs_irrep_mat,                  &
                                        x2c_cb_pq_to_uq,                   &
                                        x2c_cb_uq_to_pq,                   &
                                        x2c_iqmult_trip_q,                 &
                                        x2c_bs_to_fs,                      &
                                        x2c_pointer_quat,                  &
                                        num_nuclei,                        &
                                        type_nuclei,                       &
                                        print_lvl)
#endif
       case(2)
!        4c-Fock-Dirac Hamiltonian...

!        out: 4c-fock operator in SA-AO basis -> h1_saao_bas
         call get_4cdf_saao_basis(h1_saao_bas,        &
                                  scrmat1,            &
                                  scrmat2,            &
                                  naosh_ls,           &
                                  naosh_all,          &
                                  norb_dim,           &
                                  nfsym,              &
                                  nz,                 &
                                  nzt,                &
                                  ipq_off_in,         &
                                  linear_sym,         &
                                  mdirac,             &
                                  x2c_file_scr,       &
                                  print_lvl)

       case(3)

#ifdef PRG_DIRAC
!        Free-particle matrix h_fp
         allocate(wrk(lwrk))
         call reset_x2c_cb_onefck(.false.)
         call freemt(h1_saao_bas,print_lvl,wrk,lwrk)
         call reset_x2c_cb_onefck(.true.)
         deallocate(wrk)
#endif
        case(4)
  !        bare-nucleus h1 Dirac-Hamiltonian...
#ifdef PRG_DIRAC
          allocate(wrk(lwrk))
          call reset_x2c_cb_onefck(.false.)
          call onefck(h1_saao_bas,print_lvl,wrk,lwrk)
          call reset_x2c_cb_onefck(.true.)
          deallocate(wrk)
#endif

          integral_types = x2c_intgen_4c
#ifdef MOD_XAMFI
  !       ... + F^{4c,2e}[D^4c_(+)atoms]
          call x2c_add_2e_from_atomic_dens(h1_saao_bas,                       &
                                           naosh_ls,                          &
                                           nz,                                &
                                           x2c_bs_irrep_mat,                  &
                                           x2c_cb_pq_to_uq,                   &
                                           x2c_cb_uq_to_pq,                   &
                                           x2c_iqmult_trip_q,                 &
                                           x2c_bs_to_fs,                      &
                                           x2c_pointer_quat,                  &
                                           num_nuclei,                        &
                                           type_nuclei,                       &
                                           4,                                 & ! xcmode
                                           integral_types,                    & ! integral flag
                                           lwrk,                              & ! lwrk
                                           isDFT,                             & ! add XC potential (DFT)
                                           print_lvl                          &
                                          )
#endif

       case default
         print *, ' get_h1_4c_saao_basis, unkown selection code defining_h1_mat=',&
         defining_h1mat
         call quit('get_h1_4c_saao_basis, unkown selection code defining_h1_mat')
       end select

  end subroutine get_h1_4c_saao_basis

!**********************************************************************
  subroutine get_4cdf_saao_basis(h1_saao_bas,            &
                                 scrmat1,                &
                                 scrmat2,                &
                                 naosh_ls,               &
                                 naosh_all,              &
                                 norb_dim,               &
                                 nfsym,                  &
                                 nz,                     &
                                 nzt,                    &
                                 ipq_off_in,             &
                                 linear_sym,             &
                                 mdirac,                 &
                                 x2c_file_scr,           &
                                 print_lvl)
!**********************************************************************
!
!    purpose: get the 4c-fock operator in symmetry-adapted (SA)-AO basis
!
!             if it is already available on file x2cscr 
!             (for example when it was prepared as basic_h1mat) just read it; 
!             otherwise obtain it from DFFOCK in orthonormal MO basis (LS sorted) 
!             and perform the following mo2ao transformation sequence:
!                                                  if(linear_sym)
!             4c-fock_onmo (LS sorted) --> (4c-fock_onmo_nonlin LS sorted) --> 4c-fock_saao
!
!----------------------------------------------------------------------
     use x2cmod_cfg
     real(8), intent(inout)       :: h1_saao_bas(*)
     real(8), intent(inout)       :: scrmat1(*)
     real(8), intent(inout)       :: scrmat2(*)
     integer, intent(in)          :: nfsym
     integer, intent(in)          :: naosh_ls
     integer, intent(in)          :: naosh_all(nfsym)
     integer, intent(in)          :: norb_dim(nfsym)
     integer, intent(in)          :: nz
     integer, intent(in)          :: nzt
     integer, intent(in)          :: ipq_off_in(4,0:7)
     logical, intent(in)          :: linear_sym
     logical, intent(in)          :: mdirac
     integer, intent(in)          :: x2c_file_scr
     integer, intent(in)          :: print_lvl
!----------------------------------------------------------------------
     integer                      :: ioff_osh1
     integer                      :: ioff_ao_o_sh1
     integer                      :: ioff_ao_ao_sh1
     integer                      :: i
     integer                      :: norb1_f
     integer                      :: nbas1_f
     character (len=6)            :: op_mat_on_file
     character (len=12)           :: flabel
     character (len=28)           :: debug_string
     logical                      :: fndlab12
!**********************************************************************

      rewind(x2c_file_scr)
      if(fndlab12('h14csAO   11',x2c_file_scr))then

        call dzero(h1_saao_bas,naosh_ls*naosh_ls*nz)
        ioff_ao_ao_sh1 = 1
        do i = 1, nfsym

!         read the 4c-fock operator in SA-AO basis from file (label == h1_4cAO)
          nbas1_f = naosh_all(i)
          write(flabel,'(a7,i4,i1)') 'h14csAO',1,i

          if(nbas1_f > 0) call x2c_read(flabel,h1_saao_bas(ioff_ao_ao_sh1),nbas1_f**2 * nz,x2c_file_scr)
          ioff_ao_ao_sh1 = ioff_ao_ao_sh1 + nbas1_f**2 * nz
        end do

      else
!       step 1: read 4c-fock_onmo (LS sorted) from file op_mat_on_file
!       --------------------------------------------------------------
!       note: the following if clause saves one dcopy call and might therefore 
!             be worthwhile. the 4c-fock_onmo matrix is expected to reside in 
!             scrmat1 for "step 4" below.

        op_mat_on_file = 'DFFOCK'

        if(linear_sym)then
          call read_4cfock_operator_x2c(h1_saao_bas,       &
                                        op_mat_on_file,    &
                                        print_lvl)
        else
          call read_4cfock_operator_x2c(scrmat1,           &
                                        op_mat_on_file,    &
                                        print_lvl)
        end if

        if(x2c_fock_saao_basis) goto 100 ! debug option - can be activated in the input

!       step 2: 4c-fock_onmo (LS sorted) --> 4c-fock_onmo_nonlin LS sorted
!       ------------------------------------------------------------------
        if(linear_sym)then
          ioff_osh1 = 1
          do i = 1, nfsym
            norb1_f = norb_dim(i)
            if(norb1_f > 0)then

              write(flabel,'(a7,i4,i1)') '4cMOMOl',1,i
              call x2c_read(flabel,scrmat2(ioff_osh1),norb1_f**2 * nz,x2c_file_scr)

!             debug print
              if(print_lvl > 2)then
                write(debug_string,'(a27,i1)') 'x2c - 4c-fock_onmol- fsym = ',i
                call print_x2cmat(h1_saao_bas(ioff_osh1),norb1_f,norb1_f,nz,ipq_off_in,debug_string,6)
                write(debug_string,'(a27,i1)') 'x2c - on_onMOlin2nl- fsym = ',i
                call print_x2cmat(scrmat2(ioff_osh1),norb1_f,norb1_f,nz,ipq_off_in,debug_string,6)
              end if

!             scrmat1 (alias 4c-fock_onmo_nonlin) = scrmat2 * h1_saao_bas (alias 4c-fock_onmo) * scrmat2^+
              call generic_interface_ao2mo_mo2ao(scrmat1(ioff_osh1),        &
                                                 h1_saao_bas(ioff_osh1),    &
                                                 scrmat2(ioff_osh1),        &
                                                 norb1_f,                   &
                                                 norb1_f,                   &
                                                 norb1_f,                   &
                                                 nz,                        &
                                                 nz,                        &
                                                 nz,                        &
                                                 ipq_off_in,                &
                                                 'MOAO',                    &
                                                 'S',                       &
                                                 print_lvl)
!             debug print
              if(print_lvl > 2)then
                write(debug_string,'(a27,i1)') 'x2c - 4c-fock_onNOL- fsym = ',i
                call print_x2cmat(scrmat1(ioff_osh1),norb1_f,norb1_f,nz,ipq_off_in,debug_string,6)
              end if
            end if
!           update offsets for matrices
            ioff_osh1     = ioff_osh1     + norb1_f**2 * nz
          end do
        end if

!       step 3: get half-transformed overlap matrix Sao
!       -----------------------------------------------
        ioff_ao_o_sh1 = 1
        do i = 1, nfsym
          nbas1_f = naosh_all(i)
          norb1_f = norb_dim(i)
          write(flabel,'(a7,i4,i1)') '4cAOMOo',1,i
          if(nbas1_f * norb1_f > 0) call x2c_read(flabel,h1_saao_bas(ioff_ao_o_sh1), &
                                          nbas1_f * norb1_f * nzt, x2c_file_scr)
          ioff_ao_o_sh1 = ioff_ao_o_sh1 + nbas1_f * norb1_f * nzt
        end do

        call get_saomo_x2c(                &
                           scrmat2,        &
                           h1_saao_bas,    &
                           naosh_ls,       &
                           naosh_all,      &
                           norb_dim,       &
                           ioff_aomat_x,   &
                           nfsym,          &
                           nz,             &
                           nzt,            &
                           ipq_off_in,     &
                           mdirac,         &
                           print_lvl)

!       step 4: 4c-fock_onmo_nonlin LS sorted --> 4c-fock_saao
!       ------------------------------------------------------
 100    continue
        ioff_ao_ao_sh1 = 1
        ioff_osh1      = 1
        ioff_ao_o_sh1  = 1
        do i = 1, nfsym
          nbas1_f = naosh_all(i)
          norb1_f = norb_dim(i)
          if(norb1_f > 0)then

!           save on file x2cscr for use in post-X2Cmmf steps like (e)amf
            write(flabel,'(a7,i4,i1)') '4cSaomo',1,i
            call x2c_write(flabel,scrmat2(ioff_ao_o_sh1),nbas1_f*norb1_f * nzt,x2c_file_scr)

!           h1_saao_bas = scrmat2 * scrmat1 * scrmat2^+
            if(.not.x2c_fock_saao_basis)then
            call generic_interface_ao2mo_mo2ao(h1_saao_bas(ioff_ao_ao_sh1),    &
                                               scrmat1(ioff_osh1),             &
                                               scrmat2(ioff_ao_o_sh1),         &
                                               nbas1_f,                        &
                                               norb1_f,                        &
                                               nbas1_f,                        &
                                               nz,                             &
                                               nz,                             &
                                               nzt,                            &
                                               ipq_off_in,                     &
                                               'MOAO',                         &
                                               'S',                            &
                                               print_lvl)
            end if


          end if

!         save on file x2cscr for re-use as defining h1
          write(flabel,'(a7,i4,i1)') 'h14csAO',1,i
          call x2c_write(flabel,h1_saao_bas(ioff_ao_ao_sh1),nbas1_f**2 * nz,x2c_file_scr)
!         debug print
          if(print_lvl > 2)then
            write(debug_string,'(a27,i1)') 'x2c - 4c-fock_AO   - fsym = ',i
            call print_x2cmat(h1_saao_bas(ioff_ao_ao_sh1),nbas1_f,nbas1_f,nz,ipq_off_in,debug_string,6)
          end if

          ioff_ao_ao_sh1 = ioff_ao_ao_sh1 + nbas1_f * nbas1_f * nz
          ioff_osh1      = ioff_osh1      + norb1_f * norb1_f * nz
          ioff_ao_o_sh1  = ioff_ao_o_sh1  + nbas1_f * norb1_f * nzt
        end do

        call dzero(scrmat1,naosh_ls**2*nz)
        call dzero(scrmat2,naosh_ls**2*nz)

!
!       step 5: read in the "SL resorted" transformation matrix (the calling programs expect it to reside in scrmat2
!       ------------------------------------------------------------------------------------------------------------
        ioff_ao_o_sh1 = 1
        do i = 1, nfsym
          nbas1_f = naosh_all(i)
          norb1_f = norb_dim(i)
          write(flabel,'(a7,i4,i1)') '4cAOMOr',1,i
          if(nbas1_f * norb1_f > 0) call x2c_read(flabel,scrmat2(ioff_ao_o_sh1),nbas1_f*norb1_f*nzt,x2c_file_scr)
          ioff_ao_o_sh1 = ioff_ao_o_sh1 + nbas1_f * norb1_f * nzt
        end do

      end if ! if fndlab12('h1_4cAO   11') --> h1_4cAO already exists on file x2cscr

  end subroutine get_4cdf_saao_basis

!**********************************************************************
  subroutine get_h1_4c_on_basis(h1_onmo_bas,            &
                                h1_saao_bas,            &
                                tra_mat_4cAOMO,         &
                                naosh_ls,               &
                                nbas1_f,                &
                                norb1_f,                &
                                nz1,                    &
                                nz2,                    &
                                nzt,                    &
                                ipq_off_in,             &
                                spinfree,               &
                                boson_info,             &
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
     integer, intent(in)    :: boson_info(*)
     integer, intent(in)    :: nbas1_f
     integer, intent(in)    :: norb1_f
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: nz1
     integer, intent(in)    :: nz2
     integer, intent(in)    :: nzt
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(in)    :: print_lvl
     logical, intent(in)    :: spinfree
!----------------------------------------------------------------------

!      h1_4cON = V+ * h1_4cao * V
!      h1_onmo_bas = tra_mat_4cAOMO+ * h1_saao_bas * tra_mat_4cAOMO
       call generic_interface_ao2mo_mo2ao(h1_saao_bas,           &
                                          h1_onmo_bas,           &
                                          tra_mat_4cAOMO,        &
                                          nbas1_f,               &
                                          norb1_f,               &
                                          naosh_ls,              &
                                          nz1,                   &
                                          nz2,                   &
                                          nzt,                   &
                                          ipq_off_in,            &
                                          'AOMO',                &
                                          'S',                   &
                                          print_lvl) 

!#ifdef x2c_spinfree_before_4c_to_2c
!       test the spinfree-before scheme by Miro and HJAAJ
        if(spinfree)then
          call make_spinfree_h1_onmo(                            &
                                     h1_onmo_bas,                & 
                                     norb1_f,                    & 
                                     norb1_f,                    & 
                                     norb1_f,                    & 
                                     norb1_f,                    & 
                                     boson_info,                 &
                                     boson_info,                 &
                                     nz1,                        &
                                     print_lvl                   &
                                    )
        end if
!#endif

  end subroutine get_h1_4c_on_basis
!**********************************************************************

end module x2c_def_h1_4c
