!dirac_copyright_start
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
!dirac_copyright_end

module dirx2c

! stefan: - this module contains all necessary functionality 
!           to perform the transformation from a four-component 
!           to a two-component framework 
!
!           written by sknecht, jan 2010 - july 2010

! space for some notes on the pure 1-step approach: 
!
!
!
! ---------------------------------------------------------------------------------------
!                            "roadmap" of the X2C code
! ---------------------------------------------------------------------------------------
!
!   x2c_driver: top-level routine looping over molecular/atomic fragments
!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       1. do fragment = 1 , all_fragments
!
!           --> x2c_tra: for each molecular/atomic fragments:
!    
!               --> a. set up a defining matrix 4c-h1 in orthonormal MO basis
! 
!               --> b. get decoupling matrix "R" in orthonormal MO basis
!                   +  get picture-change transformation matrix U in orthonormal MO basis   
!
!               --> c. construct picture-change transformation matrix U in DIRAC-sorted SA-AO basis
! 
!               --> d. picture-change transform the "basic" 4c-h1 matrix in DIRAC-sorted SA-AO basis
!                      note: the defining matrix 4c-h1 and the "basic" 4c-h1 matrix 
!                            may or may not be identical.
!
!               --> e. (add AMFI contributions to the bare-nucleus 1e Dirac-Hamiltonian in case 
!                      h1 == bare-nucleus 1e Dirac-Hamiltonian (or extended Huckel guess) .and. AMFI is desired)
!
!       end do
!
!       2. (x2c_combine_aba_x2c: approximate a "full" picture-change transformation matrix U from all fragments)
!
!       3. pctrafo_all_op_driver: picture-change transform 
!
!           a. list of 1e property operator 
!
!           b. bare-nucleus corrected full 1e Dirac-Hamiltonian
!
!
!
!   successful X2C procedure: renew common blocks
!   +++++++++++++++++++++++++++++++++++++++++++++
!
! ---------------------------------------------------------------------------------------
!                            end of "roadmap" of the X2C code
! ---------------------------------------------------------------------------------------
!
  use memory_allocator
! x2c modules
  use x2cmod_cfg
  use x2c_cb_interface
  use x2c_def_h1_4c
  use x2c_get_h1inf2c_mo
  use x2c_decoupling_mat
  use x2c_h1inf2c_mo2ao
  use x2c_amfi_interface
  use x2c_pct_ao
  use x2c_fio
  use x2c_utils
  use atomic_x2c
  use x2c_h1mat_base
  use x2c_pct_mo_coefficients
  use x2c_import_exportC1

  implicit none

  public x2c_tra
  public x2c_driver

  private

  save

  real(8), parameter :: x2c_d0       =  0.0D0
  real(8), parameter :: x2c_dummy(2) =  0.0D0

contains

!**********************************************************************
  subroutine x2c_driver(work,                   &
                        lwork,                  &
                        atom_by_atom_x2c,       &
                        x2c_nuccent,            &
                        return_h1_inf_2c)
!**********************************************************************
!
!    purpose: driver routine of the X2C algorithm
!
!             for a detailed overview of the code structure the programmer 
!             may be refered to the "figure" at the beginning of this file.
!             A in-depth theoretical description of the molecular X2C (one-step) 
!             procedure can be found in reference 
!             "M. Ilias and T. Saue, J. Chem. Phys., 126, 064102 (2007)". 
!
!             the driver routine splits into three major parts:
!                 1. loop over (molecular fragments)/ "all centers" for decoupling
!                 2. approximate a combined picture-change transformation matrix U_comb
!                 3. picture-change transform all relevant operator
!
!**********************************************************************
     integer, intent(in)    :: lwork
     real(8), intent(inout) :: work(lwork)
     integer, intent(in)    :: x2c_nuccent        ! total number of nuclei/molecular fragments
     logical, intent(in)    :: atom_by_atom_x2c   ! atom-by-atom or molecular x2c run
     logical, intent(in)    :: return_h1_inf_2c   ! return transformed h1_2c^++
!----------------------------------------------------------------------
     real(8), pointer       :: wrk(:)
     real(8), allocatable   :: xmat(:)
     real(8), allocatable   :: scr1mat(:)
     real(8), allocatable   :: scr2mat(:)
     real(8), allocatable   :: eigvl(:)
     integer                :: prop_idummy(2)
     character (len=16)     :: prop_cdummy(2)
     integer                :: lwrk
     integer                :: is_allocated_dyn
     integer                :: i,j
     integer                :: is_center_x2c
     integer                :: dim1_fsym_Sao
     integer                :: dim2_fsym_Sao
     integer                :: final_ham_lvl_x2c  ! level of final hamiltonian stored on global file X2CMAT
     logical                :: op_2c_ll
     logical                :: lin_aomo_mat_req
     character (len=12)     :: flabel             ! label for matrices on file
!**********************************************************************

!      ------------------------------------------------------------------------------
!      step 0: initialize and print header
!      ------------------------------------------------------------------------------

!      atom-by-atom X2C run
       if( atom_by_atom_x2c )then
         call quit(' *** error in x2c_driver: atom-by-atom X2C code is not activated yet. ***')
       end if

!      print header
       call hello_dirx2c()

!      initialize interface to f77 work array
       call set_wrk_x2c(work,lwork)

!      initialize variable(s)
       final_ham_lvl_x2c = 0
       is_allocated_dyn  = 0
 
!      open file X2CMAT
       open(x2c_funit,file='X2CMAT',status='replace',form='unformatted',    &
            access='sequential',action='readwrite',position='rewind')
       call x2c_write(flabel,x2c_dummy,-1,x2c_funit)

!      ------------------------------------------------------------------------------
!      step 1: loop over (molecular fragments)/ "all centers"
!      ------------------------------------------------------------------------------

       do is_center_x2c = 1, x2c_nuccent


!        open scratch file x2cscr (data for previous fragment is erased if atom-by-atom X2C)
         open(x2c_fscru,file='x2cscr',status='replace',form='unformatted',  &
              access='sequential',action="readwrite",position='rewind')
         call x2c_write(flabel,x2c_dummy,-1,x2c_fscru)

!        open scratch file X2Camfi_scr (data for previous fragment is erased if atom-by-atom X2C)
         open(x2c_famfi,file='X2Camfi_scr',status='replace',form='unformatted',  &
              access='sequential',action="readwrite",position='rewind')
         call x2c_write(flabel,x2c_dummy,-1,x2c_famfi)

!        ------------------------------------------------------------------------------
!        step 1.0: initialize dimensions in Dirac-sorted SA-AO basis and symmetry settings 
!        ------------------------------------------------------------------------------
         
         call init_x2c_cb(.true.,               &
                          atom_by_atom_x2c,     &
                          is_center_x2c,        &
                          prop_idummy,          &
                          prop_idummy,          &
                          prop_cdummy)

!        set logical value for matrices (Lowdin and others) to be adapted for linear symmetry
         lin_aomo_mat_req = ((x2c_is_defining_h1mat == 2).and.x2c_linsym)

!        allocate necessary matrices and arrays - xmat for step 1.1 + 1.2
         call alloc(xmat   ,nr_ao_total_x2c**2 * nr_quat)
         xmat    = x2c_d0
         is_allocated_dyn  = nr_ao_total_x2c**2 * nr_quat

!        ------------------------------------------------------------------------------
!        step 1.1: compute 4c-Sao_LL overlap matrix for fragment "is_center_x2c"
!        ------------------------------------------------------------------------------

!        open temporary scratch file for 4c-Sao_LL overlap matrix
         open(x2c_Sao_scrF,file='X2C_Sao_LL_scr',status='replace',form='unformatted',    &
              access='sequential',action="readwrite",position='rewind')
         call x2c_write(flabel,x2c_dummy,-1,x2c_Sao_scrF)

         dim1_fsym_Sao = nr_ao_all(1)
         dim2_fsym_Sao = 0
         if(nr_fsym == 2) dim2_fsym_Sao = nr_ao_all(nr_fsym)

         call get_Sao_LL_mat(xmat,dim1_fsym_Sao,nr_fsym,dim2_fsym_Sao)
!                                                       ^ (compute+write if >= 0)
!        ------------------------------------------------------------------------------
!        step 1.2: get LOWDIN matrix for fragment "is_center_x2c".
!                  1.2.a. get the 4c-Lowdin matrix
!                  1.2.b. get the 2c-Lowdin matrix (only needed if the defining Hamiltonian is the 4c-Fock operator)
!        ------------------------------------------------------------------------------

!        step 1.2.a. get the 4c-Lowdin matrix
!        ------------------------------------
         lwrk = len_wrk_f77_x2c - is_allocated_dyn
         call sel_wrk_x2c(wrk,lwrk)
         call lowdin(xmat,wrk,lwrk,x2c_prt)
         call desel_wrk_x2c(wrk,lwrk)

!        allocate remaining required matrices and arrays: TODO: check if dimensions can actually be reduced...
         j                               = nr_ao_total_x2c**2 * nr_quat
         if(x2c_is_defining_h1mat < 1) j = n2bastq_dim_x2c
         call alloc(scr2mat,nr_ao_total_x2c**2 * nr_quat)
         call alloc(scr1mat,j)
         call alloc(eigvl  ,nr_ao_total_x2c)

         scr1mat = x2c_d0
         scr2mat = x2c_d0
         eigvl   = x2c_d0

!        set maximum allowed length of wrk array allocations
         is_allocated_dyn = is_allocated_dyn + nr_ao_total_x2c + nr_ao_total_x2c**2 * nr_quat + j
         lwrk             = len_wrk_f77_x2c  - is_allocated_dyn

!        step 1.2.b. get the 2c-Lowdin matrix
!        ------------------------------------
         if(x2c_is_defining_h1mat == 2)then
!          construct the two-component Lowdin matrix
!          initialize lowdin dimensions
           call init_x2c_cb(.false.,           &
                            atom_by_atom_x2c,  &
                            is_center_x2c,     &
                            prop_idummy,       &
                            prop_idummy,       &
                            prop_cdummy)

!          open temporary scratch file for 4c-Lowdin matrix
           open(x2c_lowdin_scrF,file='X2C_lowdin_scr',status='replace',form='unformatted',    &
                access='sequential',action="readwrite",position='rewind')
           call x2c_write(flabel,x2c_dummy,-1,x2c_lowdin_scrF)

!          store the 4c-Lowdin matrix temporarily
           call x2c_write('x2c-4clowdin',xmat,nr_ao_total_x2c*nr_ao_total_x2c,x2c_lowdin_scrF)

           call x2c_lowd2c(xmat,               &
                           scr1mat,            &
                           scr2mat,            &
                           eigvl,              &
                           dim_pshell,         &
                           dim_eshell,         &
                           nr_ao_l,            &
                           nr_ao_s,            &
                           nr_ao_all,          &
                           nr_mo_lw_l,         &
                           nr_mo_lw_s,         &
                           nr_mo_lw_all,       &
                           nr_ao_total_x2c,    &
                           nr_ao_large_x2c,    &
                           nr_fsym,            &
                           nr_quat,            &
                           lwrk,               &
                           is_center_x2c,      &
                           x2c_funit,          &
                           lin_aomo_mat_req)         

!          recover the 4c-Lowdin matrix
           call x2c_read('x2c-4clowdin',xmat,nr_ao_total_x2c*nr_ao_total_x2c,x2c_lowdin_scrF)
           close(x2c_lowdin_scrF,status="delete")

         end if

!        ------------------------------------------------------------------------------
!        step 1.3: construct modified dirac equation in ON basis --> retrieve AOtoMO transformation matrices.
!        ------------------------------------------------------------------------------
!
!              in : 4c-LOWDIN matrix -> xmat
!              out: AOtoMO transformation matrices -> xmat, scr1mat
!                                                   + (scr2mat, if x2c_is_defining_h1mat == 2 .and. x2c_linsym == .true.) 
!        ------------------------------------------------------------------------------


         call sel_wrk_x2c(wrk,lwrk)

         call modham(scr1mat,                  &
                     xmat,                     &
                     scr2mat,                  &
                     eigvl,                    &
                     .true.,                   &
                     lin_aomo_mat_req,         &
                     x2c_is_defining_h1mat,    &
                     wrk,                      &
                     lwrk)

!        nullify pointer to work array and release matrices
         call desel_wrk_x2c(wrk,lwrk)
         call dealloc(eigvl)
         if(x2c_is_defining_h1mat < 1)then
           call dealloc(scr1mat)
!          reallocate with full size
           call alloc(scr1mat,nr_ao_total_x2c**2 * nr_quat)
           scr1mat = x2c_d0
         end if
!
!        recalculate the new maximum allowed length of wrk array allocations to be used in the subroutines.
         is_allocated_dyn = is_allocated_dyn - nr_ao_total_x2c - j + nr_ao_total_x2c**2 * nr_quat
         lwrk             = len_wrk_f77_x2c  - is_allocated_dyn

!        ------------------------------------------------------------------------------
!        step 1.4: transfer dimensions on common blocks concerning the orthonormal MO basis -
!                  these dimensions have been evaluated in MODHAM. 
!        ------------------------------------------------------------------------------
         call init_x2c_cb(.true.,               &
                          atom_by_atom_x2c,     &
                          is_center_x2c,        &
                          prop_idummy,          &
                          prop_idummy,          &
                          prop_cdummy)

!        ------------------------------------------------------------------------------
!        step 1.5: obtain 
!                  a. decoupling matrix "R" in orthonormal MO basis
!                  b. picture-change transformation matrix "U" in Dirac-sorted SA-AO basis
!                  c. h1_2c_inf (bare-nucleus 1e Dirac-Hamiltonian) in Dirac-sorted SA-AO basis
!        ------------------------------------------------------------------------------

         call x2c_tra(xmat,                   &
                      scr1mat,                &
                      scr2mat,                &
                      nr_ao_large_x2c,        &
                      nr_ao_total_x2c,        &
                      nr_ao_all,              &
                      nr_ao_l,                &
                      dim_oshell,             &
                      dim_eshell,             &
                      dim_pshell,             &
                      fullomat_dim,           &
                      fulleomat_dim,          &
                      ioff_aomat_x,           &
                      nr_fsym,                &
                      nr_quat,                &
                      x2c_cb_pq_to_uq,        &
                      is_center_x2c,          &
                      x2c_funit,              &
                      x2c_fscru,              &
                      x2c_famfi,              &
                      x2c_is_defining_h1mat,  &
                      final_ham_lvl_x2c,      &
                      x2c_add_amfi,           &
                      x2c_bs_to_fs,           &
                      lwrk,                   &
                      x2c_prt,                &
                      atom_by_atom_x2c,       &
                      x2c_linsym,             &
                      return_h1_inf_2c,       & 
                      x2c_do_spherical)

!        free memory
         call dealloc(scr1mat)
         call dealloc(scr2mat)
         call dealloc(xmat)

         close(x2c_Sao_scrF,status="delete")

!        --------------------------------------------------------------------------------------------------
!        step 1.6: add AMFI contributions in case
!                  h1 == bare-nucleus 1e Dirac-Hamiltonian (or extended Huckel guess) .and. AMFI is desired
!        --------------------------------------------------------------------------------------------------
         if(x2c_add_amfi >0)then
           call x2c_amfi_driver(final_ham_lvl_x2c,                 &
                                x2c_add_amfi,                      &
                                nr_ao_total_x2c,                   &
                                nr_ao_large_x2c,                   &
                                nr_ao_all,                         &
                                nr_ao_l,                           &
                                nr_fsym,                           &
                                nr_quat,                           &
                                x2c_bs_irrep_mat,                  &
                                x2c_cb_pq_to_uq,                   &
                                x2c_cb_uq_to_pq,                   &
                                x2c_iqmult_trip_q,                 &
                                x2c_bs_to_fs,                      &
                                x2c_pointer_quat,                  &
                                x2c_pointer_quat_op,               &
                                x2c_max_quant_num,                 &
                                num_nuclei,                        &
                                x2c_tot_charge_for_mfsum,          &
                                x2c_amfi_order,                    &
                                len_wrk_f77_x2c,                   &
                                x2c_do_spherical,                  &
                                x2c_famfi,                         &
                                x2c_funit,                         &
                                x2c_cspeed,                        &
                                x2c_prt)
         end if

         close(x2c_famfi,status="delete")

!      ------
       end do ! end of loop over (molecular fragments)/ "all centers"
!      ------

!      -------------------------------------------------------
!      step 2: atom-by-atom X2C run: construct full pct-matrix
!      -------------------------------------------------------
       if(atom_by_atom_x2c.and.(x2c_nuccent > 1))then

         call x2c_atomic_x2c_driver()

       end if

!      ------------------------------------------------------------------------------------------
!      step 4: picture-change transform all (remaining) operator (in DIRAC-sorted SA-AO basis):
!              a. if(appropriate) bare-nucleus corrected 1e Dirac-Hamiltonian
!              b. if(appropriate) 1-electron/2-electron Fock matrices
!              c. list of 1e-property operators
!
!              in : i.  on file X2CMAT   (fh: x2c_funit):
!                      (combined) picture-change transformation matrix "U"
!                   ii. on file AOPROPER (fh: 1int_4c_file): 
!                       all property integrals in Hermit-sorted AO basis
!
!              out: all operator (listed above) in picture-change transformed form on file
!      -----------------------------------------------------------------------------------
       call pctrafo_all_op_driver(nr_ao_total_x2c,               &
                                  nr_ao_large_x2c,               &
                                  nr_ao_all,                     &
                                  nr_ao_l,                       &
                                  nr_fsym,                       &
                                  nr_quat,                       &
                                  x2c_cb_pq_to_uq,               &
                                  ioff_aomat_x,                  &
                                  x2c_funit,                     &
                                  fh_1int_4c,                    &
                                  x2c_bs_to_fs,                  &
                                  all_prp_op,                    &
                                  len_wrk_f77_x2c,               &
                                  x2c_do_spherical,              &
                                  x2c_prep_bnccorr,              &
                                  x2c_1e_system,                 &
                                  x2c_is_defining_h1mat,         &
                                  nr_2e_fock_matrices,           &
                                  file_name_1e_fock_matrix,      &
                                  file_name_2e_fock_matrices,    &
                                  x2c_prt)
!      ---------------------------------------------------------------------------------------------------------------
!      step 3: picture-change transform the 4c-MO coefficients (defining
!      operator matrix is the 4c-Dirac Fock operator)
!      ---------------------------------------------------------------------------------------------------------------
       if(x2c_is_defining_h1mat == 2)then
         call x2c_pctrafo_mo_coefficients_driver('DFCOEF',        &
                                                 nr_cmo_q,        &
                                                 nr_ao_total_x2c, &
                                                 nr_ao_large_x2c, &
                                                 nr_ao_all,       &
                                                 nr_ao_l,         &
                                                 dim_oshell,      &
                                                 dim_eshell,      &
                                                 dim_pshell,      &
                                                 nr_fsym,         &
                                                 nr_quat,         &
                                                 x2c_cb_pq_to_uq, &
                                                 x2c_funit,       &
                                                 x2c_prt)
       end if
!      --------------------------------------------------------------------------------------------
!      step 5: export picture-change transformation matrix to AO basis with nosym (C1 double group)
!      --------------------------------------------------------------------------------------------
       if(x2cmod_h1dirty)then
         call x2c_export_pctmat2C1()
       end if

!      ----------------------------------------------------------------------------------------
!      step 6: renew dimensions on common blocks, i.e. "remove" references to positronic shells
!      ----------------------------------------------------------------------------------------
       call renew_x2c_cb_orb_shell_dim(final_ham_lvl_x2c)

       close(x2c_funit,status="keep")

!      print goodbye 
       call goodbye_dirx2c()

  end subroutine

!**********************************************************************
  subroutine x2c_tra(xmat,                 &
                     scr1mat,              &
                     scr2mat,              &
                     nr_ao_L,              &
                     naosh_ls,             &
                     naosh_all,            &
                     naosh_L,              &
                     norb_dim,             &
                     nesh_dim,             &
                     npsh_dim,             &
                     h1mat_dim_4c,         &
                     momat_orbesh_dim_4c,  &
                     ioff_aomx,            &
                     nfsym,                &
                     nz,                   &
                     ipq_off,              &
                     is_fragment,          &
                     x2c_file_glb,         &
                     x2c_file_scr,         &
                     x2c_file_amf,         &
                     defining_h1mat,       &
                     is_final_ham_lvl,     &
                     add_amfi_contrib,     &
                     op_bs_to_fs,          &
                     is_len_wrk,           &
                     print_lvl,            &
                     fragment_x2c,         &
                     linear_sym,           &
                     return_h1_inf_2c,     &
                     spherical_on)
!**********************************************************************
!
!    purpose: X2C transformation driver routine
!
!             requires: 4c-Lowdin matrix; 
!                       value of "defining_h1mat" determines the defining 4c-h1 matrix
!                       wrt which the decoupling step is performed.
!
!                       defining_h1mat == 0: bare-nucleus one-electron Dirac-Hamiltonian
!                                      == 1: (extended) Huckel start guess
!                                      == 2: converged 4c-Fock-Dirac operator (not implemented yet)
!
!             delivers: a. decoupling matrix "R" in orthonormal MO basis
!                       b. picture-change transformation matrix "U" in orthonormal MO basis and Dirac-sorted SA-AO basis
!                       c. h1_2c_inf in Dirac-sorted SA-AO basis
!
!    for more information see also the equation (9), (10), (23) and (24) in the
!    reference: M. Ilias and T. Saue, J. Chem. Phys., 126, 064102 (2007)
!
!**********************************************************************
     real(8), intent(inout) :: xmat(*)
     real(8), intent(inout) :: scr1mat(*)
     real(8), intent(inout) :: scr2mat(*)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: nr_ao_L
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: naosh_L(nfsym)
     integer, intent(in)    :: norb_dim(nfsym)
     integer, intent(in)    :: nesh_dim(nfsym)
     integer, intent(in)    :: npsh_dim(nfsym)
     integer, intent(in)    :: h1mat_dim_4c
     integer, intent(in)    :: momat_orbesh_dim_4c
     integer, intent(in)    :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)    :: ipq_off(4,0:7)
     integer, intent(in)    :: is_fragment
     integer, intent(in)    :: x2c_file_glb
     integer, intent(in)    :: x2c_file_scr
     integer, intent(in)    :: x2c_file_amf
     integer, intent(in)    :: defining_h1mat
     integer, intent(inout) :: is_final_ham_lvl
     integer, intent(in)    :: add_amfi_contrib
     integer, intent(in)    :: op_bs_to_fs(0:7,1:2)
     integer, intent(in)    :: is_len_wrk
     integer, intent(in)    :: print_lvl
     logical, intent(in)    :: fragment_x2c
     logical, intent(in)    :: linear_sym
     logical, intent(in)    :: return_h1_inf_2c
     logical, intent(in)    :: spherical_on
!----------------------------------------------------------------------
     real(8), allocatable   :: eigvl(:)
     integer                :: i,j
     integer                :: lwrk
     integer                :: is_fragment_pctmat
     integer                :: basic_h1mat        ! defining parameter for the basic 4c-h1 which will be transformed to 2c-h1
     character (len=12)     :: flabel             ! label for matrices on file
!**********************************************************************
!

!      initialize maximum allowed length of wrk array allocations
       lwrk = is_len_wrk

!      step 1: save RKB AOtoMO transformation matrices to file x2cscr
!
!              in : AOtoMO transformation matrices
!              out: ["SL resorted"] matrix -> scr2mat
!              on file x2cscr: 
!                 a. AOtoMO transformation matrix ["SL resorted"] under the label 4cAOMOr
!                 b. AOtoMO transformation matrix under the label 4cAOMOo
!                 c. MOtoMO non-linear -> linear symmetry transformation matrix under the label 4cMOMOl
!      ------------------------------------------------------------------------------
       call store_aomo_trafo_matrices(xmat,                &
                                      scr1mat,             &
                                      scr2mat,             &
                                      norb_dim,            &
                                      naosh_all,           &
                                      nfsym,               &
                                      nz,                  &
                                      defining_h1mat,      &
                                      is_fragment,         &
                                      x2c_file_scr,        &
                                      linear_sym)

!      step 2: store the hamiltonian for which the 4c-->2c transformation 
!              should be done on file (== basic "h1_4c" in orthonormal MO basis) h1_4cON
!               in: RKB AOtoMO transformation matrix ["SL resorted"] -> scr2mat
!              out: - (nothing in matrices)
!                   on file x2cscr: basic 4c-h1 in orthonormal basis under label "h1_4cON"
!      ------------------------------------------------------------------------------
!      the basic "h1_4c" will always be the bare-nucleus h1 hamiltonian except 
!      for the "4c-Dirac-Fock operator as basis" option. sknecht+tsaue: march 2011
                               basic_h1mat = 0
       if(defining_h1mat == 2) basic_h1mat = 2
       call x2c_get_h1_mat_base(xmat,                  &
                                scr1mat,               &
                                scr2mat,               &
                                naosh_ls,              &
                                naosh_all,             &
                                norb_dim,              &
                                ioff_aomx,             &
                                nfsym,                 &
                                nz,                    &
                                ipq_off,               &
                                basic_h1mat,           &
                                lwrk,                  &
                                is_fragment,           &
                                x2c_file_scr,          &
                                linear_sym,            &
                                print_lvl)
       
!      step 3: get defining matrix (== "h1-defining" in orthonormal MO basis)
!              in : ["SL resorted"] matrix -> scr2mat
!              out: defining 4c-h1 in orthonormal basis -> xmat
!      ------------------------------------------------------------------------------
       call x2c_get_h1_defining_mat(xmat,                  &
                                    scr1mat,               &
                                    scr2mat,               &
                                    naosh_ls,              &
                                    naosh_all,             &
                                    norb_dim,              &
                                    ioff_aomx,             &
                                    nfsym,                 &
                                    nz,                    &
                                    ipq_off,               &
                                    defining_h1mat,        &
                                    lwrk,                  &
                                    is_fragment,           &
                                    x2c_file_scr,          &
                                    linear_sym,            &
                                    print_lvl)
 
!      allocate scratch matrix with appropriate dimension
       call alloc(eigvl  ,max(naosh_ls,momat_orbesh_dim_4c,h1mat_dim_4c))
       eigvl   = x2c_d0

!      set new maximum allowed length of wrk array allocations
       lwrk = lwrk - max(naosh_ls,momat_orbesh_dim_4c,h1mat_dim_4c)

!      step 4: get the transformed basic operator "h1_4c" --> "h1_infinite_2c" in orthonormal basis
!
!              in : defining  4c-h1 in orthonormal basis -> xmat
!                   on file x2cscr in orthonormal basis: the basic operator "h1_4c" (label == h1_4cON)
!              out: the basic 2c-h1 in orthonormal basis -> xmat
!                   on file x2cscr in orthonormal basis: pct-matrix U, "R"-matrix for nuclei is_fragment
!      ------------------------------------------------------------------------------
       call get_h1_2c_onmo_bas(xmat,                      &
                               scr1mat,                   &
                               scr2mat,                   &
                               eigvl,                     &
                               norb_dim,                  &
                               nesh_dim,                  &
                               npsh_dim,                  &
                               nfsym,                     &
                               nz,                        &
                               ipq_off,                   &
                               is_fragment,               &
                               x2c_file_scr,              &
                               defining_h1mat,            &
                               print_lvl)

!      step 5: get the transformed operator "h1_infinite_2c" in SA-AO basis
!
!              in : 2c-h1 in orthonormal MO basis -> xmat
!              out: 2c-h1 in SA-AO basis -> xmat + *
!                   
!              * == on file X2CMAT     : under the label h12cAOn (written in h1_2c^++ format)
!                         or
!              * == on file X2Camfi_scr: under the label h12cAOa (written in 4c h1_2c^++ format) 
!                                        in preparation for AMFI contributions.
!                                              
!      ------------------------------------------------------------------------------
       call get_h1_inf_2c_saao_bas(xmat,                    &
                                   eigvl,                   &
                                   scr1mat,                 &
                                   scr2mat,                 &
                                   norb_dim,                &
                                   nesh_dim,                &
                                   naosh_L,                 &
                                   naosh_all,               &
                                   nr_ao_L,                 &
                                   naosh_ls,                &
                                   ioff_aomx,               &
                                   nfsym,                   &
                                   nz,                      &
                                   ipq_off,                 &
                                   is_fragment,             &
                                   is_final_ham_lvl,        &
                                   add_amfi_contrib,        &
                                   lwrk,                    &
                                   x2c_file_scr,            &
                                   x2c_file_glb,            &
                                   x2c_file_amf,            &
                                   op_bs_to_fs,             &
                                   spherical_on,            &
                                   fragment_x2c,            &
                                   basic_h1mat,             &
                                   print_lvl)

!      step 6: construct the picture-change transformation matrix in SA-AO basis
!
!              in : -
!              out: -                         
!                   on file X2CMAT: picture-change transformation matrix U in SA-AO basis 
!                                   under the label pctmtAO
!      ------------------------------------------------------------------------------
       if(.not. fragment_x2c) then 
         is_fragment_pctmat = 0
       else
         is_fragment_pctmat = is_fragment
       end if

       call construct_pctmat_saao_bas(xmat,                &
                                      scr1mat,             &
                                      scr2mat,             &
                                      eigvl,               &
                                      norb_dim,            &
                                      nesh_dim,            &
                                      naosh_all,           &
                                      naosh_L,             &
                                      nfsym,               &
                                      nz,                  &
                                      ipq_off,             &
                                      is_fragment_pctmat,  &
                                      is_fragment,         &
                                      x2c_file_scr,        &
                                      x2c_file_glb,        &
                                      print_lvl)

!      free scratch memory
       call dealloc(eigvl)

       close(x2c_file_scr,status="delete")

  end subroutine

!**********************************************************************

end module
