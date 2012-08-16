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
! module containing the functionality to perform the picture-change
! transformation of a four-component (effective) one-electron operator 
! in DIRAC-sorted SA-AO basis to an infinite-order two-component operator 
! following the X2C formalism. 
!
! written by sknecht july 2010
!
module x2c_pct_ao

  use x2c_fio
  use x2c_cb_interface, only:            &
      init_x2c_cb,                       &
      reset_x2c_cb_onefck
  use x2c_utils, only:                   &
      pick_LL_block_saao_bas,            &
      read_fock_matrices_saao_basis_x2c, &
      write_fock_matrices_saao_basis_x2c,&
      make_two_e_fock_so_exclusive,      &
      print_x2cmat
#ifdef MOD_INTEREST
  use module_interest_interface
#endif

  implicit none

  public pctrafo_all_op_driver

  private

  real(8), parameter :: val_d1      =  1.0d0
  real(8), parameter :: val_d0      =  0.0d0
  real(8), parameter :: val_dm1     = -1.0d0

contains

!**********************************************************************
  subroutine pctrafo_all_op_driver(naosh_ls,                     &
                                   nr_ao_L,                      &
                                   naosh_all,                    &
                                   naosh_L,                      &
                                   nfsym,                        &
                                   nz,                           &
                                   ipq_off_in,                   &
                                   ioff_aomx,                    &
                                   x2c_file_glb,                 &
                                   int1_4c_file,                 &
                                   op_bs_to_fs,                  &
                                   num_prp_op,                   &
                                   is_len_wrk,                   &
                                   spherical_on,                 &
                                   mdirac,                       &
                                   prep_bncorr,                  &
                                   is_1e_system,                 &
                                   defining_h1mat,               &
                                   nr_2e_fock_matrices,          &
                                   file_name_1e_fock_matrix,     &
                                   file_name_2e_fock_matrices,   &
                                   print_lvl)   
!**********************************************************************
!
!    purpose: driver routine for the picture-change transformation of 
!             all 1e property operator f_4c listed in the pre-defined 
!             database and (upon request) the bare-nucleus corrected
!             1e Dirac Hamiltonian.
!
!----------------------------------------------------------------------
     integer, intent(in)             :: naosh_ls
     integer, intent(in)             :: nr_ao_L
     integer, intent(in)             :: naosh_all(nfsym)
     integer, intent(in)             :: naosh_L(nfsym)
     integer, intent(in)             :: nfsym
     integer, intent(in)             :: nz
     integer, intent(in)             :: ipq_off_in(4,0:7)
     integer, intent(in)             :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)             :: x2c_file_glb
     integer, intent(in)             :: int1_4c_file
     integer, intent(in)             :: op_bs_to_fs(0:7,1:2)
     integer, intent(in)             :: num_prp_op
     integer, intent(in)             :: is_len_wrk
     integer, intent(in)             :: defining_h1mat
     integer, intent(in)             :: nr_2e_fock_matrices
     integer, intent(in)             :: print_lvl
     logical, intent(in)             :: spherical_on
     logical, intent(in)             :: mdirac
     logical, intent(in)             :: prep_bncorr
     logical, intent(in)             :: is_1e_system
     character (len=6), intent(in)   :: file_name_1e_fock_matrix
     character (len=6), intent(in)   :: file_name_2e_fock_matrices
!----------------------------------------------------------------------
     logical                         :: op_2c_ll
     logical                         :: prop_op_nopct
     integer, allocatable            :: prop_op_dgsym(:)
     integer, allocatable            :: prop_op_trsym(:)
     integer                         :: is_fragment_tot
     character (len=16), allocatable :: prop_op_name(:)
!**********************************************************************

!      initialize the fragment value - after the combination of all fragment 
!      pct-matrices the new default value is 0.
       is_fragment_tot = 0

!      step 1: picture-change transform the bare-nucleus corrected 1e Dirac-Hamiltonian
!      --------------------------------------------------------------------------------
       if(prep_bncorr)then
         call pctrafo_h1bnccorr_4cto2c_saao_bas(naosh_ls,            &
                                                nr_ao_L,             &
                                                naosh_all,           &
                                                naosh_L,             &
                                                nfsym,               &
                                                nz,                  &
                                                ipq_off_in,          &
                                                ioff_aomx,           &
                                                is_fragment_tot,     &
                                                x2c_file_glb,        &
                                                op_bs_to_fs,         &
                                                is_len_wrk,          &
                                                mdirac,              &
                                                print_lvl)   
       end if

!      step 2: picture-change transform the 1e-/2e-Fock matrices
!      ---------------------------------------------------------
       if(defining_h1mat == 2 .and. .not. is_1e_system)then
!#ifdef not_needed_at_all
         call pctrafo_fockmatrices_4cto2c_saao_bas(naosh_ls,                   &
                                                   nr_ao_L,                    &
                                                   naosh_all,                  &
                                                   naosh_L,                    &
                                                   nfsym,                      &
                                                   nz,                         &
                                                   ipq_off_in,                 &
                                                   ioff_aomx,                  &
                                                   is_fragment_tot,            &
                                                   x2c_file_glb,               &
                                                   op_bs_to_fs,                &
                                                   nr_2e_fock_matrices,        &
                                                   file_name_1e_fock_matrix,   &
                                                   file_name_2e_fock_matrices, &
                                                   print_lvl)
!#else
!        the picture-change transformation of the 1e/2e-Fock matrices in AO basis
!        is superfluous - for completeness we delete the files here.
!        stefan - june 2011
         open(99,file=file_name_2e_fock_matrices,status='old',form='unformatted',      &
              access='sequential',action="readwrite",position='rewind')
         close(99,status='delete')
         open(99,file=file_name_1e_fock_matrix,status='old',form='unformatted',        &
              access='sequential',action="readwrite",position='rewind')
         close(99,status='delete')
!#endif
       end if

!      step 3: picture-change transform the list of 1e property operator
!      -------------------------------------------------------------------
       allocate(prop_op_name(num_prp_op*16))
       allocate(prop_op_dgsym(num_prp_op))
       allocate(prop_op_trsym(num_prp_op))

!      initialize property operator attributes
       call init_x2c_cb(.false.,           &
                        prop_op_nopct,     &
                        prop_op_dgsym,     &
                        prop_op_trsym,     &
                        prop_op_name)

!      default for all 1e property operator: write in "op_2c^++" format on file
       op_2c_ll = .true.
       call pctrafo_op_4cto2c_saao_bas(prop_op_dgsym,       &
                                       prop_op_trsym,       &
                                       prop_op_name,        &
                                       naosh_ls,            &
                                       nr_ao_L,             &
                                       naosh_all,           &
                                       naosh_L,             &
                                       nfsym,               &
                                       nz,                  &
                                       ipq_off_in,          &
                                       ioff_aomx,           &
                                       is_fragment_tot,     &
                                       x2c_file_glb,        &
                                       int1_4c_file,        &
                                       op_bs_to_fs,         &
                                       num_prp_op,          &
                                       is_len_wrk,          &
                                       op_2c_ll,            &
                                       prop_op_nopct,       &
                                       spherical_on,        &
                                       print_lvl)

       deallocate(prop_op_dgsym)
       deallocate(prop_op_trsym)
       deallocate(prop_op_name)

  end subroutine pctrafo_all_op_driver

!**********************************************************************
  subroutine pctrafo_op_4cto2c_saao_bas(prop_op_dgsym,       &
                                        prop_op_trsym,       &
                                        prop_op_name,        &
                                        naosh_ls,            &
                                        nr_ao_L,             &
                                        naosh_all,           &
                                        naosh_L,             &
                                        nfsym,               &
                                        nz,                  &
                                        ipq_off_in,          &
                                        ioff_aomx,           &
                                        is_fragment,         &
                                        x2c_file_glb,        &
                                        int1_4c_file,        &
                                        op_bs_to_fs,         &
                                        num_prp_op,          &
                                        is_len_wrk,          &
                                        op_2c_ll,            &
                                        prop_op_nopct,       &
                                        spherical_on,        &
                                        print_lvl)   
!**********************************************************************
!
!    purpose: picture-change transform all 1e property operator f_4c listed 
!             in the pre-defined database.
!             Since the property integrals are sorted by default in the 
!             Hermit format a resorting step to Dirac format precedes 
!             the picture-change transformation (pct).
!
!             equation: f_2c^++ = U+ * f_4c * U
!
!----------------------------------------------------------------------
     integer, intent(in)            :: prop_op_dgsym(num_prp_op)
     integer, intent(in)            :: prop_op_trsym(num_prp_op)
     integer, intent(in)            :: naosh_ls
     integer, intent(in)            :: nr_ao_L
     integer, intent(in)            :: naosh_all(nfsym)
     integer, intent(in)            :: naosh_L(nfsym)
     integer, intent(in)            :: nfsym
     integer, intent(in)            :: nz
     integer, intent(in)            :: ipq_off_in(4,0:7)
     integer, intent(in)            :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)            :: is_fragment
     integer, intent(in)            :: x2c_file_glb
     integer, intent(in)            :: int1_4c_file
     integer, intent(in)            :: op_bs_to_fs(0:7,1:2)
     integer, intent(in)            :: num_prp_op
     integer, intent(in)            :: is_len_wrk
     integer, intent(in)            :: print_lvl
     logical, intent(in)            :: op_2c_ll
     logical, intent(in)            :: prop_op_nopct
     logical, intent(in)            :: spherical_on
     character (len=16), intent(in) :: prop_op_name(num_prp_op)
!----------------------------------------------------------------------
     real(8), allocatable           :: wrk(:)
     real(8), allocatable           :: scr_mat1(:)
     real(8), allocatable           :: scr_mat2(:)
     real(8), allocatable           :: pctmat(:)
     integer                        :: i
     integer                        :: j
     integer                        :: lwrk
     integer                        :: nr_sph
     integer                        :: op_fer_rep
     integer                        :: op_mat_dim
     integer                        :: ioff_pct_mat
     logical                        :: do_pct_prp_op
     logical                        :: aoprop_avail
     logical, allocatable           :: quat_array_filled_tmp(:)
     character (len=12)             :: flabel
     character (len=4)              :: numstr
!**********************************************************************


!      test if any one-particle operator integrals have been created 
       inquire(file="AOPROPER",exist=aoprop_avail)
       if(.not.aoprop_avail) return

!      initialize offset
       ioff_pct_mat = 1

       j = 0
       do i = 1, nfsym
!        print'(a,2i6)', 'naosh_all(i), naosh_L(i)  =',naosh_all(i),naosh_L(i)
         j = j + naosh_all(i)*naosh_L(i)*nz
       end do
!      allocate pct-matrix
       allocate(pctmat(j))
       pctmat = 0

!      step 1: read picture-change transformation matrix U in SA-AO basis from global file
!      ---------------------------------------------------------------------------------

!      print '(a,i0)', 'nfsym       == ',nfsym
!      print '(a,i0)', 'nz          == ',nz
!      print '(a,i0)', 'is_fragment == ',is_fragment

       do i = 1, nfsym
         write(flabel,'(a7,i4,i1)') 'pctmtAO',is_fragment,i
         if(naosh_all(i)*naosh_L(i) > 0) call x2c_read(flabel,pctmat(ioff_pct_mat),naosh_all(i)*naosh_L(i)*nz,     &
                                                       x2c_file_glb)
         ioff_pct_mat  = ioff_pct_mat + naosh_all(i) * naosh_L(i)* nz
       end do

       open(int1_4c_file,file='AOPROPER',status='old',form='unformatted', &
            action="read",position='rewind') 

!      step 2: perform pct in DIRAC-sorted SA-AO basis
!      -----------------------------------------------

!      allocate logical array required for resorting of the 4c-operator from
!      Hermit-sorted AO basis to DIRAC-sorted SA-AO basis.
       allocate(quat_array_filled_tmp(nz))

!      allocate two matrices for pct-procedure of each property operator
!      idea: for optimal reuse of the matrices allocate with full dimensionality
       allocate(scr_mat1(naosh_ls**2 * nz))
       allocate(scr_mat2(naosh_ls**2 * nz))
       scr_mat1        = 0
       scr_mat2        = 0

!      loop over all pre-defined property operator f
       do i = 1, num_prp_op

!        define required operator attributes
         op_fer_rep = prop_op_dgsym(i) - 1
         op_mat_dim = naosh_ls**2
         if(abs(prop_op_trsym(i)) > 0) op_mat_dim = naosh_ls*(naosh_ls+1)/2


!        step 2a: get 4c-operator f_4c in SA-AO basis
!        ------------------------------------------------------------------------------
         call lset(nz,.true.,quat_array_filled_tmp)

         call prpmao(int1_4c_file,                         &
                     i,                                    & 
                     op_mat_dim,                           &
                     scr_mat1,                             &
                     scr_mat2,                             &
                     quat_array_filled_tmp,                &
                     print_lvl)                           

         if(abs(prop_op_trsym(i)) > 0)then 
           scr_mat2 = 0
           call dcopy(op_mat_dim*nz,scr_mat1,1,scr_mat2,1)
         end if

         call sort2sao_no_work(op_mat_dim,                 &
                               scr_mat2,                   &
                               scr_mat1,                   &
                               quat_array_filled_tmp,      &
     &                         op_fer_rep,                 &
                               prop_op_trsym(i),           &
                               print_lvl)

!        step 2b: perform the pct: f_2c^++ == U+ * f_4c * U in DIRAC-sorted SA-AO basis
!        --------------------------------------------------------------------------------
!                 in : f_4c    --> scr_mat1; U --> pctmat
!                 out: f_2c^++ --> scr_mat1

!        decide if a pct is needed based on the database in "decidepctr"
         call decidepctr(do_pct_prp_op,prop_op_name(i),print_lvl)
         scr_mat2 = 0
 
         if(prop_op_nopct) do_pct_prp_op = .false.

         if(do_pct_prp_op)then

           call perform_pct_saao_bas(scr_mat1,            &
                                     pctmat,              &
                                     scr_mat2,            &
                                     naosh_ls,            &
                                     nr_ao_L,             &
                                     naosh_all,           &
                                     naosh_L,             &
                                     nfsym,               &
                                     nz,                  &
                                     ipq_off_in,          &
                                     ioff_aomx,           &
                                     op_fer_rep,          &
                                     op_bs_to_fs,         &
                                     op_2c_ll,            &
                                     print_lvl)

         else 

           call pick_LL_block_saao_bas(scr_mat1,            &
                                       scr_mat2,            &
                                       naosh_L,             &
                                       naosh_all,           &
                                       naosh_ls,            &
                                       nr_ao_L,             &
                                       nfsym,               &
                                       nz,                  &
                                       op_fer_rep,          &
                                       op_bs_to_fs)
           call dcopy(nr_ao_L**2 * nz,scr_mat2,1,scr_mat1,1)
         end if

!        step 2c: store f_2c^++ in DIRAC-sorted SA-AO basis on global file X2CMAT
!        ------------------------------------------------------------------------------

!        define unique file label for each property operator f
         call num2str(i,numstr)
         write(flabel,'(a8,a4)') 'prpint2c',numstr

         if(spherical_on)then

           lwrk = is_len_wrk
           allocate(wrk(lwrk))
           call sph_iotc(scr_mat1,scr_mat2,nr_sph,print_lvl,wrk,lwrk)
           deallocate(wrk)
           

           call x2c_write(flabel,scr_mat2,nr_sph**2  * nz,x2c_file_glb)
         else
           call x2c_write(flabel,scr_mat1,nr_ao_L**2 * nz,x2c_file_glb)
         end if

!        reset logical array
         call lset(nz,.false.,quat_array_filled_tmp)

       end do

       deallocate(quat_array_filled_tmp)
       deallocate(scr_mat2)
       deallocate(scr_mat1)
       deallocate(pctmat)

       close(int1_4c_file,status='keep')

  end subroutine pctrafo_op_4cto2c_saao_bas

!**********************************************************************
  subroutine pctrafo_h1bnccorr_4cto2c_saao_bas(naosh_ls,             &
                                               nr_ao_L,              &
                                               naosh_all,            &
                                               naosh_L,              &
                                               nfsym,                &
                                               nz,                   &
                                               ipq_off_in,           &
                                               ioff_aomx,            &
                                               is_fragment,          &
                                               x2c_file_glb,         &
                                               op_bs_to_fs,          &
                                               is_len_wrk,           &
                                               mdirac,               &
                                               print_lvl)
!**********************************************************************
!
!    purpose: picture-change transform the bare-nucleus corrected 
!             1e Dirac Hamiltonian and store it on file. 
!
!             This Hamiltonian is used in the initial SCF step whereas
!             the default Hamiltonian (including AMFI contributions) 
!             will be considered in all subsequent SCF steps. 
!
!             equation: h_{2c_D_bnc_corr}^++ = U+ * h_{4c_D_bnc_corr} * U
!
!----------------------------------------------------------------------
     integer, intent(in)            :: naosh_ls
     integer, intent(in)            :: nr_ao_L
     integer, intent(in)            :: naosh_all(nfsym)
     integer, intent(in)            :: naosh_L(nfsym)
     integer, intent(in)            :: nfsym
     integer, intent(in)            :: nz
     integer, intent(in)            :: ipq_off_in(4,0:7)
     integer, intent(in)            :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)            :: is_fragment
     integer, intent(in)            :: x2c_file_glb
     integer, intent(in)            :: op_bs_to_fs(0:7,1:2)
     integer, intent(in)            :: is_len_wrk
     integer, intent(in)            :: print_lvl
     logical, intent(in)            :: mdirac
!----------------------------------------------------------------------
     real(8), allocatable           :: wrk(:)
     real(8), allocatable           :: scr_mat1(:)
     real(8), allocatable           :: scr_mat2(:)
     real(8), allocatable           :: pctmat(:)
     integer                        :: i
     integer                        :: j
     integer                        :: lwrk
     integer                        :: op_fer_rep
     integer                        :: ioff_pct_mat
     character (len=12)             :: flabel
!**********************************************************************

       lwrk = is_len_wrk

!      allocate first required scratch matrix
       allocate(scr_mat1(naosh_ls**2 * nz))
       scr_mat1 = 0

!      step 1: get the 4c bare-nucleus corrected 1e Dirac Hamiltonian in SA-AO basis
!      -----------------------------------------------------------------------------
       call reset_x2c_cb_onefck(.false.)

       if(mdirac)then
#ifdef MOD_INTEREST
         scr_mat1(1:naosh_ls**2*nz) = 0.0d0
         call process_kinetic_energy_integrals_4comp( scr_mat1 )
         call process_nuclear_attraction_integrals_4comp( scr_mat1 ) 
#endif
         allocate(wrk(lwrk))
       else
#ifdef PRG_DIRAC
         allocate(wrk(lwrk))
         call onefck(scr_mat1,print_lvl,wrk,lwrk)
#endif
       end if

       call bncorr(scr_mat1,wrk,lwrk,print_lvl)

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr_mat1,naosh_ls,naosh_ls,nz,ipq_off_in(1,0),'x2c - 4ch1bnc',6)
       end if

!      reset common block and WORK array
       deallocate(wrk)
       call reset_x2c_cb_onefck(.true.)

!      step 2: read picture-change transformation matrix in SA-AO basis from global file
!      ---------------------------------------------------------------------------------
       j = 0
       do i = 1, nfsym
         j = j + naosh_all(i)*naosh_L(i)*nz
       end do

!      allocate pct-matrix (2nd required matrix)
       allocate(pctmat(j))
       pctmat = 0

!      initialize offset
       ioff_pct_mat = 1

       do i = 1, nfsym
         write(flabel,'(a7,i4,i1)') 'pctmtAO',is_fragment,i
         if(naosh_all(i)*naosh_L(i) > 0) call x2c_read(flabel,pctmat(ioff_pct_mat),naosh_all(i)*naosh_L(i)*nz,     &
                                                       x2c_file_glb)
         ioff_pct_mat  = ioff_pct_mat + naosh_all(i) * naosh_L(i)* nz
       end do

!      step 3: perform the pct: op_2c^++ == U+ * op_4c * U in DIRAC-sorted SA-AO basis
!      -------------------------------------------------------------------------------
!              in : op_4c    --> scr_mat1; U --> pctmat
!              out: op_2c^++ --> scr_mat1

!      allocate third required matrix
       allocate(scr_mat2(naosh_ls**2 * nz))
       scr_mat2   = 0
       op_fer_rep = 0

       call perform_pct_saao_bas(scr_mat1,            &
                                 pctmat,              &
                                 scr_mat2,            &
                                 naosh_ls,            &
                                 nr_ao_L,             &
                                 naosh_all,           &
                                 naosh_L,             &
                                 nfsym,               &
                                 nz,                  &
                                 ipq_off_in,          &
                                 ioff_aomx,           &
                                 op_fer_rep,          &
                                 op_bs_to_fs,         &
                                 .true.,              &
                                 print_lvl)

!      store h_{2c_D_bnc_corr}^++ in DIRAC-sorted SA-AO basis on global file X2CMAT
       i = 0  ! one-step storage of "gerade+ungerade" matrix parts
       write(flabel,'(a11,i1)') 'h1bncc2cAOl',i
       call x2c_write(flabel,scr_mat1,nr_ao_L**2 * nz,x2c_file_glb)

       if(print_lvl > 2)then
         call print_x2cmat(scr_mat1,nr_ao_L,nr_ao_L,nz,ipq_off_in(1,0),'x2c - bare-corrected',6)
       end if

!#define blubb_blubb_so
#ifdef blubb_blubb_so
       scr_mat2   = 0
       scr_mat1   = 0
       op_fer_rep = 0
       open(99,file='blubb_oo',status='old',form='unformatted',&
            access='sequential',action="readwrite",position='rewind')
       read(99) (scr_mat1(i),i=1, naosh_ls**2*nz)
       close(99,status='keep')

       call perform_pct_saao_bas(scr_mat1,            &
                                 pctmat,              &
                                 scr_mat2,            &
                                 naosh_ls,            &
                                 nr_ao_L,             &
                                 naosh_all,           &
                                 naosh_L,             &
                                 nfsym,               &
                                 nz,                  &
                                 ipq_off_in,          &
                                 ioff_aomx,           &
                                 op_fer_rep,          &
                                 op_bs_to_fs,         &
                                 .true.,              &
                                 print_lvl)

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr_mat1,nr_ao_L,nr_ao_L,nz,ipq_off_in(1,0),'x2c - 2fock-tra',6)
       end if
#endif

       deallocate(scr_mat2)
       deallocate(pctmat)
       deallocate(scr_mat1)

  end subroutine pctrafo_h1bnccorr_4cto2c_saao_bas

!**********************************************************************
  subroutine pctrafo_fockmatrices_4cto2c_saao_bas(naosh_ls,                    &
                                                  nr_ao_L,                     &
                                                  naosh_all,                   &
                                                  naosh_L,                     &
                                                  nfsym,                       &
                                                  nz,                          &
                                                  ipq_off_in,                  &
                                                  ioff_aomx,                   &
                                                  is_fragment,                 &
                                                  x2c_file_glb,                &
                                                  op_bs_to_fs,                 &
                                                  nr_2e_fock_matrices,         &
                                                  file_name_1e_fock_matrix,    &
                                                  file_name_2e_fock_matrices,  &
                                                  print_lvl)
!**********************************************************************
!
!    purpose: picture-change transform the 2-electron Fock matrices of 
!             a converged 4c-SCF run and store them on file.
!
!             This transformation is needed for the molecular-mean field
!             approach within the X2C formalism. 
!             For more information see:
!             -------------------------------------------------------------------
!             J. Sikkema et al., JCP, 131, (2009) 124116. Table I:
!             Hamiltonian scheme --> ^2DC(G)^M (U_mol) 
!             "the converged Fock operator defines the transformation to 2C basis".
!             -------------------------------------------------------------------
!
!             basic equation applied in this subroutine: 
!                      f_{2c}^++ = U+ * f_{4c} * U
!----------------------------------------------------------------------
     integer, intent(in)            :: naosh_ls
     integer, intent(in)            :: nr_ao_L
     integer, intent(in)            :: naosh_all(nfsym)
     integer, intent(in)            :: naosh_L(nfsym)
     integer, intent(in)            :: nfsym
     integer, intent(in)            :: nz
     integer, intent(in)            :: ipq_off_in(4,0:7)
     integer, intent(in)            :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)            :: is_fragment
     integer, intent(in)            :: x2c_file_glb
     integer, intent(in)            :: op_bs_to_fs(0:7,1:2)
     integer, intent(in)            :: nr_2e_fock_matrices
     integer, intent(in)            :: print_lvl
     character (len=6), intent(in)  :: file_name_1e_fock_matrix
     character (len=6), intent(in)  :: file_name_2e_fock_matrices
!----------------------------------------------------------------------
     real(8), allocatable           :: scr_mat1(:)
     real(8), allocatable           :: scr_mat2(:)
     real(8), allocatable           :: pctmat(:)
     integer                        :: i
     integer                        :: j
     integer                        :: fh
     integer                        :: fock_mat_offset
     integer                        :: op_fer_rep
     integer                        :: ioff_pct_mat
     character (len=12)             :: flabel
!**********************************************************************

!      step 1: read the picture-change transformation matrix U in SA-AO basis from global file
!      ---------------------------------------------------------------------------------------
       j = 0
       do i = 1, nfsym
         j = j + naosh_all(i)*naosh_L(i)*nz
       end do
!      allocate pct-matrix
       allocate(pctmat(j))
       pctmat = 0

!      initialize offset
       ioff_pct_mat = 1
       do i = 1, nfsym
         write(flabel,'(a7,i4,i1)') 'pctmtAO',is_fragment,i
         call x2c_read(flabel,pctmat(ioff_pct_mat),naosh_all(i)*naosh_L(i)*nz,x2c_file_glb)
         ioff_pct_mat  = ioff_pct_mat + naosh_all(i) * naosh_L(i)* nz
       end do

!      allocate scratch matrices
       allocate(scr_mat1(naosh_ls**2 * nz * nr_2e_fock_matrices))
       allocate(scr_mat2(naosh_ls**2 * nz))
       scr_mat2   = 0
       scr_mat1 = 0

!      step 2: transform the 2-electron fock matrices
!      ----------------------------------------------

!      step 2a: read the 4c 2-electron fock matrices from file
!      -------------------------------------------------------
       fh = 99
       open(fh,file=file_name_2e_fock_matrices,status='old',form='unformatted',  &
            access='sequential',action="readwrite",position='rewind')
       call read_fock_matrices_saao_basis_x2c(scr_mat1,              &
                                              naosh_ls**2*nz,        &
                                              nr_2e_fock_matrices,   &
                                              fh,                    &
                                              print_lvl)

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr_mat1,naosh_ls,naosh_ls,nz,ipq_off_in(1,0),'x2c - 4cfock2',6)
       end if

!      step 2b: perform the pct: f_2c^++ == U+ * f_4c * U
!      --------------------------------------------------
!              in : f_4c    --> scr_mat1; U --> pctmat
!              out: f_2c^++ --> scr_mat1

       op_fer_rep = 0

!      initialize offset for fock matrix array
       fock_mat_offset = 1

       do i = 1, nr_2e_fock_matrices

         call perform_pct_saao_bas(scr_mat1(fock_mat_offset),    &
                                   pctmat,                       &
                                   scr_mat2,                     &
                                   naosh_ls,                     &
                                   nr_ao_L,                      &
                                   naosh_all,                    &
                                   naosh_L,                      &
                                   nfsym,                        &
                                   nz,                           &
                                   ipq_off_in,                   &
                                   ioff_aomx,                    &
                                   op_fer_rep,                   &
                                   op_bs_to_fs,                  &
                                   .true.,                       &
                                   print_lvl)

!        keep track of correct offset
         fock_mat_offset = fock_mat_offset + i * naosh_ls**2 * nz

       end do

!      step 2c: store the picture-change transformed fock matrices on file
!      -------------------------------------------------------------------
       call write_fock_matrices_saao_basis_x2c(scr_mat1,              &
                                               naosh_ls**2*nz,        &
                                               nr_2e_fock_matrices,   &
                                               fh,                    &
                                               print_lvl)
       close(fh,status='keep')

!      step 3: transform the 1-electron fock matrix
!      --------------------------------------------

!      step 2a: read the 4c 1-electron fock matrix from file
!      -----------------------------------------------------
       scr_mat1 = 0
       scr_mat2 = 0
       fh = 99
       open(fh,file=file_name_1e_fock_matrix,status='old',form='unformatted',  &
            access='sequential',action="readwrite",position='rewind')
       call read_fock_matrices_saao_basis_x2c(scr_mat1,              &
                                              naosh_ls**2*nz,        &
                                              1,                     &
                                              fh,                    &
                                              print_lvl)

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr_mat1,naosh_ls,naosh_ls,nz,ipq_off_in(1,0),'x2c - 4cfock1',6)
       end if

!      step 3b: perform the pct: f_2c^++ == U+ * f_4c * U
!      --------------------------------------------------
!              in : f_4c    --> scr_mat1; U --> pctmat
!              out: f_2c^++ --> scr_mat1

       op_fer_rep = 0

!      initialize offset for fock matrix array
       fock_mat_offset = 1

       call perform_pct_saao_bas(scr_mat1(fock_mat_offset),    &
                                 pctmat,                       &
                                 scr_mat2,                     &
                                 naosh_ls,                     &
                                 nr_ao_L,                      &
                                 naosh_all,                    &
                                 naosh_L,                      &
                                 nfsym,                        &
                                 nz,                           &
                                 ipq_off_in,                   &
                                 ioff_aomx,                    &
                                 op_fer_rep,                   &
                                 op_bs_to_fs,                  &
                                 .true.,                       &
                                 print_lvl)

!      step 3c: store the picture-change transformed fock matrix on file
!      -----------------------------------------------------------------
       call write_fock_matrices_saao_basis_x2c(scr_mat1,              &
                                               naosh_ls**2*nz,        &
                                               1,                     &
                                               fh,                    &
                                               print_lvl)
       close(fh,status='keep')

#define blubb_blubb_aooso
#ifdef blubb_blubb_aooso
       scr_mat2   = 0
       scr_mat1   = 0
       op_fer_rep = 0
       open(99,file='blubb_oo',status='old',form='unformatted',&
            access='sequential',action="readwrite",position='rewind')
       read(99) (scr_mat1(i),i=1, naosh_ls**2*nz)
       close(99,status='keep')

       call make_two_e_fock_so_exclusive(scr_mat1,nr_ao_L,nr_ao_L,nz)
!      debug print
!      if(print_lvl > 2)then
         call print_x2cmat(scr_mat1,naosh_ls,naosh_ls,nz,ipq_off_in(1,0),'x2c - 2fock-before tra',6)
!      end if
       call perform_pct_saao_bas(scr_mat1(fock_mat_offset),    &
                                 pctmat,                       &
                                 scr_mat2,                     &
                                 naosh_ls,                     &
                                 nr_ao_L,                      &
                                 naosh_all,                    &
                                 naosh_L,                      &
                                 nfsym,                        &
                                 nz,                           &
                                 ipq_off_in,                   &
                                 ioff_aomx,                    &
                                 op_fer_rep,                   &
                                 op_bs_to_fs,                  &
                                 .true.,                       &
                                 print_lvl)
!      debug print
!      if(print_lvl > 2)then
         call print_x2cmat(scr_mat1,nr_ao_L,nr_ao_L,nz,ipq_off_in(1,0),'x2c - 2fock-tra',6)
!      end if
       open(99,file='blubb_oo',status='replace',form='unformatted',&
            access='sequential',action="readwrite",position='rewind')
       write(99) (scr_mat1(i),i=1, nr_ao_L**2*4)
       close(99,status='keep')
#undef blubb_blubb_aooso
#endif
       deallocate(scr_mat2)
       deallocate(scr_mat1)
       deallocate(pctmat)

  end subroutine pctrafo_fockmatrices_4cto2c_saao_bas

!**********************************************************************
  subroutine perform_pct_saao_bas(op_matrix_2c_ll,     &
                                  pctmat,              &
                                  scr_mat1,            &
                                  naosh_ls,            &
                                  nr_ao_L,             &
                                  naosh_all,           &
                                  naosh_L,             &
                                  nfsym,               &
                                  nz,                  &
                                  ipq_off_in,          &
                                  ioff_aomx,           &
                                  op_fer_rep,          &
                                  op_bs_to_fs,         &
                                  op_2c_ll,            &
                                  print_lvl)   
!**********************************************************************
!
!    purpose: perform the picture-change transformation of a four-component
!             (property) operator X_4c using the total pct-matrix U:
!
!             ===================================
!             equation: X_2c^++ = U^+ * X_4c * U
!             ===================================
!
!             if op_2c_ll == .true. : return only the upper op_2c^++ block (remaining parts zeroed out)
!             if op_2c_ll == .false.: return full four-component matrix block containing the
!                                     pct-transformed op_2c^++ block
!
!             input : 4c-operator matrix    --> op_matrix_2c_ll; pct-matrix --> pctmat 
!             output: 2c^++-operator matrix --> op_matrix_2c_ll
!
!----------------------------------------------------------------------
     real(8), intent(inout) :: op_matrix_2c_ll(naosh_ls**2,nz)
     real(8), intent(inout) :: scr_mat1(naosh_ls**2,nz)
     real(8), intent(in)    :: pctmat(*)
     integer, intent(in)    :: naosh_ls
     integer, intent(in)    :: nr_ao_L
     integer, intent(in)    :: naosh_all(nfsym)
     integer, intent(in)    :: naosh_L(nfsym)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(in)    :: ipq_off_in(4,0:7)
     integer, intent(in)    :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)    :: print_lvl
     integer, intent(in)    :: op_fer_rep
     integer, intent(in)    :: op_bs_to_fs(0:7,1:2)
     logical, intent(in)    :: op_2c_ll
!----------------------------------------------------------------------
     integer                :: i
     integer                :: k
     integer                :: ik
     integer                :: nbas1_f_i
     integer                :: nbasL_f_i
     integer                :: nbas1_f_ik
     integer                :: nbasL_f_ik
     integer                :: op_sym
     integer                :: ioff_pct_mat_dag
     integer                :: ioff_pct_mat
!**********************************************************************

!      initialize offset for pctmat+
       ioff_pct_mat_dag = 1

!      step 1: perform pct in DIRAC-sorted SA-AO basis
!      ------------------------------------------------------------------------------

!      operator: map fermion to boson symmetry
       op_sym = op_bs_to_fs(op_fer_rep,1)

       do i = 1, nfsym

         ik = mod(i+op_sym,2) + 1

!        set dimensions for matrices used in the subroutine
         nbas1_f_i  = naosh_all(i)
         nbasL_f_i  = naosh_L(i)
         nbas1_f_ik = naosh_all(ik)
         nbasL_f_ik = naosh_L(ik)

         if(naosh_L(i) > 0 .and. naosh_L(ik) > 0)then

!          determine offset for pctmat
           ioff_pct_mat = 1
           if(ik == 2) ioff_pct_mat = 1 + naosh_L(1) * naosh_all(1) * nz

!          perform picture-change transformation: 
!          X_2c = pctmat+ * X_4c * pctmat
!          scr_mat1 = pctmat+ * op_matrix_2c_ll * pctmat
           call qtrans90('AOMO','S',val_d0,                      &
                         nbas1_f_i,                              &
                         nbas1_f_ik,                             &
                         nbasL_f_i,                              &
                         nbasL_f_ik,                             &
                         op_matrix_2c_ll(1+ioff_aomx(i,ik),1),   &
                         naosh_ls,                               &
                         naosh_ls,                               &
                         nz,                                     &
                         ipq_off_in(1,op_fer_rep),               &
                         scr_mat1(1+ioff_aomx(i,ik),1),          &
                         naosh_ls,                               &
                         naosh_ls,                               &
                         nz,                                     &
                         ipq_off_in(1,op_fer_rep),               &
                         pctmat(ioff_pct_mat_dag),               &
                         nbas1_f_i,                              &
                         nbasL_f_i,                              &
                         nz,                                     &
                         ipq_off_in(1,0),                        &
                         pctmat(ioff_pct_mat),                   &
                         nbas1_f_ik,                             &
                         nbasL_f_ik,                             &
                         nz,                                     &
                         ipq_off_in(1,0),                        &
                         print_lvl)
          
         end if

!        update offset for pctmat+
         ioff_pct_mat_dag  = ioff_pct_mat_dag + naosh_L(i) * naosh_all(i) * nz

       end do

!      step 2: return pct-operator in DIRAC-sorted SA-AO basis on array op_matrix_2c_ll
!      ------------------------------------------------------------------------------
!              case a. op_2c_ll == .true. : return upper op_2c++ block (remaining parts zeroed out)
!              case b. op_2c_ll == .false.: return full 4c block with upper op_2c++ block containing the pct-operator

       if(op_2c_ll)then
         op_matrix_2c_ll = 0
         call pick_LL_block_saao_bas(scr_mat1,            &
                                     op_matrix_2c_ll,     &
                                     naosh_L,             &
                                     naosh_all,           &
                                     naosh_ls,            &
                                     nr_ao_L,             &
                                     nfsym,               &
                                     nz,                  &
                                     op_fer_rep,          &
                                     op_bs_to_fs)
       else 
         call dcopy(naosh_ls**2 * nz,scr_mat1,1,op_matrix_2c_ll,1)
       end if

  end subroutine perform_pct_saao_bas

!**********************************************************************

end module x2c_pct_ao
