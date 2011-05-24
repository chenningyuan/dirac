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
!
!
! module containing the functionality to perform the picture-change
! transformation of the four-component MO coefficients following the X2C formalism. 
!
! written by sknecht may 2011
!
module x2c_pct_mo_coefficients

  use memory_allocator
  use x2c_fio
  use x2c_utils, only:                 &
      print_x2cmat

  implicit none

  public x2c_pctrafo_mo_coefficients_driver

  private

  real(8), parameter :: val_d1      =  1.0d0
  real(8), parameter :: val_d0      =  0.0d0
  real(8), parameter :: val_dm1     = -1.0d0

contains

!**********************************************************************
  subroutine x2c_pctrafo_mo_coefficients_driver(mo_coeff_file_name,  &
                                                naosh_all,           &
                                                naosh_L,             &
                                                norb_dim,            &
                                                nesh_dim,            &
                                                npsh_dim,            &
                                                nfsym,               &
                                                nz,                  &
                                                ipq_off_in,          &
                                                x2c_file_glb,        &
                                                print_lvl)   
!**********************************************************************
!
!    purpose: driver routine for the picture-change transformation of 
!             the 4c-MO coefficients. 
!             this path of the X2C procedure is only active when the 
!             4c-fock operator has been chosen as the defining hamiltonian 
!             for the decoupling step.
!             --> this module is part of the molecular-mean-field option within 
!             the X2C framework.
!
!----------------------------------------------------------------------
     integer, intent(in)             :: nfsym
     integer, intent(in)             :: nz
     integer, intent(in)             :: naosh_all(nfsym)
     integer, intent(in)             :: naosh_L(nfsym)
     integer, intent(in)             :: norb_dim(nfsym)
     integer, intent(in)             :: nesh_dim(nfsym)
     integer, intent(in)             :: npsh_dim(nfsym)
     integer, intent(in)             :: ipq_off_in(4,0:7)
     integer, intent(in)             :: x2c_file_glb
     character (len= 6), intent(in)  :: mo_coeff_file_name
     integer, intent(in)             :: print_lvl
!----------------------------------------------------------------------
     real(8), allocatable            :: pctmat(:)
     real(8), allocatable            :: mo_coeff_4c(:)
     real(8), allocatable            :: mo_coeff_2c(:)
     real(8), allocatable            :: eigen_values(:)
     integer, allocatable            :: boson_irreps_jz_eigenvalues(:)
     integer                         :: i, j
     integer                         :: is_fragment_tot
!**********************************************************************

!      initialize the fragment value - in the mmf-case it has (by definition) to be zero.
       is_fragment_tot = 0

!      step 0: set dimensions for allocation of scratch space + allocate
!      -----------------------------------------------------------------
       i = 1
       call alloc(pctmat,i)
       call alloc(mo_coeff_4c,i)
       call alloc(mo_coeff_2c,i)
       call alloc(eigen_values,i)
       call alloc(boson_irreps_jz_eigenvalues,i)

!      step 1: picture-change transform the MO-coefficients
!      ----------------------------------------------------


!      step 2: save the picture-change transformed MO-coefficients
!      -----------------------------------------------------------

!      step 3: release memory
!      ----------------------
       call dealloc(boson_irreps_jz_eigenvalues)
       call dealloc(eigen_values)
       call dealloc(mo_coeff_2c)
       call dealloc(mo_coeff_4c)
       call dealloc(pctmat)

  end subroutine x2c_pctrafo_mo_coefficients_driver

!**********************************************************************
  subroutine x2c_pctrafo_4cto2c_mo_coefficients()
!**********************************************************************
#ifdef undef
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
     logical, intent(in)            :: spherical_on
     character (len=16), intent(in) :: prop_op_name(num_prp_op)
!----------------------------------------------------------------------
     real(8), pointer               :: wrk(:)
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
     logical, allocatable           :: quat_array_filled_tmp(:)
     character (len=12)             :: flabel
     character (len=4)              :: numstr
!**********************************************************************

!      initialize offset
       ioff_pct_mat = 1

       j = 0
       do i = 1, nfsym
!        print'(a,2i6)', 'naosh_all(i), naosh_L(i)  =',naosh_all(i),naosh_L(i)
         j = j + naosh_all(i)*naosh_L(i)*nz
       end do
!      allocate pct-matrix
       call alloc(pctmat,j)
       pctmat = 0

!      step 1: read picture-change transformation matrix U in SA-AO basis from global file
!      ---------------------------------------------------------------------------------

!      print '(a,i0)', 'nfsym       == ',nfsym
!      print '(a,i0)', 'nz          == ',nz
!      print '(a,i0)', 'is_fragment == ',is_fragment

       do i = 1, nfsym
         write(flabel,'(a7,i4,i1)') 'pctmtAO',is_fragment,i
         call x2c_read(flabel,pctmat(ioff_pct_mat),naosh_all(i)*naosh_L(i)*nz,x2c_file_glb)
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
       call alloc(scr_mat1,naosh_ls**2 * nz)
       call alloc(scr_mat2,naosh_ls**2 * nz)
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
           call sel_wrk_x2c(wrk,lwrk)
           call sph_iotc(scr_mat1,scr_mat2,nr_sph,print_lvl,wrk,lwrk)
           call desel_wrk_x2c(wrk,lwrk)

           call x2c_write(flabel,scr_mat2,nr_sph**2  * nz,x2c_file_glb)
         else
           call x2c_write(flabel,scr_mat1,nr_ao_L**2 * nz,x2c_file_glb)
         end if

!        reset logical array
         call lset(nz,.false.,quat_array_filled_tmp)

       end do

       deallocate(quat_array_filled_tmp)
       call dealloc(scr_mat2)
       call dealloc(scr_mat1)
       call dealloc(pctmat)

       close(int1_4c_file,status='keep')
#endif
  end subroutine x2c_pctrafo_4cto2c_mo_coefficients
!**********************************************************************

end module x2c_pct_mo_coefficients
