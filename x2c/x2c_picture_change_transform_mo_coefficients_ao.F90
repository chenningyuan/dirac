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
  use x2c_cb_interface, only:          &
      renew_x2c_cb_orb_shell_dim

  implicit none

  public x2c_pctrafo_mo_coefficients_driver
  public x2c_molecular_mean_field_mos

  private

  real(8), parameter :: val_d1      =  1.0d0
  real(8), parameter :: val_d0      =  0.0d0
  real(8), parameter :: val_dm1     = -1.0d0

contains

!**********************************************************************
  subroutine x2c_pctrafo_mo_coefficients_driver(mo_coeff_file_name,  &
                                                pct_mocoeff_strategy,&
                                                nr_cmo_q,            &
                                                naosh_ls,            &
                                                nr_ao_L,             &
                                                naosh_all,           &
                                                naosh_L,             &
                                                ioff_aomx,           &
                                                dimension_for_2mat,  &
                                                norb_dim,            &
                                                nesh_dim,            &
                                                npsh_dim,            &
                                                nfsym,               &
                                                nz,                  &
                                                ipq_off_in,          &
                                                x2c_file_glb,        &
                                                linear_sym,          &
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
     integer, intent(in)             :: pct_mocoeff_strategy
     integer, intent(in)             :: nfsym
     integer, intent(in)             :: nz
     integer, intent(in)             :: nr_cmo_q
     integer, intent(in)             :: naosh_ls
     integer, intent(in)             :: nr_ao_L
     integer, intent(in)             :: naosh_all(nfsym)
     integer, intent(in)             :: naosh_L(nfsym)
     integer, intent(in)             :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)             :: dimension_for_2mat
     integer, intent(in)             :: norb_dim(nfsym)
     integer, intent(in)             :: nesh_dim(nfsym)
     integer, intent(in)             :: npsh_dim(nfsym)
     integer, intent(in)             :: ipq_off_in(4,0:7)
     integer, intent(in)             :: x2c_file_glb
     logical, intent(in)             :: linear_sym
     character (len= 6), intent(in)  :: mo_coeff_file_name
     integer, intent(in)             :: print_lvl
!----------------------------------------------------------------------
     real(8), allocatable            :: scratch1(:)
     real(8), allocatable            :: scratch2(:)
     real(8), allocatable            :: scratch3(:)
     real(8), allocatable            :: scratch4(:)
     integer, allocatable            :: scratch5(:)
     integer                         :: i, j, k, l
     integer                         :: is_fragment_tot
!**********************************************************************

!      initialize the fragment value - in the mmf-case it has (by definition) to be zero.
       is_fragment_tot = 0
       j               = 0
       k               = 0
       l               = 0

!      step 1: set dimensions for allocation of scratch space + allocate 
!      (for simplicity we keep it independent from the strategy chosen below)
!      ----------------------------------------------------------------------
       do i = 1, nfsym
         k = k + norb_dim(i)
       end do
       j = naosh_ls*naosh_ls*nz
       l = max(naosh_ls*naosh_ls*nz,k*4*nz+8) ! see rsjaco and qdiag memory requests 
                                              ! for further explanation of this odd dimensioning

       call alloc(scratch1,dimension_for_2mat)
       call alloc(scratch2,j)
       call alloc(scratch3,l)
       call alloc(scratch4,k)
       call alloc(scratch5,k)
       scratch1 = 0       
       scratch2 = 0
       scratch3 = 0
       scratch4 = 0
       scratch5 = 0

       select case(pct_mocoeff_strategy)

         case(1) 

!        path 1: directly picture-change transform the MO-coefficients
!        -------------------------------------------------------------
!        this branch is currently inactive as there seems to be a problem with the 
!        current implementation scheme. stefan - june 2011
!        objects handled inside this routine:
!        in:
!            - 4c-MO coefficients residing on DFCOEF
!            - picture-change-transformation matrix U in SA-AO basis
!        out:
!            - 2c-MO coefficients residing on DFCOEF
         call x2c_pctrafo_4cto2c_mo_coefficients(mo_coeff_file_name,           &
                                                 scratch1,                     &
                                                 scratch2,                     &
                                                 scratch3,                     &
                                                 scratch4,                     &
                                                 scratch5,                     &
                                                 nr_cmo_q,                     &
                                                 naosh_ls,                     &
                                                 nr_ao_L,                      &
                                                 naosh_all,                    &
                                                 naosh_L,                      &
                                                 norb_dim,                     &
                                                 nesh_dim,                     &
                                                 npsh_dim,                     &
                                                 nfsym,                        &
                                                 nz,                           &
                                                 ipq_off_in,                   &
                                                 x2c_file_glb,                 &
                                                 is_fragment_tot,              &
                                                 print_lvl)   
         case(2) 

!        path 2: diagonalize the picture-change transformed 2c-Fock operator to obtain the pc-transformed MO-coefficients 
!        ----------------------------------------------------------------------------------------------------------------
!        objects handled inside this routine:
!        in:
!            - picture-change-transformed Fock-operator in SA-AO basis residing on file X2CMAT
!            - Lowdin-transformation matrix in 2c-form                 residing on file X2CMAT
!        out:
!            - 2c-MO coefficients in orthonormal basis residing on X2CMAT
!            - 2c-eigenvalues                          residing on X2CMAT
!            - 2c-boson irrep / j_z information        residing on X2CMAT
         call x2c_diagonalize_2c_fock_operator_onmo(scratch1,                     &
                                                    scratch2,                     &
                                                    scratch3,                     &
                                                    scratch4,                     &
                                                    scratch5,                     &
                                                    naosh_ls,                     &
                                                    nr_ao_L,                      &
                                                    naosh_all,                    &
                                                    naosh_L,                      &
                                                    ioff_aomx,                    &
                                                    norb_dim,                     &
                                                    nesh_dim,                     &
                                                    nfsym,                        &
                                                    nz,                           &
                                                    ipq_off_in,                   &
                                                    l,                            &
                                                    x2c_file_glb,                 &
                                                    is_fragment_tot,              &
                                                    linear_sym,                   &
                                                    print_lvl)   
       end select

!      step 2: release memory
!      ----------------------
       call dealloc(scratch5)
       call dealloc(scratch4)
       call dealloc(scratch3)
       call dealloc(scratch2)
       call dealloc(scratch1)

  end subroutine x2c_pctrafo_mo_coefficients_driver

!**********************************************************************
  subroutine x2c_pctrafo_4cto2c_mo_coefficients(mo_coeff_file_name,          &
                                                pctmat,                      &
                                                mo_coeff_4c,                 &
                                                mo_coeff_2c,                 &
                                                eigen_values,                &
                                                boson_irreps_jz_eigenvalues, &
                                                nr_cmo_q,                    &
                                                naosh_ls,                    &
                                                nr_ao_L,                     &
                                                naosh_all,                   &
                                                naosh_L,                     &
                                                norb_dim,                    &
                                                nesh_dim,                    &
                                                npsh_dim,                    &
                                                nfsym,                       &
                                                nz,                          &
                                                ipq_off_in,                  &
                                                x2c_file_glb,                &
                                                is_fragment,                 &
                                                print_lvl)   
!**********************************************************************
!
!    purpose: 
!
!             equation: mo_coeff_2c^LL = U+ * mo_coeff_4c * U
!
!----------------------------------------------------------------------
     real(8), intent(inout)         :: eigen_values(*)
     real(8), intent(inout)         :: mo_coeff_2c(*)
     real(8), intent(inout)         :: mo_coeff_4c(*)
     real(8), intent(inout)         :: pctmat(*)
     integer, intent(inout)         :: boson_irreps_jz_eigenvalues(*)
     integer, intent(in)            :: nr_cmo_q
     integer, intent(in)            :: nfsym
     integer, intent(in)            :: nz
     integer, intent(in)            :: naosh_ls
     integer, intent(in)            :: nr_ao_L
     integer, intent(in)            :: naosh_all(nfsym)
     integer, intent(in)            :: naosh_L(nfsym)
     integer, intent(in)            :: norb_dim(nfsym)
     integer, intent(in)            :: nesh_dim(nfsym)
     integer, intent(in)            :: npsh_dim(nfsym)
     integer, intent(in)            :: ipq_off_in(4,0:7)
     integer, intent(in)            :: is_fragment
     integer, intent(in)            :: x2c_file_glb
     integer, intent(in)            :: print_lvl
     character (len= 6), intent(in) :: mo_coeff_file_name
!----------------------------------------------------------------------
     integer                        :: i, j
     integer                        :: ioff_pct_mat
     integer                        :: norb1_f
     integer                        :: n1esh_f
     integer                        :: n1psh_f
     integer                        :: nbasL_f
     integer                        :: nbas1_f
     integer                        :: ioff_mo_coeff_4c
     integer                        :: ioff_mo_coeff_2c
     integer                        :: offset_4c      
     integer                        :: offset_2c
     integer                        :: ndim_2c_mo
     integer                        :: ndim_2c_ev
     integer                        :: ndim_2c_bj
     integer                        :: fh
     real(8)                        :: total_energy
     character (len=12)             :: flabel
     character (len=24)             :: fdate
     character (len=74)             :: info_text
     integer, allocatable           :: dimension_info(:,:)
!**********************************************************************

!      initialize offset
       ioff_pct_mat = 1

!      step 1: read picture-change transformation matrix U in SA-AO basis from global file
!      ---------------------------------------------------------------------------------

       do i = 1, nfsym
         write(flabel,'(a7,i4,i1)') 'pctmtAO',is_fragment,i
         call x2c_read(flabel,pctmat(ioff_pct_mat),naosh_all(i)*naosh_L(i)*nz,x2c_file_glb)
         ioff_pct_mat  = ioff_pct_mat + naosh_all(i) * naosh_L(i)* nz
       end do

!      debug print
       if(print_lvl > 2)then
         ioff_pct_mat = 1
         do i = 1, nfsym
           write(info_text,'(a27,i1)') 'x2c - pctmat_saao - xcmo - ',i
           call print_x2cmat(pctmat(ioff_pct_mat),naosh_all(i),naosh_L(i),nz,ipq_off_in,info_text,6)
           ioff_pct_mat  = ioff_pct_mat + naosh_all(i) * naosh_L(i)* nz
         end do
       end if

!      read MO coefficients
       fh = 99
       open(fh,file=mo_coeff_file_name,status='old',form='unformatted',access='sequential',position='rewind') 
       call reacmo_no_work(fh,mo_coeff_file_name,mo_coeff_4c,eigen_values,boson_irreps_jz_eigenvalues,total_energy,15)
       close(fh,status='delete')

!      step 3: perform pct of MO coefficients: MO_coeff_2c^++ == U+ * MO_coeff_4c
!      --------------------------------------------------------------------------
       ioff_pct_mat     = 1
       ioff_mo_coeff_4c = 1
       ioff_mo_coeff_2c = 1

       do i = 1, nfsym

         norb1_f          = norb_dim(i)
         n1psh_f          = npsh_dim(i)
         n1esh_f          = nesh_dim(i)
         nbas1_f          = naosh_all(i)
         nbasL_f          = naosh_L(i)
         ioff_mo_coeff_4c = ioff_mo_coeff_4c + nbas1_f*n1psh_f*nz

         if(nesh_dim(i) > 0)then
           call qgemm(nbasL_f,                                  &
                      n1esh_f,                                  &
                      nbas1_f,                                  &
                      val_d1,                                   &
                      'H','N',                                  &
                      ipq_off_in,                               &
                      pctmat(ioff_pct_mat),                     &
                      nbas1_f,                                  &
                      nbasL_f,                                  &
                      nz,                                       &
                      'N','N',                                  &
                      ipq_off_in,                               &
                      mo_coeff_4c(ioff_mo_coeff_4c),            &
                      nbas1_f,                                  &
                      n1esh_f,                                  &
                      nz,                                       &
                      val_d0,                                   &
                      ipq_off_in,                               &
                      mo_coeff_2c(ioff_mo_coeff_2c),            &
!                     nbasL_f,                                  &
                      nbas1_f,                                  &
                      n1esh_f,                                  &
                      nz)

!          zero LS part of coefficients
           do j = 1, n1esh_f
             call dzero(mo_coeff_2c(ioff_mo_coeff_2c+((j-1)*nbas1_f+nbasL_f)),nbas1_f-nbasL_f)
           end do
           if(print_lvl > 2)then
             print *, 'dimensions: nbas1_f, nbasL_f, n1esh_f, n1psh_f: ',nbas1_f, nbasL_f, n1esh_f, n1psh_f
             print *, '4c-MO coefficients (all) ==> '
             write(info_text,'(a27,i1)') 'x2c - 4c-MO (all) - xcmo - ',i
             call print_x2cmat(mo_coeff_4c,nbas1_f,norb1_f,nz,ipq_off_in,info_text,6)
       
             print *, '2c-MO coefficients ==> '
             write(info_text,'(a27,i1)') 'x2c - 2c-MO (all) - xcmo - ',i
             call print_x2cmat(mo_coeff_2c(ioff_mo_coeff_2c),nbas1_f,n1esh_f,nz,ipq_off_in,info_text,6)
!            call print_x2cmat(mo_coeff_2c(ioff_mo_coeff_2c),nbasL_f,n1esh_f,nz,ipq_off_in,info_text,6)
           end if
         end if

!        set new offsets in matrices
         ioff_mo_coeff_2c =  ioff_mo_coeff_2c + nbas1_f*n1esh_f*nz
!        ioff_mo_coeff_2c =  ioff_mo_coeff_2c + nbasL_f*n1esh_f*nz
         ioff_mo_coeff_4c =  ioff_mo_coeff_4c + nbas1_f*n1esh_f*nz
         ioff_pct_mat     =  ioff_pct_mat     + nbas1_f*nbasL_f*nz
       end do

!      step 4: write 2c MO coefficients back to file
!      ---------------------------------------------
 
       call alloc(dimension_info,3,nfsym)
!      if(.false.)then
       if(.true.)then
       print *, 'compress...'
!      4a. compress arrays to large-components
       offset_4c        = 1
       offset_2c        = 1
       ioff_mo_coeff_2c = 1
       do i = 1, nfsym

         nbas1_f = naosh_all(i)
         nbasL_f = naosh_L(i)

         call dcopy(nesh_dim(i),eigen_values(offset_4c+npsh_dim(i)),1,mo_coeff_4c,1)
         call dcopy(nesh_dim(i),mo_coeff_4c,1,eigen_values(offset_2c),1)
         call icopy(nesh_dim(i),boson_irreps_jz_eigenvalues(offset_4c+npsh_dim(i)),1,  &
                                boson_irreps_jz_eigenvalues(offset_2c),1)
         ioff_mo_coeff_4c = 1

!        do j = 1, nesh_dim(i)
!          call dcopy(nbasL_f,mo_coeff_2c(ioff_mo_coeff_2c+((j-1)*nbas1_f)),1,         &
!                             mo_coeff_4c(ioff_mo_coeff_4c+((j-1)*nbasL_f)),1)
!          call dcopy(nbasL_f,mo_coeff_4c(ioff_mo_coeff_4c+((j-1)*nbasL_f)),1,         &
!                             mo_coeff_2c(ioff_mo_coeff_2c+((j-1)*nbasL_f)),1)
!        end do

         offset_2c           = offset_2c + nesh_dim(i)
         offset_4c           = offset_4c + norb_dim(i)
         ioff_mo_coeff_2c    = ioff_mo_coeff_2c + nbas1_f*nesh_dim(i)*nz

         dimension_info(1,i) = 0
         dimension_info(2,i) = nesh_dim(i)
         dimension_info(3,i) = naosh_all(i)
!        dimension_info(3,i) = naosh_L(i)
       end do
       ndim_2c_mo = ioff_mo_coeff_2c - 1
       ndim_2c_ev = 0
       ndim_2c_bj = 0
       do i = 1, nfsym
         ndim_2c_ev = ndim_2c_ev + nesh_dim(i) 
         ndim_2c_bj = ndim_2c_bj + nesh_dim(i)
       end do
       print *, 'compress done...'
       else
       print *, 'compress...'
!      4a. compress arrays to large-components
       offset_4c        = 1
       offset_2c        = 1
       ioff_mo_coeff_2c = 1
       do i = 1, nfsym

         nbas1_f = naosh_all(i)
         nbasL_f = naosh_L(i)

         call dcopy(nesh_dim(i),eigen_values(offset_4c+npsh_dim(i)),1,mo_coeff_4c,1)
         call dcopy(nesh_dim(i),mo_coeff_4c,1,eigen_values(offset_2c),1)
         call icopy(nesh_dim(i),boson_irreps_jz_eigenvalues(offset_4c+npsh_dim(i)),1,  &
                                boson_irreps_jz_eigenvalues(offset_2c),1)
         ioff_mo_coeff_4c = 1

!        do j = 1, nesh_dim(i)
!          call dcopy(nbasL_f,mo_coeff_2c(ioff_mo_coeff_2c+((j-1)*nbas1_f)),1,         &
!                             mo_coeff_4c(ioff_mo_coeff_4c+((j-1)*nbasL_f)),1)
!          call dcopy(nbasL_f,mo_coeff_4c(ioff_mo_coeff_4c+((j-1)*nbasL_f)),1,         &
!                             mo_coeff_2c(ioff_mo_coeff_2c+((j-1)*nbasL_f)),1)
!        end do

         offset_2c           = offset_2c + nesh_dim(i)
         offset_4c           = offset_4c + norb_dim(i)
         ioff_mo_coeff_2c    = ioff_mo_coeff_2c + nbasL_f*nesh_dim(i)*nz

         dimension_info(1,i) = 0
         dimension_info(2,i) = nesh_dim(i)
         dimension_info(3,i) = naosh_L(i)
       end do
         ndim_2c_mo = ioff_mo_coeff_2c - 1
         ndim_2c_ev = 0
         ndim_2c_bj = 0
         do i = 1, nfsym
           ndim_2c_ev = ndim_2c_ev + nesh_dim(i) 
           ndim_2c_bj = ndim_2c_bj + nesh_dim(i)
         end do
       print *, 'ndim_2c_mo is: ',ndim_2c_mo
       end if
       print *, 'compress done...'

       fh = 99
       open(fh,file=mo_coeff_file_name,status='new',form='unformatted',access='sequential',position='rewind') 

       info_text(1:50)  = ' X2Cmod - molecular-mean field 2c MOs              '
       info_text(51:74) = fdate()
       write(fh) info_text,nfsym,((dimension_info(i,j),i = 1,3),j=1,nfsym),total_energy
       write(fh) (mo_coeff_2c(i),i=1,ndim_2c_mo)
       write(fh) (eigen_values(i),i=1,ndim_2c_ev)
       write(fh) (boson_irreps_jz_eigenvalues(i),i=1,ndim_2c_bj)

       close(fh,status='keep')
       call dealloc(dimension_info)

  end subroutine x2c_pctrafo_4cto2c_mo_coefficients
!**********************************************************************

  subroutine x2c_diagonalize_2c_fock_operator_onmo(xmat,                     &
                                                   scratch2,                 &
                                                   scratch3,                 &
                                                   scratch4,                 &
                                                   scratch5,                 &
                                                   naosh_ls,                 &
                                                   nr_ao_L,                  &
                                                   naosh_all,                &
                                                   naosh_L,                  &
                                                   ioff_aomx,                &
                                                   norb_dim,                 &
                                                   nesh_dim,                 &
                                                   nfsym,                    &
                                                   nz,                       &
                                                   ipq_off_in,               &
                                                   lwork,                    &
                                                   x2c_file_glb,             &
                                                   is_fragment,              &
                                                   linear_sym,               &
                                                   print_lvl)   
!**********************************************************************
!
!    purpose: 
!
!----------------------------------------------------------------------
     real(8), intent(inout)         :: xmat(*)
     real(8), intent(inout)         :: scratch2(*)
     real(8), intent(inout)         :: scratch3(*)
     real(8), intent(inout)         :: scratch4(*)
     integer, intent(inout)         :: scratch5(*)
     integer, intent(in)            :: nfsym
     integer, intent(in)            :: nz
     integer, intent(in)            :: naosh_ls
     integer, intent(in)            :: nr_ao_L
     integer, intent(in)            :: naosh_all(nfsym)
     integer, intent(in)            :: naosh_L(nfsym)
     integer, intent(in)            :: ioff_aomx(nfsym,nfsym)
     integer, intent(in)            :: norb_dim(nfsym)
     integer, intent(in)            :: nesh_dim(nfsym)
     integer, intent(in)            :: ipq_off_in(4,0:7)
     integer, intent(inout)         :: lwork
     integer, intent(in)            :: is_fragment
     integer, intent(in)            :: x2c_file_glb
     integer, intent(in)            :: print_lvl
     logical, intent(in)            :: linear_sym
!----------------------------------------------------------------------
     integer                        :: i, j, l
     integer                        :: norb1_f
     integer                        :: n1esh_f
     integer                        :: nbasL_f
     integer                        :: nbas1_f
     integer                        :: offset_1
     integer                        :: offset_2
     integer                        :: nz_lowdin = 1
     integer                        :: fh = 99
     integer                        :: nborb_dim(4,2,0:2)
     logical                        :: fndlab12
     character (len=12)             :: flabel
     character (len=74)             :: info_text
     real(8), allocatable           :: original_eigen_values(:)
     real(8)                        :: dummy
!**********************************************************************

!      step 1: read from global file: 
!              a. the 2c-Fock operator in SA-AO basis
!              b. the 2c-Lowdin transformation matrix
!      -----------------------------------------------------------------
!      a. 2c-Fock operator
       do i = 1, nfsym
         nbas1_f  = naosh_all(i)
         if(naosh_all(i) > 0)then
           write(flabel,'(a7,i4,i1)') 'h12cAOn',1,i
           call x2c_read(flabel,xmat(1+ioff_aomx(i,i)),nbas1_f*naosh_ls*nz,x2c_file_glb)
         end if
       end do

!      b. 2c-Lowdin transformation matrix
       write(flabel,'(a7,i4,i1)') '2cLWAOg',1,0
       if(linear_sym) write(flabel,'(a7,i4,i1)') '2cLWAOl',1,0
       call x2c_read(flabel,scratch3,(naosh_all(1)*nesh_dim(1)+naosh_all(2)*nesh_dim(2)),x2c_file_glb)

!      debug print
       if(print_lvl > 2)then
         offset_1 = 1
         do i = 1, nfsym

           write(info_text,'(a27,i1)') 'x2c - 2c-fock-saao - fsym =',i
           call print_x2cmat(xmat(1+ioff_aomx(i,i)),naosh_all(i),naosh_ls,nz,ipq_off_in,info_text,6,.true.)

           write(info_text,'(a27,i1)') 'x2c - 2c-Lowdin    - fsym =',i
           call print_x2cmat(scratch3(offset_1),naosh_all(i),nesh_dim(i),1,ipq_off_in,info_text,6)

           offset_1  = offset_1 + naosh_all(i) * nesh_dim(i)
         end do
       end if

!      step 2: perform the transformation to 2c-onMO basis by means of the Lowdin transformation matrix
!      ------------------------------------------------------------------------------------------------
       offset_1 = 1
       offset_2 = 1

       do i = 1, nfsym

         norb1_f          = norb_dim(i)
         n1esh_f          = nesh_dim(i)
         nbas1_f          = naosh_all(i)
         nbasL_f          = naosh_L(i)

         if(norb_dim(i) > 0)then
           call qtrans90('AOMO','S',val_d0,      &
                         nbasL_f,                &
                         nbasL_f,                &
                         n1esh_f,                &
                         n1esh_f,                &
                         xmat(1+ioff_aomx(i,i)), &
                         naosh_ls,               &
                         naosh_ls,               &
                         nz,                     &
                         ipq_off_in,             &
                         scratch2(offset_2),     &
                         n1esh_f,                &
                         n1esh_f,                &
                         nz,                     &
                         ipq_off_in,             &
                         scratch3(offset_1),     &
                         nbas1_f,                &
                         n1esh_f,                &
                         nz_lowdin,              &
                         ipq_off_in,             &
                         scratch3(offset_1),     &
                         nbas1_f,                &
                         n1esh_f,                &
                         nz_lowdin,              &
                         ipq_off_in,             &
                         print_lvl)

!          debug print
           if(print_lvl > 2)then
             write(info_text,'(a27,i1)') 'x2c - 2c-fock-onMO - fsym =',i
             call print_x2cmat(scratch2(offset_2),n1esh_f,n1esh_f,nz,ipq_off_in,info_text,6)
           end if

         end if

!        set new offsets in matrices
         offset_1 = offset_1 + nbas1_f * n1esh_f
         offset_2 = offset_2 + n1esh_f * n1esh_f * nz
       end do


!      step 3: read original eigenvalues (to be kept) from DFCOEF
!      ----------------------------------------------------------
       l = 0
       do i = 1, nfsym
         l = l + norb_dim(i)
       end do
       call alloc(original_eigen_values,l)
       open(fh,file='DFCOEF',status='old',form='unformatted',access='sequential',position='rewind')
       call reacmo_no_work(fh,'DFCOEF',xmat,original_eigen_values,scratch5,dummy,4)
       close(fh,status='keep')
       offset_1 = 1
       offset_2 = 1
       do i = 1, nfsym
         call dcopy(nesh_dim(i),original_eigen_values(offset_1+norb_dim(i)-nesh_dim(i)),1,scratch4(offset_2),1) 
         call dcopy(nesh_dim(i),scratch4(offset_2),1,original_eigen_values(offset_2),1)
         offset_1 = offset_1 + norb_dim(i)
         offset_2 = offset_2 + nesh_dim(i)
       end do
       

!      step 4: initialize sub-block data (if linear_sym) required for block-diagonalization inside DFDIAG
!              reset common blocks with references to negative energy solutions deleted
!      --------------------------------------------------------------------------------------------------
       call renew_x2c_cb_orb_shell_dim(-137)
       
!      step 5: diagonalize the 2c-Fock matrix in orthonormal reference MO-basis
!      ------------------------------------------------------------------------
!          in: 2c-Fock matrix in onMO basis => scratch2
!
!         out: eigenvectors of the 2c-Fock matrix in onMO basis => xmat
!              eigenvalues of the 2c-Fock matrix                => scratch4
!              boson irreps/ <j_z>-blocking info                => scratch5

       offset_1 = 1
       call dfdiag(scratch2,scratch4,scratch5,xmat,.TRUE.,scratch3,offset_1,lwork)

!      step 6. save original eigenvalues, boson irreps and new MO coefficients in orthonormal reference MO basis on file
!      -----------------------------------------------------------------------------------------------------------------
       offset_1 = 1
       offset_2 = 1
       do i = 1, nfsym

!        a. eigenvalues
         write(flabel,'(a7,i4,i1)') 'eig-val',1,i
         call x2c_write(flabel,original_eigen_values(offset_1),nesh_dim(i),x2c_file_glb)
         if(print_lvl > 2)then
           write(info_text,'(a27,i1)') 'x2c - 2c-eigenvalues fsym =',i
           call print_x2cmat(scratch4(offset_1),nesh_dim(i),1,1,ipq_off_in,info_text,6)
           write(info_text,'(a27,i1)') 'x2c - orig. eigenval fsym =',i
           call print_x2cmat(original_eigen_values(offset_1),nesh_dim(i),1,1,ipq_off_in,info_text,6)
         end if

!        b. boson irreps / <j_z>-blocking information
         rewind x2c_file_glb
         write(flabel,'(a7,i4,i1)') 'bos-irr',1,i
         if(fndlab12('EOFLABEL-x2c',x2c_file_glb))then
           backspace(x2c_file_glb)
           call newlab12(flabel,x2c_file_glb,6)
           call writi(x2c_file_glb,nesh_dim(i),scratch5(offset_1))
           call newlab12('EOFLABEL-x2c',x2c_file_glb,6)
         end if

!        c. 2c-MO coefficients in MO reference basis
         write(flabel,'(a7,i4,i1)') 'cmo-2cM',1,i
         call x2c_write(flabel,xmat(offset_2),nesh_dim(i)*nesh_dim(i)*nz,x2c_file_glb)
!        debug print
         if(print_lvl > 2)then
           write(info_text,'(a27,i1)') 'x2c - 2c-cmo-onMO  - fsym =',i
           call print_x2cmat(xmat(offset_2),n1esh_f,n1esh_f,nz,ipq_off_in,info_text,6)
         end if

!        keep track of proper offsets
         offset_1  = offset_1 + nesh_dim(i)
         offset_2  = offset_2 + nesh_dim(i) * nesh_dim(i) * nz
       end do

       call dealloc(original_eigen_values)
       
  end subroutine x2c_diagonalize_2c_fock_operator_onmo
!**********************************************************************

  subroutine x2c_molecular_mean_field_mos(transformation_matrix_mo2ao,      &
                                          xmat,                             &
                                          mo_coeff_2c_ao,                   &
                                          eigen_values,                     &
                                          i2tmt,                            &
                                          i2tmot,                           &
                                          icmoq,                            &
                                          ntmo,                             &
                                          nfbas,                            &
                                          nesh_dim,                         &
                                          nz,                               &
                                          nzt,                              &
                                          nfsym,                            &
                                          nesh_dim_total,                   &
                                          ipq_off_in,                       &
                                          print_lvl)
!**********************************************************************
!
!    purpose: renew the DFCOEF file with the new picture-change-transformed MO
!             coefficients including the boson / <j_z>-blocking information 
!             and eigenvalues array.
!
!             note: the new DFCOEF still has references to the L+S AO-basis 
!             as we have not reset the AO-basis common blocks, only the MO-basis 
!             common blocks. 
!
!----------------------------------------------------------------------
     real(8), intent(inout)      :: transformation_matrix_mo2ao(*)
     real(8), intent(inout)      :: mo_coeff_2c_ao(*)
     real(8), intent(inout)      :: xmat(*)
     real(8), intent(inout)      :: eigen_values(*)
     integer, intent(in)         :: nfsym
     integer, intent(in)         :: i2tmt(nfsym)
     integer, intent(in)         :: i2tmot(nfsym)
     integer, intent(in)         :: icmoq(nfsym)
     integer, intent(in)         :: ntmo(nfsym)
     integer, intent(in)         :: nfbas(nfsym,0:2)
     integer, intent(in)         :: nesh_dim(nfsym)
     integer, intent(in)         :: nz
     integer, intent(in)         :: nzt
     integer, intent(in)         :: nesh_dim_total
     integer, intent(in)         :: ipq_off_in(4,0:7)
     integer, intent(in)         :: print_lvl
!----------------------------------------------------------------------
     real(8)                     :: total_energy
     integer                     :: i, j
     integer                     :: ndim_2c_mo_ao
     integer                     :: ndim_2c_ev
     integer                     :: ndim_2c_bj
     integer                     :: offset_1
     integer                     :: offset_2
     integer                     :: fh = 99
     character (len=12)          :: flabel
     character (len=24)          :: fdate
     character (len=74)          :: info_text
     logical                     :: fndlab12
     integer, allocatable        :: dimension_info(:,:)
     integer, allocatable        :: boson_irreps_jz_eigenvalues(:)
!----------------------------------------------------------------------

       call alloc(dimension_info,3,nfsym)
       call alloc(boson_irreps_jz_eigenvalues,nesh_dim_total)

!      read eigenvalues, boson irreps and new MO coefficients in orthonormal basis from file X2CMAT
       open(fh,file='X2CMAT',status='old',form='unformatted',access='sequential',action='read',position='rewind')

       offset_1 = 1
       offset_2 = 1
       do i = 1, nfsym

!        a. eigenvalues
         write(flabel,'(a7,i4,i1)') 'eig-val',1,i
         call x2c_read(flabel,eigen_values(offset_1),nesh_dim(i),fh)

!        b. boson irreps / <j_z>-blocking information
         rewind fh
         write(flabel,'(a7,i4,i1)') 'bos-irr',1,i
         if(fndlab12(flabel,fh))then
           call readi(fh,nesh_dim(i),boson_irreps_jz_eigenvalues(offset_1))
         end if

!        c. 2c-MO coefficients in MO reference basis
         write(flabel,'(a7,i4,i1)') 'cmo-2cM',1,i
         call x2c_read(flabel,xmat(offset_2),nesh_dim(i)*nesh_dim(i)*nz,fh)

!        debug print
         if(print_lvl > 2)then
           write(info_text,'(a27,i1)') 'x2c - 2c-cmo-onMO r- fsym =',i
           call print_x2cmat(xmat(offset_2),nesh_dim(i),nesh_dim(i),nz,ipq_off_in,info_text,6)
         end if

!        keep track of proper offsets
         offset_1            = offset_1 + nesh_dim(i)
         offset_2            = offset_2 + nesh_dim(i) * nesh_dim(i) * nz

!        write dimension info for the new DFCOEF
         dimension_info(1,i) = 0
         dimension_info(2,i) = nesh_dim(i)
         dimension_info(3,i) = nfbas(i,0)
       end do
       close(fh,status='keep')
!      set dimensions for the boson irreps / <j_z>-blocking information + eigenvalues arrays
       ndim_2c_bj = offset_1 - 1
       ndim_2c_ev = offset_1 - 1


!      read the title line and total energy (read code -1) from the old DFCOEF
       open(fh,file='DFCOEF',status='old',form='unformatted',access='sequential',position='rewind')
       call REACMO_no_work(fh,'DFCOEF',mo_coeff_2c_ao,eigen_values,boson_irreps_jz_eigenvalues,total_energy,-1)
       close(fh,status='keep')

!      backtransform the 2c-MO coefficients from the reference MO basis to AO basis
       ndim_2c_mo_ao = 0
       do i = 1, nfsym
         call bcktr1(mo_coeff_2c_ao(icmoq(i)+1),             &
                     nfbas(i,0),                             &
                     nesh_dim(i),                            &
                     xmat(i2tmot(i)+1),                      &
                     ntmo(i),                                &
                     ntmo(i),                                &
                     ntmo(i),                                &
                     nz,                                     &
                     ntmo(i),                                &
                     1,                                      &
                     nfbas(i,0),                             &
                     transformation_matrix_mo2ao(i2tmt(i)+1),&
                     nfbas(i,0),                             &
                     ntmo(i),                                &
                     nzt,                                    &
                     print_lvl)
         ndim_2c_mo_ao = ndim_2c_mo_ao + nfbas(i,0) * nesh_dim(i) * nz
       end do

!      write the new DFCOEF file
       open(fh,file='DFCOEF',status='replace',form='unformatted',access='sequential',position='rewind')
       info_text(1:50)  = ' X2Cmod - molecular-mean field 2c MOs              '
       info_text(51:74) = fdate()
       write(fh) info_text,nfsym,((dimension_info(i,j),i = 1,3),j=1,nfsym),total_energy
       write(fh) (mo_coeff_2c_ao(i),i=1,ndim_2c_mo_ao)
       write(fh) (eigen_values(i),i=1,ndim_2c_ev)
       write(fh) (boson_irreps_jz_eigenvalues(i),i=1,ndim_2c_bj)
       close(fh,status='keep')

       call dealloc(boson_irreps_jz_eigenvalues)
       call dealloc(dimension_info)

  end subroutine x2c_molecular_mean_field_mos

!**********************************************************************

end module x2c_pct_mo_coefficients
