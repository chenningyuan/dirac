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
                                                nr_cmo_q,            &
                                                naosh_ls,            &
                                                nr_ao_L,             &
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
     integer, intent(in)             :: nr_cmo_q
     integer, intent(in)             :: naosh_ls
     integer, intent(in)             :: nr_ao_L
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
     integer                         :: i, j, k
     integer                         :: is_fragment_tot
!**********************************************************************

!      initialize the fragment value - in the mmf-case it has (by definition) to be zero.
       is_fragment_tot = 0
       j               = 0
       k               = 0

!      step 0: set dimensions for allocation of scratch space + allocate
!      -----------------------------------------------------------------
       do i = 1, nfsym
         j = j + naosh_all(i)*naosh_L(i)*nz
         k = k + norb_dim(i)
       end do
       call alloc(pctmat,j)
       call alloc(mo_coeff_4c,nr_cmo_q)
       call alloc(mo_coeff_2c,nr_cmo_q)
       call alloc(eigen_values,k)
       call alloc(boson_irreps_jz_eigenvalues,k)
       pctmat                      = 0       
       mo_coeff_4c                 = 0
       mo_coeff_2c                 = 0
       eigen_values                = 0
       boson_irreps_jz_eigenvalues = 0

!      step 1: picture-change transform the MO-coefficients
!      ----------------------------------------------------
       call x2c_pctrafo_4cto2c_mo_coefficients(mo_coeff_file_name,           &
                                               eigen_values,                 &
                                               mo_coeff_2c,                  &
                                               mo_coeff_4c,                  &
                                               pctmat,                       &
                                               boson_irreps_jz_eigenvalues,  &
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
  subroutine x2c_pctrafo_4cto2c_mo_coefficients(mo_coeff_file_name,          &
                                                eigen_values,                &
                                                mo_coeff_2c,                 &
                                                mo_coeff_4c,                 &
                                                pctmat,                      &
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
     character (len=4)              :: numstr
     integer, allocatable           :: dimension_info(:,:)
!#define debug_mmf
#ifdef debug_mmf
     real(8), allocatable           :: tmp_buf(:)
     real(8), allocatable           :: tmp_buf2(:)
#endif
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
           call qgemm(nbas1_f,                                  &
                      n1esh_f,                                  &
                      nbasL_f,                                  &
                      val_d1,                                   &
                      'N','N',                                  &
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
           end if
         end if

!        set new offsets in matrices
         ioff_mo_coeff_2c =  ioff_mo_coeff_2c + nbas1_f*n1esh_f*nz
         ioff_mo_coeff_4c =  ioff_mo_coeff_4c + nbas1_f*n1esh_f*nz
         ioff_pct_mat     =  ioff_pct_mat     + nbas1_f*nbasL_f*nz
       end do

!      step 4: write 2c MO coefficients back to file
!      ---------------------------------------------
 
       call alloc(dimension_info,3,nfsym)
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
       print *, 'do not compress...'
!      4a. compress arrays to large-components
       ndim_2c_mo = ioff_mo_coeff_4c - 1
       ndim_2c_ev = 0
       ndim_2c_bj = 0
       do i = 1, nfsym
         dimension_info(1,i) = 0
         dimension_info(2,i) = nesh_dim(i)
         dimension_info(3,i) = naosh_all(i)
         ndim_2c_ev = ndim_2c_ev + nesh_dim(i) 
         ndim_2c_bj = ndim_2c_bj + nesh_dim(i)
       end do
       end if

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

#ifdef debug_mmf
       i = naosh_ls**2*nz*50
       allocate(tmp_buf(i))
       allocate(tmp_buf2(i))
       CALL NRMSTMO(mo_coeff_2c,mo_coeff_4c,tmp_buf,tmp_buf2,i,print_lvl,.FALSE.)
       deallocate(tmp_buf)
       deallocate(tmp_buf2)
#endif

#undef debug_mmf

  end subroutine x2c_pctrafo_4cto2c_mo_coefficients
!**********************************************************************

end module x2c_pct_mo_coefficients
