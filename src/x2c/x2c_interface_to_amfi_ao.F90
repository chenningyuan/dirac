!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
! module containing the driver routines for the addition of atomic-mean-field
! one- and two-electron spin-orbit corrections to the hamiltonian in DIRAC-sorted 
! SA-AO basis. 
!
! written by sknecht july 2010
!
module x2c_amfi_interface

  use memory_allocator
  use x2c_fio
  use x2c_utils, only:                 &
      pick_LL_block_saao_bas,          &
      print_x2cmat
  use x2c_cb_interface, only:          &
      sel_wrk_x2c,                     &
      desel_wrk_x2c

  implicit none

  public x2c_amfi_driver

  private

  real(8), parameter :: val_d1      =  1.0d0
  real(8), parameter :: val_d0      =  0.0d0
  real(8), parameter :: val_dm1     = -1.0d0

contains

!**********************************************************************
  subroutine x2c_amfi_driver(is_final_ham_lvl,        &
                             x2c_add_amfi,            &
                             naosh_ls,                &
                             nr_ao_L,                 &
                             naosh_all,               &
                             naosh_L,                 &
                             nfsym,                   &
                             nz,                      &
                             bs_irrep_mat,            &
                             ipq_off_in,              &
                             iqp_off_in,              &
                             iqmult_trip_q,           &
                             op_bs_to_fs,             &
                             quat_pointer,            &
                             quat_pointer_op,         &
                             max_quant_num,           &
                             nuclei_totnum,           &
                             tot_charge_for_mfsum,    &
                             amfi_order_contrib,      &
                             is_len_wrk,              &
                             spherical_on,            &
                             x2c_file_amf,            &
                             x2c_file_glb,            &
                             cspeed,                  &
                             print_lvl)
!**********************************************************************
!
!    purpose: driver routine pointing to the appropriate AMFI interface 
!             routine according to the value x2c_add_amfi.
!             The resulting hamiltonian h1+AMFI is stored on the global file 
!             X2CMAT and the "hamiltonian level" is set to its appropriate 
!             value (see read_1fock_x2c/dirx2c_utils.F90 for a detailed list 
!             of "hamiltonian levels").
!
!----------------------------------------------------------------------
     real(8), intent(in)     :: cspeed
     integer, intent(inout)  :: is_final_ham_lvl
     integer, intent(in)     :: x2c_add_amfi
     integer, intent(in)     :: naosh_ls
     integer, intent(in)     :: nr_ao_L
     integer, intent(in)     :: naosh_all(nfsym)
     integer, intent(in)     :: naosh_L(nfsym)
     integer, intent(in)     :: nfsym
     integer, intent(in)     :: nz
     integer, intent(in)     :: ipq_off_in(4,0:7)
     integer, intent(in)     :: iqp_off_in(4,0:7)
     integer, intent(in)     :: op_bs_to_fs(0:7,1:2)
     integer, intent(in)     :: bs_irrep_mat(4,0:7)
     integer, intent(in)     :: iqmult_trip_q(4,4,4)
     integer, intent(in)     :: quat_pointer(0:7,2)
     integer, intent(in)     :: quat_pointer_op(0:7)
     integer, intent(in)     :: x2c_file_amf
     integer, intent(in)     :: x2c_file_glb
     integer, intent(in)     :: is_len_wrk
     integer, intent(in)     :: max_quant_num
     integer, intent(in)     :: nuclei_totnum
     integer, intent(in)     :: tot_charge_for_mfsum
     integer, intent(in)     :: amfi_order_contrib
     integer, intent(in)     :: print_lvl
     logical, intent(in)     :: spherical_on
!----------------------------------------------------------------------
!**********************************************************************
 
       if(x2c_add_amfi == 1)then
!        interface to "old" AMFI (by B. Schimmelpfennig)
         call x2c_interface2old_amfi(naosh_ls,             &
                                     nr_ao_L,              &
                                     naosh_all,            &
                                     naosh_L,              &
                                     nfsym,                &
                                     nz,                   &
                                     bs_irrep_mat,         &
                                     ipq_off_in,           &
                                     iqp_off_in,           &
                                     iqmult_trip_q,        &
                                     op_bs_to_fs,          &
                                     quat_pointer,         &
                                     quat_pointer_op,      &
                                     cspeed,               &
                                     max_quant_num,        &
                                     nuclei_totnum,        &
                                     tot_charge_for_mfsum, &
                                     amfi_order_contrib,   &
                                     is_len_wrk,           &
                                     spherical_on,         &
                                     x2c_file_amf,         &
                                     x2c_file_glb,         &
                                     print_lvl)
         is_final_ham_lvl = -3
       else if(x2c_add_amfi == 2)then
!        interface to new AMFI
!        call x2c_interface2new_amfi
!        is_final_ham_lvl = TODO
       else if(x2c_add_amfi == 3)then
!        interface to new AMFI for atom-by-atom X2C
!        call x2c_interface2new_amfi_abax2c
!        is_final_ham_lvl = TODO
       end if

  end subroutine x2c_amfi_driver

!**********************************************************************
  subroutine x2c_interface2old_amfi(naosh_ls,              &
                                    nr_ao_L,               &
                                    naosh_all,             &
                                    naosh_L,               &
                                    nfsym,                 &
                                    nz,                    &
                                    bs_irrep_mat,          &
                                    ipq_off_in,            &
                                    iqp_off_in,            &
                                    iqmult_trip_q,         &
                                    op_bs_to_fs,           &
                                    quat_pointer,          &
                                    quat_pointer_op,       &
                                    cspeed,                &
                                    max_quant_num,         &
                                    nuclei_totnum,         &
                                    tot_charge_for_mfsum,  &
                                    amfi_order_contrib,    &
                                    is_len_wrk,            &
                                    spherical_on,          &
                                    x2c_file_amf,          &
                                    x2c_file_glb,          &
                                    print_lvl)   
!**********************************************************************
!
!    purpose: interface routine to the old AMFI program by B. Schimmelpfennig. 
!             According to the value of "amfi_order_contrib" two-electron
!             spin-orbit contributions are added to the X2C Hamiltonian.
!
!             the outline of this routine may be sketched as follows:
!
!             1. read h1_++_saao_cartesian_gaussian from file X2Camfi_scr (in ++ form but 4c-format retained)
!            /  a. insert quaternion phases
!   "see"   /   b. transform h1 from Dirac-sorted AO (SA-AO) to Hermit-sorted SO-AO basis
!  -------  \   c. transform h1 from Hermit-sorted SO-AO basis to simple AO basis
!            \  d. transform h1 from cartesian-gaussian AO basis to spherical-gaussian AO basis
!             2. add AMFI contributions
!            /  e. transform h1 from spherical-gaussian AO basis to cartesian-gaussian AO basis
!   "saw"   /   f. transform h1 from simple AO basis to Hermit-sorted SO-AO basis
!  -------  \   g. transform h1 from Hermit-sorted SO-AO basis to Dirac-sorted AO (SA-AO) basis
!            \  h. insert quaternion phases
!             3. pick out LL block
!             4. save h1_++_LL_saao_cartesian_gaussian to file X2CMAT (only LL-block retained)

!    note: sphgen modifies the common blocks for 
!          cartesian --> spherical gaussian transformations;
!          a call should therefore be followed by a reset call 
!          or rewrite the code such that it works without these common blocks
!          (taken from the notes by Andre Gomes and Luuk Visscher).
!
!----------------------------------------------------------------------
     real(8), intent(in)             :: cspeed
     integer, intent(in)             :: naosh_ls
     integer, intent(in)             :: nr_ao_L
     integer, intent(in)             :: naosh_all(nfsym)
     integer, intent(in)             :: naosh_L(nfsym)
     integer, intent(in)             :: nfsym
     integer, intent(in)             :: nz
     integer, intent(in)             :: ipq_off_in(4,0:7)
     integer, intent(in)             :: iqp_off_in(4,0:7)
     integer, intent(in)             :: op_bs_to_fs(0:7,1:2)
     integer, intent(in)             :: bs_irrep_mat(4,0:7)
     integer, intent(in)             :: iqmult_trip_q(4,4,4)
     integer, intent(in)             :: quat_pointer(0:7,2)
     integer, intent(in)             :: quat_pointer_op(0:7)
     integer, intent(in)             :: x2c_file_amf
     integer, intent(in)             :: x2c_file_glb
     integer, intent(in)             :: is_len_wrk
     integer, intent(in)             :: max_quant_num
     integer, intent(in)             :: nuclei_totnum
     integer, intent(in)             :: tot_charge_for_mfsum
     integer, intent(in)             :: amfi_order_contrib
     integer, intent(in)             :: print_lvl
     logical, intent(in)             :: spherical_on
!--------------------------------------------------------------
     real(8), pointer                :: wrk(:)
     real(8), allocatable            :: scr1_mat(:)
     real(8), allocatable            :: scr2_mat(:)
     real(8), allocatable            :: rcha_pnuc(:)
     integer, allocatable            :: itmp_mat1(:)
     integer, allocatable            :: itmp_mat2(:)
     integer, allocatable            :: itmp_mat3(:)
     integer, allocatable            :: itmp_mat4(:)
     integer, allocatable            :: itmp_mat5(:)
     integer, allocatable            :: itmp_mat6(:)
     integer, allocatable            :: itmp_mat7(:)
     integer, allocatable            :: itmp_mat8(:)
     integer, allocatable            :: itmp_mat9(:)
     character (len=4), allocatable  :: clab_int(:)
     integer                         :: sph_fun_L
     integer                         :: nr_sph
     integer                         :: i
     integer                         :: j
     integer                         :: iq
     integer                         :: ipq
     integer                         :: iz
     integer                         :: irepd
     integer                         :: nsphcm
     integer                         :: wf_components
     integer                         :: lwrk
     character (len=12)              :: flabel
!**********************************************************************
 
!      allocate first scratch matrix
       call alloc(scr1_mat,naosh_ls*naosh_ls*4)
       scr1_mat = 0

!      -------------------------------------------
!      step 1: read h1 in Dirac-sorted SA-AO basis
!      -------------------------------------------

       i = 0
       j = 0
       write(flabel,'(a7,i4,i1)') 'h12cAOa',j,i
       call x2c_read(flabel,scr1_mat,naosh_ls*naosh_ls*nz,x2c_file_amf)

!      a. insert quaternion phase factors
!      ----------------------------------
       if(nz < 4)then
         do iz = 1, nz                                                
           iq = ipq_off_in(iz,0)
           call q2bphase('F',iq,1,scr1_mat(1+(naosh_ls*naosh_ls*(iz-1))))
         end do           
       end if          

!      b. transform hamiltonian integral elements from DIRAC-sorted AO to Hermit-sorted AO basis  
!      -----------------------------------------------------------------------------------------
       call bstobu_no_work(scr1_mat,nz)

!      allocate second scratch matrix
       call alloc(scr2_mat,naosh_ls*naosh_ls*4)
       scr2_mat = 0

!      c. transform from SO-AO basis to AO-basis
!      -----------------------------------------
       i = op_bs_to_fs(0,1)
       do iz = 1, 4
         irepd = bs_irrep_mat(iz,0)
         iq    = iqmult_trip_q(1,quat_pointer(irepd,i),iz)
         ipq   = iqp_off_in(iq,0)
         call mtsoao(scr1_mat(1+(naosh_ls*naosh_ls*(ipq-1))),    &
                     scr2_mat(1+(naosh_ls*naosh_ls*(iz-1))),     &
                     naosh_ls,                                   &
                     irepd,                                      &
                     print_lvl)
       end do

!      release scratch space
       call dealloc(scr1_mat)

       lwrk = is_len_wrk
       call sel_wrk_x2c(wrk,lwrk)

       call sphgen(1,2,.false.,wrk,lwrk,print_lvl)

!      re-allocate scratch matrix
       call alloc(scr1_mat,naosh_ls*naosh_ls*4)
       scr1_mat = 0

!      d. transform from cartesian Gaussian AO basis to spherical Gaussian basis
!      -------------------------------------------------------------------------
!        in : h1_cartesian -> scr2_mat; wf_components == 0: transform large+small component
!        out: h1_spherical -> scr1_mat
       wf_components = 0 
       do iz = 1, 4
         call mtaocs(scr2_mat(1+(naosh_ls*naosh_ls*(iz-1))),     &
                     naosh_ls,                                   &
                     scr1_mat(1+(naosh_ls*naosh_ls*(iz-1))),     &
                     naosh_ls,                                   &
                     wf_components,                              &
                     wrk,                                        &
                     lwrk,                                       &
                     print_lvl)
       end do

!      determine the total number of spherical functions for the large (L => 1) component
       sph_fun_L = nsphcm(1)

       allocate(clab_int(max_quant_num**2))

!      get proper labels (including signs) for s,p,d,f,... functions
       call sphlab(max_quant_num-1,clab_int)

!      release scratch space
       call dealloc(scr2_mat)

!      allocate temporary matrices for AMFI module
       call alloc(itmp_mat1,sph_fun_L)
       call alloc(itmp_mat2,sph_fun_L)
       call alloc(itmp_mat3,sph_fun_L)
       call alloc(itmp_mat4,sph_fun_L)
       call alloc(itmp_mat5,sph_fun_L)
       call alloc(itmp_mat6,sph_fun_L)
       call alloc(itmp_mat7,sph_fun_L*2)
       call alloc(itmp_mat8,nuclei_totnum)
       call alloc(itmp_mat9,nuclei_totnum)
       call alloc(rcha_pnuc,nuclei_totnum)

!      ------------------------------
!      step 2. add AMFI contributions
!      ------------------------------
!        in : h1_spherical        -> scr1_mat
!        out: h1_spherical + AMFI -> scr1_mat

       call amfiin(scr1_mat,               &
                   naosh_ls,               &
                   itmp_mat2,              &
                   itmp_mat1,              &
                   itmp_mat3,              &
                   itmp_mat8,              &
                   itmp_mat4,              &
                   itmp_mat5,              &
                   itmp_mat6,              &
                   itmp_mat9,              &
                   itmp_mat7,              &
                   rcha_pnuc,              &
                   sph_fun_L,              &
                   clab_int,               &
                   .false.,                &
                   .false.,                &
                   tot_charge_for_mfsum,   &
                   amfi_order_contrib,     &
                   quat_pointer_op,        &
                   cspeed,                 &
                   print_lvl,              &
                   wrk,                    &
                   lwrk)


!      release scratch memory
       deallocate(clab_int)
       call dealloc(itmp_mat1)
       call dealloc(itmp_mat2)
       call dealloc(itmp_mat3)
       call dealloc(itmp_mat4)
       call dealloc(itmp_mat5)
       call dealloc(itmp_mat6)
       call dealloc(itmp_mat7)
       call dealloc(itmp_mat8)
       call dealloc(itmp_mat9)
       call dealloc(rcha_pnuc)


!      set up the proper factors for the backtransformation
       call sphgen(1,2,.true.,wrk,lwrk,print_lvl)

!      re-allocate scratch space
       call alloc(scr2_mat,naosh_ls*naosh_ls*4)
       scr2_mat = 0

!      e. transform from spherical Gaussian AO basis to cartesian Gaussian AO basis
!      ----------------------------------------------------------------------------
!        in : h1_spherical -> scr1_mat
!        out: h1_cartesian -> scr2_mat
       do iz = 1, 4
         call mtaosc(scr2_mat(1+(naosh_ls*naosh_ls*(iz-1))),     &
                     naosh_ls,                                   &
     &               scr1_mat(1+(naosh_ls*naosh_ls*(iz-1))),     &
                     naosh_ls,                                   &
                     wf_components,                              &
                     wrk,                                        &
                     lwrk,                                       &
                     print_lvl)
       end do

!      reset the common blocks (see note in the header)
       call sphgen(1,2,.false.,wrk,lwrk,print_lvl)

       call desel_wrk_x2c(wrk,lwrk)

!      f. transform from AO-basis to SO-AO basis
!      -----------------------------------------
!        in : h1_AO    -> scr2_mat
!        out: h1_SO-AO -> scr1_mat
       scr1_mat = 0
       i = op_bs_to_fs(0,1)
       do iz = 1, 4
         irepd = bs_irrep_mat(iz,0)
         iq    = iqmult_trip_q(1,quat_pointer(irepd,i),iz)
         ipq   = iqp_off_in(iq,0)
         call mtaoso(scr2_mat(1+(naosh_ls*naosh_ls*(iz-1))),     &
                     scr1_mat(1+(naosh_ls*naosh_ls*(ipq-1))),    &
                     naosh_ls,                                   &
                     irepd,                                      &
                     print_lvl)
       end do

!      release scratch space
       call dealloc(scr2_mat)

!      g. transform from Hermit-sorted AO to DIRAC-sorted AO basis
!      -----------------------------------------------------------
       call butobs_no_work(scr1_mat,nz)

!      h. insert quaternion phase factors
!      ----------------------------------
       if(nz < 4)then
         do iz = 1, nz                                                
           iq = ipq_off_in(iz,0)
           call q2bphase('F',iq,1,scr1_mat(1+(naosh_ls*naosh_ls*(iz-1))))
         end do           
       end if          

!      re-allocate scratch space
       call alloc(scr2_mat,naosh_ls*naosh_ls*nz)
       scr2_mat = 0

!      --------------------------------------------------------
!      step 3. select LL-block of h1_++_saao_cartesian_gaussian
!      --------------------------------------------------------

       call pick_LL_block_saao_bas(scr1_mat,            &
                                   scr2_mat,            &
                                   naosh_L,             &
                                   naosh_all,           &
                                   naosh_ls,            &
                                   nr_ao_L,             &
                                   nfsym,               &
                                   nz,                  &
                                   0,                   &
                                   op_bs_to_fs)

!      -------------------------------------------------------------
!      step 4. write h1_++_LL_saao_cartesian_gaussian to file X2CMAT 
!      -------------------------------------------------------------

       i = 0
       j = 0
       write(flabel,'(a7,i4,i1)') 'h12cAOA',j,i

       if(spherical_on)then

         lwrk = is_len_wrk
         call sel_wrk_x2c(wrk,lwrk)
         call sph_iotc(scr2_mat,scr1_mat,nr_sph,print_lvl,wrk,lwrk)
         call desel_wrk_x2c(wrk,lwrk)

         call x2c_write(flabel,scr1_mat,nr_sph**2  * nz,x2c_file_glb)
       else

!        debug print
         if(print_lvl > 2)then
           call print_x2cmat(scr2_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - h12c+AMFI',6)
         end if

         call x2c_write(flabel,scr2_mat,nr_ao_L**2 * nz,x2c_file_glb)
       end if

!      release scratch space
       call dealloc(scr1_mat)
       call dealloc(scr2_mat)

  end subroutine x2c_interface2old_amfi

!**********************************************************************

end module x2c_amfi_interface
