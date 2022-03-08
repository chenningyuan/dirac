#define DEBUG_SOC
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
! module containing the driver routines for the addition of atomic-mean-field
! two-electron spin-orbit corrections to the hamiltonian in DIRAC-sorted SA-AO basis. 
!
! written by sknecht july 2010
!
module x2c_2e_soc_interface

  use x2c_fio
  use x2c_utils, only:                 &
      print_x2cmat
  use picture_change_operations
#ifdef MOD_XAMFI
  use xamfi_driver
  use x_fock, only: get_x_fock
  use gab_settings, only: recompute_Gab
#endif
  use x2cmod_cfg, only: x2c_intgen_2c
  use xcint_main, only: set_xc_energy_delta

  implicit none

#include "cbirea.h"
#include "dcbgen.h"
#include "dcbcls.h"

  public x2c_add_2e_so_corrections
  public x2c_add_2e_diagonal_atomic
  public x2c_add_2e_from_atomic_dens

  private

contains

!**********************************************************************
  subroutine x2c_add_2e_so_corrections(is_final_ham_lvl,        &
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
                                       nuclei_type,             &
                                       tot_charge_for_mfsum,    &
                                       amfi_order_contrib,      &
                                       is_len_wrk,              &
                                       spherical_on,            &
                                       x2c_file_amf,            &
                                       x2c_file_glb,            &
                                       cspeed,                  &
                                       ioff_aomat_x,            &
                                       isDFT,                   &
                                       print_lvl)
!**********************************************************************
!
!    purpose: driver routine pointing to the appropriate Dirac-2e-SOC/external 
!             AMFI interface routine according to the value x2c_add_amfi.
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
     integer, intent(in)     :: ioff_aomat_x(nfsym,nfsym)
     integer, intent(in)     :: x2c_file_amf
     integer, intent(in)     :: x2c_file_glb
     integer, intent(in)     :: is_len_wrk
     integer, intent(in)     :: max_quant_num
     integer, intent(in)     :: nuclei_totnum
     integer, intent(in)     :: nuclei_type
     integer, intent(in)     :: tot_charge_for_mfsum
     integer, intent(in)     :: amfi_order_contrib
     integer, intent(in)     :: print_lvl
     logical, intent(in)     :: spherical_on
     logical, intent(in)     :: isDFT
!----------------------------------------------------------------------
     integer                 :: i
!**********************************************************************
 
       if(x2c_add_amfi == 1)then
!        interface to "old" AMFI (AMFI was written by the late B. Schimmelpfennig)
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
       else
#ifdef MOD_XAMFI
         call x2c_interface2aoosoc(naosh_ls,             &
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
                                   cspeed,               &
                                   max_quant_num,        &
                                   nuclei_totnum,        &
                                   nuclei_type,          &
                                   tot_charge_for_mfsum, &
                                   amfi_order_contrib,   &
                                   is_len_wrk,           &
                                   spherical_on,         &
                                   x2c_file_amf,         &
                                   x2c_file_glb,         &
                                   (x2c_add_amfi == 3),  &
                                   ioff_aomat_x,         &
                                   isDFT,                &
                                   print_lvl)
#else
         call quit('unknown 2e-SOC interface in x2c module')
#endif
       end if

       is_final_ham_lvl = -3

  end subroutine x2c_add_2e_so_corrections

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
!             1. prepare for spherical harmonics
!             2. calculate AMFI contributions
!            /  a. transform from spherical-gaussian AO basis to cartesian-gaussian AO basis
!   "saw"   /   b. transform from simple AO basis to Hermit-sorted SO-AO basis
!  -------  \   c. transform from Hermit-sorted SO-AO basis to Dirac-sorted AO (SA-AO) basis
!            \  d. insert quaternion phases
!             3. read h1_++_saao_cartesian_gaussian from file X2Camfi_scr (in ++ form but 4c-format retained)
!             4. pick out LL block
!             5. save h1_++_LL_saao_cartesian_gaussian to file X2CMAT (only LL-block retained)

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
     real(8), allocatable            :: wrk(:)
     real(8), allocatable            :: scr1_mat(:)
     real(8), allocatable            :: scr2_mat(:)
     real(8), allocatable            :: scr3_mat(:)
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
 
!      ---------------------------------------
!      step 1. prepare for spherical harmonics
!      ---------------------------------------

       lwrk = is_len_wrk
       allocate(wrk(lwrk))

       call sphgen(1,2,.false.,wrk,lwrk,print_lvl)

!      determine the total number of spherical functions for the large (L => 1) component
       sph_fun_L = nsphcm(1)

       allocate(clab_int(max_quant_num**2))

!      get proper labels (including signs) for s,p,d,f,... functions
       call sphlab(max_quant_num-1,clab_int)


!      allocate temporary matrices for AMFI module
       allocate(itmp_mat1(sph_fun_L)    )
       allocate(itmp_mat2(sph_fun_L)    )
       allocate(itmp_mat3(sph_fun_L)    )
       allocate(itmp_mat4(sph_fun_L)    )
       allocate(itmp_mat5(sph_fun_L)    )
       allocate(itmp_mat6(sph_fun_L)    )
       allocate(itmp_mat7(sph_fun_L*2)  )
       allocate(itmp_mat8(nuclei_totnum))
       allocate(itmp_mat9(nuclei_totnum))
       allocate(rcha_pnuc(nuclei_totnum))

!      ------------------------------
!      step 2. add AMFI contributions
!      ------------------------------
!        out: AMFI (spherical) -> scr1_mat

!      allocate first scratch matrix
       allocate(scr1_mat(naosh_ls*naosh_ls*4))
       scr1_mat = 0.0d0

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
       deallocate(itmp_mat1)
       deallocate(itmp_mat2)
       deallocate(itmp_mat3)
       deallocate(itmp_mat4)
       deallocate(itmp_mat5)
       deallocate(itmp_mat6)
       deallocate(itmp_mat7)
       deallocate(itmp_mat8)
       deallocate(itmp_mat9)
       deallocate(rcha_pnuc)


!      set up the proper factors for the backtransformation
       call sphgen(1,2,.true.,wrk,lwrk,print_lvl)

!      allocate scratch space
       allocate(scr2_mat(naosh_ls*naosh_ls*4),scr3_mat(naosh_ls*naosh_ls*4))
       scr2_mat = 0; scr3_mat = 0

!      a. transform from spherical Gaussian AO basis to cartesian Gaussian AO basis
!      ----------------------------------------------------------------------------
!        in : AMFI contributions -> scr1_mat
!        out: AMFI contributions -> scr2_mat
!             wf_components = 0: large+small
       wf_components = 0 
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

       deallocate(wrk)

!      b. transform from AO-basis to SO-AO basis
!      -----------------------------------------
!        in : h1_AO    -> scr2_mat
!        out: h1_SO-AO -> scr1_mat
       scr1_mat = 0.0d0
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

!      c. transform from Hermit-sorted AO to DIRAC-sorted AO basis
!      -----------------------------------------------------------
       call butobs_no_work(scr1_mat,nz)

!      d. insert quaternion phase factors
!      ----------------------------------
       if(nz < 4)then
         do iz = 1, nz                                                
           iq = ipq_off_in(iz,0)
           call q2bphase('F',iq,1,scr1_mat(1+(naosh_ls*naosh_ls*(iz-1))))
         end do           
       end if          

       call pick_LL_block_saao_bas(scr1_mat,            &
                                   scr3_mat,            &
                                   naosh_L,             &
                                   naosh_all,           &
                                   naosh_ls,            &
                                   nr_ao_L,             &
                                   nfsym,               &
                                   nz,                  &
                                   0,                   &
                                   op_bs_to_fs)

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr3_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - AMFI contribution',6)
       end if

#ifdef DEBUG_SOC
       open(99,file='soc-contributions',status='replace',form='formatted',  &
       access='sequential',action="readwrite",position='rewind')

       call print_x2cmat(scr1_mat,naosh_ls,naosh_ls,nz,ipq_off_in,'x2c - AMFI contribution',99)
       close(99,status='keep')
#endif

!      -------------------------------------------------------------------------------
!      step 3. read h1_2c_saao_cartesian_gaussian from file and add AMFI contributions 
!      -------------------------------------------------------------------------------

       i = 0
       j = 0
       write(flabel,'(a7,i4,i1)') 'h12cAOa',1,i
       scr2_mat = 0.0d0
       call x2c_read(flabel,scr2_mat,naosh_ls*naosh_ls*nz,x2c_file_amf)

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr2_mat,naosh_ls,naosh_ls,nz,ipq_off_in,'x2c - h1 wo AMFI contribution',6)
       end if


!      add amfi contributions to h1_2c
       call daxpy(naosh_ls*naosh_ls*nz,1.0d0,scr2_mat,1,scr1_mat,1)

       scr2_mat = 0.0d0

!      --------------------------------------------------------
!      step 4. select LL-block of h1_++_saao_cartesian_gaussian
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
!      step 5. write h1_++_LL_saao_cartesian_gaussian to file X2CMAT 
!      -------------------------------------------------------------

       i = 0
       j = 0
       write(flabel,'(a7,i4,i1)') 'h12cAOA',1,i

       if(spherical_on)then

         lwrk = is_len_wrk
         allocate(wrk(lwrk))
         call sph_iotc(scr2_mat,scr1_mat,nr_sph,print_lvl,wrk,lwrk)
         deallocate(wrk)

         call x2c_write(flabel,scr1_mat,nr_sph**2  * nz,x2c_file_glb)
       else

!        debug print
         if(print_lvl > 2)then
           call print_x2cmat(scr2_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - h12c+AMFI',6)
         end if

         call x2c_write(flabel,scr2_mat,nr_ao_L**2 * nz,x2c_file_glb)
       end if

!      release scratch space
       deallocate(scr1_mat,scr2_mat,scr3_mat)

  end subroutine x2c_interface2old_amfi

!**********************************************************************

  subroutine x2c_interface2aoosoc(naosh_ls,              &
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
                                  cspeed,                &
                                  max_quant_num,         &
                                  nuclei_totnum,         &
                                  nuclei_type,           &
                                  tot_charge_for_mfsum,  &
                                  amfi_order_contrib,    &
                                  is_len_wrk,            &
                                  spherical_on,          &
                                  x2c_file_amf,          &
                                  x2c_file_glb,          &
                                  DOeamf,                &
                                  ioff_aomat_x,          &
                                  isDFT,                 &
                                  print_lvl)   
!**********************************************************************
!
!    purpose: interface to the (extended) atomic-mean-field picture-change error-correction (PCE)
!             module. The PCE corrections comprise two-electron scalar- as well as spin-orbit  
!             contributions.
!             This work is published in XXX (2022).
!
!              amf: atomic mean-field (Algorithm 1 in XXX)
!             eamf: extended atomic mean-field (Algorithm 2 in XXX)
!
!             the outline of this routine may be sketched as follows:
!      IF
!             eamf:
!             =====
!
!      ELSE
!              amf:
!             =====
!             for ATOM in MOLECULE
!             do
!                    retrieve amf corrections from file (assumed to be in C1*)
!                    add atomic block in molecular matrix 
!             done
!                                                                           
!             transform the molecular matrix (C1* symmetry) to actual symmetry
!   "saw"   /   a. transform from simple AO basis to Hermit-sorted SO-AO basis
!  -------  \   b. transform from Hermit-sorted SO-AO basis to Dirac-sorted AO (SA-AO) basis
!            \  c. insert quaternion phases
!
!      END IF
!             
!             common:
!             =======

!             1. read h1_++_saao_cartesian_gaussian from file X2Camfi_scr (in ++ form but the 4c-format was kept)
!             2. pick out ++ block
!             3. add (e)amf corrections to h1 and save the final h1_++_LL_saao_cartesian_gaussian to file X2CMAT 
!                (only the ++ block will be kept)
!
!             rewritten by S. Knecht, Jan 2022 (first draft from Sep 2014?)
!
!----------------------------------------------------------------------
     real*8 , intent(in)             :: cspeed
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
     integer, intent(in)             :: ioff_aomat_x(nfsym,nfsym)
     integer, intent(in)             :: x2c_file_amf
     integer, intent(in)             :: x2c_file_glb
     integer, intent(in)             :: is_len_wrk
     integer, intent(in)             :: max_quant_num
     integer, intent(in)             :: nuclei_totnum
     integer, intent(in)             :: nuclei_type
     integer, intent(in)             :: tot_charge_for_mfsum
     integer, intent(in)             :: amfi_order_contrib
     integer, intent(in)             :: print_lvl
     logical, intent(in)             :: spherical_on
     logical, intent(in)             :: DOeamf
     logical, intent(in)             :: isDFT
!--------------------------------------------------------------
     real*8 , allocatable            :: scr1_mat(:)
     real*8 , allocatable            :: scr2_mat(:)
     real*8 , allocatable            :: scr3_mat(:)
     real*8 , allocatable            :: scr4_mat(:)
     real*8 , allocatable            :: wrk(:)
     real*8                          :: delta_xc_e
     integer, allocatable            :: aooso_info(:)
     integer                         :: i
     integer                         :: j
     integer                         :: iq
     integer                         :: ipq
     integer                         :: iz
     integer                         :: irepd
     integer                         :: nsphcm
     integer                         :: nr_sph
     integer                         :: nrows_frag
     integer                         :: ncols_frag
     integer                         :: center_total
     integer                         :: lwrk
     integer                         :: wf_components
     integer                         :: ioff_U
     integer                         :: op_fer_rep
     integer                         :: integral_types
     INTEGER, PARAMETER              :: nzc1 = 4
     character (len=12)              :: flabel
!**********************************************************************
 
!      allocate scratch space
       allocate(scr1_mat(naosh_ls*naosh_ls*nzc1),scr2_mat(naosh_ls*naosh_ls*nzc1),&
                scr3_mat(naosh_ls*naosh_ls*nzc1))
       scr1_mat = 0; scr2_mat = 0; scr3_mat = 0

       lwrk = is_len_wrk;

       if(.not.DOeamf)then ! amf

!             ----------------------------------------------------------------------
!             read and distribute amf contributions for symmetry independent centers
!             ----------------------------------------------------------------------
!             out: atomic-wise amf corrections in molecular matrix -> scr1_mat

              allocate(aooso_info(naosh_ls))

              center_total = 1
              scr1_mat     = 0

              !> a. read contributions for symmetry independent centers
              do i = 1, nuclei_type

!             count the number of basis functions for center_total
              call labcount(ncols_frag,aooso_info,nr_ao_L,1,-1,center_total,1,-1)
              nrows_frag = ncols_frag

#ifdef MOD_XAMFI
              !> read X-AMFI contributions (assumed to be exported in C1 symmetry)
              call put_xamfi_correction(                &
                                        scr1_mat,       &
                                        nr_ao_L,        &
                                        nrows_frag,     &
                                        ncols_frag,     &
                                        nzc1,           &
                                        i,              &
                                        center_total,   &
                                        print_lvl       &
                                       )
#else
              call quit('X-AMFI not avaliable in this version.')
#endif
              end do

              !> set offsets in atomic basis to two-component mode
              IPREAD  = IPREAD - 10000
              TWOCOMP = .TRUE.
              RDINPC  = .FALSE.
              CALL READIN(.FALSE.)
              CALL SETDC1(IPREAD)
              call getlab(ipread)

!             b. transform from AO-basis (C1* symm sorted) to SO-AO basis (symmetry-ordered)
!             ------------------------------------------------------------------------------
!             in : xamfi_AO    -> scr1_mat
!             out: xamfi-SO-AO -> scr3_mat
              i = op_bs_to_fs(0,1)
              do iz = 1, nzc1
              irepd = bs_irrep_mat(iz,0)
              iq    = iqmult_trip_q(1,quat_pointer(irepd,i),iz)
              ipq   = iqp_off_in(iq,0)
              call mtaoso(scr1_mat(1+(nr_ao_L*nr_ao_L*(iz-1))),     &
                          scr3_mat(1+(nr_ao_L*nr_ao_L*(ipq-1))),    &
                          nr_ao_L,                                  &
                          irepd,                                    &
                          print_lvl)
              end do

!             c. transform from Hermit-sorted AO to DIRAC-sorted AO basis
!             -----------------------------------------------------------
              call butobs_no_work(scr3_mat,nz)

!             d. insert quaternion phase factors
!             ----------------------------------
              if(nz < 4)then
                 do iz = 1, nz
                     iq = ipq_off_in(iz,0)
                     call q2bphase('F',iq,1,scr3_mat(1+(nr_ao_L*nr_ao_L*(iz-1))))
                 end do
              end if

              !> restore offsets in atomic basis to four-component mode (for now)
              TWOCOMP = .FALSE.
              CALL READIN(.FALSE.)
              CALL SETDC1(IPREAD)
              call getlab(ipread)

              IPREAD = IPREAD + 10000

!             debug print
              if(print_lvl > 2)then
                call print_x2cmat(scr3_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - amf contributions',6)
              end if
              deallocate(aooso_info)

       else ! eamf

              !> our correction term is \Delta F^{2c,2e} = { U+ F^{4c,2e}[D^{4c}_(+)atoms] U }^{++} -  F^{2c,2e}[D^{2c}_(+)atoms]

              !> to save memory we start with the second term: F^{2c,2e}[D^{2c}_(+)atoms]
!                                                              --------------------------
              !> set offsets in atomic basis to two-component mode
              IPREAD  = IPREAD - 10000
              TWOCOMP = .TRUE.; TWOCOMPBSS = .TRUE.
              RDINPC  = .FALSE.
              CALL READIN(.FALSE.)
              CALL SETDC1(IPREAD)
              call getlab(ipread)

              allocate(scr4_mat(nr_ao_L*nr_ao_L*nz)); scr4_mat = 0
              integral_types = x2c_intgen_2c
#ifdef MOD_XAMFI
              !> since we calculated a Fock matrix in 4c mode to compute the defining X2C Hamiltonian
              !> we need to enforce the recomputation of the GAB matrix in 2c mode
              recompute_Gab = .true.;

              ! returns F^{2c,2e}[D^{2c}_(+)atoms] in scr4_mat
              call x2c_add_2e_from_atomic_dens(scr4_mat,                          &
                                               nr_ao_L,                           &
                                               nz,                                &
                                               bs_irrep_mat,                      &
                                               ipq_off_in,                        &
                                               iqp_off_in,                        &
                                               iqmult_trip_q,                     &
                                               op_bs_to_fs,                       &
                                               quat_pointer,                      &
                                               nuclei_totnum,                     &
                                               nuclei_type,                       &
                                               2,                                 & ! Xc-mode
                                               integral_types,                    & ! integral flag
                                               lwrk,                              & ! lwrk
                                               isDFT,                             & ! add XC potential (DFT)
                                               print_lvl)
              !> ... and reset (we will from here on continue in 2c mode)
              recompute_Gab = .false.;
#else
              call quit('X-AMFI not avaliable in this version.')
#endif

!             debug print
              if(print_lvl > 2)then
                call print_x2cmat(scr4_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - eamf term 2 contributions',6)
              end if

              !> restore offsets in atomic basis to four-component mode (for now)
              TWOCOMP = .FALSE.; TWOCOMPBSS = .FALSE.
              CALL READIN(.FALSE.)
              CALL SETDC1(IPREAD)
              call getlab(ipread)

              IPREAD = IPREAD + 10000

              !> first term: { U+ F^{4c,2e}[D^{4c}_(+)atoms] U }^{++}

  !           !> read the U matrix
              ioff_U = 1
              do i = 1, nfsym
                write(flabel,'(a11,i1)') 'pctmtAO   1',i
                if(naosh_all(i) * naosh_L(i) > 0)then
                  call x2c_read(flabel,scr2_mat(ioff_U),naosh_all(i) * naosh_L(i)*nz,x2c_file_glb)
                end if
                ioff_U = ioff_U + naosh_all(i) * naosh_L(i) * nz
              end do

              !> read F^{4c,2e}[D^{4c}_(+)atoms]
              open(103,file='eamfX2C_F4c',status='old',form='unformatted',  &
                   access='sequential',action="readwrite",position='rewind')
              read(103) scr3_mat(1:naosh_ls*naosh_ls*nz)
              close(103,status='keep')

              !> picture-change transformation { U+ F^{4c,2e}[D^{4c}_(+)atoms] U }^{++}
              op_fer_rep = 0

              call perform_pct_saao_bas(scr3_mat,                     & ! on output: { U+ F^{4c,2e}[D^{4c}_(+)atoms] U }^{++}
                                        scr2_mat,                     & ! U
                                        scr1_mat,                     & ! scratch
                                        naosh_ls,                     &
                                        nr_ao_L,                      &
                                        naosh_all,                    &
                                        naosh_L,                      &
                                        nfsym,                        &
                                        nz,                           &
                                        ipq_off_in,                   &
                                        ioff_aomat_x,                 &
                                        op_fer_rep,                   &
                                        op_bs_to_fs,                  &
                                        .true.,                       &
                                        print_lvl)

!             debug print
              if(print_lvl > 2)then
                call print_x2cmat(scr3_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - eamf term 1 contributions',6)
              end if

              !> - second term + first term => complete eamf PCE correction
              call daxpy(nr_ao_L*nr_ao_L*nz,-1.0d0,scr4_mat,1,scr3_mat,1);
              deallocate(scr4_mat);
!             debug print
              if(print_lvl > 2)then
                call print_x2cmat(scr3_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - eamf contributions',6)
              end if

       end if

       !> collect correction term for XC energy
       if(isDFT)then
         delta_xc_e   = 0.0d0
         center_total = 1
         do i = 1, nuclei_type
           delta_xc_e  = delta_xc_e + collect_delta_xc_energy(i,              &
                                                              center_total,   &
                                                              print_lvl       &
                                                             )
         end do
         call set_xc_energy_delta(delta_xc_e)
       end if

#define DEBUG_SOC
#ifdef DEBUG_SOC
       open(99,file='soc-contributions',status='replace',form='formatted',  &
       access='sequential',action="readwrite",position='rewind')

       call print_x2cmat(scr3_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - X-AMFI contribution',99)
       close(99,status='keep')
#undef DEBUG_SOC
#endif

       open(99,file='XAMFI-ss-soc-contributions',status='replace',form='unformatted',  &
       access='sequential',action="readwrite",position='rewind')

       write(99) scr3_mat(1:nr_ao_L*nr_ao_L*nz)
       close(99,status='keep') 

!      -------------------------------------------------------------------------
!      read h1_2c_saao_cartesian_gaussian from file and add aoosoc contributions 
!      -------------------------------------------------------------------------

       scr2_mat = 0; scr1_mat = 0;

       i = 0; j = 0
       write(flabel,'(a7,i4,i1)') 'h12cAOa',1,i
       call x2c_read(flabel,scr1_mat,naosh_ls*naosh_ls*nz,x2c_file_amf)

!      -----------------------------------------------------------
!      select ++ block of h1_++_saao_cartesian_gaussian -> h1_{2c}
!      -----------------------------------------------------------
!          in: scr1_mat
!         out: scr2_mat
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
!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr2_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - h1 w/o (e)amf contribution',6)
       end if

!      -----------------------------------
!      add (e)amf contributions to h1_{2c}
!      -----------------------------------
       call daxpy(nr_ao_L*nr_ao_L*nz,1.0d0,scr3_mat,1,scr2_mat,1)

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr2_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - h1 + (e)amf contributions',6)
       end if

!      ----------------------------
!      write h1_{2c} to file X2CMAT 
!      ----------------------------
       i = 0
       j = 0
       write(flabel,'(a7,i4,i1)') 'h12cAOA',1,i

       if(spherical_on)then

         allocate(wrk(lwrk))
         call sph_iotc(scr2_mat,scr1_mat,nr_sph,print_lvl,wrk,lwrk)
         deallocate(wrk)

         call x2c_write(flabel,scr1_mat,nr_sph**2  * nz,x2c_file_glb)
       else

!        debug print
         if(print_lvl > 2)then
           call print_x2cmat(scr2_mat,nr_ao_L,nr_ao_L,nz,ipq_off_in,'x2c - h12c+amf',6)
         end if

         call x2c_write(flabel,scr2_mat,nr_ao_L**2 * nz,x2c_file_glb)
       end if

!      release scratch space
       deallocate(scr1_mat,scr2_mat,scr3_mat)

  end subroutine x2c_interface2aoosoc

!**********************************************************************

  subroutine x2c_add_2e_diagonal_atomic(&
       h1_4c_defining,        &
       naosh_ls,              &
       nfsym,                 &
       nz,                    &
       bs_irrep_mat,          &
       ipq_off_in,            &
       iqp_off_in,            &
       iqmult_trip_q,         &
       op_bs_to_fs,           &
       quat_pointer,          &
       nuclei_totnum,         &
       nuclei_type,           &
       print_lvl)   
!**********************************************************************
!
!    purpose: interface routine to add 2e-mean field atomic blocks to h1
!
!             the outline of this routine may be sketched as follows:
!
!             1. retrieve 2e-mean field atomic block from file
!   "saw"   /   a. transform from simple AO basis to Hermit-sorted SO-AO basis
!  -------  \   b. transform from Hermit-sorted SO-AO basis to Dirac-sorted AO (SA-AO) basis
!            \  c. insert quaternion phases
!             2. add atomic block to h1_4c
! 
!
!----------------------------------------------------------------------
real*8,  INTENT(INOUT)          :: h1_4c_defining(*)
integer, intent(in)             :: naosh_ls
integer, intent(in)             :: nfsym
integer, intent(in)             :: nz
integer, intent(in)             :: ipq_off_in(4,0:7)
integer, intent(in)             :: iqp_off_in(4,0:7)
integer, intent(in)             :: op_bs_to_fs(0:7,1:2)
integer, intent(in)             :: bs_irrep_mat(4,0:7)
integer, intent(in)             :: iqmult_trip_q(4,4,4)
integer, intent(in)             :: quat_pointer(0:7,2)
integer, intent(in)             :: nuclei_totnum
integer, intent(in)             :: nuclei_type
integer, intent(in)             :: print_lvl
!--------------------------------------------------------------
real(8), allocatable            :: scr1_mat(:), scr2_mat(:)
integer, allocatable            :: aooso_info(:)
integer                         :: i
integer                         :: j
integer                         :: iq
integer                         :: ipq
integer                         :: iz
integer                         :: irepd
integer                         :: nsphcm
integer                         :: nr_sph
integer                         :: nrows_frag
integer                         :: ncols_frag
integer                         :: center_total
integer                         :: wf_components
INTEGER, PARAMETER              :: nzc1 = 4
character (len=12)              :: flabel
!**********************************************************************

!      allocate scratch space
allocate(scr1_mat(naosh_ls*naosh_ls*nzc1),scr2_mat(naosh_ls*naosh_ls*nzc1)); scr1_mat = 0; scr2_mat = 0

!      debug print
if(print_lvl > 2)then
call print_x2cmat(h1_4c_defining,naosh_ls,naosh_ls,nz,ipq_off_in,'x2c - defining h1_4c wo aoosoc contribution',6)
end if

!      ---------------------------------------------------------------------------------
!      step 1. read and distribute aoosoc contributions for symmetry independent centers
!      ---------------------------------------------------------------------------------
!        out: atomic 2e-mean field blocks  -> scr1_mat

allocate(aooso_info(naosh_ls))

center_total = 1
scr1_mat     = 0

!> a. read contributions for symmetry independent centers
do i = 1, nuclei_type

!        count the number of basis functions for center_total
call labcount(ncols_frag,aooso_info,naosh_ls,1,-1,center_total,1,-1)
nrows_frag = ncols_frag

#ifdef MOD_XAMFI
!> read 2e-mean-field contributions (assumed to be exported in C1 symmetry)
call put_atomic_2e_mean_field_correction(                &
        scr1_mat,       &
        naosh_ls,       &
        nrows_frag,     &
        ncols_frag,     &
        nzc1,           &
        i,              &
        center_total,   &
        print_lvl       &
       )
#else
call quit('X-AMFI module not avaliable in this version.')
#endif
end do

!      b. transform from AO-basis (C1* symm sorted) to SO-AO basis (symmetry-ordered)
!      ------------------------------------------------------------------------------
!        in : 2e_atomic_AO    -> scr1_mat
!        out: 2e_atomic-SO-AO -> scr2_mat
i = op_bs_to_fs(0,1)
do iz = 1, nzc1
irepd = bs_irrep_mat(iz,0)
iq    = iqmult_trip_q(1,quat_pointer(irepd,i),iz)
ipq   = iqp_off_in(iq,0)
call mtaoso(scr1_mat(1+(naosh_ls*naosh_ls*(iz-1))),     &
scr2_mat(1+(naosh_ls*naosh_ls*(ipq-1))),                &
naosh_ls,                                               &
irepd,                                                  &
print_lvl)
end do

!      c. transform from Hermit-sorted AO to DIRAC-sorted AO basis
!      -----------------------------------------------------------
call butobs_no_work(scr2_mat,nz)

!      d. insert quaternion phase factors
!      ----------------------------------
if(nz < 4)then
do iz = 1, nz
iq = ipq_off_in(iz,0)
call q2bphase('F',iq,1,scr2_mat(1+(naosh_ls*naosh_ls*(iz-1))))
end do
end if

!      debug print
if(print_lvl > 2)then
call print_x2cmat(scr2_mat,naosh_ls,naosh_ls,nz,ipq_off_in,'x2c - 2e-diagonal contributions',6)
end if

deallocate(aooso_info)

!      -----------------------------------------
!      step 2. add XAMFI-SO contributions to h1_2c
!      -----------------------------------------

call daxpy(naosh_ls*naosh_ls*nz,1.0d0,scr2_mat,1,h1_4c_defining,1)

!      debug print
if(print_lvl > 2)then
call print_x2cmat(h1_4c_defining,naosh_ls,naosh_ls,nz,ipq_off_in,'x2c - h1 + diagonal contribution',6)
end if

! release scratch space
deallocate(scr1_mat,scr2_mat)

end subroutine x2c_add_2e_diagonal_atomic
!**********************************************************************

subroutine x2c_add_2e_from_atomic_dens(&
       h1_Xc,                 &
       naosh_ls,              &
       nz,                    &
       bs_irrep_mat,          &
       ipq_off_in,            &
       iqp_off_in,            &
       iqmult_trip_q,         &
       op_bs_to_fs,           &
       quat_pointer,          &
       nuclei_totnum,         &
       nuclei_type,           &
       xcmode,                &
       integral_types,        &
       lwrk,                  &
       isDFT,                 &
       print_lvl)
!**********************************************************************
!
!    purpose: interface routine to add 2e-molecular mean field contributions to an arbitrary h1 (could be 0)
!             calculated from a superposition of Xc [X=2,4] atomic densities
!
!             the outline of this routine may be sketched as follows:
!
!             1. retrieve atomic densities D^Xc_atom from file
!             2. form molecular density D^Xc_(+)atoms as superposition of atomic densities: D^Xc_(+)atoms = (+) D^Xc_atom
!   "saw"   /   a. transform from simple AO basis to Hermit-sorted SO-AO basis
!  -------  \   b. transform from Hermit-sorted SO-AO basis to Dirac-sorted AO (SA-AO) basis
!            \  c. insert quaternion phases
!             4. compute Fock matrix F^{Xc,2e}[D^Xc_(+)atoms]
!             5. add F^{Xc,2e}D^Xc_(+)atoms] to h1_Xc
!
!
!----------------------------------------------------------------------
real*8,  INTENT(INOUT)          :: h1_Xc(*)
integer, intent(in)             :: naosh_ls
integer, intent(in)             :: nz
integer, intent(in)             :: ipq_off_in(4,0:7)
integer, intent(in)             :: iqp_off_in(4,0:7)
integer, intent(in)             :: op_bs_to_fs(0:7,1:2)
integer, intent(in)             :: bs_irrep_mat(4,0:7)
integer, intent(in)             :: iqmult_trip_q(4,4,4)
integer, intent(in)             :: quat_pointer(0:7,2)
integer, intent(in)             :: nuclei_totnum
integer, intent(in)             :: nuclei_type
integer, intent(in)             :: xcmode
integer, intent(in)             :: integral_types
integer, intent(in)             :: lwrk
integer, intent(in)             :: print_lvl
logical, intent(in)             :: isDFT
!--------------------------------------------------------------
real*8 , allocatable            :: scr1_mat(:), scr2_mat(:)
integer, allocatable            :: aooso_info(:)
integer                         :: i, lwrk_
integer                         :: j
integer                         :: iq
integer                         :: ipq
integer                         :: iz
integer                         :: irepd
integer                         :: nsphcm
integer                         :: nr_sph
integer                         :: nrows_frag
integer                         :: ncols_frag
integer                         :: center_total
integer                         :: wf_components
INTEGER, PARAMETER              :: nzc1 = 4
character (len=12)              :: flabel
logical                         :: TWOCOMPBSS_
!**********************************************************************

!      allocate scratch space
allocate(scr1_mat(naosh_ls*naosh_ls*nzc1),scr2_mat(naosh_ls*naosh_ls*nzc1)); scr1_mat = 0; scr2_mat = 0

!      debug print
if(print_lvl > 2)then
call print_x2cmat(h1_Xc,naosh_ls,naosh_ls,nz,ipq_off_in,'x2c - incoming h1_Xc w/o eamf contributions',6)
end if

!      ---------------------------------------------------------------------------------
!      step 1. read and distribute atomic density contributions for symmetry independent centers
!      ---------------------------------------------------------------------------------
!        out: molecular density from atomic blocks  -> scr1_mat

allocate(aooso_info(naosh_ls))

center_total = 1
scr1_mat     = 0

!> a. read contributions for symmetry independent centers
do i = 1, nuclei_type

!        count the number of basis functions for center_total
call labcount(ncols_frag,aooso_info,naosh_ls,1,-1,center_total,1,-1)

if(ncols_frag == 0) cycle

nrows_frag = ncols_frag

#ifdef MOD_XAMFI
!> read atomic density (assumed to be exported in C1 symmetry)
call put_atomic_density_correction(                &
        scr1_mat,       &
        naosh_ls,       &
        nrows_frag,     &
        ncols_frag,     &
        nzc1,           &
        i,              &
        center_total,   &
        xcmode,         &
        print_lvl       &
       )
#else
call quit('X-AMFI module not avaliable in this version.')
#endif
end do

!      b. transform from AO-basis (C1* symm sorted) to SO-AO basis (symmetry-ordered)
!      ------------------------------------------------------------------------------
!        in : DM_AO    -> scr1_mat
!        out: DM-SO-AO -> scr2_mat
i = op_bs_to_fs(0,1)
do iz = 1, nzc1
irepd = bs_irrep_mat(iz,0)
iq    = iqmult_trip_q(1,quat_pointer(irepd,i),iz)
ipq   = iqp_off_in(iq,0)
call dtaoso(scr1_mat(1+(naosh_ls*naosh_ls*(iz-1))),     &
scr2_mat(1+(naosh_ls*naosh_ls*(ipq-1))),                &
naosh_ls,                                               &
irepd,                                                  &
print_lvl)
end do

!      c. transform from Hermit-sorted AO to DIRAC-sorted AO basis
!      -----------------------------------------------------------
call butobs_no_work(scr2_mat,nz)

!      d. insert quaternion phase factors
!      ----------------------------------
if(nz < 4)then
do iz = 1, nz
iq = ipq_off_in(iz,0)
call q2bphase('D',iq,1,scr2_mat(1+(naosh_ls*naosh_ls*(iz-1))))
end do
end if

!      debug print
if(print_lvl > 2)then
call print_x2cmat(scr2_mat,naosh_ls,naosh_ls,nz,ipq_off_in,'x2c - mol DM assembled from atomic DMs',6)
end if

deallocate(aooso_info)

!      ----------------------------------------------------
!      step 3. compute Fock matrix F^{Xc,2e}[D^Xc_(+)atoms]
!      ----------------------------------------------------
lwrk_ = lwrk;
TWOCOMPBSS_ = TWOCOMPBSS
if(xcmode ==4)then
       TWOCOMPBSS = .false.
else
       TWOCOMPBSS = .true.
end if

call get_x_fock(                        &
                scr2_mat,               & ! DMAT
                scr1_mat,               & ! FOCK
                naosh_ls,               &
                naosh_ls,               &
                1,                      & ! one Fock matrix (remember that we use an [approximate] FOCC density matrix)
                integral_types,         & ! integral flag: which type of integrals to include
                nz,                     &
                isDFT,                  & !> if true: add XC potential
                lwrk_,                  &
                print_lvl               &
      )
TWOCOMPBSS_ = TWOCOMPBSS
!      ------------------------------------------------------
!      step 4. add Fock matrix F^{Xc,2e}[D^Xc_(+)atoms] to h1
!      ------------------------------------------------------
call daxpy(naosh_ls*naosh_ls*nz,1.0d0,scr1_mat,1,h1_Xc,1)

!      debug print
if(print_lvl > 2)then
call print_x2cmat(h1_Xc,naosh_ls,naosh_ls,nz,ipq_off_in,&
                  'x2c - h1 + F^{Xc,2e}[D^Xc_(+)atoms] contribution',6)
end if

! save F^{Xc,2e}[D^Xc_(+)atoms] on file for possible reuse
if(xcmode == 4)then
open(103,file='eamfX2C_F4c',status='replace',form='unformatted',  &
     access='sequential',action="readwrite",position='rewind')
write(103) scr1_mat(1:naosh_ls*naosh_ls*nz)
close(103,status='keep')
end if
! release scratch space
deallocate(scr1_mat,scr2_mat)

end subroutine x2c_add_2e_from_atomic_dens

!**********************************************************************


end module x2c_2e_soc_interface
