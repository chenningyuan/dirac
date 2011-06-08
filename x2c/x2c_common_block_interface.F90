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
!
module x2c_cb_interface

! stefan (following the advice of Radovan):  
!          this module is the interface between the x2c code and common blocks
!          ideally no common blocks should be included in the x2c code
!          outside this module to keep the interface clean
!          hide as much as you can for public - it will save us lots of trouble in future

  use x2cmod_cfg

  implicit none

  public init_x2c_cb
  public set_wrk_x2c
  public sel_wrk_x2c
  public desel_wrk_x2c
  public reset_x2c_cb_lin_lwdmat
  public reset_x2c_cb_onefck
  public renew_x2c_cb_orb_shell_dim

  private

  real(8), pointer :: f77_wrk_x2c(:)      ! f77 work array pointer
  integer, public  :: len_wrk_f77_x2c = 0 ! length of the work array
  integer          :: mem_wrk_use_x2c = 0 ! amount of work array used
  integer          :: mem_wrk_cnt_x2c = 0 ! counts

contains

  subroutine init_x2c_cb(init_x2c_cb_blocks,   &
                         prop_op_dgsym,        &
                         prop_op_trsym,        &
                         prop_op_name)

      use quaternion_algebra
#include "dcbbas.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "maxaqn.h"
#include "ccom.h"
#include "dcbprp.h"
#include "dcbxpr.h"
#include "dcbdhf.h"
#include "priunit.h"

!   ----------------------------------------------------------------------------
    logical, intent(in)                         :: init_x2c_cb_blocks
    integer, optional, intent(out)              :: prop_op_dgsym(*)
    integer, optional, intent(out)              :: prop_op_trsym(*)
    character (len=16), optional, intent(out)   :: prop_op_name(*)
!   ----------------------------------------------------------------------------
!   local variables
!   ----------------------------------------------------------------------------
    integer             :: i, j, k, l, m
!   ----------------------------------------------------------------------------

       if(init_x2c_cb_blocks)then
         x2c_prt               = iprham
         if(x2cmod_debug)      &
         x2c_prt               = 3
         nr_quat               = nz
         nr_fsym               = nfsym
         mfsym                 = min(nfsym,2)

!        set value for amfi module (if at all)
         if(noamfi)then
           x2c_add_amfi = -1
         end if

!        common block variables for AMFI(old +new(?))
         if(x2c_add_amfi > 0)then
           x2c_tot_charge_for_mfsum = imfch
           x2c_amfi_order           = isorder_amfi_x2c
           x2c_max_quant_num        = mxqn
         end if

!        set file unit on common block in dcbgen.h
         lux2c = x2c_funit

!
!        initialize integer variables
         num_nuclei      = nucdep
!                         [  L+S  ] (ifsym==1) + [  L+S  ] (ifsym==2)
         nr_ao_total_x2c = ntbas(0)
!                         [  L    ] (ifsym==1) + [  L    ] (ifsym==2)
         nr_ao_large_x2c = ntbas(1)
!                         [    S  ] (ifsym==1) + [    S  ] (ifsym==2)
         nr_ao_small_x2c = ntbas(2)
!
         n2bastq_dim_x2c = n2bastq
!
         nr_cmo_q        = ncmotq
!
         fullemat_dim    = 0
         fullomat_dim    = 0
         fulleomat_dim   = 0
         do i = 1, nr_fsym
           j               =  nesh(i)
           k               =  norb(i)
           l               =  npsh(i)
           dim_pshell(i)   =  l
           dim_eshell(i)   =  j
           dim_eshell2(i)  =  j*2
           dim_e2shell(i)  =  j**2
           dim_e2shellq(i) = (j**2)*nr_quat
           dim_oshell(i)   =  k
           dim_o2shellq(i) = (k**2)*nr_quat
           fullemat_dim    = (j**2)*nr_quat + fullemat_dim 
           fullomat_dim    = (k**2)*nr_quat + fullomat_dim 
           fulleomat_dim   =  k* j *nr_quat + fulleomat_dim 
!                            [  L+S  ]
           nr_ao_all(i)    = nfbas(i,0)
!                            [  L    ]
           nr_ao_l(i)      = nfbas(i,1)
!                            [    S  ]
           nr_ao_s(i)      = nfbas(i,2)
           do m = 1, nr_fsym
             ioff_aomat_x(m,i) = i2basx(m,i)
           end do
         end do

         x2c_cb_pq_to_uq(1:4, 0:7)        = ipqtoq(1:4, 0:7)
         x2c_cb_uq_to_pq(1:4, 0:7)        = iqtopq(1:4, 0:7)
         x2c_bs_to_fs(0:7, 1:2)           = jbtof(0:7, 1:2)
         x2c_pointer_quat(0:7, 1:2)       = jqbas(0:7, 1:2)
         x2c_bs_irrep_mat(1:4, 0:7)       = irqmat(1:4, 0:7)
         x2c_iqmult_trip_q(1:4, 1:4, 1:4) = iqmult(1:4, 1:4, 1:4)
         x2c_pointer_quat_op(0:7)         = jm4pos(0:7)

         fh_1int_4c = lu1int
         all_prp_op = nprps

!        initialize logicals
         x2c_linsym       = linear
         x2c_do_spherical = dosphe
         x2c_prep_bnccorr = bncron.and..not.onesys
         x2c_1e_system    = onesys
  
!        initialize real*8
         x2c_cspeed = cval

!        set defining matrix wrt which the pct-mtx will be derived
         if(x2c_4c_fock_mtx_defh1)then
!          converged 4c-Fock-Dirac operator
           x2c_is_defining_h1mat      =  2
           scf_iter_counter           = -1
           diis_counter               = -1
           nr_2e_fock_matrices        = nfmat
           file_name_1e_fock_matrix   = 'DFFCK1'
           file_name_2e_fock_matrices = 'DFFCK2'
           x2c_prep_bnccorr           = .false.
         else if(x2c_free_part_mtx_defh1)then
!          free-particle matrix
           x2c_is_defining_h1mat  = 3
!        else if(dohuckel)then
!          (extended) Huckel start guess
!          x2c_is_defining_h1mat  = 1
         else
!          default: bare-nucleus one-electron Dirac-Hamiltonian
           x2c_is_defining_h1mat  = 0
         end if
! 
       else 

         if(present(prop_op_dgsym) .and. present(prop_op_trsym) .and.  present(prop_op_name))then
           do i = 1, nprps
             prop_op_dgsym(i) = iprpsym(i)
             prop_op_trsym(i) = iprptim(i)
             prop_op_name(i)  = prpnam(i)
           end do
         else 
           do i = 1, nr_fsym
!            orbital dimensions after LOWDIN (needed for a proper handling of 2c-LOWDIN matrix)
             nr_mo_lw_all(i) = nforb(i,0)
             nr_mo_lw_l(i)   = nforb(i,1)
             nr_mo_lw_s(i)   = nforb(i,2)
           end do
         end if

       end if ! initialize x2c information from Dirac common blocks

  end subroutine

!**********************************************************************
  subroutine reset_x2c_cb_lin_lwdmat(set_cb_blocks,        &
                                     npsh_dim,             &
                                     nesh_dim,             &
                                     nfsym,                &
                                     nz,                   &
                                     nborb_dim)
#include "dcborb.h"
!**********************************************************************
!   
!    purpose: reset orbital dimensions on common block in dcborb.h 
!             for "call linsym" inside the X2C module
!
!**********************************************************************
     integer, intent(in)    :: set_cb_blocks
     integer, intent(in)    :: npsh_dim(nfsym)
     integer, intent(in)    :: nesh_dim(nfsym)
     integer, intent(in)    :: nfsym
     integer, intent(in)    :: nz
     integer, intent(inout) :: nborb_dim(4,2,0:2)
!-----------------------------------------------------------------------
     integer                :: i
     integer                :: j
     integer                :: k
!**********************************************************************

!      initialize
       j = 4/nz

!      modify common block
       if(set_cb_blocks == 1)then

         do i = 1, nfsym
           npsh(i) = 0
           norb(i) = nesh_dim(i) + npsh(i)
           do k = 1, j
             nborb_dim(k,i,2) = nborb(k,i,2)
             nborb(k,i,2)     = 0
             nborb(k,i,0)     = nborb(k,i,1) + nborb(k,i,2)
           end do
         end do

!      restore common block
       else

         do i = 1, nfsym
           npsh(i) = npsh_dim(i)
           norb(i) = nesh_dim(i) + npsh(i)
           do k = 1, j
             nborb(k,i,2) = nborb_dim(k,i,2)
             nborb(k,i,0) = nborb(k,i,1) + nborb(k,i,2)
           end do
         end do
       end if

  end subroutine reset_x2c_cb_lin_lwdmat

!**********************************************************************
  subroutine reset_x2c_cb_onefck(x2c_run_vl_tmp)

  use x2cmod_cfg

#include "dcbham.h"
!**********************************************************************
!   
!    purpose: reset X2C flags on modules/common block (dcbham.h) 
!             for the "call onefck" inside the X2C module
!
!**********************************************************************
     logical, intent(in)    :: x2c_run_vl_tmp
!-----------------------------------------------------------------------
!**********************************************************************

       x2cmod_x2c = x2c_run_vl_tmp
!      modify common block
       X2C        = x2c_run_vl_tmp
!      BSS        = x2c_run_vl_tmp
!      TWOCOMPBSS = x2c_run_vl_tmp

  end subroutine reset_x2c_cb_onefck

!**********************************************************************
  subroutine renew_x2c_cb_orb_shell_dim(ham_lvl_x2c)

  use x2cmod_cfg
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dgroup.h"
!**********************************************************************
!   
!    purpose: 1. after the X2C procedure:
!             - renew orbital dimensions on common block(s) in 
!               dcborb.h and dcbbas.h. 
!             - renew logicals, parameters and scaling-factors related to 
!               small-component handlings on common block(s) in 
!               dcbgen.h, dcbham.h, dcbdhf.h and dgroup.h.

!               This is equivalent to a "removal" of the positronic shells.
!
!             2. set up the Hamiltonian flag "i2cofk" which points to the 
!                final Hamiltonian integrals present on the file X2CMAT.
!             
!
!**********************************************************************
     integer, intent(in)    :: ham_lvl_x2c
!-----------------------------------------------------------------------
     integer                :: i
     integer                :: isym
     integer                :: nsym
!**********************************************************************

!      step 1: modify shell dimension on common block
!      ----------------------------------------------
       nsym = 4/nz
       do i=1, nfsym
         npsh(i) = 0
         norb(i) = nesh(i) + npsh(i)
         ntmo(i) = norb(i)
         do isym = 1, nsym
           nborb(isym,i,2) = 0
           nborb(isym,i,0) = nborb(isym,i,1) + nborb(isym,i,2)
         end do
       end do

!      step 2: recalculate offsets and globally used pointers
!      ----------------------------------------------
       call setdc2(0)
!      if(x2c_4c_fock_mtx_defh1) call setdc1(0)

!      step 3: reset integral default, overlap metric and # of transformation matrices
!      -------------------------------------------------------------------------------
       intdef_save = intdef
       ssmtrc_save = ssmtrc
       intgen = 1
       intdef = 1
       ssmtrc = 0.0d0

!      # number of wave function components (1: large; 2: large + small) 
       mc = 1

!      # of transformation matrices == 1 --> Lowdin ao2mo
       nzt = 1

!      final hamiltonian level for matrix stored on file X2CMAT
       i2cofk = ham_lvl_x2c

  end subroutine renew_x2c_cb_orb_shell_dim
!**********************************************************************
!-----------------------------------------------------------------------
! rip from openrsp/dirac_interface.F90 / prp/london_reorth.F90

!**********************************************************************
  subroutine set_wrk_x2c(wrk,lwork)
!**********************************************************************
     real(8),target     :: wrk(:)
     integer            :: lwork
!**********************************************************************

     len_wrk_f77_x2c =  lwork
     f77_wrk_x2c     => wrk

  end subroutine set_wrk_x2c

!**********************************************************************
  subroutine sel_wrk_x2c(wrk, mem_req)
!**********************************************************************
     real(8), pointer    :: wrk(:)
     integer, intent(in) :: mem_req
!**********************************************************************

     wrk => f77_wrk_x2c
     mem_wrk_use_x2c = mem_req
     mem_wrk_cnt_x2c = mem_wrk_cnt_x2c + 1
     if(mem_wrk_cnt_x2c > 1)                               &
        call quit(' *** dirac/dirx2c.F90 - sel_wrk_x2c: counter > 1 not implemented ***')

  end subroutine sel_wrk_x2c

!**********************************************************************
  subroutine desel_wrk_x2c(wrk, mem_req)
!**********************************************************************
     real(8), pointer    :: wrk(:)
     integer, intent(in) :: mem_req
!**********************************************************************

     nullify(wrk)
     mem_wrk_use_x2c = 0
     mem_wrk_cnt_x2c = mem_wrk_cnt_x2c - 1

  end subroutine desel_wrk_x2c
!-----------------------------------------------------------------------
! end of rip
end module
