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
!
module module_aoosoc_init

! stefan:
!          this module is the interface between the aoosoc code and common blocks
!          ideally no common blocks should be included in the aoosoc code
!          outside this module to keep the interface clean
!          hide as much as you can for public - it will save us lots of trouble in future

  use atomic_oo_order_so_correction_cfg
  use module_aoosoc_pfg

  implicit none

  public aoosoc_init

  private

contains

  subroutine aoosoc_init()

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
!   local variables
!   ----------------------------------------------------------------------------
    integer             :: i, j, k, l, m
    integer             :: ier
!   ----------------------------------------------------------------------------

      aoo_prt               = iprham
      if(aoomod_debug)      &
      aoo_prt               = 3
      nr_quat               = nz
      nzt_aoo               = nzt
      if(mdirac) nzt_aoo    = nz ! default: 1
      nr_fsym               = nfsym
      mfsym                 = min(nfsym,2)

!     initialize integer variables
      aoo_num_nuclei_i = nucind
      aoo_num_nuclei_d = nucdep
!                       [  L+S  ] (ifsym==1) + [  L+S  ] (ifsym==2)
      nr_ao_total_aoo  = ntbas(0)
      nrows            = ntbas(0)
      ncols            = ntbas(0)
!                       [  L    ] (ifsym==1) + [  L    ] (ifsym==2)
      nr_ao_large_aoo  = ntbas(1)
!                       [    S  ] (ifsym==1) + [    S  ] (ifsym==2)
      nr_ao_small_aoo  = ntbas(2)
!
      n2bastq_dim_aoo  = n2bastq
!
      nr_cmo_q         = ncmotq
!
      do i = 1, nr_fsym
        j               =  nesh(i)
        k               =  norb(i)
        l               =  npsh(i)
        dim_pshell(i)   =  l
        dim_eshell(i)   =  j
        dim_oshell(i)   =  k
!                         [  L+S  ]
        nr_ao_all(i)    = nfbas(i,0)
!                         [  L    ]
        nr_ao_l(i)      = nfbas(i,1)
!                         [    S  ]
        nr_ao_s(i)      = nfbas(i,2)

        do m = 1, nr_fsym
          ioff_aomat_x(m,i) = i2basx(m,i)
        end do
      end do

      aoo_cb_pq_to_uq(1:4, 0:7)        = ipqtoq(1:4, 0:7)
      aoo_cb_uq_to_pq(1:4, 0:7)        = iqtopq(1:4, 0:7)
      aoo_bs_to_fs(0:7, 1:2)           = jbtof(0:7, 1:2)
      aoo_pointer_quat(0:7, 1:2)       = jqbas(0:7, 1:2)
      aoo_bs_irrep_mat(1:4, 0:7)       = irqmat(1:4, 0:7)
      aoo_iqmult_trip_q(1:4, 1:4, 1:4) = iqmult(1:4, 1:4, 1:4)
      aoo_pointer_quat_op(0:7)         = jm4pos(0:7)
      aoo_ihqmat(1:4,-1:1)             = ihqmat(1:4,-1:1)

!     initialize logicals
      aoo_mdirac       = mdirac
  
!     initialize real*8
      aoo_cspeed = cval

      if(mc > 1)then
!       4c-component picture
        nr_2e_fock_matrices        = nfmat
        aoo_intflg                 = intflg
      end if

      j = nucind
!     j = 0
!     do i = 1, nucind
!       do k = 1, nucdeg(i)
!         j = j + 1
!       end do
!     end do

!    active fock matrix factors
     if(nr_2e_fock_matrices > 1)then 
       allocate( aoo_dfopen(nr_2e_fock_matrices), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: aoo_dfopen(:)'
       aoo_dfopen(1:(nr_2e_fock_matrices-1)) = df(1:(nr_2e_fock_matrices-1))
     end if

    !> interface atomic/molecular data
       allocate( atom(j), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: atom(:)'
       do i=1,size(atom)
         atom(i) = type_atom( charge(i), &
                              cord(1,i), &
                              cord(2,i), &
                              cord(3,i)  &   
                            )
       end do
       allocate( aoo_nont(mxatom), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: aoo_nont(:)'
       allocate( aoo_nucdeg(mxcent), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: aoo_nucdeg(:)'
       aoo_mxatom           = mxatom
       aoo_mxcent           = mxcent
       aoo_nontyp           = nontyp
       aoo_nont(1:mxatom)   = nont(1:mxatom)
       aoo_nucdeg(1:mxcent) = nucdeg(1:mxcent)
    !> end

! 
  end subroutine
!**********************************************************************

end module
