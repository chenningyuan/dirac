!dirac_copyright_start
!      Copyright (c) 2021 by the authors of DIRAC.
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
!
module xamfi_environment

! stefan:
!          this module is the interface between the aoosoc code and common blocks
!          ideally no common blocks should be included in the aoosoc code
!          outside this module to keep the interface clean
!          hide as much as you can for public - it will save us lots of trouble in future

  use xamfi_global_parameters
  use xamfi_internal_parameters

  implicit none

  public xamfi_init
  public xamfi_finalize
  public xamfi_set_cb

  private

contains

  subroutine xamfi_init()

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
      n2bbasxq_dim_aoo = n2bbasxq
!
      nr_cmo_q         = ncmotq
!
      aoo_n2tmt        = n2tmt
      lutmat_aoo       = lutmat

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

        nr_tmo(i)       = ntmo(i)
        ioff_tmt(i)     = i2tmt(i)
        ioff_tmot(i)    = i2tmot(i)

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

!     initialize data if we are still in the 4c-component picture
      if(mc > 1)then
                    nr_2e_fock_matrices        = nfmat
        if(aoofocc) nr_2e_fock_matrices        = 1
        aoo_intflg                             = intflg
        aoo_nopen                              = nopen
      end if

      aoo_dfopen(0) = df(0)
!     active density matrix factors
      if(aoo_nopen >=1) aoo_dfopen(1:aoo_nopen) = df(1:aoo_nopen)

      j = nucind
!     j = 0
!     do i = 1, nucind
!       do k = 1, nucdeg(i)
!         j = j + 1
!       end do
!     end do

    !> interface atomic/molecular data
       allocate( atom(j), stat=ier ); if( ier.ne.0 )stop ' aoosoc init: Error in allocation: atom(:)'
       do i=1,size(atom)
         atom(i) = type_atom( charge(i), &
                              cord(1,i), &
                              cord(2,i), &
                              cord(3,i)  &   
                            )
       end do
       allocate( aoo_nont(mxatom), stat=ier ); if( ier.ne.0 )stop ' aoosoc init: Error in allocation: aoo_nont(:)'
       allocate( aoo_nucdeg(mxcent), stat=ier ); if( ier.ne.0 )stop ' aoosoc init: Error in allocation: aoo_nucdeg(:)'
       aoo_mxatom           = mxatom
       aoo_mxcent           = mxcent
       aoo_nontyp           = nontyp
       aoo_nont(1:mxatom)   = nont(1:mxatom)
       aoo_nucdeg(1:mxcent) = nucdeg(1:mxcent)
    !> end

! 
  end subroutine xamfi_init
!**********************************************************************

  subroutine xamfi_finalize()

     deallocate(atom)
     deallocate(aoo_nont)
     deallocate(aoo_nucdeg)

  end subroutine xamfi_finalize
!**********************************************************************

  subroutine xamfi_set_cb(myway)

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

    character(len=6), intent(in) :: myway
!   ----------------------------------------------------------------------------
!   local variables
!   ----------------------------------------------------------------------------
    integer             :: i, j
    character(len=1)    :: c
!   ----------------------------------------------------------------------------

    c = myway(1:1)
    select case(c)
      case ('s')
            
            !print *, 'saving cb...',myway(5:6)
            open(99,file='xamfi-'//myway(5:6),status='unknown',form='unformatted',   &
            access='sequential',action='readwrite',position='rewind')

            nzt_aoo_save            = nzt
            if(mdirac) nzt_aoo_save = nz ! default: 1
      
            write(99)   nz, nzt_aoo_save,ntbas(0),ntbas(0),ntbas(0),ntbas(1),ntbas(2),&
                        n2bastq, n2bbasxq, ncmot, n2tmt
            !print *, '   stored values are ',nz, nzt,ntbas(0),ntbas(0),ntbas(0),ntbas(1),ntbas(2),&
            !n2bastq, n2bbasxq, ncmot, n2tmt
            do i= 1, nr_fsym
                  write(99) &
                        npsh(i),&
                        nesh(i),&
                        norb(i),&
                        nfbas(i,0),&
                        nfbas(i,1),&
                        nfbas(i,2),&
                        ntmo(i),&
                        i2tmt(i),&
                        i2tmot(i)
                  do j = 1, nr_fsym
                        write(99)&
                        i2basx(i,j)
                  end do
            end do
            close(99,status="keep")

      case ('g')
            !print *, 'retrieving cb... '//myway(5:6)
            open(99,file='xamfi-'//myway(5:6),status='old',form='unformatted',     &
            access='sequential',action='readwrite',position='rewind')

            read(99)    nz, nzt,ntbas(0),nrows,ncols,ntbas(1),ntbas(2),&
                        n2bastq, n2bbasxq, ncmot, n2tmt
                  ncmot = 0
                  do i= 1, nr_fsym
                        read(99) &
                              npsh(i),&
                              nesh(i),&
                              norb(i),&
                              nfbas(i,0),&
                              nfbas(i,1),&
                              nfbas(i,2),&
                              ntmo(i),&
                              i2tmt(i),&
                              i2tmot(i)
                        do j = 1, nr_fsym
                              read(99)&
                              i2basx(i,j)
                        end do
                        ICMO(i)  = NCMOT
                        ICMOQ(i) = NCMOT*NZ
                        NCMO(i)  = NFBAS(i,0)*NORB(i)
                        NCMOQ(i) = NCMO(i)*NZ
                        NCMOT    = NCMOT + NCMO(i)
                  end do
                  ncmotq = ncmot * nz
                  !print *, ' restored values are ',nz, nzt,ntbas(0),nrows,ncols,ntbas(1),ntbas(2),&
                  !n2bastq, n2bbasxq, ncmot, n2tmt
            close(99,status="keep")
            if(mdirac) nzt = nz ! default: 1

      end select

  end subroutine xamfi_set_cb
!**********************************************************************
end module
