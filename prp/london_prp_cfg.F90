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

module london_prp_cfg

      use quaternion_algebra
#ifdef VAR_MPI
  use dirac_mpi
#endif

  implicit none

  public common_blocks_vars
  public broadcast_prp_lao_variables


  save


! parameters:
  integer, parameter, public :: fcklon_unit = 70
  integer, parameter, public :: twolon_unit = 69
  integer, parameter, public :: p_lupri = 6
  integer, parameter, public :: max_prp_operators = 400
  integer, parameter, public :: mxqn = 10
  integer, parameter, public :: mxaqn = mxqn*(mxqn+1)/2 
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (INT_STAR8)
  integer, parameter, public :: p_irat = 1
#else
  integer, parameter, public :: p_irat = 2
#endif

  integer, public :: int1_unit !should also be a parameter
  integer, public :: prp_iprint
  integer, public :: nr_ao_total
  integer, public :: p_n2bbasxq
  integer, public :: p_n2bbasx
  integer, public :: nr_qu
  integer, public :: p_norb(1:2)
  integer, public :: p_norbt
  integer, public :: p_n2orbxq
  integer, public :: p_nfsym 
  integer, public :: p_nbsym
  integer, public :: p_nfbas(1:2, 0:2)
  integer, public :: p_icmoq(1:2)
  integer, public :: p_i2basx(1:2, 1:2)
  integer, public :: p_i2orbx(1:2, 1:2)
  integer, public :: p_ipqtoq(1:4, 0:7)
  integer, public :: p_iqtopq(1:4, 0:7)
  integer, public :: p_jbtof(0:7, 1:2)
  integer, public :: ip_conmat(1:3)
  integer, public :: nmr_int
  integer, public :: p_iprptim(1:max_prp_operators)
  integer, public :: p_iprpsym(1:max_prp_operators) 
  logical, public :: dopct_rntlon = .false.
  logical, public :: sum_ep_rn
  integer, public :: p_maxrep
  integer, public :: p_kmax
  integer, public :: p_khk(mxqn)
  integer, public :: p_khkt(mxqn)
  integer, public :: p_nhkt(mxqn)
  integer, public :: p_nhtyp
  integer, public :: p_nlrgbl
  integer, public :: p_nsmlbl
  logical, public :: p_first1
  logical, public :: p_first2
  logical, public :: p_first3
  integer, public :: p_iprnta
  integer, public :: p_iprntb
  integer, public :: p_iprntc
  integer, public :: p_iprntd
  logical, public :: p_tktime
  integer, public :: p_isymax(3, 2)
  integer, public :: p_irqmat(4, 0:7)
  integer, public :: p_iqmult(4, 4, 4)
  integer, public :: p_jqbas(0:7, 2)
  integer, public :: p_ihqmat(4, -1:1)
  integer, public :: p_ipar(0:7)

  integer, parameter, public :: p_noptyp = 21
  integer, public :: p_jm4(3, p_noptyp)
  integer, public :: p_jm4pos(0:7)

! new variables
  integer, public :: nr_xyz_comp
  logical, public :: p_levyle
!
  integer, public :: p_nlrgsh




  private
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "dcbfir.h"
#include "dcbham.h"
#include "ccom.h"
#include "iratdef.h"
#include "shells.h"
#include "symmet.h"
#include "cbihr2.h"



#include "infpar.h"
#ifdef VAR_MPI
#include "mpif.h"
  dimension     istat(mpi_status_size)
#endif /* ifdef VAR_MPI */


contains


  subroutine common_blocks_vars()
!-----------------------------------------------------------------------
    integer :: i, j, k
!-----------------------------------------------------------------------

    int1_unit          = lu1int
    dopct_rntlon       = dopct
    sum_ep_rn          = doeprn
    prp_iprint         = iprprp 

! atomic orbitals:
    nr_ao_total        = ntbas(0)

    p_n2bbasxq         = n2bbasxq  ! = nr_ao_total*nr_ao_total*nr_qu
    p_n2bbasx          = n2bbasx   ! = nr_ao_total*nr_ao_total

    nr_qu              = nz

! molecular orbitals:
    p_norb(1:2)        = norb(1:2)
    p_norbt            = norbt
    p_n2orbxq          = n2orbxq

    p_nfsym            = nfsym
    p_nbsym            = nbsym
    
!
    p_nfbas = 0
    do i = 1, 2
      do j = 0, 2
        p_nfbas(i, j)  = nfbas(i, j)
      end do
    end do

!icmoq(2):
    p_icmoq(1:2) = icmoq(1:2)

!i2basx(2, 2) & i2orbx(2, 2)
    p_i2basx(1:2, 1:2) = i2basx(1:2, 1:2)
    p_i2orbx(1:2, 1:2) = i2orbx(1:2, 1:2)

! ipqtoq & iqtopq
    p_ipqtoq = 0
    p_iqtopq = 0
    do i = 1, 4
      do j = 0, 7 
        p_ipqtoq(i, j) = ipqtoq(i, j)
        p_iqtopq(i, j) = iqtopq(i, j)
      end do
    end do

! jbtof(0:7, 2)
    p_jbtof = 0
    do j = 1, 2
      do i = 0, 7 
        p_jbtof(i, j) = jbtof(i, j)
      end do
    end do

! connection matrices:
    ip_conmat = 0
    do i = 1, 3
      ip_conmat(i)     = ipcon(i)
    end do

! intnmr = illint + 2*islint + 4*issint
    nmr_int            = intnmr

!iprptim, iprpsym:
    p_iprptim = 0
    p_iprpsym = 0
    do i = 1, max_prp_operators
      p_iprptim(i)  = iprptim(i)     !(+1,-1)
      p_iprpsym(i) = iprpsym(i)      !(1 - 8)
    end do

    p_maxrep = maxrep
    p_kmax   = kmax
    p_nhtyp  = nhtyp
    do i = 1, mxqn
      p_khk(i)  = khk(i)
      p_khkt(i) = khkt(i)
      p_nhkt(i) = nhkt(i)
    end do

    p_nlrgbl = nlrgbl
    p_nsmlbl = nsmlbl

    p_first1 = first1
    p_first2 = first2
    p_first3 = first3
 

    p_iprnta = iprnta
    p_iprntb = iprntb
    p_iprntc = iprntc
    p_iprntd = iprntd
  
    p_tktime = tktime

    p_isymax = 0
    do i = 1, 3
      do j = 1, 2
        p_isymax(i, j) = isymax(i, j)
      end do
    end do

    p_irqmat = 0
    do i = 1, 4
      do j = 0, 7
        p_irqmat(i, j) = irqmat(i, j)
      end do
    end do

    p_iqmult = 0
    do i = 1, 4
      do j = 1, 4
        do k = 1, 4
          p_iqmult(i, j, k) = iqmult(i, j, k)
        end do
      end do
    end do 

    p_jqbas = 0
    do i = 0, 7
      do j = 1, 2
        p_jqbas(i, j) = jqbas(i, j)
      end do
    end do

    p_ihqmat = 0
    do i = 1, 4
      do j = -1, 1
        p_ihqmat(i, j) = ihqmat(i, j)
      end do
    end do

    if (bxlao .or. bylao .or. bzlao) then
      nr_xyz_comp = 1
    else
      nr_xyz_comp = 3
    end if


!    do i = 0, 7
    do i = 0, maxrep
      p_ipar(i) = ipar(i)
    end do

    p_nlrgsh = nlrgsh
    p_levyle = levyle


    p_jm4 = 0
    do i = 1, p_noptyp
      do j = 1, 3
        p_jm4(j, i) = jm4(j, i)
      end do
    end do

    p_jm4pos = 0
    do i = 0, 7
      p_jm4pos(i) = jm4pos(i)
    end do


  end subroutine


  subroutine broadcast_prp_lao_variables()

#ifdef VAR_MPI
    if (parcal) then
      call dirac_mpi_bcast(int1_unit   )
      call dirac_mpi_bcast(prp_iprint  )
      call dirac_mpi_bcast(nr_ao_total )
      call dirac_mpi_bcast(p_n2bbasxq  )
      call dirac_mpi_bcast(p_n2bbasx   )
      call dirac_mpi_bcast(nr_qu       )
      call dirac_mpi_bcast(p_norb      )
      call dirac_mpi_bcast(p_norbt     )
      call dirac_mpi_bcast(p_n2orbxq   )
      call dirac_mpi_bcast(p_nfsym     )
      call dirac_mpi_bcast(p_nbsym     )
      call dirac_mpi_bcast(p_nfbas     )
      call dirac_mpi_bcast(p_icmoq     )
      call dirac_mpi_bcast(p_i2basx    )
      call dirac_mpi_bcast(p_i2orbx    )
      call dirac_mpi_bcast(p_ipqtoq    )
      call dirac_mpi_bcast(p_iqtopq    )
      call dirac_mpi_bcast(p_jbtof     )
      call dirac_mpi_bcast(ip_conmat   )
      call dirac_mpi_bcast(nmr_int     )
      call dirac_mpi_bcast(p_iprptim   )
      call dirac_mpi_bcast(p_iprpsym   )
      call dirac_mpi_bcast(dopct_rntlon)
      call dirac_mpi_bcast(sum_ep_rn   )
      call dirac_mpi_bcast(p_maxrep    )
      call dirac_mpi_bcast(p_kmax      )
      call dirac_mpi_bcast(p_khk       )
      call dirac_mpi_bcast(p_khkt      )
      call dirac_mpi_bcast(p_nhkt      )
      call dirac_mpi_bcast(p_nhtyp     )
      call dirac_mpi_bcast(p_nlrgbl    )
      call dirac_mpi_bcast(p_nsmlbl    )
      call dirac_mpi_bcast(p_first1    )
      call dirac_mpi_bcast(p_first2    )
      call dirac_mpi_bcast(p_first3    )
      call dirac_mpi_bcast(p_iprnta    )
      call dirac_mpi_bcast(p_iprntb    )
      call dirac_mpi_bcast(p_iprntc    )
      call dirac_mpi_bcast(p_iprntd    )
      call dirac_mpi_bcast(p_tktime    )
      call dirac_mpi_bcast(p_isymax    )
      call dirac_mpi_bcast(p_irqmat    )
      call dirac_mpi_bcast(p_iqmult    )
      call dirac_mpi_bcast(p_jqbas     )
      call dirac_mpi_bcast(p_ihqmat    )
      call dirac_mpi_bcast(p_ipar      )
      call dirac_mpi_bcast(p_jm4       )
      call dirac_mpi_bcast(p_jm4pos    )
      call dirac_mpi_bcast(nr_xyz_comp )
      call dirac_mpi_bcast(p_levyle    )
      call dirac_mpi_bcast(p_nlrgsh    )
!parameters:
      call dirac_mpi_bcast(p_noptyp    )
      call dirac_mpi_bcast(fcklon_unit )
      call dirac_mpi_bcast(twolon_unit )
      call dirac_mpi_bcast(p_lupri     )
      call dirac_mpi_bcast(max_prp_operators)
      call dirac_mpi_bcast(mxqn        )
      call dirac_mpi_bcast(mxaqn       )
      call dirac_mpi_bcast(p_irat      )
    end if
#endif

  end subroutine

end module
