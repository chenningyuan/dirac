module interface_for_pammag

  implicit none

  public common_blocks_vars

  save

  integer, public :: int1_unit
  integer, public :: fcklon_unit = 70
  integer, public :: print_prp_level
  integer, public :: nr_ao_total
  integer, public :: nr_ao_large
  integer, public :: nr_ao_small
  integer, public :: t_n2bbasxq
  integer, public :: nr_qu
  integer, public :: t_norb(1:2)
  integer, public :: t_norbt
  integer, public :: t_n2orbxq
  integer, public :: nr_fermion_irr
  integer, public :: nr_boson_irr
  integer, public :: nfbas_fi(1:2, 0:2)
  integer, public :: cmoq_off(1:2)
  integer, public :: basx_off(1:2, 1:2)
  integer, public :: orbx_off(1:2, 1:2)
  integer, public :: prp_cb_pq_to_uq(1:4, 0:7)
  integer, public :: prp_cb_uq_to_pq(1:4, 0:7)
  integer, public :: t_jbtof(0:7, 1:2)
  integer, public :: ip_conmat(1:3)
  integer, public :: nmr_int

  integer, parameter, public :: max_prp_operators = 400
  integer, public :: prp_time_sym(1:max_prp_operators)
  integer, public :: prp_boson_sym(1:max_prp_operators) 

  logical, public :: dopct_rntlon = .false.
  logical, public :: sum_ep_rn

  private

contains

  subroutine common_blocks_vars()
!-----------------------------------------------------------------------
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
#include "dummy.h"
#include "dcbfir.h"
#include "dcbham.h"
!-----------------------------------------------------------------------
    integer :: i, j
!-----------------------------------------------------------------------

    int1_unit          = lu1int
    fcklon_unit        = lufckl
    dopct_rntlon       = dopct
    sum_ep_rn          = doeprn

    print_prp_level    = iprprp 

! atomic orbitals:
    nr_ao_total        = ntbas(0)
    nr_ao_large        = ntbas(1)
    nr_ao_small        = ntbas(2)
    t_n2bbasxq         = n2bbasxq  ! = nr_ao_total*nr_ao_total*nr_qu

    nr_qu              = nz

! molecular orbitals:
    t_norb(1:2)        = norb(1:2)
    t_norbt            = norbt
    t_n2orbxq          = n2orbxq

!
    nr_fermion_irr     = nfsym
    nr_boson_irr       = nbsym
    
!
    nfbas_fi = 0
    do i = 1, 2
      do j = 0, 2
        nfbas_fi(i, j)  = nfbas(i, j)
      end do
    end do

    !nr_boson_irr(2, 0:2)    = nbbas(2, 0:2)
    
!icmoq(2):
    cmoq_off(1:2) = icmoq(1:2)

!i2basx(2, 2) & i2orbx(2, 2)
    basx_off(1:2, 1:2) = i2basx(1:2, 1:2)
    orbx_off(1:2, 1:2) = i2orbx(1:2, 1:2)


! ipqtoq & iqtopq
    prp_cb_pq_to_uq = 0
    prp_cb_uq_to_pq = 0
    do i = 1, 4
      do j = 0, 7 
        prp_cb_pq_to_uq(i, j) = ipqtoq(i, j)
        prp_cb_uq_to_pq(i, j) = iqtopq(i, j)
      end do
    end do

! jbtof(0:7, 2)
    t_jbtof = 0
    do j = 1, 2
      do i = 0, 7 
        t_jbtof(i, j) = jbtof(i, j)
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
    prp_time_sym = 0
    prp_boson_sym = 0
    do i = 1, max_prp_operators
    !do i = 1, maxprps
      prp_time_sym(i)  = iprptim(i) !(+1,-1)
      prp_boson_sym(i) = iprpsym(i) !(1 - 8)
    end do

  end subroutine

end module
