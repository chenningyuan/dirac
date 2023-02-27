!dirac_copyright_start
!      Copyright (c) 2019 by the authors of DIRAC.
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
! module containing the functionality to renew AOMOMAT as final step of the 
! X2Cmmf procedure
!
! written by sknecht july 2020
!
module x2c_renew_AOMOMAT

  use common_matvec_op, only: insert_elm_a_b

  implicit none

  public x2c_molecular_mean_field_AOMOMAT

contains

!**********************************************************************
  subroutine x2c_molecular_mean_field_AOMOMAT(&
                  fh    ,                     &
                  tm1   ,                     &
                  tm2   ,                     &
                  tmX   ,                     &
                  lr1   ,                     &
                  lc1   ,                     &
                  lr2   ,                     &
                  lc2   ,                     &
                  lrX   ,                     &
                  lcX   ,                     &
                  nfsym ,                     &
                  nzt   ,                     &
                  i2tmt ,                     &
                  n2tmt ,                     &
                  lupri ,                     &
                  debug                       &
                  )  
!**********************************************************************
!
!    purpose: renew AOMOMAT.
!             --> this module is part of the molecular-mean-field option within 
!             the X2C framework.
!
!----------------------------------------------------------------------
     integer, intent(in)             :: nzt
     integer, intent(in)             :: nfsym
     integer, intent(in)             :: lr1(1:nfsym)
     integer, intent(in)             :: lc1(1:nfsym)
     integer, intent(in)             :: lr2(1:nfsym)
     integer, intent(in)             :: lc2(1:nfsym)
     integer, intent(in)             :: lrX(1:nfsym)
     integer, intent(in)             :: lcX(1:nfsym)
     integer, intent(inout)          :: i2tmt(1:nfsym)
     integer, INTENT(INOUT)          :: n2tmt
     integer, intent(in)             :: fh
     integer, intent(in)             :: lupri
     logical, intent(in)             :: debug
     real*8 , INTENT(IN)             :: tm1(*)
     real*8 , INTENT(IN)             :: tm2(*)
     real*8 , INTENT(INOUT)          :: tmX(*)     
!----------------------------------------------------------------------
     integer                         :: i, j, k, l
!**********************************************************************
      
      !> TM1
      j = 1; k = 1
      do i = 1, nfsym
        call insert_elm_a_b(tm1(j),lr1(i),lc1(i),nzt,tmX(k),lrX(i),lcX(i),nzt)
        j = j + lr1(i)*lc1(i)*nzt
        k = k + lrX(i)*lcX(i)*nzt
      end do
      if(debug) &
          write(lupri,*) ' fh is ',fh, ' dim is ', k-1, 'first 10 elms are ',tmX(1:10)
      CALL WRITT(fh,k-1,tmX)

      !> TM2
      j = 1; k = 1
      do i = 1, nfsym
        call insert_elm_a_b(tm2(j),lr2(i),lc2(i),nzt,tmX(k),lrX(i),lcX(i),nzt)
        j = j + lr2(i)*lc2(i)*nzt
        k = k + lrX(i)*lcX(i)*nzt
      end do
      if(debug) &
          write(lupri,*) ' fh is ',fh, ' dim is ', k-1, 'first 10 elms are ',tmX(1:10)
      CALL WRITT(fh,k-1,tmX)
      
      !> renew dimensions
      N2TMT = 0
      do i = 1, nfsym
        I2TMT(I)     = N2TMT
        N2TMT        = N2TMT + lrX(i)*lcX(i)*NZT
      end do

    end subroutine x2c_molecular_mean_field_AOMOMAT
!**********************************************************************

end module x2c_renew_AOMOMAT
