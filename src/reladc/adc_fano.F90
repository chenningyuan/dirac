MODULE adc_fano

CONTAINS

  SUBROUTINE fanoadc(ckks,bufhp,ladc,krep,itapadc,iw)
!
!  Purpose: Calculate decay widths using the FanoADC method.
!           This is the controlling module. All routines are to be found
!           in adc_fano_routines.
!           This module is supposed to be tidied up. No other routines are
!           to be added.
!
!  Author:  Elke Fasshauer

    use adc_cfg
    use adc_fano_routines
    use memory_allocator

    IMPLICIT NONE

!
!  Common Blocks
!
#include "../relccsd/symm.inc"

!
!  Data dictionary: Calling Variables
!
    REAL(8), INTENT(IN)   :: ckks(*), bufhp(*)
    INTEGER, INTENT(IN)   :: krep
    INTEGER, INTENT(IN)   :: ladc
    INTEGER, INTENT(IN)   :: itapadc
    INTEGER, INTENT(IN)   :: iw

!
!  Data dictionary: Local Variables (and to be shared with subroutines
!
    INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: bigicra

!
!  Code
!
    WRITE(iw,*)
    WRITE(iw,*) '*********************************************'
    WRITE(iw,*) '*                 FANOADC                   *'
    WRITE(iw,*) '*       written by Elke Fasshauer           *'
    WRITE(iw,*) '*********************************************'
    WRITE(iw,*) 

    ! Create lookup table bigicra
    CALL alloc(bigicra,9,NVOOT(krep))
    CALL make_bigicra(krep,bigicra,iw)

    ! Create sorted h/2h1p block
    CALL fano_write_hhp(ckks,bufhp,ladc,krep,bigicra,   &
                             itapadc,iw)

    CALL dealloc(bigicra)

  END SUBROUTINE fanoadc

END MODULE adc_fano
