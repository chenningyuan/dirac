* Set KR array according to number of functions in fermion
* symmetries.
*
      subroutine kkr_set(KR,IDIM,ISYDI,NIRR,INVSYM)
*
#include "implicit.inc"
*
      dimension KR(-IDIM:IDIM),ISYDI(NIRR)
*
      NTESTL = 0
*
      ICT = 1
      IV = 1
      do IRR = 1,1+INVSYM,1
        IRRR = IRR
        if (INVSYM.eq.1.and.NIRR.gt.2.and.IRR.eq.2) IRRR = IRR + 1
        do I = 1,ISYDI(IRRR),1
          KR(IV) = ICT
          KR(-IV) = ICT + ISYDI(IRRR)
          ICT = ICT + 1
          IV = IV + 1
        end do
        ICT = ICT + ISYDI(IRR)
      end do
*
      if (NTESTL.ge.1) then
        write(6,*) ' KKR has set KR array :'
        do ID = 1,IDIM,1
          write(6,'(3X,2I4)') KR(-ID),KR(ID)
        end do
      end if
*
      return
      end
