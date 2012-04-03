MODULE ARGOS_ONEINT_NUC
CONTAINS

SUBROUTINE ONEINT_NUC_ISEQJS(is,js,ic,ng,ncc,mcu,msu,ica,icb,nfct,ngw)
! center combinations for is=js, Y.C.Park
! ica(NUC_I/J, JS/IS, SYM/SYMUNQ)
! icb(NUC_I/J, JS/IS, SYM/SYMUNQ)
!        ^       ^       ^
  USE RECP_OUTPUT
  IMPLICIT NONE
#include "inc_print.h"
! vaiables : global
  INTEGER :: is             ! DO is->ns
  INTEGER :: js             ! DO js->is
  INTEGER :: ic             ! CENTER_I nc(is)
  INTEGER :: ng             ! SYM_N : symmetry_number from ngen?
  INTEGER :: ncc            ! number of symmetry unique center?
  INTEGER :: mcu, msu       ! from common(parmi)
  INTEGER :: ica(mcu,msu,*)
  INTEGER :: icb(4,24,*)
  INTEGER :: nfct(*)
  INTEGER :: ngw(*)
! vaiables : local
  INTEGER :: NUC_J          ! NUC_CENTER_J (jc)
  INTEGER :: SYM_I          ! SYM_I : symmetry_ith (ig)
  INTEGER :: ICA_MIN
  INTEGER :: igw, igwu
  INTEGER :: NUC_SYMUNQ_I   ! NUC_CENTER_SYM_UNIQUE_I (icc)
  INTEGER :: TEMP

      CALL OUTPUT_LOC('ONEINT_CENTER_ISEQJS','E')
      do 140 NUC_J = 1, ic
      igwu=0
      do 132 SYM_I = 1, ng
        if(ica(ic,is,SYM_I)   .ne.ic .and. &
           ica(NUC_J,is,SYM_I).ne.ic) goto 132
        ICA_MIN=min(ica(ic,is,SYM_I),ica(NUC_J,is,SYM_I))
 
        do 112 NUC_SYMUNQ_I=ncc,1,-1
          if(ICA_MIN-icb(2,1,NUC_SYMUNQ_I)) 112,136,116
  112   continue
 
  116   if(ICA_MIN.ne.NUC_J) go to 132
 
        do igw=igwu,1,-1
          if(ica(ic,is,SYM_I)   .eq.icb(1,igw,ncc+1).and. &
             ica(NUC_J,is,SYM_I).eq.icb(2,igw,ncc+1)) goto 132
        enddo
 
        igwu=igwu+1
        icb(1,igwu,ncc+1)=ica(ic,is,SYM_I)
        icb(2,igwu,ncc+1)=ica(NUC_J,is,SYM_I)
  132 continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc)=igwu
      go to 140
  136 nfct(NUC_SYMUNQ_I) = nfct(NUC_SYMUNQ_I) + ic
  140 continue

      do NUC_SYMUNQ_I = 1, ncc-1
        nfct(NUC_SYMUNQ_I) = nfct(NUC_SYMUNQ_I)/2
      enddo

      CALL OUTPUT_LOC('ONEINT_CENTER_ISEQJS','X')
END SUBROUTINE ONEINT_NUC_ISEQJS


SUBROUTINE ONEINT_NUC_ISGTJS(is,js,nc,ng,ic,mcu,msu,ica,icb,ncc,nfct)
! Center combinations for is>js
! Y. C. Park
!
! ica(CENTER_I/J, JS/IS, SYM/SYMUNQ)
! icb(CENTER_I/J, JS/IS, SYM/SYMUNQ)
!        ^          ^       ^
  USE RECP_OUTPUT
  IMPLICIT NONE
#include "inc_print.h"
! vaiables : global
  INTEGER :: is
  INTEGER :: js
  INTEGER :: nc(*)
  INTEGER :: ng             ! SYM_N : symmetry_number from ngen?
  INTEGER :: ic
  INTEGER :: mcu, msu       ! from common(parmi)
  INTEGER :: ica(mcu,msu,*)
  INTEGER :: icb(4,24,*)
  INTEGER :: ncc            ! number of symmetry unique center?
  INTEGER :: nfct(*)
! vaiables : local
  INTEGER :: NUC_J         ! CENTER_J (jc)
  INTEGER :: SYM_I          ! SYM_I : symmetry_ith (ig)
  INTEGER :: NUC_SYMUNQ_I  ! CENT_SYM_UNIQUE_I (icc)

      CALL OUTPUT_LOC('ONEINT_CENTER_ISGTJS','E')

      do 240 NUC_J = 1, nc(js)
      do 232 SYM_I = 1, ng
         if(ica(ic,is,SYM_I) .eq. ic) then
           do 212 NUC_SYMUNQ_I=ncc,1,-1
              if(ica(NUC_J,js,SYM_I) - icb(2,1,NUC_SYMUNQ_I)) 212,236,232
  212      continue
         endif
  232 continue
      ncc = ncc + 1
      icb(1,1,ncc) = ic
      icb(2,1,ncc) = NUC_J
      nfct(ncc)    = ic
      goto 240
  236 nfct(NUC_SYMUNQ_I) = nfct(NUC_SYMUNQ_I) + ic
  240 continue

      CALL OUTPUT_LOC('ONEINT_CENTER_ISGTJS','X')

END SUBROUTINE ONEINT_NUC_ISGTJS

END MODULE ARGOS_ONEINT_NUC
