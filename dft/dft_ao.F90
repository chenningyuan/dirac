module dft_ao

  use dft_cb_interface

  implicit none


contains

         subroutine dft_getsos(gao, gao1, &
                         corx, cory, corz, &
                         buf, ndim, &
                         nder)

      real(8) :: GAO(ndim), GAO1(ndim,3), &
     &          BUF(ndim)
      real(8) :: corx, cory, corz
      integer :: ndim
      integer :: nder, j
!
!     First get values for the Ao-basis
      CALL dft_GETAOS(GAO,GAO1,CORX,CORY,CORZ, &
     &            ndim,NDER)

      if (nr_boson_sym == 1) return
!
!     Functions
      CALL WTAOSO(GAO,BUF,ndim,1,0)
      CALL VTBUBS(BUF,GAO,1)
!
!     First derivatives
      IF(NDER.GE.1) THEN
        DO J = 1,3
          CALL WTAOSO(GAO1(1,J),BUF,ndim,1,0)
          CALL VTBUBS(BUF,GAO1(1,J),1)
        ENDDO
      ENDIF

      end subroutine



      SUBROUTINE dft_GETAOS(GAO,GAO1,CORPX,CORPY,CORPZ, &
     &                  ndim,NDER)
!
!     T. Helgaker sep 99
!
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
!   
      real(8) :: GAO(ndim), GAO1(ndim,3)


      integer :: ishbegin, ishend, iadr, ndim, nder
      integer :: isymop, ishela
      real(8) :: pax, pay, paz, corpx, corpy,corpz

!
#include "onecom.h"
#include "lmns.h"
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
!
!
!
      gao = 0.0d0
      if (nder > 0) then
        gao1 = 0.0d0
      end if
      IADR = 1
       ISHBEGIN=1
       ISHEND  =KMAX

!     if Levy-Leblond or IOTC
!     don't work on small component shells
      IF (nr_ls_components .EQ. 1) ISHEND = NLRGSH

      DO 100 ISHELA = ISHBEGIN,ISHEND
         NHKTA = NHKT(ISHELA)
         KHKTA = KHKT(ISHELA)
         KCKTA = KCKT(ISHELA)
         SPHRA = SPHR(ISHELA)
         IF (NDER.GT.1 .OR. NHKTA.GT.2) THEN
            CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)
         END IF
         NUCA   = NUCO(ISHELA)
         NUMCFA = NUMCF(ISHELA)
         JSTA   = JSTRT(ISHELA)
!
!tsaue: Loop over symmetry dependent centers   
!
         MULA   = ISTBAO(ISHELA)
         DO ISYMOP = 0,MAXOPR
         IF (iand(ISYMOP,MULA) .EQ. 0) THEN
           PAX = CORPX-PT(iand(ISYMAX(1,1),ISYMOP))*CENT(ISHELA,1,1)
           PAY = CORPY-PT(iand(ISYMAX(2,1),ISYMOP))*CENT(ISHELA,2,1)
           PAZ = CORPZ-PT(iand(ISYMAX(3,1),ISYMOP))*CENT(ISHELA,3,1)
           IF (NDER.EQ.0) CALL dft_GETGAO(GAO(IADR),PAX,PAY,PAZ)
           IF (NDER.GT.0) THEN
             CALL dft_GETGA1(GAO(IADR), &
     &                   GAO1(IADR,1),GAO1(IADR,2),GAO1(IADR,3), &
     &                   PAX,PAY,PAZ)
           END IF
           IADR = IADR + KHKTA
        ENDIF
        ENDDO
  100 CONTINUE
      end subroutine
      SUBROUTINE dft_GETGAO(GAO,PAX,PAY,PAZ)
!
!     T. Helgaker Sep 99
!
#include "aovec.h"
#include "maxorb.h"
#include "maxaqn.h"
!
      real(8) :: GAO(KCKTA), CAO(MXAQN)
      integer :: ISTR, IEND, IPRIMA, ICOMPA
      real(8) :: pax, pay, paz, pa2, ga
!
#include "lmns.h"
#include "onecom.h"
#include "primit.h"
#include "sphtrm.h"
!
      ISTR = JSTA + 1
      IEND = JSTA + NUCA
      PA2 = PAX**2 + PAY**2 + PAZ**2
      IF (KCKTA .EQ. 1) THEN
         GAO(1) = PRICCF(ISTR,NUMCFA)*DEXP(-PRIEXP(ISTR)*PA2)
         DO 100 IPRIMA = ISTR + 1, IEND
            GAO(1) = GAO(1) + PRICCF(IPRIMA,NUMCFA) &
     &                   *DEXP(-PRIEXP(IPRIMA)*PA2)
  100    CONTINUE
      ELSE IF (KCKTA .EQ. 3) THEN
         GA = PRICCF(ISTR,NUMCFA)*DEXP(-PRIEXP(ISTR)*PA2)
         GAO(1) = PAX*GA
         GAO(2) = PAY*GA
         GAO(3) = PAZ*GA
         DO 200 IPRIMA = ISTR + 1, IEND
            GA = PRICCF(IPRIMA,NUMCFA)*DEXP(-PRIEXP(IPRIMA)*PA2)
            GAO(1) = GAO(1) + PAX*GA
            GAO(2) = GAO(2) + PAY*GA
            GAO(3) = GAO(3) + PAZ*GA
  200    CONTINUE
      ELSE
         IF (SPHRA) THEN
            GA = PRICCF(ISTR,NUMCFA)*DEXP(-PRIEXP(ISTR)*PA2)
            DO 300 ICOMPA = 1,KCKTA
               CAO(ICOMPA) = (PAX**LVALUA(ICOMPA))* &
     &                       (PAY**MVALUA(ICOMPA))* &
     &                       (PAZ**NVALUA(ICOMPA))*GA
  300       CONTINUE
            DO 400 IPRIMA = ISTR + 1, IEND
               GA = PRICCF(IPRIMA,NUMCFA)*DEXP(-PRIEXP(IPRIMA)*PA2)
               DO 500 ICOMPA = 1,KCKTA
                  CAO(ICOMPA) = CAO(ICOMPA) + (PAX**LVALUA(ICOMPA))* &
     &                                        (PAY**MVALUA(ICOMPA))* &
     &                                        (PAZ**NVALUA(ICOMPA))*GA
  500          CONTINUE
  400       CONTINUE
            CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAO,1,KCKTA,GAO,1, &
     &                KHKTA,KHKTA,KCKTA,1)
         ELSE
            GA = PRICCF(ISTR,NUMCFA)*DEXP(-PRIEXP(ISTR)*PA2)
            DO 600 ICOMPA = 1,KCKTA
               GAO(ICOMPA) = (PAX**LVALUA(ICOMPA))* &
     &                       (PAY**MVALUA(ICOMPA))* &
     &                       (PAZ**NVALUA(ICOMPA))*GA
  600       CONTINUE
            DO 700 IPRIMA = ISTR + 1, IEND
               GA = PRICCF(IPRIMA,NUMCFA)*DEXP(-PRIEXP(IPRIMA)*PA2)
               DO 800 ICOMPA = 1,KCKTA
                  GAO(ICOMPA) = GAO(ICOMPA) + (PAX**LVALUA(ICOMPA))* &
     &                                        (PAY**MVALUA(ICOMPA))* &
     &                                        (PAZ**NVALUA(ICOMPA))*GA
  800          CONTINUE
  700       CONTINUE
         END IF
      END IF
      end subroutine
      SUBROUTINE dft_GETGA1(GAO,GAOX,GAOY,GAOZ,PAX,PAY,PAZ)
!
!     T. Helgaker Sep 99
!
#include "aovec.h"
#include "maxorb.h"
#include "maxaqn.h"
!
      LOGICAL FIRST
      real(8) :: GAO(KCKTA), GAOX(KCKTA), GAOY(KCKTA), GAOZ(KCKTA), &
     &          CAO(MXAQN), CAOX(MXAQN), CAOY(MXAQN), CAOZ(MXAQN)
      integer :: l, m, n, ICOMPA, iprima
      real(8) :: TGA, tgx, tgy, tgz
      real(8) :: pax, pay, paz ,pa2, p0
      real(8) :: gax, gay, gaz ,ga
      real(8) :: alpha, talph, tp0
!
#include "lmns.h"
#include "onecom.h"
#include "primit.h"
#include "sphtrm.h"
!
!
      PA2 = PAX**2 + PAY**2 + PAZ**2
      FIRST = .TRUE.
      DO 100 IPRIMA = JSTA + 1, JSTA + NUCA
         ALPHA = PRIEXP(IPRIMA)
         TALPH = -2.0d0*ALPHA
         GA = PRICCF(IPRIMA,NUMCFA)*DEXP(-ALPHA*PA2)
         IF (KCKTA .EQ. 1) THEN
            TGA = TALPH*GA
            IF (FIRST) THEN
               GAO(1)  = GA
               GAOX(1) = PAX*TGA
               GAOY(1) = PAY*TGA
               GAOZ(1) = PAZ*TGA
            ELSE
               GAO(1)  = GAO(1)  + GA
               GAOX(1) = GAOX(1) + PAX*TGA
               GAOY(1) = GAOY(1) + PAY*TGA
               GAOZ(1) = GAOZ(1) + PAZ*TGA
            END IF
         ELSE IF (KCKTA .EQ. 3) THEN
            TGA = TALPH*GA 
            TGX = PAX*TGA 
            TGY = PAY*TGA 
            TGZ = PAZ*TGA
            IF (FIRST) THEN
               GAO(1)  = PAX*GA
               GAO(2)  = PAY*GA
               GAO(3)  = PAZ*GA
               GAOX(1) = PAX*TGX + GA
               GAOX(2) = PAY*TGX
               GAOX(3) = PAZ*TGX
               GAOY(1) = PAX*TGY
               GAOY(2) = PAY*TGY + GA
               GAOY(3) = PAZ*TGY
               GAOZ(1) = PAX*TGZ
               GAOZ(2) = PAY*TGZ
               GAOZ(3) = PAZ*TGZ + GA
            ELSE
               GAO(1)  = GAO(1)  + PAX*GA
               GAO(2)  = GAO(2)  + PAY*GA
               GAO(3)  = GAO(3)  + PAZ*GA
               GAOX(1) = GAOX(1) + PAX*TGX + GA
               GAOX(2) = GAOX(2) + PAY*TGX
               GAOX(3) = GAOX(3) + PAZ*TGX
               GAOY(1) = GAOY(1) + PAX*TGY
               GAOY(2) = GAOY(2) + PAY*TGY + GA
               GAOY(3) = GAOY(3) + PAZ*TGY
               GAOZ(1) = GAOZ(1) + PAX*TGZ
               GAOZ(2) = GAOZ(2) + PAY*TGZ
               GAOZ(3) = GAOZ(3) + PAZ*TGZ + GA
            END IF
         ELSE
            DO 200 ICOMPA = 1,KCKTA
               L = LVALUA(ICOMPA)
               M = MVALUA(ICOMPA)
               N = NVALUA(ICOMPA)
               P0  = (PAX**L)*(PAY**M)*(PAZ**N)
!
!              up contributions 
!
               TP0 = TALPH*P0
               GAX = PAX*TP0
               GAY = PAY*TP0
               GAZ = PAZ*TP0
!
!              down contributions
!
               IF(L.GT.0)GAX=GAX+dfloat(L)*(PAX**(L-1))*(PAY**M)*(PAZ**N)
               IF(M.GT.0)GAY=GAY+dfloat(M)*(PAX**L)*(PAY**(M-1))*(PAZ**N)
               IF(N.GT.0)GAZ=GAZ+dfloat(N)*(PAX**L)*(PAY**M)*(PAZ**(N-1))
!
               IF (SPHRA) THEN
                  IF (FIRST) THEN
                     CAO (ICOMPA) =  P0*GA
                     CAOX(ICOMPA) = GAX*GA
                     CAOY(ICOMPA) = GAY*GA
                     CAOZ(ICOMPA) = GAZ*GA
                  ELSE
                     CAO (ICOMPA) = CAO (ICOMPA) +  P0*GA
                     CAOX(ICOMPA) = CAOX(ICOMPA) + GAX*GA
                     CAOY(ICOMPA) = CAOY(ICOMPA) + GAY*GA
                     CAOZ(ICOMPA) = CAOZ(ICOMPA) + GAZ*GA
                  END IF
               ELSE
                  IF (FIRST) THEN
                     GAO (ICOMPA) =   P0*GA
                     GAOX(ICOMPA) =  GAX*GA
                     GAOY(ICOMPA) =  GAY*GA
                     GAOZ(ICOMPA) =  GAZ*GA
                  ELSE
                     GAO (ICOMPA) = GAO (ICOMPA) +  P0*GA
                     GAOX(ICOMPA) = GAOX(ICOMPA) + GAX*GA
                     GAOY(ICOMPA) = GAOY(ICOMPA) + GAY*GA
                     GAOZ(ICOMPA) = GAOZ(ICOMPA) + GAZ*GA
                  END IF
               END IF
  200       CONTINUE
         END IF
         FIRST = .FALSE.
  100 CONTINUE
      IF (SPHRA) THEN
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAO ,1,KCKTA,GAO ,1, &
     &              KHKTA,KHKTA,KCKTA,1)
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAOX,1,KCKTA,GAOX,1, &
     &              KHKTA,KHKTA,KCKTA,1)
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAOY,1,KCKTA,GAOY,1, &
     &              KHKTA,KHKTA,KCKTA,1)
          CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,CAOZ,1,KCKTA,GAOZ,1, &
     &              KHKTA,KHKTA,KCKTA,1)
      END IF
      end subroutine
end module
