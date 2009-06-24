module dft_ao

! * gao, gao1, ... are not dzeroed for more performance
!   if you change the code make sure that elements of these vectors
!   are not read before they are written for the first time
!   note that in Levy-Leblond or X2C4 mode some elements will
!   contain garbage but they are not (= should not be) used
! * fixme: many things that do not change during a calculation
!   are re-evaluated in each grid point (prefactors, symmetry info)
! * todo: exploit mirror symmetries when evaluating gao1 components

  use dft_cb_interface

  implicit none

  public dft_getsos

  private

contains

  subroutine dft_getsos(gao,    &
                        gao1,   &
                        corx,   &
                        cory,   &
                        corz,   &
                        buffer, &
                        ndim,   &
                        nder)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: gao(ndim)
    real(8), intent(inout) :: gao1(ndim, 3)
    real(8), intent(in)    :: corx
    real(8), intent(in)    :: cory
    real(8), intent(in)    :: corz
    real(8), intent(inout) :: buffer(ndim, 3)
    integer, intent(in)    :: ndim
    integer, intent(in)    :: nder
!   ----------------------------------------------------------------------------
    integer                :: ixyz
!   ----------------------------------------------------------------------------

!   get values in ao basis
    call dft_getaos(gao,  &
                    gao1, &
                    corx, &
                    cory, &
                    corz, &
                    ndim, &
                    nder)

!   without symmetry skip the rest
    if (nr_boson_sym == 1) return

!   transform to so
    call wtaoso(gao, buffer, ndim, 1, 0)
!   transform to boson sorted only if there is a small component
    if (nr_ao_total > nr_ao_large) then
      call vtbubs(buffer, gao, 1)
    else
      call dcopy(ndim, buffer, 1, gao, 1)
    end if

    if (nder == 1) then
      do ixyz = 1, 3
!       transform to so
        call wtaoso(gao1(1, ixyz), buffer, ndim, 1, 0)
!       transform to boson sorted only if there is a small component
        if (nr_ao_total > nr_ao_large) then
          call vtbubs(buffer, gao1(1, ixyz), 1)
        else
          call dcopy(ndim, buffer, 1, gao1(1, ixyz), 1)
        end if
      end do
    end if

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

      IADR = 1
      ISHBEGIN=1
      ISHEND  =KMAX

!     if Levy-Leblond or IOTC
!     don't work on small component shells
!     this will make problems if you calculate the relativistic current density
!     expression using Levy-Leblond, in this case comment this out
      IF (nr_ls_components .EQ. 1) ISHEND = NLRGSH

      DO 100 ISHELA = ISHBEGIN,ISHEND

!        NHKT is l
         NHKTA = NHKT(ISHELA)

!        in principle the degeneracy among spherical Gaussians
!        KHKT = 2*(l - 1) + 1
!        but gives cartesian degeneracy
         KHKTA = KHKT(ISHELA)

!        corresponding cartesian degeneracy
!        KCKT = l*(l + 1)/2
         KCKTA = KCKT(ISHELA)

!        spherical or not (i think)
         SPHRA = SPHR(ISHELA)

         IF (NDER.GT.1 .OR. NHKTA.GT.2) THEN
            CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)
         END IF

!        nr of uncotracted shells (i think)
         NUCA   = NUCO(ISHELA)

!        contracted shell nr
         NUMCFA = NUMCF(ISHELA)

!        offset for the shell
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
