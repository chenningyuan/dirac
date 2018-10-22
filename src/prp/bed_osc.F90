      SUBROUTINE GENOSC(K,FOSC,ATMPPF,OMEGA)
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbxpp.h"
      real*8 :: FOSC(MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      integer, allocatable :: indices(:)
      integer :: j,k,k2,iterm,jterm
      iterm = 0
      jterm = 0
      K2=K+K
      WRITE(6,'(/A,I4)') 'GENOSC: Rotational weights of order: ',K2
      CALL PRSYMB(6,'-',60,0)
      WRITE(6,'(18X,100A1)') (' ',J=1,K2+2),'P','Q'
!.....Accumulate oscillator strengths
      FOSC=0.0D0
      IF(K.EQ.0) THEN
        call rotfac(indices,0,iterm,jterm,FOSC,ATMPPF,OMEGA)
      ELSE
        allocate(indices(k2))
        call bedloop(1,k2,indices,iterm,jterm,FOSC,ATMPPF,OMEGA)
        deallocate(indices)
     ENDIF
!....Add final factors for oscillator strengths
     DO ISYM = 1,NBSYM
       NEXC = KEXCNV(ISYM)
       DO IEXC = 1,NEXC
         FOSC(IEXC,ISYM)=FOSC(IEXC,ISYM)*(OMEGA(IEXC,ISYM)/CVAL)**K2*2.0D0/OMEGA(IEXC,ISYM)
       ENDDO
     ENDDO   
     write(6,*) 'Total/non-zero terms: ',ITERM,JTERM
     CALL PRSYMB(6,'-',60,0)
    END SUBROUTINE GENOSC

      recursive subroutine bedloop(m,k2,indices,iterm,jterm,FOSC,ATMPPF,OMEGA)
#include "dgroup.h"
#include "dcbxpp.h"
      real*8 :: FOSC(MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      integer :: indices(k2),ind(3),k2
      if(m==k2) then
        do i=1,3
           indices(m) = i
           call rotfac(indices,k2,iterm,jterm,FOSC,ATMPPF,OMEGA)
        enddo
        return
      else
        do i=1,3
          indices(m) = i
          call bedloop(m+1,k2,indices,iterm,jterm,FOSC,ATMPPF,OMEGA)
        enddo
      endif 
      END
      
      SUBROUTINE ROTFAC(indices,k2,iterm,jterm,FOSC,ATMPPF,OMEGA)
#include "pgroup.h"
#include "dgroup.h"
#include "dcbxpp.h"
#include "dcbxpr.h"      
      real*8 :: FOSC(MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      integer :: indices(k2),ind(3)
      CHARACTER*10 DENOM
      PARAMETER(D2=2.0D0,D1=1.0D0)
      K=K2/2
      DO IP = 1,3
        DO IQ = 1,3
          iterm = iterm + 1
          IND = 0
          DO L = 1,K2
            J = INDICES(L)
            IND(J) = IND(J) + 1
          ENDDO
          IND(IP) = IND(IP) + 1
          IND(IQ) = IND(IQ) + 1               
          IF(MOD(IND(1),2)==1) CYCLE
          IF(MOD(IND(2),2)==1) CYCLE
! Since M is even, there is no need to test on IND(3)              
          JTERM = JTERM + 1
          WROT=-GAMMA(REAL(IND(1)+1)/2)*GAMMA(REAL(IND(2)+1)/2)*GAMMA(REAL(IND(3)+1)/2)/&
                GAMMA(REAL(IND(1)+IND(2)+IND(3)+3)/2)
          IF(IP.EQ.IQ) THEN
            IND(IP) = IND(IP)-2
            WROT = WROT+GAMMA(REAL(IND(1)+1)/2)*GAMMA(REAL(IND(2)+1)/2)*GAMMA(REAL(IND(3)+1)/2)/&
                        GAMMA(REAL(IND(1)+IND(2)+IND(3)+3)/2)
          ENDIF
          WROT=WROT/4.0D0/ACOS(-1.0D0)
          FAC = GAMMA(REAL(K2+5))/GAMMA(REAL(K2/2+3))/2**(K2/2+1)
          IFAC= NINT(FAC)
          WRITE(DENOM,'(A,I8,A)') '/',IFAC,':'
          WRITE(6,'(F6.1,A10,2X,18A1)') WROT*IFAC,DENOM,&
               (CHAR(87+INDICES(J)),J=1,K2),' ',' ',CHAR(87+IP),CHAR(87+IQ)
          DO M = 0,K
            FAC = D2*(-D1)**M
            IF(M.EQ.0) FAC = D1
!...........Left operator.. order (K+M)
            IORD  = K+M
            INDAP = IBDOFF(IORD) + IP
            IF(IORD.GT.0) THEN
              IND   = 0
              DO L = 1,K+M
                J      = INDICES(L)
                IND(J) = IND(J) + 1
              ENDDO
              IA    = IORD - IND(1)
              IB    = IND(3)
              INDAP = INDAP + 3*(IB+(IA+1)*IA/2)
            ENDIF
            INDPRA  = LPPAPU(INDAP)
            ISYM    = IPRPSYM(INDPRA)
!...........Right operator: order (K-M)
            IORD  = K-M
            INDBP = IBDOFF(IORD) + IQ
            IF(IORD.GT.0) THEN
              IND   = 0
              DO L = K+M+1,K2
                J      = INDICES(L)
                IND(J) = IND(J) + 1
              ENDDO
              IA = IORD - IND(1)
              IB = IND(3)
              INDBP = INDBP + 3*(IB+(IA+1)*IA/2)
            ENDIF
            INDPRB  = LPPAPU(INDBP)
            ISYMB   = IPRPSYM(INDPRB)
            IF(ISYMB.NE.ISYM) CALL QUIT('ROTFAC: Symmetry problem !')
!...........Form product
            WRITE(6,'(2(2X,I3,2X,A16,2X,A3))') K+M,PRPNAM(INDPRA),REP(ISYM-1),K-M,PRPNAM(INDPRB),REP(ISYMB-1)
            NEXC=KEXCNV(ISYM)
            DO IEXC = 1,NEXC
               FOSC(IEXC,ISYM) = FOSC(IEXC,ISYM) &
                    + FAC*WROT*ATMPPF(IEXC,INDAP,ISYM)*ATMPPF(IEXC,INDBP,ISYM)
!               WRITE(6,*)  'TEST..',FAC,WROT,ATMPPF(IEXC,INDAP,ISYM),ATMPPF(IEXC,INDBP,ISYM)  
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    END SUBROUTINE ROTFAC
      
