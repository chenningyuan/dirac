      SUBROUTINE GENOSC_BED(K,FOSC,ATMPPF,OMEGA)
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbxpp.h"
      real*8 :: FOSC(MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      integer, allocatable :: indices(:)
      integer :: j,k,k2,iterm,jterm
      iterm = 0
      jterm = 0
      K2=K+K
      WRITE(6,'(/A,I4)') 'GENOSC_BED: Rotational weights of order: ',K2
      CALL PRSYMB(6,'-',60,0)
      WRITE(6,'(18X,100A1)') (' ',J=1,K2+2),'P','Q'
!.....Accumulate oscillator strengths
      FOSC=0.0D0
      IF(K.EQ.0) THEN
        call rotfac_bed(indices,0,iterm,jterm,FOSC,ATMPPF,OMEGA)
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
    END SUBROUTINE GENOSC_BED

      recursive subroutine bedloop(m,k2,indices,iterm,jterm,FOSC,ATMPPF,OMEGA)
#include "dgroup.h"
#include "dcbxpp.h"
      real*8 :: FOSC(MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      integer :: indices(k2),ind(3),k2
      if(m==k2) then
        do i=1,3
           indices(m) = i
           call rotfac_bed(indices,k2,iterm,jterm,FOSC,ATMPPF,OMEGA)
        enddo
        return
      else
        do i=1,3
          indices(m) = i
          call bedloop(m+1,k2,indices,iterm,jterm,FOSC,ATMPPF,OMEGA)
        enddo
      endif 
      END SUBROUTINE BEDLOOP
      
      SUBROUTINE ROTFAC_BED(indices,k2,iterm,jterm,FOSC,ATMPPF,OMEGA)
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
            IF(ISYMB.NE.ISYM) CALL QUIT('ROTFAC_BED: Symmetry problem !')
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
    END SUBROUTINE ROTFAC_BED
      
      SUBROUTINE GENOSC_MPOL(K,FOSC,ATMPPF,OMEGA)
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbxpp.h"
      real*8 :: FOSC(4,MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      integer, allocatable :: indices(:)
      integer :: j,k,k2,iterm,jterm
      iterm = 0
      jterm = 0
      K2=K+K
      WRITE(6,'(/A,I4)') 'GENOSC_MPOL: Rotational weights of order: ',K2
      CALL PRSYMB(6,'-',60,0)
      WRITE(6,'(18X,100A1)') (' ',J=1,K2+2),'P','Q'
!.....Accumulate oscillator strengths
      FOSC=0.0D0
      IF(K.EQ.0) THEN
        call rotfac_mpol(indices,0,iterm,jterm,FOSC,ATMPPF,OMEGA)
      ELSE
        allocate(indices(k2))
        call mpolloop(1,k2,indices,iterm,jterm,FOSC,ATMPPF,OMEGA)
        deallocate(indices)
     ENDIF
!....Add final factors for oscillator strengths
     DO ISYM = 1,NBSYM
       NEXC = KEXCNV(ISYM)
       DO IEXC = 1,NEXC
         FAC = (OMEGA(IEXC,ISYM)/CVAL)**K2*2.0D0
         FOSC(2,IEXC,ISYM)=FOSC(2,IEXC,ISYM)*FAC*OMEGA(IEXC,ISYM)
         FOSC(3,IEXC,ISYM)=FOSC(3,IEXC,ISYM)*FAC
         FOSC(4,IEXC,ISYM)=FOSC(4,IEXC,ISYM)*FAC/OMEGA(IEXC,ISYM)
         FOSC(1,IEXC,ISYM)=FOSC(2,IEXC,ISYM)+FOSC(3,IEXC,ISYM)+FOSC(4,IEXC,ISYM)
       ENDDO
     ENDDO   
     write(6,*) 'Total/non-zero terms: ',ITERM,JTERM
     CALL PRSYMB(6,'-',60,0)
    END SUBROUTINE GENOSC_MPOL

      recursive subroutine mpolloop(m,k2,indices,iterm,jterm,FOSC,ATMPPF,OMEGA)
#include "dgroup.h"
#include "dcbxpp.h"
      real*8 :: FOSC(4,MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      integer :: indices(k2),ind(3),k2
      if(m==k2) then
        do i=1,3
           indices(m) = i
           call rotfac_mpol(indices,k2,iterm,jterm,FOSC,ATMPPF,OMEGA)
        enddo
        return
      else
        do i=1,3
          indices(m) = i
          call mpolloop(m+1,k2,indices,iterm,jterm,FOSC,ATMPPF,OMEGA)
        enddo
      endif 
      END SUBROUTINE MPOLLOOP
      
      SUBROUTINE ROTFAC_MPOL(indices,k2,iterm,jterm,FOSC,ATMPPF,OMEGA)
#include "pgroup.h"
#include "dgroup.h"
#include "dcbxpp.h"
#include "dcbxpr.h"      
      real*8 :: FOSC(4,MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
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
!...........Left EPOL operator: order (K+M+1)
            IORD    = K + M + 1
            IND     = 0
            IND(IP) = IND(IP) + 1
            DO L = 1,K+M
              J      = INDICES(L)
              IND(J) = IND(J) + 1
            ENDDO
            IA      = IORD - IND(1)
            IB      = IND(3)
            INDAPE  = IEPOFF(IORD) + IB + 1 + (IA+1)*IA/2
            INDPEL  = LPPAPU(INDAPE)
            ISYMEL  = IPRPSYM(INDPEL)
            FACEL   = D1/GAMMA(REAL(IORD+1))
!...........Right EPOL operator: order (K-M+1)
            IORD    = K - M + 1
            IND     = 0
            IND(IQ) = IND(IQ) + 1
            DO L = K+M+1,K2
              J      = INDICES(L)
              IND(J) = IND(J) + 1
            ENDDO
            IA = IORD - IND(1)
            IB = IND(3)
            INDBPE  = IEPOFF(IORD) + IB + 1 + (IA+1)*IA/2
            INDPER  = LPPAPU(INDBPE)
            ISYMER  = IPRPSYM(INDPER)
            FACER   = D1/GAMMA(REAL(IORD+1))
!...........Left MPOL operator: order (K+M)
            INDAPM  = -1
            IORD    = K + M
            IF(IORD.GT.0) THEN    ! No magnetic monopole....
              LM      = IORD - 1
              IND     = 0
              CALL LEVI_CEVITA(INDICES(K+M),IP,IR,IJKL)
              IF(IJKL.NE.0) THEN
                INDAPM  = IMPOFF(IORD) + IR
                IF(LM.GT.0) THEN
                  IND   = 0
                  DO L = 1,K+M-1
                    J      = INDICES(L)
                    IND(J) = IND(J) + 1
                  ENDDO
                  IA     = LM - IND(1)
                  IB     = IND(3)
                  INDAPM = INDAPM + 3*(IB+(IA+1)*IA/2)
                ENDIF
                INDPML  = LPPAPU(INDAPM)
                ISYMML  = IPRPSYM(INDPML)
                FACML = IJKL*D1/GAMMA(REAL(IORD+1))
              ENDIF
            ENDIF
!...........Right MPOL operator: order (K+M)
            INDBPM  = -1
            IORD    = K - M
            IF(IORD.GT.0) THEN    ! No magnetic monopole....
              LM      = IORD - 1
              IND     = 0
              CALL LEVI_CEVITA(INDICES(K2),IQ,IR,IJKR)
              IF(IJKR.NE.0) THEN
                INDBPM   = IMPOFF(IORD) + IR
                IF(LM.GT.0) THEN
                  IND   = 0
                  DO L = K+M+1,K2-1
                    J      = INDICES(L)
                    IND(J) = IND(J) + 1
                  ENDDO
                  IA     = LM - IND(1)
                  IB     = IND(3)
                  INDBPM = INDBPM + 3*(IB+(IA+1)*IA/2)
                ENDIF
                INDPMR  = LPPAPU(INDBPM)
                ISYMMR  = IPRPSYM(INDPMR)
                FACMR   = IJKR*D1/GAMMA(REAL(IORD+1))
              ENDIF
            ENDIF
!...........Form product
            IF(INDAPM.GT.0) THEN            
               WRITE(6,'(A,2(2X,I3,2X,A16,2X,A3),I3)') 'Left : ',&
                    K+M+1,PRPNAM(INDPEL),REP(ISYMEL-1),K+M,PRPNAM(INDPML),REP(ISYMML-1),IJKL
            ELSE
              WRITE(6,'(A,(2X,I3,2X,A16,2X,A3))')  'Left : ',K+M+1,PRPNAM(INDPEL),REP(ISYMEL-1)
            ENDIF
            IF(INDBPM.GT.0) THEN
               WRITE(6,'(A,2(2X,I3,2X,A16,2X,A3),I3)') 'Right: ',&
                    K-M+1,PRPNAM(INDPER),REP(ISYMER-1),K-M,PRPNAM(INDPMR),REP(ISYMMR-1),IJKR
            ELSE
              WRITE(6,'(A,(2X,I3,2X,A16,2X,A3))')  'Right: ',K-M+1,PRPNAM(INDPER),REP(ISYMER-1)
            ENDIF
            ISYM=ISYMEL           
            NEXC=KEXCNV(ISYM)
            IF(ISYMER.NE.ISYMEL.AND.ISYMML.NE.ISYMEL.AND.ISYMMR.NE.ISYMEL) CALL QUIT('ROTFAC_MPOL: Symmetry problem !')
!           Q-Q part
            FACQQ =  FAC*WROT*FACEL*FACER
            DO IEXC = 1,NEXC
               FOSC(2,IEXC,ISYM) = FOSC(2,IEXC,ISYM) &
                    + FACQQ*ATMPPF(IEXC,INDAPE,ISYM)*ATMPPF(IEXC,INDBPE,ISYM)
!               WRITE(6,*)  'TEST..',FAC,WROT,ATMPPF(IEXC,INDAP,ISYM),ATMPPF(IEXC,INDBP,ISYM)  
            ENDDO
!           Q-M part
            IF(INDBPM.GT.0) THEN            
              FACQM = -FAC*WROT*FACEL*FACMR
              DO IEXC = 1,NEXC
                 FOSC(3,IEXC,ISYM) = FOSC(3,IEXC,ISYM) &
                      + FACQM*ATMPPF(IEXC,INDAPE,ISYM)*ATMPPF(IEXC,INDBPM,ISYM)
!               WRITE(6,*)  'TEST..',FAC,WROT,ATMPPF(IEXC,INDAP,ISYM),ATMPPF(IEXC,INDBP,ISYM)  
              ENDDO
            ENDIF
!           M-Q part
            IF(INDAPM.GT.0) THEN            
              FACMQ = -FAC*WROT*FACML*FACER
              DO IEXC = 1,NEXC
                 FOSC(3,IEXC,ISYM) = FOSC(3,IEXC,ISYM) &
                      + FACQM*ATMPPF(IEXC,INDAPM,ISYM)*ATMPPF(IEXC,INDBPE,ISYM)
!               WRITE(6,*)  'TEST..',FAC,WROT,ATMPPF(IEXC,INDAP,ISYM),ATMPPF(IEXC,INDBP,ISYM)  
              ENDDO
            ENDIF
!           M-M part
            IF(INDAPM.GT.0.AND.INDBPM.GT.0) THEN            
              FACMM =  FAC*WROT*FACML*FACMR
              DO IEXC = 1,NEXC
                 FOSC(4,IEXC,ISYM) = FOSC(4,IEXC,ISYM) &
                      + FACMM*ATMPPF(IEXC,INDAPM,ISYM)*ATMPPF(IEXC,INDBPM,ISYM)
!               WRITE(6,*)  'TEST..',FAC,WROT,ATMPPF(IEXC,INDAP,ISYM),ATMPPF(IEXC,INDBP,ISYM)  
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      END SUBROUTINE ROTFAC_MPOL
      
