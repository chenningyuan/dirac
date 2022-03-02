      SUBROUTINE GENOSC_BED(K,FOSC,ATMPPF,OMEGA)
!***********************************************************************
!
!     Calculate rotationally averaged oscillator strengths
!
!***********************************************************************
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbxpp.h"
      real*8 :: FOSC(MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      integer, allocatable :: indices(:)
      integer :: j,k,k2,iterm,jterm
      iterm = 0
      jterm = 0
      K2=K+K
      IF(IPRXPP.GE.4) THEN
        WRITE(6,'(/A,I4)') 'GENOSC_BED: Rotational weights of order: ',K2
        CALL PRSYMB(6,'-',60,0)
        WRITE(6,'(18X,100A1)') (' ',J=1,K2+2),'P','Q'
      ENDIF
!.....Accumulate oscillator strengths
      FOSC=0.0D0
      IF(K.EQ.0) THEN
        allocate(indices(1))
        call rotfac_bed(indices,0,iterm,jterm,FOSC,ATMPPF,OMEGA)
        deallocate(indices)
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
      IF(IPRXPP.GE.4) THEN
        WRITE(6,*) 'Total/non-zero terms: ',ITERM,JTERM
        CALL PRSYMB(6,'-',60,0)
      ENDIF
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
#include "implicit.h"        
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
          IF(IPRXPP.GE.4) THEN
            IFAC= NINT(FAC)
            WRITE(DENOM,'(A,I8,A)') '/',IFAC,':'
            WRITE(6,'(F12.1,A10,2X,18A1)') WROT*IFAC,DENOM,&
                 (CHAR(87+INDICES(J)),J=1,K2),' ',' ',CHAR(87+IP),CHAR(87+IQ)
          ENDIF
          DO M = 0,K
            FAC = D2*(-D1)**M
            IF(M==0) FAC = D1
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
            IF(IPRXPP.GE.5) THEN
              WRITE(6,'(2(2X,I3,2X,A16,2X,A3))') K+M,PRPNAM(INDPRA),REP(ISYM-1),K-M,PRPNAM(INDPRB),REP(ISYMB-1)
            ENDIF 
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


      SUBROUTINE GENOSC_VPOL(K,FOSC,ATMPPF,OMEGA)
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbxpp.h"
      real*8 :: FOSC(4,MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      integer, allocatable :: indices(:)
      integer :: j,k,k2,iterm,jterm
      iterm = 0
      jterm = 0
      K2=K+K
      IF(IPRXPP.GE.4) THEN
        WRITE(6,'(/A,I4)') 'GENOSC_VPOL: Rotational weights of order: ',K2
        CALL PRSYMB(6,'-',60,0)
        WRITE(6,'(18X,100A1)') (' ',J=1,K2+2),'P','Q'
      ENDIF
!.....Accumulate oscillator strengths
      FOSC=0.0D0
      IF(K.EQ.0) THEN
        allocate(indices(1))
        call rotfac_vpol(indices,0,iterm,jterm,FOSC,ATMPPF,OMEGA)
        deallocate(indices)
      ELSE
        allocate(indices(k2))
        call vpolloop(1,k2,indices,iterm,jterm,FOSC,ATMPPF,OMEGA)
        deallocate(indices)
     ENDIF
!....Add final factors for oscillator strengths
     DO ISYM = 1,NBSYM
       NEXC = KEXCNV(ISYM)
       DO IEXC = 1,NEXC
         FAC = (OMEGA(IEXC,ISYM)/CVAL)**K2*2.0D0
         FOSC(2,IEXC,ISYM)=FOSC(2,IEXC,ISYM)*FAC/OMEGA(IEXC,ISYM)
         FOSC(3,IEXC,ISYM)=FOSC(3,IEXC,ISYM)*FAC/OMEGA(IEXC,ISYM)
         FOSC(4,IEXC,ISYM)=FOSC(4,IEXC,ISYM)*FAC/OMEGA(IEXC,ISYM)
         FOSC(1,IEXC,ISYM)=FOSC(2,IEXC,ISYM)+FOSC(3,IEXC,ISYM)+FOSC(4,IEXC,ISYM)
       ENDDO
     ENDDO
     IF(IPRXPP.GE.4) THEN
       WRITE(6,*) 'Total/non-zero terms: ',ITERM,JTERM
       CALL PRSYMB(6,'-',60,0)
     ENDIF
     END SUBROUTINE GENOSC_VPOL

      recursive subroutine vpolloop(m,k2,indices,iterm,jterm,FOSC,ATMPPF,OMEGA)
#include "dgroup.h"
#include "dcbxpp.h"
      real*8 :: FOSC(4,MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      integer :: indices(k2),ind(3),k2
      if(m==k2) then
        do i=1,3
           indices(m) = i
           call rotfac_vpol(indices,k2,iterm,jterm,FOSC,ATMPPF,OMEGA)
        enddo
        return
      else
        do i=1,3
          indices(m) = i
          call vpolloop(m+1,k2,indices,iterm,jterm,FOSC,ATMPPF,OMEGA)
        enddo
      endif 
      END SUBROUTINE VPOLLOOP
      
      SUBROUTINE ROTFAC_VPOL(indices,k2,iterm,jterm,FOSC,ATMPPF,OMEGA)
#include "implicit.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbxpp.h"
#include "dcbxpr.h"      
      real*8 :: FOSC(4,MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      integer :: indices(k2),ind(3),IOFF(3,2)
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
          IF(IPRXPP.GE.4) THEN
            IFAC= NINT(FAC)
            WRITE(DENOM,'(A,I8,A)') '/',IFAC,':'
            WRITE(6,'(F12.1,A10,2X,18A1)') WROT*IFAC,DENOM,&
                 (CHAR(87+INDICES(J)),J=1,K2),' ',' ',CHAR(87+IP),CHAR(87+IQ)
          ENDIF
          DO M = 0,K
            FAC = D2*(-D1)**M
            IF(M.EQ.0) FAC = D1
!...........Left EPOL operator: order (K+M+1)
            IORD    = K + M + 1
            IF (IORD.EQ.1) THEN
              INDAPE = IBDOFF(0)+IP
            ELSE
              IND     = 0
              DO L = 1,K+M-1
                J      = INDICES(L)
                IND(J) = IND(J) + 1
              ENDDO
              LE = IORD - 2
              IA      = LE - IND(1)
              IB      = IND(3)
              IF (IP.EQ.INDICES(K+M)) THEN 
                INDAPE = IVPOFF(IORD) + IP + 9*(IB+(IA+1)*IA/2)
              ELSE  
                CALL LEVI_CEVITA(IP,INDICES(K+M),IR,IJKL)             
                IF (IORD.EQ.2) THEN
                  IOFF = RESHAPE((/3,1,2,2,3,1/),SHAPE(IOFF)) 
                ELSE 
                  IOFF = RESHAPE((/5,1,3,4,6,2/),SHAPE(IOFF))
                ENDIF 
                IF (IJKL.EQ.1) THEN
!                  IOFF = [5,1,3]
                  INDAPE = IVPOFF(IORD)+IOFF(IP,1) + 3 + 9*(IB+(IA+1)*IA/2)
                ELSEIF(IJKL.EQ.-1) THEN 
!                  IOFF = [4,6,2]
                  INDAPE = IVPOFF(IORD)+IOFF(IP,2) + 3 + 9*(IB+(IA+1)*IA/2)
                ENDIF
              ENDIF 
            ENDIF 
              INDPEL  = LPPAPU(INDAPE)
              ISYMEL  = IPRPSYM(INDPEL)
              FACEL   = D1/GAMMA(REAL(IORD+1))
!...........Right EPOL operator: order (K-M+1)
            IORD    = K - M + 1
            IF (IORD.EQ.1) THEN
              INDBPE = IBDOFF(0)+IQ
            ELSE
              IND     = 0
              DO L = K+M+1,K2-1
                J      = INDICES(L)
                IND(J) = IND(J) + 1
              ENDDO
              LE = IORD - 2
              IA = LE - IND(1)
              IB = IND(3)
              IF (IQ.EQ.INDICES(K2)) THEN 
                INDBPE = IVPOFF(IORD) + IQ + 9*(IB+(IA+1)*IA/2)
              ELSE  
                  CALL LEVI_CEVITA(IQ,INDICES(K2),IR,IJKL)             
                  IF (IORD.EQ.2) THEN
                    IOFF = RESHAPE((/3,1,2,2,3,1/),SHAPE(IOFF)) 
                  ELSE 
                    IOFF = RESHAPE((/5,1,3,4,6,2/),SHAPE(IOFF))
                  ENDIF 
                  IF (IJKL.EQ.1) THEN
!                    IOFF = [5,1,3]
                    INDBPE = IVPOFF(IORD)+IOFF(IQ,1) + 3 + 9*(IB+(IA+1)*IA/2)
                  ELSEIF(IJKL.EQ.-1) THEN 
!                    IOFF = [4,6,2]
                    INDBPE = IVPOFF(IORD)+IOFF(IQ,2) + 3 + 9*(IB+(IA+1)*IA/2)
                  ENDIF
                ENDIF 
            ENDIF 
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
! original code
!            IF(IPRXPP.GE.5) THEN            
!              IF(INDAPM.GT.0) THEN            
!                 WRITE(6,'(A,2(2X,I3,2X,A16,2X,A3),I3)') 'Left : ',&
!                      K+M+1,PRPNAM(INDPEL),REP(ISYMEL-1),K+M,PRPNAM(INDPML),REP(ISYMML-1),IJKL
!              ELSE
!                WRITE(6,'(A,(2X,I3,2X,A16,2X,A3))')  'Left : ',K+M+1,PRPNAM(INDPEL),REP(ISYMEL-1)
!              ENDIF
!              IF(INDBPM.GT.0) THEN
!                 WRITE(6,'(A,2(2X,I3,2X,A16,2X,A3),I3)') 'Right: ',&
!                      K-M+1,PRPNAM(INDPER),REP(ISYMER-1),K-M,PRPNAM(INDPMR),REP(ISYMMR-1),IJKR
!              ELSE
!                WRITE(6,'(A,(2X,I3,2X,A16,2X,A3))')  'Right: ',K-M+1,PRPNAM(INDPER),REP(ISYMER-1)
!              ENDIF
!            ENDIF

            IF(IPRXPP.GE.5) THEN   
              WRITE(6,*)'INDICES: ', INDICES        
              WRITE(6,'(A3,I2,1X,A3,I2)')'IP=',IP,'IQ=',IQ
              IF(INDAPM.GT.0) THEN            
                 WRITE(6,'(A,2(2X,I3,2X,A16,2X,A3),I3)') 'Left : ',&
                      K+M+1,PRPNAM(INDPEL),REP(ISYMEL-1),K+M,PRPNAM(INDPML),REP(ISYMML-1),IJKL
              ELSE
                WRITE(6,'(A,(2X,I3,2X,A16,2X,A3))')  'Left : ',K+M+1,PRPNAM(INDPEL)&
                ,REP(ISYMEL-1)
              ENDIF
              IF(INDBPM.GT.0) THEN
                 WRITE(6,'(A,2(2X,I3,2X,A16,2X,A3),I3)') 'Right: ',&
                      K-M+1,PRPNAM(INDPER),REP(ISYMER-1),K-M,PRPNAM(INDPMR),REP(ISYMMR-1),IJKR            
              ELSE
                WRITE(6,'(A,(2X,I3,2X,A16,2X,A3))')  'Right: ',K-M+1,PRPNAM(INDPER)&
                ,REP(ISYMER-1)
              ENDIF
            ENDIF


            ISYM=ISYMEL           
            NEXC=KEXCNV(ISYM)
            IF(ISYMER.NE.ISYMEL.AND.ISYMML.NE.ISYMEL.AND.ISYMMR.NE.ISYMEL) CALL QUIT('ROTFAC_VPOL: Symmetry problem !')
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
                      + FACMQ*ATMPPF(IEXC,INDAPM,ISYM)*ATMPPF(IEXC,INDBPE,ISYM)
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
      END SUBROUTINE ROTFAC_VPOL
     
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
      IF(IPRXPP.GE.4) THEN
        WRITE(6,'(/A,I4)') 'GENOSC_MPOL: Rotational weights of order: ',K2
        CALL PRSYMB(6,'-',60,0)
        WRITE(6,'(18X,100A1)') (' ',J=1,K2+2),'P','Q'
      ENDIF
!.....Accumulate oscillator strengths
      FOSC=0.0D0
      IF(K.EQ.0) THEN
        allocate(indices(1))
        call rotfac_mpol(indices,0,iterm,jterm,FOSC,ATMPPF,OMEGA)
        deallocate(indices)
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
     IF(IPRXPP.GE.4) THEN
       WRITE(6,*) 'Total/non-zero terms: ',ITERM,JTERM
       CALL PRSYMB(6,'-',60,0)
     ENDIF
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
#include "implicit.h"
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
          IF(IPRXPP.GE.4) THEN
            IFAC= NINT(FAC)
            WRITE(DENOM,'(A,I8,A)') '/',IFAC,':'
            WRITE(6,'(F12.1,A10,2X,18A1)') WROT*IFAC,DENOM,&
                 (CHAR(87+INDICES(J)),J=1,K2),' ',' ',CHAR(87+IP),CHAR(87+IQ)
          ENDIF
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
            IF(IPRXPP.GE.5) THEN            
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
                      + FACMQ*ATMPPF(IEXC,INDAPM,ISYM)*ATMPPF(IEXC,INDBPE,ISYM)
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
      
      SUBROUTINE OSC_BED(K,FOSC,ATMPPF,OMEGA)
!***********************************************************************
!
!     Calculated oscillator strength for oriented sample
!
!***********************************************************************
#include "implicit.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbxpp.h"
      real*8              :: FOSC(MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      real*8,allocatable  :: TMOML(:,:),TMOMR(:,:)
      PARAMETER(D2=2.0D0,D1=1.0D0,D0=0.0D0)
      allocate(tmoml(MAXEXC,NBSYM))
      allocate(tmomr(MAXEXC,NBSYM))      
      K2=K+K
      DO M = 0,K
        FAC = D2*(-D1)**M
        IF(M==0) FAC = D1
!.......Accumulate left operator .. order (K+M)
        IORD  = K+M
        TMOML = D0
        CALL accBEDMOM(TMOML,ATMPPF,IORD)
!.......Accumulate right operator .. order (K-M)
        IORD  = K-M
        TMOMR = D0
        CALL accBEDMOM(TMOMR,ATMPPF,IORD)
!.......Accumulate oscillator strengths
        DO ISYM = 1,NBSYM        
          NEXC=KEXCNV(ISYM)
          DO IEXC=1,NEXC
 ! WRITE(6,*) "IEXC ", IEXC, "ISYM ", ISYM 
 ! WRITE(6,*)  "BED L", TMOML(IEXC,ISYM), "R ", TMOMR(IEXC,ISYM)


             FOSC(IEXC,ISYM) = FOSC(IEXC,ISYM) + FAC*TMOML(IEXC,ISYM)*TMOMR(IEXC,ISYM)
!             WRITE(6,*) 'TEST..',K,M,IEXC,ISYM,FOSC(IEXC,ISYM),TMOML(IEXC,ISYM),TMOMR(IEXC,ISYM)
          ENDDO
        ENDDO
      ENDDO
      deallocate(tmoml)
      deallocate(tmomr)      
!....Add final factors for oscillator strengths
      DO ISYM = 1,NBSYM
        NEXC = KEXCNV(ISYM)
        DO IEXC = 1,NEXC
          FOSC(IEXC,ISYM)=FOSC(IEXC,ISYM)*(OMEGA(IEXC,ISYM)/CVAL)**K2*D2/OMEGA(IEXC,ISYM)
        ENDDO
      ENDDO   
      END SUBROUTINE OSC_BED

      SUBROUTINE accBEDMOM(TMOM,ATMPPF,IORDER)
#include "implicit.h"        
#include "dgroup.h"
#include "dcbxpr.h"      
#include "dcbxpp.h"
      real*8  :: TMOM(MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM)
      INDAP = IBDOFF(IORDER)
      IBKA = 1     ! Binomial coefficient (IORDER A)
      DO IA = 0,IORDER
        IX = IORDER - IA
        IBAB = 1     ! Binomial coefficient (A B)
        DO IB = 0,IA
          IY = IA - IB
          IZ = IB
          TFAC=(UWAVE(1)**IX)*(UWAVE(2)**IY)*(UWAVE(3)**IZ)*IBKA*IBAB
          DO IP = 1,3
            INDAP   = INDAP + 1
            INDPRA  = LPPAPU(INDAP)
            ISYM    = IPRPSYM(INDPRA)
            NEXC    = KEXCNV(ISYM)
            DO IEXC=1,NEXC
              TMOM(IEXC,ISYM) = TMOM(IEXC,ISYM)+UPOL(IP,1)*TFAC*ATMPPF(IEXC,INDAP,ISYM)
!      write(6,*) "INDAP;IX;IY;IZ;IP;TMOM",INDAP," ",IX," ",IY," ",IZ," ",IP," ",TMOM(IEXC,ISYM)
            ENDDO
          ENDDO
          IBAB = IBAB*IY/(IB+1)
       ENDDO
       IBKA = IBKA*IX/(IA+1)
      ENDDO
      RETURN
      END SUBROUTINE accBEDMOM

      SUBROUTINE accVMOM(TMOM,ATMPPF,IORDER)
#include "implicit.h"        
#include "dgroup.h"
#include "dcbxpr.h"      
#include "dcbxpp.h"
      INTEGER CARPOW(3,3), INDU(3)
      real*8  :: TMOM(MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM)


      CARPOW = 0.0D0    
      DO I = 1,3
        CARPOW(I,I) = 1
      ENDDO
      FAC = 1.0D0
      DO I = 2,IORDER
        FAC = FAC/I
      ENDDO

      INDAP = IVPOFF(IORDER)
      IBKA  = 1
      DO IA = 0,IORDER - 2
        IX    = IORDER - 2 - IA
        IBAB  = 1
        DO IB = 0,IA
          IY    = IA - IB 
          IZ    = IB 

          DO IP = 1, 3
            INDAP  = INDAP + 1 
            INDPRA = LPPAPU(INDAP)
            ISYM   = IPRPSYM(INDPRA)
            NEXC   = KEXCNV(ISYM)
            JX = IX + CARPOW(1,IP)
            JY = IY + CARPOW(2,IP)
            JZ = IZ + CARPOW(3,IP)
            TFAC  = (UWAVE(1)**JX)*(UWAVE(2)**JY)*(UWAVE(3)**JZ)*IBKA*IBAB*FAC
            DO IEXC = 1,NEXC
              TMOM(IEXC,ISYM) = TMOM(IEXC,ISYM)+UPOL(IP,1)*TFAC*ATMPPF(IEXC,INDAP,ISYM)
!write(6,*) "INDAP;NX;NY;NZ;INDU(IP);ATMPPF",INDAP," ",NX," ",NY," ",NZ," ",INDU(IP)," ",ATMPPF(IEXC,INDAP,ISYM)
!write(6,*) "INDAP;NX;NY;NZ;TFAC;UPOL;ATMPPF",INDAP," ",NX," ",NY," ",NZ," ",TFAC," ",UPOL(IP,1)," ",ATMPPF(IEXC,INDAP,ISYM) 
            ENDDO
          ENDDO

          DO IP = 1, 3
            INDAP  = INDAP + 1 
            INDPRA = LPPAPU(INDAP)
            ISYM   = IPRPSYM(INDPRA)
            NEXC   = KEXCNV(ISYM)
            INDU = [2,3,1]
            KX = CARPOW(2,IP)
            KY = CARPOW(3,IP)
            KZ = CARPOW(1,IP)   
            JX = IX + CARPOW(2,IP)
            JY = IY + CARPOW(3,IP)
            JZ = IZ + CARPOW(1,IP)
            TFAC  = (UWAVE(1)**JX)*(UWAVE(2)**JY)*(UWAVE(3)**JZ)*IBKA*IBAB*FAC
            DO IEXC = 1,NEXC
              TMOM(IEXC,ISYM) = TMOM(IEXC,ISYM)+UPOL(INDU(IP),1)*TFAC*ATMPPF(IEXC,INDAP,ISYM)
!write(6,*) "INDAP;NX;NY;NZ;INDU(IP);ATMPPF",INDAP," ",NX," ",NY," ",NZ," ",INDU(IP)," ",ATMPPF(IEXC,INDAP,ISYM)   
!write(6,*) "INDAP;NX;NY;NZ;TFAC;UPOL;ATMPPF",INDAP," ",NX," ",NY," ",NZ," ",TFAC," ",UPOL(INDU(IP),1)," ",ATMPPF(IEXC,INDAP,ISYM) 
            ENDDO
            IF (IORDER.GT.2) THEN
              INDAP = INDAP + 1
            ENDIF
              TFAC  = (UWAVE(1)**IX)*(UWAVE(2)**IY)*(UWAVE(3)**IZ)*IBKA*IBAB*FAC
              TFAC = TFAC*(UPOL(1,1)**KX)*(UPOL(2,1)**KY)*(UPOL(3,1)**KZ)
            DO IEXC = 1,NEXC
              TMOM(IEXC,ISYM) = TMOM(IEXC,ISYM)+UWAVE(INDU(IP))*TFAC*ATMPPF(IEXC,INDAP,ISYM)
!write(6,*) "INDAP;NX;NY;NZ;TFAC;UWAVE;ATMPPF",INDAP," ",NX," ",NY," ",NZ," ",TFAC," ",UWAVE(INDU(IP))," ",ATMPPF(IEXC,INDAP,ISYM) 

            ENDDO
          ENDDO
          IBAB = IBAB*IY/(IB+1)
        ENDDO
        IBKA = IBKA*IX/(IA+1)
      ENDDO


!      DO ISYM = 1,2
!        NEXC   = KEXCNV(ISYM)
!        DO IEXC = 1,NEXC
!
!        write(6,*) "excitation ", iexc, " order ", IORDER,  " sym ", ISYM, "  ", TMOM(IEXC,ISYM) 
!      ENDDO 
!      ENDDO


      RETURN

      END SUBROUTINE accVMOM


      SUBROUTINE OSC_VPOL(K,FOSC,ATMPPF,OMEGA)
#include "implicit.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbxpp.h"
      real*8              :: FOSC(4,MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      real*8,allocatable  :: TEMOML(:,:),TMMOML(:,:),TEMOMR(:,:),TMMOMR(:,:),TMOMTESTL(:,:),TMOMTESTR(:,:)
      PARAMETER(D2=2.0D0,D1=1.0D0,D0=0.0D0)
      allocate(temoml(MAXEXC,NBSYM))
      allocate(tmmoml(MAXEXC,NBSYM))      
      allocate(temomr(MAXEXC,NBSYM))
      allocate(tmmomr(MAXEXC,NBSYM))
      !allocate(TMOMTESTR(MAXEXC,NBSYM))
      !allocate(TMOMTESTL(MAXEXC,NBSYM))             
      K2=K+K
      DO M = 0,K
        FAC = D2*(-D1)**M
        IF(M==0) FAC = D1
!.......Accumulate left operator .. order (K+M)
        IORD  = K+M
!.........Electric multipole        
          TEMOML = D0
          IF (IORD.EQ.0) THEN
            CALL accBEDMOM(TEMOML,ATMPPF,0)
          ELSE
            CALL accVMOM(TEMOML,ATMPPF,IORD+1)
          ENDIF
!.........Magnetic multipole        
          TMMOML = D0
          CALL accMMOM(TMMOML,ATMPPF,IORD)
!         DEBUG
         ! CALL accBEDMOM(TMOMTESTL,ATMPPF,IORD)


!.......Accumulate right operator .. order (K-M)
        IORD  = K-M
!.........Electric multipole        
          TEMOMR = D0
          IF (IORD.EQ.0) THEN
            CALL accBEDMOM(TEMOMR,ATMPPF,0)
          ELSE
            CALL accVMOM(TEMOMR,ATMPPF,IORD+1)
          ENDIF
!.........Magnetic multipole        
          TMMOMR = D0
          CALL accMMOM(TMMOMR,ATMPPF,IORD)
!         DEBUG
         ! CALL accBEDMOM(TMOMTESTR,ATMPPF,IORD)

!.......Accumulate oscillator strengths
        DO ISYM = 1,NBSYM        
          NEXC=KEXCNV(ISYM)
          DO IEXC=1,NEXC
!  WRITE(6,*) "IEXC ", IEXC, "ISYM ", ISYM 
!  WRITE(6,*)  "VPOLL ", TEMOML(IEXC,ISYM), "MAGNL ", TMMOML(IEXC,ISYM)!, "BED ", TMOMTESTL(IEXC,ISYM)
!  WRITE(6,*)  "VPOLR ", TEMOMR(IEXC,ISYM), "MAGNR ", TMMOMR(IEXC,ISYM)
!  write(6,*)  "SUM ",TEMOML(IEXC,ISYM)-TMMOML(IEXC,ISYM)  
!  !WRITE(6,*) "DEBUG ", "VPOL ", TEMOMR(IEXC,ISYM), "MAGN ", TMMOMR(IEXC,ISYM), "BED ", TMOMTESTR(IEXC,ISYM)
!  !write(6,*) "SUM ",TEMOMR(IEXC,ISYM)+TMMOMR(IEXC,ISYM)  


             FOSC(2,IEXC,ISYM) = FOSC(2,IEXC,ISYM) + FAC*TEMOML(IEXC,ISYM)*TEMOMR(IEXC,ISYM)
             FOSC(3,IEXC,ISYM) = FOSC(3,IEXC,ISYM) - FAC*(TEMOML(IEXC,ISYM)*TMMOMR(IEXC,ISYM)&
                                                   +      TMMOML(IEXC,ISYM)*TEMOMR(IEXC,ISYM))
             FOSC(4,IEXC,ISYM) = FOSC(4,IEXC,ISYM) + FAC*TMMOML(IEXC,ISYM)*TMMOMR(IEXC,ISYM)
           ENDDO
        ENDDO
      ENDDO
      deallocate(temoml)
      deallocate(tmmoml)
      deallocate(temomr)      
      deallocate(tmmomr)
      !deallocate(TMOMTESTL)
      !deallocate(TMOMTESTR)
      !....Add final factors for oscillator strengths
      DO ISYM = 1,NBSYM
        NEXC = KEXCNV(ISYM)
        DO IEXC = 1,NEXC
          FAC = (OMEGA(IEXC,ISYM)/CVAL)**K2*2.0D0
          FOSC(2,IEXC,ISYM)=FOSC(2,IEXC,ISYM)*FAC/OMEGA(IEXC,ISYM)
          FOSC(3,IEXC,ISYM)=FOSC(3,IEXC,ISYM)*FAC/OMEGA(IEXC,ISYM)
          FOSC(4,IEXC,ISYM)=FOSC(4,IEXC,ISYM)*FAC/OMEGA(IEXC,ISYM)
          FOSC(1,IEXC,ISYM)=FOSC(2,IEXC,ISYM)+FOSC(3,IEXC,ISYM)+FOSC(4,IEXC,ISYM)
        ENDDO
      ENDDO   
      END SUBROUTINE OSC_VPOL





      SUBROUTINE OSC_MPOL(K,FOSC,ATMPPF,OMEGA)
#include "implicit.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbxpp.h"
      real*8              :: FOSC(4,MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM),OMEGA(MAXEXC,NBSYM)
      real*8,allocatable  :: TEMOML(:,:),TMMOML(:,:),TEMOMR(:,:),TMMOMR(:,:)
      PARAMETER(D2=2.0D0,D1=1.0D0,D0=0.0D0)
      allocate(temoml(MAXEXC,NBSYM))
      allocate(tmmoml(MAXEXC,NBSYM))      
      allocate(temomr(MAXEXC,NBSYM))
      allocate(tmmomr(MAXEXC,NBSYM))            
      K2=K+K
      DO M = 0,K
        FAC = D2*(-D1)**M
        IF(M==0) FAC = D1
!.......Accumulate left operator .. order (K+M)
        IORD  = K+M
!.........Electric multipole        
          TEMOML = D0
          CALL accEMOM(TEMOML,ATMPPF,IORD+1)
!.........Magnetic multipole        
          TMMOML = D0
          CALL accMMOM(TMMOML,ATMPPF,IORD)
!.......Accumulate right operator .. order (K-M)
        IORD  = K-M
!.........Electric multipole        
          TEMOMR = D0
          CALL accEMOM(TEMOMR,ATMPPF,IORD+1)
!.........Magnetic multipole        
          TMMOMR = D0
          CALL accMMOM(TMMOMR,ATMPPF,IORD)
!.......Accumulate oscillator strengths
        DO ISYM = 1,NBSYM        
          NEXC=KEXCNV(ISYM)
          DO IEXC=1,NEXC
             FOSC(2,IEXC,ISYM) = FOSC(2,IEXC,ISYM) + FAC*TEMOML(IEXC,ISYM)*TEMOMR(IEXC,ISYM)
             FOSC(3,IEXC,ISYM) = FOSC(3,IEXC,ISYM) - FAC*(TEMOML(IEXC,ISYM)*TMMOMR(IEXC,ISYM)&
                                                   +      TMMOML(IEXC,ISYM)*TEMOMR(IEXC,ISYM))
             FOSC(4,IEXC,ISYM) = FOSC(4,IEXC,ISYM) + FAC*TMMOML(IEXC,ISYM)*TMMOMR(IEXC,ISYM)
           ENDDO
        ENDDO
      ENDDO
      deallocate(temoml)
      deallocate(tmmoml)
      deallocate(temomr)      
      deallocate(tmmomr)
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
      END SUBROUTINE OSC_MPOL

      SUBROUTINE accMMOM(TMOM,ATMPPF,IORDER)
!     Accumulate magnetic multipole contributions to order IORDER        
#include "implicit.h"        
#include "dgroup.h"
#include "dcbxpr.h"      
#include "dcbxpp.h"
      real*8  :: TMOM(MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM)
      IF(IORDER.EQ.0) RETURN  
      FAC = 1.0D0
      DO I = 2,IORDER
        FAC = FAC/I
      ENDDO
      IORDM = IORDER - 1
      DO IP = 1,3
        INDAP = IMPOFF(IORDER)
        IBKA = 1     ! Binomial coefficient (IORDM A)
        DO IA = 0,IORDM
          IX = IORDM - IA
          IBAB = 1     ! Binomial coefficient (A B)
          DO IB = 0,IA
            IY = IA - IB
            IZ = IB
            TFAC=FAC*UPOL(IP,1)*(UWAVE(1)**IX)*(UWAVE(2)**IY)*(UWAVE(3)**IZ)*IBKA*IBAB
            DO IR = 1,3
              INDAP   = INDAP + 1
              CALL LEVI_CEVITA(IP,IR,IS,IPRS)
              IF(IPRS.NE.0) THEN
                INDPRA  = LPPAPU(INDAP)
                ISYM    = IPRPSYM(INDPRA)
                NEXC    = KEXCNV(ISYM)
                DO IEXC=1,NEXC
                  TMOM(IEXC,ISYM) = TMOM(IEXC,ISYM)+IPRS*UWAVE(IS)*TFAC*ATMPPF(IEXC,INDAP,ISYM)
!               WRITE(6,*) 'TEST2..',IORDER,IX,IY,IZ,IEXC,ISYM,TMOM(IEXC,ISYM),UPOL(IP,1),TFAC,ATMPPF(IEXC,INDAP,ISYM)
                ENDDO
              ENDIF 
            ENDDO
            IBAB = IBAB*IY/(IB+1)
          ENDDO
          IBKA = IBKA*IX/(IA+1)
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE accMMOM

      SUBROUTINE accEMOM(TMOM,ATMPPF,IORDER)
!     Accumulate electric multipole contributions to order IORDER        
#include "implicit.h"        
#include "dgroup.h"
#include "dcbxpr.h"      
#include "dcbxpp.h"
      real*8  :: TMOM(MAXEXC,NBSYM), ATMPPF(MAXEXC,NPPAPT,NBSYM)
      integer :: IND(3)  
      FAC = 1.0D0
      DO I = 2,IORDER
        FAC = FAC/I
      ENDDO
      IORDM = IORDER - 1
      IBKA = 1     ! Binomial coefficient (IORDM A)
      DO IA = 0,IORDM
        IX = IORDM - IA
        IBAB = 1     ! Binomial coefficient (A B)
        DO IB = 0,IA
          IY = IA - IB
          IZ = IB
          TFAC=FAC*(UWAVE(1)**IX)*(UWAVE(2)**IY)*(UWAVE(3)**IZ)*IBKA*IBAB
          DO IP = 1,3
            IND(1)  = IX
            IND(2)  = IY 
            IND(3)  = IZ
            IND(IP) = IND(IP) + 1
            IA1     = IORDER - IND(1)
            IB1     = IND(3)
            INDAP   = IEPOFF(IORDER) + IB1 + 1 + (IA1+1)*IA1/2
            INDPRA  = LPPAPU(INDAP)
            ISYM    = IPRPSYM(INDPRA)
            NEXC    = KEXCNV(ISYM)
            DO IEXC=1,NEXC
              TMOM(IEXC,ISYM) = TMOM(IEXC,ISYM)+UPOL(IP,1)*TFAC*ATMPPF(IEXC,INDAP,ISYM)


!              WRITE(6,*) 'TEST2..',IORDER,IX,IY,IZ,IEXC,ISYM,TMOM(IEXC,ISYM),UPOL(IP,1),TFAC,ATMPPF(IEXC,INDAP,ISYM)
            ENDDO
          ENDDO
          IBAB = IBAB*IY/(IB+1)
        ENDDO
        IBKA = IBKA*IX/(IA+1)
      ENDDO
      RETURN
      END SUBROUTINE accEMOM
      
      
