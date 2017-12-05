C  /* Deck getbedf */                                                  
      SUBROUTINE GET_BED(OMEGA,ISSYM,CMO,IBEIG,EVECR,IBTYP,
     &                   IBCVC,IBEVC,IBPVC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Build representation matrix of full light-matter interaction
C     Scalar integrals are sorted on symmetries (max. 8)
C
C     Upon entering this routine solution vectors have been written
C     to unformatted file PAMXVC, to be opened with unit number LUXVC
C     References are:
C     1. Nanna Holmgaard List, Joanna Kauczor, Trond Saue,
C        Hans Jørgen Aagaard Jensen and Patrick Norman, J. Chem. Phys. 142 (2015) 244111
C                               
C                                                                      
C***********************************************************************
      use orbital_rotation_indices
C     
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0,DP5=0.5D0)
#include "dcbgen.h"
#include "dcbbas.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbxrs.h"
#include "dcbxpp.h"
#include "orgcom.h"
#include "dcbxpr.h"
#include "dcbprl.h"            
#include "dcbibn.h"
#include "dcbwav.h"
C#include "cbiher.h"
      LOGICAL DOINT(2,2)
      CHARACTER*6 MXFORM, FMT
      INTEGER LWORK, NCOMP, NOPTP,IBEIG(*)
      DIMENSION OMEGA(MAXEXC),ISSYM(MAXEXC),CMO(*),WORK(*)
      real(8), allocatable              :: PRPINT(:,:),BBUF(:)
      real(8), allocatable              :: TMOM(:,:,:),OSC(:,:)
      character (len=8),  allocatable   :: LABINT(:)
      character (len=5)                 :: TPAR(3)
      integer, allocatable              :: INTREP(:),IPSYM2INT(:,:)
      integer, allocatable              :: IVECS(:)
      real(8)             :: PBUF(3,2),EVECR(NREDM,*)
      integer             :: IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*)
      real(8), allocatable              :: GPOE(:,:),GPOP(:,:)
      real(8), allocatable              :: XPOE(:,:),XPOP(:,:)
      real(8), allocatable              :: TMOM_P(:,:),TMOM_DER(:,:,:,:)
      real(8), allocatable              :: TMOM_DE(:,:),TMOM_DOM(:,:)
      real(8), allocatable              :: OSC_P(:),OSC_DER(:,:,:)
      real(8), allocatable              :: OSC_DE(:),OSC_DOM(:)
      logical, allocatable              :: first(:)
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE
C
#include "ibtfun.h"
      NOP = 2
C     
C     Allocate placeholders for integral labels and irreps 
C
      allocate(LABINT(8))
      allocate(INTREP(8))
      allocate(first(nz))
C
C     Integral directives
C
      IPRINT=0
      CALL PR1DIR('CXIKR   ',INTTYP,NPRPINT,INTREP,
     &            ANTI,SQUARE,DUMMY,LABINT,TRIANG,TRASPH,
     &            SOLVNT,0,.FALSE.,0,NBAST,NELMNT,IPRINT)
C     NPRPINT is by construction an even number,
C     the first and second halfs corresponding to the real (eq. 44) and imaginary parts (eq. 45)
      IF(NPRPLBL+NPRPINT.GT.MAXPRPLBL) THEN
        CALL QUIT('GET_BED:Redimension MAXPRPLBL !')
      ENDIF
      DO ICOMP = 1,NPRPINT
        PRPLBL (NPRPLBL+ICOMP) = LABINT(ICOMP)
        IPRLREP(NPRPLBL+ICOMP) = INTREP(ICOMP)
        PDOINT (NPRPLBL+ICOMP) = '0++0'
        IPRLTYP(NPRPLBL+ICOMP) =  1        
      ENDDO
      DOINT(1,1)=.FALSE.
      DOINT(2,1)=.TRUE.
      DOINT(1,2)=.TRUE.
      DOINT(2,2)=.FALSE.
C.....allocate array for holding integrals (on lower triangular form)
      allocate(PRPINT(NELMNT,NPRPINT))
      IORDER=0
      NPQUAD=0
C      
C     Find what symmetries are available
C
      allocate(IPSYM2INT(NBSYM,NOP))
      NCOMP=NPRPINT/2
      CALL ICOPY(2*NBSYM,-1,0,IPSYM2INT,1)
      DO I = 1,NCOMP
        ISYM = INTREP(I) + 1
        IPSYM2INT(ISYM,1) = I
        II = NCOMP + I
        ISYM = INTREP(II) + 1
        IPSYM2INT(ISYM,2) = II
      ENDDO
C
C     Operator setup
C
      IREP   =  JSYMOP-1
      ITIM   = -JTIMOP  ! Time antisymmetric operator
      NEXCNV = KEXCNV(JSYMOP)
      IF (NEXCNV .EQ. 0) RETURN
      IOPTYP  = 8                ! iA.P           * dot-product of alpha and vector operator
      NMAT    = MCMP(IOPTYP)     ! number of components, should be three
      IF((NPRPS+NOP).GT.MAXPRPS) THEN
        WRITE(LUPRI,'(A,I5)') 
     &   ' ** ERROR in GET_BED **  Redimension MAXPRPS to ',NPRPS+NOP
        CALL QUIT('GET_BED: Too may properties !')
      ENDIF
      DO IOP = 1,NOP
        WRITE(PRPNAM(NPRPS+IOP)(1:8),'(A3,I1,A1,A3)')
     &        'BED',IOP,'_',REP(IREP)
        IPRPSYM(NPRPS+IOP) = JSYMOP ! Total symmetry of operator
        IPRPTIM(NPRPS+IOP) = ITIM
        DO IMAT = 1,NMAT
          IAIND = JM4(IMAT,IOPTYP)
          IAREP = JM4REP(IAIND)
          IBSYM = IBTXOR(IAREP,IREP)+1
          IF(IPSYM2INT(IBSYM,IOP).EQ.-1) THEN
            PBUF(IMAT,IOP) = D0
          ELSE
            PBUF(IMAT,IOP)          = UPOL(IMAT) ! factor 0.5 not included
            IPRPLBL(IMAT,NPRPS+IOP) = IPSYM2INT(IBSYM,IOP)+NPRPLBL
          ENDIF
        ENDDO
        IF(DSUM(NMAT,PBUF(1,IOP),1).EQ.D0) THEN
           IPRPTYP(NPRPS+IOP) = -1  ! Means that there is no contribution
        ELSE
           IPRPTYP(NPRPS+IOP) = IOPTYP
        ENDIF
      ENDDO
      CALL PRSYMB(LUPRI,'=',80,0)
      WRITE(LUPRI,'(/A,I2,A,A3/)')
     &     ' *** Excitations of boson symmetry',JSYMOP,' : ',REP(IREP)
C
C     Prepare solution vectors
C     ========================
C
C      - orbital (e-e) part      
C
      NPAR = 0
      IF(NZXOPE.GT.0) THEN
         allocate(XPOE(NZXOPEQ,NEXCNV))
         allocate(BBUF(NZXOPEQ))
         allocate(IVECS(NESIM))
         NPAR = NPAR + 1
         TPAR(NPAR) = ' e-e '
         OPEN(LUBOE,FILE='PAMBOE',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPEQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBENDX,XPOE,EVECR,NEXCNV,IBTYP,
     &        IBEVC,IVECS,BBUF)
         CLOSE(LUBOE,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
C
C     - orbital (e-p) part      
C
      IF(NZXOPP.GT.0) THEN
         allocate(XPOP(NZXOPPQ,NEXCNV))
         allocate(BBUF(NZXOPPQ))
         allocate(IVECS(NPSIM))         
         NPAR = NPAR + 1
         TPAR(NPAR) = ' e-p '
         OPEN(LUBOP,FILE='PAMBOP',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPPQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBPNDX,XPOP,EVECR,NEXCNV,IBTYP,
     &        IBPVC,IVECS,BBUF)
         CLOSE(LUBOP,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
      TPAR(NPAR+1) = 'Total'
C      
C     Loop over excitation energies
C     =============================
C
      allocate(TMOM(NOP,NPAR+1,NEXCNV))
      allocate(OSC(NPAR+1,NEXCNV))
      IF(NZXOPE.GT.0) allocate(GPOE(NZXOPEQ,NOP))
      IF(NZXOPP.GT.0) allocate(GPOP(NZXOPPQ,NOP))      
      IF(DOBEDDER) THEN
         IH = 4
         H  = BEDSTP
         allocate(TMOM_P(NOP,-IH:IH))
         allocate(TMOM_DER(NOP,4,3,NEXCNV))
         allocate(TMOM_DOM(NOP,4))
         allocate(TMOM_DE(NOP,4))                  
         allocate(OSC_P(-IH:IH))
         allocate(OSC_DOM(4))
         allocate(OSC_DE(4))
         allocate(OSC_DER(4,3,NEXCNV))
      ELSE
         IH = 0
      ENDIF
      WRITE(6,'(A18,10X,A9,8X,A5)')
     &   'Transition moments','Freq(au)','Value'
      DO 10 IEXC = 1,NEXCNV
        IF(ABS(OMEGA(IEXC)).LT.1.0D-3) THEN
           WRITE(LUPRI,*)
     &     'GET_BED: Skipping excitation. Frequency too small: ',
     &     OMEGA(IEXC)
           GOTO 10
        ENDIF
C
C       Length of wave vector
C 
        FAC=OMEGA(IEXC)/CVAL
        DO IOP = 1,2
           WRITE(PRPNAM(NPRPS+IOP)(9:16),'(A4,I4)') '_exc',IEXC
           DO IMAT = 1,NMAT
              FACPRP(IMAT,NPRPS+IOP) = PBUF(IMAT,IOP)/FAC
           ENDDO
        ENDDO
C
C       Debug output
C
        IF(IPRXPP.GE.2) THEN
          DO IOP = 1,NOP
            WRITE(6,*) PRPNAM(NPRPS+IOP),IPRPTYP(NPRPS+IOP),FAC
            IF(IPRPTYP(NPRPS+IOP).EQ.-1) THEN
              WRITE(6,*) ' - skipped'
            ELSE
              DO IMAT = 1,NMAT
                II = IPRPLBL(IMAT,NPRPS+IOP)
                WRITE(6,*) PRPLBL(II),FACPRP(IMAT,NPRPS+IOP)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
        DO K = -IH,IH
C     
C         Generate integrals
C
          IF(BEDCHK) THEN
             CALL DCOPY(3,D0,0,WAVEVEC,1)
          ELSEIF(DOBEDDER) THEN
            DO I = 1,3
              WAVEVEC(I) = K*H*FAC*UWAVE(I)
            ENDDO
          ELSE
            DO I = 1,3
              WAVEVEC(I) = FAC*UWAVE(I)
            ENDDO
          ENDIF
C       
          CALL PR1IN2(WORK(KFREE),LFREE,PRPINT,INTTYP,
     &       NPRPINT,INTREP,ANTI,SQUARE,DUMMY,LABINT,TRIANG,
     &       TRASPH,.FALSE.,SOLVNT,0,0,DUMMY,
     &       0,NBAST,NELMNT,IPRINT,DOINT,.FALSE.)
          IF(IPRXPP.GE.5) THEN
            DO ICOMP = 1,NPRPINT
              WRITE(6,*) '* Integrals of ',LABINT(ICOMP),
     &        '. Symmetry: ',REP(INTREP(ICOMP))
              CALL OUTPAK(PRPINT(1,ICOMP),NBAST,1,LUPRI)
            ENDDO
          ENDIF
C
          OSC(NPAR+1,IEXC)    = D0
          TMOM(1,NPAR+1,IEXC) = D0
          TMOM(2,NPAR+1,IEXC) = D0
          IPAR = 0
C
C         Gradient - orbital (e-e) part
C         ====================================
          IF(NZXOPE.GT.0) THEN
            IPAR = IPAR + 1
            TMP = D0
            DO IOP = 1,NOP
              INDXPR = NPRPS+IOP
              IF(IPRPTYP(INDXPR).NE.-1) THEN
                CALL GPGET(INDXPR,GPOE(1,IOP),
     &               get_orbital_rotation_indices_pp(),
     &               NZXOPE,CMO,IBEIG,JBENDX,.FALSE.,PRPINT,
     &               WORK,KFREE,LFREE,IPRXPP)
                TMOM(IOP,IPAR,IEXC)=D2*
     &               DDOT(NZXOPEQ,GPOE(1,IOP),1,XPOE(1,IEXC),1)
              ELSE
                TMOM(IOP,IPAR,IEXC)=D0
              ENDIF
              TMOM(IOP,NPAR+1,IEXC)=TMOM(IOP,NPAR+1,IEXC)
     &                             +TMOM(IOP,IPAR,IEXC)
              TMP = TMP + TMOM(IOP,IPAR,IEXC)*TMOM(IOP,IPAR,IEXC)
            ENDDO
            OSC(IPAR,IEXC) = D2*OMEGA(IEXC)*TMP
            OSC(NPAR+1,IEXC) = OSC(NPAR+1,IEXC) + OSC(IPAR,IEXC)
          ENDIF
C
C         Gradient - orbital (e-p) part
C         ====================================
          IF(NZXOPP.GT.0) THEN       
            IPAR = IPAR + 1
            TMP = D0
            DO IOP = 1,NOP
              INDXPR = NPRPS+IOP
              IF(IPRPTYP(INDXPR).NE.-1) THEN
                CALL GPGET(INDXPR,GPOP(1,IOP),
     &               get_orbital_rotation_indices_pn(),
     &               NZXOPP,CMO,IBEIG,JBPNDX,.FALSE.,PRPINT,
     &               WORK,KFREE,LFREE,IPRXPP)
                TMOM(IOP,IPAR,IEXC)=D2*
     &               DDOT(NZXOPPQ,GPOP(1,IOP),1,XPOP(1,IEXC),1)
              ELSE
                TMOM(IOP,IPAR,IEXC)=D0
              ENDIF
              TMOM(IOP,NPAR+1,IEXC)=TMOM(IOP,NPAR+1,IEXC)
     &                             +TMOM(IOP,IPAR,IEXC)
              TMP = TMP + TMOM(IOP,IPAR,IEXC)*TMOM(IOP,IPAR,IEXC)
            ENDDO
            OSC(IPAR,IEXC)   = D2*OMEGA(IEXC)*TMP
            OSC(NPAR+1,IEXC) = OSC(NPAR+1,IEXC) + OSC(IPAR,IEXC)
          ENDIF
          IF(DOBEDDER) THEN
            DO IOP = 1,NOP
              TMOM_P(IOP,K) = TMOM(IOP,NPAR+1,IEXC)
           ENDDO
           OSC_P(K) = OSC(NPAR+1,IEXC)           
          ENDIF
        ENDDO
        IF(DOBEDDER) THEN
C.........Prepare derivatives of transition moments
          DO IOP = 1,NOP
            TMOM_DOM(IOP,1)=(TMOM_P(IOP,+1)-TMOM_P(IOP,-1))/D2
            TMOM_DOM(IOP,2)=((TMOM_P(IOP,+2)-TMOM_P(IOP,-2))
     &          -2.0D0*(TMOM_P(IOP,+1)-TMOM_P(IOP,-1)))/D2
            TMOM_DOM(IOP,3)=((TMOM_P(IOP,+3)-TMOM_P(IOP,-3))
     &           -4.0D0*(TMOM_P(IOP,+2)-TMOM_P(IOP,-2))
     &           +5.0D0*(TMOM_P(IOP,+1)-TMOM_P(IOP,-1)))/D2
            TMOM_DOM(IOP,4)=((TMOM_P(IOP,+4)-TMOM_P(IOP,-4))
     &           -6.0D0*(TMOM_P(IOP,+3)-TMOM_P(IOP,-3))
     &           +14.0D0*(TMOM_P(IOP,+2)-TMOM_P(IOP,-2))
     &           -14.0D0*(TMOM_P(IOP,+1)-TMOM_P(IOP,-1)))/D2
            TMOM_DE(IOP,1)=(TMOM_P(IOP,+1)+TMOM_P(IOP,-1))
     &           -2.0D0*TMOM_P(IOP,0)
            TMOM_DE(IOP,2)=(TMOM_P(IOP,+2)+TMOM_P(IOP,-2))
     &           -4.0D0*(TMOM_P(IOP,+1)+TMOM_P(IOP,-1))
     &           +6.0D0*TMOM_P(IOP,0)
            TMOM_DE(IOP,3)=(TMOM_P(IOP,+3)+TMOM_P(IOP,-3))
     &           -6.0D0*(TMOM_P(IOP,+2)+TMOM_P(IOP,-2))
     &           +15.0D0*(TMOM_P(IOP,+1)+TMOM_P(IOP,-1))
     &           -20.0D0*TMOM_P(IOP,0)
            TMOM_DE(IOP,4)=(TMOM_P(IOP,+4)+TMOM_P(IOP,-4))
     &           -8.0D0*(TMOM_P(IOP,+3)+TMOM_P(IOP,-3))
     &           +28.0D0*(TMOM_P(IOP,+2)+TMOM_P(IOP,-2))
     &           -56.0D0*(TMOM_P(IOP,+1)+TMOM_P(IOP,-1))
     &           +70.0*TMOM_P(IOP,0)
            WRITE(6,*) '* Transition moments --- tabulated points:',
     &        IEXC,OMEGA(IEXC),IOP
            DO K = -IH,IH
               WRITE(6,*) K*H*FAC,TMOM_P(IOP,K)
            ENDDO
C...........Zeroth order            
            TMOM(IOP,NPAR+1,IEXC)=TMOM_P(IOP,0)
C...........First derivative
            TFAC=D1/H
            DFAC=TFAC
            TMOM_DER(IOP,1,1,IEXC)=TMOM_DOM(IOP,1)
            TMOM_DER(IOP,1,2,IEXC)=TMOM_DER(IOP,1,1,IEXC)
     &                            -TMOM_DOM(IOP,2)/6.0D0
            TMOM_DER(IOP,1,3,IEXC)=TMOM_DER(IOP,1,2,IEXC)
     &           +TMOM_DOM(IOP,3)/30.0D0
            DO I = 1,3
              TMOM_DER(IOP,1,I,IEXC)=DFAC*TMOM_DER(IOP,1,I,IEXC)
            ENDDO
C...........Second derivative
            DFAC=DFAC*TFAC/2.0D0
            TMOM_DER(IOP,2,1,IEXC)=TMOM_DE(IOP,1)
            TMOM_DER(IOP,2,2,IEXC)=TMOM_DER(IOP,2,1,IEXC)
     &                            -TMOM_DE(IOP,2)/12.0D0
            TMOM_DER(IOP,2,3,IEXC)=TMOM_DER(IOP,2,2,IEXC)
     &                            +TMOM_DE(IOP,3)/90.0D0
            DO I = 1,3
              TMOM_DER(IOP,2,I,IEXC)=DFAC*TMOM_DER(IOP,2,I,IEXC)
            ENDDO
C...........Third derivative
            DFAC=DFAC*TFAC/3.0D0
            TMOM_DER(IOP,3,1,IEXC)=TMOM_DOM(IOP,2)
            TMOM_DER(IOP,3,2,IEXC)=TMOM_DER(IOP,3,1,IEXC)
     &                            -TMOM_DOM(IOP,3)/4.0D0
            TMOM_DER(IOP,3,3,IEXC)=TMOM_DER(IOP,3,2,IEXC)
     &                            +TMOM_DOM(IOP,4)*7.0D0/120.0D0
            DO I = 1,3
              TMOM_DER(IOP,3,I,IEXC)=DFAC*TMOM_DER(IOP,3,I,IEXC)
            ENDDO
C...........Fourth derivative
            DFAC=DFAC*TFAC/4.0D0
            TMOM_DER(IOP,4,1,IEXC)=TMOM_DE(IOP,2)
            TMOM_DER(IOP,4,2,IEXC)=TMOM_DER(IOP,4,1,IEXC)
     &                            -TMOM_DE(IOP,3)/6.0D0
            TMOM_DER(IOP,4,3,IEXC)=TMOM_DER(IOP,4,2,IEXC)
     &                            +TMOM_DE(IOP,4)*7.0D0/240.0D0
            DO I = 1,3
              TMOM_DER(IOP,4,I,IEXC)=DFAC*TMOM_DER(IOP,4,I,IEXC)
            ENDDO
          ENDDO
C.........Prepare derivatives of oscillator strenghts         
          OSC_DOM(1)=(OSC_P(+1)-OSC_P(-1))/D2
          OSC_DOM(2)=((OSC_P(+2)-OSC_P(-2))
     &              -2.0D0*(OSC_P(+1)-OSC_P(-1)))/D2
          OSC_DOM(3)=((OSC_P(+3)-OSC_P(-3))-4.0D0*(OSC_P(+2)-OSC_P(-2))
     &              +5.0D0*(OSC_P(+1)-OSC_P(-1)))/D2
          OSC_DOM(4)=((OSC_P(+4)-OSC_P(-4))-6.0D0*(OSC_P(+3)-OSC_P(-3))
     &              +14.0D0*(OSC_P(+2)-OSC_P(-2))
     &              -14.0D0*(OSC_P(+1)-OSC_P(-1)))/D2
          OSC_DE(1)=(OSC_P(+1)+OSC_P(-1))-2.0D0*OSC_P(0)
          OSC_DE(2)=(OSC_P(+2)+OSC_P(-2))-4.0D0*(OSC_P(+1)+OSC_P(-1))
     &             +6.0D0*OSC_P(0)
          OSC_DE(3)=(OSC_P(+3)+OSC_P(-3))-6.0D0*(OSC_P(+2)+OSC_P(-2))
     &             +15.0D0*(OSC_P(+1)+OSC_P(-1))
     &             -20.0D0*OSC_P(0)
          OSC_DE(4)=(OSC_P(+4)+OSC_P(-4))-8.0D0*(OSC_P(+3)+OSC_P(-3))
     &             +28.0D0*(OSC_P(+2)+OSC_P(-2))
     &             -56.0D0*(OSC_P(+1)+OSC_P(-1))+70.0*OSC_P(0)
            WRITE(6,*) '* Oscillator strengths --- tabulated points:',
     &        IEXC,OMEGA(IEXC)
            DO K = -IH,IH
               WRITE(6,*) K*H*FAC,OSC_P(K)
            ENDDO
C.........Zeroth order
          OSC(NPAR+1,IEXC)= OSC_P(0)
C.........First derivative
          TFAC=D1/H
          DFAC=TFAC
          OSC_DER(1,1,IEXC)=OSC_DOM(1)
          OSC_DER(1,2,IEXC)=OSC_DER(1,1,IEXC)-OSC_DOM(2)/6.0D0
          OSC_DER(1,3,IEXC)=OSC_DER(1,2,IEXC)+OSC_DOM(3)/30.0D0
          DO I = 1,3
            OSC_DER(1,I,IEXC)=DFAC*OSC_DER(1,I,IEXC)
          ENDDO
C.........Second derivative
          DFAC=DFAC*TFAC/2.0D0
          OSC_DER(2,1,IEXC)=OSC_DE(1)
          OSC_DER(2,2,IEXC)=OSC_DER(2,1,IEXC)-OSC_DE(2)/12.0D0
          OSC_DER(2,3,IEXC)=OSC_DER(2,2,IEXC)+OSC_DE(3)/90.0D0
          DO I = 1,3
            OSC_DER(2,I,IEXC)=DFAC*OSC_DER(2,I,IEXC)
          ENDDO
C.........Third derivative
          DFAC=DFAC*TFAC/3.0D0
          OSC_DER(3,1,IEXC)=OSC_DOM(2)
          OSC_DER(3,2,IEXC)=OSC_DER(3,1,IEXC)-OSC_DOM(3)/4.0D0
          OSC_DER(3,3,IEXC)=OSC_DER(3,2,IEXC)+OSC_DOM(4)*7.0D0/120.0D0
          DO I = 1,3
            OSC_DER(3,I,IEXC)=DFAC*OSC_DER(3,I,IEXC)
          ENDDO
C.........Fourth derivative
          DFAC=DFAC*TFAC/4.0D0
          OSC_DER(4,1,IEXC)=OSC_DE(2)
          OSC_DER(4,2,IEXC)=OSC_DER(4,1,IEXC)-OSC_DE(3)/6.0D0
          OSC_DER(4,3,IEXC)=OSC_DER(4,2,IEXC)+OSC_DE(4)*7.0D0/240.0D0
          DO I = 1,3
            OSC_DER(4,I,IEXC)=DFAC*OSC_DER(4,I,IEXC)
          ENDDO
       ENDIF
 10   CONTINUE
C
C     Print section
C
      IF(DOBEDDER) THEN
         WRITE(6,*)
     &   '* Transition moments expansion --- real part (sin [kr]):'
         DO IP = 1,3
           WRITE(6,*) '   Stencil: ',IP
           WRITE(6,'(5X,A,14X,6(A6,9X))') 'Freq.','Tmom0 ','Tmom1 ',
     &        'Tmom2 ','Tmom3 ','Tmom4 ','Ttotal'
           DO IEXC = 1,NEXCNV
             TMP=TMOM(2,NPAR+1,IEXC)
             DO IDER = 1,4
               TMP = TMP + TMOM_DER(2,IDER,IP,IEXC)
             ENDDO
           WRITE(6,'(I4,1X,F18.10,6(1X,E14.8))')
     &            IEXC,OMEGA(IEXC),TMOM(2,NPAR+1,IEXC),
     &           (TMOM_DER(2,IDER,IP,IEXC),IDER=1,4),TMP
           ENDDO
         ENDDO
         WRITE(6,*)
     &   '* Transition moments expansion --- imag. part (cos[kr]) :'
         DO IP = 1,3
           WRITE(6,*) '   Stencil: ',IP
           WRITE(6,'(5X,A,14X,6(A6,9X))') 'Freq.','Tmom0 ','Tmom1 ',
     &        'Tmom2 ','Tmom3 ','Tmom4 ','Ttotal'
           DO IEXC = 1,NEXCNV
             TMP=TMOM(1,NPAR+1,IEXC)
             DO IDER = 1,4
               TMP = TMP + TMOM_DER(1,IDER,IP,IEXC)
             ENDDO
             WRITE(6,'(I4,1X,F18.10,6(1X,E14.8))')
     &            IEXC,OMEGA(IEXC),-TMOM(1,NPAR+1,IEXC),
     &           (-TMOM_DER(1,IDER,IP,IEXC),IDER=1,4),-TMP
           ENDDO
         ENDDO
         WRITE(6,*) '* Oscillator strengths expansion:'
           DO IP = 1,3
             WRITE(6,*) '   Stencil: ',IP
             WRITE(6,'(5X,A,14X,6(A6,9X))') 'Freq.','Osc0  ','Osc1  ',
     &        'Osc2  ','Osc3  ','Osc4  ','Ototal'
             DO IEXC = 1,NEXCNV
               TMP=OSC(NPAR+1,IEXC)
               DO IDER = 1,4
                 TMP = TMP + OSC_DER(IDER,IP,IEXC)
             ENDDO
             WRITE(6,'(I4,1X,F18.10,6(1X,E14.8))')
     &            IEXC,OMEGA(IEXC),OSC(NPAR+1,IEXC),
     &           (OSC_DER(IDER,IP,IEXC),IDER=1,4),TMP
           ENDDO
         ENDDO
       ELSE
         WRITE(6,*) '* Transition moments:'
         WRITE(6,'(A4,1X,A18,14X,3(3X,A5,9X))')     
     &  'Exc.','Frequency (au)    ',(TPAR(IPAR),IPAR=1,NPAR+1)
         DO IEXC = 1,NEXCNV
           TMP=MAX(ABS(TMOM(1,NPAR+1,IEXC)),ABS(TMOM(2,NPAR+1,IEXC)))
           FMT=MXFORM(TMP,14)
           WRITE(6,'(I4,1X,F18.10,A,3(3X,'//FMT//'))')
     &     IEXC,OMEGA(IEXC),' real part (sin [kr]): ',
     &     (TMOM(2,IPAR,IEXC),IPAR = 1,NPAR+1)
           WRITE(6,'(23X,A,3(3X,'//FMT//'))')
     &     ' imag. part (cos[kr]) : ',
     &     (-TMOM(1,IPAR,IEXC),IPAR = 1,NPAR+1)
        ENDDO
        WRITE(6,*) '* Oscillator strengths:'
        WRITE(6,'(A4,1X,A18,3(3X,A5,9X))')     
     &  'Exc.','Frequency (au)    ',(TPAR(IPAR),IPAR=1,NPAR+1)
        DO IEXC = 1,NEXCNV
          WRITE(6,'(I4,1X,F18.10,3(3X,'//FMT//'))')
     &     IEXC,OMEGA(IEXC),
     &     (OSC(IPAR,IEXC),IPAR = 1,NPAR+1)
        ENDDO
      ENDIF
C     
C     Memory deallocation
      IF(NZXOPE.GT.0) THEN
        deallocate(GPOE)
        deallocate(XPOE)
      ENDIF
      IF(NZXOPP.GT.0) THEN
        deallocate(GPOP)
        deallocate(XPOP)         
      ENDIF
      deallocate(TMOM)
      deallocate(OSC)
      IF(DOBEDDER) THEN
         deallocate(TMOM_P)
         deallocate(TMOM_DER)                           
         deallocate(TMOM_DOM)
         deallocate(TMOM_DE)
         deallocate(OSC_P)
         deallocate(OSC_DER)
         deallocate(OSC_DOM)
         deallocate(OSC_DE)                  
      ENDIF
      deallocate(PRPINT)
      deallocate(LABINT)
      deallocate(INTREP)
      deallocate(IPSYM2INT)
      deallocate(first)
C     
      RETURN
      END
C
C  /* Deck geted_length */                                                  
      SUBROUTINE GET_ED_LENGTH(OMEGA,ISSYM,CMO,IBEIG,EVECR,IBTYP,
     &                   IBCVC,IBEVC,IBPVC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Test routine for electric dipole approximation in length gauge.
C     Scalar integrals are sorted on symmetries (max. 8)
C
C     Upon entering this routine solution vectors have been written
C     to unformatted file PAMXVC, to be opened with unit number LUXVC
C     References are:
C     1. Nanna Holmgaard List, Joanna Kauczor, Trond Saue,
C        Hans Jørgen Aagaard Jensen and Patrick Norman, J. Chem. Phys. 142 (2015) 244111
C                               
C                                                                      
C***********************************************************************
      use orbital_rotation_indices
C     
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0=0.0D0,D2=2.0D0,D1=1.0D0)
#include "dcbgen.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "pgroup.h"
#include "dgroup.h"
#include "nuclei.h"      
#include "dcbxrs.h"
#include "dcbxpp.h"
#include "orgcom.h"
#include "dcbxpr.h"
#include "dcbprl.h"            
#include "dcbibn.h"
#include "dcbwav.h"
C#include "cbiher.h"
      LOGICAL DOINT(2,2)
      INTEGER LWORK, NCOMP, NOPTP,IBEIG(*)
      DIMENSION OMEGA(MAXEXC),ISSYM(MAXEXC),CMO(*),WORK(*)
      real(8), allocatable              :: PRPINT(:,:),BBUF(:)
      character (len=8),  allocatable   :: LABINT(:)
      integer, allocatable              :: INTREP(:)
      integer, allocatable              :: IVECS(:)
      real(8)                           :: EVECR(NREDM,*)
      integer             :: IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*)
      real(8), allocatable              :: GPOE(:,:),GPOP(:,:)
      real(8), allocatable              :: XPOE(:,:),XPOP(:,:)
      logical, allocatable              :: first(:),doatom(:)
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE
C
#include "ibtfun.h"
C     
C     Allocate placeholders for integral labels and irreps 
C
      allocate(LABINT(3))
      allocate(INTREP(3))
      allocate(first(nz))
C
C     Integral directives
C
      IPRINT=0
      allocate(doatom(nucind))
      CALL PR1DIR('DIPLEN  ',INTTYP,NPRPINT,INTREP,
     &            ANTI,SQUARE,DUMMY,LABINT,TRIANG,TRASPH,
     &            SOLVNT,0,DOATOM,NUCIND,NBAST,NELMNT,IPRINT)
C     NPRPINT is by construction an even number,
C     the first and second halfs corresponding to the real (eq. 44) and imaginary parts (eq. 45)
      IF(NPRPLBL+NPRPINT.GT.MAXPRPLBL) THEN
        CALL QUIT('GET_BED:Redimension MAXPRPLBL !')
      ENDIF
      DO ICOMP = 1,NPRPINT
        PRPLBL (NPRPLBL+ICOMP) = LABINT(ICOMP)
        IPRLREP(NPRPLBL+ICOMP) = INTREP(ICOMP)
        PDOINT (NPRPLBL+ICOMP) = '+00+'
        IPRLTYP(NPRPLBL+ICOMP) =  1        
      ENDDO
      DOINT(1,1)=.TRUE.
      DOINT(2,1)=.FALSE.
      DOINT(1,2)=.FALSE.
      DOINT(2,2)=.TRUE.
C.....allocate array for holding integrals (on lower triangular form)
      allocate(PRPINT(NELMNT,NPRPINT))
      IORDER=0
      NPQUAD=0
C
C     Operator setup
C
      NOP = 3
      IREP   = JSYMOP-1
      NEXCNV = KEXCNV(JSYMOP)
      IF (NEXCNV .EQ. 0) RETURN
      IOPTYP  = 1                ! iA.P           * dot-product of alpha and vector operator
      NMAT    = MCMP(IOPTYP)     ! number of components, should be three
      IF((NPRPS+NOP).GT.MAXPRPS) THEN
        WRITE(LUPRI,'(A,I5)') 
     &        ' ** ERROR in GET_ED_LENGTH **  Redimension MAXPRPS to ',
     &        NPRPS+NOP
        CALL QUIT('GET_ED_LENGTH: Too may properties !')
      ENDIF
      DO IOP = 1,NOP
        WRITE(PRPNAM(NPRPS+IOP),'(A4,A3,A1,A8)')
     &        'ED__',REP(IREP),'_',LABINT(IOP)
        IPRPSYM(NPRPS+IOP) = JSYMOP ! Total symmetry of operator
        IPRPTIM(NPRPS+IOP) = JTIMOP ! Time antisymmetric operator
        IF(INTREP(IOP).EQ.IREP) THEN
           IPRPTYP(NPRPS+IOP)   = IOPTYP
           IPRPLBL(1,NPRPS+IOP) = IOP + NPRPLBL
           FACPRP(1,NPRPS+IOP)  = D1
        ELSE
           IPRPTYP(NPRPS+IOP) = -1
           FACPRP(1,NPRPS+IOP) = D0
        ENDIF
      ENDDO
      CALL PRSYMB(LUPRI,'=',80,0)
      WRITE(LUPRI,'(/A,I2,A,A3/)')
     &     ' *** Excitations of boson symmetry',JSYMOP,' : ',REP(IREP)
C
C     Prepare solution vectors
C     ========================
C
C      - orbital (e-e) part      
C
      ITIM = JTIMOP
      IF(NZXOPE.GT.0) THEN
         allocate(XPOE(NZXOPEQ,NEXCNV))
         allocate(BBUF(NZXOPEQ))
         allocate(IVECS(NESIM))         
         OPEN(LUBOE,FILE='PAMBOE',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPEQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBENDX,XPOE,EVECR,NEXCNV,IBTYP,
     &        IBEVC,IVECS,BBUF)
         CLOSE(LUBOE,STATUS = 'KEEP')
         deallocate(BBUF)

         deallocate(IVECS)         
      ENDIF
C
C     - orbital (e-e) part      
C
      IF(NZXOPP.GT.0) THEN
         allocate(XPOP(NZXOPPQ,NEXCNV))
         allocate(BBUF(NZXOPPQ))
         allocate(IVECS(NPSIM))         
         OPEN(LUBOP,FILE='PAMBOP',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPPQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBPNDX,XPOP,EVECR,NEXCNV,IBTYP,
     &        IBPVC,IVECS,BBUF)
         CLOSE(LUBOP,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
C      
C     Loop over excitation energies
C     =============================
C     
      IF(NZXOPE.GT.0) allocate(GPOE(NZXOPEQ,NOP))
      IF(NZXOPP.GT.0) allocate(GPOP(NZXOPPQ,NOP))      
      WRITE(6,'(A18,10X,A9,8X,A5)')
     &   'Transition moments','Freq(au)','Value'
      DO 10 IEXC = 1,KEXCNV(JSYMOP)
        IF(ABS(OMEGA(IEXC)).LT.1.0D-3) THEN
           WRITE(LUPRI,*)
     &     'GET_ED_LENGTH: Skipping excitation. Frequency too small: ',
     &     OMEGA(IEXC)
           GOTO 10
        ENDIF
C
C       Debug output
C
        IF(IPRXPP.GE.2) THEN
          DO IOP = 1,NOP
            WRITE(6,*) PRPNAM(NPRPS+IOP),IPRPTYP(NPRPS+IOP)
            IF(IPRPTYP(NPRPS+IOP).EQ.-1) THEN
              WRITE(6,*) ' - skipped'
            ELSE
              DO IMAT = 1,1
                II = IPRPLBL(IMAT,NPRPS+IOP)
                WRITE(6,*) PRPLBL(II),FACPRP(IMAT,NPRPS+IOP)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
C     
C       Generate integrals
C
        CALL PR1IN2(WORK(KFREE),LFREE,PRPINT,INTTYP,
     &       NPRPINT,INTREP,ANTI,SQUARE,DUMMY,LABINT,TRIANG,
     &       TRASPH,.FALSE.,SOLVNT,0,0,DOATOM,
     &       NUCIND,NBAST,NELMNT,IPRINT,DOINT,.FALSE.)
        IF(IPRXPP.GE.5) THEN
          DO ICOMP = 1,NPRPINT
            WRITE(6,*) '* Integrals of ',LABINT(ICOMP),
     &        '. Symmetry: ',REP(INTREP(NCOMP))
            CALL OUTPAK(PRPINT(1,ICOMP),NBAST,1,LUPRI)
          ENDDO
        ENDIF
C
C     Gradient - orbital (e-e) part
C     ====================================
        IF(NZXOPE.GT.0) THEN       
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
               CALL GPGET(INDXPR,GPOE(1,IOP),
     &              get_orbital_rotation_indices_pp(),
     &              NZXOPE,CMO,IBEIG,JBENDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
               TMP=D2*DDOT(NZXOPEQ,GPOE(1,IOP),1,XPOE(1,IEXC),1)
               WRITE(6,'(A,1X,A16,F18.10,3X,F14.10)')
     &            'e-e :',PRPNAM(NPRPS+IOP),OMEGA(IEXC),TMP
            ENDIF
          ENDDO
        ENDIF
C
C     Gradient - orbital (e-p) part
C     ====================================
        IF(NZXOPP.GT.0) THEN       
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
               CALL GPGET(INDXPR,GPOP(1,IOP),
     &              get_orbital_rotation_indices_pn(),
     &              NZXOPP,CMO,IBEIG,JBPNDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
               TMP=D2*DDOT(NZXOPPQ,GPOP(1,IOP),1,XPOP(1,IEXC),1)
               WRITE(6,'(A,1X,A16,F18.10,3X,F14.10)')
     &            'e-p :',PRPNAM(NPRPS+IOP),OMEGA(IEXC),TMP
            ENDIF
          ENDDO
        ENDIF
 10   CONTINUE
C
C     Memory deallocation
      IF(NZXOPE.GT.0) THEN
        deallocate(GPOE)
        deallocate(XPOE)
      ENDIF
      IF(NZXOPP.GT.0) THEN
        deallocate(GPOP)
        deallocate(XPOP)         
      ENDIF
      deallocate(PRPINT)
      deallocate(LABINT)
      deallocate(INTREP)
      deallocate(first)
      deallocate(doatom)
C     
      RETURN
      END
C  /* Deck geted_length */                                                  
      SUBROUTINE GET_ED_VELOCITY(OMEGA,ISSYM,CMO,IBEIG,EVECR,IBTYP,
     &                   IBCVC,IBEVC,IBPVC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Test routine for electric dipole approximation in velocity gauge.
C     Scalar integrals are sorted on symmetries (max. 8)
C
C     Upon entering this routine solution vectors have been written
C     to unformatted file PAMXVC, to be opened with unit number LUXVC
C     References are:
C     1. Nanna Holmgaard List, Joanna Kauczor, Trond Saue,
C        Hans Jørgen Aagaard Jensen and Patrick Norman, J. Chem. Phys. 142 (2015) 244111
C                               
C                                                                      
C***********************************************************************
      use orbital_rotation_indices
C     
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0=0.0D0,D2=2.0D0,D1=1.0D0)
#include "dcbgen.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "chrxyz.h"
#include "pgroup.h"
#include "dgroup.h"
#include "nuclei.h"      
#include "dcbxrs.h"
#include "dcbxpp.h"
#include "orgcom.h"
#include "dcbxpr.h"
#include "dcbprl.h"            
#include "dcbibn.h"
#include "dcbwav.h"
C#include "cbiher.h"
      LOGICAL DOINT(2,2)
      INTEGER LWORK, NCOMP, NOPTP,IBEIG(*)
      DIMENSION OMEGA(MAXEXC),ISSYM(MAXEXC),CMO(*),WORK(*)
      real(8), allocatable              :: PRPINT(:,:),BBUF(:)
      character (len=8),  allocatable   :: LABINT(:)
      integer, allocatable              :: INTREP(:)
      integer, allocatable              :: IVECS(:)
      real(8)                           :: EVECR(NREDM,*)
      integer             :: IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*)
      real(8), allocatable              :: GPOE(:,:),GPOP(:,:)
      real(8), allocatable              :: XPOE(:,:),XPOP(:,:)
      logical, allocatable              :: first(:),doatom(:)
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE
C
#include "ibtfun.h"
      NOP = 3
      NEXCNV = KEXCNV(JSYMOP)
      IF (NEXCNV .EQ. 0) RETURN
C     
C     Allocate placeholders for integral labels and irreps 
C
      allocate(LABINT(1))
      allocate(INTREP(1))
      allocate(first(nz))
C
C     Integral directives
C
      IPRINT=0
      allocate(doatom(nucind))
      CALL PR1DIR('OVERLAP ',INTTYP,NPRPINT,INTREP,
     &            ANTI,SQUARE,DUMMY,LABINT,TRIANG,TRASPH,
     &            SOLVNT,0,DOATOM,NUCIND,NBAST,NELMNT,IPRINT)
      IF(NPRPLBL+NPRPINT.GT.MAXPRPLBL) THEN
        CALL QUIT('GET_ED_VELOCITY:Redimension MAXPRPLBL !')
      ENDIF
      DO ICOMP = 1,NPRPINT
        PRPLBL (NPRPLBL+ICOMP) = LABINT(ICOMP)
        IPRLREP(NPRPLBL+ICOMP) = INTREP(ICOMP)
        PDOINT (NPRPLBL+ICOMP) = '0++0'
        IPRLTYP(NPRPLBL+ICOMP) =  1        
      ENDDO
      DOINT(1,1)=.FALSE.
      DOINT(2,1)=.TRUE.
      DOINT(1,2)=.TRUE.
      DOINT(2,2)=.FALSE.
C.....allocate array for holding integrals (on lower triangular form)
      allocate(PRPINT(NELMNT,NPRPINT))
      IORDER=0
      NPQUAD=0
C
C     Operator setup
C
      NOP    = 3
      IREP   = JSYMOP-1
      ITIM   = -JTIMOP ! Time antisymmetric operator
      IF((NPRPS+NOP).GT.MAXPRPS) THEN
        WRITE(LUPRI,'(A,A,I5)') 
     &        ' ** ERROR in GET_ED_VELOCITY **  ',
     &        'Redimension MAXPRPS to ',
     &        NPRPS+NOP
        CALL QUIT('GET_ED_LENGTH: Too may properties !')
      ENDIF
      DO IOP = 1,NOP
        WRITE(PRPNAM(NPRPS+IOP),'(A4,A3,2A1,A7)')
     &        'ED__',REP(IREP),'_',CHRXYZ(IOP),'dipvel_'
        IOPTYP             = 1 + IOP
        IAIND = JM4(1,IOPTYP)
        IAREP = JM4REP(IAIND)
        IPRPSYM(NPRPS+IOP) = JSYMOP ! Total symmetry of operator
        IPRPTIM(NPRPS+IOP) = ITIM   
        IF(IAREP.EQ.IREP) THEN
           IPRPTYP(NPRPS+IOP)   = IOPTYP
           IPRPLBL(1,NPRPS+IOP) = 1 + NPRPLBL
           FACPRP(1,NPRPS+IOP)  = CVAL
        ELSE
           IPRPTYP(NPRPS+IOP)  = -1
           FACPRP(1,NPRPS+IOP) = D0
        ENDIF
      ENDDO
      CALL PRSYMB(LUPRI,'=',80,0)
      WRITE(LUPRI,'(/A,I2,A,A3/)')
     &     ' *** Excitations of boson symmetry',JSYMOP,' : ',REP(IREP)
C
C     Prepare solution vectors
C     ========================
C
C      - orbital (e-e) part      
C
      IF(NZXOPE.GT.0) THEN
         allocate(XPOE(NZXOPEQ,NEXCNV))
         allocate(BBUF(NZXOPEQ))
         allocate(IVECS(NESIM))         
         OPEN(LUBOE,FILE='PAMBOE',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPEQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBENDX,XPOE,EVECR,NEXCNV,IBTYP,
     &        IBEVC,IVECS,BBUF)
         CLOSE(LUBOE,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
C
C     - orbital (e-e) part      
C
      IF(NZXOPP.GT.0) THEN
         allocate(XPOP(NZXOPPQ,NEXCNV))
         allocate(BBUF(NZXOPPQ))
         allocate(IVECS(NPSIM))         
         OPEN(LUBOP,FILE='PAMBOP',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPPQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBPNDX,XPOP,EVECR,NEXCNV,IBTYP,
     &        IBPVC,IVECS,BBUF)
         CLOSE(LUBOP,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
C      
C     Loop over excitation energies
C     =============================
C     
      IF(NZXOPE.GT.0) allocate(GPOE(NZXOPEQ,NOP))
      IF(NZXOPP.GT.0) allocate(GPOP(NZXOPPQ,NOP))      
      WRITE(6,'(A18,10X,A9,8X,A5)')
     &   'Transition moments','Freq(au)','Value'
      DO 10 IEXC = 1,KEXCNV(JSYMOP)
        IF(ABS(OMEGA(IEXC)).LT.1.0D-3) THEN
           WRITE(LUPRI,*)
     &     'GET_ED_LENGTH: Skipping excitation. Frequency too small: ',
     &     OMEGA(IEXC)
           GOTO 10
        ENDIF
C
C       Debug output
C
        IF(IPRXPP.GE.2) THEN
          DO IOP = 1,NOP
            WRITE(6,*) PRPNAM(NPRPS+IOP),IPRPTYP(NPRPS+IOP)
            IF(IPRPTYP(NPRPS+IOP).EQ.-1) THEN
              WRITE(6,*) ' - skipped'
            ELSE
              DO IMAT = 1,1
                II = IPRPLBL(IMAT,NPRPS+IOP)
                WRITE(6,*) PRPLBL(II),FACPRP(IMAT,NPRPS+IOP)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
C     
C       Generate integrals
C
        CALL PR1IN2(WORK(KFREE),LFREE,PRPINT,INTTYP,
     &       NPRPINT,INTREP,ANTI,SQUARE,DUMMY,LABINT,TRIANG,
     &       TRASPH,.FALSE.,SOLVNT,0,0,DOATOM,
     &       NUCIND,NBAST,NELMNT,IPRINT,DOINT,.FALSE.)
        IF(IPRXPP.GE.5) THEN
          DO ICOMP = 1,NPRPINT
            WRITE(6,*) '* Integrals of ',LABINT(ICOMP),
     &        '. Symmetry: ',REP(INTREP(ICOMP))
            CALL OUTPAK(PRPINT(1,ICOMP),NBAST,1,LUPRI)
          ENDDO
        ENDIF
C
C     Gradient - orbital (e-e) part
C     ====================================
        IF(NZXOPE.GT.0) THEN       
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
               CALL GPGET(INDXPR,GPOE(1,IOP),
     &              get_orbital_rotation_indices_pp(),
     &              NZXOPE,CMO,IBEIG,JBENDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
               TMP=D2*DDOT(NZXOPEQ,GPOE(1,IOP),1,XPOE(1,IEXC),1)
               WRITE(6,'(A,1X,A16,F18.10,3X,F14.10)')
     &            'e-e :',PRPNAM(NPRPS+IOP),OMEGA(IEXC),TMP
            ENDIF
          ENDDO
        ENDIF
C
C     Gradient - orbital (e-e) part
C     ====================================
        IF(NZXOPP.GT.0) THEN       
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
               CALL GPGET(INDXPR,GPOP(1,IOP),
     &              get_orbital_rotation_indices_pn(),
     &              NZXOPP,CMO,IBEIG,JBPNDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
               TMP=D2*DDOT(NZXOPPQ,GPOP(1,IOP),1,XPOP(1,IEXC),1)
               WRITE(6,'(A,1X,A16,F18.10,3X,F14.10)')
     &            'e-p :',PRPNAM(NPRPS+IOP),OMEGA(IEXC),TMP
            ENDIF
          ENDDO
        ENDIF
 10   CONTINUE
C
C     Memory deallocation
      IF(NZXOPE.GT.0) THEN
        deallocate(GPOE)
        deallocate(XPOE)
      ENDIF
      IF(NZXOPP.GT.0) THEN
        deallocate(GPOP)
        deallocate(XPOP)         
      ENDIF
      deallocate(PRPINT)
      deallocate(LABINT)
      deallocate(INTREP)
      deallocate(first)
      deallocate(doatom)
C     
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck bedleb */                                                   C
      SUBROUTINE BEDLEB(OMEGA,ISSYM,CMO,IBEIG,EVECR,IBTYP,
     &                   IBCVC,IBEVC,IBPVC,WORK,KFREE,LFREE)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C***********************************************************************
      use orbital_rotation_indices
C     
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0,DP5=0.5D0)
#include "dcbgen.h"
#include "chrxyz.h"
#include "dcbbas.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbxrs.h"
#include "dcbxpp.h"
#include "orgcom.h"
#include "dcbxpr.h"
#include "dcbprl.h"            
#include "dcbibn.h"
#include "dcbwav.h"
C#include "cbiher.h"
      LOGICAL DOINT(2,2)
      CHARACTER*6 MXFORM, FMT
      INTEGER LWORK, NCOMP, NOPTP,IBEIG(*)
      DIMENSION OMEGA(MAXEXC),ISSYM(MAXEXC),CMO(*),WORK(*)
      real(8), allocatable              :: PRPINT(:,:),BBUF(:)
      real(8), allocatable              :: TMOM(:,:,:,:),OSC(:,:)
      character (len=8),  allocatable   :: LABINT(:)
      character (len=5)                 :: TPAR(3)
      integer, allocatable              :: INTREP(:),IPSYM2INT(:,:)
      integer, allocatable              :: IVECS(:)
      real(8)             :: PBUF(3,2),EEWT(6),EVECR(NREDM,*)
      integer             :: IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*)
      real(8), allocatable              :: GPOE(:,:),GPOP(:,:)
      real(8), allocatable              :: XPOE(:,:),XPOP(:,:)
      real(8), allocatable              :: X(:),Y(:),Z(:),WT(:)      
      logical, allocatable              :: first(:)
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE
C
#include "ibtfun.h"
      NEXCNV = KEXCNV(JSYMOP)
      IF (NEXCNV .EQ. 0) RETURN
      NPOINTS = NLEB(NROTAV)
C     
C     Allocate placeholders for integral labels and irreps 
C
      allocate(LABINT(8))
      allocate(INTREP(8))
      allocate(first(nz))
C
C     Integral directives
C
      IPRINT=0
      CALL PR1DIR('CXIKR   ',INTTYP,NPRPINT,INTREP,
     &            ANTI,SQUARE,DUMMY,LABINT,TRIANG,TRASPH,
     &            SOLVNT,0,.FALSE.,0,NBAST,NELMNT,IPRINT)
C     NPRPINT is by construction an even number,
C     the first and second halfs corresponding to the real (eq. 44) and imaginary parts (eq. 45)
      IF(NPRPLBL+NPRPINT.GT.MAXPRPLBL) THEN
        CALL QUIT('GET_BED:Redimension MAXPRPLBL !')
      ENDIF
      DO ICOMP = 1,NPRPINT
        PRPLBL (NPRPLBL+ICOMP) = LABINT(ICOMP)
        IPRLREP(NPRPLBL+ICOMP) = INTREP(ICOMP)
        PDOINT (NPRPLBL+ICOMP) = '0++0'
        IPRLTYP(NPRPLBL+ICOMP) =  1        
      ENDDO
      DOINT(1,1)=.FALSE.
      DOINT(2,1)=.TRUE.
      DOINT(1,2)=.TRUE.
      DOINT(2,2)=.FALSE.
C.....allocate array for holding integrals (on lower triangular form)
      allocate(PRPINT(NELMNT,NPRPINT))
      IORDER=0
      NPQUAD=0
C      
C     Find what symmetries are available
C
      allocate(IPSYM2INT(NBSYM,2))
      NCOMP=NPRPINT/2
      CALL ICOPY(2*NBSYM,-1,0,IPSYM2INT,1)
      DO I = 1,NCOMP
        ISYM = INTREP(I) + 1
        IPSYM2INT(ISYM,1) = I
        II = NCOMP + I
        ISYM = INTREP(II) + 1
        IPSYM2INT(ISYM,2) = II
      ENDDO
C
C     Operator setup
C
      NOP    = 3
      IREP   = JSYMOP-1
      ITIM   = -JTIMOP ! Time antisymmetric operator
      IF((NPRPS+2*NOP).GT.MAXPRPS) THEN
        WRITE(LUPRI,'(A,I5)') 
     &        ' ** ERROR in BEDLEB **  Redimension MAXPRPS to ',
     &        NPRPS+2*NOP
        CALL QUIT('BEDLEB: Too may properties !')
      ENDIF
      II = 0
      DO IC = 1,2
        DO IOP = 1,NOP
          II = II + 1
          IOPTYP = 1 + IOP
          IAIND  = JM4(1,IOPTYP)
          IAREP  = JM4REP(IAIND)
          IBSYM  = IBTXOR(IAREP,IREP)+1
          WRITE(PRPNAM(NPRPS+II),'(A3,A1,I1,A7)')
     &        'BED',CHRXYZ(IC),IC,REP(IREP)
          IPRPSYM(NPRPS+II) = JSYMOP ! Total symmetry of operator
          IPRPTIM(NPRPS+II) = ITIM   
          IF(IPSYM2INT(IBSYM,IC).EQ.-1) THEN
            PBUF(IOP,IC) = D0
            IPRPTYP(NPRPS+II) = -1 ! Means that there is no contribution
          ELSE
            PBUF(IOP,IC)         = D1 ! factor 0.5 not included
            IPRPTYP(NPRPS+II)    = IOPTYP
            IPRPLBL(1,NPRPS+II)  = IPSYM2INT(IBSYM,IC)+NPRPLBL
          ENDIF
        ENDDO
      ENDDO 
      CALL PRSYMB(LUPRI,'=',80,0)
      WRITE(LUPRI,'(/A,I2,A,A3/)')
     &     ' *** Excitations of boson symmetry',JSYMOP,' : ',REP(IREP)
      WRITE(6,*) 'BED rotational average. '
      WRITE(6,*) '  Lmax    = ',NROTAV
      WRITE(6,*) '  Npoints = ',NPOINTS
C
C     Prepare solution vectors
C     ========================
C
C      - orbital (e-e) part      
C
      NPAR = 0
      IF(NZXOPE.GT.0) THEN
         allocate(XPOE(NZXOPEQ,NEXCNV))
         allocate(BBUF(NZXOPEQ))
         allocate(IVECS(NESIM))
         NPAR = NPAR + 1
         TPAR(NPAR) = ' e-e '
         OPEN(LUBOE,FILE='PAMBOE',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPEQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBENDX,XPOE,EVECR,NEXCNV,IBTYP,
     &        IBEVC,IVECS,BBUF)
         CLOSE(LUBOE,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
C
C     - orbital (e-p) part      
C
      IF(NZXOPP.GT.0) THEN
         allocate(XPOP(NZXOPPQ,NEXCNV))
         allocate(BBUF(NZXOPPQ))
         allocate(IVECS(NPSIM))         
         NPAR = NPAR + 1
         TPAR(NPAR) = ' e-p '
         OPEN(LUBOP,FILE='PAMBOP',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPPQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBPNDX,XPOP,EVECR,NEXCNV,IBTYP,
     &        IBPVC,IVECS,BBUF)
         CLOSE(LUBOP,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
      TPAR(NPAR+1) = 'Total'
C     
C     Loop over quadrature points
C     =============================
C
      allocate(TMOM(NOP,2,NPAR,NEXCNV))
      allocate(OSC(NPAR+1,NEXCNV))
      CALL DCOPY((NPAR+1)*NEXCNV,D0,0,OSC,1)
      IF(NZXOPE.GT.0) allocate(GPOE(NZXOPEQ,NOP))
      IF(NZXOPP.GT.0) allocate(GPOP(NZXOPPQ,NOP))      
      allocate( X(NPOINTS))
      allocate( Y(NPOINTS))
      allocate( Z(NPOINTS))
      allocate(WT(NPOINTS))
      CALL LEBEDEV(NROTAV,X,Y,Z,WT,NCOUNT)
      IF(NCOUNT.NE.NPOINTS) CALL QUIT('BEDLEB: Grid inconsistency.')
      DO IP = 1,NPOINTS
C.......Set up chi-average over polarization directions        
        EEWT(1) =  DP5*(Y(IP)*Y(IP)+Z(IP)*Z(IP)) ! (1,1)
        EEWT(2) = -X(IP)*Y(IP)                   ! (2,1) ! factor 0.5 cancelled by triangular sum
        EEWT(3) = -X(IP)*Z(IP)                   ! (3,1) ! factor 0.5 cancelled by triangular sum
        EEWT(4) =  DP5*(X(IP)*X(IP)+Z(IP)*Z(IP)) ! (2,2)
        EEWT(5) = -Y(IP)*Z(IP)                   ! (3,2) ! factor 0.5 cancelled by triangular sum
        EEWT(6) =  DP5*(X(IP)*X(IP)+Y(IP)*Y(IP)) ! (3,3)
        DO 10 IEXC = 1,NEXCNV
          IF(ABS(OMEGA(IEXC)).LT.1.0D-3) THEN
             WRITE(LUPRI,*)
     &       'BEDLEB: Skipping excitation. Frequency too small: ',
     &       OMEGA(IEXC)
             GOTO 10
          ENDIF
C
C         Length of wave vector
C
          FAC=OMEGA(IEXC)/CVAL
          II = 0
          DO IC = 1,2
            DO IOP = 1,NOP
              II = II + 1
              WRITE(PRPNAM(NPRPS+II)(9:16),'(A4,I4)') '_exc',IEXC
              FACPRP(1,NPRPS+II) = PBUF(IOP,IC)/FAC
            ENDDO
          ENDDO 
C     
C         Generate integrals
C
          IF(BEDCHK) THEN
            CALL DCOPY(3,D0,0,WAVEVEC,1)
          ELSE
            WAVEVEC(1) = FAC*X(IP)
            WAVEVEC(2) = FAC*Y(IP)
            WAVEVEC(3) = FAC*Z(IP)
          ENDIF
          CALL PR1IN2(WORK(KFREE),LFREE,PRPINT,INTTYP,
     &       NPRPINT,INTREP,ANTI,SQUARE,DUMMY,LABINT,TRIANG,
     &       TRASPH,.FALSE.,SOLVNT,0,0,DUMMY,
     &       0,NBAST,NELMNT,IPRINT,DOINT,.FALSE.)
          IF(IPRXPP.GE.5) THEN
            DO ICOMP = 1,NPRPINT
              WRITE(6,*) '* Integrals of ',LABINT(ICOMP),
     &          '. Symmetry: ',REP(INTREP(ICOMP))
              CALL OUTPAK(PRPINT(1,ICOMP),NBAST,1,LUPRI)
            ENDDO
          ENDIF
C
       IPAR = 0
C
C     Gradient - orbital (e-e) part
C     ====================================
       IF(NZXOPE.GT.0) THEN
          IPAR = IPAR + 1
          II = 0
          TMP = D0
          DO IC = 1,2
            DO IOP = 1,NOP
              II = II + 1
              INDXPR = NPRPS+II
              IF(IPRPTYP(INDXPR).NE.-1) THEN
                CALL GPGET(INDXPR,GPOE(1,IOP),
     &              get_orbital_rotation_indices_pp(),
     &              NZXOPE,CMO,IBEIG,JBENDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
                TMOM(IOP,IC,IPAR,IEXC)=D2*
     &              DDOT(NZXOPEQ,GPOE(1,IOP),1,XPOE(1,IEXC),1)
              ELSE
                TMOM(IOP,IC,IPAR,IEXC) = D0
              ENDIF
            ENDDO
            TMP = TMP
     &     + EEWT(1)*TMOM(1,IC,IPAR,IEXC)*TMOM(1,IC,IPAR,IEXC)
     &     + EEWT(2)*TMOM(2,IC,IPAR,IEXC)*TMOM(1,IC,IPAR,IEXC)
     &     + EEWT(3)*TMOM(3,IC,IPAR,IEXC)*TMOM(1,IC,IPAR,IEXC)
     &     + EEWT(4)*TMOM(2,IC,IPAR,IEXC)*TMOM(2,IC,IPAR,IEXC)
     &     + EEWT(5)*TMOM(3,IC,IPAR,IEXC)*TMOM(2,IC,IPAR,IEXC)
     &     + EEWT(6)*TMOM(3,IC,IPAR,IEXC)*TMOM(3,IC,IPAR,IEXC)
          ENDDO
          OSC(IPAR,IEXC) = OSC(IPAR,IEXC) + WT(IP)*TMP
       ENDIF
C
C     Gradient - orbital (e-p) part
C     ====================================
       IF(NZXOPP.GT.0) THEN
          IPAR = IPAR + 1
          II = 0
          TMP = D0
          DO IC = 1,2
            DO IOP = 1,NOP
              II = II + 1
              INDXPR = NPRPS+II
              IF(IPRPTYP(INDXPR).NE.-1) THEN
                CALL GPGET(INDXPR,GPOP(1,IOP),
     &              get_orbital_rotation_indices_pn(),
     &              NZXOPP,CMO,IBEIG,JBPNDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
                TMOM(IOP,IC,IPAR,IEXC)=D2*
     &              DDOT(NZXOPPQ,GPOP(1,IOP),1,XPOP(1,IEXC),1)
              ELSE
                TMOM(IOP,IC,IPAR,IEXC) = D0
              ENDIF
            ENDDO
            TMP = TMP
     &     + EEWT(1)*TMOM(1,IC,IPAR,IEXC)*TMOM(1,IC,IPAR,IEXC)
     &     + EEWT(2)*TMOM(2,IC,IPAR,IEXC)*TMOM(1,IC,IPAR,IEXC)
     &     + EEWT(3)*TMOM(3,IC,IPAR,IEXC)*TMOM(1,IC,IPAR,IEXC)
     &     + EEWT(4)*TMOM(2,IC,IPAR,IEXC)*TMOM(2,IC,IPAR,IEXC)
     &     + EEWT(5)*TMOM(3,IC,IPAR,IEXC)*TMOM(2,IC,IPAR,IEXC)
     &     + EEWT(6)*TMOM(3,IC,IPAR,IEXC)*TMOM(3,IC,IPAR,IEXC)
          ENDDO
          OSC(IPAR,IEXC) = OSC(IPAR,IEXC) + WT(IP)*TMP
        ENDIF
   10   CONTINUE
      ENDDO
C
C     Final values for oscillator strengths
C
      DO IEXC = 1,NEXCNV
         OSC(NPAR+1,IEXC) = D0
         DO IPAR = 1,NPAR
           OSC(IPAR  ,IEXC) = D2*OMEGA(IEXC)*OSC(IPAR,IEXC)
           OSC(NPAR+1,IEXC) = OSC(NPAR+1,IEXC) + OSC(IPAR,IEXC)
        ENDDO
      ENDDO
C
C     Print section
C
      WRITE(6,*) '* Isotropic oscillator strengths:'
      WRITE(6,'(A4,1X,A18,3(3X,A5,9X))')     
     &  'Exc.','Frequency (au)    ',(TPAR(IPAR),IPAR=1,NPAR+1)
      DO IEXC = 1,NEXCNV
         FMT=MXFORM(OSC(NPAR+1,IEXC),14)
         WRITE(6,'(I4,1X,F18.10,3(3X,'//FMT//'))')
     &     IEXC,OMEGA(IEXC),
     &     (OSC(IPAR,IEXC),IPAR = 1,NPAR+1)
      ENDDO
C            
C     Memory deallocation
      IF(NZXOPE.GT.0) THEN
        deallocate(GPOE)
        deallocate(XPOE)
      ENDIF
      IF(NZXOPP.GT.0) THEN
        deallocate(GPOP)
        deallocate(XPOP)         
      ENDIF
      deallocate(TMOM)
      deallocate(OSC)
      deallocate(PRPINT)
      deallocate(LABINT)
      deallocate(INTREP)
      deallocate(IPSYM2INT)
      deallocate(first)
C     
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCC
