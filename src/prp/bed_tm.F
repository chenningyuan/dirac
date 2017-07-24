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
      PARAMETER(D0=0.0D0,D2=2.0D0,DP5=0.5D0)
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
      INTEGER LWORK, NCOMP, NOPTP,IBEIG(*)
      DIMENSION OMEGA(MAXEXC),ISSYM(MAXEXC),CMO(*),WORK(*)
      real(8), allocatable              :: PRPINT(:,:),BBUF(:)
      real(8), allocatable              :: TMOM(:,:,:)
      character (len=8),  allocatable   :: LABINT(:)
      integer, allocatable              :: INTREP(:),IPSYM2INT(:,:)
      integer, allocatable              :: IVECS(:)
      real(8)             :: PBUF(3,2),EVECR(NREDM,*)
      integer             :: IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*)
      real(8), allocatable              :: GPOE(:,:),GPOP(:,:)
      real(8), allocatable              :: XPOE(:,:),XPOP(:,:)
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
            PBUF(IMAT,IOP)          = UPOL(IMAT) ! facot 0.5 not included
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
         allocate(XPOE(NZXOPPQ,NEXCNV))
         allocate(BBUF(NZXOPPQ))
         allocate(IVECS(NPSIM))         
         NPAR = NPAR + 1
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
      allocate(TMOM(NOP,NVAR,NEXCNV))
      IF(NZXOPE.GT.0) allocate(GPOE(NZXOPEQ,NOP))
      IF(NZXOPP.GT.0) allocate(GPOP(NZXOPPQ,NOP))      
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
C     
C       Generate integrals
C
        DO I = 1,3
          WAVEVEC(I) = FAC*UWAVE(I)
        ENDDO
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
       IVAR = 0
C
C     Gradient - orbital (e-e) part
C     ====================================
       IF(NZXOPE.GT.0) THEN
          IVAR = IVAR + 1
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
               CALL GPGET(INDXPR,GPOE(1,IOP),
     &              get_orbital_rotation_indices_pp(),
     &              NZXOPE,CMO,IBEIG,JBENDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
               TMOM(IOP,IVAR,IEXC)=D2*
     &              DDOT(NZXOPEQ,GPOE(1,IOP),1,XPOE(1,IEXC),1)
            ELSE
               TMOM(IOP,IVAR,IEXC)=D0
            ENDIF
            WRITE(6,*) 'e-e TMOM:',PRPNAM(NPRPS+IOP),OMEGA(IEXC),
     &              TMOM(IOP,IVAR,IEXC)
           ENDDO
        ENDIF
        IF(NZXOPP.GT.0) THEN       
          IVAR = IVAR + 1
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
               CALL GPGET(INDXPR,GPOP(1,IOP),
     &              get_orbital_rotation_indices_pn(),
     &              NZXOPP,CMO,IBEIG,JBPNDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
               TMOM(IOP,IVAR,IEXC)=D2*
     &              DDOT(NZXOPPQ,GPOP(1,IOP),1,XPOP(1,IEXC),1)
               WRITE(6,*) 'e-p TMOM:',PRPNAM(NPRPS+IOP),OMEGA(IEXC),
     &              TMOM(IOP,IVAR,IEXC)
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
      deallocate(TMOM)
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
         allocate(XPOE(NZXOPPQ,NEXCNV))
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
               WRITE(6,*) 'e-e TMOM:',PRPNAM(NPRPS+IOP),OMEGA(IEXC),TMP
            ENDIF
          ENDDO
        ENDIF
        IF(NZXOPP.GT.0) THEN       
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
               CALL GPGET(INDXPR,GPOP(1,IOP),
     &              get_orbital_rotation_indices_pn(),
     &              NZXOPP,CMO,IBEIG,JBPNDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
               TMP=D2*DDOT(NZXOPPQ,GPOP(1,IOP),1,XPOP(1,IEXC),1)
               WRITE(6,*) 'e-p TMOM:',PRPNAM(NPRPS+IOP),TMP
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
      NEXCNV = KEXCNV(JSYMOP)
      IF (NEXCNV .EQ. 0) RETURN
      IF((NPRPS+NOP).GT.MAXPRPS) THEN
        WRITE(LUPRI,'(A,I5)') 
     &        ' ** ERROR in GET_ED_LENGTH **  Redimension MAXPRPS to ',
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
         allocate(XPOE(NZXOPPQ,NEXCNV))
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
               WRITE(6,*) 'e-e TMOM:',PRPNAM(NPRPS+IOP),OMEGA(IEXC),TMP
            ENDIF
          ENDDO
        ENDIF
        IF(NZXOPP.GT.0) THEN       
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
               CALL GPGET(INDXPR,GPOP(1,IOP),
     &              get_orbital_rotation_indices_pn(),
     &              NZXOPP,CMO,IBEIG,JBPNDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
               TMP=D2*DDOT(NZXOPPQ,GPOP(1,IOP),1,XPOP(1,IEXC),1)
               WRITE(6,*) 'e-p TMOM:',PRPNAM(NPRPS+IOP),TMP
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
C  /* Deck bedleb */                                                     C
      SUBROUTINE BEDLEB(OMEGA,ISSYM,IPRINT,WORK, KFREE,LFREE)
C      SUBROUTINE BEDLEB(X,Y,Z,WT,NSCHEME,NCOUNT,NATOM,RNODE,RBRAGG)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h"
#include "priunit.h"
#include "pi.h"
C
#include "dcbgen.h"
#include "dcbbas.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbxpp.h"
#include "orgcom.h"
C
      PARAMETER (D0=0.0D0,DP5=0.5D0,D2=2.0D0)
      PARAMETER (MAXNUM=1455)
      DIMENSION X(MAXNUM),Y(MAXNUM),Z(MAXNUM),WT(MAXNUM)
      DIMENSION OMEGA(MAXEXC,NBSYM),ISSYM(MAXEXC,NBSYM),WORK(*)
      REAL(8):: KVEC(3,MAXNUM),E1VEC(3),E2VEC(3),FAC,TMPNRM
      REAL(8):: E1E1(3,3),E1E2(3,3),OSC(MAXEXC),ECD(MAXEXC)

C      IF LEBEDEV THEN
        CALL LD0014(X,Y,Z,WT,NCOUNT)
C      END IF
      FAC=OMEGA(1,1)/CVAL
C
      DO IK=1,NCOUNT
C   
C     Creating wave vector for the given Lebedev point or as specified
C     in the input      
C
        KVEC(1,IK)=FAC * X(IK) 
        KVEC(2,IK)=FAC * Y(IK) 
        KVEC(3,IK)=FAC * Z(IK) 
CC       
C       CALL INTEGRAL ROUTINE
C
C.......Weights for chi-average over pol-directions
C       Linear absorption 
C       <E1VEC,E1VEC>_chi average
        E1E1(1,1)=DP5*(Y(IK)*Y(IK)+Z(IK)*Z(IK))
        E1E1(1,2)=-DP5*X(IK)*Y(IK)
        E1E1(1,3)=-DP5*X(IK)*Z(IK)
        E1E1(2,1)=E1E1(1,2)
        E1E1(2,2)=DP5*(X(IK)*X(IK)+Z(IK)*Z(IK))
        E1E1(2,3)=-DP5*Y(IK)*Z(IK)
        E1E1(3,1)=E1E1(1,3)
        E1E1(3,2)=E1E1(2,3)
        E1E1(3,3)=DP5*(X(IK)*X(IK)+Y(IK)*Y(IK))
C
C     AVERAGE OVER POL-DIRECTION      
      WRITE(LUPRI,*) 'VEC', X(IK), Y(IK), Z(IK), WT(IK)
      DO I=1,3
        WRITE(LUPRI,*)  E1E1(I,1), E1E1(I,2), E1E1(I,3)
      END DO
C
C      DO I=1,3
C        DO J=I,3
C          IF I.EQ.J THEN
C            OSC(IEXC)=OSC(IEXC)+WT(IK)*E1E1(I,I)
C     &                         *(TMR(I)*TMR(I)+TMI(I)*TMI(I))
C          ELSE
C            OSC(IEXC)=OSC(IEXC)+D2*WT(IK)*E1E1(I,J)
C     &                         *(TMR(I)*TMR(J)-TMI(I),TMI(J))
C            ECD(IEXC)=ECD(IEXC)+WT(IK)*E1E2(I,J)
C     &                         *(TMI(I)*TMR(J)+TMR(I)*TMI(J))
C          END IF 
C        END DO 
C      END DO
C     
      END DO
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCC
