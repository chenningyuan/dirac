      PROGRAM TWOFIT
C***********************************************************************
C
C     Utility program for diatomics
C     This program will do a quartic fit of potential curve,
C     then calculate spectroscopic constants:
C       - equilibrium distance
C       - harmonic frequency
C       - anharmonic constant
C       - mean displacement in harmonic ground state
C     With property curves also given the program will
C     calculate vibrational effects on the corresponding
C     expectation value
C
C     Written by T. Saue April 2001
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <codata.h>
#include <pi.h>
      PARAMETER(AMUKG = 1.66054020D-27,ANGM = 1.0D-10)
      PARAMETER(XTM = XTANG*ANGM,XTCM=2.19474625D+05)
      PARAMETER(DTOL = 1.0D-4,DM1 = -1.0D0,D1 = 1.0D0,D2 = 2.0D0)
      PARAMETER(D2PI = D2*PI,D4=4.0D0)
      PARAMETER(NITER=10)
      PARAMETER(MORDER=10,MP=50,MXPRP=5)
      PARAMETER(MDIM=MORDER+1,MNDIM=(MDIM*(MDIM+1)/2))
      DIMENSION A(MNDIM),B(MDIM),IBUF(MDIM)
      DIMENSION X(MP),Y(MP),P(MP,MXPRP)
      CHARACTER POTFIL*6,RESFIL*10,REPLY*1
C
C     Heading
C
      CALL TITLER('TWOFIT for diatomics : Written by T. Saue ','*',110)
      WRITE(LUPRI,'(A)')
     &  'Select one of the following:',
     &  '  1. Spectroscopic constants.',
     &  '  2. Spectroscopic constants + properties.'
      READ(LUCMD,*) I
      IF(I.EQ.1) THEN
        WRITE(LUPRI,'(A)') 
     &    'Name input file with potential curve (A6)'
        NPRP = 0
      ELSEIF(I.EQ.2) THEN
        WRITE(LUPRI,'(A)') 
     &    'Name input file with potential/property curve (A6)'
        WRITE(LUPRI,'(A)') 
     &    'Give number of property curves'
        READ(LUCMD,*) NPRP
      ELSE
        STOP 'Unknown option !'
      ENDIF
C
C     Read points
C     ===========
C
      IUNIT = 1
      ILOGG = 2
      READ(LUCMD,'(A)') POTFIL
      WRITE(LUPRI,'(A)') 'Note that atomic units should used !'
      OPEN(IUNIT,NAME=POTFIL,STATUS='OLD',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL',ERR=20)
      RESFIL = POTFIL//'.res'
      OPEN(ILOGG,NAME=RESFIL,STATUS='UNKNOWN',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL')
      GOTO 30
 20   CONTINUE
      STOP 'Error opening POTFIL'
 30   CONTINUE
      NPOINTS = 0
      DO I = 1,MP
        READ(IUNIT,*,END=40) X(I),Y(I),(P(I,J),J=1,NPRP)
        NP = NP + 1
      ENDDO
 40   CONTINUE
CTROND
      DO I = 1,NP
       X(I) = X(I)/XTANG
      ENDDO
      WRITE(LUPRI,'(A,I5)') '* Number of points read:',NP
      IF(NP.LT.5) THEN
        WRITE(LUPRI,'(A)') '* You must have at least five points !'
        GOTO 666
      ENDIF
C
C     Get masses
C     ==========
C
 10   CONTINUE
      WRITE(LUPRI,'(A)')
     &  'Select one of the following:',
     &  '  1. Select masses of the most abundant isotopes.',
     &  '  2. Employ user-defined atomic masses.'
      READ(LUCMD,*) I
      IF(I.EQ.1) THEN
        WRITE(LUPRI,'(A)') '* Give charge of atom A:'
        READ(LUCMD,*) IZ
        CALL GETMSS(IZ,AM,ABUND,INFO)
        IF(INFO.EQ.0) THEN
          WRITE(LUPRI,'(A,F12.4)') 
     &    '* Mass     :', AM,
     &    '* Abundance:', ABUND
        ELSE
          WRITE(LUPRI,'(A)') 'Routine failed. Give mass in Daltons:'
          READ(LUCMD,*) AM
        ENDIF
        WRITE(LUPRI,'(A)') '* Give charge of atom B :'
        READ(LUCMD,*) IZ
        CALL GETMSS(IZ,BM,ABUND,INFO)
        IF(INFO.EQ.0) THEN
          WRITE(LUPRI,'(A,F12.4)') 
     &    '* Mass     :', BM,
     &    '* Abundance:', ABUND
        ELSE
          WRITE(LUPRI,'(A)') 'Routine failed. Give mass in Daltons:'
          READ(LUCMD,*) BM
        ENDIF
      ELSEIF(I.EQ.2) THEN
        WRITE(LUPRI,'(A)') '* Give mass of atom A in Daltons:'
        READ(LUCMD,*) AM
        WRITE(LUPRI,'(A)') '* Give mass of atom B in Daltons:'
        READ(LUCMD,*) BM
      ELSE
        WRITE(LUPRI,'(A)') ' You stupid fart ! Try again !'
        GOTO 10
      ENDIF
      UM = D1/(D1/AM + D1/BM)
      WRITE(ILOGG,'(A)') '* MASSES:'
      WRITE(ILOGG,'(3X,A,F8.4)') 'Atom A      :',AM
      WRITE(ILOGG,'(3X,A,F8.4)') 'Atom B      :',BM
      WRITE(ILOGG,'(3X,A,F8.4)') 'Reduced mass:',UM
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      UM = XFAMU*UM
C
C     Polynomial fit
C     ==============
C
      WRITE(ILOGG,'(A,A6)') '* Points read from ',POTFIL
      WRITE(ILOGG,'(F6.3,4X,F15.6)') (X(I),Y(I),I=1,NP)
 60   CONTINUE
      WRITE(LUPRI,'(A)') 'Polynomial fit: Give order of polynomial'
      READ(LUCMD,*) NORDER
      IF(NORDER.LT.4) THEN
        WRITE(LUPRI,'(A)') 
     &    '* This program requires at least a quartic fit.'
        GOTO 666
      ENDIF
      IF(NP.LE.NORDER) THEN
        WRITE(LUPRI,'(A,I3)') 
     +    'Too few points for a polynomial fit of order',NORDER
        GOTO 666
      ENDIF
      IF(NORDER.GT.MORDER) THEN
        WRITE(LUPRI,'(A)') 'Order to large ....'
        GOTO 666
      ENDIF
      NDIM = NORDER + 1
      CALL CURVE(NORDER,NP,A,B,X,Y,VAR,IBUF)
C* Estimate fit:
      WRITE(*,'(/A,E9.4/)') '* Fit :  ',VAR      
      WRITE(LUPRI,'(A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(LUPRI,'(A)') '* Coefficients:'
      WRITE(LUPRI,'(A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',B(I),I=1,NDIM)  
      WRITE(LUPRI,'(/A,E9.4/)') '* Fit :  ',VAR      
C
      WRITE(ILOGG,'(72A1)') ('=',I=1,72)
      WRITE(ILOGG,'(A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(ILOGG,'(A)') '* Coefficients:'
      WRITE(ILOGG,'(A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',
     &      B(I),I=1,NDIM)  
      WRITE(ILOGG,'(/A,E9.4/)') '* Fit :  ',VAR      
C
C     Find minimum
C     ============
C
C* Find smallest y-value and take corresponding x-value as start value
C* for Newton-Raphsons method      
      IYMIN = IDMIN(NP,Y,1)
      XMIN = X(IYMIN)
      CALL NEWRAP(NORDER,DTOL,B,XMIN,NITER,IERR)
      IF(IERR.EQ.1) THEN
        WRITE(LUPRI,'(A)') 
     &    'Did not find minimumn: Newtons method failed'
        WRITE(ILOGG,'(A)') 'Newtons method failed'
        GOTO 666
      ELSE
        YMIN = POLVAL(NORDER,B,XMIN)
        WRITE(LUPRI,'(A,F17.10,A,F17.10,A)') '* Local minimum:  (',
     +        XMIN*XTANG,',',YMIN,')'
        WRITE(ILOGG,'(72A1)') ('-',I=1,72)
        WRITE(ILOGG,'(A,F17.10,A,F17.10,A)') '* Local minimum:  (',
     +        XMIN*XTANG,',',YMIN,')'
      ENDIF
C
C     Harmonic constant
C     ==================
C
      FORCE = POL2DER(NORDER,B,XMIN)
      FORCE2 = FORCE*XTJ/XTM/XTM
      OMEGA = SQRT(FORCE/UM)
      FREQ  = OMEGA/D2PI
      FREQ2 = OMEGA*XTHZ
      FREQ  = SQRT(FORCE/UM)/(SQRT(AMUKG)*D2PI)
      WRITE(LUPRI,'(A,1P,E14.4,A)') '* Force constant:',FORCE2, 'N/m'
      WRITE(LUPRI,'(A,1P,E10.3,A,3X,E10.3,A)') '* Frequency:',
     &     OMEGA*XTHZ,'Hz',OMEGA*XTCM,'cm-1'
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      WRITE(ILOGG,'(A,F14.4,A)') '* Force constant:',FORCE2, 'N/m'
      WRITE(ILOGG,'(A,1P,E10.3,A,3X,E10.3,A)') '* Frequency:',
     &     OMEGA*XTHZ,'Hz',OMEGA*XTCM,'cm-1'
C      
C     Effective geometry (Newton-Raphson search)
C     ==========================================
C
      MAXIT = 10
      ACC=1.0D-6
      XX = XMIN
      UMS= SQRT(UM)
      AA = D1/D4/UMS
      DO IT = 1,MAXIT
        V1 = POLNDER(NORDER,B,XX,1)
        V2 = POLNDER(NORDER,B,XX,2)
        V2S= SQRT(V2)
        V3 = POLNDER(NORDER,B,XX,3)
        V4 = POLNDER(NORDER,B,XX,4)
        F0 = V1+ AA*V3/V2S
        F1 = V2 + AA*V4/V2S - AA*V3*V3/V2/V2S/D2
        DX = F0/F1
        XX = XX - DX
        IF(ABS(DX).LT.ACC) GOTO 70
      ENDDO
 70   CONTINUE
      WRITE(LUPRI,'(A,F17.10,A,F17.10,A)') 
     &    '* Effective geometry: ',XX*XTANG
 666  CONTINUE
      END
C
