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
      PARAMETER(D2PI = D2*PI,D4=4.0D0,CENT=1.0D2)
      PARAMETER(NITER=10)
      PARAMETER(MORDER=10,MP=50,MXPRP=5)
      PARAMETER(MDIM=MORDER+1,MNDIM=(MDIM*(MDIM+1)/2))
      DIMENSION A(MNDIM),B(MDIM),BP(MDIM),IBUF(MDIM)
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
C
C     Read points
C     ===========
C
      IUNIT = 1
      ILOGG = 2
      READ(LUCMD,*) I
      IF(I.EQ.1) THEN
        WRITE(LUPRI,'(A)') 
     &    'Name input file with potential curve (A6)'
        READ(LUCMD,'(A)') POTFIL
        NPRP = 0
      ELSEIF(I.EQ.2) THEN
        WRITE(LUPRI,'(A)') 
     &    'Name input file with potential/property curve (A6)'
        READ(LUCMD,'(A)') POTFIL
        WRITE(LUPRI,'(A)') 
     &    'Give number of property curves'
        READ(LUCMD,*) NPRP
        IF(NPRP.GT.MXPRP) THEN
          WRITE(LUPRI,'(A,I5)') 
     &       '* Maximum number of properties  : ',MXPRP,
     &       '* Requested number of properties: ',NPRP
          STOP 'Too many properties !'
        ENDIF
      ELSE
        STOP 'Unknown option !'
      ENDIF
      WRITE(LUPRI,'(A)') 
     &  'Note that bond lengths should be given in Angstroms',
     &  ' and all other quantities in atomic units !'
      OPEN(IUNIT,NAME=POTFIL,STATUS='OLD',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL',ERR=20)
      RESFIL = POTFIL//'.res'
      OPEN(ILOGG,NAME=RESFIL,STATUS='UNKNOWN',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL')
      GOTO 30
 20   CONTINUE
      STOP 'Error opening POTFIL'
 30   CONTINUE
      NP = 0
      IF(NPRP.GT.0) THEN
        DO I = 1,MP
          READ(IUNIT,*,END=40) X(I),Y(I),(P(I,J),J=1,NPRP)
          NP = NP + 1
        ENDDO
      ELSE
        DO I = 1,MP
          READ(IUNIT,*,END=40) X(I),Y(I)
          NP = NP + 1
        ENDDO
      ENDIF
 40   CONTINUE
      WRITE(LUPRI,'(A,I5)') '* Number of points read:',NP
      WRITE(ILOGG,'(A,A6)') '* Points read from ',POTFIL
      WRITE(ILOGG,'(3X,F6.3,4X,E18.8)') (X(I),Y(I),I=1,NP)
      IF(NP.LT.5) THEN
        WRITE(LUPRI,'(A)') '* You must have at least five points !'
        STOP 'Too few points for a quartic fit !'
      ENDIF
      DO I = 1,NP
       X(I) = X(I)/XTANG
      ENDDO
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
 60   CONTINUE
      WRITE(LUPRI,'(/A/)') '** SPECTROSCOPIC CONSTANTS:'
      WRITE(ILOGG,'(/A/)') '** SPECTROSCOPIC CONSTANTS:'
      NORDER = 4
CTROND      WRITE(LUPRI,'(A)') 'Polynomial fit: Give order of polynomial'
CTROND      READ(LUCMD,*) NORDER
      IF(NORDER.LT.4) THEN
        WRITE(LUPRI,'(A)') 
     &    '* This program requires at least a quartic fit.'
        STOP 'Too low order fit !'
      ENDIF
      IF(NP.LE.NORDER) THEN
        WRITE(LUPRI,'(A,I3)') 
     +    'Too few points for a polynomial fit of order',NORDER
        STOP 'Too few points !'
      ENDIF
      IF(NORDER.GT.MORDER) THEN
        STOP 'Order to large ....'
      ENDIF
      NDIM = NORDER + 1
      CALL CURVE(NORDER,NP,A,B,X,Y,VAR,IBUF)
C* Estimate fit:
      WRITE(LUPRI,'(3X,A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(LUPRI,'(3X,A)') '* Coefficients:'
      WRITE(LUPRI,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',
     &      B(I),I=1,NDIM)  
      WRITE(LUPRI,'(3X,A,E9.4)') '* Fit :  ',VAR      
C
      WRITE(ILOGG,'(72A1)') ('=',I=1,72)
      WRITE(ILOGG,'(3X,A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(ILOGG,'(3X,A)') '* Coefficients:'
      WRITE(ILOGG,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',
     &      B(I),I=1,NDIM)  
      WRITE(ILOGG,'(3X,A,E9.4)') '* Fit :  ',VAR      
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
        WRITE(ILOGG,'(A)') 'Newtons method failed'
        STOP 'Did not find minimum. Newtons method failed...'
      ELSE
        YMIN = POLVAL(NORDER,B,XMIN)
        TEMP = XMIN*XTANG
        WRITE(LUPRI,'(A,F18.5,A)') 
     &    '* Local minimum   :',TEMP,' Angstroms'
        WRITE(ILOGG,'(A,F18.5,A)') 
     &    '* Local minimum   :',TEMP,' Angstroms'
        WRITE(LUPRI,'(A,1P,E18.10,A)') 
     &    '* Expected energy :',YMIN,' Hartrees'
        WRITE(ILOGG,'(A,1P,E18.10,A)') 
     &    '* Expected energy :',YMIN,' Hartrees'
      ENDIF
C
C     Harmonic constant
C     ==================
C
      FORCE = POL2DER(NORDER,B,XMIN)
      TEMP  = FORCE*XTJ/XTM/XTM
      WRITE(LUPRI,'(A,1P,E18.4,A)') 
     &   '* Force constant  :',TEMP, ' N/m'
      WRITE(ILOGG,'(A,1P,E18.4,A)') 
     &   '* Force constant  :',TEMP, ' N/m'
      OMEGA = SQRT(FORCE/UM)
      WRITE(LUPRI,'(A,1P,E10.3,A)') 
     &   '* Frequency       :',OMEGA*XTHZ,' Hz',
     &   '                   ',OMEGA*XTCM,' cm-1'
      WRITE(ILOGG,'(A,1P,E10.3,A)') 
     &   '* Frequency       :',OMEGA*XTHZ,' Hz',
     &   '                   ',OMEGA*XTCM,' cm-1'
      TEMP = D1/SQRT(D2*UM*OMEGA)
      WRITE(LUPRI,'(A,F10.4,A)') 
     &   '* Mean displacement in harmonic ground state:',TEMP,
     &   ' Angstroms'
      WRITE(ILOGG,'(A,F10.4,A)') 
     &   '* Mean displacement in harmonic ground state:',TEMP,
     &   ' Angstroms'
C
C     Anharmonic constant
C     ===================
C
      V3 = POLNDER(NORDER,B,XMIN,3)
      V4 = POLNDER(NORDER,B,XMIN,4)
      TEMP = D1/UM/OMEGA
      T2   = TEMP*TEMP
      WXE = (5.0D0*TEMP*T2*V3*V3/48.0D0/OMEGA)-(V4*T2/16.0D0)
      WRITE(LUPRI,'(A,1P,E10.3,A)') 
     &   '* Omega*x_e       :',WXE*XTCM,' cm-1'
      WRITE(ILOGG,'(A,1P,E10.3,A)') 
     &   '* Omega*x_e       :',WXE*XTCM,' cm-1'
      DE = OMEGA*OMEGA/WXE/D4
      WRITE(LUPRI,'(A,F8.4,A)') '* Estimated D_e:',DE*XTEV, 'eV'
      WRITE(LUPRI,'(A,F8.4,A)') '* Estimated D_e:',DE*XKJMOL, 'kJ/mol'
      WRITE(ILOGG,'(A,F8.4,A)') '* Estimated D_e:',DE*XTEV, 'eV'
      D0 = DE-(OMEGA/D2)+WXE/D4
      WRITE(LUPRI,'(A,F8.4,A)') '* Estimated dissociation energy:',
     &       D0*XTEV, 'eV'
      WRITE(LUPRI,'(A,F8.4,A)') '* Estimated dissociation energy:',
     &       D0*XKJMOL, 'kJ/mol'
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      IF(NPRP.EQ.0) GOTO 666
C      
C     Effective geometry (Newton-Raphson search)
C     ==========================================
C
      WRITE(LUPRI,'(/A/)') '** PROPERTIES **'
      WRITE(ILOGG,'(/A/)') '** PROPERTIES **'
      MAXIT = 10
      ACC=1.0D-6
      XEFF = XMIN
      UMS= SQRT(UM)
      AA = D1/D4/UMS
      DO IT = 1,MAXIT
        V1 = POLNDER(NORDER,B,XEFF,1)
        V2 = POLNDER(NORDER,B,XEFF,2)
        V2S= SQRT(V2)
        V3 = POLNDER(NORDER,B,XEFF,3)
        V4 = POLNDER(NORDER,B,XEFF,4)
        F0 = V1+ AA*V3/V2S
        F1 = V2 + AA*V4/V2S - AA*V3*V3/V2/V2S/D2
        DX = F0/F1
        XEFF = XEFF - DX
        IF(ABS(DX).LT.ACC) GOTO 70
      ENDDO
      STOP 'Newtons method failed for effective geometry !'
 70   CONTINUE
      TEMP = XEFF*XTANG
      WRITE(LUPRI,'(A,F18.4,A)') 
     &    '* Effective geometry  :',TEMP,' Angstroms'
      WRITE(ILOGG,'(A,F18.4,A)') 
     &    '* Effective geometry  :',TEMP,' Angstroms'
      FORCE = POL2DER(NORDER,B,XEFF)
      OMEFF = SQRT(FORCE/UM)
      WRITE(LUPRI,'(A,1P,E18.3,A)') 
     &   '* Effective frequency  :',OMEFF*XTHZ,' Hz',
     &   '                        ',OMEFF*XTCM,' cm-1'
      WRITE(ILOGG,'(A,1P,E18.3,A)') 
     &   '* Effective frequency  :',OMEFF*XTHZ,' Hz',
     &   '                       ',OMEFF*XTCM,' cm-1'
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
C     Zero-point vibrational averages
      DO IPRP = 1,NPRP
        WRITE(LUPRI,'(/A,I3/)') '** PROPERTY NO. ',IPRP
        WRITE(ILOGG,'(/A,I3/)') '** PROPERTY NO. ',IPRP
C       Do polynomial fit
        CALL CURVE(NORDER,NP,A,BP,X,P(1,IPRP),VAR,IBUF)
C* Estimate fit:
        WRITE(*,'(/A,E9.4/)') '* Fit :  ',VAR      
        WRITE(LUPRI,'(A,I3)') '* Polynomial fit of order:',NORDER
        WRITE(LUPRI,'(A)') '* Coefficients:'
        WRITE(LUPRI,'(A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',
     &        BP(I),I=1,NDIM)  
        WRITE(LUPRI,'(A,E9.4)') '* Fit :  ',VAR      
C
        WRITE(ILOGG,'(72A1)') ('=',I=1,72)
        WRITE(ILOGG,'(A,I3)') '* Polynomial fit of order:',NORDER
        WRITE(ILOGG,'(A)') '* Coefficients:'
        WRITE(ILOGG,'(A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',
     &        BP(I),I=1,NDIM)  
        WRITE(ILOGG,'(A,E9.4)') '* Fit :  ',VAR      
C
        PMIN = POLVAL(NORDER,BP,XMIN)
        WRITE(LUPRI,'(A,1P,E18.10)') 
     &    '* Expectation value at local minimum     : ', PMIN
        WRITE(ILOGG,'(A,1P,E18.10)') 
     &    '* Expectation value at local minimum     : ', PMIN
        PEFF = POLVAL(NORDER,BP,XEFF)
        WRITE(LUPRI,'(A,1P,E18.10)') 
     &    '* Expectation value at effective geometry: ', PEFF
        WRITE(ILOGG,'(A,1P,E18.10)') 
     &    '* Expectation value at effective geometry: ', PEFF
        P2EFF = POL2DER(NORDER,BP,XEFF)
        PZPV = PEFF + (P2EFF/D4/UM/OMEFF)
        WRITE(LUPRI,'(A,1P,E18.10)') 
     &    '* Zero-point vibrational average         : ', PZPV
        WRITE(ILOGG,'(A,1P,E18.10)') 
     &    '* Zero-point vibrational average         : ', PZPV
        PVIB = CENT*(PZPV-PMIN)/PZPV
        WRITE(LUPRI,'(A,F6.2,A)') 
     &   '* Vibrational effect                      : ',PVIB, ' \%'
        WRITE(ILOGG,'(A,F6.2,A)') 
     &   '* Vibrational effect                      : ',PVIB, ' \%'
      ENDDO
 666  CONTINUE
      END
C
