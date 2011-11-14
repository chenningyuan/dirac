!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

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
      use codata
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(AMUKG = 1.66054020D-27,ANGM = 1.0D-10)
      PARAMETER(XTM = XTANG*ANGM,XTCM=2.19474625D+05)
      PARAMETER(DM1 = -1.0D0,D1 = 1.0D0,D2 = 2.0D0)
      PARAMETER(D2PI = D2*PI,D4=4.0D0,CENT=1.0D2)
      PARAMETER(NITER=10)
      PARAMETER(MORDER=16,MP=50,MXPRP=5)
      PARAMETER(MDIM=MORDER+1)
      DIMENSION A(MP,MDIM),B(MP),C(MDIM),D(MDIM)
      DIMENSION X(MP),Y(MP),P(MP,MXPRP)
      CHARACTER POTFIL*40,REPLY*1
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
     &    'Name of input file with potential curve with '//
     &    '" R E(R)" values (A40)'
        READ(LUCMD,'(A)') POTFIL
        NPRP = 0
      ELSEIF(I.EQ.2) THEN
        WRITE(LUPRI,'(A)') 
     &    'Name of input file with potential/property curve (A40)'
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
     &  'Select one of the following:',
     &  '  1. Bond lengths in Angstroms.',
     &  '  2. Bond lengths in atomic units.',
     &  ' (Note that all other quantities only in atomic units !)'
      READ(LUCMD,*) IANG
      i=LNBLNK(POTFIL)
      OPEN(IUNIT,FILE=POTFIL,STATUS='OLD',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL',ERR=20)
      OPEN(ILOGG,FILE=POTFIL(1:i)//'.twofit',STATUS='UNKNOWN',
     &     FORM='FORMATTED',ACCESS='SEQUENTIAL')
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
      WRITE(ILOGG,'(A,A40)') '* Points read from ',POTFIL
      WRITE(ILOGG,'(3X,F6.3,4X,F20.10)') (X(I),Y(I),I=1,NP)
C     Convert x-values to atomic units
      IF (IANG .EQ. 1) THEN
        DO I = 1,NP
         X(I) = X(I)/XTANG
        ENDDO
      ELSE IF (IANG .NE. 2) THEN
        STOP 'Unknown bond length unit option!'
      END IF
C
C     Polynomial fit
C     ==============
C
 60   CONTINUE
      WRITE(LUPRI,'(/A/)') '** SPECTROSCOPIC CONSTANTS:'
      WRITE(ILOGG,'(/A/)') '** SPECTROSCOPIC CONSTANTS:'
      CALL MODEFIT(NP,X,Y,A,B,C,D,MORDER,NORDER,ILOGG)
C
C     Find minimum
C     ============
C
      CALL MODEMIN(XOPT,NP,X,Y,C,NORDER,NITER,ILOGG)
C
C     Get masses
C     ==========
C
 10   CONTINUE
      WRITE(LUPRI,'(A)')
     &  'Select one of the following:',
     &  '  1. Select masses of the most abundant isotopes.',
     &  '  2. Employ user-defined atomic masses.',
     &  '  3. Normal modes: Give reduced mass.'
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
      ELSEIF(I.EQ.3) THEN
        WRITE(LUPRI,'(A)') '* Give reduced mass in Daltons:'
        READ(LUCMD,*) UM
        GOTO 50
      ELSE
        WRITE(LUPRI,'(A)') ' You stupid fart ! Try again !'
        GOTO 10
      ENDIF
      UM = D1/(D1/AM + D1/BM)
      WRITE(ILOGG,'(A)') '* MASSES:'
      WRITE(ILOGG,'(3X,A,F8.4)') 'Atom A      :',AM
      WRITE(ILOGG,'(3X,A,F8.4)') 'Atom B      :',BM
 50   CONTINUE
      WRITE(ILOGG,'(3X,A,F8.4)') 'Reduced mass:',UM
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      UM = XFAMU*UM
C
C     Harmonic constant
C     ==================
C
      CALL MODEHARM(OMEGA,V2,XOPT,UM,C,NP,NORDER,ILOGG)
      CALL MODEDISP(X,NP,XOPT,UM,OMEGA,ILOGG)
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      IF(NORDER.LT.4) GOTO 666 
C
C     Anharmonic constant
C     ===================
C
      CALL MODEANHARM(WXE,OMEGA,XOPT,UM,C,NP,NORDER,ILOGG)
      CALL MODEDISS(OMEGA,WXE,ILOGG)
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
C
      WRITE(LUPRI,*) 'Do you want to calculate spectroscopic constants',
     &   ' with other isotope(s) (y/n)?'
      READ(LUCMD,*) REPLY
      IF(REPLY.EQ.'y') GO TO 10
C
      WRITE(LUPRI,*) 'Do you want to calculate spectroscopic constants',
     &   ' with another polynomial order (y/n)?'
      READ(LUCMD,*) REPLY
      IF(REPLY.EQ.'y') GO TO 60
C
C     Anharmonic constant at non-stationary geometry
C     ==============================================
C
      WRITE(LUPRI,*) 'Do you want to calculate spectroscopic constants',
     &   ' at another geometry (y/n)?'
      READ(LUCMD,*) REPLY
      IF(REPLY.EQ.'y') THEN
        WRITE(LUPRI,*) 'Give geometry in Angstroms '
        READ(5,*) XNEW
        WRITE(ILOGG,'(A,F9.3,A)') 
     &    'Spectroscopic at non-stationary point ',
     &                  XNEW, 'Angstrom '
        XNEW = XNEW/XTANG
        CALL MODEHARM(OMEGAN,V2N,XNEW,UM,C,NP,NORDER,ILOGG)
        CALL MODEANHARM(WXEN,OMEGAN,XNEW,UM,C,NP,NORDER,ILOGG)
      ENDIF
      IF(NPRP.EQ.0) GOTO 666
C      
C     Effective geometry (Newton-Raphson search)
C     ==========================================
C
      WRITE(LUPRI,'(/A/)') '** PROPERTIES **'
      WRITE(ILOGG,'(/A/)') '** PROPERTIES **'
      XEFF = XOPT - (V3/(4.0D0*OMEGA*(OMEGA*UM)*(OMEGA*UM)))
      TEMP = XEFF*XTANG
      WRITE(LUPRI,'(A,F18.4,A)') 
     &    '* Effective geometry  :',TEMP,' Angstroms'
      WRITE(ILOGG,'(A,F18.4,A)') 
     &    '* Effective geometry  :',TEMP,' Angstroms'
      FORCE = POL2DER(NORDER,C,XEFF)
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
        CALL POLSVD(NDIM,NP,A,B,X,P(1,IPRP),C,D,CHISQ,ISKIP)
C* Estimate fit:
        WRITE(*,'(/A,E9.4/)') '* Chi square :  ',CHISQ
        WRITE(LUPRI,'(A,I3)') '* Polynomial fit of order:',NORDER
        WRITE(LUPRI,'(A)') '* Coefficients:'
        WRITE(LUPRI,'(A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',
     &        C(I),I=1,NDIM)  
        WRITE(LUPRI,'(A,E9.4)') '* Chi square :  ',CHISQ
        WRITE(LUPRI,'(3X,A,I3)') 
     &    '* Number of singularities(SVD): ',ISKIP
C
        WRITE(ILOGG,'(72A1)') ('=',I=1,72)
        WRITE(ILOGG,'(A,I3)') '* Polynomial fit of order:',NORDER
        WRITE(ILOGG,'(A)') '* Coefficients:'
        WRITE(ILOGG,'(A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',
     &        C(I),I=1,NDIM)  
        WRITE(ILOGG,'(A,E9.4)') '* Chi square :  ',CHISQ
        WRITE(ILOGG,'(3X,A,I3)') 
     &    '* Number of singularities(SVD): ',ISKIP
C
        PMIN = POLVAL(NORDER,C,XOPT)
        WRITE(LUPRI,'(A,1P,E18.10)') 
     &    '* Expectation value at local minimum     : ', PMIN
        WRITE(ILOGG,'(A,1P,E18.10)') 
     &    '* Expectation value at local minimum     : ', PMIN
        PEFF = POLVAL(NORDER,C,XEFF)
        WRITE(LUPRI,'(A,1P,E18.10)') 
     &    '* Expectation value at effective geometry: ', PEFF
        WRITE(ILOGG,'(A,1P,E18.10)') 
     &    '* Expectation value at effective geometry: ', PEFF
        P2EFF = POL2DER(NORDER,C,XEFF)
        PZPV = PEFF + (P2EFF/D4/UM/OMEFF)
        WRITE(LUPRI,'(A,1P,E18.10)') 
     &    '* Zero-point vibrational average         : ', PZPV
        WRITE(ILOGG,'(A,1P,E18.10)') 
     &    '* Zero-point vibrational average         : ', PZPV
        PVIB = CENT*(PZPV-PMIN)/PZPV
        WRITE(LUPRI,'(A,F6.2,A)') 
     &   '* Vibrational effect                      : ',PVIB, ' \\%'
        WRITE(ILOGG,'(A,F6.2,A)') 
     &   '* Vibrational effect                      : ',PVIB, ' \\%'
      ENDDO
 666  CONTINUE
      END
