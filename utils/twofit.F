C 
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C 
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
      PARAMETER(AMUKG = 1.66054020D-27,ANGM = 1.0D-10)
      PARAMETER(XTM = XTANG*ANGM,XTCM=2.19474625D+05)
      PARAMETER(DTOL = 1.0D-4,DM1 = -1.0D0,D1 = 1.0D0,D2 = 2.0D0)
      PARAMETER(D2PI = D2*PI,D4=4.0D0,CENT=1.0D2)
      PARAMETER(NITER=10)
      PARAMETER(MORDER=16,MP=50,MXPRP=5)
      PARAMETER(MDIM=MORDER+1)
      DIMENSION A(MP,MDIM),B(MP),C(MDIM),D(MDIM)
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
      OPEN(IUNIT,FILE=POTFIL,STATUS='OLD',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL',ERR=20)
      RESFIL = POTFIL//'.res'
      OPEN(ILOGG,FILE=RESFIL,STATUS='UNKNOWN',FORM='FORMATTED',
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
      WRITE(ILOGG,'(3X,F6.3,4X,F20.10)') (X(I),Y(I),I=1,NP)
      DO I = 1,NP
       X(I) = X(I)/XTANG
      ENDDO
C
C     Polynomial fit
C     ==============
C
 60   CONTINUE
      WRITE(LUPRI,'(/A/)') '** SPECTROSCOPIC CONSTANTS:'
      WRITE(ILOGG,'(/A/)') '** SPECTROSCOPIC CONSTANTS:'
      WRITE(LUPRI,'(A)') 
     &  'Polynomial fit: Give order of polynomial '
      READ(LUCMD,*) NORDER
      IF(NORDER.LT.2) THEN
        WRITE(LUPRI,'(A)') 
     &    '* This program requires at least a quartic fit.'
        STOP 'Too low order fit !'
      ENDIF
      IF(NORDER.LT.4) THEN
        WRITE(LUPRI,'(A)') 
     &    '* WARNING: a quartic fit is needed for full functionality.'
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
      CALL POLSVD(NDIM,NP,A,B,X,Y,C,D,CHISQ,ISKIP)
C* Estimate fit:
      WRITE(LUPRI,'(3X,A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(LUPRI,'(3X,A)') '* Coefficients:'
      WRITE(LUPRI,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',
     &      C(I),I=1,NDIM)  
      WRITE(LUPRI,'(3X,A,E9.4)') '* Chi square :  ',CHISQ
      WRITE(LUPRI,'(3X,A,E9.4)') '* Chi square per point:  ',CHISQ/NP
      WRITE(LUPRI,'(3X,A,I3)') '* Number of singularities(SVD): ',ISKIP
C
      WRITE(ILOGG,'(72A1)') ('=',I=1,72)
      WRITE(ILOGG,'(3X,A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(ILOGG,'(3X,A)') '* Coefficients:'
      WRITE(ILOGG,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',
     &      C(I),I=1,NDIM)  
      DO I = 1,NP
        YP = (POLVAL(NORDER,C,X(I))-Y(I))/Y(I)
        WRITE(LUPRI,'(3X,F6.3,4X,E12.4)') XTANG*X(I),YP
      ENDDO
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square :  ',CHISQ
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square per point:  ',CHISQ/NP
      WRITE(ILOGG,'(3X,A,I3)') '* Number of singularities(SVD): ',ISKIP
C
C     Find minimum
C     ============
C
C* Find smallest y-value and take corresponding x-value as start value
C* for Newton-Raphsons method      
      IYOPT = IDMIN(NP,Y,1)
      XOPT = X(IYOPT)
      CALL NEWRAP(NORDER,DTOL,C,XOPT,NITER,IERR)
      IF(IERR.EQ.1) THEN
        WRITE(ILOGG,'(A)') 'Newtons method failed'
        STOP 'Did not find minimum. Newtons method failed...'
      ELSE
        YOPT = POLVAL(NORDER,C,XOPT)
        TEMP = XOPT*XTANG
        WRITE(LUPRI,'(A,F18.5,A)') 
     &    '* Local minimum   :',TEMP,' Angstroms'
        WRITE(ILOGG,'(A,F18.5,A)') 
     &    '* Local minimum   :',TEMP,' Angstroms'
        WRITE(LUPRI,'(A,1P,E18.10,A)') 
     &    '* Expected energy :',YOPT,' Hartrees'
        WRITE(ILOGG,'(A,1P,E18.10,A)') 
     &    '* Expected energy :',YOPT,' Hartrees'
      ENDIF
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
      V2 = POL2DER(NORDER,C,XOPT)
      TEMP  = V2*XTJ/XTM/XTM
      WRITE(LUPRI,'(A,1P,E18.4,A)') 
     &   '* Force constant  :',TEMP, ' N/m'
      WRITE(ILOGG,'(A,1P,E18.4,A)') 
     &   '* Force constant  :',TEMP, ' N/m'
      OMEGA = SQRT(V2/UM)
      WRITE(LUPRI,'(A,1P,E12.5,A)') 
     &   '* Frequency       :',OMEGA*XTHZ,' Hz',
     &   '                   ',OMEGA*XTCM,' cm-1'
      WRITE(ILOGG,'(A,1P,E12.5,A)') 
     &   '* Frequency       :',OMEGA*XTHZ,' Hz',
     &   '                   ',OMEGA*XTCM,' cm-1'
      XDEV=D1/SQRT(D2*UM*OMEGA)
      XMIN=XOPT-XDEV
      XMAX=XOPT+XDEV
      TEMP=XDEV*XTANG
      WRITE(LUPRI,'(A,F10.4,A/A,A1,F10.4,A1,F10.4,A1)') 
     &   '* Mean displacement in harmonic ground state:',TEMP,
     &   ' Angstroms',
     &   '  corresponding to interval: ',
     &   '[',XMIN*XTANG,',',XMAX*XTANG,']'
      WRITE(ILOGG,'(A,F10.4,A/A,A1,F10.4,A1,F10.4,A1)') 
     &   '* Mean displacement in harmonic ground state:',TEMP,
     &   ' Angstroms',
     &   '  corresponding to interval: ',
     &   '[',XMIN*XTANG,',',XMAX*XTANG,']'
      IOUT=0
      DO I = 1,NP
        IF(X(I).LT.XMIN.OR.X(I).GT.XMAX) IOUT=IOUT+1
      ENDDO
      IF(IOUT.GT.0) THEN
        WRITE(LUPRI,'(A,I5,A/A)')
     &   '* WARNING : ',IOUT,' points lie outside interval...',
     &   'This may reduce the quality of your spectroscopic constants.'
        WRITE(ILOGG,'(A,I5,A/A)')
     &   '* WARNING : ',IOUT,' points lie outside interval...',
     &   'This may reduce the quality of your spectroscopic constants.'
      ENDIF
      IF(NORDER.LT.4) GOTO 666 
C
C     Anharmonic constant
C     ===================
C
      V3 = POLNDER(NORDER,C,XOPT,3)
      V4 = POLNDER(NORDER,C,XOPT,4)
      TEMP = 7.0D0*V3*V3/V2/9.0D0
      T2 = 100.0D0*(TEMP - V4)/V4
      WRITE(ILOGG,'(A,D20.10)')
     &   '* Third derivative                   : ', V3,
     &   '* Fourth derivative                  : ', V4
      WRITE(ILOGG,'(A,D20.10, F9.2,A)')
     &   '* Fourth derivative (Morse estimate) : ', TEMP,T2,' \\%'
      TEMP = D1/UM/OMEGA
      T2   = TEMP*TEMP
      WXE = (5.0D0*TEMP*T2*V3*V3/48.0D0/OMEGA)-(V4*T2/16.0D0)
      WRITE(LUPRI,'(A,1P,E12.5,A)') 
     &   '* Omega*x_e       :',WXE*XTCM,' cm-1'
      WRITE(ILOGG,'(A,1P,E12.5,A)') 
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
        V2N = POL2DER(NORDER,C,XNEW)
        TEMP  = V2N*XTJ/XTM/XTM
        WRITE(LUPRI,'(A,1P,E18.4,A)') 
     &   '* Force constant  :',TEMP, ' N/m'
        WRITE(ILOGG,'(A,1P,E18.4,A)') 
     &   '* Force constant  :',TEMP, ' N/m'
        OMEGA = SQRT(V2N/UM)
        WRITE(LUPRI,'(A,1P,E12.5,A)') 
     &   '* Frequency       :',OMEGA*XTHZ,' Hz',
     &   '                   ',OMEGA*XTCM,' cm-1'
        WRITE(ILOGG,'(A,1P,E12.5,A)') 
     &   '* Frequency       :',OMEGA*XTHZ,' Hz',
     &   '                   ',OMEGA*XTCM,' cm-1'
        V3N = POLNDER(NORDER,C,XNEW,3)
        V4N = POLNDER(NORDER,C,XNEW,4)
        TEMP = 7.0D0*V3N*V3N/V2N/9.0D0
        T2 = 100.0D0*(TEMP - V4N)/V4N
        WRITE(ILOGG,'(A,D20.10)')
     &   '* Third derivative                   : ', V3N,
     &   '* Fourth derivative                  : ', V4N
        WRITE(ILOGG,'(A,D20.10, F9.2,A)')
     &   '* Fourth derivative (Morse estimate) : ', TEMP,T2,' \\%'
        TEMP = D1/UM/OMEGA
        T2   = TEMP*TEMP
        WXE = (5.0D0*TEMP*T2*V3N*V3N/48.0D0/OMEGA)-(V4N*T2/16.0D0)
        WRITE(LUPRI,'(A,1P,E12.5,A)') 
     &   '* Omega*x_e       :',WXE*XTCM,' cm-1'
        WRITE(ILOGG,'(A,1P,E12.5,A)') 
     &   '* Omega*x_e       :',WXE*XTCM,' cm-1'
      ELSEIF(REPLY.EQ.'n') THEN
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
C
