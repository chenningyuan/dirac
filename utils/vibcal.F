!
!      Copyright (c) 2008 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.sdu.dk
!
      PROGRAM VIBCAL
C********************************************************************
C
C     For a given normal mode this program will calculate the
C     expectation value for a property in the five lowest
C     vibrational states
C
C     Written by T. Saue Feb 13 2007
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "codata.h"
      CHARACTER POTFIL*30,PRPFIL*30,REPLY*1
C      PARAMETER(V2=4.6822365486E-01,V3=-1.1951834335E+00,U=9.703)
C      PARAMETER(V4=6.2052157257E+00)
      PARAMETER(DP5=0.5D0,D1=1.0D0,DP16=0.0625D0,XTCM=2.19474625D+05)
      DIMENSION B(0:4),H(0:4)
      CALL TITLER('VIBCAL for properties : Written by T. Saue ','*',110)
      WRITE(LUPRI,'(A)') 
     &    'Name input file (A30) with reduced mass ',
     &    'and force constants V2,V3,V4'
      READ(LUCMD,'(A)') POTFIL
      WRITE(LUPRI,'(A)') 
     &    'Name input file (A30) with property derivatives P0,P1,P2'
      READ(LUCMD,'(A)') PRPFIL
      IPOT = 1
      OPEN(IPOT,FILE=POTFIL,STATUS='OLD',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL',ERR=20)
      i=LNBLNK(POTFIL)
      IPRP = 2
      OPEN(IPRP,FILE=PRPFIL,STATUS='UNKNOWN',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL',ERR=40)
      j=LNBLNK(PRPFIL)
      ILOGG = 3
      OPEN(ILOGG,FILE=POTFIL(1:i)//'_'//PRPFIL(1:j)//'.vib',
     &     STATUS='UNKNOWN',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL')
C
      READ(IPOT,*,ERR=30,END=30) U,V2,V3,V4
      WRITE(ILOGG,*) 'Reduced mass: ',U
      WRITE(ILOGG,*) '2nd potential derivative: ', V2
      WRITE(ILOGG,*) '3rd potential derivative: ', V3
      WRITE(ILOGG,*) '4th potential derivative: ', V4
C
      READ(IPRP,*,ERR=50,END=50) P0,P1,P2
      WRITE(ILOGG,*) '0th property derivative: ', P0
      WRITE(ILOGG,*) '1st property derivative: ', P1
      WRITE(ILOGG,*) '2nd property derivative: ', P2
      UM = U*XFAMU
      OMEGA = SQRT(V2/UM)
      FAC=D1/UM/OMEGA
      A0 = P0
      A1 = DP5*P2*FAC
      A2 = -DP5*P1*V3*FAC*FAC/OMEGA
      A3 = -DP16*P2*V4*FAC*FAC*FAC/OMEGA
      WRITE(ILOGG,'(A,A)')
     &  '    n  Term 1 (P0)     Term 2 (P2)     Term 3 (P1 V3) ',
     &  ' Term 4 (P2 V4)  Harmonic cont.  Anharm. cont.   Total value'
      DO N = 0,4
         V = N + DP5
         W = N*N+N+DP5
         HAR = A0 + A1*V
         ANH = A2*V + A3*W
         TOT = HAR + ANH 
         B(N) = TOT         
         H(N) = HAR
         WRITE(ILOGG,'(I5,7E16.8)') N,A0,A1*V,A2*V,A3*W,HAR,ANH,TOT
      ENDDO
      DIF = B(1)-B(0)
      DIF2 = DIF*XTHZ*1000.0D0
      WRITE(ILOGG,'(A,E16.8,A,F10.1,A)') 'Difference :',DIF, ' au = ',
     &     DIF2, ' mHz'
      WRITE(ILOGG,'(A,F10.1,A)') 'PV shift:',2.0D0*DIF2,' mHz'
      WRITE(ILOGG,'(A)') 'Harmonic model:'
      DIF = H(1)-H(0)
      DIF2 = DIF*XTHZ*1000.0D0
      WRITE(ILOGG,'(A,E16.8,A,F10.1,A)') 'Difference :',DIF, ' au = ',
     &     DIF2, ' mHz'
      WRITE(ILOGG,'(A,F10.1,A)') 'PV shift:',2.0D0*DIF2,' mHz'
      GOTO 100
 20   STOP 'Error opening POTFIL '
 30   STOP 'Error reading POTFIL '
 40   STOP 'Error opening PRPFIL '
 50   STOP 'Error reading PRPFIL '
 100  CONTINUE
      END

