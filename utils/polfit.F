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
      PROGRAM POLFIT
#include <implicit.h>
#include <priunit.h>
      PARAMETER(MORDER=16,MP=50)
      PARAMETER(MDIM=MORDER+1)
      PARAMETER(DTOL = 1.0D-4,DM1 = -1.0D0,D1 = 1.0D0)
      PARAMETER(NITER=10)
      DIMENSION A(MP,MDIM),B(MDIM),C(MDIM),D(MDIM)
      DIMENSION X(MP),Y(MP)
      CHARACTER RESFIL*6,REPLY*1,POTFIL*6
      IUNIT = 1
      ILOGG = 2
      WRITE(LUPRI,'(A)') '* Name output file (A6)'
      READ(LUCMD,'(A6)') RESFIL
      OPEN(ILOGG,FILE=RESFIL,STATUS='UNKNOWN',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL')
      WRITE(LUPRI,'(A)') 'Name input file with potential curve(A6)'
      READ(LUCMD,'(A)') POTFIL
      OPEN(IUNIT,FILE=POTFIL,STATUS='OLD',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL',ERR=5)
      GOTO 6
 5    CONTINUE
      STOP 'Error opening POTFIL'
 6    CONTINUE
      NP = 0
   10 CONTINUE  
C Read points
      NP = NP + 1
      IF(NP.GT.MP) STOP 'Too many points !!!'
      READ(IUNIT,*,END=20) X(NP),Y(NP)
      GOTO 10
 20   CONTINUE
      NP = NP-1
      WRITE(LUPRI,'(A,I5)') '* Number of points read:',NP
      WRITE(ILOGG,'(A,A6)') '* Points read from ',POTFIL
      WRITE(ILOGG,'(3X,F6.3,4X,F20.10)') (X(I),Y(I),I=1,NP)
 30   CONTINUE
      WRITE(*,'(A)') 'Give order of polynomial:'
      READ(LUCMD,*) NORDER
      WRITE(*,'(I3)') NORDER
      IF(NP.LE.NORDER) THEN
        WRITE(*,'(A,I3)') 
     +    'Too few points for a polynomial fit of order',NORDER
        WRITE(*,'(A)') '  Try again.'
        GOTO 30
      ENDIF
      NDIM = NORDER + 1
      CALL POLSVD(NDIM,NP,A,B,X,Y,C,D,CHISQ,ISKIP)
C      CALL CURVE(NORDER,NP,A,B,X,Y,VAR,IBUF)
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
        WRITE(LUPRI,'(3X,F6.3,4X,E12.4)') X(I),YP
      ENDDO
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square :  ',CHISQ
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square per point:  ',CHISQ/NP
      WRITE(ILOGG,'(3X,A,I3)') '* Number of singularities(SVD): ',ISKIP
      WRITE(LUPRI,'(A)') '* Do you want to find min. of curve (y/n) ?'
      READ(LUCMD,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 666
C* Find smallest y-value and take corresponding x-value as start value
C* for Newton-Raphsons method      
C      XVAL = X(MIND(NP,Y))
      IYMIN = IDMIN(NP,Y,1)
      XMIN = X(IYMIN)
      CALL NEWRAP(NORDER,DTOL,C,XMIN,NITER,IERR)
C* Give coordinates of local minimum
      IF(IERR.EQ.1) THEN
        WRITE(*,'(A)') 'Newtons method failed'
      ELSE
        YMIN = POLVAL(NORDER,C,XMIN)
        WRITE(LUPRI,'(A,F17.10,A,F17.10,A)') '* Local minimum:  (',
     +        XMIN,',',YMIN,')'
        WRITE(ILOGG,'(A,F17.10,A,F17.10,A)') '* Local minimum:  (',
     +        XMIN,',',YMIN,')'
      ENDIF
C* Calculate derivatives:
      WRITE(LUPRI,'(A)') 
     &  '* Do you want derivatives of curve at minimum(y/n) ?'
      READ(LUCMD,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 666
      TMP = POLNDER(NORDER,C,XMIN,2)
      WRITE(LUPRI,'(A,E18.10)') '2nd derivative at minimum: ',TMP
      WRITE(ILOGG,'(A,E18.10)') '2nd derivative at minimum: ',TMP
      TMP = POLNDER(NORDER,C,XMIN,3)
      WRITE(LUPRI,'(A,E18.10)') '3rd derivative at minimum: ',TMP
      WRITE(ILOGG,'(A,E18.10)') '3rd derivative at minimum: ',TMP
      TMP = POLNDER(NORDER,C,XMIN,4)
      WRITE(LUPRI,'(A,E18.10)') '4th derivative at minimum: ',TMP
      WRITE(ILOGG,'(A,E18.10)') '4th derivative at minimum: ',TMP
 666  CONTINUE
      WRITE(LUPRI,'(A)') '* Do you want a new polynomial fit ?(y/n)'
      READ(LUCMD,'(A1)') REPLY
      IF(REPLY.EQ.'y') GOTO 30
      END
C
C
