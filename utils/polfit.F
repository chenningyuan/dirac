      PROGRAM POLFIT
#include <implicit.h>
#include <priunit.h>
      PARAMETER(MORDER=10,MP=50)
      PARAMETER(MDIM=MORDER+1,MNDIM=(MDIM*(MDIM+1)/2))
      PARAMETER(DTOL = 1.0D-4,DM1 = -1.0D0,D1 = 1.0D0)
      PARAMETER(NITER=10)
      DIMENSION A(MNDIM),B(MDIM),IBUF(MDIM)
      DIMENSION X(MP),Y(MP)
      CHARACTER RESFIL*6,REPLY*1,POTFIL*6
      IUNIT = 1
      ILOGG = 2
      WRITE(LUPRI,'(A)') '* Name output file (A6)'
      READ(LUCMD,'(A6)') RESFIL
      OPEN(ILOGG,NAME=RESFIL,STATUS='UNKNOWN',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL')
      WRITE(LUPRI,'(A)') 'Name input file with potential curve(A6)'
      READ(LUCMD,'(A)') POTFIL
      OPEN(IUNIT,NAME=POTFIL,STATUS='OLD',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL',ERR=5)
      GOTO 6
 5    CONTINUE
      STOP 'Error opening POTFIL'
 6    CONTINUE
      NPOINTS = 0
   10 CONTINUE  
C Read points
      NPOINTS = NPOINTS + 1
      IF(NPOINTS.GT.MP) STOP 'Too many points !!!'
      READ(IUNIT,*,END=20) X(NPOINTS),Y(NPOINTS)
      GOTO 10
 20   CONTINUE
      NPOINTS = NPOINTS-1
      WRITE(LUPRI,*) '* Points read:',NPOINTS
 30   CONTINUE
      WRITE(*,'(A)') 'Give order of polynomial:'
      READ(LUCMD,*) NORDER
      WRITE(*,'(I3)') NORDER
      IF(NPOINTS.LE.NORDER) THEN
        WRITE(*,'(A,I3)') 
     +    'Too few points for a polynomial fit of order',NORDER
        WRITE(*,'(A)') '  Try again.'
        GOTO 30
      ENDIF
      NDIM = NORDER + 1
      CALL CURVE(NORDER,NPOINTS,A,B,X,Y,VAR,IBUF)
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
      WRITE(LUPRI,'(A)') '* Do you want to find min. of curve (y/n) ?'
      READ(LUCMD,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 666
C* Find smallest y-value and take corresponding x-value as start value
C* for Newton-Raphsons method      
C      XVAL = X(MIND(NPOINTS,Y))
      IYMIN = IDMIN(NPOINTS,Y,1)
      XMIN = X(IYMIN)
      CALL NEWRAP(NORDER,DTOL,B,XMIN,NITER,IERR)
C* Give coordinates of local minimum
      IF(IERR.EQ.1) THEN
        WRITE(*,'(A)') 'Newtons method failed'
      ELSE
        YMIN = POLVAL(NORDER,B,XMIN)
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
      WRITE(ILOGG,'(I5,E18.10)') (I,POLNDER(NORDER,B,XMIN,I),I=2,4)
 666  CONTINUE
      WRITE(LUPRI,'(A)') '* Do you want a new polynomial fit ?(y/n)'
      READ(LUCMD,'(A1)') REPLY
      IF(REPLY.EQ.'y') GOTO 30
      END
C
C
