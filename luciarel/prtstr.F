      SUBROUTINE PRTSTR(ISTR,NEL,NSTR)
*
* Print NSTR strings each containing NEL electrons
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ISTR(NEL,NSTR)
*
      DO JSTR = 1, NSTR
        WRITE(6,'(1H0,A,I6,A,4X,10(2X,I4),/,(1H ,19X,10(2X,I4)))' )
     &   ' String ',JSTR,' : ',(ISTR(IEL,JSTR),IEL=1,NEL)
      END DO
*
      RETURN
      END
