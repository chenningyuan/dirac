      SUBROUTINE PRTITL(LINES)
*
* Print title cards
*
      CHARACTER*72 LINES
      DIMENSION LINES(3)
      CHARACTER*80 STARS
*
      STARS(1:1) = ' '
      DO 80 I = 2, 80
        STARS(I:I) = '*'
   80 CONTINUE
      WRITE(6,'(A)') STARS
      WRITE(6,'(3A)') ' *  ',LINES(1),'  *'
      WRITE(6,'(3A)') ' *  ',LINES(2),'  *'
      WRITE(6,'(3A)') ' *  ',LINES(3),'  *'
      WRITE(6,'(A)') STARS
*
      RETURN
      END
