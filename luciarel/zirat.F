      SUBROUTINE ZIRAT
*
* Ratio between real and integer
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      include 'SysDef.inc'
#include "machine.inc"
*. Output
#include "irat.inc"
*. For CRAY or explicit 64 bit architecture : IRAT = 1
c     IF(MACHINE(1:4).EQ.'CRAY'.OR.MACHINE(1:2).EQ.'64') THEN
c       IRAT = 1
c     ELSE
*. For all other architectures : 2 integers per real
c       IRAT = 2
c     END IF
*
CTF No variability in thsi version.
*   We always use RTOI = 1, but in some MOLCAS interface routines,
*   RTOI is set to 2.
c     IRAT = RTOI
      IRAT = 1
*
      WRITE(6,*)
      WRITE(6,'(1H ,6X,A,I2)')
     &'Ratio between Integer and Real word length ', IRAT
      WRITE(6,*)
*
      RETURN
      END
