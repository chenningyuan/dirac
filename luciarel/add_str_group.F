      SUBROUTINE ADD_STR_GROUP(NSTADD,IOFADD,ISTADD,NSTB,NSTA,
     &                         ISTRING,IELOF,NELADD,NELTOT)
*
* Part of assembling strings in individual types to
* super group of strings
*
*. Copying strings belonging to a given type to supergroup of strings
*
* Jeppe Olsen, for once improving performance of LUCIA
*
*.Input
* =====
* NSTADD : Number of strings to be added
* IOFADD : First string to be added
* ISTADD : Strings to be added
* NSTB   : Number of strings belonging to lower gasspaces
* NSTA   : Number of strings belonging to higher gasspaces
* ISTRING: Supergroup of strings under construction
* IELOF  : Place of first electron to be added
* NELADD : Number of electrons to be added
* NELTOT : Total number of electrons
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ISTADD(NELADD,*)
*. Input and output
      DIMENSION ISTRING(NELTOT,*)
*
      IF(NSTA.GT.1) THEN
        DO IISTR = 1,NSTADD
*. Address of A(1,IISTR,1)
*. A(I(after),Igas,I(before))
          IOFF1 = (IISTR-1)*NSTA + 1
          DO ISTB = 1, NSTB
*. Address of A(1,IISTR,ISTB)
            IOFF2 = IOFF1 + (ISTB-1)*NSTADD*NSTA
            DO ISTA = 1, NSTA
              DO JEL = 1, NELADD
                ISTRING(IELOF-1+JEL,IOFF2-1+ISTA)
     &        = ISTADD(JEL,IOFADD-1+IISTR)
              END DO
            END DO
          END DO
        END DO
      ELSE IF (NSTA .EQ. 1 ) THEN
*. Address of A(1,IISTR,1)
*. A(I(after),Igas,I(before))
        DO ISTB = 1, NSTB
          IOFF0 = (ISTB-1)*NSTADD
          DO IISTR = 1,NSTADD
*. Address of A(1,IISTR,ISTB)
            IOFF2 = IISTR  + IOFF0
            DO JEL = 1, NELADD
              ISTRING(IELOF-1+JEL,IOFF2)
     &      = ISTADD(JEL,IOFADD-1+IISTR)
            END DO
          END DO
        END DO
      END IF
*
      NTEST = 0000
      if (NTEST.ge.1000) then
        write(6,*) 'ADD_STR_GROUP speaking: Input string array'
        call iwrtma(ISTADD,NELADD,IOFADD-1+NSTADD,
     &                     NELADD,IOFADD-1+NSTADD)
      end if
*
      RETURN
      END
