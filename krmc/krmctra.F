C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rtractl */
      SUBROUTINE RTRACTL(ITRLVL,WORK,LWORK)
C*****************************************************************************
C
C     Input : ITRLVL - level of integrals transformation (see below)
C
C     Output: 4-index transformed integrals on disk.
C
C      
C     ITRLVL = -1: Transform all integrals (pq|rs)
C                     (which is more than we need, but for a quick start
C                      this is good)
C     ITRLVL =  0: Transform all integrals (uv|xy) for CI & MC energy
C     ITRLVL =  1: (pv|xy) for CI & MC gradients
C     ITRLVL =  2: (pq|xy) for Newton-Raphson and NEO CI & MC
C
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbopt.h>
#include <dgroup.h>
#include <dcborb.h>
C
      CHARACTER SZTRA(4,2)*80
      DIMENSION WORK(*)
C
#include <memint.h>
C
C
      CALL QENTER('RTRACTL')
C
      IF ( IPROPT .GT. 5 ) THEN
         CALL HEADER('Output from RTRACTL',-1)
         WRITE(LUPRI,'(A,I2)')
     $        '* Transformation level (ITRLVL) = ',ITRLVL
         WRITE(LUPRI,'(3A,I1,A)')
     $        '* Integrals in transformation = ',OPT_INTTYP,
     $        ' (INTFLG=',IOPT_INTFLG,')'
      END IF
C
      IF ( ITRLVL .LT. -1 .OR. ITRLVL .GT. 4 ) THEN
         WRITE(LUPRI,'(2(/,A),I2)')
     $        '*** ERROR in RTRACTL ***',
     $        'Illegal transformation level (ITRLVL) = ',ITRLVL
         CALL QUIT('*** ERROR in RTRACTL ***')
      END IF
C
      IF ( ITRLVL .GE. 0 ) THEN
         NGENERAL = ITRLVL
      ELSE
         NGENERAL = 4
      END IF
      DO I = 1, NGENERAL
         DO J = 1, NFSYM
            IF ( NESH (J) .GT. 0 ) THEN
               IF ( NPSH (J) .GT. 0 ) THEN
                  WRITE(SZTRA(I,J),'(I4,A,I4)')
     $                 -NPSH(J), '..-1, 1..', NESH(J)
               ELSE
                  WRITE(SZTRA(I,J),'(A,I4)')
     $                 '1..', NESH(J)
               END IF
            ELSE
               IF ( NPSH ( J ) .GT. 0 ) THEN
                  WRITE(SZTRA(I,J),'(I4,A)')
     $                 -NPSH(J), '..-1'
               ELSE
                  SZTRA(I,J) = ' '
               END IF
            END IF
         END DO
      END DO
      DO I = NGENERAL+1, 4
         DO J = 1, NFSYM
            IF ( NASH (J) .GT. 0 ) THEN
               WRITE(SZTRA(I,J),'(A,I4)')
     $              '1..', NASH(J)
            ELSE
               SZTRA(I,J) = ' '
            END IF
         END DO
      END DO
C
C     Print orbital strings
C
      IF ( IPROPT .GE. 5 ) THEN
         DO I = 1, 4
            DO J = 1,NFSYM
               WRITE(LUPRI,'(A,I1,A,I1,2A)')
     $              '* Orbital string ',I,
     $              ' for irrep ',J,' = ',SZTRA(I,J)
            END DO
         END DO
      END IF
      write(6,*) 'nopqfq ',opt_nopfq
C
C     Transfer information to moltra common blocks
C
      CALL SETUP_DCBTRA(SZTRA,IOPT_INTFLG,OPT_NOPFQ,0)
C
C     Call 4-index transformation driver
C
C
      CALL PAMTRA(WORK(KFREE),LFREE)
C
C
      CALL QEXIT('RTRACTL')
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rtractl */
      SUBROUTINE SETUP_DCBTRA(SZ,INTFLG,LNOPAIR,IPRINT)
C*****************************************************************************
C
C     Input : SZ     - orbitals strings
C             IPRINT - print level
C             INTFLG - integral classes to include in transformation
C
C     Output: none
C
C     Sets the common blocks used in the 4-index transformation.
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER ( DM1 = -1.0D00 )
C
      CHARACTER SZ(4,2)*(*)
      LOGICAL LNOPAIR
C
#include <dcbtra.h>
C     
C     /DCBTC/
C
      DO I = 1, 4
         DO J = 1, 2
            TRA4_INDSTR(I,J) = SZ(I,J)
         END DO
      END DO
      DO I = 1,2
         TRA_INDSTR(I) = ' '
         TRA_CORSTR(I) = ' '
         DO J = 1,2
            TRA2_INDSTR(I,J) = ' '
         END DO
      END DO
C
C     /DCBTRI/
C
      ITRA_INTFLG = INTFLG
      ITRA_INTFL2 = INTFLG
      ITRA_INTFL4 = INTFLG
      IPRTRA = IPRINT
      ISTRAT = 4
C
C     /DCBTRL/
C
      TRA_ANTIS = .FALSE.
      TRA_MSOUT  = .FALSE.
      NO2IND = .TRUE.
      NO4IND = .FALSE.
      NOPAIR = LNOPAIR
      PRPTRA = .FALSE.
      RCORBS = .FALSE.
      FKRMC = .TRUE.
C      
C     /DCBTRR/
C
      THROUT = 1.0D-16
      SCRTRA = DM1
C
C
      RETURN
      END
