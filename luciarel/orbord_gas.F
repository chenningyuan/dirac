      SUBROUTINE ORBORD_GAS(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB,NGSOBT,
     &                  NOCOBS,NTOOBS,NTOOB,
     &                  IREOST,IREOTS,ISFTO,ITFSO,IBSO,
     &                  NTSOB,IBTSOB,ITSOB,NOBPTS,IOBPTS,
     &                  ISFSO,ITFTO,NOBPT,IPRNT,IORD,ICOMP)
*
* Obtain Reordering arrays for orbitals
* ( See note below for assumed ordering )
*
*
* GAS version
*
* =====
* Input
* =====
*  NSMOB  : Number of orbital symmetries
*  MXPOBS : Max number of orbital symmetries allowed by program
*  MXPNGAS: Max number of GAS spaces allowed by program
*  NGAS   : Number of GAS spaces
*  NGSOB  : Number of GAS orbitals per symmetry and space
*  NGSOBT : Number of GAS orbitals per space
*  NOCOBS : Number of occupied orbitals per symmetry
*  NTOOBS : Number of orbitals per symmetry,all types
*  IORD   : Order in which orbital symmetries are filled
*
* ======
* Output
* ======
*  IREOST : Reordering array symmetry => type
*  IREOTS : Reordering array type     => symmetry
*  ISFTO  : Symmetry array for type ordered orbitals
*  ITFSO  : Type array for symmetry ordered orbitals( not activated )
*  IBSO   : First orbital of given symmetry ( symmetry ordered )
*  NOBPTS : Number of orbitals per subtype and symmetry
*  IOBPTS : Offsets for orbitals of given subtype and symmetry
*           ordered according to input integrals
*
* ISFSO  : Symmetry of orbitals, symmetry ordereing
* ITFTO  : Type of orbital, type ordering
*
* Jeppe Olsen, Winter 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION NGSOB(MXPOBS,MXPNGAS),NOCOBS(*),NTOOBS(*)
      DIMENSION NGSOBT(MXPNGAS)
      DIMENSION IORD(MXPOBS,*)
*. Output
      DIMENSION IREOST(*),IREOTS(*),ISFTO(*),ITFSO(*),IBSO(*)
      DIMENSION ISFSO(*),ITFTO(*)
      DIMENSION NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*)
      DIMENSION NOBPT(MXPNGAS)

* ==========================
* Note on order of orbitals
* ==========================
*
* The orbitals are supposed to be imported ordered symmetry-type
* ordered as
*
* Loop over symmetries of orbitals
*  Loop over GAS spaces
*   Loop over orbitals of this sym and GAS
*   End of Loop over orbitals
*  End of Loop over Gas spaces
* End of loop over symmetries
*
* Internally the orbitals are reordered to type symmetry order
* where the outer loop is over types and the inner loop is
* over symmetries, i.e.
*
* Loop over GAS spaces
*  Loop over symmetries of orbitals
*   Loop over orbitals of this sym and GAS
*   End of Loop over orbitals
*  End of loop over symmetries
* End of Loop over Gas spaces
*
*. 1:  Construct ISFTO, ITFTO, IREOST,IREOTS,NOBPTS,IOBPTS
*
      NTEST = 0
      NTEST = MAX(IPRNT,NTEST)
*
      if (NTEST.ge.10) then
        WRITE(6,*) ' ICOMP = ', ICOMP
        WRITE(6,*) ' Four elements of relevant IORD row '
        WRITE(6,*)
     &   IORD(ICOMP,1),IORD(ICOMP,2),
     &   IORD(ICOMP,3),IORD(ICOMP,4)
      end if
*
      ITSOFF = 1
      DO IGAS = 1, NGAS
        DO IISYM = 1, NSMOB
          ISYM = IORD(ICOMP,IISYM)
          IF(IISYM.EQ.1) THEN
            IBSSM = 1
          ELSE
            IBSSM = IBSSM + NTOOBS(IORD(ICOMP,IISYM-1))
          END IF
          NPREV = 0
          DO JGAS = 1, IGAS-1
            NPREV = NPREV + NGSOB(ISYM,JGAS)
          END DO
          IADD = 0
          NOBPTS(IGAS,ISYM) = NGSOB(ISYM,IGAS)
          IOBPTS(IGAS,ISYM) = ITSOFF
          DO IORB = ITSOFF,ITSOFF+NGSOB(ISYM,IGAS)-1
            IADD = IADD + 1
            IREOTS(IORB) = IBSSM-1+NPREV+IADD
            IREOST(IBSSM-1+NPREV+IADD) = IORB
            ITFTO(IORB) = IGAS
            ISFTO(IORB) = ISYM
          END DO
          ITSOFF = ITSOFF + NGSOB(ISYM,IGAS)
        END DO
      END DO
*
* 2 : ISFSO,ITFSO
*
      ISTOFF = 1
      DO IISYM = 1, NSMOB
        ISYM = IORD(ICOMP,IISYM)
        DO IGAS = 1, NGAS
          DO IORB = ISTOFF,ISTOFF+NGSOB(ISYM,IGAS)-1
            ISFSO(IORB) = ISYM
            ITFSO(IORB) = IGAS
          END DO
          ISTOFF = ISTOFF + NGSOB(ISYM,IGAS)
        END DO
      END DO
*
* 3 IBSO, NOBPT
*
      IOFF = 1
      DO IISM = 1, NSMOB
       ISM = IORD(ICOMP,IISM)
       IBSO(ISM) = IOFF
       IOFF = IOFF + NTOOBS(ISM)
      END DO
      DO IGAS = 1, NGAS
        NOBPT(IGAS) = NGSOBT(IGAS)
      END DO
*
      IF( NTEST .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' Output from ORBORD '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        WRITE(6,*) ' Symmetry of orbitals , type ordered '
        CALL IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' Symmetry => type reordering array '
        CALL IWRTMA(IREOST,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' Type => symmetry reordering array '
        CALL IWRTMA(IREOTS,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' IBSO array '
        CALL IWRTMA(IBSO,1,NSMOB,1,NSMOB)
*
        WRITE(6,*) ' NOBPTS '
        CALL IWRTMA(NOBPTS,NGAS,NSMOB,MXPNGAS,MXPOBS)
        WRITE(6,*) ' NOBPT '
        CALL IWRTMA(NOBPT,NGAS,1,MXPNGAS,1)
        WRITE(6,*) ' IOBPTS '
        CALL IWRTMA(IOBPTS,NGAS,NSMOB,MXPNGAS,MXPOBS)
*
        WRITE(6,*) ' ISFTO array : '
        CALL IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' ITFSO array : '
        CALL IWRTMA(ITFSO,1,NTOOB,1,NTOOB)
*
        WRITE(6,*) ' ISFSO array : '
        CALL IWRTMA(ISFSO,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' ITFTO array : '
        CALL IWRTMA(ITFTO,1,NTOOB,1,NTOOB)
      END IF
*

      RETURN
      END
