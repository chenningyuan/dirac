      SUBROUTINE NATORB(RHO1,NSMOB,NTOPSM,NACPSM,NINPSM,
     &                  ISTOB,XNAT,RHO1SM,OCCNUM,
     &                  NACOB,SCR,IPRDEN)
*
* Obtain natural orbitals in symmetry blocks
*
* Jeppe Olsen, June 1994
*              Modification, Oct 94
*              Last modification, Feb. 1998 (reorder deg eigenvalues)
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION RHO1(NACOB,NACOB)
      DIMENSION ISTOB(*)
      DIMENSION NTOPSM(NSMOB), NACPSM(NSMOB),NINPSM(NSMOB)
*. Output
      DIMENSION RHO1SM(*),OCCNUM(*),XNAT(*)
*. Scratch ( Largest symmetry block )
      DIMENSION SCR(*)
*
      NTESTL = 0
      NTEST = MAX(NTESTL,IPRDEN)
*. IOBOFF : Offset for active orbitals in symmetry order
      DO ISMOB = 1, NSMOB
        IF(ISMOB.EQ.1) THEN
          IOBOFF = NINPSM(1)+1
          IMTOFF = 1
        ELSE
          IOBOFF =
     &    IOBOFF + NTOPSM(ISMOB-1)-NINPSM(ISMOB-1)+NINPSM(ISMOB)
          IMTOFF = IMTOFF + NACPSM(ISMOB-1)**2
        END IF
        LOB = NACPSM(ISMOB)
*
*. Extract symmetry block of density matrix
*
        DO IOB = IOBOFF,IOBOFF + LOB-1
           DO JOB = IOBOFF,IOBOFF + LOB-1
*. Corresponding type indeces
             IOBP = ISTOB(IOB)
             JOBP = ISTOB(JOB)
             RHO1SM(IMTOFF-1+(JOB-IOBOFF)*LOB+IOB-IOBOFF+1)
     &     = RHO1(IOBP,JOBP)
           END DO
        END DO
*
        IF(NTEST.GE.2 ) THEN
          WRITE(6,*)
          WRITE(6,*) ' Density matrix for symmetry  = ', ISMOB
          WRITE(6,*) ' ======================================='
          WRITE(6,*)
          CALL WRTMAT(RHO1SM(IMTOFF),LOB,LOB,LOB,LOB)
        END IF
*. Pack and diagonalize
        CALL TRIPAK(RHO1SM(IMTOFF),SCR,1,LOB,LOB)
        ONEM = -1.0D0
*. scale with -1 to get highest occupation numbers as first eigenvectors
        CALL SCALVE(SCR,ONEM,LOB*(LOB+1)/2)
        CALL EIGEN(SCR,XNAT(IMTOFF),LOB,0,1)
*
        DO  I = 1, LOB
          OCCNUM(IOBOFF-1+I) = - SCR(I*(I+1)/2)
        END DO
*. Order the degenerate eigenvalues so diagonal terms are maximized
        TESTY = 1.0D-11
        DO IOB = 2, LOB
          IF(ABS(OCCNUM(IOBOFF-1+IOB)-OCCNUM(IOBOFF-2+IOB))
     &       .LE.TESTY) THEN
            XII   = ABS(XNAT(IMTOFF-1+(IOB-1)  *LOB+IOB  ))
            XI1I1 = ABS(XNAT(IMTOFF-1+(IOB-1-1)*LOB+IOB-1))
            XII1  = ABS(XNAT(IMTOFF-1+(IOB-1-1)*LOB+IOB  ))
            XI1I  = ABS(XNAT(IMTOFF-1+(IOB-1)  *LOB+IOB-1))
*
            IF( XI1I.GT.XII.AND.XII1.GT.XI1I1 ) THEN
*. interchange orbital IOB and IOB -1
              CALL SWAPVE(XNAT(IMTOFF+(IOB-1)*LOB),
     &                    XNAT(IMTOFF+(IOB-1-1)*LOB),LOB)
              SS = OCCNUM(IOBOFF-1+IOB-1)
              OCCNUM(IOBOFF-1+IOB-1) = OCCNUM(IOBOFF-1+IOB)
              OCCNUM(IOBOFF-1+IOB)   = SS
              write(6,*) ' Orbitals interchanged ',
     &        IOBOFF-1+IOB,IOBOFF-2+IOB
            END IF
          END IF
        END DO
*
        IF(NTEST.GE.1) THEN
          WRITE(6,*)
          WRITE(6,*)
     &    ' Natural occupation numbers for symmetry = ', ISMOB
          WRITE(6,*)
     &    ' ==================================================='
          WRITE(6,*)
          CALL WRTMAT(OCCNUM(IOBOFF),1,LOB,1,LOB)
          IF(NTEST.GE.2 ) THEN
            WRITE(6,*)
            WRITE(6,*) ' Corresponding Eigenvectors '
            WRITE(6,*)
            CALL WRTMAT(XNAT(IMTOFF),LOB,LOB,LOB,LOB)
          END IF
        END IF
      END DO
*. ( End of loop over orbital symmetries )
*
      RETURN
      END
