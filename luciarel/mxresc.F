      SUBROUTINE MXRESC(IAB,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &                  NSMST,NSTFSMSPGP,
     &                  MXPNSMST,
     &                  NSMOB,MXPTOB,NTPOB,NTSOB,NTESTG,MXPKA,
     &                  NEL1234,
     &                  MXCJ,MXCIJA,
     &                  MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
*
* Find largest dimension of matrix C(Ka,Ib,J)
* Find largest dimension of matrix C(ij,Ka,Ib)
* Find largest dimension of matrix C(ij,Ia,Kb)
* Find largest dimension of matrix C(ij,Ka,Kb)
*
* Largest block of single excitations MXSXBL

*. Input
* IAB :allowed combination of alpha and beta supergroups
* IOCPTA : Number of first active alpha supergroup
* IOCPTB : Number of first active beta  supergroup
* NOCTPA : Number of active alpha supergroups
* NOCTPB : Number of active alpha supergroups

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IAB(NOCTPA,NOCTPB)
      DIMENSION NSTFSMSPGP(MXPNSMST,*)
      DIMENSION NTSOB(MXPTOB,NSMOB)
      DIMENSION NEL1234(MXPTOB,*)
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.100) WRITE(6,*) ' MXRESC : MXPKA ', MXPKA
*
* matrix C(j,Ka,Ib)
*
*. Note : Only done for alpha-strings, problems when transposing
*         constructing C(J,Ia,Kb)
      MXCJ = 0
      MXADKBLK = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(j,Ka,Ib)'
        write(6,*)
      end if
      DO 100 IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO 200 IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100)
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXB = 0
            DO 210 ISM = 1, NSMST
              MXB =MAX(MXB,NSTFSMSPGP(ISM,IBTPABS))
  210       CONTINUE
            IF(NTEST.GE.100) WRITE(6,*) ' MXB = ', MXB
            DO 300 IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,1,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO 310 KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
  310           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                MXKAO = MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &          MXKA= MXPKA
                MXSOB = 0
                DO 320 ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
  320           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKAO)
                LCJBLK = MXSOB*MXKA*MXB
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
  300       CONTINUE
          END IF
  200   CONTINUE
  100 CONTINUE
*
* matrix C(j,Ia,Kb)
*
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(j,Ia,Kb)'
        write(6,*)
      end if
*
      DO IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100)
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXA = 0
            DO ISM = 1, NSMST
              MXA =MAX(MXA,NSTFSMSPGP(ISM,IATPABS))
            END DO
            IF(NTEST.GE.100) WRITE(6,*) ' MXA = ', MXA
            DO IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IBTPABS,1,IOBTP,1,KBTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KBTP ',IOBTP,KBTP
              IF(KBTP.GT.0) THEN
                MXKB = 0
                DO KSM = 1, NSMST
                  MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                MXKBO = MXKB
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &          MXKB= MXPKA
                MXSOB = 0
                DO ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKBO)
                LCJBLK = MXSOB*MXKB*MXB
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
            END DO
          END IF
        END DO
      END DO
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXADKBLK,MXCJ ', MXADKBLK,MXCJ
      END IF
*
* matrix C(ij,Ka,Ib)
* both Ka and Ib blocked
*
      MXCIJA = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(ij,Ka,Ib)'
        write(6,*)
      end if
*
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA -1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIB = 0
            DO  ISM = 1, NSMST
              MXIB = MAX(MXIB,NSTFSMSPGP(ISM,IBTPABS))
            END DO
            IF(MXIB.GT.MXPKA) MXIB = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,1,IOBTP,1,K1ATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1ATP ',IOBTP,K1ATP
              IF(K1ATP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP(K1ATP,1,JOBTP,1,KATP)
                  IF(KATP.GT.0) THEN
                    MXKA = 0
                    DO KSM = 1, NSMST
                      MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                    IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &              MXKA= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKA*MXIB
                    MXCIJA = MAX(MXCIJA,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXCIJA ', MXCIJA
      END IF
*
*
* matrix C(ij,Ia,Kb)
* both Ka and Ib blocked
*
      IF(NTEST.GE.100) WRITE(6,*) ' MXCIJB under development '
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(ij,Ia,Kb)'
        write(6,*)
      end if
*
      MXCIJB = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIA = 0
            DO  ISM = 1, NSMST
              MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS ))
            END DO
            IF(MXIA.GT.MXPKA) MXIA = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IBTP
              CALL NEWTYP(IBTPABS,1,IOBTP,1,K1BTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1BTP ',IOBTP,K1BTP
              IF(K1BTP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP(K1BTP,1,JOBTP,1,KBTP)
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXIA
                    MXCIJB = MAX(MXCIJB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXCIJB ', MXCIJB
      END IF
*
*
* matrix C(ij,Ka,Kb)
* both Ka and Kb blocked
*
      MXCIJAB = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(ij,Ka,Kb)'
        write(6,*)
      end if
*
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            DO  IOBTP = 1, NTPOB
*. type of Ka string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,1,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) MXKA= MXPKA
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from IBTP
                  CALL NEWTYP(IBTPABS,1,JOBTP,1,KBTP)
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXKA
                    MXCIJAB = MAX(MXCIJAB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
*
* Largest block of single excitations :
* Strings of given type and sym, orbitals of given type and sym
*
* Largest block of creations : a+i !kstring> where K string is
* obtained as single annihilations
      MXSXBL = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Largest block of creations a+i |kstring>'
        write(6,*) 'alpha strings'
        write(6,*)
      end if
*
*. For alpha strings :
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        MXIA = 0
        DO  ISM = 1, NSMST
          MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP(IATPABS,1,JOBTP,1,KATP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KATP ',JOBTP,KATP
          IF(KATP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IATPABS))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP(KATP,2,IOBTP,1,JATP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
*
*. For beta  strings :
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Largest block of creations a+i |kstring>'
        write(6,*) 'beta strings'
        write(6,*)
      end if
*
      DO  IBTP = 1, NOCTPB
        IBTPABS = IBTP + IOCTPB - 1
        MXIB = 0
        DO  ISM = 1, NSMST
          MXIB = MAX(MXIB,NSTFSMSPGP(ISM,IBTPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP(IBTPABS,1,JOBTP,1,KBTP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KBTP ',JOBTP,KBTP
          IF(KBTP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IBTP))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP(KBTP,2,IOBTP,1,JBTP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC: MXSXBL : ', MXSXBL
      END IF
*
      RETURN
      END
