      SUBROUTINE DIST_2EINTS(NLIST,IBU,ISM,IPRTRA)
*
* Distribute nonrelativistic integrals into
* double group format
*
*---------------------------------------------
* Jeppe Olsen and Timo Fleig, 1997/98
*---------------------------------------------
*
* Distribute DIRAC integrals from symmetry reduced list
* on WORK(KINT2) to double group lists
*
*---------------------------------------------
* Timo Fleig, October 2000
*---------------------------------------------
*
* Input :
*
* NLIST : Number of two electron integral lists
*
* IBU :
*
* IBU(INDEX,ILIST) = 1 => Index INDEX  in list ILIST is   barred
*                  = 2 => Index INDEX  in list ILIST is unbarred
*
* ISM :
*
* ISM(1,ILIST) = 1 => Symmetry between index 1 and 2 in list ILIST
*              = 0 => No Symmetry between index 1 and 2 in list ILIST
*
* ISM(2,ILIST) = 1 => Symmetry between index 3 and 4 nd  in list ILIST
*              = 0 => No Symmetry between index 3 and 4  in list ILIST
*
* ISM(3,ILIST) = 1 => Symmetry between index (12) and (34)  in list ILIST
*              = 0 => No Symmetry between  index (12) and (34) in  ILIST
*
* Note : at the moment all symmetry flags are supposed to be zero !!
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "symm.inc"
#include "cintfo.inc"
#include "glbbas.inc"
*
      INTEGER IBU(4,*), ISM(3,*)
*
      dimension NREOSM(2,MXPNGAS)
*
      NTESTL = 0
      NTEST = max(IPRTRA,NTESTL)
*
      IPNTGRP_SAVE = PNTGRP
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIST2E')
*. Block holding largest block of integrals
      MXTSOB_L = 0
      DO IOBTP = 1, NGAS
        DO IOBSM = 1, NSMOB
          MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IOBTP,IOBSM))
        END DO
      END DO
      INTBLK = MXTSOB_L**4
      CALL MEMMAN(KLINT2,INTBLK,'ADDL  ',2,'KLINT2')
*. Loop over block of normal integrals, fetch them and dispose them
      DO LSYM = 1, NSMOB/2
       DO KSYM = 1, NSMOB/2
        DO JSYM = 1, NSMOB/2
*. Point group symmetry
          PNTGRP = 1
          CALL SYMCOM(3,1,LSYM,KSYM,LKSYM)
          CALL SYMCOM(2,1,LKSYM,JISYM,1)
          CALL SYMCOM(2,1,JSYM,ISYM,JISYM)
*. ISYM, JSYM, KSYM and LSYM are now known
          IF(NTEST.GE.3) THEN
            WRITE(6,*) ' ISYM, JSYM, KSYM, LSYM',
     &      ISYM,JSYM,KSYM,LSYM
          END IF
          PNTGRP = IPNTGRP_SAVE
*. Total number of orbitals of these  symmetries
          NITOT = NGSSHT(ISYM)
          NJTOT = NGSSHT(JSYM)
          NKTOT = NGSSHT(KSYM)
          NLTOT = NGSSHT(LSYM)
          IF(NTEST.GE.3) THEN
            WRITE(6,*) ' NITOT is  ',NITOT
            WRITE(6,*) ' NJTOT is  ',NJTOT
            WRITE(6,*) ' NKTOT is  ',NKTOT
            WRITE(6,*) ' NLTOT is  ',NLTOT
          END IF
*. Corresponding double group symmetries
          ISYM_U = IMOSF_SP(ISYM,1)
          ISYM_B = IMOSF_SP(ISYM,2)
          JSYM_U = IMOSF_SP(JSYM,1)
          JSYM_B = IMOSF_SP(JSYM,2)
          KSYM_U = IMOSF_SP(KSYM,1)
          KSYM_B = IMOSF_SP(KSYM,2)
          LSYM_U = IMOSF_SP(LSYM,1)
          LSYM_B = IMOSF_SP(LSYM,2)
*
          IF(NTEST.GE.3) THEN
            WRITE(6,*) ' ISYM,ISYM_B,ISYM_U', ISYM,ISYM_B,ISYM_U
            WRITE(6,*) ' JSYM,JSYM_B,JSYM_U', JSYM,JSYM_B,JSYM_U
            WRITE(6,*) ' KSYM,KSYM_B,KSYM_U', KSYM,KSYM_B,KSYM_U
            WRITE(6,*) ' LSYM,LSYM_B,LSYM_U', LSYM,LSYM_B,LSYM_U
          END IF
*
          NGSOBCT = 0
*
*
          DO LGAS = 1, NGAS
           NL = NGSSH(LSYM,LGAS)
           DO LLGAS = 1, LGAS
             IF(LLGAS.EQ.1) THEN
               LOFF = 1
             ELSE
               LOFF = LOFF + NGSSH(LSYM,LLGAS)
             END IF
           END DO
*
*
           DO KGAS = 1, NGAS
            NK = NGSSH(KSYM,KGAS)
            DO KKGAS = 1, KGAS
              IF(KKGAS.EQ.1) THEN
                KOFF = 1
              ELSE
                KOFF = KOFF + NGSSH(KSYM,KKGAS)
              END IF
            END DO
*
*
            DO JGAS = 1, NGAS
             NJ = NGSSH(JSYM,JGAS)
             DO JJGAS = 1, JGAS
               IF(JJGAS.EQ.1) THEN
                JOFF = 1
              ELSE
                JOFF = JOFF + NGSSH(JSYM,JJGAS)
              END IF
             END DO
*
*
             DO IGAS = 1, NGAS
              NI = NGSSH(ISYM,IGAS)
              DO IIGAS = 1, IGAS
                IF(IIGAS.EQ.1) THEN
                  IOFF = 1
                ELSE
                  IOFF = IOFF + NGSSH(ISYM,IIGAS)
                END IF
              END DO
*
*
              IGSOBCT = 0
              JGSOBCT = 0
              KGSOBCT = 0
              LGSOBCT = 0
*
              do IIGAS = 1,IGAS-1,1
                 IGSOBCT = IGSOBCT + NGSSH(ISYM,IIGAS)
              end do
*
              do JJGAS = 1,JGAS-1,1
                 JGSOBCT = JGSOBCT + NGSSH(JSYM,JJGAS)
              end do
              JGSOBCT = JGSOBCT * NITOT
*
              do KKGAS = 1,KGAS-1,1
                 KGSOBCT = KGSOBCT + NGSSH(KSYM,KKGAS)
              end do
              KGSOBCT = KGSOBCT * NITOT * NJTOT
*
              do LLGAS = 1,LGAS-1,1
                 LGSOBCT = LGSOBCT + NGSSH(LSYM,LLGAS)
              end do
              LGSOBCT = LGSOBCT * NITOT * NJTOT * NKTOT
*
              NGSOBCT = IGSOBCT + JGSOBCT + KGSOBCT + LGSOBCT
*
*
              if (NI.eq.0.or.NJ.eq.0.or.NK.eq.0.or.NL.eq.0)
     &           goto 666
*
              if (NTEST.ge.3) then
                write(6,'(A23,4I4)') ' IGAS, JGAS, KGAS, LGAS ',
     &          IGAS, JGAS, KGAS, LGAS
                write(6,'(A12,4I5)') ' NI,NJ,NK,NL', NI,NJ,NK,NL
c               write(6,*) 'IGSOBCT,JGSOBCT,KGSOBCT,LGSOBCT ='
c    &                     ,IGSOBCT,JGSOBCT,KGSOBCT,LGSOBCT
c               write(6,*) 'NGSOBCT =  ',NGSOBCT
              end if
*
*. Fetch nonrelativistic or scalar relativistic integrals
*  They are spin-free, that's for sure!
              if (NTEST.ge.2) then
                 write(6,*) 'Input for GETINTN:'
                 write(6,*) 'ISYM =  ',ISYM
                 write(6,*) 'JSYM =  ',JSYM
                 write(6,*) 'KSYM =  ',KSYM
                 write(6,*) 'LSYM =  ',LSYM
              end if
*
              CALL GETINTN(WORK(KLINT2),
     &           IGAS,ISYM,JGAS,JSYM,KGAS,KSYM,LGAS,LSYM,
     &           0,0,0,0,IRI)
C                  GETINTN(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
C    &                    IXCHNG,IKSM,JLSM,ICOUL,IRI)
*
*(IJ!KL) is in WORK(KLINT2) as full (IK!JL) matrix
*
              if (NTEST.ge.2) then
                 write(6,*) 'What has GETINTN fetched?'
                 call wrtmat(WORK(KLINT2),NI*NK*NJ*NL,1,NI*NK*NJ*NL,1)
              end if
CTF
*. Distribute to integral lists
              DO JLIST = 1, NLIST
*. Spin integration : index 1,2 and index 3,4 should be identical
*    (well, actaully rather the barred or unbarred type should be
*     identical ...)
*     This is implemented for later generalization
CTF
                IF(IBU(1,JLIST).EQ.IBU(2,JLIST).AND.
     &             IBU(3,JLIST).EQ.IBU(4,JLIST)) THEN
CTF
* This means that the third list is accessed all the time, as it is
* the only one where all indices have the same IBU value (2), meaning
* all four indices are unbarred, corresponding to NR integrals.
CTF
*. nonvanishing spin integration
*. Offset to given symmetryblock in JLIST
                IF(IBU(1,JLIST).EQ.1) THEN
                 ISYM_D = ISYM_B
                ELSE
                 ISYM_D = ISYM_U
                END IF
*
                IF(IBU(2,JLIST).EQ.1) THEN
                 JSYM_D = JSYM_B
                ELSE
                 JSYM_D = JSYM_U
                END IF
*
                IF(IBU(3,JLIST).EQ.1) THEN
                 KSYM_D = KSYM_B
                ELSE
                 KSYM_D = KSYM_U
                END IF
*
                IF(IBU(4,JLIST).EQ.1) THEN
                 LSYM_D = LSYM_B
                ELSE
                 LSYM_D = LSYM_U
                END IF
*
                if (NTEST.ge.3) then
                   write(6,*) 'ISYM_D =  ',ISYM_D
                   write(6,*) 'JSYM_D =  ',JSYM_D
                   write(6,*) 'KSYM_D =  ',KSYM_D
                   write(6,*) 'LSYM_D =  ',LSYM_D
                end if
*
                IOFF_DBG = IOFF_DBG_LIST(
     &                     JLIST,ISYM_D,JSYM_D,KSYM_D,LSYM_D)
                if (NTEST.ge.3) write(6,*)
     &          'IOFF_DBG is  ',IOFF_DBG
*. And now : Distribute to double group list
                DO L = LOFF, LOFF+NL-1
                 DO K = KOFF, KOFF+NK-1
                  DO J = JOFF, JOFF+NJ-1
                   DO I = IOFF, IOFF+NI-1
*
                    IREL = I - IOFF + 1
                    JREL = J - JOFF + 1
                    KREL = K - KOFF + 1
                    LREL = L - LOFF + 1
*
*. Address in input list (has form (IK!JL))
                    IJKL_IN =
     &               ((LREL-1)*NJ+JREL-1)*NI*NK+(KREL-1)*NI+IREL
*
*. Address in output list (has usual structure (IJ!KL))
                    IJKL_OUT = IOFF_DBG - 1
     &                + ((LREL-1)*NKTOT+KREL-1)*NITOT*NJTOT
     &                + (JREL-1)*NITOT+IREL
     &                + NGSOBCT
*. and now
                    WORK(KINT2D(1,JLIST)-1+IJKL_OUT) =
     &              WORK(KLINT2-1+IJKL_IN)
                    if (NTEST.ge.5) then
                      write(6,*) 'IJKL_IN = ',IJKL_IN
                      write(6,*) 'WORK(KLINT2-1+IJKL_IN) ='
                      write(6,*) '   ',WORK(KLINT2-1+IJKL_IN)
                      write(6,*) 'IJKL_OUT = ',IJKL_OUT
                      write(6,*) 'WORK(KINT2D(1,JLIST)-1+IJKL_OUT) ='
                      write(6,*) '   ',WORK(KINT2D(1,JLIST)-1+IJKL_OUT)
                    end if
                   END DO
                  END DO
                 END DO
                END DO
*               ^ End of loop over orbital of given TS
               END IF
*              ^ End of check for nonvanishing spin integration
              END DO
*             ^End of loop over output lists
 666          continue
             END DO
            END DO
           END DO
          END DO
*          ^ End of loop over gasspaces
        END DO
       END DO
      END DO
*     ^ End of loop over symmetries
*( Oh Yes, I do, I do I do)
*. Vanishing integral lists
      DO JLIST = 1, NLIST
       LENGTH = L2D(JLIST)
       ZERO = 0.0D0
* 1 : all imaginary terms
       IOFF = KINT2D(2,JLIST)
C?     WRITE(6,*) ' JLIST, KINT2D(2,JLIST)',JLIST,  KINT2D(2,JLIST)
       CALL SETVEC(WORK(IOFF),ZERO,LENGTH)
* 2 : zero by spin integration
       IF(.NOT.(IBU(1,JLIST).EQ.IBU(2,JLIST).AND.
     &          IBU(3,JLIST).EQ.IBU(4,JLIST))     ) THEN
C?     WRITE(6,*) 'JLIST,  KINT2D(1,JLIST)',JLIST,  KINT2D(1,JLIST)
         IOFF = KINT2D(1,JLIST)
         CALL SETVEC(WORK(IOFF),ZERO,LENGTH)
       END IF
      END DO
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*)
        WRITE(6,*) ' Integral lists over spinors(real part) '
        WRITE(6,*) ' ---------------------------------------'
        WRITE(6,*) ' Distributed 2-electron integrals       '
        WRITE(6,*)
        DO ILIST = 1, NLIST
          WRITE(6,*)
          LEN_T = L2D(ILIST)
          IOFF = KINT2D(1,ILIST)
          write(6,*) 'List and length :'
          WRITE(6,'(I9)') ILIST,LEN_T
          WRITE(6,*)
          CALL WRTMAT(WORK(IOFF),LEN_T,1,LEN_T,1)
        END DO
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'DIST2E')
*
      RETURN
      END
