      SUBROUTINE ORBINF(LUOUT,IPRNT)
*
* Obtain information about orbitals from shell information
*
* =====
* input
* =====
* Shell and symmetry information in /LUCINP/
*
* ======
* Output
* ======
* Orbital information in /ORBINP/
*
* Jeppe Olsen, Winter of 1991
*              Updated, July 97 : Double group symmetry added

*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER CITYP
#include "mxpdim.inc"
#include "lucinp.inc"
#include "cgas.inc"
*
#include "orbinp.inc"
#include "symm.inc"
#include "csm.inc"
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
*. Double group
      IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
        IDBG = 1
*. number of symmetries in point group sym
        NSMSH = NSMOB/2
      ELSE
        IDBG = 0
        NSMSH = NSMOB
      END IF
*. Total number of shells per symmetry
      DO ISHSM = 1, NSMSH
        LORB = 0
        DO JGAS = 1, NGAS
          LORB = LORB + NGSSH(ISHSM,JGAS)
        END DO
        NGSSHT(ISHSM) = LORB
      END DO
*
CTF   DO ISMSH = 1, NSMSH
CTF     IF(ISMSH.EQ.1) THEN
CTF      IBSH(ISMSH) = 1
CTF     ELSE
CTF      IBSH(ISMSH) = IBSH(ISMSH-1)+NGSSHT(ISMSH)
CTF     END IF
CTF   END DO
* There is something fishy here; try this:
      IBSH(1) = 1
      do ISMSH = 2,NSMSH,1
         IVAL = 1
         do I = 1,ISMSH-1,1
            IVAL = IVAL + NGSSHT(I)
         end do
         IBSH(ISMSH) = IVAL
      end do
*
C?    WRITE(6,*) ' IBSH in ORBINP '
C?    CALL IWRTMA(IBSH,1,NSMSH,1,NSMSH)

*
* Point group => double group information if required
* already done. Info on symm.inc
      call grinfo(PNTGRP)
*                                              *
* 1 : From shell format to orbital format *
*
      CALL OSPIR(NOSPIR,IOSPIR,PNTGRP,NIRREP,IPRNT,
     &           NORB_SET)
*
* 2 : Shell information to orbital information for each group of orbital
*
*
* ===============
*     GAS case
* ===============
*
      DO IGAS = 1, NGAS
*. Shell => orbitals for each GAS space
        CALL SHTOOB(NGSSH(1,IGAS),NIRREP,MXPOBS,NSMOB,NOSPIR,
     &              IMOSF_SP(1,1),NGSOB(1,IGAS),NGSOBT(IGAS) )
      END DO
*
      IF(IDBG.EQ.1) THEN
*. A second orbital set must be generated corresponding to beta
* spinors. The corresponding orbital information is stored in
* NGSOBT2 NGSOBT2
*. Orbital => spinors for beta spintors
C?      WRITE(6,*) ' Generate beta spinors '
C?      WRITE(6,*) ' NIRR_DG ', NIRR_DG
C       DO IRREP = 1, NIRR_DG
C         IOSPIR(1,IRREP) = IOSPIR(2,IRREP)
C       END DO
        DO IGAS = 1, NGAS
*. Shell => orbitals for each GAS space
          CALL SHTOOB(NGSSH(1,IGAS),NIRREP,MXPOBS,NSMOB,NOSPIR,
     &         IMOSF_SP(1,2),NGSOB2(1,IGAS),NGSOBT2(IGAS))
C    &                IOSPIR,NGSOB2(1,IGAS),NGSOBT2(IGAS))
        END DO
      END IF
*
*  ========================================================
*. Number of inactive, active, occupied , deleted orbitals
*  ========================================================
*
*
* current inactive and deleted orbitals are not identified so
      IGSINA = 0
      IGSDEL = 0
*
      CALL ISETVC(NTOOBS,0,NSMOB)
      CALL ISETVC(NOCOBS,0,NSMOB)
      CALL ISETVC(NACOBS,0,NSMOB)
*
      NTOOB = 0
      NACOB = 0
      NOCOB = 0
      DO IGAS = 1, NGAS
*. Inactive orbitals
        IF(IGAS.EQ.IGSINA) THEN
          CALL ICOPVE(NGSOB(1,IGAS),NINOBS,NSMOB)
          NINOB = NGSOBT(IGAS)
        END IF
*. Deleted orbitals
        IF(IGAS.EQ.IGSDEL) THEN
          CALL ICOPVE(NGSOB(1,IGAS),NDEOBS,NSMOB)
          NDEOB = NGSOBT(IGAS)
        END IF
*. Add to total number of orbitals
        CALL IVCSUM(NTOOBS,NTOOBS,NGSOB(1,IGAS),1,1,NSMOB)
        NTOOB = NTOOB + NGSOBT(IGAS)
*. Add to occupied orbitals
        IF(IGAS.NE.IGSDEL) THEN
          CALL IVCSUM(NOCOBS,NOCOBS,NGSOB(1,IGAS),1,1,NSMOB)
          NOCOB = NOCOB + NGSOBT(IGAS)
        END IF
*. Add to active orbitals
        IF(IGAS.NE.IGSINA.AND.IGAS.NE.IGSDEL) THEN
          CALL IVCSUM(NACOBS,NACOBS,NGSOB(1,IGAS),1,1,NSMOB)
          NACOB = NACOB + NGSOBT(IGAS)
        END IF
      END DO
* ===============================================
*. Well, report back
* ===============================================
      IF(NTEST.GT.0) THEN
        WRITE(LUOUT,*)
        IF(IDBG.EQ.0) THEN
          WRITE(LUOUT,*) ' Number of orbitals per symmetry :'
          WRITE(LUOUT,*) ' ================================='
        ELSE
          WRITE(LUOUT,*) ' Number of up spinors per symmetry :'
          WRITE(LUOUT,*) ' =================================='
        END IF
        WRITE(LUOUT,*)
        WRITE(LUOUT,'(1H ,A,10I4,A)')
     &  '            Symmetry  ',(I,I = 1,NSMOB)
        WRITE(LUOUT,'(1H ,A,2X,10A)')
     &  '           ========== ',('====',I = 1,NSMOB)
        DO IGAS = 1, NGAS
          WRITE(LUOUT,'(1H      ,A,I3,7X,A,10I4,8X,I3)')
     &    '   GAS',IGAS,'      ',(NGSOB(I,IGAS),I=1,NSMOB),
     &    NGSOBT(IGAS)
        END DO
*
        WRITE(LUOUT,*) ' Total number of orbitals ', NTOOB
        WRITE(LUOUT,*) ' Total number of occupied orbitals ', NOCOB
      END IF
*. Offsets for orbitals of given symmetry
      ITOOBS(1) = 1
      DO  ISMOB = 2, NSMOB
        ITOOBS(ISMOB) = ITOOBS(ISMOB-1)+NTOOBS(ISMOB-1)
      END DO
*
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Offsets for orbital of given symmetry '
        CALL IWRTMA(ITOOBS,1,NSMOB,1,NSMOB)
      END IF
*
      IF(IDBG.EQ.1) THEN
*. Repete the above for construction of beta spinors
* current inactive and deleted orbitals are not identified so
       IGSINA = 0
       IGSDEL = 0
*
       CALL ISETVC(NTOOBS2,0,NSMOB)
       CALL ISETVC(NOCOBS2,0,NSMOB)
       CALL ISETVC(NACOBS2,0,NSMOB)
*
       NTOOB2 = 0
       NACOB2 = 0
       NOCOB2 = 0
       DO IGAS = 1, NGAS
*. Inactive orbitals
         IF(IGAS.EQ.IGSINA) THEN
           CALL ICOPVE(NGSOB2(1,IGAS),NINOBS2,NSMOB)
           NINOB2 = NGSOBT2(IGAS)
         END IF
*. Deleted orbitals
         IF(IGAS.EQ.IGSDEL) THEN
           CALL ICOPVE(NGSOB2(1,IGAS),NDEOBS2,NSMOB)
           NDEOB2 = NGSOBT(IGAS)
         END IF
*. Add to total number of orbitals
         CALL IVCSUM(NTOOBS2,NTOOBS2,NGSOB2(1,IGAS),1,1,NSMOB)
         NTOOB2 = NTOOB2 + NGSOBT2(IGAS)
*. Add to occupied orbitals
         IF(IGAS.NE.IGSDEL) THEN
           CALL IVCSUM(NOCOBS2,NOCOBS2,NGSOB2(1,IGAS),1,1,NSMOB)
           NOCOB2 = NOCOB2 + NGSOBT2(IGAS)
         END IF
*. Add to active orbitals
         IF(IGAS.NE.IGSINA.AND.IGAS.NE.IGSDEL) THEN
           CALL IVCSUM(NACOBS2,NACOBS2,NGSOB2(1,IGAS),1,1,NSMOB)
           NACOB2 = NACOB2 + NGSOBT2(IGAS)
         END IF
       END DO
* ===============================================
*. Well, report back
* ===============================================
       IF(NTEST.GT.0) THEN
         WRITE(LUOUT,*)
         WRITE(LUOUT,*) ' Number of down spinors per symmetry :'
         WRITE(LUOUT,*) ' ===================================='
         WRITE(LUOUT,*)
         WRITE(LUOUT,'(1H ,A,10I4,A)')
     &   '            Symmetry  ',(I,I = 1,NSMOB)
         WRITE(LUOUT,'(1H ,A,2X,10A)')
     &   '           ========== ',('====',I = 1,NSMOB)
         DO IGAS = 1, NGAS
           WRITE(LUOUT,'(1H      ,A,I3,7X,A,10I4,8X,I3)')
     &     '   GAS',IGAS,'      ',(NGSOB2(I,IGAS),I=1,NSMOB),
     &     NGSOBT2(IGAS)
         END DO
         WRITE(LUOUT,*) ' Total number of orbitals ', NTOOB2
         WRITE(LUOUT,*) ' Total number of occupied orbitals ', NOCOB2
       END IF
*. Offsets for orbitals of given symmetry
       ITOOBS2(1) = 1
       DO  ISMOB = 2, NSMOB
         ITOOBS2(ISMOB) = ITOOBS2(ISMOB-1)+NTOOBS2(ISMOB-1)
       END DO
*
       IF(NTEST.GT.0) THEN
         WRITE(6,*) ' Offsets for orbital of given symmetry '
         CALL IWRTMA(ITOOBS2,1,NSMOB,1,NSMOB)
       END IF
      END IF
*.Information about spinor symmetries
      write(6,*)
        WRITE(LUOUT,'(A)')
     &  '  Spin-up Spinors in Double group symmetry  '
        WRITE(LUOUT,'(A)')
     &  ' =========================================== '
        WRITE(LUOUT,*)
        WRITE(LUOUT,'(1H ,A,10I4,A)')
     &    '                Irrep ',(I,I = 1,NSMOB)
        WRITE(LUOUT,'(1H ,A,2X,10A,A)')
     &    '                ===== ',('====',I = 1,NSMOB )
        DO IGAS = 1, NGAS
          WRITE(LUOUT,'(A,I2,A,10I4,6X,2I6)')
     &    '        GAS',IGAS,'          ',
     &    (NGSOB(IRREP,IGAS),IRREP = 1, NSMOB )
        END DO
*
        WRITE(LUOUT,'(A)')
     &  '  Spin-down Spinors in Double group symmetry  '
        WRITE(LUOUT,'(A)')
     &  ' =========================================== '
        WRITE(LUOUT,*)
        WRITE(LUOUT,'(1H ,A,10I4,A)')
     &    '                Irrep ',(I,I = 1,NSMOB)
        WRITE(LUOUT,'(1H ,A,2X,10A,A)')
     &    '                ===== ',('====',I = 1,NSMOB )
        DO IGAS = 1, NGAS
          WRITE(LUOUT,'(A,I2,A,10I4,6X,2I6)')
     &    '        GAS',IGAS,'          ',
     &    (NGSOB2(IRREP,IGAS),IRREP = 1, NSMOB )
        END DO

********************************************
*                                          *
* Part 2 : Reordering arrays for orbitals  *
*                                          *
********************************************
        CALL ORBORD_GAS(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB,NGSOBT,
     &       NOCOBS,NTOOBS,NTOOB,
     &       IREOST,IREOTS,ISMFTO,ITPFSO,
     &       IBSO,NTSOB,IBTSOB,ITSOB,NOBPTS,IOBPTS,
     &       ISMFSO,ITPFTO,NOBPT,IPRNT,IOSPIR,1)
      IF(IDBG.EQ.1) THEN
*. Reorder array for down spinors
        CALL ORBORD_GAS(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB2,NGSOBT2,
     &       NOCOBS2,NTOOBS2,NTOOB2,
     &       IREOST2,IREOTS2,ISMFTO2,ITPFSO2,
     &       IBSO2,NTSOB2,IBTSOB2,ITSOB2,NOBPTS2,IOBPTS2,
     &       ISMFSO2,ITPFTO2,NOBPT,IPRNT,IOSPIR,2)
        END IF
*
      RETURN
      END
