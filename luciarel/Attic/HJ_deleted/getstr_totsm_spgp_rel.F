      SUBROUTINE GETSTR_TOTSM_SPGP_REL(ISTRTP,ISPGRP,ISPGRPSM,
     &                                 NEL,NSTR,ISTR,
     &                                 NORBT,IDOREO,IZ,IREO)
*
* Obtain all super-strings of given total symmetry and given
* occupation in each GAS space
*
*.If  IDOREO .NE. 0 THEN reordering array : lexical => actual order is obtained
*
* Nomenclature of the day : superstring : string in complete
*                           orbital space, product of strings in
*                           each GAS space
* =====
* Input
* =====
*
* ISTRTP  : Type of of superstrings ( alpha => 1, beta => 2 )
* ISPGRP :  supergroup number, (relative to start of this type )
* ISPGRPSM : Total symmetry of superstrings
* NEL : Number of electrons
* IZ  : Reverse lexical ordering array for this supergroup
*
*
* ======
* Output
* ======
*
* NSTR : Number of superstrings generated
* ISTR : Occupation of superstring
* IREO : Reorder array ( if IDOREO.NE.0)
*
*
* Jeppe Olsen, July 1995
*
* Double group modifications (april 98 ) : Separate treatment of alpha and beta
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "strbas.inc"
#include "csm.inc"
      INTEGER IZ(NORBT,NEL)
*. output
      INTEGER ISTR(NEL,*), IREO(*)
*. Local scratch
      INTEGER NELFGS(MXPNGAS),ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
*
      common/LOCSYMS/NIRR_DG_LOC
*
      CALL QENTER('GETST')
*
      IF(MOD(ISTRTP,2).EQ.0) THEN
*. Barred spinor
       IUB = 2
      ELSE
       IUB = 1
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================== '
        WRITE(6,*) ' Welcome to GETSTR_TOTSM_SPGP '
        WRITE(6,*) ' ============================== '
        WRITE(6,*)
        WRITE(6,'(A,3I3)')
     & ' Strings to be obtained : Type, supergroup, symmetry ',
     &   ISTRTP,ISPGRP,ISPGRPSM
        WRITE(6,*)
        WRITE(6,*) ' IUB = ', IUB
        IF(IDOREO.NE.0) THEN
          WRITE(6,*) ' IZ matrix '
          CALL IWRTMA(IZ,NORBT,NEL,NORBT,NEL)
        END IF
      END IF
*. Absolut number of this supergroup
      ISPGRPA = IBSPGPFTP(ISTRTP) - 1 + ISPGRP
*. Occupation per gasspace
*. Largest occupied space
      NGASL = 0
*. Largest and lowest symmetries active in each GAS space
      DO IGAS = 1, NGAS
        ITPFGS(IGAS) = ISPGPFTP(IGAS,ISPGRPA)
        NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
        IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
      IF(NGASL.EQ.0) NGASL = 1
*. Number of strings per GAS space and offsets for strings of given sym
      DO IGAS = 1, NGAS
        IF(IUB.EQ.1) THEN
         CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                NNSTSGP(1,IGAS))
         CALL ICOPVE2(WORK(KISTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                IISTSGP(1,IGAS))
        ELSE
         CALL ICOPVE2(WORK(KNSTSGP2(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                NNSTSGP(1,IGAS))
         CALL ICOPVE2(WORK(KISTSGP2(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                IISTSGP(1,IGAS))
        END IF
*
      END DO
*
* Initialize and set MINVAL and MAXVAL
      IZERO = 0
      DO IGAS = 1, NGAS
        MINVAL(IGAS) = IZERO
        MAXVAL(IGAS) = IZERO
        DO ISMST =1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) then
            MAXVAL(IGAS) = ISMST
          end if
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) then
            MINVAL(IGAS) = ISMST
          end if
        END DO
      END DO
      if (NTEST.ge.300) then
        write(6,*) 'NSMST :  ',NSMST
        do IGAS = 1, NGAS
          do ISMST =1, NSMST
            write(6,*) 'NNSTSGP(',ISMST,IGAS,') : ',NNSTSGP(ISMST,IGAS)
          end do
        end do
      end if
*
CTF
* If no strings are found in a given GAS, we set the symmetry
* value in MINVAL to 1, just to prevent crashes (only a test,
* waiting for more elegant solution ...)
      IONE = 1
      do IGAS = 1,NGAS-1,1
        if (MINVAL(IGAS).eq.0) then
          MINVAL(IGAS) = IONE
        end if
      end do
CTF
*
* Largest and lowest active symmetries for each GAS space
      IF(NTEST.GE.200) THEN
         WRITE(6,*) ' Type of each GAS space '
         CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
         WRITE(6,*) ' Number of elecs per GAS space '
         CALL IWRTMA(NELFGS,1,NGAS,1,NGAS)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' =================== '
        WRITE(6,*) ' MINVAL and MAXVAL : '
        WRITE(6,*) ' =================== '
        DO I = 1, NGAS
          WRITE(6,'(2I5)') MINVAL(I),MAXVAL(I)
        END DO
      END IF
*
*. Loop over symmetries of each GAS
*
      MAXLEX = 0
      IFIRST = 1
      ISTRBS = 1
      NGASCT = NGAS
 1000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
C        NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         do IGS=1,NGAS,1
            if (ISMFGS(IGS).gt.NIRR_DG_LOC) then
              if (NTEST.ge.1) then
                write(6,*) 'Problem in getstr_totsm_spgp'
                write(6,*) 'NIRR_DG_LOC = ',NIRR_DG_LOC
                write(6,*) 'ISMFGS out of range.'
                write(6,*) 'ISMFGS(',IGS,') : ',ISMFGS(IGS)
                write(6,*) 'setting symmetry to 1'
              end if
              ISMFGS(IGS) = 1
            end if
         end do
         IF(NONEW.NE.0) GOTO 1001
        END IF
        IFIRST = 0
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL-1 spaces '
          CALL IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
        END IF
*. Symmetry of NGASL -1 spaces given, symmetry of total space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          if (NTEST.ge.300) then
            write(6,*) 'ISMFGS(',IGAS,') in getstr_totsm_spgp',
     &                  ISMFGS(IGAS)
          end if
C?        write(6,*) ' ISTSMM1 : ', ISTSMM1
          CALL SYMCOM_REL(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
C              SYMCOM(ITASK,IOBJ,I1,I2,I12)
          ISTSMM1 = JSTSMM1
        END DO
*. required sym of SPACE NGASL
        if (NTEST.ge.300) then
          write(6,*) 'ISMGSN in getstr_totsm_spgp',ISMGSN
        end if
        CALL SYMCOM_REL(2,1,ISTSMM1,ISMGSN,ISPGRPSM)
        ISMFGS(NGASL) = ISMGSN
*
        DO IGAS = NGASL+1,NGAS
          ISMFGS(IGAS) = 1
        END DO
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' Next symmetry distribution '
          CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        END IF
*. Obtain all strings of this symmetry
CT      CALL QENTER('GASSM')
        CALL GETSTRN_GASSM_SPGP_REL(ISMFGS,ITPFGS,ISTR(1,ISTRBS),
     &                              NSTR,NEL,
     &                              NNSTSGP,IISTSGP,IUB)
CT      CALL QEXIT('GASSM')
*. Reorder Info : Lexical => actual number
        IF(IDOREO.NE.0) THEN
*. Lexical number of NEL electrons
*. Can be made smart by using common factor for first NGAS-1 spaces
          DO JSTR = ISTRBS, ISTRBS+NSTR-1
            LEX = 1
            DO IEL = 1, NEL
              LEX = LEX + IZ(ISTR(IEL,JSTR),IEL)
            END DO
C           WRITE(6,*) ' string '
C           CALL IWRTMA(ISTR(1,JSTR),1,NEL,1,NEL)
C           WRITE(6,*) ' JSTR and LEX ', JSTR,LEX
*
            MAXLEX = MAX(MAXLEX,LEX)
            IREO(LEX) = JSTR
          END DO
        END IF
*
        ISTRBS = ISTRBS + NSTR
*. ready for next symmetry distribution
        if (NTEST.ge.1) then
          write(6,*) 'Loop status in getstr_totsm_spgp'
          write(6,*) 'NGAS-1 is  ',NGAS-1
        end if
        NGASCT = NGASCT - 1
        if (NGAS-1.ne.0) then
c         if (NGASCT.eq.0) goto 1001
          GOTO 1000
        else
          goto 1001
        end if
 1001 CONTINUE
*. End of loop over symmetry distributions
      NSTR = ISTRBS - 1
*. Check for zeroes
       I_DO_ZERO_CHECK = 0
       IF(I_DO_ZERO_CHECK.EQ.1) THEN
         DO JSTR = 1, NSTR
           DO JEL = 1, NEL
             IF(ISTR(JEL,JSTR).EQ.0) THEN
                WRITE(6,*) ' Occupation string with zeroes '
                CALL IWRTMA(ISTR(1,JSTR),1,NEL,1,NEL)
                Call Abend2(       ' Occupation string with zeroes ' )
              END IF
            END DO
         END DO
       END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings generated ', NSTR
        WRITE(6,*)
        WRITE(6,*) ' Strings : '
        WRITE(6,*)
        if (NSTR.ne.0.and.NEL.ne.0) CALL PRTSTR(ISTR,NEL,NSTR)
*
        IF(IDOREO.NE.0) THEN
          WRITE(6,*) 'Largest Lexical number obtained ', MAXLEX
          WRITE(6,*) ' Reorder array '
          CALL IWRTMA(IREO,1,NSTR,1,NSTR)
        END IF
      END IF
*
 1002 continue
*
      CALL QEXIT('GETST')
      RETURN
      END
