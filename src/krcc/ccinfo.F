      SUBROUTINE ZIRAT_KRCC
*
* Ratio between real and integer
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
#include "irat.inc"
*
      IRAT = 1
*
      WRITE(6,*)
      WRITE(6,'(1H ,6X,A,I2)')
     &'Ratio between Integer and Real word length ', IRAT
      WRITE(6,*)
*
      RETURN
      END
*
***********************************************************************
      SUBROUTINE ORBINF_KRCC(LUOUT)
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
* Much in this routine should be deleted!! Lasse 2011
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
*. Well we first need to get the info from dirac common blocks
      CALL SYM_INF_DIRAC
  
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
*.
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
* !!!!May not be used?!!!!
      call grinfo_KRCC(PNTGRP)
*                                              *
* 1 : From shell format to orbital format *
*
      CALL OSPIR_KRCC(NOSPIR,IOSPIR,PNTGRP,NIRREP,
     &               NORB_SET)
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
        CALL SHTOOB_KRCC(NGSSH(1,IGAS),NIRREP,MXPOBS,NSMOB,NOSPIR,
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
          CALL SHTOOB_KRCC(NGSSH(1,IGAS),NIRREP,MXPOBS,NSMOB,NOSPIR,
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
     &    '   GAS',IGAS,'      ',(NGSOB(I,IGAS),I=1,NSMOB)
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
     &     '   GAS',IGAS,'      ',(NGSOB2(I,IGAS),I=1,NSMOB)
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
      if (NTEST.ge.1) then
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
      end if
*
********************************************
*                                          *
* Part 2 : Reordering arrays for orbitals  *
*                                          *
********************************************
        CALL ORBORD_GAS_KRCC(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB,NGSOBT,
     &       NOCOBS,NTOOBS,NTOOB,
     &       IREOST,IREOTS,ISMFTO,ITPFSO,
     &       IBSO,NTSOB,IBTSOB,ITSOB,NOBPTS,IOBPTS,
     &       ISMFSO,ITPFTO,NOBPT,IOSPIR,1)
      IF(IDBG.EQ.1) THEN
*. Reorder array for down spinors
        CALL ORBORD_GAS_KRCC(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB2,NGSOBT2,
     &       NOCOBS2,NTOOBS2,NTOOB2,
     &       IREOST2,IREOTS2,ISMFTO2,ITPFSO2,
     &       IBSO2,NTSOB2,IBTSOB2,ITSOB2,NOBPTS2,IOBPTS2,
     &       ISMFSO2,ITPFTO2,NOBPT,IOSPIR,2)
        END IF
*
      RETURN
      END
*
*  This subroutine generates all group info needed for
*   - integral transformation
*   - CI
*   - Spinor rotation
*
      subroutine grinfo_KRCC(PNTGRP)
      implicit real*8 (A-H,O-Z)
*
      integer PNTGRP
*
#include "mxpdim.inc"
#include "symm.inc"
*
      if (PNTGRP.eq.9) DOUGRP = 4
      if (PNTGRP.eq.8) DOUGRP = 5
      if (PNTGRP.eq.7) DOUGRP = 6
      if (PNTGRP.eq.6) DOUGRP = 7
      if (PNTGRP.eq.5) DOUGRP = 8
*
      if (DOUGRP.eq.1) then
         NIRR_DG = 10
         NIRR_PN = 8
         IHASDG = 4
         PNTGRP = 1
      else if (DOUGRP.eq.2) then
         NIRR_DG = 5
         NIRR_PN = 4
         IHASDG = 5
         PNTGRP = 0
         write(6,*) 'Point group not defined in LUCIA '
         write(6,*) 'watch out!'
      else if (DOUGRP.eq.3) then
         NIRR_DG = 5
         NIRR_PN = 4
         IHASDG = 5
         PNTGRP = 0
         write(6,*) 'Point group not defined in LUCIA '
         write(6,*) 'watch out!'
      else if (DOUGRP.eq.4) then
         NIRR_DG = 8
         NIRR_PN = 4
         IHASDG = DOUGRP
         PNTGRP = 9
      else if (DOUGRP.eq.5.or.DOUGRP.eq.6) then
         NIRR_DG = 4
         NIRR_PN = 2
         IHASDG = DOUGRP
         if (DOUGRP.eq.5) PNTGRP = 8
         if (DOUGRP.eq.6) PNTGRP = 7
      else if (DOUGRP.eq.7) then
         NIRR_DG = 4
         NIRR_PN = 2
         IHASDG = DOUGRP
         PNTGRP = 6
      else if (DOUGRP.eq.8) then
         NIRR_DG = 2
         NIRR_PN = 1
         IHASDG = DOUGRP
         PNTGRP = 5
      else
         write(*,*) 'Double group out of range!'
         write(*,*) 'Double group is  ',DOUGRP
         write(*,*) 'Picasso will put down his brush.'
         Call Abend()
      end if
*
*  Set up info for highest abelian double subgroup
*
      if (DOUGRP.ne.IHASDG) then
         if (DOUGRP.eq.1) then
            NSDGIRR = 8
            NSPGIRR = 4
         else if (DOUGRP.eq.2.or.DOUGRP.eq.3) then
            NSDGIRR = 4
            NSPGIRR = 2
         end if
      end if
*
      MXPDBGDIM = 8
      MXPDBGDIM2 = MXPDBGDIM ** 2
*
      IZERO = 0
      CALL ISETVC(IDBGMULT,IZERO,MXPDBGDIM2)
      CALL ISETVC(IMOSF_SP,IZERO,MXPDBGDIM)
      CALL ISETVC(ISPSF_MO,IZERO,MXPDBGDIM)
      IF(PNTGRP.EQ.5) THEN
*. C1 point group
        DBGSYM(1:6) ='C1    '
*
        IPGMULT(1,1) = 1
        IDBGMULT(1,1) = 1
        IDBGMULT(1,2) = 2
        IDBGMULT(2,1) = 2
        IDBGMULT(2,2) = 1
*. All spinors have symmetry 2, all orbitals symmetry 1
        IMOSF_SP(1,1) = 2
        IMOSF_SP(1,2) = 2
*
        ISPSF_MO(2,1) = 1
        ISPSF_MO(2,2) = 1
*. Spinors and their cc have identical symmetries
        IADJSYM(1) = 1
        IADJSYM(2) = 2
* inverse elements
        INVELM(1) = 1
        INVELM(2) = 2
      ELSE IF(PNTGRP.EQ.6 ) THEN
*. Ci pointgroup
        DBGSYM(1:6) ='Ci    '
        IPGMULT(1,1) = 1
        IPGMULT(1,2) = 2
        IPGMULT(2,1) = 2
        IPGMULT(2,2) = 1
*
        IDBGMULT(1,1) = 1
        IDBGMULT(1,2) = 2
        IDBGMULT(1,3) = 3
        IDBGMULT(1,4) = 4
        IDBGMULT(2,1) = 2
        IDBGMULT(2,2) = 1
        IDBGMULT(2,3) = 4
        IDBGMULT(2,4) = 3
        IDBGMULT(3,1) = 3
        IDBGMULT(3,2) = 4
        IDBGMULT(3,3) = 1
        IDBGMULT(3,4) = 2
        IDBGMULT(4,1) = 4
        IDBGMULT(4,2) = 3
        IDBGMULT(4,3) = 2
        IDBGMULT(4,4) = 1
*. MO have symmetries 1,2 spinors symmetries 3,4 alpha,beta have symm 3
        IMOSF_SP(1,1) = 3
        IMOSF_SP(1,2) = 3
        IMOSF_SP(2,1) = 4
        IMOSF_SP(2,2) = 4
*
        ISPSF_MO(3,1) = 1
        ISPSF_MO(3,2) = 1
        ISPSF_MO(4,1) = 2
        ISPSF_MO(4,2) = 2
*. Functions and adjoints have identical symmetries
        IADJSYM(1) = 1
        IADJSYM(2) = 2
        IADJSYM(3) = 3
        IADJSYM(4) = 4
* inverse elements
        INVELM(1) = 1
        INVELM(2) = 2
        INVELM(3) = 3
        INVELM(4) = 4
      ELSE IF(PNTGRP.EQ.7.OR.PNTGRP.EQ.8) THEN
*. CS and C2 point groups
        IF(PNTGRP.EQ.7) THEN
          DBGSYM(1:6) ='Cs    '
        ELSE
          DBGSYM(1:6) ='C2    '
        END IF
*
        IPGMULT(1,1) = 1
        IPGMULT(1,2) = 2
        IPGMULT(2,1) = 2
        IPGMULT(2,2) = 1
*
        IDBGMULT(1,1) = 1
        IDBGMULT(1,2) = 2
        IDBGMULT(1,3) = 3
        IDBGMULT(1,4) = 4
        IDBGMULT(2,1) = 2
        IDBGMULT(2,2) = 1
        IDBGMULT(2,3) = 4
        IDBGMULT(2,4) = 3
        IDBGMULT(3,1) = 3
        IDBGMULT(3,2) = 4
        IDBGMULT(3,3) = 2
        IDBGMULT(3,4) = 1
        IDBGMULT(4,1) = 4
        IDBGMULT(4,2) = 3
        IDBGMULT(4,3) = 1
        IDBGMULT(4,4) = 2
*. Orbitals have sym 1,2, spinors 3,4, alpha 3, beta 4
        IMOSF_SP(1,1) = 3
        IMOSF_SP(1,2) = 4
        IMOSF_SP(2,1) = 4
        IMOSF_SP(2,2) = 3
*
        ISPSF_MO(3,1) = 1
        ISPSF_MO(3,2) = 2
        ISPSF_MO(4,1) = 2
        ISPSF_MO(4,2) = 1
*. Adjoints
        IADJSYM(1) = 1
        IADJSYM(2) = 2
        IADJSYM(3) = 4
        IADJSYM(4) = 3
* inverse elements
        INVELM(1) = 1
        INVELM(2) = 2
        INVELM(3) = 4
        INVELM(4) = 3
      ELSE IF(PNTGRP.EQ.9) THEN
*. C2H
        DBGSYM(1:6) ='C2h   '
*
        IADJSYM(1) = 1
        IADJSYM(2) = 2
        IADJSYM(3) = 3
        IADJSYM(4) = 4
        IADJSYM(5) = 6
        IADJSYM(6) = 5
        IADJSYM(7) = 8
        IADJSYM(8) = 7
*
        INVELM(1) = 1
        INVELM(2) = 2
        INVELM(3) = 3
        INVELM(4) = 4
        INVELM(5) = 6
        INVELM(6) = 5
        INVELM(7) = 8
        INVELM(8) = 7
*
        ISPSF_MO(5,1) = 1
        ISPSF_MO(6,2) = 1
        ISPSF_MO(6,1) = 2
        ISPSF_MO(5,2) = 2
        ISPSF_MO(7,1) = 3
        ISPSF_MO(8,2) = 3
        ISPSF_MO(8,1) = 4
        ISPSF_MO(7,2) = 4
*
        IMOSF_SP(1,1) = 5
        IMOSF_SP(1,2) = 6
        IMOSF_SP(2,1) = 6
        IMOSF_SP(2,2) = 5
        IMOSF_SP(3,1) = 7
        IMOSF_SP(3,2) = 8
        IMOSF_SP(4,1) = 8
        IMOSF_SP(4,2) = 7
*
        do JJ=1,8,1
           do II=1,8,1
              if (II.le.4.and.JJ.le.4) then
                 if (II.eq.JJ) then
                    IDBGMULT(II,JJ) = 1
                 else if ((II+JJ).eq.3.or.(II+JJ).eq.7) then
                    IDBGMULT(II,JJ) = 2
                 else if ((II+JJ).eq.5) then
                    IDBGMULT(II,JJ) = 4
                 else
                    IDBGMULT(II,JJ) = 3
                 end if
              else if (II.gt.4.and.JJ.gt.4) then
                 if (II.eq.JJ) then
                    IDBGMULT(II,JJ) = 2
                 else if ((II+JJ).eq.13) then
                    IDBGMULT(II,JJ) = 3
                 else if ((II+JJ).eq.11.or.(II+JJ).eq.15) then
                    IDBGMULT(II,JJ) = 1
                 else
                    IDBGMULT(II,JJ) = 4
                 end if
              else
                 if ((II+JJ).eq.9) then
                    IDBGMULT(II,JJ) = 8
                 else if ((II+JJ).eq.7.or.(II+JJ).eq.11) then
                    IDBGMULT(II,JJ) = 6
                 else if ((II+JJ).eq.6.or.(II+JJ).eq.12) then
                    IDBGMULT(II,JJ) = 5
                 else if (II.eq.6.and.JJ.eq.2) then
                    IDBGMULT(II,JJ) = 5
                 else if (II.eq.2.and.JJ.eq.6) then
                    IDBGMULT(II,JJ) = 5
                 else if (II.eq.3.and.JJ.eq.7) then
                    IDBGMULT(II,JJ) = 5
                 else if (II.eq.7.and.JJ.eq.3) then
                    IDBGMULT(II,JJ) = 5
                 else
                    IDBGMULT(II,JJ) = 7
                 end if
              end if
           end do
        end do
        do JJ=1,4,1
           do II=1,4,1
              if (II.eq.JJ) then
                 IPGMULT(II,JJ) = 1
              else if ((II+JJ).eq.3.or.(II+JJ).eq.7) then
                 IPGMULT(II,JJ) = 2
              else if ((II+JJ).eq.5) then
                 IPGMULT(II,JJ) = 4
              else
                 IPGMULT(II,JJ) = 3
              end if
           end do
        end do
        continue
      END IF
*
      NTEST = 00
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Symmetry Information from GET_DBGSYM'
        WRITE(6,*) ' ===================================='
        WRITE(6,*)
        WRITE(6,'(A,A)') ' Double group : ', DBGSYM
        WRITE(6,*)
        WRITE(6,'(A,I2)') ' Number of irreps ',NIRR_DG
        WRITE(6,*)
        WRITE(6,*) ' Multiplication table '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        DO IREP = 1, NIRR_DG
          WRITE(6,'(10X,8I3)') (IDBGMULT(IREP,JREP),JREP=1,NIRR_DG)
        END DO
*
        WRITE(6,*)
        WRITE(6,*) ' orbital X spin => spinor '
        WRITE(6,*) ' ========================='
        WRITE(6,*)
        DO IOBSM = 1, NIRR_DG/2
         DO ISF = 1, 2
          WRITE(6,'(3X,I3,3X,I3,5X,I3)')
     &    IOBSM, ISF, IMOSF_SP(IOBSM,ISF)
         END DO
        END DO
*
        WRITE(6,*)
        WRITE(6,*) ' Spinor  X spin => orbital'
        WRITE(6,*) ' ========================='
        WRITE(6,*)
        DO ISPSM =  NIRR_DG/2+1, NIRR_DG
         DO ISF = 1, 2
          WRITE(6,'(3X,I3,3X,I3,5X,I3)')
     &    ISPSM, ISF, ISPSF_MO(ISPSM,ISF)
         END DO
        END DO
*
        WRITE(6,*)
        WRITE(6,*) ' Adjoint symmetry array '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        WRITE(6,'(10X,8I3)') (IADJSYM(IREP),IREP=1,NIRR_DG)
*
        WRITE(6,*) ' Inverse elements '
        WRITE(6,*) ' ================='
        WRITE(6,'(10X,8I3)') (INVELM(IREP),IREP=1,NIRR_DG)
      END IF
*
      return
      end
*
      SUBROUTINE OSPIR_KRCC(NOSPIR,IOSPIR,PNTGRP,NIRREP,
     &                     NORB_SET)
*
* Number and symmetries of orbitals corresponding to a given shell
*
* =====
* Input
* =====
*
*   PNTGRP  : type of pointgroup
*         = 1 => D2h or a subgroup of D2H
*         = 2 => C inf v
*         = 3 => D inf h
*         = 4 => O 3
*         = 5 to 9 => different double groups :
*             9 :   C2H
*             8 :   C2
*             7 :   CS
*             6 :   CI
*             5 :   C1
*   NIRREP : Number of irreducible representations per point group
*   MXPIRR : Largest allowed number of shell irreps
*   MXPOBS : Largest allowed number of orbital symmetries
*
* ======
* Output
* ======
*
*   NOSPIR : Number of orbital symmetries per irrep
*   IOSPIR : Orbital symmetries corresponding to a given irrep
*
* Jeppe Olsen , Winter of 1991
*              Updated, July 97 : Double group symmetry added
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
      INTEGER PNTGRP
*. Double group info
#include "symm.inc"
*. Output
      DIMENSION NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR)
*
*. In general -except for double groups we have just a single
*  orbitals set
      NORB_SET = 1
      IF(PNTGRP.GE.5) NORB_SET = 2
      if (PNTGRP.eq.1) then
*=====
*.D2h
*=====
        NSMOB = 0
        DO 10 IRREP = 1, 8
          NOSPIR(IRREP) = 1
          IOSPIR(1,IRREP) = IRREP
   10   CONTINUE
      ELSE IF(PNTGRP.EQ.2) THEN
* =========
*. C inf V
* =========
* orbital symmetry is numbered as IML - MNMLOB + 1
        MNMLOB = -(NIRREP-1)
        DO 20 IRREP = 1, NIRREP
*.Irrep I contains orbitals with ML = -(IRREP-1),+(IRREP-1)
          IF(IRREP.EQ.1) THEN
            NOSPIR(IRREP) = 1
            IOSPIR(1,IRREP) = IRREP - 1 - MNMLOB + 1
          ELSE
            NOSPIR(IRREP) = 2
            IOSPIR(1,IRREP) = -(IRREP - 1) - MNMLOB + 1
            IOSPIR(2,IRREP) =  (IRREP - 1) - MNMLOB + 1
          END IF
   20   CONTINUE
      ELSE IF(PNTGRP.EQ.3) THEN
* ========
*. D inf H
* ========
* orbital symmetry is numbered as (PARITY-1) * NMLOB + IML - MNMLOB + 1
        MXMLOB =  NIRREP/2-1
        MNMLOB = -MXMLOB
        NMLOB =   NIRREP - 1
        IRREP = 0
        DO 35 IPARI = 1, 2
          IADD = (IPARI-1)*NMLOB
          DO 30 ML = 0,MXMLOB
            IRREP = IRREP + 1
            IF(ML.EQ.0) THEN
              NOSPIR(IRREP) = 1
              IOSPIR(1,IRREP) = IADD + ML - MNMLOB + 1
            ELSE
              NOSPIR(IRREP) = 2
              IOSPIR(1,IRREP) = IADD - ML - MNMLOB + 1
              IOSPIR(2,IRREP) = IADD + ML - MNMLOB + 1
            END IF
   30     CONTINUE
   35   CONTINUE

      ELSE IF(PNTGRP.EQ.4) THEN
* =====
*. O 3
* =====
* orbital symmetry is numbered as (PARITY-1) * NMLOB + IML - MNMLOB + 1
        MXMLOB =  NIRREP/2-1
        MNMLOB = -MXMLOB
        NMLOB =   NIRREP - 1
        DO 45 L = 0, NIRREP - 1
          IF(MOD(L,2).EQ.0) THEN
            IPARI = 1
          ELSE
            IPARI = 2
          END IF
          IADD = (IPARI-1)*NMLOB
          IRREP = L + 1
          NOSPIR(IRREP) = 2 * L + 1
          ICOMP = 0
          DO 40 ML = MNMLOB,MXMLOB
            ICOMP = ICOMP + 1
            IOSPIR(ICOMP,IRREP) = IADD + ML - MNMLOB + 1
   40     CONTINUE
   45   CONTINUE
      ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
* ========================
*. Double group symmetry
* ========================
*
* We are handling the transformation from orbitals to spinors
* by having two sets of orbitals. This turns out to be
* significant more convenient ( we hope ) than having
* a single set with double dimension
        NORB_SET = 2
        DO  IRREP = 1, NIRR_DG
          NOSPIR(IRREP) = 1
        END DO
*. IOSPIR become the symmetry of the spinors
        DO  IRREP = 1, NIRR_DG
*. The first set corresponds to alpha spin functions
          IOSPIR(1,IRREP) = IMOSF_SP(IRREP,1)
*. The second set corresponds to beta spin functions
          IOSPIR(2,IRREP) = IMOSF_SP(IRREP,2)
        END DO
*. The first NSMOB/2 symmetries are now the fermion irreps, whereas
*  the boson irreps are missing. Add these for completeness
        DO ISMOB = 1, NIRR_PN
         IOSPIR(1,NIRR_PN+ISMOB) = ISMOB
         IOSPIR(2,NIRR_PN+ISMOB) = ISMOB
        END DO
      ELSE
        WRITE(6,*) ' Sorry  PNTGRP out of range , PNTGRP = ', PNTGRP
        WRITE(6,*) ' OSPIR fatally wounded '
        Call Abend1( 5 )
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' OSPIR speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' Number of orbital sets ',NORB_SET
        WRITE(6,*) ' Number of orbitals per irrep '
*
        IF(.NOT.(PNTGRP.GE.5.AND.PNTGRP.LE.9)) THEN
*. Info for orbitals
          CALL IWRTMA(NOSPIR,1,NIRREP,1,NIRREP)
          WRITE(6,*) ' Orbital symmetries per irrep '
          DO 100 IRREP = 1, NIRREP
          CALL IWRTMA(IOSPIR(1,IRREP),1,NOSPIR(IRREP),1,NOSPIR(IRREP))
  100     CONTINUE
        ELSE
*. Info for spinors
          WRITE(6,*) ' IOSPIR for alpha spinors '
          WRITE(6,'(10I4)') (IOSPIR(1,I),I=1,NIRR_DG)
          WRITE(6,*) ' IOSPIR for beta spinors '
          WRITE(6,'(10I4)') (IOSPIR(2,I),I=1,NIRR_DG)
        END IF
      END IF
*
      RETURN
      END
*
      SUBROUTINE SHTOOB_KRCC(NSHPIR,NIRREP,MXPOBS,NSMOB,NOSPIR,IOSPIR,
     &                      NOBPS,NOB)
*
* Number of shells per irrep => Number of orbitals per symmetry
*
* =====
* Input
* =====
*
*  NSHPIR : Number of shells per irrep
*  NIRREP : Number of irreps
*  MXPOBS : Largest allowed number of orbitals symmetries
*  NSMOB  : Number of orbital symmetries
*  NOSPIR : Number of orbital symmetries per irrep
*  IOSPIR : Orbital symmetries per irrep
*
* ======
* Output
* ======
*  NOBPS  : Number of orbitals per symmetry
*  NOB    : Number of orbitals
*
* Jeppe Olsen, Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION NSHPIR(*),NOSPIR(*),IOSPIR(*)
*. Output
      DIMENSION NOBPS(*)
      CALL ISETVC(NOBPS,0,NSMOB)
      NOB = 0
      JSM = 0
      DO 100 IRREP = 1, NIRREP
C         IISM = IOSPIR(ISM,IRREP)
          IISM  = IOSPIR(IRREP)
          if (IISM.ne.0) then
             JSM = JSM + 1
             NOBPS(IISM) = NOBPS(IISM) + NSHPIR(IRREP)
          end if
          NOB = NOB + NSHPIR(IRREP)
  100 CONTINUE
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
         WRITE(6,*) ' SHTOOB Speaking '
         WRITE(6,*) ' =============== '
         WRITE(6,*) ' Number of orbitals obtained ', NOB
         WRITE(6,*) ' Number of orbitals per symmetry '
         CALL IWRTMA(NOBPS,1,NSMOB,1,NSMOB)
      END IF
*
      RETURN
      END
*
      SUBROUTINE ORBORD_GAS_KRCC(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB,NGSOBT,
     &                  NOCOBS,NTOOBS,NTOOB,
     &                  IREOST,IREOTS,ISFTO,ITFSO,IBSO,
     &                  NTSOB,IBTSOB,ITSOB,NOBPTS,IOBPTS,
     &                  ISFSO,ITFTO,NOBPT,IORD,ICOMP)
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
*
      SUBROUTINE DETTYP_KRCC
*
* Allowed set of determinants : Distribution of
* electrons in alpha(up) -and beta(down) spinors
*
* Defined in terms of
* NMS2VAL : Number of MS2 values (MK2 values )
* MS2VAL  : The allowed MS2 values
*
* Jeppe Olsen, July 97
* Timo Fleig, August '97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
#include "mxpdim.inc"
#include "crun.inc"
#include "lucinp.inc"
*. Output
#include "cstate.inc"
*
      NTEST = 00
*
      if (NTEST.ge.10) then
        write(6,*)
        write(6,*) ' -------------------------------- '
        write(6,*) '   DETTYP speaking  '
        write(6,*) ' -------------------------------- '
        write(6,*)
        write(6,'(A17,I4,I4)') 'MK2REF,MK2DEL are',MK2REF,MK2DEL
      end if
*
      IRELACI = 1
      IDBGCI = 1
      IF(IDBGCI.EQ.0.AND.IRELACI.EQ.0) THEN
*. Single MS2 space
        NMS2VAL = 1
        MS2VAL(1) = MK2REF
      ELSE
*. Potentially several MS2/MK2 values
        MK2MAX = MK2REF + MK2DEL
*. MK2MAX gives highest MK2 value for possibly coupled strings (?)
        if (MK2MAX.gt.NACTEL) then
          MK2MAX = NACTEL
        else
          call evenodd(IEVODMX,MK2MAX)
          call evenodd(IEVODN,NACTEL)
          if (IEVODN.eq.1.and.IEVODMX.eq.2) MK2MAX = MK2MAX - 1
        end if
        MK2MIN = MK2REF - MK2DEL
        IF(MK2MIN.LT.-NACTEL) MK2MIN = -NACTEL
* Pairing ? Ensure symmetry between +MK and - MK spaces
        IPAIR = 1
        call evenodd(IEVOD,NACTEL)
        IF(IPAIR.EQ.1) THEN
          IF(MK2MIN.GT.0) THEN
            MK2MINP = -MK2MIN
          ELSE
            if (IEVOD.eq.1) then
               MK2MIN = 1
               MK2MINP = -1
            else if (IEVOD.eq.2) then
               MK2MIN = 0
               MK2MINP = -2
            else
               write(6,*) 'Problem in evenodd'
               Call Abend2( 'Quitting in dettyp.' )
            end if
          END IF
          MK2MAXP = -MK2MAX
        ELSE
          if (IEVOD.eq.1) then
             MK2MINP = -1
             MK2MAXP = 1
          else if (IEVOD.eq.2) then
             MK2MINP = -2
             MK2MAXP = 0
          else
             write(6,*) 'Problem in evenodd'
             Call Abend2( 'Quitting in dettyp.' )
          end if
        END IF
*. And construct the MS2 spaces
        NMS2VAL = 0
        DO MS2 = MK2MIN,MK2MAX,2
          NMS2VAL = NMS2VAL + 1
          MS2VAL(NMS2VAL) = MS2
        END DO
        DO MS2 = MK2MAXP,MK2MINP,2                           !sic
          NMS2VAL = NMS2VAL+1
          MS2VAL(NMS2VAL) = MS2
        END DO
      END IF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,'(A,I2)') ' Number of MS2/MK2 values', NMS2VAL
        WRITE(6,*)
        WRITE(6,'(A)') ' Included MS2/MK2 values '
        WRITE(6,'(A)') ' ======================= '
        DO IMS2 = 1, NMS2VAL
          WRITE(6,'(1H ,I3)') MS2VAL(IMS2)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE SYMINF_KRCC
*
* Information about number of symmetries
*
* Input : /LUCINP/,/ORBINP
* Output : /CSM/,/CSMPRO/
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
#include "mxpdim.inc"
#include "lucinp.inc"
*

*. Output
* NSMSX : number of symmetries of single excitations
* NSMDX : Number of symmetries of double excitations
* NSMST : Number of symmetries of strings
* NSMCI : NUmber of symmetries of CI spaces
* ITSSX : Total symmetrix single excitation
* ITSDX : Total symmetrix double excitation
#include "csm.inc"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
* ADASX : symmetry of orbs i and i => symmetry of a+iaj
* ASXAD : symmetry of orb j and excit a+iaj => symmetry of i
* ADSXA : symmetry of orb i and excit a+iaj => symmetry of j
*
* SXSXDX : Symmetry of two single excitations
*          => symmetry of double  excitation
* SXDXSX : Symmetry of single excitation and double excitation
*          => symmetry of single  excitation
      NTEST = 0
*.
      IF(PNTGRP.EQ.1) THEN
* =====
* D 2 h
* =====
        CALL ZSYM1_KRCC(NIRREP,NTEST)
      ELSE IF(PNTGRP.EQ.2) THEN
* ========
* C inf V
* ========
        CALL ZNONAB_KRCC(0,MAXML,NSMOB,NTEST)
        CALL ZSYM2_KRCC(NTEST)
      ELSE IF(PNTGRP.EQ.3.OR.PNTGRP.EQ.4) THEN
* ===========
* D inf H O3
* ===========
        CALL ZNONAB_KRCC(1,MAXML,NSMOB,NTEST)
        CALL ZSYM2_KRCC(NTEST)
      ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
* ===========
* Double group
* ============
*
        CALL ZSYMDG_KRCC(NTEST)
      ELSE
        WRITE(6,*) ' You are too early , sorry '
        WRITE(6,*) ' Illegal PNTGRP in SYMINF ',PNTGRP
        Call Abend1( 11 )
      END IF
*
      RETURN
      END
*
      SUBROUTINE ZSYM1_KRCC(NIRREP,IPRNT)
*
* Number of symmetries for d2h
* Symmetry connecting arrays
* ( trivial, written for compatibility with higher point groups)
*
      INTEGER SYMPRO(8,8)
      DATA  SYMPRO/1,2,3,4,5,6,7,8,
     &             2,1,4,3,6,5,8,7,
     &             3,4,1,2,7,8,5,6,
     &             4,3,2,1,8,7,6,5,
     &             5,6,7,8,1,2,3,4,
     &             6,5,8,7,2,1,4,3,
     &             7,8,5,6,3,4,1,2,
     &             8,7,6,5,4,3,2,1 /
#include "csm.inc"
*
C     PARAMETER ( MXPOBS = 20 )
#include "mxpdim.inc"
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

      NSMSX = NIRREP
      NSMDX = NIRREP
      NSMST = NIRREP
      NSMCI = NIRREP
      ITSSX = 1

*
C     COPMT2(AIN,AOUT,NINR,NINC,NOUTR,NOUTC,IZERO)
      CALL ICPMT2(SYMPRO,ADASX,8,8,MXPOBS,MXPOBS,1)
      CALL ICPMT2(SYMPRO,ADSXA,8,8,MXPOBS,2*MXPOBS,1)
      CALL ICPMT2(SYMPRO,ASXAD,8,8,MXPOBS,2*MXPOBS,1)
      CALL ICPMT2(SYMPRO,SXSXDX,8,8,2*MXPOBS,2*MXPOBS,1)
      CALL ICPMT2(SYMPRO,SXDXSX,8,8,2*MXPOBS,4*MXPOBS,1)
*
      RETURN
      END
*
      SUBROUTINE ZSYM2_KRCC(IPRNT)
*
* Symmetry connecting arrays
*
* ======
*. Input
* ======
*
#include "mxpdim.inc"
#include "lucinp.inc"
*./NONAB/
      LOGICAL INVCNT
      COMMON/NONABR/ INVCNT,NIG,NORASM(MXPOBS),
     &              MNMLOB,MXMLOB,NMLOB,
     &              MXMLST,MNMLST,NMLST,
     &              NMLSX ,MNMLSX,MXMLSX,
     &              MNMLCI,MXMLCI,NMLCI,
     &              MXMLDX,MNMLDX,NMLDX
#include "csm.inc"
*
* ======
*.Output
* ======
*
*./CSMPRD/
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*
**. ADASX,ASXAD,ADSXA
      CALL ISETVC(ADASX,0,MXPOBS**2)
      CALL ISETVC(ASXAD,0,2*MXPOBS**2)
      CALL ISETVC(ADSXA,0,2*MXPOBS**2)
*
      DO 100 ISM = 1, NSMOB
C       MLSM_REL(IML,IPARI,ISM,TYPE,IWAY)
        CALL MLSM_KRCC(IML,IPARI,ISM,'OB',2)
        DO 90 JSM = 1, NSMOB
          CALL MLSM_KRCC(JML,JPARI,JSM,'OB',2)
*.a+ i a j symmetry
          IJML = IML - JML
          IF((IPARI.EQ.1.AND.JPARI.EQ.1).OR.
     &       (IPARI.EQ.2.AND.JPARI.EQ.2)) THEN
            IJPARI = 1
          ELSE
            IJPARI = 2
          END IF
          IJSM = (IJPARI-1)*NMLSX + IJML - MNMLSX + 1
          ADASX(ISM,JSM) = IJSM
          ASXAD(JSM,IJSM) = ISM
          ADSXA(ISM,IJSM) = JSM
   90   CONTINUE
  100 CONTINUE
*.SXSXDX,SXDXSX
      DO 200 ISX = 1, NSMSX
C       MLSM_REL(IML,IPARI,ISM,TYPE,IWAY)
        CALL MLSM_KRCC(IML,IPARI,ISX,'SX',2)
        DO 190 JSX = 1, NSMSX
          CALL MLSM_KRCC(JML,JPARI,JSX,'SX',2)
          IF((IPARI.EQ.1.AND.JPARI.EQ.1).OR.
     &       (IPARI.EQ.2.AND.JPARI.EQ.2)) THEN
            IJPARI = 1
          ELSE
            IJPARI = 2
          END IF
          IJML = IML + JML
          IJSM = (IJPARI-1)*NMLDX+IJML - MNMLDX + 1
          SXSXDX(ISX,JSX) = IJSM
          SXDXSX(ISX,IJSM) = JSX
  190   CONTINUE
  200 CONTINUE
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' ADASX '
         WRITE(6,*) ' ===== '
         CALL IWRTMA(ADASX,NSMOB,NSMOB,MXPOBS,MXPOBS)
         WRITE(6,*) ' ASXAD '
         WRITE(6,*) ' ===== '
         CALL IWRTMA(ASXAD,NSMOB,NSMSX,MXPOBS,2*MXPOBS)
         WRITE(6,*) ' ADSXA '
         WRITE(6,*) ' ===== '
         CALL IWRTMA(ADSXA,NSMOB,NSMSX,MXPOBS,2*MXPOBS)
         WRITE(6,*) ' SXSXDX'
         WRITE(6,*) ' ======'
         CALL IWRTMA(SXSXDX,NSMSX,NSMSX,2*MXPOBS,2*MXPOBS)
         WRITE(6,*) ' SXDXSX'
         WRITE(6,*) ' ======'
         CALL IWRTMA(SXDXSX,NSMSX,NSMDX,2*MXPOBS,4*MXPOBS)
      END IF
*
      RETURN
      END
*
      SUBROUTINE MLSM_KRCC(IML,IPARI,ISM,TYPE,IWAY)
*
* Transfer between ML,IPARI notation and compound notation ISM
*
* IWAY = 1 : IML,IPARI => Compound
* IWAY = 2 : IML,IPARI <= Compound
*
* TYPE : 'SX','OB','ST','DX','CI'
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
      CHARACTER*2 TYPE
*./NONAB/
      LOGICAL INVCNT
C     PARAMETER (MXPOBS = 20 )
      COMMON/NONABR/ INVCNT,NIG,NORASM(MXPOBS),
     &              MNMLOB,MXMLOB,NMLOB,
     &              MXMLST,MNMLST,NMLST,
     &              NMLSX ,MNMLSX,MXMLSX,
     &              MNMLCI,MXMLCI,NMLCI,
     &              MXMLDX,MNMLDX,NMLDX
#include "csm.inc"
*
*.(Tired of warnings from 3090 Compiler so )
* (
      NML = 0
      MXML= 0
      MNML= 0
*             )
      IF(TYPE.EQ.'OB') THEN
        NML = NMLOB
        MXML = MXMLOB
        MNML = MNMLOB
      ELSE IF(TYPE.EQ.'SX') THEN
        NML = NMLSX
        MXML = MXMLSX
        MNML = MNMLSX
      ELSE IF(TYPE.EQ.'DX') THEN
        NML = NMLDX
        MXML = MXMLDX
        MNML = MNMLDX
      ELSE IF(TYPE.EQ.'ST') THEN
        NML = NMLST
        MXML = MXMLST
        MNML = MNMLST
      ELSE IF(TYPE.EQ.'CI') THEN
        NML = NMLCI
        MXML = MXMLCI
        MNML = MNMLCI
      END IF
*
      IF(IWAY.EQ.1) THEN
C        ISM = (IPARI-1)*NML + MNML - 1
         ISM = (IPARI-1)*NML + IML - MNML + 1
      ELSE IF(IWAY.EQ.2) THEN
        IF(ISM.GT.NML) THEN
          IPARI = 2
          IML = ISM - NML + MNML - 1
        ELSE
          IPARI = 1
          IML = ISM       + MNML - 1
        END IF
      ELSE
        WRITE(6,*) ' Error in MLSM , IWAY = ' ,IWAY
        WRITE(6,*) ' MLSM stop !!! '
        Call Abend1( 20 )
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,A)') ' MLSM speaking ,type= ',TYPE
        WRITE(6,'(A,3I4)') ' IML IPARI ISM ',IML,IPARI,ISM
      END IF
*
      RETURN
      END
*
      SUBROUTINE ZNONAB_KRCC(INVCEN,MAXMLO,NSMOB,IPRNT)
*
*
* ============================
* Set up common block /NONAB/
* ============================
*
*========
* Input :
*========
*      INVCNT :inversion center is present(1), absent(0)
*      MAXMLO : Largest ML value of orbitals
*      NSMOB  : Number of symmetries of orbitals
*      Contents of common block /STRINP/,/ORBINP/
*=========
* output :
*=========
*======================
* Orbital Information
*======================
*      NORASM : Number of orbitals per abelian symmetry
*      MNMLOB : Smallest ML of orbitals
*      MXMLOB : largest ML of orbitals
*      NMLOB  : number of ML values for orbitals
*
*======================
* String Information
*======================
*      MNMLST : smallest ML of any strings
*      MXMLST : largest ML of any strings
*      NMLST  : Number of ML values of strings
*      NSMST  : Number of symmetries of strings
*
*==============================
* Single excitation Information
*==============================
*      MNMLSX : SMALLEST ML OF SINGLE EXCITATION
*      MXMLSX : LARGEST ML OF SINGLE EXCITATIONS
*      NMLSX  : NUMBER OF ML VALUES FOR SINGLE EXCITATIONS
*      NSMSX  : NUMBER OF SYMMETRIES FOR SINGLE EXCITATIONS
*=============================================
* External configurations (upto 4 electrons )
*=============================================
*      MNMLXT : SMALLEST ML OF External configurations
*      MXMLSX : LARGEST ML OF external configurations
*      NMLXT  : NUMBER OF ML VALUES FOR ext. configurations
*
* =============
* General input
* =============
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "strinp.inc"
* =======
*. Output
* =======
*./NONAB/
      LOGICAL INVCNT
      COMMON/NONABR/ INVCNT,NIG,NORASM(MXPOBS),
     &              MNMLOB,MXMLOB,NMLOB,
     &              MXMLST,MNMLST,NMLST,
     &              NMLSX ,MNMLSX,MXMLSX,
     &              MNMLCI,MXMLCI,NMLCI,
     &              MXMLDX,MNMLDX,NMLDX
*
      NTEST = 0
      NTEST = MAX(IPRNT,NTEST)
*. Inversion symmetry
      IF( INVCEN .EQ. 0 ) THEN
        INVCNT = .FALSE.
        NIG = 1
      ELSE
        INVCNT = .TRUE.
        NIG = 2
      END IF
*
** 1 : Information about orbitals
*
      MXMLOB = MAXMLO
      MNMLOB =-MAXMLO
      NMLOB = MXMLOB - MNMLOB + 1
*. Number of orbitals per symmetry
      DO 10 ISYM = 1, NSMOB
        NORASM(ISYM) = IFREQ(ISMFTO,ISYM,NACOB)
   10 CONTINUE
      IF( NTEST.GE. 2 ) THEN
        WRITE(6,*) ' NORASM '
        CALL IWRTMA(NORASM,1,NSMOB,1,NSMOB)
        WRITE(6,*) ' MNMLOB,MXMLOB ',MNMLOB,MXMLOB
        WRITE(6,*) ' NMLOB, NSMOB ',NMLOB,NSMOB
      END IF
*
**  2. Information about strings
*
      MXMLST = 0
      MNMLST = 0
      DO 50 ITYPE = 1, NSTTYP
        IEL = NELEC(ITYPE)
*
        MXMLTP = 0
        DO 40 IML = MXMLOB,MNMLOB,-1
          IORB = NORASM(IML-MNMLOB+1)
          IF(INVCNT) IORB = IORB + NORASM(NMLOB+IML-MNMLOB+1)
          IEL2 = MIN(IORB,IEL)
          MXMLTP = MXMLTP + IEL2*IML
          IEL = IEL - IEL2
   40   CONTINUE
        MXMLST = MAX(MXMLST,MXMLTP)
*
        MNMLTP = 0
        IEL = NELEC(ITYPE)
        DO 45 IML = MNMLOB,MXMLOB
          IORB = NORASM(IML-MNMLOB+1)
          IF(INVCNT) IORB = IORB + NORASM(NMLOB+IML-MNMLOB+1)
          IEL2 = MIN(IORB,IEL)
          MNMLTP = MNMLTP + IEL2*IML
          IEL = IEL - IEL2
   45   CONTINUE
        MNMLST = MIN(MNMLST,MNMLTP)
   50 CONTINUE
*
      NMLST  = MXMLST - MNMLST + 1
      NSMST  = NIG * NMLST
*
      IF( NTEST .GE. 2 ) THEN
        WRITE(6,*) ' MXMLST,MNMLST,NSMST'
        WRITE(6,*)   MXMLST,MNMLST,NSMST
      END IF
*
** 3. Information about single excitations
*
      MNMLSX = MNMLOB - MXMLOB
      MXMLSX = MXMLOB - MNMLOB
      NMLSX  = MXMLSX - MNMLSX +1
      NSMSX  = NIG * NMLSX

      IF( NTEST .GE.2 ) THEN
        WRITE(6,*) ' NMLSX,NSMSX,MNMLSX ',NMLSX,NSMSX,MNMLSX
      END IF
*
** 4 : External configurations(double excitations)
*
      MXMLDX = 4*MAXMLO
      MNMLDX = -4*MAXMLO
      NMLDX  = MXMLDX - MNMLDX + 1
      NSMDX  = NIG * NMLDX
      IF( NTEST .GE.2 ) THEN
        WRITE(6,*) ' NMLDX,NSMDX,MNMLDX ',NMLDX,NSMDX,MNMLDX
      END IF
*
** 5 : Determinants
*
      MXMLCI =  2*MXMLST + MXMLDX
      MNMLCI = - MXMLCI
      NMLCI = 2 * MXMLCI + 1
      NSMCI = NIG * NMLCI
*
*.6 Total symmetrix single excitation and external
*
      ITSSX = 0 - MNMLSX + 1
      ITSDX = 0 - MNMLDX + 1

      IF ( NTEST .GE. 1 ) THEN
        WRITE(6,*)
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of orbitals     .. ', NSMOB
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of strings      .. ', NSMST
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of single excit. . ', NSMSX
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of double excit. . ', NSMDX
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of determinants .. ', NSMCI
        WRITE(6,*)
*
        WRITE(6,*) ' Total symmetric single excitation .. ',ITSSX
        WRITE(6,*) ' Total symmetric double excitation .. ',ITSDX
      END IF
*
      RETURN
      END
*
      SUBROUTINE ZSYMDG_KRCC(IPRNT)
*
* Number of symmetries for d2h
* Symmetry connecting arrays
*
*
#include "mxpdim.inc"
#include "csm.inc"
#include "symm.inc"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*
      NSMSX = NIRR_DG
      NSMDX = NIRR_DG
      NSMST = NIRR_DG
      NSMCI = NIRR_DG
      ITSSX = 1

*
      IZERO = 0
      CALL ISETVC(ADASX,IZERO,MXPOBS**2)
      CALL ISETVC(ASXAD,IZERO,2*MXPOBS**2)
      CALL ISETVC(ADSXA,IZERO,2*MXPOBS**2)
*
      DO IREP = 1, NIRR_DG
        DO JREP = 1, NIRR_DG
*A+ A => SX
          ADASX(IREP,JREP) = IDBGMULT(IREP,IADJSYM(JREP))
        END DO
      END DO
*
      DO IREP = 1, NIRR_DG
        DO JREP = 1, NIRR_DG
          DO IJSM = 1, NIRR_DG
            IF(ADASX(IREP,JREP).EQ.IJSM) THEN
              ADSXA(IREP,IJSM) = JREP
              ASXAD(JREP,IJSM) = IREP
            END IF
          END DO
        END DO
      END DO

      CALL ICPMT2(IDBGMULT,SXSXDX,8,8,2*MXPOBS,2*MXPOBS,1)
      CALL ICPMT2(IDBGMULT,SXDXSX,8,8,2*MXPOBS,4*MXPOBS,1)
*
      RETURN
      END
*
      SUBROUTINE OCCLS_KRCC(IWAY,NOCCLS,IOCCLS,NEL,NGAS,IGSMIN,IGSMAX)
*
* IWAY = 1 :
* obtain NOCCLS =
* Number of allowed ways of distributing the orbitals in the
* active spaces
*
* IWAY = 2 :
* OBTAIN NOCCLS and
* IOCCLS = allowed distributions of electrons
*
*
*
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION IGSMIN(NGAS),IGSMAX(NGAS)
*. Output
      DIMENSION  IOCCLS(NGAS,*)
*. Local scratch
#include "mxpdim.inc"
      DIMENSION IOCA(MXPNGAS),IOC(MXPNGAS)
*
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
         WRITE(6,*)  ' OCCLS in action '
         WRITE(6,*) ' =================='
         WRITE(6,*) ' NGAS NEL ', NGAS,NEL
      END IF
*
      ISKIP = 1
      NOCCLS = 0
*. start with smallest allowed number
      DO IGAS = 1, NGAS
        IOCA(IGAS) = IGSMIN(IGAS)
      END DO
      NONEW = 0
      IFIRST = 1
*. Loop over possible occupations
 1000 CONTINUE
        IF(IFIRST.EQ.0) THEN
*. Next accumulated occupation
          CALL NXTNUM3(IOCA,NGAS,IGSMIN,IGSMAX,NONEW)
        END IF
        IF(NONEW.EQ.0) THEN
*. ensure that IOCA corresponds to an accumulating occupation,
*. i.e. a non-decreasing sequence
        IF(ISKIP.EQ.1) THEN
          KGAS = 0
          DO IGAS = 2, NGAS
            IF(IOCA(IGAS-1).GT.IOCA(IGAS)) KGAS = IGAS
          END DO
          IF(KGAS .NE. 0 ) THEN
            DO IGAS = 1, KGAS-1
              IOCA(IGAS) = IGSMIN(IGAS)
            END DO
            IOCA(KGAS) = IOCA(KGAS)+1
          END IF
        END IF
C?      WRITE(6,*) ' Another accumulated occupation: '
C?      CALL IWRTMA(IOCA,1,NGAS,1,NGAS)
*. corresponding occupation of each active space
        NEGA=0
        DO IGAS = 1, NGAS
          IF(IGAS.EQ.1) THEN
            IOC(IGAS) = IOCA(IGAS)
          ELSE
            IOC(IGAS) = IOCA(IGAS)-IOCA(IGAS-1)
            IF(IOC(IGAS).LT.0) NEGA = 1
          END IF
        END DO
C?      WRITE(6,*) ' Another occupation: '
C?      CALL IWRTMA(IOC,1,NGAS,1,NGAS)
        IFIRST = 0
*. Correct number of electrons
        IEL = IELSUM(IOC,NGAS)
        IF(IEL.EQ.NEL.AND.NEGA.EQ.0) THEN
          NOCCLS = NOCCLS + 1
          IF(IWAY.EQ.2) THEN
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' Another allowed class : '
              CALL IWRTMA(IOC,1,NGAS,1,NGAS)
            END IF
            CALL ICOPVE(IOC,IOCCLS(1,NOCCLS),NGAS)
          END IF
        END IF
      END IF
      IF(NONEW.EQ.0) GOTO 1000
*
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Number of Allowed occupation classes ', NOCCLS
         IF(IWAY.EQ.2.AND.NTEST.GE.20) THEN
           WRITE(6,*) ' Occupation classes '
           CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
         END IF
      END IF
*
      RETURN
      END
*
      SUBROUTINE INTDIM_KRCC
*
* Number of integrals and storage mode
*
*. Last modifications : Jan 98, relativistic integrals added
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* =====
*.Input
* =====
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "crun.inc"
#include "cgas.inc"
*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),
     &              SXDXSX(2*MXPOBS,4*MXPOBS)
*. Local arrays
      INTEGER NNORB(MXNDGIRR,4)
*
* =======
*. Output
* =======
*. Atmost 10 double group integral lists
#include "cintfo.inc"
*
      NTEST = 0
*
      if (NTEST.ge.2)
     &    WRITE(6,'(A25,I3)') ' INTDIM speaking, IDBG = ', IDBG
*
* =================================
* integrals in point group symmetry
* =================================
*
      IPNTGRP_SAVE = PNTGRP
*. Force the use of point group symmetry
      IF(IDBG.EQ.1) PNTGRP = 1
*.1 : Number of one-electron integrals
      NINT1 =  NSXFSM2(NSMSH,NGSSHT,NGSSHT,ITSSX,1,IPRNT)
C              NSXFSM2(NSMOB,NO1PS,NO2PS,ISXSM,ISYM,IPRNT)
*.2 : Number of two-electron integrals
*. Full eightfold symmetry can be used
      I12S = 1
      I34S = 1
      I1234S = 1
      NINT2 = NDXFSM2(NSMSH,NGSSHT,NGSSHT,NGSSHT,NGSSHT,
     &                I12S,I34S,I1234S,IPRNT )
C     NDXFSM2(NSMOB,NO1PS,NO2PS,NO3PS,NO4PS,
C    &         IS12,IS34,IS1234,IPRNT)
*
* =======================
* Double group integrals
* =======================
*
      IF(IDBG.EQ.0) THEN
        NL1D = 0
        NL2D = 0
      ELSE
        PNTGRP = IPNTGRP_SAVE
*. Number of double group integral lists
        CALL DBG_INT_LISTS_KRCC(NL1D,IUB1D,NL2D,IUB2D,ISM2D,IDENSI,
     &                     ISPINFREE)
*
        DO IL1 = 1, NL1D
          DO INDEX = 1, 2
            IF(IUB1D(INDEX,IL1).EQ.1) THEN
              CALL ICOPVE(NTOOBS,NNORB(1,INDEX),NSMOB)
            ELSE
              CALL ICOPVE(NTOOBS2,NNORB(1,INDEX),NSMOB)
            END IF
          END DO
          L1D(IL1) = NSXFSM2(NSMOB,NNORB(1,1),NNORB(1,2),ITSSX,1,IPRNT)
        END DO
*
        DO IL2 = 1, NL2D
          DO INDEX = 1, 4
            IF(IUB2D(INDEX,IL2).EQ.1) THEN
              CALL ICOPVE(NTOOBS,NNORB(1,INDEX),NSMOB)
            ELSE
              CALL ICOPVE(NTOOBS2,NNORB(1,INDEX),NSMOB)
            END IF
          END DO
          L2D(IL2) = NDXFSM2(NSMOB,
     &               NNORB(1,1),NNORB(1,2),NNORB(1,3),NNORB(1,4),
     &               ISM2D(1,IL2),ISM2D(2,IL2),ISM2D(3,IL2),
     &               IPRNT )
        END DO
      END IF
*. Largest block of two-electron integrals
      LMAXE2LST = NINT2
      DO IL2 = 1, NL2D
        LMAXE2LST = MAX(LMAXE2LST,L2D(IL2))
      END DO
*
      write(6,*)
      write(6,'(A,I10)')
     &   '   The largest block of two-electron integrals is',LMAXE2LST
*
      NBINT1 = NSMOB
      NBINT2 = NSMOB ** 3
*
      RETURN
      END
*
      SUBROUTINE DBG_INT_LISTS_KRCC(NL1,IUB1,NL2,IUB2,ISM2,IDENSI,
     &                         ISPINFREE)
*
* Number of double group integral lists to be constructed
*
*     Jeppe Olsen, Feb. 98
*
* Reordered lists for convenience: We only need the list
* (uu|uu) in an initial or plain SOCI calculation. So this
* is going to be list number 1. Allocation modified as well.
* Modifications in GETINT_DBG incorporated (and elsewhere).
* We are furthermore using real wave functions. This implies
* That complex conjugation symmetry can also be used in the
* quaternion matrix groups C1 and Ci, because their irrep
* basis functions are not generally complex. But if we
* generalize the code to e.g. atomic complex 4- (or 2-)spinors
* or the like, these things have to be distinguished.
* So we will have 3 cases:
*
*   1) Initial or plain SOCI: 1 list
*   2) SOCI + spinor optimization in some way (INS or MCSCF): 3 lists
*   3) Complex wave function: 6 lists
*
* Here, time-reversal symmetry has already been accounted for in
* all cases, so these are non-redundant integral lists.
* At the moment, the first two cases are assumed to be implemented,
* so we will distinguish them by the input keyword DENSI, i.e.
* whether density matrices are to be calculated or not.
*
*     Timo Fleig, June 1999
*
* revised for DIRAC environment
*     Timo Fleig, December 2000
*
      IMPLICIT REAL*8(A-H,O-Z)
*.output
      INTEGER IUB1(2,10),IUB2(4,10),ISM2(3,10)
*
*. ======================================
*. 1e Lists : all four lists constructed
*. ======================================
*
      NL1 = 4
*. 1 : unbarred, unbarred
      IUB1(1,1) = 1
      IUB1(2,1) = 1
*. 2 : unbarred, barred
      IUB1(1,2) = 1
      IUB1(2,2) = 2
*. 3 : barred, unbarred
      IUB1(1,3) = 2
      IUB1(2,3) = 1
*. 4 : barred, barred
      IUB1(1,4) = 2
      IUB1(2,4) = 2
*
*. =========
*. 2e Lists
*. =========
*
*  For the moment, no complex wave function:
      IWFCOMP = 0
*
* The lists of integrals are
*
* 1 : (u u | u u)
*--------------------
* 2 : (u b | u u)
* 3 : (b u | u b)
* 4 : (b u | b u)
*--------------------
* 5 : (u b | u b)
* 6 : (b u | u u)
*--------------------
*
*=====================
* List 1 ( u u ! u u )
*=====================
*
      IUB2(1,1) = 1
      IUB2(2,1) = 1
      IUB2(3,1) = 1
      IUB2(4,1) = 1
*. No Symmetry between 12 and 34
      ISM2(1,1) = 0
      ISM2(2,1) = 0
      ISM2(3,1) = 0
*
      if (IDENSI.ne.0.or.ISPINFREE.eq.0) then
*
*=====================
* List 2 ( u b ! u u )
*=====================
*
        IUB2(1,2) = 1
        IUB2(2,2) = 2
        IUB2(3,2) = 1
        IUB2(4,2) = 1
*. No Symmetry between 12 and 34
        ISM2(1,2) = 0
        ISM2(2,2) = 0
        ISM2(3,2) = 0
*
*=====================
* List 3 ( b u ! u b )
*=====================
*
        IUB2(1,3) = 2
        IUB2(2,3) = 1
        IUB2(3,3) = 1
        IUB2(4,3) = 2
*. No Symmetry between 12 and 34
        ISM2(1,3) = 0
        ISM2(2,3) = 0
        ISM2(3,3) = 0
*
*=====================
* List 4 ( b u ! b u )
*=====================
*
        IUB2(1,4) = 2
        IUB2(2,4) = 1
        IUB2(3,4) = 2
        IUB2(4,4) = 1
*. No Symmetry between 12 and 34
        ISM2(1,4) = 0
        ISM2(2,4) = 0
        ISM2(3,4) = 0
*
        if (IWFCOMP.eq.1) then
*
*=====================
* List 5 ( u b ! u b )
*=====================
*
          IUB2(1,5) = 1
          IUB2(2,5) = 2
          IUB2(3,5) = 1
          IUB2(4,5) = 2
*. Symmetry between 12 and 34
          ISM2(1,5) = 0
          ISM2(2,5) = 0
          ISM2(3,5) = 0
CTEMP     ISM2(3,5) = 1
*
*=====================
* List 6 ( b u ! u u )
*=====================
*
          IUB2(1,6) = 2
          IUB2(2,6) = 1
          IUB2(3,6) = 1
          IUB2(4,6) = 1
*. No Symmetry between 12 and 34
          ISM2(1,6) = 0
          ISM2(2,6) = 0
          ISM2(3,6) = 0
*
        end if
      end if
*
      if (IDENSI.eq.0.and.ISPINFREE.eq.1) then
        NL2 = 1
      else
        NL2 = 4
        if (IWFCOMP.eq.1) NL2 = 6
      end if
*
      RETURN
      END
*
      SUBROUTINE LCISPC_KRCC(WORK,KFREE,LFREE,LCSBLK)
*
* Number of dets and combinations
* per symmetry for each type of internal space
*
* Jeppe Olsen , Winter 1994/1995 ( woops !)
*               July 1997, multiple Ms spaces allowed
*
* GAS VERSION
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
* ===================
*.Input common blocks
* ===================
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "strbas.inc"
#include "csm.inc"
#include "stinf.inc"
#include "cgas.inc"
#include "gasstr.inc"
*
* ====================
*. Output common block : XISPSM is calculated
* ====================
*
#include "cicisp.inc"
*
      DIMENSION WORK(*)
C
C
C#include "memint.h"
*
      NTEST = 200
*
      CALL QENTER('LCISP')
*. Number of spaces
      NICISP = NCMBSPC
*
*.Local memory
*
*.  Largest NOCTPA*NOCTPB
      print*,'entering LCISPC_KRCC'
      CALL MEMCHK_KRCC(WORK)
      MXNOCAB = 0
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
        NOCAB = NOCTYP(IATP)*NOCTYP(IBTP)
        MXNOCAB = MAX(MXNOCAB,NOCAB)
      END DO
C     WRITE(6,*) ' MXNOCAB = ',MXNOCAB
      print*,'MXNOCAB',MXNOCAB
      CALL MEMGET('REAL',KLIOIO,MXNOCAB,WORK,KFREE,LFREE)
      print*,'KFREE,LFREE',KFREE,LFREE
      CALL TEST_ARRAY(WORK(KLIOIO),MXNOCAB)
      CALL MEMCHK_KRCC(WORK)
      print*,'after KLIOIO'
*
      print*,'NSMST',NSMST
      CALL MEMGET('REAL',KLBLTP,NSMST,WORK,KFREE,LFREE)
      print*,'after KLBLTP'
*. Obtain array giving symmetry of sigma v reflection times string
*. symmetry.
*. Array defining symmetry combinations of internal strings
*. Number of internal dets for each symmetry
      CALL SMOST_KRCC(NSMST,NSMCI,MXPCSM,ISMOST)
*
      MXSB = 0
* MXSOOB_AS added by Lasse 20-01-04
      MXSOOB_AS = 0
      MXSOOB = 0
      CALL MEMCHK_KRCC(WORK)
*. Loop over MS2 spaces
      print*,'NMS2VAL,NICISP',NMS2VAL,NICISP
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
*
        NOCTPA =  NOCTYP(IATP)
        NOCTPB =  NOCTYP(IBTP)
*
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
*
C?      WRITE(6,*) ' IMS2 IATP IBTP NOCTPA NOCTPB '
C?      WRITE(6,*)   IMS2,IATP,IBTP,NOCTPA,NOCTPB
*
        DO ICI = 1, NICISP
*. allowed combination of types
          CALL IAIBCM_GAS(LCMBSPC(ICI),ICMBSPC(1,ICI),
     &                    IGSOCCX,NOCTPA,
     &                    NOCTPB,ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &                    NELFGP,MXPNGAS,NGAS,WORK(KLIOIO),0)
      CALL MEMCHK_KRCC(WORK)
*
          DO  ISYM = 1, NSMCI
            CALL ZBLTP_KRCC(ISMOST(1,ISYM),NSMST,IDC,WORK(KLBLTP))
            CALL NGASDT_KRCC(IGSOCCX(1,1,ICI),IGSOCCX(1,2,ICI),
     &           NGAS,ISYM,
     &           NSMST,NOCTPA,NOCTPB,WORK(KNSTSO(IATP)),
     &           WORK(KNSTSO2(IBTP)),
     &           ISPGPFTP(1,IBSPGPFTP(IATP)),
     &           ISPGPFTP(1,IBSPGPFTP(IBTP)),MXPNGAS,NELFGP,
     &           NCOMB,XNCOMB,MXS,MXSOO,WORK(KLBLTP),NTTSBL,
     &           LCOL,WORK(KLIOIO),MXSOO_AS)
*
            IF(IMS2.EQ.1) THEN
              XISPSM(ISYM,ICI) = XNCOMB
              NBLKIC(ISYM,ICI) = NTTSBL
              LCOLIC(ISYM,ICI) = LCOL
            ELSE
              XISPSM(ISYM,ICI) = XISPSM(ISYM,ICI) + XNCOMB
              NBLKIC(ISYM,ICI) = NBLKIC(ISYM,ICI) + NTTSBL
              LCOLIC(ISYM,ICI) = LCOLIC(ISYM,ICI) + LCOL
            END IF
            MXSOOB = MAX(MXSOOB,MXSOO)
            MXSB = MAX(MXSB,MXS)
            MXSOOB_AS = MAX(MXSOO_AS,MXSOOB_AS)
          END  DO
*         ^ End of loop over symmetries
        END DO
*       ^ End of loop over CI spaces
      END DO
*     ^ End of loop over MS2 spaces
*
* Check for minimum size of resolution matrices:
      if (MXSB.gt.LCSBLK) then
        write(6,*) ' Exceeds space for resolution matrices: ',LCSBLK
        write(6,*) ' Increase LCSBLK (in krcc_inp) to at least ',MXSB
        call quit ('Stop in LCISPC_REL.')
      end if
*
      if (NTEST.ge.1) then
        WRITE(6,*)
     &  ' Number of internal combinations per symmetry '
        WRITE(6,*)
     &  ' =========================================== '
*
        DO 200 ICI = 1, NCMBSPC
          WRITE(6,*) ' CI space ', ICI
          CALL WRTMAT(XISPSM(1,ICI),1,NSMCI,1,NSMCI)
  200   CONTINUE
        WRITE(6,*) ' Largest symmetry block           ',MXSB
        WRITE(6,*) ' Largest Symmetry-type-type block ',MXSOOB
*
        WRITE(6,*) ' Number of TTS subblocks per CI expansion '
        WRITE(6,*) ' ======================================== '
*
        DO ICI = 1,  NCMBSPC
          WRITE(6,*) ' Internal CI space ', ICI
          CALL IWRTMA(NBLKIC(1,ICI),1,NSMCI,1,NSMCI)
        END DO
      end if
*
*. Largest number of BLOCKS in a CI expansion
      MXNTTS = 0
      DO ICI = 1,NCMBSPC
       DO ISM =1, NSMCI
        MXNTTS = MAX(MXNTTS,NBLKIC(ISM,ICI))
       END DO
      END DO
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Largest number of blocks in CI expansion',MXNTTS
        WRITE(6,*) ' Number of columns per CI expansion '
        WRITE(6,*) ' =================================== '
*
        DO ICI = 1,NCMBSPC
          WRITE(6,*) ' Internal CI space ', ICI
          CALL IWRTMA(LCOLIC(1,ICI),1,NSMCI,1,NSMCI)
        END DO
      END IF
*
* Release memory
      CALL MEMREL('LCISP',WORK,KLIOIO,KLIOIO,KFREE,LFREE)
*
      CALL QEXIT('LCISP')
*
      RETURN
      END
*
      SUBROUTINE SMOST_KRCC(NSMST,NSMCI,MXPCSM,ISMOST)
*
* ISMOST(ISYM,ITOTSM) : Symmetry of an internal state is ITOTSM
*                       if symmetry of 1 string is ISYM, the
*                       symmetry of the other string is
*                       ISMOST(ISYM,ITOTSM)
*
* Jeppe Olsen , Spring of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ISMOST(MXPCSM,MXPCSM)
*
      DO 1000 ITOTSM = 1, NSMCI
       DO 900 ISTSM  = 1, NSMST
C            SYMCOM_REL(ITASK,IOBJ,I1,I2,I12)
        CALL SYMCOM_KRCC(2,1,ISTSM,JSTSM,ITOTSM)
        ISMOST(ISTSM,ITOTSM) = JSTSM
  900  CONTINUE
 1000 CONTINUE
*
      NTEST = 0
      IF( NTEST.NE. 0 ) THEN
        WRITE(6,*) ' ==============='
        WRITE(6,*) ' Info from SMOST '
        WRITE(6,*) ' ==============='
        DO 1010 ITOTSM = 1, NSMCI
          WRITE(6,*) ' ISMOST array for ITOTSM = ', ITOTSM
          CALL IWRTMA(ISMOST(1,ITOTSM),1,NSMST,1,NSMST)
 1010   CONTINUE
      END IF
*
      RETURN
      END
*
*
* Codes for general symmetry handling
*
*                - ZSTINF : generate /STINF/ info on strings and mapping
*                - MEMSTR : allocates memory for string information
*                - WEIGHT : Weights for strings
*                - NSTRSO : Number of strings per sym and class
*                - ZBASE  : offset arrays for strings
*                - ZSMCL  : symmetry and class for each string
*                - GENSTR : Generate strings ordered by sym and class
*                - MEMEXT : Memory for external blocks
*
      SUBROUTINE ZBLTP_KRCC(ISMOST,MAXSYM,IDC,ICBLTP)
*
* Generate vector ICBLTP giving type of each block
*
* Simplified version : IDC = 3, IDC = 4 path eliminated
*
*
* ICBLTP gives type of symmetry block :
* = 0 : symmetry block is not included
* = 1 : symmetry block is included , all OO types
* = 2 : symmetry block is included , lower OO types
*
*. Input
      DIMENSION ISMOST(*)
*. Output
      DIMENSION ICBLTP(*)
*
      DO 100 IASYM = 1, MAXSYM
*
        IBSYM = ISMOST(IASYM)
        IF(IBSYM .EQ. 0 ) GOTO 100
        IF(IDC.EQ.2.AND.IBSYM.GT.IASYM) THEN
*.Symmetry block excluded
          ICBLTP(IASYM) = 0
        ELSE IF((IDC.EQ.2.AND.IASYM.GT.IBSYM).OR.IDC.EQ.1) THEN
*.Complete symmetry block included
          ICBLTP(IASYM) = 1
        ELSE
*.Lower half  symmetry block included
          ICBLTP(IASYM) = 2
        END IF
  100 CONTINUE
*
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Block type of symmetry blocks '
         CALL IWRTMA(ICBLTP,1,MAXSYM,1,MAXSYM)
      END IF
*
      RETURN
      END
*
      SUBROUTINE NGASDT_KRCC(IOCCMN,IOCCMX,NGAS,ITOTSM,
     &                      NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,
     &                      IAOCC,IBOCC,MXPNGAS,NELFGP,
     &                      NCOMB,XNCOMB,MXSB,MXSOOB,
     &                      IBLTP,NTTSBL,LCOL,IOCOC,MXSOOB_AS)
*
* Number of combinations with symmetry ITOTSM and
* occupation between IOCCMN and IOCCMX
*
* In view of the limited range of I*4, the number of dets
* is returned as integer and  real*8
*
* MXSB is largest UNPACKED symmetry block
* MXSOOB is largest UNPACKED symmetry-type-type block
* NTTSBL is number of TTS blocks in vector
* LCOL is the sum of the number of columns in each block
*
*
* Winter 94/95
* Winter 05 MXSOOB_AS added, Lasse
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Allowed combinations of alpha and beta types
      INTEGER IOCOC(NOCTPA,NOCTPB)
*. Occupation constraints
      DIMENSION IOCCMN(NGAS),IOCCMX(NGAS)
*. Occupation of alpha and beta strings
      DIMENSION IAOCC(MXPNGAS,*),IBOCC(MXPNGAS,*)
*. Number of strings per supergroup and symmetry
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*),NELFGP(*)
*. block types
      DIMENSION IBLTP(*)
*
      CALL QENTER('NGASD')
      NTEST = 000
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' NGASDT speaking'
        WRITE(6,*) ' ==============='
        WRITE(6,*) ' NGAS NOCTPA,NOCTPB ',NGAS,NOCTPA,NOCTPB
        WRITE(6,*) ' ITOTSM ', ITOTSM
        WRITE(6,*) ' Upper and lower occupation constraints'
        CALL IWRTMA(IOCCMN,1,NGAS,1,NGAS)
        CALL IWRTMA(IOCCMX,1,NGAS,1,NGAS)
        WRITE(6,*) ' IOCOC matrix '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        WRITE(6,*) ' Number of alpha and beta strings '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
*
      MXSB = 0
      MXSOOB = 0
      NCOMB = 0
      XNCOMB = 0.0D0
      NTTSBL = 0
      LCOL = 0
      MXSOOB_AS = 0
      DO 300 IASM = 1, NSMST
        IF(IBLTP(IASM).EQ.0) GOTO 300
        CALL SYMCOM_KRCC(2,1,IASM,IBSM,ITOTSM)
C       write(6,*) ' Iasm ibsm ', iasm,ibsm
        LSB = 0
        IF(IBSM.NE.0) THEN
          IF(IBLTP(IASM).EQ.2) THEN
            ISYM = 1
          ELSE
            ISYM = 0
          END IF
          DO 200 IATP = 1, NOCTPA
           IF(ISYM.EQ.1) THEN
C Lasse change
C            MXBTP = IATP
             MXBTP = MIN(IATP,NOCTPB)
           ELSE
             MXBTP = NOCTPB
           END IF
           LTSSA = NSSOA(IASM,IATP)
           XLTSSA = FLOAT(LTSSA)
           DO 100 IBTP = 1, MXBTP
*
             IF(NTEST.GE.10) THEN
               WRITE(6,*) ' Alpha super group and beta super group'
               CALL IWRTMA(IAOCC(1,IATP),1,NGAS,1,NGAS)
               CALL IWRTMA(IBOCC(1,IBTP),1,NGAS,1,NGAS)
             END IF
*
             IF(IOCOC(IATP,IBTP).EQ.1) THEN
*. Size of unpacked block
               LTTSUP =  LTSSA*NSSOB(IBSM,IBTP)
C              print*,'LTTSUP,LTSSA,NSSOB(IBSM,IBTP),IBSM,IBTP',
C    &                 LTTSUP,LTSSA,NSSOB(IBSM,IBTP),IBSM,IBTP
*. Size of packed block
               IF(ISYM.EQ.0.OR.IATP.NE.IBTP) THEN
                 LTTSBL = LTSSA*NSSOB(IBSM,IBTP)
                 XNCOMB = XNCOMB + XLTSSA*FLOAT(NSSOB(IBSM,IBTP))
               ELSE
                 LTTSBL = LTSSA *(LTSSA +1)/2
                 XNCOMB = XNCOMB + XLTSSA*(XLTSSA+1.0D0)*0.5
               END IF
               NCOMB = NCOMB + LTTSBL
               LSB = LSB + LTTSUP
               MXSOOB = MAX(MXSOOB,LTTSUP)
               NTTSBL = NTTSBL + 1
               LCOL = LCOL + NSSOB(IBSM,IBTP)
             END IF
  100      CONTINUE
  200     CONTINUE
          MXSB = MAX(MXSB,LSB)
        END IF
            MXSOOB_AS = MAX(MXSOOB_AS,LSB)
  300 CONTINUE
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' NGASDT : NCOMB XNCOMB ,NTTSBL',
     &               NCOMB,XNCOMB,NTTSBL
      END IF
*
      CALL QEXIT('NGASD')
*
      RETURN
      END
*
      SUBROUTINE GET_GASOCC_KRCC(JCMBSPC)
*
* Does what the name implies. Though it only tells if a GAS is occupied
* or not (excitations will be made into it or not) for a given JCMBSPC
*
* 1 = occupied
* 0 = unoccupied
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*
C     print*,'enter GET_GASOCC_KRCC'
C     print*,'JCMBSPC',JCMBSPC
      NTEST = 0
      NELTOT2 = IGSOCCX(NGAS,2,JCMBSPC)
C     print*,'NELTOT2',NELTOT2
      DO IGAS=1,NGAS
        IDIFF = IGSOCCX(IGAS,2,JCMBSPC) - IGSOCCX(IGAS,1,JCMBSPC)
        IF(IDIFF.eq.0.and.IGAS.ne.1) THEN
          IF(NGASOCC(IGAS-1).eq.0.or.
     &      IGSOCCX(IGAS-1,1,JCMBSPC).eq.NELTOT2) THEN
            NGASOCC(IGAS) = 0
          ELSE
            NGASOCC(IGAS) = 1
          END IF
        ELSE IF(IDIFF.eq.0.and.IGAS.eq.1) THEN
          IF(NELTOT2.le.2*NGSOBT(IGAS)) THEN
C            NGASOCC(IGAS) = 0
C          ELSE
            NGASOCC(IGAS) = 1
          END IF
        ELSE
          NGASOCC(IGAS) = 1
        END IF
        IF(NTEST.ge.1) THEN
          WRITE(6,*) 'IGAS,NGASOCC',IGAS,NGASOCC(IGAS)
        END IF
      END DO
* Check if its ok (if not probably a bug in this routine)
      DO IGAS=2,NGAS
        IDIFF = NGASOCC(IGAS-1) - NGASOCC(IGAS)
        IF(IDIFF.lt.0) THEN
          CALL ABEND2('BUG IN GET_GASOCC')
        END IF
      END DO
      RETURN
      END
*
      SUBROUTINE SETUP_FOR_itrctl_KRCC
C
C Due to some silly renaming this routine is here
C 
      implicit real*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "symm.inc"
#include "cgas.inc"
*  Setup for later
      NOBPIRR = 0
      IOFPGSM(1) = 0
      do IRR=1,NIRR_DG,1
         do IGAS=1,NGAS,1
            NOBPIRR = NOBPIRR + NGSSH(IRR,IGAS)
         end do
         if (IRR.lt.NIRR_PN) IOFPGSM(IRR+1) = IOFPGSM(IRR) + NOBPIRR
         ITOTDM(IRR) = NOBPIRR
         SYDI(IRR) = ITOTDM(IRR)
         NOBPIRR = 0
      end do
*
      END
*
      subroutine sspir_KRCC
*
      implicit real*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "symm.inc"
#include "coeff.inc"
#include "crun.inc"
*
*
* Number and symmetries of orbitals corresponding to a given shell
*
* =====
* Input
* =====
*
*   DOUGRP  : type of double group
*   NIRREP : Number of irreducible representations per double group
*   MXPIRR : Largest allowed number of shell irreps
*   MXPSBS : Largest allowed number of spinor symmetries
*   MXNSPIR : Maximum number of spatial irreps
*   MXNSPIPDMGRP: Maximum number of spinors per dimension group
*
* ======
* Output
* ======
*
*
* Jeppe Olsen , Winter of 1991
*
* modified: Timo Fleig, Fall of 1996
*
      integer ISSPIR,NSSPIR
*. Output
*
*  Scratch
      do IRREP=1,NIRR_DG,1
         NBSPFSM(IRREP) = 0
         NUSPFSM(IRREP) = 0
      end do
*
      if (DOUGRP.eq.1) then
*=======
* D2h(*)
*=======
        NSMOB = 0
      else if (DOUGRP.eq.2.or.DOUGRP.eq.3) then
         write(*,*) 'These cases have not yet been implemented.'
         Call Abend1( 5 )
*=======
* C2v(*)
*=======
      else if (DOUGRP.eq.4) then
*=======
* C2h(*)
*=======
         SYMDIM(1) = ITOTDM(1) + ITOTDM(2)
         SYMDIM(2) = ITOTDM(3) + ITOTDM(4)
         SYMDIM_F(1) = NORB_F(1) + NORB_F(2)
         SYMDIM_F(2) = NORB_F(3) + NORB_F(4)
         DDIM = (SYMDIM(1))*2 + (SYMDIM(2))*2
         NBSPFSM(5) = ITOTDM(2)
         NBSPFSM(6) = ITOTDM(1)
         NBSPFSM(7) = ITOTDM(4)
         NBSPFSM(8) = ITOTDM(3)
         NUSPFSM(5) = ITOTDM(1)
         NUSPFSM(6) = ITOTDM(2)
         NUSPFSM(7) = ITOTDM(3)
         NUSPFSM(8) = ITOTDM(4)
*
         do I=1,2,1
            ITRIADB(I) = (SYMDIM(I)**2 + SYMDIM(I)) / 2
            ITRIADB_F(I) = (SYMDIM_F(I)**2 + SYMDIM_F(I)) / 2
         end do
         do I=1,2,1
            ITRIADBSUM_F = ITRIADBSUM_F + ITRIADB_F(I)
         end do
      else if (DOUGRP.eq.5.or.DOUGRP.eq.6) then
*=======
* C2(*) and Cs(*)
*=======
         SYMDIM(1) = ITOTDM(1) + ITOTDM(2)
         SYMDIM_F(1) = NORB_F(1) + NORB_F(2)
         DDIM = (SYMDIM(1))*2
         NBSPFSM(3) = ITOTDM(2)
         NBSPFSM(4) = ITOTDM(1)
         NUSPFSM(3) = ITOTDM(1)
         NUSPFSM(4) = ITOTDM(2)
         ITRIADB(1) = (SYMDIM(1)**2 + SYMDIM(1)) / 2
         ITRIADB_F(1) = (SYMDIM_F(1)**2 + SYMDIM_F(1)) / 2
         ITRIADBSUM_F = ITRIADB_F(1)
      else if (DOUGRP.eq.7) then
*=======
* Ci(*)
*=======
         SYMDIM(1) = ITOTDM(1)*2
         SYMDIM(2) = ITOTDM(2)*2
         SYMDIM_F(1) = NORB_F(1)*2
         SYMDIM_F(2) = NORB_F(2)*2
         DDIM = SYMDIM(1) + SYMDIM(2)
         NBSPFSM(3) = ITOTDM(1)
         NBSPFSM(4) = ITOTDM(2)
         NUSPFSM(3) = ITOTDM(1)
         NUSPFSM(4) = ITOTDM(2)
         ITRIADB(1) = (SYMDIM(1)**2 + SYMDIM(1)) / 2
         ITRIADBSUM = ITRIADB(1)
         ITRIADB_F(1) = (SYMDIM_F(1)**2 + SYMDIM_F(1)) / 2
         ITRIADBSUM_F = ITRIADB_F(1)
      else if (DOUGRP.eq.8) then
*=======
* C1(*)
*=======
         SYMDIM(1) = ITOTDM(1)*2
         SYMDIM_F(1) = NORB_F(1)*2
         NBSPFSM(2) = ITOTDM(1)
         NUSPFSM(2) = ITOTDM(1)
         ITRIADB(1) = (SYMDIM(1)**2 + SYMDIM(1)) / 2
         ITRIADBSUM = ITRIADB(1)
         ITRIADB_F(1) = (SYMDIM_F(1)**2 + SYMDIM_F(1)) / 2
         ITRIADBSUM_F = ITRIADB_F(1)
      else
         write(*,*) ' Sorry  DOUGRP out of range , DOUGRP = ', DOUGRP
         write(*,*) ' SSPIR fatally wounded '
         Call Abend1( 5 )
      end if
*
*  Determine spinor info and symmetry within the double group
*
Cjan      call spiinf(LU1)
*
Cjan      write(*,*) 'This is sspir speaking:'
Cjan      do I=1,DDIM,1
Cjan         write(*,*) 'Spinor ',I,' type is   ',SPINOR(I)
Cjan         write(*,*) 'Spinor symmetry is        ',SPISYM(I)
Cjan      end do
*
      return
      end
*
      SUBROUTINE IAIBCM_KRCC(ICISPC,IATP_TF,IBTP_TF,IAIB)
*
* obtain allowed combination of alpha- and beta- supergroups
* for CI space ICISPC
*
* Master for IAIBCM_GAS
*
*      Jeppe Olsen, august 1995
*
*  Corrected for rel. GASCI
*       Timo Fleig, 1999
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "stinf.inc"
#include "strinp.inc"
*. Output
      INTEGER IAIB(*)
*
      NAEL = NELEC(IATP_TF)
      NBEL = NELEC(IBTP_TF)
*
      NOCTPA = NOCTYP(IATP_TF)
      NOCTPB = NOCTYP(IBTP_TF)
*
      IOCTPA = IBSPGPFTP(IATP_TF)
      IOCTPB = IBSPGPFTP(IBTP_TF)
*
C?    write(6,*) ' IAIB ::::::'
C?    write(6,*) ' LCMBSPC, ICISPC, ICMBSPC '
C?    WRITE(6,*) ICISPC,  LCMBSPC(ICISPC)
C?    WRITE(6,*) (ICMBSPC(II,ICISPC),II=1, LCMBSPC(ICISPC))

      CALL IAIBCM_GAS(LCMBSPC(ICISPC),ICMBSPC(1,ICISPC),
     &                IGSOCCX,NOCTPA,
     &                NOCTPB,ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &                NELFGP,MXPNGAS,NGAS,IAIB,0)
*
      RETURN
      END
*
      SUBROUTINE PART_CIV3_KRCC(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                     IOFFBLK,IOFFBTC,
     &                     NSMST,MXLNG,IOCOC,ISMOST,
     &                     NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP,
     &                     ITTSS_ORD,ISIMSYM)
*
* Partition a CI vector with given MS2 into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1 the complete ci vector is constructed
*
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type , relative to start of supergroup
*   IBATCH(2,*) : Beta type , relative to start of supergroup
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expanded form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, Jan. 1998
*
* Modified from PART_CIV2 for relativistic purposes
*
* Compared to PART_CIV2 the changes are
* 1 : The type numbers are absolute, no need for additional
*     reference to offset for given type
* 2 : Blocks filled from a given offset
*
* ISIMSYM added 2006 Lasse
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) '     PART_CIV3     '
        WRITE(6,*) ' =================='
        WRITE(6,*) ' IDC = ', IDC
        WRITE(6,*)
        WRITE(6,*) ' NSSOA array ( input ) '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        WRITE(6,*) ' NSSOB array ( input ) '
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
        WRITE(6,*) 'ISIMSYM, ICOMP = ',ISIMSYM,ICOMP
      END IF
*
*. block  zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = IOFFBTC-1
      IBLOCK = IOFFBLK-1
      IFINI = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
*. Next block : Order is currently : IB, IA, ISM  (leftmost inner loop )
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
*. Size of TT block ( all symmetries)
      LBLOCK_AS = 0
      IF(ISIMSYM.EQ.1 .AND. ISM. EQ. 1 ) THEN
        DO IASM = 1, NSMST
          IBSM = ISMOST(IASM)
          NSTA = NSSOA(IASM,IA)
          NSTB = NSSOB(IBSM,IB)
          IF(IBLTP(IASM).EQ.0) GOTO 99
          IF(IBLTP(IASM).EQ.2.AND.IA.LT.IB) GOTO 99
          LBLOCK_AS = LBLOCK_AS + NSTA*NSTB
   99   CONTINUE
        END DO
        INCLUDE = 0
C?      WRITE(6,*) ' IA IB LBLOCK_AS', IA,IB, LBLOCK_AS
        IF(LENGTH+LBLOCK_AS.LE.MXLNG.OR.ICOMP.EQ.1) INCLUDE = 1
      END IF
C MOVED IF(IBLTP(ISM).EQ.0) GOTO 1000
C MOVED IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
C MOVED IF(IOCOC(IA,IB).EQ.0) GOTO 1000 with ISIMSYM
C?    write(6,*) ' PART_CIV3 IDC IBLTP ', IDC,IBLTP(ISM)
*. can this block be included
      IBSM = ISMOST(ISM)
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
*
        IF(ISIMSYM.EQ.0) THEN
          INCLUDE = 0
          IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) INCLUDE = 1
        END IF
*
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK
        IF(INCLUDE.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(6,*) ' Not enough scratch space to include a single Block'
        WRITE(6,*) ' Since I cannot procede I will stop '
        WRITE(6,*) ' Insufficient buffer detected in PART_CIV3'
        write(6,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(6,*) ' Alter GAS space of raise Buffer from ', MXLNG
        stop ' In PART_CIV3. '
      ELSE
*. This batch is finished, goto next batch
        GOTO 2000
      END IF
 2002 CONTINUE
*
      IF(NTEST.NE.0) THEN
        WRITE(6,*) 'Output from PART_CIV3'
        WRITE(6,*) '====================='
        WRITE(6,*)
        WRITE(6,*) ' Number of added batches ', NBATCH-IOFFBTC+1
        DO JBATCH = IOFFBTC, NBATCH
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*)
          WRITE(6,*) '      Number of blocks included ', LBATCH(JBATCH)
          WRITE(6,*) '      TTSS and offsets and lenghts of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(6,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE REF_DET_KRCC
*
*     Routine to choose reference determinant based on input
*     This is just setting NELEC(1) and NELEC(2) and
*     the number of electrons in each GAS for the reference
*     Lasse
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "strinp.inc"
#include "lucinp.inc"
#include "crun.inc"
#include "cgas.inc"
*
      NTEST = 0
      NAELRF  = (MK2REF + NACTEL)/2
      NBELRF  = NACTEL - NAELRF
      NELEC(1) = NAELRF
      NELEC(2) = NBELRF
*
      NELECTOT = 0
      DO IGAS=1,NGAS
        IF(IGAS.EQ.1) THEN
          IREF_ELEC(IGAS) = MIN(NACTEL,2*NGSOBT(IGAS)) 
          NELECTOT = NELECTOT + IREF_ELEC(IGAS)
        ELSE IF(NELECTOT.EQ.NACTEL) THEN
          IREF_ELEC(IGAS) = 0
        ELSE
          IRESTELEC = NACTEL - NELECTOT
          IREF_ELEC(IGAS) = MIN(IRESTELEC,2*NGSOBT(IGAS))
          NELECTOT = NELECTOT + IREF_ELEC(IGAS)
        END IF
      END DO
*
      NELECATOT = 0
      NELECBTOT = 0
      DO IGAS=1,NGAS
        IF(IGAS.EQ.1) THEN
          IREF_AB_ELEC(1,IGAS) = MIN(NELEC(1),NGSOBT(IGAS))
          IREF_AB_ELEC(2,IGAS) = MIN(NELEC(2),NGSOBT(IGAS))
          NELECATOT = NELECATOT + IREF_AB_ELEC(1,IGAS)
          NELECBTOT = NELECBTOT + IREF_AB_ELEC(2,IGAS)
        ELSE IF(NELECATOT.EQ.NELEC(1).AND.NELECBTOT.EQ.NELEC(2)) THEN
          IREF_AB_ELEC(1,IGAS) = 0
          IREF_AB_ELEC(2,IGAS) = 0
        ELSE IF(NELECATOT.EQ.NELEC(1)) THEN
          IREF_AB_ELEC(1,IGAS) = 0
          IREF_AB_ELEC(2,IGAS) = MIN(NELEC(2),NGSOBT(IGAS))
          NELECBTOT = NELECBTOT + IREF_AB_ELEC(2,IGAS)
        ELSE IF(NELECBTOT.EQ.NELEC(2)) THEN
          IREF_AB_ELEC(1,IGAS) = MIN(NELEC(1),NGSOBT(IGAS))
          IREF_AB_ELEC(2,IGAS) = 0
          NELECATOT = NELECATOT + IREF_AB_ELEC(1,IGAS)
        ELSE
          IRESTAELEC = NELEC(1) - NELECATOT
          IRESTBELEC = NELEC(2) - NELECBTOT
          IREF_AB_ELEC(1,IGAS) = MIN(IRESTAELEC,NGSOBT(IGAS))
          IREF_AB_ELEC(2,IGAS) = MIN(IRESTBELEC,NGSOBT(IGAS))
          NELECATOT = NELECATOT + IREF_AB_ELEC(1,IGAS)
          NELECBTOT = NELECBTOT + IREF_AB_ELEC(2,IGAS)
        END IF
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info on reference determinant '
        WRITE(6,*) ' Number of active electrons    ',NACTEL
        WRITE(6,*) ' MK2 projection                ',MK2REF
        WRITE(6,*) ' Unbarred number of electrons  ',NELEC(1)
        WRITE(6,*) ' Barred number of electrons    ',NELEC(2)
        WRITE(6,*) '                               '
        WRITE(6,*) '==============================='
        WRITE(6,*) '                               '
        WRITE(6,*) ' Number of electrons in GAS    '
        DO IGAS = 1,NGAS
        WRITE(6,*) ' GAS ',IGAS,' N electrons ',IREF_ELEC(IGAS)
        END DO
        WRITE(6,*) '                               '
        WRITE(6,*) ' Number of unbarred/barred     '
        DO IGAS = 1,NGAS
        WRITE(6,*) ' GAS ',IGAS,' N Unbarred  ',IREF_AB_ELEC(1,IGAS)
        WRITE(6,*) ' GAS ',IGAS,' N Barred    ',IREF_AB_ELEC(2,IGAS)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE CC_AC_SPACES_KRCC(ISPC,IREFTYP)
*
*  Divide orbital spaces  into
*  Hole spaces     : Only annihilation allowed
*  particle spaces : Only creation allowed
*  valence  spaces : Both annihilation and creation allowed
*
* Result is stored in IHPVGAS
*
* Division based upon occupation in CI space ISPC
* Used for Coupled Cluster Calculations
*
* Find type of reference state
*
* IREFTYP = 1 : Closed Shell Hartree-Fock state
* IREFTYP = 2 : Highspin open shell state
* IREFTYP = 3 : More general reference space (CAS, RAS, GAS ..)
*
* Jeppe Olsen, Summer of 98 ( not much of a summer !)
*
*
*
* Major changes for KRCC with only one calc by Lasse 2010
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "strinp.inc"
*
      NEL_REF = NELEC(1) + NELEC(2)
C     WRITE(6,*) ' NELEC(1), NELEC(2) ', NELEC(1),NELEC(2)
*
*. To get rid of annoying and incorrect compiler warnings
      NEL_MAX = 0
*
      NHOLE = 0
      NVAL  = 0
*
* A more recent code without assuming that the hole spaces
* occur first ( for the QD project)
      NELECTOT = 0
      DO IGAS = 1, NGAS
*. Minimum number of electrons in this space
C       IF(IGAS.EQ.1) THEN
C         NEL_MIN = IGSOCCX(1,1,ISPC)
C       ELSE
C         NEL_MIN = IGSOCCX(IGAS,1,ISPC)-IGSOCCX(IGAS-1,2,ISPC)
C         NEL_MIN = MAX(0,NEL_MIN)
C       END IF
C. Largest number of electrons in this space
C       IF(IGAS.EQ.1) THEN
C         NEL_MAX = IGSOCCX(1,2,ISPC)
C       ELSE
C         NEL_MAX = IGSOCCX(IGAS,2,ISPC)-IGSOCCX(IGAS-1,1,ISPC)
C         print*,'IGSOCCX(IGAS,2,ISPC),IGSOCCX(IGAS-1,1,ISPC)',
C    &            IGSOCCX(IGAS,2,ISPC),IGSOCCX(IGAS-1,1,ISPC)
C         NEL_MAX = MIN(NEL_MAX,2*NGSOBT(IGAS))
C       END IF
C       WRITE(6,*) ' IGAS, NEL_MAX, NEL_MIN = ',
C    &               IGAS, NEL_MAX, NEL_MIN
C       IF(NEL_MAX.EQ.0) THEN
C. particle space
C        IHPVGAS(IGAS) = 2
C       ELSE IF (NEL_MIN.EQ.2*NGSOBT(IGAS)) THEN
C. Hole space
C        IHPVGAS(IGAS) = 1
C        NHOLE = NHOLE + NGSOBT(IGAS)
C       ELSE
C. Valence space
C        IHPVGAS(IGAS) = 3
C        NVAL = NVAL + NGSOBT(IGAS)
C       END IF
        IF(IREF_ELEC(IGAS).EQ.2*NGSOBT(IGAS)) THEN
C. Hole space
          IHPVGAS(IGAS) = 1
          NHOLE = NHOLE + NGSOBT(IGAS)
        ELSE IF(IREF_ELEC(IGAS).EQ.0) THEN
C. particle space
          IHPVGAS(IGAS) = 2
        ELSE
C. Valence space
          IHPVGAS(IGAS) = 3
          NVAL = NVAL + NGSOBT(IGAS)
        END IF
      END DO
*
      NEL_AL = NELEC(1)
      NEL_BE = NELEC(2)
C      print*,'NEL_AL,NEL_BE',NEL_AL,NEL_BE
      IF(NEL_AL.EQ.NHOLE.AND.NEL_BE.EQ.NHOLE) THEN
*. Closed shell Hartree-Fock
        IREFTYP = 1
      ELSE IF(NEL_AL.EQ.NHOLE.AND.NEL_BE.EQ.NHOLE+NVAL.OR.
     &        NEL_BE.EQ.NHOLE.AND.NEL_AL.EQ.NHOLE+NVAL) THEN
*. Highspin openshell
        IREFTYP = 2
      ELSE
*. More general, not analyzed in detail p.t.
        IREFTYP = 3
      END IF
*
*. IHPVGAS_AB : Hole, particle, valence for alpha and beta orbitals
*
*. Differs from IHPVGAS for high spin open shell case
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.1) THEN
          IHPVGAS_AB(IGAS,1) = 1
          IHPVGAS_AB(IGAS,2) = 1
        ELSE IF (IHPVGAS(IGAS).EQ.2) THEN
          IHPVGAS_AB(IGAS,1) = 2
          IHPVGAS_AB(IGAS,2) = 2
        ELSE IF (IHPVGAS(IGAS).EQ.3) THEN
          IF(IREFTYP.EQ.2.AND.NEL_AL.GT.NEL_BE) THEN
            IHPVGAS_AB(IGAS,1) = 1
            IHPVGAS_AB(IGAS,2) = 2
          ELSE IF(IREFTYP.EQ.2.AND.NEL_AL.LT.NEL_BE) THEN
            IHPVGAS_AB(IGAS,1) = 2
            IHPVGAS_AB(IGAS,2) = 1
          ELSE
            IHPVGAS_AB(IGAS,1) = 3
            IHPVGAS_AB(IGAS,2) = 3
          END IF
        END IF
      END DO
*

      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CC division of orbitals '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
        WRITE(6,*) ' Hole =>1, Part=>2, Val=>3 '
        WRITE(6,*)
        CALL IWRTMA(IHPVGAS,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' CC division for alpha-spinorbitals'
        CALL IWRTMA(IHPVGAS_AB(1,1),1,NGAS,1,NGAS)
*
        WRITE(6,*) ' CC division for beta-spinorbitals'
        CALL IWRTMA(IHPVGAS_AB(1,2),1,NGAS,1,NGAS)
*
        IF(IREFTYP.EQ.1) THEN
          WRITE(6,*) ' Reference state is closed shell single SD'
        ELSE IF( IREFTYP.EQ.2) THEN
          WRITE(6,*)
     &    ' Reference state is high spin open shell single SD'
        ELSE
          WRITE(6,*)
     &    ' Reference state is general multireference state '
        END IF

      END IF
*
      RETURN
      END
*
      SUBROUTINE ACAC_EXC_TYP_KRCC(IAAEXC,MX_AAEXC)
*
* Information about active-active excitations
*
* IAAEXC_TYP = 0 => No active-active excitations (Closed shell HF)
*              1 => Active alpha => Active beta (High-spin, Ms = Max)
*              2 => Active beta  => Active alpha(High-spin, Ms = Min)
*              3 => All types of active-active   (General CAS)
*
* Identification of double occupied/valence orbitals is based on
* info in IHPVGAS
*
* Used for Coupled Cluster Calculations
*
* Find type of reference state
*
* Jeppe Olsen, March of 2000
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "strinp.inc"
*
      NTEST = 00
*
*. Number of hole and valence orbitals
      NHOLE = 0
      NVAL  = 0
      DO IGAS = 1, NGAS
*
       IF(IHPVGAS(IGAS).EQ.1) THEN
*. hole space
          NHOLE = NHOLE + NGSOBT(IGAS)
       ELSE IF(IHPVGAS(IGAS).EQ.3) THEN
*. Valence space
          NVAL = NVAL + NGSOBT(IGAS)
       END IF
      END DO
*
      NEL_AL = NELEC(1)
      NEL_BE = NELEC(2)
C Add Lasse
C Make sure highest excitation does not surpass number of electrons
      IF(NHOLE+NVAL.GT.NEL_AL+NEL_BE) THEN
        NVAL = NEL_AL+NEL_BE-NHOLE
      END IF
C End Lasse
      IF(NEL_AL.EQ.NHOLE.AND.NEL_BE.EQ.NHOLE) THEN
*. Closed shell Hartree-Fock
        IAAEXC = 0
        MX_AAEXC = 0
      ELSE IF(NEL_BE.EQ.NHOLE.AND.NEL_AL.EQ.NHOLE+NVAL) THEN
*. High spin open shell case with Max MS
        IAAEXC = 1
        MX_AAEXC = NVAL
      ELSE IF(NEL_AL.EQ.NHOLE.AND.NEL_BE.EQ.NHOLE+NVAL ) THEN
*. High spin open shell case with Min MS
        IAAEXC = 2
        MX_AAEXC = NVAL
      ELSE
*. More general, not analyzed in detail p.t. , assumed CAS
        IAAEXC = 3
        MX_AAEXC = NVAL
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) ' Allowed types of active-active excitations '
        WRITE(6,*) ' ========================================== '
        WRITE(6,*)
        IF(IAAEXC.EQ.0) THEN
          WRITE(6,*) ' No excitations in active orbital space '
        ELSE IF(IAAEXC.EQ.1) THEN
          WRITE(6,*) ' Alpha => Beta active excitations '
        ELSE IF (IAAEXC.EQ.2) THEN
          WRITE(6,*) ' Beta => alpha active excitations '
        ELSE IF (IAAEXC.EQ.3) THEN
          WRITE(6,*) ' All active-active excitations '
        END IF
*
        WRITE(6,'(A,I5)')
     &  ' Largest excitation level for active-active excitations',
     &    MX_AAEXC
      END IF
*
      RETURN
      END
*
      SUBROUTINE OCCLSE_REF(IOCCLS)
*
* Lasse I think this works
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "strinp.inc"
*
      INTEGER :: IOCCLS(NGAS)
*
      DO IGAS=1,NGAS
        IOCCLS(IGAS) = IREF_ELEC(IGAS)
      END DO
*
      RETURN
      END
*
      SUBROUTINE OCCLSE_KRCC(IWAY,NOCCLS,IOCCLS,NEL,ICISPC,
     &                  I_DO_BASSPC,IBASSPC,NOBPT)
*
* IWAY = 1 :
* obtain NOCCLS =
* Number of allowed ways of distributing the orbitals in
* CI space ICISPC
* active spaces
*
* IWAY = 2 :
* OBTAIN NOCCLS and
* IOCCLS = the allowed distributions of electrons
*
* Extended OCCLS : Allows the use of CI spaces obtained
*                  as combinations of occupation spaces
*
*
*
* Jeppe Olsen, October 2000
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "cgas.inc"
*. Input
      DIMENSION NOBPT(NGAS)
*. Output
      DIMENSION IOCCLS(NGAS,*)
      DIMENSION IBASSPC(*)
*. Local scratch
      INTEGER IOCMIN(MXPNGAS),IOCMAX(MXPNGAS),IREFOCCLSE(MXPNGAS)
      DIMENSION IOCA(MXPNGAS),IOC(MXPNGAS)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
         WRITE(6,*)  ' OCCLS in action '
         WRITE(6,*) ' =================='
         WRITE(6,*) ' ICISPC = ', ICISPC
      END IF
*. Largest and smallest accumulated occupations
      DO IGAS = 1, NGAS
        KMIN = IGSOCCX(IGAS,1,ICMBSPC(1,ICISPC))
        KMAX = IGSOCCX(IGAS,2,ICMBSPC(1,ICISPC))
        DO JSPC = 2, LCMBSPC(ICISPC)
          LSPC = ICMBSPC(JSPC,ICISPC)
          KMIN = MIN(KMIN, IGSOCCX(IGAS,1,LSPC))
          KMAX = MAX(KMAX, IGSOCCX(IGAS,2,LSPC))
        END DO
        IOCMIN(IGAS) = KMIN
        IOCMAX(IGAS) = KMAX
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' IOCMIN and IOCMAX from OCCLSE '
        CALL IWRTMA(IOCMIN,1,NGAS,1,NGAS)
        CALL IWRTMA(IOCMAX,1,NGAS,1,NGAS)
      END IF
*
      NOCCLS = 0
*. start with smallest allowed number
      DO IGAS = 1, NGAS
*. Smallest allowed occ in this GASspace
        IOCA(IGAS) =  IOCMIN(IGAS)
      END DO
      NONEW = 0
      IFIRST = 1
*. Loop over possible occupations
 1000 CONTINUE
        IF(IFIRST.EQ.0) THEN
*. Next accumulated occupation
          CALL NXTNUM3(IOCA,NGAS,IOCMIN,IOCMAX,NONEW)
        END IF
        IFIRST = 0
        IF(NONEW.EQ.0) THEN
*. ensure that IOCA corresponds to an accumulating occupation,
*. i.e. a non-decreasing sequence
          KGAS = 0
          DO IGAS = 2, NGAS
            IF(IOCA(IGAS-1).GT.IOCA(IGAS)) KGAS = IGAS
          END DO
* Test below looks odd? see line 59 as well Lasse
* this does not appear to be enough to ensure a decreasing
* sequence except if NGAS=1,2. But we will not make silly inputs
          IF(KGAS .NE. 0 ) THEN
            DO IGAS = 1, KGAS-1
              IOCA(IGAS) = IOCMIN(IGAS)
            END DO
            IOCA(KGAS) = IOCA(KGAS)+1
          END IF
C?      WRITE(6,*) ' Another accumulated occupation: '
C?      CALL IWRTMA(IOCA,1,NGAS,1,NGAS)
*. corresponding occupation of each active space
        NEGA=0
        IM_TO_STUFFED = 0
        DO IGAS = 1, NGAS
          IF(IGAS.EQ.1) THEN
            IOC(IGAS) = IOCA(IGAS)
          ELSE
            IOC(IGAS) = IOCA(IGAS)-IOCA(IGAS-1)
            IF(IOC(IGAS).LT.0) NEGA = 1
            IF(IOC(IGAS).GT.2*NOBPT(IGAS)) IM_TO_STUFFED = 1
          END IF
        END DO
C?      WRITE(6,*) ' Another occupation: '
C?      CALL IWRTMA(IOC,1,NGAS,1,NGAS)
*. Correct number of electrons ?
        IEL_TOT = IELSUM(IOC,NGAS)
*. Belongs to one of the accumulated spaces ?
        IN_SPC = 0
        DO JSPC = 1,  LCMBSPC(ICISPC)
          JCISPC = ICMBSPC(JSPC,ICISPC)
C?        WRITE(6,*) ' JSPC, JCISPC = ', JSPC, JCISPC
          INSPC_LOC = 1
          DO IGAS = 1, NGAS
            IEL = IOCA(IGAS)
            IF(IEL.LT. IGSOCCX(IGAS,1,JCISPC).OR.
     &         IEL.GT. IGSOCCX(IGAS,2,JCISPC)    ) INSPC_LOC = 0
          END DO
C?        WRITE(6,*) ' JSPC, INSPC_LOC =', JSPC, INSPC_LOC
          IF(INSPC_LOC.EQ.1) IN_SPC = 1
        END DO
*
C?      WRITE(6,*) ' IEL_TOT, NEGA, IN_SPC = ',
C?   &               IEL_TOT, NEGA, IN_SPC
        IF(IEL_TOT.EQ.NEL.AND.NEGA.EQ.0.AND.IM_TO_STUFFED.EQ.0.AND.
     &     IN_SPC.EQ.1) THEN
C
C This is a fast routine so check here if this occupation is mixing core and valence
C At the moment only one core is assumed! Lasse
C
C A possiblitity to do mixed core-valence in cas space is now allowed
C Assumption at the moment is core space 1 and cas space 3
C
          IF(INOCOVA.EQ.0) THEN
            NOCCLS = NOCCLS + 1
            IF(IWAY.EQ.2) THEN
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) ' Another allowed class : '
                CALL IWRTMA(IOC,1,NGAS,1,NGAS)
              END IF
              CALL ICOPVE(IOC,IOCCLS(1,NOCCLS),NGAS)
*
              IF(I_DO_BASSPC.EQ.1) THEN
              IBASSPC(NOCCLS) = IBASSPC_FOR_CLS(IOC)
              END IF
*
            END IF
          ELSE IF(INOCOVA.EQ.1) THEN !Do the check
C Find max mixed core valence excitation
            IMAXMIX = IOCMAX(1) - IOCMIN(1)
C Find reference OCCLSE
            IREFOCCLSE(1) = IOCMAX(1)
            DO IGAS= 2,NGAS-1
              IREFOCCLSE(IGAS) = IOCMAX(IGAS) - IOCMAX(IGAS-1)
            END DO
            IREFOCCLSE(NGAS) = 0
C Find current mixed excitation
            IMIX = IREFOCCLSE(1) - IOC(1)
            IF(IMIX.GE.1) THEN
              DO IGAS = 2,NGAS
                IF(IREFOCCLSE(IGAS).EQ.0) THEN
                  EXIT
                ELSE
                  IMIX = IMIX + IREFOCCLSE(IGAS) - IOC(IGAS)
                END IF
              END DO
            ELSE
              IMIX = 0
            END IF
C           print*,'IMIX',IMIX
            IF(IMAXMIX.GE.IMIX) THEN
              IPASS = 1
            ELSE
              IPASS = 0
            END IF
            IF(IPASS.EQ.1) THEN
              NOCCLS = NOCCLS + 1
              IF(IWAY.EQ.2) THEN
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' Another allowed class : '
                  CALL IWRTMA(IOC,1,NGAS,1,NGAS)
                END IF
                CALL ICOPVE(IOC,IOCCLS(1,NOCCLS),NGAS)
*
                IF(I_DO_BASSPC.EQ.1) THEN
                  IBASSPC(NOCCLS) = IBASSPC_FOR_CLS(IOC)
                END IF
*
              END IF
            END IF
          ELSE
C Here the mixed core-valence correlation in cas comes in
C Find max mixed core valence excitation
            IMAXMIX  = IOCMAX(1) - IOCMIN(1)
            IMAXMIX2 = IOCMAX(3) - IOCMIN(3)
C Find reference OCCLSE
            IREFOCCLSE(1) = IOCMAX(1)
            DO IGAS= 2,NGAS-1
              IREFOCCLSE(IGAS) = IOCMAX(IGAS) - IOCMAX(IGAS-1)
            END DO
            IREFOCCLSE(NGAS) = 0
C Find current mixed excitation
            IMIX = IREFOCCLSE(1) - IOC(1)
            IF(IMIX.GE.1) THEN
              DO IGAS = 2,NGAS
                IF(IREFOCCLSE(IGAS).EQ.0) THEN
                  IF(IMIX.GT.IMAXMIX) THEN
C Here the 3. space is skipped in counting for IMIX
                    DO IIGAS = IGAS+1,NGAS
                      IF(IOC(IIGAS).NE.0) THEN
C In case of 4+ spaces (and then virtuals are split)
                        IMIX = IMAXMIX + 1
                        EXIT
                      ELSE
C This will be an allowed excitation
                        IMIX = IMAXMIX
                      END IF
                    END DO
                  ELSE
                    EXIT
                  END IF
                ELSE
                  IMIX = IMIX + IREFOCCLSE(IGAS) - IOC(IGAS)
                END IF
              END DO
            ELSE
              IMIX = 0
            END IF
C           print*,'IMIX',IMIX
            IF(IMAXMIX.GE.IMIX) THEN
              IPASS = 1
            ELSE
              IPASS = 0
            END IF
            IF(IPASS.EQ.1) THEN
              NOCCLS = NOCCLS + 1
              IF(IWAY.EQ.2) THEN
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' Another allowed class : '
                  CALL IWRTMA(IOC,1,NGAS,1,NGAS)
                END IF
                CALL ICOPVE(IOC,IOCCLS(1,NOCCLS),NGAS)
*
                IF(I_DO_BASSPC.EQ.1) THEN
                  IBASSPC(NOCCLS) = IBASSPC_FOR_CLS(IOC)
                END IF
*
              END IF
            END IF
          END IF
        END IF
      END IF
      IF(NONEW.EQ.0) GOTO 1000
*
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Number of Allowed occupation classes ', NOCCLS
         IF(IWAY.EQ.2.AND.NTEST.GE.20) THEN
           WRITE(6,*) ' Occupation classes : '
           WRITE(6,*) ' ===================='
           WRITE(6,*)
           WRITE(6,*) ' Class    Occupation in GASpaces '
           WRITE(6,*) ' ================================'
           DO I = 1, NOCCLS
             WRITE(6,'(1H ,I5,3X,16I3)')
     &       I, (IOCCLS(IGAS,I),IGAS=1, NGAS)
           END DO
C          CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
         END IF
      END IF
*
      RETURN
      END
*
      SUBROUTINE TP_OBEX2_KRCC(NOCCLS,NEL,NGASX,IOBEX_TP,LCOBEX_TP,
     &                   LAOBEX_TP,
     &                   IOCCLS,IOCCLS_REF,MX_NCREA,MX_NANNI,
     &                   MX_EXC_LEVEL,IEXTP_TO_OCCLS,MX_AAEXC,IFLAG,
     &                   NOBEX_TP)
*
* Obtain the orbital excitation types needed to generate occupation classes
* in IOCCLS from IOCCLS_REF
*
* Jeppe Olsen, Updated version of TP_OBEX, March 2000, Still on the train
*              April 2001         unit operator added as excitation
*                                 NOBEX_TP + 1. The arrays
*                                 IOBEX_TP, LCOBEX, LAOBEX, IEXTP_TO_OCCLS
*                                 should this be dimensioned with NOBEX_TP+1
*
*
* Allows active-active excitations. Active orbitals are assumed to
* be in a single orbital space,
*
* If IFLAG = 1, only the number of orbital excitation types is generated
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*. Input
      INTEGER IOCCLS(NGAS,*), IOCCLS_REF(NGAS)
*. Output
      INTEGER IOBEX_TP(2*NGAS,*),  LCOBEX_TP(*), LAOBEX_TP(*)
      INTEGER IEXTP_TO_OCCLS(*)
*. Local scratch
      DIMENSION ICREA(MXPNGAS),IANNI(MXPNGAS)
*. Number of active orbital spaces
      NACT_SPC = 0
      IACT_SPC = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.3) THEN
          NACT_SPC = NACT_SPC + 1
          IACT_SPC = IGAS
        END IF
      END DO
*
* This have been checked several places, Lasse
*
      IF(NACT_SPC.GT.1) THEN
        WRITE(6,*) ' TP_OBEX2 in problems '
        WRITE(6,*) ' More than one active orbital spaces '
        WRITE(6,*) ' NACT_SPC = ',  NACT_SPC
        STOP ' TP_OBEX2 :  More than one active orbital spaces '
      END IF
C?    WRITE(6,*) ' TP_OBEX2 : IACT_SPC,NACT_SPC',IACT_SPC,NACT_SPC
C?    WRITE(6,*) ' TP_OBEX2 : MX_AAEXC ', MX_AAEXC
*
* The orbital excitation operator IEXTP is  organized as
*
* LCOBEX(IEXTP) : Number of creation operators
* LAOBEX(IEXTP) : Number of annihilation operators
* IOBEX_TP(1 - NGAS, IEXTP) : Number of creation operators per gassspace
* IOBEX(NGAS+1  -  2*NGAS, , IEXTP) : Number of annihilation operators
*
* IEXTP_TO_OCCLS is map from orbital excitation type to occupation
* class for CI coefficients
      NTEST = 00
      MX_NCREA = 0
      MX_NANNI = 0
      NOBEX_TP = 0
      MX_EXC_LEVEL = 0
      JREFCLS = 0
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'TP_OBEX speaking '
        WRITE(6,*) ' Reference occupation class '
        CALL IWRTMA(IOCCLS_REF,1,NGAS,1,NGAS)
      END IF
*
      IZERO = 0
      DO JOCCLS = 1, NOCCLS
        NANNI = 0
        NCREA = 0
        CALL ISETVC(ICREA,IZERO,NGAS)
        CALL ISETVC(IANNI,IZERO,NGAS)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Excited occupation class '
          CALL IWRTMA(IOCCLS(1,JOCCLS),1,NGAS,1,NGAS)
        END IF
        DO IGAS = 1, NGAS
          IF(IOCCLS(IGAS,JOCCLS).GT.IOCCLS_REF(IGAS)) THEN
            ICREA(IGAS) = IOCCLS(IGAS,JOCCLS) - IOCCLS_REF(IGAS)
            NCREA = NCREA + ICREA(IGAS)
          ELSE IF (IOCCLS(IGAS,JOCCLS).LT.IOCCLS_REF(IGAS)) THEN
            IANNI(IGAS) = -( IOCCLS(IGAS,JOCCLS) - IOCCLS_REF(IGAS))
            NANNI = NANNI + IANNI(IGAS)
          END IF
        END DO
        IF(NANNI.EQ.0.AND.NCREA.EQ.0) JREFCLS = JOCCLS
*. Add active -active excitation
        DO IAA_EXC = 0, MX_AAEXC
         ITSOKAY = 1
         IF(IAA_EXC.GT.0) THEN
*. Can another active-active excitation be added ?
          IF(IANNI(IACT_SPC)+1.LE.MX_AAEXC.AND.
     &       ICREA(IACT_SPC)+1.LE.MX_AAEXC) THEN
              ITSOKAY = 1
          ELSE
              ITSOKAY = 0
          END IF
         END IF
         IF(ITSOKAY.EQ.1) THEN
*
          IF(IAA_EXC.GT.0) THEN
            ICREA(IACT_SPC)  = ICREA(IACT_SPC) + 1
            IANNI(IACT_SPC)  = IANNI(IACT_SPC) + 1
            NCREA = NCREA + 1
            NANNI = NANNI + 1
          END IF
*
          MX_EXC_LEVEL = MAX(MX_EXC_LEVEL,NCREA)
          IF(NCREA+NANNI.NE.0) THEN
            NOBEX_TP = NOBEX_TP + 1
            IF(IFLAG.NE.1) THEN
              LCOBEX_TP(NOBEX_TP) = NCREA
              LAOBEX_TP(NOBEX_TP) = NANNI
C             print*,'NCREA,NANNI',NCREA,NANNI
              IEXTP_TO_OCCLS(NOBEX_TP) = JOCCLS
              CALL ICOPVE(ICREA,IOBEX_TP(1,NOBEX_TP),NGAS )
              CALL ICOPVE(IANNI,IOBEX_TP(NGAS+1,NOBEX_TP),NGAS )
            END IF
          END IF
*
         END IF
        END DO
*       ^ End of loop over active-active excitations
      END DO
*
      IF(IFLAG.NE.1) THEN
*. Add unit operator as excition NOBEX_TP + 1
        LCOBEX_TP(NOBEX_TP + 1) = 0
        LAOBEX_TP(NOBEX_TP + 1) = 0
        IEXTP_TO_OCCLS(NOBEX_TP+1) = JREFCLS
        IZERO = 0
        CALL ISETVC(IOBEX_TP(1,NOBEX_TP+1),IZERO,NGAS)
        CALL ISETVC(IOBEX_TP(NGAS+1,NOBEX_TP+1),IZERO,NGAS)
      END IF
*
      IF(NTEST.GE.3) THEN
        WRITE(6,*) ' Largest excitation level : ', MX_EXC_LEVEL
        WRITE(6,*)
        WRITE(6,*) ' Number of types of orbital excitations ', NOBEX_TP
        WRITE(6,*)
*
        IF(IFLAG.NE.1) THEN
          WRITE(6,*) ' Creation part,  Annihilation  part '
          WRITE(6,*) ' ==================================='
          DO IOBEX = 1, NOBEX_TP+1
            WRITE(6,'(16I4,16I4)')
     &      (IOBEX_TP(I,IOBEX),I=1, NGAS),
     &      (IOBEX_TP(NGAS+I,IOBEX),I=1, NGAS)
          END DO
*
          WRITE(6,*) ' Orbital excitation type to occupation class '
          CALL IWRTMA(IEXTP_TO_OCCLS,1,NOBEX_TP+1,1,NOBEX_TP+1)
        END IF
*
      END IF
*
      RETURN
      END
*
      SUBROUTINE OBEX_TO_SPOBEX_KRCC(IFLAG,IOBEX_TP,LCOBEX_TP,LAOBEX_TP,
     &                          NOBEX_TP,ISPOBEX_TP,NSPOBEX_TP,NGAS,
     &                          NOBPT,MS2TOT,IAAEXC,IACT_SPC,
     &                          ISOX_TO_OX,MXSPOX,NSOX_FOR_OX,
     &                          IBSOX_FOR_OX,ISOX_FOR_OX,NAEL,NBEL,
     &                          IREFSPC,IRECOM,WORK,KFREE,LFREE)
*
* Orbital excitation types => Spin-orbital excitations
*
* IFLAG = 1 : Just Number of spinorbital excitations
* IFLAG = 2 : Also the actual excitations
*
*. Jeppe Olsen, Summer of 99
*               Updated with Active-Active excitations, March 2000
*               SOX<=>OX mapping added, Spring 2001
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cstate.inc"
*
      DIMENSION WORK(*)
*. Input
      INTEGER IOBEX_TP(2*NGAS, NOBEX_TP), NOBPT(*)
      INTEGER LCOBEX_TP(NOBEX_TP), LAOBEX_TP(NOBEX_TP)
*. Output : Creation in alpha, Creation in beta,
*           Annihilation in alpha, annihilation in beta
      INTEGER ISPOBEX_TP(4*NGAS,*)
*. Spin orbital type => orbital type
      INTEGER ISOX_TO_OX(*)
*. Number of spinorbital excitations for given orbital excitation
      INTEGER NSOX_FOR_OX(*)
*. And the actual spinorbital excitations for given orbital excitation
      INTEGER ISOX_FOR_OX(*)
*. Offset to given orbital type in ISOX_FOR_OX
      INTEGER IBSOX_FOR_OX(*)
*. Local scratch for occ of reference alpha and beta
      INTEGER IREF_AL(MXPNGAS),IREF_BE(MXPNGAS)
*
      NTEST = 00
*
      IDUM = 0
      CALL MEMCHK_KRCC(WORK)
*. Largest number of creation and annihilation operators
      MX_CREA = IMNMX(LCOBEX_TP,NOBEX_TP,2)
      MX_ANNI = IMNMX(LAOBEX_TP,NOBEX_TP,2)
*. Largest possible block of spincombination
*. Nup + Ndown = MX_OPS, Nup - Ndown = Ms2_TOT
*. Ms2_TOT needs changing from in- or outside .Lasse
      MX_OPS = MX_CREA + MX_ANNI
      MX_UP = (MX_OPS+MS2TOT)/2
      MX_OPSD2 = MX_OPS/2
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NOBEX_TP = ', NOBEX_TP
        WRITE(6,*) ' MX_CREA, MX_ANNI = ', MX_CREA,MX_ANNI
        WRITE(6,*) ' MX_OPS, MX_UP = ',  MX_OPS, MX_UP
        WRITE(6,*) ' MS2TOT, NAEL, NBEL = ', MS2TOT, NAEL, NBEL
        WRITE(6,*) ' IREFSPC = ',IREFSPC
        WRITE(6,*) ' NMS2VAL = ',NMS2VAL
      END IF
*
C     MX_SOEX_BLK = IBION(MX_OPS,MX_OPSD2)
Cori  MX_SOEX_BLK = IBION(MX_OPS,MX_UP)
* all MS2 values needed. Lasse
      MX_SOEX_BLK = 0
      DO I = 1, NMS2VAL
        TEMP = ABS(MS2VAL(I))
        MX_SOEX_BLK_TEMP = IBION(MX_OPS/2,TEMP)
        MX_SOEX_BLK = MX_SOEX_BLK + MX_SOEX_BLK_TEMP
      END DO
Cori      CALL MEMMAR(KLSOEX ,MX_SOEX_BLK*MX_OPS,'ADDL  ',2,'SOEXBL')
C
C
C     CHECK IF THIS IS OK FOR OPEN SHELL!!!!!
C
C
      CALL MEMGET('REAL',KLSOEX,MX_SOEX_BLK*MX_OPS*NMS2VAL*2*40,
     &            WORK,KFREE,LFREE)
      CALL MEMCHK_KRCC(WORK)
*.Alpha and beta-occupations for reference space
*. Notice : IREFSPC is not used, info in IHPVGAS is used !
*. if condition in use here
*. since IREFSPC=1 is hardwired in lucia_gencc
      IF(IREFSPC.NE.0) THEN
        CALL GET_REF_ALBE_OCC_KRCC(IREFSPC,IREF_AL,IREF_BE,IMS2_REF)
      END IF
*.
      NSPOBEX_TP = 0
      DO ICREA = 0, MX_CREA
        DO IANNI = 0, MX_ANNI
*. Any operators with this number  of creation and annihilation operators ?
          NCOMP = 0
          DO JOBEX_TP = 1, NOBEX_TP
            IF(LCOBEX_TP(JOBEX_TP).EQ.ICREA .AND.
     &         LAOBEX_TP(JOBEX_TP).EQ.IANNI ) NCOMP = NCOMP + 1
C            print*,'LCOBEX_TP,LAOBEX_TP,JOBEX_TP,ICREA,IANNI',
C     &      LCOBEX_TP(JOBEX_TP),LAOBEX_TP(JOBEX_TP),JOBEX_TP,ICREA,IANNI
          END DO
*
C          print*,'NCOMP',NCOMP
C          NCOMP = 0
          IF(NCOMP.NE.0) THEN
*. ICREA creation ops and IANNI annihilation strings with total
*. spin projection MSTOT
*. Changes needs to be made from here
            DO I = 1,NMS2VAL
C              IMS2TOT = MS2VAL(I)
C              CALL AB_COMP_FOR_OBOP_REL(ICREA,IANNI,IMS2TOT,IMS2_REF,
C     &                 NABCOMP,WORK(KLSOEX),NAEL,NBEL,IFLAG)
              IMS2TOT = MS2VAL(I)
C              CALL EVENODD(ITEST,IMS2TOT)
C              IF(ITEST.EQ.1) CYCLE
C             print*,'IMS2TOT',IMS2TOT
C              if(j.eq.1) then
              IMS2TOT = IMS2TOT - IMS2_REF
C              else
C              IMS2TOT = IMS2TOT - IMS2_REF
C Notice following line need to be looked at in restricted/adapted case
C              IALBES = IREF_AL
C              IREF_AL = IREF_BE
C              IREF_BE = IALBES
C              end if
*
* Notice IMS2_REF is not used in AB_COMP_FOR_OBOP_REL since added on outside
*
              CALL AB_COMP_FOR_OBOP_KRCC(ICREA,IANNI,IMS2TOT,IMS2_REF,
     &                 NABCOMP,WORK(KLSOEX),NAEL,NBEL,IFLAG)
              IF(NTEST.GE.100) THEN
              WRITE(6,*) ' NABCOMP = ', NABCOMP
              WRITE(6,*) ' IMS2TOT = ', IMS2TOT
              WRITE(6,*) ' IMS2_REF = ', IMS2_REF
              END IF
              DO JOBEX_TP = 1, NOBEX_TP
                IF(LCOBEX_TP(JOBEX_TP).EQ.ICREA .AND.
     &             LAOBEX_TP(JOBEX_TP).EQ.IANNI .AND.
     &             (MXSPOX.EQ.0.OR.ICREA.LE.MXSPOX)  ) THEN
*. Find the number of the above spincombinations that can be included
*  for this orbital excitaion
                   NOPEN = ICREA + IANNI
*
                   CALL ACT_SPOBEX_OBEX_KRCC(IFLAG,
     &                IOBEX_TP(1,JOBEX_TP),
     &                NGAS,NOPEN,NABCOMP,WORK(KLSOEX),NSPCOMP_ACT,
     &                ISPOBEX_TP(1,NSPOBEX_TP+1),NOBPT,
     &                IACT_SPC,IAAEXC,IREF_AL,IREF_BE,
     &                IREFSPC,NAEL,NBEL,IRECOM)
*
                   IF(IFLAG.EQ.2) THEN
C?    WRITE(6,*) ' First SPOBEX_TP (a) '
C?    CALL IWRTMA(ISPOBEX_TP(1,1),4*NGAS,1,4*NGAS)
                     IOFF = NSPOBEX_TP+1
C                     WRITE(6,*) ' IOFF, NSPCOMP_ACT = ',
C     &                            IOFF, NSPCOMP_ACT
                     CALL ISETVC(ISOX_TO_OX(IOFF),JOBEX_TP,
     &                           NSPCOMP_ACT)
C                     print*,'ISOX_TO_OX,IOFF',ISOX_TO_OX(IOFF),IOFF
C?    WRITE(6,*) ' First SPOBEX_TP (b) '
C?    CALL IWRTMA(ISPOBEX_TP(1,1),4*NGAS,1,4*NGAS)
                   END IF
*
                   NSPOBEX_TP =  NSPOBEX_TP + NSPCOMP_ACT
                END IF
              END DO
            END DO
          END IF
*         ^ End if nonvanishing number of operators
        END DO
      END DO
*     ^ End of loop over creation and annihilation operators
C?    WRITE(6,*) ' Memchk after construction of spobex '
C?    WRITE(6,*) ' Memcheck passed '
*
C?    WRITE(6,*) ' First SPOBEX_TP'
C?    CALL IWRTMA(ISPOBEX_TP(1,1),4*NGAS,1,4*NGAS)
*
      IF(IFLAG.EQ.2) THEN
*. OX => SOX lists
*. Number of sox per ox
        IZERO = 0
        CALL ISETVC(NSOX_FOR_OX,IZERO,NOBEX_TP)
        DO ISOX = 1, NSPOBEX_TP
          IOX = ISOX_TO_OX(ISOX)
          NSOX_FOR_OX(IOX) = NSOX_FOR_OX(IOX)+1
        END DO
*. offsets for sox's with given ox
        DO IOX = 1, NOBEX_TP
         IF(IOX.EQ.1) THEN
           IBSOX_FOR_OX(IOX) = 1
         ELSE
           IBSOX_FOR_OX(IOX) = IBSOX_FOR_OX(IOX-1) + NSOX_FOR_OX(IOX-1)
         END IF
        END DO
*. and the sox for a given ox
        CALL ISETVC(NSOX_FOR_OX,IZERO,NOBEX_TP)
        DO ISOX = 1, NSPOBEX_TP
          IOX = ISOX_TO_OX(ISOX)
          NSOX_FOR_OX(IOX) = NSOX_FOR_OX(IOX)+1
          IADR = IBSOX_FOR_OX(IOX) -1 + NSOX_FOR_OX(IOX)
          ISOX_FOR_OX(IADR) = ISOX
        END DO
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' ***************************************** '
        WRITE(6,*) ' Information about spinorbital excitations '
        WRITE(6,*) ' ***************************************** '
        WRITE(6,*)
        WRITE(6,*) ' Total number ', NSPOBEX_TP
        IF(IFLAG.EQ.2) THEN
          DO JSPCOMP_ACT = 1, NSPOBEX_TP
            WRITE(6,*)
            WRITE(6,*) ' Included spinorbitalexcitation ', JSPCOMP_ACT
            WRITE(6,'(A,16I4)')
     &      ' Creation of alpha     :',
     &      (ISPOBEX_TP(I+0*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)')
     &      ' Creation of beta      :',
     &      (ISPOBEX_TP(I+1*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)')
     &      ' Annihilation of alpha :',
     &      (ISPOBEX_TP(I+2*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)')
     &      ' Annihilation of beta  :',
     &      (ISPOBEX_TP(I+3*NGAS,JSPCOMP_ACT),I=1,NGAS)
          END DO
*
          WRITE(6,*) ' OX => SOX mapping '
          DO IOX = 1, NOBEX_TP
            N = NSOX_FOR_OX(IOX)
            IB = IBSOX_FOR_OX(IOX)
            WRITE(6,*)
            WRITE(6,*) ' ===================='
            WRITE(6,*) ' Info for OX : ', IOX
            WRITE(6,*) ' ===================='
            WRITE(6,*) ' Number of SOXS = ', N
            WRITE(6,*)
            WRITE(6,*) ' SOXS : '
            CALL IWRTMA(ISOX_FOR_OX(IB),1,N,1,N)
          END DO
*
        END IF
      END IF
* 
      CALL MEMREL('O_TO_S',WORK,KLSOEX,KLSOEX,KFREE,LFREE)
      RETURN
      END
*
      SUBROUTINE SPOBEX_FOR_OCCLS(
     &           IEXTP_TO_OCCLS,NOCCLS,ISOX_TO_OX,NSOX,
     &           NSOX_FOR_OCCLS,ISOX_FOR_OCCLS,IBSOX_FOR_OCCLS)
*
* Obtain spin-orbital excitation corresponding to given occupation class
*
* Jeppe Olsen, April 24 in Kerkrade The Netherlands
*
* Note : in input NSOX should be the extended number including
*        the unit operator
#include "implicit.inc"
*.  Input
      INTEGER IEXTP_TO_OCCLS(NOCCLS),ISOX_TO_OX(NSOX)
*. Output
      INTEGER NSOX_FOR_OCCLS(NOCCLS),IBSOX_FOR_OCCLS(NOCCLS)
      INTEGER ISOX_FOR_OCCLS(NSOX)
*
      IB = 1
      DO IOCCLS = 1, NOCCLS
C?      WRITE(6,*) ' IOCCLS ', IOCCLS
        IBSOX_FOR_OCCLS(IOCCLS) = IB
        N = 0
        DO ISOX = 1, NSOX
C?        WRITE(6,*) ' ISOX, ISOX_TO... ', ISOX,
C?   &                 ISOX_TO_OX(ISOX)
          IF(IEXTP_TO_OCCLS(ISOX_TO_OX(ISOX)).EQ.IOCCLS) THEN
             N = N + 1
             L =  IB - 1 + N
             ISOX_FOR_OCCLS(L) = ISOX
          END IF
        END DO
        IB = IB + N
        NSOX_FOR_OCCLS(IOCCLS) = N
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' spin-orbital excitations to occupations '
        DO IOCCLS = 1, NOCCLS
          WRITE(6,*) ' Occupation class : ', IOCCLS
          N = NSOX_FOR_OCCLS(IOCCLS)
          IB =  IBSOX_FOR_OCCLS(IOCCLS)
          WRITE(6,*) ' Number of spin orbital orbital types : ', N
          WRITE(6,*) ' The spin-orbital excitations'
          CALL IWRTMA(ISOX_FOR_OCCLS(IB),1,N,1,N)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE SPOBEXTP_PAIRS_KRCC(
     &           NSPOBEX_TP,ISPOBEX_TP,NGAS,ISPOBEX_PAIRS)
*
* Find pairs of spin-orbital excitaion types related to each other
* by spin-flip
*
* Jeppe Olsen, July11 2001
*
* Operators can in openshell cases not always be related by spin-flip
* since a certain configuration is choosen as fermi vacuum for the unrestricted case.
* They are therefore labelled with their own minus
* A possiblility to add operator have been added but currently not in use!
* Added by Lasse 2007
*
#include "implicit.inc"
#include "mxpdim.inc"
*. Input
C      DIMENSION ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
      DIMENSION ISPOBEX_TP(4*NGAS,2*NSPOBEX_TP)
*. Output
C      INTEGER ISPOBEX_PAIRS(NSPOBEX_TP)
      INTEGER ISPOBEX_PAIRS(2*NSPOBEX_TP)
*. Local scratch
      INTEGER IABFLIP(MXPNGAS*4)
*
      NSPOBEX_TPX = NSPOBEX_TP
      IZERO = 0
C      CALL ISETVC(ISPOBEX_PAIRS,IZERO,NSPOBEX_TP)
      CALL ISETVC(ISPOBEX_PAIRS,IZERO,2*NSPOBEX_TP)
*
      DO ITP = 1, NSPOBEX_TP
      IF(ISPOBEX_PAIRS(ITP).EQ.0) THEN
*. Perform spinflip on type ITP
        CALL ABFLIP_SPOXTP(ISPOBEX_TP(1,ITP),IABFLIP,NGAS)
*. Find address of spinflipped type in input list
        DO JTP = 1, NSPOBEX_TP
        IF(ISPOBEX_PAIRS(JTP).EQ.0) THEN
          IDIFF = 0
          DO IOP = 1, 4*NGAS
            IF(IABFLIP(IOP).NE.ISPOBEX_TP(IOP,JTP)) IDIFF = IDIFF + 1
          END DO
          IF(IDIFF.EQ.0) THEN
           ISPOBEX_PAIRS(ITP) = JTP
           ISPOBEX_PAIRS(JTP) = ITP
          END IF
        END IF
* Lasse added
        IF(JTP.EQ.NSPOBEX_TP.AND.ISPOBEX_PAIRS(ITP).EQ.0) THEN
         IADD = 0
         IF(IADD.EQ.0) THEN
          ISPOBEX_PAIRS(ITP) = - ITP
         ELSE
* Add new operator and increase NSPOBEX_TP
          NSPOBEX_TPX = NSPOBEX_TPX + 1
          CALL ICOPVE(IABFLIP,ISPOBEX_TP(1,NSPOBEX_TPX),4*NGAS)
          ISPOBEX_PAIRS(ITP) = NSPOBEX_TPX
          ISPOBEX_PAIRS(NSPOBEX_TPX) = ITP
         END IF
        END IF
* End Lasse addition
        END DO
*       ^ End of loop over JTP
      END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Table of paired spin-orbital excitations '
C        CALL IWRTMA(ISPOBEX_PAIRS,1,NSPOBEX_TP,1,NSPOBEX_TP)
        CALL IWRTMA(ISPOBEX_PAIRS,1,2*NSPOBEX_TP,1,2*NSPOBEX_TP)
      END IF
*
      RETURN
      END
*
      SUBROUTINE SPOBEX_TO_ABOBEX_CC(ISPOBEX_TP,NSPOBEX_TP,NGAS,
     &           IFLAG,NAOBEX_TP,NBOBEX_TP,IAOBEX_TP,IBOBEX_TP)
*
* Split spin-orbital excitations into alpha and beta-orbital excitations
*
* IFLAG = 1 : Find only number of alpha- and beta- orbital excitations
*
* Jeppe Olsen, July 2000
*
#include "implicit.inc"
*. Input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
*. Output
      INTEGER IAOBEX_TP(2*NGAS,*),IBOBEX_TP(2*NGAS,*)
*
C?    WRITE(6,*) ' SPOBEX_TO_ABOBEX_CC : NSPOBEX_TP = ', NSPOBEX_TP
      DO IAB = 1, 2
        LEX_TP = 0
        DO JSPOBEX_TP = 1, NSPOBEX_TP
*. Has this a or b excitation been observed before
          I_AM_OLD_HAT = 0
          DO KSPOBEX_TP = 1, JSPOBEX_TP-1
            I_DENTICAL = 1
            DO ICA = 1, 2
              IEXP_OFF = (IAB-1)*NGAS + (ICA-1)*2*NGAS + 1
              DO IGAS = 1, NGAS
                IF(ISPOBEX_TP(IGAS+IEXP_OFF-1,JSPOBEX_TP).NE.
     &             ISPOBEX_TP(IGAS+IEXP_OFF-1,KSPOBEX_TP)    )
     &             I_DENTICAL = 0
              END DO
            END DO
            IF(I_DENTICAL.EQ.1) I_AM_OLD_HAT = 1
          END DO
          IF(I_AM_OLD_HAT.EQ.0) THEN
            LEX_TP = LEX_TP + 1
            IF(IFLAG.EQ.0) THEN
              DO ICA = 1, 2
               IEXP_OFF =  (IAB-1)*NGAS + (ICA-1)*2*NGAS + 1
               IAB_OFF  =  (ICA-1)*NGAS + 1
               IF(IAB.EQ.1) THEN
                 CALL ICOPVE(ISPOBEX_TP(IEXP_OFF,JSPOBEX_TP),
     &                        IAOBEX_TP(IAB_OFF,LEX_TP),NGAS)
               ELSE IF (IAB.EQ.2) THEN
                 CALL ICOPVE(ISPOBEX_TP(IEXP_OFF,JSPOBEX_TP),
     &                        IBOBEX_TP(IAB_OFF,LEX_TP),NGAS)
               END IF
              END DO
            END IF
          END IF
*         ^ End of I_AM_OLD_HAT = 0
        END DO
*       ^ End of loop over JSPOBEX_TP
      IF(IAB.EQ.1) THEN
        NAOBEX_TP = LEX_TP
      ELSE
        NBOBEX_TP = LEX_TP
      END IF
      END DO
*     ^ End of loop over IAB
*
      NTEST = 00
      IF(NTEST.GE.3) THEN
       WRITE(6,*) ' Number of alpha-excitation operators ', NAOBEX_TP
       WRITE(6,*) ' Number of beta-excitation operators ',  NBOBEX_TP
      END IF
      IF(NTEST.GE.5.AND.IFLAG.EQ.0) THEN
        WRITE(6,*) ' Alpha-excitation operators : '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        DO JEX_TP = 1, NAOBEX_TP
         WRITE(6,*)
         WRITE(6,*) ' alphaorbitalexcitation ', JEX_TP
         WRITE(6,'(A,16I4)')
     &   ' Creation      :',  (IAOBEX_TP(I+0*NGAS,JEX_TP),I=1,NGAS)
         WRITE(6,'(A,16I4)')
     &   ' Annihilation  :',  (IAOBEX_TP(I+1*NGAS,JEX_TP),I=1,NGAS)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' beta-excitation operators : '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        DO JEX_TP = 1, NBOBEX_TP
         WRITE(6,*)
         WRITE(6,*) ' betaorbitalexcitation ', JEX_TP
         WRITE(6,'(A,16I4)')
     &   ' Creation      :',  (IBOBEX_TP(I+0*NGAS,JEX_TP),I=1,NGAS)
         WRITE(6,'(A,16I4)')
     &   ' Annihilation  :',  (IBOBEX_TP(I+1*NGAS,JEX_TP),I=1,NGAS)
        END DO
      END IF
*     ^ End of print is active
*
      RETURN
      END
*
      SUBROUTINE LEN_GENOP_STR_MAP_KRCC(NGENOP,IGENOP,
     &                             NSPGRP,ISPGRP,NOBPT,
     &                             NGAS,MAXLEN)
*
* A set of general operators, IGENOP, and a set of supergroups, ISPGRP,
* are defined. Obtain max.length of mappings from resolution strings
* to istrings.
*
* The resolution strings KSTR are assumed to be inserted between the
* creation and annihilation operators
*
*   <ISPGRP!ICREA!KSTR> <KSTR!IANNI!ISPGRP'>
* So the mappings from KSTR to ISTR are always creation mappings
*
* Jeppe Olsen, July 2000 ( At summerschool, HNIE)
*
#include "implicit.inc"
#include "mxpdim.inc"
*
*. Specific input
      INTEGER IGENOP(2*NGAS,*), ISPGRP(MXPNGAS,*), NOBPT(NGAS)
*. Local scratch
      INTEGER KSPGRP(MXPNGAS), IGENOP_EXP(MXPLCCOP)
*
      MAXLEN = 0
      IONE = 1
      IMONE = -1
      DO ICA = 1, 2
      DO JSPGRP = 1, NSPGRP
*. ICREA/IANNI(dag)*!KSTR> = !ISTR>
       DO JGENOP = 1, NGENOP
*. Occupation of KSPGRP
         IF(ICA.EQ.1) THEN
           CALL IVCSUM(KSPGRP,ISPGRP(1,JSPGRP),IGENOP(1,JGENOP),
     &                 IONE,IMONE,NGAS)
         ELSE
           CALL IVCSUM(KSPGRP,ISPGRP(1,JSPGRP),IGENOP(NGAS+1,JGENOP),
     &                 IONE,IMONE,NGAS)
         END IF
*. Is KSPGRP a correct supergroup( all occ larger than zero)
         I_AM_OKAY = 1
         DO IGAS = 1, NGAS
           IF(KSPGRP(IGAS).LT.0) I_AM_OKAY = 0
           IF(KSPGRP(IGAS).GT.NOBPT(IGAS)) I_AM_OKAY = 0
         END DO
         IF(I_AM_OKAY.EQ.1) THEN
*
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' Initial occ. of K supergroup '
             CALL IWRTMA(KSPGRP,1,NGAS,1,NGAS)
             WRITE(6,*) ' Active part of IGENOP'
             IF(ICA.EQ.1) THEN
               CALL IWRTMA(IGENOP(1,JGENOP),1,NGAS,1,NGAS)
             ELSE
               CALL IWRTMA(IGENOP(NGAS+1,JGENOP),1,NGAS,1,NGAS)
             END IF
           END IF
*
*. IGENOP in expanded form
C  REF_OP(IOPGAS,IOP,NOP,NGAS,IWAY)
          LEN = 0
          IF(ICA.EQ.1) THEN
            CALL REF_OP_KRCC(IGENOP(1,JGENOP),IGENOP_EXP,NOP,
     &                  NGAS,1,MXPLCCOP)
          ELSE
            CALL REF_OP_KRCC(IGENOP(NGAS+1,JGENOP),IGENOP_EXP,NOP,
     &                  NGAS,1,MXPLCCOP)
          END IF
          DO IOP = 1, NOP
            IOPTP = IGENOP_EXP(IOP)
*. Number of Kstrings
C                  NST_FOR_OCC(NELEC_PER_GAS,NORB_PER_GAS,NGAS)
            NSTR = NST_FOR_OCC(KSPGRP,NOBPT,NGAS)
            NORB = NOBPT(IOPTP)
            LEN = LEN + NSTR*NORB
*. Update KSPGRP
            KSPGRP(IOPTP) = KSPGRP(IOPTP) + 1
          END DO
          MAXLEN = MAX(MAXLEN,LEN)
         END IF
       END DO
      END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Max length of KSTR => ISTR map = ', MAXLEN
      END IF
*
      RETURN
      END
*
      SUBROUTINE REF_OP_KRCC(IOPGAS,IOP,NOP,NGAS,IWAY,MXNOP)
* An operatorstring may be specifed as
*
* IOPGAS : Number of operators per GASspace
* IOP    : GASpace of each operator
*
* Transform between these two form
*
* Iway = 1 : IOPGAS => IOP
* Iway = 2 : IOP    => IOPGAS
*
* Jeppe Olsen, Summer of 99
*
#include "implicit.inc"
*. Input/Output
      INTEGER IOPGAS(NGAS),IOP(MXNOP)
*
      IF(IWAY.EQ.1) THEN
        JOP = 0
        DO JGAS = 1, NGAS
          LJGAS = IOPGAS(JGAS)
          DO JJOP = 1, LJGAS
            JOP = JOP + 1
            IOP(JOP) = JGAS
          END DO
        END DO
        NOP = JOP
        if (NOP.gt.MXNOP) then
          write(6,*) 'NOP too large in REF_OP. NOP =',NOP
          write(6,*) 'Increase MXNOP or find bug !!'
        end if
      ELSE
        DO JGAS = 0, NGAS
          JOP = 0
          DO JJOP = 1, NOP
            IF(IOP(JJOP).EQ.JGAS) JOP = JOP +1
          END DO
          IOPGAS(JGAS) = JOP
        END DO
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.1) THEN
         WRITE(6,*) ' IOPGAS => IOP '
        ELSE
         WRITE(6,*) ' IOP => IOPGAS '
        END IF
        WRITE(6,*) ' IOPGAS and IOP '
        CALL IWRTMA(IOPGAS,1,NGAS,1,NGAS)
        CALL IWRTMA(IOP,1,NOP,1,NOP)
      END IF
*
      RETURN
      END
*
      SUBROUTINE STUPID_ROUTINE(ITSOSO_TP,NTSOSO_TP,NGAS)
#include "implicit.inc"
      INTEGER ITSOSO_TP(4*NGAS,NTSOSO_TP)
      DO J=1,4*NGAS
      DO I = 1,NTSOSO_TP
        WRITE(6,*) 'ITSOSO_TP',ITSOSO_TP(J,I)
      END DO
      END DO
      END
*
      SUBROUTINE IDIM_TCC_KRCC(ITSOSO_TP,NTSOSO_TP,ISYM,
     &           MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &           LTSOSO_TP,IBTSOSO_TP,IDIM_T,
     &           MX_TBLK_AS,IFTONE,IFTTWO,ILTTWO,N1ELINT,N2ELINT)
*
* Dimension of T operators in spin-orbital basis
*
* Largest number of strings of given sym in T(ICA,ICB,IAA,IAB)
* i.e. largest block of ICA, ICB,IAA,IAB of given sym
*
* Size of block required to hold above blocks
*
* Jeppe Olsen, Summer of 99
*
* Double group version
*
#include "implicit.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "csm.inc"
#include "symm.inc"
*. Specific input
      INTEGER ITSOSO_TP(4*NGAS,NTSOSO_TP)
*. Local scratch
      INTEGER ICA_GRP(MXPNGAS), ICB_GRP(MXPNGAS)
      INTEGER IAA_GRP(MXPNGAS), IAB_GRP(MXPNGAS)
*.Output : Length of each type and offset for each type
      INTEGER LTSOSO_TP(NTSOSO_TP),IBTSOSO_TP(NTSOSO_TP)
*
      NTEST = 0003
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' IDIM_TCC : ISYM = ', ISYM
      END IF
*
      LENGTH = 0
      LENGTH2 = 0
      MX_ST_TSOSO = 0
      MX_ST_TSOSO_BLK = 0
      MX_TBLK = 0
      MX_SBSTR = 0
      MX_TBLK_AS = 0
*
      DO ITSS_TP = 1, NTSOSO_TP
       IF(NTEST.GE.100) WRITE(6,*) ' ITSS_TP = ', ITSS_TP
       IF(NTEST.GE.100) CALL WRT_SPOX_TP_CC_KRCC(ITSOSO_TP(1,ITSS_TP),1)
*. Occupation to group translation
       CALL OCC_TO_GRP_KRCC(ITSOSO_TP(1+0*NGAS,ITSS_TP),ICA_GRP,1)
       CALL OCC_TO_GRP_KRCC(ITSOSO_TP(1+1*NGAS,ITSS_TP),ICB_GRP,1)
       CALL OCC_TO_GRP_KRCC(ITSOSO_TP(1+2*NGAS,ITSS_TP),IAA_GRP,1)
       CALL OCC_TO_GRP_KRCC(ITSOSO_TP(1+3*NGAS,ITSS_TP),IAB_GRP,1)
*
       if (NTEST.ge.50) then
         write(6,*) ' Output from OCC_TO_GRP : '
         write(6,*) ' ICA_GRP : ',(ICA_GRP(I),I=1,NGAS,1)
         write(6,*) ' ICA_GRP : ',(ICB_GRP(I),I=1,NGAS,1)
         write(6,*) ' ICA_GRP : ',(IAA_GRP(I),I=1,NGAS,1)
         write(6,*) ' ICA_GRP : ',(IAB_GRP(I),I=1,NGAS,1)
       end if
*
       NEL_CA = IELSUM(ITSOSO_TP(1+0*NGAS,ITSS_TP),NGAS)
       NEL_CB = IELSUM(ITSOSO_TP(1+1*NGAS,ITSS_TP),NGAS)
       NEL_AA = IELSUM(ITSOSO_TP(1+2*NGAS,ITSS_TP),NGAS)
       NEL_AB = IELSUM(ITSOSO_TP(1+3*NGAS,ITSS_TP),NGAS)
*
       IBTSOSO_TP(ITSS_TP) = LENGTH + 1
       DO I_CR_SM = 1, NSMST
*. symmetry of annihilation strings
        I_AN_SM = IDBGMULT(ISYM,INVELM(I_CR_SM))
*. symmetry of creation strings corresponding to annihilation strings ...
        I_AN_SM = IADJSYM(I_AN_SM)
        DO I_CR_AL_SM = 1, NSMST
          DO II_AN_AL_SM = 1, NSMST
             I_AN_AL_SM = IADJSYM(II_AN_AL_SM)
*
             I_CR_BE_SM = IDBGMULT(INVELM(I_CR_AL_SM),I_CR_SM)
             I_AN_BE_SM = IDBGMULT(INVELM(I_AN_AL_SM),I_AN_SM)
*
             IUB = 1
             CALL NST_SPGRP_KRCC(NGAS,ICA_GRP,I_CR_AL_SM,NSMST,
     &                          LEN_CA,NDIST_CA,IUB)
             IUB = 2
             CALL NST_SPGRP_KRCC(NGAS,ICB_GRP,I_CR_BE_SM,NSMST,
     &                          LEN_CB,NDIST_CB,IUB)
             IUB = 1
             CALL NST_SPGRP_KRCC(NGAS,IAA_GRP,I_AN_AL_SM,NSMST,
     &                          LEN_AA,NDIST_AA,IUB)
             IUB = 2
             CALL NST_SPGRP_KRCC(NGAS,IAB_GRP,I_AN_BE_SM,NSMST,
     &                          LEN_AB,NDIST_AB,IUB)
*
             IF(NTEST.GE.1000) THEN
               WRITE(6,'(A,4I5)') ' sym of CA CB AA AB',
     &         I_CR_AL_SM, I_CR_BE_SM, I_AN_AL_SM, I_AN_BE_SM
               WRITE(6,'(A,4I5)') ' LEN_CA, LEN_CB, LEN_AA, LEN_AB ',
     &                      LEN_CA, LEN_CB, LEN_AA, LEN_AB
             END IF
*
             LENGTH = LENGTH + LEN_CA*LEN_CB*LEN_AA*LEN_AB
             LENGTH2 = LENGTH2 + LEN_CA*LEN_CB*LEN_AA*LEN_AB
             if (NTEST.ge.1000) then
               write(6,*) 'LENGTH ',LEN_CA*LEN_CB*LEN_AA*LEN_AB
               write(6,*) 'LENGTH2 ',LENGTH2
             end if
*
             MX_ST_TSOSO =
     &       MAX(MX_ST_TSOSO,LEN_CA,LEN_CB,LEN_AA,LEN_AB)
             MX_ST_TSOSO_BLK =
     &       MAX(MX_ST_TSOSO_BLK,LEN_CA*NEL_CA,LEN_CB*NEL_CB,
     &                           LEN_AA*NEL_CA,LEN_AB*NEL_CB)
             MX_TBLK = MAX(MX_TBLK,LEN_CA*LEN_CB*LEN_AA*LEN_AB)
             MX_SBSTR = MAX(MX_SBSTR,LEN_CA,LEN_CB,LEN_AA,LEN_AB)
C            print*,'MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,MX_SBSTR',
C    &               MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,MX_SBSTR
          END DO
        END DO
       END DO
C     IF(ILTTWO.NE.-1) THEN
      LTSOSO_TP(ITSS_TP) = LENGTH2
      MX_TBLK_AS = MAX(MX_TBLK_AS,LENGTH2)
C     END IF
      END DO
*
      IDIM_T = LENGTH
*
* Number of non-spin-flip one- and two electron integrals
* (for diagonal):
* If statement added as temp fix by Lasse since IDIM_TCC_DBG is called
* in SIGDEN_CC_REL where there is no need for N1ELINT and N2ELINT.
* Used ILTTWO as a switch since it is not in use.
      IF(ILTTWO.NE.-1) THEN
      N1ELINT = IBTSOSO_TP(IFTTWO) - IBTSOSO_TP(IFTONE)
      N2ELINT = LENGTH - N1ELINT
      END IF
*
      IF(NTEST.GE.3) THEN
        WRITE(6,*) ' Number of blocks ',NTSOSO_TP
        WRITE(6,*) ' Number of T-coefficients ', LENGTH
        WRITE(6,*) ' Largest symmetry block of T ', MX_TBLK
        WRITE(6,*) ' Largest block of T ', MX_TBLK_AS
        WRITE(6,*) ' Largest substring : ', MX_SBSTR
        IF(ILTTWO.NE.-1) THEN
          WRITE(6,*) ' Number of 1e intergrals ', N1ELINT
          WRITE(6,*) ' Number of 2e intergrals ', N2ELINT
        END IF
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
     &  ' Largest number of strings of given sym and type ',
     &               MX_ST_TSOSO
        WRITE(6,*)
     & ' Largest block of strings of given sym and type ',
     &               MX_ST_TSOSO_BLK
        WRITE(6,*) ' Number of elements per block '
        CALL IWRTMA(LTSOSO_TP,1,NTSOSO_TP,1,NTSOSO_TP)
        WRITE(6,*) ' Offset for each block '
        CALL IWRTMA(IBTSOSO_TP,1,NTSOSO_TP,1,NTSOSO_TP)
      END IF
*
      RETURN
      END
*
      SUBROUTINE SET_INTERMEDIATES_MASTER(WORK,KFREE,LFREE)
*
* Driver routine for setting up the intermediates 
* Intermediate operators will be setup according to (P,H)
* There are 5 possible intermediates (when assuming P>H)
* namely (1,2),(1,1),(1,0),(0,2),(0,1)
* Can be extended if so desire...
*
* Lasse 2010
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "crun.inc"
#include "lucinp.inc"
#include "ctcc.inc"
#include "ctccp.inc"
* for real or complex groups
#include "../include/dgroup.h"
*
      DIMENSION WORK(*)
*
      print*,'say hi'
*
* Will first havve to ensure that all possible ABOBEX are accounted for.
* This may not be the case for real and complex groups
* Will do this seperately for alpha and beta types
      IF(NZ.LE.2) THEN
* First alpha
        INI = 0
        CALL FIND_ALL_ABOBEX(INI,NGAS,WORK(KLAOBEX_CC),
     &                       NAOBEX_TP,IDUMMY,LENA_KRCC)
*
        CALL MEMGET('INTE',KLAOBEX_KRCC,2*NGAS*LENA_KRCC,
     &              WORK,KFREE,LFREE)
        INI = 1
        CALL FIND_ALL_ABOBEX(INI,NGAS,WORK(KLAOBEX_CC),
     &                       NAOBEX_TP,WORK(KLAOBEX_KRCC),LENA_KRCC)
        CALL MEMCHK_KRCC(WORK)
* Now beta
        INI = 0
        CALL FIND_ALL_ABOBEX(INI,NGAS,WORK(KLBOBEX_CC),
     &                       NBOBEX_TP,IDUMMY,LENB_KRCC)
*
        CALL MEMGET('INTE',KLBOBEX_KRCC,2*NGAS*LENB_KRCC,
     &              WORK,KFREE,LFREE)
        INI = 1
        CALL FIND_ALL_ABOBEX(INI,NGAS,WORK(KLBOBEX_CC),
     &                       NBOBEX_TP,WORK(KLBOBEX_KRCC),LENB_KRCC)
        CALL MEMCHK_KRCC(WORK)
      END IF 
*
* Values
*
      IZERO = 0
*
* Setup (1,0)
*
      IP = 1
      IH = 0
      IKRFLIPMAX = MIN(2*MK2DEL,MX_EXC_LEVEL)
      IKRFLIPMAX = MIN(2+MK2DEL,IKRFLIPMAX)
      IHOLEATOT = NGAS+1
      IHOLEBTOT = NGAS+1
      IPARTATOT = NGAS+1
      IPARTBTOT = NGAS+1
*
* Allocate a few dummy arrays
*
      CALL MEMGET('INTE',KIHOLEA,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIHOLEB,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEA,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEB,NGAS**2,WORK,KFREE,LFREE)
*
* Set vectors
*
      CALL ISETVC(WORK(KIHOLEA),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIHOLEB),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIPARTICLEA),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIPARTICLEB),IZERO,NGAS**2)
*
* Initialize
*
      IINI = 0
*
      CALL SET_INTMEDIATES(IINI,IP,IH,NINTER10,IDUMMY,IKRFLIPMAX,
     &                     IHOLEATOT,IHOLEBTOT,
     &                     WORK(KIHOLEA),WORK(KIHOLEB),
     &                     IPARTATOT,IPARTBTOT,
     &                     WORK(KIPARTICLEA),WORK(KIPARTICLEB),
     &                     LENA_KRCC,LENB_KRCC,
     &                     WORK(KLAOBEX_KRCC),WORK(KLBOBEX_KRCC))
*
* Remove dummy arrays
*
      CALL MEMCHK_KRCC(WORK)
      CALL MEMREL('INT10',WORK,KIHOLEA,KIHOLEA,KFREE,LFREE)
*
* Allocate permanent array
*
      CALL MEMGET('INTE',KINTM10,4*NGAS*NINTER10,WORK,KFREE,LFREE)
*
* Allocate a few dummy arrays
*
      CALL MEMGET('INTE',KIHOLEA,IHOLEATOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIHOLEB,IHOLEBTOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEA,IPARTATOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEB,IPARTBTOT*NGAS,WORK,KFREE,LFREE)
*
* Set vectors
*
      CALL ISETVC(WORK(KIHOLEA),IZERO,IHOLEATOT*NGAS)
      CALL ISETVC(WORK(KIHOLEB),IZERO,IHOLEBTOT*NGAS)
      CALL ISETVC(WORK(KIPARTICLEA),IZERO,IPARTATOT*NGAS)
      CALL ISETVC(WORK(KIPARTICLEB),IZERO,IPARTBTOT*NGAS)
*
* Find the 01 intermediates
*
      CALL MEMCHK_KRCC(WORK)
      IINI = 1
*
      CALL SET_INTMEDIATES(IINI,IP,IH,NINTER10,WORK(KINTM10),IKRFLIPMAX,
     &                     IHOLEATOT,IHOLEBTOT,
     &                     WORK(KIHOLEA),WORK(KIHOLEB),
     &                     IPARTATOT,IPARTBTOT,
     &                     WORK(KIPARTICLEA),WORK(KIPARTICLEB),
     &                     LENA_KRCC,LENB_KRCC,
     &                     WORK(KLAOBEX_KRCC),WORK(KLBOBEX_KRCC))
*
* Remove dummy arrays
*
      CALL MEMCHK_KRCC(WORK)
      CALL MEMREL('INT10',WORK,KIHOLEA,KIHOLEA,KFREE,LFREE)
*
* Setup (0,1)
*
      IP = 0
      IH = 1
      IKRFLIPMAX = MIN(2*MK2DEL,MX_EXC_LEVEL)
      IHOLEATOT = NGAS+1
      IHOLEBTOT = NGAS+1
      IPARTATOT = NGAS+1
      IPARTBTOT = NGAS+1
*
* Allocate a few dummy arrays
*
      CALL MEMGET('INTE',KIHOLEA,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIHOLEB,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEA,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEB,NGAS**2,WORK,KFREE,LFREE)
*
* Set vectors
*
      CALL ISETVC(WORK(KIHOLEA),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIHOLEB),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIPARTICLEA),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIPARTICLEB),IZERO,NGAS**2)
*
* Initialize
*
      IINI = 0
*
      CALL SET_INTMEDIATES(IINI,IP,IH,NINTER01,IDUMMY,IKRFLIPMAX,
     &                     IHOLEATOT,IHOLEBTOT,
     &                     WORK(KIHOLEA),WORK(KIHOLEB),
     &                     IPARTATOT,IPARTBTOT,
     &                     WORK(KIPARTICLEA),WORK(KIPARTICLEB),
     &                     LENA_KRCC,LENB_KRCC,
     &                     WORK(KLAOBEX_KRCC),WORK(KLBOBEX_KRCC))
*
* Remove dummy arrays
*
      CALL MEMCHK_KRCC(WORK)
      CALL MEMREL('INT01',WORK,KIHOLEA,KIHOLEA,KFREE,LFREE)
*
* Allocate permanent array
*
      CALL MEMGET('INTE',KINTM01,4*NGAS*NINTER01,WORK,KFREE,LFREE)
*
* Allocate a few dummy arrays
*
      CALL MEMGET('INTE',KIHOLEA,IHOLEATOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIHOLEB,IHOLEBTOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEA,IPARTATOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEB,IPARTBTOT*NGAS,WORK,KFREE,LFREE)
*
* Set vectors
*
      CALL ISETVC(WORK(KIHOLEA),IZERO,IHOLEATOT*NGAS)
      CALL ISETVC(WORK(KIHOLEB),IZERO,IHOLEBTOT*NGAS)
      CALL ISETVC(WORK(KIPARTICLEA),IZERO,IPARTATOT*NGAS)
      CALL ISETVC(WORK(KIPARTICLEB),IZERO,IPARTBTOT*NGAS)
*
* Find the 01 intermediates
*
      CALL MEMCHK_KRCC(WORK)
      IINI = 1
*
      CALL SET_INTMEDIATES(IINI,IP,IH,NINTER01,WORK(KINTM01),IKRFLIPMAX,
     &                     IHOLEATOT,IHOLEBTOT,
     &                     WORK(KIHOLEA),WORK(KIHOLEB),
     &                     IPARTATOT,IPARTBTOT,
     &                     WORK(KIPARTICLEA),WORK(KIPARTICLEB),
     &                     LENA_KRCC,LENB_KRCC,
     &                     WORK(KLAOBEX_KRCC),WORK(KLBOBEX_KRCC))
*
* Remove dummy arrays
*
      CALL MEMCHK_KRCC(WORK)
      CALL MEMREL('INT01',WORK,KIHOLEA,KIHOLEA,KFREE,LFREE)
*
* Setup (1,1)
*
      IP = 1
      IH = 1
      IKRFLIPMAX = MIN(2*MK2DEL,MX_EXC_LEVEL)
      IKRFLIPMAX = MIN(2+1*MK2DEL,IKRFLIPMAX)
      IHOLEATOT = NGAS+1
      IHOLEBTOT = NGAS+1
      IPARTATOT = NGAS+1
      IPARTBTOT = NGAS+1
*
* Allocate a few dummy arrays
*
      CALL MEMGET('INTE',KIHOLEA,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIHOLEB,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEA,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEB,NGAS**2,WORK,KFREE,LFREE)
*
* Set vectors
*
      CALL ISETVC(WORK(KIHOLEA),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIHOLEB),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIPARTICLEA),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIPARTICLEB),IZERO,NGAS**2)
*
* Initialize
*
      IINI = 0
*
      CALL SET_INTMEDIATES(IINI,IP,IH,NINTER11,IDUMMY,IKRFLIPMAX,
     &                     IHOLEATOT,IHOLEBTOT,
     &                     WORK(KIHOLEA),WORK(KIHOLEB),
     &                     IPARTATOT,IPARTBTOT,
     &                     WORK(KIPARTICLEA),WORK(KIPARTICLEB),
     &                     LENA_KRCC,LENB_KRCC,
     &                     WORK(KLAOBEX_KRCC),WORK(KLBOBEX_KRCC))
*
* Remove dummy arrays
*
      CALL MEMCHK_KRCC(WORK)
      CALL MEMREL('INT11',WORK,KIHOLEA,KIHOLEA,KFREE,LFREE)
*
* Allocate permanent array
*
      CALL MEMGET('INTE',KINTM11,4*NGAS*NINTER11,WORK,KFREE,LFREE)
*
* Allocate a few dummy arrays
*
      CALL MEMGET('INTE',KIHOLEA,IHOLEATOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIHOLEB,IHOLEBTOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEA,IPARTATOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEB,IPARTBTOT*NGAS,WORK,KFREE,LFREE)
*
* Set vectors
*
      CALL ISETVC(WORK(KIHOLEA),IZERO,IHOLEATOT*NGAS)
      CALL ISETVC(WORK(KIHOLEB),IZERO,IHOLEBTOT*NGAS)
      CALL ISETVC(WORK(KIPARTICLEA),IZERO,IPARTATOT*NGAS)
      CALL ISETVC(WORK(KIPARTICLEB),IZERO,IPARTBTOT*NGAS)
*
* Find the 11 intermediates
*
      CALL MEMCHK_KRCC(WORK)
      IINI = 1
*
      CALL SET_INTMEDIATES(IINI,IP,IH,NINTER11,WORK(KINTM11),IKRFLIPMAX,
     &                     IHOLEATOT,IHOLEBTOT,
     &                     WORK(KIHOLEA),WORK(KIHOLEB),
     &                     IPARTATOT,IPARTBTOT,
     &                     WORK(KIPARTICLEA),WORK(KIPARTICLEB),
     &                     LENA_KRCC,LENB_KRCC,
     &                     WORK(KLAOBEX_KRCC),WORK(KLBOBEX_KRCC))
*
* Remove dummy arrays
*
      CALL MEMCHK_KRCC(WORK)
      CALL MEMREL('INT11',WORK,KIHOLEA,KIHOLEA,KFREE,LFREE)
*
* Setup (2,0)
*
      IP = 2
      IH = 0
      IKRFLIPMAX = MIN(3*MK2DEL,MX_EXC_LEVEL)
      IKRFLIPMAX = MIN(2+2*MK2DEL,IKRFLIPMAX)
      IHOLEATOT = NGAS+1
      IHOLEBTOT = NGAS+1
      IPARTATOT = NGAS+1
      IPARTBTOT = NGAS+1
*
* Allocate a few dummy arrays
*
      CALL MEMGET('INTE',KIHOLEA,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIHOLEB,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEA,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEB,NGAS**2,WORK,KFREE,LFREE)
*
* Set vectors
*
      CALL ISETVC(WORK(KIHOLEA),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIHOLEB),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIPARTICLEA),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIPARTICLEB),IZERO,NGAS**2)
*
* Initialize
*
      IINI = 0
*
      CALL SET_INTMEDIATES(IINI,IP,IH,NINTER20,IDUMMY,IKRFLIPMAX,
     &                     IHOLEATOT,IHOLEBTOT,
     &                     WORK(KIHOLEA),WORK(KIHOLEB),
     &                     IPARTATOT,IPARTBTOT,
     &                     WORK(KIPARTICLEA),WORK(KIPARTICLEB),
     &                     LENA_KRCC,LENB_KRCC,
     &                     WORK(KLAOBEX_KRCC),WORK(KLBOBEX_KRCC))
*
* Remove dummy arrays
*
      CALL MEMCHK_KRCC(WORK)
      CALL MEMREL('INT20',WORK,KIHOLEA,KIHOLEA,KFREE,LFREE)
*
* Allocate permanent array
*
      CALL MEMGET('INTE',KINTM20,4*NGAS*NINTER20,WORK,KFREE,LFREE)
*
* Allocate a few dummy arrays
*
      CALL MEMGET('INTE',KIHOLEA,IHOLEATOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIHOLEB,IHOLEBTOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEA,IPARTATOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEB,IPARTBTOT*NGAS,WORK,KFREE,LFREE)
*
* Set vectors
*
      CALL ISETVC(WORK(KIHOLEA),IZERO,IHOLEATOT*NGAS)
      CALL ISETVC(WORK(KIHOLEB),IZERO,IHOLEBTOT*NGAS)
      CALL ISETVC(WORK(KIPARTICLEA),IZERO,IPARTATOT*NGAS)
      CALL ISETVC(WORK(KIPARTICLEB),IZERO,IPARTBTOT*NGAS)
*
* Find the 20 intermediates
*
      CALL MEMCHK_KRCC(WORK)
      IINI = 1
*
      CALL SET_INTMEDIATES(IINI,IP,IH,NINTER20,WORK(KINTM20),IKRFLIPMAX,
     &                     IHOLEATOT,IHOLEBTOT,
     &                     WORK(KIHOLEA),WORK(KIHOLEB),
     &                     IPARTATOT,IPARTBTOT,
     &                     WORK(KIPARTICLEA),WORK(KIPARTICLEB),
     &                     LENA_KRCC,LENB_KRCC,
     &                     WORK(KLAOBEX_KRCC),WORK(KLBOBEX_KRCC))
*
* Remove dummy arrays
*
      CALL MEMCHK_KRCC(WORK)
      CALL MEMREL('INT20',WORK,KIHOLEA,KIHOLEA,KFREE,LFREE)
*
* Setup (0,2)
*
      IP = 0
      IH = 2
      IKRFLIPMAX = MIN(3*MK2DEL,MX_EXC_LEVEL)
      IKRFLIPMAX = MIN(2+2*MK2DEL,IKRFLIPMAX)
      IHOLEATOT = NGAS+1
      IHOLEBTOT = NGAS+1
      IPARTATOT = NGAS+1
      IPARTBTOT = NGAS+1
*
* Allocate a few dummy arrays
*
      CALL MEMGET('INTE',KIHOLEA,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIHOLEB,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEA,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEB,NGAS**2,WORK,KFREE,LFREE)
*
* Set vectors
*
      CALL ISETVC(WORK(KIHOLEA),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIHOLEB),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIPARTICLEA),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIPARTICLEB),IZERO,NGAS**2)
*
* Initialize
*
      IINI = 0
*
      CALL SET_INTMEDIATES(IINI,IP,IH,NINTER02,IDUMMY,IKRFLIPMAX,
     &                     IHOLEATOT,IHOLEBTOT,
     &                     WORK(KIHOLEA),WORK(KIHOLEB),
     &                     IPARTATOT,IPARTBTOT,
     &                     WORK(KIPARTICLEA),WORK(KIPARTICLEB),
     &                     LENA_KRCC,LENB_KRCC,
     &                     WORK(KLAOBEX_KRCC),WORK(KLBOBEX_KRCC))
*
* Remove dummy arrays
*
      CALL MEMCHK_KRCC(WORK)
      CALL MEMREL('INT02',WORK,KIHOLEA,KIHOLEA,KFREE,LFREE)
*
* Allocate permanent array
*
      CALL MEMGET('INTE',KINTM02,4*NGAS*NINTER02,WORK,KFREE,LFREE)
*
* Allocate a few dummy arrays
*
      CALL MEMGET('INTE',KIHOLEA,IHOLEATOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIHOLEB,IHOLEBTOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEA,IPARTATOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEB,IPARTBTOT*NGAS,WORK,KFREE,LFREE)
*
* Set vectors
*
      CALL ISETVC(WORK(KIHOLEA),IZERO,IHOLEATOT*NGAS)
      CALL ISETVC(WORK(KIHOLEB),IZERO,IHOLEBTOT*NGAS)
      CALL ISETVC(WORK(KIPARTICLEA),IZERO,IPARTATOT*NGAS)
      CALL ISETVC(WORK(KIPARTICLEB),IZERO,IPARTBTOT*NGAS)
*
* Find the 02 intermediates
*
      CALL MEMCHK_KRCC(WORK)
      IINI = 1
*
      CALL SET_INTMEDIATES(IINI,IP,IH,NINTER02,WORK(KINTM02),IKRFLIPMAX,
     &                     IHOLEATOT,IHOLEBTOT,
     &                     WORK(KIHOLEA),WORK(KIHOLEB),
     &                     IPARTATOT,IPARTBTOT,
     &                     WORK(KIPARTICLEA),WORK(KIPARTICLEB),
     &                     LENA_KRCC,LENB_KRCC,
     &                     WORK(KLAOBEX_KRCC),WORK(KLBOBEX_KRCC))
*
* Remove dummy arrays
*
      CALL MEMCHK_KRCC(WORK)
      CALL MEMREL('INT02',WORK,KIHOLEA,KIHOLEA,KFREE,LFREE)
*
* Setup (1,2)
*
      IP = 1
      IH = 2
      IKRFLIPMAX = MIN(2*MK2DEL,MX_EXC_LEVEL)
      IKRFLIPMAX = MIN(2+1*MK2DEL,IKRFLIPMAX)
      IHOLEATOT = NGAS+1
      IHOLEBTOT = NGAS+1
      IPARTATOT = NGAS+1
      IPARTBTOT = NGAS+1
*
* Allocate a few dummy arrays
*
      CALL MEMGET('INTE',KIHOLEA,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIHOLEB,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEA,NGAS**2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEB,NGAS**2,WORK,KFREE,LFREE)
*
* Set vectors
*
      CALL ISETVC(WORK(KIHOLEA),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIHOLEB),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIPARTICLEA),IZERO,NGAS**2)
      CALL ISETVC(WORK(KIPARTICLEB),IZERO,NGAS**2)
*
* Initialize
*
      IINI = 0
*
      CALL SET_INTMEDIATES(IINI,IP,IH,NINTER12,IDUMMY,IKRFLIPMAX,
     &                     IHOLEATOT,IHOLEBTOT,
     &                     WORK(KIHOLEA),WORK(KIHOLEB),
     &                     IPARTATOT,IPARTBTOT,
     &                     WORK(KIPARTICLEA),WORK(KIPARTICLEB),
     &                     LENA_KRCC,LENB_KRCC,
     &                     WORK(KLAOBEX_KRCC),WORK(KLBOBEX_KRCC))
*    
* Remove dummy arrays      
*    
      CALL MEMCHK_KRCC(WORK)
      CALL MEMREL('INT12',WORK,KIHOLEA,KIHOLEA,KFREE,LFREE)
*    
* Allocate permanent array
*
      CALL MEMGET('INTE',KINTM12,4*NGAS*NINTER12,WORK,KFREE,LFREE)
*     
* Allocate a few dummy arrays
*
      CALL MEMGET('INTE',KIHOLEA,IHOLEATOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIHOLEB,IHOLEBTOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEA,IPARTATOT*NGAS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPARTICLEB,IPARTBTOT*NGAS,WORK,KFREE,LFREE)
*
* Set vectors
*     
      CALL ISETVC(WORK(KIHOLEA),IZERO,IHOLEATOT*NGAS)
      CALL ISETVC(WORK(KIHOLEB),IZERO,IHOLEBTOT*NGAS)
      CALL ISETVC(WORK(KIPARTICLEA),IZERO,IPARTATOT*NGAS)
      CALL ISETVC(WORK(KIPARTICLEB),IZERO,IPARTBTOT*NGAS)
*
* Find the 12 intermediates
*     
      CALL MEMCHK_KRCC(WORK)
      IINI = 1
*     
      CALL SET_INTMEDIATES(IINI,IP,IH,NINTER12,WORK(KINTM12),IKRFLIPMAX,
     &                     IHOLEATOT,IHOLEBTOT,
     &                     WORK(KIHOLEA),WORK(KIHOLEB),
     &                     IPARTATOT,IPARTBTOT,
     &                     WORK(KIPARTICLEA),WORK(KIPARTICLEB),
     &                     LENA_KRCC,LENB_KRCC,
     &                     WORK(KLAOBEX_KRCC),WORK(KLBOBEX_KRCC))
*    
* Remove dummy arrays
*
      CALL MEMCHK_KRCC(WORK)
      CALL MEMREL('INT12',WORK,KIHOLEA,KIHOLEA,KFREE,LFREE)
*
      RETURN
      END
*
      SUBROUTINE SET_INTMEDIATES(IINI,IP,IH,NINTER,INTER,IKRFLIPMAX,
     &                           IHOLEATOT,IHOLEBTOT,
     &                           IHOLEA,IHOLEB,
     &                           IPARTATOT,IPARTBTOT,
     &                           IPARTICLEA,IPARTICLEB,
     &                           NAOBEX_TP,NBOBEX_TP,
     &                           IAOBEX_TP,IBOBEX_TP)
*
* Will find the actual intermediates depending on IH and IP
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
* Output
      INTEGER INTER(4*NGAS,*)
      INTEGER IHOLEA(NGAS,*),IHOLEB(NGAS,*)
      INTEGER IPARTICLEA(NGAS,*),IPARTICLEB(NGAS,*)
* Input
      INTEGER IAOBEX_TP(2*NGAS,NAOBEX_TP),IBOBEX_TP(2*NGAS,NBOBEX_TP)
*
      NTEST = 00
*
* Reset things
*
      IHOLEATOT = 0
      IHOLEBTOT = 0
      IPARTATOT = 0
      IPARTBTOT = 0
      IHOLEAIDENT = 0
      IHOLEBIDENT = 0
      IPARTAIDENT = 0
      IPARTBIDENT = 0
      NINTER = 0
      IALLOWIDENT = 0
C     print*,'start'
C     print*,IHOLEA(1,1),IHOLEA(2,1)
C     print*,IHOLEB(1,1),IHOLEB(2,1)
*
* Now distribute IH number of holes in IHOLEATOT + IHOLEBTOT number of GAS
      IF(IH.GE.1) THEN
        DO IGAS = 1,NGAS
          IF(IHPVGAS_AB(IGAS,1).EQ.1) THEN
* place an alpha hole in IGAS
            IHOLEATOT = IHOLEATOT + 1
            IHOLEA(IGAS,IHOLEATOT) = 1
          END IF
          IF(IHPVGAS_AB(IGAS,2).EQ.1) THEN
* place an beta hole in IGAS
            IHOLEBTOT = IHOLEBTOT + 1
            IHOLEB(IGAS,IHOLEBTOT) = 1
          END IF
* Add identity operator for holes in the valence space
* These are taken from T 
          IF(IHPVGAS_AB(IGAS,1).EQ.3) THEN
* place an alpha hole in IGAS
            IHOLEATOT = IHOLEATOT + 1
            IHOLEA(IGAS,IHOLEATOT) = 0
            IHOLEAIDENT = IHOLEATOT
          END IF
          IF(IHPVGAS_AB(IGAS,2).EQ.3) THEN
* place an beta hole in IGAS
            IHOLEBTOT = IHOLEBTOT + 1
            IHOLEB(IGAS,IHOLEBTOT) = 0
            IHOLEBIDENT = IHOLEBTOT
          END IF
        END DO
      END IF
*
* Now distribute IP number of particles in IPARTTOT number of GAS
      IF(IP.GE.1) THEN
        DO IGAS = 1,NGAS
          IF(IHPVGAS_AB(IGAS,1).EQ.2) THEN
* place an alpha particle in IGAS
            IPARTATOT = IPARTATOT + 1
            IPARTICLEA(IGAS,IPARTATOT) = 1
          END IF
          IF(IHPVGAS_AB(IGAS,2).EQ.2) THEN
* place an beta particle in IGAS
            IPARTBTOT = IPARTBTOT + 1
            IPARTICLEB(IGAS,IPARTBTOT) = 1
          END IF
* Add identity operator for particles in the valence space
* These are taken from T
          IF(IHPVGAS_AB(IGAS,1).EQ.3) THEN
* place an alpha hole in IGAS
            IPARTATOT = IPARTATOT + 1
            IPARTICLEA(IGAS,IPARTATOT) = 0
            IPARTAIDENT = IPARTATOT
          END IF
          IF(IHPVGAS_AB(IGAS,2).EQ.3) THEN
* place an beta hole in IGAS
            IPARTBTOT = IPARTBTOT + 1
            IPARTICLEB(IGAS,IPARTBTOT) = 0
            IPARTBIDENT = IPARTBTOT
          END IF
        END DO
      END IF
*
* Now create the intermediate supergroups
C     print*,'NAOBEX_TP,NBOBEX_TP',NAOBEX_TP,NBOBEX_TP
C     print*,'IPARTATOT,IPARTBTOT',IPARTATOT,IPARTBTOT
C     print*,IPARTICLEA(1,1),IPARTICLEA(2,1)
C     print*,IPARTICLEB(1,1),IPARTICLEB(2,1)
C     print*,'IHOLEATOT,IHOLEBTOT',IHOLEATOT,IHOLEBTOT
C     print*,IHOLEA(1,1),IHOLEA(2,1)
C     print*,IHOLEB(1,1),IHOLEB(2,1)
      DO IAOBEX =1,NAOBEX_TP
        DO IBOBEX =1,NBOBEX_TP
* List the five cases of (P,H)
C     print*,IPARTICLEA(1,1),IPARTICLEA(2,1)
C     print*,IPARTICLEB(1,1),IPARTICLEB(2,1)
C     print*,IHOLEA(1,1),IHOLEA(2,1)
C     print*,IHOLEB(1,1),IHOLEB(2,1)
          IF(IH.EQ.0.AND.IP.EQ.1) THEN
            CALL IH0IP1(NGAS,IINI,IP,IH,NINTER,INTER,IKRFLIPMAX,
     &                  IAOBEX_TP(1,IAOBEX),IBOBEX_TP(1,IBOBEX),
     &                  IPARTATOT,IPARTBTOT,IPARTICLEA,IPARTICLEB,
     &                  IPARTAIDENT,IPARTBIDENT,
     &                  IALLOWIDENT)
          ELSE IF(IH.EQ.1.AND.IP.EQ.0) THEN
            CALL IH1IP0(NGAS,IINI,IP,IH,NINTER,INTER,IKRFLIPMAX,
     &                  IAOBEX_TP(1,IAOBEX),IBOBEX_TP(1,IBOBEX),
     &                  IHOLEATOT,IHOLEBTOT,IHOLEA,IHOLEB,
     &                  IHOLEAIDENT,IHOLEBIDENT,
     &                  IALLOWIDENT)
          ELSE IF(IH.EQ.1.AND.IP.EQ.1) THEN
            CALL IH1IP1(NGAS,IINI,IP,IH,NINTER,INTER,IKRFLIPMAX,
     &                  IAOBEX_TP(1,IAOBEX),IBOBEX_TP(1,IBOBEX),
     &                  IPARTATOT,IPARTBTOT,IPARTICLEA,IPARTICLEB,
     &                  IPARTAIDENT,IPARTBIDENT,
     &                  IHOLEATOT,IHOLEBTOT,IHOLEA,IHOLEB,
     &                  IHOLEAIDENT,IHOLEBIDENT,
     &                  IALLOWIDENT)
          ELSE IF(IH.EQ.2.AND.IP.EQ.0) THEN
            CALL IH2IP0(NGAS,IINI,IP,IH,NINTER,INTER,IKRFLIPMAX,
     &                  IAOBEX_TP(1,IAOBEX),IBOBEX_TP(1,IBOBEX),
     &                  IHOLEATOT,IHOLEBTOT,IHOLEA,IHOLEB,
     &                  IHOLEAIDENT,IHOLEBIDENT,
     &                  IALLOWIDENT)
          ELSE IF(IH.EQ.0.AND.IP.EQ.2) THEN
            CALL IH0IP2(NGAS,IINI,IP,IH,NINTER,INTER,IKRFLIPMAX,
     &                  IAOBEX_TP(1,IAOBEX),IBOBEX_TP(1,IBOBEX),
     &                  IPARTATOT,IPARTBTOT,IPARTICLEA,IPARTICLEB,
     &                  IPARTAIDENT,IPARTBIDENT,
     &                  IALLOWIDENT)
          ELSE IF(IH.EQ.2.AND.IP.EQ.1) THEN
            CALL IH2IP1(NGAS,IINI,IP,IH,NINTER,INTER,IKRFLIPMAX,
     &                  IAOBEX_TP(1,IAOBEX),IBOBEX_TP(1,IBOBEX),
     &                  IPARTATOT,IPARTBTOT,IPARTICLEA,IPARTICLEB,
     &                  IPARTAIDENT,IPARTBIDENT,
     &                  IHOLEATOT,IHOLEBTOT,IHOLEA,IHOLEB,
     &                  IHOLEAIDENT,IHOLEBIDENT,
     &                  IALLOWIDENT)
          ELSE
            STOP "HOW DID YOU MANAGE TO GET HERE??? "
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Information for intemediate ',IH,IP
        WRITE(6,*) ' =============================== '
        WRITE(6,*) ' Holes and particles '
        WRITE(6,*) ' IHOLEATOT ',IHOLEATOT
        DO I =1,IHOLEATOT
          WRITE(6,*) ' ALPHA HOLE STRING ',I
          WRITE(6,*) (IHOLEA(IGAS,I),IGAS=1,NGAS)
        END DO
        WRITE(6,*) ' IHOLEBTOT ',IHOLEBTOT
        DO I =1,IHOLEBTOT
          WRITE(6,*) ' BETA HOLE STRING ',I
          WRITE(6,*) (IHOLEB(IGAS,I),IGAS=1,NGAS)
        END DO
        WRITE(6,*) ' IPARTATOT ',IPARTATOT
        DO I =1,IPARTATOT
          WRITE(6,*) ' ALPHA PARTICLE STRING ',I
          WRITE(6,*) (IPARTICLEA(IGAS,I),IGAS=1,NGAS)
        END DO
        WRITE(6,*) ' IPARTBTOT ',IPARTBTOT
        DO I =1,IPARTBTOT
          WRITE(6,*) ' BETA PARTICLE STRING ',I
          WRITE(6,*) (IPARTICLEB(IGAS,I),IGAS=1,NGAS)
        END DO
        WRITE(6,*) ' =============================== '
        WRITE(6,*) ' Placement of identity operators '
        WRITE(6,*) ' IHOLEAIDENT ',IHOLEAIDENT
        WRITE(6,*) ' IHOLEBIDENT ',IHOLEBIDENT
        WRITE(6,*) ' IPARTAIDENT ',IPARTAIDENT
        WRITE(6,*) ' IPARTAIDENT ',IPARTBIDENT
        WRITE(6,*) ' =============================== '
        WRITE(6,*) ' Number of intermediates ',NINTER
        IF(IINI.EQ.1) THEN
        CALL WRT_SPOX_TP_CC_KRCC(INTER,NINTER)
        END IF
      END IF
*
      RETURN
      END
*
      SUBROUTINE IH2IP0(NGAS,IINI,IP,IH,NINTER,INTER,IKRFLIPMAX,
     &                  IAOBEX_TP,IBOBEX_TP,
     &                  IHOLEATOT,IHOLEBTOT,IHOLEA,IHOLEB,
     &                  IHOLEAIDENT,IHOLEBIDENT,IALLOWIDENT)
* For 2 holes and 0 particle
* Will only allow one identity operator to pass though
      IMPLICIT REAL*8(A-H,O-Z)
* Output    
      INTEGER INTER(4*NGAS,*)
      INTEGER IHOLEA(NGAS,IHOLEATOT),IHOLEB(NGAS,IHOLEBTOT)
* Input   
      INTEGER IAOBEX_TP(2*NGAS),IBOBEX_TP(2*NGAS)
* Scratch
      INTEGER ITEMP(4*NGAS) 
*           
      IF(IALLOWIDENT.EQ.0) THEN
        IDENTHA = 0
        IDENTHB = 0
      ELSE         
        IDENTHA = 1
        IDENTHB = 1
      END IF
* hole  A B
      DO IAH = 1,IHOLEATOT
*
        IF(IAH.EQ.IHOLEAIDENT) THEN
          IF(IDENTHA.EQ.1) THEN
            CYCLE
          ELSE
            IDENTHA = 1
          END IF
        END IF
*
        DO IBH = 1,IHOLEBTOT
*
          IF(IBH.EQ.IHOLEBIDENT) THEN
            IF(IDENTHB.EQ.1) THEN
              CYCLE
            ELSE
              IDENTHB = 1
            END IF
          END IF
*
          DO IGAS =1,NGAS
* Creator alpha
            ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) +
     &                           IHOLEA(IGAS,IAH)
* Creator Beta
            ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS) +
     &                           IHOLEB(IGAS,IBH)
* Annihilator Alpha
            ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) 
* Annihilator Beta
            ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS)
          END DO
C         CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
          CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
          IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
            NINTER = NINTER + 1
            IF(IINI.EQ.1) THEN
* Store super group
              CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C             CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
            END IF
          END IF
        END DO
      END DO
* Hole A A
      DO IAH = 1,IHOLEATOT
*
        IF(IAH.EQ.IHOLEAIDENT) THEN
          IF(IDENTHA.EQ.1) THEN
            CYCLE
          ELSE
            IDENTHA = 1
          END IF
        END IF
*
        DO IAH2 = 1,IHOLEATOT
*
          IF(IAH2.EQ.IHOLEAIDENT) THEN
            IF(IDENTHA.EQ.1) THEN
              CYCLE
            ELSE
              IDENTHA = 1
            END IF
          END IF
*
          DO IGAS =1,NGAS
* Creator alpha
            ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) +
     &                           IHOLEA(IGAS,IAH)       +
     &                           IHOLEA(IGAS,IAH2)
* Creator Beta
            ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS) 
* Annihilator Alpha
            ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) 
* Annihilator Beta
            ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS)
          END DO
C         CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
          CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
          IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
            NINTER = NINTER + 1
            IF(IINI.EQ.1) THEN
* Store super group
              CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C             CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
            END IF
          END IF
        END DO
      END DO
* Hole B B
      DO IBH = 1,IHOLEBTOT
*
        IF(IBH.EQ.IHOLEBIDENT) THEN
          IF(IDENTHB.EQ.1) THEN
            CYCLE
          ELSE
            IDENTHB = 1
          END IF        
        END IF          
*
        DO IBH2 = 1,IHOLEBTOT
*
          IF(IBH2.EQ.IHOLEBIDENT) THEN
            IF(IDENTHB.EQ.1) THEN
              CYCLE
            ELSE
              IDENTHB = 1
            END IF        
          END IF          
*
          DO IGAS =1,NGAS
* Creator alpha
            ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) 
* Creator Beta
            ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS) +
     &                           IHOLEB(IGAS,IBH)       +
     &                           IHOLEB(IGAS,IBH2)
* Annihilator Alpha
            ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) 
* Annihilator Beta
            ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS)
          END DO
C         CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
          CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
          IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
            NINTER = NINTER + 1
            IF(IINI.EQ.1) THEN
* Store super group
              CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C             CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
            END IF
          END IF
        END DO
      END DO 
*
      RETURN
      END
*
      SUBROUTINE IH0IP2(NGAS,IINI,IP,IH,NINTER,INTER,IKRFLIPMAX,
     &                  IAOBEX_TP,IBOBEX_TP,
     &                  IPARTATOT,IPARTBTOT,IPARTICLEA,IPARTICLEB,
     &                  IPARTAIDENT,IPARTBIDENT,
     &                  IALLOWIDENT)
* For 0 holes and 2 particles
* Will only allow one identity operator to pass though
      IMPLICIT REAL*8(A-H,O-Z)
* Output    
      INTEGER INTER(4*NGAS,*)
      INTEGER IPARTICLEA(NGAS,IPARTATOT),IPARTICLEB(NGAS,IPARTBTOT)
* Input   
      INTEGER IAOBEX_TP(2*NGAS),IBOBEX_TP(2*NGAS)
* Scratch
      INTEGER ITEMP(4*NGAS) 
*           
      IF(IALLOWIDENT.EQ.0) THEN
        IDENTPA = 0
        IDENTPB = 0
      ELSE         
        IDENTPA = 1
        IDENTPB = 1
      END IF
* hole  A B
      DO IAP = 1,IPARTATOT
*
        IF(IAP.EQ.IPARTAIDENT) THEN
          IF(IDENTPA.EQ.1) THEN
            CYCLE
          ELSE
            IDENTPA = 1
          END IF
        END IF
*
        DO IBP = 1,IPARTBTOT
*
          IF(IBP.EQ.IPARTBIDENT) THEN
            IF(IDENTPB.EQ.1) THEN
              CYCLE
            ELSE
              IDENTPB = 1
            END IF
          END IF
*
          DO IGAS =1,NGAS
* Creator alpha
            ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS)
* Creator Beta
            ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS)
* Annihilator Alpha
            ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) +
     &                           IPARTICLEA(IGAS,IAP)
* Annihilator Beta
            ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS) +
     &                           IPARTICLEB(IGAS,IBP)
          END DO
C         CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
          CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
          IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
            NINTER = NINTER + 1
            IF(IINI.EQ.1) THEN
* Store super group
              CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C             CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
            END IF
          END IF
        END DO
      END DO
* Hole A A
      DO IAP = 1,IPARTATOT
*
        IF(IAP.EQ.IPARTAIDENT) THEN
          IF(IDENTPA.EQ.1) THEN
            CYCLE
          ELSE
            IDENTPA = 1
          END IF
        END IF
*
        DO IAP2 = 1,IPARTATOT
*
          IF(IAP2.EQ.IPARTAIDENT) THEN
            IF(IDENTPA.EQ.1) THEN
              CYCLE
            ELSE
              IDENTPA = 1
            END IF
          END IF
*
          DO IGAS =1,NGAS
* Creator alpha
            ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS)
* Creator Beta
            ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS) 
* Annihilator Alpha
            ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) +
     &                           IPARTICLEA(IGAS,IAP)   +
     &                           IPARTICLEA(IGAS,IAP2)
* Annihilator Beta
            ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS)
          END DO
C         CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
          CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
          IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
            NINTER = NINTER + 1
            IF(IINI.EQ.1) THEN
* Store super group
              CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C             CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
            END IF
          END IF
        END DO
      END DO
* Hole B B
      DO IBP = 1,IPARTBTOT
*
        IF(IBP.EQ.IPARTBIDENT) THEN
          IF(IDENTPB.EQ.1) THEN
            CYCLE
          ELSE
            IDENTPB = 1
          END IF        
        END IF          
*
        DO IBP2 = 1,IPARTBTOT
*
          IF(IBP2.EQ.IPARTBIDENT) THEN
            IF(IDENTPB.EQ.1) THEN
              CYCLE
            ELSE
              IDENTPB = 1
            END IF        
          END IF          
*
          DO IGAS =1,NGAS
* Creator alpha
            ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) 
* Creator Beta
            ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS)
* Annihilator Alpha
            ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) 
* Annihilator Beta
            ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS) +
     &                           IPARTICLEB(IGAS,IBP)   +
     &                           IPARTICLEB(IGAS,IBP2)
          END DO
C         CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
          CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
          IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
            NINTER = NINTER + 1
            IF(IINI.EQ.1) THEN
* Store super group
              CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C             CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
            END IF
          END IF
        END DO
      END DO 
*
      RETURN
      END
*
      SUBROUTINE IH2IP1(NGAS,IINI,IP,IH,NINTER,INTER,IKRFLIPMAX,
     &                  IAOBEX_TP,IBOBEX_TP,
     &                  IPARTATOT,IPARTBTOT,IPARTICLEA,IPARTICLEB,
     &                  IPARTAIDENT,IPARTBIDENT,
     &                  IHOLEATOT,IHOLEBTOT,IHOLEA,IHOLEB,
     &                  IHOLEAIDENT,IHOLEBIDENT,IALLOWIDENT)
* For 1 holes and 1 particle
* Will only allow one identity operator to pass though for given
* ha,hb,pa,pb
      IMPLICIT REAL*8(A-H,O-Z)
* Output
      INTEGER INTER(4*NGAS,*)
      INTEGER IHOLEA(NGAS,IHOLEATOT),IHOLEB(NGAS,IHOLEBTOT)
      INTEGER IPARTICLEA(NGAS,IPARTATOT),IPARTICLEB(NGAS,IPARTBTOT)
* Input
      INTEGER IAOBEX_TP(2*NGAS),IBOBEX_TP(2*NGAS)
* Scratch
      INTEGER ITEMP(4*NGAS)
*
      IF(IALLOWIDENT.EQ.0) THEN
        IDENTP  = 0
        IDENTHA = 0
        IDENTHB = 0
      ELSE
        IDENTP  = 1
        IDENTHA = 1
        IDENTHB = 1
      END IF
* hole  A B part A and B
      DO IAH = 1,IHOLEATOT
*
        IF(IAH.EQ.IHOLEAIDENT) THEN
          IF(IDENTHA.EQ.1) THEN
            CYCLE
          ELSE
            IDENTHA = 1
          END IF
        END IF
*
        DO IBH = 1,IHOLEBTOT
*
          IF(IBH.EQ.IHOLEBIDENT) THEN
            IF(IDENTHB.EQ.1) THEN
              CYCLE
            ELSE
              IDENTHB = 1
            END IF
          END IF
* Part A
          DO IAP = 1,IPARTATOT
*
            IF(IAP.EQ.IPARTAIDENT) THEN
              IF(IDENTP.EQ.1) THEN
                CYCLE
              ELSE
                IDENTP = 1
              END IF
            END IF
*
            DO IGAS =1,NGAS
* Creator alpha
              ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) +
     &                             IHOLEA(IGAS,IAH)
* Creator Beta
              ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS) +
     &                             IHOLEB(IGAS,IBH)
* Annihilator Alpha
              ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) +
     &                             IPARTICLEA(IGAS,IAP)
* Annihilator Beta
              ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS)
            END DO
C           CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
            CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
            IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
              NINTER = NINTER + 1
              IF(IINI.EQ.1) THEN
* Store super group
                CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C               CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
              END IF
            END IF
          END DO
* Part B
          DO IBP = 1,IPARTBTOT
*
            IF(IBP.EQ.IPARTBIDENT) THEN
              IF(IDENTP.EQ.1) THEN
                CYCLE
              ELSE
                IDENTP = 1
              END IF
            END IF
*
            DO IGAS =1,NGAS
* Creator alpha
              ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) +
     &                             IHOLEA(IGAS,IAH)
* Creator Beta
              ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS) +
     &                             IHOLEB(IGAS,IBH)
* Annihilator Alpha
              ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS)
* Annihilator Beta
              ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS) +
     &                             IPARTICLEB(IGAS,IBP)
            END DO
C           CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
            CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
            IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
              NINTER = NINTER + 1
              IF(IINI.EQ.1) THEN
* Store super group
                CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C               CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
              END IF
            END IF
          END DO
        END DO
      END DO
* Hole A A Part A and B
      DO IAH = 1,IHOLEATOT
*
        IF(IAH.EQ.IHOLEAIDENT) THEN
          IF(IDENTHA.EQ.1) THEN
            CYCLE
          ELSE
            IDENTHA = 1
          END IF
        END IF
*
        DO IAH2 = 1,IHOLEATOT
*
          IF(IAH2.EQ.IHOLEAIDENT) THEN
            IF(IDENTHA.EQ.1) THEN
              CYCLE
            ELSE
              IDENTHA = 1
            END IF
          END IF
* Part A
          DO IAP = 1,IPARTATOT
*
            IF(IAP.EQ.IPARTAIDENT) THEN
              IF(IDENTP.EQ.1) THEN
                CYCLE
              ELSE
                IDENTP = 1
              END IF
            END IF
*
            DO IGAS =1,NGAS
* Creator alpha
              ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) +
     &                             IHOLEA(IGAS,IAH)       +
     &                             IHOLEA(IGAS,IAH2)
* Creator Beta
              ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS) 
* Annihilator Alpha
              ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) +
     &                             IPARTICLEA(IGAS,IAP)
* Annihilator Beta
              ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS)
            END DO
C           CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
            CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
            IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
              NINTER = NINTER + 1
              IF(IINI.EQ.1) THEN
* Store super group
                CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C               CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
              END IF
            END IF
          END DO
* Part B
          DO IBP = 1,IPARTBTOT
*
            IF(IBP.EQ.IPARTBIDENT) THEN
              IF(IDENTP.EQ.1) THEN
                CYCLE
              ELSE
                IDENTP = 1
              END IF
            END IF
*
            DO IGAS =1,NGAS
* Creator alpha
              ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) +
     &                             IHOLEA(IGAS,IAH)       +
     &                             IHOLEA(IGAS,IAH2)
* Creator Beta
              ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS) 
* Annihilator Alpha
              ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS)
* Annihilator Beta
              ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS) +
     &                             IPARTICLEB(IGAS,IBP)
            END DO
C           CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
            CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
            IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
              NINTER = NINTER + 1
              IF(IINI.EQ.1) THEN
* Store super group
                CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C               CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
              END IF
            END IF
          END DO
        END DO
      END DO
* Hole B B Part A and B
      DO IBH = 1,IHOLEBTOT
*
        IF(IBH.EQ.IHOLEBIDENT) THEN
          IF(IDENTHB.EQ.1) THEN
            CYCLE
          ELSE
            IDENTHB = 1
          END IF
        END IF
*
        DO IBH2 = 1,IHOLEBTOT
*
          IF(IBH2.EQ.IHOLEBIDENT) THEN
            IF(IDENTHB.EQ.1) THEN
              CYCLE
            ELSE
              IDENTHB = 1
            END IF
          END IF
* Part A
          DO IAP = 1,IPARTATOT
*
            IF(IAP.EQ.IPARTAIDENT) THEN
              IF(IDENTP.EQ.1) THEN
                CYCLE
              ELSE
                IDENTP = 1
              END IF
            END IF
*
            DO IGAS =1,NGAS
* Creator alpha
              ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) 
* Creator Beta
              ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS) +
     &                             IHOLEB(IGAS,IBH)       +
     &                             IHOLEB(IGAS,IBH2)
* Annihilator Alpha
              ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) +
     &                             IPARTICLEA(IGAS,IAP)
* Annihilator Beta
              ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS)
            END DO
C           CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
            CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
            IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
              NINTER = NINTER + 1
              IF(IINI.EQ.1) THEN
* Store super group
                CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C               CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
              END IF
            END IF
          END DO
* Part B
          DO IBP = 1,IPARTBTOT
*
            IF(IBP.EQ.IPARTBIDENT) THEN
              IF(IDENTP.EQ.1) THEN
                CYCLE
              ELSE
                IDENTP = 1
              END IF
            END IF
*
            DO IGAS =1,NGAS
* Creator alpha
              ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) 
* Creator Beta
              ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS) +
     &                             IHOLEB(IGAS,IBH)       +
     &                             IHOLEB(IGAS,IBH2)
* Annihilator Alpha
              ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS)
* Annihilator Beta
              ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS) +
     &                             IPARTICLEB(IGAS,IBP)
            END DO
C           CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
            CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
            IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
              NINTER = NINTER + 1
              IF(IINI.EQ.1) THEN
* Store super group
                CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C               CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
              END IF
            END IF
          END DO
        END DO
      END DO
*
      RETURN
      END
*
      SUBROUTINE IH1IP1(NGAS,IINI,IP,IH,NINTER,INTER,IKRFLIPMAX,
     &                  IAOBEX_TP,IBOBEX_TP,
     &                  IPARTATOT,IPARTBTOT,IPARTICLEA,IPARTICLEB,
     &                  IPARTAIDENT,IPARTBIDENT,
     &                  IHOLEATOT,IHOLEBTOT,IHOLEA,IHOLEB,
     &                  IHOLEAIDENT,IHOLEBIDENT,IALLOWIDENT)
* For 1 holes and 1 particle
* Will only allow one identity operator to pass though for given
* ha,hb,pa,pb
      IMPLICIT REAL*8(A-H,O-Z)
* Output
      INTEGER INTER(4*NGAS,*)
      INTEGER IHOLEA(NGAS,IHOLEATOT),IHOLEB(NGAS,IHOLEBTOT)
      INTEGER IPARTICLEA(NGAS,IPARTATOT),IPARTICLEB(NGAS,IPARTBTOT)
* Input
      INTEGER IAOBEX_TP(2*NGAS),IBOBEX_TP(2*NGAS)
* Scratch
      INTEGER ITEMP(4*NGAS)
*
      IF(IALLOWIDENT.EQ.0) THEN
        IDENTPA = 0
        IDENTPB = 0
        IDENTHA = 0
        IDENTHB = 0
      ELSE              
        IDENTPA = 1
        IDENTPB = 1
        IDENTHA = 1
        IDENTHB = 1
      END IF            
*
      DO IAH = 1,IHOLEATOT
*
        IF(IAH.EQ.IHOLEAIDENT) THEN
          IF(IDENTHA.EQ.1) THEN
            CYCLE
          ELSE
            IDENTHA = 1
          END IF
        END IF
*
        DO IAP = 1,IPARTATOT
*
          IF(IAP.EQ.IPARTAIDENT) THEN
            IF(IDENTPA.EQ.1) THEN
              CYCLE
            ELSE
              IDENTPA = 1
            END IF
          END IF
*
          DO IGAS =1,NGAS
* Creator alpha
            ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) +
     &                           IHOLEA(IGAS,IAH)
* Creator Beta
            ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS)
* Annihilator Alpha
            ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) +
     &                           IPARTICLEA(IGAS,IAP)
* Annihilator Beta
            ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS)
          END DO
C         CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
          CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
          IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
            NINTER = NINTER + 1
            IF(IINI.EQ.1) THEN
* Store super group
              CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C             CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
            END IF
          END IF
        END DO
      END DO
*
      DO IAH = 1,IHOLEATOT
*
        IF(IAH.EQ.IHOLEAIDENT) THEN
          IF(IDENTHA.EQ.1) THEN
            CYCLE
          ELSE
            IDENTHA = 1 
          END IF
        END IF
*
        DO IBP = 1,IPARTBTOT
*
          IF(IBP.EQ.IPARTBIDENT) THEN
            IF(IDENTPB.EQ.1) THEN
              CYCLE
            ELSE
              IDENTPB = 1
            END IF
          END IF
*
          DO IGAS =1,NGAS
* Creator alpha
            ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) +
     &                           IHOLEA(IGAS,IAH)
* Creator Beta
            ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS)
* Annihilator Alpha
            ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) 
* Annihilator Beta
            ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS) +
     &                           IPARTICLEB(IGAS,IBP)
          END DO
C         CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
          CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
          IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
            NINTER = NINTER + 1
            IF(IINI.EQ.1) THEN
* Store super group
              CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C             CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
            END IF
          END IF
        END DO                   
      END DO
*
      DO IBH = 1,IHOLEBTOT
*
        IF(IBH.EQ.IHOLEBIDENT) THEN
          IF(IDENTHB.EQ.1) THEN
            CYCLE
          ELSE
            IDENTHB = 1
          END IF
        END IF
*
        DO IAP = 1,IPARTATOT
*    
          IF(IAP.EQ.IPARTAIDENT) THEN
            IF(IDENTPA.EQ.1) THEN
              CYCLE
            ELSE
              IDENTPA = 1
            END IF
          END IF 
*       
          DO IGAS =1,NGAS
* Creator alpha         
            ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) 
* Creator Beta
            ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS) +
     &                           IHOLEB(IGAS,IBH)
* Annihilator Alpha
            ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) +
     &                           IPARTICLEA(IGAS,IAP)
* Annihilator Beta
            ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS)
          END DO
C         CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
          CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
          IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
            NINTER = NINTER + 1
            IF(IINI.EQ.1) THEN
* Store super group
              CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C             CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
            END IF
          END IF
        END DO
      END DO
*
      DO IBH = 1,IHOLEBTOT
*
        IF(IBH.EQ.IHOLEBIDENT) THEN
          IF(IDENTHB.EQ.1) THEN  
            CYCLE
          ELSE
            IDENTHB = 1
          END IF
        END IF
*
        DO IBP = 1,IPARTBTOT
*    
          IF(IBP.EQ.IPARTBIDENT) THEN
            IF(IDENTPB.EQ.1) THEN
              CYCLE
            ELSE
              IDENTPB = 1
            END IF
          END IF 
*       
          DO IGAS =1,NGAS
* Creator alpha         
            ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) 
* Creator Beta
            ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS) +
     &                           IHOLEB(IGAS,IBH)
* Annihilator Alpha
            ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) 
* Annihilator Beta
            ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS) +
     &                           IPARTICLEB(IGAS,IBP)
          END DO
C         CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
          CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
          IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
            NINTER = NINTER + 1
            IF(IINI.EQ.1) THEN
* Store super group
              CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C             CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
            END IF
          END IF
        END DO
      END DO
*
      RETURN
      END
*
      SUBROUTINE IH1IP0(NGAS,IINI,IP,IH,NINTER,INTER,IKRFLIPMAX,
     &                  IAOBEX_TP,IBOBEX_TP,
     &                  IHOLEATOT,IHOLEBTOT,IHOLEA,IHOLEB,
     &                  IHOLEAIDENT,IHOLEBIDENT,IALLOWIDENT)
* For 1 holes and 0 particle
* Will only allow one identity operator to pass though
      IMPLICIT REAL*8(A-H,O-Z)
* Output
      INTEGER INTER(4*NGAS,*)
      INTEGER IHOLEA(NGAS,IHOLEATOT),IHOLEB(NGAS,IHOLEBTOT)
* Input
      INTEGER IAOBEX_TP(2*NGAS),IBOBEX_TP(2*NGAS)
* Scratch
      INTEGER ITEMP(4*NGAS)
*
      IF(IALLOWIDENT.EQ.0) THEN
        IDENT = 0
      ELSE
        IDENT = 1
      END IF
*
      DO IA = 1,IHOLEATOT
        IF(IA.EQ.IHOLEAIDENT) THEN
          IF(IDENT.EQ.1) THEN
            CYCLE
          ELSE
            IDENT = 1
          END IF
        END IF
        DO IGAS =1,NGAS
* Creator alpha
          ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS) +
     &                         IHOLEA(IGAS,IA)
* Creator Beta
          ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS)
* Annihilator Alpha
          ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) 
* Annihilator Beta
          ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS)
        END DO
C       CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
        CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
        IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
          NINTER = NINTER + 1
          IF(IINI.EQ.1) THEN
* Store super group
            CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C           CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
          END IF
        END IF
      END DO
*
      DO IB = 1,IHOLEBTOT
        IF(IB.EQ.IHOLEBIDENT) THEN
          IF(IDENT.EQ.1) THEN
            CYCLE
          ELSE
            IDENT = 1
          END IF
        END IF
        DO IGAS =1,NGAS
* Creator alpha
            ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS)
* Creator Beta
            ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS) +
     &                           IHOLEB(IGAS,IB)
* Annihilator Alpha
            ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS)
* Annihilator Beta
            ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS)
        END DO
C       CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
        CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
        IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
          NINTER = NINTER + 1
          IF(IINI.EQ.1) THEN
* Store super group
            CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C           CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
          END IF
        END IF
      END DO
*
      RETURN
      END
*
      SUBROUTINE IH0IP1(NGAS,IINI,IP,IH,NINTER,INTER,IKRFLIPMAX,
     &                  IAOBEX_TP,IBOBEX_TP,
     &                  IPARTATOT,IPARTBTOT,IPARTICLEA,IPARTICLEB,
     &                  IPARTAIDENT,IPARTBIDENT,IALLOWIDENT)
* For 0 holes and 1 particle
* Will only allow one identity operator to pass though
      IMPLICIT REAL*8(A-H,O-Z)
* Output
      INTEGER INTER(4*NGAS,*)
      INTEGER IPARTICLEA(NGAS,IPARTATOT),IPARTICLEB(NGAS,IPARTBTOT)
* Input
      INTEGER IAOBEX_TP(2*NGAS),IBOBEX_TP(2*NGAS)
* Scratch
      INTEGER ITEMP(4*NGAS)
*
      IF(IALLOWIDENT.EQ.0) THEN
        IDENT = 0
      ELSE
        IDENT = 1
      END IF
*
      DO IA = 1,IPARTATOT
        IF(IA.EQ.IPARTAIDENT) THEN
          IF(IDENT.EQ.1) THEN
            CYCLE
          ELSE
            IDENT = 1
          END IF
        END IF
        DO IGAS =1,NGAS
* Creator alpha
          ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS)
* Creator Beta
          ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS)
* Annihilator Alpha
          ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS) +
     &                         IPARTICLEA(IGAS,IA)
* Annihilator Beta
          ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS)
        END DO
C       CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
        CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
        IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
          NINTER = NINTER + 1
          IF(IINI.EQ.1) THEN
* Store super group
            CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C           CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
          END IF
        END IF
      END DO
*
      DO IB = 1,IPARTBTOT
        IF(IB.EQ.IPARTBIDENT) THEN
          IF(IDENT.EQ.1) THEN
            CYCLE
          ELSE
            IDENT = 1
          END IF
        END IF
        DO IGAS =1,NGAS
* Creator alpha
            ITEMP(IGAS+0*NGAS) = IAOBEX_TP(IGAS+0*NGAS)
* Creator Beta
            ITEMP(IGAS+1*NGAS) = IBOBEX_TP(IGAS+0*NGAS)
* Annihilator Alpha
            ITEMP(IGAS+2*NGAS) = IAOBEX_TP(IGAS+1*NGAS)
* Annihilator Beta
            ITEMP(IGAS+3*NGAS) = IBOBEX_TP(IGAS+1*NGAS) +
     &                           IPARTICLEB(IGAS,IB)
        END DO
C       CALL WRT_SPOX_TP_CC_KRCC(ITEMP,1)
        CALL CANITBE(ITEMP,IANSWER,IKRFLIPMAX,IH,IP)
        IF(IANSWER.EQ.1) THEN
* We have found an intemediate super group
          NINTER = NINTER + 1
          IF(IINI.EQ.1) THEN
* Store super group
            CALL ICOPVE(ITEMP,INTER(1,NINTER),4*NGAS)
C           CALL WRT_SPOX_TP_CC_KRCC(INTER(1,NINTER),1)
          END IF
        END IF
      END DO
*
      RETURN
      END
*
      SUBROUTINE CANITBE(IOP,IANSWER,IKRFLIPMAX,IH,IP)
*
* Checks if an operator can exist
* We need ::
* Conservation of particle number
* Not higher particle rank than T
* Stay inside GAS
* Correct number of holes and particles
* Cannot exceede Kramers flip 
* Notice this is passed down since a restriction 
* on the Kramers Flip for T will not give the same 
* restriction on the way for an Intermediate operator
* In general will the max flip for an intermediate be
* NELEC >= KF INTER >= (N+1)*T 
* Where N is the number of commutators since
* any comtraction has Delta MK=0 so the total
* kramers flip of the resultant operator will 
* just be the sum of the two ooperators contracted
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "cgas.inc"
#include "crun.inc"
#include "dgroup.h"
*
      INTEGER IOP(4*NGAS)
*
      NTEST = 00
*
      IANSWER = 0
*
* Conservation of particle number
      ICONS = 0
      DO IGAS =1,NGAS
        ICONS = IOP(IGAS+0*NGAS) + 
     &          IOP(IGAS+1*NGAS) -
     &          IOP(IGAS+2*NGAS) -
     &          IOP(IGAS+3*NGAS) +
     &          ICONS
      END DO
      IF(ICONS.NE.0) RETURN
      IF(NTEST.EQ.100) WRITE(6,*) " we have conservation "
*
* Particle rank cannot excede max excitation level
      IPART = 0
      DO IGAS =1,NGAS
        IPART = IOP(IGAS+0*NGAS) + 
     &          IOP(IGAS+1*NGAS) +
     &          IOP(IGAS+2*NGAS) +
     &          IOP(IGAS+3*NGAS) +
     &          IPART
      END DO
      IF(IPART.GT.2*MX_EXC_LEVEL) RETURN
      IF(NTEST.EQ.100) WRITE(6,*) "not too highly excited",IPART
*
* Must contain the correct number of holes and particles
      NHOLE = 0
      NPART = 0
      DO IGAS =1,NGAS
        NHOLEAADD = IOP(IGAS+0*NGAS) 
        NHOLEBADD = IOP(IGAS+1*NGAS)
        NPARTAADD = IOP(IGAS+2*NGAS)
        NPARTBADD = IOP(IGAS+3*NGAS)
        DO IAB = 1,2
          IF(IHPVGAS_AB(IGAS,IAB).EQ.1) THEN
            IF(IAB.EQ.1) THEN
              NHOLE = NHOLE + NHOLEAADD
            ELSE
              NHOLE = NHOLE + NHOLEBADD
            END IF 
          ELSE IF(IHPVGAS_AB(IGAS,IAB).EQ.2) THEN
            IF(IAB.EQ.1) THEN
              NPART = NPART + NPARTAADD
            ELSE
              NPART = NPART + NPARTBADD
            END IF
          ELSE IF(IHPVGAS_AB(IGAS,IAB).EQ.3) THEN
            IF(IAB.EQ.1) THEN
              NHOLE = NHOLE + NHOLEAADD
              NPART = NPART + NPARTAADD
            ELSE
              NHOLE = NHOLE + NHOLEBADD
              NPART = NPART + NPARTBADD
            END IF
          ELSE
            STOP 'THERE IS NO WAY YOU COULD GET HERE'
          END IF
        END DO
      END DO
      IF(IH.NE.NHOLE.AND.IP.NE.NPART) RETURN
      IF(NTEST.EQ.100) WRITE(6,*) " right number of holes ",NHOLE,NPART
      
*
* Is this an allowed operator in the GAS.
* For all indices yet to be contracted rules are different
*
* Should still check if valence space (consider this)
      DO IGAS=1,NGAS
        ICOMP = 0
        DO J=1,NFSYM
          ICOMP = ICOMP + NGSSH(J,IGAS)
        END DO
        IF(IHPVGAS(IGAS).NE.1) THEN
          IF(IOP(IGAS+0*NGAS).GT.ICOMP) RETURN
          IF(IOP(IGAS+1*NGAS).GT.ICOMP) RETURN
        END IF
        IF(IHPVGAS(IGAS).NE.2) THEN
          IF(IOP(IGAS+2*NGAS).GT.ICOMP) RETURN
          IF(IOP(IGAS+3*NGAS).GT.ICOMP) RETURN
        END IF
      END DO
      IF(NTEST.EQ.100) WRITE(6,*) " operator allowed for GAS"
*
* Kramers flipping of operator
      IKRFLIP = 0
      DO IGAS =1,NGAS
        IKRFLIP = IOP(IGAS+0*NGAS) -
     &            IOP(IGAS+1*NGAS) -
     &            IOP(IGAS+2*NGAS) +
     &            IOP(IGAS+3*NGAS) +
     &            IKRFLIP
      END DO
      IF(ABS(IKRFLIP).GT.IKRFLIPMAX) RETURN
      IF(NTEST.EQ.100) WRITE(6,*) ' I flip ok ',IKRFLIP/2
*
* Even or Odd Kramers flip
      IF(NZ.LT.4) THEN
        IKRFLIP = IKRFLIP/2
        CALL EVENODD(ITEST,ABS(IKRFLIP))
        IF(ITEST.EQ.1) RETURN
        IF(NTEST.EQ.100) WRITE(6,*) ' even KRAMERS FLIP'
      END IF
*
* An operator has been found
      IANSWER = 1
      IF(NTEST.EQ.100) WRITE(6,*) ' An operator has been found '
*
      RETURN
      END
*
      SUBROUTINE INI_CC_AMP_KRCC(CC)
*
* Initialize Coupled Cluster amplitudes
*
* IFORM = 1 => Set to zero
* IFORM = 2 => Read in from LU_CCAMP
*
* Jeppe Olsen, Summer of 98
*
* Made automated restart from LU_CCAMP (CCAMP) if it is there
* Also checks for formatted or unformatted file
*
* Lasse, Summer of 08
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "crun.inc"
*(    ^contains number of single and double excitations and N_CC_AMP)
#include "clunit.inc"
*. Amplitudes to be defined
      DIMENSION CC(*)
      LOGICAL YN
      CHARACTER(LEN=2) UNFORM
*
* for real and complex
      IF(IRECOM.EQ.1) THEN
        IMULTFAC = 1
      ELSE
        IMULTFAC = 2
      END IF
*
* Automated restart if LU_CCAMP (CCAMP) exicts
      INQUIRE(FILE='CCAMP',EXIST=YN)
      IF(YN) THEN
        IFORM = 2
        INQUIRE(FILE='CCAMP',UNFORMATTED=UNFORM)
        IF(UNFORM(1:2).EQ.'UN') THEN
          I_FORMATTED = 0
        ELSE
          I_FORMATTED = 1
        END IF
      ELSE
        IFORM = 1
      END IF
*
      IF(IFORM.EQ.1) THEN
        ZERO = 0.0D0
        CALL SETVEC(CC,ZERO,IMULTFAC*N_CC_AMP)
        WRITE(6,*) ' Initial set of amplitudes set to zero '
C       DO I=1,N_CC_AMP
C         CC(I) = I
C         print*,'CC,I',CC(I),I
C       END DO
*
      ELSE IF(IFORM.EQ.2) THEN
        ZERO = 0.0D0
        CALL SETVEC(CC,ZERO,IMULTFAC*N_CC_AMP)
        WRITE(6,*) ' Reading in CC amplitudes from ', LU_CCAMP
        IF(I_FORMATTED.EQ.1) THEN
*. Formatted
          OPEN(Unit=LU_CCAMP,File='CCAMP',Status='UNKNOWN',
     &     Form='FORMATTED')
          CALL REWINO(LU_CCAMP)
          READ(LU_CCAMP,*) N_CC_AMPP
          N_CC_AMP_READ = MIN(N_CC_AMPP,N_CC_AMP)
          DO I = 1, IMULTFAC*N_CC_AMP_READ
            READ(LU_CCAMP,*) CC(I)
          END DO
        ELSE
*. Unformatted
          OPEN(Unit=LU_CCAMP,File='CCAMP',Status='UNKNOWN',
     &     Form='UNFORMATTED')
          CALL REWINO(LU_CCAMP)
          READ(LU_CCAMP) N_CC_AMPP
          N_CC_AMP_READ = MIN(N_CC_AMPP,N_CC_AMP)
          READ(LU_CCAMP) (CC(I),I=1,IMULTFAC*N_CC_AMP_READ)
        END IF
          CLOSE(LU_CCAMP)
      ELSE
        WRITE(6,*) ' Unknown parameter in INI_CC_AMP ',IFORM
        STOP       ' Unknown parameter in INI_CC_AMP '
      END IF
*
      IF(IFORM.EQ.2) THEN
        WRITE(6,*) ' A set of amplitudes have been found '
        WRITE(6,*) ' This is a restarted calculation '
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Initial set of amplitudes '
        WRITE(6,*) ' ========================= '
        WRITE(6,*)
        CALL WRT_CC_VEC(CC,6)
      END IF
*
      RETURN
      END
*
      SUBROUTINE OCC_TO_GRP_KRCC(IOCC,IGRP,IWAY)
*. Translate between occupation and group labels of supergroup
*
* IWAY = 1 => OCC to group
*      = 2 => Group to occ
*
* Number of strings in supergroup with given symmetry and number of
* electrons in each space
*
* Jeppe Olsen, March 1999
*
#include "implicit.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Local scratch
C     INTEGER JGRP(MXPNGAS)
*. Specific input/output
      INTEGER IOCC(NGAS),IGRP(NGAS)
*
      NTEST = 00
*
      IF(IWAY.EQ.1) THEN
*. Occupation => Group number
        DO IOBTP = 1, NGAS
          JJGRP = 0
          DO KGRP = IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
            IF(NELFGP(KGRP).EQ.IOCC(IOBTP)) JJGRP = KGRP
          END DO
          IGRP(IOBTP) = JJGRP
*
          IF(JJGRP.EQ.0.or.NTEST.ge.10) THEN
C          WRITE(6,*) ' Group not included in list '
           WRITE(6,*) ' Input occupations : '
           CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
           write(6,*) 'Current orbital type :',IOBTP
           write(6,*) 'IBGPSTR array :'
           call iwrtma(IBGPSTR,1,NGAS,1,NGAS)
           write(6,*) 'NGPSTR array :'
           call iwrtma(NGPSTR,1,NGAS,1,NGAS)
           write(6,*) 'Group range ',
     &                 IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
           write(6,*) 'List of electron occupations: '
           do K=IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
             write(6,*) 'NELFGP(',K,') =', NELFGP(K)
           end do
           if (JJGRP.eq.0) STOP ' Group not included in list '
          END IF
        END DO
*
      ELSE
*. Group => Occupation
        DO IOBTP = 1, NGAS
          IOCC(IOBTP) = NELFGP(IGRP(IOBTP))
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation and corresponding group array '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
        CALL IWRTMA(IGRP,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
*
      SUBROUTINE OCC_TO_GRP_CC_KRCC(IOCC,IGRP,IWAY)
*. Translate between occupation and group labels of supergroup
*
* IWAY = 1 => OCC to group
*      = 2 => Group to occ
*
* Number of strings in supergroup with given symmetry and number of
* electrons in each space
*
* Jeppe Olsen, March 1999
*
#include "implicit.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Local scratch
C     INTEGER JGRP(MXPNGAS)
*. Specific input/output
      INTEGER IOCC(NGAS),IGRP(NGAS)
*
      NTEST = 00
*
      IF(IWAY.EQ.1) THEN
*. Occupation => Group number
        DO IOBTP = 1, NGAS
          JJGRP = 0
          DO KGRP = IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
            IF(NELFGP(KGRP).EQ.IOCC(IOBTP)) JJGRP = KGRP
          END DO
          IGRP(IOBTP) = JJGRP
*
          IF(JJGRP.EQ.0) THEN
           WRITE(6,*) ' Group not included in list '
           WRITE(6,*) ' GAS space with problem : ', IOBTP
           WRITE(6,*) ' Input occupations : '
           CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
           print*,'KGRP range ',IBGPSTR(IOBTP), IBGPSTR(IOBTP)
     &                          + NGPSTR(IOBTP)-1
           STOP       ' Group not included in list '
          END IF
        END DO
*
      ELSE
*. Group => Occupation
        DO IOBTP = 1, NGAS
          IOCC(IOBTP) = NELFGP(IGRP(IOBTP))
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation and corresponding group array '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
        CALL IWRTMA(IGRP,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
*
      SUBROUTINE NST_SPGRP_KRCC(NIGRP,IGRP,ISM_TOT,NSMST,NSTRIN,NDIST,
     &                         IAB)
*
* Number of strings for given combination of groups and
* symmetry.
*
*. Input
*
*
*   NIGRP : Number of active groups
*   IGRP : The active groups
*   ISM_TOT : Total symmetry of supergroup
*   NSMST   : Number of string symmetries
*   IAB : = 1 => unbarred operators
*         = 2 => barred operators
*
*. Output
*
*  NSTRIN : Number of strings with symmetry ISM_TOT
*  NDIST  : Number of symmetry distributions
*
* Jeppe Olsen, September 1997
*
*. Double group version
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
*. Specific Input
      DIMENSION IGRP(NIGRP)
*. General input
#include "gasstr.inc"
C     DIMENSION NSTSGP(NSMST,*)
*. Scratch
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' =========================='
        WRITE(6,*) ' NST_SPGRP_DBG is speaking '
        WRITE(6,*) ' =========================='
*
        WRITE(6,*) ' Supergroup in action : '
        WRITE(6,'(A,I3  )') ' Number of active spaces ', NIGRP
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
      END IF
*
      IF(NIGRP.EQ.0) THEN
        IF(ISM_TOT.EQ.1) THEN
          LENGTH = 1
        ELSE
          LENGTH = 0
        END IF
        GOTO 1001
      END IF
*
*. Set up min and max values for symmetries
      CALL MINMAX_FOR_SYM_DIST_KRCC(NIGRP,IGRP,MNSM,MXSM,NDISTX,IAB)
*. Loop over symmetry distributions
      IFIRST = 1
      LENGTH = 0
      NDIST = 0
*. Last group with symmetry differing from total symmetric
      NIGRPL = 1
      DO JGRP = 1, NIGRP
        IF(MXSM(JGRP).GT.1) NIGRPL = JGRP
      END DO
*. Number of strings in groups after NIGRPL
      NSTRL = 1
      DO JGRP = NIGRPL+1,NIGRP
        IF(IAB.EQ.1) THEN
          NSTRL = NSTRL*NSTFSMGP(1,IGRP(JGRP))
        ELSE
          NSTRL = NSTRL*NSTFSMGP2(1,IGRP(JGRP))
        END IF
      END DO
 1000 CONTINUE
*. Next symmetry distribution
* \/  NEXT_SYM_DISTR_REL should work for double group symmetry
        CALL NEXT_SYM_DISTR_KRCC(NIGRPL,MNSM,MXSM,ISM,ISM_TOT,IFIRST,
     &                          NONEW)
        IF(NONEW.EQ.0) THEN
          LDIST = NSTRL
          DO JGRP = 1, NIGRPL
            IF(IAB.EQ.1) THEN
              LDIST = LDIST*NSTFSMGP(ISM(JGRP),IGRP(JGRP))
            ELSE
              LDIST = LDIST*NSTFSMGP2(ISM(JGRP),IGRP(JGRP))
            END IF
          END DO
          LENGTH = LENGTH + LDIST
          NDIST = NDIST + 1
      GOTO 1000
        END IF
*
 1001 CONTINUE
      NSTRIN = LENGTH
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings obtained ', LENGTH
        WRITE(6,*) ' Number of symmetry-distributions',NDIST
      END IF
*
      RETURN
      END
*
      SUBROUTINE MINMAX_FOR_SYM_DIST_KRCC(NIGRP,IGRP,MNVAL,MXVAL,NDIST,
     &                                   IAB)
*
* A combination of NIGRP groups are given (IGRP)
*. Find MIN and MAX for symmetry in each group
*
* Jeppe Olsen, September 1997
*              April 1998     From  MINMAX_SM_GP
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Include blocks
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "csm.inc"
*. Input
      DIMENSION IGRP(NIGRP)
*.Output
      DIMENSION MNVAL(NIGRP),MXVAL(NIGRP)
*. Local scratch
      DIMENSION LSMGP(MXPOBS,MXPNGAS)
*
      NTEST = 0000
      IF(NTEST.GE.100) WRITE(6,*) ' >> Entering MINMAX_... <<'
*
      DO JGRP = 1, NIGRP
C Requires setting up of  MINMAX_SM_GP and  MINMAX_SM_GP2
        IF(IAB.EQ.1) THEN
          MNVAL(JGRP) = MINMAX_SM_GP(1,IGRP(JGRP))
          MXVAL(JGRP) = MINMAX_SM_GP(2,IGRP(JGRP))
        ELSE
          MNVAL(JGRP) = MINMAX_SM_GP2(1,IGRP(JGRP))
          MXVAL(JGRP) = MINMAX_SM_GP2(2,IGRP(JGRP))
        END IF
      END DO
*. Total number of symmetry distributions
      NDIST = 1
      DO JGRP = 1, NIGRP
        NDIST = NDIST*(MXVAL(JGRP)-MNVAL(JGRP)+1)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Group combination : '
        WRITE(6,'(5X,10I3)') (IGRP(JGRP),JGRP=1, NIGRP)
        WRITE(6,*)
        WRITE(6,*) ' Group Minsym Maxsym'
        WRITE(6,*) ' ==================='
        DO JGRP = 1, NIGRP
          WRITE(6,'(3I6)') IGRP(JGRP),MNVAL(JGRP),MXVAL(JGRP)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' Total number of distributions', NDIST
      END IF
*
      IF(NTEST.GE.1000) WRITE(6,*) ' >> Leaving MINMAX_... <<'
*
      RETURN
      END
*
      SUBROUTINE NEXT_SYM_DISTR_KRCC(NGAS,MINVAL,MAXVAL,
     &                              ISYM,ISYM_TOT,IFIRST,NONEW)
*
* Obtain next distribution of symmetries with given total
* Symmetry.
*
* Loop over first NGAS-1 spaces are performed, and the symmetry
* of the last space is then fixed by the required total sym
*
* Jeppe Olsen, Sept 97
* Obtain next distribution of symmetries with given total
* Symmetry.
*
* Loop over first NGAS-1 spaces are performed, and the symmetry
* of the last space is then fixed by the required total sym
*
* Jeppe Olsen, Sept 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION MINVAL(NGAS),MAXVAL(NGAS)
*. Input and output
      DIMENSION ISYM(NGAS)
*. Symmetry of first NGAS -1 spaces
      IF(IFIRST.EQ.1) THEN
        DO IGAS = 1, NGAS-1
          ISYM(IGAS) = MINVAL(IGAS)
        END DO
        IFIRST = 0
        NONEW = 0
      ELSE
C            NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
        CALL NXTNUM3(ISYM,NGAS-1,MINVAL,MAXVAL,NONEW)
      END IF
*. Symmetry of last space
      IF(NONEW.EQ.0) THEN
        JSYM = 1
        DO IGAS = 1, NGAS-1
          CALL SYMCOM_KRCC(3,0,JSYM,ISYM(IGAS),KSYM)
          JSYM = KSYM
        END DO
        CALL SYMCOM_KRCC(2,0,JSYM,ISYM(NGAS),ISYM_TOT)
C      write(6,*) ' after SYMCOM : JSYM, ISYM(NGAS), ISYM_TOT'
C      write(6,*) JSYM,ISYM(NGAS),ISYM_TOT
*
C       IF(MINVAL(NGAS).GT.ISYM(NGAS).OR.
C    &     MAXVAL(NGAS).LT.ISYM(NGAS)    )NONEW = 1
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(NONEW.EQ.1) THEN
         WRITE(6,*) ' No new symmetry distributions '
        ELSE
         WRITE(6,*) ' Next symmetry distribution '
         CALL IWRTMA(ISYM,1,NGAS,1,NGAS)
        END IF
      END IF
*
      RETURN
      END
*
      SUBROUTINE SYMCOM_KRCC(ITASK,IOBJ,I1,I2,I12)
*
* Symmetries I1,I2,I12 are related as
* I1 I2 = I12
* IF(ITASK = 1 ) I2 and I12 are known, find I1
* IF(ITASK = 2 ) I1 and I12 are known, find I2
* IF(ITASK = 3 ) I1 and I2 are known , find I12
*
* IOBJ = 1 : I1,I2 are strings I12 determinant
* ( Other things can follow )
* IOBJ = 2 : I1,I2,I3 are externals
* IOBJ = 3 : I1 is an external, I2,I3 are dets
* IOBJ = 4 : I1 is orbital, I2 is string,l, I12 is string
* IOBJ = 5 : I1 is single excitation, I2 is string,l, I12 is string
* IOBJ = 6 : I1 is orbital, I2 is Orbital I12 is single excitation
*
* If obtained symmetry I1 or I2 is outside bounds,
* zero is returned.
*
* Jeppe Olsen , Spring of 1991
*
* ================
*. Driver routine
* ================
#include "implicit.inc"
#include "mxpdim.inc"
#include "lucinp.inc"
*
      NTEST = 0000
*
      if (NTEST.ge.2000) then
        write(6,*)
        write(6,*) 'symcom test:'
        write(6,*) 'ITASK : ',ITASK
        write(6,*) 'IOBJ  : ',IOBJ
        write(6,*) 'I1    : ',I1
        write(6,*) 'I2    : ',I2
        write(6,*) 'I12   : ',I12
      end if
*
      IF(PNTGRP.EQ.1) THEN
        CALL SYMCM1(ITASK,IOBJ,I1,I2,I12)
      ELSE IF(PNTGRP.GE.2.AND.PNTGRP.LE.4) THEN
        CALL SYMCM2(ITASK,IOBJ,I1,I2,I12)
      ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
C       CALL SYMCMDG_KRCC(ITASK,IOBJ,I1,I2,I12)
        CALL SYMCMDG_KRCC(ITASK,I1,I2,I12)
      ELSE
        WRITE(6,*) ' PNTGRP parameter out of bounds ', PNTGRP
        WRITE(6,*) ' Enforced stop in SYMCOM '
        Call Abend1( 11 )
      END IF
*
      RETURN
      END
*
      SUBROUTINE SYMCMDG_KRCC(ITASK,I1,I2,I12)
*
* Symmetries I1,I2,I12 are related as
* I1*I2 = 12
* IF(ITASK = 1 ) I2 and I12 are known, find I1
* IF(ITASK = 2 ) I1 and I12 are known, find I1
* IF(ITASK = 3 ) I1 and I2 are known , find I12
*
* Double group version, written for compatibility with general symmetry
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "symm.inc"
*
      IF(ITASK.EQ.1) THEN
        I1  = IDBGMULT(INVELM(I2),I12)
      ELSE IF(ITASK.EQ.2) THEN
        I2  = IDBGMULT(INVELM(I1),I12)
      ELSE IF (ITASK.EQ.3) THEN
        I12 = IDBGMULT(I1,I2)
      END IF
*
      RETURN
      END
*
      SUBROUTINE DIAG_EXC_KRCC(NCA,NCB,NAA,NAB,NGAS,IDIAG)
*
* Check if alpha and betaexcitation parts are identical
*
*. Output
* IDIAG = 0 : Not diagonal in alpha, beta
* IDIAG = 1 : Diag in alpha, beta
*
*
#include "implicit.inc"
*. Input
      INTEGER NCA(NGAS),NCB(NGAS),NAA(NGAS),NAB(NGAS)
*
      IDIAG = 1
      DO IGAS = 1, NGAS
        IF(NCA(IGAS).NE.NCB(IGAS)) IDIAG = 0
        IF(NAA(IGAS).NE.NAB(IGAS)) IDIAG = 0
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NCA, NCB, NAA, NAB = '
        CALL IWRTMA(NCA,1,NGAS,1,NGAS)
        CALL IWRTMA(NCB,1,NGAS,1,NGAS)
        CALL IWRTMA(NAA,1,NGAS,1,NGAS)
        CALL IWRTMA(NAB,1,NGAS,1,NGAS)
        WRITE(6,*) ' IDIAG = ', IDIAG
      END IF
*
      RETURN
      END
*
      SUBROUTINE AB_COMP_FOR_OBOP_KRCC(NCREA,NANNI,MS2TOT,IMS2_REF,
     &                                NCOMP,IABCOMP,NAEL,NBEL,IFLAG)
*
* Spinprojection components of orbital operator containing
* NCREA creation operators and NANNI annihilation operators.
*
* (Modified SPNCOM routine )
*
* Jeppe Olsen,  Summer of 99
*
#include "implicit.inc"
      INTEGER ADD
#include "mxpdim.inc"
*. Output
      DIMENSION IABCOMP(NCREA+NANNI,*)
*. Local scratch
      DIMENSION IWORK(2**MXOPEN)
*
      NTEST = 000
C     WRITE(6,*) ' IFLAG, NAEL, NBEL at start of AB... ',
C    &             IFLAG, NAEL, NBEL
      NCOMP=0
*
* combinations are considered as binary numbers,1=alpha,0=beta
*
      NOPEN = NCREA+NANNI
      MX=2 ** NOPEN
C      print*,'NOPEN,MX,NCREA,NANNI,MXPORB',NOPEN,MX,NCREA,NANNI,MXPORB
* Lasse additional check start
      IF(MX.LE.2**MXOPEN) THEN
        CALL ISETVC(IWORK,0,MX)
      ELSE
        WRITE(6,*) 'MX,2**MXOPEN',MX,2**MXOPEN
        WRITE(6,*) 'AB_COMP_FOR_OBOP_REL IN TROUBLE'
        WRITE(6,*) 'PLEASE INCREASE MXOPEN'
        CALL ABEND('AB_COMP_FOR_OBOP_REL IN TROUBLE')
      END IF
* Lasse additional check end
*
* Loop over all possible binary numbers
      DO 200 I=1,MX
C.. 1 : NEXT BINARY NUMBER
        ADD=1
        J=0
  190   CONTINUE
        J=J+1
* lasse addition to prevent random 1 in iwork(j)
* in- or outside the array
        IF(J.GT.2**MXOPEN.OR. J.GT.MX) THEN
          ADD=0
          stop 'outside array'
          GOTO 200
        ELSE IF(IWORK(J).EQ.1) THEN
          IWORK(J)=0
        ELSE IF(IWORK(J).EQ.0) THEN
          IWORK(J)=1
          ADD=0
        END IF
* lasse addition end
        IF( ADD .EQ. 1 ) GOTO 190
C.. 2 :  CORRECT SPIN PROJECTION ?
        NUP=0
        DO J=1,NCREA
          NUP=NUP+IWORK(J)
        END DO
        DO J = NCREA+1, NCREA+NANNI
          IF(IWORK(J).EQ.0) NUP = NUP + 1
        END DO
        NDOWN = NOPEN - NUP
        MS2 = NUP - NDOWN !+ IMS2_REF
C        print*,'MS2,NUP,NDOWN',MS2,NUP,NDOWN
*. Number of alpha and beta annihilations must not be greater than
*       the number of electrons with this spinprojection
        MALPHA = 0
        MBETA = 0
        DO J = NCREA + 1, NCREA + NANNI
          IF(IWORK(J).EQ.1) THEN
            MALPHA = MALPHA + 1
          ELSE
            MBETA = MBETA + 1
          END IF
        END DO
C        print*,'MALPHA,MBETA',MALPHA,MBETA
C
* Seems to be only place needing change since only place to find ms2tot
* Though nothing seems nessecary
        IF(MS2.EQ.MS2TOT) THEN
C          print*,'MS = MS2TOT'
          IF(MALPHA.LE.NAEL) THEN
C            print*,'MALPHA LE NAEL'
            IF(MBETA.LE.NBEL) THEN
C              stop 'print'
              NCOMP = NCOMP + 1
              CALL ICOPVE(IWORK,IABCOMP(1,NCOMP),NOPEN)
            END IF
          END IF
         END IF
C          NCOMP = NCOMP + 1
C          WRITE(6,*) ' ICOMP, MALPHA, MBETA',
C     &                 ICOMP, MALPHA, MBETA
C          CALL ICOPVE(IWORK,IABCOMP(1,NCOMP),NOPEN)
C        END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ICOMP, NAEL, NBEL, MALPHA, MBETA,MS2,MS2TOT',
     &               ICOMP, NAEL, NBEL, MALPHA, MBETA,MS2,MS2TOT
      END IF
C
  200 CONTINUE
C
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Generation of spincomponents of orbitaloperator'
         WRITE(6,*)
         WRITE(6,*) ' NCREA, NANNI, MS2TOT', NCREA, NANNI,MS2TOT
         WRITE(6,*) ' Number of terms generated ', NCOMP
         WRITE(6,*)
         IF(IFLAG.NE.1) THEN
           WRITE(6,*) ' The operators (alpha = 1, beta = 0)'
           WRITE(6,*) ' ==================================='
           WRITE(6,*)
           DO 20 J=1,NCOMP
             WRITE(6,1020) J,(IABCOMP(K,J),K=1,NOPEN)
  20       CONTINUE
         END IF
 1020    FORMAT(1H0,I5,2X,30I2,/,(1H ,7X,30I2))
      END IF
*
      RETURN
      END
*
      SUBROUTINE ACT_SPOBEX_OBEX_KRCC(IFLAG,IOBEX_TP,NGAS,NOPEN,
     &           NSPCOMP_F,ISPCOMP_F,NSPCOMP_ACT,ISPCOMP_ACT,
     &           NOBPT,IACT_SPC,IAAEXC,IREF_AL,IREF_BE,
     &           IREFSPC,NAEL,NBEL,IRECOM)
*
* A set of spincomponents, ISCOMP_F  of an orbital excitation
* is given. Find the number of spinorbital excitations consistent
* with the given orbital excitation IOBEX_TP
*
* Jeppe Olsen, Summer of 99
*
* Spin combinations added Jan 2. 2000 (Sitting in the kitchen,
*                                      Jette preparing dinner )
*. Active-active excitations added March 2000, In the train
*. Test for excitation * Reference is nonvanishing, April 2001
*
* A spincomponent can be excluded of the following reasons :
* 1 : the number of alpha or beta operators in a given group is
*     larger than the number of orbitals
* 2 : When there are several creation (annihilation) operators
*     in a given orbital group, each combination of alpha and
*     beta operators is included only once.
*     For example if there are three operators belonging to a
*     given orbital space we include : a a a, a a b, a b b, b b b
*     (and not b a a, a b a, b a b, b b a )
* 4 : IF IAAEXC = 1 or 2, only selected active crea/anni of
*     active orbitals are allowed
*
*
* If IFLAG = 1, then only the number of spinorbital excitation operators is
*               returned, not the actual operators.
*
#include "implicit.inc"
#include "mxpdim.inc"
*. Local scratch
      INTEGER NCA(MXPNGAS),NCB(MXPNGAS),NAA(MXPNGAS),NAB(MXPNGAS)
*.Input
      INTEGER IOBEX_TP(2*NGAS)
      INTEGER ISPCOMP_F(NOPEN,NSPCOMP_F)
      INTEGER NOBPT(*)
*. Alpha and beta occupation of reference
      INTEGER IREF_AL(*),IREF_BE(*)
*.Output
      INTEGER ISPCOMP_ACT(4*NGAS,*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ACT_SPOEX....: Orbital excitation '
        CALL IWRTMA(IOBEX_TP,1,2*NGAS,1,2*NGAS)
      END IF
C?    WRITE(6,*) ' ACT.., alpha and beta of ref .. '
C?    CALL IWRTMA(IREF_AL,1,NGAS,1,NGAS)
C?    CALL IWRTMA(IREF_BE,1,NGAS,1,NGAS)
C?    WRITE(6,*) '  ACT_SPOBEX_OBEX, IREFSPC = ', IREFSPC
*
      NSPCOMP_ACT = 0
      IZERO = 0
      DO JSPCOMP = 1, NSPCOMP_F
*
        CALL ISETVC(NCA,IZERO,NGAS)
        CALL ISETVC(NCB,IZERO,NGAS)
        CALL ISETVC(NAA,IZERO,NGAS)
        CALL ISETVC(NAB,IZERO,NGAS)
*
        I_AM_OKAY = 1
        IOP = 1
        DO ICA = 1, 2
          IF(ICA.EQ.1) THEN
            IB = 1
          ELSE
            IB = NGAS + 1
          END IF
*. Check on creation operators
          DO IGAS = 1, NGAS
C?          WRITE(6,*) ' IGAS = ', IGAS
            LOP = NOBPT(IGAS)
C?          WRITE(6,*) ' LOP = ', LOP
            IF(IOBEX_TP(IB-1+IGAS).GT.0) THEN
              NOPFGAS = IOBEX_TP(IB-1+IGAS)
              IBA = 0
              LALPHA = 0
              LBETA = 0
              DO JOP = IOP, IOP + NOPFGAS - 1
C?              WRITE(6,*) ' JOP = ', JOP
                IF(ISPCOMP_F(JOP,JSPCOMP).EQ.1) THEN
                  LALPHA = LALPHA + 1
                  IF(LBETA.NE.0) IBA = 1
                ELSE
                  LBETA = LBETA + 1
                END IF
              END DO
              IF(ICA.EQ.1) THEN
                NCA(IGAS) = LALPHA
                NCB(IGAS) = LBETA
              ELSE
                NAA(IGAS) = LALPHA
                NAB(IGAS) = LBETA
              END IF
              IF(LALPHA.GT.LOP.OR.LBETA.GT.LOP.OR.IBA.EQ.1) THEN
                I_AM_OKAY = 0
              END IF
C             IOP = IOP + LOP
              IOP = IOP + NOPFGAS
            END IF
*           ^ End of there is operators in this space
          END DO
*         ^ End of loop over Gas spaces
        END DO
*       ^ End of loop over crea/anni
*
        I_EXCLUDE = 0
        IF(IAAEXC.EQ.1.AND.
     &     (NCA(IACT_SPC).NE.0.OR.NAB(IACT_SPC).NE.0))I_EXCLUDE = 1
        IF(IAAEXC.EQ.2.AND.
     &     (NCB(IACT_SPC).NE.0.OR.NAA(IACT_SPC).NE.0))I_EXCLUDE = 1
*. Test that annihilation on reference state is not vanishing
C        IREFSPCS = IREFSPC
C        IREFSPC = 0
        IF(IREFSPC.NE.0) THEN
*. There is a well defined reference space from which the
*. excitation are applied, check for vanishing excitations
          DO IGAS = 1, NGAS
            IF(NAA(IGAS).GT.IREF_AL(IGAS) ) I_EXCLUDE = 1
            IF(NAB(IGAS).GT.IREF_BE(IGAS) ) I_EXCLUDE = 1
          END DO
          DO IGAS = 1, NGAS
Cori            IF(NCA(IGAS)-NAA(IGAS)+IREF_AL(IGAS).GT.NOBPT(IGAS))
Cori     &      I_EXCLUDE = 1
Cori            IF(NCB(IGAS)-NAB(IGAS)+IREF_BE(IGAS).GT.NOBPT(IGAS))
Cori     &      I_EXCLUDE = 1
C removed since we will always have kramers pairs! and the limiting factor
C is the number of barred or unbarred in a space.
C And T will always be from hole to particle
            IF(NCA(IGAS)+IREF_AL(IGAS).GT.NOBPT(IGAS))
     &      I_EXCLUDE = 1
            IF(NCB(IGAS)+IREF_BE(IGAS).GT.NOBPT(IGAS))
     &      I_EXCLUDE = 1
          END DO
        ELSE
*. No well defined reference space, just check that excitation in
*. principle is feasible
          DO IGAS = 1, NGAS
            IF(NAA(IGAS).GT.MIN(NOBPT(IGAS),NAEL) ) I_EXCLUDE = 1
            IF(NAB(IGAS).GT.MIN(NOBPT(IGAS),NBEL) ) I_EXCLUDE = 1
          END DO
          DO IGAS = 1, NGAS
            IF(NCA(IGAS).GT.NOBPT(IGAS))I_EXCLUDE = 1
            IF(NCB(IGAS).GT.NOBPT(IGAS))I_EXCLUDE = 1
          END DO
        END IF
C        IREFSPC = IREFSPCS
*
* New check by lasse for the Kramers flip
*
        IF(I_AM_OKAY.EQ.1.AND.I_EXCLUDE.EQ.0) THEN
* No need to check if it allready not okay
          IMK2 = 0
          DO IGAS = 1,NGAS
            IMK2 = NCA(IGAS) - NCB(IGAS) - NAA(IGAS) + NAB(IGAS) + IMK2
          END DO
          IMK2 = IMK2/2
          IF(IRECOM.EQ.1) THEN
            CALL EVENODD(ITEST,ABS(IMK2))
            IF(ITEST.EQ.1) I_AM_OKAY = 0
          END IF
        END IF

*
        IF(I_AM_OKAY.EQ.1.AND.I_EXCLUDE.EQ.0) THEN
*. All test passed, welcome to the included
          NSPCOMP_ACT = NSPCOMP_ACT+1
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Spincombination = ', JSPCOMP
            CALL IWRTMA(ISPCOMP_F(1,JSPCOMP),1,NOPEN,1,NOPEN)
            WRITE(6,*) ' I_AM_OKAY = ', I_AM_OKAY
            WRITE(6,*) ' NAA and NAB '
            CALL IWRTMA(NAA,1,NGAS,1,NGAS)
            CALL IWRTMA(NAB,1,NGAS,1,NGAS)
          END IF
          IF(IFLAG.EQ.2) THEN
            CALL ICOPVE(NCA,ISPCOMP_ACT(1+0*NGAS,NSPCOMP_ACT),NGAS)
            CALL ICOPVE(NCB,ISPCOMP_ACT(1+1*NGAS,NSPCOMP_ACT),NGAS)
            CALL ICOPVE(NAA,ISPCOMP_ACT(1+2*NGAS,NSPCOMP_ACT),NGAS)
            CALL ICOPVE(NAB,ISPCOMP_ACT(1+3*NGAS,NSPCOMP_ACT),NGAS)
          END IF
        END IF
*       ^ End of I_AM_OKAY = 1
      END DO
*     ^ End of loop over spin combinations
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Generating allowed spincombinations '
        WRITE(6,*)
     &  ' Input per GASpace : creation and annihilation'
        CALL IWRTMA(IOBEX_TP(1)     ,1,NGAS,1,NGAS)
        CALL IWRTMA(IOBEX_TP(1+NGAS),1,NGAS,1,NGAS)
        WRITE(6,*)
     &  ' Number of included spincombinations ', NSPCOMP_ACT
        IF(IFLAG.EQ.2) THEN
          DO JSPCOMP_ACT = 1, NSPCOMP_ACT
            WRITE(6,*)
            WRITE(6,*) ' Included spinorbitalexcitation ', JSPCOMP_ACT
            WRITE(6,'(A,16I4)')
     &      ' Creation of alpha     :',
     &      (ISPCOMP_ACT(I+0*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)')
     &      ' Creation of beta      :',
     &      (ISPCOMP_ACT(I+1*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)')
     &      ' Annihilation of alpha :',
     &      (ISPCOMP_ACT(I+2*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)')
     &      ' Annihilation of beta  :',
     &      (ISPCOMP_ACT(I+3*NGAS,JSPCOMP_ACT),I=1,NGAS)
          END DO
        END IF
      END IF
*     ^ End if NTEST .ge. 100
*
      RETURN
      END
*
      SUBROUTINE SYM_INF_DIRAC
*************************************************************
*                                                           *
* Routine to transfer info from dirac common blocks to KRCC *
*                                                           *
*************************************************************
#include "implicit.h"
C Dirac common blocks
#include "pgroup.h"
C KRCC common blocks
#include "mxpdim.inc"
#include "lucinp.inc" 
C
C     Double group
      IF (GROUP.eq.'C1 ') THEN
         PNTGRP = 5
         NIRREP = 2
         IRIQ = 2
      ELSE IF (GROUP.eq.'Ci ') THEN
         PNTGRP = 6
         NIRREP = 4
         IRIQ = 2
      ELSE IF (GROUP.eq.'Cs ') THEN
         PNTGRP = 7
         NIRREP = 4
         IRIQ = 2
      ELSE IF (GROUP.eq.'C2 ') THEN
         PNTGRP = 8
         NIRREP = 4
         IRIQ = 2
      ELSE IF (GROUP.eq.'C2h') THEN
         PNTGRP = 9
         NIRREP = 8
         IRIQ = 2
      ELSE IF (GROUP.eq.'C2v') THEN
         PNTGRP = 8
         NIRREP = 4
         IRIQ = 1
      ELSE IF (GROUP.eq.'D2 ') THEN
         PNTGRP = 8
         NIRREP = 4
         IRIQ = 1
      ELSE IF (GROUP.eq.'D2h') THEN
         PNTGRP = 9
         NIRREP = 8
         IRIQ = 1
      ELSE
         CALL QUIT('Illegal double group.')
      END IF
      NSMOB = NIRREP
      RETURN
      END
*
      SUBROUTINE FIND_MUB(IEXCITE,IKRFLIP,IVEC)
*
* Find the possible MUB values for a given Kramers flip
*
#include "implicit.inc"
*
      INTEGER IVEC(0:IEXCITE,0:IKRFLIP,-IEXCITE:IEXCITE)
*
      NTEST = 00
*
      DO I = 0,IEXCITE
        DO J = 0,MIN(I,IKRFLIP)
          IMUBMAX = I-J
          DO K = -IMUBMAX,IMUBMAX,2
            IVEC(I,J,K) = K
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
      WRITE(6,*) ' All possible Mub values for given N and Mk'
        DO I = 0,IEXCITE
          DO J = 0,MIN(I,IKRFLIP)
            IMUBMAX = I-J
            WRITE(6,*) ' N ',I,' Mk ',J
            WRITE(6,*) ' Mub Range '
            DO K = -IMUBMAX,IMUBMAX,2
              WRITE(6,*) K
            END DO
          END DO
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE SORT_OP_3_INDEX(IOP,IDIMOP,NGAS,
     &                           IEXCITE,IKRFLIP,
     &                           IVECS,IVECF,
     &                           IVECH,IHELP_OP,IHELP_IDX)
*
* Will sort the operator according to N,Mk,Mub and return the sorted
* operator along with position of and operator with given  N,Mk,Mub in
* new sorting
*
#include "implicit.inc"
*
      INTEGER IVECS(0:IEXCITE,-IKRFLIP:IKRFLIP,-IEXCITE:IEXCITE)
      INTEGER IVECF(0:IEXCITE,-IKRFLIP:IKRFLIP,-IEXCITE:IEXCITE)
      INTEGER IVECH(0:IEXCITE,-IKRFLIP:IKRFLIP,-IEXCITE:IEXCITE)
      INTEGER IOP(4*NGAS,IDIMOP)
      INTEGER IHELP_OP(4*NGAS,IDIMOP),IHELP_IDX(IDIMOP)
*
      NTEST = 00
*
* Find the number of operators with given N,Mk,Mub
*
C     print*,'NGAS,IDIMOP',NGAS,IDIMOP
      DO I =1,IDIMOP
        CALL ANALYZE_OP(IOP(1,I),NGAS,N,MK,MUB)
        IVECS(N,MK,MUB) = IVECS(N,MK,MUB) + 1
C       print*,'IVECS(N,MK,MUB),N,MK,MUB',IVECS(N,MK,MUB),N,MK,MUB
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of operators of given N,Mk,Mub '
        DO I = 0,IEXCITE
          DO J = -IKRFLIP,IKRFLIP
            DO K = -IEXCITE,IEXCITE
              IF(IVECS(I,J,K).GT.0) THEN
                WRITE(6,*) IVECS(I,J,K),' Operators with ',
     &                     'N ',I,' Mk',J,' Mub ',K
              END IF
            END DO
          END DO
        END DO
      END IF
*
* Want to store operators as N,-Mk to Mk,-Mub to Mub at the moment!!!
* Therefore find the range where an operator can be stored
*
      ISUM = 0
*
      DO I =0,IEXCITE
        DO J = -IKRFLIP,IKRFLIP
          DO K = -IEXCITE,IEXCITE
            ISUM = ISUM + IVECS(I,J,K)
            IVECF(I,J,K) = ISUM 
            IVECS(I,J,K) = ISUM - IVECS(I,J,K) + 1
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        IDIFF = 0
        WRITE(6,*) ' Start and finish for Operators '
        DO I =0,IEXCITE
          DO J = -IKRFLIP,IKRFLIP
            DO K = -IEXCITE,IEXCITE
              IF(IVECF(I,J,K)-IVECS(I,J,K).GE.0) THEN
                WRITE(6,*) ' Operators with ',
     &                     'N ',I,' Mk',J,' Mub ',K,
     &                     ' begins ',IVECS(I,J,K),
     &                     ' ends ',IVECF(I,J,K)
                IDIFF = IDIFF + 1
              END IF
            END DO
          END DO
        END DO
        WRITE(6,*) ' Number of different operators ',IDIFF
      END IF
*
* Now find the new place for the operators and copy to temp array
*
      DO I =1,IDIMOP
        CALL ANALYZE_OP(IOP(1,I),NGAS,N,MK,MUB)
        IHELP_IDX(I) = IVECS(N,MK,MUB) + IVECH(N,MK,MUB)
        CALL ICOPVE(IOP(1,I),IHELP_OP(1,IHELP_IDX(I)),4*NGAS)
        IVECH(N,MK,MUB) = IVECH(N,MK,MUB) + 1
      END DO
*
      IF(NTEST.GE.100) THEN
        DO I =1,IDIMOP
          WRITE(6,*) ' Original order ',I,' New order ',IHELP_IDX(I)
        END DO
      END IF
*
* Copy operators from temp array back again
*
      DO I =1,IDIMOP
        CALL ICOPVE(IHELP_OP(1,I),IOP(1,I),4*NGAS)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' New order of operators '
        CALL WRT_SPOX_TP_CC_KRCC(IOP,IDIMOP)
      END IF
*
      RETURN
      END
*
      SUBROUTINE ANALYZE_OP(IOP,NGAS,N,MK,MUB)
*
* Will give N,MK,MUB of any operator
*
#include "implicit.inc"
*
      INTEGER IOP(4*NGAS)
*
      NTEST = 00
*
      NCA = 0
      NCB = 0
      NAA = 0
      NAB = 0
*
* Find number of creator and annihilators of alpha and beta
*
      DO IGAS =1,NGAS
* Creator alpha
        NCA = NCA + IOP(IGAS+0*NGAS) 
* Creator Beta
        NCB = NCB + IOP(IGAS+1*NGAS) 
* Annihilator Alpha
        NAA = NAA + IOP(IGAS+2*NGAS)
* Annihilator Beta
        NAB = NAB + IOP(IGAS+3*NGAS)
      END DO
*
* Particle rank
*
      N = (NCA + NCB + NAA + NAB)/2
*
* Kramers flip
*
      MK = (NCA - NCB - NAA + NAB)/2
*
* Mub
*
      MUB = (NCA - NCB + NAA - NAB)/2
*
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Information from ANALYZE_OP '
        CALL WRT_SPOX_TP_CC_KRCC(IOP,1)
        WRITE(6,*) ' Analysis results '
        WRITE(6,*) ' Particle rank = ',N
        WRITE(6,*) ' Kramers flip  = ',MK
        WRITE(6,*) ' MUB           = ',MUB
        IF(NTEST.EQ.1000) THEN
          WRITE(6,*) ' Number of GAS = ',NGAS
          WRITE(6,*) ' Number of CA  = ',NCA
          WRITE(6,*) ' Number of CB  = ',NCB
          WRITE(6,*) ' Number of AA  = ',NAA
          WRITE(6,*) ' Number of AB  = ',NAB
        END IF
      END IF
*
      RETURN
      END
