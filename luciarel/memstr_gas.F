      SUBROUTINE MEMSTR_GAS
*
*
* Construct pointers for saving information about strings and
* their mappings
*
* GAS version
*
*========
* Input :
*========
* Number and groups of strings defined by /GASSTR/
* Symmetry information stored in         /CSM/
* String information stored in           /STINF/
*=========
* Output
*=========
* Pointers stored in common block /STRBAS/
*
* Jeppe Olsen , Winter of 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "strbas.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "stinf.inc"

*
*. Start of string information
      CALL MEMMAN(KSTINF,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
      NTEST = 0
      IF(NTEST.NE.0)
     &WRITE(6,*) ' First word with string information',KSTINF
*
*.  Offsets for occupation and reorder array of strings
*
      DO IGRP = 1, NGRP
        NSTRIN = NSTFGP(IGRP)
        LSTRIN = NSTRIN*NELFGP(IGRP)
        CALL MEMMAN(KOCSTR(IGRP),LSTRIN,'ADDS  ',1,'OCSTR ')
        CALL MEMMAN(KSTREO(IGRP),NSTRIN,'ADDS  ',1,'STREO ')
        IF(PNTGRP.GE.5) THEN
          CALL MEMMAN(KOCSTR2(IGRP),LSTRIN,'ADDS  ',1,'OCSTR2')
          CALL MEMMAN(KSTREO2(IGRP),NSTRIN,'ADDS  ',1,'STREO2')
        END IF
*. For spin down strings
      END DO
*
*. Number of strings per symmetry and offset for strings of given sym
*. for groups
*
      CALL MEMMAN(KNSTSGP(1),NSMST*NGRP,'ADDS  ',1,'NSTSGP')
      CALL MEMMAN(KISTSGP(1),NSMST*NGRP,'ADDS  ',1,'ISTSGP')
      IF(PNTGRP.GE.5) THEN
        CALL MEMMAN(KNSTSGP2(1),NSMST*NGRP,'ADDS  ',1,'NSTSGP')
        CALL MEMMAN(KISTSGP2(1),NSMST*NGRP,'ADDS  ',1,'ISTSGP')
      END IF


*
*. Number of strings per symmetry and offset for strings of given sym
*. for types
*
      DO  ITP  = 1, NSTTP
        CALL MEMMAN(KNSTSO(ITP),NSPGPFTP(ITP)*NSMST,'ADDS  ',1,'NSTSO ')
        CALL MEMMAN(KISTSO(ITP),NSPGPFTP(ITP)*NSMST,'ADDS  ',1,'ISTSO ')
      END DO
      IF(PNTGRP.GE.5) THEN
       DO ITP = 1, NSTTP
         LEN_T = NSPGPFTP(ITP)*NSMST
         CALL MEMMAN(KNSTSO2(ITP),LEN_T,'ADDS  ',1,'NSTSO2')
         CALL MEMMAN(KISTSO2(ITP),LEN_T,'ADDS  ',1,'ISTSO2')
       END DO
      END IF
*
**. Lexical adressing of arrays : use array indeces for complete active space
*
*. Not in use so
      DO  IGRP = 1, NGRP
        CALL MEMMAN(KZ(IGRP),NACOB*NELFGP(IGRP),'ADDS  ',1,'Zmat  ')
      END DO
*
*. Mappings between different groups
*
      DO  IGRP = 1, NGRP
        IEL = NELFGP(IGRP)
        IGAS = IGSFGP(IGRP)
        IORB = NOBPT(IGAS)
        ISTRIN = NSTFGP(IGRP)
*. Use full orbital notation
*. Note Only creation maps are currently installed,
*. Generalize when/if annihilation maps are included
        IF(ISTAC(IGRP,2).NE.0) THEN
          LENGTH = IORB*ISTRIN
          CALL MEMMAN(KSTSTM(IGRP,1),LENGTH,'ADDS  ',1,'ORBMAP')
          CALL MEMMAN(KSTSTM(IGRP,2),LENGTH,'ADDS  ',1,'STRMAP')
          IF(PNTGRP.GE.5) THEN
            CALL MEMMAN(KSTSTM2(IGRP,1),LENGTH,'ADDS  ',1,'ORBMAP')
            CALL MEMMAN(KSTSTM2(IGRP,2),LENGTH,'ADDS  ',1,'STRMAP')
          END IF
        ELSE
          KSTSTM(IGRP,1) = -1
          KSTSTM(IGRP,2) = -1
          IF(PNTGRP.GE.5) THEN
            KSTSTM2(IGRP,1) = -1
            KSTSTM2(IGRP,2) = -1
          END IF
        END IF
      END DO
*
*. Symmetry of conjugated orbitals and orbital excitations
*
*     KCOBSM,KNIFSJ,KIFSJ,KIFSJO
      CALL MEMMAN(KCOBSM,NACOB,'ADDS  ',1,'Cobsm ')
      CALL MEMMAN(KNIFSJ,NACOB*NSMSX,'ADDS  ',1,'Nifsj ')
      CALL MEMMAN(KIFSJ,NACOB**2 ,'ADDS  ',1,'Ifsj  ')
      CALL MEMMAN(KIFSJO,NACOB*NSMSX,'ADDS  ',1,'Ifsjo ')
*
*. Symmetry of excitation connecting  strings of given symmetry
*
      CALL MEMMAN(KSTSTX,NSMST*NSMST,'ADDS  ',1,'Ststx ')
*
*. Occupation classes
*
      CALL MEMMAN(KIOCLS,NMXOCCLS*NGAS,'ADDS  ',1,'IOCLS ')
*. Annihilation/Creation map of supergroup types
      CALL MEMMAN(KSPGPAN,NTSPGP*NGAS,'ADDS  ',1,'SPGPAN')
      CALL MEMMAN(KSPGPCR,NTSPGP*NGAS,'ADDS  ',1,'SPGPCR')
*
*. Last word of string information
*
      CALL MEMMAN(KSTINE,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
      IF(NTEST.NE.0)
     &WRITE(6,*) ' Last word with string information',KSTINE-1
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' R*8 words used for storing Strinformation ',
     &               KSTINE - KSTINF
      END IF
*
      WRITE(6,*) ' Memory map '
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'PRINT',IDUM,'IDUMMY')
      RETURN
      END
