      SUBROUTINE MEMSTR_GAS_REL
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
      call memmar(KSTINF,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
      NTEST = 0
      IF(NTEST.NE.0)
     &WRITE(6,*) ' First word with string information',KSTINF
*
*.  Offsets for occupation and reorder array of strings
*
      DO IGRP = 1, NGRP
        NSTRIN = NSTFGP(IGRP)
        LSTRIN = NSTRIN*NELFGP(IGRP)
        call memmar(KOCSTR(IGRP),LSTRIN,'ADDS  ',1,'OCSTR ')
        call memmar(KSTREO(IGRP),NSTRIN,'ADDS  ',1,'STREO ')
        IF(PNTGRP.GE.5) THEN
        call memmar(KOCSTR2(IGRP),LSTRIN,'ADDS  ',1,'OCSTR2')
        call memmar(KSTREO2(IGRP),NSTRIN,'ADDS  ',1,'STREO2')
        END IF
*. For spin down strings
      END DO
*
*. Number of strings per symmetry and offset for strings of given sym
*. for groups
*
      call memmar(KNSTSGP(1),NSMST*NGRP,'ADDS  ',1,'NSTSGP')
      call memmar(KISTSGP(1),NSMST*NGRP,'ADDS  ',1,'ISTSGP')
      IF(PNTGRP.GE.5) THEN
        call memmar(KNSTSGP2(1),NSMST*NGRP,'ADDS  ',1,'NSTSGP')
        call memmar(KISTSGP2(1),NSMST*NGRP,'ADDS  ',1,'ISTSGP')
      END IF


*
*. Number of strings per symmetry and offset for strings of given sym
*. for types
*
      DO  ITP  = 1, NSTTP
        call memmar(KNSTSO(ITP),NSPGPFTP(ITP)*NSMST,'ADDS  ',1,'NSTSO ')
        call memmar(KISTSO(ITP),NSPGPFTP(ITP)*NSMST,'ADDS  ',1,'ISTSO ')
      END DO
      IF(PNTGRP.GE.5) THEN
       DO ITP = 1, NSTTP
         LEN_T = NSPGPFTP(ITP)*NSMST
         call memmar(KNSTSO2(ITP),LEN_T,'ADDS  ',1,'NSTSO2')
         call memmar(KISTSO2(ITP),LEN_T,'ADDS  ',1,'ISTSO2')
       END DO
      END IF
*
**. Lexical adressing of arrays : use array indeces for complete active space
*
*. Not in use so
      DO  IGRP = 1, NGRP
        call memmar(KZ(IGRP),NACOB*NELFGP(IGRP),'ADDS  ',1,'Zmat  ')
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
          call memmar(KSTSTM(IGRP,1),LENGTH,'ADDS  ',1,'ORBMAP')
          call memmar(KSTSTM(IGRP,2),LENGTH,'ADDS  ',1,'STRMAP')
          IF(PNTGRP.GE.5) THEN
            call memmar(KSTSTM2(IGRP,1),LENGTH,'ADDS  ',1,'ORBMAP')
            call memmar(KSTSTM2(IGRP,2),LENGTH,'ADDS  ',1,'STRMAP')
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
      call memmar(KCOBSM,NACOB,'ADDS  ',1,'Cobsm ')
      call memmar(KNIFSJ,NACOB*NSMSX,'ADDS  ',1,'Nifsj ')
      call memmar(KIFSJ,NACOB**2 ,'ADDS  ',1,'Ifsj  ')
      call memmar(KIFSJO,NACOB*NSMSX,'ADDS  ',1,'Ifsjo ')
*
*. Symmetry of excitation connecting  strings of given symmetry
*
      call memmar(KSTSTX,NSMST*NSMST,'ADDS  ',1,'Ststx ')
*
*. Occupation classes
*
      call memmar(KIOCLS,NMXOCCLS*NGAS,'ADDS  ',1,'IOCLS ')
*. Annihilation/Creation map of supergroup types
      call memmar(KSPGPAN,NTSPGP*NGAS,'ADDS  ',1,'SPGPAN')
      call memmar(KSPGPCR,NTSPGP*NGAS,'ADDS  ',1,'SPGPCR')
*
*. Last word of string information
*
      call memmar(KSTINE,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
      IF(NTEST.NE.0)
     &WRITE(6,*) ' Last word with string information',KSTINE-1
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' R*8 words used for storing Strinformation ',
     &               KSTINE - KSTINF
      END IF
*
      WRITE(6,*) ' Memory map '
      IDUM = 0
      call memmar(IDUM,IDUM,'PRINT',IDUM,'IDUMMY')
      RETURN
      END
