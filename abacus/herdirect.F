C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dirac, a relativistic ab initio electronic structure program, 
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : herdirect.F
C SOURCE  : $Source$
C
C
#ifdef UNDEF
/* Comdeck herdirect_log */
941206-hjaaj
HERPAR: new routine for defining parallel hermit calc. of Fock matr.
HERFCK: call PARDRV if PARALL true
        changed STATIS to KPRINT.gt.0 (STATIS missing in newest
        aba.cdk from oslo)
941205-from oslo:
HERFCK: also IPRINT = IPRFCK for PARALL true (was = 0)
941103-hjaaj:
implemented ISYMDM(NDMAT), IFCTYP(NDMAT)
HERDI1,HERDI2: excluded; not used, are also incomplete
ABARUN: moved to here from abadrv.u (needed for HERMIT)
931026-hjaaj: sether now calls READIN; + consistency check
HERFCK: new parameter (print level in CBITWO)
931025-hjaaj: sirdirn04.u (i.e. HERFCK,SETHER) added to herdirect.u
GETFCK now calls HERFCK
AOLABE renamed to AOLAB4, generalized, and moved to GP.
GETFCK is now only internal test routine; HERFCK is called from the
outside (SIRIUS and RESPONS).
921201-hjaaj:
   new Sirius input routine for Hermit parameters ? an extended sether?
   also direct calculation of solvent integrals ?
#endif
C  /* Deck tstdir */
      SUBROUTINE TSTDIR(FMAT,DMAT,NDMAT,WORK,LWORK,IPRINT)
C
C     Written by Henrik Koch and Trygve Helgaker 27-November-1991.
C
C     PURPOSE : Driver routine for the test calculation of the
C               two-electron part of the fock matrices.
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      PARAMETER (MAXDIS=500)
      DIMENSION DMAT(*), FMAT(*), WORK(LWORK)
#include <inforb.h>
C
      CALL TIMER('START ',TIMSTR,TIMEND)
      KGMAT = 1
      KINDX = KGMAT + NBAST*NBAST*MAXDIS
      KEND  = KINDX + (MAXDIS + 1)/IRAT
      LFREE = LWORK - KEND
      IF (KEND .GT. LWORK) CALL STOPIT('TSTDIR',' ',KEND,LWORK)
C
      NUMDIS = -1
    1 CONTINUE
C
C     Get distributions
C
      CALL GETDIS(WORK(KGMAT),WORK(KINDX),NUMDIS,MAXDIS,WORK(KEND),
     &            LFREE)
      IF (NUMDIS .EQ. -1) GOTO 2
C
C     Add contributions to Fock matrix
C
      CALL DISFCK(FMAT,DMAT,NDMAT,WORK(KGMAT),WORK(KINDX),NUMDIS)
C
      GOTO 1
C
    2 CONTINUE
C
C     Fock matrix has now been constructed
C
C------------------------------------------------------
C     Write out densities and associated fock matrices.
C------------------------------------------------------
C
      IF (IPRINT.GT.2) THEN
         CALL HEADER('Density and Fock matrices from '//
     &               'distributions in TSTDIR',-1)
         DO 300 I = 1, NDMAT
            ISTR = NBAST*NBAST*(I - 1) + 1
            WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
            CALL OUTPUT(DMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  300    CONTINUE
      END IF
      CALL TIMER('GETDIS',TIMSTR,TIMEND)
      RETURN
      END
C  /* Deck disfck */
      SUBROUTINE DISFCK(FMAT,DMAT,NDMAT,GMAT,INDEX,NUMDIS)
C
C     Written by Henrik Koch and Trygve Helgaker 27-November-1991.
C
C
C     PURPOSE : Calculate fock matrices using distributions.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ONE=1.0D00, HALF=0.5D00, FOURTH=0.25D00)
      INTEGER   P,Q,R,S
      DIMENSION FMAT(NBAST,NBAST,NDMAT), DMAT(NBAST,NBAST,NDMAT)
      DIMENSION GMAT(NBAST,NBAST,NUMDIS), INDEX(NUMDIS)
#include <inforb.h>
C
      DO 100 I = 1,NUMDIS
C
         CALL UNPKIJ(INDEX(I),R,S)
C
         IF (R .EQ. S) THEN
            FACCOU = HALF
            FACEXC = FOURTH
         ELSE
            FACCOU = ONE
            FACEXC = HALF
         END IF
C
         DO 110 J = 1,NDMAT
C
            VALUE = FACCOU*DDOT(NBAST*NBAST,DMAT(1,1,J),1,GMAT(1,1,I),1)
            FMAT(R,S,J) = FMAT(R,S,J) + VALUE
            FMAT(S,R,J) = FMAT(S,R,J) + VALUE
C
            DO 120 Q = 1,NBAST
C
               FMAT(R,Q,J) = FMAT(R,Q,J)
     &                  - FACEXC*DDOT(NBAST,DMAT(1,S,J),1,GMAT(1,Q,I),1)
               FMAT(S,Q,J) = FMAT(S,Q,J)
     &                  - FACEXC*DDOT(NBAST,DMAT(1,R,J),1,GMAT(1,Q,I),1)
C
  120       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck getdis */
      SUBROUTINE GETDIS(GMAT,INDXAB,NUMDIS,MAXDIS,WORK,LWORK)
C
C     Written by Henrik Koch and Trygve Helgaker 26-November-1991.
C
C     PURPOSE : Driver routine for the calculation of the two-electron
C               distributions (**|cd).
C               The distributions are stored as full squares and without
C               symmetry reduction.
C
C                      *CALL CBITWO   from TWOINP
C                      *CALL HRUNIT   from HERMIT
C                      *CALL INFORB   from sirius
C                      *CALL GAMCOM   from READIN
C                      *CALL CCOM     from READIN
C                      *CALL PINCOM   from READIN
C                      *CALL PRIMIT   from READIN
C                      *CALL XYZPOW   from READIN
C                      *CALL SHELLS   from READIN
C                      *CALL NUCLEI   from READIN
C                      *CALL SYMMET   from READIN
C                      *CALL DORPS    from abatro.u
C                       ( only DOREPS and DOCOOR used)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <dummy.h>
      LOGICAL ABA
      DIMENSION GMAT(NBAST,NBAST,MAXDIS), INDXAB(MAXDIS)
      DIMENSION WORK(LWORK)
#include <inforb.h>
#include <cbitwo.h>
#include <dorps.h>
#include <ccom.h>
C
      CALL ABARUN(ABA)
C
      I2TYP = 0
C
C-------------------------------------------------------------
C     Setup information for the two-electron integralroutines.
C-------------------------------------------------------------
C
      IF (.NOT.ABA) THEN
         NODV = .FALSE.
         NOPV = .FALSE.
         IF (NASHT .EQ. 0) NODV = .TRUE.
         IF (NASHT .LT. 2) NOPV = .TRUE.
         IPRALL =  0
         DO 100 I = 0,7
            DOREPS(I) = .TRUE.
  100    CONTINUE
         DO 110 I = 1,MXCENT
            DOCOOR(1,I) = .TRUE.
            DOCOOR(2,I) = .TRUE.
            DOCOOR(3,I) = .TRUE.
  110    CONTINUE
         IPRNTA = 0
         IPRNTB = 0
         IPRNTC = 0
         IPRNTD = 0
         RETUR  = .FALSE.
         NOCONT = .FALSE.
         TKTIME = .FALSE.
      END IF
C
      MAXDIF = 0
      ITYPE  = 4
      CALL PAOVEC(WORK,LWORK,0,IPRALL)
C
C----------------------------
C     Calculate distributions
C----------------------------
C
      CALL TIMER('START ',TIMSTR,TIMEND)
      CALL TWOINT(WORK,LWORK,DUMMY,DUMMY,NDMAT,IDUMMY,IDUMMY,GMAT,
     &            INDXAB,NUMDIS,MAXDIS,ITYPE,MAXDIF,0,NODV,NOPV,NOCONT,
     &            TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,
     &            IDUMMY,I2TYP,IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,IDUMMY,
     &            .FALSE.)
      CALL TIMER('TWOINT',TIMSTR,TIMEND)
C
C--------------------------------------
C     Write out integrals distributions
C--------------------------------------
C
      IF (IPRINT.GT.3) THEN
         CALL HEADER('Distribution matrices in GETDIS',-1)
         DO 300 I = 1, NUMDIS
            CALL UNPKIJ(INDXAB(I),IA,IB)
            WRITE(LUPRI,'(//,1X,A,I4)')' Distribution matrix No.',I
            WRITE(LUPRI,'(1X,A,2I4)')  ' Orbital indices:       ',IA,IB
            CALL OUTPUT(GMAT(1,1,I),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  300    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck getfck */
#if defined (VAR_GETFCK)
C
C     Obsolete
C
      SUBROUTINE GETFCK(FMAT,DMAT,NDMAT,WORK,LWORK,DIRECT,TRIPLE)
C
C     Written by Henrik Koch and Trygve Helgaker 22-November-1991.
C
C     PURPOSE : Driver routine for the calculation of the two-electron
C               part of the fock matrices.
C               We assume the densities and fock matrices are full
C               squares and without symmtry reduction .
C               The subroutine offers the possibility of calculating
C               the fock matrices directly or by reading the integrals
C               from disk.
C
C               DIRECT = .TRUE.   Direct calculation using abacus.
C               DIRECT = .FALSE.  Read integrals from AOTWOINT.
C
C               TRIPLE = .TRUE.   Triplet fock-matrix.
C               TRIPLE = .FALSE.  Singlet fock-matrix.
C
C               Common blocks needed to be initialized before this
C               subroutine is called:
C
C                      *CALL CBITWO   from TWOINP
C                      *CALL HRUNIT   from HERMIT
C                      *CALL INFORB   from sirius
C                      *CALL GAMCOM   from READIN
C                      *CALL CCOM     from READIN
C                      *CALL PINCOM   from READIN
C                      *CALL PRIMIT   from READIN
C                      *CALL XYZPOW   from READIN
C                      *CALL SHELLS   from READIN
C                      *CALL NUCLEI   from READIN
C                      *CALL SYMMET   from READIN
C                      *CALL DORPS    from abatro.u
C                       ( only DOREPS and DOCOOR used)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <dummy.h>
      LOGICAL   DIRECT, TRIPLE
      DIMENSION FMAT(*), DMAT(*)
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  INFORB: NBAST,?
C  CBITWO: IPRINT,?
C
#include <inforb.h>
#include <cbitwo.h>
#include <dorps.h>
#include <hrunit.h>
#include <symmet.h>
C
      IF (DIRECT) THEN
         JPRINT = 2
         CALL HERFCK(FMAT,DMAT,NDMAT,WORK,LWORK,TRIPLE,JPRINT)
      ELSE
         IF (TRIPLE) THEN
            CALL QUIT('GETFCK error: triplet not implemented yet!')
         END IF
C
         LUINTA = 11
         OPEN (LUINTA,STATUS='OLD',FORM='UNFORMATTED',
     *         FILE='AOTWOINT')
C
         KINT  = 1
         KIINT = KINT  + 600
         KIND4 = KIINT + (601+1)/IRAT
         KEND  = KIND4 + (4*600+1)/IRAT
C
         LENINT = IRAT*600 + 600 + 1
C
         REWIND LUINTA
         CALL MOLLAB('BASTWOEL',LUINTA,6)
C
    1    CONTINUE
C
         CALL READI(LUINTA,LENINT,WORK(KINT))
         CALL AOLAB4(WORK(KIINT),600,WORK(KIND4),NUMINT)
         IF (NUMINT .EQ. 0)  GOTO 1
         IF (NUMINT .EQ. -1) GOTO 2
C
         CALL FCKDIR(FMAT,DMAT,NDMAT,WORK(KINT),WORK(KIND4),
     &               NUMINT,IPRINT)
C
         GOTO 1
C
    2    CONTINUE
C
      ENDIF
C
C------------------------------------------------------
C     Write out densities and associated fock matrices.
C------------------------------------------------------
C
C     IF (IPRINT.GT.2) THEN
         CALL HEADER('Density and Fock matrices in GETFCK',-1)
         DO 300 I = 1, NDMAT
            ISTR = NBAST*NBAST*(I - 1) + 1
            WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
            CALL OUTPUT(DMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  300    CONTINUE
C     END IF
C
      RETURN
      END
#endif
C  /* Deck fcktes */
#if !defined (VAR_GETFCK)
      SUBROUTINE FCKTES(WORK,LWORK,MAXDIF,NODV,NOPV,NOCONT,TKTIME,
     &                  IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR)
C     April 1995: FCKTES below must be changed to call new Fock matrix
C     routines instead of GETFCK.
      CALL QUIT('FCKTES not implemented in this version!')
      END
#else
      SUBROUTINE FCKTES(WORK,LWORK,MAXDIF,NODV,NOPV,NOCONT,TKTIME,
     &                  IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0)
      LOGICAL NODV, NOPV, NOCONT, TKTIME, RETUR, DIRECT
      DIMENSION WORK(LWORK)
#include <inforb.h>
C
#if !defined (PRG_DIRAC)
      IF (NODV) THEN
         NDMAT = 1
      ELSE
         NDMAT = 2
      END IF
      KFMAT = 1
      KDMAT = KFMAT + N2BASX*NDMAT
      KCMO  = KDMAT + N2BASX*NDMAT
      KDV   = KCMO  + NCMOT
      KDTSO = KDV   + NNASHX
      KDASO = KDTSO + NNBAST
      KLAST = KDASO + NNBAST
      IF (KLAST .GT. LWORK) CALL STOPIT('FCKTES',' ',KLAST,LWORK)
      CALL ONEDSF(WORK(KCMO),WORK(KDV),WORK(KDMAT),
     &            WORK(KDMAT+N2BASX),IPRINT,NODC,NODV)
      CALL DAXPY(N2BASX,-D1,WORK(KDMAT+N2BASX),1,WORK(KDMAT),1)
      CALL RDONEL('ONEHAMIL',.TRUE.,WORK(KDTSO),NNBAST)
      CALL SQDENS(WORK(KDTSO),WORK(KFMAT),IPRINT)
      IF (.NOT.NODV) CALL DZERO(WORK(KFMAT+N2BASX),N2BASX)
C
      DIRECT = .TRUE.
      KLAST = KCMO
      LWRK  = LWORK - KLAST + 1
      CALL TIMER('START ',TIMSTR,TIMEND)
      CALL GETFCK(WORK(KFMAT),WORK(KDMAT),NDMAT,WORK(KLAST),LWRK,DIRECT,
     *          .FALSE.)
      CALL TIMER('GETFCK',TIMSTR,TIMEND)
      if (.true.) return
C
C     Test distributions
C
      KFMAT = 1
      KDMAT = KFMAT + N2BASX*NDMAT
      KCMO  = KDMAT + N2BASX*NDMAT
      KDV   = KCMO  + NCMOT
      KDTSO = KDV   + NNASHX
      KDASO = KDTSO + NNBAST
      KLAST = KDASO + NNBAST
      IF (KLAST .GT. LWORK) CALL STOPIT('FCKTES',' ',KLAST,LWORK)
      CALL ONEDSF(WORK(KCMO),WORK(KDV),WORK(KDMAT),WORK(KDMAT+N2BASX),
     &            IPRINT,NODC,NODV)
      CALL DAXPY(N2BASX,-D1,WORK(KDMAT+N2BASX),1,WORK(KDMAT),1)
      CALL RDONEL('ONEHAMIL',.TRUE.,WORK(KDTSO),NNBAST)
      CALL SQDENS(WORK(KDTSO),WORK(KFMAT),IPRINT)
      IF (.NOT.NODV) CALL DZERO(WORK(KFMAT+N2BASX),N2BASX)
C
      KLAST = KCMO
      LWRK  = LWORK - KLAST + 1
      CALL TSTDIR(WORK(KFMAT),WORK(KDMAT),NDMAT,WORK(KLAST),LWRK,IPRINT)
#endif
      RETURN
      END
#endif
C  /* Deck sqdens */
      SUBROUTINE SQDENS(DTSO,DMAT,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      DIMENSION DTSO(*), DMAT(NBAST,NBAST)
#include <inforb.h>
#include <symmet.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from SQDENS','*',103)
      CALL DZERO(DMAT,NBAST*NBAST)
      ISOFF = 0
      DO 100 IREP = 1, NSYM
         DO 200 I = 1, NBAS(IREP)
            INDA = IBAS(IREP) + I
            DO 300 J = 1, I
               INDB = IBAS(IREP) + J
               ISOFF = ISOFF + 1
               DMAT(INDA,INDB) = DMAT(INDA,INDB) + DTSO(ISOFF)
               DMAT(INDB,INDA) = DMAT(INDA,INDB)
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Squared matrix in SQDENS',-1)
         CALL OUTPUT(DMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck sether */
      SUBROUTINE SETHER(JPRINT,NEWGEO,RELCAL,WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  CBIREA: LUMLCL
C  SYMMET: MAXREP,NAOS(8))
C  INFORB: NSYM,NBAS(8)
C
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <cbirea.h>
#include <symmet.h>
#include <inforb.h>
C
      LOGICAL SFIRST, RELCAL, TSTINP, HERMIT, NEWGEO
      SAVE    SFIRST
      DATA    SFIRST /.TRUE./
C
C     Initialize /CBIREA/ and call readin
C
      IF (SFIRST .OR. NEWGEO) THEN
         TSTINP = .FALSE.
         CALL REAINI(JPRINT,RELCAL,TSTINP)
         LUMLCL = 56
C
         IPREAD = JPRINT
         HERMIT = .FALSE.
         CALL READIN(WORK,LWORK,HERMIT)
C
C        consistency check
C
         IF(.NOT.RELCAL) THEN
           IF (MAXREP+1 .NE. NSYM) THEN
            WRITE (LUPRI,'(2(/A,I5))')
     &      ' SETHER fatal error: NSYM  from SIRIUS is',NSYM,
     &      '                  MAXREP+1 from READIN is',MAXREP+1
            CALL QUIT('SETHER error: NSYM .ne. MAXREP+1')
           END IF
           NERR = 0
           DO 100 ISYM = 1,NSYM
              IF (NAOS(ISYM) .NE. NBAS(ISYM)) NERR = NERR + 1
  100      CONTINUE
           IF (NERR .GT. 0) THEN
              WRITE (LUPRI,'(/A/A/)')
     &        ' SETHER fatal error: NBAS(:) from LUONEL',
     &        '   is different from NAOS(:) from READIN.'
              WRITE (LUPRI,'(A,8I5)') ' NBAS(:) =',(NBAS(I),I=1,NSYM)
              WRITE (LUPRI,'(A,8I5)') ' NAOS(:) =',(NAOS(I),I=1,NSYM)
              CALL QUIT('SETHER error: NBAS(:) .ne. NAOS(:)')
           END IF
         ELSE
C
C        In relativistic calculation we must generate what is needed in INFORB
C
           CALL RELINF
         ENDIF
         SFIRST = .FALSE.
      END IF
      RETURN
      END
C  /* Deck herfck */
      SUBROUTINE HRFCK1(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,PARFCK,IPRFCK,
     &                  WORK,LWORK)
C
C     Written by Henrik Koch and Trygve Helgaker 22-November-1991.
C
C     PURPOSE : Driver routine for direct calculation of the
C               two-electron part of the fock matrices.
C               We assume the densities and fock matrices are full
C               squares and without symmetry reduction .
C
C               Common blocks needed to be initialized before this
C               subroutine is called:
C
C                      *CALL CBITWO   from TWOINP
C                      *CALL HRUNIT   from HERMIT
C                      *CALL INFORB   from sirius
C                      *CALL GAMCOM   from READIN
C                      *CALL CCOM     from READIN
C                      *CALL PINCOM   from READIN
C                      *CALL PRIMIT   from READIN
C                      *CALL XYZPOW   from READIN
C                      *CALL SHELLS   from READIN
C                      *CALL NUCLEI   from READIN
C                      *CALL SYMMET   from READIN
C                      *CALL DORPS    from abatro.u
C                       ( only DOREPS and DOCOOR used)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <dummy.h>
      LOGICAL   PARFCK, ABA
      DIMENSION FMAT(N2BASX,NDMAT), DMAT(N2BASX,NDMAT),
     &          ISYMDM(NDMAT), IFCTYP(NDMAT)
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  INFORB: N2BASX,NBAST,?
C  CBITWO: IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,?
C  INFPAR: KPRINT,MAXNOD,MAXTSK,?
C
#include <abainf.h>
#include <inforb.h>
#include <cbitwo.h>
#include <dorps.h>
#include <hrunit.h>
#include <symmet.h>
#include <infpar.h>
#include <ccom.h>
C
      I2TYP = 0
C
C     Inquire if this is an ABACUS run or not:
C
      CALL ABARUN(ABA)
C
C
C----------------------------------------------------------------
C     Setup information for the two-electron integral routines.
C----------------------------------------------------------------
C
      IF (.NOT.ABA) THEN
         LUONEL = 16
         LUINTA = 11
         LUSUPM = 17
         LUDASP = 24
         LUSOL  = 31
C
         NODV = (NASHT .EQ. 0)
         NOPV = (NASHT .LT. 2)
C
C        Remember to change this!!!!!!!
C
         IPRALL =  0
C
         DO 100 I = 0,7
            DOREPS(I) = .TRUE.
  100    CONTINUE
         DO 110 I = 1,MXCENT
            DOCOOR(1,I) = .TRUE.
            DOCOOR(2,I) = .TRUE.
            DOCOOR(3,I) = .TRUE.
  110    CONTINUE
C
         IPRINT = IPRFCK
         IPRNTA = 0
         IPRNTB = 0
         IPRNTC = 0
         IPRNTD = 0
         RETUR  = .FALSE.
         NOCONT = .FALSE.
         TKTIME = .FALSE.
      END IF
C
      ITYPE  = 3
      MAXDER = 0
      CALL PAOVEC(WORK,LWORK,0,IPRALL)
C
C     Transform density to AO basis
C
      KDAO  = 1
      KLAST = KDAO + N2BASX*NDMAT
      IF(KLAST.GT.LWORK) CALL STOPIT('HERFCK','DSOTAO',KLAST,LWORK)
      LWRK  = LWORK - KLAST + 1
      DO 200 IMAT = 1,NDMAT
         JDAO = KDAO + (IMAT-1)*N2BASX
         CALL DSOTAO(DMAT(1,IMAT),WORK(JDAO),NBAST,ISYMDM(IMAT),IPRINT)
  200 CONTINUE
C
C--------------------------------
C     Calculate fock matrices.
C--------------------------------
C
#if defined (VAR_MPI) 
      IF (PARFCK .OR. PARALL) THEN
         IATOM = 0
         KNSTAT = KLAST
         KINDEX = KNSTAT + MAXNOD
         KSORT  = KINDEX + MAXTSK
         KTIMES = KSORT  + MAXTSK + 2
         KLAST  = KTIMES + MAXTSK
         LWRK   = LWORK - KLAST
C
         IF (KLAST .GT. LWORK) THEN
            CALL STOPIT('HERFCK','PARDRV',KLAST,LWORK)
         END IF
C
         IF (KPRINT .GT. 0) CALL GETTIM(CPU1,WALL1)
C
         CALL PARDRV(WORK(KLAST),LWRK,FMAT,WORK(KDAO),NDMAT,ISYMDM,
     &               IFCTYP,ITYPE,MAXDER,IATOM,NODV,NOPV,NOCONT,
     &               TKTIME,RETUR,WORK(KNSTAT),WORK(KINDEX),
     &               WORK(KSORT),WORK(KTIMES))
C
C        Symmetrize skeleton Fock matrices
C
         CALL SKLFCK(FMAT,WORK(KLAST),LWRK,IPRINT,.TRUE.,.FALSE.,
     &               .FALSE.,.FALSE.,NODV,MAXDER,.FALSE.,NDMAT,
     &               ISYMDM,IFCTYP,IATOM)
C
         IF (KPRINT .GT. 0) THEN
            CALL GETTIM(CPU2,WALL2)
            CPU    = CPU2 - CPU1
            WALL   = WALL2 - WALL1
            TOTWAL = TOTWAL + WALL
C
            IMINS  = INT(TOTWAL)/60
            IHOURS = IMINS/60
            IMINS  = IMINS - 60*IHOURS
            ISECS  = NINT(TOTWAL) - 3600*IHOURS - 60*IMINS
C
            WRITE(6,'(A,I5.2,A,I2.2,A,I2.2)')
     &           '>>>> Total wall time used in PARDRV so far   :',
     &              IHOURS,':',IMINS,':',ISECS
            IF (KPRINT .GT. 5) THEN
               WRITE(LUPRI,'(A,F11.2,/A,F11.2)')
     &           '>>>> CPU  time used in PARDRV last iteration :',
     &             CPU,
     &           '>>>> Wall time used in PARDRV last iteration :',
     &             WALL
            END IF
         END IF
      ELSE
#endif
         CALL TIMER('START ',TIMSTR,TIMEND)
         CALL TWOINT(WORK(KLAST),LWRK,FMAT,WORK(KDAO),NDMAT,ISYMDM,
     &               IFCTYP,DUMMY,IDUMMY,NUMDIS,1,ITYPE,MAXDER,0,
     &               NODV,NOPV,NOCONT,TKTIME,IPRINT,IPRNTA,
     &               IPRNTB,IPRNTC,IPRNTD,RETUR,IDUMMY,I2TYP,
     &               IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,IDUMMY,.FALSE.)
         CALL TIMER('TWOINT',TIMSTR,TIMEND)
#if defined (VAR_MPI)
      END IF
#endif
C
C------------------------------------------------------
C     Write out densities and associated fock matrices.
C------------------------------------------------------
C
      IF (IPRINT.GT.4) THEN
         CALL HEADER('Density and Fock matrices in HERFCK',-1)
         DO 300 I = 1, NDMAT
            WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
            CALL OUTPUT(DMAT(1,I),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(1,I),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  300    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck aindex */
      SUBROUTINE AINDEX(ISHELA,NAINTS,INDEXA,DOINDX,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <aovec.h>
#include <mxorb.h>
#include <mxcent.h>
      LOGICAL DOINDX
      DIMENSION INDEXA(*)
#include <blocks.h>
#include <symmet.h>
#include <ibtfun.h>
C
      IORB = 0
      MULA   = ISTBSH(ISHELA)
      DO 100 IORBA = 1, NORBSH(ISHELA)
         DO 200 ICOMPA = 1, KHKTSH(ISHELA)
            ITYNA = ISYMAO(NHKTSH(ISHELA),ICOMPA)
            IADR  = IORBSB(IORBSH(ISHELA,1))
     &            + KHKTSH(ISHELA)*(IORBA - 1) + ICOMPA
            DO 300 IREPA = 0, MAXREP
               IF (IBTAND(MULA,IBTXOR(IREPA,ITYNA)) .EQ. 0) THEN
                  IORB = IORB + 1
                  IF (DOINDX) INDEXA(IORB) = IPTSYM(IADR,IREPA)
               END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      NAINTS = IORB
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from AINDEX',2)
         IF (.NOT.DOINDX) THEN
            WRITE (LUPRI,'(2X,A,I5)') ' Number of a orbitals:',NAINTS
         ELSE
            WRITE (LUPRI,'(2X,I3,A,8I5/,(40X,8I5))')
     &          NAINTS,' distributions in this TWOINT call:',
     &            (INDEXA(I),I=1,NAINTS)
         END IF
      END IF
      RETURN
      END
C  /* Deck abarun */
      SUBROUTINE ABARUN(RUNABA)
#include <implicit.h>
      LOGICAL ABA, RUNABA
      SAVE ABA
      DATA ABA /.FALSE./
      RUNABA = ABA
      RETURN
      ENTRY ABASET
      ABA = .TRUE.
      RETURN
      END
C  /* Deck herdi1 */
#if defined (VAR_HERDI1)
----> 941103-hjaaj: HERDI1 is incomplete(?) and is never called
----> therefore I have excluded HERDI1
      SUBROUTINE HERDI1(WORK,LWORK,IPRFCK)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <aovec.h>
#include <dummy.h>
      LOGICAL   ABA
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  INFORB: N2BASX,NBAST,?
C  CBITWO: IPRINT,?
C
#include <inforb.h>
#include <cbitwo.h>
#include <dorps.h>
#include <hrunit.h>
#include <symmet.h>
#include <blocks.h>
C
C     Inquire if this is an ABACUS run or not:
C
      CALL ABARUN(ABA)
C
C----------------------------------------------------------------
C     Setup information for the two-electron integralroutines.
C----------------------------------------------------------------
C
      LUINTA = 42
      IF (.NOT.ABA) THEN
         LUONEL = 16
         LUSUPM = 17
         LUDASP = 24
         LUSOL  = 31
C
         NODV = (NASHT .EQ. 0)
         NOPV = (NASHT .LT. 2)
C
C        Remember to change this!!!!!!!
C
         IPRALL =  0
C
         DO 100 I = 0,7
            DOREPS(I) = .TRUE.
  100    CONTINUE
         DO 110 I = 1,MXCENT
            DOCOOR(1,I) = .TRUE.
            DOCOOR(2,I) = .TRUE.
            DOCOOR(3,I) = .TRUE.
  110    CONTINUE
C
         IPRINT = IPRFCK
         IPRNTA = 0
         IPRNTB = 0
         IPRNTC = 0
         IPRNTD = 0
         RETUR  = .FALSE.
         NOCONT = .FALSE.
         TKTIME = .FALSE.
      END IF
C
      OPEN (LUINTA,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE='AOTWODIS')
      REWIND LUINTA
C
      CALL PAOVEC(WORK,LWORK,0,IPRALL)
C
      RETURN
      END
#endif
C  /* Deck herdi2 */
#if defined (VAR_HERDI2)
----> 941103-hjaaj: HERDI2 is incomplete and is never called
----> therefore I have excluded HERDI2
      SUBROUTINE HERDI2(WORK,LWORK,INDEXA,ISHELA,NAINTS,IPRFCK)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <aovec.h>
#include <dummy.h>
      LOGICAL DOINDX
      DIMENSION INDEXA(MXCORB)
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  INFORB: N2BASX,NBAST,?
C  CBITWO: IPRINT,?
C
#include <inforb.h>
#include <cbitwo.h>
#include <dorps.h>
#include <hrunit.h>
#include <symmet.h>
#include <blocks.h>
#include <ccom.h>
C
      I2TYP = 0
C
      ITYPE  = 7
      MAXDIF = 0
      JPRINT = 6
      DOINDX = .TRUE.
      CALL PAOVEC(0,IPRALL)
      CALL AINDEX(ISHELA,NAINTS,INDEXA,DOINDX,JPRINT)
      CALL TIMER('START ',TIMSTR,TIMEND)
      CALL TWOINT(WORK,LWORK,DUMMY,DUMMY,NDMAT,ISYMDM,
     &            IFCTYP,DUMMY,IDUMMY,NUMDIS,1,ITYPE,MAXDIF,0,NODV,
     &            NOPV,NOCONT,TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,
     &            IPRNTD,RETUR,ISHELA,I2TYP,IDUMMY,IDUMMY,
     &            DUMMY,DUMMY,DUMMY,IDUMMY,.FALSE.)
      CALL TIMER('TWOINT',TIMSTR,TIMEND)
C
      RETURN
      END
#endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck RELINF */
      SUBROUTINE RELINF
C*****************************************************************************
C
C     Generate INFORB for DIRAC
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
C
#include <nuclei.h>
#include <inforb.h>
#include <symmet.h>
C
      NSYM  = MAXREP + 1
      NBAST = NBASIS
      N2BASX = NBAST*NBAST
      DO I = 1,NSYM
        NBAS(I) = NAOS(I)
      ENDDO
C
      RETURN
      END
