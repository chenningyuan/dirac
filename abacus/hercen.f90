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
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : hercen.F
C
C
C  /* Deck ceninp */
      SUBROUTINE CENINP(LUINFO,WORK,LWORK,NONTYP,NSYMOP,NONT,IQM,NBLCK,
     &                  JCO,BASREF,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,ISGEN,
     &                  KATOM,KANG,KSETS,KBLOCK,KPRIM,ANG,ZMTTST)
C*****************************************************************************
C
C       Loop over expansion centers and read coordinates and basis set information.
C
C       Revised version of BASINP, L. Visscher, december 2003.
C
C       Descripption of some of the arguments :
C       LUINFO : File from which the input is read (MOLECULE.INP, the .mol file in pam)
C       NONTYP : Number of different atom types defined by the user.
C       NONTYP : Number of different atom types defined by the user.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0)
#include "cbirea.h"
      DIMENSION WORK(LWORK)
      LOGICAL SEG,ANG,NOORBTS,NOFITFS,ZMTTST,CNTBAS
      DIMENSION NONT(KATOM),IQM(KATOM,KSETS),JCO(KANG,KATOM,KSETS),
     &          NBLCK(KATOM,KSETS),NUC(KBLOCK,KSETS),NRC(KBLOCK,KSETS),
     &          SEG(KBLOCK,KSETS),
     &          ALPHA(KPRIM,KBLOCK,KSETS),
     &          CPRIM(KPRIM,KPRIM,KBLOCK,KSETS),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS),ISGEN(KBLOCK)
      CHARACTER SPDCAR*1,BSET*5,FBSET*6,ZMATL*1,BASFIL*20,BSKEYWORD*14
      CHARACTER*80 BASREF(10,KATOM,KSETS)
      LOGICAL GENERATE
      DIMENSION IBLCK(MXBSETS)
#include "molinp.h"
#include "nuclei.h"
#include "huckel.h"
#include "gencon.h"
#include "dcbgen.h"
#include "ibtfun.h"
      CNTBAS = .FALSE.
      NUCIND = 0
      NHTYP  = 0
      CALL IZERO (NBLCK,KATOM*NSETS)
C
C     Run over atomic centers
C     =======================
C
      DO 100 I = 1, NONTYP
C        Define some offsets
         INUC = NUCIND   ! start of information for nuclei
         DO ISET = 1, NSETS
            IBLCK(ISET) = NBLCK(I,ISET) + 1 ! start of orbital blocks
         ENDDO
C
C*****************************************************************************
C        Card 6: ATOMIC DATA
C         Q     - nuclear charge
C         NONT  - number of symmetry independent atoms
C         IQM   - maximum angular quantum number (s=1,p=2,d=3 etc.)
C         NHTYP - maximum angular quantum number for ALL orbitals
C         JCO   - number of AO-blocks for each l-value
C*****************************************************************************
C
         CALL RDLINE(LUINFO)
C
C        Test if Z-matrix for first line (I.eq.1):
C
         IF (I .EQ. 1 .AND. .NOT. DIRAC) THEN
           READ (MLINE(NMLINE),'(A1)') ZMATL
           IF (ZMATL .NE. ' ') THEN
              ZMTTST = .TRUE.
              WRITE (LUPRI,'(/A)')
     &        'Atoms and their geometry is read in Z-matrix format.'
              IF (.NOT. BASIS) THEN
                 WRITE (LUPRI,'(/A)') '*** READI1 ERROR:'
                 CALL QUIT('Z-matrix only implemented with BASIS')
              END IF
           END IF
         END IF
C
C       Read/process data on centers/blocks
C       ===================================
C
        QEXP = D0
        IF(.NOT.ZMTTST) THEN
          IF (DIRAC) THEN
            ADDSTO = .FALSE.
            READ (MLINE(NMLINE),'(BN,1X,F9.0,I5,F20.5)')
     &            Q,NONT(I),QEXP
          ELSE IF (.NOT.ZMTTST) THEN
            NSETS = 1
            IF (ADDSTO) NSETS = 2
            BSET  = 'Basis'
            READ (MLINE(NMLINE),'(BN,1X,F9.0,25I5)') Q,NONT(I),
     &                           IQM(I,1), (JCO(J,I,1), J=1, IQM(I,1))
            IF (NINT(Q) .GT. 18) THEN
              NSETS = 1
              ADDSTO = .FALSE.
              WRITE (LUPRI,'(/A)') ' No Huckel start guess due to '//
     &             'nuclear charge larger than 18'
            END IF
          ENDIF
        ENDIF
C
C       Read/process data on symmetry independent centers
C       =================================================
C
        IF (ZMTTST) THEN
           CALL ZMAT(NONT,KATOM)
           Q = CHARGE(NUCIND)
        ELSE
           IF (INPTST) THEN
              WRITE(LUPRI,'(/A,I5)') 'Atomic type: ',I
              WRITE(LUPRI,'(3X,A,F4.1)') 'Nuclear charge:  ',Q
              WRITE(LUPRI,'(3X,A,I5)') 'Number of symmetry '//
     &             'independent centres: ',NONT(I)
           ENDIF
C
C          Read coordinates of the centers of that belong to this atom type
C
           CALL CNTINP(LUINFO,NONT(I),Q,QEXP,NSYMOP,ANG)
C
C          Loop over the basis sets that are to be associated with this atom type.
C
           IF (DIRAC) THEN
              NSETS = 4
           ELSE
              NSETS = 3
           ENDIF
           DO ISET = 1, NSETS
              KLAST = 1
              LWRK = LWORK
C             IBLOCK gives the first free block for this atom type and set.
              IBLOCK = IBLCK(ISET)
              CALL READ_BASIS (LUINFO,WORK(KLAST),LWRK,
     &             IQM(I,ISET),NBLCK(I,ISET),JCO(1,I,ISET),
     &             BASREF(1,I,ISET),
     &             NUC(IBLOCK,ISET),NRC(IBLOCK,ISET),SEG(IBLOCK,ISET),
     &             ALPHA(1,IBLOCK,ISET),CPRIM(1,1,IBLOCK,ISET),
     &             CPRIMU(1,1,IBLOCK,ISET),ISGEN,
     &             KANG,ISET,KBLOCK,KPRIM,GENERATE)
C
C             A set might be generated from a previous set for this atom type.
C
              IF (GENERATE) THEN
C
                  IBLOCK1 = IBLCK(1)
                  IF (ISET.EQ.2) THEN
C                    When uncontracted we generate the small component basis
C                    using the kinetic balance relation
C                    ...or...
C                    Make uncontracted small component basis
C                    set from kinetic balance.
C                    =======================================================
                     CALL KINBAL(I,IBLOCK1,IBLOCK,
     &                    IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &                    NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,INPTST,
     &                    CNTBAS)
                     BASREF(1,I,ISET) = ' Derived from large component'
                  ELSE
C                    TO DO : automatic generation of fit set ?
                  ENDIF
              ENDIF
           ENDDO
        ENDIF
C
C     End of loop over unique atom types
  100 CONTINUE
C
      RETURN
      END
C  /* Deck READ_BASIS */
      SUBROUTINE READ_BASIS (LUINFO,WORK,LWORK,IQM,NBLCK,
     &           JCO,BASREF,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,ISGEN,
     &           KANG,ISET,KBLOCK,KPRIM,GENERATE)
C*****************************************************************************
C
C       Read an expansion set for a given center
C       Note : NBLCK is the total number of blocks for this type of basis set
C              NBLOCK is the number of blocks for this type and atom
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (DSM = 1.0D-30)
#include "cbirea.h"
      DIMENSION WORK(LWORK)
      LOGICAL SEG,NOORBTS,NOFITFS,CNTBAS
      DIMENSION JCO(KANG),NUC(KBLOCK),NRC(KBLOCK),
     &          SEG(KBLOCK),
     &          ALPHA(KPRIM,KBLOCK),
     &          CPRIM(KPRIM,KPRIM,KBLOCK),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK),ISGEN(KBLOCK)
      CHARACTER SPDCAR*1,BSET*5,BASFIL*20,BSKEYWORD*14
      CHARACTER*80 BASREF(10)
      LOGICAL GENERATE
#include "molinp.h"
#include "nuclei.h"
#include "huckel.h"
#include "gencon.h"
#include "dcbgen.h"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Additional scratch arrays needed for MULTIBASIS option
C
      PARAMETER (MAXBSFIL=5,KNG=MXAQN,KBLK=500,KPRM=30)
C
      INTEGER IQMT
      DIMENSION JCOT(KNG,MAXBSFIL)
      DIMENSION IBLKAGNMOM(KBLK)
      DIMENSION NUCT(KBLK)
      DIMENSION NRCT(KBLK)
      DIMENSION SEGT(KBLK)
      DIMENSION ALPHAT(KPRM,KBLK)
      DIMENSION CPRIMT(KPRM,KPRM,KBLK)
      DIMENSION CPRIMUT(KPRM,KPRM,KBLK)
      DIMENSION ISGENT(KBLK)
      CHARACTER*80 BASREFT(10,KATOM,MAXBSFIL)
      CHARACTER*20 BSFIL(MAXBSFIL)
      LOGICAL SEGT
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
#include "ibtfun.h"
      GENERATE = .FALSE.
      CNTBAS = .FALSE.
      NHTYP  = 0
      IQM = 0
C
C     Initialize reference strings for this atom type.
C
      BASREF(1) = ' Reference not found in input'
      DO IREF = 2, 10
         BASREF(IREF) = ' Not initialized'
      ENDDO
C
C     Read the basis set. The input differs in DALTON/DIRAC due to
C     the optional reading of small component basis sets. The code 
C     switch starts here with first the processing of the DIRAC format.
C
C     Construct basis set for dummy centers if requested
C     ==================================================
C
      NQ = NINT(Q + 0.01D0)
      IF ((LCMMAX .GT. -1) .AND. (NQ .EQ. 0)) THEN
        WRITE(LUPRI,'(/A)') 'BASINP: CM basis not implemented !'
        CALL QUIT('BASINP: CM basis not implemented in Dirac !')
C
C     Get basis set from basis set library
C     ===================================
C
      ELSE IF (BASIS) THEN
        WRITE(LUPRI,'(/A/A)')
     &    'BASINP: "BASIS" in first line not allowed in Dirac.',
     &    'Use "LARGE BASIS <basis>" for each atom type instead.'
        CALL QUIT('BASINP: "BASIS" in first line not allowed in Dirac!')
C
      ELSE
C
C     Read basis set info from file MOLECULE.INP
C     ==========================================
C
C        The input is handled by different routines, depending on the
C         keywords that are given. The first gives the type of the set
C         of functions that is read in. This shoold correspond to the set
C         parameter ISET that is given in the call to BASINP. Valid
C         names are :
C
C         ISET     BSKEYWORD   Type of functions
C           1      LARGE       Large component basis functions
C           2      SMALL       Small component basis functions
C           3      FTSET       Primary set for density fitting (total density matrix)
C           3      FTLRG       Primary set for density fitting (LL density matrix)
C           4      FTSML       Secundary set for density fitting (SS density matrix)
C          (note that FTLRG is an alias for FTSET at the moment)
C
C         The next keyword gives the type of basis set. Valid options are :
C
C         IBSFLAG  BSKEYWORD   Type of functions
C           1      (blank)     Read from input file in standard format
C           1      INTGRL      Read from input file in standard format
C           1      EXPLICIT    Read from input file in standard format
C           2      MOLFBAS     Read from specified basisfile in MOLFDIR format
C           3      BASIS       Read from specified basisfile in DALTON  format
C           4      GEOM        Read from input in even tempered format (alpha beta N IQM JCO())
C           4      EVENTEMP    Read from input in even tempered format (alpha beta N IQN JCO())
C           4      WELLTEMP    Read from input in well tempered format (alpha beta gamma delta N IQM JCO())
C           4      FAMILY      Read from input in family basis format (N IQM JCO())
C           4      DUALFAMILY  Read from input in dual family basis format (N IQM JCO())
C           5      MULTIBASIS  Read from multiple basis sets in DALTON format (used to augment sets)
C           6      KINBAL      Generate the small component basis by kinetic balance
C
         CALL RDLINE(LUINFO)
         READ (MLINE(NMLINE),'(A5,X,A14)') BSET,BSKEYWORD
C
         SELECT CASE (BSET)
         CASE ("LARGE")
              JSET = 1
         CASE ("SMALL")
              JSET = 2
         CASE ("FTSET")
              JSET = 3
         CASE ("FTLRG")
              JSET = 3
         CASE ("FTSML")
              JSET = 4
         CASE DEFAULT
              JSET = 0
         END SELECT
C
C        Check the correct order of the input. We do supply defaults,
C        but want to detect trivial mistakes before doing so.
C
         IF (JSET.EQ.0.AND.ISET.EQ.1) THEN
C           The large component should always be specified
            WRITE(LUPRI,'(A,I5,A,/,A)')'Expected LARGE keyword, read',
     $      ' line no. :',NMLINE,':',MLINE(NMLINE)
            CALL QUIT('Error reading keyword in BASINP')
         ELSEIF (JSET.LT.ISET) THEN
             WRITE(LUPRI,'(A,A)')'Wrong order of keywords, specify',
     $       ' them in the order LARGE SMALL FTSET/FTLRG FTSML'
            CALL QUIT ('BASINP : Keywords in wrong order')
         ENDIF
C
C        Initialize basis set flag
C
         SELECT CASE (BSKEYWORD)
         CASE ("INTGRL")
              IBSFLAG = 1
              READ (MLINE(NMLINE),'(12X,12I5)')
     &              IQM,(JCO(K),K=1,IQM)
         CASE ("EXPLICIT")
              IBSFLAG = 1
              READ (MLINE(NMLINE),'(14X,12I5)')
     &             IQM,(JCO(K),K=1,IQM)
         CASE ("MOLFBAS")
              IBSFLAG = 2
              READ (MLINE(NMLINE),'(14X,A20)') BASFIL
         CASE ("BASIS")
              IBSFLAG = 3
              READ (MLINE(NMLINE),'(12X,A20)') BASFIL
         CASE ("WELLTEMP")
              IBSFLAG = 4
              READ (MLINE(NMLINE)(16:),*,ERR=9000)
     $             (FAMPAR(K),K=1,4),NFAMEXP(1),
     $             IQM, (JCO(K),K=1,IQM)
              CALL GENFAMEXP(LUINFO,0)
         CASE ("EVENTEMP")
              IBSFLAG = 4
              READ (MLINE(NMLINE)(16:),*,ERR=9000)
     $             (FAMPAR(K),K=1,2),NFAMEXP(1),
     $             IQM, (JCO(K),K=1,IQM)
              CALL GENFAMEXP(LUINFO,0)
         CASE ("GEOM")
              IBSFLAG = 4
              READ (MLINE(NMLINE)(11:),*,ERR=9000)
     $             (FAMPAR(K),K=1,2),NFAMEXP(1),
     $             IQM, (JCO(K),K=1,IQM)
              FAMPAR(3) = D0
              FAMPAR(4) = D0
              CALL GENFAMEXP(LUINFO,0)
         CASE ("FAMILY")
              IBSFLAG = 4
              READ (MLINE(NMLINE),'(12X,12I5)',ERR=9000)
     $             NFAMEXP(1), IQM, (JCO(K),K=1,IQM)
              CALL GENFAMEXP(LUINFO,1)
         CASE ("DUALFAMILY")
              IBSFLAG = 4
               READ (MLINE(NMLINE),'(16X,12I5)',ERR=9000)
     $              NFAMEXP(1), NFAMEXP(2),
     &              IQM, (JCO(K),K=1,IQM)
              CALL GENFAMEXP(LUINFO,2)
         CASE ("MULTIBASIS")
              IBSFLAG = 5
              READ (MLINE(NMLINE)(18:),*,ERR=9000)
     &               NFILES,(BSFIL(K),K=1,NFILES)
              IF ( NFILES .EQ. 0 ) THEN
                 WRITE(*,*) 
     &             'PLEASE SPECIFY THE NUMBER OF BASIS SET FILES'
                 GOTO 9000
              END IF
         CASE ("KINBAL")
              IF (ISET.EQ.2) THEN
                 IBSFLAG = 6
              ELSE
                 CALL QUIT ('KINBAL only defined for SC basis')
              ENDIF
         CASE DEFAULT
C            The default depends on the set that we are reading.
             SELECT CASE (ISET)
             CASE (1)
C               Assume it is like 'LARGE   2    1....'
                IBSFLAG = 1
                READ (MLINE(NMLINE),'(5X,12I5)',ERR=9000)
     &               IQM,(JCO(K),K=1,IQM)
             CASE (2)
C               Assume kinetic balance, the user has omittted this line
                IBSFLAG = 6
                NMLINE = NMLINE - 1
             CASE (3,4)
C               Assume that no fit set was given, the user has omittted this line
                IBSFLAG = 0
                NMLINE = NMLINE - 1
             END SELECT
         END SELECT
C
C        We know what should be read, call the appropriate reader.
C
C        Compute number of storage blocks remaining before we start
         KAOVEC = KBLOCK - IBLOCK + 1
C
         SELECT CASE (IBSFLAG)
         CASE (0)
C
C            No basis set given, no need to read.
C            We have not set IQM yet in this case, do this here.
             IQM = 0
C
         CASE (1)
C
C            Read explicitly given basis set.
C
             CALL GTOINP(LUINFO,IQM,JCO,
     &            NUC,NRC,SEG,
     &            ALPHA,CPRIM,CPRIMU,ISGEN,
     &            NBLOCK,KAOVEC,KPRIM)
             NBLCK = NBLCK + NBLOCK
             BASREF(1) = 
     &       "Basis set typed explicitly in input file "
         CASE (2)
C
C            Molfdir basis set
C
             KBFPR = 1
             KBFCO = KBFPR + KPRIM
             KLAST = KBFCO + KPRIM*KPRIM
             IF (KLAST.GT.LWORK)
     &          CALL STOPIT('BASINP','BFGINP',KLAST,LWORK)
             DIRCON(1) = .TRUE.
             CALL BFGINP(I,ISET,WORK(KBFPR),WORK(KBFCO),
     &           IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &           NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,BASFIL,
     &           BASREF)
         CASE (3)
C
C            Dalton basis set
C
             BASNAM = BASFIL
             CALL BASLIB(IQM,JCO,NUC,NRC,SEG,
     &            ALPHA,CPRIM,CPRIMU,ISGEN,
     &            NBLOCK,ISGEN,KAOVEC,KPRIM,Q,1,DSM,
     &            BASREF)
             NBLCK = NBLCK + NBLOCK
         CASE (4)
C
C            Generate well-temperered/family basis set
C
             CALL FAMBAS(LUINFO,IQM,JCO,NUC,NRC,SEG,
     &            ALPHA,CPRIM,CPRIMU,ISGEN,
     &            NBLOCK,KAOVEC,KPRIM)
             NBLCK = NBLCK + NBLOCK
             BASREF(1) = 
     &   "Well-tempered basis set typed explicitly in input file "
C
         CASE (5)
C
C TODO (10-12-2003, LV) : Isolate this code in a subroutine.
C
C reading of basis sets from a number of files, including when the files
C do not necessarily have the same maximum and/or minimum orbital 
C angular momenta (last modif. 11/09/2003 by gomes) 
C
C we will move through the basis set files using the dummy variable iblock1
C 
                IBLOCK1 = 1
                MAXIQM  = 0
C
C we make sure we don't have any junk on the temporary variables 
C
                CALL IZERO(JCOT,KNG*NFILES)
                CALL IZERO(NUCT,KBLK)
                CALL IZERO(NRCT,KBLK)
C
                BASREF(1) =
     &          "Basis set read from multiple files"
C
                DO 463 M=1,NFILES
                   BASNAM = BSFIL(M)
                   KAOVEC  = KBLOCK + 1 - IBLOCK1 
C
                   BASREF(3*(M-1)+2) = " " 
                   BASREF(3*(M-1)+3) = " ref. in file "//BSFIL(M) 
C 
                   CALL BASLIB(IQMT,JCOT(1,M),
     &                  NUCT(IBLOCK1),
     &                  NRCT(IBLOCK1),SEGT(IBLOCK1),
     &                  ALPHAT(1,IBLOCK1),
     &                  CPRIMT(1,1,IBLOCK1),
     &                  CPRIMUT(1,1,IBLOCK1),
     &                  NBLOCK,ISGENT(IBLOCK1),KAOVEC,KPRIM,Q,1,DSM,
     &                  BASREFT(1,I,M))

                   BASREF(3*(M-1)+4) = BASREFT(1,I,M)
C
C to help in the reordering, we create a table mapping the orbital
C angular momentum of a block to the number of the block in the
C temp storage 
C 
                   IDBLK = IBLOCK1
                   IDANG = 1
                   DO 464 IRR=1,IQMT
                      DO 465 IRS=1,JCOT(IRR,M)
                         IBLKAGNMOM(IDBLK) = IDANG
                         IDBLK = IDBLK + 1
  465                 CONTINUE 
                      IDANG = IDANG + 1
  464              CONTINUE 
C
                   IF ( IQMT .GT. MAXIQM ) THEN
                      MAXIQM = IQMT
                   ENDIF
                   IBLOCK1 = IBLOCK1 + NBLOCK
  463           CONTINUE
C
C after reading all files, we reorder the data wrt orbital angular 
C momentum. in this reordering variable kt is the "final" index of a given block
C
                IQM = MAXIQM
                KT = IBLOCK
                DO 466 IQ=1,MAXIQM
                DO 467 IR=1,IBLOCK1 - 1
C
C in filling up the permanent storage with data, care should be taken to skip 
C zero-sized blocks that arise when we have only parts of a basis
C set in one file, as when the polarization functions are in one file and
C the scf set is in another...
C
                     IF ( IBLKAGNMOM(IR) .EQ. IQ ) THEN
                       IF ( NUCT(IR) .GT. 0 ) THEN

                            JCO(IQ)= JCO(IQ)+ 1
                            NUC(KT)  = NUC(KT) + NUCT(IR)
                            NRC(KT)  = NRC(KT) + NRCT(IR)
                            SEG(KT)  = SEGT(IR)
                            ISGEN(KT)  = ISGENT(IR)

                            DO 468 IX=1, NUCT(IR)
                               ALPHA(IX,KT)      =  ALPHAT(IX,IR)
                               DO 468 JX=1, NRCT(IR) 
                                CPRIM(IX,JX,KT) =  CPRIMT(IX,JX,IR)
                                CPRIMU(IX,JX,KT)= CPRIMUT(IX,JX,IR)
  468                       CONTINUE
                            KT = KT + 1
                         ENDIF 
                      ENDIF
  467              CONTINUE
  466           CONTINUE
C
C finally, we update nblck(i,1) with the total number of non-zero sized
C blocks read, subtracting one only because on the loops above kt was 
C incremented after the last block was put in place ...
C
               NBLCK = KT - 1 
C
C ... and we are done
C
         CASE (6)
C
C           Kinetic balance is special as it requires information
C           about the other basis that we do not have at this level.
C           Do this at a higher level.
C
            GENERATE = .TRUE.
         END SELECT
C
C        Check whether the basis set does not exceed the limit in angular momentum
C
         NOORBTS    = IQM.EQ.0
         IF (.NOT.NOORBTS) THEN
             NHTYP     = MAX(NHTYP,IQM)
             IF(NHTYP.GT.MXQN) GOTO 5000
         ENDIF
C
C        Check if basis set is contracted
C        ============================================
C
         DO JBLCK = 1, NBLOCK
            IF ((NUC(JBLCK).NE.NRC(JBLCK)).AND.
     &          (NRC(JBLCK).GT.0)) CNTBAS = .TRUE.
         END DO
         IF ( CNTBAS .AND.
     $        (IBSFLAG .EQ. 3 .OR. IBSFLAG .EQ. 1 )) THEN
            DIRCON = .TRUE.
         END IF
C
C     End of IF statement choosing between global set and individualy specified sets
      ENDIF
C
      RETURN
C
C     Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(6X,A,I3,3A/9X,2A/9X,2(A,I3),A)')
     &     '*  Input specifies highest orbital of atomic type ',
     &     I,' as "',SPDCAR(NHTYP-1),'".',
     &     ' Highest allowed orbital in this version: ',
     &     SPDCAR(MXQN-1),
     &     ' Increase MXQN from',MXQN,' to',NHTYP/2,' and recompile.'
        CALL QUIT('Too high angular specified in input.')
 5001 CONTINUE
        WRITE (LUPRI,'(6X,A,I3,3A/9X,2A/9X,2(A,I3),A)')
     &     '*  Input specifies highest fit function of atomic type ',
     &     I,' as "',SPDCAR(NHTYPF-1),'".',
     &     ' Highest allowed fit function in this version: ',
     &     SPDCAR(2*MXQN-2),
     &     ' Increase MXQN from',MXQN,' to',NHTYPF,' and recompile.'
        CALL QUIT('Too high angular specified in input.')
 5200 CONTINUE
        CALL QUIT('Error reading keyword SMALL or SBFIL.')
 9000 CONTINUE
      WRITE(LUPRI,'(A,/,A,I3,/,A)')
     $     '*** ERROR in READIN ***',
     &     'Error reading line ',NMLINE,
     &     MLINE(NMLINE)
         CALL QUIT('*** ERROR in READIN ***')
      END
