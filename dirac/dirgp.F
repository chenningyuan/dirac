C 
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
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
C
C VERSION : $Revision: 1.45.2.1 $
C DATE    : $Date: 2004/09/01 08:08:47 $
C FILE    : dirgp.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pcmoin */
      SUBROUTINE PCMOIN(WORK,LWORK)
C*****************************************************************************
C
C     Read coefficients from formatted DFPCMO and
C     write to unformatted DFCOEF
C
C     Written by T.Saue Feb 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
C
      DIMENSION WORK(LWORK)
      CALL QENTER('PCMOIN')
#include "memint.h"
C
C     Memory allocation
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBG ,NORBT   ,WORK,KFREE,LFREE)
C
      CALL PCMOI1(WORK(KCMO),WORK(KEIG),WORK(KIBG),WORK(KFREE),LFREE)
C
C     Memory deallocation
      CALL MEMREL('PCMOIN',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('PCMOIN')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pcmout */
      SUBROUTINE PCMOUT(WORK,LWORK)
C*****************************************************************************
C
C     Read coefficients from unformatted DFCOEF and
C     write to formatted DFPCMO
C
C     Written by T.Saue Feb 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
C
      DIMENSION WORK(LWORK)
      CALL QENTER('PCMOUT')
#include "memint.h"
C
C     Memory allocation
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBG ,NORBT   ,WORK,KFREE,LFREE)
C
      CALL PCMOU1(WORK(KCMO),WORK(KEIG),WORK(KIBG),WORK(KFREE),LFREE)
C
C     Memory deallocation
      CALL MEMREL('PCMOUT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('PCMOUT')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wricmo */
      SUBROUTINE WRICMO(IUNIT,CMO,EIG,IBEIG,TOTERG)
C*****************************************************************************
C
C     Write DHF-coefficients and eigenvalues to unformatted file
C
C     IDIM(1,J) - number of positronic solutions
C     IDIM(2,J) - number of electronic solutions
C     IDIM(3,J) - number of AO-basis functions
C
C     Written by T.Saue Sept 1 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbgen.h"
      CHARACTER TEXT*74
      DIMENSION CMO(NCMOTQ),EIG(NORBT),IBEIG(NORBT),IDIM(3,2)
C
      DO I = 1,NFSYM
        IDIM(1,I) = NPSH(I)
        IDIM(2,I) = NESH(I)
        IDIM(3,I) = NFBAS(I,0)
      ENDDO
      TEXT(1:50) = TITLE
      CALL GTINFO(TEXT(51:74))
      REWIND IUNIT
      WRITE(IUNIT) TEXT,NFSYM,((IDIM(I,J),I = 1,3),J=1,NFSYM),
     &             TOTERG
      WRITE(IUNIT) CMO
      WRITE(IUNIT) EIG
      WRITE(IUNIT) IBEIG
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reacmo */
      SUBROUTINE REACMO(IUNIT,CMOFIL,CMO,EIG,IBEIG,TOTERG,IOPT,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     Read DHF-coefficients and eigenvalues from unformatted file
C     Read option is provided from bit-packed IOPT:
C       0001 - give restart info
C       0010 - read coefficients
C       0100 - read eigenvalues
C       1000 - read boson irrep identification
C
C     IDIM(1,J) - number of positronic solutions
C     IDIM(2,J) - number of electronic solutions
C     IDIM(3,J) - number of AO-basis functions
C
C     Written by T.Saue Sept 1 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
      LOGICAL INFO,LBIT,TEST
      CHARACTER TEXT*74,CMOFIL*6,FMT*6,MXFORM*6
      DIMENSION CMO(NCMOTQ),EIG(NORBT),IBEIG(NORBT),WORK(*)
      DIMENSION IDIM(3,2)
C
      IF(IOPT.EQ.0) RETURN
      REWIND IUNIT
C
C     Read title line
C     ===============
C
      READ (IUNIT,END=10,ERR=20) TEXT,NSYM,
     &       ((IDIM(I,J),I = 1,3),J=1,NSYM),TOTERG
C
C     Compatibility check for NFSYM and NFBAS
C     =======================================
C
      IERR = 0
      IF(NSYM .NE. NFSYM) THEN
         IERR = IERR + 1
         INFO = .TRUE.
         WRITE(LUPRI,'(/A/A,2I5)')
     &      'REACMO FATAL ERROR: Incompatible number of fermion ircops',
     &      ' from file and in this run:',NSYM,NFSYM
      ELSE
         DO IFRP = 1,NFSYM
           IF(IDIM(3,IFRP).NE.NFBAS(IFRP,0)) IERR = IERR + 1
         END DO
         IF (IERR.GT.0) THEN
           INFO = .TRUE.
           WRITE(LUPRI,'(/A/A,2(/A,I2,2I5))')
     &     'REACMO FATAL ERROR. Incompatible number of basis functions',
     &     ' from file and in this run:', (' Fermion ircop no.',IFRP,
     &     IDIM(3,IFRP),NFBAS(IFRP,0),IFRP=1,NFSYM)
         END IF
      END IF
C
C     File info
C     =========
C
      IF(LBIT(IOPT,1)) THEN
        FMT = MXFORM(TOTERG,20)
        WRITE(LUPRI,'(/3A,'//FMT//')')
     &   '* REACMO: Coefficients read from file ',CMOFIL,
     &   ' - Total energy: ',TOTERG
        WRITE(LUPRI,'(2A)') '* Heading :',TEXT
      ENDIF
C
      IF (IERR.GT.0) CALL QUIT(
     &'REACMO: old DFCOEF file is not consistent with this calculation')
C
      IF(IOPT.LT.2) RETURN
C
C     Compatibility tests for number of orbitals
C     ==========================================
C
      TEST = .FALSE.
      NC = 0
      NO = 0
      DO IFRP = 1,NFSYM
        IDIM(3,IFRP) = IDIM(1,IFRP)+IDIM(2,IFRP)
        IF(IDIM(3,IFRP).NE.NORB(IFRP)) TEST=.TRUE.
        NC = NC + IDIM(3,IFRP)*NFBAS(IFRP,0)*NZ
        NO = NO + IDIM(3,IFRP)
      ENDDO
C
C     Read coefficients
C     =================
C
      IF(LBIT(IOPT,2)) THEN
        IF(TEST) THEN
          CALL MEMGET('REAL',KC,NC,WORK,KFREE,LFREE)
          CALL READT(IUNIT,NC,WORK(KC))
          JC = KC 
          DO IFRP = 1,NFSYM
            CALL ADACMO(NFBAS(IFRP,0),NZ,
     &         CMO(ICMOQ(IFRP)+1),NORB(IFRP),NPSH(IFRP),NESH(IFRP),
     &         WORK(JC),IDIM(3,IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
            JC = JC + IDIM(3,IFRP)*NFBAS(IFRP,0)*NZ
          ENDDO
          CALL MEMREL('REACMO.cf',WORK,KC,KC,KFREE,LFREE)
        ELSE
          READ (IUNIT,END=30,ERR=40) CMO
        ENDIF
        IF(IOPT.LT.4) RETURN
      ELSE
        READ(IUNIT)
      ENDIF
C
C     Read eigenvalues
C     ================
C
      IF(LBIT(IOPT,3)) THEN
        IF(TEST) THEN
          CALL MEMGET('REAL',KO,NO,WORK,KFREE,LFREE)
          CALL READT(IUNIT,NO,WORK(KO))
          JO = KO 
          DO IFRP = 1,NFSYM
            CALL ADAEIG(
     &         EIG(IORB(IFRP)+1),NORB(IFRP),NPSH(IFRP),NESH(IFRP),
     &         WORK(JO),IDIM(3,IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
            JO = JO + IDIM(1,IFRP)+IDIM(2,IFRP)
          ENDDO
          CALL MEMREL('REACMO.eig',WORK,KO,KO,KFREE,LFREE)
        ELSE
          READ (IUNIT,END=50,ERR=60) EIG
        ENDIF
        IF(IOPT.LT.8) RETURN
      ELSE
        READ(IUNIT)
      ENDIF
C
C     Read single group irrep (spinfree and Levy-Leblond calcs)
C     =========================================================
C
      IF(LBIT(IOPT,4)) THEN
        IF(TEST) THEN
          CALL MEMGET('INTE',KO,NO,WORK,KFREE,LFREE)
          CALL READI(IUNIT,NO,WORK(KO))
          CALL ADAIBP(IBEIG,WORK(KO),IDIM)
          CALL MEMREL('REACMO.ibos',WORK,KO,KO,KFREE,LFREE)
        ELSE
          READ (IUNIT,END=70,ERR=80) IBEIG
        ENDIF
        IF(IOPT.LT.16) RETURN
      ELSE
        READ(IUNIT)
      ENDIF
C
      RETURN
C
C     Error messages
C
 10   CONTINUE
      CALL QUIT('REACMO: END OF FILE reading TEXT')
 20   CONTINUE
      CALL QUIT('REACMO: ERROR reading TEXT')
 30   CONTINUE
      CALL QUIT('REACMO: END OF FILE reading coefficients')
 40   CONTINUE
      CALL QUIT('REACMO: ERROR reading coefficients')
 50   CONTINUE
      CALL QUIT('REACMO: END OF FILE reading eigenvalues')
 60   CONTINUE
c     CALL QUIT('REACMO: ERROR reading eigenvalues')
      WRITE(LUPRI,'(/A)') 'REACMO: WARNING !!! No eigenvalues read'
      RETURN
 70   CONTINUE
      CALL QUIT('REACMO: END OF FILE reading irrep identification')
 80   CONTINUE
      CALL QUIT('REACMO: ERROR reading irrep identification')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pcmoi1 */
      SUBROUTINE PCMOI1(CMO,EIG,IBEIG,WORK,LWORK)
C*****************************************************************************
C
C     Read DHF-coefficients and eigenvalues from FORMATTED file
C
C     Written by T.Saue Feb 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
      LOGICAL TOBE
      CHARACTER TEXT*74
      DIMENSION CMO(NCMOTQ),EIG(NORBT),IBEIG(NORBT),WORK(*)
      DIMENSION IDIM(3,2)
C
C     Get coefficients from FDPCMO
C
      INQUIRE(FILE='DFPCMO',EXIST=TOBE)
      IF(.NOT.TOBE) THEN
          WRITE(LUPRI,'(A)') 
     &      'PCMOIN: Formatted file DFPCMO not found !'
          CALL QUIT('PCMOIN: Formatted file DFPCMO not found !')
      ENDIF
      OPEN(LUCOEF,FILE='DFPCMO',STATUS='OLD',FORM='FORMATTED',
     +     ACCESS='SEQUENTIAL')
      READ (LUCOEF,'(A74)') TEXT
      READ (LUCOEF,'(I4,6I10)') NSYM,((IDIM(I,J),I = 1,3),J=1,NSYM)
      READ (LUCOEF,'(D24.16)') TOTERG
      READ (LUCOEF,'(6F22.16)') (CMO(I),I=1,NCMOTQ)
      READ (LUCOEF,'(6D22.12)') (EIG(I),I=1,NORBT)
      READ (LUCOEF,'(66I2)') (IBEIG(I),I=1,NORBT)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Write coefficients to DFCOEF
C
      CALL OPNFIL(LUCOEF,'DFCOEF','UNKNOWN','PCMOIN')        
      CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,TOTERG)
      CLOSE(LUCOEF,STATUS='KEEP')
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)')
     &    '* PCMOIN: Coefficients read from formatted DFPCMO ',
     &    '          and written to unformatted DFCOEF'
      CALL PRSYMB(LUPRI,'=',75,0)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pcmou1 */
      SUBROUTINE PCMOU1(CMO,EIG,IBEIG,WORK,LWORK)
C*****************************************************************************
C
C     Write DHF-coefficients and eigenvalues from FORMATTED file
C
C     Written by T.Saue Feb 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
      LOGICAL TOBE
      CHARACTER TEXT*74
      DIMENSION CMO(NCMOTQ),EIG(NORBT),IBEIG(NORBT),WORK(*)
      DIMENSION IDIM(3,2)
#include "memint.h"
C
C     Get coefficients from DFCOEF
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) THEN
          WRITE(LUPRI,'(A)') 
     &      'PCMOUT: Unformatted file DFCOEF not found !'
          CALL QUIT('PCMOUT: Unformatted file DFCOEF not found !')
      ENDIF
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PCMOUT')        
      CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,14,
     &            WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Write coefficients to DFPCMO
C
      DO I = 1,NFSYM
        IDIM(1,I) = NPSH(I)
        IDIM(2,I) = NESH(I)
        IDIM(3,I) = NFBAS(I,0)
      ENDDO
      OPEN(LUCOEF,FILE='DFPCMO',STATUS='UNKNOWN',
     &     FORM='FORMATTED',ACCESS='SEQUENTIAL')
      TEXT(1:50) = TITLE
      CALL GTINFO(TEXT(51:74))
      WRITE(LUCOEF,'(A74)') TEXT
      WRITE(LUCOEF,'(I4,6I10)') NFSYM,((IDIM(I,J),I = 1,3),J=1,NFSYM)
      WRITE(LUCOEF,'(D24.16)') TOTERG
      WRITE(LUCOEF,'(6F22.16)') (CMO(I),I=1,NCMOTQ)
      WRITE(LUCOEF,'(6D22.12)') (EIG(I),I=1,NORBT)
      WRITE(LUCOEF,'(66I2)') (IBEIG(I),I=1,NORBT)
      CLOSE(LUCOEF,STATUS='KEEP')
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)')
     &    '* PCMOUT: Coefficients read from unformatted DFCOEF ',
     &    '          and written to formatted DFPCMO'
      CALL PRSYMB(LUPRI,'=',75,0)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrifck */
      SUBROUTINE WRIFCK(IUNIT,FOCK,N)
C*****************************************************************************
C
C     Write N AO - Fock matrices to unformatted file
C
C     Written by T.Saue Sept 1 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbdhf.h"
      CHARACTER TEXT*74
      DIMENSION FOCK(N2BBASXQ*N)
C
      TEXT(1:50) = TITLE
      CALL GTINFO(TEXT(51:74))
      REWIND IUNIT
      WRITE(IUNIT) TEXT
      WRITE(IUNIT) N2BBASXQ,N
      WRITE(IUNIT) FOCK
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reafck */
      SUBROUTINE REAFCK(IUNIT,FOCK,INFO,N)
C*****************************************************************************
C
C     Read AO Fock matrix from unformatted file
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbdhf.h"
      LOGICAL INFO
      CHARACTER TEXT*74, FILNAM*20
      DIMENSION FOCK(N2BBASXQ*N)
C
      REWIND IUNIT
      READ(IUNIT,END=10,ERR=20) TEXT
      READ(IUNIT) M2BBASXQ,M
      IF (N2BBASXQ .NE. M2BBASXQ .OR. N .GT. M) INFO = .TRUE.
      IF(INFO) THEN
        INQUIRE(UNIT=IUNIT, NAME=FILNAM)
        WRITE(LUPRI,'(/2A)')
     &   '* REAFCK: Fock matrix read from file ',FILNAM
        WRITE(LUPRI,'(2A)') '* Heading :',TEXT
      ENDIF
      IF (N2BBASXQ .NE. M2BBASXQ .OR. N .GT. M) THEN
         WRITE(LUPRI,'(/A,2(/A,2I10))') '*** FATAL ERROR ***',
     &      'Fock matrix size in this run and size on file: ',
     &      N2BBASXQ,M2BBASXQ,
     &      'Number of Fock matrices in this run and on file:',N,M
         CALL QUIT('REAFCK: inconsistent dimensions on file'//FILNAM)
      END IF
      READ (IUNIT,END=10,ERR=20) FOCK
C
      RETURN
C
C     Error messages
C
 10   CONTINUE
C     End of file
      CALL QUIT('REAFCK: End of file')
 20   CONTINUE
C     Error reading file
 30   CONTINUE
      CALL QUIT('REAFCK: ERROR reading file')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wridns */
      SUBROUTINE WRIDNS(IUNIT,DMAT)
C*****************************************************************************
C
C     Write SO - density matrix to unformatted file
C
C     Written by T.Saue Oct 12 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbdhf.h"
      CHARACTER TEXT*74
      DIMENSION DMAT(N2BBASXQ*NFMAT)
C
      TEXT(1:50) = TITLE
      CALL GTINFO(TEXT(51:74))
      REWIND IUNIT
      WRITE(IUNIT) TEXT
      WRITE(IUNIT) DMAT
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck readns */
      SUBROUTINE READNS(IUNIT,DMAT)
C*****************************************************************************
C
C     Read SO density matrix from unformatted file
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcbdhf.h"
      CHARACTER TEXT*74
      DIMENSION DMAT(N2BBASXQ*NFMAT)
C
      REWIND IUNIT
      READ(IUNIT,END=10,ERR=20) TEXT
      READ (IUNIT,END=30,ERR=40) DMAT
C
      RETURN
C
C     Error messages
C
 10   CONTINUE
      CALL QUIT('READNS: END reading TEXT')
 20   CONTINUE
      CALL QUIT('READNS: ERROR reading TEXT')
 30   CONTINUE
      CALL QUIT('READNS: END reading record 2')
 40   CONTINUE
      CALL QUIT('READNS: ERROR reading record 2')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck numlst */
      SUBROUTINE NUMLST(A,ILIST,NLIST,IRNGE1,IRNGE2,IFRP,IP,WORK,LWORK)
C*****************************************************************************
C
C     Original purpose was to read the numbers in the string A which has
C     the format : '3..6, 8, -27..2, 123, 14..15'.
C     
C     !! Number 0 (Zero) NOT counted or stored !!
C
C     Input: 
C           IF(IP.NE.0) : Give list of numbers in array
C                         ILIST(NLIST).
C           ELSE : Only IP as OUTPUT
C
C     Output:
C           IP : number of elements in the list A.
C                Duplicates only counted once.
C
C     IRNGE1..IRNGE2 is range of allowed values in A.
C     ITMPAR   : Temporary storage
C     If the string starts with 'all', all elements in the range enters.
C
C                           Jon K. Laerdahl 30.12.96
C     Minor polish T.Saue Jan 7 1997

C     Additional feature (LV, 15-5-2001) : Determine list based on energy or
C     other criteria. This is more user-friendly for standard cases.
C
C     This option is toggled by starting the string by energy, followed by the
C     lower treshold, upper treshold and minimal gap (to prevent cutting through
C     quasi-degenerate orbitals).
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dcbgen.h"
C
      CHARACTER A*(*)
      DIMENSION ILIST(NLIST),WORK(LWORK)
      LOGICAL TOBE
C
C     Memory initialization
C     =====================
C
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK
C
C     Distinguish between the two branches
C
      NA = LEN(A)
      IF (INDEX(A(1:NA),'energy') .EQ. 0) THEN
         CALL MEMGET('INTE',KITM,IRNGE2-IRNGE1+1,WORK,KFREE,LFREE)
         CALL NUMLS1(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,WORK(KITM))
      ELSE
         CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KIBG,NORBT,WORK,KFREE,LFREE)
C
C        Get eigenvalues from DFCOEF
C
         INQUIRE(FILE='DFCOEF',EXIST=TOBE)
         IF(.NOT.TOBE) THEN
C
C           DFCOEF does not (yet) exist, put IP = - 666, which makes 
C           NUMLS2 return with an empty list
C
            IP = - 666
            CALL NUMLS2(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,WORK(KEIG))
         ELSE
            CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PCMOUT')
            CALL REACMO(LUCOEF,'DFCOEF',DUM,WORK(KEIG),WORK(KIBG),DUM,4,
     &                  WORK,KFREE,LFREE)
            CLOSE(LUCOEF,STATUS='KEEP')
            IOFF = IORB(IFRP)+NPSH(IFRP)
            CALL NUMLS2(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,WORK(KEIG+IOFF))
         ENDIF
      ENDIF
C
      CALL MEMREL('NUMLST',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck numls1 */
      SUBROUTINE NUMLS1(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,ITMPAR)
C*****************************************************************************
C
C     Reads the numbers in the string A which has
C     the format : '3..6, 8, -27..2, 123, 14..15'.
C     
C     !! Number 0 (Zero) NOT counted or stored !!
C
C     Input: 
C           IP > 0: Give list of numbers in array
C                         ILIST(NLIST).
C           IP = 0: Only IP as OUTPUT
C           IP < 0: Give list of numbers in array
C                         ILIST(NLIST) - but the order is maintained,
C                         ie, 1-8,10,9  .NE.  1-10
C
C     Output:
C           IP : number of elements in the list A.
C                Duplicates only counted once.
C
C     IRNGE1..IRNGE2 is range of allowed values in A.
C     ITMPAR   : Temporary storage
C     If the string starts with 'all', all elements in the range enters.w
C
C                           Jon K. Laerdahl 30.12.96
C     Minor polish T.Saue Jan 7 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (MAXLB = 72)
C
      CHARACTER A*(*),B*(MAXLB)
      CHARACTER*1 TEGN
      DIMENSION ITMPAR(IRNGE1:IRNGE2),ILIST(NLIST)
C
C     Initialize
C
      NUM1 = 0
      NUM2 = 0
      INO  = 0
      NA   = LEN(A)
      IF(NA.GE.3) THEN
        IF( INDEX(A(1:NA),'all') .GT. 0) THEN
          IC = 0
          IF(IP.EQ.0) THEN
            DO 10 I=IRNGE1,IRNGE2
              IF (I.EQ.0) GOTO 10
              IC = IC + 1
C              IF (IC.GT.NLIST) GOTO 930
 10         CONTINUE
          ELSE
            DO 11 I=IRNGE1,IRNGE2
              IF (I.EQ.0) GOTO 11
              IC = IC + 1
              IF (IC.GT.NLIST) GOTO 930
              ILIST(IC) = I
 11         CONTINUE
          ENDIF
          IP = IC
          RETURN
        ENDIF
      ENDIF
      IPOINT = ICHAR('.')
      ISPACE = ICHAR(' ')
      ISIGN  = ICHAR('-')
      ICOMMA = ICHAR(',')
      INULL  = ICHAR('0')
      ININE  = ICHAR('9')
      IOOH   = ICHAR('o')
      CALL IZERO(ITMPAR,IRNGE2-IRNGE1+1)
      B = ' '
C
C     Take out blanks
C
      NB = 0
      DO 20 I = 1,NA
         TEGN = A(I:I)
         IF (TEGN.EQ.' ') GOTO 20
         NB = NB + 1
         IF(NB.GT.MAXLB) GOTO 940
         B(NB:NB) = TEGN
 20   CONTINUE
      I = 1
 100  CONTINUE
      NEG1 = 1
      NEG2 = 1
      INT = ICHAR(B(I:I))
      IF (INT.EQ.ISIGN) THEN
         NEG1 = -1
         I = I+1
         INT = ICHAR(B(I:I))
      ENDIF
      IF ((INT.GE.INULL).AND.(INT.LE.ININE) .OR.
     &    (INT.EQ.IOOH) ) THEN
         IF (INT .EQ. IOOH) THEN
            I = I + 1
            TEGN = B(I:I)
            IF (TEGN .NE. 'o') GOTO 950
            IF (NEG1 .EQ. -1) THEN
               NUM1 = -IRNGE1
            ELSE
               NUM1 = IRNGE2
            END IF
            I = I + 1
         ELSE
C           ... found a number
            NUM1 = INT-INULL
            CALL READNM(I,NUM1,B)
            NUM1 = NUM1*NEG1
         END IF
         TEGN = B(I:I)
         IF (TEGN.EQ.'.') THEN
C        ... checking for range
            I = I+1
            TEGN = B(I:I)
            IF (TEGN.NE.'.') THEN
               GOTO 900
            ELSE
C           ... find end of range
               I = I+1
               INT = ICHAR(B(I:I))
               IF (INT.EQ.ISIGN) THEN
                  IF (NEG1.GT.0) GOTO 900
                  NEG2 = -1
                  I = I+1
                  INT = ICHAR(B(I:I))
               ENDIF
               IF ((INT.GE.INULL).AND.(INT.LE.ININE)) THEN
C              ... read end of range
                  NUM2 = INT-INULL
                  CALL READNM(I,NUM2,B)
                  NUM2 = NUM2*NEG2
               ELSE IF (INT .EQ. IOOH) THEN
                  I = I + 1
                  TEGN = B(I:I)
                  IF (TEGN .NE. 'o') GOTO 950
                  IF (NEG2 .EQ. -1) THEN
                     NUM2 = -IRNGE1
                  ELSE
                     NUM2 = IRNGE2
                  END IF
                  I = I + 1
               ELSE 
                  GOTO 900
               ENDIF
            ENDIF
            IF (NUM1.GT.NUM2) GOTO 900
            IF ((NUM1.GE.IRNGE1).AND.(NUM2.LE.IRNGE2)) THEN
               DO 200 J=NUM1,NUM2
                  INO = INO + 1
                  ITMPAR(J) = INO
 200           CONTINUE
            ELSE
               WRITE(LUPRI,'(A,I5,A,I5)') 'Range:', NUM1,'..',NUM2
               GOTO 920
            ENDIF
         ELSEIF ((TEGN.EQ.',').OR.(TEGN.EQ.' ')) THEN
C           ... single number
            IF ((NUM1.GE.IRNGE1).AND.(NUM1.LE.IRNGE2)) THEN
               INO = INO + 1
               ITMPAR(NUM1) = INO
            ELSE
               WRITE(LUPRI,'(A,I5)') 'Value:', NUM1
               GOTO 920
            ENDIF
         ELSE 
            GOTO 910
         ENDIF
      ELSEIF ((INT.NE.ISPACE).AND.(INT.NE.IPOINT)) THEN
         GOTO 910
      ENDIF
      I = I+1
      IF (I.LE.NB) GOTO 100
C     
C     Count number of strings
C
      ICNT = 0
      DO 300 I=IRNGE1,IRNGE2
C
C     Do NOT count ZERO...
C
         IF (I.EQ.0) GOTO 300
         IF (ITMPAR(I).NE.0) ICNT = ICNT+1
 300  CONTINUE
C
      IF(IP.NE.0) THEN
         CALL IZERO(ILIST,NLIST)
         IC = 0
         DO 400 I=IRNGE1,IRNGE2
C
C     Do NOT store ZERO...
C
            IF (I.EQ.0) GOTO 400
            IF (ITMPAR(I).NE.0) THEN
               IC = IC + 1
               IF (IC.GT.NLIST) GOTO 930
               IF (IP.GT.0) THEN
                  ILIST(IC) = I
               ELSE
                  ILIST(ITMPAR(I)) = I
               END IF
            ENDIF
 400     CONTINUE
      ENDIF
      IP = ICNT
C
      RETURN
 900  CONTINUE
      WRITE(LUPRI,'(//A,I3//3A/9X,200A1)')
     &     'NUMLS1: Wrong syntax in input at character ',I,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS1: Wrong Syntax in Input')
 910  CONTINUE
      WRITE(LUPRI,'(/A,A1,A,I4,A,I3//3A/9X,200A1)')
     &   'NUMLS1: Wrong Character in Input: ',CHAR(INT),
     &   ' ( ASCII code',INT,') at character',I,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS1: Wrong Character in Input')
 920  CONTINUE
      WRITE(LUPRI,'(/A,2I6//3A/9X,200A1)')
     &     'Out of Range in NUMLS1. Range:',IRNGE1,IRNGE2,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS1: Out of Range!')
 930  CONTINUE
      WRITE(LUPRI,'(/A,I7//3A/9X,200A1)')
     &     'Out of Range in NUMLS1. NLIST:',NLIST,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS1: Out of Range in ITMPAR') 
 940  CONTINUE
      WRITE(LUPRI,'(/A//3A/9X,200A1)')
     &     'NUMLS1: Character string beyond 72.',
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS1: Character string beyond 72') 
 950  CONTINUE
      WRITE(LUPRI,'(/A//3A/9X,200A1)')
     &     'NUMLS1: oo must be two consecutive characters.',
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS1: oo must be two consecutive characters.')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck numls2 */
      SUBROUTINE NUMLS2(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,EIG)
C*****************************************************************************
C
C     Gives the list of oribitals that fulfill the energy criteria
C
C     Input: 
C           IP = 0: Only IP as OUTPUT
C           IP > 0: Give also list
C           IP < 0: Exit after printing the criteria 
C
C     Output:
C           IP : number of elements in the list A.
C
C     IRNGE1..IRNGE2 is range of allowed values in A.
C     EIG   : List of eigenvalues
C
C     L. Visscher, May 2001
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      CHARACTER A*(*)
      DIMENSION EIG(*)
      DIMENSION ILIST(NLIST)
      LOGICAL CHMIN,CHMAX
C     
C     Determine begin and end of the part of the string that contains
C     the energy values
C
      NA = LEN(A)
      MA = INDEX(A(1:NA),'energy') + 6
C     
      READ (A(MA:NA),*,END=900,ERR=900) ELOW,EHIGH,EGAP
C
      IF (IP.EQ.-666) THEN
         IP = 0
         RETURN
      ENDIF
C
      WRITE (LUPRI,1000) ELOW,EHIGH,EGAP
C
C     Check that the specified ELOW and EHIGH are separated by a 
C     large enough gap from the next lower/higher energies
C
   1  EMIN = 1.E16
      EMAX = - 1.E16
C
      DO I = IRNGE1, IRNGE2
        EI = EIG(I)
        IF (EI.GE.ELOW.AND.EI.LE.EHIGH) THEN
           IF (EMIN.GT.EI) EMIN = EI
           IF (EMAX.LT.EI) EMAX = EI
        ENDIF
      ENDDO
C
      CHMIN = .FALSE.
      CHMAX = .FALSE.
C
      DO I = IRNGE1, IRNGE2
        EI = EIG(I)
        IF (EI.LT.EMIN.AND.EI.GE.EMIN-EGAP) THEN
           IF (ELOW.GT.EI) ELOW = EI
           CHMIN = .TRUE.
        ENDIF
        IF (EI.GT.EMAX.AND.EI.LE.EMAX+EGAP) THEN
           IF (EHIGH.LT.EI) EHIGH = EI
           CHMAX = .TRUE.
        ENDIF
      ENDDO
C
      IF (CHMIN) WRITE (LUPRI,1001) 'minimum',ELOW
      IF (CHMAX) WRITE (LUPRI,1001) 'maximum',EHIGH
C
      IF (CHMIN.OR.CHMAX) GOTO 1
C
C     Determine the orbitals that fall within this energy range
C
      ICNT = 0
      DO I = IRNGE1, IRNGE2
        EI = EIG(I)
        IF (EI.GE.ELOW.AND.EI.LE.EHIGH) THEN
           ICNT = ICNT + 1
           IF (IP.NE.0) ILIST(ICNT) = I
        ENDIF
      ENDDO
C
      IP = ICNT
      RETURN
C
 900  CONTINUE
      WRITE(LUPRI,'(//A)')
     & 'NUMLS2: Error reading energy criteria, the string is ',A(1:NA)
      CALL QUIT('NUMLS2: Error in Input')
 1000 FORMAT (5X,'Energy selection of active orbitals :',
     & F10.2,' < Eps. < ',F8.2,' with a mininum gap of ',F8.4,' au.')
 1001 FORMAT (/' Changed ',A7,' energy for active orbital range to ',
     & F10.4,' because of near-degeneracies.')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck readnm */
      SUBROUTINE READNM(I,NUM,A)
      CHARACTER A*(*),TEGN*1
C
      INULL  = ICHAR('0')
      ININE  = ICHAR('9')
 10   I = I+1
      TEGN = A(I:I)
      INT = ICHAR(TEGN)
      IF ((INT.GE.INULL).AND.(INT.LE.ININE)) THEN
         NUM = NUM*10+(INT-INULL)
         GO TO 10
      ELSE
         RETURN
      ENDIF
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C   /* Deck wrtiar */
      SUBROUTINE WRTIAR(ILIST,IP)
#include "implicit.h"
#include "priunit.h"
      DIMENSION ILIST(IP)
C
      WRITE(LUPRI,'(10I5)') (ILIST(I), I=1,IP)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck EIGEXT */
      SUBROUTINE EIGEXT(EIG,IORBCL,ITYP,WORK,KFREE,LFREE)
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
C
      DIMENSION EIG(*),IORBCL(2,NORBT),WORK(*)
      KFRSAV = KFREE
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBEIG,NORBT ,WORK,KFREE,LFREE)
      LUBUF = 22
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','EIGEXT')
      CALL REACMO(LUBUF,'DFCOEF',WORK(KCMO),WORK(KEIG),WORK(KIBEIG),
     &            TOTERG,14,WORK,KFREE,LFREE)
      CLOSE(LUBUF,STATUS='KEEP')
      DO I = 1,NORBT
      IF(IORBCL(1,I).EQ.ITYP) THEN
        IND      = IORBCL(2,I) 
        EIG(IND) = WORK(KEIG+I-1)
      ENDIF
      ENDDO
      CALL MEMREL('EIGEXT',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Adacmo */
      SUBROUTINE ADACMO(NB,NZ,CF1,NO1,NP1,NE1,CF2,NO2,NP2,NE2)
C***********************************************************************
C
C     Adapt coefficients CF2 to current dimensions CF1
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION CF1(NB,NO1,NZ),CF2(NB,NO2,NZ)
C     Positronic part
      IV1 = 1
      IV2 = 1
      NC = MIN(NP1,NP2)
      NR = NP1-NP2
      IF(NR.GT.0) THEN
        NN = NB*NR
        DO IZ = 1,NZ
          CALL DZERO(CF1(1,IV1,IZ),NN)
        ENDDO
        IV1 = IV1 + NR
      ELSE
        IV2 = 1-NR
      ENDIF
      NN = NB*NC
      DO IZ = 1,NZ
        CALL DCOPY(NN,CF2(1,IV2,IZ),1,CF1(1,IV1,IZ),1)        
      ENDDO
C     Electronic part
      IV1 = NP1 + 1
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      NN = NB*NC
      DO IZ = 1,NZ
        CALL DCOPY(NN,CF2(1,IV2,IZ),1,CF1(1,IV1,IZ),1)        
      ENDDO
      IF(NR.GT.0) THEN
        IV1 = IV1 + NC
        NN  = NB*NR
        DO IZ = 1,NZ
          CALL DZERO(CF1(1,IV1,IZ),NN)
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Adaeig */
      SUBROUTINE ADAEIG(EIG1,NO1,NP1,NE1,EIG2,NO2,NP2,NE2)
C***********************************************************************
C
C     Adapt eigenvalues EIG2 to current dimensions EIG1
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION EIG1(NO1),EIG2(NO2)
C     Positronic part
      IV1 = 1
      IV2 = 1
      NC = MIN(NP1,NP2)
      NR = NP1-NP2
      IF(NR.GT.0) THEN
        CALL DZERO(EIG1(IV1),NR)
        IV1 = IV1 + NR
      ELSE
        IV2 = 1-NR
      ENDIF
      CALL DCOPY(NC,EIG2(IV2),1,EIG1(IV1),1)        
C     Electronic part
      IV1 = NP1 + 1
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      CALL DCOPY(NC,EIG2(IV2),1,EIG1(IV1),1)        
      IF(NR.GT.0) THEN
        IV1 = IV1 + NC
        CALL DZERO(EIG1(IV1),NR)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Adaibp */
      SUBROUTINE ADAIBP(IBO1,IBO2,IDIM)
C***********************************************************************
C
C     Pre-routine to adaibo
C     Written by T.Saue Feb 1 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION IBO1(*),IBO2(*),IDIM(3,2)
      I2 = 1
      DO IFRP = 1,NFSYM
        CALL ADAIBO(
     &         IBO1(IORB(IFRP)+1),NORB(IFRP),NPSH(IFRP),NESH(IFRP),
     &         IBO2(I2),IDIM(3,IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
          I2 = I2 + IDIM(2,IFRP)+IDIM(1,IFRP)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Adaibo */
      SUBROUTINE ADAIBO(IBO1,NO1,NP1,NE1,IBO2,NO2,NP2,NE2)
C***********************************************************************
C
C     Adapt boson irrep array IBO2 to current dimensions IBO1
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION IBO1(NO1),IBO2(NO2)
C     Positronic part
      IV1 = 1
      IV2 = 1
      NC = MIN(NP1,NP2)
      NR = NP1-NP2
      IF(NR.GT.0) THEN
        CALL IZERO(IBO1(IV1),NR)
        IV1 = IV1 + NR
      ELSE
        IV2 = 1-NR
      ENDIF
      CALL ICOPY(NC,IBO2(IV2),1,IBO1(IV1),1)        
C     Electronic part
      IV1 = NP1 + 1
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      CALL ICOPY(NC,IBO2(IV2),1,IBO1(IV1),1)        
      IF(NR.GT.0) THEN
        IV1 = IV1 + NC
        CALL IZERO(IBO1(IV1),NR)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Nucbas */
      SUBROUTINE NUCBAS(NUCORB,AA,IPRINT)
C***********************************************************************
C
C     Extract basis information for all centers
C
C     Written by T.Saue March 12 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "dummy.h"
C
#include "nuclei.h"
#include "ccom.h"
#include "shells.h"
#include "primit.h"
      CHARACTER SPDCAR*1
      DIMENSION NUCORB(NHTYP,2,NUCIND),AA(2,NHTYP,2,NUCIND)
C
C     Initialize
C
      NDIM = 2*NHTYP*NUCIND
      CALL IZERO(NUCORB,NDIM)
C
      JCENT = 0
      JPRIM = -1
      JC = -1
      DO ISHELL = 1,KMAX
        ICENT = NCENT(ISHELL)
        IF(ICENT.NE.JCENT) THEN
          JCENT = ICENT
          JLVAL = 0
        ENDIF
        IC = LCLASS(ISHELL)
        IF(IC.NE.JC) THEN
          JC    = IC
          JLVAL = 0
        ENDIF
        ILVAL = NHKT(ISHELL)
        IF(ILVAL.NE.JLVAL) THEN
          JLVAL = ILVAL
          NUCORB(ILVAL,IC,ICENT) = 0
          AA(1,ILVAL,IC,ICENT)=D0
          AA(2,ILVAL,IC,ICENT)=DUMMY
        ENDIF
        NUCORB(ILVAL,IC,ICENT)=NUCORB(ILVAL,IC,ICENT)+1
        IPRIM = JSTRT(ISHELL)
        IF(IPRIM.NE.JPRIM) THEN
          JPRIM = IPRIM
          NPRIM = NUCO(ISHELL)
          DO IEXP = 1,NPRIM
            A=PRIEXP(IPRIM+IEXP)
            AA(1,ILVAL,IC,ICENT)=MAX(AA(1,ILVAL,IC,ICENT),A)
            AA(2,ILVAL,IC,ICENT)=MIN(AA(2,ILVAL,IC,ICENT),A)
          ENDDO
        ENDIF
      ENDDO
C
      IF(IPRINT.GE.2) THEN
         CALL HEADER('NUCORB:Basis set information:',-1)
        DO I = 1,NUCIND
          WRITE(LUPRI,'(/A,A4,A/)') '*** Center: ',NAMN(I),' ***'
          WRITE(LUPRI,'(2X,A)') '* Large components:'
          IC = 1
          DO LL = 1,NHTYP
          IF(NUCORB(LL,IC,I).GT.0) THEN
            L=LL-1
            WRITE(LUPRI,'(3X,A1,A,I6,2(3X,A,E12.5))')
     &      SPDCAR(L),'-orbitals:',NUCORB(LL,IC,I),
     &      'Alpha_H :',AA(1,LL,IC,I),
     &      'Alpha_L :',AA(2,LL,IC,I)
          ENDIF
          ENDDO
          WRITE(LUPRI,'(2X,A)') '* Small components:'
          IC = 2
          DO LL = 1,NHTYP
          IF(NUCORB(LL,IC,I).GT.0) THEN
            L=LL-1
            WRITE(LUPRI,'(3X,A1,A,I6,2(3X,A,E12.5))')
     &      SPDCAR(L),'-orbitals:',NUCORB(LL,IC,I),
     &      'Alpha_H: ',AA(1,LL,IC,I),
     &      'Alpha_L: ',AA(2,LL,IC,I)
          ENDIF
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck acmout */
      SUBROUTINE ACMOUT(WORK,LWORK)
C*****************************************************************************
C
C     Read coefficients from unformatted DFCOEF, 
C     take out symmetry and write to unformatted DFACMO 
C
C     Written by T.Saue May 26 2003
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
C
      DIMENSION WORK(LWORK)
      CALL QENTER('ACMOUT')
#include "memint.h"
C
C     Memory allocation
      CALL MEMGET('REAL',KCMO1,N2BBASX*4,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO2,N2BBASX*4,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG1,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG2,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBG1,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBG2,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINUC,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,NTBAS(0),WORK,KFREE,LFREE)
C
      CALL ACMOU1(WORK(KCMO1),WORK(KEIG1),WORK(KIBG1),
     &            WORK(KCMO2),WORK(KEIG2),WORK(KIBG2),
     &            WORK(KINUC),WORK(KIBUF),WORK(KFREE),LFREE)
C
C     Memory deallocation
      CALL MEMREL('ACMOUT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('ACMOUT')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck acmou1 */
      SUBROUTINE ACMOU1(CMO1,EIG1,IBEIG1,CMO2,EIG2,IBEIG2,INUC,
     &                  IBUF,WORK,LWORK)
C*****************************************************************************
C
C     Read coefficients from unformatted DFCOEF, 
C     take out symmetry and write to unformatted DFACMO 
C
C     Written by T.Saue May 26 2003
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER(DM1 = -1.0D0)
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcblab.h"
#include "dgroup.h"
#include "nuclei.h"
      LOGICAL TOBE
      CHARACTER TEXT*74
      DIMENSION CMO1(*),EIG1(NORBT),IBEIG1(NORBT),
     &          CMO2(*),EIG2(NORBT),IBEIG2(NORBT),
     &          INUC(*),IBUF(*),WORK(LWORK)
      DIMENSION IDIM(3)
#include "memint.h"
C
C     Find correct ordering of nuclei
C
      CALL LSQMA1(INUC)
      I1 = 0
      I2 = 0
      DO IC = 1,2
        DO ICENT = 1,NUCDEP
          DO IB = 1,NTBAS(IC)
          IF(INUC(I1+IB).EQ.ICENT) THEN
            I2 = I2 + 1
            IBUF(I2) = I1+IB
          ENDIF
          ENDDO
        ENDDO
        I1 = I1 + NTBAS(IC)
      ENDDO
C
C     Get coefficients from DFCOEF --> stored in CMO1
C     ===============================================
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) THEN
          WRITE(LUPRI,'(A)') 
     &      'ACMOUT: Unformatted file DFCOEF not found !'
          CALL QUIT('ACMOUT: Unformatted file DFCOEF not found !')
      ENDIF
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PCMOUT')        
      CALL REACMO(LUCOEF,'DFCOEF',CMO1,EIG1,IBEIG1,TOTERG,14,
     &            WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Transform to AO-basis --> stored in CMO2
C     ========================================
C
      KC = 1
      DO IFRP = 1,NFSYM
      IF(NORB(IFRP).GT.0) THEN
        IF(NBSYM.EQ.1) THEN
          CALL CFSORT(CMO2(KC),CMO1(ICMOQ(IFRP)+1),NORB(IFRP),
     &                IFRP,NFBAS(IFRP,0))
        ELSE
          NCBF2 = NFBAS(IFRP,0)*NORB(IFRP)*4
          CALL MEMGET('REAL',KCBF2,NCBF2,WORK,KFREE,LFREE)
          CALL CFSORT(WORK(KCBF2),CMO1(ICMOQ(IFRP)+1),NORB(IFRP),
     &                IFRP,NFBAS(IFRP,0))
          CALL WTSOAO(WORK(KCBF2),CMO2(KC),IFRP,
     &                    NFBAS(IFRP,0),NTBAS(0),NORB(IFRP)*4,0)
          CALL MEMREL('ACMOU1.AO',WORK,KCBF2,KCBF2,KFREE,LFREE)
        ENDIF
        KC = KC + NTBAS(0)*NORB(IFRP)*4
      ENDIF
      ENDDO
C
C     Get correct row ordering (NBAS) --> stored in CMO1     
C     ==================================================
C
      I1 = 1
      I2 = 1
      DO IFRP = 1,NFSYM
        CALL ACMOU2(CMO1(I1),NORBT,CMO2(I2),NORB(IFRP),
     &              NTBAS(0),IBUF)
        I1 = I1 + NTBAS(0)*NORB(IFRP)
        I2 = I2 + NTBAS(0)*NORB(IFRP)*4
      ENDDO
C
C     In CFSORT the BR matrix is scaled by -1. Scale back
C     ===================================================
C
      I1 = 1 + NTBAS(0)*NORBT*2
      CALL DSCAL(NTBAS(0)*NORBT,DM1,CMO1(I1),1)
C
C     Get correct column ordering
C     ===========================
C 
      CALL INDEXX(NORBT,EIG1,IBUF)
C     Sort eigenvalues and boson irrep info
      DO J = 1,NORBT
        EIG2(J) = EIG1(IBUF(J))
        IBEIG2(J) = IBEIG1(IBUF(J))
        CALL ACMOU3(CMO1,CMO2,IBUF(J),J,NTBAS(0),NORBT)
      ENDDO
C
C     Write coefficients etc. to file
C     ===============================
C
      CALL OPNFIL(LUCOEF,'DFACMO','NEW','ACMOUT')
      DO I = 1,2
        IDIM(I) = 0
      ENDDO
      DO IFRP = 1,NFSYM
        IDIM(1) = IDIM(1) + NPSH(IFRP)  
        IDIM(1) = IDIM(1) + NESH(IFRP)  
      ENDDO
      IDIM(3) = NTBAS(0)
      TEXT(1:50) = TITLE
      CALL GTINFO(TEXT(51:74))
      WRITE(LUCOEF) TEXT,1,(IDIM(I),I = 1,3),TOTERG
      NDIM = NTBAS(0)*NORBT*4
      CALL WRITT(LUCOEF,NDIM,CMO2)
      CALL WRITT(LUCOEF,NORBT,EIG2)
      CALL WRITI(LUCOEF,NORBT,IBEIG2)
      CLOSE(LUCOEF,STATUS='KEEP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck acmou2 */
      SUBROUTINE ACMOU2(CMO1,NMO1,CMO2,NMO2,NBAS,IBUF)
C***********************************************************************
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION CMO1(NBAS,NMO1,4),CMO2(NBAS,NMO2,4),IBUF(*)
C
      DO IZ = 1,4
        DO J = 1,NMO2
          DO I = 1,NBAS
            CMO1(I,J,IZ) = CMO2(IBUF(I),J,IZ)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck acmou3 */
      SUBROUTINE ACMOU3(CMO1,CMO2,J1,J2,NBAS,NORB)
C***********************************************************************
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION CMO1(NBAS,NORB,4),CMO2(NBAS,NORB,4)
C
      DO IZ = 1,4
        CALL DCOPY(NBAS,CMO1(1,J1,IZ),1,CMO2(1,J2,IZ),1)        
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Reaprj */
      SUBROUTINE REAPRJ(LUPRJ,NP,NR,KRVEC,KBVEC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Read DFPRJC, file from projection analysis
C     Written by T. Saue, Oct 21 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbprj.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION NR(NREFS,NFSYM),NP(2,0:2),WORK(*)
      NCDIM = 0
      NBDIM = 0
      READ(LUPRJ) MREFS
      IF(NREFS.GT.MREFS) CALL QUIT('REAPRJ:MREFS.NE.NREFS !')
      DO IFRP = 1,NFSYM
C       Dimension....nref.......nmol......
        READ(LUPRJ) NP(IFRP,1),NP(IFRP,2),(NR(I,IFRP),I=1,NREFS)
        NP(IFRP,0) = NP(IFRP,1)+NP(IFRP,2)
        NCDIM = NCDIM + NFBAS(IFRP,0)*NP(IFRP,0)*NZ
        NBDIM = NBDIM + NP(IFRP,0)*NP(IFRP,2)*NZ
      ENDDO
      CALL MEMGET('REAL',KRVEC,NCDIM,WORK,KFREE,LFREE)       
      CALL MEMGET('REAL',KBVEC,NBDIM,WORK,KFREE,LFREE)       
      CALL READT(LUPRJ,NCDIM,WORK(KRVEC))
      CALL READT(LUPRJ,NBDIM,WORK(KBVEC))
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck numls3 */
      SUBROUTINE NUMLS3(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,NZERO)
C*****************************************************************************
C
C     Reads the numbers in the string A which has
C     the format : '3..6, 8, -27..2, 123, 14..15'.
C     
C     !! Number 0 (Zero) IS counted AND stored, and returned as number
C
C     Input: 
C           IP.EQ.0: Only IP as OUTPUT
C           IP.NE.0: Give list of numbers in array ILIST(NLIST).
C
C     Output:
C           IP : number of elements in the list A.
C
C     IRNGE1..IRNGE2 is range of allowed values in A.
C
C     Based on the routine NUMLS1 written by Jon K. Laerdahl 30.12.96
C     Written by T.Saue 2003
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
C
      PARAMETER (MAXLB = 72)
C
      CHARACTER A*(*),B*(MAXLB)
      CHARACTER*1 TEGN
      DIMENSION ILIST(NLIST)
C
C     Initialize
C
      NZERO = 0
      NUM1  = 0
      NUM2  = 0
      INO   = 0
      NA = LEN(A)
      IF(NA.GE.3) THEN
        IF( INDEX(A(1:NA),'all') .GT. 0) THEN
          IC = 0
          IF(IP.EQ.0) THEN
            DO 10 I=IRNGE1,IRNGE2
              IF (I.EQ.0) GOTO 10
              IC = IC + 1
 10         CONTINUE
          ELSE
            DO 11 I=IRNGE1,IRNGE2
              IF (I.EQ.0) GOTO 11
              IC = IC + 1
              IF (IC.GT.NLIST) GOTO 930
              ILIST(IC) = I
 11         CONTINUE
          ENDIF
          IP = IC
          RETURN
        ENDIF
      ENDIF
      IPOINT = ICHAR('.')
      ISPACE = ICHAR(' ')
      ISIGN  = ICHAR('-')
      ICOMMA = ICHAR(',')
      INULL  = ICHAR('0')
      ININE  = ICHAR('9')
      IOOH   = ICHAR('o')
      B = ' '
C
C     Take out blanks
C
      NB = 0
      DO 20 I = 1,NA
         TEGN = A(I:I)
         IF (TEGN.EQ.' ') GOTO 20
         NB = NB + 1
         IF(NB.GT.MAXLB) GOTO 940
         B(NB:NB) = TEGN
 20   CONTINUE
      I = 1
 100  CONTINUE
      NEG1 = 1
      NEG2 = 1
      INT = ICHAR(B(I:I))
      IF (INT.EQ.ISIGN) THEN
         NEG1 = -1
         I = I+1
         INT = ICHAR(B(I:I))
      ENDIF
      IF ((INT.GE.INULL).AND.(INT.LE.ININE) .OR.
     &    (INT.EQ.IOOH) ) THEN
         IF (INT .EQ. IOOH) THEN
C           ... "infinity"
            I = I + 1
            TEGN = B(I:I)
            IF (TEGN .NE. 'o') GOTO 950
            IF (NEG1 .EQ. -1) THEN
               NUM1 = -IRNGE1
            ELSE
               NUM1 = IRNGE2
            END IF
            I = I + 1
         ELSE
C           ... found a number
            NUM1 = INT-INULL
            CALL READNM(I,NUM1,B)
            NUM1 = NUM1*NEG1
         END IF
         TEGN = B(I:I)
         IF (TEGN.EQ.'.') THEN
C        ... checking for range
            I = I+1
            TEGN = B(I:I)
            IF (TEGN.NE.'.') THEN
               GOTO 900
            ELSE
C           ... find end of range
               I = I+1
               INT = ICHAR(B(I:I))
               IF (INT.EQ.ISIGN) THEN
                  IF (NEG1.GT.0) GOTO 900
                  NEG2 = -1
                  I = I+1
                  INT = ICHAR(B(I:I))
               ENDIF
               IF ((INT.GE.INULL).AND.(INT.LE.ININE)) THEN
C              ... read end of range
                  NUM2 = INT-INULL
                  CALL READNM(I,NUM2,B)
                  NUM2 = NUM2*NEG2
               ELSE IF (INT .EQ. IOOH) THEN
                  I = I + 1
                  TEGN = B(I:I)
                  IF (TEGN .NE. 'o') GOTO 950
                  IF (NEG2 .EQ. -1) THEN
                     NUM2 = -IRNGE1
                  ELSE
                     NUM2 = IRNGE2
                  END IF
                  I = I + 1
               ELSE 
                  GOTO 900
               ENDIF
            ENDIF
            IF (NUM1.GT.NUM2) GOTO 900
            IF ((NUM1.GE.IRNGE1).AND.(NUM2.LE.IRNGE2)) THEN
               DO 200 J=NUM1,NUM2
                  INO = INO + 1
                  IF(J.EQ.0) NZERO=NZERO+1
                  IF(IP.NE.0) ILIST(INO)=J
 200           CONTINUE
            ELSE
               WRITE(LUPRI,'(A,I5,A,I5)') 'Range:', NUM1,'..',NUM2
               GOTO 920
            ENDIF
         ELSEIF ((TEGN.EQ.',').OR.(TEGN.EQ.' ')) THEN
C           ... single number
            IF ((NUM1.GE.IRNGE1).AND.(NUM1.LE.IRNGE2)) THEN
               INO = INO + 1
               IF(NUM1.EQ.0) NZERO=NZERO+1
               IF(IP.NE.0) ILIST(INO)=NUM1
            ELSE
               WRITE(LUPRI,'(A,I5)') 'Value:', NUM1
               GOTO 920
            ENDIF
         ELSE 
            GOTO 910
         ENDIF
      ELSEIF ((INT.NE.ISPACE).AND.(INT.NE.IPOINT)) THEN
         GOTO 910
      ENDIF
      I = I+1
      IF (I.LE.NB) GOTO 100
      IP = INO
      RETURN
 900  CONTINUE
      WRITE(LUPRI,'(//A,I3//3A/9X,200A1)')
     &     'NUMLS3: Wrong syntax in input at character ',I,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS3: Wrong Syntax in Input')
 910  CONTINUE
      WRITE(LUPRI,'(/A,A1,A,I4,A,I3//3A/9X,200A1)')
     &   'NUMLS3: Wrong Character in Input: ',CHAR(INT),
     &   ' ( ASCII code',INT,') at character',I,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS3: Wrong Character in Input')
 920  CONTINUE
      WRITE(LUPRI,'(/A,2I6//3A/9X,200A1)')
     &     'Out of Range in NUMLS3. Range:',IRNGE1,IRNGE2,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS3: Out of Range!')
 930  CONTINUE
      WRITE(LUPRI,'(/A,I7//3A/9X,200A1)')
     &     'Out of Range in NUMLS3. NLIST:',NLIST,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS3: Out of Range in ITMPAR') 
 940  CONTINUE
      WRITE(LUPRI,'(/A//3A/9X,200A1)')
     &     'NUMLS3: Character string beyond 72.',
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS3: Character string beyond 72') 
 950  CONTINUE
      WRITE(LUPRI,'(/A//3A/9X,200A1)')
     &     'NUMLS3: oo must be two consecutive characters.',
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS3: oo must be two consecutive characters.')
      END
