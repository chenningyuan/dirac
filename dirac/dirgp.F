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
C FILE    : dirgp.F
C SOURCE  : $Source$
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
#include <implicit.h>
#include <priunit.h>
#include <dcbgen.h>
#include <dcbbas.h>
#include <dcborb.h>
C
      DIMENSION WORK(LWORK)
      CALL QENTER('PCMOIN')
#include <memint.h>
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
#include <implicit.h>
#include <priunit.h>
#include <dcbgen.h>
#include <dcbbas.h>
#include <dcborb.h>
C
      DIMENSION WORK(LWORK)
      CALL QENTER('PCMOUT')
#include <memint.h>
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
#include <implicit.h>
#include <priunit.h>
#include <dgroup.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dcbgen.h>
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reacmo */
      SUBROUTINE REACMO(IUNIT,CMOFIL,CMO,EIG,IBEIG,TOTERG,IOPT,
     &                   WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Read DHF-coefficients and eigenvalues from unformatted file
C     IOPT:
C       0 - read only title
C       1 - read coefficients and eigenvalues
C       2 - read coefficients, and eigenvalues and irrep identification
C       3 - read coefficients and eigenvalues - give restart info
C       4 - read coefficients only - give info
C       5 - read only eigenvalues
C       6 - read coefficients only
C
C     IDIM(1,J) - number of positronic solutions
C     IDIM(2,J) - number of electronic solutions
C     IDIM(3,J) - number of AO-basis functions
C
C     Written by T.Saue Sept 1 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
      LOGICAL INFO,LBIT,TEST
      CHARACTER TEXT*74,CMOFIL*6,FMT*6,MXFORM*6
      DIMENSION CMO(NCMOTQ),EIG(NORBT),IBEIG(NORBT),WORK(*)
      DIMENSION IDIM(3,2)
C
      INFO = IOPT.EQ.3.OR.IOPT.EQ.0.OR.IOPT.EQ.4
      REWIND IUNIT
C
      READ (IUNIT,END=10,ERR=20) TEXT,NSYM,
     &       ((IDIM(I,J),I = 1,3),J=1,NSYM),TOTERG
C
      IF(INFO) THEN
        FMT = MXFORM(TOTERG,20)
        WRITE(LUPRI,'(/A,I5,3X,A,'//FMT//')')
     &   'REACMO: Coefficients read from unit',IUNIT,
     &   '- Total energy: ',TOTERG
        WRITE(LUPRI,'(2A)') '* Heading :',TEXT
      ENDIF
      IF(IOPT.EQ.0) RETURN
C
C     Compatibility test
      IF(NSYM.NE.NFSYM) 
     &   CALL QUIT('REACMO: Incompatible number of fermion ircops')
      TEST = .FALSE.
      NC = 0
      NO = 0
      DO IFRP = 1,NFSYM
        IF(IDIM(3,IFRP).NE.NFBAS(IFRP,0)) 
     &    CALL QUIT('REACMO: Incompatible number of basis functions')
        IDIM(3,IFRP) = IDIM(1,IFRP)+IDIM(2,IFRP)
        IF(IDIM(3,IFRP).NE.NORB(IFRP)) TEST=.TRUE.
        NC = NC + IDIM(3,IFRP)*NFBAS(IFRP,0)*NZ
        NO = NO + IDIM(3,IFRP)
      ENDDO
C
C     Read coefficients
      IF(IOPT.EQ.5) THEN
        READ(IUNIT)
      ELSEIF(TEST) THEN
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
      IF(IOPT.EQ.4 .OR. IOPT.EQ.6) RETURN
C     Read eigenvalues
      IF(TEST) THEN
        CALL MEMGET('REAL',KO,NO,WORK,KFREE,LFREE)
        CALL READT(IUNIT,NO,WORK(KO))
        JO = KO 
        DO IFRP = 1,NFSYM
          CALL ADAEIG(
     &         EIG(IORB(IFRP)+1),NORB(IFRP),NPSH(IFRP),NESH(IFRP),
     &         WORK(JO),IDIM(3,IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
          JC = JC + IDIM(3,IFRP)
        ENDDO
        CALL MEMREL('REACMO.eig',WORK,KO,KO,KFREE,LFREE)
      ELSE
        READ (IUNIT,END=50,ERR=60) EIG
      ENDIF
      IF(IOPT.NE.2) RETURN
C     Read single group irrep (spinfree and Levy-Leblond calcs)
      IF(TEST) THEN
        CALL MEMGET('INTE',KO,NO,WORK,KFREE,LFREE)
        CALL READI(IUNIT,NO,WORK(KO))
        CALL ADAIBP(IBEIG,WORK(KO),IDIM)
        CALL MEMREL('REACMO.ibos',WORK,KO,KO,KFREE,LFREE)
      ELSE
        READ (IUNIT,END=70,ERR=80) IBEIG
      ENDIF
C
      RETURN
C
C     Error messages
C
 10   CONTINUE
      CALL QUIT('REACMO: END reading TEXT')
 20   CONTINUE
      CALL QUIT('REACMO: ERROR reading TEXT')
 30   CONTINUE
      CALL QUIT('REACMO: END reading coefficients')
 40   CONTINUE
      CALL QUIT('REACMO: ERROR reading coefficients')
 50   CONTINUE
      CALL QUIT('REACMO: END reading eigenvalues')
 60   CONTINUE
c     CALL QUIT('REACMO: ERROR reading eigenvalues')
      WRITE(LUPRI,'(A)') 'WARNING !!! No eigenvalues read'
      RETURN
 70   CONTINUE
      CALL QUIT('REACMO: END reading irrep identification')
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
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcbdhf.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
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
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcbdhf.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
      LOGICAL TOBE
      CHARACTER TEXT*74
      DIMENSION CMO(NCMOTQ),EIG(NORBT),IBEIG(NORBT),WORK(*)
      DIMENSION IDIM(3,2)
#include <memint.h>
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
      CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,2,
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
C     Write AO - Fock matrix to unformatted file
C
C     Written by T.Saue Sept 1 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbgen.h>
#include <dcbdhf.h>
      CHARACTER TEXT*74
      DIMENSION FOCK(N2BBASXQ*N)
C
      TEXT(1:50) = TITLE
      CALL GTINFO(TEXT(51:74))
      REWIND IUNIT
      WRITE(IUNIT) TEXT
      WRITE(IUNIT) FOCK
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reafck */
      SUBROUTINE REAFCK(IUNIT,FOCK,RESTART,N)
C*****************************************************************************
C
C     Read AO Fock matrix from unformatted file
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbdhf.h>
      LOGICAL RESTART
      CHARACTER TEXT*74
      DIMENSION FOCK(N2BBASXQ*N)
C
      REWIND IUNIT
      READ(IUNIT,END=10,ERR=20) TEXT
      IF(RESTART) THEN
        CALL HEADER('REAFCK: Restart from two-electron Fock matrix',-1)
        WRITE(LUPRI,'(A,I5)')
     &   '* Fock matrix read from unit ',IUNIT
        WRITE(LUPRI,'(2A)') '* Heading :',TEXT
      ENDIF
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
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dcbgen.h>
#include <dcbdhf.h>
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
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dcbdhf.h>
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
      SUBROUTINE NUMLST(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,WORK,LWORK)
C*****************************************************************************
C
C     Reads the numbers in the string A which has
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
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      CHARACTER A*(*)
      DIMENSION ILIST(NLIST),WORK(LWORK)
C
C#include <memint.h>
C
C     Memory initialization
C     =====================
C
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK
C
      CALL MEMGET('INTE',KITM,IRNGE2-IRNGE1+1,WORK,KFREE,LFREE)
      CALL NUMLS1(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,WORK(KITM))
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
#include <implicit.h>
#include <priunit.h>
      CHARACTER A*(*),B*72
      CHARACTER*1 TEGN
      DIMENSION ITMPAR(IRNGE1:IRNGE2),ILIST(NLIST)
C
C     Initialize
C
      NUM1 = 0
      NUM2 = 0
      INO  = 0
      NA   = LEN(A)
      IF(NA.GT.72) GOTO 940
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
            WRITE(LUPRI,'(A,A1,A1,I3,A1)')
     &        'NUMLST: Wrong Character in Input: ',CHAR(INT),'(',INT,')'
            GOTO 910
         ENDIF
      ELSEIF ((INT.NE.ISPACE).AND.(INT.NE.IPOINT)) THEN
         WRITE(LUPRI,'(A,A1,A1,I3,A1)')
     &      'NUMLST: Wrong Character in Input: ',CHAR(INT),'(',INT,')'
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
      CALL QUIT('NUMLST: Wrong Syntax in Input')
 910  CONTINUE
      CALL QUIT('NUMLST: Wrong Characters in Input')
 920  CONTINUE
      WRITE(LUPRI,'(A,2I6)') 'Out of Range in NUMLST. Range:',
     &     IRNGE1,IRNGE2
      CALL QUIT('NUMLST: Out of Range!')
 930  CONTINUE
      WRITE(LUPRI,'(A,I7)') 'Out of Range in NUMLST. NLIST:',
     &     NLIST
      CALL QUIT('NUMLST: Out of Range in ITMPAR') 
 940  CONTINUE
      WRITE(LUPRI,'(A)') 'Character string beyond 72'
      CALL QUIT('NUMLST: Character string beyond 72') 
 950  CONTINUE
      CALL QUIT('NUMLST: oo must be two consecutive characters.')
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
#include <implicit.h>
#include <priunit.h>
      DIMENSION ILIST(IP)
C
      WRITE(LUPRI,'(10I5)') (ILIST(I), I=1,IP)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION POLVAL(NORDER,B,XVAL)
#include <implicit.h>
      DIMENSION B(NORDER+1)
      POLVAL = B(1)
      XBUF = 1
      DO 10 I = 2,(NORDER+1)
        XBUF = XBUF*XVAL
        POLVAL = POLVAL + XBUF*B(I)
   10 CONTINUE     
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck EIGEXT */
      SUBROUTINE EIGEXT(EIG,IORBCL,ITYP,WORK,KFREE,LFREE)
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
C
      DIMENSION EIG(*),IORBCL(2,NORBT),WORK(*)
      KFRSAV = KFREE
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBEIG,NORBT ,WORK,KFREE,LFREE)
      IOPT  = 2
      LUBUF = 22
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','EIGEXT')
      CALL REACMO(LUBUF,'DFCOEF',WORK(KCMO),WORK(KEIG),WORK(KIBEIG),
     &            TOTERG,IOPT,WORK,KFREE,LFREE)
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
#include <implicit.h>
#include <priunit.h>
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
      IV1 = IV1 + NC
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      NN = NB*NC
      DO IZ = 1,NZ
        CALL DCOPY(NN,CF2(1,IV2,IZ),1,CF1(1,IV1,IZ),1)        
      ENDDO
      IF(NR.GT.0) THEN
        NN = NB*NR
        DO IZ = 1,NZ
          CALL DZERO(CF1(1,IV1,IZ),NN)
        ENDDO
        IV1 = IV1 + NR
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
#include <implicit.h>
#include <priunit.h>
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
      IV1 = IV1 + NC
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      CALL DCOPY(NC,EIG2(IV2),1,EIG1(IV1),1)        
      IF(NR.GT.0) THEN
        CALL DZERO(EIG1(IV1),NR)
        IV1 = IV1 + NR
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
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dgroup.h>
      DIMENSION IBO1(*),IBO2(*),IDIM(3,2)
      I2 = 1
      DO IFRP = 1,NFSYM
        CALL ADAIBO(
     &         IBO1(IORB(IFRP)+1),NORB(IFRP),NPSH(IFRP),NESH(IFRP),
     &         IBO2(I2),IDIM(3,IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
          I2 = I2 + IDIM(3,IFRP)
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
#include <implicit.h>
#include <priunit.h>
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
      IV1 = IV1 + NC
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      CALL ICOPY(NC,IBO2(IV2),1,IBO1(IV1),1)        
      IF(NR.GT.0) THEN
        CALL IZERO(IBO1(IV1),NR)
        IV1 = IV1 + NR
      ENDIF
C
      RETURN
      END
