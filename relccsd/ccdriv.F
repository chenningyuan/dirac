C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : ccdriv.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCHEADER
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Write header for Relativistic Coupled Cluster Program
C     Keep track of all the modifications made in the program.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER DATEX*10,TIMEX*8
C
C---------------Executable code--------------------------------------
C
      WRITE (IW,1000) "2.0 ","October 1997"
      CALL DAYTIME (DATEX,TIMEX)
      WRITE (IW,1001) DATEX,TIMEX
C
 1000 FORMAT (//" Relativistic Coupled Cluster program RELCCSD"
     &/" Version number ",A
     &/" Version date ",A
     &//" Written by :",/" Luuk Visscher"/
     &/" NASA Ames Research Center     "
     &/" Rijks Universiteit Groningen  "
     &/" Odense Universitet            ")
 1001 FORMAT (//" Today is :",T15,A10/" The time is :",T17,A8)
C
      RETURN
C---------------Record of modifications made in the code-------------
C     11-12-97 : MP2 gradients works, CCSD Lambda equations written
C--------------------------------------------------------------------
C     October 1997 : major restructuring
C     Split the code in energy, gradient and second derivative modules
C     Final version should contain energies and 1st and 2nd order
C     properties at all the available levels of theory
C     Improved memory management (thanks to H-J Aa Jensens routines)
C--------------------------------------------------------------------
C     Aug  6 1997 :  
C     Implemented excitation energy module. 
C--------------------------------------------------------------------
C     Apr 16 1997 :  
C     Interfaced with DIRAC, reorganized code. 
C--------------------------------------------------------------------
C     Jul 10 1996 :  
C     Added RPA module, changed output format, 
C     added NOCCSD (MP2) option
C     Note that the RPA is not very optimized : needs complex 
C     arithmetics for completely imaginary operators, does not
C     consider symmetry of the operators
C--------------------------------------------------------------------
C     Oct 5 1995 :  
C     Corrected errors in complex arithmetic, implemented -T correction
C     Frozen option was also implemented somewhat earlier
C--------------------------------------------------------------------
C     Aug 13 1995 :  
C     Implemented triples correction in symmetry-adapted version
C--------------------------------------------------------------------
C     Jun 29 1995 :  
C     Dynamic memory allocation & incore files for CRAY
C--------------------------------------------------------------------
C     Jan 10 1995 : 
C     Real/complex arithmetic possible
C     Symmetry packed version
C--------------------------------------------------------------------
C     Oct 19 1994 : 
C     Writes integrals out on disk
C     Otherwise identical to previous version
C--------------------------------------------------------------------
C     Sep 22 1994 : First numbered version of program 
C     Keeps all integrals in core.
C     Does not use any symmetry.
C     But is vectorized by using matrix multiplies where ever possible
C--------------------------------------------------------------------
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCINPT (IN,IW,MIOBUF)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Reads input for relativistic Coupled Cluster program
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "inpt.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL REWINP
C
      NAMELIST/RELCCSD/ NELEC,NFROZ,CARITH,
     &                  DOSORT,DOENER,DOGRAD,DOLINR,DOEXC,
     &                  MWORD,DEBUG,TIMING,IPRNT,INTERFACE,
     &                  MIOBUF
      NAMELIST/CCSORT/ SRTALL
      NAMELIST/CCENER/ DOMP2,DOCCSD,DOCCSDT,MAXIT,MAXDIM,NTOL
      NAMELIST/CCGRAD/ DOMP2G,DOCCSDG,DOCCSDTG,NEOPER,NAMEE,
     &                 MAXIT,MAXDIM,NTOL
      NAMELIST/CCLINR/ NAOPER,NBOPER,NFREQ,NAMEA,NAMEB,EFREQ,
     &                 MAXIT,MAXDIM,NTOL
      NAMELIST/CCEXC / NEXC,MAXIT,MAXDIM,NTOL
      NAMELIST/CCEND / REWINP
      DOUBLE PRECISION D0
      PARAMETER (D0 = 0.0D00)
C
C---------------Executable code--------------------------------------
C
C     Read number of occupied spinors from input
C
      MIOBUF = 0
      REWIND (IN)
      DO I = 1, 16
         NELEC(I) = 0
         NFROZ(I) = 0
      ENDDO
      IPRNT = 0
      DOSORT  = .TRUE.
      DOENER  = .TRUE.
      DOGRAD  = .FALSE.
      DOLINR  = .FALSE.
      DOEXC   = .FALSE.
      DEBUG  = .FALSE.
      TIMING = .FALSE.
      MWORD = 0
      READ (IN,RELCCSD,END=1,ERR=1)
      GOTO 2
    1 CALL QUIT('Error exit while reading namelist RELCCSD')
    2 CONTINUE
C
C     Input for the energy calculation.
C
      MAXDIM = 8
      MAXIT = 30
      NTOL = 12
      DOMP2   = .TRUE.
      DOCCSD  = .FALSE.
      DOCCSDT = .FALSE.
      REWIND (IN)
      READ (IN,CCENER,END=4,ERR=3)
      GOTO 4
    3 CALL QUIT('Error exit while reading namelist CCENER')
    4 CONTINUE
      MXDIMCC = MAXDIM
      MXITCC = MAXIT
      NTOLCC = NTOL
C
C     Input for the gradient calculation.
C
      DOMP2G   = .FALSE.
      DOCCSDG  = .FALSE.
      DOCCSDTG = .FALSE.
      NEOPER = 1
      NAMEE(1) = 'ZDIPLEN '
      MAXDIM = 8
      MAXIT = 30
      NTOL = 12
      REWIND (IN)
      READ (IN,CCGRAD,END=6,ERR=5)
      GOTO 6
    5 CALL QUIT('Error exit while reading namelist CCGRAD')
    6 CONTINUE
      MXDIMGR = MAXDIM
      MXITGR = MAXIT
      NTOLGR = NTOL
C
C     Input for the linear response calculation.
C
      NAOPER = 1
      NBOPER = 1
      NFREQ = 1
      DO I = 1, 16
         NEXC(I) = 1
      ENDDO
      NAMEA(1) = 'ZDIPLEN '
      NAMEB(1) = 'ZDIPLEN '
      DO IFREQ = 1, MAXFR
         EFREQ(IFREQ) = D0
      ENDDO
      MAXDIM = 8
      MAXIT = 30
      NTOL = 12
      REWIND (IN)
      READ (IN,CCLINR,END=8,ERR=7)
      GOTO 8
    7 CALL QUIT('Error exit while reading namelist CCLINR')
    8 CONTINUE
      MXDIMRP = MAXDIM
      MXITRP = MAXIT
      NTOLRP = NTOL
C
C     Input for the excitation energy calculation.
C
      DO I = 1, 16
         NEXC(I) = 1
      ENDDO
      MAXDIM = 8
      MAXIT = 30
      NTOL = 12
      REWIND (IN)
      READ (IN,CCEXC,END=10,ERR=9)
      GOTO 10
    9 CALL QUIT('Error exit while reading namelist CCEXC')
   10 CONTINUE
      MXDIMEX = MAXDIM
      MXITEX = MAXIT
      NTOLEX = NTOL
C
C     Position the file at the right spot when this routine is
C     called as a part of DIRAC
C
      REWIND (IN)
      REWINP = .FALSE.
      READ (IN,CCEND,END=100,ERR=99)
      GOTO 100
   99 CALL QUIT('Error exit while reading namelist CCEND')
  100 CONTINUE
      IF (REWINP) REWIND (IN)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if defined (CC_ALONE)
      SUBROUTINE CCMAIN (MAXCORE)
#else
      SUBROUTINE CCMAIN (CC,MAXCORE)
#endif
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Driver for Relativistic Coupled Cluster program
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "iobuf.inc"
      INCLUDE "diis.inc"
      INCLUDE "inpt.inc"
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
C     CC is the static memory  (when called as subroutine)
C     CC is the dynamic memory (when used as stand-alone code)
C
C
#if defined (CC_ALONE)
#if defined (SYS_CRAY) || defined (SYS_AIX) || defined (SYS_IRIX) || defined (SYS_HPUX) || defined (SYS_T3D)
      POINTER (ICC,CC(1))
      REAL*8 CC
#else
      check dynamic memory allocation
#endif
#else
      REAL*8 CC(MAXCORE)
#endif
      PARAMETER (MAXSP=1000)
      INTEGER IREPSP(MAXSP*3),IREPSPI(MAXSP*16)
      REAL*8 EPS(MAXSP),EPS1(MAXSP)
      LOGICAL ALLOC
C
C---------------Executable code--------------------------------------
C
C     Memory and I/O setup section
C     ----------------------------
C
      MAXCOR = MAXCORE
      MXCORR = 0
C
C     Read information from MOLFDIR
C
      OPEN (MRCONEE,FILE='MRCONEE',FORM='UNFORMATTED')
      READ (MRCONEE) NSP,BREIT,ECORE
      CLOSE (MRCONEE)
C
      IF (NSP.GT.MAXSP) THEN
         PRINT*," INCREASE MAXSP IN CCMAIN TO",NSP
         CALL QUIT("NSP.GT.MAXSP")
      ENDIF
C
C     Read input
C
      CALL CCINPT (IN,IW,MIOBUF)
C
C     Number of integer words for REAL*8 and 
C     number of REAL*8 words per variable (2 for complex groups, 
C     1 for real)
C
      IRW = INTOWP(1)
      IF (CARITH) THEN
         RCW = 2
      ELSE
         RCW = 1
      ENDIF
C
      IF (MWORD.GT.0) MAXCOR = MWORD * 1024 * 1024
      MWORD = MAXCOR / (1024 * 1024)
      IF (MOD(MAXCOR,1024*1024).NE.0) MWORD = MWORD + 1
C
C  Check the user supplied memory limit with the static maximum.
C
      IERR = 0
#if !defined (CC_ALONE)
      IF (MAXCOR.GT.MAXCORE) IERR = 1
#endif
      IF (IERR.NE.0) THEN
         MWORD = MAXCORE / (1024 * 1024)
         WRITE (IW,1006) MAXCORE,MWORD
         CALL QUIT("DECREASE MWORD")
      ENDIF
C
C     In case of dynamic allocation :
C     Allow for overhead : subtract 800 kWord
C     Allow for file buffers (Cray only) : NCCFIL files with 20 blocks each
C
#if defined (CC_ALONE)
         MAXCOR = MAXCOR - 800 * 1024
#if defined (SYS_CRAY) || defined (SYS_T3D)
         MAXCOR = MAXCOR - NCCFIL * 20 * 512
#endif
#endif
C
      WRITE (IW,1000) MAXCOR
C
C  Calculate minimum size of I/O buffers 
C  Records defined in integer*4 words, minimum size is one record
C
      NREC = 1
      MINBUF = SEC2I(NREC)
      MINSEC = IRW * RCW * NREC
C
      WRITE (IW,1005) MINBUF*RCW
C
C     General setup
C     -------------
C
      CALL CCDINI
C
C  Set up symmetry tables
C
      CALL SYMTAB (NELEC,NFROZ,NSP,EPS1,EPS,IREPSP,IREPSPI)
C
      NELECT = 0
      DO I = 1, 16
         NELECT = NELECT + NELEC(I)
      ENDDO
C
      WRITE (IW,1004) NELECT,NSP-NELECT,CARITH,DOSORT,DOENER,DOGRAD,
     &                DOLINR,DOEXC,DEBUG,TIMING,IPRNT,MWORD,INTERFACE
C
C     Check the memory requirements of the active modules
C
      ICALC = 1
      ALLOC = .FALSE.
      KFREE = 1
      LFREE = MAXCOR
      CALL MEMOPT('REAL',FOO,NFOO*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',FVO,NFVO*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',FVV,NFVV*RCW,CC,KFREE,LFREE,ALLOC)
C
      MXCORR1 = 0
      CALL CCSETI (ICALC,CC,KFREE,LFREE,MXCORR1,MREMAIN,EPS,
     &             NSP,IREPSP,IREPSPI,CC(FOO),CC(FVO),CC(FVV))
      IF (DOSORT) MXCORR = MAX0(MXCORR,MXCORR1)

      MXCORR1 = 0
      CALL CCENER (ICALC,CC,KFREE,LFREE,MXCORR1,EPS,
     &             CC(FOO),CC(FVO),CC(FVV))
      IF (DOENER) MXCORR = MAX0(MXCORR,MXCORR1)
C
      MXCORR1 = 0
      CALL CCGRAD (ICALC,NSP,IREPSPI,CC,KFREE,LFREE,MXCORR1,EPS,
     &             CC(FOO),CC(FVO),CC(FVV))
      IF (DOGRAD) MXCORR = MAX0(MXCORR,MXCORR1)
C
      MXCORR1 = 0
      CALL CCHESS (ICALC,CC,KFREE,LFREE,MXCORR1,EPS,
     &             NSP,IREPSPI,CC(FOO),CC(FVO),CC(FVV))
      IF (DOLINR) MXCORR = MAX0(MXCORR,MXCORR1)
C
      MXCORR1 = 0
      CALL CCEXC (ICALC,CC,KFREE,LFREE,MXCORR1,EPS,
     &             NSP,IREPSPI,CC(FOO),CC(FVO),CC(FVV))
      IF (DOEXC) MXCORR = MAX0(MXCORR,MXCORR1)
C
      IF (MXCORR.GT.MAXCOR) THEN
#if defined (CC_ALONE)
         MXCORR = MXCORR+800*1024
#if defined (SYS_CRAY) || defined (SYS_T3D)
         MXCORR = MXCORR+NCCFIL*20*512
#endif
#endif
         MWORD = MXCORR / (1024 * 1024)
         IF (MOD(MXCORR,1024*1024).NE.0) MWORD = MWORD + 1
         WRITE (IW,1002) MXCORR,MWORD
         CALL QUIT("INCREASE MWORD")
      ENDIF
C
C  Allocate the memory
C
#if defined (CC_ALONE)
      IERR = 0
      IABORT = 0
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL HPALLOC (ICC,MXCORR,IERR,IABORT)
#endif
#if defined (SYS_AIX)
      ICC =  MALLOC (%VAL(MXCORR*8))
      IF (ICC.EQ.0) IERR = 1
#endif
#if defined (SYS_IRIX) || defined (SYS_HPUX)
      ICC =  MALLOC (MXCORR*8)
      IF (ICC.EQ.0) IERR = 1
#endif
      IF (IERR.NE.0) THEN
          PRINT*, " Problems allocating ",MXCORR," words"
          CALL QUIT("MEMORY PROBLEMS")
      ENDIF
#endif
C
      WRITE (IW,1001) MXCORR
C
C  Use remaining free memory as extra buffer for integral and
C  coefficient files.
C  At present only used in Cray version.
C
#if defined (CC_ALONE)
      MREMAIN = MAXCOR - MXCORR
#else
#if defined (SYS_CRAY) || defined (SYS_T3D)
      MREMAIN = NCCFIL*20*512
      IF(MIOBUF.NE.0) THEN
        IF(MIOBUF.LT.MREMAIN) THEN
          PRINT*, "MIOBUF too small; need at least ", MREMAIN
          CALL QUIT("Too small IO-buffer")
        ENDIF
        MREMAIN = MIOBUF
      ENDIF
#endif
#endif
C
      NINTSB = MAX0(MREMAIN/512,20)
      NINTSB  = MIN0(MXDIMCC*(NV3*RCW/512+FVO*RCW/512+2),NINTSB)
CTROND      CALL RFILE(ITAMPL,NINTSB)
CTROND      MREMAIN = MREMAIN - (NINTSB-20)*512
C
      NINTSB = MAX0(MREMAIN/512,20)
      NINTSB  = MIN0(MXDIMCC*(NV3*RCW/512+FVO*RCW/512+2),NINTSB)
CTROND      CALL RFILE(ITAMPD,NINTSB)
CTROND      MREMAIN = MREMAIN - (NINTSB-20)*512
C
C     Start the calculations
C     ----------------------
C
      ICALC = 2
      ALLOC = .TRUE.
      KFREE = 1
      LFREE = MXCORR
      CALL MEMOPT('REAL',FOO,NFOO*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',FVO,NFVO*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',FVV,NFVV*RCW,CC,KFREE,LFREE,ALLOC)
C
      IF (DOSORT) CALL CCSETI (ICALC,CC,KFREE,LFREE,MXCORR,MREMAIN,EPS,
     &            NSP,IREPSP,IREPSPI,CC(FOO),CC(FVO),CC(FVV))
      CALL FLSHFO (IW)
C
      IF (DOENER) CALL CCENER (ICALC,CC,KFREE,LFREE,MXCORR,EPS,
     &            CC(FOO),CC(FVO),CC(FVV))
      CALL FLSHFO (IW)
C
      IF (DOGRAD) CALL CCGRAD (ICALC,NSP,IREPSPI,CC,KFREE,LFREE,
     &            MXCORR,EPS,CC(FOO),CC(FVO),CC(FVV))
      CALL FLSHFO (IW)
C
      IF (DOLINR) CALL CCHESS (ICALC,CC,KFREE,LFREE,MXCORR,EPS,
     &                         NSP,IREPSPI,CC(FOO),CC(FVO),CC(FVV))
      CALL FLSHFO (IW)
C
      IF (DOEXC ) CALL CCEXC (ICALC,CC,KFREE,LFREE,MXCORR,EPS,
     &                        NSP,IREPSPI,CC(FOO),CC(FVO),CC(FVV))
      CALL FLSHFO (IW)
C
C  FINISHED; RETURN AFTER CLOSING FILES
C
      JCODE = 3
      CALL RCLOSE(ITAPT+0,JCODE)
      JCODE = 3
      CALL RCLOSE(ITAPT+1,JCODE)
      JCODE = 3
      CALL RCLOSE(ITAPT+2,JCODE)
      JCODE = 3
      CALL RCLOSE(ITAPT+3,JCODE)
      JCODE = 3
      CALL RCLOSE(ITAPT+4,JCODE)
      JCODE = 3
      CALL RCLOSE(ITAPT+5,JCODE)
      JCODE = 3
      CALL RCLOSE(ITAPT+6,JCODE)
C
      CALL CCSUMM
C
      RETURN
C
1000  FORMAT (/" Maximum core available :",T50,I10," 8-byte words")
1001  FORMAT (" Maximum core used for active modules :",
     &        T50,I10," 8-byte words")
1002  FORMAT (" Not enough core, increase to at least",
     &        T50,I10," 8-byte words"//" Make Mword equal to ",I5)
1003  FORMAT (/" Size of 3- and 4-virtual batches :",T50,I10
     &        /" Number of passes :",T50,I10)
1004  FORMAT (//" Number of electrons :",T40,I5
     &        /" Number of virtual spinors :",T40,I5
     &        /" Complex arithmetic mode :",T44,L1
     &        /" Do integral sorting     :",T44,L1
     &        /" Do energy calculation   :",T44,L1
     &        /" Do gradient calculation :",T44,L1
     &        /" Do response calculation :",T44,L1
     &        /" Calculate exc. energies :",T44,L1
     &        /" Debug information       :",T44,L1
     &        /" Timing information      :",T44,L1
     &        /" Print level          :",T40,I5
     &        /" Memory limit (MWord) :",T40,I5
     &        /" Interface used       :",T40,A10//)
1005  FORMAT (" Size of I/O buffers :",T50,I10," 8-byte words")
1006  FORMAT (" Too much core specified, maximum is",
     &        T50,I10," 8-byte words"//" Make Mword equal to ",I5)
1008  FORMAT (" Amount of memory used :",T50,F10.1," megawords")
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SYMTAB (NELEC,NFROZ,NSP,EPS1,EPS,IREPSP,IREPSPI)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Read symmetry information from MOLFDIR
C     Generate symmetry pointer arrays
C
C---------------Routines called----------------------------------------
C
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER IREPSP(NSP,3),IREPSPI(NSP,8,2),NELEC(16),NFROZ(16)
      REAL*8 EPS1(*),EPS(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "files.inc"
      INCLUDE "symm.inc"
      INCLUDE "freeze.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER*14 REPNT(16)
      INTEGER NVIRT(16)
      INTEGER NOOOT2(8),NVVVT2(8)
C
C---------------Executable code--------------------------------------
C
C     Read symmetry information
C
      OPEN (MRCONEE,FILE='MRCONEE',FORM='UNFORMATTED')
      READ (MRCONEE) 
      READ (MRCONEE) NSYMRP,(REPNT(IRP),IRP=1,NSYMRP)
      READ (MRCONEE) NREP,(REPNA(IRP),IRP=1,2*NREP)
      READ (MRCONEE) ((MULTB(I,J,1),I=1,2*NREP),J=1,2*NREP)
      READ (MRCONEE) (IREPSP(ISP,1),IREPSP(ISP,2),EPS1(ISP),ISP=1,NSP)
      CLOSE (MRCONEE,STATUS='KEEP')
C
C     Right adjust representation names
C
      DO IRP = 1, NSYMRP
         REPN(IRP) = REPNT(IRP)(1:4)
         CALL BLANKR (4,REPN(IRP))
      ENDDO
      DO IRP = 1, NREP
         CALL BLANKR (4,REPNA(IRP))
      ENDDO
C
C     Make inverse of multiplication table
C
      DO JREP = 1, 2*NREP
         DO IREP = 1, 2*NREP
            KREP = MULTB(IREP,JREP,1)
            MULTB(IREP,KREP,2) = JREP
         ENDDO
      ENDDO
C
C     Get rid of offset for boson reps
C
      DO JREP = 1, NREP
         DO IREP = 1, NREP
            MULTB(IREP,JREP,1) = MULTB(IREP,JREP,1) - NREP
            MULTB(IREP,JREP,2) = MULTB(IREP,JREP,2) - NREP
         ENDDO
      ENDDO
      DO JREP = NREP + 1, 2*NREP
         DO IREP = NREP + 1, 2*NREP
            MULTB(IREP,JREP,1) = MULTB(IREP,JREP,1) - NREP
            MULTB(IREP,JREP,2) = MULTB(IREP,JREP,2) - NREP
         ENDDO
      ENDDO
C
C     Initialize dimension arrays
C
      DO IREP = 1, NREP
         NO(IREP) = 0
         NV(IREP) = 0
         NC(IREP) = 0
         MOO(IREP) = 0
         MVO(IREP) = 0
         MVV(IREP) = 0
         NOO(IREP) = 0
         NOV(IREP) = 0
         NVO(IREP) = 0
         NVV(IREP) = 0
         NOOT(IREP) = 0
         NVVT(IREP) = 0
         NOOOT(IREP) = 0
         NVOOT(IREP) = 0
         NVVOT(IREP) = 0
         NOOVT(IREP) = 0
         NVVVT(IREP) = 0
         NVOO (IREP) = 0
         NOVO (IREP) = 0
         NOVVT(IREP) = 0
         NOOOT2(IREP) = 0
         NVVVT2(IREP) = 0
      ENDDO
C
C     Convert input information to abelian subgroup
C
      NFROZEN = 0
      DO IRP = 1, NSYMRP
         NVIRT(IRP) = 0
         IN  = 0
         DO I = 1, NSP
            IF (IREPSP(I,1).EQ.IRP) THEN
               IN = IN + 1
               IREP = IREPSP(I,2)
               IF (IN.LE.NELEC(IRP)) THEN
                  NO(IREP) = NO(IREP) + 1
                  IREPSP(I,3) = NO(IREP)
                  IREPSPI(NO(IREP),IREP,1) = I
                  IF (NO(IREP).GT.MAXFROZ) CALL QUIT("INCREASE MAXFROZ")
                  IF (IN.LE.NFROZ(IRP)) THEN
                     NC(IREP) = NC(IREP) + 1
                     IFROZ(NO(IREP),IREP) =  NC(IREP)
                     NFROZEN = NFROZEN + 1
                  ELSE
                     IFROZ(NO(IREP),IREP) =  0
                  ENDIF
               ELSE
                  NVIRT(IRP) = NVIRT(IRP) + 1
                  NV(IREP) = NV(IREP) + 1
                  IREPSP(I,3) = - NV(IREP)
                  IREPSPI(NV(IREP),IREP,2) = I
               ENDIF
            ENDIF
         ENDDO
      ENDDO
C
      WRITE (IW,1000) (REPN(I),I=1,NSYMRP)
      WRITE (IW,1002) (NELEC(I),I=1,NSYMRP)
      WRITE (IW,1002) (NVIRT(I),I=1,NSYMRP)
      WRITE (IW,1001) (REPNA(I),I=1,NREP)
      WRITE (IW,1002) (NO(I),I=1,NREP)
      WRITE (IW,1002) (NV(I),I=1,NREP)
C
      IF (NFROZEN.GT.0) THEN
         WRITE (IW,1003)
         WRITE (IW,1000) (REPN(I),I=1,NSYMRP)
         WRITE (IW,1002) (NFROZ(I),I=1,NSYMRP)
         WRITE (IW,1001) (REPNA(I),I=1,NREP)
         WRITE (IW,1002) (NC(I),I=1,NREP)
      ENDIF
C
C     Set up dimension arrays
C
      DO JREP = 1, NREP
         DO IREP = 1, NREP
            IJREP = MULTB(IREP,JREP,1)
            IIOV(IREP,JREP) = NOV(IJREP)
            IIVO(IREP,JREP) = NVO(IJREP)
            IIVV(IREP,JREP) = NVV(IJREP)
            NOO(IJREP) = NOO(IJREP) + NO(IREP)*NO(JREP)
            NOV(IJREP) = NOV(IJREP) + NO(IREP)*NV(JREP)
            NVO(IJREP) = NVO(IJREP) + NV(IREP)*NO(JREP)
            NVV(IJREP) = NVV(IJREP) + NV(IREP)*NV(JREP)
            IJREP = MULTB(IREP,JREP,2)
            JJOO(IREP,JREP) = MOO(IJREP)
            JJVO(IREP,JREP) = MVO(IJREP)
            JJVV(IREP,JREP) = MVV(IJREP)
            MOO(IJREP) = MOO(IJREP) + NO(IREP)*NO(JREP)
            MVO(IJREP) = MVO(IJREP) + NV(IREP)*NO(JREP)
            MVV(IJREP) = MVV(IJREP) + NV(IREP)*NV(JREP)
         ENDDO
      ENDDO
C
      DO JREP = 1, NREP
         IJREP = MULTB(JREP,JREP,1)
         IIOOT(JREP,JREP) = NOOT(IJREP)
         IIVVT(JREP,JREP) = NVVT(IJREP)
         NOOT(IJREP) = NOOT(IJREP) + NO(JREP)*(NO(JREP)-1)/2
         NVVT(IJREP) = NVVT(IJREP) + NV(JREP)*(NV(JREP)-1)/2
         DO IREP = JREP + 1, NREP
            IJREP = MULTB(IREP,JREP,1)
            IIOOT(IREP,JREP) = NOOT(IJREP)
            IIVVT(IREP,JREP) = NVVT(IJREP)
            NOOT(IJREP) = NOOT(IJREP) + NO(IREP)*NO(JREP)
            NVVT(IJREP) = NVVT(IJREP) + NV(IREP)*NV(JREP)
         ENDDO
      ENDDO
C
      DO KREP = 1, NREP
         DO IJREP = 1, NREP
            IJKREP = MULTB(KREP,IJREP+NREP,2)
            KKOOOT(IJREP,KREP) = NOOOT(IJKREP)
            KKVVOT(IJREP,KREP) = NVVOT(IJKREP)
            KKOOVT(IJREP,KREP) = NOOVT(IJKREP)
            KKVOO (IJREP,KREP) = NVOO (IJKREP)
            NOOOT(IJKREP) = NOOOT(IJKREP) + NOOT(IJREP)*NO(KREP)
            NVVOT(IJKREP) = NVVOT(IJKREP) + NVVT(IJREP)*NO(KREP)
            NOOVT(IJKREP) = NOOVT(IJKREP) + NOOT(IJREP)*NV(KREP)
            NVVVT(IJKREP) = NVVVT(IJKREP) + NVVT(IJREP)*NV(KREP)
            NVOO (IJKREP) = NVOO (IJKREP) + NVO (IJREP)*NO(KREP)
         ENDDO
      ENDDO
C   
      DO IJREP = 1, NREP
         DO KREP = 1, NREP
            IJKREP = MULTB(KREP,IJREP+NREP,2)
            LLOVO (KREP,IJREP) = NOVO (IJKREP)
            LLOOOT(KREP,IJREP) = NOOOT(IJKREP)
            LLVOOT(KREP,IJREP) = NVOOT(IJKREP)
            LLOVVT(KREP,IJREP) = NOVVT(IJKREP)
            LLVVVT(KREP,IJREP) = NVVVT2(IJKREP)
            NOVO (IJKREP) = NOVO (IJKREP) + NVO(IJREP)*NO(KREP)
            NVOOT(IJKREP) = NVOOT(IJKREP) + NOOT(IJREP)*NV(KREP)
            NOVVT(IJKREP) = NOVVT(IJKREP) + NVVT(IJREP)*NO(KREP)
            NOOOT2(IJKREP) = NOOOT2(IJKREP) + NOOT(IJREP)*NO(KREP)
            NVVVT2(IJKREP) = NVVVT2(IJKREP) + NVVT(IJREP)*NV(KREP)
         ENDDO
      ENDDO
C
      NFOO = MOO(1)
      NFVO = MVO(1)
      NFVV = MVV(1)
      NDIMT1 = 0
      NDIMT2 = 0
      DO IREP = 1, NREP
         NDIMT1 = NDIMT1 + NV(IREP)*NO(IREP)
         NDIMT2 = NDIMT2 + NVVT(IREP)*NOOT(IREP)
      ENDDO
C
C     Set up pointer arrays : IXXXX for <IJ,KL:KLREP> ordening
C                             JXXXX for (IK|JL:JLREP) ordening
C                             KXXXX for <IJK,L:LREP> ordening
C                             LXXXX for <I,JKL:JKLREP> ordening
C
      IO(1) = 0
      IV(1) = 0
C
      DO IREP = 1, NREP
         IO(IREP+1) = IO(IREP) + NO(IREP)
         IV(IREP+1) = IV(IREP) + NV(IREP)
      ENDDO
C
      IOO(1) = 0
      IVO(1) = 0
      IVV(1) = 0
      IOOOOTT(1) = 0
      IVOOOT (1) = 0
      IVVOOTT(1) = 0
      IVOVO  (1) = 0
      IVOVVT (1) = 0
      IVVVVTT(1) = 0
      IVVOO  (1) = 0
      IVVOOT (1) = 0
C
      NVVVOTMX = 1
      DO IREP = 1, NREP
         IOO    (IREP+1) = IOO    (IREP) + NO  (IREP)*NO  (IREP)
         IVO    (IREP+1) = IVO    (IREP) + NV  (IREP)*NO  (IREP)
         IVV    (IREP+1) = IVV    (IREP) + NV  (IREP)*NV  (IREP)
         IOOOOTT(IREP+1) = IOOOOTT(IREP) + NOOT(IREP)*NOOT(IREP)
         IVOOOT (IREP+1) = IVOOOT (IREP) + NVO (IREP)*NOOT(IREP)
         IVVOOTT(IREP+1) = IVVOOTT(IREP) + NVVT(IREP)*NOOT(IREP)
         IVOVO  (IREP+1) = IVOVO  (IREP) + NVO (IREP)*NVO (IREP)
         IVOVVT (IREP+1) = IVOVVT (IREP) + NVO (IREP)*NVVT(IREP)
         IVVVVTT(IREP+1) = IVVVVTT(IREP) + NVVT(IREP)*NVVT(IREP)
         IVVOO  (IREP+1) = IVVOO  (IREP) + NVV (IREP)*NOO (IREP)
         IVVOOT (IREP+1) = IVVOOT (IREP) + NVV (IREP)*NOOT(IREP)
         IF (NVVT(IREP)*NVO(IREP).GT.NVVVOTMX) 
     &         NVVVOTMX=NVVT(IREP)*NVO(IREP)
      ENDDO
C
      NV1 = IOOOOTT(NREP+1)
      NV2 = IVOOOT (NREP+1)
      NV3 = IVVOOTT(NREP+1)
      NV4 = IVOVO  (NREP+1)
      NV5 = IVOVVT (NREP+1)
      NV6 = IVVVVTT(NREP+1)
C
      JVOXX(1) = 0
      JOOOO(1) = 0
      JVOOO(1) = 0
      JVVOO(1) = 0
      JVOVO(1) = 0
      JVVVO(1) = 0
      J2VOVO(1) = 0
      DO IREP = 1, NREP
         JVOXX(IREP+1) = JVOXX(IREP) + MVO(IREP)
         JOOOO(IREP+1) = JOOOO(IREP) + MOO(IREP)*MOO(IREP)
         JVOOO(IREP+1) = JVOOO(IREP) + MVO(IREP)*MOO(IREP)
         JVVOO(IREP+1) = JVVOO(IREP) + MVV(IREP)*MOO(IREP)
         JVOVO(IREP+1) = JVOVO(IREP) + MVO(IREP)*MVO(IREP)
         JVVVO(IREP+1) = JVVVO(IREP) + MVV(IREP)*MVO(IREP)
         JREP = MULTB(IREP+NREP,1+NREP,2)
         J2VOVO(IREP+1) = J2VOVO(IREP) + MVO(IREP)*MVO(JREP)
      ENDDO
      NDIMX = JVOXX(NREP+1)
C
      KVOOV (1) = 0
      KVOOO (1) = 0
      KOOOOT(1) = 0
      KVVOOT(1) = 0
      KOOVOT(1) = 0
      KOOOVT(1) = 0
      KOOVVT(1) = 0
      KVVVOT(1) = 0
      DO IREP = 1, NREP
         KVOOV  (IREP+1) = KVOOV  (IREP) + NVOO (IREP)*NV(IREP)
         KVOOO  (IREP+1) = KVOOO  (IREP) + NVOO (IREP)*NO(IREP)
         KOOOOT (IREP+1) = KOOOOT (IREP) + NOOOT(IREP)*NO(IREP)
         KVVOOT (IREP+1) = KVVOOT (IREP) + NVVOT(IREP)*NO(IREP)
         KOOVOT (IREP+1) = KOOVOT (IREP) + NOOVT(IREP)*NO(IREP)
         KOOOVT (IREP+1) = KOOOVT (IREP) + NOOOT(IREP)*NV(IREP)
         KOOVVT (IREP+1) = KOOVVT (IREP) + NVOOT(IREP)*NV(IREP)
         KVVVOT (IREP+1) = KVVVOT (IREP) + NVVVT(IREP)*NO(IREP)
      ENDDO
C
      LOOVO (1) = 0
      LVOVO (1) = 0
      LOOOOT(1) = 0
      LVVOOT(1) = 0
      LOVOOT(1) = 0
      LVOOOT(1) = 0
      LOVVOT(1) = 0
      DO IREP = 1, NREP
         LOOVO  (IREP+1) = LOOVO  (IREP) + NO(IREP)*NOVO (IREP)
         LVOVO  (IREP+1) = LVOVO  (IREP) + NV(IREP)*NOVO (IREP)
         LOOOOT (IREP+1) = LOOOOT (IREP) + NO(IREP)*NOOOT(IREP)
         LVVOOT (IREP+1) = LVVOOT (IREP) + NV(IREP)*NVOOT(IREP)
         LOVOOT (IREP+1) = LOVOOT (IREP) + NO(IREP)*NVOOT(IREP)
         LVOOOT (IREP+1) = LVOOOT (IREP) + NV(IREP)*NOOOT(IREP)
         LOVVOT (IREP+1) = LOVVOT (IREP) + NO(IREP)*NVVOT(IREP)
      ENDDO
C
C     3-index arrays are used in t3corr : make index arrays
C
      DO IJKREP = 1, NREP
       IVVVT(1,IJKREP) = 0
       DO JKREP = 1, NREP
        IREP = MULTB(JKREP+NREP,IJKREP,2)
        IVVVT(JKREP+1,IJKREP)=IVVVT(JKREP,IJKREP)+NV(IREP)*NVVT(JKREP)
       ENDDO
       JVVVT(1,IJKREP) = 0
       DO KREP = 1, NREP
        IJREP = MULTB(KREP,IJKREP,2)
        JVVVT(KREP+1,IJKREP)=JVVVT(KREP,IJKREP)+NVVT(IJREP)*NV(KREP)
       ENDDO
      ENDDO
      IVVVTMX = 1
      DO IJKREP = 1, NREP
      IF (IVVVT(NREP+1,IJKREP).GT.IVVVTMX) IVVVTMX=IVVVT(NREP+1,IJKREP)
      ENDDO
C
C
C     Reorder the spinor eigenvalues
C
      DO ISP = 1, NSP
         IREP = IREPSP(ISP,2)
         II = IREPSP(ISP,3)
         IF (II.GT.0) THEN
            I = IO(IREP) + II
         ELSE
            I = IO(NREP+1) + IV(IREP) - II
         ENDIF
         EPS(I) = EPS1(ISP)
      ENDDO
C
 1000 FORMAT (/' Configuration in highest pointgroup'//16(1X,A4))
 1001 FORMAT (//' Configuration in abelian subgroup'//8(1X,A4))
 1002 FORMAT (16(1X,I4))
 1003 FORMAT (//' List of frozen spinors (amplitudes zeroed out)')
      RETURN 
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCSETI (ICALC,CC,KFREE,LFREE,MXCORR,MREMAIN,EPS,
     &                   NSP,IREPSP,IREPSPI,FOO,FVO,FVV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Driver to generate the sorted integral files
C
C     Menu driven routine. Depending on value of ICALC
C
C       1) Calculate memory requirements for the different modules
C       2) Produce desired integral files or rebuild Fock matrix
C          - OOOO, VOOO, VVOO and VOVO files
C          - VOVV file
C          - VVVV file
C          - Fock matrix
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 CC(*)
      REAL*8 EPS(*),FOO(*),FVO(*),FVV(*)
      INTEGER IREPSP(NSP,3),IREPSPI(NSP,8,2)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "iobuf.inc"
      INCLUDE "diis.inc"
      INCLUDE "inpt.inc"
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL ALLOC
C
C---------------Executable code--------------------------------------
C
      ALLOC = ICALC.GE.2
C
C  Allocate core for integral sorting step
C
      KFRSAV = KFREE
      LFRSAV = LFREE
      CALL MEMOPT('INTE',KR,2*NSP+2,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('INTE',INDEX,NSP*NSP*4,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('INTE',INDK,NSP*NSP,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('INTE',INDL,NSP*NSP,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',CBUF,NSP*NSP*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('INTE',IBUF,NSP*NSP*7*16,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BUF1,NSP*NSP,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BUF2,NSP*NSP,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',V1 ,NV1*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',V1A,NV1*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',V2 ,NV2*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',V2A,NV2*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',V3 ,NV3*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',V3A,NV3*RCW,CC,KFREE,LFREE,ALLOC)
      MXCORR = MAX0(MXCORR,KFREE)
C
C     We cannot reuse the memory when the bound checking is active
C     Here is room for improvement if the memory requirements become
C     large for this step. Start after VVOO arrays.
C
C     We must be able to store at least a NV4 size array, check that.
C     Since we use this array as work in MAKEFM we also need to store the
C     first representation of the sorted (non-triangular) OOOO and VOOO ints
C
C     The maximum necessary size is
C
      NVT = MAX0(NV5,NV6)
      NVT = MAX0(NVT,NV4)
      NVT = MAX0(NVT,JOOOO(2))
      NVT = MAX0(NVT,JVOOO(2))
C
C     The remaining space is
C
      NVF = (LFREE-2) / (2 * RCW )
C
      IF (NVT.LT.NVF) THEN
         NBUF = NVT / MINBUF + 1
      ELSE
         NBUF = NVF / MINBUF
      ENDIF
      IF (NBUF.LT.1) NBUF = 1
      NVT = NBUF * MINBUF
C
      CALL MEMOPT('REAL',VT ,NVT*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',VTA,NVT*RCW,CC,KFREE,LFREE,ALLOC)
      V4  = VT
      V4A = VTA
      MXCORR = MAX0(MXCORR,KFREE)
C
      IF (ICALC.EQ.1) THEN
         WRITE (IW,1001) "reading and sorting integrals",MXCORR
         KFREE = KFRSAV
         LFREE = LFRSAV
         RETURN
      ENDIF
C
C  Read integrals & make fock matrices
C
      CALL RDINTS(IPRNT,NSP,CC(INDEX),CC(KR),CC(INDK),CC(INDL),
     &            CC(CBUF),CC(IBUF),NVT,CC(BUF1),CC(BUF2),IREPSP,
     &            IREPSPI,CC(V1),CC(V1A),
     &            CC(V2),CC(V2A),CC(V3),CC(V3A),CC(V4),CC(V4A),
     &            CC(VT),CC(VTA),CC(VT),CC(VTA),EPS,FOO,FVO,
     &            FVV,BREIT,DOCCSD,INTERFACE,MREMAIN)
C
      CALL MEMREL('CCSETI',CC,1,KFRSAV,KFREE,LFREE)
C
      RETURN
1001  FORMAT (" Core used for ",A," :",T50,I10," 8-byte words")
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCENER (ICALC,CC,KFREE,LFREE,MXCORR,EPS,FOO,FVO,FVV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Driver for energy calculations
C
C     Menu driven routine. Depending on value of ICALC
C
C       1) Calculate memory requirements for the different modules
C       2) Calculate energy and wave function at the
C          - SCF level
C          - MP2 level
C          - CCSD level
C          - CCSD(T) level
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "iobuf.inc"
      INCLUDE "diis.inc"
      INCLUDE "inpt.inc"
      INCLUDE "symm.inc"
      INCLUDE "eqns.inc"
      INCLUDE "results.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 CC(*)
      REAL*8 EPS(*),FOO(*),FVO(*),FVV(*)
      REAL*8 TOL
      LOGICAL ALLOC
C
C---------------Executable code--------------------------------------
C
      ALLOC = ICALC.GE.2
      IF (ICALC.LE.0) RETURN
C
C Allocate space for the calculation of the amplitudes
C
      KFRSAV = KFREE
      LFRSAV = LFREE
      CALL MEMOPT('REAL',T1,NDIMT1*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',T2,NDIMT2*RCW,CC,KFREE,LFREE,ALLOC)
      KFRSAV1 = KFREE
      CALL MEMOPT('REAL',S1,NDIMT1*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',S2,NDIMT2*RCW,CC,KFREE,LFREE,ALLOC)
C
C This is all we need to calculate the MP2 energy
C
      IF (DOMP2) MXCORR = KFREE
C
      CALL MEMOPT('REAL',AZ,NV1*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',H ,NV4*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',HO,NFOO*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',HV,NFVV*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',HOV,NFVO*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',GOO,NFOO*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',GVV,NFVV*RCW,CC,KFREE,LFREE,ALLOC)
      MAXDIM = MXDIMCC
      CALL MEMOPT('REAL',CCC,MAXDIM,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BB1,(MAXDIM+1)*(MAXDIM+2),CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BB2,(MAXDIM+1)*(MAXDIM+2),CC,KFREE,LFREE,ALLOC)
      NBUF1 = MAX0(JOOOO(NREP+1),JVOOO(NREP+1),
     &             JVVOO(NREP+1),JVOVO(NREP+1))
      CALL MEMOPT('REAL',BUF1,NBUF1*RCW,CC,KFREE,LFREE,ALLOC)
C
C     Calculate size of 2nd and 3rd work arrays
C
      NV4IRP = 1
      DO IRP = 1, NREP
         NV4IRP = MAX0(NV4IRP,NVV(IRP)*NVVT(IRP))
      ENDDO
      NBUF2 = (LFREE - 2) / (2 * RCW)
      IF (NBUF2.GT.NV4IRP) NBUF2 = NV4IRP
      NBUF2 = MAX0(NBUF1,NBUF2)
      NBUF3 = NBUF2
C
      IF (NBUF3.LT.JVVVO(2)) THEN
C
C        We need to increase the size of BUF3 for the in-core sort
C        in GVINTM. Move some memory from BUF2 to BUF3.
C
         NBUF3D = JVVVO(2) - NBUF3
         NBUF3  = NBUF3 + NBUF3D
         NBUF2  = NBUF2 - NBUF3D
C        Check minimum size of BUF2
         NBUF2 = MAX0(NBUF1,NBUF2)
      ENDIF
C
C     Allocate the work arrays
C
      CALL MEMOPT('REAL',BUF2,NBUF2*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BUF3,NBUF3*RCW,CC,KFREE,LFREE,ALLOC)
C        
      IF (ICALC.EQ.1) WRITE (IW,1001) "calculating amplitudes",KFREE
C
      IF (DOCCSD) THEN
         MXCORR = MAX0(MXCORR,KFREE)
      ENDIF
C
C  Calculate requirements for triples correction.
C  Use MEMOPT later.
C
      TOP = S1
      VOOO  = TOP
      T2S   = VOOO + NV2*RCW + RCW
      T2T   = T2S + IVVOOT(NREP+1)*RCW + RCW
      VVOO  = T2T + LOVVOT(NREP+1)*RCW + RCW
      VVVO  = VVOO + NDIMT2*RCW + RCW
      BUF1B = VVVO + KVVVOT(NREP+1)*RCW + RCW
      BUF2B = BUF1B + IVVVTMX*RCW + RCW
      BUF3B = BUF2B + IVVVTMX*RCW + RCW
      BUF4B = BUF3B + IVVVTMX*RCW + RCW
      TOP = BUF4B + MAX0(NVVVOTMX,NDIMT2)*RCW + RCW
      IF(ICALC.EQ.1) WRITE (IW,1001) "in core evaluation of triples",TOP
      IF (DOCCSDT) MXCORR = MAX0(MXCORR,TOP)
C
      IF (ICALC.EQ.1) THEN
         KFREE = KFRSAV
         LFREE = LFRSAV
         RETURN
      ENDIF
C
      WRITE (IW,1007) DOMP2,DOCCSD,DOCCSDT
      EQNS = 'AMPLTD'
C
C  Solve relativistic MP2 equations
C
      IF (DOMP2) THEN
         IF (IMSTAT(4).EQ.3) THEN
            WRITE (IW,1010)
            GOTO 10
         ENDIF
C
         CALL MP2EQN(EPS,FVO,CC(T1),CC(S1),CC(T2),CC(S2))
C
C        This is a restart point: save the amplitudes and update status.
C
         CALL PUTAMPT(CC(T1),CC(T2))
         IMSTAT(4) = 3
         CALL CCDUMP
   10    CONTINUE
      ENDIF
C
C  Solve relativistic CCSD equations
C
   11 IF (DOCCSD) THEN
         IF (IMSTAT(5).EQ.3) THEN
            WRITE (IW,1011)
            GOTO 20
         ENDIF
         WRITE (IW,1008) MXITCC,MXDIMCC,10.0D0**(-NTOLCC)
         MAXDIM = MXDIMCC
         MAXIT  = MXITCC
         TOL = 10.0D0**(-NTOLCC)
C
C        Initialize T1 and T2 amplitudes
C
         IF (IMSTAT(5).EQ.2) THEN
C
C           Use restart amplitudes
C
            CALL GETAMPT (CC(T1),CC(T2))
C
         ELSEIF (IMSTAT(4).EQ.3) THEN
C
C           Use MP2 amplitudes
C
            CALL GETAMPT (CC(T1),CC(T2))
C
         ELSE   
C
C           Use standard initialization
C
            CALL XCOPY (NDIMT1,FVO,1,CC(T1),1)
            CALL GETVVOO (CC(T2))
            CALL DENOM (EPS,CC(T1),CC(T2),CC(T1),CC(T2))
            CALL ZCORE (CC(T1),CC(T2))
C
         ENDIF
C
         CALL PUTAMPT(CC(T1),CC(T2))
         IMSTAT(5) = 2
         CALL CCDUMP
C
         CALL CCEQN(EPS,FOO,FVO,FVV,CC(T1),CC(S1),
     &           CC(T2),CC(S2),CC(AZ),DUM,CC(H),CC(HO),CC(HV),CC(HOV),
     &           CC(GOO),CC(GVV),CC(BUF1),CC(BUF2),CC(BUF3),
     &           CC(CCC),CC(BB1),CC(BB2),
     &           TOL,NBUF2,NBUF3)
C
C        This is a restart point: save the amplitudes and update status.
C
         CALL PUTAMPT(CC(S1),CC(S2))
         IMSTAT(5) = 3
         CALL CCDUMP
   20    CONTINUE
      ENDIF
C
C  Calculate T3 corrections
C
C
      IF (DOCCSDT) THEN
         IF (IMSTAT(6).EQ.3) THEN
            WRITE (IW,1012)
            GOTO 30
         ENDIF
         CALL MEMREL('CCENER',CC,1,KFRSAV1,KFREE,LFREE)
         CALL MEMOPT('REAL',VOOO,NV2*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',T2S,IVVOOT(NREP+1)*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',T2T,LOVVOT(NREP+1)*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',VVOO,NDIMT2*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',VVVO,KVVVOT(NREP+1)*RCW,
     &               CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',BUF1,IVVVTMX*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',BUF2,IVVVTMX*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',BUF3,IVVVTMX*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',BUF4,MAX0(NVVVOTMX,NDIMT2)*RCW,
     &               CC,KFREE,LFREE,ALLOC)
C
C        Check if CCSD amplitudes are available
C
         IF (IMSTAT(5).EQ.3) THEN
            CALL GETAMPT (CC(T1),CC(T2))
         ELSEIF (IMSTAT(4).EQ.3) THEN
            WRITE (IW,1006)
            CALL GETAMPT (CC(T1),CC(T2))
         ELSE
            DOCCSD = .TRUE.
            GOTO 11
         ENDIF
C
         CALL T3CORR (FVO,EPS,CC(T1),CC(T2),CC(T2S),CC(T2T),
     &                CC(BUF1),CC(BUF2),CC(BUF3),CC(BUF4),
     &                CC(VVVO),CC(VVOO),CC(VOOO))
C
C     This is a restart point : update status.
C
         IMSTAT(6) = 3
         CALL CCDUMP
   30    CONTINUE
      ENDIF
C
      CALL MEMREL('CCENER',CC,1,KFRSAV,KFREE,LFREE)
C
      RETURN
1001  FORMAT (" Core used for ",A," :",T50,I10," 8-byte words")
1007  FORMAT (//" Energy calculations"
     &        /" MP2 module active :",T44,L1
     &        /" CCSD module active :",T44,L1
     &        /" CCSD(T) module active :",T44,L1)
1008  FORMAT (//" CCSD options :"
     &        /" Maximum number of iterations :",T40,I5
     &        /" Maximum size of DIIS space :",T40,I5
     &        /" Convergence criterium :",T38,E7.1)
1006  FORMAT (/' WARNING : Calculating CCSD(T) energies with MP2',
     &         ' amplitudes')
1010  FORMAT (/' Skipping MP2 calculation')
1011  FORMAT (/' Skipping CCSD calculation')
1012  FORMAT (/' Skipping CCSD(T) calculation')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCGRAD (ICALC,NSP,IREPSPI,CC,KFREE,LFREE,MXCORR,
     &                   EPS,FOO,FVO,FVV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Driver for gradient calculations
C
C     Menu driven routine. Depending on value of ICALC
C
C       1) Calculate memory requirements for the different modules
C       2) Calculate first order properties at the
C          - SCF level
C          - MP2 level
C          - CCSD level
C          - CCSD(T) level
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "iobuf.inc"
      INCLUDE "diis.inc"
      INCLUDE "inpt.inc"
      INCLUDE "eqns.inc"
      INCLUDE "results.inc"
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 CC(*)
      REAL*8 EPS(*),FOO(*),FVO(*),FVV(*)
      REAL*8 TOL
      COMPLEX*16 EPHASE,CEXPA
      LOGICAL ALLOC
C
C---------------Executable code--------------------------------------
C
      if (doccsdg) then
         print*
         print*,' CCSD gradient not yet implemented,'
     &         ,' calculate MP2 gradient instead'
         print*
         doccsdg = .false.
         domp2g = .true.
      endif
      ALLOC = ICALC.GE.2
      IF (ICALC.LE.0) RETURN
C
C The DIIS algorithm is used for the Z-vector and Lamda-equations
C Allocate the workarrays only once with the maximum dimensions
C
      MAXIT  = MXITGR
      MAXDIM = MXDIMGR
      IF (DOCCSDG) MAXIT = MAX0(MAXIT,MXITCC)
      IF (DOCCSDG) MAXDIM = MAX0(MAXDIM,MXDIMCC)
C
C Allocate space for the calculation of the gradient
C
      KFRSAV = KFREE
      LFRSAV = LFREE
      CALL MEMOPT('REAL',T1,NDIMT1*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',T2,NDIMT2*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',S1,NDIMT1*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',S2,NDIMT2*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',DOO,NFOO*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',DVV,NFVV*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',XOO,NFOO*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',XVV,NFVV*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',DVO,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',DOV,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',D2VO,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',D2OV,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',XVO,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',XOV,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',EVO,NEOPER*NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',EOV,NEOPER*NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',EOO,NFOO*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',EVV,NFVV*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',DG ,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      NBUF1 = MAX0(JOOOO(NREP+1),JVOOO(NREP+1),
     &             JVVOO(NREP+1),JVOVO(NREP+1))
      CALL MEMOPT('REAL',CCC,MAXDIM,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BB1,(MAXDIM+1)*(MAXDIM+2),CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BB2,(MAXDIM+1)*(MAXDIM+2),CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BUF1,NBUF1*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',CBUF,2*NSP*NSP,CC,KFREE,LFREE,ALLOC)
C
C     Calculate size of 2nd and 3rd work arrays
C
      NV4IRP = 1
      DO IRP = 1, NREP
         NV4IRP = MAX0(NV4IRP,NVV(IRP)*NVVT(IRP))
      ENDDO
      NBUF2 = (LFREE - 2) / (2 * RCW)
      IF (NBUF2.GT.NV4IRP) NBUF2 = NV4IRP
      NBUF2 = MAX0(NBUF1,NBUF2)
      NBUF3 = NBUF2
C
      IF (NBUF3.LT.JVVVO(2)) THEN
C
C        We need to increase the size of BUF3 for the in-core sort
C        in GVINTM. Move some memory from BUF2 to BUF3.
C
         NBUF3D = JVVVO(2) - NBUF3
         NBUF3  = NBUF3 + NBUF3D
         NBUF2  = NBUF2 - NBUF3D
C        Check minimum size of BUF2
         NBUF2 = MAX0(NBUF1,NBUF2)
      ENDIF
C
C     Allocate the work arrays
C
      CALL MEMOPT('REAL',BUF2,NBUF2*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BUF3,NBUF3*RCW,CC,KFREE,LFREE,ALLOC)
C
      IF (DOCCSDG) THEN
         CALL MEMOPT('REAL',AZ,NV1*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',H ,NV4*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',HO,NFOO*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',HV,NFVV*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',HOV,NFVO*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',GOO,NFOO*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',GVV,NFVV*RCW,CC,KFREE,LFREE,ALLOC)
         CALL MEMOPT('REAL',CZ,NV2*RCW,CC,KFREE,LFREE,ALLOC)
      ENDIF
C
      IF (DOCCSDG) MXCORR = KFREE
C
      IF (ICALC.EQ.1) THEN
         WRITE (IW,1001) "Gradient calculation",KFREE
         KFREE = KFRSAV
         LFREE = LFRSAV
         RETURN
      ENDIF
C
      WRITE (IW,1007) DOMP2G,DOCCSDG,DOCCSDTG,
     &                MXITGR,MXDIMGR,10.0D0**(-NTOLGR)
      EQNS = 'LAMBDA'
C
C     Start with SCF (reference determinant) expectation values
C
      WRITE (IW,1002) '    SCF'
      CALL XCOPY (NFOO,A0,0,CC(DOO),1)
      CALL XCOPY (NFVV,A0,0,CC(DVV),1)
      CALL XCOPY (NDIMX,A0,0,CC(DVO),1)
      CALL XCOPY (NDIMX,A0,0,CC(DOV),1)
C
      DO EOPER = 1, NEOPER
         EOFF = (EOPER-1)*NDIMX*RCW+1
         CALL RDPROP(IPRNT,NAMEE(EOPER),NSP,LTR,CC(CBUF),IREPSPI,
     &               EPHASE,CC(EVO),CC(EOV),CC(EOO),CC(EVV))
         CALL EXPVAL(CC(DOO),CC(DVV),CC(DVO),CC(DOV),
     &               CC(EOO),CC(EVV),CC(EVO),CC(EOV),
     &               EPHASE,CEXPA)
         WRITE (IW,1004) NAMEE(EOPER),CEXPA
      ENDDO
C
C     Get MP2 relaxed density matrix
C
      IF (DOMP2G) THEN
C
C        Solve relativistic MP2 lambda equations
C
         CALL MP2EQN(EPS,FVO,CC(T1),CC(S1),CC(T2),CC(S2))
C
C        Get amplitude part of density matrix
C
         CALL CCDENSC (1,CC(T1),CC(S1),CC(T2),CC(S2),
     &                 CC(BUF1),CC(BUF2),CC(BUF3),CC(DOO),CC(DVV),
     &                 CC(DVO),CC(DOV))
C
C        Solve Z-vector equations
C
         TOL = 10.0D0**(-NTOLGR)
         MAXIT  = MXITGR
         MAXDIM = MXDIMGR
         CALL CCDENSZ (1,NSP,IREPSPI,EPS,CC(T1),CC(S1),CC(T2),CC(S2),
     &                 CC(BUF1),CC(BUF2),CC(BUF3),NBUF2,NBUF3,CC(CBUF),
     &                 CC(XVO),CC(XOV),CC(DOO),CC(DVV),CC(DVO),CC(DOV),
     &                 CC(EVO),CC(EOV),CC(D2VO),CC(D2OV),CC(DG),
     &                 CC(CCC),CC(BB1),CC(BB2),TOL)
C
         WRITE (IW,1002) '    MP2'
         DO EOPER = 1, NEOPER
            EOFF = (EOPER-1)*NDIMX*RCW+1
            CALL RDPROP(IPRNT,NAMEE(EOPER),NSP,LTR,CC(CBUF),IREPSPI,
     &                  EPHASE,CC(EVO),CC(EOV),CC(EOO),CC(EVV))
            CALL EXPVAL(CC(DOO),CC(DVV),CC(DVO),CC(DOV),
     &                  CC(EOO),CC(EVV),CC(EVO),CC(EOV),
     &                  EPHASE,CEXPA)
            WRITE (IW,1004) NAMEE(EOPER),CEXPA
         ENDDO
C
      ENDIF
C
C     Get CCSD relaxed density matrix
C
      IF (DOCCSDG) THEN
C
C        Solve relativistic CCSD lambda equations
C
C        The CCSD energy calculation should be completed
C
         IF (IMSTAT(5).NE.3) THEN
            WRITE (IW,1009) 
            RETURN
         ENDIF
C
C        No restart: start by copying T1 and T2 amplitdues to L1 and L2
C
         CALL GETAMPT (CC(T1),CC(T2))
C
         WRITE (IW,1008) MXITCC,MXDIMCC,10.0D0**(-NTOLCC)
         MAXDIM = MXDIMCC
         MAXIT  = MXITCC
         TOL = 10.0D0**(-NTOLCC)
         NINTSB = 20
         CALL RFILE(ITAPL,NINTSB)
C
         CALL CCEQN(EPS,FOO,FVO,FVV,CC(T1),CC(S1),
     &           CC(T2),CC(S2),CC(AZ),CC(CZ),CC(H),CC(HO),CC(HV),
     &           CC(HOV),CC(GOO),CC(GVV),CC(BUF1),CC(BUF2),CC(BUF3),
     &           CC(CCC),CC(BB1),CC(BB2),
     &           TOL,NBUF2,NBUF3)
C
C        Put T1 and T2 amplitudes back in place
C        L1 and L2 amplitudes are in S1, S2
C
         CALL GETAMPT(CC(T1),CC(T2))
C
C        Get amplitude part of density matrix
C
         CALL CCDENSC (1,CC(T1),CC(S1),CC(T2),CC(S2),
     &                 CC(BUF1),CC(BUF2),CC(BUF3),CC(DOO),CC(DVV),
     &                 CC(DVO),CC(DOV))
C
C        Solve Z-vector equations
C
         TOL = 10.0D0**(-NTOLGR)
         MAXIT  = MXITGR
         MAXDIM = MXDIMGR
c        CALL CCDENSZ (1,NSP,IREPSPI,EPS,CC(T1),CC(S1),CC(T2),CC(S2),
c    &                 CC(BUF1),CC(BUF2),CC(BUF3),NBUF2,NBUF3,CC(CBUF),
c    &                 CC(XVO),CC(XOV),CC(DOO),CC(DVV),CC(DVO),CC(DOV),
c    &                 CC(EVO),CC(EOV),CC(D2VO),CC(D2OV),CC(DG),
c    &                 CC(CCC),CC(BB1),CC(BB2),TOL)
C
         WRITE (IW,1002) '   CCSD'
         DO EOPER = 1, NEOPER
            EOFF = (EOPER-1)*NDIMX*RCW+1
            CALL RDPROP(IPRNT,NAMEE(EOPER),NSP,LTR,CC(CBUF),IREPSPI,
     &                  EPHASE,CC(EVO),CC(EOV),CC(EOO),CC(EVV))
            CALL EXPVAL(CC(DOO),CC(DVV),CC(DVO),CC(DOV),
     &                  CC(EOO),CC(EVV),CC(EVO),CC(EOV),
     &                  EPHASE,CEXPA)
            WRITE (IW,1004) NAMEE(EOPER),CEXPA
         ENDDO
C
      ENDIF
C
      CALL MEMREL('CCGRAD',CC,1,KFRSAV,KFREE,LFREE)
C
      RETURN
1001  FORMAT (" Core used for ",A," :",T50,I10," 8-byte words")
1007  FORMAT (//" Gradient calculations"
     &        /" MP2 module active :",T44,L1
     &        /" CCSD module active :",T44,L1
     &        /" CCSD(T) module active :",T44,L1
     &        //" Z-vector equations :"
     &        /" Maximum number of iterations :",T40,I5
     &        /" Maximum size of DIIS space :",T40,I5
     &        /" Convergence criterium :",T38,E7.1)
1002  FORMAT(//,1X,A7,' Gradient results',//3X,'Expectation value',13X,
     & 'Real part',8X,'Imaginary part')
1004  FORMAT('  <',A8,'> ',T30,2F20.13)
1008  FORMAT (//" CCSD lambda-equations :"
     &        /" Maximum number of iterations :",T40,I5
     &        /" Maximum size of DIIS space :",T40,I5
     &        /" Convergence criterium :",T38,E7.1)
1009  FORMAT (/" WARNING : Cannot solve CCSD lambda equations :",
     &         " T1 and T2 amplitudes are not available "/
     &         " Check that you specified DOENER and DOCCSD !")
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCHESS (ICALC,CC,KFREE,LFREE,MXCORR,EPS,
     &                   NSP,IREPSPI,FOO,FVO,FVV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Driver for hessian calculations
C
C     Menu driven routine. Depending on value of ICALC
C
C       1) Calculate memory requirements for the different modules
C       2) Calculate second order properties at the
C          - SCF (RPA) level
C          - MP2 level
C          - CCSD level
C          - CCSD(T) level
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 CC(*)
      REAL*8 EPS(*),FOO(*),FVO(*),FVV(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "iobuf.inc"
      INCLUDE "diis.inc"
      INCLUDE "inpt.inc"
      INCLUDE "symm.inc"
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 APHASE(MAXOP),BPHASE
      REAL*8 TOL
      LOGICAL BDONE
      LOGICAL ALLOC
C
C---------------Executable code--------------------------------------
C
      ALLOC = ICALC.GE.2
      IF (ICALC.LE.0) RETURN
C
C Allocate space for the calculation of the hessian
C
      KFRSAV = KFREE
      LFRSAV = LFREE
      CALL MEMOPT('REAL',X1,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',X2,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',Y1,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',Y2,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',AVO,NAOPER*NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',AOV,NAOPER*NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BVO,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BOV,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',DG ,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      MAXDIM = MXDIMRP
      CALL MEMOPT('REAL',CCC,MAXDIM,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BB1,(MAXDIM+1)*(MAXDIM+2),CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BB2,(MAXDIM+1)*(MAXDIM+2),CC,KFREE,LFREE,ALLOC)
      NBUF1 = MAX0(JVVOO(NREP+1),JVOVO(NREP+1))
      CALL MEMOPT('REAL',BUF1,NBUF1*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BUF2,NBUF1*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',CBUF,2*NSP*NSP,CC,KFREE,LFREE,ALLOC)
      MXCORR = KFREE
C
      IF (ICALC.LE.1) THEN
         WRITE (IW,1001) "RPA calculation",KFREE
         KFREE = KFRSAV
         LFREE = LFRSAV
         RETURN
      ENDIF
C
C  Solve relativistic Random Phase Approximation equations
C
      WRITE (IW,1009) MXITRP,MXDIMRP,10.0D0**(-NTOLRP)
      MAXIT  = MXITRP
      TOL = 10.0D0**(-NTOLRP)
      MAXDIM = MXDIMRP
      DO AOPER = 1, NAOPER
         AOFF = (AOPER-1)*NDIMX*RCW
         CALL RDPROP(IPRNT,NAMEA(AOPER),NSP,LFA,CC(CBUF),IREPSPI,
     &               APHASE(AOPER),CC(AVO+AOFF),CC(AOV+AOFF),DUM,DUM)
      ENDDO
      DO BOPER = 1, NBOPER
         BDONE = .FALSE.
         CALL RDPROP(IPRNT,NAMEB(BOPER),NSP,LFA,CC(CBUF),IREPSPI,
     &               BPHASE,CC(BVO),CC(BOV),DUM,DUM)
         DO IFREQ = 1, NFREQ
            IF (IPRNT.GE.1) WRITE (IW,1010) NAMEB(BOPER),EFREQ(IFREQ)
            CALL RPAA
     &           (NAOPER,NAMEA,NAMEB(BOPER),APHASE,BPHASE,BDONE,EPS,
     &           EFREQ(IFREQ),CC(AVO),CC(AOV),CC(BVO),CC(BOV),CC(X1),
     &           CC(X2),CC(Y1),CC(Y2),CC(DG),CC(BUF1),CC(BUF2),
     &           CC(CCC),CC(BB1),CC(BB2),TOL,IPRNT,DEBUG,
     &           TIMING)
            CALL FLSHFO(IW)
         ENDDO
      ENDDO
C
      RETURN
1001  FORMAT (" Core used for ",A," :",T50,I10," 8-byte words")
1009  FORMAT (//" Propagator (RPA) module"
     &        /" Maximum number of iterations :",T40,I5
     &        /" Maximum size of DIIS space :",T40,I5
     &        /" Convergence criterium :",T38,E7.1)
1010  FORMAT (//" Random Phase Approximation : ",
     &        /" Calculate propagators of ",A8," (B)",
     &        " at frequency :",F15.5)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCEXC (ICALC,CC,KFREE,LFREE,MXCORR,EPS,
     &                  NSP,IREPSPI,FOO,FVO,FVV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Driver for excited state calculations
C
C     Menu driven routine. Depending on value of ICALC
C
C       1) Calculate memory requirements for the different modules
C       2) Calculate excitation energies at the
C          - SCF (RPA) level
C          - MP2 level
C          - CCSD level
C          - CCSD(T) level
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 CC(*)
      REAL*8 EPS(*),FOO(*),FVO(*),FVV(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "iobuf.inc"
      INCLUDE "diis.inc"
      INCLUDE "inpt.inc"
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 APHASE(MAXOP),BPHASE
      REAL*8 TOL
      LOGICAL BDONE
      LOGICAL ALLOC
C
C---------------Executable code--------------------------------------
C
      ALLOC = ICALC.GE.2
      IF (ICALC.LE.0) RETURN
C
C Allocate space for the excitation energy calculation
C
      KFRSAV = KFREE
      LFRSAV = LFREE
      CALL MEMOPT('REAL',X1,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',X2,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',Y1,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',Y2,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',DG ,NDIMX*RCW,CC,KFREE,LFREE,ALLOC)
      MAXDIM = MAX0(6,2*MXDIMRP+2)
      CALL MEMOPT('REAL',FREQR,MAXDIM,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',FREQI,MAXDIM,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',FREQD,MAXDIM,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',CCC,MAXDIM*MAXDIM,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BB1,MAXDIM*MAXDIM,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BB2,MAXDIM*MAXDIM,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BB3,MAXDIM*MAXDIM,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BB4,MAXDIM*MAXDIM,CC,KFREE,LFREE,ALLOC)
      NBUF1 = MAX0(JVVOO(NREP+1),JVOVO(NREP+1))
      CALL MEMOPT('REAL',BUF1,NBUF1*RCW,CC,KFREE,LFREE,ALLOC)
      CALL MEMOPT('REAL',BUF2,NBUF1*RCW,CC,KFREE,LFREE,ALLOC)
      MXCORR = KFREE
C
      IF (ICALC.LE.1) THEN
         KFREE = KFRSAV
         LFREE = LFRSAV
         WRITE (IW,1001) "excitation energies",KFREE
         RETURN
      ENDIF
C
C  Solve relativistic Random Phase Approximation equations
C
      WRITE (IW,1009) MXITRP,MXDIMRP,10.0D0**(-NTOLRP)
      MAXIT  = MXITRP
      TOL = 10.0D0**(-NTOLRP)
C
      DO IREP = 1, NREP
        IF (NEXC(IREP).NE.0) THEN
          WRITE (IW,1011) NEXC(IREP),REPNA(IREP+NREP)
          CALL RPAB
     &        (EPS,CC(X1),CC(X2),CC(Y1),CC(Y2),CC(DG),CC(BUF1),
     &        CC(BUF2),CC(CCC),CC(BB1),CC(BB2),CC(BB3),CC(BB4),
     &        CC(FREQR),CC(FREQI),CC(FREQD),
     &        NEXC(IREP),IREP,TOL,IPRNT,DEBUG,TIMING)
          CALL FLSHFO(IW)
        ENDIF
      ENDDO
C
      RETURN
1001  FORMAT (" Core used for ",A," :",T50,I10," 8-byte words")
1009  FORMAT (//" Excitation energy module (RPA)"
     &        /" Maximum number of iterations :",T40,I5
     &        /" Maximum size of reduced space :",T40,I5
     &        /" Convergence criterium :",T38,E7.1)
1011  FORMAT (//" Random Phase Approximation : ",
     &        /" Calculate",I3," excitation energies for irrep ",A4)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCDINI
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Initialize the dump file and status information.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "results.inc"
      INCLUDE "inpt.inc"
      INCLUDE "symm.inc"
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      DIMENSION BUFF(RSRCLEN)
      LOGICAL TOBE
C
C---------------Executable code--------------------------------------
C
C     Check the status of the restart file
C
      IRECL = IRECLEN (0,RSRCLEN,0)
      INQUIRE (FILE='MCCRES',EXIST=TOBE)
      IF (TOBE) THEN
         OPEN (MCCRES,FILE='MCCRES',ACCESS='DIRECT',STATUS='OLD',
     &         RECL=IRECL)
C
C         Read the input information from the previous run
C
C         Not implemented yet
C
C         Read the status information
C
          READ (MCCRES,REC=2) IMSTAT
C
C         Read the energies (read the whole common block)
C
          READ (MCCRES,REC=3) BUFF
          CALL DCOPY(7,BUFF,1,ETOT,1)
      ELSE
         OPEN (MCCRES,FILE='MCCRES',ACCESS='DIRECT',STATUS='NEW',
     &         RECL=IRECL)
         CALL DCOPY(7,A0,0,ETOT,1)
         DO I = 1, NSTATS
            IMSTAT(I) = 0
         ENDDO
      ENDIF
C
C     Initialize the status information
C     When restarting the previous status is retained if the module
C     was on the menu.
C
      IF (DOSORT.AND.IMSTAT(1).EQ.0)             IMSTAT(1) = 1
      IF (DOSORT.AND.SRTALL.AND.IMSTAT(2).EQ.0)  IMSTAT(2) = 1
      IF (DOSORT.AND.IMSTAT(3).EQ.0)             IMSTAT(3) = 1
      IF (DOENER.AND.DOMP2.AND.IMSTAT(1).EQ.0)   IMSTAT(4) = 1
      IF (DOENER.AND.DOCCSD.AND.IMSTAT(1).EQ.0)  IMSTAT(5) = 1
      IF (DOENER.AND.DOCCSDT.AND.IMSTAT(1).EQ.0) IMSTAT(6) = 1
C
C     Initialize the record information
C
      IRECIN = 1 ! Input for the current run
      IRECST = 2 ! Status of the current run
      IRECEN = 3 ! Energies and other results
      IRECT1 = 4 ! T1 amplitudes
      NRECT1 = INTOWP(NDIMT1*RCW) / RSRCLEN + 1
      IRECT2 = IRECT1 + NRECT1 ! T2 amplitudes
      NRECT2 = INTOWP(NDIMT2*RCW) / RSRCLEN + 1
C
      CLOSE (MCCRES,STATUS='KEEP')
C
      CALL CCDUMP
C
      RETURN
 1000 FORMAT (3X,'Step : ',A,' done in previous run; will be skipped')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCDUMP
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Dump all the information in the input and restart block
C     to the begin of the restart file
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "results.inc"
      INCLUDE "inpt.inc"
C
C---------------Local variables--------------------------------------
C
      DIMENSION BUFF(RSRCLEN)
C
C---------------Executable code--------------------------------------
C
      IRECL = IRECLEN (0,RSRCLEN,0)
      OPEN (MCCRES,FILE='MCCRES',ACCESS='DIRECT',STATUS='OLD',
     &      RECL=IRECL)
C
C     Dump the input information
C
C     Not implemented yet
C
C     Dump the status information
C
      WRITE (MCCRES,REC=2) IMSTAT
C
C     Dump the energies (write the whole common block)
C
      CALL DCOPY(7,ETOT,1,BUFF,1)
      WRITE (MCCRES,REC=3) BUFF
C
      CLOSE (MCCRES,STATUS='KEEP')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCSUMM
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Write summary of the results, based on the restart file
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "results.inc"
      INCLUDE "inpt.inc"
C
C---------------Local variables--------------------------------------
C
      DIMENSION BUFF(RSRCLEN)
      CHARACTER DATEX*10,TIMEX*8
      CHARACTER*30 STATCH(0:4)
      DATA STATCH/'Never asked for              ',
     &            'On menu but not yet called   ',
     &            'Started, restartable         ',
     &            'Completed, restartable       ',
     &            'Completed, not restartable   '/
C
C---------------Executable code--------------------------------------
C
      IRECL = IRECLEN (0,RSRCLEN,0)
      OPEN (MCCRES,FILE='MCCRES',ACCESS='DIRECT',STATUS='OLD',
     &      RECL=IRECL)
C
C     Read the status information
C
      READ (MCCRES,REC=2) IMSTAT
C
C     Read the energies (read the whole common block)
C
      READ (MCCRES,REC=3) BUFF
      CALL DCOPY(7,BUFF,1,ETOT,1)
C
      CLOSE (MCCRES,STATUS='KEEP')
C
      CALL DAYTIME (DATEX,TIMEX)
C
      WRITE (IW,1000) ('-',I=1,80)
C
C     Write date & time
C
      WRITE (IW,1001) DATEX,TIMEX
C
C     Write status information
C
      WRITE (IW,1002)
C
      WRITE (IW,1021) STATCH(IMSTAT(1))
      WRITE (IW,1022) STATCH(IMSTAT(2))
      WRITE (IW,1023) STATCH(IMSTAT(3))
      WRITE (IW,1024) STATCH(IMSTAT(4))
      WRITE (IW,1025) STATCH(IMSTAT(5))
      WRITE (IW,1026) STATCH(IMSTAT(6))
C
C     Write energies
C
      WRITE (IW,1003) 
      WRITE (IW,1070) ESCF
      WRITE (IW,1072) EMP2
      WRITE (IW,1073) ECCSD
      WRITE (IW,1074) ET1
      WRITE (IW,1075) ET2
      WRITE (IW,1076) ET3
      WRITE (IW,1080) ESCF+EMP2
      WRITE (IW,1081) ESCF+ECCSD
      WRITE (IW,1082) ESCF+ECCSD+ET1
      WRITE (IW,1083) ESCF+ECCSD+ET1+ET2
      WRITE (IW,1084) ESCF+ECCSD+ET1+ET3
C
      WRITE (IW,1000) ('-',I=1,80)
C
 1000 FORMAT (//80A)
 1001 FORMAT (//" Today is :",T15,A10/" The time is :",T17,A8)
 1002 FORMAT (/' Status of the calculations')
 1003 FORMAT (/' Overview of calculated energies')
 1021 FORMAT (' Integral sort # 1 :',T40,A30)
 1022 FORMAT (' Integral sort # 2 :',T40,A30)
 1023 FORMAT (' Fock matrix build :',T40,A30)
 1024 FORMAT (' MP2 energy calculation :',T40,A30)
 1025 FORMAT (' CCSD energy calculation :',T40,A30)
 1026 FORMAT (' CCSD(T) energy calculation :',T40,A30)
 1070 FORMAT('@ SCF energy :',T40,F25.15)
 1072 FORMAT('@ MP2 correlation energy :',T40,F25.15)
 1073 FORMAT('@ CCSD correlation energy :',T40,F25.15)
 1074 FORMAT('@ 4th order triples correction :',T40,F25.15)
 1075 FORMAT('@ 5th order triples (T) correction :',T40,F25.15)
 1076 FORMAT('@ 5th order triples -T  correction :',T40,F25.15)
 1080 FORMAT('@ Total MP2 energy :',T40,F25.15)
 1081 FORMAT('@ Total CCSD energy :',T40,F25.15)
 1082 FORMAT('@ Total CCSD+T  energy :',T40,F25.15)
 1083 FORMAT('@ Total CCSD(T) energy :',T40,F25.15)
 1084 FORMAT('@ Total CCSD-T  energy :',T40,F25.15)
C
      RETURN
      END
