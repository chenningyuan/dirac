C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : ccgenp.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BLANKR (N,A)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Circular shift of character array such that the blanks are to the
C     right side of the array.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N
      CHARACTER*1 A(N)
C
C---------------Common Blocks--------------------------------------
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IS = 0
      DO I = N, 1, -1
         IF (A(I).NE.' ') GOTO 2
         IS = IS + 1
      ENDDO
C
    2 IF (IS.EQ.0) RETURN
C
      DO I = N-IS, 1, -1
            A(I+IS) = A(I)
      ENDDO
      DO I = 1, IS
         A(I) = ' '
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CONJUGA (N,X,INCX)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Take complex conjugate of vector X
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N,INCX
      COMPLEX*16 X(*)
C
C---------------Common Blocks--------------------------------------
C
C---------------Local variables--------------------------------------
C
      INTEGER II
C
C---------------Executable code--------------------------------------
C
      II = 1
      DO I = 1, N
         X(II) = DCONJG(X(II))
C64B     X(II) =  CONJG(X(II))
         II = II + INCX
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      Subroutine CpuUsed (Sec)
C
      Real*8 Sec
C
#if defined (SYS_AIX)  || defined (SYS_IRIX) || defined (SYS_LINUX) \
 || defined (SYS_HPUX)
      Sec = MCLOCK()/100.D0
#elif defined (SYS_CRAY) || defined (SYS_T3D)
      Sec = Second()
#else
      Does MCLOCK exist on your computer???
#endif
C
      Return
      End
C
      Real*8 Function CpuDelta()
C
      Real*8 Sec1,Sec0
      Common/CPU/Sec0,Sec1
      Data Sec0 /0.d0/
C
      Call CPUused(Sec1)
      CPUdelta = Sec1 - Sec0
      Sec0 = Sec1
C
      Return
      End
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      Subroutine DayTime (Datum,Tijd)
#if defined (SYS_HPUX)
      CHARACTER FDATE*24
#endif
      Character*10 Datum,Tijd*8,String*26
      external time
      Datum = '          '
      Tijd  = '        '
#if defined (SYS_CRAY) || defined (SYS_T3D)
      Call Date(Datum)
      Call Clock(Tijd)
#elif defined (SYS_HPUX)
      STRING(1:24) = FDATE()
      DATUM = STRING(1:10)
      TIJD  = STRING(12:19)
#elif defined (SYS_IRIX) || defined (SYS_LINUX) || defined (SYS_HPUX)
      Call Date(Datum)
      Call Time(Tijd)
#elif defined (SYS_AIX) 
      Call fdate_ (String)
      Datum = Strinf(1:10)
      Tijd  = String(12:19)
#else
      Date and time on your computer???
#endif
      Return
      End
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      INTEGER FUNCTION INTOWP(N)
#if defined (SYS_CRAY) || defined (SYS_T3D)
      INTOWP = N
#else
      INTOWP=2*N
#endif
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      INTEGER FUNCTION SEC2I(N)
C
C     BUFF in the common block IOBUFF should have a length of
C     at least two times the value used here.
C
      SEC2I=512*N
C64B  SEC2I=512*N
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RSAVE(ITAPE,BLOCKS)
C
C FLUSH THE FILE
C
      IMPLICIT INTEGER (A-Z)
C
      JCODE = 3
      CALL RCLOSE (ITAPE,JCODE)
      CALL RFILE  (ITAPE,BLOCKS)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RFILE(ITAPE,BLOCKS)
C
C USE THIS FOR ALL UNFORMATTED FILES
C
      IMPLICIT INTEGER (A-Z)
C
      CHARACTER*4 DN,FNAM*8
C
      COMMON /POINTR/ POS(100)
      COMMON /SECT/ SECTOR,MAXLN
      COMMON/ FLNAME/ DN(100) 
C
      IF (ITAPE.GT.99.OR.ITAPE.LT.1) THEN
      WRITE(6,*) ' INVALID TAPE NUMBER GIVEN: ITAPE=' ,ITAPE
      CALL QUIT('ccgenp.F line 184')
      ENDIF
      SECTOR = SEC2I(1)
      MAXLN = SECTOR*1024*1024
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IF (BLOCKS.EQ.0) BLOCKS = 20
      IF (BLOCKS.EQ.1) BLOCKS = 2
      ISTATS = 1
      CALL WOPEN(ITAPE,BLOCKS,ISTATS,IERR)
      IF(IERR .NE. 0) THEN
          WRITE(6,*) ' PROBLEMS OPENING FILE ',ITAPE,' IERR =',IERR
          CALL QUIT('ccgenp line 195')
      END IF
#else
      FNAM = DN(ITAPE)
      FNAM(5:5)=CHAR(0)
      SIZE = 0
      STATUS = 1
      CALL OPENC(ITAPE,FNAM,SIZE,STATUS)
#endif
C
C PLACE THE POINTER OF ITAPE AT THE BEGINNING
C
      POS(ITAPE) = 1
C
C RETURN
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RGETSA(ITAPE,IADR)
C
C THIS ROUTINE RETURNS THE CURRENT SECTOR ADDRESS FOR FILES USING
C THE SREAD/SWRIT AND/OR RREAD/RWRIT IO ROUTINES.
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /POINTR/ POS(100)
      COMMON /SECT/ SECTOR,MAXLN
C
      IPOS = POS(ITAPE)
      IREC = (IPOS-1)/SECTOR + 1
      TEST = SECTOR*(IREC-1) + 1
      IF(IPOS.NE.TEST) THEN
      WRITE(*,*) ' ERROR ENCOUNTERED IN RGETSA FOR FILE ',ITAPE
      WRITE(*,*) ' IPOS,TEST = ',IPOS,TEST
      CALL QUIT('ccgenp.F line  231')
      END IF
C
      IADR=IREC
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RSETSA(ITAPE,IADR)
C
C THIS ROUTINE SETS THE POINTER OF A UTILITIES FILE TO THE DESIRED
C SECTOR ADDRESS.
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /POINTR/ POS(100)
      COMMON /SECT/ SECTOR,MAXLN
C
      IPOS = SEC2I(IADR-1) + 1
      POS(ITAPE) = IPOS
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WWRITW(ITAPE,ARRAY,NLEN,FWORD,LWORD)
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /POINTR/ POS(100)
      COMMON /SECT/ SECTOR,MAXLN
C
      DIMENSION ARRAY(NLEN)
C
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL PUTWA(ITAPE,ARRAY,FWORD,NLEN,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN WWRITW FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 269')
      END IF
#else
      OFF = FWORD - 1
      CALL WRABSF(ITAPE,ARRAY,NLEN,OFF)
#endif
C
C UPDATE LWORD AND THEN RETURN
C
      LWORD = FWORD + NLEN
      POS(ITAPE) = LWORD
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WREADW(ITAPE,ARRAY,NLEN,FWORD,LWORD)
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /POINTR/ POS(100)
      COMMON /SECT/ SECTOR,MAXLN
C
      DIMENSION ARRAY(NLEN)
C
#if defined (SYS_CRAY) || defined(SYS_T3D)
      CALL GETWA(ITAPE,ARRAY,FWORD,NLEN,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN WREADW FOR FILE ',ITAPE,FWORD,NLEN
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 298')
      END IF
#else
      OFF = FWORD - 1
      CALL WRABSF(ITAPE,ARRAY,NLEN,OFF)
#endif
C
C UPDATE LWORD AND THEN RETURN
C
      LWORD = FWORD + NLEN
      POS(ITAPE) = LWORD
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RREAD(ITAPE,ARRAY,NLEN,IREC)
C
C THIS ROUTINE READS NLEN INTEGER WORDS INTO ARRAY STARTING AT
C SECTOR IREC.
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /POINTR/ POS(100)
      COMMON /SECT/ SECTOR,MAXLN
      COMMON /IOBUFF/ BUFF(1024)
C
      DIMENSION ARRAY(NLEN)
C
      IF (NLEN.EQ.0) RETURN
      FWORD = (IREC-1)*SECTOR + 1
      NGRP = NLEN/MAXLN
      IF(MAXLN*NGRP.NE.NLEN) NGRP = NGRP + 1
C
      IF(NGRP.EQ.1) THEN
C
      NW = (NLEN/SECTOR)*SECTOR
      OFF = FWORD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IERR = 0
      IF(NW.GT.0) CALL GETWA(ITAPE,ARRAY,FWORD,NW,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN RREAD FOR FILE ',ITAPE,FWORD,NW
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 341')
      END IF
#else
      IF(NW.GT.0) CALL RDABSF(ITAPE,ARRAY,NW,OFF)
#endif
C
      EXTRA = MOD(NLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
      FWD = FWORD + NW
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL GETWA(ITAPE,BUFF,FWD,SECTOR,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN RREAD FOR FILE ',ITAPE,FWD,NW
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 356')
      END IF
#else
      CALL RDABSF(ITAPE,BUFF,SECTOR,OFF)
#endif
      DO 12 I = 1,EXTRA
      ARRAY(NW+I) = BUFF(I)
   12 CONTINUE
      END IF
C
      ELSE
C
      LLEN = MOD(NLEN,MAXLN)
      IF(LLEN.EQ.0) LLEN = MAXLN
      FWD = FWORD
C
      DO 11 IGRP = 1,NGRP-1
      IPT = (IGRP-1)*MAXLN + 1
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL GETWA(ITAPE,ARRAY(IPT),FWD,MAXLN,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN RREAD FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 380')
      END IF
#else
      CALL RDABSF(ITAPE,ARRAY(IPT),MAXLN,OFF)
#endif
      FWD = FWD + MAXLN
  11  CONTINUE
C
      NW = (LLEN/SECTOR)*SECTOR
      IPT = (NGRP-1)*MAXLN + 1
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IF(NW.GT.0) CALL GETWA(ITAPE,ARRAY(IPT),FWD,NW,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN RREAD FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 396')
      END IF
#else
      IF(NW.GT.0) CALL RDABSF(ITAPE,ARRAY(IPT),NW,OFF)
#endif
C
      EXTRA = MOD(LLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
      FWD = FWD + NW
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL GETWA(ITAPE,BUFF,FWD,SECTOR,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN RREAD FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 411')
      END IF
#else
      CALL RDABSF(ITAPE,BUFF,SECTOR,OFF)
#endif
      IPT = IPT + NW - 1
      DO 13 I = 1,EXTRA
      ARRAY(IPT+I) = BUFF(I)
   13 CONTINUE
      END IF
C
      END IF      
C
C WE MUST POSITION THE FILE POINTER AT THE BEGINNING OF THE NEXT SECTOR
C
      TEST = NLEN/SECTOR
      IF(SECTOR*TEST .NE. NLEN) TEST = TEST + 1
      TESTL = SECTOR*TEST
      POS(ITAPE) = FWORD + TESTL
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RWRIT(ITAPE,ARRAY,NLEN,IREC)
C
C THIS ROUTINE WRITES NLEN INTEGER WORDS FROM ARRAY TO FILE ITAPE
C STARTING AT SECTOR IREC.
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /POINTR/ POS(100)
      COMMON /SECT/ SECTOR,MAXLN
      COMMON /IOBUFF/ BUFF(1024)
C
      DIMENSION ARRAY(NLEN)
C
      IF (NLEN.EQ.0) RETURN
      FWORD = (IREC-1)*SECTOR + 1
      NGRP = NLEN/MAXLN
      IF(MAXLN*NGRP.NE.NLEN) NGRP = NGRP + 1
C
      IF(NGRP.EQ.1) THEN
      NW = (NLEN/SECTOR)*SECTOR
      OFF = FWORD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IERR = 0
      IF(NW.GT.0) CALL PUTWA(ITAPE,ARRAY,FWORD,NW,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN RWRIT FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 461')
      END IF
#else
      CALL WRABSF(ITAPE,ARRAY,NW,OFF)
#endif
C
      EXTRA = MOD(NLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
      DO 12 I = 1,EXTRA
      BUFF(I) = ARRAY(NW+I)
   12 CONTINUE
      FWD = FWORD + NW
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL PUTWA(ITAPE,BUFF,FWD,SECTOR,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN RWRIT FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 479')
      END IF
#else
      CALL WRABSF(ITAPE,BUFF,SECTOR,OFF)
#endif
      END IF
C
      ELSE
C
      LLEN = MOD(NLEN,MAXLN)
      IF(LLEN.EQ.0) LLEN = MAXLN
      FWD = FWORD
C
      DO 11 IGRP = 1,NGRP-1
      IPT = (IGRP-1)*MAXLN + 1
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL PUTWA(ITAPE,ARRAY(IPT),FWD,MAXLN,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN RWRIT FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 500')
      END IF
#else
      CALL WRABSF(ITAPE,ARRAY(IPT),MAXLN,OFF)
#endif
      FWD = FWD + MAXLN
  11  CONTINUE
C
      NW = (LLEN/SECTOR)*SECTOR
      IPT = (NGRP-1)*MAXLN + 1
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IF(NW.GT.0) CALL PUTWA(ITAPE,ARRAY(IPT),FWD,NW,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN RWRIT FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 516')
      END IF
#else
      IF(NW.GT.0) CALL WRABSF(ITAPE,ARRAY(IPT),NW,OFF)
#endif
C
      EXTRA = MOD(LLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
      IPT = IPT + NW - 1
      DO 13 I = 1,EXTRA
      BUFF(I) = ARRAY(IPT+I)
   13 CONTINUE
      FWD = FWD + NW
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL PUTWA(ITAPE,BUFF,FWD,SECTOR,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN RWRIT FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 535')
      END IF
#else
      CALL WRABSF(ITAPE,BUFF,SECTOR,OFF)
#endif
      END IF
C
      END IF      
C
C WE MUST POSITION THE FILE POINTER AT THE BEGINNING OF THE NEXT SECTOR
C
      TEST = NLEN/SECTOR
      IF(SECTOR*TEST .NE. NLEN) TEST = TEST + 1
      TESTL = SECTOR*TEST
      POS(ITAPE) = FWORD + TESTL
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SREAD(ITAPE,ARRAY,NLEN)
C
C THIS ROUTINE READS NLEN INTEGER WORDS INTO ARRAY STARTING AT
C STARTING AT THE CURRENT POINTER LOCATION.
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /POINTR/ POS(100)
      COMMON /SECT/ SECTOR,MAXLN
      COMMON /IOBUFF/ BUFF(1024)
C
      DIMENSION ARRAY(NLEN)
C
      FWORD = POS(ITAPE)
      NGRP = NLEN/MAXLN
      IF(MAXLN*NGRP.NE.NLEN) NGRP = NGRP + 1
C
      IF(NGRP.EQ.1) THEN
C
      NW = (NLEN/SECTOR)*SECTOR
      OFF = FWORD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IERR = 0
      IF(NW.GT.0) CALL GETWA(ITAPE,ARRAY,FWORD,NW,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN SREAD FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 581')
      END IF
#else
      IF(NW.GT.0) CALL RDABSF(ITAPE,ARRAY,NW,OFF)
#endif
C
      EXTRA = MOD(NLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
      FWD = FWORD + NW
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL GETWA(ITAPE,BUFF,FWD,SECTOR,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN SREAD FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 596')
      END IF
#else
      CALL RDABSF(ITAPE,BUFF,SECTOR,OFF)
#endif
      DO 12 I = 1,EXTRA
      ARRAY(NW+I) = BUFF(I)
   12 CONTINUE
      END IF
C
      ELSE
C
      LLEN = MOD(NLEN,MAXLN)
      IF(LLEN.EQ.0) LLEN = MAXLN
      FWD = FWORD
C
      DO 11 IGRP = 1,NGRP-1
      IPT = (IGRP-1)*MAXLN + 1
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL GETWA(ITAPE,ARRAY(IPT),FWD,MAXLN,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN SREAD FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 620')
      END IF
#else
      CALL RDABSF(ITAPE,ARRAY(IPT),MAXLN,OFF)
#endif
      FWD = FWD + MAXLN
  11  CONTINUE
C
      NW = (LLEN/SECTOR)*SECTOR
      IPT = (NGRP-1)*MAXLN + 1
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IF(NW.GT.0) CALL GETWA(ITAPE,ARRAY(IPT),FWD,NW,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN SREAD FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 636')
      END IF
#else
      IF(NW.GT.0) CALL RDABSF(ITAPE,ARRAY(IPT),NW,OFF)
#endif
C
      EXTRA = MOD(LLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
      FWD = FWD + NW
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL GETWA(ITAPE,BUFF,FWD,SECTOR,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN SREAD FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 651')
      END IF
#else
      CALL RDABSF(ITAPE,BUFF,SECTOR,OFF)
#endif
      IPT = IPT + NW - 1
      DO 13 I = 1,EXTRA
      ARRAY(IPT+I) = BUFF(I)
   13 CONTINUE
      END IF
C
      END IF      
C
C WE MUST POSITION THE FILE POINTER AT THE BEGINNING OF THE NEXT SECTOR
C
      TEST = NLEN/SECTOR
      IF(SECTOR*TEST .NE. NLEN) TEST = TEST + 1
      TESTL = SECTOR*TEST
      POS(ITAPE) = FWORD + TESTL
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SWRIT(ITAPE,ARRAY,NLEN)
C
C THIS ROUTINE WRITES NLEN INTEGER WORDS FROM ARRAY TO FILE ITAPE
C STARTING AT THE CURRENT POINTER LOCATION.
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /POINTR/ POS(100)
      COMMON /SECT/ SECTOR,MAXLN
      COMMON /IOBUFF/ BUFF(1024)
C
      DIMENSION ARRAY(NLEN)
C
      FWORD = POS(ITAPE)
      NGRP = NLEN/MAXLN
      IF(MAXLN*NGRP.NE.NLEN) NGRP = NGRP + 1
C
      IF(NGRP.EQ.1) THEN
      NW = (NLEN/SECTOR)*SECTOR
      OFF = FWORD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IERR = 0
      IF(NW.GT.0) CALL PUTWA(ITAPE,ARRAY,FWORD,NW,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN SWRIT FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 700')
      END IF
#else
      IF(NW.GT.0) CALL WRABSF(ITAPE,ARRAY,NW,OFF)
#endif
C
      EXTRA = MOD(NLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
      DO 12 I = 1,EXTRA
      BUFF(I) = ARRAY(NW+I)
   12 CONTINUE
      FWD = FWORD + NW
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL PUTWA(ITAPE,BUFF,FWD,SECTOR,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN SWRIT FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 718')
      END IF
#else
      CALL WRABSF(ITAPE,BUFF,SECTOR,OFF)
#endif
      END IF
C
      ELSE
C
      LLEN = MOD(NLEN,MAXLN)
      IF(LLEN.EQ.0) LLEN = MAXLN
      FWD = FWORD
C
      DO 11 IGRP = 1,NGRP-1
      IPT = (IGRP-1)*MAXLN + 1
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL PUTWA(ITAPE,ARRAY(IPT),FWD,MAXLN,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN SWRIT FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 739')
      END IF
#else
      CALL WRABSF(ITAPE,ARRAY(IPT),MAXLN,OFF)
#endif
      FWD = FWD + MAXLN
  11  CONTINUE
C
      NW = (LLEN/SECTOR)*SECTOR
      IPT = (NGRP-1)*MAXLN + 1
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IF(NW.GT.0) CALL PUTWA(ITAPE,ARRAY(IPT),FWD,NW,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN SWRIT FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 755')
      END IF
#else
      IF(NW.GT.0) CALL WRABSF(ITAPE,ARRAY(IPT),NW,OFF)
#endif
C
      EXTRA = MOD(LLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
      IPT = IPT + NW - 1
      DO 13 I = 1,EXTRA
      BUFF(I) = ARRAY(IPT+I)
   13 CONTINUE
      FWD = FWD + NW
      OFF = FWD - 1
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CALL PUTWA(ITAPE,BUFF,FWD,SECTOR,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR IN SWRIT FOR FILE ',ITAPE
         WRITE(*,*) ' IERR = ',IERR
         CALL QUIT('ccgenp.F line 774')
      END IF
#else
      CALL WRABSF(ITAPE,BUFF,SECTOR,OFF)
#endif
      END IF
C
      END IF      
C
C WE MUST POSITION THE FILE POINTER AT THE BEGINNING OF THE NEXT SECTOR
C
      TEST = NLEN/SECTOR
      IF(SECTOR*TEST .NE. NLEN) TEST = TEST + 1
      TESTL = SECTOR*TEST
      POS(ITAPE) = FWORD + TESTL
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SREW(ITAPE)
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /POINTR/ POS(100)
C
C PLACE THE POINTER OF ITAPE AT THE BEGINNING
C
C
C UPDATE THE POSITION ARRAY
C
      POS(ITAPE) = 1
C
C RETURN
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RCLOSE(ITAPE,JCODE)
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*8 FNAM
      CHARACTER*4 DN
C
      COMMON /POINTR/ POS(100)
      common /flname/ dn(100)
C
C JCODE = 4     CLOSE AND DELETE FILE
C JCODE = 3     CLOSE AND SAVE FILE
C
      IF(JCODE.NE.3.AND.JCODE.NE.4) THEN
         WRITE(*,*) ' INVALID JCODE IN RCLOSE,  JCODE = ',JCODE
         WRITE(*,*) ' FILE ',ITAPE,'  CLOSED AND SAVED.'
         JCODE = 3
      END IF
C
      IF(ITAPE.EQ.6) CALL QUIT(' YOU CANNOT CLOSE A FILE ON UNIT 6')
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IERR = 0
      CALL WCLOSE(ITAPE,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,*) ' ERROR ENCOUNTERED CLOSING FILE',ITAPE
         WRITE(*,*) ' IERR,JCODE = ',IERR,JCODE
      END IF
#else
      CALL CLOSEC(ITAPE)
#endif
C
      IF(JCODE.EQ.4) THEN
        FNAM = DN(ITAPE)
        FNAM(5:5)=CHAR(0)
#if defined (SYS_CRAY) || defined (SYS_T3D)
        IF (ITAPE.LT.10) WRITE (FNAM,'(A5,I1,A1)') 'fort.',ITAPE,char(0)
        IF (ITAPE.GT.10) WRITE (FNAM,'(A5,I2,A1)') 'fort.',ITAPE,char(0)
        CALL PXFUNLINK(FNAM,0,STAT)
#else
        STAT = UNLINK(FNAM)
#endif
        IF(STAT.NE.0) THEN
          WRITE(6,*) '  PROBLEMS DELETING FILE ',ITAPE,'  STAT = ',STAT
          WRITE(6,*) ' DN = ',DN(ITAPE)
        END IF
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      BLOCK DATA FLNAMES
      CHARACTER*4 DN
      COMMON/ FLNAME/ DN(100) 
      DATA DN /'ft1 ','ft2 ','ft3 ','ft4 ','ft5 ','ft6 ','ft7 ','ft8 ',
     &  'ft9 ',
     &  'ft10','ft11','ft12','ft13','ft14','ft15','ft16','ft17','ft18',
     &  'ft19','ft20','ft21','ft22','ft23','ft24','ft25','ft26','ft27',
     &  'ft28','ft29','ft30','ft31','ft32','ft33','ft34','ft35','ft36',
     &  'ft37','ft38','ft39','ft40','ft41','ft42','ft43','ft44','ft45',
     &  'ft46','ft47','ft48','ft49','ft50','ft51','ft52','ft53','ft54',
     &  'ft55','ft56','ft57','ft58','ft59','ft60','ft61','ft62','ft63',
     &  'ft64','ft65','ft66','ft67','ft68','ft69','ft70','ft71','ft72',
     &  'ft73','ft74','ft75','ft76','ft77','ft78','ft79','ft80','ft81',
     &  'ft82','ft83','ft84','ft85','ft86','ft87','ft88','ft89','ft90',
     &  'ft91','ft92','ft93','ft94','ft95','ft96','ft97','ft98','ft99',
     &  'help'/
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FLIN(A,IDIM,IN,IM,DET)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     LINEAR SIMULTANEOUS EQUATION
C
C     A(IN*IN) * X(IN*IM) = B(IN*IM)
C
C     A & B SHOULD BE STORED ON A(IN*(IN+IM))
C     SOLUTION X WILL BE STORED ON B PART IN DIMENSION A.
C
      DIMENSION A(IDIM,*)
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C64B  PARAMETER (ZERO=0.0,ONE=1.0)
C
      N=IN
      NR=IM
      JMAX=N+NR
      SIGN=ONE
C M IS THE STAGE OF ELIMINATION
      DO 49 M=1,N
      TEMP=ZERO
      DO 41 I=M,N
      IF(M.GT.1)A(I,M)=A(I,M)-DOTT(A(I,1),IDIM,A(1,M),1,M-1)
      AVAL=A(I,M)
      IF(ABS(AVAL).LE.TEMP)GOTO 41
C64B  IF( ABS(AVAL).LE.TEMP)GOTO 41
      TEMP=ABS(AVAL)
C64B  TEMP= ABS(AVAL)
      IMAX=I
 41   CONTINUE
      IF(TEMP.LE.ZERO)GOTO 999
      IF(IMAX.EQ.M)GOTO 45
      SIGN=-SIGN
      DO 44 J=1,JMAX
      STOR=A(M,J)
      A(M,J)=A(IMAX,J)
      A(IMAX,J)=STOR
 44   CONTINUE
 45   CONTINUE
      JJ=M+1
      IF(JJ.GT.JMAX)GOTO 49
      IF(M.GT.1)GOTO 47
      DO 46 J=JJ,JMAX
      A(1,J)=A(1,J)/A(1,1)
 46   CONTINUE
      D=A(1,1)
      GOTO 49
 47   CONTINUE
      DO 48 J=JJ,JMAX
      A(M,J)=(A(M,J)-DOTT(A(M,1),IDIM,A(1,J),1,M-1))/A(M,M)
 48   CONTINUE
      D=D*A(M,M)
 49   CONTINUE
      IF(NR.EQ.0) RETURN
      DO 59 I=1,NR
      NPI=N+I
      DO 58 K=2,N
      J=N+1-K
      A(J,NPI)=A(J,NPI)-DOTT(A(J,J+1),IDIM,A(J+1,NPI),1,K-1)
 58   CONTINUE
 59   CONTINUE
C***  IF(ABS(D).GE.1.0D+36) D=1.0D+36
C***  IF(ABS(D).LE.1.0D-36) D=1.0D-36
      DET=D*SIGN
      RETURN
C ON ZERO PIVOT, SET DET=0.AND RETURN TO CALLING PROGRAM NOV 1972
 999  DET=ZERO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION DOTT(A,NA,B,NB,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1)
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C64B  PARAMETER (ZERO=0.0,ONE=1.0)
C
      IAPT=1
      IBPT=1
      D=ZERO
      DO 10 I=1,N
      D=D+A(IAPT)*B(IBPT)
      IAPT=IAPT+NA
      IBPT=IBPT+NB
 10   CONTINUE
      DOTT=D
      RETURN
      END
