C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : ccints.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RDINTS(IPRNT,NSP,INDEX,KR,INDK,INDL,CBUF,IBUF,NVT,
     &                  BUF1,BUF2,IREPSP,IREPSPI,
     &                  VOOOO,VOOOOA,VVOOO,VVOOOA,VVVOO,VVVOOA,
     &                  VVOVO,VVOVOA,VVOVV,VVOVVA,VVVVV,VVVVVA,
     &                  EPS,FOO,FVO,FVV,BREIT,SRTALL,INTERFACE,
     &                  MREMAIN)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Reads 1-electron and 2-electron integrals
C     Constructs Fock matrix elements
C     This version reads from files created by TMOONE & ROTRAN or DIRAC
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NVT,NSP
      INTEGER KR(-NSP:NSP),INDK(*),INDL(*)
      INTEGER INDEX(NSP,NSP,4),IBUF(16,6,*)
      INTEGER IREPSP,IREPSPI
      REAL*8 FOO(*),FVO(*),FVV(*)
      REAL*8 VOOOO(RCW,*),VOOOOA(RCW,*),VVOOO(RCW,*),VVOOOA(RCW,*)
      REAL*8 VVVOO(RCW,*),VVVOOA(RCW,*),VVOVO(RCW,*),VVOVOA(RCW,*)
      REAL*8 VVOVV(RCW,*),VVOVVA(RCW,*),VVVVV(RCW,*),VVVVVA(RCW,*)
      REAL*8 BUF1(*),BUF2(*),CBUF(RCW,*)
      REAL*8 EPS(*)
      LOGICAL BREIT,SRTALL
      CHARACTER*10 INTERFACE
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
      INCLUDE "results.inc"
      INCLUDE "files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 REPNUC,E0,E1,EPSDIF
      CHARACTER*10 DATEX,TIMEX*8
      INTEGER IDUM
      LOGICAL LDUM
C
C---------------Executable code--------------------------------------
C
      call memchk ('Entering rdints',foo(0),1)
      CALL MKINDEX (NSP,IREPSP,IREPSPI,INDEX)
C
      NSPQ = NSP*NSP
C
C  Check buffer size for passing of VOVV and VVVV integrals
C
      IF (MOD(NVT*IRW*RCW,SEC2I(1)).NE.0)
     &    CALL QUIT("INVALID BUFFER SIZE IN RDINTS")
C
C  GET THE NUCLEAR REPULSION ENERGY
C
      OPEN (MRCONEE,FILE='MRCONEE',FORM='UNFORMATTED')
      READ (MRCONEE) IDUM,LDUM,REPNUC
C
C  READ ONE ELECTRON INTEGRALS
C
      READ (MRCONEE) 
      READ (MRCONEE) 
      READ (MRCONEE) 
      READ (MRCONEE) 
      READ (MRCONEE) (BUF1(IJ),BUF2(IJ),IJ=1,NSPQ)
      CLOSE (MRCONEE)
C  
C  ORDER AND PUT IN FOCK MATRIX
C
      IJ = 0
      DO J = 1, NSP
         DO I = 1, NSP
            IJ = IJ + 1
            ITYP = INDEX(I,J,1)
            IND = INDEX(I,J,2)
            IF (IND.NE.0) THEN
            IF (ITYP.EQ.1) THEN
               IF (CARITH) THEN
                  FOO(RCW*IND-1) = BUF1(IJ)
                  FOO(RCW*IND)   = - BUF2(IJ)
               ELSE
                  FOO(IND) = BUF1(IJ)
               ENDIF
            ELSEIF (ITYP.EQ.2) THEN
               IF (CARITH) THEN
                  FVO(RCW*IND-1) = BUF1(IJ)
                  FVO(RCW*IND)   = - BUF2(IJ)
               ELSE
                  FVO(IND) = BUF1(IJ)
               ENDIF
            ELSEIF (ITYP.EQ.4) THEN
               IF (CARITH) THEN
                  FVV(RCW*IND-1) = BUF1(IJ)
                  FVV(RCW*IND)   = - BUF2(IJ)
               ELSE
                  FVV(IND) = BUF1(IJ)
               ENDIF
            ENDIF
            ENDIF
         ENDDO
      ENDDO
C
C  INITIALIZE READING OF COULOMB INTEGRALS
C
      OPEN(MDINT,FILE='MDCINT',FORM='UNFORMATTED')
      REWIND(MDINT)
      IF (INTERFACE(1:7).EQ.'MOLFDIR') THEN
         READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,
     &    (IDUM,I=1,4*NKR),(KR(I),KR(-I),I=1,NKR)
      ELSEIF (INTERFACE(1:5).EQ.'DIRAC') THEN
         READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,
     &                                    (KR(I),KR(-I),I=1,NKR)
      ELSE
         WRITE (IW,'(//1X,3A)') 'Interface ',INTERFACE,' not defined'
         CALL QUIT('Adapt RDINTS')
      ENDIF
      IF (2*NKR.NE.NSP) 
     &   CALL QUIT("INCONSISTENT MRCONEE AND MDCINT FILES")
      ITRSGN = 1
C
C  Initialize writing to sorted integral files
C  Use extra memory as buffer space
C
      NINTSB0 = MAX0(MREMAIN/512,20)
      NINTSB0 = MIN0((NFOO+NFVO+NFVV)*RCW/512+1,NINTSB0)
      CALL RFILE(ITAPT,NINTSB0)
      MREMAIN = MREMAIN - (NINTSB0-20)*512
C
      NINTSB1 = MAX0(MREMAIN/512,20)
      NINTSB1 = MIN0(NV1*RCW/512+1,NINTSB1)
      CALL RFILE(ITAPT+1,NINTSB1)
      MREMAIN = MREMAIN - (NINTSB1-20)*512
C
      NINTSB2 = MAX0(MREMAIN/512,20)
      NINTSB2 = MIN0(NV2*RCW/512+1,NINTSB2)
      CALL RFILE(ITAPT+2,NINTSB2)
      MREMAIN = MREMAIN - (NINTSB2-20)*512
C
      NINTSB3 = MAX0(MREMAIN/512,20)
      NINTSB3 = MIN0(NV3*RCW/512+1,NINTSB3)
      CALL RFILE(ITAPT+3,NINTSB3)
      MREMAIN = MREMAIN - (NINTSB3-20)*512
C
      NINTSB4 = MAX0(MREMAIN/512,20)
      NINTSB4 = MIN0(NV4*RCW/512+1,NINTSB4)
      CALL RFILE(ITAPT+4,NINTSB4)
      MREMAIN = MREMAIN - (NINTSB4-20)*512
C
      IF (IMSTAT(1).GT.1) THEN
         WRITE (IW,1010)
         GOTO 150
      ENDIF
C
C     RE-ENTRY POINT IN CASE THE GAUNT INTERACTION HAS TO BE ADDED
C
    1 ICOUNT = 0
C
C     GET THE OOOO, VOOO, VVOO AND VOVO INTEGRALS
C
      CALL XCOPY (NV1,A0,0,VOOOO,1)
      CALL XCOPY (NV2,A0,0,VVOOO,1)
      CALL XCOPY (NV3,A0,0,VVVOO,1)
      CALL XCOPY (NV4,A0,0,VVOVO,1)
      CALL XCOPY (NV1,A0,0,VOOOOA,1)
      CALL XCOPY (NV2,A0,0,VVOOOA,1)
      CALL XCOPY (NV3,A0,0,VVVOOA,1)
      CALL XCOPY (NV4,A0,0,VVOVOA,1)
C
  100 CALL GET2INTS (MDINT,KR,NSP,INDEX,ITRSGN,INDK,INDL,
     &               NZ,CBUF,IBUF)
      IF (NZ.LT.0) GOTO 110
      ICOUNT = ICOUNT + NZ
      DO 102 I = 1, NZ
         DO 101 IP = 1, 16
            ITYP = IBUF(IP,1,I)
            IREP = IBUF(IP,2,I) 
            IND1 = IBUF(IP,3,I)
            IND2 = IBUF(IP,4,I)
            ITRS = IBUF(IP,5,I)
            JTRS = IBUF(IP,6,I)
            IF (ITYP.EQ.1) THEN         
               IND = IOOOOTT(IREP) + (IND1-1)*NOOT(IREP) + IND2
               VOOOO(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VOOOO(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.-1) THEN         
               IND = IOOOOTT(IREP) + (IND1-1)*NOOT(IREP) + IND2
               VOOOOA(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VOOOOA(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.2) THEN       
               IND = IVOOOT(IREP) + (IND1-1)*NVO(IREP) + IND2
               VVOOO(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VVOOO(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.-2) THEN
               IND = IVOOOT(IREP) + (IND1-1)*NVO(IREP) + IND2
               VVOOOA(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VVOOOA(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.4) THEN
               IND = IVVOOTT(IREP) + (IND1-1)*NVVT(IREP) + IND2
               VVVOO(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VVVOO(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.-4) THEN
               IND = IVVOOTT(IREP) + (IND1-1)*NVVT(IREP) + IND2
               VVVOOA(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VVVOOA(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.6) THEN
               IND = IVOVO(IREP) + (IND1-1)*NVO(IREP) + IND2
               VVOVO(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VVOVO(2,IND) = CBUF(2,I) * JTRS
            ELSEIF (ITYP.EQ.-6) THEN
               IND = IVOVO(IREP) + (IND1-1)*NVO(IREP) + IND2
               VVOVOA(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) VVOVOA(2,IND) = CBUF(2,I) * JTRS
            ENDIF
  101    CONTINUE
  102 CONTINUE
      GOTO 100
C
C     ADD THE DIRECT AND EXCHANGE INTEGRALS
C
  110 CALL XAXPY (NV1,-A1,VOOOOA,1,VOOOO,1)
      CALL XAXPY (NV2,-A1,VVOOOA,1,VVOOO,1)
      CALL XAXPY (NV3,-A1,VVVOOA,1,VVVOO,1)
      CALL XAXPY (NV4,-A1,VVOVOA,1,VVOVO,1)
C
C     ADD BREIT AND COULOMB INTEGRALS
C
      IF (BREIT.AND.ITRSGN.LT.0) THEN
         CALL RREAD(ITAPT+1,VOOOOA,NV1*RCW*IRW,1)
         CALL RREAD(ITAPT+2,VVOOOA,NV2*RCW*IRW,1)
         CALL RREAD(ITAPT+3,VVVOOA,NV3*RCW*IRW,1)
         CALL RREAD(ITAPT+4,VVOVOA,NV4*RCW*IRW,1)
         CALL XAXPY(NV1,A1,VOOOOA,1,VOOOO,1)
         CALL XAXPY(NV2,A1,VVOOOA,1,VVOOO,1)
         CALL XAXPY(NV3,A1,VVVOOA,1,VVVOO,1)
         CALL XAXPY(NV4,A1,VVOVOA,1,VVOVO,1)
      ENDIF
      CALL RWRIT(ITAPT+1,VOOOO,NV1*RCW*IRW,1)
      CALL RWRIT(ITAPT+2,VVOOO,NV2*RCW*IRW,1)
      CALL RWRIT(ITAPT+3,VVVOO,NV3*RCW*IRW,1)
      CALL RWRIT(ITAPT+4,VVOVO,NV4*RCW*IRW,1)
C
C     This is a restart point : save the files and update status.
C
  150 CONTINUE
      CALL RSAVE(ITAPT+1,NINTSB1)
      CALL RSAVE(ITAPT+2,NINTSB2)
      CALL RSAVE(ITAPT+3,NINTSB3)
      CALL RSAVE(ITAPT+4,NINTSB4)
      IMSTAT(1) = 3
      CALL CCDUMP
C
      IF (.NOT.SRTALL) GOTO 600
C
      NINTSB5 = MAX0(MREMAIN/512,20)
      NINTSB5 = MIN0(NV5*RCW/512+1,NINTSB5)
      CALL RFILE(ITAPT+5,NINTSB5)
      MREMAIN = MREMAIN - (NINTSB5-20)*512
C
      NINTSB6 = 20
      CALL RFILE(ITAPT+6,NINTSB6)
C
      IF (IMSTAT(2).GT.1) THEN
         WRITE (IW,1011)
         GOTO 550
      ENDIF
C
C     GET THE VOVV INTEGRALS
C
      NPASS = (NV5-1) / NVT + 1
      NSEC5 = 1
      N5L = 0
      WRITE (IW,1005) NPASS,'VOVV'
C
      DO IPASS = 1, NPASS
      N5F = N5L + 1
      N5L = MIN0(N5L+NVT,NV5)
      N5 = N5L - N5F + 1
      NSEC = N5*IRW*RCW / SEC2I(1)
      IF (MOD(N5*IRW*RCW,SEC2I(1)).NE.0) NSEC = NSEC + 1
C
      REWIND (MDINT)
      READ (MDINT)
C
      CALL XCOPY (N5,A0,0,VVOVV,1)
      CALL XCOPY (N5,A0,0,VVOVVA,1)
C
  200 CALL GET2INTS (MDINT,KR,NSP,INDEX,ITRSGN,INDK,INDL,
     &               NZ,CBUF,IBUF)
      IF (NZ.LT.0) GOTO 210
      DO 202 I = 1, NZ
         DO 201 IP = 1, 16
            ITYP = IBUF(IP,1,I)
            IREP = IBUF(IP,2,I) 
            IND1 = IBUF(IP,3,I)
            IND2 = IBUF(IP,4,I)
            ITRS = IBUF(IP,5,I)
            JTRS = IBUF(IP,6,I)
            IF (ITYP.EQ.14) THEN       
               IND = IVOVVT(IREP) + (IND1-1)*NVO(IREP) + IND2
               IF (IND.GE.N5F.AND.IND.LE.N5L) THEN
                  IND = IND - N5F + 1
                  VVOVV(1,IND) = CBUF(1,I) * ITRS
                  IF (CARITH) VVOVV(2,IND) = CBUF(2,I) * JTRS
               ENDIF
            ELSEIF (ITYP.EQ.-14) THEN
               IND = IVOVVT(IREP) + (IND1-1)*NVO(IREP) + IND2
               IF (IND.GE.N5F.AND.IND.LE.N5L) THEN
                  IND = IND - N5F + 1
                  VVOVVA(1,IND) = CBUF(1,I) * ITRS
                  IF (CARITH) VVOVVA(2,IND) = CBUF(2,I) * JTRS
               ENDIF
            ENDIF
  201    CONTINUE
  202 CONTINUE
      GOTO 200
C
C     ADD THE DIRECT AND EXCHANGE INTEGRALS
C
  210 CALL XAXPY (N5,-A1,VVOVVA,1,VVOVV,1)
C
C     ADD BREIT AND COULOMB INTEGRALS
C
      IF (BREIT.AND.ITRSGN.LT.0) THEN
         CALL RREAD(ITAPT+5,VVOVVA,N5*RCW*IRW,NSEC5)
         CALL XAXPY(N5,A1,VVOVVA,1,VVOVV,1)
      ENDIF
      CALL RWRIT(ITAPT+5,VVOVV,N5*RCW*IRW,NSEC5)
      NSEC5 = NSEC5 + NSEC
      ENDDO
C
C     GET THE VVVV INTEGRALS
C
      NPASS = (NV6-1) / NVT + 1
      NSEC6 = 1
      N6L = 0
      WRITE (IW,1005) NPASS,'VVVV'
C
      DO IPASS = 1, NPASS
      N6F = N6L + 1
      N6L = MIN0(N6L+NVT,NV6)
      N6 = N6L - N6F + 1
      NSEC = N6*IRW*RCW / SEC2I(1)
      IF (MOD(N6*IRW*RCW,SEC2I(1)).NE.0) NSEC = NSEC + 1
C
      REWIND (MDINT)
      READ (MDINT)
C
      CALL XCOPY (N6,A0,0,VVVVV,1)
      CALL XCOPY (N6,A0,0,VVVVVA,1)
C
  300 CALL GET2INTS (MDINT,KR,NSP,INDEX,ITRSGN,INDK,INDL,
     &               NZ,CBUF,IBUF)
      IF (NZ.LT.0) GOTO 310
      DO 302 I = 1, NZ
         DO 301 IP = 1, 16
            ITYP = IBUF(IP,1,I)
            IREP = IBUF(IP,2,I) 
            IND1 = IBUF(IP,3,I)
            IND2 = IBUF(IP,4,I)
            ITRS = IBUF(IP,5,I)
            JTRS = IBUF(IP,6,I)
            IF (ITYP.EQ.16) THEN       
               IND = IVVVVTT(IREP) + (IND1-1)*NVVT(IREP) + IND2
               IF (IND.GE.N6F.AND.IND.LE.N6L) THEN
                  IND = IND - N6F + 1
                  VVVVV(1,IND) = CBUF(1,I) * ITRS
                  IF (CARITH) VVVVV(2,IND) = CBUF(2,I) * JTRS
               ENDIF
            ELSEIF (ITYP.EQ.-16) THEN
               IND = IVVVVTT(IREP) + (IND1-1)*NVVT(IREP) + IND2
               IF (IND.GE.N6F.AND.IND.LE.N6L) THEN
                  IND = IND - N6F + 1
                  VVVVVA(1,IND) = CBUF(1,I) * ITRS
                  IF (CARITH) VVVVVA(2,IND) = CBUF(2,I) * JTRS
               ENDIF
            ENDIF
  301    CONTINUE
  302 CONTINUE
      GOTO 300
C
C     ADD THE DIRECT AND EXCHANGE INTEGRALS
C
  310 CALL XAXPY (N6,-A1,VVVVVA,1,VVVVV,1)
C
C     ADD BREIT AND COULOMB INTEGRALS
C
      IF (BREIT.AND.ITRSGN.LT.0) THEN
         CALL RREAD(ITAPT+6,VVVVVA,N6*RCW*IRW,NSEC6)
         CALL XAXPY(N6,A1,VVVVVA,1,VVVVV,1)
      ENDIF
      CALL RWRIT(ITAPT+6,VVVVV,N6*RCW*IRW,NSEC6)
      NSEC6 = NSEC6 + NSEC
      ENDDO
C
  500 CLOSE(MDINT)
C
C  NOW PASS THE BREIT INTEGRAL FILE IF DESIRED
C
      IF (BREIT.AND.ITRSGN.EQ.1) THEN
         WRITE (IW,1000) 'Coulomb',DATEX,TIMEX,ICOUNT
         OPEN(MDINT,FILE='MDBINT',FORM='UNFORMATTED')
         REWIND(MDINT)
         READ (MDINT,ERR=10001,END=10001) DATEX,TIMEX,NKR2
         IF (NKR2.NE.NKR) GOTO 10003
         ITRSGN = -1
         GOTO 1
      ELSEIF (BREIT) THEN
         WRITE (IW,1000) 'Breit',DATEX,TIMEX,ICOUNT
      ELSE
         WRITE (IW,1000) 'Coulomb',DATEX,TIMEX,ICOUNT
      ENDIF
C
C     This is a restart point : save the files and update status.
C
  550 CONTINUE
      CALL RSAVE(ITAPT+5,NINTSB5)
      CALL RSAVE(ITAPT+6,NINTSB6)
      IMSTAT(2) = 3
      CALL CCDUMP
C
  600 IF (SRTALL) THEN
      WRITE (IW,1004) 'OOOO',NV1,'VOOO',NV2,'VVOO',NV3,'VOVO',NV4,
     &                'VOVV',NV5,'VVVV',NV6
      ELSE
      WRITE (IW,1004) 'OOOO',NV1,'VOOO',NV2,'VVOO',NV3,'VOVO',NV4
      ENDIF
C
C  FORM FOCK MATRIX; ADD TWO-ELECTRON PART
C  CALCULATE THE SCF (OR REFERENCE) ENERGY
C
      CALL MAKEFM (VVOVO,VVOVOA,FOO,FVO,FVV,E0,E1)
C
C     Write the Fock matrix to file
C
      NSEC0 = 1
      CALL RWRIT(ITAPT+0,FVO,NFVO*RCW*IRW,NSEC0)
      NSEC = NFVO*IRW*RCW / SEC2I(1)
      IF (MOD(NFVO*IRW*RCW,SEC2I(1)).NE.0) NSEC = NSEC + 1
      NSEC0 = NSEC0 + NSEC
      CALL RWRIT(ITAPT+0,FOO,NFOO*RCW*IRW,NSEC0)
      NSEC = NFOO*IRW*RCW / SEC2I(1)
      IF (MOD(NFOO*IRW*RCW,SEC2I(1)).NE.0) NSEC = NSEC + 1
      NSEC0 = NSEC0 + NSEC
      CALL RWRIT(ITAPT+0,FVV,NFVV*RCW*IRW,NSEC0)
C
C     This is a restart point : save the files and update status.
C
      CALL RSAVE(ITAPT+0,NINTSB0)
      IMSTAT(3) = 3
      CALL CCDUMP
C
      WRITE (IW,1001)
      II = 1
      I  = 1
      DO IREP = 1, NREP
         DO J = 1, NO(IREP)
            EPSDIF = DABS(EPS(I)-FOO(II))
C64B        EPSDIF =  ABS(EPS(I)-FOO(II))
            IF (DABS(FOO(II)).GT.1.D0) EPSDIF = EPSDIF / DABS(FOO(II))
C64B        IF ( ABS(FOO(II)).GT.1.0 ) EPSDIF = EPSDIF /  ABS(FOO(II))
            IF (EPSDIF.GT.1.D-7.OR.IPRNT.GE.1)
     &         WRITE (IW,1002) I,REPNA(IREP),EPS(I),FOO(II)
            EPS(I) = FOO(II)
            I  = I  + 1
            II = II + (NO(IREP) + 1) * RCW
         ENDDO
         II = II - NO(IREP) * RCW
      ENDDO
      II = 1
      DO IREP = 1, NREP
         DO J = 1, NV(IREP)
            EPSDIF = DABS(EPS(I)-FVV(II))
C64B        EPSDIF =  ABS(EPS(I)-FVV(II))
            IF (DABS(FVV(II)).GT.1.D0) EPSDIF = EPSDIF / DABS(FVV(II))
C64B        IF ( ABS(FVV(II)).GT.1.0 ) EPSDIF = EPSDIF /  ABS(FVV(II))
            IF (EPSDIF.GT.1.D-7.OR.IPRNT.GE.2)
     &         WRITE (IW,1002) I,REPNA(IREP),EPS(I),FVV(II)
            EPS(I) = FVV(II)
            I  = I  + 1
            II = II + (NV(IREP) + 1) * RCW
         ENDDO
         II = II - NV(IREP) * RCW
      ENDDO
C
      ESCF = E0 + E1 
      WRITE(IW,1003) REPNUC,E0,E1,ESCF,ESCF+REPNUC
      ESCF = ESCF + REPNUC
      ETOT = ESCF
C
 1000 FORMAT (/' Reading ',A7,' integrals :'/' File date :',T18,A10
     &/' File time :',T20,A8,/' # of integrals',I8)
 1005 FORMAT (/' Need ',I6,' passes to sort ',A4,' integrals')
 1004 FORMAT (/' Expanded and sorted integrals to 6 types :'/
     & (' Type ',A4,' :',I10,' integrals'))
 1001 FORMAT (//' Checking the orbital energies, the program uses ',
     & 'the diagonal elements of the'/' reconstructed Fock matrix.',
     & ' Differences with the reference orbital energies'/
     & ' are given if above a treshold or if iprnt > 1'
     &//' Spinor  Abelian Rep.',9X,'Energy   Recalc. Energy')
 1002 FORMAT (2x,I5,6X,A4,4X,2F16.10)
 1003 FORMAT (//' Nuclear repulsion + core energy :',T40,F25.15
     &/' Zero order electronic energy :',T40,F25.15 
     &/' First order electronic energy :',T40,F25.15 
     &/' Electronic energy :',T40,F25.15 
     &/' SCF energy :',T40,F25.15) 
 1010 FORMAT (/' Skipping integral sorting step 1')
 1011 FORMAT (/' Skipping integral sorting step 2')
 1012 FORMAT (/' Skipping Fock-matrix build')
      RETURN
10000 CALL QUIT('ERROR READING HEADER MDCINT')
10001 CALL QUIT('ERROR READING HEADER MDBINT')
10003 CALL QUIT('MDCINT AND MDBINT FILES ARE INCOMPATIBLE')
10010 CALL QUIT('ERROR READING INTEGRALS FROM MDCINT')
10011 CALL QUIT('ERROR READING INTEGRALS FROM MDBINT')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RDPROP(IPRNT,NAMEA,NSP,RDDIAG,PROP,IREPSPI,APHASE,
     &                  AVO,AOV,AOO,AVV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Reads property integral files for operator A
C     If (RDDIAG) it reads also the diagonal blocks (OO and VV)
C     This version reads from the file created by PRTRAN
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NSP,IREPSPI(NSP,8,2)
      REAL*8 AVO(*),AOV(*),AOO(*),AVV(*)
      REAL*8 PROP(2,NSP,NSP)
      REAL*8 APHASE(2)
      CHARACTER*8 NAMEA
      LOGICAL RDDIAG,PHASE
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
      INCLUDE "files.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER*32 ACHAR
      REAL*8 AVOI,AOVI,AMAXR,AMAXI,ANORM
C
C---------------Executable code--------------------------------------
C
      OPEN (MDPROP,FILE='MDPROP',FORM='UNFORMATTED')
    1 READ (MDPROP,END=100,ERR=100) ACHAR
      IF (ACHAR(1:8).NE.'********'.OR.ACHAR(25:32).NE.NAMEA) GOTO 1
      IF (IPRNT.GE.1) WRITE (IW,1000) NAMEA,ACHAR(9:16),ACHAR(17:24)
      READ (MDPROP) PROP
      CLOSE (MDPROP,STATUS='KEEP')
C
      APHASE(1) = D1
      APHASE(2) = D0
      PHASE = .FALSE.
C
    2 AMAXI = D0
      AAMAX = 0
      IIMAX = 0
      AI = 1
      DO AIREP = 1, NREP
         DO IREP = 1, NREP
            AREP = MULTB (IREP,AIREP+NREP,2)
            AREP = MULTB (AREP,NREP+1,2)
            DO I = 1, NO(IREP)
               II = IREPSPI(I,IREP,1)
               DO A = 1, NV(AREP)
                  AA = IREPSPI(A,AREP,2)
                  IF (CARITH) THEN
                      AVO(RCW*AI-1) = PROP(1,AA,II)
                      AVO(RCW*AI)   = PROP(2,AA,II)
                      AOV(RCW*AI-1) = PROP(1,II,AA)
                      AOV(RCW*AI)   = PROP(2,II,AA)
                   ELSE
                      AVO(AI) = PROP(1,AA,II)*APHASE(1) -
     &                          PROP(2,AA,II)*APHASE(2)
                      AOV(AI) = PROP(1,II,AA)*APHASE(1) -
     &                          PROP(2,II,AA)*APHASE(2)
                      AVOI = PROP(1,AA,II)*APHASE(2) +
     &                       PROP(2,AA,II)*APHASE(1)
                      AOVI = PROP(1,II,AA)*APHASE(2) +
     &                       PROP(2,II,AA)*APHASE(1)
                      IF (ABS(AVOI).GT.ABS(AMAXI)) THEN
                         AMAXR = AVO(AI) 
                         AMAXI = AVOI
                         AAMAX = AA
                         IIMAX = II
                      ENDIF
                      IF (ABS(AOVI).GT.ABS(AMAXI)) THEN
                         AMAXR = AOV(AI) 
                         AMAXI = AOVI 
                         AAMAX = II
                         IIMAX = AA
                      ENDIF
                   ENDIF
                   AI = AI + 1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      IF (.NOT.CARITH.AND.ABS(AMAXI).GT.ACCUR) THEN
         IF (PHASE) THEN
            WRITE (IW,1010) AA,II,AMAXI
            CALL QUIT('USE COMPLEX ARITHMETICS')
         ELSE
            APHASE(1) = AMAXR
            APHASE(2) = AMAXI 
            ANORM = DSQRT(AMAXR*AMAXR+AMAXI*AMAXI)
C64B        ANORM =  SQRT(AMAXR*AMAXR+AMAXI*AMAXI)
            APHASE(1) = APHASE(1) / ANORM
            APHASE(2) = APHASE(2) / ANORM
            PHASE = .TRUE.
            IF (IPRNT.GE.1) WRITE (6,1020) NAMEA,APHASE
            GOTO 2
         ENDIF
      ENDIF
C
      IF (.NOT.RDDIAG) RETURN
C
      IJ = 1
      DO JREP = 1, NREP
         IREP = JREP
         DO J = 1, NO(JREP)
            JJ = IREPSPI(J,JREP,1)
            DO I = 1, NO(IREP)
               II = IREPSPI(I,IREP,1)
               IF (CARITH) THEN
                  AOO(RCW*IJ-1) = PROP(1,II,JJ)
                  AOO(RCW*IJ)   = PROP(2,II,JJ)
               ELSE
                  AOO(IJ) = PROP(1,II,JJ)*APHASE(1) -
     &                      PROP(2,II,JJ)*APHASE(2)
               ENDIF
               IJ = IJ + 1
            ENDDO
         ENDDO
      ENDDO
C
      AB = 1
      DO BREP = 1, NREP
         AREP = BREP
         DO B = 1, NV(BREP)
            BB = IREPSPI(B,BREP,2)
            DO A = 1, NV(AREP)
               AA = IREPSPI(A,AREP,2)
               IF (CARITH) THEN
                  AVV(RCW*AB-1) = PROP(1,AA,BB)
                  AVV(RCW*AB)   = PROP(2,AA,BB)
               ELSE
                  AVV(AB) = PROP(1,AA,BB)*APHASE(1) -
     &                      PROP(2,AA,BB)*APHASE(2)
               ENDIF
               AB = AB + 1
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
  100 WRITE (6,*) ' Property ',NAMEA,' not found on file MDPROP'
      CALL QUIT(' Property integrals missing')
 1000 FORMAT (/' Read integral type ',A8,' created ',
     & A8,' storage info : ',A8)
 1010 FORMAT (/' Largest imaginary part of matrix element',2I5,F10.2)
 1020 FORMAT (/' Property ',A8,' scaled with phase factor ',2F10.4)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET2INTS (MDINT,KR,NSPO,INDEX,ITRSGN,INDK,INDL,
     &                     NZ,CBUF,IBUF)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Reads block of 2-electron integrals from MOLFDIR file MDINT
C     Expands Kramers' reduced list into full list.
C     Reorder from (ij|kl) to <ik||jl>
C     Expanded (16 permutations) indices are written in array IBUF :
C     IBUF(IPERM,1,INZ) : Type (OOOO, VOOO, VVOO, VOVO, VVVO, VVVV)
C     IBUF(IPERM,2,INZ) : Number of combined irrep
C     IBUF(IPERM,3,INZ) : Index of ket pair
C     IBUF(IPERM,4,INZ) : Index of bra pair
C     IBUF(IPERM,5,INZ) : Index of exchange ket pair
C     IBUF(IPERM,6,INZ) : Sign under time-reversion (real part)
C     IBUF(IPERM,7,INZ) : Sign under time-reversion and conjugations 
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER KR(-NSPO:NSPO),INDK(*),INDL(*)
      REAL*8 CBUF(RCW,*)
      INTEGER IBUF(16,6,*)
      INTEGER INDEX(NSPO,NSPO,4)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER IKR,JKR,KKR,LKR
      REAL*8 RDUM
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN 
         READ (MDINT,END=10010,ERR=10010) IKR,JKR,NZ,
     &        (INDK(INZ),INDL(INZ),INZ=1,NZ),
     &        (CBUF(1,INZ),CBUF(2,INZ),INZ=1,NZ) 
      ELSE
         READ (MDINT,END=10010,ERR=10010) IKR,JKR,NZ,
     &        (INDK(INZ),INDL(INZ),INZ=1,NZ),
     &        (CBUF(1,INZ),RDUM,INZ=1,NZ) 
      ENDIF
      IF (IKR.EQ.0) THEN
         NZ = -1
         RETURN
      ENDIF
      IJSIGN = SIGN(1,IKR*JKR) * ITRSGN
C
C     Expand to the full list
C
      DO 50 INZ = 1, NZ
         KKR = INDK(INZ)
         LKR = INDL(INZ)
         I0 = KR(IKR)
         J0 = KR(JKR)
         IPERM = 0 
         DO 40 ITR1 = 1, 2
            K0 = KR(KKR)
            L0 = KR(LKR)
            DO 30 ITR2 = 1, 2
               I1 = I0
               J1 = J0
               K1 = K0
               L1 = L0
               DO 20 IPI = 1, 2
                  I = I1
                  J = J1
                  K = K1
                  L = L1
                  DO 10 ICC = 1, 2
                     IPERM = IPERM + 1
                     JLIND = INDEX(J,L,4)
                     IKIND = INDEX(I,K,4)
                     IF (JLIND.EQ.0.OR.I.EQ.K) THEN ! Non-canonical integral
                        IBUF(IPERM,1,INZ) = 0
                     ELSEIF (IKIND.NE.0) THEN  ! Contribution to direct part
                        ITYP = (INDEX(J,L,1)-1)*4 + INDEX(I,K,1)
                        IBUF(IPERM,1,INZ) = ITYP
                        IBUF(IPERM,2,INZ) = INDEX(J,L,3)
                        IBUF(IPERM,3,INZ) = INDEX(J,L,4)
                        IBUF(IPERM,4,INZ) = INDEX(I,K,4)
                     ELSE       ! Contribution to exchange part
                        ITYP = (INDEX(J,L,1)-1)*4 + INDEX(K,I,1)
                        ITYP = - ITYP
                        IBUF(IPERM,1,INZ) = ITYP
                        IBUF(IPERM,2,INZ) = INDEX(J,L,3)
                        IBUF(IPERM,3,INZ) = INDEX(J,L,4)
                        IBUF(IPERM,4,INZ) = INDEX(K,I,4)
                     ENDIF
                     J = I1
                     I = J1
                     L = K1
                     K = L1
   10             CONTINUE
                  I1 = K0
                  J1 = L0
                  K1 = I0
                  L1 = J0
   20          CONTINUE
               K0 = KR(-LKR)
               L0 = KR(-KKR)
   30       CONTINUE
            I0 = KR(-JKR)
            J0 = KR(-IKR)
   40    CONTINUE
   50 CONTINUE
C
      DO INZ = 1, NZ
         DO IPERM = 1, 8
            IBUF(IPERM,5,INZ) = 1
         ENDDO
         DO IPERM = 9, 16
            IBUF(IPERM,5,INZ) = IJSIGN
         ENDDO
      ENDDO
C
      DO INZ = 1, NZ
         KKR = INDK(INZ)
         LKR = INDL(INZ)
         KLSIGN = SIGN(1,KKR*LKR) * ITRSGN
         DO IPERM = 5, 8
            IBUF(IPERM,5,INZ) = IBUF(IPERM,5,INZ) * KLSIGN
         ENDDO
         DO IPERM = 13, 16
            IBUF(IPERM,5,INZ) = IBUF(IPERM,5,INZ) * KLSIGN
         ENDDO
      ENDDO
C
      DO INZ = 1, NZ
         DO IPERM = 1, 15, 2
            IBUF(IPERM,6,INZ) = IBUF(IPERM,5,INZ)
            IBUF(IPERM+1,6,INZ) = - IBUF(IPERM+1,5,INZ)
         ENDDO
      ENDDO
C
      RETURN
10010 CALL QUIT('ERROR READING 2-E INTEGRALS')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKEFM (BUF1,BUF2,FOO,FVO,FVV,E0,E1)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Generates Fock matrix. Calculates E0 and E1
C
C---------------Routines called----------------------------------------
C
C     GETOOOO, GETVOOO, GETVOVO, SRT1TT4, SRT1ST4, SRT1SS4
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 E0,E1
      REAL*8 FOO(*),FVO(*),FVV(*)
      REAL*8 BUF1(*),BUF2(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
      INCLUDE "symm.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 XSUM
C
C---------------Executable code--------------------------------------
C
      E0 = AR0
      E1 = AR0
C
      CALL GETOOOO (BUF1)
      CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NO,NO,NO,NO,MOO,JOOOO,JJOO,JJOO,
     &              BUF1,BUF2)
      KK = 1
      DO KREP = 1, NREP
         DO K = 1, NO(KREP)
C---------------------------------------------------------------------
C  F(IJ) = F(IJ) + V(IJ,KK)
C---------------------------------------------------------------------
            CALL XAXPY(NFOO,A1,BUF2(KK),1,FOO,1)
C---------------------------------------------------------------------
C  E1 = E1 - V(II,KK)
C---------------------------------------------------------------------
            II = 0
            DO IREP = 1, NREP
               E1 = E1 - XSUM(NO(IREP),BUF2(II+KK),NO(IREP)+1)
               II = II + NO(IREP) * NO(IREP) * RCW
            ENDDO
            KK = KK + (NO(KREP)+1) * NFOO * RCW
         ENDDO
         KK = KK - NO(KREP) * NFOO * RCW
      ENDDO
      E1 = E1 * DP5
C---------------------------------------------------------------------
C  E0 = E0 + FOO(I,I)
C---------------------------------------------------------------------
      II = 1
      DO IREP = 1, NREP
         E0 = E0 + XSUM(NO(IREP),FOO(II),NO(IREP)+1)
         II = II + NO(IREP) * NO(IREP) * RCW
      ENDDO
C
      CALL GETVOOO (BUF1)
      CALL SRT1ST4 (NREP,MULTB,LTR,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &              BUF1,BUF2)
      KK = 1
      DO KREP = 1, NREP
         DO K = 1, NO(KREP)
C---------------------------------------------------------------------
C  F(AI) = F(AI) + V(AI,KK)
C---------------------------------------------------------------------
            CALL XAXPY(NFVO,A1,BUF2(KK),1,FVO,1)
            KK = KK + (NO(KREP)+1) * NFVO * RCW
         ENDDO
         KK = KK - NO(KREP) * NFVO * RCW
      ENDDO
C
      CALL GETVOVO (BUF1)
      CALL SRT1SS4 (NREP,MULTB,NV,NO,NV,NO,MVV,JVVOO,JJVV,JJOO,
     &              BUF1,BUF2)
      KK = 1
      DO KREP = 1, NREP
         DO K = 1, NO(KREP)
C---------------------------------------------------------------------
C  F(AB) = F(AB) + V(AB,KK)
C---------------------------------------------------------------------
            CALL XAXPY(NFVV,A1,BUF2(KK),1,FVV,1)
            KK = KK + (NO(KREP)+1) * NFVV * RCW
         ENDDO
         KK = KK - NO(KREP) * NFVV * RCW
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MKINDEX (NSP,IREPSP,IREPSPI,INDEX)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Set up index array for initial integral sort
C     Index(i,j,1) : type (oo,vo,ov,vv)
C     Index(i,j,2) : position of fm(i,j) in appropriate part of fock matrix
C     Index(i,j,3) : symmetry of ij pair
C     Index(i,j,4) : position of ij pair in batch of this symmetry
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
      INTEGER IREPSP(NSP,3),IREPSPI(NSP,8,2),INDEX(NSP,NSP,4)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "symm.inc"
C
C---------------Local variables--------------------------------------
C
C
C---------------Executable code--------------------------------------
C
      DO J = 1, NSP
         IF (IREPSP(J,3).GT.0) THEN
            JJ = 0 
         ELSE 
            JJ = 2
         ENDIF
         DO I = 1, NSP
            IF (IREPSP(I,3).GT.0) THEN
               II = 1 
            ELSE 
               II = 2
            ENDIF
            INDEX(I,J,1) = II + JJ
            INDEX(I,J,2) = 0
            INDEX(I,J,3) = 0
            INDEX(I,J,4) = 0
         ENDDO
      ENDDO
C
      IJ = 0
      DO IREP = 1, NREP
         DO JJ = 1, NO(IREP)
            J = IREPSPI(JJ,IREP,1)
            DO II = 1, NO(IREP)
               I = IREPSPI(II,IREP,1)
               IJ = IJ + 1
               INDEX(I,J,2) = IJ
            ENDDO
         ENDDO
      ENDDO
C
      IJ = 0
      DO IREP = 1, NREP
         DO JJ = 1, NO(IREP)
            J = IREPSPI(JJ,IREP,1)
            DO II = 1, NV(IREP)
               I = IREPSPI(II,IREP,2)
               IJ = IJ + 1
               INDEX(I,J,2) = IJ
            ENDDO
         ENDDO
      ENDDO
C
      IJ = 0
      DO IREP = 1, NREP
         DO JJ = 1, NV(IREP)
            J = IREPSPI(JJ,IREP,2)
            DO II = 1, NO(IREP)
               I = IREPSPI(II,IREP,1)
               IJ = IJ + 1
               INDEX(I,J,2) = IJ
            ENDDO
         ENDDO
      ENDDO
C
      IJ = 0
      DO IREP = 1, NREP
         DO JJ = 1, NV(IREP)
            J = IREPSPI(JJ,IREP,2)
            DO II = 1, NV(IREP)
               I = IREPSPI(II,IREP,2)
               IJ = IJ + 1
               INDEX(I,J,2) = IJ
            ENDDO
         ENDDO
      ENDDO
C
      DO IJREP = 1, NREP
         IJ = 0
         DO JREP = 1, NREP
            DO 10 IREP = JREP, NREP
               IF (MULTB(IREP,JREP,1).NE.IJREP) GOTO 10
               DO JJ = 1, NO(JREP)
                  J = IREPSPI(JJ,JREP,1)
                  IMIN = 1
                  IF (IREP.EQ.JREP) IMIN = JJ + 1
                  DO II = IMIN, NO(IREP)
                     I = IREPSPI(II,IREP,1)
                     IJ = IJ + 1
                     INDEX(I,J,3) = IJREP
                     INDEX(I,J,4) = IJ
                  ENDDO
               ENDDO
 10         CONTINUE
         ENDDO
      ENDDO
C
      DO IJREP = 1, NREP
         IJ = 0
         DO JREP = 1, NREP
            DO 20 IREP = 1, NREP
               IF (MULTB(IREP,JREP,1).NE.IJREP) GOTO 20
               DO JJ = 1, NO(JREP)
                  J = IREPSPI(JJ,JREP,1)
                  DO II = 1, NV(IREP)
                     I = IREPSPI(II,IREP,2)
                     IJ = IJ + 1
                     INDEX(I,J,3) = IJREP
                     INDEX(I,J,4) = IJ
                  ENDDO
               ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
C
      DO IJREP = 1, NREP
         IJ = 0
         DO JREP = 1, NREP
            DO 40 IREP = JREP, NREP
               IF (MULTB(IREP,JREP,1).NE.IJREP) GOTO 40
               DO JJ = 1, NV(JREP)
                  J = IREPSPI(JJ,JREP,2)
                  IMIN = 1
                  IF (IREP.EQ.JREP) IMIN = JJ + 1
                  DO II = IMIN, NV(IREP)
                     I = IREPSPI(II,IREP,2)
                     IJ = IJ + 1
                     INDEX(I,J,3) = IJREP
                     INDEX(I,J,4) = IJ
                  ENDDO
               ENDDO
 40         CONTINUE
         ENDDO
      ENDDO
C
      RETURN 
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
