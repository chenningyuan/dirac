*  Read DIRAC-COULOMB two-electron integrals
*  (Later also BREIT integrals)
*  (GAUNT integrals can be included in COULOMB type ints.,
*   because they have the same symmetry.)
*
      subroutine rdtwo_dirac(RKLR,RKLI,LUTWO,DIRACTWO,
     &                       NDIM,ISPINFREE,MORB,
     &                       KR,IOFF12,NUMINT,ITPLOC,INDK,INDL)
*
      implicit real*8 (A-H,O-Z)
*
#include "input.inc"
*
      character*7 DIRACTWO
      character*10 DATEX,DATEXB,TIMEX*8,TIMEXB*8
      logical REALAR
c     logical BREIT
*
      dimension RKLR((2*NDIM)**4),RKLI((2*NDIM)**4)
      dimension KR(-NDIM:NDIM),INDK((2*NDIM)**4),INDL((2*NDIM)**4),
     &          IOFF12(MORB,MORB,4),NUMINT(MORB,MORB,4),
     &          ITPLOC(MORB,MORB)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
      REALAR = .TRUE.
      TRESINT = 1.D-12
      NINT = 0
      KLS = 0
      KLF = M2
      MM = M1 * M2
*
      open (LUTWO,FILE=DIRACTWO,FORM='UNFORMATTED')
      rewind (LUTWO)
      read (LUTWO) datex,timex,nkr,
     & (kr(i),kr(-i),i=1,nkr)
*
      if (2*NKR.ne.MORB) then
        write(6,*) 'MRCONEE and MDCINT are incompatible.'
        stop 'Quitting.'
      end if
      write(6,*) 'Integrals have been calculated on the'
      write(6,'(A12,X,A)') '            ',DATEX,' at ',TIMEX
      if (NTEST.ge.5) then
        write(6,*)
        write(6,*) 'Info from MDCINT file'
        write(6,*) 'NKR is ',NKR
        write(6,*) 'KR array:'
        do I=1,NKR,1
          write(6,'(2I8)') KR(-I),KR(I)
        end do
      end if
C
C     -------------------------
C     READ IN COULOMB INTEGRALS
C     -------------------------
C
* Initialize
      ISYMOFF = 0
      ITYPE = 1
      IKR_SAVE = 0
      JKR_SAVE = 0
      do M=1,MORB,1
        do N=1,MORB,1
          ITPLOC(M,N) = 0
          do IT=1,4,1
            NUMINT(M,N,IT) = 0
          end do
        end do
      end do
*
      rewind (LUTWO)
      read (LUTWO)
*
* Read a block of integrals
  10  continue
      read (LUTWO) IKR,JKR,NZ,
     &             (INDK(INZ+ISYMOFF),INDL(INZ+ISYMOFF),INZ=1,NZ,1),
     &             (RKLR(INZ+ISYMOFF),RKLI(INZ+ISYMOFF),INZ=1,NZ,1)
*
* End of file 
* (FIXME: No multiple files implemented, yet!)
      if (IKR.eq.0) then
        write(6,*) 'End of file MDCINT.'
        goto 20
      end if
*
      ITPLOC(KR(IKR),KR(JKR)) = ITPLOC(KR(IKR),KR(JKR)) + 1
      IOFF12(KR(IKR),KR(JKR),ITPLOC(KR(IKR),KR(JKR))) = ISYMOFF
      NUMINT(KR(IKR),KR(JKR),ITPLOC(KR(IKR),KR(JKR))) = NZ
*
      if (NTEST.ge.10) then
        write(6,'(A,2I6)') 'IKR,JKR : ',IKR,JKR
        write(6,'(A,3I4,A,I6)') 
     &        'IOFF12(',KR(IKR),KR(JKR),ITPLOC(KR(IKR),KR(JKR)),') = ',
     &         IOFF12(KR(IKR),KR(JKR),ITPLOC(KR(IKR),KR(JKR)))
        write(6,'(A,3I4,A,I6)')
     &        'NUMINT(',KR(IKR),KR(JKR),ITPLOC(KR(IKR),KR(JKR)),') = ',
     &         NUMINT(KR(IKR),KR(JKR),ITPLOC(KR(IKR),KR(JKR)))
        write(6,'(A,X,2I6)') 'Integral block: ',KR(IKR),KR(JKR)
        do INT=1,NZ,1
          write(6,'(A,2I8,X,2F18.12)') 
     &      'Indices K,L ',KR(INDK(INT+ISYMOFF)),KR(INDL(INT+ISYMOFF)),
     &                     RKLR(INT+ISYMOFF),RKLI(INT+ISYMOFF)
        end do
      end if
*
      ISYMOFF = ISYMOFF + NZ
      IKR_SAVE = IKR
      JKR_SAVE = JKR
      goto 10
*
   20 if (NTEST.ge.5) then
        write(6,*)
        write(6,*) '===================='
        write(6,*) 'rdtwo_dirac speaking'
        write(6,*) '===================='
        write(6,*)
        write(6,*) ' Index pair, int type, offset, number of integrals:'
        do IKPAIR=1,-1,-2
          do II=IKPAIR,MORB*IKPAIR/2,IKPAIR
            do JKPAIR=1,-1,-2
              do IJ=JKPAIR,MORB*JKPAIR/2,JKPAIR
                do IT=1,ITPLOC(KR(II),KR(IJ)),1
                  if (NUMINT(KR(II),KR(IJ),IT).gt.0) then
                    write(6,'(1X,5I5)') KR(II),KR(IJ),IT,
     &                                  IOFF12(KR(II),KR(IJ),IT),
     &                                  NUMINT(KR(II),KR(IJ),IT)
                    do INT=1,NUMINT(KR(II),KR(IJ),IT),1
                      write(6,'(5X,2I4,2F18.10)') 
     &                 KR(INDK(INT+IOFF12(KR(II),KR(IJ),IT))),
     &                 KR(INDL(INT+IOFF12(KR(II),KR(IJ),IT))),
     &                 RKLR(INT+IOFF12(KR(II),KR(IJ),IT)),
     &                 RKLI(INT+IOFF12(KR(II),KR(IJ),IT))
                    end do
                  end if
                end do
              end do
            end do
          end do
        end do
      end if


c     DO INZ = 1, NZ
c        IF (ABS(RKLI(INZ)).GT.TRESINT) REALAR = .FALSE.
c     ENDDO
c     MINT = MINT + NZ
c     If (ikr.EQ.0) GoTo 40 ! End of File
c     i = indmo(kr(ikr))
c     j = indmo(kr(jkr))
c     itr = indmo(kr(-ikr))
c     jtr = indmo(kr(-jkr))
c     SignIJ = SIGN(1,ikr) * SIGN(1,jkr)
c     DO INZ = 1, NZ
c        kkr = indk(inz)
c        lkr = indl(inz)
c        k = indmo(kr(kkr))
c        l = indmo(kr(lkr))
c        ktr = indmo(kr(-kkr))
c        ltr = indmo(kr(-lkr))
c        SignKL = SIGN(1,kkr) * SIGN(1,lkr)
C--> Original integral
c        IF (I.NE.0.AND.J.NE.0.AND.K.NE.0.AND.L.NE.0) THEN
c        TWOR = RKLR(INZ)
c        TWOI = RKLI(INZ)
c        IJ = (I-1) * NORB + J
c        KL = (K-1) * NORB + L
c        JI = (J-1) * NORB + I
c        LK = (L-1) * NORB + K
c        IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
c           TIR(IJ+(KL-KLS-1)*M1) = TWOR
c           TII(IJ+(KL-KLS-1)*M1) = TWOI
c        ENDIF
c        IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
c           TIR(JI+(LK-KLS-1)*M1) = TWOR
c           TII(JI+(LK-KLS-1)*M1) = - TWOI
c        ENDIF
c        IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
c           TIR(KL+(IJ-KLS-1)*M1) = TWOR
c           TII(KL+(IJ-KLS-1)*M1) = TWOI
c        ENDIF
c        IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
c           TIR(LK+(JI-KLS-1)*M1) = TWOR
c           TII(LK+(JI-KLS-1)*M1) = - TWOI
c        ENDIF
c        ENDIF ! check on deleted spinors
c     ENDDO
C
C
c  40 IF (MINT.EQ.0) WRITE (6,1010) 'MDCINT'
C
c     IF (.NOT.BREIT) GOTO 61
C
C     -----------------------
C     READ IN BREIT INTEGRALS
C     -----------------------
C
c  50 read (mdbint,ERR=10011,END=10011) ikr,jkr,nz,
c    &             (indk(inz),indl(inz),inz=1,nz),
c    &             (rklr(inz),rkli(inz),inz=1,nz)
c     DO INZ = 1, NZ
c        IF (ABS(RKLI(INZ)).GT.TRESINT) REALAR = .FALSE.
c     ENDDO
c     MINTB = MINTB + NZ
c     If (ikr.EQ.0) GoTo 60 ! End of File
c     i = indmo(kr(ikr))
c     j = indmo(kr(jkr))
c     itr = indmo(kr(-ikr))
c     jtr = indmo(kr(-jkr))
c     SignIJ = - SIGN(1,ikr) * SIGN(1,jkr)
c     DO INZ = 1, NZ
c        kkr = indk(inz)
c        lkr = indl(inz)
c        k = indmo(kr(kkr))
c        l = indmo(kr(lkr))
c        ktr = indmo(kr(-kkr))
c        ltr = indmo(kr(-lkr))
c        SignKL = - SIGN(1,kkr) * SIGN(1,lkr)
c--> Original integral
c        IF (I.NE.0.AND.J.NE.0.AND.K.NE.0.AND.L.NE.0) THEN
c        TWOR = RKLR(INZ)
c        TWOI = RKLI(INZ)
c        IJ = (I-1) * NORB + J
c        KL = (K-1) * NORB + L
c        JI = (J-1) * NORB + I
c        LK = (L-1) * NORB + K
c        IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
c           TIR(IJ+(KL-KLS-1)*M1+MM) = TWOR
c           TII(IJ+(KL-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
c           TIR(JI+(LK-KLS-1)*M1+MM) = TWOR
c           TII(JI+(LK-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
c           TIR(KL+(IJ-KLS-1)*M1+MM) = TWOR
c           TII(KL+(IJ-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
c           TIR(LK+(JI-KLS-1)*M1+MM) = TWOR
c           TII(LK+(JI-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        ENDIF ! check on deleted spinors
c--> Time reversal left
c        IF (ITR.NE.0.AND.JTR.NE.0.AND.K.NE.0.AND.L.NE.0) THEN
c        TWOR = RKLR(INZ) * SignIJ
c        TWOI = RKLI(INZ) * SignIJ
c        IJ = (JTR-1) * NORB + ITR
c        KL = (K-1) * NORB + L
c        JI = (ITR-1) * NORB + JTR
c        LK = (L-1) * NORB + K
c        IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
c           TIR(IJ+(KL-KLS-1)*M1+MM) = TWOR
c           TII(IJ+(KL-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
c           TIR(JI+(LK-KLS-1)*M1+MM) = TWOR
c           TII(JI+(LK-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
c           TIR(KL+(IJ-KLS-1)*M1+MM) = TWOR
c           TII(KL+(IJ-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
c           TIR(LK+(JI-KLS-1)*M1+MM) = TWOR
c           TII(LK+(JI-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        ENDIF ! check on deleted spinors
c--> Time reversal right
c        IF (I.NE.0.AND.J.NE.0.AND.KTR.NE.0.AND.LTR.NE.0) THEN
c        TWOR = RKLR(INZ) * SignKL
c        TWOI = RKLI(INZ) * SignKL
c        IJ = (I-1) * NORB + J
c        KL = (LTR-1) * NORB + KTR
c        JI = (J-1) * NORB + I
c        LK = (KTR-1) * NORB + LTR
c        IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
c           TIR(IJ+(KL-KLS-1)*M1+MM) = TWOR
c           TII(IJ+(KL-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
c           TIR(JI+(LK-KLS-1)*M1+MM) = TWOR
c           TII(JI+(LK-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
c           TIR(KL+(IJ-KLS-1)*M1+MM) = TWOR
c           TII(KL+(IJ-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
c           TIR(LK+(JI-KLS-1)*M1+MM) = TWOR
c           TII(LK+(JI-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        ENDIF ! check on deleted spinors
c--> Time reversal both
c        IF (ITR.NE.0.AND.JTR.NE.0.AND.KTR.NE.0.AND.LTR.NE.0) THEN
c        TWOR = RKLR(INZ) * SignIJ * SignKL
c        TWOI = RKLI(INZ) * SignIJ * SignKL
c        IJ = (JTR-1) * NORB + ITR
c        KL = (LTR-1) * NORB + KTR
c        JI = (ITR-1) * NORB + JTR
c        LK = (KTR-1) * NORB + LTR
c        IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
c           TIR(IJ+(KL-KLS-1)*M1+MM) = TWOR
c           TII(IJ+(KL-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
c           TIR(JI+(LK-KLS-1)*M1+MM) = TWOR
c           TII(JI+(LK-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
c           TIR(KL+(IJ-KLS-1)*M1+MM) = TWOR
c           TII(KL+(IJ-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
c           TIR(LK+(JI-KLS-1)*M1+MM) = TWOR
c           TII(LK+(JI-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        ENDIF ! check on deleted spinors
c     ENDDO
c
c     GOTO 50
c
c  60 IF (MINTB.EQ.0) WRITE (6,1010) 'MDBINT'
c
c     -----------------------------------------
c     ADD BREIT INTEGRALS AND COULOMB INTEGRALS
c     -----------------------------------------
c
c     DO IJKL = 1, (KLF-KLS) * M1
c        TIR(IJKL) = TIR(IJKL) + TIR(IJKL+MM)
c        TII(IJKL) = TII(IJKL) + TII(IJKL+MM)
c     ENDDO
c
   61 CONTINUE
C
C     ----------------------------------
C     WRITE 2-ELECTRON INTEGRALS TO FILE
C     ----------------------------------
C
c     IJKL = 0
c     DO 200 KL = KLS + 1, KLF
c     TIR (KL+(KL-KLS-1)*M1) = 0.5D0 * TIR(KL+(KL-KLS-1)*M1)
c     TII (KL+(KL-KLS-1)*M1) = 0.5D0 * TII(KL+(KL-KLS-1)*M1)
c     NKL = 0
c     DO 100 IJ = 1, M1
c        IJKL = IJKL + 1
c        IF (ABS(TIR(IJKL))+ABS(TII(IJKL)).GT.TRESINT) THEN
c           NKL = NKL + 1
c           INDK(NKL) = IJ
c           RKLR(NKL) = TIR(IJKL)
c           RKLI(NKL) = TII(IJKL)
c         ENDIF
c 100     CONTINUE
c     NINT = NINT + NKL
c     IF (REALAR) THEN
c        WRITE (MRCTWOE) KL,NKL,(INDK(IX),IX=1,NKL),(RKLR(IX),IX=1,NKL)
c     ELSE
c        WRITE (MRCTWOE) KL,NKL,(INDK(IX),IX=1,NKL),(RKLR(IX),IX=1,NKL),
c    &                   (RKLI(IX),IX=1,NKL)
c     ENDIF
c 200 CONTINUE
c     KLS = KLF
c     KLF = MIN0 (M1, KLF+M2)
c     IF (KLS.NE.M1) GOTO 10
c     IF (.NOT.BREIT) THEN
c        WRITE (6,1000) MINT,DATEX,TIMEX,NINT
c     ELSE
c        WRITE(6,1001) MINT,DATEX,TIMEX,MINTB,DATEXB,TIMEXB,NINT
c     ENDIF
c     IF (NINT.GT.0) THEN
c        WRITE (MRCTWOE) -NINT,0
c     ELSE
c        WRITE (MRCTWOE) -1,0
c     ENDIF
c     IF (REALAR) THEN
c        WRITE (6,1002)
c     ELSE
c        WRITE (6,1003)
c     ENDIF
C
c     CLOSE (MDCINT,STATUS='KEEP')
c     IF (BREIT) CLOSE (MDBINT,STATUS='KEEP')
      RETURN
 1000 FORMAT (
     &/I8,' Unique Coulomb integrals read from MDCINT (',A10,1X,A8,')'
     &/I8,' Coulomb integrals written to MRCTWOE')
 1001 FORMAT (
     &/I8,' Unique Coulomb integrals read from MDCINT (',A10,1X,A8,')'
     &/I8,' Unique Breit integrals read from MDBINT (',A10,1X,A8,')'
     &/I8,' Combined integrals written to MRCTWOE')
 1002 FORMAT (/' Real arithmetic will be used in EXPAND and TWOLIN')
 1003 FORMAT (/' Complex arithmetic will be used')
 1010 FORMAT (//' CAUTION : No non-zero integrals found on ',A//)
      END
