      SUBROUTINE GETINM(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,NORBMOL,NIRR_PN)
*
* Obtain integrals XINT(IK,JL) = (IJ!KL) for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "clunit.inc"
*.LUNDIO
      Parameter ( mxBatch = 106  )
      Parameter ( mxSyBlk = 666  )
      Common / LundIO / LuTr2,lTr2Rec,iDAdr(mxBatch),nBatch(mxSyBlk)
      dimension NORBMOL(NIRR_PN)
*.Output
      DIMENSION XINT(*)
*.Local
      Parameter ( lBuf    = 9600 )
      Dimension Scr(lBuf)
      Logical iSymj,kSyml,iSymk,jSyml,ijSymkl,ikSymjl
      Logical ijPerm,klPerm,ijklPerm
*
      iOrb=NOBPTS(ITP,ISM)
      jOrb=NOBPTS(JTP,JSM)
      kOrb=NOBPTS(KTP,KSM)
      lOrb=NOBPTS(LTP,LSM)
      iOff=IOBPTS(ITP,ISM)
      jOff=IOBPTS(JTP,JSM)
      kOff=IOBPTS(KTP,KSM)
      lOff=IOBPTS(LTP,LSM)
*
      IAD13 = 0
      NLONGER = 2
      call findhigh(NORBMOL,NIRR_PN,MAXDMMO)
      write(6,*) 'Highest symmetry dimension is  ',MAXDMMO
      KBUFF = MAXDMMO**4
*
*     Collect Coulomb terms
*
      ijPerm=.false.
      klPerm=.false.
      ijklPerm=.false.
      If ( iSm.gt.jSm ) Then
        iSym=iSm
        jSym=jSm
      Else
        iSym=jSm
        jSym=iSm
        ijPerm=.true.
      End If
      ijBlk=jSym+iSym*(iSym-1)/2
      If ( kSm.gt.lSm ) Then
        kSym=kSm
        lSym=lSm
      Else
        kSym=lSm
        lSym=kSm
        klPerm=.true.
      End If
      klBlk=lSym+kSym*(kSym-1)/2
      If ( klBlk.gt.ijBlk ) Then
        iTemp=iSym
        iSym=kSym
        kSym=iTemp
        iTemp=jSym
        jSym=lSym
        lSym=iTemp
        iTemp=ijBlk
        ijBlk=klBlk
        klBlk=iTemp
        ijklPerm=.true.
      End If
C?    print *,' GetInt reports to you'
C?    print *,' iSym,jSym,kSym,lSym',iSym,jSym,kSym,lSym
      iSyBlk=klBlk+ijBlk*(ijBlk-1)/2
      iBatch=nBatch(iSyBlk)
C?    print *,' iSyBlk,iBatch',iSyBlk,iBatch
      iSymj=iSym.eq.jSym
      kSyml=kSym.eq.lSym
      iSymk=iSym.eq.kSym
      jSyml=jSym.eq.lSym
      ijSymkl=iSymj.and.kSyml
      ikSymjl=iSymk.and.jSyml
      itOrb=NTOOBS(iSym)
      jtOrb=NTOOBS(jSym)
      ktOrb=NTOOBS(kSym)
      ltOrb=NTOOBS(lSym)
C?    print *,' itOrb,jtOrb,ktOrb,ltOrb',itOrb,jtOrb,ktOrb,ltOrb
      ijPairs=itOrb*jtOrb
      klPairs=ktOrb*ltOrb
      nInts=ijPairs*klPairs
      If ( ikSymjl ) Then
        If ( iSymj ) Then
          ijPairs=itOrb*(itOrb+1)/2
          klPairs=ijPairs
          nInts=ijPairs*(klPairs+1)/2
        Else
          ijPairs=itOrb*jtOrb
          klPairs=ijPairs
          nInts=ijPairs*(klPairs+1)/2
        End If
      Else If ( ijSymkl ) Then
        ijPairs=itOrb*(itOrb+1)/2
        klPairs=ktOrb*(ktOrb+1)/2
        nInts=ijPairs*klPairs
      End If
C?    print *,' ijPairs,klPairs',ijPairs,klPairs
      iRecOld=-1
      iInt=0
C?    print *,' Start loop over Coulomb integrals'
C?    print *,' ijkl,ij,kl,i,j,k,l'
      Do lJeppe=lOff,lOff+lOrb-1
        jMin=jOff
        If ( JLSM.ne.0 ) jMin=lJeppe
        Do jJeppe=jMin,jOff+jOrb-1
          Do kJeppe=kOff,kOff+kOrb-1
            iMin = iOff
            If(IKSM.ne.0) iMin = kJeppe
            Do iJeppe=iMin,iOff+iOrb-1
              If ( ijPerm ) Then
                ii=jJeppe
                jj=iJeppe
              Else
                ii=iJeppe
                jj=jJeppe
              End If
              If ( klPerm ) Then
                kk=lJeppe
                ll=kJeppe
              Else
                kk=kJeppe
                ll=lJeppe
              End If
              If ( ijklPerm ) Then
                i=kk
                j=ll
                k=ii
                l=jj
              Else
                i=ii
                j=jj
                k=kk
                l=ll
              End If
              If ( iSymj ) Then
                If ( i.gt.j ) Then
                  ij=j+i*(i-1)/2
                Else
                  ij=i+j*(j-1)/2
                End If
              Else
                ij=j+(i-1)*jtOrb
              End If
              ijOff=ij+(ij-1)*(ij-2)/2-1
              If ( kSyml ) Then
                If ( k.gt.l ) Then
                  kl=l+k*(k-1)/2
                Else
                  kl=k+l*(l-1)/2
                End If
              Else
                kl=l+(k-1)*ltOrb
              End If
              klOff=kl+(kl-1)*(kl-2)/2-1
              If ( ikSymjl ) Then
                If ( ij.gt.kl ) Then
                  ijkl=ij+(kl-1)*ijPairs-klOff
                Else
                  ijkl=kl+(ij-1)*klPairs-ijOff
                End If
              Else
                ijkl=ij+(kl-1)*ijPairs
              End If
C?            print '(7i6)',ijkl,ij,kl,i,j,k,l
              iRec=(ijkl-1)/lTr2Rec
              If ( iRec.eq.iRecOld ) then
                ijkl=ijkl-iRec*lTr2Rec
              Else
                iDisk=iDAdr(iBatch)
                Do iSkip=1,iRec
                  Call DaFile(LuTr2,0,Scr,2*lTr2Rec,iDisk)
                End Do
                Call DaFile(LU2INT,2,Scr,KBUFF*NLONGER,IAD13)
                ijkl=ijkl-iRec*lTr2Rec
                iRecOld=iRec
              End If
              iInt=iInt+1
              XInt(iInt)=Scr(ijkl)
            End Do
          End Do
        End Do
      End Do
*
*     Collect Exchange terms
*
      If ( IXCHNG.ne.0 ) Then
        ijPerm=.false.
        klPerm=.false.
        ijklPerm=.false.
        If ( iSm.gt.lSm ) Then
          iSym=iSm
          jSym=lSm
        Else
          iSym=lSm
          jSym=iSm
          ijPerm =.True.
        End If
        ijBlk=jSym+iSym*(iSym-1)/2
        If ( kSm.gt.jSm ) Then
          kSym=kSm
          lSym=jSm
        Else
          kSym=jSm
          lSym=kSm
          klPerm =.True.
        End If
        klBlk=lSym+kSym*(kSym-1)/2
        If ( klBlk.gt.ijBlk ) Then
          iTemp=iSym
          iSym=kSym
          kSym=iTemp
          iTemp=jSym
          jSym=lSym
          lSym=iTemp
          iTemp=ijBlk
          ijBlk=klBlk
          klBlk=iTemp
          ijklPerm = .True.
        End If
        iSyBlk=klBlk+ijBlk*(ijBlk-1)/2
        iBatch=nBatch(iSyBlk)
        iSymj=iSym.eq.jSym
        kSyml=kSym.eq.lSym
        iSymk=iSym.eq.kSym
        jSyml=jSym.eq.lSym
        ijSymkl=iSymj.and.kSyml
        ikSymjl=iSymk.and.jSyml
        itOrb=NTOOBS(iSym)
        jtOrb=NTOOBS(jSym)
        ktOrb=NTOOBS(kSym)
        ltOrb=NTOOBS(lSym)
        ijPairs=itOrb*jtOrb
        klPairs=ktOrb*ltOrb
        nInts=ijPairs*klPairs
        If ( ikSymjl ) Then
          If ( iSymj ) Then
            ijPairs=itOrb*(itOrb+1)/2
            klPairs=ijPairs
            nInts=ijPairs*(klPairs+1)/2
          Else
            ijPairs=itOrb*jtOrb
            klPairs=ijPairs
            nInts=ijPairs*(klPairs+1)/2
          End If
        Else If ( ijSymkl ) Then
          ijPairs=itOrb*(itOrb+1)/2
          klPairs=ktOrb*(ktOrb+1)/2
          nInts=ijPairs*klPairs
        End If
        iRecOld=-1
        iInt=0
C?      print *,' iSym,jSym,kSym,lSym',iSym,jSym,kSym,lSym
C?      print *,' iSyBlk,iBatch',iSyBlk,iBatch
C?      print *,' itOrb,jtOrb,ktOrb,ltOrb',itOrb,jtOrb,ktOrb,ltOrb
C?      print *,' ijPairs,klPairs',ijPairs,klPairs
C?      print *,' Start loop over Exchange integrals'
C?      print *,' ijkl,ij,kl,i,j,k,l'
        Do lJeppe=lOff,lOff+lOrb-1
          jMin=jOff
          If ( JLSM.ne.0 ) jMin=lJeppe
          Do jJeppe=jMin,jOff+jOrb-1
            Do kJeppe=kOff,kOff+kOrb-1
              iMin = iOff
              If(IKSM.ne.0) iMin = kJeppe
              Do iJeppe=iMin,iOff+iOrb-1
              If ( ijPerm ) Then
                ii=lJeppe
                ll=iJeppe
              Else
                ii=iJeppe
                ll=lJeppe
              End If
              If ( klPerm ) Then
                kk=jJeppe
                jj=kJeppe
              Else
                kk=kJeppe
                jj=jJeppe
              End If
              If ( ijklPerm ) Then
                i=kk
                l=jj
                k=ii
                j=ll
              Else
                i=ii
                j=jj
                k=kk
                l=ll
              End If
                If ( iSymj ) Then
                  If ( i.gt.l ) Then
                    ij=l+i*(i-1)/2
                  Else
                    ij=i+l*(l-1)/2
                  End If
                Else
                  ij=l+(i-1)*jtOrb
                End If
                ijOff=ij+(ij-1)*(ij-2)/2-1
                If ( kSyml ) Then
                  If ( k.gt.j ) Then
                    kl=j+k*(k-1)/2
                  Else
                    kl=k+j*(j-1)/2
                  End If
                Else
                  kl=j+(k-1)*ltOrb
                End If
                klOff=kl+(kl-1)*(kl-2)/2-1
                If ( ikSymjl ) Then
                  If ( ij.gt.kl ) Then
                    ijkl=ij+(kl-1)*ijPairs-klOff
                  Else
                    ijkl=kl+(ij-1)*klPairs-ijOff
                  End If
                Else
                  ijkl=ij+(kl-1)*ijPairs
                End If
C?              print '(7i6)',ijkl,ij,kl,i,j,k,l
                iRec=(ijkl-1)/lTr2Rec
                If ( iRec.eq.iRecOld ) then
                  ijkl=ijkl-iRec*lTr2Rec
                Else
                  iDisk=iDAdr(iBatch)
                  Do iSkip=1,iRec
                    Call DaFile(LU2INT,0,Scr,2*lTr2Rec,iDisk)
                  End Do
                  Call DaFile(LU2INT,2,Scr,KBUFF*NLONGER,IAD13)
                  ijkl=ijkl-iRec*lTr2Rec
                  iRecOld=iRec
                End If
                iInt=iInt+1
                XInt(iInt)=XInt(iInt)-Scr(ijkl)
              End Do
            End Do
          End Do
        End Do
      End If
*
      Return
      End
