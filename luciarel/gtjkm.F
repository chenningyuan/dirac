      Subroutine GTJKM(RJ,RK)
*
*     Gather all integrals RJ(I,J) = (II!JJ)
*     Gather all integrals RK(I,J) = (IJ!IJ)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "orbinp.inc"
#include "lucinp.inc"
*.ORBINP
C     PARAMETER(MXPORB=500)
C     PARAMETER(MXPIRR = 20)
C     PARAMETER(MXPOBS = 20)
C     PARAMETER(MXPR4T = 10)
CFUE Addition required to common /ORBINP/
C     PARAMETER(MXTYP  = 10)
C     PARAMETER(MXSYM  = 8 )
CFUE
C     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPR4T,MXPOBS),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB),IBSO(MXPOBS),
C    &              NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB),
C    &              NOBPTS(6+MXPR4T,MXPOBS),IOBPTS(6+MXPR4T,MXPOBS),
C    &              ITOOBS(MXPOBS)
*./LUCINP/ ( NSMOB used )
C     INTEGER PNTGRP,EXTSPC
C     COMMON/LUCINP/PNTGRP,NIRREP,NSMCMP,MAXML,MAXL,
C    &              INTSPC,EXTSPC,NRSSH(MXPIRR,3),
C    &              MNRS1R,MXRS1R,MNRS3R,MXRS3R,NACTEL,
C    &              NSMOB,NRS0SH(1,MXPIRR),NRS4SH(MXPR4T,MXPIRR),
C    &              MXR4TP, MXHR0,MXER4,
C    &              NINASH(MXPIRR),
C    &              INTXCI,NDELSH(MXPIRR),MNRS10,MXRS30
C    Addition required to common /ORBINP/
CFUE
CThis is a safty line to kill the compiler. Note that you have to
Cto some changes
*.LUNDIO
      Parameter ( mxBatch = 106  )
      Parameter ( mxSyBlk = 666  )
      Common / LundIO / LuTr2,lTr2Rec,iDAdr(mxBatch),nBatch(mxSyBlk)
*.Output
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
*.Local
      Parameter ( lBuf    = 9600 )
      Dimension Scr(lBuf)
*
      Do iSym=1,nSmOb
        itOrb=NTOOBS(iSym)
        itOff=ITOOBS(iSym)
        iiBlk=iSym*(iSym+1)/2
        iiPairs=itOrb*(itOrb+1)/2
        Do jSym=1,iSym
          jtOrb=NTOOBS(jSym)
          jtOff=ITOOBS(jSym)
          jjBlk=jSym*(jSym+1)/2
          jjPairs=jtOrb*(jtOrb+1)/2
          ijPairs=itOrb*jtOrb
          If ( iSym.eq.jSym ) ijPairs=jtOrb+itOrb*(itOrb-1)/2
          ijBlk=jSym+iSym*(iSym-1)/2
*
*     collect all RJ(iOrb,jOrb)=(II,JJ)
*
          iRecOld=-1
          iSyBlk=jjBlk+iiBlk*(iiBlk-1)/2
          iBatch=nBatch(iSyBlk)
          iDisk=iDAdr(iBatch)
          nInts=iiPairs*jjPairs
          Do i=1,itOrb
            ii=i*(i+1)/2
*JOS
            MaxJ = jtOrb
            If(Isym.eq.Jsym) MaxJ = i
            Do j=1,MaxJ
*JOE
              jj=j*(j+1)/2
              iijj=ii+(jj-1)*iiPairs
*JOS
              If ( Isym.eq.Jsym ) Then
                jjOff=jj+(jj-1)*(jj-2)/2-1
                iijj =  iijj - jjOff
              End If
*JOE
              iRec=(iijj-1)/lTr2Rec
              If ( iRec.eq.iRecOld ) then
                iijj=iijj-iRec*lTr2Rec
              Else
                iDisk=iDAdr(iBatch)
                Do iSkip=1,iRec
                  Call DaFile(LuTr2,0,Scr,2*lTr2Rec,iDisk)
                End Do
                Call DaFile(LuTr2,2,Scr,2*lTr2Rec,iDisk)
                iijj=iijj-iRec*lTr2Rec
                iRecOld=iRec
              End If
              RJ(i+itOff-1,j+jtOff-1)=Scr(iijj)
              RJ(j+jtOff-1,i+itOff-1)=Scr(iijj)
            End Do
          End Do
*
*     collect all RK(iOrb,jOrb)=(IJ,IJ)
*
          iRecOld=-1
          iSyBlk=ijBlk*(ijBlk+1)/2
          iBatch=nBatch(iSyBlk)
          iDisk=iDAdr(iBatch)
          nInts=ijPairs*(ijPairs+1)/2
          ij=0
          Do i=1,itOrb
            jMax=jtOrb
            If ( iSym.eq.jSym ) jMax=i
            Do j=1,jMax
              ij=ij+1
              ijOff=ij+(ij-1)*(ij-2)/2-1
              ijij=ij+(ij-1)*ijPairs-ijOff
              iRec=(ijij-1)/lTr2Rec
              If ( iRec.eq.iRecOld ) then
                ijij=ijij-iRec*lTr2Rec
              Else
                iDisk=iDAdr(iBatch)
                Do iSkip=1,iRec
                  Call DaFile(LuTr2,0,Scr,2*lTr2Rec,iDisk)
                End Do
                Call DaFile(LuTr2,2,Scr,2*lTr2Rec,iDisk)
                ijij=ijij-iRec*lTr2Rec
                iRecOld=iRec
              End If
              RK(i+itOff-1,j+jtOff-1)=Scr(ijij)
              RK(j+jtOff-1,i+itOff-1)=Scr(ijij)
            End Do
          End Do
*
        End Do
      End Do
*
      RETURN
      END
