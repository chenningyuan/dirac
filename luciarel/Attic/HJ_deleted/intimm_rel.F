*
      SUBROUTINE INTIMM_REL(XINT,MAXSYM,IDIM,IPRT)
*
* Import all two electron integrals from MOTRA 2e-file
*
*
      IMPLICIT REAL*8(A-H,O-Z)

      Include 'SysDef.inc'

*.ORBINP
#include "mxpdim.inc"
#include "orbinp.inc"
*.LUNDIO
      Parameter ( mxBatch = 106  )
      Parameter ( mxSyBlk = 666  )
      Common / LundIOR / LuTr2,lTr2Rec,iDAdr(mxBatch),nBatch(mxSyBlk)
*
      dimension IDIM(MAXSYM)
*.Output
      DIMENSION XINT(*)
*. For testing
      Ntestl = 0
      Ntest = max(NTESTL,IPRT)
*
      if (NTEST.ge.1) then
        write(6,*) '--------------------------------'
        write(6,*) ' Information from INTIMM        '
        write(6,*) '--------------------------------'
      end if
*
      Ioff = 1
      Do 101 Ism = 1, Maxsym
        Do 102 Jsm = 1,Ism
          Do 103 Ksm = 1, ISm
            If(Ism .eq. Ksm ) Then
             LsmMX = Jsm
            Else
             LsmMX = Ksm
            End if
            DO 104 Lsm = 1, LsmMX
            If ( ieor(iSm-1,jSm-1).ne.ieor(kSm-1,lSm-1) )  goto 104
*
              IJbl = Ism*(Ism-1)/2 + Jsm
              KLbl = Ksm*(Ksm-1)/2 + Lsm
              IJKLbl = IJbl*(IJbl-1)/2 + KLbl
*
              NIORB = IDIM(ISM)
              NJORB = IDIM(JSM)
              NKORB = IDIM(KSM)
              NLORB = IDIM(LSM)
*
              If(Ism.Eq.Jsm) Then
                nIJ = NIorb*(NIorb+1)/2
              Else
                nIJ = NIORB*NJORB
              End If
              If(Ksm.Eq.Lsm) Then
               nKL = nKorb*(nKorb+1)/2
              Else
                nKL = nKorb*nLorb
              End If
*
              If(Ism .Eq. Ksm .And. Jsm .Eq. Lsm ) Then
                nIJKL= nIJ*(nIJ+1)/2
              Else
                nIJKL= nIJ*nKL
              End if
*
              If ( Ntest .Ne. 0 ) then
                Write(6,*) ' Ism Jsm Ksm Lsm ', Ism,Jsm,Ksm,Lsm
                Write(6,*) ' Ioff, nIJKL ', Ioff,nIJKL
              End if
              Ibatch = nbatch(IJKLbl)
              iDisk=iDAdr(iBatch)
              If(Ntest.Ne.0) Write(6,*)
     &        ' IJKLbl iBatch iDisk ', IJKLbl,iBatch,iDisk
*. Loop over records
              nRec = nIJKL/lTr2Rec
              If(Nrec*lTr2Rec.Lt. nIJKL) nRec = nRec + 1
              Ioffo = Ioff
              Do 50 IRec = 1, Nrec
                If ( IRec .Ne. Nrec ) Then
                   Nintrc = lTr2Rec
                Else
                   Nintrc = nIJKL -(Nrec-1)*Ltr2Rec
                End if
                Call Dafile(LuTr2,2,Xint(Ioff),RtoI*Nintrc,Idisk)
                Ioff = Ioff + Nintrc
   50        Continue
*
             If( Ntest .Ne. 0 ) then
               Write(6,*) ' Integral block '
               If(Ism.Eq.Ksm .And. Jsm. Eq. Lsm ) Then
                 Call Prsym(Xint(Ioffo),nIJ)
               Else
                 Call Wrtmat(Xint(Ioffo),nIJ,nKL,nIJ,nKL)
               End if
             End if
*
  104       Continue
  103     Continue
  102   Continue
  101 Continue
*
      Return
      End
