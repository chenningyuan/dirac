      Subroutine MkLundIO
************************************************************************
*                                                                      *
*     Purpose:                                                         *
*     Initialize the Common /LundIO/                                   *
*                                                                      *
*     Calling parameters: none                                         *
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************
*
      Parameter ( mxBatch = 106  )
      Parameter ( mxSyBlk = 666  )
      Parameter ( lBlk    = 9600 )
      Parameter ( LuTwo   = 13   )
      Common / LundIO / LuTr2,lTr2Rec,iDAdr(mxBatch),nBatch(mxSyBlk)
*----------------------------------------------------------------------*
*     Start procedure:                                                 *
*     First set the unit number, record length and open file           *
*     per symmetry element                                             *
*----------------------------------------------------------------------*
      LuTr2=LuTwo
      lTr2Rec=lBlk
      Call DaName(luTr2,'TRAINT')
*----------------------------------------------------------------------*
*     Load the table of disk adresses                                  *
*----------------------------------------------------------------------*
      Call DaFile(LuTr2,2,iDAdr,mxBatch,iDisk)
*----------------------------------------------------------------------*
*     Generate the symmetry block to batch number translation table    *
*----------------------------------------------------------------------*
      iBatch=0
      Do iSym=1,8
        Do jSym=1,iSym
          Do kSym=1,iSym
            mxlSym=kSym
            If ( kSym.eq.iSym ) mxlSym=jSym
            Do lSym=1,mxlSym
              If ( ieor(iSym-1,jSym-1).eq.ieor(kSym-1,lSym-1) ) Then
                ijPair=jSym+iSym*(iSym-1)/2
                klPair=lSym+kSym*(kSym-1)/2
                iSyBlk=klPair+ijPair*(ijPair-1)/2
                iBatch=iBatch+1
                nBatch(iSyBlk)=iBatch
              End If
            End Do
          End Do
        End Do
      End Do
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
      Return
      End
