      Subroutine Getobs2(ECORE,NAOS_ENV,NMOS_ENV)
************************************************************************
*                                                                      *
*     Purpose:                                                         *
*     Open one-electron integral file
*.    Read in orbital information  and core energy
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************
*
      Implicit Real*8 (A-H,O-Z)
*
      Include 'SysDef.inc'
*
      Parameter( mxOrb = 510 )
      Parameter( LuOne = 14  )
*
      DIMENSION NAOS_ENV(*),NMOS_ENV(*)
*
      COMMON/MOLOBS/
     & IOList(64),iToc(64),nBas(8),nOrb(8),nFro(8),nDel(8),
     & Nsym
      Character*4 Name(2,mxOrb)
*----------------------------------------------------------------------*
*     Start procedure:                                                 *
*     open the transformed one-electron integral file                  *
*----------------------------------------------------------------------*
      LuX = LuOne
      Call DaName(LuX,'TRAONE')
*----------------------------------------------------------------------*
*     Set up the scatter/gather list and read table of contents        *
*----------------------------------------------------------------------*
      iDisk=0
      Call GSList(IOLIST,8,iToc,64,ECORE,RtoI,nSym,1,
     &            nBas,8,nOrb,8,nFro,8,nDel,8,Name,2*4*mxOrb/ItoB)
      Call DaFile(LuX,4,IOList,iDum,iDisk)
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
*
* Copy to external arrays
*
      CALL ICOPVE(NORB,NMOS_ENV,NSYM)
      CALL ICOPVE(NBAS,NAOS_ENV,NSYM)
*
      Return
      End
