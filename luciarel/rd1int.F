*
      subroutine rd1int(ONEINT,FNONEAO,LUONEAO,NATOMS,IFOCKEX)
************************************************************************
*                                                                      *
*     Objective: Read the header of the one-electron integral file     *
*                Extract also symmetry and basis set information.      *
*                In addition read the overlap, the nuclear attraction  *
*                and kinetic integrals.                                *
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************
*                                                                      *
*     modified for LUCIAREL adapted to MOLCAS (4.1, 4.3 running)       *
*                                                                      *
*     Timo Fleig, 16.05.99    (another sunny day in may; a sun-day)    *
*                              may-day ?
*                 25-27.05.99 (cont'd, more sunny days in may ...)     *
*                 30.11.99    (read nucl. rep. en. when INPORB ..)     *
*                             (A fucking rainy and windy november day) *
*                                                                      *
************************************************************************
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "symm.inc"
#include "onel.inc"
#include "input.inc"
*
      common/CECORE/ECORE
*
      character*6 FNONEAO
      character*8 OneLbl,ComLbl
      character*72 HEADER
*
      dimension NBAS_SCF(MXNSPIR)
      dimension COOR_SCF(3,MXNATOMS),ONEINT(ITRIASUM)
      dimension HONE_SCF(MXPORB*MXPORB)
*
      NTESTL = 0
      NTEST = max(IPRINT,NTESTL)
*
*  Initialize
      do IINT=1,ITRIASUM,1
        ONEINT(IINT) = 0.0
      end do
*
      Call qEnter('rd1int')
*----------------------------------------------------------------------*
*     Start procedure                                                  *
*     Open the one electron integral file                              *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      Call OpnOne(iRC,iOpt,FNONEAO,LUONEAO)
      If ( iRc.ne.0 ) Goto 990
*----------------------------------------------------------------------*
*     Read one-electron integral file header etc.                      *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      iComp=0
      iSyLbl=1
      OneLbl='Title   '
      Call RdOne(iRc,iOpt,OneLbl,iComp,HEADER,iSyLbl)
      If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read no.of symm. species                                         *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      iComp=0
      iSyLbl=1
      OneLbl='nSym    '
      Call RdOne(iRc,iOpt,OneLbl,iComp,NSYM_SCF,iSyLbl)
      if (NSYM_SCF.ne.NIRR_PN) then
        write(6,*) 'Number of symms on comfile does not'
        write(6,*) 'equal those from input.'
        Call Abend2( 'severe error, quitting.' )
      else
        write(6,*) 'Number of pg symmetries: ',NSYM_SCF
      end if
      If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read symm. oper per symm. species                                *
*----------------------------------------------------------------------*
c     iRc=-1
c     iOpt=0
c     iComp=0
c     iSyLbl=1
c     OneLbl='SymOp   '
c     Call RdOne(iRc,iOpt,OneLbl,iComp,IOPER_SCF,iSyLbl)
c     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read no. of basis functions per symm. species                    *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      iComp=0
      iSyLbl=1
      OneLbl='nBas    '
      Call RdOne(iRc,iOpt,OneLbl,iComp,NBAS_SCF,iSyLbl)
      do ISYM=1,NIRR_PN,1
        write(6,*) 'Dimension of symm ',ISYM,' is ',NBAS_SCF(ISYM)
        if (SYDI(ISYM).ne.NBAS_SCF(ISYM)) then
          write(6,*) 'Dimensions of spatial symmetries corrupted'
          write(6,'(A5,I1,A4,I3)')
     &             'SYDI(',ISYM,') = ',SYDI(ISYM)
          write(6,'(A9,I1,A4,I3)')
     &             'NBAS_SCF(',ISYM,') = ',NBAS_SCF(ISYM)
          Call Abend2( 'fatal error in rd1int' )
        end if
      end do
      If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read no. of basis functions per symm. species                    *
*----------------------------------------------------------------------*
c     iRc=-1
c     iOpt=0
c     iComp=0
c     iSyLbl=1
c     OneLbl='AtomType'
c     Call RdOne(iRc,iOpt,OneLbl,iComp,BsLbl,iSyLbl)
c     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read no. of unique atoms in the system                           *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      iComp=0
      iSyLbl=1
      OneLbl='nAtoms  '
      Call RdOne(iRc,iOpt,OneLbl,iComp,NATOMS,iSyLbl)
      write(6,*) 'Number of unique atoms:  ',NATOMS
      if (NATOMS.gt.MXNATOMS) then
        write(6,*) 'Number of atoms is too large.'
        write(6,*) 'Adjust MXNATOMS in mxpdim.inc or reconsider.'
        Call Abend2( 'Quitting.' )
      end if
      If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read atom labels                                                 *
*----------------------------------------------------------------------*
c     iRc=-1
c     iOpt=0
c     iComp=0
c     iSyLbl=1
c     OneLbl='AtomLbl '
c     Call RdOne(iRc,iOpt,OneLbl,iComp,AtLbl,iSyLbl)
c     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read coordinates of atoms                                        *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      iComp=0
      iSyLbl=1
      OneLbl='Coord   '
      Call RdOne(iRc,iOpt,OneLbl,iComp,COOR_SCF,iSyLbl)
      do IATM=1,NATOMS,1
        write(6,'(A9,I3,A1,X,3F10.6)') 'Atom no. ',IATM,':',
     &                                 (COOR_SCF(IC,IATM),IC=1,3,1)
      end do
      If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read nuclear repulsion energy                                    *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      iComp=0
      iSyLbl=1
      OneLbl='PotNuc  '
      Call RdOne(iRc,iOpt,OneLbl,iComp,POTNUC_SCF,iSyLbl)
      write(6,'(A,X,F18.15)') 'Nuclear repulsion energy: ',POTNUC_SCF
      If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read overlap integrals                                           *
*----------------------------------------------------------------------*
c     iRc=-1
c     iOpt=6
c     iComp=1
c     iSyLbl=1
c     OneLbl='Mltpl  0'
c     Call RdOne(iRc,iOpt,OneLbl,iComp,HONE_SCF,iSyLbl)
c     write(6,*)
c     write(6,*) 'Overlap integrals:'
c     do IINT=1,ITRIASUM,1
c       write(6,'(A,X,I4,A,F18.14)')
c    &        'Symm orb integral ',IINT,':',HONE_SCF(IINT)
c     end do
c     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read core Hamiltonian                                            *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=6
      iComp=1
      iSyLbl=1
      OneLbl='OneHam  '
      Call RdOne(iRc,iOpt,OneLbl,iComp,HONE_SCF,iSyLbl)
      write(6,*)
      write(6,*) 'Reading core Hamiltonian'
      write(6,*) 'from MOLCAS file ONEINT ...'
      write(6,*)
*  Copy to SO integral storage array
      do IINT=1,ITRIASUM,1
        ONEINT(IINT) = ONEINT(IINT) + HONE_SCF(IINT)
      end do
*  Check
      if (NTEST.ge.1) then
        do IINT=1,ITRIASUM,1
          write(6,'(A,X,I4,A,F20.12)')
     &          'Symm orb integral ',IINT,':',HONE_SCF(IINT)
        end do
      end if
      If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read nuclear attraction integrals                                *
*----------------------------------------------------------------------*
c     iRc=-1
c     iOpt=6
c     iComp=1
c     iSyLbl=1
c     OneLbl='Attract '
c     Call RdOne(iRc,iOpt,OneLbl,iComp,HONE_SCF,iSyLbl)
c     write(6,*)
c     write(6,*) 'Nuclear attraction integrals:'
c     do IINT=1,ITRIASUM,1
c       write(6,'(A,X,I4,A,F20.12)')
c    &        'Symm orb integral ',IINT,':',HONE_SCF(IINT)
c     end do
c     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read kinetic energy integrals                                    *
*----------------------------------------------------------------------*
c     iRc=-1
c     iOpt=6
c     iComp=1
c     iSyLbl=1
c     OneLbl='Kinetic '
c     Call RdOne(iRc,iOpt,OneLbl,iComp,HONE_SCF,iSyLbl)
c     write(6,*)
c     write(6,*) 'kinetic energy integrals:'
c     do IINT=1,ITRIASUM,1
c       write(6,'(A,X,I4,A,F20.12)')
c    &        'Symm orb integral ',IINT,':',HONE_SCF(IINT)
c     end do
c     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     If this is a perturbative reaction field calculation then        *
*     modify the one-electron Hamiltonian by the reaction field and    *
*     the nuclear attraction by the cavity self-energy                 *
*----------------------------------------------------------------------*
c     If ( iRFpert.ne.0 ) then
c        nTemp=0
c        Do iSym=1,nSym
c           nTemp=nTemp+nBas(iSym)*(nBas(iSym)+1)/2
c        End Do
c        Call GetMem('RFFLD','Allo','Real',lTemp,nTemp)
c        iRc=-1
c        iOpt=0
c        Call OpnCom(iRc,iOpt,FnCom,LuCom)
c        If ( iRc.ne.0 ) Goto 993
c        iRc=-1
c        iOpt=0
c        ComLbl='ERFself'
c        Call RdCom(iRc,iOpt,ComLbl,ERFself,2)
c        PotNuc=PotNuc+ERFself
c        If ( iRc.ne.0 ) Goto 994
c        iRc=-1
c        iOpt=0
c        ComLbl='RFFLD'
c        Call RdCom(iRc,iOpt,ComLbl,Work(lTemp),2*nTemp)
c        Call Daxpy(nTemp,1.0D0,Work(lTemp),1,Hone,1)
c        If ( iRc.ne.0 ) Goto 994
c        iRc=-1
c        iOpt=0
c        Call ClsCom(iRc,iOpt)
c        If ( iRc.ne.0 ) Goto 995
c        Call GetMem('RFFLD','Free','Real',lTemp,nTemp)
c     End If
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      call clsone(iRc,iOpt)
      If ( iRc.ne.0 ) Goto 992
*
c     nTot=0
c     nTot1=0
c     nTot2=0
c     n2max=0
c     Do iSym=1,nSym
c       iBas=nBas(iSym)
c       nTot=nTot+iBas
c       nTot1=nTot1+iBas*(iBas+1)/2
c       nTot2=nTot2+iBas*iBas
c       n2max=Max(n2max,iBas*iBas)
c     End Do
      Call qExit('Rd1Int')
*
      if (IFOCKEX.eq.0.and.ORBFILE) then
        ECORE = POTNUC_SCF
        write(6,*) 'Passing POTNUC to ECORE'
      end if
*
      write(6,*) 'file opened successfully, returning ...'
      Return
*----------------------------------------------------------------------*
*     Error Exit                                                       *
*----------------------------------------------------------------------*
990   Call ErrMsg('Rd1Int',001)
      Call Quit(16)
991   Call StkA(0,OneLbl)
      Call ErrMsg('Rd1Int',002)
      Call Quit(16)
992   Call ErrMsg('Rd1Int',003)
      Call Quit(16)
993   Call ErrMsg('Rd1Int',004)
      Call Quit(16)
994   Call StkA(0,ComLbl)
      Call ErrMsg('Rd1Int',005)
      Call Quit(16)
995   Call ErrMsg('Rd1Int',006)
      Call Quit(16)
      End
