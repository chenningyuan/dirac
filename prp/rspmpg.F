C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpxlr_mpg */
      SUBROUTINE PRPXLR_MPG(DMO,XMO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by C.V. Larsen
C     Last revision J. van Stralen - June 2002
C
C     PURPOSE: Driver routine for the calculation of the V-O and P-O 
C              part of the D^2 matrix. It starts by setting up parameters
C              that the normal response program gets from the input file 
C              or default settings.This includes setting up an operator.
C
C     INPUT:   DMO - second order density matrix 
C              XMO - the Lagrangian
C
C     REMARKS: -
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbxlr.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dummy.h"
      dimension work(*)
      DIMENSION DMO(*)
      DIMENSION XMO(*)
      logical orbdef(3)
      dimension pfac(3)
      character plabel(3)*8
      character tmp(4)*4

      call qenter('PRPXLR_MPG')
C
C     Initialization - from xlrinp
C     ============================
C
      illint = 1
      islint = 1
      issint = 0
C
C     Initialize /XCBLR/
C     ==================
C
      iprxlr    = 0            !general printlevel
      thclr     = 1.d-7         !maximum norm of solutionerror vector
      thrcom    = 1.0d-6         !threshold for compression
      resxlr    = 1.d3
      cnvxlr(1) = dummy
      cnvxlr(2) = 1.d-2
      MAXRM     = 200           
      itrxlr    = 100
      intxlr    = illint+2*islint+4*issint
      if(levyle) intxlr = illint
      itrilr(1) = 1
      itrilr(2) = 1
      nlrapt    = 0             !number of A operators
      nlrbpt    = 0             !number of B operators
      call izero(nlrap,8)
      call izero(nlrbp,8)
      nbfreq    = 1             !number of frequencies
      bfreq(1)  = 0.d0
      allcmb    = .false.       !generate all operator combinations
      xlr_skipep= .false.       !do include e-p orbital rotations
      xlr_skipee= .false.       !do include e-e orbital rotations
c     qedrun    = .TRUE.        !not implemented
      xlr_lsfg(1) = .true.        !Generate G-matrix in MO basis
      xlr_lsfg(2) = .true.        !Generate G-matrix in MO basis
      xlrnrm    = .true.        !normalize solutionvector
      xLRDIH    = .TRUE.        !pre-conditioning of trail vectors
C                                
      orbdef(1) = .true.
      orbdef(2) = .true.
      orbdef(3) = .true.
      triab     = .false.       !force triangularity
      XLR_XVCFIL    = 'NORXVC'
      orbxlr    = .false.       !is this an orbital response problem
      xstern    = .false.
C
C     set work array to zero
C     ======================
C
      CALL DZERO(WORK,LFREE)
C
C     Initialize VIRTUAL orbital strings
C     ==================================
C
C     Virtual input strings
C
      do I=1,nfsym
         call num2str(nocc(I),tmp(1))
         call num2str(npsh(I),tmp(2))
         call num2str(nocc(I)+1,tmp(3))
         call num2str(nesh(I),tmp(4))
         xlr_indstr(1,I)='1..'//tmp(1)
         xlr_indstr(2,I)='0'
         if(levyle)then
            xlr_indstr(3,I)=tmp(3)//'..'//tmp(4)
         else
            xlr_indstr(3,I)='-'//tmp(2)//'..-1,'//tmp(3)//'..'//tmp(4)
         endif
      enddo
C
C     Initialize the operators needed
C     ===============================
C
      indxpr=-1
      isyxpr=-1
      itrxpr=-1
C
      iptyp=1
      npcomp=1
C
C     ---Right side operator
C
      pfac(1)=-1.d0
      plabel(1)='OVERLAP'
      write(lupri,'(a,6i4)') 'XPRIND',iptyp,npcomp,indxpr,isyxpr,itrxpr,
     &     iprxlr
      call xprind('Lagrangian rhs. ',iptyp,npcomp,pfac,plabel,
     &     indxpr,isyxpr,itrxpr,work(kfree),lfree,iprxlr)
      write(lupri,'(a,6i4)') 'XPRIND',iptyp,npcomp,indxpr,isyxpr,itrxpr,
     &     iprxlr
      call op1ind('NLRBPT',ind1op,llrbpu,nlrbpt,indxpr,maxllr)
C
C     Allocate memory
C     ===============
C     
C{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
      ifrp = 1
C{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
      nrspf=nlrbpt*nbfreq
      if(orbxlr)then
         call memget('INTE',klro,2*norb(ifrp),work,kfree,lfree)
         call deforb(work(klro),ncls_xlr,3,xlr_indstr,
     &        work,kfree,lfree,iprxlr)
         ni=0
         do i=1,nfsym
            ni=ni+ncls_xlr(1,1,i)+ncls_xlr(2,1,i)
         enddo
      else
         ni=0
         klrorb=kfree
      endif
      if(xlr_skipee)then
         nee=0
      else
         nee=ni+1
      endif
      if(xlr_skipep)then
         nep=0
      else
         nep=ni+1
      endif
      nci=0
      n2lrf=nlrapt*nrspf*(1+nee+nep+nci)
      call memget('REAL',klrf,n2lrf,work,kfree,lfree)
      call dzero(work(klrf),n2lrf)
      call memget('REAL',kcnv,nrspf,work,kfree,lfree)
C
C
      CALL FLSHFO(LUPRI)
C
      call prpxl1_mpg(work(klrf),work(kcnv),work(klro),DMO,XMO,
     &                work,kfree,lfree)
C
      call qexit('PRPXLR_MPG')
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpxl1_mpg */
      SUBROUTINE PRPXL1_MPG(ABLRF,RCNV,IORBCL,DMO,XMO,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Written by C.V. Larsen
C     Last revision J. van Stralen - June 2002
C
C     INPUT:   ABLRF  - something for the linear response function,
C                       not needed for the MP2 gradient
C              RCNV   - total norm of residual relative to norm of
C                       solution vector
C              IORBCL - something for the linear response function,
C                       not needed for the MP2 gradient
C              XMO   -  the Lagrangian
C     
C     REMARKS: -
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "aovec.h"
#include "mxorb.h"
#include "infpar.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbxlr.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbxrs.h"
C
#include "blocks.h"
      parameter(gpthrs=1.d-9)
      dimension ablrf(nfreq,nlrapt,nlrbpt,*)
      dimension rcnv(nbfreq,nlrbpt)
      dimension iorbcl(2,*)
      DIMENSION DMO(*)
      DIMENSION XMO(*)
      DIMENSION WORK(*)
C
      call qenter('PRPXL1_MPG')
C
      kfrsav=kfree
C
C     Initialize
C     ==========
C
      call setxlr_mpg(work,kfree,lfree)
C
C     Linear response function
C     ========================
C
C
C     Dimensions of reduced system
C     ----------------------------
C
      kxlr=kfree
      nevecr=nredm*nexsim
      call memget('INTE',kibtyp,2*nredm,work,kfree,lfree)
      call memget('INTE',kibcvc,  nredm,work,kfree,lfree)
      call memget('INTE',kibevc,  nredm,work,kfree,lfree)
      call memget('INTE',kibpvc,  nredm,work,kfree,lfree)
      call memget('REAL',kevalr, nexsim,work,kfree,lfree)
      call memget('REAL',kevecr, nevecr,work,kfree,lfree)
C
C     Transfer frequencies (actually only one frequency needed 0)
C     -----------------------------------------------------------
C
      call dcopy(nfreq,bfreq,1,work(kevalr),1)
C
C     Loop over operators
C     -------------------
C
      do jopsy=1,nfsym
         nfc=nbsym/nfsym
         do is=1,nfc
            jsymop=jfsym(is,jopsy)
            nopb=nlrbp(jsymop)
            jopb=jlrbp(jsymop)
            do 30 iopb=1,nopb
               indbp=llbpsu(jopb+iopb)
               indprb=llrbpu(indbp)
               jtimop=iprptim(indprb)
C
C     Initialize orbital parameters
C     =============================
C
               call xrspar(xlr_indstr,xlr_skipee,xlr_skipep,
     &                     work,kfree,lfree,IPRINT)
C
C     Property gradient
C     =================
C
               kgrad=kfree
               call memget('REAL',kgpci,nzconfq,work,kfree,lfree)
               call memget('REAL',kgpoe,nzxopeq,work,kfree,lfree)
               call memget('REAL',kgpop,nzxoppq,work,kfree,lfree)
               call pamgrd_mpg(indprb,work(kgpoe),
     &                         work(kgpop),jopsy,work(kzconf),
     &                         work(kzxope),work(kzxopp),nzconf,nzxope,
     &                         nzxopp,XMO,work,kfree,lfree,iprxlr)
C
C     Check norm of gradient
C     ======================
C
               gpcnrm=dnrm2(nzconfq,work(kgpci),1)
               gpenrm=dnrm2(nzxopeq,work(kgpoe),1)
               gppnrm=dnrm2(nzxoppq,work(kgpop),1)
               gptnrm=sqrt(gpcnrm**2+gpenrm**2+gpprnm**2)
c               write(lupri,'(a,f18.10)') 'Gnorm',gptnrm
               if(gptnrm.lt.gpthrs)then
                  goto 30
               endif
C
C     Solve response equation
C     =======================
C
               ncred=0
               nered=0
               npred=0
               nzred=ncred+nered+npred
               RSPLAB = PRPNAM(INDPRB)
               call xrsctl(work(kgpci),work(kgpoe),work(kgpop),
     &                     work(kibtyp),work(kibcvc),work(kibevc),
     &                     work(kibpvc),rcnv(1,indbp),work(kevalr),
     &                     work(kevecr),work,kfree,lfree)
               call memrel('PRPGR1.XRSCTL',work,1,kgrad,kfree,lfree)
C
C     Evaluate result
C     ===============
C
      KEIG = KFREE
               call getlrf_mpg(ablrf,work(kibtyp),work(kibcvc),
     &                         work(kibevc),work(kibpvc),RCNV(1,INDBP),
     &                         work(kevalr),work(kevecr),iorbcl,
     &                         work(keig),DMO,work,kfree,lfree)
 30         continue
         enddo
      enddo
      call memrel('PRPGR1.ERR',work,1,kxlr,kfree,lfree)
C
C     Final remarks
C     =============
C
      close(luxvc,status='KEEP')
      call memrel('PRPGR1.FINAL',work,1,kfrsav,kfree,lfree)
C
      call qexit('PRPXL1_MPG')
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setxlr_mpg */
      SUBROUTINE SETXLR_MPG(WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by C.V. Larsen
C     Last revision J. van Stralen - June 2002
C
C     Initialize variables/Open necessary files etc.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"

#include "dgroup.h"
#include "dcbxlr.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbxrs.h"
#include "cbihr2.h"

      logical tobe
      dimension work(*)
      call qenter('SETXLR_MPG')
      kfrsav=kfree
C
C     Sort operators on symmetries
C     ============================
C
      call op1srt('LLBPSU',llrbpu,llbpsu,nlrbpt,nlrbp,jlrbp,
     &             work(kfree),lfree)
C
C     Transfer information
C     ====================
C
      lineq     = .true.
      lsvcfg(1) = xlr_lsfg(1)
      lsvcfg(2) = xlr_lsfg(2)
      tknorm    = xlrnrm
      diaghe    = xlrdih
      iprxrs    = iprxlr
      thcxrs    = thclr
      resfac    = resxlr
      maxitr    = itrxlr
      nfreq     = nbfreq
      nexsim    = nfreq
      nexstv    = nfreq
      nexcnv    = nfreq
      nredm     = maxrm
      n2redm    = maxrm*maxrm
      loffty    = 0
      static    = .true.
      cnvint(1) = cnvxlr(1)
      cnvint(2) = cnvxlr(2)
      itrint(1) = itrilr(1)
      itrint(2) = itrilr(2)
      intdef    = intxlr
      sternh    = xstern
C
C     Set printlevel in two-electron routines
C     =======================================
C
      iprtwo=0
C
C     initialize file for solution vector
C     ===================================
C
      luxvc=29
      call opnfil(luxvc,'PAMXVC','UNKNOWN','PRPGRA')
      WRITE(LUXVC) 'END_OF_THIS_FILE'
C
      call qexit('SETXLR_MPG')
C
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamgrd_mpg */
      SUBROUTINE PAMGRD_MPG(INDXPR,GPOE,GPOP,JOPSY,JCONF,JXOPE,
     &                      JXOPP,NCONF,NXOPE,NXOPP,XMO,WORK,KFREE,
     &                      LFREE,IPRINT)
C***********************************************************************
C
C     Written by C.V. Larsen
C     Last revision J. van Stralen - June 2002
C
C     PURPOSE:
C
C     INPUT:   NXOPE  - # of e-e orbital rotations
C              NXOPP  - # of e-p orbital rotations
C
C     REMARKS: -
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbxpr.h"
      dimension gpoe(*)
      dimension gpop(*)
      dimension jconf(*)
      dimension jxope(*)
      dimension jxopp(*)
      DIMENSION XMO(*)
      dimension work(*)
      
      call qenter('PAMGRD_MPG')
C
C     Consistency check
C     =================
C
      isym=iprpsym(indxpr)
      iopsy=jbtof(isym-1,1)
      if(iopsy.ne.jopsy)then
         write(lupri,'(a)') 'Incorrect fermion ircop of property ',
     &        prpnam(indxpr)
         call quit('*** ERROR(PAMGRD_MPG) Incorrect fermion ircop')
      endif
C
C     Orbital part
C     ============
C
      nxopt=nxope+nxopp
      if(nxopt.gt.0)then
         call grdorb_mpg(indxpr,gpoe,gpop,jxope,jxopp,nxope,nxopp,
     &                   XMO,work,kfree,lfree,iprint)
      endif

      call qexit('PAMGRD_MPG')
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck grdorb_mpg */
      SUBROUTINE GRDORB_MPG(INDXPR,GPOE,GPOP,JXOPE,JXOPP,NXOPE,
     &                      NXOPP,XMO,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Written by C.V. Larsen
C     Last revision J. van Stralen - June 2002
C
C     PURPOSE: Construct property gradient from the Lagrangian
C
C              e[1]_Q*(ma)=L_ma
C     
C              m - secondary (positronic and virtual)
C              a - occupied
C     
C     INPUT:   GPOE   - Property gradient,V-O part (this is a vector)
C              GPOP   - Property gradient,V-P part (this is a vector)
C              JXOPE  - Orbital rotation index, V-O
C              JXOPP  - Orbital rotation index, P-O
C              NXOPE  - # of occupied - virtual rotations
C              NXOPP  - # of occupied - positronic rotations
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbxpr.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
      dimension gpoe(*)
      dimension gpop(*)
      dimension jxope(*)
      dimension jxopp(*)
      DIMENSION XMO(NORBT,NORBT,NZ)
      dimension work(*)
C
      call qenter('GRDORB_MPG')
C
#if defined (SYS_IBM)
      open(lu1int,status='OLD',form='UNFORMATTED',file='APROPER')
#else
      open(lu1int,status='OLD',form='UNFORMATTED',file='AOPROPER')
#endif
C
      close(lu1int,status='KEEP')
C
C     Construct gradient
C     ==================
C
C     + ---------------+
C     | electronic part|
C     + ---------------+
C
      if(nxope.gt.0)then
         call xopgat(gpoe,jxope,nxope,XMO,norbt,nz)
         call dscal(nz*nxope,-D2,gpoe,1)
      endif
C
C     +-----------------+
C     | positronic part |
C     +-----------------+
C
      if(nxopp.gt.0)then
         call xopgat(gpop,jxopp,nxopp,XMO,norbt,nz)
         call dscal(nz*nxopp,-D2,gpop,1)
      endif
C     ====================
C
C     Printsection
C     ============
C     
c      IF(IPRINT.GE.4) THEN
        CALL HEADER('Output from GRDORB_MPG',-1)
        IREP = IPRPSYM(INDXPR)-1
        IF(NXOPE.GT.0) THEN
          WRITE(LUPRI,'(A,A16)') 
     &    '* Gradient (e-e) of property ',PRPNAM(INDXPR)
          CALL PRQMAT(GPOE,NXOPE,1,NXOPE,1,NZ,IPQTOQ(1,IREP),
     &                LUPRI)
        ENDIF
        IF(NXOPP.GT.0) THEN
          WRITE(LUPRI,'(A,A16)') 
     &    '* Gradient (e-p) of property ',PRPNAM(INDXPR)
          CALL PRQMAT(GPOP,NXOPP,1,NXOPP,1,NZ,IPQTOQ(1,IREP),
     &                LUPRI)
        ENDIF
c      ENDIF
C
      call qexit('GRDORB_MPG')
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck num2str */
      SUBROUTINE NUM2STR(NUM,RESSTR)
C***********************************************************************
C     
C     Written by C.V. Larsen
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      character str(4)*1
      character resstr*4
C
      CALL QENTER('NUM2STR')
C
      ifac=1000
      istr=1
      inum=num
 10   ich=inum/ifac
      str(istr)=char(ich+48)
      istr=istr+1
      inum=mod(inum,ifac)
      ifac=ifac/10
      if(ifac.ge.1)goto 10
C
      resstr=str(1)//str(2)//str(3)//str(4)
c      write(lupri,'(2a)') 'Num2Str: ',resstr
C
      CALL QEXIT('NUM2STR')
C
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getlrf_mpg */
      SUBROUTINE GETLRF_MPG(ABLRF,IBTYP,IBCVC,IBEVC,IBPVC,RCNV,EVALR,
     &                      EVECR,IORBCL,EIG,DMO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by C.V. Larsen
C     Last revision J. van Stralen - June 2002
C
C     PURPOSE: Retrieves the result of the response equation and puts it
C              in the VO and PO block of the 2nd order density matrix
C
C     INPUT:   ABLRF  - linear response function,
C                       not needed for the MP2 gradient
C              IBTYP  - type of the trial vector
C
C                  trial vectors will be of three kinds:
C
C                  C-type: contains only configurational parameters
C                              (WE DON'T NEED THIS ONE NOW!)
C                  E-type: contains only e-e rotations
C                  P-type: contains only e-p rotations
C              IBCVC  - pointer from C-type trial vectors to IBTYP
C              IBEVC  - pointer from E-type trial vectors to IBTYP
C              IBPVC  - pointer from P-type trial vectors to IBTYP
C              RCNV   - total norm of residual relative to norm of
C                       solution vector
C              EVALR  - eigenvalues of the reduced response equation
C              EVECR  - eigenvectors of the reduced response equation
C              IORBCL - something for the linear response function,
C                       not needed for the MP2 gradient
C              EIG    - eigenvalues of the mo's
C              DMO    - the D^2 matrix
C
C     REMARKS: -
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "dcbxrs.h"
#include "dcbxlr.h"
#include "dcbibn.h"
#include "dcborb.h"
#include "dgroup.h"
#include "../moltra/dcbtr3.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
      logical first
      integer a
      dimension ablrf(nfreq,nlrapt,nlrbpt,*)
      dimension ibtyp(2,*)
      dimension ibcvc(*)
      dimension ibevc(*)
      dimension iorbcl(2,*)
      dimension ibpvc(*)
      dimension evalr(*)
      dimension evecr(nredm,*)
      dimension eig(*)
      DIMENSION DMO(NORBT,NORBT,NZ)
      dimension work(*)
C
      CALL QENTER('GETLRF_MPG')
      kfrsav=kfree
C
C     Evaluate solution vectors for all frequencies
C     =============================================
C
      nesim=0
      npsim=0
      ncsim=0
      ifac=1
      if(.not.static) ifac=2
C
C     Get coefficients
C     ================
C
      call memget('REAL',kcmo,ncmotq,work,kfree,lfree)
      iopt=4
      lubuf=22
      call opnfil(lubuf,'DFCOEF','OLD','GETLRF')
      call reacmo(lubuf,'DFCOEF',work(kcmo),dum,dum,toterg,iopt)
      close(lubuf,status='KEEP')
C
C     Solution vector (e-e) part
C     ==========================
C
      if(nzxope.gt.0)then
         write(lupri,'(a)') 'e-e solutions'
c        PRINT*, 'in getlrf_mpg'
         nesim=nfreq*ifac
         kee=kfree
         call memget('REAL',kgvec,      nzxopeq,work,kfree,lfree)
         call memget('REAL',kxvec,nzxopeq*nfreq,work,kfree,lfree)
         call memget('REAL',kbbuf,      nzxopeq,work,kfree,lfree)
         CALL MEMGET('INTE',KIVEC,        NESIM,WORK,KFREE,LFREE)
         call getlr1(ablrf(1,1,1,ioff+1),jbendx,evecr,evalr,RCNV,
     &               ibtyp,ibevc,
     &               work(kgvec),work(kxvec),work(kbbuf),nzxope,
     &               work(kcmo),iorbcl,eig,work(kzxope),WORK(KIVEC),
     &               work,kfree,lfree)
C
C     +-----------------------------+
C     | Result can be isolated here |
C     +-----------------------------+
C
C
         K = 0
         DO IZ=1,NZ
            IO1 = 0
            IO2 = 0
            DO IFRP=1,NFSYM
               IO1 = IO1 + NPSH(IFRP)
               IO2 = IO2 + NPSH(IFRP) + NOCC(IFRP)
               DO I = 1, NOCC(IFRP)
                  II = I + IO1
                  CALL DCOPY(NSSH(IFRP),WORK(KXVEC+K),   1,
     &                                  DMO(1+IO2,II,IZ),1)
                  K = K + NSSH(IFRP)
               ENDDO
               IO1 = IO1 - NPSH(IFRP) + NORB(IFRP)
               IO2 = IO2 - NPSH(IFRP) - NOCC(IFRP) + NORB(IFRP)
            ENDDO
         ENDDO
         CALL MEMREL('GETLRF.e-e',WORK,1,KEE,KFREE,LFREE)
      endif
C
         IO1 = 0
         IO2 = 0
         DO IFRP = 1, NFSYM
            IO1 = IO1 + NPSH(IFRP)
            IO2 = IO2 + NPSH(IFRP) + NOCC(IFRP)
            DO IZ=1,NZ
               DO I=1,NOCC(IFRP)
                  II = I + IO1
                  DO A=1,NSSH(IFRP)
                     IA = A + IO2
                     write(lupri,'(2i3,f18.10)') IA,II,DMO(IA,II,IZ)
                  ENDDO
               ENDDO
            ENDDO
            IO1 = IO1 - NPSH(IFRP) + NORB(IFRP)
            IO2 = IO2 - NPSH(IFRP) - NOCC(IFRP) + NORB(IFRP)
         ENDDO
C
C     Solution vector (e-p) part
C     ==========================
C
      if(nzxopp.gt.0)then
         write(lupri,'(a)') 'e-p solutions'
         npsim=nfreq*ifac
         call memget('REAL',kxvec,nzxoppq*nfreq,work,kfree,lfree)
         call memget('REAL',kgvec,nzxoppq      ,work,kfree,lfree)
         call memget('REAL',kbbuf,nzxoppq      ,work,kfree,lfree)
         CALL MEMGET('INTE',KIVEC,NPSIM        ,WORK,KFREE,LFREE)
         call getlr1(ablrf(1,1,1,ioff+1),jbpndx,evecr,evalr,RCNV,ibtyp,
     &               ibpvc,work(kgvec),work(kxvec),work(kbbuf),nzxopp,
     &               work(kcmo),iorbcl,eig,work(kzxopp),WORK(KIVEC),
     &               work,kfree,lfree)
C
C    +---------------------------+
C    |Result can be isolated here|
C    +---------------------------+
C
         K = 0
         DO IZ=1,NZ
            IO1 = 0
            IO2 = 0
            DO IFRP=1,NFSYM
               IO1 = IO1 + NPSH(IFRP)
               DO I = 1, NOCC(IFRP)
                  II = I + IO1
                  CALL DCOPY(NPSH(IFRP),WORK(KXVEC+K),   1,
     &                                  DMO(1+IO2,II,IZ),1)
                  K = K + NPSH(IFRP)
               ENDDO
               IO1 = IO1 - NPSH(IFRP) + NORB(IFRP)
               IO2 = IO2 + NORB(IFRP)
            ENDDO
         ENDDO
      endif
C
c        IO1 = 0
c        IO2 = 0
c        DO IFRP = 1, NFSYM
c           IO1 = IO1 + NPSH(IFRP)
c           DO IZ=1,NZ
c              DO I=1,NOCC(IFRP)
c                 II = I + IO1
c                 DO A=1,NPSH(IFRP)
c                    IA = A + IO2
c                    write(lupri,'(2i3,f18.10)') IA,II,DMO(IA,II,IZ)
c                 ENDDO
c              ENDDO
c           ENDDO
c           IO1 = IO1 - NPSH(IFRP) + NORB(IFRP)
c           IO2 = IO2 + NORB(IFRP)
c        ENDDO
      CALL MEMREL('GETLRF.end',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL MEMCHK('einde rspmpg',WORK,1)
      CALL QEXIT('GETLRF_MPG')
C
      RETURN
      END
