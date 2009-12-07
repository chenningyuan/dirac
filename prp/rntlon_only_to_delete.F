      subroutine rntlonnew(cmo,ibeig,work,kfree,lfree)
!-----------------------------------------------------------------------
!like xrssep subroutine (pamrvc.F)
#define debug_me
!
#include "implicit.h"
#include "priunit.h"  
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "dummy.h"
#include "dcbfir.h"
#include "dcbham.h"

      PARAMETER (D1=1.0D0,D0=0.0D0)

      DIMENSION CMO(*),IBEIG(*),WORK(*)
      DIMENSION IREP(3), ISYM(3)

      LOGICAL EPMIX,HERM,AH,DOPCTRA_SAVE

!-----------------------------------------------------------------------
      KFRSAV = KFREE
      IPRINT = IPRPRP
      EPMIX = DOEPRN
      DOPCTRA_SAVE=DOPCT
      DOPCT=.FALSE.

! get memory for the final result:
! we allocate 3*N2BBASXQ as we'll reuse it inside called routines
! for the result we need only 3*N2ORBXQ
      call memget('REAL',kh0tb1,3*n2bbasxq,work,kfree,lfree)

! calculate {T^B, g(0)_add}
      call get_gtb(work(kh0tb1),cmo,ibeig,work,kfree,lfree)

! calculate {T^B, F(0)}
! contribution from {T^B, F(0)} is added on top of contribution {T^B, g(0)_add}
      call get_ftb(work(kh0tb1),cmo,ibeig,work,kfree,lfree)

      dopct=dopctra_save
      call memrel('rntlonnew',work,1,kfrsav,kfree,lfree)
      return
      end

!!!!!! do nothing 
      subroutine do_nothing(MAT,WORK,KFREE,LFREE)
        write(*, *) 'alamakota'
      end
!!!!!! do nothing


      subroutine get_gtb(gmo,cmo,ibeig,work,kfree,lfree)
!-----------------------------------------------------------------------
#define debug_me
! subroutine like gmolitx (pamrvc.F)
!we don't need all common blocks, but anyway this is temporary:
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "dummy.h"
#include "dcbfir.h"
#include "dcbham.h"

      PARAMETER (D1=1.0D0,D0=0.0D0)

      DIMENSION CMO(NCMOTQ),IBEIG(NTBAS(0)),WORK(*)
      DIMENSION ISYMOP(3),IHRMOP(3),IFCKOP(3)
      DIMENSION IREP(3),IOPSY(3),ISYM(3), ITIM(3)

      DIMENSION GMO(*)

      LOGICAL EPMIX,HERM,AH,DOPCTRA_SAVE
!-----------------------------------------------------------------------

      kfrsav = kfree

      call memget('REAL',kgao,3*n2bbasxq,work,kfree,lfree)
      call get_gtb_1(gmo,work(kgao),cmo,ibeig,work,kfree,lfree)
      call memrel('get_gtb',work,1,kfrsav,kfree,lfree)
      return
      end


      subroutine get_gtb_1(gmo,gao,cmo,ibeig,work,kfree,lfree)
!-----------------------------------------------------------------------
#define debug_me
!this flag does not work now
!#define use_gmolitx2
! subroutine like gmolitx1
      use dirac_cfg
#ifdef MOD_DFT
      use dft_cfg
      use dft_main
      use dft_temp
#endif /* MOD_DFT */
!we don't need all common blocks, but anyway this is temporary:
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "dummy.h"
#include "dcbfir.h"
#include "dcbham.h"

      PARAMETER (D1=1.0D0,D0=0.0D0)

      DIMENSION CMO(NCMOTQ),IBEIG(NTBAS(0)),WORK(*)
      DIMENSION ISYMOP(3),IHRMOP(3),IFCKOP(3)
      DIMENSION IREP(3),IOPSY(3),ISYM(3), ITIM(3)

      DIMENSION gmo(n2orbxq,*), gao(n2bbasxq,*)

      LOGICAL EPMIX,HERM,AH,DOPCTRA_SAVE

#ifdef MOD_DFT
      type(dft_options_t) :: dft_options
#endif
!-----------------------------------------------------------------------
      KFRSAV = KFREE
      IPRINT = IPRPRP
      EPMIX = DOEPRN

!flags first, todo: use SETGMT (prp/pamrvc.F) for setting these flags?
!now copy'n'pasted from old rntlon
      DO I = 1, 3
        INDXPR    = IPCON(I)
        ITIM      = IPRPTIM(INDXPR)
        ISYM(I)   = IPRPSYM(INDXPR)
        ISYMOP(I) = ISYM(I)
        IHRMOP(I) = -1
        IFCKOP(I) = JBTOF(ISYM(I)-1,1)
        IOPSY(I)  = JBTOF(ISYM(I)-1,1)
        IREP(I)   = ISYM(I) - 1
      ENDDO

!calculate modified densities: (todo: adapt and call dmolitx?)
!we reuse memory and write modified densities to gmo:
      call get_mod_den(gmo,cmo,ibeig,work,kfree,lfree)

!we need npos:
!todo: see STEX_GET_NPOS routine in prp/pamstex.F - maybe call it here?
      call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. ,.TRUE. /))
      call SetIntTaskArrayDimension(NPOS,PARCAL)
      if (NPOS.GT.0) THEN
         CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
      else
         KPOS = KFREE
      endif
    
      nfmtx = 3

      call twofck(isymop,
     &            ihrmop,
     &            ifckop,
     &            gao,           !3 Fock matrices(n2bbasxq each)
     &            gmo,           !3 modified density matrices(n2orbxq each)
     &            nfmtx,         !=3
     &            work(kpos),
     &            intnmr,
     &            iprint,
     &            work(kfree),
     &            lfree) 

!!todo here will come dft contribution (copy'n'pasted, not tested)
!#ifdef MOD_DFT
!        if (dirac_cfg_dft_calculation) then
!
!          call memget('REAL', kbuf, n2bbasxq, work, kfree, lfree)
!!unperturbed density matrix in work(kbuf):
!          call genden(work(kbuf), cmo, 1, iprint)
!
!          call reset_dft_options(dft_options)
!          dft_options%lr = .true.
!
!          call dftexc(work(kbuf),
!     &                gmo,
!     &                gao,
!     &                nfmtx,
!     &                nfmtx,
!     &                isymop,
!     &                ihrmop,
!     &                dft_options,
!     &                iprint)
!
!!          call memrel('get_gtb_1.dft', work, 1, kbuf, kfree, lfree)
!
!        end if
!#endif /* MOD_DFT */

!now in gao we have {T(1), g(0)_add} in AO basis
#ifdef debug_me
        do icomp=1,3
          CALL HEADER(
     &    'RNTLON: {T(1),g(0)_add}_ic in AO basis ',-1)
          WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
          CALL PRQMAT(gao(1,icomp),NTBAS(0),NTBAS(0),
     &          NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
        enddo
#endif /* debug_me */

! do AO-MO transformation
! todo: this ao-mo transformation is done in GMOLITX2 (pamrvc.F) using qbtrans
! (qbtrans is sth that Miro wanted to have - see his comments in old rntlon)

      call dzero(gmo,3*n2bbasxq)

      do icomp=1,3

#ifndef use_gmolitx2

       DO I1 = 1, NFSYM
          I2 = MOD(I1+IOPSY(ICOMP),2) + 1
          CALL QTRANS('AOMO','S',D0,
     &     NFBAS(I1,0),NFBAS(I2,0),
     &     NORB(I1),NORB(I2),
     &     gao(1,icomp),                   !result in AO basis
     &     NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP(ICOMP)),
     &     gmo(1,icomp),                   !space for new result in MO basis
     &     NORBT,NORBT,NZ,IPQTOQ(1,IREP(ICOMP)),
     &     CMO(1+ICMOQ(I1)),
     &     NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &     CMO(1+ICMOQ(I2)),
     &     NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &     WORK(KFREE),LFREE,IPRINT)
       ENDDO

#else /*ifdef use_gmolitx2*/
!gosia: it does not work properly even though it looks really similar
!to 'call qtrans' above
!todo: check what is different
!        call gmolitx2(gmo(1,icomp),
!     &                gao(1,icomp),
!     &                cmo,
!     &                ibeig,
!     &                work(kfree),
!     &                lfree)
#endif /* use_gmolitx2 */

!now in gmo we have {T(1), g(0)_add} in MO basis:
#ifdef debug_me
         CALL HEADER(
     &'RNTLON: Total transf. term, {T(1),g(0)_add}_ic, in MO:',-1)
         WRITE(LUPRI,'(2X,A,I2)') '*** component :',ICOMP
         CALL PRQMAT(gmo(1,icomp),NORBT,NORBT,
     &               NORBT,NORBT,NZ,
     &            IPQTOQ(1,IREP(ICOMP)),LUPRI)
#endif /* debug_me */
      enddo !do icomp = 1, 3

      call memrel('get_gtb_1',work,1,kfrsav,kfree,lfree)
      return
      end


      subroutine get_mod_den(dtbmat,cmo,ibeig,work,kfree,lfree)
!-----------------------------------------------------------------------
!like dmolitx (pamrvc.F)
#define debug_me
!we don't need all common blocks, but anyway this is temporary:
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "dummy.h"
#include "dcbfir.h"
#include "dcbham.h"
!
      PARAMETER (D1=1.0D0,D0=0.0D0)

      DIMENSION CMO(NCMOTQ),IBEIG(NTBAS(0)),WORK(*)
      DIMENSION ISYMOP(3),IHRMOP(3),IFCKOP(3)
      DIMENSION IREP(3),IOPSY(3),ISYM(3), ITIM(3)
      DIMENSION KF1TR(3),KF2TR(3),KMDM(3)
!
      DIMENSION DTBMAT(*)

      LOGICAL EPMIX,HERM,AH,DOPCTRA_SAVE
!-----------------------------------------------------------------------

      KFRSAV = KFREE

      call memget('REAL',ktmat,n2orbxq,work,kfree,lfree)
      call get_mod_den_1(dtbmat,work(ktmat),cmo,ibeig,work,kfree,lfree)
      CALL MEMREL('get_mod_den',WORK,1,KFRSAV,KFREE,LFREE)
      return
      end


      subroutine get_mod_den_1(dtbmat,tmat,cmo,ibeig,work,kfree,lfree)
!-----------------------------------------------------------------------
!like dmolitx1 (pamrvc.F)
!choose one:
#define use_getmdm
!dmoli3 cannot yet be used (comments in pamrvc.F, line 819)
!#define use_dmoli3
#define debug_me
!we don't need all common blocks, but anyway this is temporary:
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "dummy.h"
#include "dcbfir.h"
#include "dcbham.h"
!
      PARAMETER (D1=1.0D0,D0=0.0D0)

      DIMENSION CMO(NCMOTQ),IBEIG(NTBAS(0)),WORK(*)
      DIMENSION ISYMOP(3),IHRMOP(3),IFCKOP(3)
      DIMENSION IREP(3),IOPSY(3),ISYM(3), ITIM(3)
      DIMENSION KF1TR(3),KF2TR(3),KMDM(3)
!
      DIMENSION DTBMAT(N2BBASXQ,*) 
      dimension tmat(n2orbxq) 
! we take n2orbxq for tmat as tmat is used in a do-loop 

      LOGICAL EPMIX,HERM,AH,DOPCTRA_SAVE
!-----------------------------------------------------------------------

      KFRSAV = KFREE
      IPRINT = IPRPRP
      EPMIX = DOEPRN

      call opnfil(lu1int,'AOPROPER','OLD','GET_MOD_DEN')

#ifdef use_getmdm
!inside getmdm we have: dmat(n2bbasx, nz) and tmat(norbt,norbt,nz)
!the same inside dmoli3...

      call memget('REAL',kmdmc,n2bbasxq,work,kfree,lfree)

      DO ICOMP = 1, 3

        INDXPR       = IPCON(ICOMP)
        ITIM(ICOMP)  = IPRPTIM(INDXPR)
        ISYM(ICOMP)  = IPRPSYM(INDXPR)
        IOPSY(ICOMP) = JBTOF(ISYM(ICOMP)-1,1)
        IREP(ICOMP)  = ISYM(ICOMP) - 1 

        call dzero(tmat, n2orbxq)
 
        CALL PRPMAT(INDXPR,IOPSY,TMAT,CMO,
     &              IBEIG,ICMOQ,NORB,WORK,KFREE,LFREE,IPRPRP)

#ifdef debug_me
        CALL HEADER(
     & 'RNTLON: T(1)_ic connection matrix in MO basis (get_mod_den_1)',
     & -1)
        WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
        CALL PRQMAT(TMAT,NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
#endif /* debug_me */

!was in old rntlon
        ITYP = 1
        IH = -1

        call getmdm(epmix,
     &              isym(icomp),
     &              work(kmdmc),
     &              cmo,
     &              ibeig,
     &              tmat,
     &              work,
     &              kfree,
     &              lfree,
     &              iprprp)

#ifdef debug_me
         CALL HEADER(
     &'RNTLON: Total SYMMETRIZED modified density matrix',-1)
        WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
         CALL PRQMAT(WORK(KMDMC),NTBAS(0),NTBAS(0),
     &     NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
         CALL FLSHFO(LUPRI)
#endif /* ifdef debug_me */

        call dcopy(n2bbasxq,work(kmdmc),1,dtbmat(1,icomp),1)
        call dzero(work(kmdmc),n2bbasxq)

      END DO
#endif /* use_getmdm */

#ifdef use_dmoli3

      call memget('REAL',kmdmc,n2bbasxq,work,kfree,lfree)
      ITYP = 1
      IH = -1
      do icomp = 1, 3
   
        INDXPR       = IPCON(ICOMP)
        ITIM(ICOMP)  = IPRPTIM(INDXPR)
        ISYM(ICOMP)  = IPRPSYM(INDXPR)
        IOPSY(ICOMP) = JBTOF(ISYM(ICOMP)-1,1)
        IREP(ICOMP)  = ISYM(ICOMP) - 1 

        CALL DZERO(TMAT, N2ORBXQ)
 
        CALL PRPMAT(INDXPR,IOPSY,TMAT,CMO,
     &              IBEIG,ICMOQ,NORB,WORK,KFREE,LFREE,IPRPRP)

#ifdef debug_me
        CALL HEADER(
     & 'RNTLON: T(1)_ic connection matrix in MO basis (call dmoli3)',
     & -1)
        WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
        CALL PRQMAT(TMAT,NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
#endif /* debug_me */

        call dmoli3(ityp,ihrmop(icomp),work(kmdmc),cmo,tmat,
     &              work,kfree,lfree)

        call dcopy(n2bbasxq,work(kmdmc),1,dtbmat(1,icomp),1)
        call dzero(work(kmdmc),n2bbasxq)
      enddo
#endif /* use_dmoli3 */
!     
      CLOSE(LU1INT, STATUS='KEEP')

      CALL MEMREL('GET_MOD_DEN',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END 


      subroutine get_ftb(ftb,cmo,ibeig,work,kfree,lfree)
!-----------------------------------------------------------------------
!ftb contains G contrib.
#define debug_me
!
!we don't need all common blocks, but anyway this is temporary:
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "dummy.h"
#include "dcbfir.h"
#include "dcbham.h"
!
! in ttra we have: tmat(norbt,norbt,nz), hmat and thmat - the same
! in fmoli3: fmo, fbuf, bmat also: (norbt, norbt, nz)
!
      PARAMETER (D1=1.0D0,D0=0.0D0)

      DIMENSION CMO(NCMOTQ),IBEIG(NTBAS(0)),WORK(*)
      DIMENSION ISYMOP(3),IHRMOP(3),IFCKOP(3)
      DIMENSION IREP(3),IOPSY(3),ISYM(3), ITIM(3)

      DIMENSION ftb(*)

      LOGICAL EPMIX,HERM,AH,DOPCTRA_SAVE
!-----------------------------------------------------------------------
      KFRSAV = KFREE
      IPRINT = IPRPRP
      EPMIX = DOEPRN

      call memget('REAL',ktmat,n2bbasxq,work,kfree,lfree)
      call memget('REAL',kfmo,n2bbasxq,work,kfree,lfree)

!get unperturbed Fock in MO basis:
!todo: it would be less expensive to read Fock from files
      call getfck(work(kfmo),iprprp,work,kfree,lfree) 

#ifdef debug_me
        CALL HEADER(
     &'RNTLON: Total Fock matrix in MO basis:',-1)
        CALL PRQMAT(WORK(KFMO),NORBT,NORBT,NORBT,
     &            NORBT,NZ,IPQTOQ(1,0),LUPRI)
#endif /* debug_me */

      call get_ftb_1(ftb,work(kfmo),work(ktmat),
     &               cmo,ibeig,work,kfree,lfree)

      call memrel('get_ftb',work,1,kfrsav,kfree,lfree)
      return
      end


      subroutine get_ftb_1(ftb,fmat,tmat,cmo,ibeig,work,kfree,lfree)
!-----------------------------------------------------------------------
!like fmoli2
!#define use_ttra
#define use_fmoli3
#define debug_me
!
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "dummy.h"
#include "dcbfir.h"
#include "dcbham.h"
!
! in ttra we have: tmat(norbt,norbt,nz), hmat and thmat - the same
! in fmoli3: fmo, fbuf, bmat also: (norbt, norbt, nz)
!
      PARAMETER (D1=1.0D0,D0=0.0D0)

      DIMENSION CMO(NCMOTQ),IBEIG(NTBAS(0)),WORK(*)
      DIMENSION ISYMOP(3),IHRMOP(3),IFCKOP(3)
      DIMENSION IREP(3),IOPSY(3),ISYM(3), ITIM(3)
      DIMENSION KF1TR(3),KF2TR(3),KMDM(3)
!unperturbed Fock in FMAT, result in FTB
      DIMENSION FMAT(*)  
      DIMENSION FTB(N2ORBXQ,*)
      dimension tmat(*)

      LOGICAL EPMIX,HERM,AH,DOPCTRA_SAVE
!-----------------------------------------------------------------------
      KFRSAV = KFREE
      IPRINT = IPRPRP
      EPMIX = DOEPRN

#ifdef use_ttra
! flags for ttra:
    !old:
    !ityp = 1       ! e-e only
    !herm = .false. !antisymmetric
    !isymh = 1      ! do {T,F} and not [T,F]
    !ah    = .true. !then take -1 for hermitian conjugate, else take +1

      OPEN(LUFCKL,FILE='FCKLON',FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=8*N2ORBXQ,STATUS='UNKNOWN')

      CALL OPNFIL(LU1INT,'AOPROPER','OLD','RNTLON')

      DOPCTRA_SAVE=DOPCT
      DOPCT=.FALSE.

      call memget('REAL',ktemp,n2orbxq,work,kfree,lfree)

      DO ICOMP = 1, 3
       INDXPR          =  IPCON(ICOMP)
       ITIM            =  IPRPTIM(INDXPR)
       ISYM(ICOMP)     =  IPRPSYM(INDXPR)
       IOPSY(ICOMP)    =  JBTOF(ISYM(ICOMP)-1,1)
       IREP(ICOMP)     =  ISYM(ICOMP) - 1

       ITYP = 1
       HERM  = .FALSE.

        CALL PRPMAT(INDXPR,IOPSY,tmat,CMO,
     &              IBEIG,ICMOQ,NORB,WORK,KFREE,LFREE,IPRPRP)

#ifdef debug_me
        CALL HEADER(
     &  'RNTLON: T(1)_ic connection matrix in MO basis (get_ftb_1)',-1)
        WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
        CALL PRQMAT(tmat,NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
#endif /* debug_me */

       ISYMH = 1
       AH = .TRUE.
       call dzero(work(ktemp),n2orbxq)
       CALL TTRA(EPMIX,AH,ISYM(ICOMP),ISYMH,tmat,
     &           FMAT,work(ktemp),IPRPRP)

!now in work(ktemp) we have {T(1), F(0)} in MO basis:
#ifdef debug_me
        CALL HEADER(
     &  'RNTLON: Transformed {T(1),F(0)}_icomp in MO basis ',-1)
        WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
        CALL PRQMAT(WORK(ktemp),NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
#endif /* debug_me */

       CALL DAXPY(N2ORBXQ,D1,WORK(ktemp),1,FTB(1,ICOMP),1)
!ftb contains all now 
#ifdef debug_me
        CALL HEADER(
     &'RNTLON: Total transformed term, {T(1),H(0)}_comp, in MO:',-1)
        WRITE(LUPRI,'(2X,A,I1)') '*** Component x-1/y-2/z-3:',ICOMP
        CALL PRQMAT(ftb(1,ICOMP),
     &              NORBT,NORBT,NORBT,NORBT,NZ,
     &            IPQTOQ(1,IREP(ICOMP)),LUPRI)
#endif /* debug_me */
!so write to FCKLON:
       CALL WRTDAC(LUFCKL,N2ORBXQ,ftb(1,icomp),ICOMP)

#ifdef debug_me
        CALL HEADER(
     &  'RNTLON: Transformed {T(1),F(0)}_icomp in MO basis ',-1)
        WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
        CALL PRQMAT(WORK(ktemp),NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
#endif /* debug_me */

      ENDDO

      close(lufckl,status='keep')

#endif /* use_ttra */

#ifdef use_fmoli3
! flags for fmoli3:

      OPEN(LUFCKL,FILE='FCKLON',FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=8*N2ORBXQ,STATUS='UNKNOWN')

      CALL OPNFIL(LU1INT,'AOPROPER','OLD','RNTLON')

      DOPCTRA_SAVE=DOPCT
      DOPCT=.FALSE.

      call memget('REAL',ktemp,n2orbxq,work,kfree,lfree)

      DO ICOMP = 1, 3
       INDXPR          =  IPCON(ICOMP)
       ITIM            =  IPRPTIM(INDXPR)
       ISYM(ICOMP)     =  IPRPSYM(INDXPR)
       IOPSY(ICOMP)    =  JBTOF(ISYM(ICOMP)-1,1)
       IREP(ICOMP)     =  ISYM(ICOMP) - 1

       ITYP = 1
       HERM  = .FALSE.

        CALL PRPMAT(INDXPR,IOPSY,tmat,CMO,
     &              IBEIG,ICMOQ,NORB,WORK,KFREE,LFREE,IPRPRP)

#ifdef debug_me
        CALL HEADER(
     &  'RNTLON: T(1)_ic connection matrix in MO basis (get_ftb_1)',-1)
        WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
        CALL PRQMAT(tmat,NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
#endif /* debug_me */

       ISYMH = 1
       AH = .TRUE.
       call dzero(work(ktemp),n2orbxq)

       call fmoli3(ityp,work(ktemp),tmat,fmat,work,kfree,lfree)

!now in work(ktemp) we have {T(1), F(0)} in MO basis:
#ifdef debug_me
        CALL HEADER(
     &  'RNTLON: Transformed {T(1),F(0)}_icomp in MO basis ',-1)
        WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
        CALL PRQMAT(WORK(ktemp),NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
#endif /* debug_me */

       CALL DAXPY(N2ORBXQ,D1,WORK(ktemp),1,FTB(1,ICOMP),1)
!ftb contains all now 
#ifdef debug_me
        CALL HEADER(
     &'RNTLON: Total transformed term, {T(1),H(0)}_comp, in MO:',-1)
        WRITE(LUPRI,'(2X,A,I1)') '*** Component x-1/y-2/z-3:',ICOMP
        CALL PRQMAT(ftb(1,ICOMP),
     &              NORBT,NORBT,NORBT,NORBT,NZ,
     &            IPQTOQ(1,IREP(ICOMP)),LUPRI)
#endif /* debug_me */
!so write to FCKLON:
       CALL WRTDAC(LUFCKL,N2ORBXQ,ftb(1,icomp),ICOMP)

#ifdef debug_me
        CALL HEADER(
     &  'RNTLON: Transformed {T(1),F(0)}_icomp in MO basis ',-1)
        WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
        CALL PRQMAT(WORK(ktemp),NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
#endif /* debug_me */

!       CALL DAXPY(N2ORBXQ,D1,WORK(KBUF),1,FTB(1,ICOMP),1)

      ENDDO

      close(lufckl,status='keep')



#endif /* use_fmoli3 */

      CLOSE(LU1INT,STATUS='KEEP')
      CALL MEMREL('DO_FTB',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END

