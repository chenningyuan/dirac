      SUBROUTINE IIII(IM,INN)
      implicit real*8 (A-H,O-Z)
      DIMENSION IM(INN)
      CALL ISETVC(IM,INN,INN)
      DO I=1,INN
        PRINT*,'IM(I),I',IM(I),I
      END DO
      END
      SUBROUTINE INTIM_KRCC(INTSPACE,MCSCF,
     &                   WORK,KFREE,LFREE)
*
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
      integer POI(30)
      logical MCSCF
      common /WORKPIC/IPOI
*
#include "mxpdim.inc"
#include "symm.inc"
#include "crun.inc"
*
      DIMENSION WORK(*)
C
C
C#include "memint.h"
*
Commmented out by Lasse since i do not think i need it on krcc
C     CALL INTIM_KRCC_POINTERS
* for global purposes
      INTITER = 1
C Stupid transfer of information
      CALL SETUP_FOR_itrctl_KRCC
*  Determine double group information
      CALL SSPIR_KRCC
*
*  Move to integral transformation section:
*  Define type of transformation
*  Here for now: complete two-component trans for CI
*
C Assuming this position (at least until a new one is taken)
      ITRA_TYPE = 2
*
      call itrctl_KRCC(INTSPACE,ITRA_TYPE,INTITER,
     &            ISPINFREE,NORB_F,ISYM_T,
     &            F1,F2,IBOSYM,MCSCF,CIRUN,
     &            WORK,KFREE,LFREE)
*
      return
      end
*
*  Control section for integral transformation
*
      subroutine itrctl_KRCC(INTSPACE,ITRA_TYPE,INTITER,
     &                  ISPINFREE,NORB_FULL,
     &                  ISYM_T,
     &                  F1,F2,IBOSYM,MCSCF,CIRUN,
     &                  WORK,KFREE,LFREE)
*
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "symm.inc"
#include "glbbas.inc"
#include "clunit.inc"
#include "ctcc.inc"
#include "ctccp.inc"
#include "cgas.inc"
#include "dgroup.h"
*
      character*1 COORD
      logical MCSCF
      dimension INTEG(1000)
      dimension NORB_FULL(MXPIRR)
*
* ITRA_TYPE       Defines the type of transformation
*---------------------------------------------------
* ITRA_TYPE = 1 : Normal nonrelativistic transformation
* ITRA_TYPE = 2 : 2/4-component relativistic transformation for CI,
*                 all nonredundant integral lists
* ITRA_TYPE = 3 : 2/4-component relativistic transformation for CI,
*                 Excluding integrals with an odd number of bars
* ITRA_TYPE = 4 : 2/4-component relativistic transformation for MCSCF,
*                 or density matrix spinor optimization
*                 all nonredundant integral lists
* ITRA_TYPE = 5 : 2/4-component relativistic transformation for MCSCF,
*                 or density matrix spinor optimization
*                 Excluding integrals with an odd number of bars
*
*
C Will always assume dirac enviroment, sigden route and transformation of 4 or higher 
C                                                                 Lasse
C
*
      DIMENSION WORK(*)
C
C
C#include "memint.h"
*
      DIRACONE = 'MRCONEE'
      DIRACTWO = 'MDCINT'
C assuming INTIMP = 7 Lasse
      INTIMP = 7
      if (MCSCF) INTIMP = 7
*
      call qenter('ITRCTRL')
*
      NTESTL = 00
      NTEST = 0
*
      if (NTEST.ge.1) then
        WRITE(6,*) ' *************************************************'
        WRITE(6,*) ' *                                               *'
        WRITE(6,*) ' *  Transformer man                              *'
        WRITE(6,*) ' *  You run the show                             *'
        WRITE(6,*) ' *  Remote control                               *'
        WRITE(6,*) ' *  Direct the action with the push of a button  *'
        WRITE(6,*) ' *                                               *'
        WRITE(6,*) ' *                               Neil Young      *'
        WRITE(6,*) ' *                                               *'
        WRITE(6,*) ' *************************************************'
      end if
*
*  CI or MCSCF transformation?
*  Note: We will do the one-particle transformation always with
*        the full orbital space, it doesn't kill us. Afterwards,
*        a routine will select a matrix of integrals over the
*        active spaces for the GAS routines.
*
      NTRIAEL = ((SYMDIM(1)+SYMDIM(2))**2
     &          +(SYMDIM(1)+SYMDIM(2))) /2
*
      NTRSM_FULL = 0
      NDIM_TOT = 0
      do ISYM=1,NIRR_PN,1
        NDIM_TOT = NDIM_TOT + SYDI(ISYM)
      end do
*
      call isetvc(SYMDIM_F,0,NIRR_DG)
      if (DOUGRP.lt.7.and.DOUGRP.ge.4) then
        I = 1
        do IRR=1,NIRR_PN,1
          if (IRR.eq.3.or.IRR.eq.5.or.IRR.eq.7) I = I + 1
          SYMDIM_F(I) = SYMDIM_F(I) + NORB_FULL(IRR)
        end do
      else if (DOUGRP.eq.7) then
        SYMDIM_F(1) = 2*NORB_FULL(1)
        SYMDIM_F(2) = 2*NORB_FULL(2)
      else if (DOUGRP.eq.8) then
        SYMDIM_F(1) = 2*NORB_FULL(1)
      end if
*
*
*  1) Static memory
*
*  Transformed integrals, 32 lists
      print*,'Static memory'
      print*,'MXPIRR,MXPNGAS',MXPIRR,MXPNGAS
      print*,'NIRR_DG,NIRR_PN',NIRR_DG,NIRR_PN
      call findhigh(NGSSH,MXPIRR*MXPNGAS,MXGSDIM)
      if (NTEST.ge.1) write(6,*) 'Highest GAS dimension is  ',MXGSDIM
      call findhigh(SYMDIM(1),NIRR_DG,IMXDDIM)
      if (NTEST.ge.1) write(6,*) 'Sohing else high is IMXDDIM ',IMXDDIM
C NORB_FULL is for externalle generated FOCK matrix and not needed here. Can be deleted along with SYMDIM_F
      call findhigh(NORB_FULL(1),NIRR_PN,IMXDIM_F)
      if (NTEST.ge.1) write(6,*) 'Thing else high is IMXDIM_F ',IMXDIM_F
      call findhigh(SYDI(1),NIRR_DG,IMXDIM)
      if (NTEST.ge.1) write(6,*) 'Something else high is IMXDIM ',IMXDIM
      INTSPACE = (2*MXGSDIM)**2*NGAS**2*32
      ITRALINE = (2*MXGSDIM)**2*NGAS**2
*
*
      CALL MEMGET('INTE',KBOSYM,2*NDIM_TOT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KISREO,2*NDIM_TOT,WORK,KFREE,LFREE)
*
* Allocate for required integrals/densities only
C Should depend on real or complex groups! 
      IRIMLT = MIN(NZ,2)
* !! Still enforced
      NALLINT = LEN_T_VEC * IRIMLT
      print*,'IRIMLT,NZ,NALLINT,LEN_T_VEC',IRIMLT,NZ,NALLINT,LEN_T_VEC
      write(6,*) ' Allocating for ',NALLINT/IRIMLT,' integrals.'
      write(6,*) ' Allocated amount of memory ',NALLINT
      IF(IRIMLT.EQ.1) THEN
        write(6,*) ' only real integrals '
      ELSE
        write(6,*) ' real and imaginary integrals '
      END IF
      CALL MEMGET('REAL',KT_CC,NALLINT,WORK,KFREE,LFREE)
*
*
************************************************
* Set mark for flushing local memory
************************************************
c     call memmar(KDUM,IDUM,'MARK  ',IDUM,'ITRCTL')
*
*  2) Local memory
*  Some index and GAS info
*  Scratch space, HELPMAT2
        NUM_ONEI = (2*NDIM_TOT)**2
        NUM_TWEI = (2*NDIM_TOT)**4
      CALL MEMGET('INTE',KKR,(2*NDIM_TOT+1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINVREO,2*NDIM_TOT+1,WORK,KFREE,LFREE)
*
*  for C2H reordering:
      IHELEN = 0
      do I=2,3,1
         IHELEN = IHELEN + (ITOTDM(I)**2 + ITOTDM(I)) / 2
      end do
*
      write(6,*) 'IOPER,ITRA_TYPE',
     &            IOPER,ITRA_TYPE
*
*******************************************************
*  Reading integrals from DIRAC
*******************************************************
*
C     print*,'check before READ_DIRAC_FOR_KRCC'
C     CALL MEMCHK_KRCC(WORK)
C     CALL READ_DIRAC_FOR_KRCC(WORK,KFREE,LFREE) !to start with
C     print*,'check after READ_DIRAC_FOR_KRCC'
C     CALL MEMCHK_KRCC(WORK)
*
*******************************************************
*  Integral reading section
*******************************************************
*
*
*******************************************************
*  Distribution and transformation section
*******************************************************
*
*
*
          if (NTEST.ge.1) then
            write(6,*)
            write(6,*) '********************************************'
            write(6,*) '* Distributing MOLTRA transformed          *'
            write(6,*) '* dirac one- (and two-) electron integrals *'
            write(6,*) '* using input integral reordering          *'
            write(6,*) '********************************************'
          end if
*
          IRIOFF = (2*NDIM_TOT)**2
          NDIM_TOT2 = NDIM_TOT*2
C     CALL MEMCHK_KRCC(WORK)
C     CALL STUPID_ROUTINE(WORK(KLSOBEX),NSPOBEX_TP,NGAS)
      CALL MEMGET('INTE',KISTR_CA,MX_ST_TSOSO_BLK_MX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KISTR_CB,MX_ST_TSOSO_BLK_MX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KISTR_AA,MX_ST_TSOSO_BLK_MX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KISTR_AB,MX_ST_TSOSO_BLK_MX,WORK,KFREE,LFREE)
*
C     CALL MEMCHK_KRCC(WORK)
       call dist_ints_sigden_KRCC(INTIMP,ISPINFREE,ISYM_T,
     &                            WORK(KLSOBEX),WORK(KT_CC),
     &                            WORK(KISTR_CA),WORK(KISTR_CB),
     &                            WORK(KISTR_AA),WORK(KISTR_AB),
     &                            NDIM_TOT2,WORK(KFCACM),WORK(KH2AC),
     &                            WORK(KBOSYM),
     &                            WORK(KLIBSOBEX),WORK(KINC_HX),
     &                            WORK(KLABEXTP),NTEST,
     &                            WORK,KFREE,LFREE)
*
        write(6,*)
        write(6,*) ' Number of integrals: ',NALLINT
        write(6,*) ' (including imaginary parts) '
C       stop 'after dist_ints_sigden_KRCC'
*
*  Flush local memory to mark set at beginning of do loop
C     call memmar(KDUM,IDUM,'FLUSM ',IDUM,'ITRCTL')
      CALL MEMREL('ITRCTRL',WORK,1,KKR,KFREE,LFREE)
*
      call qexit('ITRCTRL')
      return
      end
*
******************************************************************
*                                                                *
*  Another silly routine by Lasse                                *
*                                                                *
******************************************************************
*
      SUBROUTINE SUM_ARRAY
*
*     Array ISUMORB will be used in GET_INT2
*
#include "mxpdim.inc"
#include "cgas.inc"
*
      IMAXORB = 0
      DO IGAS = 1, NGAS
        IORB = NGSOBT(IGAS)
        IF(IORB.GT.IMAXORB) THEN
          IMAXORB = IORB
        END IF
      END DO
*
*     Now the sum
*
      ISUMORB(1) = 0
      ISUM = 0
      DO I = 1, IMAXORB-1
       ISUM = ISUM + I
       ISUMORB(I+1) = ISUM
      END DO
*
      RETURN
      END
*
*  Import one-electron integrals from DIRAC environment
*
      subroutine rdone_dirac_KRCC(ONEINT,LUONE,DIRACONE,MORB,NDIM,
     &                           ISPINFREE,IBOSYM,IDUMAR,DUMAR,IPRNT)
*
      implicit real*8 (A-H,O-Z)
*
*
      character*14 CHA(8)
      common/CECORER/ECORE
*
      character*7 DIRACONE
      logical BREIT
*
      dimension IDUMAR(*),IBOSYM(*),NSTR(2)
      dimension ONEINT((2*NDIM)**2*2),DUMAR(*)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRNT)
*
      open(LUONE,FILE=DIRACONE,FORM='UNFORMATTED')
      REWIND(LUONE)
      read(LUONE) MORB,BREIT,ECORE_DIR,NFSYM,NZ
      read(LUONE) NSYMRP,(CHA(I),I=1,NSYMRP),(IDUMAR(I),I=1,NSYMRP),
     &                                         (NSTR(I),I=1,NFSYM)
      read(LUONE)
      read(LUONE)
      read(LUONE) (IDUMAR(ID),IDUMAR(ID),DUMAR(ID),ID=1,MORB,1),
     &            (IBOSYM(I),I=1,MORB,1)
      if (NTEST.ge.3) then
        write(6,*)
        write(6,*) 'Testing read MRCONEE'
        write(6,*) 'MORB = ',MORB
        write(6,*) 'BREIT = ',BREIT
        write(6,*) 'ECORE_DIR = ',ECORE_DIR
        write(6,*) 'Number of spinors per fermion symmetry:'
        write(6,'(2I5)') (NSTR(I),I=1,NFSYM)
        write(6,*) 'Read boson symmetries of spinors:'
        do IB=1,MORB,1
          write(6,'(A11,I5,A3,I12)') 'Spinor no. ',IB,' : ',IBOSYM(IB)
        end do
      end if
*
      if (2*NDIM.ne.MORB) then
        write(6,*) 'Orbital dimensions of '
        write(6,*) 'DIRAC and input do not match.'
        write(6,*) 'DIRAC : ',MORB,' spinors'
        write(6,*) 'Input : ',NDIM,' orbitals'
        call abend2('Quitting.')
      end if
*
* Read the integrals
      NDIMSQ = MORB*MORB
      if (ISPINFREE.eq.1) then
        read(LUONE) (ONEINT(I),TDUM,I=1,MORB*MORB)
      else if (ISPINFREE.eq.0) then
        read(LUONE) (ONEINT(I),ONEINT(I+NDIMSQ),I=1,MORB*MORB)
*
* Set IBOSYM array so that U->alpha and B->beta spinors
*              To get correct Fermion symmetry labeling in
*              LUCIAREL we need to specify same(!) gerade boson
*              irrep for all gerade orbitals and same ungerade
*              boson irrep for all ungerade orbitals (note: which
*              gerade or ungerade boson irrep we choose doesn't
*              matter).
        IOFFFS = 0
        IRROFF = 0
        do IFSYM = 1,NFSYM,1
          do ISP = 1+IOFFFS,NSTR(IFSYM)+IOFFFS,1
            IBOSYM(ISP) = 0 + IRROFF
          end do
*
          if (NTEST.ge.3) then
            write(6,*) ' Reset IBOSYM array, ferm symm. :',IFSYM
            write(6,*) (IBOSYM(I),I=1+IOFFFS,NSTR(IFSYM)+IOFFFS,1)
          end if
*
          IOFFFS = NSTR(IFSYM)
          IRROFF = 2
          if (NZ.eq.1) IRROFF = 1
        end do
*
      end if
*
      if (NTEST.ge.5) then
        write(6,*)
        MORB_L = MORB
        write(6,*) 'Real One-electron ints from dirac:'
        call wrtmat(ONEINT(1),MORB_L,MORB_L,MORB_L,MORB_L)
        if (ISPINFREE.eq.0) then
          write(6,*) 'Imaginary One-electron ints from dirac:'
          call wrtmat(ONEINT(1+NDIMSQ),MORB_L,MORB_L,MORB_L,MORB_L)
        end if
      end if
*
      ECORE = ECORE_DIR
      close(LUONE)
*
      return
      end
*
      SUBROUTINE READ_DIRAC_FOR_KRCC(WORK,KFREE,LFREE)
*
* Chaotic routine written by lasse to get the needed integrals 
* for the integral sorting step. Probably a lot of redundancies 
* since taken without any knowledge from som routines in krmc.
* Don't think just do
* Much is a copy paste from KRCI_CALC and KRCI_CALC1
*
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "thrzer.h"
C
#include "ipoist8.inc"
#include "maxorb.h"
!hope we don't need it
#include "../krmc/dcbnrt.h"
#include "dcbpsi.h"
#include "dcbkrmc_itinfo.h"
#include "dcbgen.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "frame.h"
#include "dgroup.h"
C#include "inforb.h"
#include "krciprop.h"
#include "mxpdim.inc"
#include "ctcc.inc"
#include "ctccp.inc"
#include "cgas.inc"
#include "glbbas.inc"
C
      LOGICAL   WFCONV, DO4INDEX, TOBE, EX
      CHARACTER SZTRA(2,4)*80
C     DIMENSION CMO(*), CREF(*)
C     CHARACTER WFTYP*(*), SECTID*12, CPUTID*12, WALLTID*12
C
*
      DIMENSION WORK(*)
C
C
      NTEST = 0
C#include "memint.h"
*
      CALL QENTER('KRCC_INT')
C
C From KRCI_CALC
C
      CALL OPNFIL(LUKRMC,'KRMCSCF','UNKNOWN','PSIOPT')
      CALL NEWLAB('*KRMCSCF',LUKRMC,LUPRI)
      CALL SETDCBIDX()
      CALL RSETWOP('KR-CI',WORK,KFREE,LFREE)
      CALL FLSHFO(LUPRI)

C
C From KRCI_CALC_1
C
C
C
C     ***************************************************
C     *** Initialize some parameters ********************
C     ***************************************************
C
      DO I = 1,2
         IOPT_ISYMOP(I) = 1
         IOPT_IFCKOP(I) = 1
         IOPT_IHRMOP(I) = 1
      END DO
      JTRLVL = 0
      ITMAC  = 0
      DO4INDEX = ( MCTYPE .GE. JMCMIN )
C     ... MCTYPE .lt. JMCMIN is SCF
      EMCOLD = D0
      WFCONV = .FALSE.
C
      CALL DZERO(GNORM,5)
C
C     Integrals for KR-CC calculation. Which are actually needed
C
      print*,'AT THE MOMENT HARDWIRED TO LL LS AND SS INTEGRALS'
      IOPT_INTDEF = 7!0111
      IOPT_INTBUF = 7!0111
      CALL INTCON(IOPT_INTFLG,IOPT_INTBUF,IOPT_INTDEF,GNORM(5),
     &            OPT_CNVINT,ITMAC,IOPT_ITRINT,OPT_INTTYP)
*     Mo's from DFCOEF
      EX = .FALSE.
      print*,'NCMOTQ',NCMOTQ
      CALL MEMGET('REAL',   KCMO  ,NCMOTQ  ,WORK,KFREE,LFREE)
         INQUIRE ( FILE = 'DFCOEF', EXIST = EX )
         IF ( EX ) THEN
            CALL OPNFIL(LUCOEF,'DFCOEF','OLD','RREADMO')
            CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUMMY,DUMMY,DUMMY,2,
     &                  WORK,KFREE,LFREE)
            CLOSE(LUCOEF,STATUS='KEEP')
            IF (NTEST.GE. 0) THEN
               WRITE(LUPRI,'(/A)')
     &              ' (ROPTST)  Start orbitals read from' //
     &              ' file DFCOEF'
            END IF
         END IF
      CALL RGMCMO(WORK(KCMO),WORK(KFREE),LFREE)
* Stupid routine to see if the 4-index transformation should be done
      CALL AN_IDIOT(NO4IDX)
* Transformation 
      CALL INISZT(JTRLVL,SZTRA,LNOROT,NOROT)
C      CALL RTRACTL(CMO,SZTRA,WORK(KFREE),LFREE)
      CALL MEMGET('INTE',KIBEIG,NORBT,WORK,KFREE,LFREE)
      CALL IZERO(WORK(KIBEIG),NORBT)
C
      IF(NO4IDX.EQ.0) THEN
      CALL RTRACTL(WORK(KCMO),WORK(KIBEIG),SZTRA,WORK(KFREE),LFREE)
C      CALL RCIST(CMO,CREF,WORK,KFREE,LFREE)
C no need for these calls. Lasse 2012
C     CALL MEMGET('REAL',KCREF ,0       ,WORK,KFREE,LFREE)
C     CALL RCIST(WORK(KCMO),WORK(KCREF),WORK,KFREE,LFREE)
      CALL MEMCHK('cemory2',WORK,KCMO)
      END IF
*
*
* Calculate NASHT
*
c     NASHT = 0
c     print*,'NGAS,NFSYM',NGAS,NFSYM
c     DO I=1,NGAS
c      DO J =1,NFSYM
c       NASHT = NASHT + NGSSH(J,I)
c       print*,'NGSSH(J,I)',NGSSH(J,I)
c      END DO
c     END DO
*
      N_FCACM = (2*NASHT)*(2*NASHT)
C
C       Calculate core Hamiltonian:
C       ---------------------------
C
C       The core Hamiltonian is the active-active part of FC,
C       a.k.a. FCAC.
C
      MZ = MIN(NZ,2)
c     N2ASHXQ = NASHT*NASHT* NZ
C     print*,'MZ,N_FCACM,N2ASHXQ,N2ORBXQ',MZ,N_FCACM,N2ASHXQ,N2ORBXQ
      CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFCAC,N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFC,N2ORBXQ,WORK,KFREE,LFREE)
      CALL rGETFC(CWRTFO_MAT,WORK(KCMO),WORK(KFC),ECORE,
     &            WORK(KFREE),LFREE)
*
      IPROPT = 000
      CALL RGETAC(WORK(KFC),WORK(KFCAC),IPROPT)
*
c     NASH(1) = NASHT
C     NASH(3) =9
C     NASH(6) = 9
C     NASH(7) = 9
      CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPROPT)
      CALL MEMREL('after rGETFC',WORK,1,KFCAC,KFREE,LFREE)
      LEN_ALL_INT = 0
c     NNASHX = NASHT * (NASHT+1) / 2
      LEN_ALL_INT = NASHT*NASHT*NNASHX*NZ*3
      CALL MEMGET('REAL', KH2AC, LEN_ALL_INT, WORK, KFREE, LFREE)

C     orbital symmetry info array
      call memget('INTE',kibeig, norbt,work,kfree,lfree)
      call izero(work(kibeig),norbt)

!     read integrals into WORK(KH2AC)
      call RGETH2(dummy,work(kh2ac),dummy,work(kibeig),.false.,
     &            .true.,.true.,work(kfree),lfree)
C      CALL RGETH21_KRCC(.FALSE.,.TRUE.,DUMMY,DUMMY,WORK(KH2AC),
C    &               WORK(KFREE),LFREE,IPROPT)
*
      CALL QEXIT('KRCC_INT')
      print*,'exit KRCC_INT'
*
      RETURN
      END
*
      SUBROUTINE AN_IDIOT(I)
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "crun.inc"
*
      I = NO4IDX
*
      RETURN
      END
*
* Driver routine for distribution of DIRAC integrals to
* double group GAS lists
*   Routine using excitation operator sequence for integral
*   storage definition (SIGDEN route)
*
*----------------------------------------------
*  Timo Fleig, Mar. 2001 (Ellen Tabea is here!)
*  Jeppe Olsen modifies it. (Apr. 2001)
*  In complete form, June 2001, Timo Fleig
*
*  Implemented real groups (DIRAC) treated like complex groups
*  with zero imaginary parts.
*  Timo Fleig, August 25, 2003
*
*  Added fetching of integrals for inactive Fock-matrix 05 Lasse
*  Removed them again 08 since KT_CC is now used instead. Lasse
* 
*  Added array IOFFINTTYPE for integral fetching in effecient CC code. Lasse
*
* Demolished and revised for KRCC Lasse 2010
* Now all integrals are coulomb - exchange
*
* Made an ugly fix when reading in the integrals from the ADUCCA code
* since there are changes in the signs from the old to the new code
* which i do not know how to handle - Lasse 2013
*
*----------------------------------------------
*
      subroutine dist_ints_sigden_KRCC(INTIMP,ISPINFREE,ISYM_T,ISOBEX,T,
     &                            ISTR_CA,ISTR_CB,ISTR_AA,ISTR_AB,
     &                            NASHT,F1,F2,
     &                            IBOSYM,IBSPOBEX_TP,INCLUDED_HX,
     &                            LABEL,IPRNT,
     &                            WORK,KFREE,LFREE)
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "ctcc.inc"
#include "symm.inc"
#include "cgas.inc"
#include "orbinp.inc"
#include "glbbas.inc"
#include "clunit.inc"
*
#include "dgroup.h"
*
      character*6 CIRUN
      INTEGER   INCLUDED_HX(NSPOBEX_TP)
      dimension NICA_STR(MXPNGAS),NICB_STR(MXPNGAS),
     &          NIAA_STR(MXPNGAS),NIAB_STR(MXPNGAS)
      dimension ISOBEX(NGAS,4,NSPOBEX_TP),
     &          IBSPOBEX_TP(NSPOBEX_TP),LABEL(4,NSPOBEX_TP)
      dimension I_ADD(4),INDEX(4),IND_SAVE(4)
*
* Group notation of the strings:
      dimension IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS),
     &          IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
      dimension T(NALLINT)
*
* Space for the construction of SUBSTRINGS
* dimensioned outside as integers MX_ST_TSOSO_BLK_MX
      dimension ISTR_CA(*),ISTR_CB(*),ISTR_AA(*),ISTR_AB(*)
* Always remember this when work is passed down
      DIMENSION WORK(*)
*
      NTESTL = 1
C     NTEST = max(NTESTL,IPRNT)
      NTEST = 1
*
* Constants
      IZERO = 0
      ZERO = 0.D0
*
* Double group type
      if (NZ.eq.1) then
        IDGTYP = 1
        write(6,*)
        write(6,*) ' NOTICE from DIST_INTS_SIGDEN: '
        write(6,*) ' Imaginary parts treated, but zeroed out. '
      else if (NZ.eq.2) then
        IDGTYP = 2
      else if (NZ.eq.4) then
        IDGTYP = 3
      end if
*
* Check for inversion symmetry of group
      INVSYM = 0
      if (DOUGRP.eq.7.or.DOUGRP.eq.4.or.DOUGRP.eq.1) INVSYM = 1
*
* Get reordering info depending on environment
C     print*,'INVSYM,DOUGRP,INTIMP',INVSYM,DOUGRP,INTIMP
C       print*,'before reo_info_dirac_KRCC'
C       print*,'NDIM_TOT,NIRR_PN,SYMDIM,IBOSYM',
C    &          NDIM_TOT,NIRR_PN,SYMDIM,IBOSYM
C       print*,'SYDI,DOUGRP,INVSYM,IDGTYP,ISPINFREE',
C    &          SYDI,DOUGRP,INVSYM,IDGTYP,ISPINFREE
C       print*,'KISREO,KINVREO,IPRNT',KISREO,KINVREO,IPRNT
C       CALL MEMCHK_KRCC(WORK)
        call reo_info_dirac_KRCC(NDIM_TOT,NIRR_PN,SYMDIM,IBOSYM,
     &                      SYDI,DOUGRP,WORK(KISREO),WORK(KINVREO),
     &                      INVSYM,IDGTYP,ISPINFREE,IPRNT)
        if (ISPINFREE.eq.1.or.IDGTYP.eq.1) then
          IRILP = 1
        else
          IRILP = 2
        end if
      IRILP = MIN(2,NZ)
* !! Still enforced
*
* Open ASCII unit for checking integral sorting information
      if (NTEST.ge.1) then
        LUINGR = LUSC81
        open(unit=LUINGR,file='INTEGRALS',
     &       status='REPLACE',form='FORMATTED',iostat=ios)
        print*,'ios,LUSC81,LUINGR',ios,LUSC81,LUINGR
        if(ios.ne.0) then
          call abend2('problem with file integrals')
        end if
        rewind(LUINGR)
      end if
*
* File for saving index sequence
* changed from status UNKNOWN to REPLACE to accomedate
* for several calculations
      LUINSV = LUSC82
      open(unit=LUINSV,file='INDICES',
     &     status='REPLACE',form='FORMATTED',iostat=ios)
      if(ios.ne.0) then
        call abend2('problem with file indices')
      end if
      rewind(LUINSV)
*
* Ugly fix to see if the integrals are there from the ARDUCCA code
*
      IMPORT_ARDUCCA = 0
      CALL IMPORT_ARDUCCA_INTEGRALS(T,IMPORT_ARDUCCA)
      IF(IMPORT_ARDUCCA.EQ.1) THEN
        WRITE(6,*) ' INTEGRALS TAKEN FROM ARDUCCA'
      END IF
C       CALL MEMCHK_KRCC(WORK)
*
* Loop over real/imaginary part
      ICOUNT = 1
      ICOUNT2 = 1
      IISUM = 0
      INAC = 0
      INBC = 0
      INAA = 0
      INBA = 0
      ICOUNTRELTYPE = 0
      do IRI=1,IRILP,1
        if (NTEST.ge.1) then
          write(6,*)
          if (IRI.eq.1.and.NTEST.ge.2) write(6,*) 'REAL LOOP'
          if (IRI.eq.1) write(LUINGR,'(A)') 'Real  integrals.'
          if (IRI.eq.2.and.NTEST.ge.2) write(6,*) 'IMAG LOOP'
          if (IRI.eq.2) write(LUINGR,'(A)') 'Imag. integrals.'
        end if
        ITPSV = 0
        if (NTEST.ge.1) write(LUINGR,'(A)') 'INTEGRAL FILE'
C       CALL MEMCHK_KRCC(WORK)
*
* Loop over excitation types
        do IXTP=1,NSPOBEX_TP,1
*
* Check for correct pointers to type blocks
          if (IRI.eq.1.and.ICOUNT.ne.IBSPOBEX_TP(IXTP)) then
            write(6,*) 'Pointers to types do not match numbers'
            write(6,*) 'of integrals generated.'
            write(6,*) 'LABEL array: ',(LABEL(I,IXTP),I=1,4,1)
            write(6,*) 'int.no.,offset: ',ICOUNT,IBSPOBEX_TP(IXTP)
            write(6,*) 'Current type: ',IXTP
            call abend2('Quitting.')
          end if
          if (NTEST.ge.2) then
            write(6,*) '======================'
            write(6,*) 'Processing type ',IXTP
            write(6,*) '======================'
            write(6,*) 'ISOBEX partition:'
            call iwrtma(ISOBEX(1,1,IXTP),NGAS,4,NGAS,4)
          end if
C       CALL MEMCHK_KRCC(WORK)
C
* write statement about next integral type (to file)
C
C Notice at the moment TEST must be at least 1 for newccv to work. Lasse
C
            if (NTEST.ge.1) THEN
              call integral_type_KRCC(ISOBEX(1,1,IXTP),
     &                           NGAS,ITPSV,LUINGR,
     &                           ISUM,NAC,NBC,NAA,NBA)
              CALL ICHANGE_INT_TYPE(ISUM,NAC,NBC,NAA,NBA,
     &                              IISUM,INAC,INBC,INAA,INBA,
     &                              ICHANGE)
              IF(ICHANGE.EQ.1) THEN
                ICOUNTRELTYPE = ICOUNTRELTYPE + 1
                IRELTYPE(ICOUNTRELTYPE) = IXTP
C             print*,'IRELTYPE(ICOUNTRELTYPE)',IRELTYPE(ICOUNTRELTYPE)
C             print*,'ICOUNTRELTYPE',ICOUNTRELTYPE
              END IF
            END IF
            if (NTEST.ge.1)
     &        write(LUINGR,'(A,I5)') 'Type ',IXTP
* New addition by Lasse 12-05-10
C         print*,'INCLUDED_HX(IXTP),IXTP',INCLUDED_HX(IXTP),IXTP
          IF(INCLUDED_HX(IXTP).EQ.1) cycle
*
*. Change from Occupation to group notation of strings
          call occ_to_grp_KRCC(ISOBEX(1,1,IXTP),IGRP_CA,1)
          call occ_to_grp_KRCC(ISOBEX(1,2,IXTP),IGRP_CB,1)
          call occ_to_grp_KRCC(ISOBEX(1,3,IXTP),IGRP_AA,1)
          call occ_to_grp_KRCC(ISOBEX(1,4,IXTP),IGRP_AB,1)
          if (NTEST.ge.5.and.IRI.eq.1) then
            write(6,*) 'Group information:'
            write(6,'(A,4I6)') 'IGRP_CA : ',(IGRP_CA(I),I=1,NGAS,1)
            write(6,'(A,4I6)') 'IGRP_CB : ',(IGRP_CB(I),I=1,NGAS,1)
            write(6,'(A,4I6)') 'IGRP_AA : ',(IGRP_AA(I),I=1,NGAS,1)
            write(6,'(A,4I6)') 'IGRP_AB : ',(IGRP_AB(I),I=1,NGAS,1)
          end if
*
* Corresponding number of electrons
          NEL_CA = ielsum(ISOBEX(1,1,IXTP),NGAS)
          NEL_CB = ielsum(ISOBEX(1,2,IXTP),NGAS)
          NEL_AA = ielsum(ISOBEX(1,3,IXTP),NGAS)
          NEL_AB = ielsum(ISOBEX(1,4,IXTP),NGAS)
          if (NTEST.ge.2.and.IRI.eq.1) then
            write(6,*) 'NEL_CA = ',NEL_CA
            write(6,*) 'NEL_CB = ',NEL_CB
            write(6,*) 'NEL_AA = ',NEL_AA
            write(6,*) 'NEL_AB = ',NEL_AB
          end if
*
* Symmetry loops
* (Annihilator symm. is assumed in terms of corresponding creators!)
*  Symmetry restrictions according to number of electrons
          NEL_C = NEL_CA + NEL_CB
          NEL_A = NEL_AA + NEL_AB
          call evenodd(IEONELC,NEL_C)
          if (IEONELC.eq.1) then
            ISMCSTA = NIRR_PN + 1
            ISMCEND = NIRR_DG
          else if (IEONELC.eq.2) then
            ISMCSTA = 1
            ISMCEND = NIRR_PN
          end if
*
          call evenodd(IEONELCA,NEL_CA)
          if (IEONELCA.eq.1) then
            ISMCASTA = NIRR_PN + 1
            ISMCAEND = NIRR_DG
          else if (IEONELCA.eq.2) then
            ISMCASTA = 1
            ISMCAEND = NIRR_PN
          end if
*
          if (NTEST.ge.5.and.IRI.eq.1) then
            write(6,*) 'C  symmetry loop: ',ISMCSTA,ISMCEND
            write(6,*) 'CA symmetry loop: ',ISMCASTA,ISMCAEND
          end if
*
* Symmetry loops
* Symmetry of an elementary second-quantized operator or a string of such
* operators is defined referring to CREATION operators (particle). Then,
* consequently, annihilation operator symmetries are complex conjugate.
          do ISM_C=ISMCSTA,ISMCEND,1
            ISM_A = IDBGMULT(ISM_C,IADJSYM(ISYM_T))
            do ISM_CA=ISMCASTA,ISMCAEND,1
              ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
              do ISM_AA=1,NIRR_DG,1
                ISM_AB = IDBGMULT(ISM_A,INVELM(ISM_AA))
                if (NTEST.ge.2000.and.IRI.eq.1) then
                  write(6,*)
                  write(6,*)
                  write(6,*) 'Printing symmetries:'
                  write(6,*) 'Creator symm.      :',ISM_C
                  write(6,*) 'Annihilator symm.  :',ISM_A
                  write(6,*) 'Creator alpha symm.:',ISM_CA
                  write(6,*) 'Creator beta symm. :',ISM_CB
                  write(6,*) 'Annihi. alpha symm.:',ISM_AA
                  write(6,*) 'Annihi. beta symm. :',ISM_AB
                  write(6,*) 'Hamiltonian symm.  :',ISYM_T
                  write(6,*) 'Adj. Ham. symm.    :',IADJSYM(ISYM_T)
                  write(6,*) 'more',INVELM(ISM_CA),INVELM(ISM_AA)
                end if
C               CALL MEMCHK_KRCC(WORK)
*
*. Generate strings of these symmetries
*. obtain strings
                IUB = 1
                call getstr2_totsm_spgp_KRCC(IUB,IGRP_CA,NGAS,ISM_CA,
     &                                      NEL_CA,NSTR_CA,ISTR_CA,
     &                                      NORBT,0,IDUM,IDUM,
     &                                      WORK,KFREE,LFREE)
                IUB = 2
                call getstr2_totsm_spgp_KRCC(IUB,IGRP_CB,NGAS,ISM_CB,
     &                                      NEL_CB,NSTR_CB,ISTR_CB,
     &                                      NORBT,0,IDUM,IDUM,
     &                                      WORK,KFREE,LFREE)
                IUB = 1
                call getstr2_totsm_spgp_KRCC(IUB,IGRP_AA,NGAS,ISM_AA,
     &                                      NEL_AA,NSTR_AA,ISTR_AA,
     &                                      NORBT,0,IDUM,IDUM,
     &                                      WORK,KFREE,LFREE)
                IUB = 2
                call getstr2_totsm_spgp_KRCC(IUB,IGRP_AB,NGAS,ISM_AB,
     &                                      NEL_AB,NSTR_AB,ISTR_AB,
     &                                      NORBT,0,IDUM,IDUM,
     &                                      WORK,KFREE,LFREE)
*
* Loops over strings for the given symmetries and types
* Loop over T elements as matrix T(I_CA,I_CB,I_AA,I_AB)
*  Determine possible connections
C               print*,'NSTR_AB,NSTR_AA,NSTR_CB,NSTR_CA',
C    &                  NSTR_AB,NSTR_AA,NSTR_CB,NSTR_CA
                if ((NSTR_CA.ge.1.or.NSTR_CB.ge.1).and.
     &              (NSTR_AA.ge.1.or.NSTR_AB.ge.1)) then
*
                  do I_AB=1,NSTR_AB,1
                    do I_AA=1,NSTR_AA,1
                      do I_CB=1,NSTR_CB,1
                        do I_CA=1,NSTR_CA,1
                          IF(IMPORT_ARDUCCA.EQ.0) THEN
*
* Notes:
*  1) Standard notation in second quantization for 2-el. ints is:
*     a+I a+K_ aL aJ_ (IJ|KL) = - a+I a+K_ aJ_ aL (IJ|KL)
*     The code proceeds according to the right hand side. The
*     signs are located in WORK(KSIGNNHX) for each excitation
*     class. So first indices 2 and 3 are interchanged, then
*     3 and 4 (implicit in the code).
*  2) Integral classes with 2 alpha or 2 beta creators (or both)
*     are treated differently: (due to exchange character)
*     a+I a+K aL aJ [(IJ|KL) - (IL|KJ)]
*  3) In the full spin-dependent formalism, extra classes of
*     integrals are required.
*
*. Indices of STRING I_CA are stored in IOCC_CA and first element is
*  IOCC_CA(1+NEL_CA*(I_CA-1))
*
C                         print*,'NEL_C',NEL_C
C                         print*,'NEL_CA,NEL_CB',NEL_CA,NEL_CB
C                         print*,'NEL_AA,NEL_AB',NEL_AA,NEL_AB
* One-electron part
                          if (NEL_C.eq.1) then
                            IXCH = 0
*
                            if (NEL_CA.eq.1.and.NEL_CB.eq.0) then
                              INDEX(1)
     &                              = ISTR_CA(1+NEL_CA*(I_CA-1))
                              INDEX(2) = 0
                            else if (NEL_CA.eq.0.and.NEL_CB.eq.1) then
                              INDEX(1) = 0
                              INDEX(2)
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                            else
                              write(6,*)
     &                              'Impossible number of creators.'
                              write(6,*)
     &                              'NEL_CA,NEL_CB = ',NEL_CA,NEL_CB
                              call abend2('Quitting.')
                            end if
*
                            if (NEL_AA.eq.1.and.NEL_AB.eq.0) then
                              INDEX(3)
     &                              = ISTR_AA(1+NEL_AA*(I_AA-1))
                              INDEX(4) = 0
                            else if (NEL_AA.eq.0.and.NEL_AB.eq.1) then
                              INDEX(3) = 0
                              INDEX(4)
     &                              = ISTR_AB(1+NEL_AB*(I_AB-1))
                            else
                              write(6,*)
     &                              'Impossible number of annihil.'
                              write(6,*)
     &                              'NEL_AA,NEL_AB = ',NEL_AA,NEL_AB
                              call abend2('Quitting.')
                            end if
* Two-electron part
                          else if (NEL_C.eq.2) then
*
* Will now make exchange on all integrals. Lasse
*
                            if (NEL_CA.eq.2.and.NEL_CB.eq.0) then
                              IXCH = 1
                              do IEL=1,NEL_CA,1
                                INDEX(IEL) =ISTR_CA(IEL+NEL_CA*(I_CA-1))
                              end do
                            else if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
                              IXCH = 1
                              INDEX(1) = ISTR_CA(1+NEL_CA*(I_CA-1))
                              INDEX(3) = ISTR_CB(1+NEL_CB*(I_CB-1))
                            else if (NEL_CA.eq.0.and.NEL_CB.eq.2) then
                              IXCH = 1
                              do IEL=1,NEL_CB,1
                                INDEX(IEL) =ISTR_CB(IEL+NEL_CB*(I_CB-1))
                              end do
                            end if
*
                            if (NEL_AA.eq.2.and.NEL_AB.eq.0) then
                              IXCH = 1
                              do IEL=1,NEL_AA,1
                               INDEX(IEL+2)=ISTR_AA(IEL+NEL_AA*(I_AA-1))
                              end do
                            else if (NEL_AA.eq.1.and.NEL_AB.eq.1) then
                              IXCH = 1
                              INDEX(2) = ISTR_AA(1+NEL_AA*(I_AA-1))
                              INDEX(4) = ISTR_AB(1+NEL_AB*(I_AB-1))
                            else if (NEL_AA.eq.0.and.NEL_AB.eq.2) then
                              IXCH = 1
                              do IEL=1,NEL_AB,1
                                INDEX(IEL+2) =
     &                          ISTR_AB(IEL+NEL_AB*(I_AB-1))
                              end do
                            else
                              write(6,*)
     &                              'Impossible number of annihil.'
                              write(6,*)
     &                              'NEL_AA,NEL_AB = ',NEL_AA,NEL_AB
                              call abend2('Quitting.')
                            end if
*
                          else
                            write(6,*) 'Impossible number of creators.'
                            write(6,*) 'NEL_C = ',NEL_C
                            call abend2('Quitting.')
                          end if
*
*  Stored integrals (DIRAC) refer to symmetry-type ordering.
*  Internal ordering is type-symmetry, so rename:
                          IR_IND = 0
                          JR_IND = 0
                          KR_IND = 0
                          LR_IND = 0
C     DO IJ =1,4
C     print*,'IREOTS(INDEX(IJ)),INDEX(IJ),IJ',
C    &        IREOTS(INDEX(IJ)),INDEX(IJ),IJ
C     END DO
                          if (INDEX(1).ne.IZERO)
     &                              IR_IND = IREOTS(INDEX(1))
                          if (INDEX(2).ne.IZERO)
     &                              JR_IND = IREOTS(INDEX(2))
                          if (INDEX(3).ne.IZERO)
     &                              KR_IND = IREOTS(INDEX(3))
                          if (INDEX(4).ne.IZERO)
     &                              LR_IND = IREOTS(INDEX(4))
                          INDEX(1) = IR_IND
                          INDEX(2) = JR_IND
                          INDEX(3) = KR_IND
                          INDEX(4) = LR_IND
C      print*,'I,J,K,L',IR_IND,JR_IND,KR_IND,LR_IND
*
* Index reordering and second integral in exchange case:
* (Reason: CACA string generated e.g. for two alpha excitations,
*          but this does not correspond to index ordering in integral!)
*  Apply only to 2-electron integrals:
*
                          YXINT = ZERO
                          YINT = ZERO
C                         print*,'writing out IXCH',IXCH
                          if (IXCH.eq.1) then
* Index sequence reads:
*  (IJ|KL) (which is not the operator sequence for the desired int!)
* Reorder for Coulomb integral:
*  (IJ|KL) --> (IK|LJ)
* Now create exchange integral:
*  (IJ|KL) --> (IJ|LK)
* (Overall, indices L and J are interchanged and a "-" introduced,
*  due to I+ K+ L J = - I+ K+ J L
                            IF(NEL_CA.NE.1) THEN
* Coulomb:
                            I_SV = INDEX(1)
                            L_SV = INDEX(2)
                            K_SV = INDEX(3)
                            J_SV = INDEX(4)
                            INDEX(1) = I_SV
                            INDEX(2) = K_SV
                            INDEX(3) = L_SV
                            INDEX(4) = J_SV
* Exchange:
                            I_ADD(1) = I_SV!LABEL(1,IXTP)
                            I_ADD(2) = J_SV!LABEL(4,IXTP)
                            I_ADD(3) = L_SV!LABEL(3,IXTP)
                            I_ADD(4) = K_SV!LABEL(2,IXTP)
                            do IRUN = 1,4,1
                              I_ADD(IRUN) =
     &                              I_ADD(IRUN) *
     &                              LABEL(IRUN,IXTP)
C                             print*,'I_ADD(IRUN),IRUN',I_ADD(IRUN),IRUN
                            end do
                            ELSE
C           write(6,'(A,6X,4I3)') '1 INDEX array: ',(I_ADD(ID),ID=1,4,1)
                            I_ADD(1) = INDEX(1)*LABEL(1,IXTP)
                            I_ADD(4) = INDEX(2)*LABEL(2,IXTP)
                            I_ADD(3) = INDEX(3)*LABEL(3,IXTP)
                            I_ADD(2) = INDEX(4)*LABEL(4,IXTP)
C           write(6,'(A,6X,4I3)') '2 INDEX array: ',(I_ADD(ID),ID=1,4,1)
                            END IF
C                           print*,'IRI,NZ',IRI,NZ
                            if (IRI.eq.1.or.(IRI.eq.2.and.NZ.ne.1)) then
                                call getint_sigden_KRCC(YXINT,I_ADD,IRI,
     &                                    NDIM_TOT,NASHT,F1,F2,IPRNT,
     &                                    WORK,KFREE,LFREE)
                            end if
                          end if
*
*
*  Modify if integrals are over relativistic spinors:
                            if ((INTIMP.eq.6.or.INTIMP.eq.7).and.
     &                           ISPINFREE.eq.0) then
                              if (INDEX(1).ne.0.and.
     &                            INDEX(2).ne.0) then
                                do IRUN = 1,4,1
                                    INDEX(IRUN) =
     &                                  INDEX(IRUN) *
     &                                  LABEL(IRUN,IXTP)
                                end do
*
* Correct internal ordering of two-bar exchange type integrals:
* (This might be implemented in some more elegant way later on.)
                                if (INDEX(1).gt.0.and.
     &                              INDEX(2).lt.0.and.
     &                              INDEX(3).lt.0.and.
     &                              INDEX(4).gt.0) then
                                  J_SV = INDEX(2)
                                  L_SV = INDEX(4)
                                  INDEX(2) = -L_SV
                                  INDEX(4) = -J_SV
                                end if
                              end if
*
*                           Write INDEX sequence to line in file
*                           for use in diagonal routine (saving memory):
                            write(LUINSV,'(4I4)') (INDEX(JJ),JJ=1,4)
*
C                           print*,'IRI,NZ',IRI,NZ
                            if (IRI.eq.1.or.(IRI.eq.2.and.NZ.ne.1)) then
C           write(6,'(A,6X,4I3)') '3 INDEX array: ',(INDEX(ID),ID=1,4,1)
                              call getint_sigden_KRCC(YINT,INDEX(1),
     &                                  IRI,NDIM_TOT,NASHT,F1,F2,IPRNT,
     &                                  WORK,KFREE,LFREE)
                            end if
*
C                       print*,'ICOUNT,YINT,YXINT',ICOUNT,YINT,YXINT
                            T(ICOUNT) = YINT - YXINT
*
* FIXME !! Uncommented sign here!! 
C                          print*,'IXCH',IXCH
                            if (IXCH.eq.1.AND.NEL_CA.NE.1) 
     &                          T(ICOUNT) = - T(ICOUNT)
* FIXME !! Uncommented sign here!!
*
                            if (NTEST.ge.1) then
                              if (NTEST.ge.10)
     &                            write(6,'(A,I4,A,4I3,A,1F18.12)')
     &                                  'INTEGRAL no. ',ICOUNT,'(',
     &                                   (INDEX(M),M=1,4),
     &                                  ') Value: ',T(ICOUNT)
*                           write integral to file (for debugging):
                              if (IXCH.eq.0) then
                                write(LUINGR,'(I4,4X,4I3,1F18.10)')
     &                                ICOUNT,INDEX(1),INDEX(2),
     &                                       INDEX(3),INDEX(4),T(ICOUNT)
                              else
                           write(LUINGR,'(I4,4X,4I3,X,A,X,4I3,1F18.10)')
     &                                ICOUNT,INDEX(1),INDEX(2),
     &                                       INDEX(3),INDEX(4),'-',
     &                             I_ADD(1),I_ADD(2),I_ADD(3),I_ADD(4),
     &                                       T(ICOUNT)
                              end if
                            end if
*
*                           Write INDEX sequence to line in file
*                           for use in diagonal routine (saving memory):
                            write(LUINSV,'(4I4)') (INDEX(JJ),JJ=1,4)
                          end if
                          END IF
                          ICOUNT = ICOUNT + 1
                          ICOUNT2 = ICOUNT2 + 1
                        end do
*                       ^ End loop over alpha creation strings
                          IF(IRI.EQ.2) THEN
                            IXTPCOUNT = IXTP+NSPOBEX_TP
                          ELSE
                            IXTPCOUNT = IXTP
                          END IF
                          IF(IXTPCOUNT.EQ.1) THEN
                          IOFFINTTYPE(IXTPCOUNT) = 1
                          IOFFINTTYPE(IXTPCOUNT+1) = ICOUNT
                          ELSE IF(IXTPCOUNT.LT.2*NSPOBEX_TP) THEN
                          IOFFINTTYPE(IXTPCOUNT+1) = ICOUNT
                          END IF
C     print*,'IXTPCOUNT,ICOUNT,IXTP',IXTPCOUNT,ICOUNT,IXTP
                      end do
*                     ^ End loop over beta creation strings
                    end do
*                   ^ End loop over alpha annihilation strings
                  end do
*                 ^ End loop over beta annihilation strings
                end if
*               ^ End if connection condition
              end do
*             ^ End loop over alpha annihilator symmetry
            end do
*           ^ End loop over alpha creator symmetry
          end do
*         ^ End loop over creation part symmetry
        end do
*       ^ End loop over excitation types
      end do
*     ^ End loop over real/imaginary
      if (NTEST.ge.1.and.(CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 '))
     &   close(LUINGR)
      close(LUINSV)
      NALLINT = ICOUNT - 1
*
C Since IOFFINTTYPE may not be correctly set when there is only one
C Kramers pair in a gas space we need to insure that there is always
C an increasing number for the next HTYPE in IOFFINTTYPE
      CALL NXT_NUM_GE_NUM_MINUS_ONE(IOFFINTTYPE,2*NSPOBEX_TP)
*
C     CALL MEMCHK_KRCC(WORK)
      NTEST = 100
      if (NTEST.ge.1) then
        write(6,*) ' DIST_INTS_SIGDEN reporting   : '
        write(6,*) '  Total number of integrals   : ',NALLINT
      end if
      IF(NTEST.GE.100) THEN
        DO I=1,10
          WRITE(6,*) ' Change in integral type ',IRELTYPE(I)
        END DO
        DO I=1,NSPOBEX_TP
          WRITE(6,*) ' INTEGRAL OFFSET ',IOFFINTTYPE(I),' FOR HTYPE ',I
        END DO
      END IF
      ICHECK = 1
      IF(ICHECK.EQ.1) THEN
        open(unit=54,file='INTE',
     &       status='REPLACE',form='FORMATTED',iostat=ios)
        if(ios.ne.0) then
          call abend2('problem with file integrals')
        end if
        rewind(54)
        DO I=1,NALLINT
        WRITE(54,*) T(I)
        END DO
      END IF
*
      return
      end
*
       SUBROUTINE ICHANGE_INT_TYPE(ISUM,NAC,NBC,NAA,NBA,
     &                            IISUM,INAC,INBC,INAA,INBA,ICHANGE)
C
C Routine to see if there has been a change in the type of integral
C
#include "implicit.inc"
C from 1 to
      NTEST = 00
      ICHANGE = 0
      IF(ISUM.NE.IISUM) ICHANGE = 1
C
      IF(NAC.NE.INAC.OR.NBC.NE.INBC.OR.NAA.NE.INAA.OR.NBA.NE.INBA) THEN
        ICHANGE = 1
      END IF
C reset
      INAC = NAC
      INBC = NBC
      INAA = NAA
      INBA = NBA
      IISUM = ISUM
C
      IF(NTEST.GE.100) THEN
        IF(ICHANGE.GE.1) THEN
          WRITE(6,*) 'A new type has been found'
          WRITE(6,*) 'NAC,NBC,NAA,NBA',NAC,NBC,NAA,NBA
          WRITE(6,*) 'INAC,INBC,INAA,INBA',INAC,INBC,INAA,INBA
        ELSE
          WRITE(6,*) 'No new type has been found'
          WRITE(6,*) 'NAC,NBC,NAA,NBA',NAC,NBC,NAA,NBA
        END IF
      END IF
      RETURN
      END
*
* Determine matching dimensions of boson functions,
* and reordering arrays for DIRAC environment.
*
      subroutine reo_info_dirac_KRCC(ITOTDIM,NIRR_PN,SYMDIM,IBOSYM,
     &                          SYDI,DOUGRP,ISREO,INVREO,INVSYM,
     &                          IDGTYP,ISPINFREE,IPRNT)
*
#include "implicit.h"
#include "mxpdim.inc"
*
      integer SYDI(NIRR_PN),SYMDIM(NIRR_PN),DOUGRP
      dimension ISREO(2*ITOTDIM),INVREO(-ITOTDIM:ITOTDIM),
     &          IBOSYM(2*ITOTDIM)
*
      NTESTL = 00
      NTEST = max(IPRNT,NTESTL)
*
* Check for matching dimensions of boson symmetry functions
      if (ISPINFREE.eq.1) then
        do IRR=1,NIRR_PN,1
          IBCT = 0
          IOFFBS = 0
          do ISP=1,ITOTDIM,1
            if (IBOSYM(ISP+IOFFBS)+1.eq.IRR) IBCT = IBCT + 1
            if (INVSYM.eq.1.and.ISP.eq.SYMDIM(1)/2) IOFFBS = SYMDIM(1)/2
          end do
*
          if (IBCT.ne.SYDI(IRR)) then
            write(6,*)
            write(6,*) 'Problem with IBSPI array.'
            write(6,*) 'Boson symmetry dimensions do not match'
            write(6,*) 'with GAS matrix input.'
            write(6,*) 'Symmetry ',IRR
            write(6,*) 'No. of input functions ',SYDI(IRR)
            write(6,*) 'No. of functions on IBSPI array ',IBCT
            call abend('Quitting.')
          end if
        end do
        write(6,*)
        write(6,*) 'Boson symm. dimension check passed.'
*
      else if (ISPINFREE.eq.0) then
*
* Set IBOSYM array so that U->alpha and B->beta spinors
*              To get correct Fermion symmetry labeling in
*              LUCIAREL we need to specify same(!) gerade boson
*              irrep for all gerade orbitals and same ungerade
*              boson irrep for all ungerade orbitals (note: which
*              gerade or ungerade boson irrep we choose doesn't
*              matter).
        IOFFFS = 0
        IRROFF = 0
        do IFSYM = 1,NIRR_PN,1
          do ISP = 1+IOFFFS,2*SYMDIM(IFSYM)+IOFFFS,1
            IBOSYM(ISP) = 0 + IRROFF
          end do
*
          if (NTEST.ge.3) then
            write(6,*) ' Reset IBOSYM array, ferm symm. :',IFSYM
            write(6,*) (IBOSYM(I),I=1+IOFFFS,2*SYMDIM(IFSYM)+IOFFFS,1)
          end if
*
          IOFFFS = 2*SYMDIM(IFSYM)
          IRROFF = 2
          if (IDGTYP.eq.3) IRROFF = 1
        end do
      else
        write(6,*) 'ISPINFREE = ',ISPINFREE
        write(6,*) ' has illegal value.'
        stop 'Quitting.'
      end if
*
*  Create symmetry reordering array for this list
*  ... and the inverted array
      ICT = 1
      JCT = 1
      IOFFBS = 0
      do INSYM=1,NIRR_PN,1
        if (INSYM.le.2) then
          IDIMGRP = 1
        else
          IDIMGRP = 2
        end if
        IDIV = 1
        if (IDGTYP.eq.3) then
          IDIMGRP = INSYM
          IDIV = 2
        end if
*
* ... continue with offsets and skipping of even spinor
*     symmetries (these have redundant, i.e. Kramers-paired
*     classification
* For the moment, we hardwire the irrep order of C2h to the
*  generators XY and XYZ (necessarily in this order !!).
*  I.e., the irreps are classified
*  as Ag, Bg, Au, Bu. This should be made more flexible later on.
        ILEN = SYMDIM(IDIMGRP) / IDIV
        if (INSYM.ge.3) IOFFBS = 2*SYMDIM(1)
        if (DOUGRP.eq.4) then
*  INSYM_SEA corresponds to the irrep ordering generated by
*  the choice of symmetry generators in the integral program.
          if (INSYM.eq.1) INSYM_SEA = 1
          if (INSYM.eq.2) INSYM_SEA = 2
          if (INSYM.eq.3) INSYM_SEA = 3
          if (INSYM.eq.4) INSYM_SEA = 4
        else
          INSYM_SEA = INSYM
        end if
        if (INSYM.eq.3) ICT = 1
        if (NTEST.ge.10) then
          write(6,*) 'Current run       : ',INSYM
          write(6,*) 'Seeking symmetry  : ',INSYM_SEA
          write(6,*) 'Doing IND from ',1+IOFFBS, 'to',IOFFBS+ILEN
        end if
        do IND=1+IOFFBS,IOFFBS+ILEN,1
          if ((IBOSYM(IND)+1).eq.INSYM_SEA) then
            ISREO(IND) = ICT
            IND_REL = IND
            if (INSYM.ge.3) IND_REL = IND - SYMDIM(1)
            INVREO(JCT) = IND_REL
            INVREO(-JCT) = -IND_REL
            JCT = JCT + 1
            ICT = ICT + 1
          end if
        end do
      end do
*  ISREO(orbital number DIRAC) = boson symmetry ordered orbital number
*        for one-electron part
*  INVREO(boson symmetry ordered orbital number) = orbital number DIRAC
*        for two-electron part
*
*
      if (NTEST.ge.3) then
        write(6,*) 'Symmetry reordering array DIRAC:'
        do IGRP=1,NIRR_PN/2,1
          if (IGRP.eq.1) then
            write(6,*) 'First dimension group:'
            call iwrtma(ISREO(1),1,SYMDIM(IGRP),1,SYMDIM(IGRP))
          else if (IGRP.eq.2) then
            IFEROFF = 2*SYMDIM(1)
            write(6,*) 'Second dimension group:'
            call iwrtma(ISREO(1+IFEROFF),1,SYMDIM(IGRP),1,SYMDIM(IGRP))
            write(6,*) 'Ferm. index offset for second IDIMGRP: ',IFEROFF
          end if
        end do
        write(6,*) 'Inverse reordering array (One Kramers branch):'
        call iwrtma(INVREO(1),1,ITOTDIM,1,ITOTDIM)
      end if
*
      return
      end
*
      subroutine integral_type_KRCC(ISOB,NGAS,ITPSV,LUINGR,
     &                              ISUM,NAC,NBC,NAA,NBA)
*
#include "implicit.h"
      character*11 CELEC2
      character*13 CELEC
      dimension ISOB(NGAS,4)
      ISUM = 0
      do J=1,4,1
        do I=1,NGAS,1
          ISUM = ISUM + ISOB(I,J)
        end do
      end do
      if (ISUM.eq.2) CELEC = 'One-electron '
      if (ISUM.eq.4) CELEC = 'Two-electron '
      NAC = 0
      NBC = 0
      NAA = 0
      NBA = 0
      do I=1,NGAS,1
        NAC = NAC + ISOB(I,1)
        NBC = NBC + ISOB(I,2)
        NAA = NAA + ISOB(I,3)
        NBA = NBA + ISOB(I,4)
      end do
      if (NAC.eq.NAA.and.NBC.eq.NBA) CELEC2 = 'spin-free  '
      if (NAC.ne.NAA.or.NBC.ne.NBA)  CELEC2 = 'spin-orbit '
      if (ISUM.ne.ITPSV)
     &     write(LUINGR,'(A,A,A)') CELEC,CELEC2,'integrals.'
      ITPSV = ISUM
*
      NTEST = 00
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Operator in integral_type_KRCC '
        CALL WRT_SPOX_TP_CC_KRCC(ISOB,1)
        WRITE(6,*) ' NAC,NBC,NAA,NBA ',NAC,NBC,NAA,NBA
      END IF
*
      return
      end
*
* Fetch integral from integral list according to SIGDEN
* excitation class formalism - driver routine
*
      subroutine getint_sigden_KRCC(YINT,INDEX,IRI,ITOTDIM,
     &                         NASHT,F1,F2,IPRNT,
     &                         WORK,KFREE,LFREE)
*
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "glbbas.inc"
#include "crun.inc"
#include "symm.inc"
*
      logical NZ3FORM
      dimension INDEX(4)
      DIMENSION WORK(*)
*
* New (memory saving) (NZ,3) integral interface (if true)
      NZ3FORM = .TRUE.
*
C     CALL MEMCHK_KRCC(WORK)
      
!
!         this is the default integral scheme of today's LUCIAREL.
!
            call kkr_set_KRCC(WORK(KKR),ITOTDIM,SYDI,NIRR_PN,INVSYM)
            NNASHX = ((NASHT/2) * ((NASHT/2)+1)) / 2
            call getind_nz3_KRCC(YINT,INDEX,IRI,ITOTDIM,WORK(KINVREO),
     &                      WORK(KKR),NASHT,NNASHX,F1,F2,IPRNT)
*
C     CALL MEMCHK_KRCC(WORK)
      end
*
* Set KR array according to number of functions in fermion
* symmetries.
*
      subroutine kkr_set_KRCC(KR,IDIM,ISYDI,NIRR,INVSYM)
*
#include "implicit.inc"
*
      dimension KR(-IDIM:IDIM),ISYDI(NIRR)
*
      NTESTL = 0
*
      ICT = 1
      IV = 1
      do IRR = 1,1+INVSYM,1
        IRRR = IRR
        if (INVSYM.eq.1.and.NIRR.gt.2.and.IRR.eq.2) IRRR = IRR + 1
        do I = 1,ISYDI(IRRR),1
          KR(IV) = ICT
          KR(-IV) = ICT + ISYDI(IRRR)
          ICT = ICT + 1
          IV = IV + 1
        end do
        ICT = ICT + ISYDI(IRR)
      end do
*
      if (NTESTL.ge.1) then
        write(6,*) ' KKR has set KR array :'
        do ID = 1,IDIM,1
          write(6,'(3X,2I4)') KR(-ID),KR(ID)
        end do
      end if
*
      return
      end
*
* Determine 1- or 2- particle integral
* Fetch integral.
*
*=========================
* Timo Fleig, Oct 20, 2004
*=========================
*
      subroutine getind_nz3_KRCC(YINT,INDEX,IRI,IDIM,INVREO,
     &                      KR,NASHT,NNASHX,DMFC,H2AC,IPRNT)
*
#include "implicit.inc"
#include "dgroup.h"
*
C     dimension DMFC(*),H2AC(*)
      dimension DMFC(*),H2AC(NASHT/2,NASHT/2,NNASHX,NZ,*)
      dimension ICSAVE(4),INDEX(4),INVREO(-IDIM:IDIM),KR(-IDIM:IDIM),
     &          INDEX_LOC(4)
*
      do IL = 1,4,1
        INDEX_LOC(IL) = INDEX(IL)
      end do
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRNT)
*
      if (IRI.eq.1) then
        MZ = 1
      else if (IRI.eq.2) then
        MZ = NZ
      end if
*
      if (NTEST.ge.2) then
        write(6,*)
        write(6,'(A,6X,4I3)') '  INDEX array: ',(INDEX(I),I=1,4,1)
      end if
      if (NTEST.ge.5) then
        write(6,*)
        write(6,*) ' GETIND_NZ3 : '
        write(6,*) '  NZ,MZ : ',NZ,MZ
        write(6,*) '  No. of active orbitals : ',NASHT/2
        write(6,*) '  Compound index dimension: ',NNASHX
      end if
*
      if (NTEST.ge.20) then
        write(6,*)
        write(6,*) '1-el. integral array in getind_nz3:'
        call wrtmat(DMFC,IDIM*2,IDIM*2,IDIM*2,IDIM*2)
        write(6,*)
        write(6,*) '2-el. integrals in getind_nz3, (NZ,3) format:'
        call prdnz3(H2AC(1,1,1,1,1),NASHT/2,NNASHX,NZ,IPQTOQ(1,0),6)
      end if
*
      NCT = 0
      do IC=1,4,1
        if (INDEX_LOC(IC).ne.0) then
          NCT = NCT + 1
          ICSAVE(NCT) = IC
        end if
      end do
*
      if (NCT.eq.2) then
        if (NTEST.ge.10) write(6,*) '1-electron integral.'
        I = KR(INVREO(INDEX_LOC(ICSAVE(1))))
        J = KR(INVREO(INDEX_LOC(ICSAVE(2))))
* Compute input address
        IIN =  (IRI-1)*NASHT*NASHT +
     &          (J-1) *NASHT +
     &           I
        if (NTEST.ge.5) then
          write(6,*) 'Index pair ',I,J,' address :',IIN
        end if
        YINT = DMFC(IIN)
*
* Sign for imaginary 1-electron integrals:
* (According to h_(IJ) = -h_(I_J_).
*  Reference are the unbarred spinors.)
        if (IRI.eq.2.and.INDEX_LOC(1).eq.0) YINT = - YINT
*
        if (NTEST.ge.2) write(6,'(2I3,2X,1F18.10)') I,J,YINT
*
      else if (NCT.eq.4) then
*
* Interface to (NZ,3) format for 2-electron integrals
*
* Use Kramers symmetry/particle symmetry to generate a stored
* index sequence on (NZ,3) format:
*  NZ3TP = 1     (pq|rs)
*  NZ3TP = 2     (pq_|rs_)
*  NZ3TP = 3     (p_q|rs_)
          IRIFAC = 1
          ISEQTP = 0
          if (INDEX_LOC(1).lt.0.and.INDEX_LOC(2).lt.0.and.
     &        INDEX_LOC(3).lt.0.and.INDEX_LOC(4).lt.0) then
            call index_exch(INDEX_LOC(1),INDEX_LOC(2),1,-1,1)
            call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,-1,1)
            NZ3TP = 1
            ISEQTP = 1
          else if (INDEX_LOC(1).gt.0.and.INDEX_LOC(2).gt.0.and.
     &             INDEX_LOC(3).gt.0.and.INDEX_LOC(4).gt.0) then
            NZ3TP = 1
            ISEQTP = 2
          else if (INDEX_LOC(1).gt.0.and.INDEX_LOC(2).gt.0.and.
     &             INDEX_LOC(3).lt.0.and.INDEX_LOC(4).lt.0) then
            call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,-1,1)
            NZ3TP = 1
            ISEQTP = 3
          else if (INDEX_LOC(1).gt.0.and.INDEX_LOC(2).lt.0.and.
     &             INDEX_LOC(3).lt.0.and.INDEX_LOC(4).gt.0) then
            call index_exch(INDEX_LOC(1),INDEX_LOC(3),1,1,1)
            call index_exch(INDEX_LOC(2),INDEX_LOC(4),1,1,1)
            INDEX_LOC(1) = -INDEX_LOC(1)
            INDEX_LOC(4) = -INDEX_LOC(4)
            NZ3TP = 3
            ISEQTP = 4
          else if (INDEX_LOC(1).gt.0.and.INDEX_LOC(2).lt.0.and.
     &             INDEX_LOC(3).gt.0.and.INDEX_LOC(4).lt.0) then
            INDEX_LOC(2) = -INDEX_LOC(2)
            INDEX_LOC(4) = -INDEX_LOC(4)
            NZ3TP = 2
            ISEQTP = 5
          else if (INDEX_LOC(1).lt.0.and.INDEX_LOC(2).gt.0.and.
     &             INDEX_LOC(3).lt.0.and.INDEX_LOC(4).gt.0) then
            if (IRI.eq.2) IRIFAC = -1
            call index_exch(INDEX_LOC(1),INDEX_LOC(2),1,1,1)
            call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,1,1)
            INDEX_LOC(2) = -INDEX_LOC(2)
            INDEX_LOC(4) = -INDEX_LOC(4)
            NZ3TP = 2
            ISEQTP = 6
          end if
          if (NTEST.ge.2) write(6,'(A,X,4I3)')
     &       '  After 3 list step:', (INDEX_LOC(I),I=1,4,1)
*
* Find smallest index and its position in Mulliken notation
          call findlow(INDEX_LOC,4,ILOW,IPOS)
*
* Transfer smallest to fourth index by using
*   particle symmetry
*   Kramers symmetry
*   complex conjugation symmetry (IRIFAC2 = -1 if imaginary part)
*   CAUTION HERE:
*    The sign has to be derived from the original index sequence!
*    Thus the below sign arrangements may look queer.
*    They are nevertheless correct.
*
          IRIFAC2 = 1
          if (IPOS.eq.1) then
            call index_exch(INDEX_LOC(1),INDEX_LOC(3),1,1,1)
            call index_exch(INDEX_LOC(2),INDEX_LOC(4),1,1,1)
            if (IRI.eq.2.and.NZ3TP.ne.2) IRIFAC2 = -1
            call index_exch(INDEX_LOC(1),INDEX_LOC(2),1,1,1)
            call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,1,1)
          else if (IPOS.eq.2) then
            if (IRI.eq.2.and.NZ3TP.eq.3) IRIFAC2 = -1
            call index_exch(INDEX_LOC(1),INDEX_LOC(3),1,1,1)
            call index_exch(INDEX_LOC(2),INDEX_LOC(4),1,1,1)
          else if (IPOS.eq.3) then
            if (IRI.eq.2.and.NZ3TP.eq.1) IRIFAC2 = -1
            call index_exch(INDEX_LOC(1),INDEX_LOC(2),1,1,1)
            call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,1,1)
          end if
          if (NTEST.ge.2) write(6,'(A,X,4I3)')
     &       '  After posit. step:', (INDEX_LOC(I),I=1,4,1)
*
* Determine address in (NZ,3) format
          NASHTH = NASHT/2
*
* Determine compound index KL: NDIM2
          NDIM2 = 0
          do IC = 1,INDEX_LOC(3),1
            if (IC.eq.INDEX_LOC(3)) then
              JCEND = INDEX_LOC(4)
            else if (IC.ne.INDEX_LOC(3)) then
              JCEND = IC
            end if
            do JC = 1,JCEND,1
              NDIM2 = NDIM2 + 1
            end do
          end do
*
          IIN = (NZ3TP-1)        *MZ*NNASHX*NASHTH*NASHTH +
     &          (MZ-1)           *NNASHX*NASHTH*NASHTH +
     &          (NDIM2-1)        *NASHTH*NASHTH +
     &          (INDEX_LOC(1)-1) *NASHTH +
     &           INDEX_LOC(2)
*
C         YINT = H2AC(IIN) * IRIFAC * IRIFAC2
          YINT = H2AC(INDEX_LOC(1),INDEX_LOC(2),NDIM2,MZ,NZ3TP)
     &                * IRIFAC * IRIFAC2
          if (NTEST.ge.2) write(6,'(A,4I3,2X,4I3,2X,1F18.10)')
     &                    'MZ,NZ3TP,IRIFAC,IRIFAC2  integral',
     &                     MZ,NZ3TP,IRIFAC,IRIFAC2,
     &                    (INDEX_LOC(ID),ID=1,4,1),YINT
*
      else
        write(6,*) 'No QFT implemented YET in LUCIAREL.'
        write(6,*) 'you specified ',NCT,' indices.'
        call abend2('Quitting.')
      end if
*
      return
      end
*
      SUBROUTINE RGETH21_KRCC(DOFQ,DOH2AC,PV,FQ,H2AC,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Calculate FQ and H2AC.
C
C     Input:
C        DOFQ  - calculate the auxiliary Fock matrix FQ
C        DOH2AC - collect active two-electron integrals.
C        PV - two-electron density matrix
C
C     Output:
C        FQ - auxiliary Fock matrix
C        (H2AC written to disk)
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : L. Visscher - Mar 15 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "maxorb.h"
#include "dcbidx.h"
#include "dcbham.h"
C
      DIMENSION FQ(NORBT,NASHT,NZ)
      DIMENSION PV(NASHT,NASHT,NNASHX,NZ,3)
      DIMENSION H2AC(NASHT,NASHT,NNASHX,NZ,3)
      DIMENSION WORK(*)
      LOGICAL DOFQ, DOH2AC
C     index arrays
      DIMENSION NFPCK12(2), NFPCK34(2)
      DIMENSION NSTR(2,0:2,4)
      DIMENSION IKLB(2)
C     strings
      CHARACTER DATEX*8,TIMEX*8,FNAM*5,FNODE*10,LFNAM*10
      logical fopen
C
#include "memint.h"
C
      CALL QENTER('RGETH21')
C
C     Read information from integral file
C     -----------------------------------
C
      IPRINT = 10
      LUINT = 1
      OPEN (LUINT, FILE = '4INDINFO', STATUS = 'OLD',
     &     FORM = 'UNFORMATTED' )
C
      WRITE(LUPRI,*)
      WRITE(LUPRI,'(A)') ' Reading from file : 4INDINFO '
C
C     Read time & date for creation
C
      READ (LUINT) DATEX, TIMEX
C
      IF ( IPRINT .GE. 5 ) WRITE(LUPRI,'(4A)')
     &     ' (RGETH2) Time and date read from 4INDINFO: ',
     &     TIMEX, ' ',DATEX
C
      READ (LUINT) NFPCK12, NFPCK12T, NFPCK34, NFPCK34T, NSTR
C
      CALL MEMGET('INTE',KDXKR12,2*NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB12,2*NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB34,2*NFPCK34T,WORK,KFREE,LFREE)
C
C     I and J indices
      CALL READI(LUINT,NFPCK12T*2,WORK(KDXKR12))
C     Boson irrep for IJ indices
      CALL READI(LUINT,NFPCK12T*2,WORK(KDXB12))
C     Boson irrep for KL indices
      CALL READI(LUINT,NFPCK34T*2,WORK(KDXB34))
C
      CLOSE (LUINT, STATUS = 'KEEP')
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Open integral files.
C     --------------------
C
      IIKLB = 0
      LUINT = 56
C     LUINT = 56 + IREPIJ, IREPIJ = 1,NFSYM are actually used
      DO IREPIJ = 1, NFSYM
         IKLB(IREPIJ) = IIKLB
         IRECL = IRECLEN(NCLASS * NFPCK34 (IREPIJ), 0, 0)
         IRECL = MAX(IRECL,IRECLEN(1,0,0))
         WRITE (FNAM,'(A4,I1)') '4IND',IREPIJ
         FNODE = LFNAM(FNAM)
         WRITE(LUPRI,'(A,A)') ' Reading from file : ',FNODE
         OPEN (LUINT+IREPIJ,FILE=FNODE,STATUS='OLD',
     &        ACCESS='DIRECT',RECL=IRECL)
         IIKLB = IIKLB + NFPCK34 (IREPIJ)
      END DO
C
C     Allocate memory for integrals (H2XY)
C
      CALL MEMGET('REAL',KHMAT,N2ORBX*NCLASS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2XY,N2ORBX*NCLASS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2YX,N2ORBX*NCLASS,WORK,KFREE,LFREE)
c
#ifdef MOD_DEBUG
      write(6,*) '(RGETH2) IDXKR12'
      do i = 1, nfpck12t
         write(6,*) I,
     &        IGETINT(WORK(KDXKR12),1+(I-1)*2),
     &        IGETINT(WORK(KDXKR12),I*2)
      end do
#endif
c
C
C     Read integrals (x>y)
C     --------------------
C
C     zero arrays
C
      IF (DOFQ) CALL DZERO(FQ,NASHT*NORBT*NZ)
      IF (DOH2AC) CALL DZERO(H2AC,NASHT*NASHT*NNASHX*NZ*3)
C
C     Loop over records in integrals files.
C
      NXY = 0
      DO 10 IREPIJ = 1, NFSYM
         IREC = 0
         DO 20 IXY = 1, NFPCK12(IREPIJ)
            NXY = NXY + 1
            IREC = IREC + 1
            CALL GTXY(NXY,IX,IY,IXU,IYU,WORK(KDXKR12))
            IREPIJB = IGETINT(WORK(KDXB12),NXY*2-1)
#ifdef MOD_DEBUG
            write(6,*) 'x & y: ',ix,iy
            write(6,*) 'ative : ',ixu,iyu
            write(6,*) ' '
#endif
C           We only need active-active distributions.
C
            IF (IXU .LT. 0 .OR. IYU .LT. 0) GOTO 20
C
C           ...also skip if X<Y
C
            IF (IXU .LT. IYU) GOTO 20
C
            IUU = MAX(IXU,IYU) * ( MAX(IXU,IYU) - 1 ) / 2
     &           + MIN(IXU,IYU)
C
C           Read integrals
C
            NKL = NFPCK34(IREPIJ)
C
            READ (LUINT + IREPIJ, REC=IREC,
     &           ERR=8000, IOSTAT=IOSVAL)
     &           (WORK(KH2XY+I-1),I=1,NKL*NCLASS)
#ifdef MOD_DEBUG
            ii = 0
            do i = 1, 4 * nz
               write(6,*) 'class ',i,' before transformation'
               call output(work(kh2xy+ii),1,norbt,1,norbt,
     &              norbt,norbt,1,lupri)
               ii = ii + nkl
            end do
#endif
C
C           Reorder integrals to full matrix
C           --------------------------------
C
            CALL DZERO(WORK(KH2YX),N2ORBX*NCLASS)
C
            IREPABR = 0
            DO IREPAB = 0, NBSYM - 1
               IF (IREPIJ .EQ. JBTOF(IREPAB,1)) THEN
                  IREPABR = IREPABR + 1
                  DO IZ2 = 1, NZ
                     IQ2 = IPQTOQ(IZ2,IREPAB)
C
C                    Transform symmetry packed
C                    integrals to full matrix.
C                    (take care of diff. integral transformations etc.)
C
                     CALL INTB2X(WORK(KH2XY),IZ2,IREPABR,IREPAB,
     &                    WORK(KH2YX),IREPIJ,NKL)
C                    CALL INTB2X(HT,IZ2,IREPABR,IREPAB,HX,IREPIJ,NKL)
C
C                    Transform boson-symmetry array to
C                    full
C
C                    (FIXME: write me)
C
                     IF (SPINFR) THEN
                        write(6,*) 'SPINFR not implemented in RGETH2'
                        CALL QUIT('SPINFR not implemented in RGETH2')
                     END IF
C
                  END DO
               END IF
            END DO
C
C           Transform to Dirac (NZ,3) format
C           --------------------------------
C
            CALL CL2NZ3(WORK(KH2YX),WORK(KHMAT),N2ORBX,IREPIJ,
     &           IREPIJB,WORK(KDXB34),IKLB)
C           CALL CL2NZ3(VCL,VNZ3,NELEM,IREPIJ,IREPIJB,IREPKLB,IKLB)
C
C           Get (gu) and (ug) blocks + H2AC
C           -------------------------------
C
            CALL X2GU(WORK(KHMAT),WORK(KH2XY),WORK(KH2YX),
     &           H2AC,IUU,IXU .EQ. IYU, DOFQ, DOH2AC)
C           CALL X2GU(HX,H2XY,H2YX,H2AC,IUU,ONLYXY,DOFQ,DOH2AC)
C
#ifdef MOD_DEBUG
            do i = 1, 3 * nz
               write(6,*) 'xy matrix no ',i
               call output(work(kh2xy + (i-1)*norbt*nasht),
     &              1,norbt,1,nasht,norbt,nasht,1,lupri)
            end do
            if ( ix .ne. iy) then
               do i = 1, 3 * nz
                  write(6,*) 'yx matrix no ',i
                  call output(work(kh2yx + (i-1)*norbt*nasht),
     &                 1,norbt,1,nasht,norbt,nasht,1,lupri)
               end do
            end if
#endif
C
C           Do the gazillion matrix multiplies to construct
C           the auxiliary Fock matrix FQ.
C
            IF (DOFQ) THEN
               CALL RADDFQ(FQ,WORK(KH2XY),WORK(KH2YX),PV,
     &              IUU,IXU,IYU,IREPIJ,IPRINT)
            END IF
C
 20      CONTINUE
 10   CONTINUE
C
C     Close integral files
C     --------------------
C
      DO IFSYM = 1, NFSYM
         CLOSE (LUINT + IFSYM, STATUS = 'KEEP')
      END DO
C
      IF (IPRINT .GE. 10) THEN
         IF (DOH2AC) THEN
            WRITE(LUPRI,'(A)') ' (RGETH21) H2AC matrix'
            CALL PRDNZ3(H2AC,NASHT,NNASHX,NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
         IF (DOFQ) THEN
            WRITE(LUPRI,'(A)') ' (RGETH21) FQ matrix'
            CALL PRQMAT(FQ,NORBT,NASHT,NORBT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
      END IF
C
      CALL QEXIT('RGETH21')
      RETURN
C
C     Read errors
C
 8000 CONTINUE
      WRITE(LUPRI,'(/A)')
     &     '*** ERROR in RGETH2 *** Read error'
 8100 CONTINUE
      WRITE(LUPRI,'(1X,A,I4)')
     &     'IREPIJ = ',IREPIJ,
     &     'IREC   = ',IREC,
     &     'IOSTAT = ',IOSVAL,
     &     'I      = ',I,
     &     'IX     = ',IX,
     &     'IY     = ',IY,
     &     'NXY    = ',NXY
      CALL QUIT('*** ERROR in RGTETH2 *** Error reading integral file')
      END
*
      SUBROUTINE IMPORT_ARDUCCA_INTEGRALS(T,IMPORT_ARDUCCA)
*
* Import sorted integrals from the ARDUCCA code
*
* Made automated restart from LU_ARDUCCA (KRCC_INTEGRALS) if it is there
* Also checks for formatted or unformatted file
*
* Lasse,2013 
      IMPLICIT REAL*8(A-H,O-Z)
#include "clunit.inc"
*. Integrals to be defined
      DIMENSION T(*)
      LOGICAL YN
      CHARACTER(LEN=2) UNFORM
*
* Automated restart if LU_ARDUCCA (KRCC_INTEGRALS) exits
      INQUIRE(FILE='KRCC_INTEGRALS',EXIST=YN)
      IF(YN) THEN
        INQUIRE(FILE='KRCC_INTEGRALS',UNFORMATTED=UNFORM)
        IF(UNFORM(1:2).EQ.'UN') THEN
          I_FORMATTED = 0
        ELSE
          I_FORMATTED = 1
        END IF
* Now we know that we improt integrals
        IMPORT_ARDUCCA = 1
      ELSE
        IMPORT_ARDUCCA = 0
      END IF
          I_FORMATTED = 1
      print*,'I_FORMATTED',I_FORMATTED
*
      IF(IMPORT_ARDUCCA.EQ.1) THEN
        WRITE(6,*) ' Reading in integrals from ', LU_ARDUCCA
        IF(I_FORMATTED.EQ.1) THEN
*. Formatted
          OPEN(Unit=LU_ARDUCCA,File='KRCC_INTEGRALS',
     &     Status='UNKNOWN',Form='FORMATTED')
          CALL REWINO(LU_ARDUCCA)
          READ(LU_ARDUCCA,*) NINT
          READ(LU_ARDUCCA,*) (T(I),I=1,NINT)
C         DO I = 1, NINT
C           READ(LU_ARDUCCA,*) T(I)
C         END DO
        ELSE
*. Unformatted
          OPEN(Unit=LU_ARDUCCA,File='KRCC_INTEGRALS',
     &     Status='UNKNOWN',Form='UNFORMATTED',IOSTAT=IOS)
        IF(IOS.NE.0) THEN
          call abend2('problem with file krcc_integrals')
        END IF
          CALL REWINO(LU_ARDUCCA)
          READ(LU_ARDUCCA) NINT
          READ(LU_ARDUCCA) (T(I),I=1,NINT)
C         print*,'NINT',NINT
C         DO I = 1, NINT
C           READ(LU_ARDUCCA) T(I)
C           print*,'T,I',T(I),I
C         END DO
        END IF
        CLOSE(LU_ARDUCCA)
      END IF
*
      IF(IMPORT_ARDUCCA.EQ.1) THEN
        WRITE(6,*) ' A set of integrals have been found '
        WRITE(6,*) ' This is a restarted calculation '
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Initial set of integrals '
        WRITE(6,*) ' ======================== '
        WRITE(6,*)
        WRITE(6,*) ' Number of integrals read ',NINT
        WRITE(6,*)
        DO I=1,NINT
          WRITE(6,*) T(I)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE DIFFF(VEC1,VEC2,NDIM)
*
#include "implicit.inc"
      DIMENSION VEC1(NDIM), VEC2(NDIM)
*
      ITOOMUCH = 0
      AERROR = 1.0E-10
      DO I=1,NDIM
        DIFF = VEC1(I)-VEC2(I)
        print*,'DIFFFF',DIFF,VEC1(I),VEC2(I)
        IF(ABS(DIFF).GE.AERROR) THEN
         ITOOMUCH = ITOOMUCH + 1
        END IF
      END DO
      IF(ITOOMUCH.GE.1) THEN
       print*,'SOMETHING IS WRONG',ITOOMUCH
       stop 'loll'
      END IF
      END 
