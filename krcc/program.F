C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PAMKRCC(WORK,LWORK)
C***********************************************************************
C
C Getting information from Dirac and setting up strings etc.
C No calculation of offset. Will just use dirac wrk or allocate
C
C Written by Lasse Soerensen 2010
C
C
C common blocks needs cleaning!!!!!!!!!!!!!!!!
C
C***********************************************************************
#include "implicit.inc"
*. Parameters for dimensioning
#include "mxpdim.inc"
C will probably have to delete most of these again
*.File numbers
#include "clunit.inc"
*.Print flags
#include "cprnt.inc"
#include "lucinp.inc"
#include "cstate.inc"
#include "crun.inc"
#include "cicisp.inc"
#include "oper.inc"
#include "cgas.inc"
C#include "chktyp.inc"
#include "glbbas.inc"
C#include "wrkspc.inc"
C added by lasse
#include "cintfo.inc"
C#include "iratdef.h"
#include "ipoist8.inc"
C
C#include "memint.h"
      DIMENSION WORK(LWORK)
C
C
#include "memint.h"
C
      CALL QENTER('KRCC ')
*. Ratio beteeen real and integer word length
      CALL ZIRAT_KRCC !ok
* First make an entry on memory to look for
      IONE = 1
      CALL MEMGET('INTE',KLOOK,IONE,WORK,KFREE,LFREE)
* First do transformation and get integrals from DIRAC
      CALL READ_DIRAC_FOR_KRCC(WORK,KFREE,LFREE)
*. Assign diskunits 
!!!!needs updating!!!!
!!!!needs updating!!!!
!!!!needs updating!!!!
!!!!needs updating!!!!
      CALL DISKUN_KRCC
*. Header 
      CALL WRTHD_KRCC(LUOUT) !ok
      CALL MEMCHK('LL1',WORK,KLOOK)
*. Ratio beteeen real and integer word length
C      CALL ZIRAT_KRCC !ok
*. From shells to orbitals
      CALL ORBINF_KRCC(LUOUT,IPRORB)
*. Types of determinants to be included
      CALL DETTYP_KRCC(IPRSTR)
*. Number of string types
      CALL MEMCHK('LL2',WORK,KLOOK)
      CALL STRTYP_GAS_KRCC(IPRSTR)
*. Symmetry information
      CALL SYMINF_KRCC(IPRORB)
*. Number of integrals
      CALL INTDIM_KRCC(IPRTRA)
*. Allocate local and global arrays 
!!!!use allocate here!!!!
!!!!needs updating!!!!
!!!!needs updating!!!!
!!!!needs updating!!!!
!!!!needs updating!!!!
      CALL ALLOC_KRCC(WORK,KFREE,LFREE)
      CALL MEMCHK('LULU6',WORK,KLOOK)
*. Internal string information (stored in WORK, bases in /STRBAS/)
      CALL STRINF_GAS_KRCC(WORK,KFREE,LFREE,IPRSTR)
      CALL MEMCHK('LULU5',WORK,KLOOK)
*. Internal subspaces
      CALL LCISPC_KRCC(WORK,KFREE,LFREE,LCSBLK,IPRCIX)
      CALL MEMCHK('LULU1',WORK,KLOOK)
*. Divide orbital spaces into inactive/active/secondary!!!!not sure to include!!!!
C      CALL GASSPC_KRCC
* Find occupation of GAS
C Notice JCMBSPC is hardwired
      JCMBSPC = 1
      CALL GET_GASOCC_KRCC(JCMBSPC)
      CALL MEMCHK('LULU3',WORK,KLOOK)
* Get Hamiltonian strings. 
* Notice a similar routine is called later which does the actual setup for the cc runs. 
* This is just for the integral sorting
C      CALL SET_HOP_DBG_KRCC(4,1,WORK,LWORK)
C assumes totally symmetric hamiltonian
!!!!!May need updating to elminate obsolete operators!!!!!
!!!!!May need updating to elminate obsolete operators!!!!!
!!!!!May need updating to elminate obsolete operators!!!!!
!!!!!May need updating to elminate obsolete operators!!!!!
      ISYM_T = 1
      CALL SET_HOP_DBG_FOR_NEWCCV(4,1,WORK,KFREE,LFREE)
C     CALL SET_HOP_DBG_KRCC(4,1,WORK,KFREE,LFREE)
      CALL MEMCHK('LULU2',WORK,KLOOK)
C      CALL MEMCHK_KRCC(WORK)
      print*,'after SET_HOP_DBG_KRCC'
*. Read in integrals !!!!Notice non-matching call!!!!
C       call picasso(INTSPACE,ISO1,IPRTRA,
C    &               NL2D,IUB2D,ISM2D)
      CALL INTIM_KRCC(INTSPACE,MCSCF,WORK,KFREE,LFREE)
      CALL MEMCHK('LUL22',WORK,KLOOK)
*
*******************************************************
C
C Here we need to allocate some vectors and other stuff
C
C We will only assume that we can do CC calculations
C
*******************************************************
*
      CALL KRCC_GENCC(IREFSM,ITSPC,JCMBSPC,IPRDIA,
     &                I_TRANS_WF,II_RES_EXC,
     &                II_DO_CMPCCI,ECORE_KRMC,
     &                WORK,KFREE,LFREE)


*
      CALL QEXIT('KRCC ')
*
      RETURN
      END
*
       SUBROUTINE KRCC_GENCC(ISM,ISPC,ICISPC,IPRNT,
     &                       I_TRANS_WF,II_RES_EXC,II_DO_CMPCCI,
     &                       ECORE_KRMC,
     &                       WORK,KFREE,LFREE)
*
* Master routine for general coupled cluster calc with LUCIA
* - not any more
*      Lasse
*
* Now master routine for general relativistic mr coupled cluster
*
*
#include "implicit.inc"
      REAL*8 INPRDD
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "strinp.inc"
#include "cgas.inc"
#include "gasstr.inc"
C#include "wrkspc.inc"
#include "orbinp.inc"
#include "glbbas.inc"
#include "cands.inc"
#include "crun.inc"
#include "clunit.inc"
#include "cc_exc.inc"
#include "cecore.inc"
#include "cstate.inc"
#include "cprnt.inc"
#include "csm.inc"
#include "ctcc.inc"
#include "oper.inc"
#include "cintfo.inc"
*
* for real or complex groups
#include "../include/dgroup.h"
      COMMON/CMXCJ/MXCJ
      CHARACTER*6 CCTYPE
      DIMENSION WORK(*)
*
*************************************************************
*                                                           *
* Things not yet appliable (or never will be?)              *
*                                                           *
* I_TRANS_WF, II_DO_CMPCCI, ICCSOLVE, II_DO_CMPCCI,         *
* II_RES_EXC, MXSPOX                          *
*                                                           *
*                                                 Lasse     *
*                                                           *
*                                                           *
*************************************************************
*CALL MEMGET('INTE',KLIHIND,4*NSPOBEX_TP,WORK,KFREE,LFREE)
*CALL MEMGET('REAL',KSIGNNHX,NSPOBEX_TP,WORK,KFREE,LFREE)
*
* Hardwirerings
***************************
* Hardwirering IHTYPE for now Lasse used in get_hx_rela
      IHTYPE = 4
* NOTICE HARDWIRERING OF ISPC SINCE ONLY ONE SPACE WILL NOW BE AS INPUT
      ISPC = ICISPC
* def
      JCMBSPC=ICISPC
      I_DO_CC = 1
      ECORE = ECORE_KRMC
      ISIMSYM = 0
*
* Check if real or complex/quarternion algebra is used
      IF(NZ.EQ.1) THEN
      WRITE(6,*) ' Using real algebra '
      IRECOM = 1
      IMULTFAC = 1
      WRITE(6,*) ' Using diagonal + subspace jacobian'
      I_DO_SBSPJA = 1
      ELSE
      WRITE(6,*) ' Using complex algebra '
      IRECOM = 2
      IMULTFAC = 2
      WRITE(6,*) ' Using diagonal approximation to Jacobian'
      I_DO_SBSPJA = 0
      WRITE(6,*) ' Complex algebra still only possible with old CC'
C     IF(I_DO_NEWCCV.NE.0) THEN
        WRITE(6,*) ' Complex algebra still only possible with old CC '
        WRITE(6,*) ' Can at the moment only be used with CCONCI '
        CALL ABEND()
C     END IF
      END IF
*. General CC so
      CCTYPE(1:6) = 'GEN_CC'
* ADDED BY LASSE SINCE ONLY FORM HERE ALLOWED
      CCFORM(1:3) = 'TCC'
*
      WRITE(6,*)
      WRITE(6,*) ' ========================================'
      WRITE(6,*) ' General Coupled Cluster section entered '
      WRITE(6,*) ' ========================================'
      WRITE(6,*)
      WRITE(6,*) ' Reference space ',ICISPC
      WRITE(6,*) ' Space defining T-operators ', ISPC
      WRITE(6,*) ' I_TRANS_WF  = ',  I_TRANS_WF
      WRITE(6,*) ' II_DO_CMPCCI = ', II_DO_CMPCCI
*
      NTEST = 0
      IPRNT = MAX(IPRNT,NTEST)
* Taken print level from CI. Lasse
      IPRCC = IPRNT
*
      I12 = 2
      IDUM = 0
*
* Routine added for integral fetching will just set up an array needed for NEWCCV
* Always on now
C     IF(I_DO_NEWCCV.EQ.1) THEN
        CALL SUM_ARRAY
C     END IF
* Dimensioning
* ADDED FOR REL CALC
C     CALL Z_BLKFO_KRCC(IDC,NMS2VAL,
C    &                 ISPC,ISM,KPCLBT,KPCLEBT,KPCI1BT,KPCIBT,KPCBLTP,
C    &                 NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,NBAT_MS2,
C    &                 IBAT_MS2,
C    &                 WORK,KFREE,LFREE)
      CALL MEMCHK_KRCC(WORK)
*
*
*. Space for BMATRIX and scratch for DIIS
C For space savings Lasse
C Hardwirering of optimizer
      ICCSOLVE = 1
      IF(ICCSOLVE.EQ.2) THEN
*
        MXDIM = MAX(MAXIT,MXCIV_ORIG)+1
        LENB = MXDIM*(MXDIM+1)/2
        LENSCR = 3*(MXDIM+1)**2 + MXDIM + 1
        CALL MEMMAR(KBMAT,LENB,'ADDL  ',2,'BM_DIS')
        CALL MEMMAR(KSCDIIS,LENSCR,'ADDL  ',2,'SC_DIS')
      ELSE
        KBMAT = 0
        KSDIIS = 0
      END IF
*
      IF(I_DO_SBSPJA.EQ.1) THEN
        MXDIM = MAX(MAXIT,MXCIV_ORIG)
*. Space for working with Subpspace Jacobian :
        LSBSPJA = 5*MXDIM**2 + 2*MXDIM
      ELSE
        LSBSPJA = 0
      END IF
      CALL MEMGET('REAL',KSBSPJA,LSBSPJA,WORK,KFREE,LFREE)
      CALL MEMCHK_KRCC(WORK)
*
*
*. Will now choose a reference determinant based on input!
*. All possible mk2 projections possible. Lasse
*
      CALL REF_DET_KRCC
*
* End the addition of determinant
*
      ITSPC  = ISPC
      IETSPC = ICISPC
*. The reference space is the first space
      IREFSPC = 1
*
* Find the type of reference state
* ================================
*
*. Divide orbitals into HOLE ( Double occupied in reference space,
*  PARTICLES( unoccupied in reference) and Valence orbitals
*  (not complete occupied or completely unoccupied)
*
*   IREFTYP = 1 => CLOSED Shell HARTREE-FOCK,
*   IREFTYP = 2 => High spin open shell single det state
*   IREFTYP = 3 => Cas state or more general multireference state
      CALL CC_AC_SPACES_KRCC(IREFSPC,IREFTYP)
      WRITE(6,*) ' IREFTYP after call to CC_AC_SPACES ', IREFTYP
*. Number of active orbital spaces
      NACT_SPC = 0
      IACT_SPC = 0
*. Number of hole spaces
      NHOL_SPC = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.3) THEN
          NACT_SPC = NACT_SPC + 1
          IACT_SPC = IGAS
        END IF
        IF(IHPVGAS(IGAS).EQ.1) THEN
          NHOL_SPC = NHOL_SPC + 1
        END IF
      END DO
      IF(IPRNT.GE.5)
     &WRITE(6,*) ' Number of hole spaces = ', NHOL_SPC
*
      IF(NACT_SPC.GT.1) THEN
        WRITE(6,*) ' GENCC  in problems '
        WRITE(6,*) ' More than one active orbital spaces '
        WRITE(6,*) ' NACT_SPC = ',  NACT_SPC
        CALL ABEND2(' GENCC :  More than one active orbital spaces ')
      END IF
      IF(IPRNT.GE.5)
     &WRITE(6,*) ' GENCC : IACT_SPC,NACT_SPC',IACT_SPC,NACT_SPC
*
*. Info on active-active excitation types
      CALL ACAC_EXC_TYP_KRCC(IAAEXC_TYP,MX_AAEXC,IPRCC)
C     print*,'IAAEXC_TYP,MX_AAEXC,IPRCC',IAAEXC_TYP,MX_AAEXC,IPRCC
      CALL MEMCHK_KRCC(WORK)
*. pt active rotations are eliminated for CAS
      IF(IPRNT.GE.5)
     &WRITE(6,*) ' GENCC : MX_AAEXC ', MX_AAEXC
*
* 1 : Find the set of T amplitudes.
*     These are defined as the set of excitation operators needed
*     to obtain the CI space ISPC by exciting determinants in
*     the reference space
*
*. Orbital excitation operators
*. Number of types of orbital excitation operators and length to store
*  these types
*
*. Number of orbital excitations equals the number of occupation
* classes in ISPC ( minus the reference space, which should be
* included )
*
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NEL = NAEL + NBEL
*
      ICSPC = ICISPC
      ISSPC = ICISPC
*. Routines use complete TT blocks so
C     ISIMSYM = 1
*
* ========================
* info for reference space
* THIS NEEDS CHANGING!!!
* Could bouild a dummy like reference on ref_det
* Assumes from above a single reference occupation space
* ========================
*
      NOCCLS_REF = 1
      CALL MEMGET('INTE',KLOCCLS_REF,NGAS,WORK,KFREE,LFREE)
      CALL OCCLSE_REF(WORK(KLOCCLS_REF))
*
*. Make sure that there is just a single occupation space
C     CALL OCCLSE_KRCC(1,NOCCLS_REF,IOCCLS,NEL,IREFSPC,0,0,NOBPT)
C     IF(NOCCLS_REF.NE.1) THEN
C       WRITE(6,*) ' Problem in general CC '
C       WRITE(6,*)
C    &  ' Reference space is not a single occupation space'
C       STOP
C    &  ' Reference space is not a single occupation space'
C     END IF
*. and the reference occupation space
C     CALL MEMMAR(KLOCCLS_REF,NGAS,'ADDL  ',1,'OCC_RF')
C     CALL MEMGET('INTE',KLOCCLS_REF,NGAS,WORK,KFREE,LFREE)
C     CALL OCCLSE_KRCC(2,NOCCLS_REF,WORK(KLOCCLS_REF),NEL,IREFSPC,0,0,
C    &                NOBPT)
*
* ====================================
* Info for space defining excitations
* ====================================
*
*. Number
      CALL OCCLSE_KRCC(1,NOCCLS,IOCCLS,NEL,ISPC,0,0,NOBPT)
*. And the occupation classes
      CALL MEMGET('INTE',KLOCCLS,NOCCLS*NGAS,WORK,KFREE,LFREE)
      CALL OCCLSE_KRCC(2,NOCCLS,WORK(KLOCCLS),NEL,ISPC,0,0,NOBPT)
*. Number of occupation classes for T-operators
      NTOCCLS = NOCCLS
*. It could be an idea to check that reference space is included
*
* ========================
* Orbital excitation types
* ========================
*
*. Number of excitation types
      IFLAG = 1
      IDUM = 1
      CALL MEMCHK_KRCC(WORK)
      CALL TP_OBEX2_KRCC(NOCCLS,NEL,NGAS,WORK(IDUM),
     &             WORK(IDUM),WORK(IDUM),
     &             WORK(KLOCCLS),WORK(KLOCCLS_REF),MX_NCREA,MX_NANNI,
     &             MX_EXC_LEVEL,WORK(IDUM),MX_AAEXC,IFLAG,
     &             NOBEX_TP,IPRCC)
*
      IF(IPRNT.GE.5)
     &WRITE(6,*) ' NOBEX_TP,MX_EXC_LEVEL = ', NOBEX_TP,MX_EXC_LEVEL
*. And the actual orbital excitations
*.  An orbital excition operator is defined by
*   1 : Number of creation operators
*   2 : Number of annihilation operators
*   3 : The actual creation and annihilation operators
*. The number of orbital excitations is increased by one to include
*. excitations within the reference space
      NOBEX_TPE = NOBEX_TP+1
      CALL MEMGET('INTE',KLCOBEX_TP,NOBEX_TPE,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLAOBEX_TP,NOBEX_TPE,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KOBEX_TP,NOBEX_TPE*2*NGAS,WORK,KFREE,LFREE)
*. Excitation type => Original occupation class
      CALL MEMGET('INTE',KEX_TO_OC,NOBEX_TPE,WORK,KFREE,LFREE)
      CALL MEMCHK_KRCC(WORK)
      IFLAG = 0
      CALL TP_OBEX2_KRCC(NOCCLS,NEL,NGAS,WORK(KOBEX_TP),
     &             WORK(KLCOBEX_TP),WORK(KLAOBEX_TP),
     &             WORK(KLOCCLS),WORK(KLOCCLS_REF),MX_NCREA,MX_NANNI,
     &             MX_EXC_LEVEL,WORK(KEX_TO_OC),MX_AAEXC,IFLAG,
     &             NOBEX_TP,IPRCC)
      CALL MEMCHK_KRCC(WORK)
*
* =======================
* Spinorbital excitations
* =======================
*
*. Spin combinations of CC excitations : Currently we assume that
*. The T-operator is a singlet, can 'easily' be changed. Done Lasse
      WRITE(6,*) ' MSCOMB_CC in LUCIA_GENCC ', MSCOMB_CC
*
*. Notice : The first time in OBEX_TO_SPOBEX we always use MSCOMB_CC = 0.
*. This may lead to the allocation of too much space for
*. spinorbital excitations, but MSCOMB_CC = 1, requires access
*. to WORK(KLSOBEX) which has not been defined
*
*. Largest spin-orbital excitation level
*
* MXSPOX set from outside and not in use
*
* NOTICE MXSPOX is here hardwired (it is the largest allowed excitation)
      MXSPOX = 0
      IF(MXSPOX.NE.0) THEN
        MXSPOX_L = MXSPOX
      ELSE
        MXSPOX_L = MX_EXC_LEVEL
      END IF
C     WRITE(6,*) ' MXSPOX, MXSPOX_L, MX_EXC_LEVEL = ',
C    &             MXSPOX, MXSPOX_L, MX_EXC_LEVEL
      IZERO = 0
*
* Notice MS2TOT == 0 set from this side NOBPT,0,IZERO (it is the 0)
* It is taken from common block on inside
      CALL OBEX_TO_SPOBEX_KRCC(1,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,IDUMMY,NSPOBEX_TP_CC,NGAS,
     &     NOBPT,0,IZERO ,IAAEXC_TYP,IACT_SPC,IPRCC,IDUMMY,
     &     MXSPOX_L,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),NAEL,NBEL,IREFSPC,
     &     IRECOM,WORK,KFREE,LFREE)
*. Extended number of spin-orbital excitations : Include
*. unit operator as last spinorbital excitation operator
C allowing ample space for all operator that is yet to be found
C Notice OBEX_TO_SPOBEX_REL has been changed significantly inside. Lasse
      NSPOBEX_TPE = NSPOBEX_TP_CC + 1
      IF(IPRNT.GE.5) WRITE(6,*) ' NSPOBEX_TP_CC = ', NSPOBEX_TP_CC
C      print*,'NSPOBEX_TP_CC,NSPOBEX_TPE',NSPOBEX_TP_CC,NSPOBEX_TPE
*. And the actual spinorbital excitation operators
      CALL MEMGET('INTE',KLSOBEX_CC,4*NGAS*NSPOBEX_TPE,WORK,KFREE,LFREE)
*. Map spin-orbital exc type => orbital exc type
      CALL MEMGET('INTE',KLSOX_TO_OX,NSPOBEX_TPE,WORK,KFREE,LFREE)
*. First SOX of given OX ( including zero operator )
      CALL MEMGET('INTE',KIBSOX_FOR_OX,NOBEX_TP+1,WORK,KFREE,LFREE)
*. Number of SOX's for given OX
      CALL MEMGET('INTE',KNSOX_FOR_OX,NOBEX_TP+1,WORK,KFREE,LFREE)
*. SOX for given OX
      CALL MEMGET('INTE',KISOX_FOR_OX,NSPOBEX_TPE,WORK,KFREE,LFREE)
      CALL OBEX_TO_SPOBEX_KRCC(2,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,WORK(KLSOBEX_CC),NSPOBEX_TP_CC,
     &     NGAS,NOBPT,0,MSCOMB_CC,IAAEXC_TYP,IACT_SPC,
     &     IPRCC,WORK(KLSOX_TO_OX),MXSPOX_L,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),NAEL,NBEL,IREFSPC,
     &     IRECOM,WORK,KFREE,LFREE)
      NSPOBEX_TPE = NSPOBEX_TP_CC + 1
      CALL MEMCHK_KRCC(WORK)
*. Add unit-operator as last spinorbital excitation
      IZERO = 0
      CALL ISTVC3(WORK(KLSOBEX_CC),(NSPOBEX_TPE-1)*4*NGAS+1,IZERO,
     &            4*NGAS)
      IF(IPRNT.GE.5) THEN
        WRITE(6,*) ' Extended list of spin-orbital excitations : '
        CALL WRT_SPOX_TP_CC_KRCC(WORK(KLSOBEX_CC),NSPOBEX_TPE)
      END IF
      CALL ISTVC3(WORK(KLSOX_TO_OX),NSPOBEX_TPE,NOBEX_TP+1,1)
*. Mapping spinorbital excitations => occupation classes
      CALL MEMGET('INTE',KIBSOX_FOR_OCCLS,NOCCLS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNSOX_FOR_OCCLS,NOCCLS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KISOX_FOR_OCCLS,NSPOBEX_TPE,WORK,KFREE,LFREE)
      CALL SPOBEX_FOR_OCCLS(WORK(KEX_TO_OC),NOCCLS,WORK(KLSOX_TO_OX),
     &     NSPOBEX_TPE,WORK(KNSOX_FOR_OCCLS),WORK(KISOX_FOR_OCCLS),
     &     WORK(KIBSOX_FOR_OCCLS))
      CALL MEMCHK_KRCC(WORK)
*
*. Frozen spin-orbital excitation types
*******
* setting NFRZ_CC = 0 Lasse, FRZ_SPOBEX not yet appliable
*******
C     NFRZ_CC = 0
C     CALL FRZ_SPOBEX(WORK(KLSPOBEX_FRZ),WORK(KLCOBEX_TP),NSPOBEX_TP_CC,
C    &                WORK(KLSOX_TO_OX),IFRZ_CC_AR,NFRZ_CC)
C     IZERO = 0
C     CALL ISTVC3(WORK(KLSPOBEX_FRZ),NSPOBEX_TP_CC+1,IZERO,1)
*. Spin-orbital excitation types related by spin-flip
      CALL MEMGET('INTE',KLSPOBEX_AB,2*NSPOBEX_TP_CC+1,WORK,KFREE,LFREE)
C NEXT ARRAY ONLY ALLOCATED FOR THE PROGAM TO RUN LASSE
      CALL MEMGET('INTE',KLSPOBEX22,2*NSPOBEX_TP_CC+1,WORK,KFREE,LFREE)
      CALL SPOBEXTP_PAIRS_KRCC(NSPOBEX_TP_CC+1,WORK(KLSOBEX_CC),NGAS,
     &                    WORK(KLSPOBEX_AB))
      IF(NTEST.GE.5) WRITE(6,*) ' After SPOBEXTP_PAIRS'
*************
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Extended list of spin-orbital excitations : '
        CALL WRT_SPOX_TP_CC_KRCC(WORK(KLSOBEX_CC),2*NSPOBEX_TPE)
      END IF
**************
*. Alpha- and beta-excitations constituting the spinorbital excitations
*. Number
      CALL MEMCHK_KRCC(WORK)
      CALL SPOBEX_TO_ABOBEX_CC(WORK(KLSOBEX_CC),NSPOBEX_TP_CC,NGAS,
     &     1,NAOBEX_TP,NBOBEX_TP,IDUMMY,IDUMMY)
*. And the alpha-and beta-excitations
      LENA = 2*NGAS*NAOBEX_TP
      LENB = 2*NGAS*NBOBEX_TP
      CALL MEMCHK_KRCC(WORK)
      print*,'LENA,LENB,KFREE,LFREE',LENA,LENB,KFREE,LFREE
      CALL MEMCHK_KRCC(WORK)
      CALL MEMGET('INTE',KLAOBEX_CC,LENA,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLBOBEX_CC,LENB,WORK,KFREE,LFREE)
      CALL SPOBEX_TO_ABOBEX_CC(WORK(KLSOBEX_CC),NSPOBEX_TP_CC,NGAS,
     &     0,NAOBEX_TP,NBOBEX_TP,WORK(KLAOBEX_CC),WORK(KLBOBEX_CC))
      IF(NTEST.GE.10) WRITE(6,*) ' After SPOBEX_TO_ABOBEX_CC '
*. Max dimensions of CCOP !KSTR> = !ISTR> maps
*. For alpha excitations
      IATP = 1
      IOCTPA = IBSPGPFTP(IATP)
      NOCTPA = NSPGPFTP(IATP)
*
      CALL LEN_GENOP_STR_MAP_KRCC(
     &     NAOBEX_TP,WORK(KLAOBEX_CC),NOCTPA,NELFSPGP(1,IOCTPA),
     &     NOBPT,NGAS,MAXLENA)
      IBTP = 2
      IOCTPB = IBSPGPFTP(IBTP)
      NOCTPB = NSPGPFTP(IBTP)
      CALL LEN_GENOP_STR_MAP_KRCC(
     &     NBOBEX_TP,WORK(KLBOBEX_CC),NOCTPB,NELFSPGP(1,IOCTPB),
     &     NOBPT,NGAS,MAXLENB)
      IF(NTEST.GE.10) WRITE(6,*) ' After LEN_GENOP_STR_MAP_CC '
      MAXLEN_I1 = MAX(MAXLENA,MAXLENB)
      IF(IPRNT.GE.5) WRITE(6,*) ' MAXLEN_I1 = ', MAXLEN_I1
*
* Max Dimension of spinorbital excitation operators
*
      CALL MEMGET('INTE',KLLSOBEX_CC,NSPOBEX_TPE,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLIBSOBEX_CC,NSPOBEX_TPE,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLSPOBEX_AC_CC,NSPOBEX_TPE,WORK,KFREE,LFREE)
*. KLSPOBEX_AC allocated for later and transfered via common blocks
*. ALl spinorbital excitations are initially active
*
            IONE = 1
            CALL ISETVC(WORK(KLSPOBEX_AC_CC),IONE,NSPOBEX_TPE)
*
CLasse      MX_ST_TSOSO_MX = 0
CLasse      MX_ST_TSOSO_BLK_MX = 0
CLasse      MX_TBLK_MX = 0
CLasse      MX_TBLK_AS_MX = 0
CLasse      LEN_T_VEC_MX = 0
*
C      print*,'IFTTWO,ILTTWO ',IFTTWO,ILTTWO
CLasse    DO ICCAMP_SM = 1, NSMST
CLasse        print*,'ICCAMP_SM ',ICCAMP_SM
*
COri        CALL IDIM_TCC_DBG(WORK(KLSOBEX),NSPOBEX_TP,ICCAMP_SM,
COri     &       MX_ST_TSOSOL,MX_ST_TSOSO_BLKL,MX_TBLKL,
COri     &       WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VECL,
COri     &       MSCOMB_CC,MX_TBLK_AS,
COri     &       WORK(KISOX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
COri     &       NTCONF,IPRCC)
CLasse        CALL IDIM_TCC_DBG(WORK(KLSOBEX),NSPOBEX_TP,ICCAMP_SM,
CLasse     &       MX_ST_TSOSOL,MX_ST_TSOSO_BLKL,MX_TBLKL,
CLasse     &       WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VECL,
CLasse     &       MSCOMB_CC,MX_TBLK_AS,
CLasse     &       IFTONE,IFTTWO,ILTTWO,N1ELINT,N2ELINT,
CLasse     &       IPRCC)
*
CLasse        MX_ST_TSOSO_MX = MAX(MX_ST_TSOSO_MX,MX_ST_TSOSOL)
CLasse        MX_ST_TSOSO_BLK_MX = MAX(MX_ST_TSOSO_BLK_MX,MX_ST_TSOSO_BLKL)
CLasse        MX_TBLK_MX = MAX(MX_TBLK_MX,MX_TBLKL)
CLasse        MX_TBLK_AS_MX = MAX(MX_TBLK_AS_MX,MX_TBLK_AS)
CLasse        LEN_T_VEC_MX = MAX(LEN_T_VEC_MX, LEN_T_VECL)
C        print*,'N1ELINT,N2ELINT ',N1ELINT,N2ELINT
*
CLasse      END DO
      IF(IPRNT.GE.5) WRITE(6,*) ' MX_TBLK_AS_MX = ', MX_TBLK_AS_MX
      IF(IPRNT.GE.5) WRITE(6,*) ' LEN_T_VEC_MX = ', LEN_T_VEC_MX
C      stop 'after IDIM'
*
*
* Cluster operator setup is now finshed
* Find the number of holes and particles to be contracted for a given H
* operator
*
      CALL HTYPE_PH_CONTRACTIONS_MASTER(WORK,KFREE,LFREE)
*
* Start setting up intermediates
*
      CALL SET_INTERMEDIATES_MASTER(WORK,KFREE,LFREE)
*
* Sort Cluster operator and intermediates according to 
* N, delta Mk, Mub
* for easy solution of equations
*
      CALL SORT_OP_MASTER(WORK,KFREE,LFREE)
*
* Not so pretty temp solution (need to dimension T after sort!)
*
*. And dimensions for symmetry 1
      ITOP_SM = 1
*
*
* NOT SURE IDIM_TCC_KRCC WORKS AS IT SHOULD!!!! LASSE
* Should be moved since T is resorted!!!!! Lasse 2011 dont ignore
C     CALL IDIM_TCC_KRCC(WORK(KLSOBEX_CC),NSPOBEX_TP_CC,ITOP_SM,
      CALL IDIM_TCC_KRCC(WORK(KLSOBEX_CC),NSPOBEX_TPE,ITOP_SM,
     &     MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &     WORK(KLLSOBEX_CC),WORK(KLIBSOBEX_CC),LEN_T_VEC_CC,
     &     MSCOMB_CC,MX_SBSTR,
     &     IFTONE,IFTTWO,-1,N1ELINT,N2ELINT,
     &     IPRCC)
*
      IF(NTEST.GE.100) WRITE(6,*) ' After IDIM_TCC section '
      MX_TBLK_MX = MX_TBLK
      MX_ST_TSOSO_MX =  MX_ST_TSOSO
C temp change in size
      MX_ST_TSOSO_BLK_MX = 2*MX_ST_TSOSO_BLK
      N_CC_AMP = LEN_T_VEC_CC
*. Hard wire info for unit operator stored as last spinorbital excitation
C  ISTVC2(IVEC,IBASE,IFACT,NDIM)
      IONE = 1
C
      CALL ISTVC3(WORK(KLLSOBEX_CC),NSPOBEX_TPE,IONE,1)
      N_CC_AMPP1 = N_CC_AMP + 1
*. Why is this done???
C     CALL ISTVC3(WORK(KLIBSOBEX_CC),NSPOBEX_TPE,N_CC_AMPP1,1)
      IF(IPRNT.GE.5)
     &WRITE(6,*) ' LUCIA_GENCC : N_CC_AMP = ', N_CC_AMP
*
*. Pointers to KVEC1 and KVEC2, transferred through GLBBAS
      CALL GET_3BLKS_GCC_REL(KVEC1,KVEC2,KVEC3,MXCJ,
     &                       WORK,KFREE,LFREE)
       KVEC1P = KVEC1
       KVEC2P = KVEC2
       print*,'KVEC1,KVEC2,KVEC1P,KVEC2P',KVEC1,KVEC2,KVEC1P,KVEC2P
       WRITE(6,*) ' KVEC3 after GET_3BLKS.. ', KVEC3
*. and two CC vectors , extra element for unexcited SD at end of vectors
       N_SD_INT = 1
* next line added by Lasse not sure if correct but thinks so
       LEN_T_VEC_MX = LEN_T_VEC_CC
*
       LENNY = LEN_T_VEC_MX + N_SD_INT
       print*,'LENNY',LENNY
      CALL MEMGET('REAL',KCC1,IMULTFAC*LENNY,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCC2,IMULTFAC*LENNY,WORK,KFREE,LFREE)
*
* End not so pretty temp solution
*
*
* Now that operators are found and sorted it is time to solve the cc
* equations. This will be done by first finding all possible
* contractions and then eliminating the redundant ones and in the
* process calculate permutation factors and signs for contractions
*
      CALL CONTRACTION_SOLVER(WORK,KFREE,LFREE)
*
* Now to dimension and setup arrays for intermediates.
* This should later be moved to CC_VEC_FNC_KRCC to save space
*
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC)
*
* ============================================================
* Initialize : Restart, set coefs to zero, start from CI coefs
* ============================================================
*
C Now with automated restart. Lasse
      CALL INI_CC_AMP_KRCC(WORK(KCC1))
      CALL MEMCHK_KRCC(WORK)
*
* ==========================================
*. Transfer control to optimization routine
* ==========================================
*
*. New optimizer
C     CALL OPTIM_CC2_KRCC(WORK(KCC1),WORK(KCC2),WORK(KDIA),
C    &                   WORK(KVEC1),WORK(KVEC2),CCTYPE,ECC,
C    &                   WORK,KFREE,LFREE)
C     STOP 'CC ENERGY HAS BEEN OPTIMIZED'
*
* ===========================
* Calculate density matrices
* ===========================
*
      IDENSI = 0
      IF(IDENSI.NE.0) THEN
*. Three more vectors
        CALL MEMMAR(KCC3,LEN_T_VEC_MX,'ADDL  ',2,'CC3_VE')
        CALL MEMMAR(KCC4,LEN_T_VEC_MX,'ADDL  ',2,'CC4_VE')
        CALL MEMMAR(KCC5,LEN_T_VEC_MX,'ADDL  ',2,'CC5_VE')
C       CALL CC_DENSI(WORK(KCC1),WORK(KCC2),WORK(KDIA),
C    &       WORK(KCC3),WORK(KCC4),WORK(KCC5),
C    &       WORK(KVEC1),WORK(KVEC2) )
C     CC_DENSI(CCAMP,CCLAG,CCDIA,CCVEC1,CCVEC2,CCVEC3,
C    &           VEC1,VEC2)
      END IF
*
* =========================================
*. Coupled cluster first order properties
* =========================================
*
      NPROP = 0
      IF(NPROP.GT.0) THEN
*. Calculate properties
        KLDUM = 1
C       CALL ONE_EL_PROP(1,0,WORK(KLDUM))
      END IF
*
* ======================
*. Analyze T- ampitudes
* ======================
*
      WRITE(6,*) ' KVEC3 before ANA_GCC ', KVEC3
      CALL ANA_GENCC_CC_KRCC(WORK(KCC1),1,
     &                       WORK,KFREE,LFREE)
      WRITE(6,*) ' KVEC3 after ANA_GCC ', KVEC3
      STOP 'AFTER ANALYZE'
*
*
* =========================
* Expectation value energy
* =========================
*
*.  Exp(t) !ref>
      I_DO_CC_EXP = 0
      IF(I_DO_CC_EXP.EQ.1) THEN
      MX_TERM = 100
      ICC_EXC = 1
      XCONV = 1.0D-20
      I_USE_SIMTRH = 0
*
      ICSPC = IETSPC
      ISSPC = IETSPC
*
C     CALL EXPT_REF(LUC,LUSC1,LUHC,LUSC2,LUSC3,XCONV,MX_TERM,
C    &             WORK(KVEC1),WORK(KVEC2),CCTYPE)
*. H Exp T !ref>
      ICC_EXC = 0
C?    WRITE(6,*) ' Input to MV7 '
C?    CALL WRTVCD(WORK(KVEC1),LUSC1,1,-1)
*. Regenerate inactive Fock matrix - different definition of
*. inactive Fock matrix may have been used in CC calculation
      CALL COPVEC(WORK(KINT1O),WORK(KINT1),NINT1)
      I_USE_SIMTRH = 0
      CALL FI(WORK(KINT1),ECCX,1)
      ECORE = ECORE_INI

      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUSC1,LUHC)
C?    WRITE(6,*) ' Output from MV7 '
C?    CALL WRTVCD(WORK(KVEC1),LUHC,1,-1)
*. E = <ref! exp(T)+ H exp(T)!ref>/<ref! exp(T)+exp(T)!ref>
      LBLK = -1
      CHC  = INPRDD(WORK(KVEC1),WORK(KVEC2),LUSC1,LUHC  ,1,LBLK)
      CHHC = INPRDD(WORK(KVEC1),WORK(KVEC2),LUHC,LUHC  ,1,LBLK)
      CC   = INPRDD(WORK(KVEC1),WORK(KVEC2),LUSC1,LUSC1 ,1,LBLK)
      WRITE(6,*)
      WRITE(6,*) ' Energy as expectation value : '
      WRITE(6,*) ' ============================================='
      WRITE(6,*)
      WRITE(6,'(5X,A,E25.12)')
     &'  <ref! exp(T)+ H exp(T)!ref> (- Ecore)= ', CHC
      WRITE(6,'(5X,A,E25.12)')
     &'  <ref! exp(T)+   exp(T)!ref> = ', CC
      WRITE(6,'(5X,A,E25.12)')
     &' Expectation value coupled cluster energy   = ', CHC/CC + ECORE
      WRITE(6,*)
      WRITE(6,*) ' (<CC|H^2|CC> - <CC|H|CC>^2/<CC!CC>)/<CC!CC>  = ',
     &              CHHC/CC - (CHC/CC)**2
*
      WRITE(6,*) 'I_TRANS_WF (2) = ', I_TRANS_WF
      IF(I_TRANS_WF.EQ.1) THEN
        WRITE(6,*) ' cc wf is normalized and transferred to LUC '
        XNORM = SQRT(CC)
        FACTOR = 1.0D0/XNORM
        CALL SCLVCD(LUSC1,LUC,FACTOR,WORK(KVEC1),1,-1)
      END IF
*
      IF(IDENSI.NE.0) THEN
*. Normalize wavefunction
        WRITE(6,*) ' cc wf is normalized and transferred to LUC '
        XNORM = SQRT(CC)
        FACTOR = 1.0D0/XNORM
        CALL SCLVCD(LUSC1,LUSC2,FACTOR,WORK(KVEC1),1,-1)
*. Calculate density matrices in standard expectation value formulation
        LBLK = -1
        CALL COPVCD(LUSC2,LUSC1,WORK(KVEC1),1,LBLK)
        CALL DENSI2(IDENSI,WORK(KRHO1),WORK(KRHO2),
     &       WORK(KVEC1),WORK(KVEC2),LUSC1,LUSC2,EXPS2,
     &       0,WORK(KSRHO1))
      END IF
*.    ^ End if density matrices should be calculated
      IF(NPROP.GT.0) THEN
*. Calculate properties
        KLDUM = 1
        CALL ONE_EL_PROP(1,0,WORK(KLDUM))
      END IF
*.    ^ End if properties should be calculated
      END IF
*     ^ End of expectation values are calculated
C     IF(I_DO_CC3.EQ.1) THEN
C       CALL CC3_JACO(WORK(KCC1),WORK(KVEC1),WORK(KVEC2),
C    &                WORK(KCC2))
C     END IF
*
*. Test Jacobian by explicit calculation
*
      I_TEST_JACO = 0
      IF(I_TEST_JACO.EQ.1) THEN
        WRITE(6,*) ' Control will be transferred to TEST_JACO'
        call quit('subroutine JACO_TEST missing')
!       CALL JACO_TEST(WORK(KCC1),WORK(KVEC1),WORK(KVEC2),
!    &                 WORK(KCC2))
      END IF
*
*
* =======================
* CC excitation energies
* =======================
*
      I_DO_CC_EXC_E = 0
      IF(I_DO_CC_EXC_E.EQ.1) THEN
C
C Linear response works only for real groups. Lasse
C
*. Three more vectors
        CALL MEMMAR(KCC3,LEN_T_VEC_MX,'ADDL  ',2,'CC3_VE')
        CALL MEMMAR(KCC4,LEN_T_VEC_MX,'ADDL  ',2,'CC4_VE')
        CALL MEMMAR(KCC5,LEN_T_VEC_MX,'ADDL  ',2,'CC5_VE')
*
C       CALL CC_EXC_E_REL(WORK(KVEC1),WORK(KVEC2),WORK(KCC1),
C    &                WORK(KCC2),WORK(KCC3),WORK(KCC4),
C    &                WORK(KDIA),II_RES_EXC,ECC)
      END IF
*
      WRITE(6,*)  ' Returning from LUCIA_GENCC_REL '
      CALL MEMMAR(IDUM,IDUM,'FLUSM ',IDUM,'GENCC ')
*
      RETURN
      END


