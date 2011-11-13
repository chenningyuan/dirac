!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl

***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GASCI_REL(F1,F2,INTSPACE,ISM,ISPC,IPRNT,EREF,
     &                     ECORE_KRMC,CIRUN,IPROCLIST,IGROUPLIST)
*
      use memory_allocator
*
* CI optimization in GAS space number ISPC for symmetry ISM
*
* Jeppe Olsen, Winter of 1995
*              Winter of 1998 : Modified for relativistic calculations
*
* Stefan Knecht, Jan 2010: FIXME: introduce use of alloc routines if 
*                                 they can handle int*8 allocations 
*                                 within the ilp32 model.
*                                 For now, we are using the generic 
*                                 allocate/deallocate routines
*
      IMPLICIT REAL*8(A-H,O-Z)
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "ipoist8.inc"
*
      EXTERNAL MV7_REL
*
      character*6 CIRUN
*
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "clunit.inc"
#include "csm.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "oper.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "ctcc.inc"
#include "csfbas.inc"
#include "cintfo.inc"
*
*. Common block for communicating with sigma
#include "cands.inc"
#include "spinfo.inc"
*  nat. orb. occ. num.
#include "noccn_inf.inc"
#include "integrals_off.inc"
#include "dgroup.h"
#include "maxorb.h"
#include "../krmc/dcbbos.h"
*     /CECORE/
      COMMON/CECORER/ECORE
      DIMENSION F1(*), F2(*), IPROCLIST(*), IGROUPLIST(*)
      INTEGER*8 L_MAX_DENS_BLK
      LOGICAL CPIC_PASS_D, DIAGONAL_ALLOC
C 
C     some local scratch
      real*8,    allocatable :: scratch_dav1(:)
      real*8,    allocatable :: scratch_dav2(:)
      integer,   allocatable :: iscratch_dav1(:)
C     memory pointer for TT-block distribution among CPU's
C     memory pointer for TT-block offsets
C     memory pointer for T coefficients
C     memory pointer for T diag block
C     memory pointer for C1 vector block
C     memory pointer for C2 vector block
#if defined (VAR_MPI2)
#include "krmc_shmem.h"
#include "mxdim_mpi2.h"
      INTEGER(KIND=MPI_OFFSET_KIND) ITEST_OFF
#if defined VAR_SHMEM
      INTEGER(KIND=MPI_OFFSET_KIND) IT_TTOL
      POINTER (MY_TTPL_PTR, IT_TTPL(IDUMMY_BUFF_SZ))
      POINTER (MY_TTOL_PTR, IT_TTOL(IDUMMY_BUFF_SZ))
      POINTER (MY_T_PTR   , T_BUFF(IDUMMY_BUFF_SZ))
      POINTER (MY_T_D_PTR , T_BUFF_D(IDUMMY_BUFF_SZ))
      POINTER (MY_XT_PTR  , XT_BUFF(IDUMMY_BUFF_SZ))
      POINTER (MY_VEC1_PTR, VEC1_BUFF(IDUMMY_BUFF_SZ))
      POINTER (MY_VEC2_PTR, VEC2_BUFF(IDUMMY_BUFF_SZ))
      POINTER (MY_DL_PTR  , IT_DLIST_TMP(IDUMMY_BUFF_SZ))
#else
      integer,   allocatable :: IT_TTPL(:)
      integer,   allocatable :: IT_DLIST_TMP(:)
      integer*8, allocatable :: IT_TTOL(:)
      real*8,    allocatable :: T_BUFF(:)
      real*8,    allocatable :: T_BUFF_D(:)
#endif /* VAR_SHMEM */
#endif /* VAR_MPI2 */
C
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'GASCI ')
C
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRNT)
C     set some variables ...
      MIN2  = - 2
      NZERO =   0
      IIONE =   1
C     CPIC_PASS_D controls whether the T_DIAG part is 
C     stored on T_BUFF_D or not
C     DIAGONAL_ALLOC keeps track of allocation of T_BUFF_D
      CPIC_PASS_D    = .FALSE.
      DIAGONAL_ALLOC = .FALSE.
C
CSK   NTEST = 3
      IF(NTEST.GE.3) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ====================================='
        WRITE(LUWRT,*) ' Control has been transferred to GASCI'
        WRITE(LUWRT,*) ' ====================================='
        WRITE(LUWRT,*)
        WRITE(LUWRT,'(A)') '  A few pertinent data : '
        WRITE(LUWRT,*)
        WRITE(LUWRT,'(A,I2)') '  CI space         ',ISPC
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Number of GAS spaces included ',LCMBSPC(ISPC)
        WRITE(LUWRT,'(A,10I3)') ' GAS spaces included           ',
     &               (ICMBSPC(II,ISPC),II=1,LCMBSPC(ISPC))
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Occupation constraints : '
        WRITE(LUWRT,*) '========================= '
        WRITE(LUWRT,*)
        WRITE(LUWRT,*)
        DO JJGASSPC = 1, LCMBSPC(ISPC)
         JGASSPC = ICMBSPC(JJGASSPC,ISPC)
        WRITE(LUWRT,*)
     &  ' Gas space  Min acc. occupation Max acc. occupation '
        WRITE(LUWRT,*)
     &  ' ================================================== '
        DO IGAS = 1, NGAS
          WRITE(LUWRT,'(3X,I2,13X,I3,16X,I3)') IGAS,
     &     IGSOCCX(IGAS,1,JGASSPC),IGSOCCX(IGAS,2,JGASSPC)
        END DO
        END DO
C
       END IF
CSK    NTEST = 0
C
      ECORE = ECORE_KRMC
      NDET = XISPSM(ISM,ISPC)
      write(LUWRT,*) ' Number of determinants/combinations  ',NDET
C     transfer NDET to common block LUCIPARREL
      L_COMBI     = NDET
      L_COMBI_MAX = NDET
C
C     real or complex eigenvectors
      IF(IRELACI.EQ.0) THEN
C       real
        IRC = 1
      ELSE
        IRC = 2
C       well, generally complex, but
        if (NZ.eq.1) then
C         quaternion real
          IRC = 1
        end if
      END IF
      if (IRELACI.eq.1.and.NZ.eq.1) then
        write(LUWRT,*) ' Quaternion real relativistic CI '
      end if
C     control variable for general file structure
      LBLK = -1
C
C     transfer symmetry of C and sigma to CANDS
      ICSM = ISM
      ISSM = ISM
      ICSPC = ISPC
      ISSPC = ISPC
C
C     transfer to LUCIPARREL
      IRC_SAVE = IRC
C
      IF(NOCSF.EQ.1) THEN
        NVAR = NDET
      ELSE
        NVAR = NCSASM(ISM)
      END IF
      IF (IPRNT.GE.5) WRITE(LUWRT,*) '  NVAR in GASCI ', NVAR
C
#if defined (VAR_MPI2)
      IF( IT_SHL .ge. 0 .or. SPLIT_IJKL )THEN
C
C       allocate T-block processor and offset list
        MY_TTPL_LEN = NSPOBEX_TP * IRC
        MY_TTOL_LEN = NMPROC
#if defined VAR_SHMEM
        MY_TTOL_LEN = NSPOBEX_TP * IRC
        CALL MPIXMEM_ALLOC( MY_TTPL_PTR, 2, MY_TTPL_LEN,
     &                      MPI_INFO_NULL, .FALSE.)
        CALL MPIXMEM_ALLOC( MY_TTOL_PTR, 3, MY_TTOL_LEN,
     &                      MPI_INFO_NULL, .FALSE.)
        CALL IMINO8(IT_TTOL,MY_TTOL_LEN)
#else
        allocate(IT_TTPL(MY_TTPL_LEN))
        allocate(IT_TTOL(MY_TTOL_LEN))
        IT_TTPL =  0
        IT_TTOL = -1
#endif
      END IF
#endif
C
C     partition CI vector - for all MS2 spaces into batches
      CALL Z_BLKFO_REL_1(IDC,NMS2VAL,ICSPC,ICSM,KPCIBT,KPCBLTP,
     &                   KNODELIST,KBLOCKLIST,KRCCTOS,
     &                   NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,NBAT_MS2,
     &                   IBAT_MS2,1,0,IDUMMY,IRC,IT_TTPL,
     &                   IGROUPLIST)
C
#if defined (VAR_MPI2)
C
C     define length of T_DIAG block and store on common block 
C     variable MY_T_LEN_D
      LEN_T_D_BUFF = (3*(NTOOB**2)) + (2*NTOOB)
      MY_T_LEN_D   = LEN_T_D_BUFF
C     
C     INFO about T block distribution on IT_TTPL
C     ==> if (shared memory route active) then
C            allocate T buff, (read (ij|kl) in), open T window
C         end if
C     <==
C
      IF( IT_SHL .ge. 0 .or. SPLIT_IJKL )THEN
C
C       calculate length of T_BUFF
        MY_T_LEN     = IGET_MY_T_LEN(IT_TTPL,NSPOBEX_TP,
     &                               N1ELINT,IRC)
        LEN_T_BUFF   = MY_T_LEN
        IF( SPLIT_IJKL ) THEN
C         real part of T array
          LEN_T_BUFF = MY_T_LEN/IRC
        END IF
csk     write(luwrt,*) ' LEN_T_BUFF, MY_T_LEN, MY_T_LEN_D',
csk  &                   LEN_T_BUFF, MY_T_LEN, MY_T_LEN_D
csk     CALL IWRTMAMN(IT_TTPL,1,NSPOBEX_TP,1,NSPOBEX_TP,LUWRT)
C
C       FIXME: provide information if available
C
        IS_DUMMY = 0
C
C       allocate T_BUFF array (replaces use of WORK(KT_CC))
#if defined VAR_SHMEM
        CALL GET_MEM_T( IS_DUMMY )
        CALL GET_MEM_T_D( IS_DUMMY )
C       zero arrays
        CALL DZERO8(T_BUFF,MY_T_LEN)
        CALL DZERO8(T_BUFF_D,MY_T_LEN_D)
#else
        allocate(T_BUFF(MY_T_LEN))
        allocate(T_BUFF_D(MY_T_LEN_D))
        DIAGONAL_ALLOC = .TRUE.
        T_BUFF   =  0
        T_BUFF_D =  0
#endif
C
        IF( SPLIT_IJKL ) THEN
          CALL MPI_GATHER(MY_T_LEN,1,MPI_INTEGER8,IT_TTOL,1,
     &                    MPI_INTEGER8,MASTER,MPI_COMM_WORLD,IERR)
        END IF
C
        WRITE(LUWRT,'(/2X,A)')
     & '=============================================================='
        WRITE(LUWRT,'(2X,A/)')
     & '==>  allocation of process unique T (integral) elements    <=='
        DO II = 1, NMPROC
        WRITE(LUWRT,'(2X,A,1X,I4,1X,A,I12,1X,A)')
     & 'process',II-1,'allocates for ',IT_TTOL(II),'elements          '
        END DO
      
        WRITE(LUWRT,'(2X,A/)')
     & '=============================================================='

        IF(SPLIT_IJKL) deallocate(IT_TTOL)
C
C       read in all (active) blocked integrals (T blocks)
C       ==> integrals required for diagonal H on T_BUFF_D
C
        IF (CIRUN .ne. 'ANALYZ' ) THEN
C           FIXME: 11-02-08 SK: definition of IBU2D?
          call picasso(INTSPACE,ISO1,IPRNT, .TRUE. ,
     &                 NL2D,IBU2D,ISM2D,F1,F2,T_BUFF,
     &                 T_BUFF_D,IT_TTPL,IBOSYM,CIRUN,.FALSE.)
          CPIC_PASS_D = .TRUE.
        END IF
C
csk     CALL WRTMATMN8(T_BUFF,1,MY_T_LEN,1,MY_T_LEN,LUWRT)
csk     CALL WRTMATMN8(T_BUFF_D,1,MY_T_LEN_D,1,MY_T_LEN_D,LUWRT)
C
#if defined VAR_SHMEM
        IF( IT_SHL .ge. 0 )THEN
C
C         open memory window (collective call for communicator
C         group MY_NEWCOMM_SM)
C
          CALL GET_T_WIN( T_BUFF, MYNEW_COMM_SM, MY_T_LEN )
C
          CALL CALC_WIN_OFFSET(IT_TTPL,IT_TTOL,NSPOBEX_TP,IRC)
C
csk        WRITE(LUWRT,*) ' printing offset array IT_TTOL'
csk        CALL IWRTMAMN8(IT_TTOL,1,MY_TTOL_LEN,1,MY_TTOL_LEN,LUWRT)
        END IF
#endif
C
C       transfer IS_LENGHT_TT to common block LUCIPAR_I8
C       == max. T block length
C
        IS_LENGTH_TT = 0
        IS_LENGTH_TT = NCALC_MAX_BLK(IRC)
C
      END IF
C     ... read-in reordered integrals for the old-school way.
      IF( .NOT. CPIC_PASS_D .and. REORD_IJKL ) THEN
C
         IF( .NOT. DIAGONAL_ALLOC ) allocate(T_BUFF_D(MY_T_LEN_D))
         DIAGONAL_ALLOC = .TRUE. 
C
CSK      CALL DZERO8(T_BUFF,MY_T_LEN)
         IF (CIRUN .ne. 'ANALYZ' ) THEN
           call picasso(INTSPACE,ISO1,IPRNT, .TRUE. ,
     &                  NL2D,IBU2D,ISM2D,F1,F2,T_BUFF,
     &                  T_BUFF_D,IT_TTPL,IBOSYM,CIRUN,.FALSE.)
           CPIC_PASS_D = .TRUE.
         END IF
      END IF
#endif
C
      CALL Z_BLKFO_REL_2(IDC,NMS2VAL,ICSPC,ICSM,KPCIBT,KPCBLTP,
     &                   NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,NBAT_MS2,
     &                   IBAT_MS2,1,0,IDUMMY,IRC)
C
#if defined (VAR_MPI2)
C
      IF ( CIRUN .eq. 'DIAG  ' .or. CIRUN .eq. 'CIINII' .or.
     &     CIRUN .eq. 'KR-CI ' .or.
     &     CIRUN .eq. 'SIGMA ' .or. CIRUN .eq. 'SIGMAD' .or. 
     &     CIRUN .eq. 'DENS1 ' .or. CIRUN .eq. 'DENS2 ' ) THEN
C
C       dummy file array allocation for 'remaining' CPUs
        IF( IAM_NOT_INV .eq. 0 ) THEN
          IALL_LU1 = 1
          IALL_LU2 = 1
          IALL_LU3 = 1
          IALL_LU4 = 1
          IALL_LU5 = 1
          IALL_LUC = NUM_BLOCKS2
        END IF
C
        CALL CALC_OFF_MPI_FILE('ILUC  ',IGROUPLIST,WORK(KBLOCKLIST),
     &                         WORK(KNODELIST),1)
C
        CALL CALC_OFF_MPI_FILE('ILU1  ',IGROUPLIST,WORK(KBLOCKLIST),
     &                         WORK(KNODELIST),MXCIV_CI+NROOT)
C
        CALL CALC_OFF_MPI_FILE('ILU2  ',IGROUPLIST,WORK(KBLOCKLIST),
     &                         WORK(KNODELIST),MXCIV_CI+NROOT)
C
C       initialize general offset for ILUC
        MY_LUC_OFF = 0
C
C       length for allocation of file arrays
        IALL_LU1 = ( MXCIV_CI + NROOT ) * MY_ACT_BLK2
        IALL_LU2 = ( MXCIV_CI + NROOT ) * MY_ACT_BLK2
        IALL_LUC = NUM_BLOCKS2
C
C       allocate file arrays
        CALL MEMMAR(KILU1LIST,IALL_LU1,'ADDS  ',1,'LU1LST')
        CALL MEMMAR(KILU2LIST,IALL_LU2,'ADDS  ',1,'LU2LST')
        CALL MEMMAR(KILUCLIST,IALL_LUC,'ADDS  ',1,'LUCLST')
C       initialize ...
        CALL IZERO(WORK(KILU1LIST),IALL_LU1)
        CALL IZERO(WORK(KILU2LIST),IALL_LU2)
        CALL IZERO(WORK(KILUCLIST),IALL_LUC)
C
        IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'SIGMAD' .or. 
     &      CIRUN .eq. 'DIAG  ' .or. CIRUN .eq. 'KR-CI ' ) THEN 
C
          CALL CALC_OFF_MPI_FILE('LUDIA ',IGROUPLIST,WORK(KBLOCKLIST),
     &                           WORK(KNODELIST),1)
C
          IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'KR-CI ' ) THEN
C
            CALL CALC_OFF_MPI_FILE('ILU3  ',IGROUPLIST,WORK(KBLOCKLIST),
     &                             WORK(KNODELIST),NROOT)
C
            CALL CALC_OFF_MPI_FILE('ILU4  ',IGROUPLIST,WORK(KBLOCKLIST),
     &                             WORK(KNODELIST),1)
C
            CALL CALC_OFF_MPI_FILE('ILU5  ',IGROUPLIST,WORK(KBLOCKLIST),
     &                             WORK(KNODELIST),1)
C
C           some more file arays ... 
            IALL_LU3 = ( NROOT ) * MY_ACT_BLK2
            IALL_LU4 = ( 1 )     * MY_ACT_BLK2
            IALL_LU5 = ( 1 )     * MY_ACT_BLK2
C           allocate ...
            CALL MEMMAR(KILU3LIST,IALL_LU3,'ADDS  ',1,'LU3LST')
            CALL MEMMAR(KILU4LIST,IALL_LU4,'ADDS  ',1,'LU4LST')
            CALL MEMMAR(KILU5LIST,IALL_LU5,'ADDS  ',1,'LU5LST')
C           initialize ...
            CALL IZERO(WORK(KILU3LIST),IALL_LU3)
            CALL IZERO(WORK(KILU4LIST),IALL_LU4)
            CALL IZERO(WORK(KILU5LIST),IALL_LU5)
C
          END IF
C         ^ CIRUN == CIINII or KR-CI
        END IF            
C       ^ CIRUN == CIINII .or. SIGMAD .or. DIAG .or. KR-CI
      END IF
C     ^ CIRUN == DIAG .or. CIRUN == CIINII .or. CIRUN == SIGMA ...
#endif
C
C     allocate memory : 2 vector blocks 
C     store veclen for  1 vector/scratch block
      CALL GET_3BLKS_REL(KVEC1,KVEC2,LSCR2,CIRUN)
C
C     CI diagonal - if required
*+++++++++++++++++++++++++++++++++++++++++++++++++++
      if (CIRUN.eq.'DIAG  '.or.CIRUN.eq.'CIINII'.or.
     &    CIRUN.eq.'KR-CI '.or.
     &    CIRUN.eq.'SIGMA '.or.CIRUN.eq.'SIGMAD') then
*
        if (CIRUN.ne.'SIGMA ') then
*+++++++++++++++++++++++++++++++++++++++++++++++++++
*
      IF(IDIAG.EQ.2) THEN
        LUDIA = LUSC1
      END IF
#if defined (VAR_MPI2)
      IF( USE_EX_IDIA ) THEN
        WRITE(LUWRT,*) '  saving CPU/WALL-time in this MPI run   '
        WRITE(LUWRT,*) '  using existing diagonal file DIAPAR.1. '
        WRITE(LUWRT,*) '                                         '
        GOTO 55
      END IF
#endif
      IF(.NOT.(IDIAG.EQ.2.AND.IRESTR.EQ.1)) THEN
        ECOREP = ECORE
C
        IF(ICISTR.GE.2) CALL REWINE(LUDIA,-1)
        IRECO = 1
        CALL GASDIAT_REL(WORK(KVEC1),LUDIA,ECOREP,ICISTR,IHAM12,
     &               WORK(KPCBLTP),NBLOCK,WORK(KPCIBT),
     &               KPICI1,INTSPACE,T_BUFF,T_BUFF_D,IDENSI,
     &               NIRREP,ISPINFREE,
     &               ENVIRO,SIGDEN_ROUTE,IRECO,WORK(KNODELIST),
     &               CPIC_PASS_D)
C       end of real part
        CALL ITODS(-1,1,-1,LUDIA)
        IF (IRC.EQ.2) THEN
C         diagonal once again- for the imaginary part
          IRECO = 2
          CALL GASDIAT_REL(WORK(KVEC1),LUDIA,ECOREP,ICISTR,IHAM12,
     &                 WORK(KPCBLTP),NBLOCK,WORK(KPCIBT),
     &                 KPICI1,INTSPACE,T_BUFF,T_BUFF_D,IDENSI,
     &                 NIRREP,ISPINFREE,
     &                 ENVIRO,SIGDEN_ROUTE,IRECO,WORK(KNODELIST),
     &                 CPIC_PASS_D)
C       end of imaginary  part
        CALL ITODS(-1,1,-1,LUDIA)
        END IF
        IF(IPRCIX.GE.2) WRITE(LUWRT,*) ' Diagonal constructed  '
        if (IMEMCK.gt.0) then
          write(LUWRT,*) ' Memory map after GASDIAT_REL in GASCI_REL'
          IDUM = 0
          call memmar(KDUM,IDUM,'PRINT ',IDUM,'IDUMMY')
        end if
      ELSE
         WRITE(LUWRT,*) ' Diagonal not calculated '
      END IF
c
#if defined (VAR_MPI2)
C
 55   CONTINUE
C
      NPTEST_VAR = 00
      IF( NPTEST_VAR .ge. 10 ) THEN
      IRECO = IRC_SAVE
      ITEST_OFF = 0
      DO ILPRECO = 1, IRECO
        ITEST_OFF = MY_DIA_OFF + MY_VEC1_IOFF * ( ILPRECO - 1 )
        DO IIBLK = 1, NUM_BLOCKS
          ILEN = 0
          CALL GET_BLOCK_PROC(WORK(KNODELIST),IIBLK,IPROC)
          IF( MYPROC .eq. IPROC ) THEN
            CALL GET_BLOCK_LENGTH(WORK(KBLOCKLIST),IIBLK,ILEN)
            CALL MPI_FILE_READ_AT(IDIA,ITEST_OFF,WORK(KVEC1),ILEN,
     &                            MPI_REAL8,ISTAT,IERR)
            WRITE(LUWRT,*) ' Diagonal elements',ILEN
            CALL WRTMATMN(WORK(KVEC1),1,ILEN,1,ILEN,LUWRT)
            ITEST_OFF = ITEST_OFF + ILEN
          END IF
        END DO
      END DO
      
      ENDIF
      NPTEST_VAR = 00
#endif
c
*+++++++++++++++++++++++++++++++++++++++++++++++++++
        end if
        if (CIRUN.ne.'DIAG  ') then
*+++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     transfer control to optimization routine
c
      MINST = 1
      IF(IRESTR.EQ.0) THEN
        INICI = 0
      ELSE
        INICI = -1
      END IF
      NPRDET = 0
C     this is a CI and not a perturbation calculation
      IPERTOP = 0
C
C     restart from previous vector(s)
C
#if defined (VAR_MPI2)
C
      IF( INICI .lt. 0 )THEN
        if ( CIRUN.eq.'SIGMAD'.or.CIRUN.eq.'SIGMA ') GOTO 88 
          WRITE(LUWRT,'(/A,A)') '  MPI restart file will be'//
     &    ' prepared from C-vector file: KRCI_CVECS.',SYMFLABEL
C
         IF( MYPROC .eq. MASTER ) CALL REWINE(LUC,-1)
C
         CALL COP_REST_VEC_REL(WORK(KVEC1),LUC,ILU1,MY_LU1_OFF,
     &                         WORK(KILU1LIST),
     &                         WORK(KNODELIST),WORK(KBLOCKLIST),
     &                         NBLOCK,NROOT,IRC)
C
        WRITE(LUWRT,'(A/)') '  file conversion from standard I/O'//
     &                      ' to MPI-IO finished.'
C
      END IF
 88   CONTINUE
C
#endif
C
C
*+++++++++++++++++++++++++++++++++++++++++++++++++++
          if (CIRUN.eq.'SIGMAD'.or.CIRUN.eq.'SIGMA ') then
*+++++++++++++++++++++++++++++++++++++++++++++++++++
C
C         C-vector resides on MASTER in file LUC
C
C                        ============
C                        PARALLEL RUN
C                        ============
C
C         NODES recieve the appropriate C-blocks
C
C         MASTER recieves the corresponding SIGMA-blocks to 
C         write on LUHC
C
C         IF( SIGMAD )   MASTER recieves the corresponding diagonal
C                        elements to write on LUDIA
C
#if defined (VAR_MPI2)
C
            IF( MYPROC .eq. MASTER ) CALL REWINE(LUC,-1)
C
            CALL COP_REST_VEC_REL(WORK(KVEC1),LUC,ILU1,MY_LU1_OFF,
     &                            WORK(KILU1LIST),
     &                            WORK(KNODELIST),WORK(KBLOCKLIST),
     &                            NBLOCK,IIONE,IRC)
C
C            Arrays for partitioning of the CI vector
C            (*_MS2 arrays are in CSTATE )
             call z_blkfo_rel_par(IDC,NMS2VAL,ISPC,ICSM,
     &                            KLVLBT,KLVLEBT,KLVI1BT,KLVIBT,KDUM,
     &                            NBATCH_LD,NM_BLK_LD,NBLK_MS2_LD,
     &                            IBLK_MS2_LD,NBAT_MS2_LD,IBAT_MS2_LD,
     &                            NZERO,IIONE,WORK(KNODELIST))
C
C           set offset for sigma-file
            JVEC_SF = 0
C
C           erase old LUCLIST
            CALL IZERO(WORK(KILUCLIST),IALL_LUC)
C
            IF( IRC .eq. 1 ) THEN
C             real quaternion double group
              CALL COPVCD_PP_CC_B(ILU1,ILUC,WORK(KVEC1),NBATCH_LD,
     &                            WORK(KLVLBT),WORK(KLVLEBT),
     &                            WORK(KLVI1BT),WORK(KLVIBT),
     &                            MY_LU1_OFF,MY_LUC_OFF,
     &                            WORK(KILU1LIST),WORK(KILUCLIST),
     &                            WORK(KBLOCKLIST),NZERO)
            ELSE
C             complex quaternion double group
              CALL COPVCD_PP_CC_B_C(ILU1,ILUC,WORK(KVEC1),NBATCH_LD,
     &                              WORK(KLVLBT),WORK(KLVLEBT),
     &                              WORK(KLVI1BT),WORK(KLVIBT),
     &                              MY_LU1_OFF,MY_LUC_OFF,
     &                              WORK(KILU1LIST),WORK(KILUCLIST),
     &                              WORK(KBLOCKLIST),NZERO)
            END IF
C
#endif
C           start calculation: sigma = H x C
C
************************************************************************
            call sigden_ctrl(WORK(KVEC1),WORK(KVEC2),
#if defined (VAR_MPI2)
     &                       ILUC,ILU2,
#else
     &                       LUC,LUHC,
#endif
     &                       WORK(KT_CC),T_BUFF,1
#if defined (VAR_MPI2)
     &                      ,WORK(KILUCLIST),WORK(KILU2LIST),
     &                       WORK(KBLOCKLIST),WORK(KNODELIST),
     &                       IGROUPLIST,IPROCLIST,WORK(KRCCTOS),
     &                       IT_TTPL,IT_TTOL
     &                      ,WORK(KBLOCKLIST),WORK(KNODELIST)
#endif
     &                       )
************************************************************************
C
#if defined (VAR_MPI2)
C
            IF( IAM_NOT_INV .eq. 0 ) GOTO 405
C
            IF( MYPROC .eq. MASTER )  CALL REWINE(LUHC,-1)
C
            LBLK = -1
C
            CALL COPVCD_PAR_BDRIV_REL(ILU2,LUHC,WORK(KVEC1),
     &                                WORK(KNODELIST),NUM_BLOCKS,
     &                                WORK(KBLOCKLIST),
     &                                MPI_COMM_WORLD,LBLK,IIONE,
     &                                WORK(KILU2LIST),MY_LU2_OFF,IRC)
C
C
            IF( CIRUN .eq. 'SIGMAD' ) THEN
C
C             save H diagonal on MASTER
C
              IF( MYPROC .eq. MASTER ) CALL REWINE(LUDIA,LBLK)
C             H is a hermitian matrix...
              CALL COPVCD_PAR_BDRIV_DIA_REL2(IDIA,LUDIA,WORK(KVEC1),
     &                                       WORK(KNODELIST),NUM_BLOCKS,
     &                                       WORK(KBLOCKLIST),
     &                                       MPI_COMM_WORLD,LBLK,IIONE,
     &                                       MY_DIA_OFF)
            END IF
C           ^ CIRUN .eq. SIGMAD
 405    CONTINUE
C
#endif
          else
 
            if (NROOT.gt.MXROOT) then
              WRITE(LUWRT,'(A,I4)') ' Too many roots. MXROOT is ',MXROOT
              WRITE(LUWRT,'(A,I4)') ' Increase MXROOT to ',NROOT
              CALL QUIT(' *** ERROR in GASCI_REL: too many roots ***')
            end if
C
C           allocate some scratch space
            NFINDM = min(NVAR,4*NROOT)
            call alloc(scratch_dav1, 4*NFINDM+1)
            call alloc(iscratch_dav1,5*NFINDM)
            call alloc(scratch_dav2, NROOT)
C
            CALL CIEIG5_REL(MV7_REL,INICI,EROOT,WORK(KVEC1),WORK(KVEC2),
     &                      MINST,LUDIA,LUC,LUHC,LUSC1,LUSC2,LUSC3,
     &                      LUSC34,LUSC35,NVAR,NBLK,NROOT,MXCIV_CI,
     &                      MAXIT,MXROOT,LUCIVI,IPRNT,WORK(KSBEVC),
     &                      NPRDET,WORK(KSBIDT),MXP1,MXP2,MXQ,
     &                      scratch_dav1,iscratch_dav1,scratch_dav2,
     &                      WORK(KT_CC),T_BUFF,ECORE,ICISTR,LBLK,
     &                      IDIAG,THRES_E,NBATCH,WORK(KPCLBT),
     &                      WORK(KPCLEBT),WORK(KPCIBT),INIDEG,
     &                      E_THRE,C_THRE,E_CONV,C_CONV,ICLSSEL,NOCCLS,
     &                      IRC,IDCOMH,IHAM12,SIGDEN_ROUTE,
     &                      WORK(KBLOCKLIST),WORK(KNODELIST),
     &                      IPROCLIST,IGROUPLIST
#if defined (VAR_MPI2)
     &                      ,WORK(KILU1LIST),WORK(KILU2LIST),
     &                      WORK(KILU3LIST),WORK(KILU4LIST),
     &                      WORK(KILU5LIST),WORK(KILUCLIST),
     &                      WORK(KRCCTOS),IT_TTPL,IT_TTOL
#endif
     &                      )
C
C           save energy of the lowest eigenstate
            EREF = EROOT(IRFROOT)
C
C           release local scratch memory
            call dealloc(scratch_dav2)
            call dealloc(iscratch_dav1)
            call dealloc(scratch_dav1)
C
#if defined (VAR_MPI2)
C
C           delete big scratch files when doing large-scale CI run
C           (CIRUN .eq. KR-CI) before collecting all vectors
            CALL SETUNITS_PAR_CLOSE_REL(CIRUN,1)
C
            IF( IAM_NOT_INV .eq. 0 ) GOTO 505
C
            LBLK = -1
C           copy c-vectors from nodes and master back to the master
C           purpose: master does the ana-stuff
C
C           NOTE: only C-vectors will be saved.
C
            IF( MYPROC .eq. MASTER ) CALL REWINE(LUC,-1)
            CALL REWINE(LUSC41,-1)
            DO JROOT = 1, NROOT
              CALL COPVCD_PAR_BDRIV_REL(ILU1,LUSC41,WORK(KVEC1),
     &                                  WORK(KNODELIST),NUM_BLOCKS,
     &                                  WORK(KBLOCKLIST),
     &                                  MPI_COMM_WORLD,LBLK,JROOT,
     &                                  WORK(KILU1LIST),MY_LU1_OFF,IRC)
            END DO
            CALL SETUNITS_PAR_CLOSE_REL(CIRUN,2)
            IF(MYPROC.EQ.MASTER) THEN
              CALL REWINE(LUSC41,-1)
              CALL REWINE(LUC,-1)
              DO JROOT = 1, NROOT
                CALL copvcdc(LUSC41,LUC,WORK(KVEC1),0,IRC,LBLK)
              ENDDO
              CALL REWINE(LUC,-1)
            END IF
            NPTESTVAR = 00
            IF( NPTESTVAR .ge. 10 )THEN
              IF(MYPROC.EQ.MASTER) THEN
                DO IVEC = 1, NROOT
                  WRITE(LUWRT,*) '  final new C vector',IVEC
                  CALL WRTVCDC(WORK(KVEC1),LUC,0,IRC,LBLK)
                END DO
              END IF
            END IF
            NPTESTVAR = 00
C
 505        CONTINUE
            IF( IAM_NOT_INV .eq. 0 ) 
     &      CALL SETUNITS_PAR_CLOSE_REL(CIRUN,2)
#endif
C
*+++++++++++++++++++++++++++++++++++++++++++++++++++
          end if
C         ^ CIRUN ne SIGMAD or SIGMA
        end if
C       ^ CIRUN ne DIAG
      else if (CIRUN.eq.'ANALYZ') then
C     analysis of optimized wave function
*+++++++++++++++++++++++++++++++++++++++++++++++++++
C
#if defined (VAR_MPI2)
C
C     only MASTER analyzes the optimized wave function ...
C
      IF( MYPROC .ne. MASTER ) GOTO 888
C
#endif
      LBLK = -1
      call REWINE(LUC,-1)
      NPTESTVAR = 00
      IF( NPTESTVAR .ge. 10 )THEN
        IF(MYPROC.EQ.MASTER) THEN
          DO IVEC = 1, NROOT
            WRITE(LUWRT,*) '  final new C vector before ana',IVEC
            CALL WRTVCDC(WORK(KVEC1),LUC,0,IRC,LBLK)
          END DO
        END IF
      END IF
      NPTESTVAR = 00
      call REWINE(LUC,-1)
      do JROOT=1,NROOT,1
        call REWINE(LUSC61,-1)
        call copvcdc(LUC,LUSC61,WORK(KVEC2),0,IRC,LBLK)
        call REWINE(LUSC61,-1)
        WRITE(LUWRT,'(/A,I3)')
     &  ' ***************************************************'
        WRITE(LUWRT,'(A,I3)')
     &  ' Analysis of Optimized Wave Function for ROOT = ',JROOT
        WRITE(LUWRT,'(A,I3)')
     &  ' ***************************************************'
C
        call gasana_rel(WORK(KVEC1),NBLOCK,WORK(KPCIBT),
     &                  WORK(KPCBLTP),LUSC61,LUSC37,
     &                  ICISTR,IRC)
      end do
      call REWINE(LUC,-1) 
      NPTESTVAR = 00
      IF( NPTESTVAR .ge. 10 )THEN
        IF(MYPROC.EQ.MASTER) THEN
          DO IVEC = 1, NROOT
            WRITE(LUWRT,*) '  final new C vector after ana',IVEC
            CALL WRTVCDC(WORK(KVEC1),LUC,0,IRC,LBLK)
          END DO
        END IF
      END IF
      NPTESTVAR = 00

 888  CONTINUE 
#if defined (VAR_MPI2)
C     let the slaves wait... it is only a short way to go so why not?
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#endif
C
*+++++++++++++++++++++++++++++++++++++++++++++++++++
      else if (CIRUN.eq.'DENS1 '.or.CIRUN.eq.'DENS2 ') then
C
C     analysis of density and occupation
C
*+++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     real or complex eigenvectors
      IF(IRELACI.EQ.0) THEN
C       real
        IRC = 1
      ELSE
C       well, complex
        IRC = 2
        if (NZ.eq.1) then
C         quaternion real
          IRC = 1
        end if
      END IF
C
      LBLK = -1
C
      if (IDENSI.ne.0) then
        if (IDENSI.eq.1.or.IDENSI.eq.2) then
          write(LUWRT,*) 'IDENSI is ',IDENSI
          write(LUWRT,*) 
     &    'I am not allowed to calculate non-relativistic'
          write(LUWRT,*) ' density matrices YET.'
          write(LUWRT,*) ' Please ask the experts for instructions'
          Call Abend2( 'LUCIAREL quitting' )
        end if
#if defined (VAR_MPI2)
C
C       arrays for partitioning of the CI vector
C       (*_MS2 arrays are in CSTATE )
        call z_blkfo_rel_par(IDC,NMS2VAL,ISPC,ICSM,
     &                       KLVLBT,KLVLEBT,KLVI1BT,KLVIBT,KDUM,
     &                       NBATCH_LD,NM_BLK_LD,NBLK_MS2_LD,
     &                       IBLK_MS2_LD,NBAT_MS2_LD,IBAT_MS2_LD,
     &                       0,1,WORK(KNODELIST))
C
#endif
C       
C       calculation of densities for natural orbital occ. numbers
C       ---------------------------------------------------------
        IF( NOOCCN_DENS_RUN ) THEN
C
          IF( MYPROC .eq. MASTER ) THEN
            CALL REWINE(LUSC61,LBLK)
            CALL REWINE(LUSC62,LBLK)
            CALL SKPVCD(LUC,IRC * (IDENSI_STATE - 1),WORK(KVEC1),1,LBLK)
          END IF 
#if defined (VAR_MPI2)
          IVEC_TMP = 1
C         right hand-side vector
          CALL COP_REST_VEC_REL(WORK(KVEC1),LUC,ILU1,MY_LU1_OFF,
     &                          WORK(KILU1LIST),
     &                          WORK(KNODELIST),WORK(KBLOCKLIST),
     &                          NUM_BLOCKS,IVEC_TMP,IRC)
C          
          IF( MYPROC .eq. MASTER) THEN
            CALL SKPVCD(LUC,IRC * (IDENSI_STATE - 1),WORK(KVEC1),1,LBLK)
          END IF
C         left hand-side vector
          CALL COP_REST_VEC_REL(WORK(KVEC1),LUC,ILU2,MY_LU2_OFF,
     &                          WORK(KILU2LIST),
     &                          WORK(KNODELIST),WORK(KBLOCKLIST),
     &                          NUM_BLOCKS,IVEC_TMP,IRC)
C       files are in place now... say hello
        WRITE(LUWRT,*) ' '
        WRITE(LUWRT,'(12X,A)') ' ==================================== '
        WRITE(LUWRT,'(12X,A)') ' PARALLEL 1p DENSITY MAT. CALCULATION '
        WRITE(LUWRT,'(12X,A)') ' ==================================== '
        WRITE(LUWRT,*) ' '
C
C
C
#else
C         right hand-side vector to LUSC61
          CALL COPVCDC(LUC,LUSC61,WORK(KVEC1),0,IRC,LBLK)
C
          CALL SKPVCD(LUC,IRC * (IDENSI_STATE - 1),WORK(KVEC1),1,LBLK)
C         left hand-side vector on LUSC62
          CALL COPVCDC(LUC,LUSC62,WORK(KVEC1),0,IRC,LBLK)
#endif
C
          GOTO 895
C
        END IF
C
C       prepare vector files for pure dens. matrix calc. run (parallel)
C       ---------------------------------------------------------------
#if defined (VAR_MPI2)
        IF( MYPROC .eq. MASTER ) THEN
          CALL REWINE(LUSC61,LBLK)
          CALL REWINE(LUSC62,LBLK)
        END IF
C       right hand-side vector on LUSC61
        CALL COP_REST_VEC_REL(WORK(KVEC1),LUSC61,ILU1,MY_LU1_OFF,
     &                        WORK(KILU1LIST),
     &                        WORK(KNODELIST),WORK(KBLOCKLIST),
     &                        NUM_BLOCKS,NROOT,IRC)
C       left hand-side vector on LUSC62
        CALL COP_REST_VEC_REL(WORK(KVEC1),LUSC62,ILU2,MY_LU2_OFF,
     &                        WORK(KILU2LIST),
     &                        WORK(KNODELIST),WORK(KBLOCKLIST),
     &                        NUM_BLOCKS,NROOT,IRC)
C
C       files are in place now... say hello
        WRITE(LUWRT,*) ' '
        WRITE(LUWRT,'(12X,A)') ' =================================== '
        WRITE(LUWRT,'(12X,A)') ' PARALLEL DENSITY MATRIX CALCULATION '
        WRITE(LUWRT,'(12X,A)') ' =================================== '
        WRITE(LUWRT,*) ' '
C
#endif
C
 895    CONTINUE
C
C       one-body density matrix
C       spin-dependent first-order density matrix
C
        JROOT = 1
        JROOT_SCR = 0
C
        IF( NOOCCN_DENS_RUN ) THEN
           JROOT_SCR = IDENSI_STATE 
        ELSE
           JROOT_SCR = 1
        END IF
C
        if (IRELACI.eq.1) then
          INONREL = 0
        else
          INONREL = 1
        end if
        if (IDENSI.eq.3) then
          IDEN12 = 1
          write(LUWRT,*) ' '
          write(LUWRT,'(2X,A)') 
     &    ' Computing spin-dependent one-particle density'
          write(LUWRT,*) ' '
        else if (IDENSI.eq.4) then
          IDEN12 = 2
          write(LUWRT,*) ' '
          write(LUWRT,'(2X,A)') 
     &   ' Computing spin-dependent one- and two-particle densities.'
          write(LUWRT,*) ' '
        else
          write(LUWRT,*) 'IDENSI =  ',IDENSI
          Call Abend2( 'Problemo in gasci. Quitting.' )
        end if
C
C
        if (IPRDEN .ge. 500 .and. .not. NOOCCN_DENS_RUN) then
          IF( MYPROC .eq. MASTER )THEN
            write(LUWRT,*) 'Current vector',JROOT_SCR,' on unit',LUSC61
            call wrtvcdc(WORK(KVEC1),LUSC61,1,IRC,LBLK)
            write(LUWRT,*) 'Current vector ',JROOT_SCR,' on unit',LUSC62
            call wrtvcdc(WORK(KVEC2),LUSC62,1,IRC,LBLK)
          END IF
        end if
C
        if (SIGDEN_ROUTE.eq.'TRADIT') then
          call dbg_rho(WORK(KVEC1),WORK(KVEC2),INONREL,IDEN12,JROOT)
        else if (SIGDEN_ROUTE.eq.'SIGDEN') then
C
          if (CIRUN.eq.'DENS1 ') then
            NTERMS = N1ELINT
          else if (CIRUN.eq.'DENS2 ') then
            NTERMS = NALLINT
          end if
C
C         erase T amplitudes from sigma calculation
C
          IF( IT_SHL .ge. 0 .or. SPLIT_IJKL )THEN
            call DZERO8(T_BUFF,MY_T_LEN)
          ELSE
            call DZERO8(WORK(KT_CC),NALLINT)
          END IF
C
csk    WRITE(LUWRT,*) ' start: MY T (BUFF) array'
csk    IF( SPLIT_IJKL )THEN
csk      CALL WRTMATMN8(T_BUFF,1,MY_T_LEN,1,MY_T_LEN,LUWRT)
csk    ELSE
csk      CALL WRTMATMN8(WORK(KT_CC),1,NALLINT,1,NALLINT,LUWRT)
csk    END IF
C
#if defined (VAR_MPI2)
C         erase old LUCLIST
          CALL IZERO(WORK(KILUCLIST),IALL_LUC)
C         copy right-hand vector to scratch file
          IF( IRC .eq. 1 ) THEN
C           real quaternion double group
            CALL COPVCD_PP_CC_B(ILU1,ILUC,WORK(KVEC1),NBATCH_LD,
     &                          WORK(KLVLBT),
     &                          WORK(KLVLEBT),WORK(KLVI1BT),
     &                          WORK(KLVIBT),MY_LU1_OFF,MY_LUC_OFF,
     &                          WORK(KILU1LIST),WORK(KILUCLIST),
     &                          WORK(KBLOCKLIST), JROOT - 1 )
          ELSE
C           complex quaternion double group
            CALL COPVCD_PP_CC_B_C(ILU1,ILUC,WORK(KVEC1),NBATCH_LD,
     &                            WORK(KLVLBT),WORK(KLVLEBT),
     &                            WORK(KLVI1BT),WORK(KLVIBT),
     &                            MY_LU1_OFF,MY_LUC_OFF,
     &                            WORK(KILU1LIST),WORK(KILUCLIST),
     &                            WORK(KBLOCKLIST), JROOT - 1 )
          END IF
C
#endif
C
C         set offset for left-hand vector on ILU2 
C         (parallel calculation)
C
          JVEC_SF = JROOT - 1
C
C         compute density matrix elements
C         -------------------------------
C
*=======================================================================
          call sigden_ctrl(WORK(KVEC1),WORK(KVEC2),
#if defined (VAR_MPI2)
     &                     ILUC,ILU2,
#else
     &                     LUSC61,LUSC62,
#endif
     &                     WORK(KT_CC),T_BUFF,2
#if defined (VAR_MPI2)
     &                    ,WORK(KILUCLIST),WORK(KILU2LIST),
     &                     WORK(KBLOCKLIST),WORK(KNODELIST),
     &                     IGROUPLIST,IPROCLIST,WORK(KRCCTOS),
     &                     IT_TTPL,IT_TTOL
     &                    ,WORK(KBLOCKLIST),WORK(KNODELIST)
#endif
     &                   )
*=======================================================================
C
C         file INDICES.x
C
C         IF( .NOT. SHARED_M ) THEN
C           read indices from cpu-specific file INDICES.'MYPROC'
C           ( written in dist_ints_sigden --> interface_r.F )
C         ELSE
C           read indices from node-master specific file 
C           ( written in dist_ints_sigden --> interface_r.F )
C         ENDIF
C
C         SK-30-11-2007: since only LUCIAREL is parallel at present
C                        it is sufficient to allocate the 4-INDEX-array
C                        on the real (!) master
C                         
C         FIXME!!!       HAS TO BE CHANGED AS SOON AS THE COMPLETE
C                        KR-MCSCF IS RUNNING IN PARALLEL!
C
          IF( MYPROC .eq. MASTER )THEN
            IF( NOOCCN_DENS_RUN ) THEN
              WRITE(LUWRT,*) ' '
              WRITE(LUWRT,*) ' Allocating 4-index array of length',
     &                         NTERMS*4
              WRITE(LUWRT,*) ' '
C
              call memmar_i8(KINDEX,NTERMS*4,'ADDL  ',1,'KINDEX')
              call IZERO8(WORK(KINDEX),NTERMS*4)
C
            ELSE
              WRITE(LUWRT,*) ' Allocating 4-index array of length',
     &                         NALLINT*4
              WRITE(LUWRT,*) ' '
C
              call memmar_i8(KINDEX,NALLINT*4,'ADDL  ',1,'KINDEX')
              call IZERO8(WORK(KINDEX),NALLINT*4)
C
            END IF
C
            call read_indices(WORK(KINDEX),NTERMS)
C
          END IF
C         ^ MASTER !!!
          IPRDEN = 1
C
#if defined (VAR_MPI2)
C
csk       WRITE(LUWRT,*) ' MY index array '
csk       DO IJINDEX = 1, NTERMS 
csk       WRITE(LUWRT,'(A,I4)') 'index: ',WORK(KINDEX+IJINDEX-1)
csk       END DO
csk    WRITE(LUWRT,*) ' MY T (BUFF) array'
csk    IF( SPLIT_IJKL )THEN
csk      CALL WRTMATMN8(T_BUFF,1,MY_T_LEN,1,MY_T_LEN,LUWRT)
csk    ELSE
csk      CALL WRTMATMN8(WORK(KT_CC),1,NALLINT,1,NALLINT,LUWRT)
csk    END IF
          ISHARE_DENSM = 0
C
C         FIXME: control ISHARE_DENSM from the outside
C
C         case 1: only MASTER will get all density matrix elements
C             --> ISHARE_DENSM = 0
C         case 2: density matrix elements on all CPUs
C             --> ISHARE_DENSM = 1
C
          IF( IT_SHL .lt. 0 )THEN
C
            L_MAX_DENS_BLK = 0
            L_MAX_DENS_BLK = NCALC_MAX_BLK(IRC)
C
            IF( NOOCCN_DENS_RUN .or. SPLIT_IJKL )THEN
              CALL MEMMAR_I8(KMAX_DENS_BLK,L_MAX_DENS_BLK,'ADDL  ',
     &                       1,'KDMBLK')
              CALL DZERO8(WORK(KMAX_DENS_BLK),L_MAX_DENS_BLK)
            END IF
C
            IF( NOOCCN_DENS_RUN ) THEN
              IF( SPLIT_IJKL ) THEN
                CALL COLLECT_DENS_MAT_SPLIT(WORK(KMAX_DENS_BLK),T_BUFF,
     &                                      IT_TTPL,IRC)
              ELSE 
              CALL COLLECT_DENS_MAT(WORK(KMAX_DENS_BLK),WORK(KT_CC),
     &                               IRC,ISHARE_DENSM)
              END IF
            ELSE
              IF( SPLIT_IJKL ) THEN
                CALL COLLECT_DENS_MAT_SPLIT(F2,T_BUFF,IT_TTPL,IRC)
              ELSE

                CALL COLLECT_DENS_MAT(F2,WORK(KT_CC),IRC,ISHARE_DENSM)
              END IF
            END IF
C           ^ NO run?
C         MPI one-sided based T split
          ELSE 
#if defined VAR_SHMEM
C
C           make sure each CPU gets here before starting T block 
C           collection
C
            CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
            IF( IT_SHL .eq. 0 )THEN
C
               MY_XT_BUFF_LEN = NCALC_MAX_BLK(IRC)
C
               CALL MPIXMEM_ALLOC(MY_XT_PTR,1,MY_XT_BUFF_LEN,
     &                            MPI_INFO_NULL,.FALSE.)
               CALL DZERO8(XT_BUFF,MY_XT_BUFF_LEN)
               CALL COLLECT_DENS_MAT_SM_LOC(T_BUFF,XT_BUFF,IT_TTPL,
     &                                      IT_TTOL,IRC,ISHARE_DENSM)
               CALL MPIXMEM_FREE(XT_BUFF)
C
            END IF
C
C           parallel MCSCF: provide WORK(KINDEX) to all CPU's
C
#endif /* VAR_SHMEM */
          END IF
C         ^ IT_SHL .lt. 0
C
#endif
C
          IPRDEN = 00
C
          IF (IPRDEN .ge. 10 ) THEN
            IF( MYPROC .eq. MASTER .and. ( IT_SHL .lt. 0 
     &          .and. .NOT. SPLIT_IJKL) ) THEN
              write(LUWRT,*) '======================================='
              write(LUWRT,'(A,I4)') ' Density matrix for root ',
     &                                JROOT_SCR
              write(LUWRT,*) '======================================='
C
              call wrtdens(WORK(KT_CC),WORK(KINDEX),WORK(KLIHIND),
     &                     WORK(KSIGNNHX),NTOOB2,WORK(KLABEXTP),
     &                     WORK(KLOP_REO),WORK(KSIOPREO),IRC)
            END IF
          END IF
C
          IPRDEN = 00
#if defined (VAR_MPI2)
C
C         parallel MCSCF: use all CPU's to transfer the density 
C         matrices...
C
          IF( SPLIT_IJKL )THEN
C
C           collect information about density T blocks and their 
C           "contributors". Choose process with the lowest tag 
C           (MYPROC) as sender for the master (receiver). 
C           IT_TTPL --> gather --> sort --> write back on IT_TTPL
C           previous information on IT_TTPL is lost!
C
C           NOTE: --> parallel MCSCF: receiver may be all processes 
C           or a subgroup depending on the MCSCF integral distribution.  
C           That is, why we already use MPI_ALLGATHER!
C
            MY_DL_TMP_LEN = NSPOBEX_TP*NMPROC*IRC_SAVE
C
            allocate(IT_DLIST_TMP(MY_DL_TMP_LEN))
            IT_DLIST_TMP = 0
            CALL MPI_ALLGATHER(IT_TTPL,NSPOBEX_TP*IRC_SAVE,
     &                         my_MPI_INTEGER,IT_DLIST_TMP,
     &                         NSPOBEX_TP*IRC_SAVE,my_MPI_INTEGER,
     &                         MPI_COMM_WORLD,IERR)
C
            CALL FIND_SENDER(IT_TTPL,IT_DLIST_TMP,NSPOBEX_TP)
C
csk         WRITE(LUWRT,*) ' AFTER SENDER matching IT_TTPL'
csk         CALL IWRTMAMN(IT_TTPL,1,NSPOBEX_TP*IRC,1,
csk  &                    NSPOBEX_TP*IRC,LUWRT)
C
C           release scratch memory
            deallocate(IT_DLIST_TMP)
          ELSE
            IF( MYPROC .ne. MASTER ) GOTO 919
          END IF
#endif
C
C         transfer density matrices to KRMC format
C         ----------------------------------------
C
C         local or global T distribution
          IF( IT_SHL .ge. 0 )THEN
#if defined (VAR_MPI2)
            IF( NOOCCN_DENS_RUN )THEN
C
              call trnsfdens_no_sm(F1,T_BUFF,IT_TTPL,IT_TTOL,
     &                             WORK(KINDEX),
     &                             WORK(KLABEXTP),WORK(KLIHIND),
     &                             WORK(KSIGNNHX),WORK(KLOP_REO),
     &                             WORK(KSIOPREO),IRC,NTOOB2,IPRDEN)
C
            ELSE
C
              CALL DZERO(F2,LF2_ZERO)
C
C             FIXME: WRITE ME!!!
C
csk           call trnsfdens_sm(F1,F2,T_BUFF,IT_TTPL,IT_TTOL,
csk     &                       WORK(KINDEX),
csk     &                       WORK(KLABEXTP),WORK(KLIHIND),
csk     &                       WORK(KSIGNNHX),WORK(KLOP_REO),
csk     &                       WORK(KSIOPREO),IRC,NTOOB2,IPRDEN)
              CALL Abend2( ' IMPLEMENT ME! ' )
            END IF
#endif
C         T split / no T distribution
          ELSE
            IF( NOOCCN_DENS_RUN )THEN
              IF( SPLIT_IJKL )THEN
#if defined (VAR_MPI2)
                call trnsfdens_no_split(F1,T_BUFF,WORK(KMAX_DENS_BLK),
     &                            IT_TTPL,WORK(KINDEX),WORK(KLABEXTP),
     &                            WORK(KLIHIND),
     &                            WORK(KSIGNNHX),WORK(KLOP_REO),
     &                            WORK(KSIOPREO),IRC,NTOOB2,IPRDEN)
#endif
              ELSE
                call trnsfdens_no(F1,WORK(KT_CC),WORK(KINDEX),
     &                            WORK(KLABEXTP),WORK(KLIHIND),
     &                            WORK(KSIGNNHX),WORK(KLOP_REO),
     &                            WORK(KSIOPREO),IRC,NTOOB2,IPRDEN)
              END IF
            ELSE
C             ... two-particle density matrix
              IF( IDEN12 .eq. 2 ) CALL DZERO(F2,LF2_ZERO)
C
              IF( SPLIT_IJKL )THEN
#if defined (VAR_MPI2)
                call trnsfdens_split(F1,F2,T_BUFF,WORK(KMAX_DENS_BLK),
     &                        IT_TTPL,
     &                        WORK(KINDEX),WORK(KLABEXTP),WORK(KLIHIND),
     &                        WORK(KSIGNNHX),WORK(KLOP_REO),
     &                        WORK(KSIOPREO),IRC,NTOOB2,IPRDEN)
#endif
              ELSE
csk             IF( IDEN12 .eq. 1 ) THEN
csk               write(LUWRT,*) ' '
csk               write(LUWRT,*) '    this is my INDEX array        '
csk               write(LUWRT,*) ' '
csk               CALL IWRTMAMN8(WORK(KINDEX),1,NTERMS,1,NTERMS,LUWRT)
csk             END IF
                call trnsfdens(F1,F2,WORK(KT_CC),
     &                        WORK(KINDEX),WORK(KLABEXTP),WORK(KLIHIND),
     &                        WORK(KSIGNNHX),WORK(KLOP_REO),
     &                        WORK(KSIOPREO),IRC,NTOOB2,IDEN12,IPRDEN)
              END IF
            END IF
          END IF
C
csk       WRITE(LUWRT,*) ' MY updated densities'
csk       WRITE(LUWRT,*) '           F1        '
csk       CALL WRTMATMN(F1,1,N1ELINT,1,N1ELINT,LUWRT)
C         IF (.NOT. NOOCCN_DENS_RUN )THEN
C           WRITE(LUWRT,*) '           F2        '
C           CALL WRTMATMN8(F2,1,N2ELINT,1,N2ELINT,LUWRT)
C         END IF
C
C
#if defined (VAR_MPI2)  
C
C        all slaves are waiting for the master who is transferring the
C        densities...
C
 919     CONTINUE
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#endif        
C
        else
          write(LUWRT,*) 'SIGDEN_ROUTE with unallowed value in gasci.'
          write(LUWRT,*) 'SIGDEN_ROUTE ',SIGDEN_ROUTE
        end if
C
      end if
C
*+++++++++++++++++++++++++++++++++++++++++++++++++++
      else
        WRITE(LUWRT,*) 'CIRUN has unallowed value :',CIRUN
        CALL QUIT(' *** ERROR in GASCI *** Unallowed CIRUN type.')
      end if
C     ^ CIRUN cases for MCSCF calls.
*+++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
#if defined (VAR_MPI2)
      IF( IT_SHL .ge. 0 ) THEN
#if defined VAR_SHMEM
C
C       de-allocate T array and destroy T memory window
        CALL MPIXMEM_FREE(IT_TTPL)
        CALL MPIXMEM_FREE(IT_TTOL)
        CALL REL_MEM_T_WIN(T_BUFF)
        IF ( CPIC_PASS_D ) MPIXMEM_FREE(T_BUFF_D)
#endif
      ELSE IF (SPLIT_IJKL ) THEN
        deallocate(IT_TTPL)
        deallocate(T_BUFF)
      END IF
#if !defined VAR_SHMEM
      IF( DIAGONAL_ALLOC ) deallocate(T_BUFF_D)
#endif
C
C     reset partitioning of CI vector to 'sequential' fashion
C     (*_MS2 arrays are in CSTATE )
      CALL Z_BLKFO_REL(IDC,NMS2VAL,
     &                 ISPC,ICSM,KPCLBT,KPCLEBT,KPCI1BT,KPCIBT,KPCBLTP,
     &                 NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,NBAT_MS2,
     &                 IBAT_MS2,IDUMMY,.FALSE.,0,0,IDUMMY)
C
#endif
      write(LUWRT,'(/2A/)') '  I have completed task  ',CIRUN
C     release memory
      IDUM = 0
      call memmar(KDUM,  IDUM,  'FLUSM ',IDUM,'GASCI ')
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE INPCTL_REL(LUIN,LUOUT)
*
* Read and check input for LUCIA
*
*.1 : Position input file after  line containing &LUCIA &END
C      CALL SLASK
*.2 : read input in
      CALL READIN_REL(LUIN,LUOUT)
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
#ifdef LUCIAREL_STANDALONE
      PROGRAM LUCIAREL
#else
      SUBROUTINE PAMLUCIAR(WRK,LWRK)
#endif
*
* L U C I A R E L
*
*
* CI for program for :FCI
*                     RASCI
*                     GAS GAS GAS GAS GAS GAS
*
* Written by Jeppe Olsen , winter of 1991
*                          GAS version in action summer of 95
*
* General relativistic formalism for 4-/2-component MCSCF
*                          Timo Fleig + Jeppe Olsen
*                          2000 - 2001
*
* Excitation operator formalism for sigma and density implemented.
*                          Timo Fleig + Jeppe Olsen
*                          2001 - 2002
*
* Implementation of relativistic double group CI,
*                          Timo Fleig + Jeppe Olsen
*                          1995 - 1999
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Parameters for dimensioning
#include "mxpdim.inc"
*.Memory
#include "wrkspc.inc"
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
#include "glbbas.inc"
#include "csm.inc"
#include "cintfo.inc"
#include "picasso.inc"
*
      common/glbref/ISPC,ISM,ITERGLB
*
*
#ifdef LUCIAREL_STANDALONE
*.
********************************************
********************************************
* Initialize BIGGY
      do I=1,MXPWRD,1
         WORK(I) = 0.0D0
      end do
#else
      DIMENSION WRK(LWRK)
C
C     Compute offset from WRK to WORK
C
      CALL COMPOFF(WRK,WORK,K_OFFSET)
C
      MXPWRD = LWRK
      KADD = MXPWRD
      call memmar(K_OFFSET,KADD,'INI   ',IDUMMY,'DUMMY')
C
#endif
C FIXME !!!
      write(6,*) 'CAUTION !'
      write(6,*) 'LUCIAREL as standalone main program needs attention.'
      call quit ('Standalone version not fully implemented.')
C FIXME !!!
********************************************
********************************************
      CALL QENTER('LUCIAR')
*. Assign diskunits
      CALL DISKUN_REL
*. Header
      CALL WRTHD_REL(LUOUT)
      WRITE(LUOUT,'(10X,A,I12)')
     &'Dimension of R*8 workspace : ', MXPWRD
* ======================================================
* 1 : .Read input,insert defaults and cross check input
* ======================================================
#ifdef LUCIAREL_STANDALONE
      CALL INPCTL_REL(LUIN,LUOUT)
#else
      OPEN(LUIN,FILE='DIRAC.INP')
      CALL INPCTL_REL(LUIN,LUOUT)
      CLOSE(LUIN, STATUS='KEEP')
#endif
*. Ratio beteeen real and integer word length
      CALL ZIRAT_REL
*. From shells to orbitals
      CALL ORBINF_REL(LUOUT,IPRORB)
*. Types of determinants to be included
      CALL DETTYP(IPRSTR)
      print*,'ISO1,ISO2 ',ISO1,ISO2
*. Number of string types
      CALL STRTYP_GAS_REL(IPRSTR)
*. Symmetry information
      CALL SYMINF_REL(IPRORB)
*. Number of integrals
      CALL INTDIM_REL(IPRTRA)
*. Static memory, initialization and allocation
#ifdef LUCIAREL_STANDALONE
      K_OFFSET = 1
      KADD = MXPWRD
      call memmar(K_OFFSET,KADD,'INI   ',IDUMMY,'DUMMY')
#endif
*. Allocate work space, general
      CALL ALLOC_REL
*. Internal string information (stored in WORK, bases in /STRBAS/)
      CALL STRINF_GAS_REL(WORK(1),IPRSTR)
*. Internal subspaces
      CALL LCISPC_REL(LCSBLK,IPRCIX)
*. Space for Hamiltonian in form suitable for SIGDEN
*. Type of Hamiltonian : 2 e Delta MK = 0 + all 1 e ops
      if (SIGDEN_ROUTE.eq.'SIGDEN') then
        if (ISO1.eq.0.and.ISO2.eq.0) then
          IHTYPE = 1
        else if (ISO1.eq.1.and.ISO2.eq.0) then
          IHTYPE = 2
        else if (ISO1.eq.1.and.ISO2.eq.1.and.ENVIRO(1:6).ne.'DIRAC ') 
     &       then
          IHTYPE = 3
        else if (ISO1.eq.1.and.ISO2.eq.1.and.ENVIRO(1:6).eq.'DIRAC ') 
     &       then
          IHTYPE = 4
        else
          write(6,*) 'Setting of ISO1,ISO2 not allowed.'
          write(6,*) 'ISO1,ISO2 : ',ISO1,ISO2
          call abend2('Quitting.')
        end if
        ISYM_T = 1
        NO_TTWO = 0
        CALL SET_HOP_DBG(IHTYPE,NO_TTWO,ISYM_T,IPRHAM_CI)
      END IF
*. Symmetry of reference
      IF(PNTGRP.GT.1.and.PNTGRP.lt.5)
     &           CALL MLSM_REL(IREFSM,IREFPA,IREFSM,'CI',1)
      IF(NOINT.EQ.1) THEN
        WRITE(6,*) ' End of calculation without integrals'
        Call Abend2( ' End of calculation without integrals' )
      END IF
*
*====================================================
*  Iterative natural spinor CI method
*  ----------------------------------
*   (Default is one iteration, so normal CI)
*====================================================
*
      do INATIT=1,NATITER,1
*
*  Generate symmetry info and the like
        if (INATIT.eq.1) call luc2pic(ISO1,IPRTRA)
*
*. Generate and store pointers for transformed integrals
*  (and integrals !)
        IF(NOINT.EQ.0.and.INATIT.eq.1) THEN
          CALL INTIM_REL(IPRTRA)
        ELSE
          WRITE(6,*)
     &      ' No integrals imported (or > 1st NS iteration) '
        END IF
*
*. Now call PICASSO
*. PICASSO will read the integrals (in initial run)
*.              read SPAO coefficients (in initial run)
*.              and do the transformation
*.              of all integrals
*
C
C       FIXME: SK-26-11-2007 implement use of T_BUFF for
C              LUCIARELs-standalone version
C
        call picasso(INTSPACE,ISO1,IPRTRA,
     &               NL2D,IUB2D,ISM2D,T_BUFF )
        write(6,*) 'Integrals transformed and/or distributed.'
*
*
* ======================
* CI in GAS spaces
* ======================
*
        ISTOSPC = 0
        IF(IRESTR.EQ.1) ISTOSPC = 1
*
*. Restarted Lambda calculations need special attention
*  since restart is realized in in second calculation
*  First calculation is used to establish H0
*
        IRESTR_ORIG=IRESTR
        IF(IRESTR.EQ.1.AND.XLAMBDA.NE.1.0D0) THEN
*
          write(6,*) 'THIS ROUTE IS STILL FORBIDDEN !!'
          write(6,'(A16,I2,1F3.1)')
     &          'IRESTR,XLAMBDA =',IRESTR,XLAMBDA
          Call Abend2( 'Quitting.' )
*
        END IF
*       ^ End of special handling of restarted calc with lambda modified op.
*
        ISKIPEI_INI = ISKIPEI
        DO JCMBSPC = 1, NCMBSPC
          ISKIPEI = 0
*
*. Perform GASCI/GASPT calculation for internal space ICISPC
*
          WRITE(6,*)
          WRITE(6,'(A)')
     &    ' =============================================='
          WRITE(6,'(A,I3)')
     &    ' Information on calculations in space ', JCMBSPC
          WRITE(6,'(A)')
     &    ' =============================================='
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &    ' Number of calculation in this CI space ', NSEQCI(JCMBSPC)
*
          IF(IRESTR.EQ.1.AND.XLAMBDA.NE.1.0D0.AND.JCMBSPC.EQ.1) THEN
            WRITE(6,*) ' Remember No restart in calc 1 (Lambda calc)'
            IRESTR = 0
            Call Abend2( 'Quitting.' )
          END IF
*
          I_EXPAND = 1
          IF(XLAMBDA.NE.1.0D0 .AND.JCMBSPC.GT.1) THEN
            WRITE(6,*) ' =================================='
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) ' =================================='
            IF(JCMBSPC.EQ.2) THEN
*. Modify the operator
              WRITE(6,*) ' Operator will be modified '
              WRITE(6,*) ' This route is still forbidden! '
              Call Abend2( 'Quitting.' )
CTF           CALL GENH1(XLAMBDA)
CTF           CALL SCLH2(XLAMBDA)
            END IF
          END IF
          IF(JCMBSPC.EQ.2.AND.IRESTR_ORIG.EQ.1
     &                   .AND.XLAMBDA.NE.1.0D0)  THEN
            WRITE(6,*) ' Restart vectors will be copied to LUC'
            WRITE(6,*) ' CI will restart with vectors from prev. calc'
            IF(ICISTR.EQ.1) THEN
              LBLK = XISPSM(IREFSM,2)
            ELSE
              LBLK = -1
            END IF
            call memmar(KDUM,IDUM,'MARK ',IDUM,'VCSAVE')
            IF(ICISTR.EQ.1) THEN
              LBLOCK = XISPSM(IREFSM,2)
            ELSE IF (ICISTR.EQ.2) THEN
              LBLOCK = MXSB
            ELSE IF (ICISTR.EQ.3) THEN
              LBLOCK = MXSOOB
            END IF
            call memmar(KVEC1,LBLOCK,'ADDS  ',2,'VEC1  ')
            CALL REWINE(LUC,-1)
            CALL REWINE(LUSC39,-1)
            DO JROOT = 1, NROOT
C             CALL COPVCD_REL(LUIN,LUOUT,SEGMNT,IREW,LBLK)
              CALL COPVCD_REL(LUSC39,LUC,WORK(KVEC1),0,LBLK)
            END DO
C?        WRITE(6,*) ' Back copied first vector on LUC'
C?        CALL WRTVCD(WORK(KVEC1),LUC,1,LBLK)
            call memmar(KDUM,IDUM,'FLUSM',IDUM,'VCSAVE')
*. No expansion should follow
            I_EXPAND = 0
*. But normal restart
            IRESTR = 1
          END IF
*
          DO JSEQ = 1,NSEQCI(JCMBSPC)
            IF((ISEQCI(JSEQ,JCMBSPC).GT.0
     &          .OR.ISEQCI(JSEQ,JCMBSPC).EQ.-5)
     &          .AND. JCMBSPC.NE.1.AND.ISTOSPC.NE.JCMBSPC
     &          .AND.I_EXPAND.EQ.1)                   THEN
*. Restart from previous spaces ( Assuming a progressing sequence :
*  spaces are just added, not subtracted )
*  ( Used for perturbation calculations and CI calculations )
*
              LUIN = LUC
              LUOUT = LUSC1
              IF(ICISTR.EQ.1) THEN
                LBLK = XISPSM(IREFSM,JCMBSPC)
              ELSE
                LBLK = -1
              END IF
              WRITE(6,*) ' LBLK = ', LBLK
              CALL EXPCIV_REL(IREFSM,ISTOSPC,LUIN,JCMBSPC,LUOUT,
     &                        LBLK,LUSC2,
     &                        NROOT,1,IDC,IPRDIA)
*. Last space where vectors were stored
              ISTOSPC = JCMBSPC
              ISKIPEI = ISKIPEI_INI
*. Expanded vector will be used as initial vector in the
*. zero space calculation. Tell next CI to restart from
*. CI vectors
              IRESTR = 1
            END IF
            IF(ISEQCI(JSEQ,JCMBSPC).GT.0) THEN
              IF(JSEQ.EQ.1.AND.JCMBSPC.EQ.2.AND.IRST2.EQ.0) THEN
*. No restart from previous vectors - IRST2 has been set to zero
                IRESTR = 0
                WRITE(6,*) ' No restart from previous vectors'
              END IF
*. Good old normal CI !!!!
*. do CI in space JCMBPSC
              MAXIT = ISEQCI(JSEQ,JCMBSPC)
*
              ISPC = JCMBSPC
              ISM = IREFSM
              ITERGLB = 1
*
*
*
              write(6,*) 'Standalone version of LUCIAREL not activated.'
              call abend('LUCIAREL')
C             CALL GASCI_REL(INTSPACE,IREFSM,JCMBSPC,IPRDIA,EREF)
*
*
*
*. Last space where vectors were stored
              ISTOSPC = JCMBSPC
CC            ELSE IF (ISEQCI(JSEQ,JCMBSPC).GT. -5 ) THEN
CC              WRITE(6,'(A,I3)') ' Vector free calculation at level ',
CC     &                            -ISEQCI(JSEQ,JCMBSPC)
CC              LEVEL =  -ISEQCI(JSEQ,JCMBSPC)
CC              IF(IC1DSC.LE.0) THEN
CC                LU1EFF = 0
CC              ELSE
CC                LU1EFF = LUHC
CC              END IF
CC              WRITE(6,*) ' LU1EFF in MAIN ', LU1EFF
CC              CALL DIRDIR(JCMBSPC-1,JCMBSPC,IREFSM,LUC,LEVEL,EREF,
CC     &                    LUSC2,LUSC3,LU1EFF)
CC            ELSE IF (ISEQCI(JSEQ,JCMBSPC).LE.-5) THEN
CC              WRITE(6,'(A)') ' Perturbation calculation '
CC              CALL PERTCTL(IREFSM,JCMBSPC,EREF)
*. also perturbation expansion of EKT ??
CC              IF(IPTEKT.EQ.1) THEN
CC                WRITE(6,*) ' Perturbation expansion of EKT '
CC                WRITE(6,*) ' ============================= '
C                    PTEKT(LUN,LU0,N,NEKT,LPTEKT,ISM,ISPC)
CC                CALL PTEKT(LUC,LUSC36,NPERT,IREFSM,JCMBSPC)
CC              END IF
            END IF
          END DO
*
        END DO
*       ^ End loop over CI spaces
*
        if (NOMOFL.eq.0.or.IFINMO.eq.6) THEN
*. Create set of spinors from densities
c         WRITE(6,*) ' I am going to call FINMO '
c         CALL FINMO_REL
        end if
      end do
*     ^ End loop over natural spinor iterations (default=1)
*
      CALL QEXIT('LUCIAR')
#ifdef LUCIAREL_STANDALONE
      Call Abend2( ' I am home from  the loops ' )
#endif
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE READIN_REL(LUIN,LUOUT)
*
*
* File is supposed to be positioned at first record of input
* The end of the input stream is identified by END OF INPUT
* Unless MOLCS is specified,
* All keywords are initiated by a point ., while comments are
* initiated by a *.
*
* The keywords are can broadly be divided into two types
*  1 : Keywords describing CI calculation to be carried out
*  2 : Keywords describing how CI optimization should be performed
*
*
* All input parameter concerning CI space are saved in /LUCIN1/
* All input concerning actual CI vectors are save in /CSTATE/
* All input paramters concerning run are saved in /CRUN/
*
* Since the keywords are read in from one pass over input file,
* the keywords must be in logical order.For example, the number
* of irreducible representations (irreps) must be give before
* the number of shells per irrep
*
********************************************
* Jeppe Olsen, Spring of 1991
*
*  Lots of additions and modifications.
*    Timo Fleig, 1995 - 2000
********************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
      PARAMETER(MXPLNC = 72 )
      CHARACTER*72 TITLEC
      CHARACTER*72 CARD
      CHARACTER*72 CARD1
      CHARACTER*72 LASTCARD
      COMMON/CTITLER/ TITLEC(3)
      CHARACTER*6 KEYWOR
      PARAMETER(MXPKW = 98)
      DIMENSION KEYWOR(MXPKW)
      DIMENSION ISETKW(MXPKW)
*
      DATA KEYWOR/'TITLE ','PNTGRP','NIRREP','INTSPC','EXTSPC',
     &     'NACTEL','INACT ','CORE  ','RAS1  ','RAS2  ',
     &     'RAS3  ','MXSCTP','SECOND','REFSPC','INTSEL',
     &     'MS2   ','MULTS ','IREFSM','ROOTS ','IDIAG ',
     &     'MAXIT ','EXPHAM','RESTRT','INTIMP','INCORE',
     &     'DELETE','MSCOMB','MLCOMB','IPRSTR','IPRCIX',
     &     'IPRORB','IPRDIA','MXCIV ','CISTOR','NOCSF ',
     &     'IPRXT ','NOINT ','DMPINT','RESDIM','CJKAIB',
     &     'INIREF','RESTRF','IPROCC','MOCAA ','MOCAB ',
     &     'ECORE ','PERTU ','APRREF','APRZER','GASSH ',
     &     'GASSPC','CMBSPC','CICONV','SEQUEN','EXTKOP',
     &     'MACHIN','C1DSC ','H0SPC ','H0FRM ','RFROOT',
     &     'H0EX  ','INIDEG','LAMBDA','LCSBLK','IPRDEN',
     &     'NOMOFL','ECHO  ','FINORB','E_THRE','C_THRE',
     &     'E_CONV','C_CONV','CLSSEL','DENSI ','PTEKT ',
     &     'H0ROOT','NORST2','SKIPEI','XYZSYM','PROPER',
     &     'DBGCI ','RELACI','SO-OP ','MK2REF','MK2DEL',
     &     'TRANS ','SP<-MO','SP<-RE','IPRTRA','IPRSIG',
     &     'COMHAM','ONEPAR','FOCKEX','INPORB','IPRWFC',
     &     'NATSPI','SIGDEN','IPRHAM'/
*
*
      INTEGER CITYP
*.Largest allowed number of allowed irreps for orbs
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "cstate.inc"
#include "crun.inc"
#include "cprnt.inc"
#include "cgas.inc"
#include "oper.inc"
#include "machine.inc"
#include "cands.inc"
*     for picasso
#include "input.inc"
*
#include "krmcluci_inf.h"
*
*./CECORE/
      COMMON/CECORER/ECORE
*. Flag for compatibility with normal MOLCAS input format
      MOLCS = 1
*
      NERROR = 0
      NWARN = 0
      EXTSPC = 0
      IECHO = 0
      LASTCARD = '      '
* Default is spin-orbit relativistic calculation:
* Default is Kramers symmetry active
*  (general modification independent of environment when
*   properties are implemented and spinor basis works)
      KRASYM = .TRUE.
      ISPINFREE = 0
*****************************************************************
*                                                               *
* Part 1 : Read in Keywords and perform some preliminary checks *
*                                                               *
*****************************************************************
*
*. Defaults for pointgroup and number of irreps must be set here
*. Default point group D2H
      PNTGRP = 1
*. Default number of irreps
      NIRREP = 8
      NSMOB = 8
*
      CALL ISETVC(ISETKW,0,MXPKW)
*
 1000 CONTINUE
*. Next potential keyword
      READ(LUIN,'(A)') CARD
*. Left-position nonblank characters in CARD
      CALL LFTPOS(CARD,MXPLNC)
*. Change to upper case
C          UPPCAS(LINE,LENGTH)
      CALL UPPCAS(CARD,MXPLNC)
C?      WRITE(6,'(A,A)') ' Upper case form ', CARD
      IF(CARD(1:1).EQ.'*'.OR.CARD(1:1).EQ.'#'.OR.
     &     CARD(1:1).EQ.'!') THEN
*. Skip comment cards
         GOTO 999
*. End of input card
      ELSE IF(CARD(1:5).EQ.'ENDOF') THEN
         GOTO 1001
      ELSE IF(MOLCS.EQ.0.AND.CARD(1:1).NE.'.') THEN
*. Line out of context
         WRITE(LUOUT,'(A)') ' Warning, card out of context : '
         WRITE(LUOUT,'(A)') CARD
         NWARN = NWARN + 1
      ELSE IF(MOLCS.EQ.1.OR.CARD(1:1).EQ.'.') THEN
         IF(MOLCS.EQ.1) THEN
*. Move characters one place to right
            DO 1286 ICHAR = 7,2,-1
               CARD(ICHAR:ICHAR) = CARD(ICHAR-1:ICHAR-1)
 1286       CONTINUE
            CARD(1:1) = ' '
         END IF
*. A keyword has been identified, match with possible keywords
         IF(CARD(2:6).EQ.'TITLE' ) THEN
*
* =========================
*.Keyword 1 :  TITLE cards
* =========================
*
*. Three title cards
            ISETKW(1) = 1
            DO 20 IC = 1, 3
               READ(LUIN,'(A)') TITLEC(IC)
 20         CONTINUE
            GOTO 999
         END IF
*
*
*================================================
*. Keyword 2 : <POINTG> : Point group of orbitals
*================================================
*
* Possible point groups : D2H,CINFV,DINFH,O3
         IF(CARD(2:4).EQ.'D2H'   .OR.
     &        CARD(2:6).EQ.'CINFV' .OR.
     &        CARD(2:6).EQ.'DINFH' .OR.
     &        CARD(2:3).EQ.'O3'    .OR.
     &        CARD(2:4).EQ.'DB-'   ) THEN
*
            ISETKW(2) = 1
            IF(CARD(2:4).EQ.'D2H') THEN
               PNTGRP = 1
            ELSE IF(CARD(2:6).EQ.'CINFV') THEN
               PNTGRP = 2
            ELSE IF(CARD(2:6).EQ.'DINFH') THEN
               PNTGRP = 3
            ELSE IF(CARD(2:3).EQ.'O3') THEN
               PNTGRP = 4
            ELSE IF(CARD(2:4).EQ.'DB-') THEN
*. Double group symmetry, of type ??
               PNTGRP = 5
*. What double group
*
* Point group         Code (PNTGRP)
*  C2H          =>          9
*  C2           =>          8
*  CS           =>          7
*  CI           =>          6
*  C1           =>          5
               IF(CARD(5:7).EQ.'C2H') THEN
                  PNTGRP = 9
                  NIRREP = 8
*. The individual spatial orbitals are the shells,
*. so the number of shell symmetries equal just half the number
*. of irreps
               ELSE IF(CARD(5:6).EQ.'C2') THEN
                  PNTGRP = 8
                  NIRREP = 4
               ELSE IF(CARD(5:6).EQ.'CS') THEN
                  PNTGRP = 7
                  NIRREP = 4
               ELSE IF(CARD(5:6).EQ.'CI') THEN
                  PNTGRP = 6
                  NIRREP = 4
               ELSE IF(CARD(5:6).EQ.'C1') THEN
                  PNTGRP = 5
                  NIRREP = 2
               ELSE
                  WRITE(6,*) ' Unknown double group '
                  WRITE(LUOUT,'(A)') CARD
                  ISETKW(2) = -1
                  NERROR = NERROR + 1
                  NIRREP = 0
               END IF
*. The individual spatial orbitals are the shells,
*. so the number of shell symmetries equal just half the number
*. of irreps
               NSMSH = NIRREP/2
*. The spinor symmetries will be NIRREP/2+1 to NIRREP,
*. so the number of spinor (orbital) symmetries are set as
               NSMOB = NIRREP
            END IF
            GOTO 999
         END IF
*
c        IF(CARD(2:7).EQ.'NIRREP') THEN
*
*=====================================================
*. Keyword 3 : <NIRREP> : Number of irreps of orbitals
*=====================================================
*
* Number of irreducible representations in point group
* D2h              : 1,2,4,8
* C inf H, D inf H : largest ML
* O3               : Largest L
* DBG              :  1,2,4,8 ( Dimension of all irreps)
*.D2h or subgroup
* ===============
c           IF(PNTGRP.EQ.1) THEN
c              READ(LUIN,*) NIRREP
c              NSMCMP = NIRREP
c              NSMOB  = NIRREP
c              NSMSH =  NIRREP
c
c              ISETKW(3) = 1
*.Dimensions 3,5,6,7,8 are not allowed
c              IF(NIRREP.EQ.3.OR.(NIRREP.GT.4.AND.NIRREP.LT.8)) THEN
c                 WRITE(LUOUT,*) ' Input error : NIRREP = ', NIRREP
c                 WRITE(LUOUT,*) ' Allowed values of NIRREP :1,2,4,8'
c                 NERROR = NERROR + 1
c                 ISETKW(3) = -1
c              END IF
*. Zero values used for other pointgroups
c              MAXML  = -1
c              MAXL   = -1
c             INVCNT = -1
c          ELSE IF (PNTGRP.EQ.2) THEN
*. Cinf V
* =======
c             READ(LUIN,*) MAXML
c             ISETKW (3) = 1
c             IF(MAXML.LT.0) THEN
c                WRITE(LUOUT,*)
c    &           ' Largest ML values of shells must be atleast be zero '
c                WRITE(LUOUT,*) ' MAXML from input :' ,MAXML
c                NERROR = NERROR + 1
c                ISETKW(3) = -1
c             END IF
c             NIRREP =  MAXML + 1
c             NSMCMP = 2 * MAXML + 1
c             NSMOB = NSMCMP
c             INVCNT = 0
c             MAXL = -1
c          ELSE IF (PNTGRP.EQ.3) THEN
*. Dinf H
* =======
c             READ(LUIN,*) MAXML
c             ISETKW (3) = 1
c             IF(MAXML.LT.0) THEN
c                WRITE(LUOUT,*)
c    &           ' Largest ML values of shells must be atleast be zero '
c                WRITE(LUOUT,*) ' MAXML from input : ',MAXML
c                NERROR = NERROR + 1
c                ISETKW(3) = -1
c             END IF
c             NIRREP = 2 * ( MAXML + 1)
c             NSMCMP = 2 * ( 2*MAXML + 1 )
c             NSMOB = NSMCMP
c             INVCNT = 1
c             MAXL  = -1
c          ELSE IF (PNTGRP.EQ.4) THEN
*. O 3
* =======
c             READ(LUIN,*) MAXL
c             ISETKW (3) = 1
c             IF(MAXL.LT.0) THEN
c                WRITE(LUOUT,*)
c    &            ' Largest L values of shells must be atleast be zero '
c                WRITE(LUOUT,*) ' MAXL from input : ' , MAXL
c                NERROR = NERROR + 1
c                ISETKW(3) = -1
c             END IF
c             MAXML = MAXL
c             NIRREP = MAXL + 1
c             NSMCOM = 2 * (2 *MAXML + 1 )
c             NSMOB = NSMCMP
c             INVCNT = 1
c          ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
*. Double group
* ===============
*. Well the dimension of the point group is given by the
*. specified double group but for general compatibility
*.we read the number and skip it
c             READ(LUIN,*) NIRREPX
c             WRITE(6,*) ' The supplied number of IRREPS'
c             WRITE(6,*) ' was read but is not used '
c             ISETKW(3) = 1
c          END IF
c          IF(ISETKW(3).EQ.-1)
c    &          WRITE(LUOUT,*) ' .NIRREP input incorrect !! . '
c          GOTO 999
c       END IF
*
* ================================================
*. Keyword 4 : INTSPC : Type of internal CI space
* ================================================
*
        IF(CARD(2:4).EQ.'CAS'.OR.CARD(2:4).EQ.'FCI'.OR.
     &       (CARD(2:4).EQ.'RAS'.AND.CARD(5:5).EQ.' ')) THEN
           ISETKW(4) = 1
           IF(CARD(2:4).EQ.'CAS'.OR. CARD(2:4).EQ.'FCI' ) THEN
              INTSPC = 1
           ELSE IF (CARD(2:4).EQ.'RAS') THEN
              INTSPC = 2
*. Limits on allowed number of electrons in RASI and RAS III
              READ(LUIN,*) MNRS10,MXRS30
           END IF
           GOTO 999
        END IF
*
*===========================
* Keyword 5 : External space
*===========================
*
        IF(CARD(2:7).EQ.'EXTSPC') THEN
           ISETKW(5) = 1
           READ(LUIN,'(A)') CARD1
           CALL LFTPOS(CARD1,MXPLNC)
           CALL UPPCAS(CARD1,MXPLNC)
           IF(CARD1(1:4).EQ.'NONE') THEN
              EXTSPC = 0
              MXER4 = 0
              MXHR0 = 0
           ELSE IF
     &        (CARD1(1:4).EQ.'CORE'.AND.CARD1(5:10).EQ.'SECOND') THEN
              EXTSPC = 3
              READ(LUIN,*) MXHR0,MXER4
           ELSE IF
     &        (CARD1(1:6).EQ.'SECOND'.AND.CARD1(7:10).EQ.'CORE') THEN
              EXTSPC = 3
              READ(LUIN,*) MXER4,MXHR0
           ELSE IF(CARD1(1:4).EQ.'CORE') THEN
              EXTSPC = 1
              READ(LUIN,*) MXHR0
           ELSE IF(CARD1(1:6).EQ.'SECOND') THEN
              EXTSPC = 2
              READ(LUIN,*) MXER4
           ELSE
              ISETKW(5) = - 1
              WRITE(6,*) ' Illegal card for EXTSPC : '
              WRITE(6,'(A)') CARD1
              NERROR = NERROR + 1
           END IF
           GOTO 999
        END IF
*
* =============================================
* Keyword 6 NACTEL : Number of active electrons
* =============================================
*
        IF(CARD(2:7).EQ.'NACTEL') THEN
           READ(LUIN,*) NACTEL
           ISETKW(6) = 1
           IF(NACTEL.LT.0) THEN
              WRITE(LUOUT,*)
     &          ' ERROR : Illegal number of active electrons ', NACTEL
              ISETKW(6) = -1
              NERROR = NERROR + 1
           end if
           if (NACTEL.eq.1) then
              IELCI = 1
           else
              IELCI = 2
           END IF
           GOTO 999
        END IF
*==================
* 7 : Inactive shells
*==================
         IF(CARD(2:7).EQ.'INACTI') THEN
            READ(LUIN,*) (NINASH(IRREP),IRREP = 1, NIRREP)
            ISETKW(7) = 1
            GOTO 999
         END IF
*=================================
* 8 : Core shells ( = RAS0 shells)
*==================================
         IF(CARD(2:5).EQ.'CORE') THEN
            READ(LUIN,*) (NRS0SH(1,IRREP),IRREP = 1, NIRREP)
            ISETKW(8) = 1
            EXTSPC = EXTSPC + 1
            GOTO 999
         END IF
*===========
* 9 : RAS 1
*===========
         IF(CARD(2:5).EQ.'RAS1') THEN
*.Number of RAS 1 shells per irrep
            READ(LUIN,*) (NRSSH(IRREP,1),IRREP = 1, NIRREP)
*.Smallest allowed number of electrons in RAS 1
C!         READ(LUIN,*) MNER10
            ISETKW(9) = 1
            GOTO 999
         END IF
*===========
* 10 : RAS 2
*===========
         IF(CARD(2:5).EQ.'RAS2'.OR.CARD(2:7).EQ.'ACTIVE') THEN
            READ(LUIN,*) (NRSSH(IRREP,2),IRREP = 1, NIRREP)
            ISETKW(10) = 1
            GOTO 999
         END IF
*===========
* 11 : RAS 3
*===========
         IF(CARD(2:5).EQ.'RAS3') THEN
            ISETKW(11) = 1
*.Number of RAS 3 shells per irrep
            READ(LUIN,*) (NRSSH(IRREP,3),IRREP = 1, NIRREP)
*.Largest allowed number of electrons in RAS III
C!         READ(LUIN,*) MXER30
            GOTO 999
         END IF
* ==================================================
* 12 : RAS 4, subdivision of shells of given irrep
* ==================================================
*. Transferred to read in of RAS1
* =================================================
* 13 : Number of shells in secondary space per type
* =================================================
         IF(CARD(2:7).EQ.'SECOND') THEN
*. Number of subdivisions of RAS4
            READ(LUIN,*) MXR4TP
*. Number of secondary irreps per type
            DO 30 ITYPE = 1, MXR4TP
               READ(LUIN,*) (NRS4SH(IRREP,ITYPE),IRREP = 1, NIRREP)
 30         CONTINUE
            ISETKW(13) = 1
*.Largest allowed number of electron in secondary space
            GOTO 999
         END IF
* =========================
* 14 : Reference space
* =========================
*
* For RAS reference spaces,the reference space can be reduced compared
* to the zero order space by introducing different constraints on
* the number of allowed electrons in RAS1 and RAS3. The reference space
* must however be a subset of the zero order space.
*
         IF(CARD(2:7).EQ.'REFSPC') THEN
            READ(LUIN,*) MNRS1R,MXRS3R
            ISETKW(14) = 1
            IF(MNRS1R.LT.MNRS10.OR.MXRS3R.GT.MXRS30) THEN
               WRITE(LUOUT,*)
     &              ' Reference space larger than zero order space'
               NERROR = NERROR + 1
               ISETKW(14) = -1
            END IF
            GOTO 999
         END IF
*
* =========================================================
* 15 : selection of active configurations in internal space
* =========================================================
*
         IF(CARD(2:7).EQ.'INTSEL' ) THEN
            ISETKW(15) = 1
*
            READ(LUIN,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
            CALL UPPCAS(CARD1,MXPLNC)
*
            IF(CARD1(1:4).EQ.'NONE') THEN
*. All internals are included
               INTSEL = 0
            ELSE IF(CARD1(1:6).EQ.'INDTST') THEN
*. Include coeffcients larger than CTHRES or having energy contributions
*. larger than ETHREA
               INTSEL = 1
               READ(LUIN,*)  CTHRES,ETHRES
            ELSE IF(CARD1(1:6).EQ.'TOTTST') THEN
*. Obtain CTHRES of the total wavefinction and ETHRES of the total
*. energy
               INTSEL = 2
               READ(LUIN,*)  CTHRES,ETHRES
            ELSE IF(CARD1(1:6).EQ.'INDWCN') THEN
*. Include all configutations with reference weights larger than
*. a given threshold in reference CI
               INTSEL = 3
               READ(LUIN,*) XWCNF
            ELSE IF(CARD1(1:6).EQ.'TOTWCN') THEN
*. Include the largest configurations so all a given fraction
*. of the Zero order reference is included in the CI
               INTSEL = 4
               READ(LUIN,*) XWCNF
            ELSE
               WRITE(LUOUT,*) ' Error : Illegal card in INTSEL :'
               WRITE(LUOUT,'(A)') CARD1
               ISETKW(15) = - 1
               NERROR = NERROR + 1
            END IF
            GOTO 999
         END IF
*===============================
* 16 : Two times spin projection
* ==============================
         IF(CARD(2:4).EQ.'MS2') THEN
            ISETKW(16) = 1
            READ(LUIN,*) MS2
            GOTO 999
         END IF
*========================
* 17 : spin multiplicity
* =======================
c        IF(CARD(2:6).EQ.'MULTS') THEN
c           ISETKW(17) = 1
c           READ(LUIN,*) MULTS
c           GOTO 999
c        END IF
*========================
* 18 : Reference symmetry
* =======================
         IF(CARD(2:7).EQ.'IREFSM') THEN
            ISETKW(18) = 1
            IF(PNTGRP.EQ.1) THEN
               READ(LUIN,*) IREFSM
            ELSE IF(PNTGRP.EQ.2) THEN
               READ(LUIN,*) IREFML
            ELSE IF (PNTGRP.EQ.3) THEN
               READ(LUIN,*) IREFML,IREFPA
               IF(IREFPA.EQ.-1) IREFPA = 2
            ELSE IF (PNTGRP.EQ.4) THEN
               READ(LUIN,*) IREFL,IREFML,IREFPA
               IF(IREFPA.EQ.-1) IREFPA = 2
            ELSE IF (PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
               READ(LUIN,*) IREFSM
            END IF
            GOTO 999
         END IF
*==========================
* 19 : Roots to be obtained
* =========================
         IF(CARD(2:6).EQ.'ROOTS') THEN
            ISETKW(19) = 1
            READ(LUIN,*) NROOT
            do I=1,NROOT,1
               IROOT(I) = I
            end do
            GOTO 999
         END IF
*===============================
* 20 : Diagonalization algorithm   : .MEGACI , .TERACI
*===============================
         IF(CARD(2:7).EQ.'MEGACI') THEN
            ISETKW(20) = 1
            IDIAG = 1
            GOTO 999
         ELSE IF(CARD(2:7).EQ.'TERACI') THEN
            ISETKW(20) = 1
            IDIAG = 2
            GOTO 999
         END IF
*==================================
* 21 : Explicit hamilton matrix   : MXP1,MXP2,MXQ
*==================================
         IF(CARD(2:7).EQ.'EXPHAM') THEN
            READ(LUIN,*) MXP1,MXP2,MXQ
            ISETKW(21) = 1
            GOTO 999
         END IF
*===================================================
* 22 : Largest allowed number of Iterations per root : MAXIT
*===================================================
         IF(CARD(2:6).EQ.'MAXIT') THEN
            ISETKW(22) = 1
            READ(LUIN,*) MAXIT
            GOTO 999
         END IF
*====================
* 23 : Restart option
*====================
         IF(CARD(2:7).EQ.'RESTRT') THEN
            ISETKW(23) = 1
            IRESTR = 1
            GOTO 999
         END IF
*========================================
* 24 Import of integrals and environment
*========================================
         IF(CARD(2:7).EQ.'MOLCAS'.OR.CARD(2:7).EQ.'ENV=MO') THEN
*. Integrals imported from MOLCAS
            INTIMP = 1
            ISETKW(24) = 1
            ENVIRO(1:6) = 'MOLCAS'
            GOTO 999
         ELSE IF(CARD(2:6).EQ.'LUCAS')THEN
*. Integrals imported from LUCAS
            INTIMP = 2
            ISETKW(24) = 1
            ENVIRO(1:6) = 'LUCAS '
            GOTO 999
         ELSE IF(CARD(2:7).EQ.'FMINSM'.OR.
     &           CARD(2:7).EQ.'ENV=FM') THEN
*. Integrals read formatted in, only integrals differing from
*. zero by symmetry are  included
            INTIMP = 3
            ISETKW(24) = 1
            ENVIRO(1:6) = 'LUCIA '
            GOTO 999
         ELSE IF(CARD(2:7).EQ.'SIRIUS'.OR.CARD(2:7).EQ.'DALTON'
     &           .OR.CARD(2:7).EQ.'ENV=DA') THEN
*. Integrals imported from SIRIUS/DALTON
            INTIMP = 5
C         write(6,*) ' Sirius Flag activated '
            ISETKW(24) = 1
            ENVIRO(1:6) = 'DALTON'
            GOTO 999
         else if (CARD(2:7).eq.'DIRAC ') then
           INTIMP = 6
           ISETKW(24) = 1
           ENVIRO(1:6) = 'DIRAC '
           goto 999
*
* Integral import variable 7 is reserved for DIRAC MCSCF runs !!
*
         ELSE IF(CARD(2:7).EQ.'ENV=NO' ) THEN
*. No program environment, integrals, coefs will just be set to zero
            ENVIRO(1:6) = 'NONE  '
            ISETKW(24) = 1
            GOTO 999
Ctimo Added MOLECULE SWEDEN information
         else if(CARD(2:7).eq.'MOLECU') then
*. Integrals imported from MOLECULE
            INTIMP = 8
            ISETKW(24) = 1
            ENVIRO(1:6) = 'MOLECU'
            goto 999
         END IF
* ===============================
* 25 :INCORE option for integrals
* ==============================
         IF(CARD(2:7).EQ.'INCORE') THEN
            ISETKW(25) = 1
            INCORE = 1
            GOTO 999
         END IF
* ===================
* 26 : Deleted shells
* ===================
         IF(CARD(2:7).EQ.'DELETE') THEN
            ISETKW(26) = 1
            READ(LUIN,*) (NDELSH(IRREP),IRREP= 1, NIRREP)
            GOTO 999
         END IF
* ===================
* 27: Ms combinations
* ===================
         IF(CARD(2:7).EQ.'MSCOMB') THEN
            ISETKW(27) = 1
            READ(LUIN,*) PSSIGN
            IF(.NOT.(PSSIGN.EQ.-1.0D0.OR.PSSIGN.EQ.1.0D0)) THEN
              WRITE(LUOUT,*)' Illegal Spin combination factor ',PSSIGN
              NERROR = NERROR + 1
            END IF
            GOTO 999
         END IF
* ===================
* 28: Ml combinations
* ===================
         IF(CARD(2:7).EQ.'MLCOMB') THEN
            ISETKW(28) = 1
            READ(LUIN,*) PLSIGN
            IF(.NOT.(PLSIGN.EQ.-1.0D0.OR.PLSIGN.EQ.1.0D0)) THEN
               WRITE(LUOUT,*)' Illegal ml combination factor ',PLSIGN
               NERROR = NERROR + 1
            END IF
            GOTO 999
         END IF
* ======================================
* 29 : Print flag for string information
* ======================================
         IF(CARD(2:7).EQ.'IPRSTR') THEN
            ISETKW(29) = 1
            READ(LUIN,*) IPRSTR
            GOTO 999
         END IF
* ========================================
* 30 : Print flag for CI space information
* ========================================
         IF(CARD(2:7).EQ.'IPRCIX') THEN
            ISETKW(30) = 1
            READ(LUIN,*) IPRCIX
            GOTO 999
         END IF
* =======================================
* 31 : Print flag for orbital information
* =======================================
         IF(CARD(2:7).EQ.'IPRORB') THEN
            ISETKW(31) = 1
            READ(LUIN,*) IPRORB
            GOTO 999
         END IF
* ===============================================
* 32 : Print flag for diagonalization information
* ===============================================
         IF(CARD(2:7).EQ.'IPRDIA') THEN
            ISETKW(32) = 1
            READ(LUIN,*) IPRDIA
            GOTO 999
         END IF
* ===============================================
* 36 : Print flag for Externals
* ===============================================
         IF(CARD(2:6).EQ.'IPRXT') THEN
            ISETKW(36) = 1
            READ(LUIN,*) IPRXT
            GOTO 999
         END IF
* =====================================
* 43 : Print occupation of lowest Dets
* =====================================
         IF(CARD(2:7).EQ.'IPROCC') THEN
            ISETKW(43) = 1
            READ(LUIN,*) IPROCC
            GOTO 999
         END IF
* ====================================
* 65 : Print level for density matrices
* ====================================
         IF(CARD(2:7).EQ.'IPRDEN') THEN
            ISETKW(65) = 1
            READ(LUIN,*) IPRDEN
            GOTO 999
         END IF
*=========================================================
* 33 : Largest allowed number of Vectors in diagonalization
*=========================================================
         IF(CARD(2:6).EQ.'MXCIV') THEN
            ISETKW(33) = 1
            READ(LUIN,*) MXCIV_CI
            GOTO 999
         END IF
* =============================
* 34 : Storage mode for vectors
* =============================
*  Default: 3 symm-occ blocks held in core:
         ICISTR = 3
*
c        IF(CARD(2:7).EQ.'CISTOR')THEN
c           ISETKW(34) = 1
c           READ(LUIN,*) ICISTR
c           GOTO 999
c        END IF
* ================================
* 35 : Do not employ CSF expansion
* ================================
         IF(CARD(2:6).EQ.'NOCSF') THEN
            ISETKW(35) = 1
            NOCSF = 1
            GOTO 999
         END IF
* ================================
* 37 : Do not read in integrals
* ================================
         IF(CARD(2:6).EQ.'NOINT') THEN
            ISETKW(37) = 1
            NOINT = 1
            GOTO 999
         END IF
* ================================
* 38 : Dump integrals in formatted form on file 90
* ================================
         IF(CARD(2:7).EQ.'DMPINT') THEN
            ISETKW(38) = 1
            IDMPIN = 1
            GOTO 999
         END IF
* ================================
* 39 : Define dimension of resolution matrices
* ================================
         IF(CARD(2:7).EQ.'RESDIM') THEN
            ISETKW(39) = 1
            READ(LUIN,*) MXINKA
            GOTO 999
         END IF
* ====================================================================
* 40 : Use CJKAIB matrices as intermediate matrices in alpha-beta-loop
* ====================================================================
         IF(CARD(2:7).EQ.'CJKAIB') THEN
            ISETKW(40) = 1
            ICJKAIB = 1
            GOTO 999
         END IF
* ====================================================================
* 44 : Use Minimal operatioon count method for alpha-alpha and beta-beta
* ====================================================================
         IF(CARD(2:6).EQ.'MOCAA') THEN
            ISETKW(44) = 1
            MOCAA = 1
            GOTO 999
         END IF
* ====================================================================
* 45 : Use Minimal operation count method for alpha-beta
* ====================================================================
         IF(CARD(2:6).EQ.'MOCAB') THEN
            ISETKW(45) = 1
            MOCAB = 1
            GOTO 999
         END IF
* ====================================================================
* 41 : Initial CI in reference space
* ====================================================================
         IF(CARD(2:7).EQ.'INIREF') THEN
            ISETKW(41) = 1
            INIREF  = 1
            GOTO 999
         END IF
* ====================================================================
* 42 : Restart from reference CI expansion
* ====================================================================
         IF(CARD(2:7).EQ.'RESTRF') THEN
            ISETKW(42) = 1
            IRESTRF = 1
*. Flag that restart will be used for zero space calculation
            ISETKW(23) = 1
            IRESTR = 1
            GOTO 999
         END IF
* ====================================================================
* 46 : Read in of core energy
* ====================================================================
         IF(CARD(2:6).EQ.'ECORE') THEN
            ISETKW(46) = 1
            READ(LUIN,*) ECORE
            GOTO 999
         END IF
*
* =====================================================================
* 47 : Use Perturbation theory for zero order space
* =====================================================================
*
         IF(CARD(2:6).EQ.'PERTU') THEN
*
*. Perturbation theory : Three parameters to be specified :
*
*      1 : Max order of correction vectors required
*      2 : Type of partitioning ( H0 )
*          Current choices : MP, EN, H0READ
*      3 : zero order energy : E0=EX ( use exact energy of reference state )
*                              E0=AV ( Use expectation value of H0 )
*                              E0=RE ; Readin zero order energy in
            ISETKW(47) = 1
            IPERT = 1
*. Number of correction vectors
            READ(LUIN,*) NPERT
*. Moeller-Plesset or Epstein-Nesbet partitioning
            READ(LUIN,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
            CALL UPPCAS(CARD1,MXPLNC)
C?       WRITE(6,'(A)') CARD1
*
            IF(CARD1(1:2) .EQ. 'MP' ) THEN
               MPORENP = 1
               IPART = 1
            ELSE  IF(CARD1(1:2) .EQ. 'EN' ) THEN
               MPORENP = 2
               IPART = 2
            ELSE IF(CARD1(1:6).EQ.'H0READ' ) THEN
*. Read in one body hamiltonian
               MPORENP = 0
               IPART = 3
*.
            ELSE
               WRITE(LUOUT,*) ' Unknown partitioning '
               WRITE(LUOUT,'(A)') CARD1
               NERROR = NERROR + 1
            END IF
* Zero order energy :
            READ(LUIN,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
            CALL UPPCAS(CARD1,MXPLNC)
C?       WRITE(6,'(A)') CARD1
*
            IF(CARD1(1:5).EQ.'E0=AV') THEN
               IE0AVEX = 1
            ELSE IF(CARD1(1:5).EQ.'E0=EX') THEN
               IE0AVEX = 2
            ELSE IF(CARD1(1:5).EQ.'E0=RE') THEN
               IE0AVEX = 3
               READ(5,*) E0READ
               WRITE(6,*) ' Zero order energy =',E0READ
            ELSE
               WRITE(6,*) ' Unknown form of zero order energy '
               WRITE(LUOUT,'(A)') CARD1
               NERROR = NERROR + 1
            END IF
*
            GOTO 999
         END IF
*
* =====================================================================
* 48 : Approximate Hamiltonian in reference space
* =====================================================================
*
         IF(CARD(2:7).EQ.'APRREF') THEN
            ISETKW(48) = 1
            READ(LUIN,*)  MNRS1RE,MXRS3RE
*. Moeller-Plesset or Epstein-Nesbet partitioning
            READ(LUIN,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
*. Change to upper case
            CALL UPPCAS(CARD1,MXPLNC)
*
            IF(CARD1(1:2) .EQ. 'MP' ) THEN
               MPORENR = 1
            ELSE  IF(CARD1(1:2) .EQ. 'EN' ) THEN
               MPORENR = 2
            ELSE
               WRITE(LUOUT) ' Unknown partitioning '
               WRITE(LUOUT,'(A)') CARD1
               NERROR = NERROR + 1
            END IF
*
            IAPRREF = 1
            GOTO 999
         END IF
*
* =====================================================================
* 49 : Approximate Hamiltonian in zero order space
* =====================================================================
*
         IF(CARD(2:7).EQ.'APRZER') THEN
            ISETKW(49) = 1
            READ(LUIN,*)  MNRS1ZE,MXRS3ZE
*. Moeller-Plesset or Epstein-Nesbet partitioning
            READ(LUIN,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
*. Change to upper case
            CALL UPPCAS(CARD1,MXPLNC)
*
            IF(CARD1(1:2) .EQ. 'MP' ) THEN
               MPORENZ = 1
            ELSE  IF(CARD1(1:2) .EQ. 'EN' ) THEN
               MPORENZ = 2
            ELSE
               WRITE(LUOUT) ' Unknown partitioning '
               WRITE(LUOUT,'(A)') CARD1
               NERROR = NERROR + 1
            END IF
            IAPRZER = 1
            GOTO 999
         END IF
*
* =====================================================================
* 50 : Generalized active space concept invoked, orbital spaces
* =====================================================================
*
         IF(CARD(2:6).EQ.'GASSH') THEN
            Write(6,*) ' READIN, GASSH, NIRREP = ', NIRREP
*. A note when double group symmetry is being used.
*. The number of irreps in the double group is twice
* the number of irreps in the spatial point group
* (Dear Jeppe, this is only true for the quaternionic and
*  the complex matrix groups !! Timo)
* We read symmetry info for the orbitals, not the
* resulting spinors.
*. Generalized active space in use
            IDOGAS = 1
            call isetvc(NGSSH,0,MXPIRR*MXPNGAS)
            READ(LUIN,*) NGAS
            DO IGAS = 1, NGAS
*. Number of irreps per GAS
               READ(LUIN,*) (NGSSH(IRREP,IGAS),IRREP = 1, NSMSH)
            END DO
            ISETKW(50) = 1
            GOTO 999
         END IF
*
* =====================================================================
* 51 : Generalized active space occupation restrictions
* =====================================================================
*
         IF(CARD(2:7).EQ.'GASSPC') THEN
*. Orbital constraints in gas spaces
*. GASSH must have been defined vefore, check this
            IF(ISETKW(50).EQ.0) THEN
               WRITE(6,*) ' Dear User'
               WRITE(6,*)
               WRITE(6,*) ' GASSH must be specified before GASSPC'
               WRITE(6,*)
     &          ' Else I do not know about the number of orbital spaces'
               WRITE(6,*) ' So I will stop '
               Call Abend2( 'READIN: put GASSH before GASSPC' )
            END IF
            IDOGAS = 1
            ISETKW(51) = 1
*. Number of oribtal spaces
            READ(LUIN,*) NCISPC
            DO ISPC = 1, NCISPC
*. Upper and lower limits for each orbital space
               DO IGAS = 1, NGAS
                  READ(LUIN,*)  IGSOCCX(IGAS,1,ISPC),
     &                 IGSOCCX(IGAS,2,ISPC)
               END DO
            END DO
            GOTO 999
         END IF
*
* =====================================================================
* 52 : Calculations will be performed in combination of different GAS spaces
* =====================================================================
*
         IF(CARD(2:7).EQ.'CMBSPC') THEN
*. Calculations will be performed in combination of different GAS spaces
            IDOGAS = 1
            ISETKW(52) = 1
*. Check if SEQUEN have been specified.
            IF(ISETKW(54).EQ.1) THEN
             WRITE(6,*) ' Dear user '
             WRITE(6,*)
             WRITE(6,*)' SEQUEN flag has been specified before CMBSPC'
             WRITE(6,*)' This confuses me and makes me wonder what the'
             WRITE(6,*)' meaning of everything is. '
             WRITE(6,*)' Please ensure that CMBSPC is given before'
             write(6,*)' SEQUEN'
             WRITE(6,*)
             WRITE(6,*)'                                  Lucia  '
             WRITE(6,*)
             Call Abend2( 'READIN: Specify CMBSPC before SEQUEN' )
            END IF
*. Number of combination spaces
            READ(LUIN,*) NCMBSPC
            DO JCMBSPC  = 1, NCMBSPC
*. Number of gas spaces in this space
              READ(LUIN,*) LLCMBSPC
              LCMBSPC(JCMBSPC) = LLCMBSPC
*. Gasspaces included
              READ(LUIN,*) (ICMBSPC(IGASSPC,JCMBSPC),IGASSPC=1,LLCMBSPC)
            END DO
            GOTO 999
         END IF
*
* ===============================================
* 53 : Energy convergence of CI
* ===============================================
*
         IF(CARD(2:7).EQ.'CICONV') THEN
*. Energy convergence of CI
            READ(LUIN,*) THRES_E
            ISETKW(53) = 1
            GOTO 999
         END IF
*
* =========================================
* 54 : SEQUEN KEYWORD
* =========================================
*
         IF(CARD(2:7).EQ.'SEQUEN') THEN
            IF(ISETKW(52).EQ.0) THEN
*. Combination spaces were not explicitly defined,
*. assume each gas space is a comb space
              NCMBSPC = NCISPC
            END IF
            DO JCMBSPC = 1, NCMBSPC
             READ(LUIN,*) NSEQCI(JCMBSPC)
             IF(NSEQCI(JCMBSPC).GT.0)
     &       READ(LUIN,*) (ISEQCI(ICI,JCMBSPC),ICI = 1, NSEQCI(JCMBSPC))
            END DO
            ISETKW(54) = 1
            GOTO 999
         END IF
*
* =====================================================================
* 55 : Call EXTENDED KOOPMANS' THEOREM ROUTINE
* =====================================================================
*
         IF(CARD(2:7).EQ.'EXTKOP') THEN
*. Oh yes, we will do it !
            IEXTKOP = 1
            ISETKW(55) = 1
            GOTO 999
         END IF
*
* ==========================
* 56 : What's your engine ?
* ==========================
*
*  Default setting for machine:
         ARCHI = 'SGO'
*
c        IF(CARD(2:7).EQ.'MACHIN') THEN
c           ISETKW(56) = 1
c           READ(LUIN,'(A)') CARD1
c           CALL LFTPOS(CARD1,MXPLNC)
*. Change to upper case
c           CALL UPPCAS(CARD1,MXPLNC)
c           MACHINE(1:6) = CARD1(1:6)
c           WRITE(6,'(A,A)') ' Machine = ', MACHINE
c           if (MACHINE(1:3).eq.'SGI') then
c              ARCHI = 'SGO'
c           else if (MACHINE(1:3).eq.'IBM') then
c              ARCHI = 'RS6'
c           else if (MACHINE(1:2).eq.'HP') then
c              ARCHI = 'SGO'
c           end if
c           GOTO 999
c        END IF
*
* ==========================================================
* 57 : Save first order correction to wavefunction on DISC?
* ==========================================================
*
* ( For vector free calculations )
*
         IF(CARD(2:6).EQ.'C1DSC') THEN
            ISETKW(57) = 1
            IC1DSC = 1
            GOTO 999
         END IF
*
* ==========================================================
*.58 :  Specify subspaces in which perturbation is nonvanishing
* ==========================================================
*
         IF(CARD(2:6).EQ.'H0SPC') THEN
*. Ensure that number of GASSPACES have been defined
            IF(ISETKW(50).EQ.0) THEN
               WRITE(6,*) ' Dear User'
               WRITE(6,*)
               WRITE(6,*) ' GASSH must be specified before H0SPC'
               WRITE(6,*)
     &          ' Else I do not know about the number of orbital spaces'
               WRITE(6,*) ' So I will stop '
               Call Abend2( 'READIN: put GASSH before H0SPC ' )
            END IF
            READ(LUIN,*) NPTSPC
            IF(NPTSPC.GT.MXPPTSPC) THEN
*
               WRITE(LUOUT,*) ' To many perturbation spaces '
               WRITE(LUOUT,*)
     &              ' raise MXPPTSPC from ', MXPPTSPC ,' to ',NPTSPC
               Call Abend2( 'NPTSPC>MXPPTSPC in READIN ' )
            END IF
*
            IH0SPC = 1
            DO JPTSPC = 1, NPTSPC
*. Number of occupation spaces in this subspace
               DO JGAS = 1, NGAS
                  READ(LUIN,*)
     &                IOCPTSPC(1,JGAS,JPTSPC),IOCPTSPC(2,JGAS,JPTSPC)
               END DO
            END DO
            ISETKW(58) = 2
            GOTO 999
         END IF
*
* ============================================
*.59 :  Specify Type of H0 for each subspace
* ============================================
*
         IF(CARD(2:6).EQ.'H0FRM') THEN
*. Ensure that number of Perturbation subspaces have been defined
            IF(ISETKW(58).EQ.0) THEN
               WRITE(6,*) ' Dear User'
               WRITE(6,*)
               WRITE(6,*) ' H0SPC must be specified before H0FRM'
               WRITE(6,*)
     &              ' Else I do not know about the number of spaces'
               WRITE(6,*) ' So I will stop '
               Call Abend2( 'READIN: put H0SPC before H0FRM ' )
            END IF
*. Type of perturbation in this subspace
*
* 1 => Diagonal MP
* 2 => EN
* 3 => Nondiagonal MP
* 4 => Exact Hamiltonian
* 5 => Nondiagonal FI+FA + exact in orbital subspaces
*
            DO JPTSPC = 1, NPTSPC
               READ(LUIN,*) IH0INSPC(JPTSPC)
            END DO
            ISETKW(59) = 2
            GOTO 999
         END IF
*
* =============================================
* 60 : Reference root for Perturbation theory
* =============================================
*
         IF(CARD(2:7).EQ.'RFROOT') THEN
            ISETKW(60) = 1
            READ(LUIN,*) IRFROOT
C        WRITE(6,*) ' Reference Root = ',IRFROOT
            GOTO 999
         END IF
*
* ======================================================
* 61 : Orbital spaces in which Exact Hamiltonian is used
* ======================================================
*
         IF(CARD(2:5).EQ.'H0EX') THEN
            ISETKW(61) = 1
            READ(LUIN,*)  NH0EXSPC
            READ(LUIN,*) (IH0EXSPC(I),I=1, NH0EXSPC)
C?       WRITE(6,*) ' Keyword : H0EX activated '
C?       WRITE(6,*) '  NH0EXSPC ',  NH0EXSPC
C?       WRITE(6,*) (IH0EXSPC(I),I=1, NH0EXSPC)
            GOTO 999
         END IF
*
* ================================================
* 62 : Treatment of degenerencies of initial guess
* ================================================
*
         IF(CARD(2:7).EQ.'INIDEG') THEN
            ISETKW(62) = 1
            READ(LUIN,*) INIDEG
            GOTO 999
         END IF
*
* ========================================================
* 63 : Use modified Hamilton operator in CI optimization
* ========================================================
*
         IF(CARD(2:7).EQ.'LAMBDA') THEN
            ISETKW(63) = 1
            READ(LUIN,*) XLAMBDA
            GOTO 999
         END IF
*
* =============================================================
* 64 : Length of smallest block for batch of C an Sigma vectors
* =============================================================
*
         IF(CARD(2:7).EQ.'LCSBLK') THEN
            ISETKW(64) = 1
            READ(LUIN,*) LCSBLK
            GOTO 999
         END IF
*
*
* =============================================================
* 66 : No MO-AO file
* =============================================================
*
         IF(CARD(2:7).EQ.'NOMOFL') THEN
C?      write(6,*) ' NOMOFL flag active '
*. No MO-AO file
            NOMOFL = 1
            ISETKW(66) = 1
            GOTO 999
         END IF
*
*
* =============================================================
* 67 : ECHO the following keywords
* =============================================================
*
         IF(CARD(2:5).EQ.'ECHO') THEN
            IECHO = 1
            ISETKW(67) = 1
            GOTO 999
         END IF
*
*
* ====================
* 68 : Final orbitals
* ====================
*
*. Should be specified after NIRREP, I have not added the
* test!!
         IF(CARD(2:7).EQ.'FINORB') THEN
*. Type of final orbitals
            READ(LUIN,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
*. Change to upper case
            CALL UPPCAS(CARD1,MXPLNC)
*
            WRITE(6,'(A,A)')
     &           ' Type of final orbitals ',CARD1
            ISETKW(68) = 1
*
            IF(CARD1(1:5).EQ.'NATUR') THEN
*. Natural orbitals
               IFINMO = 1
               read(LUIN,*) NATITER
               write(6,*) 'Natural Orbitals. ',NATITER,'iterations'
            ELSE IF ( CARD(1:5).EQ.'CANON' ) THEN
*. Canonical orbitals
               IFINMO = 2
            ELSE IF ( CARD1(1:6).EQ.'PS_NAT') THEN
*. Pseudo natural orbitals
               IFINMO = 3
            ELSE IF ( CARD1(1:6) .EQ. 'PS_CAN') THEN
*. Pseudo canonical orbitals
               IFINMO = 4
            ELSE IF (CARD1(1:6) .EQ. 'PS_NC') THEN
*. Pseudo natural-canonical orbitals
               IFINMO = 5
*. requires input of subshells in which to define
*. Pseudo-natural orbitals
               READ(LUIN,*) NPSSPC
               DO IPSSPC = 1, NPSSPC
                  READ(LUIN,*) (NPSSH(IRREP,IPSSPC),IRREP=1,NIRREP)
               END DO
            else if (CARD1(1:6).eq.'NATSPI') then
*. Number of iterations in nat. orbital/spinor calculation
               read(LUIN,*) NATITER
               write(6,*) 'Natural Spinors. ',NATITER,'iterations'
               IFINMO = 6
            ELSE
*. Unidentified type of final orbitals
               WRITE(6,*) ' Unidentified type of final orbitals'
               WRITE(6,'(A,A)') '  you suggested : ', CARD1
               WRITE(6,*)
               WRITE(6,*) ' Allowed types of final orbitals'
               WRITE(6,*) ' ==============================='
               WRITE(6,*)
               WRITE(6,*) '     NATUR'
               WRITE(6,*) '     CANON'
               WRITE(6,*) '     PS_NAT'
               WRITE(6,*) '     PS_CAN'
               WRITE(6,*) '     PS_NC'
               WRITE(6,*) '     NATSPI'
               NERROR = NERROR + 1
               ISETKW(68) = - 1
            END IF
            GOTO 999
         END IF
*
*
* ===================================================================
* 69 : Threshold on second order energy corrections, individual coefs
* ===================================================================
*
         IF(CARD(2:7).EQ.'E_THRE') THEN
            READ(LUIN,*) E_THRE
            ISETKW(69) = 1
            GOTO 999
         END IF
*
*
* =======================================================================
* 70 : Threshold on first order wavefunction corrections,individual coefs
* =======================================================================
*
         IF(CARD(2:7).EQ.'C_THRE') THEN
            READ(LUIN,*) C_THRE
            ISETKW(70) = 1
            GOTO 999
         END IF
*
* ===================================================================
* 71 : Threshold on second order energy corrections, Total Threshold
* ===================================================================
*
         IF(CARD(2:7).EQ.'E_CONV') THEN
            READ(LUIN,*) E_CONV
            ISETKW(71) = 1
            GOTO 999
         END IF
*
*
* =======================================================================
* 72 : Threshold on first order wavefunction corrections,Total Threshold
* =======================================================================
*
         IF(CARD(2:7).EQ.'C_CONV') THEN
            READ(LUIN,*) C_CONV
            ISETKW(72) = 1
            GOTO 999
         END IF
*
*
* ===============================
* 73 : Selection of classes
* ===============================
*
         IF(CARD(2:7).EQ.'CLSSEL') THEN
            ICLSSEL = 1
            ISETKW(73) = 1
            GOTO 999
         END IF
*
*
* =====================================
* 74 : Calculation of density matrices
* ======================================
*
         IF(CARD(2:6).EQ.'DENSI') THEN
            READ(LUIN,*) IDENSI
*. IDENSI = 0 =>  No calculation of density matrices
*  IDENSI = 1 =>  Calculation of one- body density matrix
*  IDENSI = 2 =>  Calculation of one- and two-body density matrices
*  IDENSI = 3 =>  Calculation of spin-dependent one-body density matrices
*  IDENSI = 4 =>  Calculation of spin-dependent one- and two-body
*                                density matrices
            ISETKW(74) = 1
            GOTO 999
         END IF
*
*
* =====================================
* 75 : Perturbation expansion of EKT
* ======================================
*
         IF(CARD(2:6).EQ.'PTEKT') THEN
            IPTEKT = 1
*. Number of EKT to be analyzed, atmost 20
            READ(5,*)  NPTEKT
            IF(NPTEKT.GT.20) THEN
               WRITE(6,*) ' Atmost 20 perturbation expansions'
               Call Abend2( ' NPTEKT in .PTEKT to Large ' )
            END IF
*. orbital and symmetry for zero order solution
            DO JEKT = 1, NPTEKT
               READ(5,*) LPTEKT(1,JEKT),LPTEKT(2,JEKT)
            END DO
            ISETKW(75) = 1
            WRITE(6,*) ' NPTEKT = ', NPTEKT
            WRITE(6,*) ' LPTEKT = ',LPTEKT(1,1),LPTEKT(2,1)
            GOTO 999
         END IF
*
* =================================================
* 76 : Root used to define Zero order Hamiltonian
* =================================================
*
         IF(CARD(2:7).EQ.'H0ROOT') THEN
            ISETKW(76) = 1
            READ(LUIN,*) IH0ROOT
C        WRITE(6,*) ' Reference Root = ',IH0ROOT
            GOTO 999
         END IF
*
* ======================================
* 77 : No restart in CI calculation 2
* =====================================
*
         IF(CARD(2:7).EQ.'NORST2') THEN
            ISETKW(77) = 1
            IRST2 =  0
            WRITE(6,*) ' NORST2 flag read '
            GOTO 999
         END IF
*
* =====================================================
* 78 : Skip initial evaluation of energy from CI calc 2
* ====================================================
*
         IF(CARD(2:7).EQ.'SKIPEI') THEN
            ISETKW(78) = 1
            ISKIPEI =  1
            WRITE(6,*) ' SKIPEI flag set  '
            GOTO 999
         END IF
*
* =================================================================
* 79 : Symmetry of X, Y and Z - Yes it could be obtained from files
* ================================================================
*
         IF(CARD(2:7).EQ.'XYZSYM') THEN
            ISETKW(79) = 1
            READ(LUIN,*) (IXYZSYM(I),I=1,3)
            WRITE(6,*) 'IXYZSYM', (IXYZSYM(I),I=1,3)
            GOTO 999
         END IF
*
* ==============================================
* 80 : One-electron properties to be calculated
* ==============================================
*
         IF(CARD(2:7).EQ.'PROPER') THEN
            ISETKW(80) = 1
            READ(LUIN,*) NPROP
            DO IPROP = 1, NPROP
               READ(LUIN,'(A)') CARD1
               CALL LFTPOS(CARD1,MXPLNC)
               CALL UPPCAS(CARD1,MXPLNC)
               PROPER(IPROP)=CARD1(1:6)
               WRITE(6,'(A,A)') ' Property to be calculated ',
     &              PROPER(IPROP)
            END DO
            GOTO 999
         END IF
*
* =================================
* 81 : CI in double group symmetry
* =================================
*
         IF(CARD(2:6).EQ.'DBGCI') THEN
            ISETKW(81) = 1
            IDBGCI = 1
            WRITE(6,*) ' Double group CI flag set '
            GOTO 999
         END IF
*
* ======================
* 82 : Relativistic CI
* ======================
*
         IF(CARD(2:7).EQ.'RELACI') THEN
            ISETKW(82) = 1
            IRELACI = 1
            WRITE(6,*) ' RELACI  flag set '
            GOTO 999
         END IF
*
* ================================================
* 83 : Type of spin-orbit operator to be included
* ================================================
*
         IF(CARD(2:6).EQ.'SO-OP') THEN
            READ(LUIN,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
            CALL UPPCAS(CARD1,MXPLNC)
            IF(CARD1(1:3).EQ.'1+2'.OR.CARD1(1:7).EQ.'ONE+TWO'.OR.
     &         CARD1(1:3).EQ.'2+1'.OR.CARD1(1:7).EQ.'TWO+ONE')THEN
               ISO1 = 1
               ISO2 = 1
               ISETKW(83) = 1
               GOTO 999
            ELSE IF(CARD1(1:1).EQ.'1'.OR.CARD1(1:3).EQ.'ONE') THEN
               ISO1 = 1
               ISO2 = 0
               ISETKW(83) = 1
               GOTO 999
            ELSE IF(CARD1(1:1).EQ.'2'.OR.CARD1(1:3).EQ.'TWO') THEN
               ISO1 = 0
               ISO2 = 1
               ISETKW(83) = 1
               GOTO 999
            else if(CARD1(1:1).EQ.'0'.or.CARD(1:6).eq.'SPINFR') then
CTF
* Added for test runs without SO integrals. And DIRAC spinfree
*  formalism.
               ISPINFREE = 1
               ISO1 = 0
               ISO2 = 0
               ISETKW(83) = 1
               goto 999
CTF
            ELSE
               WRITE(6,*) ' Unidentified form of SO-operator'
               WRITE(6,'(A,A)') '  you suggested : ', CARD1
               WRITE(6,*) ' Allowed forms : '
               WRITE(6,*)
     &         ' One+two, 1+2  => both one- and two-electron operators'
               WRITE(6,*)
     &         ' One, 1        =>      one-electron operators'
               WRITE(6,*)
     &         ' Two, 2        =>      two-electron operators'
               NERROR = NERROR + 1
               ISO1 = 0
               ISO2 = 0
               ISETKW(83) = -1
            END IF
         END IF
*
* =================================
* 84 : Reference MK2 value
* =================================
*
         IF(CARD(2:7).EQ.'MK2REF') THEN
            ISETKW(84) = 1
            READ(LUIN,*) MK2REF_CI
            GOTO 999
         END IF
*
* =================================
* 85 : Delta MK2 value
* =================================
*
         IF(CARD(2:7).EQ.'MK2DEL') THEN
            ISETKW(85) = 1
            READ(LUIN,*) MK2DEL_CI
            GOTO 999
         END IF
*
* ============================================
* 86 : Initialize with integral transformation
* ============================================
*
         IF(CARD(2:6).EQ.'TRANS') THEN
            ISETKW(86) = 1
            READ(LUIN,*) ITRANS
            GOTO 999
         END IF
*
* =================================================
* 87 : Spinor expansion from MO expansion
* =================================================
*
         IF(CARD(2:7).EQ.'SP<-MO') THEN
            ISETKW(87) = 1
            ISPEXP_INI = 1
            GOTO 999
         END IF
*
* =================================================
* 88 : Spinor expansion read in
* =================================================
*
         IF(CARD(2:7).EQ.'SP<-RE') THEN
            ISETKW(88) = 1
            ISPEXP_INI = 8
            GOTO 999
         END IF
*
* ================================================
* 89 : Print flag for transformation information
* ================================================
*
         if(CARD(2:7).eq.'IPRTRA') then
            ISETKW(89) = 1
            read(LUIN,*) IPRTRA
            goto 999
         end if
*
* ================================================
* 90 : Print flag for sigma vector information
* ================================================
*
         if(CARD(2:7).eq.'IPRSIG') then
            ISETKW(90) = 1
            read(LUIN,*) IPRSIG
            goto 999
         end if
*
* ========================================================
* 91 : COMHAM: Test run diagonalizing complete Hamiltonian
* ========================================================
*
         if(CARD(2:7).eq.'COMHAM') then
            ISETKW(91) = 1
            IDCOMH = 1
            goto 999
         end if
*
* ========================================================
* 92 : ONEPAR: Do full Hamiltonian or only one-particle op.
* ========================================================
*
         if(CARD(2:7).eq.'PARHAM') then
            ISETKW(92) = 1
            read(LUIN,*) IPARHAM
            goto 999
         end if
*
* ==================================================================
* 93 : FOCKEX: Do CI with externally generated Fock matrix integrals
* ==================================================================
*
         if(CARD(2:7).eq.'FOCKEX') then
            ISETKW(93) = 1
            IFOCKEX = 1
            read(LUIN,*) (NORB_F(ISYM),ISYM=1,NSMSH,1)
            write(6,*)
            goto 999
         end if
*
* ======================================================================
* 94 : INPORB: Read MO vector from formatted ASCII file (MOLCAS: INPORB)
* ======================================================================
*
         if(CARD(2:7).eq.'INPORB') then
            ISETKW(94) = 1
            ORBFILE = .true.
            goto 999
         end if
*
* ======================================================================
* 95 : IPRWFC: Print more or less detailed info about CI expansions 
*              (opt.)
* ======================================================================
*
         if(CARD(2:7).eq.'IPRWFC') then
            ISETKW(95) = 1
            read(LUIN,*) IPRWFC
            goto 999
         end if
*
* ======================================================================
* 97 : SIGDEN: New route for sigma vectors and densities 
*              (one set of rout.)
* ======================================================================
*
         if(CARD(2:7).eq.'SIGDEN') then
           ISETKW(97) = 1
           goto 999
         end if
*
            read(LUIN,*) IPRHAM_CI
* ======================================================================
* 98 : IPRHAM: Print flag for Hamiltonian definition
* ======================================================================
*
         if (CARD(2:7).eq.'IPRHAM') then
           ISETKW(98) = 1
           goto 999
         end if
*
*. KEYWORD was not identified
*
         WRITE(6,*)
     &       '  ****  Error, unidentified KEYWORD in READIN   **** '
         WRITE(6,*)
         WRITE(6,*) ' Last line read  '
         WRITE(6,*) ' ================'
         WRITE(LUOUT,'(10X,A)') CARD
         WRITE(6,*)
         WRITE(6,*) ' Preceeding KEYWORD'
         WRITE(6,*) ' ==================='
         WRITE(LUOUT,'(10X,A)') LASTCARD
         NERROR = NERROR + 1
*
      END IF
 999  CONTINUE
      IF(IECHO.EQ.1)
     &     WRITE(6,'(A,A)') ' processed KEYWORD/COMMENT  : ', CARD
*. Save previous keyword
      LASTCARD(1:72) = CARD(1:72)
      GOTO 1000
*.End of loop over KEYWORDS
 1001 CONTINUE
*
      IF(NERROR.NE.0) THEN
         WRITE(LUOUT,'(A)')
     &        ' Run will be aborted due to input errors '
         WRITE(LUOUT,'(A,I9)')
     &        ' Number of input errors detected in READIN ', NERROR
*
         WRITE(6,*) ' The following keywords were correctly identified'
         WRITE(6,*) ' ================================================'
         DO  IENTRY = 1, MXPKW
            IF(ISETKW(IENTRY).EQ.+1)
     &           WRITE(LUOUT,'(10X,A)') KEYWOR(IENTRY)
         END DO
         WRITE(6,*)
         WRITE(LUOUT,*)
         WRITE(LUOUT,*)
         WRITE(LUOUT,*)
         WRITE(LUOUT,*)
     &        '     An expert is a man who has made all the mistakes,'
         WRITE(LUOUT,*)
     &        '     which can be made, in a very narrow field        '
         WRITE(LUOUT,*)
     &        '                                                      '
         WRITE(LUOUT,*)
     &        '                                      Niels Bohr      '

         Call Abend1( 5 )
      END IF
*
**********************************************************************
*                                                                    *
* Part 2 : Insert defaults for missing optional keywords             *
*          and print error messages for missing mandatory keywords   *
*                                                                    *
**********************************************************************
*
      NMISS = 0
*
*.1 : Default title
*
      IF(ISETKW(1).EQ.0) THEN
         TITLEC(1) =
     &        ' Some molecule or some atom                           '
         TITLEC(2) =
     &        ' Some type of CI expansion                            '
         TITLEC(3) =
     &        ' Some user who is too lazy to supply a TITLE          '
         ISETKW(1) = 2
      END IF
*
*.2  Missing pointgroup ( has actually been defaulted )
*
      IF(ISETKW(2).EQ.0) THEN
         PNTGRP = 1
         ISETKW(2) = 2
      END IF
*
*.3 Missing number of irreps, allowed for D2h and DBG, illegal else
*
      IF(ISETKW(3).EQ.1.AND.PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
*
*. Double group symmetry in action
*
*. NIRREP was defined, but it is actually known from
*. Double group label. Use the latter value !
*. ( we computers do not really trust you human beings,
*.   You talk so much about the errors we do, but that
*.   is nothing compared to the amount of errors you do
*    - and we are polite enough not to complain)
         IF(PNTGRP.EQ.5) THEN
            NIRREP = 2
         ELSE IF(PNTGRP.EQ.6.OR.PNTGRP.EQ.7.OR.PNTGRP.EQ.8) THEN
            NIRREP = 4
         ELSE IF(PNTGRP.EQ.9) THEN
            NIRREP = 8
         END IF
         ISETKW(3) = 2
      END IF
*
      IF(ISETKW(3).EQ.0) THEN
         IF(PNTGRP .EQ.1 ) THEN
*. Repeat default
            NIRREP = 8
            NSMCMP = NIRREP
            NSMOB  = NIRREP
            ISETKW(3) = 2
         ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
*. Double group
            IF(PNTGRP.EQ.5) THEN
               NIRREP = 2
            ELSE IF(PNTGRP.EQ.6.OR.PNTGRP.EQ.7.OR.PNTGRP.EQ.8) THEN
               NIRREP = 4
            ELSE IF(PNTGRP.EQ.9) THEN
               NIRREP = 8
            END IF
            ISETKW(3) = 2
         ELSE
*. Number of irreps is mandatory for CINV,DINFH,O3
            NMISS = NMISS + 1
            WRITE(LUOUT,*)
     &     '  Input error ! .NIRREP must be specified for CinV,DinH,O3'
         END IF
      END IF
*
* 4 : Internal CI expansion
*UNUSED
*.Default is CAS
      IF(ISETKW(50) .EQ. 0 ) THEN
         IF(ISETKW(4).EQ.0) THEN
            INTSPC = 1
            ISETKW(4) = 2
*. If a RAS1 or a RAS 3 space has been defined, RAS must have
*  been specified
            IF(ISETKW(9).EQ.1.OR.ISETKW(11).EQ.1) THEN
               ISETKW(4) = 0
               NMISS = NMISS + 1
               WRITE(LUOUT,*)
     &     '  Input error ! .RAS must be specified when .RAS1 or .RAS3'
               WRITE(LUOUT,*)
     &              '                 has been activated '

            END IF
         END IF
      ELSE IF (ISETKW(50) .EQ. 0 ) THEN
*. Gas expansion
         INTSPC = 3
      END IF
*
* 6 : Number of active electrons
*
*. Mandatory
      if (ISETKW(6).eq.0) then
        NACTEL = IGSOCCX(NGAS,2,1)
        if (NACTEL.eq.1) then
          IELCI = 1
        else
          IELCI = 2
        end if
        ISETKW(6) = 2
      end if
*
* 7 : Inactive orbitals
*UNUSED
      IF(ISETKW(7).EQ.0) THEN
         CALL ISETVC(NINASH,0,NIRREP)
         ISETKW(7) = 0
      END IF
*
* 8 : Core orbitals, only of interest if EXTSPC .ne. 0
*UNUSED
      IF(ISETKW(8).EQ.0) THEN
         CALL ISETVC(NRS0SH,0,NIRREP)
         MNHR0 = 0
         IF(EXTSPC.EQ.0) THEN
            ISETKW(8) = 3
         ELSE
            ISETKW(8) = 2
         END IF
      END IF
*
* 9 : RAS 1 orbitals
*UNUSED
      IF(ISETKW(9).EQ.0) THEN
         CALL ISETVC(NRSSH(1,1),0,NIRREP)
         IF(INTSPC.EQ.1) THEN
            ISETKW(9) = 3
         ELSE
            ISETKW(9) = 2
         END IF
      END IF
*
* 10 : RAS 2 orbitals
*UNUSED
      IF(ISETKW(10).EQ.0) THEN
         CALL ISETVC(NRSSH(1,2),0,NIRREP)
         ISETKW(10) = 2
      END IF
*
* 11 : RAS 3 orbitals
*UNUSED
      IMLCR3 = 0
      IF(ISETKW(11).EQ.0) THEN
         CALL ISETVC(NRSSH(1,3),0,NIRREP)
         IF(MOLCS.EQ.1.AND.INTSPC.EQ.2) THEN
*. Use information from one-electron integral file to obtain
* default
            IMLCR3 = 1
         END IF
         IF(INTSPC.EQ.1) THEN
            ISETKW(11) = 3
         ELSE
            ISETKW(11) = 2
         END IF
      END IF
*
* 12 : Partitioning of secondary space ( default 1 set in SECOND)
*
C     IF(ISETKW(12).EQ.0.OR.ISETKW(12).EQ.2) THEN
C       MXR4TP = 1
C       IF(EXTSPC.EQ.0) THEN
C         ISETKW(12) = 3
C       ELSE
C         ISETKW(12) = 2
C       END IF
C     END IF
*
* 13 : Secondary space
*UNUSED
      IF(ISETKW(13).EQ.0) THEN
         MXR4TP = 1
         DO 210 IR4TP = 1, MXR4TP
            CALL ISETVC(NRS4SH(1,IR4TP),0,NIRREP)
 210     CONTINUE
         MXER4 = 0
         IF(EXTSPC.EQ.0) THEN
            ISETKW(13) = 3
         ELSE
            ISETKW(13) = 2
         END IF
      END IF
*
* 14 : occupation restrictions for Reference space
*UNUSED
      IF(ISETKW(14).EQ.0) THEN
         MNRS1R = MNRS10
         MXRS3R = MXRS30
         IF(EXTSPC.EQ.0.OR.INTSPC.EQ.1) THEN
            ISETKW(14) = 3
         ELSE
            ISETKW(14) = 2
         END IF
      END IF
*
* 15 : Selection of active configurations
*UNUSED
      IF(ISETKW(15).EQ.0) THEN
*. Standard is no selection
         INTSEL = 0
      END IF
*
* 16 : Two times spin projection
*UNUSED (only if not double group CI and not relativistic CI)
      IF(ISETKW(16).EQ.0) THEN
        MS2 = MK2REF_CI
        ISETKW(16) = 2
      END IF
*
* 17 : Spin multiplicity
*UNUSED
C     IF(ISETKW(17).EQ.0) THEN
C       WRITE(LUOUT,*)
C    &  '  Input error ! .MULTS must be specified '
C       NMISS = NMISS + 1
C     END IF
*
* 18 : Reference symmetry
*
      IF(ISETKW(18).EQ.0) THEN
        WRITE(LUOUT,*)
     &       '  Input error ! .IREFSM must be specified '
        NMISS = NMISS + 1
      else
        if (IREFSM.gt.NIRREP) then
          write(6,*) 'Fatal error in input.'
          write(6,*) 'IREFSM larger than NIRREP.'
          write(6,*) 'IREFSM = ',IREFSM,' ,NIRREP = ',NIRREP
          Call Abend2( 'Please reconsider.' )
        end if
      END IF
*
* 19 : Roots to be optimized
*
      IF(ISETKW(19).EQ.0) THEN
        WRITE(LUOUT,*)
     &       '  Input error ! .ROOTS must be specified '
        NMISS = NMISS + 1
      else
        if (NROOT.gt.MXROOT) then
          write(6,*) 'Fatal error in input.'
          write(6,*) 'Number of roots to be obtained is',NROOT
          write(6,*) 'But max. no. of roots is ',MXROOT
          Call Abend2( 'Try again.' )
        end if
      END IF
*
* 20 : Diagonalization routine
*
      IF(ISETKW(20).EQ.0) THEN
*. Standard is currently MINDV4
         IDIAG = 1
         ISETKW(20) = 2
      END IF
*
* 21 : Explicit Hamiltonian
*
      IF(ISETKW(21).EQ.0) THEN
*. Default is no explicit Hamiltonian
         MXP1 = 0
         MXP2 = 0
         MXQ  = 0
         ISETKW(21) = 2
      END IF
*
* 22 : Largest allowed number of CI iterations per root
*
      IF(ISETKW(22).EQ.0) THEN
*. Default is no explicit Hamiltonian
         MAXIT = 40
         ISETKW(22) = 2
      END IF
*
* 23 : Restart option
*
      IF(ISETKW(23).EQ.0) THEN
*. Default is no explicit Hamiltonian
         IRESTR = 0
         ISETKW(23) = 3
      END IF
*
* 24 : Integral import
*
      IF(ISETKW(24).EQ.0) THEN
*. Default is - from NOV26 : Dalton
         INTIMP = 5
         ENVIRO(1:6) = 'DALTON'
         ISETKW(24) = 2
      END IF
*
* 25 : INCORE option for integrals
*
      IF(ISETKW(25).EQ.0) THEN
         IF(EXTSPC.EQ.0 ) THEN
            INCORE = 1
         ELSE
            INCORE = 0
         END IF
         ISETKW(25) = 2
C
C       IF(INTEXP.EQ.0) THEN
C         ISETKW(25) = 3
C       ELSE
C         ISETKW(25) = 3
C       END IF

      END IF
*
* 26 : DELETEd shells
*UNUSED
      IF(ISETKW(26) .EQ. 0 ) THEN
*. If CAS + Active have been set or RAS + Ras3 have been set,
*. obtain for MOLCAS Interface from number of basis functions
         IF(INTSPC.EQ.1.OR.
     &        (INTSPC.EQ.2.AND.ISETKW(11).EQ.1)) THEN
            IMLCR3 = 2
         ELSE
            CALL ISETVC(NDELSH,0,NIRREP)
         END IF
         ISETKW(26) = 2
      END IF
*
* 27 : Ms combinations
*
      IF(ISETKW(27).EQ.0) THEN
         PSSIGN = 0.0D0
         ISETKW(27) = 2
      ELSE IF(MS2.NE.0) THEN
       WRITE(LUOUT,*) ' Spin combinations only allowed with MS2 = 0'
       WRITE(LUOUT,*) ' Your value of MS2 = ',MS2, ' differs from zero'
       WRITE(LUOUT,*) ' LUCIA will neglect your nice suggestion '
       WRITE(LUOUT,*)  ' to use spin combinations '
       PSSIGN = 0.0D0
       ISETKW(27) = 2
      END IF
*
* 28 : Ml combinations
*
      IF(ISETKW(28).EQ.0) THEN
         PLSIGN = 0.0D0
         ISETKW(28) = 2
      ELSE IF(PNTGRP.EQ.1) THEN
         WRITE(LUOUT,*) ' Ml combinations not allowed with d2h '
         WRITE(LUOUT,*) ' LUCIA will neglect your nice suggestion '
         WRITE(LUOUT,*)  ' to use ML combinations '
         PLSIGN = 0.0D0
         ISETKW(28) = 2
      ELSE IF(IREFML.NE.0) THEN
         WRITE(LUOUT,*) ' ML combinations only allowed with ML = 0'
         WRITE(LUOUT,*)
     &        ' Your value of IREFML = ',IREFML, ' differs from zero'
         WRITE(LUOUT,*) ' LUCIA will neglect your nice suggestion '
         WRITE(LUOUT,*)  ' to use ML combinations '
         PLSIGN = 0.0D0
         ISETKW(28) = 2
      END IF
      IF(PSSIGN.EQ.0.0D0.AND.PLSIGN.EQ.0.0D0) THEN
         IDC = 1
      ELSE IF(PSSIGN.NE.0.0D0.AND.PLSIGN.EQ.0.0D0) THEN
         IDC = 2
      ELSE IF(PSSIGN.EQ.0.0D0.AND.PLSIGN.NE.0.0D0) THEN
         IDC = 3
      ELSE IF(PSSIGN.NE.0.0D0.AND.PLSIGN.NE.0.0D0) THEN
         IDC = 4
      END IF
C?    WRITE(6,* ) ' TEST readin IDC = ', IDC
*
* 29 : print flag for string information
*
      IF(ISETKW(29).EQ.0) THEN
         IPRSTR = 0
         ISETKW(29) = 2
      END IF
*
* 30 : print flag for string information
*
      IF(ISETKW(30).EQ.0) THEN
         IPRCIX = 0
         ISETKW(30) = 2
      END IF
*
* 31 : print flag for orbital information
*
      IF(ISETKW(31).EQ.0) THEN
         IPRORB = 0
         ISETKW(31) = 2
      END IF
*
* 32 : print flag for diagonalization information
*
      IF(ISETKW(32).EQ.0) THEN
         IPRDIA = 3
         ISETKW(32) = 2
      END IF
*
* 36 : print flag for External blocks
*
      IF(ISETKW(36).EQ.0) THEN
         IPRXT  = 0
         ISETKW(36) = 2
      END IF
*
* 43 : Print occupation of lowest SD's / configurations
*
      IF(ISETKW(43).EQ.0) THEN
         IPROCC = 0
         ISETKW(43) = 2
      END IF
*
* 65 : Print level for densities, default is to print
*      natural occupation numbers only
*
      IF(ISETKW(65).EQ.0) THEN
         IPRDEN = 1
         ISETKW(65) = 2
      END IF
*
* 89 : print flag for transformation information
*
      IF(ISETKW(89).EQ.0) THEN
         IPRTRA = 0
         ISETKW(89) = 2
      END IF
*
* 90 : print flag for sigma vector information
*
      IF(ISETKW(90).EQ.0) THEN
         IPRSIG = 0
         ISETKW(90) = 2
      END IF
*
* 91 : Test run with complete Hamiltonian
*
      if (ISETKW(91).eq.0) then
         IDCOMH = 0
         ISETKW(91) = 2
      end if
*
* 33 : Number of CI vectors in subspace
*
      IF(ISETKW(33).EQ.0) THEN
* default is 3 vectors per root
         MXCIV_CI = 3 * NROOT
         ISETKW(33) = 2
      END IF
      IF(ISETKW(33).EQ.1.AND.MXCIV_CI .LT.2*NROOT) THEN
         WRITE(LUOUT,*)
     &    '   The number of vectors is increased to 2*NROOT = ',2*NROOT
         MXCIV_CI = 2 * NROOT
      END IF
*
* 34 : CI storage mode
*
      IF(ISETKW(34).EQ.0) THEN
        IF(IDIAG.EQ.1) THEN
*: Default is three symmetry blocks and two vectors
c         ICISTR = 1
          ISETKW(34) = 2
        ELSE
*.Default is three symmetry blocks
c         ICISTR = 2
          ISETKW(34) = 2
        END IF
      END IF
      IF(IDIAG.EQ.2.AND.ICISTR.EQ.1) THEN
*. complete vectors should not be used together with PICO
      WRITE(LUOUT,*)
     &'    You have suggested the use of two complete vectors in core'
      WRITE(LUOUT,*)
     &'    together with the TERACI option '
      WRITE(LUOUT,*)
     &'    Although this could be an interesting suggestion '
      WRITE(LUOUT,*)
     &'    I allow myself to reduce the storage mode to 3 sym. blocks '
      ICISTR = 2
      END IF
*
* 35 : Employ CSF expansion ?
*
      IF(ISETKW(35).EQ.0) THEN
        NOCSF = 0
        ISETKW(35) = 2
      END IF
* CSF expansion must only be used when two vectors are stored in CORE
      IF(NOCSF.EQ.0.AND.ICISTR.NE.1) THEN
        WRITE(LUOUT,*)
     &  '   Do not use CSF expansions for blocked CI calculations ! '
        WRITE(LUOUT,*)
     &  '   Employ the NOCSF option or set CISTOR to 1 ! '
        WRITE(LUOUT,*)
     &  '   I will stop and wait for your decision  ! '
        Call Abend1( 11 )
      END IF
*
* 37 : Avoid any readin of integrals ( useful for obtaining
*      size of CI expansion etc.
*
      IF(ISETKW(37).EQ.0 ) THEN
        NOINT = 0
        ISETKW(37) = 2
      END IF
*
* 38 : Dump integrals in formatted form : Default is No
*UNUSED
      IF(ISETKW(38).EQ.0) THEN
        IDMPIN = 0
        ISETKW(38) = 2
      END IF
*
* 39 : Explicitly dimension of dimension of block of resolution strings
*
      IF(ISETKW(39).EQ.0) THEN
        MXINKA = 500
        ISETKW(39) = 2
      END IF
*
* 40 : Use CJKAKB intermediate matrices in alpha-beta loop,
*      Default is  YES !!!!!
*UNUSED
      IF(ISETKW(40).EQ.0) THEN
        ICJKAIB = 1
        ISETKW(40) = 2
      END IF
*
*  41 : Initial CI in reference space, default is : No
*UNUSED
      IF(ISETKW(41).EQ.0) THEN
         ISETKW(41) = 2
         INIREF = 0
      END IF
*
*  42 : Restart with CI in reference space
*UNUSED
      IF(ISETKW(42).EQ.0) THEN
         ISETKW(42) = 2
         IRESTRF = 0
      END IF
*
*  44 : Use MOC method for alpha-alpha loop, default is NO !
*
      IF(ISETKW(44).EQ.0) THEN
         ISETKW(44) = 2
         MOCAA = 0
      END IF
*
*  45 : Use MOC method for alpha-beta loop, default is NO !
*
      IF(ISETKW(45).EQ.0) THEN
         ISETKW(45) = 2
         MOCAB = 0
      END IF
*
* 46 : Core energy : Default is 0 / MOLCAS : Value read in !
*
      IF(ISETKW(46).EQ.0) THEN
         ISETKW(46) = 2
         ECORE = 0.0D0
      END IF
*
* 47 : Use perturbation theory for zero order space . def is no !
*UNUSED
      IF(ISETKW(47).EQ.0) THEN
        IPERT = 0
        NPERT = 0
        ISETKW(47) = 2
*. Else ensure that a CI in reference space is performed
      ELSE
        INIREF = 1
      END IF
*
*
*. 48 : Approximate Hamiltonian in reference space : NO !!
*UNUSED
      IF(ISETKW(48).EQ.0) THEN
        IAPRREF = 0
        MNRS1RE = MNRS1R
        MXRS3RE = MXRS3R
        ISETKW(48) = 2
      END IF
*
*. 49 : Approximate Hamiltonian in zero order space : NO !!
*UNUSED
      IF(ISETKW(49).EQ.0) THEN
        IAPRZER = 0
        MNRS1ZE = MNRS10
        MXRS3ZE = MXRS30
        ISETKW(49) = 2
      END IF
*
* 50 : GAS expansion, default is : NO !!!
*UNUSED
      IF(ISETKW(50).EQ.0) THEN
        IDOGAS = 0
        ISETKW(50) = 2
      END IF
*
* 52 : Combination of gasspaces : Default is just to take each  space
*      By itself
*
      IF(ISETKW(52).EQ.0) THEN
        NCMBSPC = NCISPC
        DO ICISPC = 1, NCISPC
          LCMBSPC(ICISPC) = 1
          ICMBSPC(1,ICISPC) = ICISPC
        END DO
        ISETKW(52) = 2
      END IF
*
* 53 : Convergence threshold for CI
*
      IF(ISETKW(53).EQ.0) THEN
        THRES_E = 1.0D-10
        ISETKW(53) = 2
      END IF
*
* 54 : General sequencer : default is just straight sequence
*      of CI with default number of iterations
      IF(ISETKW(54).EQ.0) THEN
        DO JCMBSPC = 1, NCMBSPC
          NSEQCI(JCMBSPC) = 1
          ISEQCI(1,JCMBSPC) = MAXIT
        END DO
        ISETKW(54) = 2
      END IF
*
* 55 : EKT calculation : Default is no
*UNUSED
      IF(ISETKW(55).EQ.0) THEN
        IEXTKOP = 0
        ISETKW(55) = 2
      END IF
*
*. 56 : Default Machine : Good old IBM machine
*UNUSED
      IF(ISETKW(56).EQ.0) THEN
c       MACHINE(1:3) = 'IBM'
c       ISETKW(56) = 2
      END IF
*
* 57 : Allow first order correction to be saved on DISC
*UNUSED     (For vector free calculations )
*     Default is : NO !!
      IF(ISETKW(57).EQ.0) THEN
        IC1DSC = 0
        ISETKW(57) = 2
      END IF
*
* 58 : Restrictions on interactions of perturbation
*UNUSED
*. Default is : no
      IF(ISETKW(58).EQ.0) THEN
        IH0SPC = 0
        ISETKW(58) = 2
      END IF
*
* 59 : Type of perturbation in subspaces spaces
*UNUSED
* Default is specified by IPART from keyword PERTU
      IF(ISETKW(59).EQ.0) THEN
       ISETKW(59) = 2
       IF(IH0SPC.NE.0) THEN
         DO JPTSPC = 1, NPTSPC
           IH0INSPC(JPTSPC) = IPART
         END DO
       END IF
      END IF
*
* 60 : Reference Root, default is NROOT
*
*. Should be less or equal to NROOT
      IF(ISETKW(60).EQ.1) THEN
        IF(IRFROOT.GT.NROOT) THEN
          WRITE(6,*) ' Reference root (RFROOT) larger '
          WRITE(6,*) ' than total number of roots (NROOT) '
          WRITE(6,*) ' CHANGE NROOT or RFROOT '
          NMISS = NMISS + 1
        END IF
      END IF

      IF(ISETKW(60).EQ.0) THEN
       ISETKW(60) = 2
       IRFROOT = NROOT
      END IF
*
* 61 : H0EX : Orbital spaces in which exaxt Hamiltonian is used
*      No default
*.
*. Is H0EX required ( Has H0FRM = 5 been used )
      IUSED = 0
      IF(ISETKW(59).EQ.1) THEN
         IUSED = 0
         DO JPTSPC = 1, NPTSPC
           IF( IH0INSPC(JPTSPC) .EQ. 5 ) IUSED = 1
         END DO
       END IF
       IF(IUSED.EQ.0.AND.ISETKW(61).EQ.0) THEN
*. No exact spaces included and none have been defined !
         NH0EXSPC = 0
         IH0EXSPC(1) = -1
       END IF
       IF(IUSED.EQ.1.AND.ISETKW(61).EQ.0) THEN
*. Needed, but not supplied
          WRITE(6,*) ' You have specified that zero order operator'
          WRITE(6,*) ' Include exact Hamilton operator in subspace'
          WRITE(6,*) ' You should then also supply Keyword H0EX '
          NMISS = NMISS + 1
       END IF
*
*. If perturbation theory will be invoked be sure that the
*. form of perturbation theory has been specified through
* KEYWORD PERTU ( number 47 as you maybe know )
      IDOPERT = 0
      DO JCMBSPC = 1, NCMBSPC
        DO JSEQCI = 1, NSEQCI(JCMBSPC)
          IF(ISEQCI(JSEQCI,JCMBSPC).EQ.-5 ) IDOPERT = 1
        END DO
      END DO
*
      IF(IDOPERT.EQ.1 .AND. IPERT.EQ.0) THEN
        WRITE(6,*) ' Perturbation theory will be used '
        WRITE(6,*) ' Please specify form through PERTU keyword '
        NMISS = NMISS + 1
      END IF
*
*. 62 : Default Handling of degenrences of initial CI vectors
*.      Default is : No action
*
      IF(ISETKW(62).EQ.0) THEN
        INIDEG = 0
        ISETKW(62) = 2
      END IF
*
*. 63 : Use F + Lambda(H-F) as operator instead of H
*.      Default is : No i.e Lambda = 1
*
      IF(ISETKW(63).EQ.0) THEN
        XLAMBDA = 1.0D0
        ISETKW(63) = 2
      END IF
*
*. 64 : Smallest block in batch of C and sigma
*.      Default is zero
*
      IF(ISETKW(64).EQ.0) THEN
        LCSBLK = 100000
        ISETKW(64) = 2
      END IF
*
*. 66 : NO MO file : Default is access to MO-AO file
*
      IF(ISETKW(66).EQ.0) THEN
        NOMOFL = 0
        ISETKW(66) = 2
      END IF
*
*. 68 : Type of natural orbitals, default is natural orbitals
*
      IF(ISETKW(68).EQ.0) THEN
        ISETKW(68) = 2
        IFINMO = 1
        NATITER = 1
      END IF
*
*. 69 : Default Threshold for individual energy correction = 0.0
*UNUSED
      IF(ISETKW(69).EQ.0) THEN
        E_THRE = 0.0D0
        ISETKW(69) = 2
      END IF
*
*. 70 : Default Threshold for wave individual function corrections = 0.0
*UNUSED
      IF(ISETKW(70).EQ.0) THEN
        C_THRE = 0.0D0
        ISETKW(70) = 2
      END IF
*
*. 71 : Default Threshold for total energy corrections = 0.0
*UNUSED
      IF(ISETKW(71).EQ.0) THEN
        E_CONV = 0.0D0
        ISETKW(71) = 2
      END IF
*
*. 72 : Default Threshold for total wave function correction = 0.0
*UNUSED
      IF(ISETKW(72).EQ.0) THEN
        C_CONV = 0.0D0
        ISETKW(72) = 2
      END IF
*
*. 73 : Perform Class selection : Default if still NO
*
      IF(ISETKW(73).EQ.0) THEN
        ICLSSEL = 0
        ISETKW(73) = 2
      END IF
*
* 74 : Calculation of density matrices : Default is
*       calculation of one-body density
*
      IF(ISETKW(74).EQ.0) THEN
        IDENSI = 0
        ISETKW(74) = 2
      END IF
*. If IDENSI was set to zero and properties were requested
*  overwrite input to obtain 1-el matrix
      IF(IDENSI.EQ.0.AND.ISETKW(80).EQ.1) THEN
        WRITE(6,*) ' You have specified calculation of'
        WRITE(6,*) ' one-electron properties, and this'
        WRITE(6,*) ' requires the calculation of the '
        WRITE(6,*) ' one-electron density. '
        WRITE(6,*)
        WRITE(6,*) ' You have, however, specified IDENSI=0'
        WRITE(6,*) ' corresponding to no densities'
        WRITE(6,*)
        WRITE(6,*) ' I will allow myself to modify your'
        WRITE(6,*) ' input to allow calculation of the '
        WRITE(6,*) ' one-electron densities, so property'
        WRITE(6,*) ' calculation can proceed as planned '
        WRITE(6,*)
        WRITE(6,*)                        ' Lucia '
*. and do it
C       IDENSI = 1
        write(6,*) 'Ok, but I will change IDENSI to 3 '
        write(6,*) 'meaning that spin-dependent densities '
        write(6,*) 'will be calculated. If you do not want '
        write(6,*) 'this, set DENSI to 0 in input.'
        write(6,*) 'After all, one-electron properties cannot '
        write(6,*) 'be determined by LUCIAREL currently anyway.'
        write(6,*)
        write(6,*) '                        Luciarel '
        IDENSI = 3
      END IF
*
* 75 : Perturbation expansion of EKT, default is no
*UNUSED
      IF(ISETKW(75).EQ.0) THEN
        IPTEKT = 0
        ISETKW(75) = 2
      END IF
*
* 76 : Root for zero order operator , default is NROOT
*UNUSED
*. Should be less or equal to NROOT
      IF(ISETKW(76).EQ.1) THEN
        IF(IH0ROOT.GT.NROOT) THEN
          WRITE(6,*) ' Zero order operator root (H0ROOT) larger '
          WRITE(6,*) ' than total number of roots (NROOT) '
          WRITE(6,*) ' CHANGE NROOT or H0ROOT '
          NMISS = NMISS + 1
        END IF
      END IF
      IF(ISETKW(76).EQ.0) THEN
       ISETKW(76) = 2
       IH0ROOT = NROOT
      END IF
*
* 77 : NO restart from previous vectors in calculation 2
*      Default is NO NO, i.e. restart in calculation 2
*
      IF(ISETKW(77).EQ.0) THEN
        IRST2 = 1
        ISETKW(77) = 2
      END IF
*
* 78 : NO skipping  of initial energy evaluations
*
      IF(ISETKW(78).EQ.0) THEN
        ISKIPEI = 0
        ISETKW(78) = 2
      END IF
*
* 79 : Symmetry of x,y,z - needed for property calculations
*
      IF(ISETKW(79).EQ.0) THEN
*. Problematic if Properties should be calculated
       IF(ISETKW(80).EQ.1) THEN
         WRITE(6,*) ' Symmetry of X,Y,Z has not been given'
         WRITE(6,*) ' You have to specify this for property calc'
         WRITE(6,*) ' Please add KEYWORD XYZSYM '
         NMISS = NMISS + 1
         ISETKW(79) = -1
       ELSE
*. Is not needed, just supply zeroes
         DO ICOMP = 1, 3
           IXYZSYM(ICOMP) = 0
         END DO
         ISETKW(79) = 2
       END IF
      END IF
*
* 80 : Property calculation, default is no
*
      IF(ISETKW(80).EQ.0) THEN
        NPROP = 0
        ISETKW(80) = 2
      END IF
*
* 81 : Double group CI, default is still no
*      corresponding to a socalled no-no pair CI calculation
*
      IF(ISETKW(81).EQ.0) THEN
        IDBGCI = 0
      ELSE
*. Double group CI in action
*. Check to ensure that a double group was specified
*. for symmetry
        IF(.NOT.(PNTGRP.GE.5.AND.PNTGRP.LE.9)) THEN
          WRITE(6,*) ' You specified a double group CI calculation'
          WRITE(6,*) ' But you did not give a double group '
          WRITE(6,*)
          WRITE(6,*) ' (DB_C1, DB_CI,DB_CS,DB_C2,DB_C2H ) '
          WRITE(6,*)
          WRITE(6,*) ' as symmetry group. Please correct this'
          NERROR = NERROR + 1
          ISETKW(81) = -1
        END IF
      END IF
*
* 82 : Relativistic calculation : Default is no
*
      IF(ISETKW(82).EQ.0) THEN
        IRELACI = 0
        ISETKW(82) = 2
        WRITE(6,*) ' RELACI flag has not been set.'
        WRITE(6,*) ' I can do non-rel. CI in principle,'
        WRITE(6,*) ' but not this way.'
        WRITE(6,*) ' Set SO-OP to 0 instead.'
        WRITE(6,*)
        WRITE(6,*) '                   Luciarel '
        Call Abend('Quitting')
      END IF
*. If relativistic calculation ensure that DBGCI is turned on
      IF(IRELACI.EQ.1.AND.IDBGCI.EQ.0) THEN
        WRITE(6,*) ' Current code performs only '
        WRITE(6,*) ' relativistic CI in double groups '
        WRITE(6,*)
        WRITE(6,*) ' You have specified the RELACI flag '
        WRITE(6,*) ' but not DBGCI '
        WRITE(6,*)
        WRITE(6,*) ' Please add DBGCI and ensure that '
        WRITE(6,*) ' a double group has been specified as'
        WRITE(6,*) ' symmetrygroup'
        Call Abend('Quitting')
      END IF
*
* 83 : Type of spin orbit operator
*
      IF(ISETKW(83).EQ.0) THEN
*. Must be included for relativistic calculation
        IF(IRELACI.EQ. 1 ) THEN
          WRITE(6,*) ' Type of spin-orbit operator '
          WRITE(6,*) ' ( Keyword SO-OP ) not specified'
          WRITE(6,*) ' This is required for relativistic calc'
          WRITE(6,*)
          WRITE(6,*) ' Please add and return '
          call abend('Quitting')
        END IF
      END IF
*
* 84  MK2REF , no default, required for double group/rel calculation
*
      IF(ISETKW(84).EQ.0) THEN
         IF(IDBGCI.EQ.0.OR.IRELACI.EQ.0) THEN
           MK2REF_CI = 0
           ISETKW(84) = 2
         ELSE
           WRITE(6,*) ' Keyword MK2REF_CI not defined '
           WRITE(6,*) ' This is required for double group CI'
           WRITE(6,*) ' and/or relativistic CI '
           WRITE(6,*) ' Please add MK2REF_CI and return '
           ISETKW(84) = -1
           NERROR = NERROR + 1
         END IF
      END IF
*
* 85  MK2DEL , no default, required for double group/rel calculation
*
      IF(ISETKW(85).EQ.0) THEN
         IF(IDBGCI.EQ.0.OR.IRELACI.EQ.0) THEN
           MK2DEL_CI = 0
           ISETKW(85) = 2
         ELSE
           WRITE(6,*) ' Keyword MK2DEL_CI not defined '
           WRITE(6,*) ' This is required for double group CI'
           WRITE(6,*) ' and/or relativistic CI '
           WRITE(6,*) ' Please add MK2DEL_CI and return '
           ISETKW(85) = -1
           NERROR = NERROR + 1
         END IF
      END IF
*
* 86  Initial Integral transformation, default is no
*
      IF(ISETKW(86).EQ.0) THEN
         ITRANS = 0
         ISETKW(86) = 2
      END IF
*
* 87  Initial spinor-expansion matrix, default is readin - if relativistic
*
      IF(ISETKW(87).EQ.0.AND.ISETKW(88).EQ.0) THEN
        IF(IRELACI.NE.0) THEN
          if (ENVIRO(1:6).ne.'DIRAC ') then
            ISPEXP_INI = 2
            write(6,*) ' Sorry, but spinor expansion has to be read'
            write(6,*) ' from MO expansion. '
            write(6,*) ' Set SP<-MO in input.'
            write(6,*)
            write(6,*) ' SP<-RE true spinor expansion from start '
            write(6,*) ' not possible yet.'
            Call Abend('Quitting')
          else
            ISPEXP_INI = 3
          end if
        ELSE
          ISPEXP_INI = 1
        END IF
        ISETKW(87) = 2
        ISETKW(88) = 2
      END IF
*
* 92 : One- or two-particle (full) Hamiltonian
*
      if (ISETKW(92).eq.0) then
         IHAM12 = 2
         ISETKW(92) = 2
      else
        if (IPARHAM.eq.1) IHAM12 = 1
        if (IPARHAM.eq.2) IHAM12 = 2
        if (IPARHAM.eq.3) IHAM12 = 3
      end if
*
* 93 : Calculation with externally determined Fock matrix
*
      if (ISETKW(93).eq.0) then
         IFOCKEX = 0
         ISETKW(93) = 2
      end if
*
* 94 : MOLCAS orbitals from JOBIPH or INPORB file
*
      if (ISETKW(94).eq.0) then
        ORBFILE = .false.
        ISETKW(94) = 2
      else
        if (INTIMP.ne.1) then
          write(LUOUT,*)
          write(LUOUT,*) ' MOLCAS environment required for '
          write(LUOUT,*) ' INPORB file. Please specify first! '
          Call Abend2( 'Quitting.' )
        end if
      end if
*
* 95 : print flag for wave function information
*
      IF(ISETKW(95).EQ.0) THEN
         IPRHAM_CI = 0
         ISETKW(95) = 2
      END IF
*
* 97 : Selected route for sigma vectors and densities
*
      IF(ISETKW(97).EQ.0) THEN
        SIGDEN_ROUTE = 'TRADIT'
        ISETKW(97) = 2
      else if (ISETKW(97).EQ.1) then
        SIGDEN_ROUTE = 'SIGDEN'
      END IF
*
* 98 : print flag for hamiltonian information
*
      IF(ISETKW(98).EQ.0) THEN
         IPRHAM_CI = 0
         ISETKW(98) = 2
      END IF
*
*
*. Thresholds only active in connection with IDIAG = 2,
*. Check and maybe iisue a warning
      IF(IDIAG.EQ.2) THEN
*. Check to ensure that zero or two thresholds were  set,
        IF(ISETKW(69).NE.ISETKW(70)) THEN
          WRITE(LUOUT,*)
     &    ' Only a single threshold (E_THRE or C_THRE) '
          WRITE(LUOUT,*)
     &    ' on individual determinants given. '
          WRITE(LUOUT,*)
     &    ' One of the thresholds vanishes therefore and '
          WRITE(LUOUT,*)
     &    ' all determinants will therefore be included  '
          WRITE(LUOUT,*)
          WRITE(LUOUT,*) '                   Warns '
          WRITE(LUOUT,*)
          WRITE(LUOUT,*) '                   LUCIA  '
        END IF
      ELSE
*. Good old diagonalization, thrsholds not active
        IF(ISETKW(69).EQ.1.OR.ISETKW(70).EQ.1) THEN
          WRITE(LUOUT,*)
     &    ' Thresholds on selection of individual coefficients '
          WRITE(LUOUT,*)
     &    ' are only active in connection with keyword TERACI  '
          WRITE(LUOUT,*)
          WRITE(LUOUT,*) '                   Warns '
          WRITE(LUOUT,*)
          WRITE(LUOUT,*) '                   LUCIA  '
        END IF
      END IF
*
      IF(ISETKW(50).EQ.1.AND. ISETKW(51).EQ.0) THEN
* Number of GAS shells given but no occupations !!
        WRITE(6,*) ' GAS calculation (GASSH specified)'
        WRITE(6,*) ' But no Occupation constraints (GASSPC) '
        WRITE(6,*)
        WRITE(6,*) ' Please add GASSPC '
        NMISS = NMISS + 1
      END IF
*
*
*
      IF(NMISS.NE.0) THEN
        WRITE(LUOUT,'(A,I9)')
     &  ' Number of missing required keyword ', NMISS
        WRITE(LUOUT,*)
     &  ' You have wounded me I give up '
        WRITE(LUOUT,*)
        WRITE(LUOUT,*)
        WRITE(LUOUT,*)
        WRITE(LUOUT,*)
     & '     An expert is a man who has made all the mistakes,'
        WRITE(LUOUT,*)
     &  '     which can be made, in a very narrow field        '
        WRITE(LUOUT,*)
     &  '                                                      '
        WRITE(LUOUT,*)
     &  '                                      Niels Bohr      '
        Call Abend('Quitting')
      END IF
*
*. Open one-electron file to obtain core energy and
*. dimensions of orbital space
*  This is not used at the moment. Further calls in routine
*  are commented out. Data are determined from input.
      IF(NOINT.EQ.0) THEN
        CALL GET_ORB_DIM_ENV_REL(ECORE_ENV)
        IF(ISETKW(46).EQ.2) ECORE = ECORE_ENV
      ELSE
        WRITE(6,*) ' GETOBS not called since NOINT = 1'
        ECORE = 0.0D0
      END IF

************************************************************
*                                                          *
* Part 3 : Print input                                     *
*                                                          *
************************************************************
*
      WRITE(LUOUT,*)
      WRITE(LUOUT,*) '******************'
      WRITE(LUOUT,*) '*  Title of run  *'
      WRITE(LUOUT,*) '******************'
      WRITE(LUOUT,*)
      CALL PRTITL(TITLEC)
      WRITE(LUOUT,*)
*
*. Machine in use
c     WRITE(6,'(A,A)') '    Machine in use : ', MACHINE
*
*. Type of calculation
*
      WRITE(LUOUT,*)
      WRITE(LUOUT,*) '************************'
      WRITE(LUOUT,*) '*  Type of calculation *'
      WRITE(LUOUT,*) '************************'
      WRITE(LUOUT,*)
      IF(IDBGCI.EQ.0.AND.IRELACI.EQ.0) THEN
        WRITE(LUOUT,'(A)')
     &  '     Non-relativistic calculation in spatial pointgroup'
      ELSE IF(IDBGCI.EQ.1.AND.IRELACI.EQ.0) THEN
        WRITE(LUOUT,'(A)')
     &  '     Non-relativistic calculation in double group '
      ELSE IF(IDBGCI.EQ.1.AND.IRELACI.EQ.1) THEN
        WRITE(LUOUT,'(A)')
     &  '     Relativistic calculation in double group '
      ELSE IF(IDBGCI.EQ.0.AND.IRELACI.EQ.1) THEN
        WRITE(LUOUT,'(A)')
     &  '     Relativistic calculation in spatial pointgroup '
        WRITE(LUOUT,'(A)')
     &  '     This is dangerous. You better think about it! '
      END IF
      IF(IRELACI.EQ.1) THEN
        IF(ISO1.EQ.1.AND.ISO2.EQ.1) then
          if (ENVIRO(1:6).eq.'DIRAC ') then
            WRITE(LUOUT,'(A)')
     &  '         Dirac-Coulomb operator '
          else
            WRITE(LUOUT,'(A)')
     &  '         One- and two- electron spin-orbit operator '
          end if
        end if
        IF(ISO1.EQ.1.AND.ISO2.EQ.0) then
          if (ENVIRO(1:6).ne.'DIRAC ') then
            WRITE(LUOUT,'(A)')
     &  '         One-electron spin-orbit (mean-field) operator '
          else
            WRITE(LUOUT,'(A)')
     &  '         Dirac-Coulomb Hamiltonian        '
          end if 
        end if
        if (ISO1.eq.0.and.ISO2.eq.0.and.ENVIRO(1:6).eq.'DIRAC ') then
          write(LUOUT,'(A)')
     &  '         Dirac-Coulomb Hamiltonian        '
          write(LUOUT,'(A)')
     &  '         Spinfree formalism               '
        end if
        if (ISO1.EQ.0.AND.ISO2.EQ.1) then
          if (ENVIRO(1:6).eq.'DIRAC ') then
            write(6,*) 'Setting of relativistic operators not'
            write(6,*) 'allowed in DIRAC environment.'
            write(6,*) 'This is relativistic theory!'
            call Abend2 ('Input error.')
          end if
          WRITE(LUOUT,'(A)')
     &  '         Two-electron spin-orbit operator '
        end if
      END IF
      if (IFOCKEX.eq.1) then
        WRITE(LUOUT,'(A)')
     &  '         Fock matrix gen. externally for inactive orbitals '
        WRITE(LUOUT,'(A43,8I4)')
     &  '         Full orbital basis in symmetries: ',
     &            (NORB_F(ISYM),ISYM=1,NSMSH,1)
      end if
      if (INTIMP.eq.1) then
        if (ORBFILE) then
          WRITE(LUOUT,'(A)')
     &    '         MO vector read from formatted ASCII file INPORB '
      else
          WRITE(LUOUT,'(A)')
     &    '         MO vector read from RASSCF file JOBIPH          '
        end if
      end if
*
*. Type of reference state
      WRITE(LUOUT,*)
      WRITE(LUOUT,*) '*************************************'
      WRITE(LUOUT,*) '*  Symmetry and spin of CI vectors  *'
      WRITE(LUOUT,*) '*************************************'
      WRITE(LUOUT,*)
*. Point group
      IF(PNTGRP.EQ.1) THEN
        WRITE(LUOUT,'(A)')
     &  '     Point group ............ D2H'
      ELSE IF(PNTGRP.EQ.2) THEN
        WRITE(LUOUT,'(A)')
     &  '     Point group ............ C inf v'
      ELSE IF(PNTGRP.EQ.3) THEN
        WRITE(LUOUT,'(A)')
     &  '     Point group ............ D inf h'
      ELSE IF(PNTGRP.EQ.4) THEN
        WRITE(LUOUT,'(A)')
     &  '     Point group ............ O3'
      ELSE IF(PNTGRP.EQ.5) THEN
        WRITE(LUOUT,'(A)')
     &  '     Double group ........... C1'
      ELSE IF(PNTGRP.EQ.6) THEN
        WRITE(LUOUT,'(A)')
     &  '     Double group ........... Ci'
      ELSE IF(PNTGRP.EQ.7) THEN
        WRITE(LUOUT,'(A)')
     &  '     Double group ........... Cs'
      ELSE IF(PNTGRP.EQ.8) THEN
        WRITE(LUOUT,'(A)')
     &  '     Double group ........... C2'
      ELSE IF(PNTGRP.EQ.9) THEN
        WRITE(LUOUT,'(A)')
     &  '     Double group ........... C2H'
      END IF
*.Spatial symmetry
      IF(PNTGRP.EQ.1) THEN
        WRITE(LUOUT,'(A,I1)')
     &  '     Spatial symmetry ....... ', IREFSM
      ELSE IF(PNTGRP.EQ.2) THEN
        WRITE(LUOUT,'(A,I1)')
     &  '     ML value ............... ', IREFML
      ELSE IF(PNTGRP.EQ.3) THEN
        WRITE(LUOUT,'(A,I1)')
     &  '     ML value ............... ', IREFML
        IF(IREFPA.EQ.1) WRITE(LUOUT,'(A)')
     &  '     Parity   ..............  Gerade'
        IF(IREFPA.EQ.2) WRITE(LUOUT,'(A)')
     &  '     Parity   ..............  Ungerade'
      ELSE IF(PNTGRP.EQ.4) THEN
        WRITE(LUOUT,'(A,I1)')
     &  '     L  value ............... ', IREFL
        WRITE(LUOUT,'(A,I1)')
     &  '     ML value ............... ', IREFML
        IF(IREFPA.EQ.1) WRITE(LUOUT,'(A)')
     &  '     Parity   ..............  Gerade'
        IF(IREFPA.EQ.2) WRITE(LUOUT,'(A)')
     &  '     Parity   ..............  Ungerade'
       ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
        WRITE(LUOUT,'(A,I1)')
     &  '     Double group symmetry... ', IREFSM
      END IF
*.Spin
      WRITE(LUOUT,*)
      IF(IRELACI.EQ.0.AND.IDBGCI.EQ.0) THEN
        WRITE(LUOUT,'(A,I2)')
     &    '     2 times spinprojection  ', MS2
c       WRITE(LUOUT,'(A,I2)')
c    &    '     Spin multiplicity ....  ', MULTS
      ELSE
        WRITE(LUOUT,'(A,I2)')
     &    '     2*reference value of MK ', MK2REF_CI
        WRITE(LUOUT,'(A,I2,A,I2)')
     &    '     Allowed interval of MK2 ',
     &          MK2REF_CI-MK2DEL_CI,' to ', MK2REF_CI+MK2DEL_CI
        WRITE(6,*)
      END IF
*.Number of active electrons
      WRITE(LUOUT,'(A,I2)')
     &  '     Active electrons .....  ', NACTEL
      WRITE(LUOUT,*)
      WRITE(LUOUT,*) '*********************************************'
      WRITE(LUOUT,*) '*  Shell spaces and occupation constraints  *'
      WRITE(LUOUT,*) '********************************************* '
      WRITE(LUOUT,*)
*
      IF(IDOGAS.EQ.0) THEN
*
*. NOT a GAS expansion, NOT IN USE !!!!!
*
      WRITE(LUOUT,'(A,10I4)')
     &  '                Irrep ',(I,I = 1,NSMSH)
      WRITE(LUOUT,'(A,2X,10A)')
     &  '                ===== ',('====',I = 1,NSMSH)
*
*. Inactive
      IF(ISETKW(7).EQ.1) THEN
        WRITE(LUOUT,'(A,10I4)')
     &  '     Inactive         ',(NINASH(I),I=1,NSMSH)
      END IF
*. Core
      IF(ISETKW(8).EQ.1) THEN
        WRITE(LUOUT,'(A,10I4)')
     &  '     Core             ',(NRS0SH(I,1),I=1,NSMSH)
      END IF
*. RAS1
      IF(ISETKW(9).EQ.1) THEN
        WRITE(LUOUT,'(A,10I4)')
     &  '     Ras1             ',(NRSSH(I,1),I=1,NSMSH)
      END IF
*. RAS2/ACTIVE
      IF(ISETKW(10).EQ.1) THEN
        IF(INTSPC.EQ.1) THEN
          WRITE(LUOUT,'(A,10I4)')
     &  '     Active           ',(NRSSH(I,2),I=1,NSMSH)
        ELSE IF(INTSPC.EQ.2) THEN
          WRITE(LUOUT,'(A,10I4)')
     &  '     Ras2             ',(NRSSH(I,2),I=1,NSMSH)
        END IF
      END IF
*. RAS3
      IF(ISETKW(11).EQ.1) THEN
        WRITE(LUOUT,'(A,10I4)')
     &  '     Ras3             ',(NRSSH(I,3),I=1,NSMSH)
      END IF
      IF(INTSPC.EQ.2.AND.IMLCR3.EQ.1) WRITE(LUOUT,'(A)')
     &  '     ( RAS 3 space supplied by courtesy of TRAONE )'
*. Secondary space
      IF(ISETKW(13).EQ.1) THEN
        DO 310 ITP = 1,MXR4TP
          WRITE(LUOUT,'(A,I2,A,10I4)')
     &  '     Secondary',ITP,'      ',(NRS4SH(I,ITP),I=1,NSMSH)
  310   CONTINUE
      END IF
*. Deleted space
      IF(ISETKW(26).EQ.1) THEN
        WRITE(LUOUT,'(A,10I4)')
     &  '     Deleted          ',(NDELSH(I),I=1,NSMSH)
      END IF
      IF(IMLCR3.EQ.2) WRITE(LUOUT,'(A)')
     &  '     ( Deleted shells supplied by courtesy of TRAONE )'
*. Occupation constraints

      IF(INTSPC.EQ.2) THEN
*. Zero order space
        WRITE(LUOUT,*)
        IFRST = 1
        IF(ISETKW(9).EQ.1) THEN
         WRITE(LUOUT,'(A,I2)')
     &  '     Zero order space, Minimal occupation in Ras1 .. ',MNRS10
         IFRST = 0
        END IF
        IF(ISETKW(11).EQ.1) THEN
          IF(IFRST.EQ.1) THEN
            WRITE(LUOUT,'(A,I2)')
     &  '     Zero order space, Maximal occupation in Ras3 .. ',MXRS30
          ELSE
            WRITE(LUOUT,'(A,I2)')
     &  '                       Maximal occupation in Ras3 .. ',MXRS30
          END IF
        END IF
*.Reference space
        WRITE(LUOUT,*)
        IFRST = 1
        IF(ISETKW(9).EQ.1) THEN
         WRITE(LUOUT,'(A,I2)')
     &  '     Reference  space, Minimal occupation in Ras1 .. ',MNRS1R
         IFRST = 0
        END IF
        IF(ISETKW(11).EQ.1) THEN
          IF(IFRST.EQ.1) THEN
            WRITE(LUOUT,'(A,I2)')
     &  '     Reference  space, Maximal occupation in Ras3 .. ',MXRS3R
          ELSE
            WRITE(LUOUT,'(A,I2)')
     &  '                       Maximal occupation in Ras3 .. ',MXRS3R
          END IF
        END IF
      END IF
*.Core space
      WRITE(LUOUT,*)
      IF(ISETKW(8).EQ.1) THEN
        WRITE(LUOUT,'(A,I2)')
     &  '     Largest number of excitations out of core .....   ',MXHR0
      END IF
*.Secondary space
      IF(ISETKW(13).EQ.1) THEN
        WRITE(LUOUT,'(A,I2)')
     &  '     Largest number of excitations to secondary space  ',MXER4
      END IF

      ELSE
*
*. GAS space
*
      WRITE(LUOUT,*)
      WRITE(LUOUT,*) ' *************************'
      WRITE(LUOUT,*) ' Generalized active space '
      WRITE(LUOUT,*) ' *************************'
      WRITE(LUOUT,*)
      IF(IRELACI.EQ.0) THEN
        WRITE(LUOUT,'(A)') ' Orbital subspaces:'
        WRITE(LUOUT,'(A)') ' ================== '
      ELSE
        WRITE(LUOUT,'(A)')
     &  ' Orbital subspaces in molecular pointgroup'
        WRITE(LUOUT,'(A)')
     &  ' ========================================='
      END IF
      WRITE(LUOUT,'(A,10I4,A)')
     &  '                Irrep ',(I,I = 1,NSMSH)
      WRITE(LUOUT,'(A,2X,10A,A)')
     &  '                ===== ',('====',I = 1,NSMSH )
      DO IGAS = 1, NGAS
        WRITE(LUOUT,'(A,I2,A,10I4,6X,2I6)')
     &  '        GAS',IGAS,'          ',
     &  (NGSSH(IRREP,IGAS),IRREP = 1, NSMSH )
      END DO
*
      IF(IRELACI.NE.0) THEN
        WRITE(6,*)
     &' ( Information about spinors in double group comes later)'
      END IF
*
      WRITE(LUOUT,'(A,I3)')
     & ' Number of Occupation spaces : ',NCISPC
      WRITE(LUOUT,*)
      DO ICISPC = 1, NCISPC
      WRITE(LUOUT,'(A,I3)')
     &' Bounds on accumulated occupations for subspace : ',ICISPC
      WRITE(LUOUT,'(A)')
     & ' ====================================================== '
      WRITE(LUOUT,'(A)')
      WRITE(LUOUT,'(A)') '         Min. occ    Max. occ '
      WRITE(LUOUT,'(A)') '         ========    ======== '
      DO IGAS = 1, NGAS
        WRITE(LUOUT,'(A,I2,3X,I3,9X,I3)')
     &  '   GAS',IGAS,IGSOCCX(IGAS,1,ICISPC),IGSOCCX(IGAS,2,ICISPC)
      END DO
      END DO
*
      WRITE(6,*)
      WRITE(6,'(A,I3)')
     &' Number of CI spaces included : ', NCMBSPC
      WRITE(6,*)
      DO JCMBSPC = 1, NCMBSPC
        WRITE(6,*)
        WRITE(6,'(A,I3)') ' Information about CI space ', JCMBSPC
        WRITE(6,'(A)')    ' =================================='
        WRITE(6,'(A,I3)')
     &  ' Number of GAS spaces included ', LCMBSPC(JCMBSPC)
        WRITE(6,'(A,10I3)')
     &  ' GAS spaces included           ',
     &    (ICMBSPC(II,JCMBSPC),II=1,LCMBSPC(JCMBSPC))
*
        WRITE(6,'(A,I3)') ' Number of calculations in this CI space ',
     &  NSEQCI(JCMBSPC)
        WRITE(6,'(A)')   '  Calculations in this subspace '
        WRITE(6,'(A)')   '  =============================='
        DO JSEQ = 1, NSEQCI(JCMBSPC)
          IF(ISEQCI(JSEQ,JCMBSPC).LT.0.AND.ISEQCI(JSEQ,JCMBSPC).GT.-5)
     &    THEN
            WRITE(6,'(A,I3)')
     &      '       Vector free calculation at level ',
     &      -ISEQCI(JSEQ,JCMBSPC)
          ELSE IF( ISEQCI(JSEQ,JCMBSPC).GT.0 ) THEN
           WRITE(6,'(A,I3)')
     &      '       Normal CI with max. iterations = ',
     &      ISEQCI(JSEQ,JCMBSPC)
          ELSE IF (ISEQCI(JSEQ,JCMBSPC).EQ.-5 ) THEN
           WRITE(6,'(A,I3)')
     &      '       Perturbation calculation '
          END IF
        END DO
      END DO
      WRITE(6,*)
      END IF
*
      IF(XLAMBDA.NE.1.0D0) THEN
        WRITE(6,*)
        WRITE(6,'(A,F13.8)')
     &  ' Modifed operator H(l) = l*F + l*(H-F) used with l =',XLAMBDA
        IF(IRESTR.EQ.0) THEN
        WRITE(6,*)
     &  ' Notice : This madness starts  in second calculation'
        ELSE
         WRITE(6,*) ' You have specified a calculation with modified '
         WRITE(6,*) ' Hamiltonian (the LAMBDA option) and RESTART '
         WRITE(6,*) ' so this is what I will do '
         WRITE(6,*)
         WRITE(6,*) '   1:) Perform CI in space 1 to obtain Hamiltonian'
         WRITE(6,*) '       (no RESTART in this space )'
         WRITE(6,*) '   2:) CI calculation in space 2  with '
         WRITE(6,*) '       modified Hamiltonian and RESTART from LU21'
         WRITE(6,*) ' Space 2 should therefore correspond to the'
         WRITE(6,*) ' restarted calculation '
       END IF
      END IF
*
      WRITE(LUOUT,*)
      WRITE(LUOUT,*) '***********'
      WRITE(LUOUT,*) '*  Roots  *'
      WRITE(LUOUT,*) '*********** '
      WRITE(LUOUT,*)
      WRITE(LUOUT,'(A,I3)')
     &  '     Number of roots to be obtained ', NROOT
      WRITE(LUOUT,'(A,(20I3))')
     &  '     Roots to be obtained ', (IROOT(I),I=1, NROOT )
*
      WRITE(LUOUT,*)
      WRITE(LUOUT,*) '**************************'
      WRITE(LUOUT,*) '*  Run time definitions  *'
      WRITE(LUOUT,*) '************************** '
      WRITE(LUOUT,*)
*. Program environment
      WRITE(6,'(A,A6)') '      Program environment... ', ENVIRO
*. Integral transformation
      IF(ITRANS.NE.0) THEN
        WRITE(LUOUT,'(A,I1)')
     &  '     Initial integral transformation, level = ', ITRANS
      END IF
*. Initial spinor transformation
      IF(ISPEXP_INI.EQ.1) THEN
        WRITE(LUOUT,'(A)')
     &  '     Initial spinor expansion obtained from MO expansion'
      ELSE IF (ISPEXP_INI.EQ.2) THEN
        WRITE(LUOUT,'(A)')
     &  '     Initial spinor expansion read in '
      else if (ISPEXP_INI.eq.3) then
        WRITE(LUOUT,'(A)')
     &  '     Reading integrals over transformed spinors. '
      END IF
*. No. of particle Hamiltonian
      if (IPARHAM.eq.1) then
        WRITE(LUOUT,'(A)')
     &  '     Only 1-particle Hamiltonian under consideration '
      else if (IPARHAM.eq.2) then
        WRITE(LUOUT,'(A)')
     &  '     Full 1- and 2-particle Hamiltonian active       '
      else if (IPARHAM.eq.3) then
        WRITE(LUOUT,'(A)')
     &  '     Only 2-particle Hamiltonian under consideration '
      end if
*. Integral import
      IF(NOINT.EQ.1) THEN
        WRITE(LUOUT,'(A)')
     &  '     No integrals will be read in       '
      ELSE IF(NOINT.EQ.0) THEN
      IF(INTIMP.EQ.1) THEN
        WRITE(LUOUT,'(A)')
     &  '     Integrals read in in MOLCAS format '
      ELSE IF(INTIMP.EQ.2) THEN
        WRITE(LUOUT,'(A)')
     &  '     Integrals read in in LUCAS format '
      ELSE IF(INTIMP.EQ.3) THEN
        WRITE(LUOUT,'(A)')
     &  '     Integrals read in in formatted form (E22.15) ',
     &  '      From unit 13'
        WRITE(LUOUT,'(A)')
     &  '     All integrals of correct symmetry combination read in'
      ELSE IF(INTIMP.EQ.5) THEN
        WRITE(LUOUT,'(A)')
     &  '     Integrals read in in DALTON/SIRIUS format '
      else if (INTIMP.eq.6) then
        WRITE(LUOUT,'(A)')
     &  '     Integrals read in in DIRAC format '
        WRITE(LUOUT,'(A)')
     &  '       Dirac-Coulomb Hamiltonian curr. assumed '
      else if (INTIMP.eq.8) then
        WRITE(LUOUT,'(A)')
     &  '     Integrals read in in MOLECULE-SWEDEN format '
        WRITE(LUOUT,'(A)')
     &  '     one-electron ints from unit 2, transformation '
        WRITE(LUOUT,'(A)')
     &  '     spin-orbit ints from SOMFINT_*, trafo '
        WRITE(LUOUT,'(A)')
     &  '     two-electron ints from unit 13, transformed '
      END IF
*. Integral storage
      IF(INCORE.EQ.1) WRITE(LUOUT,'(A)')
     &  '     All integrals stored in core'
      END IF
      WRITE(LUOUT,*)
* ( END IF for NOINT
*. CSF or SD expansion
      IF(NOCSF.EQ.0) THEN
        WRITE(LUOUT,'(A)')
     &  '     CI optimization performed with CSF''s '
      ELSE
        WRITE(LUOUT,'(A)')
     &  '     CI optimization performed with SD''s '
      END IF
*. Ms,Ml combinations
      IF(ISETKW(27).EQ.1) THEN
        WRITE(LUOUT,'(A,F8.3)')
     &  '     Spin combinations used with sign ',PSSIGN
      END IF
      IF(ISETKW(28).EQ.1) THEN
        WRITE(LUOUT,'(A,F8.3)')
     &  '     ML   combinations used with sign ',PLSIGN
      END IF
*. Initial approximation to vectors
      WRITE(LUOUT,*)
      IF(IRESTR.EQ.1.AND.IRESTRF.EQ.0) THEN
         WRITE(LUOUT,'(A)')
     &  '     Restarted calculation '
      ELSE IF(IRESTRF.EQ.1) THEN
         WRITE(LUOUT,'(A)')
     &  '     Restarted calculation from REFERENCE space expansion'
      ELSE
         IF(MXP1.NE.0) THEN
           WRITE(LUOUT,'(A)')
     &  '     Initial vectors obtained from explicit Hamiltonian'
         ELSE IF(MXP1.EQ.0) THEN
           WRITE(LUOUT,'(A)')
     &  '     Initial vectors obtained from diagonal'
         END IF
      END IF
*. Handling of degeneracies of initial vectors
      IF(INIDEG.EQ.1) THEN
        WRITE(LUOUT,'(A)')
     &  '     Symmetric combination of degenerate initial vectors'
      ELSE IF (INIDEG.EQ.-1) THEN
        WRITE(LUOUT,'(A)')
     &  '     Antisymmetric combination of degenerate initial vectors'
      ELSE IF (INIDEG.EQ.0) THEN
        WRITE(LUOUT,'(A)')
     &  '     No combination of degenerate initial vectors'
      END IF
*. CI storage mode
      WRITE(6,*)
      if (ISETKW(97).eq.1) then
        WRITE(6,*)
     &  '     SIGDEN route for sigma vectors and densities invoked   '
      end if
      IF(ICISTR.EQ.1) THEN
        WRITE(6,*)
     &  '     3 symmetry blocks and two vectors will be held in core '
      ELSE IF( ICISTR.EQ.2) THEN
        WRITE(6,*)
     &  '     3 symmetry blocks will be held in core '
      ELSE IF( ICISTR.EQ.3) THEN
        WRITE(6,*)
     &  '     3 symmetry-occ-occ blocks will be held in core '
      END IF
      IF(LCSBLK.NE.0) WRITE(6,'(A,I10)')
     &  '      Smallest allowed size of sigma- and C-batch ',LCSBLK
      WRITE(LUOUT,'(A,I4)')
     &  '     Dimension of block of resolution strings ', MXINKA
      IF(ICJKAIB.EQ.1) THEN
        WRITE(6,*)
     &  '     C(Ka,Jb,j) matrices will be used in alpha-beta loop '
      END IF
      IF(MOCAA.NE.0) WRITE(LUOUT,'(A,I4)')
     &  '     MOC method used for alpha-alpha+beta-beta loop '
      IF(MOCAB.NE.0) WRITE(LUOUT,'(A,I4)')
     &  '     MOC method used for alpha-beta loop            '
*
      WRITE(LUOUT,*)
      IF(IDENSI.EQ.0) THEN
        WRITE(LUOUT,'(A)')
     &  '     No calculation of density matrices  '
      ELSE IF(IDENSI.EQ.1) THEN
        WRITE(LUOUT,'(A)')
     &  '     One-body density matrix calculated '
      ELSE IF(IDENSI.EQ.2) THEN
        WRITE(LUOUT,'(A)')
     &  '     One- and two-body density matrices  calculated '
      else if (IDENSI.eq.3) then
        WRITE(LUOUT,'(A)')
     &  '     Spin-dependent one-body density matrix calculated '
      else if (IDENSI.eq.4) then
        WRITE(LUOUT,'(A)')
     &  '     Spin-dependent one- and two-body densities
     &        calculated '
      else
        WRITE(LUOUT,'(A,I15)')
     &  '     Unknown code for density matrices :',IDENSI
      END IF
*
*
*. Diagonalization information
      WRITE(LUOUT,'(A)')
     &  '     CI diagonalization : '
      WRITE(LUOUT,'(A)')
     &  '     ==================== '
*. Subspace Hamiltonian
      IF(MXP1+MXP2+MXQ .EQ.0) THEN
        WRITE(LUOUT,'(A)')
     &  '        No subspace Hamiltonian '
      ELSE
        WRITE(LUOUT,'(A,3I4)')
     &  '        Dimensions of subspace Hamiltonian ',MXP1,MXP2,MXQ
      END IF
*. Diagonalizer
      IF(IDIAG.EQ.1.AND.ICISTR.EQ.1) THEN
        WRITE(LUOUT,'(A)')
     &    '        Diagonalizer : MINDV4 '
      ELSE IF(IDIAG.EQ.1.AND.ICISTR.GE.2) THEN
        WRITE(LUOUT,'(A)')
     &    '        Diagonalizer : MICDV4 '
      ELSE IF(IDIAG.EQ.1.AND.ICISTR.GE.3) THEN
        WRITE(LUOUT,'(A)')
     &    '        Diagonalizer : CMICDV '
      ELSE IF(IDIAG.EQ.2) THEN
      WRITE(LUOUT,'(A)')
     &  '        Diagonalizer : PICO2  '
      END IF
*. No restart in CI calc 2
      IF(IRST2.EQ.0) THEN
      WRITE(LUOUT,'(A)')
     &  '        No restart from previous vectors in second calc '
      END IF
      IF(ISKIPEI.EQ.1) THEN
      WRITE(LUOUT,'(A)')
     &  '        Initial energy evaluations skipped after first calc'
      WRITE(LUOUT,'(A)')
     &  '        (Only active in connection with TERACI )'
      END IF
*
*
*
*. Number of iterations
      WRITE(LUOUT,'(A,I2)')
     &  '        Allowed number of iterations    ',MAXIT
*. Number of CI vectors in subspace
      WRITE(LUOUT,'(A,I2)')
     &  '        Allowed Dimension of CI subspace',MXCIV_CI
*
      WRITE(LUOUT,'(A,E11.5)')
     &  '        Convergence threshold for energy',THRES_E
      WRITE(6,*)
      IF(IDIAG.EQ.2) THEN
        WRITE(LUOUT,'(A,E11.5)')
     &   '        Individual second order energy threshold ',E_THRE
        WRITE(LUOUT,'(A,E11.5)')
     &   '        Individual first order wavefunction threshold ',C_THRE
        IF(ICLSSEL.EQ.1) THEN
         WRITE(LUOUT,*)
         WRITE(LUOUT,'(A)')
     &   '         Class selection will be performed : '
         WRITE(LUOUT,'(A)')
     &   '         =================================== '
         WRITE(LUOUT,'(A,E11.5)')
     &    '          Total second order energy threshold ',E_CONV
         WRITE(LUOUT,'(A,E11.5)')
     &    '          Total first order wavefunction threshold ',C_CONV
        ELSE
         WRITE(LUOUT,'(A)')
     &'            No class selection in iterative procedure '
        END IF
      END IF
*
      IF(IPERT.NE.0) THEN
        WRITE(LUOUT,'(A)')
     &    '     Perturbation calculation'
        WRITE(LUOUT,'(A)')
     &  '     ======================= '
        WRITE(6,*)
     &  '        Root Choosen as zero order state ', IRFROOT
        WRITE(6,*)
     &  '        Root used for zero order operator ', IH0ROOT
        IF(IE0AVEX.EQ.1) THEN
          WRITE(6,*)
     &  '        Expectation value of H0 used as zero order energy '
        ELSE IF( IE0AVEX.EQ.2) THEN
          WRITE(6,*)
     &  '        Exact energy of reference used as zero order energy'
        END IF
        WRITE(6,*)
     &  '        Correction vectors obtained through  order ', NPERT
        IF(IH0SPC.EQ.0) THEN
        WRITE(6,*)
     &  '        No restrictions on perturbation interactions '
        ELSE
        WRITE(6,*)
     &  '        Perturbation restricted to interactions in subspaces'
        END IF
*
        WRITE(6,*)
        IF(IH0SPC.NE.0) THEN
        WRITE(6,*)
     &  '        Number of perturbation subspaces ', NPTSPC
        WRITE(6,*)
        WRITE(6,*)
     &  '        ======================== '
        WRITE(6,*)
     &  '        Perturbation subspaces : '
        WRITE(6,*)
     &  '        ======================== '
        DO JPTSPC = 1, NPTSPC
COLD      WRITE(LUOUT,'(A)')
COLD &     ' ====================================================== '
          WRITE(LUOUT,'(A)')
          WRITE(LUOUT,'(7X,A)') '         Min. occ    Max. occ '
          WRITE(LUOUT,'(7X,A)') '         ========    ======== '
          DO IGAS = 1, NGAS
            WRITE(LUOUT,'(7X,A,I2,3X,I3,9X,I3)')
     &      '   GAS',IGAS,IOCPTSPC(1,IGAS,JPTSPC)
     &                   ,IOCPTSPC(2,IGAS,JPTSPC)
          END DO
        END DO
*
        WRITE(6,*)
        WRITE(6,'(7X,A)') ' ========================================'
        WRITE(6,'(7X,A)') ' Approximate Hamiltonian in CI subspaces '
        WRITE(6,'(7X,A)') ' ========================================'
        WRITE(6,'(7X,A)')
        WRITE(6,'(7X,A)') '    Subspace          H(apr)   '
        WRITE(6,'(7X,A)') '  ============================='
        WRITE(6,'(7X,A)')
        DO JPTSPC = 1, NPTSPC
          IF(IH0INSPC(JPTSPC).EQ.1) THEN
            WRITE(LUOUT,'(12X,I3,8X,A)')
     &      JPTSPC, ' Diagonal Fock operator '
          ELSE IF(IH0INSPC(JPTSPC).EQ.2) THEN
            WRITE(LUOUT,'(12X,I3,8X,A)')
     &      JPTSPC, ' Epstein-Nesbet operator'
          ELSE IF(IH0INSPC(JPTSPC).EQ.3) THEN
            WRITE(LUOUT,'(12X,I3,8X,A)')
     &      JPTSPC, ' Nondiagonal Fock operator '
          ELSE IF(IH0INSPC(JPTSPC).EQ.4) THEN
            WRITE(LUOUT,'(12X,I3,8X,A)')
     &      JPTSPC, ' Complete Hamiltonian  '
          ELSE IF(IH0INSPC(JPTSPC).EQ.5) THEN
            WRITE(LUOUT,'(12X,I3,8X,A)')
     &      JPTSPC, ' Mix of Fock and Exact operator '
          END IF
         END DO
         IF(ISETKW(61).GT.0) THEN
           WRITE(6,*)
           WRITE(6,'(7X,A)')
     &     ' Orbital subspaces where exact Hamiltonian is used : '
           WRITE(6,'(7X,A)')
     &      '===================================================='
           WRITE(6,*)
           WRITE(LUOUT,'(10X,10(2X,I3))') (IH0EXSPC(I),I=1, NH0EXSPC)
           WRITE(6,*)
         END IF
*
       END IF
       END IF
*
       WRITE(6,*)
       IF(NPROP.EQ.0) THEN
         WRITE(6,*) '        No calculation of properties'
       ELSE
         WRITE(6,'(7X,A,I3)')
     &   ' Number of properties to be calculated', NPROP
         WRITE(6,*)
         WRITE(6,'(9X,A)')    ' Properties : '
         WRITE(6,'(9X,A)')   ' ============='
         DO IPROP = 1, NPROP
           WRITE(6,'(16X,A)') PROPER(IPROP)
         END DO
       END IF
*
       WRITE(6,*)
       WRITE(6,*)
       IF(IEXTKOP.EQ.0) THEN
         WRITE(6,'(7X,A)') ' No extended Koopmans'' calculations '
       ELSE
         WRITE(6,'(7X,A)') ' Extended Koopmans'' calculations '
       END IF
*
      IF(NOMOFL.EQ.0) THEN
        WRITE(6,*)
        WRITE(6,'(7X,A)') ' Final orbitals :'
        WRITE(6,'(7X,A)') ' ================'
        WRITE(6,*)
        IF(IFINMO.EQ.1) THEN
          WRITE(6,'(10X,A)') ' Natural orbitals'
        ELSE IF(IFINMO.EQ.2) THEN
          WRITE(6,'(10X,A)') ' Canonical orbitals'
        ELSE IF(IFINMO.EQ.3) THEN
          WRITE(6,'(10X,A)') ' Pseudo-natural orbitals'
          WRITE(6,'(10X,A)')
     &   ' (Density matrix diagonalized in orbital subspaces )'
        ELSE IF(IFINMO.EQ.4) THEN
          WRITE(6,'(10X,A)') ' Pseudo-canonical orbitals'
          WRITE(6,'(10X,A)')
     &   ' (FI+FA  diagonalized in orbital subspaces )'
         ELSE IF (IFINMO .EQ. 5 ) THEN
          WRITE(6,'(10X,A)')
     &   ' Pseudo-natural-canonical orbitals (sic)'
          WRITE(6,'(10X,A)')
     &   ' (Pseudo natural orbitals are first obtained'
          WRITE(6,'(10X,A)')
     &   '  by diagonalizing density matrix in orbital subpspaces.'
          WRITE(6,'(10X,A)')
     &   '  FI+FA is transformed to this basis, and the transformed'
          WRITE(6,'(10X,A)')
     &   '  matrix is block diagonalized) '
          WRITE(6,*)
          WRITE(6,'(10X,A)')
     &   ' Orbital spaces in which transformed FIFA is diagonalized'
          WRITE(6,'(10X,A)')
     &   ' ========================================================'
          DO IPSSPC = 1, NPSSPC
            WRITE(LUOUT,'(A,I2,A,10I4,6X,2I6)')
     &      '     SPACE',IPSSPC,'          ',
     &     (NPSSH(IRREP,IPSSPC),IRREP = 1, NIRREP)
          END DO
        END IF
      END IF
*
*. Print levels
*
      WRITE(LUOUT,*)
      WRITE(LUOUT,'(A)')  '     Print levels : '
      IF(ISETKW(29).EQ.2) THEN
        WRITE(LUOUT,'(A,I3)')
     &  '      Default print level for string    information = ', IPRSTR
      ELSE
        WRITE(LUOUT,'(A,I3)')
     &  '      Raised  print level for string    information = ', IPRSTR
      END IF
      IF(ISETKW(95).EQ.2) THEN
        WRITE(LUOUT,'(A,I3)')
     &  '      Default print level for wave function info    = ', IPRWFC
      ELSE
        WRITE(LUOUT,'(A,I3)')
     &  '      Raised  print level for wave function info    = ', IPRWFC
      END IF
      IF(ISETKW(30).EQ.2) THEN
        WRITE(LUOUT,'(A,I3)')
     &  '      Default print level for CI space  information = ', IPRCIX
      ELSE
        WRITE(LUOUT,'(A,I3)')
     &  '      Raised  print level for CI space  information = ', IPRCIX
      END IF
      IF(ISETKW(31).EQ.2) THEN
        WRITE(LUOUT,'(A,I3)')
     &  '      Default print level for orbital   information = ', IPRORB
      ELSE
        WRITE(LUOUT,'(A,I3)')
     &  '      Raised  print level for orbital   information = ', IPRORB
      END IF
      IF(ISETKW(65).EQ.2) THEN
        WRITE(LUOUT,'(A,I3)')
     &  '      Default print level for density matrix        = ', IPRDEN
      ELSE
        WRITE(LUOUT,'(A,I3)')
     &  '      Raised  print level for density matrix        = ', IPRDEN
      END IF
      IF(ISETKW(32).EQ.2) THEN
        WRITE(LUOUT,'(A,I3)')
     &  '      Default print level for iterative information = ', IPRDIA
      ELSE
        WRITE(LUOUT,'(A,I3)')
     &  '      Raised  print level for iterative information = ', IPRDIA
      END IF
      IF(ISETKW(36).EQ.2) THEN
        WRITE(LUOUT,'(A,I3)')
     &  '      Default print level for External blocks       = ', IPRXT
      ELSE
        WRITE(LUOUT,'(A,I3)')
     &  '      Raised  print level for External blocks       = ', IPRXT
      END IF
      IF(ISETKW(89).EQ.2) THEN
        WRITE(LUOUT,'(A,I3)')
     &  '      Default print level for transformation info   = ', IPRTRA
      ELSE
        WRITE(LUOUT,'(A,I3)')
     &  '      Raised  print level for transformation info   = ', IPRTRA
      END IF
      IF(ISETKW(90).EQ.2) THEN
        WRITE(LUOUT,'(A,I3)')
     &  '      Default print level for sigma vector info   = ', IPRSIG
      ELSE
        WRITE(LUOUT,'(A,I3)')
     &  '      Raised  print level for sigma vector info   = ', IPRSIG
      END IF
*
      IF(IPROCC.NE.0) THEN
        WRITE(LUOUT,'(A,I3)')
     &  '      Occupation of SD''s/ configurations           = ',IPROCC
      END IF
*
      IF(MOLCS.EQ.1) WRITE(LUOUT,'(A,E18.9)')
     &  '      Core energy : ', ECORE
*
      IF(IDMPIN.EQ.1) THEN
        WRITE(LUOUT,'(A)')
        WRITE(6,*) ' Integrals written in formatted form (E22.15)'
        WRITE(6,*) ' on file 90 '
      END IF
*
      NIRR_DG_LOC = NIRREP
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE WRTHD_REL(LUOUT)
*
* Introduce LUCIAREL
*
      WRITE(6,'(8X,A)')
     &'****************************************************************'
      WRITE(6,'(8X,A)')
     &'*                                                              *'
      WRITE(6,'(8X,A)')
     &'*           Welcome to LUCIAREL, a direct relativistic         *'
      WRITE(6,'(8X,A)')
     &'*                                CI program                    *'
      WRITE(6,'(8X,A)')
     &'*                                                              *'
      WRITE(6,'(8X,A)')
     &'* LUCIA written by Jeppe Olsen, University of Lund, Sweden     *'
      WRITE(6,'(8X,A)')
     &'* Relativistic CI version by Timo Fleig and Jeppe Olsen        *'
      WRITE(6,'(8X,A)')
     &'*    - Double group symmetry implemented: JO and TF            *'
      WRITE(6,'(8X,A)')
     &'*    - Time-reversal symmetry implemented: JO and TF           *'
      WRITE(6,'(8X,A)')
     &'*    - Integral transformation section (one-el, spin-orbit): TF*'
      WRITE(6,'(8X,A)')
     &'*    - Integral transformation section (two-el): JO and TF     *'
      WRITE(6,'(8X,A)')
     &'*                                                              *'
      WRITE(6,'(8X,A)')
     &'* V 2.0                                                        *'
      WRITE(6,'(8X,A)')
     &'* Version of November 2000                                     *'
      WRITE(6,'(8X,A)')
     &'*   1.0       R-GAS-CI                                         *'
      WRITE(6,'(8X,A)')
     &'*   1.1       MOLCAS version                                   *'
      WRITE(6,'(8X,A)')
     &'*   1.2       External Fock matrices from MOTRA in use         *'
      WRITE(6,'(8X,A)')
     &'*   1.3       Density matrices can be obtained                 *'
      WRITE(6,'(8X,A)')
     &'*             (Implement C_1, C_i, and C_2h)                   *'
      WRITE(6,'(8X,A)')
     &'*    1.3.1     Integral transformation section generalized     *'
      WRITE(6,'(8X,A)')
     &'*   1.4       Resolution matrix batching active                *'
      WRITE(6,'(8X,A)')
     &'*   1.5       Analysis of optimized CI vectors possible        *'
      WRITE(6,'(8X,A)')
     &'*   1.6       Quaternionic and complex double groups active    *'
      WRITE(6,'(8X,A)')
     &'*   1.7       MOLCAS 5.0 version (input section changed)       *'
      WRITE(6,'(8X,A)')
     &'*    1.7.1     Abend calls for parallel runs included          *'
      WRITE(6,'(8X,A)')
     &'*    1.7.2     SEWARD-AMFI calls implemented for MOLCAS4.3 - 5.*'
      WRITE(6,'(8X,A)')
     &'*   1.8       DIRAC 4-component version                        *'
      WRITE(6,'(8X,A)')
     &'*    1.8.1     Spinfree formalism (Dyall, 1993) implemented    *'
      WRITE(6,'(8X,A)')
     &'*    1.8.2     Full spin-dependent formalism part. implemented *'
      WRITE(6,'(8X,A)')
     &'*                    ------------                              *'
      WRITE(6,'(8X,A)')
     &'*   2.0       Version for relativistic MCSCF                   *'
      WRITE(6,'(8X,A)')
     &'*              General sigma and density (from CC) implemented *'
      WRITE(6,'(8X,A)')
     &'****************************************************************'
*
      WRITE(6,*)
      WRITE(6,'(8X,A)')
     &' In case of trouble please contact : '
      WRITE(6,*)
      WRITE(6,'(12X,2A)')
     &' Timo Fleig, Theoretical Chemistry, Faculty of Sciences '
      WRITE(6,'(12X,2A)')
     &' Free University of Amsterdam                   '
      WRITE(6,'(12X,2A)')
     &' De Boelelaan 1083, NL-1081 HV Amsterdam, NL    '
      WRITE(6,'(12X,2A)')
     &' Telephone : +31 20 44 47627 '
      WRITE(6,'(12X,2A)')
     &' E-mail : fleig@chem.vu.nl '
      WRITE(6,*)
     &' or : '
      WRITE(6,*)
      WRITE(6,'(12X,2A)')
     &' Jeppe Olsen, Teoretisk Kemi, Aarhus Universitet '
      WRITE(6,'(12X,2A)')
     &' Langelandsgade 140, DK-8000 Aarhus C '
      WRITE(6,'(12X,2A)')
     &' E-mail : jeppe@signe.teokem.lu.se  or'
      WRITE(6,'(12X,2A)')
     &'          jeppe@aausgi.uni-c.dk'
      WRITE(6,*)
      WRITE(6,'(8X,A)')
      write(6,*)
      write(6,*)
      write(6,*) '-----------------------------------'
      write(6,*) '|  Stability on systems:          |'
      write(6,*) '-----------------------------------'
      write(6,*)
      write(6,*) ' SGI (Octane, Origin 2000) '
      write(6,*) ' IBM (RS6000)              '
      write(6,*) ' REDHAT LINUX              '
      write(6,*) ' SuSE Linux                '
      write(6,*)
      write(6,*)
*
      RETURN
      END
