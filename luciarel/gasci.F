      SUBROUTINE GASCI(F1,F2,INTSPACE,ISM,ISPC,IPRNT,EREF,CIRUN)
*
* CI optimization in GAS space number ISPC for symmetry ISM
*
* Jeppe Olsen, Winter of 1995
*              Winter of 1998 : Modified for relativistic calculations
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      EXTERNAL MV7
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
*
#include "cgas.inc"
#include "lucinp.inc"
#include "ctcc.inc"
*
      COMMON/CSFBAS/KDFTP,KCFTP,KDTOC,KICONF(MXCNSM),KTPFCN(MXCNSM),
     &              KICTS(MXCNSM),KSCTS(MXCNSM),KCNFCN(MXCNSM)
*
*. Common block for communicating with sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
      COMMON/SPINFO/MULTSP,MS2P,
     &              MINOP,MAXOP,NTYP,NDPCNT(MXPCTP),NCPCNT(MXPCTP),
     &              NCNATS(MXPCTP,MXPCSM),NDTASM(MXPCSM),NCSASM(MXPCSM),
     &              NCNASM(MXPCSM)
*./CECORE/
      COMMON/CECORE/ECORE
*
      dimension SCR2(MXROOT)
      dimension F1(*),F2(*)
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GASCI ')
*
      NTESTL = 0010
      NTEST = MAX(NTESTL,IPRNT)
*
      IF(NTEST.GE.3) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================================='
        WRITE(6,*) ' Control has been transferred to GASCI'
        WRITE(6,*) ' ====================================='
        WRITE(6,*)
        WRITE(6,'(A)') '  A few pertinent data : '
        WRITE(6,*)
        WRITE(6,'(A,I2)') '  CI space         ',ISPC
        WRITE(6,*)
        WRITE(6,*) ' Number of GAS spaces included ',LCMBSPC(ISPC)
        WRITE(6,'(A,10I3)') ' GAS spaces included           ',
     &               (ICMBSPC(II,ISPC),II=1,LCMBSPC(ISPC))
        WRITE(6,*)
        WRITE(6,*) ' Occupation constraints : '
        WRITE(6,*) '========================= '
        WRITE(6,*)
        WRITE(6,*)
        DO JJGASSPC = 1, LCMBSPC(ISPC)
         JGASSPC = ICMBSPC(JJGASSPC,ISPC)
        WRITE(6,*)
     &  ' Gas space  Min acc. occupation Max acc. occupation '
        WRITE(6,*)
     &  ' ================================================== '
        DO IGAS = 1, NGAS
          WRITE(6,'(3X,I2,13X,I3,16X,I3)') IGAS,
     &     IGSOCCX(IGAS,1,JGASSPC),IGSOCCX(IGAS,2,JGASSPC)
        END DO
        END DO
*
*
       END IF
*
      NDET = XISPSM(ISM,ISPC)
*. real or complex eigenvectors
      IF(IRELACI.EQ.0) THEN
*. real
        IRC = 1
      ELSE
*. Well, complex
        IRC = 2
      END IF
      WRITE(6,*) ' Number of determinants/combinations  ',NDET
*
*.Transfer to CANDS
      ICSM = ISM
      ISSM = ISM
      ICSPC = ISPC
      ISSPC = ISPC
*
      IF(NOCSF.EQ.1) THEN
        NVAR = NDET
      ELSE
        NVAR = NCSASM(ISM)
      END  IF
      IF(IPRNT.GE.5)
     &WRITE(6,*) '  NVAR in GASCI ', NVAR
*
*. Allocate memory for diagonalization
      IF(ICISTR.EQ.1) THEN
        LBLOCK = NDET
      ELSE IF (ICISTR.EQ.2) THEN
        LBLOCK = MXSB
      ELSE IF (ICISTR.EQ.3) THEN
        LBLOCK = MXSOOB
      END IF
*
      LBLOCK = MAX(LBLOCK,LCSBLK)
      CALL MEMMAN(KVEC1,LBLOCK,'ADDS  ',2,'VEC1  ')
      CALL MEMMAN(KVEC2,LBLOCK,'ADDS  ',2,'VEC2  ')
*. Partition CI vector - for all MS2 spaces into batches
* (*_MS2 arrays are in CSTATE )
      CALL Z_BLKFO(IDC,NMS2VAL,
     &             ISPC,ISM,KLCLBT,KLCLEBT,KLI1BT,KLCIBT,KLCBLTP,
     &             NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,NBAT_MS2,IBAT_MS2)
C          Z_BLKFO(ISPC,ISM,KPCLBT,KPCLEBT,
C    &             KPCI1BT,KPCIBT,KPCBLTP,NBATCH,NBLOCK,
C    &             NBLK_MS2,IBLK_MS2,NBAT_MS2,IBAT_MS2)
*
*. Class divisions of dets
      IF(ICLSSEL.EQ.1) THEN
        WRITE(6,*) ' Info for class selection constructed'
        WRITE(6,*) ' Not implemented for relativistic version '
        Call Abend2
     & ( ' GASCI : CLSSEL not implemented for relativistic version' )
      END IF
*. Allocate memory : 2 vector blocks + 1 vector/scratch blocks
      I_USE_SBLOCK=1
      IF(IDIAG.EQ.2.OR.I_USE_SBLOCK.EQ.1) THEN
        CALL GET_3BLKS(KVEC1,KVEC2,KVEC3,LSCR2)
*. LSCR2 : Will be used for resolution matrices, save
         LCJRES = LSCR2
      END IF
*
*. CI diagonal - if required
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++
      if (CIRUN.eq.'DIAG  '.or.CIRUN.eq.'CIINII'.or.
     &    CIRUN.eq.'SIGMA '.or.CIRUN.eq.'SIGMAD') then
*
        if (CIRUN.ne.'SIGMA ') then
*+++++++++++++++++++++++++++++++++++++++++++++++++++
*
      IF(IDIAG.EQ.2) THEN
        LUDIA = LUSC1
      END IF
      WRITE(6,*) ' Diagonal will be written to unit ',LUDIA
      IF(.NOT.(IDIAG.EQ.2.AND.IRESTR.EQ.1)) THEN
        ECOREP = ECORE
        IF(ICISTR.GE.2) CALL REWINO(LUDIA)
        IRECO = 1
        CALL GASDIAT(WORK(KVEC1),LUDIA,ECOREP,ICISTR,IHAM12,
     &               WORK(KLCBLTP),NBLOCK,WORK(KLCIBT),
     &               KPICI1,INTSPACE,IDENSI,NIRREP,ISPINFREE,ENVIRO,
     &               SIGDEN_ROUTE,IRECO)
*. End of real part
        CALL ITODS(-1,1,-1,LUDIA)
        IF(IRC.EQ.2) THEN
* Diagonal once again- for the imaginary part
          IRECO = 2
          CALL GASDIAT(WORK(KVEC1),LUDIA,ECOREP,ICISTR,IHAM12,
     &                 WORK(KLCBLTP),NBLOCK,WORK(KLCIBT),
     &                 KPICI1,INTSPACE,IDENSI,NIRREP,ISPINFREE,ENVIRO,
     &                 SIGDEN_ROUTE,IRECO)
*. End of imaginary  part
        CALL ITODS(-1,1,-1,LUDIA)
        END IF
        IF(IPRCIX.GE.2) WRITE(6,*) ' Diagonal constructed  '
      ELSE
         WRITE(6,*) ' Diagonal not calculated '
      END IF
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++
        end if
        if (CIRUN.ne.'DIAG  ') then
*+++++++++++++++++++++++++++++++++++++++++++++++++++
*
*. Transfer control to optimization routine
*
      MINST = 1
      IF(IRESTR.EQ.0) THEN
        INICI = 0
      ELSE
        INICI = -1
      END IF
      NPRDET = 0
      LBLK = -1
*. This is a CI and not a perturbation calculation
      IPERTOP = 0
*. real or complex eigenvectors
      IF(IRELACI.EQ.0) THEN
*. real
        IRC = 1
      ELSE
*. Well, complex
        IRC = 2
      END IF
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++
          if (CIRUN.eq.'SIGMAD'.or.CIRUN.eq.'SIGMA ') then
*+++++++++++++++++++++++++++++++++++++++++++++++++++
            call sigden_ctrl(WORK(KVEC1),WORK(KVEC2),LUC,LUHC,
     &                       WORK(KT_CC),1)
            print*,'C VECTOR'
            call wrtvcdc(WORK(KVEC1),LUC,1,-1)
            print*,'SIGMA VECTOR'
            call wrtvcdc(WORK(KVEC2),LUHC,1,-1)
          else
*
      CALL CIEIG5(MV7,INICI,EROOT,WORK(KVEC1),WORK(KVEC2),
     &     MINST,LUDIA,LUC,LUHC,LUSC1,LUSC2,LUSC3,LUSC34,LUSC35,
     &     NVAR,NBLK,NROOT,MXCIV,MAXIT,MXROOT,
     &     LUCIVI,IPRNT,WORK(KSBEVC),
     &     NPRDET,WORK(KH0),WORK(KSBIDT),
     &     MXP1,MXP2,MXQ,WORK(KH0SCR),SCR2,
     &     WORK(KT_CC),
     &     ECORE,ICISTR,LBLK,
     &     IDIAG,WORK(KVEC3),THRES_E,
     &     NBATCH,WORK(KLCLBT),WORK(KLCLEBT),
CTF  &     WORK(KLCI1BT),
     &     WORK(KLCIBT),INIDEG,E_THRE,C_THRE,
     &     E_CONV,C_CONV,ICLSSEL,
CTF  &     WORK(KLBLKCLS),
     &     NOCCLS,
CTF  &     WORK(KLCLSC),WORK(KLCLSE), WORK(KLCLSCT),WORK(KLCLSET),
CTF  &     WORK(KLCLSA),WORK(KLCLSL),WORK(KLBLKA),
CTF  &     WORK(KLCLSD),WORK(KLCLSDT),ISKIPEI,WORK(KLC2B),WORK(KLCLSA2),
     &     LBLOCK,IRC,IDCOMH,IHAM12,SIGDEN_ROUTE)
      EREF = EROOT(IRFROOT)
CTF   WRITE(6,*) ' Memcheck after CIEIG5'
CTF   CALL LMEMCHK
*+++++++++++++++++++++++++++++++++++++++++++++++++++
          end if
        end if
      else if (CIRUN.eq.'ANAL  ') then
* Analysis of optimized wave function
*+++++++++++++++++++++++++++++++++++++++++++++++++++
*
      if (IPRWFC.ge.1) then
        call rewino(LUC)
        do JROOT=1,NROOT,1
          call rewino(LUSC61)
          call copvcdc(LUC,LUSC61,WORK(KVEC2),0,LBLK)
          call rewino(LUSC61)
          WRITE(6,*)
          WRITE(6,'(1H ,A,I3)')
     &    ' ***************************************************'
          WRITE(6,'(1H ,A,I3)')
     &    ' Analysis of Optimized Wave Function for ROOT = ',JROOT
          WRITE(6,'(1H ,A,I3)')
     &    ' ***************************************************'
*
          call gasana(WORK(KVEC1),NBLOCK,WORK(KLCIBT),
     &                WORK(KLCBLTP),LUSC61,LUSC37,
     &                KVEC3)
        end do
      end if
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++
      else if (CIRUN.eq.'DENS1 '.or.CIRUN.eq.'DENS2 ') then
* Analysis of density and occupation
*+++++++++++++++++++++++++++++++++++++++++++++++++++
      LBLK = -1
      if (IDENSI.ne.0) then
        if (IDENSI.eq.1.or.IDENSI.eq.2) then
          write(6,*) 'IDENSI is ',IDENSI
          write(6,*) 'I am not allowed to calculate non-relativistic'
          write(6,*) 'density matrices YET.'
          write(6,*) 'See readin.f or ask the experts for instructions'
          Call Abend2( 'LUCIAREL quitting' )
        end if
*
* One-body density matrix
* Spin-dependent first-order density matrix
*
C?      write(6,*) 'LUC is  ',LUC
C?      write(6,*) 'LUHC is  ',LUHC
C?      write(6,*) 'LUSC1 is  ',LUSC1
        CALL REWINO(LUC)
        CALL REWINO(LUHC)
        CALL REWINO(LUSC1)
*
*  Start the loops over approximate eigenvectors
        do JROOT = 1,NROOT,1
          WRITE(6,*)
          WRITE(6,'(1H ,A,I3)')
     &    ' **************************************************'
          WRITE(6,'(1H ,A,I3)')
     &    ' Analysis of Density and occupation for ROOT = ',JROOT
          WRITE(6,'(1H ,A,I3)')
     &    ' **************************************************'
CTF
*  New density routines for relativistic (and somewhen also nonrel.)
*  calculations. DBG_RHO is called and vectors are read from disk
*  internally.
*
          if (IRELACI.eq.1) then
            INONREL = 0
          else
            INONREL = 1
          end if
          if (IDENSI.eq.3) then
            IDEN12 = 1
            write(6,*)
            write(6,*) '   Computing spin-dependent one-particle '
            write(6,*) '   density'
          else if (IDENSI.eq.4) then
            IDEN12 = 2
            write(6,*)
            write(6,*) '  Computing spin-dependent'
            write(6,*) '  one- and two-particle densities.'
            write(6,*) '  I hope you have set SIGDEN_ROUTE = SIGDEN!'
          else
            write(6,*) 'IDENSI =  ',IDENSI
            Call Abend2( 'Problemo in gasci. Quitting.' )
          end if
*
* Pick current root from multi-root vector on LUC
* (I guess SKPVCD for skipping a given number of vectors on
*  file should be used here for convenience ...)
          if (JROOT.ge.2) then
C
C  MCSCF: C vectors are placed to files outside in KR_LUCI_IO.
C
c           call copvcdc(LUC,LUHC,WORK(KVEC2),1,LBLK)
c           do IRT=1,JROOT-2,1
c             call copvcdc(LUC,LUHC,WORK(KVEC2),0,LBLK)
c           end do
c           rewind(LUSC61)
c           rewind(LUSC62)
c           call copvcdc(LUC,LUSC61,WORK(KVEC2),0,LBLK)
c           call copvcdc(LUSC61,LUSC62,WORK(KVEC2),1,LBLK)
c         else if (JROOT.eq.1) then
c           call copvcdc(LUC,LUSC61,WORK(KVEC2),1,LBLK)
c           call copvcdc(LUC,LUSC62,WORK(KVEC2),1,LBLK)
          end if
*
          if (IPRDEN.ge.500) then
            write(6,*) 'Current vector ',JROOT,'  on unit ',LUSC61
            call wrtvcdc(WORK(KVEC1),LUSC61,1,LBLK)
            write(6,*) 'Current vector ',JROOT,'  on unit ',LUSC62
            call wrtvcdc(WORK(KVEC2),LUSC62,1,LBLK)
          end if
*
          if (SIGDEN_ROUTE.eq.'TRADIT') then
            call dbg_rho(WORK(KVEC1),WORK(KVEC2),INONREL,IDEN12,JROOT)
          else if (SIGDEN_ROUTE.eq.'SIGDEN') then
* Erase T amplitudes from sigma calculation
            ZERO = 0.D0
            call setvec(WORK(KT_CC),ZERO,NINT_CC)
*
* Compute density matrix elements
            call sigden_ctrl(WORK(KVEC1),WORK(KVEC2),LUSC61,LUSC62,
     &                       WORK(KT_CC),2)
*
            if (IPRDEN.ge.10) then
              write(6,*) '======================================='
              write(6,'(A,I4)') ' Density matrix for root ',JROOT
              write(6,*) '======================================='
              IRC_D = 2
              call wrtdens(WORK(KT_CC),WORK(KINDEX),WORK(KLIHIND),
     &                     WORK(KSIGNNHX),WORK(KLOP_REO),IRC_D)
            end if
*
* Transfer density matrices to KRMC format
            IRC_D = 2
            call trnsfdens(F1,F2,WORK(KT_CC),
     &                     WORK(KINDEX),WORK(KLABEXTP),WORK(KLIHIND),
     &                     WORK(KSIGNNHX),WORK(KLOP_REO),WORK(KSIOPREO),
     &                     IRC_D,NTOOB2,IPRDEN)
*
          else
            write(6,*) 'SIGDEN_ROUTE with unallowed value in gasci.'
            write(6,*) 'SIGDEN_ROUTE ',SIGDEN_ROUTE
          end if
*
CTF     WRITE(6,*)
CTF     IF(ICISTR.EQ.1) THEN
CTF       CALL FRMDSC(WORK(KVEC1),NVAR,-1,LUC,IMZERO,IAMPACK)
CTF       CALL COPVEC(WORK(KVEC1),WORK(KVEC2),NVAR)
CTF       CALL DENSI2(IDENSI,WORK(KRHO1),
CTF  &           WORK(KRHO2),
CTF  &           WORK(KVEC1),WORK(KVEC2),0,0)
CTF     ELSE
CTF       IF(IDENSI.GE.1.OR.NROOT.GT.1) THEN
CTF          CALL REWINO(LUSC1)
CTF          CALL COPVCD(LUC,LUSC1,WORK(KVEC1),0,LBLK)
CTF          CALL COPVCD(LUSC1,LUHC,WORK(KVEC1),1,LBLK)
CTF       END IF
CTF       CALL DENSI2(IDENSI,WORK(KRHO1),
CTF  &                WORK(KRHO2),
CTF  &                WORK(KVEC1),WORK(KVEC2),LUSC1,LUHC)
C         IF(NROOT.GT.1) THEN
C           CALL REWINO(LUSC1)
C           CALL GASANA(WORK(KVEC1),NBLOCK,WORK(KLCIBT),WORK(KLCBLTP),
C    &                  LUSC1,ICISTR)
C         ELSE
C           CALL REWINO(LUC)
C           CALL GASANA(WORK(KVEC1),NBLOCK,WORK(KLCIBT),WORK(KLCBLTP),
C    &                  LUC,ICISTR)
C           END IF
CTF       END IF
C         IF(NPROP.GT.0) CALL ONE_EL_PROP
        END DO
      end if
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++
      else
        write(6,*) 'CIRUN has unallowed value :',CIRUN
        call quit('Stopping in GASCI')
      end if
*     ^ CIRUN cases for MCSCF calls.
*+++++++++++++++++++++++++++++++++++++++++++++++++++
*
      print*,'I have completed my task'
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUM,'GASCI')
      RETURN
      END
