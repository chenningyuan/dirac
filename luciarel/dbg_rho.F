      subroutine dbg_rho(C,HC,INONREL,IDEN12,JROOT)
*
* Outer routine for generation of density matrix
* for double group CI.
*
* The two vectors are supposed to reside on LUSC62 (C) and LUSC61 (S)
* C and HC are work space partitions called WORK(KVEC1) and
* WORK(KVEC2), respectively, above (in gasci).
*
* Generalized to allow for relativistic densities, T.F. + J.O., may 98
* cont'd ... fall 1998, T.F.
* cont'd ... august 1999, Timo Fleig
* cont'd ... december 1999, Timo Fleig
* cont'd ... june 2000, Timo Fleig
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*.Definition of c and sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC

#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "oper.inc"
#include "symm.inc"
#include "clunit.inc"
*
      common/glbref/ISPC,ISM,ITERGLB
*
      character*3 FILESTAT
*
      dimension EIGVDEN(MXNSPIPDMGRP)
      dimension NGSOB_HE(MXPIRR*MXPNGAS),
     &          NGSOB2_HE(MXPIRR*MXPNGAS)
*
      NTESTL = 100
      NTEST = max(NTESTL,IPRDEN)
*
      CALL QENTER('DBG_RH')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DBG_RH')
*
*. Arrays for partitioning of sigma
* (MS2 arrays are in CSTATE )
      CALL Z_BLKFO(IDC,NMS2VAL,
     &             ISPC,ISM,KLSLBT,KLSLEBT,KLSI1BT,KLSIBT,
     &             KLSBLTP,NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,
     &             NBAT_MS2,IBAT_MS2)
*
*  Erase density matrices
      do IC=1,2,1
        do IB=1,2,1
          do IA=1,2,1
            do IER=1,NTOOB**2,1
              WORK(KRHO1_DBG(IA,IB,IC,JROOT)+IER-1) = 0.D0
            end do
          end do
        end do
      end do
*
      CALL DASSG3(WORK(KRHO1_DBG(2,2,1,JROOT)),
     &            WORK(KRHO1_DBG(2,2,2,JROOT)),
     &            WORK(KRHO1_DBG(2,1,1,JROOT)),
     &            WORK(KRHO1_DBG(2,1,2,JROOT)),
     &            WORK(KRHO1_DBG(1,1,1,JROOT)),
     &            WORK(KRHO1_DBG(1,1,2,JROOT)),
     &            WORK(KRHO1_DBG(1,2,1,JROOT)),
     &            WORK(KRHO1_DBG(1,2,2,JROOT)),
     &            C,HC,WORK(KVEC3),NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &            WORK(KLSI1BT),WORK(KLSIBT),
     &            INONREL,IDEN12,IFINMO,IPRDEN)
*
*    WORK(KRHO1_DBG(A,B,C,D)), A=1,2  first  index barred,unbarred
*                              B=1,2  second index barred, unbarred
*                              C=1,2  real, imaginary part of density
*                              D=JROOT density for root JROOT
*
*  Check determined density matrices
      if (NTEST.ge.00) then
       write(6,*) 'Checking calculated density matrices for this root:'
       do I=1,2,1
c       I = 2
        do J=1,2,1
         do K=1,2,1
          write(6,*)
          write(6,*) 'Checking address on WORK array:'
          write(6,'(A10,I1,A1,I1,A1,I1,A1,I2,A2,I9)')
     &            'KRHO1_DBG(',I,',',J,',',K,',',JROOT,')=',
     &             KRHO1_DBG(I,J,K,JROOT)
          write(6,*)
          if (I.eq.1.and.J.eq.1.and.K.eq.1) write(6,*) 'BB den, real'
          if (I.eq.1.and.J.eq.1.and.K.eq.2) write(6,*) 'BB den, imag'
          if (I.eq.1.and.J.eq.2.and.K.eq.1) write(6,*) 'BU den, real'
          if (I.eq.1.and.J.eq.2.and.K.eq.2) write(6,*) 'BU den, imag'
          if (I.eq.2.and.J.eq.1.and.K.eq.1) write(6,*) 'UB den, real'
          if (I.eq.2.and.J.eq.1.and.K.eq.2) write(6,*) 'UB den, imag'
          if (I.eq.2.and.J.eq.2.and.K.eq.1) write(6,*) 'UU den, real'
          if (I.eq.2.and.J.eq.2.and.K.eq.2) write(6,*) 'UU den, imag'
          call wrtmat(WORK(KRHO1_DBG(I,J,K,JROOT)),
     &                NTOOB,NTOOB,NTOOB,NTOOB)
          write(6,*)
         end do
        end do
       end do
      end if
*
* Arrange density matrices for diagonalization
* and do the job. Write natural spinors to files.
*
      if (JROOT.eq.1) then
*  Currently we determine natural spinors only for the
*  lowest root. On the long run, the spinors should be
*  averaged over a selection of states (components).
        if (INATIT.eq.1) then
          FILESTAT = 'NEW'
        else if (INATIT.ge.2) then
          FILESTAT = 'OLD'
        end if
        open(unit=LUMOUTR,file='NATSPIR',
     &       status='UNKNOWN',form='FORMATTED')
        open(unit=LUMOUTI,file='NATSPII',
     &       status='UNKNOWN',form='FORMATTED')
        call rewino(LUMOUTR)
        call rewino(LUMOUTI)
        ISUM = 0
        ISTORE = 0
        do ISYM=1,NIRR_PN,1
          do IGS=1,NGAS,1
            ISUM = ISUM + NGSSH(ISYM,IGS)
          end do
          if (ISUM.gt.ISTORE) ISTORE = ISUM
          ISUM = 0
        end do
        NFHELEN = ISTORE**2
        ISTODIM = (2*ISTORE)**2
c       ISTOTRI = ((2*ISTORE)**2 + 2*ISTORE) / 2
        ISTOTRI = (ISTORE**2 + ISTORE) / 2
        call memman(KNBLK,NFHELEN,'ADDL',2,'TRHEMA')
        call memman(KREADEN,ISTODIM*2,'ADDL',2,'READEN')
        call memman(KREATRI,ISTOTRI*2,'ADDL',2,'REATRI')
        if (NTEST.ge.20) then
          write(6,'(A24,I8,A4)')
     &            'Alloc. full help matrix ',NFHELEN,' R*8'
          write(6,'(A24,I8,A4)')
     &            'Alloc. storing dimension',ISTODIM,' R*8'
          write(6,'(A24,I8,A4)')
     &            'Alloc. storing triangle ',ISTOTRI,' R*8'
        end if
*
        do ISYM=1,NIRR_PN,1
          ISPSM = ISYM + NIRR_PN
          if (NTEST.ge.10) then
            write(6,*) 
            write(6,*) ' ========================================'
            write(6,*) '  Determining natural spinors, symm.',ISPSM
            write(6,*) ' ========================================'
          end if
          if (ISYM.le.2) then
            IGRP = 1
            if (DOUGRP.eq.7.and.ISYM.eq.2) IGRP = 2
          else if (ISYM.ge.3) then
            IGRP = 2
          end if
          ITRIALEN = (SYMDIM(IGRP)**2 + SYMDIM(IGRP)) / 2
          ISCRLEN = (2*SYMDIM(IGRP)+1)*SYMDIM(IGRP)
          call memman(KNATSPI,(SYMDIM(IGRP)*2)**2,
     &                'ADDL',2,'NATSPI')
          call memman(KNSPIR,SYMDIM(IGRP)**2,'ADDL',2,'NSPINR')
          call memman(KNSPII,SYMDIM(IGRP)**2,'ADDL',2,'NSPINI')
          call memman(KSCRATCH,ISCRLEN,'ADDL',2,'SCRALE')
          if (NTEST.ge.20) then
            write(6,'(A24,I8,A4)')
     &           'Alloc. nat. spinor array',SYMDIM(IGRP)*2,' R*8'
            write(6,'(A24,I8,A4)')
     &           'Alloc. real NS array    ',SYMDIM(IGRP)**2,' R*8'
            write(6,'(A24,I8,A4)')
     &           'Alloc. imag. NS array   ',SYMDIM(IGRP)**2,' R*8'
          end if
*
          do IRI=1,2,1
            call arr_dens(WORK(KRHO1_DBG(2,2,IRI,JROOT)),
     &                    WORK(KRHO1_DBG(2,1,IRI,JROOT)),
     &                    WORK(KRHO1_DBG(1,2,IRI,JROOT)),
     &                    WORK(KRHO1_DBG(1,1,IRI,JROOT)),
     &                    NGSOB,NGSOB2,
     &                    NGSOB_HE,NGSOB2_HE,
     &                    WORK(KNBLK),
     &                    WORK(KREADEN+(IRI-1)*ISTODIM),
     &                    DOUGRP,SYMDIM(IGRP),
     &                    NTOOB,ISPSM,NIRR_PN,IGRP,
     &                    ISPSF_MO,ITOTDM,NGAS,IPRDEN)
            call fold(WORK(KREATRI+(IRI-1)*ISTOTRI),
     &                WORK(KREADEN+(IRI-1)*ISTODIM),
     &                SYMDIM(IGRP),1,ITRIALEN)
          end do
*
*  Diagonalize one-particle density matrix
*
          if (NTEST.ge.100) then
            write(6,*)
            write(6,*) 'Folded real density matrix:'
            call wrtmat(WORK(KREATRI),ISTOTRI,1,ISTOTRI,1)
            write(6,*)
            write(6,*) 'Folded imaginary density matrix:'
            call wrtmat(WORK(KREATRI+ISTOTRI),ISTOTRI,1,
     &                  ISTOTRI,1)
          end if
          call ceigen(WORK(KREATRI),WORK(KREATRI+ISTOTRI),
     &                SYMDIM(IGRP),EIGVDEN,WORK(KNATSPI),
     &                WORK(KSCRATCH))
*
          if (NTEST.ge.100) then
            write(6,*) '--------------------------------'
            write(6,*) 'Writing output natural spinors.'
            write(6,*) '--------------------------------'
            call wrtmat(WORK(KNATSPI),(SYMDIM(IGRP)*2)**2,1,
     &                   (SYMDIM(IGRP)*2)**2,1)
            write(6,*) '--------------------------------'
            write(6,*) '--------------------------------'
          end if
*
          if (NTEST.ge.50) then
            write(6,*) 'Natural spinor occupation numbers'
            write(6,*) 'Spinor symmetry: ',ISPSM
            call wrtmat(EIGVDEN,SYMDIM(IGRP),1,SYMDIM(IGRP),1)
          end if
*
*  Write this symmetry to file LUMOUT
*
          call wrt_natspi(WORK(KNATSPI),SYMDIM(IGRP),
     &                    LUMOUTR,LUMOUTI,ISPSM)
*
          write(6,*) 'Checking memory in SPISYM loop ...'
          call lmemchk
*
        end do
*       ^ End loop over spinor symmetries
*
        close (unit=LUMOUTR)
        close (unit=LUMOUTI)
*
      end if
*     ^ End if which root is being treated (to be modified)
*
*. Eliminate local memory
      CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'DBG_RH')
*
      CALL QEXIT('DBG_RH')
*
      RETURN
      END
