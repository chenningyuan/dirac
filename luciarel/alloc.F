      SUBROUTINE ALLOC
*
* Dimensions and
* Allocation of static memory
*
* =====
* Input
* =====
*
* KFREE : Pointer to first element of free space
* Information in /LUCINP/,/ORBINP/,/CSYM/
*
* ======
* Output
* ======
* KFREE : First array of free space after allocation of
*         static memory
* /GLBBAS/,/CDIM/
*
*
* =======
* Version
* =======
*
* Modified Jan 1997
* SIGDEN_ROUTE SIGDEN added by Timo Fleig, May-June 2001

*. Input
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cstate.inc"
#include "csm.inc"
#include "crun.inc"
#include "cprnt.inc"
#include "cintfo.inc"
#include "ctcc.inc"
*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

*.Output
#include "glbbas.inc"
CTF adding global info
      common/glbref/ISPC/ISM/ITERGLB
*.1 : One electron integrals( Complete matrix allocated )
      CALL MEMMAN(KINT1,NTOOB ** 2,'ADDS  ',2,'INT1   ')
*.1.1 : Inactive fock matrix
      CALL MEMMAN(KFI  ,NTOOB ** 2,'ADDS  ',2,'FI     ')
*.2 : Two electron integrals
      IF(NOINT.EQ.0.AND.INCORE.EQ.1) THEN
         CALL MEMMAN(KINT2,NINT2,'ADDS  ',2,'INT2   ')
      END IF
*. Pointers to symmetry block of integrals
      CALL MEMMAN(KPINT1,NBINT1,'ADDS  ',2,'PINT1   ')
      CALL MEMMAN(KPINT2,NBINT2,'ADDS  ',2,'PINT2   ')
*. Symmetry of last index as a function of initial index
      CALL MEMMAN(KLSM1,NBINT1,'ADDS  ',2,'LSM1   ')
      CALL MEMMAN(KLSM2,NBINT2,'ADDS  ',2,'LSM2   ')
*
* =============================
*. For double group integrals
* =============================
* Caution!
* Double group calculation can also be non-relativistic!
* (Well, at least technically!)
      IF(IDBG.EQ.1) THEN
        if (SIGDEN_ROUTE.eq.'TRADIT') then
*. 1e integrals
          DO IL1 = 1, NL1D
*. Real part
            CALL MEMMAN(KINT1D(1,IL1),NTOOB ** 2,'ADDS  ',2,'INT1D ')
*. Imaginary part
            CALL MEMMAN(KINT1D(2,IL1),NTOOB ** 2,'ADDS  ',2,'INT1D ')
*. and a pointer
            CALL MEMMAN(KPINT1D(IL1),NBINT1,'ADDS  ',2,'PINT1D')
          END DO
          CALL MEMMAN(KPINT1,NBINT1,'ADDS  ',2,'PINT1   ')
          CALL MEMMAN(KPINT2,NBINT2,'ADDS  ',2,'PINT2   ')
*. 2e integrals
          DO IL2 = 1, NL2D
*. Real part
            LENGTH = L2D(IL2)
C?          write(6,*) 'Current LENGTH is  ',LENGTH,'IL2 is  ',IL2
            CALL MEMMAN(KINT2D(1,IL2),LENGTH,'ADDS  ',2,'INT2D ')
*. Imaginary part
            CALL MEMMAN(KINT2D(2,IL2),LENGTH,'ADDS  ',2,'INT2D ')
*. and a pointer
            CALL MEMMAN(KPINT2D(IL2),NBINT2,'ADDS  ',2,'PINT2D')
          END DO
*
*  Allocation of density matrices
*  See readin.f for explanations
          if (IDENSI.eq.1.or.IDENSI.eq.2) then
            call MEMMAN(KRHO1,NTOOB ** 2,'ADDS  ',2,'RHO1  ')
            if (IDENSI.eq.2) then
              LRHO2 = NTOOB**2*(NTOOB**2+1)/2
              call MEMMAN(KRHO2,LRHO2     ,'ADDS  ',2,'RHO2  ')
              write(6,*) 'Message from alloc:'
              write(6,*) 'These density routines are not operating yet.'
              write(6,*) 'Change DENSI in input or ask the experts.'
              Call Abend2( 'LUCIAREL is quitting.' )
            end if
          else if (IDENSI.eq.3.or.IDENSI.eq.4) then
* We will calculate UU and UB, real and imaginary parts
            LRHO1 = NTOOB ** 2
            do IROOTCT = 1,NROOT,1
              do IRI = 1,2,1
                do JBU = 1,2,1
                  do IBU = 1,2,1
* Redundant density blocks!
c                   IBU = 2
                    call memman(KRHO1_DBG(IBU,JBU,IRI,IROOTCT),
     &                   LRHO1,'ADDS  ',2,'RHO_DB')
                  end do
                end do
              end do
            end do
            if (IDENSI.eq.4) then
* Still to be allocated !!
              write(6,*) 'Spin-dependent two-particle densities '
              write(6,*) 'not yet debugged. '
              Call Abend2( 'I would better quit.' )
            end if
          end if
*
* We set the pointer for the T array here as well. There is no
* explicit access to the array in this case, so no memory faults
* should occur.
          call memman(KT_CC,1,'ADDS  ',2,'KT_CC ')
        else if (SIGDEN_ROUTE.eq.'SIGDEN') then
*
* Allocate the amplitude array T only with general sigden
* structure for integrals or densities.
          NINT_CC = 0
          do IL1=1,NL1D,1
            NINT_CC = NINT_CC + 2 * NTOOB**2
          end do
*
          do IL2=1,NL2D,1
            NINT_CC = NINT_CC + 2 * L2D(IL2)
          end do
*
          call memman(KT_CC,NINT_CC,'ADDS  ',2,'KT_CC ')
C FIXME: Add symmetry savings
          write(6,*) 'NOTE:'
          write(6,*) 'Allocation should be modified to account '
          write(6,*) 'for possible symmetry savings.'
        end if
      end if
*
*.4. Integrals (ij!kk) -(ik!kj), pointer and space
      CALL MEMMAN(KPNIJ,NTOOB ** 2,'ADDS  ',1,'KPNIJ ')
      CALL MEMMAN(KIJKK,NTOOB**2 *(NTOOB+1) / 2, 'ADDS  '
     &            ,2,'KIJKK ')
CTF  Allocate memory for initial guess for CI vector
      IF(MXP1+MXP2+MXQ. NE. 0 .OR. IPROCC.NE.0 ) THEN
*. Allocate memory for explicit hamiltonian and roots
        NSBDET = MXP1 + MXP2 + MXQ
        NSBDETP = MAX(NSBDET,IPROCC)
        MXP = MXP1 + MXP2
        CALL MEMMAN(KSBEVC,NSBDETP*NROOT,'ADDS  ',2,'KSBEVC')
        CALL MEMMAN(KSBEVL,       NROOT,'ADDS  ',2,'KSBEVL')
*. KSBIDT must be able to hold list of subspace dets / print dets
        LSCR = NSBDETP+1
        CALL MEMMAN(KSBIDT,LSCR,'ADDS  ',1,'KSBIDT')
*. ( NSBDET + 1 elements allocated since NSBDET+1 values are
* obtained in order to check for degenerencies )
CADDB
        CALL MEMMAN(KSBCNF,LSCR        ,'ADDS  ',1,'KSBCNF')
CADDE
        CALL MEMMAN(KSBIA ,LSCR        ,'ADDS  ',1,'KSBIA ')
        CALL MEMMAN(KSBIB ,LSCR        ,'ADDS  ',1,'KSBIB ')
*
        LH0 = MXP*(MXP+1)/2 + MXP1*MXQ + MXQ
        CALL MEMMAN(KH0,LH0,'ADDS  ',2,'KH0   ')
*. Scratch space for manipulating with H0
        LH0SCR =
     &  MAX(6*NSBDETP,4*NSBDETP+4*NOCOB,
     &      MXP1*(MXP1+1)/2+MXP1**2)
        CALL MEMMAN(KH0SCR,LH0SCR      ,'ADDS  ',2,'KH0SCR')
      ELSE
        KSBEVC = 1
        KSBEVL = 1
        KSBIDT = 1
        KSBCNF = 1
        KH0    = 1
        KH0SCR = 1
      END IF
*. indices for pair of orbitals symmetry ordered
*
      CALL MEMMAN(KINH1,NTOOB*NTOOB,'ADDS  ',2,'KINH1  ')
*
      IF (NOMOFL.EQ.0) THEN
        write(6,*)
        WRITE(6,*)
     &  ' Alloc : General allocation for SPAO trafo '
        WRITE(6,*)
     &  '         still to be done !! (now locally)'
*. MO-AO transformation
CTF     CALL MEMMAN(KMOAOIN,NTOOB ** 2,'ADDS  ',2,'MOAOIN')
CTF     CALL MEMMAN(KMOAOUT,NTOOB ** 2,'ADDS  ',2,'MOAOUT')
*. Spinor-spin atomic orbital transformation
CTF     CALL MEMMAN(KSPSAOIN,4*NTOOB**2,'ADDS  ',2,'SPOAIN')
CTF     CALL MEMMAN(KSPSAOUT,4*NTOOB**2,'ADDS  ',2,'SPOAUT')
*. and their offsets
CTF     CALL MEMMAN(KISPSAOIN,4,'ADDS  ',1,'ISPOAI')
CTF     CALL MEMMAN(KISPSAOUT,4,'ADDS  ',1,'ISPOAU')
*
      ELSE
        KMOAOIN = 0
        KMOAOUT = 0
        KSPSAOIN = 0
        KSPSAOUT = 0
        KISPSAOIN = 0
        KISPSAOUT = 0
      END IF
*. (more to come )
      CALL MEMMAN(IDUMMY,IDUMMU,'CHECK',IDUMMY,'Dummy ')
      RETURN
      END
