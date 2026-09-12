!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

!#define DEBUG_SOC

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck onegen */
      SUBROUTINE ONEGEN(WORK,LWORK)
C*****************************************************************************
C
C     This routine goes through the list of one-electron integrals specified
C     in the COMMON block XCBCLS.
C
C     Call other routines for integrals generation.
C
C     Called from: MINPRP
C                  PAMSET
C                  LSQMAT
C
C     Written by T.Saue - May 1996
C
C*****************************************************************************
! checking
      use interface_ao_specific
      use fde_mod

      use mxcent_mod
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
      DIMENSION WORK(LWORK)
      CHARACTER COMB*4
      LOGICAL TEST
      type(fde_import) :: itmp
C
#include "dcbham.h"
#include "dcbcls.h"
#include "cbihr1.h"
#include "cbiher.h"
#include "nuclei.h"

      CHARACTER CPUTID*12, WALLTID*12, SECTID*12

      CALL QENTER('ONEGEN')
C.....With modified Hamiltonian this routine is skipped for now
      IF(MDIRAC) GOTO 20
      !mi: always print nuclear contribution to the dipole moment (costs nothing)
      TEST=(DSUM(NUCIND,CHARGE,1).GT.D0)
      IF(TEST) CALL DIPNUC(1,.FALSE.)
C     Dump AO information in case of numerical integration      
      call interface_ao_write()
C
C     Loop through list of class labels; calculate only the active ones(CLSCAL)
C     =========================================================================
C
      IF (NOSMLV) CALL NOSSNUCATT
      IF (IPRONE .GT. 2) CALL GETTIM(CPU0,WALL0)
      DO 10 I = 1,NPRPCLS
      IF(.NOT.CLSCAL(I)) GOTO 10
        IF (IPRONE.GE.2) WRITE(LUPRI,'(/2A,2X,A)')
     &       '* ONEGEN: class operator and class comb ',
     &       CLSINT(I), CLSCMB(I)
Chj     CLSCMB and PDOINT is now '+', '-', or '0',
Chj     corresponding to a factor of +1, -1, or 0 on that block.
C       The four blocks (AO basis): LL, SL, LS, SS
        DO J = 1,4
          IF (CLSCMB(I)(J:J) .EQ. '0') THEN
             COMB(J:J) = 'F'
          ELSE
             COMB(J:J) = 'T'
          END IF
        END DO
C
C     Nuclear attraction integrals (molecular field)
C
      IF    (CLSINT(I).EQ.'MOLFLD ') THEN
        IF(NOSMLV) COMB(4:4) = 'F'
        CALL NUCATT(WORK,LWORK,COMB,IPRONE)
C
C     Effective electronic density on nuclei
C
      ELSEIF (CLSINT(I).EQ.'EFFDEN ') THEN
        CALL EFFECTIVE_DENS_NUC(WORK,LWORK,COMB)
C
C     Self-energy(3component)
C
      ELSEIF (CLSINT(I)(3:6).EQ.'FGMG') THEN
        CALL NUMINT(I,COMB,3,'AN',IPRONE) 
C
C     Self-energy(Flambaum)
C
      ELSEIF (CLSINT(I)(1:4).EQ.'SEFG') THEN
        IF(NOSMLV) COMB(4:4) = 'F'
        CALL NUMINT(I,COMB,0,'SY',IPRONE)       
C
C     Self-energy(Pyyko)
C
      ELSEIF (CLSINT(I)(1:4).EQ.'SEPZ') THEN
        IF(NOSMLV) COMB(4:4) = 'F'
        CALL NUMINT(I,COMB,0,'SY',IPRONE)       
C
C     Uehling potential for vacuum polarization
C
      ELSEIF (CLSINT(I).EQ.'UEHLING') THEN
C
        IF(NOSMLV) COMB(4:4) = 'F'
        CALL NUMINT(I,COMB,0,'SY',IPRONE)       
C
C
C     Contribution from (static) embedding potential
C
      ELSEIF (CLSINT(i).EQ.'FDEVEMB') THEN
        IF(NOSMLV) COMB(4:4) = 'F'
        call fde_get_import_info(itmp)
        if (itmp%im_vemb) then
            CALL fde_calculate_emb_pot_mat(COMB,IPRONE)
         endif
C
C     Contribution from beta-matrix
C
      ELSEIF(CLSINT(I).EQ.'BETAMAT') THEN
        CALL BETAMT(COMB,WORK,LWORK,IPRONE)
C
C     Contribution from CAP integrals
C
      ELSEIF(CLSINT(I).EQ.'CAP_RE '.OR.CLSINT(I).EQ.'CAP_IM ' .OR.
     &       CLSINT(I).EQ.'CAP_OVL'.OR.CLSINT(I).EQ.'CAPD1R ') THEN
CMI     ... this serves only for 'catching' CAP... labels because
CMI     these integrals are not computed in the abacus part but are from external source
       CONTINUE

CMI  Contribution from Heff_EDM integrals
!     ELSEIF (CLSINT(I).EQ.'EDM') THEN
!       CALL EDM_EFF_OPER(WORK,LWORK,IPRONE)

      ELSE
C
C Process all other integrals, which are coded in the abacus part !
C
        TEST = ALLATM.AND.
     &         ((CLSINT(I).EQ.'NEFIELD').OR.
     &          (CLSINT(I).EQ.'ELFGRDC').OR.
     &          (CLSINT(I).EQ.'RM1N1H ').OR.
     &          (CLSINT(I).EQ.'NUCSNLO').OR.
     &          (CLSINT(I).EQ.'FERMI C').OR.
     &          (CLSINT(I).EQ.'PVC    ').OR.
     &          (CLSINT(I).EQ.'NSTCGO '))
        IF(TEST) THEN
          ALLATM = .FALSE.
          NPATOM = 1
          DO K = 1,NUCIND
            IPATOM(1)=K
            CALL OP1GEN(CLSINT(I),COMB,IORDCL(I),WORK,LWORK,IPRONE)
          ENDDO
          ALLATM = .TRUE.
          IPATOM(1) = 0
          NPATOM    = 0
        ELSE
          CALL OP1GEN(CLSINT(I),COMB,IORDCL(I),WORK,LWORK,IPRONE)
        ENDIF
      ENDIF
      CLSCAL(I) = .FALSE.
   10 CONTINUE
C
      IF (IPRONE .GT. 2) THEN
         CALL GETTIM(CPUE,WALLE)
         CPUTID  = SECTID(CPUE-CPU0)
         WALLTID = SECTID(WALLE-WALL0)
         WRITE(LUPRI,'(/4A)')
     &  'Total time used in ONEGEN (CPU) ',CPUTID,' and (WALL) ',WALLTID
      END IF
 20   CONTINUE
      CALL QEXIT('ONEGEN')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck numint_ang */
      SUBROUTINE NUMINT(IOP,COMB,NXYZ,MATINF,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      
C On iuput: IOP   ... Index of operator: e.g. CLSINT(IOP)='UEHLING'
C           IXYZ  ... Index of XYZ, which is used in angular roop
C                     ( WT8 = WT8 * ANG(IY,IXYZ) )
C           NXYZ  ... 0 for scalar operator, 3 for vector operator
C           MATINF... SY : Symmetric
C                     AN : Anti-symmetrix
C                     OF : Off-diagonal
C      
C The numerical integration for constructing AO basis matrix
C
C Called from: ONEGEN
C
C A. Sunaga Nov. 2018
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use num_grid_cfg
      use interface_ao
      use dirac_ao_eval
      use dirac_ao_eval_sub
      
      use primitives_mod
      use maxaqn_mod
      use mxcent_mod
      use symmet_mod
      use nucdata_mod
#include "implicit.h" 
#include "priunit.h" 
#include "pi.h"
C
#include "nuclei.h"      
#include "ibtfun.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbgen.h"
#include "dcbcls.h"
C      
      LOGICAL ADDSML
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4,MATINF*2
      CHARACTER RADC*2,ANGC*2,RADANGC*6
      PARAMETER (D0=0.0D0,D1=1.0D0,D2=2.0D0,DM1=-1.0D0)      
      real(8), allocatable :: ANG(:,:),GAO(:,:),GAO2(:,:),BUFGAO(:,:)
      real(8), allocatable :: BUF(:),TC(:), NUMIR(:),R(:)
      real(8), allocatable :: INTAO(:),INTAOTRI(:)
      real(8) C(3,MXCENT)
C
      CHARACTER CPUTID*12, WALLTID*12, SECTID*12, TITLE_PRINT*50      
C
      CALL QENTER('NUMINT')
      CALL GETTIM(CPU0,WALL0)
      IF(num_grid_cfg_numint_print)THEN
         CALL NUMINTOP !Welcome message
      ENDIF      
C      
      NCUTZ = num_grid_cfg_minz
      INTANG=num_grid_cfg_ANGMIN
      m = NLEB(INTANG)
      NDIM2 = NTBAS(0)*NTBAS(0)
      allocate( ANG      (m,4)   )
      allocate( TC       (NDIM2) )
      allocate( INTAO    (NDIM2) )
      allocate( INTAOTRI (NDIM2) )
C      
      CALL DZERO(TC,NDIM2)
C      
      CALL LEBEDEV(INTANG,ANG(1,1),ANG(1,2),ANG(1,3),ANG(1,4),n)
      CALL DSCAL(m,4.0D0*PI,ANG(1,4),1)
! DSCAl is corresponding to the below part in SLEB/dftgrf.F
!      DO J=1,NCOUNT
!         WT(J)=WT(J)*4.0D0*PI
!     ENDDO
C
      IF(NXYZ.EQ.0)THEN
        INIXYZ = 0
        MAXXYZ = 0
      ELSEIF(NXYZ.EQ.3)THEN
        INIXYZ = 1
        MAXXYZ = 3
      ELSEIF(NXYZ.LT.0)THEN
        INIXYZ = -NXYZ
        MAXXYZ = -NXYZ
      ENDIF
C      
      DO 5 IXYZ = INIXYZ, MAXXYZ
C
C     Label for the operator
C         
      IF(IXYZ.EQ.0)THEN
!     This transformation is important, because LABEL must be 8, while CLSINT(IOP) is 7
        LABEL = CLSINT(IOP)
        WRITE(RTNLBL(2)(3:4),'(I2)') 1
      ELSE   
        LABEL(2:8) = CLSINT(IOP)
        WRITE(RTNLBL(2)(3:4),'(I2)') ISYMAX(IXYZ,1)+1
        IF(IXYZ.EQ.1)THEN
          LABEL(1:1) = 'X'
        ELSEIF(IXYZ.EQ.2)THEN
          LABEL(1:1) = 'Y'
        ELSEIF(IXYZ.EQ.3)THEN
          LABEL(1:1) = 'Z'
        ELSE
          CALL QUIT('NUMINT_ANG: The initial character of LABEL' //
     &               'should be X/Y/Z')
        ENDIF
      ENDIF 
C         
      CALL DZERO(INTAO   ,NDIM2)
      CALL DZERO(INTAOTRI,NDIM2)
      NAT = 0
      ISYMOP = 0
      DO 10 IATOM = 1, NUCIND !change
        IF(.NOT.NOORBT(IATOM)) THEN 
         IF(NINT(CHARGE(IATOM)).GE.NCUTZ .OR.
     &      num_grid_cfg_allatom.EQV. .TRUE.) THEN
         NAT = NAT + 1
         GEXP = GNUEXP(IATOM)
         NDEG = MULT(ISTBNU(IATOM))
         C(1,IATOM)= PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,IATOM)
         C(2,IATOM)= PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,IATOM)
         C(3,IATOM)= PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,IATOM)          
C         
C   SET DEFAULT VALUES OF PARAMETERS FOR CALCULATION
C
         MXNP = 4000            !same to GRASP
         CALL MEMSET_NUM(CHARGE(IATOM),IATOM,MXNP,RNT,H,EPH,GRDC,
     &                    num_grid_cfg_grasp,IPRINT)
C
         allocate( NUMIR(MXNP) )
         allocate( R(MXNP)     )
C         
         IF(num_grid_cfg_grasp)THEN
           R(1) = RNT
           DO J = 2,MXNP
              R(J) = R(J-1)*EPH
           ENDDO
C     Print out the information
           IF(num_grid_cfg_numint_print)THEN
             CALL NUMINFO(ISYMOP,IATOM,CHARGE(IATOM),RNT,R(MXNP),H,MXNP)
           ENDIF           
C     R**2 is multiplized in NUMPOT
           CALL NUMPOT(IOP,NUMIR,R,CHARGE(IATOM),GEXP,MXNP,IPRINT)
C Additional R is needed for GRASP grid (QUAD2)
           DO I = 1,MXNP
              NUMIR(I) = NUMIR(I) * R(I)
           ENDDO
C
         ELSE  ! num_grid_cfg_grasp: DIRAC grid start from here
           R(1) = RNT
           DO IR = 2,MXNP
             R(IR) = (R(IR-1)+GRDC)*EPH-GRDC
           ENDDO
C     Print out the information
           IF(num_grid_cfg_numint_print)THEN
             CALL NUMINFO(ISYMOP,IATOM,CHARGE(IATOM),RNT,R(MXNP),H,MXNP)
           ENDIF           
C R**2 is multiplized in NUMPOT
           CALL NUMPOT(IOP,NUMIR,R,CHARGE(IATOM),GEXP,MXNP,IPRINT)
           W = (RNT+GRDC)*H 
           NUMIR(1) = NUMIR(1) * W
           DO IR = 2,MXNP
             W = (R(IR)+GRDC)*H
             NUMIR(IR) = NUMIR(IR) * W
           ENDDO
C          
         ENDIF   !num_grid_cfg_grasp
C     Numerical integration
C     Here, INTAOTRI is a buffer
         CALL INTSUM(INTAOTRI,NDIM2,NUMIR,R,MXNP,ANG,m,IXYZ,
     &    C(1,IATOM),num_grid_cfg_grasp,MATINF,IATOM,LABEL,IPRINT)
C     writing the matrix in output
         IF(IPRINT.GE.12)THEN
           WRITE(TITLE_PRINT,'(5A)')
     &       'NUMINT: ',LABEL,'Symmetrized AO-matrix of ',
     &       NAMN(IATOM),"atom"
           CALL HEADER(TITLE_PRINT,-1)
           CALL OUTPUT(INTAOTRI,1,NTBAS(0),1,NTBAS(0),
     &                 NTBAS(0),NTBAS(0),-1,LUPRI) 
         ENDIF   

C       Scale weights with multiplicity of center
         IF(NDEG.GT.0) THEN
            FAC = dble(NDEG) 
            CALL DSCAL(NDIM2,FAC,INTAOTRI,1)
         ENDIF
         
C     Here, INTAOTRI is a buffer
         CALL DAXPY(NDIM2,D1,INTAOTRI,1,INTAO,1)
         CALL DZERO(INTAOTRI,NDIM2)
         deallocate(NUMIR)
         deallocate(R)
         ENDIF    ! IF(ZAN(NAT).GE.NCUTZ .OR.num_grid_cfg_allatom.EQV. .TRUE.)
        ELSE      ! Ghost center case
         IF(num_grid_cfg_numint_print)THEN
            WRITE(LUPRI,'(A,A)')" - For ", NAMN(IATOM)
            WRITE(LUPRI,'(A,/)')"     effQED potential is skipped"
         ENDIF  
        ENDIF   ! IF(NOORBT(IATOM))  
        IF(IATOM.EQ.NUCIND)THEN
          num_grid_cfg_numint_print = .FALSE.
        ENDIF
10    CONTINUE   ! IATOM = 1, NUCIND, ISYMOP = 0, MXOPR 
C
      IF(IPRINT.GE.6)THEN
        WRITE(TITLE_PRINT,'(5A)')
     &    'NUMINT: ',LABEL,'Total AO-matrix'
        CALL HEADER(TITLE_PRINT,-1)
        CALL OUTPUT(INTAO,1,NTBAS(0),1,NTBAS(0),
     &              NTBAS(0),NTBAS(0),-1,LUPRI) 
      ENDIF  
C
C    Transformation from full-matrix to triangle-matrix
C      
      CALL MTFULTRI(INTAO,INTAOTRI,NTBAS(0),IPRINT)
C     
C     Generate integral labels
C
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = MATINF
      IF(MATINF.EQ.'OF') RTNLBL(2)(1:2) = 'SY' !OFF diagonal but non-anti
      RTNLBL(2)(5:8) = COMB
C
C     Write integrals to file
C
      CALL WRTPRO(INTAOTRI,NNBBASX,LABEL,RTNLBL,IPRINT)
C      
C     debug print
C      
      NRAD = int(-log10(num_grid_cfg_RADINT))
 5    CONTINUE !IXYZ
      deallocate(ANG)
      deallocate(INTAO)
      deallocate(INTAOTRI)
C
      CALL GETTIM(CPUE,WALLE)
      CPUTID  = SECTID(CPUE-CPU0)
      WALLTID = SECTID(WALLE-WALL0)
      WRITE(LUPRI,'(6A)')
     &     '  Total time used in ',CLSINT(IOP),' (CPU) ',CPUTID,
     &     ' and (WALL) ',WALLTID    
      CALL QEXIT('NUMINT')
      RETURN         
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck NUMPOT */
      SUBROUTINE NUMPOT(IOP,NUMIR,R,Z,GEXP,MXNP,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Getting numerical potential for operators
C
C     On input:
C     R(MXNP)... Radial grid
C     Z......... Nuclear charge       
C     GEXP...... 3/(2 * <R^2_nuc>)      
C      
C     On output:
C     NUMIR ... Numerical potential
C
C     Called from: NUMINT
C
C     A. Sunaga AuG. 2019      
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use mxcent_mod
#include "implicit.h" 
#include "priunit.h"
#include "nuclei.h"
#include "dcbcls.h"
C
      REAL*8 NUMIR(MXNP), R(MXNP)
C     
      IF(GAUNUC)THEN
         NUCTYP = 2
      ELSE   
         NUCTYP = 1
      ENDIF   
C
      IF(CLSINT(IOP).EQ.'SEFGMG')THEN
        CALL INTMG(NUMIR,R,Z,GEXP,MXNP,NUCTYP)
        DO IR = 1,MXNP
         NUMIR(IR) = NUMIR(IR) * R(IR) * R(IR) 
        ENDDO
      ELSEIF(CLSINT(IOP).EQ.'UEHLING')THEN
        CALL VCMPOL(NUMIR,R,Z,GEXP,MXNP,NNUC,NUCTYP,1)
        DO IR = 1,MXNP
          NUMIR(IR) = NUMIR(IR) * R(IR)
        ENDDO
      ELSEIF(CLSINT(IOP).EQ.'SEFGHF')THEN
        CALL INTHF (NUMIR,R,Z,GEXP,MXNP,NUCTYP)
        DO IR = 1,MXNP
          NUMIR(IR) = NUMIR(IR) * R(IR) * R(IR) 
        ENDDO
      ELSEIF(CLSINT(IOP).EQ.'SEFGLF')THEN
        CALL INTLF (NUMIR,R,Z,MXNP             )
        DO IR = 1,MXNP
          NUMIR(IR) = NUMIR(IR) * R(IR) * R(IR) 
        ENDDO
      ELSEIF(CLSINT(IOP).EQ.'SEPZ'  )THEN
        CALL INTPZ (NUMIR,R,Z,MXNP             )
        DO IR = 1,MXNP
          NUMIR(IR) = NUMIR(IR) * R(IR) * R(IR) 
        ENDDO
      ENDIF
C
      IF(IPRINT.GE.2)THEN
        WRITE(LUPRI,'(2X,/A)')
     &        '**Numerical potential (not operator) at each radial grid'
        WRITE(LUPRI,'(2A)')"  Potential: ",CLSINT(IOP)
        write(LUPRI,'(A)')'   R(Bohr),           Potential(a.u.)'
        DO IR = 1,MXNP
          write(LUPRI,'(A,E14.8e2,E20.8e2)')
     &      "  ",R(IR),-NUMIR(IR)/(R(IR)*R(IR))
        ENDDO
      ENDIF                  
      RETURN
      END
      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck intsum */
      SUBROUTINE INTSUM(INTSO,NDIM2,VWRR,R,MXNP,ANG,m,IXYZ,C,
     &                   GRASP,MAT,IATOM,LABEL,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     The integration of radial and angular integral   
C      
C     On input:
C     NDIM2     ... NDIM2 = NTBAS(0)*NTBAS(0)
C     VWRR      ... r-dependent part (Potential and weight)
C                     DIRAC grid: VWRR = V(r) * W(r) * r * r 
C                     GRASP grid: VWRR = V(r) * r    * r * r 
C                   V(r) ... Target potential
C                   W(r) ... Weight for radial integration using Lindh's grid
C     R         ... Radial grid      
C     MXNP      ... Number of radial grid
C     ANG       ... Angular grid (Lebedev grid)
C                   Please multiply 4pi to ANG(1,4) before this routine.
C                   ANG(iy,IXYZ): iy is the grid,
C                   IXYZ=1,3 is XYZ, IXYZ=4 is weight
C     m         ... Number of angular grid
C                    e.g. m = NLEB(num_grid_cfg_ANGMIN)
C     IXYZ      ... Additional angular momentum
C     C         ... Center of IATOM: 
C                   It is obtained by CALL SET_GRID_INFO(NAT,NUM,ZAN,C)
C     GRASP     ... TRUE : GRASP grid
C                   FALSE: DIRAC grid
C     MAT       ... SY : Symmetric
C                   AN : Anti-symmetrix
C                   OF : Off-diagonal
C
C     On output:
C     INTSO ... AO matrix at SO basis after the integration + previous INTSO     
C
C The numerical integration for constructing AO basis matrix
C
C Called from: NUMINT_ANG, NUMINT
C
C A. Sunaga Jul. 2019
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use num_grid_cfg
      use interface_ao
      use dirac_ao_eval
      use dirac_ao_eval_sub
      use xc_max_block_length

      use primitives_mod
      use shells_mod
      use maxaqn_mod
      use mxcent_mod
      use symmet_mod
#include "implicit.h" 
#include "priunit.h" 
#include "dcbbas.h"
C
#include "dgroup.h"
#include "dcbham.h"
C      
#include "nuclei.h"
      REAL*8 R(MXNP), ANG(m,4), INTSO(NDIM2), VWRR(MXNP), C(3,1)
C
      CHARACTER MAT*2, LABEL*8
      CHARACTER TITLE_PRINT*100
      LOGICAL GRASP,TRIMAT,ANTI
      PARAMETER (D0=0.0D0,D1=1.0D0,DM1=-1.0D0,D2=2.0D0)      
      REAL*8    ALPHA(2)
      INTEGER*8 IA(2), IB(2), ITC(2), IOFA(2), IOFB(2)
      real(8), allocatable :: GAO(:,:),GAO2(:,:),GMAX(:,:),BUFGAO(:,:)
      real(8), allocatable :: TC(:),TD(:,:),INTAOBUF(:)  
      real(8), allocatable :: ANGXYZ(:),XYZ8(:,:)
      real(8), allocatable :: pxyz(:,:),p2(:)
      real(8), parameter   :: bitstring_parity(0:7) =
     &    (/1.0d0, -1.0d0, -1.0d0, 1.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0/)     
C    
      allocate( TC       (NDIM2) )
      allocate( INTAOBUF (NDIM2) )
      allocate( ANGXYZ   (m)     )
C
      SCREEN = num_grid_cfg_screen
      nblock_size = min(max_block_length, m)
      nsize3      = nblock_size*3
      allocate( XYZ8(nblock_size,3) )
      allocate( pxyz(nblock_size,3) )
      allocate( p2  (nblock_size  ) )
C
      call dirac_ao_eval_init(0,0,.FALSE.)       
      allocate( GMAX(0:nr_boson_ireps-1, nr_shells) )
C
      CALL DZERO(INTAOBUF,NDIM2)
      IF(IXYZ.EQ.0)THEN
        CALL DCOPY(m,ANG(1,4),1,ANGXYZ,1)
        MTXYZ = 0 
      ELSE  
        IF (IXYZ.LE.3)THEN
          DO iy = 1, m
            ANGXYZ(iy) = ANG(iy,4)*ANG(iy,IXYZ)
          ENDDO
          MTXYZ = ISYMAX(IXYZ,2)
       ELSE
         write(*,*)'IXYZ=',IXYZ
         write(*,*)'is wrong number!' !debug
         stop
        ENDIF
      ENDIF
C
      IF(MAT.EQ.'SY')THEN
         IA(1) = 1
         IB(1) = 1
         IA(2) = 2
         IB(2) = 2
         ITC(1) = 1
         ITC(2) = NTBAS(1)*(NTBAS(0)+1) + 1
         IOFA(1) = 1
         IOFB(1) = 1
         IOFA(2) = 1 + NTBAS(1)
         IOFB(2) = 1 + NTBAS(1)
         ALPHA(1) = D1
         ALPHA(2) = D1
      ELSE !ANTI or OFFD
         IA(1) = 2
         IB(1) = 1
         IA(2) = 1
         IB(2) = 2
         ITC(1) = NTBAS(1)          + 1 
         ITC(2) = NTBAS(1)*NTBAS(0) + 1
         IOFA(1) = 1 + NTBAS(1)
         IOFB(1) = 1
         IOFA(2) = 1
         IOFB(2) = 1 + NTBAS(1)
         ALPHA(1) = D1
         ALPHA(2) = D1
         IF(MC.EQ.1)THEN
           IA(1) =  1 
           ITC(1) =  1 
           IOFA(1) = 1
         ENDIF
      ENDIF
C
Ciaa  The next line must be removed if the magnetic part of the Ginges-Flambaum self-energy effective Hamiltonian
C     is used with i.beta.alpha dot SEFGMG, instead of the current ALPHADOT version (currently we use: i.alpha dot SEFGMG).
C     If next line is removed, then ALPHADOT must be replaced by BEALDOT in the inputs.
      IF(MAT.EQ.'AN')ALPHA(1) = DM1
C
C      
      nr_points_on_this_proc=MXNP*m
C
C     Determination of the maximum value for each (isym,ishell)
C     for cutoff of the numerical integration. The flow is
C      
C     1. Several AOs are picked up in the numerical integration region
C       do IR = 1, MXNP_ODD, (MXNP_ODD-1)/2
C     2. Getting the GMAX(isym,ishell), which is the largest AO at (isym,ishell)
C      
C     3. If the GMAX is smaller than the SCREEN, the get AO step is skipped.
C        at all radial grid.      
C         IF(GMAX(isym,ishell).GE.SCREEN)THEN
C      
      CALL DZERO(GMAX, nr_boson_ireps*nr_shells)
      MXNP_ODD = MXNP - mod(MXNP,2)
      do IR = 1, MXNP_ODD, (MXNP_ODD-1)/2  
        IOF_block = 1
C                 
C     In the case of m>100, a special treatment is needed.
C
        m_NUM_MAX = INT(m/100) + 1               
C
C    GET AO's (from dirac_ao_eval.F90, but adapted for this routine)
C
        IOFM = 1
        DO m_NUM = 1, m_NUM_MAX
          nblock_length = min(max_block_length, m - IOF_block+ 1)
          nblock_length = min(nblock_length, m)
          allocate( GAO2(nblock_length, NTBAS(0)) )
          CALL DZERO(GAO2,nblock_length*NTBAS(0))
C         The spherical grid (R(IR), ANG(IOFM,JXYZ)) is transfromed to
C         Cartesian grid (pxyz(nblock_size,3))
C         to use get_cartesian_ao_g0 routine
          DO JXYZ = 1,3
             XYZ8(:,JXYZ) = C(JXYZ,1)
             CALL DAXPY(nblock_length,R(IR),ANG(IOFM,JXYZ),1,
     &                               XYZ8(   1,JXYZ),1)
          ENDDO
C
          i = 1
          do ishell = 1, nr_shells
            do isym = 0, nr_boson_ireps - 1
              if (iand(isym, shell_stabilizing_irep(ishell)) == 0)then
               DO JXYZ = 1, 3
                pxyz(:,JXYZ) = DM1 
     &             * bitstring_parity(iand(irep_of_axes(JXYZ), isym))
     &             * center_xyz(shell_center(ishell), JXYZ)
               ENDDO
C               
               CALL DAXPY(nsize3,D1,XYZ8,1,pxyz,1)
C
               do k = 1, nblock_length
                  p2(k) = DDOT(3,pxyz(k,1),nblock_size,
     &                           pxyz(k,1),nblock_size) 
               end do
C
               call get_cartesian_ao_g0(nblock_length, 
     &                                    ishell,       
     &                                    pxyz(1,1),           
     &                                    pxyz(1,2),           
     &                                    pxyz(1,3),           
     &                                    p2,
     &                                    GAO2(1, i))
C              
               DO k = 1, nblock_length
                 DO ICOMP = 1,KHKT(ishell)
                 GMAX(isym,ishell)
     &          = MAX( ABS(GAO2(k,i+ICOMP-1)), ABS(GMAX(isym,ishell)) )
                 ENDDO
               ENDDO
               i = i + cartesian_deg(ishell)
              end if
            end do  ! isym
          end do               ! ishell
          IOFM = IOFM + nblock_length
          IOF_block = IOF_block + nblock_length
          deallocate(GAO2)
        ENDDO                  ! m_NUM
      end do                ! IR         
C
      IF(GRASP)THEN
        allocate( TD(NDIM2,MXNP) )         
        do IR = 1, MXNP
          CALL DZERO(TC,      NDIM2)
          IOF_block = 1
          m_NUM_MAX = INT(m/100) + 1              
          allocate( GAO   (m, NTBAS(0)) )
          allocate( BUFGAO(m, NTBAS(0)) )
          CALL DZERO(GAO, m*NTBAS(0))
          CALL DZERO(BUFGAO,m*NTBAS(0))
C
          call dirac_ao_eval_init(0,0,.FALSE.) 
C
C    GET AO's
C
          IOFM = 1
          DO m_NUM = 1, m_NUM_MAX
C                 
C     In the case of m>100, a speciAL treatment is needed.
C     m is the number of angular grid             
C
            nblock_length = min(max_block_length, m - IOF_block+ 1)
            nblock_length = min(nblock_length, m)
            nblock3 = nblock_length*3
C              
            allocate( GAO2(nblock_length, NTBAS(0)) )              
            CALL DZERO(GAO2,nblock_length*NTBAS(0))
            DO JXYZ = 1,3
               XYZ8(:,JXYZ) = C(JXYZ,1)
               CALL DAXPY(nblock_length,R(IR),ANG(IOFM,JXYZ),1,
     &                                 XYZ8(1   ,JXYZ),1)
            ENDDO   
             
            i = 1
            do ishell = 1, nr_shells
              do isym = 0, nr_boson_ireps - 1
               if (iand(isym, shell_stabilizing_irep(ishell)) == 0)then                 
                IF(GMAX(isym,ishell).GE.SCREEN)THEN   
                 DO JXYZ = 1, 3
                  pxyz(:,JXYZ) = DM1 
     &               * bitstring_parity(iand(irep_of_axes(JXYZ), isym))
     &               * center_xyz(shell_center(ishell), JXYZ)
                 ENDDO 
C
                 CALL DAXPY(nsize3,D1,XYZ8,1,pxyz,1)
C
                 do k = 1, nblock_length
                    p2(k) = DDOT(3,pxyz(k,1),nblock_size,
     &                             pxyz(k,1),nblock_size) 
                 end do
C      
                 call get_cartesian_ao_g0(nblock_length, 
     &                                      ishell,       
     &                                      pxyz(1,1),           
     &                                      pxyz(1,2),           
     &                                      pxyz(1,3),           
     &                                      p2,
     &                                      GAO2(1, i))
                ENDIF  ! GMAX cutoff  
                i = i + cartesian_deg(ishell)
               end if
              end do  ! isym
            end do               ! ishell
            DO J = 1, NTBAS(0)
               CALL DCOPY(nblock_length,GAO2(1,J),1,GAO(IOFM,J),1)
            ENDDO
            IOFM = IOFM + nblock_length
            deallocate(GAO2)
            IOF_block = IOF_block + nblock_length
          ENDDO                 ! m_NUM          
C          
          DO  iy = 1, m
            CALL DAXPY(NTBAS(0),ANGXYZ(iy),GAO(iy,1),m,BUFGAO(iy,1),m)
          ENDDO
          CALL DSCAL(NTBAS(0)*m,VWRR(IR),BUFGAO(1,1),1)
          DO IC = 1, MC
           CALL DGEMM('T','N',NTBAS(IA(IC)),NTBAS(IB(IC)),m,
     &                 ALPHA(IC),BUFGAO(1,IOFA(IC)),m,GAO(1,IOFB(IC)),
     &                 m,D0,TC(ITC(IC)),NTBAS(0))
          ENDDO 
          CALL DCOPY(NDIM2,TC,1,TD(1,IR),1)
          deallocate(GAO)
          deallocate(BUFGAO)
        end do                ! IR
          CALL QUAD2(INTAOBUF,TD,MXNP,NDIM2)
        deallocate(TD)
      ELSE  !DIRAC grid from here
        do IR = 1, MXNP  
          CALL DZERO(TC,      NDIM2)
          IOF_block = 1
C                 
C     In the case of m>100, a speciAL treatment is needed.
C
          m_NUM_MAX = INT(m/100) + 1               
          allocate( GAO    (m, NTBAS(0)) ) 
          allocate( BUFGAO (m, NTBAS(0)) )
          CALL DZERO(GAO, m*NTBAS(0))
          CALL DZERO(BUFGAO,m*NTBAS(0))
C
          call dirac_ao_eval_init(0,0,.FALSE.) 
C
C    GET AO's
C
          IOFM = 1
          DO m_NUM = 1, m_NUM_MAX
            nblock_length = min(max_block_length, m - IOF_block+ 1)
            nblock_length = min(nblock_length, m)
C     This allocation is needed, because the size of GAO2
C     in get_cartesian_ao_g0 must be nblock_length 
            allocate( GAO2(nblock_length, NTBAS(0)) )   
            CALL DZERO(GAO2,nblock_length*NTBAS(0))
            DO JXYZ = 1,3
               XYZ8(:,JXYZ) = C(JXYZ,1)
               CALL DAXPY(nblock_length,R(IR),ANG(IOFM,JXYZ),1,
     &                                 XYZ8(1   ,JXYZ),1)
            ENDDO
C
            i = 1
            do ishell = 1, nr_shells
              do isym = 0, nr_boson_ireps - 1
               if (iand(isym, shell_stabilizing_irep(ishell)) == 0)then                 
                IF(GMAX(isym,ishell).GE.SCREEN)THEN   
                 DO JXYZ = 1, 3
                  pxyz(:,JXYZ) = DM1 
     &               * bitstring_parity(iand(irep_of_axes(JXYZ), isym))
     &               * center_xyz(shell_center(ishell), JXYZ)
                 ENDDO 
C
                 CALL DAXPY(nsize3,D1,XYZ8,1,pxyz,1)
C
                 do k = 1, nblock_length
                    p2(k) = DDOT(3,pxyz(k,1),nblock_size,
     &                             pxyz(k,1),nblock_size) 
                 end do
C      
                 call get_cartesian_ao_g0(nblock_length, 
     &                                      ishell,       
     &                                      pxyz(1,1),           
     &                                      pxyz(1,2),           
     &                                      pxyz(1,3),           
     &                                      p2,
     &                                      GAO2(1, i))
                end if
                i = i + cartesian_deg(ishell)
               ENDIF  ! GMAX
              end do  ! isym
            end do               ! ishell
C     Copy a set of 100 AO (GAO2) to total AO array (GAO)            
            DO J = 1, NTBAS(0)
              CALL DCOPY(nblock_length,GAO2(1,J),1,GAO(IOFM,J),1)
            ENDDO
            IOFM = IOFM + nblock_length
            IOF_block = IOF_block + nblock_length
            deallocate(GAO2)
          ENDDO                  ! m_NUM
C     Multiplying the factor for angular part ANGXYZ. BUFGAO = ANGXYZ*GAO
          DO  iy = 1, m
            CALL DAXPY(NTBAS(0),ANGXYZ(iy),GAO(iy,1),m,BUFGAO(iy,1),m)
          ENDDO
          DO IC = 1, MC
           CALL DGEMM('T','N',NTBAS(IA(IC)),NTBAS(IB(IC)),m,
     &                 ALPHA(IC),BUFGAO(1,IOFA(IC)),m,GAO(1,IOFB(IC)),
     &                 m,D0,TC(ITC(IC)),NTBAS(0))
          ENDDO
C sum  weight*r^2*V(r)*Basis^2        
          CALL DAXPY(NDIM2,VWRR(IR),TC,1,INTAOBUF,1)
          deallocate(GAO)
          deallocate(BUFGAO)
        end do                ! IR
      ENDIF                     ! GRASP or DIRAC
      IF(IPRINT.GE.6)THEN
        WRITE(TITLE_PRINT,'(5A)')
     &     'INTSUM: ',LABEL,'Before-symmetrized AO-matrix of ',
     &     NAMN(IATOM),"atom"
           CALL HEADER(TITLE_PRINT,-1)
           CALL OUTPUT(INTAOBUF,1,NTBAS(0),1,NTBAS(0),
     &                 NTBAS(0),NTBAS(0),-1,LUPRI) 
      ENDIF   
C
      CALL MTAOSO(INTAOBUF,INTSO,NTBAS(0),MTXYZ,IPRINT)
C
      IF(IPRINT.GE.6)THEN
        WRITE(TITLE_PRINT,'(5A)')
     &     'INTSUM: ',LABEL,'After-symmetrized AO-matrix of ',
     &     NAMN(IATOM),"atom"
           CALL HEADER(TITLE_PRINT,-1)
           CALL OUTPUT(INTSO,1,NTBAS(0),1,NTBAS(0),
     &                 NTBAS(0),NTBAS(0),-1,LUPRI) 
      ENDIF   
C
      deallocate( TC      )
      deallocate( INTAOBUF)
      deallocate( XYZ8    )
      deallocate( pxyz    )
      deallocate( p2      )
      deallocate( GMAX    )
      deallocate( ANGXYZ  )
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtfultri */
      SUBROUTINE MTFULTRI(FUL,TRI,NDIM,IPRINT)
C***********************************************************************
C
C     Matrix transform:
C     full-matrix (FUL) to lower-triangle matrix (TRI)
C
C     Written by A.Sunaga - Aug 2019
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C#include "dcbbas.h"
      DIMENSION FUL(NDIM*NDIM), TRI(NDIM*NDIM)
C
      DO J = 1, NDIM
        IOFF = 0
        IOFT = (J-1)*J/2
        DO IT = 1, J
           TRI(IT+IOFT) = FUL(J+IOFF)
           IOFF = IOFF + NDIM
        ENDDO
      ENDDO
C      
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE QUAD2(RESULT,TE,NR,NDIM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Numerical radial integral
C
C Called from: NMINT2
C
C copied from /QUAD/grasp.F/GRASP
C
C     Modified to DIRAC A. Sunaga Nov. 2018       
C     Modification: e.g. TE is two dimension 
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
C
      PARAMETER( W7 = 1.0D0/90.0D0,
     &           W6 = -4.0D0*W7,  W5 =   7.0D0*W7,
     &           W4 = 30.0D0*W7,  W3 =  44.0D0*W7,
     &           W2 = 68.0D0*W7,  W1 = 112.0D0*W7)
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      PARAMETER(EPSI = 1.0D-35,H=1.0D-2) !-35
!      REAL(8) TE(NR)      
!      DIMENSION TE(NR)      
      DIMENSION TE(NDIM,NR), RESULT(NDIM)
C     
C   FIND THE STARTING POINT FOR THE INTEGRATION. THIS IS THE POINT
C   WHERE THE INTEGRAND TE FIRST BECOMES SIGNIFICANT, AND IS CHOSEN
C   TO PREVENT DIVIDING BY ZERO.
C
      DO 5  II = 1,NDIM
C
C     In the atomic case (GRASP), 'DO 10 J = 1,"NR"' is fine, because
C     normally TE(II,nuclear region) >> TE(II,NR).          
C     In the case of molecule, a situation TE(II,nuclear region) < TE(II,NR)
C     can occur. So, Ayaki changed from 'DO 10 J = 1,"NR"' to 'DO 10 J = 1,"NR-4"'
C     "4" is the minumum number for JJ+3 =< NR        
C         
      DO 10 J = 1,NR-4
        IF (ABS (TE(II,J)) .GT. EPSI) GO TO 20
   10 CONTINUE
      RESULT(II) = ZERO
      GO TO 5
C
C   CALCULATE THE INTEGRAL USING EXTENDED SIMPSONS RULE, NEGLECTING THE
C   CONTRIBUTION FROM THE STARTING POINT. CHOOSE THE STARTING POINT
C   R(JJ) TO GIVE AN ODD NUMBER OF INTEGRATION POINTS.
C   THE CRAY VERSION OF THE CODE USES THE VECTOR FUNCTION SSUM TO
C   PERFORM THE SUMMATION.
C
   20 SUM = ZERO
!      write(*,*)"TE1,J,NR",TE(J),J,NR !sunaga
      NP = NR-J+1
      ND2 = NP/2
      N1D2 = (NP+1)/2
      K = ND2-N1D2+1
      JJ = J+K
      JJ1 = JJ+1
      DO 30 I = JJ1,NR-3,2
         SUM = SUM + W1*TE(II,I) + W2*TE(II,I+1)
 30   CONTINUE
      SUM = SUM + W1*TE(II,NR-1)
      SUM = SUM + W4*( TE(II,JJ)   + TE(II,NR)   )
     &          + W5*( TE(II,JJ1)  + TE(II,NR-1) )
     &          + W6*( TE(II,JJ+2) + TE(II,NR-2) )
     &          + W7*( TE(II,JJ+3) + TE(II,NR-3) )
C      ENDIF
C
C   ADD THE CONTRIBUTION FROM 0 TO R(JJ), PLUS THE CONTRIBUTION
C   TO SIMPSONS RULE FROM TE(JJ).
C
      COP = TE(II,J+1)/TE(II,J)
      IF (COP .LE. ONE) THEN
         RPX = ZERO
      ELSE
         RPX = TE(II,JJ)/LOG (COP)
      ENDIF
      RESULT(II) = H*(SUM+RPX)
    5 CONTINUE
      RETURN
C
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck memset */
      SUBROUTINE MEMSET_NUM(Z,IATOM,MXNP,RNT,H,EPH,GRDC,GRASP,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C     i) Estimating memories, by calculating the uehling potential of
C     point charge
C
C     ii) Deciding some values for numerical integral
C
C     On iuput: Z...Nuclear charge of the atom
C      
C      
C     On output: GRDC - Grid for numerical integral
C                RNT...initial grid      
C                H.....step size
C                EPH = exp(H)
C      
C     
C     Called from: NUMINT
C
C    Modified for DIRAC, A. Sunaga Dec. 2018
C  
************************************************************************
      use num_grid_cfg

      use maxaqn_mod
      use mxcent_mod
      use ccom_mod, only: NHTYP
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "nuclei.h"
#include "dcbham.h"      
#include "dcbgen.h"
C     
      LOGICAL ADDSML,GRASP,UEHCUT
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0)     
      DIMENSION NUCORB(NHTYP,2,NUCIND),AA(2,NHTYP,2,NUCIND)
C
      UEHCUT = num_grid_cfg_uehcut
C      
      N = MXNP
      BZ = 0.074D0 + 0.35D0*Z/CVAL
      CVAL3 = CVAL*CVAL*CVAL
      PREFAC = BZ * Z*Z*Z*Z / CVAL3
CUEHLING      
      PI3 = THREE*PI
      WD = CVAL+CVAL
      FACTOR = -(Z+Z)/(PI3*CVAL)
      IF(GRASP)THEN 
C
C     SET UP POTENTIAL FOR POINT NUCLEUS
C
        H     = 1.0D-2
        EPH   = EXP(H)
        RNT   = 1.0D-06/Z       
        RTEMP = RNT
        DO 1 I = 1,N
           X = PREFAC*EXP(-Z*RTEMP)
           IF(X.LT.1.0D-50)THEN
              TEMP = X
              MXNP = I
              exit
           ENDIF   
           RTEMP = RTEMP*EPH
C
c     Ayaki
c     As far as I check GRASP, the number of grid (MXNP) must be odd number,
c     for the integral of Vacuum polarization using Ne and Hg atoms.         
c     If it is even number, the value of QUAD
c     does not agree with the one of GRASP (Maximu error in each grid is 10e-5).
C           
 1      CONTINUE
C
        IF(UEHCUT)THEN
          RTEMP = RNT
          DO J = 1,N
             X = WD*RTEMP
             IF(X.GT.1.0D+02)THEN !same to FUNK
                TEMP = FACTOR*FUNK (X,1)              
                MXNP = J
                exit
             ENDIF   
             RTEMP = RTEMP*EPH
          ENDDO
        ENDIF        
      ELSE !num_grid_cfg_grasp
C
C     Get information about basis sets
C
        CALL NUCBAS(NUCORB,AA,IPRINT)
        ADDSML = num_grid_cfg_force_4c_grid.OR.
     &              (MC.GT.1.AND.(.NOT.LEVYLE))
        RADDER = num_grid_cfg_RADINT
        CALL LMGGRID(H,RNT,RH,RADDER,IATOM,ADDSML,NUCORB,AA,IPRINT)
C
        EPH   = EXP(H)
        GRDC  = RNT/(EPH-ONE)
        RTEMP = RNT
C     
C     SET UP POTENTIAL FOR POINT NUCLEUS
C
        DO 2 I = 1,N
           X = PREFAC*EXP(-Z*RTEMP)
           IF(X.LT.1.0D-50)THEN
              TEMP = X
              MXNP = I 
             exit
           ENDIF 
           RTEMP = (RTEMP+GRDC)*EPH-GRDC
 2      CONTINUE
C
        IF(UEHCUT)THEN
          RTEMP = RNT
          DO J = 1,N
             X = WD*RTEMP
             IF(X.GT.1.0D+02)THEN !same to FUNK
                TEMP = FACTOR*FUNK (X,1)              
                MXNP = J
                exit
             ENDIF   
             RTEMP = (RTEMP+GRDC)*EPH-GRDC
          ENDDO
        ENDIF
      ENDIF    !num_grid_cfg_grasp
C
      IF (IPRINT.GE.6) THEN
         WRITE(LUPRI,'(2X,A)')'**Truncation of radial grid'
         WRITE(LUPRI,'(2X,A,E8.2e2,A,/)')
     &    'The value of outer grid ',TEMP,' is less than 1.0D-50'
      ENDIF   
C      
      RETURN
C
      END
************************************************************************
*
      SUBROUTINE VCMPOL(TB,R,Z,GEXP,N,NNUC,NUCTYP,NVPORD)
*      
*     NUCTYP:
*     1  Point
*     2  Gauss
*     3  Fermi
*     Ayaki: NUCTYP is optimized to the definition of the DIRAC code      
*
*     NVPORD:
*     1  second (Uehling)
*     2  Forth order   
*     
*   ----------------   SECTION 11   SUBPROGRAM 10  -----------------   *
*                                                                      *
*   THIS  ROUTINE CONTROLS THE SETTING UP OF THE VACUUM POLARIZATION   *
*   POTENTIAL FOR THE GIVEN NUCLEAR CHARGE DISTRIBUTION AT EACH GRID   *
*   POINT USING THE ANALYTIC FUNCTIONS DEFINED BY  L WAYNE FULLERTON   *
*   AND  G  A  RINKER  JR  IN  PHYS REV A 13 (1976) 1283. ACTUALLY R   *
*   TIMES  THEIR  POTENTIAL  IS  TABULATED, AS WE SHALL ONLY REQUIRE   *
*   THIS WHEN WE INTEGRATE. THE POTENTIAL IS  ACCUMULATED  IN  ARRAY   *
*   TB(I),I = 1,N WHICH IS IN COMMON BLOCK /TATB/.                     *
*                                                                      *
*   SUBROUTINES CALLED: VAC2, VAC4, VACUSR, NCHARG.                    *
*                                                                      *
*                                         LAST REVISION: 08 FEB 1988   *
*                                                                      *
************************************************************************
#include "implicit.h"
C      
      DIMENSION TB(N), R(N)  
C     
C   SECOND-ORDER VACUUM POLARISATION
C
      IF(NVPORD.EQ.1)THEN
        CALL VAC2(TB,R,Z,GEXP,N,NNUC,NUCTYP)
      ENDIF   
C
C   FOURTH-ORDER VACUUM POLARIZATION ... not implemented yet
C
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE VAC2(TB,R,Z,GEXP,N,NNUC,NUCTYP)
*                                                                      *
*   ----------------   SECTION 11   SUBPROGRAM 11  -----------------   *
*                                                                      *
*   THIS ROUTINE SETS UP THE SECOND-ORDER VACUUM POLARIZATION POTEN-   *
*   TIAL  USING  EQUATIONS (1), (4) AND (6) OF FULLERTON AND RINKER.   *
*   THE POTENTIAL  IS ACCUMULATED IN ARRAY TB(I),I = 1,N WHICH IS IN   *
*   COMMON BLOCK /TATB/.                                               *
*                                                                      *
*   SUBROUTINES CALLED: FUNK, QUAD.                                    *
*                                                                      *
*                                         LAST REVISION: 02 OCT 1989   *
*                                                                      *
************************************************************************
      use num_grid_cfg

#include "implicit.h"
#include "pi.h"
#include "dcbgen.h"
C     
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0)     
      DIMENSION TB(N),R(N)
      real(8), allocatable :: ZCHARG(:),TA(:),RG(:)
C      
      PI3 = THREE*PI
      C2 = CVAL*CVAL
      EPS = 1.0D-05
      WD = CVAL+CVAL
      FACTOR = -(Z+Z)/(PI3*CVAL)
C
C   SET UP POTENTIAL FOR POINT NUCLEUS
C
      DO 1 I = 1,N
         X = WD*R(I)
         TB(I) = FACTOR*FUNK (X,1)
    1 CONTINUE
      IF (NUCTYP .EQ. 1) RETURN
C
C   OTHER DISTRIBUTIONS (EQUATION (6))
C
    2 WB = (C2+C2)/(Z*THREE)
      IF(num_grid_cfg_grasp)THEN
         NG = N
         allocate(RG(NG))
         CALL DCOPY(NG,R,1,RG,1)
      ELSE
        NG = 4000  
        CALL MEMSET_NUM(Z,IATOM,NG,RNT,H,EPH,GRDC,
     &                   .TRUE.,1)
        allocate(RG(NG))
        RG(1) = RNT
        DO J = 2,NG
           RG(J) = RG(J-1)*EPH
        ENDDO        
      ENDIF   
!
      allocate(TA(NG))
      allocate(ZCHARG(NG))
      CALL DZERO(ZCHARG,NG)
      CALL DZERO(TA,NG)

C     Nuclear charge distribution (e.g., Gauss)
      CALL NCHARG(ZCHARG,RG,Z,GEXP,NG,NNUC,NUCTYP)
!      
      DO 3 I = 1,NG
        TA(I) = ZCHARG(I)*(RG(I)*RG(I)*RG(I))
    3 CONTINUE
      CALL QUAD (FAC3,TA,NG)
      FAC3 = FAC3*FACTOR*WB
      DO 4 I = 1,NG
        TA(I) = TA(I)*(RG(I)*RG(I))
    4 CONTINUE
      CALL QUAD (FAC5,TA,NG)
      FAC5 = FAC5*FACTOR*C2*WB/5.0D00
C
C   TABULATE POTENTIALS 
C
    5 DO 6 I = 1,N
         X = WD*R(I)
         TB(I) = TB(I)+FAC3*FUNK (X,3)+FAC5*FUNK (X,5) 
    6 CONTINUE
C
C   TABULATE POTENTIAL BY INTEGRATING EQUATION (4)
C
      FACTOR = -ONE/(PI3*C2*TWO)
      DO 8 K = 1,N
         XK = WD*R(K)
         DO 7 I = 1,NNUC
            XI = WD*RG(I)
            XM = ABS (XK-XI)
            XP = XK+XI
            TA(I) = (FUNK (XM,0)-FUNK (XP,0))*ZCHARG(I)
    7    CONTINUE
         CALL QUAD (X,TA,NG)
         X = X*FACTOR
         IF (X .EQ. ZERO) GO TO 11
         IF (ABS ((TB(K)-X)/X) .LT. EPS) RETURN
         TB(K) = X
    8 CONTINUE
      RETURN
 11   DO 12 I = K,N
         TB(I) = ZERO
 12   CONTINUE 
      deallocate(RG)
      deallocate(TA)
      deallocate(ZCHARG)
      RETURN
C
      END

************************************************************************
*                                                                      * 
C /* Deck FUNK */
      FUNCTION FUNK (X,N)
*                                                                      *
*   ----------------   SECTION 11   SUBPROGRAM 12  -----------------   *
*                                                                      *
*   THIS  FUNCTION  EVALUATES THE KN(X) FUNCTIONS USING THE ANALYTIC   *
*   FUNCTIONS DEFINED IN TABLES 1 AND 3 OF FULLERTON AND RINKER.       *
*                                                                      *
*   NO SUBROUTINES CALLED.                                             *
*                                                                      *
*                                         LAST REVISION: 06 AUG 1987   *
*                                                                      *
************************************************************************
#include "implicit.h"
C
      DIMENSION P(10,4),B(2,4),C(3,4),D(5,4),E(5,4),NP(4)
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0)
C     
#include "priunit.h"
C
      DATA P/ 8.8357293375D-01, -2.8259817381D-01,
     &       -5.8904879578D-01,  1.2500133434D-01,
     &       -3.2729913852D-02,  8.2888574511D-03,
     &       -1.0327765800D-05,  6.3643668900D-05,
     &        0.0D00, 0.0D00,
     &       -7.1740181754D-01,  1.1780972274D00,
     &       -3.7499963087D-01,  1.3089675530D-01,
     &       -3.8258286439D-02, -2.4297287300D-05,
     &       -3.5920148670D-04, -1.7170090700D-05,
     &        0.0D00, 0.0D00,
     &        9.9999999987D-01,  1.9770200000D-08,
     &       -7.5000050190D-01,  7.8540306316D-01,
     &       -3.4988601655D-01,  6.4596333000D-05,
     &       -9.8189080747D-03,  8.6513145800D-05,
     &       -2.3969236620D-04,  0.0D00,
     &        6.0000000002D 00, -6.4305200000D-08,
     &        2.1049413000D-06, -2.6711271500D-05,
     &       -1.3705236152D-01, -6.3476104090D-04,
     &       -7.8739801501D-02, -1.9641740173D-03,
     &       -3.4752369349D-03, -7.3145316220D-04/
C
      DATA B/-3.19999594323D+02,  2.53900995981D 00,
     &       -6.40514843293D+01,  7.11722714285D-01,
     &        5.19010136460D+03,  8.28495496200D+01,
     &        3.18150793824D+02,  4.33898867347D+01/
C
      DATA C/-3.19999594333D+02,  2.53901020662D+00,
     &        0.0D 00,            6.40514843287D+01,
     &       -7.11722686403D-01,  8.042207748D-04,
     &        2.76805406060D+04, -3.27039477790D+02,
     &        0.0D00,             8.48402116837D+02,
     &       -2.56939867765D+01,  3.20844906346D-01/
C
      DATA D/ 5.018065179D+00,  7.151891262D+01,
     &        2.116209929D+02,  3.140327478D+01,
     &       -1.0D00,           2.172386409D+02,
     &        1.643364528D+03,  2.122244512D+03,
     &       -4.512004044D+01,  1.0D00,
     &        8.540770444D 00,  6.076242766D+01,
     &        9.714630584D+01,  3.154973593D+01,
     &        1.0D00,           5.9243015865D-01,
     &        2.0596312871D00,  3.7785190424D 00,
     &        3.5614853214D00,  1.0D00/
C
      DATA E/ 2.669207401D+00,  5.172549669D+01,
     &        2.969809720D+02,  5.364324164D+02,
     &        1.535335924D+02,  1.155589983D+02,
     &        1.292191441D+03,  3.831198012D+03,
     &        2.904410075D+03,  0.000D00,
     &        4.543392478D00,   3.514920169D+01,
     &        6.019668656D+01,  8.468839579D00,
     &        0.0D00,           3.1511867816D-01,
     &        3.473245222D-01,  3.8791936870D-02,
     &       -1.3059741497D-03, 0.0D00/
C
      DATA NP/8,8,9,10/
C
      IF (X .EQ. ZERO) GO TO 11
      IF ((N .LT. 0) .OR. (N .EQ. 2) .OR.
     &    (N .EQ. 4) .OR. (N .GT. 5)) GO TO 98
      IF (N-3) 1,2,3
    1 K = N+1
      XN = ONE
      GO TO 4
    2 K = N
      XN = ONE/(X*X)
!      IF(N.eq.3)THEN
!         write(*,*)"BSUM,LOGX,CSUM",BSUM,LOG(X),CSUM
!         write(*,*)
!         write(*,*)"X*X,XN",X*X,XN
!      ENDIF   
      GO TO 4
    3 K = N-1
      XN = ONE/(X**4.0D0)
    4 IF (X .GT. ONE) GO TO 9
C
C   CALCULATE FUNCTION FOR X < = 1
C
      NN = NP(K)
      SUM = ZERO
      DO 5 I = 1,NN
         SUM = SUM+P(I,K)*XN
         XN = XN*X
    5 CONTINUE
      X2 = X*X
      BSUM = B(1,K)+X2*(B(2,K)+X2)
      CSUM = C(1,K)+X2*(C(2,K)+X2*C(3,K))
      GO TO (6,8,7,8), K
    6 BSUM = BSUM*X
      GO TO 8
    7 BSUM = BSUM*X2
    8 CONTINUE
!      IF(N.eq.3)THEN
!         write(*,*)"BSUM,LOGX,CSUM",BSUM,LOG(X),CSUM
!         write(*,*)
!         write(*,*)"SUM1",SUM
!      ENDIF   
      SUM = SUM+BSUM*LOG (X)/CSUM
!      IF(N.eq.3)THEN
!         write(*,*)"BSUM,LOGX,CSUM",BSUM,LOG(X),CSUM
!         write(*,*)"SUM2",SUM
!      ENDIF   
      FUNK = SUM
      RETURN
C
C   CALCULATE FUNCTION FOR X > 1
C
    9 XN = ONE
      DSUM = ZERO
      ESUM = ZERO
      IF (X .GT. 1.0D+02) GO TO 12
      DO 10 I = 1,5
         DSUM = DSUM+D(I,K)*XN
         ESUM = ESUM+E(I,K)*XN
         XN = XN/X
   10 CONTINUE
      XM = -X
      SUM = DSUM*EXP (XM)/(ESUM*SQRT (X**3))
      FUNK = SUM
      RETURN
   11 IF (N .NE. 0) GO TO 99
      FUNK = P(1,1)
      RETURN
   12 FUNK = ZERO
      RETURN
C
C   ERROR SECTION
C
   98 WRITE (LUPRI,300)
      WRITE (LUPRI,302)
      STOP 'Error sentense'
   99 WRITE (LUPRI,300)
      WRITE (LUPRI,301)
      STOP 'Error sentense'
C
  300 FORMAT (///' ***** ERROR IN FUNK *****')
  301 FORMAT (/' ATTEMPT TO CALCULATE FUNK(0,N) FOR N > 0')
  302 FORMAT (/' ATTEMPT TO CALCULATE FUNK(X,N) FOR N VALUES OTHER THAN'
     & ,' 0, 1, 3 AND 5.')
C
      END
************************************************************************
*                                                                      *
      SUBROUTINE NCHARG(TA,R,Z,GEXP,N,NNUC,NUCTYP)
*                                                                      *
*   ----------------   SECTION 11   SUBPROGRAM 19   ----------------   *
*                                                                      *
*   EVALUATE  THE normalized NUCLEAR CHARGE DISTRIBUTION (TZ).         *
*                                                                      *
*   (Fermi is UNnormalized distribution (10/04/2019)       *
*                                                                      *
*   NUCTYP  =  1 : Point                                               *
*              2 : GAUSSIAN                                            *
*              3 : FERMI                                               *
*                                                                      *
*   NO SUBROUTINES CALLED.                                             *
*                                                                      *
*                                         LAST REVISION: 08 FEB 1988   *
*                                                                      *
*   Adapted to DIRAC A. Sunaga Apr. 2019                               *
*                                                                      *
************************************************************************
      use num_grid_cfg
      use codata
      
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0)     
      DIMENSION TA(N),R(N)
C
      IF (GEXP .LE. ZERO) GO TO 99
      IF (NUCTYP .EQ. 1) GO TO 200
C
C   GAUSSIAN CHARGE DISTRIBUTION
C
      IF (NUCTYP .EQ. 2)THEN
        NNUC = N
        WB = GEXP                 !THREE/(TWO*RAD*RAD)
        WA = SQRT(WB)
        ZNORM = 4.0D0*Z
        ZNORM = ZNORM*WB*WA/SQRT(PI) ! ZNORM = Z * (xi/pi)^3/2 * 4pi
        DO 6 I = 1,N
           TA(I) = ( EXP(-WB*R(I)*R(I)) ) * R(I)*R(I)        
  6     CONTINUE
C Normalization
  8     CALL DSCAL(NNUC,ZNORM,TA,1)
     
C
C   FERMI CHARGE DISTRIBUTION
C
      ELSEIF(NUCTYP.EQ.3)THEN
C Ayaki: Normalization has not been implemented yet
        GO TO 9  
        TFM = 2.30D00
        T = TFM*1.0D-5/xtang  !xtang = 0.5291772083d0. See gp/codata.F90
        SKI = T/(4.0D0*LOG (3.0d0))      
        IF (SKI .LE. ZERO) GO TO 98
        write(*,*)"NNUC",NNUC
        DO 4 I = 1,NNUC
           TA(I) = (R(I)*R(I))/(ONE+EXP ((R(I)-RAD)/SKI))
!           write(*,*)"R,RAD,SKI",R(I),RAD,SKI,TA(I)
    4   CONTINUE
      ENDIF

      RETURN
C      
 200  write(*,*)'Point charge should not use NCHARGE routine.
     &           Something is wrong.'
      stop
 9    write(*,*)'Fermi has not been made in NCHARGE'
      stop

C      
C     ERROR SECTION
C
   97 WRITE (LUPRI,300)
      WRITE (LUPRI,302)
      STOP '111902'
   98 WRITE (LUPRI,300)
      WRITE (LUPRI,303)
      STOP '111903'
   99 WRITE (LUPRI,300)
      WRITE (LUPRI,304)
      STOP '111904'
C
  300 FORMAT (///' ***** ERROR IN NCHARG *****')
  301 FORMAT (/' NO CODE FOR NUCTYP  =  ',I2)
  302 FORMAT (/' SIX-PARAMETER DISTRIBUTIONS NOT IMPLEMENTED')
  303 FORMAT (/' SKIN THICKNESS MUST BE POSITIVE')
  304 FORMAT (/' NUCLEAR RADIUS MUST BE POSITIVE')
C
      END
C
*****************************************************************************
*
      SUBROUTINE INTMG(TB,R,Z,GEXP,N,NUCTYP)
*
*     PHI_MG CONTRIBUTION (7) of Flambaum et al. Phys. Rev. A 72 052115 (2005)
*     Magnetic contribution
*
*     Called from: NUMINT
*      
*     Modified to DIRAC A. Sunaga Mar. 2019      
*           
*****************************************************************************
#include "implicit.h"
#include "dcbgen.h"
#include "pi.h"
C     
      DIMENSION TB(N),R(N)
      PARAMETER(ONE=1.0D0,TWO=2.0D0)
C
      IF(NUCTYP.EQ.1)THEN
*     POINT
        POTIR=Z 
        DO 10 I = 1,N
          DIFFNUCPOT = POTIR/(R(I)*R(I))
          CALL INTMG2(TB(I),R(I),POTIR,DIFFNUCPOT)
 10     CONTINUE
      ELSE   
*     GAUSS
        XI = GEXP              !1.5/(RAD*RAD)
        WA = SQRT(XI)
        DO 20 I = 1,N
          POTIR = Z*ERF(WA*R(I))
          DIFFNUCPOT = POTIR/(R(I)*R(I))
     &                 - TWO*Z*SQRT(XI)/(SQRTPI*R(I))*EXP(-XI*R(I)*R(I))
          CALL INTMG2(TB(I),R(I),POTIR,DIFFNUCPOT)
 20     CONTINUE
      ENDIF  
C      
      RETURN
      END
*****************************************************************************
*
      SUBROUTINE INTMG2(TBI,RI,POTIR,DIFFNUCPOT)
*
*     PHI_MG CONTRIBUTION (7) of Flambaum et al. Phys. Rev. A 72 052115 (2005)
*     Magnetic contribution
*
*     Called from: NUMINT
*      
*     Modified to DIRAC A. Sunaga Mar. 2019      
*           
*****************************************************************************
#include "implicit.h"
#include "dcbgen.h"
#include "pi.h"
C     
      PARAMETER(NTMAX = 10000,ONE=1.0D0,TWO=2.0D0)
      real(8), allocatable :: TT(:),TTX(:)      
C
      allocate(TT(NTMAX))
      allocate(TTX(NTMAX))
C
      PREFAC = ONE/(4.0D0*CVAL*CVAL*PI)
C
C     PREFAC is different from GRASP.
C     In GRASP, large*small + small*large are calculated.
C     In DIRAC, only triangle matrix (large*small) are calculated.
C     From this, we have to remove TWO from GRASP's source.
C      
      CALL G1TGRID(TT,TTX,NTMAX,NT,RI)
      CALL QUAD_T_US2(TT,TTX,NTMAX,NT,RESULT)

*     radial derivative of the nuclear potential
      G1A = RESULT*DIFFNUCPOT
      G1B = DIFFNUCPOT
*
      CALL G2TGRID(TT,TTX,NTMAX,NT,RI)
      CALL QUAD_T_US2(TT,TTX,NTMAX,NT,RESULT)
      G2 = TWO*RESULT*CVAL*POTIR/RI
      TBI = PREFAC*(G1A - G1B + G2) 
C      
      deallocate(TT)
      deallocate(TTX)
C      
      RETURN
      END
************************************************************************
*                                                                      *
         SUBROUTINE G1TGRID (F,T,NTMAX,NT,R)
*                                                                      *
*   ----------------   SECTION 11   SUBPROGRAM 27   ----------------   *
*                                                                      *
*   CREATE T-GRID FOR THE RADIATIVE POTENTIALS                         *
*   SUBROUTINES CALLED:   PHI_G1_INNER                                 *
*                                                                      *
*                                         LAST REVISION: 15 DEC 2008   *
*                                                                      *
************************************************************************
*
#include <implicit.h>
*
      PARAMETER (ONE = 1.0D0, TEN = 10.0D0)
      INTEGER NT,NTMAX
      DOUBLEPRECISION T(NTMAX)
      DOUBLEPRECISION F(NTMAX)
      TWELVE  = 12.0D0
      HUNDRED = 100.0D0
*
      DO 1 J = 1, NTMAX
        T(J) = ONE + TEN**(-(TWELVE - DBLE(J)*ONE/HUNDRED))
        F(J) = PHI_G1_INNER(T(J),R)
*
        IF (ABS(F(J)) .LE. 1.0D-32) THEN
           NT = J
           GOTO 2
        END IF
*
    1 CONTINUE
    2 RETURN
*
      END
************************************************************************
*                                                                      *
      FUNCTION PHI_G1_INNER (T,R)
*                                                                      *
*   ----------------   SECTION 11   SUBPROGRAM 22  -----------------   *
*                                                                      *
*   CALCULATE THE FIRST PART OF THE INTEGRAL-KERNAL FOR THE RADIATIVE  *
*   PHI_G POTENTIAL NO SUBROUTINES CALLED.                             *
*                                                                      *
*                                         LAST REVISION: 12 DEC 2007   *
*                                                                      *
************************************************************************
*
#include "implicit.h"
#include "dcbgen.h"
*      
      PARAMETER (ONE = 1.0D0,TWO=2.0D0)
*     
* CALCULATE THE T DEPENDENT PART
*
      X1 = ONE/(T*T*SQRT(T*T - ONE))
      X2 = EXP(-TWO*T*R*CVAL)
*
      PHI_G1_INNER = X1*X2
*
      RETURN
*
      END      
************************************************************************
*                                                                      *
         SUBROUTINE QUAD_T_US2 (F,T,NTMAX,NT,RESULT)
*                                                                      *
*   ----------------   SECTION 11   SUBPROGRAM 24   ----------------   *
*                                                                      *
*   SIMPSON RULE FOR AN NONUNIFORM GRID                                *
*   NO SUBROUTINES CALLED.                                             *
*                                                                      *
*                                         LAST REVISION: 12 DEC 2007   *
*                                                                      *
************************************************************************
*
#include "implicit.h"
*
      PARAMETER( W7 = 1.0D0/90.0D0,
     &           W6 = -4.0D0*W7,  W5 =   7.0D0*W7,
     &           W4 = 30.0D0*W7,  W3 =  44.0D0*W7,
     &           W2 = 68.0D0*W7,  W1 = 112.0D0*W7)
*
      DOUBLEPRECISION T(NTMAX)
      DOUBLEPRECISION F(NTMAX)
*
*   FIND THE STARTING POINT FOR THE INTEGRATION. THIS IS THE POINT
*   WHERE THE INTEGRAND TA FIRST BECOMES SIGNIFICANT, AND IS CHOSEN
*   TO PREVENT DIVIDING BY ZERO.
*
      SUM = 0.0D0
      DO 1 J = 1,NT-3,2
         J1 = J
         J2 = J + 1
         J3 = J + 2
         F1=F(J1)
         F2=F(J2)
         F3=F(J3)
         T1=T(J1)
         T2=T(J2)
         T3=T(J3)
         X1 = F2*(T1-T3)*(T1-T3)
         X2 = F1*(T2-T3)*(2.0D0*T1-3*T2+T3)
         X3 = F3*(T1-T2)*(T1-3.0D0*T2+2.0D0*T3)
         SUM = SUM - (T1-T3)*(X1+X2-X3)/(6.0D0*(T1-T2)*(T2-T3))
    1 CONTINUE
      RESULT = SUM
      RETURN
*
      END
************************************************************************
*                                                                      *
         SUBROUTINE G2TGRID (F,T,NTMAX,NT,R)
*                                                                      *
*   ----------------   SECTION 11   SUBPROGRAM 28   ----------------   *
*                                                                      *
*   CREATE T-GRID FOR THE RADIATIVE POTENTIALS                         *
*   SUBROUTINES CALLED:   PHI_G2_INNER                                 *
*                                                                      *
*                                         LAST REVISION: 15 JAN 2008   *
*                                                                      *
************************************************************************
*
#include <implicit.h>
*
      PARAMETER (ONE = 1.0D0,TEN = 10.0D0)
      INTEGER NT,NTMAX
      DOUBLEPRECISION T(NTMAX)
      DOUBLEPRECISION F(NTMAX)
      TWELVE  = 12.0D0
      HUNDRED = 100.0D0
*
      DO 1 J = 1, NTMAX
        T(J) = ONE + TEN**(-(TWELVE - DBLE(J)*ONE/HUNDRED))
        F(J) = PHI_G2_INNER(T(J),R)
        IF (ABS(F(J)) .LE. 1.0D-32) THEN
           NT = J
           GOTO 2
        END IF
*
    1 CONTINUE
    2 RETURN
*
      END
************************************************************************
*                                                                      *
      FUNCTION PHI_G2_INNER (T,R)
*                                                                      *
*   ----------------   SECTION 11   SUBPROGRAM 23  -----------------   *
*                                                                      *
*   CALCULATE THE SECOND PART OF THE INTEGRAL-KERNAL FOR THE           *
*   RADIATIVE PHI_G POTENTIAL NO SUBROUTINES CALLED.                   *
*                                                                      *
*                                         LAST REVISION: 12 DEC 2007   *
*                                                                      *
************************************************************************
*
#include "implicit.h"
#include "dcbgen.h"
*
      PARAMETER (ONE = 1.0D0,TWO=2.0D0)
*     CALCULATE THE T DEPENDENT PART
*
      X1 = ONE/(T*SQRT(T*T - ONE))
      X2 = EXP(-TWO*T*R*CVAL)
*
      PHI_G2_INNER = X1*X2
*
      RETURN
*
      END      
*****************************************************************************
*
      SUBROUTINE INTLF(TB,R,Z,N)
*
*     PHI_L CONTRIBUTION (9) of Flambaum et al. Phys. Rev. A 72 052115 (2005)
*     Low flequency contribution
*
*     Called from: NUMINT
*      
*     Modified to DIRAC A. Sunaga Mar. 2019      
*           
*****************************************************************************
#include "implicit.h"
#include "dcbgen.h"
C      
      DIMENSION TB(N),R(N)
C
      BZ = 0.074D0 + 0.35D0*Z/CVAL
      CVAL3 = CVAL*CVAL*CVAL
      PREFAC = BZ * Z*Z*Z*Z / CVAL3
      DO 1 I = 1,N
         TB(I) = EXP(-Z*R(I))
 1    CONTINUE
      CALL DSCAL(N,PREFAC,TB,1)
C      
      RETURN
      END
*****************************************************************************
*
      SUBROUTINE INTHF(TB,R,Z,GEXP,N,NUCTYP)
*
*     PHI_L CONTRIBUTION (10) of Flambaum et al. Phys. Rev. A 72 052115 (2005)
*     High flequency contribution
*
*     Called from: NUMINT
*      
*     Modified to DIRAC A. Sunaga Mar. 2019      
*           
*****************************************************************************
#include "implicit.h"
#include "dcbgen.h"
#include "pi.h"
C     
      DIMENSION TB(N),R(N)
      PARAMETER(ONE = 1.0D0)
      PARAMETER(NTMAX = 10000)
      real(8), allocatable :: TT(:),TTX(:),AZTEST(:)      
C
      allocate(TT(NTMAX))
      allocate(TTX(NTMAX))
      allocate(AZTEST(120))
C     
      X1  = (Z-80.0D0)/CVAL
      X2  = X1 * X1
      X3  = X2 * X1
      X4  = X3 * X1
      ZC7 = 0.07D0 * Z*Z/(CVAL*CVAL)      
      AZ1COEF = 1.071D0 - 1.976D0*X2 - 2.128D0*X3 + 0.169D0*X4
      PREFAC  = ONE/(CVAL*PI)
      IF(NUCTYP.EQ.2)THEN
        XI = GEXP !1.5D0/(RAD*RAD)
        WA = SQRT(XI)
      ENDIF
C
      CALL DCOPY(N,R,1,TB,1)
      CALL DSCAL(N,CVAL,TB,1) 
C
      DO 5 I = 1,N
        IF (NUCTYP .EQ. 1) THEN
*     POINT
          POT = Z/R(I) 
        ELSE
*     GAUSS
          POT = Z*ERF(WA*R(I))/R(I)
        ENDIF
c     hth  GRID FOR T-INTEGRATION
        CALL FTGRID(TT,TTX,NTMAX,NT,R(I),Z)
chth  T-INTEGRATION
        CALL QUAD_T_US2(TT,TTX,NTMAX,NT,RESULT)
        AZ1 = AZ1COEF * TB(I)
        AZ2 = TB(I) + ZC7
        AZ  = AZ1/AZ2
        TB(I) = PREFAC*RESULT*POT*AZ
        
 5    CONTINUE
C
C     It does not agree with GRASP, because the default of GRASP is n-dependent HF
C     In order to compare with GRASP, we have to use n-INdependent HF term in GRASP
      deallocate(TT)
      deallocate(TTX)
      deallocate(AZTEST)
C      
      RETURN
      END
C
************************************************************************
*                                                                      *
         SUBROUTINE FTGRID (F,T,NTMAX,NT,R,Z)
*                                                                      *
*   ----------------   SECTION 11   SUBPROGRAM 26   ----------------   *
*                                                                      *
*   CREATE T-GRID FOR THE RADIATIVE POTENTIALS                         *
*   SUBROUTINES CALLED:   PHI_F_INNER                                  *
*                                                                      *
*                                         LAST REVISION: 12 DEC 2007   *
*                                                                      *
************************************************************************
*
#include <implicit.h>
*
      INTEGER NT,NTMAX,J
      PARAMETER(ONE = 1.0D0, TEN = 10.D0)
      DOUBLEPRECISION T(NTMAX)
      DOUBLEPRECISION F(NTMAX)
      TWELVE  = 12.0D0
      HUNDRED = 100.0D0
*
      DO 1 J = 1, NTMAX
        T(J) = ONE + TEN**(-(TWELVE - DBLE(J)*ONE/HUNDRED))
        F(J) = PHI_F_INNER(T(J),R,Z)
        IF (ABS(F(J)) .LE. 1.0D-32) THEN
           NT = J
           GOTO 2
        END IF
*
    1 CONTINUE
    2 RETURN
*
      END
************************************************************************
*                                                                      *
      FUNCTION PHI_F_INNER (T,R,Z)
*                                                                      *
*   ----------------   SECTION 11   SUBPROGRAM 21  -----------------   *
*                                                                      *
*   CALCULATE THE INTEGRAL-KERNAL FOR THE RADIATIVE PHI_F POTENTIAL    *
*   NO SUBROUTINES CALLED.                                             *
*                                                                      *
*                                         LAST REVISION: 12 DEC 2007   *
*                                                                      *
************************************************************************
*
#include "implicit.h"
#include "dcbgen.h"
*      
      PARAMETER(ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0, FOUR = 4.0D0,
     &          HALF = 0.5D0)
*     
* CALCULATE THE T DEPENDENT PART
*
      X1 = ONE/SQRT(T*T-ONE)
      X2 = ONE - ONE/(TWO*T*T)
      X3 = LOG(T*T - ONE)
      X4 = FOUR*LOG(CVAL/Z + HALF)
      X5 = THREE/TWO
      X6 = ONE/(T*T)
      X7 = EXP(-TWO*T*R*CVAL)
*
      PHI_F_INNER = X1*(X2*(X3 + X4) - X5 + X6)*X7
*
      RETURN
*
      END
      
*****************************************************************************
*
      SUBROUTINE INTPZ(TB,R,Z,N)
*
*     Model potential of Pekka and Zhao of J Phys B 36 1469 (2004)
*
*     Called from: NUMINT
*      
*     Modified to DIRAC A. Sunaga Mar. 2019      
*           
*****************************************************************************
#include "implicit.h"
C     
      DIMENSION TB(N),R(N)
C
      B0 = -48.6166D0
      B1 = 1.53666D0
      B2 = 0.0301129D0
      BETA0 = -12751.3D0
      BETA1 = 916.038D0
      BETA2 = 5.7797D0
      BFIT = B0 + B1*Z + B2*Z*Z
      BETAFIT = BETA0 + BETA1*Z + BETA2*Z*Z
C
      DO 1 I = 1,N
         TB(I) = EXP(-BETAFIT*R(I)*R(I))
 1    CONTINUE
      CALL DSCAL(N,BFIT,TB,1)
C
C     change
C     error sentence for 29<Z<83
C      
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE QUAD(RESULT,TE,NR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Numerical radial integral for one dimension
C
C copied from /QUAD/grasp.F/GRASP
C Adapted to DIRAC A. Sunaga Nov. 2018       
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
C
      PARAMETER( W7 = 1.0D0/90.0D0,
     &           W6 = -4.0D0*W7,  W5 =   7.0D0*W7,
     &           W4 = 30.0D0*W7,  W3 =  44.0D0*W7,
     &           W2 = 68.0D0*W7,  W1 = 112.0D0*W7)
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      PARAMETER(H=1.0D-2,EPSI = 1.0D-35) !-35
!      REAL(8) TE(NR)      
!      DIMENSION TE(NR)      
      DIMENSION TE(*)
C     
C   FIND THE STARTING POINT FOR THE INTEGRATION. THIS IS THE POINT
C   WHERE THE INTEGRAND TE FIRST BECOMES SIGNIFICANT, AND IS CHOSEN
C   TO PREVENT DIVIDING BY ZERO.
C
      DO 1 J = 1,NR
        IF (ABS (TE(J)) .GT. EPSI) GO TO 2
    1 CONTINUE
      RESULT = ZERO
      RETURN
C
C   CALCULATE THE INTEGRAL USING EXTENDED SIMPSONS RULE, NEGLECTING THE
C   CONTRIBUTION FROM THE STARTING POINT. CHOOSE THE STARTING POINT
C   R(JJ) TO GIVE AN ODD NUMBER OF INTEGRATION POINTS.
C   THE CRAY VERSION OF THE CODE USES THE VECTOR FUNCTION SSUM TO
C   PERFORM THE SUMMATION.
C
    2 SUM = ZERO
      NP = NR-J+1
      ND2 = NP/2
      N1D2 = (NP+1)/2
      K = ND2-N1D2+1
      JJ = J+K
      JJ1 = JJ+1
      DO 3 I = JJ1,NR-3,2
         SUM = SUM + W1*TE(I) + W2*TE(I+1)
    3 CONTINUE
      SUM = SUM + W1*TE(NR-1)
      SUM = SUM + W4*( TE(JJ)   + TE(NR)   )
     &          + W5*( TE(JJ1)  + TE(NR-1) )
     &          + W6*( TE(JJ+2) + TE(NR-2) )
     &          + W7*( TE(JJ+3) + TE(NR-3) )
C
C   ADD THE CONTRIBUTION FROM 0 TO R(JJ), PLUS THE CONTRIBUTION
C   TO SIMPSONS RULE FROM TE(JJ).
C
      COP = TE(J+1)/TE(J)
      IF (COP .LE. ONE) THEN
         RPX = ZERO
      ELSE
         RPX = TE(JJ)/LOG (COP)
      ENDIF
      RESULT = H*(SUM+RPX)
      RETURN
C
      END      
************************************************************************
*                                                                      *
C     /* Deck XKONKN */
      FUNCTION XKONKN (X,N,TOL)
*                                                                      *
*   ----------------   SECTION 09   SUBPROGRAM 10   ----------------   *
*                                                                      *
*   EVALUATE THE SUM OF THE SERIES (-X)**(K+1)/K**N FROM K = 1         *
*   UNTIL THE PRESCRIBED TOLERANCE IS MET.                             *
*                                                                      *
*   NO SUBROUTINES CALLED                                              *
*                                                                      *
*                                         LAST REVISION: 30 AUG 1987   *
*                                                                      *
************************************************************************
#include "implicit.h"
#include "priunit.h"
C     
      PARAMETER(MXXK = 100) ! 100 in GRASP      
      DIMENSION DELTA(MXXK)
      PARAMETER(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0)
C
C   SUMMATION FOR K = 1,...,6
C
      XK = ONE
      FKN = ONE
      SUM = XK
      DO 1 K = 2,6
         XK = -XK*X
         FKN = DBLE (K**N)
         TERM = XK/FKN
         SUM = SUM+TERM
         IF (ABS (TERM/SUM) .LT. TOL) GO TO 4
    1 CONTINUE
C
C   REMAINDER OF SUM
C
      K = 6
      M = 7
    2 K = K+1
      KM = K-M
      IF (KM .GT. MXXK-1) GO TO 99
      XK = -XK*X
      FKN = DBLE (K**N)
      DELTA(KM+1) = XK*HALF/FKN
      DO 3 J = KM,1,-1
         DELTA(J) = HALF*(DELTA(J)+DELTA(J+1))
    3 CONTINUE
      SUM = SUM+DELTA(1)
      IF (ABS (DELTA(1)/SUM) .GT. TOL) GO TO 2
C
    4 XKONKN = SUM*X
      RETURN
C
   99 WRITE (LUPRI,300) MXXK

      STOP '091001'
C
  300 FORMAT(/' ***** ERROR IN XKONKN *****',//' OVERFLOW IN ARRAY ',
     & 'DELTA. RECOMPILE WITH VALUE OF SYMBOL XK GREATER THAN ', I3 )
C
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck numintop */
      SUBROUTINE NUMINTOP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      
C Welcome message of numerical integration code
C
C Called from: NUMINT
C
C A. Sunaga Aug. 2020
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use num_grid_cfg
#include "priunit.h"
C      
       CALL TITLER
     &   ('Numerical integration for one-electron operators','*',199)
       WRITE(LUPRI,100)
 100   FORMAT(
     &    3X,'General structure: ',/,
     &    3X,'    Ayaki Sunaga and Trond Saue ',//,
     &    3X,'effQED: ',/,
     &    3X,'    Ayaki Sunaga and Trond Saue',//,
     &    1X,79('*'))
       
       WRITE(LUPRI,200)
 200   FORMAT(/' This module is published in:',//
     &    4X,'effQED: A. Sunaga, M. Salman and T. Saue',/,
     &    20X,'J. Chem. Phys. 157 (2022) 164101 ',/,
     &    1X,79('*'),/)

       WRITE(LUPRI,'(A)')
     &   " ***** SET-UP for NUMERICAL INTEGRATION *****"
       write(lupri, '(a)') ' - radial quadrature according to'
C       
       IF(num_grid_cfg_grasp)THEN
         write(lupri, '(a)')
     &'     the GRASP code.'
       ELSE
         write(lupri, '(a)')
     &'     R. Lindh, et al, Theor. Chem. Acc. 106, 178 (2001).'
         write(lupri, '(a, e12.5)')
     &'     precision of radial quadrature set to',
     &           num_grid_cfg_radint
       ENDIF
C
       write(lupri, '(a)')
     &' - radial integration is truncated by'
       IF(num_grid_cfg_uehcut)THEN
         write(lupri, '(a)')
     &'     Uehling potential '
       ELSE
         write(lupri, '(a)')
     &"     Flambaum-Ginges's low-frequency term"
       ENDIF 
C       
       write(lupri, '(a)')
     &' - angular quadrature according to'
       write(lupri, '(a, i2)')
     &'     the Lebedev scheme exact up to order L = ',
     &      num_grid_cfg_ANGQED
       WRITE(LUPRI,'(/,A)')
     &" ***** SET-UP for EACH ATOM *****"
C       
      RETURN         
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck numinfo */
      SUBROUTINE NUMINFO(ISYMOP,IATOM,Z,RNT,RMAX,H,MXNP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      
C Printing the information for numerical integration 
C
C Called from: NUMINT, NUMINT_ANG
C
C A. Sunaga Nov. 2020
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use num_grid_cfg

      use mxcent_mod
      use nucdata_mod, only: NUCMASS
#include "implicit.h"
#include "priunit.h"
#include "nuclei.h"
C
      IF(ISYMOP.NE.0)RETURN
C      
      WRITE(LUPRI,'(A,A)')
     &" - For ", NAMN(IATOM)
      WRITE(LUPRI,'(A,I3,A,E8.2e2,A,E8.2e2)')
     &"     mass = ",INT(NUCMASS(INT(Z))),
     &", Rmin = ", RNT, ", Rmax = ",RMAX
      WRITE(LUPRI,'(A,E8.2e2,A,I4)')
     &"     stepsize = ",H,", # of radial grid = ",MXNP 
C
      IF(IATOM.EQ.NUCIND)THEN
         num_grid_cfg_numint_print = .FALSE.
         WRITE(LUPRI,'(A)')         
      ENDIF
      RETURN
      END       
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EFFECTIVE_DENS_NUC(WORK,LWORK,COMB)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Generates integrals used for the calculation of the effective electronic
C density at nuclei.
C
C Called from: ONEGEN
C
C Written by Robert van Meer, Sept. 2010, last edited Mar. 2011
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use mxcent_mod
#include "implicit.h"
#include "priunit.h"
#include "dcbham.h"
#include "nuclei.h"
#include "dcbbas.h"
#include "dcbxpr.h"
C     GEXPEFFDE2 contains nuclear gaussian exponents related to size change
C     EXPPLC = Exponent placeholder
      REAL*8 GEXPEFFDE2, EXPPLC, PREFACTOR
      CHARACTER RTNLBL(2)*8,COMB*4,OMITVNUC(2)*4,
     &          LABINT(15*MXCENT)*8
      LOGICAl DOINT(4), DOATOM(NUCIND)
      INTEGER ITYP
      real*8 prefac2(4)
      DIMENSION WORK(LWORK), GEXPEFFDE2(4), INTREP(9*MXCENT)
#include "memint.h"
      CALL QENTER('EFFECTIVE_DENS_NUC')

C     Memory allocation
      CALL MEMGET('REAL',KONEIN,4*NNBBASX,WORK,KFREE,LFREE)

      IF (BSS.OR.x2C.OR.ZORA.OR.LEVYLE) THEN
        WRITE(LUPRI,*) 'EFFECTIVE_DENS_NUC: BSS, x2C, ZORA and LEVY'//
     &  '-LEBLOND modes not supported'
        CALL QUIT('EFFECTIVE_DENS_NUC: not supported for current'//
     &  ' Hamiltonian')
      END IF

      READ(COMB,'(4L1)',ERR=1000) (DOINT(I),I=1,4)
      OMITVNUC(1) = 'FFFF'
      OMITVNUC(2) = 'FFFF'
C     Fetch and set labels/variables (-99 is used to turn on .SELECT)
      CALL SETATM(DOATOM,NATOM,-99)
      CALL EFFDE2TYP(NOPTYP,INTREP,LABINT,DOATOM,NATOM)
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = 'SY'
      WRITE(RTNLBL(2)(3:4),'(I2)') 1
      RTNLBL(2)(5:8) = COMB
      ITYP = 0

C     Get the integrals, exponent is changed for a single nucleus
C     and reset to the original value after each calculation
      DO K = 1,NUCIND
        IF (DOATOM(K)) THEN
C         Generate nuclear exponents for nuclear size displacments and prefactor
          CALL EFF_DEN_NUC_EXP(GEXPEFFDE2, K, PREFACTOR)
          EXPPLC = GNUEXP(K)
C         Calculate property matrices for nuclear displacements
          DO L = 1, 4
            GNUEXP(K) = GEXPEFFDE2(L)
C           Calculate the integrals
            CALL NUCAT1(WORK(KONEIN+(L-1)*NNBBASX),WORK(KFREE),LFREE,
     &                  DOINT,OMITVNUC,IPRHAM)
          ENDDO
          GNUEXP(K) = EXPPLC
C         Create total property matrix with correct prefactor
          CALL MEMGET('REAL',KHELP,NNBBASX,WORK,KFREE,LFREE)
          call dzero(work(khelp),nnbbasx)
          prefac2(1) = +1.0d0
          prefac2(2) = -8.0d0
          prefac2(3) = +8.0d0
          prefac2(4) = -1.0d0
          do l = 1, 4
          call daxpy(NNBBASX,PREFACTOR*prefac2(l),
     &               work(konein+(l-1)*nnbbasx),1,
     &               WORK(khelp),1)
          end do
C         Store matrix
          ITYP = ITYP + 1
          CALL WRTPRO(WORK(khelp),NNBBASX,LABINT(ITYP),RTNLBL,IPRHAM)
          CALL MEMREL('EFFECTIVE_DENS_NUC',WORK,KWORK,KHELP,KFREE,LFREE)
        ENDIF
      ENDDO

C     Memory deallocation
      CALL MEMREL('EFFECTIVE_DENS_NUC',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('EFFECTIVE_DENS_NUC')
      RETURN

 1000 CONTINUE
C     Not able to read DOINT information from COMB
      WRITE(LUPRI,'(A,A)') 'EFFECTIVE_DENS_NUC: Not able to read COMB ='
     & ,COMB
      CALL QUIT('EFFECTIVE_DENS_NUC: Not able to read COMB')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EFF_DEN_NUC_EXP(GEXPEFFDE2, NUCLEUS, PREFACTOR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Fills Mossbauer exponent array with exponent values
C
C Called from: MOSSBAUER_ISOMER
C
C Written by Robert van Meer, Sept. 2010 last edited Mar. 2011
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use mxcent_mod
#include "implicit.h"
#include "nuclei.h"
#include "pi.h"
      REAL*8 GEXPEFFDE2, RADIUS
      DIMENSION GEXPEFFDE2(4)
C     Generate nuclear radius from the gaussian exponent
      RADIUS = SQRT(3D0 / (2D0 * GNUEXP(NUCLEUS)) )
C     Generate new gaussian parameters with dR=0.005R
      GEXPEFFDE2(1) = 3D0 / (2D0 * 0.990D0 * 0.990D0 * RADIUS * RADIUS)
      GEXPEFFDE2(2) = 3D0 / (2D0 * 0.995D0 * 0.995D0 * RADIUS * RADIUS)
      GEXPEFFDE2(3) = 3D0 / (2D0 * 1.005D0 * 1.005D0 * RADIUS * RADIUS)
      GEXPEFFDE2(4) = 3D0 / (2D0 * 1.010D0 * 1.010D0 * RADIUS * RADIUS)
C     Prefactor = 3 / (48*Z*Pi*R*dR)
      PREFACTOR = 12.5D0 / (RADIUS * RADIUS * CHARGE(NUCLEUS) * PI)
!     write(lupri,*) 'gnuexp, radius, pi, prefac, new exp',
!    & GNUEXP(NUCLEUS), radius, pi, prefactor,GEXPEFFDE2(1:4)
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck op1gen */
      SUBROUTINE OP1GEN(WORD,COMB,IORDER,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Driver for the generation of one-electron integrals
C     to be placed on the file CHECKPOINT
C
C     On input:  WORD     ...
C                COMB     ... integral characteritics
C                IORDER   ...
C
C     Called from ONEGENdirone.F
C
C     Written by Trond Saue - May 1996
C     Last revisions:  May 19 1996
C                  MI(+TEC): minor changes for LAO based RPA, dec.2002
C
C*****************************************************************************
      use primitives_mod
      use maxaqn_mod
      use mxcent_mod
      use symmet_mod
#include "implicit.h"
#include "priunit.h"
C
#include "dummy.h"
#include "nuclei.h"
#include "cbiher.h"
      LOGICAL DOINT(4),TOFILE
      CHARACTER WORD*7,COMB*4,LABINT(5*MXCOOR)*8,MTFORM*6
      DIMENSION WORK(LWORK)
C
      CALL QENTER('OP1GEN')
#include "memint.h"
C     Write integrals to file
      TOFILE   = .TRUE.
      NCOMP    = 0
      MAXTYP = 5*MXCOOR
      IF ((WORD.EQ.'ELFGRDC').OR.(WORD.EQ.'ELFGRDS')) THEN
         NINTAD = 9*NUCIND*(MAXREP+1)
      ELSE IF (WORD.EQ.'DSO    ') THEN
CMI ... to be fixed for the future LAO implementation...
         NINTAD = 9*NUCDEP*NUCDEP
      ELSE
         NINTAD = MAXTYP
      END IF
      READ(COMB,'(4L1)',ERR=1000) (DOINT(I),I=1,4)
C
C     Memory allocation
C
      CALL MEMGET('INTE',KINTRP,MAXTYP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINTAD,NINTAD,WORK,KFREE,LFREE)
C
C  extra assignement of characteristics for very narrow group of integrals
C
      NCOMP = 0
      IF (WORD.EQ.'S1MAGR ') THEN
         TRIANG = .FALSE.
         MTFORM = 'SQUARE'
      ELSE
         TRIANG = .TRUE.
         MTFORM = 'TRIANG'
      END IF
C
C     Generate integrals in abacus
C
C     ... calling of ../abacus/her1pro.F
      CALL PR1INT_1(WORK(KFREE),LFREE,WORK(KINTRP),WORK(KINTAD),
     &            LABINT,WORD,IORDER,NPQUAD,TRIANG,
     &            PROPRI,IPRINT,DUMMY,NCOMP,TOFILE,MTFORM,
     &            DOINT)
C
C     Memory deallocation
C
      CALL MEMREL('OP1GEN',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('OP1GEN')
C
      RETURN
 1000 CONTINUE
C
C    Not able to read DOINT information from COMB
C
      WRITE(LUPRI,'(A,A)') 'OP1GEN: Class label:',WORD
      WRITE(LUPRI,'(A,A)') 'OP1GEN ERROR: Not able to read COMB =',COMB
      CALL QUIT('OP1GEN: Not able to read COMB')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck betamt */
      SUBROUTINE BETAMT(COMB,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Driver for the generation of beta-matrix integrals
C     to be placed on file CHECKPOINT
C
C     Written by Trond Saue - May 1996
C     Last revision May 19 1996
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dummy.h"
      LOGICAL DOINT(4)
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4
      DIMENSION WORK(LWORK)
#include "dcbbas.h"
C
      CALL QENTER('BETAMT')
#include "memint.h"
C
C     Initialization
C
      NCOMP    = 1
      READ(COMB,'(4L1)',ERR=1000) (DOINT(I),I=1,4)
C
C     Memory allocation
C
      CALL MEMGET('REAL',KONEIN,NNBBASX,WORK,KFREE,LFREE)
C
C     Generate integrals
C
      CALL RGETINT('OVERLAP',WORK(KONEIN),DOINT,NCOMP,IDUMMY,IDUMMY,
     &            'TRIANG',.FALSE.,IPRINT,WORK(KFREE),LFREE)
C
      IF(IPRINT.GE.6) THEN
        CALL HEADER('Small component overlap integrals',-1)
        CALL OUTPAK(WORK(KONEIN),NTBAS(0),1,LUPRI)
      ENDIF
C
C     Generate integral labels
C
      LABEL = 'BETAMAT '
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = 'SY'
      WRITE(RTNLBL(2)(3:4),'(I2)') 1
      RTNLBL(2)(5:8) = COMB
C
C     Write integrals to file
C
      CALL WRTPRO(WORK(KONEIN),NNBBASX,LABEL,RTNLBL,IPRINT)
C
C     Memory deallocation
C
      CALL MEMREL('BETAMT',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('BETAMT')
C
      RETURN
 1000 CONTINUE
C
C     Not able to read DOINT information from COMB
C
      WRITE(LUPRI,'(A,A)') 'BETAMT ERROR: Not able to read COMB =',COMB
      CALL QUIT('BETAMT: Not able to read COMB')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck nucatt */
      SUBROUTINE NUCATT(WORK,LWORK,COMB,IPRINT)
C*****************************************************************************
C
C     Generate nuclear attraction integrals
C
C     IPRINT is IPRONE
C
C     Called from: ONEGEN
C
C     Written by T.Saue May 1995
C     Last revision May 19 1996 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      LOGICAl DOINT(4)
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4,OMITVNUC(2)*4
      DIMENSION WORK(LWORK)
#include "dcbbas.h"
C
      CALL QENTER('NUCATT')
#include "memint.h"
C     Memory allocation
      CALL MEMGET('REAL',KONEIN,NNBBASX,WORK,KFREE,LFREE)
C
C     Generate integrals
C
      READ(COMB,'(4L1)',ERR=1000) (DOINT(I),I=1,4)
      OMITVNUC(1) = 'FFFF'
      OMITVNUC(2) = 'FFFF'
      CALL NUCAT1(WORK(KONEIN),WORK(KFREE),LFREE,
     &            DOINT,OMITVNUC,IPRINT)
      IF(IPRINT.GE.4) THEN
        CALL HEADER('Nuclear attraction Integrals',-1)
        CALL OUTPAK(WORK(KONEIN),NTBAS(0),1,LUPRI)
      ENDIF
C
C     Generate integral labels
C
      LABEL = 'MOLFIELD'
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = 'SY'
      WRITE(RTNLBL(2)(3:4),'(I2)') 1
      RTNLBL(2)(5:8) = COMB
C
C     Write integrals to file
C
      CALL WRTPRO(WORK(KONEIN),NNBBASX,LABEL,RTNLBL,IPRINT)
C
C     Memory deallocation
      CALL MEMREL('NUCATT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('NUCATT')
      RETURN
 1000 CONTINUE
C
C     Not able to read DOINT information from COMB
C
      WRITE(LUPRI,'(A,A)') 'NUCATT: Not able to read COMB =',COMB
      CALL QUIT('NUCATT: Not able to read COMB')
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck nucat1 */
      SUBROUTINE NUCAT1(OP1INT,WORK,LWORK,DOINT,OMITVNUC,IPRINT)
      use mxcent_mod
#include "implicit.h"
#include "dummy.h"
C
      LOGICAL DOINT(2,2)
      CHARACTER*4 OMITVNUC(2)
      DIMENSION OP1INT(NNBBASX),WORK(LWORK)
#include "dcbbas.h"
#include "nuclei.h"
#include "dcbham.h"
C
      CALL QENTER('NUCAT1')
#include "memint.h"
C
      CALL MEMGET('REAL',KSTHMA,3*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFACIN,2*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      LWRK   = LWORK - KFREE + 1
C
      CALL ONEDR1(WORK(KSTHMA),DUMMY,DUMMY,WORK(KFACIN),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KNCENT),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),DUMMY,WORK(KFREE),LWRK,
     &            IPRINT,.FALSE.,0,.TRUE.,.TRUE.,.FALSE.,NTBAS(0),
     &            NNBBASX,NNBBAST,DOINT,OMITVNUC,.FALSE.)
C
C     Symmetry unpack nuclear attraction integrals
C
      KNUCS = KSTHMA + NNBBASX
C
C      *********************************************
C      **                                         **
C      **      this part adds ECP AREP            **
C      **    nuc-elec attraction integrals.       **
C      **                                         **
C      *********************************************
C
      IF (ECPCALC) then
         OPEN(59,FORM='UNFORMATTED',FILE='RECP_INT_C')
         CALL MEMGET('REAL',KARPP,NNBBAST,WORK,KFREE,LFREE)
         CALL READT(59,NNBBAST,WORK(KARPP))
         close(59)
         DO I=0,NNBBAST-1
            WORK(KNUCS+I)=WORK(KNUCS+I)+WORK(KARPP+I)
         ENDDO
         CALL MEMREL('NUCAT1.ECPCALC',WORK,KWORK,KARPP,KFREE,LFREE)
      ENDIF
C
      CALL SYMUPK(WORK(KNUCS),OP1INT,1,NNBBASX)
C
C     Memory deallocation
C
      CALL MEMREL('NUCAT1.ONEDR1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('NUCAT1')
      RETURN
C
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EDM_EFF_OPER(WORK,LWORK,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  prepare imaginary, antisymmetric KINENERG integrals into LS/SL blocks,
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
      DIMENSION WORK(LWORK)
!     ... local variables
      LOGICAl DOINT(4)
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4,OMITVNUC(2)*4
#include "memint.h"

C     Memory allocation
      CALL MEMGET('REAL',K_EDM,NNBBASX,WORK,KFREE,LFREE)

      write(LUPRI,'(/,2X,A,/)') '....integral routine for EDM '

C     Memory deallocation
      CALL MEMREL('EDM',WORK,KWORK,KWORK,KFREE,LFREE)
      !CALL QUIT('EDM end')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck onefck */
      SUBROUTINE ONEFCK(FOCK,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C   PURPOSE:
C ===========
C     Generate one-electron Fock matrix, using the one-electron operators
C     specified in COMMON block CBIHAM.
C
C     It can be (in SA-AO basis):
C         i) Dirac 1-el. operator
C        ii) BSS/IOTC two-component operator (either in the 4comp. or in the reduced 2comp. framework)
C       iii) Non-relativistic operator (in reduced two-component framework)
C
C     Routine is called from many places.
C
C     NOTE: If BSS=.true., get BSS integrals (TWOCOMP does not matter later)
C              BSS=.false., TWOCOMP does matter (Dirac or nonrelativistic Ham.)
C
C     Written by Trond Saue May 1996
C     Last revision: May 19 1996
C
C MI&HJAaJ (2003): If BSS=.TRUE., read the TWO-COMPONENT
C                   relativistic bare nuclei (in AO basis) from the file
C                   BSSMAT (label='H1AO_DK')
C
C   MI/Strasbg,2005 added new labels containing
C        various modified atomic integrals of H2c AO
C   MI/Strasbg,Febr. 2006 H2cAO Hamiltonian reading placed into separate routine
C
C*****************************************************************************
      use x2cmod_cfg
      use x2c_utils, only:
     &    read_1fock_x2c
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbgen.h"
      DIMENSION FOCK(*),WORK(LWORK)
C
      CALL QENTER('ONEFCK')
#include "memint.h"

CMI   ... control output ...
      IF (IPRHAM.GE.4) THEN
        IF (BSS .or. x2c) THEN
         if(x2c)then
          CALL HEADER(
     &    'ONEFCK: One-electron Fock matrix from - '//
     &    ' X2C 2c Hamiltonian ',-1)
         end if
         IF (TWOCOMPBSS) THEN
          CALL HEADER(
     &    'ONEFCK: One-electron Fock matrix from - '//
     &    ' BSS 2c Hamiltonian (in pure two-comp mode, '//
     &    ' saving memory)',-1)
         ELSE
          CALL HEADER(
     &    'ONEFCK: One-electron Fock matrix from - '//
     &    ' BSS 2c Hamiltonian (in four-component mode, '//
     &    'wasting memory)',-1)
         ENDIF
        ELSE
          IF (TWOCOMP) THEN
         CALL HEADER(
     &   'ONEFCK: One-electron Fock matrix - nonrelativistic'//
     &   '  Hamiltonian (BSS=.false.,TWOCOMP=.true.)',-1)
          ELSEIF(MDIRAC) THEN
         CALL HEADER(
     &   'ONEFCK: One-electron Fock matrix from the'//
     &   ' modified Dirac Hamiltonian (BSS=.false.,TWOCOMP=.false.)',-1)
          ELSE
         CALL HEADER(
     &   'ONEFCK: One-electron Fock matrix from the'//
     &   ' Dirac Hamiltonian (BSS=.false.,TWOCOMP=.false.)',-1)
          ENDIF
       ENDIF
      ENDIF

      if(bss)then
!MI If BSS is active read one-electron DKn(BSS) AO integrals from the
!   BSSMAT file (lu=LUBSS) according to the parameter I2COFK
       CALL GETH2CAO(I2COFK,FOCK,WORK(KFREE),LFREE,IPRHAM)
      else if(x2c.and..not.x2c_4c_fock_mtx_defh1)then
        call read_1fock_x2c(fock,n2bbasx,nz,i2cofk)
      else
C
C     Memory allocation
C
cisl-b    getting NZ*N2BBASX instead of NZ*NNBBASX
C
C     ************************************************
C     ****** Construct one-electron Fock matrix ******
C     ************************************************
C
       CALL MEMGET('LOGI',KFIRST,NZ        ,WORK,KFREE,LFREE)
       CALL ONEFC1(FOCK,WORK(KFREE),LFREE)
       CALL MEMREL('ONEFCK.ONEFC1',WORK,KWORK,KWORK,KFREE,LFREE)

      end if
C
C     Output section
C     ==============
C
!#define blubb_blubb_so
#ifdef blubb_blubb_so
      iprham_save = iprham
      iprham      = 4
#endif
      IF(IPRHAM.GE.4) THEN
        IF (BSS) THEN
         IF (TWOCOMPBSS) THEN
          CALL HEADER(
     &    'ONEFCK: One-electron Fock (symm.blocked SA-AO)'//
     &    ' matrix from - '//
     &    ' BSS 2c Hamiltonian (in pure two-comp mode, '//
     &    ' saving memory)',-1)
         ELSE
          CALL HEADER(
     &    'ONEFCK: One-electron Fock (symm.blocked SA-AO)'//
     &    'matrix from - '//
     &    ' BSS 2c Hamiltonian (LL blocks in four-component mode, '//
     &    'wasting memory)',-1)
         ENDIF
        ELSE IF(x2c)THEN
          CALL HEADER(
     &    'ONEFCK: One-electron Fock (symm.blocked SA-AO)'//
     &    ' matrix from - X2C 2c Hamiltonian)',-1)
        ELSE
          IF (TWOCOMP) THEN
         CALL HEADER(
     &   'ONEFCK: One-electron Fock matrix - nonrelativistic'//
     &   '  Hamiltonian (BSS=.false.,TWOCOMP=.true.)',-1)
          ELSEIF(MDIRAC) THEN
         CALL HEADER(
     &   'ONEFCK: One-electron Fock matrix from the'//
     &   ' modified Dirac Hamiltonian (BSS=.false.,TWOCOMP=.false.)',-1)
          ELSE
         CALL HEADER(
     &   'ONEFCK: One-electron Fock matrix from the'//
     &   ' Dirac Hamiltonian (BSS=.false.,TWOCOMP=.false.)',-1)
          ENDIF
        ENDIF
        DO 10 I = 1,NFSYM
          IF(NFBAS(I,0).EQ.0) GOTO 10
          WRITE(LUPRI,'(A,I1,A,I1)')
     &     '* Fermion ircop no.',I,'/',NFSYM
          CALL PRQMAT(FOCK(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   10   CONTINUE
      ENDIF
#ifdef blubb_blubb_so
      iprham = iprham_save
#undef blubb_blubb_so
#endif
C
      CALL QEXIT('ONEFCK')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck onefc1 */
      SUBROUTINE ONEFC1(FOCK,WORK,LWORK)
C*****************************************************************************
C
C     Generate one-electron Fock matrix, using the one-electron operators
C     specified in COMMON block CBIHAM.
C
C     Called from: ONEFCK
C
C     Written by Trond Saue May 1996
C     Last revision: May 19 1996
C                  MI / Febr.2006
C
C*****************************************************************************
      use dirac_cfg
      use fde_mod

#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "dummy.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dcbxpr.h"
#include "dcbprl.h"
C
      type(fde_import) :: itmp
      real(8), allocatable :: ftri(:,:),op1int(:)
      logical, allocatable :: first(:)
      LOGICAL TOBE
      DIMENSION FOCK(N2BBASX,NZ),WORK(LWORK)
C
      CALL QENTER('ONEFCK1')
#include "memint.h"
C
C     Initialize
C
      OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &      FILE='AOPROPER')
      allocate(ftri(nnbbasx,nz))
      allocate(op1int(nnbbasx))
      allocate(first(nz))
      CALL LSET(NZ,.TRUE.,FIRST)
C
C     Nuclear attraction integrals
C     ============================
C
      IOP = IPMOLFLD
      CALL PRPMAO(LU1INT,IOP,.TRUE.,WORK,.TRUE.,NNBBASX,FTRI,OP1INT,
     &            FIRST,IPRHAM)
      IF (SMLV1C) THEN
         CALL CATCORR(FTRI,DUMMY,DUMMY,WORK,LFREE,IPRHAM)
         CALL PR1CEX1 (FTRI,FIRST,FOCK,WORK,KFREE,LFREE,IPRHAM)
      ENDIF

      if (dirac_cfg_fde) then
C     Embedding potential integrals, if FDE is requested
C     ==================================================
C
C aspg: i wonder whether or not there should be more control to allow
C       for the possibility of not including the static potential to
C       the one-electron matrix prior to the scf...
         call fde_get_import_info(itmp)
         if (itmp%im_vemb) then
            IOP = IPVEMB0
            CALL PRPMAO(LU1INT,IOP,.TRUE.,WORK,.TRUE.,NNBBASX,FTRI,
     &      OP1INT,FIRST,IPRHAM)
         endif
      ENDIF

C
C     Kinetic energy
C     ==============
C
CMI/Febr.2006 For TWOCOMP=.true. this is the exact nonrelativistic kinetic energy
C
      IOP = IPKINERG
      CALL PRPMAO(LU1INT,IOP,.TRUE.,WORK,.TRUE.,NNBBASX,FTRI,
     &            OP1INT,FIRST,IPRHAM)

C
C     Beta matrix
C     ===========
C
CMI/Febr.2006 Beta matrix is not initialized in OP1INI for the two-component mode !

      IF (.NOT.TWOCOMP) THEN
        IOP = IPBETAMT
        CALL PRPMAO(LU1INT,IOP,.TRUE.,WORK,.TRUE.,NNBBASX,FTRI,
     &              OP1INT,FIRST,IPRHAM)
      ENDIF
C
C     Energy shift
C     ============
C
      IF (.NOT.TWOCOMP) THEN
        IOP = IPESHIFT
        IF(IOP.NE.-1) THEN 
          CALL PRPMAO(LU1INT,IOP,.TRUE.,WORK,.TRUE.,NNBBASX,FTRI,
     &          OP1INT,FIRST,IPRHAM)
        ENDIF
      ENDIF

C
C     Finite field perturbations:
C     Loop through  the list of one-electron operators in CBIHAM;
C     ======================================================
C     ...skip Fock-matrices for now ...
      TOBE = .FALSE.
      DO I = 1,N1OPER
        IOP = IPR1OP(I)
        IOPTYP = IPRPTYP(IOP)
        IF(IOPTYP.EQ.0) THEN
C.........Fock matrix
          TOBE = .TRUE.
        ELSE
           CALL PRPMAO(LU1INT,IOP,.TRUE.,WORK,.TRUE.,NNBBASX,FTRI,
     &                 OP1INT,FIRST,IPRHAM)
        ENDIF
      ENDDO
C
C
C     Make full matrix
C     ================
C
      deallocate(op1int)
      CALL TRI2SQ(FOCK,FTRI,FIRST,0,1)
      deallocate(ftri)
c
c     now add SO integrals
c     SOPP integrals
c     ==============
c
      if (ECPCALC) then
         allocate(op1int(n2bbasx))
         IOP = IPSOPP
         CALL PRPMAO(LU1INT,IOP,.TRUE.,WORK,.FALSE.,N2BBASX,FOCK,
     &          OP1INT,FIRST,IPRHAM)
         deallocate(op1int)
      endif
C
C     Reindex to sorted basis
C     =======================
C
      CALL BUTOBS(FOCK,NZ,WORK(KFREE),LFREE)
C
C     Add atomic shift projectors
C     ---------------------------
C
      IF(DOASHIFT) CALL ASHIFT(FOCK,IPRHAM,WORK,KFREE,LFREE)      
C
C     Add effective QED potentials
C     ----------------------------      
C
      IF(DOEFFQED) CALL EFFQED(FOCK,WORK,KFREE,LFREE)
C
C     Add Fock matrices when requested
C
      IF(TOBE) THEN
        allocate(op1int(n2bbasx))
        DO I = 1,N1OPER
          IOP = IPR1OP(I)
          IOPTYP = IPRPTYP(IOP)
          IF(IOPTYP.EQ.0) THEN
            KLBL  = IPRPLBL(1,IOP)
            INQUIRE(FILE=PRPLBL(KLBL),EXIST=FIRST(1))
            IF(FIRST(1)) THEN
              CALL OPNFIL(2,PRPLBL(KLBL),'OLD','ONEFC1')
              CALL REAFCK(2,OP1INT,.FALSE.,1)
              CLOSE(2,STATUS = 'KEEP')
              PFAC = FACPRP(1,IOP)
              CALL DAXPY(N2BBASXQ,PFAC,OP1INT,1,FOCK,1)
              CALL LSET(NZ,.FALSE.,FIRST)
            ELSE
              WRITE(LUPRI,'(A,A8,A,A16,A)')
     &        'Property file ',PRPLBL(KLBL),' for operator ',
     &        PRPNAM(IOP),' not found !'
              CALL QUIT('ONEFC1: Fock matrix file not found !')
            ENDIF
          ENDIF
        ENDDO
        deallocate(op1int)
      ENDIF
      deallocate(first)
      CALL QEXIT('ONEFCK1')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpmao */
      SUBROUTINE PRPMAO(LU1INT,IOP,ONFILE,PRPINT,TRIANG,NDIM,PMAT,
     &                  OP1INT,FIRST,IPRINT)
C*****************************************************************************
C
C     Generates  lower triangular/square  matrix of one-electron integrals
C     corresponding to property IOP of COMMON block XCBPRP
C
C     On input:  LU1INT    - ONFILE: lu of integral file; CACHE: offset
C                IOP     - index of the property
C                NDIM    - dimension of array (either triangular or square)
C                FIRST   - array, must be set to true
C                IPRINT  - print level (very usefull for control printouts)
C
C     On output: PMAT - accumulated property matrix (Hermit sorting)
C
C
C     Written by Trond Saue May 1996
C     Last revision: May 19 1996 - tsaue
C     Dec.2005, MI - generalized for square properties (NDIM = NNBBASX or N2BBASX)
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,PTOL = 1.0D-12)
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbxpr.h"
#include "dcbprl.h"
C
!      IPRINT=16 !sunaga      
      CHARACTER MXFORM*6, FMT*6, RTNLBL(2)*8,JCHR*1
      LOGICAL FIRST(*), DOINT(4), TRIANG, ONFILE
      DIMENSION PMAT(NDIM,*),OP1INT(NDIM),PRPINT(NDIM,*),IDOINT(4)
C
C     Extract information about property
C     ==================================
C
      IOPTYP = IPRPTYP(IOP)
      IF(IOPTYP.EQ.0) THEN
        CALL QUIT('PRPMAO: No Fock matrix read here !')
      ENDIF
      IHREF = 1
      IF(.NOT.TRIANG) IHREF = 0
      IF(IPRINT.GE.5) THEN
        CALL HEADER('Output from PRPMAO',-1)
        WRITE(LUPRI,*)
     & 'Generating property matrix for the property: ',PRPNAM(IOP)
        IF (TRIANG) THEN
          WRITE(LUPRI,*) '.... it is a LOWER TRIANGULAR matrix '
        ELSE
          WRITE(LUPRI,*) '.... it is a SQUARE matrix '
        ENDIF
      ENDIF
C
C     Loop over components
C     ====================
C
      IBRPOP = IPRPSYM(IOP)-1
      NCMP   = MCMP(IOPTYP)
      DO 10 I = 1,NCMP
        PFAC  = FACPRP(I,IOP)
        IF(ABS(PFAC).LT.PTOL) GOTO 10 ! Skip this component of factor is very small
C
C       Extract information about component
C
C       SCALAR OPERATOR:
C       ----------------
C       KLBL  - pointer to list of labels
C       IH    - scalar operator is symmetric(1)/antisymmetric(2) about diagonal
C       IBREP - boson irrep of scalar operator
C
C       MATRIX OPERATOR:
C       ----------------
C       IM4  - pointer to matrix operator
C       IFAC - sign of coefficient in linear combination
C       IFAS - sign of matrix operator
C         Anti-symmetric scalar operators are provided with an imaginary phase
C         in order to retain Hermiticity
C       IMAT - quaternion position of matrix operator
C
        KLBL = IPRPLBL(I,IOP)
        IH   = IPRLTYP(KLBL)
        IF(IABS(IH).NE.IHREF) THEN
           WRITE(6,*) 'PRPMAO ERROR: Mismatch square/triangular. IH=',IH
          CALL QUIT('PRPMAO: mismatch square/triangular')
        ENDIF
C
C       Determine what matrix the operator enters
C       =========================================
C
C       Is this a LL+SS operator (IC=1) or a LS+SL operator (IC=2) ?
        IF(PDOINT(KLBL)(2:3) .EQ. '00') THEN
          IC = 1
        ELSE IF((PDOINT(KLBL)(1:1) .EQ. '0') .AND.
     &          (PDOINT(KLBL)(4:4) .EQ. '0')) THEN
          IC = 2
        ELSE
          write(lupri,'(3A)') 'PRPMAO: illegal PDOINT = "',PDOINT,'"'
          CALL QUIT('PRPMAO: illegal PDOINT!')
        ENDIF
        IBRP = IPRLREP(KLBL)
        IM4  = JM4(I,IOPTYP)
        IZ   = JM4POS(IM4)
        IQ   = IQMULT(1,JQBAS(IBRP,IC),IZ)
        IQP  = IQTOPQ(IQ,IBRPOP)
! check with data printout
        IF(IQP.LE.0) THEN
          WRITE(LUPRI,'(/A)') '***** PRPMAO ERROR *****'
          WRITE(LUPRI,'(A,I3,A,I2,3X,A,A8)')
     &     'Component nr. ',I,'/',NCMP,' label:',PRPLBL(KLBL)
          WRITE(LUPRI,'(A,I2,A,I2,A,I2,A)')
     &    'IQP(',IQP,')=IQTOPQ(',IQ,',',IBRPOP,').LE.0 !'
          CALL WRIXPR(0,IOP)
          WRITE(LUPRI,'(A,I5)') 'Q: Scalar :',JQBAS(IBRP,IC)
          WRITE(LUPRI,'(A,I5)') 'Q: Matrix :',IZ
          WRITE(LUPRI,'(A,9I3)')
     &        'IM4,IBRP,IH,IZ,IQ,IBRPOP,IQP,IC,IOPTYP :',
     &         IM4,IBRP,IH,IZ,IQ,IBRPOP,IQP,IC,IOPTYP
          WRITE(LUPRI,'(A)') 'IQTOPQ(1:4,0:NBSYM-1):'
          DO IZ = 1,4
             WRITE(LUPRI,'(8I5)') (IQTOPQ(IZ,J),J=0,NBSYM-1)
          END DO
          WRITE(LUPRI,'(A)') 'IPQTOQ(1:NZ,0:NBSYM-1):'
          DO IZ = 1,NZ
             WRITE(LUPRI,'(8I5)') (IPQTOQ(IZ,J),J=0,NBSYM-1)
          END DO
          CALL QUIT('PRPMAO: IQP.LE.0 !!!')
        END IF
        IF(IQP.GT.NZ) THEN
          WRITE(LUPRI,'(/A)') '***** PRPMAO ERROR *****'
          WRITE(LUPRI,'(A,I5,3X,A8)') 'Component nr.',I,PRPLBL(KLBL)
          WRITE(LUPRI,'(A,I2,A,I2,A)') 'IQP(',IQP,').GT.NZ(',
     &          NZ,')'
          CALL WRIXPR(0,IOP)
          WRITE(LUPRI,'(A,I5)') 'Q: Scalar :',JQBAS(IBRP,IC)
          WRITE(LUPRI,'(A,I5)') 'Q: Matrix :',IZ
          CALL QUIT('PRPMAO: IQP.GT.NZ !!!')
        ENDIF
        IF(IPRINT.GE.5) THEN
          WRITE(LUPRI,'(A,I5,3X,A8)') 'Component nr.',I,PRPLBL(KLBL)
          FMT = MXFORM(PFAC,20)
          WRITE(LUPRI,'(A,4X,'//FMT//',A,A4)')
     &      'Factor: ',PFAC,' and PDOINT:',PDOINT(KLBL)
          WRITE(LUPRI,'(A,I5)') 'Integrals enter matrix :',IQP
          WRITE(LUPRI,'(A,L1,A,I2)')
     &    'FIRST = ',FIRST(IQP),' and IC(1 or 2):',IC
        ENDIF
C =============================================================
C               Get integrals from the CHECKPOINT file
C =============================================================
        IF(ONFILE) THEN
          CALL PRPREA(LU1INT,PRPLBL(KLBL),RTNLBL,PDOINT(KLBL),
     &          OP1INT,NDIM,IPRINT)
        ELSE
           IINT=KLBL-NPRPLBL
           CALL DCOPY(NDIM,PRPINT(1,IINT),1,OP1INT,1)
        ENDIF  
        IF(IPRINT.GE.10) THEN
          WRITE(LUPRI,'(A,I3)')
     &    '* Integrals from CHECKPOINT, size of matrix ',NTBAS(0)
          IF (NTBAS(0).EQ.NTBAS(1)) THEN
           WRITE(LUPRI,'(2X,A)')
     &   '...two-component framework (only LL block is present) '
          ELSE
           WRITE(LUPRI,'(2X,A)')
     &   '...four-component framework  '
          ENDIF
          IF (TRIANG) THEN
            CALL OUTPAK(OP1INT,NTBAS(0),-1,LUPRI)
          ELSE
            CALL OUTPUT(OP1INT,1,NTBAS(0),1,NTBAS(0),
     &                  NTBAS(0),NTBAS(0),-1,LUPRI)
          ENDIF
        ENDIF
C
C       Determine factor to scale matrix with
C       =====================================
C
        IFAC  = JCOM(I,IOPTYP)
        PFAC  = FACPRP(I,IOP)*IFAC
        IF(PFAC.NE.D1) THEN
          CALL DSCAL(NDIM,PFAC,OP1INT,1)
          IF(IPRINT.GE.15) THEN
            WRITE(LUPRI,'(/A,D13.6)')
     &    '* scaled integrals with the factor=',PFAC
           IF (TRIANG) THEN
             CALL OUTPAK(OP1INT,NTBAS(0),-1,LUPRI)
           ELSE
             CALL OUTPUT(OP1INT,1,NTBAS(0),1,NTBAS(0),
     &                   NTBAS(0),NTBAS(0),-1,LUPRI)
           ENDIF
          ENDIF
        ENDIF
C
C       Phase transformation/PDOINT factor if necessary.
C       The four blocks (AO basis): LL, SL, LS, SS
C
Chj     PDOINT is now '+', '-', or '0',
Chj     corresponding to a factor of +1, -1, or 0 on that block.
        NM1 = 0
        DO J = 1,4
           JCHR = PDOINT(KLBL)(J:J)
           IF (JCHR .EQ. '0') THEN
              IDOINT(J) = 0
              DOINT(J) = .FALSE.
           ELSE IF (JCHR .EQ. '+') THEN
              IDOINT(J) = 1
              DOINT(J) = .TRUE.
           ELSE IF (JCHR .EQ. '-') THEN
              NM1 = 1
              IDOINT(J) = -1
              DOINT(J) = .TRUE.
           ELSE
              CALL QUIT('PRPMAO: Illegal factor in PDOINT')
           END IF
        END DO
C
        IF (NZ.LT.4 .OR. NM1.GT.0) THEN
          IF (TRIANG) THEN
            CALL QNPHASE(IZ,IBRP,OP1INT,IDOINT)
          ELSE
            CALL QMPHASE('S',IZ,IBRP,OP1INT,DOINT)
          ENDIF

          IF(IPRINT.GE.15) THEN
             WRITE(LUPRI,'(A/3A)')
     &          '* Phase and/or PDOINT transformed integrals:',
     &          '  ( PDOINT = ',PDOINT(KLBL),' )'
             IF (TRIANG) THEN
               CALL OUTPAK(OP1INT,NTBAS(0),-1,LUPRI)
             ELSE
               CALL OUTPUT(OP1INT,1,NTBAS(0),1,NTBAS(0),
     &              NTBAS(0),NTBAS(0),-1,LUPRI)
             ENDIF
          ENDIF
        ENDIF
C
C       Determine where integrals go and add them on
C
        IF(FIRST(IQP)) THEN
          CALL DCOPY(NDIM,OP1INT,1,PMAT(1,IQP),1)
          FIRST(IQP) = .FALSE.
        ELSE
C         ... accumulate the property ...
          CALL DAXPY(NDIM,D1,OP1INT,1,PMAT(1,IQP),1)
        ENDIF
 10   CONTINUE

      IF(IPRINT.GE.10) THEN
        DO IZ = 1, NZ
          WRITE(LUPRI,'(A,I2,A)') '* Property matrix (IZ=',IZ,')'
          IF (FIRST(IZ)) THEN
             WRITE(LUPRI,'(/A)') '  --- zero matrix ---'
          ELSE
           IF (TRIANG) THEN
             CALL OUTPAK(PMAT(1,IZ),NTBAS(0),-1,LUPRI)
           ELSE
             CALL OUTPUT(PMAT(1,IZ),1,NTBAS(0),1,NTBAS(0),
     &                   NTBAS(0),NTBAS(0),-1,LUPRI)
           ENDIF
          END IF
        ENDDO
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRPMSAO(INDXPR,PRPAO,ONFILE,PRPINT,FIRST,
     &                   WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     PURPOSE: Get 2c/4c SA-AO matrix of one-electron operator
C     corresponding to property INDXPR of COMMON block XCBPRP
C
C     On input:
C    =============
C                INDXPR - index of the operator in the list
C                IDIMENS - dimension of arrays (triangular/square)
C                PBUF,OP1INT,FIRST - arrays to be used inside
C
C         PBUF, OP1INT have size of IDIMENS.
C
C    On output:
C   =============
C                PRPAO - 2c/4c property matrix in SA-AO basis
C                PBUF(IDIMENS,NZ) - property matrix in BU AO basis
C
C NOTE: In the TWOCOM.and.BSS=.true. case one needs only PRPAO array !
C       All other input parameters in this case can be dummy ...
C
C     Written by Miro ILIAS/Strasbourg, March 2006
C
C*****************************************************************************
!     use x2c_fio
      use quaternion_algebra

      use mxcent_mod
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0, PTOL = 1.0D-12)
C
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbprp.h"
#include "dcbxpr.h"
#include "dcbprl.h"
C
      LOGICAL TRIANGULAR,TOBE,YESPCTRA, MAKEPCTRA,FIRST(NZ),ONFILE
      DIMENSION PRPAO(N2BBASX,NZ),PRPINT(*),WORK(LWORK)
      CHARACTER PRPLABEL*8, prplab_x2c*12, NUMSTR*4
C
      CALL QENTER('PRPMSAO')
#include "memint.h"
C
C     Extract symmetry information - ISYM, IREP, ITIM
C     =============================
C
      IOPTYP     = IPRPTYP(INDXPR)
      ISYM       = IPRPSYM(INDXPR)
      IREP       = ISYM - 1
      ITIM       = IPRPTIM(INDXPR)
      IPRINT2    = MAX(IPRPRP,IPRINT)
      KLBL       = IPRPLBL(1,INDXPR)
      TRIANGULAR = ABS(IPRLTYP(KLBL)).EQ.1
      IF (IPRINT2.GE.10) THEN
       CALL HEADER('*** Output from PRPMSAO ***',-1)
       write(lupri,'(2X,A,I3,A,A)')
     & '* dealing with the operator with INDXPR=',
     & INDXPR,' name=',PRPNAM(INDXPR)
       write(lupri,'(2X,A,I1,A,I1,A,I2)')
     & '* operator ISYM=',ISYM,' IREP=',
     & IREP,' ITIM=',ITIM
       IF (TRIANGULAR) THEN
        write(lupri,'(2X,a)') '* this is a TRIANGULAR operator'
       ELSE
        write(lupri,'(2X,A)') '* this is a SQUARE operator'
       ENDIF
      ENDIF
C ========================================================================================
C    ... in this case read picture change transformed property integrals from BSSMAT file
C      which are there under the name  'P2C_'//NUMSTR
C
C     X2C module: picture change transformed property integrals reside
C                 on file X2CMAT
C ========================================================================================
      IF((TWOCOMP.AND.BSS).or.x2c) THEN
        IF(.NOT.ONFILE) CALL QUIT('PRPMSAO:2c not implemented yet !')
        IPRINT2 = MAX(IPRHAM,IPRINT)
        IF (N2BBASXQ.NE.(NTBAS(1)*NTBAS(1)*NZ)) THEN
          write(lupri,*) 'PRPMSAO error output:'
          write(lupri,*) 'N2BBASXQ=',N2BBASXQ
          write(lupri,*) 'NTBAS(1)=',NTBAS(1)
          write(lupri,*) 'NTBAS(1)*NTBAS(1)*NZ=',
     &                    NTBAS(1)*NTBAS(1)*NZ
          CALL QUIT('PRPMSAO: N2BBASXQ.ne.N2BBASXQ_L!!!')
        ENDIF
CMI    ...the label will contain the INDXPR operator index in the string form !
        CALL NUM2STR(INDXPR,NUMSTR)
        if(x2c)then
          write(prplab_x2c,'(a8,a4)') 'prpint2c',numstr
          open(lux2c,file='X2CMAT',status='old',form='unformatted',
     &        access='sequential',action="read",position='rewind')
          call x2c_read(prplab_x2c,prpao,n2bbasxq,lux2c)
          close(lux2c,status="keep")
        else
          PRPLABEL = 'P2C_'//NUMSTR
          CALL RFBSSMAT(PRPLABEL,.FALSE.,'PRPMSAO',
     &                 PRPAO,N2BBASXQ,LUBSS,IPRINT2)
          IF (IPRINT2.GE.6) THEN
            write(lupri,'(2x,A,A,A)')
     &      'PRPMSAO: 2comp. property record ',PRPLABEL,
     &      ' was read from BSSMAT file.'
          ENDIF
        end if ! if(x2c)then
        CALL LSET(NZ,.FALSE.,FIRST)
      ELSEIF(IOPTYP.EQ.0) THEN
C.......Fock matrix read
        KLBL  = IPRPLBL(1,INDXPR)
        INQUIRE(FILE=PRPLBL(KLBL),EXIST=TOBE)
        IF(TOBE) THEN
          CALL OPNFIL(2,PRPLBL(KLBL),'OLD','PRPMA1')
          CALL REAFCK(2,PRPAO,.FALSE.,1)
          CLOSE(2,STATUS = 'KEEP')
          PFAC = FACPRP(1,INDXPR)
          IF(PFAC.NE.D1) CALL DSCAL(N2BBASXQ,PFAC,PRPAO,1)
          CALL LSET(NZ,.FALSE.,FIRST)
        ELSE
          WRITE(LUPRI,'(A,A8,A,A16,A)')
     &      'Property file ',PRPLBL(KLBL),' for operator ',
     &      PRPNAM(INDXPR),' not found !'
          CALL QUIT('PRPMSAO: Property file not found !')
        ENDIF
      ELSE  ! IF((TWOCOMP.AND.BSS).or.x2c) THEN
CMI    -----  get the full four-component operator -------
        CALL LSET(NZ,.TRUE.,FIRST)
C     ... do resort PBUF to DIRAC sorting
        IPRINT2 = MAX(IPRPRP,IPRINT)
        IF (TRIANGULAR) THEN
          IDIMENS=NNBBASX
          CALL MEMGET('REAL',KPBUF  ,NZ*IDIMENS,WORK,KFREE,LFREE)
          CALL PRPMAO(LU1INT,INDXPR,ONFILE,PRPINT,TRIANGULAR,IDIMENS,
     &                WORK(KPBUF),PRPAO,FIRST,IPRINT2)
          CALL TRI2SQ(PRPAO,WORK(KPBUF),FIRST,IREP,ITIM)
          CALL MEMREL('PRPMSAO.tri2sq',WORK,KWORK,KPBUF,KFREE,LFREE)
        ELSE
          IDIMENS=N2BBASX
          CALL MEMGET('REAL',KONEIN,N2BBASX,WORK,KFREE,LFREE)
          CALL PRPMAO(LU1INT,INDXPR,ONFILE,PRPINT,TRIANGULAR,IDIMENS,
     &                PRPAO,WORK(KONEIN),FIRST,IPRINT2)
          CALL MEMREL('PRPMSAO.prpmao',WORK,KWORK,KONEIN,KFREE,LFREE)
          DO IZ = 1,NZ
            IF(FIRST(IZ)) CALL DZERO(PRPAO(1,IZ),N2BBASX)
          ENDDO
        ENDIF
        DO IZ = 1,NZ
          IF(.NOT.FIRST(IZ)) THEN
            CALL BUTOBS(PRPAO(1,IZ),1,WORK(KFREE),LFREE)
          ENDIF
        ENDDO
CMI   ----  Do the picture change transformation -----
        IF (BSS) THEN
         IF(.NOT.ONFILE) CALL QUIT('PRPMSAO: BSS not implemented yet !')
         IPRINT1 = MAX(IPRHAM,IPRPRP)
C     ... decide if the picture change transformation is needed
         CALL DECIDEPCTR(YESPCTRA,PRPNAM(INDXPR),IPRINT1)
C     ... if needed, do the picture change transformation of the PRPAO
         IF (YESPCTRA) THEN
          MAKEPCTRA = DOPCT .AND. .NOT.NOPCT
C  ... allocate array for the picture change transformation matrix
          ISIZE = NFBAS(1,0)*NFBAS(1,1)*NZ
          IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NFBAS(2,1)*NZ
          CALL MEMGET('REAL',KPCTM,ISIZE,WORK,KFREE,LFREE)
C     ... array for the 2-component picture change transformed operator
          CALL MEMGET('REAL',KPCTO,N2BBASXQ,WORK,KFREE,LFREE)

          CALL DOPCTRA(PRPAO,WORK(KPCTM),WORK(KPCTO),ISIZE,MAKEPCTRA,
     &                  PRPNAM(INDXPR),ISYM,IREP,ITIM,
     &                 'S',WORK(KFREE),LFREE,IPRINT1)
          CALL MEMREL('DOPCTRA',WORK,KWORK,KPCTM,KFREE,LFREE)
CMI       ... this is needed for the picture change transformed operator
          CALL LSET(NZ,.FALSE.,FIRST)
         ENDIF
        ENDIF
      ENDIF

      CALL QEXIT('PRPMSAO')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck decidepctr */
      SUBROUTINE DECIDEPCTR(DOTRANS,PROPNAME,IPRINT)
C*****************************************************************************
C
C  Purpose: Decides if the property operator (PROPNAME from the operator list)
C           has to be picture change transformed (returns logical DOTRANS).
C
C  Called from: PRPMSAO
C               PROP2BSS
C
C  Written by Miro ILIAS/Strasbourg/March 2006
C  Last modfications: Miro Ilias/Tel Aviv/November 2007
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"

      LOGICAL DOTRANS
      CHARACTER*16  PROPNAME, TPROPNAME

      CALL QENTER('DECIDEPCTR')

      DOTRANS = .TRUE.

      TPROPNAME = trim(PROPNAME)

C ... here is list of operators which DO not have to be picture change tranformed
C ... LL-forms of some operators (Beta matrix, Kinetic energy...) are not usable, however,
C ... some are usable (ZANGMOM)

      IF (TPROPNAME.EQ.'Overlap matrix'   .OR.  !
     &    TPROPNAME.EQ.'z-spin'           .OR.  ! overlap * 1/2
     &    TPROPNAME.EQ.'Nuc. attraction'  .OR.
     &    TPROPNAME.EQ.'Beta matrix'      .OR.
     &    TPROPNAME.EQ.'Kinetic energy'   .OR.
     &    TPROPNAME.EQ.'Orbital z-moment' .OR.  ! ... ZANGMOM
     &    TPROPNAME.EQ.'z-angular mom'    .OR.  ! ... ZANGMOM
     &    TPROPNAME.EQ.'XANGMOM'          .OR.  ! ... XANGMOM
     &    TPROPNAME.EQ.'YANGMOM'          .OR.  ! ... YANGMOM
     &    TPROPNAME.EQ.'ZANGMOM'          .OR.  ! ... ZANGMOM
!    &    TPROPNAME.EQ.'Spin-orbit matri' .OR.  ! ... spin-orbit matrix used for kappa matrix
     &    TPROPNAME.EQ.'d|S>_dBX'         .OR.
     &    TPROPNAME.EQ.'d|S>_dBY'         .OR.
     &    TPROPNAME.EQ.'d|S>_dBZ'         .OR.
     &    TPROPNAME.EQ.'dS_dBX'           .OR.
     &    TPROPNAME.EQ.'dS_dBY'           .OR.
     &    TPROPNAME.EQ.'dS_dBZ'           .OR.
     &    TPROPNAME.EQ.'dS_dB2XX'         .OR.
     &    TPROPNAME.EQ.'dS_dB2XY'         .OR.
     &    TPROPNAME.EQ.'dS_dB2XZ'         .OR.
     &    TPROPNAME.EQ.'dS_dB2YY'         .OR.
     &    TPROPNAME.EQ.'dS_dB2YZ'         .OR.
     &    TPROPNAME.EQ.'dS_dB2ZZ'
     &      ) DOTRANS = .FALSE.

      IF (IPRINT.GE.4) THEN
       IF (DOTRANS) THEN
         WRITE(lupri,'(/2X,3A)')
     &   'DECIDEPCTR: operator of entering name ',PROPNAME,
     &   ' HAS to be picture change transformed.'
       ELSE
         WRITE(lupri,'(/2X,3A)')
     &   'DECIDEPCTR: operator of entering name ',PROPNAME,
     &   ' DOES NOT have to be picture change transformed.'
       ENDIF
       CALL FLSHFO(LUPRI)
      ENDIF

      CALL QEXIT('DECIDEPCTR')
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tri2sq */
      SUBROUTINE TRI2SQ(PMAT,PTRI,FIRST,IREP,ITIM)
C*****************************************************************************
C
C    PURPOSE: Expand lower triangular matrix to squre form
C
C    Written by Miro ILIAS/Strasbourg/March 2006
C
C*****************************************************************************
#include "implicit.h"
#include "dgroup.h"
#include "dcbbas.h"
      DIMENSION PMAT(N2BBASX,NZ),PTRI(NNBBASX,NZ)
      LOGICAL   FIRST(NZ)
C
      CALL QENTER('SORT2SAO')

      DO IZ = 1,NZ
        IF(FIRST(IZ)) THEN
C ... zero out unused quaternionic blocks of property matrices...
          CALL DZERO(PMAT(1,IZ),N2BBASX)
        ELSE
          IQ = IPQTOQ(IZ,IREP)
          IH = IHQMAT(IQ,ITIM)
          IF(IH.EQ.1) THEN
            CALL DSPTSI(NTBAS(0),PTRI(1,IZ),PMAT(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DAPTGE(NTBAS(0),PTRI(1,IZ),PMAT(1,IZ))
          ENDIF
        ENDIF
      END DO
      CALL QEXIT('SORT2SAO')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lowdin */
      SUBROUTINE LOWDIN(VMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     LOWDIN and LOWDI1 routines
C     ===========================
C
C     PURPOSE:
C        Construction of the transformation matrix for Lowdin's
C        canonical orthonormalization
C        NB! Not the symmetrical orthonormalization, which
C            often has been called the "Lowdin orthonormalization"
C
C     METHOD:
C        In this subroutine V+ is computed. As the basis set is
C        real it reduces to the transpose:
C
C           V+ = TRN(V) = t(s**(-0.5))
C
C     Written by T.Saue
C     Last revision: Nov 23 1994 - tsaue
C     Last revision: Jan.2006 modified for pure nonrel.case...
C
C*****************************************************************************
      use primitives_mod
      use aovec_mod
      use maxaqn_mod
#include "implicit.h"
C
#include "dcbgen.h"
#include "dcbbas.h"
      LOGICAL LBIT
      DIMENSION VMAT(*),WORK(LWORK)
C
      CALL QENTER('LOWDIN')
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET('REAL',KSMAT,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,MXBBAS ,WORK,KFREE,LFREE)
C.....generate IRKB array if needed
      IF(.NOT.TWOCOMP.AND.LBIT(ISPHTR,2)) THEN
        CALL MEMGET('INTE',KIRKB,NPSHEL,WORK,KFREE,LFREE)
      ELSE
        KIRKB = KFREE
      ENDIF
C
C     Make Lowdin canonical matrix
C     ============================
C
      CALL LOWDI1(VMAT,WORK(KSMAT),WORK(KEIG),WORK(KIRKB),
     &            WORK(KFREE),LFREE,IPRINT)
C
C     Memory deallocation
C
      CALL MEMREL('LOWDIN.LOWDI1',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
      CALL QEXIT('LOWDIN')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lowdi1 */
      SUBROUTINE LOWDI1(VMAT,SMAT,EIG,IRKB,WORK,LWORK,IPRINT)
C*****************************************************************************
C     See SUBROUTINE LOWDIN for explanation
C*****************************************************************************
      use maxaqn_mod
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
      LOGICAL LBIT
      CHARACTER COMP(2)*1
      DIMENSION VMAT(*),SMAT(*),EIG(*),IRKB(*),WORK(LWORK)
C
      CALL QENTER('LOWDI1')
      IF(IPRINT.GE.5) CALL TITLER('Output from LOWDI1','*',103)
      COMP(1) = 'L'
      COMP(2) = 'S'
C
C     Get overlap matrix packed on fermion irreps
C     ===========================================
C
      CALL GTOVLT(SMAT,D1,IPRINT)
C      
C
C     Generate IRKB array for use in RKBLOW
C     =====================================
C
      IF(LBIT(ISPHTR,2).AND..NOT.(TWOCOMP.OR.MDIRAC)) THEN
        CALL RKBREL(IRKB,IPRINT)
      ENDIF
C
C     Initialize Lowdin canonical matrix
C     ========================
C
      CALL DZERO(VMAT,N2BBASX)
      NBRP = 4/NZ
      ISOFF = 1
      IVOFF = 1
      NLINDEP=0
      CALL HEADER('Generating Lowdin canonical matrix:',-1)
      DO 10 IFRP = 1,NFSYM
        IF(NFBAS(IFRP,0).EQ.0) GOTO 10
        JSOFF = ISOFF
        JVOFF = IVOFF
C
C       Large components
C       ----------------
C
        IC = 1
        NFORB(IFRP,IC) = 0
        IF(MDIRAC) THEN
          IP = IFRP
        ELSE
          IP = MOD(IFRP+IC,NFSYM) + 1
        ENDIF
        DO JSYM = 1,NBRP
          NBORB(JSYM,IFRP,IC) = 0
          IREP = JFSYM(JSYM,IP) - 1
          IF(NBBAS(IREP,IC).GT.0) THEN
            IF(LBIT(ISPHTR,IC)) THEN
              CALL SPHLOW(IREP,IC,SMAT(JSOFF),VMAT(JVOFF),
     &                    NFBAS(IFRP,0),NEFF,NLDP,STOL(IC),EIG,
     &                    WORK,LWORK,IPRINT)
              NDEL  = NBBAS(IREP,IC) - NEFF
              NLINDEP=NLINDEP+NLDP
            ELSE
              CALL LOWGEN(SMAT(JSOFF),NBBAS(IREP,IC),
     &                    VMAT(JVOFF),NEFF,NFBAS(IFRP,0),
     &                    STOL(IC),IPRINT,EIG,WORK,LWORK)
              NDEL  = NBBAS(IREP,IC) - NEFF
              NLINDEP=NLINDEP+NDEL              
              WRITE(LUPRI,'(3X,A1,3X,A3,3X,A,I10,A,E8.2)')
     &          COMP(IC),REP(IREP),'* Deleted: ',NDEL,
     &          ' *Smin: ',EIG(NBBAS(IREP,IC))
            ENDIF
C
C           Zero columns of Lowdin canonical matrix that are to be deleted
C
            IF(NDEL.GT.0) THEN
              J = JVOFF + NFBAS(IFRP,0)*NEFF
              DO I = 1,NDEL
                CALL DZERO(VMAT(J),NBBAS(IREP,IC))
                J = J + NFBAS(IFRP,0)
              ENDDO
            ENDIF
            JSOFF = JSOFF + (NFBAS(IFRP,0)+1)*NBBAS(IREP,IC)
            JVOFF = JVOFF + NFBAS(IFRP,0)*NEFF + NBBAS(IREP,IC)
            NBORB(JSYM,IFRP,IC) = NEFF
            NFORB(IFRP,IC) = NFORB(IFRP,IC) + NEFF
          ENDIF
        ENDDO
C
C       Small components
C       ----------------
C
        IC = 2
        NFORB(IFRP,IC) = 0
        IF(MDIRAC) THEN
          IP = IFRP
        ELSE
          IP = MOD(IFRP+IC,NFSYM) + 1
        ENDIF
        DO JSYM = 1,NBRP
          NBORB(JSYM,IFRP,IC) = 0
          IREP = JFSYM(JSYM,IP) - 1
          IF(NBBAS(IREP,IC).GT.0) THEN
            IF(LBIT(ISPHTR,IC)) THEN
              IF(MDIRAC) THEN
                CALL SPHLOW(IREP,IC,SMAT(JSOFF),VMAT(JVOFF),
     &                      NFBAS(IFRP,0),NEFF,NLDP,STOL(IC),EIG,
     &                      WORK,LWORK,IPRINT)
              ELSE
                CALL RKBLOW(IREP,IC,SMAT(JSOFF),VMAT(JVOFF),
     &                      NFBAS(IFRP,0),NEFF,NLDP,STOL(IC),EIG,
     &                      IRKB,WORK,LWORK,IPRINT)
              ENDIF
              NDEL    = NBBAS(IREP,IC) - NEFF
              NLINDEP = NLINDEP+NLDP              
            ELSE
              CALL LOWGEN(SMAT(JSOFF),NBBAS(IREP,IC),
     &                    VMAT(JVOFF),NEFF,NFBAS(IFRP,0),
     &                    STOL(IC),IPRINT,EIG,WORK,LWORK)
              NDEL  = NBBAS(IREP,IC) - NEFF
              NLINDEP=NLINDEP+NDEL
              WRITE(LUPRI,'(3X,A1,3X,A3,3X,A,I10,38X,A,E10.2)')
     &           COMP(IC),REP(IREP),'* Deleted: ',NDEL,
     &          ' *Smin: ',EIG(NBBAS(IREP,IC))
            ENDIF
C
C           Zero columns of Lowdin canonical matrix that are to be deleted
C
            IF(NDEL.GT.0) THEN
              J = JVOFF + NFBAS(IFRP,0)*NEFF
              DO I = 1,NDEL
                CALL DZERO(VMAT(J),NBBAS(IREP,IC))
                J = J + NFBAS(IFRP,0)
              ENDDO
            ENDIF
            NBORB(JSYM,IFRP,IC) = NEFF
            NFORB(IFRP,IC) = NFORB(IFRP,IC) + NEFF
            JSOFF = JSOFF + (NFBAS(IFRP,0)+1)*NBBAS(IREP,IC)
            JVOFF = JVOFF + NFBAS(IFRP,0)*NEFF+NBBAS(IREP,IC)
          ENDIF
          NBORB(JSYM,IFRP,0) = NBORB(JSYM,IFRP,1) + NBORB(JSYM,IFRP,2)
        ENDDO
        NFORB(IFRP,0) = NFORB(IFRP,1) + NFORB(IFRP,2)
        IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &    'LOWDI1: Lowdin canonical matrix (in ferm.symm.blocks) ',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** fermion corep ',IFRP,'/',NFSYM
           CALL PRQMAT(VMAT(IVOFF),NFBAS(IFRP,0),NFORB(IFRP,0),
     &           NFBAS(IFRP,0),NFORB(IFRP,0),1,IPQTOQ(1,0),LUPRI)
        END IF
        ISOFF = ISOFF + NFBAS(IFRP,0)*NFBAS(IFRP,0)
        IVOFF = IVOFF + NFBAS(IFRP,0)*NFORB(IFRP,0)
 10   CONTINUE
      IF(NLINDEP.GT.0) THEN
        WRITE(LUPRI,'(A,I0,A)') '*** WARNING *** : ',NLINDEP,
     &    ' functions deleted due to numerical linear dependence.'
      ENDIF

      CALL QEXIT('LOWDI1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgetint */
      SUBROUTINE RGETINT(CLASLB,OP1INT,DOINT,NCOMP,IORDER,NPQUAD,
     &                  MXFORM,TOFILE,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for getting one-electron property integrals
C     Written by T.Saue October 1994
C     Last revision : Nov 15 1994
C
C*****************************************************************************
      use primitives_mod
      use shells_mod
      use maxaqn_mod
      use mxcent_mod
      use symmet_mod
#include "implicit.h"
C
#include "nuclei.h"
      DIMENSION OP1INT(*),WORK(LWORK)
      CHARACTER CLASLB*7,MXFORM*6,LABINT(3*MXCOOR)*8
      LOGICAL DOINT(2,2),PROPRI,TRIANG,TOFILE
C
      CALL QENTER('RGETINT')
#include "memint.h"
C
      MAXTYP = 3*MXCOOR
      IF ((CLASLB.EQ.'ELFGRDC').OR.(CLASLB.EQ.'ELFGRDS')) THEN
         NINTAD = 9*NUCIND*(MAXREP+1)
      ELSE
         NINTAD = MAXTYP
      END IF
      PROPRI   = .FALSE.
      TRIANG   = .TRUE.
C
C     Memory allocation
C
      CALL MEMGET('INTE',KINTRP,MAXTYP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINTAD,NINTAD,WORK,KFREE,LFREE)
C
C     Calculate integrals in abacus part
C
      CALL PR1INT_1(WORK(KFREE),LFREE,WORK(KINTRP),WORK(KINTAD),
     &            LABINT,CLASLB,IORDER,NPQUAD,TRIANG,
     &            PROPRI,IPRINT,OP1INT,NCOMP,TOFILE,MXFORM,DOINT)
C
C     Memory deallocation
C
      CALL MEMREL('RGETINT.PR1INT_1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('RGETINT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vfsort */
      SUBROUTINE VFSORT(NLEFF,NSEFF,VMAT,VBUF,WORK,LWORK)
C*****************************************************************************
C
C     Sort Lowdin canonical matrix on fermion symmetries
C
C     Written by T.Saue Sep 12 1995
C     Last revision: Sep 12 1995
C*****************************************************************************
#include "implicit.h"
C Used from COMMON blocks:
C  DCBBAS: NTBAS(0)
C
#include "dcbbas.h"
      DIMENSION VMAT(*),VBUF(*),WORK(LWORK)
      CALL QENTER('VFSORT')
#include "memint.h"
C
C     Memory allocation
      CALL MEMGET('INTE',KSORT,NTBAS(0),WORK,KFREE,LFREE)
C
      CALL VFSOR1(NLEFF,NSEFF,VMAT,VBUF,WORK(KSORT))
C     Memory release
      CALL MEMREL('VFSORT.VFSOR1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('VFSORT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vfsor1 */
      SUBROUTINE VFSOR1(NLEFF,NSEFF,VMAT,VBUF,ISORT)
C*****************************************************************************
C
C     Sort Lowdin canonical matrix on fermion symmetries
C
C     Written by T.Saue Sep 12 1995
C     Last revision: Sep 12 1995
C*****************************************************************************
#include "implicit.h"
#include "dcbbas.h"
#include "dcborb.h"
      DIMENSION VBUF(NTBAS(0),*)
      DIMENSION VMAT(*),ISORT(*),NVEC(2)
C
C
C     Fermion symmetry 1
C     ==================
C
      CALL IZERO(NFORB,6)
      IOFF    = 1
      N2      = 0
      NVEC(1) = NLEFF
      NVEC(2) = NSEFF
      JJ = 0
      DO IC = 1,2
        DO J = 1,NVEC(IC)
          JJ = JJ + 1
          IND = IDAMAX(NTBAS(0),VBUF(1,JJ),1)
          IF(IND.LE.NFBAS(1,0)) THEN
            NFORB(1,IC) = NFORB(1,IC) + 1
            CALL DCOPY(NFBAS(1,0),VBUF(1,JJ),1,VMAT(IOFF),1)
            IOFF = IOFF + NFBAS(1,0)
          ELSE
            N2          = N2 + 1
            NFORB(2,IC) = NFORB(2,IC) + 1
            ISORT(N2)   = JJ
          ENDIF
        ENDDO
      ENDDO
C
C     Fermion  irrep 2
C     ================
C
      DO J = 1,N2
        CALL DCOPY(NFBAS(2,0),VBUF(NFBAS(1,0)+1,ISORT(J)),1,
     &             VMAT(IOFF),1)
        IOFF = IOFF + NFBAS(2,0)
      ENDDO
      NFORB(1,0) = NFORB(1,1)+NFORB(1,2)
      NFORB(2,0) = NFORB(2,1)+NFORB(2,2)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck freemt */
      SUBROUTINE FREEMT(FOCK,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for the construction of the one-electron free particle matrix
C     in SO-basis
C     Trond Saue, October 1995 - Oslo
C     Last revision: May 21 1996 - tsaue
C
C*****************************************************************************
#include "implicit.h"
C
      DIMENSION FOCK(*),WORK(LWORK)
#include "dcbbas.h"
#include "dgroup.h"
C
      CALL QENTER('FREEMT')
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET2('REAL','FTRI',KFTRI ,NZ*NNBBASX,WORK,KFREE,LFREE)
C
C     ****** Construct Fock matrix ******
C
      CALL FREEM1(FOCK,WORK(KFTRI),WORK(KFREE),LFREE,IPRINT)

C     Memory deallocation
      CALL MEMREL('FREEMT.FREEM1',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('FREEMT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck freem1 */
      SUBROUTINE FREEM1(FOCK,FTRI,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Driver for the construction of the one-electron free particle matrix
C     in SO-basis
C     Trond Saue, October 1995 - Oslo
C
C-----------------------------------------------------------------------------
C
C     The free particle matrix has the following form in the C1 case:
C
C             La      Sa      Lb      Sb
C     La  |   0    | -icDz |   0   | -icD- |
C     Sa  | -icDz  |  -C2  | -icD- |   0   |
C     Lb  |   0    | -icD+ |   0   | -icDz |
C     Sb  | -icD+  |   0   | -icDz |  -C2  |
C
C
C     The integrals:
C        S  - overlap integrals                  symmetric
C        C2 - (2c^2)S                            symmetric
C        D  - dipole velocity integrals          antisymmetric
C            D+ = Dx + iDy
C            D- = Dx - iDy
C     The ordering of components is chosen to reveal the following
C     matrix structure which is due to time reversal symmetry:
C
C     |  A  |  B  |
C     | -B* |  A* |
C
C     In this routine the full matrices A and B are constructed:
C
C            A(real)   |   A(imag)   |    B(real)   |  B(imag)
C     ------------------------------------------------------------
C           La      Sa     La    Sa  |  Lb     Sb     Lb     Sb
C     La |   0  |   0  |   0  | -cDz |    0  | -cDy |   0  | -cDx |
C     Sa |   0  |  -C2 | -cDz |   0  |  -cDy |   0  | -cDx |   0  |
C
C     It should be noted that Ar is symmetric, whereas Ai,Br and Bi
C     are antisymmetric.
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbgen.h"
C
      LOGICAL FIRST(4)
      DIMENSION FOCK(N2BBASX,NZ),FTRI(NNBBASX,NZ),
     &          WORK(LWORK)
      CALL QENTER('FREEM1')

C
C     Initialize
C
      FIRST(1:4) = .TRUE.
C
C     Beta matrix (FOCK is used as scratch space)
C     ===========
C
      IOP = IPBETAMT
      CALL PRPMAO(LU1INT,IOP,.TRUE.,WORK,.TRUE.,NNBBASX,FTRI,
     &            FOCK,FIRST,IPRHAM)
C
C     Kinetic energy (FOCK is used as scratch space)
C     ==============
C
      IOP = IPKINERG
      CALL PRPMAO(LU1INT,IOP,.TRUE.,WORK,.TRUE.,NNBBASX,FTRI,
     &            FOCK,FIRST,IPRHAM)
C
C     Make full matrix
C     ================
C
      DO IZ = 1,NZ
        IF(FIRST(IZ)) THEN
          CALL DZERO(FOCK(1,IZ),N2BBASX)
        ELSE
          IQ = IPQTOQ(IZ,0)
          IH = IHQMAT(IQ,1)
          IF(IH.EQ.1) THEN
            CALL DSPTSI(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DAPTGE(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ENDIF
        ENDIF
      ENDDO
      IF(IPRHAM.GE.5) THEN
        CALL HEADER(
     &   'FREEM1:Unsorted total one-electron Fock matrix',-1)
        CALL PRQMAT(FOCK,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &              IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Reindex to sorted basis
C     =======================
C
      CALL BUTOBS(FOCK,NZ,WORK,LWORK)
C
C     Output section
C     ==============
C
      IF(IPRHAM.GE.4) THEN
        CALL HEADER(
     &  'FREEM1: One-electron (symm.blocked) free particle matrix',-1)
        DO 10 I = 1,NFSYM
          IF(NFBAS(I,0).EQ.0) GOTO 10
          WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(FOCK(I2BASX(I,I)+1,1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   10   CONTINUE
      ENDIF
C
      CALL QEXIT('FREEM1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gmotra */
      SUBROUTINE GMOTRA(GOFAST)
C*****************************************************************************
C
C PURPOSE:
C     This matrix provides initial transformation matrix from
C     SO- to MO-basis. With no restrictions on kinetic balance
C     the Lowdin canonical matrix is used. With restricted balance, the
C     Lowdin canonical matrix is used to diagonalize the free-particle
C     matrix. Unphysical solutions are then deleted from the set
C     of vectors and the remaining vectors used as transformation matrix.
C
C MI,HJAaJ/2002-3
C     For the two-component relativistic calculations
C     apply the Barysz-Sadlej-Snijders transformation
C     on the Dirac one-electron Hamiltonian.
C     Save its elements and the picture change
C     transformation matrix as well.
C MI/2006
C     For the two-component nonrelativistic
C     calculations run the "URKBAL" branch as only
C     Lowdin canonical matrix is needed. It is assumed that small
C     components were eliminated in HERINP.
C
C     Called from: PSISCF (when BSS2DC or DC2BSS mode is activated)
C                  PAMSET (at the begining of the DIRAC run, and second time
C                          when pure two-comp. scheme+BSS is wished)
C
C     Written by T.Saue, October 1995
C     Last revision May 21 1996 - tsaue
C                    2006 - M.Ilias
C
C*****************************************************************************
      use memory_allocator
      use checkpoint
      use dircmo
      use x2cmod_cfg
      use dirx2c
      use x2c_utils, only:
     &    write_fock_matrices_saao_basis_x2c,
     &    read_fock_matrices_saao_basis_x2c

      use x2c_pct_mo_coefficients, only:
     &    x2c_molecular_mean_field_mos
      use x2c_cb_interface, only:
     &    renew_x2c_cb_orb_shell_dim
      use x2c_renew_AOMOMAT
      use quaternion_algebra
      use num_grid_gen
      use interface_mo_specific
      use interface_ao_specific
      use dirac_cfg, only: dirac_cfg_dft_calculation
#ifdef MOD_XAMFI
      use x_fock
      use xamfi_global_parameters
      use xamfi_internal_parameters
      use xamfi_environment
      use xamfi_utils, only:
     &    dump_normalized_2e_soc,
     &    read_xamfi_general_matrix,
     &    dump_delta_xc_energy,
     &    goodbye_xamfi,
     &    dump_Xc_c1_dens
      use xcint_main, only: get_xc_energy
      use aobasis_collector, only: deallocate_aobasis_container
#endif
!
      use maxaqn_mod
      use mxcent_mod
      use ccom_mod, only: DOSPHE

#include "implicit.h"
#ifndef MOD_XAMFI
      logical :: aoomod   = .false.
#endif
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
C Used from COMMON blocks:
C   DCBBAS: N2BBASX
C   PAMIOU: LUTMAT
#include "dcbgen.h"
#include "dcbcls.h"
#include "dcbham.h"
#include "cbihr1.h"
#include "cbirea.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
! needed for restart from X2CMAT
#include "dgroup.h"
#ifdef MOD_XAMFI
#include "nuclei.h"
#endif
!x2cmmf
#include "siripc.h"
      real(8), allocatable :: WORK(:)
      logical, intent(in)  :: GOFAST
!     local scratch
      logical              :: TOBE, FNDLAB, FNDLAB12, gmotra_case3
      logical              :: is_2c_rel_run, lsymad
      integer              :: is_defining_h1mat
      integer              :: ndim_lowmat,idum
      integer              :: lulmat = 99
      character(len=12)    :: flabel
      real(8), allocatable :: tmat_mmf(:)
      real(8), allocatable :: tbuf_mmf(:)
      real(8), allocatable :: cmo_mmf(:)
      real(8), allocatable :: scratch_mmf(:)
#ifdef MOD_XAMFI
      real(8), allocatable :: scr1_aoosoc(:)
      real(8), allocatable :: dmat_aoosoc(:)
      logical              :: no1esystem_aoosoc, aoc_save
      logical              :: x2c_save
      real*8               :: delta_xc_energy_aoosoc
#endif
      real(8)              :: timstr(2), timend(2), dum

      CALL QENTER('GMOTRA')

!     set local 2c-relativistic run flag
      is_2c_rel_run = bss.or.x2c
!     set file unit number for X2CMAT
      lux2c = 40

      IF(GOFAST)then
        if (.NOT.is_2c_rel_run) GOTO 10
        if (x2c) goto 7
      end if

      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in GMOTRA')
      KFRSAV = KFREE
C
!     logical flag for keeping track of whether case 3 in gmotra has
!     been entered
      gmotra_case3 = .false.

C     By default we do not partition the space. This may be
C     modified if we have linear symmetry or approximate the
C     Hamiltonian : Levy-Leblond, Spinfree, ..
C
      SUB_BL = .FALSE.
C
C     Make one-center projection matrix
C     =================================
C
      IF (ONECAP.AND.(INTV1C.EQ.1.AND.ICTLV1C(2).EQ.1)) THEN
        CALL MK1CEX (WORK,KFREE,LFREE,IPRGEN)
      ENDIF

!     we use N2BBASXQ instead of N2BBASX for BSS...to be checked
      CALL MEMGET2('REAL','VMAT',KVMAT,N2BBASXQ,WORK,KFREE,LFREE)
C
C     Get Lowdin canonical orthonormalization matrix
C     ==============================================

      inquire(exist=tobe,file='LOWDMAT')
!     delete existing LOWDIN matrix file when performing a geometry
!     optimization - it could be the one from a previous geometry.
!     stefan - june 2011
      if(tobe.and.optimi)then
        call opnfil(lulmat,'LOWDMAT','OLD','gmodel')
        close(lulmat,status='delete')
        tobe = .false.
      end if

      call opnfil(lulmat,'LOWDMAT','UNKNOWN','gmotra')
      rewind lulmat

      if(tobe)then
        CALL READI(LULMAT,1,NDIM_LOWMAT)
        CALL READT(LULMAT,NDIM_LOWMAT,WORK(KVMAT))
      else
        CALL TIMER2('START ',TIMSTR,TIMEND)
        CALL LOWDIN(WORK(KVMAT),WORK(KFREE),LFREE,IPRHAM)
        CALL TIMER2('Lwdn_a',TIMSTR,TIMEND)
!       save LOWDIN matrix on file
        NDIM_LOWMAT = 0
        DO I = 1, NFSYM
          NDIM_LOWMAT = NDIM_LOWMAT + NFORB(I,0)*NFBAS(I,0)
        END DO
        CALL WRITI(LULMAT,1,NDIM_LOWMAT)
        CALL WRITT(LULMAT,NDIM_LOWMAT,WORK(KVMAT))
      end if

      close(lulmat,status='keep')

!     release memory when we invoke the x2c module
      IF(X2C)THEN
        CALL MEMREL('X2CFRE',WORK,KWORK,KWORK,KFREE,LFREE)
        CALL DEALLOC(WORK)
      end if

#ifdef MOD_XAMFI
!     initialize atomic oo-order 2e-SOC module if enabled
      if(aoomod) call xamfi_init()
#endif

C     Unrestricted kinetic balance; use Lowdin canonical matrix
C     =========================================================
C
C     LV : The unrestricted procedure should be revised so that
C     the large components are also sorted on boson irreps !!!
C
      IF(.NOT.URKBAL.AND.NTBAS(2).EQ.0.AND.
     &   .NOT.(TWOCOMP.OR.is_2c_rel_run)) THEN
        URKBAL = .TRUE.
        WRITE(LUPRI,'(A)')
     &  '*** WARNING *** : No small bases. Setting URKBAL = T !'
      ENDIF
C
C     ** Branching point **
C=======================================================================
C
!       ---------------------------------------------------
!       Case 1: Unrestricted kinetic balance
!               or 2-component non-relativistic calculation
!       ---------------------------------------------------
 7    CONTINUE
      IF(URKBAL.OR.(TWOCOMP.AND..NOT.is_2c_rel_run))then
        NZT   = 1
        KTMAT = KVMAT
        DO IFRP = 1,NFSYM
          NESH(IFRP)  = NFORB(IFRP,1)
          NPSH(IFRP)  = NFORB(IFRP,2)
          NORB(IFRP)  = NESH(IFRP) + NPSH(IFRP)
          NTMO(IFRP)  = NORB(IFRP)
        ENDDO
        IF (TWOCOMP) THEN
          SSMTRC = D0
          MC     = 1
        ENDIF
        IF (IPRHAM.GE.5) THEN
          write(lupri,*) 'GMOTRA: branch URKBAL'
          write(lupri,*) 'NESH(IFRP),NPSH(IFRP)',
     &     (NESH(IFRP),NPSH(IFRP),IFRP=1,NFSYM)
        ENDIF
!       --------------------------------------------
!       Case 2: relativistic 2-component Hamiltonian
!       --------------------------------------------
      ELSE IF(is_2c_rel_run)THEN
!
        CALL TIMER2('START ',TIMSTR,TIMEND)

        TOBE=.FALSE.
        IF(GOFAST) THEN
          if(bss)then
            inquire(file='BSSMAT',exist=tobe)
          else if(x2c)then
            inquire(file='X2CMAT',exist=tobe)
          end if
        ENDIF

!       case 2.1 - step a: restart 2c-run from file
!       -------------------------------------------
        IF(TOBE) THEN
          if(bss)then
            CALL BSS_RESTART(WORK(KVMAT),IPRHAM,WORK,KFREE,LFREE)  
          else if(x2c)then
            CALL X2C_RESTART(IPRHAM)
            if(gofast) goto 10
          endif
!
!       case 2.2 - step a: generate 2c-Hamiltonian
!       ------------------------------------------
        ELSE
          if(bss)then
!           bss module (if asked for it can also perform an X2C transformation)
            CALL MAKE_H2C(WORK(KVMAT),WORK(KFREE),LFREE)
          else if(x2c)then
!           todo: put in an interface routine to be called inside x2c_driver
            if(mdirac)then
              NZT            = 1
              NPSH(1:nfsym)  = NFORB(1:nfsym,2) ! check!!!
            else
              NZT            = NZ
              NPSH(1:nfsym)  = NFORB(1:nfsym,1) ! check!!!
            end if
            DO IFRP = 1,NFSYM
              NESH(IFRP)  = NFORB(IFRP,1)
              NORB(IFRP)  = NESH(IFRP)    +  NPSH(IFRP)
              NTMO(IFRP)  = NORB(IFRP)
            END DO
!           optimized x2c module
            call x2c_main()
!           re-allocate
            call alloc(WORK,LWORK,id='WORK in GMOTRA - part 2')
#include "memint.h"

          end if

        ENDIF ! restart X2C or re-do X2C switch

        CALL TIMER2('mk_h2c',TIMSTR,TIMEND)

!       case 2 - step b: generate Lowdin canonical orthonormalization matrix in 2c-mode
!     -------------------------------------------------------------------------------
        if(.not.(twocompbss.or.x2c).or.x2cmod_mmf)then


          if(x2cmod_mmf)then
            CALL MEMGET2('REAL','TMAT',KTMAT,N2BBASXQ,WORK,KFREE,LFREE)
            call dzero(work(KTMAT),N2BBASXQ)
            open(lux2c,file='X2CMAT',status='old',form='unformatted',
     &           access='sequential',action='readwrite',
     &           position='rewind')
            ndim_lowdinmat   = 0
            do i = 1, nfsym
              ndim_lowdinmat = ndim_lowdinmat + nfbas(i,0)*nesh(i)
            end do

            if(linear)then
              gmotra_case3 = .true.
            else
              call x2c_read('2cLWAOg   10',work(KTMAT),
     &                      ndim_lowdinmat,lux2c)
            end if

            close(lux2c,status='keep')
            sub_bl = .false.
          else
!           WORK(KTMAT) carries already the large-electronic part of the Lowdin canonical matrix
            KTMAT = KVMAT
          end if
        else
C
          if(.not.x2c)
     &      CALL MEMREL('BSS.2c',WORK,KWORK,KWORK,KFREE,LFREE)

          TWOCOMP = .TRUE.
          RDINPC = .FALSE.
          IPREAD = IPREAD - 5
          CALL READIN(.FALSE.)
          !> stknecht: take care of 2c-mmf-sf for a linear molecule
          if(spinfr .and. linear) linear = .false.
          if(spinfr .and. atomic) atomic = .false.
          CALL SETDC1(IPREAD)
!         set AO-labels for large components - in x2c_main initialized for large and small components
          call getlab(ipread)
          IPREAD = IPREAD + 5
C.........reset NEWPRP to re-initialize AO matrices on CHECKPOINT
          NEWPRP = .TRUE.
          CALL LSET(NPRPCLS,.TRUE.,CLSCAL)
          CALL ONEGEN(WORK,LWORK)

          CALL MEMGET2('REAL','TMAT',KTMAT,N2BBASXQ,WORK,KFREE,LFREE)
C.........do canonical orthonormalization
          CALL TIMER2('START ',TIMSTR,TIMEND)
          CALL LOWDIN(WORK(KTMAT),WORK(KFREE),LFREE,IPRHAM)
          CALL TIMER2('Lwdn_b',TIMSTR,TIMEND)

!         save new 2c-LOWDIN matrix on file
          call opnfil(lulmat,'LOWDMAT','UNKNOWN','gmotra')
          rewind lulmat
          ndim_lowmat = 0
          do i = 1, nfsym
            ndim_lowmat = ndim_lowmat + nforb(i,0)*nfbas(i,0)
          end do
          call writi(lulmat,1,ndim_lowmat)
          call writt(lulmat,ndim_lowmat,work(kvmat))
          close(lulmat,status='keep')

!         reset orbital dimensions
          NZT = 1
          DO IFRP = 1,NFSYM
            NESH(IFRP)  = NFORB(IFRP,1)
            NPSH(IFRP)  = NFORB(IFRP,2)
            NORB(IFRP)  = NESH(IFRP) + NPSH(IFRP)
            NTMO(IFRP)  = NORB(IFRP)
          ENDDO
          SSMTRC = D0
          MC = 1

        end if

        IF (IPRHAM.GE.5) THEN
          write(lupri,*) 'GMOTRA: relativistic 2c-branch'
          write(lupri,*) 'NESH(IFRP),NPSH(IFRP)',
     &     (NESH(IFRP),NPSH(IFRP),IFRP=1,NFSYM)
        ENDIF
!
!     -------------------------------------------------------------------
!     Case 3: Transformation matrix with modified Dirac equation embedded
!     -------------------------------------------------------------------
      ELSE
        CALL MEMGET2('REAL','TMAT',KTMAT,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','EIG' ,KEIG ,NTBAS(0),WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','TBUF',KTBUF,N2BBASXQ,WORK,KFREE,LFREE)

!       defining h1 matrix control variable (in the X2C module) -
!       set here for consistency to the maximum value
        is_defining_h1mat = 4
        lsymad = .false.
        if(linear.or.atomic) lsymad = .true.

!        in: LOWDIN matrix                   --> KVMAT
!
!       out: AO2MO                           --> KTMAT
!       out: AO2MO SL resorted               --> KVMAT
!       out: MO2MO linear symmetry adapted   --> KTBUF (if linear == .true.; nothing otherwise)
        CALL MODHAM(WORK(KTMAT),WORK(KVMAT),WORK(KTBUF),WORK(KEIG),
     &              .not.(VEXTPJ.OR.FREEPJ.OR.LEVYLE),
     &              lsymad,
     &              is_defining_h1mat,.true.,WORK(KFREE),LFREE)

        CALL MEMREL('MODHAM',WORK,KEIG,KEIG,KFREE,LFREE)

        IF (IPRHAM.GE.5) THEN
          write(lupri,*) 'GMOTRA: branch MODHAM'
          write(lupri,*) 'NESH(IFRP),NPSH(IFRP)',
     &     (NESH(IFRP),NPSH(IFRP),IFRP=1,NFSYM)
        ENDIF

        gmotra_case3 = .true.

      ENDIF
C
C     ** End of branching **
C=======================================================================
C
      CALL MEMGET2('REAL','TBUF',KTBUF,N2BBASXQ,WORK,KFREE,LFREE)

      N2TMT   = 0
      N2TMOTQ = 0
      DO IFRP = 1,NFSYM
        NESHMF(IFRP) = NESH(IFRP)
        NPSHMF(IFRP) = NPSH(IFRP)
        NISHMF(IFRP) = NISH(IFRP)
        NOCCMF(IFRP) = NOCC(IFRP)
        I2TMT(IFRP)  = N2TMT
        N2TMT        = N2TMT + NFBAS(IFRP,0)*NORB(IFRP)*NZT
        I2TMOT(IFRP) = N2TMOTQ
        N2TMO(IFRP)  = NTMO(IFRP)*NTMO(IFRP)
        N2TMOTQ      = N2TMOTQ + N2TMO(IFRP)*NZ
      ENDDO

#ifdef MOD_XAMFI
      !> store common block variables while still being in 4c mode
      if(aoomod)then
        if (x2c) then
          call xamfi_set_cb("set-2c")
        else
          call xamfi_set_cb("set-4c")
        endif
      endif
#endif
C
C     Copy information from NBORB to the subblock array if we have not
C     partioned the orbital space otherwise.
C
      IF (.NOT.SUB_BL) CALL INISUB
      IF (SPINFR.or.(linear.and.x2c.and.x2cmod_mmf)) THEN
         SUB_BL = .TRUE.
      ENDIF
      if(linear) then
        if(gmotra_case3) then
          if(x2c.and.x2cmod_mmf)then
            open(lux2c,file='X2CMAT',status='old',form='unformatted',
     &           access='sequential',action='readwrite',
     &           position='rewind')
            call x2c_read('2cLWAOl   10',work(KTMAT),
     7                    ndim_lowdinmat,lux2c)
            close(lux2c,status='keep')
          else
            open(99,file='AOMOlin',status='old',form='unformatted',
     &          access='sequential',action='readwrite',
     &          position='rewind')

            !> rec # 1: U
            read(99)
            !> rec # 2: VU
            call readt(99,N2TMT,work(KTMAT))
!           file is needed in the 4c-2c transformation after the 4c-SCF
!           in order to ensure numerical stability (avoiding double
!           diagonalization, etc.)
            if(x2cmod_mmf)then
              close(99,status='keep')
            else
              close(99,status='delete')
            end if
          end if
        else 
C.........SUB_BL will be set true if linear symmetry works
          CALL DCOPY(N2TMT,WORK(KTMAT),1,WORK(KTBUF),1)
!         in:  KTMAT + KTBUF
!         out: KTMAT
          CALL MEMGET2('REAL','EIG'  ,KEIG  ,NTBAS(0),WORK,KFREE,LFREE)
          CALL MEMGET2('REAL','TBUF2',KTBUF2,N2BBASXQ,WORK,KFREE,LFREE)

!         in: AO2MO                         --> KTBUF
!
!         out: AO2MO linear symmetry adapted --> KTMAT
!         out: AO2MO                         --> KTBUF
!         out: MO2MO linear symmetry adapted --> KTBUF2 (discarded if not x2c-mmf)
          IF(ATOMIC)THEN
            CALL GATMSM(WORK(KTMAT),WORK(KTBUF),WORK(KTBUF2),WORK(KEIG),
     &                  WORK(KFREE),LFREE)
          ELSE
            CALL GLINSM(WORK(KTMAT),WORK(KTBUF),WORK(KTBUF2),WORK(KEIG),
     &                WORK(KFREE),LFREE)
          ENDIF
          CALL MEMREL('MODHAM',WORK,KEIG,KEIG,KFREE,LFREE)
        end if ! if(gmotra_case3) then
      end if       ! if(linear) then
C
C     Construct right index transformed overlap matrix Saomo = Sao V
C     Since V(dagger) Sao V = I, Saomo is the inverse of V(dagger)
C
!      in: AO2MO (linear symmetry adapted if linsym == true) --> KTMAT
!     out: right index transformed overlap matrix Saomo      --> KTBUF
      CALL MKSAOMO(WORK(KTBUF),WORK(KTMAT),.FALSE.,
     &             IPRHAM,WORK(KFREE),LFREE)
C
C     If requested, reduce transformation matrix by projection with
C     respect to a set of fragment orbitals
C
      IF(PROJEC) THEN
        CALL PSIPRJ(WORK(KTMAT),WORK(KTBUF),IPRONE,WORK,KFREE,LFREE)
        CALL SETDC2(0)
      ENDIF
C
C     Write transformation matrices
      CALL OPNFIL(LUTMAT,'AOMOMAT','UNKNOWN','GMOTRA')
!     AO2MO transformation matrix (linear symmetry adapted if linsym == .true.)
      CALL WRITT(LUTMAT,N2TMT,WORK(KTMAT))
!     right index transformed overlap matrix Saomo
      CALL WRITT(LUTMAT,N2TMT,WORK(KTBUF))
!     orbital dimensions, etc --> "runtime" data
      CALL ORBDAT

      CLOSE(LUTMAT,STATUS='KEEP')
#ifdef MOD_XAMFI
      if(.not.x2c.and.x2cmod_mmf.and.aoomod.and.nucdep==1
     &   .and.nucind==1)then
       CALL MEMGET2('REAL','SCR1'  ,KTMATaoo,N2TMT,WORK,KFREE,LFREE)
       CALL MEMGET2('REAL','SCR2'  ,KTBUFaoo,N2TMT,WORK,KFREE,LFREE)

       CALL OPNFIL(LUTMAT,'AOMOMAT','OLD','GMOTRA ')
       CALL READT(LUTMAT,N2TMT,work(ktmataoo))
       CALL READT(LUTMAT,N2TMT,work(ktbufaoo))
       CLOSE(LUTMAT,STATUS='KEEP')

       CALL OPNFIL(LUTMAT,'AOMOMATaoo','UNKNOWN','GMOTRA')
       CALL WRITT(LUTMAT,N2TMT,WORK(ktmataoo))
       CALL WRITT(LUTMAT,N2TMT,WORK(ktbufaoo))
       CALL ORBDAT
       CLOSE(LUTMAT,STATUS='KEEP')
      end if
#endif
      if(x2c.and.x2cmod_mmf)then
        if(aoomod)then
#ifdef MOD_XAMFI
          if(nucdep == 1 .and. nucind == 1)then
!           release memory
            CALL MEMREL('GMOTRA - aoosoc',WORK,KWORK,KWORK,KFREE,LFREE)
            call dealloc(work)
!           generate 2e-SO+scalar corrections in an atomic run
!>          ensure to not calculate any spurious 2e-SO+scalar corrections for actual 1e-systems
            call rmolchr(ichrg) ! sum of charge of all constituents, here only one atom
            no1esystem_aoosoc = .false.
            no1esystem_aoosoc = (ichrg - kcharg > 1)
            if(onesys) no1esystem_aoosoc = .false.
            call xamfi_atomcorr(delta_xc_energy_aoosoc,
     &                          no1esystem_aoosoc)
          end if ! check for atomic run
#else
        call quit('X-AMFI not available in this version')
#endif
        end if ! aoomod
        
        if(x2c_2c_mmf_mos)then
 
!         release memory
          if(allocated(work))then
            CALL MEMREL('GMOTRA - mmfmos',WORK,KWORK,KWORK,KFREE,LFREE)
            call dealloc(work)
          end if
!     
!         perform the final steps for the molecular mean-field procedure within the X2C-module
!         ------------------------------------------------------------------------------------

          !> step 1
          !> ****** 
!         !> purpose: transform the 2c-MO-coefficients from the orthonormal reference basis to AO basis
!                     and renew the coefficient file DFCOEF
!                     FIXME: a potential issue could be frozen orbitals - test me!
          allocate(tmat_mmf(n2bbasxq))
          allocate(tbuf_mmf(n2bbasxq))
          allocate(scratch_mmf(norbt))
          allocate(cmo_mmf(ncmotq))
          cmo_mmf     = 0
          scratch_mmf = 0
          tbuf_mmf    = 0
          tmat_mmf    = 0

          CALL OPNFIL(LUTMAT,'AOMOMAT','UNKNOWN','GMOTRA')
          CALL READT(LUTMAT,N2TMT,tmat_mmf)
          CLOSE(LUTMAT,STATUS='KEEP')

          call x2c_molecular_mean_field_mos(tmat_mmf,tbuf_mmf,cmo_mmf,
     &                                      scratch_mmf,
     &                                      i2tmt,i2tmot,icmoq,ntmo,
     &                                      nfbas(1,0),nfbas(1,1),nesh,
     &                                      nz,nzt,nfsym,
     &                                      nesht,ipqtoq,iprham)
          deallocate(cmo_mmf)
          deallocate(scratch_mmf)
          deallocate(tbuf_mmf)
          deallocate(tmat_mmf)

          !> step 2
          !> ******
          !> pupose: renew AOMOMAT with matrices of proper dimensions

          allocate(tmat_mmf(n2tmt),tbuf_mmf(n2tmt),scratch_mmf(n2tmt))
          tbuf_mmf = 0; tmat_mmf = 0; scratch_mmf = 0

          CALL OPNFIL(LUTMAT,'AOMOMAT','UNKNOWN','GMOTRA')
          CALL READT(LUTMAT,N2TMT,tmat_mmf)
          CALL READT(LUTMAT,N2TMT,tbuf_mmf)
          CLOSE(LUTMAT,STATUS='DELETE')

          CALL OPNFIL(LUTMAT,'AOMOMAT','UNKNOWN','GMOTRA')

          call x2c_molecular_mean_field_AOMOMAT(
     &            fh    = LUTMAT,
     &            tm1   = tmat_mmf,
     &            tm2   = tbuf_mmf,
     &            tmX   = scratch_mmf,
     &            lr1   = nfbas(1:nfsym,0),
     &            lc1   = norb(1:nfsym),
     &            lr2   = nfbas(1:nfsym,0),
     &            lc2   = norb(1:nfsym),
     &            lrX   = nfbas(1:nfsym,1),
     &            lcX   = nesh(1:nfsym),
     &            nfsym = nfsym,
     &            nzt   = nzt,
     &            I2TMT = i2tmt(1:nfsym),
     &            N2TMT = N2TMT,
     &            lupri = lupri,
     &            debug = iprham > 2
     &            )
          !> NOTE: ORBDAT is called after the reset of dimensions below and eventually AOMOMAT is closed
          deallocate(scratch_mmf,tbuf_mmf,tmat_mmf)

!         re-allocate
          call alloc(WORK,LWORK,id='WORK in GMOTRA - mmfmos')
#include "memint.h"

!         set intflag to (LL|LL) integrals only...
          call set_intflg(1)
          TWOCOMP = .TRUE.
          RDINPC = .FALSE.
          IPREAD = IPREAD - 5
          CALL READIN(.FALSE.)
          !> stknecht: take care of 2c-mmf-sf for a linear molecule
          if(spinfr .and. linear) linear = .false.
          if(spinfr .and. atomic) atomic = .false.
          CALL SETDC1(IPREAD)
!         set AO-labels for large components - in x2c_main initialized for large and small components
          call getlab(ipread)
          IPREAD = IPREAD + 5
!.........reset NEWPRP to re-initialize AO matrices on CHECKPOINT
          NEWPRP = .TRUE.
          CALL LSET(NPRPCLS,.TRUE.,CLSCAL)
          CALL ONEGEN(WORK(kfree),lfree)
          call renew_x2c_cb_orb_shell_dim(3)
C
C         Prepare for parallel calculation (set inforb)
C         ==============================================
C
          CALL RELINF
          NEWGEO = .TRUE.

!         reset 2e-common blocks in Hermit
          CALL PAOVEC(WORK(kfree),lfree,0,0)
!         reset DFDENS
          inquire(file='DFDENS',exist=tobe)
          if(tobe) then
            lublubb = 99
            call opnfil(lublubb,'DFDENS','OLD','GMOTRA')
            close(lublubb,status='DELETE')
          end if
!         delete and reset DFT grid
          call delete_num_grid(); call reset_num_grid()
!         reset interface-mo
          call interface_mo_write()
!         reset interface-ao
          call interface_ao_write()


          !> final step: store orbital dimensions, etc --> "runtime" data on AOMOMAT
          CALL ORBDAT
          CLOSE(LUTMAT,STATUS='KEEP')

!         release memory
          CALL MEMREL('GMOTRA - mmfmos-2',WORK,KWORK,KWORK,KFREE,LFREE)
          call dealloc(work)

!         !> post-SCF X2Cmmf:                   calculate H1_mmf_X2C_AO and save it on DFFCK1_mmf
!         !> atomic 2e-SO + scalar corrections: calculate G_A[2c] and save it on D2FCK_mmf

#ifdef MOD_XAMFI
          if(aoomod .and. aoofocc .and. aoc)then
            aoc_save   = aoc
            nfmat_save = nfmat
            !> implies focc ...
            aoc        = .false.
            nfmat      = 1
          end if
#endif

          len_fmat  = ntbas(0)**2*nz
          len_fmats = len_fmat*nfmat

          allocate(cmo_mmf(ncmotq))
          allocate(tbuf_mmf(len_fmats))
          cmo_mmf     = 0
          tbuf_mmf    = 0
          !> read 2c-mos
          call reacmo_new(cmo=cmo_mmf)

          !> calculate inactive (and possibly active) 2c density matrix
          call denmat(tbuf_mmf,cmo_mmf,iprham)

#ifdef MOD_XAMFI
          if(aooeamf)then
            !> save atomic density matrix (of atom(1) ... == atomic fragment)
            allocate(dmat_aoosoc(len_fmats)); dmat_aoosoc = 0
            call dcopy(len_fmats,tbuf_mmf,1,dmat_aoosoc,1)
            call dump_Xc_c1_dens(
     &                           dmat_aoosoc,
     &                           ntbas(0),
     &                           ntbas(0),
     &                           nz,
     &                           nfmat,
     &                           nopen,
     &                           df(0:nopen),
     &                           'aooDM_2c',
     &                           nint(charge(1)),
     &                           102,
     &                           .true.
     &                           )
            aooeamf= .false.; deallocate(dmat_aoosoc);
          endif
#endif

          allocate(tmat_mmf(len_fmats))
          tmat_mmf    = 0

#ifdef MOD_XAMFI
          if(.not.onesys)
     &    call get_scso_2e_fock(dmat                 = tbuf_mmf,
     &                          cmo                  = cmo_mmf,
     &                          fmat                 = tmat_mmf,
     &                          nrows                = ntbas(0),
     &                          ncols                = ntbas(0),
     &                          nr_2e_fock_matrices  = nfmat,
     &                          intflg               = intflg,
     &                          nz                   = nz,
     &                          ipqtoq               = ipqtoq,
     &                          x_nopen              = nopen,
     &                          x_dfopen             = df(0:nopen),
     &                          x_n2tmt              = n2tmt,
     &                          nzt_x                = nzt,
     &                          nr_fsym              = nfsym,
     &                          N2BBASXQ_dim_x       = n2bbasxq,
     &                          nr_tmo               = ntmo,
     &                          nr_ao_all            = nfbas(1:nfsym,0),
     &                          nr_mo_all            = norb(1:nfsym),
     &                          ioff_tmot            = i2tmot,
     &                          ioff_tmt             = i2tmt,
     &                          ioff_aomat_x         = i2basx,
     &                          myAOMOMAT            = 'AOMOMAT   ',
     &                          isDFT                =
     &                          dirac_cfg_dft_calculation,
     &                          isLoA                =
     &                          (linear.or.atomic),
     &                          xcmode               = 2,
     &                          lwork_ext            = lfree,
     &                          print_lvl            = iprham
     &                         )

          !> subtract XC_energy_2c
          if(dirac_cfg_dft_calculation)then
            delta_xc_energy_aoosoc =
     &      delta_xc_energy_aoosoc - get_xc_energy()
          end if
#endif

          !> here we split...
          !> a. either construct H1_mmf_X2C_AO for post-SCF X2Cmmf calculation ...
          if(.not.aoomod)then

            open(lux2c,file='X2CMAT',status='old',form='unformatted',
     &           access='sequential',action='readwrite',
     &           position='rewind')

!           2c-Fock operator, picture-change transformed from 4c-Fock operator ("++"-part)
            write(flabel,'(a7,i4,i1)') 'h12cAOn',1,0
            call x2c_read(flabel,tbuf_mmf,len_fmat,lux2c)
            call daxpy(len_fmat,-1.0d0,tmat_mmf,1,tbuf_mmf,1)
            close(lux2c, status='keep')

!           save H1_mmf_X2C_AO
            open(lux2c,file='DFFCK_mmf',status='unknown',
     &           access='sequential',action='readwrite',
     &           form='unformatted',position='rewind')
            call write_fock_matrices_saao_basis_x2c(
     &                                              tbuf_mmf,
     &                                              len_fmat, 1,'Fock',
     &                                              lux2c, 0)
            close(lux2c, status='keep')

          !> ... or
          !> b. save the G_A[2c] contribution to calculate the final 2e-SO+scalar atomic correction below
          else !
            open(lux2c,file='D2FCK_mmf',status='unknown',
     &           access='sequential',action='readwrite',
     &           form='unformatted',position='rewind')
            call write_fock_matrices_saao_basis_x2c(
     &                                              tmat_mmf,
     &                                              len_fmat, 1,'F[2]',
     &                                              lux2c, 0)
            close(lux2c, status='keep')

#ifdef MOD_XAMFI
            !> restore original settings
            if(aoofocc .and. aoc)then
              aoc        = aoc_save
              nfmat      = nfmat_save
            end if
#endif

          end if

          deallocate(tmat_mmf)
          deallocate(cmo_mmf)
          deallocate(tbuf_mmf)
#ifdef MOD_XAMFI
        !>                                                              ~
        !>  final step for atomic 2e-spin-orbit and scalar corrections: G_A^{++} = {U+ G_A[4c] U}^{++} - G_A[2c]

        if(aoomod .and. no1esystem_aoosoc)then
           if(nucdep == 1 .and. nucind == 1)then
            call xamfi_atomcorr2(delta_xc_energy_aoosoc,
     &                          no1esystem_aoosoc)              
          end if ! check for atomic run
        end if ! aoomod
#endif

        end if ! x2c_2c_mmf_mos == .true.

      end if ! x2c + x2cmod_mmf == .true.

!     final output
!     ------------
      CALL TITLER('Orbital dimensions','*',125)

      IF ( NFSYM .EQ. 1 ) THEN
         WRITE(LUPRI,'(A,I5)')
     &        'No. of positive energy orbitals (NESH): ',NESH(1)
         WRITE(LUPRI,'(A,I5)')
     &        'No. of negative energy orbitals (NPSH): ',NPSH(1)
         WRITE(LUPRI,'(A,I5)')
     &        'Total no. of orbitals           (NORB): ',NORB(1)
      ELSE
         WRITE(LUPRI,'(35X,5(A5,1X,I1,1X))')
     &        ('Irrep',I,I=1,NFSYM),'Sum'
         J = 0
         DO I = 1, NFSYM
            J = J + NESH(I)
         END DO
         WRITE(LUPRI,'(A,5(I5,3X))')
     &        'No. of electronic orbitals (NESH): ',
     &        (NESH(I),I=1,NFSYM),J
         J = 0
         DO I = 1, NFSYM
            J = J + NPSH(I)
         END DO
         WRITE(LUPRI,'(A,5(I5,3X))')
     &        'No. of positronic orbitals (NPSH): ',
     &        (NPSH(I),I=1,NFSYM),J
         J = 0
         DO I = 1, NFSYM
            J = J + NORB(I)
         END DO
         WRITE(LUPRI,'(A,5(I5,3X)/)')
     &        'Total no. of orbitals      (NORB): ',
     &        (NORB(I),I=1,NFSYM),J
      END IF
C
C     Memory deallocation
      if(allocated(work))then
        CALL MEMREL('GMOTRA',WORK,KWORK,KWORK,KFREE,LFREE)
        call dealloc(work)
      end if

 10   CONTINUE

#ifdef MOD_XAMFI
      if(aoomod) call xamfi_finalize()
#endif

      CALL QEXIT('GMOTRA')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck modham */
      SUBROUTINE MODHAM(TMAT,VMAT,TBUF,EIG,SLGROUP,LSYMAD,
     &                  is_defining_h1mat,print_header,WORK,LWORK)
C*****************************************************************************
C
C     Generate transformation matrix to MO-basis with the
C     quaternion modified Dirac equation embedded.
C     Ref: Visscher and Saue, JCP 113 (2000) 3996.
C
C     and for S-L resorting of columns (flag SLGROUP).
C
C     On input:
C   =============
C
C     VMAT - 4comp Lowdin transformation matrix for canonical orthonormalization
C            possibly with sperical and RKB transformations embedded.
C            (is destroyed!)
C
C     SLGROUP - if .true. (=Trond's 2comp.approach-one step BSS transformation)
C                resort columns of the transformation matrix into "S-L" order);
C               routine SLSORT does nothing if conditions for transformations are not met
C
C     LSYMAD  - if .true., apply the LINSYM procedure for adaptation
C               of the transformation matrix for the linear symmetry - needed in connection
C               to 4c->2c transition (.DO4C2C keyword)
C
C     is_defining_h1mat - control parameter for the defining h1 matrix
C                         in the X2C module wrt which the decoupling will be performed.
C
C     On output:
C   ==============
C
C  TMAT - total 4c AO -> MO transformation matrix
C
C  VMAT - if SLGROUP=.true., then it contains the standard (ALWAYS non-linear symmetry) transformation
C         matrix with SL-resorted columns
C
C  TBUF - if linear symmetry (LSYMAD) is on, it contains the transformation matrix from non-linear symmetry MOs
C         to the MOs of linear symmetry
C
C   Written by T.Saue, October 1995
C   Modified by L. Visscher, Jan 20, 1998 to include Levy-Leblond
C            M. Ilias, 2005, to include LSYMAD, SLGROUP
C            S. Knecht, 2010, to include is_defining_h1mat
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, DM2 = -2.0D0, D0 = 0.0D0)
C
C Used from COMMON blocks
C  CBIPAM:
C  DCBBAS:
C  DCBORB:
C  DGROUP:
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "cbihr1.h"
C
      DIMENSION TMAT(*),TBUF(*),VMAT(*),EIG(*),WORK(LWORK)
      LOGICAL, intent(in) :: SLGROUP,LSYMAD
      logical, intent(in) :: print_header
      integer is_defining_h1mat
      CALL QENTER('MODHAM')
#include "memint.h"
C
CMI  ... check of keywords..
      IF (SLGROUP.AND.(VEXTPJ.OR.FREEPJ.OR.LEVYLE))
     &  CALL QUIT('MODHAM: SLGROUP does not correspond with'//
     &  ' VEXTPJ/FREEPJ/LEVYLE flags!')
!
      FAC = DM2*CVAL*CVAL
      NSYM = 4/NZ
      if(.not.MDIRAC)THEN
C
C       Generate free particle (or bare nucleus) matrix in Lowdin-basis
C       ===============================================================
C
        IF(VEXTPJ) THEN
          CALL ONEFCK(TBUF,IPRONE,WORK(KFREE),LFREE)
        ELSE
          CALL FREEMT(TBUF,IPRONE,WORK(KFREE),LFREE)
        ENDIF
        IOFT = 1
        IOFV = 1
        DO 10 I = 1,NFSYM
          IF(NFORB(I,0).EQ.0) GOTO 10
CMI        ... transform to the Lowdin basis "Theta"
          CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                              NFORB(I,0),NFORB(I,0),
     &         TBUF(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &         TMAT(IOFT),NFORB(I,0),NFORB(I,0),NZ,IPQTOQ(1,0),
     &         VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &         VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &         WORK,LWORK,IPRGEN)
          IF (IPRHAM.GE.6) THEN
             IF(VEXTPJ) THEN
               CALL HEADER('MODHAM: '//
     &'Bare nucleus Dirac matrix in orthonormal canonical Lowdin basis'
     &,-1)
             ELSE
               CALL HEADER('MODHAM: '//
     &'Free particle Dirac matrix in orthonormal canonical Lowdin basis'
     &,-1)
             ENDIF
             WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
             CALL PRQMAT(TMAT(IOFT),NFORB(I,0),NFORB(I,0),NFORB(I,0),
     &                   NFORB(I,0),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
          IOFT = IOFT + NFORB(I,0)*NFORB(I,0)*NZ
          IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)
   10   CONTINUE
C
        NZT = NZ
      END IF

      IF(LEVYLE) THEN
C
C       Generate transformation matrix for the Levy-Leblond equation
C       ============================================================
C
        CALL LL_EQ(TMAT,TBUF,VMAT)

      ELSE
C
C       The modified Dirac equation at the matrix level (Hermit):
C       =========================================================

        CALL MOD_DIREQ(TMAT,TBUF,VMAT,EIG,WORK(KFREE),LFREE)
      ENDIF
C
C     Print section
C
      if(print_header)then
        CALL HEADER('Output from MODHAM',-1)
        IF(LEVYLE.OR.VEXTPJ.OR.FREEPJ) THEN
          WRITE(LUPRI,'(A)') ' * All positronic solutions deleted !'
        ELSE
          WRITE(LUPRI,'(A)') ' * Applied strict kinetic balance !'
        ENDIF
      end if
      IF(IPRHAM.GE.5) THEN
        CALL HEADER(
     &    'MODHAM: Actual transformation matrix 4c AO->MO ',-1)
        IMAT = 1
        DO 80 I = 1,NFSYM
        IF(NORB(I).GT.0) THEN
          WRITE(LUPRI,'(A,I1,A,I1)')
     &    '* Fermion ircop no.',I,'/',NFSYM
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
        ENDIF
 80     CONTINUE
      ENDIF

CMI ...  For the the BSS.RKB (ie 1step approach,X2c) resort columns of transformation matrix so that
CMI     they are in the order "S1 S2 ... S(nbsym) L1 L2 ... L(nbsym)" instead of default "S1 L1 S2 L2 ..."
!Miro:  Resorting is also applied for 4c DC_RKB Hamiltonians
      IF(SLGROUP)THEN
        ISIZE = NFBAS(1,0)*NORB(1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NORB(2)*NZ
!       is_defining_h1mat is the control variable for the defining h1
!       matrix wrt which the decoupling is performed in the X2C module.
!       see the following lsymad procedure for more comments - SK Aug 2010
        open(99,file='AOMOSLR',status='replace',form='unformatted',
     &       access='sequential',action='readwrite',
     &       position='rewind')
        if(is_defining_h1mat.gt.3.or.
     &    (is_defining_h1mat.eq.2.and.LSYMAD))then
!       if(is_defining_h1mat.gt.3)then
          CALL MEMGET2('REAL','BUF',KBUF,ISIZE,WORK,KFREE,LFREE)
!         in:  AO2MO transformation matrix             --> TMAT
!         out: AO2MO transformation matrix             --> TMAT
!              AO2MO transformation matrix SL resorted --> WORK(KBUF)
          CALL SLSORT(TMAT,WORK(KBUF),IPRHAM,WORK(KFREE),LFREE)
          call writt(99,ISIZE, WORK(KBUF))
          !> keep KBUF as we will need it later...
!CTROND
        else
          call dzero(VMAT,n2bbasxq)
!         in:  AO2MO transformation matrix             --> TMAT
!         out: AO2MO transformation matrix             --> TMAT
!              AO2MO transformation matrix SL resorted --> VMAT
          CALL SLSORT(TMAT,VMAT,IPRHAM,WORK(KFREE),LFREE)
          call writt(99,ISIZE, vmat)
!CTROND
        end if
        close(99,status='keep')
      ENDIF

!     adapt the MO part of the AO2MO transformation matrix wrt linear symmetry by
!     using the eigenvectors of <j_z> in MO basis as transformation matrix
!     ---------------------------------------------------------------------------
      open(99,file='AOMOlin',status='replace',form='unformatted',
     &     access='sequential',action='readwrite',
     &     position='rewind')
!     remember: in linear symmetry NZ ==> 1
      ISIZE  = NFBAS(1,0)*NORB(1)*NZ
      ISIZE2 = NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2)then
        ISIZE  = ISIZE  + NFBAS(2,0) * NORB(2)*NZ
        ISIZE2 = ISIZE2 + NORB(2)    * NORB(2)*NZ
      END IF

      IF(LSYMAD)THEN
        CALL MEMGET2('REAL','BUF1',KBUF1,ISIZE,WORK,KFREE,LFREE)
        CALL DCOPY(ISIZE,TMAT,1,VMAT,1)
        CALL DCOPY(ISIZE,TMAT,1,WORK(KBUF1),1)
        IF(IPRHAM.GE.5) THEN
          CALL HEADER(
     &    'MODHAM: BEFORE LINSYM final transformation'//
     &    ' matrix 4c AO->MO, VMAT(copyied from TMAT) ',-1)
          IMAT = 1
          DO I = 1,NFSYM
            IF(NORB(I).GT.0) THEN
              WRITE(LUPRI,'(/A,I1,A,I1)')
     &        '* Fermion ircop no.',I,'/',NFSYM
              CALL PRQMAT(VMAT(IMAT),NFBAS(I,0),NORB(I),
     &                    NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
              IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
            ENDIF
          ENDDO
        ENDIF

        CALL DZERO(TMAT,N2BBASXQ)
        CALL DZERO(TBUF,N2BBASXQ)
        IF (.NOT.SUB_BL) CALL INISUB
        IF (SPINFR) THEN
          SUB_BL = .TRUE.
        ENDIF

!       in:  VMAT

        IF(ATOMIC)THEN
!       out: AO2MO atomic symmetry adapted (kappa/mj) --> TMAT
!       out: AO2MO atomic symmetry adapted (kappa   ) --> VMAT
!       out: MO2MO atomic symmetry adapted WRONG .... --> TBUF (discarded if not x2c-mmf)
          CALL GATMSM(TMAT,VMAT,TBUF,EIG,WORK(KFREE),LFREE)
        ELSE
!       out: AO2MO linear symmetry adapted --> TMAT
!       out: AO2MO                         --> VMAT
!       out: MO2MO linear symmetry adapted --> TBUF (discarded if not x2c-mmf)
          CALL GLINSM(TMAT,VMAT,TBUF,EIG,WORK(KFREE),LFREE)
        ENDIF
!       save matrices on file to be re-used in the x2c-module
        if(.not.spinfr)then
!         MO2MO lin
          call writt(99,ISIZE2,TBUF)
!         AO2MO lin
          call writt(99,ISIZE, TMAT)
        end if
!       AO2MO (on backup buffer space)
        call writt(99,ISIZE, WORK(KBUF1))

!       by default: return the unmodified AO2MO transformation matrix
        CALL DCOPY(ISIZE,WORK(KBUF1),1,TMAT,1)

        IF(IPRHAM.GE.5) THEN
          IMAT = 1
          IBUF = 1
          DO I = 1,NFSYM
            IF(NORB(I).GT.0) THEN
                  CALL HEADER(
     &    'MODHAM: AFTER LINSYM final transformation'//
     &    ' matrix 4c AO->MO, TMAT(output)',-1)
              WRITE(LUPRI,'(A,I1,A,I1)')
     &        '* Fermion ircop no.',I,'/',NFSYM
              CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                    NFBAS(I,0),NORB(I),1,IPQTOQ(1,0),LUPRI)
              WRITE(LUPRI,'(/5X,A/)')
     &    'MODHAM: 4c ON  RKB "theta" -> LSYM "theta" transformation'//
     &        '  matrix (TBUF).'
              WRITE(LUPRI,'(A,I1,A,I1)')
     &        '* Fermion ircop no.',I,'/',NFSYM
              CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),
     &                    NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
              IMAT = IMAT + NFBAS(I,0)*NORB(I)*1
              IBUF = IBUF + NORB(I)*NORB(I)*1
            ENDIF
          ENDDO
        ENDIF
        CALL MEMREL('SLGROU1',WORK,KBUF1,KBUF1,KFREE,LFREE)
      else
!       AO2MO transformation matrix
        call writt(99,ISIZE, TMAT)
      endif
      close(99,status='keep')

      IF(SLGROUP.and.(is_defining_h1mat.gt.3.or.
     &    (is_defining_h1mat.eq.2.and.LSYMAD)))THEN
C  ...  copy the SL resorted transformation matrix to be exported into VMAT ...
        CALL DCOPY(ISIZE,WORK(KBUF),1,VMAT,1)
        IF(IPRHAM.GE.5) THEN
          CALL HEADER(
     &    'MODHAM: "SL" resorted transformation matrix'//
     &    ' 4c AO->MO, VMAT(ouput) ',-1)
          IBUF = 1
          DO I = 1, NFSYM
            IF(NORB(I).GT.0) THEN
              WRITE(LUPRI,'(A,I1,A,I1)')
     &        '* Fermion ircop no.',I,'/',NFSYM
              CALL PRQMAT(VMAT(IBUF),NFBAS(I,0),NORB(I),
     &                    NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
              IBUF = IBUF + (NFBAS(I,0)*NORB(I)*NZ)
            ENDIF
          ENDDO
        ENDIF
        CALL MEMREL('SLGROU2',WORK,KBUF,KBUF,KFREE,LFREE)
      ENDIF

      CALL FLSHFO(LUPRI)
      CALL QEXIT('MODHAM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prprea */
      SUBROUTINE PRPREA(LU1INT,PRPLBL,RTNLBL,PCOMB,PRPINT,NDIM,IPRINT)
      use labeled_storage
C***********************************************************************
C
C     Read property integrals from file with the label PRPLBL
C     and class blocks defined by COMB
C     Return integrals in PRPINT and return labels in RTNLBL
C
C     Called from: PRPMAO - general routine to get AO integrals from file
C                  GTOVLX - to get overlap integrals
C                  GTOVLT - to get symmetry packed overlap integrals
C
C     Written by T.Saue, January 1995
C     Last revision: May 19 1996 - tsaue
C                    MI/March 2006 - more print-out
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      LOGICAL FNDLB3, FINDLABEL
      CHARACTER PRPLBL*8, RTNLBL(2)*8, PCOMB*4, COMB*4
      DIMENSION PRPINT(NDIM)
      character*(:), allocatable :: label
      logical                    :: label_present
      CALL QENTER('PRPREA')
C
Chj   PDOINT is now '+', '-', or '0',
Chj   corresponding to a factor of +1, -1, or 0 on that block.
C     The four blocks (AO basis): LL, SL, LS, SS
      DO J = 1,4
        IF (PCOMB(J:J) .EQ. '0' .OR. PCOMB(J:J) .EQ. 'F') THEN
           COMB(J:J) = 'F'
        ELSE
           COMB(J:J) = 'T'
        END IF
      END DO
C
      IF (IPRINT.GE.12) THEN
       WRITE(LUPRI,'(A,I4)')
     & 'PRPREA: >>> Going to read the property '//PRPLBL//
     & ' with COMB '//COMB//
     & ' from the file with LU1INT = ',LU1INT
      ENDIF
      label = '/result/operators/ao_matrices/'//
     &        PRPLBL//COMB
      call lab_query(dirac_datafiles(1),label,label_present)
      if (label_present) then
        call lab_read(dirac_datafiles(1),label,rdata=prpint(1:ndim))
        IF (IPRINT.GE.17) THEN
           WRITE(LUPRI,'(A,I4)')
     &    'PRPREA: >>> Property '//PRPLBL//' '//COMB//
     &    'succesfully read from '//dirac_datafiles(1)%name
        ENDIF
      ELSE
        WRITE(LUPRI,'(A,I4)')
     &    'PRPREA: >>> Property '//PRPLBL//' '//COMB//
     &                     ' not found on '//dirac_datafiles(1)%name
        CALL QUIT('PRPREA:Property '//PRPLBL//
     &   ' of COMB '//COMB//' not found')
      ENDIF
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('PRPREA')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fndlb3 */
      LOGICAL FUNCTION FNDLB3(SRCLBL,COMB,RTNLBL,LU,IPRINT)
C***********************************************************************
C
C  18-Jan-1995 tsaue
C  1. Find MOLECULE label and block COMB on unformatted file LU
C       a)Cannot find label: return FALSE
C       b)Finds label      : continues

CMI: 12/12/2008: added printout for betted debugg
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      LOGICAL FNDLB2
      CHARACTER SRCLBL*8, RTNLBL(2)*8, COMB*4
      FNDLB3 = .FALSE.
   10 CONTINUE
      IF(FNDLB2(SRCLBL,RTNLBL,LU)) THEN
C     Found label SRCLBL on LU; check COMB
        IF(RTNLBL(2)(5:8).NE.COMB) THEN
CMI       ... control print out ...
          IF (IPRINT.GE.30) THEN
            WRITE(LUPRI,'(2X,A,A)') 'FNDLB3: searching for SRCLBL=',
     &      SRCLBL
            WRITE(LUPRI,'(2X,A,A4,A,A4)') 'FNDLB3: '//
     &      'RTNLBL(2)(5:8)=',RTNLBL(2)(5:8),'.ne. COMB=',COMB
          ENDIF
          GOTO 10
        ENDIF
        FNDLB3 = .TRUE.
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GTOVLX(SMAT,SSMETR)
C*****************************************************************************
C
C     Get overlap matrix from file CHECKPOINT
C     The SS block is scaled by the factor SSMETR
C     No symmetry packing
C
C     Written by T.Saue May 26 1996
C     Last revision: Jan 20 1998 - L. Visscher
C
C*****************************************************************************
      use memory_allocator
#include "implicit.h"
      PARAMETER(D1 = 1.0D0)
C
#include "dcbham.h"
#include "dcbgen.h"
#include "cbihr1.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION SMAT(N2BBASX)
      CHARACTER*8 RTNLBL(2)
      real(8), allocatable :: sbuf(:)
      CALL QENTER('GTOVLX')
C
C
C       Allocate integral buffer
C       ========================
C
        call alloc(sbuf, NNBBASX)
C
C       Read integral from file
C       =======================
C
        CALL PRPREA(LU1INT,'OVERLAP ',RTNLBL,'TFFT',sbuf,NNBBASX,
     &              IPRONE)
C
C       Make full matrix
C       ================
C
        CALL DSPTSI(NTBAS(0),SBUF,SMAT)
C
C       Memory deallocation
C       ===================
C
        call dealloc(sbuf)
C
C       Reindex to sorted basis
C       =======================
C
        CALL BUTOBS_no_work(SMAT,1)
C
C       Scale SS-block
C       ==============
C
        IF(SSMETR.NE.D1) THEN
          DO I = 1,NFSYM
             JJ = I2BASX(I,I)+(NFBAS(I,1)*NTBAS(0))+NFBAS(I,1)+1
             DO J = 1, NFBAS(I,2)
                CALL DSCAL (NFBAS(I,2),SSMETR,SMAT(JJ),1)
                JJ = JJ + NTBAS(0)
             ENDDO
          ENDDO
        ENDIF
C
      CALL QEXIT('GTOVLX')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gtovlt */
      SUBROUTINE GTOVLT(SMAT,SSMETR,IPRINT)
C*****************************************************************************
C
C     Get overlap matrix from file CHECKPOINT
C     The SS block is scaled by the factor SSMETR
C     Packed on fermion irreps
C
C     Written by T.Saue May 26 1996
C     Last revision: Jan 20 1998 - L. Visscher
C
CMI/Sept2005 The SMAT must be allocated at least to the N2BAST size !
C
C*****************************************************************************
      use primitives_mod
      use maxaqn_mod
      use mxcent_mod
      use symmet_mod
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "cbihr1.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION SMAT(*)
      CHARACTER*8 RTNLBL(2)
      real(8), allocatable :: stri(:)
      real(8), allocatable :: sbuf(:)
      CALL QENTER('GTOVLT')
C
C     Allocate integral buffer
C     ========================
C
      CALL DZERO(SMAT,N2BAST)
      NBRP = 4/NZ

        allocate(sbuf(n2bbasx))
        allocate(stri(nnbbasx))
C
C       Read integral from file
C       =======================
C
        CALL PRPREA(LU1INT,'OVERLAP ',RTNLBL,'TFFT',STRI,NNBBASX,
     &              IPRONE)
C
C       Make full matrix
C       ================
C
        CALL DSPTSI(NTBAS(0),STRI,SBUF)
C
        CALL DZERO(SMAT,N2BAST)
        NBRP     = 4/NZ
        IOFF = 0
        IF (TWOCOMP) THEN
           IC_MAX = 1
        ELSE
           IC_MAX = 2
        ENDIF
        DO I = 1,NFSYM
          DO IC = 1,IC_MAX
            IP    = MOD(I+IC,NFSYM) + 1
            DO JSYM = 1,NBRP
              ISYM           = JFSYM(JSYM,IP)
              IREP           = ISYM - 1
              IF(NBBAS(IREP,IC).GT.0) THEN
                ISOFF = IOFF + (IBBAS(IREP,IC)-IBAS(I))*(NFBAS(I,0)+1)
                IUOFF = ICOS(ISYM,IC)*(NTBAS(0)+1)
                DO J = 1,NBBAS(IREP,IC)
                  CALL DCOPY(NBBAS(IREP,IC),SBUF(1+IUOFF),1,
     &                                      SMAT(1+ISOFF)    ,1)
                  IUOFF = IUOFF + NTBAS(0)
                  ISOFF = ISOFF + NFBAS(I,0)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
          IOFF = IOFF + N2BAS(I)
        ENDDO

        deallocate(stri)
        deallocate(sbuf)
C
C     Scale SS-block
C     ==============
C
      IF(SSMETR.NE.D1 .AND. .NOT.TWOCOMP) THEN
        IOFF = 1
        DO I = 1,NFSYM
          JJ = IOFF + NFBAS(I,1)*(NFBAS(I,0)+1)
          DO J = 1, NFBAS(I,2)
            CALL DSCAL (NFBAS(I,2),SSMETR,SMAT(JJ),1)
            JJ = JJ + NFBAS(I,0)
          ENDDO
          IOFF = IOFF + N2BAS(I)
        ENDDO
      ENDIF
C
C     Print section
C     =============
C
      IF(IPRINT.GE.5) THEN
        CALL HEADER('GTOVLT: Overlap (metric) matrix:',-1)
        IOFF = 1
        DO I = 1,NFSYM
        IF(NFBAS(I,0).GT.0) THEN
          WRITE(LUPRI,'(A,I1,A,I1)')
     &     '* Fermion ircop no.',I,'/',NFSYM
          CALL PRQMAT(SMAT(IOFF),NFBAS(I,0),NFBAS(I,0),
     &                NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),LUPRI)
          IOFF = IOFF + N2BAS(I)
        ENDIF
        ENDDO
      ENDIF
C
      CALL QEXIT('GTOVLT')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reno */
      SUBROUTINE RENO(LDM,NDIM,SMAT)
#include "implicit.h"
      DIMENSION SMAT(LDM,NDIM)
      DO J = 1,NDIM
        DO I = 1,NDIM
          SMAT(I,J) = SMAT(I,J)/SQRT(SMAT(I,I)*SMAT(J,J))
        ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck sphlow */
      SUBROUTINE SPHLOW(IREP,IC,SMAT,VMAT,LDM,NEFF,NLDP,STOL,
     &                  EIG,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate Lowdin's canoncial orthonormalization matrix with
C     cartesian to spherical transformation embedded.
C
C***********************************************************************
#include "implicit.h"
C
C
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),
     &          WORK(*)
C
      CALL QENTER('SPHLOW')
#include "memint.h"
      NCAR = NBBAS(IREP,IC)
      NMAT = NCAR*NCAR
      CALL MEMGET2('REAL','TMAT',KTMAT,NMAT,WORK,KFREE,LFREE)
      CALL SPHLO1(IREP,IC,NCAR,SMAT,VMAT,LDM,NEFF,NLDP,STOL,EIG,
     &            WORK(KTMAT),WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('SPHLOW',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('SPHLOW')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck sphlo1 */
      SUBROUTINE SPHLO1(IREP,IC,NCAR,SMAT,VMAT,LDM,NEFF,NLDP,STOL,EIG,
     &                  TMAT,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate Lowdin's canoncial orthonormalization matrix with
C     cartesian to spherical transformation embedded.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
#include "dgroup.h"
#include "pgroup.h"
#include "dcbham.h"
      LOGICAL DONRM
      CHARACTER COMP(2)*1
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),TMAT(NCAR,NCAR),WORK(*)
      DATA COMP/'L','S'/
C
C**************************************************
C***  Generate spherical transformation matrix  ***
C**************************************************
C
      DONRM = .TRUE.
      CALL SPHCAR(TMAT,IREP,IC,NCAR,NSPH,
     &            DONRM,IPRINT,WORK,LWORK)
C
C*********************************************************
C***  Transform overlap matrix to spherical basis and
C***  generate Lowdin canonical orthonromalization matrix
C*********************************************************
C
      CALL QTRANS('AOMO','S',D0,NCAR,NCAR,NSPH,NSPH,
     &            SMAT,LDM ,NCAR,1,IPQTOQ(1,0),
     &            VMAT,LDM ,NSPH,1,IPQTOQ(1,0),
     &            TMAT,NCAR,NSPH,1,IPQTOQ(1,0),
     &            TMAT,NCAR,NSPH,1,IPQTOQ(1,0),
     &            WORK,LWORK,IPRINT)
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A,1X,A3,2X,A1)')
     &    'SPHLOW: Overlap matrix in solid harmonics for',
     &    REP(IREP),COMP(IC)
        CALL OUTPUT(VMAT,1,NSPH,1,NSPH,LDM,NSPH,-1,LUPRI)
      ENDIF
      CALL LOWGEN(VMAT,NSPH,SMAT,NEFF,LDM,STOL,IPRINT,
     &            EIG,WORK,LWORK)
      CALL DGEMM('N','N',NCAR,NEFF,NSPH,D1,TMAT,NCAR,
     &           SMAT,LDM,D0,VMAT,LDM)
      NA   = NCAR - NSPH
      NLDP = NSPH - NEFF
      ND   = NA + NLDP
      WRITE(LUPRI,'(3X,A1,3X,A3,3X,3(A,I10),A,E8.2)')
     &           COMP(IC),REP(IREP),'* Deleted: ',ND,
     &           '(Proj: ',NA,', Lindep: ',NLDP,') Smin: ',
     &           EIG(NSPH)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck rkblow */
      SUBROUTINE RKBLOW(IREP,IC,SMAT,VMAT,LDM,NEFF,NLDP,STOL,
     &                  EIG,IRKB,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate Lowdin's canoncial orthonormalization matrix with
C     the restricted kinetic balance (RKB) transformation embedded.
C     See SUBROUTINE RKBLO1 for further explanations.
C
C***********************************************************************
#include "implicit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),IRKB(*),WORK(*)
C
      CALL QENTER('RKBLOW')
#include "memint.h"
      NCAR = NBBAS(IREP,IC)
      NMAT = NCAR*NCAR
      CALL MEMGET2('REAL','TMAT',KTMAT,NMAT,WORK,KFREE,LFREE)
      CALL RKBLO1(IREP,IC,NCAR,SMAT,VMAT,LDM,NEFF,NLDP,STOL,EIG,
     &            WORK(KTMAT),IRKB,WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('RKBLOW',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('RKBLOW')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck rkblo1 */
      SUBROUTINE RKBLO1(IREP,IC,NCAR,SMAT,VMAT,LDM,NEFF,NLDP,STOL,EIG,
     &                  TMAT,IRKB,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate Lowdin's canoncial orthonormalization matrix with
C     the restricted kinetic balance (RKB) transformation embedded.
C
C     For each center and L value a transformation matrix is extracted
C     from spherical transformation matric CSP according to allowed symmetries.
C       IFUN(*,1) - points to Cartesian component in full list
C       IFUN(*,2) - points to spherical component in full list
C       IFUN(*,3) - points to Cartesian comp. in symmetry-reduced list
C     The symmetry-reduced spherical transformation matrix is
C     stored in CRED
C     For each L value
C       IRED(1,*) - gives number of symmetry-reduced Cartesians KXYZ
C       IRED(2,*) - gives number of symmetry-reduced sphericals KLM
C       IRED(3,*) - gives offsets to each L value
C     For each L value CRED has the dimension CRED(KXYZ,KLM)
C     The RKB part is adjoined to CRED and has the dimension
C        CRED(KRKB,KLMM2) where KRKB = KXYZ+KXYZM2
C     KXYZM2 and KLMM2 are the number of symmetry-reduced Cartesian
C     and spherical components, respectively, of the L-2 partner.
C
C

C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
C
#include "dgroup.h"
#include "pgroup.h"
      LOGICAL DONRM
      CHARACTER COMP(2)*1
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),TMAT(NCAR,NCAR),
     &          IRKB(*),WORK(*)
      DATA COMP/'L','S'/
C
C**************************************************
C***  Generate spherical transformation matrix  ***
C**************************************************
C
      DONRM = .TRUE.
      CALL RKBCAR(TMAT,IREP,IC,NCAR,NSPH,
     &            DONRM,IRKB,IPRINT,WORK,LWORK)
C
C*********************************************************
C***  Transform overlap matrix to spherical basis and
C***  generate Lowdin canonical orthonormalization matrix
C*********************************************************
C
      CALL QTRANS('AOMO','S',D0,NCAR,NCAR,NSPH,NSPH,
     &            SMAT,LDM ,NCAR,1,IPQTOQ(1,0),
     &            VMAT,LDM ,NSPH,1,IPQTOQ(1,0),
     &            TMAT,NCAR,NSPH,1,IPQTOQ(1,0),
     &            TMAT,NCAR,NSPH,1,IPQTOQ(1,0),
     &            WORK,LWORK,IPRINT)
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A,1X,A3,2X,A1)')
     &    'RKBLOW: Overlap matrix in modified solid harmonics for',
     &    REP(IREP),COMP(IC)
        CALL OUTPUT(VMAT,1,NSPH,1,NSPH,LDM,NSPH,-1,LUPRI)
      ENDIF
      CALL LOWGEN(VMAT,NSPH,SMAT,NEFF,LDM,STOL,IPRINT,
     &     EIG,WORK,LWORK)
      CALL DGEMM('N','N',NCAR,NEFF,NSPH,D1,TMAT,NCAR,
     &     SMAT,LDM,D0,VMAT,LDM)
      NA   = NCAR - NSPH
      NLDP = NSPH - NEFF
      ND   = NA + NLDP
      WRITE(LUPRI,'(3X,A1,3X,A3,3X,3(A,I10),A,E8.2)')
     &     COMP(IC),REP(IREP),'* Deleted: ',ND,
     &     '(Proj: ',NA,', Lindep: ',NLDP,') Smin: ',
     &     EIG(NSPH)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lowgen */
      SUBROUTINE LOWGEN(SMAT,NS,VMAT,NV,LDM,SSTOL,IPRINT,
     &                  EIG,WORK,LWORK)
C*****************************************************************************
C
C     Generate Lowdin matrix
C     (Lowdin's canonical orthonormalization)
C
C     V_ij = U_ij/SQRT(s_j) where U_ij is the matrix that diagonalizes
C     the overlap matrix S and s_j is eigenvalue j.
C     Linear dependence is removed by eliminating columns of V
C     corresponding to eigenvalues of the overlap matrix below
C     the given threshold.
C
C     Written by T.Saue Sep 11 1995
C     Last revision Apr 1 1998 - tsaue
C                    Febr2006 - MI
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dummy.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),WORK(LWORK)
C
#include "memint.h"
      CALL QENTER('LOWGEN')

CMI .. mover before diagonalization which can destroy SMAT
      IF(IPRINT.GE.6) THEN
        CALL HEADER('LOWGEN: Lowdin input (overlap) matrix:',-1)
        WRITE(LUPRI,'(3X,A,2I5)') 'size :',NS,NS
        CALL OUTPUT(SMAT,1,NS,1,NS,LDM,NS,-1,LUPRI)
      ENDIF

#ifdef MOD_MATLAB_LOG
      call matexport_text('% Entering LOWGEN')
      call matexport_double2('S',SMAT,LDM,NS,NS)
#endif
C
C     Diagonalize overlap matrix
C
      IF (.NOT.DOQJACO) THEN
         IF (LOWJACO) THEN
            CALL RSJACO(LDM,NS,NS,SMAT,EIG,1,-1,0,VMAT)
         ELSE
            CALL QDIAG(1,NS,SMAT,LDM,NS,EIG,1,VMAT,LDM,NS,WORK,LWORK
     $           ,IERR)
c    QDIAG gives the reverse value ordering from RSJACO, so fix that here.
c    It would be better to rewrite the code below instead..
            CALL ORDER3(VMAT,EIG,LDM,NS,NS,-1)
            IF (IERR.NE.0) THEN
               CALL QUIT('Diagonalization of overlap matrix failed')
            ENDIF
         ENDIF
      ELSE
CMI   ... (quaternion) Jacobi diagonalization routine
        NSS=NS*NS
        CALL MEMGET2('REAL','SMAT',KSMAT,NSS,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','VMAT',KVMAT,NSS,WORK,KFREE,LFREE)
        IPOS=0
        DO II=1,NS
        DO JJ=1,NS
          WORK(KSMAT+IPOS)=SMAT(JJ,II) ! needed NxN matrix for input
          IPOS=IPOS+1
        ENDDO
        ENDDO
        CALL QJACOBI(WORK(KSMAT),WORK(KVMAT),NS,1,0,IDUMMY,
     &               .TRUE.,IPRINT)
        IPOS=0
        DO II=1,NS
          EIG(II)=WORK(KSMAT+(NS*(II-1)+II-1)) ! extract eigenvalues
        DO JJ=1,NS
          VMAT(JJ,II)=WORK(KVMAT+IPOS) ! extract eigenvectors
          IPOS=IPOS+1
        ENDDO
        ENDDO
        CALL ORDER3(VMAT,EIG,LDM,NS,NS,-1) ! resort
        CALL MEMREL('LOWGEN',WORK,KSMAT,KSMAT,KFREE,LFREE)
      ENDIF
!    ... ensure V_ij = U_ij/SQRT(s_j)
      NV = 0
      DO I = 1,NS
        IF(EIG(I).LE.SSTOL) GOTO 10
        NV  = NV + 1
        FAC = D1/SQRT(EIG(I))
        CALL DSCAL(NS,FAC,VMAT(1,I),1)
      ENDDO
 10   CONTINUE
C
C     Print section
C
      IF(IPRINT.GE.2) THEN
          WRITE(LUPRI,'(7X,A,1P,D8.1)')
     &   'LOWGEN: Smallest ".LINDEP" test value of a kept orbital:',
     &   EIG(NV)
      ENDIF

C
      IF(IPRINT.GE.2) THEN
        CALL HEADER('LOWGEN: Eigenvalues in Lowdin diagonalization:',-1)
        DO I = 1,NS
          WRITE(6,*) I,EIG(I)
        ENDDO
      ENDIF
C
      IF(IPRINT.GE.6) THEN
        CALL HEADER(
     &  'LOWGEN: Canonical Lowdin matrix (i.e. eigenvectors):',-1)
        CALL OUTPUT(VMAT,1,NS,1,NV,LDM,NS,-1,LUPRI)
      ENDIF

#ifdef MOD_MATLAB_LOG
      call matexport_double2('V',VMAT,LDM,NS,NV)
      call matexport_text('% Leaving LOWGEN')
#endif

      CALL QEXIT('LOWGEN')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lowgenq */
      SUBROUTINE LOWGENQ(NZ,N,NEFF,SMAT,LRS,LCS,VMAT,LRV,LCV,
     &                   SSTOL,IPRINT,EIG,IPQ,WORK,LWORK)
C*****************************************************************************
C
C     Generate Lowdin matrix
C     (Lowdin's canonical orthonormalization)
C
C     V_ij = U_ij/SQRT(s_j) where U_ij is the matrix that diagonalizes
C     the overlap matrix S and s_j is eigenvalue j.
C     Linear dependence is removed by eliminating columns of V
C     corresponding to eigenvalues of the overlap matrix below
C     the given threshold.
C
C     Written by T. Saue Sep 2013
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
      DIMENSION SMAT(LRS,LCS,NZ),VMAT(LRV,LCV,NZ),EIG(*),IPQ(4),
     &          WORK(LWORK)
C
#include "memint.h"
      CALL QENTER('LOWGENQ')
      IF(IPRINT.GE.6) THEN
        CALL HEADER('LOWGENQ: Lowdin input (overlap) matrix:',-1)
        CALL PRQMAT(SMAT,LRS,LCS,N,N,NZ,IPQ,LUPRI)
      ENDIF
C
C     Diagonalize overlap matrix; scale with minus one to change order of eigenvalues
C
      NSDIM = LRS*LCS*NZ
      CALL DSCAL(NSDIM,DM1,SMAT,1)
      CALL QDIAG(NZ,N,SMAT,LRS,LCS,EIG,1,VMAT,LRV,LCV,
     &           WORK,LWORK,IERR)
      CALL DSCAL(N,DM1,EIG,1)
!    ... ensure V_ij = U_ij/SQRT(s_j)
      NEFF = 0
      DO I = 1,N
        IF(EIG(I).LE.SSTOL) GOTO 10
        NEFF  = NEFF + 1
        FAC = D1/SQRT(EIG(I))
        DO IZ = 1,NZ
          CALL DSCAL(N,FAC,VMAT(1,I,IZ),1)
        ENDDO
      ENDDO
 10   CONTINUE
      CALL QEXIT('LOWGENQ')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sym_ortho */
      SUBROUTINE SYM_ORTHO(NZ,N,NEFF,SMAT,LRS,LCS,VMAT,LRV,LCV,
     &                   SSTOL,IPRINT,EIG,IPQ,WORK,LWORK)
C*****************************************************************************
C
C     Generate S^{-1/2} matrix for symmetric orthonormalization
C
C     Linear dependence is removed by eliminating columns of V
C     corresponding to eigenvalues of the overlap matrix below
C     the given threshold.
C
C     Written by T. Saue Sep 2013
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0,D0=0.0D0)
C
      DIMENSION SMAT(LRS,LCS,NZ),VMAT(LRV,LCV,NZ),EIG(*),IPQ(4),
     &          WORK(LWORK)
C
#include "memint.h"
      CALL QENTER('SYM_ORTHO')
      IF(IPRINT.GE.6) THEN
        CALL HEADER('SYM_ORTHO: Lowdin input (overlap) matrix:',-1)
        CALL PRQMAT(SMAT,LRS,LCS,N,N,NZ,IPQ,LUPRI)
      ENDIF
C
C     Diagonalize overlap matrix; scale with minus one to change order of eigenvalues
C
      NSDIM = LRS*LCS*NZ
      CALL DSCAL(NSDIM,DM1,SMAT,1)
      CALL QDIAG(NZ,N,SMAT,LRS,LCS,EIG,1,VMAT,LRV,LCV,
     &           WORK,LWORK,IERR)
      CALL DSCAL(N,DM1,EIG,1)
!    ... ensure V_ij = U_ij/SQRT(s_j)
      NEFF = 0
      IF(EIG(1).LE.SSTOL) GOTO 10
      NEFF = NEFF + 1
      ALPHA = D1/SQRT(EIG(1))
      CALL DENST1(SMAT,LRS,LCS,NZ,ALPHA,D0,
     &              VMAT,LRV,LCV,1,1,N)
      DO I = 2,N
        IF(EIG(I).LE.SSTOL) GOTO 10
        NEFF  = NEFF + 1
        ALPHA = D1/SQRT(EIG(I))
        CALL DENST1(SMAT,LRS,LCS,NZ,ALPHA,D1,
     &              VMAT,LRV,LCV,I,1,N)
      ENDDO
 10   CONTINUE
      DO IZ = 1,NZ
        DO J = 1,NEFF
          CALL DCOPY(N,SMAT(1,J,IZ),1,VMAT(1,J,IZ),1)
        ENDDO
      ENDDO
      CALL QEXIT('SYM_ORTHO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tscorr */
      SUBROUTINE TSCORR(FMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Generate SS nuclear attraction integrals to mimic
C     two-electron SS-integrals
C
C     Written by T.Saue July 1997
C     Last revision July 15 1997 - tsaue
C
C*****************************************************************************
      use mxcent_mod
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
      DIMENSION WORK(LWORK),CBUF(2),GBUF(2),FMAT(*)
      LOGICAL DOINT(2,2)
      CHARACTER OMITVNUC(2)*4
#include "dgroup.h"
#include "dcbbas.h"
#include "nuclei.h"
C
      CALL QENTER('TSCORR')
#include "memint.h"
C     Memory allocation
      CALL MEMGET2('REAL','ONEMT',KONEMT,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','ONEIN',KONEIN,NNBBASX,WORK,KFREE,LFREE)
C
C     Buffer info
C
      DO I = 1,NUCIND
        CBUF(I) = CHARGE(I)
        GBUF(I) = GNUEXP(I)
        GNUEXP(I) = 0.0D0
      ENDDO
      call quit('TSCORR not properly implemented in this version')
C     HJAaJ: this is only for a specific datomic!!! (TeSr ??)
C
      CHARGE(1) = -0.1812D0
      CHARGE(2) = -0.0832D0
C
C     Generate integrals
C
      DOINT(1,1) = .FALSE.
      DOINT(2,2) = .TRUE.
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.

      OMITVNUC(1) = 'TTTT'
      OMITVNUC(2) = 'FFFF'
      CALL NUCAT1(WORK(KONEIN),WORK(KFREE),LFREE,
     &            DOINT,OMITVNUC,IPRINT)
C
      CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
      CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
      IF(IPRINT.GE.6) THEN
        CALL HEADER('SS nuclear attraction integrals',-1)
          CALL PRQMAT(WORK(KONEMT),NTBAS(0),NTBAS(0),
     &                NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
      ENDIF
      CALL DAXPY(N2BBASX,D1,WORK(KONEMT),1,FMAT,1)
C
C     Memory deallocation
      CALL MEMREL('TSCORR',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Buffer retrieve
C
      DO I = 1,NUCIND
        CHARGE(I) = CBUF(I)
        GNUEXP(I) = GBUF(I)
      ENDDO
C
      CALL QEXIT('TSCORR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bncorr */
      SUBROUTINE BNCORR(FMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     BNCORR: Generate screened bare nucleus matrix
C     Jun 2003 Jesper Kielberg Pedersen (based on OLDBNCR)
C
C     Description :
C    ===============
C     Extension of BNCORR that in a simple way takes the shell-structure
C     into account. For each center in a molecule we estimate the screening
C     neglected in the bare nucleus approximation by a sum of the
C     contributions from each shells (n) in the atoms :
C           sum_(n) <X_A|(-Z(n)*G(n,a,r_C))/r_C|X_B>
C           X=L,S ; G(a,r)=exp(-a*r^2)
C
C     On input/ouput:     FMAT - assumed one-electron (2c/4c bare nucleus) Fock matrix (input),
C    =================    screening bare nucleus added (ouput)
C
C     Called from: PREDHF/dirscf.F
C                  RH1DIAG/../krmc/krmcopt.F
C
C      In the case of two-component mode do the picture change transformation
C      of the four-compo. correction matrix. In pure two-component mode add only
C      LL block.
C
C    Last modifications: March 2006/MI,Strasbourg
C
C*****************************************************************************
!     use x2c_fio, only: x2c_read
      use primitives_mod
      use maxaqn_mod
      use mxcent_mod
      use symmet_mod
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0)
C
      LOGICAL DOINT(2,2), FINISH, ISONBSSMAT, WEHAVE
      CHARACTER*4 OMITVNUC
      CHARACTER*8 BNCRLAB
      REAL*8    WORK(LWORK),FMAT(*)
      logical   fndlab12
      character*12 flabel
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "cbihr1.h"
#include "cbisol.h"

      real(8)  :: TIMSTR(2), TIMEND(2)
C
      CALL QENTER('BNCORR')
c
c     Start timer
c
      CALL TIMER2('START',TIMSTR,TIMEND)
c
#include "memint.h"

      WEHAVE = .FALSE.

      if(bss)then
        WEHAVE = ISONBSSMAT('BNCR2C_2','BNCORR',LUBSS,IPRHAM).OR.
     &           ISONBSSMAT('BNCR2C_4','BNCORR',LUBSS,IPRHAM)

      else if(x2c)then

       open(lux2c,file='X2CMAT',status='old',form='unformatted',
     &      access='sequential',action='read',position='rewind')

       write(flabel,'(a11,i1)') 'h1bncc2cAOl',0
       if(fndlab12(flabel,lux2c)) then
         wehave = .true.
       else
         close(lux2c,status='keep')
       end if
      end if

      IF (WEHAVE) THEN
CMI   ... read 2c screened bare nucleus from the BSSMAT file in 4c or 2c metric
        if(bss)then
          if(twocomp)then
            BNCRLAB='BNCR2C_2'
          else
            BNCRLAB='BNCR2C_4'
          end if
          CALL RFBSSMAT(BNCRLAB,.FALSE.,'BNCORR',
     &                  FMAT,N2BBASXQ,LUBSS,IPRHAM)
        else if(x2c)then
          call x2c_read(flabel,fmat,n2bbasx * nz,lux2c)
          close(lux2c,status='keep')
        end if
      ELSE

        if(mdirac)then
        write(lupri,*)
     &  ' bare-nucleus correction not implemented yet for MDIRAC'
        write(lupri,*)
     &  ' program will continue with the default 1-el Hamiltonian'
          goto 999
        end if

CMI   ... continue the standard way
C     Memory allocation
      CALL MEMGET2('REAL','ONEMT',KONEMT, N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','ONEIN',KONEIN, NNBBASX,WORK,KFREE,LFREE)
      LFCBA = NUCDEP**3
      CALL MEMGET2('REAL','FCBA' ,KFCBA ,  LFCBA ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','COORC',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','SIGNC',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','NCENT',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','JSYMC',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','JCENT',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','GEXP' ,KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','XGEXP',KXGEXP,  NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','CHRG' ,KCHRG ,  NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','ESHEL',KESHEL,7*NUCIND,WORK,KFREE,LFREE)
!     allocate(eshell_vec(7*NUCIND))
      LEN_WRK   = LWORK - KFREE + 1
      CALL DZERO(WORK(KCHRG),NUCIND)
      CALL DZERO(WORK(KONEMT),N2BBASX)
      CALL DZERO(WORK(KESHEL),7*NUCIND)
C
C     Initialization (calculate LL and SS)
C
      DOINT(1,1) = .TRUE.
      IF (NOSMLV.OR.BSS.OR.x2c) THEN
         DOINT(2,2) = .FALSE.
         IF (IPRINT.GE.1) WRITE(LUPRI,'(/A)')
     &     '  BNCORR: SS contribution to the bare nuclei omitted !'
      ELSE
         DOINT(2,2) = .TRUE.
      END IF
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.
c
c   ... find max number of occupied shells
c
      Z = D0
      DO I = 1,NUCIND
         IF (.not.NOORBT(I)) Z = MAX(Z,CHARGE(I)) ! skip point charges
      ENDDO

C    ... find maximum number of shells we need to run over
C        and return number in NSHELL
      CALL BNCOR1(GNU,ZEFF,0,Z,NSHELL,0,DUMMY,NUCIND,.TRUE.)

c   ... loop over shells
c       in the first shell (N=1) we only place one electron to avoid
c       self-interaction. This electron is delocalized (to some degree)
c       in the zeroth run (N=0) by using a diffuse Gaussian (see BNCOR1)
c
      IF (IPRINT .GE. 1) THEN
         WRITE(LUPRI,'(/4X,A/4X,A//14X,A/14X,A)')
     &' *** BNCORR for < X_A | (Qeff*G(a,r_C))/r_C | X_B > ***',
     &'     X=L,S ; G(a,r_C)=exp(-a*r_C^2)',
     &      ' Nucleus   shell charge=Qeff  Exponent(a)',
     &      '-----------------------------------------'
      ENDIF
      DO N = 0,NSHELL
c     ... In Slater's rules we write the electron configuration as
c         [1s] [2s2p] [3s3p][3d] [4s4p][4d][4f]... so for each shell we have
c          1     1        2            3           loops
         IF (N.LE.1) THEN
            NSLOOP = 1
         ELSE
            NSLOOP = N - 1
         ENDIF
c        ... we will not go beyond f-orbitals
         NSLOOP = MIN(NSLOOP,3)
         DO L = 1,NSLOOP
c
         ZTEST = D0
         DO I = 1,NUCIND
            IF (NOORBT(I) .OR. CHARGE(I).EQ.0.0D0) THEN ! skip point charges and floating orbitals
               WORK(KXGEXP-1+I) = 0.0D0
               WORK(KCHRG -1+I) = 0.0D0
            ELSE
c           ... Get charges and exponents for the electrons in this group
               CALL BNCOR1(GNU,ZEFF,N,CHARGE(I),L,I,WORK(KESHEL),
     &              NUCIND,.FALSE.)
               ZTEST = MAX(ZTEST,ZEFF)
               WORK(KXGEXP-1+I) =  GNU
               WORK(KCHRG -1+I) = -ZEFF
            END IF
         ENDDO ! I = 1,NUCIND
c
         IF(IPRINT .GE. 1) THEN
            IF (N .LE. 1) THEN
               WRITE(LUPRI,'(14X,A,I2,A)')
     &         '-- Subshell n =',N,', l = 0'
            ELSE IF (L .EQ. 1) THEN
               WRITE(LUPRI,'(14X,A,I2,A)')
     &         '-- Subshell n =',N,', l = 0,1'
            ELSE
               WRITE(LUPRI,'(14X,A,I2,A,I2)')
     &         '-- Subshell n =',N,', l =',L
            END IF
            DO I = 1,NUCIND
               IF (WORK(KCHRG-1+I) .NE. D0)
     &         WRITE(LUPRI,'(15X,A5,F15.2,F15.2)')
     &         NAMN(I),WORK(KCHRG-1+I),WORK(KXGEXP-1+I)
            ENDDO
         ENDIF
c
         IF (ZTEST .NE. D0) THEN
C
C     Omit selected terms depending on OMITVNUC (See CATNUC for details)
C     In this version one-center terms are calculated with the same
C     charges/exponents as multi-center terms.
C
         OMITVNUC = 'FFFF'
c
c     ... Make list of nuclei.
c
         CALL CATNUC(NUCDEP,NCENTC,WORK(KCHRG),WORK(KFCBA),WORK(KGEXP),
     &            WORK(KXGEXP),WORK(KCOORC), WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),OMITVNUC,IPRONE)
c
c    ... Calculate integrals (both LL and SS)
c
         CALL CATDR1(WORK(KONEMT),DUMMY,NCENTC,WORK(KFCBA),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KFREE),LEN_WRK,IPRONE,
     &            .FALSE.,0,NNBBASX,DOINT)
c
         END IF ! if (ztest .ne. d0)
c  ... end loop within shells
         ENDDO ! DO L = 1,NSLOOP
c  ... end loop over shells
      ENDDO ! DO N = 0,NSHELL

      IF (IPRINT .GE. 1) WRITE(LUPRI,'(14X,A/)')
     &      '-----------------------------------------'
C
C     Symmetry unpack integrals
C
      CALL SYMUPK(WORK(KONEMT),WORK(KONEIN),1,NNBBASX)
      CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
      CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
C
      IF(IPRINT.GE.11) THEN
         CALL HEADER('BNCORR correction for DC FMAT',-1)
         CALL OUTPUT(WORK(KONEMT),1,NTBAS(0),1,NTBAS(0),
     &        NTBAS(0),NTBAS(0),-1,LUPRI)
      END IF

C
C     Add correction to (entering) Fock-matrix, FMAT
C
      CALL DAXPY(N2BBASX,D1,WORK(KONEMT),1,FMAT,1)
      IF(IPRINT.GE.11) THEN
         CALL HEADER('BNCORR: screening corrected bare nucleus FMAT',-1)
         CALL PRQMAT(FMAT,NTBAS(0),NTBAS(0),
     &               NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
      END IF
C
C     Memory deallocation
      CALL MEMREL('BNCORR',WORK,KWORK,KWORK,KFREE,LFREE)

      ENDIF

C
c     Stop timer
c
 999  IF (IPRINT.GE.2)
     &  CALL TIMER2('>>> Time used in BNCORR is',TIMSTR,TIMEND)
      CALL QEXIT('BNCORR')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bncor1 */
      SUBROUTINE BNCOR1(GNU,ZEFF,N,Z,L,INUC,SHELE,NNUC,MAXSHELL)
C*****************************************************************************
C
C     BNCOR1 :  Find effective charges and Gaussian exponent
C               (from Slater rules) for a given shell (N) for
C               a given atom (Z)
C
C     Jun 2003 Jesper Kielberg Pedersen (used for BNCORR)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0,
     &          D6 = 6.0D0, D8 = 8.0D0, D10=10.0D0, D14=14.0D0,
     &          D18=18.0D0)
      LOGICAL MAXSHELL
      DIMENSION STO1G(8), XN_EFF(8), SHELE(NNUC,8)
      SAVE STO1G
C     STO1G(n) is the exponent of the 1s-Gaussian that
C     best fits the ns Slater-type atomic orbitals of exponent 1.
C     (Robert F. Stewart, JCP 52(1),431 (Jan 1970))
C     This article provides exponent up to the n=5 shell. The exponents
C     for the last two shells have been generated from a fit :
C     (y = 0.28196 x^-1.5532)
      DATA STO1G /0.2709D0, 0.1012D0, 0.05297D0,0.03265D0,0.02217D0,
     &            0.01744D0,0.01373D0,0.01116D0/
!     DATA XN_EFF/1.0D0, 2.0D0, 3.0D0, 3.7D0, 4.0D0, 4.2D0, 4.4D0/
      DATA XN_EFF/1.0D0, 2.0D0, 3.0D0, 3.7D0, 4.4D0, 5.1D0, 5.8D0, 
     & 6.5D0/
C     ... Slater's n_eff according to wikipedia for n=1:6,
C         I extrapolated to n=7 /hjaaj Sep 2011.
C     ... Empirical modification Feb. 2012 hjaaj: XeF6 did not
C         converge with the Sep 2011 values; on the other hand
C         we had an example in Sep 2011 which benefitted from the
C         n_eff values.
C
C     Table of s,p,d,f,g  orbital occupations of the elements
C     (from which the electron-configuration can be deduced)
      DIMENSION NOCC(138*5)
c     DIMENSION NOCC(127*5)
      SAVE NOCC
      DATA NOCC
C     Periode 1
     &          / 1, 0, 0, 0, 0, 2, 0, 0, 0, 0,
C     Periode 2
     &            3, 0, 0, 0, 0, 4, 0, 0, 0, 0,
     &            4, 1, 0, 0, 0, 4, 2, 0, 0, 0,
     &            4, 3, 0, 0, 0, 4, 4, 0, 0, 0, 
     &            4, 5, 0, 0, 0, 4, 6, 0, 0, 0,
C     Periode 3
     &            5, 6, 0, 0, 0, 6, 6, 0, 0, 0,
     &            6, 7, 0, 0, 0, 6, 8, 0, 0, 0,
     &            6, 9, 0, 0, 0, 6,10, 0, 0, 0,
     &            6,11, 0, 0, 0, 6,12, 0, 0, 0,
C     Periode 4
     &            7,12, 0, 0, 0, 8,12, 0, 0, 0,
     &            8,12, 1, 0, 0, 8,12, 2, 0, 0,
     &            8,12, 3, 0, 0, 8,12, 4, 0, 0,
     &            8,12, 5, 0, 0, 8,12, 6, 0, 0,
     &            8,12, 7, 0, 0, 8,12, 8, 0, 0,
     &            8,12, 9, 0, 0, 8,12,10, 0, 0,
     &            8,13,10, 0, 0, 8,14,10, 0, 0,
     &            8,15,10, 0, 0, 8,16,10, 0, 0,
     &            8,17,10, 0, 0, 8,18,10, 0, 0,
C     Periode 5 (5s block)
     &            9,18,10, 0, 0, 10,18,10, 0, 0,
C     (4d block)
     &           10,18,11, 0, 0, 10,18,12, 0, 0, 
     &           10,18,13, 0, 0, 10,18,14, 0, 0, 
     &           10,18,15, 0, 0, 10,18,16, 0, 0,
     &           10,18,17, 0, 0, 10,18,18, 0, 0, 
     &           10,18,19, 0, 0, 10,18,20, 0, 0,
C     (5p block)
     &           10,19,20, 0, 0, 10,20,20, 0, 0, 
     &           10,21,20, 0, 0, 10,22,20, 0, 0, 
     &           10,23,20, 0, 0, 10,24,20, 0, 0,
C     Periode 6 (6s block)
     &           11,24,20, 0, 0, 12,24,20, 0, 0,
     &           12,24,21, 0, 0,
C    --> Lanthanides (4f block)
     &           12,24,21, 1, 0, 12,24,20, 3, 0, 
     &           12,24,20, 4, 0, 12,24,20, 5, 0, 
     &           12,24,20, 6, 0, 12,24,20, 7, 0,
     &           12,24,20, 8, 0, 12,24,20, 9, 0, 
     &           12,24,20,10, 0, 12,24,20,11, 0, 
     &           12,24,20,12, 0, 12,24,20,13, 0,
     &           12,24,20,14, 0, 12,24,21,14, 0,
C    <-- Lanthanides; Hf -> Hg (5d block)
     &           12,24,22,14, 0, 12,24,23,14, 0,
     &           12,24,24,14, 0, 12,24,25,14, 0,
     &           12,24,26,14, 0, 12,24,27,14, 0,
     &           12,24,28,14, 0, 12,24,29,14, 0, 
     &           12,24,30,14, 0,
C     Tl -> Rn (6p block)
     &           12,25,30,14, 0, 12,26,30,14, 0, 
     &           12,27,30,14, 0, 12,28,30,14, 0, 
     &           12,29,30,14, 0, 12,30,30,14, 0,
C     Periode 7 (7s block)
     &           13,30,30,14, 0, 14,30,30,14, 0,
     &           14,30,31,14, 0,
C    --> Actinides
     &           14,30,32,14, 0, 14,30,31,16, 0, 
     &           14,30,31,17, 0, 14,30,31,18, 0, 
     &           14,30,30,20, 0, 14,30,30,21, 0,
     &           14,30,31,21, 0, 14,30,30,23, 0, 
     &           14,30,30,24, 0, 14,30,30,25, 0, 
     &           14,30,30,26, 0, 14,30,30,27, 0,
     &           14,30,30,28, 0, 14,30,31,28, 0,
C    <-- Actinides; (6d block and 7p block)
     &           14,30,32,28, 0, 14,30,33,28, 0,
     &           14,30,34,28, 0, 14,30,35,28, 0, 
     &           14,30,36,28, 0, 14,30,37,28, 0,
     &           13,30,39,28, 0, 13,30,40,28, 0, 
     &           14,30,40,28, 0, 14,31,40,28, 0, 
     &           14,32,40,28, 0, 14,33,40,28, 0,
     &           14,34,40,28, 0, 14,35,40,28, 0, 
     &           14,36,40,28, 0,
C     Periode 8 (8s block)
     &           15,36,40,28,0,  16,36,40,28,0,
     &           16,37,40,28,0,
C    --> superlanthanides (122-138)
     &           16,37,41,28,0,  16,37,41,29,0,
     &           16,37,41,30,0,  16,37,41,30,1,
     &           16,37,41,30,2,  16,38,40,30,3,
     &           16,38,40,30,4,  16,38,40,30,5,
     &           16,38,40,30,6,  16,38,40,30,7,
     &           16,38,40,30,8,  16,38,40,31,8,
     &           16,38,40,32,8,  16,38,40,32,9,
     &           16,38,40,32,10, 16,38,41,31,11,
     &           16,38,41,31,12/
      CALL QENTER('BNCOR1')
c
      IF (NINT(Z) .GT. 138) THEN
         WRITE(LUPRI,'(2(A/))') ' FATAL : The bare nucleus '//
     &   'correction is only implemented for elemtents 1-138',
     &   '         You can disable BNCORR with the .NOBNCR keyword '//
     &   'in *SCF'
         CALL QUIT('BNCORR : Unsupported element!')
      ENDIF

      ZEFF = D0
      GNU  = D0

C     ... If this is a ghost-atom (f.ex for Counter-Poise
C         correction) we return (charge is zero in BNCORR)
      IF (Z .EQ. D0) GOTO 9000

      IF (MAXSHELL) THEN
C    ... just find maximum number of shells we need to run over
C        and return number in L-variable (used in BNCORR)
         NST = NOCC((NINT(Z)-1)*5 + 1)
         L = NST/2
         IF(MOD(NST,2).NE.0) L = L + 1
         GOTO 9000
      ENDIF
C
      IF (N.EQ.0) THEN
c     ...Always make 1 electron semi-delocalized
         ZEFF = D1
         GNU  = 0.4D0
         SHELE(INUC,1) = SHELE(INUC,1) + D1
         GOTO 9000
      ENDIF
c     ... If Hydrogen we are done
      IF (Z.EQ.D1) GOTO 9000
c     ... special case for the (2-1) 1s electron
      IF (N.EQ.1) THEN
         ZEFF = D1
         SEXP = Z - 0.35
         IF (Z.EQ.2) SEXP = SEXP + 0.05
         SHELE(INUC,1) = SHELE(INUC,1) + D1
         GOTO 20
      ENDIF
c     ... Generel case. We are past H and He and have thus accounted
c         for two electrons.
c     ... Find numbers of s,p,d,f,g-electrons for element (Z)
         NST=0
         NPT=0
         NDT=0
         NFT=0
         NGT=0
      IF (Z.GE.1.0D0) THEN ! miro: prevent NOCC of out-of-bounds
         NST = NOCC((NINT(Z)-1)*5 + 1)
         NPT = NOCC((NINT(Z)-1)*5 + 2)
         NDT = NOCC((NINT(Z)-1)*5 + 3)
         NFT = NOCC((NINT(Z)-1)*5 + 4)
         NGT = NOCC((NINT(Z)-1)*5 + 5)
      ENDIF
C
c     ... Now figure out the occupation used for the Slater rules (N >= 2) :
c      [1s][2s2p][3s3p][3d][4s4p][4d][4f]...
c
c     =================================
c     How many s-electrons in this shell
c     =================================
c
      NS = NST - (N-1)*2
c
      IF (NS .LE. 0) THEN
         NS = 0
      ELSEIF (NS.GE.2) THEN
         NS = 2
      ENDIF
c
c     =================================
c     How many p-electrons in this shell
c     =================================
c
      NP = NPT - (N-2)*6
c
      IF (NP .LE. 0) THEN
         NP = 0
      ELSEIF (NP.GE.6) THEN
         NP = 6
      ENDIF
c
c     =================================
c     How many d-electrons in this shell
c     =================================
c
      ND = NDT - (N-3)*10
c
      IF (ND .LE. 0) THEN
         ND = 0
      ELSEIF (ND.GE.10) THEN
         ND = 10
      ENDIF
c
c     =================================
c     How many f-electron in this shell
c     =================================
c
      NF = NFT - (N-4)*14
c
      IF (NF .LE. 0) THEN
         NF = 0
      ELSEIF (NF.GE.14) THEN
         NF = 14
      ENDIF
c     =================================
c     How many g-electron in this shell
c     =================================
c
      NG = NGT - (N-5)*18
c
      IF (NG .LE. 0) THEN
         NG = 0
      ELSEIF (NG.GE.18) THEN
         NF = 18
      ENDIF
      IF (L.EQ.1) ZEFF = NS + NP
      IF (L.EQ.2) ZEFF = ND
      IF (L.EQ.3) ZEFF = NF
      IF (L.EQ.4) ZEFF = NG
      SHELE(INUC,N) = SHELE(INUC,N) + ZEFF
      IF (ZEFF.EQ.D0) GO TO 9000
C
C     ... Finally use Slater's rules to estimate screening
C         (and orbital exponent)
C
      IF (L.EQ.1) THEN
c        [nsnp]
         EMN = D0
         DO I = 1,N-2
            EMN = EMN + SHELE(INUC,I)
         ENDDO
         SEXP = (Z - ((ZEFF-1.0)*0.35D0 + SHELE(INUC,N-1)*0.85D0 + EMN))
         SEXP = SEXP / XN_EFF(N)
      ELSE
c        [nd],[nf], [ng]
         EMN = D0
         DO I = 1,N-1
            EMN = EMN + SHELE(INUC,I)
         ENDDO
         SEXP = (Z - ((ZEFF-1.0D0)*0.35D0 + EMN)) / XN_EFF(N)
      ENDIF
c
   20 CONTINUE
c     ... Convert STO to GTO using 1 Gaussian (STO1G)
      GNU = D2*SEXP*SEXP*STO1G(N)
!     write (lupri,*) 'BNCOR1: N,L,STO1G(N),SEXP,GNU',
!    &                         N,L,STO1G(N),SEXP,GNU
C     ... Make sure we dont get more delocalized than the first
C         1s electron.
C     GNU = MAX(GNU,0.50D0)
 9000 CONTINUE
      CALL QEXIT ('BNCOR1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck eonecap */
      SUBROUTINE EONECAP(FMAT,DMAT,WORK,LWORK,IQ)
C****************************************************************************
C     PURPOSE:
C        Calculate correction to the nuclear attraction energy
C        from the one-center models
C
C      E = D(IJ)(*)F(IJ)
C
C     Sep 2001 : Jesper Kielberg Pedersen
C     Based on ERGCAL
C     Revised Juli 2002 : jkp
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D2 = 2.0D0)
C
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
      CHARACTER FMT*6,MXFORM*6
      DIMENSION FMAT(*), DMAT(*),WORK(LWORK)
C
      CALL QENTER('EONECAP')
C
      IF(INTV1C.EQ.3.OR.(INTV1C.EQ.2.AND.ICTLV1C(2).EQ.0)) THEN
         IF (NOPEN .GE. 1) THEN
            CALL QUIT('EONECAP not implemented for open shell '//
     &                'for this ONECAP model yet.')
         END IF
      END IF
C
      ECORR = D2*DDOT(N2BBASX,DMAT,1,FMAT,1)
              ! Factor D2 comes from Time-Reversal Symmetry.
C
      CALL HEADER('Output from EONECAP',-1)
      FMT = MXFORM(ABS(ECORR),20)
      IF (IQ .EQ. 1) THEN
         WRITE(LUPRI,'(A,'//FMT//')')
     &  'ONECAP correction to LL-Nuc att.  :  ',ECORR
      ELSE IF (IQ .EQ. 2) THEN
         WRITE(LUPRI,'(A,'//FMT//')')
     &  'ONECAP correction to SS-Nuc att.  :  ',ECORR
      ELSE IF (IQ .EQ. 3) THEN
         WRITE(LUPRI,'(A,'//FMT//')')
     &  'Estimation of error on the energy in this model :  ',ECORR
      END IF
C
      CALL QEXIT('EONECAP')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiprj */
      SUBROUTINE PSIPRJ(TMAT,TINV,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Eliminate a set of fragment orbitals from the variational space;
C     this is done by projecting them out of the MO transformation matrix
C     On input:
C       TMAT contains the MO transformation matrix T (canonical orthogonalization)
C       TINV contains the inverse of T^{dagger}, that is, TINV = ST
C     Written by T.Saue March 26 1999
C
C***********************************************************************
#include "implicit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbdhf.h"
C
      LOGICAL TOBE
      DIMENSION TMAT(*),TINV(*),WORK(*)
      DIMENSION NSTR(2,0:2,MFRAG),KVEC(2,MFRAG),NVECS(2)
      DIMENSION IOFF(2)
C
      KFRSAV = KFREE
      CALL TITLER('Reduction of variational space by projection',
     &     '*',125)
C     
C     Set dimensions for coefficients
C     ===============================
C
      NORBT = 0
      NCMOT  = 0
      DO IFRP = 1,NFSYM
        IORB(IFRP)  = NORBT
        NORBT       = NORBT  + NORB(IFRP)
        ICMO(IFRP)  = NCMOT
        ICMOQ(IFRP) = NCMOT*NZ
        NCMO(IFRP)  = NFBAS(IFRP,0)*NORB(IFRP)
        NCMOQ(IFRP) = NCMO(IFRP)*NZ
        NCMOT       = NCMOT + NCMO(IFRP)
      ENDDO
      NCMOTQ  = NCMOT*NZ
C
C     Find total number of fragment orbitals
C     =======================================
C
C     NFRAG - number of fragments
C     NSTR(IFRP,IC,IFRAG)
C     - number of fragment orbitals of fermion ircop IFRP for fragment IFRAG
C        IC = 0 total
C        IC = 1 positive energy
C        IC = 2 negative energy
C     NVECS(IFRP) - total number of fragment orbitals in fermion ircop IFRP
C     NVECT       - total number of fragment orbtials
C
      NVECT = 0
C.....loop over fermion ircops
      DO IFRP = 1,NFSYM
        NVECS(IFRP) = 0
C.......loop over fragments
        DO IFRAG = 1,NFRAG
          NSTR(IFRP,1,IFRAG) = NESH(IFRP)
          NSTR(IFRP,2,IFRAG) = NPSH(IFRP)
          CALL ORBNUM(VCPROJ(IFRP,IFRAG),IFRP,KVEC(IFRP,IFRAG),
     &                NSTR(1,0,IFRAG),WORK,KFREE,LFREE)
          NVECS(IFRP) = NVECS(IFRP) + NSTR(IFRP,0,IFRAG)
        ENDDO
        NVECT = NVECT + NVECS(IFRP)
      ENDDO
      IF(NVECT.EQ.0) GOTO 30
C     
C     Calculate dimensions of fragment coefficients array
C     ====================================================
C
C     KQ - start adress for coefficient array
C     KE - start adress for eigenvalue array
C     KI - start adress for boson irreps info
C
      NCDIM = 0
      NEDIM = 0
      DO IFRP = 1,NFSYM
        NCDIM = NCDIM + NFBAS(IFRP,0)*NVECS(IFRP)*NZ
        NEDIM = NEDIM + NVECS(IFRP)
      ENDDO
      CALL MEMGET2('REAL','Q',KQ,NCDIM,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','E',KE,NEDIM,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','I',KI,NEDIM,WORK,KFREE,LFREE)
      CALL IZERO(WORK(KI),NEDIM)
C
C     Allocate buffer space for coefficients, eigenvalues
C     and boson irrep information
C     ===================================================
C
      CALL MEMGET2('REAL','CMO',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','EIG',KEIG,NORBT   ,WORK,KFREE,LFREE)
      IF(SUB_BL) THEN
        IOPT = 15
        CALL MEMGET2('REAL','IBEIG',KIBE,NORBT,WORK,KFREE,LFREE)
      ELSE
        IOPT = 7
        KIBE = KFREE
      ENDIF
C
C     Select fragment orbitals
C     =========================
C
      DO IFRP = 1,NFSYM
        IOFF(IFRP) = 0
      ENDDO
      CALL SELFRAG(WORK(KQ),NVECS,WORK(KE),WORK(KI),
     &             NFRAG,NSTR,PRJFIL,PROOWN,NPRJNUC,IOPT,KVEC,
     &             WORK(KCMO),WORK(KEIG),WORK(KIBE),IOFF,
     &             KRMC_FLG,.FALSE.,WORK,KFREE,LFREE)
C
C     Collect fragment info
C      
      CALL MEMGET2('INTE','SEL',KSEL,2*NVECT,WORK,KFREE,LFREE)
      CALL FRAGINFO(WORK(KSEL),NVECT,KVEC,NFRAG,NSTR,
     &              WORK,KFREE,LFREE)
C
C     Throw away the full set of coefficients and eigenvalues
C     =======================================================
C
      CALL MEMREL('PSIPRJ.cfs',WORK,KFRSAV,KCMO,KFREE,LFREE)
C
C     Reduce variational space by projection using the fragment orbitals
C     ===================================================================
C
      CALL PSIPR1(NVECS,WORK(KQ),WORK(KE),WORK(KI),TMAT,TINV,
     &            NSTR,WORK(KSEL),IPRINT,WORK,KFREE,LFREE)
 30   CONTINUE
      CALL MEMREL('PSIPRJ.final',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mk1cex */
      SUBROUTINE MK1CEX (WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C     PURPOSE:
C        Construct one-center expansion coefficients
C
C        Written by L.Visscher
C        Last revision: Aug 16 2001 - lv
C*****************************************************************************
      use maxaqn_mod
      use mxcent_mod
#include "implicit.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
#include "dcbgen.h"
#include "dcbbas.h"
#include "nuclei.h"
      DIMENSION WORK(*)
C
      CALL QENTER('MK1CEX')
      KSAVE = KFREE
C
      LDM = NTBAS(0)
      LDM2 = LDM * LDM
      CALL MEMGET2('REAL','SMAT1' ,KSMAT1 ,LDM2,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','SMAT2' ,KSMAT2 ,LDM2,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','INDAO' ,KINDAO ,LDM ,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','INDAOS',KINDAOS,LDM ,WORK,KFREE,LFREE)
      NIND = 2 * NUCDEP + 2
      CALL MEMGET('INTE',KINDNAO,NIND,WORK,KFREE,LFREE)
C
C     Get overlap matrix
C     ==================
C
      CALL GTOVLX(WORK(KSMAT1),D1)
C
C     Backtransform overlap matrix to AO-basis
C     ========================================
C
      CALL MTBSBU(WORK(KSMAT1),WORK(KSMAT2))
      CALL MTSOAO(WORK(KSMAT2),WORK(KSMAT1),NTBAS(0),0,IPRINT)
C
C     Get index arrays to sort on expansion centers (nuclei)
C     ======================================================
C
      CALL INDAOC (WORK(KINDAO),WORK(KINDAOS),WORK(KINDNAO))
C
C     Sort the overlap matrix to this order
C     =====================================
C
      CALL REORMT (WORK(KSMAT1),WORK(KSMAT2),LDM,LDM,
     &             WORK(KINDAO),WORK(KINDAO))
C
C     Construct the projection operator and put it in SMAT1
C     =====================================================
C
      CALL MK1CEX2 (LDM,WORK(KSMAT2),WORK(KSMAT1),WORK(KINDNAO),
     &              WORK,KFREE,LFREE,IPRINT)
C
C     Sort the projection operator back to the standard AO order
C     ========================================================
C
      CALL REORMT (WORK(KSMAT1),WORK(KSMAT2),LDM,LDM,
     &             WORK(KINDAOS),WORK(KINDAOS))
C
C     Write the projection operator to file for later use
C     ===================================================
C
      CALL OPNFIL(LUPMAT,'DFPMAT','UNKNOWN','MK1CEX')
      CALL WRITT(LUPMAT,N2BBASX,WORK(KSMAT2))
      CLOSE(LUPMAT,STATUS='KEEP')
C
      CALL MEMREL('MK1CEX',WORK,KSAVE,KSAVE,KFREE,LFREE)
C
      CALL QEXIT('MK1CEX')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mk1cex2 */
      SUBROUTINE MK1CEX2 (LDM,STOT,PTOT,INDNAO,WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C     PURPOSE:
C        Construct one-center expansion coefficients
C
C        Written by L.Visscher
C        Last revision: Aug 16 2001 - lv
C*****************************************************************************
      use maxaqn_mod
      use mxcent_mod
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "nuclei.h"
      DIMENSION STOT(LDM,LDM),PTOT(LDM,LDM)
      DIMENSION WORK(*)
      DIMENSION INDNAO(2,*)
C
      CALL QENTER('MK1CEX2')
      KSAVE = KFREE
      CALL DZERO (PTOT,LDM*LDM)
C
C     Loop over large/small component
C     Only the small part is considered because the
C     one-center approximation is not used in the
C     large component.
C
      DO IC = 2, 2
         NB0 = INDNAO(IC,1)+1
C
C        Loop over expansion centers
C
         DO IEXPC = 1, NUCDEP
           NB1 = INDNAO(IC,IEXPC)+1
           NB2 = INDNAO(IC,IEXPC+1)
           NLOC = NB2-NB1+1
           NV = NLOC*LDM
           NW = NLOC*NLOC
           CALL MEMGET2('REAL','EIG',KEIG,NLOC,WORK,KFREE,LFREE)
           CALL MEMGET2('REAL','V'  ,KV,NV,WORK,KFREE,LFREE)
           CALL MEMGET2('REAL','W'  ,KW,NW,WORK,KFREE,LFREE)
C
C          Make the local projection operator
C
           CALL MK1CEX3 (STOT(NB1,NB0),STOT(NB1,NB1),PTOT(NB1,NB0),LDM,
     &                   NTBAS(IC),NLOC,
     &                   WORK(KV),WORK(KW),WORK(KEIG),STOL(IC),
     &                   IPRINT,WORK,KFREE,LFREE)
           CALL MEMREL('MK1CEX2',WORK,KEIG,KEIG,KFREE,LFREE)
C
C          This operator contains also the projection on the functions of the same
C          center which should be (close to) the unit matrix. Make this part zero.
C
C          TODO : check that this part is indeed the unit matrix.
C
           DO I = NB1, NB2
              CALL DZERO(PTOT(NB1,I),NLOC)
           ENDDO
         ENDDO
      ENDDO
C
      IF(IPRINT.GE.5) THEN
        CALL HEADER('Projection operator in AOC basis',-1)
        CALL OUTPUT(PTOT,1,LDM,1,LDM,LDM,LDM,-1,LUPRI)
      ENDIF
C
      CALL MEMREL('MK1CEX2',WORK,KSAVE,KSAVE,KFREE,LFREE)
      CALL QEXIT('MK1CEX2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mk1cex3 */
      SUBROUTINE MK1CEX3 (STOT,SLOC,PTOT,LDM,NTOT,NLOC,V,W,EIG,STOL,
     &                    IPRINT,WORK,KFREE,LFREE)
C*****************************************************************************
C     PURPOSE:
C        Make projection operator for this center.
C        We reduce the basis if the overlap matrix on this
C        center is nearly singular. This may happen as we work in
C        the untransformed non-orthogonal AO basis. STOL should be
C        consistent with the values used in LOWDIN !!
C
C        On input :
C           STOT : The complete overlap matrix
C           SLOC : The overlap matrix for the current center
C           LDM : Leading dimension (should be NBAST)
C           NTOT : The total number of functions to be projected
C                  (now all large or small functions in the molecule
C                   but generalization to domains is easy)
C           NLOC : The number of functions on the current center
C           STOL : Linear dependency treshold
C
C        On output :
C           PTOT : The updated projection operator matrix (added the
C                  rows for the current center)
C
C        NB : STOT and PTOT should be positioned such that the current
C             center is the first row !
C
C        Scratch arrays : V, W, EIG
C
C        Written by L.Visscher
C        Last revision: Aug 16 2001 - lv
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
      DIMENSION WORK(*)
      DIMENSION STOT(LDM,*),PTOT(LDM,*)
      DIMENSION SLOC(LDM,NLOC),V(LDM,NLOC),W(NLOC,NLOC),EIG(NLOC)
C
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Local overlap matrix',-1)
        CALL OUTPUT(SLOC,1,NLOC,1,NLOC,LDM,NLOC,-1,LUPRI)
      ENDIF
C
C     Form (Sloc)**-1/2 * U
C
      CALL LOWGEN(SLOC,NLOC,V,NV,LDM,STOL,IPRINT,EIG,WORK(KFREE),LFREE)
C
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Local Lowdin canonical matrix',-1)
        CALL OUTPUT(V,1,NLOC,1,NV,LDM,NLOC,-1,LUPRI)
      ENDIF
C
C     Form (Sloc)**-1/2 * U * U+ * (Sloc)**-1/2 = (Sloc)**-1
C
      CALL DGEMM('N','T',NLOC,NLOC,NV,D1,V,LDM,V,LDM,D0,W,NLOC)
C
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Local inverse overlap matrix',-1)
        CALL OUTPUT(W,1,NLOC,1,NLOC,NLOC,NLOC,-1,LUPRI)
      ENDIF
C
C     Make (Sloc)**-1  * Stot : the projection operator
C
      CALL DGEMM('N','N',NLOC,NTOT,NLOC,D1,W,NLOC,STOT,LDM,D0,PTOT,LDM)
C
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Local projection operator',-1)
        CALL OUTPUT(PTOT,1,NLOC,1,NTOT,LDM,NLOC,-1,LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PR1CEX1 */
      SUBROUTINE PR1CEX1 (FTRI,FIRST,FOCK,WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C     Correct for neglect of multicenter integrals in symmetry packed SO matrix
C
C     Written by L. Visscher, august 2001
C     Last revision: Aug 16, 2001
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbham.h"
C
      LOGICAL FIRST(NZ)
      DIMENSION FOCK(N2BBASX,NZ),FTRI(NNBBASX,NZ),
     &          WORK(*)
C
C     Make full matrix
C     ================
C
      DO IZ = 1,NZ
        IF(FIRST(IZ)) THEN
          CALL DZERO(FOCK(1,IZ),N2BBASX)
        ELSE
          IQ = IPQTOQ(IZ,0)
          IH = IHQMAT(IQ,1)
          IF(IH.EQ.1) THEN
            CALL DSPTSI(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DAPTGE(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ENDIF
        ENDIF
      ENDDO
C
C     Allocate work space to store the projection operator and
C     a work array.
C
      CALL MEMGET2('REAL','PROJ',KPROJ,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','BUF1',KBUF1,N2BBASX*4,WORK,KFREE,LFREE)
C
C     Get the projection operator
C
      CALL OPNFIL(LUPMAT,'DFPMAT','OLD','PR1CEX1')
      CALL READT(LUPMAT,N2BBASX,WORK(KPROJ))
      CLOSE(LUPMAT,STATUS='KEEP')
C
C     Add contributions from missing multicenter integrals
C     ====================================================
C
      IREP = 0
      CALL PR1CEX2 (NTBAS(0),IREP,FOCK,WORK(KPROJ),WORK(KBUF1),IPRINT)
C
C     Release the memory
C
      CALL MEMREL('PR1CEX1',WORK,KPROJ,KPROJ,KFREE,LFREE)
C
C     Pack matrix
C     ===========
C
      DO IZ = 1,NZ
        IF(.NOT.FIRST(IZ)) THEN
          IQ = IPQTOQ(IZ,0)
          IH = IHQMAT(IQ,1)
          IF(IH.EQ.1) THEN
            CALL DSITSP(NTBAS(0),FOCK(1,IZ),FTRI(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DGETAP(NTBAS(0),FOCK(1,IZ),FTRI(1,IZ))
          ENDIF
        ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PR1CEX2 */
      SUBROUTINE PR1CEX2 (NB,IREP,FOCK,PROJ,BUF1,IPRINT)
C***********************************************************************
C     Correct for neglect of multicenter integrals in SO matrix
C
C     A basis function |q> can be projected on the set of functions
C     |r> on center A using the projection operators P(r,q)
C     made in MK1CEX1.
C
C     We use this to correct the input matrix F(p,q) that contains
C     only the one-center contributions (F1) because the multicenter
C     contributions (F2) were neglected.
C
C     F(p,q) = F1(p,q) + F2(p,q)
C
C     FR(p,q) = F1(p,r) * P(r,q)
C     FL(p,q) = P+(p,r) * F1(r,q)
C
C     F(p,q) is approximated by F1(p,q) + 1/2 (FL(p,q) + FR(p,q))
C     We use the average of the left and right projected matrices to
C     preserve the permutation symmetry.
C
C     Written by L. Visscher, august 2001
C     Last revision: Aug 16, 2001
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
C
      PARAMETER (D0=0.D0, D1=1.D0, HALF=0.5D0)
      DIMENSION FOCK(NB,NB,NZ),PROJ(NB,NB),BUF1(NB,NB,4)
C
C     Transform Fock matrix to AO-basis
C     ==================================
C
      IF (IPRINT.GE.5) THEN
         CALL HEADER('PR1CEX2: Fock matrix in SO basis',-1)
         CALL PRQMAT(FOCK,NB,NB,NB,NB,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      IPAR = JBTOF(IREP,1)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL MTSOAO(FOCK(1,1,IPQ),BUF1(1,1,IZ),NB,IREPD,IPRINT)
      ENDDO
C
      DO IZ = 1, 4
         IF (IPRINT.GE.8) THEN
            CALL HEADER('PR1CEX2: Fock matrix in AO basis',-1)
            CALL OUTPUT(BUF1(1,1,IZ),1,NB,1,NB,NB,NB,-1,LUPRI)
         ENDIF
C
C        Construct the right projected fock matrix and store it in FOCK
C
         CALL DGEMM('N','N',NB,NB,NB,D1,BUF1(1,1,IZ),NB,PROJ,NB,
     &               D0,FOCK,NB)
C
C        Construct and add the left projected fock matrix
C
         CALL DGEMM('T','N',NB,NB,NB,D1,PROJ,NB,BUF1(1,1,IZ),NB,
     &               D1,FOCK,NB)
C
C        Add the correction to the Fock matrix
C
         IF (IPRINT.GE.10) THEN
            CALL HEADER('PR1CEX2 : Projection operator',-1)
            CALL OUTPUT(PROJ,1,NB,1,NB,NB,NB,-1,LUPRI)
            CALL HEADER('PR1CEX2: Correction to Fock matrix (*2)',-1)
            CALL OUTPUT(FOCK,1,NB,1,NB,NB,NB,-1,LUPRI)
         ENDIF
C
         NB2 = NB * NB
         CALL DAXPY (NB2,HALF,FOCK,1,BUF1(1,1,IZ),1)
C
         IF (IPRINT.GE.8) THEN
            CALL HEADER('PR1CEX2: Corrected Fock matrix',-1)
            CALL OUTPUT(BUF1(1,1,IZ),1,NB,1,NB,NB,NB,-1,LUPRI)
         ENDIF
      ENDDO
C
C     Transform Fock matrix to SO-basis
C     =================================
C
      CALL DZERO(FOCK,NB2*NZ)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL MTAOSO(BUF1(1,1,IZ),FOCK(1,1,IPQ),NB,IREPD,IPRINT)
      ENDDO
C
      IF (IPRINT.GE.5) THEN
         CALL HEADER('PR1CEX2: Corrected Fock matrix in SO basis',-1)
         CALL PRQMAT(FOCK,NB,NB,NB,NB,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PR1CEX3 */
      SUBROUTINE PR1CEX3 (NB,IREP,DENS,PROJ,BUF1,IPRINT)
C***********************************************************************
C     Correct for neglect of multicenter integrals in SO matrix
C     Note that this routines resembles PR1CEX2 that is used for
C     Fock matrix (contributions). This one is good for densities
C     where the transpose of P needs to be used and the transformation
C     to AO basis requires a different routine.
C
C     A basis function |q> can be projected on the set of functions
C     |r> on center A using the projection operators P(r,q)
C     made in MK1CEX1.
C
C     We use this to correct the input matrix D(q,p) to map multicenter
C     contributions to the calculated one-center integrals
C
C     D(q,p) = D1(q,p) + D2(q,p)
C
C     DR(q,p) = D1(q,r) * P+(r,p)
C     DL(q,p) = P(q,r) * D1(r,p)
C
C     D(q,p) is approximated by F1(q,p) + 1/2 (DL(q,p) + DR(q,p))
C     We use the average of the left and right projected matrices to
C     preserve the permutation symmetry.
C
C     Written by L. Visscher, august 2001
C     Last revision: Aug 17, 2001
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
C
      PARAMETER (D0=0.D0, D1=1.D0, HALF=0.5D0)
      DIMENSION DENS(NB,NB,NZ),PROJ(NB,NB),BUF1(NB,NB,4)
C
C     Transform density matrix to AO-basis
C     =====================================
C
      IF (IPRINT.GE.5) THEN
         CALL HEADER('PR1CEX3: Density matrix in SO basis',-1)
         CALL PRQMAT(DENS,NB,NB,NB,NB,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      IPAR = JBTOF(IREP,1)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL DTSOAO(DENS(1,1,IPQ),BUF1(1,1,IZ),NB,IREPD,IPRINT)
      ENDDO
C
      DO IZ = 1, 4
         IF (IPRINT.GE.8) THEN
            CALL HEADER('PR1CEX3: Density matrix in AO basis',-1)
            CALL OUTPUT(BUF1(1,1,IZ),1,NB,1,NB,NB,NB,-1,LUPRI)
         ENDIF
C
C        Construct the right projected fock matrix and store it in DENS
C
         CALL DGEMM('N','T',NB,NB,NB,D1,BUF1(1,1,IZ),NB,PROJ,NB,
     &               D0,DENS,NB)
C
C        Construct and add the left projected fock matrix
C
         CALL DGEMM('N','N',NB,NB,NB,D1,PROJ,NB,BUF1(1,1,IZ),NB,
     &               D1,DENS,NB)
C
C        Add the correction to the density matrix
C
         IF (IPRINT.GE.10) THEN
            CALL HEADER('PR1CEX3 : Projection operator',-1)
            CALL OUTPUT(PROJ,1,NB,1,NB,NB,NB,-1,LUPRI)
            CALL HEADER('PR1CEX3: Correction to density matrix (*2)',-1)
            CALL OUTPUT(DENS,1,NB,1,NB,NB,NB,-1,LUPRI)
         ENDIF
C
         NB2 = NB * NB
         CALL DAXPY (NB2,HALF,DENS,1,BUF1(1,1,IZ),1)
C
         IF (IPRINT.GE.8) THEN
            CALL HEADER('PR1CEX3: Corrected density matrix',-1)
            CALL OUTPUT(BUF1(1,1,IZ),1,NB,1,NB,NB,NB,-1,LUPRI)
         ENDIF
      ENDDO
C
C     Transform density matrix to SO-basis
C     ====================================
C
      CALL DZERO(DENS,NB2*NZ)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL DTAOSO(BUF1(1,1,IZ),DENS(1,1,IPQ),NB,IREPD,IPRINT)
      ENDDO
C
      IF (IPRINT.GE.5) THEN
         CALL HEADER('PR1CEX3: Corrected density matrix in SO basis',-1)
         CALL PRQMAT(DENS,NB,NB,NB,NB,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck glinsm */
      SUBROUTINE GLINSM(TMAT,VMAT,TBUF,EIG,WORK,LWORK)
C*****************************************************************************
C
C     This routine makes an additional transformation from the
C     matrix supplied by GMOTRA to a basis that is diagonal in
C     one component of angular momentum (j in general, l for spinfree
C     and Levy-Leblond). This can be used for systems that posess
C     linear symmetry.
C
C     Called from GMOTRA
C
C     Written by L. Visscher october/november 2003
C     Modified by T. Saue Jan 12 2007
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbham.h"
#include "dgroup.h"
      DIMENSION WORK(LWORK),TMAT(*),VMAT(*), TBUF(*), EIG(*)
      CALL QENTER('GLINSM')

#include "memint.h"
      IF (URKBAL) THEN
         WRITE (LUPRI,1000) "unrestricted kinetic balance"
         LINEAR = .FALSE.
         CALL QEXIT('GLINSM')
         RETURN
      ENDIF
C
      IF (SPINFR) THEN
         WRITE (LUPRI,1000) "spinfree"
         LINEAR = .FALSE.
         CALL QEXIT('GLINSM')
         RETURN
      ENDIF
C
      IF (ZORA) THEN
         WRITE (LUPRI,1000) "ZORA"
         LINEAR = .FALSE.
         CALL QEXIT('GLINSM')
         RETURN
      ENDIF
C
      CALL LINSYM(TMAT,VMAT,TBUF,EIG,WORK(KFREE),LFREE)
      SUB_BL = .TRUE.

      IF (IPRHAM.GE.1) THEN
        WRITE(LUPRI,'(/2X,A/)')
     &'* TMAT and Saomo were adapted for the linear symmetry in GLINSM.'
      ENDIF
      CALL FLSHFO(LUPRI)

      CALL QEXIT('GLINSM')
      RETURN
C
 1000 FORMAT (//" WARNING : linear symmetry is not implemented for ",
     &A," yet" /"           continue calculation in lower symmetry")
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gatmsm */
      SUBROUTINE GATMSM(TMAT,VMAT,TBUF,EIG,WORK,LWORK)
C*****************************************************************************
C
C     This routine makes an additional transformation from the
C     matrix supplied by GMOTRA to a basis that is diagonal in
C     one component of angular momentum (Kappa and mj).
C     This can be used for systems that posess atomic symmetry.
C
C     Routine based on GLINSM
C      
C     Called from GMOTRA
C
C     Written by A. Sunaga July 2018
C
C*****************************************************************************
      use mxcent_mod
#include "implicit.h"
#include "priunit.h"
#include "dcbham.h"
#include "dgroup.h"
#include "nuclei.h"
      DIMENSION WORK(LWORK),TMAT(*),VMAT(*), TBUF(*), EIG(*)
      CALL QENTER('GATMSM')

#include "memint.h"
      IF (URKBAL) THEN
         WRITE (LUPRI,1000) "unrestricted kinetic balance"
         LINEAR = .FALSE.
         ATOMIC = .FALSE.
         CALL QEXIT('GATMSM')
         RETURN
      ENDIF
C
      IF (SPINFR) THEN
         WRITE (LUPRI,1000) "spinfree"
         LINEAR = .FALSE.
         ATOMIC = .FALSE.
         CALL QEXIT('GATMSM')
         RETURN
      ENDIF
C
      IF (ZORA) THEN
         WRITE (LUPRI,1000) "ZORA"
         LINEAR = .FALSE.
         ATOMIC = .FALSE.
         CALL QEXIT('GATMSM')
         RETURN
      ENDIF
C
      IF(IPRHAM == 6 .and. atomic)THEN
        ATOMIC = .FALSE.
        IPRHAM = 5
      END IF
C
C The following part is needed to turn off ATOMIC
C in the BSSE calculation, that is
C   - NINT(CHARGE(IATOM)).EQ.0 ... ghost center
C   - NOORBT(IATOM).EQV. .FALSE. ... ghost with basis set
C However, when a single ghost center (NUCIND.EQ.1), 
C ATOMIC is not turned off.
C
      IF(NUCIND.GT.1)THEN
        DO IATOM = 1, NUCIND
          IF(NINT(CHARGE(IATOM)).EQ.0)THEN
          IF(NOORBT(IATOM).EQV. .FALSE.)THEN
            ATOMIC = .FALSE.         
          ENDIF
          ENDIF
        ENDDO
      ENDIF    
      IF(ATOMIC)THEN
       CALL ATMSYM(TMAT,VMAT,TBUF,EIG,WORK(KFREE),LFREE)
       IF (IPRHAM.GE.1) THEN
        WRITE(LUPRI,'(/2X,A/)')
     &'* TMAT and Saomo were adapted for the Atomic symmetry in GATMSM.'
       ENDIF
      ELSE
         CALL LINSYM(TMAT,VMAT,TBUF,EIG,WORK(KFREE),LFREE)         
      ENDIF
C
      SUB_BL = .TRUE.

      CALL FLSHFO(LUPRI)

      CALL QEXIT('GATMSM')
      RETURN
C
 1000 FORMAT (//" WARNING : atomic symmetry is not implemented for ",
     &A," yet" /"           continue calculation in lower symmetry")
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck linsym */
      SUBROUTINE LINSYM(TMAT,VMAT,TBUF,EIG,WORK,LWORK)
C*****************************************************************************
C
C     Construct matrix Jz(x,y,z) and transform to eigenvectors of this operator
C
C     Routine based on MODHAM and parts of the ZORA diagonalization routine
C
C     Called from GLINSM, LOWD2C
C
C     On input:  VMAT 4c AO-MO RKB transformation matrix
C
C     On output: TMAT - new 4c AO-MO RKB transf.matrix adapted for the linear symmetry
C                TBUF - MO to MO transformation matrix adapted for linear symmetry...
C
C     Written by L.Visscher, october 2003
C     Last changes: M.Ilias, Strasbourg, 2005
C                  M.Ilias, Tel Aviv, April 2009 - extended with QJACOBI
!                 MI, July 2014 - option to print out J_z matrices for
!                 diagonalization
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0, D1=1.D0, D2=2.D0)
C
C Used from COMMON blocks
C  CBIPAM:
C  DCBBAS:
C  DCBORB:
C  DGROUP:
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dummy.h"
C
      DIMENSION TMAT(*),TBUF(*),VMAT(*),EIG(*),WORK(LWORK)
      DIMENSION IQBF(4),IQW(4)
      LOGICAL :: DOJACO_SAVE
      character*100 :: hfile_format,hfile_name

      CALL QENTER('LINSYM')

#include "memint.h"
C
      CALL HEADER('Output from LINSYM',-1)
      NBRP = 4 / NZ
C
C     Construct mj matrix in AO basis
C
      CALL AMJMAT(TBUF,IPRHAM,WORK,LWORK)
C
C     Transform to MO-basis
C     =====================
C
C Note : It's tempting to use I2ORBT to find the start of a block, but this is not
C possible as this array is not yet set (is done in dirset). Took me half an hour to
C figure this one out (LV).
C
      IOFT = 1
      IOFV = 1
      DO 10 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 10
        IOFTB = I2BASX(I,I) + 1
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &            NORB(I),NORB(I),
     &            TBUF(IOFTB),NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &            TMAT(IOFT),NORB(I),NORB(I),1,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(I,0),NORB(I),1,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(I,0),NORB(I),1,IPQTOQ(1,0),
     &            WORK,LWORK,IPRHAM)

        IF (IPRHAM.GE.12) THEN

         CALL HEADER(
     &    'LINSYM: Entering AO->MO transformation matrix',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IOFV),NFBAS(I,0),NORB(I),
     &   NFBAS(I,0),NORB(I),1,IPQTOQ(1,0),LUPRI)

         CALL HEADER('LINSYM: MJ matrix in AO basis',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IOFTB),NFBAS(I,0),NFBAS(I,0),
     &   NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)

         CALL HEADER('LINSYM: MJ matrix in MO basis',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IOFT),NORB(I),NORB(I),
     &   NORB(I),NORB(I),1,IPQTOQ(1,0),LUPRI)

         CALL FLSHFO(LUPRI)

        ENDIF

        IOFT = IOFT + NORB(I)*NORB(I)
        IOFV = IOFV + NORB(I)*NFBAS(I,0)
   10 CONTINUE

C
C     Initialize TBUF to zero (we will fill in parts later and its easier to initialize
C     the rest here). Also initialize the array that counts the number of functions for
C     each Mj value.
C
      CALL DZERO (TBUF,N2TMT)

      CALL IZERO (NORB_SUB,MAX_SUB_BL*3*2)
      WRITE (LUPRI,1000)

!     stefan - june 2011: save LL block of <j_z>-eigenvectors (temporary) on file
!                         this block is needed to transform the 2c-LOWDIN matrix
!                         in the x2c-mmf approach (thus ensuring numerical stability by not re-doing this step).
      if(do4c2c)then
        open(13,file='JzMOLLb',status='replace',form='unformatted',
     &       access='sequential',action='readwrite',
     &       position='rewind')
      end if

C
      IOFT = 1
      DO 20 I = 1,NFSYM
C
        IF(NORB(I).EQ.0) GOTO 20
C
C       Set dimensions
C
        NL = NESH(I)
        NS = NPSH(I)
        NT = NORB(I)
        NLL = NL * NL
        NSS = NS * NS
        NTT = NT * NT
        IEIGP = IORB(I) + 1
        IEIGE = IEIGP + NS
C
C       Allocate memory for the diagonalization. We want to
C       have the LL and SS blocks in separate blocks.
C
        CALL MEMGET2('REAL','FLL',KFLL,NLL,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','VLL',KVLL,NLL,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','FSS',KFSS,NSS,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','VSS',KVSS,NSS,WORK,KFREE,LFREE)
C
C       Extract LL block of Jz matrix
C
        CALL EXTRSB (1,I,NBORB,4,NBRP,NBORB,4,NBRP,0,0,1,1,1,NZ,
     &               WORK(KFLL),NL,NL,1,
     &               TMAT(IOFT),NT,NT,NZ)
C     
C       Find eigenvalues of Jz
C
        IF (NL.GT.0) THEN
          IF (IPRHAM.GE.6) THEN
           CALL HEADER('LINSYM: Jz_LL to be diagonalized',-1)
           WRITE(LUPRI,'(2X,A,L1)') 'DOJACO=',DOJACO
           WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(WORK(KFLL),NL,NL,NL,NL,1,
     &          IPQTOQ(1,0),LUPRI)
          ENDIF

          IF (WRITE_Jz_MATRIX) THEN
! get the file name first, containing symmetry
            hfile_format="(A13,A7,I1,A1,I1)"
            write(hfile_name,hfile_format)
     &       "Jz_LL_matrix.","fermirp",I,"-",NFSYM
! open the unique file and write the matrix
             luhm=95
             open(luhm,file=trim(hfile_name),access="sequential",
     &         form="formatted",status="unknown")
               write(luhm,*) NL,1 ! this is pure real matrix (NZ=1)
               ic=0
               do ii=1,NL
               do jj=1,NL
                 write(luhm,*) ii,jj,WORK(KFLL+ic)
                 ic = ic + 1
               enddo
               enddo
             close(luhm,status="keep")
          ENDIF

          IF (.NOT.DOQJACO) THEN
            IF(DOJACO) THEN
              CALL RSJACO(NL,NL,NL,WORK(KFLL),EIG(IEIGE),1,1,0,
     &                  WORK(KVLL))
            ELSE
               CALL QDIAG(1,NL,WORK(KFLL),NL,NL,EIG(IEIGE),1,
     &                WORK(KVLL),NL,NL,WORK(KFREE),LFREE,IERR)
              IF (IERR.NE.0) THEN
                WRITE(LUPRI,'(2X,A)')
     &         'LINSYM: Erroneous QDIAG diagonalization of Jz_LL!'
               CALL QUIT(
     &        'LINSYM: Erroneous QDIAG diagonalization of Jz_LL !')
              ENDIF
            ENDIF
          ELSE
            CALL QJACOBI(WORK(KFLL),WORK(KVLL),NL,1,0,
     &                   IDUMMY,.TRUE.,IPRHAM)
            DO II=1,NL
              EIG(IEIGE+II-1)=WORK(KFLL+NL*(II-1)+II-1) ! store eigenvalues
            ENDDO
          ENDIF
          IF (IPRHAM.GE.5) THEN
           CALL HEADER('LINSYM: Jz_LL eigenvalues / eigenvectors',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,'/',NFSYM
           DO II = 1,NL
             WRITE(LUPRI,*) II,EIG(IEIGE+II-1)
           ENDDO
           CALL PRQMAT(WORK(KVLL),NL,NL,NL,NL,1,
     &          IPQTOQ(1,0),LUPRI)
           CALL FLSHFO(LUPRI)
          ENDIF
        ENDIF
C
C       Extract SS block of Jz matrix
C
        CALL EXTRSB (1,I,NBORB,4,NBRP,NBORB,4,NBRP,0,0,2,2,1,NZ,
     &               WORK(KFSS),NS,NS,1,
     &               TMAT(IOFT),NT,NT,NZ)
C
        IF (NS.GT.0) THEN
          IF (IPRHAM.GE.6) THEN
           CALL HEADER('LINSYM: Jz_SS to be diagonalized',-1)
           WRITE(LUPRI,'(2X,A,L1)') 'DOJACO=',DOJACO
           WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(WORK(KFSS),NS,NS,NS,NS,1,
     &          IPQTOQ(1,0),LUPRI)
          ENDIF

          IF (WRITE_Jz_MATRIX) THEN
! get the file name first, containing symmetry
            hfile_format="(A13,A7,I1,A1,I1)"
            write(hfile_name,hfile_format)
     &       "Jz_SS_matrix.","fermirp",I,"-",NFSYM
             print *,trim(hfile_name)
! open the unique file and write the matrix
             luhm=95
             open(luhm,file=trim(hfile_name),access="sequential",
     &         form="formatted",status="unknown")
               write(luhm,*) NL,1 ! this is pure real matrix (NZ=1)
               ic=0
               do ii=1,NS
               do jj=1,NS
                 write(luhm,*) ii,jj,WORK(KFSS+ic)
                 ic = ic + 1
               enddo
               enddo
             close(luhm,status="keep")
          ENDIF

          IF (.NOT.DOQJACO) THEN
            IF(DOJACO) THEN
              CALL RSJACO(NS,NS,NS,WORK(KFSS),EIG(IEIGP),1,1,0,
     &                  WORK(KVSS))
            ELSE
              CALL QDIAG(1,NS,WORK(KFSS),NS,NS,
     &                EIG(IEIGP),1,
     &                WORK(KVSS),NS,NS,
     &                WORK(KFREE),LFREE,IERR)
              IF (IERR.NE.0) THEN
                WRITE(LUPRI,'(2X,A)')
     &         'LINSYM: Erroneous QDIAG diagonalization of Jz_SS!'
               CALL QUIT(
     &        'LINSYM: Erroneous QDIAG diagonalization of Jz_SS !')
              ENDIF
            ENDIF
          ELSE
            CALL QJACOBI(WORK(KFSS),WORK(KVSS),NS,1,0,IDUMMY,
     &                   .TRUE.,IPRHAM)
            DO II=1,NS
              EIG(IEIGP+II-1)=WORK(KFSS+NS*(II-1)+II-1) ! store eigenvalues
            ENDDO
          ENDIF
          IF (IPRHAM.GE.5) THEN
           CALL HEADER('LINSYM: Jz_SS eigenvalues/eigenvectors',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &          '*** Fermion corep ',I,'/',NFSYM
           DO II = 1,NS
             WRITE(LUPRI,*) II,EIG(IEIGP+II-1)
           ENDDO
           CALL PRQMAT(WORK(KVSS),NS,NS,NS,NS,1,
     &          IPQTOQ(1,0),LUPRI)
           CALL FLSHFO(LUPRI)
          ENDIF
        ENDIF
C
C       Count the Mj-values.
C       We will either find a negative or a positive value for each Mj.
C
        MJMAX = 0
        DO IX = IEIGE, IEIGE+NL-1
           MJ = NINT(D2*EIG(IX))
           MJABS = (ABS(MJ)+1)/2
           IF (MJABS.GT.MAX_SUB_BL) THEN
             write(lupri,*) 'MAX_SUB_BL:',MAX_SUB_BL,
     &       ' actual MJ:',MJABS,' max value of MJ:',
     &       (ABS( NINT(D2*EIG(IEIGE+NL-1)))+1)/2
             CALL QUIT ("increase max_sub_bl")
           ENDIF
!        miro: for spin-free ensure the MJABS index not to be out-of-bounds
           IF (MJABS.GT.0) THEN
             ID_SUB_BL(MJABS,I)  = MJ
             NORB_SUB(MJABS,I,1) = NORB_SUB(MJABS,I,1) + 1
             NORB_SUB(MJABS,I,0) = NORB_SUB(MJABS,I,0) + 1
             MJMAX = MAX(MJABS,MJMAX)
           ENDIF
        ENDDO
C
        DO IX = IEIGP, IEIGP+NS-1
           MJ = NINT(D2*EIG(IX))
           MJABS = (ABS(MJ)+1)/2
           IF (MJABS.GT.MAX_SUB_BL) THEN
             write(lupri,*) 'MAX_SUB_BL:',MAX_SUB_BL,
     &        ' max value of MJ:',
     &       (ABS( NINT(D2*EIG(IEIGE+NL-1)))+1)/2
             CALL QUIT ("increase max_sub_bl")
           ENDIF
!        miro: for spin-free avoid MJABS index out-of-bounds
           IF (MJABS.GT.0) THEN
             ID_SUB_BL(MJABS,I)  = MJ
             NORB_SUB(MJABS,I,2) = NORB_SUB(MJABS,I,2) + 1
             NORB_SUB(MJABS,I,0) = NORB_SUB(MJABS,I,0) + 1
             MJMAX = MAX(MJABS,MJMAX)
           ENDIF
        ENDDO
C
        N_SUB_BL(I) = MJMAX
C
C       Write information to output
C
        DO MJABS = 1, MJMAX
           WRITE (LUPRI,1010) 3-2*I,ID_SUB_BL(MJABS,I),
     &                       (NORB_SUB(MJABS,I,IC),IC=0,2)
        ENDDO
C
C       We need to be careful now. In making NORB_SUB we chose to
C       use the absolute value of MJ because this gives the most logical
C       block structure. The eigenvectors and solutions are, however,
C       still ordered by the signed value. Reorder them first before
C       continuing.
C
        DO IX = IEIGP, IEIGP+NS-1
           EIG(IX) = ABS(EIG(IX))
        ENDDO
C
C       The SELBOS routine does just the reordering that we want,
C       use work(kfss) as scratch, and allocate a few other small
C       scratch arrays. The IBEIG part will not be used but we will
C       initialize it to prevent problems with uninitialized data while
C       sorting it.
C
        CALL MEMGET2('INTE','INDX' ,KINDX,NS,WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IDUM' ,KIDUM,NS,WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IBEIG',KIBEI,NS,WORK,KFREE,LFREE)
        MJ = 0
        CALL ICOPY (NS,MJ,0,WORK(KIBEI),1)
        CALL SELBOS (I,NS,WORK(KINDX),EIG(IEIGP),
     &               WORK(KIBEI),WORK(KVSS),
     &               WORK(KFSS),WORK(KIDUM),IPRHAM)
        CALL MEMREL('LINSYM',WORK,KINDX,KINDX,KFREE,LFREE)
C
C       Same thing for the positive energy solutions
C
        DO IX = IEIGE, IEIGE+NL-1
           EIG(IX) = ABS(EIG(IX))
        ENDDO
        CALL MEMGET2('INTE','INDX' ,KINDX,NL,WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IDUM' ,KIDUM,NL,WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IBEIG',KIBEI,NL,WORK,KFREE,LFREE)
        MJ = 0
        CALL ICOPY (NL,MJ,0,WORK(KIBEI),1)
        CALL SELBOS (I,NL,WORK(KINDX),EIG(IEIGE),
     &               WORK(KIBEI),WORK(KVLL),
     &               WORK(KFLL),WORK(KIDUM),IPRHAM)
        CALL MEMREL('LINSYM',WORK,KINDX,KINDX,KFREE,LFREE)
C
C       Put the MO vectors in TBUF
C
C       We now have the solutions ordered on on the right index
C       by absolute mj-value and at the left index on boson irreps.
C       Keep this ordering when inserting into the full matrix so
C       that the new transformation matrix (generated below) will
C       give an ordering on MJ values.
C

        CALL DZERO(TBUF(IOFT),NT*NT)

        CALL EXTRSB (2,I,NBORB,4,NBRP,
     &               NORB_SUB,MAX_SUB_BL,N_SUB_BL(I),0,0,1,1,1,1,
     &               WORK(KVLL),NL,NL,1,
     &               TBUF(IOFT),NT,NT,1)

        if(do4c2c)then ! x2c-mmf approach
          call writt(13,NL**2,WORK(KVLL))
        end if

        IF (IPRHAM.GE.6) THEN
         CALL HEADER('LINSYM: MO vectors in TBUF,LL',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IOFT),NT,NT,NT,NT,1,
     &   IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)

         CALL HEADER('LINSYM: MO vectors in buffer space,LL',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(WORK(KVLL),NL,NL,NL,NL,1,
     &   IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
        ENDIF

        CALL EXTRSB (2,I,NBORB,4,NBRP,
     &               NORB_SUB,MAX_SUB_BL,N_SUB_BL(I),0,0,2,2,1,1,
     &               WORK(KVSS),NS,NS,1,
     &               TBUF(IOFT),NT,NT,1)

        IF (IPRHAM.GE.6) THEN
         CALL HEADER('LINSYM: MO vectors in TBUF,SS',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IOFT),NT,NT,NT,NT,1,
     &   IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
        ENDIF
C
C       Done; release the extra memory used for diagonalization
C
        CALL MEMREL('LINSYM',WORK,KFLL,KFLL,KFREE,LFREE)
C
C       Update offset
C
        IOFT = IOFT + NORB(I)*NORB(I)
C
   20 CONTINUE

      if(do4c2c)then ! x2c-mmf approach
        close(13,status='keep')
      end if
C
C     Generate transformation matrix in non-orthogonal basis
C     ======================================================
C
      IOFT = 1
      IOFU = 1
      IOFV = 1
      DO 30 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 30
        IF (IPRHAM.GE.5) THEN
         CALL HEADER(
     &'LINSYM: Transformation matrix in orthon. basis (MO->MO),TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IOFU),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),1,IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
        ENDIF
C
C       Backtransform all solutions
C       ===========================
C
        CALL BCKTRA(TMAT(IOFT),NFBAS(I,0),NORB(I),
     &              TBUF(IOFU),NORB(I),NORB(I),
     &              NORB(I),1,
     &              NORB(I),1,NFBAS(I,0),
     &              VMAT(IOFV),NFBAS(I,0),NORB(I),1,
     &              IPRHAM)

        IF (IPRHAM.GE.5) THEN
         CALL HEADER(
     &'LINSYM: Final (symm.blocked) transformation'//
     &' matrix (AO->MO),TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IOFT),NFBAS(I,0),NORB(I),
     &         NFBAS(I,0),NORB(I),1,IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
        ENDIF

        IOFT = IOFT + NFBAS(I,0)*NORB(I)
        IOFU = IOFU + NORB(I)*NORB(I)
        IOFV = IOFV + NORB(I)*NFBAS(I,0)
 30   CONTINUE
C.....Also initialize NTMO_SUB which will be the number of orbitals in each subblock in orthonormal basis
      CALL ICOPY(MAX_SUB_BL*3*2,NORB_SUB,1,NTMO_SUB,1)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('LINSYM')
      RETURN
 1000 FORMAT (/' Parity   MJ  Functions(total) Functions(LC)'
     &         ' Functions(SC)')
 1010 FORMAT (T5,I2,T9,I3,'/2',T20,I5,T33,I5,T46,I5)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck atmsym */
      SUBROUTINE ATMSYM(TMAT,VMAT,TBUF,EIG,WORK,LWORK)
C*****************************************************************************
C
C     Generate transformation to atomic supersymmetry
C     This will lead to orbitals with well-defined kappa and m_j.
C     The routine first builds the matrix representation of the kappa operator
C     in orthonormal basis and diagonalizes it in a second step.
C     The transformation to orthonormal basis is then updated with the
C     transformation to kappa blocks.
C
C     In the ensuing step, for each kappa block, the matrix representation of the j_z operator
C     is diagonalized. The transformation to orthonormal basis is then updated with the
C     transformation to (kappa, mj) - blocks
C
C     The kappa quantum number can be expressed as
C       kappa = a(j+1/2), a = 2(l-j)
C     The degeneracy of mj-values is therefore
C       2j+1 = 2|kappa|
C     However, since we do Kramers-restricted calculations, we only keep one value of mj for each |mj|
C     and so we need to reserve |kappa| mj-values for each kappa.      
C
C     Note that large component functions partner with small component functions of opposite sign of kappa.
C
C     In the normal atomic case, that is, with inversion symmetry (NFSYM=2), the kappa values will be ordered as
C         gerade: -1,+2,-3,...
C       ungerade: +1,-2,+3,-4 ,...
C     and so for each inversion symmetry the offset to the subblocks associated with a given kappa will be |kappa|(|kappa|-1)/2.
C
C     Atomic supersymmetry is also activated in the presence of one or more ghost atoms lowering symmetry to linear.
C     In that case, inversion symmetry is lost, and kappa values will be ordered as
C       -1,+1,-2,+2,-3,+3,-4....
C     Offsets are then
C       |kappa|(|kappa|-1) for negative kappa
C       |kappa|x|kappa|    for positive kappa    
C
C     With atomic supersymmetry each symmetry subblock are characterized by a double index. They will be packed as a single superindex.
C     We recall that a superindex IJ refers to upper triangular matrices, that is
C     (1,1)  (1,2) (1,3) (1,4) ..
C            (2,2) (2,3) (2,4) ..       
C                  (3,3) (3,4) ..
C                        (4,4) ..
C                              ..
C     From this it will be seen that a superindex may be expressed as
C       IJ = J(J-1)/2 + I
C
C     In the atomic case we have
C       |kappa|   1    2    3    4
C                1/2  1/2  1/2  1/2
C                    -3/2 -3/2 -3/2      
C                          5/2  5/2
C                              -7/2
C
C     In the above matrix, |kappa| is the column index, whereas rows relate to m_j.
C     It will be noted that m_j has alternating sign; this is a consequence of the quaternion
C     symmetry scheme of DIRAC. We may encode m_j with respect to the row index as
C
C             m_j = (1/2)(2I-1)(-1)^(I+1)
C
C     For kappa, we note that in the usual atomic case, the kappa values are distributed as
C       IFRP=1 (gerade  ):  -1,+2,-3,...
C       IFRP=2 (ungerade):  +1,-2,+3,...
C
C     Routine based on LINSYM
C     
C     Called from GLINSM, LOWD2C
C
C     On input:  VMAT 4c AO-MO RKB transformation matrix
C
C     On output: TMAT - new 4c AO-MO RKB transf.matrix adapted for atomic supersymmetry
C                TBUF - MO to MO transformation matrix adapted for atomic supersymmetry
C
C     Written by A.Sunaga, July 2018
C     Bugfixes by S. Knecht, May 2021
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0, D1=1.D0, D2=2.D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dummy.h"
#include "dcbdhf.h"            
#include "dcbibt.h"      
C
      DIMENSION TMAT(*),TBUF(*),VMAT(*),EIG(*),WORK(LWORK)
      DIMENSION IQBF(4),IQW(4),KPMAX(2)
      REAL*8, allocatable :: TMPA(:),TMPB(:)
      LOGICAL :: DOJACO_SAVE
      character*100 :: hfile_format,hfile_name

      CALL QENTER('ATMSYM')

#include "memint.h"
C

      CALL HEADER('Output from ATMSYM',-1)
      NBRP = 4 / NZ
C     
C     Construct kappa matrix in AO basis
C
      CALL AMKMAT(TBUF,IPRHAM,WORK,LWORK)
C
C     Transform to MO-basis
C     =====================
C
      IOFT = 1
      IOFV = 1
      N2TMT = 0
      N2TMOTQ = 0
      DO 10 IFRP = 1,NFSYM
        IF(NORB(IFRP).EQ.0) GOTO 10
        I2TMT(IFRP)  = N2TMT
        N2TMT        = N2TMT + NFBAS(IFRP,0)*NORB(IFRP)*NZT
        I2TMOT(IFRP) = N2TMOTQ
        N2TMO(IFRP)  = NTMO(IFRP)*NTMO(IFRP)
        N2TMOTQ      = N2TMOTQ + N2TMO(IFRP)*NZ
C        
        IOFTB = I2BASX(IFRP,IFRP) + 1        
        CALL QTRANS('AOMO','S',D0,NFBAS(IFRP,0),NFBAS(IFRP,0),
     &            NORB(IFRP),NORB(IFRP),
     &            TBUF(IOFTB),NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &            TMAT(IOFT),NORB(IFRP),NORB(IFRP),1,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(IFRP,0),NORB(IFRP),1,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(IFRP,0),NORB(IFRP),1,IPQTOQ(1,0),
     &            WORK,LWORK,IPRHAM)

        IF (IPRHAM.GE.10) THEN

         CALL HEADER(
     &    'ATMSYM: Entering AO->MO transformation matrix',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',IFRP,'/',NFSYM
         CALL PRQMAT(VMAT(IOFV),NFBAS(IFRP,0),NORB(IFRP),
     &   NFBAS(IFRP,0),NORB(IFRP),1,IPQTOQ(1,0),LUPRI)

         CALL HEADER('ATMSYM: K matrix in AO basis',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',IFRP,'/',NFSYM
         CALL PRQMAT(TBUF(IOFTB),NFBAS(IFRP,0),NFBAS(IFRP,0),
     &   NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)

         CALL HEADER('ATMSYM: K matrix in MO basis',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',IFRP,'/',NFSYM
         CALL PRQMAT(TMAT(IOFT),NORB(IFRP),NORB(IFRP),
     &   NORB(IFRP),NORB(IFRP),1,IPQTOQ(1,0),LUPRI)

         CALL FLSHFO(LUPRI)

        ENDIF

        IOFT = IOFT + NORB(IFRP)*NORB(IFRP)
        IOFV = IOFV + NORB(IFRP)*NFBAS(IFRP,0)
   10 CONTINUE

C
C     Initialize TBUF to zero (we will fill in parts later and its easier to initialize
C     the rest here). Also initialize the array that counts the number of functions for
C     each k value.
C
      CALL DZERO (TBUF,N2TMT)
      CALL IZERO (ID_SUB_BL,MAX_SUB_BL*2)
      CALL IZERO (NORB_SUB,MAX_SUB_BL*3*2)
      IF (IPRHAM.GE.2) THEN
        WRITE (LUPRI,1000)
      ENDIF
!     stefan - june 2011: save LL block of <j_z>-eigenvectors (temporary) on file
!                         this block is needed to transform the 2c-LOWDIN matrix
!                         in the x2c-mmf approach (thus ensuring numerical stability by not re-doing this step).
      if(do4c2c)then
        open(13,file='KMJMOLLb',status='replace',form='unformatted',
     &       access='sequential',action='readwrite',
     &       position='rewind')
      end if
      
C
      IOFT = 1
      CALL IZERO (N_SUB_BL,2)
      DO 20 IFRP = 1,NFSYM
        IF(NORB(IFRP).EQ.0) GOTO 20
C
C       Set dimensions
C
        NL     = NESH(IFRP)
        NS     = NPSH(IFRP)
        NT     = NORB(IFRP)
        NLL    = NL * NL
        NSS    = NS * NS
        NTT    = NT * NT
        IEIGP  = IORB(IFRP) + 1
        IEIGE  = IEIGP + NS
        KPMAXL = 0
        KPMAXS = 0

C
C       Allocate memory for the diagonalization. We want to
C       have the LL and SS blocks in separate blocks.
C
        CALL MEMGET('REAL',KFLL,NLL,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFSS,NSS,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KVLL,NLL,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KVSS,NSS,WORK,KFREE,LFREE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       LL Extract and Diagonalization
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       Extract LL block of Kappa matrix
C
        CALL EXTRSB(1,IFRP,NBORB,4,NBRP,NBORB,4,NBRP,0,0,1,1,1,NZ,
     &              WORK(KFLL),NL,NL,1,
     &              TMAT(IOFT),NT,NT,NZ)
C
C       Find eigenvalues of Kappa
C
        IF (NL.GT.0) THEN
          IF (IPRHAM.GE.6) THEN
            CALL HEADER('ATMSYM: Kappa_LL to be diagonalized',-1)
            WRITE(LUPRI,'(2X,A,L1)') 'DOJACO=',DOJACO
            WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',IFRP,'/',NFSYM
            CALL PRQMAT(WORK(KFLL),NL,NL,NL,NL,1,
     &          IPQTOQ(1,0),LUPRI)
          ENDIF
          IF (WRITE_K_MATRIX) THEN
! get the file name first, containing symmetry
            hfile_format="(A13,A7,I1,A1,I1)"
            write(hfile_name,hfile_format)
     &       "Kappa_LL_matrix.","fermirp",IFRP,"-",NFSYM
! open the unique file and write the matrix
            luhm=95
            open(luhm,file=trim(hfile_name),access="sequential",
     &          form="formatted",status="unknown")
            write(luhm,*) NL,1 ! this is pure real matrix (NZ=1)
            ic=0
            do ii=1,NL
              do jj=1,NL
                write(luhm,*) ii,jj,WORK(KFLL+ic)
                ic = ic + 1
              enddo
            enddo
            close(luhm,status="keep")
          ENDIF

          IF (.NOT.DOQJACO) THEN
            IF(DOJACO) THEN
              CALL RSJACO(NL,NL,NL,WORK(KFLL),EIG(IEIGE),1,1,0,
     &                  WORK(KVLL))
            ELSE
              CALL QDIAG(1,NL,WORK(KFLL),NL,NL,
     &                EIG(IEIGE),1,
     &                WORK(KVLL),NL,NL,
     &                WORK(KFREE),LFREE,IERR)
              IF (IERR.NE.0) THEN
                WRITE(LUPRI,'(2X,A)')
     &         'ATMSYM: Errorneous QDIAG diagonalization of Kappa_LL!'
               CALL QUIT(
     &        'ATMSYM: Errorneous QDIAG diagonalization of Kappa_LL !')
              ENDIF
            ENDIF
          ELSE
            CALL QJACOBI(WORK(KFLL),WORK(KVLL),NL,1,0,
     &                   IDUMMY,.TRUE.,IPRHAM)
            DO II=1,NL
              EIG(IEIGE+II-1)=WORK(KFLL+NL*(II-1)+II-1) ! store eigenvalues
            ENDDO
          ENDIF
          KPMAXL = NINT(MAX(ABS(EIG(IEIGE)),ABS(EIG(IEIGE+NL-1))))
          IF (IPRHAM.GE.3) THEN
            CALL HEADER('ATMSYM: Kappa_LL eigenvalues',-1)
            WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &       '*** Fermion corep ',IFRP,'/',NFSYM
            WRITE(LUPRI,'(A,I2)') ' * Max ABS(kappa_large):', KPMAXL
            DO II = 1,NL
              WRITE(LUPRI,*) II,EIG(IEIGE+II-1)
            ENDDO
            CALL FLSHFO(LUPRI)
          ENDIF
          IF (IPRHAM.GE.5) THEN
            CALL HEADER('ATMSYM: Kappa_LL eigenvectors',-1)
            WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',IFRP,'/',NFSYM
            CALL PRQMAT(WORK(KVLL),NL,NL,NL,NL,1,
     &                  IPQTOQ(1,0),LUPRI)
            CALL FLSHFO(LUPRI)
          ENDIF
        ENDIF  ! IF (NL.GT.0) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SS Extract and Diagonalization
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C       Extract SS block of Kappa matrix
C
        CALL EXTRSB (1,IFRP,NBORB,4,NBRP,NBORB,4,NBRP,0,0,2,2,1,NZ,
     &               WORK(KFSS),NS,NS,1,
     &               TMAT(IOFT),NT,NT,NZ)
C
        IF (NS.GT.0) THEN
          IF (IPRHAM.GE.6) THEN
            CALL HEADER('ATMSYM: Kappa_SS to be diagonalized',-1)
            WRITE(LUPRI,'(2X,A,L1)') 'DOJACO=',DOJACO
            WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',IFRP,'/',NFSYM
            CALL PRQMAT(WORK(KFSS),NS,NS,NS,NS,1,
     &                  IPQTOQ(1,0),LUPRI)
 9        ENDIF

          IF (WRITE_K_MATRIX) THEN
! get the file name first, containing symmetry
            hfile_format="(A13,A7,I1,A1,I1)"
            write(hfile_name,hfile_format)
     &      "Kappa_SS_matrix.","fermirp",IFRP,"-",NFSYM
            print *,trim(hfile_name)
! open the unique file and write the matrix
            luhm=95
            open(luhm,file=trim(hfile_name),access="sequential",
     &        form="formatted",status="unknown")
            write(luhm,*) NL,1 ! this is pure real matrix (NZ=1)
            ic=0
            do ii=1,NS
              do jj=1,NS
                write(luhm,*) ii,jj,WORK(KFSS+ic)
                ic = ic + 1
              enddo
            enddo
            close(luhm,status="keep")
          ENDIF
          IF (.NOT.DOQJACO) THEN
            IF(DOJACO) THEN
              CALL RSJACO(NS,NS,NS,WORK(KFSS),EIG(IEIGP),1,1,0,
     &                  WORK(KVSS))
            ELSE
              CALL QDIAG(1,NS,WORK(KFSS),NS,NS,
     &                EIG(IEIGP),1,
     &                WORK(KVSS),NS,NS,
     &              WORK(KFREE),LFREE,IERR)
              IF (IERR.NE.0) THEN
                WRITE(LUPRI,'(2X,A)')
     &         'ATMSYM: Errorneous QDIAG diagonalization of Kappa_SS!'
               CALL QUIT(
     &        'ATMSYM: Errorneous QDIAG diagonalization of Kappa_SS !')
              ENDIF
            ENDIF
          ELSE
            CALL QJACOBI(WORK(KFSS),WORK(KVSS),NS,1,0,IDUMMY,
     &                   .TRUE.,IPRHAM)
            DO II=1,NS
              EIG(IEIGP+II-1)=  WORK(KFSS+NS*(II-1)+II-1) ! store eigenvalues
            ENDDO
         ENDIF
          KPMAXS=NINT(MAX(ABS(EIG(IEIGP)),ABS(EIG(IEIGP+NS-1))))
          IF (IPRHAM.GE.3) THEN
            CALL HEADER('ATMSYM: Kappa_SS eigenvalues',-1)
            WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &           '*** Fermion corep ',IFRP,'/',NFSYM
            WRITE(LUPRI,'(A,I2)') ' * Max ABS(kappa_small):', KPMAXS
            DO II = 1,NS
              WRITE(LUPRI,*) II,EIG(IEIGP+II-1)
            ENDDO
            CALL FLSHFO(LUPRI)
          ENDIF
          IF (IPRHAM.GE.5) THEN
            CALL HEADER('ATMSYM: Kappa_SS eigenvectors',-1)
            WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &           '*** Fermion corep ',IFRP,'/',NFSYM
            CALL PRQMAT(WORK(KVSS),NS,NS,NS,NS,1,
     &           IPQTOQ(1,0),LUPRI)
            CALL FLSHFO(LUPRI)
          ENDIF
        ENDIF  !IF (NS.GT.0) THEN
C        
C       Find maximal kappa value
C       This determines the number of sub-blocks for this fermion irrep
C       Note that we already anticipating the further division into mj-blocks for each kappa        
C
        KPMAX(IFRP) = MAX(KPMAXL,KPMAXS)
        IF(NFSYM.EQ.2)THEN
          MAXSUB = KPMAX(IFRP)*(KPMAX(IFRP)+1)/2
        ELSE
          MAXSUB = KPMAX(IFRP)*KPMAX(IFRP)
        ENDIF         
        IF (MAXSUB.GT.MAX_SUB_BL) THEN
          write(lupri,*) 'MAX_SUB_BL:',MAX_SUB_BL,
     &       ' max value of Kappa:',KPMAX(IFRP),
     &       ' increase max_sub_bl to ', MAXSUB        
          CALL QUIT ("increase max_sub_bl")
        ENDIF      
        N_SUB_BL(IFRP) = MAXSUB
C
C       We now have kappa eigenvalues and corresponding eigenvectors
C       However, since the diagonalization routines orders on
C       increasing eigenvalue, we shall use SELBOS to get the right ordering
C     Note how eigenvalues are modified to get the desired ordering
C        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       LL counting and sorting
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
C      (KPABS,KPADD)=(1,1),(2,2),(3,4),(4,7)
C       Discrete KPASS is needed to allocate the corresponding Mj           
        KPBUF = 0
        KPADD = 0
        DO IX = IEIGE, IEIGE+NL-1
           KP = NINT(EIG(IX))
           IF (EIG(IX).GT.0) THEN
              EIG(IX) = EIG(IX)**D2 + EIG(IX)
           ELSE
              EIG(IX) = EIG(IX)**D2
           ENDIF
           IF(KP.NE.KPBUF) THEN
              KPBUF = KP
              KPABS = ABS(KP)
              IF(NFSYM.eq.2)THEN 
                 KPADD = KPABS*(KPABS-1)/2 + 1
              ELSE
                 KPADD = KPABS*(KPABS-1) + 1
                 IF(KP.GT.0) KPADD = KPADD + KP
              ENDIF
              ID_SUB_BL(KPADD,IFRP) = KP
           ENDIF
           NORB_SUB(KPADD,IFRP,1) = NORB_SUB(KPADD,IFRP,1) + 1
           NORB_SUB(KPADD,IFRP,0) = NORB_SUB(KPADD,IFRP,0) + 1
        ENDDO
        CALL MEMGET('INTE',KINDX,NL,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIDUM,NL,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBEI,NL,WORK,KFREE,LFREE)
        KP = 0
        CALL ICOPY (NL,KP,0,WORK(KIBEI),1) ! IBEIG is formally not needed, but initialized to avoid problems
        CALL SELBOS (IFRP,NL,WORK(KINDX),EIG(IEIGE),
     &               WORK(KIBEI),WORK(KVLL),
     &               WORK(KFLL),WORK(KIDUM),IPRHAM)
        CALL MEMREL('ATMSYM',WORK,KINDX,KINDX,KFREE,LFREE)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SS counting and sorting
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        KPBUF = 0
        KPADD = 0
        DO IX = IEIGP, IEIGP+NS-1
           KP = -NINT(EIG(IX))     ! Note the opposite sign of kappa 
           IF (EIG(IX).GT.0) THEN
              EIG(IX) = EIG(IX)**D2
           ELSE
              EIG(IX) = EIG(IX)**D2 - EIG(IX) ! Notice sign change to get right ordering
           ENDIF
           IF(KP.NE.KPBUF) THEN
              KPBUF = KP
              KPABS = ABS(KP)
              IF(NFSYM.eq.2)THEN 
                 KPADD = KPABS*(KPABS-1)/2 + 1
              ELSE
                 KPADD = KPABS*(KPABS-1) + 1
                 IF(KP.GT.0) KPADD = KPADD + KP
              ENDIF  
              ID_SUB_BL(KPADD,IFRP) = KP
           ENDIF
           NORB_SUB(KPADD,IFRP,2) = NORB_SUB(KPADD,IFRP,2) + 1
           NORB_SUB(KPADD,IFRP,0) = NORB_SUB(KPADD,IFRP,0) + 1
        ENDDO
        CALL MEMGET('INTE',KINDX,NS,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIDUM,NS,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBEI,NS,WORK,KFREE,LFREE)
        KP = 0
        CALL ICOPY (NS,KP,0,WORK(KIBEI),1) ! IBEIG is formally not needed, but initialized to avoid problems
        CALL SELBOS (IFRP,NS,WORK(KINDX),EIG(IEIGP),
     &               WORK(KIBEI),WORK(KVSS),
     &               WORK(KFSS),WORK(KIDUM),IPRHAM)
        CALL MEMREL('ATMSYM',WORK,KINDX,KINDX,KFREE,LFREE)
C
C     Write information to output
C
      IF (IPRHAM.GE.2) THEN
        DO KPADD = 1, MAXSUB
           WRITE (LUPRI,1010) 3-2*IFRP,ID_SUB_BL(KPADD,IFRP),
     &                       (NORB_SUB(KPADD,IFRP,IC),IC=0,2)
        ENDDO
      ENDIF  
C            
C       Put the MO vectors in TBUF
C
C       We now have the solutions ordered on on the right index
C       by absolute mj-value and at the left index on boson irreps.
C       Keep this ordering when inserting into the full matrix so
C       that the new transformation matrix (generated below) will
C       give an ordering on KP values.
C        
        CALL DZERO(TBUF(IOFT),NT*NT)

        CALL EXTRSB (2,IFRP,NBORB,4,NBRP,
     &               NORB_SUB,MAX_SUB_BL,N_SUB_BL(IFRP),0,0,1,1,1,1,
     &               WORK(KVLL),NL,NL,1,
     &               TBUF(IOFT),NT,NT,1)
C
        IF (IPRHAM.GE.5) THEN
         CALL HEADER('ATMSYM: MO vectors in TBUF,LL',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',IFRP,'/',NFSYM
         CALL PRQMAT(TBUF(IOFT),NT,NT,NT,NT,1,
     &   IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)

         CALL HEADER('ATMSYM: MO vectors in buffer space,LL',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',IFRP,'/',NFSYM
         CALL PRQMAT(WORK(KVLL),NL,NL,NL,NL,1,
     &   IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
        ENDIF

        CALL EXTRSB (2,IFRP,NBORB,4,NBRP,
     &               NORB_SUB,MAX_SUB_BL,N_SUB_BL(IFRP),0,0,2,2,1,1,
     &               WORK(KVSS),NS,NS,1,
     &               TBUF(IOFT),NT,NT,1)

        IF (IPRHAM.GE.5) THEN
         CALL HEADER('ATMSYM: MO vectors in TBUF,SS',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',IFRP,'/',NFSYM
         CALL PRQMAT(TBUF(IOFT),NT,NT,NT,NT,1,
     &   IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
        ENDIF
C
C       Done; release the extra memory used for diagonalization
C
        CALL MEMREL('ATMSYM',WORK,KFLL,KFLL,KFREE,LFREE)
C
C       Update offset
C
        IOFT = IOFT + NORB(IFRP)*NORB(IFRP)
C        
 20   CONTINUE ! DO 20 IFRP = 1,NFSYM
C         

C
C     Generate transformation matrix in non-orthogonal basis
C     ======================================================
C
      IOFT = 1
      IOFU = 1
      IOFV = 1
      DO 30 IFRP = 1,NFSYM
        IF(NORB(IFRP).EQ.0) GOTO 30
        IF (IPRHAM.GE.5) THEN
         CALL HEADER(
     &'ATMSYM: Transformation matrix in orthon. basis (MO->MO),TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',IFRP,'/',NFSYM
         CALL PRQMAT(TBUF(IOFU),NORB(IFRP),NORB(IFRP),NORB(IFRP),
     &                 NORB(IFRP),1,IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
        ENDIF
C
C       Backtransform all solutions
C       ===========================
C
        CALL BCKTRA(TMAT(IOFT),NFBAS(IFRP,0),NORB(IFRP),
     &              TBUF(IOFU),NORB(IFRP),NORB(IFRP),
     &              NORB(IFRP),1,
     &              NORB(IFRP),1,NFBAS(IFRP,0),
     &              VMAT(IOFV),NFBAS(IFRP,0),NORB(IFRP),1,
     &              IPRHAM)
C
        IF (IPRHAM.GE.5) THEN
         CALL HEADER(
     &'ATMSYM: TMAT (symm.blocked) transformation'//
     &' matrix (AO->MO) aka VU/kappa,TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',IFRP,'/',NFSYM
         CALL PRQMAT(TMAT(IOFT),NFBAS(IFRP,0),NORB(IFRP),
     &         NFBAS(IFRP,0),NORB(IFRP),1,IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
        ENDIF

        IOFT = IOFT + NFBAS(IFRP,0)*NORB(IFRP)
        IOFU = IOFU + NORB(IFRP)*NORB(IFRP)
        IOFV = IOFV + NORB(IFRP)*NFBAS(IFRP,0)
 30   CONTINUE !  DO 30 IFRP = 1,NFSYM

      CALL ICOPY(MAX_SUB_BL*3*2,NORB_SUB,1,NTMO_SUB,1)  
      !> copy kappa-adapted 4c AO->atomic MO transformation matrix from TMAT to VMAT
      !> so that VMAT becomes VU, i.e. V => VU
      CALL DCOPY(N2TMT,TMAT,1,VMAT,1)
      IF(IPRHAM.GE.5) THEN
            CALL HEADER(
     &    'dcopy: AFTER LINSYM final transformation'//
     &    ' matrix 4c AO->LINSYM MO, TMAT(output)',-1)
        IMAT = 1
        IBUF = 1
        DO IFRP = 1,NFSYM
            IF(NORB(IFRP).GT.0) THEN
                WRITE(LUPRI,'(A,I1,A,I1)')
     &      '* Fermion ircop no.',IFRP,'/',NFSYM
                CALL PRQMAT(TMAT(IMAT),NFBAS(IFRP,0),NORB(IFRP),
     &                  NFBAS(IFRP,0),NORB(IFRP),1,IPQTOQ(1,0),LUPRI)
                WRITE(LUPRI,'(A,I1,A,I1)')
     &      '* Fermion ircop no.',IFRP,'/',NFSYM
                CALL PRQMAT(VMAT(IMAT),NFBAS(IFRP,0),NORB(IFRP),
     &                  NFBAS(IFRP,0),NORB(IFRP),1,IPQTOQ(1,0),LUPRI)
                WRITE(LUPRI,'(/5X,A/)')
                IMAT = IMAT + NFBAS(IFRP,0)*NORB(IFRP)*1
                IBUF = IBUF + NORB(IFRP)*NORB(IFRP)*1
            ENDIF
         ENDDO
      ENDIF
    
C       
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C.....Here we start MJ-part
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Construct mj matrix in AO basis
      allocate(TMPA(NTBAS(0)*NTBAS(0)))
      CALL DZERO(TMPA,NTBAS(0)*NTBAS(0))      
      CALL AMJMAT(TMPA,IPRHAM,WORK,LWORK) 
C
C     Output    
C  
      WRITE (LUPRI,'(/A)')
     &' Parity  Kappa  MJ  Functions(total) Functions(LC) Functions(SC)'
C      
      IOFT    = 1
      IOFU    = 1
      IOFVS   = 1
      IOFVL   = 1
      INVMAX  = 3-NFSYM  ! INV=2 for NFSYM=1 and INV=1 for NFSYM=2
      DO 40 IFRP = 1,NFSYM
C
C        output for ABS(Kappa)=1,
C        which is NOT included in the below loop (IKAP = 2,KPMAX(IFRP))
C
       IF(NFSYM.EQ.2)THEN
         WRITE (LUPRI,1020) 3-2*IFRP,(-1)**IFRP,1,
     &          (NORB_SUB(1,IFRP,IC),IC=0,2)
         KPADD = 2
         NSKIP = NORB_SUB(1,IFRP,0)
       ELSE
         WRITE (LUPRI,1020) 3-2*IFRP,-1,1,(NORB_SUB(1,IFRP,IC),IC=0,2)
         WRITE (LUPRI,1020) 3-2*IFRP, 1,1,(NORB_SUB(2,IFRP,IC),IC=0,2)
         KPADD = 3
         NSKIP = NORB_SUB(1,IFRP,0)+NORB_SUB(2,IFRP,0)
       ENDIF
       NORB1 = NFBAS(IFRP,0)*NSKIP
C     
       IOFVS = IOFVS + NORB1
       IOFVL = IOFVL + NORB1
       IOFT  = IOFT  + NORB1
       IF(NORB(IFRP).EQ.0) GOTO 40
C
C........allocate KFLL,KVLL,KFSS,KVSS once and for all
C
       NLL=0
       NSS=0
       KPBUF = KPADD
       DO IKAP = 2,KPMAX(IFRP)
       DO INV=1,INVMAX          
         NLL   = MAX(NLL,NORB_SUB(KPBUF,IFRP,1)*NORB_SUB(KPBUF,IFRP,1))
         NSS   = MAX(NSS,NORB_SUB(KPBUF,IFRP,2)*NORB_SUB(KPBUF,IFRP,2))
         KPBUF = KPBUF + IKAP      
       ENDDO
       ENDDO
       CALL MEMGET('REAL',KFLL,NLL,WORK,KFREE,LFREE)
       CALL MEMGET('REAL',KFSS,NSS,WORK,KFREE,LFREE)
       CALL MEMGET('REAL',KVLL,NLL,WORK,KFREE,LFREE)
       CALL MEMGET('REAL',KVSS,NSS,WORK,KFREE,LFREE)
C
       IF(do4c2c)THEN 
C       make U_Mj only for kappa = 1,-1, because U_Mj is diagonal at Kappa = -1,1 
C       kappa = not 1,-1 will be made later       
         allocate(TMPB(NORB(IFRP)*NORB(IFRP)))       
         CALL DZERO(TMPB,NORB(IFRP)*NORB(IFRP))
         CALL SBUNIT(1,NSKIP,NORB(IFRP),TMPB)  
       ENDIF
C      
C........Here we transform the MJ-matrix for particular kappa; look TSUBBL
C       
       DO IKAP = 2,KPMAX(IFRP)
       KPOFF=IKAP*(IKAP-1)/2
       DO INV=1,INVMAX ! IF (NFSYM.EQ.1) we cycle twice to get both signs of kappa
         NL = NORB_SUB(KPADD,IFRP,1)
         NS = NORB_SUB(KPADD,IFRP,2)
         IEIGP = IORB(IFRP) + 1
         IEIGE = IEIGP + NS
         IF(NFSYM.EQ.2)THEN
           ISKP = (-1)**(IKAP+IFRP+1)
         ELSE
           ISKP = (-1)**INV
         ENDIF
C           
CCCCCCCCCCCCCCCCCC         
C Large AO-MO transformation (diagonal in kappa because of VMAT => VU)
CCCCCCCCCCCCCCCCCC         
         IF (NL.GT.0) THEN
           CALL TSBBL2(IFRP,1,KPADD,TMPA,WORK(KFLL),VMAT,VMAT,
     &                 WORK(KFREE),LFREE,IPRHAM) 
CCCCCCCCCCCCCCCCCCCCCCC 
C Large diagonalization
CCCCCCCCCCCCCCCCCCCCCCC         
C
           IF (IPRHAM.GE.6) THEN
             CALL HEADER('ATMSYM: Mj_LL to be diagonalized',-1)
             WRITE(LUPRI,'(2X,A,L1)') 'DOJACO=',DOJACO
             WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &       '*** Fermion corep ',IFRP,'/',NFSYM
             CALL PRQMAT(WORK(KFLL),NL,NL,NL,NL,1,
     &                   IPQTOQ(1,0),LUPRI)
           ENDIF
C        
           IF (.NOT.DOQJACO) THEN
             IF(DOJACO) THEN
               CALL RSJACO(NL,NL,NL,WORK(KFLL),EIG(IEIGE),1,1,0,
     &                       WORK(KVLL))
             ELSE
               CALL QDIAG(1,NL,WORK(KFLL),NL,NL,
     &                    EIG(IEIGE),1,
     &                    WORK(KVLL),NL,NL,
     &                    WORK(KFREE),LFREE,IERR)
               IF (IERR.NE.0) THEN
                 WRITE(LUPRI,'(2X,A)')
     &           'ATMSYM: Erroneous QDIAG diagonalization of Mj_LL!'
                 CALL QUIT(
     &          'ATMSYM: Erroneous QDIAG diagonalization of Mj_LL !')
               ENDIF
             ENDIF
           ELSE
             CALL QJACOBI(WORK(KFLL),WORK(KVLL),NL,1,0,
     &                    IDUMMY,.TRUE.,IPRHAM)
             DO II=1,NL
               EIG(IEIGE+II-1)=WORK(KFLL+NL*(II-1)+II-1) ! store eigenvalues
             ENDDO
           ENDIF
           MJMAXL = NINT( MAX( ABS(EIG(IEIGE)), ABS(EIG(IEIGE+NL-1)) ) )
           IF (IPRHAM.GE.3) THEN
             CALL HEADER('ATMSYM: Mj_LL eigenvalues',-1)
             WRITE(LUPRI,'(A,I2)') ' * Kappa:', ID_SUB_BL(KPADD,IFRP)
             WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &       '*** Fermion corep ',IFRP,'/',NFSYM
             DO II = 1,NL
               WRITE(LUPRI,*) II,EIG(IEIGE+II-1)
             ENDDO
             CALL FLSHFO(LUPRI)
           ENDIF
           IF (IPRHAM.GE.5) THEN
             CALL HEADER('ATMSYM: Mj_LL eigenvectors',-1)
             WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &       '*** Fermion corep ',IFRP,'/',NFSYM
             CALL PRQMAT(WORK(KVLL),NL,NL,NL,NL,1,
     &            IPQTOQ(1,0),LUPRI)
             CALL FLSHFO(LUPRI)
           ENDIF
         ENDIF !NL.GT.0
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC         
C Small AO-MO transformation (diagonal in kappa because of VMAT => VU)
CCCCCCCCCCCCCCCCCCCCCCCCCCCC         
         IF (NS.GT.0) THEN
           CALL TSBBL2(IFRP,2,KPADD,TMPA,WORK(KFSS),VMAT,VMAT,
     &                 WORK(KFREE),LFREE,IPRHAM) 
CCCCCCCCCCCCCCCCCCCCCCC         
C Small diagonalization
CCCCCCCCCCCCCCCCCCCCCCC         
C
           IF (IPRHAM.GE.6) THEN
             CALL HEADER('ATMSYM: Mj_SS to be diagonalized',-1)
             WRITE(LUPRI,'(2X,A,L1)') 'DOJACO=',DOJACO
             WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &       '*** Fermion corep ',IFRP,'/',NFSYM
             CALL PRQMAT(WORK(KFSS),NS,NS,NS,NS,1,
     &                IPQTOQ(1,0),LUPRI)
           ENDIF
C     
           IF (.NOT.DOQJACO) THEN
             IF(DOJACO) THEN
               CALL RSJACO(NS,NS,NS,WORK(KFSS),EIG(IEIGP),1,1,0,
     &                     WORK(KVSS))
             ELSE
               CALL QDIAG(1,NS,WORK(KFSS),NS,NS,
     &                    EIG(IEIGP),1,
     &                    WORK(KVSS),NS,NS,
     &                    WORK(KFREE),LFREE,IERR)
               IF (IERR.NE.0) THEN
                 WRITE(LUPRI,'(2X,A)')
     &           'ATMSYM: Erroneous QDIAG diagonalization of Mj_SS!'
                 CALL QUIT(
     &           'ATMSYM: Erroneous QDIAG diagonalization of Mj_SS !')
               ENDIF
             ENDIF
           ELSE
             CALL QJACOBI(WORK(KFSS),WORK(KVSS),NS,1,0,IDUMMY,
     &                   .TRUE.,IPRHAM)
             DO II=1,NS
               EIG(IEIGP+II-1)=  WORK(KFSS+NS*(II-1)+II-1) ! store eigenvalues
             ENDDO
           ENDIF
           MJMAXS=NINT(MAX(ABS(EIG(IEIGP)),ABS(EIG(IEIGP+NS-1))))
           IF (IPRHAM.GE.3) THEN
             CALL HEADER('ATMSYM: Mj_SS eigenvalues',-1)
             WRITE(LUPRI,'(A,I2)') ' * Kappa:', ID_SUB_BL(KPADD,IFRP)
             WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &         '*** Fermion corep ',IFRP,'/',NFSYM
             DO II = 1,NS
               WRITE(LUPRI,*) II,EIG(IEIGP+II-1)
             ENDDO
             CALL FLSHFO(LUPRI)
           ENDIF
           IF (IPRHAM.GE.5) THEN
             CALL HEADER('ATMSYM: Mj_SS eigenvectors',-1)
             WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &        '*** Fermion corep ',IFRP,'/',NFSYM
             CALL PRQMAT(WORK(KVSS),NS,NS,NS,NS,1,
     &                   IPQTOQ(1,0),LUPRI)
             CALL FLSHFO(LUPRI)
           ENDIF
C      
           MJBUF = 0
           DO IX = IEIGP, IEIGP+NS-1
             MJ = NINT(D2*EIG(IX))
             IF(MJ.NE.MJBUF) THEN
               MJBUF = MJ
               MJABS = (ABS(MJ)+1)/2
               ID_SUB_BL(KPADD+MJABS-1,IFRP)=ISKP*(KPOFF+MJABS)
             ENDIF
           ENDDO        
C
           DO IX = IEIGP, IEIGP+NS-1
             EIG(IX) = ABS(EIG(IX))
           ENDDO
C
C          The SELBOS routine does just the reordering that we want,
C          use work(kfss) as scratch, and allocate a few other small
C          scratch arrays. The IBEIG part will not be used but we will
C          initialize it to prevent problems with uninitialized data while
C          sorting it.
C
           CALL MEMGET('INTE',KINDX,NS,WORK,KFREE,LFREE)
           CALL MEMGET('INTE',KIDUM,NS,WORK,KFREE,LFREE)
           CALL MEMGET('INTE',KIBEI,NS,WORK,KFREE,LFREE)
           MJ = 0
           CALL ICOPY (NS,MJ,0,WORK(KIBEI),1)
           CALL SELBOS (IFRP,NS,WORK(KINDX),EIG(IEIGP),
     &                  WORK(KIBEI),WORK(KVSS),
     &                  WORK(KFSS),WORK(KIDUM),IPRHAM)
           CALL MEMREL('ATMSYM',WORK,KINDX,KINDX,KFREE,LFREE)
         ENDIF                     !NS.GT.0
C
C        Same thing for large component solutions
C
         MJBUF = 0
         DO IX = IEIGE, IEIGE+NL-1
           MJ = NINT(D2*EIG(IX))
           IF(MJ.NE.MJBUF) THEN
             MJBUF = MJ
             MJABS = (ABS(MJ)+1)/2
             ID_SUB_BL(KPADD+MJABS-1,IFRP)=ISKP*(KPOFF+MJABS)             
           ENDIF
         ENDDO        
C
         DO IX = IEIGE, IEIGE+NL-1
           EIG(IX) = ABS(EIG(IX))
         ENDDO
         CALL MEMGET('INTE',KINDX,NL,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KIDUM,NL,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KIBEI,NL,WORK,KFREE,LFREE)
         MJ = 0
         CALL ICOPY (NL,MJ,0,WORK(KIBEI),1)
         CALL SELBOS (IFRP,NL,WORK(KINDX),EIG(IEIGE),
     &                WORK(KIBEI),WORK(KVLL),
     &                WORK(KFLL),WORK(KIDUM),IPRHAM)
         CALL MEMREL('ATMSYM',WORK,KINDX,KINDX,KFREE,LFREE)
C
         DO IC = 0,MC !MC = 1 for L only and 2 for L+S
           NDIM = NORB_SUB(KPADD,IFRP,IC)/IKAP 
           DO IK = 1,IKAP
             NORB_SUB(KPADD+IK-1,IFRP,IC) = NDIM
           ENDDO
         ENDDO
C
C     Write information to output
         DO IK = 1,IKAP
           ID = ID_SUB_BL(KPADD+IK-1,IFRP)
           IF(NORB_SUB(KPADD+IK-1,IFRP,0).NE.0) THEN
             CALL ATOMIC_ID(ID,KP,J,MJ,LL)
             WRITE (LUPRI,1020) 3-2*IFRP,KP,MJ,
     &             (NORB_SUB(KPADD+IK-1,IFRP,IC),IC=0,MC)
           ENDIF
         ENDDO
C
         IOFS = 0
         IOFL = 0
         IOFVL = IOFVL + NFBAS(IFRP,0)*NS


         DO IK = 1,IKAP
           IOFK = KPADD+IK-1
           IF(NS.GT.0)THEN
             CALL BCKTRA(TMAT(IOFT),NFBAS(IFRP,0),
     &                   NORB_SUB(IOFK,IFRP,2),
     &                   WORK(KVSS+IOFS),NS,NORB_SUB(IOFK,IFRP,2),
     &                   NS,1,
     &                   NORB_SUB(IOFK,IFRP,2),1,NFBAS(IFRP,0),
     &                   VMAT(IOFVS),NFBAS(IFRP,0),NS,1,
     &                   IPRHAM)
             IF(do4c2c)THEN
C... Insert WORK(KFSS) to TMPB (U_Mj) for each Mj
                CALL INSRTBL(IFRP,2,KPADD,IOFK,
     &                     NTMO_SUB,NORB_SUB,NORB_SUB(IOFK,IFRP,2),
     &                     NTMO_SUB(KPADD,IFRP,2),NORB(IFRP),
     &                     WORK(KFSS+IOFS),TMPB(1),IPRHAM)
             ENDIF
           ENDIF           
C
C...  Update small here: : WORK(KVSS+IOFS),VMAT(IOFV),TMAT(IOFV)
C
           IOFS = IOFS + NS*NORB_SUB(IOFK,IFRP,2)
           IOFT = IOFT + NFBAS(IFRP,0)*NORB_SUB(IOFK,IFRP,2)
           CALL BCKTRA(TMAT(IOFT),NFBAS(IFRP,0),
     &                 NORB_SUB(IOFK,IFRP,1),
     &                 WORK(KVLL+IOFL),NL,NORB_SUB(IOFK,IFRP,1),
     &                 NL,1,
     &                 NORB_SUB(IOFK,IFRP,1),1,NFBAS(IFRP,0),
     &                 VMAT(IOFVL),NFBAS(IFRP,0),NL,1,
     &                 IPRHAM)
           IF(do4c2c)THEN
C... Insert WORK(KFLL) to TMPB (U_Mj) for each Mj
             CALL INSRTBL(IFRP,1,KPADD,IOFK,
     &                 NTMO_SUB,NORB_SUB,NORB_SUB(IOFK,IFRP,1),
     &                 NTMO_SUB(KPADD,IFRP,1),NORB(IFRP),
     &                 WORK(KFLL+IOFL),TMPB(1),IPRHAM)
           ENDIF
C...  Update large here: : WORK(KVLL+IOFL),VMAT(IOFV),TMAT(IOFV)
           IOFL = IOFL + NL*NORB_SUB(IOFK,IFRP,1)
           IOFT = IOFT + NFBAS(IFRP,0)*NORB_SUB(IOFK,IFRP,1)
         ENDDO  ! IK = 1,IKAP
C     
         IOFVS = IOFVS + NFBAS(IFRP,0)*(NL+NS)
         IF(NS.NE.0)THEN 
           IOFVL = IOFVL + NFBAS(IFRP,0)*NS
         ELSE
           IOFVL = IOFVL + NFBAS(IFRP,0)*NORB_SUB(KPADD,IFRP,1)*IKAP
         ENDIF
         KPADD = KPADD + IKAP
       ENDDO                     ! INV=1,INVMAX
       ENDDO                     ! IKAP = 2,KPMAX(IFRP)
C
       CALL MEMREL('ATMSYM',WORK,KFLL,KFLL,KFREE,LFREE)
       IF(IFRP.EQ.1)IBUF2=IOFT 
       if(do4c2c)then            ! x2c-mmf approach
C... Making U (= U_kappa * U_Mj); U=TMPA, U_kappa=TBUF, U_Mj=TMPB           
        CALL DGEMM('N','N',NORB(IFRP),NORB(IFRP),NORB(IFRP),D1,
     &             TBUF(IOFU),NORB(IFRP),TMPB(1),NORB(IFRP),
     &             D0,TMPA(1),NORB(IFRP))
C... Since U should be in TBUF, copy TMPA to TBUF
        CALL DCOPY(NORB(IFRP)*NORB(IFRP),TMPA,1,TBUF(IOFU),1)
C... Extract only LL block (TMPA) from total U (TBUF)
        CALL EXTRSB (1,IFRP,NBORB,4,NBRP,
     &               NORB_SUB,MAX_SUB_BL,N_SUB_BL(IFRP),0,0,1,1,1,1,
     &               TMPA(1),NESH(IFRP),NESH(IFRP),1,
     &               TBUF(IOFU),NORB(IFRP),NORB(IFRP),1)
        call writt(13,NESH(IFRP)*NESH(IFRP),TMPA(1))
C
        IF(IPRHAM.GE.5)THEN        
         CALL HEADER(
     &    'ATMSYM: Unitary transformation'//
     &    ' matrix Kappa->Mj MO, TMPB',-1)         
         WRITE(LUPRI,'(A,I1,A,I1)')
     &   '* Fermion ircop no.',IFRP,'/',NFSYM
         CALL PRQMAT(TMPB(1),NORB(IFRP),NORB(IFRP),
     &        NORB(IFRP),NORB(IFRP),1,IPQTOQ(1,0),LUPRI)
C
         CALL HEADER(
     &    'ATMSYM: Final Unitary transformation'//
     &    ' matrix MO->MO (only LL), TMPA',-1)         
         WRITE(LUPRI,'(A,I1,A,I1)')
     &   '* Fermion ircop no.',IFRP,'/',NFSYM
         CALL PRQMAT(TMPA(1),NESH(IFRP),NESH(IFRP),
     &        NESH(IFRP),NESH(IFRP),1,IPQTOQ(1,0),LUPRI)
       ENDIF 
        deallocate(TMPB)
       end if !do4c2c
       IOFU = IOFU + NORB(IFRP)*NORB(IFRP)
C
C
 40   CONTINUE !IFRP
C.....Also initialize NTMO_SUB which will be the number of orbitals
C      in each subblock in orthonormal basis      
      CALL ICOPY(MAX_SUB_BL*3*2,NORB_SUB,1,NTMO_SUB,1)
      deallocate(TMPA)
      if(do4c2c)then ! x2c-mmf approach
        close(13,status='keep')
      end if
C     
C
C... Debug sentences
C      
      IOFT = 1
      IOFU = 1
      IOFV = 1
      IMAT = 1
      IBUF = 1 
      IF (IPRHAM.GE.5) THEN
            DO 50 IFRP = 1,NFSYM
            IF(NORB(IFRP).EQ.0) GOTO 50      
      
C
            CALL HEADER(
     &'ATMSYM: Final unitary transformation matrix'// 
     &' in orthon. basis U (MO->MO),TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',IFRP,'/',NFSYM
         CALL PRQMAT(TBUF(IOFU),NORB(IFRP),NORB(IFRP),NORB(IFRP),
     &                 NORB(IFRP),1,IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
C
         CALL HEADER(
     &'ATMSYM: TMAT (symm.blocked) transformation'//
     &' matrix (AO->MO) aka VU TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',IFRP,'/',NFSYM
         CALL PRQMAT(TMAT(IOFT),NFBAS(IFRP,0),NORB(IFRP),
     &         NFBAS(IFRP,0),NORB(IFRP),1,IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
C
         CALL HEADER(
     &    'ATMSYM: VMAT final transformation'//
     &    ' matrix 4c AO->Kappa MO aka V',-1)
         WRITE(LUPRI,'(A,I1,A,I1)')
     &   '* Fermion ircop no.',IFRP,'/',NFSYM
         CALL PRQMAT(VMAT(IMAT),NFBAS(IFRP,0),NORB(IFRP),
     &               NFBAS(IFRP,0),NORB(IFRP),1,IPQTOQ(1,0),LUPRI)
         WRITE(LUPRI,'(/5X,A/)')         
C
            IOFT = IOFT + NFBAS(IFRP,0)*NORB(IFRP)
            IOFU = IOFU + NORB(IFRP)*NORB(IFRP)
            IOFV = IOFV + NORB(IFRP)*NFBAS(IFRP,0)
            IMAT = IMAT + NFBAS(IFRP,0)*NORB(IFRP)*1
            IBUF = IBUF + NORB(IFRP)*NORB(IFRP)*1
C
 50         CONTINUE
      ENDIF            
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('ATMSYM')
      RETURN
 1000 FORMAT (/' Parity  Kappa  Functions(total) Functions(LC)'//
     &       ' Functions(SC)')
 1010 FORMAT (T5,I2,T9,I3,T20,I5,T33,I5,T46,I5)
 1020 FORMAT (T5,I2,T9,I3,T13,I4,'/2',T25,I5,T39,I5,T53,I5)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck amjmat */
      SUBROUTINE AMJMAT(FOCK,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Generate one-electron matrix of Jz
C     IPRINT = IPRHAM
C
C     Called from LINSYM
C
C     Written by Luuk Visscher, oct 1999
C
C*****************************************************************************
#include "implicit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION FOCK(*),WORK(LWORK)
C
      CALL QENTER('AMJMAT')
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET('REAL',KFTRI ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN,   NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST,NZ        ,WORK,KFREE,LFREE)
C
C     *************************************************
C     ****** Construct one-electron matrix of Jz ******
C     *************************************************
C
      CALL AMJMA1(FOCK,WORK(KFTRI),WORK(KONEIN),WORK(KFIRST),
     &            WORK,KFREE,LFREE,IPRINT)
C
C     Memory deallocation
      CALL MEMREL('AMJMAT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('AMJMAT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck amjma1 */
      SUBROUTINE AMJMA1(FOCK,FTRI,OP1INT,FIRST,WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C
C   Generate one-electron matrix of Jz
C
C  Written by Luuk Visscher, oct 1999
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbham.h"
C
      PARAMETER (HALF=0.5D0, D1=1.0D0)
      LOGICAL FIRST
      DIMENSION FOCK(*),FTRI(NNBBASX),OP1INT(NNBBASX),
     &          WORK(*)
      CALL QENTER('AMJMA1')
C
C     Sanity check : linear symmetry requires use of real groups (C2v or D2h)
C
      IF (NZ.NE.1) CALL QUIT ("Linear symmetry requires NZ=1")
C
C     Get integrals of component of the z-angular momentum that we want to
C     quantize on.
C
      FIRST = .TRUE.
      IOP =  IPANGMOM
      CALL PRPMAO(LU1INT,IOP,.TRUE.,WORK,.TRUE.,NNBBASX,FTRI,
     &         OP1INT,FIRST,IPRHAM)
      CALL DSPTSI(NTBAS(0),FTRI,FOCK)
C
C     Get integrals of the spin-z (sigma_z) operator
C     ===============================================
C
      FIRST = .TRUE.
      IOP =  IPSPNMOM
      CALL PRPMAO(LU1INT,IOP,.TRUE.,WORK,.TRUE.,NNBBASX,FTRI,
     &        OP1INT,FIRST,IPRHAM)
      CALL MEMGET('REAL',KBUF,N2BBASX,WORK,KFREE,LFREE)
      CALL DSPTSI(NTBAS(0),FTRI,WORK(KBUF))
C
C     Make the Jz operator (in case of spinfree we make Lz)
C     =====================================================
C
      IF(.NOT.SPINFR) CALL DAXPY (N2BBASX,D1,WORK(KBUF),1,FOCK,1)
C
C     Release memory
C
      CALL MEMREL('AMJMA1',WORK,KBUF,KBUF,KFREE,LFREE)
C
C     Reindex to sorted basis
C     =======================
C
      CALL BUTOBS(FOCK,1,WORK(KFREE),LFREE)
C
C     Output section
C     ==============
C
      IF(IPRHAM.GE.4) THEN
        CALL HEADER(
     &  'AMJMA1: One-electron MJ matrix in symm.blocked SA-AO basis',-1)
        DO 10 I = 1,NFSYM
          IF(NFBAS(I,0).EQ.0) GOTO 10
          WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(FOCK(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
   10   CONTINUE
      ENDIF
C
      CALL QEXIT('AMJMA1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck amkmat */
      SUBROUTINE AMKMAT(FOCK,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Generate one-electron matrix of Kappa
C     IPRINT = IPRHAM
C
C     Routine based on AMJMAT
C
C     Called from ATMSYM
C
C     Written by A. Sunaga, July 2018
C
C*****************************************************************************
#include "implicit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION FOCK(*),WORK(LWORK)
C
      CALL QENTER('AMKMAT')
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET('REAL',KFTRI ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN,   NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST,NZ        ,WORK,KFREE,LFREE)
C
C     *************************************************
C     ****** Construct one-electron matrix of Kappa ***
C     *************************************************
C
      CALL AMKMA1(FOCK,WORK(KFTRI),WORK(KONEIN),WORK(KFIRST),
     &            WORK,KFREE,LFREE,IPRINT)
C
C     Memory deallocation
      CALL MEMREL('AMKMAT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('AMKMAT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck amkma1 */
      SUBROUTINE AMKMA1(FOCK,FTRI,OP1INT,FIRST,WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C
C  Generate one-electron matrix of Kappa
C
C  Written by A. Sunaga, July 2018
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbham.h"
C
      PARAMETER (HALF=0.5D0, D1=1.0D0, DM1=-1.0D0)
      LOGICAL FIRST
      DIMENSION FOCK(*),FTRI(NNBBASX),OP1INT(NNBBASX),
     &          WORK(*)
      CALL QENTER('AMKMA1')
C
C     Sanity check : ATOMIC symmetry requires use of real groups (C2v or D2h)
C
      IF (NZ.NE.1) CALL QUIT ("Atomic symmetry requires NZ=1")
C
C     Get integrals of component of the Spin-orbit operatar
C
      FIRST = .TRUE.
      IOP =  IPSPNORB
      CALL PRPMAO(LU1INT,IOP,.TRUE.,WORK,.TRUE.,NNBBASX,FTRI,
     &        OP1INT,FIRST,IPRHAM)
      CALL DSPTSI(NTBAS(0),FTRI,FOCK)
C
C     Get integrals of the Overlap integral
C     ===============================================
C
      FIRST = .TRUE.
      IOP =  IPOVRLAP
      CALL PRPMAO(LU1INT,IOP,.TRUE.,WORK,.TRUE.,NNBBASX,FTRI,
     &         OP1INT,FIRST,IPRHAM)
      CALL MEMGET('REAL',KBUF,N2BBASX,WORK,KFREE,LFREE)
      CALL DSPTSI(NTBAS(0),FTRI,WORK(KBUF))
C
C     Make the Kappa operator (in case of spinfree ..)
C     =====================================================
C
      IF(.NOT.SPINFR) CALL DAXPY (N2BBASX,DM1,WORK(KBUF),1,FOCK,1)
C
C     Release memory
C
      CALL MEMREL('AMKMA1',WORK,KBUF,KBUF,KFREE,LFREE)
C
C     Reindex to sorted basis
C     =======================
C
      CALL BUTOBS(FOCK,1,WORK(KFREE),LFREE)
C
C     Output section
C     ==============
C
      IF(IPRHAM.GE.4) THEN
        CALL HEADER(
     &  'AMKMA1: One-electron K matrix in symm.blocked SA-AO basis',-1)
        DO 10 I = 1,NFSYM
          IF(NFBAS(I,0).EQ.0) GOTO 10
          WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(FOCK(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
   10   CONTINUE
      ENDIF
C
      CALL QEXIT('AMKMA1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck sbunit */
      SUBROUTINE SBUNIT(NINI,NSUB,NTOT,FFUL)
C***********************************************************************
C 
C     Making a sub unit matrix 
C
C     A Sunaga 2020
C
C***********************************************************************
#include "implicit.h"
      DIMENSION FFUL(NTOT,NTOT)
C
      DO I = NINI,NSUB
         FFUL(I,I) = 1.0D0
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck insrtbl */
      SUBROUTINE INSRTBL(IFRP,IC,KPSUB,MJSUB,KPBLCK,MJBLCK,
     &                   NINSRT,NSUB,NFUL,
     &                   FSUB,FFUL,IPRNT)
C***********************************************************************
C 
C     Insert a subblock to a full matrix
C     Moving pointers of row and column of FFUl separately       
C
C     IFRP.... Fermion irreps     
C     IC...... 1 - Large; 2 - Small
C      
C     FSUB..... subblock
C               Eigenvector for specified (Kappa,Mj)
C     NSUB..... Number of Kappa-MO basis
C     NINSRT... Number of eigenvectors (Kappa,Mj)
C      
C     FFUL... full matrix
C     NFUL... size of FFUL      
C      
C     A Sunaga 2021 Jun
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
C
      DIMENSION KPBLCK(MAX_SUB_BL,2,0:2),MJBLCK(MAX_SUB_BL,2,0:2),
     &          FSUB(NSUB,NINSRT),FFUL(NFUL,NFUL)
C     
C... pointer for Kappa
C      
      IF(KPSUB.GE.1)THEN
        IKP = 0 
        DO IS = 1, KPSUB-1
           IKP = IKP + KPBLCK(IS,IFRP,0)
        ENDDO   
      ELSE
        IKP = 0
      ENDIF
C      
      IF(IC.EQ.1)THEN
        KPBLCK_BUF = NSUB
      ELSEIF(IC.EQ.2)THEN
        KPBLCK_BUF = 0
      ENDIF
      ITMKP = IKP+KPBLCK_BUF !final pointer
C
C... pointer for Mj
C      
      IF(MJSUB.GE.1)THEN
        IMJ = 0 
        DO IS = 1, MJSUB-1
           IMJ = IMJ + MJBLCK(IS,IFRP,0)
        ENDDO   
      ELSE
        IMJ = 0
      ENDIF
C
      IF(MJBLCK(MJSUB,IFRP,IC).GT.0) THEN
        IF(IC.EQ.1)THEN
          MJBLCK_BUF = NINSRT
        ELSEIF(IC.EQ.2)THEN
          MJBLCK_BUF = 0
        ENDIF
        ITMMJ = IMJ + MJBLCK_BUF !final pointer
C
        DO I = 1,NINSRT
           CALL DCOPY(NSUB,FSUB(1,I),1,FFUL(ITMKP+1,ITMMJ+I),1)
        ENDDO
      ENDIF 
C
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Wrifrz */
      SUBROUTINE WRIFRZ(IUNIT,IFRP,CSEL,ESEL,ISEL,IVEC,BUF,IBUF,
     &                  NVEC,NDEL)
C***********************************************************************
C
C     Dump frozen orbitals to file.
C     Also adjust orbital parameters accordingly      
C     Written by T. Saue May 2003
C
C***********************************************************************
#include "implicit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dcbham.h"
      DIMENSION CSEL(*),ESEL(*),ISEL(*),
     &          BUF(*),IVEC(*),IBUF(*),NDEL(2)
C
C     Get index array for frozen orbitals;
C
      NTOT=1
      NTOP = NESH(IFRP) + NDEL(1)
      NLOW = NPSH(IFRP) + NDEL(2)
      CALL NUMLS3(VCFROZ(IFRP),IVEC,NFBAS(IFRP,0),
     &            -NLOW,NTOP,NTOT,NZERO)
C     
C     Consistency check
C
      IF(NFRO(IFRP).NE.(NTOT-NZERO)) THEN
        CALL QUIT('WRIFRZ:Inconsistency !')
      ENDIF
C
C     Selected orbitals that are not frozen have index zero;
C     they must be taken out and pointers made to the remaining ones
C     Write compressed index array to file
C
C     We recompute NASHMFT as this value may change if there are frozen open shells
C
      NE = 0
      NP = 0
      NASHMFT = 0
      II = 0
      DO I = 1,NTOT
      IP = IVEC(I)
      IF(IP.NE.0) THEN
        II = II + 1
        IBUF(II)=I
C       Positronic position
        IF(IP.LT.0.AND.IP.GE.(-NLOW)) THEN
          IVEC(II)=1+IP
          NP = NP + 1
          IF (SUB_BL) THEN
            DO ISUB = 1,N_SUB_BL(IFRP)
            IF(ISEL(I).EQ.ID_SUB_BL(ISUB,IFRP)) THEN
              NORB_SUB(ISUB,IFRP,2)=NORB_SUB(ISUB,IFRP,2)+1
              NORB_SUB(ISUB,IFRP,0)=NORB_SUB(ISUB,IFRP,0)+1
            ENDIF
            ENDDO
          ENDIF
C       Electronic position
        ELSEIF (IP.GT.0.AND.IP.LE.NTOP) THEN
          IVEC(II)=IP
          NE = NE + 1
          IF (SUB_BL) THEN
            DO ISUB = 1,N_SUB_BL(IFRP)
            IF(ISEL(I).EQ.ID_SUB_BL(ISUB,IFRP)) THEN             
              NORB_SUB(ISUB,IFRP,1)=NORB_SUB(ISUB,IFRP,1)+1
              NORB_SUB(ISUB,IFRP,0)=NORB_SUB(ISUB,IFRP,0)+1
            ENDIF
            ENDDO
          ENDIF
C         Check if orbital belongs to closed or open shells.
C         If so, subtract from occupation number
          IF(IP.GT.0.AND.IP.LE.NISH(IFRP)) THEN
            NISHMF(IFRP) = NISHMF(IFRP) - 1
          ELSE
            I1 = NISH(IFRP)
            DO IOPEN = 1,NOPEN
              I2 = I1 + NACSH(IFRP,IOPEN)
              IF(IP.GT.I1.AND.IP.LE.I2) THEN
                NACSHMF(IFRP,IOPEN) = NACSHMF(IFRP,IOPEN) - 1
                GOTO 20
              ENDIF
              I1 = I2
            ENDDO
 20         CONTINUE
          ENDIF
        ENDIF
      ENDIF
      ENDDO
      NOCCMF(IFRP) = NISHMF(IFRP)
      DO IOPEN = 1,NOPEN
        NOCCMF(IFRP) = NOCCMF(IFRP) + NACSHMF(IFRP,IOPEN)
        NASHMFT = NASHMFT + NACSHMF(IFRP,IOPEN)
      ENDDO
C
C     Final adjustments
C
      NPSH(IFRP) = NPSH(IFRP) + NP
      NESH(IFRP) = NESH(IFRP) + NE
      DO I = 1,NFRO(IFRP)
        IVEC(I) = IVEC(I) + NPSH(IFRP)
      ENDDO
C
C     Write index array to file
C
      CALL WRITI(IUNIT,NFRO(IFRP),IVEC)
C
C     Write coefficients and eigenvalues to file
C
      IDIM = NFBAS(IFRP,0)*NVEC
      JDIM = NFBAS(IFRP,0)*NZ
      DO II = 1,NFRO(IFRP)
        I = IBUF(II)
        IOFF = NFBAS(IFRP,0)*(I-1)+1
        JOFF = 1
        DO IZ = 1,NZ
          CALL DCOPY(NFBAS(IFRP,0),CSEL(IOFF),1,BUF(JOFF),1)
          IOFF = IOFF + IDIM
          JOFF = JOFF + NFBAS(IFRP,0)
        ENDDO
        CALL WRITT(IUNIT,JDIM,BUF)
 99     WRITE(IUNIT) ESEL(I),ISEL(I)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Inisub */
      SUBROUTINE INISUB
C***********************************************************************
C
C     Copy information from NBORB to the subblock array if we have not
C     partioned the orbital space otherwise.
C
C     Written by Luuk and copied by Trond Oct 19 2004
C
C***********************************************************************
#include "implicit.h"
C
#include "dcborb.h"
#include "dgroup.h"
      DO IFRP = 1, NFSYM
        N_SUB_BL(IFRP) = 4 / NZ
        DO IB = 1, N_SUB_BL(IFRP)
          NORB_SUB(IB,IFRP,0) = NBORB(IB,IFRP,0)
          NORB_SUB(IB,IFRP,1) = NBORB(IB,IFRP,1)
          NORB_SUB(IB,IFRP,2) = NBORB(IB,IFRP,2)
          IREP = JFSYM(IB,IFRP) - 1
          ID_SUB_BL(IB,IFRP) = IREP
        ENDDO
      ENDDO
C.....Also initialize NTMO_SUB which will be the number of orbitals in each subblock in orthonormal basis
      CALL ICOPY(MAX_SUB_BL*3*2,NORB_SUB,1,NTMO_SUB,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck DELPOS */
      SUBROUTINE DELPOS(TMAT,TMOD,FMAT,EIG,IBEIG,VMAT,WORK,KFREE,LFREE)
C***********************************************************************
C
C
C     Deleting negative-energy solutions of:
C     ======================================
C       FREEPJ - free particle matrix
C       VEXTPJ - bare nucleus matrix
C
C     On input:
C       TMAT(NFORB,NFORB,IFRP)
C            - free particle/bare nucleus matrix
C              in orthonormal canonical Lowdin basis
C       TMOD(NFORB,NORB,IFRP)
C            - transformation matrix with modified Dirac equation
C              embedded, in orthonormal basis
C       FMAT - empty
C     On output:
C       TMAT - transformation matrix with negative-energy solutions deleted
C     Written by T.Saue, October 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "dgroup.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
      DIMENSION TMAT(*),TMOD(*),FMAT(*),EIG(*),IBEIG(*),VMAT(*),WORK(*)
C
      NSYM = 4/NZ
C
C     Transform free particle/bare nucleus matrix to modified Dirac form;
C
      IOFT = 1
      IOFV = 1
      DO 10 I = 1,NFSYM
        IF(NTMO(I).EQ.0) GOTO 10
        CALL QTRANS('AOMO','S',D0,NFORB(I,0),NFORB(I,0),
     &       NTMO(I),NTMO(I),
     &       TMAT(IOFT),NFORB(I,0),NFORB(I,0),NZ,IPQTOQ(1,0),
     &       FMAT(1+I2TMOT(I)),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &       TMOD(1+I2TMT(I)),NFORB(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &       TMOD(1+I2TMT(I)),NFORB(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)
        IF (IPRHAM.GE.6) THEN
          IF(VEXTPJ) THEN
            CALL HEADER('DELPOS: '//
     &      'Bare nucleus Dirac matrix in mod. Dirac eq. MO-basis',-1)
          ELSE
            CALL HEADER('DELPOS: '//
     &      'Free particle Dirac matrix in mod. Dirac eq. MO-basis',-1)
          ENDIF
          CALL PRQMAT(FMAT(1+I2TMOT(I)),NTMO(I),NTMO(I),
     &              NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
C       get full transformation as well !
        CALL BCKTR1(TMAT(1+I2ORBT(I)),NFBAS(I,0),NTMO(I),
     &              TMOD(1+I2TMT(I)),NFORB(I,0),NTMO(I),
     &              NFORB(I,0),NZ,NTMO(I),1,NFBAS(I,0),
     &              VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,
     &              IPRHAM)
        IF (IPRHAM.GE.6) THEN
          CALL HEADER('DELPOS: Transformation matrix',-1)
          CALL PRQMAT(TMAT(1+I2ORBT(I)),NFBAS(I,0),NTMO(I),
     &                NFBAS(I,0),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
        IOFT = IOFT + NFORB(I,0)*NFORB(I,0)*NZ
        IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)
 10   CONTINUE
C
C     Diagonalize: eigenvectors are now stored in TMOD
C
      CALL DFDIAG (FMAT,EIG,IBEIG,TMOD,.FALSE.,WORK,KFREE,LFREE)
C
C     Backtransform only electronic solutions
C
      IBUF = 1
      DO 20 I = 1,NFSYM
        IF(NESH(I).EQ.0) GOTO 20
          IF(IPRHAM.GE.5) THEN
            CALL HEADER('DELPOS: Eigenvalues',-1)
            WRITE(LUPRI,'(I5,E16.8)') (J,EIG(IORB(I)+J),J=1,NTMO(I))
          ENDIF
          CALL BCKTR1(FMAT(IBUF),NFBAS(I,0),NESH(I),
     &                TMOD(1+I2TMOT(I)),NTMO(I),NTMO(I),NTMO(I),NZ,
     &                NESH(I),(NPSH(I)+1),NFBAS(I,0),
     &                TMAT(1+I2ORBT(I)),NFBAS(I,0),NTMO(I),NZT,IPRHAM)
          NPSH(I) = 0
          NORB(I) = NESH(I)
          NBUF = NFBAS(I,0)*NESH(I)*NZ
          CALL DCOPY(NBUF,FMAT(IBUF),1,TMAT(IBUF),1)
          IBUF = IBUF + NBUF
          NTMO(I) = NORB(I)
C         Finish by adjusting the number of positron solutions to zero
          DO ISYM = 1, NSYM
            NBORB(ISYM,I,2) = 0
            NBORB(ISYM,I,0) = NBORB(ISYM,I,1) + NBORB(ISYM,I,2)
          ENDDO
 20     CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck LL_EQ */
      SUBROUTINE LL_EQ(TMAT,TBUF,VMAT)
C***********************************************************************
C
C     The Levy-Leblond equation:
C     ==========================
C
C       Set up basis for Levy-Leblond equation as:
C
C              [                (Psi(L), 0 ]
C        Psi = [                           ]
C              [(2mc)^-1 sigma.p Psi(L), 0 ]
C
C     Reference:
C       L. Visscher and T. Saue, J. Chem. Phys. 113(2000) 3996
C       "Approximate relativistic electronic structure methods based on
C       the quaternion modified Dirac equation"
C
C     On input:
C       TMAT - free particle matrix in Lowdin orthonormal basis
C       VMAT - Lowdin canonical transformation matrix
C
C     On output:
C       TMAT - transformation matrix with Levy-Leblond equation embedded
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, D2 = 2.0D0, D0 = 0.0D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
C
      DIMENSION TMAT(*),TBUF(*),VMAT(*)
C
      FAC  = D2*CVAL*CVAL
      FAC  = D1/FAC
      NSYM = 4/NZ
C
      IBUF  = 1
      IMAT  = 1
      DO 10 IFRP = 1,NFSYM
        NESH(IFRP)  = NFORB(IFRP,1)
        NPSH(IFRP)  = 0
        NORB(IFRP)  = NESH(IFRP) + NPSH(IFRP)
        NTMO(IFRP)  = NORB(IFRP)
        IF(NORB(IFRP).EQ.0) GOTO 10
        CALL DZERO(TBUF(IBUF),NFORB(IFRP,0)*NORB(IFRP)*NZT)
C.......set LL block to unit matrix
        JJ = IBUF
        DO J = 1, NFORB(IFRP,1)
          TBUF(JJ) = D1
          JJ = JJ + NFORB(IFRP,0) + 1
        ENDDO
C.......set SL block to (2mc)^-1 sigma.p
        DO IZ = 1, NZ
          JSL = IBUF + NFORB(IFRP,1)
          ISL = IMAT + NFORB(IFRP,1)
          DO J = 1, NFORB(IFRP,1)
            CALL DCOPY (NFORB(IFRP,2),TMAT(ISL),1,TBUF(JSL),1)
            CALL DSCAL (NFORB(IFRP,2),FAC,TBUF(JSL),1)
            ISL = ISL + NFORB(IFRP,0)
            JSL = JSL + NFORB(IFRP,0)
          ENDDO
          IF (IPRHAM.GE.6) THEN
           CALL HEADER('Transformation matrix in orthon. basis',-1)
           CALL PRQMAT(TBUF(IBUF),NFORB(IFRP,0),NORB(IFRP),
     &                 NFORB(IFRP,0),NORB(IFRP),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
          IMAT = IMAT + NFORB(IFRP,0) * NFORB(IFRP,0)
          IBUF = IBUF + NFORB(IFRP,0) * NORB(IFRP)
        ENDDO
C.......finish by adjusting the number of positron solutions to zero
        DO ISYM = 1, NSYM
          NBORB(ISYM,IFRP,2) = 0
          NBORB(ISYM,IFRP,0) = NBORB(ISYM,IFRP,1) + NBORB(ISYM,IFRP,2)
        ENDDO
   10 CONTINUE
C
C    Generate transformation matrix
C    ==============================
C
      IBUF = 1
      IMAT = 1
      IOFV = 1
      DO 20 IFRP = 1,NFSYM
        IF(NORB(IFRP).EQ.0) GOTO 20
        CALL BCKTR1(TMAT(IMAT),NFBAS(IFRP,0),NFORB(IFRP,1),
     &              TBUF(IBUF),NFORB(IFRP,0),NFORB(IFRP,1),
     &              NFORB(IFRP,0),NZ,
     &              NFORB(IFRP,1),1,NFBAS(IFRP,0),
     &              VMAT(IOFV),NFBAS(IFRP,0),NFORB(IFRP,0),1,
     &              IPRGEN)
        IMAT = IMAT + NFBAS(IFRP,0)*NORB(IFRP)*NZ
        IBUF = IBUF + NFORB(IFRP,0)*NORB(IFRP)*NZ
        IOFV = IOFV + NFORB(IFRP,0)*NFBAS(IFRP,0)
   20   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MOD_DIREQ */
      SUBROUTINE MOD_DIREQ(TMAT,TBUF,VMAT,EIG,WORK,LWORK)
C***********************************************************************
C
C       The modified Dirac eqation:
C       ===========================
C
C       Basis as in Levy-Leblond but without fixed ratio
C
C
C               [                (Psi(L),                      0  ]
C         Psi = [                                                 ]
C               [                      0, (2mc)^-1 sigma.p Psi(L) ]
C
C     Reference:
C       L. Visscher and T. Saue, J. Chem. Phys. 113(2000) 3996
C       "Approximate relativistic electronic structure methods based on
C       the quaternion modified Dirac equation"
C
C     On input:
C       TMAT - free particle matrix in Lowdin orthonormal basis
C       VMAT - Lowdin canonical transformation matrix
C
C     On output:
C       TMAT - transformation matrix with modified Dirac equation embedded
C
C     Note that there is a difference in boson irrep ordering
C     between TMAT and the Lowdin canonical transformation matrix VMAT:
!
!
C     After transformation with VMAT the boson irrep ordering
C     in a given fermion ircop is:
C       L1 - L2 - L3 .. S1 - S2 - S3 - ..
C     whereas with TMAT the ordering is:
C       S1 - L1 - S2 - L2 - S3 - L3 ..
C
C     where L and S refer to the large and small components,
C     respectively.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, DM2 = -2.0D0, D0 = 0.0D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "cbihr1.h"
      DIMENSION TMAT(*),TBUF(*),VMAT(*),EIG(*),WORK(*)
#include "memint.h"

      NSYM = 4/NZ
      IMAT = 1
      IBUF = 1
      DO 10 I = 1,NFSYM
        NESH(I)  = NFORB(I,1)
        NPSH(I)  = NFORB(I,1)
        NORB(I)  = NESH(I) + NPSH(I)
        NTMO(I)  = NORB(I)
        IF(NORB(I).EQ.0) GOTO 10
C
C       Prepare small component part
C       ----------------------------
C
        LWRK = MXBBAS*MXBBAS*NSYM
        CALL MEMGET2('REAL','SMT',KSMT,LWRK,WORK,KFREE,LFREE)
        LWRK = NFORB(I,1) * NFORB(I,1)
        CALL MEMGET2('REAL','VMT',KVMT,LWRK,WORK,KFREE,LFREE)
        IVMT  = KVMT
        IBOL  = 0
        NEFF  = 0
        DO ISYM = 1,NSYM
        IF(NBORB(ISYM,I,1).GT.0) THEN
C
C         Make (sigma.p)(sigma.p) = p2
C
          IREP = JFSYM(ISYM,I) - 1
          JSL = IMAT + IBOL * NFORB(I,0) + NFORB(I,1)
          JLS = IMAT + NFORB(I,1) * NFORB(I,0) + IBOL
          CALL QGEMM(NBORB(ISYM,I,1),NBORB(ISYM,I,1),NFORB(I,2),D1,
     &          'N','N',IPQTOQ(1,0),TMAT(JLS),NFORB(I,0),NFORB(I,0),NZ,
     &          'N','N',IPQTOQ(1,0),TMAT(JSL),NFORB(I,0),NFORB(I,0),NZ,
     &          D0,IPQTOQ(1,0),WORK(KSMT),
     &          NBORB(ISYM,I,1),NBORB(ISYM,I,1),1)
          IF (IPRHAM.GE.7) THEN
            WRITE(LUPRI,'(A,A3)') ' * Boson irrep: ',REP(IREP)
            CALL HEADER('Kinetic energy matrix in orthon. basis',-1)
            CALL PRQMAT(WORK(KSMT),NBORB(ISYM,I,1),NBORB(ISYM,I,1),
     &                  NBORB(ISYM,I,1),NBORB(ISYM,I,1),
     &                  1,IPQTOQ(1,0),LUPRI)
          ENDIF
          CALL LOWGEN(WORK(KSMT),NBORB(ISYM,I,1),WORK(IVMT),NEFFI,
     &                NBORB(ISYM,I,1),STOL(2),IPRONE,EIG,
     &                WORK(KFREE),LFREE)
          IF (NEFFI.LT.NBORB(ISYM,I,1)) WRITE(LUPRI,'(A)')
     &            '*** WARNING *** Linear dependence in RKB basis'
          NEFF = NEFF + NEFFI
          NBORB(ISYM,I,2) = NEFFI
          IF (IPRHAM.GE.7) THEN
            CALL HEADER('Ext.transf. matrix in orthon. basis',-1)
            CALL PRQMAT(WORK(IVMT),NBORB(ISYM,I,1),NBORB(ISYM,I,2),
     &                  NBORB(ISYM,I,1),NBORB(ISYM,I,2),
     &                  1,IPQTOQ(1,0),LUPRI)
           ENDIF
          IVMT = IVMT + NBORB(ISYM,I,1)*NBORB(ISYM,I,2)
          IBOL = IBOL + NBORB(ISYM,I,1)
        ELSE
          NBORB(ISYM,I,2) = 0
        ENDIF
        NBORB(ISYM,I,0) = NBORB(ISYM,I,1) + NBORB(ISYM,I,2)
        ENDDO
        NPSH(I)  = NEFF
        NORB(I)  = NESH(I) + NPSH(I)
        NTMO(I)  = NORB(I)
C
        CALL DZERO(TBUF(IBUF),NFORB(I,0)*NORB(I)*NZ)
        IBO  = 0
        IBOL = 0
        IVMT = KVMT
        DO ISYM = 1,NSYM
        IF(NBORB(ISYM,I,1).GT.0) THEN
          IREP = JFSYM(ISYM,I) - 1
C
C         Small component part
C
          JBUF = IBUF + NFORB(I,0)*IBO + NFORB(I,1)
          JSL  = IMAT + IBOL * NFORB(I,0) + NFORB(I,1)
          CALL QGEMM(NFORB(I,2),NBORB(ISYM,I,2),NBORB(ISYM,I,1),D1,
     &         'N','N',IPQTOQ(1,0),TMAT(JSL),NFORB(I,0),NFORB(I,0),NZ,
     &         'N','N',IPQTOQ(1,0),WORK(IVMT),
     &                 NBORB(ISYM,I,1),NBORB(ISYM,I,2),1,
     &         D0,IPQTOQ(1,0),TBUF(JBUF),NFORB(I,0),NORB(I),NZ)
          IBO = IBO + NBORB(ISYM,I,2)
C
C         Large component part
C
          JBUF = IBUF + NFORB(I,0)*IBO + IBOL
          DO J = 1, NBORB(ISYM,I,1)
            TBUF(JBUF) = D1
            JBUF = JBUF + NFORB(I,0) + 1
          ENDDO
          IBO = IBO + NBORB(ISYM,I,1)
C
C         Update offsets
C
          IBOL = IBOL + NBORB(ISYM,I,1)
          IVMT = IVMT + NBORB(ISYM,I,1)*NBORB(ISYM,I,2)
        ENDIF
        ENDDO
        IF (IPRHAM.GE.6) THEN
         CALL HEADER('Transformation matrix in orthon. basis',-1)
         CALL PRQMAT(TBUF(IBUF),NFORB(I,0),NORB(I),
     &                 NFORB(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
        IMAT = IMAT + NFORB(I,0)*NFORB(I,0)*NZ
        IBUF = IBUF + NFORB(I,0)*NORB(I)*NZ
        CALL MEMREL('MOD_DIREQ',WORK,KSMT,KSMT,KFREE,LFREE)
   10 CONTINUE
C
      IF(VEXTPJ.OR.FREEPJ) THEN
C     These variables are later set by SETDC2, but needed here...
        N2TMT  = 0
        NORBT  = 0
        N2ORBT = 0
        N2TMOTQ = 0
        DO IFRP = 1,NFSYM
          I2TMT(IFRP)  = N2TMT
          N2TMT        = N2TMT + NFORB(IFRP,0)*NORB(IFRP)*NZT
          IORB(IFRP)   = NORBT
          NORBT        = NORBT  + NORB(IFRP)
          I2ORBT(IFRP) = N2ORBT*NZ
          N2ORBT       = N2ORBT + NFBAS(IFRP,0)*NORB(IFRP)
          I2TMOT(IFRP) = N2TMOTQ
          N2TMO(IFRP)  = NTMO(IFRP)*NTMO(IFRP)
          N2TMOTQ      = N2TMOTQ + NTMO(IFRP)*NTMO(IFRP)*NZ
        ENDDO
        CALL MEMGET2('REAL','BUF'  ,KBUF,N2ORBT*NZ,WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IBEIG',KBEIG,NTBAS(0),WORK,KFREE,LFREE)
        CALL DELPOS(TMAT,TBUF,WORK(KBUF),EIG,WORK(KBEIG),
     &              VMAT,WORK,KFREE,LFREE)
        CALL MEMREL('MOD_DIREQ.delpos',WORK,KBUF,KBUF,KFREE,LFREE)
      ELSE
C
C     Generate transformation matrix
C     ==============================
C
        IMAT = 1
        IBUF = 1
        IOFV = 1
        DO 20 I = 1,NFSYM
          IF(NFORB(I,0).EQ.0) GOTO 20
          CALL BCKTR1(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                TBUF(IBUF),NFORB(I,0),NORB(I),
     &                NFORB(I,0),NZ,
     &                NORB(I),1,NFBAS(I,0),
     &                VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,
     &                IPRGEN)
          IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
          IBUF = IBUF + NFORB(I,0)*NORB(I)*NZ
          IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)
   20   CONTINUE
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck psipr1 */
      SUBROUTINE PSIPR1(NVECS,CSEL,ESEL,ISEL,TMAT,TINV,NSTR,
     &                  INFOSEL,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C      Reduce variational space by projection using the fragment orbitals
C      On input:
C        NVECS(ifrp) is the number of fragment orbitals in ifrp
C        CSEL, ESEL, ISEL contains fragment orbitals (coefficients), 
C           their energies and, if available, supersymmetry labels
C        TMAT is the AOMO transformation matrix T (canonical orhogonalization)
C        TINV is the inverse of T^{dagger}, that is TINV = ST, where S is the
C           AO overlap matrix
C
C      On input, for a given fermion ircop, the AOMO transformation matrix
C      has the dimension
C         TMAT(nfbas(ifrp,0),norb(ifr),nz)
C      On output, the dimension has been modified to
C         TMAT(nfbas(ifrp,0),ntmo(ifr),nz)
C   
C
C      Written by Trond Saue long time ago....
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
C
      LOGICAL TOBE
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbdhf.h"
      integer, allocatable :: NVECS_SUB(:)
      DIMENSION CSEL(*),ESEL(*),ISEL(*),TMAT(*),TINV(*),
     &          INFOSEL(2,*),WORK(*)
      DIMENSION NVECS(2),NQ(2),NSTR(2,0:2,MFRAG)
C
      CALL QENTER('PSIPR1')
      KFRSAV = KFREE
      NFROT = 0
      DO IFRP = 1,NFSYM
        NFROT = NFROT + NFRO(IFRP)
      ENDDO
      IF(NFROT.GT.0) CALL OPNFIL(LUCOEF,'DFFROZ','UNKNOWN','PSIPR1')
C.....For simplification , set subblock information, also 
C     when there is no subblock ... allows single code
      IF(SUB_BL) THEN
        NZMO = 1 ! in the spinfree and linear case arrays in orthonormal basis are real
      ELSE
        NZMO = NZ
        DO IFRP = 1,NFSYM
          N_SUB_BL(IFRP) = 1
          NORB_SUB(1,IFRP,0) = NORB(IFRP)
          ID_SUB_BL(1,IFRP)  = 0
        ENDDO
      ENDIF
      TOBE    = .FALSE.
      N2TMT   = 0
      N2TMOTQ = 0      
      ICOFF   = 0
      IEOFF   = 0
      IVOFF   = 0
      ITEFF = 1
      DO 10 IFRP =  1,NFSYM
        WRITE(LUPRI,'(A,A3,A2)') '* Fermion ircop ',FREP(IFRP),' :'
        NQ(1) = 0
        NQ(2) = 0
        IF(IPRINT.GE.5) THEN
          IF(NVECS(IFRP).GT.0) THEN
            WRITE(6,*) '* Energies and symmetry labels of coefficients:'
            DO IVEC = 1,NVECS(IFRP)
              WRITE(6,*) IVEC,ESEL(IVOFF+IVEC),ISEL(IVOFF+IVEC)
            ENDDO
            WRITE(6,*)
     &        '  Input selected coefficients: '
            CALL PRQMAT(CSEL(ICOFF+1),
     &              NFBAS(IFRP,0),NVECS(IFRP),
     &              NFBAS(IFRP,0),NVECS(IFRP),
     &              NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
        ENDIF
        allocate(NVECS_SUB(N_SUB_BL(IFRP)))
C.......First find how the selected vectors are distributed on subblocks
        NREF = 0
        DO ISUB = 1,N_SUB_BL(IFRP)
        IF(NORB_SUB(ISUB,IFRP,0).GT.0) THEN
C.........Find orbitals in this subblock
          NVECS_SUB(ISUB) = 0
          DO IVEC = 1,NVECS(IFRP)
          IF(ISEL(IVEC+IVOFF).EQ.ID_SUB_BL(ISUB,IFRP)) THEN
            NVECS_SUB(ISUB) = NVECS_SUB(ISUB) + 1
            IF(INFOSEL(1,IVOFF+IVEC).GT.0) THEN
C..............Positive energy solution                  
               NTMO_SUB(ISUB,IFRP,1) = NTMO_SUB(ISUB,IFRP,1) - 1
               NQ(1) = NQ(1) + 1
            ELSE
C..............Negative-energy solution                  
               NTMO_SUB(ISUB,IFRP,2) = NTMO_SUB(ISUB,IFRP,2) - 1
               NQ(2)= NQ(2) + 1
            ENDIF
          ENDIF
          ENDDO
          NREF = NREF+NORB_SUB(ISUB,IFRP,0)*NVECS_SUB(ISUB)*NZMO
        ENDIF
        ENDDO
C.......Allocate an array to hold all vectors
        CALL MEMGET2('REAL','REF',KREF,NREF,WORK,KFREE,LFREE)
C.......Transform vector of each subblock to orthonormal basis
        IBORB = 0
        IREF  = KREF
        DO ISUB = 1,N_SUB_BL(IFRP)
        IF(NVECS_SUB(ISUB).GT.0) THEN
C.........Transform vector of this subblock to orthonormal basis
          ITM  = I2TMT(IFRP) + NFBAS(IFRP,0)*IBORB + 1 ! old offset for transformation
          IMO  = ICOFF + 1
          JREF = IREF
          DO IVEC = 1,NVECS(IFRP)
            IF(ISEL(IVEC+IVOFF).EQ.ID_SUB_BL(ISUB,IFRP)) THEN
               CALL QGEMM(NORB_SUB(ISUB,IFRP,0),1,NFBAS(IFRP,0),D1,
     &           'H','N',IPQTOQ(1,0),TINV(ITM),
     &                   NFBAS(IFRP,0),NORB(IFRP),NZT,
     &           'N','N',IPQTOQ(1,0),CSEL(IMO),
     &                   NFBAS(IFRP,0),NVECS(IFRP),NZ,
     &           D0,IPQTOQ(1,0),WORK(JREF),NORB_SUB(ISUB,IFRP,0),
     &           NVECS_SUB(ISUB),NZMO)
               JREF = JREF + NORB_SUB(ISUB,IFRP,0)
            ENDIF
            IMO = IMO + NFBAS(IFRP,0)
          ENDDO
          IF(IPRINT.GE.5) THEN
             WRITE(6,*)
     &       '  Transformed coefficients for subblock ',ISUB
              CALL PRQMAT(WORK(IREF),
     &              NORB_SUB(ISUB,IFRP,0),NVECS_SUB(ISUB),
     &              NORB_SUB(ISUB,IFRP,0),NVECS_SUB(ISUB),
     &              NZMO,IPQTOQ(1,0),LUPRI)
          ENDIF
          IREF = IREF + NORB_SUB(ISUB,IFRP,0)*NVECS_SUB(ISUB)*NZMO
        ENDIF             
        IBORB = IBORB + NORB_SUB(ISUB,IFRP,0)
        ENDDO
C.......We have now transformed all vectors of this irrep to orthonormal basis
C.......and overwrite TINV with TMAT        
        IOFF = I2TMT(IFRP)  + 1
        NDIM = NFBAS(IFRP,0)*NORB(IFRP)
        DO IZ = 1,NZ
          CALL DCOPY(NDIM,TMAT(IOFF),1,TINV(IOFF),1)
          IOFF = IOFF + NDIM
        ENDDO
C.......We now start accumulating the new transformation matrix for each subblock
        NTMO(IFRP) = NORB(IFRP) - NVECS(IFRP)
        WRITE(LUPRI,'(2X,A,I6)')
     &     '- Total number of orbitals deleted: ',NVECS(IFRP),
     &     '       - positive energy solutions: ',NQ(1),
     &     '       - negative energy solutions: ',NQ(2)
        NESH(IFRP)   = NESH(IFRP)-NQ(1) 
        NPSH(IFRP)   = NPSH(IFRP)-NQ(2) 
        NESHMF(IFRP) = NESH(IFRP)
        NPSHMF(IFRP) = NPSH(IFRP)          
C
        IBORB      = 0
        IREF       = KREF
        DO ISUB = 1,N_SUB_BL(IFRP)
        IF(NORB_SUB(ISUB,IFRP,0).GT.0) THEN
          ITM   = I2TMT(IFRP) + NFBAS(IFRP,0)*IBORB + 1 ! old offset for transformation
          NEFF  = NORB_SUB(ISUB,IFRP,0) - NVECS_SUB(ISUB)
C.........New matrix
          IF(NVECS_SUB(ISUB).GT.0) THEN
            IF(SUB_BL) THEN
              ID = ID_SUB_BL(ISUB,IFRP)
              IF (SPINFR) THEN
                WRITE(LUPRI,'(/2A)') '* Boson symmetry ',REP(ID)
              ELSE
                WRITE(LUPRI,'(/A,I4,3A,I2,A)')
     &              '* Block',ISUB,' in ',FREP(IFRP),
     &              ':  Omega = ',ABS(ID),'/2'
              ENDIF
            ENDIF
            NPRJ = NORB_SUB(ISUB,IFRP,0)*NORB_SUB(ISUB,IFRP,0)*NZMO
            CALL MEMGET2('REAL','BUF',KBUF,NPRJ,WORK,KFREE,LFREE)
C...........Do projection
            CALL PRJORT(WORK(KBUF),WORK(IREF),NORB_SUB(ISUB,IFRP,0),
     &           NVECS_SUB(ISUB),NZMO,WORK,KFREE,LFREE)
            IF(IPRINT.GE.5) THEN
            WRITE(6,*)
     &         'PSIPR1: Eigenvectors from PRJORT for subblock ', ISUB
            CALL PRQMAT(WORK(KBUF),
     &              NORB_SUB(ISUB,IFRP,0),NORB_SUB(ISUB,IFRP,0),
     &              NORB_SUB(ISUB,IFRP,0),NORB_SUB(ISUB,IFRP,0),
     &              NZMO,IPQTOQ(1,0),LUPRI)
            ENDIF            
C...........Backtransform new transformation matrix (TINV now contains this part of TMAT)
           CALL QGEMM(NFBAS(IFRP,0),NEFF,NORB_SUB(ISUB,IFRP,0),D1,
     &     'N','N',IPQTOQ(1,0),TINV(ITM),NFBAS(IFRP,0),NORB(IFRP),NZT,
     &     'N','N',IPQTOQ(1,0),WORK(KBUF),NORB_SUB(ISUB,IFRP,0),
     &     NORB_SUB(ISUB,IFRP,0),NZMO,D0,
     &     IPQTOQ(1,0),TMAT(ITEFF),NFBAS(IFRP,0),NTMO(IFRP),NZ)
            CALL MEMREL('PSIPRJ.prj',WORK,KFRSAV,KBUF,KFREE,LFREE)
            IREF  = IREF + NORB_SUB(ISUB,IFRP,0)*NVECS_SUB(ISUB)*NZMO
            TOBE = .TRUE.
          ELSE
            IOFF = ITM
            JOFF = ITEFF
            NDIM = NFBAS(IFRP,0)*NORB_SUB(ISUB,IFRP,0)
            DO IZ = 1,NZ
              CALL DCOPY(NDIM,TINV(IOFF),1,TMAT(JOFF),1)
              IOFF = IOFF + NFBAS(IFRP,0)*NORB(IFRP)
              JOFF = JOFF + NFBAS(IFRP,0)*NTMO(IFRP)
            ENDDO
          ENDIF ! IF(NVECS_SUB(ISUB).GT.0) THEN
          IF(IPRINT.GE.5) THEN
            WRITE(6,*) '* PSIPR1: Transformation matrix..',
     &               ID_SUB_BL(ISUB,IFRP)
            CALL PRQMAT(TMAT(ITEFF),
     &            NFBAS(IFRP,0),NEFF,
     &            NFBAS(IFRP,0),NTMO(IFRP),
     &            NZ,IPQTOQ(1,0),LUPRI)
          ENDIF            
          IBORB = IBORB + NORB_SUB(ISUB,IFRP,0)
          ITEFF = ITEFF + NFBAS(IFRP,0)*NEFF   ! new offset for transformation
          NTMO_SUB(ISUB,IFRP,0) = NTMO_SUB(ISUB,IFRP,1)
     &                          + NTMO_SUB(ISUB,IFRP,2)
          NORB_SUB(ISUB,IFRP,0) = NTMO_SUB(ISUB,IFRP,0) 
          NORB_SUB(ISUB,IFRP,1) = NTMO_SUB(ISUB,IFRP,1)
          NORB_SUB(ISUB,IFRP,2) = NTMO_SUB(ISUB,IFRP,2)            
        ENDIF
        ENDDO                    ! DO ISUB = 1,N_SUB_BL(IFRP)
        deallocate(NVECS_SUB)
C
C       Save frozen orbitals to file
C
        IF(NFRO(IFRP).GT.0) THEN
          NBUF = NFBAS(IFRP,0)*NZ
          CALL MEMGET2('INTE','SEL',KSEL,NVECS(IFRP),WORK,KFREE,LFREE)
          CALL MEMGET2('INTE','IND',KIND,NFRO(IFRP),WORK,KFREE,LFREE)
          CALL MEMGET2('REAL','BUF',KBUF,NBUF,WORK,KFREE,LFREE)
          CALL WRIFRZ(LUCOEF,IFRP,CSEL(ICOFF+1),ESEL(IVOFF+1),
     &         ISEL(IVOFF+1),WORK(KSEL),WORK(KBUF),WORK(KIND),
     &         NVECS(IFRP),NQ)
          CALL MEMREL('PSIPRJ.froz',WORK,KFRSAV,KBUF,KFREE,LFREE)
       ENDIF
       NORB(IFRP)   = NTMO(IFRP) + NFRO(IFRP)
       I2TMT(IFRP)  = N2TMT
       N2TMT        = N2TMT + NFBAS(IFRP,0)*NTMO(IFRP)*NZT
       I2TMOT(IFRP) = N2TMOTQ
       N2TMO(IFRP)  = NTMO(IFRP)*NTMO(IFRP)
       N2TMOTQ      = N2TMOTQ + N2TMO(IFRP)*NZ
       ICOFF = ICOFF + NFBAS(IFRP,0)*NVECS(IFRP)*NZ  ! update offset for selected coefficients
       IEOFF = IEOFF + NVECS(IFRP) ! update offset for energies and supersymmetries
 10   CONTINUE
      IF(NFROT.GT.0) CLOSE(LUCOEF,STATUS='KEEP')
      IF(IPRHAM.GE.5) THEN
        CALL HEADER('PSIPRJ: Final transformation matrix ',-1)
        DO 20 IFRP = 1,NFSYM
        IF(NTMO(IFRP).GT.0) THEN
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',IFRP
          CALL PRQMAT(TMAT(I2TMT(IFRP)+1),NFBAS(IFRP,0),NTMO(IFRP),
     &                NFBAS(IFRP,0),NTMO(IFRP),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
 20     CONTINUE
      ENDIF
      IF(TOBE) THEN
C
C       Construct right index transformed overlap matrix Saomo = Sao T
C       Since T(dagger) Sao T = I, Saomo is the inverse of T(dagger)
C
        CALL MKSAOMO(TINV,TMAT,.FALSE.,IPRHAM,WORK(KFREE),LFREE)
      ENDIF
C
      CALL QEXIT('PSIPR1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Get_Overlap_AO */
      SUBROUTINE GET_OVERLAP_AO(SMAT,SSMTRC,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Get overlap matrix in AO-basis
C     Written by Trond Saue June 4 2012
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
      DIMENSION SMAT(*),WORK(*)
      KFRSAV=KFREE
C.....Get overlap matrix in SO-basis (no symmetry packing)
      CALL MEMGET2('REAL','SBUF',KSBUF,N2BBASX,WORK,KFREE,LFREE)
      CALL GTOVLX(WORK(KSBUF),SSMTRC)
C.....Transform to AO-basis and release memory
      CALL MTSOAO(WORK(KSBUF),SMAT,NTBAS(0),0,IPRINT)
      CALL MEMREL('GET_OVERLAP_AO',WORK,1,KFRSAV,KFREE,LFREE)
C.....Print section
      IF(IPRINT.GE.5) THEN
        WRITE(6,*) 'GET_OVERLAP_AO: Overlap matrix in AO-basis'
        CALL OUTPUT(SMAT,1,NTBAS(0),1,NTBAS(0),NTBAS(0),NTBAS(0),
     &                   1,LUPRI)
      ENDIF
C      
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
      SUBROUTINE SET_INTFLG(new_intflg)
C***********************************************************************      
      use mxcent_mod
#include "implicit.h"
      integer, intent(in) :: new_intflg
#include "mxgas.h"
#include "dcbcosc.h"
#include "dcbdhf.h"
#include "dcbgrd.h"
#include "dcbkrci.h"
#include "dcbmp2.h"
#include "dcbopt.h"
#include "../prp/dcbxqr.h"
      ICOSC_INTDEF = new_intflg
      ICOSC_INTFLG = new_intflg
      ICOSC_INTBUF = new_intflg

      INTDEF       = new_intflg
      INTDEF_SAVE  = new_intflg
      INTFLG       = new_intflg
      INTBUF       = new_intflg

      IGRD_INTFLG  = new_intflg
      IGRD_INTBUF  = new_intflg

      iKRCI_INTDEF = new_intflg
      iKRCI_INTFLG = new_intflg
      iKRCI_INTBUF = new_intflg

      MP2_INTFLG   = new_intflg

      iOPT_INTDEF = new_intflg
      iOPT_INTFLG = new_intflg
      iOPT_INTBUF = new_intflg

      intxqr      = new_intflg

!     intxpp      = new_intflg

      call set_intfl1(new_intflg)
      call set_intfl2(new_intflg)

      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
      SUBROUTINE SET_INTFL1(NEW_INTFLG)
C***********************************************************************      
      use dcbnmr_mod, only: INTNMR
      implicit none
      integer, intent(in) :: new_intflg
#include "mxgas.h"
#include "dcbxrs.h"
#include "dcbkrmc.h"
#include "dcbxlr.h"
      INTXRS       = new_intflg
      INTDEF       = new_intflg
      INTFLG       = new_intflg
      INTBUF       = new_intflg

      iKRMC_INTDEF = new_intflg
      iKRMC_INTFLG = new_intflg
      iKRMC_INTBUF = new_intflg

      INTXLR       = new_intflg

      INTNMR       = new_intflg

      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
      SUBROUTINE SET_INTFL2(NEW_INTFLG)
C***********************************************************************      
      implicit none
      integer, intent(in) :: new_intflg
#include "dcbgen.h"

!     INTGEN       = new_intflg
      ILLDIR       = new_intflg
      ISLDIR       = 0
      ISSDIR       = 0
      IGTDIR       = 0
!     IDFLAG       = 1

      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck nossnucatt */
      SUBROUTINE NOSSNUCATT
C***********************************************************************
C     
C     Zero SS-block of potentials
C
C***********************************************************************
#include "dcbham.h"
#include "dcbxpr.h"      
#include "dcbprl.h"      
C.....Nuclear attraction integrals
      IF(IPMOLFLD.GT.0) THEN
         ILBL  = IPRPLBL(1,IPMOLFLD)
           PDOINT(ILBL)(4:4) = '0'
       ENDIF
C.....Static embedding potential
      IF(IPVEMB0.GT.0) THEN
         ILBL  = IPRPLBL(1,IPVEMB0)
           PDOINT(ILBL)(4:4) = '0'
       ENDIF
       END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prjort */
      SUBROUTINE PRJORT(VMO,COEF,NORB,NVEC,NZ,WORK,KFREE,LFREE)
C***********************************************************************
C     This routine projects vectors COEF out of the identity matrix
C     Presently, this is done by first forming the density matrix 
C     of vectors COEF and then diagonalize it.
C     Provided that eigenvalues are given in ascending order, then
C     the first NORB-NVEC columns of the resulting eigenvectors have
C     zero eigenvalues and provide the orthogonal complement.
C     A possible better approach is to do a "fat" QR decomposition.
C     Thanks to Ulf Ekstrom for suggestions.
C
C     Written by Trond Saue
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
      DIMENSION VMO(NORB,NORB,NZ),COEF(NORB,NVEC,NZ),WORK(*)
C.....Form density matrix
      KFRSAV = KFREE
      CALL MEMGET2('REAL','DMAT',KDMAT,NORB*NORB*NZ,WORK,KFREE,LFREE)
      CALL DENST1(WORK(KDMAT),NORB,NORB,NZ,D1,D0,
     &     COEF,NORB,NVEC,1,NVEC,NORB)
      CALL MEMGET2('REAL','EIG',KEIG,NORB,WORK,KFREE,LFREE)
      CALL QDIAG(NZ,NORB,WORK(KDMAT),NORB,NORB,
     &           WORK(KEIG),1,VMO,NORB,
     &           NORB,WORK(KFREE),LFREE,IERR)
      WRITE(6,*) ' Eigenvalues from PRJORT: ',NVEC
      CALL OUTPUT(WORK(KEIG),1,NORB,1,1,
     &            NORB,NORB,-1,LUPRI)
      CALL MEMREL('PRJORT',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE Add_Screening_Potential(FMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Replacement for BNCORR, add screening potential by reading in fitted local
C     potential from GRASP DFT calculations. 
C
C     Lucas Visscher, July 2016
C
C     For each center in a molecule we estimate the screening
C     neglected in the bare nucleus approximation by a sum of the
C     contributions from Gaussian functions :
C           sum_(n) <X_A|(c(n)*G(n,a,r_C))/r_C|X_B>
C           X=L,S ; G(a,r)=exp(-a*r^2)
C
C     On input/ouput:     FMAT - assumed one-electron (2c/4c bare nucleus) Fock matrix (input),
C    =================    screening bare nucleus added (ouput)
C
C     Called from: PREDHF/dirscf.F
C                  RH1DIAG/../krmc/krmcopt.F
C
C
C*****************************************************************************
      use primitives_mod
      use maxaqn_mod
      use mxcent_mod
      use symmet_mod
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0)
C
      LOGICAL DOINT(2,2), FINISH
      DIMENSION WORK(LWORK),FMAT(*)
      real*8, allocatable  :: ONEMT(:),ONEIN(:),FCBA(:,:,:),COORC(:,:)
      real*8, allocatable  :: SIGNC(:,:),GEXP(:)
      integer, allocatable :: NCENT(:),JSYMC(:),JCENT(:)
      real*8, allocatable  :: fit_exp(:),fit_coef(:)
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "cbihr1.h"
#include "cbisol.h"
C
#include "ibtfun.h"
c
#include "memint.h"

        if(mdirac)then
        write(lupri,*)
     &  ' bare-nucleus correction not implemented yet for MDIRAC'
        write(lupri,*)
     &  ' program will continue with the default 1-el Hamiltonian'
          return
        end if

C     Memory allocation
      Allocate (ONEMT(N2BBASX))
      Allocate (ONEIN(N2BBASX))
      Allocate (COORC(3,NUCDEP))
      Allocate (SIGNC(3,NUCDEP))
      Allocate (GEXP(NUCDEP))
      Allocate (NCENT(NUCDEP))
      Allocate (JSYMC(NUCDEP))
      Allocate (JCENT(NUCDEP))

      ONEMT = D0
C
C     Initialization (calculate LL and SS)
C
      DOINT(1,1) = .TRUE.
      IF (NOSMLV.OR.BSS.OR.x2c) THEN
         DOINT(2,2) = .FALSE.
         IF (IPRINT.GE.1) WRITE(LUPRI,'(/A)')
     &     '  SCRPOT: SS contribution to the bare nuclei omitted !'
      ELSE
         DOINT(2,2) = .TRUE.
      END IF
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.
c

      DO I = 1,NUCIND
         IF (NOORBT(I) .OR. CHARGE(I).EQ.0.0D0) THEN ! skip point charges and floating orbitals
            CYCLE
         ELSE
c        ... Get charges and exponents for the electrons in this group
            n_charge = nint(CHARGE(I))
            ! do no use atom number in IZATOM(I) here
            ! Reason: if ECP we need the reduced charge to get a reasonable start guess ...
            n_exp = -1 ! signal the routine to provide the number of exponents only
            CALL Read_screening_parameters
     &           (n_charge,n_exp,DUMMY,DUMMY)
            Allocate (fit_exp(n_exp))
            Allocate (fit_coef(n_exp))
            CALL Read_screening_parameters 
     &           (n_charge,n_exp,fit_exp,fit_coef)
         END IF

c        The driver routines below cannot take a generally contracted potential, so we need to loop. 
c        if this becomes a bottleneck: rewrite the catdr1 routine....

         MXCENTC = 0
         MULC   = ISTBNU(I)
         CORCX0 = CORD(1,I)
         CORCY0 = CORD(2,I)
         CORCZ0 = CORD(3,I)
c
c        ... Count number of symmetry-dependent nuclei so we can allocate FCBA
c
         DO ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,MULC) .EQ. 0) THEN
               MXCENTC = MXCENTC + 1
            END IF
         end do
         Allocate (FCBA(MXCENTC,NUCDEP,NUCDEP))

         DO N = 1, n_exp
            NCENTC = 0
            CHARG1 = - fit_coef(N)
            DO ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,MULC) .EQ. 0) THEN
                  NCENTC = NCENTC + 1
                  JSYMC(NCENTC)   = ISYMOP
                  JCENT(NCENTC)   = I
                  SIGNC(1,NCENTC) = PT(IBTAND(ISYMAX(1,1),ISYMOP))
                  SIGNC(2,NCENTC) = PT(IBTAND(ISYMAX(2,1),ISYMOP))
                  SIGNC(3,NCENTC) = PT(IBTAND(ISYMAX(3,1),ISYMOP))
                  COORC(1,NCENTC) = SIGNC(1,NCENTC)*CORCX0
                  COORC(2,NCENTC) = SIGNC(2,NCENTC)*CORCY0
                  COORC(3,NCENTC) = SIGNC(3,NCENTC)*CORCZ0
                  GEXP(NCENTC)    = fit_exp(N)
                  CALL DCOPY(NUCDEP*NUCDEP,CHARG1,0,
     &                       FCBA(NCENTC,1,1),MXCENTC)
               END IF
            end do
c
c       ... Calculate integrals (both LL and SS)
c
            CALL CATDR1(ONEMT,DUMMY,NCENTC,FCBA,
     &                  COORC,SIGNC,JSYMC,
     &                  JCENT,GEXP,WORK,LWORK,IPRONE,
     &                  .FALSE.,0,NNBBASX,DOINT)
         ENDDO ! Loop over exponents
         DeAllocate (fit_exp)
         DeAllocate (fit_coef)
         DeAllocate (FCBA)
      ENDDO ! Loop over nuclei
C
C     Symmetry unpack integrals
C
      CALL SYMUPK(ONEMT,ONEIN,1,NNBBASX)
      CALL DSPTSI(NTBAS(0),ONEIN,ONEMT)
      CALL BUTOBS(ONEMT,1,WORK,LWORK)
 
C     Add correction to (entering) Fock-matrix, FMAT

      CALL DAXPY(N2BBASX,D1,ONEMT,1,FMAT,1)

      IF(IPRINT.GE.11) THEN
         CALL HEADER('Screening potential for DC FMAT',-1)
         CALL OUTPUT(ONEMT,1,NTBAS(0),1,NTBAS(0),
     &        NTBAS(0),NTBAS(0),-1,LUPRI)
      END IF
C
C     Memory deallocation
      DeAllocate (ONEMT)
      DeAllocate (ONEIN)
      DeAllocate (COORC)
      DeAllocate (SIGNC)
      DeAllocate (GEXP)
      DeAllocate (NCENT)
      DeAllocate (JSYMC)
      DeAllocate (JCENT)

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE READ_SCREENING_PARAMETERS (N_CHARGE,N_EXP,
     &     FIT_EXP,FIT_COEF)
C***********************************************************************      
      implicit none
      integer,intent (in)    :: n_charge
      integer,intent (inout) :: n_exp
      real*8,intent (out)    :: fit_exp(*),fit_coef(*)

      integer                :: i,intorb, lupot, lu_not_used
      logical                :: newel, early_exit, emsl_type
      integer, parameter     :: iprint=0

      lupot  = lu_not_used(91)
      early_exit = n_exp .eq. -1 ! Use this to query the number of exponents
      call find_element ('SCRPOT',n_charge,lupot,iprint,emsl_type)
      call find_pos (newel,n_exp,intorb,lupot,iprint,emsl_type)

      if (.not. early_exit) then
         do i = 1, n_exp
            read (lupot,*) fit_exp(i),fit_coef(i)
         end do
      end if

      close (lupot, STATUS = 'KEEP')

      end
#ifdef DEBUG_SOC
#undef DEBUG_SOC
#endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck ashift */
      SUBROUTINE ASHIFT(FOCK,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Build atomic shift operator
C
C     Written by Trond Saue
C
C***********************************************************************
      use mxcent_mod
#include "implicit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
C     
C
#include "dcbham.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcb_ashift.h"      
      DIMENSION FOCK(*),WORK(*)
      CHARACTER TEXT*74
      KFRSAV = KFREE
      NDIM=N2BBASX*NZ
      CALL MEMGET2('REAL','KBUF',KBUF,NDIM,WORK,KFREE,LFREE)      
      CALL ATOMIC_DENS(WORK(KBUF),MAXAREF,MAXAVAL,AREFFIL,AVECREF,
     &     ASHIFT_NVAL,ASHIFT_PAR,IPRINT,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','KSMAT',KSMAT,N2BAST,WORK,KFREE,LFREE)
      CALL GTOVLT(WORK(KSMAT),SSMTRC,IPRINT)
      DO I = 1,NFSYM
         CALL QTRANS('AOMO','S',D1,
     &              NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),
     &              WORK(I2BASX(I,I)+KBUF),NTBAS(0),NTBAS(0),
     &                 NZ,IPQTOQ(1,0),
     &              FOCK(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),
     &                 NZ,IPQTOQ(1,0),
     &              WORK(I2BAST(I)+KSMAT),NFBAS(I,0),NFBAS(I,0),1,1,
     &              WORK(I2BAST(I)+KSMAT),NFBAS(I,0),NFBAS(I,0),1,1,      
     &              WORK(KFREE),LFREE,IPRINT)
      ENDDO
      CALL MEMREL('ASHIFT',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BSS_RESTART(VMAT,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Restart BSS run from file
C
C     BSSMAT available:
C     1) open BSSMAT and remove all property integrals
C     2) do picture change transformation of integrals on CHECKPOINT
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      logical :: FNDLAB
      integer :: IPRINT
      DIMENSION VMAT(*),WORK(*)
      WRITE(LUPRI,'(A,A)') '* Info: DIRAC found BSSMAT ',
     &      '(BSS stuff) ... skipping generation of H2C'
      LUBSS = 40
      CALL OPNFIL(LUBSS,'BSSMAT','OLD','BSS_RESTART')
      IF (FNDLAB('P2C_0001',LUBSS)) THEN
        BACKSPACE(LUBSS)
        CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
      ENDIF
      CLOSE(LUBSS,STATUS='KEEP')
      CALL PROP2BSS(VMAT,IPRINT,WORK(KFREE),LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE X2C_RESTART(IPRINT)
C***********************************************************************
C
C     Restart X2C run from file
C
C     X2CMAT available:
C     1) open X2CMAT and remove all property integrals
C     2) do picture change transformation of integrals on CHECKPOINT
C***********************************************************************
      use memory_allocator
      use checkpoint      
      use x2cmod_cfg      
      use x2c_pct_ao, only:
     &    pctrafo_all_op_driver
      use quaternion_algebra

      use maxaqn_mod
      use ccom_mod
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbdhf.h"
#include "dcbcls.h"
#include "dcbxpr.h"            
      real(8), allocatable :: WORK(:)
      logical              :: x2c_AOA,x2c_Aon
      integer              :: iprint
C      
      write(lupri,'(a,a)') '* Info: DIRAC found X2CMAT ',
     &    '(X2C stuff) ... skipping generation of H2C'
      open(lux2c,file='X2CMAT',status='old',form='unformatted',
     &             access='sequential',action='readwrite',
     &             position='rewind')
!     find Hamiltonian on file property integral... and step back
!     I2COFK     - integer number specifying which 2comp. Hamiltonian (in AO basis)
!                  is to be read in the ONEFCK routine
      if(.not.x2cmod_fragment_x2c)then
        i2cofk =  0
!       check for hamiltonian integrals on file X2CMAT and set i2cofk
        call checkpoint_query(
     &'/result/hamiltonian/x2c/ao_matrices/h12cAOA   10',x2c_AOA)
        call checkpoint_query(
     &'/result/hamiltonian/x2c/ao_matrices/h12cAOn   10',x2c_AOn)
        if(x2c_AOA)then ! with AMFI
          if(noamfi)then
            write(lupri,'(a)')
     &' *** restart error: try to read h1+amfi and noamfi == .true. ***'
            call quit(
     &'*** error in gmotra: mismatch between X2CMAT data and input ***')
          end if
!         set i2cofk to |3|: since we restart, it does not matter
!         how the AMFI part was derived.
          i2cofk = 3
        else if(x2c_AOn)then !without AMFI
          i2cofk = 1
        else
          if(x2cmod_mmf_restart)then
            i2cofk = -2
          else
            write(lupri,'(a)')
     &'*** restart error: no hamiltonian integrals stored on X2CMAT ***'
            call quit(
     &'*** error in gmotra: no h1 data on file X2CMAT ***')
          end if
        end if
      else
!     nothing to do here for now: x2cmod_fragment_x2c is always false.
      end if
      if(.not.x2cmod_skip_op_pct)then
        call pctrafo_all_op_driver(ntbas(0),
     &                             ntbas(1),
     &                             nfbas(1,0),
     &                             nfbas(1,1),
     &                             nfsym,
     &                             nz,
     &                             irqmat,
     &                             ipqtoq,
     &                             iqtopq,
     &                             iqmult,
     &                             jqbas,
     &                             iqdef,
     &                             i2basx,
     &                             lux2c,
     &                             lu1int,
     &                             jbtof,
     &                             nprps,
     &                             lwork,
     &                             dosphe,
     &                             mdirac,
     &                             .false.,
     &                             onesys,
     &                             -1,
     &                             df,
     &                             -1,
     &                             'undef ',
     &                             'undef ',
     &                             iprint)
      end if
      close(lux2c,status='keep')
!     allocate
      call legacy_lwork_get(LWORK)
      call alloc(WORK,LWORK,id='WORK in GMOTRA - part 2')
#include "memint.h"
      TWOCOMP = .TRUE.
      RDINPC = .FALSE.
      IPREAD = IPREAD - 5
      CALL READIN(.FALSE.)
!     > stknecht: take care of 2c-mmf-sf for a linear molecule
      if(spinfr .and. linear) linear = .false.
      if(spinfr .and. atomic) atomic = .false.          
      CALL SETDC1(IPREAD)
!     set AO-labels for large components - in x2c_main initialized for large and small components
      CALL GETLAB(IPREAD)
      IPREAD = IPREAD + 5
C.......reset NEWPRP to re-initialize AO matrices on CHECKPOINT
      NEWPRP = .TRUE.
      CALL LSET(NPRPCLS,.TRUE.,CLSCAL)
      CALL ONEGEN(WORK,LWORK)
      call dealloc(work)
      RETURN
      END
#ifdef MOD_XAMFI      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE xamfi_atomcorr(delta_xc_energy_aoosoc,
     &                          no1esystem_aoosoc)
      use memory_allocator
      use xamfi_environment
      use xamfi_global_parameters
      use xamfi_driver, only: get_xamfi_correction
      use xcint_main  , only: get_xc_energy      
      use dircmo
      use dirac_cfg, only: dirac_cfg_dft_calculation      
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"            
#include "dcbdhf.h"
#include "dgroup.h"      
      real(8), allocatable :: dmat_aoosoc(:)
      real(8), allocatable :: cmo_aoosoc(:)
      real(8), allocatable :: pctm_aoosoc(:)      
      real*8               :: delta_xc_energy_aoosoc
      logical              :: x2c_save      
      logical              :: no1esystem_aoosoc, aoc_save
      character(len=12)    :: flabel
      integer              :: lwork
      integer              :: nfmat_save,ioff_pct_mat,ichrg
!      integer              :: lux2c=40
!      lux2c=40
!>    temporarily switch on 4c mode
      call xamfi_set_cb("get-4c")
!> if requested (.FOCC): density matrix for 2e-SO+scalar corrections
!> always with fractional occupation for open-shell atoms
      if(aoofocc .and. aoc)then
        aoc_save   = aoc
        nfmat_save = nfmat
        aoc        = .false.
        nfmat      = 1
      end if
!> turn off x2c to properly call twofck in 4c mode
      x2c_save = x2c
      x2c      = .false.
!> generate density matrix - incoming variable to the aoosoc module
      allocate(dmat_aoosoc(ntbas(0)**2*nz*nfmat))
      allocate(cmo_aoosoc(ncmotq))
      dmat_aoosoc = 0.0d0
      cmo_aoosoc  = 0.0d0
!     read coefficients
      call reacmo_new(cmo=cmo_aoosoc)
!     generate density matrix
      call denmat(dmat_aoosoc,cmo_aoosoc,iprham)
      allocate(pctm_aoosoc(ntbas(0)**2*4))
!>    read the U matrix
      ioff_pct_mat = 1
      do i = 1, nfsym
        write(flabel,'(a11,i1)') 'pctmtAO   1',i
        if(nfbas(i,0) * nfbas(i,1) > 0)then
          call x2c_read(flabel,pctm_aoosoc(ioff_pct_mat),
     &                  nfbas(i,0) * nfbas(i,1)*nz,lux2c)
        end if
        ioff_pct_mat = ioff_pct_mat + nfbas(i,0) * nfbas(i,1) * nz
      end do
      close(lux2c,status='keep')
      call legacy_lwork_get(lwork)      
      lwork_local = lwork-ntbas(0)**2*nz*nfmat-ntbas(0)**2*nz-ncmotq
      call get_xamfi_correction(
     &                          dmat_aoosoc,
     &                          cmo_aoosoc,
     &                          pctm_aoosoc,
     &                          no1esystem_aoosoc,
     &                          dirac_cfg_dft_calculation,
     &                          (linear.or.atomic),
     &                          lwork_local,
     &                          iprham
     &                          )
!>    save XC_energy_4c
      delta_xc_energy_aoosoc     = 0.0d0
      if(dirac_cfg_dft_calculation)then
        delta_xc_energy_aoosoc =
     &  delta_xc_energy_aoosoc + get_xc_energy()
      end if
!>    reset common block variables to 2c mode
      call xamfi_set_cb("get-2c")
!> restore logicals
      if(aoofocc)then
        aoc   = aoc_save
        nfmat = nfmat_save
      end if
      x2c = x2c_save
!>    release memory      
      deallocate(pctm_aoosoc)
      deallocate(cmo_aoosoc)
      deallocate(dmat_aoosoc)

      end subroutine xamfi_atomcorr
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE xamfi_atomcorr2(delta_xc_energy_aoosoc,
     &     no1esystem_aoosoc)
      use dircmo      
      use dirac_cfg, only: dirac_cfg_dft_calculation
      use x2c_utils, only:
     &     read_fock_matrices_saao_basis_x2c
      use xamfi_global_parameters
      use xamfi_internal_parameters
      use xamfi_environment      
      use xamfi_utils, only:
     &    dump_normalized_2e_soc,
     &    read_xamfi_general_matrix,
     &    dump_delta_xc_energy,
     &    goodbye_xamfi
      use aobasis_collector, only: deallocate_aobasis_container

      use mxcent_mod
#include "dcbbas.h"
#include "dcbham.h"      
#include "dgroup.h"      
#include "nuclei.h"
      real(8), allocatable :: tmat_mmf(:)
      real(8), allocatable :: tbuf_mmf(:)
      real(8), allocatable :: cmo_mmf(:)
      real(8), allocatable :: scratch_mmf(:),loop_factor(:)
      real*8               :: delta_xc_energy_aoosoc
      real*8               :: xxxBLUBB
      integer              :: len_fmat      

      len_fmat  = ntbas(0)**2*nz
      allocate(tmat_mmf(len_fmat))
      allocate(tbuf_mmf(len_fmat))
      allocate(loop_factor(0:1))
      loop_factor = 0.0d0
      open(lux2c,file='D2FCK_mmf',status='unknown',
     &           access='sequential',action='readwrite',
     &           form='unformatted',position='rewind')
!>    read G_A[2c]
      call read_fock_matrices_saao_basis_x2c(
     &                                       tmat_mmf,
     &                                       len_fmat, 1,
     &                                       lux2c, 0,
     &                                       loop_factor)
      deallocate(loop_factor)
      close(lux2c, status='keep')
#ifdef DEBUG_SOC
      CALL HEADER('X-AMFI - G_A[2c]',-1)
      CALL PRQMAT(tmat_mmf,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &            IPQTOQ(1,0),LUPRI)
#endif
!>    read {U+ G_A[4c] U}^{++}
      call read_xamfi_general_matrix(
     &                               tbuf_mmf,
     &                               ntbas(0),
     &                               ntbas(0),
     &                               nz,
     &                               'aoo2esoX',
     &                               nint(charge(1)),
     &                               lux2c,
     &                               1 !> atomic center id
     &                              )
#ifdef DEBUG_SOC
      CALL HEADER('X-AMFI - {U+ G_A[4c] U}^{++}',-1)
      CALL PRQMAT(tbuf_mmf,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &            IPQTOQ(1,0),LUPRI)
#endif
!> G_A^{++} = {U+ G_A[4c] U}^{++} - G_A[2c]
      call daxpy(len_fmat,-1.0d0,tmat_mmf,1,tbuf_mmf,1)
#ifdef DEBUG_SOC
      CALL HEADER('X-AMFI - 2e-normalized (final)',-1)
      CALL PRQMAT(tbuf_mmf,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &            IPQTOQ(1,0),LUPRI)
#undef DEBUG_SOC
#endif
      call dump_normalized_2e_soc(
     &                            tbuf_mmf,
     &                            ntbas(0),
     &                            ntbas(0),
     &                            nz,
     &                            'aoo2esoc',
     &                            nint(charge(1)),
     &                            lux2c,
     &                            .true.
     &                           )
      if(dirac_cfg_dft_calculation)then
#ifdef DEBUG_SOC
        WRITE(LUPRI,'(/a,f12.9)')
     &        'X-AMFI - DELTA XC energy: ',delta_xc_energy_aoosoc
#undef DEBUG_SOC
#endif
        tbuf_mmf = 0; tmat_mmf = 0
!>      read {U+ V_A[4c] U}^{++}
        call read_xamfi_general_matrix(
     &                                 tbuf_mmf,
     &                                 ntbas(0),
     &                                 ntbas(0),
     &                                 nz,
     &                                 'amftVxc2',
     &                                 nint(charge(1)),
     &                                 103,
     &                                 1 !> atomic center id
     &                                )
        open(103,file='amfVxc2',status='old',
     &       access='sequential',action='readwrite',
     &       form='unformatted',position='rewind')
!>      read V_A[2c]
        read(103) tmat_mmf(1:ntbas(0)*ntbas(0)*nz)
        close(103,status="delete")
#ifdef DEBUG_SOC
        CALL HEADER('X-AMFI - Vxc-tilde             ',-1)
        CALL PRQMAT(tbuf_mmf,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &              NZ,IPQTOQ(1,0),LUPRI)
        CALL HEADER('X-AMFI - Vxc2                  ',-1)
        CALL PRQMAT(tmat_mmf,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &            NZ,IPQTOQ(1,0),LUPRI)
#endif
!>      V_A^{++} = {U+ V_A[4c] U}^{++} - V_A[2c]
        call daxpy(len_fmat,-1.0d0,tmat_mmf,1,tbuf_mmf,1)
#ifdef DEBUG_SOC
        CALL HEADER('X-AMFI - Vxc-normalized (final)',-1)
        CALL PRQMAT(tbuf_mmf,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &              NZ,IPQTOQ(1,0),LUPRI)
        allocate(cmo_mmf(ncmotq))
        allocate(scratch_mmf(len_fmats))
        cmo_mmf     = 0
        scratch_mmf = 0
!>      read 2c-mos
        call reacmo_new(cmo=cmo_mmf)
!>      calculate inactive (and possibly active) density matrix
        call denmat(scratch_mmf,cmo_mmf,iprham)
        xxxBLUBB = ddot(len_fmat,scratch_mmf,1,tbuf_mmf,1)
        WRITE(LUPRI,'(/a,f12.9)')
     &        'X-AMFI - DELTA XCmat (PCE) energy: ',xxxBLUBB
        deallocate(scratch_mmf)
        deallocate(cmo_mmf)
#undef DEBUG_SOC
#endif
        call dump_normalized_2e_soc(
     &                              tbuf_mmf,
     &                              ntbas(0),
     &                              ntbas(0),
     &                              nz,
     &                              'aooVxcMA',
     &                              nint(charge(1)),
     &                              103,
     &                              .true.
     &                              )
        call dump_delta_xc_energy(
     &                            delta_xc_energy_aoosoc,
     &                            'aooDLXCE',
     &                            nint(charge(1)),
     &                            lux2c
     &                           )
      end if ! dirac_cfg_dft_calculation
      call goodbye_xamfi(atom(1)%charge)
      call deallocate_aobasis_container()
!>    release memory
      deallocate(tbuf_mmf)
      deallocate(tmat_mmf)
      
      end subroutine xamfi_atomcorr2
#endif      
! -- end of dirone.F --
