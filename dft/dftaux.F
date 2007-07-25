C
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C 
C /* Deck dft_linres */
#if defined (MOD_DNF)
Caspg: extra argument to keep track of b0, used for debugging of df-vxc-lr
      SUBROUTINE DFT_LR(DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,GAO,GAO1,
     &                  RHO,RHOGRD,RHG,WGHT,DOGGA,BUF,B0)
#else
      SUBROUTINE DFT_LR(DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,GAO,GAO1,
     &                  RHO,RHOGRD,RHG,WGHT,DOGGA,BUF)
#endif /* MOD_DNF */
C
C     T. Helgaker sep 99/oct 00
C     adapted to DIRAC by T.Saue May 3 2001
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0,DP5 = 0.5D0,D2=2.0D0)
C
#include "inforb.h"
#include "nuclei.h"
#if defined (PRG_DIRAC)
#include "dgroup.h"
#else
#include "wrkrsp.h"
#include "dftinf.h"
#endif
#include "dftcom.h"
#include "functionals.h"
C
      INTEGER A, B
      LOGICAL DOGGA, TRPLET, LBIT
      DIMENSION DTRMAT(NBAST,NBAST,NZ,NXCMAT),
     &          EXCMAT(NBAST,NBAST,NZ,NXCMAT),
     &          ISYMOP(NXCMAT),IHRMOP(NXCMAT),
     &          GAO(NBAST), GAO1(NBAST,3),
     &          RHG(3),VT(3), VXC(DERVS2), BUF(*)
      DIMENSION B3(3)
C
      DO IM = 1,NXCMAT
      IREPOP = ISYMOP(IM)-1
      IF(IHRMOP(IM).EQ.1) THEN
C
C       Get perturbed density
C       =====================
C
        CALL GETRHO(B0,IREPOP,GAO,DTRMAT(1,1,1,IM),BUF)
#if defined (MOD_DNF)
c       write (*,*) 'debug, in dft_lr: b0',b0
#endif /* MOD_DNF */
C
C       ***************
C       ***** LDA *****
C       ***************
C
        IF(.NOT.DOGGA) THEN
          IF(DABS(B0).GT.DFTHRL) THEN
C
C           Calculated VT
C
            CALL DFTPOT1(VXC, WGHT,RHO, 0D0)
            SXC = VXC(FRR)*B0
C
C           Linear transformation
C
            CALL DOEXC1(SXC,IREPOP,EXCMAT(1,1,1,IM),GAO)
          END IF
        ELSE
C
C         ***************
C         ***** GGA *****
C         ***************
C
          CALL GETGRHO(B3,IREPOP,GAO,GAO1,DTRMAT(1,1,1,IM),BUF)
          BMAX = MAX(ABS(B0),ABS(B3(1)),ABS(B3(2)),ABS(B3(3)))
C
          IF (BMAX.GT.DFTHRL) THEN
            CALL DFTPOT1(VXC,WGHT,RHO,RHOGRD)
            BR = (B3(1)*RHG(1)+ B3(2)*RHG(2)+ B3(3)*RHG(3))/RHOGRD
            FACR = VXC(FRZ)*B0 + (VXC(FZZ)-VXC(FZ0)/RHOGRD)*BR
            SXC = VXC(FRR)*B0 + VXC(FRZ)*BR
            DO I = 1,3
              VT(I) = (VXC(FZ0)*B3(I)+FACR*RHG(I))/RHOGRD
            ENDDO
#if defined (MOD_DNF)
c        write (*,*) 'debug: '
c        write (*,'(A3,E20.12)') '       b3(i,j,k)',b3(1),b3(2),b3(3)
c        write (*,'(A3,E20.12)') '      rhg(i,j,k)',rhg(1),rhg(2),rhg(3)
c        write (*,'(A3,E20.12)') '       vt(i,j,k)',vt(1),vt(2),vt(3)
c        write (*,'(A2,E20.12)') ' sxc,b0,br,f1,f2',sxc,b0,br,
c    &                              VXC(FRR)*B0,VXC(FRZ)*BR
#endif /* MOD_DNF */
C           Linear transformation
C
            CALL DOEXC2(SXC,VT,IREPOP,EXCMAT(1,1,1,IM),GAO,GAO1)
CTROND            CALL DOLRG(VT,IREPOP,EXCMAT(1,1,1,IM),GAO,GAO1)
          END IF
        END IF
      ENDIF
      ENDDO
      RETURN
      END
      subroutine dft_qr(dtrmat,excmat,isymop,ihrmop,gao,gao1,
     &                  rho,nrmgrd,grdrho,wght,dogga,buf)
!***********************************************************************
!
!     Calculate the DFT contributions the quadratic response.
!
!     Input:
!     dtrmat(*,*,*,1) --- B density matrix
!     dtrmat(*,*,*,2) --- C density matrix
!     dtrmat(*,*,*,3) --- BC density matrix
!     isymop          --- operator symmetry
!     ihrmop          --- operator hermiticity
!     gao             --- 
!     gao1            --- 
!     rho             --- \rho in current grid point
!     nrmgrd          --- |\nabla\rho| in current grid point
!                         (corresponds to rhogrd in LR)
!     grdrho          --- \nabla\rho in current grid point
!                         (corresponds to rhg in LR)
!     wght            --- 
!     dogga           --- LDA or GGA
!     buf             --- scratch memory
!
!     Output:
!     excmat(*,*,*,1) --- G^{[B]C} matrix
!     excmat(*,*,*,2) --- G^{[C]B} matrix
!     excmat(*,*,*,3) --- H^{[BC]} matrix
!
!     Written by johhe Jan 2007
!
!***********************************************************************
!
      implicit none
!
#include "mxcent.h"
#include "inforb.h"
#if defined (PRG_DIRAC)
#include "dgroup.h"
#else
#include "wrkrsp.h"
#include "dftinf.h"
#endif
#include "dftcom.h"
#include "functionals.h"
!
!     I/O VARIABLES
!
      double precision dtrmat(nbast,nbast,nz,3)
      double precision excmat(nbast,nbast,nz,3)
      integer isymop(3), ihrmop(3)
      double precision gao(nbast), gao1(nbast,3)
      double precision rho, nrmgrd, grdrho(3), wght, buf(*)
      logical dogga
!
!     INTERNAL VARIABLES
!
!     Symmetry of the different operators
      integer irepbop, irepcop, irepbcop
!
!     One and two times modified densities w.r.t. perturbations B and C
!     (rho_B corresponds to b0 in LR)
      double precision rho_B, rho_C, rho_BC
!
!     Gradients of the modified densities
!     (grdrho_B corresponds to b3 in LR)
      double precision grdrho_B(3), grdrho_C(3), grdrho_BC(3)
!
!     \zeta^X=|\nabla\rho^X|^2 for the different modified densities
!     (zeta_B corresponds to br in LR)
      double precision zeta_B, zeta_C, zeta_BC
!
!     xc-potential (functional derivatives of e_{xc} w.r.t. \rho and \zeta)
      double precision vxc(dervs3)
!
!     s and v in the equation
!       f_{pq} = s \Omega_{pq} + \mathbf{v} \cdot \nabla\Omega_{pq}
      double precision sxc, vt(3)
!
!     Pre-factors associated with different parts of \mathbf{v}
      double precision facr, facb, facc
!
!     To check if a point needs to be included or not.
      double precision bmax, cmax, bcmax
!
!     Loop index
      integer i
!
!     Calculate the ireps for the different perturbing operators
      irepbop  = isymop(1)-1
      irepcop  = isymop(2)-1
      irepbcop = isymop(3)-1
!
!     Get perturbed densities
!     =======================
!
!     Modified density w.r.t. B
!!    if(ihrmop(1).eq.1) then
         call getrho(rho_B,irepbop,gao,dtrmat(1,1,1,1),buf)
!!    else
!!       call quit('DFT-QR general case not implemented!')
!!    end if
!
!     Modified density w.r.t. C
!!    if(ihrmop(2).eq.1) then
         call getrho(rho_C,irepcop,gao,dtrmat(1,1,1,2),buf)
!!    else
!!       call quit('DFT-QR general case not implemented!')
!!    end if
!
!     Two times modified density w.r.t. B and C
!!    if(ihrmop(3).eq.1) then
         call getrho(rho_BC,irepbcop,gao,dtrmat(1,1,1,3),buf)
!!    else
!!       call quit('DFT-QR general case not implemented!')
!!    end if
!
!       ***************
!       ***** LDA *****
!       ***************
!
      if(.not.dogga) then
!
!     Get the functional derivatives
         call dftpot2(vxc,wght,rho,0D0)
!
!     Constructing g^{[C]B} 
         if(dabs(rho_B).gt.dfthrl) then 
!
!           Calculate sxc
            sxc = vxc(frr)*rho_B
!
!           Linear transformation
            call doexc1(sxc,irepbop,excmat(1,1,1,1),gao)
         end if
!
!     Constructing g^{[B]C}
         if(dabs(rho_C).gt.dfthrl) then 
!
!           Calculate sxc
            sxc = vxc(frr)*rho_C
!
!           Linear transformation
            call doexc1(sxc,irepcop,excmat(1,1,1,2),gao)
         END IF
!
!     Constructing h^{[BC]}
         if(dabs(rho_BC).gt.dfthrl) then 
!
!           Calculate sxc
            sxc = 2.0D0*vxc(frrr)*rho_B*rho_C + vxc(frr)*rho_BC
!
!           Linear transformation
            call doexc1(sxc,irepbcop,excmat(1,1,1,3),gao)
         end if
!
      else  
!
!         ***************
!         ***** GGA *****
!         ***************
!
!     Get gradients of perturbed densities
!     ====================================
!
!     Gradients of the different modified densities
      if(ihrmop(1) .eq. 1) then
        call getgrho(grdrho_b,irepbop,
     &               gao,gao1,
     &               dtrmat(1,1,1,1),buf)
      else
#if defined (MOD_CDFT)
        CALL GET_GR_NH(.TRUE.,
     &                 grdrho_b,
     &                 irepbop,
     &                 GAO,
     &                 GAO1,
     &                 DTRMAT(1,1,1,1),
     &                 BUF)

#else
        call quit('DFT-QR general case not implemented!')
#endif
      endif
      if(ihrmop(2) .eq. 1) then
        call getgrho(grdrho_c,irepcop,
     &               gao,gao1,
     &               dtrmat(1,1,1,2),buf)
      else
#if defined (MOD_CDFT)
        CALL GET_GR_NH(.TRUE.,
     &                 grdrho_c,
     &                 irepcop,
     &                 GAO,
     &                 GAO1,
     &                 DTRMAT(1,1,1,2),
     &                 BUF)

#else
        call quit('DFT-QR general case not implemented!')
#endif
      endif
      if(ihrmop(3) .eq. 1) then
        call getgrho(grdrho_bc,irepbcop,
     &               gao,gao1,
     &               dtrmat(1,1,1,3),buf)
      else
#if defined (MOD_CDFT)
        CALL GET_GR_NH(.TRUE.,
     &                 grdrho_bc,
     &                 irepbcop,
     &                 GAO,
     &                 GAO1,
     &                 DTRMAT(1,1,1,3),
     &                 BUF)

#else
        call quit('DFT-QR general case not implemented!')
#endif
      endif
!
!     \zeta^B = \frac{\partial\zeta}{\partial\epsilon_B}
         zeta_B = 2.0D0*(grdrho_B(1)*grdrho(1) + grdrho_B(2)*grdrho(2) +
     &                   grdrho_B(3)*grdrho(3))
!
!     \zeta^B = \frac{\partial\zeta}{\partial\epsilon_C}
         zeta_C = 2.0D0*(grdrho_C(1)*grdrho(1) + grdrho_C(2)*grdrho(2) +
     &                   grdrho_C(3)*grdrho(3))
!
!     \zeta^{BC} = \frac{\partial^2\zeta}{\partial\epsilon_B\partial\epsilon_C}
!     Note the factor 4 which comes from that \rho^{BC} is already symmetrized,
!     and, hence, a factor 2 is enough on the first line.
         zeta_BC = 2.0D0*(grdrho_BC(1)*grdrho(1) + 
     &                    grdrho_BC(2)*grdrho(2) +
     &                    grdrho_BC(3)*grdrho(3)) +
     &             4.0d0*(grdrho_B(1)*grdrho_C(1) +
     &                    grdrho_B(2)*grdrho_C(2) +
     &                    grdrho_B(3)*grdrho_C(3))
!
!     Get the functional derivatives
         call dftpot2(vxc,wght,rho,nrmgrd)
!
!     Constructing g^{[B]C}
!     
         bmax = max(abs(rho_B),
     &              abs(grdrho_B(1)),abs(grdrho_B(2)),abs(grdrho_B(3)))
         if (bmax.gt.dfthrl) then
!
!           Calculate sxc
            sxc  = vxc(frr)*rho_B + vxc(frz)*zeta_B
!     
!           Calculate vt
            facr = vxc(frz)*rho_B + vxc(fzz)*zeta_B
            do i=1,3
               vt(i) = 2.0D0 * (facr*grdrho(i) + vxc(fz0)*grdrho_B(i))
            end do
!
!           Linear transformation
            call doexc2(sxc,vt,irepbop,excmat(1,1,1,1),gao,gao1)
!
         end if
!
!     Constructing g^{[C]B}
!     
         cmax = max(abs(rho_C),
     &              abs(grdrho_C(1)),abs(grdrho_C(2)),abs(grdrho_C(3)))
         if (cmax.gt.dfthrl) then
!
!           Calculate sxc
            sxc  = vxc(frr)*rho_C + vxc(frz)*zeta_C
!     
!           Calculate vt
            facr = vxc(frz)*rho_C + vxc(fzz)*zeta_C
            do I=1,3
               vt(i) = 2.0D0 * (facr*grdrho(i) + vxc(fz0)*grdrho_C(i))
            end do
!
!           Linear transformation
            call doexc2(sxc,vt,irepcop,excmat(1,1,1,2),gao,gao1)
!
         end if
!
!     Constructing h^{[BC]}
!     
         bcmax = max(bmax,cmax,abs(rho_BC),
     &            abs(grdrho_BC(1)),abs(grdrho_BC(2)),abs(grdrho_BC(3)))
         if (bcmax.gt.dfthrl) then
!
!           Calculate sxc
            sxc =  2.0D0*vxc(frrr)*rho_B*rho_C +
     &             2.0d0*vxc(frrz)*(rho_B*zeta_C+rho_C*zeta_B) +
     &             2.0D0*vxc(frzz)*zeta_B*zeta_C+vxc(frr)*rho_BC + 
     &             vxc(frz)*zeta_BC
!     
!           Calculate vt
            facr = 2.0D0*vxc(fzzz)*zeta_B*zeta_C +
     &             2.0d0*vxc(frzz)*(rho_B*zeta_C+rho_C*zeta_B) +
     &             2.0D0*vxc(frrz)*rho_B*rho_C + vxc(frz)*rho_BC +
     &             vxc(fzz)*zeta_BC
            facb = 2.0d0*vxc(frz)*rho_C + 2.0d0*vxc(fzz)*zeta_C
            facc = 2.0d0*vxc(frz)*rho_B + 2.0d0*vxc(fzz)*zeta_B
            do i=1,3
               vt(i) = 2.0D0 * ( facr*grdrho(i) + 
     &                           facb*grdrho_B(i) + facc*grdrho_C(i) +
     &                           vxc(fz0) * grdrho_BC(i) )
            end do
!
!           Linear transformation
            call doexc2(sxc,vt,irepbcop,excmat(1,1,1,3),gao,gao1)
!
         end if
!
      end if
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dolrg */
      SUBROUTINE DOLRG(VT,IREPOP,EXCMAT,GAO,GAO1)
C***********************************************************************
C
C     Construct Kohn-Sham matrix
C
C     Written by T.Saue and T.Helgaker and P.Salek
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
C
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAO(NBAST),GAO1(NBAST,3),EXCMAT(NBAST,NBAST)
#else
      DIMENSION GAO(NTBAS(0)),GAO1(NTBAS(0),3),EXCMAT(NTBAS(0),NTBAS(0))
#endif
      DIMENSION VT(0:3)
C
#include "ibtfun.h"
C
      DO 10 IRP1 = 0,NBSYM-1
        DO 20 IRP2 = 0,NBSYM-1
          IF(IBTXOR(IRP1,IRP2).NE.IREPOP) GOTO 20
          DO 30 IC = 1,MC
            NB1 = NBBAS(IRP1,IC)
            NB2 = NBBAS(IRP2,IC)
            JB1 = IBBAS(IRP1,IC)
            JB2 = IBBAS(IRP2,IC)
            DO J = 1, NB2
               G0 = GAO(JB2+J)
               GX = GAO1(JB2+J,1)
               GY = GAO1(JB2+J,2)
               GZ = GAO1(JB2+J,3)
               DO I = 1, NB1
                  A0 = G0*GAO(JB1+I)
                  AX = GX*GAO(JB1+I) + G0*GAO1(JB1+I,1)
                  AY = GY*GAO(JB1+I) + G0*GAO1(JB1+I,2)
                  AZ = GZ*GAO(JB1+I) + G0*GAO1(JB1+I,3)
                  EXCMAT(JB1+I,JB2+J) = EXCMAT(JB1+I,JB2+J)
     &                 + VT(0)*A0 + VT(1)*AX
     &                 + VT(2)*AY + VT(3)*AZ
              END DO
            ENDDO
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
C
      END
c
      SUBROUTINE FORTWRT(str, slength)
#include "priunit.h"
      INTEGER unit, slength
      CHARACTER str(slength)
      WRITE(LUPRI, *) str
      END

      SUBROUTINE DFTSETCAM(W,BETA)
#include "implicit.h"
#include "dcbham.h"
      HFXMU = W
      HFXATT = BETA 
      END
C /* Deck dft_df_lr */
#if defined (MOD_DNF) 
      SUBROUTINE fit_DFT_LR(pCFIT,pVXCint,NXCMAT,ISYMOP,IHRMOP,
     &                      GAO,GAO1,RHO,RHOGRD,RHG,WGHT,DOGGA,B0,
     &                      nskipped)
C
C adapted from DFT_LR for density fitting
C andre gomes and radovan bast, january 2007 
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0,DP5 = 0.5D0,D2=2.0D0)
C
#include "inforb.h"
#include "nuclei.h"
#if defined (PRG_DIRAC)
#include "dgroup.h"
#else
#include "wrkrsp.h"
#include "dftinf.h"
#endif
#include "dftcom.h"
#include "functionals.h"
C
      INTEGER A, B
      LOGICAL DOGGA, TRPLET, LBIT
      DIMENSION pCFIT(*),pVXCint(*),ISYMOP(NXCMAT),IHRMOP(NXCMAT),
     &          GAO(*), GAO1(*),
     &          RHG(3),VT(3), VXC(DERVS2)
      DIMENSION B3(3)
C
      if (NXCMAT.GT.1) 
     &   CALL QUIT ('DF-VXC not working for multiple densities yet...')

      DO IM = 1,NXCMAT
      IREPOP = ISYMOP(IM)-1
      IF(IHRMOP(IM).EQ.1) THEN
C
C       Get perturbed density
C       =====================
C
        B0 = DDOT(NBASISAUX,GAO,1,pCFIT,1)
C
C       ***************
C       ***** LDA *****
C       ***************
C
        IF(.NOT.DOGGA) THEN
          IF(DABS(B0).GT.DFTHRL) THEN
C
C            Calculated VT
C
             CALL DFTPOT1(VXC, WGHT,RHO, 0D0)
             SXC = VXC(FRR)*B0
C
C            Linear transformation, saving the integrals
C
             CALL fit_DOEXC(DOGGA,SXC,DUMMY,GAO,DUMMY,
     &                      pVXCint,NBASISAUX)
          else
             nskipped=nskipped+1
          END IF
        ELSE
C
C         ***************
C         ***** GGA *****
C         ***************
C
C  getting gradient from perturbed fitted density
C
           DO IXYZ = 1,3
              IOFF = NBASISAUX * (IXYZ - 1) + 1
              b3(ixyz) = 
     &              D1*DDOT(NBASISAUX,GAO1(IOFF),1,pCFIT,1)
           ENDDO
           BMAX = MAX(ABS(B0),ABS(B3(1)),ABS(B3(2)),ABS(B3(3)))
C
           IF (BMAX.GT.DFTHRL) THEN
              CALL DFTPOT1(VXC,WGHT,RHO,RHOGRD)
              BR = (B3(1)*RHG(1)+ B3(2)*RHG(2)+ B3(3)*RHG(3))/RHOGRD
              FACR = VXC(FRZ)*B0 + (VXC(FZZ)-VXC(FZ0)/RHOGRD)*BR
              SXC = VXC(FRR)*B0 + VXC(FRZ)*BR
              DO I = 1,3
                 VT(I) = (VXC(FZ0)*B3(I)+FACR*RHG(I))/RHOGRD
              ENDDO
c        write (*,*) 'debug: '
c        write (*,'(A3,E20.12)') '       b3(i,j,k)',b3(1),b3(2),b3(3)
c        write (*,'(A3,E20.12)') '      rhg(i,j,k)',rhg(1),rhg(2),rhg(3)
c        write (*,'(A3,E20.12)') '       vt(i,j,k)',vt(1),vt(2),vt(3)
c        write (*,'(A2,E20.12)') ' sxc,b0,br,f1,f2',sxc,b0,br,
c    &                              VXC(FRR)*B0,VXC(FRZ)*BR
C
C           Linear transformation, saving the integrals
C
              CALL fit_DOEXC(DOGGA,SXC,VT,GAO,GAO1,pVXCint,NBASISAUX)
          else
             nskipped=nskipped+1
          END IF
        END IF

      ENDIF
      ENDDO
      RETURN
      END
#endif /* MOD_DNF */
