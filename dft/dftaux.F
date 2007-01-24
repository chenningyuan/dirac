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
      SUBROUTINE DFT_LINRES(DTRMAT,EXCMAT,NXCMAT,ISYMOP,IHRMOP,GAO,GAO1,
     &                      RHO,RHOGRD,RHG,WGHT,DOGGA,BUF)
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
#include "cdft.h"
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

C           Linear transformation
C
            CALL DOEXC2(SXC,VT,IREPOP,EXCMAT(1,1,1,IM),GAO,GAO1)
CTROND            CALL DOLRG(VT,IREPOP,EXCMAT(1,1,1,IM),GAO,GAO1)
          END IF
        END IF
#if defined (MOD_CDFT)
      ELSEIF(IHRMOP(IM) .EQ. -1) THEN

        IF(LBIT(LBIT_CDFT,1)) THEN !spin density contribution to linear response
          CALL SDFT_LR(DTRMAT(1,1,1,IM),
     &                 EXCMAT(1,1,1,IM),
     &                 IREPOP,
     &                 GAO,GAO1,
     &                 RHO,RHG,RHOGRD,
     &                 WGHT,DOGGA,BUF)
        ENDIF

        IF(LBIT(LBIT_CDFT,2)) THEN !current density contribution to linear response
          CALL CDFT_LR(DTRMAT(1,1,1,IM),
     &                 EXCMAT(1,1,1,IM),
     &                 IREPOP,
     &                 GAO,GAO1,
     &                 RHO,RHG,RHOGRD,
     &                 WGHT,DOGGA,BUF)
        ENDIF

      ELSE
        CALL QUIT('DTRMAT has to be hermitian or antihermitian')
#endif
      ENDIF
      ENDDO
      RETURN
      END
      subroutine dft_quadres(dtrmat,excmat,isymop,ihrmop,gao,gao1,
     &                       rho,rhogrd,rhg,wght,dogga,buf)
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
!     rhogrd          --- |\nabla\rho| in current grid point
!     rhg             --- \nabla\rho in current grid point
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
      integer irepbop, irepcop, irepbcop
      double precision gao(nbast), gao1(nbast,3)
      double precision rho, rhogrd,  rhg(3), wght, buf(*)
      logical dogga
!
!     INTERNAL VARIABLES
!
      double precision b3(3), c3(3), bc3(3)
      double precision b0, br, c0, cr, bc0, bcr
      double precision sxc, vt(3)
      double precision vxc(dervs3)
      double precision facr, facb, facc, bmax, cmax, bcmax
      integer i
!
!     Calculate the ireps for the different perturbing operators
      irepbop = isymop(1)-1
      irepcop = isymop(2)-1
      irepbcop = isymop(3)-1
!
!     Get perturbed densities
!     =======================
!
!     B-op density
      if(ihrmop(1).eq.1) then
         call getrho(b0,irepbop,gao,dtrmat(1,1,1,1),buf)
      else
         call quit('DFT-QR general case not implemented!')
      end if
!
!     C-op density
      if(ihrmop(2).eq.1) then
         call getrho(c0,irepcop,gao,dtrmat(1,1,1,2),buf)
      else
         call quit('DFT-QR general case not implemented!')
      end if
!
!     BC-op density
      if(ihrmop(3).eq.1) then
         call getrho(bc0,irepbcop,gao,dtrmat(1,1,1,3),buf)
      else
         call quit('DFT-QR general case not implemented!')
      end if
!
!       ***************
!       ***** LDA *****
!       ***************
!
      if(.not.dogga) then
!
!     Possibly check sizes of B0, C0, and BC0, if all small, skip this.
!     Get the functional derivatives
         call dftpot2(vxc,wght,rho,0D0)
!
!     Constructing g^{[C]B} 
         if(dabs(b0).gt.dfthrl) then 
!
!           Calculate sxc
            sxc = vxc(frr)*b0
!
!           Linear transformation
            call doexc1(sxc,irepbop,excmat(1,1,1,1),gao)
         end if
!
!     Constructing g^{[B]C}
         if(dabs(c0).gt.dfthrl) then 
!
!           Calculate sxc
            sxc = vxc(frr)*c0
!
!           Linear transformation
            call doexc1(sxc,irepcop,excmat(1,1,1,2),gao)
         END IF
!
!     Constructing h^{[BC]}
         if(dabs(bc0).gt.dfthrl) then 
!
!           Calculate sxc
            sxc = 2.0D0*vxc(frrr)*b0*c0 + vxc(frr)*bc0
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
!     B-, C, and BC-op gradient densityíes
         call getgrho(b3, irepbop, gao,gao1,dtrmat(1,1,1,1),buf)
         call getgrho(c3, irepcop, gao,gao1,dtrmat(1,1,1,2),buf)
         call getgrho(bc3,irepbcop,gao,gao1,dtrmat(1,1,1,3),buf)
!
!     \frac{\partial\zeta}{\partial\epsilon_B}
         br = 2.0D0*(b3(1)*rhg(1)+ b3(2)*rhg(2)+ b3(3)*rhg(3))
!
!     \frac{\partial\zeta}{\partial\epsilon_C}
         cr = 2.0D0*(c3(1)*rhg(1)+ c3(2)*rhg(2)+ c3(3)*rhg(3))
!
!     \frac{\partial^2\zeta}{\partial\epsilon_B\partial\epsilon_C}
!     Not the extra factor 2 on the second line!!!
         bcr = 2.0D0*(bc3(1)*rhg(1)+ bc3(2)*rhg(2)+ bc3(3)*rhg(3)) +
     &         4.0d0*(b3(1)*c3(1) + b3(2)*c3(2) + b3(3)*c3(3))
!
!     Possibly check sizes of B0, C0, and BC0, if all small, skip this.
   !     Get the functional derivatives
         call dftpot2(vxc,wght,rho,rhogrd)
!
!     Constructing g^{[B]C}
!     
         bmax = max(abs(b0),abs(b3(1)),abs(b3(2)),abs(b3(3)))
         if (bmax.gt.dfthrl) then
!
!           Calculate sxc
            sxc  = vxc(frr)*b0 + vxc(frz)*br
!     
!           Calculate vt
            facr = vxc(frz)*b0 + vxc(fzz)*br
            do i=1,3
               vt(i) = 2.0D0 * ( facr*rhg(i) + vxc(fz0)*b3(i) )
            end do
!
!           Linear transformation
            call doexc2(sxc,vt,irepbop,excmat(1,1,1,1),gao,gao1)
!
         end if
!
!     Constructing g^{[C]B}
!     
         cmax = max(abs(c0),abs(c3(1)),abs(c3(2)),abs(c3(3)))
         if (cmax.gt.dfthrl) then
!
!           Calculate sxc
            sxc  = vxc(frr)*c0 + vxc(frz)*cr
!     
!           Calculate vt
            facr = vxc(frz)*c0 + vxc(fzz)*cr
            do I=1,3
               vt(i) = 2.0D0 * ( facr*rhg(i) + vxc(fz0)*c3(i) )
            end do
!
!           Linear transformation
            call doexc2(sxc,vt,irepcop,excmat(1,1,1,2),gao,gao1)
!
         end if
!
!     Constructing h^{[BC]}
!     
         bcmax = max(abs(bc0),abs(bc3(1)),abs(bc3(2)),abs(bc3(3)),
     &               bmax,cmax)
         if (bcmax.gt.dfthrl) then
!
!           Calculate sxc
            sxc =  2.0D0*vxc(frrr)*b0*c0+2.0d0*vxc(frrz)*(b0*cr+c0*br)+
     &             2.0D0*vxc(frzz)*br*cr+vxc(frr)*bc0 + 
     &             vxc(frz)*bcr
!     
!           Calculate vt
            facr = 2.0D0*vxc(fzzz)*br*cr+2.0d0*vxc(frzz)*(b0*cr+c0*br) +
     &             2.0D0*vxc(frrz)*b0*c0 + vxc(frz)*bc0 +
     &             vxc(fzz)*bcr
            facb = 2.0d0*vxc(frz)*c0 + 2.0d0*vxc(fzz)*cr
            facc = 2.0d0*vxc(frz)*b0 + 2.0d0*vxc(fzz)*br
            do i=1,3
               vt(i) = 2.0D0 * ( facr*rhg(i) + facb*b3(i) + facc*c3(i) +
     &                           vxc(fz0) * bc3(i) )
            end do
!
!           Linear transformation
            call doexc2(sxc,vt,irepbcop,excmat(1,1,1,3),gao,gao1)
!
         end if
!
      end if
      return
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
