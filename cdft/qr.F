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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      INTEGER FUNCTION I_DERIVATIVE(I)
C=======================================================================
C     describe me ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                               last revision: june 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

          IF(I .EQ.     1) THEN
            I_DERIVATIVE =  1
      ELSEIF(I .EQ.   100) THEN
            I_DERIVATIVE =  2
      ELSEIF(I .EQ. 10000) THEN
            I_DERIVATIVE =  3
      ELSEIF(I .EQ.     2) THEN
            I_DERIVATIVE =  4
      ELSEIF(I .EQ.    20) THEN
            I_DERIVATIVE =  5
      ELSEIF(I .EQ.   200) THEN
            I_DERIVATIVE =  6
      ELSEIF(I .EQ.  2000) THEN
            I_DERIVATIVE =  7
      ELSEIF(I .EQ.   101) THEN
            I_DERIVATIVE =  8
      ELSEIF(I .EQ.  1010) THEN
            I_DERIVATIVE =  9
      ELSEIF(I .EQ. 10001) THEN
            I_DERIVATIVE = 10
      ELSEIF(I .EQ. 10100) THEN
            I_DERIVATIVE = 11
      ELSEIF(I .EQ.     3) THEN
            I_DERIVATIVE = 12
      ELSEIF(I .EQ.   300) THEN
            I_DERIVATIVE = 13
      ELSEIF(I .EQ.    21) THEN
            I_DERIVATIVE = 14
      ELSEIF(I .EQ.   201) THEN
            I_DERIVATIVE = 15
      ELSEIF(I .EQ.  2001) THEN
            I_DERIVATIVE = 16
      ELSEIF(I .EQ.   102) THEN
            I_DERIVATIVE = 17
      ELSEIF(I .EQ.   120) THEN
            I_DERIVATIVE = 18
      ELSEIF(I .EQ.  1011) THEN
            I_DERIVATIVE = 19
      ELSEIF(I .EQ.  1110) THEN
            I_DERIVATIVE = 20
      ELSEIF(I .EQ.  2100) THEN
            I_DERIVATIVE = 21
      ELSE
            I_DERIVATIVE =  99
      ENDIF

      RETURN
      END
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     /* qr_c */
      subroutine qr_c(ivar,
     &                jvar,
     &                kvar,
     &                der,
     &                ixyz_start,
     &                var,
     &                var_pq,
     &                use_var)
c=======================================================================
c     describe me ...
c-----------------------------------------------------------------------
c                                                           radovan bast
c                                               last revision: june 2007
c=======================================================================

#include "implicit.h"
#include "priunit.h"

      parameter (d0 = 0.0d0,d1 = 1.0d0,d2 = 2.0d0,d3 = 3.0d0,d4 = 4.0d0)

      logical   use_var(5)
      dimension var   (3,0:3,5),
     &          var_pq(3,0:3,5)

      if(use_var(ivar) .and. use_var(jvar) .and. use_var(kvar)) then
       do ixyz = ixyz_start,3
        var_pq(ixyz,0,ivar) = var_pq(ixyz,0,ivar)
     &                      + var(ixyz,1,jvar)*var(ixyz,2,kvar)*der
     &                      + var(ixyz,2,jvar)*var(ixyz,1,kvar)*der
       enddo
      endif

      return
      end
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     /* qr_b */
      subroutine qr_b(ivar,
     &                jvar,
     &                der,
     &                ixyz_start,
     &                var,
     &                var_pq,
     &                use_var)
c=======================================================================
c     describe me ...
c-----------------------------------------------------------------------
c                                                           radovan bast
c                                               last revision: june 2007
c=======================================================================

#include "implicit.h"
#include "priunit.h"

      parameter (d0 = 0.0d0,d1 = 1.0d0,d2 = 2.0d0,d3 = 3.0d0,d4 = 4.0d0)

      logical   use_var(5)
      dimension var   (3,0:3,5),
     &          var_pq(3,0:3,5)

      if(use_var(ivar) .and. use_var(jvar)) then
       do ixyz = ixyz_start,3
        var_pq(ixyz,0,ivar) = var_pq(ixyz,0,ivar) + var(ixyz,3,jvar)*der
        var_pq(ixyz,1,ivar) = var_pq(ixyz,1,ivar) + var(ixyz,2,jvar)*der
        var_pq(ixyz,2,ivar) = var_pq(ixyz,2,ivar) + var(ixyz,1,jvar)*der
       enddo
      endif

      return
      end
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     /* qr_bx */
      subroutine qr_bx(ivar,
     &                 jvar,
     &                 der,
     &                 ixyz_start,
     &                 var,
     &                 var_pq,
     &                 use_var)
c=======================================================================
c     describe me ...
c-----------------------------------------------------------------------
c                                                           radovan bast
c                                               last revision: june 2007
c=======================================================================

#include "implicit.h"
#include "priunit.h"

      parameter (d0 = 0.0d0,d1 = 1.0d0,d2 = 2.0d0,d3 = 3.0d0,d4 = 4.0d0)

      logical   use_var(5)
      dimension var   (3,0:3,5),
     &          var_pq(3,0:3,5)

      if(use_var(ivar) .and. use_var(jvar)) then
       do ixyz = ixyz_start,3
        var_pq(ixyz,1,ivar) = var_pq(ixyz,1,ivar) + var(ixyz,2,jvar)*der
        var_pq(ixyz,2,ivar) = var_pq(ixyz,2,ivar) + var(ixyz,1,jvar)*der
       enddo
      endif

      return
      end
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     /* dft_response */
      subroutine dft_response(add_r,
     &                        add_s,
     &                        add_gr,
     &                        add_gs,
     &                        dtrmat,
     &                        excmat,
     &                        isymop,
     &                        ihrmop,
     &                        gao,
     &                        gao1,
     &                        r_0,
     &                        gr_0_norm,
     &                        gr_0,
     &                        w,
     &                        buf)
c=======================================================================
c     describe me ...
c-----------------------------------------------------------------------
c                                                           radovan bast
c                                               last revision: june 2007
c=======================================================================

#include "implicit.h"
#include "priunit.h"

      parameter (d0 = 0.0d0,d1 = 1.0d0,d2 = 2.0d0,d3 = 3.0d0,d4 = 4.0d0)
      parameter (nr_vars = 5)
!     parameter (nr_vars = 6)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dcbham.h"
#include "dgroup.h"
#include "inforb.h"
#include "cdft.h"
#include "dftcom.h"
#include "../dft/functionals.h"

      logical   add_r,
     &          add_s,
     &          add_gr,
     &          add_gs,
     &          lbit,
     &          use_var(nr_vars),
     &          var_is_product(nr_vars)
      dimension dtrmat(nbast,nbast,nz,3),
     &          excmat(nbast,nbast,nz,3),
     &          isymop(3),
     &          ihrmop(3),
     &          gao(*),
     &          gao1(*),
     &          gr_0(3),
     &          buf(*),
!
     &          gr(  3,0:3),
     &          gs(3,3,0:3),
!
     &          var   (3,0:3,nr_vars),
     &          var_pq(3,0:3,nr_vars),
!
     &          vxc_vector(0:3,3,0:3,2),
!
     &          vxc(99)

#include "ibtfun.h"

      do ivar = 1,nr_vars
        use_var(ivar)        = .false.
        var_is_product(ivar) = .false.
      enddo

      if(add_r)  use_var(1) = .true.
      if(add_s)  use_var(2) = .true.
      if(add_gr) use_var(3) = .true.
      if(add_gs) use_var(4) = .true.
      if(add_gs) use_var(5) = .true.
!                use_var(6) = .true.

      var_is_product(3) = .true.
      var_is_product(4) = .true.
      var_is_product(5) = .true.

      call dzero(gr,4*3)
      call dzero(gs,4*3*3)

      call dzero(var,   nr_vars*4*4)
      call dzero(var_pq,nr_vars*4*4)

      call dzero(vxc_vector,4*3*4*2)

      ixyz_start = 1
      if(lbit(lbit_cdft,3)) ixyz_start = 3 !collinear approximation

      do ibc = 1,3
        if(use_var(1)) then
          call get_r(.true.,
     &               var(3,ibc,1),
     &               isymop(ibc)-1,
     &               gao,
     &               dtrmat(1,1,1,ibc),
     &               buf)
        endif
        if(use_var(2)) then
          call get_m(.true.,
     &               var(1,ibc,2),
     &               isymop(ibc)-1,
     &               gao,
     &               dtrmat(1,1,1,ibc),
     &               buf)
        endif
        if(use_var(3)) then
          call get_gr_nh(.true.,
     &                   gr(1,ibc),
     &                   isymop(ibc)-1,
     &                   gao,
     &                   gao1,
     &                   dtrmat(1,1,1,ibc),
     &                   buf)
        endif
        if(use_var(4)) then
          call get_m_dervs_nh(.true.,
     &                       gs(1,1,ibc),
     &                       isymop(ibc)-1,
     &                       gao,
     &                       gao1,
     &                       dtrmat(1,1,1,ibc),
     &                       buf)
        endif
      enddo

      gr(1,0) = gr_0(1)
      gr(2,0) = gr_0(2)
      gr(3,0) = gr_0(3)

!     Z = \nabla R \cdot \nabla R
      ivar = 3
      if(use_var(ivar)) then
        ixyz = 3
          var(ixyz,1,ivar) = ddot(3,gr(1,1),1,gr(1,0),1)*d2
          var(ixyz,2,ivar) = ddot(3,gr(1,2),1,gr(1,0),1)*d2
          var(ixyz,3,ivar) = ddot(3,gr(1,3),1,gr(1,0),1)*d2
     &                     + ddot(3,gr(1,1),1,gr(1,2),1)*d4
      endif

!     Y = \nabla R \cdot \nabla S
      ivar = 4
      if(use_var(ivar)) then
        do ixyz = ixyz_start,3
          var(ixyz,1,ivar) = ddot(3,gs(1,ixyz,1),1,gr(1,0),1)
          var(ixyz,2,ivar) = ddot(3,gs(1,ixyz,2),1,gr(1,0),1)
          var(ixyz,3,ivar) = ddot(3,gs(1,ixyz,3),1,gr(1,0),1)
     &                     + ddot(3,gs(1,ixyz,1),1,gr(1,2),1)*d2
     &                     + ddot(3,gs(1,ixyz,2),1,gr(1,1),1)*d2
        enddo
      endif

!     X = \nabla S \cdot \nabla S
      ivar = 5
      if(use_var(ivar)) then
        do ixyz = ixyz_start,3
          var(ixyz,1,ivar) = d0
          var(ixyz,2,ivar) = d0
          var(ixyz,3,ivar) = ddot(3,gs(1,ixyz,1),1,gs(1,ixyz,2),1)*d4
        enddo
      endif

!     M = "S" \cdot "S"
!     ivar = 6
!     if(use_var(ivar)) then
!       do ixyz = ixyz_start,3
!         var(ixyz,1,ivar) = d0
!         var(ixyz,2,ivar) = d0
!         var(ixyz,3,ivar) = var(ixyz,1,2)*var(ixyz,2,2)*d4
!       enddo
!     endif

      if(add_r .and.       add_gr) call  dftpot2(vxc,w,r_0,gr_0_norm)
      if(add_r .and. .not. add_gr) call  dftpot2(vxc,w,r_0,d0)
      if(add_s .and.       add_gs) call sdftpot2(vxc,w,r_0,gr_0_norm)
      if(add_s .and. .not. add_gs) call sdftpot2(vxc,w,r_0,d0)

#if defined brute_force_it
!     if i_derivative() doesn't find a non-zero derivative
!     it will return the pointer 99
      vxc(99) = d0

      do ixyz = ixyz_start,3

!       this way i can use the loop over x, y, z
!       also for purely collinear variables (R and Z)
        do ibc = 1,3
          var(ixyz,ibc,1) = var(3,ibc,1)
          var(ixyz,ibc,3) = var(3,ibc,3)
        enddo

        do ivar = 1,nr_vars
          if(use_var(ivar)) then
            ichecksum =             10**(ivar-1)
            der = vxc(i_derivative(ichecksum))

            var_pq(ixyz,3,ivar) = der

            do jvar = 1,nr_vars
              if(use_var(jvar)) then
                jchecksum = ichecksum + 10**(jvar-1)
                der = vxc(i_derivative(jchecksum))

                var_pq(ixyz,2,ivar) = var_pq(ixyz,2,ivar)
     &                              + var(ixyz,1,jvar)*der
                var_pq(ixyz,1,ivar) = var_pq(ixyz,1,ivar)
     &                              + var(ixyz,2,jvar)*der
                var_pq(ixyz,0,ivar) = var_pq(ixyz,0,ivar)
     &                              + var(ixyz,3,jvar)*der

                do kvar = 1,nr_vars
                  if(use_var(kvar)) then
                    kchecksum = jchecksum + 10**(kvar-1)
                    der = vxc(i_derivative(kchecksum))

                    var_pq(ixyz,0,ivar) = var_pq(ixyz,0,ivar)
     &                           + var(ixyz,1,jvar)*var(ixyz,2,kvar)*der
     &                           + var(ixyz,2,jvar)*var(ixyz,1,kvar)*der

                  endif
                enddo
              endif
            enddo
          endif
        enddo
      enddo
#endif

      do ixyz = ixyz_start,3

!       this way i can use the loop over x, y, z
!       also for purely collinear variables (R and Z)
        do ibc = 1,3
          var(ixyz,ibc,1) = var(3,ibc,1)
          var(ixyz,ibc,3) = var(3,ibc,3)
        enddo
      enddo

!     a(Z)
      ivar = 3
      ixyz = 3
        var_pq(ixyz,3,ivar) = vxc(fZ0)

!     a(X)
      ivar = 5
      do ixyz = ixyz_start,3
        var_pq(ixyz,3,ivar) = vxc(fX)
      enddo

!     b(R,R)
      call qr_b(1,1,vxc(fRR),3,var,var_pq,use_var)

!     b(S,S)
      call qr_b(2,2,vxc(fSS),ixyz_start,var,var_pq,use_var)

!     a(M)
!     ivar = 6
!     do ixyz = ixyz_start,3
!       var_pq(ixyz,3,ivar) = vxc(fSS)/d2
!     enddo

!     b(Z,Z)
      call qr_b(3,3,vxc(fZZ),3,var,var_pq,use_var)

!     b(Y,Y)
      call qr_b(4,4,vxc(fYY),ixyz_start,var,var_pq,use_var)

!     b(R,Z)
      call qr_b(1,3,vxc(fRZ),3,var,var_pq,use_var)
      call qr_b(3,1,vxc(fRZ),3,var,var_pq,use_var)

!     b(S,Y)
      call qr_b(2,4,vxc(fSY),ixyz_start,var,var_pq,use_var)
      call qr_b(4,2,vxc(fSY),ixyz_start,var,var_pq,use_var)

!     b'(X,R)
      call qr_bx(5,1,vxc(fRX),ixyz_start,var,var_pq,use_var)

!     b'(X,Z)
      call qr_bx(5,3,vxc(fZX),ixyz_start,var,var_pq,use_var)

!     c(R,R,R)
      call qr_c(1,1,1,vxc(fRRR),3,var,var_pq,use_var)

!     c(Z,Z,Z)
      call qr_c(3,3,3,vxc(fZZZ),3,var,var_pq,use_var)

!     c(R,R,Z)
      call qr_c(1,1,3,vxc(fRRZ),3,var,var_pq,use_var)
      call qr_c(1,3,1,vxc(fRRZ),3,var,var_pq,use_var)
      call qr_c(3,1,1,vxc(fRRZ),3,var,var_pq,use_var)

!     c(R,S,S)
      call qr_c(1,2,2,vxc(fRSS),ixyz_start,var,var_pq,use_var)
      call qr_c(2,1,2,vxc(fRSS),ixyz_start,var,var_pq,use_var)
      call qr_c(2,2,1,vxc(fRSS),ixyz_start,var,var_pq,use_var)

!     b(n,M)
!     call qr_b(1,6,vxc(fRSS)/d2,ixyz_start,var,var_pq,use_var)
!     call qr_b(6,1,vxc(fRSS)/d2,ixyz_start,var,var_pq,use_var)

!     c(S,S,Z)
      call qr_c(2,2,3,vxc(fSSZ),ixyz_start,var,var_pq,use_var)
      call qr_c(2,3,2,vxc(fSSZ),ixyz_start,var,var_pq,use_var)
      call qr_c(3,2,2,vxc(fSSZ),ixyz_start,var,var_pq,use_var)

!     c(R,Z,Z)
      call qr_c(1,3,3,vxc(fRZZ),3,var,var_pq,use_var)
      call qr_c(3,1,3,vxc(fRZZ),3,var,var_pq,use_var)
      call qr_c(3,3,1,vxc(fRZZ),3,var,var_pq,use_var)

!     c(R,Y,Y)
      call qr_c(1,4,4,vxc(fRYY),ixyz_start,var,var_pq,use_var)
      call qr_c(4,1,4,vxc(fRYY),ixyz_start,var,var_pq,use_var)
      call qr_c(4,4,1,vxc(fRYY),ixyz_start,var,var_pq,use_var)

!     c(Z,Y,Y)
      call qr_c(3,4,4,vxc(fZYY),ixyz_start,var,var_pq,use_var)
      call qr_c(4,3,4,vxc(fZYY),ixyz_start,var,var_pq,use_var)
      call qr_c(4,4,3,vxc(fZYY),ixyz_start,var,var_pq,use_var)

!     c(R,S,Y)
      call qr_c(1,2,4,vxc(fRSY),ixyz_start,var,var_pq,use_var)
      call qr_c(1,4,2,vxc(fRSY),ixyz_start,var,var_pq,use_var)
      call qr_c(2,1,4,vxc(fRSY),ixyz_start,var,var_pq,use_var)
      call qr_c(2,4,1,vxc(fRSY),ixyz_start,var,var_pq,use_var)
      call qr_c(4,1,2,vxc(fRSY),ixyz_start,var,var_pq,use_var)
      call qr_c(4,2,1,vxc(fRSY),ixyz_start,var,var_pq,use_var)

!     c(S,Z,Y)
      call qr_c(2,3,4,vxc(fSZY),ixyz_start,var,var_pq,use_var)
      call qr_c(2,4,3,vxc(fSZY),ixyz_start,var,var_pq,use_var)
      call qr_c(3,2,4,vxc(fSZY),ixyz_start,var,var_pq,use_var)
      call qr_c(3,4,2,vxc(fSZY),ixyz_start,var,var_pq,use_var)
      call qr_c(4,2,3,vxc(fSZY),ixyz_start,var,var_pq,use_var)
      call qr_c(4,3,2,vxc(fSZY),ixyz_start,var,var_pq,use_var)

!     R
!     contributes to the h+ part
      ivar = 1
      if(use_var(ivar)) then
        ixyz = 3
          jxyz = 0
            do ibc = 1,3
              vxc_vector(jxyz,ixyz,ibc,1) = var_pq(ixyz,3-ibc,ivar)
            enddo
      endif

!     S
!     contributes to the h- part
      ivar = 2
      if(use_var(ivar)) then
        do ixyz = ixyz_start,3
          jxyz = 0
            do ibc = 1,3
              vxc_vector(jxyz,ixyz,ibc,2) = var_pq(ixyz,3-ibc,ivar)
            enddo
        enddo
      endif

!     Z = \nabla R \cdot \nabla R
!     contributes to the h+ part
      ivar = 3
      if(use_var(ivar)) then
        ixyz = 3
          do jxyz = 1,3
            do ibc = 1,3
              vxc_vector(jxyz,ixyz,ibc,1) = vxc_vector(jxyz,ixyz,ibc,1)
     &                         + var_pq(ixyz,3-ibc,ivar)*gr(jxyz,0  )*d2
     &                         + var_pq(ixyz,3    ,ivar)*gr(jxyz,ibc)*d2
            enddo
            ibc = 3
              vxc_vector(jxyz,ixyz,ibc,1) = vxc_vector(jxyz,ixyz,ibc,1)
     &                         + var_pq(ixyz,1    ,ivar)*gr(jxyz,1  )*d4
     &                         + var_pq(ixyz,2    ,ivar)*gr(jxyz,2  )*d4
          enddo
      endif

!     Y = \nabla R \cdot \nabla S
!     contributes to both h+ and h-
      ivar = 4
      if(use_var(ivar)) then
        do ixyz = ixyz_start,3
          do jxyz = 1,3
            do ibc = 1,3
              vxc_vector(jxyz,ixyz,ibc,2) = vxc_vector(jxyz,ixyz,ibc,2)
     &                         + var_pq(ixyz,3-ibc,ivar)*gr(jxyz,0  )
     &                         + var_pq(ixyz,3    ,ivar)*gr(jxyz,ibc)
              vxc_vector(jxyz,ixyz,ibc,1) = vxc_vector(jxyz,ixyz,ibc,1)
     &                        + var_pq(ixyz,3,ivar)*gs(jxyz,ixyz,ibc)
            enddo
            ibc = 3
              vxc_vector(jxyz,ixyz,ibc,2) = vxc_vector(jxyz,ixyz,ibc,2)
     &                         + var_pq(ixyz,1    ,ivar)*gr(jxyz,1  )*d2
     &                         + var_pq(ixyz,2    ,ivar)*gr(jxyz,2  )*d2
              vxc_vector(jxyz,ixyz,ibc,1) = vxc_vector(jxyz,ixyz,ibc,1)
     &                        + var_pq(ixyz,1,ivar)*gs(jxyz,ixyz,1  )*d2
     &                        + var_pq(ixyz,2,ivar)*gs(jxyz,ixyz,2  )*d2
          enddo
        enddo
      endif

!     X = \nabla S \cdot \nabla S
!     contributes to the h- part
      ivar = 5
      if(use_var(ivar)) then
        do ixyz = ixyz_start,3
          do jxyz = 1,3
            do ibc = 1,3
              vxc_vector(jxyz,ixyz,ibc,2) = vxc_vector(jxyz,ixyz,ibc,2)
     &                        + var_pq(ixyz,3,ivar)*gs(jxyz,ixyz,ibc)*d2
            enddo
            ibc = 3
              vxc_vector(jxyz,ixyz,ibc,2) = vxc_vector(jxyz,ixyz,ibc,2)
     &                        + var_pq(ixyz,1,ivar)*gs(jxyz,ixyz,1  )*d4
     &                        + var_pq(ixyz,2,ivar)*gs(jxyz,ixyz,2  )*d4
          enddo
        enddo
      endif

!     M = "S" \cdot "S"
!     contributes to the h- part
!     ivar = 6
!     if(use_var(ivar)) then
!       do ixyz = ixyz_start,3
!         jxyz = 0
!           do ibc = 1,3
!             vxc_vector(jxyz,ixyz,ibc,2) = vxc_vector(jxyz,ixyz,ibc,2)
!    &                        + var_pq(ixyz,3,ivar)*var(ixyz,ibc,2)*d2
!           enddo
!           ibc = 3
!             vxc_vector(jxyz,ixyz,ibc,2) = vxc_vector(jxyz,ixyz,ibc,2)
!    &                        + var_pq(ixyz,1,ivar)*var(ixyz,1,2)*d4
!    &                        + var_pq(ixyz,2,ivar)*var(ixyz,2,2)*d4
!       enddo
!     endif

      do ibc = 1,3
        worth_it = max(
     &                 dabs(vxc_vector(0,3,ibc,1)),
     &                 dabs(vxc_vector(1,3,ibc,1)),
     &                 dabs(vxc_vector(2,3,ibc,1)),
     &                 dabs(vxc_vector(3,3,ibc,1))
     &                )
        if(worth_it .gt. dfthrl) then
          call doexc2(vxc_vector(0,3,ibc,1),
     &                vxc_vector(1,3,ibc,1),
     &                isymop(ibc)-1,
     &                excmat(1,1,1,ibc),
     &                gao,
     &                gao1)
        endif

        do ixyz = ixyz_start,3
          worth_it = max(
     &                   dabs(vxc_vector(0,ixyz,ibc,2)),
     &                   dabs(vxc_vector(1,ixyz,ibc,2)),
     &                   dabs(vxc_vector(2,ixyz,ibc,2)),
     &                   dabs(vxc_vector(3,ixyz,ibc,2))
     &                  )
          if(worth_it .gt. dfthrl) then

!           scaling only the h- part to get symmetrization right
!           the h+ part is scaled inside doexc2
            vxc_vector(1,ixyz,ibc,2) = vxc_vector(1,ixyz,ibc,2)*d2
            vxc_vector(2,ixyz,ibc,2) = vxc_vector(2,ixyz,ibc,2)*d2
            vxc_vector(3,ixyz,ibc,2) = vxc_vector(3,ixyz,ibc,2)*d2

            call exc_m(add_gs,
     &                 vxc_vector(0,ixyz,ibc,2),
     &                 vxc_vector(1,ixyz,ibc,2),
     &                 ixyz,
     &                 isymop(ibc),
     &                 excmat(1,1,1,ibc),
     &                 gao,
     &                 gao1)
          endif
        enddo
      enddo

      return
      end
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     /* dft_qr_slim */
      subroutine dft_qr_slim(add_grad,
     &                       add_spin,
     &                       dtrmat,
     &                       excmat,
     &                       isymop,
     &                       ihrmop,
     &                       gao,
     &                       gao1,
     &                       r_0,
     &                       gr_0_norm,
     &                       gr_0,
     &                       w,
     &                       buf)
c=======================================================================
c     describe me ...
c-----------------------------------------------------------------------
c                                                           radovan bast
c                                               last revision: june 2007
c=======================================================================

#include "implicit.h"
#include "priunit.h"

      parameter (d0 = 0.0d0,d1 = 1.0d0,d2 = 2.0d0,d3 = 3.0d0,d4 = 4.0d0)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dcbham.h"
#include "dgroup.h"
#include "inforb.h"
#include "cdft.h"
#include "dftcom.h"
#include "../dft/functionals.h"

      logical   add_grad,
     &          add_spin,
     &          lbit
      dimension dtrmat(nbast,nbast,nz,3),
     &          excmat(nbast,nbast,nz,3),
     &          isymop(3),
     &          ihrmop(3),
     &          gao(*),
     &          gao1(*),
     &          buf(*),
!
     &          gr_0 (  3),
     &          gr_b (  3),
     &          gr_c (  3),
!
     &          s_b  (  3),
     &          s_c  (  3),
!
     &          gs_b (3,3),
     &          gs_c (3,3),
!
     &          Y_b  (  3),
     &          Y_c  (  3),
     &          Y_bc (  3),
!
     &          X_bc (  3),
!
     &          v_bc_h_Y_0(3),
     &          v_bc_h_Y_b(3),
     &          v_bc_h_Y_c(3),
     &          v_bc_h_X_b(3),
     &          v_bc_h_X_c(3),
!
!               here we construct (see thesis for notation):
!
!               K_{pq,0}^[BC] = r_0^[BC] \Omega_{pq} + q_{0j}^[BC] \nabla \Omega_{pq}
!               K_{pq,z}^[BC] = r_z^[BC] \Omega_{pq} + q_{zj}^[BC] \nabla \Omega_{pq}
!               K_{pq,y}^[BC] = r_y^[BC] \Omega_{pq} + q_{yj}^[BC] \nabla \Omega_{pq}
!               K_{pq,x}^[BC] = r_x^[BC] \Omega_{pq} + q_{xj}^[BC] \nabla \Omega_{pq}
!
!               all the rest is done elsewhere
!
!               scalar h+ part: r_0^[BC]
!               s_bc_h
!               vector h+ part: q_{0j}^[BC] with j = x,y,z
     &          v_bc_h (  3),
!               scalar h- part: r_k^[BC]    with k = x,y,z
     &          s_bc_ah(  3),
!               vector h- part: q_{kj}^[BC] with k = x,y,z and j = x,y,z
     &          v_bc_ah(3,3),
!
     &          vxc(99)

#include "ibtfun.h"

        call getrho(r_b,isymop(1)-1,gao,
     &            dtrmat(1,1,1,1),buf)
        call getrho(r_c,isymop(2)-1,gao,
     &            dtrmat(1,1,1,2),buf)

      if(add_spin) then
        call get_m(.true.,s_b,isymop(1)-1,gao,
     &                  dtrmat(1,1,1,1),buf)
        call get_m(.true.,s_c,isymop(2)-1,gao,
     &                  dtrmat(1,1,1,2),buf)
      endif

      if(add_grad) then
        call get_gr_nh(.true.,gr_b,isymop(1)-1,gao,gao1,
     &                       dtrmat(1,1,1,1),buf)
        call get_gr_nh(.true.,gr_c,isymop(2)-1,gao,gao1,
     &                       dtrmat(1,1,1,2),buf)

!       Z = \nabla R \cdot \nabla R
        Z_b  = d2*ddot(3,gr_0,1,gr_b,1)
        Z_c  = d2*ddot(3,gr_0,1,gr_c,1)
        Z_bc = d2*ddot(3,gr_b,1,gr_c,1)
      endif

      if(add_grad .and. add_spin) then
        call get_m_dervs_nh(.true.,gs_b,isymop(1)-1,gao,gao1,
     &                            dtrmat(1,1,1,1),buf)
        call get_m_dervs_nh(.true.,gs_c,isymop(2)-1,gao,gao1,
     &                            dtrmat(1,1,1,2),buf)

!       Y = \nabla R \cdot \nabla S
        do k = k_start,3
          Y_b (k) = ddot(3,gr_0,1,gs_b(1,k),1)
          Y_c (k) = ddot(3,gr_0,1,gs_c(1,k),1)
          Y_bc(k) = ddot(3,gr_b,1,gs_c(1,k),1)
     &            + ddot(3,gr_c,1,gs_b(1,k),1)
        enddo

!       X = \nabla S \cdot \nabla S
        do k = k_start,3
          X_bc(k) = d2*ddot(3,gs_b(1,k),1,gs_c(1,k),1)
        enddo
      endif

      if(add_grad) then
        call dftpot2(vxc,w,r_0,gr_0_norm)
      else
        call dftpot2(vxc,w,r_0,d0)
      endif

      if(add_spin) then
        if(add_grad) then
          call sdftpot2(vxc,w,r_0,gr_0_norm)
        else
          call sdftpot2(vxc,w,r_0,d0)
        endif
      endif

      k_start = 1
      if(lbit(lbit_cdft,3)) k_start = 3 !collinear approximation

!-------------------------------------------------------------------------------
!     collect s_bc_h

      s_bc_h = d0

      s_bc_h = s_bc_h + vxc(fRRR)*(r_b*r_c)

      if(add_grad) then
        s_bc_h = s_bc_h + vxc(fRRZ)*(r_b*Z_c + Z_b*r_c)
     &                  + vxc(fRZZ)*(Z_b*Z_c)
     &                  + vxc(fRZ)*Z_bc
      endif

      if(add_spin) then
        do k = k_start,3
          s_bc_h = s_bc_h + vxc(fRSS)*(s_b(k)*s_c(k))

          if(add_grad) then
            s_bc_h = s_bc_h + vxc(fRSY)*(s_b(k)*Y_c(k) + Y_b(k)*s_c(k))
     &                      + vxc(fRYY)*(Y_b(k)*Y_c(k))
     &                      + vxc(fRX)*X_bc(k)
          endif
        enddo
      endif

!-------------------------------------------------------------------------------
!     collect s_bc_ah

      call dzero(s_bc_ah,3)

      if(add_spin) then
        do k = k_start,3
          s_bc_ah(k) = s_bc_ah(k) + vxc(fRSS)*(r_b*s_c(k) + s_b(k)*r_c)

          if(add_grad) then
           s_bc_ah(k) = s_bc_ah(k) + vxc(fRSY)*(r_b*Y_c(k) + Y_b(k)*r_c)
     &                             + vxc(fSSZ)*(Z_b*s_c(k) + s_b(k)*Z_c)
     &                             + vxc(fSZY)*(Z_b*Y_c(k) + Y_b(k)*Z_c)
     &                             + vxc(fSY)*Y_bc(k)
          endif
        enddo
      endif

!-------------------------------------------------------------------------------
!     collect vector components

      call dzero(v_bc_h,3)
      call dzero(v_bc_ah,9)

      if(add_grad) then

        v_bc_h_Z_0 = d0
        v_bc_h_Z_b = d0
        v_bc_h_Z_c = d0

        v_bc_h_Z_0 = v_bc_h_Z_0
     &             + vxc(fZZZ)*(Z_b*Z_c)
     &             + vxc(fRZZ)*(r_b*Z_c + Z_b*r_c)
     &             + vxc(fRRZ)*(r_b*r_c)

        v_bc_h_Z_0 = v_bc_h_Z_0
     &             + vxc(fZZ)*Z_bc
        v_bc_h_Z_b = v_bc_h_Z_b
     &             + vxc(fRZ)*r_c
     &             + vxc(fZZ)*Z_c
        v_bc_h_Z_c = v_bc_h_Z_c
     &             + vxc(fRZ)*r_b
     &             + vxc(fZZ)*Z_b

        if(add_spin) then

          call dzero(v_bc_h_Y_0,3)
          call dzero(v_bc_h_Y_b,3)
          call dzero(v_bc_h_Y_c,3)

          call dzero(v_bc_h_X_b,3)
          call dzero(v_bc_h_X_c,3)

          do k = k_start,3
            v_bc_h_Z_0  = v_bc_h_Z_0
     &                  + vxc(fZYY)*(Y_b(k)*Y_c(k))
     &                  + vxc(fSZY)*(s_b(k)*Y_c(k) + Y_b(k)*s_c(k))
     &                  + vxc(fSSZ)*(s_b(k)*s_c(k))

            v_bc_h_Y_0(k) = v_bc_h_Y_0(k)
     &                    + vxc(fYY)*Y_bc(k)
            v_bc_h_Y_b(k) = v_bc_h_Y_b(k)
     &                    + vxc(fSY)*s_c(k)
     &                    + vxc(fYY)*Y_c(k)
            v_bc_h_Y_c(k) = v_bc_h_Y_c(k)
     &                    + vxc(fSY)*s_b(k)
     &                    + vxc(fYY)*Y_b(k)

            v_bc_h_X_b(k) = v_bc_h_X_b(k)
     &                    + vxc(fRX)*r_c
     &                    + vxc(fZX)*Z_c
            v_bc_h_X_c(k) = v_bc_h_X_c(k)
     &                    + vxc(fRX)*r_b
     &                    + vxc(fZX)*Z_b

            v_bc_h_Y_0(k) = v_bc_h_Y_0(k)
     &                    + vxc(fRSY)*(r_b*s_c(k) + s_b(k)*r_c)
     &                    + vxc(fRYY)*(r_b*Y_c(k) + Y_b(k)*r_c)
     &                    + vxc(fSZY)*(Z_b*s_c(k) + s_b(k)*Z_c)
     &                    + vxc(fZYY)*(Z_b*Y_c(k) + Y_b(k)*Z_c)
          enddo
        endif

        do j = 1,3
          v_bc_h(j) = v_bc_h(j) + v_bc_h_Z_0*gr_0(j)*d2
     &                          + v_bc_h_Z_b*gr_b(j)*d2
     &                          + v_bc_h_Z_c*gr_c(j)*d2
          if(add_spin) then
            do k = k_start,3
              v_bc_ah(j,k) = v_bc_ah(j,k) + v_bc_h_Y_0(k)*gr_0(j)
     &                                    + v_bc_h_Y_b(k)*gr_b(j)
     &                                    + v_bc_h_Y_c(k)*gr_c(j)

              v_bc_h(j) = v_bc_h(j) + v_bc_h_Y_b(k)*gs_b(j,k)
     &                              + v_bc_h_Y_c(k)*gs_c(j,k)

              v_bc_ah(j,k) = v_bc_ah(j,k) + v_bc_h_X_b(k)*gs_b(j,k)*d2
     &                                    + v_bc_h_X_c(k)*gs_c(j,k)*d2
            enddo
          endif
        enddo
      endif

!-------------------------------------------------------------------------------

      s_bc_h =     d2*s_bc_h
      call dscal(3,d2,v_bc_h,1)

      worth_it = max(
     &               dabs(s_bc_h),
     &               dabs(v_bc_h(1)),
     &               dabs(v_bc_h(2)),
     &               dabs(v_bc_h(3))
     &              )
      if(worth_it .gt. dfthrl) then
        call doexc2(s_bc_h,
     &              v_bc_h,
     &              isymop(3)-1,
     &              excmat(1,1,1,3),
     &              gao,
     &              gao1)
      endif

      if(.not. add_spin) goto 999

      call dscal(3,d2,s_bc_ah,1)
      call dscal(9,d2,v_bc_ah,1)

!     scaling the vector h- part v_bc_ah to get symmetrization right
!     the h+ part v_bc_h is scaled inside doexc2
      call dscal(9,d2,v_bc_ah,1)

      do k = k_start,3
        worth_it = max(
     &                 dabs(s_bc_ah(k)),
     &                 dabs(v_bc_ah(1,k)),
     &                 dabs(v_bc_ah(2,k)),
     &                 dabs(v_bc_ah(3,k))
     &                )
        if(worth_it .gt. dfthrl) then
          call exc_m(add_grad,
     &               s_bc_ah(k),
     &               v_bc_ah(1,k),
     &               k,
     &               isymop(3),
     &               excmat(1,1,1,3),
     &               gao,
     &               gao1)
        endif
      enddo

 999  continue

      return
      end
