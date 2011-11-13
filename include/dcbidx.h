! FILE    : include/dcbidx.h
! VERSION : $Revision: 7317 $
! DATE    : $Date: 2008-08-07 16:47:12 +0200 (Thu, 07 Aug 2008) $
!
! Needs: maxorb.h, maxash.h
!
! PURPOSE: index arrays for the krmc module (and prp module when krmc)
!
! Originally designed by HJAaJ 2001 based on similar data structures
! in SIRIUS in DALTON.
!
!     ISMO  : Fermion symmetry of orbital (1 for gerade, 2 for ungerade)
!     IOBTYP: type of orbital:
!        - frozen orbital         : 0
!        - positronic orbital     : 1
!        - inactive orbital       : 2
!        - active orbital         : 3
!        - secondary orbital      : 4
!
!     IDXa2b index vector gives index of "a" to "b" renumbering
!     where "a" and "b" are one of:
!       G : general orbital
!       U : active orbital
!       E : electronic orbital ("positive energy" orbital)
!       P : positronic orbital ("negative energy" orbital)
!       I : inactive orbital (doubly occupied orbital)
!       F : frozen inactive orbital
!       S : secondary orbital (empty electronic orbital)
!       V : active secondary orbital (non-deleted secondary orbitals)
!       D : deleted empty orbitals
!     IDXT2G(:,1:4) :
!       general orbital index of each of the four "transformation" orbitals
!       for the four indices in the integral transformation (1 2 | 3 4)
!     NIDX*, IIDX* :
!       info related to the four indices in the integral transformation
!
!     
      INTEGER JTFRO, JTPOSI, JTINAC, JTACT, JTSEC
      PARAMETER (JTFRO=0, JTPOSI=1, JTINAC=2, JTACT=3, JTSEC=4)

      INTEGER IFSMO(MAXORB) , IOBTYP(MAXORB), IDXG2U(MAXORB),           &
     &        IDXU2G(MAXASH), IDXE2G(MAXORB), IDXG2E(MAXORB),           &
     &        IDXP2G(MAXORB), IDXG2P(MAXORB), IDXI2G(MAXORB),           &
     &        IDXG2I(MAXORB), IDXS2G(MAXORB), IDXG2S(MAXORB),           &
     &        IDXV2G(MAXORB), IDXG2V(MAXORB), IDXG2D(MAXORB),           &
     &        IDXD2G(MAXORB), IDXT2G(MAXORB, 4),                        &
     &        NIDXT(4), IIDX1(2), NIDX1(2), IIDX2(2), NIDX2(2),         &
     &                  IIDX3(2), NIDX3(2), IIDX4(2), NIDX4(2)
      COMMON /DCBIDX/ IFSMO,  IOBTYP, IDXG2U, IDXU2G, IDXE2G, IDXG2E,   &
     &                IDXP2G, IDXG2P, IDXI2G, IDXG2I, IDXS2G, IDXG2S,   &
     &                IDXV2G, IDXG2V, IDXG2D, IDXD2G, IDXT2G,           &
     &                NIDXT, IIDX1, NIDX1, IIDX2, NIDX2,                &
     &                       IIDX3, NIDX3, IIDX4, NIDX4
! end of include/dcbidx.h
