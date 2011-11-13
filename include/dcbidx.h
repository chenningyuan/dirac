! VERSION : $Revision: 1.5 $
! DATE    : $Date: 2004/05/26 13:06:54 $
!
!     needs: maxorb.h, maxash.h
!
!     ISMO  : Fermion symmetry of orbital
!     IOBTYP: type of orbital:
!        - frozen orbital    : 0
!        - positronic orbital: 1
!        - inactive orbital  : 2
!        - active orbital    : 3
!        - secondary orbital : 4
!     IDXG2U: index vector: general orbital to active orbital
!     IDXU2G: index vector: active orbital to general orbital
!
!     IDXE2U: index vector: electron orbital to active orbital
!     IDXU2E: index vector: active orbital to electronic orbital
!     
      INTEGER JTFRO, JTPOSI, JTINAC, JTACT, JTSEC
      PARAMETER (JTFRO = 0, JTPOSI = 1, JTINAC = 2, JTACT = 3, JTSEC =  &
     &           4)

      INTEGER IFSMO(MAXORB), IOBTYP(MAXORB), IDXG2U(MAXORB),            &
     &        IDXU2G(MAXASH), IDXE2G(MAXORB), IDXG2E(MAXORB),           &
     &        IDXP2G(MAXORB), IDXG2P(MAXORB), IDXI2G(MAXORB),           &
     &        IDXG2I(MAXORB), IDXV2G(MAXORB), IDXG2V(MAXORB),           &
     &        IDXT2G(MAXORB, 4), NIDXT(4), IIDX1(2), NIDX1(2), IIDX2(2),&
     &        NIDX2(2), IIDX3(2), NIDX3(2), IIDX4(2), NIDX4(2)
      COMMON /DCBIDX/ IFSMO, IOBTYP, IDXG2U, IDXU2G, IDXE2G, IDXG2E,    &
     &                IDXP2G, IDXG2P, IDXI2G, IDXG2I, IDXV2G, IDXG2V,   &
     &                IDXT2G, NIDXT, IIDX1, NIDX1, IIDX2, NIDX2, IIDX3, &
     &                NIDX3, IIDX4, NIDX4
