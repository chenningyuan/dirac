C VERSION : $Revision$
C DATE    : $Date$
C SOURCE  : $Source$
C
C     needs: maxorb.h, maxash.h
C
C     ISMO  : Fermion symmetry of orbital
C     IOBTYP: type of orbital:
C        - frozen orbital    : 0
C        - positronic orbital: 1
c        - inactive orbital  : 2
c        - active orbital    : 3
c        - secondary orbital : 4
C     IDXG2U: index vector: general orbital to active orbital
C     IDXU2G: index vector: active orbital to general orbital
C
C     IDXE2U: index vector: electron orbital to active orbital
C     IDXU2E: index vector: active orbital to electronic orbital
C     
      PARAMETER (
     &     JTFRO = 0,
     &     JTPOSI = 1,
     &     JTINAC = 2,
     &     JTACT = 3,
     &     JTSEC = 4)
C
      PARAMETER (
     &     JDHF = 0,
     &     JODHF = 1,
     &     JCAS = 2,
     &     JGAS = 3)
C
      COMMON /DCBIDX/
     &     IFSMO(MAXORB),
     &     IOBTYP(MAXORB),
     &     IDXG2U(MAXORB),
     &     IDXU2G(MAXASH),
     &     IDXE2G(MAXORB),
     &     IDXG2E(MAXORB),
     &     IDXP2G(MAXORB),
     &     IDXG2P(MAXORB),
     &     IDXI2G(MAXORB),
     &     IDXG2I(MAXORB),
     &     IDXV2G(MAXORB),
     &     IDXG2V(MAXORB),
     &     IDXT2G(MAXORB,4), NIDXT(4),
     &     IIDX1(2), NIDX1(2), 
     &     IIDX2(2), NIDX2(2),
     &     IIDX3(2), NIDX3(2), 
     &     IIDX4(2), NIDX4(2)
