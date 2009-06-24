!
!MI... some description would be usefull....
!
! radovan: i agree Miro ...
!          i am guessing because there is hardly any documentation:
!          NHKT is l
!          in Dalton KHKT would be the degeneracy among spherical Gaussians
!          KHKT = 2*(l - 1) + 1
!          and KCKT the corresponding cartesian degeneracy
!          KCKT = l*(l + 1)/2
!          but in DIRAC with cartesian basis both give the same (cartesian) degeneracy
!          JSTRT is an offset for the shell
!
      LOGICAL SHARE(MXSHEL), SEGM(MXSHEL), SPHR(MXSHEL)
      REAL*8 CENT(MXSHEL, 3, 2)
      INTEGER NHKT(MXSHEL), KHKT(MXSHEL), KCKT(MXSHEL), ISTBAO(MXSHEL), &
     &        NUCO(MXSHEL), JSTRT(MXSHEL), NSTRT(MXSHEL), NCENT(MXSHEL),&
     &        NRCO(MXSHEL), NUMCF(MXSHEL), NBCH(MXSHEL), KSTRT(MXSHEL), &
     &        LCLASS(MXSHEL), IPTSHL(MXSHEL), NUMCFT(MXSHEL), KMAX,     &
     &        NLRGSH, NSMLSH, NORBS, KMAXTOT,                           &
     &        MST(MXSHEL), MBSID(MXSHEL), MBIDBT(MXSHEL),               &
     &        NAUXSH, KMAXAUX, NORBAUX, IFSTSHAUX
      COMMON /SHELLS/ CENT, NHKT, KHKT, KCKT, ISTBAO, NUCO, JSTRT,      &
     &                NSTRT, NCENT, SHARE, NRCO, NUMCF, NBCH, KSTRT,    &
     &                SEGM, LCLASS, IPTSHL, NUMCFT, SPHR, KMAX, NLRGSH, &
     &                NSMLSH, NORBS, KMAXTOT,                           &
     &                MST, MBSID, MBIDBT, NAUXSH, KMAXAUX, NORBAUX,     &
     &                IFSTSHAUX
!     MBSID and MBIDBT have been added for multiple basis sets (WK/UniKA/31-10-2002).
