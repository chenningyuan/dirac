!
!MI... some description would be usefull....
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
