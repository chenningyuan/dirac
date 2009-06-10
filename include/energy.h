      REAL*8          ENERKE,                                           &
     &                ENEREE,                                           &
     &                ENERNA,                                           &
     &                ENERNN,                                           &
#ifndef PRG_DIRAC
!radovan: Hessian not implemented in DIRAC
!         if you implement it, please allocate dynamically
!         and don't put the matrices in here
     &                HESFS1(MXCOOR, MXCOOR),                           &
     &                HESFS2(MXCOOR, MXCOOR),                           &
     &                HESREL(MXCOOR, MXCOOR),                           &
     &                HESRL2(MXCOOR, MXCOOR),                           &
     &                HESRL3(MXCOOR, MXCOOR),                           &
     &                HESSEE(MXCOOR, MXCOOR),                           &
     &                HESSKE(MXCOOR, MXCOOR),                           &
     &                HESSNA(MXCOOR, MXCOOR),                           &
     &                HESSNN(MXCOOR, MXCOOR),                           &
#endif /* ifndef PRG_DIRAC */
     &                GRADEE(MXCOOR),                                   &
     &                GRADFS(MXCOOR),                                   &
     &                GRADKE(MXCOOR),                                   &
     &                GRADNA(MXCOOR),                                   &
     &                GRADNN(MXCOOR)
      COMMON /ENERGY/ ENERKE,                                           &
     &                ENEREE,                                           &
     &                ENERNA,                                           &
     &                ENERNN,                                           &
#ifndef PRG_DIRAC
!radovan: Hessian not implemented in DIRAC
!         if you implement it, please allocate dynamically
!         and don't put the matrices in here
     &                HESFS1,                                           &
     &                HESFS2,                                           &
     &                HESREL,                                           &
     &                HESRL2,                                           &
     &                HESRL3,                                           &
     &                HESSEE,                                           &
     &                HESSKE,                                           &
     &                HESSNA,                                           &
     &                HESSNN,                                           &
#endif /* ifndef PRG_DIRAC */
     &                GRADEE,                                           &
     &                GRADFS,                                           &
     &                GRADKE,                                           &
     &                GRADNA,                                           &
     &                GRADNN
