!
!     PROPERTIES:
!
!       IPRPTYP - indicates kind of property operator
!         0. FOCKMAT         * Fock-matrix
!         1. P               * diagonal times scalar operator
!         2. i[alpha_x]P     * x-component of alpha times scalar operator
!         3. i[alpha_y]P     * y-component of alpha times scalar operator
!         4. i[alpha_z]P     * z-component of alpha times scalar operator
!         5. i[alpha x P]_x  * vector product of alpha and vector operator,
!                              x-component
!         6. i[alpha x P]_y  * vector product of alpha and vector operator,
!                              y-component
!         7. i[alpha x P]_z  * vector product of alpha and vector operator,
!                              z-component
!         8. i[alpha.P]      * dot-product of alpha and vector operator
!         9. [gamma5]P       * gamma5 times scalar operator
!        10. i[Sigma_x]P     * x-component of Sigma times scalar operator
!        11. i[Sigma_y]P     * y-component of Sigma times scalar operator
!        12. i[Sigma_z]P     * z-component of Sigma times scalar operator
!        13. i[betaSigma_x]P * x-component of beta Sigma times scalar operator
!        14. i[betaSigma_y]P * y-component of beta Sigma times scalar operator
!        15. i[betaSigma_z]P * z-component of beta Sigma times scalar operator
!        16. i[betaalpha_x]P * x-component of beta alpha times scalar operator
!        17. i[betaalpha_y]P * y-component of beta alpha times scalar operator
!        18. i[betaalpha_z]P * z-component of beta alpha times scalar operator
!        19. [beta]P         * beta times scalar operator
!        20. i[Sigma.P]      * dot-product of Sigma and vector operator
!        21. [betagamma5]P   * beta gamma5 times scalar operator
!        22. i[beta alpha.P] * dot-product of beta alpha and vector operator
!        23. i[betagamma5]P  * i.beta gamma5 times scalar operator
!
!       IPRPLBL - pointers to property labels
!       IPRPSYM - boson symmetry of operator  (1 - 8)
!       IPRPTIM - symmetry under time reversal (+1,-1)
!       PRPNAM  - user specified name of operator
!
      INTEGER MAXPRPS
      PARAMETER (MAXPRPS = 5328)

      CHARACTER*16 PRPNAM(MAXPRPS)
      REAL*8 FACPRP(3, MAXPRPS)
      INTEGER NPRPS, IPRPTYP(MAXPRPS), IPRPLBL(3, MAXPRPS),             &
     &        IPRPSYM(MAXPRPS), IPRPTIM(MAXPRPS)
      COMMON /XCBPRP/ FACPRP, NPRPS, IPRPTYP, IPRPLBL, IPRPSYM, IPRPTIM

      COMMON /XCBPRPC/ PRPNAM
