C
C     PROPERTIES:
C
C       IPRPTYP  - indicates kind of property operator
C         1. P             * scalar operator
C         2. [alpha_x]p    * x-component of alpha times scalar operator
C         3. [alpha_y]p    * y-component of alpha times scalar operator
C         4. [alpha_z]p    * z-component of alpha times scalar operator
C         5. [Alpha x P]_x * vector product of alpha and vector operator,
C                            x-component    
C         6. [Alpha x P]_y * vector product of alpha and vector operator,
C                            y-component    
C         7. [Alpha x P]_z * vector product of alpha and vector operator,
C                            z-component    
C         8. A.P           * dot-product of alpha and vector operator
C
C       IPRPLBL - pointers to property labels
C       IPRPSYM - boson symmetry of operator  (1 - 8)
C       IPRPTIM - symmetry under time reversal (1,2)
C       PRPNAM  - user specified name of operator
C
      PARAMETER ( MAXPRPS = 200 )
      CHARACTER*16  PRPNAM
      COMMON/XCBPRP/FACPRP(3,MAXPRPS),
     &              NPRPS,IPRPTYP(MAXPRPS),IPRPLBL(3,MAXPRPS),
     &              IPRPSYM(MAXPRPS),IPRPTIM(MAXPRPS)
      COMMON/XCBPRPC/PRPNAM(MAXPRPS)
