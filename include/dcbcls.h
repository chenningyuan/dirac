C
C     CLASS LABELS
C
C     CLSINT - keyword for a class of integrals
C     CLSCMB - indicates what blocks to compute : LL,SL,LS,SS
C     CLSCAL - indicates whether class is active, that is to be calculated
C
      PARAMETER(MAXCLS = 15)
      LOGICAL CLSCAL
      CHARACTER CLSINT*7,CLSCMB*4
      COMMON/XCBCLS/NPRPCLS,IORDCL(MAXCLS),CLSCAL(MAXCLS)
      COMMON/CLSLBL/CLSINT(MAXCLS),CLSCMB(MAXCLS)
