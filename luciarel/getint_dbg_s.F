      SUBROUTINE GETINT_DBG_S(XINT,ITP,ISM,IBU,JTP,JSM,JBU,
     &                        KTP,KSM,KBU,LTP,LSM,LBU,
     &                        IXCHNG,IKSM,JLSM,ICOUL,
     &                        INTLST,IJKLOF,IREOINDX,SCALE)
*
* Obtain integrals
*
*     ICOUL = 0 :
*                  XINT(IK,JL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 1 :
*                  XINT(IJ,KL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
* Storing for ICOUL = 1 not working if IKSM or JLSM .ne. 0
*
* Version for integrals stored in INTLST
*
* Double group integral version. At the moment no
* permutation symmetry
*
* INDXREO gives mapping between input and output indeces
* output  index => input  index
*
* Input index : Index in complete integral list
* Output index : Index in outblock block
*
* all symmetries etc refers to data for the output block
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "lucinp.inc"
*. Reordering of indeces
      INTEGER IREOINDX(4)
*. INDXREO(I) gives the
*. Integral list
      REAL*8 INTLST(*)
      Dimension IJKLof(NsmOB,NsmOb,NsmOB)
*.Output
      DIMENSION XINT(*)
*. A bit of local scratch
      INTEGER INSM(4),IOUTSM(4),ISCR(4)
*
      NTEST = 000
*
      if (NTEST.ge.10) then
         WRITE(6,*) ' Godmorgen dette er GETINT_DBG_S '
         WRITE(6,*) ' IREOINDX : ', (IREOINDX(I),I=1,4)
         WRITE(6,*) ' ISM, JSM, KSM, LSM',ISM, JSM, KSM, LSM
         WRITE(6,*) ' IBU, JBU, KBU, LBU',IBU, JBU, KBU, LBU
         WRITE(6,*) ' SCALE = ', SCALE
      end if
*
      IF(IBU.EQ.2) THEN
         IORB = NOBPTS(ITP,ISM)
         ITORB = NTOOBS(ISM)
      ELSE
         IORB = NOBPTS2(ITP,ISM)
         ITORB = NTOOBS2(ISM)
      END IF
*
      IF(JBU.EQ.2) THEN
         JORB = NOBPTS(JTP,JSM)
         JTORB = NTOOBS(JSM)
      ELSE
         JORB = NOBPTS2(JTP,JSM)
         JTORB = NTOOBS2(JSM)
      END IF
*
      IF(KBU.EQ.2) THEN
         KORB = NOBPTS(KTP,KSM)
         KTORB = NTOOBS(KSM)
      ELSE
         KORB = NOBPTS2(KTP,KSM)
         KTORB = NTOOBS2(KSM)
      END IF
*
      IF(LBU.EQ.2) THEN
         LORB = NOBPTS(LTP,LSM)
         LTORB = NTOOBS(LSM)
      ELSE
         LORB = NOBPTS2(LTP,LSM)
         LTORB = NTOOBS2(LSM)
      END IF
*
      IF(NTEST.GE.100) THEN
         WRITE(6,*)
     &   ' GETINT_DGB_S: IORB,JORB,KORB,LORB',IORB,JORB,KORB,LORB
         WRITE(6,*)
     &   '  ITORB, JTORB, KTORB, LTORB' ,
     &        ITORB, JTORB, KTORB, LTORB
      END IF
*. Offsets relative to start of all orbitals, symmetry ordered
      IF(IBU.EQ.2) THEN
         IOFF = IBSO(ISM)
         IOFF0 = IOFF
         DO IITP = 1, ITP -1
            IOFF = IOFF + NOBPTS(IITP,ISM)
         END DO
      ELSE
         IOFF = IBSO2(ISM)
         IOFF0 = IOFF
         DO IITP = 1, ITP -1
            IOFF = IOFF + NOBPTS2(IITP,ISM)
         END DO
      END IF
*
      IF(JBU.EQ.2) THEN
         JOFF = IBSO(JSM)
         JOFF0 = JOFF
         DO JJTP = 1, JTP -1
            JOFF = JOFF + NOBPTS(JJTP,JSM)
         END DO
      ELSE
         JOFF = IBSO2(JSM)
         JOFF0 = JOFF
         DO JJTP = 1, JTP -1
            JOFF = JOFF + NOBPTS2(JJTP,JSM)
         END DO
      END IF
*
      IF(KBU.EQ.2) THEN
         KOFF = IBSO(KSM)
         KOFF0 = KOFF
         DO KKTP = 1, KTP -1
            KOFF = KOFF + NOBPTS(KKTP,KSM)
         END DO
      ELSE
         KOFF = IBSO2(KSM)
         KOFF0 = KOFF
         DO KKTP = 1, KTP -1
            KOFF = KOFF + NOBPTS2(KKTP,KSM)
         END DO
      END IF
*
      IF(LBU.EQ.2) THEN
         LOFF = IBSO(LSM)
         LOFF0 = LOFF
         DO LLTP = 1, LTP -1
            LOFF = LOFF + NOBPTS(LLTP,LSM)
         END DO
      ELSE
         LOFF = IBSO2(LSM)
         LOFF0 = LOFF
         DO LLTP = 1, LTP -1
            LOFF = LOFF + NOBPTS2(LLTP,LSM)
         END DO
      END IF
*
      IF(NTEST.GE.100) THEN
         WRITE(6,*)
     &  ' GETINT_DBG_S: IOFF,JOFF,KOFF,LOFF',
     &                  IOFF,JOFF,KOFF,LOFF
         WRITE(6,*)
     &  ' GETINT_DBG_S: IOFF0,JOFF0,KOFF0,LOFF0',
     &                  IOFF0,JOFF0,KOFF0,LOFF0
      END IF
*
*. Symmetry of input block
      ISCR(IREOINDX(1)) = ISM
      ISCR(IREOINDX(2)) = JSM
      ISCR(IREOINDX(3)) = KSM
      ISCR(IREOINDX(4)) = LSM
*
      ISMI = ISCR(1)
      JSMI = ISCR(2)
      KSMI = ISCR(3)
      LSMI = ISCR(4)
*
      if (NTEST.ge.100) then
        WRITE(6,*) ' ISMI = ', ISMI
        WRITE(6,*) ' JSMI = ', JSMI
        WRITE(6,*) ' KSMI = ', KSMI
        WRITE(6,*) ' LSMI = ', LSMI
      end if
*
*. Input dimensions
      ISCR(IREOINDX(1)) = ITORB
      ISCR(IREOINDX(2)) = JTORB
      ISCR(IREOINDX(3)) = KTORB
      ISCR(IREOINDX(4)) = LTORB
*
      NI_IN = ISCR(1)
      NJ_IN = ISCR(2)
      NK_IN = ISCR(3)
      NL_IN = ISCR(4)
*
      if (NTEST.ge.100) then
        WRITE(6,*) ' NI_IN = ', NI_IN
        WRITE(6,*) ' NJ_IN = ', NJ_IN
        WRITE(6,*) ' NK_IN = ', NK_IN
        WRITE(6,*) ' NL_IN = ', NL_IN
      end if
*
*
*     Collect Coulomb terms
*
*
      if (NTEST.ge.100) WRITE(6,*) ' Before IBLOFF '
*
      IBLOFF = IJKLOF(ISMI,JSMI,KSMI)
*
      if (NTEST.ge.100) then
         WRITE(6,*) ' IBLOFF = ', IBLOFF
         WRITE(6,*) ' Block of two-electron integrals in input lists '
         CALL WRTMAT(INTLST(IBLOFF),NI_IN*NJ_IN,NK_IN*NL_IN,
     &               NI_IN*NJ_IN,NK_IN*NL_IN )
      end if
*
      IINT = 0
*. Loop over indeces of output blocks
      Do L=LOFF,LOFF+LORB-1
         if (NTEST.ge.100) write(6,*) 'Doing L= ',L
         JMIN = JOFF
         If (JLSM.ne.0) JMIN = L
         Do J=JMIN,JOFF+JORB-1
            if (NTEST.ge.100) write(6,*) 'Doing J= ',J
*
            Do K=KOFF,KOFF+KORB-1
               if (NTEST.ge.100) write(6,*) 'Doing K= ',K
               IMIN = IOFF
               If (IKSM.ne.0) IMIN = K
               IF (ICOUL.EQ.1) THEN
*. Address before integral (1,j!k,l)
                IINT = (L-LOFF)*JORB*KORB*IORB + (K-KOFF)*JORB*IORB
     &                  + (J-JOFF)*IORB
               END IF
*
*. Block (ISM JSM ! KSM LSM ) with (Ism,jsm) > (ksm,lsm)
               Do I=IMIN,IOFF+IORB-1
                  if (NTEST.ge.100) write(6,*) 'Doing I= ',I
*. 4 output indeces are now known, corresponding input indeces
                  ISCR(IREOINDX(1)) = I - IOFF0 + 1
                  ISCR(IREOINDX(2)) = J - JOFF0 + 1
                  ISCR(IREOINDX(3)) = K - KOFF0 + 1
                  ISCR(IREOINDX(4)) = L - LOFF0 + 1
                  I_IN = ISCR(1)
                  J_IN = ISCR(2)
                  K_IN = ISCR(3)
                  L_IN = ISCR(4)
*. (I J K L ) in input list
                  IJKL_IN = (L_IN-1)*NK_IN*NJ_IN*NI_IN
     &                 +    (K_IN-1)      *NJ_IN*NI_IN
     &                 +    (J_IN-1)            *NI_IN
     &                 +     I_IN
*
                  IJKL = IBLOFF-1+ IJKL_IN
                  IINT = IINT + 1
*
                  if (NTEST.ge.100) then
                     WRITE(6,*) ' IJKL_IN' , IJKL_IN
                     WRITE(6,*) ' IINT and IJKL ',IINT, IJKL
                  end if
*
                  XINT(IINT) = SCALE*INTLST(IJKL)
*
                  if (NTEST.ge.100) then
                     write(6,*) 'Processed INTLST(',IJKL,')
     &                          to XINT(',IINT,')'
                  end if
               End Do
*
            End Do
         End Do
      End Do
*
*     Collect Exchange terms
*
      If ( IXCHNG.ne.0 ) Then
*
*
         IBLOFF = IJKLOF(ISMI,LSMI,KSMI)
*
         IINT = 0
         Do L=LOFF,LOFF+LORB-1
            JMIN = JOFF
            If (JLSM.ne.0) JMIN = L
            Do J=JMIN,JOFF+JORB-1
               Do K=KOFF,KOFF+KORB-1
*
                  IF(ICOUL.EQ.1)  THEN
*. Address before integral (1,j!k,l)
                   IINT = (L-LOFF)*JORB*KORB*IORB + (K-KOFF)*JORB*IORB
     &                     + (J-JOFF)*IORB
                  END IF
*
                  If(IKSM.ne.0) IMIN = K
                  Do I=IMIN,IOFF+IORB-1
*. 4 output indeces are now known, corresponding input indeces
                     ISCR(IREOINDX(1)) = I - IOFF0 + 1
                     ISCR(IREOINDX(2)) = J - JOFF0 + 1
                     ISCR(IREOINDX(3)) = K - KOFF0 + 1
                     ISCR(IREOINDX(4)) = L - LOFF0 + 1
                     I_IN = ISCR(1)
                     J_IN = ISCR(2)
                     K_IN = ISCR(3)
                     L_IN = ISCR(4)
*. (I L !K J) in output list
                     ILKJ0 = (J_IN-1)*NK_IN*NL_IN*NI_IN
     &                     + (K_IN-1)      *NL_IN*NI_IN
     &                     + (L_IN-1)            *NI_IN
     &                     +  I_IN
                     ILKJ = ILKJ0 + IBLOFF - 1
                     IINT = IINT + 1
                     XINT(IINT)=XINT(IINT)-SCALE*INTLST(ILKJ)
                  End Do
*
               End Do
            End Do
         End Do
      End If
*
      IF(NTEST.GE.100) THEN
        IF(ICOUL.EQ.1) THEN
          WRITE(6,*) ' output block in output format '
          CALL WRTMAT(XINT,IORB*JORB,KORB*LORB,IORB*JORB,KORB*LORB)
        END IF
      END IF
*
      Return
      End
