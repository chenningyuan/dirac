C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dirac, a relativistic ab initio electronic structure program, 
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dirtra.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getlab */
      SUBROUTINE GETLAB(WORK,LWORK,IPRINT)
C*****************************************************************************
C
C      Generate primitive labels for basis functions
C
C      Written by T.Saue - November 1994
C      Last revision: Nov 23 1994 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
C
#include <dcblab.h>
#include <symmet.h>
#include <aosotr.h>
#include <dgroup.h>
#include <pgroup.h>
#include <nuclei.h>
#include <shells.h>
#include <dcbbas.h>
#include <ccom.h>
      CHARACTER CLS(2)*1,DEG*2
      DIMENSION WORK(LWORK)
C
#include <dcbibt.h>
      CALL QENTER('GETLAB')
#include <memint.h>
      CLS(1) = 'L'
      CLS(2) = 'S'
C
C     **********************************************
C     ***** L A B E L S  for   A O - B A S I S *****
C     **********************************************
C
C
C     Make labels in the order they appear in HERMIT
C     ==============================================
C
      NLAB  = 0
      NAORB = 0
      DO 10 ISHELL = 1,NLRGSH
        LVAL  = NHKT(ISHELL)
        ICENT = NCENT(ISHELL)
        ICLS  = LCLASS(ISHELL)
        NDEG  = NUCDEG(ICENT)
        DO 20 IDEG = 1,NDEG
          IOFF = NHKOFF(LVAL)
          DO 30 ICOMP = 1,KHK(LVAL)
            ITYP  = IOFF + ICOMP
            NAORB = NAORB + 1
            JLAB  = IPACK(ICLS,ICENT,IDEG,ITYP)
C
C           Check if labels are already defined
C           ==================================
            DO 40 ILAB = 1,NLAB
            IF(JLAB.EQ.IATTR(ILAB,1)) THEN
              IPLAB(NAORB,1) = ILAB
              GOTO 30
            ENDIF
C
 40         CONTINUE
            NLAB          = NLAB + 1
            IF(NLAB.GT.MAXLAB) GOTO 999
            DEG = '  '
            IF(NDEG.GT.1) WRITE(DEG,'(I2)') IDEG
            IATTR(NLAB,1) = JLAB
            PLABEL(NLAB,1) = CLS(ICLS)//' '//NAMN(ICENT)(1:3)//
     &          DEG//' '//GTOTYP(ITYP)
            IPLAB(NAORB,1) = NLAB
 30       CONTINUE
 20     CONTINUE
 10   CONTINUE
      NPLAB(1) = NLAB
      DO 50 ISHELL = NLRGSH+1,KMAX
        LVAL  = NHKT(ISHELL)
        ICENT = NCENT(ISHELL)
        ICLS  = LCLASS(ISHELL)
        NDEG  = NUCDEG(ICENT)
        DO 60 IDEG = 1,NDEG
          IOFF = NHKOFF(LVAL)
          DO 70 ICOMP = 1,KHK(LVAL)
            ITYP  = IOFF + ICOMP
            NAORB = NAORB + 1
            JLAB  = IPACK(ICLS,ICENT,IDEG,ITYP)
C
C           Check if labels is already defined
C           ==================================
            DO 80 ILAB = 1,NLAB
            IF(JLAB.EQ.IATTR(ILAB,1)) THEN
              IPLAB(NAORB,1) = ILAB
              GOTO 70
            ENDIF
 80         CONTINUE
C
C           Define new label
C           ================
            NLAB          = NLAB + 1
            IF(NLAB.GT.MAXLAB) GOTO 999
            DEG = '  '
            IF(NDEG.GT.1) WRITE(DEG,'(I2)') IDEG
            IATTR(NLAB,1) = JLAB
            PLABEL(NLAB,1) = CLS(ICLS)//' '//NAMN(ICENT)(1:3)//
     &         DEG//' '//GTOTYP(ITYP)
            IPLAB(NAORB,1) = NLAB
 70       CONTINUE
 60     CONTINUE
 50   CONTINUE
      NPLAB(0) = NLAB
      NPLAB(2) = NPLAB(0) - NPLAB(1)
C
C     Print section
C     =============
C
      IF(IPRINT.GE.1) THEN
        CALL HEADER('GETLAB: AO-labels',-1)
        WRITE(LUPRI,'(3X,A,I5)') '* Large components:',NPLAB(1)
        WRITE(LUPRI,'(6(3X,I3,2X,A12))') (I,PLABEL(I,1),I=1,NPLAB(1))
        IOFF = NPLAB(1)
        WRITE(LUPRI,'(3X,A,I5)') '* Small components:',NPLAB(2)
        WRITE(LUPRI,'(6(3X,I3,2X,A12))')
     &      (I,PLABEL(I,1),I=1+IOFF,NPLAB(0))
        IF(IPRINT.GE.7) THEN
          CALL PRIVEC('IPLAB ',IPLAB(1,1),NAORB)
        ENDIF
      ENDIF
C
C     **********************************************
C     ***** L A B E L S  for   S O - B A S I S *****
C     **********************************************
C
      NLAB = 0
      DO IC = 1,2
        DO IRP = 0,MAXREP
          DO 100 I = 1,NBBAS(IRP,IC)
            IU    = ICOS(IRP+1,IC) + I
            IS    = IBBAS(IRP,IC)  + I
            ICENT = IPCEN(IU)
            ITYP  = IPTYP(IU)
            JLAB  = IPACK(IC,ICENT,IRP,ITYP)
            DO ILAB = 1,NLAB
            IF(JLAB.EQ.IATTR(ILAB,2)) THEN
              IPLAB(IS,2) = ILAB
              GOTO 100
            ENDIF
            ENDDO
            NLAB = NLAB + 1
            IATTR(NLAB,2)  = JLAB
            PLABEL(NLAB,2) = CLS(IC)//' '//REP(IRP)
     &              //NAMN(ICENT)(1:3)//GTOTYP(ITYP)
            IPLAB(IS,2)    = NLAB
 100      CONTINUE
        ENDDO
      ENDDO      
C
C     Print section
C     =============
C
      IF((IPRINT.GE.1).AND.(NBSYM.GT.1)) THEN
        CALL HEADER('GETLAB: SO-labels',-1)
        WRITE(LUPRI,'(3X,A,I5)') '* Large components:',NPLAB(1)
        WRITE(LUPRI,'(6(3X,I3,2X,A12))') (I,PLABEL(I,2),I=1,NPLAB(1))
        IOFF = NPLAB(1)
        WRITE(LUPRI,'(3X,A,I5)') '* Small components:',NPLAB(2)
        WRITE(LUPRI,'(6(3X,I3,2X,A12))')
     &      (I,PLABEL(I,2),I=1+IOFF,NPLAB(0))
         IF(IPRINT.GE.7) THEN
           CALL PRIVEC('IPLAB ',IPLAB(1,2),NAORB)
        ENDIF
      ENDIF
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('GETLAB')
      RETURN
C
  999 CONTINUE
      WRITE (LUPRI,'(/A,I6/A,I6)')
     *   ' >>> PLABEL error, '//
     *   'no. of labels',NLAB,
     *   '                  '//
     *   ' current maximum number                  ',MAXLAB
      CALL QUIT('GETLAB:Too many labels !!!')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qnphase */
      SUBROUTINE QNPHASE(IM,JSYMOP,ONEINT,DOINT)
C*****************************************************************************
C
C     Quaternion phase insertion in lower triangular matrix
C       Phase: (e_K1*)(e_IM)(e_K2)
C
C     Theory: T.Saue &  H.J.Aa.Jensen
C     Written by T.Saue October 1994
C     Last revision: June 25 1996 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0)
C
#include <dgroup.h>
#include <dcbqal.h>
#include <dcbbas.h>
      LOGICAL DOINT(2,2),DIAG,SAME
      DIMENSION ONEINT(NNBBASX)
      EXTERNAL BDQALGEBRA
C
#include <ibtfun.h>
C
      IL = 2
      IR = 1
      I1 = 0
      DO 10 IRP1 = 0,NBSYM-1
        DO 20 IC1 = 1,2
          NB1 = NBBAS(IRP1,IC1)
          K1  = JQBAS(IRP1,IC1)
          I2 = 0
          DO 30 IRP2 = 0,IRP1
            SAME = IRP1.EQ.IRP2
            ICMX = 2
            IF(SAME) ICMX = IC1
            DO 40 IC2 = 1,ICMX
              NB2 = NBBAS(IRP2,IC2)
              K2  = JQBAS(IRP2,IC2)
              DIAG = SAME.AND.(IC1.EQ.IC2)
              IF(DOINT(IC1,IC2).AND.(IBTXOR(IRP1,IRP2).EQ.JSYMOP)) THEN
                KFAC = IQPHASE(K1,K2,IM)*IQPH(K1,IL)*IQPH(K2,IR)
                IF(KFAC.EQ.-1) THEN
                  DO 50 J1 = 1,NB1
                    JMX = NB2
                    IF(DIAG) JMX = J1
                    JOFF = (I1+J1)*(I1+J1-1)/2 + I2
                    DO 60 J2 = 1,JMX
                      ONEINT(JOFF+J2) = -ONEINT(JOFF+J2)
   60               CONTINUE
   50             CONTINUE
                ENDIF
              ENDIF
              I2 = I2 + NB2
   40       CONTINUE
   30     CONTINUE
          I1   = I1 + NB1
   20   CONTINUE
   10 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck q2phase */
      SUBROUTINE Q2PHASE(TYP,IM,JFAC,AMAT)
C*****************************************************************************
C
C     Quaternion phase insertion in square matrix
C       TYP = F:  (e_K1)(e_IM)(e_K2*)   as Fock matrix
C       TYP = D:  (e_K1*)(e_IM)(e_K2)   as density matrix
C
C     Written by T.Saue, January 1995
C     Last revision: June 25 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(DM1 = -1.0D0)
C
#include <mxorb.h>
#include <maxaqn.h>
#include <mxcent.h>
C
#include <symmet.h>
#include <dcbbas.h>
#include <dcbqal.h>
#include <dgroup.h>
      CHARACTER TYP*1
      DIMENSION AMAT(*)
      EXTERNAL BDQALGEBRA
C
      IF(TYP.NE.'F'.AND.TYP.NE.'D') THEN
        CALL QUIT('Q2PHASE: Unknown TYP !')
      ENDIF
      DO 10 IRP1 = 0,NBSYM-1
        ISY1 = IRP1 + 1
        DO 20 IC1 = 1,2
          NB1 = NBBAS(IRP1,IC1)
          IF(NB1.EQ.0) GOTO 20
          K1  = JQBAS(IRP1,IC1)
          DO 30 IRP2 = 0,NBSYM-1
            ISY2 = IRP2 + 1
            DO 40 IC2 = 1,2
              NB2 = NBBAS(IRP2,IC2)
              IF(NB2.EQ.0) GOTO 40
              K2   = JQBAS(IRP2,IC2)
              IF    (TYP.EQ.'F') THEN
                KFAC = JFAC*IQPHASE(K1,K2,IM)*IQPH(K1,1)*IQPH(K2,2)
              ELSEIF(TYP.EQ.'D') THEN
                IQ   = IQMULT(K1,K2,IM)
                KFAC = JFAC*IQPHASE(K1,K2,IQ)*IQPH(K1,2)*IQPH(K2,1)
              ENDIF
              IF(KFAC.EQ.-1) THEN
                IOFF = I2COSX(ISY1,ISY2,IC1,IC2)
                DO 50 ICOL = 1,NB2
                  CALL DSCAL(NB1,DM1,AMAT(IOFF+1),1)
                  IOFF = IOFF + NTBAS(0)
   50           CONTINUE
              ENDIF
   40       CONTINUE
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rpgen */
      SUBROUTINE RPGEN(IREPA)
C*****************************************************************************
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <symmet.h>
      DIMENSION IREPA(*)
      DO 10 ISYM = 1,NBSYM
        IREP = ISYM - 1
        DO 20 IC = 1,2
          CALL ICOPY(NCOS(ISYM,IC),IREP,0,IREPA(ICOS(ISYM,IC)+1),1)
 20     CONTINUE
 10   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtqoso */
      SUBROUTINE MTQOSO(AQO,ASO,IPRINT)
C*****************************************************************************
C
C     Matrix Transform: QO-basis to SO-basis
C     Written by Trond Saue
C     Last revision: Jan 7 1998 - jth LINUX
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
C
#include <dcblab.h>
#include <dgroup.h>
#include <dcbqal.h>
#include <dcbbas.h>
C  Used from COMMON blocks:
C     MXORB :  MXCORB
C     LABELS:  IPLAB,IATTR
C     DGROUP:  NZ,IQDEF
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AQO(NBAST,NBAST,NZ),ASO(NBAST,NBAST,4)
#else
      DIMENSION AQO(NTBAS(0),NTBAS(0),NZ),ASO(NTBAS(0),NTBAS(0),4)
#endif
      EXTERNAL BDQALGEBRA
#include <dcbibt.h>
      CALL DZERO(ASO,N2BBASX*4)
      DO K = 1,NZ
        IK = IPQTOQ(K,0)
        DO J = 1,NTBAS(0)
          IQJ  = JQBAS(IGET(IATTR(IPLAB(J,1),1)),
     &                 JGET(IATTR(IPLAB(J,1),1)))
          DO I = 1,NTBAS(0)
             IQI  = JQBAS(IGET(IATTR(IPLAB(I,1),1)),
     &                    JGET(IATTR(IPLAB(I,1),1)))
             IM   = IQMULT(IQI,IQJ,IK)
             ASO(I,J,IM) = (IQPHASE(IQI,IQJ,IK)*
     &          IQPH(IQI,1)*IQPH(IQJ,2))*AQO(I,J,K)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtsoqo */
      SUBROUTINE MTSOQO(ASO,AQO,IPRINT)
C*****************************************************************************
C
C     Matrix transform: SO-basis to QO-basis
C     Written by Trond Saue
C     Last revision: Jan 7 1998 -jth LINUX
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
C
#include <dcblab.h>
#include <dgroup.h>
#include <dcbqal.h>
#include <dcbbas.h>
C  Used from COMMON blocks:
C     MXORB :  MXCORB
C     LABELS:  IPLAB,IATTR
C     DGROUP:  NZ,IQDEF
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AQO(NBAST,NBAST,NZ),ASO(NBAST,NBAST,4)
#else
      DIMENSION AQO(NTBAS(0),NTBAS(0),NZ),ASO(NTBAS(0),NTBAS(0),4)
#endif
      EXTERNAL BDQALGEBRA
#include <dcbibt.h>
      CALL DZERO(AQO,N2BBASX*NZ)
      DO M = 1,NZ
        DO J = 1,NTBAS(0)
          IQJ  = JQBAS(IGET(IATTR(IPLAB(J,1),1)),
     &                 JGET(IATTR(IPLAB(J,1),1)))
          DO I = 1,NTBAS(0)
             IQI  = JQBAS(IGET(IATTR(IPLAB(I,1),1)),
     &                    JGET(IATTR(IPLAB(I,1),1)))
             IM    = IPQTOQ(M,0)
             IK    = IQMULT(IQI,IQJ,IM)
             AQO(I,J,M) = (IQPHASE(IQI,IQJ,IK)*IQPH(I,2)*IQPH(J,1))
     &                    *ASO(I,J,IK)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtqoso */
      SUBROUTINE VTQOSO(AQO,ASO,IPRINT)
C*****************************************************************************
C
C     Vector Transform: QO-basis to SO-basis
C     Last revision: Jan 7 1998 - jth LINUX
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
C
#include <dcblab.h>
#include <dcbqal.h>
#include <dgroup.h>
#include <dcbbas.h>
C  Used from COMMON blocks:
C     MXORB :  MXCORB
C     LABELS:  IPLAB,IATTR
C     DGROUP:  NZ,IQDEF
C     DCBBAS:  NTBAS(0),N2BBASX
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AQO(NBAST,NBAST,NZ),ASO(NBAST,NBAST,4)
#else
      DIMENSION AQO(NTBAS(0),NTBAS(0),NZ),ASO(NTBAS(0),NTBAS(0),4)
#endif
      EXTERNAL BDQALGEBRA
#include <dcbibt.h>
      CALL DZERO(ASO,N2BBASX*4)
      DO K = 1,NZ
        IK = IPQTOQ(K,0)
        DO J = 1,NTBAS(0)
          IQJ  = JQBAS(IGET(IATTR(IPLAB(J,1),1)),
     &                 JGET(IATTR(IPLAB(J,1),1)))
          DO I = 1,NTBAS(0)
             IQI  = JQBAS(IGET(IATTR(IPLAB(I,1),1)),
     &                    JGET(IATTR(IPLAB(I,1),1)))
             IM   = IQMULT(IQI,IQJ,IK)
             ASO(I,J,IM) = (IQPHASE(IQI,IQJ,IK)*
     &          IQPH(IQI,1)*IQPH(IQJ,2))*AQO(I,J,K)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtsoqo */
      SUBROUTINE VTSOQO(ASO,AQO,IPRINT)
C*****************************************************************************
C
C     Vector transform: SO-basis to QO-basis
C     Written by Trond Saue
C     Last revision: Jan 7 1998 - jth LINUX
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
C
#include <dcblab.h>
#include <dgroup.h>
#include <dcbqal.h>
#include <dcbbas.h>
C  Used from COMMON blocks:
C     MXORB :  MXCORB
C     LABELS:  IPLAB,IATTR
C     DGROUP:  NZ,IQDEF
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AQO(NBAST,NBAST,NZ),ASO(NBAST,NBAST,4)
#else
      DIMENSION AQO(NTBAS(0),NTBAS(0),NZ),ASO(NTBAS(0),NTBAS(0),4)
#endif
      EXTERNAL BDQALGEBRA
#include <dcbibt.h>
      CALL DZERO(AQO,N2BBASX*NZ)
      DO M = 1,NZ
        DO J = 1,NTBAS(0)
          IQJ  = JQBAS(IGET(IATTR(IPLAB(J,1),1)),
     &                 JGET(IATTR(IPLAB(J,1),1)))
          DO I = 1,NTBAS(0)
             IQI  = JQBAS(IGET(IATTR(IPLAB(I,1),1)),
     &                    JGET(IATTR(IPLAB(I,1),1)))
             IM    = IPQTOQ(M,0)
             IK    = IQMULT(IQI,IQJ,IM)
             AQO(I,J,M) = (IQPHASE(IQI,IQJ,IK)*IQPH(I,2)*IQPH(J,1))
     &                   *ASO(I,J,IK)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtsoao */
      SUBROUTINE MTSOAO(ASO,AAO,NBAST,IREPDM,IPRINT)
C*****************************************************************************
C
C     Matrix transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1=1.0D0)
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <pgroup.h>
#include <ibtfun.h>
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine MTSOAO',-1)
C
C     Loop over all irreps in molecule
C
      IR(1) = IREPDM
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(2) = IBTXOR(IREPDM,IFXYZ)
      ISTRA = 1
      CALL DZERO(AAO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
            AN     = D1/FMULT(MULA)
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA = AN*
     &              PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
CTOND               FACA = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     BN     = D1/FMULT(MULB)
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = BN*PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        AAO(INDA,INDB) = AAO(INDA,INDB)
     &                                 + FACA*FACB*ASO(I,J)
                     END IF
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('OUTPUT: Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtaoso */
      SUBROUTINE MTAOSO(AAO,ASO,NBAST,IREPDM,IPRINT)
C*****************************************************************************
C
C     Matrix transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1=1.0D0)
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <pgroup.h>
#include <ibtfun.h>
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine MTAOSO',-1)
C
C     Loop over all irreps in molecule
C
      IR(1) = IREPDM
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(2) = IBTXOR(IREPDM,IFXYZ)
      ISTRA = 1
C      CALL DZERO(AAO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
            AN     = D1/FMULT(MULA)
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     BN     = D1/FMULT(MULB)
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        ASO(I,J) = ASO(I,J)
     &                                 + FACA*FACB*AAO(INDA,INDB)
                     END IF
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT:  Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('OUTPUT: Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dtsoao */
      SUBROUTINE DTSOAO(ASO,AAO,NBAST,IREPDM,IPRINT)
C*****************************************************************************
C
C     Matrix transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1=1.0D0)
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <pgroup.h>
#include <ibtfun.h>
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine DTSOAO',-1)
C
C     Loop over all irreps in molecule
C
      IR(1) = IREPDM
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(2) = IBTXOR(IREPDM,IFXYZ)
      ISTRA = 1
      CALL DZERO(AAO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        AAO(INDA,INDB) = AAO(INDA,INDB)
     &                                 + FACA*FACB*ASO(I,J)
                     END IF
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('OUTPUT:Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dtaoso */
      SUBROUTINE DTAOSO(AAO,ASO,NBAST,IREPDM,IPRINT)
C*****************************************************************************
C
C     Matrix transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1=1.0D0)
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <pgroup.h>
#include <ibtfun.h>
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine DTAOSO',-1)
C
C     Loop over all irreps in molecule
C
      IR(1) = IREPDM
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(2) = IBTXOR(IREPDM,IFXYZ)
      ISTRA = 1
C      CALL DZERO(AAO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            AN     = D1/FMULT(MULA)
            INDA   = KSTRT(IA) + NA - KHKTA
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA = AN*
     &              PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     BN     = D1/FMULT(MULB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = BN*PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        ASO(I,J) = ASO(I,J)
     &                                 + FACA*FACB*AAO(INDA,INDB)
                     END IF
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('OUTPUT:Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ntaoso */
      SUBROUTINE NTAOSO(AAO,ASO,NBAST,IREPM,IPRINT)
C*****************************************************************************
C
C     Matrix Transform:  AO-basis to SO-basis
C     This routine is based on HERMIT routine DAOTSO (880418  PRT)
C     Note that the SO-matrix is NOT initialized !!!
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sept 10 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0)
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <aosotr.h>
#include <pgroup.h>
#include <ibtfun.h>
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine MTAOSO',-1)
C
C     Loop over all irreps in molecule
C
      IFXYZ  = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(1) = IREPM
      IR(2) = IBTXOR(IREPM,IFXYZ)
      ISTRA = 1
C      CALL DZERO(ASO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            JKBA   = JTRAN(I)
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            MULA   = ISTBAO(IA)
            AN     = D1/FMULT(MULA)
            DO 300 IKBA = 1,JKBA
               INDA = ITRAN(I,IKBA)
               FACA = CTRAN(I,IKBA)
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     JKBB = JTRAN(J)
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     MULB   = ISTBAO(IB)
                     BN     = D1/FMULT(MULB)
                     DO 600 IKBB = 1,JKBB
                        INDB = ITRAN(J,IKBB)
                        FACB = CTRAN(J,IKBB)
                        ASO(I,J) = ASO(I,J)
     &                                 + FACA*FACB*AAO(INDA,INDB)
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtsoao */
      SUBROUTINE VTSOAO(ASO,AAO,NBAST,NVEC,IPRINT)
C*****************************************************************************
C
C     Vector transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      DIMENSION ASO(NBAST,NVEC), AAO(NBAST,NVEC)
      DIMENSION IR(2)
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <pgroup.h>
#include <ibtfun.h>
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine VTSOAO',-1)
C
C     Loop over all irreps in molecule
C
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      ISTRA = 1
      CALL DZERO(AAO,NBAST*NVEC)
      DO 100 IREPA = 0, MAXREP
        NORBA = NAOS(IREPA+1)
        DO 200 I = ISTRA,ISTRA + NORBA - 1
          IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
          NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
          NHKTA  = NHKT(IA)
          KHKTA  = KHKT(IA)
          MULA   = ISTBAO(IA)
          INDA   = KSTRT(IA) + NA - KHKTA
C
C         Loop over symmetry dependent centers
C
          DO 300 ISYMA = 0, MAXOPR
          IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
            INDA = INDA + KHKTA
            FAC  = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
            CALL DAXPY(NVEC,FAC,ASO(I,1),NBAST,AAO(INDA,1),NBAST)
          ENDIF
  300     CONTINUE
  200   CONTINUE
        ISTRA = ISTRA + NORBA
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Vectors in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NVEC,NBAST,NVEC,1,LUPRI)
         CALL HEADER('INPUT: Vectors in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NVEC,NBAST,NVEC,1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wtsoao */
      SUBROUTINE WTSOAO(ASO,AAO,IFRP,KNFBAS,NBAST,NVEC,IPRINT)
C*****************************************************************************
C
C     Vector transform: SO-basis to AO-basis (in sorted basis !)
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Jan 7 1998 (jth: NFBAS(IFRP,0) now parameter)
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
C
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <pgroup.h>
#include <dgroup.h>
#include <dcbbas.h>
C
C     KNFBAS = NFBAS(IFRP,0)
C     NBAST  = NTBAS(0)
C
      DIMENSION ASO(KNFBAS,NVEC), AAO(NBAST,NVEC)
      DIMENSION IR(2)
C
#include <ibtfun.h>
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine WTSOAO',-1)
C
C     Loop over all irreps in molecule
C
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      CALL DZERO(AAO,NTBAS(0)*NVEC)
      DO 100 IREPA = 0, MAXREP
        DO IC = 1,2
        NORBA = NCOS(IREPA+1,IC)
        IF(JBTOF(IREPA,IC).EQ.IFRP) THEN
          DO 200 I = 1,NORBA
            IU     = ICOS(IREPA+1,IC) + I
            IS     = IBBAS(IREPA,IC) - IBAS(IFRP) + I
            IA     = IBTAND(IBTSHR(IPIND(IU),16),65535)
            NA     = IBTAND(IBTSHR(IPIND(IU), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
C
C           Loop over symmetry dependent centers
C
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
              INDA = INDA + KHKTA
              FAC  = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
              CALL DAXPY(NVEC,FAC,ASO(IS,1),NFBAS(IFRP,0),
     &                            AAO(INDA,1),NTBAS(0))
            ENDIF
  300       CONTINUE
  200     CONTINUE
        ENDIF
        ENDDO
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Vectors in SO basis',-1)
         CALL OUTPUT(ASO,1,NFBAS(IFRP,0),1,NVEC,
     &               NFBAS(IFRP,0),NVEC,1,LUPRI)
         CALL HEADER('INPUT: Vectors in AO basis',-1)
         CALL OUTPUT(AAO,1,NTBAS(0),1,NVEC,NTBAS(0),NVEC,1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtaoso */
      SUBROUTINE VTAOSO(AAO,ASO,NBAST,NVEC,IPRINT)
C*****************************************************************************
C
C     Vector Transform:  AO-basis to SO-basis
C     This routine is based on HERMIT routine DAOTSO (880418  PRT)
C     Note that the SO-matrix is NOT initialized !!!
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sept 10 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0)
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      DIMENSION ASO(NBAST,NVEC), AAO(NBAST,NVEC)
      DIMENSION IR(2)
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <pgroup.h>
#include <ibtfun.h>
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine VTAOSO',-1)
C
C     Loop over all irreps in molecule
C
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      ISTRA = 1
      CALL DZERO(ASO,NBAST*NVEC)
      DO 100 IREPA = 0, MAXREP
        NORBA = NAOS(IREPA+1)
        DO 200 I = ISTRA,ISTRA + NORBA - 1
          IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
          NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
          NHKTA  = NHKT(IA)
          KHKTA  = KHKT(IA)
          MULA   = ISTBAO(IA)
          AN     = D1/FMULT(MULA)
          INDA   = KSTRT(IA) + NA - KHKTA
C
C         Loop over symmetry dependent centers
C
          DO 300 ISYMA = 0, MAXOPR
          IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
            INDA = INDA + KHKTA
            FAC  = AN*PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
            CALL DAXPY(NVEC,FAC,AAO(INDA,1),NBAST,ASO(I,1),NBAST)
          ENDIF
  300     CONTINUE
  200   CONTINUE
        ISTRA = ISTRA + NORBA
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Vectors in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NVEC,NBAST,NVEC,1,LUPRI)
         CALL HEADER('OUTPUT: Vectors in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NVEC,NBAST,NVEC,1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bstobu */
      SUBROUTINE BSTOBU(AMAT,MZ,WORK,LWORK)
C*****************************************************************************
C
C     Matrix transform:  from sorted basis to unsorted basis
C
C     Written by T.Saue Mar 2 1995
C     Last revision : May 6 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <symmet.h>
C
      DIMENSION AMAT(*),WORK(LWORK)
      CALL QENTER('BSTOBU')
#include <memint.h>
C
C     Memory allocation
C
      CALL MEMGET('REAL',KBUF,MZ*N2BBASX,WORK,KFREE,LFREE)
      CALL DCOPY(MZ*N2BBASX,AMAT,1,WORK(KBUF),1)
      IZOFF = 0
      DO IZ = 1,MZ
        DO JC = 1,2
          DO JREP = 0,MAXREP
            DO J = 0,NBBAS(JREP,JC)-1
              JSOFF = IZOFF + (IBBAS(JREP  ,JC)+J)*NTBAS(0)
              JUOFF = IZOFF + (ICOS (JREP+1,JC)+J)*NTBAS(0)
              DO IC = 1,2
                DO IREP = 0,MAXREP
                  IUOFF = JUOFF + ICOS (IREP+1,IC)
                  ISOFF = JSOFF + IBBAS(IREP  ,IC)
                  IF(NBBAS(IREP,IC).GT.0)
     &              CALL DCOPY(NBBAS(IREP,IC),WORK(KBUF+ISOFF),1,
     &                                        AMAT(IUOFF+1)   ,1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IZOFF = IZOFF + N2BBASX
      ENDDO
C
C     Memory deallocation
C
      CALL MEMREL('BSTOBU',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT ('BSTOBU')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck butobs */
      SUBROUTINE BUTOBS(AMAT,MZ,WORK,LWORK)
C*****************************************************************************
C
C     Matrix transform:  from sorted basis to unsorted basis
C
C     Written by T.Saue Mar 2 1995
C     Last revision : May 6 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <symmet.h>
C
      DIMENSION AMAT(*),WORK(LWORK)
      CALL QENTER('BSTOBU')
#include <memint.h>
C

      CALL MEMGET('REAL',KBUF,MZ*N2BBASX,WORK,KFREE,LFREE)
      CALL DCOPY(MZ*N2BBASX,AMAT,1,WORK(KBUF),1)
      IZOFF = 0
      DO IZ = 1,MZ
        DO JC = 1,2
          DO JREP = 0,MAXREP
            DO J = 0,NBBAS(JREP,JC)-1
              JSOFF = IZOFF + (IBBAS(JREP  ,JC)+J)*NTBAS(0)
              JUOFF = IZOFF + (ICOS (JREP+1,JC)+J)*NTBAS(0)
              DO IC = 1,2
                DO IREP = 0,MAXREP
                  IUOFF = JUOFF + ICOS (IREP+1,IC)
                  ISOFF = JSOFF + IBBAS(IREP  ,IC)
                  IF(NBBAS(IREP,IC).GT.0) 
     &              CALL DCOPY(NBBAS(IREP,IC),WORK(KBUF+IUOFF),1,
     &                                        AMAT(ISOFF+1)   ,1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IZOFF = IZOFF + N2BBASX
      ENDDO
C
      CALL MEMREL('BSTOBU',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT ('BSTOBU')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtbubs */
      SUBROUTINE VTBUBS(VBU,VBS,NVEC)
C*****************************************************************************
C
C     Vector transform:  from unsorted basis to sorted basis
C
C     Written by T.Saue Mar 2 1995
C     Last revision : Jan 7 1998 - jth LINUX
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <symmet.h>
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION VBU(NBAST,NVEC),VBS(NBAST,NVEC)
#else
      DIMENSION VBU(NTBAS(0),NVEC),VBS(NTBAS(0),NVEC)
#endif
      CALL QENTER('VTBUBS')
C
      DO J = 1,NVEC
        DO IC = 1,2
          DO IREP = 0,MAXREP
            IUOFF = ICOS (IREP+1,IC)
            ISOFF = IBBAS(IREP  ,IC)
            IF(NBBAS(IREP,IC).GT.0) 
     &        CALL DCOPY(NBBAS(IREP,IC),VBU(IUOFF+1,J),1,
     &                                  VBS(ISOFF+1,J),1)
          ENDDO
        ENDDO
      ENDDO
C
      CALL QEXIT ('VTBUBS')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtbsbu */
      SUBROUTINE VTBSBU(VBS,VBU,NVEC)
C*****************************************************************************
C
C     Vector transform:  from sorted basis to unsorted basis
C
C     Written by T.Saue Mar 2 1995
C     Last revision : Jan 7 1998 - jth LINUX
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <symmet.h>
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION VBS(NBAST,NVEC),VBU(NBAST,NVEC)
#else
      DIMENSION VBS(NTBAS(0),NVEC),VBU(NTBAS(0),NVEC)
#endif
      CALL QENTER('VTBSBU')
      DO J = 1,NVEC
        DO IC = 1,2
          DO IREP = 0,MAXREP
            IUOFF = ICOS (IREP+1,IC)
            ISOFF = IBBAS(IREP  ,IC)
            IF(NBBAS(IREP,IC).GT.0)
     &        CALL DCOPY(NBBAS(IREP,IC),VBS(ISOFF+1,J),1,
     &                                  VBU(IUOFF+1,J),1)
          ENDDO
        ENDDO
      ENDDO
      CALL QEXIT ('VTBSBU')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtbsbu */
      SUBROUTINE MTBSBU(AMS,AMU)
C*****************************************************************************
C
C     Matrix transform:  from sorted basis to unsorted basis
C
C     Written by T.Saue Mar 2 1995
C     Last revision : July 17 1997 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <symmet.h>
C
      DIMENSION AMS(*),AMU(*)
      CALL QENTER('MTBSBU')
C
      DO JC = 1,2
        DO JREP = 0,MAXREP
          DO J = 0,NBBAS(JREP,JC)-1
            JSOFF = (IBBAS(JREP  ,JC)+J)*NTBAS(0)
            JUOFF = (ICOS (JREP+1,JC)+J)*NTBAS(0)
            DO IC = 1,2
              DO IREP = 0,MAXREP
                IUOFF = JUOFF + ICOS (IREP+1,IC)
                ISOFF = JSOFF + IBBAS(IREP  ,IC)
                IF(NBBAS(IREP,IC).GT.0)
     &            CALL DCOPY(NBBAS(IREP,IC),AMS(ISOFF+1),1,
     &                                      AMU(IUOFF+1),1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL QEXIT ('MTBSBU')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtbubs */
      SUBROUTINE MTBUBS(AMU,AMS)
C*****************************************************************************
C
C     Matrix transform:  from unsorted basis to sorted basis
C
C     Written by T.Saue Mar 2 1995
C     Last revision : July 17 1997 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <symmet.h>
C
      DIMENSION AMU(*),AMS(*)
      CALL QENTER('MTBUBS')
      DO JC = 1,2
        DO JREP = 0,MAXREP
          DO J = 0,NBBAS(JREP,JC)-1
            JSOFF = (IBBAS(JREP  ,JC)+J)*NTBAS(0)
            JUOFF = (ICOS (JREP+1,JC)+J)*NTBAS(0)
            DO IC = 1,2
              DO IREP = 0,MAXREP
                IUOFF = JUOFF + ICOS (IREP+1,IC)
                ISOFF = JSOFF + IBBAS(IREP  ,IC)
                IF(NBBAS(IREP,IC).GT.0)
     &            CALL DCOPY(NBBAS(IREP,IC),AMU(IUOFF+1),1,
     &                                      AMS(ISOFF+1),1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL QEXIT ('MTBUBS')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sphgen */
      SUBROUTINE SPHGEN(WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Generate coefficients for transforming from Cartesian to
C     spherical components
C
C     This is just a driver exploiting all the goodies in SPHCOM of
C     the HERMIT package
C
C     Written by T.Saue Sep 19 1995
C     Last revision: Sep 19 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0)
#include <maxaqn.h>
C
C Used from COMMON blocks
C  SPHCOM: CSP,ISPADR
C
#include <dcbsph.h>

      DIMENSION WORK(LWORK)
C
      CALL QENTER('SPHGEN')
#include <memint.h>
      NCSPT = (MXQN*(MXQN + 1)/2)*(2*MXQN - 1)
      CALL DZERO(CSP,NCSP)
C
      CALL MEMGET('REAL',KCSPT,NCSPT,WORK,KFREE,LFREE)
C
C     MORDER eq 0 : M_l order 0, +1, -1, +2, -2, ..., +LVAL, -LVAL
C     MORDER ne 0 : M_l order -LVAL, ..., -1, 0, 1, ..., +LVAL
C
C     MINTEG eq 1 : the spherical components are normalized such
C                   that all coefficients are integers.
      MORDER = 1
      MINTEG = 1
      IOFF = 1
      DO 10 LVAL = 0, MXQN - 1
        NLM            = 2*LVAL + 1
        NXYZ           = (LVAL+1)*(LVAL+2)/2
        ISPADR(LVAL+1) = IOFF
        IOFF           = IOFF + NXYZ*NXYZ
        IF (LVAL .EQ. 0) THEN
          CSP(1)       = D1
        ELSE IF (LVAL .EQ. 1) THEN
          CSP (2)      = D1
          CSP (6)      = D1
          CSP(10)      = D1
        ELSE
          CALL SPHCOM(LVAL,WORK(KCSPT),NLM,NXYZ,MORDER,MINTEG,
     &                WORK(KFREE),LFREE,IPRINT)
C
C         Normalize matrix
C
        IF(LVAL.EQ.3) THEN
C
C     Jon K. Laerdahl (16.12.96):
C         Normalization NOT necessary
C         Makes RKB transformation more complicated
C         --> Remove Normalization
C
C          JOFF = KCSPT
C1
C          FAC = D1/SQRT(24.0000D0)
C          CALL DSCAL(NXYZ,FAC,WORK(JOFF),1)
C          JOFF = JOFF + NXYZ
C2
C          JOFF = JOFF + NXYZ
C3
C          FAC = D1/SQRT(40.0D0)
C          CALL DSCAL(NXYZ,FAC,WORK(JOFF),1)
C          JOFF = JOFF + NXYZ
C4
C          FAC = D1/SQRT(6.666666666670D0)
C          CALL DSCAL(NXYZ,FAC,WORK(JOFF),1)
C          JOFF = JOFF + NXYZ
C5
C          FAC = D1/SQRT(40.0D0)
C          CALL DSCAL(NXYZ,FAC,WORK(JOFF),1)
C          JOFF = JOFF + NXYZ
C6
C          FAC = D1/SQRT(4.0D0)
C          CALL DSCAL(NXYZ,FAC,WORK(JOFF),1)
C          JOFF = JOFF + NXYZ
C7
C          FAC = D1/SQRT(2.66666667D0)
C          CALL DSCAL(NXYZ,FAC,WORK(JOFF),1)
C          JOFF = JOFF + NXYZ
          CALL MTRSP(NXYZ,NLM,WORK(KCSPT),NXYZ,CSP(ISPADR(LVAL+1)),NLM)
        ELSE
C          JOFF = KCSPT
C          DO K = 1,NLM
C            FAC = DDOT(NXYZ,WORK(JOFF),1,WORK(JOFF),1)
C            FAC = D1/SQRT(FAC)
C            CALL DSCAL(NXYZ,FAC,WORK(JOFF),1)
C            JOFF = JOFF + NXYZ
C          ENDDO
          CALL MTRSP(NXYZ,NLM,WORK(KCSPT),NXYZ,CSP(ISPADR(LVAL+1)),NLM)
        END IF
        END IF
 10   CONTINUE
C
C     Generate backtransformation as well
C
      DO L = 1,MXQN
        IOFF = ISPADR(L)
        NLM  = 2*L- 1
        NXYZ = (L*(L+1))/2
        CALL MTRSP(NLM,NXYZ,CSP(IOFF),NLM,SPC(IOFF),NXYZ)
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.4) THEN
        CALL HEADER('SPHGEN: Cartesian transformation matrix:',-1)
        DO L = 1,MXQN
          LVAL = L-1
          NLM            = 2*LVAL + 1
          NXYZ           = (LVAL+1)*(LVAL+2)/2
          WRITE (LUPRI,'(1X,A,I5)')
     &     ' Coefficients for angular quantum number ',LVAL
          CALL OUTPUT(CSP(ISPADR(L)),1,NLM,1,NXYZ,NLM,NXYZ,1,LUPRI)
        ENDDO
        CALL HEADER('SPHGEN: Spherical transformation matrix:',-1)
        DO L = 1,MXQN
          LVAL = L-1
          NLM            = 2*LVAL + 1
          NXYZ           = (LVAL+1)*(LVAL+2)/2
          WRITE (LUPRI,'(1X,A,I5)')
     &     ' Coefficients for angular quantum number ',LVAL
          CALL OUTPUT(SPC(ISPADR(L)),1,NXYZ,1,NLM,NXYZ,NLM,1,LUPRI)
        ENDDO
      ENDIF
      CALL MEMREL('SPHGEN.SPHGE1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('SPHGEN')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtaocs */
      SUBROUTINE VTAOCS(AOV,LDA,CSV,LDC,NVEC,LCOMP,IPRINT)
C*****************************************************************************
C
C     Transform NVEC vectors from Cartesian AO basis to Spherical
C     Components
C
C     Options:
C     LCOMP = 0   Transform all components
C     LCOMP = 1   Transform only large components
C     LCOMP = 2   Transform only small components
C
C     NSPH - number of spherical components
C
C     Written by T.Saue Sep 19 1995
C     Last revision Sep 19 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
C
C Used from COMMON blocks:
C    SHELLS: KMAX
C    SPHCOM:
C    NUCLEI: NLARGE,NSMALL,NBASIS
#include <nuclei.h>
#include <shells.h>
#include <dcbsph.h>
      DIMENSION AOV(LDA,NVEC),CSV(LDC,NVEC)
C
C     Set limits
C
      NSPH = NSPHCM(LCOMP)
      IF    (LCOMP.EQ.0) THEN
C       All components
        ISH1 = 1
        ISH2 = KMAX
      ELSEIF(LCOMP.EQ.1) THEN
C       Large components
        ISH1 = 1
        ISH2 = NLRGSH
      ELSEIF(LCOMP.EQ.2) THEN
C       Small components
        ISH1 = NLRGSH+1
        ISH2 = KMAX
      ELSE
        WRITE(LUPRI,'(A,I5)') '* Unknown LCOMP = ',LCOMP
        CALL QUIT('VTAOCS: Unknown component LCOMP !!!')
      ENDIF
C
C     Transform to spherical components
C
      IOFF = 1
      JOFF = 1
      DO ISHELL = ISH1,ISH2
        L     = NHKT(ISHELL)
        NLM   = 2*(L-1) + 1
        NXYZ  = L*(L+1)/2
        NDEG   = NUCDEG(NCENT(ISHELL))
        DO IDEG = 1,NDEG
          CALL DGEMM('N','N',NLM,NVEC,NXYZ,D1,CSP(ISPADR(L)),NLM,
     &               AOV(IOFF,1),LDA,D0,CSV(JOFF,1),LDC)
          IOFF = IOFF + NXYZ
          JOFF = JOFF + NLM
        ENDDO
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.5) THEN
        CALL HEADER('VTAOCS: Transformed vectors:',-1)
        CALL OUTPUT(CSV,1,NSPH,1,NVEC,LDC,NVEC,1,LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtaocs */
      SUBROUTINE MTAOCS(AOM,LDA,CSM,LDC,LCOMP,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Transform NVEC vectors from Cartesian AO basis to Spherical
C     Components
C
C     Options:
C     LCOMP = 0   Transform all components
C     LCOMP = 1   Transform only large components
C     LCOMP = 2   Transform only small components
C
C     NSPH - total number of spherical components
C
C     Written by T.Saue Sep 19 1995
C     Last revision Sep 19 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
#include <mxcent.h>
C
#include <nuclei.h>
      DIMENSION AOM(LDA,*),CSM(LDC,*),WORK(LWORK)
      CALL QENTER('MTAOCS')
C
C     Set limits
C
      NSPH = NSPHCM(LCOMP)
      IF    (LCOMP.EQ.0) THEN
C       All components
        NVEC = NBASIS
      ELSEIF(LCOMP.EQ.1) THEN
C       Large components
        NVEC = NLARGE
      ELSEIF(LCOMP.EQ.2) THEN
C       Small components
        NVEC = NSMALL
      ELSE
        WRITE(LUPRI,'(A,I5)') '* Unknown LCOMP = ',LCOMP
        CALL QUIT('MTAOCS: Unknown component LCOMP !!!')
      ENDIF
C
C     Allocate work space
C
#include <memint.h>
      NWRK = NVEC*NSPH
      CALL MEMGET('REAL',KW1,NWRK,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KW2,NWRK,WORK,KFREE,LFREE)
C
C     First transform rows
C
      CALL VTAOCS(AOM,LDA,WORK(KW1),NSPH,NVEC,LCOMP,0)
      IF(IPRINT.GE.6) THEN
        CALL HEADER('MTAOCS: Half-transformed matrix:',-1)
        CALL OUTPUT(WORK(KW1),1,NSPH,1,NVEC,NSPH,NVEC,1,LUPRI)
      ENDIF
C
C     Transpose and transform columns
C
      CALL MTRSP(NSPH,NVEC,WORK(KW1),NSPH,WORK(KW2),NVEC)
      CALL VTAOCS(WORK(KW2),NVEC,WORK(KW1),NSPH,NSPH,LCOMP,0)
C
C     Backtranspose to obtain final matrix
C
      CALL MTRSP(NSPH,NSPH,WORK(KW1),NSPH,CSM,LDC)
C
C     Memory deallocation
      CALL MEMREL('MTAOCS',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Print section
C
      IF(IPRINT.GE.5) THEN
        CALL HEADER('MTAOCS: Transformed matrix:',-1)
        CALL OUTPUT(CSM,1,NSPH,1,NSPH,LDC,NSPH,1,LUPRI)
      ENDIF
      CALL QEXIT ('MTAOCS')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vtcsao */
      SUBROUTINE VTCSAO(CSV,LDC,AOV,LDA,NVEC,LCOMP,IPRINT)
C*****************************************************************************
C
C     Transform NVEC vectors from Spherical Components to Cartesian
C     AO basis.
C     NOTE: The transformation matrix is generally rectanuglar,
C           so the extra components are lost !!!!
C
C     Options:
C     LCOMP = 0   Transform all components
C     LCOMP = 1   Transform only large components
C     LCOMP = 2   Transform only small components
C
C     NSPH - number of spherical components
C
C     Written by T.Saue Sep 19 1995
C     Last revision Sep 19 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
C
C Used from COMMON blocks:
C    SHELLS: KMAX
C    SPHCOM:
C    NUCLEI: NLARGE,NSMALL,NBASIS
#include <nuclei.h>
#include <shells.h>
#include <dcbsph.h>
      DIMENSION AOV(LDA,NVEC),CSV(LDC,NVEC)
C
C     Set limits
C
      NSPH = NSPHCM(LCOMP)
      IF    (LCOMP.EQ.0) THEN
C       All components
        ISH1 = 1
        ISH2 = KMAX
      ELSEIF(LCOMP.EQ.1) THEN
C       Large components
        ISH1 = 1
        ISH2 = NLRGSH
      ELSEIF(LCOMP.EQ.2) THEN
C       Small components
        ISH1 = NLRGSH+1
        ISH2 = KMAX
      ELSE
        WRITE(LUPRI,'(A,I5)') '* Unknown LCOMP = ',LCOMP
        CALL QUIT('VTCSAO: Unknown component LCOMP !!!')
      ENDIF
C
C     Transform to spherical components
C
      IOFF = 1
      JOFF = 1
      DO ISHELL = ISH1,ISH2
        L     = NHKT(ISHELL)
        NLM   = 2*(L-1) + 1
        NXYZ  = L*(L+1)/2
        NDEG   = NUCDEG(NCENT(ISHELL))
        DO IDEG = 1,NDEG
          CALL DGEMM('T','N',NXYZ,NVEC,NLM,D1,CSP(ISPADR(L)),NLM,
     &               CSV(JOFF,1),LDC,D0,AOV(IOFF,1),LDA)
          IOFF = IOFF + NXYZ
          JOFF = JOFF + NLM
        ENDDO
      ENDDO
      IF(IPRINT.GE.5) THEN
        NDIM = IOFF-1
        CALL HEADER('VTCSAO: Transformed vectors:',-1)
        CALL OUTPUT(AOV,1,NDIM,1,NVEC,LDA,NVEC,1,LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck nsphcm */
      FUNCTION NSPHCM(LCOMP)
C*****************************************************************************
C
C     Count the number of spherical functions
C
C     Options:
C     LCOMP = 0   All components
C     LCOMP = 1   Only large components
C     LCOMP = 2   Only small components
C
C     Written by T.Saue Sep 19 1995
C     Last revision Sep 19 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
#include <mxorb.h>
#include <mxcent.h>
C
C Used from COMMON blocks:
C    SHELLS: KMAX,NHKT
C    NUCLEI: NLARGE,NSMALL,NBASIS
#include <nuclei.h>
#include <shells.h>
C
C     Set limits
C
      IF    (LCOMP.EQ.0) THEN
C       All components
        ISH1 = 1
        ISH2 = KMAX
      ELSEIF(LCOMP.EQ.1) THEN
C       Large components
        ISH1 = 1
        ISH2 = NLRGSH
      ELSEIF(LCOMP.EQ.2) THEN
C       Small components
        ISH1 = NLRGSH+1
        ISH2 = KMAX
      ELSE
        WRITE(LUPRI,'(A,I5)') '* Unknown LCOMP = ',LCOMP
        CALL QUIT('NSPHCM: Unknown component LCOMP !!!')
      ENDIF
C
C     Count the number of spherical functions
C
      NSPHCM = 0
      DO ISHELL = ISH1,ISH2
        L      = NHKT(ISHELL)
        NDEG   = NUCDEG(NCENT(ISHELL))
        NSPHCM = NSPHCM + (2*(L-1) + 1)*NDEG
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sphtst */
      SUBROUTINE SPHTST(WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Check normalization of tranformation matrix
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
C
#include <dcbsph.h>
#include <dcbbas.h>
      DIMENSION CSPSUM(((MXQN+1)*MXQN)/2),WORK(LWORK)
#include <memint.h>
C
      CALL HEADER('Output from SPHTST',-1)
      CALL SPHGEN(WORK,LWORK,IPRINT)
C
C     Check normalization
C
      DO L = 1,MXQN
        LVAL = L-1
        NLM            = 2*LVAL + 1
        NXYZ           = (LVAL+1)*(LVAL+2)/2
        IOFF = ISPADR(L)
        DO K = 1,NXYZ
          CSPSUM(K) = SQRT(DDOT(NLM,CSP(IOFF),1,CSP(IOFF),1))
          IOFF = IOFF + NLM
        ENDDO
        WRITE(LUPRI,'(I5,5X,40F9.4)') (L-1),(CSPSUM(I),I=1,NXYZ)
      ENDDO
      NSPH = NSPHCM(0)
      WRITE(LUPRI,'(A,I5)') 'NSPH: ',NSPH
      CALL MEMGET('REAL',KA,NTBAS(0)*NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KB,NSPH*NTBAS(0),WORK,KFREE,LFREE)
      CALL SPHTS1(WORK(KA),WORK(KB),NSPH,WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('SPHTST',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sphts1 */
      SUBROUTINE SPHTS1(AMAT,BMAT,NSPH,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION AMAT(NBAST,NBAST),BMAT(NSPH,NBAST),WORK(LWORK)
#else
      DIMENSION AMAT(NTBAS(0),NTBAS(0)),BMAT(NSPH,NTBAS(0)),WORK(LWORK)
#endif
C
C     GENERATE MATRIX
C
      OFF = 0.0D0
      DO J = 1,NTBAS(0)
        DO I = 1,NTBAS(0)
          OFF = OFF + 0.1D0
          AMAT(I,J) = OFF
        ENDDO
      ENDDO
      CALL HEADER('SPHTS1: Test matrix',-1)
      CALL OUTPUT(AMAT,1,NTBAS(0),1,NTBAS(0),NTBAS(0),NTBAS(0),
     &            1,LUPRI)
      CALL VTAOCS(AMAT,NTBAS(0),BMAT,NSPH,NTBAS(0),0,IPRINT)
      CALL VTCSAO(BMAT,NSPH,AMAT,NTBAS(0),NTBAS(0),0,IPRINT)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mmsoao */
      SUBROUTINE MMSOAO(ASO,AAO,NBAST,IPRINT)
C*****************************************************************************
C
C     Matrix transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1=1.0D0)
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <pgroup.h>
#include <ibtfun.h>
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine MTSOAO',-1)
C
C     Loop over all irreps in molecule
C
      ISTRA = 1
      CALL DZERO(AAO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
            AN     = D1/FMULT(MULA)
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA =
     &              PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
C                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     BN     = D1/FMULT(MULB)
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        AAO(INDA,INDB) = AAO(INDA,INDB)
     &                                 + FACA*FACB*ASO(I,J)
                     END IF
  600                CONTINUE
  500             CONTINUE
C                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('OUTPUT: Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qtrans */
      SUBROUTINE QTRANS(TYP,TREV,FADD,NRAO,NCAO,NRMO,NCMO,
     &                  FAO,LRAO,LCAO,NZAO,IQAO,
     &                  FMO,LRMO,LCMO,NZMO,IQMO,
     &                  TM1,LR1,LC1,NZTM1,IQTM1,
     &                  TM2,LR2,LC2,NZTM2,IQTM2,
     &                  WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     This routine performs unitary transformations indicated by TYP:
C
C       TYP = AOMO : AO-to-MO-transformation  FMO = (C+)FAOC
C       TYP = MOAO : MO-to-AO-transformation  FAO = CFMO(C+)
C
C     ADDF gives the possibility of adding the results
C
C     TREV indicates symmetry of F-matrix under timereversal
C
C       TREV = 'S' - symmetric     (t =  1)
C       TREV = 'A' - anti-symmtric (t = -1)
C
C     The unitary C-matrix has a time-symmetric structure.
C
C     The AO-to-MO-transformation can be set up quaternionically as
C
C       FMO = [(Ca+)-t(CbT)j][(FAOa)+(FAOb)j][(Ca)+(Cb)j]
C
C       F(NBAS,NBAS,NZ) --->  F(NORB,NORB,NZ)
C
C     and the MO-to-Ao-transformation can be set up as:
C
C       FAO = [(Ca)+t(Cb)j][(FMOa)+(FMOb)j][(Ca+)-(Cb)j]
C
C       F(NORB,NORB,NZ) --->  F(NBAS,NBAS,NZ)
C
C     Written by T.Saue,January 1995
C     Last revision: 29 January 1997 : Luuk Visscher
C****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
      CHARACTER TREV*1,TRA*1,TYP*4
      DIMENSION FAO(LRAO,LCAO,*),FMO(LRMO,LCMO,*),
     &          TM1(LR1,LC1,*),TM2(LR2,LC2,*),
     &          IQAO(*),IQMO(*),IQTM1(*),IQTM2(*),
     &          WORK(LWORK)
      DIMENSION IQW(4)
C
      CALL QENTER('QTRANS')
#include <memint.h>
C
C     Add the results on to the contents of the MO matrix ?
C
c     print*, 'The AO matrix'
c     CALL PRQMAT(FAO,NRAO,NCAO,LRAO,LCAO,NZAO,IQAO,LUPRI)
c     print*, 'The transformation matrix TM1'
c     CALL PRQMAT(TM1,NRAO,NCMO,LR1,LC1,NZTM1,IQTM1,LUPRI)
c     print*, 'The transformation matrix TM2'
c     CALL PRQMAT(TM2,NRAO,NCMO,LR2,LC2,NZTM2,IQTM2,LUPRI)
C
C
C     Time reversal symmetry of F-matrix
C
      IF    (TREV.EQ.'S') THEN
        TRA = 'N'
      ELSEIF(TREV.EQ.'A') THEN
        TRA = 'I'
      ELSE
        CALL QUIT('QTRANS: Unknown keyword TREV '//TREV)
      ENDIF
C
C     *******************************************************
C     *****  AO - to - MO -transformation : F' = (C+)FC *****
C     *******************************************************
C
      IF(TYP.EQ.'AOMO') THEN
C
C       Look for the most efficient transformation
C       ------------------------------------------
C
        IF(NCAO.GE .NRAO) THEN
C
C         Perform the two steps:
C           1. W  = FC
C           2. F' = (C+)W
C
          CALL IQPACK(IQAO,NZAO,IQTM2,NZTM2,IQW,NZW)
          LWRK = NRAO*NCMO*NZW
          CALL MEMGET('REAL',KWRK,LWRK,WORK,KFREE,LFREE)
C
C         First part of transformation: W = FC
C         ------------------------------------
C
          CALL QGEMM(NRAO,NCMO,NCAO,D1,
     &             'N','N',IQAO,FAO,LRAO,LCAO,NZAO,
     &             'N','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &        D0,IQW,WORK(KWRK),NRAO,NCMO,NZW)
C
C         Second part of transformation: H = (C+)W
C         ----------------------------------------
C
          CALL IQPACK(IQW,NZW,IQTM1,NZTM1,IQMO,NZMO)
C          IF(NZBF.NE.NZMO) CALL QUIT('QTRANS: NZBF.NE.NZMO!')
          CALL QGEMM(NRMO,NCMO,NRAO,D1,
     &       'H',TRA,IQTM1,TM1 ,LR1,LC1,NZTM1,
     &       'N','N',IQW ,WORK(KWRK),NRAO,NCMO,NZW,
     &       FADD,IQMO,FMO ,LRMO,LCMO,NZMO)
          CALL MEMREL('QTRANS.AOMO',WORK,KWORK,KWORK,KFREE,LFREE)
        ELSE
C
C         Perform the two steps:
C           1. W  = (C+)F
C           2. F' = WC
C
          CALL IQPACK(IQAO,NZAO,IQTM1,NZTM1,IQW,NZW)
          LWRK = NRMO*NCAO*NZW
          CALL MEMGET('REAL',KWRK,LWRK,WORK,KFREE,LFREE)
C
C         First part of transformation: W = (C+)F
C         ----------------------------------------
C
          CALL QGEMM(NRMO,NCAO,NRAO,D1,
     &       'H',TRA,IQTM1,TM1,LR1,LC1,NZTM1,
     &       'N','N',IQAO,FAO,LRAO,LCAO,NZAO,
     &        D0,IQW,WORK(KWRK),NRMO,NCAO,NZW)
C
C         Second part of transformation: H = WC
C         -------------------------------------
C
          CALL IQPACK(IQW,NZW,IQTM2,NZTM2,IQMO,NZMO)
C          IF(NZBF.NE.NZMO) CALL QUIT('QTRANS: NZBF.NE.NZMO!')
          CALL QGEMM(NRMO,NCMO,NCAO,D1,
     &       'N','N',IQW ,WORK(KWRK),NRMO,NCAO,NZW,
     &       'N','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       FADD,IQMO,FMO ,LRMO,LCMO,NZMO)
          CALL MEMREL('QTRANS.AOMO',WORK,KWORK,KWORK,KFREE,LFREE)
C
        ENDIF
C        IF(IPRINT.GE.5) THEN
C          CALL HEADER('QTRANS: AO-to-MO-transformed Fock matrix',-1)
C          CALL PRQMAT(FMO,NRMO,NCMO,LRMO,LCMO,NZMO,IPQTOQ(1,0,1),LUPRI)
C        ENDIF
C
C     *******************************************************
C     *****  MO - to - AO -transformation : F' = CF(C+) *****
C     *******************************************************
C
      ELSEIF(TYP.EQ.'MOAO') THEN
C
C       Look for the most efficient transformation
C       ------------------------------------------
C
        IF(NCMO.GE.NRMO) THEN
          CALL IQPACK(IQMO,NZMO,IQTM2,NZTM2,IQW,NZW)
          LWRK = NRMO*NCAO*NZW
          CALL MEMGET('REAL',KWRK,LWRK,WORK,KFREE,LFREE)
C
C         First part of transformation: W = F(C+)
C
          CALL QGEMM(NRMO,NCAO,NCMO,D1,
     &       'N','N',IQMO,FMO,LRMO,LCMO,NZMO,
     &       'H','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       D0,IQW,WORK(KWRK),NRMO,NCAO,NZW)
C
C         Second part of transformation: H = CW
C
          CALL IQPACK(IQW,NZW,IQTM1,NZTM1,IQAO,NZAO)
C          IF(NZBF.NE.NZAO) CALL QUIT('QTRANS: NZBF.NE.NZAO!')
          CALL QGEMM(NRAO,NCAO,NRMO,D1,
     &       'N',TRA,IQTM1,TM1 ,LR1,LC1,NZTM1,
     &       'N','N',IQW ,WORK(KWRK),NRMO,NCAO,NZW,
     &       FADD,IQAO,FAO,LRAO,LCAO,NZAO)
          CALL MEMREL('QTRANS.MOAO',WORK,KWORK,KWORK,KFREE,LFREE)
        ELSE
          CALL IQPACK(IQMO,NZMO,IQTM1,NZTM1,IQW,NZW)
          LWRK = NRAO*NCMO*NZW
          CALL MEMGET('REAL',KWRK,LWRK,WORK,KFREE,LFREE)
C
C         First part of transformation: W = CF
C
          CALL QGEMM(NRAO,NCMO,NRMO,D1,
     &       'N',TRA,IQTM1,TM1 ,LR1,LC1,NZTM1,
     &       'N','N',IQMO,FMO,LRMO,LCMO,NZMO,
     &       D0,IQW,WORK(KWRK),NRAO,NCMO,NZW)
C
C         Second part of transformation: H = CW
C
          CALL IQPACK(IQW,NZW,IQTM2,NZTM2,IQAO,NZAO)
C          IF(NZBF.NE.NZAO) CALL QUIT('QTRANS: NZBF.NE.NZAO!')
          CALL QGEMM(NRAO,NCAO,NCMO,D1,
     &       'N','N',IQW ,WORK(KWRK),NRAO,NCMO,NZW,
     &       'H','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       FADD,IQAO,FAO,LRAO,LCAO,NZAO)
          CALL MEMREL('QTRANS.MOAO',WORK,KWORK,KWORK,KFREE,LFREE)
        ENDIF
C       IF(IPRINT.GE.5) THEN
C          CALL HEADER('QTRANS: MO-to-AO-transformed Fock matrix',-1)
C          CALL PRQMAT(FAO,NRAO,NCAO,LRAO,LCAO,NZAO,IPQTOQ(1,0,1),LUPRI)
C        ENDIF
      ELSE
        CALL QUIT('FTRANS: Unknown keyword TYP '//TYP)
      ENDIF
C
      CALL QEXIT('QTRANS')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bcktra */
      SUBROUTINE BCKTRA(COFT,LRT,LCT,COFU,LRU,LCU,NORB,NZ,IQC,
     &                  NVEC,JVEC,NBAS,TMAT,LR1,LC1,NZT,IQT,IPRINT)
C*****************************************************************************
C
C     PURPOSE:
C
C     METHOD:
C
C     INPUT:
C	NORB - number of vectors to be tranformed
C	COFU - array of untransformed vectors of dimension (NORB,NORB,NZ)
C	COFT - array of backtransformed vectors of dimension (NBAS,NVEC,NZ)
C       TMAT - transformation matrix of dimension (NBAS,NORB,NZT)
C
C	Written by T.Saue January 1995
C	LAST VERSION: June 27 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0)
C
      DIMENSION COFU(LRU,LCU,NZ),COFT(LRT,LCT,NZ),
     &          TMAT(LR1,LC1,NZT),IQC(NZ),IQT(NZT)
C
      IF(NVEC.EQ.0) RETURN
      IF(IPRINT.GE.5) THEN
        CALL HEADER('BCKTRA:Nontransformed coefficients',-1)
        WRITE(LUPRI,'(2(3X,A,I5))')
     +    'Start index: ',JVEC,'No. of vectors:',NVEC
        CALL PRQMAT(COFU,NORB,NORB,LRU,LCU,NZ,IQC,LUPRI)
        CALL HEADER('BCKTRA:Transformation matrix:',-1)
        CALL PRQMAT(TMAT,NBAS,NORB,LR1,LC1,NZT,IQC,LUPRI)
      ENDIF
      CALL QGEMM(NBAS,NVEC,NORB,D1,
     &           'N','N',IQT,TMAT,LR1,LC1,NZT,
     &           'N','N',IQC,COFU(1,JVEC,1),LRU,LCU,NZ,
     &                  D0,IQC,COFT,LRT,LCT,NZ)
      IF(IPRINT.GE.5) THEN
        CALL HEADER('BCKTRA:Backtransformed coefficients',-1)
        WRITE(LUPRI,'(2(3X,A,I5))')
     +    'Start index: ',JVEC,'No. of vectors:',NVEC
        CALL PRQMAT(COFT,NBAS,NVEC,LRT,LCT,NZ,IQC,LUPRI)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck delsph */
      SUBROUTINE DELSPH(AMAT,LCOMP,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Delete totally symmetric contributions of spherical transformation
C     from LL or SS-block of matrix
C
C     LCOMP  = 1   LL-block
C     LCOMP  = 2   SS-block
C
C     Written June 20 1996 - Trond Saue
C     Last revision June 20 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <mxcent.h>
C
#include <nuclei.h>
C
      DIMENSION AMAT(*),WORK(LWORK)
C
      CALL QENTER('DELSPH')
#include <memint.h>
      IF    (LCOMP.EQ.1) THEN
C     LL-block
        NDIM = NLARGE
        IOFF = 1
      ELSEIF(LCOMP.EQ.2) THEN
C     SS-block
        NDIM = NSMALL
        IOFF = NBASIS*NLARGE+NLARGE+1
      ELSE
        WRITE(LUPRI,'(A,I5)') 
     &    'DELSPH: Unknown option ',LCOMP
        CALL QUIT('DELSPH: Unknown option LCOMP !')
      ENDIF
C
      NSPH = NSPHCM(LCOMP)
      NWRK = NDIM*NSPH
      WRITE(LUPRI,'(A,I6)') '* DELSPH: NSPH',NSPH
      CALL OUTPUT(AMAT,1,NLARGE,1,NLARGE,NLARGE,NLARGE,1,LUPRI)
C     Allocate buffers
      CALL MEMGET('REAL',KW1,NWRK,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KW2,NWRK,WORK,KFREE,LFREE)
C     First transform rows
      CALL VTAOCS(AMAT(IOFF),NBASIS,WORK(KW1),NSPH,
     &              NDIM,LCOMP,IPRINT)
C     ...then transpose
      CALL MTRSP(NSPH,NDIM,WORK(KW1),NSPH,WORK(KW2),NDIM)
C     ...and transform columns
      CALL VTAOCS(WORK(KW2),NDIM,WORK(KW1),NSPH,
     &              NSPH,LCOMP,IPRINT)
C     Backtransform columns
      CALL VTCSAO(WORK(KW1),NSPH,WORK(KW2),NDIM,
     &              NSPH,LCOMP,IPRINT)
C     ...then backtranspose
      CALL MTRSP(NDIM,NSPH,WORK(KW2),NDIM,WORK(KW1),NSPH)
C     ...and finally backtransform rows
      CALL VTCSAO(WORK(KW1),NSPH,AMAT(IOFF),NBASIS,
     &              NDIM,LCOMP,IPRINT)
C
C
C     Print section
C     =============
C
C      IF(IPRINT.GE.5) THEN
C        CALL HEADER('Output from DELSPH',-1)
C        CALL PRQMAT(AMAT(IOFF),NBASIS,NBASIS,NBASIS,NBASIS,1,
C     &              IPQTOQ(1,0,1),LUPRI)
C      ENDIF
      CALL QEXIT('DELSPH')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mksphc */
      SUBROUTINE MKSPHC(WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     This routine generates an orthogonal transformation from Cartesian
C     to solid harmonics
C
C     Written by T.Saue Aug 15 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <maxaqn.h>
C
      DIMENSION WORK(LWORK)
      CALL QENTER('MKSPHC')
#include <memint.h>
      NIJK = MXQN*MXQN*MXQN
      MXYZ = (MXQN*(MXQN+1))/2      
      NFUN = (MXQN*(MXQN+1)*(MXQN+2))/6
      CALL MEMGET('INTE',KINADR,3*NFUN,WORK,KFREE,LFREE)            
      CALL MEMGET('INTE',KINIJK,  NIJK,WORK,KFREE,LFREE)            
      CALL MEMGET('INTE',KIX,     MXYZ,WORK,KFREE,LFREE)            
      CALL MEMGET('INTE',KJY,     MXYZ,WORK,KFREE,LFREE)            
      CALL MEMGET('INTE',KKZ,     MXYZ,WORK,KFREE,LFREE)            
      CALL MEMGET('INTE',KLOFF,   MXQN,WORK,KFREE,LFREE)            
      CALL MEMGET('REAL',KDFAC,   MXQN,WORK,KFREE,LFREE)            
      CALL MKSPH1(WORK(KINADR),WORK(KINIJK),
     &            WORK(KIX),WORK(KJY),WORK(KKZ),WORK(KLOFF),
     &            WORK(KDFAC),
     &            WORK(KFREE),LFREE,IPRINT)
      CALL QEXIT('MKSPHC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mksph1 */
      SUBROUTINE MKSPH1(INADR,INIJK,IX,JY,KZ,LOFF,DFAC,
     &                  WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     This routine generates an orthogonal transformation from Cartesian
C     to solid harmonics
C
C     Written by T.Saue Aug 15 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0)
C
#include <maxaqn.h>
C
#include <sphtrm.h>
      DIMENSION INADR(3,*),INIJK(0:MXQN-1,0:MXQN-1,0:MXQN-1),
     &          IX(*),JY(*),KZ(*),LOFF(*),DFAC(*),WORK(*)
C
      CALL QENTER('MKSPH1')
C
      IND  = 0
      IOFF = 1
      FAC  = D1
      CALL DZERO(CSP,NCSP)
      DO LVAL = 0,MXQN-1
        NLM            = 2*LVAL + 1
        NXYZ           = (LVAL+1)*(LVAL+2)/2
        N2XYZ          = NXYZ*NXYZ
        ISPADR(LVAL+1) = IOFF
        IOFF           = IOFF + N2XYZ
        LOFF(LVAL+1)   = IND
        FAC            = FAC*NLM
        DFAC(LVAL+1)   = FAC
        CALL LMNVAL(LVAL+1,NXYZ,IX,JY,KZ)
        DO K = 1,NXYZ
          IND = IND + 1
          INADR(1,IND) = IX(K)
          INADR(2,IND) = JY(K)
          INADR(3,IND) = KZ(K)
        ENDDO
        IF (LVAL .EQ. 0) THEN
          CSP(1)       = D1
        ELSE IF (LVAL .EQ. 1) THEN
          CSP (2)      = D1
          CSP (6)      = D1
          CSP(10)      = D1
        ELSE
          CALL MKSPH2(LVAL,INADR,INIJK,LOFF,DFAC,WORK,LWORK,IPRINT)
        ENDIF
      ENDDO
C
      CALL QEXIT('MKSPH1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mksph2 */
      SUBROUTINE MKSPH2(LVAL,INADR,INIJK,LOFF,DFAC,
     &                  WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     This routine generates an orthogonal transformation from Cartesian
C     to solid harmonics for a given l-value.
C
C     Written by T.Saue Aug 15 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0)
C
#include <maxaqn.h>
C
      DIMENSION INADR(3,*),INIJK(0:MXQN-1,0:MXQN-1,0:MXQN-1),
     &          LOFF(*),DFAC(*),WORK(*)
C
      CALL QENTER('MKSPH2')
#include <memint.h>
      NXYZ           = (LVAL+1)*(LVAL+2)/2
      NNXYZ          = (NXYZ*(NXYZ+1))/2
      N2XYZ          = NXYZ*NXYZ
      CALL MEMGET('REAL',KSXYZ,NNXYZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOEF,N2XYZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF ,NXYZ ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,NXYZ ,WORK,KFREE,LFREE)
      CALL MKSPH3(LVAL,INADR,INIJK,LOFF,DFAC,WORK(KSXYZ),
     &            WORK(KCOEF),WORK(KBUF),WORK(KIBUF),
     &                  WORK(KFREE),LFREE,IPRINT)
      CALL QEXIT('MKSPH2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mksph3 */
      SUBROUTINE MKSPH3(LVAL,INADR,INIJK,LOFF,DFAC,
     &                  SXYZ,COEF,BUF,IBUF,
     &                  WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     This routine generates an orthogonal transformation from Cartesian
C     to solid harmonics for a given l-value.
C
C     Written by T.Saue Aug 15 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0)
C
#include <maxaqn.h>
C
#include <ccom.h>
      DIMENSION INADR(3,*),INIJK(0:MXQN-1,0:MXQN-1,0:MXQN-1),
     &          LOFF(*),DFAC(*),SXYZ(*),COEF(*),BUF(*),
     &          IBUF(*),WORK(*)
C
      NXYZ = (LVAL+1)*(LVAL+2)/2
C
C     Make overlap matrix in minimal basis
C
      CALL MKSXYZ(LVAL,SXYZ,INADR,LOFF,DFAC,IPRINT)      
C
C     Diagonalize minimal overlap matrix
C
      CALL DUNIT(COEF,NXYZ)
      CALL JACO(SXYZ,COEF,NXYZ,NXYZ,NXYZ,BUF,IBUF)
      IF(IPRINT.GE.4) THEN
        WRITE(LUPRI,'(A)') 
     &    '* Eigenvalues of minimal overlap matrix'
        II = 0
        DO I = 1,NXYZ
          II = (I*(I+1))/2
          WRITE(LUPRI,'(A,I5,F10.5)') '*EIG :',II,SXYZ(II)
        ENDDO
        WRITE(LUPRI,'(A)') 
     &    '* Eigenvectors of minimal overlap matrix'
        CALL OUTPUT(COEF,1,NXYZ,1,NXYZ,NXYZ,NXYZ,1,LUPRI)
      
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mksxyz */
      SUBROUTINE MKSXYZ(LVAL,SXYZ,INADR,LOFF,DFAC,IPRINT)
C*****************************************************************************
C
C     Make lower triangular row-packed overlap matrix in minimal 
C     Cartesian Gaussian basis for angular momentum value L
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0=0.0D0,D1 = 1.0D0)
C
      DIMENSION INADR(3,*),LOFF(*),DFAC(*),SXYZ(*)
C
      NXYZ = (LVAL+1)*(LVAL+2)/2
      IOFF = LOFF(LVAL+1)
      IND  = 0
      DO I = 1,NXYZ
        DO J = 1,I
          IND = IND + 1
          SBUF = D1
          DO K = 1,3
            II = INADR(K,IOFF+I) + INADR(K,IOFF+J)
            IF(MOD(II,2).EQ.1) THEN
              SBUF = D0
            ELSEIF(II.GT.0) THEN
              II = II/2
              SBUF = SBUF*DFAC(II)
            ENDIF
          ENDDO
          SXYZ(IND) = SBUF
        ENDDO
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.4) THEN
        CALL HEADER('Output from MKSXYZ',-1)
        WRITE(LUPRI,'(A,I5)') '* L -value:',LVAL
        WRITE(LUPRI,'(A)') '* Minimal overlap matrix'
        CALL OUTPAK(SXYZ,NXYZ,NXYZ,LUPRI)
      ENDIF      
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Getsph */
      SUBROUTINE GETSPH(CSPH,ICSP,LMAX,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Get spherical transformation matrices up to LMAX
C     Written by T.Saue May 1 1998
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      PARAMETER(D1 = 1.0D0)
#include <ccom.h>
#include <sphtrm.h>
#include <symmet.h>
#include <pgroup.h>
      DIMENSION CSPH(*),ICSP(*),WORK(*)
C
C     MORDER eq 0 : M_l order 0, +1, -1, +2, -2, ..., +LVAL, -LVAL
C     MORDER ne 0 : M_l order -LVAL, ..., -1, 0, 1, ..., +LVAL
C
C     MINTEG eq 1 : the spherical components are normalized such
C                   that all coefficients are integers.
      MORDER = 1
      MINTEG = 1
      IOFF   = 0
      DO LL = 1, LMAX
        L    = LL - 1
        NLM  = 2*L + 1
        NXYZ = (LL*(LL+1))/2
        CALL SPHCOM(L,CSPH(IOFF+1),NLM,NXYZ,MORDER,MINTEG,
     &                WORK,LWORK,IPRINT)
        ICSP(LL) = IOFF
        IOFF     = IOFF + NLM*NXYZ
      ENDDO
      RETURN
C     
      END
