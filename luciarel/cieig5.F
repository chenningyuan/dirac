      SUBROUTINE CIEIG5(MV7,INICI,EROOT,VEC1,
     &           VEC2,MINST,LUDIA,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &           NDIM,NBLK,NROOT,MAXVEC,MXCIIT,MXROOT,LUINCI,
     &           IPRT,PEIGVC,NPRDET,H0,IPNTR,NP1,NP2,NQ,H0SCR,
     &           SCR2,T_CC,
     &           EIGSHF,ICISTR,LBLK,IDIAG,VEC3,THRES_E,
     &           NBATCH,LBATCHB,LBATCHE,
CTF  &           LBLOCK,
     &           IBLOCK,INIDEG,E_THRE,C_THRE,E_CONV,C_CONV,ICLSSEL,
CTF  &           IBLK_TO_CLS,
     &           NCLS,
CTF  &           CLS_C,CLS_E,CLS_CT,CLS_ET,CLS_A,ICLS_L,BLKS_A,
CTF  &           CLS_DEL,CLS_DELT,ISKIPEI,I2BLK,ICLS_A2,
     &           MXLNG,IRC,IDCOMH,IHAM12,SIGDEN_ROUTE)
*
* Master routine for CI diagonalization
*
* Modified to handle PQ - preconditioner , May 1990
* PICO,MICDV4 added spring of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION INIDET(100)
      PARAMETER( LLWRK = 70000 )
      COMMON/SCR/SCR1(LLWRK),ISCR1(LLWRK)
*. Output from Subspace dagonalization
      DIMENSION H0(*),IPNTR(*),H0SCR(*),PEIGVC(*),SCR2(*)
*.
      DIMENSION EROOT(NROOT)
      IF( IPRT.GT. 1 )  WRITE(6,'(/A)')
     &'          *** information from ci diagonalization  ***'
      WRITE(6,*)
*
      NTEST=0
*
*               ====================================
** 1 :               INITIAL VARIATIONAL SUBSPACE
*               ====================================
*
      IF( INICI .EQ. 0 ) THEN
        IF(NPRDET .EQ. 0 ) THEN
* ==================================================
*. Initial guess from lowest elements of CI diagonal
* ==================================================
* In order to treat degeneracies, the lowest 4 * NROOT elements are
*.obtained

          write(6,*) ' CIEIG5 NDIM NROOT ', NDIM,NROOT
          NFINDM = MIN(NDIM,4*NROOT)
          CALL FNDMND(LUDIA,LBLK,VEC1,NFINDM,NFINDA,ISCR1(1+2*NFINDM),
     &                SCR1(1+2*NFINDM),ISCR1,SCR1,IPRT)
C              FNDMND(LU,LBLK,SEGMNT,NSUBMX,NSUB,ISCR,SCR,ISCAT,SUBVAL
          CALL REWINO(LU1)
          IBASE = 1
          TEST = 1.0D-10
*
          DO 100 IROOT = 1, NROOT
*. Number of degenerate elements
            NDEG = 1
            XVAL = SCR1(IBASE)
   90       CONTINUE
            IF(IBASE-1+NDEG+1.LE.NFINDA) THEN
              IF (ABS(SCR1(IBASE-1+NDEG+1)-XVAL).LE.TEST) THEN
                NDEG = NDEG + 1
                GOTO 90
              END IF
            END IF
            WRITE(6,*) ' IROOT NDEG ', IROOT,NDEG
*
            IF(INIDEG.EQ.0.AND.NDEG.GT.1) THEN
              WRITE(6,*) ' WARNING WARNING WARNING WARNING ! '
              WRITE(6,*) ' DEGENERATE INITIAL VECTORS FOR CI '
              WRITE(6,*) ' I AM NOT ALLOWED TO TAKE THIS INTO '
              WRITE(6,*) ' CONSIDERATION SINCE  INIDEG = 0 '
              WRITE(6,*)
              WRITE(6,*) ' I hope you know what you are doing '
              NDEG = 1
            END IF

*. Initial guess in compressed form in SCR1
            SCALE = 1.0D0/SQRT(DFLOAT(NDEG))
            DO 250 II = 1,NDEG
*. Anti symmetric combination
              IF(INIDEG.EQ.-1) THEN
                SCR1(II) = (-1.0D0)**II * SCALE
*. Symmetric combination
              ELSE IF (INIDEG.EQ.1.OR.INIDEG.EQ.0) THEN
                SCR1(II) =  SCALE
              END IF
  250       CONTINUE
C                WRSVCD(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM)
* Real part
            CALL REWINE(LUDIA,LBLK)
c           ISCR1(1) = 16
            CALL WRSVCD(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,NDEG,NDIM,
     &                  LUDIA)
*. Imaginary part
            CALL WRSVCD(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,0,NDIM,
     &                  LUDIA)
            IBASE = IBASE + NDEG
  100     CONTINUE
        ELSE
* =====================================
*. Initial approximations are in PEIGVC
* =====================================
CTF  Currently never used !!!
          CALL REWINO(LU1)
          DO 1984 IROOT = 1, NROOT
C                WRSVCD(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM)
            CALL WRSVCD(LU1,LBLK,VEC1,IPNTR,
     &           PEIGVC((IROOT-1)*NPRDET+1),NPRDET,NDIM,LUDIA)
 1984     CONTINUE
        END IF
      END IF
*
* ======================================================
* Initial CI vectors are already on file LU1, do nothing
* ======================================================
*
      IF( INICI .LT. 0 ) THEN
*. Vectors assumed already in LU1
        IF(IPRT.GT. 1 )
     &  WRITE(6,*)' Initial CI vector assumed in place '
      END IF
*
*                 ========================
* 2 :                  Diagonalization
*                 ========================
*
      CALL QENTER('CIEIG5')
* Inverse iteration modified Davidson with 2 vectors in core
      IF(IPRT .GE. 2 ) THEN
         WRITE(6,*)
         WRITE(6,'(A,I3)')
     &   '  Number of roots to be converged..  ',NROOT
         WRITE(6,*)
         WRITE(6,'(A,I3)')
     &   '  Largest allowed number of vectors..',MAXVEC
         WRITE(6,*)
         WRITE(6,'(A,I3)')
     &   '  Allowed number of CI iterations  ..',MXCIIT
      END IF
*
      KRNRM = 1
      IF(IRC.EQ.1) THEN
       IFAC = 1
      ELSE
       IFAC = 2
      END IF
      KEIG = KRNRM + (MXCIIT+1)*NROOT
      KFIN = KEIG  + (MXCIIT+1)*2*MAXVEC
      KAPROJ = KFIN + NROOT
      IF(IRC.EQ.2) THEN
        KAPROJI = KAPROJ +  MAXVEC*(MAXVEC+1)/2
      ELSE
        KAPROJI = KAPROJ
      END IF
      KAVEC = KAPROJI + MAXVEC*(MAXVEC+1)/2
      KWORK = KAVEC + MAXVEC ** 2*IFAC*IFAC
      KLFREE = KWORK + MAXVEC*(MAXVEC+1)*IFAC*IFAC
      IF( IPRT .GE. 5 ) THEN
         WRITE(6,*) ' KRNRM KEIG KFIN KAPROJ KAVEC KWORK KLFREE '
         WRITE(6,'(6I4)')  KRNRM,KEIG,KFIN,KAPROJ,KAVEC,KWORK,KLFREE
      END IF
      IF( KLFREE-1 .GT. LLWRK) THEN
         WRITE(6,'(A,2I7)' )
     &   ' Not enough memory in CIEIG5 : needed and available ',
     &   KLFREE-1, LLWRK
         Call Abend2( ' insufficient memory in cieig5 ' )
      END IF
*
      IF(IRC.EQ.1) THEN
*. Real numbers only
       IF(IDIAG.EQ.1.AND.ICISTR.EQ.1) THEN
*. Routine using two complete vectors in core
         IOLSEN = 1
         IPICO = 0
         IF(MXCIIT.NE.0) THEN
         CALL MINDV4(MV7,VEC1,VEC2,LU1,LU2,SCR1(KRNRM),SCR1(KEIG),
     &              EROOT     ,MXCIIT,NDIM,LU3,LUDIA,NROOT,
     &              MAXVEC,NROOT,SCR1(KAPROJ),SCR1(KAVEC),
     &              SCR1(KWORK) ,IPRT,
     &              NPRDET,H0,IPNTR,NP1,NP2,NQ,H0SCR,EIGSHF,
     &              IOLSEN,IPICO)
         ELSE
*. No iterations, set energy to 0
           DO IROOT = 1, NROOT
             EROOT(IROOT) = 0.0D0
           END DO
         END IF
       ELSE IF(IDIAG.EQ.1.AND.ICISTR.GE.2) THEN
*. Routine using two blocks in core
*
        I_USE_MICDV6 = 1
         IF(I_USE_MICDV6.EQ.1) THEN
C       IF(.NOT.(MAXVEC.EQ.2.AND.NROOT.EQ.1)) THEN
         WRITE(6,*) ' MICDV6 will be called '
         CALL MICDV6(VEC1,VEC2,LU1,LU2,SCR1(KRNRM),SCR1(KEIG),
     &              EROOT     ,MXCIIT,NDIM,LU3,LU4,LU5,LU6,LU7,
     &              LUDIA,NROOT,
     &              MAXVEC,NROOT,SCR1(KAPROJ),SCR1(KAVEC),
     &              SCR1(KWORK) ,IPRT,
     &              NPRDET,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &              THRES_E)
        ELSE
*. Special version for one root 2 vectors
          WRITE(6,*) ' MICDV5 will be called '
         CALL MICDV5(VEC1,VEC2,LU1,LU2,SCR1(KRNRM),SCR1(KEIG),
     &              EROOT     ,MXCIIT,NDIM,LU3,LU4,LU5,LUDIA,NROOT,
     &              MAXVEC,NROOT,SCR1(KAPROJ),SCR1(KAVEC),
     &              SCR1(KWORK) ,IPRT,
     &              NPRDET,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &              THRES_E)
        END IF
       ELSE IF (IDIAG.EQ.2)  THEN
*.Routine using two vector segments and three files
CTF      NSUB = 0
CTF      WRITE(6,*) ' New PICO3'
CTF      WRITE(6,*) ' New PICO3'
CTF      WRITE(6,*) ' New PICO3'
CTF      WRITE(6,*) ' New PICO3'
CTF      WRITE(6,*) ' New PICO3'
CTF      WRITE(6,*) ' New PICO3'
CTF      WRITE(6,*) ' New PICO3'
CTF      WRITE(6,*) ' New PICO3'
CTF      WRITE(6,*) ' New PICO3'
CTF      WRITE(6,*) ' New PICO3'
CTF      WRITE(6,*) ' New PICO3'
CTF      WRITE(6,*) ' ISKIPEI before call to PICO2',ISKIPEI
CTF
CTF      CALL PICO3(VEC1,VEC2,LU1,LU2,LU3,LU4,SCR1(KRNRM),SCR1(KEIG),
CTF  &        EROOT,MXCIIT,NBATCH,LBATCHB,LBATCHE,LBLOCK,IBLOCK,
CTF  &        IPRT,
CTF  &        NPRDET,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
CTF  &        THRES_E,E_THRE,C_THRE,E_CONV,C_CONV,ICLSSEL,
CTF  &        IBLK_TO_CLS,NCLS,CLS_C,CLS_E,CLS_CT,CLS_ET,
CTF  &        CLS_A,ICLS_L,BLKS_A,CLS_DEL,CLS_DELT,ISKIPEI,I2BLK,VEC3,
CTF  &        ICLS_A2,MXLNG)
CTF
         continue
       END IF
*
      ELSE
*. Complex numbers may be encountered
*
        CALL CMICDV(MV7,VEC1,VEC2,LU1,LU2,SCR1(KRNRM),SCR1(KEIG),
     &     SCR2(1),MXCIIT,
     &     NDIM,LU3,LU4,LU5,LU6,LUDIA,NROOT,MAXVEC,NROOT,
     &     T_CC,
     &     SCR1(KAPROJI),SCR1(KAPROJ),SCR1(KAVEC),SCR1(KWORK),IPRT,
     &     NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &     0,0,IDCOMH,IHAM12,SIGDEN_ROUTE)
*
C      SUBROUTINE CMICDV(MV7,VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,
C     &                  NVAR,LU3,LU4,LU5,LU6,LUDIA,NROOT,MAXVEC,NINVEC,
C     &                  T_CC,
C     &                  APROJI,APROJR,AVEC,WORK,IPRT,
C     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
C     &                  IPREDZN,INONREL,IDCOMH,IHAM12,SIGDEN_ROUTE)
C
      END IF
*     ^ End if switch real/complex numbers
      ENOT = SCR1(KFIN)
*
      call qexit('CIEIG5')
      RETURN
      END
