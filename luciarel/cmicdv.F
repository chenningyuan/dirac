      SUBROUTINE CMICDV(MV7,VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,
     &                  NVAR,LU3,LU4,LU5,LU6,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  T_CC,
     &                  APROJI,APROJR,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  IPREDZN,INONREL,IDCOMH,IHAM12,SIGDEN_ROUTE)
*
* Multi root Davidson routine for complex Hermitian
* matrices , requires two blocks in core
*
* May 1994, Jeppe Olsen
*
* Initial version without general preconditioner
* or inverse iteration corrections
*
*
* Input :
* =======
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU3,LU4   : Scatch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
*
* The vectors are in general multi record vectors.
* The real and imaginary parts are separated by end of vector marks
* (-1) in my simple world
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
      REAL * 8   INPROD,INPRDD
      DIMENSION RNRM(NROOT,MAXIT),EIG(2*MAXVEC,MAXIT)
      DIMENSION APROJI(*),APROJR(*),AVEC(*),WORK(*)
      DIMENSION H0(*),IPNTR(1)
      DIMENSION H0SCR(*)
*
* Dimensioning required of local vectors
*      APROJR : MAXVEC*(MAXVEC+1)/2
*      APROJI : MAXVEC*(MAXVEC+1)/2
*      AVEC   : (2*MAXVEC) ** 2
*      WORK   : 2*MAXVEC*(2*MAXVEC+1)/2
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
      DIMENSION FINEIG(*)
      character*6 SIGDEN_ROUTE
      LOGICAL CONVER,RTCNV(100)
      LOGICAL EXSTOP
      EXTERNAL MV7
*
      call qenter('CMICDV')
*
      NTEST = 0000
      IPRT = MAX(NTEST,IPRT)
*
      write(6,*)
      write(6,*) '----------------------------------------------'
      write(6,*) ' I´ll take you to a place '
      write(6,*) ' Where we shall find our'
      write(6,*) '  ROOTS, BLOODY ROOTS '
      write(6,*)
      write(6,*) '               Sepultura '
      write(6,*) '----------------------------------------------'
      write(6,*)
      write(6,*) ' Print flag, Nuclear repulsion energy'
      write(6,*) ' IPRT, EIGSHF ',IPRT,EIGSHF
*
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(6,*) ' Sorry CMICDV wounded , MAXVEC .LT. 2*NROOT '
         WRITE(6,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(6,*) ' Raise MXCIV to be at least 2 * Nroot '
         WRITE(6,*) ' Enforced stop in CMICDV '
         Call Abend1( 20 )
       END IF
*
       CONVLIM = 1.0D-5
       CONVER = .FALSE.
*
       IF(MAXIT.EQ.0) THEN
         WRITE(6,*) ' MAXIT = 0, so just a fast visit '
         RETURN
       END IF
*
       if (IPRT.ge.1) then
         WRITE(6,*) ' Diagonal in CMICDV'
         WRITE(6,*) ' =================='
         CALL REWINE(LUDIA,LBLK)
         CALL WRTVCDC(VEC1,LUDIA,1,LBLK)
       end if
*
       if (IDCOMH.eq.0) then
          I_DO_COMHAM = 0
       else if (IDCOMH.eq.1) then
          I_DO_COMHAM = 1
       else
          write(6,*) 'Problemo in CMICDV.'
          Call Abend2( 'Quitting.' )
       end if
*
       IF(I_DO_COMHAM .EQ.1) THEN
         write(6,*)
         write(6,*) ' ===================================== '
         WRITE(6,*) ' Control will be transferred to COMHAM '
         write(6,*) ' ===================================== '
         WRITE(6,*) ' NVAR = ', NVAR
         CALL COMHAM(NVAR,T_CC,SIGDEN_ROUTE)
         Call Abend2( ' After COMHAM  ' )
       END IF

*
* ===================
*.Initial iteration
* ===================
       ITER = 1
       NVEC = 0
       IF( MAXIT .GE. 1 ) THEN
*
         IF(IPRT.GE.1000) THEN
           WRITE(6,*) ' Initial vectors '
           WRITE(6,*) ' ================'
           CALL REWINE(LU1,LBLK)
           DO IVEC = 1, NINVEC
             WRITE(6,*) ' Initial vector ...  ', IVEC
             CALL WRTVCDC(VEC1,LU1,0,LBLK)
           END DO
         END IF
*
         CALL REWINE(LU1,LBLK)
         CALL REWINE(LU2,LBLK)
         DO 10 IVEC = 1,NINVEC
           CALL REWINE(LU3,LBLK)
           CALL REWINE(LU4,LBLK)
*. Copy to scratch unit LU3
           CALL COPVCDC(LU1,LU3,VEC1,0,LBLK)
*
*--------------------------------------------------------------
* Generalization to SIGDEN routine block implemented.
*--------------------------------------------------------------
           if (SIGDEN_ROUTE.eq.'TRADIT') then
             call mv7(VEC1,VEC2,LU3,LU4,INONREL)
           else if (SIGDEN_ROUTE.eq.'SIGDEN') then
             call sigden_ctrl(VEC1,VEC2,LU3,LU4,T_CC,1)
           else
             write(6,*) 'Route for sigma and densities not allowed.'
             write(6,*) 'SIGDEN_ROUTE is ',SIGDEN_ROUTE
             call abend ('Quitting.')
           end if
*
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
           CALL REWINE(LU4,LBLK)
           CALL COPVCDC(LU4,LU2,VEC1,0,LBLK)
*. Projected matrix
           DO 8 JVEC = 1, IVEC
             IJ = IVEC*(IVEC-1)/2 + JVEC
             CALL INPRDDC(VEC1,VEC2,LU2,LU3,JVEC,1,1,LBLK,XR,XI)
C                 INPRDDC(VEC1,VEC2,LU1,LU2,IOFF1,IOFF2,IREW,LBLK,XR,XI)
             APROJR(IJ) = XR
             APROJI(IJ) = -XI
    8      CONTINUE
   10    CONTINUE
*
         IF(IPRT.GE.1000) THEN
           WRITE(6,*) ' Initial sigma vectors '
           WRITE(6,*) ' ======================'
           CALL REWINE(LU2,LBLK)
           DO IVEC = 1, NINVEC
             WRITE(6,*) ' Initial vector ...  ', IVEC
             CALL WRTVCDC(VEC1,LU2,0,LBLK)
           END DO
         END IF
*
         IF( IPRT .GE.10.OR.IPREDZN.NE.0 ) THEN
           WRITE(6,*) ' Initial projected matrix, real part'
           CALL PRSYM(APROJR,NINVEC)
           WRITE(6,*) ' Initial projected matrix, imag part'
           CALL PRSYM(APROJI,NINVEC)
          END IF
*. Diagonalize initial projected matrix
C             CEIGEN(HR,HI,NDIM,EIG,C,SCR
         print*,'NINVEC FOR CEIGEN ',NINVEC
         CALL CEIGEN(APROJR,APROJI,NINVEC,EIG(1,1),AVEC,WORK)
*
         IF(IPRT .GE. 3.OR.IPREDZN.NE.0 ) THEN
           WRITE(6,'(A,I4)') ' Eigenvalues of initial iteration '
           WRITE(6,'(F22.13)')
     &     ( EIG(IROOT,1)+EIGSHF,IROOT=1,NROOT)
         END IF
         NVEC = NINVEC
       END IF
*
       IF (MAXIT.EQ.1.OR.IPREDZN.EQ.1) GOTO  901
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
       write(6,*)
       WRITE(6,*) ' Info from iteration .... ', ITER
       WRITE(6,*) ' _________________________'
*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = H*X - EIGAPR*X
*
       IADD = 0
       CONVER = .TRUE.
       DO 100 IROOT = 1, NROOT
         EIGAPR = EIG(IROOT,ITER)
C        EIGAPR = EIG(IROOT,ITER-1)
*. ==========================
*. Real part of HX - EIGAPR*X
*. ==========================
         DO IVEC = 1, NVEC
           WORK(2*(IVEC-1)+1) = AVEC((IROOT-1)*2*NVEC+IVEC)
           WORK(2*(IVEC-1)+2) =-AVEC((IROOT-1)*2*NVEC+NVEC+IVEC)
         END DO
*. Real part, Hx
         CALL MVCSMD(LU2,WORK,LU3,LU6,VEC1,VEC2,2*NVEC,1,LBLK)
C             MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
*. Real part X
         CALL MVCSMD(LU1,WORK,LU5,LU6,VEC1,VEC2,2*NVEC,1,LBLK)
*
C?       WRITE(6,*) 
C?   &   ' Current approximation to real part of eigenvector'
C?       CALL WRTVCD(VEC1,LU5,1,LBLK)
*. Real part of Hx -Ex on LU4
         EIGAPRM = - EIGAPR
         ONE = 1.0D0
         CALL REWINE(LU3,LBLK)
         CALL REWINE(LU5,LBLK)
         CALL REWINE(LU4,LBLK)
         CALL VECSMD(VEC1,VEC2,ONE,EIGAPRM,LU3,LU5,LU4,0,LBLK)
*. ==========================
*. Imag part of HX - EIGAPR*X
*. ==========================
         DO IVEC = 1, NVEC
           WORK(2*(IVEC-1)+1) = AVEC((IROOT-1)*2*NVEC+NVEC+IVEC)
           WORK(2*(IVEC-1)+2) = AVEC((IROOT-1)*2*NVEC     +IVEC)
         END DO
*. Imag part of HX
         CALL MVCSMD(LU2,WORK,LU3,LU6,VEC1,VEC2,2*NVEC,1,LBLK)
C               MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
*. Imag part X
         CALL MVCSMD(LU1,WORK,LU5,LU6,VEC1,VEC2,2*NVEC,1,LBLK)
*. Imag part of Hx -Ex as second vector on LU4
         EIGAPRM = - EIGAPR
         ONE = 1.0D0
         CALL REWINE(LU3,LBLK)
         CALL REWINE(LU5,LBLK)
         CALL VECSMD(VEC1,VEC2,ONE,EIGAPRM,LU3,LU5,LU4,0,LBLK)
         IF ( IPRT  .GE. 600 ) THEN
           WRITE(6,*) '  ( HX - EX )'
           CALL WRTVCDC(VEC1,LU4,1,LBLK)
         END IF
*  Strange place to put convergence but ....
         CALL REWINE(LU4,LBLK)
         RNORM = SQRT( INPRDD(VEC1,VEC1,LU4,LU4,0,LBLK)
     &         +       INPRDD(VEC1,VEC1,LU4,LU4,0,LBLK) )
         RNRM(IROOT,ITER) = RNORM
         WRITE(6,'(A19,7X,I3,3X,1E18.13,3X,1F19.12)')
     &     ' Iter RNORM EIGAPR ', ITER,RNORM,EIGAPR+EIGSHF
*
* Flush output buffer
         CALL FLSHFO(6)
*
* Set logical for convergence results
         IF(RNORM.LT. CONVLIM ) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            CONVER = .FALSE.
         END IF
* I suppose this is only a check for the necessity of defining
* a new search direction. Using intermediate ITERIM instead ...
CTF      ITER = ITER + 1
         ITERIM = ITER + 1
         IF( ITERIM .GT. MAXIT) GOTO 100
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new directio
* =====================================================================
         IF( .NOT. RTCNV(IROOT) ) THEN
           IADD = IADD + 1
*. Inverted diagonal times (HX-EX) on LU3
*. Real part
           CALL DMTVCD(VEC1,VEC2,LUDIA,LU4,LU3,-EIGAPR,1,1,LBLK)
*. Imag part
           CALL REWINE(LUDIA,LBLK)
           CALL DMTVCD(VEC1,VEC2,LUDIA,LU4,LU3,-EIGAPR,0,1,LBLK)
           IF ( IPRT  .GE. 600) THEN
             WRITE(6,*) '(D-E)-1 *( HX - EX ) '
             CALL WRTVCDC(VEC1,LU3,1,LBLK)
           END IF
* =====================================================================
*. 1.3 Orthogonalize to all previous vectors
* =====================================================================
           CALL REWINE(LU3,LBLK)
           XNORM =  INPRDD(VEC1,VEC1,LU3,LU3,0,LBLK)
     &           +  INPRDD(VEC1,VEC1,LU3,LU3,0,LBLK)
           XNORMO = XNORM
           DO JVEC = 1,NVEC+IADD-1
C    9            INPRDDC(VEC1,VEC2,LU1,LU2,IOFF1,IOFF2,IREW,LBLK,XR,XI)
             CALL INPRDDC(VEC1,VEC2,LU1,LU3,JVEC,1,1,LBLK,
     &                    XR,XI)
             WORK(2*(JVEC-1)+1) = XR
             WORK(2*(JVEC-1)+2) = XI
             XNORM = XNORM - XR**2 - XI**2
           END DO
           THRSLDP = 1.0D-10
           IF(XNORM/XNORMO .LE. THRSLDP) THEN
*. Linear dependency
             WRITE(6,*)
     &       ' CMICDV : Vector eliminated due to linear dependence '
             IADD = IADD -1
             GOTO 100
           END IF
           SCALE = 1.0D0/SQRT(XNORM)
C          write(6,*) ' Scale, scale**2*xnorm',scale,
C    &                scale**2*xnorm
           CALL SCALVE(WORK,SCALE,2*(NVEC+IADD-1) )
*. Real part of orthogonalized vector on LU5
           DO JVEC = 1, NVEC+IADD-1
             WORK(2*(JVEC-1)+2) =
     &      -WORK(2*(JVEC-1)+2)
           END DO
           CALL MVCSMD(LU1,WORK,LU4,LU5,VEC1,VEC2,2*(NVEC+IADD-1),
     &                 1,LBLK)
C               MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
           CALL REWINE(LU3,LBLK)
           CALL REWINE(LU4,LBLK)
           CALL REWINE(LU5,LBLK)
           ONEM = -1.0D0
           CALL VECSMD(VEC1,VEC2,SCALE,ONEM,LU3,LU4,LU5,0,LBLK)
*. Imag part of othogonal vector
           DO JVEC = 1, NVEC+IADD-1
             BUF = WORK(2*(JVEC-1)+1)
             WORK(2*(JVEC-1)+1) = - WORK(2*(JVEC-1)+2)
             WORK(2*(JVEC-1)+2) = BUF
           END DO
           CALL MVCSMD(LU1,WORK,LU4,LU6,VEC1,VEC2,2*(NVEC+IADD-1),
     &                 1,LBLK)
C               MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
           CALL SKPVCD(LU3,1,VEC1,1,LBLK)
           CALL REWINE(LU4,LBLK)
           ONEM = -1.0D0
           CALL VECSMD(VEC1,VEC2,SCALE,ONEM,LU3,LU4,LU5,0,LBLK)
           IF ( IPRT  .GE. 600 ) THEN
             WRITE(6,*) '   Orthonormalized (D-E)-1 *( HX - EX ) '
             CALL WRTVCDC(VEC1,LU5,1,LBLK)
           END IF
C               SKPVCD(LU,NVEC,SEGMNT,IREW,LBLK)
           CALL SKPVCD(LU1,2*(NVEC+IADD-1),VEC1,1,LBLK)
           CALL REWINE(LU5,LBLK)
           CALL COPVCDC(LU5,LU1,VEC1,0,LBLK)
         END IF
  100 CONTINUE
CTF
      ITER = ITER + 1
CTF
      IF( CONVER ) THEN
        ITER = ITER - 1
        GOTO  901
      END IF
*
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
* =======================================================
**  2 : Optimal combination of new and old directions
* =======================================================
*
*  2.1: Multiply new directions with matrix
      CALL SKPVCD(LU1,2*NVEC,VEC1,1,LBLK)
      CALL SKPVCD(LU2,2*NVEC,VEC1,1,LBLK)
      DO IVEC = 1, IADD
        CALL REWINE(LU3,LBLK)
        CALL COPVCDC(LU1,LU3,VEC1,0,LBLK)
        CALL REWINE(LU4,LBLK)
        IF(IPRT.GE.10000) THEN
          WRITE(6,*) ' Input vector to mV7 '
          CALL WRTVCDC(VEC1,LU3,1,LBLK)
        END IF
*
* Generalization to SIGDEN routine block implemented.
        if (SIGDEN_ROUTE.eq.'TRADIT') then
          call mv7(VEC1,VEC2,LU3,LU4,INONREL)
        else if (SIGDEN_ROUTE.eq.'SIGDEN') then
          call sigden_ctrl(VEC1,VEC2,LU3,LU4,T_CC,1)
        end if
*
        IF(IPRT.GE.10000) THEN
          WRITE(6,*) ' Output vector from MV7 '
          CALL WRTVCDC(VEC1,LU4,1,LBLK)
        END IF
        CALL REWINE(LU4,LBLK)
        CALL COPVCDC(LU4,LU2,VEC1,0,LBLK)
*. Augment projected matrix
        DO  JVEC = 1, NVEC+IVEC
          IJ = IVEC*(IVEC-1)/2 + JVEC
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          CALL INPRDDC(VEC1,VEC2,LU3,LU2,1,JVEC,1,LBLK,XR,XI)
C               INPRDDC(VEC1,VEC2,LU1,LU2,IOFF1,IOFF2,IREW,LBLK,XR,XI)
          APROJR(IJ) = XR
          APROJI(IJ) = XI
        END DO
      END DO
*
*. Diagonalize projected matrix
      NVEC = NVEC + IADD
      IF( IPRT  .GE. 10 ) THEN
        WRITE(6,*) ' Projected matrices before CEIGEN '
        CALL PRSYM(APROJR,NVEC)
        CALL PRSYM(APROJI,NVEC)
      END IF
      write(6,*) 'Dimension of projected matrix:',NVEC
      CALL CEIGEN(APROJR,APROJI,NVEC,EIG(1,ITER),AVEC,WORK)
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F22.13)')
     &   ( EIG(IROOT,ITER)+EIGSHF,IROOT=1,NROOT)
CJO      call towave(nroot,eig(1,iter))
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJR,NVEC)
        CALL PRSYM(APROJI,NVEC)
        CALL WRTMAT(AVEC,2*NVEC,NROOT,2*NVEC,NROOT)
      END IF
  901 CONTINUE
*
* =============================================
**  reset or assemble eigenvectors
* =============================================
*
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER.
     &   OR.IPREDZN.EQ.1)THEN
*. Reset C vectors
        CALL REWINE( LU5,LBLK)
        DO IROOT = 1, NROOT
*. Real part of vector IROOT on LU3
          DO JVEC = 1, NVEC
           WORK((JVEC-1)*2+1) = AVEC((IROOT-1)*2*NVEC+JVEC)
           WORK((JVEC-1)*2+2) =-AVEC((IROOT-1)*2*NVEC+NVEC+JVEC)
          END DO
          CALL MVCSMD(LU1,WORK,
     &    LU3,LU4,VEC1,VEC2,2*NVEC,1,LBLK)
*. Imaginary part of vector on LU6
          DO JVEC = 1, NVEC
           WORK((JVEC-1)*2+1) = AVEC((IROOT-1)*2*NVEC+NVEC+JVEC)
           WORK((JVEC-1)*2+2) = AVEC((IROOT-1)*2*NVEC+JVEC)
          END DO
          CALL MVCSMD(LU1,WORK,
     &    LU6,LU4,VEC1,VEC2,2*NVEC,1,LBLK)
          XNORM = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
     &          + INPRDD(VEC1,VEC1,LU6,LU6,1,LBLK)
          SCALE  = 1.0D0/SQRT(XNORM)
C          write(6,*) ' Scale, scale**2*xnorm',scale,
C    &                scale**2*xnorm
          WORK(2*NROOT*NVEC+IROOT) = SCALE
*. Copy current vector to LU5
C?        WRITE(6,*) ' Current vector on LU3, IROOT = ', IROOT
C?        CALL WRTVCD(VEC1,LU3,1,LBLK)
*
          CALL REWINE(LU3,LBLK)
          CALL SCLVCD(LU3,LU5,SCALE,VEC1,0,LBLK)
*
C?        WRITE(6,*) ' Current vector on LU6, IROOT = ', IROOT
C?        CALL WRTVCD(VEC1,LU6,1,LBLK)
*
          CALL REWINE(LU6,LBLK)
          CALL SCLVCD(LU6,LU5,SCALE,VEC1,0,LBLK)
CE        IF(IPRT.GT.100) THEN
CE           WRITE(6,*) ' Normalized vector '
CE           call rewine(LU5,LBLK)
CE           CALL WRTVCDC(VEC1,LU5,1,LBLK)
CE        END IF
        END DO
*. Transfer C vectors to LU1
*
        if (IPRT.ge.1000) then
          WRITE(6,*) ' Initial vectors on LU1 '
          CALL REWINE(LU1,LBLK)
          DO IVEC = 1, NROOT
            WRITE(6,*) ' Vector ', IVEC
            CALL WRTVCDC(VEC1,LU1,0,LBLK)
            WRITE(6,*) ' Writing of vector completed '
          END DO
*
          WRITE(6,*) ' Vectors on LU5 '
          CALL REWINE(LU5,LBLK)
          DO IVEC = 1, NROOT
            WRITE(6,*) ' Vector ', IVEC
            CALL WRTVCDC(VEC1,LU5,0,LBLK)
            WRITE(6,*) ' Writing of vector completed '
          END DO
        end if
*
        CALL REWINE(LU1,LBLK)
        CALL REWINE(LU5,LBLK)
        DO IVEC = 1,NROOT
          CALL COPVCDC(LU5,LU1,VEC1,0,LBLK)
C?        WRITE(6,*) ' Vector copied : ', IVEC
        END DO
        IF(IPRT.GE.1000) THEN
          WRITE(6,*) ' New reset C-vectors '
          WRITE(6,*) ' =================== '
          CALL REWINO(LU1)
          DO IVEC = 1, NROOT
             CALL WRTVCDC(VEC1,LU1,0,LBLK)
          END DO
        END IF
*. corresponding sigma vectors
        CALL REWINE(LU5,LBLK)
        DO IROOT = 1, NROOT
*. Real part of vector IROOT on LU3
          DO JVEC = 1, NVEC
           WORK((JVEC-1)*2+1) = AVEC((IROOT-1)*2*NVEC+JVEC)
           WORK((JVEC-1)*2+2) =-AVEC((IROOT-1)*2*NVEC+NVEC+JVEC)
          END DO
          CALL MVCSMD(LU2,WORK,
     &    LU3,LU4,VEC1,VEC2,2*NVEC,1,LBLK)
*. Imaginary part of vector on LU6
          DO JVEC = 1, NVEC
           WORK((JVEC-1)*2+1) = AVEC((IROOT-1)*2*NVEC+NVEC+JVEC)
           WORK((JVEC-1)*2+2) = AVEC((IROOT-1)*2*NVEC+JVEC)
          END DO
          CALL MVCSMD(LU2,WORK,
     &    LU6,LU4,VEC1,VEC2,2*NVEC,1,LBLK)
          SCALE = WORK(2*NROOT*NVEC+IROOT)
*. Copy to LU5
          CALL REWINE(LU3,LBLK)
          CALL SCLVCD(LU3,LU5,SCALE,VEC1,0,LBLK)
          CALL REWINE(LU6,LBLK)
          CALL SCLVCD(LU6,LU5,SCALE,VEC1,0,LBLK)
        END DO
*. Transfer HC vectors to LU2
        CALL REWINE(LU2,LBLK)
        CALL REWINE(LU5,LBLK)
        DO IVEC = 1,NROOT
          CALL COPVCDC(LU5,LU2,VEC1,0,LBLK)
        END DO
        IF(IPRT.GE.1000) THEN
          WRITE(6,*) ' New reset sigma -vectors '
          WRITE(6,*) ' =================== '
          CALL REWINO(LU2)
          DO IVEC = 1, NROOT
             CALL WRTVCDC(VEC1,LU2,0,LBLK)
          END DO
        END IF
        NVEC = NROOT
*. Reset projected matrices
        CALL SETVEC(AVEC,0.0D0,2*NVEC**2)
        CALL SETVEC(APROJR,0.0D0,NVEC*(NVEC+1)/2)
        CALL SETVEC(APROJI,0.0D0,NVEC*(NVEC+1)/2)
        DO IROOT = 1,NROOT
          AVEC((IROOT-1)*2*NROOT+IROOT) = 1.0D0
          APROJR(IROOT*(IROOT+1)/2 ) = EIG(IROOT,ITER)
        END DO
*
      END IF
cbs
      EXSTOP =.FALSE. 
      INQUIRE(FILE='STOPPER',EXIST=EXSTOP)
      if (EXSTOP) then
      write(6,*) 'I will stop now, like requested '
      goto 1001
      endif
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER.AND.IPREDZN.EQ.0) THEN
        GOTO 1000
      ELSE IF(IPREDZN.EQ.1) THEN
        WRITE(6,*) ' PREDIA vectors written to file LUCVECT'
        RETURN
      END IF
 1001 CONTINUE

* ( End of loop over iterations )
*
      IF(IPREDZN.EQ.0) THEN
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 0 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
C!?      ITER = ITER - 1
         IF (IPRT .GE. 0 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
      END IF
      END IF
*
      IF ( IPRT .GE. 0 ) THEN
        CALL REWINE(LU1,LBLK)
        CALL REWINE(LU2,LBLK)
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          FINEIG(IROOT) = EIG(IROOT,ITER)
          WRITE(6,1190) FINEIG(IROOT)+EIGSHF
 1190     FORMAT(' The final approximation to eigenvalue ',F22.10)
          IF(IPRT.GE.10) THEN
            WRITE(6,1201)
 1201       FORMAT(1H0,'The final approximation to eigenvector')
            CALL WRTVCDC(VEC1,LU1,0,LBLK)
 1200       FORMAT(1H0,'The final sigma vector                ')
            CALL WRTVCDC(VEC1,LU2,0,LBLK)
          END IF
          WRITE(6,1300)
 1300     FORMAT(1H0,' Summary of iterations ',/,1H
     +          ,' ----------------------')
          WRITE(6,1310)
 1310     FORMAT
     &    (1H0,' Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(6,1340) I,EIG(IROOT,I)+EIGSHF,RNRM(IROOT,I)
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      ELSE
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(IROOT,ITER)+EIGSHF
 1601   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(IROOT,ITER),
     &                 RNRM(IROOT,1)/RNRM(IROOT,ITER)
 1607   CONTINUE
      END IF
C
      call qexit('CMICDV')
*
* Delete scratch units
      close(unit=LU3,status='DELETE')
      close(unit=LU4,status='DELETE')
      close(unit=LU5,status='DELETE')
      close(unit=LU6,status='DELETE')
*
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
