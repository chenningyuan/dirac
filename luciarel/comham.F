*
      SUBROUTINE COMHAM(NVAR,T_CC,SIGDEN_ROUTE)
*
* Construct and diagonalize complete Hamiltonian
*
*
* LUDIA needed for giving form of vectors
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cprnt.inc"
#include "clunit.inc"
*
      character*6 SIGDEN_ROUTE
*
      common/CECORE/ECORE
*
      PARAMETER(MXDIM = 700)
      DIMENSION VEC1(MXDIM),VEC2(MXDIM),VEC3(MXDIM)
      DIMENSION SCR1((2*MXDIM)**2),SCR2((2*MXDIM)**2)
*
      NTEST = 050
*
      write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(6,*) 'MXDIM in COMHAM set to ',MXDIM
      write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
*
      LUIN = LUSC92
      LUOUT = LUSC93
      LBLK = -1
*
      LENH = NVAR**2
      call memman(KHR,LENH,'ADDS  ',2,'HREAL  ')
      call memman(KHI,LENH,'ADDS  ',2,'HIMAG  ')
*
*. Loop over elements in real part
      DO I = 1, NVAR
        WRITE(6,*) ' COMHAM : I = ', I
        CALL FLSHFO(6)
* Real part
        CALL REWINE(LUDIA,LBLK)
        CALL REWINE(LUIN,LBLK)
        ONE = 1.0D0
        CALL WRSVCD(LUIN,LBLK,VEC1,I,ONE,1,NVAR,LUDIA)
* Imaginary part
        CALL WRSVCD(LUIN,LBLK,VEC1,I,ONE,0,NVAR,LUDIA)
* Do the sigma
        CALL REWINE(LUOUT,LBLK)
        if (SIGDEN_ROUTE.eq.'TRADIT') then
          call mv7(VEC1,VEC2,LUIN,LUOUT,0)
        else if (SIGDEN_ROUTE.eq.'SIGDEN') then
          call sigden_ctrl(VEC1,VEC2,LUIN,LUOUT,T_CC,1)
        else
          write(6,*) 'Route for sigma and densities not allowed.'
          write(6,*) 'SIGDEN_ROUTE is ',SIGDEN_ROUTE
          call abend ('Quitting.')
        end if
*.
        CALL REWINE(LUOUT,LBLK)
        CALL MERGE_VEC(LUOUT,WORK(KHR+(I-1)*NVAR))
        CALL MERGE_VEC(LUOUT,WORK(KHI+(I-1)*NVAR))
C            MERGE_VEC(LU,VEC)
C       WRITE(6,*) I,' th Sigma vector read in in COMHAM '
C       CALL WRTMAT(WORK(KHR+(I-1)*NVAR),1,NVAR,1,NVAR)
      END DO
*
      if (NTEST.ge.50) then
         WRITE(6,*) ' Real and Imaginary parts of Hamiltonian '
         CALL WRTMAT(WORK(KHR),NVAR,NVAR,NVAR,NVAR)
         WRITE(6,*)
         WRITE(6,*)
         WRITE(6,*)
         CALL WRTMAT(WORK(KHI),NVAR,NVAR,NVAR,NVAR)
      end if
*
*. Check imag part > Antihermitian
      XMXDIFF = 1.0D-11
      DO I = 1, NVAR
        DO J = 1, I
          IF(ABS(WORK(KHI+(I-1)*NVAR+J-1)+WORK(KHI+(J-1)*NVAR+I-1)).
     &           GT.XMXDIFF) THEN
            Write(6,*)
            WRITE(6,*)
     &      ' Imag matrix not antihermitean, I,J, H(I,J),H(J,I), DIFF'
            WRITE(6,*) I,J,WORK(KHI+(I-1)*NVAR+J-1),
     &                     WORK(KHI+(J-1)*NVAR+I-1),
     &                     WORK(KHI+(I-1)*NVAR+J-1)
     &                     +WORK(KHI+(J-1)*NVAR+I-1)
          END IF
        END DO
      END DO
*. Check real part : Hermitian ?
      XMXDIFF = 1.0D-11
      DO I = 1, NVAR
        DO J = 1, I
          IF(ABS(WORK(KHR+(I-1)*NVAR+J-1)-WORK(KHR+(J-1)*NVAR+I-1))
     &                .GT.XMXDIFF) THEN
            Write(6,*)
            WRITE(6,*)
     &      ' Real matrix not hermitean, I,J, H(I,J),H(J,I), DIFF'
            WRITE(6,*) I,J,WORK(KHR+(I-1)*NVAR+J-1),
     &                     WORK(KHR+(J-1)*NVAR+I-1),
     &                     WORK(KHR+(I-1)*NVAR+J-1)
     &                     -WORK(KHR+(J-1)*NVAR+I-1)
          END IF
        END DO
      END DO
*
*. Diagonalize
*
*. But first reform to lower half form
      CALL TRIPAK(WORK(KHR),SCR1,1,NVAR,NVAR)
      CALL COPVEC(SCR1,WORK(KHR),NVAR*(NVAR+1)/2)
*. Read in two-electron terms from unit 30
*. Kept for historical (sentimental) reasons
C     IDIR2E = 1
C     IF(IDIR2E.EQ.1) THEN
C       LU30 =
C       REWIND (LU30)
C       DO IDIM = 1, NVAR
C         READ(LU30,*) (SCR1(IDIM*(IDIM-1)/2+JCOL),JCOL=1,IDIM)
C       END DO
C       WRITE(6,*) ' 2 el  terms read in '
C       CALL PRSYM(SCR1,NVAR)
*. And add
C       call VECSUM(HR,HR,SCR1,1.0D0,1.0D0,NVAR*(NVAR+1)/2)
C       write(6,*) ' Updated HR '
C       CALL PRSYM(HR,NVAR)
C     END IF
C          TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
      CALL TRIPAK(WORK(KHI),SCR1,1,NVAR,NVAR)
      CALL COPVEC(SCR1,WORK(KHI),NVAR*(NVAR+1)/2)
*
      write(6,*) 'NVAR : ',NVAR
      NVAR_SAVE = NVAR
      CALL CEIGEN(WORK(KHR),WORK(KHI),NVAR,VEC1,SCR1,SCR2)
C          CEIGEN(HR,HI,NDIM,EIGVAL,EIGVEC,SCR)
      write(6,*) 'NVAR_SAVE : ',NVAR_SAVE
*
*  Add nuclear repulsion energy
      DO I = 1, NVAR_SAVE
         VEC1(I) = VEC1(I) + ECORE
      END DO
*
      WRITE(6,*) ' Lowest Eigenvalue ', VEC1(1)
      WRITE(6,*) ' All eigenvalues '
      DO I = 1, NVAR_SAVE
         WRITE(6,'(X,I4,2X,1F18.10)') I,VEC1(I)
      END DO
*
*  Print MXANAVEC eigenvectors
      if (IPRWFC.ge.1) then
        do IANA=1,MXANAVEC,1
          write(6,*)
          write(6,'(A,X,I4)') ' Eigenvector number ... ',IANA
        end do
      end if
*
      close(unit=LUIN,status='DELETE')
      close(unit=LUOUT,status='DELETE')
*
      RETURN
      END
