      subroutine reomoluc(SYMINT_M,ISM,JSM,KSM,LSM,
     &                    IOFF,JOFF,KOFF,LOFF,
     &                    NIORB,NJORB,NKORB,NLORB,INLEN,
     &                    SYMINT_L)
*
*
* Transfer integrals of given symmetry block from Molecule to
* internal LUCIA order
*
*
*. Input :
* SYMINT_M : Input block of MOLECULE_SWEDEN INTEGRAL
*. Output
* SYMINT_L : Complete list of LUCIA integrals
      implicit real*8 (A-H,O-Z)
*
*. Input
      dimension SYMINT_M(*)
*. Output
      DIMENSION SYMINT_L(*)
*
      NTEST = 0000
*
      if (NTEST.ge.1000) then
         write(6,*) 'This is reomoluc speaking'
         write(6,*) 'ISM = ',ISM
         write(6,*) '     NIORB = ',NIORB
         write(6,*) 'JSM = ',JSM
         write(6,*) '     NJORB = ',NJORB
         write(6,*) 'KSM = ',KSM
         write(6,*) '     NKORB = ',NKORB
         write(6,*) 'LSM = ',LSM
         write(6,*) '     NLORB = ',NLORB
         write(6,*) 'Input integrals:'
         write(6,*) 'Length of current block: ',INLEN
         do I=1,INLEN,1
            write(6,*) SYMINT_M(I)
         end do
      end if
*
*. Index restrictions
*
      IF(ISM.EQ.JSM) THEN
        IJ_RES = 1
      ELSE
        IJ_RES = 0
      END IF
*
      IF(KSM.EQ.LSM) THEN
       KL_RES = 1
      ELSE
       KL_RES = 0
      END IF
*
      IF(ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
       IJKL_RES = 1
      ELSE
       IJKL_RES = 0
      END IF
*
* K
*  L
*   I
*    J
*
      IJKL_M = 0
      DO K = KOFF, KOFF + NKORB - 1
        IF(KL_RES.EQ.1) THEN
          L_LAST = K
        ELSE
          L_LAST = LOFF + NLORB - 1
        END IF
        DO L = LOFF, L_LAST
          IF(IJKL_RES.EQ.1) THEN
            I_FIRST = K
          ELSE
            I_FIRST = IOFF
          END IF
          DO I = I_FIRST, IOFF + NIORB - 1
            IF(IJ_RES.EQ.1) THEN
              J_LAST = I
            ELSE
              J_LAST = JOFF + NJORB - 1
            END IF
            IF(IJKL_RES.EQ.1.AND.I.EQ.K) THEN
              J_FIRST = L
            ELSE
              J_FIRST = JOFF
            END IF
            DO J = J_FIRST, J_LAST
*. Well, go back to relative ordering
              IR = I - IOFF + 1
              JR = J - JOFF + 1
              KR = K - KOFF + 1
              LR = L - LOFF + 1
              if (NTEST.ge.1000) then
                write(6,'(A,4I4)') 'IR JR KR LR  ',IR,JR,KR,LR
                write(6,'(A,4I4)') 'I J K L  ',I,J,K,L
                write(6,'(A,4I4)') 'IOFF JOFF KOFF LOFF  ',
     &                              IOFF,JOFF,KOFF,LOFF
              end if
              IJKL_L = I2EAD2(ISM,IR,JSM,JR,KSM,KR,LSM,LR)
              IJKL_M = IJKL_M + 1
              SYMINT_L(IJKL_L) = SYMINT_M(IJKL_M)
              if (NTEST.ge.1000) then
                 WRITE(6,'(A,6I4,E22.12)')
     &           ' I J K L IJKL_M IJKL_L VAL', I,J,K,L,
     &           IJKL_M,IJKL_L, SYMINT_M(IJKL_M)
              end if
            end do
          end do
        end do
      end do
*
      if (NTEST.ge.1000) then
         write(6,*) ' First IJKL_M integrals of SYMINT '
         CALL WRTMAT(SYMINT_L,IJKL_M-1,1,IJKL_M-1,1)
      end if
      return
      end
