      SUBROUTINE GET_INT(XOUT,XINT,TYPE,ISEA,JSEA,KSEA,LSEA,
     &                           N12,IKFAC,IKPERM,
     &                           INDEX,NSTART,NEND,IPRNT,NALPHA,INT)
*
#include "implicit.inc"
#include "ipoist8.inc"
*
* XOUT the integral being found
* XINT all integrals
* ISEA,JSEA,KSEA,LSEA indices ijkl
* TYPE of integral
* NSTART,NEND Start of loop until end (over which integrals) 
* INDEX contains indices of integrals
* N12 all integrals
* NALPHA needed for signs when NALPHA eq 0 or 4 for first COUL loop
      character*4 TYPE
*
      dimension IKPERM(4),
     &          INDEX(4,*)
      dimension XINT(N12)
*
      NTESTL = 100
      NTEST = max(NTESTL,IPRNT)
*
      IERROR = 1
      MONE = -1.0D0
C      N2B = IOFF2(I2DIATP)
C     do INT = N2B,N12,1
      do INT = NSTART+1,NEND,1
        if (TYPE.eq.'COUL') then
          I = INDEX(IKPERM(1),INT)
          J = INDEX(IKPERM(2),INT)
          K = INDEX(IKPERM(3),INT) * IKFAC
          L = INDEX(IKPERM(4),INT) * IKFAC
        if (I.eq.ISEA.and.
     &      J.eq.JSEA.and.
     &      K.eq.KSEA.and.
     &      L.eq.LSEA) then
          if (NTEST.ge.100) then
C           write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
C    &            'Diagonal int type ',TYPE,
C    &             INT,(INDEX(II,INT),II=1,4),XINT(INT)
            print*,'I,J,K,L',I,J,K,L
            write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
     &            'int type ',TYPE,
     &             INT,(INDEX(II,INT),II=1,4),
     &             XINT(INT)
            write(6,*)
          end if
C         XOUT = XINT(INT)
          IF(NALPHA.EQ.4.OR.NALPHA.EQ.0) THEN
            XOUT = MONE * XINT(INT)
          ELSE
            XOUT = XINT(INT)
          END IF
          IERROR = 0
* and exit.
          goto 10
        end if
        else if (TYPE.eq.'EXC1') then
          I = INDEX(IKPERM(1),INT)
          J = INDEX(IKPERM(4),INT) * IKFAC
          K = INDEX(IKPERM(2),INT)
          L = INDEX(IKPERM(3),INT) * IKFAC
        if (I.eq.ISEA.and.
     &      J.eq.JSEA.and.
     &      K.eq.KSEA.and.
     &      L.eq.LSEA) then
          if (NTEST.ge.100) then
C           write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
C    &            'Diagonal int type ',TYPE,
C    &             INT,(INDEX(II,INT),II=1,4),XINT(INT)
            print*,'I,J,K,L',I,J,K,L
            write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
     &            'int type ',TYPE,
     &             INT,(INDEX(II,INT),II=1,4),
     &             XINT(INT)
            write(6,*)
          end if
C         XOUT = XINT(INT)
          XOUT = XINT(INT)
          IERROR = 0
* and exit.
          goto 10
        end if
        else if (TYPE.eq.'EXC2') then
C          I = INDEX(IKPERM(1),INT)
C          J = INDEX(IKPERM(4),INT)
C          K = INDEX(IKPERM(2),INT) * IKFAC
C          L = INDEX(IKPERM(3),INT) * IKFAC
          I = INDEX(IKPERM(1),INT)
          J = INDEX(IKPERM(4),INT)
          K = INDEX(IKPERM(3),INT) * IKFAC
          L = INDEX(IKPERM(2),INT) * IKFAC
* Rename: External ST ordering, internal TS ordering:
        if (I.eq.ISEA.and.
     &      J.eq.JSEA.and.
     &      K.eq.KSEA.and.
     &      L.eq.LSEA) then
          if (NTEST.ge.100) then
C           write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
C    &            'Diagonal int type ',TYPE,
C    &             INT,(INDEX(II,INT),II=1,4),XINT(INT)
            print*,'I,J,K,L',I,J,K,L
            write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
     &            'int type ',TYPE,
     &             INT,(INDEX(II,INT),II=1,4),
     &             XINT(INT)
            write(6,*)
          end if
C         XOUT = XINT(INT)
          XOUT = XINT(INT)
          IERROR = 0
* and exit.
          goto 10
        end if
        else if (TYPE.eq.'SPIN') then
          I = INDEX(IKPERM(1),INT)
          J = INDEX(IKPERM(2),INT) * IKFAC
          K = INDEX(IKPERM(3),INT) 
          L = INDEX(IKPERM(4),INT) * IKFAC
        if (I.eq.ISEA.and.
     &      J.eq.JSEA.and.
     &      K.eq.KSEA.and.
     &      L.eq.LSEA) then
          if (NTEST.ge.100) then
C           write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
C    &            'Diagonal int type ',TYPE,
C    &             INT,(INDEX(II,INT),II=1,4),XINT(INT)
            print*,'I,J,K,L',I,J,K,L
            write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
     &            'int type ',TYPE,
     &             INT,(INDEX(II,INT),II=1,4),
     &             XINT(INT)
            write(6,*)
          end if
C         XOUT = XINT(INT)
C          stop 'spin-orbit'
C          IF(NALPHA.EQ.4.OR.NALPHA.EQ.0) THEN
C            XOUT = MONE * XINT(INT)
C          ELSE
          print*,'getting spin-orbit'
            XOUT = MONE * XINT(INT)
C          END IF
          IERROR = 0
* and exit.
          goto 10
        end if
        end if
C
C temp fix for problems with na or nb = 4 with sign and index
C
        IF(INT.EQ.NEND.AND.IERROR.EQ.1.AND.TYPE.eq.'COUL') THEN
C          INT = 0
          do INT2 = NSTART+1,NEND,1
          I = INDEX(IKPERM(1),INT2)
          J = INDEX(IKPERM(4),INT2) * IKFAC
          K = INDEX(IKPERM(3),INT2) 
          L = INDEX(IKPERM(2),INT2) * IKFAC
        if (I.eq.ISEA.and.
     &      J.eq.JSEA.and.
     &      K.eq.KSEA.and.
     &      L.eq.LSEA) then
          if (NTEST.ge.100) then
C           write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
C    &            'Diagonal int type ',TYPE,
C    &             INT,(INDEX(II,INT),II=1,4),XINT(INT)
            print*,'I,J,K,L',I,J,K,L
            write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
     &            'int type ',TYPE,
     &             INT2,(INDEX(II,INT2),II=1,4),
     &             XINT(INT2)
            write(6,*)
            print*,'IERROR was 1'
          end if
C         XOUT = XINT(INT)
          XOUT = XINT(INT2)
          IERROR = 0
* and exit.
          goto 10
          end if
          end do
        END IF
        IF(INT2.GE.NEND.AND.IERROR.EQ.1.AND.TYPE.eq.'COUL') THEN
C          INT = 0
          do INT3 = NSTART+1,NEND,1
          I = INDEX(IKPERM(3),INT3)
          J = INDEX(IKPERM(4),INT3) 
          K = INDEX(IKPERM(1),INT3) * IKFAC
          L = INDEX(IKPERM(2),INT3) * IKFAC
        if (I.eq.ISEA.and.
     &      J.eq.JSEA.and.
     &      K.eq.KSEA.and.
     &      L.eq.LSEA) then
          if (NTEST.ge.100) then
C           write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
C    &            'Diagonal int type ',TYPE,
C    &             INT,(INDEX(II,INT),II=1,4),XINT(INT)
            print*,'I,J,K,L',I,J,K,L
            write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
     &            'int type ',TYPE,
     &             INT3,(INDEX(II,INT3),II=1,4),
     &             XINT(INT3)
            write(6,*)
            print*,'IERROR was 1'
          end if
C         XOUT = XINT(INT)
          IF(NALPHA.EQ.4.OR.NALPHA.EQ.0) THEN
            XOUT = MONE * XINT(INT3)
          ELSE
            XOUT = XINT(INT3)
          END IF
          IERROR = 0
* and exit.
          goto 10
          end if
          end do
        END IF
        IF(INT3.GE.NEND.AND.IERROR.EQ.1.AND.TYPE.eq.'COUL') THEN
C          INT = 0
          do INT4 = NSTART+1,NEND,1
          I = INDEX(IKPERM(3),INT4)
          J = INDEX(IKPERM(2),INT4) * IKFAC
          K = INDEX(IKPERM(1),INT4)
          L = INDEX(IKPERM(4),INT4) * IKFAC
        if (I.eq.ISEA.and.
     &      J.eq.JSEA.and.
     &      K.eq.KSEA.and.
     &      L.eq.LSEA) then
          if (NTEST.ge.100) then
C           write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
C    &            'Diagonal int type ',TYPE,
C    &             INT,(INDEX(II,INT),II=1,4),XINT(INT)
            print*,'I,J,K,L',I,J,K,L
            write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
     &            'int type ',TYPE,
     &             INT4,(INDEX(II,INT4),II=1,4),
     &             XINT(INT4)
            write(6,*)
            print*,'IERROR was 1'
          end if
C         XOUT = XINT(INT)
          XOUT = XINT(INT4)
          IERROR = 0
* and exit.
          goto 10
          end if
          end do
        END IF
        IF(INT.EQ.NEND.AND.IERROR.EQ.1.AND.TYPE.eq.'SPIN') THEN
C          INT = 0
          do INT2 = NSTART+1,NEND,1
          I = INDEX(IKPERM(1),INT2)
          J = INDEX(IKPERM(4),INT2) * IKFAC
          K = INDEX(IKPERM(3),INT2)
          L = INDEX(IKPERM(2),INT2) * IKFAC
        if (I.eq.ISEA.and.
     &      J.eq.JSEA.and.
     &      K.eq.KSEA.and.
     &      L.eq.LSEA) then
          if (NTEST.ge.100) then
C           write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
C    &            'Diagonal int type ',TYPE,
C    &             INT,(INDEX(II,INT),II=1,4),XINT(INT)
            print*,'I,J,K,L',I,J,K,L
            write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
     &            'int type ',TYPE,
     &             INT2,(INDEX(II,INT2),II=1,4),
     &             XINT(INT2)
            write(6,*)
            print*,'IERROR was 1'
          end if
C         XOUT = XINT(INT)
C          stop 'spin-orbit 2'
          print*,'getting 2 spin-orbit'
          XOUT = XINT(INT2)
          IERROR = 0
* and exit.
          goto 10
          end if
          end do
        END IF
        IF(INT2.GE.NEND.AND.IERROR.EQ.1.AND.TYPE.eq.'SPIN') THEN
C          INT = 0
          do INT3 = NSTART+1,NEND,1
          I = INDEX(IKPERM(3),INT3)
          J = INDEX(IKPERM(4),INT3) * IKFAC
          K = INDEX(IKPERM(1),INT3)
          L = INDEX(IKPERM(2),INT3) * IKFAC
        if (I.eq.ISEA.and.
     &      J.eq.JSEA.and.
     &      K.eq.KSEA.and.
     &      L.eq.LSEA) then
          if (NTEST.ge.100) then
C           write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
C    &            'Diagonal int type ',TYPE,
C    &             INT,(INDEX(II,INT),II=1,4),XINT(INT)
            print*,'I,J,K,L',I,J,K,L
            write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
     &            'int type ',TYPE,
     &             INT3,(INDEX(II,INT3),II=1,4),
     &             XINT(INT3)
            write(6,*)
            print*,'IERROR was 1'
          end if
C         XOUT = XINT(INT)
C          stop 'spin-orbit 2'
          print*,'getting 3 spin-orbit'
CFIRST          XOUT = XINT(INT3)
          XOUT = MONE * XINT(INT3)
          IERROR = 0
* and exit.
          goto 10
          end if
          end do
        END IF
        IF(INT3.GE.NEND.AND.IERROR.EQ.1.AND.TYPE.eq.'SPIN') THEN
C          INT = 0
          do INT4 = NSTART+1,NEND,1
          I = INDEX(IKPERM(3),INT4)
          J = INDEX(IKPERM(2),INT4) * IKFAC
          K = INDEX(IKPERM(1),INT4)
          L = INDEX(IKPERM(4),INT4) * IKFAC
        if (I.eq.ISEA.and.
     &      J.eq.JSEA.and.
     &      K.eq.KSEA.and.
     &      L.eq.LSEA) then
          if (NTEST.ge.100) then
C           write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
C    &            'Diagonal int type ',TYPE,
C    &             INT,(INDEX(II,INT),II=1,4),XINT(INT)
            print*,'I,J,K,L',I,J,K,L
            write(6,'(A,X,A,I5,3X,4I3,1F18.10)')
     &            'int type ',TYPE,
     &             INT4,(INDEX(II,INT4),II=1,4),
     &             XINT(INT4)
            write(6,*)
            print*,'IERROR was 1'
          end if
C         XOUT = XINT(INT)
C          stop 'spin-orbit 2'
          print*,'getting 4 spin-orbit'
CFIRST          XOUT = MONE * XINT(INT4)
          XOUT = XINT(INT4)
          IERROR = 0
* and exit.
          goto 10
          end if
          end do
        END IF
      end do
      if (IERROR.eq.1) then
        write(6,*) 'No integral found in geth2dia_sigden.'
C       write(6,*) 'N2B,N12 = ',N2B,N12
        write(6,*) 'Type of integral = ',TYPE
        write(6,*) 'Indices ',I,J,K,L
        write(6,*) 'Indices input',ISEA,JSEA,KSEA,LSEA
        stop 'Quitting.'
      end if
*
  10  continue
C set counter to int for debug purposes
      print*,'INT,INT2,INT3,INT4',INT,INT2,INT3,INT4
      IF(INT.GE.NEND) THEN
        INT = INT2
        IF(INT2.GE.NEND) THEN
          INT = INT3
          IF(INT3.GE.NEND) THEN
            INT = INT4
          END IF
        END IF
      END IF
      return
      end
*
      SUBROUTINE GET_INT2(XOUT,XINT,ISEA,JSEA,KSEA,LSEA,
     &                                   NALPHA)
*
* Routine to fetch integral one by one from CC code.
* Done by calculating an offset on XINT, accieved by arranging indices
* Integral storage mode is assumed the default for luciarel hence integral indexing
* 
C  Plan is as follows
C 1. Identify which GAS space each index belong to
C 2. Make sure that this integral type will be there (only for NALPHA = 0 or 4)
C 3. Identify TYPE of integral
C 4. Calculate integral offset 
C 5. Fetch integral
C
C Could perhaps just pass down IHTYPE instead of finding it?
C ISUMORB and ISUBLENGHT could be replaced by Gauss's formula
C Parts can be reduced to what is found in FIND_TYPE_SPIN_ORBIT and thereby made faster
C
C ISEA = C1
C JSEA = A1
C KSEA = C2
C LSEA = A2
C
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*. Local scratch
      INTEGER I_INDX_TO_GAS(4), IMARKONEORBGAS(NGAS), ISUBLENGHT(NGAS+1)
      INTEGER IORB(NGAS+1), ITYPE2ALPHA(3), IPOS2ALPHA(3)
      DIMENSION XINT(*)
*
* Set some local variables
*
      NTEST = 00
      SIGN2 = 1.0D0
* 1. Identify the GAS and mark those with one orbital
      IORBCOUNT = 0
      ILENGHT = 0
      IEND = 0
      JEND = 0
      KEND = 0
      LEND = 0
      ISWAP = 0
      ISUBLENGHT(1) = 0
      IORB(1) = 0
      DO IGAS = 1,NGAS
        IORBCOUNT = IORBCOUNT + NGSOBT(IGAS)
        IF(IEND.EQ.0) THEN
          IF(IORBCOUNT.GE.ISEA) THEN
            I_INDX_TO_GAS(1) = IGAS
            IEND = 1
          END IF
        END IF
        IF(JEND.EQ.0) THEN
          IF(IORBCOUNT.GE.JSEA) THEN
            I_INDX_TO_GAS(2) = IGAS
            JEND = 1
          END IF
        END IF
        IF(KEND.EQ.0) THEN
          IF(IORBCOUNT.GE.KSEA) THEN
            I_INDX_TO_GAS(3) = IGAS
            KEND = 1
          END IF
        END IF
        IF(LEND.EQ.0) THEN
          IF(IORBCOUNT.GE.LSEA) THEN
            I_INDX_TO_GAS(4) = IGAS
            LEND = 1
          END IF
        END IF
        IF(NGSOBT(IGAS).EQ.1) IMARKONEORBGAS(IGAS) = 1
        IMARKONEORBGAS(IGAS) = 0
        ILENGHT = ILENGHT + IGAS
        ISUBLENGHT(IGAS+1) =  ILENGHT
        IORB(IGAS+1) = IORBCOUNT
      END DO
*
* 2. Finding the index of this type of integral in luciarel format
*    We here want the following condition on indices l>j and k>i for NALPHA 4 or 0
*    NALPHA 2 is not sorted like this
*
C      WRITE(6,*) 'ISEA,JSEA,KSEA,LSEA,start = ',ISEA,JSEA,KSEA,LSEA
      IF(NALPHA.EQ.4.OR.NALPHA.EQ.0.OR.ISPINORBIT.EQ.1) THEN
        IF(ABS(ISEA).GT.ABS(KSEA)) THEN
          IF(ABS(JSEA).GT.ABS(LSEA)) THEN
C swap all 4 indices and I_INDX_TO_GAS
C           print*,'enter 1'
            ISWAP = I_INDX_TO_GAS(1)
            I_INDX_TO_GAS(1) = I_INDX_TO_GAS(3)
            I_INDX_TO_GAS(3) = ISWAP
            ISWAP = I_INDX_TO_GAS(2)
            I_INDX_TO_GAS(2) = I_INDX_TO_GAS(4)
            I_INDX_TO_GAS(4) = ISWAP
            ISWAP = ISEA
            ISEA = KSEA
            KSEA = ISWAP
            ISWAP = JSEA
            JSEA = LSEA
            LSEA = ISWAP
            SIGN2 = -1.0D0 
          ELSE
C swap only 2 indices and I_INDX_TO_GAS
C           print*,'enter 2'
            ISWAP = I_INDX_TO_GAS(1)
            I_INDX_TO_GAS(1) = I_INDX_TO_GAS(3)
            I_INDX_TO_GAS(3) = ISWAP
            ISWAP = ISEA
            ISEA = KSEA
            KSEA = ISWAP
            SIGN2 = 1.0D0 !?
          END IF
        ELSE IF(ABS(JSEA).GT.ABS(LSEA)) THEN
C swap only 2 indices and I_INDX_TO_GAS
C           print*,'enter 3'
            ISWAP = I_INDX_TO_GAS(2)
            I_INDX_TO_GAS(2) = I_INDX_TO_GAS(4)
            I_INDX_TO_GAS(4) = ISWAP
            ISWAP = JSEA
            JSEA = LSEA
            LSEA = ISWAP
            SIGN2 = 1.0D0
        ELSE
C           print*,'enter 4'
            SIGN2 = -1.0D0
        END IF
      END IF
*
* 3. Identify relative type of integral and find offset of type and integral from lists
* Remember to check if there is only one kramers pair in gas shell. Done in step 1
* Should also add relative offset of type. Done
      IF(NALPHA.EQ.4.OR.NALPHA.EQ.0.OR.ISPINORBIT.EQ.1) THEN
        ITYPE = ISUBLENGHT(I_INDX_TO_GAS(4)) + I_INDX_TO_GAS(2)-1 +
     & ILENGHT*(ISUBLENGHT(I_INDX_TO_GAS(3)) + I_INDX_TO_GAS(1)-1)
*
        IF(NALPHA.EQ.4) THEN
          ITYPE = ITYPE + IRELTYPE(3)
C        print*,'ITYPE ABS',ITYPE
        ELSE IF(NALPHA.EQ.0) THEN
C          print*,'ITYPE,IRELTYPE(4)',ITYPE,IRELTYPE(4)
          ITYPE = ITYPE + IRELTYPE(4)
        ELSE
          ITYPE = ITYPE + IRELTYPE(8)
        END IF
      ELSE
C
C Here we need coulomb exchange2 and spin-orbit. Spin-orbit is like above
C Spin-orbit is here a special case since it is sorted like for NALPHA 4 or 0
C Will just make a small subroutine that goes like above to find it (if it is there!)
C Placed at the very end since no duplication of indices is then needed
        ITYPE =           (I_INDX_TO_GAS(2)-1) + 
     &          NGAS    * (I_INDX_TO_GAS(4)-1) + 
     &          NGAS**2 * (I_INDX_TO_GAS(3)-1) +
     &          NGAS**3 * (I_INDX_TO_GAS(1)-1)
        ITYPE2ALPHA(1) = ITYPE + IRELTYPE(5)
        ITYPE =           (I_INDX_TO_GAS(4)-1) +
     &          NGAS    * (I_INDX_TO_GAS(2)-1) +
     &          NGAS**2 * (I_INDX_TO_GAS(3)-1) +
     &          NGAS**3 * (I_INDX_TO_GAS(1)-1)
        ITYPE2ALPHA(2) = ITYPE + IRELTYPE(5) +
     &                  (IRELTYPE(6)-IRELTYPE(5))/2
C        print*,'ITYPE2ALPHA',ITYPE2ALPHA(1),ITYPE2ALPHA(2)
      END IF
* 4. Calculate position of integral
* First calculate position of integral in type then add offset for type
C Number of integrals can be found from array
      IF(NALPHA.EQ.4.OR.NALPHA.EQ.0.OR.ISPINORBIT.EQ.1) THEN
C Four possibilties
        IF(I_INDX_TO_GAS(1).EQ.I_INDX_TO_GAS(3).AND.
     &     I_INDX_TO_GAS(2).EQ.I_INDX_TO_GAS(4)) THEN
          IPOS = ISEA - IORB(I_INDX_TO_GAS(1)) - 1 +
     &           ISUMORB(KSEA - IORB(I_INDX_TO_GAS(3)) - 1) +
     &          (JSEA - IORB(I_INDX_TO_GAS(2)) - 1) * 
     &           ISUMORB(NGSOBT(I_INDX_TO_GAS(3))) +
     &           ISUMORB(LSEA - IORB(I_INDX_TO_GAS(4)) - 1) *
     &           ISUMORB(NGSOBT(I_INDX_TO_GAS(3)))
C          print*,'JSEA,LSEA,ISEA,KSEA',JSEA,LSEA,ISEA,KSEA
C          print*,'IORB(I_INDX_TO_GAS(2)-1)',IORB(I_INDX_TO_GAS(2)-1)
C          print*,'IPOS for 1',IPOS
          IPOS = IPOS + IOFFINTTYPE(ITYPE)
C          print*,'IPOS,IOFFINTTYPE(ITYPE) for 1',IPOS,IOFFINTTYPE(ITYPE)
        ELSE IF (I_INDX_TO_GAS(1).NE.I_INDX_TO_GAS(3).AND.
     &     I_INDX_TO_GAS(2).NE.I_INDX_TO_GAS(4)) THEN
          IPOS = (JSEA - IORB(I_INDX_TO_GAS(2))-1) * 
     &                   NGSOBT(I_INDX_TO_GAS(1)) * 
     &                   NGSOBT(I_INDX_TO_GAS(3)) * 
     &                   NGSOBT(I_INDX_TO_GAS(4)) +
     &           (LSEA - IORB(I_INDX_TO_GAS(4))-1) *
     &                   NGSOBT(I_INDX_TO_GAS(1)) *
     &                   NGSOBT(I_INDX_TO_GAS(3)) +
     &           (ISEA - IORB(I_INDX_TO_GAS(1))-1) *
     &                   NGSOBT(I_INDX_TO_GAS(3)) +
     &            KSEA - IORB(I_INDX_TO_GAS(3))-1
C          print*,'JSEA,LSEA,ISEA,KSEA',JSEA,LSEA,ISEA,KSEA
C          print*,'IORB(I_INDX_TO_GAS(2)-1)',IORB(I_INDX_TO_GAS(2)-1)
C          print*,'IPOS',IPOS
          IPOS = IPOS + IOFFINTTYPE(ITYPE)
C          print*,'IPOS,IOFFINTTYPE(ITYPE) for 2',IPOS,IOFFINTTYPE(ITYPE)
C          print*,'
        ELSE IF (I_INDX_TO_GAS(1).EQ.I_INDX_TO_GAS(3)) THEN
          IPOS = ISEA - IORB(I_INDX_TO_GAS(1)) - 1 +
     &           ISUMORB(KSEA - IORB(I_INDX_TO_GAS(3)) - 1) +
     &           (JSEA - IORB(I_INDX_TO_GAS(2))-1) *
     &                   NGSOBT(I_INDX_TO_GAS(4))  *
     &           ISUMORB(NGSOBT(I_INDX_TO_GAS(3))) +
     &           (LSEA - IORB(I_INDX_TO_GAS(4))-1) *
     &           ISUMORB(NGSOBT(I_INDX_TO_GAS(3)))
C          print*,'IPOS for 3',IPOS
          IPOS = IPOS + IOFFINTTYPE(ITYPE)
C          print*,'IPOS,IOFFINTTYPE(ITYPE) for 3',IPOS,IOFFINTTYPE(ITYPE)
        ELSE IF (I_INDX_TO_GAS(2).EQ.I_INDX_TO_GAS(4)) THEN
          IPOS = (ISEA - IORB(I_INDX_TO_GAS(1))-1) *
     &                   NGSOBT(I_INDX_TO_GAS(3)) +
     &            KSEA - IORB(I_INDX_TO_GAS(3))-1 +
     &           (JSEA - IORB(I_INDX_TO_GAS(2))-1) *
     &                   NGSOBT(I_INDX_TO_GAS(1)) *
     &                   NGSOBT(I_INDX_TO_GAS(3)) +
     &           ISUMORB(LSEA - IORB(I_INDX_TO_GAS(4)) - 1) *
     &                   NGSOBT(I_INDX_TO_GAS(1)) *
     &                   NGSOBT(I_INDX_TO_GAS(3)) 
C          print*,'IPOS 4',IPOS
          IPOS = IPOS + IOFFINTTYPE(ITYPE)
C          print*,'IPOS,IOFFINTTYPE(ITYPE) 4',IPOS,IOFFINTTYPE(ITYPE)
        END IF
      ELSE
          IPOS = (LSEA - IORB(I_INDX_TO_GAS(4))-1) *
     &                        NGSOBT(I_INDX_TO_GAS(1)) *
     &                        NGSOBT(I_INDX_TO_GAS(3)) *
     &                        NGSOBT(I_INDX_TO_GAS(2)) +
     &           (JSEA - IORB(I_INDX_TO_GAS(2))-1) *
     &                        NGSOBT(I_INDX_TO_GAS(1)) *
     &                        NGSOBT(I_INDX_TO_GAS(3)) +
     &           (KSEA - IORB(I_INDX_TO_GAS(3))-1) *
     &                        NGSOBT(I_INDX_TO_GAS(1)) +
     &            ISEA - IORB(I_INDX_TO_GAS(1))-1
C        print*,'ISEA,JSEA,KSEA,LSEA',ISEA,JSEA,KSEA,LSEA
        IPOS2ALPHA(1) = IPOS + IOFFINTTYPE(ITYPE2ALPHA(1))
          IPOS = (LSEA - IORB(I_INDX_TO_GAS(4))-1) *
     &                        NGSOBT(I_INDX_TO_GAS(1)) *
     &                        NGSOBT(I_INDX_TO_GAS(3)) *
     &                        NGSOBT(I_INDX_TO_GAS(2)) +
     &           (JSEA - IORB(I_INDX_TO_GAS(2))-1) *
     &                        NGSOBT(I_INDX_TO_GAS(1)) *
     &                        NGSOBT(I_INDX_TO_GAS(3)) +
     &           (KSEA - IORB(I_INDX_TO_GAS(3))-1) *
     &                        NGSOBT(I_INDX_TO_GAS(1)) +
     &            ISEA - IORB(I_INDX_TO_GAS(1))-1
        IPOS2ALPHA(2) = IPOS + IOFFINTTYPE(ITYPE2ALPHA(2))
C        print*,'IPOS',IPOS
C see if there is a spin-orbit contribution and find it!
        IF(ISPINORBIT.EQ.2) THEN
        IF(ISEA.NE.KSEA.AND.JSEA.NE.LSEA) THEN
          CALL FIND_TYPE_SPIN_ORBIT(ISEA,JSEA,KSEA,LSEA,I_INDX_TO_GAS,
     &                              ILENGHT,ISUBLENGHT,IORB,
     &                              SIGNS,IPOS2ALPHA(3))
          IPOS2ALPHA(3) = 0
        END IF
        ELSE
          IPOS2ALPHA(3) = 0
        END IF
C        print*,'IPOS2ALPHA',IPOS2ALPHA(1),IPOS2ALPHA(2),IPOS2ALPHA(3)
      END IF
* 5. Fetch integral from list
C
C      print*,'ipos',ipos
C      print*,'i am here'
      IF(NALPHA.EQ.4.OR.NALPHA.EQ.0.OR.ISPINORBIT.EQ.1) THEN
        XOUT = SIGN2*XINT(IPOS)
      ELSE
        IF(IPOS2ALPHA(3).EQ.0) THEN
          XOUT = XINT(IPOS2ALPHA(1)) -
     &           XINT(IPOS2ALPHA(2))
        ELSE
          XOUT = XINT(IPOS2ALPHA(1)) - 
     &           XINT(IPOS2ALPHA(2)) +
     &    SIGNS* XINT(IPOS2ALPHA(3))
        END IF
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'IORBCOUNT = ', IORBCOUNT 
        WRITE(6,*) 'ILENGHT = ',ILENGHT
        DO IGAS = 1,NGAS
          WRITE(6,*) 'ISUBLENGHT(IGAS) = ', ISUBLENGHT(IGAS)
        END DO
        DO IGAS = 1,NGAS
          WRITE(6,*) 'IORB(IGAS) = ', IORB(IGAS)
        END DO
        DO I = 1,4
          WRITE(6,*) 'I_INDX_TO_GAS(I) = ',I_INDX_TO_GAS(I)
        END DO
        WRITE(6,*) 'ISEA,JSEA,KSEA,LSEA = ',ISEA,JSEA,KSEA,LSEA
        IF(NALPHA.EQ.4.OR.NALPHA.EQ.0) THEN
          WRITE(6,*) 'ITYPE = ',ITYPE
        ELSE
          DO I=1,3
            WRITE(6,*) 'ITYPE2ALPHA(I) = ',ITYPE2ALPHA(I)
          END DO
        END IF
        IF(NALPHA.EQ.4.OR.NALPHA.EQ.0) THEN
          WRITE(6,*) 'IPOS = ',IPOS
        ELSE
          DO I=1,3
            WRITE(6,*) 'IPOS2ALPHA(I) = ',IPOS2ALPHA(I)
          END DO
        END IF
        WRITE(6,*) 'XOUT = ',XOUT
      END IF
      RETURN
      END
***********************************************************************
*                                                                     *
*     Routine to find the type related to a block                     *
*     of spin-orbit integrals and the index of the integral           *
*                                           Lasse                     *
***********************************************************************
*
      SUBROUTINE FIND_TYPE_SPIN_ORBIT(ISEA,JSEA,KSEA,LSEA,I_INDX_TO_GAS,
     &                                ILENGHT,ISUBLENGHT,IORB,
     &                                SIGNS,IPOS)
*
* Does exactly what GET_INT2 does for finding IHTYPE for NALPHA 4 or 0
* but for NALPHA 2. Hence permutes indices and finds type and then index.
* Will only fetch from the first block of 2-electron spin-orbit integrals.
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*. Local scratch
      INTEGER I_INDX_TO_GAS(4), ISUBLENGHT(NGAS+1), IORB(NGAS+1)
*
* 2. Finding the index of this type of integral in luciarel format
*    We here want the following condition on indices l>j and k>i
*
      IF(ABS(ISEA).GT.ABS(KSEA)) THEN
        IF(ABS(JSEA).GT.ABS(LSEA)) THEN
C swap all 4 indices and I_INDX_TO_GAS
C          print*,'enter 1 spin'
          ISWAP = I_INDX_TO_GAS(1)
          I_INDX_TO_GAS(1) = I_INDX_TO_GAS(3)
          I_INDX_TO_GAS(3) = ISWAP
          ISWAP = I_INDX_TO_GAS(2)
          I_INDX_TO_GAS(2) = I_INDX_TO_GAS(4)
          I_INDX_TO_GAS(4) = ISWAP
          ISWAP = ISEA
          ISEA = KSEA
          KSEA = ISWAP
          ISWAP = JSEA
          JSEA = LSEA
          LSEA = ISWAP
          SIGNS = -1.0D0 !?
        ELSE
C swap only 2 indices and I_INDX_TO_GAS
C          print*,'enter 2 spin'
          ISWAP = I_INDX_TO_GAS(1)
          I_INDX_TO_GAS(1) = I_INDX_TO_GAS(3)
          I_INDX_TO_GAS(3) = ISWAP
          ISWAP = ISEA
          ISEA = KSEA
          KSEA = ISWAP
          SIGNS = 1.0D0 !?
        END IF
      ELSE IF(ABS(JSEA).GT.ABS(LSEA)) THEN
C swap only 2 indices and I_INDX_TO_GAS
C          print*,'enter 3 spin'
          ISWAP = I_INDX_TO_GAS(2)
          I_INDX_TO_GAS(2) = I_INDX_TO_GAS(4)
          I_INDX_TO_GAS(4) = ISWAP
          ISWAP = JSEA
          JSEA = LSEA
          LSEA = ISWAP
          SIGNS = 1.0D0
      ELSE
C        print*,'enter 4 spin'
        SIGNS = -1.0D0
      END IF
*
* 3. Identify relative type of integral and find offset of type and integral from lists
* Remember to check if there is only one kramers pair in gas shell. Done in step 1
* Should also add relative offset of type. Done
      ITYPE = ISUBLENGHT(I_INDX_TO_GAS(4)) + I_INDX_TO_GAS(2)-1 +
     &ILENGHT*(ISUBLENGHT(I_INDX_TO_GAS(3)) + I_INDX_TO_GAS(1)-1)
*
      ITYPE = ITYPE + IRELTYPE(8)
* 4. Calculate position of integral
* First calculate position of integral in type then add offset for type
C Four possibilties
      IF(I_INDX_TO_GAS(1).EQ.I_INDX_TO_GAS(3).AND.
     &   I_INDX_TO_GAS(2).EQ.I_INDX_TO_GAS(4)) THEN
         IPOS = ISEA - IORB(I_INDX_TO_GAS(1)) - 1 +
     &          ISUMORB(KSEA - IORB(I_INDX_TO_GAS(3)) - 1) +
     &         (JSEA - IORB(I_INDX_TO_GAS(2)) - 1) *
     &          ISUMORB(NGSOBT(I_INDX_TO_GAS(3))) +
     &          ISUMORB(LSEA - IORB(I_INDX_TO_GAS(4)) - 1) *
     &          ISUMORB(NGSOBT(I_INDX_TO_GAS(3)))
C          print*,'IPOS for 1',IPOS
         IPOS = IPOS + IOFFINTTYPE(ITYPE)
C          print*,'IPOS,IOFFINTTYPE(ITYPE) for 1',IPOS,IOFFINTTYPE(ITYPE)
      ELSE IF (I_INDX_TO_GAS(1).NE.I_INDX_TO_GAS(3).AND.
     &         I_INDX_TO_GAS(2).NE.I_INDX_TO_GAS(4)) THEN
         IPOS = (JSEA - IORB(I_INDX_TO_GAS(2))-1) *
     &                  NGSOBT(I_INDX_TO_GAS(1)) *
     &                  NGSOBT(I_INDX_TO_GAS(3)) *
     &                  NGSOBT(I_INDX_TO_GAS(4)) +
     &          (LSEA - IORB(I_INDX_TO_GAS(4))-1) *
     &                  NGSOBT(I_INDX_TO_GAS(1)) *
     &                  NGSOBT(I_INDX_TO_GAS(3)) +
     &          (ISEA - IORB(I_INDX_TO_GAS(1))-1) *
     &                   NGSOBT(I_INDX_TO_GAS(3)) +
     &           KSEA - IORB(I_INDX_TO_GAS(3))-1
C          print*,'IPOS',IPOS
          IPOS = IPOS + IOFFINTTYPE(ITYPE)
C          print*,'IPOS,IOFFINTTYPE(ITYPE)',IPOS,IOFFINTTYPE(ITYPE)
      ELSE IF (I_INDX_TO_GAS(1).EQ.I_INDX_TO_GAS(3)) THEN
         IPOS = ISEA - IORB(I_INDX_TO_GAS(1)) - 1 +
     &          ISUMORB(KSEA - IORB(I_INDX_TO_GAS(3)) - 1) +
     &          (JSEA - IORB(I_INDX_TO_GAS(2))-1) *
     &                  NGSOBT(I_INDX_TO_GAS(4))  *
     &          ISUMORB(NGSOBT(I_INDX_TO_GAS(3))) +
     &          (LSEA - IORB(I_INDX_TO_GAS(4))-1) *
     &          ISUMORB(NGSOBT(I_INDX_TO_GAS(3)))
C         print*,'IPOS for 3',IPOS
         IPOS = IPOS + IOFFINTTYPE(ITYPE)
C         print*,'IPOS,IOFFINTTYPE(ITYPE) for 3',IPOS,IOFFINTTYPE(ITYPE)
      ELSE IF (I_INDX_TO_GAS(2).EQ.I_INDX_TO_GAS(4)) THEN
         IPOS = (ISEA - IORB(I_INDX_TO_GAS(1))-1) *
     &                  NGSOBT(I_INDX_TO_GAS(3)) +
     &           KSEA - IORB(I_INDX_TO_GAS(3))-1 +
     &          (JSEA - IORB(I_INDX_TO_GAS(2))-1) *
     &                  NGSOBT(I_INDX_TO_GAS(1)) *
     &                  NGSOBT(I_INDX_TO_GAS(3)) +
     &          ISUMORB(LSEA - IORB(I_INDX_TO_GAS(4)) - 1) *
     &                  NGSOBT(I_INDX_TO_GAS(1)) *
     &                  NGSOBT(I_INDX_TO_GAS(3))
C         print*,'IPOS 4',IPOS
         IPOS = IPOS + IOFFINTTYPE(ITYPE)
C         print*,'IPOS,IOFFINTTYPE(ITYPE) 4',IPOS,IOFFINTTYPE(ITYPE)
      END IF
*
      RETURN
      END
