      SUBROUTINE GET_DIA_2E(UUUUC,UUUUE,UBBUE,IDENSI,ISPINFREE,
     &                      ISPSF_MO,
     &                      ENVIRO,SIGDEN_ROUTE,NIRR_PN,IPRNT)
*
* Obtain 2electron integrals needed for
* constructing the determinant diagonal
*
*===============================================
* Jeppe Olsen, Timo Fleig
*  SIGDEN route implemented, 
*     Timo Fleig, August 2001
*===============================================
*
* UUUUC(I,J) = ( I I  ! J  J)
* UUUUE(I,J) = ( I J  ! J  I)
* UBBUE(I,J) = ( I J_ ! J_ I)
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION UUUUC(NTOOB,NTOOB),UUUUE(NTOOB,NTOOB)
      DIMENSION UBBUE(NTOOB,NTOOB)
*
      character*6 ENVIRO,SIGDEN_ROUTE
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "ctcc.inc"
*
      dimension INDGATH(MXPORB),ISPSF_MO(MXNDGIRR,2),IKPERM(4)
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT)
*
      if (NTEST.ge.2) then
        write(6,*)
      ZERO = 0.D0
        write(6,*) ' //////////////////////////////////// '
        write(6,*) '  Entering get_dia_2e                 '
        write(6,*) ' //////////////////////////////////// '
      end if
*
* UUUUC and UUUUE obtained from integral list 3
*
* ( I I ! J J )
*
      if (NTEST.ge.2) then
        write(6,*)
        write(6,*) ' Getting diagonal Coulomb ints. (II|JJ):'
      end if
*
      if (SIGDEN_ROUTE.eq.'TRADIT') then
        DO ISM = NSMOB/2+1, NSMOB
          DO JSM = NSMOB/2+1, NSMOB
*. Pointer to first integral in list ( ISM ISM ! JSM JSM)
            IOFF_DBG = IOFF_DBG_LIST(1,ISM,ISM,JSM,JSM)
            NI = NTOOBS(ISM)
            NJ = NTOOBS(JSM)
            IOFF = IBSO(ISM)
            JOFF = IBSO(JSM)
            if (NTEST.ge.50) then
              write(6,*) 'I symmetry is ',ISM
              write(6,*) 'J symmetry is ',JSM
              write(6,*) 'LIST offset is ',IOFF_DBG
              write(6,*) 'I Orbitals in this symmetry (all GAS) ',NI
              write(6,*) 'J Orbitals in this symmetry (all GAS) ',NJ
              write(6,*) 'I offset for this symmetry ',IOFF
              write(6,*) 'J offset for this symmetry ',JOFF
            end if
            DO I = 1, NI
              DO J = 1, NJ
*. address of (II!JJ)
                IJ = IOFF_DBG-1 + ((J-1)*NJ+(J-1))*NI*NI
     &                          +  (I-1)*NI+I
                UUUUC(I-1+IOFF,J-1+JOFF) =
     &          WORK(KINT2D(1,1)-1+IJ)
                IADRIN = KINT2D(1,1)-1+IJ
                if (NTEST.ge.100) then
                   write(6,*) 'selected address:  ',IADRIN
                   write(6,*) '  selected integral: ',WORK(IADRIN)
                end if
              END DO
            END DO
*           ^ End of loop over orbitals
          END DO
        END DO
*       ^ End of loop over symmetries
*
      else if (SIGDEN_ROUTE.eq.'SIGDEN') then
        IKFAC = 1
        IKPERM(1) = 1
        IKPERM(2) = 2
        IKPERM(3) = 3
        IKPERM(4) = 4
        if (ENVIRO(1:6).eq.'DIRAC '.and.ISPINFREE.eq.0) then
* Storage mode for integrals contains "-" sign for beta/barred indices.
          IKFAC = -1
* Use Kramers symmetry on second index pair.
          IKPERM(1) = 1
          IKPERM(2) = 2
          IKPERM(3) = 4
          IKPERM(4) = 3
        end if  
        N12INTS = N1ELINT + N2ELINT
        do ISM = NSMOB/2+1, NSMOB
          do JSM = NSMOB/2+1, NSMOB
            NI = NTOOBS(ISM)
            NJ = NTOOBS(JSM)
            IOFF = IBSO(ISM)
            JOFF = IBSO(JSM)
            ISM_B = ISPSF_MO(ISM,1)
            JSM_B = ISPSF_MO(JSM,1)
            call orbind_gssm(ISM_B,NIRR_PN,INDGATH,NI,IPRNT)
            call orbind_gssm(JSM_B,NIRR_PN,INDGATH(1+NI),NJ,IPRNT)
            do I = 1,NI,1
              do J = 1,NJ,1
                if (INDGATH(I).ne.0.and.INDGATH(J+NI).ne.0) then
                  call geth2dia_sigden(UUUUC(I-1+IOFF,J-1+JOFF),
     &                                 WORK(KT_CC),WORK(KINDEX),
     &                                 INDGATH(I),INDGATH(J+NI),
     &                                 'COUL',WORK(KLIBSOBEX),N12INTS,
     &                                 I2DIATP,IKFAC,IKPERM,NTEST)
                else
                  UUUUC(I-1+IOFF,J-1+JOFF) = ZERO
                end if
              end do
            end do
          end do
        end do
      end if
*
*
* ( I J ! J I )
      if (NTEST.ge.2) then
        write(6,*)
        write(6,*) ' Getting diagonal exchange ints. (IJ|JI):'
      end if
*
      if (SIGDEN_ROUTE.eq.'TRADIT') then
        DO ISM =  NSMOB/2+1, NSMOB
          DO JSM  = NSMOB/2+1, NSMOB
*. Pointer to first integral in list ( ISM JSM ! JSM ISM)
            IOFF_DBG = IOFF_DBG_LIST(1,ISM,JSM,JSM,ISM)
            NI = NTOOBS(ISM)
            NJ = NTOOBS(JSM)
            IOFF = IBSO(ISM)
            JOFF = IBSO(JSM)
            if (NTEST.ge.50) then
              write(6,*) 'I symmetry is ',ISM
              write(6,*) 'J symmetry is ',JSM
              write(6,*) 'LIST offset is ',IOFF_DBG
              write(6,*) 'I Orbitals in this symmetry (all GAS) ',NI
              write(6,*) 'J Orbitals in this symmetry (all GAS) ',NJ
              write(6,*) 'I offset for this symmetry ',IOFF
              write(6,*) 'J offset for this symmetry ',JOFF
            end if
            DO I = 1, NI
             DO J = 1, NJ
*. address of (IJ!JI)
               IJ = IOFF_DBG-1 + ((I-1)*NJ+(J-1))*NJ*NI
     &                         +  (J-1)*NI+I
               UUUUE(I-1+IOFF,J-1+JOFF) =
     &         WORK(KINT2D(1,1)-1+IJ)
               IADRIN = KINT2D(1,1)-1+IJ
               if (NTEST.ge.100) then
                  write(6,*) 'selected address:  ',IADRIN
                  write(6,*) '  selected integral: ',WORK(IADRIN)
               end if
             END DO
            END DO
*           ^ End of loop over orbitals
          END DO
        END DO
*
      else if (SIGDEN_ROUTE.eq.'SIGDEN') then
        IKFAC = 1
        IKPERM(1) = 1
        IKPERM(2) = 2
        IKPERM(3) = 3
        IKPERM(4) = 4
        if (ENVIRO(1:6).eq.'DIRAC '.and.ISPINFREE.eq.0) then
          IKFAC = -1
* Use Kramers symmetry on second index pair.
          IKPERM(1) = 1
          IKPERM(2) = 2
          IKPERM(3) = 4
          IKPERM(4) = 3
        end if  
        do ISM = NSMOB/2+1, NSMOB
          do JSM = NSMOB/2+1, NSMOB
            NI = NTOOBS(ISM)
            NJ = NTOOBS(JSM)
            IOFF = IBSO(ISM)
            JOFF = IBSO(JSM)
            ISM_B = ISPSF_MO(ISM,1)
            JSM_B = ISPSF_MO(JSM,1)
            call orbind_gssm(ISM_B,NIRR_PN,INDGATH,NI,IPRNT)
            call orbind_gssm(JSM_B,NIRR_PN,INDGATH(1+NI),NJ,IPRNT)
            do I = 1,NI,1
              do J = 1,NJ,1
                if (INDGATH(I).ne.0.and.INDGATH(J+NI).ne.0) then
                  call geth2dia_sigden(UUUUE(I-1+IOFF,J-1+JOFF),
     &                                 WORK(KT_CC),WORK(KINDEX),
     &                                 INDGATH(I),INDGATH(J+NI),
     &                                 'EXC1',WORK(KLIBSOBEX),N12INTS,
     &                                 I2DIATP,IKFAC,IKPERM,NTEST)
                else
                  UUUUE(I-1+IOFF,J-1+JOFF) = ZERO
                end if
              end do
            end do
          end do
        end do
      end if
*
* ( I J_ ! J_ I)
      if (NTEST.ge.2) then
        write(6,*)
        write(6,*) ' Getting diagonal exchange ints. (IJ_|J_I):'
      end if
*
*( obtained from list 3 )
      if ((ENVIRO(1:6).eq.'DIRAC '.and.ISPINFREE.eq.0).or.
     &    (ENVIRO(1:6).eq.'MOLCAS'.and.IDENSI.ge.1)) then
        if (SIGDEN_ROUTE.eq.'TRADIT') then
          DO ISM =  NSMOB/2+1, NSMOB
            DO JSM  = NSMOB/2+1, NSMOB
*. Pointer to first integral in list ( ISM JSM ! JSM ISM)
              IOFF_DBG = IOFF_DBG_LIST(3,ISM,JSM,JSM,ISM)
              NI = NTOOBS(ISM)
              NJ = NTOOBS(JSM)
              IOFF = IBSO(ISM)
              JOFF = IBSO(JSM)
              if (NTEST.ge.50) then
                write(6,*) 'I symmetry is ',ISM
                write(6,*) 'J symmetry is ',JSM
                write(6,*) 'LIST offset is ',IOFF_DBG
                write(6,*) 'I Orbitals in this symmetry (all GAS) ',NI
                write(6,*) 'J Orbitals in this symmetry (all GAS) ',NJ
                write(6,*) 'I offset for this symmetry ',IOFF
                write(6,*) 'J offset for this symmetry ',JOFF
              end if
              DO I = 1, NI
               DO J = 1, NJ
*. address of (IJ_!J_I)
                 IJ = IOFF_DBG-1 + ((I-1)*NJ+(J-1))*NJ*NI
     &                           +  (J-1)*NI+I
                 UBBUE(I-1+IOFF,J-1+JOFF) =
     &           WORK(KINT2D(1,3)-1+IJ)
                 IADRIN = KINT2D(1,3)-1+IJ
                 if (NTEST.ge.100) then
                   write(6,*) 'selected address:  ',IADRIN
                   write(6,*) '  selected integral: ',WORK(IADRIN)
                 end if
               END DO
              END DO
*             ^ End of loop over orbitals
            END DO
          END DO
*         ^ end loop over symmetries
*
        else if (SIGDEN_ROUTE.eq.'SIGDEN') then
          IKFAC = 1
          IKPERM(1) = 1
          IKPERM(2) = 2
          IKPERM(3) = 3
          IKPERM(4) = 4
          if (ENVIRO(1:6).eq.'DIRAC '.and.ISPINFREE.eq.0) then
            IKFAC = -1
            IKPERM(1) = 1
            IKPERM(2) = 2
            IKPERM(3) = 3
            IKPERM(4) = 4
          end if  
          do ISM = NSMOB/2+1, NSMOB
            do JSM = NSMOB/2+1, NSMOB
              NI = NTOOBS(ISM)
              NJ = NTOOBS2(JSM)
              IOFF = IBSO(ISM)
              JOFF = IBSO2(JSM)
              ISM_B = ISPSF_MO(ISM,1)
              JSM_B = ISPSF_MO(JSM,2)
              call orbind_gssm(ISM_B,NIRR_PN,INDGATH,NI,IPRNT)
              call orbind_gssm(JSM_B,NIRR_PN,INDGATH(1+NI),NJ,IPRNT)
              do I = 1,NI,1
                do J = 1,NJ,1
                  if (INDGATH(I).ne.0.and.INDGATH(J+NI).ne.0) then
                    call geth2dia_sigden(UBBUE(I-1+IOFF,J-1+JOFF),
     &                                   WORK(KT_CC),WORK(KINDEX),
     &                                   INDGATH(I),INDGATH(J+NI),
     &                                   'EXC2',WORK(KLIBSOBEX),N12INTS,
     &                                   I2DIATP,IKFAC,IKPERM,NTEST)
                  else
                    UBBUE(I-1+IOFF,J-1+JOFF) = ZERO
                  end if
                end do
              end do
            end do
          end do
        end if
      end if
*     ^ end if non-vanishing integral lists
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Two-electron integrals for CI diagonal'
        WRITE(6,*)
        WRITE(6,*) ' (I I ! J J) matrix '
        WRITE(6,*)
        CALL WRTMAT(UUUUC,NTOOB,NTOOB,NTOOB,NTOOB)
        WRITE(6,*)
        WRITE(6,*) ' (I J ! J I) matrix '
        WRITE(6,*)
        CALL WRTMAT(UUUUE,NTOOB,NTOOB,NTOOB,NTOOB)
        if (IDENSI.ge.1.or.ISPINFREE.eq.0) then
          WRITE(6,*)
          WRITE(6,*) ' (I J_! J_I) matrix '
          WRITE(6,*)
          CALL WRTMAT(UBBUE,NTOOB,NTOOB,NTOOB,NTOOB)
        end if
      END IF
*
      RETURN
      END
