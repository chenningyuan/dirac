      SUBROUTINE GASANAS(C,LUC,NBL_C,NSSOA,NSSOB,
     &                 NOCTPA,NOCTPB,
     &                 MXPNGAS,IOCTPA,IOCTPB,NBLOCK,IBLOCK,
     &                 THRES,MAXTRM,NAEL,NBEL,
     &                 IASTR,IBSTR,IBLTP,NSMST,IUSLAB,
     &                 IOBLAB,NCPMT,WCPMT,NELFSPGP,NOCCLS,
     &                 NGAS,IOCCLS,
     &                 IASPGPTP,IBSPGPTP,
     &                 ICISTR,NORB,IPRNCIV,IPRNT)
*
* Analyze CI vector :
*
*      1) Print atmost MAXTRM  combinations with coefficients
*         larger than THRES
*         Currently the corresponding dets are not GIVEN !!
*
*      2) Number of coefficients in given range
*
*      3) Number of coefficients in given range for given
*         occupation class
*
* Jeppe Olsen , Jan. 1989 ,
*               Aug 1995 : GAS version
*               May 1997 : BLOCK driven
*
*´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
* Modified for relativistic CI calculations
*
*        Timo Fleig, Feb. 2000
*´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
*
*
*. If IUSLAB  differs from zero Character*6 array IOBLAB is used to identify
*  Orbitals
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(*)
*. General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION IBLTP(*)
      CHARACTER*6 IOBLAB(*)
      DIMENSION NELFSPGP(MXPNGAS,*)
      DIMENSION IOCCLS(NGAS,NOCCLS)
*. Specific input
      DIMENSION IBLOCK(8,NBL_C)
*. Output
      DIMENSION NCPMT(10,NOCCLS)
      DIMENSION WCPMT(10,NOCCLS)
*
      NTESTL = 0000
      NTEST = max(IPRNT,NTESTL)
*
      IF(IUSLAB.NE.0) THEN
         WRITE(6,*)
         WRITE(6,*)
     & ' Labels for orbitals are of the type n l ml starting with n = 1'
         WRITE(6,*)
     & ' so the user should not be  alarmed by labels like 1 f+3 '
         WRITE(6,*)
      END IF

      MINPRT = 0
      ITRM = 0
      IDET = 0
      IIDET = 0
      ILOOP = 0
      NCIVAR = 0
      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
      CNORM = 0.0D0
 2001 CONTINUE
      IF( ICISTR .GE. 2 ) CALL REWINO(LUC)
      IIDET = 0
      ILOOP = ILOOP + 1
      IF ( ILOOP  .EQ. 1 ) THEN
         XMAX = 1.0D0
         XMIN = 1.0D0/SQRT(10.0D0)
      ELSE
         XMAX = XMIN
         XMIN = XMIN/SQRT(10.0D0)
      END IF
      IF(XMIN .LT. THRES  ) XMIN =  THRES
      IF(IPRNCIV.EQ.1) THEN
*. Print in one shot
         XMAX = 3006.1956
         XMIN = 0.0D0
      END IF
      IDET = 0
*
      WRITE(6,*)
      WRITE(6,'(A,E10.4,A,E10.4)')
     & '  Printout of coefficients in interval  ',XMIN,' to ',XMAX
      WRITE(6,'(A)')
     & '  ========================================================='
      WRITE(6,*)
*
      if (NTEST.ge.10) then
         WRITE(6,'(A,I6)') ' GASANAS : NBL_C = ',NBL_C
         WRITE(6,'(A,I6)') ' GASANAS : IASPGPTP = ',IASPGPTP
         WRITE(6,'(A,I6)') ' GASANAS : IBSPGPTP = ',IBSPGPTP
      end if
*
      DO JBLOCK = 1, NBL_C
*
         IATP = IBLOCK(1,JBLOCK)
         IBTP = IBLOCK(2,JBLOCK)
         IASM = IBLOCK(3,JBLOCK)
         IBSM = IBLOCK(4,JBLOCK)
         if (NTEST.ge.5) write(6,'(A,X,I6,I6,I6,I6)')
     &        ' IATP,IBTP,IASM,IBSM = ',IATP,IBTP,IASM,IBSM
*. Obtain alpha strings of sym IASM and type IATP
         IDUM = 0
         CALL GETSTR_TOTSM_SPGP(IASPGPTP,IATP,IASM,NAEL,
     &        NASTR1,IASTR,
     &        NORB,0,IDUM,IDUM)
*
*. Obtain beta strings of sym IBSM and type IBTP
         IDUM = 0
         CALL GETSTR_TOTSM_SPGP(IBSPGPTP,IBTP,IBSM,NBEL,
     &        NBSTR1,IBSTR,
     &        NORB,0,IDUM,IDUM)
*
         IF(IBLTP(IASM).EQ.2) THEN
            IRESTR = 1
         ELSE
            IRESTR = 0
         END IF
*
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
         if (NTEST.ge.5) write(6,'(A,X,I9,I9)')
     &        ' NIA,NIB = ',NIA,NIB
*
         IMZERO = 0
         if (ICISTR.ge.2) then
*. Read in a Type-Type-symmetry block
*? Perhaps use GSTTBL here?
            CALL IFRMDS(IDET,1,-1,LUC)
            CALL FRMDSC(C,IDET,-1,LUC,IMZERO,IAMPACK)
            if (NTEST.ge.10) then
              write(6,'(A,I8)') ' Number of elements readin ',IDET
            end if
            if (NTEST.ge.1000) then
              write(6,*)
              write(6,*) 'Partition of C vector is: '
              call wrtmat(C,IDET,1,IDET,1)
            end if
            IDET = 0
         END IF
C?       write(6,*) ' IMZERO in gasanas:  ',IMZERO
         IF(IMZERO.NE.1) THEN
*
            IBBAS = 1
            IABAS = 1
*
            DO IB = IBBAS,IBBAS+NIB-1
               IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
                  MINIA = IB - IBBAS + IABAS
               ELSE
                  MINIA = IABAS
               END IF
               DO  IA = MINIA,IABAS+NIA-1
*
                  IF (ILOOP .EQ. 1 ) NCIVAR = NCIVAR + 1
                  IDET = IDET + 1
                  IF (XMAX .GE. ABS(C(IDET)) .AND.
     &                 ABS(C(IDET)).GT. XMIN ) THEN
                     ITRM = ITRM + 1
                     IIDET = IIDET + 1
                     IF( ITRM .LE. MAXTRM ) THEN
                        CNORM = CNORM + C(IDET) ** 2
                        WRITE(6,'(A)')
                        WRITE(6,'(A)')
     &                       '                 =================== '
                        WRITE(6,*)
*
                        WRITE(6,'(A,I8,A,E14.8)')
     &                   '  Coefficient of combination ',IDET,' is ',
     &                       C(IDET)
                        WRITE(6,'(A)')
     &                   '  Corresponding alpha - and beta string '
                        IF(IUSLAB.EQ.0) THEN
                           WRITE(6,'(4X,10I4)')
     &                          (IASTR(IEL,IA),IEL = 1, NAEL )
                           WRITE(6,'(4X,10I4)')
     &                          (IBSTR(IEL,IB),IEL = 1, NBEL )
                        ELSE
                           WRITE(6,'(4X,10(1X,A6))')
     &                      (IOBLAB(IASTR(IEL,IA)),IEL = 1, NAEL )
                           WRITE(6,'(4X,10(1X,A6))')
     &                      (IOBLAB(IBSTR(IEL,IB)),IEL = 1, NBEL )
                        END IF
                     END IF
                  END IF
               END DO
*     ^ End of loop over alpha strings
            END DO
*     ^ End of loop over beta strings
         END IF
*     ^ End of if statement for nonvanishing blocks
      END DO
*     ^ End of loop over blocks
      IF(IIDET .EQ. 0 ) WRITE(6,*) '   ( no coefficients )'
      IF( XMIN .GT. THRES .AND. ILOOP .LE. 30 ) GOTO 2001
*
      WRITE(6,'(A,E15.8)')
     &     '  Norm of printed CI vector .. ', CNORM
*
*.Size of CI coefficients
*
*
      IDET = 0
      IF(ICISTR .GE. 2 ) CALL REWINO(LUC)
      CALL ISETVC(NCPMT,0    ,10*NOCCLS)
      CALL SETVEC(WCPMT,0.0D0,10*NOCCLS)
      DO JBLOCK = 1, NBL_C
         IATP = IBLOCK(1,JBLOCK)
         IBTP = IBLOCK(2,JBLOCK)
         IASM = IBLOCK(3,JBLOCK)
         IBSM = IBLOCK(4,JBLOCK)
*
         IF(IBLTP(IASM).EQ.2) THEN
            IRESTR = 1
         ELSE
            IRESTR = 0
         END IF
*. Occupation class corresponding to given occupation
         JOCCLS = 0
         DO JJOCCLS = 1, NOCCLS
            IM_THE_ONE = 1
            DO IGAS = 1, NGAS
               IF(NELFSPGP(IGAS,IATP-1+IOCTPA)+
     &          NELFSPGP(IGAS,IBTP-1+IOCTPB).NE.IOCCLS(IGAS,JJOCCLS))
     &          IM_THE_ONE = 0
            END DO
            IF(IM_THE_ONE .EQ. 1 ) JOCCLS = JJOCCLS
         END DO
*
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
*
         IMZERO = 0
         IF( ICISTR.GE.2 ) THEN
*. Read in a Type-Type-symmetry block
            CALL IFRMDS(IDET,1,-1,LUC)
            CALL FRMDSC(C,IDET,-1,LUC,IMZERO,IAMPACK)
            IDET = 0
         END IF
         IF(IMZERO.EQ.0) THEN

            DO IB = IBBAS,IBBAS+NIB-1
               IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
                  MINIA = IB - IBBAS + IABAS
               ELSE
                  MINIA = IABAS
               END IF
               DO IA = MINIA,IABAS+NIA-1
                  IDET = IDET + 1
                  DO IPOT = 1, 10
                     IF(10.0D0 ** (-IPOT+1).GE.ABS(C(IDET)).AND.
     &                 ABS(C(IDET)).GT. 10.0D0 ** ( - IPOT )) THEN
                       NCPMT(IPOT,JOCCLS)=
     &                       NCPMT(IPOT,JOCCLS)+ 1
                       WCPMT(IPOT,JOCCLS)=
     &                       WCPMT(IPOT,JOCCLS)+ C(IDET) ** 2
                     END IF
                  END DO
*     ^ End of loop over powers of ten
               END DO
*     ^ End of loop over alpha strings
            END DO
*     ^ End of loop over beta strings
         END IF
*     ^ End of test for novanishing blocks
      END DO
*     ^ End of loop over blocks
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') '   Magnitude of CI coefficients '
      WRITE(6,'(A)') '  =============================='
      WRITE(6,'(A)')
      WACC = 0.0D0
      NACC = 0
      DO 300 IPOT = 1, 10
         W = 0.0D0
         N = 0
         DO 290 JOCCLS = 1, NOCCLS
            N = N + NCPMT(IPOT,JOCCLS)
            W = W + WCPMT(IPOT,JOCCLS)
 290     CONTINUE
         WACC = WACC + W
         NACC = NACC + N
         WRITE(6,'(A,I2,A,I2,3X,I9,X,E15.8,3X,E15.8)')
     &        '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC
 300  CONTINUE
*
      WRITE(6,*) ' Number of coefficients less than  10-11',
     &     ' IS  ',NCIVAR - NACC
*
      IF(NOCCLS.NE.1) THEN
         WRITE(6,'(A)')
         WRITE(6,'(A)')
     &    '   Magnitude of CI coefficients for each excitation level '
         WRITE(6,'(A)')
     &    '  ========================================================='
         WRITE(6,'(A)')
         DO 400 JOCCLS = 1, NOCCLS
            N = 0
            DO 380 IPOT = 1, 10
               N = N + NCPMT(IPOT,JOCCLS)
 380        CONTINUE
            IF(N .NE. 0 ) THEN
               WRITE(6,*)
               WRITE(6,'(A,15I3)')'       Occupation of active sets :',
     &              (IOCCLS(IGAS,JOCCLS),IGAS=1, NGAS)
               WRITE(6,'(A,I9)')
     &          '         Number of coefficients larger than 10-11 ', N
               WRITE(6,*)
               WACC = 0.0D0
               DO 370 IPOT = 1, 10
                  N =  NCPMT(IPOT,JOCCLS)
                  W =  WCPMT(IPOT,JOCCLS)
                  WACC = WACC + W
                  WRITE(6,'(A,I2,A,I2,3X,I9,1X,E15.8,3X,E15.8)')
     &                 '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC
 370           CONTINUE
            END IF
 400     CONTINUE
*
*. Total weight and number of dets per excitation level
*
         WRITE(6,'(A)')
         WRITE(6,'(A)')
     &        '   Total weight and number of SD''s (> 10 ** -11 )  : '
         WRITE(6,'(A)')
     &        '  ================================================='
         WRITE(6,'(A)')
       WRITE(6,*) '        N      Weight      Acc. Weight   Occupation '
       WRITE(6,*) ' ==================================================='
         WACC = 0.0D0
         DO 500 JOCCLS = 1, NOCCLS
            N = 0
            W = 0.0D0
            DO 480 IPOT = 1, 10
               N = N + NCPMT(IPOT,JOCCLS)
               W = W + WCPMT(IPOT,JOCCLS)
 480        CONTINUE
            WACC = WACC + W
            IF(N .NE. 0 ) THEN
               WRITE(6,'(1X,I9,3X,E9.4,7X,E9.4,2X,16(1X,I2))')
     &              N,W,WACC,(IOCCLS(IGAS,JOCCLS),IGAS=1,NGAS)
            END IF
 500     CONTINUE
      END IF
*
      RETURN
      END
