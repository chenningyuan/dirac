      SUBROUTINE GASDIAS_REL(IMK2NUM,NMS2VAL,MS2VAL,NAEL,
     &           IASTR,NBEL,IBSTR,
     &           NORB,DIAG,NSMST,NIRREP,
     &           H1A,H1B,XA,XB,SCR,RJ,RK,
     &           RK2,NSSOA,NSSOB,LUDIA,ECORE,
     &           PLSIGN,PSSIGN,IPRNT,NTOOB,ICISTR,RJKAA,IHAM12,
     &           IBLTP,NBLOCK,IBLKFO,IASPGPTP,IBSPGPTP,
     &           IREOTS,NGAS,NIRR_PN,MXTSOB,ISPINFREE)
*
* Calculate determinant diagonal
* Turbo-ras version
*
* Driven by IBLKFO, May 97
*
* ========================
* General symmetry version
* ========================
*
* Jeppe Olsen, July 1995, GAS version
*
* Fixed bug in first exchange term.    Timo Fleig, 07.Aug.98
* Fixed bug in third integral term.    Timo and Jeppe, 25.Aug.98
* Correct sym-type ordering for all integrals implemented.
*                                      Timo Fleig, fall 2000
*
* IHAM12 = 1 => only one-body part
*        = 2 =>      one+two-body part
*
* separate alpha and beta 1-electron integrals  May 1998
      IMPLICIT REAL*8           (A-H,O-Z)
*
*.General input
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*)
      DIMENSION H1A(NORB), H1B(NORB)
*. Specific input
      DIMENSION IBLTP(*),IBLKFO(8,NBLOCK)
      dimension MS2VAL(NMS2VAL)
      dimension IREOTS(NORB)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION RK2(NTOOB,NTOOB)
*
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION RJKAA(*)
*. Output
      DIMENSION DIAG(*)
*
      call qenter('GASDIAS')
*
      NTEST =  0000
      NTEST = MAX(NTEST,IPRNT)
*
      IF(PSSIGN.EQ.-1.0D0) THEN
         XADD = 1000000.0
      ELSE
         XADD = 0.0D0
      END IF
*
      IF (NTEST.GE.20) THEN
        if (IMK2NUM.eq.1) then
          WRITE(6,*)
     &    ' Diagonal one electron integrals, unbarred'
          CALL WRTMAT(H1A,1,NORB,1,NORB)
          WRITE(6,*)
     &    ' Diagonal one electron integrals,   barred'
          CALL WRTMAT(H1B,1,NORB,1,NORB)
          WRITE(6,*) ' Core energy ', ECORE
          IF(IHAM12.EQ.2) THEN
            WRITE(6,*) ' Coulomb and exchange integrals '
            CALL WRTMAT(RJ,NTOOB,NTOOB,NTOOB,NTOOB)
            WRITE(6,*)
            CALL WRTMAT(RK,NTOOB,NTOOB,NTOOB,NTOOB)
            if (IDENSI.ge.1.or.ISPINFREE.eq.0) then
              WRITE(6,*)
              CALL WRTMAT(RK2,NTOOB,NTOOB,NTOOB,NTOOB)
            end if
          END IF
          write(6,*)
          write(6,*) 'IREOTS (T-S --> S-T) array'
          call iwrtma(IREOTS,1,NORB,1,NORB)
        end if
*
        WRITE(6,*) ' TTSS for Blocks '
        DO IBLOCK = 1, NBLOCK
          WRITE(6,'(10X,4I3,2I8)') (IBLKFO(II,IBLOCK),II=1,4)
        END DO

        WRITE(6,*) ' Number of alpha strings , type 1 '
        CALL IWRTMA(NSSOA,1,NSMST,1,NSMST)
        WRITE(6,*) ' Number of beta strings , type 1 '
        CALL IWRTMA(NSSOB,1,NSMST,1,NSMST)
*
      END IF
*
*  Diagonal elements according to Timo's formula
*
*   DIAG(IDET)
*          = sum_I (NIU*HII + NIB*HI_I_)
*            + sum_KM 1/2 * NKU*NMU [(KK|MM)     - (KM|MK)    ]
*                   + 1/2 * NKB*NMB [(K_K_|M_M_) - (K_M_|M_K_)]
*                   +       NKU*NMB [(KK|M_M_)   - (KM_|M_K)  ]
*
      IDET = 0
      ITDET = 0
*
      DO IBLK = 1, NBLOCK
        if (NTEST.ge.50) then
          write(6,*) '-----------------------------------'
          write(6,'(A14,I8,A13)') 
     &               ' Block number ',IBLK,' in diagonal '
          write(6,*) '-----------------------------------'
        end if
*
        IATP = IBLKFO(1,IBLK)
        IBTP = IBLKFO(2,IBLK)
        IASM = IBLKFO(3,IBLK)
        IBSM = IBLKFO(4,IBLK)
        if (NTEST.ge.10)
     &    WRITE(6,'(A21,4I4)') 
     &    ' IATP IBTP IASM IBSM ',IATP,IBTP,IASM,IBSM
*
        IF(IBLTP(IASM).EQ.2) THEN
          IREST1 = 1
        ELSE
          IREST1 = 0
        END IF
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP_REL(IASPGPTP,IATP,IASM,NAEL,NASTR1,IASTR,
     &                         NORB,0,IDUM,IDUM)
*
*. Obtain beta strings of sym IBSM and type IBTP
        CALL GETSTR_TOTSM_SPGP_REL(IBSPGPTP,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                         NORB,0,IDUM,IDUM)
*
        IBSTRT = 1
        IBSTOP =  NSSOB(IBSM,IBTP)
        DO IB = IBSTRT,IBSTOP
C?        WRITE(6,*)  ' IB = ', IB
          IBREL = IB - IBSTRT + 1
*
          IF(IREST1.EQ.1.AND.IATP.EQ.IBTP) THEN
            IASTRT =  IB
          ELSE
            IASTRT = 1
          END IF
          IASTOP = NSSOA(IASM,IATP)
C?        WRITE(6,*)' IASTRT,IASTOP', IASTRT,IASTOP
*
          DO IA = IASTRT,IASTOP
C?          WRITE(6,*)  ' IA = ', IA
            IDET = IDET + 1
            ITDET = ITDET + 1
*
            X = ECORE
*
            DO IEL = 1, NAEL
* Occupation of alpha spinors
              X = X + H1A(IREOTS(IASTR(IEL,IA)))
              if (NTEST.ge.1000) then
                 write(6,'(A,I3,A,2X,1F18.10)') 
     &               'Adding H1A(',IREOTS(IASTR(IEL,IA)),')',
     &                       H1A(IREOTS(IASTR(IEL,IA)))
              end if
            END DO
*
            DO IEL = 1, NBEL
* Occupation of beta spinors
              X = X + H1B(IREOTS(IBSTR(IEL,IB)))
              if (NTEST.ge.1000) then
                 write(6,'(A,I3,A,2X,1F18.10)') 
     &               'Adding H1B(',IREOTS(IBSTR(IEL,IB)),')',
     &                       H1B(IREOTS(IBSTR(IEL,IB)))
              end if
            END DO
*
            DO IEL = 1, NAEL
              DO JEL = 1, NAEL
                X = X + 0.5*
     &            RJ(IREOTS(IASTR(IEL,IA)),IREOTS(IASTR(JEL,IA)))
     &                - 0.5*
     &            RK(IREOTS(IASTR(IEL,IA)),IREOTS(IASTR(JEL,IA)))
                  if (NTEST.ge.1000) then
                     write(6,'(2X,A,2I3,A,1F18.10,3X,A,2I3,A,1F18.10)') 
     &                          'Adding 0.5 J(',
     &                          IREOTS(IASTR(IEL,IA)),
     &                          IREOTS(IASTR(JEL,IA)),
     &                          ') = ',
     &            0.5*RJ(IREOTS(IASTR(IEL,IA)),IREOTS(IASTR(JEL,IA))),
     &                          '- 0.5 K (',
     &                          IREOTS(IASTR(IEL,IA)),
     &                          IREOTS(IASTR(JEL,IA)),') =',
     &            0.5*RK(IREOTS(IASTR(IEL,IA)),IREOTS(IASTR(JEL,IA)))
                  end if
              END DO
            END DO
*
            DO IEL = 1, NBEL
              DO JEL = 1, NBEL
               X = X + 0.5*
     &            RJ(IREOTS(IBSTR(IEL,IB)),IREOTS(IBSTR(JEL,IB)))
     &               - 0.5*
     &            RK(IREOTS(IBSTR(IEL,IB)),IREOTS(IBSTR(JEL,IB)))
                  if (NTEST.ge.1000) then
                     write(6,'(2X,A,2I3,A,1F18.10,3X,A,2I3,A,1F18.10)') 
     &                          'Adding 0.5 J(',
     &                          IREOTS(IBSTR(IEL,IB)),
     &                          IREOTS(IBSTR(JEL,IB)),
     &                          ') = ',
     &            0.5*RJ(IREOTS(IBSTR(IEL,IB)),IREOTS(IBSTR(JEL,IB))),
     &                          '- 0.5 K (',
     &                          IREOTS(IBSTR(IEL,IB)),
     &                          IREOTS(IBSTR(JEL,IB)),') =',
     &            0.5*RK(IREOTS(IBSTR(IEL,IB)),IREOTS(IBSTR(JEL,IB)))
                  end if
              END DO
            END DO
*
            DO IEL = 1, NAEL
              DO JEL = 1, NBEL
                X = X + 
     &            RJ(IREOTS(IASTR(IEL,IA)),IREOTS(IBSTR(JEL,IB)))
                if (IDENSI.ge.1.or.ISPINFREE.eq.0) then
                  X = X - 
     &            RK2(IREOTS(IASTR(IEL,IA)),IREOTS(IBSTR(JEL,IB)))
                  if (NTEST.ge.1000) then
                     write(6,'(2X,A,2I3,A,1F18.10,3X,A,2I3,A,1F18.10)') 
     &                          'Adding J(',
     &                          IREOTS(IASTR(IEL,IA)),
     &                          -IREOTS(IBSTR(JEL,IB)),
     &                          ') = ',
     &            RJ(IREOTS(IASTR(IEL,IA)),IREOTS(IBSTR(JEL,IB))),
     &                          '- K2 (',
     &                          IREOTS(IASTR(IEL,IA)),
     &                          -IREOTS(IBSTR(JEL,IB)),') =',
     &            RK2(IREOTS(IASTR(IEL,IA)),IREOTS(IBSTR(JEL,IB)))
                  end if
                end if
              END DO
            END DO
*
            if (NTEST.ge.100) then
               write(6,*) 'Element ITDET and value X:'
               write(6,*) ITDET,X
               WRITE(6,*) ' alpha and beta string : '
               CALL IWRTMA(IASTR(1,IA),1,NAEL,1,NAEL)
               CALL IWRTMA(IBSTR(1,IB),1,NBEL,1,NBEL)
               WRITE(6,*)
            end if
            DIAG(IDET) = X
            IF(IB.EQ.IA) DIAG(IDET) = DIAG(IDET) + XADD
          END DO
*         ^ End of loop over alpha strings
        END DO
*       ^ End of loop over betastrings
*. Yet a RAS block of the diagonal has been constructed
2345    continue
        IF(ICISTR.GE.2) THEN
          IF(NTEST.GE.100) THEN
            write(6,*) ' number of diagonal elements to disc ',IDET
            CALL WRTMAT(DIAG,1,IDET,1,IDET)
          END IF
          CALL ITODS(IDET,1,-1,LUDIA)
          CALL TODSC(DIAG,IDET,-1,LUDIA)
          IDET = 0
        END IF
      END DO
*     ^ End of loop over blocks
      if (NTEST.ge.500) then
         WRITE(6,*)
     &   ' Number of diagonal elements generated ',ITDET
      end if
*
      IF (NTEST.GE.1000) THEN
        WRITE(6,*) ' CIDIAGONAL FIXME. (write statement)'
c       CALL WRTMAT(DIAG(1),ITDET,1,ITDET,1)
      END IF
*
      call qexit('GASDIAS')
      RETURN
      END
