C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rtractl */
      SUBROUTINE RTRACTL(ITRLVL,WORK,LWORK)
C*****************************************************************************
C
C     Input : ITRLVL - level of integrals transformation (see below)
C
C     Output: 4-index transformed integrals on disk.
C
C      
C     ITRLVL = -1: Transform all integrals (pq|rs)
C                     (which is more than we need, but for a quick start
C                      this is good)
C     ITRLVL =  0: Transform all integrals (uv|xy) for CI & MC energy
C     ITRLVL =  1: (pv|xy) for CI & MC gradients
C     ITRLVL =  2: (pq|xy) for Newton-Raphson and NEO CI & MC
C
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbopt.h>
#include <dgroup.h>
#include <dcborb.h>
C
      CHARACTER SZTRA(4,2)*80
      DIMENSION WORK(*)
C
#include <memint.h>
C
C
      CALL QENTER('RTRACTL')
C
      IF ( IPROPT .GT. 5 ) THEN
         CALL HEADER('Output from RTRACTL',-1)
         WRITE(LUPRI,'(A,I2)')
     $        '* Transformation level (ITRLVL) = ',ITRLVL
         WRITE(LUPRI,'(3A,I1,A)')
     $        '* Integrals in transformation = ',OPT_INTTYP,
     $        ' (INTFLG=',IOPT_INTFLG,')'
      END IF
C
      IF ( ITRLVL .LT. -1 .OR. ITRLVL .GT. 4 .OR. ITRLVL .EQ. 2 ) THEN
         WRITE(LUPRI,'(2(/,A),I2)')
     $        '*** ERROR in RTRACTL ***',
     $        'Illegal transformation level (ITRLVL) = ',ITRLVL
         CALL QUIT('*** ERROR in RTRACTL ***')
      END IF
C
C     Initialize all indeces with active orbitals.
C
      IF (ITRLVL .EQ. 1 .OR. ITRLVL .EQ. 0) THEN
         DO I = 1,4
            DO J = 1, NFSYM
               IF ( NASH (J) .GT. 0 ) THEN
                  WRITE(SZTRA(I,J),'(A,I4)')
     $                 '1..', NASH(J)
               ELSE
                  SZTRA(I,J) = ' '
               END IF
            END DO
         END DO
         IF (ITRLVL .EQ. 1 ) THEN
C           (uu|gu) integrals
            DO J = 1, NFSYM
               IF ( NESH (J) .GT. 0 ) THEN
                  IF ( NPSH (J) .GT. 0 ) THEN
                     WRITE(SZTRA(3,J),'(I4,A,I4)')
     $                 -NPSH(J), '..-1, 1..', NESH(J)
                  ELSE
                     WRITE(SZTRA(3,J),'(A,I4)')
     $                    '1..', NESH(J)
                  END IF
               ELSE
                  IF ( NPSH ( J ) .GT. 0 ) THEN
                     WRITE(SZTRA(3,J),'(I4,A)')
     $                    -NPSH(J), '..-1'
                  ELSE
                     SZTRA(3,J) = ' '
                  END IF
               END IF
            END DO
         END IF
      ELSE IF (ITRLVL .EQ. 2 ) THEN
C        cannot be implemented with the current structure...
      ELSE 
C        initialize all indeces with general orbitals
         DO I = 1, 4
            DO J = 1, NFSYM
               IF ( NESH (J) .GT. 0 ) THEN
                  IF ( NPSH (J) .GT. 0 ) THEN
                     WRITE(SZTRA(I,J),'(I4,A,I4)')
     $                    -NPSH(J), '..-1, 1..', NESH(J)
                  ELSE
                     WRITE(SZTRA(I,J),'(A,I4)')
     $                    '1..', NESH(J)
                  END IF
               ELSE
                  IF ( NPSH ( J ) .GT. 0 ) THEN
                     WRITE(SZTRA(I,J),'(I4,A)')
     $                    -NPSH(J), '..-1'
                  ELSE
                     SZTRA(I,J) = ' '
                  END IF
               END IF
            END DO
         END DO
         IF (ITRLVL .EQ. 3) THEN
C           (gu|gg) integrals
            DO J = 1, NFSYM
               IF ( NASH (J) .GT. 0 ) THEN
                  WRITE(SZTRA(2,J),'(A,I4)')
     $                 '1..', NASH(J)
               ELSE
                  SZTRA(2,J) = ' '
               END IF
            END DO
         END IF
      END IF
C
C     Print orbital strings
C
      IF ( IPROPT .GE. 5 ) THEN
         DO I = 1, 4
            DO J = 1,NFSYM
               WRITE(LUPRI,'(A,I1,A,I1,2A)')
     $              '* Orbital string ',I,
     $              ' for irrep ',J,' = ',SZTRA(I,J)
            END DO
         END DO
      END IF
      write(6,*) 'nopqfq ',opt_nopfq
C
C     Transfer information to moltra common blocks
C
      CALL SETUP_DCBTRA(SZTRA,IOPT_INTFLG,OPT_NOPFQ,0)
C
C     Call 4-index transformation driver
C
C
      CALL PAMTRA(WORK(KFREE),LFREE)
C
C
      CALL QEXIT('RTRACTL')
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rtractl */
      SUBROUTINE SETUP_DCBTRA(SZ,INTFLG,LNOPAIR,IPRINT)
C*****************************************************************************
C
C     Input : SZ     - orbitals strings
C             IPRINT - print level
C             INTFLG - integral classes to include in transformation
C
C     Output: none
C
C     Sets the common blocks used in the 4-index transformation.
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER ( DM1 = -1.0D00 )
C
      CHARACTER SZ(4,2)*(*)
      LOGICAL LNOPAIR
C
#include <dcbtra.h>
C     
C     /DCBTC/
C
      DO I = 1, 4
         DO J = 1, 2
            TRA4_INDSTR(I,J) = SZ(I,J)
         END DO
      END DO
      DO I = 1,2
         TRA_INDSTR(I) = ' '
         TRA_CORSTR(I) = ' '
         DO J = 1,2
            TRA2_INDSTR(I,J) = ' '
         END DO
      END DO
C
C     /DCBTRI/
C
      ITRA_INTFLG = INTFLG
      ITRA_INTFL2 = INTFLG
      ITRA_INTFL4 = INTFLG
      IPRTRA = IPRINT
      ISTRAT = 4
C
C     /DCBTRL/
C
      TRA_ANTIS = .FALSE.
      TRA_MSOUT  = .FALSE.
      NO2IND = .TRUE.
      NO4IND = .FALSE.
      NOPAIR = LNOPAIR
      PRPTRA = .FALSE.
      RCORBS = .FALSE.
      FKRMC = .TRUE.
      NOMDCINT = .TRUE.
C      
C     /DCBTRR/
C
      THROUT = 1.0D-16
      SCRTRA = DM1
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_rgeth2 */
      SUBROUTINE REAINTS(VUUUU,VGUUU,VGGUU,VGUGU,NDIMU,NDIMG,
     &     IFLAG,ITRLVL,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Read integrals and distribute in VUUUU, VGUUU, VGGUU and VGUGU.
C
C     IFLAG: bit 1: generate VUUUU
C            bit 2: generate VGUUU
C            bit 3: generate VGGUU
C            bit 4: generate VGUGU
C
C     Input:
C
C     Output:
C        VUUUU - integrals of type (uu|uu) 
C        VGUUU - integrals of type (gu|uu) 
C        VGGUU - integrals of type (gg|uu) 
C        VGUGU - integrals of type (gu|gu)
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbopt.h>
#include <dgroup.h>
C
      DIMENSION VUUUU(NDIMU,NDIMU,NDIMU,NDIMU,*)
      DIMENSION VGUUU(NDIMG,NDIMU,NDIMU,NDIMU,*)
      DIMENSION VGGUU(NDIMG,NDIMG,NDIMU,NDIMU,*)
      DIMENSION VGUGU(NDIMG,NDIMU,NDIMG,NDIMU,*)
      DIMENSION WORK(*)
C     index arrays
      DIMENSION NFPCK12(2), NFPCK34(2)
      DIMENSION NSTR(2,0:2,4)
      DIMENSION IDXG2U(1000)
C     logicals
      LOGICAL LBIT
      LOGICAL LUUUU,LGUUU,LGGUU,LGUGU
      LOGICAL LDZERO
C     strings
      CHARACTER CDATE*8,TIMEX*8,FNAM*5,FNODE*10,LFNAM*10
C
      DIMENSION ICLFAC(16), ICLFAC4(4)
C
      DATA ICLFAC4 / +1,-1,-1,+1 /
C
#include <sdpre.h>
#include <memint.h>
C
      CALL QENTER('REAINTS')
C
      write(6,*) 'iprnt',iprnt
      IPRINT = 4
C
C     Initialize ICLFAC
C
      CALL IZERO(ICLFAC,16)
      IF (NZ .EQ. 1) THEN
         DO I = 1, 4
            ICLFAC(I) = ICLFAC4(I)
         END DO
      ELSE IF (NZ .EQ. 2) THEN
         DO I = 0,3
            ICLFAC(2*I+1) = ICLFAC4(I+1)
            ICLFAC(2*I+2) = ICLFAC4(I+1)
         END DO
      ELSE
         DO I = 0,3
            ICLFAC(4*I+1) = ICLFAC4(I+1)
            ICLFAC(4*I+2) = ICLFAC4(I+1)
            ICLFAC(4*I+3) = ICLFAC4(I+1)
            ICLFAC(4*I+4) = ICLFAC4(I+1)
         END DO
      END IF
C
      IF (IFLAG .LT. 0) THEN
         LDZERO = .FALSE.
         IFLAG = -IFLAG
      ELSE
         LDZERO = .TRUE.
      END IF
C
      IF ( NDIMG .GE. 1000 ) 
     &     CALL QUIT('Internal error in REAINTS (NDIMG>1000)')
C
      LUUUU = LBIT(IFLAG,1)
      LGUUU = LBIT(IFLAG,2)
      LGGUU = LBIT(IFLAG,3)
      LGUGU = LBIT(IFLAG,4)
C
      IF (IPRINT .GE. 3) THEN
         CALL HEADER('Output from REAINTS',-1)
         WRITE(LUPRI,'(1X,A,I2)') '- integral flag (IFLAG): ',IFLAG
         WRITE(LUPRI,'(1X,A,I2)')
     &        '- integral transformation level (ITRLVL): ',ITRLVL
         IF (LUUUU)
     &        WRITE(LUPRI,'(1X,A)')
     &        '- distribute (uu|uu) integrals into VUUUU'
         IF (LGUUU)
     &        WRITE(LUPRI,'(1X,A)')
     &        '- distribute (gu|uu) integrals into VGUUU'
         IF (LGGUU)
     &        WRITE(LUPRI,'(1X,A)')
     &        '- distribute (gg|uu) integrals into VGGUU'
         IF (LGUGU)
     &        WRITE(LUPRI,'(1X,A)')
     &        '- distribute (gu|gu) integrals into VGUGU'
      END IF
C
C
C     Zero out integral arrays and check for consistency
C
      IF (LUUUU .AND. LDZERO) CALL DZERO(VUUUU, NDIMU**4 * NZ**2 )
      IF (LGUUU) THEN
         IF ( ITRLVL .GT. -1 .AND. ITRLVL .LT. 1 ) THEN
            WRITE(LUPRI,'(A,/,A,I3,A)')
     &           '*** ERROR in REAINTS ***',
     &           'Integral level = ',ITRLVL,
     &           ' but (gu|uu) integrals requested!'
            CALL QUIT('*** ERROR in REAINTS ***')
         END IF
         IF (LDZERO) CALL DZERO(VGUUU, NDIMG * NDIMU**3 * NZ**2 )
      END IF
      IF (LGGUU .OR. LGUGU) THEN
         IF ( ITRLVL .GT. -1 .AND. ITRLVL .LT. 2 ) THEN
            WRITE(LUPRI,'(A,/,A,I3,A)')
     &           '*** ERROR in REAINTS ***',
     &           'Integral level = ',ITRLVL,
     &           ' but (gg|uu) and/or (gu|gu) integrals requested!'
            CALL QUIT('*** ERROR in REAINTS ***')
         END IF
         IF (LDZERO) THEN
            CALL DZERO(VGGUU, NDIMG**2 * NDIMU**2 * NZ**2 )
            CALL DZERO(VGUGU, NDIMG**2 * NDIMU**2 * NZ**2 )
         END IF
      END IF
C
C     Make index vector general index --> active index
C
      IF ( ITRLVL .GT. 0 ) THEN
         NGENERAL = ITRLVL
      ELSE
         NGENERAL = 4
      END IF
      CALL MKIDXG2U(IDXG2U)
#ifdef UNDEF
      write(6,*) 'idxg2u = ',(IDXG2U(I),I=1,NDIMG)
#endif
C
C     Read integral info from MDINFO
C
      LUINT = 1 
      OPEN (LUINT, FILE = 'MDINFO', STATUS = 'OLD',
     &     FORM = 'UNFORMATTED' )
C
C     Read time & date for creation
C
      READ (LUINT) DATEX, TIMEX
C
      IF ( IPRINT .GE. 3 ) WRITE(LUPRI,'(4A)')
     &     '- time and date read from MDINFO: ',
     &     TIMEX, ' ',DATEX
C
      READ (LUINT) NFPCK12, NFPCK12T
      READ (LUINT) NFPCK34, NFPCK34T
      READ (LUINT) NSTR
C     
      CALL MEMGET('INTE',KDXKR12,2*NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB12,NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB34,NFPCK34T,WORK,KFREE,LFREE)
C
C     I and J indeces
      CALL READI(LUINT,NFPCK12T*2,WORK(KDXKR12))
#ifdef UNDEF      
      do i = 1, 2 * nfpck12t
         write(6,*) i,IGETINT(WORK(KDXKR12),I)
      end do
#endif
C     Boson irrep for IJ indeces
      CALL READI(LUINT,NFPCK12T,WORK(KDXB12))
C     Boson irrep for KL indeces
      CALL READI(LUINT,NFPCK34T,WORK(KDXB34))
C
      CLOSE (LUINT, STATUS = 'KEEP')
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Open integral files.
C
      LUINT = 0
      DO IREPIJ = 1, NFSYM
         IRECL = IRECLEN(NCLASS * NFPCK34 (IREPIJ), 0, 0)
         IRECL = MAX0(IRECL,IRECLEN(1,0,0))
         WRITE (FNAM,'(A4,I1)') '4IND',IREPIJ
         FNODE = LFNAM(FNAM)
         OPEN (LUINT+IREPIJ,FILE=FNODE,STATUS='OLD',
     &        ACCESS='DIRECT',RECL=IRECL)
      END DO
C
C     Read and distribute integrals.
C
C     - IJ is triangular.
C     - KL is not triangular.
C
      NBSYMP = NBSYM / NFSYM
C
      NIJOFF = 0
      NKLOFF = 0
      DO IREPIJ = 1, NFSYM
#ifdef UNDEF      
         write(6,*) 'REAINT IREPIJ = ',IREPIJ
#endif
         NIJ = NFPCK12(IREPIJ)
         NKL = NFPCK34(IREPIJ)
#ifdef UNDEF      
         write(6,*) 'REAINT NIJ = ',NIJ
         write(6,*) 'REAINT NKL = ',NKL
#endif
         CALL MEMGET('REAL',KHMAT,NKL*NZ*NZ*NBSYMP,WORK,KFREE,LFREE)
         DO IJ = 1, NIJ
#ifdef UNDEF      
            write(6,*) '  REAINT IJ = ', IJ
#endif
            IKR = IGETINT(WORK(KDXKR12),2*(IJ-1)+1 + NIJOFF)
            JKR = IGETINT(WORK(KDXKR12),2*IJ + NIJOFF)
            IREPIJB = IGETINT(WORK(KDXB12),IJ + NIJOFF)
#ifdef UNDEF      
            write(6,*) '  REAINT IKR = ', IKR
            write(6,*) '  REAINT JKR = ', JKR
            write(6,*) '  REAINT IREPIJB = ', IREPIJB
            call flshfo(6)
#endif
            READ (LUINT+IREPIJ,REC=IJ) (WORK(KHMAT+I-1),I=1,NKL*NCLASS)
#ifdef UNDEF      
            write(6,*) 'int ',(work(khmat+i-1),i=1,nkl*nclass)
#endif
C
            NZKL = 0
            IKR1 = IKR
            JKR1 = JKR
            KL = 0
            LL = 0
            DO IREPL = 1, NFSYM
               IREPK = MOD(IREPL + IREPIJ, 2 ) + 1
               DO L = 1, NSTR(IREPL,0,4)
                  LL = LL + 1
                  DO K = 1, NSTR(IREPK,0,3)
                     IF (IREPK .EQ. 1) THEN
                        KK = K
                     ELSE
                        KK = K + NSTR(1,0,3)
                     END IF
                     KL = KL + 1
C
c#ifdef UNDEF
                     WRITE(6,'(A,I2,A,I2,A,I2,A,I2,A,4(3X,D14.8))')
     &                    '    INT(',IKR,',',JKR,',',KK,',',LL,') = ',
     &                    (WORK(KHMAT + KL + (I-1)*NKL - 1),I=1,NCLASS)
c#endif
C
C                    Get active indeces from general/active indeces.
C
                     IF (ITRLVL .EQ. 0) THEN
                        IAC = IKR
                        JAC = JKR
                        KAC = KK
                        LAC = LL
                     ELSE IF (ITRLVL .EQ. 1) THEN
                        IAC = IKR
                        JAC = JKR
                        KAC = IDXG2U(KK)
                        LAC = LL
                     ELSE IF (ITRLVL .EQ. 2) THEN
C                       ...not implemented...
                        call quit('ralle ralle')
                        IAC = IDXG2U(IKR)
                        JAC = IDXG2U(JKR)
                        KAC = KK
                        LAC = LL
                     ELSE IF (ITRLVL .EQ. 3) THEN
                        IAC = IDXG2U(IKR)
                        JAC = JKR
                        KAC = IDXG2U(KK)
                        LAC = IDXG2U(LL)
                     ELSE 
                        IAC = IDXG2U(IKR)
                        JAC = IDXG2U(JKR)
                        KAC = IDXG2U(KK)
                        LAC = IDXG2U(LL)
                     END IF
C
C                    Distribute integrals
C                    
                     DO I = 1, NCLASS
                        IF (LUUUU .AND.
     &                       IAC .GT. 0 .AND. JAC .GT. 0 .AND.
     &                       KAC .GT. 0 .AND. LAC .GT. 0) THEN
                           VUUUU(KAC,LAC,IAC,JAC,I) =
     &                          WORK(KHMAT + KL + (I-1)*NKL - 1)
                           VUUUU(KAC,LAC,JAC,IAC,I) =
     &                          ICLFAC(I) *
     &                          WORK(KHMAT + KL + (I-1)*NKL - 1)
                        END IF
                        IF (LGUUU) THEN
                           IF (JAC .GT. 0 .AND. IAC .GT. 0
     &                          .AND. LAC .GT. 0) THEN
                              VGUUU(KK,LAC,IAC,JAC,I) = 
     &                             WORK(KHMAT + KL + (I-1)*NKL - 1)
                           END IF
                           IF (IAC .GT. 0 .AND. LAC .GT. 0
     &                          .AND. JAC .GT. 0) THEN
                              VGUUU(KK,LAC,JAC,IAC,I) =
     &                             ICLFAC(I) *
     &                             WORK(KHMAT + KL + (I-1)*NKL - 1)
                           END IF
                        END IF
                        IF (LGGUU) THEN
                           IF (IAC .GT. 0 .AND. IAC .GT. 0) THEN
                              VGGUU(KK,LL,IAC,JAC,I) = 
     &                             WORK(KHMAT + KL + (I-1)*NKL - 1)
                              VGGUU(KK,LL,JAC,IAC,I) = 
     &                             ICLFAC(I) *
     &                             WORK(KHMAT + KL + (I-1)*NKL - 1)
                           END IF
                        END IF
                        IF (LGUGU) THEN
                           IF (JAC .GT. 0 .AND. LAC .GT. 0) THEN
                              VGUGU(KK,LAC,IKR,JAC,I) = 
     &                             WORK(KHMAT + KL + (I-1)*NKL - 1)
                           END IF
                           IF (IAC .GT. 0 .AND. LAC .GT. 0) THEN
                              VGUGU(KK,LAC,JKR,IAC,I) =
     &                             ICLFAC(I) *
     &                             WORK(KHMAT + KL + (I-1)*NKL - 1)
                           END IF
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C     Close files
C
      DO IREPIJ = 1,NFSYM
         CLOSE (LUINT+IREPIJ, STATUS = 'KEEP' )
      END DO
C
 9000 FORMAT(1P,1X,4I3,4X,4(D14.7,1X),3(/,16X,4(D14.7,1X)))
C
      CALL MEMREL('REAINTS',WORK,KWORK,KWORK,KFREE,LFREE)
C
      IF (IPRINT .GE. 3) THEN
         IF (LUUUU) THEN
            CALL HEADER('VUUUU matrix',-1)
            DO I = 1,NDIMU
               DO J = 1,NDIMU
                  DO K = 1,NDIMU
                     DO L = 1,NDIMU
                        WRITE(LUPRI,9000)
     &                       I,J,K,L,(VUUUU(I,J,K,L,M),M=1,NCLASS)
                     END DO
                  END DO
               END DO
            END DO
         END IF
         IF (LGUUU) THEN
            CALL HEADER('VGUUUU matrix',-1)
            DO I = 1,NDIMG
               DO J = 1,NDIMU
                  DO K = 1,NDIMU
                     DO L = 1,NDIMU
                        WRITE(LUPRI,9000)
     &                       I,J,K,L,(VGUUU(I,J,K,L,M),M=1,NCLASS)
                     END DO
                  END DO
               END DO
            END DO
         END IF
         IF (LGGUU) THEN
            CALL HEADER('VGGUU matrix',-1)
            DO I = 1,NDIMG
               DO J = 1,NDIMG
                  DO K = 1,NDIMU
                     DO L = 1,NDIMU
                        WRITE(LUPRI,9000)
     &                       I,J,K,L,(VGGUU(I,J,K,L,M),M=1,NCLASS)
                     END DO
                  END DO
               END DO
            END DO
         END IF
         IF (LGUGU) THEN
            CALL HEADER('VGUGU matrix',-1)
            DO I = 1,NDIMG
               DO J = 1,NDIMU
                  DO K = 1,NDIMG
                     DO L = 1,NDIMU
                        WRITE(LUPRI,9000)
     &                       I,J,K,L,(VGUGU(I,J,K,L,M),M=1,NCLASS)
                     END DO
                  END DO
               END DO
            END DO
         END IF
      END IF
C
      CALL QEXIT('REAINTS')
C
      RETURN
10000 STOP 'END READING MDCINT'
10010 STOP 'ERROR READING MDCINT'
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck igetint */
      INTEGER FUNCTION IGETINT(IX,I)
C*****************************************************************************
C
C     Return element I of array IX
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
C
      DIMENSION IX(*)
C
      IGETINT = IX(I)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mkidxg2u */
      SUBROUTINE MKIDXG2U(IDXG2U)
C*****************************************************************************
C
C     Make index vector pointing from a general index to the 
C     corresponding active index.
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
C
      DIMENSION IDXG2U(*)
C
#include <dcborb.h>
#include <dgroup.h>
#include <dcbopt.h>
C
      IF ( OPT_NOPFQ ) THEN
C
C        No positronic indeces
C
         II = 0
         JJ = 0
         DO IFSYM = 1,NFSYM
            DO I = 1,NISH(IFSYM)
               II = II + 1
               IDXG2U(II) = -1
            END DO
            DO I = 1,NASH(IFSYM)
               II = II + 1
               JJ = JJ + 1
               IDXG2U(II) = JJ
            END DO
            DO I = 1,NSSH(IFSYM)
               II = II + 1
               IDXG2U(II) = -1
            END DO
         END DO
      ELSE
C
C        With positronic indeces
C
         II = 0
         JJ = 0
         DO IFSYM = 1,NFSYM
            DO I = 1,NISH(IFSYM)+NPSH(IFSYM)
               II = II + 1
               IDXG2U(II) = -1
            END DO
            DO I = 1,NASH(IFSYM)
               II = II + 1
               JJ = JJ + 1
               IDXG2U(II) = JJ
            END DO
            DO I = 1,NSSH(IFSYM)
               II = II + 1
               IDXG2U(II) = -1
            END DO
         END DO
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck D2MUUUU */
      SUBROUTINE D2MUUUU(VDUUUU,VMUUUU,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Transform a UUUU matrix in Dirac format to a UUUU matrix
C     in Molfdir format.
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
C
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION VDUUUU(NASHT,NASHT,NASHT,NASHT,NZ,NZ,*)
      DIMENSION VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
C
      DIMENSION NFPCK12(2), NFPCK34(2)
      DIMENSION NSTR(2,0:2,4)
C
#include <memint.h>
C
      CALL READ4IDX(NSTR,NFPCK12,NFPCK34,
     &     KDXKR12,KDXB12,KDXB34,WORK,KFREE,LFREE)
C
      CALL DZERO(VMUUUU,(2*NASHT)**2 * MIN0(2,NZ))
C
      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NSTR(IFSYM,0,1)
            II = II + 1
            JJ = JJ + 1
            IKRA(II) = JJ
            IKRB(II) = JJ + NSTR(IFSYM)
         ENDDO
         JJ = JJ + NSTR(IFSYM)
      ENDDO
C
      NIJOFF = 0
      NKLOFF = 0
      DO IREPIJ = 1, NFSYM
         NIJ = NFPCK12(IREPIJ)
         NKL = NFPCK34(IREPIJ)
            CALL D2MUUUU1(WORK(KFREE),LFREE,IPRINT,
     &        NSTR(1,0,3),NSTR(1,0,4),
     &        IREPIJ,NIJ,WORK(KDXKR12),NIJOFF,
     &        WORK(KDXB12),WORK(KDXB34),NKLOFF,
     &        VDUUUU,VMUUUU)
         NIJOFF = NIJOFF + NIJ
         NKLOFF = NKLOFF + NKL
      ENDDO
C
      CALL MEMREL('REAINTS',WORK,KWORK,KWORK,KFREE,LFREE)
C
 9000 FORMAT(1P,1X,4I3,4X,2(D14.7,1X))
C
      IF (IPRINT .GE. 0) THEN
         CALL HEADER('VUUUU matrix in molfdir format',-1)
         DO I = 1,2*NASHT
            DO J = 1,2*NASHT
               DO K = 1,2*NASHT
                  DO L = 1,2*NASHT
                     WRITE(LUPRI,9000)
     &                    I,J,K,L,(VMUUUU(I,J,K,L,IZ),IZ=1,MIN0(2,NZ))
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      RETURN
      END 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck d2muuuu1*/
      SUBROUTINE D2MUUUU1(WORK,LWORK,IPRINT,NSTR3,NSTR4,
     & IREPIJ,NIJ,INDXKR,NIJOFF,INDXB12,INDXB34,NKLOFF,
     &     VDUUUU,VMUUUU)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Driver for D2MUUUU2.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
C
      INTEGER   NSTR3(2),NSTR4(2)
      DIMENSION WORK(*)
      DIMENSION INDXKR(2,*),INDXB12(*),INDXB34(*)
      DIMENSION VDUUUU(*), VMUUUU(*)
C
#include <memint.h>
C
      DO IJ = 1, NIJ
         IKR = INDXKR(1,IJ+NIJOFF)
         JKR = INDXKR(2,IJ+NIJOFF)
         IREPIJB = INDXB12(IJ+NIJOFF)
         CALL D2MUUUU2(IPRINT,IREPIJ,IREPIJB,INDXB34(1+NKLOFF),IJ,
     &        IKR,JKR,NSTR3,NSTR4,NKL,VDUUUU,VMUUUU,
     &        WORK(KFREE),LFREE)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck D2MUUUU2 */
      SUBROUTINE D2MUUUU2(IPRINT,IREPIJ,IREPIJB,IREPKLB,IJ,IKR,JKR,
     &                   NSTR3,NSTR4,NKL,VDUUUU,VMUUUU,WORK,LWORK)
C
C     Symmetrize the final integrals.
C
C     Written by Joern Thyssen, Nov 3 2000
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcborb.h>
C
      DIMENSION WORK(*)
      INTEGER NSTR3(2),NSTR4(2)
      DIMENSION VDUUUU(NASHT,NASHT,NASHT,NASHT,NZ,NZ,*)
      DIMENSION VMUUUU(*)
      DIMENSION IREPKLB(NKL)
C
#include <memint.h>
C
      CALL MEMGET('REAL',KFTMAT,NKL*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDK,NKL,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDL,NKL,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRKL,2*NKL,WORK,KFREE,LFREE)
C
C     Real and complex classes
C     They are always written as complex for MOLFDIR
C
      DO ICLASS = 1, 4
         IREPABR = 0
         DO IREPAB = 0, NBSYM-1
            IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
               IREPABR = IREPABR + 1
               DO IZ1 = 1, NZ
                  IQ1 = IPQTOQ(IZ1,IREPAB)
                  DO IZ2 = 1, NZ
                     IQ2 = IPQTOQ(IZ2,IREPAB)
                     IFAC = INTCL(ICLASS,IQ1,IQ2)
                     IF (IFAC.NE.0) THEN
                        FACTR = IFAC
                        CALL DAXPY (NKL,FACTR,
     &                       VDUUUU(1,1,IKR,JKR,IZ2,IZ1,IREPABR),
     &                       1, WORK(KFTMAT),1)
                     ENDIF
                     IF (NZ.GT.1) THEN
                        IFAC=INTCL(ICLASS+4,IQ1,IQ2)
                        IF (IFAC.NE.0) THEN
                           FACTR = IFAC
                           CALL DAXPY (NKL,FACTR,
     &                          VDUUUU(1,1,IKR,JKR,IZ2,IZ1,IREPABR),1,
     &                          WORK(KFTMAT+NKL),1)
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
 100           CONTINUE
            ENDIF
         ENDDO
         CALL D2MUUUU3(IPRINT,ICLASS,IREPIJ,NSTR3,NSTR4,IKR,JKR,
     &        WORK(KFTMAT),NKL,
     &        WORK(KINDK),WORK(KINDL),
     &        VMUUUU)
      ENDDO
C
C     Quaternion classes
C
      IF (NZ.EQ.4) THEN
         DO ICLASS = 9, 12
            CALL DZERO(WORK(KFTMAT),2*NKL)
            DO IREPAB = 0, NBSYM-1
               IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
                  DO IZ1 = 1, NZ
                     IQ1 = IPQTOQ(IZ1,IREPAB)
                     DO IZ2 = 1, NZ
                        IQ2 = IPQTOQ(IZ2,IREPAB)
                        IFAC = INTCL(ICLASS,IQ1,IQ2)
                        IF (IFAC.NE.0) THEN
                            FACTR = IFAC
                            CALL DAXPY (NKL,FACTR,
     &                           VDUUUU(1,1,IKR,JKR,IZ2,IZ1,IREPABR),1,
     &                           WORK(KFTMAT),1)
                        ENDIF
                        IFAC=INTCL(ICLASS+4,IQ1,IQ2)
                        IF (IFAC.NE.0) THEN
                           FACTR = IFAC
                           CALL DAXPY (NKL,FACTR,
     &                          VDUUUU(1,1,IKR,JKR,IZ2,IZ1,IREPABR),1,
     &                          WORK(KFTMAT+NKL),1)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
            CALL D2MUUUU3(IPRINT,ICLASS,IREPIJ,NSTR3,NSTR4,IKR,JKR,
     &           WORK(KFTMAT),NKL,
     &           WORK(KINDK),WORK(KINDL),
     &           VMUUUU)
         ENDDO
      ENDIF
C
      CALL MEMREL('D2MUUUU2',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck d2muuuu3 */
      SUBROUTINE D2MUUUU3(IPRINT,ICLASS,IREPIJ,NSTR3,NSTR4,IKR,JKR,FMAT,
     &     NFM,INDK,INDL,VMUUUU)
C
C     blah blah
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
#include <dcborb.h>
      DIMENSION NSTR3(2),NSTR4(2)
      DIMENSION FMAT(NFM,2),INDK(*),INDL(*)
      DIMENSION VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,2)
C
      IKR1 = IKR
      JKR1 = JKR
      KL = 0
      LL = 0
      DO IREPL = 1, NFSYM
         IREPK = MOD(IREPL+IREPIJ,2) + 1
         DO L = 1, NSTR4(IREPL)
            LL = LL + 1
            DO K = 1, NSTR3(IREPK)
               IF (IREPK.EQ.1) THEN
                  KK = K
               ELSE
                  KK = K + NSTR3(1)
               ENDIF
               KL = KL + 1
               INDK(KL) = KK
               INDL(KL) = LL
            ENDDO
         ENDDO
      ENDDO
C
      IF (ICLASS.EQ.1) THEN
         CONTINUE
      ELSEIF (ICLASS.EQ.2) THEN
         DO INZ = 1, NZKL
            INDK(INZ) = - INDK(INZ)
            INDL(INZ) = - INDL(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.3) THEN
         JKR1 = - JKR1
         DO INZ = 1, NZKL
            INDL(INZ) = - INDL(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.4) THEN
         JKR1 = - JKR1
         DO INZ = 1, NZKL
            INDK(INZ) = - INDK(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.9) THEN
         DO INZ = 1, NZKL
            INDK(INZ) = - INDK(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.10) THEN
         DO INZ = 1, NZKL
            INDL(INZ) = - INDL(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.11) THEN
         IKR1 = - IKR1
      ELSEIF (ICLASS.EQ.12) THEN
         JKR1 = - JKR1
      ELSE
         CALL QUIT ('Class not defined in MDINTB')
      ENDIF
C
C
C     Distribute integrals in VMUUUU
C
      I   = KR(IKR)
      ITR = KR(-IKR)
      J   = KR(JKR)
      ITR = KR(-JKR)
      SIGNIJ = SIGN(1,IKR) * SIGN(1,JKR)
      DO INZ = 1, NZKL
         KKR = INDK(INZ)
         K = KR(KKR)
         KTR = KR(-KKR)
         LKR = INDL(INZ)
         L = KR(LKR)
         LTR = KR(-LKR)
         SignKL = SIGN(1,kkr) * SIGN(1,lkr)
C--> Original integral plus time-reversed partners
         VMUUUU(I,J,K,L,1) = FMAT(INZ,1)
         VMUUUU(I,J,K,L,2) = FMAT(INZ,2)
C
         VMUUUU(JTR,ITR,K,L,1) = SIGNIJ * FMAT(INZ,1)
         VMUUUU(JTR,ITR,K,L,2) = SIGNIJ * FMAT(INZ,2)
C
         VMUUUU(I,J,LTR,KTR,1) = SIGNKL * FMAT(INZ,1)
         VMUUUU(I,J,LTR,KTR,2) = SIGNKL * FMAT(INZ,2)
C
         VMUUUU(JTR,ITR,LTR,KTR,1) = SIGNIJ * SIGNKL * FMAT(INZ,1)
         VMUUUU(JTR,ITR,LTR,KTR,2) = SIGNIJ * SIGNKL * FMAT(INZ,2)
C--> Complex conjugate plus time-reversed partners
         VMUUUU(J,I,L,K,1) = FMAT(INZ,1)
         VMUUUU(J,I,L,K,2) = - FMAT(INZ,2)
C
         VMUUUU(ITR,JTR,L,K,1) = SIGNIJ * FMAT(INZ,1)
         VMUUUU(ITR,JTR,L,K,2) = - SIGNIJ * FMAT(INZ,2)
C
         VMUUUU(J,I,KTR,LTR,1) = SIGNKL * FMAT(INZ,1)
         VMUUUU(J,I,KTR,LTR,2) = - SIGNKL * FMAT(INZ,2)
C
         VMUUUU(ITR,JTR,KTR,LTR,1) = SIGNIJ * SIGNKL * FMAT(INZ,1)
         VMUUUU(ITR,JTR,KTR,LTR,2) = - SIGNIJ * SIGNKL * FMAT(INZ,2)
C
C--> Particle interchanged plus time-reversed partners
C
         VMUUUU(K,L,I,J,1) = FMAT(INZ,1)
         VMUUUU(K,L,I,J,2) = FMAT(INZ,2)
C
         VMUUUU(LTR,KTR,I,J,1) = SIGNKL * FMAT(INZ,1)
         VMUUUU(LTR,KTR,I,J,2) = SIGNKL * FMAT(INZ,2)
C
         VMUUUU(K,L,JTR,ITR,1) = SIGNIJ * FMAT(INZ,1)
         VMUUUU(K,L,JTR,ITR,2) = SIGNIJ * FMAT(INZ,2)
C
         VMUUUU(LTR,KTR,JTR,ITR,1) = SIGNKL * SIGNIJ * FMAT(INZ,1)
         VMUUUU(LTR,KTR,JTR,ITR,2) = SIGNKL * SIGNIJ * FMAT(INZ,2)
C
C--> Particle interchanged and complex conjugated plus time-reversed partners
C
         VMUUUU(L,K,J,I,1) = FMAT(INZ,1)
         VMUUUU(L,K,J,I,2) = - FMAT(INZ,2)
C
         VMUUUU(KTR,LTR,J,I,1) = SIGNKL * FMAT(INZ,1)
         VMUUUU(KTR,LTR,J,I,2) = - SIGNKL * FMAT(INZ,2)
C
         VMUUUU(L,K,ITR,JTR,1) = SIGNIJ * FMAT(INZ,1)
         VMUUUU(L,K,ITR,JTR,2) = - SIGNIJ * FMAT(INZ,2)
C
         VMUUUU(KTR,LTR,ITR,JTR,1) = SIGNKL * SIGNIJ * FMAT(INZ,1)
         VMUUUU(KTR,LTR,ITR,JTR,2) = - SIGNKL * SIGNIJ * FMAT(INZ,2)
C
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck read4idx */
      SUBROUTINE READ4IDX(NSTR,NFPCK12,NFPCK34,
     &     KDXKR12,KDXB12,KDXB34,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Read NSTR, IDXKR12, and other index vectors from MDINFO.
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
C
      DIMENSION NSTR(2,0:2,4)
      DIMENSION NFPCK12(2),NFPCK34(2)
      CHARACTER DATEX*8, TIMEX*8
C
C
C     Read integral info from MDINFO
C
      LUINT = 1 
      OPEN (LUINT, FILE = 'MDINFO', STATUS = 'OLD',
     &     FORM = 'UNFORMATTED' )
C
C     Read time & date for creation
C
      READ (LUINT) DATEX, TIMEX
C
      READ (LUINT) NFPCK12, NFPCK12T
      READ (LUINT) NFPCK34, NFPCK34T
      READ (LUINT) NSTR
C     
      CALL MEMGET('INTE',KDXKR12,2*NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB12,NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB34,NFPCK34T,WORK,KFREE,LFREE)
C
C     I and J indeces
      CALL READI(LUINT,NFPCK12T*2,WORK(KDXKR12))
C     Boson irrep for IJ indeces
      CALL READI(LUINT,NFPCK12T,WORK(KDXB12))
C     Boson irrep for KL indeces
      CALL READI(LUINT,NFPCK34T,WORK(KDXB34))
C
      CLOSE (LUINT, STATUS = 'KEEP')
C
      RETURN
      END




