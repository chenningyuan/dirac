      SUBROUTINE GET_HX_RELA(IFLAG,IHTYPE,NOBPT,NHX,
     &                       IFTONE,IFTTWO,ILTTWO,IHX_LOC,IHX,IHINDX,
     &                       SIGN_NHX,I2DIATP,LABEL,IOP_REO,SIGN_OPREO,
     &                       IPRNT)
*
* Obtain number (IFLAG = 1) 
* or number and types (IFLAG.NE.1) of excitations in 
* relativistic Hamiltonian 
*
*--------------------------------
* Jeppe Olsen and Timo Fleig
*   March/April 2001
* Orbital and electron number limitations
*   Timo Fleig, June 2001
* Updated to exclude unwanted components, Oct 2001
*   JO and TF
* Included LABEL specification of excitation type for DIRAC
*   (and other interfaces)
*--------------------------------
*
* IHTYPE = 1 : Include only Delta MK = 0 terms
* IHTYPE = 2 : Include Delta Mk = 0 for two-electron part 
*              and complete one-electron part
* IHTYPE = 3 : Include all terms for relativisitic Hamiltonian
*              in a basis of non-relativistic orbitals
* IHTYPE = 4 : Include all terms for relativisitic Hamiltonian
*              in a basis of relativistic spinors
*
* LABEL : 1  for unbarred (alpha)
*         -1 for barred   (beta)
*
* IOP_REO : Reordering array for output density elements.
*           Relate definition of operator in this routine with
*           general output CA AA CB AB (or corresponding scheme
*           in DMK .ne. 0 cases). Not used in 1-el. case, but
*           due to NHX type counting need to be defined.
*
* SIGN_OPREO : Sign for this operator reordering
*          1.0 or -1.0 : Non-redundant density and sign
*          0.0         : Redundant density type
*
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*. input
      dimension NOBPT(NGAS)
*. output 
      INTEGER IHX(NGAS,4,*),IHINDX(4,*),LABEL(4,*),IOP_REO(4,*)
      dimension SIGN_NHX(*),SIGN_OPREO(*)
*. Local scratch
      INTEGER IHX_LOC(4*NGAS)
*. An excitation is stored as : 
*  Alpha creation,beta creation,alpha annihilation,beta annihilation.
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT)
*
      IZERO = 0
      ZERO = 0.D0
      ONE = 1.0D0
*
      NHX = 0
*
* Determine maximum number of electrons per GAS first:
* FIXME  : Hardwired only one type of GASSPC calculation for the
*          moment. This would require other changes as well.
      MXNELGS(1) = IGSOCCX(1,2,1)
      do IGAS=2,NGAS,1
        MXNELGS(IGAS) = IGSOCCX(IGAS,2,1) - IGSOCCX(IGAS-1,1,1)
        MXNELGS(IGAS) = min(MXNELGS(IGAS),2*NOBPT(IGAS))
      end do
*
*. Delta MK = 0 terms
*
* a+i alpha a j alpha
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          if (NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1) then
              CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
              IHX_LOC(IGAS+(1-1)*NGAS) = 1
              IHX_LOC(JGAS+(3-1)*NGAS) = 1
*. Is this excitation type included ?
              CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN 
                WRITE(6,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  SIGN_NHX(NHX) = ONE
                  LABEL(IHINDX(1,NHX),NHX) = 1
                  LABEL(IHINDX(2,NHX),NHX) = 1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
                end if
              END IF
            end if
          end if
        END DO
      END DO
*
* a+i beta a j beta
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          if (NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1) then
              CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
              IHX_LOC(IGAS+(2-1)*NGAS) = 1
              IHX_LOC(JGAS+(4-1)*NGAS) = 1
              CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN 
                WRITE(6,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  SIGN_NHX(NHX) = ONE
                  LABEL(IHINDX(1,NHX),NHX) = -1
                  LABEL(IHINDX(2,NHX),NHX) = -1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
                end if
              END IF
            end if
          end if
        END DO
      END DO
*
* For general purposes: First type with 
*                       1-electron and 2-electron integrals:
      IFTONE = 1
      IFTTWO = NHX + 1
*
*. a+i alpha  a+j alpha  a k alpha  a l alpha, ordered so i.ge.j, k.ge.l
      DO IGAS = 1, NGAS
        DO JGAS = 1, IGAS 
          DO KGAS = 1, NGAS
            DO LGAS = 1, KGAS
              if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1))) then
                CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                IHX_LOC(IGAS+(1-1)*NGAS) = 1
                IHX_LOC(JGAS+(1-1)*NGAS) = IHX_LOC(JGAS+(1-1)*NGAS) + 1
                IHX_LOC(KGAS+(3-1)*NGAS) = 1
                IHX_LOC(LGAS+(3-1)*NGAS) = IHX_LOC(LGAS+(3-1)*NGAS) + 1
C               IHX(IGAS,1,NHX) = 1
C               IHX(JGAS,1,NHX) = IHX(JGAS,1,NHX)+1
* means IGAS=JGAS, so the number of creators of this type is 2!
C               IHX(KGAS,3,NHX) = 1
C               IHX(LGAS,3,NHX) = IHX(LGAS,3,NHX)+1
                CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                IF (INCLUDED.EQ.0) THEN 
                  WRITE(6,*) ' Hamiltonian type excluded '
                  CALL WRT_SPOX_TP(IHX_LOC,1)
                ELSE
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 2
                    IHINDX(3,NHX) = 3
                    IHINDX(4,NHX) = 4
                    SIGN_NHX(NHX) = ONE
                    LABEL(IHINDX(1,NHX),NHX) = 1
                    LABEL(IHINDX(2,NHX),NHX) = 1
                    LABEL(IHINDX(3,NHX),NHX) = 1
                    LABEL(IHINDX(4,NHX),NHX) = 1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 2
                    IOP_REO(3,NHX) = 3
                    IOP_REO(4,NHX) = 4
                    SIGN_OPREO(NHX) = ONE
                  end if
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
*
*. a+i beta a+j beta a k beta a l beta, ordered so i.ge.j, k.ge.l
      DO IGAS = 1, NGAS
        DO JGAS = 1, IGAS 
          DO KGAS = 1, NGAS
            DO LGAS = 1, KGAS
              if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &              NOBPT(IGAS).ge.2).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &              NOBPT(KGAS).ge.2)).or.
     &            ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                               MXNELGS(JGAS).ge.1).and.
     &             (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                               MXNELGS(LGAS).ge.1))) then
                CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                IHX_LOC(IGAS+(2-1)*NGAS) = 1
                IHX_LOC(JGAS+(2-1)*NGAS) = IHX_LOC(JGAS+(2-1)*NGAS) + 1
                IHX_LOC(KGAS+(4-1)*NGAS) = 1
                IHX_LOC(LGAS+(4-1)*NGAS) = IHX_LOC(LGAS+(4-1)*NGAS) + 1
                CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                IF (INCLUDED.EQ.0) THEN 
                  WRITE(6,*) ' Hamiltonian type excluded '
                  CALL WRT_SPOX_TP(IHX_LOC,1)
                ELSE
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 2
                    IHINDX(3,NHX) = 3
                    IHINDX(4,NHX) = 4
                    SIGN_NHX(NHX) = ONE
                    LABEL(IHINDX(1,NHX),NHX) = -1
                    LABEL(IHINDX(2,NHX),NHX) = -1
                    LABEL(IHINDX(3,NHX),NHX) = -1
                    LABEL(IHINDX(4,NHX),NHX) = -1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 2
                    IOP_REO(3,NHX) = 3
                    IOP_REO(4,NHX) = 4
                    SIGN_OPREO(NHX) = ONE
                  end if
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
*
* First type of integrals for 2-electron diagonal
      I2DIATP = NHX + 1
*
*. a+i alpha  a+j beta  a k beta  a l alpha
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          DO KGAS = 1, NGAS
            DO LGAS = 1, NGAS 
              if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &            MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &            NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &            NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1) then
                CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                IHX_LOC(IGAS+(1-1)*NGAS) = 1
                IHX_LOC(JGAS+(2-1)*NGAS) = 1
                IHX_LOC(KGAS+(4-1)*NGAS) = 1
                IHX_LOC(LGAS+(3-1)*NGAS) = 1
                CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                IF (INCLUDED.EQ.0) THEN 
                  WRITE(6,*) ' Hamiltonian type excluded '
                  CALL WRT_SPOX_TP(IHX_LOC,1)
                ELSE
                  NHX = NHX + 1
                  if (IFLAG.ne.1) then
                    CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  +a  a  +b  b
* so reorder accordingly:
                    IHINDX(1,NHX) = 1
                    IHINDX(2,NHX) = 3
                    IHINDX(3,NHX) = 4
                    IHINDX(4,NHX) = 2
                    SIGN_NHX(NHX) = -1.D0 * ONE
                    LABEL(IHINDX(1,NHX),NHX) = 1
                    LABEL(IHINDX(2,NHX),NHX) = -1
                    LABEL(IHINDX(3,NHX),NHX) = -1
                    LABEL(IHINDX(4,NHX),NHX) = 1
                    IOP_REO(1,NHX) = 1
                    IOP_REO(2,NHX) = 4
                    IOP_REO(3,NHX) = 3
                    IOP_REO(4,NHX) = 2
                    SIGN_OPREO(NHX) = -1.D0 * ONE
                  end if
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
*
      if (IHTYPE.ge.4) then
*
*. a+i alpha  a+j beta  a k alpha  a l beta
*   Generates type (ub|bu) of integrals.
*
        do IGAS = 1,NGAS,1
          do JGAS = 1,NGAS,1
            do KGAS = 1,NGAS,1
              do LGAS = 1,NGAS,1
                if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &              MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &              NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &              NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1) then
                  call isetvc(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = 1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = 1
                  call caab_included(IHX_LOC,INCLUDED)
                  if (INCLUDED.eq.0) then
                    write(6,*) 'Hamiltonian type excluded '
                    call wrt_spox_tp(IHX_LOC,1)
                  else
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      call icopve(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  +a  b  +b  a
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 2
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = -1
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 4
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 3
                      SIGN_OPREO(NHX) = ZERO
                    end if
                  end if
                end if
              end do
            end do
          end do
        end do
*
      end if
*
* Last type of non-spin-flip integrals
      ILTTWO = NHX
*
* one-electron spin orbit
*
      IF (IHTYPE.GT.2) THEN
*
* a+i alpha  j beta
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &          NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
              CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
              IHX_LOC(IGAS+(1-1)*NGAS) = 1
              IHX_LOC(JGAS+(4-1)*NGAS) = 1
              CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN 
                WRITE(6,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  SIGN_NHX(NHX) = ONE
                  LABEL(1,NHX) = 1
                  LABEL(2,NHX) = -1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
                end if
              END IF
            END IF
          END DO
        END DO
*
* a+i beta a j alpha 
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            if (MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &          NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1) then
              CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
              IHX_LOC(IGAS+(2-1)*NGAS) = 1
              IHX_LOC(JGAS+(3-1)*NGAS) = 1
              CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
              IF (INCLUDED.EQ.0) THEN 
                WRITE(6,*) ' Hamiltonian type excluded '
                CALL WRT_SPOX_TP(IHX_LOC,1)
              ELSE
                NHX = NHX + 1
                if (IFLAG.ne.1) then
                  CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
                  SIGN_NHX(NHX) = ONE
                  IHINDX(1,NHX) = 1
                  IHINDX(2,NHX) = 2
                  IHINDX(3,NHX) = 0
                  IHINDX(4,NHX) = 0
                  LABEL(1,NHX) = -1
                  LABEL(2,NHX) = 1
                  LABEL(3,NHX) = 0
                  LABEL(4,NHX) = 0
                  IOP_REO(1,NHX) = 1
                  IOP_REO(2,NHX) = 2
                  IOP_REO(3,NHX) = 0
                  IOP_REO(4,NHX) = 0
                  SIGN_OPREO(NHX) = ONE
                end if
              END IF
            END IF
          END DO
*    quaternionic matrix groups only!
        END DO
      END IF
*
* Two-electron terms with Delta mk .ne. 0
*
      if (IHTYPE.eq.3.or.IHTYPE.eq.5) THEN
*
**********************************************************
*  Spin-dependent spinor basis. E.g. DIRAC environment   *
*    (ub|uu), (bu|uu), (bu|bb), (ub|bb) integral classes *
**********************************************************
*
*. a+i alpha a+j alpha a k alpha a l beta, i .ge. j
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, NGAS 
                if ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &               MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &               NOBPT(IGAS).ge.2.and.NOBPT(KGAS).ge.1.and.
     &               NOBPT(LGAS).ge.1).or.
     &               (IGAS.ne.JGAS.and.
     &                MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &                MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1)) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(1-1)*NGAS) = IHX_LOC(JGAS+(1-1)*NGAS)+1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = 1
                  CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN 
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  +a  a  b  a+
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 2
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = 1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = -1 
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 3
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = -1.D0 * ONE
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
*. a+i beta a+j alpha a k alpha a l alpha, k .ge. l
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS 
                if ((KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &               MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &               NOBPT(KGAS).ge.2.and.NOBPT(IGAS).ge.1.and.
     &               NOBPT(JGAS).ge.1).or.
     &               (KGAS.ne.LGAS.and.
     &                MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &                MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1)) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = 1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(3-1)*NGAS) = IHX_LOC(LGAS+(3-1)*NGAS)+1
                  CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN 
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  +b  a  a  a+
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 2
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = -1
                      LABEL(IHINDX(2,NHX),NHX) = 1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = 1 
                      IOP_REO(1,NHX) = 3
                      IOP_REO(2,NHX) = 1
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = ONE
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
*. a+i beta  a+j beta  a k beta  a l alpha , i .ge. j
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, NGAS 
                if ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &               MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &               NOBPT(IGAS).ge.2.and.NOBPT(KGAS).ge.1.and.
     &               NOBPT(LGAS).ge.1).or.
     &               (IGAS.ne.JGAS.and.
     &                MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &                MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1)) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(2-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = IHX_LOC(JGAS+(2-1)*NGAS)+1
                  IHX_LOC(KGAS+(4-1)*NGAS) = 1
                  IHX_LOC(LGAS+(3-1)*NGAS) = 1
                  CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN 
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  +b  a  b  b+
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 4
                      IHINDX(3,NHX) = 3
                      IHINDX(4,NHX) = 2
C  ???                SIGN_NHX(NHX) = -1 * ONE
                      SIGN_NHX(NHX) = ONE
                      LABEL(IHINDX(1,NHX),NHX) = -1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = -1
                      LABEL(IHINDX(4,NHX),NHX) = 1 
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 3
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = ZERO
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
*. a+i alpha  a+j beta  a k beta  a l beta , k .ge. l
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS 
                if ((KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &               MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &               NOBPT(KGAS).ge.2.and.NOBPT(IGAS).ge.1.and.
     &               NOBPT(JGAS).ge.1).or.
     &               (KGAS.ne.LGAS.and.
     &                MXNELGS(IGAS).ge.1.and.MXNELGS(JGAS).ge.1.and.
     &                MXNELGS(KGAS).ge.1.and.MXNELGS(LGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1)) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = 1
                  IHX_LOC(KGAS+(4-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = IHX_LOC(LGAS+(4-1)*NGAS)+1
                  CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN 
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  +a  b  +b  b
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 3
                      IHINDX(3,NHX) = 4
                      IHINDX(4,NHX) = 2
                      SIGN_NHX(NHX) = ONE  !!??
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = -1
                      LABEL(IHINDX(4,NHX),NHX) = -1 
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 3
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = ZERO
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
      end if
*
**********************************************************
*  Spin-dependent spinor basis. E.g. DIRAC environment   *
*      (bu|bu), (ub|ub) integral classes                 *
**********************************************************
*
      if (IHTYPE.eq.4) then
*
*. a+i alpha  a+j alpha  a k beta  a l beta, i.ge.j, k.ge.l
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS
                if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1))) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(1-1)*NGAS) = 1
                  IHX_LOC(JGAS+(1-1)*NGAS) = IHX_LOC(JGAS+(1-1)*NGAS)+1
                  IHX_LOC(KGAS+(4-1)*NGAS) = 1
                  IHX_LOC(LGAS+(4-1)*NGAS) = IHX_LOC(LGAS+(4-1)*NGAS)+1
                  CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN 
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  +a  b  +a  b 
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 2
                      IHINDX(3,NHX) = 3
                      IHINDX(4,NHX) = 4
                      SIGN_NHX(NHX) = 1 * ONE
                      LABEL(IHINDX(1,NHX),NHX) = 1
                      LABEL(IHINDX(2,NHX),NHX) = 1
                      LABEL(IHINDX(3,NHX),NHX) = -1
                      LABEL(IHINDX(4,NHX),NHX) = -1 
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 3
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = -1.D0 * ONE
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
*
*. a+i beta  a+j beta  a k alpha  a l alpha, i.ge.j, k.ge.l
        DO IGAS = 1, NGAS
          DO JGAS = 1, IGAS
            DO KGAS = 1, NGAS
              DO LGAS = 1, KGAS
                if (((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.eq.JGAS.and.MXNELGS(IGAS).ge.2.and.
     &                NOBPT(IGAS).ge.2).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.eq.LGAS.and.MXNELGS(KGAS).ge.2.and.
     &                NOBPT(KGAS).ge.2)).or.
     &              ((IGAS.ne.JGAS.and.MXNELGS(IGAS).ge.1.and.
     &                NOBPT(IGAS).ge.1.and.NOBPT(JGAS).ge.1.and.
     &                                 MXNELGS(JGAS).ge.1).and.
     &               (KGAS.ne.LGAS.and.MXNELGS(KGAS).ge.1.and.
     &                NOBPT(KGAS).ge.1.and.NOBPT(LGAS).ge.1.and.
     &                                 MXNELGS(LGAS).ge.1))) then
                  CALL ISETVC(IHX_LOC,IZERO,4*NGAS)
                  IHX_LOC(IGAS+(2-1)*NGAS) = 1
                  IHX_LOC(JGAS+(2-1)*NGAS) = IHX_LOC(JGAS+(2-1)*NGAS)+1
                  IHX_LOC(KGAS+(3-1)*NGAS) = 1
                  IHX_LOC(LGAS+(3-1)*NGAS) = IHX_LOC(LGAS+(3-1)*NGAS)+1
                  CALL CAAB_INCLUDED(IHX_LOC,INCLUDED)
                  IF (INCLUDED.EQ.0) THEN 
                    WRITE(6,*) ' Hamiltonian type excluded '
                    CALL WRT_SPOX_TP(IHX_LOC,1)
                  ELSE
                    NHX = NHX + 1
                    if (IFLAG.ne.1) then
                      CALL ICOPVE(IHX_LOC,IHX(1,1,NHX),4*NGAS)
* This excitation corresponds to the integral
*  +b  a  +b  a
* so reorder accordingly:
                      IHINDX(1,NHX) = 1
                      IHINDX(2,NHX) = 2
                      IHINDX(3,NHX) = 3
                      IHINDX(4,NHX) = 4
                      SIGN_NHX(NHX) = 1 * ONE
                      LABEL(IHINDX(1,NHX),NHX) = -1
                      LABEL(IHINDX(2,NHX),NHX) = -1
                      LABEL(IHINDX(3,NHX),NHX) = 1
                      LABEL(IHINDX(4,NHX),NHX) = 1 
                      IOP_REO(1,NHX) = 1
                      IOP_REO(2,NHX) = 3
                      IOP_REO(3,NHX) = 2
                      IOP_REO(4,NHX) = 4
                      SIGN_OPREO(NHX) = -1.D0 * ONE
                    end if
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
      end if
*
      if (IFLAG.ne.1.and.NTEST.ge.50) then
        do I=1,NHX
          write(6,*) 'I,SIGN_NHX ',I,SIGN_NHX(I)
        end do
      end if
*
      if (IFLAG.ne.1.and.NTEST.ge.100) then
        WRITE(6,*) ' Information from  GET_HX_RELA '
        WRITE(6,*)
        write(6,*) ' Max. no of electrons per GAS:'
        call iwrtma(MXNELGS,1,NGAS,1,NGAS)
        WRITE(6,*) ' The excitation types '
        CALL WRT_SPOX_TP(IHX,NHX)
C            WRT_SPOX_TP(IEX_TP,NEX_TP)
        else
        WRITE(6,*) ' Number of excitation types ', NHX
      end if
*
      RETURN
      END 
