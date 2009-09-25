

      SUBROUTINE SCATTER_BVEC(BVEC, &
     &                        ICOUNTER, &
     &                        TYP, &
     &                        BMAT, &
     &                        IZ, &
     &                        IFSYM1, &
     &                        IFSYM2, &
     &                        STRING)
!=======================================================================
!     scatter solution vector BVEC onto BMAT
!-----------------------------------------------------------------------
!                                                           radovan bast
!                                              last revision: april 2007
!=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"

      CHARACTER TYP*2
      DIMENSION BVEC(*), &
     &          BMAT(NORBT,NORBT,NZ), &
     &          STRING(*), &
     &          NEI(2), &
     &          NES(2), &
     &          NP(2)

      DO I = 1,2
        NEI(I) = NISH(I)
        NES(I) = NESH(I) - NISH(I)
        NP(I)  = NPSH(I)
      ENDDO

      IOFF = (NESH(1) + NPSH(1))*(IFSYM1 - 1) + NP(IFSYM1)

      II1 = IOFF + 1
      II2 = II1 + NEI(IFSYM1) - 1
      IS1 = 1
      IS1 = IS1 + (NESH(1) + NPSH(1))*(IFSYM2 - 1)
      IS2 = IS1 + NP(IFSYM2)  - 1
      IS3 = IS2 + NEI(IFSYM2) + 1
      IS4 = IS3 + NES(IFSYM2) - 1

      DO II = II1,II2

!fixme
!       factor = string(ii - ioff + nesh(1)*(ifsym1-1))
        factor = 1.0d0
        IF(TYP(1:1) .EQ. 'P') THEN
          DO IS = IS1,IS2
            BMAT(IS,II,IZ) = BVEC(ICOUNTER)*FACTOR
            ICOUNTER = ICOUNTER + 1
          ENDDO
        ENDIF
        IF(TYP(1:1) .EQ. 'E') THEN
          DO IS = IS3,IS4
            BMAT(IS,II,IZ) = BVEC(ICOUNTER)*FACTOR
            ICOUNTER = ICOUNTER + 1
          ENDDO
        ENDIF

      ENDDO

      END
      SUBROUTINE SETUP_FO_DMAT_AO(CMO, &
     &                            IXVECTOR, &
     &                            D, &
     &                            STRING, D_mo)
!=======================================================================
!     description in here ...
!-----------------------------------------------------------------------
!                                                           radovan bast
!                                              last revision: april 2007
!=======================================================================

      use memory_allocator
      use visual_cfg
      use matrix_defop
      use dirac_interface

#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "pgroup.h"
#include "dcbxrs.h"

      CHARACTER LAB1*16,&
     &          LAB2*16,&
     &          TYP1*2,&
     &          TYP2*2,&
     &          LBUF*16,&
     &          PBUF*1
      LOGICAL   FILE_EXISTS
      type(matrix) :: D, D_mo
      DIMENSION CMO(*),&
     &          STRING(*)

      real(8), allocatable :: bvec(:)
      real(8), allocatable :: bmat(:)

      integer, allocatable :: orbital_rotation_indices(:, :)
      real(8), allocatable :: response_vector(:, :, :)
      integer              :: length
      character(2)         :: vector_type
      real(8)              :: ih_factor

      INQUIRE(FILE = visual_cfg_xvector_file(IXVECTOR),EXIST = FILE_EXISTS)
      IF(.NOT. FILE_EXISTS) THEN
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A,A6,A)') 'solution vector file ', &
     &    visual_cfg_XVECTOR_FILE(IXVECTOR),' not found'
        CALL QUIT('ERROR in VISUAL')
      ENDIF
      LUNIT_XVECTOR = 17
      OPEN(LUNIT_XVECTOR, &
     &     FILE = visual_cfg_xvector_file(IXVECTOR), &
     &     STATUS = 'OLD', &
     &     ACCESS = 'SEQUENTIAL', &
     &     FORM = 'UNFORMATTED')

!     get number of vectors on file

      NR_VECTORS = 0
      REWIND LUNIT_XVECTOR
 1    CONTINUE
      READ(LUNIT_XVECTOR) LAB1
      IF(LAB1 .EQ. 'END_OF_THIS_FILE') THEN
        GOTO 2
      ELSE
        NR_VECTORS = NR_VECTORS + 1
        READ(LUNIT_XVECTOR)
        GOTO 1
      ENDIF
 2    CONTINUE

!     errors

      IF(visual_cfg_NR_xVECTORS .EQ. 0) THEN
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'NR_VECTORS .EQ. 0'
        CALL QUIT('ERROR in VISUAL')
      ENDIF

      IF(visual_cfg_IWHICH_XVECTOR(IXVECTOR) .GT. NR_VECTORS) THEN
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'IWHICH_XVECTOR(IXVECTOR) .GT. NR_VECTORS'
        CALL QUIT('ERROR in VISUAL')
      ENDIF

      IF(visual_cfg_IWHICH_XVECTOR(IXVECTOR) .LT. 1) THEN
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'IWHICH_XVECTOR(IXVECTOR) .LT. 1'
        CALL QUIT('ERROR in VISUAL')
      ENDIF

!     move to the right record

      IVECTOR = 0
      REWIND LUNIT_XVECTOR
 3    CONTINUE
      READ(LUNIT_XVECTOR)
      IVECTOR = IVECTOR + 1
      IF(IVECTOR .EQ. visual_cfg_IWHICH_XVECTOR(IXVECTOR)) THEN
        BACKSPACE(LUNIT_XVECTOR)
        GOTO 4
      ELSE
        READ(LUNIT_XVECTOR)
        GOTO 3
      ENDIF
 4    CONTINUE

!     get type and dimensions of record IWHICH_XVECTOR(IXVECTOR)

      READ(LUNIT_XVECTOR) LAB1, LAB2, TYP1, TYP2, FREQ1, FREQ2, &
     &                    JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2, &
     &                    UNUSED, LEN, INTFLG, ERGRSP, NBAS, NORBS

      IF(TYP1(2:2) .EQ. '+') IHERMITICITY =  1
      IF(TYP1(2:2) .EQ. '-') IHERMITICITY = -1

      JSYMOP                 = JSYMOP1
      visual_cfg_irep_xvector(IXVECTOR) = JSYMOP - 1
      IREPOP = JSYMOP - 1
      JOPSY  = JBTOF(visual_cfg_IREP_XVECTOR(IXVECTOR),1)
      NBDIM  = LEN/NZ
      length = nbdim 
      D_mo%irep   = jsymop - 1
      D_mo%ih_sym = ihermiticity

    select case (typ1(1:1))
      case ('E')
        vector_type = 'pp'
      case ('P')
        vector_type = 'pn'
      case default
        call quit('unexpected typ1(1:1)')
    end select

    call alloc(orbital_rotation_indices, 2, length)

    select case (vector_type)
      case ('pp')
        call set_orbital_rotation_indices(D_mo%irep, &
                                          length,    &
                                          2,         &
                                          3,         &
                                          orbital_rotation_indices)
      case ('pn')
        call set_orbital_rotation_indices(D_mo%irep, &
                                          length,    &
                                          2,         &
                                          1,         &
                                          orbital_rotation_indices)
      case default
        call quit('unexpected vector_type')
    end select

 



!     allocate memory

      call alloc(bvec, nbdim*nz)
      call alloc(bmat, n2orbxq)

!     read vector

      CALL DZERO(BVEC,NBDIM*NZ)
      CALL READT(LUNIT_XVECTOR,NBDIM*NZ,BVEC)

      CLOSE(LUNIT_XVECTOR,STATUS = 'KEEP')

!     scatter vector onto MO density matrix


    select case (D_mo%ih_sym)
      case (+1)
        ih_factor =  1.0d0
      case (-1)
        ih_factor = -1.0d0
      case default
        call quit('unexpected D_mo%ih_sym')
    end select

    call scatter_vector(length,                   &
                        orbital_rotation_indices, &
                        ih_factor,                &
                        bvec,                     &
                        D_mo%elms,                &
                        D_mo%irep)

    call dealloc(orbital_rotation_indices)


!     CALL DZERO(BMAT,N2ORBXQ)
!     ICOUNTER = 1

!     my static electric perturbations are ungerade
!               magnetic                     gerade
!
!     i don't know codewise how to extract gerade/ungerade
!     from the general operator irrep
!
!     with frequency-dependent magnetic perturbations
!     the code will fail when calculating with symmetry
!     but for the moment it does its job ...

!     IGERADE = -IHERMITICITY
!     IF(DABS(FREQ1) .GT. 1.0D-9) IGERADE = -1

!     DO IZ = 1,NZ
!       IF(NFSYM .EQ. 1) THEN
!         CALL SCATTER_BVEC(BVEC, &
!    &                      ICOUNTER, &
!    &                      TYP1, &
!    &                      BMAT, &
!    &                      IZ, &
!    &                      1,1, &
!    &                      STRING)
!       ELSE
!         IF(IGERADE .EQ. -1) THEN
!           CALL SCATTER_BVEC(BVEC, &
!    &                        ICOUNTER, &
!    &                        TYP1, &
!    &                        BMAT, &
!    &                        IZ, &
!    &                        1,2, &
!    &                        STRING)
!           CALL SCATTER_BVEC(BVEC, &
!    &                        ICOUNTER, &
!    &                        TYP1, &
!    &                        BMAT, &
!    &                        IZ, &
!    &                        2,1, &
!    &                        STRING)
!         ELSE
!           CALL SCATTER_BVEC(BVEC, &
!    &                        ICOUNTER, &
!    &                        TYP1, &
!    &                        BMAT, &
!    &                        IZ, &
!    &                        1,1, &
!    &                        STRING)
!           CALL SCATTER_BVEC(BVEC, &
!    &                        ICOUNTER, &
!    &                        TYP1, &
!    &                        BMAT, &
!    &                        IZ, &
!    &                        2,2, &
!    &                        STRING)
!         ENDIF
!       ENDIF
!     ENDDO
!     CALL DSCAL(N2ORBXQ,-D1,BMAT,1)
!     call dealloc(bvec)

!     DO IZ = 1,NZ
!       IOFF = NORBT*NORBT*(IZ - 1)
!       CALL FULMAT('S',NORBT,NORBT,BMAT(1+IOFF)) !symmetrize
!     ENDDO

!     get AO density matrix

!     CALL DZERO(D%elms,N2BBASXQ)
!         IF(TYP1(1:1) .EQ. 'E') THEN
!       CALL DMOLI3(0,IHERMITICITY,D%elms,CMO,BMAT)
!     ELSEIF(TYP1(1:1) .EQ. 'P') THEN
!       CALL DMOLI3(1,IHERMITICITY,D%elms,CMO,BMAT)
!     ENDIF
!     call dealloc(bmat)

!     IF(NZ .LT. 4) THEN
!       DO IZ = 1,NZ
!         IQ = IPQTOQ(IZ,visual_cfg_IREP_XVECTOR(IXVECTOR))
!         IOFF = NTBAS(0)*NTBAS(0)*(IZ-1)
!         CALL Q2BPHASE('D',IQ,1,D%elms(1+IOFF))
!       ENDDO
!     ENDIF

!     CALL DSCAL(N2BBASXQ,D2,D%elms,1)

      END
