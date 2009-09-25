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

      implicit none

#include "priunit.h"
#include "dcbqal.h"

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
      real(8) :: CMO(*)
      real(8) :: STRING(*)

      real(8), allocatable :: bvec(:)
      real(8), allocatable :: bmat(:)

      integer, allocatable :: orbital_rotation_indices(:, :)
      real(8), allocatable :: response_vector(:, :)
      integer              :: length
      character(2)         :: vector_type
      real(8)              :: ih_factor

      integer              :: nr_vectors, LUNIT_XVECTOR
      integer              :: ixvector, NBAS, NORBS, len
      integer              :: JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2, IVECTOR
      real(8)              :: FREQ1, FREQ2, ERGRSP, unused

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


      visual_cfg_irep_xvector(IXVECTOR) = JSYMOP1 - 1
      length  = LEN/NZ
      D_mo%irep   = JSYMOP1 - 1



    select case (typ1(1:1))
      case ('E')
        vector_type = 'pp'
      case ('P')
        vector_type = 'pn'
      case default
        call quit('unexpected typ1(1:1)')
    end select

    select case (typ1(2:2))
      case ('+')
        D_mo%ih_sym =  1
        ih_factor   =  1.0d0
      case ('-')
        D_mo%ih_sym = -1
        ih_factor   = -1.0d0
      case default
        call quit('unexpected typ1(2:2)')
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

!   allocate response vector
    call alloc(response_vector, length, nz)
    response_vector = 0.0d0

!   read response vector
    call readt(lunit_xvector, length*nz, response_vector)
    close(lunit_xvector, status = 'keep')

    call scatter_vector(length,                   &
                        orbital_rotation_indices, &
                        ih_factor,                &
                        response_vector,          &
                        D_mo%elms,                &
                        D_mo%irep)

    call dealloc(orbital_rotation_indices)
    call dealloc(response_vector)

    end subroutine
