MODULE mm_memory_manager_mod

! Module to record memory allocations.
! Note it only records total allocations, not
! peak allocation.

   USE mm_global_paras_mod

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: mm_allocate,        &
             mm_init_mem_man,    &
             mm_close_mem_man

   INTERFACE mm_allocate
      MODULE PROCEDURE allocate_1
      MODULE PROCEDURE allocate_2
      MODULE PROCEDURE allocate_3
      MODULE PROCEDURE allocate_4
      MODULE PROCEDURE allocate_5
   END INTERFACE

   INTEGER(INTK), PARAMETER :: Frealk = 8
   INTEGER(INTK), PARAMETER :: Fraw_mm_paras = 92
   INTEGER(INTK), PARAMETER :: Fbox_mm_data = 4
   INTEGER(INTK), PARAMETER :: Fbox_mm_paras = 76

   INTEGER(INTK), PARAMETER :: MBytes = 1048576
   REAL(REALK), SAVE :: mem_use(NMEMDIVS)

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_close_mem_man

      IMPLICIT NONE
      INTEGER(INTK) :: i
      LOGICAL, PARAMETER :: PRINT_DEBUG_MEMORY=.FALSE.

      IF (PRINT_DEBUG_MEMORY) THEN 
         mem_use(:) = mem_use(:)/Mbytes
         
         WRITE(LUPRI,'(/,A)')      "  Cumulative memory allocation / MB"
         WRITE(LUPRI,'(A,/)')      "  ---------------------------------"
         DO i = 1, NMEMDIVS
            WRITE(LUPRI,'(3X,A,"  =",F11.2)') NSPACE(i), mem_use(i)
         END DO
         WRITE(LUPRI,*)
      END IF

   END SUBROUTINE mm_close_mem_man

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_mem_man

      IMPLICIT NONE
      mem_use = zero

   END SUBROUTINE mm_init_mem_man

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_1(name,array,a)

      IMPLICIT NONE
      REAL(REALK),   POINTER    :: array(:) 
      INTEGER(INTK), INTENT(IN) :: a, name
      INTEGER(INTK) :: error

      ALLOCATE(array(a),STAT=error)
      IF (error /= 0) CALL QUIT('allocation error in allocate_1')
      mem_use(name) = mem_use(name) + a*Frealk

   END SUBROUTINE allocate_1

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_2(name,array,a,b)

      IMPLICIT NONE
      REAL(REALK),   POINTER    :: array(:,:) 
      INTEGER(INTK), INTENT(IN) :: a,b, name
      INTEGER(INTK) :: error

      ALLOCATE(array(a,b),STAT=error)
      IF (error /= 0) CALL QUIT('allocation error in allocate_2')
      mem_use(name) = mem_use(name) + a*b*Frealk

   END SUBROUTINE allocate_2

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_3(name,array,a)

      IMPLICIT NONE
      TYPE(raw_mm_paras), POINTER    :: array(:) 
      INTEGER(INTK),      INTENT(IN) :: a, name
      INTEGER(INTK) :: error

      ALLOCATE(array(a),STAT=error)
      IF (error /= 0) CALL QUIT('allocation error in allocate_3')
      mem_use(name) = mem_use(name) + a*Fraw_mm_paras

   END SUBROUTINE allocate_3

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_4(name,array,a)

      IMPLICIT NONE
      TYPE(box_mm_data), POINTER    :: array(:) 
      INTEGER(INTK),     INTENT(IN) :: a, name
      INTEGER(INTK) :: error

      ALLOCATE(array(a),STAT=error)
      IF (error /= 0) CALL QUIT('allocation error in allocate_4')
      mem_use(name) = mem_use(name) + a*Fbox_mm_data

   END SUBROUTINE allocate_4

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_5(name,array,a)

      IMPLICIT NONE
      TYPE(box_mm_paras), POINTER    :: array(:) 
      INTEGER(INTK),      INTENT(IN) :: a, name
      INTEGER(INTK) :: error

      ALLOCATE(array(a),STAT=error)
      IF (error /= 0) CALL QUIT('allocation error in allocate_5')
      mem_use(name) = mem_use(name) + a*Fbox_mm_paras

   END SUBROUTINE allocate_5

!-------------------------------------------------------------------------------

END MODULE mm_memory_manager_mod
