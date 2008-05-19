module aoresponse_main_m

  use kinds_m
! use memory_allocator
  use dirac_matrix_operations_m
  use dirac_interface_m
  use propint
  use pertdens
  use proptest
! use aoresponse_input_m

  implicit none

! public something

! private

contains

  subroutine test_matrix_stat()

    type(matrix) :: A, B, C

                            call matrix_stat('entering subroutine   ')
    call mat_init(A, 5, 5); call matrix_stat('call mat_init(A, 5, 5)')
    A = A                 ; call matrix_stat('A = A                 ')
    B = 2.0d0*A           ; call matrix_stat('B = 2.0d0*A           ')
    C = A + B             ; call matrix_stat('C = A + B             ')
    C = C/2.0d0           ; call matrix_stat('C = C/2.0d0           ')

    call test_matrix_stat_sub()
    call test_matrix_stat_sub()
    call test_matrix_stat_sub()

    call mat_free(A)      ; call matrix_stat('call mat_free(A)      ')
    B = 0                 ; call matrix_stat('B = 0                 ')
    C = 0                 ; call matrix_stat('C = 0                 ')
                            call matrix_stat('leaving subroutine    ')

  end subroutine test_matrix_stat

  subroutine test_matrix_stat_sub()

    type(matrix) :: A, B, C

                            call matrix_stat('entering subroutine   ')
    call mat_init(A, 5, 5); call matrix_stat('call mat_init(A, 5, 5)')
    A = A                 ; call matrix_stat('A = A                 ')
    B = 2.0d0*A           ; call matrix_stat('B = 2.0d0*A           ')
    C = A + B             ; call matrix_stat('C = A + B             ')
    C = C/2.0d0           ; call matrix_stat('C = C/2.0d0           ')


    call mat_free(A)      ; call matrix_stat('call mat_free(A)      ')
    B = 0                 ; call matrix_stat('B = 0                 ')
    C = 0                 ; call matrix_stat('C = 0                 ')
                            call matrix_stat('leaving subroutine    ')

  end subroutine test_matrix_stat_sub







  subroutine aoresponse_driver()

#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

    type(matrix)               :: S
    type(matrix)               :: C
    type(matrix)               :: D_mo
    type(matrix)               :: D
    type(matrix)               :: F

    real(kind_dp), allocatable :: freq(:)
    integer           :: nr_freq = 1

    allocate(freq(nr_freq))
    freq = 0.0d0

!!!!!!!!!!
    freq(1) = 0.00d0
    freq(1) = 0.02d0
!!!!!!!!!!

    write(*, *) 'first follows'
!   call matrixstats()
    call get_S(S)
    write(*, *) 'second follows'
!   call matrixstats()

    call get_C(C)
!   call matrixstats()
    call get_D_mo(D_mo)
!   call matrixstats()

    D = C*D_mo*dag(C)
    D_mo = 0
!   call matrixstats()
    C    = 0
!   call matrixstats()

    D = D/2.0d0
    D%is_density_matrix = .true.


    call get_F(F)

    call dipmom(S, D, F)
!   call polari(S, D, F, freq(1))

    write(*, *) 'compare: -4.828661928957 at 0.00'
    write(*, *) 'compare: -4.833040214626 at 0.02'

    S    = 0
    D    = 0
    F    = 0

    deallocate(freq)

  end subroutine aoresponse_driver

end module aoresponse_main_m
