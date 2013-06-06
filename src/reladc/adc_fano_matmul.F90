MODULE adc_fano_matmul

CONTAINS

  SUBROUTINE fano_create_t_mom(fn_adcdiag,fn_adcmat,ladc,nr_in,nr_fin, &
                               e2h1p,t_mom,iw)
!
! Purpose: 
!

    use adc_fano_exchange
    use memory_allocator

    IMPLICIT NONE

!
! Calling Variables
!
    CHARACTER(6), INTENT(IN)                        :: fn_adcdiag,fn_adcmat
    INTEGER, INTENT(IN)                             :: ladc
    INTEGER, INTENT(IN)                             :: nr_in, nr_fin
    REAL(8), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: e2h1p
    REAL(8), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: t_mom
    INTEGER, INTENT(IN)                             :: iw

!
! Local Variables
!

!
! Code
!

    CALL dealloc(in_evecs)
    CALL dealloc(fin_evecs)
    CALL dealloc(fin_energies)

  END SUBROUTINE fano_create_t_mom

END MODULE adc_fano_matmul
