MODULE adc_fano_matmul

CONTAINS

  SUBROUTINE fano_create_t_mom(fn_adcdiag,fn_adcmat,ladc,nr_in,nr_fin, &
                               fhmat,fhdiag,intbuf,e2h1p,t_mom,iw)
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
    INTEGER, INTENT(IN)                             :: intbuf
    INTEGER, INTENT(IN)                             :: fhmat, fhdiag
    REAL(8), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: e2h1p
    REAL(8), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: t_mom
    INTEGER, INTENT(IN)                             :: iw

!
! Local Variables
!
    REAL(8), ALLOCATABLE, DIMENSION(:)    :: diag
    INTEGER                               :: ixx
    INTEGER                               :: nbufs
    INTEGER                               :: nact
    REAL(8), ALLOCATABLE, DIMENSION(:)    :: buf
    INTEGER, ALLOCATABLE, DIMENSION(:)    :: ioi, ioj
    INTEGER                               :: idummy
    REAL(8)                               :: buffer
    INTEGER                               :: irow, jcol
    INTEGER                               :: i, j, k
    INTEGER                               :: counter
    REAL(8), ALLOCATABLE, DIMENSION(:)    :: vec_start
    REAL(8), ALLOCATABLE, DIMENSION(:)    :: vec1, vec2


!
! Code
!

    CALL alloc(buf,intbuf)
    CALL alloc(ioi,intbuf)
    CALL alloc(ioj,intbuf)

    CALL alloc(vec_start,ladc) !to be deallocated
    CALL alloc(vec1,ladc)
    CALL alloc(vec2,ladc)

    CALL alloc(e2h1p,nr_in_evecs*nr_fin)
    CALL alloc(t_mom,nr_in_evecs*nr_fin)

    vec_start = 0
    counter   = 0

    OPEN(fhmat,FILE=fn_adcmat,FORM='unformatted',STATUS='unknown')
    OPEN(fhdiag,FILE=fn_adcdiag,FORM='unformatted',STATUS='unknown')

    DO i = 1, nr_in_evecs
      vec1 = 0
      vec_start(1:nr_in) = in_evecs(1:nr_in,i)

      CALL alloc(diag,ladc)
      REWIND(fhdiag)
      REWIND(fhmat)
      READ(fhdiag) (diag(ixx), ixx=1, ladc), nbufs
      DO j = 1, ladc
        vec1(j) = vec1(j) + diag(j) * vec_start(j)
      END DO
      CALL dealloc(diag)

      DO k = 1, nbufs
        READ(fhmat) (buf(ixx),ixx=1,intbuf), &
                    (ioi(ixx),ixx=1,intbuf), &
                    (ioj(ixx),ixx=1,intbuf), &
                    nact, idummy
        DO j = 1, nact
          irow   = ioi(j)
          jcol   = ioj(j)
          buffer = buf(j)

          vec1(irow) = vec1(irow) + buffer * vec_start(jcol)
          vec1(jcol) = vec1(jcol) + buffer * vec_start(irow)
        END DO
      END DO
!
! Now we have the initial state vector multiplied with the sorted ADC matrix.
! We will further multiply this vector to each and every final state vector.
!
      DO j = 1, nr_fin
        vec2 = 0
        counter = counter + 1
        vec2(ladc-nr_fin+1:ladc) = fin_evecs(1:nr_fin,j)
        e2h1p(counter) = fin_energies(j)
        t_mom(counter) = DOT_PRODUCT(vec1,vec2)
      END DO
    END DO

    CALL dealloc(buf)
    CALL dealloc(ioi)
    CALL dealloc(ioj)

    CALL dealloc(vec_start)
    CALL dealloc(vec1)
    CALL dealloc(vec2)

    CALL dealloc(in_evecs)
    CALL dealloc(fin_evecs)
    CALL dealloc(fin_energies)

!
! To be removed as soon as Dirac has its own Stieltjes routine
!
    CALL dealloc(e2h1p)
    CALL dealloc(t_mom)

  END SUBROUTINE fano_create_t_mom

END MODULE adc_fano_matmul
