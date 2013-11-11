MODULE adc_fano_matmul

CONTAINS

  SUBROUTINE fano_create_t_mom(fn_adcdiag,fn_adcmat,ladc,nr_in,nr_fin, &
                               fhmat,fhdiag,intbuf,e2h1p,t_mom,nbufs,iw)
!
! Purpose: Do the final multiplication of the initial wave function
!          with the coupling block and the final state wavefunctions.
!

    use adc_fano_exchange
    use memory_allocator
    use codata
    use adc_cfg

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
    INTEGER, INTENT(IN)                             :: nbufs
    INTEGER, INTENT(IN)                             :: iw

!
! Local Variables
!
    REAL(8), ALLOCATABLE, DIMENSION(:)    :: diag
    INTEGER                               :: ixx
!    INTEGER                               :: nbufs
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

    CHARACTER(3), PARAMETER               :: nmspec = 'FPS'
    CHARACTER(11)                         :: fn_trans
    INTEGER                               :: fh_trans
    INTEGER                               :: ierror
    INTEGER                               :: channel
    INTEGER                               :: first, last, off


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

    OPEN(fhmat,FILE=fn_adcmat,FORM='unformatted',STATUS='unknown')
!    OPEN(fhdiag,FILE=fn_adcdiag,FORM='unformatted',STATUS='unknown')

    DO i = 1, nr_in_evecs
      counter = 0
      vec1 = 0
      vec_start(1:nr_in) = in_evecs(1:nr_in,i)
!      WRITE(iw,*) 'starting vector'
!      WRITE(iw,'(1X,ES14.7)') (vec_start(ixx), ixx=1,ladc)

! Choose file name of the pseudo-spectrum
      IF(i.GT.9) THEN
        WRITE(fn_trans,'(a3,a1,i2)') nmspec,'.',i
      ELSE
        WRITE(fn_trans,'(a3,a2,i1)') nmspec,'.0',i
      END IF


!      CALL alloc(diag,ladc)
!      REWIND(fhdiag)
      REWIND(fhmat)
!      READ(fhdiag) (diag(ixx), ixx=1, ladc), nbufs
!      DO j = 1, ladc
!        vec1(j) = vec1(j) + diag(j) * vec_start(j)
!      END DO
!      CALL dealloc(diag)
!      !WRITE(iw,*) 'Vector with diagonal'
!      !WRITE(iw,'(1X,ES14.7)') (vec1(ixx), ixx=1,ladc)
!      WRITE(iw,*) 'nbufs =', nbufs

      DO k = 1, nbufs
        READ(fhmat) (buf(ixx),ixx=1,intbuf), &
                    (ioi(ixx),ixx=1,intbuf), &
                    (ioj(ixx),ixx=1,intbuf), &
                    nact, idummy
!        WRITE(iw,*) 'nact = ', nact
        DO j = 1, nact
          irow   = ioi(j)
          jcol   = ioj(j)
          buffer = buf(j)
!          WRITE(iw,*) 'buffer = ', buffer

          vec1(irow) = vec1(irow) + buffer * vec_start(jcol)
          vec1(jcol) = vec1(jcol) + buffer * vec_start(irow)
        END DO
      END DO
!      WRITE(iw,*) 'Vector with total sorted matrix'
!      WRITE(iw,'(1X,ES14.7)') (vec1(ixx), ixx=1,ladc)
      
!
! Now we have the initial state vector multiplied with the sorted ADC matrix.
! We will further multiply this vector to each and every final state vector.
!
      DO j = 1, nr_fin
        vec2 = 0
        counter = counter + 1
        vec2(ladc-nr_fin+1:ladc) = fin_evecs(1:nr_fin,j)
        !WRITE(iw,*) 'Final state vector'
        !WRITE(iw,'(1X,ES14.7)') (vec2(ixx), ixx=1,ladc)
        e2h1p(counter) = fin_energies(j)
        ! This part is ony necessary, when the initial and final states
        ! are non-orthogonal, which they are not at the moment!
        !t_mom(counter) = DOT_PRODUCT(vec1,vec2) &
        !                 - DOT_PRODUCT(vec_start,vec2) * in_energy(i)
        t_mom(counter) = DOT_PRODUCT(vec1,vec2)
      END DO

      WRITE(iw,*)
      WRITE(iw,*) '-------------------------------------------------------------'
      WRITE(iw,'(1X,A45,I4)') 'Pseudo-spectrum for initial state number ', i
      WRITE(iw,'(1X,A10,F10.6,A3,A20,F10.6)') 'energy: ', in_energy(1,i),&
                                              'eV', 'polestrength: ', in_energy(2,i)
      WRITE(iw,*) '-------------------------------------------------------------'
      WRITE(iw,*)
      WRITE(iw,'(1X,2A25)') 'Energies', 'Transition moments'
      WRITE(iw,'(1X,2ES25.15)') (e2h1p(ixx), 2*pi*t_mom(ixx)*t_mom(ixx), ixx=1,counter)

    END DO


!
! To be removed as soon as Dirac has its own Stieltjes routine
!
    fh_trans = fhdiag
    OPEN(fh_trans,FILE=fn_trans,STATUS='UNKNOWN',FORM='FORMATTED',    &
         ACTION='WRITE',iostat=ierror)
    WRITE(fh_trans,'(1X,2ES25.15)') (e2h1p(ixx), 2*pi*t_mom(ixx)*t_mom(ixx), ixx=1,counter)
    CLOSE(fh_trans)

!
! Starting the partial decay widths
!
!

    WRITE(iw,*) 
    WRITE(iw,*) ' Starting creation of partial decay width pseudo-spectra'
    WRITE(iw,*)

    DO i = 1, nr_in_evecs
      last = 0
      vec_start(1:nr_in) = in_evecs(1:nr_in,i)
!      WRITE(iw,*) 'starting vector'
!      WRITE(iw,'(1X,ES14.7)') (vec_start(ixx), ixx=1,ladc)

      DO channel = 1, reladc_fano_nrgroups
        counter = 0
        vec1  = 0
        first = last + 1
        last  = first + nr2h1p(channel) - 1

        REWIND(fhmat)

! Choose file name of the pseudo-spectrum
        IF(i.GT.9) THEN
          WRITE(fn_trans,'(a3,a1,i2,a1,a4)') nmspec,'.',i,'.',reladc_fano_labels(channel)
        ELSE
          WRITE(fn_trans,'(a3,a2,i1,a1,a4)') nmspec,'.0',i,'.',reladc_fano_labels(channel)
        END IF


        DO k = 1, nbufs
          READ(fhmat) (buf(ixx),ixx=1,intbuf), &
                      (ioi(ixx),ixx=1,intbuf), &
                      (ioj(ixx),ixx=1,intbuf), &
                      nact, idummy
!          WRITE(iw,*) 'nact = ', nact
          DO j = 1, nact
            irow   = ioi(j)
            jcol   = ioj(j)
            buffer = buf(j)
!            WRITE(iw,*) 'buffer = ', buffer

            IF (first <= jcol .AND. jcol <= last) THEN
              vec1(irow) = vec1(irow) + buffer * vec_start(jcol)
            END IF
            IF (first <= irow .AND. irow <= last) THEN
              vec1(jcol) = vec1(jcol) + buffer * vec_start(irow)
            END IF
          END DO
        END DO
!        IF (channel == 2) THEN
!          WRITE(iw,*) 'Vector with total sorted matrix'
!          WRITE(iw,'(1X,ES14.7)') (vec1(ixx), ixx=1,ladc)
!        END IF
      
!
! Now we have the initial state vector multiplied with the sorted ADC matrix.
! We will further multiply this vector to each and every final state vector.
!
        DO j = 1, nr_fin
          vec2 = 0
          counter = counter + 1
          vec2(ladc-nr_fin+1:ladc) = fin_evecs(1:nr_fin,j)
          !IF (j == 1) THEN
          !  WRITE(iw,*) 'Final state vector'
          !  WRITE(iw,'(1X,ES14.7)') (vec2(ixx), ixx=1,ladc)
          !END IF
          e2h1p(counter) = fin_energies(j)
          ! This part is ony necessary, when the initial and final states
          ! are non-orthogonal, which they are not at the moment!
          !t_mom(counter) = DOT_PRODUCT(vec1,vec2) &
          !                 - DOT_PRODUCT(vec_start,vec2) * in_energy(i)
          t_mom(counter) = DOT_PRODUCT(vec1,vec2)
!          WRITE(iw,*) 't_mom: ', t_mom(counter)
        END DO
        
        WRITE(iw,*)
        WRITE(iw,*) 'name: ', fn_trans
        WRITE(iw,*) '-------------------------------------------------------------'
        WRITE(iw,'(1X,A45,I4)') 'Pseudo-spectrum for initial state number ', i
        WRITE(iw,'(1X,A10,F10.6,A3,A20,F10.6)') 'energy: ', in_energy(1,i),&
                                                'eV', 'polestrength: ', in_energy(2,i)
        WRITE(iw,*) '-------------------------------------------------------------'
        WRITE(iw,*)
        WRITE(iw,'(1X,2A25)') 'Energies', 'Transition moments'
        WRITE(iw,'(1X,2ES25.15)') (e2h1p(ixx), 2*pi*t_mom(ixx)*t_mom(ixx), ixx=1,counter)
        WRITE(iw,*)

!
! This procedure only makes sense, when those final vectors not
! associated with the chosen final states are not considered in
! the calculation of the transition moments.
!
!
! To be removed as soon as Dirac has its own Stieltjes routine
!
        OPEN(fh_trans,FILE=fn_trans,STATUS='UNKNOWN',FORM='FORMATTED',    &
             ACTION='WRITE',iostat=ierror)
        WRITE(fh_trans,'(1X,2ES25.15)') (e2h1p(ixx), 2*pi*t_mom(ixx)*t_mom(ixx), ixx=1,counter)
        CLOSE(fh_trans)

      END DO ! end of loop over channels
    END DO   ! end of loop over initial states

    CALL dealloc(buf)
    CALL dealloc(ioi)
    CALL dealloc(ioj)

    CALL dealloc(vec_start)
    CALL dealloc(vec1)
    CALL dealloc(vec2)

    CALL dealloc(in_evecs)
    CALL dealloc(in_energy)
    CALL dealloc(fin_evecs)
    CALL dealloc(fin_energies)


    CALL dealloc(e2h1p)
    CALL dealloc(t_mom)

  END SUBROUTINE fano_create_t_mom



!
!-----------------------------------------------------------------
!
  SUBROUTINE cfano_create_t_mom(fn_adcdiag,fn_adcmat,ladc,nr_in,nr_fin, &
                               fhmat,fhdiag,intbuf,e2h1p,t_mom,nbufs,iw)
!
! Purpose: Do the final multiplication of the initial wave function
!          with the coupling block and the final state wavefunctions.
!

    use adc_fano_exchange
    use memory_allocator
    use codata

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
    INTEGER, INTENT(IN)                             :: nbufs
    INTEGER, INTENT(IN)                             :: iw

!
! Local Variables
!
    COMPLEX(8), ALLOCATABLE, DIMENSION(:)    :: diag
    INTEGER                               :: ixx
!    INTEGER                               :: nbufs
    INTEGER                               :: nact
    COMPLEX(8), ALLOCATABLE, DIMENSION(:)    :: buf
    INTEGER, ALLOCATABLE, DIMENSION(:)    :: ioi, ioj
    INTEGER                               :: idummy
    COMPLEX(8)                               :: buffer
    INTEGER                               :: irow, jcol
    INTEGER                               :: i, j, k
    INTEGER                               :: counter
    COMPLEX(8), ALLOCATABLE, DIMENSION(:)    :: vec_start
    COMPLEX(8), ALLOCATABLE, DIMENSION(:)    :: vec1, vec2

    CHARACTER(6),PARAMETER                :: fn_trans = 'FPSPEC'
    INTEGER                               :: fh_trans
    INTEGER                               :: ierror


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
!    OPEN(fhdiag,FILE=fn_adcdiag,FORM='unformatted',STATUS='unknown')

    DO i = 1, nr_in_evecs
      vec1 = 0
      vec_start(1:nr_in) = cin_evecs(1:nr_in,i)
!      WRITE(iw,*) 'starting vector'
!      WRITE(iw,'(1X,ES14.7)') (vec_start(ixx), ixx=1,ladc)

!      CALL alloc(diag,ladc)
!      REWIND(fhdiag)
      REWIND(fhmat)
!      READ(fhdiag) (diag(ixx), ixx=1, ladc), nbufs
!      DO j = 1, ladc
!        vec1(j) = vec1(j) + diag(j) * vec_start(j)
!      END DO
!      CALL dealloc(diag)
!      !WRITE(iw,*) 'Vector with diagonal'
!      !WRITE(iw,'(1X,ES14.7)') (vec1(ixx), ixx=1,ladc)
!      WRITE(iw,*) 'nbufs =', nbufs

      DO k = 1, nbufs
        READ(fhmat) (buf(ixx),ixx=1,intbuf), &
                    (ioi(ixx),ixx=1,intbuf), &
                    (ioj(ixx),ixx=1,intbuf), &
                    nact, idummy
!        WRITE(iw,*) 'nact = ', nact
        DO j = 1, nact
          irow   = ioi(j)
          jcol   = ioj(j)
          buffer = buf(j) ! lower triangle is stored
!          WRITE(iw,*) 'buffer = ', buffer

          vec1(irow) = vec1(irow) + CONJG(buffer) * vec_start(jcol)
          vec1(jcol) = vec1(jcol) + buffer * vec_start(irow)
        END DO
      END DO
!      WRITE(iw,*) 'Vector with total sorted matrix'
!      WRITE(iw,'(1X,ES14.7)') (vec1(ixx), ixx=1,ladc)
      
!
! Now we have the initial state vector multiplied with the sorted ADC matrix.
! We will further multiply this vector to each and every final state vector.
!
      DO j = 1, nr_fin
        vec2 = 0
        counter = counter + 1
        vec2(ladc-nr_fin+1:ladc) = cfin_evecs(1:nr_fin,j)
        !WRITE(iw,*) 'Final state vector'
        !WRITE(iw,'(1X,ES14.7)') (vec2(ixx), ixx=1,ladc)
        e2h1p(counter) = fin_energies(j)
        ! This part is ony necessary, when the initial and final states
        ! are non-orthogonal, which they are not at the moment!
        !t_mom(counter) = DOT_PRODUCT(vec1,vec2) &
        !                 - DOT_PRODUCT(vec_start,vec2) * in_energy(i)
        t_mom(counter) = DOT_PRODUCT(vec1,vec2)
      END DO
    END DO
    WRITE(iw,'(1X,2A25)') 'Energies', 'Transition moments'
    WRITE(iw,'(1X,2ES25.15)') (e2h1p(ixx), 2*pi*t_mom(ixx)*t_mom(ixx), ixx=1,counter)

    CALL dealloc(buf)
    CALL dealloc(ioi)
    CALL dealloc(ioj)

    CALL dealloc(vec_start)
    CALL dealloc(vec1)
    CALL dealloc(vec2)

    CALL dealloc(in_evecs)
    !CALL dealloc(in_energy)
    CALL dealloc(fin_evecs)
    CALL dealloc(fin_energies)

!
! To be removed as soon as Dirac has its own Stieltjes routine
!
    fh_trans = fhdiag
    OPEN(fh_trans,FILE=fn_trans,STATUS='UNKNOWN',FORM='FORMATTED',    &
         ACTION='WRITE',iostat=ierror)
    WRITE(fh_trans,'(1X,2ES25.15)') (e2h1p(ixx), 2*pi*t_mom(ixx)*t_mom(ixx), ixx=1,counter)
    CLOSE(fh_trans)

    CALL dealloc(e2h1p)
    CALL dealloc(t_mom)

  END SUBROUTINE cfano_create_t_mom
END MODULE adc_fano_matmul
