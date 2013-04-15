MODULE adc_fano

!
!  Data dictionary: common variables and arrays for the Fano ADC
!  
  INTEGER, PARAMETER                   :: fano_intbuf = 5 * 1024 * 1024 !as in ADC
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: adcbuf   !alloc in fano_write_hhp
                                                   !dealloc in fano_make_hphp
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: adcioi, adcioj !"
  INTEGER                              :: adcnbufs !number of buffers for ADC matrix
  INTEGER                              :: adcib    !book keeping number for adcbuf

  INTEGER                              :: itapfano !same as itapadc, written in a module
                                                   !file handle basis
  INTEGER                              :: fin_max=1000000  !number of final state configs
                                                   !needs to be read in later
                                               
  
  

CONTAINS

  SUBROUTINE fano_write_hhp(ckks,bufhp,ladc,crep,itape2,iw)

!
!  Purpose: Calculate the coupling block of the ADC matrix and sort
!           the matrix elements according to their description of
!           initial or final states. Write out the lower triangle
!           of the matrix for the initial state description and write
!           the coupling block of the sorted ADC matrix.
!
!           -----                -----
!     1h    |   |             1h |   |
!           -----                -----
!           |   |   ------>      |   |
!     fin   ~~~~~                |   |
!           |   |            fin ~~~~~
!           -----                -----
!
!  Author: Elke Fasshauer
!

    use adc_cfg
    use memory_allocator

    IMPLICIT NONE

!
!  Common Blocks
!
#include "../relccsd/symm.inc"
!#include "../relccsd/files.inc" !brauche ich hoffentlich nicht

!
!  Data dictionary: Calling Variables
!
    REAL(8)    :: ckks(*), bufhp(*)
    INTEGER    :: ladc              ! length of adc matrix
    INTEGER    :: crep              ! current representation
    INTEGER    :: itape2            ! file handling number of ADC matrix
    INTEGER    :: iw                ! file handling number of output file

!
!  Data dictionary: Local Variables
!
    REAL(8),ALLOCATABLE,DIMENSION(:)   :: single_col   !single column of ADC matrix
    REAL(8),ALLOCATABLE,DIMENSION(:)   :: inbuf        !buffer for initial states
    INTEGER,ALLOCATABLE,DIMENSION(:)   :: inioi, inioj !indices of initials state matrix
    INTEGER                            :: innbufs
    INTEGER                            :: inib
    REAL(8),ALLOCATABLE,DIMENSION(:,:) :: temp_array   !array to temporarily store fin part
    INTEGER                            :: temp_row     !row of temp_array
    INTEGER                            :: finrow       !row of final states
    INTEGER                            :: infh         !file handling number in matrix
    INTEGER                            :: adcrow, adccol
    INTEGER                            :: inrow, incol
    INTEGER                            :: k, ixx, fincount ! counters
    LOGICAL                            :: is_fin = .FALSE.
    INTEGER, PARAMETER                 :: jdummy = 0

!
!  Data dictionary: Variables taken from WRITE_HHP, not initialized there
!
    INTEGER                            :: off1, off2, off3
    INTEGER                            :: nok, najl

!
!  allocate arrays
!
    CALL alloc(single_col,ladc)
    CALL alloc(inbuf,fano_intbuf)
    CALL alloc(inioi,fano_intbuf)
    CALL alloc(inioj,fano_intbuf)
    CALL alloc(adcbuf,fano_intbuf)
    CALL alloc(adcioi,fano_intbuf)
    CALL alloc(adcioj,fano_intbuf)
    CALL alloc(temp_array,fin_max,3)

! ---------------------------------
!  Code
! ---------------------------------

!  clear output buffers
    single_col = 0
    inbuf      = 0
    inioi      = 0
    inioj      = 0
    adcbuf     = 0
    adcioi     = 0
    adcioj     = 0
    temp_array = 0

    inrow      = 0
    finrow     = 0
    fincount   = 0

    infh       = itapfano + 6

    WRITE(iw,*) 'Elke: Here we go!'

!                                                                       
!  write out h/h and h/2h1p buffers                                     
!
    nok  = no(crep)
    najl = nvoot(crep)

    off1 = 1
    off3 = 1

    one_column:DO k = 1, nok
      adccol = k
      off2   = 1
      CALL XCOPY(nok,ckks(off1),1,single_col(off2),1)
      off2   = off2 + nok
      CALL XCOPY(najl,bufhp(off3),1,single_col(off2),1)

      DO adcrow = adccol + 1, ladc
        incol = adccol
        !CALL subroutine
        is_fin = .TRUE.
        is_final:IF (is_fin) THEN
          finrow = finrow + 1
          IF (single_col(adcrow) /= 0.0D0) THEN
            fincount = fincount + 1
            temp_array(fincount,1) = single_col(adcrow)
            temp_array(fincount,2) = finrow
            temp_array(fincount,3) = adccol
          END IF
        ELSE is_final
          inrow = inrow + 1
          IF (single_col(adcrow) /= 0.0D0) THEN
            ! Handle the matrix of initial states
            inib   = inib + 1
            inbuf(inib) = single_col(adcrow)
            inioi(inib) = inrow
            inioj(inib) = adccol
            IF (inib == fano_intbuf) THEN
              innbufs = innbufs + 1
              WRITE(infh) (-inbuf(ixx), ixx=1, fano_intbuf),     &
                          (inioi(ixx), ixx=1, fano_intbuf),      &
                          (inioj(ixx), ixx=1, fano_intbuf),      &
                          fano_intbuf, jdummy
            END IF
            ! We write out the sorted ADC matrix
            adcib  = adcib + 1
            adcbuf(adcib) = single_col(adcrow)
            adcioi(adcib) = inrow
            adcioj(adcib) = adccol
            IF (adcib == fano_intbuf) THEN
              adcnbufs = adcnbufs + 1
              WRITE(itape2) (-adcbuf(ixx), ixx=1, fano_intbuf),  &
                            (adcioi(ixx), ixx=1, fano_intbuf),   &
                            (adcioj(ixx), ixx=1, fano_intbuf),   &
                            fano_intbuf, jdummy
              adcib = 0
            END IF
          END IF
        END IF is_final
      END DO
!
!  One coulmn has been processed, now the final states are added to
!  the end of the corresponding column
!
      DO temp_row = 1, fincount
        adcib = adcib + 1
        adcbuf(adcib) = temp_array(temp_row,1)
        adcioi(adcib) = temp_array(temp_row,2) + inrow
        adcioj(adcib) = temp_array(temp_row,3)
        IF (adcib == fano_intbuf) THEN
          WRITE(itape2) (-adcbuf(ixx), ixx=1, fano_intbuf),  &
                        (adcioi(ixx), ixx=1, fano_intbuf),   &
                        (adcioj(ixx), ixx=1, fano_intbuf),   &
                        fano_intbuf, jdummy
          adcib = 0
        END IF
      END DO 

      inrow = 0
      finrow = 0
      fincount = 0

      off1 = off1 + nok
      off3 = off3 + najl

    END DO one_column

!
!  deallocate local arrays
!
    CALL dealloc(single_col)
    CALL dealloc(inbuf)
    CALL dealloc(inioi)
    CALL dealloc(inioj)
    CALL dealloc(temp_array)

  END SUBROUTINE fano_write_hhp

END MODULE adc_fano
