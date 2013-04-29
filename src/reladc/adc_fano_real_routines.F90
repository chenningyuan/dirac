MODULE adc_fano_real_routines

!
!  Data dictionary: common variables and arrays for the Fano ADC
!  
  INTEGER, PARAMETER                   :: fano_intbuf = 2 * 1024 * 1024 !as in ADC
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: adcbuf   !alloc in fano_write_hhp
                                                   !dealloc in fano_make_hphp
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: inbuf
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: adcioi, adcioj !"
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: inioi, inioj !"
  INTEGER                              :: adcnbufs !number of buffers for ADC matrix
  INTEGER                              :: innbufs
  INTEGER                              :: adcib    !book keeping number for adcbuf

  INTEGER                              :: itapfano !same as itapadc, written in a module
                                                   !file handle basis
  INTEGER                              :: adcfh    ! file handle for ADC matrix
  INTEGER                              :: infh
  INTEGER                              :: fin_max  !number of final state configs
                                               
  
  

CONTAINS

  SUBROUTINE fano_write_hhp(ckks,bufhp,ladc,crep,bigicra,itapadc,iw)

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
    use adc_fano_routines
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
    INTEGER, INTENT(IN)    :: crep              ! current representation
    REAL(8), INTENT(IN), DIMENSION(no(crep)*no(crep))       :: ckks
    REAL(8), INTENT(IN), DIMENSION(nvoot(crep)*nvoot(crep)) :: bufhp
    INTEGER, INTENT(IN)    :: ladc              ! length of adc matrix
    INTEGER, DIMENSION(9,nvoot(crep)), INTENT(IN) :: bigicra
    INTEGER, INTENT(IN)    :: itapadc           ! file handling number of ADC matrix
    INTEGER, INTENT(IN)    :: iw                ! file handling number of output file

!
!  Data dictionary: Local Variables
!
    REAL(8),ALLOCATABLE,DIMENSION(:)   :: single_col   !single column of ADC matrix
    INTEGER                            :: inib
    REAL(8),ALLOCATABLE,DIMENSION(:,:) :: temp_array   !array to temporarily store fin part
    INTEGER                            :: temp_row     !row of temp_array
    INTEGER                            :: finrow       !row of final states
    INTEGER                            :: adcrow, adccol
    INTEGER                            :: inrow, incol
    INTEGER                            :: k, ixx, fincount ! counters
    LOGICAL                            :: is_in  = .TRUE.
    LOGICAL                            :: is_coupl = .FALSE.
    LOGICAL                            :: is_fin = .FALSE.
    INTEGER, PARAMETER                 :: jdummy = 0

!
!  Data dictionary: Variables taken from WRITE_HHP, not initialized there
!
    INTEGER                            :: off1, off2, off3
    INTEGER                            :: nok, najl
    INTEGER                            :: dummy

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

    itapfano   = itapadc
    adcfh      = itapfano + 6 !Files need to be opened
    infh       = itapfano + 7

    WRITE(iw,*) 'Elke: Here we go!'

!                                                                       
!  write out h/h and h/2h1p buffers                                     
!
    nok  = no(crep)
    najl = nvoot(crep)

    inib  = 0
    adcib = 0

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
        IF (adcrow-nok > 0) THEN
          !WRITE(iw,*) 'adcrow = ', adcrow
          !WRITE(iw,*) 'nok    = ', nok
          !WRITE(iw,*) 'adcrow-nok = ', adcrow - nok
          CALL determine_block(adcrow-nok,adccol,bigicra,crep,is_in,&
                               is_coupl,is_fin,1)
        !  dummy = det_nr_configs(1,nvoot(crep),bigicra,crep,1,iw)
        ELSE
          is_in = .TRUE.
        END IF
        is_initial:IF (is_in) THEN
          inrow = inrow + 1
          IF (single_col(adcrow) /= 0.0D0) THEN
            ! Handle the matrix of initial states
            inib   = inib + 1
            inbuf(inib) = single_col(adcrow)
            inioi(inib) = inrow
            inioj(inib) = adccol
            !WRITE(iw,*) 'adcrow     : ', adcrow
            !WRITE(iw,*) 'in number  : ', inbuf(inib)
            IF (inib == fano_intbuf) THEN
              CALL rwrite_to_disk(inbuf,inioi,inioj,fano_intbuf, &
                                  innbufs,inib,infh)
            END IF
            !Todo: write out rest of in matrix after all loops
            ! We write out the sorted ADC matrix
            adcib  = adcib + 1
            adcbuf(adcib) = single_col(adcrow)
            adcioi(adcib) = inrow
            adcioj(adcib) = adccol
            IF (adcib == fano_intbuf) THEN
              CALL rwrite_to_disk(adcbuf,adcioi,adcioj,fano_intbuf,  &
                                  adcnbufs,adcib,adcfh)
            END IF
          END IF
        ELSE is_initial
          finrow = finrow + 1
          IF (single_col(adcrow) /= 0.0D0) THEN
            fincount = fincount + 1
            temp_array(fincount,1) = single_col(adcrow)
            temp_array(fincount,2) = finrow
            temp_array(fincount,3) = adccol
            !WRITE(iw,*) 'adcrow     : ', adcrow
            !WRITE(iw,*) 'fin number : ', temp_array(fincount,1)
          END IF
        END IF is_initial
      END DO
!
!  One column has been processed, now the final states are added to
!  the end of the corresponding column
!
      DO temp_row = 1, fincount
        adcib = adcib + 1
        adcbuf(adcib) = temp_array(temp_row,1)
        adcioi(adcib) = temp_array(temp_row,2) + inrow
        adcioj(adcib) = temp_array(temp_row,3)
        IF (adcib == fano_intbuf) THEN
          CALL rwrite_to_disk(adcbuf,adcioi,adcioj,fano_intbuf,adcnbufs,  &
                adcib,adcfh)
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
    CALL dealloc(temp_array)

  END SUBROUTINE fano_write_hhp




!
!---------------------------------------------------------------------
!
  SUBROUTINE fano_makehphp(ckks,eajl,oooo,vovo,icra,ladc,krep,  &
                           finalpos,itape1,itape2,iw)

!
!  Purpose: Create all elements of the ADC Sat Block and write them to the
!           different Fano matrices
!

    use adc_cfg
    use adc_fano_routines
    use memory_allocator

    IMPLICIT NONE

!
!  Common Blocks
!
#include "../relccsd/symm.inc"

!
!  Data dictionary: Calling Variables
!
    INTEGER,INTENT(IN)                               :: krep
    REAL(8),INTENT(IN),DIMENSION(no(krep)*no(krep))  :: ckks
    REAL(8),INTENT(IN)                               :: eajl(*)
    REAL(8),INTENT(IN)                               :: oooo(*), vovo(*)
    INTEGER,INTENT(IN),DIMENSION(6,nvoot(krep))      :: icra
    INTEGER,INTENT(IN)                               :: ladc
    INTEGER,INTENT(IN),DIMENSION(fin_max)            :: finalpos
    INTEGER,INTENT(IN)                               :: itape1, itape2
    INTEGER,INTENT(IN)                               :: iw
    

!
!  Data dictionary: Local Variables
!
    INTEGER                                :: nok, najl
    REAL(8)                                :: temp
    REAL(8),ALLOCATABLE,DIMENSION(:)       :: tempe !temporary array for diag
    REAL(8),ALLOCATABLE,DIMENSION(:)       :: temp_fin
    INTEGER                                :: irep
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:)   :: oolt
    REAL(8),ALLOCATABLE,DIMENSION(:)       :: diag
    INTEGER                                :: fincount
    INTEGER                                :: finrow, fincol
    REAL(8),ALLOCATABLE,DIMENSION(:)       :: finbuf
    INTEGER,ALLOCATABLE,DIMENSION(:)       :: finioi, finioj
    INTEGER                                :: irow, jcol, i, j, k
    INTEGER                                :: counter
    INTEGER                                :: off1, ioff 

!
!  Code
!
    WRITE(iw,*) 'FanoADC: Entering ADC Satellite block'
    nok  = no(krep) ! = mxno
    najl = nvoot(krep)

    CALL alloc(diag,nok+najl, id='FanoADC: Array of diagonal')

    ! Fill in the diagonal of the h/h block
    off1 = 1
    DO k = 1, nok
      diag(k) = ckks(off1 + k - 1)
      off1 = off1 + nok
    END DO
    
    tempe = 0
    counter = 0
    CALL alloc(tempe,najl, id='FanoADC: temporary array of diagonal elements')
    CALL alloc(temp_fin, fin_max)
    ! add the 2h1p diag sorted elements
    diag_2h1p:DO k = 1, najl
      DO fincount = 1, fin_max
        temp = finalpos(fincount)
        IF (temp == k) THEN
          temp_fin(fincount) = eajl(k)
          CYCLE diag_2h1p
        END IF  
      END DO
      counter = counter + 1
      tempe(counter) = eajl(k)
      !WRITE(iw,*) 'counter = ', counter
      !WRITE(iw,*) 'k       = ', k
    END DO diag_2h1p

    ! add final state diagonal
    ! todo: write several diagonal files from temp_fin, tempe and total
    tempe(najl-fin_max+1:najl) = temp_fin(1:fin_max)
    diag(nok+1:ladc)           = tempe(1:najl)

    CALL dealloc(temp_fin)
    CALL dealloc(tempe)

    WRITE(iw,*) 'Diagonal is written to temporary array'

    offdiag:IF(reladc_adclevel >= 2) THEN
      CALL GETOOOO(oooo)
      CALL GETVOVO(vovo)

      ! construct oolt lookup table for KCOCC
      ALLOCATE(oolt(nok,nok,nrep))
      oolt = 0
      DO irep = 1, nrep
        ioff = 1
        DO j = 1, no(irep)
          DO i = j + 1, no(irep)
            oolt(i,j,irep) = ioff
            ioff = ioff + 1
          END DO
        END DO
      END DO

      DEALLOCATE(oolt)
    END IF offdiag

    CALL dealloc(diag)

    CALL dealloc(adcbuf)
    CALL dealloc(adcioi)
    CALL dealloc(adcioj)
    CALL dealloc(inbuf)
    CALL dealloc(inioi)
    CALL dealloc(inioj)

  END SUBROUTINE fano_makehphp


!
!---------------------------------------------------------------------
!
  SUBROUTINE rwrite_to_disk(buf,ioi,ioj,intbuf,nbufs,ib,fh)
!
!  Purpose:  to write the three arrays buf, ioi, ioj of length intbuf
!            to disk
!
    IMPLICIT NONE
!
!  Data dictionary: Calling Variables
!
    INTEGER, INTENT(IN)                     :: intbuf, fh
    INTEGER, INTENT(INOUT)                  :: nbufs, ib
    REAL(8), DIMENSION(intbuf), INTENT(IN)  :: buf
    INTEGER, DIMENSION(intbuf), INTENT(IN)  :: ioi, ioj
!
!  Data dictionary: Local Variables
!
    INTEGER :: ixx
    INTEGER :: jdummy

!
!  Code
!

    jdummy = 0

    nbufs = nbufs + 1

    WRITE(fh) (-buf(ixx), ixx = 1, intbuf),    &
              (ioi(ixx), ixx = 1, intbuf),     &
              (ioj(ixx), ixx = 1, intbuf),     &
              intbuf, jdummy
    ib = 0

  END SUBROUTINE rwrite_to_disk


!
!--------------------------------------------------------------------
!
  SUBROUTINE rwrite_rest_to_disk(buf,ioi,ioj,intbuf,nbufs,ib,fh)
!
!  Purpose:  write the rest of the matrix to disk
!
    IMPLICIT NONE
!
!  Data dictionary: Calling variables
!
    INTEGER, INTENT(IN)                     :: intbuf, fh
    INTEGER, INTENT(INOUT)                  :: nbufs, ib
    REAL(8), DIMENSION(intbuf), INTENT(IN)  :: buf
    INTEGER, DIMENSION(intbuf), INTENT(IN)  :: ioi, ioj
!
!  Data dictionary: Local Variables
!
    INTEGER :: ixx
    INTEGER :: jdummy

!
!  Code
!

    jdummy = 0

    nbufs = nbufs + 1

    WRITE(fh) (-buf(ixx), ixx = 1, intbuf),   &
              (ioi(ixx), ixx = 1, intbuf),    &
              (ioj(ixx), ixx = 1, intbuf),    &
              ib, jdummy
    ib = 0

  END SUBROUTINE rwrite_rest_to_disk
    


END MODULE adc_fano_real_routines
