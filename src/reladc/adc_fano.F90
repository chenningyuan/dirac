MODULE adc_fano

!
!  Data dictionary: common variables and arrays for the Fano ADC
!  
  INTEGER                              :: fano_intbuf = 5 * 1024 * 1024 !as in ADC
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: adcbuf   !alloc in fano_write_hhp
                                                   !dealloc in fano_make_hphp
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: adcioi, adcioj !"
  INTEGER                              :: adcnbufs !number of buffers for ADC matrix
  INTEGER                              :: adcib    !book keeping number for adcbuf
                                               
  
  

CONTAINS

  SUBROUTINE fano_write_hhp(ckks,bufhp,ladc,crep,itape2)

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
!  Author: EF
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

!
!  Data dictionary: Local Variables
!
    REAL(8),ALLOCATABLE,DIMENSION(:)   :: single_col   !single column of ADC matrix
    REAL(8),ALLOCATABLE,DIMENSION(:)   :: inbuf        !buffer for initial states
    INTEGER,ALLOCATABLE,DIMENSION(:)   :: inioi, inioj !indices of initials state matrix
    INTEGER                            :: innbufs
    INTEGER                            :: inib
    INTEGER                            :: fhin         !file handling number in matrix

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

!
!  Code
!
   single_col = 0
   inbuf      = 0
   inioi      = 0
   adcbuf     = 0
   adcioi     = 0
   adcioj     = 0

!
!  deallocate local arrays
!
    CALL dealloc(single_col)
    CALL dealloc(inbuf)
    CALL dealloc(inioi)
    CALL dealloc(inioj)

  END SUBROUTINE fano_write_hhp

END MODULE adc_fano
