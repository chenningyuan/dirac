MODULE adc_fano_exchange

REAL(8), ALLOCATABLE, DIMENSION(:,:)  :: in_evecs !alloc in lvecanaly_r, dealloc in fano
INTEGER                               :: nr_in_evecs
REAL(8), ALLOCATABLE, DIMENSION(:)    :: in_energy

REAL(8), ALLOCATABLE, DIMENSION(:)    :: fin_energies
REAL(8), ALLOCATABLE, DIMENSION(:,:)  :: fin_evecs

!
! Special arrays for complex algebra
!
COMPLEX(8), ALLOCATABLE, DIMENSION(:,:)  :: cfin_evecs

END MODULE adc_fano_exchange
