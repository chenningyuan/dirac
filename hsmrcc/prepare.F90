subroutine HSCCINP
!
! called from PSIINP
!
#include "priunit.h"

 write(LUPRI,*) 'in HSCCINP'


end subroutine HSCCINP


  subroutine HSMRCC_Dirac_input


  end subroutine HSMRCC_Dirac_input

!
!
!
  subroutine print_HSMRCC_Dirac_info
#include "priunit.h"

  write(LUPRI,9001)

9001 format(                                                                   ,/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|        I N T E R F A C E                                      |',/, &
            '|      PROGRAM FOR EFFECTIVE AND INTERMEDIATE HAMILTONIAN       |',/, &
            '|     HILBERT SPACE MULTI-REFERENCE COUPLED CLUSTER METHODS     |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                                                               |',/, &
            '|           Shamasundar. K. R. and Sourav Pal                   |',/, &
            '|           Theoretical Chemistry Group                         |',/, &
            '|           Physical Chemistry Division                         |',/, &
            '|           National Chemical Laboratory (NCL)                  |',/, &
            '|           Pashon Road, Pune - 411 008                         |',/, &
            '|           INDIA                                               |',/, &
            '|                                                               |',/, &
            '|                       AND                                     |',/, &
            '|                                                               |',/, &
            '|           Ephraim Eliav and Uzi Kaldor                        |',/, &
            '|           School of Chemistry                                 |',/, &
            '|           Tel Aviv University                                 |',/, &
            '|           Tel Aviv, 69978                                     |',/, &
            '|           ISRAEL                                              |',/, &
            '|                                                               |',/, &
            '|    Note : This version is ported for opteron machines at      |',/, &
            '|           University of Waterloo and Tel Aviv University.     |',/, &
            '|                                                               |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                INFORMATION ON CURRENT RUN                     |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
           )


  end subroutine print_HSMRCC_Dirac_info


!
! prepare file for the gamess-like input file
!
! $GMS/"syminf.f90"
  subroutine make_gamess_syminf

  end subroutine make_gamess_syminf

!--------------------------------------------------
module pointers

integer, dimension(:,:,:), pointer:: IA
integer, dimension(:), pointer:: IVA
real*8, dimension(:,:,:), pointer:: RA
logical, dimension(:,:,:), pointer:: LA

end module pointers
!--------------------------------------------------

subroutine check_pointers
! checks content of all filled arrays, which are absorbed in pointers module ....
use pointers
#include "implicit.h"
#include "priunit.h"

 do i=1,4

  write(LUPRI,'(A,i2,a,i4)') 'vector IVA(',i,')=',IVA(i)

 do j=1,4
 do iz=1,2

  write(LUPRI,'(A,i2,a,i2,a,i2,a,i4)') 'IA(',i,',',j,',',iz,')=',IA(i,j,iz)
  write(LUPRI,'(A,i2,a,i2,a,i2,a,f12.6)') 'RA(',i,',',j,',',iz,')=',RA(i,j,iz)
  write(LUPRI,'(A,i2,a,i2,a,i2,a,l1)') 'LW(',i,',',j,',',iz,')=',LA(i,j,iz)

 enddo
 enddo
 enddo

end subroutine check_pointers

subroutine HSCC_STUFF(WORK,LWORK)
use Memory
#include "implicit.h"
#include "priunit.h"
!integer,intent(in) :: LWORK
real*8, intent(in) :: WORK(LWORK)

integer (kind=4) :: ipos1,ipos2,ipos3,ipos4

write(LUPRI,*)
write(LUPRI,*) 'in HSCC_STUFF !! LWORK=',LWORK
call StickMemory(WORK,LWORK)

N=4;NNZ=2; ISIZE=N*N*NNZ
call mymemget('INTE',ipos1,ISIZE)
call mymemget('REAL',ipos2,ISIZE)
call mymemget('LOGI',ipos3,ISIZE)

call mymemget('INTE',ipos4,ISIZE)

call assign_pointers77( PW(ipos1), PW(ipos2), PW(ipos3), PW(ipos4),N,N,NNZ)

call check_pointers

call UnstickMemory

end subroutine HSCC_STUFF


