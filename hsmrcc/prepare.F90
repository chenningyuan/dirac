subroutine HSCCINP
!
!
!
#include "priunit.h"

  write(LUPRI,*) 'in HSCCINP'
!
!  $HSMR 
!    NCOR=1   ! 1 core MO
!    NACT=2 
!    NFRZV=8
! $END

end subroutine HSCCINP


subroutine HSMRCC_Dirac_input


end subroutine HSMRCC_Dirac_input

!
!
!
  subroutine print_HSMRCC_Dirac_info
#include "priunit.h"

  write(LUPRI,9001)

9001 format('                                                                 ',/, &
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
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' )


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
return
end subroutine check_pointers

subroutine HSCC_STUFF(WORK,LWORK)

use Memory
#include "implicit.h"
#include "priunit.h"
!integer,intent(in) :: LWORK
real*8, intent(in) :: WORK(LWORK)
logical :: check_point = .false.
!logical :: check_point = .true.
integer (kind=4) :: ipos1,ipos2,ipos3,ipos4

write(LUPRI,*)
write(LUPRI,*) 'in HSCC_STUFF !! LWORK=',LWORK

  call print_HSMRCC_Dirac_info

  call StickMemory(WORK,LWORK)

if (check_point) then
 N=4;NNZ=2; ISIZE=N*N*NNZ
 call mymemget('INTE',ipos1,ISIZE)
 call mymemget('REAL',ipos2,ISIZE)
 call mymemget('LOGI',ipos3,ISIZE)
 call mymemget('INTE',ipos4,ISIZE)
 call assign_pointers77( PW(ipos1), PW(ipos2), PW(ipos3), PW(ipos4),N,N,NNZ)
 call check_pointers
endif

!call read_moltra_integrals
call read_mrconee; call print_mrconee_data

call read_mdcint ! read MDCINT file

call UnstickMemory

end subroutine HSCC_STUFF

! ----------------------------------
module io_moltra

! ... logical numbers ... 
integer (kind=2) :: MDINT=20, MRCONEE=21

CHARACTER*10 DATEX,TIMEX*8
real*8 :: ECORE

!integer (kind=4) :: NSP
integer  :: NSP, NFSYM, NZ, NORBT, NSYMRP, NELEC(16), NREP 
CHARACTER*14 REPNT(16)

integer, parameter :: MXREP=32,MAXSP=1000
character*4 REPN(16), REPNA(MXREP*2)

integer :: IREPSP(MAXSP*4),MULTB(MXREP*2,MXREP*2,2)
real*8 EPS1(MAXSP)

LOGICAL BREIT,SPFR

end module io_moltra

subroutine read_moltra_integrals
 use io_moltra
#include "priunit.h"
!logical is_mdcint, is_mrconee
open(MDINT,file="MDCINT",form='UNFORMATTED',status='old')
write(LUPRI,*) 'MDCINT open...'
! READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,  & 
!                                       (KR(I),KR(-I),I=1,NKR)
close(MDINT)
end subroutine read_moltra_integrals

subroutine read_mrconee

 use io_moltra

#include "implicit.h"
#include "priunit.h"

 logical :: is_mrc

     inquire(FILE='MRCONEE',EXIST=is_mrc)

     if (.not.is_mrc) call QUIT('mrconee file not present !')

       OPEN (MRCONEE,FILE='MRCONEE',FORM='UNFORMATTED',status = 'old')
       READ (MRCONEE,END=10,ERR=10) NSP,BREIT,ECORE, &
                                    NFSYM,NZ,SPFR,NORBT
       READ (MRCONEE,END=10,ERR=10) NSYMRP,(REPNT(IRP),IRP=1,NSYMRP)  & 
                      ,(NELEC(IRP),IRP=1,NSYMRP)
       GOTO 11

!C---- Keep OLD (MOLFDIR) format here for compatibility
   10     WRITE (LUPRI,*) "Missing information on MRCONEE, I assume"//  & 
          " this is a MOLFDIR file and try to continue.."
!         INTERFACE='MOLFDIR'
!         REWIND(MRCONEE)
!         READ (MRCONEE) NSP,BREIT,ECORE
!         READ (MRCONEE) NSYMRP,(REPNT(IRP),IRP=1,NSYMRP) 
!         DO IRP = 1, NSYMRP
!            NELEC(IRP) = 0
!         ENDDO

   11  CONTINUE
       READ (MRCONEE) NREP,(REPNA(IRP),IRP=1,2*NREP)
       READ (MRCONEE) ((MULTB(I,J,1),I=1,2*NREP),J=1,2*NREP)
!C      IREPSP is to be used as 2-dimensional array, read it like that
       READ (MRCONEE) (IREPSP(ISP),IREPSP(NSP+ISP),   & 
                       EPS1(ISP),ISP=1,NSP)
       CLOSE (MRCONEE,STATUS='KEEP')

end subroutine read_mrconee

subroutine print_mrconee_data
 use io_moltra
#include "implicit.h"
#include "priunit.h"

write(LUPRI,*) 'NREP=',NREP,' NSP=',NSP
write(LUPRI,*) 'REPNA: ',(REPNA(IRP),IRP=1,2*NREP)
write (LUPRI,*) 'EPS1: ',(EPS1(ISP),ISP=1,NSP)

end subroutine print_mrconee_data

!
!  ... read MDCINT 
!
subroutine  read_mdcint
   use io_moltra
#include "implicit.h"
#include "priunit.h"

logical :: is_mdc

inquire(file="MDCINT",EXIST=is_mdc); if (.not.is_mdc) call QUIT('No MDCINT file!')

open(MDINT,file="MDCINT",form='UNFORMATTED',status='old')
write(LUPRI,*) 'MDCINT open...'
!READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,( KR(I),KR(-I), I=1,NKR )

CLOSE (MDINT,STATUS='KEEP')

 10000 CALL QUIT('ERROR READING HEADER MDCINT')


end subroutine  read_mdcint


SUBROUTINE ASSIGN_POINTERS77(IW,RW,LW,IVW,N1,N2,NNZ)
!
!
! in Fortran77 style !!!
!
!
 use pointers
#include "implicit.h"
#include "priunit.h"
INTEGER, TARGET:: IW(N1,N2,NNZ),IVW(N1)
REAL*8,TARGET :: RW(N1,N2,NNZ)
LOGICAL, TARGET :: LW(N1,N2,NNZ)
! .. fill arrays !!!
 DO I=1,N1
      IVW(I) = I-10
   DO J=1,N2
   DO IZ=1,NNZ
      IW(I,J,IZ) = I*2+J+2*IZ-100
      RW(I,J,IZ) = DFLOAT(I+J+2*IZ)+0.250
      IF (I.NE.J) THEN
        LW(I,J,IZ) = .TRUE.
       ELSE
        LW(I,J,IZ) = .FALSE.
     ENDIF
   ENDDO
   ENDDO
 ENDDO

!
!  ...finally  assign pointers here
!
  IA => IW
  IVA => IVW
  RA => RW
  LA => LW

RETURN
END

