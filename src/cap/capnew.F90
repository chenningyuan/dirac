!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

!
! New basis set universal CAP implementation in DIRAC
!
!  Basis set functions (exponents, contraction coefficients)
! are transfered from DIRAC
!


!=======================================================================
!
! file contains pure Fortran90 classes related to CAP project
!
!=======================================================================

  module cap

    use Vectors

! private ! must be public that subroutines can be used in F77 parts !
!#include "dcbcap.h"

!  logical, private :: cap_ovlp
!   integer(kind=4), private :: iprint
!   integer(kind=4), private :: iprcap
!  type cap_int_param
!  end type cap_int_param

  ! Number of all primitive gaussians in the basis set (s,s,px,py,pz,d...)
  ! NOTE: this value is >= than sh in dircap.F as DIRAC deals with general contractios only
   ! integer(kind=4), private :: napg = 0
    !integer(kind=4), public :: napg = 0

    integer, public :: napg = 0

   ! ... simulate R_ms(NUBF,3)
   ! type(Array2Dim) :: dir_Rms, dir_ang_pow
    type(Array2Dim), save :: dir_Rms, dir_ang_pow

! ... simulate ang_pow(NUBF,3) ...

    !type(Vect) :: dir_basis_n, dir_basis_ord, dir_ord_expon, dir_ord_coe
    type(Vect), save :: dir_basis_n, dir_basis_ord, dir_ord_expon, dir_ord_coe

  contains

  subroutine cap_init
use QMatrixes

#include "implicit.h"
#include "priunit.h"

!#include "dcbcap.h"

! iprcap = IPRCAP
! cap_ovlp = CAP_OVLP

  if (IPRCAP >= 2) then
   write(LUPRI,*) 'in init routine !'
   write(LUPRI,*) 'module iprcap =',IPRCAP
   write(LUPRI,*) 'module cap_ovlp =',CAPOVLP
!   write(LUPRI,*) 'CPRIMU ',CPRIMU
!   write(LUPRI,*) 'ALPHA: ',ALPHA
!   write(LUPRI,*) 'NUC: ',NUC

  endif

! call quit('in cap_init...')

 end subroutine cap_init

subroutine test_files
use QMatrixes
#include "implicit.h"
#include "priunit.h"

#include "mxcent.h"
#include "maxorb.h"

#include "aovec.h"
#include "maxaqn.h"
#include "primit.h"
!C#include "cbirea.h"
#include "chrxyz.h"
#include "nuclei.h"
#include "symmet.h"
#include "dummy.h"
#include "infpar.h"
!  write(LUPRI,*) 'CPRIMU ',CPRIMU
  write(LUPRI,*) 'ALPHA: ',ALPHA
!   write(LUPRI,*) 'NUC: ',NUC

  call quit('end in test files...')
end subroutine test_files


!subroutine test_basis_fu
subroutine get_napg
!------------------------------------------------------------
!
! control print out of atomic basis function;
!  exponents and contraction coefficients from non-abacus
! 
! part of the DIRAC
!
!------------------------------------------------------------
!use memory
!use Vectors
!use cap

#include "implicit.h"
#include "priunit.h"

#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "dcbgen.h"
#include "cbirea.h"
#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "frame.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "lmns.h"
#include "chrsgn.h"
#include "aosotr.h"
#include "infpar.h"

#include "dcbcap.h"


logical SPHRA

 integer :: ISHELA


call QENTER('get_napg')

 
   if (IPRCAP >= 3) then
      call HEADER('get_napg',-1)

      write(lupri,*) 'NLARGE=',NLARGE,' NSMALL=',NSMALL

      write(lupri,*) 'NUCIND=',NUCIND
      write(lupri,*) 'NLRGSH=',NLRGSH
      write(lupri,*) 'NSMLSH=',NSMLSH
      write(lupri,*) 'KMAX=',KMAX
      write(lupri,*) 'MXSHEL (must be <=KMAX) =',MXSHEL
      write(lupri,*) 'NORBS=',NORBS

      write(lupri,*) '==========================================='
      write(lupri,*) 'NPSHEL=',NPSHEL
      write(lupri,*) 'NPLSH=',NPLSH
      write(lupri,*) 'NPSSH=',NPSSH

      write(lupri,*) 'charge pre nucleus:',(charge(i),i=1,nucind)

      write(lupri,*) 'nuclear degeneracy:',(nucdeg(i),i=1,nucind)

      write(lupri,*) 'Total number of centers:',NUCDEP
      write(lupri,*)
      write(lupri,*)

  endif

      napg = 0

   ! run over shells...
     ! DO ISHELL = 1, NLRGSH
      DO ISHELA = 1, KMAX

      !  write(LUPRI,*) 'nonhurt ISHELA/KMAX=',ISHELA,'/',KMAX

     !  L     = NHKT(ISHELL)
     !  ICENT = NCENT(ISHELL)
     !  NDEG  = NUCDEG(ICENT)
     !  INDX=NUMCF(ISHELL)

         NHKTA = NHKT(ISHELA)
         KHKTA = KHKT(ISHELA)
         KCKTA = KCKT(ISHELA)
         ICA   = LCLASS(ISHELA)

      !   write(LUPRI,*) 'hurt 1 ISHELA=',ISHELA

         SPHRA = SPHR(ISHELA) 

      !  write(LUPRI,*) 'hurt 1.1 ISHELA=',ISHELA

         CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)

       ! write(LUPRI,*) 'hurt 1.1.1 ISHELA=',ISHELA
         NCENTA = NCENT(ISHELA)
      !  write(LUPRI,*) 'hurt 1.1.2 ISHELA=',ISHELA

         ICENTA = NUCNUM(NCENTA,1)

      !  write(LUPRI,*) 'hurt 1.1.3 ISHELA=',ISHELA

         MULA   = ISTBAO(ISHELA)
         MULTA  = MULT(MULA)

      !  write(LUPRI,*) 'hurt 2 ISHELA=',ISHELA

         NUCA   = NUCO(ISHELA)  ! ... uncontracted shells
         NUMCFA = NUMCF(ISHELA) ! ... contracted shells

         NUCDEGA  = NUCDEG(NCENTA)       
       !     ... offset for normalized coefficients ...
         JSTA   = JSTRT(ISHELA)
       !    ... coordinates
         CORAX  = CENT(ISHELA,1,1)
         CORAY  = CENT(ISHELA,2,1)
         CORAZ  = CENT(ISHELA,3,1)

      !  write(LUPRI,*) 'hurt 3 ISHELA=',ISHELA

   !    write(LUPRI,*) 'ISHELL=',ISHELL,' L=',L
   !    write(LUPRI,*) 'ICENT=',ICENT,' NDEG=',NDEG
   !    write(LUPRI,*) 'KHKT(ISHELL)=',KHKT(ISHELL)
   !    write(LUPRI,*) 'KCKT(ISHELL)=',KCKT(ISHELL)
   !    write(LUPRI,*) 'ICA(ISHELL)=',LCLASS(ISHELL)
   !    write(LUPRI,*) 'SPHR(ISHELL)=',SPHR(ISHELL)
   !    write(LUPRI,*) 'NCENT(ISHELL)=',NCENT(ISHELL)
   !    write(LUPRI,*) 'NUCNUM(NCENT(ISHELL),1)=',NUCNUM(NCENT(ISHELL),1)
   !    write(LUPRI,*) 'MULT(ISTBAO(ISHELL))=',MULT(ISTBAO(ISHELL))
   !    write(LUPRI,*) 'NUCO(ISHELL)=',NUCO(ISHELL)
   !    write(LUPRI,*) 'NUMCF(ISHELL)=',NUMCF(ISHELL)
   !    write(LUPRI,*) 'JSTRT(ISHELL)=',JSTRT(ISHELL)

        write(LUPRI,*)
        write(LUPRI,'(2X,A,I3,A,I3,A,I2)') 'ISHELA/KMAX=',ISHELA,'/',KMAX,' large(1)/small(2) comp., ICA:',ICA
        write(LUPRI,'(2X,A,I3,A,I3)') 'NHKTA(alias L)=',NHKTA,' cartesian degeneracy (KHKTA)=',KHKTA

    !    write(LUPRI,'(2X,A,I3,A,I3,A,I3)') 'ISHELA=',ISHELA,' NUCA=',NUCA,' NUCMFA=',NUMCFA
        write(LUPRI,'(2X,A,I3,A,I3)')  & 
          'contracted shell (NUMCFA) #/ no f of uncontracted shells (NUCA)=',NUMCFA,'/',NUCA
        write(LUPRI,'(2X,A,I3,A,I3)') & 
          'symmetric center (NCENTA) / total nu of symm.nucleis (NUCIND)=',NCENTA,'/',NUCIND

        write(LUPRI,'(2X,A,L1)')      'logical SPHRA=',SPHRA 

        write(LUPRI,'(2X,A,I3,A,I3,A,F10.3,A,A4,A,I3)')  &
   'ICENTA=',ICENTA,' NCENTA=',NCENTA,' charge:',CHARGE(NCENTA),' name:',NAMN(NCENTA), 'deg:',NUCDEGA

        write(LUPRI,'(2X,A,3F10.4)')  'coordinates:',CORAX,CORAY,CORAZ

        write(LUPRI,'(2X,A)')     '........L,M,N values:'
        DO ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
         write(LUPRI,'(2X,A,3I3)') 'LVALA,MVALA,NVALA:',LVALA,MVALA,NVALA
        ENDDO
 
        write(LUPRI,'(2X,A)') 'Exponents and normalized contract.coefficients...'
        DO  IPRIMA = 1,NUCA
           JPRIMA = JSTA + IPRIMA
           CONTA = PRICCF(JPRIMA,NUMCFA)
           EXPA = PRIEXP(JPRIMA)
           write(LUPRI,'(2X,A,F16.6,A,F20.12)')  &
           ' exp:',EXPA,'  normalized contr.coeff:',CONTA

           ! that's it !!!
           !IF (CONTA <> 0.0D0)  THEN
           IF (CONTA .ne. 0.0D0)  THEN
             napg = napg + KHKTA
           ENDIF

         ENDDO

        !napg = napg + (NUCA*KHKTA)
       ! napg = napg + NUCA
      !  write(LUPRI,*) '>>> accumulated napg:',napg

      !  write(LUPRI,*) 'hurt ISHELA/KMAX:',ISHELA,'/',KMAX

      ENDDO

      write(LUPRI,*) 'determined full_size=', napg
      CALL FLSHFO(LUPRI)

! allocate space for Rms,...

call QEXIT('get_napg')
end subroutine get_napg

end module cap

subroutine CAP_new_stuff
!----------------------------------------------------------------------
!
! Called from CAP_STUFF/prepcap90.F90
!
!
!----------------------------------------------------------------------
use cap
use Vectors
#include "priunit.h"

call QENTER('CAP_new_stuff')

call HEADER('CAP_new_stuff !',-1)

!.. prepare inputs for AOINTEGRAL routine...
call get_napg

call set_vect_iprint(10); write(LUPRI,*) 'get_vect_iprint=',get_vect_iprint()

 write(LUPRI,*) 'napg=',napg

 write(LUPRI,*) 'going to allocate dir_Rms...' 
 CALL FLSHFO(LUPRI)

 dir_Rms = Make2DimArray('REAL',napg,3,'dirac R_ms 2D REAL array structure') 
 write(LUPRI,*) 'dir_Rms allocated...' 
 CALL FLSHFO(LUPRI)
 !call Fill2D_Vect_re(dir_Rms,12.24D0)
 call Fill2D_Vect(dir_Rms,12.24D0)
! call print_2Darray_attrib(dir_Rms)
 call PrintObject(dir_Rms)
 call DeAllocObject(dir_Rms)

 dir_ang_pow = Make2DimArray('INTE',napg,3,'dirac ang_pow 2D INTEGER array structure') 
 !call Fill2D_Vect_int(dir_ang_pow,7)
 call Fill2D_Vect(dir_ang_pow,-7)
 !call print_2Darray_attrib(dir_ang_pow)
 call PrintObject(dir_ang_pow)
 call ZeroObject(dir_ang_pow);  call PrintObject(dir_ang_pow)
 call DeAllocObject(dir_ang_pow)


!call quit('CAP_new_stuff!')

call QEXIT('CAP_new_stuff')
end subroutine CAP_new_stuff


