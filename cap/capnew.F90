!
! New basis set universal CAP implementation in DIRAC
!
!
!  Basis set functions (exponents, contraction coefficients)
! are transfered from DIRAC
!


subroutine CAP_new_stuff
!----------------------------------------------------------------------
!
!
!----------------------------------------------------------------------
use Memory
call QENTER('CAP_new_stuff')

 call HEADER('CAP_new_stuff !',-1)

  call test_basis_fu
 

call QEXIT('CAP_new_stuff')
end subroutine CAP_new_stuff



subroutine test_basis_fu
!------------------------------------------------------------
!
!
!
! control print out of atomic basis function; exponents and contraction coefficients from non-abacus
! 
! part of the DIRAC
!
!
!------------------------------------------------------------
use memory
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
#include "chrsgn.h"
#include "aosotr.h"
#include "infpar.h"


logical SPHRA

!call QENTER('test_oof_basis_fu')

 
  call HEADER('test_oof_basis_fu',-1)

      write(lupri,*) 'NLARGE=',NLARGE,' NSMALL=',NSMALL

      write(lupri,*) 'NUCIND=',NUCIND
      write(lupri,*) 'NLRGSH=',NLRGSH
      write(lupri,*) 'NSMLSH=',NSMLSH
      write(lupri,*) 'KMAX=',KMAX
      write(lupri,*) 'NORBS=',NORBS

      write(lupri,*) 'charge pre nucleus:',(charge(i),i=1,nucind)

      write(lupri,*) 'nuclear degeneracy:',(nucdeg(i),i=1,nucind)

      write(lupri,*) 'Total number of centers:',NUCDEP
      write(lupri,*)


   ! run over shells...
     ! DO ISHELL = 1, NLRGSH
      DO ISHELA = 1, KMAX

     !  L     = NHKT(ISHELL)
     !  ICENT = NCENT(ISHELL)
     !  NDEG  = NUCDEG(ICENT)
     !  INDX=NUMCF(ISHELL)

         NHKTA = NHKT(ISHELA)
         KHKTA = KHKT(ISHELA)
         KCKTA = KCKT(ISHELA)
         ICA   = LCLASS(ISHELA)

         SPHRA = SPHR(ISHELA) 
         CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)

         NCENTA = NCENT(ISHELA)
         ICENTA = NUCNUM(NCENTA,1)
         MULA   = ISTBAO(ISHELA)
         MULTA  = MULT(MULA)
         NUCA   = NUCO(ISHELA)
         NUMCFA = NUMCF(ISHELA)

         NUCDEGA  = NUCDEG(NCENTA)       
       !     ... offset for normalized coefficients ...
         JSTA   = JSTRT(ISHELA)
       !    ... coordinates
         CORAX  = CENT(ISHELA,1,1)
         CORAY  = CENT(ISHELA,2,1)
         CORAZ  = CENT(ISHELA,3,1)

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
        write(LUPRI,'(2X,A,I3,A,I3)') 'contracted shell (NUMCFA) #/ no f of uncotracted shells (NUCA)=',NUMCFA,'/',NUCA
        write(LUPRI,'(2X,A,I3,A,I3)') 'symmetric center (NCENTA) / total nu of symm.nucleis (NUCIND)=',NCENTA,'/',NUCIND

        write(LUPRI,'(2X,A,L1)')      'logical SPHRA=',SPHRA 

        write(LUPRI,'(2X,A,I3,A,I3,A,F10.3,A,A4,A,I3)')  &
   'ICENTA=',ICENTA,' NCENTA=',NCENTA,' charge:',CHARGE(NCENTA),' name:',NAMN(NCENTA), 'deg:',NUCDEGA

        write(LUPRI,'(2X,A,3I3)')    'L,M,N=',LVALUA,MVALUA,NVALUA
        write(LUPRI,'(2X,A,3F10.4)')  'coordinates:',CORAX,CORAY,CORAZ
 
        write(LUPRI,'(2X,A)') 'Exponents and normalized contract.coefficients...'
        DO  IPRIMA = 1,NUCA
           JPRIMA = JSTA + IPRIMA
           CONTA = PRICCF(JPRIMA,NUMCFA)
           EXPA = PRIEXP(JPRIMA)
           write(LUPRI,'(2X,A,F16.6,A,F20.12))') ' exp:',EXPA,'  normalized contr.coeff:',CONTA
         ENDDO

    !   write(LUPRI,*)

      ENDDO

!call QEXIT('test_oof_basis_fu')
end subroutine test_basis_fu

