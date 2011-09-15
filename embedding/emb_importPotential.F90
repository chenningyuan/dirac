!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

module FDE_ImportPotentials

   implicit none
   private
   
   public FDE_SaveAOPROPER_StaticEmbPot

   contains

! routines 
! 
!     the following sets up the matrix elements for ground-state embedding 
!     potential, as a so-called "numerical operator" (that is, we have
!     values of the given operator on a set of gridpoints P, and contract
!     that with the value of the orbitals (or their derivatives, whatever)
!     and the integration weight to get the final integral: 
!
!     M_ab = < phi^n_a | operator | phi^m_b > = \sum^P_i w_i phi^n_a(i) operator(i) phi^m_b(i) 
!
!     where i is the i-th grid point, w_i the integration weight, and
!     phi^{n,m}_{a,b}(i) is the value of the {n,m}-th order derivative of
!     the (contracted) basis function {a,b} at point i.
!
!     information about the ao-basis blocks active for a given operator is
!     carried in the variable COMB, as indicated up. In the case of the
!     embedding potential for FDE, we have a diagonal operator just like
!     for the molecular field.
!  

!
!
! constructs the matrix elements for the embedding potential and adds
! make_numerical_operator_matrix: calculates matrix elements for 
!                                 a numerical operator on the grid 
!
!
         subroutine FDE_SaveAOPROPER_StaticEmbPot(comb,iprint)
 
            use memory_allocator
            use numoper_integrals
!           use fde_mod
            use fde_types
            use fde_cfg
            use fde_data
            use fde_io

#include "implicit.h"
#include "priunit.h"
            LOGICAl DOINT(4)
            CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4
            character*80 string
#include "dcbbas.h"
#include "mxcent.h"
#include "nuclei.h"

            integer :: file_unit, iprint, i, irrep_vemb
            real(kind=8), allocatable :: a_oneint(:)
            character(len=60) :: file_name
            type(fde_files)   :: ftmp

            call fde_import_static
            call fde_get_files_info(ftmp)
      
            file_unit = ftmp%embpot%unit
            file_name = ftmp%embpot%name
            call fde_open_file(file_name,file_unit)

 
            label = 'FDEVEMB '
            irrep_vemb = 0  ! embedding potential is totally symmetric

            read(comb,'(4L1)',err=1000) (doint(i),i=1,4)

            call alloc(a_oneint,nnbbasx)  
            a_oneint = 0.0d0
            call NumOper_OneElOpMatrix(a_oneint,file_unit,doint,irrep_vemb)

            IF(IPRINT.GE.4) THEN
              write (string,*) 'Integrals over (numerical) operator: '//label
              call HEADER(string,-1)
              call OUTPAK(a_oneint,NTBAS(0),1,6)
            ENDIF
!
!        Generate integral labels
!
            CALL GETDAT(RTNLBL(1),RTNLBL(2))
            RTNLBL(2)(1:2) = 'SY'
            WRITE(RTNLBL(2)(3:4),'(I2)') 1
            RTNLBL(2)(5:8) = COMB
!
!        Write integrals to file
!
            call WRTPRO(a_oneint,NNBBASX,LABEL,RTNLBL,IPRINT)
!
            call fde_close_file(file_unit)
            call dealloc (a_oneint)

            return
 1000       continue
!
!     Not able to read DOINT information from COMB
!
            write(LUPRI,'(A,A)') 'FDE_SaveAOPROPER_StaticEmbPot: '// &
                                 'Not able to read COMB =',COMB
            call quit('Vemb2AOPROPER: Not able to read COMB')

         end subroutine FDE_SaveAOPROPER_StaticEmbPot 

end module FDE_ImportPotentials
