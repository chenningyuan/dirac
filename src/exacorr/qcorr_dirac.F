!
!    Collection of routines needed for exacorr that use the DIRAC common
!    block infrastructure.
!
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE Get_DIRAC_MO_coefficients (IFRP,QMO,MO_LIST,NMO)
C***********************************************************************
C
C     Fill array with subset of MO coefficients and related information
C     Reads from either DFCOEF or KRMC file
C
C     Input
C        IFRP : Fermion irrep (1 gerade, 2 ungerade)
C
C     Output
C        QMO  : MOs (coefficients, energies, etc.)
C
C
C     Written by L.Visscher May 2011
C
C
C***********************************************************************

      Use Dirac_MO

#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      Integer, Intent(in) :: MO_LIST(NMO)
      Type(MO), intent(inout) :: QMO

      Real(8), Allocatable :: CMO(:),EIG(:)
      Integer, Allocatable :: IBEIG(:)
      Integer :: NAO,IR(2),NR
      LOGICAL TOBE,TOBEK

      NAO = NFBAS(IFRP,0)
C     Allocate the arrays inside QMO
      call create (QMO,NAO,NMO,NZ)
C     Fill the index part of QMO
      do i = 1, nmo
         QMO%index(i) = mo_list(i)
      enddo

C
C     Read all coefficients from file
C
      Allocate (CMO(NCMOTQ))
      Allocate (EIG(NORBT))
      Allocate (IBEIG(NORBT))
      INQUIRE(FILE='KRMCSCF',EXIST=TOBEK)

      IF (TOBEK) THEN
         CALL OPNFIL(LUKRMC,'KRMCSCF','OLD','PAMTRA')
         JRDMO = -1
         CALL RREADMO(CMO,JRDMO,1,LUKRMC)
         IF (JRDMO .EQ. 0) THEN
            IF (IPRTRA .GE. 0) THEN
               WRITE(LUPRI,'(/A)')
     &              ' (PAMTRA)  Orbitals read from' //
     &              ' label NEWORB on file KRMCSCF'
            END IF
         END IF
         CLOSE (LUKRMC)
      ELSE
         INQUIRE(FILE='DFCOEF',EXIST=TOBE)
         IF (.NOT.TOBE) THEN
            WRITE(LUPRI,'(/A)')
     &        'PAMTRA ERROR: MO coefficient file not present'
            CALL QUIT('PAMTRA: No MO coefficient file !')
         ELSE
            CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMTRA')
            CALL REACMO_no_work(LUCOEF,'DFCOEF',CMO,EIG,
     &                  IBEIG,TOTERG,14)
            CLOSE(LUCOEF,STATUS='KEEP')
         END IF
      END IF
C
C     Get the coefficients from the full list and place them in coeff
C
      CALL SELCFS (CMO(1+ICMOQ(IFRP)),IFRP,QMO%coeff,
     &             QMO%NMO,QMO%index,
     &             0,QMO%NMO,NAO,NORB(IFRP))
C
C     Get the eigenvalues and their idenfication in terms of boson irreps (for spinfree)
C
      CALL SELEIG(EIG(1+IORB(IFRP)),IFRP,QMO%energy,
     &            QMO%index,0,QMO%NMO)
      CALL SELIBEIG(IBEIG,IORB(IFRP),IFRP,QMO%boson_irrep,
     &              QMO%index,0,QMO%NMO)
C
C     Throw away the complete set, keep only the selected sets
C
      DeAllocate (CMO)
      DeAllocate (EIG)
      DeAllocate (IBEIG)
C
      END

       integer function Number_AOShells()

!      Read information from common block to get a cleaner interface

#include "implicit.h"
#include "maxorb.h"
#include "shells.h"

       Number_AOShells = nlrgsh

       end

       integer function Number_AOs(ils)

!      Get the number of basis functions (0:total, 1:large component, 2:small component)

#include "implicit.h"
#include "dcbbas.h"

       Number_AOs = nfbas(1,ils)

       end

       subroutine Initialize_Shells (shell)

       use dirac_exatensor

#include "implicit.h"
#include "maxorb.h"
#include "shells.h"
#include "aovec.h"
#include "primit.h"

       integer ishell, iPrimStart, iPrimEnd, k
       type(basis_func_info_t) :: shell(*)

       do ishell = 1, nlrgsh
          shell(ishell)%orb_momentum = nhkt(ishell)
          shell(ishell)%atom_number  = ncent(ishell)
          shell(ishell)%n_primitives = nuco(ishell)
          nullify  (shell(ishell)%exponent)
          allocate (shell(ishell)%exponent(nuco(ishell)))
          nullify  (shell(ishell)%coefficient)
          allocate (shell(ishell)%coefficient(nuco(ishell)))
          iPrimStart = jstrt(ishell) + 1
          iPrimEnd   = jstrt(ishell) + nuco(ishell)
          k          = numcf(ishell)
         shell(ishell)%exponent = priexp(iPrimStart:iPrimEnd)
         shell(ishell)%coefficient = priccf(iPrimStart:iPrimEnd,k)
         shell(ishell)%coord = cent(ishell,1:3,1)
       end do

       return 
       end

       subroutine Write_To_Mdcint (nmo,rcw,mo_tens)

!      Utility routine to write to old format, makes it easier to test the index transform

#include "implicit.h"

       integer rcw    ! this is used to toggle between complex (so, rcw=2) and real (scalar, rcw=1) tensors
       integer iclass, nclass
       integer nmo(4)
       real(8) mo_tens(2,nmo(1)*rcw,nmo(2)*rcw,nmo(3)*rcw,nmo(4)*rcw) ! the first rcw is real/complex, the other toggles orbitals/spinors
       real(8), allocatable :: rkl(:,:)
       integer, allocatable :: indk(:),indl(:)
       integer :: nmo_start(8,4),nmo_end(8,4)

!      We use no symmetry, but may have both scalar and spin-orbit tensors as input
       irepij = 1
       iprint = 0

       if (rcw == 1) then
          nclass = 1  ! this implies spinfree
       else
          nclass = 8
       end if

!      classes are subsets of integrals, make them consistent with old format (see symfint & mdintb in moltra)
!      Class 1 : (i    j    | k    l   )
!      Class 2 : (i    j    | kbar lbar)
!      Class 3 : (i    jbar | k    lbar)
!      Class 4 : (ibar j    | k    lbar)
!      Class 9 : (i    j    | kbar l   )
!      Class 10: (i    j    | k    lbar)
!      Class 11: (ibar j    | k    l   )
!      Class 12: (i    jbar | k    l   )
!      NB: classes 5 through 8 are not used, we therefore shift the final 4 in the array below to get a contiguous index

       do i = 1, 4
          do iclass = 1, 8
             nmo_start(iclass,i) = 1
             nmo_end  (iclass,i) = nmo(i)
          end do
       end do
!      bars for index 4
       do iclass = 2, 4
          nmo_start(iclass,4) = nmo(4) + 1
          nmo_end  (iclass,4) = 2 * nmo(4)
       end do
       nmo_start(6,4) = nmo(4) + 1
       nmo_end  (6,4) = 2 * nmo(4)
!      bars for index 3
       nmo_start(2,3) = nmo(3) + 1
       nmo_start(5,3) = nmo(3) + 1
       nmo_end  (2,3) = 2 * nmo(3)
       nmo_end  (5,3) = 2 * nmo(3)
!      bars for index 2
       nmo_start(3,2) = nmo(2) + 1
       nmo_start(8,2) = nmo(2) + 1
       nmo_end  (3,2) = 2 * nmo(2)
       nmo_end  (8,2) = 2 * nmo(2)
!      bars for index 1
       nmo_start(4,1) = nmo(1) + 1
       nmo_start(7,1) = nmo(1) + 1
       nmo_end  (4,1) = 2 * nmo(1)
       nmo_end  (7,1) = 2 * nmo(1)

!      Allocate memory needed by mdintb
       nkl = nmo(1) * nmo(2)
       allocate(rkl(rcw,nkl))
       allocate(indk(nkl))
       allocate(indl(nkl))
       do iclass = 1, nclass
          if (iclass.le.4) then
             jclass = iclass
          else
             jclass = iclass + 4  ! shift iclass to obey the old convention given above
          end if
          do j = nmo_start(iclass,4), nmo_end(iclass,4)
             do i = nmo_start(iclass,3),nmo_end(iclass,3)
                call mdintb(iprint,rcw,jclass,irepij,nmo(1),nmo(2),i,j,
     &               mo_tens(1,nmo_start(iclass,1):nmo_end(iclass,1),
     &                         nmo_start(iclass,2):nmo_end(iclass,2),
     &                         i,j),nkl,indk,indl,rkl)
             end do
          end do
       end do

       deallocate(rkl)
       deallocate(indk)
       deallocate(indl)

       return
       end
