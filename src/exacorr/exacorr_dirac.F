!
!    Collection of routines needed for exacorr that use the DIRAC common
!    block infrastructure.
!
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXA_MOLTRA()
C***********************************************************************
C
C     Driver routine for exacorr MOLTRA implementation
C
C     Written by L.Visscher Jan 2017
C
C***********************************************************************

      Use exacorr_ao_to_mo

#include "implicit.h"! because common blocks may assume this
#include "../moltra/dcbtra.h"  ! to access tra4_indstr
#include "dcborb.h"  ! to access nesh (number of positive energy "electronlike" orbitals)
#include "dgroup.h"  ! to access nfsym (should be 1)

      character (len=72)   :: orbital_string(4) ! String with range of MO's to be transformed for each index
      integer              :: kvec              ! Dummy variable needed in first call to numlst
      integer, allocatable :: mo_list(:),mo_occ(:),mo_vir(:)
      integer              :: nmo(4),noccu,nvir

      if (nfsym.ne.1) then
         call quit ('exacorr can not handle symmetry')
      end if
C
C     Transfer necessary information from common blocks
      orbital_string = tra4_indstr(:,1) ! This is to be changed once we have an exacorr-specific input driver
C
C     Loop over indices to convert string information to list of orbitals
C
C     Count number of orbitals to be transformed for each index 
      ntmo = 0
      do i = 1, 4
        nmo(i) = 0
        call numlst (orbital_string(i),kvec,1,1,nesh(1),1,nmo(i))
        ntmo = ntmo + nmo(i)
      end do

C     Fill the array with orbital indices
      allocate (mo_list(ntmo))
      imo = 1
      do i = 1, 4
       call numlst
     & (orbital_string(i),mo_list(imo),nmo(i),1,nesh(1),1,nmo(i))
       imo = imo + nmo(i)
      end do
C
C     Transformed AO to MO integrals are written to MDCINT
      call mdinti (0) 
      call compute_mo_integrals (nmo, mo_list)
      call mdintf
C
      deallocate (mo_list)

      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXA_CC()
C***********************************************************************
C
C     Driver routine for exacorr CC implementation
C
C     Written by L.Visscher Jul 2017
C
C***********************************************************************

      Use exacc_cfg
      Use exacorr_cc

#include "implicit.h"! because common blocks may assume this
#include "dcborb.h"  ! to access nesh (number of positive energy "electronlike" orbitals)
#include "dgroup.h"  ! to access nfsym (should be 1)

      integer              :: kvec              ! Dummy variable needed in first call to numlst
      integer, allocatable :: mo_list(:),mo_occ(:),mo_vir(:)
      integer              :: nmo(4),noccu,nvir

      if (nfsym.ne.1) then
         call quit ('exacorr can not handle symmetry')
      end if
C
      noccu = 0
      call numlst (string_occupied,kvec,1,1,nesh(1),1,noccu)
      allocate (mo_occ(noccu))
      call numlst (string_occupied,mo_occ,noccu,1,nesh(1),1,noccu)

      nvir = 0
      call numlst (string_virtual,kvec,1,1,nesh(1),1,nvir)
      allocate (mo_vir(nvir))
      call numlst (string_virtual,mo_vir,nvir,1,nesh(1),1,nvir)

      call cc_driver (noccu,nvir,mo_occ,mo_vir)

      deallocate (mo_occ)
      deallocate (mo_vir)

      end
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

       subroutine Initialize_Shells (shell,ierr)

!      get information from common blocks about the basis set
!      raise an error for contracted sets that can not yet be handled by Interest

       use exacorr_datatypes

#include "implicit.h"
#include "maxorb.h"
#include "shells.h"
#include "aovec.h"
#include "primit.h"

       integer ishell, iPrimStart, iPrimEnd, k, ierr
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
!        In DIRAC uncontracted shells are stored as quasi-contracted, fix this.
         call compress_unc_shell (shell(ishell),ierr)
       end do

       return 
       end

       subroutine compress_unc_shell (shell,ierr)
!If basis functions were grouped together in the input, they appear as quasi-contracted functions
!with a lot of zero coefficients. Remove these zeroes and reduce the length to just one function.
!Routine will produce an error, if the function is a true contraced function (more than one non-zero coefficient)

          use exacorr_datatypes

          integer, intent(out)     :: ierr
          type(basis_func_info_t), intent(inout) :: shell

          integer i, j, non_zero
          real(8), allocatable :: exponents(:), coefficients(:)
          real(8), parameter   :: treshold=1.D-12

          allocate(exponents(shell%n_primitives))
          allocate(coefficients(shell%n_primitives))
          exponents    = shell%exponent
          coefficients = shell%coefficient

          ierr = 0
          non_zero = 0
          do i = 1, shell%n_primitives
             if (abs(coefficients(i)) > treshold) then
                non_zero = non_zero + 1
                j = i
             end if
          end do

          if (non_zero > 1) then
             ierr = 1
             return   ! leave shell untouched and return with error condition
          else if (non_zero == 0) then
             print*, " WARNING: compress_unc_shell found empty shell"
             ierr = 2
             j = 1    ! set j to the first function, to avoid complications with zero length arrays
          end if

          deallocate (shell%exponent)
          deallocate (shell%coefficient)
          allocate   (shell%exponent(1))
          allocate   (shell%coefficient(1))

          shell%n_primitives = 1
          shell%exponent     = exponents(j)
          shell%coefficient  = coefficients(j)

          deallocate (exponents)
          deallocate (coefficients)

       end subroutine compress_unc_shell

       subroutine Write_To_Mdcint (nmo,rcw,mo_tens)

!      Utility routine to write to the DIRAC moltra format, makes it easier to test the index transform
 
!      Note that the integral matrix is transposed in moltra, this is why the 3 and 4th indices are written to the left below
!      This routine should therefore be rewritten if the matrix transpose I(ij,kl) .ne. I(kl,ij)

#include "implicit.h"

       integer rcw    ! this is used to toggle between complex (so, rcw=2) and real (scalar, rcw=1) tensors
       integer iclass, nclass
       integer nmo(4)
       real(8) mo_tens(rcw,nmo(3)*rcw,nmo(4)*rcw,nmo(1)*rcw,nmo(2)*rcw) ! the first rcw is real/complex, the other toggles orbitals/spinors
       real(8), allocatable :: fmat(:,:),rkl(:,:)
       integer, allocatable :: indk(:),indl(:)
       integer :: nmo_start(8,4),nmo_end(8,4)

       iprint = 0

!      We use no symmetry, but may have both scalar and spin-orbit tensors as input
       irepij = 1
       if (rcw == 1) then
          nclass = 1  ! this implies spinfree
       else
          nclass = 8
       end if

!      We only need to write a subset of integrals, the rest will be generated by permutation and/or Kramers' symmetry
 
!      Classes are defined to have a standard subset of integrals, make these consistent with the definition in mdintb (in moltra)
!      Class 1 : (k    l    | i    j   )
!      Class 2 : (kbar lbar | i    j   )
!      Class 3 : (k    lbar | i    jbar)
!      Class 4 : (kbar l    | i    jbar)
!      Class 9 : (kbar l    | i    j   )
!      Class 10: (k    lbar | i    j   )
!      Class 11: (k    l    | ibar j   )
!      Class 12: (k    l    | i    jbar)

!      NB: classes 5 through 8 are not used, we therefore shift the final 4 classes in the array below to get a contiguous index

       do i = 1, 4
          do iclass = 1, 8
             nmo_start(iclass,i) = 1
             nmo_end  (iclass,i) = nmo(i)
          end do
       end do

!      define bars for class 2
       nmo_start(2,3) = nmo(3) + 1
       nmo_end  (2,3) = 2 * nmo(3)
       nmo_start(2,4) = nmo(4) + 1
       nmo_end  (2,4) = 2 * nmo(4)
!      define bars for class 3
       nmo_start(3,2) = nmo(2) + 1
       nmo_end  (3,2) = 2 * nmo(2)
       nmo_start(3,4) = nmo(4) + 1
       nmo_end  (3,4) = 2 * nmo(4)
!      define bars for class 4
       nmo_start(4,2) = nmo(2) + 1
       nmo_end  (4,2) = 2 * nmo(2)
       nmo_start(4,3) = nmo(3) + 1
       nmo_end  (4,3) = 2 * nmo(3)
!      define bars for class 5=9
       nmo_start(5,3) = nmo(3) + 1
       nmo_end  (5,3) = 2 * nmo(3)
!      define bars for class 6=10
       nmo_start(6,4) = nmo(4) + 1
       nmo_end  (6,4) = 2 * nmo(4)
!      define bars for class 7=11
       nmo_start(7,1) = nmo(1) + 1
       nmo_end  (7,1) = 2 * nmo(1)
!      define bars for class 8=12
       nmo_start(8,2) = nmo(2) + 1
       nmo_end  (8,2) = 2 * nmo(2)

!      Allocate memory needed by mdintb
       nkl = nmo(3) * nmo(4)
       allocate(rkl(rcw,nkl))
       allocate(indk(nkl))
       allocate(indl(nkl)) 
       allocate(fmat(nkl,rcw))

!      loop over classes to write out subsets of integrals
       do iclass = 1, nclass
          if (iclass.le.4) then
             jclass = iclass
          else
             jclass = iclass + 4  ! shift iclass to obey the convention used in mdintb (given above)
          end if

!         loop over blocks of integrals, taking i and j as outermost indices like is done in moltra
          do j = nmo_start(iclass,2), nmo_end(iclass,2)
             if (j.le.nmo(2)) then ! mdintb takes as index Kramers pairs not spinors
                jj = j
             else
                jj = j - nmo(2)
             endif
             do i = nmo_start(iclass,1),nmo_end(iclass,1)
               if (i.le.nmo(1)) then ! mdintb takes as index Kramers pairs not spinors
                  ii = i
               else
                  ii = i - nmo(1)
               endif

!              In moltra the real/imaginary index within a block is last rather than first, so we need to reorganize

!              If the particle/particle (index 12 <-> index 34) symmetry is not fulfilled we should reorganize further
!              this can be done here, but will be quite inefficient and not worth the trouble for this utility routine that
!              is to discarded later anyhow.
     
               rkl = reshape(mo_tens(1:rcw,
     &                       nmo_start(iclass,3):nmo_end(iclass,3),
     &                       nmo_start(iclass,4):nmo_end(iclass,4),i,j),
     &                       (/rcw, nmo(3)*nmo(4)/))
               fmat = transpose(rkl)

!              write a block of integrals to file. mdintb takes care of thresholding and labelling
               call mdintb(iprint,rcw,jclass,irepij,nmo(3),nmo(4),ii,jj,
     &                     fmat,nkl,indk,indl,rkl)

             end do
          end do
       end do

       deallocate(fmat)
       deallocate(rkl)
       deallocate(indk)
       deallocate(indl)

       return
       end
