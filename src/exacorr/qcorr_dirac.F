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

