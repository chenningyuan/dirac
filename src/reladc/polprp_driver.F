!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!      *****************************************************
!      ****     POLARIZATION  PROPAGATOR  MODULE        ****
!      ****                   for                       ****
!      ****           ELECTRONIC EXCITATIONS            ****
!      *****************************************************
!
!                       AUTHOR: M.Pernpointner
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PAMPOLPRP ()
!
!---------------Description-----------------------------------------------
!
!    Double group symmetry adapted and parallel
!    polarization propagator code for the calculation of electronic
!    excitations.
!    In parallel execution this is the part called by the dirac(...)
!    subroutine, the MASTER PART. 
!    It has to wake up the slaves with the correct ITASK
!    number in order to run the shared calculation.
!    Afterwards it gives itself to the shared calculation.
!
!---------------Common Blocks--------------------------------------
!
      use interface_to_mpi
#include "implicit.h"
#include "maxorb.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbpsi.h"
#include "dcbgen.h"
#include "dcbmp2no.h"
#include "infpar.h"
#include "../relccsd/ccpar.inc"
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      REAL*8 SEC,SEC0
      LOGICAL CA
      CHARACTER*10 INTFCE
#ifdef VAR_MPI
      dimension istat(df_mpi_status_size)
#endif
!
!---------------Executable code--------------------------------------
!
      CALL CPUUSED(SEC0)
      OPEN(LUCMD,FILE = 'DIRAC.INP')
      IN = LUCMD
      IW = LUPRI
!
! iface set to DIRAC6 as default, in line with having scheme 6 as default in moltra
      INTFCE = 'DIRAC6    '
      IF( IAOMOSCHEME .eq. 4 ) INTFCE = 'DIRAC     '
      print*,' Current Dirac interface: ',INTFCE
!
!  determine condition when we have complex arithmetic !
!
      CA = NZ.GE.2.AND..NOT.(SPINFR.OR.LEVYLE)
!MI   ... if CAP is on, switch to complex algebra
      IF (CAP) CA = .TRUE.
#if defined (VAR_MPI)
      MASTER = MPARID
      MYPROC = MYTID
!     Add the master node
      NMPROC = NUMNOD + 1
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif
!
!  Do a "wake up" with the correct ITASK in order to guide the slaves
!  to the corresponding POLPRP section.
!
!
!     ITASK = 12 for parallel POLPRP
!
      print*, '+++++ in PAMPOLPRP, parcal: ',PARCAL
      print*, '+++++ trying to wake up slaves.'
      IF (PARCAL) CALL DIRAC_PARCTL( 12 )

#if defined (VAR_MPI)
      IF (NMPROC .GT. 1) THEN
         call interface_mpi_BCAST_l0(CA,1,MASTER,global_communicator)
         call interface_mpi_BCAST_l0(MP2_DENRUN,1,MASTER,
     &                  global_communicator)
      END IF
!
!     if MP2-NO is enabled, tell the slaves about the correct interface
!
      IF( MP2_DENRUN )THEN
        call interface_mpi_BCAST(IAOMOSCHEME,1,MASTER,
     &                 global_communicator)
        IF( IAOMOSCHEME .eq. 4 ) INTFCE = 'DIRAC     '
        call interface_mpi_BCAST(INFO_STRING(1),62,MASTER,
     &                 global_communicator)
        call interface_mpi_BCAST(INFO_STRING(2),62,MASTER,
     &                 global_communicator)
        call interface_mpi_BCAST(INFO_STRING(3),62,MASTER,
     &                 global_communicator)
      END IF
#endif
!
!     Enter the generic POLPRP_MAIN routine (master becomes now primus inter pares)
!
      CALL POLPRP_MAIN(IN,IW,CA,INTFCE)
!
!     Return the slaves to the general menu routine.
!
!
      if(numnod.eq.0) return
!
!     Find the slaves and release them
!
      NTEST = -1
#if defined (VAR_MPI)
      DO I = 1,NUMNOD
         CALL interface_MPI_IRECV(NODE,1,df_MPI_ANY_SOURCE,20,
     &                            global_communicator,IREQ)
         CALL interface_MPI_WAIT(IREQ,ISTAT)
         CALL interface_MPI_SEND(NTEST,1,NODE,30,
     &                           global_communicator)
      ENDDO
#endif
!
!     Let the slaves sleep (if relevant for this operating system)
!
      CALL DIRAC_PARCTL( 0 )

      CALL CPUUSED(SEC)
      WRITE (IW,1000) SEC-SEC0,SEC0,SEC
      WRITE (IW,'(//80A1)') ('#',I=1,80)
      CLOSE(LUCMD)
!
 1000 FORMAT (//' CPU time (seconds) used in POLPRP:',T50,F14.4/
     & ' CPU time (seconds) used before POLPRP:',T50,F14.4/
     & ' CPU time (seconds) used in total sofar:',T50,F14.4//
     & '  --- Normal end of POLPRP Run ---')
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE POLPRP_NODE
!
!  similar to CC_NOD due to identical requirements
!
      use interface_to_mpi

#include "implicit.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dcbmp2no.h"
#include "infpar.h"
#if defined (VAR_MPI)
      dimension istat(df_mpi_status_size)
#endif
#include "priunit.h"
#include "dgroup.h"
#include "../relccsd/ccpar.inc"
#include "../relccsd/files.inc"
      CHARACTER*10 INTFCE, MP2NOF
      LOGICAL CA,CCEND
C
C     arrange for the MPI stuff and correct node number
C     to the total number of running invocations.
C
#if defined (VAR_MPI)
      MASTER = MPARID
      MYPROC = MYTID
C     Add the master node
      NMPROC = NUMNOD + 1
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif

C aspg 20080903
C iface set to DIRAC6 as default, in line with having scheme 6 as default in moltra
      INTFCE = 'DIRAC6    '
      IF( IAOMOSCHEME .eq. 4 ) INTFCE = 'DIRAC     '
!     MP2_DENRUN = .FALSE.
      print*,'+++++ current (node) interface: ',INTFCE
C
C     Get the correct arithmetic (NZ is not initialized on slaves if you do not run HF !)
C
#if defined (VAR_MPI)
      IF (NMPROC .GT. 1) THEN
         print*,'+++++ node',MYPROC,' entering mpi_bcast1.'
         call interface_mpi_BCAST_l0(CA,1,MASTER,global_communicator)
         print*,'+++++ node',MYPROC,' entering mpi_bcast2.'
         call interface_mpi_BCAST_l0(MP2_DENRUN,1,MASTER,
     &                  global_communicator)
      END IF
      IF( MP2_DENRUN )THEN
        call interface_mpi_BCAST(IAOMOSCHEME,1,MASTER,
     &                  global_communicator)
        IF( IAOMOSCHEME .eq. 4 ) INTFCE = 'DIRAC     '
      END IF
#endif
C
C     create a node-unique filename as output file. Important on
C     shared file systems. Otherwise all the output gets mingled in one
C     file
C
      CCBASF="POL_PRP.OUT"
      
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (CCFILN,'(A11,A1,I1)') CCBASF,'.',MYPROC
         LCCFIL    = 13
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (CCFILN,'(A11,A1,I2)') CCBASF,'.',MYPROC
         LCCFIL    = 14
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (CCFILN,'(A11,A1,I3)') CCBASF,'.',MYPROC
         LCCFIL    = 15
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND POLPRP_NOD Routine!")
      ENDIF
#else
      CCFILN    = CCBASF
      LCCFIL    = 11
#endif
C
C
C     open the local input file and the node specific output file.
C     The associated file handle MSLVOUT is transferred to CCMAIN.
C     Every access to the local stdout handle then automatically writes
C     to the corresponding output file.
C
      OPEN(LUCMD,FILE = 'DIRAC.INP')
      OPEN(MSLVOUT,FILE = CCFILN(1:LCCFIL))
      CALL POLPRP_MAIN(LUCMD,MSLVOUT,CA,INTFCE)
      CLOSE (MSLVOUT,STATUS='KEEP')
      CLOSE (LUCMD,STATUS='KEEP')
      write(*,*) 'Node',MYPROC,' has arrived after POLPRP_MAIN!'
C
C     The slave is finished, wait for the master
C
      CCEND = .FALSE.
#if defined (VAR_MPI)
      CALL interface_MPI_ISEND(MYTID,1,MPARID,20,
     &                         global_communicator,IREQ)
C
    1    CONTINUE
         CALL SLEEP(1)
         call interface_mpi_iprobe(MPARID,30,CCEND,global_communicator,
     &                             istat)
         IF (.NOT.CCEND) GOTO 1
C
      CALL interface_MPI_RECV(NTEST,1,MPARID,30,global_communicator)
#endif
C
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE POLPRP_MAIN (IN123,IW123,CA,INTFCE)

      use memory_allocator
      use interface_to_mpi

      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Driver for Polarization Propagator 
!
!---------------Calling variables--------------------------------------
!
      INTEGER IN123,IW123
      LOGICAL CA
      CHARACTER*10 INTFCE
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/files.inc"
#include "../relccsd/inpt.inc"
#include "../relccsd/symm.inc"
#include "../relccsd/eqns.inc"
#include "../relccsd/complex.inc"
#include "../relccsd/ccpar.inc"
#include "../relccsd/waio.h"
#include "dcbham.h"
!
!---------------Local variables--------------------------------------
!
      integer , parameter        :: MAXSP=1000
      INTEGER IREPSP(MAXSP*4),IREPSPI(MAXSP*MXREP*2)
      REAL*8 EPS(MAXSP),EPS1(MAXSP)
      REAL*8 CPU0,CPU1,CPU2,CPU3

      real(8), allocatable :: CC(:)

      integer :: ncore_elec_adc = 0  ! new argument for RDSYMI, needed in RELCCSD, not in ADC.

      integer           :: maxcor,mxcorr,mxcorr1,mword_limit
      real*8            :: total_mem_new = 0   !number of bytes, if number of megawords is given in --aw
      real*8            :: new_mword_limit  = 0

      integer     :: MAXCORE    ! contains --mw * 1000 * 1000 (number of words, decimal mega)
      integer     :: maxcore_info(2)
      character*6 :: maxcore_infc(2)
!
!---------------interface region --------------------------------------
!

      interface

        SUBROUTINE POLPRP_PAREXE(i1,i2,ra1,ia1,i3)
          integer                               :: i1,i2,i3
          real*8, dimension(:)                  :: ra1
          integer,dimension(:)                  :: ia1
        END SUBROUTINE

      end interface
!
!---------------executable code --------------------------------------
!
      maxcore_info = 0
      maxcore_infc = (/'CCSETI','POLPRP'/)

      call legacy_lwork_get(MAXCORE)

!
!---------------Executable code--------------------------------------

!
!     The INTERFACE/EQNS variable is defined in ccinpt.inc and is needed
!     by CCSETI. The default scheme is DIRAC6. This is also hardcoded
!     in the RELCCSD code.
!
      CARITH = CA
      IN=IN123
      IW=IW123
      INTERFACE = INTFCE
      EQNS='RELADC'
!
!    indicate next module run
!
      write(IW,1999)
 1999 format(//4X,90('*')/
     &   4X,21('*'),5X,'Start of Polarization Propagator Module',4X,
     &      21('*')/
     &   4X,90('*')//)
      WRITE(IW,*)
      WRITE(IW,*) ' The theory for this module is outlined in:'
      WRITE(IW,*)
      WRITE(IW,*) ' M. Pernpointner, J. Chem. Phys. 140, 084108 (2014).'
      WRITE(IW,*) ' doi: 10.1063/1.4865964'
      WRITE(IW,*)
      WRITE(IW,*)
      WRITE(IW,*)' ________________________________________________ '
      WRITE(IW,*)'|                                                |'
      WRITE(IW,*)'|       Integral preparation section             |'
      WRITE(IW,*)'|________________________________________________|'
      WRITE(IW,*)

      WRITE (IW,1011) MYPROC+1,NMPROC
 1011 FORMAT (/'---< Process ',I5,' of ',I5,'----<'/)
!
!     Initialize low-level I/O
!
      CALL WAIO_INIT (IW)
!
!     Initialize the timing.
!
      CALL CPUUSED(CPU0)
!
!     Write out memory information and
!     call new mem alloc routines
!
      WRITE (IW,*) 'Memory info:'
      WRITE (IW,*) 'Megawords (MB) given by --mw: ',
     &             MAXCORE/(1000*1000),MAXCORE / (128 * 1024)
      call allocator_get_available_mem(total_mem_new)
      WRITE (IW,*)
      WRITE (IW,*) 'Maximally allowed dynamical memory (MW,B):',
     &      total_mem_new/(8*1000*1000),total_mem_new
      maxcor = nint((total_mem_new/8.0d0),kind(8))
      mword_limit = maxcor / (1024*1024)
!
!     Read symmetry information from MRCONEE
!
      CALL RDSYMI (NSP,EPS1,IREPSP,ncore_elec_adc)
      IF (NSP.GT.MAXSP) THEN
         PRINT*,' INCREASE MAXSP IN POLPRP_MAIN TO',NSP
         CALL QUIT('NSP.GT.MAXSP')
      ENDIF
!
      IPRNT = 2   ! print orbital energies (variable in inpt.inc!)
!
!     Number of integer words for REAL*8 and 
!     number of REAL*8 words per variable (2 for complex groups, 
!     1 for real)
!
      IRW = INTOWP(1)
      IF (CARITH) THEN
         RCW = 2
      ELSE
         RCW = 1
      ENDIF
!
!  Report memory in megabytes
!
!
!     Set up symmetry tables and parallel distribution arrays
!
      CALL SYMTAB (NELEC,NFROZ,NINACT,NACT,NSP,EPS1,EPS,IREPSP,
     &             IREPSPI,.FALSE.,.TRUE.,NELEC_F1,NELEC_F2)
!
!  determine number of electrons and print general information
!
      NELECT = 0
      DO I = 1, 16
         NELECT = NELECT + NELEC(I)
      ENDDO
!
      WRITE (IW,1004) NELECT,NSP-NELECT,CARITH
1004  FORMAT (//' Number of electrons :',T40,I5
     &        /' Number of virtual spinors :',T40,I5
     &        /' Complex arithmetic mode :',T44,L1/)

!
!  initialize USEOE=T default also for POLPRP
      USEOE   = .true.
!  because it relies on CCSETI and needs the orbital
!  energies from the original HF calculation
!  stefan: i am not sure whether you always want to set USEOE = .true.
!          but let us keep it.
!    
!  in principle this allows for the use of X2Cmmf
!  but the NORECMP parameter is NOT irrelevant - stefan march 2014
!
!
      if(X2CMMF)then
!       X2Cmmf hamiltonian in use - by default always set USEOE and NORECMP to .true.
!       everything else results in nonsense energies
        NORECMP = .true.
        useoe   = .true. 
        write(IW,'(/A,/A,/A,/A,/A/)')
     &  ' ************************************************************',
     &  ' *** molecular-mean-field X2C Hamiltonan (X2Cmmf) active  ***',
     &  ' *** use of SCF orbital energies and no re-computation of ***',
     &  ' *** Fock matrix elements set explicitly                  ***',
     &  ' ************************************************************'
      else
        if(USEOE) WRITE (IW,*) 'Orbital energies from HF are used'
      end if

!
!  call integral preparation routine from RELCCSD
!  Sort 2E-Ints into the six classes OOOO,VOOO,VVOO,VOVO,VVVO,VVVV
!  
!       ICALC = 1
!       CALL CCSETI (ICALC,maxcore_info(1),EPS,NSP,IREPSP,IREPSPI)

      ICALC = 2
      CALL CCSETI (ICALC,maxcore_info(1),EPS,NSP,IREPSP,IREPSPI)
      WRITE(IW,*) 'Memory requirements of ',
     &             maxcore_infc(1),maxcore_info(1)
!
!  call *parallel* polarization propagator code
!  
      CALL POLPRP_PAREXE(IW,ITAPADC,EPS,IREPSPI,NSP)
!
!     return after closing low-level I/O files
!
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+0,JCODE)
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+1,JCODE)
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+2,JCODE)
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+3,JCODE)
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+4,JCODE)
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+5,JCODE)
      JCODE = 3
      CALL WAIO_CLOSE(ITAPT+6,JCODE)
!
!     Print I/O statistics
!
      CALL WAIO_PRST (IW)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE POLPRP_PAREXE(IW,IOBASE,EPS,IREPSPI,NSP)
!
! executed in parallel or in serial from here.
!
      use interface_to_mpi
      use polprp_cfg     !these data are visible on master node only!
                         !They contain control input for the POLPRP run
                         !and user input for the Davidson diagonalizer

 
!     IMPLICIT INTEGER(A-Z)
      IMPLICIT NONE
!
!---------------Argument Variables---------------------------------
!
      integer                           :: iw       ! I/O file handle for screen output
      integer                           :: iobase   ! I/O base file handle, all purposes
      real*8, dimension(:)              :: eps
      integer, dimension(:)             :: irepspi  ! spinor definitions
      integer                           :: nsp      ! number of active spinors
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"     !LTR/LFA are defined here (true, false)
#include  "../relccsd/ccpar.inc"
!
!---------------Local variables---------------------------------
!
      integer                 :: n1,n2
      integer                 :: ialloc
      real*8, allocatable     :: daibj(:)
      integer                 :: irploop,desrep
      CHARACTER*6             :: namematfil = 'PP_ELE'  !name of the main ADC matrix file
      CHARACTER*6             :: namecnffil = 'ADCXCF'  !name of the configuration file
      integer                 :: irecl, ireccnt, intbuf
      integer                 :: ladc
      integer                 :: lencouplbl = 0
      integer                 :: nbufs = 0
      integer                 :: signal
      integer                 :: serialindicator
      integer                 :: i
      logical                 :: propertyints
!  
!  variables for rdprop_adc integral testing
!
      integer      ::  iprnt = 1
      character*8  :: nameofprop
      real*8, allocatable   :: a_CBUF(:)
      real*8, allocatable   :: a_EVO(:), a_EOV(:), a_EOO(:), a_EVV(:)
      complex*16            :: ephase
    
!
!-----------------------interface section---------------------
!
      interface

        SUBROUTINE property_handler(i1,i2,l1,i3,ia1)
          integer                  :: i1,i2,i3
          logical                  :: l1
          integer, dimension(:)    :: ia1
        END SUBROUTINE

        SUBROUTINE propinit(i1,i2,l1,i3,ia1)
          integer                  :: i1,i2,i3
          logical                  :: l1
          integer, dimension(:)    :: ia1
        END SUBROUTINE

        SUBROUTINE BUILD_X_MAIN (ra1,ra2,i1)
          real*8, dimension(:)  :: ra1
          real*8, dimension(:)  :: ra2
          integer               :: i1
        END SUBROUTINE

        SUBROUTINE SOLVE_HPHP_R(ra1,i1,i2)
          real*8, dimension(:)  :: ra1
          integer               :: i1,i2
        END SUBROUTINE

        SUBROUTINE SOLVE_HPHP_C(ra1,i1,i2)
          real*8, dimension(:)  :: ra1
          integer               :: i1,i2
        END SUBROUTINE

        SUBROUTINE BUILD_X_COUSAT_R(h,a,i1,i2,i3,i4,l1,r2,i5)
          real*8, dimension(:)  :: h
          real*8, dimension(*)  :: a
          integer               :: i1,i2,i3,i4
          logical               :: l1
          real*8                :: r2
          integer               :: i5
        END SUBROUTINE

        SUBROUTINE BUILD_X_COUSAT_C(h,a,i1,i2,i3,i4,l1,r2,i5)
          real*8, dimension(:)  :: h
          real*8, dimension(*)  :: a
          integer               :: i1,i2,i3,i4
          logical               :: l1
          real*8                :: r2
          integer               :: i5
        END SUBROUTINE

        SUBROUTINE mpi_bxcs_r(ra1,i1,ra2,i2,i3,i4,i5,l1,r1,i6,i7)
          real*8, dimension(:)  :: ra1
          integer               :: i1,i2,i3,i4,i5,i6,i7
          real*8, dimension(*)  :: ra2
          real*8                :: r1
          logical               :: l1
        END SUBROUTINE

        SUBROUTINE mpi_bxcs_c(ra1,i1,ra2,i2,i3,i4,i5,l1,r1,i6,i7)
          real*8, dimension(:)  :: ra1
          integer               :: i1,i2,i3,i4,i5,i6,i7
          real*8, dimension(*)  :: ra2
          real*8                :: r1
          logical               :: l1
        END SUBROUTINE

        SUBROUTINE mpi_davidson_master_r(i1,i2,i3,i4,c1,c2,i5,i6,i7)
          integer               :: i1,i2,i3,i4,i5,i6,i7
          character*6           :: c1,c2
        END SUBROUTINE

        SUBROUTINE mpi_davidson_master_c(i1,i2,i3,i4,c1,c2,i5,i6,i7)
          integer               :: i1,i2,i3,i4,i5,i6,i7
          character*6           :: c1,c2
        END SUBROUTINE

        SUBROUTINE mpi_xmat_dist(i1,c1,i2,i3,i4,i5,i6)
          integer               :: i1,i2,i3,i4,i5,i6
          character*6           :: c1
        END SUBROUTINE

        SUBROUTINE mpi_davidson_helper_r (i1,c1,i2)
          integer               :: i1,i2
          character*6           :: c1
        END SUBROUTINE

        SUBROUTINE mpi_davidson_helper_c (i1,c1,i2)
          integer               :: i1,i2
          character*6           :: c1
        END SUBROUTINE

        SUBROUTINE build_tramo(i1,i2,ra1)
          integer               :: i1,i2
          real*8, dimension(:)  :: ra1
        END SUBROUTINE

      end interface
!
!---------------Execution ---------------------------------
!

!  _________________________ INTRO SECTION __________________________
! |
! |
! |
#if defined (VAR_MPI)
      if(MYPROC.EQ.MASTER) then   !in parallel: only master prints info
#endif

        WRITE(IW,*)' ________________________________________________ '
        WRITE(IW,*)'|                                                |'
        WRITE(IW,*)'|           Input control section                |'
        WRITE(IW,*)'|________________________________________________|'
        WRITE(IW,*)
        CALL print_input_info(iw,nrep) 

#if defined (VAR_MPI)
      endif
      WRITE(IW,'(A,I5)') ' Total # of nodes: ',NMPROC
      WRITE(IW,'(A,I5)') ' Current node # :  ',MYPROC
      
#endif

! |
! |
! |_________________________ END INTRO SECTION __________________________



!  _________________________ CALCULATION SECTION __________________________
! |
! |  in parallel: only master calcs ADC mat (is very fast),
! |  the nodes help with various
! |  tasks. Nodes will massively participate in Davidson.
! |
#if defined (VAR_MPI)
      if(MYPROC.EQ.MASTER) then 
#endif

        WRITE(IW,*)' ________________________________________________ '
        WRITE(IW,*)'|                                                |'
        WRITE(IW,*)'|           Calculation section                  |'
        WRITE(IW,*)'|________________________________________________|'
        WRITE(IW,*)
!
!  check if there are transformed dipole moment integrals available
!  that are required for transition moment calculations.
!  If yes fetch them and adapt them for further processing.
!
        CALL property_handler(iobase,iobase+1,propertyints,NSP,IREPSPI)
!
!  do HPHP block first (it can easily be kept in memory)
!
!   all symmetries in the HP-HP block are created at once.
!   Allocate hphp block memory
!

        n1 = jvovo(nrep+1)
        allocate(daibj(n1*rcw),stat=ialloc)
        if(ialloc.ne.0) stop 'allocation problem daibj'
        daibj=0.0d0

        CALL BUILD_X_MAIN (daibj,eps,n1) !in matrix storage mode: <||>,|><|
!
! diagonalize main block for Davidson start vectors
!
        if(carith) then
          CALL SOLVE_HPHP_C(daibj,n1,iobase)   !iobase is reusable!
        else
          CALL SOLVE_HPHP_R(daibj,n1,iobase)   !iobase is reusable!
        endif
#if defined (VAR_MPI)
      endif
#endif
!__________________________________________________
!|
!|      this part is executed by all nodes
!|

        INTBUF = 200*1024   !buffer size for matrix R/W buffers
!|
!|_________________________________________________


#if defined (VAR_MPI)
      if(MYPROC.EQ.MASTER) then 
#endif
! ______________ PRECONTRACTION of the VVVV and the VVOO integrals
! |
! |   In case of parallel execution we have to collect the VVVV
! |   integrals from the nodes. Otherwise all are on one node.
! |   Files will be written to 'VVVVXVVOO'

        CALL PRECONTRACT(iobase,iw)

! ______________ SYMMETRY LOOP FOR COUPLING + SAT BLOCK___________
! |
! |   In case of parallel execution 
! |   only master loops through the irreps and the slaves provide 
! |   helper routines.
! |

        DO IRPLOOP = 1,NREP

! check for valid inputs

          desrep = polprp_statesym(IRPLOOP)
          if(desrep.EQ.0) cycle
          if(desrep.GT.NREP) then
             WRITE(IW,*) 'Chosen Symmetry',desrep,
     &       ' does not exist and is skipped!'
             CYCLE
          endif
          if(MVO(desrep).EQ.0) then
            WRITE(iw,*) 'No states in symmetry',desrep,repna(desrep)
            CYCLE
          endif

! print current run info

          CALL PST('Starting calculation of Coupling and SAT block+')
          WRITE(IW,'(A,I4,X,A4)') 'Curr. symmetry:',desrep,repna(desrep)
          WRITE(IW,'(A,A)') 'Used main matrix file name: ',namematfil
          WRITE(IW,'(A,I14)') 'Current I/O buffer size: ',INTBUF
          WRITE(IW,'(A,I14)') 'Real/complex word: ',rcw
          WRITE(IW,*)

! generate configuration file for this symmetry

          WRITE(IW,*) 'Configuration data written to ',namecnffil
          CALL WCONDAT_X(iobase,namecnffil,irecl,desrep,ireccnt)
          WRITE(IW,*) 'Record length in configuration file:',IRECL
          WRITE(IW,*) 'Number of records (=matrix dimension):',ireccnt

! generate coupling and sat block in this symmetry
! together with the helpers simultaneously started.
! open main matrix file and construct main matrix.

#if defined (VAR_MPI)
!
!  send the slave a start signal only in case if there *is* something to do
!
          signal = 1  !start code.
          write(iw,*) 'Master indicates start.'
          call interface_mpi_bcast(signal,1,MASTER,
     &                             global_communicator)

          open(unit=iobase,file=namematfil, access='SEQUENTIAL',
     &         form='UNFORMATTED',status='NEW')
          n2 = jvovo(nrep+1)
          if(carith) then
            call mpi_bxcs_c(daibj,
     &                      n2,
     &                      eps,
     &                      desrep,
     &                      lencouplbl,
     &                      intbuf,
     &                      nbufs,
     &                      polprp_doextended,
     &                      polprp_writethr,
     &                      iobase,
     &                      iw)
          else
            call mpi_bxcs_r(daibj,
     &                      n2,
     &                      eps,
     &                      desrep,
     &                      lencouplbl,
     &                      intbuf,
     &                      nbufs,
     &                      polprp_doextended,
     &                      polprp_writethr,
     &                      iobase,
     &                      iw)
          endif
          close(iobase)  !after distribution matrix will be deleted.

!                           lencouplbl and nbufs are set in the
!                           mpi_bxcs_r routine!
#else
!
! remind: in the serial case the number of buffers (NBUFS) is known
!         after matrix construction.
!
        OPEN(UNIT=IOBASE,FILE=namematfil, ACCESS='SEQUENTIAL',
     &     FORM='UNFORMATTED',STATUS='NEW')
        serialindicator = -1
        Write(iobase) serialindicator
        IF(CARITH) THEN
          CALL BUILD_X_COUSAT_C(DAIBJ,EPS,DESREP,LENCOUPLBL,
     &                          INTBUF,NBUFS,polprp_doextended,
     &                          polprp_writethr,iobase)
          CALL XDIAG_INI_C(iobase,iobase+1,
     &                     INTBUF,NBUFS,lencouplbl + MVO(desrep))
        ELSE
          CALL BUILD_X_COUSAT_R(DAIBJ,EPS,DESREP,LENCOUPLBL,
     &                          INTBUF,NBUFS,polprp_doextended,
     &                          polprp_writethr,iobase)
          CALL XDIAG_INI_R(iobase,iobase+1,
     &                     INTBUF,NBUFS,lencouplbl + MVO(desrep))
        ENDIF
        CLOSE(iobase)

        Write(*,*)
        Write(*,*) 'Serial run information:'
        Write(*,*) 'ADC matrix has',NBUFS,' buffers of length',INTBUF
        Write(*,*)
#endif

!
! check for internal consistency (we use various ways for deriving the matrix dimension)
! 
          ladc = lencouplbl + MVO(desrep)
          if(ladc.ne.ireccnt) then
            call quit(' *** Internal error. ADC matrix dim. mismatch')
          endif

#if defined (VAR_MPI)
! Main ADC matrix built in this symmetry. We synchronize all nodes

          call interface_mpi_BARRIER(global_communicator)
          write(iw,*) 'Master node synchronized after',
     &               ' matrix construction.'

! slaves have knowledge only of the variables 'iobase', 'iw' and 'namematfil'
! all other relevant variables for matrix distribution have to be bcast.

          call interface_mpi_bcast(irploop,1,MASTER,global_communicator)
          call interface_mpi_bcast(intbuf,1,MASTER,global_communicator)
          call interface_mpi_bcast(nbufs,1,MASTER,global_communicator)
          call interface_mpi_bcast(ladc,1,MASTER,global_communicator)

!  enter matrix distribution routine (master entry)

          call mpi_xmat_dist(iobase,namematfil,irploop,
     &                     intbuf,nbufs,ladc,iw)

! We synchronize all nodes after matrix distribution

          call interface_mpi_BARRIER(global_communicator)
          write(iw,*) 'Master node synchronized after',
     &               ' matrix distribution.'

! Master calls Davidson diagonalizer

          if(carith) then
            call fulldia_c3(iobase,namematfil,intbuf,ladc)
            call mpi_davidson_master_c(iw,
     &                                 iobase,
     &                                 desrep,
     &                                 intbuf,
     &                                 namematfil,
     &                                 namecnffil,
     &                                 irecl,
     &                                 ladc,
     &                                 MVO(desrep))
          else
            call fulldia_r3(iobase,namematfil,intbuf,ladc)
            call mpi_davidson_master_r(iw,
     &                                 iobase,
     &                                 desrep,
     &                                 intbuf,
     &                                 namematfil,
     &                                 namecnffil,
     &                                 irecl,
     &                                 ladc,
     &                                 MVO(desrep))
          endif

! Master deletes its main ADC matrix file

          open(unit=iobase,file=NAMEMATFIL, access='SEQUENTIAL',
     &         form='UNFORMATTED',status='unknown')
          close(iobase,status='DELETE')
#else
!
!________________________________________________________________
!|   
!|     serial case diagonalizer call
!|   
!
          if(carith) then
            call fulldia_c3(iobase,namematfil,intbuf,ladc)
            call mpi_davidson_master_c(iw,
     &                                 iobase,
     &                                 desrep,
     &                                 intbuf,
     &                                 namematfil,
     &                                 namecnffil,
     &                                 irecl,
     &                                 ladc,
     &                                 MVO(desrep))
          else
            call fulldia_r3(iobase,namematfil,intbuf,ladc)
            call mpi_davidson_master_r(iw,
     &                                 iobase,
     &                                 desrep,
     &                                 intbuf,
     &                                 namematfil,
     &                                 namecnffil,
     &                                 irecl,
     &                                 ladc,
     &                                 MVO(desrep))
          endif

          open(unit=iobase,file=NAMEMATFIL, access='SEQUENTIAL',
     &         form='UNFORMATTED',status='unknown')
          close(iobase,status='DELETE')
#endif
!
!_________________________________________________________________
!|
!| serial or parallel/MASTER: calculate transition moments
!|

          if(propertyints) then
            write(*,*)
            write(*,*)  'Property integrals exist.'
            write(*,*)  'TM are calculated for this state.'
            write(*,*)
            write(*,*)  'Attention: temporarily commented out!'
            write(*,*)
!           call build_tramo(desrep,iobase,eps)
          else
            write(*,*) 'No transition moments for irep',desrep,
     &                 ' can be calculated.'
          endif
!|
!|________________________________________________________________

        ENDDO   !IRPLOOP

#if defined (VAR_MPI)
        signal = -1  !code for slave loop termination
        write(iw,*) 'Master indicates loop break.'
        call interface_mpi_bcast(signal,1,MASTER,
     &                           global_communicator)
#endif
        deallocate(daibj)
! |
! |
! |________________ END SYMMETRY LOOP (MASTER) _________________________

!   ________________________________________________________________
!   |
!   |
!   |     Slave execution part
!   |
!   |     Consists of corresponding helper routines mainly for heavy duty
!   |     distributed matrix multiplication
!   |
#if defined (VAR_MPI)
      else

        write(*,*) 'Node ',myproc,' is waiting for master finalizing',
     &               ' HPHP block.'

        Do while(.true.)    ! Circulating in the helper loop.
                            ! Slave gets a start signal and a finalizing signal.

          call interface_mpi_bcast(signal,1,MASTER,
     &                             global_communicator)
          if(signal.eq.-1) then
            write(iw,*) 'Node',MYPROC,' received termination signal.'
            exit  !Leave do loop
          else if (signal.eq.1) then
            write(iw,*) 'Node',MYPROC,' received start sig. ADC matrix.'
          else if (signal.eq.2) then
            write(iw,*) 'Node',MYPROC,' received start sig. PRECONTR.'
          else
            call quit('Illegal signal!')
          endif
!
!  check if precontraction is requested (only once per program execution!)
!  In this case call the int_supplier because
!  nothing else is needed for this task. Afterwards we have to reenter main cycle
!  WITHOUT DOING ANYTHING ELSE for the individual tasks per IRREP. 
!
          if(signal.eq.2) then
            call mpi_int_supplier(iw)
            cycle
          endif
!
!  call the integral completion routine for the slaves
!  serves real and complex case.
!
          call mpi_int_supplier(iw)

! Main ADC matrix built in this symmetry. We synchronize all nodes...
! ... and prepare for matrix distribution.

          call interface_mpi_BARRIER(global_communicator)
          write(*,*) 'Node',MYPROC,' synchronized after',
     &               ' matrix construction.'

! fetch unknown variables from master node

          call interface_mpi_bcast(irploop,1,MASTER,global_communicator)
          call interface_mpi_bcast(intbuf,1,MASTER,global_communicator)
          call interface_mpi_bcast(nbufs,1,MASTER,global_communicator)
          call interface_mpi_bcast(ladc,1,MASTER,global_communicator)

! and enter matrix distribution routine (slave entry)

          call mpi_xmat_dist(iobase,namematfil,irploop,
     &                     intbuf,nbufs,ladc,iw)

! We synchronize all nodes after matrix distribution

          call interface_mpi_BARRIER(global_communicator)
          write(iw,*) 'Node',MYPROC,' synchronized after',
     &               ' matrix distribution.'

! We enter the helper for Davidson diagonalisation
! required variables are known to the slaves.

          if(carith) then
            write(iw,*) 'Node',MYPROC,' enters (complex) dav. helper.'
            call mpi_davidson_helper_c(iw,namematfil,iobase)
          else
            write(iw,*) 'Node',MYPROC,' enters (real) dav. helper.'
            call mpi_davidson_helper_r(iw,namematfil,iobase)
          endif

        Enddo  !end of infinite helper loop

        write(iw,*) 'Node',MYPROC,' will return to main program' 

      endif
!   |
!   |
!   |   End of slave execution part
!   |_______________________________________________________________

#endif

! |
! |
! |_________________________ END CALCULATION SECTION __________________________

      return
      end 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE print_input_info(iw,nrep)
      use polprp_cfg     !these data are visible on master node only!

      implicit none

      integer               :: iw         ! I/O file handle for proper screen output
      integer               :: nrep       ! number of total symmetries in the calculation
!
      integer               :: i,n,icount

      write(iw,*) 'POLPRP control variables:'
      write(iw,*) '-------------------------'
      write(iw,*)
      write(iw,*) 'we do extended calculation: ',polprp_doextended
      write(iw,*) 'threshold for writing elements to ADC matrix:',
     &             polprp_writethr
      write(iw,*) 'we perform a Davidson diagonalization:',polprp_dodiag
      write(iw,*) 'The way final state symmetries are chosen is:'
!
! determine symmetries to be calculated
!
      n = 0
      do i = 1,nrep
        n = n + polprp_statesym(i)
      enddo
      if(n.eq.0) then
        write(iw,*) 'Program checks symmetries automatically'
        do i=1,nrep
          polprp_statesym(i) = i
        enddo
      else
        write(iw,*) 'User has selected the following final'
        write(iw,*) 'state symmetries:'
        icount=1
        do i=1,nrep
          if(polprp_statesym(i).ne.0) then
            write(iw,*) '   No: ',icount,
     &                  '  symmetry: ',polprp_statesym(i)
            icount = icount + 1
          endif
        enddo
      endif
      write(iw,*)
!
! report diagonalization parameters
!
      write(iw,*) 'Davidson control variables:'
      write(iw,*) '---------------------------'
      write(iw,*)

      if(polprp_dodiag) then
        WRITE(IW,*) 'Number of eigenstates to be calculated (roots)',
     &     polprp_davroots
        WRITE(IW,*) 'Size of Krylov space: ',polprp_davmaxsp
        WRITE(IW,*) 'Number of Davidson macro iterations: ',
     &     polprp_davmaxit
        WRITE(IW,'(A,E20.4)') ' Requested energy conv.: ',polprp_davconv
        IF(polprp_davooc.eqv..false.) THEN
           WRITE(IW,*) 'Trying full incore storage of Davidson vectors.'
           WRITE(IW,*) 'Will default to out-of-core if resources are'
           WRITE(IW,*) 'not available.'
        ENDIF
      endif

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRECONTRACT(iounit,iw)
!
!
!  this routine is called by master only when parallel execution takes place.
!  Then the nodes provide integrals.
!
!

#if defined (VAR_MPI)
      use interface_to_mpi
#endif

      implicit none

      integer               :: iounit         ! I/O file handle for writing to file
      integer               :: iw             ! file handle for writing to screen
!
! local variables
!
      real*8, allocatable, dimension(:)        :: vvvvz
      real*8, allocatable, dimension(:)        :: vvoo1
      real*8, allocatable, dimension(:)        :: vvoop
      integer                                  :: n1,ixx
      integer                                  :: ialloc
      integer                                  :: signal
      integer*8                                :: n8,i8
      COMPLEX*16                               :: A0,A1
      character*9                              :: preconfn='VVVVXVVOO'


#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/ccpar.inc"
!
! execution
!
      if(carith) then
        CALL PST('Performing (complex) precontraction+')
        if(rcw.ne.2) call quit ('intenal error in PRECONTRACT')
      else
        CALL PST('Performing (real) precontraction+')
        if(rcw.ne.1) call quit ('intenal error in PRECONTRACT')
      endif
      A0=(0.0D00,0.0D00)
      A1=(1.0D00,0.0D00)


#if defined (VAR_MPI)
!
!  send the slaves a start signal for int supply only
!
      signal = 2  !start code.
      write(iw,*) 'Master indicates PRECONTRACTION (signal=2)'
      call interface_mpi_bcast(signal,1,MASTER,
     &                             global_communicator)
#endif

! allocate buffers and read in integrals

      n1 = ivvoott(nrep+1)
      allocate(vvoo1(n1*rcw))
      allocate(vvoop(n1*rcw))

      n8 = ivvvvtt(nrep+1)
      allocate(vvvvz(n8*rcw),stat=ialloc)
      if(ialloc.ne.0) then
        STOP 'VVVV buffer could not be allocated.'
      else
        write(*,*) 'Master allocated',n8,' words for VVVV buffer.'
      endif

#if defined (VAR_MPI)
      do i8=1,n8*rcw
        vvvvz(i8) = -4096.0d0
      enddo
      call rdvvvv(vvvvz)  ! master reads its own chunks
      write(*,*) 'masters chunk reading ... complete.'
      call mpi_vvvv_complete(vvvvz,n8*rcw)  !.... and completes the stream.
      write(*,*) 'VVVV check in PRECONTRACT...'
      do i8=1,n8
        if(vvvvz(i8).eq. -4096.0d0)
     &    call quit('Internal error: gap in the VVVV stream!')
      enddo
      write(*,*) '.... success. No leaks.'
#else
      call rdvvvv(vvvvz)  ! master reads complete stream in serial mode
#endif
      call getvvoo(vvoo1)
!
!  do the contraction over all symmetries
!
      CALL CNTRCT('N','N',NVVT,NOOT,NVVT,A1,vvvvz,vvoo1,A0,vvoop,NREP)
      write(*,*) 'PRECONTRACTION finished.'
!
!  open corresponding file and write out contracted integrals.
!
      open(unit=iounit,file=preconfn, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='NEW')
      write(iounit) (vvoop(ixx), ixx=1,n1*rcw)
      close(iounit)
!
!  free buffers immediately 
!
      deallocate(vvoop)
      deallocate(vvoo1)
      deallocate(vvvvz)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
