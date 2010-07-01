! $Id: diremb.F 9568 2009-09-15 09:50:44Z gomes $
!
! description : this subroutine will read an embedding potential from 
!               file and construct a matrix representation of that 
!               potential in the ao basis currently in use, via numerical
!               integration. 
!
!               the routine expects that the grid file has the following
!               format:
!
!               line        one : number of grid point
!               lines to to n-1 : quadruplet of numbers, where the
!                                 first three are the (x,y,z) coords
!                                 and the last is the value of the
!                                 embedding potential at that point
!                                 times the weight of the quadrature
!                                 [= v(i)*w(i)]
!               line          n : a non-negative integer, like  -42 =)
!
!               the default name for the grid file is EMBPOT, but that
!               can be overriden via a menu option. see routine fdeinp_
!               for further details
!
!
      MODULE FDE_MOD
      USE GRID_FILE
#ifdef MOD_XML
      USE FILE, ONLY: OPEN_FILE
#endif
#ifdef MOD_XCFUN
      use xcfun
#endif
      use dft_input
      use dft_cfg

      PRIVATE
!     PUBLIC GetFDEPotential_FromFile, FDEINI,FDEINP,FDE_GRIDOUT
      PUBLIC FDE_ImportStaticPotential_FileOpen, FDE_ImportStaticPotential_FileClose
      public FDE_setupFDEgrid, FDE_testFrozendensity, FDE_cleanupFDEgrid
      public FDEINI,FDEINP,FDE_GRIDOUT,have_frozen_density
      PUBLIC DO_GRID_OUT
      public fde_grid, fde_frozen_prp, fde_grid_np, fde_nprop

      LOGICAL,SAVE :: DO_GRID_OUT=.FALSE.
      logical,save :: have_frozen_density = .false.
      logical,save :: import_vemb = .false.
      logical,save :: update_vemb = .false.
      INTEGER,SAVE :: LUFDEPOT=45,LUTESTING=46,LUGRIDOUT=47, LUFDEFRZDNS=48
      INTEGER,SAVE :: print_fde=0, fde_grid_np=0, fde_nprop=1
      CHARACTER(LEN=60),SAVE :: FDEPOT='EMBPOT',FDEDNS='FRZDNS',GRIDOUT_NAME=''
      CHARACTER(LEN=60),SAVE :: out_level='HF'
      real(kind=8), pointer :: fde_grid(:,:), fde_frozen_prp(:,:) 
      
      CONTAINS

         subroutine FDE_ImportStaticPotential_FileOpen(embpot_unit)
            integer, intent(out) :: embpot_unit
            logical              :: embpot_found

            inquire(FILE=FDEPOT,EXIST=embpot_found)

            if (embpot_found) then
                open(lufdepot,               &
                     FORM   = 'FORMATTED',   &
                     STATUS = 'UNKNOWN',     &
                     ACCESS = 'SEQUENTIAL',  &
                     FILE   =  FDEPOT)
                embpot_unit = lufdepot
             else
                call quit("File with embedding potential not found!")
             endif

         end subroutine FDE_ImportStaticPotential_FileOpen

         subroutine FDE_ImportStaticPotential_FileClose()
            close(lufdepot) 
         end subroutine FDE_ImportStaticPotential_FileClose


         subroutine FDE_ImportFrozenDensity_FileOpen(frozdens_unit)
            integer, intent(out) :: frozdens_unit
            logical              :: frozdens_found

            inquire(FILE=FDEDNS,EXIST=frozdens_found)

            if (frozdens_found) then
                open(lufdefrzdns,               &
                     FORM   = 'FORMATTED',   &
                     STATUS = 'UNKNOWN',     &
                     ACCESS = 'SEQUENTIAL',  &
                     FILE   =  FDEDNS)
                frozdens_unit = lufdefrzdns
             else
                call quit("File with Frozen Density(+gradient) not found!")
             endif

         end subroutine FDE_ImportFrozenDensity_FileOpen

         subroutine FDE_ImportFrozenDensity_FileClose()
            close(lufdefrzdns) 
         end subroutine FDE_ImportFrozenDensity_FileClose


         subroutine FDE_setupFDEgrid()
            integer :: frozdens_unit

            call FDE_ImportFrozenDensity_FileOpen(frozdens_unit)
            call read_grid(frozdens_unit,fde_grid,fde_frozen_prp)
            call FDE_ImportFrozenDensity_FileClose()

            fde_grid_np = size(fde_grid,2)
            fde_nprop   = size(fde_frozen_prp,1)
            
         end subroutine FDE_setupFDEgrid

         subroutine FDE_testFrozenDensity()
            integer      :: point = 0
            real(kind=8) :: frozen_particle_number = 0.0

            do point = 1, fde_grid_np 
               frozen_particle_number = frozen_particle_number + fde_grid(4,point)*fde_frozen_prp(2,point) 
            end do
            
            write (*,*) '  Test: imported frozen density integrates to ',frozen_particle_number
            write (*,*) '  Test: grid points processed ',fde_grid_np
            write (*,*) '  Test: number of properties  ',fde_nprop
            frozen_particle_number = 0.0

         end subroutine FDE_testFrozenDensity

         subroutine FDE_cleanupFDEgrid()
            deallocate(fde_grid)
            deallocate(fde_frozen_prp) 
         end subroutine FDE_cleanupFDEgrid

      SUBROUTINE FDEINI()

      END SUBROUTINE


      SUBROUTINE FDEINP(WORD,Echo_Input)
#include "implicit.h"
#include "priunit.h" 
#include "dcbgen.h"
      PARAMETER (NTABLE = 8)
!
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7, line*80
!
!
      logical :: Echo_Input,newdef
!
#include "memint.h"

      print_fde=0
      ICHANG = 0
      IF (WORD(1:4).eq.'*FDE') THEN
        WORD1 = WORD
        DO
          READ (LUCMD, '(A7)') WORD
          SELECT CASE (WORD(1:1))
          CASE('!','#')
            CYCLE
          CASE('.')
            ICHANG = ICHANG + 1
            SELECT CASE(WORD(1:7))
            CASE ('.OPTION')
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
              CYCLE
            CASE ('.EMBPOT')
              import_vemb = .true.
!             read filename for the embedding potential
              READ(LUCMD,'(A60)') FDEPOT 

            CASE ('.UPDATE')
              update_vemb = .true.
!             read filename for the embedding potential
              write(*,*) 'Embedding potential will be updated and saved'

            CASE ('.FRDENS')
!             read filename for the frozen density and density gradient
              READ(LUCMD,'(A60)') FDEDNS
              write(*,*) 'File with frozen density (+gradient etc) (=fdedns) is ',fdedns
              have_frozen_density = .true.

! read kinetic energy functional for non-additive term
            CASE ('.OFKEF')
#ifdef MOD_XCFUN
              read(lucmd, '(a80)') line
              call parse_functional(line, id_kin_fun)
              call set_kin_fun_alda(id_kin_fun_alda)
#else         
              call quit('XCFun not available in this version')
#endif /* ifdef MOD_XCFUN */

            CASE ('.PRINT')
!             set printlevel 
              READ(LUCMD,*) print_fde 
            CASE ('.LEVEL')
              READ(LUCMD,'(A60)') out_level
              select case(out_level)
              case('HF','MP2')
                 WRITE(LUPRI,*) 'Outputted density will be taken from: ' &
                      ,trim(out_level)
              case default
                 WRITE(LUPRI,*) 'Input for outputted density not a &
     &  x  valid calculation, HF will be used instead. Given:',out_level
              end select
            CASE ('.GRIDOU')
! WRITE TO THIS GRIDFILE
              READ(LUCMD,'(A60)') GRIDOUT_NAME
              DO_GRID_OUT=.TRUE.
              DOPRP=.TRUE.
            CASE DEFAULT
              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD, &
                 '" not recognized in FDEINP.'
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
              CALL QUIT('Illegal keyword in FDEINP.')
            END SELECT
          CASE('*')
            EXIT
          CASE DEFAULT
            WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD, &
               '" not recognized in FDEINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal prompt in FDEINP.')
          END SELECT
        END DO
      END IF

!     options are echoed
      IF(Echo_Input) THEN
        write (lupri,*) ' '
        CALL PRSYMB(LUPRI,'=',101,1)
        WRITE(LUPRI,'(1X,2A,/1X,2A)') &
           'Frozen Density Embedding to include enviroment effects', &
           ' using the fixed potential method described in:',        &
           'A.S.P. Gomes, C. R. Jacob and L. Visscher, ',            &
           'Phys. Chem. Chem. Phys. 10 (2008) 5353–5362.'
        CALL PRSYMB(LUPRI,'=',101,1)
        WRITE(LUPRI,*)
        WRITE(LUPRI,'(1X,A)') &
           ' Frozen Density Embedding Calculation set-up:'  
        CALL PRSYMB(LUPRI,'-',44,2)
        WRITE(LUPRI,*)
        WRITE(LUPRI,'(A,I5)')  ' * Print level: ',print_fde

        if (import_vemb) then
           write(LUPRI,'(A,A60)')  ' * FDE Potential read from file: ',FDEPOT
        else
           write(LUPRI,'(A,A60)')  ' * FDE Potential (re)generated internally'
           if (.not.have_frozen_density) then
              write(LUPRI,'(A,A60)')  '   Data from frozen subsystems unavailable!'
              call quit("Missing datafile FRDENS")
           endif
        endif

        if (update_vemb) then 
           write(LUPRI,'(A,A60)') ' * Embedding potential will be updated and saved'
           if (.not.have_frozen_density) then
              write(LUPRI,'(A,A60)')  '   Data from frozen subsystems unavailable!'
              call quit("Missing datafile FRDENS")
           endif
        else
           write(LUPRI,'(A,A60)') ' * Fixed embedding potential will be used' 
        endif

        if (have_frozen_density) then
           write(LUPRI,'(A,A60)') &
           ' * Density (and gradient) from frozen subsystems read from file: ',FDEDNS
           call FDE_setupFDEgrid()
        endif

        IF (DO_GRID_OUT) then
           write(LUPRI,'(2A)') &
           ' * FDE Gridfile with updated density written to:',GRIDOUT_NAME
        endif

      END IF

      END SUBROUTINE


      SUBROUTINE TIME()
!     Simple timing function
#include "priunit.h"
      INTEGER :: file
      
      INTEGER,save :: time_val1(8)=0
      INTEGER      :: time_val2(8), delta_time
      LOGICAL      :: test
!     For timing
      CHARACTER(LEN=12) :: date_tmp,time_tmp,zone_tmp

      test=.FALSE.
      DO i=1,8; if(.NOT.(time_val1(i)==0)) test=.TRUE.
      END DO

      IF (.NOT.(test)) THEN
        CALL DATE_AND_TIME(date_tmp,time_tmp,zone_tmp,time_val1)
        RETURN
      END IF

      CALL DATE_AND_TIME(date_tmp,time_tmp,zone_tmp,time_val2)

      delta_time=1*(time_val2(7)-time_val1(7)) &
               +60*(time_val2(6)-time_val1(6)) &
               +60*(time_val2(5)-time_val1(5)) &
             +3600*(time_val2(4)-time_val1(4))  

      WRITE(LUPRI,*) 'Time spent in embedding output routine: ', & 
           delta_time, ' s.'
      WRITE(LUPRI,*) 'Time spent in embedding output routine: ', &
           time_val2(4)-time_val1(4), ' d.', &
           time_val2(5)-time_val1(5), ' h.', &
           time_val2(6)-time_val1(6), ' m.', &
           time_val2(7)-time_val1(7), ' s.' 
      time_val1(:)=time_val2(:)

      END SUBROUTINE

#ifdef MOD_XML
       SUBROUTINE FDE_GRIDOUT(WORK,LWORK)
      use memory_allocator
      USE VC_GRID
      use xml_parser
!     use conversion
!      IMPLICIT NONE

#include "priunit.h"
#include "dcbbas.h"
#include "maxaqn.h"
#ifdef MOD_DNF
#include "densfit.h"
#endif
#include "onecom.h"
#include "ccom.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbgen.h"
!     Passed variables
      REAL(KIND=8) :: WORK(*)
      INTEGER :: LWORK

      type(xml_tag),pointer :: tag

!     Local variables
      character(len=132) :: csizept1, csizept2
      character(len=132) :: csizegrdrho1, csizehesrho1
      REAL(KIND=8) :: DMAT(NTBAS(0),NTBAS(0),4)
      REAL(KIND=8) :: CMO(NTBAS(0),NTBAS(0),4)
      REAL(kind=8), pointer :: points(:,:)
      REAL(kind=8), allocatable, target :: grad_rho(:,:)
      REAL(kind=8), allocatable, target :: hes_rho(:,:),vc(:),rho(:,:)

      INTEGER ::i,ierr
      INTEGER :: NDIMENSION
!     JHS THROW AWAY IN THE END?
      INTEGER :: NDER=0, DOLND=0

      REAL(KIND=8),ALLOCATABLE :: GAO(:), GAO1(:),GAO2(:),BUF(:),NCNT(:)
      REAL(KIND=8),ALLOCATABLE :: GAB1(:)

      REAL(KIND=8)::dummy,TOLS,TOLOG
      REAL(KIND=8)::gexp=0.0,ncentc=1.0,factor=1.0
      
      CALL TITLER('Output of density to XML file','*',118)
      WRITE(LUPRI,*) 'Name of the outputfile:', trim(GRIDOUT_NAME)

      CALL TIME()

      NDIMENSION=NTBAS(0)

#ifdef MOD_DNF
      IF(DF_FITVXC) NDIMENSION=MAX(NTBAS(0), NBASISAUX)
#endif

      call alloc(gao ,   ndimension)
      call alloc(gao1, 3*ndimension)
      call alloc(gao2, 6*ndimension)
      call alloc(gab1, 3*ndimension)
      call alloc(ncnt,   ndimension)
      call alloc(buf , 4*ndimension)

!     Initialize the matrices to zero.

      GAO  = 0.0
      GAO1 = 0.0
      GAO2 = 0.0
      BUF  = 0.0
      NCNT = 0.0
      GAB1 = 0.0
      DMAT = 0.0


!     JHS Construct the density matrix.
!     
!     For the case of the SCF density:
!     Get CMO from file.
      CMO(:,:,:)=0.0

      select case(out_level)
      case('MP2')
         write(lupri,*) 'Reading CC_DENSITY'
         call GET_CC_DENSITY(out_level,'AO',DMAT,ierr)
         if (ierr.ne.0) then
            write(lupri,*) 'Error in reading density matrix'
            call quit('error readin density matrix')
         end if
      case default
!     HF density used
         write(lupri,*) 'Reading DFCOEF'
         CALL OPNFIL(LUCOEF,'DFCOEF','UNKNOWN','EMBDRV')
         CALL REACMO(LUCOEF,'DFCOEF',CMO,dummy,dummy,dummy,3,WORK,1, &
              LWORK)
         if (print_fde>20) then
            WRITE(LUPRI,*) 'CMO', CMO
            WRITE(LUPRI,*) 'DMAT OLD', DMAT
         end if
!     SCF construct density matrix
         CALL DENMAT(DMAT,CMO,print_fde)
!     Correct DMAT for occupationnumbers
         DMAT(:,:,:)=2*DMAT(:,:,:)

!     TMP scaling
      end select

!     GET an id for LUGRIDOUT.
      LUGRIDOUT=OPEN_FILE(GRIDOUT_NAME)
      IF (LUGRIDOUT<0) &
           CALL QUIT('ERROR OPENING GRIDOUT FILE. FILE: ' &
           //GRIDOUT_NAME)
      IF (LUGRIDOUT==0) &
           CALL QUIT('FILE FOR GRIDOUT DOES NOT EXISTS. FILE: ' &
           //GRIDOUT_NAME)

      WRITE(LUPRI,*) 'Reading grid.'
      CALL READ_GRID(LUGRIDOUT,points)
      CLOSE(LUGRIDOUT)
      WRITE(LUPRI,*) 'Grid read from file.'

      write(lupri,*) 'Number of gridpoints:', size(points,2)

      call alloc(rho,     1,size(points,2))
      call alloc(grad_rho,3,size(points,2))
      call alloc(hes_rho, 6,size(points,2))
      call alloc(vc,      size(points,2))

      rho      = 0.0d0
      grad_rho = 0.0d0
      hes_rho  = 0.0d0
      vc       = 0.0d0

      WRITE(LUPRI,*) 'Calculating the density and its derivatives.'
      DO i=1,size(points,2)
!     JHS Get the values of the basisfunctions in every gridpoint
        CALL GETSOS(GAO,GAO1,GAO2,GAB1,NCNT,  &
                    points(1,i),points(2,i),points(3,i), &
                    BUF,NDIMENSION,2,DOLND,print_fde)
!     JHS calculate the SCF density and its derivatives. 
        CALL GETRHO(rho(1,i),0,GAO,DMAT,BUF)
        CALL GETGRHO(grad_rho(1:3,i),0,GAO,GAO1,DMAT,BUF)
        CALL GETLRHO(hes_rho(1,i),hes_rho(2,i),hes_rho(3,i), &
             hes_rho(4,i),hes_rho(5,i),hes_rho(6,i), &
             GAO,GAO1,GAO2,DMAT,BUF)
!     JHS Calculate the coulomb potential.
      END DO

      call dealloc (gao)
      call dealloc (gao1)
      call dealloc (gao2)
      call dealloc (gab1)
      call dealloc (ncnt)

      WRITE(LUPRI,*) & 
          'Calculating the coulomb potential and its derivatives.'
      CALL EMB_VC(DMAT,points,vc,print_fde, WORK,LWORK)

      WRITE(LUPRI,*) 'Writing data to xml_file:', GRIDOUT_NAME

      write (csizept1,*)     size(points,1)    ! ; trim (csizept1) 
      write (csizept2,*)     size(points,2)    ! ; trim (csizept2) 
      write (csizehesrho1,*) size(hes_rho,1)   ! ; trim (csizehesrho1) 
      write (csizegrdrho1,*) size(grad_rho,1)  ! ; trim (csizegrdrho1) 

      tag=>xml_open(GRIDOUT_NAME)
        write(lupri,*) 'tst'
        tag=>xml_tag_open(tag,'grid')
          call xml_attribute(tag,'size',csizept2)
          tag=>xml_tag_open(tag,'dataset','gridpoints')
            call xml_attribute(tag,'size',  csizept2)
            call xml_attribute(tag,'width', csizept1)
            call xml_add_data(tag,points)
          tag=>xml_tag_close(tag)

          tag=>xml_tag_open(tag,'dataset','density')
            call xml_attribute(tag,'size',  csizept2)
            call xml_attribute(tag,'width', '1')
            call xml_add_data(tag,rho)
          tag=>xml_tag_close(tag)

          tag=>xml_tag_open(tag,'dataset','gradient')
            call xml_attribute(tag,'size', csizept2)
            call xml_attribute(tag,'width',csizegrdrho1)
            call xml_add_data(tag,grad_rho)
          tag=>xml_tag_close(tag)

          tag=>xml_tag_open(tag,'dataset','hessian')
            call xml_attribute(tag,'size',  csizept2)
            call xml_attribute(tag,'width', csizehesrho1)    
            call xml_add_data(tag,hes_rho)
          tag=>xml_tag_close(tag)

          tag=>xml_tag_open(tag,'dataset','vc')
            call xml_attribute(tag,'size',  csizept2)
            call xml_attribute(tag,'width', '1')
            call xml_add_data(tag,vc)
          tag=>xml_tag_close(tag)
        tag=>xml_tag_close(tag)

        write(*,*) 'WRITING:'       
        call xml_write(GRIDOUT_NAME)
        write(*,*) 'CLOSING:'
      tag=>xml_close(GRIDOUT_NAME)
      write(*,*) 'CLOSED!'

      call dealloc(rho)
      call dealloc(grad_rho)
      call dealloc(hes_rho)
      call dealloc(vc)

      CALL TIME()

      END SUBROUTINE
#endif

!
! auxiliary routines
!
!
! in case we use dalton 2.0 or dirac without enabling ERI
! we include printdmat.
!
#ifndef MOD_ERI
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     /* PRINTDMAT */
! aspg
! taken from eri/eri2ftuv.F, added here for convenience

      SUBROUTINE PRINTDMAT( NSIM, NORBT, D )
#include "implicit.h"
#include "priunit.h"
      INTEGER :: NSIM, NORBT
      REAL(KIND=8) :: D(:)

      DO ISIM=1,NSIM
        WRITE( LUPRI, * ) 'MATRIX NUMBER : ', ISIM
        DO I=1,NORBT
          WRITE( LUPRI, 400 ) (D((I-1)*NORBT+J),J=1,NORBT)
        END DO
      END DO

 400  FORMAT(6F16.12)

      RETURN
      END SUBROUTINE
!
#endif


      END MODULE FDE_MOD
