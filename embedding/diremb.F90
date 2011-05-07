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

! $Id$
!
! description : basic fde module, handling i/o of imported and exported
!               data (potentials, densities and density gradients etc) on
!               a grid, and also for setting options (do embedding response,
!               regenerate the embedding potential, set kinetic energy 
!               functional and so on).
!
!               the module also keeps some internal variables (like the grid
!               and the imported quantities) that can be made available to
!               different parts of the code.
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
      use xc_derv

#ifdef PRG_DIRAC
!miro: introduced in order to fix g95 integer*8 problem
#define SIZE(a,b) int(size(a,b))
#endif

      PRIVATE
      PUBLIC FDE_ImportStaticPotential_FileOpen, FDE_ImportStaticPotential_FileClose
      public FDE_setupFDEgrid, FDE_testFrozendensity, FDE_cleanupFDEgrid
      public FDEINI,FDEINP,FDE_GRIDOUT,have_frozen_density, update_vemb, import_vemb
      PUBLIC DO_GRID_OUT
      public fde_grid, fde_frozen_prp, fde_grid_np, fde_nprop
      public skip_kin_contribution, skip_xc_contribution
      public get_xml_export_level
      public get_fde_nadd_fun_derv

      LOGICAL,SAVE :: DO_GRID_OUT=.FALSE.
      logical,save :: have_frozen_density = .false.
      logical,save :: import_vemb = .false.
      logical,save :: update_vemb = .false.
      logical,save :: is_fdersp_alda = .false.
      logical,save :: skip_kin_contribution = .false.
      logical,save :: skip_xc_contribution  = .false.
      INTEGER,SAVE :: LUFDEPOT=45,LUTESTING=46,LUGRIDOUT=47, LUFDEFRZDNS=48
      INTEGER,SAVE :: print_fde=0, fde_grid_np=0, fde_nprop=1
      CHARACTER(LEN=60),SAVE :: FDEPOT='EMBPOT',FDEDNS='FRZDNS',GRIDOUT_NAME='GRIDOUT'
      CHARACTER(LEN=4),SAVE :: xmlexport_dens_level='HF'
      real(kind=8), pointer :: fde_grid(:,:), fde_frozen_prp(:,:) 

#ifdef MOD_UNRELEASED
!     public fde_kef,       fde_xcf
!     public fde_kef_alda,  fde_xcf_alda
!     public fde_kef_xalda, fde_xcf_xalda

      type(functional), save :: fde_kef,       fde_xcf
      type(functional), save :: fde_kef_alda,  fde_xcf_alda
      type(functional), save :: fde_kef_xalda, fde_xcf_xalda

! the default (orbital-free) kinetic energy functional for now is the thomas-fermi one
      character(len=80), save :: kefun_fde='kin_tf'
      character(len=80), save :: xcfun_fde='lda'

      logical, save :: fde_set_functionals(20) 
#endif

      real(8) :: hfx_out, mu_out, beta_out
      
      CONTAINS

        subroutine get_xml_export_level(level)
           CHARACTER(LEN=4), intent(OUT) :: level
   
           level = xmlexport_dens_level 
        end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDE_ImportStaticPotential_FileOpen: opens file containing the embedding potential
!                                    (static) exported on a grid
!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDE_ImportStaticPotential_FileClose: closes file containing the embedding potential
!                                      (static) exported on a grid
!
         subroutine FDE_ImportStaticPotential_FileClose()
            close(lufdepot) 
         end subroutine FDE_ImportStaticPotential_FileClose

#ifdef MOD_UNRELEASED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDE_ImportFrozenDensity_FileOpen: opens file containing the exported density and
!                                    other quantities (density gradient components etc)
!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDE_ImportFrozenDensity_FileClose: closes file containing the exported density and
!                                    other quantities (density gradient components etc)
!
         subroutine FDE_ImportFrozenDensity_FileClose()
            close(lufdefrzdns) 
         end subroutine FDE_ImportFrozenDensity_FileClose


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDE_setupFDEgrid: opens the file containing the frozen density, reads in the density
!                   and whatever other quantities we have exportend (such as the 
!                   components of the density gradient etc) using read_grid (internally,
!                   it is read_grid which will allocate storage, and later we have to 
!                   free it using FDE_cleanupFDEgrid for instance). for convenience, 
!                   we also set fde_grid_np (number of grid points) and fde_nprop 
!                   (number of properties read in) 
!
         subroutine FDE_setupFDEgrid()
            integer :: frozdens_unit

            call FDE_ImportFrozenDensity_FileOpen(frozdens_unit)
            call read_grid(frozdens_unit,fde_grid,fde_frozen_prp)
            call FDE_ImportFrozenDensity_FileClose()

            fde_grid_np = size(fde_grid,2)
            fde_nprop   = size(fde_frozen_prp,1)

         end subroutine FDE_setupFDEgrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDE_testFrozenDensity: integrates the imported density, useful to check if the
!                        number of electrons does match that of the frozen subsystems
!
         subroutine FDE_testFrozenDensity()
            integer      :: point = 0
            real(kind=8) :: frozen_particle_number = 0.0

            do point = 1, fde_grid_np 
               frozen_particle_number = frozen_particle_number + fde_grid(4,point)*fde_frozen_prp(3,point) 
            end do
            
            write (*,*) '  FDE info: imported frozen density integrates to ',frozen_particle_number,' electrons'
            write (*,*) '            grid points processed ',fde_grid_np
            write (*,*) '            number of properties  ',fde_nprop
            frozen_particle_number = 0.0

         end subroutine FDE_testFrozenDensity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDE_cleanupFDEgrid: dealocates the storage used for holding the grid, and whatever
!                     properties exported at the grid points
!
         subroutine FDE_cleanupFDEgrid()
            deallocate(fde_grid)
            deallocate(fde_frozen_prp) 
         end subroutine FDE_cleanupFDEgrid
#endif /* MOD_UNRELEASED */

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDEINI: performs the initialization in dirac
!
      SUBROUTINE FDEINI()
#include "implicit.h"
#include "priunit.h"
#include "dcbham.h"
         integer :: i

#ifdef MOD_UNRELEASED
#ifdef MOD_XCFUN

         update_vemb = .true. 
         have_frozen_density = .true.

         call parse_functional(kefun_fde, fde_kef,      hfx_out, mu_out, beta_out, .false.)
         call parse_functional(kefun_fde, fde_kef_alda, hfx_out, mu_out, beta_out, .false.)
         call parse_functional(xcfun_fde, fde_xcf,      hfx_out, mu_out, beta_out, .false.)
         call parse_functional(xcfun_fde, fde_xcf_alda, hfx_out, mu_out, beta_out, .false.)

         do i = 1, 20
            fde_set_functionals(i) = .true.
         enddo

#else 
         call quit('XCFun not available in this version')
#endif /* ifdef MOD_XCFUN */
#endif
         write (lupri,*) ' '
         CALL PRSYMB(LUPRI,'=',101,1)
         write (lupri,*) ' '
         write(LUPRI,'(1X,A)') &
            ' Frozen Density Embedding Calculation defaults:'
         call FDE_setupInfo(lupri)    
         CALL PRSYMB(LUPRI,'=',101,1)
         write (lupri,*) ' '
      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDEINP: dirac-specific input reader for the module
!
      SUBROUTINE FDEINP(WORD,Echo_Input)
#include "implicit.h"
#include "priunit.h" 
#include "dcbgen.h"
#include "dcbham.h"
      PARAMETER (NTABLE = 8)
!
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7, line*80
!
!
      logical :: Echo_Input,newdef
!
#include "memint.h"

      update_vemb = .false.
      have_frozen_density = .false.

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
              update_vemb = .false.
!             read filename for the embedding potential
              READ(LUCMD,'(A60)') FDEPOT 

            CASE ('.PRINT')
!             set printlevel 
              READ(LUCMD,*) print_fde 

#ifdef MOD_UNRELEASED
            CASE ('.SKIPK')
               skip_kin_contribution = .true.
            CASE ('.SKIPX')
               skip_xc_contribution  = .true.

            CASE ('.FRDENS')
!             read filename for the frozen density and density gradient
              READ(LUCMD,'(A60)') FDEDNS
              have_frozen_density = .true.

! read kinetic energy functional for non-additive term
            CASE ('.NAXCF')
#ifdef MOD_XCFUN
              read(lucmd,'(a80)') xcfun_fde
              call parse_functional(xcfun_fde, fde_xcf, hfx_out, mu_out, beta_out, .false.)
              write(*,*) 'testing ',xcfun_fde
#else         
              call quit('XCFun not available in this version')
#endif /* ifdef MOD_XCFUN */

! read kinetic energy functional for non-additive term
            CASE ('.NAKEF')
#ifdef MOD_XCFUN
              read(lucmd,'(a80)') kefun_fde
              call parse_functional(kefun_fde, fde_kef, hfx_out, mu_out, beta_out, .false.)
#else         
              call quit('XCFun not available in this version')
#endif /* ifdef MOD_XCFUN */

! read kinetic energy functional for non-additive term
            CASE ('.RSPLDA')
#ifdef MOD_XCFUN
              is_fdersp_alda = .true.
              call set_kin_fun_alda(fde_kef_alda)
              call set_xc_fun_alda(fde_xcf_alda)
#else
              call quit('XCFun not available in this version')
#endif /* ifdef MOD_XCFUN */

            CASE ('.LEVEL')
              READ(LUCMD,'(A4)') xmlexport_dens_level
              select case(xmlexport_dens_level)
              case('DHF','MP2')
                 WRITE(LUPRI,*) 'Outputted density will be taken from: ' &
                      ,trim(xmlexport_dens_level)
              case default
                 WRITE(LUPRI,*) 'Input for outputted density not a &
     &  x  valid calculation, HF will be used instead. Given:',xmlexport_dens_level
              end select
            CASE ('.GRIDOU')
! WRITE TO THIS GRIDFILE
              READ(LUCMD,'(A60)') GRIDOUT_NAME
              DO_GRID_OUT=.TRUE.
              DOPRP=.TRUE.
#endif /* MOD_UNRELEASED */
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
         if (Echo_Input) then
           write (lupri,*) ' '
           CALL PRSYMB(LUPRI,'=',101,1)
           write (lupri,*) ' '
           write(LUPRI,'(1X,A)') &
              ' Frozen Density Embedding Calculation, user-defined setup:'  
           call FDE_setupInfo(lupri) 
           CALL PRSYMB(LUPRI,'=',101,1)
           write (lupri,*) ' '
         endif 

      END SUBROUTINE


      subroutine FDE_setupInfo(lupri)
         integer :: lupri 

         CALL PRSYMB(LUPRI,'-',44,2)
         write(LUPRI,*)
         write(LUPRI,'(A,I5)')  ' * Print level                : ',print_fde
#ifdef MOD_UNRELEASED
         write(LUPRI,'(A,A20)') ' * Non-additive KE functional : ',kefun_fde
         write(LUPRI,'(A,A20)') ' * Non-additive XC functional : ',xcfun_fde
!        write(LUPRI,'(A,A20)') '   if .ALDA is set, will use : ',kefun_alda_fde
#endif
         if (import_vemb) then
            write(LUPRI,'(A,A60)')  ' * FDE Potential to be read from file: ',FDEPOT
            write(LUPRI,*)
            write(LUPRI,'(1X,2A,/1X,2A)') &
            'Frozen Density Embedding to include enviroment effects', &
            ' using the fixed potential method described in:',        &
            'A.S.P. Gomes, C. R. Jacob and L. Visscher, ',            &
            'Phys. Chem. Chem. Phys. 10 (2008) 5353â5362.'
#ifdef MOD_UNRELEASED
         else
            write(LUPRI,'(A,A60)')  ' * FDE Potential generated from frozen, active densities'
            update_vemb = .true.
            have_frozen_density = .true.
         endif

         if (have_frozen_density) then
            write(LUPRI,'(A,A60)') &
                ' * Density (and gradient) from frozen subsystems read from file: ',FDEDNS
         endif

         IF (DO_GRID_OUT) then
            write(LUPRI,'(2A)') &
            ' * FDE Gridfile with updated density written to:',GRIDOUT_NAME
#endif
         endif
      end subroutine FDE_setupInfo

#ifdef MOD_UNRELEASED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! time: keeps track of time spent in operations within FDE_GRIDOUT
!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDE_GRIDOUT: calculates and exports the coulomb potential on a grid (taken to be
!              the same from which the density etc was imported from). the export
!              uses jetze sikkema's xml library (sitting in gp/xml in the dirac
!              tree).
! 
       SUBROUTINE FDE_GRIDOUT(WORK,LWORK)
      use memory_allocator
      USE VC_GRID
      use xml_parser

      implicit none

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
#include "dcborb.h"
#include "mxcent.h"
#include "nuclei.h"
!     Passed variables
      REAL(KIND=8) :: WORK(*)
      INTEGER :: LWORK

      type(xml_tag),pointer :: tag

!     Local variables
      character(len=132) :: csizept1, csizept2
      character(len=132) :: csizegrdrho1, csizehesrho1
      real(kind=8), allocatable :: dmat(:)
      real(kind=8), allocatable :: cmo(:)
      real(kind=8) :: test_int_density
      REAL(kind=8), pointer :: points(:,:)

      real(8), allocatable :: grad_rho(:, :), hes_rho(:, :), vc(:), rho(:, :)

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

      call alloc(dmat, n2bbasxq)
      dmat = 0.0d0
      call alloc(cmo, ncmotq)
      cmo = 0.0d0

      select case(xmlexport_dens_level)
      case('MP2')
         write(lupri,*) 'Reading CC_DENSITY'
         call GET_CC_DENSITY(xmlexport_dens_level,'AO',DMAT,ierr)
         if (ierr.ne.0) then
            write(lupri,*) 'Error in reading density matrix'
            call quit('error readin density matrix')
         end if
      case default
!     HF density used
         write(lupri,*) 'Reading DFCOEF'
         CALL OPNFIL(LUCOEF,'DFCOEF','UNKNOWN','EMBDRV')
         call reacmo_no_work(lucoef,    &
                             'DFCOEF',  &
                             cmo,       &
                             (/0.0d0/), &
                             (/0/),     &
                             dummy,     &
                             2)
         if (print_fde.gt.20) then
            WRITE(LUPRI,*) 'CMO', CMO
            WRITE(LUPRI,*) 'DMAT OLD', DMAT
         end if
!     SCF construct density matrix
      CALL DENMAT(DMAT,CMO,print_fde)
      call dealloc(cmo)

!     TMP scaling
      end select
!     Correct DMAT for occupationnumbers, both for scf and mp2...
      dmat = 2.0d0*dmat

!     GET an id for LUGRIDOUT.
      LUGRIDOUT=OPEN_FILE(GRIDOUT_NAME)
      write (*,*) 'lugridout id is:',LUGRIDOUT
      IF (LUGRIDOUT.lt.0) &
           CALL QUIT('ERROR OPENING GRIDOUT FILE. FILE: ' &
           //GRIDOUT_NAME)
      IF (LUGRIDOUT.eq.0) &
           CALL QUIT('FILE FOR GRIDOUT DOES NOT EXISTS. FILE: ' &
           //GRIDOUT_NAME)

      WRITE(LUPRI,*) 'Reading grid.'
      nullify(points)
      CALL READ_GRID(LUGRIDOUT,points)
      CLOSE(LUGRIDOUT)
      WRITE(LUPRI,*) 'Grid read from file.'

      write(lupri,*) 'Number of gridpoints:', size(points,2)

      call alloc(rho,     1,SIZE(points,2))
      call alloc(grad_rho,3,SIZE(points,2))
      call alloc(hes_rho, 6,SIZE(points,2))
      call alloc(vc,        SIZE(points,2))

      rho      = 0.0d0
      grad_rho = 0.0d0
      hes_rho  = 0.0d0
      vc       = 0.0d0

      test_int_density = 0.0d0

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
         test_int_density = test_int_density + rho(1,i)*points(4,i)
      END DO

      WRITE(LUPRI,*) ' Integrated density adds up to ',test_int_density,&
     &               ' electrons.'

      call dealloc (gao)
      call dealloc (gao1)
      call dealloc (gao2)
      call dealloc (gab1)
      call dealloc (ncnt)

      WRITE(LUPRI,*) & 
          'Calculating the coulomb potential and its derivatives.'
      CALL EMB_VC(DMAT,points,vc,print_fde, WORK,LWORK)
      call dealloc(dmat)

      WRITE(LUPRI,*) 'Writing data to xml_file:', GRIDOUT_NAME

      write (csizept1,*)     size(points,1)    ! ; trim (csizept1) 
      write (csizept2,*)     size(points,2)    ! ; trim (csizept2) 
      write (csizehesrho1,*) size(hes_rho,1)   ! ; trim (csizehesrho1) 
      write (csizegrdrho1,*) size(grad_rho,1)  ! ; trim (csizegrdrho1) 

      tag=>xml_open(GRIDOUT_NAME)
        write(lupri,*) 'tst'
        tag=>xml_tag_open(tag,'grid')
          call set_attribute(tag,'size',csizept2)
          tag=>xml_tag_open(tag,'dataset','gridpoints')
            call set_attribute(tag,'size',  csizept2)
            call set_attribute(tag,'width', csizept1)
            call xml_add_data(tag,points)
          tag=>xml_tag_close(tag)

          tag=>xml_tag_open(tag,'dataset','density')
            call set_attribute(tag,'size',  csizept2)
            call set_attribute(tag,'width', '1')
            call xml_add_data(tag,rho)
          tag=>xml_tag_close(tag)

          tag=>xml_tag_open(tag,'dataset','gradient')
            call set_attribute(tag,'size', csizept2)
            call set_attribute(tag,'width',csizegrdrho1)
            call xml_add_data(tag,grad_rho)
          tag=>xml_tag_close(tag)

          tag=>xml_tag_open(tag,'dataset','hessian')
            call set_attribute(tag,'size',  csizept2)
            call set_attribute(tag,'width', csizehesrho1)    
            call xml_add_data(tag,hes_rho)
          tag=>xml_tag_close(tag)

          tag=>xml_tag_open(tag,'dataset','vc')
            call set_attribute(tag,'size',  csizept2)
            call set_attribute(tag,'width', '1')
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
! auxiliary routines, dirac/dalton-specific
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

#endif /* MOD_UNRELEASED */


#ifdef  MOD_XCFUN

   subroutine get_fde_nadd_fun_derv(order,        &
                                   block_length, &
                                   w,            &
                                   r_0,          &
                                   r_0_t,        &
                                   z_0,          &
                                   z_0_t,        &
                                   n,            &
                                   derv,         &
                                   parallel_xc)

   use dirac_cfg
   use xc_derv
   use xc_max_block_length

!   ----------------------------------------------------------------------------
    integer,          intent(in)  :: order
    integer,          intent(in)  :: block_length
    real(8),          intent(in)  :: w(*)
    real(8),          intent(in)  :: r_0(*)
    real(8),          intent(in)  :: z_0(*)
    real(8),          intent(in)  :: r_0_t(*)
    real(8),          intent(in)  :: z_0_t(*)
    integer,          intent(in)  :: n
    real(8),          intent(out) :: derv(max_block_length, 0:n)
    logical,          intent(in)  :: parallel_xc
    real(8)                       :: dervtmp(max_block_length, 0:n)
    integer :: i, j
!   ----------------------------------------------------------------------------

    derv = 0.0d0

    if (fde_set_functionals(order)) then
       call xcfun_set_functional(fde_xcf     , order, parallel_xc) 
       call xcfun_set_functional(fde_xcf_alda, order, parallel_xc) 
       call xcfun_set_functional(fde_kef,      order, parallel_xc)
       call xcfun_set_functional(fde_kef_alda, order, parallel_xc)
       fde_set_functionals(order) = .false.
   endif

    if (.not.skip_xc_contribution) then

       call get_functional_derv(fde_xcf,       &
                                fde_xcf_alda,  &
                                fde_xcf_alda,  &
                                order,         &
                                block_length,  &
                                w,             &
                                r_0_t,         &
                                z_0_t,         &
                                n,             &
                                dervtmp,       &
                                alda_real=is_fdersp_alda, &
                                alda_imag=is_fdersp_alda)

       do i = 1, max_block_length
          do j = 0, n
             derv(i,j) = derv(i,j) + dervtmp(i,j)
          enddo
       enddo

       call get_functional_derv(fde_xcf,       &
                                fde_xcf_alda,  &
                                fde_xcf_alda,  &
                                order,         &
                                block_length,  &
                                w,             &
                                r_0,           &
                                z_0,           &
                                n,             &
                                dervtmp,       &
                                alda_real=is_fdersp_alda, &
                                alda_imag=is_fdersp_alda)

       do i = 1, max_block_length
          do j = 0, n
             derv(i,j) = derv(i,j) - dervtmp(i,j)
          enddo
       enddo

    endif

! now we construct the contribution from the kinetic energy to the response
! of the embedding potential.
! first set functional to the of-kef, do the non-additive ke contribution

    if (.not.skip_kin_contribution) then

       call get_functional_derv(fde_kef,       &
                                fde_kef_alda,  &
                                fde_kef_alda,  &
                                order,         &
                                block_length,  &
                                w,             &
                                r_0_t,         &
                                z_0_t,         &
                                n,             &
                                dervtmp,       &
                                alda_real=is_fdersp_alda, &
                                alda_imag=is_fdersp_alda)

       do i = 1, max_block_length
          do j = 0, n
             derv(i,j) = derv(i,j) + dervtmp(i,j)
          enddo
       enddo

! now subtract the active subsystem's kinetic energy contribution

       call get_functional_derv(fde_kef,       &
                                fde_kef_alda,  &
                                fde_kef_alda,  &
                                order,         &
                                block_length,  &
                                w,             &
                                r_0,           &
                                z_0,           &
                                n,             &
                                dervtmp,       &
                                alda_real=is_fdersp_alda, &
                                alda_imag=is_fdersp_alda)

       do i = 1, max_block_length
          do j = 0, n
             derv(i,j) = derv(i,j) - dervtmp(i,j)
          enddo
       enddo

    endif

! ! and restore the default dft functional, if we are to use it
! ! down the line...
!     call xcfun_set_functional(f, order, parallel_xc)
!     call xcfun_set_functional(f_alda, order, parallel_xc)
!     call xcfun_set_functional(f_xalda, order, parallel_xc)

  end subroutine get_fde_nadd_fun_derv

#endif /* ifdef MOD_XCFUN */




      END MODULE FDE_MOD
