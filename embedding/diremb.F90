! $Id: diremb.F 9568 2009-09-15 09:50:44Z gomes $
!
! subroutine  : GetFDEPotential_FromFile
!
! authors     : Andre Gomes, Christoph Jacob (VUA)
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
! arguments   :
!
!   1. fdemat : matrix to hold the embedding potential. can be the
!               one-electron fock matrix, or any other matrix one
!               wishes to use... 
!   2.  lwork : 
!   3.   work : workspace array 
!
! in the end  : the embedding potential is added to fdemat. 
!
      MODULE FDE_MOD
      USE GRID_FILE
#ifdef MOD_XML
      USE FILE, ONLY: OPEN_FILE
#endif

      PRIVATE
!     PUBLIC GetFDEPotential_FromFile, FDEINI,FDEINP,FDE_GRIDOUT
      PUBLIC FDE_ImportStaticPotential_FileOpen, FDE_ImportStaticPotential_FileClose
      public FDEINI,FDEINP,FDE_GRIDOUT
      PUBLIC DO_GRID_OUT
           
      LOGICAL,SAVE :: DO_GRID_OUT=.FALSE.
      INTEGER,SAVE :: LUFDEPOT=45,LUTESTING=46,LUGRIDOUT=47
      INTEGER,SAVE :: print_fde=0
      CHARACTER(LEN=60),SAVE :: FDEPOT='EMBPOT',GRIDOUT_NAME=''
      CHARACTER(LEN=60),SAVE :: out_level='HF'
      
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


!
! deprecated code!, being replaced by the code based on the "numoper_integrals" module,
! that writes the embedding potential matrix to AOPROPER...
!
      SUBROUTINE GetFDEPotential_FromFile(fdemat,work,lwork)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"

#ifndef PRG_DIRAC
!radovan: it needs dfthri, but not if DIRAC is used
#include "dftcom.h"
#endif

#include "ccom.h"
#include "inforb.h"
      INTEGER :: p,q
      real(kind=8), pointer :: points(:,:),vc(:)
      REAL(KIND=8) :: work(lwork)
      REAL(KIND=8),ALLOCATABLE :: GAO(:),BUF(:),NCNT(:)

#if defined (PRG_DIRAC)
      REAL(KIND=8) :: fdemat(nbast,nbast)
#else
!     maybe better to make implicit size for dalton, or in general?
!     dimension fdemat(*)
#include "symmet.h"
      REAL(KIND=8) :: fdemat(nnbasx)
#endif
      LOGICAL :: debug, embpot_found
!
      call qenter('GetFDEPotential_FromFile')
#include "memint.h"

      ndimension   = nbast 
      debug        = .false.

      if (debug) write (LUPRI,*) 'debug: before opening FDE pot file'

!     before proceeding, we check if the FDE  potential is available, 
      inquire(FILE=FDEPOT,EXIST=embpot_found)
      if (.not.embpot_found) then
         write(LUPRI,*) FDEPOT, ' file not found, we can not proceed..'
         call quit ('Fatal error: FDE potential file not present')
      endif
 
      if (debug) then
!        for testing purposes, we will use this file to dump the grid we read.
         open(lutesting,              &
              FORM   = 'FORMATTED',   &
              STATUS = 'UNKNOWN',     &
              ACCESS = 'SEQUENTIAL',  &
              FILE   =  'TSTPOT')
      endif

      open(lufdepot,                  &
             FORM   = 'FORMATTED',    &
             STATUS = 'UNKNOWN',      &
             ACCESS = 'SEQUENTIAL',   &
             FILE   =  FDEPOT)

!     we start processing the grid file by first getting the number of point in 
!     the grid. the size of the grid is then used to allocate workspace for 
!     reading the actual grid, along with the embedding potential at each point
      rewind lufdepot

      CALL READ_GRID(lufdepot,points,vc)
    
!      if (debug) then
!         call header ('1-electron matrix, before fde update',-1)
!         call PRINTDMAT(1,ndimension,fdemat)
!      endif

!     we allocate the variables to hold the value of the basis functions
!     at the grid ponts, and also some work arrays
#if defined (PRG_DIRAC)
      ALLOCATE(GAO(NDIMENSION))
      ALLOCATE(NCNT(4*NDIMENSION))
      ALLOCATE(BUF(4*NDIMENSION))
#else
      ALLOCATE(GAO(NDIMENSION*NTYPESO))
      ALLOCATE(KNCNT(NDIMENSION))
      call setupsos(0,.false.,idum1,idum2)
#endif

!     Initialize variables.
      gao=0
      NCNT=0
      BUF=0
      nder   = 0
      irepop = 0

      if (debug) write (lutesting,'(I10)') size(points)

      do i = 1, size(points,2)       ! begin numerical integration

         if (debug) write (lutesting,*) points(:,i)


#if defined (PRG_DIRAC)
!       we get the value of all gaussians at the current grid point with 
!       routines from dft
         call getsos(gao,dummy,dummy,dummy,    &   
                     ncnt, &
                     points(1,i), &
                     points(2,i), &
                     points(3,i),BUF, &
                     ndimension,0,.false.,print_fde)

!        only the upper triangular part of the one electron embedding 
!        potential matrix is calculated. 
         do p = 1, ndimension
           do q = p, ndimension
              
!            first we get the value of the orbital pair (gamma_pq) at point i
!            and multiply it by the embedding potential, and just for the
!            LL or SS blocks... 
             if ((p.le.nlarge).and.(q.gt.nlarge)) then
               gamma_pq      = 0
             else
               gamma_pq      = gao(p)*gao(q)
             endif
!            ep_pq_update  = points(4,i)*vc(i)*gamma_pq
             ep_pq_update  = vc(i)*gamma_pq

!            in dirac, we can work with the full matrix
             fdemat(p,q) = fdemat(p,q) + ep_pq_update 
             fdemat(q,p) = fdemat(p,q)

          enddo ! loop over q
        enddo ! loop over p
!
#else   /* DALTON */
!       getsos_fde is a simplification of getsos for our case. if there was
!       no symmetry, getaos could have been used directly... 
         call getsos_fde(gao,ncnt,
     &               points(1,i),
     &               points(2,i),
     &               points(3,i),
     &               work,lwork,
     &               ndimension,dfthri,iprdfe)

!         in order to treat the symmetry blocking with packed storage, 
!         i (aspg) came up with the following (remember that the matrix
!         is stored as an array with the upper/lower triangular of the
!         first irrep coming before the second and so on):
!         
!         a. "irrep offset" (irrep_offt) will tell us where in the array a 
!             given irrep begins 
!         b. "packed storage offset" (jpkds_offt)) will tell us the number 
!            of elements from the packed storage within a given irrep that
!            we've looped over so far 
!         c. "basis offset" (ibas_offt) will tell us the starting index of 
!            the basis set in the gao array for the different irreps 
          irrep_offt = 0 
          jpkds_offt = 0
          ibas_offt  = 0

          do lambda = 1, maxrep+1
             do p = 1,naos(lambda)                 
                do q = p,naos(lambda)    ! loop over the unique elements

                   gamma_pq      = gao(ibas_offt+p)
     &                           * work(ibas_offt+q)
                   ep_pq_update  = points(4,i)*vc(i)*gamma_pq

!              we do the indexing for packed array storage [ij=i+j(j-1)/2] (see
!              e.g the lapack documentation), but taking care of the offsetting
!               due to the symmetry blocking
                  pq = irrep_offt + p + q*(q-1)/2
                  jpkds_offt =  jpkds_offt + 1

                  fdemat(pq)  = fdemat(pq) + ep_pq_update 
               enddo ! loop over q
            enddo ! loop over p

            irrep_offt = irrep_offt + jpkds_offt 
            ibas_offt  = ibas_offt  + naos(lambda)
            jpkds_offt = 0

         enddo ! loop over symmetry representations
#endif /* DIRAC vs DALTON, matrix update over the grid */
      enddo ! loop over points

!     if (debug) then
!        call header ('1-electron matrix, after FDE update',-1)
!        call printdmat(1,ndimension,fdemat)
!     endif

#if defined (PRG_DIRAC)
! wonder why this crashes with dalton...
#endif

      DEALLOCATE(NCNT)
      DEALLOCATE(BUF)
      DEALLOCATE(GAO)
      
      DEALLOCATE(vc)
      DEALLOCATE(points)
      
      if (debug) close(lutesting)
      call qexit('GetFDEPotential_FromFile')

      END SUBROUTINE

      SUBROUTINE FDEINI()

      END SUBROUTINE


      SUBROUTINE FDEINP(WORD,Echo_Input)
#include "implicit.h"
#include "priunit.h" 
#include "dcbgen.h"
      PARAMETER (NTABLE = 8)
!
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
!
!
      logical :: Echo_Input,newdef
!
#include "memint.h"

#if defined (PRG_DIRAC)
      iunit = LUPRI
#else
      iunit = LUPRI
#endif

      print_fde=0
      ICHANG = 0
      IF (WORD(1:7)=='*FDE   ') THEN
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
! read filename for the embedding potential
              READ(LUCMD,'(A30)') FDEPOT 
            CASE ('.PRINT')
! set printlevel 
              READ(LUCMD,*) print_fde 
            CASE ('.LEVEL')
              READ(LUCMD,'(A30)') out_level
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
              READ(LUCMD,'(A30)') GRIDOUT_NAME
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
         write(LUPRI,'(A,A60)')  ' * FDE Potential read from file: ', &
         FDEPOT
         IF (DO_GRID_OUT) write(LUPRI,'(2A)') &
         ' * FDE Gridfile with updated density written to:',GRIDOUT_NAME
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
      REAL(kind=8), allocatable, target :: points(:,:),grad_rho(:,:)
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
        CALL GETGRHO(grad_rho(3,i),0,GAO,GAO1,DMAT,BUF)
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
!      FLUSH(LUPRI)
      CALL EMB_VC(DMAT,points,vc,print_fde, WORK,LWORK)

      WRITE(LUPRI,*) 'Writing data to xml_file:', GRIDOUT_NAME

      write (csizept1,*)     size(points,1)    ! ; trim (csizept1) 
      write (csizept2,*)     size(points,2)    ! ; trim (csizept2) 
      write (csizehesrho1,*) size(hes_rho,1)   ! ; trim (csizehesrho1) 
      write (csizegrdrho1,*) size(grad_rho,1)  ! ; trim (csizegrdrho1) 

#ifdef MOO
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
            call xml_attribute(tag,'width', csizept1)
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
#endif
      

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


#if !defined (PRG_DIRAC)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETSOS_FDE(GSO,NCNT,PX,PY,PZ,WORK,LWORK,NBAST,
     &                      DFTHRI,IPRINT)
!
!     original implementation of getsos: T. Helgaker feb 01
!     adaptation to fde case: A. Gomes, apr 07
!
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h> 
!        
      LOGICAL DOLND, DOGGA
      DIMENSION GSO(NBAST*NTYPSO), WORK(LWORK), NCNT(NBAST)
!       
#include <dftinf.h>
#include <symmet.h> 
!
#include <ibtfun.h> 
!   
      DOLND = .FALSE.
      DOGGA = .FALSE.

      IF (MAXREP.EQ.0) THEN
         CALL DFTAOS(GSO,dummy,dummy,dummy,dummy,
     &               NCNT,PX,PY,PZ,WORK,LWORK,NBAST,DOLND,
     &               DOGGA,DFTHRI,IPRINT)
      ELSE
         KGAO = 1
         KLST = KGAO + NTYPSO*NBAST 
         IF (KLST.GT.LWORK) CALL STOPIT('DFTAOS','LWORK',KLST,LWORK)
         CALL DFTAOS(WORK(KSO0),dummy,dummy,dummy,dummy,NCNT,PX,PY,PZ,
     &               WORK,LWORK,NBAST,DOLND,DOGGA,DFTHRI,IPRINT)
         CALL DFTSOS(WORK(KSO0),GSO,NBAST,NTYPSO,NCNT,IPRINT)
      END IF
      END
#endif

      END MODULE FDE_MOD
