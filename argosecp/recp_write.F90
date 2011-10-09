MODULE RECP_WRITE

CONTAINS 

! ------------------------------------------------

SUBROUTINE RECP_WRITE_ISA(ircru,jrcru,nst,mstu,nprir,ijsf,nopir,ilxyz,iapt, &
                          mpru,ipr,N_INTEGRAL,h2,esf)

! # put integrals and labels into the intermediate storage arrays.
  USE RECP_IPT
  IMPLICIT NONE
! global variables
  INTEGER :: ircru,jrcru,nst,mstu,nprir(2,mstu,*),ijsf,nopir(*),ilxyz(3,*),iapt(3)
  INTEGER :: mpru,ipr
  INTEGER :: N_INTEGRAL(3)     ! number of integral element of each type 
  REAL(8) :: h2(*)
  LOGICAL :: esf
! local variables
  INTEGER :: ist,iop,INDEX1,NPAIR,IPAIR,I,J,INT_TYPE 

  INDEX1 = 0

  do I = 1, ircru
    if(esf) jrcru = I
    do J = 1, jrcru
      do ist = 1, nst
!       setting NPAIR ?
        if ( esf .and. J.eq.I ) then
          NPAIR = nprir(1,ist,ijsf)
        else
          NPAIR = nprir(2,ist,ijsf)
        endif

        if ( (nopir(ist).ne.0) .and. (NPAIR .ne. 0) ) then
          do iop = 1, nopir(ist)
!           INT_TYPE : integral type: S,T,V = 1 / SO(x,y,z) = 1,2,3
            INT_TYPE = ilxyz(iop,ist)

            do IPAIR = 1, NPAIR
              INDEX1 = INDEX1 + 1
              iapt(INT_TYPE) = iapt(INT_TYPE) + 1
              IF (INT_TYPE.EQ.1) WRITE(64,'(x,f20.16)') h2(INDEX1)
              IF (INT_TYPE.EQ.2) WRITE(65,'(x,f20.16)') h2(INDEX1)
              IF (INT_TYPE.EQ.3) WRITE(66,'(x,f20.16)') h2(INDEX1)
              IF (INT_TYPE.EQ.1) WRITE(67,'(x,i10)') (ipr+IPAIR) 
              IF (INT_TYPE.EQ.2) WRITE(68,'(x,i10)') (ipr+IPAIR)
              IF (INT_TYPE.EQ.3) WRITE(69,'(x,i10)') (ipr+IPAIR)

              N_INTEGRAL(INT_TYPE) = N_INTEGRAL(INT_TYPE) + 1
            enddo
          enddo
        endif

        ipr = ipr + NPAIR

      enddo
    enddo
  enddo

END SUBROUTINE RECP_WRITE_ISA

! ------------------------------------------------

SUBROUTINE RECP_WRITE_MAIN(nop,iapt,ipr,ipr1,ipr2,N1MAX,L1REC,n_bfn_sym, &
                           N_INTEGRAL,MN_INTEGRAL,INT_AREP,nnam,mpru,nbft,nnbft,nst)
  USE RECP_IPT
  USE RECP_OUTPUT
! write 1-e integral arrays, Y.C. Park
  IMPLICIT NONE
! variables
#include "argos_mxvalue.h"
#include "argos_files.h"
#include "argos_error.h"
#include "inc_print.h"

! global variable
  INTEGER :: mpru, nbft, nnbft, nst   ! common block (parmi)
  INTEGER :: nnam                     ! common block (ntgr)
  INTEGER :: INT_AREP
  INTEGER :: N1MAX,  L1REC
  INTEGER :: N_INTEGRAL(3) ! number of integral element of each type
  INTEGER :: MN_INTEGRAL   ! maximum of N_INTEGRAL(i,i=1-3)
  INTEGER :: nop, iapt(3), ipr, ipr1, ipr2

! local variable
  REAL(8), ALLOCATABLE :: array_stvc(:)
  REAL(8), ALLOCATABLE :: array_so(:,:)
  REAL(8), ALLOCATABLE :: values_stvc(:)    ! to be fixed
  INTEGER, ALLOCATABLE :: labels_stvc(:,:)  ! from labels(nipv,*)

  INTEGER :: int_type     ! integral type 
  INTEGER :: ilfact = 1024
  INTEGER :: itypea, itypeb
  INTEGER :: lblop
  INTEGER :: n_bfn_sym(mstup)  ! (nbpsy or nso)
  INTEGER :: num, I,value_temp,value_index
  CHARACTER(8) :: chrtyp

  INTEGER :: itypan(7), itypbn(7)
! #             s    t    v  veff so:x so:y so:z
  data itypan / 0,   0,   0,   0,   2,   2,   2 /
  data itypbn / 0,   1,   2,   3,   0,   1,   2 /

  CALL OUTPUT_LOC('RECP_WRITE_MAIN','E')

  IF (nnam.EQ.5) THEN
     OPEN (UNIT=74,FILE='RECP_ISA1')
     OPEN (UNIT=75,FILE='RECP_ISA2')
     OPEN (UNIT=76,FILE='RECP_ISA3')
     OPEN (UNIT=77,FILE='RECP_ISB1')
     OPEN (UNIT=78,FILE='RECP_ISB2')
     OPEN (UNIT=79,FILE='RECP_ISB3')
     REWIND(74)
     REWIND(75)
     REWIND(76)
     REWIND(77)
     REWIND(78)
     REWIND(79)
  ELSE
     OPEN (UNIT=74,FILE='RECP_ISA1')
     OPEN (UNIT=77,FILE='RECP_ISB1') 
     REWIND(74)
     REWIND(77)
  ENDIF

  ALLOCATE( array_stvc(MN_INTEGRAL) )
  ALLOCATE( array_so(nbft,nbft) )
  ALLOCATE( values_stvc(MN_INTEGRAL) )
  ALLOCATE( labels_stvc(2,MN_INTEGRAL) )

  IF (PRINT_LEVEL.GE.5) PRINT *,' nop :',nop
! # write out the components of the operator arrays.
! value_index = 1
  do lblop = 1, nop

!    check that array(*) was correctly allocated.
!    if (iapt(lblop).gt.mpru)   call bummer('oneint: array(*) overflow, iapt=',iapt(lblop), faterr )

!    # check integral type
     itypea = itypan(nnam + lblop - 1)
     itypeb = itypbn(nnam + lblop - 1)
     call siftyp( itypea, itypeb, chrtyp )
     call RECP_WRITE_TYPE(itypea, itypeb, int_type, lblop,PRINT_LEVEL)

    DO I = 1, N_INTEGRAL(lblop) 
       IF (lblop.EQ.1) READ(74,'(x,f20.16)') values_stvc(I)
       IF (lblop.EQ.2) READ(75,'(x,f20.16)') values_stvc(I)
       IF (lblop.EQ.3) READ(76,'(x,f20.16)') values_stvc(I)
       IF (lblop.EQ.1) READ(77,'(x,i10)') value_temp 
       IF (lblop.EQ.2) READ(78,'(x,i10)') value_temp 
       IF (lblop.EQ.3) READ(79,'(x,i10)') value_temp
!       print *,'values_stvc(',I,')', values_stvc(I), value_temp
       labels_stvc(1,I) =     IPT_IL(value_temp)/ilfact
       labels_stvc(2,I) = mod(IPT_IL(value_temp),ilfact)
!      print *,'labels',lblop,labels_stvc(1,I),labels_stvc(2,I),I
!      value_index = value_index + 1
    ENDDO
!   # write integral to file
    call RECP_WRITE_ARRAY( nst,int_type,n_bfn_sym,nbft,labels_stvc,mpru,array_stvc,array_so, &
                            values_stvc,N_INTEGRAL,MN_INTEGRAL,lblop)
    call RECP_WRITE_FILE( int_type, N_INTEGRAL, nbft,array_stvc, array_so,INT_AREP,MN_INTEGRAL, &
                           LBLOP)

    IF(PRINT_LEVEL.GE. 5) write (nlist,'(4x,a,a)') chrtyp,' integrals were written in a file'
    IF(PRINT_LEVEL.GE.10) PRINT *,'   N_INTEGRAL(',LBLOP,')',N_INTEGRAL(LBLOP)

  enddo

  DEALLOCATE( array_stvc )
  DEALLOCATE( array_so )
  DEALLOCATE( values_stvc )
  DEALLOCATE( labels_stvc )

  IF (nnam.EQ.5) THEN
     CLOSE (UNIT=74,STATUS='DELETE')
     CLOSE (UNIT=75,STATUS='DELETE')
     CLOSE (UNIT=76,STATUS='DELETE')
     CLOSE (UNIT=77,STATUS='DELETE')
     CLOSE (UNIT=78,STATUS='DELETE')
     CLOSE (UNIT=79,STATUS='DELETE')
  ELSE
     CLOSE (UNIT=74,STATUS='DELETE')
     CLOSE (UNIT=77,STATUS='DELETE')
  ENDIF
  CALL OUTPUT_LOC('RECP_WRITE_MAIN','X')
END SUBROUTINE RECP_WRITE_MAIN

! ------------------------------------------------


SUBROUTINE RECP_WRITE_ARRAY(n_sym,int_type,n_bfn_sym,nbft,labels_stvc,mpru,array_stvc,array_so, & 
                             values_stvc,N_INTEGRAL,MN_INTEGRAL,lblop)
! From cnvrt_stvz, Y.C. Park
  USE RECP_OUTPUT
  IMPLICIT NONE
#include "argos_mxvalue.h"
#include "inc_print.h"
  integer :: n_bfn_sym(mstup) ! nbpsy or nso : number of basis functions per symmetry block
  integer :: int_type         ! integral type 
  integer :: n_sym            ! nsym or nst : number of symmetry blocks
  integer :: nbft
  integer :: mpru
  integer :: lblop
  integer :: N_INTEGRAL(3), MN_INTEGRAL
  real(8) :: array_stvc(:)   ! array_stvc(MN_INTEGRAL) 
  real(8) :: array_so(:,:)     ! array_so(nbft,nbft) 
  real(8) :: values_stvc(:)  ! values_stvc(MN_INTEGRAL)
  integer :: labels_stvc(2,*) ! labels(1:2,1:n1max) orbital label buffer. 
! variables : local
  integer :: i, j, ij
  integer :: n_tot        ! n_total
  integer :: n_n_tot      ! n! of n_total
  integer :: n_skp(8)     ! n_skip
  integer :: n_n_skp(8)   ! n! of n_skip

  integer :: i_lab, j_lab
  integer :: SYM1, SYM2, SYM_12   ! isym, jsym, ijsym

  integer :: nbfmx=50000     ! FIXME 
  integer, ALLOCATABLE :: sym_bfn(:) ! symb(1:nbft) symmetry index of each basis function 'n_bfn_sym' is needed

  integer :: nndxf
  nndxf(i) = ( i * ( i - 1 ) ) / 2     ! define nndxf

  CALL OUTPUT_LOC('RECP_WRITE_ARRAY','E')
  ALLOCATE ( sym_bfn(nbfmx) )
  DO I = 1, nbfmx
     sym_bfn(i) = 0
  ENDDO

  IF (PRINT_LEVEL .GE. 5) THEN 
     PRINT *,'   INT_TYPE       =',INT_TYPE
     PRINT *,'   N_INTEGRAL     =',N_INTEGRAL(LBLOP)
     PRINT *,'   N_SYM          =',N_SYM
     PRINT *,'   MN_INTEGRAL    =',MN_INTEGRAL
  ENDIF
! make sym_bfn 
  n_tot  = 0
  n_n_tot = 0
  do SYM1 = 1, n_sym
     n_skp(SYM1)   = n_tot  
     n_n_skp(SYM1) = n_n_tot 
     n_tot          = n_tot   + n_bfn_sym(SYM1)               
     n_n_tot        = n_n_tot + nndxf(n_bfn_sym(SYM1) + 1)  
!    print *, '1-n_bfn_sym :',SYM1,n_bfn_sym(SYM1)
     do j = (n_skp(SYM1)+1), n_tot
        sym_bfn(j)  =SYM1
     enddo
  enddo

  if (int_type .le. 4) then
! # symmetric matrix, totally symmetric operator.
!    #  initialize matrix element 
     array_stvc = 0

     do i = 1, N_INTEGRAL(1)
       i_lab  = max( labels_stvc(1,i), labels_stvc(2,i) )
       j_lab  = min( labels_stvc(1,i), labels_stvc(2,i) )
       SYM1  = sym_bfn(i_lab)
       SYM2  = sym_bfn(j_lab)
       SYM_12 = nndxf(SYM1) + SYM2
       if ( SYM1 .eq. SYM2 ) then
!         # only valid array(*) elements are referenced.
          ij  =  n_n_skp(SYM1) + nndxf(i_lab - n_skp(SYM1)) + j_lab - n_skp(SYM1)
!       | ij  :  skip for       +                             +       -              |
!       |        each symmetry                                                       | 
!         print *,'ij',ij,i_lab,j_lab,labels_stvc(1,i), labels_stvc(2,i)
          array_stvc(ij) = values_stvc(i)
       endif
     enddo

  elseif ( int_type .ge. 5 ) then
! # antisymmetric matrix, nontotally symmetric operator.
! # initialize spin-orbit matrix element 
     array_so = 0
     do i = 1, N_INTEGRAL(LBLOP)
        i_lab  = labels_stvc(1,i)
        j_lab  = labels_stvc(2,i)
        if (i_lab .lt. j_lab) then
           i_lab = labels_stvc(2,i)
           j_lab = labels_stvc(1,i)
!          the sign change upon transposition.
           values_stvc(i) = -values_stvc(i)
        endif
        array_so(i_lab,j_lab) =  values_stvc(i)
        array_so(j_lab,i_lab) = -values_stvc(i)
!       print *,'array_so',i_lab,j_lab,array_so(i_lab,j_lab)
     enddo
  endif

  DEALLOCATE ( sym_bfn )

  CALL OUTPUT_LOC('RECP_WRITE_ARRAY','X')
END SUBROUTINE RECP_WRITE_ARRAY


! ------------------------------------------------


SUBROUTINE RECP_WRITE_FILE(int_type,N_INTEGRAL,nbft,array_stvc,array_so, &
                            INT_AREP,MN_INTEGRAL,LBLOP)
! write ro file
  USE RECP_OUTPUT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER :: LBLOP
  integer :: i,j
  integer :: int_type, INT_AREP
  integer :: N_INTEGRAL(3), MN_INTEGRAL
  integer :: nbft
  real(8) :: array_stvc(:)  ! array_stvc(MN_INTEGRAL)
  real(8) :: array_so(:,:)    ! array_so(nbft,nbft)

  LOGICAL :: FILE_EXIST  = .FALSE.
  LOGICAL :: ARGOS_DEBUG = .FALSE.
  INTEGER :: FILE_NUM(7)
  CHARACTER(10) :: FILE_INT(7),FILE_DBG(7),FILE_NAME(7)
  DATA FILE_NUM / 64, 64, 64, 64, 65, 66, 67 /
  DATA FILE_INT /'RECP_INT_S','RECP_INT_T','RECP_INT_V', &
                 'RECP_INT_C','RECP_INT_X','RECP_INT_Y', &
                 'RECP_INT_Z'/
  DATA FILE_DBG /'RECP_DBG_S','RECP_DBG_T','RECP_DBG_V', &
                 'RECP_DBG_C','RECP_DBG_X','RECP_DBG_Y', &
                 'RECP_DBG_Z'/
  DATA FILE_NAME /'Overlap   ','Kinetic   ','Potential ', &
                  'RECP AVE  ','RECP SO(X)','RECP SO(Y)', &
                  'RECP SO(Z)'/

  CALL OUTPUT_LOC('RECP_WRITE_FILE','E')

! # Skip writing S,T,V integral
  IF(INT_TYPE .LE. 3) GOTO 100

! # Delete if the file exist (for geometry optimization)
  INQUIRE(FILE=FILE_INT(INT_TYPE),EXIST=FILE_EXIST)
  IF (FILE_EXIST) THEN
     WRITE(6,*) 'File ',FILE_INT(INT_TYPE),' delete..'
     OPEN (UNIT=FILE_NUM(INT_TYPE),FILE=FILE_INT(INT_TYPE),FORM='UNFORMATTED',STATUS='OLD')
     CLOSE(UNIT=FILE_NUM(INT_TYPE),STATUS='DELETE')
  ENDIF

! # Write RECP integrals to file 
  OPEN(UNIT=FILE_NUM(INT_TYPE),FILE=FILE_INT(INT_TYPE),FORM='UNFORMATTED')
  IF(INT_TYPE.EQ.4)  WRITE(FILE_NUM(INT_TYPE)),(ARRAY_STVC(I),I=1,N_INTEGRAL(LBLOP)) !yc-check
  IF(INT_TYPE.GE.5)  WRITE(FILE_NUM(INT_TYPE)), ((ARRAY_SO(J,I),J=1,NBFT),I=1,NBFT)
  CLOSE(UNIT=FILE_NUM(INT_TYPE))

! # Write SOREP as 0 when AREP calculation
  IF(INT_TYPE.EQ.4 .AND. INT_AREP.eq.1) THEN 
     ARRAY_SO = 0
     DO INT_TYPE = 5, 7  !INT_TYPE = 5~7
!       Delete if the file exist (for geometry optimization)
        INQUIRE(FILE=FILE_INT(INT_TYPE),EXIST=FILE_EXIST)
        IF(FILE_EXIST) THEN
           OPEN (UNIT=FILE_NUM(INT_TYPE),FILE=FILE_INT(INT_TYPE), FORM='UNFORMATTED',STATUS='OLD')
           CLOSE(UNIT=FILE_NUM(INT_TYPE),STATUS='DELETE')
        ENDIF
!       Make 'ARRAY_SO = 0'
        OPEN(UNIT=FILE_NUM(INT_TYPE),FILE=FILE_INT(INT_TYPE), FORM='UNFORMATTED')
        WRITE(FILE_NUM(INT_TYPE)) ,((ARRAY_SO(J,I),J=1,NBFT),I=1,NBFT) 
        CLOSE(UNIT=FILE_NUM(INT_TYPE))
     ENDDO
     INT_TYPE = 4
  ENDIF

! # Write RECP integrals to file (for debug)
  IF (RECP_DBG.GE.5) THEN
     OPEN (UNIT=FILE_NUM(INT_TYPE),FILE=FILE_DBG(INT_TYPE), FORM='FORMATTED')
     IF (INT_TYPE.EQ.4) THEN
        WRITE(FILE_NUM(INT_TYPE),'(a)') FILE_NAME(INT_TYPE) 
        DO I = 1,N_INTEGRAL(LBLOP)
           WRITE(FILE_NUM(INT_TYPE),'(i7,2x,f20.16)') I,ARRAY_STVC(I)
        ENDDO
     ELSEIF (INT_TYPE.GE.5) THEN
        WRITE(FILE_NUM(INT_TYPE),'(a)') FILE_NAME(INT_TYPE) 
        DO I = 1,NBFT
           DO J = 1,NBFT
              WRITE(FILE_NUM(INT_TYPE),'(i7,i7,2x,f20.16)') I,J,ARRAY_SO(J,I)
           ENDDO
        ENDDO
     ENDIF
     CLOSE(UNIT=FILE_NUM(INT_TYPE))
  ENDIF
 
  100  CONTINUE

  CALL OUTPUT_LOC('RECP_WRITE_FILE','X')

END SUBROUTINE RECP_WRITE_FILE


! ------------------------------------------------


SUBROUTINE RECP_WRITE_TYPE(itypea, itypeb, int_type, lblop,PRINT_LEVEL)
! check integral type 
! Y. C. Park
  IMPLICIT NONE

  INTEGER :: PRINT_LEVEL
  INTEGER :: itypea,  itypeb, int_type, lblop

  IF(PRINT_LEVEL.GE.5)  PRINT *,' * ECPINT : entering... ARGOS_WRITE_TYPE * '

! int_type = 1-4 : symmetric matrix, 
!                  totally symmetric operator
! int_type =     : symmetric matrix, 
!                  nontotally symmetric operator
! int_type = 5-7 : antisymmetric matrix, 
!                  nontotally symmetric operator
  if (itypea.eq.0) then
     if     (itypeb.eq.0) then
        int_type = 1        !   S
     elseif (itypeb.eq.1) then
        int_type = 2        !   T 
     elseif (itypeb.eq.2) then
        int_type = 3        !   V
     elseif (itypeb.eq.3) then
        int_type = 4        !   Core or Veff
     endif
  elseif (itypea.eq.1) then
     print *, 'This integral type is not supported for the now'
     stop
  elseif (itypea.eq.2) then
     if     (lblop.eq.1) then
        int_type = 5         !   SO(x)
     elseif (lblop.eq.2) then
        int_type = 6         !   SO(y)
     elseif (lblop.eq.3) then
        int_type = 7         !   SO(z)
     endif
  endif

  IF (PRINT_LEVEL.GE.5)  PRINT *,' * ECPINT : exiting...  ARGOS_WRITE_TYPE * '

END SUBROUTINE RECP_WRITE_TYPE

END MODULE RECP_WRITE
