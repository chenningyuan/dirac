MODULE ARGOS_WRITE
  IMPLICIT NONE

CONTAINS 

! ------------------------------------------------

SUBROUTINE ARGOS_WRITE_INTERMEDIATE(ircru,jrcru,nst,mstu,nprir,ijsf,nopir,ilxyz,iapt, &
                                    mpru,iarray,ipr,N_INTEGRAL,array,h2,esf)

! # put integrals and labels into the intermediate storage arrays.


! global variables
  INTEGER :: ircru, jrcru
  INTEGER :: nst
  INTEGER :: mstu
  INTEGER :: nprir(2,mstu,*)
  INTEGER :: ijsf
  INTEGER :: nopir(*)
  INTEGER :: ilxyz(3,*)
  INTEGER :: iapt(3)
  INTEGER :: mpru
  INTEGER :: iarray(mpru,*)
  INTEGER :: ipr
  INTEGER :: N_INTEGRAL(3)     ! number of integral element of each type 

  REAL(8) :: array(mpru,*)
  REAL(8) :: h2(*)
  LOGICAL :: esf
! local variables
  INTEGER :: ist
  INTEGER :: iop
  INTEGER :: INDEX1
  INTEGER :: NPAIR, IPAIR
  INTEGER :: I, J
  INTEGER :: INT_TYPE

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
              iapt(INT_TYPE)                 = iapt(INT_TYPE) + 1
              array(  iapt(INT_TYPE), INT_TYPE) = h2(INDEX1)
              iarray( iapt(INT_TYPE), INT_TYPE) = (ipr+IPAIR)
!             print *, 'ARGOS_ONEINT ',
!                      iapt(INT_TYPE),
!                      INT_TYPE,       ! S,T,V=1, Veff=1,2,3
!                      INDEX1, h2(INDEX1)
!                      ,array(iapt(INT_TYPE), INT_TYPE)
!                      ,iarray(iapt(INT_TYPE), INT_TYPE)
              N_INTEGRAL(INT_TYPE) = N_INTEGRAL(INT_TYPE) + 1
            enddo
          enddo
        endif

        ipr = ipr + NPAIR

      enddo
    enddo
  enddo

END SUBROUTINE ARGOS_WRITE_INTERMEDIATE

! ------------------------------------------------

SUBROUTINE ARGOS_WRITE_MAIN(nop,iapt,ipr,ipr1,ipr2,N1MAX,L1REC,il,array,iarray,n_bfn_sym, &
                            N_INTEGRAL,MN_INTEGRAL,INT_AREP,nnam,mpru,nbft,nnbft,nst,PRINT_LEVEL)

! write 1-e integral arrays
! Y.C. Park (Prof. Y.S.Lee's group, KAIST)

! variables
#include "argos_mxvalue.h"
#include "argos_files.h"
#include "argos_error.h"

! global variable
  INTEGER :: PRINT_LEVEL
  INTEGER :: mpru, nbft, nnbft, nst   ! common block (parmi)
  INTEGER :: nnam                     ! common block (ntgr)
  INTEGER :: INT_AREP
  INTEGER :: N1MAX,  L1REC

  INTEGER :: N_INTEGRAL(3) ! number of integral element of each type
  INTEGER :: MN_INTEGRAL   ! maximum of N_INTEGRAL(i,i=1-3)

  INTEGER :: nop, iapt(3), ipr, ipr1, ipr2
  INTEGER :: iarray(mpru,*), il(nnbft)
  REAL(8) :: array(mpru,*)

! local variable
  REAL(8), ALLOCATABLE :: array_stvc(:)
  REAL(8), ALLOCATABLE :: array_so(:,:)
  REAL(8), ALLOCATABLE :: values_stvc(:)    ! to be fixed
  INTEGER, ALLOCATABLE :: labels_stvc(:,:)  ! from labels(nipv,*)

  INTEGER :: int_type     ! integral type 
  INTEGER :: ilfact = 1024
  INTEGER :: itypea, itypeb
  INTEGER :: lblop
  INTEGER :: ibuf, ibuf_temp
  INTEGER :: n_bfn_sym(mstup)  ! (nbpsy or nso)
  INTEGER :: num
  CHARACTER(8) :: chrtyp

  INTEGER :: itypan(7), itypbn(7)
! #             s    t    v  veff so:x so:y so:z
  data itypan / 0,   0,   0,   0,   2,   2,   2 /
  data itypbn / 0,   1,   2,   3,   0,   1,   2 /

  ALLOCATE( array_stvc(MN_INTEGRAL) )
  ALLOCATE( array_so(nbft,nbft) )
  ALLOCATE( values_stvc(MN_INTEGRAL) )
  ALLOCATE( labels_stvc(2,MN_INTEGRAL) )

  IF (PRINT_LEVEL .GE. 5) THEN
     PRINT *,' * ECPINT : entering... ARGOS_WRITE'
     PRINT *,'   nop            :',nop
  ENDIF

! # write out the components of the operator arrays.

  do lblop = 1, nop

!    check that array(*) was correctly allocated.
     if (iapt(lblop).gt.mpru)   call bummer('oneint: array(*) overflow, iapt=',iapt(lblop), faterr )

!    initialize the integral record counters for this array.
     ibuf   = 0
     ibuf_temp = 0

     IF (PRINT_LEVEL.GE.10) THEN
        print *,'   after w1stup'
        print *,'   + ibuf      = ',ibuf
        print *,'   + ibuf_temp = ',ibuf_temp
     ENDIF

!    # check integral type
     itypea = itypan(nnam + lblop - 1)
     itypeb = itypbn(nnam + lblop - 1)
     call siftyp( itypea, itypeb, chrtyp )
     call argos_write_type(itypea, itypeb, int_type, lblop,PRINT_LEVEL)

!    # loop over output records.
!    # do while ( ipr1 .le. iapt(lblop) )...
     ipr1 = 1

     820 CONTINUE
     if( ipr1 .le. iapt(lblop) ) then
       num  = min( (n1max - ibuf_temp), (iapt(lblop) - ipr1 + 1) )
       ipr2 = ipr1 + num - 1
!      print information
       IF (PRINT_LEVEL.GE.10) CALL ARGOS_WRITE_INFO1(ibuf,ibuf_temp,N1MAX,lblop,iapt,ipr1,ipr2,num,N_INTEGRAL)
!      fill up the record.
       do ipr = ipr1, ipr2
          ibuf           = ibuf + 1
          ibuf_temp      = ibuf_temp + 1
          labels_stvc(1,ibuf)=    il(iarray(ipr,lblop))/ilfact
          labels_stvc(2,ibuf)=mod(il(iarray(ipr,lblop)),ilfact)
          values_stvc(ibuf)   = array(ipr,lblop)
!         # print values
!         IF (PRINT_LEVEL.GE.10) WRITE(NLIST,'(a,i6,i4,i4,f20.16)') &
!                        ' ARGOS_WRITE : ',ibuf, labels(1,ibuf),labels(2,ibuf),values_stvc(ibuf)
       enddo
       ipr1 = ipr2 + 1
       ibuf_temp = 0 !ycp wtlab1->sifew1 (ibuf=0)
       GOTO 820
    endif
 
    IF(PRINT_LEVEL.GE.10) CALL ARGOS_WRITE_INFO2(int_type,nst,lblop,N1MAX,ibuf,N_INTEGRAL) 

!   # write integral to file
    call argos_write_array( nst,int_type,n_bfn_sym,nbft,labels_stvc,mpru,array_stvc,array_so, &
                            values_stvc,N_INTEGRAL,MN_INTEGRAL,lblop,N1MAX,PRINT_LEVEL)
    call argos_write_file( int_type, N_INTEGRAL, nbft,array_stvc, array_so,INT_AREP,MN_INTEGRAL, &
                           PRINT_LEVEL,LBLOP)

    IF(PRINT_LEVEL.GE. 5) write (nlist,'(4x,a,a)') chrtyp,' integrals were written in a file'
    IF(PRINT_LEVEL.GE.10) PRINT *,'   N_INTEGRAL(',LBLOP,')',N_INTEGRAL(LBLOP)

  enddo

! print information
  IF(PRINT_LEVEL.GE. 5) PRINT *,' * ECPINT : exiting...  ARGOS_WRITE'

  DEALLOCATE( array_stvc )
  DEALLOCATE( array_so )
  DEALLOCATE( values_stvc )
  DEALLOCATE( labels_stvc )

END SUBROUTINE ARGOS_WRITE_MAIN

! ------------------------------------------------


SUBROUTINE ARGOS_WRITE_ARRAY(n_sym,int_type,n_bfn_sym,nbft,labels_stvc,mpru,array_stvc,array_so, & 
                             values_stvc,N_INTEGRAL,MN_INTEGRAL,lblop,N1MAX,PRINT_LEVEL)
! From cnvrt_stvz 
! Y. C. Park

#include "argos_mxvalue.h"

  INTEGER :: PRINT_LEVEL
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
  INTEGER :: N1MAX


! variables : local
  integer :: i, j, ij
  integer :: n_tot        ! n_total
  integer :: n_n_tot      ! n! of n_total
  integer :: n_skp(8)     ! n_skip
  integer :: n_n_skp(8)   ! n! of n_skip

  integer :: i_lab, j_lab
  integer :: i_sym, j_sym, ij_sym   ! isym, jsym, ijsym

  integer :: nbfmx=500    ! tempory 
  integer, ALLOCATABLE :: sym_bfn(:) ! symb(1:nbft) symmetry index of each basis function 'n_bfn_sym' is needed

  integer :: nndxf
  nndxf(i) = ( i * ( i - 1 ) ) / 2     ! define nndxf

  ALLOCATE ( sym_bfn(nbfmx) )

  IF (PRINT_LEVEL .GE. 5) THEN 
     PRINT *,' * ECPINT : entering... ARGOS_WRITE_ARRAY * '
     PRINT *,'   INT_TYPE       =',INT_TYPE
     PRINT *,'   N_INTEGRAL     =',N_INTEGRAL(LBLOP)
     PRINT *,'   N_SYM          =',N_SYM
     PRINT *,'   MN_INTEGRAL    =',MN_INTEGRAL
     PRINT *,'   N1MAX          =',N1MAX
  ENDIF

! make sym_bfn 
  n_tot  = 0
  n_n_tot = 0
  do i_sym = 1, n_sym
     n_skp(i_sym)   = n_tot  
     n_n_skp(i_sym) = n_n_tot 
     n_tot          = n_tot   + n_bfn_sym(i_sym)               
     n_n_tot        = n_n_tot + nndxf(n_bfn_sym(i_sym) + 1)  
!    print *, '1-n_bfn_sym :',i_sym,n_bfn_sym(i_sym)
     do j = (n_skp(i_sym)+1), n_tot
        sym_bfn(j)  =i_sym
     enddo
  enddo

  if (int_type .le. 4) then
! # symmetric matrix, totally symmetric operator.

!    #  initialize matrix element 
     array_stvc = 0

!    do i = 1, ibuf 
     do i = 1, N_INTEGRAL(1)
       i_lab  = max( labels_stvc(1,i), labels_stvc(2,i) )
       j_lab  = min( labels_stvc(1,i), labels_stvc(2,i) )
       i_sym  = sym_bfn(i_lab)
       j_sym  = sym_bfn(j_lab)
       ij_sym = nndxf(i_sym) + j_sym
       if ( i_sym .eq. j_sym ) then
!         # only valid array(*) elements are referenced.
          ij  =  n_n_skp(i_sym) + nndxf(i_lab - n_skp(i_sym)) + j_lab - n_skp(i_sym)
!       | ij  :  skip for       +                             +       -              |
!       |        each symmetry                                                       | 
          array_stvc(ij) = values_stvc(i)
       endif
     enddo

     IF (PRINT_LEVEL.GE.10) THEN 
        do i = 1, N_INTEGRAL(1)
           print *,'array_stvc',i,array_stvc(i)
        enddo
     ENDIF

  elseif ( int_type .ge. 5 ) then
! # antisymmetric matrix, nontotally symmetric operator.

!    # initialize spin-orbit matrix element 
     array_so = 0

!    do i = 1, ibuf
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

  IF (PRINT_LEVEL .GE. 5) PRINT *,' * ECPINT : exiting...  ARGOS_WRITE_ARRAY * '
  
END SUBROUTINE ARGOS_WRITE_ARRAY


! ------------------------------------------------


SUBROUTINE ARGOS_WRITE_FILE(int_type,N_INTEGRAL,nbft,array_stvc,array_so, &
                            INT_AREP,MN_INTEGRAL,PRINT_LEVEL,LBLOP)
! write ro file
! Y. C. Park

  INTEGER :: PRINT_LEVEL
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

! # start subroutine
  IF(PRINT_LEVEL.GE.5)  PRINT *,' * ECPINT : entering... ARGOS_WRITE_FILE * '

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
! IF (ARGOS_DEBUG) THEN
!    OPEN (UNIT=FILE_NUM(INT_TYPE),FILE=FILE_DBG(INT_TYPE), FORM='FORMATTED')
!    IF (INT_TYPE.EQ.4) THEN
!       WRITE(FILE_NUM(INT_TYPE),'(a)') FILE_NAME(INT_TYPE) 
!       DO I = 1,N_INTEGRAL
!          WRITE(FILE_NUM(INT_TYPE),'(i7,2x,f20.16)') I,ARRAY_STVC(I)
!       ENDDO
!    ELSEIF (INT_TYPE.GE.5) THEN
!       WRITE(FILE_NUM(INT_TYPE),'(a)') FILE_NAME(INT_TYPE) 
!       DO I = 1,NBFT
!          DO J = 1,NBFT
!             WRITE(FILE_NUM(INT_TYPE),'(i7,i7,2x,f20.16)') I,J,ARRAY_SO(J,I)
!          ENDDO
!       ENDDO
!    ENDIF
!    CLOSE(UNIT=FILE_NUM(INT_TYPE))
! ENDIF
 
  100  CONTINUE

  IF (PRINT_LEVEL.GE.5)  PRINT *,' * ECPINT : exiting...  ARGOS_WRITE_FILE * '

END SUBROUTINE ARGOS_WRITE_FILE


! ------------------------------------------------


SUBROUTINE ARGOS_WRITE_TYPE(itypea, itypeb, int_type, lblop,PRINT_LEVEL)
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

END SUBROUTINE ARGOS_WRITE_TYPE


! ------------------------------------------------


SUBROUTINE ARGOS_WRITE_INFO1(ibuf,ibuf_temp,N1MAX,lblop,iapt,ipr1,ipr2,num,N_INTEGRAL)
! Y. C. Park

  IMPLICIT NONE

  INTEGER :: ibuf, ibuf_temp
  INTEGER :: N1MAX
  INTEGER :: lblop
  INTEGER :: iapt(3) 
  INTEGER :: ipr1, ipr2
  INTEGER :: num
  INTEGER :: N_INTEGRAL(3) ! number of integral element of each type

  print *, ' * ARGOS_WRITE : in the loop'
  print *, ' ---------------------------'
! print *, '   n1max              = ',n1max
  print *, '   ibuf               = ',ibuf
  print *, '   ibuf_temp          = ',ibuf_temp
! print *, ' 1:n1max - ibuf_temp  =',n1max-ibuf_temp
  print *, ' 2:iapt(lblop)-ipr1+1 =',iapt(lblop)-ipr1+1
  print *, ' num=(min(1,2))       =', num
  print *, ' ---------------------------'
  print *, ' N_INTEGRAL         = ', N_INTEGRAL(LBLOP)
  print *, ' iapt(lblop)        = ', iapt(lblop)
  print *, ' ipr1               = ',ipr1
  print *, ' ipr2               = ',ipr2
  print *, ' ---------------------------'
  
END SUBROUTINE ARGOS_WRITE_INFO1 


! ------------------------------------------------


SUBROUTINE ARGOS_WRITE_INFO2(int_type,nst,lblop,N1MAX,ibuf,N_INTEGRAL) 
! Y. C. Park

  IMPLICIT NONE

  INTEGER :: int_type      ! integral type 
  INTEGER :: nst           !
  INTEGER :: lblop         !
  INTEGER :: N1MAX         !
  INTEGER :: ibuf          !
  INTEGER :: N_INTEGRAL(3) ! number of integral element of each type

  PRINT *,'   INT_TYPE   :', int_type
  PRINT *,'   NST        :', NST
  PRINT *,'   LBLOP      :', lblop
  PRINT *,'   N1MAX      :', N1MAX
  PRINT *,'   IBUF       :',IBUF
  PRINT *,'   N_INTEGRAL :',N_INTEGRAL(LBLOP)

END SUBROUTINE ARGOS_WRITE_INFO2 

END MODULE
