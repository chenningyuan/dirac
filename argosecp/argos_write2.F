*deck stvcz_write
      SUBROUTINE ARGOS_WRITE(
     & nop,iapt,ipr,ipr1,ipr2,num,
     & il,array,iarray,n_bfn_sym,chrtyp,
     & N_INTEGRAL,MN_INTEGRAL,int_AREP)

!     write 1-e integral arrays
!     modified from original argos
!     Y.C. Park (Prof. Y.S.Lee's group, KAIST)

!     variables
#include "argos_mxvalue.h"
#include "argos_files.h"
#include "argos_error.h"
#include "argos_print.h"

      integer    ilfact
      parameter( ilfact=1024 )
      integer    nipv
      parameter( nipv=2 )
      integer  INT_AREP

      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      integer
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir(3),      inam(5),     nnam, mdum(32)

!     # /bufout/ holds some output integral file parameters.
      integer
     & info,     ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
      common /bufout/
     & info(10), ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape

      integer       n1max
      equivalence ( n1max, info(3) )
      integer lblop
      integer    msame,   nmsame,   nomore
      parameter( msame=0, nmsame=1, nomore=2 )
      real*8     zero
      parameter( zero=0d0 )

      integer N_INTEGRAL(3) ! number of integral element of each type
      integer MN_INTEGRAL   ! maximum of N_INTEGRAL(i,i=1-3)
      integer iarray(mpru,*), il(nnbft)
      real*8 array(mpru,*),values_stvc(MN_INTEGRAL)
      real*8  array_stvc(MN_INTEGRAL)
      integer labels_stvc(2,MN_INTEGRAL) !from labels(nipv,*)
      real*8 array_so(nbft,nbft)
      character*8 chrtyp
      integer nop, iapt(3), ipr, ipr1, ipr2, num

      integer ibuf_temp
      integer n_bfn_sym(mstup)  ! (nbpsy or nso)

      integer itypan(7), itypbn(7)
!     #             s    t    v  veff so:x so:y so:z
      data itypan / 0,   0,   0,   0,   2,   2,   2 /
      data itypbn / 0,   1,   2,   3,   0,   1,   2 /

      IF (PRINT_LEVEL .GE. 5) THEN
         PRINT *,' * ECPINT : entering... ARGOS_WRITE'
         PRINT *,'   nop            :',nop
      ENDIF

!     # write out the components of the operator arrays.

      do lblop = 1, nop

!       # check that array(*) was correctly allocated.
        if ( iapt(lblop) .gt. mpru ) then
           call bummer('oneint: array(*) overflow, iapt=',
     &      iapt(lblop), faterr )
        endif

!       # initialize the integral record counters for this array.
!       call w1stup
        ibuf   = 0
        ibuf_temp = 0

        IF (PRINT_LEVEL.GE.10) THEN
           print *,'   after w1stup'
           print *,'   + ibuf      = ',ibuf
           print *,'   + ibuf_temp = ',ibuf_temp
        ENDIF

!ycp-e
!       # SIFS integral type.
        itypea = itypan(nnam + lblop - 1)
        itypeb = itypbn(nnam + lblop - 1)

!       # loop over output records.
!       # do while ( ipr1 .le. iapt(lblop) )...
        ipr1 = 1
820     if ( ipr1 .le. iapt(lblop) ) then

          num  = min( (n1max - ibuf_temp), (iapt(lblop) - ipr1 + 1) )
          ipr2 = ipr1 + num - 1

          IF (PRINT_LEVEL.GE.10) THEN
             print *, ' * ARGOS_WRITE : in the loop' 
             print *, ' ---------------------------' 
             print *, '   n1max              = ',n1max
             print *, '   ibuf               = ',ibuf
             print *, '   ibuf_temp          = ',ibuf_temp
             print *, ' 1:n1max - ibuf_temp  =',n1max-ibuf_temp
             print *, ' 2:iapt(lblop)-ipr1+1 =',iapt(lblop)-ipr1+1
             print *, ' num=(min(1,2))       =', num 
             print *, ' ---------------------------' 
             print *, ' N_INTEGRAL         = ', N_INTEGRAL(LBLOP)
             print *, ' iapt(lblop)        = ', iapt(lblop) 
             print *, ' ipr1               = ',ipr1
             print *, ' ipr2               = ',ipr2
             print *, ' ---------------------------' 
          ENDIF

!         # fill up the record.
          do ipr = ipr1, ipr2
            ibuf           = ibuf + 1
            ibuf_temp      = ibuf_temp + 1
            labels_stvc(1,ibuf)=    il(iarray(ipr,lblop))/ilfact
            labels_stvc(2,ibuf)=mod(il(iarray(ipr,lblop)),ilfact)
            values_stvc(ibuf)   = array(ipr,lblop)

!           # print values
!           IF (PRINT_LEVEL.GE.10) THEN
!              WRITE(NLIST,'(a,i6,i4,i4,f20.16)')
!    &              '   ARGOS_WRITE : ',ibuf,
!    &              labels(1,ibuf),labels(2,ibuf),values_stvc(ibuf)
!           ENDIF
          enddo

          ipr1 = ipr2 + 1
          ibuf_temp = 0 !ycp wtlab1->sifew1 (ibuf=0)

          go to 820
        endif
 
!       # check integral type
        call argos_write_type(itypea, itypeb, int_type, lblop)

        IF (PRINT_LEVEL.GE.10) THEN
           PRINT *,'   INT_TYPE =', INT_TYPE
           PRINT *,'   NST      =', NST
           PRINT *,'   LBLOP    =', LBLOP
           PRINT *,'   N1MAX    =', N1MAX
           PRINT *,'   IBUF           :',IBUF
           PRINT *,'   N_INTEGRAL     :',N_INTEGRAL
        ENDIF

        call argos_write_array(
     &       nst, int_type, n_bfn_sym, nbft, labels_stvc,n1max,
     &       mpru, array_stvc, array_so, values_stvc,
     &       N_INTEGRAL,MN_INTEGRAL, lblop)
        call argos_write_file( 
     &       int_type, N_INTEGRAL, nbft,array_stvc, array_so,
     &       INT_AREP,MN_INTEGRAL)

        IF (PRINT_LEVEL .GE. 5) THEN
           call siftyp( itypea, itypeb, chrtyp )
           write (nlist,'(4x,a,a)') chrtyp,
     &           ' integrals were written in a file'
        ENDIF

      enddo

      IF (PRINT_LEVEL .GE. 5) THEN
!        PRINT *,'   N_INTEGRAL(',LBLOP,')',N_INTEGRAL
         PRINT *,' * ECPINT : exiting...  ARGOS_WRITE'
      ENDIF

      return
      end



! ### From cnvrt_stvz ###
      subroutine argos_write_array(
     &           n_sym, int_type, n_bfn_sym,nbft,labels_stvc,n1max,
     &           mpru, array_stvc, array_so, values_stvc, 
     &           N_INTEGRAL,MN_INTEGRAL,lblop)

#include "argos_mxvalue.h"
#include "argos_print.h"

!     # /bufout/ holds some output integral file parameters.
      integer
     & info,     ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
      common /bufout/
     & info(10), ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape

      integer  n_bfn_sym(mstup)! nbpsy or nso
!                                number of basis functions 
!                                per symmetry block

      integer  int_type     ! integral type 
      integer  n_sym        ! nsym or nst
!                             number of symmetry blocks
      integer  nbft
      integer  mpru
      integer  lblop
!     real*8   array_stvc(mpru)
      integer  N_INTEGRAL(3), MN_INTEGRAL
      real*8   array_stvc(MN_INTEGRAL) 
      real*8   array_so(nbft,nbft) 

      integer       n1max
!     equivalence ( n1max, info(3) )
!     real*8        values_stvc(n1max) 
      real*8        values_stvc(MN_INTEGRAL)

      integer  labels_stvc(2,*) ! labels(1:2,1:n1max) orbital label buffer. 

!     variables : local
      integer  i, j, ij
      integer  n_tot        ! n_total
      integer  n_n_tot      ! n! of n_total
      integer  n_skp(8)     ! n_skip
      integer  n_n_skp(8)   ! n! of n_skip
      integer    nbfmx
      parameter (nbfmx=255)    ! tempory 
      integer   sym_bfn(nbfmx) ! symb(1:nbft)
!                                symmetry index of each basis function
!                                'n_bfn_sym' is needed
      integer i_lab, j_lab
      integer i_sym, j_sym, ij_sym   ! isym, jsym, ijsym

      integer nndxf
      nndxf(i) = ( i * ( i - 1 ) ) / 2     ! define nndxf

      IF (PRINT_LEVEL .GE. 5) THEN 
         PRINT *,' * ECPINT : entering... ARGOS_WRITE_ARRAY * '
         PRINT *,'   INT_TYPE       =',INT_TYPE
         PRINT *,'   N_INTEGRAL     =',N_INTEGRAL
         PRINT *,'   N_SYM          =',N_SYM
         PRINT *,'   MN_INTEGRAL    =',MN_INTEGRAL
         PRINT *,'   N1MAX          =',N1MAX
      ENDIF

!     make sym_bfn 
      n_tot  = 0
      n_n_tot = 0
      do i_sym = 1, n_sym
         n_skp(i_sym)   = n_tot  
         n_n_skp(i_sym) = n_n_tot 
         n_tot          = n_tot   + n_bfn_sym(i_sym)               
         n_n_tot        = n_n_tot + nndxf(n_bfn_sym(i_sym) + 1)  
!        print *, '1-n_bfn_sym :',i_sym,n_bfn_sym(i_sym)
         do j = (n_skp(i_sym)+1), n_tot
            sym_bfn(j)  =i_sym
         enddo
      enddo

      if (int_type .le. 4) then
!     # symmetric matrix, totally symmetric operator.

!        #  initialize matrix element 
         array_stvc = 0

!        do i = 1, ibuf 
         do i = 1, N_INTEGRAL(1)
           i_lab  = max( labels_stvc(1,i), labels_stvc(2,i) )
           j_lab  = min( labels_stvc(1,i), labels_stvc(2,i) )
           i_sym  = sym_bfn(i_lab)
           j_sym  = sym_bfn(j_lab)
           ij_sym = nndxf(i_sym) + j_sym
           if ( i_sym .eq. j_sym ) then
!             # only valid array(*) elements are referenced.
              ij  =  n_n_skp(i_sym)             ! skip for each symmetry
     &             + nndxf(i_lab - n_skp(i_sym))
     &             + j_lab
     &             - n_skp(i_sym)
              array_stvc(ij) = values_stvc(i)
           endif
         enddo

         IF (PRINT_LEVEL.GE.10) THEN 
            do i = 1, N_INTEGRAL(1)
               print *,'array_stvc',i,array_stvc(i)
            enddo
         ENDIF

      elseif ( int_type .ge. 5 ) then
!     # antisymmetric matrix, nontotally symmetric operator.

!        # initialize spin-orbit matrix element 
         array_so = 0

!        do i = 1, ibuf
         do i = 1, N_INTEGRAL(LBLOP)
            i_lab  = labels_stvc(1,i)
            j_lab  = labels_stvc(2,i)
            if (i_lab .lt. j_lab) then
               i_lab = labels_stvc(2,i)
               j_lab = labels_stvc(1,i)
!              the sign change upon transposition.
               values_stvc(i) = -values_stvc(i)
            endif
            array_so(i_lab,j_lab) =  values_stvc(i)
            array_so(j_lab,i_lab) = -values_stvc(i)
!           print *,'array_so',i_lab,j_lab,array_so(i_lab,j_lab)
         enddo

      endif

      IF (PRINT_LEVEL .GE. 5) 
     &   PRINT *,' * ECPINT : exiting...  ARGOS_WRITE_ARRAY * '
      
      return
      end




      subroutine argos_write_file( 
     &           int_type, n_integral, nbft,
     &           array_stvc, array_so, INT_AREP, MN_INTEGRAL)
!
#include "argos_print.h"

      integer i,j
      integer int_type, INT_AREP
      integer N_INTEGRAL, MN_INTEGRAL
      integer nbft
      real*8  array_stvc(MN_INTEGRAL)
!     real*8  array_stvc(*)
      real*8  array_so(nbft,nbft)

      IF (PRINT_LEVEL.GE.5) 
     &   PRINT *,' * ECPINT : entering... ARGOS_WRITE_FILE * '

!     Write integrals to file & Print

      if     (int_type .eq. 1) then  
!        write(6,'(a)')  '    integral type : S       '
!        open (59,form='unformatted',file='RECP_S')
!        write(59),(array_stvc(i),i=1,n_integral)
!        close(59)
      elseif (int_type .eq. 2) then  
!        write(6,'(a)')  '    integral type : T       '
!        open (59,form='unformatted',file='RECP_T')
!        write(59),(array_stvc(i),i=1,n_integral)
!        close(59)
      elseif (int_type .eq. 3) then  
!        write(6,'(a)')  '    integral type : V       '
!        open (59,form='unformatted',file='RECP_V')
!        write(59),(array_stvc(i),i=1,n_integral)
!        close(59)

! ### RECP_VCORE ###
      elseif (int_type .eq. 4) then 
!        debug  --------------------
         IF (ARGOS_DEBUG.GE.1) THEN 
            open (64,form='formatted',file='RECP_DBG_VCORE')
            write(64,'(a)')  '    integral type : Veff    '
            do i = 1,  n_integral
               write(64,1001) 'Veff:',i,array_stvc(i)
            enddo
            close(64)
         ENDIF
!        write --------------------
         open (64,form='unformatted',file='RECP_VCORE')
         write(64),(array_stvc(i),i=1,n_integral)
         close(64)

         IF (INT_AREP.eq.1) THEN
          array_so=0
          open(65,form='unformatted',file='RECP_SOX')
          open(66,form='unformatted',file='RECP_SOY')
          open(67,form='unformatted',file='RECP_SOZ')
          write(65),((array_so(j,i),j=1,nbft),i=1,nbft)
          write(66),((array_so(j,i),j=1,nbft),i=1,nbft)
          write(67),((array_so(j,i),j=1,nbft),i=1,nbft)
          close(65)
          close(66)
          close(67)
         ENDIF

! ### RECP_SOX ###
      elseif (int_type .eq. 5) then  
!        debug  --------------------
         IF (ARGOS_DEBUG.GE.1) THEN   
            open(65,form='formatted',file='RECP_DBG_SOX')
            write(65,'(a)')  '    integral type : SO(x)   '
            do i = 1,  nbft
               do j = 1,  nbft
                  write(65,1002) 'SO(x):',i,j,array_so(j,i)
               enddo
            enddo
            close(65)
         ENDIF
!        write --------------------
         open(65,form='unformatted',file='RECP_SOX')
         write(65),((array_so(j,i),j=1,nbft),i=1,nbft)
         close(65)

! ### RECP_SOY ###
      elseif (int_type .eq. 6) then  
!        debug  --------------------
         IF (ARGOS_DEBUG.GE.1) THEN   
            open(66,form='formatted',file='RECP_DBG_SOY')
            write(66,'(a)')  '    integral type : SO(y)   '
            do i = 1,  nbft
               do j = 1,  nbft
                  write(66,1002) 'SO(y):',i,j,array_so(j,i)
               enddo
            enddo
            close(66)
         ENDIF
!        write --------------------
         open(66,form='unformatted',file='RECP_SOY')
         write(66),((array_so(j,i),j=1,nbft),i=1,nbft)
         close(66)

! ### RECP_SOZ ###
      elseif (int_type .eq. 7) then  
!        debug  --------------------
         IF (ARGOS_DEBUG.GE.1) THEN
            open(67,form='formatted',file='RECP_DBG_SOZ')
            write(67,'(a)')  '    integral type : SO(z)   '
            do i = 1,  nbft
               do j = 1,  nbft
                  write(67,1002) 'SO(z):',i,j,array_so(j,i)
               enddo
            enddo
            close(67)
         ENDIF
!        write --------------------
         open(67,form='unformatted',file='RECP_SOZ')
         write(67),((array_so(j,i),j=1,nbft),i=1,nbft)
         close(67)
      endif

      IF (PRINT_LEVEL.GE.5) 
     &    PRINT *,' * ECPINT : exiting...  ARGOS_WRITE_FILE * '

 1001 format(a,2x,i5,2x,f20.16)
 1002 format(a,2x,i5,i5,2x,f20.16)
      
      end


      subroutine argos_write_type(itypea, itypeb, int_type, lblop)
!     check integral type 

      integer    itypea,  itypeb, int_type, lblop
      IF (PRINT_LEVEL.GE.5) 
     &   PRINT *,' * ECPINT : entering... ARGOS_WRITE_TYPE * '


!     int_type = 1-4 : symmetric matrix, 
!                      totally symmetric operator
!     int_type =     : symmetric matrix, 
!                      nontotally symmetric operator
!     int_type = 5-7 : antisymmetric matrix, 
!                      nontotally symmetric operator
      if (itypea.eq.0) then
         if     (itypeb.eq.0) then
            int_type = 1         !   S
         elseif (itypeb.eq.1) then
            int_type = 2         !   T 
         elseif (itypeb.eq.2) then
            int_type = 3         !   V
         elseif (itypeb.eq.3) then
            int_type = 4         !   Core or Veff
         endif
      elseif (itypea.eq.1) then
         print *, 'This integral type is not supported now'
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

      IF (PRINT_LEVEL.GE.5) 
     &   PRINT *,' * ECPINT : exiting...  ARGOS_WRITE_TYPE * '

      end

