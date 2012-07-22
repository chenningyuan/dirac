      subroutine csd
!---------------------------------------------------------------------
!       closed shell CC method (CCSD, CCSDT-1,...), 00 sector
!---------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use gen_mod
      use cfiles_mod
      use cint_mod
      use dat_mod
      use cij1_mod
      use cij2_mod
      use hsfs_parameters
implicit none
#include "priunit.h"
#include "csdtime.h"
#include "dbg.h"
#else
implicit none
#include "priunit.h"
#include "param_ccsd.h"
#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)
#include "cfiles.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)
#include "cint.h"
!      common /cint/ nint,idum1,vint(ndint),iint(ndint)
#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "csdtime.h"
!     common /csdtime/ csdstim,csdltim
#include "dbg.h"
!     common /dbg/ iter,idebug
#endif
character*256 :: filenm,filenm51,filenm6,filenm20
real(kind=8)::convrg,tmax,t1max,t3max,t4max, diff,df3,df4,fnew,diff1
integer::ntape4,new,itmax,it1max,id3mx,ind3mx,id4mx,ind4mx, & 
         iddf3,inddf3,idiff,ising,icnvg,inddf4,idiff1,itr1,iddf4
integer, external :: idec

!     iopt(1) - print.
!     iopt(2) - no. of iterations.
!     iopt(3) - convergence = 10**(-iopt(3))
!     iopt(4) .gt.0 - do rle every iopt(4) iterations.
!             .lt.0 - diis
!     iopt(5) - print only elements gt 10**(-iopt(5))
!     iopt(6) - start avrg.
!     iopt(8) - fock.
!     iopt(10) .ne.0 - do not save pppp.

      real*4 cputime

! open some files
      filenm(1:1)=' '
      call getenv('input',filenm)
      if (filenm(1:1).eq.' ') filenm='input'
      open(5,file=filenm)

      filenm6(1:1)=' '
      call getenv('output',filenm6)
      if (filenm6(1:1).ne.' ') open(6,file=filenm6,form='FORMATTED',    &
  status='UNKNOWN',access='SEQUENTIAL')

      filenm51(1:1)=' '
      call getenv('tape51',filenm51)
      if (filenm51(1:1).eq.' ') filenm51='tape51'
      open(51,file=filenm51,form='UNFORMATTED')

      filenm20(1:1)=' '
      call getenv('tape20',filenm20)
      if (filenm20(1:1).eq.' ') filenm20='tape20'
      open(20,file=filenm20,form='UNFORMATTED')

      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')

      call header(                                                      &
  'Standalone spin-adapted CCSD/CCSDT program,'//                       &
 ' precious heritage of the TAU group',-1)

      csdstim = cputime(real(0.))
      csdltim = csdstim

!    ... miro: important initialization (not all compilers assign zero
!    by default)
      nflcur=0
      Nstrt(1,0)=0
      nstrt(2,0)=0

! ... readin with ccsd.F routines ...
      call rdin
      call rdsym

      if (iopt(10).eq.0) then
        filenm(1:1)=' '
        call getenv('tape9',filenm)
        if (filenm(1:1).eq.' ') filenm='tape9'
        open(9,file=filenm,form='UNFORMATTED')
      endif

      ! convergence criterium 
      convrg = 10.d0**(-iopt(3))
      write (LUPRI,'('' max # of iter='',i3,''. convrg='',e10.2)') iopt(2),convrg

      call prtime("rdin",0,0)

      call ordv_csd ! sorting integrals

      iter=0

      call tinit ! initialization of amplitudes

      if (iopt(17).ne.0) then
        call files_ccsd ("t30 ",ntape4,new,1,-1,-1)
        write (6,*) 'length of three body file = ',ije1(nij1)
      end if
      if (iopt(17).ge.4) then
        call files_ccsd ("t40 ",ntape4,new,1,-1,-1)
        write (6,*) 'length of four body file = ',ije1(nij1)
      end if
      call prt  ("tc  ","tc  ",2)
      call fmax ("tc  ",tmax,itmax)
      call fmax ("t1c ",t1max,it1max)
      if (iopt(17).ne.0) then
        call fmax6 ("t30 ",t3max,id3mx,ind3mx)
        if (iopt(17).ge.4) then
          call fmax8 ("t40 ",t4max,id4mx,ind4mx)
        end if
      end if
      write(6,2400) tmax,idec(itmax),t1max,idec(it1max)
      if (iopt(17).ne.0) then
        write (6,2401) t3max,idec(id3mx),idec(ind3mx)
        if (iopt(17).ge.4) then
          write (6,2401) t4max,idec(id4mx),idec(ind4mx)
        end if
      end if

      if (iopt(4).gt.0) then
        if (iopt(17).eq.0) then
          call dbdiag  ("hp  ","hhpp","dum3","dum4",1100,1)
        else if (iopt(17).ge.4) then
          call dbdiag  ("hp  ","hhpp","t30 ","t40 ",1111,1)
        else
          call dbdiag  ("hp  ","hhpp","t30 ","dum4",1110,1)
        end if
        call wrtrle ("t1c ",1,1)    
        call wrtrle ("tc  ",2,1)    
        if (iopt(17).ne.0) then
          call wrtrle ("t30 ",3,1)    
          if (iopt(17).ge.4) then
            call wrtrle ("t40 ",4,1)    
          end if
        end if
      endif

      call t1sqr("t1c ","t1c ","tc  ","tm  ")

      call ener(0)

!     initialization
        diff=0.0d0
       diff1=0.0d0
         df3=0.0d0
         df4=0.0d0

!  ... iteration cycle ....
220   iter=iter+1

      call prtime("start iter",iter,0)
      call flush (6)
      call ordta_csd
      call calt1
      call calt2

      !call check_modules 

      if (iopt(17).ne.0) then
        call calt3
      end if
      if (iopt(17).ge.4) then
        call calt4
      end if
      call ordtb

      if (iopt(6).ne.0.and.iter.ge.mod(iopt(6),100)) then
        fnew=max(dfloat(iopt(6)/100),1.d0)
        if (fnew.le.6.5d0) then
          fnew=1.d0/fnew
        else
          fnew=fnew-5.d0
        endif
        call add("tnew",fnew/(1.d0+fnew),"tc  ",                        &
   1.d0/(1.d0+fnew),"tnew")
        call add("t1nw",fnew/(1.d0+fnew),"t1c ",                        &
   1.d0/(1.d0+fnew),"t1nw")
        if (iopt(17).ne.0) then
          call add("t3nw",fnew/(1.d0+fnew),"t30 ",                      &
     1.d0/(1.d0+fnew),"t3nw")
        end if
        if (iopt(17).ge.4) then
          call add("t4nw",fnew/(1.d0+fnew),"t40 ",                      &
     1.d0/(1.d0+fnew),"t4nw")
        end if
        write (6,'('' t averaged with fnew '',f7.3)') fnew
      endif

      if (iopt(17).ne.0) then
        call fmax6 ("t3nw",t3max,id3mx,ind3mx)
        if (iopt(17).ge.4) then
          call fmax8 ("t4nw",t4max,id4mx,ind4mx)
        end if
      end if

      call fmax("tnew",tmax,itmax)
      call fmax("t1nw",t1max,it1max)
      if (iopt(17).ge.4) then
        call dfmax4 ("t40 ","t4nw","s1  ",df4,iddf4,inddf4)
      end if
      if (iopt(17).ne.0) then
        call dfmax3 ("t30 ","t3nw","s1  ",df3,iddf3,inddf3)
      end if
      call diffmax("tc  ","tnew","s2  ",diff,idiff)
      call diffmax("t1c ","t1nw","s2  ",diff1,idiff1)

      write (LUPRI,*) ' '
      write (LUPRI,'('' iteration # '',i3,'' diffmx='',e10.3,2x,i12,    &
  e10.3,2x,i12)')                                                       &
  iter,diff,idec(idiff),diff1,idec(idiff1)

      if (iopt(17).ne.0) then
        write (LUPRI,2402) df3,idec(iddf3),idec(inddf3)
      end if
      if (iopt(17).ge.4) then
        write (LUPRI,2402) df4,idec(iddf4),idec(inddf4)
      end if
!2401  format (24x,f14.9,2x,i8,2x,i8)
2401  format (24x,d14.9,2x,i8,2x,i8)
2402  format (22x,e10.3,2x,i8,2x,i8)

      write (LUPRI,*) ' '
      write (LUPRI,2400) tmax,idec(itmax),t1max,idec(it1max)
!2400  format(' max t, t1, t3, t4=',f15.9,2x,i12,f17.9,2x,i6)
2400  format(' max t, t1, t3, t4=',d15.9,2x,i12,d17.9,2x,i6)

      if (iopt(17).ne.0) then
        write (6,2401) t3max,idec(id3mx),idec(ind3mx)
      end if
      if (iopt(17).ge.4) then
        write (6,2401) t4max,idec(id4mx),idec(ind4mx)
      end if

!     ... error output ...
      write (LUPRI,*) ' '
      if (abs(tmax).gt.1.d10.or.abs(tmax).lt.1.d-10) then
        write (*,*) ' abort. tmax out of bounds'
        call abext('tmax    ') 
      endif

! diis
      if (iopt(4).lt.0) then
        call diis(iter,icnvg,11,"t1nw","t1c ","tnew","tc  ")
        if(iter.gt.1) then
          call diffmax("tc  ","tnew","s2  ",diff,idiff)
          call diffmax("t1c ","t1nw","s2  ",diff1,idiff1)
          write (6,*) ' '
          write (6,'('' iteration # '',i3,'' diis diffmx='',e10.3,2x,   &
      i12,e10.3,2x,i12)')                                               &
      iter,diff,idec(idiff),diff1,idec(idiff1)
        endif
        if(icnvg.eq.1) go to 910
      endif
      call move ("tnew","tc  ")
      call move ("t1nw","t1c ")
      if (iopt(17).ne.0) then
        call move ("t3nw","t30 ")
      end if
      if (iopt(17).ge.4) then
        call move ("t4nw","t40 ")
      end if
      call t1sqr ("t1c ","t1c ","tc  ","tm  ")
      call ener(0)
      if (iopt(17).ge.4) then
        call prt8 ("t4nw","new t4  ",2)
      end if
      if (iopt(17).ne.0) then
        call prt6 ("t3nw","new t3  ",2)
      end if
      call prt("tnew","new t   ",2)
      call prt("t1nw","new t1  ",2)

!     write(LUPRI,*) abs(diff),abs(diff1),abs(df3),abs(df4) ! control printout

!  ----   if end of convergence, jump out
      if (max(abs(diff),abs(diff1),abs(df3),abs(df4)).lt.convrg) go to 910

      call wrtfile (51,0,"tc  ")
      call wrtfile (51,1,"t1c ")
      if (iopt(17).ne.0) then
        call wrtfile (51,2,"t30 ")
        if (iopt(17).ge.3) then
          call wrtbox_csd
        end if
        if (iopt(17).ge.4) then
          call wrtfile (51,3,"t40 ")
        end if
      end if
      if (iopt(4).gt.0) go to 260
!   no rle
      if (iter.ge.iopt(2)) go to 920
      go to 220
!   with rle
260   continue
      itr1=mod(iter-1,iopt(4))+1
      call wrtrle ("t1c ",1,itr1+1)
      call wrtrle ("tc  ",2,itr1+1)
      if (iopt(17).ne.0) then
        call wrtrle ("t30 ",3,itr1+1)
      end if
      if (iopt(17).ge.4) then
        call wrtrle ("t40 ",4,itr1+1)
      end if
      if (itr1.lt.iopt(4)) go to 220
      if (iopt(17).eq.0) then
        call rle (itr1,1100,"t1c ","tc  ","t30 ","dum4",ising)
      else if (iopt(17).ge.4) then
        call rle (itr1,1111,"t1c ","tc  ","t30 ","t40 ",ising)
      else
        call rle (itr1,1110,"t1c ","tc  ","t30 ","dum4",ising)
      end if
!  rle singular
      if (ising.eq.1) then
        if (iopt(17).eq.0) then
          call rle (itr1-1,1100,"t1c ","tc  ","dum3","dum4",ising)
        else if (iopt(17).ge.4) then
          call rle (itr1-1,1111,"t1c ","tc  ","t30 ","t40 ",ising)
        else
          call rle (itr1-1,1110,"t1c ","tc  ","t30 ","dum4",ising)
        end if
        if (ising.eq.1) then
          iopt(4)=0
          write (*,*) ' ***** rle singular, turned off'
          go to 220
        endif
      endif
      call t1sqr ("t1c ","t1c ","tc  ","tm  ")
      call wrtrle ("t1c ",1,1)
      call wrtrle ("tc  ",2,1)
      if (iopt(17).ne.0) then
        call wrtrle ("t30 ",3,1)
      end if
      if (iopt(17).ge.4) then
        call wrtrle ("t40 ",4,1)
      end if
      if (iter.ge.iopt(2)) go to 920
      go to 220
!
!  different endings
910   write (LUPRI,'('' csd  converged after iter '',i3)') iter
      call ener(2)
      write(LUPRI,*)
      go to 970
920   write (LUPRI,*) ' iteration limit reached'
      call abext('limit   ')     

970   continue

      if (iopt(17).ge.4) then
        call prt8 ("t40 ","last t4 ",0)
      end if
      if (iopt(17).ne.0) then
        call prt6 ("t30 ","last t3 ",0)
      end if
      call prt ("tc  ","last t  ",0)
      call prt ("t1c ","last t1 ",0)
      call wrtfile (51,0,"tc  ")
      call wrtfile (51,1,"t1c ")
      if (iopt(17).ne.0) then
        call wrtfile (51,2,"t30 ")
        if (iopt(17).ge.3) then
          call wrtbox_csd
        end if
        if (iopt(17).ge.4) then
          call wrtfile (51,3,"t40 ")
        end if
      end if

   !   write (LUPRI,'(/,i8,a,i8,a)')                                     &
   ! nstrt(1,nflcur),' words used in igen,',                             &
   ! nstrt(2,nflcur),' words used in gen.'
    call space_usage

  call prtime ("end     ",0,0)

  end ! of program/routine csd

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine ordv_csd
!-------------------------------------------------------------------------------
!
! va: hhpp
! vb: hphp
! vc: hhhp
! vd: hppp
!              hhpp 1324 = va1
!              hhpp 1423 = va2
!              hhpp 3412 = va4
!              hhpp 1243 = va6
!              hhpp 3421 = va8
!              hhhp 3412 = vc4
!              hhhp 3421 = vc8
!              hhhp 3124 = vc9
!              hhhp 3214 = vc10
!              hhhp 2413 = vc11
!              hppp 2413 = vd7
!              hppp 2314 = vd8
!              hppp 2134 = vd12
!              hppp 2143 = vd14
!
!-------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
 use cint_mod
 use dat_mod
 use fock_mod
 use hsfs_parameters
implicit none
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
implicit none
#include "param_ccsd.h"
#include "cint.h"
!     common /cint/ nint,idum1,vint(ndint),iint(ndint)
#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "fock.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)
#endif

      if (iopt(1).ge.5) then
        call header('Output from ordv_csd',-1)
      endif

      call readv (nint,vint,iint,ndint,nfock,fock,ifock,ndt1)

      call prtime   ("readv   ",nint,0)

      !print *,"**************************"
      call seltyp   (2233,0,"hhpp",4)
      call seltyp   (2333,0,"hppp",4)
      call seltyp   (2323,0,"hphp",4)
      call seltyp   (2222,0,"hhhh",4)
      call seltyp   (2223,0,"hhhp",4)
       !print *,"**************************"
      call seltyp4  (3333,0,"pppp",4)
      call reorder  ("hhpp","va1 ",1324,4)
      call reorder  ("hhpp","va2 ",1423,4)
      call reorder  ("hhpp","va4 ",3412,4)
      call reorder  ("hhpp","va6 ",1243,4)
      call reorder  ("hhpp","va8 ",3421,4)
      call reorder  ("hhhp","vc4 ",3412,4)
      call reorder  ("hhhp","vc8 ",3421,4)
      call reorder  ("hhhp","vc9 ",3124,4)
      call reorder  ("hhhp","vc10",3214,4)
      call reorder  ("hhhp","vc11",2413,4)
      call reorder  ("hppp","vd7 ",2413,4)
      call reorder  ("hppp","vd8 ",2314,4)
      call reorder  ("hppp","vd12",2134,4)
      call reorder  ("hppp","vd14",2143,4)
      call selt1    (22,0,"hh  ",0)
      call selt1    (23,0,"hp  ",0)
      call selt1    (33,0,"pp  ",0)
      call selt1    (32,0,"ph  ",0)
!
      if (iopt(17).ne.0) then ! for CCSDT-x 
        call reorder  ("hppp","ve1 ",2431,4)
        call reorder  ("hppp","ve2 ",2341,4)
        call reorder  ("hhpp","ve3 ",1342,4)
        call reorder  ("hhpp","ve4 ",1432,4)
      end if   
      call prtime("ordv",0,0)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine tinit
!------------------------------------------------------------------------
!  Initialize amplitudes.
!
!   tc = hhpp / eps; t1c = hp/eps
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
  use dat_mod
implicit none
#include "priunit.h"
#else
implicit none
!      implicit real*8 (a-h,o-z)
#include "priunit.h"
#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif
      if (iopt(20).ne.0) then
       call rdfile (51,0,"tc  ",-1)
       call rdfile (51,1,"t1c ",-1)
       if ((iopt(17).ne.0).and.(iopt(20).eq.1)) call rdfile (51,2,"t30 ",-1)
      else
       call move     ("hhpp","tc  ")
       call diveps   ("tc  ")
       call move     ("hp  ","t1c ")
       call diveps   ("t1c ")
      end if  

      if ((iopt(17).ne.0).and.(iopt(20).ne.1)) then
        call clear6("t30 ",222333,0)
        if (iopt(17).ge.4) then
          call clear8("t40 ",22223333,0)
        end if
      end if
end subroutine

!**************************************************************
      subroutine ener(id)
!--------------------------------------------------------------
!
!  Calculates CCSD correlation energy
!
!   id parameter : 2 - final correlation energy
!        er than 1,2 - iteration correl.energy
!   id=1 for rle, 0 regularly.
!
!--------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
#include "implicit.h"
!     implicit real*8 (a-h,o-z)
#include "priunit.h"
#include "param_ccsd.h"
#include "cfiles.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)
#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#endif
      e=0.d0
      call add("hhpp",-2.d0,"va6 ",1.d0,"s2  ")

      if (id.eq.1) then
        call t1sqr("tr1 ","tr1 ","tr  ","s3  ")
        call files_ccsd("s3  ",ntm,new,2,-1,-1)
      else
        call files_ccsd("tm  ",ntm,new,2,-1,-1)
      endif
      call files_ccsd("s2  ",nts2,new,1,-1,-1)

      nsize=ije1(nij1)
      istrt1=nstrt(2,ntm-1)
      istrt2=nstrt(2,nts2-1)
      e=-ddot(nsize,gen(istrt1+1),1,gen(istrt2+1),1)

      if (id.eq.1) then
        call files_ccsd("tr1 ",nt1,new,1,-1,-1)
      else
        call files_ccsd("t1c ",nt1,new,1,-1,-1)
      endif
      call files_ccsd("hp  ",nth,new,2,-1,-1)
      nsize=ije1(nij1)
      istrt1=nstrt(2,nt1-1)
      istrt2=nstrt(2,nth-1)

      ! ... final formula ...
      e=e+2.d0*ddot(nsize,gen(istrt1+1),1,gen(istrt2+1),1)

      if (id.eq.1) then
        write (LUPRI,'(/,2X,A,F20.14)')                                 &
   'correlation energy (rle)=',e
      else if (id.eq.2) then
        write (LUPRI,'(/,2X,A,F20.14)')                                 &
   'correlation energy(final)=',e
      else
        write (LUPRI,'(/,2X,A,F20.14)')                                 &
   'correlation energy=',e
      endif

      return
      end subroutine ener

! *****************************************************
      subroutine ordta_csd
!---------------------------------------------------------------------
!
! does reordering 
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif

      call reorder("tc  ","ta1 ",1324,4)
      call reorder("tc  ","ta2 ",1423,4)
      call reorder("tm  ","tm1 ",1324,5)
      call reorder("tm  ","tm2 ",1423,4)
      call reorder("tc  ","ta3 ",1432,4)
      call reorder("tc  ","ta4 ",3412,4)
      call reorder("tm  ","tm4 ",3412,4)
      call reorder("t1c ","t1or",1243,4)
      ! hhpp=tc
      call reorder("tc  ","ta5 ",1342,4)
      call clear("tb1 ",2323,0)
      call clear("tb4 ",3322,0)
      call clear("tb15",2323,0)
      if (iopt(17).ne.0) then
        call clear    ("tb16",2233,0)
        call clear6   ("t3b1",223332,0)
      end if
      call prtime("ordta   ",0,1)
      return
      end
! ******************************************************
  subroutine ordtb
!---------------------------------------------------------------------
!
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
    use dat_mod
implicit none
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
implicit none
#include "priunit.h"
#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif
      call diveps    ("t1nw")
      call reorder("tb1  ","s2  ",1324,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call reorder("tb4 ","s2  ",3412,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call reorder("tb15","s2  ",1342,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      if (iopt(17).ne.0) then
        call reorder   ("tb16","s1  ",1243,4)
        call add       ("s1  ",1.d0,"tnew",1.d0,"tnew")
      end if
      call sqr("tnew","s2  ")
      call diveps("tnew")
      if (iopt(17).ge.3) then
        call reord6  ("t3b1","s1  ",126345)
        call add     ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
        call sqr3    ("t3nw")
        call dveps6  ("t3nw")
      end if
      if (iopt(17).ge.4) then
        call sqr4     ("t4nw")
        call dveps8   ("t4nw")
      end if
      call prtime("ordtb   ",0,1)
      return
 end subroutine
! **********************************************************
 subroutine calt1
!-----------------------------------------------------------
!  calculate T1 amplitudes within 9 diagrams, 
! dg1n2_csd ... dg1n11_csd
!-----------------------------------------------------------
#if defined HSFS_DYN_ALLOC
  use dat_mod
implicit none
#include "priunit.h"
#else
implicit none
!     implicit real*8 (a-h,o-z)
#include "priunit.h"
#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif
      call move("hp  ","t1nw")
      call dg1n2_csd
      call dg1n3_csd
      call dg1n5_csd
      call dg1n6_csd
      call dg1n7_csd
      call dg1n8_csd
      call dg1n9_csd
      call dg1n10_csd
      call dg1n11_csd
      call prtime    ("calt1   ",0,1)
 end subroutine
! ****************************************************
 subroutine calt2
!--------------------------------------------------------------------------
! Calculates T2 aplitudes in 17 diagrams (CCSD)
!
!--------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
   use dat_mod
implicit none
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
implicit none
#include "priunit.h"
#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif
      call cmult("hhpp",0.5d0,"tnew")
      call diag4_csd
      call diag5_csd
      call diag6_csd
      call diag7_csd
      call diag8_csd
      call diag9_csd
      call diag11_csd
      call diag12_csd
      call diag15_csd
      call diag13_csd
      call diag14_csd
      call diag16_csd
      call diag21_csd
      call diag2_csd
      call diag3_csd
      call dgt1p3_csd
      call dgt1p4_csd
      if (iopt(17).ge.4) then ! CCSDTQ case...
!  diags 40,41
        call reordx   ("t40 ","s1  ",12563478)
        call reordx   ("t40 ","s2  ",12573468)
        call add      ("s1  ",0.5d0,"s2  ",-1.d0,"s1  ")
        call add      ("hhpp",2.d0,"va6 ",-1.d0,"s2  ")
        call mult8    ("s1  ","s2  ","s3  ",2)
        call add      ("s3  ",1.d0,"tnew",1.d0,"tnew")
!   diag 42
        call reordx   ("t40 ","s1  ",14672358)
        call mult8    ("s1  ","hhpp","s2  ",2)
        call add      ("s2  ",0.5d0,"tnew",1.d0,"tnew")
      end if
      call prtime    ("calt2   ",0,1)
 end subroutine calt2
! *******************************************************
      subroutine calt3
!---------------------------------------------------------------------
! Calculates T3 aplitudes (for CCSDTx like methods)
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use dat_mod
      use cfiles_mod
implicit none
#include "priunit.h"
#include "dbg.h"
#else
!     implicit real*8 (a-h,o-z)
implicit none
#include "priunit.h"
#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "cfiles.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)
#include "dbg.h"
!     common /dbg/ iter,idebug
#endif
      if (iopt(17).ge.1) then
!  dg1n13
        call reord6    ("t30 ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125634)
        call add       ("va4 ",2.d0,"va8 ",-1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",0)
        call add       ("s3  ",1.d0,"t1nw",1.d0,"t1nw")
      end if
      if (iopt(17).ge.2) then
!  diag 36
        call add       ("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",2)
        call mult      ("t1c ","s3  ","s2  ",131)
        call add       ("s2  ",1.d0,"tb16",1.d0,"tb16")
      end if
      if (iopt(17).ge.1) then
!  dg1n14
        call reord6    ("t30 ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125634)
        call add       ("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",0)
        call add       ("s3  ",1.d0,"t1nw",1.d0,"t1nw")
!  diag 31
        call add       ("vc4 ",-2.d0,"vc8 ",1.d0,"s2  ")
        call mult6     ("s2  ","s1  ","s3  ",2)
        call add       ("s3  ",1.d0,"tnew",1.d0,"tnew")
      end if
      if (iopt(17).ge.1) then
!  diag 30
        call reord6    ("t30 ","s1  ",124563)
        call add       ("ve1 ",2.d0,"ve2 ",-1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",2)
        call add       ("s3  ",1.d0,"tnew",1.d0,"tnew")
      end if
      if (iopt(17).ge.2) then
!  diag 34
        call add       ("va1 ",2.d0,"va2 ",-1.d0,"s2  ")
        call mult      ("s2  ","t1c ","s3  ",202)
! v1  Feb 95
        call add       ("s3  ",1.d0,"hp  ",1.d0,"s3  ")

        call reorder   ("s3  ","s2  ",1243,4)
        call mult6     ("s1  ","s2  ","s3  ",0)
        call add       ("s3  ",1.d0,"tnew",1.d0,"tnew")
!  diag 37
        call add       ("ve3 ",-2.d0,"ve4 ",1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",2)
        call mult      ("s3  ","t1or","s2  ",311)
        call add       ("s2  ",1.d0,"tnew",1.d0,"tnew")
      end if
      if (iopt(17).ge.1) then
!  diag 32
        call reord6    ("t30 ","s1  ",126453)
        call mult6     ("s1  ","ve1 ","s2  ",2)
        call add       ("s2  ",-1.d0,"tb16",1.d0,"tb16")
      end if
      if (iopt(17).ge.2) then
!  diag 35
        call add       ("va1 ",2.d0,"va2 ",-1.d0,"s2  ")
        call mult      ("s2  ","t1c ","s3  ",202)
! v1 Feb 95
        call add       ("s3  ",1.d0,"hp  ",1.d0,"s3  ")

        call reorder   ("s3  ","s2  ",1243,4)
        call mult6     ("s1  ","s2  ","s3  ",0)
        call add       ("s3  ",-1.d0,"tb16",1.d0,"tb16")
!  diag 39
        call mult6     ("s1  ","ve4 ","s2  ",2)
        call mult      ("s2  ","t1or","s3  ",311)
        call add       ("s3  ",1.d0,"tnew",1.d0,"tnew")
      end if
      if (iopt(17).ge.1) then
!  diag 33
        call reord6    ("t30 ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126534)
        call mult6     ("vc8 ","s1  ","s2  ",2)
        call add       ("s2  ",1.d0,"tb16",1.d0,"tb16")
      end if
      if (iopt(17).ge.2) then
!  diag 38
        call mult6     ("s1  ","va8 ","s2  ",2)
        call mult      ("t1c ","s2  ","s3  ",131)
        call add       ("s3  ",1.d0,"tb16",1.d0,"tb16")
      end if
      if (iopt(17).lt.3) then
!  dg3n1,2
        call mult6    ("tc  ","hppp","s3  ",1)
        call reorder  ("hhhp","s1  ",2143,4)
        call mult6    ("s1  ","ta3 ","s2  ",1)
        call add      ("s2  ",-1.d0,"s3  ",1.d0,"s1  ")
        call reord6   ("s1  ","t3nw",124356)
        call sqr3     ("t3nw")
        call dveps6   ("t3nw")
      else if (iopt(17).ge.3) then
! GROUP I
        call box1     
        call reorder  ("bx1 ","s1  ",3214,4)
        call mult6    ("tc  ","s1  ","t3b1",1)
        call prtime   ("I   ",0,1)
! GROUP II
        call box2    
        call reorder  ("tc  ","s1  ",3421,4)
        call mult6    ("bx2 ","s1  ","s2  ",1)
        call add      ("s2  ",1.d0,"t3b1",1.d0,"t3b1")
        call prtime   ("II  ",0,1)
! GROUP III
        call box3     
        call reord6   ("t30 ","s1  ",124563)
        call reorder  ("bx3 ","s2  ",2143,4)
        call mult6    ("s1  ","s2  ","s3  ",3)
        call add      ("s3  ",0.5d0,"t3b1",1.d0,"t3b1")
        call prtime   ("III ",0,1)
! GROUP IV
        call box4     
        call reord6   ("t30 ","s1  ",126453)
        call reorder  ("bx4 ","s2  ",1243,4)
        call mult6    ("s1  ","s2  ","s3  ",3)
        call reord6   ("s3  ","t3nw",125436)
        call prtime   ("IV  ",0,1)
! GROUP V   
        call box5     
        call reord6   ("t30 ","s1  ",126453)
        call reorder  ("bx5 ","s2  ",1243,4)
        call mult6    ("s1  ","s2  ","s3  ",3)
        call reord6   ("s3  ","s1  ",125463)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
        call prtime   ("V   ",0,1)
! GROUP VI
        call box6     
        call reord6   ("t30 ","s1  ",124356)
        call reord7   ("s1  ","s2  ",132456)
        call reord6   ("s2  ","s1  ",125634)
        call reorder  ("bx6 ","s2  ",1243,4)
        call mult6    ("s2  ","s1  ","s3  ",3)
        call reord6   ("s3  ","s1  ",123654)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
        call prtime   ("VI  ",0,1)
! GROUP VII
!  diag 1 
        call clear6   ("s1  ",222333,0)
        call t1t2     ("t1c ","tc  ","s1  ","s2  ",123456,2)
        call t1t2     ("t1c ","tc  ","s2  ","s1  ",123456,3)
        call add      ("t30 ",1.d0,"s1  ",1.d0,"s1  ")
        call mult6    ("s1  ","pppp","s2  ",3)
        call reord6   ("s2  ","s1  ",123465)
        call reord7   ("s1  ","s3  ",132456)
        call add      ("s2  ",1.d0,"s3  ",1.d0,"s1  ")
        call add      ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
!  diag 2 
        call mult6      ("t30 ","vd14","s1  ",3)
        call mult6      ("s1  ","t1or","s2  ",5)
        call add        ("s2  ",-1.d0,"t3nw",1.d0,"t3nw")
!  diag 3
        call mult6    ("t30 ","hhpp","s1  ",3)
        call mult6    ("s1  ","tm4 ","s2  ",3)
        call add      ("s2  ",0.5d0,"t3nw",1.d0,"t3nw")
        call prtime   ("VII ",0,1)
! GROUP VIII
        call box7
        call mult6    ("t30 ","bx7 ","s1  ",5)
        call add      ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
        call prtime   ("VIII",0,1)
! GROUP IX
        call box8
        call reord6   ("t30 ","s1  ",124563)
        call mult6    ("s1  ","bx8 ","s2  ",5)
        call add      ("s2  ",0.5d0,"t3b1",1.d0,"t3b1")
        call prtime   ("IX  ",0,1)
      end if
      if (iopt(17).ge.4) then
!  dg3n1
        call reordx   ("t40 ","s1  ",12356748)
        call add      ("vd7 ",2.d0,"vd8 ",-1.d0,"s2  ")
        call mult8    ("s1  ","s2  ","s3  ",3)
        call add      ("s3  ",0.5d0,"t3nw",1.d0,"t3nw")
!  dg3n2
        call reordx   ("t40 ","s1  ",12567834)
        call add      ("vc4 ",-2.d0,"vc8 ",1.d0,"s2  ")
        call mult8    ("s1  ","s2  ","s3  ",3)
        call reord6   ("s3  ","s1  ",126345)
        call add      ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
!  dg3n3
        call reordx   ("t40 ","s1  ",12457836)
        call mult8    ("s1  ","vd7 ","s2  ",3)
        call add      ("s2  ",-1.d0,"t3nw",1.d0,"t3nw")
!  dg3n4
        call reordx   ("t40 ","s1  ",12578643)
        call mult8    ("s1  ","vc4 ","s2  ",3)
        call reord6   ("s2  ","s1  ",126345)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!  dg3n5,6
        call add      ("va1 ",2.d0,"va2 ",-1.d0,"s1  ")
        call mult     ("s1  ","t1c ","s3  ",202)
        call reordx   ("t40 ","s1  ",12356748)
        call reordx   ("t40 ","s2  ",12456738)
        call add      ("s1  ",1.d0/3.d0,"s2  ",-0.5d0,"s1  ")
        call mult8    ("s1  ","s3  ","s2  ",2)
        call add      ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  dg3n7,9
        call reordx   ("t40 ","s1  ",12456387)
        call add      ("hhpp",-2.d0,"va6 ",1.d0,"s2  ")
        call mult8    ("s1  ","s2  ","s3  ",3)
        call reordx   ("t40 ","s1  ",12457386)
        call mult8    ("s1  ","hhpp","s2  ",3)
        call add      ("s3  ",0.5d0,"s2  ",1.d0,"s1  ")
        call mult6    ("s1  ","t1or","s2  ",5)
        call add      ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  dg3n8,10
        call reordx   ("t40 ","s1  ",12568743)
        call add      ("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
        call mult8    ("s1  ","s2  ","s3  ",3)
        call reordx   ("t40 ","s1  ",12578643)
        call mult8    ("s1  ","va4 ","s2  ",3)
        call add      ("s3  ",0.5d0,"s2  ",1.d0,"s1  ")
        call mult6    ("s1  ","t1c ","s2  ",5)
        call reord6   ("s2  ","s1  ",126345)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
      end if
 call prtime   ("calt3   ",0,1) 
end subroutine

! *******************************************************
 subroutine calt4
!>---------------------------------------------------------------------
!> Calculates T4 amplitudes (for CCSDTQ)
!>--------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
  use dat_mod
implicit none
#include "priunit.h"
#include "dbg.h"
#else
#include "priunit.h"
implicit none
!     implicit real*8 (a-h,o-z)
#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dbg.h"
!     common /dbg/ iter,idebug
#endif
      if (iopt(17).eq.4) then
!  dg4n1
        call mult8    ("t30 ","hppp","s1  ",1)
        call reordx   ("s1  ","t4nw",12364578)
!  dg4n2
        call reord6   ("t30 ","s1  ",124563)
        call mult8    ("s1  ","vc11","s2  ",1)
        call reordx   ("s2  ","s1  ",12683475)
        call add      ("s1  ",-0.5d0,"t4nw",0.5d0,"t4nw")
!  dg4n3
        call mult6    ("tc  ","hphp","s1  ",1)
        call mult8    ("s1  ","ta3 ","s2  ",1)
        call reordx   ("s2  ","s1  ",12463578)
        call add      ("s1  ",-1.d0,"t4nw",1.d0,"t4nw")
!  dg4n4
        call mult6    ("ta3 ","hhhh","s1  ",1)
        call mult8    ("s1  ","ta3 ","s2  ",1)
        call reordx   ("s2  ","s1  ",14563728)
        call add      ("s1  ",0.5d0,"t4nw",1.d0,"t4nw")
!  dg4n5
        call clear8   ("s1  ",22223333,0)
        call t2sqr    ("tc  ","tc  ","s1  ","s2  ",3)
        call t2sqr    ("tc  ","tc  ","s2  ","s1  ",4)
        call mult8    ("s1  ","pppp","s2  ",4)
        call reordx   ("s2  ","s1  ",12435687)
        call add      ("s1  ",1.d0,"s2  ",1.d0,"s1  ")
        call add      ("s1  ",0.25d0,"t4nw",1.d0,"t4nw")
!  dg4n6
        call mult6    ("tc  ","va1 ","s1  ",1)
        call mult8    ("s1  ","ta3 ","s2  ",1)
        call reordx   ("s2  ","s1  ",12463758)
        call add      ("s1  ",-1.d0,"t4nw",1.d0,"t4nw")
      else if (iopt(17).eq.5) then
! GROUP I
        call mult8    ("t30 ","bx1 ","s1  ",1)
        call reordx   ("s1  ","t4nw",12364587)
! GROUP II
        call reord6   ("t30 ","s1  ",124563)
        call mult8    ("s1  ","bx2 ","s2  ",1)
        call reordx   ("s2  ","s1  ",12673485)
        call add      ("s1  ",0.5d0,"t4nw",0.5d0,"t4nw")
! GROUP III
        call reordx   ("t40 ","s1  ",12356748)
        call mult8    ("s1  ","bx3 ","s2  ",4)
        call reordx   ("s2  ","s1  ",12374568)
        call add      ("s1  ",1.d0/6.d0,"t4nw",1.d0,"t4nw")
! GROUP IV
        call reordx   ("t40 ","s1  ",12356847)
        call mult8    ("s1  ","bx4 ","s2  ",4)
        call reordx   ("s2  ","s1  ",12374568)
        call add      ("s1  ",0.5d0,"t4nw",1.d0,"t4nw")
!  GROUP V
        call reordx   ("t40 ","s1  ",12356847)
        call mult8    ("s1  ","bx5 ","s2  ",4)
        call reordx   ("s2  ","s1  ",12374586)
        call add      ("s1  ",0.5d0,"t4nw",1.d0,"t4nw")
!  GROUP VI
        call reordx   ("t40 ","s1  ",12567834)
        call mult8    ("s1  ","bx6 ","s2  ",4)
        call reordx   ("s2  ","s1  ",12783456)
        call add      ("s1  ",0.5d0,"t4nw",1.d0,"t4nw")
!  GROUP VII
!   diag 1
        call t2sqr    ("tc  ","tc  ","t40 ","s1  ",3)
        call t2sqr    ("tc  ","tc  ","s1  ","s2  ",4)
        call t1t3     ("t1c ","t30 ","s2  ","s1  ",3)
        call t1t3     ("t1c ","t30 ","s1  ","s2  ",4)
        call mult8    ("s2  ","pppp","s1  ",4)
        call reordx   ("s1  ","s2  ",12435687)
        call add      ("s1  ",1.d0,"s2  ",1.d0,"s1  ")
        call add      ("s1  ",0.25d0,"t4nw",1.d0,"t4nw")
!   diag 2
        call mult8    ("t40 ","vd14","s1  ",4)
        call mult8    ("s1  ","t1or","s2  ",6)
        call add      ("s2  ",-0.5d0,"t4nw",1.d0,"t4nw")
!   diag 3
        call mult8    ("t40 ","hhpp","s1  ",4)
        call mult8    ("s1  ","tm4 ","s2  ",4)
        call add      ("s2  ",0.25d0,"t4nw",1.d0,"t4nw")
!  GROUP VIII
        call mult8    ("t40 ","bx7 ","s1  ",6)
        call add      ("s1  ",1.d0/6.d0,"t4nw",1.d0,"t4nw")
!  GROUP IX
        call reordx   ("t40 ","s1  ",12356784)
        call mult8    ("s1  ","bx8 ","s2  ",6)
        call reordx   ("s2  ","s1  ",12384567)
        call add      ("s1  ",1.d0/6.d0,"t4nw",1.d0,"t4nw")
!  GROUP XI
        call box9
        call mult8    ("bx9 ","tc  ","s1  ",1)
        call reordx   ("s1  ","s2  ",12673485)
        call add      ("s2  ",1.d0,"t4nw",1.d0,"t4nw")
!  GROUP XII
        call box10
        call mult8    ("bx10","ta3 ","s1  ",1)
        call reordx   ("s1  ","s2  ",12563478)
        call add      ("s2  ",1.d0,"t4nw",1.d0,"t4nw")
!  GROUP XIII
        call box11
        call reord6   ("t30 ","s1  ",124536)
        call mult8    ("bx11","s1  ","s2  ",5)
        call reordx   ("s2  ","s1  ",12563478)
        call add      ("s1  ",0.5d0,"t4nw",1.d0,"t4nw")
!  GROUP XIV
!   diag 1
        call reord6   ("t30 ","s1  ",124635)
        call mult6    ("s1  ","va1 ","s2  ",3)
        call mult8    ("s2  ","s1  ","s3  ",5)
        call reordx   ("s3  ","s2  ",12563478)
        call add      ("s2  ",0.5d0,"t4nw",1.d0,"t4nw")
!   diag 2
        call mult6    ("s1  ","va2 ","s2  ",3)
        call reord6   ("t30 ","s1  ",125634)
        call mult8    ("s2  ","s1  ","s3  ",5)
        call reordx   ("s3  ","s1  ",12563847)
        call add      ("s1  ",0.5d0,"t4nw",1.d0,"t4nw")
!   diag 3
        call mult6    ("t30 ","hhpp","s3  ",3)
        call reord6   ("t30 ","s1  ",124356)
        call reord7   ("s1  ","s2  ",132456)
        call reord6   ("s2  ","s1  ",125643)
        call mult8    ("s3  ","s1  ","s2  ",5)
        call reordx   ("s2  ","s1  ",12354876)
        call add      ("s1  ",0.5d0,"t4nw",1.d0,"t4nw")
      end if
      call prtime   ("calt4   ",0,1)
      return
      end
!  ******************************************************
 subroutine diag2_csd
!
! 
! called from calt2
 implicit none
#include "priunit.h"
    call mult("tc  ","pp  ","s2  ",311)
    call add("tnew",1.d0,"s2  ",1.d0,"tnew")
    call prt("s2  ","diag2_csd ",3)
    call prtime("diag2_csd ",0,1)
 end subroutine
! *******************************************************
      subroutine diag3_csd
!
!
implicit none
#include "priunit.h"
!   incl diags 3 28
!   diag 28
      call mult("t1c ","hp  ","s2  ",111)
!   diag 3
      call add("s2  ",1.d0,"hh  ",1.d0,"s2  ")
! ----------------------------------
      call mult("s2  ","ta3 ","s3  ",131)
      call add("tnew",1.d0,"s3  ",-1.d0,"tnew")
      call prt("s3  ","-diag3_csd",3)
      call prtime("diag3_csd",0,1)
      return
      end
! **********************************************
!     subroutine diag4
      subroutine diag4_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add("va1 ",2.d0,"hphp",-1.d0,"s2  ")
      call mult("ta1 ","s2  ","s1  ",222)
      call add("tb1 ",1.d0,"s1  ",1.d0,"tb1 ")
      call prt("s1  ","diag4   ",3)
      call prtime("diag4   ",0,1)
      return
      end
! *********************************************
!     subroutine diag5
      subroutine diag5_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tm2 ","va1 ","s1  ",222)
      call add("tb1 ",1.d0,"s1  ",-1.d0,"tb1 ")
      call prt("s1  ","-diag5  ",3)
      call prtime("diag5   ",0,1)
      return
      end
! *********************************************
!     subroutine diag6
      subroutine diag6_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tm2 ","hphp","s1  ",222)
      call add("tb15",1.d0,"s1  ",-1.d0,"tb15")
      call prt("s1  ","-diag6  ",3)
      call prtime("diag6   ",0,1)
      return
      end
! **************************************************
!     subroutine diag7
      subroutine diag7_csd
! most comp. demanding diagram, symmetric
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("pppp","tm  ","s2  ",222) ! "pppp"x"tm"->"s2"
      call sqr("s2  ","s3  ")
      call add("tb4 ",1.d0,"s2  ",0.5d0,"tb4 ")
      call prt("s2  ","2*diag7 ",3)
      call prtime("diag7   ",0,1)
      return
      end
! *********************************************
!     subroutine diag8
      subroutine diag8_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("hhhh","tm4 ","s2  ",222)
      call add("tnew",1.d0,"s2  ",0.5d0,"tnew")
      call prt("s2  ","2*diag8 ",3)
      call prtime("diag8   ",0,1)
      return
      end
! ******************************************************
!     subroutine diag9
      subroutine diag9_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   includes diags 9 10 20 23
      call add("tm2 ",-1.d0,"ta1 ",1.d0,"s2  ")
      call add("va1 ",2.d0,"va2 ",-1.d0,"s3  ")
      call mult("s2  ","s3  ","s1  ",222)
!   diag 20
      call add("vc8 ",-2.d0,"vc4 ",1.d0,"s2  ")
      call mult("t1or","s2  ","s3  ",131)
      call reorder("s3  ","s2  ",2143,5)
      call add("s1  ",1.d0,"s2  ",1.d0,"s1  ")
!   diag 23
      call add("vd12",2.d0,"vd14",-1.d0,"s2  ")
      call mult("t1c ","s2  ","s3  ",131)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
! ---------------------------
      call mult("s1  ","ta1 ","s2  ",222)
      call add("tb1 ",1.d0,"s2  ",1.d0,"tb1 ")
      call prt("s2  ","diag9   ",3)
      call prtime("diag9   ",0,1)
      return
      end
! **********************************************
!     subroutine diag11
      subroutine diag11_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   includes diag 11 24
      call mult("tm2 ","va1 ","s1  ",222)
!   diag 24
      call mult("t1c ","vd12","s2  ",131)
      call add("s1  ",0.5d0,"s2  ",-1.d0,"s1  ")
! --------------------------------------
      call mult("s1  ","tm2 ","s2  ",222)
      call add("tb1 ",1.d0,"s2  ",1.d0,"tb1 ")
      call prt("s2  ","diag11  ",3)
      call prtime("diag11  ",0,1)
      return
      end
! ***********************************************
!     subroutine diag12
      subroutine diag12_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   includes diag 12 25
      call mult("tm2 ","va2 ","s1  ",222)
!   diag 25
      call mult("t1c ","vd14","s2  ",131)
      call add("s1  ",0.5d0,"s2  ",-1.d0,"s1  ")
! -----------------------------------------
      call mult("s1  ","tm2 ","s2  ",222)
      call add("tb15",1.d0,"s2  ",1.d0,"tb15")
      call prt("s2  ","diag12  ",3)
      call prtime("diag12  ",0,1)
      return
      end
! **********************************************************
!     subroutine diag13
      subroutine diag13_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   include diag 13 18
      call add("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult("tm1 ","s2  ","s1  ",113)
!   diag 18
      call add("vc9 ",-2.d0,"vc10",1.d0,"s2  ")
      call mult("s2  ","t1c ","s3  ",202)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
! ----------------------------------------
      call mult("s1  ","ta3 ","s3  ",131)
      call add("tnew",1.d0,"s3  ",1.d0,"tnew")
      call prt("s3  ","diag13   ",3)
      call prtime("diag13   ",0,1)
      return
      end
! **********************************************************
!     subroutine diag14
      subroutine diag14_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   include diag 14 19
      call add("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
      call mult("tm4 ","s2  ","s1  ",113)
!   diag 19
      call add("vd7 ",2.d0,"vd8 ",-1.d0,"s2  ")
      call mult("s2  ","t1c ","s3  ",202)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
! ----------------------------------------
      call mult("tc  ","s1  ","s3  ",311)
      call add("tnew",1.d0,"s3  ",1.d0,"tnew")
      call prt("s3  ","diag14  ",3)
      call prtime("diag14  ",0,1)
      return
      end
! ******************************************************
!     subroutine diag15
      subroutine diag15_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   includes diag 15 26
      call mult("tm  ","hhpp","s1  ",222)
!   diag 26
      call mult("t1c ","hhhp","s2  ",131)
      call add("s1  ",0.5d0,"s2  ",1.d0,"s1  ")
! ------------------------------------
      call mult("s1  ","tm4 ","s2  ",222)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","diag15  ",3)
      call prtime("diag15  ",0,1)
      return
      end
! ****************************************************
!     subroutine diag16
      subroutine diag16_csd
#include "implicit.h"
#include "priunit.h"
!     implicit real*8 (a-h,o-z)
      call mult("t1c ","hppp","s2  ",131)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","diag16  ",3)
      call prtime("diag16  ",0,1)
      return
      end
! *******************************************************
      subroutine diag21_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 21 22 27 29 17
!   diag 21
      call mult("ta2 ","vc9 ","s1  ",222)
!   diag 17
      call add("s1  ",1.d0,"vc11",-1.d0,"s1  ")
!   diag 22
      call mult("vc10","ta2 ","s2  ",222)
      call reorder("s2  ","s3  ",1432,5)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
!   diag 27
      call mult("tm  ","vd14","s2  ",222)
!     diag 29
        call mult("tc  ","hp  ","s3  ",311)
        call add("s2  ",1.d0,"s3  ",1.d0,"s2  ")
!    ---------------------------
      call reorder("s2  ","s3  ",1324,5)
      call add("s1  ",1.d0,"s3  ",-1.d0,"s1  ")
! -----------------------------------
      call mult("s1  ","t1or","s2  ",311)
      call add("tb1 ",1.d0,"s2  ",1.d0,"tb1 ")
      call prt("s2  ","diag21  ",3)
      call prtime("diag21  ",0,1)
      return
      end
! *************************************************
subroutine dg1n2_csd
! called from calt1 as the 1st step
! dg1 - 1el.diagrams, n2 - second diagram in Uzi's list,
!
! cross - antisymmetric diagram, what means 2*Coulomb - Exchange (2 diagrams)
!
! diagram: h1
!
implicit none
#include "priunit.h"
  ! 2*hppp(2413) - hppp(2314) = s1, what is spin-adapted antisymmetrised integral
  call add("vd7 ",2.d0,"vd8 ",-1.d0,"s1  ")
  ! note: "tm" means tc; s3 final diagram
  call mult("tm1 ","s1  ","s3  ",113)
  ! adding results to t1 amplitudes
  call add("t1nw",1.d0,"s3  ",1.d0,"t1nw")
  call prt("s3  ","dg1n2_csd",3)
  call prtime("dg1n2_csd",0,1)
end subroutine
! **********************************************************
subroutine dg1n3_csd
! called from calt1
implicit none
#include "priunit.h"
  call add("ta1 ",2.d0,"tm2 ",-1.d0,"s2  ")
  call mult("s2  ","hp  ","s3  ",202)
  call add("t1nw",1.d0,"s3  ",1.d0,"t1nw")
  call prt("s3  ","dg1n3_csd",3)
  call prtime("dg1n3_csd",0,1)
end subroutine
! **********************************************************
!     subroutine dg1n5
      subroutine dg1n5_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("t1c ","pp  ","s2  ",111)
      call add("t1nw",1.d0,"s2  ",1.d0,"t1nw")
      call prt("s2  ","dg1n5   ",3)
      call prtime("dg1n5   ",0,1)
      return
      end
! *****************************************************************
!     subroutine dg1n6
      subroutine dg1n6_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add("vc4 ",-2.d0,"vc8 ",1.d0,"s2  ")
      call mult("s2  ","tm4 ","s3  ",113)
      call add("t1nw",1.d0,"s3  ",1.d0,"t1nw")
      call prt("s3  ","dg1n6   ",3)
      call prtime("dg1n6   ",0,1)
      return
      end
! ***************************************************************
!     subroutine dg1n7
      subroutine dg1n7_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add("va1 ",2.d0,"hphp",-1.d0,"s2  ")
      call mult("s2  ","t1c ","s3  ",202)
      call add("t1nw",1.d0,"s3  ",1.d0,"t1nw")
      call prt("s3  ","dg1n7   ",3)
      call prtime("dg1n7   ",0,1)
      return
      end
! ***********************************************************************
!     subroutine dg1n8
      subroutine dg1n8_csd
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("hh  ","t1or","s2  ",111)
      call add("t1nw",1.d0,"s2  ",-1.d0,"t1nw")
      call prt("s2  ","-dg1n8  ",3)
      call prtime("dg1n8   ",0,1)
      return
      end
! ******************************************************************
      subroutine dg1n9_csd
implicit none
#include "priunit.h"
      call add("hhpp",-2.d0,"va6 ",1.d0,"s2  ")
      call mult("tc  ","s2  ","s3  ",113)
      call mult("s3  ","t1or","s2  ",111)
      call add("t1nw",1.d0,"s2  ",1.d0,"t1nw")
      call prt("s2  ","dg1n9_csd",3)
      call prtime("dg1n9_csd",0,1)
      return
      end
! ****************************************************************
      subroutine dg1n10_csd
implicit none
#include "priunit.h"
      call add("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
      call mult("tm4 ","s2  ","s3  ",113)
      call mult("t1c ","s3  ","s2  ",111)
      call add("t1nw",1.d0,"s2  ",1.d0,"t1nw")
      call prt("s2  ","dg1n10  ",3)
      call prtime("dg1n10  ",0,1)
      return
      end
! ***************************************************************
      subroutine dg1n11_csd
implicit none
#include "priunit.h"
      call add("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult("s2  ","t1c ","s3  ",202)
      call add("ta1 ",-2.d0,"ta2 ",1.d0,"s1  ")
      call mult("s1  ","s3  ","s2  ",202)
      call add("t1nw",1.d0,"s2  ",1.d0,"t1nw")
      call prt("s2  ","dg1n11  ",3)
      call prtime("dg1n11  ",0,1)
      return
      end
! *****************************************************
      subroutine dgt1p3_csd
implicit none
#include "priunit.h"
!   compensates for overcounting t1**3
      call mult("t1c ","vd14","s2  ",131)
      call mult("t1c ","s2  ","s3  ",131)
      call mult("s3  ","t1or","s2  ",311)
      call add("tnew",1.d0,"s2  ",2.d0,"tnew")
      call prt("s2  ","dgt1p3/2",3)
      call prtime("dgt1p3  ",0,1)
      return
      end
! ****************************************
      subroutine dgt1p4_csd
implicit none
#include "priunit.h"
!   compensates for overcounting t1**4
      call mult("t1c ","hhpp","s3  ",131)
      call mult("t1c ","s3  ","s2  ",131)
      call mult("t1or","s2  ","s3  ",131)
      call mult("t1or","s3  ","s2  ",131)
      call add("tb4 ",1.d0,"s2  ",-1.d0,"tb4 ")
      call prt("s2  ","-dgt1p4 ",3)
      call prtime("dgt1p4  ",0,1)
      return
      end
!  ********************************************************************
      subroutine box1 
implicit none
#include "priunit.h"
!  diag 1
      call reorder  ("hppp","bx1 ",1342,4)
! v1 Feb 95
      call mult     ("ta5 ","ph  ","s1  ",311)
      call add      ("bx1 ",1.d0,"s1  ",-0.5d0,"bx1 ")
!  diag 2 
      call add      ("vd7 ",2.d0,"vd8 ",-1.d0,"s1  ")
      call mult     ("ta1 ","s1  ","s2  ",222)
      call add      ("s2  ",1.d0,"bx1 ",1.d0,"bx1 ")
!  diag 3 
      call mult     ("tm2 ","vd7 ","s1  ",222)
      call add      ("s1  ",-1.d0,"bx1 ",1.d0,"bx1 ")
!  diag 4 
      call mult     ("tm2 ","vd8 ","s1  ",222)
      call reorder  ("s1  ","s2  ",1324,4)
      call add      ("s2  ",-1.d0,"bx1 ",1.d0,"bx1 ")
!  diag 5 
      call mult     ("vc4 ","tm4 ","s1  ",222)
      call reorder  ("s1  ","s2  ",1342,4)
      call add      ("s2  ",1.d0,"bx1 ",1.d0,"bx1 ")
!  diag 6
      call reord6   ("t30 ","s1  ",124356)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",125634)
      call cmult    ("t1c ",0.5d0,"s2  ")
      call t1t2     ("s2  ","tc  ","s1  ","s3  ",145623,3)
      call t1t2     ("t1c ","tc  ","s3  ","s1  ",145623,2)
      call add      ("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","s3  ",2)
      call add      ("s3  ",1.d0,"bx1 ",1.d0,"bx1 ")
!  diag 7
      call reord6   ("t30 ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",126534)
      call t1t2     ("t1c ","tc  ","s1  ","s2  ",156423,1)
      call t1t2     ("t1c ","tc  ","s2  ","s1  ",156423,2)
      call t1t2     ("t1c ","tm  ","s1  ","s2  ",156423,3)
      call mult6    ("s2  ","va8 ","s1  ",2)
      call add      ("s1  ",1.d0,"bx1 ",1.d0,"bx1 ")
!  diag 8
      call mult     ("ve3 ","t1or","s1  ",311)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",-1.d0,"bx1 ",1.d0,"bx1 ")
!  diag 9
      call reorder  ("hphp","s1  ",1423,4)
      call mult     ("s1  ","t1or","s2  ",311)
      call reorder  ("s2  ","s1  ",1423,4)
      call add      ("s1  ",-1.d0,"bx1 ",1.d0,"bx1 ")
      return
      end 
!  *******************************************************************
      subroutine box2 
implicit none
#include "priunit.h"
!  diag 1
      call reorder  ("hhhp","s1  ",2143,4)
      call cmult    ("s1  ",-1.d0,"bx2 ")
! v1 Feb 95
      call mult     ("tc  ","hp  ","s1  ",311)
      call add      ("bx2 ",1.d0,"s1  ",-.5d0,"bx2 ")
!  diag 2,3   
      call add      ("vc9 ",-2.d0,"vc10",1.d0,"s1  ")
      call mult     ("ta1 ","s1  ","s2  ",222)
      call mult     ("tm2 ","vc9 ","s1  ",222)
      call add      ("s2  ",1.d0,"s1  ",1.d0,"s1  ")
      call reorder  ("s1  ","s2  ",1324,4)
      call add      ("s2  ",1.d0,"bx2 ",1.d0,"bx2 ")
!  diag 4 
      call mult     ("vc10","tm2 ","s1  ",222)
      call reorder  ("s1  ","s2  ",1342,4)
      call add      ("s2  ",1.d0,"bx2 ",1.d0,"bx2 ")
!  diag 5 
      call mult     ("tm  ","vd14","s1  ",222)
      call add      ("s1  ",-1.d0,"bx2 ",1.d0,"bx2 ")
!  diag 6 
      call reord6   ("t30 ","s1  ",124563)
      call cmult    ("t1c ",0.5d0,"s2  ")
      call t1t2     ("s2  ","tc  ","s1  ","s3  ",124563,3)
      call t1t2     ("t1c ","tc  ","s3  ","s2  ",124563,2)
      call add      ("ve3 ",-2.d0,"ve4 ",1.d0,"s1  ")
      call mult6    ("s2  ","s1  ","s3  ",2)
      call add      ("s3  ",1.d0,"bx2 ",1.d0,"bx2 ")
!  diag 7 
      call reord6   ("t30 ","s1  ",126453)
      call t1t2     ("t1c ","tc  ","s1  ","s2  ",126453,1)
      call t1t2     ("t1c ","tc  ","s2  ","s1  ",126453,2)
      call t1t2     ("t1c ","tm  ","s1  ","s2  ",126453,3)
      call mult6    ("s2  ","ve4 ","s1  ",2)
      call add      ("s1  ",1.d0,"bx2 ",1.d0,"bx2 ")
!  diag 8 
      call mult     ("va1 ","t1c ","s1  ",311)
      call reorder  ("s1  ","s2  ",1423,4)
      call add      ("s2  ",-1.d0,"bx2 ",1.d0,"bx2 ")
!  diag 9 
      call mult     ("t1c ","hphp","s1  ",131)
      call add      ("s1  ",-1.d0,"bx2 ",1.d0,"bx2 ")
!  diag 10
      call mult     ("hhhh","t1or","s1  ",311)
      call reorder  ("s1  ","s2  ",2143,4)
      call add      ("s2  ",1.d0,"bx2 ",1.d0,"bx2 ")
      return
      end 
!  *******************************************************************
      subroutine box3 
implicit none
#include "priunit.h"
!  diag 1 
      call add      ("va1 ",2.d0,"hphp",-1.d0,"bx3 ")
!  diag 2
      call add      ("vd8 ",2.d0,"vd7 ",-1.d0,"s1  ")
      call mult     ("t1c ","s1  ","s2  ",131)
      call reorder  ("s2  ","s1  ",1243,4)
      call add      ("s1  ",1.d0,"bx3 ",1.d0,"bx3 ")
!  diag 3 
      call add      ("vc8 ",-2.d0,"vc4 ",1.d0,"s1  ")
      call mult     ("s1  ","t1or","s2  ",311)
      call reorder  ("s2  ","s1  ",1432,4)
      call add      ("s1  ",1.d0,"bx3 ",1.d0,"bx3 ")
!  diag 4 
      call add      ("va1 ",2.d0,"va2 ",-1.d0,"s1  ")
      call mult     ("ta1 ","s1  ","s2  ",222)
      call add      ("s2  ",2.d0,"bx3 ",1.d0,"bx3 ")
!  diag 5 
      call add      ("va1 ",-2.d0,"va2 ",1.d0,"s1  ")
      call mult     ("tm2 ","s1  ","s2  ",222)
      call add      ("s2  ",1.d0,"bx3 ",1.d0,"bx3 ")
      return
      end
!  *******************************************************************
      subroutine box4 
implicit none
#include "priunit.h"
!  diag 1 
      call cmult    ("va1 ",-1.d0,"bx4 ")
!  diag 2
      call mult     ("t1c ","vd8 ","s1  ",131)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",-1.d0,"bx4 ",1.d0,"bx4 ")
!  diag 3 
      call mult     ("vc8 ","t1or","s1  ",311)
      call reorder  ("s1  ","s2  ",1432,4)
      call add      ("s2  ",1.d0,"bx4 ",1.d0,"bx4 ")
!  diag 4 
      call add      ("va1 ",-2.d0,"va2 ",1.d0,"s1  ")
      call mult     ("ta1 ","s1  ","s2  ",222)
      call add      ("s2  ",1.d0,"bx4 ",1.d0,"bx4 ")
!  diag 5 
      call mult     ("tm2 ","va1 ","s1  ",222)
      call add      ("s1  ",1.d0,"bx4 ",1.d0,"bx4 ")
      return
      end 
!  *******************************************************************
      subroutine box5 
implicit none
#include "priunit.h"
!  diag 1 
      call cmult    ("hphp",-1.d0,"bx5 ")
!  diag 2 
      call mult     ("t1c ","vd7 ","s1  ",131)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",-1.d0,"bx5 ",1.d0,"bx5 ")
!  diag 3
      call mult     ("vc4 ","t1or","s1  ",311)
      call reorder  ("s1  ","s2  ",1432,4)
      call add      ("s2  ",1.d0,"bx5 ",1.d0,"bx5 ")
!  diag 4 
      call mult     ("tm2 ","va2 ","s1  ",222)
      call add      ("s1  ",1.d0,"bx5 ",1.d0,"bx5 ")
      return
      end
!  *******************************************************************
      subroutine box6 
implicit none
#include "priunit.h"
!  diag 1 
      call cmult    ("hhhh",0.5d0,"bx6 ")
!  diag 2 
      call mult     ("t1c ","vc10","s1  ",131)
      call add      ("s1  ",1.d0,"bx6 ",1.d0,"bx6 ")
!  diag 3
      call mult     ("tm  ","hhpp","s1  ",222)
      call add      ("s1  ",0.5d0,"bx6 ",1.d0,"bx6 ")
      return
      end 
!  *******************************************************************
      subroutine box7
implicit none
#include "priunit.h"
!  diag 1 
      call add      ("vd7 ",2.d0,"vd8 ",-1.d0,"s1  ")
      call mult     ("s1  ","t1c ","bx7 ",202)
! v1 Feb 95
      call add      ("bx7 ",1.d0,"pp  ",1.d0,"bx7 ")
      call mult     ("t1or","ph  ","s1  ",111)
      call add      ("bx7 ",1.d0,"s1  ",-1.d0,"bx7 ")
!  diag 2 
      call add      ("va4 ",-2.d0,"va8 ",1.d0,"s1  ")
      call mult     ("tm4 ","s1  ","s2  ",113)
      call add      ("s2  ",1.d0,"bx7 ",1.d0,"bx7 ")
      return 
      end
!  *******************************************************************
      subroutine box8
implicit none
#include "priunit.h"
!  diag 1 
      call add      ("vc9 ",-2.d0,"vc10",1.d0,"s1  ")
      call mult     ("s1  ","t1c ","bx8 ",202)
! v1 Feb2
      call add      ("bx8 ",1.d0,"hh  ",1.d0,"bx8 ")
      call mult     ("t1c ","hp  ","s1  ",111)
      call add      ("bx8 ",1.d0,"s1  ",1.d0,"bx8 ")
!  diag 2 
      call add      ("hhpp",-2.d0,"va6 ",1.d0,"s1  ")
      call mult     ("tm  ","s1  ","s2  ",113)
      call add      ("s2  ",1.d0,"bx8 ",1.d0,"bx8 ")
      return
      end 
!  *******************************************************************
      subroutine box9
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  diag 1
      call reord6   ("t30 ","s1  ",124536)
      call add      ("vd7 ",2.d0,"vd8 ",-1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","bx9 ",3)
!  diag 2
      call t1t2     ("t1c ","tc  ","t30 ","s1  ",123456,2)
      call t1t2     ("t1c ","tc  ","s1  ","s2  ",123456,3)
      call reord6   ("s2  ","s1  ",124635)
      call mult6    ("s1  ","vd7 ","s2  ",3)
      call add      ("s2  ",-1.d0,"bx9 ",0.5d0,"bx9 ")
!  diag 3
      call t1t2     ("t1c ","tc  ","t30 ","s1  ",123456,2)
      call reord6   ("s1  ","s2  ",124635)
      call mult6    ("s2  ","vd8 ","s1  ",3)
      call reord6   ("s1  ","s2  ",123546)
      call add      ("s2  ",-1.d0,"bx9 ",1.d0,"bx9 ")
!  diag 4
      call reord6   ("t30 ","s1  ",124356)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",125634)
      call mult6    ("s1  ","vc4 ","s2  ",3)
      call reord6   ("s2  ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("s2  ",1.d0,"bx9 ",1.d0,"bx9 ")
!  diag 5
      call cmult    ("t1c ",0.5d0,"s1  ")
      call t1t3     ("s1  ","t30 ","t40 ","s2  ",3)
      call reordx   ("s2  ","s1  ",12568743)
      call add      ("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
      call mult8    ("s1  ","s2  ","s3  ",3)
      call add      ("s3  ",0.5d0,"bx9 ",1.d0,"bx9 ")
!  diag 6
      call cmult    ("tc  ",0.5d0,"s1  ")
      call t2sqr    ("s1  ","tm  ","t40 ","s2  ",2)
      call t1t3     ("t1c ","t30 ","s2  ","s1  ",2)
      call t1t3     ("t1c ","t30 ","s1  ","s2  ",3)
      call t1t3     ("t1c ","t30 ","s2  ","s1  ",4)
      call reordx   ("s1  ","s2  ",12578643)
      call mult8    ("s2  ","va4 ","s1  ",3)
      call add      ("s1  ",1.d0,"bx9 ",1.d0,"bx9 ")
      return
      end 
!  *********************************************************************
      subroutine box10
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  diag 1
      call reord6   ("t30 ","s1  ",124536)
      call add      ("vc9 ",-2.d0,"vc10",1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","bx10",3)
!  diag 2
      call t1t2     ("t1c ","tc  ","t30 ","s1  ",123456,2)
      call t1t2     ("t1c ","tc  ","s1  ","s2  ",123456,3)
      call reord6   ("s2  ","s1  ",124635)
      call mult6    ("s1  ","vc9 ","s2  ",3)
      call add      ("s2  ",1.d0,"bx10",0.5d0,"bx10")
!  diag 3
      call t1t2     ("t1c ","tc  ","t30 ","s1  ",123456,2)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124536)
      call mult6    ("s1  ","vc10","s2  ",3)
      call reord6   ("s2  ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124536)
      call add      ("s1  ",1.d0,"bx10",1.d0,"bx10")
!  diag 4
      call reord6   ("t30 ","s1  ",124365)
      call mult6    ("s1  ","vd12","s2  ",3)
      call reord6   ("s2  ","s1  ",123546)
      call add      ("s1  ",-1.d0,"bx10",1.d0,"bx10")
!  diag 5
      call t2sqr    ("tc  ","tc  ","t40 ","s1  ",3)
      call cmult    ("t1c ",0.5d0,"s2  ")
      call t1t3     ("s2  ","t30 ","s1  ","s3  ",3)
      call reordx   ("s3  ","s1  ",12564387)
      call add      ("hhpp",-2.d0,"va6 ",1.d0,"s2  ")
      call mult8    ("s1  ","s2  ","s3  ",3)
      call add      ("s3  ",0.5d0,"bx10",1.d0,"bx10")
!  diag 6
      call cmult    ("tc  ",0.5d0,"s1  ")
      call t2sqr    ("s1  ","tm  ","t40 ","s2  ",3)
      call t2sqr    ("tc  ","tm  ","s2  ","s1  ",2)
      call t2sqr    ("tc  ","tm  ","s1  ","s2  ",4)
      call t1t3     ("t1c ","t30 ","s2  ","s1  ",2)
      call t1t3     ("t1c ","t30 ","s1  ","s2  ",3)
      call t1t3     ("t1c ","t30 ","s2  ","s1  ",4)
      call reordx   ("s1  ","s2  ",12574386)
      call mult8    ("s2  ","hhpp","s1  ",3)
      call add      ("s1  ",1.d0,"bx10",1.d0,"bx10")
!  diag 7
      call mult6    ("tc  ","hhpp","s1  ",1)
      call reord6   ("s1  ","s2  ",124356)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",124635)
      call add      ("s2  ",-1.d0,"bx10",1.d0,"bx10")
!  diag 8
      call mult6    ("tc  ","hphp","s1  ",1)
      call reord6   ("s1  ","s2  ",123564)
      call add      ("s2  ",-1.d0,"bx10",1.d0,"bx10")
!  diag 9
      call mult6    ("ta3 ","hhhh","s1  ",1)
      call reord6   ("s1  ","s2  ",125346)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",124356)
      call add      ("s2  ",0.5d0,"bx10",1.d0,"bx10")
      return
      end 
!  *********************************************************************
      subroutine box11
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  diag 1
      call reord6   ("t30 ","s1  ",124536)
      call add      ("va1 ",2.d0,"va2 ",-1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","bx11",3)
!  diag 2
      call t1t2     ("t1c ","tc  ","t30 ","s1  ",123456,2)
      call t1t2     ("t1c ","tc  ","s1  ","s2  ",123456,3)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",124536)
      call add      ("va1 ",-2.d0,"va2 ",1.d0,"s1  ")
      call mult6    ("s2  ","s1  ","s3  ",3)
      call add      ("s3  ",1.d0,"bx11",0.5d0,"bx11")
      return
      end
!  *********************************************************************
!     subroutine wrtbox
      subroutine wrtbox_csd
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "dat.h"
!      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif
      call wrtfile (20,0,"bx1 ")
      call wrtfile (20,1,"bx2 ")
      call wrtfile (20,2,"bx3 ")
      call wrtfile (20,3,"bx4 ")
      call wrtfile (20,4,"bx5 ")
      call wrtfile (20,5,"bx6 ")
      call wrtfile (20,6,"bx7 ")
      call wrtfile (20,7,"bx8 ")
      if (iopt(17).eq.5) then
        call wrtfile (20,8,"bx9 ")
        call wrtfile (20,9,"bx10")
        call wrtfile (20,10,"bx11")
      end if
      return
      end 
!  *******************************************************************
