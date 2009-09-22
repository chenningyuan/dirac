!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine chu 
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use gen_mod
      use cij1_mod
      use cfiles_mod
      use cint_mod
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#include "csdtime.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)
#include "cfiles.h"
!     common /cint/ nint,idum1,vint(ndint),iint(ndint)
#include "cint.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /csdtime/ csdstim,csdltim
#include "csdtime.h"
#endif

      character*256 filenm,filenm51,filenm53,filenm55,filenm6,
     1              filenm20,filenm31,filenm32
c     iopt(1) - print.
c     iopt(2) - no. of iterations.
c     iopt(3) - convergence = 10**(-iopt(3))
c     iopt(4) .gt.0 - do rle every iopt(4) iterations.
c             .lt.0 - diis
c     iopt(5) - print only elements gt 10**(-iopt(5))
c     iopt(6) - start avrg.
c     iopt(8) - fock.
c     iopt(10) .ne.0 - do not save pppp.
      real*4 cputime
      filenm(1:1)=' '
      call getenv('input',filenm)
      if (filenm(1:1).eq.' ') filenm='input'
      open(5,file=filenm)
      filenm6(1:1)=' '
      call getenv('output',filenm6)
      if (filenm6(1:1).ne.' ') open(6,file=filenm6,form='FORMATTED',
     1 status='UNKNOWN',access='SEQUENTIAL')
      filenm51(1:1)=' '
      call getenv('tape51',filenm51)
      if (filenm51(1:1).eq.' ') filenm51='tape51'
      open(51,file=filenm51,form='UNFORMATTED')
      filenm53(1:1)=' '
      call getenv('tape53',filenm53)
      if (filenm53(1:1).eq.' ') filenm53='tape53'
      open(53,file=filenm53,form='UNFORMATTED')
      filenm55(1:1)=' '
      call getenv('tape55',filenm55)
      if (filenm55(1:1).eq.' ') filenm55='tape55'
      open(55,file=filenm55,form='UNFORMATTED')
      filenm20(1:1)=' '
      call getenv('tape20',filenm20)
      if (filenm20(1:1).eq.' ') filenm20='tape20'
      open(20,file=filenm20,form='UNFORMATTED')
      filenm31(1:1)=' '
      call getenv('tape31',filenm31)
      if (filenm31(1:1).eq.' ') filenm31='tape31'
      open(31,file=filenm31,form='UNFORMATTED')
      filenm32(1:1)=' '
      call getenv('tape32',filenm32)
      if (filenm32(1:1).eq.' ') filenm32='tape32'
      open(32,file=filenm32,form='UNFORMATTED')
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')
      write (*,*) '   $$$$$$$$$   chu   $$$$$$$$$  '
      csdstim = cputime(real(0.))
      csdltim = csdstim
      call rdin
      call rdsym
      if (iopt(10).eq.0) then
        filenm(1:1)=' '
        call getenv('tape9',filenm)
        if (filenm(1:1).eq.' ') filenm='tape9'
        open(9,file=filenm,form='UNFORMATTED')
      endif
      convrg = 10.**(-iopt(3))
      write (6,'('' max # of iter='',i3,''. convrg='',e10.2)') 
     1 iopt(2),convrg
      call prtime("rdin    ",0,0)
      if (iopt(9).ne.0) call rdrsp2
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      irecu=5
      call reads   (51,0,"s1  ")
      call reorder ("s1  ","tc  ",3412,4)
      call reads   (51,1,"s1  ")
      call reorder ("s1  ","t1c ",1243,4)
      call rdfile  (53,0,"thv ",-1)
      call rdfile  (53,1,"t1hv",-1)
      call rdfile  (53,2,"thvf",-1)
      if (iopt(17).ne.0) then
       ! call rdfile  (53,3,"t31 ")
        call rdfile  (53,3,"t31 ",-1)
        if (iopt(17).eq.3) then
          call reads  (51,2,"s1  ")        
          call reord6 ("s1  ","s2  ",124356)
          call reord7 ("s2  ","s1  ",312456)
          call reord6 ("s1  ","s2  ",125346)
          call reord7 ("s2  ","s1  ",132456)
          call reord6 ("s1  ","t30 ",126345)
!         call rdbox
          call rdbox_chu
        end if
      end if
      call t1sqr("t1c ","t1c ","tc  ","tcm ")
      call t1sqr("t1hv","t1c ","thv ","tvm ")
!     call ordvu
      call ordvu_chu
      call calvu0
c      call pdiag
c  pdiag eliminated for lack of io file numbers. to be reinstated
c  if possible, as significant time saving is achieved.
      iter=0
      ieof=0
!     call tuinit
      call tuinit_chu
      call prt("thu ","thu     ",2)
      call prt("thuf","thuf    ",2)
      call fmax("thu ",tmax,itmax)
      write(6,2400) tmax,idec(itmax)
      if (iopt(4).gt.0) then
        if (iopt(17).eq.0) then
          call dbdiag ("dum1","vu0 ","dum3","dum4",0100,-1)
        else
          call reord6 ("v30 ","s1  ",123546)
          call reord7 ("s1  ","s2  ",213456)
          call add    ("v30 ",1.d0,"s2  ",1.d0,"s1  ")
          call dbdiag ("dum1","vu0 ","s1  ","dum4",0110,-1)
        end if
        call wrtrle ("thu ",2,1)    
        if (iopt(17).ne.0) then
          call wrtrle ("t32 ",3,1)    
        end if
      endif
      call t1sqr("t1hv","t1hv","thu ","tm  ")
220   iter=iter+1
      call prtime("strt itr",iter,0)
      call ordtua_chu
      if (iopt(17).ne.0) then
        call caltu3_chu
      end if
      call caltu2_chu
      call ordtub_chu
      if (iopt(6).ne.0.and.iter.ge.mod(iopt(6),100)) then
        fnew=max(dfloat(iopt(6)/100),1.d0)
        if (fnew.le.6.5d0) then
          fnew=1.d0/fnew
        else
          fnew=fnew-5.d0
        endif
        call add("tnew",fnew/(1.d0+fnew),"thu ",
     1  1.d0/(1.d0+fnew),"tnew")
        if (iopt(17).ne.0) then
          call add("t3nw",fnew/(1.d0+fnew),"t32 ",
     1    1.d0/(1.d0+fnew),"t3nw")
        end if
        write (6,'('' t averaged with fnew '',f7.3)') fnew
      endif
      call fmax("tnew",tmax,itmax)
      t3max = 0.0d0 ! mi init.
      if (iopt(17).ne.0) then
        call fmax6 ("t3nw",t3max,id3mx,ind3mx)
      end if
      call diffmax("thu ","tnew","s2  ",diff,idiff)
      df3 = 0.0d0 ! mi init
      if (iopt(17).ne.0) then
        call dfmax3 ("t32 ","t3nw","s2  ",df3,iddf3,inddf3)
      end if
      write (6,*) ' '
      write (6,'('' iteration # '',i3,'' diffmx='',
     1 e10.3,2x,i12)') iter,diff,idec(idiff)
      if (iopt(17).ne.0) then
        write (6,2401) df3,idec(iddf3),idec(inddf3)
2401    format (24x,e10.3,2x,i12,2x,i12)
      end if
      write(6,2400) tmax,idec(itmax),t3max,idec(id3mx),idec(ind3mx)
2400  format(' max t,t3=',f15.9,2x,i12,2x,f15.9,2x,i6,2x,i8)
      if (abs(tmax).gt.1.d10.or.abs(tmax).lt.1.d-10) then
        write (*,*) ' abort. tmax out of bounds'
        call abext('tmax    ')      
      endif
c diis
      if (iopt(4).lt.0) then
        call diis(iter,icnvg,01,"dumm","dumm","tnew","thu ")
        if(iter.gt.1) then
          call diffmax("thu ","tnew","s2  ",diff,idiff)
          write (6,*) ' '
          write (6,'('' iteration # '',i3,'' diis diffmx='',e10.3,2x,
     1     i12)')  iter,diff,idec(idiff)
        endif
        if(icnvg.eq.1) go to 910
      endif
      call move("tnew","thu ")
      if (iopt(17).ne.0) then
        call move ("t3nw","t32 ")     
      end if
      call t1sqr("t1hv","t1hv","thu ","tm  ")
      call prt("tnew","new t   ",2)
      if (max(abs(diff),abs(df3)).lt.convrg) go to 910
      if (iopt(4).gt.0) go to 260
c   no rle
      if (iter.ge.iopt(2)) go to 920
      go to 220
c   with rle
260   continue
      itr1=mod(iter-1,iopt(4))+1
      call wrtrle ("thu ",2,itr1+1)
      if (iopt(17).ne.0) then
        call wrtrle ("t32 ",3,itr1+1)
      end if
      if (itr1.lt.iopt(4)) go to 220
      if (iopt(17).eq.0) then
        call rle  (itr1,0100,"dum1","thu ","dum3","dum4",ising)
      else
        call rle  (itr1,0110,"dum1","thu ","t32 ","dum4",ising)
      end if
c  rle singular
      if (ising.eq.1) then
        if (iopt(17).eq.0) then
          call rle  (itr1-1,0100,"dum1","thu ","dum3","dum4",ising)
        else
          call rle  (itr1-1,0110,"dum1","thu ","t32 ","dum4",ising)
        end if
        if (ising.eq.1) then
          iopt(4)=0
          write(*,*) ' ***** rle singular, turned off'
          go to 220
        endif
      endif
      call t1sqr("t1hv","t1hv","thu ","tm  ")
      call wrtrle ("thu ",2,1)
      if (iopt(17).ne.0) then
        call wrtrle ("t32 ",3,1)
      end if
      if (iter.ge.iopt(2)) go to 920
      go to 220
c
c  different endings
910   write (6,*) ' chu  converged after iter ',iter
      egr0=0.d0
      call eneru2 (egr0,1,2,max(1,iopt(1)))
      call eneru2 (egr0,2,2,max(1,iopt(1)))
      go to 970
920   write (6,*) ' iteration limit'
      call abext('limit   ')     
970   call prt("thu ","last t  ",0)
      call prt("thuf","thuf",0)
      call wrtfile(55,0,"thu ")
      call wrtfile(55,1,"thuf")
      if (iopt(17).ne.0) then
        call wrtfile  (55,2,"t32 ")
        if (iopt(17).eq.3) then
!         call wrtbox
          call wrtbox_chu
        end if
      end if
      write (*,*) nstrt(1,nflcur),' words used in igen,',
     1  nstrt(2,nflcur),' words used in gen'
      call prtime("end     ",iter,0)
!      stop
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine ordvu_chu
!------------------------------------------------------------------------
c va: pphh
c vd: phhh
c vva: gphh
c vvc: gpph
c              pphh 1324 = va1
c              pphh 1423 = va2
c              pphh 3412 = va4
c              pphh 1243 = va6
c              pphh 3421 = va8
c              phhh 2314 = vd7
c              phhh 2413 = vd11
c              phhh 2134 = vd12
c              phhh 2143 = vd14
c              gphh 1324 = vva1
c              gpph 1324 = vvc1
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cint_mod
      use fock_mod
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!     common /cint/nint,idum1,vint(ndint),iint(ndint)
#include "cint.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)
#include "fock.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      external idec
c
      call readv    (nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call prtime   ("readv   ",nint,0)
      call seltyp   (3322,0,"pphh",4)
      call seltyp   (3222,0,"phhh",4)
      call seltyp   (2222,0,"hhhh",4)
      call seltyp   (2322,1000,"gphh",4)
      call seltyp   (2232,1000,"ghph",4)
      call seltyp   (2332,1000,"gpph",4)
      call seltyp   (2222,1100,"gghh",4)
      call seltyp   (2233,1100,"ggpp",4)
      call reorder  ("pphh","va1 ",1324,4)
      call reorder  ("pphh","va2 ",1423,4)
      call reorder  ("pphh","va4 ",3412,4)
      call reorder  ("pphh","va6 ",1243,4)
      call reorder  ("pphh","va8 ",3421,4)
      call reorder  ("phhh","vd7 ",2314,4)
      call reorder  ("phhh","vd11",2413,4)
      call reorder  ("phhh","vd12",2134,4)
      call reorder  ("phhh","vd14",2143,4)
      call reorder  ("gphh","vva1",1324,4)
      call reorder  ("gpph","vvc1",1324,4)
      call selt1    (22,0,"hh  ",4)
      call selt1    (32,0,"ph  ",4)
      call selt1    (23,10,"gp  ",4)
      if (iopt(17).ne.0) then
        call seltyp   (2222,1000,"ghhh",4)
        call seltyp   (2323,1000,"gphp",4)
        call seltyp   (2223,1100,"gghp",4)
        call seltyp   (2222,1111,"gggg",4)
      end if
      call prtime("ordu    ",0,0)
      return
      end
c**************************************************************
      subroutine calvu0
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      call ordta0_chu
      call move("gghh","vu0 ")
      call cmult("vu0 ",0.5d0,"vu0 ")
      call diag3_chu
      call diag4_chu
      call diag5_chu
      call diag6_chu
      call diag8_chu
      call diag9_chu
      call diag11_chu
      call diag12_chu
      call diag13_chu
      call dgt1p3_chu
      call dgt1p4_chu
      call diag21_chu
      call diag26_chu
      call diag28_chu
      if (iopt(17).ne.0) then
c  diag 31
        call reord6    ("t31 ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125436)
        call add       ("gpph",-2.d0,"vvc1",1.d0,"s2  ")
        call mult6     ("s2  ","s1  ","s3  ",2)
        call add       ("s3  ",-1.d0,"vu0 ",1.d0,"vu0 ")
c  diag 33
        call reord6    ("t31 ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126345)
        call mult6     ("s1  ","gpph","s2  ",2)
        call reorder   ("s2  ","s1  ",1423,4)
        call add       ("s1  ",-1.d0,"vu0 ",1.d0,"vu0 ")
      end if
      if (iopt(17).ge.2) then
c    diag 36,38
        call reord6    ("t31 ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125634)
        call add       ("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",2)
        call reord6    ("t31 ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126543)
        call mult6     ("s1  ","va8 ","s2  ",2)
        call add       ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
        call mult      ("s3  ","t1hv","s1  ",311)
        call reorder   ("s1  ","s2  ",1423,4)
        call add       ("s2  ",1.d0,"vu0 ",1.d0,"vu0 ")
      end if
      call ordtb0_chu
      call sqr("vu0 ","s2  ")
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
c  dg3n1b
        call mult6      ("thv ","ghhh","s1  ",1)
        call reord6     ("s1  ","s2  ",124365)
        call reord7     ("s2  ","v30 ",132456)
c  dg3n1c
        call mult6      ("tv6 ","ghhh","s1  ",1)
        call reord6     ("s1  ","s2  ",124563)
        call reord7     ("s2  ","s1  ",132456)
        call add        ("s1  ",-1.d0,"v30 ",-1.d0,"v30 ")
c  dg3n2a
        call reorder    ("gpph","s1  ",1243,4)
        call mult6      ("tv3 ","s1  ","s2  ",1)
        call reord6     ("s2  ","s1  ",124356)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",125436)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c  dg3n2b
        call mult6      ("tv3 ","gphp","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125463)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c  dg3n2c,f1
        call reorder    ("tc  ","s1  ",1432,4)
        call mult6      ("gghp","s1  ","s3  ",1)
        call reorder    ("thv ","s1  ",2341,4)
        call mult6      ("gggg","s1  ","s2  ",1)
        call cmult      ("s2  ",-1.d0,"s1  ")
        call fold6      ("s3  ","s1  ",idum,2)
        call reord6     ("s3  ","s1  ",124356)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
      else if (iopt(17).eq.3) then
c  GROUP I
c   diag 2
        call mult6      ("thv ","bx12","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call reord7     ("s2  ","v30 ",132456)
c   diag 3
        call mult6      ("tv6 ","bx12","s1  ",1)
        call reord6     ("s1  ","s2  ",124653)
        call reord7     ("s2  ","s1  ",132456)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP II
c   diag 1
        call mult6      ("tv3 ","bx13","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125436)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 2
        call mult6      ("tv3 ","bx14","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125463)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 3
        call box24_chu  (0)
        call mult6      ("bx24","tc4 ","s1  ",1)
        call reord6     ("s1  ","s2  ",126354)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP III
c   diag 2
        call reord6     ("t31 ","s1  ",124536)
        call mult6      ("s1  ","bx15","s2  ",3)
        call reord6     ("s2  ","s1  ",125364)
        call reord7     ("s1  ","s2  ",132456)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP IV
c   diag 2
        call reord6     ("t31 ","s1  ",124635)
        call mult6      ("s1  ","bx16","s2  ",3)
        call reord6     ("s2  ","s1  ",125364)
        call reord7     ("s1  ","s2  ",132456)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 3
        call reord6     ("t31 ","s1  ",125634)
        call mult6      ("s1  ","bx16","s2  ",3)
        call reord6     ("s2  ","s1  ",125463)
        call reord7     ("s1  ","s2  ",132456)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP V
c   diag 2
        call reord6     ("t31 ","s1  ",124635)
        call mult6      ("s1  ","bx17","s2  ",3)
        call reord6     ("s2  ","s1  ",125346)
        call reord7     ("s1  ","s2  ",132456)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 3
        call reord6     ("t31 ","s1  ",125634)
        call mult6      ("s1  ","bx17","s2  ",3)
        call reord6     ("s2  ","s1  ",125643)
        call reord7     ("s1  ","s2  ",132456)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP VI
c   diag 1
        call reord6     ("t31 ","s1  ",124356)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",125634)
        call mult6      ("s1  ","bx18","s2  ",3)
        call reord6     ("s2  ","s1  ",125346)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",126345)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 2
        call box25_chu  (0)
        call reord6     ("t30 ","s1  ",126345)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",126543)
        call mult6      ("bx25","s1  ","s2  ",3)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP IX
c   diag 3
        call reord6     ("t31 ","s1  ",124563)
        call mult6      ("s1  ","bx19","s2  ",5)
        call reord6     ("s2  ","s1  ",126354)
        call reord7     ("s1  ","s2  ",132456)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
      end if
      return
      end
c ******************************************
      subroutine tuinit_chu
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif
      call move("vu0 ","thu ")
      call fold_ccsd("thu ","thuf")
      call diveps("thu ")
      call cmult ("thu ",-1.d0,"thu ")
      if (iopt(17).ne.0) then
        call move    ("v30 ","t32 ")
        call reord6  ("t32 ","s1  ",123546)
        call reord7  ("s1  ","s2  ",213456)
        call add     ("t32 ",1.d0,"s2  ",1.d0,"t32 ")
        call dveps6  ("t32 ")
        call cmult   ("t32 ",-1.d0,"t32 ")
      end if
      return
      end
c **************************************************
      subroutine ordta0_chu
!------------------------------------------------------------------------
c  tv1 = tv 1324
c  tvm1 = tvm 1324
c  tv2 = tv 1423
c  tvm2 = tvm 1423
c  tv3 = tv 1432
c  tcm4 = tcm 3412
c  t1cor = t1c 21
c  tb1 = tnew 1324
c  tb4 = tnew 3412
c  tb15 = tnew 1342
!------------------------------------------------------------------------
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call reorder("tc  ","tc1 ",1324,4)
      call reorder("tc  ","tc3 ",1423,4)
      call reorder("tc  ","tc4 ",3412,4)
      call reorder("thv ","tv1 ",1324,4)
      call reorder("thv ","tv2 ",1423,4)
      call reorder("thv ","tv3 ",1432,4)
      call reorder("thv ","tv6 ",1243,4)
      call reorder("tvm ","tvm1",1324,4)
      call reorder("tvm ","tvm2",1423,4)
      call reorder("tcm ","tcm2",1423,4)
      call reorder("tcm ","tcm4",3412,4)
      call reorder("t1c ","t1cr",1243,4)
      call reorder("t1hv","t1vr",1243,4)
      call clear("tb1 ",2222,1010)
      call clear("tb4 ",2222,11)
      call clear("tb15",2222,1010)
      call clear("thuf",2222,1111)
      call prtime("ordta0  ",0,1)
      return
      end
c ******************************************************
      subroutine ordtb0_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      !call reorder("tb1 ","s2  ",1324,4,0)
      call reorder("tb1 ","s2  ",1324,4)
      call add("vu0 ",1.d0,"s2  ",1.d0,"vu0 ")
      !call reorder("tb4 ","s2  ",3412,4,0)
      call reorder("tb4 ","s2  ",3412,4)
      call add("vu0 ",1.d0,"s2  ",1.d0,"vu0 ")
      !call reorder("tb15","s2  ",1342,4,0)
      call reorder("tb15","s2  ",1342,4)
      call add("vu0 ",1.d0,"s2  ",1.d0,"vu0 ")
      call prtime("ordtb0  ",0,1)
      return
      end
c *****************************************************
      subroutine ordtua_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      !call reorder("thu ","ta3 ",1432,4,1)
      call reorder("thu ","ta3 ",1432,4)
      !call reorder("tm  ","tm4 ",3412,4,1)
      call reorder("tm  ","tm4 ",3412,4)
      !call reorder("tm  ","tm6 ",1243,4,1)
      call reorder("tm  ","tm6 ",1243,4)
      call clear("tb1 ",2222,1010)
      call clear("tb2 ",2222,1001)
      call clear("tb4 ",2222,11)
      call clear("tb6 ",2222,1100)
      call clear("tb15",2222,1010)
      call prtime("ordta   ",0,1)
      return
      end
c ****************************************************
      subroutine caltu3_chu
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      call move          ("v30 ","t3nw")
      if (iopt(17).lt.3) then
c  dg3n1a
        call mult6       ("thu ","phhh","s1  ",1)
        call reord6      ("s1  ","s2  ",124356)
        call add         ("s2  ",-1.d0,"t3nw",1.d0,"t3nw")
      else if (iopt(17).eq.3) then
c  GROUP I
c   diag 1
        call mult6       ("thu ","bx2 ","s1  ",1)
        call reord6      ("s1  ","s2  ",126354)
        call add         ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP II
c   diag 3
        call box24_chu   (1)
        call mult6       ("bx24","tc4 ","s1  ",1)
        call reord6      ("s1  ","s2  ",126354)
        call add         ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP III
c   diag 1
        call reord6      ("t32 ","s1  ",124536)
        call mult6       ("s1  ","bx3 ","s2  ",3)
        call reord6      ("s2  ","s1  ",125346)
        call add         ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
c  GROUP IV
c   diag 1
        call reord6      ("t32 ","s1  ",124635)
        call mult6       ("s1  ","bx4 ","s2  ",3)
        call reord6      ("s2  ","s1  ",125346)
        call add         ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP V
c   diag 1
        call reord6      ("t32 ","s1  ",124635)
        call mult6       ("s1  ","bx5 ","s2  ",3)
        call reord6      ("s2  ","s1  ",125364)
        call add         ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP VI
c   diag 2
        call box25_chu   (1)
        call reord6      ("t30 ","s1  ",126345)
        call reord7      ("s1  ","s2  ",132565)
        call reord6      ("s2  ","s1  ",126543)
        call mult6       ("bx25","s1  ","s2  ",3)
        call add         ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP VII
c   diag 1
        call mult6       ("t32 ","bx6a","s1  ",3)
        call add         ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c   diag 2
        call reord6      ("t32 ","s1  ",123645)
        call mult6       ("s1  ","bx6 ","s2  ",3)
        call reord6      ("s2  ","s1  ",123564)
        call add         ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c   diag 3
        call mult6       ("thv ","hhhh","s1  ",1)
        call mult6       ("s1  ","t1hv","s2  ",5)
        call reord6      ("s2  ","s1  ",126345)
        call reord7      ("s1  ","s2  ",132456)
        call add         ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
c   diag 4
        call mult6       ("tv6 ","hhhh","s1  ",1)
        call mult6       ("s1  ","t1hv","s2  ",5)
        call reord6      ("s2  ","s1  ",126543)
        call reord7      ("s1  ","s2  ",132456)
        call add         ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP VIII
c   diag 1
        call mult6       ("t32 ","bx8 ","s1  ",5)
        call add         ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
c   diag 2
        call reord6      ("t32 ","s1  ",123564)
        call mult6       ("s1  ","bx8 ","s2  ",5)
        call reord6      ("s2  ","s1  ",123645)
        call add         ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP IX
c   diag 1
        call reord6      ("t32 ","s1  ",124563)
        call mult6       ("s1  ","bx7 ","s2  ",5)
        call reord6      ("s2  ","s1  ",126345)
        call add         ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
      end if
      call prtime      ("caltu3  ",0,0)
      return
      end
C************************************************************************
      subroutine caltu2_chu
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif
      call move("vu0 ","tnew")
      call cmult("tnew",0.5d0,"tnew")
      call diag2_chu
      call diag7_chu
      call diag14_chu
      call diag15_chu
      call diag27_chu
      if (iopt(17).ne.0) then
c  diag 30
        call reord6  ("t32 ","s1  ",124356)
        call add     ("vd14",2.d0,"vd12",-1.d0,"s2  ")
        call mult6   ("s1  ","s2  ","s3  ",2)
        call add     ("s3  ",-1.d0,"tnew",1.d0,"tnew")
c  diag 32
        call reord6  ("t32 ","s1  ",126345)
        call mult6   ("s1  ","vd14","s2  ",2)
        call add     ("s2  ",1.d0,"tb6 ",1.d0,"tb6 ")
      end if
      if (iopt(17).ge.2) then
c    diag 34,35      
        call add       ("va1 ",2.d0,"va2 ",-1.d0,"s1  ")
        call mult      ("s1  ","t1c ","s3  ",202)
c g1 Feb 95
        call add       ("s3  ",1.d0,"ph  ",1.d0,"s3  ")
 
        call reord6    ("t32 ","s1  ",124536)
        call reord6    ("t32 ","s2  ",124635)
        call add       ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
        call mult6     ("s1  ","s3  ","s2  ",0)
        call add       ("s2  ",1.d0,"tnew",1.d0,"tnew")
c    diag 37,39
        call reord6    ("t32 ","s1  ",124356)
        call add       ("pphh",-2.d0,"va6 ",1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",2)
        call reord6    ("t32 ","s1  ",126354)
        call mult6     ("s1  ","pphh","s2  ",2)
        call add       ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
        call mult      ("s3  ","t1cr","s1  ",311)
        call add       ("s1  ",1.d0,"tnew",1.d0,"tnew")
      end if
      call prtime  ("caltu2  ",0,0)
      return
      end
c ******************************************************
      subroutine ordtub_chu
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      !call reorder("tb1 ","s2  ",1324,5,1)
      call reorder("tb1 ","s2  ",1324,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      !call reorder("tb2 ","s2  ",1423,5,1)
      call reorder("tb2 ","s2  ",1423,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      !call reorder("tb4 ","s2  ",3412,5,1)
      call reorder("tb4 ","s2  ",3412,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      !call reorder("tb6 ","s2  ",1243,5,1)
      call reorder("tb6 ","s2  ",1243,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      !call reorder("tb15","s2  ",1342,5,1)
      call reorder("tb15","s2  ",1342,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call sqr("tnew","s2  ")
!     call fold("tnew","thuf")
      call fold_ccsd("tnew","thuf")
      call dfoldu_chu
      call diveps("tnew")
      call cmult("tnew",-1.d0,"tnew")
      if (iopt(17).ne.0) then
        call reord6 ("t3nw","s1  ",123546)
        call reord7 ("s1  ","s2  ",213456)
        call add    ("t3nw",1.d0,"s2  ",1.d0,"t3nw")
        call dveps6 ("t3nw")
        call cmult  ("t3nw",-1.d0,"t3nw")
      end if
      call prtime("ordtub  ",0,1)
      return
      end
c *************************************************
      subroutine diag7_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tm  ","hhhh","s2  ",222)
      call add("tnew",1.d0,"s2  ",0.5d0,"tnew")
      call prt("s2  ","diag7   ",3)
      call prtime("diag7   ",0,1)
      return
      end
c **********************************************
      subroutine diag8_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("ggpp","tcm4","s2  ",222)
      call add("vu0 ",1.d0,"s2  ",0.5d0,"vu0 ")
      call prt("s2  ","2diag8  ",3)
      call prtime("diag8   ",0,1)
      return
      end
c *************************************************
      subroutine diag3_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("gp  ","tv3 ","s1  ",131)
      call add("vu0 ",1.d0,"s1  ",1.d0,"vu0 ")
      call prt("s1  ","diag3   ",3)
      call prtime("diga3   ",0,1)
      return
      end
c **************************************************
      subroutine diag2_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("thu ","hh  ","s2  ",311)
      call add("tnew",1.d0,"s2  ",-1.d0,"tnew")
      call prt("s2  ","-diag2  ",3)
      call prtime("diga2   ",0,1)
      return
      end
c *******************************************************
      subroutine diag4_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add("vva1",2.d0,"ghph",-1.d0,"s2  ")
      call mult("s2  ","tv1 ","s1  ",222)
      call add("tb1 ",1.d0,"s1  ",1.d0,"tb1 ")
      call prt("s1  ","diag4   ",3)
      call prtime("diag4   ",0,1)
      return
      end
c ************************************************************
      subroutine diag5_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vva1","tvm2","s1  ",222)
      call add("tb1 ",1.d0,"s1  ",-1.d0,"tb1 ")
      call prt("s1  ","-diag5  ",3)
      call prtime("diag5   ",0,1)
      return
      end
c ************************************************************
      subroutine diag6_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("ghph","tvm2","s1  ",222)
      call add("tb15",1.d0,"s1  ",-1.d0,"tb15")
      call prt("s1  ","-diag6  ",3)
      call prtime("diag6   ",0,1)
      return
      end
c *******************************************************************
      subroutine diag9_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl diags 9 10 20 23
      call add("tv1 ",-1.d0,"tvm2",1.d0,"s2  ")
      call add("va1 ",-2.d0,"va2 ",1.d0,"s3  ")
      call mult("s2  ","s3  ","s1  ",222)
c   diag 20
      !call reorder("gpph","s2  ",1243,5,0)
      call reorder("gpph","s2  ",1243,5)
      !call reorder("gpph","s3  ",1342,5,0)
      call reorder("gpph","s3  ",1342,5)
      call add("s2  ",2.d0,"s3  ",-1.d0,"s2  ")
      call mult("s2  ","t1cr","s3  ",311)
      !call reorder("s3  ","s2  ",1423,5,0)
      call reorder("s3  ","s2  ",1423,5)
      call add("s1  ",1.d0,"s2  ",1.d0,"s1  ")
c   diag 23
      call add("vd12",-2.d0,"vd14",1.d0,"s2  ")
      call mult("t1hv","s2  ","s3  ",131)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
c -------------
      call mult("tv1 ","s1  ","s2  ",222)
      call add("tb1 ",1.d0,"s2  ",1.d0,"tb1 ")
      call prt("s2  ","diag9   ",3)
      call prtime("diag9   ",0,1)
      return
      end
c **************************************************
      subroutine diag11_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl diags 11 24
      call mult("tvm2","va1 ","s1  ",222)
c   diag 24
      call mult("t1hv","vd12","s2  ",131)
      call add("s1  ",0.5d0,"s2  ",1.d0,"s1  ")
c -----------------
      call mult("s1  ","tvm2","s2  ",222)
      call add("tb1 ",1.d0,"s2  ",1.d0,"tb1 ")
      call prt("s2  ","diag11  ",3)
      call prtime("diag11  ",0,1)
      return
      end
c **************************************************
      subroutine diag12_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl diags 12 25
      call mult("tvm2","va2 ","s1  ",222)
c   diag 25
      call mult("t1hv","vd14","s2  ",131)
      call add("s1  ",0.5d0,"s2  ",1.d0,"s1  ")
c ------------------
      call mult("s1  ","tvm2","s2  ",222)
      call add("tb15",1.d0,"s2  ",1.d0,"tb15")
      call prt("s2  ","diag12  ",3)
      call prtime("diag12  ",0,1)
      return
      end
c ************************************************
      subroutine diag13_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl diags 13 18
      call add("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult("tvm1","s2  ","s1  ",113)
c   diag 18
      call add("vvc1",2.d0,"gpph",-1.d0,"s2  ")
      call mult("s2  ","t1c ","s3  ",202)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
c -----------
      call mult("s1  ","tv3 ","s3  ",131)
      call add("vu0 ",1.d0,"s3  ",1.d0,"vu0 ")
      call prt("s3  ","diag13  ",3)
      call prtime("diag13  ",0,1)
      return
      end
c ****************************************************
      subroutine diag14_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl diags 14 19
      call add("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
      call mult("tcm4","s2  ","s1  ",113)
c   diag 19
      call add("vd11",-2.d0,"vd7 ",1.d0,"s3  ")
      call mult("s3  ","t1c ","s2  ",202)
      call add("s1  ",1.d0,"s2  ",1.d0,"s1  ")
c --- pdig14 to here
      call mult("thu ","s1  ","s2  ",311)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","diag14  ",3)
      call prtime("diag14  ",0,1)
      return
      end
c ****************************************************
      subroutine diag15_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tm  ","pphh","s1  ",222)
      call mult("s1  ","tcm4","s2  ",222)
      call add("tnew",1.d0,"s2  ",0.5d0,"tnew")
      call prt("s2  ","2diag15 ",3)
      call prtime("diag15  ",0,1)
      return
      end
c ****************************************************
      subroutine diag21_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl diags 21 17 22
      call mult("tv2 ","vvc1","s1  ",222)
c   diag 17
      call seltyp(2232,1100,"s2  ",5)
      !call reorder("s2  ","s3  ",2413,5,0)
      call reorder("s2  ","s3  ",2413,5)
      call add("s1  ",1.d0,"s3  ",-1.d0,"s1  ")
c   diag 22
      call mult("tv2 ","gpph","s2  ",222)
      !call reorder("s2  ","s3  ",3214,5,0)
      call reorder("s2  ","s3  ",3214,5)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
c --------------
      call mult("s1  ","t1cr","s2  ",311)
      call add("tb1 ",1.d0,"s2  ",-1.d0,"tb1 ")
      call prt("s2  ","-diag21 ",3)
      call prtime("diag21  ",0,1)
      return
      end
c ************************************************
      subroutine diag27_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl diags 27 29
c   diag 29
      call mult("thu ","ph  ","s2  ",311)
c   diag 27
      call mult("tm6 ","phhh","s1  ",222)
      !call reorder("s1  ","s3  ",1243,5,1)
      call reorder("s1  ","s3  ",1243,5)
      call add("s2  ",-1.d0,"s3  ",1.d0,"s2  ")
c ---------------
      call mult("s2  ","t1cr","s1  ",311)
      call add("tnew",1.d0,"s1  ",1.d0,"tnew")
      call prt("s1  ","diag27  ",3)
      call prtime("diag27  ",0,1)
      return
      end
c **********************************************************
      subroutine diag28_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("t1hv","ph  ","s2  ",111)
      !call reorder("thv ","s3  ",1432,5,0)
      call reorder("thv ","s3  ",1432,5)
      call mult("s2  ","s3  ","s1  ",131)
      call add("vu0 ",1.d0,"s1  ",-1.d0,"vu0 ")
      call prt("s1  ","-diag28 ",3)
      call prtime("diag28  ",0,1)
      return
      end
c *********************************************************
      subroutine diag26_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl diags 26 16
      !call reorder("tcm ","s1  ",3412,5,0)
      call reorder("tcm ","s1  ",3412,5)
      !call reorder("gpph","s2  ",1432,5,0)
      call reorder("gpph","s2  ",1432,5)
      call mult("s1  ","s2  ","s3  ",222)
c   diag 16
      call seltyp(2222,10,"s2  ",5)
      call add("s3  ",1.d0,"s2  ",1.d0,"s3  ")
c -------------
      call mult("s3  ","t1hv","s2  ",311)
      call add("tb4 ",1.d0,"s2  ",-1.d0,"tb4 ")
      call prt("s2  ","-diag26 ",3)
      call prtime("diag26  ",0,1)
      return
      end
c *********************************************************
      subroutine dgt1p3_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   compensates for overcounting
      call mult("t1hv","vd14","s2  ",131)
      call mult("t1hv","s2  ","s3  ",131)
      call mult("s3  ","t1cr","s2  ",311)
      call add("vu0 ",1.d0,"s2  ",-2.d0,"vu0 ")
      call prt("s2  ","-.5t1**3",3)
      call prtime("dgt1p3  ",0,1)
      return
      end
c ****************************************************
      subroutine dgt1p4_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   -t1**4 compensates for overcounting
      call mult("t1hv","pphh","s3  ",131)
      call mult("t1hv","s3  ","s2  ",131)
      call mult("t1cr","s2  ","s3  ",131)
      call mult("t1cr","s3  ","s2  ",131)
      call add("tb4 ",1.d0,"s2  ",-1.d0,"tb4 ")
      call prt("s2  ","-dgt1p4 ",3)
      call prtime("dgt1p4  ",0,1)
      return
      end
c *************************************************
      subroutine dfoldu_chu
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      call mult("thvf","ta3 ","s2  ",131)
c --------------
      call mult("thuf","t1vr","s3  ",311)
!     call unfold("s3  ","s1  ",3)
      call unfold_ccsd("s3  ","s1  ",3)
      call add("s2  ",1.d0,"s1  ",1.d0,"s2  ")
c --------------
      call mult("thuf","tm4 ","s3  ",222)
      call add("s2  ",1.d0,"s3  ",-0.5d0,"s2  ")
c =============
      call sqr("s2  ","s3  ")
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("tnew","tnew    ",3)
      if (iopt(17).eq.3) then
c  diag f1
        call t1t2      ("t1hv","thv ","t32 ","s1  ",123456,1)
        call t1t2      ("t1hv","thv ","s1  ","s2  ",123456,2)
        call reord7    ("s2  ","s1  ",312456)
        call reord6    ("s1  ","s2  ",124356)
        call reord7    ("s2  ","s1  ",132456)
        call reord6    ("s1  ","s2  ",125634)
        call mult6     ("thuf","s2  ","s1  ",3)
        call add       ("s1  ",-0.5d0,"t3nw",1.d0,"t3nw")
c  diag f2
        call reorder   ("thv ","s1  ",2341,4)
        call mult6     ("thuf","s1  ","s2  ",1)
        call reord6    ("s2  ","s1  ",124356)
        call cmult     ("s1  ",-1.d0,"s1  ")
        call fold6     ("t3nw","s1  ",idum,2)
c  diag f3
        call reord7    ("t32 ","s1  ",231456)
        call reord6    ("s1  ","s2  ",124563)
        call mult6     ("s2  ","thvf","s1  ",5)
        call reord6    ("s1  ","s2  ",126345)
        call reord7    ("s2  ","s1  ",312456)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
      end if
      call prtime("dfoldu  ",0,1)
      return
      end
C ********************************************************************
      subroutine rdbox_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call rdfile  (20,1,"bx2 ",-1)
      call reads   (20,2,"s1  ")
      call reorder ("s1  ","bx3 ",2143,4)
      call reads   (20,3,"s1  ")
      call reorder ("s1  ","bx4 ",2143,4)
      call reads   (20,4,"s1  ")
      call reorder ("s1  ","bx5 ",2143,4)
      call rdfile  (20,5,"bx6 ",-1)
      call rdfile  (20,6,"bx7 ",-1)
      call rdfile  (20,7,"bx8 ",-1)
      call rdfile  (31,0,"bx12",-1)
      call rdfile  (31,1,"bx13",-1)
      call rdfile  (31,2,"bx14",-1)
      call rdfile  (31,3,"bx15",-1)
      call rdfile  (31,4,"bx16",-1)
      call rdfile  (31,5,"bx17",-1)
      call rdfile  (31,6,"bx18",-1)
      call rdfile  (31,7,"bx19",-1)
      call rdfile  (31,8,"bx6a",-1)
      return
      end
C  *******************************************************************
      subroutine box24_chu (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      if (iparm.eq.0) then
c  diag 1
        call move     ("gghp","bx24")
c  diag 2
        call add      ("vvc1",-2.d0,"gpph",1.d0,"s1  ")
        call mult     ("tv1 ","s1  ","s3  ",222)
c  diag 3
        call mult     ("tvm2","vvc1","s1  ",222)
        call add      ("s1  ",1.d0,"s3  ",1.d0,"s1  ")
        call reorder  ("s1  ","s2  ",1324,4)
        call add      ("s2  ",-1.d0,"bx24",1.d0,"bx24")
c  diag 4
        call mult     ("gpph","tvm2","s1  ",222)
        call reorder  ("s1  ","s2  ",1342,4)
        call add      ("s2  ",-1.d0,"bx24",1.d0,"bx24")
c  diag 8
        call mult     ("gphh","t1hv","s1  ",311)
        call reorder  ("s1  ","s2  ",1432,4)
        call add      ("s2  ",-1.d0,"bx24",1.d0,"bx24")
c  diag 9
        call mult     ("t1hv","ghph","s1  ",131)
        call add      ("s1  ",-1.d0,"bx24",1.d0,"bx24")
c  diag 10
        call reorder  ("ggpp","s1  ",1243,4)
        call mult     ("s1  ","t1cr","s2  ",311)
        call reorder  ("s2  ","s1  ",1243,4)
        call add      ("s1  ",1.d0,"bx24",1.d0,"bx24")
      else if (iparm.eq.1) then
        call clear    ("bx24",2223,1100)
c  diag 5
        call mult     ("tm  ","vd14","s1  ",222)
        call add      ("s1  ",1.d0,"bx24",1.d0,"bx24")
c g1 Apr 95
        call mult     ("tm  ","ph  ","s1  ",311)
        call add      ("bx24",1.d0,"s1  ",-.5d0,"bx24")
c  diag 6
        call t1t2     ("t1hv","thv ","t32 ","s1  ",123456,2)
        call cmult    ("t1c ",0.5d0,"s2  ")
        call t1t2     ("s2  ","thu ","s1  ","s3  ",123456,3)
        call reord6   ("s3  ","s1  ",124356)
        call add      ("pphh",-2.d0,"va6 ",1.d0,"s2  ")
        call mult6    ("s1  ","s2  ","s3  ",2)
        call add      ("s3  ",1.d0,"bx24",1.d0,"bx24")
c  diag 7
        call t1t2     ("t1hv","thv ","t32 ","s1  ",123456,1)
        call t1t2     ("t1hv","thv ","s1  ","s2  ",123456,2)
        call t1t2     ("t1c ","tm  ","s2  ","s1  ",123456,3)
        call reord6   ("s1  ","s2  ",126354)
        call mult6    ("s2  ","pphh","s1  ",2)
        call add      ("s1  ",1.d0,"bx24",1.d0,"bx24")
      end if
      return
      end 
C  *******************************************************************
      subroutine box25_chu (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      if (iparm.eq.0) then
c  diag 1
        call cmult    ("ggpp",0.5d0,"bx25")
c  diag 2
        call mult     ("t1hv","gpph","s1  ",131)
        call add      ("s1  ",-1.d0,"bx25",1.d0,"bx25")
      else if (iparm.eq.1) then
        call clear    ("bx25",2233,1100)
c  diag 3   
        call mult     ("tm  ","pphh","s1  ",222)
        call add      ("s1  ",0.5d0,"bx25",1.d0,"bx25")
      end if
      return
      end
C  *******************************************************************
      subroutine wrtbox_chu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c
      call box24_chu(0)
      call wrtfile  (32,0,"bx24")
      call box24_chu(1)
      call reads    (32,0,"s1  ")
      call add      ("s1  ",1.d0,"bx24",1.d0,"bx24")
      call wrtfile  (32,0,"bx24")
c
      call box25_chu(0)
      call wrtfile  (32,1,"bx25")
      call box25_chu(1)
      call reads    (32,1,"s1  ")
      call add      ("s1  ",1.d0,"bx25",1.d0,"bx25")
      call wrtfile  (32,1,"bx25")
c
      return
      end
C  *******************************************************************
