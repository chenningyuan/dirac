!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine chp 
!---------------------------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
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
!     parameter (ndint=3000000)
#include "param_ccsd.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)
#include "cfiles.h"
!     common /cint/ nint,idum1,vint(ndint),iint(ndint)
#include "cint.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
! miro: synchronized with actual dimensions
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /csdtime/ csdstim,csdltim
#include "csdtime.h"
#endif

      character*256 filenm,filenm6,filenm51,filenm52,filenm53,filenm56,
     1filenm20,filenm21,filenm31,filenm41,filenm42
      real*4 cputime

c     iopt(1) - print.
c     iopt(2) - no. of iterations.
c     iopt(3) - convergence = 10**(-iopt(3))
c     iopt(4) .ne.0 - do rle every iopt(4) iterations.
c     iopt(5) - print only elements .gt. 10**(-iopt(5))
c     iopt(6) - start avrg.
c     iopt(8) - fock.
c     iopt(9) - incomplete model space
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
      filenm52(1:1)=' '
      call getenv('tape52',filenm52)
      if (filenm52(1:1).eq.' ') filenm52='tape52'
      open(52,file=filenm52,form='UNFORMATTED')
      filenm53(1:1)=' '
      call getenv('tape53',filenm53)
      if (filenm53(1:1).eq.' ') filenm53='tape53'
      open(53,file=filenm53,form='UNFORMATTED')
      filenm56(1:1)=' '
      call getenv('tape56',filenm56)
      if (filenm56(1:1).eq.' ') filenm56='tape56'
      open(56,file=filenm56,form='UNFORMATTED')
      filenm20(1:1)=' '
      call getenv('tape20',filenm20)
      if (filenm20(1:1).eq.' ') filenm20='tape20'
      open(20,file=filenm20,form='UNFORMATTED')
      filenm21(1:1)=' '
      call getenv('tape21',filenm21)
      if (filenm21(1:1).eq.' ') filenm21='tape21'
      open(21,file=filenm21,form='UNFORMATTED')
      filenm31(1:1)=' '
      call getenv('tape31',filenm31)
      if (filenm31(1:1).eq.' ') filenm31='tape31'
      open(31,file=filenm31,form='UNFORMATTED')
      filenm41(1:1)=' '
      call getenv('tape41',filenm41)
      if (filenm41(1:1).eq.' ') filenm41='tape41'
      open(41,file=filenm41,form='UNFORMATTED')
      filenm42(1:1)=' '
      call getenv('tape42',filenm42)
      if (filenm42(1:1).eq.' ') filenm42='tape42'
      open(42,file=filenm42,form='UNFORMATTED')
      filenm(1:1)=' '
      call getenv ('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')
      write (*,*) '  $$$$$$$$$$$   chp    $$$$$$$$$$  '
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
      call rdfile(51,0,"tc  ",-1)
      call rdfile(51,1,"t1c ",-1)
      if (iopt(17).ne.0) then
        call rdfile(51,2,"t3c ",-1)
      end if
      call rdfile(52,0,"tv  ",-1)
      call rdfile(52,1,"t1v ",-1)
      call rdfile(52,2,"tvf ",-1)
      if (iopt(17).ne.0) then
        call rdfile(52,3,"t3v ",-1)
      end if
!     tg is prepared by chv, must be there 
      call rdfile(53,0,"tg  ",-1)
      call rdfile(53,1,"t1g ",-1)
      call rdfile(53,2,"tgf ",-1)
      if (iopt(17).ne.0) then
        call rdfile(53,3,"t3g ",-1)
      end if
      if (iopt(17).eq.3) then
        call rdbox
      end if
      call t1sqr("t1c ","t1c ","tc  ","tcm ")
      call t1sqr("t1v ","t1c ","tv  ","tvm ")
      call reorder("t1c ","s2  ",1243,4)  
      call t1sqr("t1g ","s2  ","tg  ","tgm ")
      call ordvv_chp
      call calvv0_chp
      call pdiag
      iter=0
      call tvinit_chp
      if (iopt(4).ne.0) then
        if (iopt(17).eq.0) then
          call dbdiag ("dum1","vv0 ","dum3","dum4",0100,1)
        else
          call dbdiag ("dum1","vv0 ","v30 ","dum4",0110,1)
        end if
        call wrtrle ("thp ",2,1)    
        if (iopt(17).ne.0) then
          call wrtrle ("t3p ",3,1)     
        end if
      endif

      call prt("thp ","thp     ",2)
      call fmax("thp ",tmax,itmax)
      write(6,2400) tmax,idec(itmax)
      call t1sqr("t1v ","t1gr","thp ","tm  ")
220   iter=iter+1
      call prtime("strt itr",iter,0)
      call ordtva
      if (iopt(17).ne.0) then
        call caltv3
      end if
      call caltv2_chp
      call ordtvb
      if (iopt(6).ne.0.and.iter.ge.mod(iopt(6),100)) then
        fnew=max(dfloat(iopt(6)/100),1.d0)
        if (fnew.le.6.5d0) then
          fnew=1.d0/fnew
        else
          fnew=fnew-5.d0
        endif
        call add("tnew",fnew/(1.d0+fnew),"thp ",1.d0/(1.d0+fnew),"tnew")
        if (iopt(17).ne.0) then
          call add("t3nw",fnew/(1.d0+fnew),
     1    "t3p ",1.d0/(1.d0+fnew),"t3nw")
        end if
        write (6,'('' t averaged with fnew '',f7.3)') fnew
      endif
      call fmax("tnew",tmax,itmax)
      t3max = 0.0d0 ! mi initialization
      if (iopt(17).ne.0) then
        call fmax6 ("t3nw",t3max,id3mx,ind3mx)
      end if
      call diffmax("thp ","tnew","s2  ",diff,idiff)
      if (iopt(17).ne.0) then
        call dfmax3 ("t3p ","t3nw","s2  ",df3,iddf3,inddf3)
      end if

      write (6,*) ' '
      write (6,'('' iteration # '',i3,'' diffmx='',
     1 e10.3,2x,i12)') iter,diff
      df3 = 0.d0 ! mi initialization
      if (iopt(17).ne.0) then
        write (6,2401) df3,idec(iddf3),idec(inddf3)
2401    format (24x,e10.3,2x,i12,2x,i12)
      end if
      write(6,2400) tmax,idec(itmax)
2400  format(' max t=',f15.9,2x,i12)
      if (iopt(17).ne.0) then
        write (6,2402) t3max,idec(id3mx),idec(ind3mx)
2402    format (' max t3=',f15.9,2x,i12,2x,i12)
      end if
      if (abs(tmax).gt.1.d10.or.abs(tmax).lt.1.d-10) then
        write (*,*) ' abort. tmax out of bounds'
        call abext('tmax    ')     
      endif
      call move("tnew","thp ")
      if (iopt(17).ne.0) then
        call move("t3nw","t3p ")
      end if
      call t1sqr("t1v ","t1gr","thp ","tm  ")
      call prt("tnew","new t   ",2)
      if (max(abs(diff),abs(df3)).lt.convrg) go to 910
      if (iopt(4).ne.0) go to 260
c   no rle
      if (iter.ge.iopt(2)) go to 920
      go to 220
c   with rle
260   continue
      itr1=mod(iter-1,iopt(4))+1
      call wrtrle ("thp ",2,itr1+1)
      if (iopt(17).ne.0) then
        call wrtrle  ("t3p ",3,itr1+1)
      end if
      if (itr1.lt.iopt(4)) go to 220
      if (iopt(17).eq.0) then
        call rle (itr1,0100,"dum1","thp ","dum3","dum4",ising)
      else
        call rle (itr1,0110,"dum1","thp ","t3p ","dum4",ising)
      end if
c  rle singular
      if (ising.eq.1) then
        if (iopt(17).eq.0) then
          call rle (itr1-1,0100,"dum1","thp ","dum3","dum4",ising)
        else
          call rle (itr1-1,0110,"dum1","thp ","t3p ","dum4",ising)
        end if
        if (ising.eq.1) then
          iopt(4)=0
          write(*,*) ' ***** rle singular, turned off'
          go to 220
        endif
      endif
      call t1sqr  ("t1v ","t1gr","thp ","tm  ")
      call wrtrle ("thp ",2,1)    
      if (iopt(17).ne.0) then
        call wrtrle ("t3p ",3,1)    
      end if
      if (iter.ge.iopt(2)) go to 920
      go to 220
c
c  different endings
910   write (6,*) ' chp   converged after iter ',iter
      go to 970
920   write (6,*) ' iteration limit'
      call abext('limit   ')     

970   call prt("thp ","last t  ",0)
      call prt("thpf","last tf ",0)
      call wrtfile(56,0,"thp ")
      call wrtfile(56,1,"thpf")
      if (iopt(17).ne.0) then
        call wrtfile(56,2,"t3p ")
        call wrtbox_chp
      end if
      write (*,*) nstrt(1,nflcur),' words used in igen,',
     1  nstrt(2,nflcur),' words used in gen'
      call prtime("end     ",iter,0)
!      stop
      end
c **************************************************************************************
      subroutine ordvv_chp
!---------------------------------------------------------------------------------------
c va: hhpp
c vb: hphp
c vc: hhhp
c vd: hppp
c vva: vhpp
c vvc: vhhp
c              hhpp 1324 = va1
c              hhpp 1423 = va2
c              hhpp 3412 = va4
c              hhpp 1243 = va6
c              hhpp 3421 = va8
c              hhhp 3412 = vc4
c              hhhp 3421 = vc8
c              hhhp 3124 = vc9
c              hhhp 3214 = vc10
c              hppp 2314 = vd7
c              hppp 2413 = vd11
c              hppp 2134 = vd12
c              hppp 2143 = vd14
c              vhpp 1324 = vva1
c              vhhp 1324 = vvc1
!---------------------------------------------------------------------------------------
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
!     parameter (ndt1=1000,ndint=3000000)
#include "param_ccsd.h"
!     common /cint/nint,idum1,vint(ndint),iint(ndint)
#include "cint.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)
#include "fock.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif
c
      call readv(nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call prtime("readv",nint,0)

      call seltyp(2233,0,"hhpp",4)
      call seltyp(2323,0,"hphp",4)
      call seltyp(2223,0,"hhhp",4)
      call seltyp(2322,0,"hphh",4)
      call seltyp(2222,0,"hhhh",4)
      call seltyp(2333,0,"hppp",4)
      call seltyp(3233,1000,"vhpp",4)
      call seltyp(3323,1000,"vphp",4)
      call seltyp(3223,1000,"vhhp",4)
      call seltyp(3232,1000,"vhph",4)
      call seltyp(2322,100,"hvhh",4)
      call seltyp(2223,1000,"ghhp",4)
      call seltyp(2232,1000,"ghph",4)
      call seltyp(2323,10,"hpgp",4)
      call seltyp(2233,1000,"ghpp",4)
      call seltyp(3223,1100,"vghp",4)
      call seltyp(3233,1100,"vgpp",4)
      call seltyp(2322,110,"hvgh",4)
      call reorder("hhpp","va1 ",1324,4)
      call reorder("hhpp","va2 ",1423,4)
      call reorder("hhpp","va3 ",1243,4)
      call reorder("hhpp","va4 ",3412,4)
      call reorder("hhpp","va6 ",1243,4)
      call reorder("hhpp","va8 ",3421,4)
      call reorder("hphh","vc1 ",1423,4)
      call reorder("hhhp","vc4 ",3412,4)
      call reorder("hhhp","vc8 ",3421,4)
      call reorder("hppp","vd7 ",2314,4)
      call reorder("hppp","vd11",2413,4)
      call reorder("hppp","vd12",2134,4)
      call reorder("hppp","vd14",2143,4)
      call reorder("vhpp","vva1",1324,4)
      call reorder("vhhp","vvc2",1243,4)
      call reorder("vhhp","vvc3",1342,4)
      call reorder("ghpp","vga1",1243,4)
      call reorder("ghhp","vgb1",3124,4)
      call selt1(33,0,"pp  ",4)
      call selt1(23,0,"hp  ",4)
      call selt1(22,0,"hh  ",4)
      call selt1(32,10,"vh  ",4)
      call selt1(23,10,"gp  ",4)
      call reorder("hp  ","ph  ",1243,4)
      if (iopt(17).ne.0) then
        call seltyp  (3332,1001,"vppg",4)
      end if
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
        call seltyp  (3333,1000,"vppp",4)
        call seltyp  (3232,1111,"vgvg",4)
        call seltyp  (2222,1000,"ghhh",4)
      end if
      if (iopt(17).eq.3) then
        call seltyp4 (3333,0,"pppp",4)
        call seltyp  (2322,0010,"hpgh",4)
        call seltyp  (3322,0010,"ppgh",4)
        call seltyp  (3222,1000,"vhhh",4)
        call seltyp  (3322,1001,"vphg",4)
        call seltyp  (2333,1000,"gppp",4)
      end if
      call prtime("ordv    ",0,0)
      return
      end
c**************************************************************
      subroutine pdiag
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call pdig14b
      call pdig13a
      call pdig12a
      return
      end
c ******************************************
      subroutine tvinit_chp
!----------------------------------------------------------
!
!
!
!----------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      call move    ("vv0 ","thp ")
      call fold_ccsd ("thp ","thpf")
      call diveps  ("thp ")
      if (iopt(17).gt.0) then
        call move     ("v30 ","t3p ")
        call dveps6   ("t3p ")
      end if
      return
      end 

c **************************************************
      subroutine ordta0_chp
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder("tc  ","tc1 ",1423,4)
      call reorder("tv  ","tv1 ",1324,4)
      call reorder("tvm ","tvm2",1423,4)
      call reorder("tg  ","tg1 ",3142,4)
      call reorder("tg  ","tg2 ",3214,4)
      call reorder("tgm ","tgm2",4132,4)
      call reorder("t1c ","t1cr",1243,4)
      call reorder("t1g ","t1gr",1243,4)
      call reorder("tg  ","tgr ",3412,4)
      call clear("tb1 ",2323,110)
      call clear("tb2 ",3322,1001)
      call clear("tb4 ",3223,1010)
      call clear("tb15",3223,1100)
      call clear("thpf",3232,1111)
      call prtime("ordta0  ",0,1)
      return
      end
c ******************************************************
      subroutine calvv0_chp
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      call ordta0_chp
      call seltyp(3232,1001,"vv0 ",3)
      call diag7_chp
      call diag16_chp
      call diag17
      call diag9_chp
      call dgt1p3
      call dgt1p4
      call diag11
      call diag24
      call diag15
      call diag13b
      call diag20b
      call diag21b
      call diag26
      call diag27a
      call diag2b
      call diag3a
      call diag8
      call diag4a
      call diag4b
      call diag5a
      call diag5b
      call diag6a
      call diag22b
      call diag6b_chp
      call diag28b
      call diag14a
      call diag22a
      call diag27b
      if (iopt(17).ne.0) then
c  diags 30c,37c
        call add       ("ghpp",2.d0,"vga1",-1.d0,"boxB")
        if (iopt(17).ge.2) then
          call add     ("va4 ",-2.d0,"va8 ",1.d0,"s1  ")
          call mult    ("t1g ","s1  ","s2  ",131)
          call reorder ("s2  ","s1  ",1432,4)
          call add     ("s1  ",1.d0,"boxB",1.d0,"boxB")
        end if
        call reord6    ("t3v ","s1  ",124356)
        call mult6     ("s1  ","boxB","s2  ",2)
        call add       ("s2  ",1.d0,"vv0 ",1.d0,"vv0 ")
c  diags 31b,36b
        call reorder   ("vhhp","s1  ",1324,4)
        call add       ("vhhp",-2.d0,"s1  ",1.d0,"boxC")
        if (iopt(17).ge.2) then
          call add     ("hhpp",-2.d0,"va6 ",1.d0,"s1  ")
          call mult    ("t1v ","s1  ","s2  ",131)
          call add     ("s2  ",1.d0,"boxC",1.d0,"boxC")
        end if
        call reord6    ("t3g ","s1  ",124653)
        call mult6     ("boxC","s1  ","s2  ",2)
        call reorder   ("s2  ","s1  ",1432,4)
        call add       ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
c  diags 32c,39c
        call cmult     ("ghpp",-1.d0,"boxF")            
        if (iopt(17).ge.2) then
          call mult    ("t1g ","va4 ","s1  ",131)
          call reorder ("s1  ","s2  ",1432,4)
          call add     ("s2  ",1.d0,"boxF",1.d0,"boxF")
        end if
        call reord6    ("t3v ","s1  ",126354)
        call mult6     ("s1  ","boxF","s2  ",2)
        call add       ("s2  ",1.d0,"vv0 ",1.d0,"vv0 ")
c  diags 33b,38b
        call move      ("vhhp","boxG")            
        if (iopt(17).ge.2) then
          call mult    ("t1v ","hhpp","s1  ",131)
          call add     ("s1  ",1.d0,"boxG",1.d0,"boxG")
        end if
        call reord6    ("t3g ","s1  ",126453)
        call mult6     ("boxG","s1  ","s2  ",2)
        call reorder   ("s2  ","s1  ",1432,4)
        call add       ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
      end if
      call ordtb0_chp
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
c  dg3n1b,2h
        call mult6     ("tv  ","vga1","s1  ",1)
        call reorder   ("tg  ","s2  ",1423,4)
        call mult6     ("vhph","s2  ","s3  ",1)
        call add       ("s3  ",-1.d0,"s1  ",1.d0,"s1  ")
        call reord6    ("s1  ","v30 ",125346)
c  dg3n1c,2g
        call mult6     ("tv  ","hpgp","s1  ",1)
        call mult6     ("vhph","tg2 ","s2  ",1)
        call add       ("s2  ",-1.d0,"s1  ",1.d0,"s1  ")
        call reord6    ("s1  ","s2  ",124365)
        call reord7    ("s2  ","s1  ",132456)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c  dg3n1f,2d
        call reorder   ("tv  ","s1  ",1243,4)
        call mult6     ("s1  ","hpgp","s2  ",1)
        call mult6     ("vvc2","tg2 ","s1  ",1)
        call add       ("s1  ",-1.d0,"s2  ",1.d0,"s1  ")
        call reord6    ("s1  ","s2  ",124563)
        call reord7    ("s2  ","s1  ",132456)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c  dg3n1g,2c
        call reorder   ("tg  ","s1  ",1432,4)
        call mult6     ("vppp","s1  ","s2  ",1)
        call reorder   ("tv  ","s1  ",1432,4)
        call mult6     ("s1  ","ghhh","s3  ",1)
        call add       ("s3  ",-1.d0,"s2  ",1.d0,"s1  ")
        call reord6    ("s1  ","s2  ",126345)
        call reord7    ("s2  ","s1  ",132456)
        call reord6    ("s1  ","s2  ",126453)
        call add       ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c  dg3nf2,1h
        call reorder   ("vgvg","s1  ",1342,4)
        call reorder   ("tg  ","s2  ",4321,4)
        call mult6     ("s1  ","s2  ","s3  ",1)
        call cmult     ("s3  ",-1.d0,"s1  ")
        call reorder   ("vppg","s2  ",1342,4)
        call mult6     ("s2  ","tc  ","s3  ",1)
        call fold6     ("s3  ","s1  ",idum,2)
        call reord6    ("s3  ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125346)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c  dg3n2f
        call reorder   ("hvgh","s1  ",2134,4)
        call reorder   ("tc  ","s2  ",1432,4)
        call mult6     ("s1  ","s2  ","s3  ",1)
        call reord6    ("s3  ","s1  ",124536)
        call add       ("s1  ",-1.d0,"v30 ",1.d0,"v30 ")
c  dg3nf1
        call reorder   ("vgvg","s1  ",1243,4)
        call reorder   ("tv  ","s2  ",2341,4)
        call mult6     ("s1  ","s2  ","s3  ",1)
        call reord6    ("s3  ","s1  ",124536)
        call fold6     ("v30 ","s1  ",idum,2)
      else if (iopt(17).eq.3) then
c  GROUP I
c   diag 1
        call reads     (31,1,"s1  ")
        call mult6     ("tv  ","s1  ","s2  ",1)
        call reord6    ("s2  ","v30 ",126345)
c   diag 2
        call reads     (31,2,"s1  ")
        call mult6     ("tv  ","s1  ","s2  ",1)
        call reord6    ("s2  ","s3  ",126345)
        call reord7    ("s3  ","s2  ",132456)
        call add       ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 4
        call reorder   ("tv  ","s2  ",1243,4)
        call mult6     ("s2  ","s1  ","s3  ",1)
        call reord6    ("s3  ","s1  ",126543)
        call reord7    ("s1  ","s2  ",132456)
        call add       ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 5
        call reads     (21,0,"s1  ")
        call reorder   ("tg  ","s2  ",1432,4)
        call mult6     ("s1  ","s2  ","s3  ",1)
        call reord6    ("s3  ","s1  ",126345)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126354)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 6
        call box31     (0)
        call mult6     ("bx31","tc  ","s1  ",1)
        call reord6    ("s1  ","s2  ",125346)
        call reord7    ("s2  ","s1  ",132456)
        call reord6    ("s1  ","s2  ",125346)
        call add       ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP II
c   diag 2
        call reads     (31,0,"s1  ")
        call reorder   ("tv  ","s2  ",1432,4)
        call mult6     ("s2  ","s1  ","s3  ",1)
        call reord6    ("s3  ","s2  ",125346)
        call reord7    ("s2  ","s1  ",132456)
        call reord6    ("s1  ","s2  ",126453)
        call add       ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 3
        call reads     (21,1,"s1  ")
        call mult6     ("s1  ","tg  ","s2  ",1)
        call reord6    ("s2  ","s1  ",126543)
        call reord7    ("s1  ","s2  ",132456)
        call add       ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 4
        call box32     (0)             
        call reorder   ("tc  ","s1  ",1432,4)
        call mult6     ("bx32","s1  ","s2  ",1)
        call reord6    ("s2  ","s1  ",124536)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 5
        call reads     (21,2,"s1  ")
        call mult6     ("s1  ","tg  ","s2  ",1)
        call reord6    ("s2  ","s3  ",126345)
        call reord7    ("s3  ","s2  ",132456)
        call add       ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 6
        call reorder   ("tg  ","s2  ",1243,4)
        call mult6     ("s1  ","s2  ","s3  ",1)
        call reord6    ("s3  ","s1  ",126345)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP III
c   diag 2
        call reads     (31,3,"s1  ")
        call reord6    ("t3v ","s2  ",124563)
        call mult6     ("s2  ","s1  ","s3  ",3)
        call reord6    ("s3  ","s1  ",126354)
        call reord7    ("s1  ","s2  ",132456)
        call add       ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 3
        call reads     (21,3,"s1  ")
        call reord6    ("t3g ","s2  ",124563)
        call mult6     ("s1  ","s2  ","s3  ",3)
        call reord6    ("s3  ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126345)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP IV 
c   diag 2
        call reads     (31,4,"s1  ")
        call reord6    ("t3v ","s2  ",124653)
        call mult6     ("s2  ","s1  ","s3  ",3)
        call reord6    ("s3  ","s2  ",126354)
        call reord7    ("s2  ","s3  ",132456)
        call add       ("s3  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 4
        call reord6    ("t3v ","s2  ",125643)
        call mult6     ("s2  ","s1  ","s3  ",3)
        call reord6    ("s3  ","s1  ",126453)
        call reord7    ("s1  ","s2  ",132456)
        call add       ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 5
        call reads     (21,4,"s1  ")
        call reord6    ("t3g ","s2  ",124653)
        call mult6     ("s1  ","s2  ","s3  ",3)
        call reord6    ("s3  ","s2  ",125346)
        call reord7    ("s2  ","s3  ",132456)
        call reord6    ("s3  ","s2  ",126345)
        call add       ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 6
        call reord6    ("t3g ","s2  ",125643)
        call mult6     ("s1  ","s2  ","s3  ",3)
        call reord6    ("s3  ","s1  ",126345)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126345)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP V 
c   diag 1
        call reads     (31,5,"s1  ")
        call reord6    ("t3v ","s2  ",124653)
        call mult6     ("s2  ","s1  ","s3  ",3)
        call reord6    ("s3  ","s1  ",126354)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 5
        call reads     (21,5,"s1  ")
        call reord6    ("t3g ","s2  ",124653)
        call mult6     ("s1  ","s2  ","s3  ",3)
        call reord6    ("s3  ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126543)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 6  
        call box33     (0)
        call reord6    ("t3c ","s1  ",124635)
        call mult6     ("bx33","s1  ","s2  ",3)
        call reord6    ("s2  ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",124635)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP VI
c   diag 2
        call reads     (21,6,"s1  ")
        call mult6     ("s1  ","t3g ","s2  ",3)
        call reord6    ("s2  ","s3  ",126435)
        call reord7    ("s3  ","s2  ",132456)
        call add       ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 3
        call reord6    ("t3g ","s2  ",123654)
        call mult6     ("s1  ","s2  ","s3  ",3)
        call reord6    ("s3  ","s1  ",126435)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 4
        call mult      ("t1g ","hhhh","s1  ",131)
        call reorder   ("tv  ","s2  ",1432,4)
        call mult6     ("s2  ","s1  ","s3  ",1)
        call reord6    ("s3  ","s1  ",126345)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126453)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP VII
c   diag 5
        call reads     (31,6,"s1  ")
        call mult6     ("t3v ","s1  ","s2  ",3)
        call add       ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 6
        call reord6    ("t3v ","s2  ",123654)
        call mult6     ("s2  ","s1  ","s3  ",3)
        call reord6    ("s3  ","s1  ",123654)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c   diag 7
        call mult      ("t1v ","gppp","s1  ",131)
        call mult6     ("s1  ","tc  ","s2  ",1)
        call reord6    ("s2  ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125436)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP VIII
c   diag 2
        call reads     (31,7,"s1  ")
        call mult6     ("t3v ","s1  ","s2  ",5)
        call reord7    ("s2  ","s1  ",132456)
        call reord6    ("s1  ","s2  ",123465)
        call add       ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c  GROUP IX
c   diag 3
        call reads     (21,7,"s1  ")
        call mult6     ("t3g ","s1  ","s2  ",5)
        call reord6    ("s2  ","s1  ",126345)
        call reord7    ("s1  ","s2  ",312456)
        call reord6    ("s2  ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126534)
        call add       ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
      end if
      call prtime      ("calvv0  ",0,0)
      return
      end 
C **********************************************************************
      subroutine ordtb0_chp
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder("tb1 ","s2  ",2143,4)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call reorder("tb2 ","s2  ",1324,4)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call reorder("tb4 ","s2  ",1243,4)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call reorder("tb15","s2  ",1342,4)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call prtime("ordtb0  ",0,1)
      return
      end
c *****************************************************
      subroutine ordtva
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder("thp ","ta1 ",1243,4)
      call reorder("thp ","ta3 ",1432,4)
      call reorder("tm  ","tm4 ",2341,4)
      call reorder("tm  ","tm6 ",1243,4)
      call reorder("thp ","ta2 ",1342,4)
      call clear("tb15",3223,1100)
      call clear("tb3 ",3232,1100)
      call clear("tb4 ",3223,1010)
      call clear("thpf",3232,1111)
      call prtime("ordta   ",0,1)
      return
      end
c ****************************************************
      subroutine caltv3
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      call move        ("v30 ","t3nw")
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
c  dg3n1d
        call mult6     ("ta1 ","hppp","s1  ",1)
        call reord6    ("s1  ","s2  ",124536)
        call add       ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
c  dg3n2b
        call reorder   ("thp ","s1  ",1342,4)
        call reorder   ("hhhp","s2  ",2341,4)
        call mult6     ("s1  ","s2  ","s3  ",1)
        call reord6    ("s3  ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125346)
        call add       ("s1  ",-1.d0,"t3nw",1.d0,"t3nw")
      else if (iopt(17).eq.3) then
c  GROUP I
c   diag 3
        call mult6     ("ta1 ","bx1 ","s1  ",1)
        call reord6    ("s1  ","s2  ",124635)
        call add       ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
c   diag 6
        call box31     (1)
        call mult6     ("bx31","tc  ","s1  ",1)
        call reord6    ("s1  ","s2  ",125346)
        call reord7    ("s2  ","s1  ",132456)
        call reord6    ("s1  ","s2  ",125346)
        call add       ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP II
c   diag 1
        call mult6     ("ta3 ","bx2 ","s1  ",1)
        call reord6    ("s1  ","s2  ",125346)
        call reord7    ("s2  ","s1  ",132456)
        call reord6    ("s1  ","s2  ",125436)
        call add       ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
c   diag 4
        call box32     (1)
        call reorder   ("tc  ","s1  ",1432,4)
        call mult6     ("bx32","s1  ","s2  ",1)
        call reord6    ("s2  ","s1  ",124536)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP III
c   diag 1
        call reord6    ("t3p ","s1  ",124536)
        call mult6     ("s1  ","bx3 ","s2  ",3)
        call reord6    ("s2  ","s1  ",125346)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP IV 
c   diag 1
        call reord7    ("t3p ","s1  ",132456)
        call reord6    ("s1  ","s2  ",124536)
        call mult6     ("s2  ","bx4 ","s1  ",3)
        call reord6    ("s1  ","s2  ",125346)
        call add       ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
c   diag 3
        call reord6    ("t3p ","s1  ",125634)
        call mult6     ("s1  ","bx4 ","s2  ",3)
        call reord6    ("s2  ","s1  ",125436)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP V 
c   diag 2
        call reord7    ("t3p ","s1  ",132456)
        call reord6    ("s1  ","s2  ",124536)
        call mult6     ("s2  ","bx5 ","s1  ",3)
        call reord6    ("s1  ","s2  ",125346)
        call reord7    ("s2  ","s1  ",132456)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c   diag 3
        call reord6    ("t3p ","s1  ",125634)
        call mult6     ("s1  ","bx5 ","s2  ",3)
        call reord6    ("s2  ","s1  ",125634)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c   diag 4
        call reord7    ("t3p ","s1  ",132456)
        call reord6    ("s1  ","s2  ",125634)
        call mult6     ("s2  ","bx5 ","s1  ",3)
        call reord6    ("s1  ","s2  ",125634)
        call reord7    ("s2  ","s1  ",132456)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c   diag 6  
        call box33     (1)
        call reord6    ("t3c ","s1  ",124635)
        call mult6     ("bx33","s1  ","s2  ",3)
        call reord6    ("s2  ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",124635)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP VI
c   diag 1
        call reord6    ("t3p ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125634)
        call mult6     ("s1  ","bx6a","s2  ",3)
        call reord6    ("s2  ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126345)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP VII
c   diag 1
        call t1t2      ("t1v ","tgr ","t3p ","s1  ",123456,1)
        call t1t2      ("t1c ","thp ","s1  ","s3  ",123456,3)
        call reord6    ("s3  ","s1  ",123546)
        call mult6     ("s1  ","pppp","s2  ",3)
        call reord6    ("s2  ","s1  ",123546)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw") 
        call reord6    ("s3  ","s1  ",123564)
        call reord7    ("s1  ","s2  ",321456)
        call mult6     ("s2  ","pppp","s1  ",3)
        call reord6    ("s1  ","s2  ",123645)
        call reord7    ("s2  ","s1  ",321456)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw") 
c   diag 2
        call reord6    ("t3p ","s1  ",123546)
        call mult6     ("s1  ","vd12","s2  ",3)
        call mult6     ("s2  ","t1cr","s1  ",5)
        call reord6    ("s1  ","s2  ",123645)
        call add       ("s2  ",-1.d0,"t3nw",1.d0,"t3nw") 
c   diag 3
        call reord6    ("t3p ","s1  ",123564)
        call mult6     ("s1  ","vd12","s2  ",3)
        call mult6     ("s2  ","t1cr","s1  ",5)
        call reord6    ("s1  ","s2  ",123546)
        call add       ("s2  ",-1.d0,"t3nw",1.d0,"t3nw") 
c   diag 4
        call reord6    ("t3p ","s1  ",123546)
        call mult6     ("s1  ","hhpp","s2  ",3)
        call reorder   ("tcm ","s1  ",3412,4)
        call mult6     ("s2  ","s1  ","s3  ",3)
        call reord6    ("s3  ","s1  ",123546)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw") 
c  GROUP VIII
c   diag 1
        call mult6     ("t3p ","bx7 ","s1  ",5)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c   diag 3
        call reord6    ("t3p ","s1  ",123564)
        call mult6     ("s1  ","bx7 ","s2  ",5)
        call reord6    ("s2  ","s1  ",123645)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c  GROUP IX   
c   diag 1
        call reord6    ("t3p ","s1  ",124563)
        call mult6     ("s1  ","bx8 ","s2  ",5)
        call reord6    ("s2  ","s1  ",126345)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c   diag 2
        call reord7    ("t3p ","s1  ",132456)
        call reord6    ("s1  ","s2  ",124563)
        call mult6     ("s2  ","bx8 ","s1  ",5)
        call reord6    ("s1  ","s2  ",126345)
        call reord7    ("s2  ","s1  ",132456)
        call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
      end if
      call prtime      ("caltv3  ",0,0)
      return
      end 
C  ********************************************************************
      subroutine caltv2_chp
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      call move("vv0 ","tnew")
      call diag2
      call diag3b
      call diag6
      call diag14b
      call diag13a
      call diag12a
      if (iopt(17).ne.0) then
c  diags 30b,37b
        call add       ("vd12",2.d0,"vd14",-1.d0,"boxA")
        if (iopt(17).ge.2) then
          call add     ("va4 ",-2.d0,"va8 ",1.d0,"s1  ")
          call mult    ("t1cr","s1  ","s2  ",131)
          call reorder ("s2  ","s1  ",1423,4)
          call add     ("s1  ",1.d0,"boxA",1.d0,"boxA")
        end if
        call reord6    ("t3p ","s1  ",125364)
        call mult6     ("s1  ","boxA","s2  ",2)
        call add       ("s2  ",1.d0,"tb4 ",1.d0,"tb4 ")
c  diags 31c,36c
        call reorder   ("hhhp","s1  ",1324,4)
        call add       ("hhhp",-2.d0,"s1  ",1.d0,"boxD")
        if (iopt(17).ge.2) then
          call add     ("hhpp",-2.d0,"va6 ",1.d0,"s1  ")
          call mult    ("t1c ","s1  ","s2  ",131)
          call add     ("s2  ",1.d0,"boxD",1.d0,"boxD")
        end if
        call reord6    ("t3p ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125436)
        call mult6     ("s1  ","boxD","s2  ",2)
        call reorder   ("s2  ","s1  ",1423,4)
        call add       ("s1  ",1.d0,"tnew",1.d0,"tnew")
c  diags 32b,39b
        call cmult     ("vd12",-1.d0,"boxE")            
        if (iopt(17).ge.2) then
          call mult    ("t1cr","va4 ","s1  ",131)
          call reorder ("s1  ","s2  ",1423,4)
          call add     ("s2  ",1.d0,"boxE",1.d0,"boxE")
        end if
        call reord7    ("t3p ","s1  ",132456)
        call reord6    ("s1  ","s2  ",125364)
        call mult6     ("s2  ","boxE","s1  ",2)
        call add       ("s1  ",1.d0,"tb4 ",1.d0,"tb4 ")
c  diags 33c,38c
        call move      ("hhhp","boxH")            
        if (iopt(17).ge.2) then
          call mult    ("t1c ","hhpp","s1  ",131)
          call add     ("s1  ",1.d0,"boxH",1.d0,"boxH")
        end if
        call reord6    ("t3p ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126435)
        call mult6     ("s1  ","boxH","s2  ",2)
        call add       ("s2  ",1.d0,"tb3 ",1.d0,"tb3 ")
c  diags 34b,35b,35c
        call add       ("va1 ",2.d0,"va2 ",-1.d0,"s1  ")
        call mult      ("s1  ","t1c ","boxI",202)
        call reord6    ("t3p ","s3  ",124536)
        call reord6    ("t3p ","s1  ",126534)
        call add       ("s1  ",-1.d0,"s3  ",2.d0,"s3  ")
        call reord7    ("t3p ","s1  ",132456)
        call reord6    ("s1  ","s2  ",124536)
        call add       ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
        call mult6     ("s3  ","boxI","s1  ",0)
        call add       ("s1  ",1.d0,"tnew",1.d0,"tnew")
      end if
      call prtime    ("caltv2  ",0,0)
      return
      end
c ******************************************************
      subroutine ordtvb

#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      i=1
      call reorder("tb4 ","s2  ",1243,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
10    call reorder("tb15","s2  ",1342,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call reorder("tb3 ","s2  ",1432,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      if (i.eq.77) go to 11
      call fold_ccsd("tnew","thpf")
c     call heff2
      i=77
      call clear("tb15",3223,1100)
      call clear("tb3 ",3232,1100)
      call dfoldhp
      if (iopt(17).eq.3) then
        call dfold3
      end if
      go to 10
11    call diveps("tnew")
      if (iopt(17).ne.0) then
        call dveps6 ("t3nw")     
      end if
      call prtime("ordtvb  ",0,0)
      return
      end
c *************************************************
      subroutine diag7_chp
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call seltyp(3233,100,"s2  ",5)
      call mult("tvm ","s2  ","s3  ",222)
      call add("vv0 ",1.d0,"s3  ",1.d0,"vv0 ")
      call prt("s3  ","diag7   ",3)
      call prtime("diag7   ",0,1)
      return
      end
c ************************************************
      subroutine diag2
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("ta1 ","pp  ","s2  ",311)
      call add("tb4 ",1.d0,"s2  ",1.d0,"tb4 ")
      call prt("s2  ","diag2   ",3)
      call prtime("diag2   ",0,1)
      return
      end
c ************************************************
      subroutine diag3b
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("ta3 ","hh  ","s2  ",311)
      call add("tb3 ",1.d0,"s2  ",-1.d0,"tb3 ")
      call prt("s2  ","-diag3b ",3)
      call prtime("diag3b  ",0,1)
      return
      end
c **********************************************
      subroutine diag6
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder("hphp","s1  ",1243,5)
      call mult("ta3 ","s1  ","s2  ",222)
      call add("tb15",1.d0,"s2  ",-1.d0,"tb15")
      call prt("s2  ","-diag6  ",3)
      call prtime("diag6   ",0,1)
      return
      end
c *********************************************
      subroutine pdig14b
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  incl 14b, 19b, 29b
C  14B
      call add("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
      call reorder("tcm ","s1  ",3412,4)
      call mult("s1  ","s2  ","p14b",113)
c  19b
      call add("vd11",2.d0,"vd7 ",-1.d0,"s2  ")
      call mult("t1c ","s2  ","s1  ",022)
      call add("p14b",1.d0,"s1  ",1.d0,"p14b")
c  29b
      call mult("t1cr","ph  ","s2  ",111)
      call add("p14b",1.d0,"s2  ",1.d0,"p14b")
      call prt("p14b","pdig14b ",3)
      call prtime("pdig14b ",0,1)
      return
      end
c *********************************************
      subroutine diag14b
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("ta1 ","p14b","s2  ",311)
      call add("tb4 ",1.d0,"s2  ",1.d0,"tb4 ")
      call prt("s2  ","diag14b ",3)
      call prtime("diag14b ",0,1)
      return
      end
c **************************************************
      subroutine pdig13a
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  incl 13a 28a 28a
c  13a
      call add("hhpp",-2.d0,"va3 ",1.d0,"s2  ")
      call mult("tcm ","s2  ","p13a",113)
c  18a
      call reorder("hhhp","s2  ",1324,5)
      call add("s2  ",-2.d0,"hhhp",1.d0,"s2  ")
      call mult("t1c ","s2  ","s1  ",022)
      call add("p13a",1.d0,"s1  ",1.d0,"p13a")
c  28a
      call mult("t1c ","hp  ","s1  ",111)
      call add("p13a",1.d0,"s1  ",1.d0,"p13a")
      call prt("p13a","pdig13a ",3)
      call prtime("pdig13a ",0,1)
      return
      end
c *********************************************
      subroutine diag13a
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("ta3 ","p13a","s2  ",311)
      call add("tb3 ",1.d0,"s2  ",1.d0,"tb3 ")
      call prt("s2  ","diag13a ",3)
      call prtime("diag13a ",0,1)
      return
      end
c **************************************************
      subroutine pdig12a
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  incl 12a, 22b, 25a
c  12a
      call reorder("hhpp","s1  ",3214,5)
      call mult("tc1 ","s1  ","p12a",222)
c 22b
      call reorder("hhhp","s1  ",1432,5)
      call mult("t1cr","s1  ","s2  ",131)
      call reorder("s2  ","s1  ",2134,5)
      call add("p12a",1.d0,"s1  ",1.d0,"p12a")
c  25a
      call mult("t1c ","vd11","s2  ",131)
      call add("p12a",1.d0,"s2  ",-1.d0,"p12a")
      call prt("p12a","pdig12a ",3)
      call prtime("pdig12a ",0,1)
      return
      end
c *********************************************
      subroutine diag12a
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("ta3 ","p12a","s2  ",222)
      call add("tb15",1.d0,"s2  ",1.d0,"tb15")
      call prt("s2  ","diag12a ",3)
      call prtime("diag12a ",0,1)
      return
      end
c *********************************************
      subroutine diag4a
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder("ghph","s1  ",2143,5)
      call add("vgb1",2.d0,"s1  ",-1.d0,"s2  ")
      call mult("tv1 ","s2  ","s1  ",222)
      call add("tb2 ",1.d0,"s1  ",1.d0,"tb2 ")
      call prt("s1  ","diag4a  ",3)
      call prtime("diag4a  ",0,1)
      return
      end
c *****************************************************
      subroutine diag4b
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add("vva1",2.d0,"vphp",-1.d0,"s2  ")
      call mult("s2  ","tg1 ","s1  ",222)
      call add("tb2 ",1.d0,"s1  ",1.d0,"tb2 ")
      call prt("s1  ","diag4b  ",3)
      call prtime("diag4b  ",0,1)
      return
      end
c ********************************************
      subroutine diag5a
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tvm2","vgb1","s1  ",222)
      call add("tb2 ",1.d0,"s1  ",-1.d0,"tb2 ")
      call prt("s1  ","-diag5a ",3)
      call prtime("diag5a  ",0,1)
      return
      end
c *********************************************************
      subroutine diag5b
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vva1","tgm2","s1  ",222)
      call add("tb2 ",1.d0,"s1  ",-1.d0,"tb2 ")
      call prt("s1  ","-diag5b ",3)
      call prtime("diag5b  ",0,1)
      return
      end
c *******************************************************
      subroutine diag6a
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  only t1v*t1g product contributes to vv0
      call mult("t1v ","hphp","s1  ",131)
      call mult("s1  ","t1g ","s2  ",311)
      call add("vv0 ",1.d0,"s2  ",-1.d0,"vv0 ")
      call prt("s2  ","-diag6a ",3)
      call prtime("diag6a  ",0,1)
      return
      end
c ********************************************************
      subroutine diag22b
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  incl 22b 25a 12a. only t1v*t1g contributions to vv0.
c  22b
      call mult("t1v ","hhhp","s1  ",131)
      call mult("t1g ","s1  ","s2  ",131)
      call mult("t1cr","s2  ","s1  ",131)
      call reorder("s1  ","s2  ",3412,5)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call prt("s2  ","diag22b ",3)
      call prtime("diag22b ",0,1)
c  25a
      call mult("t1c ","vd14","s1  ",131)
      call mult("t1v ","s1  ","s2  ",131)
      call mult("s2  ","t1g ","s1  ",311)
      call add("vv0 ",1.d0,"s1  ",-1.d0,"vv0 ")
      call prt("s1  ","-diag25a",3)
      call prtime("diag25a ",0,1)
c  12a
      call mult("t1g ","va4 ","s1  ",131)
      call reorder("tc  ","s2  ",1432,5)
      call mult("s2  ","s1  ","s3  ",222)
      call mult("t1v ","s3  ","s2  ",131)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call prt("s2  ","diag12a ",3)
      call prtime("diag12a ",0,1)
      return
      end
c **********************************************************
      subroutine diag6b_chp
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder("tcm ","s1  ",1423,5)
      call mult("vghp","s1  ","s2  ",222)
      call add("tb15",1.d0,"s2  ",-1.d0,"tb15")
      call prt("s2  ","-diag6b ",3)
      call prtime("diag6b  ",0,1)
      return
      end
c *************************************************************
      subroutine diag9_chp
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl diags 9 10a 20a 23a
c   9a 10a
      call add("tg1 ",-2.d0,"tgm2",1.d0,"s2  ")
      call add("va1 ",-2.d0,"va2 ",1.d0,"s3  ")
      call mult("s2  ","s3  ","s1  ",222)
c   20a
      call reorder("hphh","s3  ",1324,5)
      call add("vc1 ",-2.d0,"s3  ",1.d0,"s2  ")
      call mult("t1g ","s2  ","s3  ",131)
      call reorder("s3  ","s2  ",2134,5)
      call add("s1  ",1.d0,"s2  ",1.d0,"s1  ")
c   23a
      call add("vga1",2.d0,"ghpp",-1.d0,"s2  ")
      call mult("t1c ","s2  ","s3  ",131)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
c --------------------
      call mult("tv1 ","s1  ","s2  ",222)
      call add("tb2 ",1.d0,"s2  ",1.d0,"tb2 ")
      call prt("s2  ","diag9   ",3)
      call prtime("diag9   ",0,1)
      return
      end
c ************************************************************
      subroutine diag20b
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl 20b, 23b
      call add("vvc2",-2.d0,"vvc3",1.d0,"s2  ")
      call mult("s2  ","t1cr","s1  ",311)
      call reorder("s1  ","s3  ",1423,5)
c   diag 23b
      call add("vd12",2.d0,"vd14",-1.d0,"s2  ")
      call mult("t1v ","s2  ","s1  ",131)
      call add("s3  ",1.d0,"s1  ",1.d0,"s3  ")
c ---------------------------------
      call mult("s3  ","tg1 ","s2  ",222)
      call add("tb2 ",1.d0,"s2  ",1.d0,"tb2 ")
      call prt("s2  ","diag20b ",3)
      call prtime("diag20b ",0,1)
      return
      end
c ***************************************************
      subroutine diag21b
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  24b
      call mult("t1v ","vd12","s2  ",131)
c  21b
      call mult("t1cr","vvc2","s3  ",131)
      call reorder("s3  ","s1  ",2134,5)
      call add("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
c ---------------
      call reorder("tg  ","s2  ",4132,4)
      call mult("s1  ","s2  ","s3  ",222)
      call add("tb2 ",1.d0,"s3  ",1.d0,"tb2 ")
      call prt("s3  ","diag21b ",3)
      call prtime("diag21b ",0,3)
      return
      end
c *************************************************
      subroutine diag26
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  incl 26a 26b
c  26a
      call mult("t1c ","vhhp","s1  ",131)
      call reorder("s1  ","s3  ",2134,5)
c  26b
      call reorder("hphh","s2  ",1342,5)
      call mult("t1v ","s2  ","s1  ",131)
c ---------------------
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
      call mult("s1  ","tgm ","s2  ",222)
      call add("tb4 ",1.d0,"s2  ",1.d0,"tb4 ")
      call prt("s2  ","diag26  ",3)
      call prtime("diag26  ",0,1)
      return
      end
c *************************************************
      subroutine diag27a
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  incl 27a 29a
      call mult("tv  ","vd14","s1  ",222)
      call mult("tv  ","hp  ","s2  ",311)
      call add("s1  ",1.d0,"s2  ",1.d0,"s1  ")
      call mult("s1  ","t1g ","s2  ",311)
      call add("vv0 ",1.d0,"s2  ",-1.d0,"vv0 ")
      call prt("s2  ","-diag27a",3)
      call prtime("diag27a ",0,1)
      return
      end
c *************************************************
      subroutine diag2b
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tv  ","gp  ","s1  ",311)
      call add("vv0 ",1.d0,"s1  ",1.d0,"vv0 ")
      call prt("s1  ","diag2b  ",3)
      call prtime("diag2b  ",0,1)
      return
      end
c *************************************************
      subroutine diag3a
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vh  ","tg2 ","s1  ",131)
      call add("vv0 ",1.d0,"s1  ",-1.d0,"vv0 ")
      call prt("s1  ","-diag3a ",3)
      call prtime("diag3a  ",0,1)
      return
      end
c **************************************************
      subroutine diag8
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("hvhh","tgm ","s1  ",222)
      call add("tb1 ",1.d0,"s1  ",1.d0,"tb1 ")
      call prt("s1  ","diag8   ",3)
      call prtime("diag8   ",0,1)
      return
      end
c **************************************************
      subroutine dgt1p3
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder("hphh","s1  ",1243,5)
      call mult("t1g ","s1  ","s2  ",131)
      call mult("t1cr","s2  ","s1  ",131)
      call mult("t1v ","s1  ","s2  ",131)
      call reorder("s2  ","s1  ",1423,5)
      call add("vv0 ",1.d0,"s1  ",-2.d0,"vv0 ")
      call prt("s1  ","-.5t1**3",3)
      call prtime("dgt1p3  ",0,1)
      return
      end
c ******************************************************
      subroutine diag11
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl diags 11a, 10b
c  11a
      call mult("tgm2","va1 ","s1  ",222)
c  10b
      call add("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult("tg1 ","s2  ","s3  ",222)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
c ---------------------
      call mult("tvm2","s1  ","s2  ",222)
      call add("tb2 ",1.d0,"s2  ",1.d0,"tb2 ")
      call prt("s2  ","diag11  ",3)
      call prtime("diag11  ",0,1)
      return
      end
c ***************************************************
      subroutine diag24
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl diags 21a, 24a
c  24
      call mult("t1c ","vga1","s1  ",131)
c   21
      call mult("t1g ","vc1 ","s2  ",131)
      call reorder("s2  ","s3  ",2134,5)
      call add("s3  ",1.d0,"s1  ",-1.d0,"s3  ")
c ----------------------
      call mult("tvm2","s3  ","s2  ",222)
      call add("tb2 ",1.d0,"s2  ",1.d0,"tb2 ")
      call prt("s2  ","diag24  ",3)
      call prtime("diag24  ",0,1)
      return
      end
c ************************************************************
      subroutine diag13b
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl 13b 18b
      call add("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call reorder("tvm ","s3  ",1324,4)
      call mult("s3  ","s2  ","s1  ",113)
c   18
      call reorder("vhhp","s3  ",1324,4)
      call add("s3  ",-2.d0,"vhhp",1.d0,"s2  ")
      call mult("s2  ","t1c ","s3  ",202)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
c ----------------------
      call mult("s1  ","tg2 ","s2  ",131)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call prt("s2  ","diag13b ",3)
      call prtime("diag13b ",0,1)
      return
      end
c ******************************************************
      subroutine diag14a
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl 14a, 19a
      call add("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
      call mult("tgm ","s2  ","s3  ",113)
c   19
      call reorder("ghpp","s1  ",1324,5)
      call reorder("ghpp","s2  ",1423,5)
      call add("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult("t1c ","s1  ","s2  ",022)
      call add("s2  ",1.d0,"s3  ",1.d0,"s2  ")
c ------------------
      call mult("tv  ","s2  ","s1  ",311)
      call add("vv0 ",1.d0,"s1  ",1.d0,"vv0 ")
      call prt("s1  ","diag14a ",3)
      call prtime("diag14a ",0,1)
      return
      end
c ******************************************************
      subroutine diag22a
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  incl 22a, 25b
      call mult("t1g ","vvc3","s1  ",131)
      call reorder("s1  ","s3  ",2134,5)
c   diag 25b
      call mult("t1v ","ghpp","s2  ",131)
      call add("s3  ",1.d0,"s2  ",-1.d0,"s3  ")
c ----------------------------------
      call mult("s3  ","tc1 ","s1  ",222)
      call add("tb15",1.d0,"s1  ",1.d0,"tb15")
      call prt("s1  ","diag22a ",3)
      call prtime("diag22a ",0,1)
      return
      end
c ************************************************
      subroutine diag27b
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tv  ","vga1","s2  ",222)
      call mult("s2  ","t1cr","s1  ",311)
      call add("tb4 ",1.d0,"s1  ",-1.d0,"tb4 ")
      call prt("s1  ","-diag27b",3)
      call prtime("diag27b ",0,1)
      return
      end
c **************************************************
      subroutine diag15
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tvm ","va3 ","s1  ",222)
      call mult("s1  ","tgm ","s2  ",222)
      call add("tb4 ",1.d0,"s2  ",1.d0,"tb4 ")
      call prt("s2  ","diag15  ",3)
      call prtime("diag15  ",0,1)
      return
      end
c ************************************************************
      subroutine diag17
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  incl 17a 17b
      call mult("hvgh","t1cr","s2  ",311)
      call add("tb1 ",1.d0,"s2  ",-1.d0,"tb1 ")
      call prt("s2  ","diag17a ",3)
      call prtime("diag17a ",0,1)
c ----
      call mult("vhph","t1g ","s2  ",311)
      call add("vv0 ",1.d0,"s2  ",-1.d0,"vv0 ")
      call prt("s2  ","diag17b ",3)
      call prtime("diag17b ",0,1)
      return
      end
c **********************************************************
      subroutine diag27
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c   incl  27, 29
c   27
      call mult("tm6 ","hppp","s1  ",222)
      call reorder("s1  ","s2  ",1243,5)
c   29a
      call mult("thp ","hp  ","s1  ",311)
      call add("s1  ",-1.d0,"s2  ",-1.d0,"s1  ")
c ------------
      call mult("s1  ","t1cr","s2  ",311)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","diag27  ",3)
      call prtime("diag27  ",0,1)
      return
      end
c ********************************************
      subroutine diag28b
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("t1v ","hp  ","s2  ",111)
      call mult("s2  ","tg2 ","s1  ",131)
      call add("vv0 ",1.d0,"s1  ",1.d0,"vv0 ")
      call prt("s1  ","diag28b ",3)
      call prtime("diag28b ",0,1)
      return
      end
c ***************************************************
      subroutine diag16_chp
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  incl 16a 16b
c  16a
      call mult("t1v ","hpgp","s2  ",131)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call prt("s2  ","diag16a ",3)
      call prtime("diag16a ",0,1)
c  16b
      call mult("t1c ","vgpp","s2  ",131)
      call add("tb1 ",1.d0,"s2  ",1.d0,"tb1 ")
      call prt("s2  ","diag16b ",3)
      call prtime("diag16b ",0,1)
      return
      end
c ************************************************
      subroutine dgt1p4
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  -2*t1**4 compensates for overcounting.
      call mult("t1v ","hhpp","s3  ",131)
      call mult("t1c ","s3  ","s2  ",131)
      call mult("t1cr","s2  ","s3  ",131)
      call mult("t1g ","s3  ","s2  ",131)
      call reorder("s2  ","s3  ",4321,5)
      call add("vv0 ",1.d0,"s3  ",-2.d0,"vv0 ")
      call prt("s3  ","-.5t1**4",3)
      call prtime("dgt1p4  ",0,1)
      return
      end
c *************************************************
      subroutine dfoldhp
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

c  f2a
      call mult("thp ","tgf ","s2  ",311)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","f2a     ",3)
c  f3a
      call reorder("thpf","s1  ",1432,5)
      call mult("s1  ","t1gr","s2  ",311)
      call unfold_ccsd("s2  ","s1  ",3)
      call add("tb3 ",1.d0,"s1  ",1.d0,"tb3 ")
      call prt("s1  ","f3a     ",3)
c  f1a
      call reorder("thpf","s3  ",1423,5)
      call mult("s3  ","tm4 ","s2  ",222)
      call add("tb15",1.d0,"s2  ",1.d0,"tb15")
      call prt("s2  ","f1a     ",3)
c  f2b
      call reorder("thp ","s1  ",2341,5)
      call mult("tvf ","s1  ","s2  ",131)
      call add("tnew",1.d0,"s2  ",-1.d0,"tnew")
      call prt("s2  ","f2b     ",3)
c  f3b
      call reorder("t1v ","s1  ",1243,4)
      call mult("s3  ","s1  ","s2  ",311)
      call unfold_ccsd("s2  ","s1  ",3)
      call add("tb15",1.d0,"s1  ",-1.d0,"tb15")
      call prt("s1  ","-f3b    ",3)
      call prtime("dfoldhp ",0,1)
      return
      end
C  ********************************************************************
      subroutine dfold3
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  diag 1
      call reord7    ("t3p ","s1  ",231456)
      call reord6    ("s1  ","s2  ",124563)
      call mult6     ("s2  ","tvf ","s1  ",5)
      call reord6    ("s1  ","s2  ",126345)
      call reord7    ("s2  ","s1  ",312456)
      call add       ("s1  ",-1.d0,"t3nw",1.d0,"t3nw")
c  diag 2
      call reord6    ("t3p ","s1  ",123465)
      call mult6     ("s1  ","tgf ","s2  ",5)
      call reord6    ("s2  ","s1  ",123465)
      call add       ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
c  diag 3
      call t1t2      ("t1v ","tgr ","t3p ","s1  ",123456,1)
      call t1t2      ("t1gr","tv  ","s1  ","s2  ",123456,2)
      call reord7    ("s2  ","s1  ",231456)
      call reord6    ("s1  ","s2  ",124635)
      call reorder   ("thpf","s1  ",1432,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord7    ("s3  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",124536)
      call add       ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
c  diag 4
      call reorder   ("thpf","s1  ",1342,4)
      call reorder   ("tg  ","s2  ",2341,4)
      call mult6     ("s1  ","s2  ","s3  ",1)
      call reord6    ("s3  ","s1  ",125346)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",126345)
      call cmult     ("s1  ",-1.d0,"s2  ")
      call fold6     ("t3nw","s2  ",idum,2)
c  diag 5
      call reorder   ("thpf","s1  ",1243,4)
      call reorder   ("tv  ","s2  ",2341,4)
      call mult6     ("s1  ","s2  ","s3  ",1)
      call reord6    ("s3  ","s1  ",124536)
      call fold6     ("t3nw","s1  ",idum,2)
      return
      end
C  *******************************************************************
      subroutine wrtbox_chp
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      call wrtfile (41,0,"boxA")
      call wrtfile (41,1,"boxB")
      call wrtfile (41,2,"boxC")
      call wrtfile (41,3,"boxD")
      call wrtfile (41,4,"boxE")
      call wrtfile (41,5,"boxF")
      call wrtfile (41,6,"boxG")
      call wrtfile (41,7,"boxH")
      call wrtfile (41,8,"boxI")
      if (iopt(17).eq.3) then
        call wrtfile (42,0,"bx31")
        call wrtfile (42,1,"bx32")
        call wrtfile (42,2,"bx33")
      end if
      return
      end 
C  *******************************************************************
      subroutine rdbox
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call rdfile (20,0,"bx1 ",-1)
      call rdfile (20,1,"bx2 ",-1)
      call rdfile (20,2,"bx3 ",-1)
      call rdfile (20,3,"bx4 ",-1)
      call rdfile (20,4,"bx5 ",-1)
      call rdfile (31,8,"bx6a",-1)
      call rdfile (20,6,"bx7 ",-1)
      call rdfile (20,7,"bx8 ",-1)
      return
      end 
C  *******************************************************************
      subroutine box31 (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      if (iparm.eq.0) then
c  diag 1
        call reorder    ("vppg","bx31",1342,4)
c  diag 2
        call reorder    ("ppgh","s1  ",3142,4)
        call reorder    ("ppgh","s2  ",3241,4)
        call add        ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
        call mult       ("tv1 ","s1  ","s2  ",222)
        call add        ("s2  ",1.d0,"bx31",1.d0,"bx31")
c  diag 3
        call reorder    ("ppgh","s1  ",3142,4)
        call mult       ("tvm2","s1  ","s2  ",222)
        call add        ("s2  ",-1.d0,"bx31",1.d0,"bx31")
c  diag 5
        call reorder    ("vhhp","s1  ",1423,4)
        call mult       ("s1  ","tgm ","s2  ",222)
        call reorder    ("s2  ","s1  ",1432,4)
        call add        ("s1  ",1.d0,"bx31",1.d0,"bx31")
c  diag 8
        call reorder    ("vhpp","s1  ",1342,4)
        call mult       ("s1  ","t1g ","s2  ",311)
        call reorder    ("s2  ","s1  ",1243,4)
        call add        ("s1  ",-1.d0,"bx31",1.d0,"bx31")
c  diag 9
        call reorder    ("vphg","s1  ",1423,4)
        call mult       ("s1  ","t1cr","s2  ",311)
        call reorder    ("s2  ","s1  ",1423,4)
        call add        ("s1  ",-1.d0,"bx31",1.d0,"bx31")
      else if (iparm.eq.1) then
c  diag 4
        call reorder    ("tm  ","s1  ",1423,4)
        call mult       ("s1  ","vd7 ","s2  ",222)
        call reorder    ("s2  ","s1  ",1324,4)
        call move       ("s1  ","bx31")
c  diag 6
        call t1t2       ("t1gr","tv  ","t3p ","s1  ",123456,2)
        call cmult      ("t1c ",0.5d0,"s2  ")
        call t1t2       ("s2  ","thp ","s1  ","s3  ",123456,3)
        call reord6     ("s3  ","s1  ",124356)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",125634)
        call add        ("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
        call mult6      ("s1  ","s2  ","s3  ",2)
        call add        ("s3  ",1.d0,"bx31",-1.d0,"bx31")
c  diag 7
        call t1t2       ("t1v ","tgr ","t3p ","s1  ",123456,1)
        call t1t2       ("t1gr","tv  ","s1  ","s2  ",123456,2)
        call t1t2       ("t1c ","tm  ","s2  ","s1  ",123456,3)
        call reord6     ("s1  ","s2  ",126345)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",126534)
        call mult6      ("s2  ","va4 ","s1  ",2)
        call add        ("s1  ",1.d0,"bx31",1.d0,"bx31")
      end if
      return
      end
C  *******************************************************************
      subroutine box32 (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      if (iparm.eq.0) then
c  diag 1
        call reorder    ("hvgh","bx32",2134,4)
c  diag 2
        call reorder    ("vhhp","s1  ",1324,4)
        call add        ("s1  ",-2.d0,"vhhp",1.d0,"s1  ")
        call reorder    ("tg  ","s2  ",1342,4)
        call mult       ("s1  ","s2  ","s3  ",222)
        call reorder    ("s3  ","s1  ",1432,4)
        call add        ("s1  ",1.d0,"bx32",-1.d0,"bx32")
c  diag 3
        call reorder    ("vhhp","s1  ",1324,4)
        call reorder    ("tgm ","s2  ",1432,4)
        call mult       ("s1  ","s2  ","s3  ",222)
        call reorder    ("s3  ","s1  ",1432,4)
        call add        ("s1  ",1.d0,"bx32",1.d0,"bx32")
c  diag 5
        call reorder    ("ppgh","s1  ",3421,4)
        call mult       ("tvm ","s1  ","s2  ",222)
        call add        ("s2  ",-1.d0,"bx32",1.d0,"bx32")
c  diag 8
        call reorder    ("hpgh","s1  ",1342,4)
        call mult       ("t1v ","s1  ","s2  ",131)
        call add        ("s2  ",-1.d0,"bx32",1.d0,"bx32")
c  diag 9
        call reorder    ("vphg","s1  ",1432,4)
        call mult       ("s1  ","t1c ","s2  ",311)
        call reorder    ("s2  ","s1  ",1423,4)
        call add        ("s1  ",-1.d0,"bx32",1.d0,"bx32")
c  diag 10
        call mult       ("vhhh","t1g ","s1  ",311)
        call reorder    ("s1  ","s2  ",1243,4)
        call add        ("s2  ",1.d0,"bx32",1.d0,"bx32")
      else if (iparm.eq.1) then
c  diag 4
        call reorder    ("tm  ","s1  ",1423,4)
        call mult       ("s1  ","hhhp","s2  ",222)
        call reorder    ("s2  ","s1  ",1324,4)
        call move       ("s1  ","bx32")
c  diag 6
        call t1t2       ("t1v ","tgr ","t3p ","s1  ",123456,1)
        call cmult      ("t1c ",0.5d0,"s2  ")
        call t1t2       ("s2  ","thp ","s1  ","s3  ",123456,3)
        call reord6     ("s3  ","s1  ",125346)
        call add        ("hhpp",-2.d0,"va6 ",1.d0,"s2  ")
        call mult6      ("s1  ","s2  ","s3  ",2)
        call add        ("s3  ",1.d0,"bx32",1.d0,"bx32")
c  diag 7
        call t1t2       ("t1v ","tgr ","t3p ","s1  ",123456,1)
        call t1t2       ("t1gr","tv  ","s1  ","s2  ",123456,2)
        call t1t2       ("t1c ","tm  ","s2  ","s1  ",123456,3)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",125346)
        call mult6      ("s1  ","hhpp","s2  ",2)
        call add        ("s2  ",1.d0,"bx32",1.d0,"bx32")
      end if
      return
      end
C  *******************************************************************
      subroutine box33 (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      if (iparm.eq.0) then
c  diag 1
        call reorder    ("vphg","bx33",1432,4)
c  diag 2
        call reorder    ("ppgh","s1  ",3412,4)
        call mult       ("t1v ","s1  ","s2  ",131)
        call add        ("s2  ",-1.d0,"bx33",-1.d0,"bx33")
c  diag 3
        call mult       ("vvc3","t1g ","s1  ",311)
        call reorder    ("s1  ","s2  ",1423,4)
        call add        ("s2  ",1.d0,"bx33",1.d0,"bx33")
      else if (iparm.eq.1) then
c  diag 4
        call reorder    ("tm  ","s1  ",1423,4)
        call mult       ("s1  ","va2 ","s2  ",222)
        call move       ("s2  ","bx33")
      end if
      return
      end
C  *******************************************************************
