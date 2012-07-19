!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine m12b
!---------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------
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
#include "dat.h"
!     common /csdtime/ csdstim,csdltim
#include "csdtime.h"
#endif

      character*256 filenm,filenm6,                                     &
 filenm51,filenm52,filenm53,filenm55,filenm56,filenm57,filenm58
      real*4 cputime

!     iopt(1) - print.
!     iopt(5) - print only elements .gt. 10**(-iopt(5))
!     iopt(8) - fock.
!     iopt(10) .ne.0 - do not write tape9
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
      filenm52(1:1)=' '
      call getenv('tape52',filenm52)
      if (filenm52(1:1).eq.' ') filenm52='tape52'
      open(52,file=filenm52,form='UNFORMATTED')
      filenm53(1:1)=' '
      call getenv('tape53',filenm53)
      if (filenm53(1:1).eq.' ') filenm53='tape53'
      open(53,file=filenm53,form='UNFORMATTED')
      filenm55(1:1)=' '
      call getenv('tape55',filenm55)
      if (filenm55(1:1).eq.' ') filenm55='tape55'
      open(55,file=filenm55,form='UNFORMATTED')
      filenm56(1:1)=' '
      call getenv('tape56',filenm56)
      if (filenm56(1:1).eq.' ') filenm56='tape56'
      open(56,file=filenm56,form='UNFORMATTED')
      filenm57(1:1)=' '
      call getenv('tape57',filenm57)
      if (filenm57(1:1).eq.' ') filenm57='tape57'
      open(57,file=filenm57,form='UNFORMATTED')
      filenm58(1:1)=' '
      call getenv('tape58',filenm58)
      if (filenm58(1:1).eq.' ') filenm58='tape58'
      open(58,file=filenm58,form='UNFORMATTED')
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')
      write (*,*) '   $$$$$$$$$   m12b   $$$$$$$$$  '
      csdstim = cputime(real(0.))
      csdltim = csdstim
      call rdin
      call rdsym
      call prtime ("rdin    ",0,0)
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      call reads   (51,0,"s1  ")
      call reorder ("s1  ","tc  ",3412,4)
      call reads   (51,1,"s1  ")
      call reorder ("s1  ","t1c ",1243,4)
      call rdfile  (52,0,"tv  ",-1)
      call rdfile  (52,1,"t1v ",-1)
      call rdfile  (52,2,"tvf ",-1)
      call reads   (56,0,"s1  ")
      call reorder ("s1  ","tp  ",4321,4)
      call reads   (56,1,"s1  ")
      call reorder ("s1  ","tpf ",4321,4)
      call reads   (57,0,"s1  ")
      call reorder ("s1  ","tq  ",3412,4)
      call reads   (57,1,"s1  ")
      call reorder ("s1  ","tqf ",3412,4)
      call rdfile  (58,1,"t12 ",-1)
      call rdfile  (53,0,"tg  ",-1)
      call rdfile  (53,1,"t1g ",-1)    
      call rdfile  (53,2,"tgf ",-1)
      call rdfile  (55,0,"thu ",-1)
      call rdfile  (55,1,"thuf",-1)
      if (iopt(17).gt.0) then
        call rdfile  (52,3,"t3v ",-1)
        call reads   (56,2,"s2  ")
        call reord6  ("s2  ","s1  ",125346)
        call reord7  ("s1  ","s2  ",312456)
        call reord6  ("s2  ","s1  ",125346)
        call reord7  ("s1  ","s2  ",132456)
        call reord6  ("s2  ","t3p ",126435)
        call reads   (57,2,"s2  ")
        call reord6  ("s2  ","s1  ",124356)
        call reord7  ("s1  ","s2  ",312456)
        call reord6  ("s2  ","s1  ",125346)
        call reord7  ("s1  ","s2  ",132456)
        call reord6  ("s2  ","t3q ",126345)
        call rdfile  (58,2,"t12p",-1)
        call rdfile  (58,3,"t12q",-1)
        call rdfile  (53,3,"t3g ",-1)
        call rdfile  (55,2,"t3hu",-1) 
      end if
      call t1sqr  ("t1c ","t1c ","tc  ","tm  ")
      call t1sqr ("t1g ","t1g ","thu ","thum")
      call reorder ("t1v ","s1  ",1243,4)
      call t1sqr ("t1g ","s1  ","tp  ","tpm ")
      call ordv_m12b
      call clear6 ("ttqf",223232,111111)
      call groupa 
      call groupb 
      call groupc 
      call groupd 
      call groupe 
      call groupf 
      call groupg 
      call grouph 
      call groupi 
      if (iopt(17).gt.0) then
        call groupq
      end if
      call ord1
      call clear6 ("ttpf",223223,111111)
      call groupj 
      call groupk 
      call groupl 
      call groupm 
      call groupn 
      call groupo 
      if (iopt(17).gt.0) then
        call groupp
      end if
      call ord2
      call move ("tgf ","thvf")
      call prt6 ("ttpf","ttpf",0)
      call prt6 ("ttqf","ttqf",0)
      egr0=0.d0
      call eneru2 (egr0,1,2,max(1,iopt(1)))
      call eneru2 (egr0,12,2,max(1,iopt(1)))
      write (6,*) nstrt(1,nflcur),'words used in igen,',                &
 nstrt(2,nflcur),' words used in gen'
      call prtime ("end ",0,0)
!     stop
      end
!  ********************************************************************
      subroutine ordv_m12b
!---------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cint_mod
      use dat_mod
      use fock_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     parameter (ndt1=1000,ndint=3000000)
#include "param_ccsd.h"
!     common /cint/ nint,idum1,vint(ndint),iint(ndint)
#include "cint.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)
#include "fock.h"
#endif

      call readv    (nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call prtime   ("readv   ",nint,0)
      call selt1    (23,0,"hp  ",4)   
      call seltyp   (3322,0,"pphh",4)
      call seltyp   (2322,0111,"hvgg",4)
      call seltyp   (2223,0010,"hhgp",4)
      call seltyp   (2232,0001,"hhpg",4)
      call seltyp   (2222,0011,"hhgg",4)
      call seltyp   (2333,0100,"hvpp",4)
      call seltyp   (2323,0110,"hvgp",4)
      call seltyp   (2332,0101,"hvpg",4)
      call seltyp   (2233,1000,"ghpp",4)
      call seltyp   (2233,0001,"hhpv",4)
      call seltyp   (2233,1010,"ghvp",4)
      call seltyp   (2233,1100,"ggpp",4)
      call reorder  ("pphh","vm1 ",3124,4)
      call reorder  ("hhgp","vm2 ",1342,4)
      call reorder  ("pphh","vm4 ",3412,4)
      call reorder  ("hhpg","vm5 ",1432,4)
      call reorder  ("pphh","vm8 ",3214,4)
      call reorder  ("pphh","vm10",3214,4)
      call fold4    ("t1c ","t1cb",1,2)
      call fold4    ("tg  ","s1  ",2,4)
      call fold4    ("s1  ","tgc ",3,4)
      call fold4    ("tg  ","tgg ",2,4)
      call fold4    ("tv  ","tvc ",2,4)
      call fold4    ("tvc ","tvd ",3,4)
      call fold4    ("thu ","thua",3,4)
      call fold4    ("tq  ","tqa ",2,4)
      call fold4    ("tv  ","s1  ",2,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","tr1 ",1243,4)
      call fold4    ("tg  ","s1  ",3,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","tr2 ",1342,4)
      call fold4    ("tq  ","s1  ",4,4)
      call reorder  ("s1  ","tr3 ",1342,4)
      call fold4    ("tp  ","s1  ",3,4)
      call reorder  ("s1  ","tr4 ",1342,4)
      call fold4    ("tg  ","s1  ",2,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","tr5 ",1243,4)
      call fold4    ("tc  ","s1  ",1,4)
      call fold4    ("s1  ","s2  ",3,4)
      call fold4    ("s2  ","s1  ",4,4)
      call reorder  ("s1  ","tr6 ",1432,4)
      call fold4    ("tc  ","s1  ",1,4)
      call fold4    ("s1  ","s2  ",3,4)
      call reorder  ("s2  ","tr7 ",1324,4)
      call fold4    ("tm  ","s1  ",1,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","tr8 ",1423,4)
      call fold4    ("tm  ","s1  ",3,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","tr9 ",3412,4)
      call fold4    ("tc  ","s1  ",1,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","tr10",1423,4)
      call fold4    ("t1c ","s1  ",2,2)
      call reorder  ("s1  ","tr11",1243,4)
      call reorder  ("tq  ","tr14",1342,4)
      call reorder  ("tpm ","tr15",1432,4)
      call fold4    ("tg  ","s1  ",3,4)
      call reorder  ("s1  ","tr16",1342,4)
      call fold4    ("tg  ","s1  ",4,4)
      call reorder  ("s1  ","tr17",1432,4)
      call fold4    ("tv  ","s1  ",4,4)
      call reorder  ("s1  ","tr18",1423,4)
      call reorder  ("tp  ","tr19",1432,4)
      call fold4    ("tp  ","s1  ",2,4)
      call reorder  ("s1  ","tr20",1243,4)
      call reorder  ("tvd ","tr21",2134,4)
      call fold4    ("tv  ","s1  ",3,4)
      call reorder  ("s1  ","tr22",1324,4)
      if (iopt(17).gt.0) then
        call seltyp    (2333,1100,"gvpp",4)
        call seltyp    (2223,0011,"hhgv",4)
        call seltyp    (2233,1001,"ghpv",4)
        call seltyp    (2232,1001,"ghpg",4)
        call seltyp    (2223,1010,"ghgp",4)
        call seltyp    (3233,1010,"vhvp",4)
        call seltyp    (3233,1001,"vhpv",4)
        call seltyp    (3223,1010,"vhgp",4)
        call seltyp    (3232,1001,"vhpg",4)
        call reorder   ("ghpg","vm11",1423,4)
        call reorder   ("vhvp","vm12",1324,4)
        call reorder   ("ghgp","vm13",1324,4)
        call reorder   ("ghvp","vm14",1324,4)
        call reorder   ("ghpv","vm15",1423,4)
        call reorder   ("hhgv","vm16",3412,4)
        call reorder   ("vhpg","vm17",1432,4)
        call reorder   ("vhgp","vm18",1342,4)
      end if
      call prtime   ("ordv",0,0)
      return
      end
!  ********************************************************************
      subroutine groupa
!---------------------------------------------------------------------
!
!
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
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif
!
!
!     GROUP A
!
!  diags hb1,2
      call add      ("vm1 ",2.d0,"vm8 ",-1.d0,"s1  ")
      call mult     ("s1  ","t1c ","s2  ",202)
      call add      ("s2  ",1.d0,"hp  ",1.d0,"hb1 ")
      call reorder  ("hb1 ","s3  ",1243,4)
!  diag a.1
      call mult     ("thua","s3  ","s1  ",311)
      call mult6    ("s1  ","tr1 ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tb1 ")
!  diag a.2
      call mult     ("tr2 ","hb1 ","s1  ",311)
      call mult6    ("s1  ","tqa ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tb2 ")
!  diag a.3
      call mult     ("tr3 ","hb1 ","s1  ",311)
      call mult6    ("tgc ","s1  ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tb3 ")
!  diag a.4
      call mult     ("tr4 ","hb1 ","s1  ",311)
      call mult6    ("s1  ","tr5 ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tb4 ")
!  diag a.5
      call mult     ("tr5 ","s3  ","s1  ",311)
      call mult6    ("s1  ","tr3 ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tb5 ")
      call prtime   ("group a ",1,0)
      return
      end
!  ********************************************************************
      subroutine groupb
!---------------------------------------------------------------------
!
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
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif
!
!
!     GROUP B
!
!     diag hbg 
!     diag hbg-a
      call move     ("hvgg","hbg ")
!  diag hbg-b
      call mult     ("hp  ","tr6 ","s1  ",131)
      !call add      ("s1  ",1.d0,"hbg ",1.d0,"hbg ",1.d0)
      call add      ("s1  ",1.d0,"hbg ",1.d0,"hbg ")
!  diag hbg-c
      call reorder  ("hhgp","s2  ",2341,4)
      call add      ("vm2 ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult     ("s1  ","tr7 ","s3  ",222)
!  diag hbg-d
      call mult     ("vm2 ","tr8 ","s1  ",222)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s1  ")
      call reorder  ("s1  ","s2  ",1324,4)
      call add      ("s2  ",1.d0,"hbg ",1.d0,"hbg ")
!  diag hbg-e
      call mult     ("vm5 ","tr8 ","s1  ",222)
      call reorder  ("s1  ","s2  ",1342,4)
      call add      ("s2  ",-1.d0,"hbg ",1.d0,"hbg ")
!  diag hbg-f
      call mult     ("hvpp","tr9 ","s1  ",222)
      call add      ("s1  ",1.d0,"hbg ",1.d0,"hbg ")
!  diag hbg-j
      call add      ("vm1 ",2.d0,"vm8 ",-1.d0,"s1  ")
      call mult     ("s1  ","t1c ","s2  ",202)
      call mult     ("s2  ","tr6 ","s1  ",131)
      call add      ("s1  ",1.d0,"hbg ",1.d0,"hbg ")
!  diag hbg-h,k,l
      call mult     ("vm1 ","tr10","s3  ",222)
      call add      ("vm1 ",2.d0,"vm10 ",-1.d0,"s1  ")
      call mult     ("s1  ","tr7 ","s2  ",222)
      call add      ("s2  ",1.d0,"s3  ",-1.d0,"s3  ")
      call reorder  ("s3  ","s1  ",1342,4)
      call reorder  ("hvpg","s2  ",1243,4)
      call add      ("s2  ",1.d0,"s1  ",1.d0,"s2  ")
      call mult     ("s2  ","tr11","s1  ",311)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",1.d0,"hbg ",1.d0,"hbg ")
!  diag hbg-i,m
      call mult     ("vm10","tr10","s1  ",222)
      call reorder  ("s1  ","s2  ",1342,4)
      call add      ("hvgp",1.d0,"s2  ",-1.d0,"s1  ")
      call mult     ("s1  ","tr11","s2  ",311)
      call add      ("s2  ",1.d0,"hbg ",1.d0,"hbg ")
!  diag hbg-g,n
      call mult     ("vm4 ","tr9 ","s2  ",222)
      call reorder  ("s2  ","s1  ",1342,4)
      call reorder  ("hhgg","s2  ",1342,4)
      call add      ("s2  ",1.d0,"s1  ",1.d0,"s2  ")
      call mult     ("s2  ","t1cb","s1  ",311)
      call reorder  ("s1  ","s2  ",1423,4)
      call add      ("s2  ",-1.d0,"hbg ",1.d0,"hbg ")
!  diag b.1 
      call reorder  ("hbg ","s1  ",2341,4)
      call reorder  ("t12 ","s2  ",2134,4)
      call mult6    ("s2  ","s1  ","s3  ",1)
      call cmult    ("s3  ",-1.d0,"tb6 ")
      call prtime   ("group b ",2,0)
      return
      end
!  ********************************************************************
      subroutine groupc
!---------------------------------------------------------------------
!
!
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
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!
!     GROUP C
!
!  diag hbu-a
      call move     ("hhgg","hbu ")
!  diag hbu-b
      call mult     ("vm4 ","tr9 ","s1  ",222)
      call add      ("s1  ",1.d0,"hbu ",1.d0,"hbu ")
!  diag hbu-c,d
      call reorder  ("hhpg","s2  ",1243,4)
      call mult     ("s2  ","tr11","s1  ",311)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",1.d0,"hbu ",1.d0,"hbu ")
      call reorder  ("s1  ","s2  ",2134,4)
      call add      ("s2  ",1.d0,"hbu ",1.d0,"hbu ")
!  diag c.1,c.2
      call fold4    ("hbu ","s1  ",1,4)
      call reorder  ("s1  ","boxc",1342,4)
!
      call reorder  ("hbu ","s1  ",3421,4)
      call mult     ("t1g ","s1  ","s2  ",131)
      call add      ("s2  ",-1.d0,"boxc",1.d0,"boxc")
!
      call mult6    ("boxc","tqa ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb2 ",1.d0,"tb2 ")
      call prtime   ("group c ",3,0)
      return
      end
!  ********************************************************************
      subroutine groupd
!---------------------------------------------------------------------
!
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
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!
!     GROUP D
!
!  diag hbp-a
      call move     ("hvgp","hbp ")
!  diag hbp-b,d
      call mult     ("vm8 ","tr8 ","s3  ",222)
      call reorder  ("hvpp","s2  ",1423,4)
      call mult     ("s2  ","tr11","s1  ",311)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
      call reorder  ("s3  ","s1  ",1342,4)
      call add      ("s1  ",-1.d0,"hbp ",1.d0,"hbp ")
!  diag hbp-c
      call mult     ("vm2 ","t1cb","s1  ",311)
      call reorder  ("s1  ","s2  ",1423,4)
      call add      ("s2  ",-1.d0,"hbp ",1.d0,"hbp ")
!  diag d.1,d.2   
      call reorder  ("hbp ","s1  ",2341,4)
      call mult     ("t1g ","s1  ","boxd",131)
!
      call fold4    ("hbp ","s1  ",1,4)
      call add      ("s1  ",1.d0,"boxd",-1.d0,"boxd")
!
      call mult6    ("boxd","tr3 ","s2  ",1)
      call add      ("s2  ",1.d0,"tb3 ",1.d0,"tb3 ")
      call prtime   ("group d ",4,0)
      return
      end
!  ********************************************************************
      subroutine groupe
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!
!     GROUP E
!
!  diag hbq-a
      call move     ("hvpg","hbq ")
!  diag hbq-b,c
      call add      ("vm1 ",2.d0,"vm8 ",-1.d0,"s1  ")
      call mult     ("s1  ","tr7 ","s2  ",222)
      call mult     ("vm1 ","tr8 ","s3  ",222)
      call add      ("s2  ",1.d0,"s3  ",-1.d0,"s3  ")
      call reorder  ("s3  ","s1  ",1324,4)
      call add      ("s1  ",1.d0,"hbq ",1.d0,"hbq ")
!  diag hbq-d
      call mult     ("vm5 ","t1cb","s1  ",311)
      call reorder  ("s1  ","s2  ",1432,4)
      call add      ("s2  ",-1.d0,"hbq ",1.d0,"hbq ")
!  diag hbq-e
      call mult     ("hvpp","tr11","s1  ",311)
      call add      ("s1  ",1.d0,"hbq ",1.d0,"hbq ")
!  diag e.1,e.4
      call reorder  ("hbq ","s1  ",4213,4)
      call mult     ("t1v ","s1  ","s3  ",131)
      call fold4    ("hbq ","s1  ",3,4)
      call reorder  ("s1  ","s2  ",3421,4)
      call add      ("s3  ",1.d0,"s2  ",-1.d0,"s1  ")
      call mult6    ("thua","s1  ","s2  ",1)
      call add      ("s2  ",1.d0,"tb1 ",1.d0,"tb1 ")
!  diag e.2,e.3,e.5,e.6
      call reorder  ("hbq ","s1  ",2431,4)
      call mult     ("t1g ","s1  ","s3  ",131)
      call fold4    ("hbq ","s1  ",1,4)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",1.d0,"s3  ",-1.d0,"s3  ")
      call mult6    ("s3  ","tr3 ","s1  ",1)
      call add      ("s1  ",1.d0,"tb5 ",1.d0,"tb5 ")
      call mult6    ("tr4 ","s3  ","s1  ",1)
      call add      ("s1  ",1.d0,"tb4 ",1.d0,"tb4 ")
      call prtime   ("group e ",5,0)
      return
      end
!  ********************************************************************
      subroutine groupf
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!
!     GROUP F
!
!  diag hbs-a
      call move     ("hhgp","hbs ")
!  diag hbs-b
      call reorder  ("pphh","s2  ",3421,4)
      call mult     ("s2  ","tr11","s1  ",311)
      call reorder  ("s1  ","s2  ",1243,4) 
      call add      ("hbs ",1.d0,"s2  ",1.d0,"hbs ")
      call reorder  ("hbs ","hbs1",3124,4)
      call reorder  ("hbs ","hbs2",3214,4)
!  diags f.1,7,12,20,21
      call fold4    ("hbs ","s1  ",2,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","s3  ",2431,4)
!
      call reorder  ("s1  ","s2  ",2134,4)
      call mult     ("s2  ","t1v ","s1  ",311)
      call reorder  ("s1  ","s2  ",1432,4)
      call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!
      call fold4    ("hbs ","s1  ",4,4)
      call reorder  ("s1  ","s2  ",4312,4)
      call mult     ("t1g ","s2  ","s1  ",131)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!
      call add      ("hbs1",2.d0,"hbs2",-1.d0,"s1  ")
      call mult     ("tr14","s1  ","s2  ",222)
      call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!    
      call mult     ("tr15","hbs1","s1  ",222)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!
      call mult6    ("tgc ","s3  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb3 ",1.d0,"tb3 ")
!  diags f.2,3,11,17,18
      call mult6    ("tr5 ","s3  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb5 ",1.d0,"tb5 ")
!  diags f.4,6,10,19
      call fold4    ("hbs ","s1  ",1,4)
      call fold4    ("s1  ","s2  ",4,4)
      call  reorder ("s2  ","boxy",1342,4)
!
      call mult     ("s1  ","t1v ","s2  ",311)
      call reorder  ("s2  ","s3  ",1342,4)
      call add      ("s3  ",1.d0,"boxy",1.d0,"boxy")
!
      call fold4    ("hbs ","s1  ",4,4)
      call reorder  ("s1  ","s2  ",3421,4)
      call mult     ("t1g ","s2  ","s1  ",131)
      call add      ("s1  ",-1.d0,"boxy",1.d0,"boxy")
!
      call mult     ("tr15","hbs2","s1  ",222)
      call reorder  ("s1  ","s2  ",1324,4)
      call add      ("s2  ",-1.d0,"boxy",1.d0,"boxy")
!
      call mult6    ("boxy","tr5 ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb4 ",1.d0,"tb4 ")
!  diags f.5,9,8,25
      call fold4    ("hbs ","s1  ",1,4)
      call fold4    ("s1  ","boxx",2,4)
!
      call fold4    ("hbs ","s1  ",1,4)
      call reorder  ("s1  ","s2  ",1432,4)
      call mult     ("s2  ","t1g ","s1  ",311)
      call reorder  ("s1  ","s2  ",1432,4)
      call add      ("s2  ",-1.d0,"boxx",1.d0,"boxx")
!
      call fold4    ("hbs ","s1  ",2,4)
      call reorder  ("s1  ","s2  ",2341,4)
      call mult     ("t1g ","s2  ","s1  ",131)
      call add      ("s1  ",-1.d0,"boxx",1.d0,"boxx")
!
      call reorder  ("hbs ","s1  ",3412,4)
      call mult     ("thum","s1  ","s2  ",222)
      call add      ("s2  ",1.d0,"boxx",1.d0,"boxx")
!
      call mult6    ("boxx","tr1 ","s1  ",1)
      call add      ("s1  ",1.d0,"tb1 ",1.d0,"tb1 ")
!  diags f.13,14,15
      call add      ("hbs1",2.d0,"hbs2",-1.d0,"s3  ")
      call mult     ("tr16","s3  ","boxw",222)
!
      call mult     ("tr17","hbs1","s1  ",222)
      call add      ("s1  ",-1.d0,"boxw",1.d0,"boxw")
!
      call mult     ("tr17","hbs2","s1  ",222)
      call reorder  ("s1  ","s2  ",1324,4)
      call add      ("s2  ",-1.d0,"boxw",1.d0,"boxw")
!
      call mult6    ("boxw","tqa ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb2 ",1.d0,"tb2 ")
!  diags f.16
      call mult     ("tr18","hbs2","s1  ",222)
      call reorder  ("s1  ","s2  ",1324,4)
      call mult6    ("thua","s2  ","s1  ",1)        
      call add      ("s1  ",1.d0,"tb1 ",1.d0,"tb1 ")
!  diags f.23,22
      call reorder  ("hbs ","s1  ",3421,4)
      call mult     ("tgg ","s1  ","s3  ",222)
      call mult6    ("tr4 ","s3  ","s2  ",1)
      call add      ("s2  ",1.d0,"tb4 ",1.d0,"tb4 ")
!
      call mult6    ("s3  ","tr3 ","s2  ",1)
      call add      ("s2  ",1.d0,"tb5 ",1.d0,"tb5 ")
!  diag f.24
      call reorder  ("hbs ","s1  ",3412,4)
      call mult     ("tgg ","s1  ","s2  ",222)
      call mult6    ("s2  ","tr3 ","s1  ",1)
      call add      ("s1  ",1.d0,"tb3 ",1.d0,"tb3 ")
      call prtime   ("group f ",6,0)
      return
      end
!  ********************************************************************
      subroutine groupg
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!
!     GROUP G
!
!     diag hbt-a
      call move     ("hvpp","hbt ")
!     diag hbt-b
      call mult     ("vm8 ","t1cb","s1  ",311)
      call reorder  ("s1  ","s2  ",1432,4)
      call add      ("s2  ",-1.d0,"hbt ",1.d0,"hbt ")
      call reorder  ("hbt ","hbt1",2314,4)
      call reorder  ("hbt ","hbt2",2413,4)
!   diags g.1,2,3,9,10
      call fold4    ("hbt ","s1  ",1,4)
      call fold4    ("s1  ","s3  ",3,4)
!
      call reorder  ("s1  ","s2  ",1243,4)
      call mult     ("s2  ","t1v ","s1  ",311)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!
      call fold4    ("hbt ","s1  ",3,4)
      call reorder  ("s1  ","s2  ",2341,4)
      call mult     ("t1g ","s2  ","s1  ",131)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!
      call add      ("hbt2",2.d0,"hbt1",-1.d0,"s1  ")
      call mult     ("tr14","s1  ","s2  ",222)
      call reorder  ("s2  ","s1  ",1324,4)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!
      call mult     ("tr15","hbt2","s1  ",222)
      call reorder  ("s1  ","s2  ",1324,4)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
!
      call mult6    ("tr2 ","s3  ","s1  ",1)
      call add      ("s1  ",1.d0,"tb2 ",1.d0,"tb2 ")
!   diag g.5,6
      call add      ("hbt2",2.d0,"hbt1",-1.d0,"s1  ")
      call mult     ("tr16","s1  ","boxg",222)
!
      call mult     ("tr17","hbt2","s1  ",222)
      call add      ("s1  ",-1.d0,"boxg",1.d0,"boxg")
!
      call reorder  ("boxg","s2  ",1324,4)
      call mult6    ("s2  ","tr3 ","s1  ",1)
      call add      ("s1  ",1.d0,"tb3 ",1.d0,"tb3 ")
!   diag g.7,8
      call mult     ("tr17","hbt1","s1  ",222)
      call reorder  ("s1  ","s3  ",1324,4)
      call cmult    ("s3  ",-1.d0,"s3  ")
!
      call mult6    ("s3  ","tr3 ","s1  ",1)
      call add      ("s1  ",1.d0,"tb5 ",1.d0,"tb5 ")
!
      call mult6    ("tr4 ","s3  ","s1  ",1)
      call add      ("s1  ",1.d0,"tb4 ",1.d0,"tb4 ")
!   diag g.4
      call reorder  ("hbt ","s1  ",2134,4)
      call mult     ("tvc ","s1  ","s2  ",222)
      call mult6    ("thua","s2  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb1 ",1.d0,"tb1 ")
      call prtime   ("group g ",7,0)
      return
      end
!  ********************************************************************
      subroutine grouph
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!
!     GROUP H
!
      call mult     ("t1g ","vm8 ","s1  ",131)
      call add      ("s1  ",-1.d0,"ghpp",1.d0,"boxp")
      call reorder  ("boxp","bp1 ",1324,4)
      call reorder  ("boxp","bp2 ",1423,4)
      call reorder  ("pphh","s2  ",3412,4)
      call mult     ("s2  ","t1v ","s1  ",311)
      call add      ("s1  ",1.d0,"hhpv",1.d0,"boxh")
!   diag h.1,2,16,19
      call mult     ("tvc ","boxp","s1  ",222)
      call reorder  ("s1  ","s3  ",3124,4)
!
      call reorder  ("boxh","s1  ",4123,4)
      call mult     ("tr17","s1  ","s2  ",222)
      call reorder  ("s2  ","s1  ",1324,4)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!
      call mult6    ("tgc ","s3  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb3 ",1.d0,"tb3 ")
!   diag h.15,17,20,23
      call mult6    ("tr5 ","s3  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb5 ",1.d0,"tb5 ")
!   diag h.13,14,18,21,22,24
      call reorder  ("boxh","s1  ",4213,4)
      call reorder  ("boxh","s2  ",4123,4)
      call add      ("s1  ",2.d0,"s2  ",-1.d0,"s2  ")
      call mult     ("tr16","s2  ","boxu",222)
!
      call mult     ("tr17","s1  ","s2  ",222)
      call add      ("s2  ",-1.d0,"boxu",1.d0,"boxu")
!
      call reorder  ("boxp","s1  ",1243,4)
      call mult     ("tvc ","s1  ","s2  ",222)
      call reorder  ("s2  ","s1  ",3214,4)
      call add      ("s1  ",1.d0,"boxu",1.d0,"boxu")
!
      call mult6    ("boxu","tr5 ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb4 ",1.d0,"tb4 ")
!   diag h.3,4,5,10,11,12
      call add      ("bp1 ",2.d0,"bp2 ",-1.d0,"s1  ")
      call mult     ("tr16","s1  ","boxv",222)
!
      call mult     ("tr17","bp1 ","s1  ",222)
      call add      ("s1  ",-1.d0,"boxv",1.d0,"boxv")
!
      call mult     ("tr17","bp2 ","s1  ",222)
      call reorder  ("s1  ","s2  ",3214,4)
      call add      ("s2  ",-1.d0,"boxv",1.d0,"boxv")
!
      call reorder  ("boxv","s2  ",1324,4)
      call mult6    ("s2  ","tr1 ","s1  ",1)
      call add      ("s1  ",1.d0,"tb1 ",1.d0,"tb1 ")
!   diag h.6,7,8,9
      call reorder  ("boxh","s1  ",4321,4)
      call mult     ("tgg ","s1  ","s3  ",222)
!
      call mult     ("tr18","bp2 ","s1  ",222)
      call reorder  ("s1  ","s2  ",3214,4)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
      call mult6    ("tr2 ","s3  ","s1  ",1)
      call add      ("s1  ",1.d0,"tb2 ",1.d0,"tb2 ")
      call prtime   ("group h ",8,0)
      return
      end
!  ********************************************************************
      subroutine groupi
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!
!     GROUP I
!
!   diag i.1,3
      call fold4    ("bp1 ","s1  ",2,4)
      call fold4    ("s1  ","s3  ",3,4)
!   diag i.2
      call reorder  ("ghvp","s1  ",1342,4)
      call mult     ("t1g ","s1  ","s2  ",131)
      call reorder  ("s2  ","s1  ",2314,4)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!   diag i.4
      call mult     ("ggpp","t1v ","s1  ",311)
      call reorder  ("s1  ","s2  ",2413,4)
      call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!   diag i.5,9
      call add      ("bp1 ",2.d0,"bp2 ",-1.d0,"s1  ")
      call mult     ("tr14","s1  ","s2  ",222)
      call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")     
!   diag i.6
      call reorder  ("ghpp","s1  ",1324,4)
      call mult     ("tr15","s1  ","s2  ",222)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
!   diag i.7,11
      call mult     ("tr15","bp2 ","s1  ",222)
      call reorder  ("s1  ","s2  ",3214,4)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
!   diag i.8,12
      call reorder  ("boxh","s1  ",4312,4)
      call mult     ("thum","s1  ","s2  ",222)
      call reorder  ("s2  ","s1  ",2314,4)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!   diag i.10
      call mult     ("t1g ","vm8 ","s1  ",131)
      call reorder  ("s1  ","s2  ",1324,4)
      call mult     ("tr19","s2  ","s1  ",222)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!
      call reorder  ("s3  ","s2  ",3124,4)
      call mult6    ("s2  ","tr6 ","s1  ",1)
      call add      ("s1  ",1.d0,"tb6 ",1.d0,"tb6 ")
      call prtime   ("group i ",9,0)
      return
      end
!  ********************************************************************
      subroutine groupq
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!  diag q1
      call fold6    ("t3v ","s1  ",2,1)
      call fold6    ("s1  ","s2  ",3,1)
      call fold6    ("s2  ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",123654)
      call mult6    ("ggpp","s2  ","s1  ",3)
      call reord6   ("s1  ","s2  ",126435)
      call add      ("s2  ",1.d0,"ttqf",1.d0,1.d0,"ttqf")
!  diag q2
      call fold6    ("t3p ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124356)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",125634)
      call mult6    ("s1  ","gvpp","s2  ",3)
      call reord6   ("s2  ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",126345)
      call add      ("s1  ",1.d0,"ttqf",1.d0,"ttqf")
!  diag q3
      call fold6    ("t3q ","s1  ",5,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124356)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",125643)
      call mult6    ("gvpp","s1  ","s2  ",3)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",123645)
      call add      ("s2  ",1.d0,"ttqf",1.d0,"ttqf")
!  diag q4
      call fold6    ("t3hu","s1  ",3,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",123645)
      call mult6    ("s1  ","vm16","s2  ",3)
      call reord6   ("s2  ","s1  ",123564)
      call add      ("s1  ",1.d0,"ttqf",1.d0,"ttqf")
!  diag q5
      call fold6    ("t12q","s1  ",3,1)
      call reord6   ("s1  ","s2  ",123546)
      call reorder  ("hhgg","s1  ",3412,4)
      call mult6    ("s2  ","s1  ","s3  ",3)
      call reord6   ("s3  ","s1  ",123546)
      call add      ("s1  ",1.d0,"ttqf",1.d0,"ttqf")
!  diag q6
      call fold6    ("t3hu","s1  ",3,1)
      call fold6    ("s1  ","s2  ",4,1)
      call reorder  ("hhgv","s1  ",4321,4)
      call mult6    ("s2  ","s1  ","s3  ",3)
      call add      ("s3  ",1.d0,"ttqf",1.d0,"ttqf")
!  diag q7
      call fold6    ("t3q ","s1  ",2,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord6   ("s2  ","s1  ",124563)
      call add      ("vm13",2.d0,"vm11",-1.d0,"s2  ")
      call mult6    ("s2  ","s1  ","s3  ",3)
      call reord7   ("s3  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",124356)
      call add      ("s2  ",1.d0,"ttqf",1.d0,"ttqf")
!  diag q8
      call fold6    ("t12q","s1  ",4,1)
      call reord6   ("s1  ","s2  ",124536)
      call add      ("vm18",2.d0,"vm17",-1.d0,"s1  ")
      call mult6    ("s2  ","s1  ","s3  ",3)
      call reord6   ("s3  ","s1  ",125346)
      call add      ("s1  ",1.d0,"ttqf",1.d0,"ttqf")
!  diag q9
      call fold6    ("t3g ","s1  ",2,1)
      call fold6    ("s1  ","s2  ",4,1)
      call fold6    ("s2  ","s1  ",5,1)
      call reord6   ("s1  ","s2  ",124563)
      call add      ("vm14",2.d0,"vm15",-1.d0,"s1  ")
      call mult6    ("s2  ","s1  ","s3  ",3)
      call reord6   ("s3  ","s1  ",125364)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("s2  ",1.d0,"ttqf",1.d0,"ttqf")
!  diag q10
      call fold6    ("t12q","s1  ",6,1)
      call reord6   ("s1  ","s2  ",125634)
      call mult6    ("s2  ","vm17","s1  ",3)
      call reord6   ("s1  ","s2  ",125634)
      call add      ("s2  ",-1.d0,"ttqf",1.d0,"ttqf")
!  diag q11
      call fold6    ("t3hu","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124653)
      call reorder  ("vhpv","s2  ",4123,4)
      call mult6    ("s1  ","s2  ","s3  ",3)
      call reord6   ("s3  ","s1  ",126354)
      call add      ("s1  ",-1.d0,"ttqf",1.d0,"ttqf")
!  diag q12
      call fold6    ("t3q ","s1  ",2,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124653)
      call mult6    ("vm11","s1  ","s2  ",3)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",124653)
      call add      ("s2  ",-1.d0,"ttqf",1.d0,"ttqf")
!  diag q13
      call fold6    ("t3p ","s1  ",3,1)
      call fold6    ("s1  ","s2  ",4,1)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",124563)
      call mult6    ("s2  ","vm11","s1  ",3)
      call reord6   ("s1  ","s2  ",125346)
      call reord7   ("s2  ","s1  ",132456)
      call add      ("s1  ",-1.d0,"ttqf",1.d0,"ttqf")
!  diag q14
      call fold6    ("t3g ","s1  ",3,1)
      call fold6    ("s1  ","s2  ",5,1)
      call fold6    ("s2  ","s1  ",6,1)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",125643)
      call mult6    ("vm15","s1  ","s2  ",3)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",124536)
      call add      ("s2  ",-1.d0,"ttqf",1.d0,"ttqf")
!  diag q15
      call fold6    ("t3p ","s1  ",3,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",125643)
      call mult6    ("s2  ","vm11","s1  ",3)
      call reord6   ("s1  ","s2  ",125634)
      call reord7   ("s2  ","s1  ",132456)
      call add      ("s1  ",-1.d0,"ttqf",1.d0,"ttqf")
!  diag q16
      call fold6    ("t12q","s1  ",6,1)
      call reord6   ("s1  ","s2  ",126534)
      call mult6    ("s2  ","vm18","s1  ",3)
      call reord6   ("s1  ","s2  ",125346)
      call add      ("s2  ",-1.d0,"ttqf",1.d0,"ttqf")
!  diag q17    
      call fold6    ("t12p","s1  ",4,1)
      !call reord6   ("s1  ","s2  ",124635,4)
      call reord6   ("s1  ","s2  ",124635)
      call mult6    ("s2  ","vm18","s1  ",3)
      call reord6   ("s1  ","s2  ",125346)
      call add      ("s2  ",-1.d0,"ttqf",1.d0,"ttqf")
!  diag q18
      call fold6    ("t3q ","s1  ",3,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",124563)
      call mult6    ("vm13","s2  ","s1  ",3)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124356)
      call add      ("s1  ",-1.d0,"ttqf",1.d0,"ttqf")
!  diag q19
      call fold6    ("t3g ","s1  ",3,1)
      call fold6    ("s1  ","s2  ",4,1)
      call fold6    ("s2  ","s1  ",5,1)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124563)
      call mult6    ("s1  ","vm14","s2  ",3)
      call reord6   ("s2  ","s1  ",125364)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("s2  ",-1.d0,"ttqf",1.d0,"ttqf")
!  diag q20
      call fold6    ("t3p ","s1  ",3,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",125643)
      call mult6    ("vm13","s2  ","s1  ",3)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124356)
      call add      ("s1  ",-1.d0,"ttqf",1.d0,"ttqf")
!  diag q21
      call fold6    ("t3g ","s1  ",3,1)
      call fold6    ("s1  ","s2  ",5,1)
      call fold6    ("s2  ","s1  ",6,1)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",125643)
      call mult6    ("s1  ","vm14","s2  ",3)
      call reord6   ("s2  ","s1  ",125364)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("s2  ",-1.d0,"ttqf",1.d0,"ttqf")
      call prtime   ("group q ",0,0)
      return
      end
!  *********************************************************************
      subroutine ord1   
!---------------------------------------------------------------------
!
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
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

      call reord6   ("tb1 ","s1  ",126345)
      call add      ("s1  ",1.d0,"ttqf",1.d0,"ttqf")
      call reord6   ("tb2 ","s1  ",124563)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124356)
      call add      ("s1  ",1.d0,"ttqf",1.d0,"ttqf")
      call reord6   ("tb3 ","s1  ",124356)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("s2  ",1.d0,"ttqf",1.d0,"ttqf")
      call reord6   ("tb4 ","s1  ",124536)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124356)
      call add      ("s1  ",1.d0,"ttqf",1.d0,"ttqf")
      call reord6   ("tb5 ","s1  ",124653)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("s2  ",1.d0,"ttqf",1.d0,"ttqf")
      call reord6   ("tb6 ","s1  ",124536)
      call add      ("s1  ",1.d0,"ttqf",1.d0,"ttqf")
      call prtime   ("ord1",0,0)
      return
      end
!  *********************************************************************
      subroutine groupj
!---------------------------------------------------------------------
!
!
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
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

! 
!     GROUP J
!
!   diag j.1
      call mult     ("tvd ","hb1 ","s1  ",311)
      call reorder  ("s1  ","s2  ",2134,4)
      call mult6    ("thua","s2  ","s1  ",1)
      call cmult    ("s1  ",-1.d0,"tba ")               
!   diag j.2
      call mult     ("tr2 ","hb1 ","s1  ",311)
      call mult6    ("s1  ","tr20","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tbd ")                 
!   diag j.3
      call mult     ("tr4 ","hb1 ","s1  ",311)
      call mult6    ("s1  ","tgc ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tbb ")                 
      call prtime   ("group j ",1,0)
      return
      end
!  ********************************************************************
      subroutine groupk
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!
!     GROUP K
!
      call mult6    ("boxc","tr20","s1  ",1)
      call add      ("s1  ",-1.d0,"tbd ",1.d0,"tbd ")
      call prtime   ("group k ",2,0)
      return
      end
!  ********************************************************************
      subroutine groupl
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!
!     GROUP L
!
!   diag l.2,4
      call mult6    ("boxd","tr4 ","tbc ",1)
!   diag l.1,3
      call fold4    ("hbp ","s1  ",4,4)
      call reorder  ("s1  ","s3  ",3421,4)
!
      call reorder  ("hbp ","s1  ",2314,4)
      call mult     ("t1v ","s1  ","s2  ",131)
      call reorder  ("s2  ","s1  ",3124,4)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!
      call mult6    ("thua","s3  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tba ",1.d0,"tba ")
      call prtime   ("group l ",3,0)
      return
      end
!  ********************************************************************
      subroutine groupm
!---------------------------------------------------------------------
!
!
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
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!
!     GROUP M
!
!   diag m.1,3,6,12
      call mult6    ("tgc ","boxy","s1  ",1)
      call add      ("s1  ",-1.d0,"tbc ",1.d0,"tbc ")
!   diag m.2,4,14,15
      call mult6    ("boxx","tr21","s1  ",1)
      call add      ("s1  ",1.d0,"tba ",1.d0,"tba ")
!   diag m.7,8,9
      call mult6    ("boxw","tr20","s1  ",1)
      call add      ("s1  ",-1.d0,"tbd ",1.d0,"tbd ")
!   diag m.10,11
      call add      ("hbs1",2.d0,"hbs2",-1.d0,"s1  ")
      call mult     ("tr22","s1  ","s3  ",222)
!
      call mult     ("tr18","hbs1","s1  ",222)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!
      call reorder  ("s3  ","s2  ",3124,4)
      call mult6    ("thua","s2  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tba ",1.d0,"tba ")
!   diag m.13
      call reorder  ("hbs ","s1  ",3412,4)
      call mult     ("tgg ","s1  ","s2  ",222)
      call mult6    ("tr4 ","s2  ","s1  ",1)
      call add      ("s1  ",1.d0,"tbb ",1.d0,"tbb ")
      call prtime   ("group m ",4,0)
      return
      end
!  ********************************************************************
      subroutine groupn
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!  
!     GROUP N
!
!   diag n.1,2,3,7
      call fold4    ("hbt ","s1  ",1,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","s3  ",4213,4)
!
      call reorder  ("s1  ","s2  ",2134,4)
      call mult     ("t1v ","s2  ","s1  ",131)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!
      call fold4    ("hbt ","s1  ",4,4)
      call reorder  ("s1  ","s2  ",4231,4)
      call mult     ("s2  ","t1g ","s1  ",311)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
!
      call mult     ("tr15","hbt1","s1  ",222)
      call reorder  ("s1  ","s2  ",2314,4)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
!
      call reorder  ("s3  ","s2  ",3214,4)
      call mult6    ("tr2 ","s2  ","s1  ",1)
      call add      ("s1  ",1.d0,"tbd ",1.d0,"tbd ")
!   diag n.5,6
      call reorder  ("boxg","s1  ",1324,4)
      call mult6    ("tr4 ","s1  ","s2  ",1)
      call add      ("s2  ",1.d0,"tbb ",1.d0,"tbb ")
!   diag n.4
      call reorder  ("hbt ","s1  ",2143,4)
      call mult     ("tvc ","s1  ","s2  ",222)
      call reorder  ("s2  ","s1  ",2134,4)
      call mult6    ("thua","s1  ","s2  ",1)
      call add      ("s2  ",-1.d0,"tba ",1.d0,"tba ")
      call prtime   ("group n ",5,0)
      return
      end
!  ********************************************************************
      subroutine groupo
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!
!     GROUP O
!
!   diag o.1,2,3,16,17,18
      call mult6    ("tgc ","boxu","s1  ",1)
      call add      ("s1  ",-1.d0,"tbc ",1.d0,"tbc ")
!   diag o.4,5,6,13,14,15
      call mult6    ("boxv","tr21","s1  ",1)
      call add      ("s1  ",1.d0,"tba ",1.d0,"tba ")
!   diag o.7,8,9,10,11,12
      call add      ("bp1 ",2.d0,"bp2 ",-1.d0,"s1  ")
      call mult     ("tr22","s1  ","s3  ",222)
!
      call mult     ("tr18","bp1 ","s1  ",222)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!
      call reorder  ("boxh","s1  ",3412,4)
      call mult     ("tgg ","s1  ","s2  ",222)
      call reorder  ("s2  ","s1  ",4213,4)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!
      call reorder  ("s3  ","s2  ",3214,4)
      call mult6    ("tr2 ","s2  ","s1  ",1)
      call add      ("s1  ",1.d0,"tbd ",1.d0,"tbd ")
      call prtime   ("group o ",6,0)
      return
      end
!  ********************************************************************
      subroutine groupp
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!
!  diag p1
      call fold6    ("t3v ","s1  ",2,1)
      call fold6    ("s1  ","s2  ",3,1)
      call fold6    ("s2  ","s1  ",4,1)
      call reorder  ("ggpp","s2  ",1243,4)
      call mult6    ("s2  ","s1  ","s3  ",3)
      call reord6   ("s3  ","s1  ",126543)
      call add      ("s1  ",0.5d0,"ttpf",1.d0,"ttpf")
!  diag p2
      call fold6    ("t3p ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",125643)
      call mult6    ("s1  ","gvpp","s2  ",3)
      call reord6   ("s2  ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",126453)
      call add      ("s1  ",1.d0,"ttpf",1.d0,"ttpf")
!  diag p3
      call fold6    ("t12p","s1  ",3,1)
      call reord6   ("s1  ","s2  ",123645)
      call reorder  ("hhgg","s1  ",3412,4)
      call mult6    ("s2  ","s1  ","s3  ",3)
      call reord6   ("s3  ","s1  ",123564)
      call add      ("s1  ",0.5d0,"ttpf",1.d0,"ttpf")
!  diag p4
      call fold6    ("t3hu","s1  ",3,1)
      call fold6    ("s1  ","s2  ",4,1)
      call mult6    ("s2  ","vm16","s1  ",3)
      call add      ("s1  ",1.d0,"ttpf",1.d0,"ttpf")
!  diag p5
      call fold6    ("t12p","s1  ",4,1)
      call reord6   ("s1  ","s2  ",124653)
      call reorder  ("hvgp","s1  ",3214,4)
      call mult6    ("s2  ","s1  ","s3  ",3)
      call reord6   ("s3  ","s1  ",126354)
      call add      ("s1  ",-1.d0,"ttpf",1.d0,"ttpf")
!  diag p6
      call fold6    ("t3g ","s1  ",2,1)
      call fold6    ("s1  ","s2  ",4,1)
      call fold6    ("s2  ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",126453)
      call mult6    ("s2  ","vm15","s1  ",3)
      call reord6   ("s1  ","s2  ",125436)    
      call reord7   ("s2  ","s1  ",132456)
      call add      ("s1  ",-1.d0,"ttpf",1.d0,"ttpf")
!  diag p7
      call fold6    ("t3p ","s1  ",2,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",125643)
      call mult6    ("vm11","s1  ","s2  ",3)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",124635)
      call add      ("s2  ",-1.d0,"ttpf",1.d0,"ttpf")
!  diag p8
      call fold6    ("t3p ","s1  ",2,1)
      call fold6    ("s1  ","s2  ",4,1)
      call reord6   ("s2  ","s1  ",124563)
      call add      ("vm13",2.d0,"vm11",-1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","s3  ",3)
      call reord6   ("s3  ","s1  ",125364)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("s2  ",1.d0,"ttpf",1.d0,"ttpf")
!  diag p9
      call fold6    ("t3hu","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord6   ("s2  ","s1  ",124563)
      call reorder  ("vhpv","s2  ",1423,4)
      call add      ("vm12",2.d0,"s2  ",-1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","s3  ",3)
      call reord6   ("s3  ","s1  ",125346)
      call add      ("s1  ",0.5d0,"ttpf",1.d0,"ttpf")
!  diag p10
      call fold6    ("t3hu","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124653)
      call mult6    ("s1  ","vm12","s2  ",3)
      call reord6   ("s2  ","s1  ",125346)
      call add      ("s1  ",-1.d0,"ttpf",1.d0,"ttpf")
!  diag p11
      call fold6    ("t3p ","s1  ",3,1)
      call fold6    ("s1  ","s2  ",4,1)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",124563)
      call mult6    ("s2  ","vm13","s1  ",3)
      call reord6   ("s1  ","s2  ",125364)
      call reord7   ("s2  ","s1  ",132456)
      call add      ("s1  ",-1.d0,"ttpf",1.d0,"ttpf")
!  diag p12
      call fold6    ("t3p ","s1  ",2,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",125643)
      call mult6    ("vm13","s1  ","s2  ",3)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",124365)
      call add      ("s2  ",-1.d0,"ttpf",1.d0,"ttpf")
!
      call prtime   ("group p ",0,0)
      return
      end 
!  *********************************************************************
      subroutine ord2  
!---------------------------------------------------------------------
!
!
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
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
#endif

!
      call reord6   ("tba ","s1  ",126345)
      call add      ("s1  ",1.d0,"ttpf",1.d0,"ttpf")
      call reord6   ("tbb ","s1  ",124563)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124356)
      call add      ("s1  ",1.d0,"ttpf",1.d0,"ttpf")
      call reord6   ("tbc ","s1  ",124356)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("s2  ",1.d0,"ttpf",1.d0,"ttpf")
      call reord6   ("tbd ","s1  ",124536)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124356)
      call add      ("s1  ",1.d0,"ttpf",1.d0,"ttpf")
      call reord6   ("ttpf","s1  ",123546)
      call reord7   ("s1  ","s2  ",213456)
      call add      ("ttpf",1.d0,"s2  ",1.d0,"ttpf")
      call prtime   ("ord2",0,0)
      return
      end
!  ********************************************************************
!  ********************************************************************
