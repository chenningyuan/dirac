!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine chx
!---------------------------------------------------------------------
!
!   (4,0)
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use cint_mod
      use dat_mod
      use gen_mod
#include "implicit.h"
#include "priunit.h"
#include "csdtime.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)
#include "cfiles.h"
!     common /cint/ nint,idum1,vint(ndint),iint(ndint)
#include "cint.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /csdtime/ csdstim,csdltim
#include "csdtime.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
#endif

      character*256 filenm,filenm2,filenm3,filenm6,                     &
 filenm51,filenm53,filenm55,filenm60,filenm62,filenm33
      real*4 cputime
!     iopt(1) - print.
!     iopt(5) - print only elements .gt. 10**(-iopt(5))
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
      filenm53(1:1)=' '
      call getenv('tape53',filenm53)
      if (filenm53(1:1).eq.' ') filenm53='tape53'
      open(53,file=filenm53,form='UNFORMATTED')
      filenm55(1:1)=' '
      call getenv('tape55',filenm55)
      if (filenm55(1:1).eq.' ') filenm55='tape55'
      open(55,file=filenm55,form='UNFORMATTED')
      filenm60(1:1)=' '
      call getenv('tape60',filenm60)
      if (filenm60(1:1).eq.' ') filenm60='tape60'
      open(60,file=filenm60,form='UNFORMATTED')
      filenm62(1:1)=' '
      call getenv('tape62',filenm62)
      if (filenm62(1:1).eq.' ') filenm62='tape62'
      open(62,file=filenm62,form='UNFORMATTED')
      filenm33(1:1)=' '
      call getenv('tape33',filenm33)
      if (filenm33(1:1).eq.' ') filenm33='tape33'
      open(33,file=filenm33,form='UNFORMATTED')
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')
      write (6,*) '   $$$$$$$$$   chx (4,0)   $$$$$$$$$  '
      csdstim=cputime(real(0.))
      csdltim=csdstim
      call rdin
      call rdsym
      call prtime("rdin    ",0,0)
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      call reads  (51,0,"s1  ")
      call reorder ("s1  ","tc  ",3412,4)
      call reads  (51,1,"s2  ")
      call reorder ("s2  ","t1c ",1243,4)
      call rdfile (53,0,"thv ",-1)
      call rdfile (53,1,"t1hv",-1)
      call rdfile (53,2,"thvf",-1)
      call rdfile (55,0,"thu ",-1)
      call rdfile (55,1,"thuf",-1)
      call rdfile (60,0,"thwf",-1)
      if (iopt(17).ge.1) then
        call rdfile (53,3,"t31 ",-1)
        call rdfile (55,2,"t32 ",-1)
        call rdfile (60,1,"t33 ",-1)
      end if
      call t1sqr    ("t1c ","t1c ","tc  ","tcm ")
      call t1sqr    ("t1hv","t1c ","thv ","tvm ")
      call t1sqr    ("t1hv","t1hv","thu ","tum ")
      call ordv_chx 
      call ordta_chx
      call diags1_chx
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
        call diags2_chx
      end if
      if (iopt(17).eq.3) then
        call diags3_chx ! mi fix
      end if
      call sqr4    ("thxf")
      egr0=0.d0
      call eneru2 (egr0,1,2,max(1,iopt(1)))                         
      call eneru2 (egr0,2,2,max(1,iopt(1)))                         
      call eneru2 (egr0,3,2,max(1,iopt(1)))                         
      call eneru2 (egr0,4,2,max(1,iopt(1)))                         
      call prt6 ("thxf","thxf",0)
      call wrtfile (62,0,"thxf")
      write (6,*) nstrt(1,nflcur),'words used in igen,',                &
 nstrt(2,nflcur),' words used in gen'
      call prtime ("end ",0,0)
!      stop
      end
!  *******************************************************************
      subroutine ordv_chx 
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
#include "param_ccsd.h"
!     common /cint/ nint,idum1,vint(ndint),iint(ndint)
#include "cint.h"
!     common /dat/  norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1) 
#include "fock.h"
#endif

      call readv    (nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call seltyp   (3322,0,"pphh",4)    
      call seltyp   (2223,0010,"hhgp",4)
      call seltyp   (2222,0011,"hhgg",4)
      call seltyp   (2222,0111,"hggg",4)   
      call seltyp   (2332,1000,"gpph",4)
      call seltyp   (2232,1001,"ghpg",4)
      call seltyp   (2223,1010,"ghgp",4)
      call seltyp   (2233,1100,"ggpp",4)
      call seltyp   (2223,1110,"gggp",4)   
      call selt1    (23,0,"hp  ",4)
      call selt1    (32,0,"ph  ",4)
      call reorder  ("pphh","va1 ",1324,4)
      call reorder  ("pphh","va2 ",1423,4)
      call reorder  ("pphh","va4 ",3412,4)
      call reorder  ("pphh","va6 ",1243,4)
      call reorder  ("pphh","va8 ",3421,4)
      call prtime   ("readv   ",nint,0)
      return
      end
!  *******************************************************************
      subroutine ordta_chx
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call fold4    ("t1c ","t1c1",2,2)
      call fold4    ("thv ","s1  ",3,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","tv7 ",1432,4)
      call prtime   ("ordta   ",0,0)
      return
      end
!  *******************************************************************
      subroutine diags1_chx
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  GROUP I
!   diag 1
      call reorder  ("hhgg","s3  ",3412,4)
!   diag 2
      call fold4    ("tcm ","s1  ",3,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","s1  ",3412,4)
      call mult     ("s1  ","va4 ","s2  ",222)
      call add      ("s2  ",0.5d0,"s3  ",0.5d0,"s3  ")
!
      call fold4    ("thu ","s1  ",3,4)
      call mult6    ("s1  ","s3  ","s2  ",1)
      call reord6   ("s2  ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","boxZ",124356)
!
!   diag 4
      call fold4    ("thv ","s1  ",3,4)
      call reorder  ("s1  ","s2  ",1324,4)
      call add      ("va1 ",-2.d0,"va2 ",1.d0,"s1  ")
      call mult     ("s2  ","s1  ","s3  ",222)
!   diag 3
      call reorder  ("ghgp","s1  ",1342,4)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!   diag 5
      call fold4    ("tvm ","s1  ",4,4)
      call reorder  ("s1  ","s2  ",1423,4)
      call mult     ("s2  ","va1 ","s1  ",222)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!
      call reorder  ("s3  ","s1  ",1243,4)
      call mult6    ("tv7 ","s1  ","s2  ",1)
      call reord6   ("s2  ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124536)
      call add      ("s1  ",1.d0,"boxZ",1.d0,"boxZ")
!
!   diag 6
      call reorder  ("ghpg","s3  ",1423,4)
!   diag 7
      call fold4    ("tvm ","s1  ",4,4)
      call reorder  ("s1  ","s2  ",1423,4)
      call mult     ("s2  ","va2 ","s1  ",222)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",1.d0,"s3  ",-1.d0,"s3  ")
!   diag 8
      call reorder  ("hhgp","s1  ",3142,4)
      call mult     ("t1hv","s1  ","s2  ",131)
      call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!
      call mult6    ("tv7 ","s3  ","s1  ",1)
      call add      ("s1  ",1.d0,"boxZ",1.d0,"boxZ")
!
!   diag 9
      call clear6   ("s1  ",223222,110101)
      call fold4    ("thu ","s2  ",3,4)
      call t1t2     ("t1c1","s2  ","s1  ","s3  ",123456,3)
      call fold4    ("thv ","s1  ",3,4)
      call fold4    ("s1  ","s2  ",4,4)
      call t1t2     ("t1hv","s2  ","s3  ","s1  ",123456,2)
      call reord6   ("s1  ","s2  ",124635)
      call reorder  ("hhgp","s1  ",2341,4)
      call mult6    ("s2  ","s1  ","s3  ",3)
      call reord6   ("s3  ","s1  ",126345)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124356)
      call add      ("s1  ",1.d0,"boxZ",1.d0,"boxZ")
!   diag 10
      call clear6   ("s1  ",233222,100111)
      call fold4    ("thv ","s2  ",3,4)
      call fold4    ("s2  ","s3  ",4,4)
      call t1t2     ("t1c1","s3  ","s1  ","s2  ",123456,3)
      call t1t2     ("t1c1","s3  ","s2  ","s1  ",123456,2)
      call reord6   ("s1  ","s2  ",124356)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",125634)
      call reorder  ("gpph","s1  ",1432,4)
      call mult6    ("s2  ","s1  ","s3  ",3)
      call reord7   ("s3  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",123546)
      call add      ("s2  ",-1.d0,"boxZ",1.d0,"boxZ")
!
      call reord6   ("boxZ","s1  ",124356)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124356)
      call fold4    ("thu ","s2  ",3,4)
      call mult8    ("s1  ","s2  ","s3  ",1)
      call reordx   ("s3  ","thxf",12673485)
!  GROUP II
!   diag 1 
      call cmult    ("ggpp",0.5d0,"s3  ")
!   diag 2
      call mult     ("tum ","pphh","s1  ",222)
      call add      ("s1  ",0.5d0,"s3  ",1.d0,"s3  ")
!   diag 3
      call mult     ("t1hv","gpph","s1  ",131)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!
      call mult6    ("s3  ","tv7 ","s1  ",1)
      call reord6   ("s1  ","s2  ",124653)
      call mult8    ("s2  ","tv7 ","s1  ",1)
      call reordx   ("s1  ","s2  ",12367548)
      call add      ("s2  ",1.d0,"thxf",1.d0,"thxf")
      call prtime   ("diags1  ",0,0)
      return
      end
!  *******************************************************************
      subroutine diags2_chx
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  diag 1
      call fold6    ("t33 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reorder  ("hggg","s1  ",2341,4)
      call mult8    ("s2  ","s1  ","s3  ",1)
      call reordx   ("s3  ","s1  ",12364578)
      call add      ("s1  ",-0.5d0,"thxf",1.d0,"thxf")
!  diag 2
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call fold6    ("s2  ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",124563)
      call mult8    ("s2  ","gggp","s1  ",1)
      call reordx   ("s1  ","s2  ",12673485)
      call add      ("s2  ",0.5d0,"thxf",1.d0,"thxf")
      call prtime   ("diags2  ",0,0)
      return
      end 
!  *******************************************************************
      subroutine diags3_chx
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  GROUP I
!   diag 1
      call rdfile   (33,0,"boxA",-1)     
! g1 Apr 95
      call fold4    ("thv  ","s1  ",3,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","s1  ",1342,4)
      call mult     ("s1  ","hp  ","s2  ",311)
      call add      ("boxA",1.d0,"s2  ",-.5d0,"boxA")

!   diag 2
      call fold6    ("t31 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord6   ("s2  ","s1  ",124356)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",125634)
      call add      ("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","s3  ",2)
      call add      ("s3  ",1.d0,"boxA",1.d0,"boxA")
!   diag 3
      call fold6    ("t31 ","s1  ",5,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",126543)
      call mult6    ("s1  ","va4 ","s2  ",2)
      call add      ("s2  ",1.d0,"boxA",1.d0,"boxA")
!
      call fold6    ("t33 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call mult8    ("s2  ","boxA","s1  ",1)
      call reordx   ("s1  ","s2  ",12364587)
      call add      ("s2  ",0.5d0,"thxf",1.d0,"thxf")
      call prtime   ("I   ",0,0)
!  GROUP II
!   diag 1
      call rdfile   (33,1,"boxB",-1)     
! g1 Apr 95
      call fold4    ("thu  ","s1  ",3,4)
      call mult     ("s1  ","ph  ","s2  ",311)
      call add      ("boxB",1.d0,"s2  ",-.5d0,"boxB")

!   diag 2
      call t1t2     ("t1c ","thu ","t32 ","s1  ",123456,3)
      call fold6    ("s1  ","s2  ",4,1)
      call reord6   ("s2  ","s1  ",124356)
      call add      ("pphh",-2.d0,"va6 ",1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","s3  ",2)
      call add      ("s3  ",1.d0,"boxB",1.d0,"boxB")
!   diag 3
      call fold6    ("t32 ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",126354)
      call mult6    ("s2  ","pphh","s1  ",2)
      call add      ("s1  ",1.d0,"boxB",1.d0,"boxB")
!
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call fold6    ("s2  ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",124563)
      call mult8    ("s2  ","boxB","s1  ",1)
      call reordx   ("s1  ","s2  ",12673485)
      call add      ("s2  ",0.5d0,"thxf",1.d0,"thxf")
      call prtime   ("II  ",0,0)
!  GROUP III
!   diag 1
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord6   ("s2  ","s1  ",124536)
      call reorder  ("hhgp","s2  ",3142,4)
      call reorder  ("hhgp","s3  ",3241,4)
      call add      ("s2  ",-2.d0,"s3  ",1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","boxC",3)
!   diag 2
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124635)
      call reorder  ("hhgp","s2  ",3142,4)
      call mult6    ("s1  ","s2  ","s3  ",3)
      call add      ("s3  ",1.d0,"boxC",0.5d0,"boxC")
!   diag 3
      call reorder  ("hhgp","s2  ",3241,4)
      call mult6    ("s1  ","s2  ","s3  ",3)
      call reord6   ("s3  ","s1  ",123546)
      call add      ("s1  ",1.d0,"boxC",1.d0,"boxC")
!   diags 4,5
      call mult     ("t1hv","pphh","s1  ",131)
      call add      ("s1  ",1.d0,"gpph",-1.d0,"s1  ")
      call fold6    ("t31 ","s2  ",4,1)
      call fold6    ("s2  ","s3  ",5,1)
      call fold6    ("s3  ","s2  ",6,1)
      call reord6   ("s2  ","s3  ",124356)
      call reord7   ("s3  ","s2  ",132456)
      call reord6   ("s2  ","s3  ",125634)
      call reorder  ("s1  ","s2  ",1432,4)
      call mult6    ("s3  ","s2  ","s1  ",3)
      call reord6   ("s1  ","s2  ",125346)
      call reord7   ("s2  ","s1  ",132456)
      call add      ("s1  ",1.d0,"boxC",1.d0,"boxC")
!
      call fold4    ("thu ","s1  ",3,4)
      call mult8    ("boxC","s1  ","s2  ",1)
      call reordx   ("s2  ","s1  ",12673485)
      call add      ("s1  ",1.d0,"thxf",1.d0,"thxf")
      call prtime   ("III ",0,0)
!  GROUP IV
!   diag 1
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord6   ("s2  ","s1  ",124536)
      call reorder  ("gpph","s2  ",1324,4)
      call add      ("s2  ",2.d0,"gpph",-1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","boxD",3)
!   diag 2
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124635)
      call reorder  ("gpph","s2  ",1324,4)
      call mult6    ("s1  ","s2  ","s3  ",3)
      call add      ("s3  ",-1.d0,"boxD",0.5d0,"boxD")
!   diag 3
      call mult6    ("s1  ","gpph","s2  ",3)
      call reord6   ("s2  ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124536)
      call add      ("s1  ",-1.d0,"boxD",1.d0,"boxD")
!   diags 4,5
      call reorder  ("hhgp","s3  ",3412,4)
      call fold4    ("t1c ","s1  ",2,2)
      call reorder  ("s1  ","s2  ",1243,4)
      call mult     ("s2  ","va4 ","s1  ",131)
      call reorder  ("s1  ","s2  ",1432,4)
      call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
      call fold6    ("t33 ","s1  ",4,1)
      call mult6    ("s1  ","s3  ","s2  ",3)
      call reord6   ("s2  ","s1  ",124536)
      call add      ("s1  ",1.d0,"boxD",1.d0,"boxD")
!
      call reord6   ("boxD","s1  ",125346)
      call mult8    ("s1  ","tv7 ","s2  ",1)
      call reordx   ("s2  ","s1  ",12364578)
      call add      ("s1  ",1.d0,"thxf",1.d0,"thxf")
      call prtime   ("IV  ",0,0)
!  GROUP V
!   diag 1
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord6   ("s2  ","s1  ",124536)
      call add      ("va1 ",2.d0,"va2 ",-1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","s3  ",3)
      call mult8    ("s3  ","s1  ","s2  ",5)
      call reordx   ("s2  ","s3  ",12563478)
      call add      ("s3  ",0.25d0,"thxf",1.d0,"thxf")
!   diag 2
      call add      ("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","s3  ",3)
      call t1t2     ("t1hv","thv ","t32 ","s1  ",123456,2)
      call t1t2     ("t1c ","thu ","s1  ","s2  ",123456,3)
      call fold6    ("s2  ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124635)
      call mult8    ("s3  ","s1  ","s2  ",5)
      call reordx   ("s2  ","s1  ",12563478)
      call add      ("s1  ",0.5d0,"thxf",1.d0,"thxf")
!   diag 3
      call cmult    ("t1hv",2.d0,"s1  ")
      call t1t2     ("s1  ","thv ","t32 ","s2  ",123456,2)
      call cmult    ("t1c ",2.d0,"s1  ")
      call t1t2     ("s1  ","thu ","s2  ","s3  ",123456,3)
      call fold6    ("s3  ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s3  ",124635)
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124635)
      call mult6    ("s1  ","va1 ","s2  ",3)
      call mult8    ("s2  ","s3  ","s1  ",5)
      call reordx   ("s1  ","s2  ",12563478)
      call add      ("s2  ",0.5d0,"thxf",1.d0,"thxf")
!   diag 4
      call cmult    ("t1hv",2.d0,"s1  ")
      call t1t2     ("s1  ","thv ","t32 ","s2  ",123456,1)
      call cmult    ("t1c ",2.d0,"s1  ")
      call t1t2     ("s1  ","thu ","s2  ","s3  ",123456,3)
      call fold6    ("s3  ","s1  ",5,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s3  ",125634)
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124635)
      call mult6    ("s1  ","va2 ","s2  ",3)
      call mult8    ("s2  ","s3  ","s1  ",5)
      call reordx   ("s1  ","s2  ",12563847)
      call add      ("s2  ",0.5d0,"thxf",1.d0,"thxf")
!   diag 5
      call fold6    ("t33 ","s1  ",4,1)
      call mult6    ("s1  ","pphh","s2  ",3)
      call fold6    ("t31 ","s1  ",4,1)
      call fold6    ("s1  ","s3  ",5,1)
      call fold6    ("s3  ","s1  ",6,1)
      call reord6   ("s1  ","s3  ",124356)
      call reord7   ("s3  ","s1  ",132456)
      call reord6   ("s1  ","s3  ",125643)
      call mult8    ("s2  ","s3  ","s1  ",5)
      call reordx   ("s1  ","s2  ",12354876)
      call add      ("s2  ",0.5d0,"thxf",1.d0,"thxf")
      call prtime   ("diags3  ",0,0)
      return
      end 
!  *******************************************************************
