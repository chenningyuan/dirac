!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine csx
!---------------------------------------------------------------------
!
! (0,4)
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
 filenm51,filenm52,filenm54,filenm59,filenm61,filenm23

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
      filenm52(1:1)=' '
      call getenv('tape52',filenm52)
      if (filenm52(1:1).eq.' ') filenm52='tape52'
      open(52,file=filenm52,form='UNFORMATTED')
      filenm54(1:1)=' '
      call getenv('tape54',filenm54)
      if (filenm54(1:1).eq.' ') filenm54='tape54'
      open(54,file=filenm54,form='UNFORMATTED')
      filenm59(1:1)=' '
      call getenv('tape59',filenm59)
      if (filenm59(1:1).eq.' ') filenm59='tape59'
      open(59,file=filenm59,form='UNFORMATTED')
      filenm61(1:1)=' '
      call getenv('tape61',filenm61)
      if (filenm61(1:1).eq.' ') filenm61='tape61'
      open(61,file=filenm61,form='UNFORMATTED')
      filenm23(1:1)=' '
      call getenv('tape23',filenm23)
      if (filenm23(1:1).eq.' ') filenm23='tape23'
      open(23,file=filenm23,form='UNFORMATTED')
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')
      write (6,*) '   $$$$$$$$$   csx   $$$$$$$$$  '
      csdstim=cputime(real(0.))
      csdltim=csdstim
      call rdin
      call rdsym
      call prtime("rdin    ",0,0)
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      call rdfile (51,0,"tc  ",-1)
      call rdfile (51,1,"t1c ",-1)
      call rdfile (52,0,"tv  ",-1)
      call rdfile (52,1,"t1v ",-1)
      call rdfile (52,2,"tvf ",-1)
      call rdfile (54,0,"tu  ",-1)
      call rdfile (54,1,"tuf ",-1)
      call rdfile (59,0,"twf ",-1)
      if (iopt(17).ge.1) then
        call rdfile (52,3,"t31 ",-1)
        call rdfile (54,2,"t32 ",-1)
        call rdfile (59,1,"t33 ",-1)
      end if
      call t1sqr    ("t1c ","t1c ","tc  ","tcm ")
      call t1sqr    ("t1v ","t1c ","tv  ","tvm ")
      call t1sqr    ("t1v ","t1v ","tu  ","tum ")
      call ordv_csx 
      call ordta_csx
      call diags1_csx
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
        call diags2_csx
      end if
      if (iopt(17).eq.3) then
        call diags3_csx
      end if
      call sqr4    ("txf ")
      egr0=0.d0
      call eneru2 (egr0,1,3,max(1,iopt(1)))                         
      call eneru2 (egr0,2,3,max(1,iopt(1)))                         
      call eneru2 (egr0,3,3,max(1,iopt(1)))                         
      call eneru2 (egr0,4,3,max(1,iopt(1)))                         
      call prt6 ("txf ","txf ",0)
      call wrtfile (61,0,"txf ")
      write (6,*) nstrt(1,nflcur),'words used in igen,',                &
 nstrt(2,nflcur),' words used in gen'
      call prtime ("end ",0,0)
!      stop
      end
!  *******************************************************************
      subroutine ordv_csx 
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

      call seltyp   (2233,0,"hhpp",4)    
      call seltyp   (3332,0010,"ppvh",4)
      call seltyp   (3333,0011,"ppvv",4)
      call seltyp   (3333,0111,"pvvv",4)   
      call seltyp   (3223,1000,"vhhp",4)
      call seltyp   (3323,1001,"vphv",4)
      call seltyp   (3332,1010,"vpvh",4)
      call seltyp   (3322,1100,"vvhh",4)
      call seltyp   (3332,1110,"vvvh",4)   
      call selt1    (32,0,"ph  ",4)
      call selt1    (23,0,"hp  ",4)
      call reorder  ("hhpp","va1 ",1324,4)
      call reorder  ("hhpp","va2 ",1423,4)
      call reorder  ("hhpp","va4 ",3412,4)
      call reorder  ("hhpp","va6 ",1243,4)
      call reorder  ("hhpp","va8 ",3421,4)
      call prtime   ("readv   ",nint,0)
      return
      end
!  *******************************************************************
      subroutine ordta_csx
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call fold4    ("t1c ","t1c1",2,2)
      call fold4    ("tv  ","s1  ",3,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","tv7 ",1432,4)
      call prtime   ("ordta   ",0,0)
      return
      end
!  *******************************************************************
      subroutine diags1_csx
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  GROUP I
!   diag 1
      call reorder  ("ppvv","s3  ",3412,4)
!   diag 2
      call fold4    ("tcm ","s1  ",3,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","s1  ",3412,4)
      call mult     ("s1  ","va4 ","s2  ",222)
      call add      ("s2  ",0.5d0,"s3  ",0.5d0,"s3  ")
!
      call fold4    ("tu  ","s1  ",3,4)
      call mult6    ("s1  ","s3  ","s2  ",1)
      call reord6   ("s2  ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","boxZ",124356)
!
!   diag 4
      call fold4    ("tv  ","s1  ",3,4)
      call reorder  ("s1  ","s2  ",1324,4)
      call add      ("va1 ",-2.d0,"va2 ",1.d0,"s1  ")
      call mult     ("s2  ","s1  ","s3  ",222)
!   diag 3
      call reorder  ("vpvh","s1  ",1342,4)
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
      call reorder  ("vphv","s3  ",1423,4)
!   diag 7
      call fold4    ("tvm ","s1  ",4,4)
      call reorder  ("s1  ","s2  ",1423,4)
      call mult     ("s2  ","va2 ","s1  ",222)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",1.d0,"s3  ",-1.d0,"s3  ")
!   diag 8
      call reorder  ("ppvh","s1  ",3142,4)
      call mult     ("t1v ","s1  ","s2  ",131)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
!
      call mult6    ("tv7 ","s3  ","s1  ",1)
      call add      ("s1  ",1.d0,"boxZ",1.d0,"boxZ")
!
!   diag 9
      call clear6   ("s1  ",332333,110101)
      call fold4    ("tu  ","s2  ",3,4)
      call t1t2     ("t1c1","s2  ","s1  ","s3  ",123456,3)
      call fold4    ("tv  ","s1  ",3,4)
      call fold4    ("s1  ","s2  ",4,4)
      call t1t2     ("t1v ","s2  ","s3  ","s1  ",123456,2)
      call reord6   ("s1  ","s2  ",124635)
      call reorder  ("ppvh","s1  ",2341,4)
      call mult6    ("s2  ","s1  ","s3  ",3)
      call reord6   ("s3  ","s1  ",126345)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124356)
      call add      ("s1  ",-1.d0,"boxZ",1.d0,"boxZ")
!   diag 10
      call clear6   ("s1  ",322333,100111)
      call fold4    ("tv  ","s2  ",3,4)
      call fold4    ("s2  ","s3  ",4,4)
      call t1t2     ("t1c1","s3  ","s1  ","s2  ",123456,3)
      call t1t2     ("t1c1","s3  ","s2  ","s1  ",123456,2)
      call reord6   ("s1  ","s2  ",124356)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",125634)
      call reorder  ("vhhp","s1  ",1432,4)
      call mult6    ("s2  ","s1  ","s3  ",3)
      call reord7   ("s3  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",123546)
      call add      ("s2  ",1.d0,"boxZ",1.d0,"boxZ")
!
      call reord6   ("boxZ","s1  ",124356)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124356)
      call fold4    ("tu  ","s2  ",3,4)
      call mult8    ("s1  ","s2  ","s3  ",1)
      call reordx   ("s3  ","txf ",12673485)
!  GROUP II
!   diag 1 
      call cmult    ("vvhh",0.5d0,"s3  ")
!   diag 2
      call mult     ("tum ","hhpp","s1  ",222)
      call add      ("s1  ",0.5d0,"s3  ",1.d0,"s3  ")
!   diag 3
      call mult     ("t1v ","vhhp","s1  ",131)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!
      call mult6    ("s3  ","tv7 ","s1  ",1)
      call reord6   ("s1  ","s2  ",124653)
      call mult8    ("s2  ","tv7 ","s1  ",1)
      call reordx   ("s1  ","s2  ",12367548)
      call add      ("s2  ",1.d0,"txf ",1.d0,"txf ")
      call prtime   ("diags1  ",0,0)
      return
      end
!  *******************************************************************
      subroutine diags2_csx
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  diag 1
      call fold6    ("t33 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reorder  ("pvvv","s1  ",2341,4)
      call mult8    ("s2  ","s1  ","s3  ",1)
      call reordx   ("s3  ","s1  ",12364578)
      call add      ("s1  ",0.5d0,"txf ",1.d0,"txf ")
!  diag 2
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call fold6    ("s2  ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",124563)
      call mult8    ("s2  ","vvvh","s1  ",1)
      call reordx   ("s1  ","s2  ",12673485)
      call add      ("s2  ",-0.5d0,"txf ",1.d0,"txf ")
      call prtime   ("diags2  ",0,0)
      return
      end 
!  *******************************************************************
      subroutine diags3_csx
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  GROUP I
!   diag 1
      call rdfile   (23,0,"boxA",-1)     
! v1 Feb 95
      call fold4    ("tv  ","s1  ",3,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","s1  ",1342,4)
      call mult     ("s1  ","ph  ","s2  ",311)
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
      call add      ("s2  ",0.5d0,"txf ",1.d0,"txf ")
      call prtime   ("I   ",0,0)
!  GROUP II
!   diag 1
      call rdfile   (23,1,"boxB",-1)     
! v1 Feb 95
      call fold4    ("tu  ","s1  ",3,4)
      call mult     ("s1  ","hp  ","s2  ",311)
      call add      ("boxB",1.d0,"s2  ",-.5d0,"boxB")
!   diag 2
      call t1t2     ("t1c ","tu  ","t32 ","s1  ",123456,3)
      call fold6    ("s1  ","s2  ",4,1)
      call reord6   ("s2  ","s1  ",124356)
      call add      ("hhpp",-2.d0,"va6 ",1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","s3  ",2)
      call add      ("s3  ",1.d0,"boxB",1.d0,"boxB")
!   diag 3
      call fold6    ("t32 ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",126354)
      call mult6    ("s2  ","hhpp","s1  ",2)
      call add      ("s1  ",1.d0,"boxB",1.d0,"boxB")
!
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call fold6    ("s2  ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",124563)
      call mult8    ("s2  ","boxB","s1  ",1)
      call reordx   ("s1  ","s2  ",12673485)
      call add      ("s2  ",0.5d0,"txf ",1.d0,"txf ")
      call prtime   ("II  ",0,0)
!  GROUP III
!   diag 1
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord6   ("s2  ","s1  ",124536)
      call reorder  ("ppvh","s2  ",3142,4)
      call reorder  ("ppvh","s3  ",3241,4)
      call add      ("s2  ",2.d0,"s3  ",-1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","boxC",3)
!   diag 2
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124635)
      call reorder  ("ppvh","s2  ",3142,4)
      call mult6    ("s1  ","s2  ","s3  ",3)
      call add      ("s3  ",-1.d0,"boxC",0.5d0,"boxC")
!   diag 3
      call reorder  ("ppvh","s2  ",3241,4)
      call mult6    ("s1  ","s2  ","s3  ",3)
      call reord6   ("s3  ","s1  ",123546)
      call add      ("s1  ",-1.d0,"boxC",1.d0,"boxC")
!   diags 4,5
      call mult     ("t1v ","hhpp","s1  ",131)
      call add      ("s1  ",1.d0,"vhhp",1.d0,"s1  ")
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
      call fold4    ("tu  ","s1  ",3,4)
      call mult8    ("boxC","s1  ","s2  ",1)
      call reordx   ("s2  ","s1  ",12673485)
      call add      ("s1  ",1.d0,"txf ",1.d0,"txf ")
      call prtime   ("III ",0,0)
!  GROUP IV
!   diag 1
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord6   ("s2  ","s1  ",124536)
      call reorder  ("vhhp","s2  ",1324,4)
      call add      ("s2  ",-2.d0,"vhhp",1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","boxD",3)
!   diag 2
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124635)
      call reorder  ("vhhp","s2  ",1324,4)
      call mult6    ("s1  ","s2  ","s3  ",3)
      call add      ("s3  ",1.d0,"boxD",0.5d0,"boxD")
!   diag 3
      call mult6    ("s1  ","vhhp","s2  ",3)
      call reord6   ("s2  ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",124536)
      call add      ("s1  ",1.d0,"boxD",1.d0,"boxD")
!   diags 4,5
      call reorder  ("ppvh","s3  ",3412,4)
      call fold4    ("t1c ","s1  ",2,2)
      call reorder  ("s1  ","s2  ",1243,4)
      call mult     ("s2  ","va4 ","s1  ",131)
      call reorder  ("s1  ","s2  ",1432,4)
      call add      ("s2  ",1.d0,"s3  ",-1.d0,"s3  ")
      call fold6    ("t33 ","s1  ",4,1)
      call mult6    ("s1  ","s3  ","s2  ",3)
      call reord6   ("s2  ","s1  ",124536)
      call add      ("s1  ",1.d0,"boxD",1.d0,"boxD")
!
      call reord6   ("boxD","s1  ",125346)
      call mult8    ("s1  ","tv7 ","s2  ",1)
      call reordx   ("s2  ","s1  ",12364578)
      call add      ("s1  ",1.d0,"txf ",1.d0,"txf ")
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
      call add      ("s3  ",0.25d0,"txf ",1.d0,"txf ")
!   diag 2
      call add      ("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","s3  ",3)
      call t1t2     ("t1v ","tv  ","t32 ","s1  ",123456,2)
      call t1t2     ("t1c ","tu  ","s1  ","s2  ",123456,3)
      call fold6    ("s2  ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124635)
      call mult8    ("s3  ","s1  ","s2  ",5)
      call reordx   ("s2  ","s1  ",12563478)
      call add      ("s1  ",0.5d0,"txf ",1.d0,"txf ")
!   diag 3
      call cmult    ("t1v ",2.d0,"s1  ")
      call t1t2     ("s1  ","tv  ","t32 ","s2  ",123456,2)
      call cmult    ("t1c ",2.d0,"s1  ")
      call t1t2     ("s1  ","tu  ","s2  ","s3  ",123456,3)
      call fold6    ("s3  ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s3  ",124635)
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124635)
      call mult6    ("s1  ","va1 ","s2  ",3)
      call mult8    ("s2  ","s3  ","s1  ",5)
      call reordx   ("s1  ","s2  ",12563478)
      call add      ("s2  ",0.5d0,"txf ",1.d0,"txf ")
!   diag 4
      call cmult    ("t1v ",2.d0,"s1  ")
      call t1t2     ("s1  ","tv  ","t32 ","s2  ",123456,1)
      call cmult    ("t1c ",2.d0,"s1  ")
      call t1t2     ("s1  ","tu  ","s2  ","s3  ",123456,3)
      call fold6    ("s3  ","s1  ",5,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s3  ",125634)
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124635)
      call mult6    ("s1  ","va2 ","s2  ",3)
      call mult8    ("s2  ","s3  ","s1  ",5)
      call reordx   ("s1  ","s2  ",12563847)
      call add      ("s2  ",0.5d0,"txf ",1.d0,"txf ")
!   diag 5
      call fold6    ("t33 ","s1  ",4,1)
      call mult6    ("s1  ","hhpp","s2  ",3)
      call fold6    ("t31 ","s1  ",4,1)
      call fold6    ("s1  ","s3  ",5,1)
      call fold6    ("s3  ","s1  ",6,1)
      call reord6   ("s1  ","s3  ",124356)
      call reord7   ("s3  ","s1  ",132456)
      call reord6   ("s1  ","s3  ",125643)
      call mult8    ("s2  ","s3  ","s1  ",5)
      call reordx   ("s1  ","s2  ",12354876)
      call add      ("s2  ",0.5d0,"txf ",1.d0,"txf ")
      call prtime   ("V   ",0,0)
      return
      end 
!  *******************************************************************
