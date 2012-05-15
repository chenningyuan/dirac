!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine chy
!---------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use cint_mod
      use dat_mod
      use fock_mod
      use idv_mod
      use prm_mod
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
!     common /dat/  norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)
#include "fock.h"
!     common /csdtime/ csdstim,csdltim
#include "csdtime.h"
!     common /idv/ iv(200),nvalo
#include "idv.h"
!     common /prm/ iperm (0:6,720)
#include "prm.h"
#endif

      dimension iord(10)
      character*256 filenm,filenm6,filenm51,filenm53,filenm55,          &
 filenm60,filenm62
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
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')
      write (6,*) '   $$$$$$$$$   chy   $$$$$$$$$  '
      csdstim=cputime(real(0.))
      csdltim=csdstim
      call rdin
      call rdsym
      call prtime  ("rdin    ",0,0)
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      call reads   (51,0,"s1  ")
      call reorder ("s1  ","tc  ",3412,4)
      call reads   (51,1,"s1  ")
      call reorder ("s1  ","t1c ",1243,4)
      call rdfile  (53,0,"thv ",-1)
      call rdfile  (53,1,"t1hv",-1)
      call rdfile  (53,2,"thvf",-1)
      call rdfile  (55,0,"thu ",-1)
      call rdfile  (55,1,"thuf",-1)
      call rdfile  (60,0,"thwf",-1)
      call rdfile  (62,0,"thxf",-1)
      if (iopt(17).ge.3) then
        call rdfile  (53,3,"t31 ",-1)
        call rdfile  (55,2,"t32 ",-1)
        call rdfile  (60,1,"t33 ",-1)
      end if
      nvalo=0
      do 20 j=1,120
        if (icv(j).eq.1) then
          nvalo=nvalo+1
          iv(nvalo)=j      
        end if
20    continue
      call ordv_chy
      call ordta_chy
      call t1sqr    ("t1c ","t1c ","tc  ","tcm ")
      call t1sqr    ("t1hv","t1c ","thv ","tvm ")
      call t1sqr    ("t1hv","t1hv","thu ","tum ")
      call diags1_chy
      if (iopt(17).eq.3) then
        call diags31_chy
        call diags32_chy
        call diags33_chy
      end if
      call fmax    ("thyf",tmax,idum)
      write (6,99) tmax
99    format (1x,'tmax,unsymmetrised=',e10.3)
      if (abs(tmax).gt.1.d-13) then
        call perms (5)
        call move    ("thyf","s1  ")
        do 30 i=2,120
          do 40 j=1,5
            iord(j)=iperm(j,i)
40        continue
          call reord10 ("s1  ","s2  ",iord)
          call reord11 ("s2  ","s3  ",iord)
          call add     ("s3  ",1.d0,"thyf",1.d0,"thyf")
          if (i/10.d0-int(i/10.d0).lt.0.05d0) then
            call prtime ("perm",i,0)
          end if
30      continue
        call fmax    ("thyf",tmax,idum)
        write (6,98) tmax
98      format (1x,'tmax,symmetrised=',e10.3)
      end if
      egr0=0.d0
      call eneru2  (egr0,1,2,max(1,iopt(1)))
      call eneru2  (egr0,2,2,max(1,iopt(1)))
      call eneru2  (egr0,3,2,max(1,iopt(1)))                         
      call eneru2  (egr0,4,2,max(1,iopt(1)))                         
      call eneru2  (egr0,5,2,max(1,iopt(1)))                         
      call wrtfile (62,1,"thyf")
      write (6,*) nstrt(1,nflcur),'words used in igen,',                &
 nstrt(2,nflcur),' words used in gen'
      call prtime ("end ",0,0)
!      stop
      end
!  *******************************************************************
      subroutine ordv_chy
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
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)
#include "fock.h"
#endif

!
      call readv    (nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call prtime   ("readv   ",nint,0)
      call seltyp   (3322,0,"pphh",4)
      call seltyp   (2233,1000,"ghpp",4)
      call seltyp   (2223,1010,"ghgp",4)
      call seltyp   (2232,1001,"ghpg",4)
      call seltyp   (2222,0011,"hhgg",4)
      call seltyp   (2223,0010,"hhgp",4)
      call seltyp   (2233,1100,"ggpp",4)
      call selt1    (23,0,"hp  ",4)
      call prtime   ("ordv",0,0)
      return
      end
!  ********************************************************************
      subroutine ordta_chy
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call fold4    ("thu ","tua ",3,4)
      call fold4    ("t1c ","s1  ",2,2)
      call reorder  ("s1  ","tw1 ",1243,4)
      call fold4    ("thv ","s1  ",3,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","tw3 ",1342,4)
      call prtime   ("ordta   ",0,0)
      return
      end
!  **********************************************************************
      subroutine diags1_chy
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      dimension iord(10)
      data iord / 1,2,6,7,5,4,10,3,8,9 /
!   diags 1,3
      call reorder ("pphh","s1  ",1342,4)
      call mult    ("s1  ","tw1 ","s2  ",311)
      call reorder ("s2  ","s1  ",2143,4)
      call reorder ("hhgp","s2  ",2431,4)
      call add     ("s1  ",1.d0,"s2  ",1.d0,"s1  ")
      call mult6   ("tua ","s1  ","s2  ",1)
      call reord6  ("s2  ","s3  ",123645)
!   diags 2,4
      call reorder ("pphh","s1  ",3124,4)
      call mult    ("t1hv","s1  ","s2  ",131)
      call reorder ("ghpp","s1  ",1243,4)
      call add     ("s2  ",1.d0,"s1  ",1.d0,"s1  ")
      call mult6   ("tw3 ","s1  ","s2  ",1)
      call reord6  ("s2  ","s1  ",124356)
      call reord7  ("s1  ","s2  ",132456)
      call add     ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!
      call diage5  ("s3  ","tua ","tw3 ","thyf",iord)
      call prtime  ("diags1  ",0,0)
      return
      end
!  *******************************************************************
      subroutine diags31_chy
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      dimension iord(10)
      data iord / 1,2,3,6,7,5,4,10,9,8 /
!   diag 1
      call reorder   ("ghgp","s1  ",1324,4)
      call mult6     ("s1  ","tw3 ","s2  ",1)   
      call reord6    ("s2  ","s1  ",124563)
      call reord7    ("s1  ","s2  ",132456)
      call cmult     ("s2  ",-1.d0,"boxA")
!   diag 3
      call reorder   ("ghpg","s1  ",1243,4)
      call mult6     ("tw3 ","s1  ","s2  ",1)
      call reord6    ("s2  ","s1  ",124365)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",-1.d0,"boxA",1.d0,"boxA")
!   diag 5
      call reorder   ("hhgg","s1  ",3421,4)
      call mult6     ("tua ","s1  ","s2  ",1)
      call add       ("s2  ",1.d0,"boxA",1.d0,"boxA")
!   diag 7
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call reord6    ("s2  ","s1  ",124536)
      call reorder   ("hhgp","s2  ",3142,4)
      call reorder   ("hhgp","s3  ",3241,4)
      call add       ("s2  ",2.d0,"s3  ",-1.d0,"s2  ")
      call mult6     ("s1  ","s2  ","s3  ",3)
      call add       ("s3  ",-0.5d0,"boxA",1.d0,"boxA")
!   diag 8
      call t1t2      ("t1hv","thv ","t32 ","s1  ",123456,2)
      call t1t2      ("t1c ","thu ","s1  ","s2  ",123456,3)
      call fold6     ("s2  ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",6,1)
      call reord6    ("s2  ","s1  ",124635)
      call reorder   ("hhgp","s2  ",3142,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call add       ("s3  ",1.d0,"boxA",1.d0,"boxA")
!   diag 9
      call reorder   ("hhgp","s2  ",3241,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125346)
      call add       ("s1  ",1.d0,"boxA",1.d0,"boxA")
!   diag 14
      call t1t2      ("t1c ","thv ","t31 ","s1  ",123456,2)
      call t1t2      ("t1c ","thv ","s1  ","s2  ",123456,3)
      call fold6     ("s2  ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call fold6     ("s2  ","s1  ",6,1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125634)
      call mult6     ("s2  ","ghpp","s1  ",3)
      call reord6    ("s1  ","s2  ",125346)
      call reord7    ("s2  ","s1  ",132456)
      call add       ("s1  ",-1.d0,"boxA",1.d0,"boxA")

!   diag 15
      call reorder   ("pphh","s1  ",3124,4)
      call reorder   ("pphh","s2  ",3214,4)
      call add       ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call mult      ("s1  ","t1c ","s3  ",202)
! g1 Apr 95
      call add       ("s3  ",1.d0,"hp  ",-1.d0,"s3  ")

      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call fold6     ("s2  ","s1  ",6,1)
      call reord6    ("s1  ","s2  ",124563)  
      call mult6     ("s2  ","s3  ","s1  ",5)
      call add       ("s1  ",0.5d0,"boxA",1.d0,"boxA")
!   diag 18
      call reorder   ("pphh","s1  ",3124,4)
      call reorder   ("pphh","s2  ",3214,4)
      call add       ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call fold4     ("thv ","s2  ",3,4)
      call reorder   ("s2  ","s3  ",1324,4)
      call mult      ("s3  ","s1  ","s2  ",222)
      call mult6     ("s2  ","tw3 ","s1  ",1)
      call reord6    ("s1  ","s2  ",124563)
      call reord7    ("s2  ","s1  ",132456)
      call add       ("s1  ",-1.d0,"boxA",1.d0,"boxA")
!   diag 20
      call reorder   ("pphh","s1  ",3241,4)
      call reorder   ("pphh","s2  ",3142,4)
      call add       ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call mult      ("tw1 ","s1  ","s2  ",131)
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s3  ",5,1)
      call reord6    ("s3  ","s1  ",124536)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call add       ("s3  ",0.5d0,"boxA",1.d0,"boxA")
!   diag 21
      call fold4     ("tvm ","s1  ",4,4)
      call reorder   ("s1  ","s2  ",1423,4)
      call reorder   ("pphh","s1  ",3124,4)
      call mult      ("s2  ","s1  ","s3  ",222)
      call mult6     ("s3  ","tw3 ","s1  ",1)
      call reord6    ("s1  ","s2  ",124563)
      call reord7    ("s2  ","s1  ",132456)
      call add       ("s1  ",1.d0,"boxA",1.d0,"boxA")
!   diag 22
      call reorder   ("pphh","s1  ",3241,4)
      call mult      ("tw1 ","s1  ","s2  ",131)
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s3  ",6,1)
      call reord6    ("s3  ","s1  ",124635)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call add       ("s3  ",1.d0,"boxA",1.d0,"boxA")
!   diag 26
      call reorder   ("pphh","s2  ",3142,4)
      call mult      ("tw1 ","s2  ","s3  ",131)
      call mult6     ("s1  ","s3  ","s2  ",3)
      call reord6    ("s2  ","s1  ",123546)
      call add       ("s1  ",1.d0,"boxA",1.d0,"boxA")
!   diag 27
      call fold4     ("tvm ","s1  ",4,4)
      call reorder   ("s1  ","s2  ",1423,4)
      call reorder   ("pphh","s1  ",3214,4)
      call mult      ("s2  ","s1  ","s3  ",222)
      call mult6     ("tw3 ","s3  ","s1  ",1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call add       ("s1  ",1.d0,"boxA",1.d0,"boxA")
!   diag 34
      call reorder   ("pphh","s1  ",3214,4)
      call mult      ("t1hv","s1  ","s2  ",131)
      call fold6     ("t31 ","s1  ",4,1)
      call fold6     ("s1  ","s3  ",5,1)
      call fold6     ("s3  ","s1  ",6,1)
      call reord6    ("s1  ","s3  ",124356)
      call reord7    ("s3  ","s1  ",132456)
      call reord6    ("s1  ","s3  ",125634)
      call mult6     ("s3  ","s2  ","s1  ",3)
      call reord6    ("s1  ","s2  ",125346)
      call reord7    ("s2  ","s1  ",132456)
      call add       ("s1  ",1.d0,"boxA",1.d0,"boxA")
!   diag 36
      call fold4     ("tcm ","s1  ",3,4)
      call fold4     ("s1  ","s2  ",4,4)
      call reorder   ("s2  ","s1  ",3421,4)
      call reorder   ("pphh","s2  ",3412,4)
      call mult      ("s1  ","s2  ","s3  ",222)
      call mult6     ("tua ","s3  ","s1  ",1)
      call add       ("s1  ",1.d0,"boxA",1.d0,"boxA")
!
      call fold6     ("t33 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call dge5n2    ("s2  ","boxA","s1  ",iord)
      call add       ("s1  ",0.5d0,"thyf",1.d0,"thyf")
      call prtime    ("diags31 ",0,0)
      return   
      end
!  ********************************************************************
      subroutine diags32_chy
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension iord(10)
      data iord / 1,2,3,6,7,5,4,10,9,8 /
!   diag 2
      call reorder   ("ghgp","s1  ",1342,4)
      call mult6     ("s1  ","tua ","s2  ",1)
      call reord6    ("s2  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",125364)
      call cmult     ("s1  ",-1.d0,"boxB")
!   diag 4
      call reorder   ("ghpg","s1  ",1432,4)
      call mult6     ("tua ","s1  ","s2  ",1)
      call reord6    ("s2  ","s1  ",124356)
      call add       ("s1  ",-1.d0,"boxB",1.d0,"boxB")
!   diag 6
      call reorder   ("ggpp","s1  ",1243,4)
      call mult6     ("tw3 ","s1  ","s2  ",1)
      call reord6    ("s2  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",125346)
      call add       ("s1  ",1.d0,"boxB",1.d0,"boxB")
!   diag 10
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call reord6    ("s2  ","s1  ",124536)
      call reorder   ("ghpp","s2  ",1342,4)
      call reorder   ("ghpp","s3  ",1432,4)
      call add       ("s2  ",-2.d0,"s3  ",1.d0,"s2  ")
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125346)
      call add       ("s1  ",-0.5d0,"boxB",1.d0,"boxB")
!   diag 11
      call t1t2      ("t1hv","thv ","t32 ","s1  ",123456,2)
      call t1t2      ("t1c ","thu ","s1  ","s2  ",123456,3)
      call fold6     ("s2  ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",6,1)
      call reord6    ("s2  ","s1  ",124635)
      call reorder   ("ghpp","s2  ",1342,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125346)
      call add       ("s1  ",-1.d0,"boxB",1.d0,"boxB")
!   diag 12
      call t1t2      ("t1hv","thv ","t32 ","s1  ",123456,1)
      call t1t2      ("t1c ","thu ","s1  ","s2  ",123456,3)
      call fold6     ("s2  ","s1  ",5,1)
      call fold6     ("s1  ","s2  ",6,1)
      call reord6    ("s2  ","s1  ",125634)
      call reorder   ("ghpp","s2  ",1432,4)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",124653)
      call add       ("s1  ",-1.d0,"boxB",1.d0,"boxB")
!   diag 13
      call t1t2      ("t1hv","thu ","t33 ","s1  ",123456,2)
      call t1t2      ("t1hv","thu ","s1  ","s2  ",123456,3)
      call fold6     ("s2  ","s1  ",4,1)
      call reorder   ("hhgp","s2  ",3412,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call add       ("s3  ",1.d0,"boxB",1.d0,"boxB")
!   diag 16
      call reorder   ("pphh","s1  ",1324,4)
      call reorder   ("pphh","s2  ",1423,4)
      call add       ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call fold4     ("thv ","s2  ",3,4)
      call reorder   ("s2  ","s3  ",1324,4)
      call mult      ("s3  ","s1  ","s2  ",222)
      call mult6     ("s2  ","tua ","s1  ",1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125364)
      call add       ("s2  ",-1.d0,"boxB",1.d0,"boxB")
!   diag 17
      call reorder   ("pphh","s1  ",1243,4)
      call add       ("s1  ",2.d0,"pphh",-1.d0,"s1  ")
      call mult      ("t1hv","s1  ","s2  ",131)
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s3  ",5,1)
      call reord6    ("s3  ","s1  ",124536)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125346)
      call add       ("s1  ",-0.5d0,"boxB",1.d0,"boxB")
!   diag 23
      call fold4     ("tvm ","s1  ",4,4)
      call reorder   ("s1  ","s2  ",1423,4)
      call reorder   ("pphh","s1  ",1324,4)
      call mult      ("s2  ","s1  ","s3  ",222)
      call mult6     ("s3  ","tua ","s1  ",1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125364)
      call add       ("s2  ",1.d0,"boxB",1.d0,"boxB")
!   diag 25
      call reorder   ("pphh","s1  ",1243,4)
      call mult      ("t1hv","s1  ","s2  ",131)
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s3  ",6,1)
      call reord6    ("s3  ","s1  ",124635)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125346)
      call add       ("s1  ",1.d0,"boxB",1.d0,"boxB")
!   diag 29
      call fold4     ("tvm ","s1  ",4,4)
      call reorder   ("s1  ","s2  ",1423,4)
      call reorder   ("pphh","s1  ",1423,4)
      call mult      ("s2  ","s1  ","s3  ",222)
      call mult6     ("tua ","s3  ","s1  ",1)
      call reord6    ("s1  ","s2  ",124356)
      call add       ("s2  ",1.d0,"boxB",1.d0,"boxB")
!   diag 30
      call mult      ("t1hv","pphh","s1  ",131)
      call fold6     ("t32 ","s2  ",5,1)
      call fold6     ("s2  ","s3  ",6,1)
      call reord6    ("s3  ","s2  ",125634)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord6    ("s3  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",124653)
      call add       ("s1  ",1.d0,"boxB",1.d0,"boxB")
!   diag 31
      call reorder   ("pphh","s1  ",1342,4)
      call mult      ("tw1 ","s1  ","s2  ",131)
      call fold6     ("t33 ","s1  ",5,1)
      call reord6    ("s1  ","s3  ",123564)
      call mult6     ("s3  ","s2  ","s1  ",3)
      call reord6    ("s1  ","s2  ",123546)
      call add       ("s2  ",1.d0,"boxB",1.d0,"boxB")
!   diag 33
      call reorder   ("pphh","s1  ",1243,4)
      call mult      ("tum ","s1  ","s2  ",222)
      call mult6     ("tw3 ","s2  ","s1  ",1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125346)
      call add       ("s2  ",1.d0,"boxB",1.d0,"boxB")
!
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call fold6     ("s2  ","s1  ",6,1)
      call reord6    ("s1  ","s2  ",124563)
      call dge5n2    ("boxB","s2  ","s1  ",iord)
      call add       ("s1  ",0.5d0,"thyf",1.d0,"thyf")
      call prtime    ("diags32 ",0,0)
      return
      end
!  *******************************************************************
      subroutine diags33_chy
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension iord1(10),iord2(10),iord3(10),iord4(10)
      data iord1 / 1,2,6,7,3,8,9,5,4,10 /
      data iord2 / 1,3,8,9,2,7,6,5,4,10 /
      data iord3 / 1,3,4,8,2,7,6,5,10,9 /
      data iord4 / 1,6,7,9,3,2,8,5,4,10 /
!   diag 19
      call reorder   ("pphh","s1  ",1324,4)
      call reorder   ("pphh","s2  ",1423,4)
      call add       ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call fold6     ("t32 ","s2  ",4,1)
      call fold6     ("s2  ","s3  ",5,1)
      call reord6    ("s3  ","s2  ",124536)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call diage5    ("s3  ","tw3 ","tua ","s1  ",iord1)
      call add       ("s1  ",0.5d0,"thyf",1.d0,"thyf")
!   diag 24
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",6,1)
      call reord6    ("s2  ","s1  ",124635)
      call reorder   ("pphh","s2  ",1324,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call diage5    ("s3  ","tw3 ","tua ","s1  ",iord1)
      call add       ("s1  ",1.d0,"thyf",1.d0,"thyf")
!   diag 28
      call fold6     ("t32 ","s1  ",5,1)
      call fold6     ("s1  ","s2  ",6,1)
      call reord6    ("s2  ","s1  ",125634)
      call reorder   ("pphh","s2  ",1423,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call diage5    ("s3  ","tw3 ","tua ","s1  ",iord2)
      call add       ("s1  ",1.d0,"thyf",1.d0,"thyf")
!   diag 32
      call fold6     ("t33 ","s1  ",5,1)
      call reord6    ("s1  ","s2  ",123546)
      call mult6     ("s2  ","pphh","s1  ",3)
      call diage5    ("s1  ","tw3 ","tw3 ","s2  ",iord3)
      call add       ("s2  ",0.5d0,"thyf",1.d0,"thyf")
!   diag 35
      call fold6     ("t31 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call fold6     ("s2  ","s1  ",6,1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125634)
      call reorder   ("pphh","s1  ",3412,4)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call diage5    ("s3  ","tua ","tua ","s1  ",iord4)
      call add       ("s1  ",0.5d0,"thyf",1.d0,"thyf")
      call prtime    ("diags33 ",0,0)
      return
      end
!  *******************************************************************
