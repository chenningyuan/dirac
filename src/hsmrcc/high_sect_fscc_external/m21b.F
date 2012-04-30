!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine m21b
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

      character*256 filenm,filenm6,
     1filenm51,filenm52,filenm53,filenm54,filenm56,filenm57,filenm58

      real*4 cputime


c     iopt(1) - print.
c     iopt(5) - print only elements .gt. 10**(-iopt(5))
c     iopt(8) - fock.
c     iopt(10) .ne.0 - do not write tape9
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
      filenm54(1:1)=' '
      call getenv('tape54',filenm54)
      if (filenm54(1:1).eq.' ') filenm54='tape54'
      open(54,file=filenm54,form='UNFORMATTED')
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
      write (*,*) '   $$$$$$$$$   m21b   $$$$$$$$$  '
      csdstim = cputime(real(0.))
      csdltim = csdstim
      call rdin
      call rdsym
      call prtime ("rdin    ",0,0)
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      call rdfile (51,0,"tc  ",-1)
      call rdfile (51,1,"t1c ",-1)
      call rdfile (52,0,"tv  ",-1)
      call rdfile (52,1,"t1v ",-1)
      call rdfile (52,2,"tvf ",-1)
      call rdfile (56,0,"tp  ",-1)
      call rdfile (56,1,"tpf ",-1)
      call rdfile (57,0,"tq  ",-1)
      call rdfile (57,1,"tqf ",-1)
      call rdfile (58,1,"t21 ",-1)
      call rdfile (53,0,"tg  ",-1)
      call rdfile (53,1,"t1g ",-1)    
      call rdfile (53,2,"tgf ",-1)
      call rdfile (54,0,"tu  ",-1)
      call rdfile (54,1,"tuf ",-1)
      call t1sqr ("t1c ","t1c ","tc  ","tm  ")
      call reorder ("t1c ","s1  ",1243,4)
      call t1sqr ("t1v ","t1v ","tu  ","tum ")
      call reorder ("t1g ","s1  ",1243,4)
      call t1sqr ("t1v ","s1  ","tp  ","tpm ")
      call ordv_m21b
      call heffq
      call heffp
      egr0=0.d0
      if (iopt(11).eq.1) then
        call clear6 ("ttpf",332332,111111)
        call clear6 ("ttqf",332233,111111)
      end if
      call eneru2 (egr0,1,3,max(1,iopt(1)))
      call eneru2 (egr0,21,3,max(1,iopt(1)))
      write (6,*) nstrt(1,nflcur),'words used in igen,',
     1nstrt(2,nflcur),' words used in gen'
      call prtime ("end ",0,0)
!     stop
      end

C  ********************************************************************
      subroutine ordv_m21b
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
      call selt1    (32,0,"ph  ",4)   
      call seltyp   (2233,0,"hhpp",4)
      call seltyp   (3233,0111,"pgvv",4)
      call seltyp   (3332,0010,"ppvh",4)
      call seltyp   (3323,0001,"pphv",4)
      call seltyp   (3333,0011,"ppvv",4)
      call seltyp   (3222,0100,"pghh",4)
      call seltyp   (3232,0110,"pgvh",4)
      call seltyp   (3223,0101,"pghv",4)
      call seltyp   (3322,1000,"vphh",4)
      call seltyp   (3322,0001,"pphg",4)
      call seltyp   (3322,1010,"vpgh",4)
      call seltyp   (3322,1100,"vvhh",4)
      call reorder  ("hhpp","vm1 ",3124,4)
      call reorder  ("ppvh","vm2 ",1342,4)
      call reorder  ("hhpp","vm4 ",3412,4)
      call reorder  ("pphv","vm5 ",1432,4)
      call reorder  ("hhpp","vm8 ",3214,4)
      call reorder  ("hhpp","vm10",3214,4)
      call fold4    ("t1c ","t1cb",1,2)
      call fold4    ("tv  ","s1  ",2,4)
      call fold4    ("s1  ","tvc ",3,4)
      call fold4    ("tv  ","tvg ",2,4)
      call fold4    ("tg  ","tgc ",2,4)
      call fold4    ("tgc ","tgd ",3,4)
      call fold4    ("tu  ","tua ",3,4)
      call fold4    ("tq  ","tqa ",2,4)
      call fold4    ("tg  ","s1  ",2,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","tr1 ",1243,4)
      call fold4    ("tv  ","s1  ",3,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","tr2 ",1342,4)
      call fold4    ("tq  ","s1  ",4,4)
      call reorder  ("s1  ","tr3 ",1342,4)
      call fold4    ("tp  ","s1  ",3,4)
      call reorder  ("s1  ","tr4 ",1342,4)
      call fold4    ("tv  ","s1  ",2,4)
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
      call fold4    ("tv  ","s1  ",3,4)
      call reorder  ("s1  ","tr16",1342,4)
      call fold4    ("tv  ","s1  ",4,4)
      call reorder  ("s1  ","tr17",1432,4)
      call fold4    ("tg  ","s1  ",4,4)
      call reorder  ("s1  ","tr18",1423,4)
      call reorder  ("tp  ","tr19",1432,4)
      call fold4    ("tp  ","s1  ",2,4)
      call reorder  ("s1  ","tr20",1243,4)
      call reorder  ("tgd ","tr21",2134,4)
      call fold4    ("tg  ","s1  ",3,4)
      call reorder  ("s1  ","tr22",1324,4)
      call prtime   ("ordv",0,0)
      return
      end
C  ********************************************************************
      subroutine heffq
!---------------------------------------------------------------------
!
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)
#include "cfiles.h"
#endif

C
C     GROUP A
C
c  diags hb1,2
      call add      ("vm1 ",2.d0,"vm8 ",-1.d0,"s1  ")
      call mult     ("s1  ","t1c ","s2  ",202)
      call add      ("s2  ",1.d0,"ph  ",1.d0,"hb1 ")
      call reorder  ("hb1 ","s3  ",1243,4)
c  diag a.1
      call mult     ("tua ","s3  ","s1  ",311)
      call mult6    ("s1  ","tr1 ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tb1 ")
c  diag a.2
      call mult     ("tr2 ","hb1 ","s1  ",311)
      call mult6    ("s1  ","tqa ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tb2 ")
c  diag a.3
      call mult     ("tr3 ","hb1 ","s1  ",311)
      call mult6    ("tvc ","s1  ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tb3 ")
c  diag a.4
      call mult     ("tr4 ","hb1 ","s1  ",311)
      call mult6    ("s1  ","tr5 ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tb4 ")
c  diag a.5
      call mult     ("tr5 ","s3  ","s1  ",311)
      call mult6    ("s1  ","tr3 ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tb5 ")
      call prtime   ("group a ",1,0)
C
C     GROUP B
C
c     diag hbv 
c     diag hbv-a
      call move     ("pgvv","hbv ")
c  diag hbv-b
      call mult     ("ph  ","tr6 ","s1  ",131)
      !call add      ("s1  ",-1.d0,"hbv ",1.d0,"hbv ",1.d0)
      call add      ("s1  ",-1.d0,"hbv ",1.d0,"hbv ")
c  diag hbv-c
      call reorder  ("ppvh","s2  ",2341,4)
      call add      ("vm2 ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult     ("s1  ","tr7 ","s3  ",222)
c  diag hbv-d
      call mult     ("vm2 ","tr8 ","s1  ",222)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s1  ")
      call reorder  ("s1  ","s2  ",1324,4)
      call add      ("s2  ",1.d0,"hbv ",1.d0,"hbv ")
c  diag hbv-e
      call mult     ("vm5 ","tr8 ","s1  ",222)
      call reorder  ("s1  ","s2  ",1342,4)
      call add      ("s2  ",-1.d0,"hbv ",1.d0,"hbv ")
c  diag hbv-f
      call mult     ("pghh","tr9 ","s1  ",222)
      call add      ("s1  ",1.d0,"hbv ",1.d0,"hbv ")
c  diag hbv-j
      call add      ("vm1 ",2.d0,"vm8 ",-1.d0,"s1  ")
      call mult     ("s1  ","t1c ","s2  ",202)
      call mult     ("s2  ","tr6 ","s1  ",131)
      call add      ("s1  ",-1.d0,"hbv ",1.d0,"hbv ")
c  diag hbv-h,k,l
      call mult     ("vm1 ","tr10","s3  ",222)
      call add      ("vm1 ",2.d0,"vm10 ",-1.d0,"s1  ")
      call mult     ("s1  ","tr7 ","s2  ",222)
      call add      ("s2  ",1.d0,"s3  ",-1.d0,"s3  ")
      call reorder  ("s3  ","s1  ",1342,4)
      call reorder  ("pghv","s2  ",1243,4)
      call add      ("s2  ",1.d0,"s1  ",1.d0,"s2  ")
      call mult     ("s2  ","tr11","s1  ",311)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",-1.d0,"hbv ",1.d0,"hbv ")
c  diag hbv-i,m
      call mult     ("vm10","tr10","s1  ",222)
      call reorder  ("s1  ","s2  ",1342,4)
      call add      ("pgvh",1.d0,"s2  ",-1.d0,"s1  ")
      call mult     ("s1  ","tr11","s2  ",311)
      call add      ("s2  ",-1.d0,"hbv ",1.d0,"hbv ")
c  diag hbv-g,n
      call mult     ("vm4 ","tr9 ","s2  ",222)
      call reorder  ("s2  ","s1  ",1342,4)
      call reorder  ("ppvv","s2  ",1342,4)
      call add      ("s2  ",1.d0,"s1  ",1.d0,"s2  ")
      call mult     ("s2  ","t1cb","s1  ",311)
      call reorder  ("s1  ","s2  ",1423,4)
      call add      ("s2  ",1.d0,"hbv ",1.d0,"hbv ")
c  diag b.1 
      call reorder  ("hbv ","s1  ",2341,4)
      call reorder  ("t21 ","s2  ",2134,4)
      call mult6    ("s2  ","s1  ","tb6 ",1)
      call prtime   ("group b ",2,0)
C
C     GROUP C
C
c  diag hbu-a
      call move     ("ppvv","hbu ")
c  diag hbu-b
      call mult     ("vm4 ","tr9 ","s1  ",222)
      call add      ("s1  ",1.d0,"hbu ",1.d0,"hbu ")
c  diag hbu-c,d
      call reorder  ("pphv","s2  ",1243,4)
      call mult     ("s2  ","tr11","s1  ",311)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",-1.d0,"hbu ",1.d0,"hbu ")
      call reorder  ("s1  ","s2  ",2134,4)
      call add      ("s2  ",-1.d0,"hbu ",1.d0,"hbu ")
c  diag c.1,c.2
      call fold4    ("hbu ","s1  ",1,4)
      call reorder  ("s1  ","boxc",1342,4)
c
      call reorder  ("hbu ","s1  ",3421,4)
      call mult     ("t1v ","s1  ","s2  ",131)
      call add      ("s2  ",1.d0,"boxc",1.d0,"boxc")
c
      call mult6    ("boxc","tqa ","s1  ",1)
      call add      ("s1  ",1.d0,"tb2 ",1.d0,"tb2 ")
      call prtime   ("group c ",3,0)
C
C     GROUP D
C
c  diag hbp-a
      call move     ("pgvh","hbp ")
c  diag hbp-b,d
      call mult     ("vm8 ","tr8 ","s3  ",222)
      call reorder  ("pghh","s2  ",1423,4)
      call mult     ("s2  ","tr11","s1  ",311)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
      call reorder  ("s3  ","s1  ",1342,4)
      call add      ("s1  ",-1.d0,"hbp ",1.d0,"hbp ")
c  diag hbp-c
      call mult     ("vm2 ","t1cb","s1  ",311)
      call reorder  ("s1  ","s2  ",1423,4)
      call add      ("s2  ",1.d0,"hbp ",1.d0,"hbp ")
c  diag d.1,d.2   
      call reorder  ("hbp ","s1  ",2341,4)
      call mult     ("t1v ","s1  ","boxd",131)
c
      call fold4    ("hbp ","s1  ",1,4)
      call add      ("s1  ",1.d0,"boxd",1.d0,"boxd")
c
      call mult6    ("boxd","tr3 ","s2  ",1)
      call add      ("s2  ",-1.d0,"tb3 ",1.d0,"tb3 ")
      call prtime   ("group d ",4,0)
C
C     GROUP E
C
c  diag hbq-a
      call move     ("pghv","hbq ")
c  diag hbq-b,c
      call add      ("vm1 ",2.d0,"vm8 ",-1.d0,"s1  ")
      call mult     ("s1  ","tr7 ","s2  ",222)
      call mult     ("vm1 ","tr8 ","s3  ",222)
      call add      ("s2  ",1.d0,"s3  ",-1.d0,"s3  ")
      call reorder  ("s3  ","s1  ",1324,4)
      call add      ("s1  ",1.d0,"hbq ",1.d0,"hbq ")
c  diag hbq-d
      call mult     ("vm5 ","t1cb","s1  ",311)
      call reorder  ("s1  ","s2  ",1432,4)
      call add      ("s2  ",1.d0,"hbq ",1.d0,"hbq ")
c  diag hbq-e
      call mult     ("pghh","tr11","s1  ",311)
      call add      ("s1  ",-1.d0,"hbq ",1.d0,"hbq ")
c  diag e.1,e.4
      call reorder  ("hbq ","s1  ",4213,4)
      call mult     ("t1g ","s1  ","s3  ",131)
      call fold4    ("hbq ","s1  ",3,4)
      call reorder  ("s1  ","s2  ",3421,4)
      call add      ("s3  ",-1.d0,"s2  ",1.d0,"s1  ")
      call mult6    ("tua ","s1  ","s2  ",1)
      call add      ("s2  ",1.d0,"tb1 ",1.d0,"tb1 ")
c  diag e.2,e.3,e.5,e.6
      call reorder  ("hbq ","s1  ",2431,4)
      call mult     ("t1v ","s1  ","s3  ",131)
      call fold4    ("hbq ","s1  ",1,4)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
      call mult6    ("s3  ","tr3 ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb5 ",1.d0,"tb5 ")
      call mult6    ("tr4 ","s3  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb4 ",1.d0,"tb4 ")
      call prtime   ("group e ",5,0)
C
C     GROUP F
C
c  diag hbs-a
      call move     ("ppvh","hbs ")
c  diag hbs-b
      call reorder  ("hhpp","s2  ",3421,4)
      call mult     ("s2  ","tr11","s1  ",311)
      call reorder  ("s1  ","s2  ",1243,4) 
      call add      ("hbs ",1.d0,"s2  ",-1.d0,"hbs ")
      call reorder  ("hbs ","hbs1",3124,4)
      call reorder  ("hbs ","hbs2",3214,4)
c  diags f.1,7,12,20,21
      call fold4    ("hbs ","s1  ",2,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","s3  ",2431,4)
c
      call reorder  ("s1  ","s2  ",2134,4)
      call mult     ("s2  ","t1g ","s1  ",311)
      call reorder  ("s1  ","s2  ",1432,4)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
c
      call fold4    ("hbs ","s1  ",4,4)
      call reorder  ("s1  ","s2  ",4312,4)
      call mult     ("t1v ","s2  ","s1  ",131)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
c
      call add      ("hbs1",2.d0,"hbs2",-1.d0,"s1  ")
      call mult     ("tr14","s1  ","s2  ",222)
      call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
c    
      call mult     ("tr15","hbs1","s1  ",222)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
c
      call mult6    ("tvc ","s3  ","s1  ",1)
      call add      ("s1  ",1.d0,"tb3 ",1.d0,"tb3 ")
c  diags f.2,3,11,17,18
      call mult6    ("tr5 ","s3  ","s1  ",1)
      call add      ("s1  ",1.d0,"tb5 ",1.d0,"tb5 ")
c  diags f.4,6,10,19
      call fold4    ("hbs ","s1  ",1,4)
      call fold4    ("s1  ","s2  ",4,4)
      call  reorder ("s2  ","boxy",1342,4)
c
      call mult     ("s1  ","t1g ","s2  ",311)
      call reorder  ("s2  ","s3  ",1342,4)
      call add      ("s3  ",-1.d0,"boxy",1.d0,"boxy")
c
      call fold4    ("hbs ","s1  ",4,4)
      call reorder  ("s1  ","s2  ",3421,4)
      call mult     ("t1v ","s2  ","s1  ",131)
      call add      ("s1  ",1.d0,"boxy",1.d0,"boxy")
c
      call mult     ("tr15","hbs2","s1  ",222)
      call reorder  ("s1  ","s2  ",1324,4)
      call add      ("s2  ",-1.d0,"boxy",1.d0,"boxy")
c
      call mult6    ("boxy","tr5 ","s1  ",1)
      call add      ("s1  ",1.d0,"tb4 ",1.d0,"tb4 ")
c  diags f.5,9,8,25
      call fold4    ("hbs ","s1  ",1,4)
      call fold4    ("s1  ","boxx",2,4)
c
      call fold4    ("hbs ","s1  ",1,4)
      call reorder  ("s1  ","s2  ",1432,4)
      call mult     ("s2  ","t1v ","s1  ",311)
      call reorder  ("s1  ","s2  ",1432,4)
      call add      ("s2  ",1.d0,"boxx",1.d0,"boxx")
c
      call fold4    ("hbs ","s1  ",2,4)
      call reorder  ("s1  ","s2  ",2341,4)
      call mult     ("t1v ","s2  ","s1  ",131)
      call add      ("s1  ",1.d0,"boxx",1.d0,"boxx")
c
      call reorder  ("hbs ","s1  ",3412,4)
      call mult     ("tum ","s1  ","s2  ",222)
      call add      ("s2  ",1.d0,"boxx",1.d0,"boxx")
c
      call mult6    ("boxx","tr1 ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb1 ",1.d0,"tb1 ")
c  diags f.13,14,15
      call add      ("hbs1",2.d0,"hbs2",-1.d0,"s3  ")
      call mult     ("tr16","s3  ","boxw",222)
c
      call mult     ("tr17","hbs1","s1  ",222)
      call add      ("s1  ",-1.d0,"boxw",1.d0,"boxw")
c
      call mult     ("tr17","hbs2","s1  ",222)
      call reorder  ("s1  ","s2  ",1324,4)
      call add      ("s2  ",-1.d0,"boxw",1.d0,"boxw")
c
      call mult6    ("boxw","tqa ","s1  ",1)
      call add      ("s1  ",1.d0,"tb2 ",1.d0,"tb2 ")
c  diag f.16
      call mult     ("tr18","hbs2","s1  ",222)
      call reorder  ("s1  ","s2  ",1324,4)
      call mult6    ("tua ","s2  ","s1  ",1)        
      call add      ("s1  ",-1.d0,"tb1 ",1.d0,"tb1 ")
c  diags f.22,23
      call reorder  ("hbs ","s1  ",3421,4)
      call mult     ("tvg ","s1  ","s3  ",222)
      call mult6    ("tr4 ","s3  ","s2  ",1)
      call add      ("s2  ",-1.d0,"tb4 ",1.d0,"tb4 ")
c
      call mult6    ("s3  ","tr3 ","s2  ",1)
      call add      ("s2  ",-1.d0,"tb5 ",1.d0,"tb5 ")
c  diag f.24
      call reorder  ("hbs ","s1  ",3412,4)
      call mult     ("tvg ","s1  ","s2  ",222)
      call mult6    ("s2  ","tr3 ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb3 ",1.d0,"tb3 ")
      call prtime   ("group f ",6,0)
C
C     GROUP G
C
c     diag hbt-a
      call move     ("pghh","hbt ")
c     diag hbt-b
      call mult     ("vm8 ","t1cb","s1  ",311)
      call reorder  ("s1  ","s2  ",1432,4)
      call add      ("s2  ",1.d0,"hbt ",1.d0,"hbt ")
      call reorder  ("hbt ","hbt1",2314,4)
      call reorder  ("hbt ","hbt2",2413,4)
c   diags g.1,2,3,9,10
      call fold4    ("hbt ","s1  ",1,4)
      call fold4    ("s1  ","s3  ",3,4)
c
      call reorder  ("s1  ","s2  ",1243,4)
      call mult     ("s2  ","t1g ","s1  ",311)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
c
      call fold4    ("hbt ","s1  ",3,4)
      call reorder  ("s1  ","s2  ",2341,4)
      call mult     ("t1v ","s2  ","s1  ",131)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
c
      call add      ("hbt2",2.d0,"hbt1",-1.d0,"s1  ")
      call mult     ("tr14","s1  ","s2  ",222)
      call reorder  ("s2  ","s1  ",1324,4)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
c
      call mult     ("tr15","hbt2","s1  ",222)
      call reorder  ("s1  ","s2  ",1324,4)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
c
      call mult6    ("tr2 ","s3  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb2 ",1.d0,"tb2 ")
c   diag g.5,6
      call add      ("hbt2",2.d0,"hbt1",-1.d0,"s1  ")
      call mult     ("tr16","s1  ","boxg",222)
c
      call mult     ("tr17","hbt2","s1  ",222)
      call add      ("s1  ",-1.d0,"boxg",1.d0,"boxg")
c
      call reorder  ("boxg","s2  ",1324,4)
      call mult6    ("s2  ","tr3 ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb3 ",1.d0,"tb3 ")
c   diag g.7,8
      call mult     ("tr17","hbt1","s1  ",222)
      call reorder  ("s1  ","s3  ",1324,4)
      call cmult    ("s3  ",-1.d0,"s3  ")
c
      call mult6    ("s3  ","tr3 ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb5 ",1.d0,"tb5 ")
c
      call mult6    ("tr4 ","s3  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb4 ",1.d0,"tb4 ")
c   diag g.4
      call reorder  ("hbt ","s1  ",2134,4)
      call mult     ("tgc ","s1  ","s2  ",222)
      call mult6    ("tua ","s2  ","s1  ",1)
      call add      ("s1  ",1.d0,"tb1 ",1.d0,"tb1 ")
      call prtime   ("group g ",7,0)
C
C     GROUP H
C
      call mult     ("t1v ","vm8 ","s1  ",131)
      call add      ("s1  ",1.d0,"vphh",1.d0,"boxp")
      call reorder  ("boxp","bp1 ",1324,4)
      call reorder  ("boxp","bp2 ",1423,4)
      call reorder  ("hhpp","s2  ",3412,4)
      call mult     ("s2  ","t1g ","s1  ",311)
      call add      ("s1  ",-1.d0,"pphg",1.d0,"boxh")
c   diag h.1,2,16,19
      call mult     ("tgc ","boxp","s1  ",222)
      call reorder  ("s1  ","s3  ",3124,4)
c
      call reorder  ("boxh","s1  ",4123,4)
      call mult     ("tr17","s1  ","s2  ",222)
      call reorder  ("s2  ","s1  ",1324,4)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
c
      call mult6    ("tvc ","s3  ","s1  ",1)
      call add      ("s1  ",1.d0,"tb3 ",1.d0,"tb3 ")
c   diag h.15,17,20,23
      call mult6    ("tr5 ","s3  ","s1  ",1)
      call add      ("s1  ",1.d0,"tb5 ",1.d0,"tb5 ")
c   diag h.13,14,18,21,22,24
      call reorder  ("boxh","s1  ",4213,4)
      call reorder  ("boxh","s2  ",4123,4)
      call add      ("s1  ",2.d0,"s2  ",-1.d0,"s2  ")
      call mult     ("tr16","s2  ","boxu",222)
c
      call mult     ("tr17","s1  ","s2  ",222)
      call add      ("s2  ",-1.d0,"boxu",1.d0,"boxu")
c
      call reorder  ("boxp","s1  ",1243,4)
      call mult     ("tgc ","s1  ","s2  ",222)
      call reorder  ("s2  ","s1  ",3214,4)
      call add      ("s1  ",1.d0,"boxu",1.d0,"boxu")
c
      call mult6    ("boxu","tr5 ","s1  ",1)
      call add      ("s1  ",1.d0,"tb4 ",1.d0,"tb4 ")
c   diag h.3,4,5,10,11,12
      call add      ("bp1 ",2.d0,"bp2 ",-1.d0,"s1  ")
      call mult     ("tr16","s1  ","boxv",222)
c
      call mult     ("tr17","bp1 ","s1  ",222)
      call add      ("s1  ",-1.d0,"boxv",1.d0,"boxv")
c
      call mult     ("tr17","bp2 ","s1  ",222)
      call reorder  ("s1  ","s2  ",3214,4)
      call add      ("s2  ",-1.d0,"boxv",1.d0,"boxv")
c
      call reorder  ("boxv","s2  ",1324,4)
      call mult6    ("s2  ","tr1 ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb1 ",1.d0,"tb1 ")
c   diag h.6,7,8,9
      call reorder  ("boxh","s1  ",4321,4)
      call mult     ("tvg ","s1  ","s3  ",222)
c
      call mult     ("tr18","bp2 ","s1  ",222)
      call reorder  ("s1  ","s2  ",3214,4)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
      call mult6    ("tr2 ","s3  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb2 ",1.d0,"tb2 ")
      call prtime   ("group h ",8,0)
C
C     GROUP I
C
c   diag i.1,3
      call fold4    ("bp1 ","s1  ",2,4)
      call fold4    ("s1  ","s3  ",3,4)
c   diag i.2
      call reorder  ("vpgh","s1  ",1342,4)
      call mult     ("t1v ","s1  ","s2  ",131)
      call reorder  ("s2  ","s1  ",2314,4)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
c   diag i.4
      call mult     ("vvhh","t1g ","s1  ",311)
      call reorder  ("s1  ","s2  ",2413,4)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
c   diag i.5,9
      call add      ("bp1 ",2.d0,"bp2 ",-1.d0,"s1  ")
      call mult     ("tr14","s1  ","s2  ",222)
      call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")     
c   diag i.6
      call reorder  ("vphh","s1  ",1324,4)
      call mult     ("tr15","s1  ","s2  ",222)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
c   diag i.7,11
      call mult     ("tr15","bp2 ","s1  ",222)
      call reorder  ("s1  ","s2  ",3214,4)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
c   diag i.8,12
      call reorder  ("boxh","s1  ",4312,4)
      call mult     ("tum ","s1  ","s2  ",222)
      call reorder  ("s2  ","s1  ",2314,4)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
c   diag i.10
      call mult     ("t1v ","vm8 ","s1  ",131)
      call reorder  ("s1  ","s2  ",1324,4)
      call mult     ("tr19","s2  ","s1  ",222)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
c
      call reorder  ("s3  ","s2  ",3124,4)
      call mult6    ("s2  ","tr6 ","s1  ",1)
      call add      ("s1  ",-1.d0,"tb6 ",1.d0,"tb6 ")
      call prtime   ("group i ",9,0)
c     
      call reord6   ("tb1 ","ttqf",126345)
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
      call prtime   ("ttqf",0,0)
      return
      end
C  *********************************************************************
      subroutine heffp
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)
#include "cfiles.h"
#endif

C 
C     GROUP J
C
c   diag j.1
      call mult     ("tgd ","hb1 ","s1  ",311)
      call reorder  ("s1  ","s2  ",2134,4)
      call mult6    ("tua ","s2  ","s1  ",1)
      call cmult    ("s1  ",-1.d0,"tba ")               
c   diag j.2
      call mult     ("tr2 ","hb1 ","s1  ",311)
      call mult6    ("s1  ","tr20","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tbd ")                 
c   diag j.3
      call mult     ("tr4 ","hb1 ","s1  ",311)
      call mult6    ("s1  ","tvc ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"tbb ")                 
      call prtime   ("group j ",1,0)
C
C     GROUP K
C
      call mult6    ("boxc","tr20","s1  ",1)
      call add      ("s1  ",1.d0,"tbd ",1.d0,"tbd ")
      call prtime   ("group k ",2,0)
C
C     GROUP L
C
c   diag l.2,4
      call mult6    ("boxd","tr4 ","s1  ",1)
      call cmult    ("s1  ",-1.d0,"tbc ")                    
c   diag l.1,3
      call fold4    ("hbp ","s1  ",4,4)
      call reorder  ("s1  ","s3  ",3421,4)
c
      call reorder  ("hbp ","s1  ",2314,4)
      call mult     ("t1g ","s1  ","s2  ",131)
      call reorder  ("s2  ","s1  ",3124,4)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
c
      call mult6    ("tua ","s3  ","s1  ",1)
      call add      ("s1  ",1.d0,"tba ",1.d0,"tba ")
      call prtime   ("group l ",3,0)
C
C     GROUP M
C
c   diag m.1,3,6,12
      call mult6    ("tvc ","boxy","s1  ",1)
      call add      ("s1  ",1.d0,"tbc ",1.d0,"tbc ")
c   diag m.2,4,14,15
      call mult6    ("boxx","tr21","s1  ",1)
      call add      ("s1  ",-1.d0,"tba ",1.d0,"tba ")
c   diag m.7,8,9
      call mult6    ("boxw","tr20","s1  ",1)
      call add      ("s1  ",1.d0,"tbd ",1.d0,"tbd ")
c   diag m.10,11
      call add      ("hbs1",2.d0,"hbs2",-1.d0,"s1  ")
      call mult     ("tr22","s1  ","s3  ",222)
c
      call mult     ("tr18","hbs1","s1  ",222)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
c
      call reorder  ("s3  ","s2  ",3124,4)
      call mult6    ("tua ","s2  ","s1  ",1)
      call add      ("s1  ",1.d0,"tba ",1.d0,"tba ")
c   diag m.13
      call reorder  ("hbs ","s1  ",3412,4)
      call mult     ("tvg ","s1  ","s2  ",222)
      call mult6    ("tr4 ","s2  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tbb ",1.d0,"tbb ")
      call prtime   ("group m ",4,0)
C  
C     GROUP N
C
c   diag n.1,2,3,7
      call fold4    ("hbt ","s1  ",1,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","s3  ",4213,4)
c
      call reorder  ("s1  ","s2  ",2134,4)
      call mult     ("t1g ","s2  ","s1  ",131)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
c
      call fold4    ("hbt ","s1  ",4,4)
      call reorder  ("s1  ","s2  ",4231,4)
      call mult     ("s2  ","t1v ","s1  ",311)
      call reorder  ("s1  ","s2  ",1243,4)
      call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
c
      call mult     ("tr15","hbt1","s1  ",222)
      call reorder  ("s1  ","s2  ",2314,4)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
c
      call reorder  ("s3  ","s2  ",3214,4)
      call mult6    ("tr2 ","s2  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tbd ",1.d0,"tbd ")
c   diag n.5,6
      call reorder  ("boxg","s1  ",1324,4)
      call mult6    ("tr4 ","s1  ","s2  ",1)
      call add      ("s2  ",-1.d0,"tbb ",1.d0,"tbb ")
c   diag n.4
      call reorder  ("hbt ","s1  ",2143,4)
      call mult     ("tgc ","s1  ","s2  ",222)
      call reorder  ("s2  ","s1  ",2134,4)
      call mult6    ("tua ","s1  ","s2  ",1)
      call add      ("s2  ",1.d0,"tba ",1.d0,"tba ")
      call prtime   ("group n ",5,0)
C
C     GROUP O
C
c   diag o.1,2,3,16,17,18
      call mult6    ("tvc ","boxu","s1  ",1)
      call add      ("s1  ",1.d0,"tbc ",1.d0,"tbc ")
c   diag o.4,5,6,13,14,15
      call mult6    ("boxv","tr21","s1  ",1)
      call add      ("s1  ",-1.d0,"tba ",1.d0,"tba ")
c   diag o.7,8,9,10,11,12
      call add      ("bp1 ",2.d0,"bp2 ",-1.d0,"s1  ")
      call mult     ("tr22","s1  ","s3  ",222)
c
      call mult     ("tr18","bp1 ","s1  ",222)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
c
      call reorder  ("boxh","s1  ",3412,4)
      call mult     ("tvg ","s1  ","s2  ",222)
      call reorder  ("s2  ","s1  ",4213,4)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
c
      call reorder  ("s3  ","s2  ",3214,4)
      call mult6    ("tr2 ","s2  ","s1  ",1)
      call add      ("s1  ",-1.d0,"tbd ",1.d0,"tbd ")
      call prtime   ("group o ",6,0)
c
      call reord6   ("tba ","ttpf",126345)
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
c
      call reord6   ("ttpf","s1  ",123546)
      call reord7   ("s1  ","s2  ",213456)
      call add      ("ttpf",1.d0,"s2  ",1.d0,"ttpf")
      call prtime   ("ttpf",0,0)
      return
      end 
C  *********************************************************************
C  *********************************************************************
