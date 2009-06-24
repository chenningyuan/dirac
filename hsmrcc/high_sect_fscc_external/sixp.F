!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine sixp
!---------------------------------------------------------------------
!
c   Reference : A.Haque, U.Kaldor; CPL 120,261 (1985)
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

      character*256 filenm,filenm2,filenm3,filenm6,
     1filenm51,filenm52,filenm53,filenm56,filenm57
      real*4 cputime

c     iopt(1) - print.
c     iopt(5) - print only elements .gt. 10**(-iopt(5))
c     iopt(8) - fock.
c     iopt(15)- -1 indicates (1,1) sector                    
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
      filenm57(1:1)=' '
      call getenv('tape57',filenm57)
      if (filenm57(1:1).eq.' ') filenm57='tape57'
      open(57,file=filenm57,form='UNFORMATTED')
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')
      write (6,*) '   $$$$$$$$$   sixp   $$$$$$$$$  '
      csdstim=cputime(real(0.))
      csdltim=csdstim
      call rdin
      call rdsym
      call prtime ("rdin    ",0,0)
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      call rdfile (51,0,"tc  ",-1)
      call rdfile (52,0,"tv  ",-1)
      call rdfile (52,1,"tvf ",-1)
      call rdfile (56,0,"thp ",-1)
      call rdfile (56,1,"thpf",-1)
      call rdfile (57,0,"thq ",-1)
      call rdfile (57,1,"thqf",-1)
      call rdfile (52,3,"t3v ",-1)
      call rdfile (53,0,"tg  ",-1)
      call rdfile (53,2,"tgf ",-1)
      call rdfile (53,3,"t3g ",-1)
      call ordv3
      call diags
      call ehp(iopt(1)+1)
      write (6,*) nstrt(1,nflcur),'words used in igen,',
     1nstrt(2,nflcur),' words used in gen'
      call prtime ("end ",0,0)
!     stop
      end
C  *********************************************************************
      subroutine ordv3
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
!     common /dat/  norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)
#include "fock.h"
#endif

      call readv     (nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call prtime    ("readv   ",nint,0)
      call seltyp    (3232,1111,"vgvg",0)
      call seltyp    (3223,1111,"vggv",0)
      call seltyp    (2322,1000,"gphh",4)
      call seltyp    (2222,0010,"hhgh",4)
      call seltyp    (2223,0,"hhhp",4)
      call seltyp    (2232,0,"hhph",4)
      call seltyp    (3223,0010,"phgp",4)
      call seltyp    (3232,0001,"phpg",4)
      call seltyp    (3233,0,"phpp",4)
      call seltyp    (3323,0001,"pphv",4)
      call seltyp    (3322,0010,"ppgh",4)
      call seltyp    (3222,1010,"vhgh",4)
      call seltyp    (3222,1001,"vhhg",4)
      call seltyp    (3223,1000,"vhhp",4)
      call seltyp    (3232,1000,"vhph",4)
      call seltyp    (3323,1010,"vpgp",4)
      call seltyp    (3322,1000,"vphh",4)
      call seltyp    (3332,1001,"vppg",4)
      call seltyp    (3333,1000,"vppp",4)
      call reorder   ("phgp","ve1 ",3421,4)
      call reorder   ("ppgh","ve2 ",3214,4)
      call reorder   ("ppgh","s1  ",3124,4)
      call add       ("s1  ",2.d0,"ve2 ",-1.d0,"ve3 ")
      call reorder   ("pphv","ve5 ",4123,4)
      call reorder   ("pphv","s1  ",4213,4)
      call add       ("ve5 ",2.d0,"s1  ",-1.d0,"ve4 ")
      call reorder   ("gphh","s1  ",1243,4)
      call add       ("gphh",-2.d0,"s1  ",1.d0,"ve6 ")
      call reorder   ("vhhp","ve18",1243,4)
      call reorder   ("vphh","ve22",1432,4)
      call reorder   ("vphh","s1  ",1342,4)
      call add       ("s1  ",-2.d0,"ve22",1.d0,"ve21")
      call reorder   ("tc  ","tc1 ",1432,4)
      call reorder   ("tv  ","tv1 ",1243,4)
      call reorder   ("tv  ","tv2 ",1342,4)
      call reorder   ("tv  ","tv3 ",2341,4)
      call reorder   ("tg  ","tg1 ",1243,4)
      call reorder   ("tg  ","tg3 ",4321,4)
      call reorder   ("thp ","thp2",1342,4)
      call prtime    ("ordv",0,0)
      return
      end 
C  *******************************************************************
      subroutine diags
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
C    
C     T3P DIAGRAMS
C
c     diags 2,15 
      call mult6     ("tv  ","ve1 ","s1  ",1)
      call mult6     ("vhph","tg1 ","s2  ",1)
      call add       ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
      call reord6    ("s1  ","t3p ",126345)
      call prtime    ("dgs 2,15",2,0)
c     diags 3,12
      call reorder   ("phpg","s3  ",4321,4)
      call mult6     ("tv  ","s3  ","s1  ",1)
      call mult6     ("vhph","tg  ","s2  ",1)
      call add       ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
      call reord6    ("s1  ","s2  ",126345)
      call reord7    ("s2  ","s1  ",132456)
      call add       ("s1  ",1.d0,"t3p ",1.d0,"t3p ")
      call prtime    ("dgs 3,12",4,0)
c     diags 4,17,f1
      call reorder   ("phpp","s2  ",2341,4)
      call reorder   ("thp ","s3  ",1243,4)
      call mult6     ("s3  ","s2  ","s1  ",1)
      call reorder   ("vhhg","s3  ",1243,4)
      call mult6     ("s3  ","tc1 ","s2  ",1)
      call add       ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
      call reorder   ("vgvg","s3  ",1243,4)
      call mult6     ("s3  ","tv3 ","s2  ",1)
      call fold6     ("s1  ","s2  ",idum,2)
      call reord6    ("s1  ","s2  ",124536)
      call add       ("s2  ",1.d0,"t3p ",1.d0,"t3p ")
      call prtime    ("dgs 4,17",7,0)
c     diags 6,11
      call reorder   ("phpg","s3  ",4321,4)
      call mult6     ("tv1 ","s3  ","s1  ",1)
      call mult6     ("ve18","tg  ","s2  ",1)
      call add       ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
      call reord6    ("s1  ","s2  ",126543)
      call reord7    ("s2  ","s1  ",132456)
      call add       ("s1  ",1.d0,"t3p ",1.d0,"t3p ")
      call prtime    ("dgs 6,11",9,0)
c     diags 7,16
      call reorder   ("vppp","s3  ",1342,4)
      call reorder   ("tg  ","s2  ",1342,4)
      call mult6     ("s3  ","s2  ","s1  ",1)
      call reorder   ("hhgh","s3  ",3124,4)
      call mult6     ("tv2 ","s3  ","s2  ",1)
      call add       ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
      call reord6    ("s1  ","s2  ",125643)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",124356)
      call add       ("s2  ",1.d0,"t3p ",1.d0,"t3p ")
      call prtime    ("dgs 7,16",11,0)
c     diags 8,13,f2
      call reorder   ("vppg","s3  ",1342,4)
      call mult6     ("s3  ","tc  ","s1  ",1)
      call reorder   ("hhhp","s3  ",2143,4)
      call mult6     ("thp2","s3  ","s2  ",1)
      call add       ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
      call reorder   ("vgvg","s3  ",1342,4)
      call mult6     ("s3  ","tg3 ","s2  ",1)
      call cmult     ("s2  ",-1.d0,"s2  ")
      call fold6     ("s1  ","s2  ",idum,2)
      call reord6    ("s1  ","s2  ",125436)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",124356)
      call add       ("s2  ",1.d0,"t3p ",1.d0,"t3p ")
      call prtime    ("dgs 8,13",14,0)
c
      call dveps6    ("t3p ")
C
C     T3Q DIAGRAMS
C 
c     diags 1,18,f4  
      call reorder   ("phpp","s3  ",2341,4)
      call mult6     ("thq ","s3  ","s1  ",1)
      call reorder   ("tc  ","s3  ",2341,4)
      call mult6     ("vhgh","s3  ","s2  ",1)
      call add       ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
      call mult6     ("vggv","tv3 ","s2  ",1)
      call fold6     ("s1  ","s2  ",idum,2)
      call reord6    ("s1  ","t3q ",124356)
      call prtime    ("dgs 1,18",17,0)
c     diags 5,14    
      call mult6     ("tv1 ","ve1 ","s1  ",1)
      call mult6     ("ve18","tg1 ","s2  ",1)
      call add       ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
      call reord6    ("s1  ","s2  ",126435)
      call add       ("s2  ",1.0d0,"t3q ",1.d0,"t3q ")
      call prtime    ("dgs 5,14",19,0)
c     diags  9,10,f3
      call reorder   ("vpgp","s3  ",1342,4)
      call mult6     ("s3  ","tc  ","s1  ",1)
      call reorder   ("thq ","s3  ",1342,4)
      call mult6     ("s3  ","hhph","s2  ",1)
      call add       ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
      call reorder   ("vggv","s3  ",1342,4)
      call mult6     ("s3  ","tg3 ","s2  ",1)
      call cmult     ("s2  ",-1.d0,"s2  ")
      call fold6     ("s1  ","s2  ",idum,2)
      call reord6    ("s1  ","s2  ",125436)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",124356)
      call add       ("s2  ",1.d0,"t3q ",1.d0,"t3q ")
      call prtime    ("dgs 9,10",22,0)
c
      call reord6    ("t3q ","s1  ",123465)
      call reord7    ("s1  ","s2  ",132456)       
      call add       ("s2  ",1.d0,"t3q ",1.d0,"t3q ")
      call dveps6    ("t3q ")
C
C     CONTRIBUTIONS TO THPF
C
c     diag p1
      call fold6     ("t3v ","s1  ",3,1)
      call fold6     ("s1  ","s2  ",4,1)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",124653) 
      call mult6     ("s2  ","ve3 ","s1  ",2)
      call add       ("s1  ",1.d0,"thpf",1.d0,"thpf")
      call prtime    ("diag p1 ",23,0)
c     diag p2
      call fold6     ("t3v ","s1  ",3,1)
      call fold6     ("s1  ","s2  ",5,1)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125463)
      call mult6     ("s2  ","ve2 ","s1  ",2)
      call add       ("s1  ",-1.d0,"thpf",1.d0,"thpf")
      call prtime    ("diag p2 ",24,0)
c     diag p3
      call fold6     ("t3g ","s1  ",3,1)
      call fold6     ("s1  ","s2  ",4,1)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",124653) 
      call mult6     ("s2  ","ve21","s1  ",2)
      call reorder   ("s1  ","s2  ",4321,4)
      call add       ("s2  ",1.d0,"thpf",1.d0,"thpf")
      call prtime    ("diag p3 ",25,0)
c     diag p4
      call fold6     ("t3g ","s1  ",3,1)
      call fold6     ("s1  ","s2  ",5,1)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125463)
      call mult6     ("s2  ","ve22","s1  ",2)
      call reorder   ("s1  ","s2  ",4321,4)
      call add       ("s2  ",1.d0,"thpf",1.d0,"thpf")
      call prtime    ("diag p4 ",26,0)
c     diag p5
      call fold6     ("t3p ","s1  ",2,1)
      call reord6    ("s1  ","s2  ",125643)
      call mult6     ("s2  ","ve4 ","s1  ",2)
      call reorder   ("s1  ","s2  ",1243,4)
      call add       ("s2  ",1.d0,"thpf",1.d0,"thpf")
      call prtime    ("diag p5 ",27,0)
c     diag p6
      call fold6     ("t3p ","s1  ",3,1)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",125643)
      call mult6     ("s1  ","ve5 ","s2  ",2)
      call reorder   ("s2  ","s1  ",1243,4)
      call add       ("s1  ",-1.d0,"thpf",1.d0,"thpf")
      call prtime    ("diag p6 ",28,0)
c     diag p7
      call fold6     ("t3p ","s1  ",4,1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125634)
      call mult6     ("s2  ","ve6 ","s1  ",2)
      call reorder   ("s1  ","s2  ",1423,4)
      call add       ("s2  ",1.d0,"thpf",1.d0,"thpf")
      call prtime    ("diag p7 ",29,0)
c     diag p8
      call fold6     ("t3p ","s1  ",6,1)
      call reord6    ("s1  ","s2  ",125346)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",126534)
      call mult6     ("s2  ","gphh","s1  ",2)
      call reorder   ("s1  ","s2  ",1432,4)
      call add       ("s2  ",1.d0,"thpf",1.d0,"thpf")
      call prtime    ("diag p8 ",30,0)
C
C     CONTRIBUTIONS TO THQF
C
c     diag q1
      call fold6     ("t3v ","s1  ",3,1)
      call fold6     ("s1  ","s2  ",6,1)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",126453)
      call mult6     ("s2  ","ve3 ","s1  ",2)
      call reorder   ("s1  ","s2  ",1243,4)
      call add       ("s2  ",1.d0,"thqf",1.d0,"thqf")
      call prtime    ("diag q1 ",31,0)
c     diag q2
      call fold6     ("t3v ","s1  ",2,1)
      call fold6     ("s1  ","s2  ",6,1)
      call reord6    ("s2  ","s1  ",126543)
      call mult6     ("s1  ","ve2 ","s2  ",2)
      call reorder   ("s2  ","s1  ",1243,4)
      call add       ("s1  ",-1.d0,"thqf",1.d0,"thqf")
      call prtime    ("diag q2 ",32,0)
c     diag q3
      call fold6     ("t3g ","s1  ",2,1)
      call fold6     ("s1  ","s2  ",5,1)
      call reord6    ("s2  ","s1  ",125463)
      call mult6     ("s1  ","ve21","s2  ",2)
      call reorder   ("s2  ","s1  ",4312,4)
      call add       ("s1  ",1.d0,"thqf",1.d0,"thqf")
      call prtime    ("diag q3 ",33,0)
c     diag q4
      call fold6     ("t3g ","s1  ",2,1)
      call fold6     ("s1  ","s2  ",6,1)
      call reord6    ("s2  ","s1  ",126543)
      call mult6     ("s1  ","ve22","s2  ",2)
      call reorder   ("s2  ","s1  ",4312,4)
      call add       ("s1  ",1.d0,"thqf",1.d0,"thqf")
      call prtime    ("diag q4 ",34,0)
c     diag q5
      call fold6     ("t3q ","s1  ",2,1)
      call reord6    ("s1  ","s2  ",124653)
      call mult6     ("s2  ","ve4 ","s1  ",2)
      call add       ("s1  ",1.d0,"thqf",1.d0,"thqf")
      call prtime    ("diag q5 ",35,0)
c     diag q6
      call fold6     ("t3p ","s1  ",3,1)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",125463)
      call mult6     ("s1  ","ve5 ","s2  ",2)
      call add       ("s2  ",-1.d0,"thqf",1.d0,"thqf")
      call prtime    ("diag q6 ",36,0)
c     diag q7
      call fold6     ("t3q ","s1  ",5,1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125634)
      call mult6     ("s2  ","ve6 ","s1  ",2)
      call reorder   ("s1  ","s2  ",1423,4)
      call add       ("s2  ",1.d0,"thqf",1.d0,"thqf")
      call prtime    ("diag q7 ",37,0)
c     diag q8
      call fold6     ("t3p ","s1  ",6,1)
      call reord6    ("s1  ","s2  ",126345)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",126543)
      call mult6     ("s2  ","gphh","s1  ",2)
      call reorder   ("s1  ","s2  ",1432,4)
      call add       ("s2  ",1.d0,"thqf",1.d0,"thqf")
      call prtime    ("diag q8 ",38,0)
      return
      end
C  *******************************************************************
