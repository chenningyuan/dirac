!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine sixh
!---------------------------------------------------------------------
!
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use cint_mod
      use dat_mod
      use idv_mod
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
!     common /idv/ iv(120),nvalo
#include "idv.h"
#endif

      character*256 filenm,filenm2,filenm6,                             &
 filenm51,filenm53,filenm55,filenm56,filenm57,filenm60,filenm62
      real*4 cputime
!     iopt(1) - print.
!     iopt(5) - print only elements .gt. 10**(-iopt(5))
!     iopt(8) - fock.
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
      filenm56(1:1)=' '
      call getenv('tape56',filenm56)
      if (filenm56(1:1).eq.' ') filenm56='tape56'
      open(56,file=filenm56,form='UNFORMATTED')
      filenm57(1:1)=' '
      call getenv('tape57',filenm57)
      if (filenm57(1:1).eq.' ') filenm57='tape57'
      open(57,file=filenm57,form='UNFORMATTED')
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
      write (6,*) '   $$$$$$$$$   sixh   $$$$$$$$$  '
      csdstim=cputime(real(0.))
      csdltim=csdstim
      call rdin
      call rdsym
      call prtime ("rdin    ",0,0)
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      if (iopt(16).eq.1) then 
        call reads    (51,0,"s1  ")
        call reorder  ("s1  ","tc  ",3412,4)
        call rdfile   (53,0,"thv ",-1)
        call rdfile   (53,2,"thvf",-1)
        if (iopt(11).ge.2) then
          call rdfile (55,0,"thu ",-1)
          call rdfile (55,1,"thuf",-1)
        end if
        if (iopt(11).ge.3) then
          call rdfile (60,0,"thwf",-1)
        end if
        if (iopt(11).ge.4) then
          call rdfile (62,0,"thxf",-1)
        end if
        if (iopt(11).ge.5) then
          call rdfile (62,1,"thyf",-1)
        end if
        if (iopt(11).eq.6) then
          call rdfile (62,2,"thzf",-1)
        end if
      else if (iopt(16).eq.2) then
        call reads    (51,0,"s1  ")
        call reorder  ("s1  ","tc  ",3412,4)
        call rdfile   (53,0,"thv ",-1)
        call rdfile   (53,2,"thvf",-1)
      end if
      if (iopt(11).ge.5) then
        nvalo=0
        do 10 j=1,120
          if (icv(j).eq.1) then
            nvalo=nvalo+1
            iv(nvalo)=j      
          end if
10      continue
      end if
      call ordv_sixh
      call ordta_sixh
      call diags1_sixh
      if (iopt(11).ge.2) then
        call diags2_sixh
      end if
      if (iopt(11).ge.3) then
        call diags3_sixh
      end if
      if (iopt(11).ge.4) then
        call diags4_sixh
      end if
      egr0=0.d0
      call eneru2 (egr0,1,2,max(1,iopt(1)))
      if (iopt(16).eq.1) then
        if (iopt(11).ge.2) then
          call eneru2 (egr0,2,2,max(1,iopt(1)))
        end if
        if (iopt(11).ge.3) then
          call eneru2 (egr0,3,2,max(1,iopt(1)))
        end if
        if (iopt(11).ge.4) then
          call eneru2 (egr0,4,2,max(1,iopt(1)))
        end if
        if (iopt(11).ge.5) then
          call eneru2 (egr0,5,2,max(1,iopt(1)))
        end if
        if (iopt(11).eq.6) then
          call eneru2 (egr0,6,2,max(1,iopt(1)))
        end if
      else if (iopt(16).eq.2) then
        call wrtfile (53,2,"thvf")
        call wrtfile (53,3,"t31 ")
      end if
      call prt6 ("t31 ","t31 ",2)
      call prtime ("end ",0,0)
!     stop
      end

!  *******************************************************************
      subroutine ordv_sixh
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

      call readv    (nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call prtime   ("readv   ",nint,0)
      call seltyp   (3332,0,"ppph",4)
      call seltyp   (3322,0,"pphh",4)
      call seltyp   (3222,0,"phhh",4)
      call seltyp   (2332,1000,"gpph",4)
      call seltyp   (2323,1000,"gphp",4)
      call seltyp   (2222,1000,"ghhh",4)
      if (iopt(11).ge.2) then
        call seltyp   (2223,0010,"hhgp",4)
        call seltyp   (2223,1100,"gghp",4)
        call seltyp   (2222,1111,"gggg",4)
      end if
      if (iopt(11).ge.3) then
        call seltyp   (2222,0011,"hhgg",4)
        call seltyp   (2223,1010,"ghgp",4)
        call seltyp   (2232,1001,"ghpg",4)
        call seltyp   (2233,1100,"ggpp",4)
      end if
      if (iopt(11).ge.4) then
        call seltyp   (2222,0111,"hggg",4)
        call seltyp   (2223,1110,"gggp",4)
      end if
      call prtime   ("ordv",0,0)
      return
      end 
!  *******************************************************************
      subroutine ordta_sixh
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
      call reorder  ("tc  ","tc3 ",1432,4)
      call reorder  ("thv ","tv3 ",1432,4)
      call reorder  ("thv ","tv6 ",1243,4)
      return
      end
!  *******************************************************************
      subroutine diags1_sixh
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
!  diags 1a,2c
      call mult6    ("thv ","phhh","s1  ",1)
      call mult6    ("gphp","tc3 ","s2  ",1)
      call add      ("s1  ",-1.d0,"s2  ",1.d0,"s1  ")
      call reord6   ("s1  ","t31 ",124356)
!  diags 1b,2b
      call mult6    ("tv6 ","phhh","s1  ",1)
      call reorder  ("gpph","s2  ",1243,4)
      call mult6    ("s2  ","tc3 ","s3  ",1)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s1  ")
      call reord6   ("s1  ","s2  ",124536)
      call add      ("s2  ",1.d0,"t31 ",1.d0,"t31 ")
!  diag 1c
      call mult6    ("ghhh","tc  ","s1  ",1)
      call reord6   ("s1  ","s2  ",124356)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",125463)
      call add      ("s2  ",-1.d0,"t31 ",1.d0,"t31 ")
!  diag 2a
      call reorder  ("ppph","s1  ",1243,4)
      call mult6    ("tv3 ","s1  ","s2  ",1)
      call reord6   ("s2  ","s1  ",124356)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",125436)
      call add      ("s1  ",1.d0,"t31 ",1.d0,"t31 ")
!
      call reord6   ("t31 ","s1  ",123465)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("t31 ",1.d0,"s2  ",1.d0,"t31 ")
      call dveps6   ("t31 ")
      call cmult    ("t31 ",-1.d0,"t31 ")
!
!  diags e1,e2
      call fold6    ("t31 ","s1  ",4,1)
      call reord6   ("s1  ","s2  ",124356)
      call reord7   ("s2  ","s1  ",132456)
      call fold6    ("t31 ","s2  ",5,1)
      call reord6   ("s2  ","s3  ",125346)
      call reord7   ("s3  ","s2  ",132456)
      call add      ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
      call reorder  ("pphh","s2  ",1243,4)
      call add      ("pphh",2.d0,"s2  ",-1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","s3  ",0)
      call add      ("s3  ",1.d0,"thvf",1.d0,"thvf")
      call prtime   ("diags1  ",0,0)
      return
      end
!  *******************************************************************
      subroutine diags2_sixh
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
!  diags 1a,2c,f1
      call mult6    ("thu ","phhh","s3  ",1)
      call mult6    ("gghp","tc3 ","s1  ",1)
      call add      ("s1  ",1.d0,"s3  ",-1.d0,"s3  ")
      call reorder  ("thv ","s1  ",2341,4)
      call mult6    ("gggg","s1  ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"s2  ")
      call fold6    ("s3  ","s2  ",idum,2)
      call reord6   ("s3  ","t32 ",124356)
!  diags 1b,2b    
      call mult6    ("thv ","ghhh","s3  ",1)
      call mult6    ("gphp","tv3 ","s1  ",1)
      call add      ("s1  ",1.d0,"s3  ",-1.d0,"s3  ")
      call reord6   ("s3  ","s1  ",124365)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("s2  ",1.d0,"t32 ",1.d0,"t32 ")
!  diags 1c,2a    
      call mult6    ("tv6 ","ghhh","s3  ",1)
      call reorder  ("gpph","s1  ",1243,4)
      call mult6    ("s1  ","tv3 ","s2  ",1)
      call add      ("s2  ",1.d0,"s3  ",-1.d0,"s3  ")
      call reord6   ("s3  ","s1  ",124563)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("s2  ",1.d0,"t32 ",1.d0,"t32 ")
!
      call reord6   ("t32 ","s1  ",123546)
      call reord7   ("s1  ","s2  ",213456)
      call add      ("t32 ",1.d0,"s2  ",1.d0,"t32 ")
      call dveps6   ("t32 ")
      call cmult    ("t32 ",-1.d0,"t32 ")
!
!  diag e1
      call fold6    ("t32 ","s1  ",4,1)
      call reord6   ("s1  ","s2  ",124356)
      call reorder  ("hhgp","s1  ",3412,4)
      call reorder  ("hhgp","s3  ",3421,4)
      call add      ("s1  ",-2.d0,"s3  ",1.d0,"s1  ")
      call mult6    ("s2  ","s1  ","dum2",2)
!  diag e2
      call fold6    ("t31 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord6   ("s2  ","s1  ",124356)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",125436)
      call reorder  ("gpph","s2  ",1324,4)
      call add      ("gpph",2.d0,"s2  ",-1.d0,"s2  ")
      call mult6    ("s2  ","s1  ","s3  ",2)
      call reorder  ("s3  ","s1  ",1243,4)
      call add      ("s1  ",1.d0,"dum2",1.d0,"dum2")
!  diag e3
      call fold6    ("t32 ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",126345)
      call reorder  ("hhgp","s1  ",3412,4)
      call mult6    ("s2  ","s1  ","s3  ",2)
      call reorder  ("s3  ","s1  ",1243,4)
      call add      ("s1  ",1.d0,"dum2",1.d0,"dum2")
!  diag e4
      call fold6    ("t31 ","s1  ",5,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",126345)
      call mult6    ("gpph","s1  ","s2  ",2)
      call reorder  ("s2  ","s1  ",1243,4)
      call add      ("s1  ",-1.d0,"dum2",1.d0,"dum2")
!
      call sqr      ("dum2","s1  ")
      call add      ("dum2",1.d0,"thuf",1.d0,"thuf")
      call prtime   ("diags2  ",0,0)
      return
      end
!  *******************************************************************
      subroutine diags3_sixh
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
!  diags 1,2,f1
      call mult6    ("thu ","ghhh","s3  ",1)
      call mult6    ("gghp","tv3 ","s1  ",1)
      call add      ("s1  ",1.d0,"s3  ",-1.d0,"s3  ")
      call reorder  ("thu ","s1  ",1432,4)
      call mult6    ("gggg","s1  ","s2  ",1)
      call cmult    ("s2  ",-1.d0,"s2  ")
      call fold6    ("s3  ","s2  ",idum,2)
      call reord6   ("s3  ","t33 ",124356)
!
      call sqr3     ("t33 ")
      call zero     ("t33 ")
      call dveps6   ("t33 ")
      call cmult    ("t33 ",-1.d0,"t33 ")
!
!  diag e1
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord6   ("s2  ","s1  ",124536)
      call reorder  ("ghgp","s2  ",1342,4)
      call reorder  ("ghpg","s3  ",1432,4)
      call add      ("s2  ",2.d0,"s3  ",-1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","s3  ",3)
      call reord6   ("s3  ","dum3",125346)
!  diag e2
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124635)
      call reorder  ("ghgp","s2  ",1342,4)
      call mult6    ("s1  ","s2  ","s3  ",3)
      call reord6   ("s3  ","s2  ",125346)
      call add      ("s2  ",-1.d0,"dum3",0.5d0,"dum3")
!  diag e3
      call reorder  ("ghpg","s2  ",1432,4)
      call mult6    ("s1  ","s2  ","s3  ",3)
      call reord6   ("s3  ","s1  ",125364)
      call add      ("s1  ",-1.d0,"dum3",1.d0,"dum3")
!  diag e4
      call fold6    ("t31 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call fold6    ("s2  ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",126345)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",126543)
      call mult6    ("ggpp","s2  ","s1  ",3)
      call add      ("s1  ",0.5d0,"dum3",1.d0,"dum3")
!  diag e5
      call fold6    ("t33 ","s1  ",4,1)
      call reorder  ("hhgg","s2  ",3412,4)
      call mult6    ("s1  ","s2  ","s3  ",3)
      call add      ("s3  ",0.5d0,"dum3",1.d0,"dum3")
!
      call sqr3     ("dum3")
      call add      ("dum3",1.d0,"thwf",1.d0,"thwf")
      call prtime   ("diags3  ",0,0)
      return
      end 
!  *******************************************************************
      subroutine diags4_sixh
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
!  diag e1
      call fold6    ("t33 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reorder  ("hggg","s1  ",2341,4)
      call mult8    ("s2  ","s1  ","s3  ",1)
      call reordx   ("s3  ","dum4",12364578)
!  diag e2
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call fold6    ("s2  ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",124563)
      call mult8    ("s2  ","gggp","s1  ",1)
      call reordx   ("s1  ","s2  ",12673485)
      call add      ("s2  ",0.5d0,"dum4",-0.5d0,"dum4")
!
      call sqr4     ("dum4")
      call add      ("dum4",1.d0,"thxf",1.d0,"thxf")
      call prtime   ("diags4  ",0,0)
      return
      end 
!  *******************************************************************
