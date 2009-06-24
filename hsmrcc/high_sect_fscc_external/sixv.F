!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine sixv
!---------------------------------------------------------------------
!
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

      character*256 filenm,filenm2,filenm6,
     1filenm51,filenm52,filenm54,filenm56,filenm57,filenm59,filenm61
      real*4 cputime
c     iopt(1) - print.
c     iopt(5) - print only elements .gt. 10**(-iopt(5))
c     iopt(8) - fock.
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
      filenm59(1:1)=' '
      call getenv('tape59',filenm59)
      if (filenm59(1:1).eq.' ') filenm59='tape59'
      open(59,file=filenm59,form='UNFORMATTED')
      filenm61(1:1)=' '
      call getenv('tape61',filenm61)
      if (filenm61(1:1).eq.' ') filenm61='tape61'
      open(61,file=filenm61,form='UNFORMATTED')
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')
      write (6,*) '   $$$$$$$$$   sixv   $$$$$$$$$  '
      csdstim=cputime(real(0.))
      csdltim=csdstim
      call rdin
      call rdsym
      call prtime ("rdin    ",0,0)
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      if (iopt(16).eq.1) then 
        call rdfile (51,0,"tc  ",-1)
        call rdfile (52,0,"tv  ",-1)
        call rdfile (52,2,"tvf ",-1)
        if (iopt(11).ge.2) then
          call rdfile (54,0,"tu  ",-1)
          call rdfile (54,1,"tuf ",-1)
        end if
        if (iopt(11).ge.3) then
          call rdfile (59,0,"twf ",-1)
        end if
        if (iopt(11).ge.4) then
          call rdfile (61,0,"txf ",-1)
        end if
        if (iopt(11).ge.5) then
          call rdfile (61,1,"tyf ",-1)
        end if
        if (iopt(11).eq.6) then
          call rdfile (61,2,"tzf ",-1)
        end if
      else if (iopt(16).eq.2) then
        call rdfile (51,0,"tc  ",-1)
        call rdfile (52,0,"tv  ",-1)
        call rdfile (52,2,"tvf ",-1)
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
      call ordv_sixv
      call ordta_sixv
      call diags1
      if (iopt(11).ge.2) then
        call diags2
      end if
      if (iopt(11).ge.3) then
        call diags3
      end if
      if (iopt(11).ge.4) then
        call diags4
      end if
      egr0=0.d0
      call eneru2 (egr0,1,3,max(1,iopt(1)))
      if (iopt(16).eq.1) then
        if (iopt(11).ge.2) then
          call eneru2 (egr0,2,3,max(1,iopt(1)))
        end if
        if (iopt(11).ge.3) then
          call eneru2 (egr0,3,3,max(1,iopt(1)))
        end if
        if (iopt(11).ge.4) then
          call eneru2 (egr0,4,3,max(1,iopt(1)))
        end if
        if (iopt(11).ge.5) then
          call eneru2 (egr0,5,3,max(1,iopt(1)))
        end if
        if (iopt(11).eq.6) then
          call eneru2 (egr0,6,3,max(1,iopt(1)))
        end if
      else if (iopt(16).eq.2) then
        call wrtfile (52,2,"tvf ")
        call wrtfile (52,3,"t31 ")
      end if
      call prtime ("end ",0,0)
!     stop
      end
C  *******************************************************************
      subroutine ordv_sixv
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
!
      call readv    (nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call prtime   ("readv   ",nint,0)
      call seltyp   (2223,0,"hhhp",4)
      call seltyp   (2233,0,"hhpp",4)
      call seltyp   (2333,0,"hppp",4)
      call seltyp   (3223,1000,"vhhp",4)
      call seltyp   (3232,1000,"vhph",4)
      call seltyp   (3333,1000,"vppp",4)
      if (iopt(11).ge.2) then
        call seltyp   (3332,0010,"ppvh",4)
        call seltyp   (3332,1100,"vvph",4)
        call seltyp   (3333,1111,"vvvv",4)
        call selt1    (23,0,"hp  ",4)
      end if
      if (iopt(11).ge.3) then
        call seltyp   (3333,0011,"ppvv",4)
        call seltyp   (3332,1010,"vpvh",4)
        call seltyp   (3323,1001,"vphv",4)
        call seltyp   (3322,1100,"vvhh",4)
        call selt1    (32,10,"vh  ",4)
        call selt1    (33,10,"vp  ",4)
      end if
      if (iopt(11).ge.4) then
        call seltyp   (3333,0111,"pvvv",4)
        call seltyp   (3332,1110,"vvvh",4)
      end if
      call prtime   ("ordv",0,0)
      return
      end 
C  *******************************************************************
      subroutine ordta_sixv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
      call reorder  ("tc  ","tc3 ",1432,4)
      call reorder  ("tv  ","tv3 ",1432,4)
      call reorder  ("tv  ","tv6 ",1243,4)
      return
      end
C  *******************************************************************
      subroutine diags1
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
c  diags 1a,2c
      call mult6    ("tv  ","hppp","s1  ",1)
      call mult6    ("vhph","tc3 ","s2  ",1)
      call add      ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
      call reord6   ("s1  ","t31 ",124356)
c  diags 1b,2b
      call mult6    ("tv6 ","hppp","s1  ",1)
      call reorder  ("vhhp","s2  ",1243,4)
      call mult6    ("s2  ","tc3 ","s3  ",1)
      call add      ("s1  ",1.d0,"s3  ",-1.d0,"s1  ")
      call reord6   ("s1  ","s2  ",124536)
      call add      ("s2  ",1.d0,"t31 ",1.d0,"t31 ")
c  diag 1c
      call mult6    ("vppp","tc  ","s1  ",1)
      call reord6   ("s1  ","s2  ",124356)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",125463)
      call add      ("s2  ",1.d0,"t31 ",1.d0,"t31 ")
c  diag 2a
      call reorder  ("hhhp","s1  ",1243,4)
      call mult6    ("tv3 ","s1  ","s2  ",1)
      call reord6   ("s2  ","s1  ",124356)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",125436)
      call add      ("s1  ",-1.d0,"t31 ",1.d0,"t31 ")
c
      call reord6   ("t31 ","s1  ",123465)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("t31 ",1.d0,"s2  ",1.d0,"t31 ")
      call dveps6   ("t31 ")
c
c  diags e1,e2
      call fold6    ("t31 ","s1  ",4,1)
      call reord6   ("s1  ","s2  ",124356)
      call reord7   ("s2  ","s1  ",132456)
      call fold6    ("t31 ","s2  ",5,1)
      call reord6   ("s2  ","s3  ",125346)
      call reord7   ("s3  ","s2  ",132456)
      call add      ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
      call reorder  ("hhpp","s2  ",1243,4)
      call add      ("hhpp",2.d0,"s2  ",-1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","s3  ",0)
      call add      ("s3  ",1.d0,"tvf ",1.d0,"tvf ")
      call prtime   ("diags1  ",0,0)
      return
      end
C  *******************************************************************
      subroutine diags2
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
c  diags 1a,2c,f1
      call mult6    ("tu  ","hppp","s3  ",1)
      call mult6    ("vvph","tc3 ","s1  ",1)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
      call reorder  ("tv  ","s1  ",2341,4)
      call mult6    ("vvvv","s1  ","s2  ",1)
      call fold6    ("s3  ","s2  ",idum,2)
      call reord6   ("s3  ","t32 ",124356)
c  diags 1b,2b    
      call mult6    ("tv  ","vppp","s3  ",1)
      call mult6    ("vhph","tv3 ","s1  ",1)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
      call reord6   ("s3  ","s1  ",124365)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("s2  ",1.d0,"t32 ",1.d0,"t32 ")
c  diags 1c,2a    
      call mult6    ("tv6 ","vppp","s3  ",1)
      call reorder  ("vhhp","s1  ",1243,4)
      call mult6    ("s1  ","tv3 ","s2  ",1)
      call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
      call reord6   ("s3  ","s1  ",124563)
      call reord7   ("s1  ","s2  ",132456)
      call add      ("s2  ",1.d0,"t32 ",1.d0,"t32 ")
c
      call reord6   ("t32 ","s1  ",123546)
      call reord7   ("s1  ","s2  ",213456)
      call add      ("t32 ",1.d0,"s2  ",1.d0,"t32 ")
      call dveps6   ("t32 ")
c
c  diag e1
      call fold6    ("t32 ","s1  ",4,1)
      call reord6   ("s1  ","s2  ",124356)
      call reorder  ("ppvh","s1  ",3412,4)
      call reorder  ("ppvh","s3  ",3421,4)
      call add      ("s1  ",2.d0,"s3  ",-1.d0,"s1  ")
      call mult6    ("s2  ","s1  ","dum2",2)
c  diag e2
      call fold6    ("t31 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord6   ("s2  ","s1  ",124356)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",125436)
      call reorder  ("vhhp","s2  ",1324,4)
      call add      ("vhhp",-2.d0,"s2  ",1.d0,"s2  ")
      call mult6    ("s2  ","s1  ","s3  ",2)
      call reorder  ("s3  ","s1  ",1243,4)
      call add      ("s1  ",1.d0,"dum2",1.d0,"dum2")
c  diag e3
      call fold6    ("t32 ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",126345)
      call reorder  ("ppvh","s1  ",3412,4)
      call mult6    ("s2  ","s1  ","s3  ",2)
      call reorder  ("s3  ","s1  ",1243,4)
      call add      ("s1  ",-1.d0,"dum2",1.d0,"dum2")
c  diag e4
      call fold6    ("t31 ","s1  ",5,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",125346)
      call reord7   ("s1  ","s2  ",132456)
      call reord6   ("s2  ","s1  ",126345)
      call mult6    ("vhhp","s1  ","s2  ",2)
      call reorder  ("s2  ","s1  ",1243,4)
      call add      ("s1  ",1.d0,"dum2",1.d0,"dum2")
CUK   V1 diags
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord6   ("s2  ","s3  ",124536)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124635)
      call add      ("s3  ",1.d0,"s1  ",-1.d0,"s3  ")
      call mult6    ("s3  ","hp  ","s1  ",0)
      call add      ("s1  ",1.d0,"dum2",1.d0,"dum2")
cuk
      call sqr      ("dum2","s1  ")
      call add      ("dum2",1.d0,"tuf ",1.d0,"tuf ")
      call prtime   ("diags2  ",0,0)
      return
      end
C  *******************************************************************
      subroutine diags3
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
c  diags 1,2,f1
      call mult6    ("tu  ","vppp","s3  ",1)
      call mult6    ("vvph","tv3 ","s1  ",1)
      call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
      call reorder  ("tu  ","s1  ",1432,4)
      call mult6    ("vvvv","s1  ","s2  ",1)
      call fold6    ("s3  ","s2  ",idum,2)
      call reord6   ("s3  ","t33 ",124356)
c
      call sqr3     ("t33 ")
      call zero     ("t33 ")
      call dveps6   ("t33 ")
c
c  diag e1
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reord6   ("s2  ","s1  ",124536)
      call reorder  ("vpvh","s2  ",1342,4)
      call reorder  ("vphv","s3  ",1432,4)
      call add      ("s2  ",2.d0,"s3  ",-1.d0,"s2  ")
      call mult6    ("s1  ","s2  ","s3  ",3)
      call reord6   ("s3  ","dum3",125346)
c  diag e2
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",6,1)
      call reord6   ("s2  ","s1  ",124635)
      call reorder  ("vpvh","s2  ",1342,4)
      call mult6    ("s1  ","s2  ","s3  ",3)
      call reord6   ("s3  ","s2  ",125346)
      call add      ("s2  ",-1.d0,"dum3",0.5d0,"dum3")
c  diag e3
      call reorder  ("vphv","s2  ",1432,4)
      call mult6    ("s1  ","s2  ","s3  ",3)
      call reord6   ("s3  ","s1  ",125364)
      call add      ("s1  ",-1.d0,"dum3",1.d0,"dum3")
c  diag e4
      call fold6    ("t31 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call fold6    ("s2  ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",126345)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",126543)
      call mult6    ("vvhh","s2  ","s1  ",3)
      call add      ("s1  ",0.5d0,"dum3",1.d0,"dum3")
c  diag e5
      call fold6    ("t33 ","s1  ",4,1)
      call reorder  ("ppvv","s2  ",3412,4)
      call mult6    ("s1  ","s2  ","s3  ",3)
      call add      ("s3  ",0.5d0,"dum3",1.d0,"dum3")
c
CUK V1 diags
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call fold6    ("s2  ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",124563)
      call mult6    ("s2  ","vh  ","s1  ",5)
      call reord6   ("s1  ","s2  ",126345)
      call add      ("s2  ",-0.5d0,"dum3",1.d0,"dum3")
c
      call fold6    ("t33 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call mult6    ("s2  ","vp  ","s1  ",5)
      call add      ("s1  ",0.5d0,"dum3",1.d0,"dum3")
cuk
      call sqr3     ("dum3")
      call add      ("dum3",1.d0,"twf ",1.d0,"twf ")
      call prtime   ("diags3  ",0,0)
      return
      end 
C  *******************************************************************
      subroutine diags4
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
c  diag e1
      call fold6    ("t33 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call reorder  ("pvvv","s1  ",2341,4)
      call mult8    ("s2  ","s1  ","s3  ",1)
      call reordx   ("s3  ","dum4",12364578)
c  diag e2
      call fold6    ("t32 ","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call fold6    ("s2  ","s1  ",6,1)
      call reord6   ("s1  ","s2  ",124563)
      call mult8    ("s2  ","vvvh","s1  ",1)
      call reordx   ("s1  ","s2  ",12673485)
      call add      ("s2  ",-0.5d0,"dum4",0.5d0,"dum4")
c
      call sqr4     ("dum4")
      call add      ("dum4",1.d0,"txf ",1.d0,"txf ")
      call prtime   ("diags4  ",0,0)
      return
      end 
C  *******************************************************************
