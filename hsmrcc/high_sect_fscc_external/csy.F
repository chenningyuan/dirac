!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine csy
!---------------------------------------------------------------------
!
!                        (0,5) sector
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
      character*256 filenm,filenm2,filenm3,filenm6,
     1filenm51,filenm52,filenm54,filenm59,filenm61
      real*4 cputime
c     iopt(1) - print.
c     iopt(5) - print only elements .gt. 10**(-iopt(5))
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
      write (6,*) '   $$$$$$$$$   csy   $$$$$$$$$  '
      csdstim=cputime(real(0.))
      csdltim=csdstim
      call rdin
      call rdsym
      call prtime  ("rdin    ",0,0)
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      call rdfile  (51,0,"tc  ",-1)
      call rdfile  (51,1,"t1c ",-1)
      call rdfile  (52,0,"tv  ",-1)
      call rdfile  (52,1,"t1v ",-1)
      call rdfile  (52,2,"tvf ",-1)
      call rdfile  (54,0,"tu  ",-1)
      call rdfile  (54,1,"tuf ",-1)
      call rdfile  (59,0,"twf ",-1)
      call rdfile  (61,0,"txf ",-1)
      if (iopt(17).ge.3) then
        call rdfile  (52,3,"t31 ",-1)
        call rdfile  (54,2,"t32 ",-1)
        call rdfile  (59,1,"t33 ",-1)
      end if

#if defined HSFS_DYN_ALLOC
      call check_modules
#endif
      nvalo=0
      do 20 j=1,120
        if (icv(j).eq.1) then
          nvalo=nvalo+1
          iv(nvalo)=j      
        end if
20    continue

      call ordv_csy
      call ordta_csy
      call t1sqr    ("t1c ","t1c ","tc  ","tcm ")
      call t1sqr    ("t1v ","t1c ","tv  ","tvm ")
      call t1sqr    ("t1v ","t1v ","tu  ","tum ")
      call diags1_csy
      if (iopt(17).eq.3) then
        call diags31_csy
        call diags32_csy
        call diags33_csy
      end if

      call fmax    ("tyf ",tmax,idum)
      write (LUPRI,99) tmax
99    format (1x,'tmax,unsymmetrised=',e10.3)

!miro: why is this restricted for abs(tmax)>1.d-13 ?
      if (abs(tmax).gt.1.d-13) then
        call perms (5)
        call move    ("tyf ","s1  ")
        do 30 i=2,120
          do 40 j=1,5
            iord(j)=iperm(j,i)
40        continue
          call reord10 ("s1  ","s2  ",iord)
          call reord11 ("s2  ","s3  ",iord)
          call add     ("s3  ",1.d0,"tyf ",1.d0,"tyf ")
          if (i/10.d0-int(i/10.d0).lt.0.05d0) then
            call prtime ("perm",i,0)
          end if
30      continue
        call fmax    ("tyf ",tmax,idum)
        write (LUPRI,98) tmax
98      format (1x,'tmax,symmetrised=',e10.3)
      end if

      egr0=0.d0
      call eneru2  (egr0,1,3,max(1,iopt(1)))
      call eneru2  (egr0,2,3,max(1,iopt(1)))
      call eneru2  (egr0,3,3,max(1,iopt(1)))                         
      call eneru2  (egr0,4,3,max(1,iopt(1)))                         
      call eneru2  (egr0,5,3,max(1,iopt(1)))                         
      call wrtfile (61,1,"tyf ")

      write (LUPRI,*) nstrt(1,nflcur),'words used in igen,',
     1nstrt(2,nflcur),' words used in gen'

      call prtime ("end ",0,0)
!      stop
      end subroutine csy
C  *******************************************************************
      subroutine ordv_csy
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

c
      call readv    (nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call prtime   ("readv   ",nint,0)
      call seltyp   (2233,0,"hhpp",4)
      call seltyp   (3322,1000,"vphh",4)
      call seltyp   (3332,1010,"vpvh",4)
      call seltyp   (3323,1001,"vphv",4)
      call seltyp   (3333,0011,"ppvv",4)
      call seltyp   (3332,0010,"ppvh",4)
      call seltyp   (3322,1100,"vvhh",4)
      call selt1    (32,0,"ph  ",4)
      call prtime   ("ordv",0,0)
      return
      end
C ********************************************************************
      subroutine ordta_csy
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call fold4    ("tu  ","tua ",3,4)
      call fold4    ("t1c ","s1  ",2,2)
      call reorder  ("s1  ","tw1 ",1243,4)
      call fold4    ("tv  ","s1  ",3,4)
      call fold4    ("s1  ","s2  ",4,4)
      call reorder  ("s2  ","tw3 ",1342,4)
      call prtime   ("ordta   ",0,0)
      return
      end 
C ********************************************************************
      subroutine diags1_csy
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension iord(10)
      data iord / 1,2,6,7,5,4,10,3,8,9 /
c   diags 1,3
      call reorder ("hhpp","s1  ",1342,4)
      call mult    ("s1  ","tw1 ","s2  ",311)
      call reorder ("s2  ","s1  ",2143,4)
      call reorder ("ppvh","s2  ",2431,4)
      call add     ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
      call mult6   ("tua ","s1  ","s2  ",1)
      call reord6  ("s2  ","s3  ",123645)
c   diags 2,4
      call reorder ("hhpp","s1  ",3124,4)
      call mult    ("t1v ","s1  ","s2  ",131)
      call reorder ("vphh","s1  ",1243,4)
      call add     ("s2  ",1.d0,"s1  ",1.d0,"s1  ")
      call mult6   ("tw3 ","s1  ","s2  ",1)
      call reord6  ("s2  ","s1  ",124356)
      call reord7  ("s1  ","s2  ",132456)
      call add     ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
c
      call diage5  ("s3  ","tua ","tw3 ","tyf ",iord)
      call prtime  ("diags1  ",0,0)
      return
      end
C  ********************************************************************
      subroutine diags31_csy
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension iord(10)
      data iord / 1,2,3,6,7,5,4,10,9,8 /
c   diag 1
      call reorder   ("vpvh","s1  ",1324,4)
      call mult6     ("s1  ","tw3 ","s2  ",1)   
      call reord6    ("s2  ","s1  ",124563)
      call reord7    ("s1  ","s2  ",132456)
      call cmult     ("s2  ",-1.d0,"boxA")
c   diag 3
      call reorder   ("vphv","s1  ",1243,4)
      call mult6     ("tw3 ","s1  ","s2  ",1)
      call reord6    ("s2  ","s1  ",124365)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",-1.d0,"boxA",1.d0,"boxA")
c   diag 5
      call reorder   ("ppvv","s1  ",3421,4)
      call mult6     ("tua ","s1  ","s2  ",1)
      call add       ("s2  ",1.d0,"boxA",1.d0,"boxA")
c   diag 7
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call reord6    ("s2  ","s1  ",124536)
      call reorder   ("ppvh","s2  ",3142,4)
      call reorder   ("ppvh","s3  ",3241,4)
      call add       ("s2  ",2.d0,"s3  ",-1.d0,"s2  ")
      call mult6     ("s1  ","s2  ","s3  ",3)
      call add       ("s3  ",0.5d0,"boxA",1.d0,"boxA")
c   diag 8
      call t1t2      ("t1v ","tv  ","t32 ","s1  ",123456,2)
      call t1t2      ("t1c ","tu  ","s1  ","s2  ",123456,3)
      call fold6     ("s2  ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",6,1)
      call reord6    ("s2  ","s1  ",124635)
      call reorder   ("ppvh","s2  ",3142,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call add       ("s3  ",-1.d0,"boxA",1.d0,"boxA")
c   diag 9
      call reorder   ("ppvh","s2  ",3241,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125346)
      call add       ("s1  ",-1.d0,"boxA",1.d0,"boxA")
c   diag 14
      call t1t2      ("t1c ","tv  ","t31 ","s1  ",123456,2)
      call t1t2      ("t1c ","tv  ","s1  ","s2  ",123456,3)
      call fold6     ("s2  ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call fold6     ("s2  ","s1  ",6,1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125634)
      call mult6     ("s2  ","vphh","s1  ",3)
      call reord6    ("s1  ","s2  ",125346)
      call reord7    ("s2  ","s1  ",132456)
      call add       ("s1  ",1.d0,"boxA",1.d0,"boxA")
c   diag 15
      call reorder   ("hhpp","s1  ",3124,4)
      call reorder   ("hhpp","s2  ",3214,4)
      call add       ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call mult      ("s1  ","t1c ","s3  ",202)
c v1 Feb 95
      call add       ("s3  ",1.d0,"ph  ",-1.d0,"s3  ")
c      
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call fold6     ("s2  ","s1  ",6,1)
      call reord6    ("s1  ","s2  ",124563)  
      call mult6     ("s2  ","s3  ","s1  ",5)
      call add       ("s1  ",0.5d0,"boxA",1.d0,"boxA")
c   diag 18
      call reorder   ("hhpp","s1  ",3124,4)
      call reorder   ("hhpp","s2  ",3214,4)
      call add       ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call fold4     ("tv  ","s2  ",3,4)
      call reorder   ("s2  ","s3  ",1324,4)
      call mult      ("s3  ","s1  ","s2  ",222)
      call mult6     ("s2  ","tw3 ","s1  ",1)
      call reord6    ("s1  ","s2  ",124563)
      call reord7    ("s2  ","s1  ",132456)
      call add       ("s1  ",-1.d0,"boxA",1.d0,"boxA")
c   diag 20
      call reorder   ("hhpp","s1  ",3241,4)
      call reorder   ("hhpp","s2  ",3142,4)
      call add       ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call mult      ("tw1 ","s1  ","s2  ",131)
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s3  ",5,1)
      call reord6    ("s3  ","s1  ",124536)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call add       ("s3  ",0.5d0,"boxA",1.d0,"boxA")
c   diag 21
      call fold4     ("tvm ","s1  ",4,4)
      call reorder   ("s1  ","s2  ",1423,4)
      call reorder   ("hhpp","s1  ",3124,4)
      call mult      ("s2  ","s1  ","s3  ",222)
      call mult6     ("s3  ","tw3 ","s1  ",1)
      call reord6    ("s1  ","s2  ",124563)
      call reord7    ("s2  ","s1  ",132456)
      call add       ("s1  ",1.d0,"boxA",1.d0,"boxA")
c   diag 22
      call reorder   ("hhpp","s1  ",3241,4)
      call mult      ("tw1 ","s1  ","s2  ",131)
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s3  ",6,1)
      call reord6    ("s3  ","s1  ",124635)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call add       ("s3  ",1.d0,"boxA",1.d0,"boxA")
c   diag 26
      call reorder   ("hhpp","s2  ",3142,4)
      call mult      ("tw1 ","s2  ","s3  ",131)
      call mult6     ("s1  ","s3  ","s2  ",3)
      call reord6    ("s2  ","s1  ",123546)
      call add       ("s1  ",1.d0,"boxA",1.d0,"boxA")
c   diag 27
      call fold4     ("tvm ","s1  ",4,4)
      call reorder   ("s1  ","s2  ",1423,4)
      call reorder   ("hhpp","s1  ",3214,4)
      call mult      ("s2  ","s1  ","s3  ",222)
      call mult6     ("tw3 ","s3  ","s1  ",1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call add       ("s1  ",1.d0,"boxA",1.d0,"boxA")
c   diag 34
      call reorder   ("hhpp","s1  ",3214,4)
      call mult      ("t1v ","s1  ","s2  ",131)
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
c   diag 36
      call fold4     ("tcm ","s1  ",3,4)
      call fold4     ("s1  ","s2  ",4,4)
      call reorder   ("s2  ","s1  ",3421,4)
      call reorder   ("hhpp","s2  ",3412,4)
      call mult      ("s1  ","s2  ","s3  ",222)
      call mult6     ("tua ","s3  ","s1  ",1)
      call add       ("s1  ",1.d0,"boxA",1.d0,"boxA")
c
      call fold6     ("t33 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call dge5n2    ("s2  ","boxA","s1  ",iord)
      call add       ("s1  ",0.5d0,"tyf ",1.d0,"tyf ")
      call prtime    ("diags31 ",0,0)
      return
      end 
C  *******************************************************************
      subroutine diags32_csy
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension iord(10)
      data iord / 1,2,3,6,7,5,4,10,9,8 /
c   diag 2
      call reorder   ("vpvh","s1  ",1342,4)
      call mult6     ("s1  ","tua ","s2  ",1)
      call reord6    ("s2  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",125364)
      call cmult     ("s1  ",-1.d0,"boxB")
c   diag 4
      call reorder   ("vphv","s1  ",1432,4)
      call mult6     ("tua ","s1  ","s2  ",1)
      call reord6    ("s2  ","s1  ",124356)
      call add       ("s1  ",-1.d0,"boxB",1.d0,"boxB")
c   diag 6
      call reorder   ("vvhh","s1  ",1243,4)
      call mult6     ("tw3 ","s1  ","s2  ",1)
      call reord6    ("s2  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",125346)
      call add       ("s1  ",1.d0,"boxB",1.d0,"boxB")
c   diag 10
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call reord6    ("s2  ","s1  ",124536)
      call reorder   ("vphh","s2  ",1342,4)
      call reorder   ("vphh","s3  ",1432,4)
      call add       ("s2  ",-2.d0,"s3  ",1.d0,"s2  ")
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125346)
      call add       ("s1  ",0.5d0,"boxB",1.d0,"boxB")
c   diag 11
      call t1t2      ("t1v ","tv  ","t32 ","s1  ",123456,2)
      call t1t2      ("t1c ","tu  ","s1  ","s2  ",123456,3)
      call fold6     ("s2  ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",6,1)
      call reord6    ("s2  ","s1  ",124635)
      call reorder   ("vphh","s2  ",1342,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125346)
      call add       ("s1  ",1.d0,"boxB",1.d0,"boxB")
c   diag 12
      call t1t2      ("t1v ","tv  ","t32 ","s1  ",123456,1)
      call t1t2      ("t1c ","tu  ","s1  ","s2  ",123456,3)
      call fold6     ("s2  ","s1  ",5,1)
      call fold6     ("s1  ","s2  ",6,1)
      call reord6    ("s2  ","s1  ",125634)
      call reorder   ("vphh","s2  ",1432,4)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",124653)
      call add       ("s1  ",1.d0,"boxB",1.d0,"boxB")
c   diag 13
      call t1t2      ("t1v ","tu  ","t33 ","s1  ",123456,2)
      call t1t2      ("t1v ","tu  ","s1  ","s2  ",123456,3)
      call fold6     ("s2  ","s1  ",4,1)
      call reorder   ("ppvh","s2  ",3412,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call add       ("s3  ",-1.d0,"boxB",1.d0,"boxB")
c   diag 16
      call reorder   ("hhpp","s1  ",1324,4)
      call reorder   ("hhpp","s2  ",1423,4)
      call add       ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call fold4     ("tv  ","s2  ",3,4)
      call reorder   ("s2  ","s3  ",1324,4)
      call mult      ("s3  ","s1  ","s2  ",222)
      call mult6     ("s2  ","tua ","s1  ",1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125364)
      call add       ("s2  ",-1.d0,"boxB",1.d0,"boxB")
c   diag 17
      call reorder   ("hhpp","s1  ",1243,4)
      call add       ("s1  ",2.d0,"hhpp",-1.d0,"s1  ")
      call mult      ("t1v ","s1  ","s2  ",131)
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s3  ",5,1)
      call reord6    ("s3  ","s1  ",124536)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125346)
      call add       ("s1  ",-0.5d0,"boxB",1.d0,"boxB")
c   diag 23
      call fold4     ("tvm ","s1  ",4,4)
      call reorder   ("s1  ","s2  ",1423,4)
      call reorder   ("hhpp","s1  ",1324,4)
      call mult      ("s2  ","s1  ","s3  ",222)
      call mult6     ("s3  ","tua ","s1  ",1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125364)
      call add       ("s2  ",1.d0,"boxB",1.d0,"boxB")
c   diag 25
      call reorder   ("hhpp","s1  ",1243,4)
      call mult      ("t1v ","s1  ","s2  ",131)
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s3  ",6,1)
      call reord6    ("s3  ","s1  ",124635)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125346)
      call add       ("s1  ",1.d0,"boxB",1.d0,"boxB")
c   diag 29
      call fold4     ("tvm ","s1  ",4,4)
      call reorder   ("s1  ","s2  ",1423,4)
      call reorder   ("hhpp","s1  ",1423,4)
      call mult      ("s2  ","s1  ","s3  ",222)
      call mult6     ("tua ","s3  ","s1  ",1)
      call reord6    ("s1  ","s2  ",124356)
      call add       ("s2  ",1.d0,"boxB",1.d0,"boxB")
c   diag 30
      call mult      ("t1v ","hhpp","s1  ",131)
      call fold6     ("t32 ","s2  ",5,1)
      call fold6     ("s2  ","s3  ",6,1)
      call reord6    ("s3  ","s2  ",125634)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord6    ("s3  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",124653)
      call add       ("s1  ",1.d0,"boxB",1.d0,"boxB")
c   diag 31
      call reorder   ("hhpp","s1  ",1342,4)
      call mult      ("tw1 ","s1  ","s2  ",131)
      call fold6     ("t33 ","s1  ",5,1)
      call reord6    ("s1  ","s3  ",123564)
      call mult6     ("s3  ","s2  ","s1  ",3)
      call reord6    ("s1  ","s2  ",123546)
      call add       ("s2  ",1.d0,"boxB",1.d0,"boxB")
c   diag 33
      call reorder   ("hhpp","s1  ",1243,4)
      call mult      ("tum ","s1  ","s2  ",222)
      call mult6     ("tw3 ","s2  ","s1  ",1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125346)
      call add       ("s2  ",1.d0,"boxB",1.d0,"boxB")
c
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call fold6     ("s2  ","s1  ",6,1)
      call reord6    ("s1  ","s2  ",124563)
      call dge5n2    ("boxB","s2  ","s1  ",iord)
      call add       ("s1  ",0.5d0,"tyf ",1.d0,"tyf ")
      call prtime    ("diags32 ",0,0)
      return
      end
C  ********************************************************************
      subroutine diags33_csy
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension iord1(10),iord2(10),iord3(10),iord4(10)
      data iord1 / 1,2,6,7,3,8,9,5,4,10 /
      data iord2 / 1,3,8,9,2,7,6,5,4,10 /
      data iord3 / 1,3,4,8,2,7,6,5,10,9 /
      data iord4 / 1,6,7,9,3,2,8,5,4,10 /
c   diag 19
      call reorder   ("hhpp","s1  ",1324,4)
      call reorder   ("hhpp","s2  ",1423,4)
      call add       ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call fold6     ("t32 ","s2  ",4,1)
      call fold6     ("s2  ","s3  ",5,1)
      call reord6    ("s3  ","s2  ",124536)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call diage5    ("s3  ","tw3 ","tua ","s1  ",iord1)
      call add       ("s1  ",0.5d0,"tyf ",1.d0,"tyf ")
c   diag 24
      call fold6     ("t32 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",6,1)
      call reord6    ("s2  ","s1  ",124635)
      call reorder   ("hhpp","s2  ",1324,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call diage5    ("s3  ","tw3 ","tua ","s1  ",iord1)
      call add       ("s1  ",1.d0,"tyf ",1.d0,"tyf ")
c   diag 28
      call fold6     ("t32 ","s1  ",5,1)
      call fold6     ("s1  ","s2  ",6,1)
      call reord6    ("s2  ","s1  ",125634)
      call reorder   ("hhpp","s2  ",1423,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call diage5    ("s3  ","tw3 ","tua ","s1  ",iord2)
      call add       ("s1  ",1.d0,"tyf ",1.d0,"tyf ")
c   diag 32
      call fold6     ("t33 ","s1  ",5,1)
      call reord6    ("s1  ","s2  ",123546)
      call mult6     ("s2  ","hhpp","s1  ",3)
      call diage5    ("s1  ","tw3 ","tw3 ","s2  ",iord3)
      call add       ("s2  ",0.5d0,"tyf ",1.d0,"tyf ")
c   diag 35
      call fold6     ("t31 ","s1  ",4,1)
      call fold6     ("s1  ","s2  ",5,1)
      call fold6     ("s2  ","s1  ",6,1)
      call reord6    ("s1  ","s2  ",124356)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",125634)
      call reorder   ("hhpp","s1  ",3412,4)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call diage5    ("s3  ","tua ","tua ","s1  ",iord4)
      call add       ("s1  ",0.5d0,"tyf ",1.d0,"tyf ")
      call prtime    ("diags33 ",0,0)
      return
      end
C  ********************************************************************
