!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine csz
!----------------------------------------------------------------------------------------
!
!  (0,6)
!
!----------------------------------------------------------------------------------------
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

      dimension iord(12)
      character*256 filenm,filenm6,                                     &
 filenm51,filenm52,filenm54,filenm59,filenm61
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
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')
      write (6,*) '   $$$$$$$$$   csz   $$$$$$$$$  '
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
      call rdfile  (61,1,"tyf ",-1)
      if (iopt(17).eq.3) then
        call rdfile  (52,3,"t31 ",-1)
        call rdfile  (54,2,"t32 ",-1)
        call rdfile  (59,1,"t33 ",-1)
      end if
      nvalo=0
      do 20 j=1,120
        if (icv(j).eq.1) then
          nvalo=nvalo+1
          iv(nvalo)=j      
        end if
20    continue
      call ordv_csz
      call ordta_csz
      call t1sqr      ("t1c ","t1c ","tc  ","tcm ")
      call t1sqr      ("t1v ","t1c ","tv  ","tvm ")
      call t1sqr      ("t1v ","t1v ","tu  ","tum ")
      call diags1_csz
      if (iopt(17).eq.3) then
        call diags31_csz
        call diags32_csz
      end if
      call fmax    ("tzf ",tmax,idum)
      write (6,99) tmax
99    format (1x,'tmax. unsymmetrised=',e10.3)
      if (abs(tmax).gt.1.d-13) then
        call perms (6)
        call move    ("tzf ","s1  ")
        do 30 i=2,720
          do 40 j=1,6
            iord(j)=iperm(j,i)
40        continue
          call reord12 ("s1  ","s2  ",iord)
          call reord13 ("s2  ","s3  ",iord)
          call add     ("s3  ",1.d0,"tzf ",1.d0,"tzf ")
          if (i/10.d0-int(i/10.d0).lt.0.05d0) then
            call prtime ("perm",i,0)
          end if
30      continue
        call fmax    ("tzf ",tmax,idum)
        write (6,98) tmax
98      format (1x,'tmax. symmetrised=',e10.3)
      end if
      egr0=0.d0
      call eneru2  (egr0,1,3,max(1,iopt(1)))                         
      call eneru2  (egr0,2,3,max(1,iopt(1)))                         
      call eneru2  (egr0,3,3,max(1,iopt(1)))                         
      call eneru2  (egr0,4,3,max(1,iopt(1)))                         
      call eneru2  (egr0,5,3,max(1,iopt(1)))                         
      call eneru2  (egr0,6,3,max(1,iopt(1)))                         
      call wrtfile (61,2,"tzf ")
      write (6,*) nstrt(1,nflcur),'words used in igen,',                &
 nstrt(2,nflcur),' words used in gen'
      call prtime ("end ",0,0)
!      stop
      end
!  *******************************************************************
      subroutine ordv_csz
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
      call seltyp   (2233,0,"hhpp",4)
      call seltyp   (3332,0010,"ppvh",4)
      call seltyp   (3333,0011,"ppvv",4)
      call seltyp   (3322,1000,"vphh",4)
      call seltyp   (3323,1001,"vphv",4)
      call seltyp   (3332,1010,"vpvh",4)
      call seltyp   (3322,1100,"vvhh",4)
      call prtime   ("ordv",0,0)
      return
      end
!  *******************************************************************
      subroutine ordta_csz
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call fold4      ("t1c ","s1  ",2,2)
      call reorder    ("s1  ","te1 ",1243,4)
      call fold4      ("tu  ","te2 ",3,4)
      call fold4      ("tv  ","s1  ",3,4)
      call fold4      ("s1  ","s2  ",4,4)
      call reorder    ("s2  ","te3 ",1342,4)
      call prtime     ("ordta   ",0,0)
      return
      end
!  ********************************************************************
      subroutine diags1_csz
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension iord(12)
      data iord / 1,7,8,3,2,9,4,5,10,6,12,11 /
!  diag 1       
      call reorder      ("hhpp","s1  ",3412,4)
      call mult6        ("te3 ","s1  ","s2  ",1)
      call diage6       ("s2  ","te2 ","te2 ","te3 ","tzf ",iord)
      call cmult        ("tzf ",0.5d0,"tzf ")
      call prtime       ("diags1   ",0,0)
      return
      end
!  ********************************************************************
      subroutine diags31_csz
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension iord(12)
      data iord / 1,2,3,7,8,9,10,6,5,4,12,11 /
!   diag 1
      call reorder    ("vpvh","s1  ",1324,4)
      call cmult      ("s1  ",-1.d0,"box1")
!   diag 5
      call reorder    ("ppvh","s1  ",3241,4)
      call mult       ("t1v ","s1  ","s2  ",131)
      call add        ("s2  ",-1.d0,"box1",1.d0,"box1")
!   diag 11
      call reorder    ("vphh","s1  ",1243,4)
      call mult       ("s1  ","te1 ","s2  ",311)
      call reorder    ("s2  ","s1  ",1423,4)
      call add        ("s1  ",1.d0,"box1",1.d0,"box1")
!   diag 17
      call reorder    ("hhpp","s1  ",3124,4)
      call reorder    ("hhpp","s2  ",3214,4)
      call add        ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call fold4      ("tv  ","s2  ",3,4)
      call reorder    ("s2  ","s3  ",1324,4)
      call mult       ("s3  ","s1  ","s2  ",222)
      call add        ("s2  ",1.d0,"box1",1.d0,"box1")
!   diag 20
      call fold4      ("tvm ","s1  ",4,4)   
      call reorder    ("s1  ","s2  ",1423,4)
      call reorder    ("hhpp","s1  ",3124,4)
      call mult       ("s2  ","s1  ","s3  ",222)
      call add        ("s3  ",1.d0,"box1",1.d0,"box1")
!
      call fold6      ("t32 ","s1  ",4,1)
      call fold6      ("s1  ","s2  ",5,1)
      call fold6      ("s2  ","s1  ",6,1)
      call reord6     ("s1  ","s2  ",124563)
      call mult8      ("s2  ","box1","s1  ",1)
      call reordx     ("s1  ","s2  ",12634758)      
      call cmult      ("s2  ",0.5d0,"boxA")
!
!   diag 2
      call reorder    ("vphv","s1  ",1423,4)
      call cmult      ("s1  ",-1.d0,"box2")
!   diag 10
      call reorder    ("ppvh","s1  ",3142,4)
      call mult       ("t1v ","s1  ","s2  ",131)
      call add        ("s2  ",-1.d0,"box2",1.d0,"box2")
!   diag 16
      call mult       ("vphh","te1 ","s1  ",311)
      call reorder    ("s1  ","s2  ",1423,4)
      call add        ("s2  ",1.d0,"box2",1.d0,"box2")
!   diag 23
      call fold4      ("tvm ","s1  ",4,4)
      call reorder    ("s1  ","s2  ",1423,4)
      call reorder    ("hhpp","s1  ",3214,4)
      call mult       ("s2  ","s1  ","s3  ",222)
      call add        ("s3  ",1.d0,"box2",1.d0,"box2")
!
      call fold6      ("t32 ","s1  ",4,1)
      call fold6      ("s1  ","s2  ",5,1)
      call fold6      ("s2  ","s1  ",6,1)
      call reord6     ("s1  ","s2  ",124563)
      call mult8      ("s2  ","box2","s1  ",1)
      call reordx     ("s1  ","s2  ",12634578)       
      call add        ("s1  ",0.5d0,"boxA",1.d0,"boxA")
!
!   diag 3
      call reorder    ("ppvv","s1  ",3421,4)
      call cmult      ("s1  ",0.5d0,"box3")
!   diag 7
      call mult       ("ppvh","te1 ","s1  ",311)
      call reorder    ("s1  ","s2  ",3421,4)
      call add        ("s2  ",-1.d0,"box3",1.d0,"box3")
!   diag 29
      call fold4      ("tcm ","s1  ",3,4)
      call fold4      ("s1  ","s2  ",4,4)
      call reorder    ("s2  ","s1  ",3412,4)
      call reorder    ("hhpp","s2  ",3421,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",0.5d0,"box3",1.d0,"box3")
!
      call fold6      ("t33 ","s1  ",4,1)
      call fold6      ("s1  ","s2  ",5,1)
      call mult8      ("s2  ","box3","s1  ",1)
      call add        ("s1  ",0.5d0,"boxA",1.d0,"boxA")
!
!   diag 6
      call reorder    ("ppvh","s1  ",3241,4)
      call mult6      ("te2 ","s1  ","s2  ",1)
      call cmult      ("s2  ",-1.d0,"box4")
!   diag 14
      call reorder    ("vphh","s1  ",1243,4)
      call mult6      ("te3 ","s1  ","s2  ",1)
      call reord6     ("s2  ","s1  ",124356)
      call reord7     ("s1  ","s2  ",132456)
      call add        ("s2  ",1.d0,"box4",1.d0,"box4")
!   diag 18
      call fold6      ("t32 ","s1  ",4,1)
      call fold6      ("s1  ","s2  ",5,1)
      call reord6     ("s2  ","s1  ",124536)
      call reorder    ("hhpp","s2  ",3124,4)
      call reorder    ("hhpp","s3  ",3214,4)
      call add        ("s2  ",-2.d0,"s3  ",1.d0,"s2  ")
      call mult6      ("s1  ","s2  ","s3  ",3)
      call add        ("s3  ",0.5d0,"box4",1.d0,"box4")
!   diag 21
      call fold6      ("t32 ","s1  ",4,1)
      call fold6      ("s1  ","s2  ",6,1)
      call reord6     ("s2  ","s1  ",124635)
      call reorder    ("hhpp","s2  ",3124,4)
      call mult6      ("s1  ","s2  ","s3  ",3)
      call add        ("s3  ",1.d0,"box4",1.d0,"box4")
!
      call mult8      ("box4","te3 ","s1  ",1)
      call reordx     ("s1  ","s2  ",12634785)       
      call add        ("s2  ",1.d0,"boxA",1.d0,"boxA")
!
!   diag 9
      call reorder    ("ppvh","s1  ",3142,4)
      call mult6      ("te2 ","s1  ","s2  ",1)
      call cmult      ("s2  ",-1.d0,"box5")
!   diag 24
      call t1t2       ("t1v ","tv  ","t32 ","s1  ",123456,2)
      call fold6      ("s1  ","s2  ",4,1)
      call fold6      ("s2  ","s1  ",6,1)
      call reord6     ("s1  ","s2  ",124635)
      call reorder    ("hhpp","s1  ",3214,4)
      call mult6      ("s2  ","s1  ","s3  ",3)
      call add        ("s3  ",1.d0,"box5",1.d0,"box5")
!
      call mult8      ("box5","te3 ","s1  ",1)
      call reordx     ("s1  ","s2  ",12638745)        
      call add        ("s2  ",1.d0,"boxA",1.d0,"boxA")
!
!   diag 28
      call t1t2       ("t1c ","tv  ","t31 ","s1  ",123456,2)
      call t1t2       ("t1c ","tv  ","s1  ","s2  ",123456,3)
      call fold6      ("s2  ","s1  ",4,1)
      call fold6      ("s1  ","s2  ",5,1)
      call fold6      ("s2  ","s1  ",6,1)
      call reord6     ("s1  ","s2  ",124356)
      call reord7     ("s2  ","s1  ",132456)
      call reord6     ("s1  ","s2  ",125643)
      call reorder    ("hhpp","s1  ",3412,4)
      call mult6      ("s2  ","s1  ","s3  ",3)
      call mult8      ("s3  ","te2 ","s1  ",1)
      call reordx     ("s1  ","s2  ",17623845)      
      call add        ("s2  ",1.d0,"boxA",1.d0,"boxA")
!
      call fold6      ("t33 ","s1  ",4,1)
      call fold6      ("s1  ","s2  ",5,1)
      call dge6n2     ("boxA","s2  ","s1  ",iord)
      call add        ("s1  ",0.5d0,"tzf ",1.d0,"tzf ")
      call prtime     ("diags31 ",0,0)
      return
      end
!  *******************************************************************
      subroutine diags32_csz
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension iord(12)
      data iord / 1,2,3,4,7,8,9,6,5,12,11,10 /
!   diag 4
      call cmult      ("vvhh",0.5d0,"s3  ")
!   diag 13
      call reorder    ("vphh","s1  ",1432,4)
      call mult       ("t1v ","s1  ","s2  ",131)
      call add        ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!   diag 27
      call mult       ("tum ","hhpp","s1  ",222)
      call add        ("s1  ",0.5d0,"s3  ",1.d0,"s3  ")
! 
      call fold6      ("t32 ","s1  ",4,1)
      call fold6      ("s1  ","s2  ",5,1)
      call fold6      ("s2  ","s1  ",6,1)
      call reord6     ("s1  ","s2  ",124563)
      call mult8      ("s2  ","s3  ","s1  ",1)
      call reordx     ("s1  ","s2  ",12763458)      
      call cmult      ("s2  ",0.5d0,"boxB")
!
!   diag 8
      call reorder    ("ppvh","s1  ",3421,4)
      call mult6      ("te2 ","s1  ","s2  ",1)
      call cmult      ("s2  ",-1.d0,"box6")
!   diag 12
      call reorder    ("vphh","s1  ",1423,4)
      call mult6      ("te3 ","s1  ","s2  ",1)
      call reord6     ("s2  ","s1  ",124356)
      call reord7     ("s1  ","s2  ",132456)
      call add        ("s2  ",1.d0,"box6",1.d0,"box6")
!   diag 19
      call fold6      ("t32 ","s1  ",4,1)
      call fold6      ("s1  ","s2  ",5,1)
      call reord6     ("s2  ","s1  ",124536)
      call reorder    ("hhpp","s2  ",1324,4)
      call reorder    ("hhpp","s3  ",1423,4)
      call add        ("s2  ",-2.d0,"s3  ",1.d0,"s2  ")
      call mult6      ("s1  ","s2  ","s3  ",3)
      call add        ("s3  ",0.5d0,"box6",1.d0,"box6")
!   diag 22
      call t1t2       ("t1v ","tv  ","t32 ","s1  ",123456,2)
      call t1t2       ("t1c ","tu  ","s1  ","s2  ",123456,3)
      call fold6      ("s2  ","s1  ",4,1)
      call fold6      ("s1  ","s2  ",6,1)
      call reord6     ("s2  ","s1  ",124635)
      call reorder    ("hhpp","s2  ",1324,4)
      call mult6      ("s1  ","s2  ","s3  ",3)
      call add        ("s3  ",1.d0,"box6",1.d0,"box6")
!
      call mult8      ("box6","te2 ","s1  ",1)
      call reordx     ("s1  ","s2  ",12673485)       
      call add        ("s2  ",1.d0,"boxB",1.d0,"boxB")
! 
!   diag 15
      call reorder    ("vphh","s1  ",1342,4)
      call mult6      ("te2 ","s1  ","s3  ",1)
!   diag 26
      call fold6      ("t33 ","s1  ",4,1)
      call mult6      ("s1  ","hhpp","s2  ",3)
      call reord6     ("s2  ","s1  ",124365)
      call add        ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!
      call mult8      ("s3  ","te3 ","s1  ",1)
      call reordx     ("s1  ","s2  ",12643875)       
      call add        ("s2  ",1.d0,"boxB",1.d0,"boxB")
!
!   diag 25
      call t1t2       ("t1v ","tv  ","t32 ","s1  ",123456,2)
      call fold6      ("s1  ","s2  ",4,1)
      call fold6      ("s2  ","s1  ",6,1)
      call reord6     ("s1  ","s2  ",124635)
      call reorder    ("hhpp","s1  ",1423,4)
      call mult6      ("s2  ","s1  ","s3  ",3)
      call mult8      ("s3  ","te2 ","s1  ",1)
      call reordx     ("s1  ","s2  ",17623485)       
      call add        ("s2  ",1.d0,"boxB",1.d0,"boxB")
!
      call fold6      ("t32 ","s1  ",4,1)
      call fold6      ("s1  ","s2  ",5,1)
      call fold6      ("s2  ","s1  ",6,1)
      call reord6     ("s1  ","s2  ",124563)
      call dge6n2     ("boxB","s2  ","s1  ",iord)
      call add        ("s1  ",0.5d0,"tzf ",1.d0,"tzf ")
      call prtime     ("diags32 ",0,0)
      return
      end
!  *******************************************************************
