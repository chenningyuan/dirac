!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine chw
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
#endif


      character*256 filenm,filenm51,filenm53,filenm55,filenm60,filenm6,
     1              filenm20,filenm31,filenm32,filenm33
      real*4 cputime

c     iopt(1) - print.
c     iopt(2) - no. of iterations.
c     iopt(3) - convergence = 10**(-iopt(3))
c     iopt(4) - do rle every iopt(4) iterations
c     iopt(5) - print only elements .gt. 10**(-iopt(5))
c     iopt(6) - start avrg.
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
      filenm20(1:1)=' '
      call getenv('tape20',filenm20)
      if (filenm20(1:1).eq.' ') filenm20='tape20'
      open(20,file=filenm20,form='UNFORMATTED')
      filenm31(1:1)=' '
      call getenv('tape31',filenm31)
      if (filenm31(1:1).eq.' ') filenm31='tape31'
      open(31,file=filenm31,form='UNFORMATTED')
      filenm32(1:1)=' '
      call getenv('tape32',filenm32)
      if (filenm32(1:1).eq.' ') filenm32='tape32'
      open(32,file=filenm32,form='UNFORMATTED')
      filenm33(1:1)=' '
      call getenv('tape33',filenm33)
      if (filenm33(1:1).eq.' ') filenm33='tape33'
      open(33,file=filenm33,form='UNFORMATTED')
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')
      write (*,*) '   $$$$$$$$$   chw   $$$$$$$$$  '
      csdstim = cputime(real(0.))
      csdltim = csdstim
      call rdin
      call rdsym
      convrg = 10.**(-iopt(3))
      write (6,'('' max # of iter='',i3,''. convrg='',e10.2)') 
     1 iopt(2),convrg
      call prtime ("rdin    ",0,0)
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      call reads   (51,0,"s1  ",-1)
      call reorder ("s1  ","tc  ",3412,4)
      call reads   (51,1,"s1  ",-1)
      call reorder ("s1  ","t1c ",1243,4)
      call rdfile  (53,0,"thv ",-1)
      call rdfile  (53,1,"t1hv",-1)
      call rdfile  (53,2,"thvf",-1)
      call rdfile  (55,0,"thu ",-1)
      call rdfile  (55,1,"thuf",-1)
      if (iopt(17).ne.0) then
        !call rdfile (53,3,"t31 ")
        call rdfile (53,3,"t31 ",-1)
        !call rdfile (55,2,"t32 ")
        call rdfile (55,2,"t32 ",-1)
        if (iopt(17).eq.3) then
          call rdbox_chw
        end if
      end if
      call t1sqr    ("t1c ","t1c ","tc  ","tcm ")
      call t1sqr    ("t1hv","t1c ","thv ","tvm ")
      call t1sqr    ("t1hv","t1hv","thu ","tum ")
      call ordv_chw
      call ordta0_chw
      if (iopt(17).ne.0) then
        call calvw0_chw
        call move     ("v30 ","t33 ")
        call sqr3     ("t33 ")
        call zero     ("t33 ")
        if (iopt(17).ne.3) then
          call dveps6   ("t33 ")
          call cmult    ("t33 ",-1.d0,"t33 ")
        else if (iopt(17).eq.3) then
          iter=0
          if (iopt(4).ne.0) then
            call dbdiag   ("dum1","dum2","t33 ","dum4",0010,1)
            call dveps6   ("t33 ")
            call cmult    ("t33 ",-1.d0,"t33 ")
            call wrtrle   ("t33 ",3,1)   
          else if (iopt(4).eq.0) then
            call dveps6   ("t33 ")
            call cmult    ("t33 ",-1.d0,"t33 ")
          end if
          call prt6 ("t33 ","t33 ",2)
          call fmax6 ("t33 ",tmax,idmax,indmax)
          write (6,2400) tmax,idec(idmax),idec(indmax)
          call prt6 ("t33 ","t33 ",2)
          call fmax6 ("t33 ",tmax,idmax,indmax)
          write (6,2400) tmax,idec(idmax),idec(indmax)
          if (iopt(4).ne.0) then
            call dbdiag   ("dum1","dum2","v30 ","dum4",0010,-1)
            call wrtrle   ("t33 ",3,1)   
          end if
          iquit=0
          do while (iquit.eq.0)
            iter=iter+1
            call prtime("strt itr",iter,0)
            call caltw3_chw 
            if (iopt(6).ne.0.and.iter.ge.mod(iopt(6),100)) then
              fnew=max(dfloat(iopt(6)/100),1.d0)
              if (fnew.le.6.5d0) then
                fnew=1.d0/fnew
              else
                fnew=fnew-5.d0
              endif
              call add ("t3nw",fnew/(1.d0+fnew),"t33 ",
     1        1.d0/(1.d0+fnew),"t3nw")
              write (6,'('' t averaged with fnew '',f7.3)') fnew
            end if
            call fmax6   ("t3nw",tmax,idmax,indmax)
            call dfmax3 ("t3nw","t33 ","s2  ",df3,iddf3,inddf3)
            write (6,*) ' '
            write (6,'('' iteration # '',i3,'' diffmx='',
     1      e10.3,2x,i12,2x,i12)') iter,df3,idec(iddf3),idec(inddf3)
            write(6,2400) tmax,idec(idmax),idec(indmax)
2400        format(' max t3=',f15.9,2x,i6,2x,i10)
            if (abs(tmax).gt.1.d10.or.abs(tmax).lt.1.d-10) then
              write (*,*) ' abort. tmax out of bounds'
              call abext('tmax    ')      
            endif
            call move ("t3nw","t33 ")   
            call prt  ("t33 ","new t33 ",2)
            if (abs(df3).lt.convrg) then      
              iquit=1
              iparm=1
            end if
            if (iter.eq.iopt(2)) then
            iquit=1
            iparm=2
            end if
            if (iquit.eq.0) then
              if (iopt(4).ne.0) then
                itr1=mod(iter-1,iopt(4))+1
                call wrtrle ("t33 ",3,itr1+1)
                if (itr1.eq.iopt(4)) then
                  call rle (itr1,0010,"dum1","dum2","t33 ","dum4",ising)
                  if (ising.eq.1) then
                    iopt(4)=0    
                    write (6,*) 'rle singular, turned off'
                  end if
                  call wrtrle ("t33 ",3,1)
                end if
              end if
            end if
          end do
          if (iparm.eq.1) then
            write (6,*) ' chw converged after iter ',iter
          else if (iparm.eq.2) then
            write (6,*) ' iteration limit'
            call abext('limit   ')     
          end if
          call prt6 ("t33 ","last t33",0)
        end if
      end if
      call endiag_chw 
      call eneru2 (egr0,1,2,max(1,iopt(1)))
      call eneru2 (egr0,2,2,max(1,iopt(1)))
      call eneru2 (egr0,3,2,max(1,iopt(1)))
      call prt ("t1hv","t1hv",0)
      call prt ("thvf","thvf",0)
      call prt ("thu ","thu ",0)
      call prt ("thuf","thuf",0)
      if (iopt(17).ne.0) then
        call prt6 ("t33 ","t33 ",0)
      end if
      call prt6 ("thwf","thwf",0)
      call wrtfile (60,0,"thwf")
      if (iopt(17).ne.0) then
        call wrtfile (60,1,"t33 ")
      end if
      write (6,*) nstrt(1,nflcur),' words used in igen,',
     1nstrt(2,nflcur),' words used in gen'
      call prtime ("end ",0,0)
!      stop
      end
C ********************************************************************
      subroutine ordv_chw
!---------------------------------------------------------------------
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
!     common /cint/nint,idum1,vint(ndint),iint(ndint)
#include "cint.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)
#include "fock.h"
#endif

      call readv(nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call prtime   ("readv   ",nint,0)
      call seltyp   (3322,0,"pphh",4)
      call seltyp   (3222,0,"phhh",4)
      call seltyp   (2222,0,"hhhh",4)
      call seltyp   (2223,0010,"hhgp",4)
      call seltyp   (2222,0011,"hhgg",4)
      call seltyp   (2222,0111,"hggg",4)
      call seltyp   (2332,1000,"gpph",4)
      call seltyp   (2232,1001,"ghpg",4)
      call seltyp   (2222,1000,"ghhh",4)
      call seltyp   (2223,1010,"ghgp",4)
      call seltyp   (2233,1100,"ggpp",4)
      call seltyp   (2223,1100,"gghp",4)
      call seltyp   (2223,1110,"gggp",4)
      call seltyp   (2222,1111,"gggg",4)
      call selt1    (32,00,"ph  ",4)
      call selt1    (23,00,"hp  ",4)
      call reorder  ("phhh","vd14",2143,4)
      call prtime   ("ordv",0,0)
      return
      end
C *********************************************************************
      subroutine ordta0_chw
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
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      call reorder    ("thv ","tv1 ",1324,4)
      call reorder    ("thv ","tv3 ",1432,4)
      call reorder    ("thu ","tu3 ",1432,4)
      call reorder    ("tcm ","tcm4",3412,4)
      call reorder    ("tum ","tum4",3412,4)
      call reorder    ("tvm ","tvm2",1423,4)
      call reorder    ("t1c ","t1cr",1243,4)
      call reorder    ("t1hv","t1vr",1243,4)
      call fold4      ("t1c ","t1c1",2,2)
      return
      end
C  ********************************************************************
      subroutine calvw0_chw  
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
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      if (iopt(17).lt.3) then
c  dg3n1,2,f1
        call mult6   ("thu ","ghhh","s3  ",1)
        call mult6   ("gghp","tv3 ","s1  ",1)
        call add     ("s1  ",1.d0,"s3  ",-1.d0,"s3  ")    
        call mult6   ("gggg","tu3 ","s2  ",1)
        call cmult   ("s2  ",-1.d0,"s1  ")
        call fold6   ("s3  ","s1  ",idum,2)
        call reord6  ("s3  ","v30 ",124356)
      else if (iopt(17).eq.3) then
c  diag 1
        call mult6   ("thu ","bx12","s1  ",1)
        call reord6  ("s1  ","v30 ",124365)
c  diag 2
        call mult6   ("bx24","tv3 ","s1  ",1)
        call reord6  ("s1  ","s2  ",124356)
        call add     ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c  diag 3
        call reord6  ("t32 ","s1  ",124536)
        call mult6   ("s1  ","bx15","s3  ",3)
c  diag 4
        call reord6  ("t32 ","s1  ",124635)
        call mult6   ("s1  ","bx16","s2  ",3)
        call add     ("s2  ",1.d0,"s3  ",0.5d0,"s1  ")
        call reord6  ("s1  ","s2  ",125346)
        call add     ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c  diag 5
        call reord6  ("t32 ","s1  ",124635)
        call mult6   ("s1  ","bx17","s2  ",3)
        call reord6  ("s2  ","s1  ",125364)
        call add     ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
c  diag 6
        call reord6  ("t31 ","s1  ",126345)
        call reord7  ("s1  ","s2  ",132456)
        call reord6  ("s2  ","s1  ",126543)
        call mult6   ("bx25","s1  ","s2  ",3)
        call add     ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
c  diag 11
        call reord6  ("t32 ","s1  ",124563)
        call mult6   ("s1  ","bx19","s2  ",5)
        call reord6  ("s2  ","s1  ",126345)
        call add     ("s1  ",0.5d0,"v30 ",1.d0,"v30 ")
      end if
      call prtime  ("calvw0  ",0,0)
      return
      end
C  **********************************************************************
      subroutine caltw3_chw
!---------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use prm_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /prm/ iperm(0:6,720)
#include "prm.h"
#endif

      call move     ("v30 ","t3nw")     
c   diag 7
      call t1t2     ("t1hv","thu ","t33 ","s1  ",123456,2)
      call t1t2     ("t1hv","thu ","s1  ","s2  ",123456,3)
      call mult6    ("s2  ","hhhh","s1  ",3)
      call add      ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
c   diag 8 
      call mult6    ("t33 ","vd14","s1  ",3)
      call mult6    ("s1  ","t1cr","s2  ",5)
      call add      ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
c   diag 9 
      call mult6    ("t33 ","pphh","s1  ",3)
      call mult6    ("s1  ","tcm4","s2  ",3)
      call add      ("s2  ",0.5d0,"t3nw",1.d0,"t3nw")
c  diag 10
      call mult6    ("t33 ","bx8 ","s1  ",5)
      call add      ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
c
      call fold6    ("t3nw","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call fold6    ("s2  ","thwf",6,1)
      call sqr3     ("thwf")
c
c   diag f1
      call reord7   ("t33 ","s1  ",231456)
      call reord6   ("s1  ","s2  ",124563)
      call mult6    ("s2  ","thvf","s1  ",5)
      call reord6   ("s1  ","s2  ",126345)
      call reord7   ("s2  ","s1  ",312456)
      call add      ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
c   diag f2
      call mult6    ("thuf","tu3 ","s1  ",1)
      call reord6   ("s1  ","s2  ",124356)
      call cmult    ("s2  ",-1.d0,"s1  ")
      call fold6    ("t3nw","s1  ",idum,2)
c   diag f3
      call t1t2     ("t1hv","thu ","t33 ","s1  ",123456,1)
      call t1t2     ("t1hv","thu ","s1  ","s2  ",123456,2)
      call reord7   ("s2  ","s1  ",312456)
      call reord6   ("s1  ","s2  ",124356)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",125634)
      call mult6    ("thuf","s2  ","s1  ",3)
      call add      ("s1  ",-0.5d0,"t3nw",1.d0,"t3nw")
c   diag f4
      call mult6    ("thwf","t1vr","s1  ",5)
      call cmult    ("s1  ",-0.5d0,"s1  ")
      call fold6    ("t3nw","s1  ",idum,2)
c   diag f5
      call mult6    ("thwf","tum4","s1  ",3)
      call cmult    ("s1  ",0.5d0,"s1  ")
      call fold6    ("t3nw","s1  ",idum,2)
c   diag f6
      call cmult    ("thu ",3.d0,"s1  ")
      call t1sqr    ("t1hv","t1hv","s1  ","s2  ")
      call t1t2     ("t1hv","s2  ","t33 ","s1  ",123456,3)
      call reord6   ("s1  ","s2  ",124356)
      call reord7   ("s2  ","s1  ",312456)
      call reord6   ("s1  ","s2  ",125346)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",126345)
      call mult6    ("thwf","s2  ","s1  ",4)
      call add      ("s1  ",1.d0/6.d0,"t3nw",1.d0,"t3nw")
c
      call zero     ("t3nw")
      call sqr3     ("t3nw")
      call dveps6   ("t3nw")
      call cmult    ("t3nw",-1.d0,"t3nw")
      call prtime   ("caltw3  ",0,0)
      return
      end
C  *********************************************************************
      subroutine endiag_chw
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
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

c  GROUP I
c   diag 1
      call reorder    ("hggg","box1",2431,4)
c g1 Apr 95
      call fold4      ("thv  ","s1  ",3,4)
      call fold4      ("s1  ","s2  ",4,4)
      call reorder    ("s2  ","s1  ",1342,4)
      call mult       ("s1  ","hp  ","s2  ",311)
      call add        ("box1",1.d0,"s2  ",-.5d0,"box1")

c   diag 2
      call fold4      ("tv1 ","s1  ",2,4)
      call reorder    ("hhgp","s2  ",3142,4)
      call reorder    ("hhgp","s3  ",3241,4)
      call add        ("s2  ",2.d0,"s3  ",-1.d0,"s2  ")
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",-1.d0,"box1",-1.d0,"box1")
c   diag 3
      call fold4      ("tvm2","s1  ",2,4)
      call reorder    ("hhgp","s2  ",3142,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",1.d0,"box1",1.d0,"box1")
c   diag 4
      call reorder    ("hhgp","s2  ",3241,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","s1  ",1324,4)
      call add        ("s1  ",1.d0,"box1",1.d0,"box1")
c   diag 5
      call reorder    ("gpph","s1  ",1432,4)
      call fold4      ("tcm ","s2  ",3,4)
      call fold4      ("s2  ","s3  ",4,4)
      call reorder    ("s3  ","s2  ",3412,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","s1  ",1342,4)
      call add        ("s1  ",-1.d0,"box1",1.d0,"box1")
c   diag 6
      call fold4      ("thv ","s1  ",3,4)
      call clear6     ("s2  ",233222,100101)
      call t1t2       ("t1c1","s1  ","s2  ","s3  ",123456,3)
      call fold4      ("s1  ","s2  ",4,4)
      call t1t2       ("t1c ","s2  ","s3  ","s1  ",123456,2)
      call reord6     ("s1  ","s2  ",124356)
      call reord7     ("s2  ","s1  ",132456)
      call reord6     ("s1  ","s2  ",126543)
      call reorder    ("pphh","s1  ",3412,4)
      call reorder    ("pphh","s3  ",3421,4)
      call add        ("s1  ",-2.d0,"s3  ",1.d0,"s1  ")
      call mult6      ("s2  ","s1  ","s3  ",2)
      call add        ("s3  ",1.d0,"box1",1.d0,"box1")
c  diag 7
      call fold4      ("tc  ","s1  ",3,4)
      call fold4      ("s1  ","s2  ",4,4)
      call clear6     ("s1  ",233222,100011)
      call t1t2       ("t1hv","s2  ","s1  ","s3  ",123456,1)
      call fold4      ("thv ","s1  ",4,4)
      call t1t2       ("t1c1","s1  ","s3  ","s2  ",123456,2)
      call fold4      ("tvm ","s1  ",4,4)
      call t1t2       ("t1c1","s1  ","s2  ","s3  ",123456,3)
      call reord6     ("s3  ","s1  ",125346)
      call reord7     ("s1  ","s2  ",132456)
      call reord6     ("s2  ","s1  ",126543)
      call reorder    ("pphh","s2  ",3412,4)
      call mult6      ("s1  ","s2  ","s3  ",2)
      call add        ("s3  ",1.d0,"box1",1.d0,"box1")
c  diag 8
      call reorder    ("ghgp","s1  ",1324,4)
      call fold4      ("t1cr","s2  ",1,2)
      call mult       ("s1  ","s2  ","s3  ",311)
      call reorder    ("s3  ","s1  ",1243,4)
      call add        ("s1  ",-1.d0,"box1",1.d0,"box1")
c  diag 9
      call reorder    ("ghpg","s1  ",1423,4)
      call mult       ("s1  ","s2  ","s3  ",311)
      call reorder    ("s3  ","s1  ",1423,4)
      call add        ("s1  ",-1.d0,"box1",1.d0,"box1")
c  diag 10
      call reorder    ("hhgg","s1  ",3421,4)
      call mult       ("t1hv","s1  ","s2  ",131)
      call add        ("s2  ",1.d0,"box1",1.d0,"box1")
c
      if (iopt(17).lt.3) then
        call fold4      ("thu ","s1  ",3,4)
        call mult6      ("s1  ","box1","s2  ",1)
        call reord6     ("s2  ","thwf",124365)
      else if (iopt(17).eq.3) then
        call wrtfile    (33,0,"box1")
      end if
c  GROUP II
c   diag 1
      call move       ("gggp","box2")
CUK hp diagram
      call fold4      ("thu ","s1  ",3,4)
      call mult       ("s1  ","ph  ","s2  ",311)
      call add        ("box2",1.d0,"s2  ",-1.d0,"box2")
cuk
c   diag 2
      call fold4      ("tv1 ","s1  ",2,4)
      call reorder    ("gpph","s2  ",1324,4)
      call add        ("s2  ",2.d0,"gpph",-1.d0,"s2  ")
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","s1  ",1324,4)
      call add        ("s1  ",1.d0,"box2",1.d0,"box2")
c   diag 3
      call fold4      ("tvm2","s1  ",2,4)
      call reorder    ("gpph","s2  ",1324,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","s2  ",1324,4)
      call add        ("s2  ",-1.d0,"box2",1.d0,"box2")
c   diag 4
      call mult       ("gpph","s1  ","s2  ",222)
      call reorder    ("s2  ","s1  ",1342,4)
      call add        ("s1  ",-1.d0,"box2",1.d0,"box2")
c   diag 5
      call reorder    ("hhgp","s1  ",3412,4)
      call mult       ("tum ","s1  ","s2  ",222)
      call add        ("s2  ",1.d0,"box2",1.d0,"box2")
c   diag 6
      call reorder    ("pphh","s1  ",1243,4)
      call add        ("s1  ",-2.d0,"pphh",1.d0,"s1  ")
      call mult       ("t1hv","s1  ","s2  ",131)
      call fold4      ("tv1 ","s1  ",2,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","s1  ",1324,4)
      call add        ("s1  ",1.d0,"box2",1.d0,"box2")
c   diag 7
      call fold4      ("thv ","s1  ",4,4)
      call clear6     ("s2  ",223222,110001)
      call t1t2       ("t1hv","s1  ","s2  ","s3  ",123456,1)
      call fold4      ("thv ","s1  ",4,4)
      call t1t2       ("t1hv","s1  ","s3  ","s2  ",123456,2)
      call t1t2       ("t1c1","tum ","s2  ","s1  ",123456,3)
      call reord6     ("s1  ","s2  ",126345)
      call mult6      ("s2  ","pphh","s1  ",2)
      call reorder    ("s1  ","s2  ",2134,4)
      call add        ("s2  ",1.d0,"box2",1.d0,"box2")
c   diag 8
      call reorder    ("ghgp","s1  ",1342,4)
      call mult       ("s1  ","t1hv","s2  ",311)
      call reorder    ("s2  ","s1  ",1423,4)
      call add        ("s1  ",-1.d0,"box2",1.d0,"box2")
c   diag 9
      call reorder    ("ghpg","s1  ",1432,4)
      call mult       ("t1hv","s1  ","s2  ",131)
      call add        ("s2  ",-1.d0,"box2",1.d0,"box2")
c   diag 10
      call fold4      ("t1cr","s1  ",1,2)
      call mult       ("ggpp","s1  ","s2  ",311)
      call reorder    ("s2  ","s1  ",2143,4)
      call add        ("s1  ",1.d0,"box2",1.d0,"box2")
c
      if (iopt(17).lt.3) then
        call fold4      ("thv ","s1  ",3,4)
        call fold4      ("s1  ","s2  ",4,4)
        call reorder    ("s2  ","s1  ",1432,4)
        call mult6      ("box2","s1  ","s2  ",1)
        call reord6     ("s2  ","s1  ",124356)
        call add        ("s1  ",1.d0,"thwf",1.d0,"thwf")
      else if (iopt(17).eq.3) then
        call wrtfile    (33,1,"box2")
      end if
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
c  diag e1
        call fold6    ("t32 ","s1  ",4,1)
        call fold6    ("s1  ","s2  ",5,1)
        call reord6   ("s2  ","s1  ",124536)
        call reorder  ("ghgp","s2  ",1342,4)
        call reorder  ("ghpg","s3  ",1432,4)
        call add      ("s2  ",2.d0,"s3  ",-1.d0,"s2  ")
        call mult6    ("s1  ","s2  ","s3  ",3)
        call reord6   ("s3  ","s1  ",125346)
        call add      ("s1  ",0.5d0,"thwf",1.d0,"thwf")
c  diag e2
        call fold6    ("t32 ","s1  ",4,1)
        call fold6    ("s1  ","s2  ",6,1)
        call reord6   ("s2  ","s1  ",124635)
        call reorder  ("ghgp","s2  ",1342,4)
        call mult6    ("s1  ","s2  ","s3  ",3)
        call reord6   ("s3  ","s2  ",125346)
        call add      ("s2  ",-1.d0,"thwf",1.d0,"thwf")
c  diag e3
        call reorder  ("ghpg","s2  ",1432,4)
        call mult6    ("s1  ","s2  ","s3  ",3)
        call reord6   ("s3  ","s1  ",125364)
        call add      ("s1  ",-1.d0,"thwf",1.d0,"thwf")
c  diag e4
        call fold6    ("t31 ","s1  ",4,1)
        call fold6    ("s1  ","s2  ",5,1)
        call fold6    ("s2  ","s1  ",6,1)
        call reord6   ("s1  ","s2  ",126345)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","s2  ",126543)
        call mult6    ("ggpp","s2  ","s1  ",3)
        call add      ("s1  ",0.5d0,"thwf",1.d0,"thwf")
c  diag e5
        call fold6    ("t33 ","s1  ",4,1)
        call reorder  ("hhgg","s2  ",3412,4)
        call mult6    ("s1  ","s2  ","s3  ",3)
        call add      ("s3  ",0.5d0,"thwf",1.d0,"thwf")
      end if
      if (iopt(17).ne.3) then
        call sqr3     ("thwf")
      end if
      call prtime     ("endiag  ",0,0)
      return
      end
C  ********************************************************************
      subroutine rdbox_chw
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call rdfile  (20,7,"bx8 ",-1)
      call rdfile  (31,0,"bx12",-1)
      call rdfile  (31,3,"bx15",-1)
      call rdfile  (31,4,"bx16",-1)
      call rdfile  (31,5,"bx17",-1)
      call rdfile  (31,7,"bx19",-1)
      call rdfile  (32,0,"bx24",-1)
      call rdfile  (32,1,"bx25",-1)
      return
      end
C  *******************************************************************
