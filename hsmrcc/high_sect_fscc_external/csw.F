      subroutine csw                     
!---------------------------------------------------------------------
!
! (0,3) sector...
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

      character*256 filenm,filenm51,filenm52,filenm54,filenm59,filenm6,
     1              filenm20,filenm21,filenm22,filenm23
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
      filenm20(1:1)=' '
      call getenv('tape20',filenm20)
      if (filenm20(1:1).eq.' ') filenm20='tape20'
      open(20,file=filenm20,form='UNFORMATTED')
      filenm21(1:1)=' '
      call getenv('tape21',filenm21)
      if (filenm21(1:1).eq.' ') filenm21='tape21'
      open(21,file=filenm21,form='UNFORMATTED')
      filenm22(1:1)=' '
      call getenv('tape22',filenm22)
      if (filenm22(1:1).eq.' ') filenm22='tape22'
      open(22,file=filenm22,form='UNFORMATTED')
      filenm23(1:1)=' '
      call getenv('tape23',filenm23)
      if (filenm23(1:1).eq.' ') filenm23='tape23'
      open(23,file=filenm23,form='UNFORMATTED')
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')

      write (LUPRI,*) '   $$$$$$$$$   csw - (0,3) sector  $$$$$$$$$  '

      write(LUPRI,*) 'iopt(17)=',iopt(17)

#if defined HSFS_DYN_ALLOC
      call check_modules
#endif
      csdstim = cputime(real(0.))
      csdltim = csdstim
      call rdin
      call rdsym
      convrg = 10.**(-iopt(3))
      write (LUPRI,'('' max # of iter='',i3,''. convrg='',e10.2)') 
     1 iopt(2),convrg
      call prtime ("rdin    ",0,0)
      if (iopt(10).eq.0) then
        filenm(1:1)=' '
        call getenv('tape9',filenm)
        if (filenm(1:1).eq.' ') filenm='tape9'
        open(9,file=filenm,form='UNFORMATTED')
      endif
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
      if (iopt(17).ne.0) then
        call rdfile (52,3,"t31 ",-1)
        call rdfile (54,2,"t32 ",-1)
        if (iopt(17).eq.3) then
          call rdbox_csw
        end if
      end if

#if defined HSFS_DYN_ALLOC
      call check_modules
#endif

      call t1sqr    ("t1c ","t1c ","tc  ","tcm ")
      call t1sqr    ("t1v ","t1c ","tv  ","tvm ")
      call t1sqr    ("t1v ","t1v ","tu  ","tum ")
      call ordv_csw
      call ordta0_csw

      ! miro: no convergence for CCSD/CCSDT, fix needed
      if (iopt(17).ne.0) then
        call calvw0_csw
        call move     ("v30 ","t33 ")
        call sqr3     ("t33 ")
        call zero     ("t33 ")

        if (iopt(17).ne.3) then
          call dveps6   ("t33 ")
        else if (iopt(17).eq.3) then
          iter=0
          if (iopt(4).ne.0) then
            call dbdiag   ("dum1","dum2","t33 ","dum4",0010,1)
            call dveps6   ("t33 ")
            call wrtrle   ("t33 ",3,1)   
          else if (iopt(4).eq.0) then
            call dveps6   ("t33 ")
          end if
          call prt6 ("t33 ","t33 ",2)
          call fmax6 ("t33 ",tmax,idmax,indmax)
          write (6,2400) tmax,idec(idmax),idec(indmax)
          iquit=0
          do while (iquit.eq.0)
            iter=iter+1
            call prtime("strt itr",iter,0)
            call flush (6) 
            call caltw3_csw 
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

!------------- iteration process ------------------------------------
            write (6,*) ' '
            write (6,'('' iteration # '',i3,'' diffmx='',
     1      e10.3,2x,i12,2x,i12)') iter,df3,idec(iddf3),idec(inddf3)
            write(6,2400) tmax,idec(idmax),idec(indmax)
2400        format(' max t3=',f15.9,2x,i6,2x,i10)


            if (abs(tmax).gt.1.d10.or.abs(tmax).lt.1.d-10) then
              write (*,*) ' abort. tmax out of bounds'
              call abext('tmax    ')      
            endif

            call wrtfile (59,0,"twf ")
            if (iopt(17).gt.0) then
              call prt6 ("t33 ","t33 ",2)
              call wrtfile (59,1,"t33 ")
            end if
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

          ! --- iterations end
          if (iparm.eq.1) then
            write (6,*) ' csw converged after iter ',iter
          else if (iparm.eq.2) then
            write (6,*) ' iteration limit'
            call abext('limit   ')     
          else
            write(LUPRI,*) 'iparm=',iparm
            write(LUPRI,*) 'dead branch for iterations...iter=',iter
          end if
          call prt6 ("t33 ","last t33",0)
        end if
      else
        !
        write(LUPRI,*) '=== no iterations ==='
        write(LUPRI,*) ' iopt(17)=',iopt(17)
      end if

#if defined HSFS_DYN_ALLOC
      call check_modules
#endif

      call endiag_csw 
      call eneru2 (egr0,1,3,max(1,iopt(1)))
      call eneru2 (egr0,2,3,max(1,iopt(1)))
      call eneru2 (egr0,3,3,max(1,iopt(1)))
      call prt6 ("twf ","twf ",0)
      call wrtfile (59,0,"twf ")
      if (iopt(17).gt.0) then
        call prt6 ("t33 ","t33 ",0)
        call wrtfile (59,1,"t33 ")
      end if

      write (LUPRI,*) nstrt(1,nflcur),' words used in igen,',
     1nstrt(2,nflcur),' words used in gen'
      call prtime ("end ",0,0)
!      stop
      end
C ********************************************************************
      subroutine ordv_csw
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
      call seltyp   (2233,0,"hhpp",4)
      call seltyp   (2333,0,"hppp",4)
      call seltyp4  (3333,0,"pppp",4)
      call seltyp   (3332,0010,"ppvh",4)
      call seltyp   (3333,0011,"ppvv",4)
      call seltyp   (3333,0111,"pvvv",4)
      call seltyp   (3223,1000,"vhhp",4)
      call seltyp   (3323,1001,"vphv",4)
      call seltyp   (3333,1000,"vppp",4)
      call seltyp   (3332,1010,"vpvh",4)
      call seltyp   (3322,1100,"vvhh",4)
      call seltyp   (3332,1100,"vvph",4)
      call seltyp   (3332,1110,"vvvh",4)
      call seltyp   (3333,1111,"vvvv",4)
      call selt1    (23,00,"hp  ",4)
      call selt1    (32,0,"ph  ",4)
      call reorder  ("hppp","vd14",2143,4)
      call prtime   ("ordv",0,0)
      return
      end
C *********************************************************************
      subroutine ordta0_csw
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

      call reorder    ("tv  ","tv1 ",1324,4)
      call reorder    ("tv  ","tv3 ",1432,4)
      call reorder    ("tu  ","tu3 ",1432,4)
      call reorder    ("tcm ","tcm4",3412,4)
      call reorder    ("tvm ","tvm2",1423,4)
      call reorder    ("tum ","tum4",3412,4)
      call reorder    ("t1c ","t1cr",1243,4)
      call reorder    ("t1v ","t1vr",1243,4)
      call fold4      ("t1c ","t1c1",2,2)
      return
      end
C  ********************************************************************
      subroutine calvw0_csw  
!---------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use dat_mod
      use cfiles_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)
#include "cfiles.h"
#endif

      if (iopt(17).lt.3) then
c  dg3n1,2,f1
        call mult6   ("tu  ","vppp","s3  ",1)
        call mult6   ("vvph","tv3 ","s1  ",1)
        call add     ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")    
        call mult6   ("vvvv","tu3 ","s2  ",1)
        call fold6   ("s3  ","s2  ",idum,2)
        call reord6  ("s3  ","v30 ",124356)
      else if (iopt(17).eq.3) then
c  diag 1
        call mult6   ("tu  ","bx12","s1  ",1)
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
      subroutine caltw3_csw
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call move     ("v30 ","t3nw")     
c   diag 7
      call t1t2     ("t1v ","tu  ","t33 ","s1  ",123456,2)
      call t1t2     ("t1v ","tu  ","s1  ","s2  ",123456,3)
      call mult6    ("s2  ","pppp","s1  ",3)
      call reord6   ("s1  ","s2  ",123465)
      call reord7   ("s2  ","s3  ",132456)
      call add      ("s1  ",1.d0,"s3  ",1.d0,"s1  ")
      call add      ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
c   diag 8 
      call mult6    ("t33 ","vd14","s1  ",3)
      call mult6    ("s1  ","t1cr","s2  ",5)
      call add      ("s2  ",-1.d0,"t3nw",1.d0,"t3nw")
c   diag 9 
      call mult6    ("t33 ","hhpp","s1  ",3)
      call mult6    ("s1  ","tcm4","s2  ",3)
      call add      ("s2  ",0.5d0,"t3nw",1.d0,"t3nw")
c  diag 10
      call mult6    ("t33 ","bx7 ","s1  ",5)
      call add      ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
c
      call fold6    ("t3nw","s1  ",4,1)
      call fold6    ("s1  ","s2  ",5,1)
      call fold6    ("s2  ","twf ",6,1)
      call sqr3     ("twf ")
c
c   diag f1
      call reord7   ("t33 ","s1  ",231456)
      call reord6   ("s1  ","s2  ",124563)
      call mult6    ("s2  ","tvf ","s1  ",5)
      call reord6   ("s1  ","s2  ",126345)
      call reord7   ("s2  ","s1  ",312456)
      call add      ("s1  ",-0.5d0,"t3nw",1.d0,"t3nw")
c   diag f2
      call mult6    ("tuf ","tu3 ","s1  ",1)
      call reord6   ("s1  ","s2  ",124356)
      call fold6    ("t3nw","s2  ",idum,2)
c   diag f3
      call t1t2     ("t1v ","tu  ","t33 ","s1  ",123456,1)
      call t1t2     ("t1v ","tu  ","s1  ","s2  ",123456,2)
      call reord7   ("s2  ","s1  ",312456)
      call reord6   ("s1  ","s2  ",124356)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",125634)
      call mult6    ("tuf ","s2  ","s1  ",3)
      call add      ("s1  ",-0.5d0,"t3nw",1.d0,"t3nw")
c   diag f4
      call mult6    ("twf ","t1vr","s1  ",5)
      call cmult    ("s1  ",0.5d0,"s1  ")
      call fold6    ("t3nw","s1  ",idum,2)
c   diag f5
      call mult6    ("twf ","tum4","s1  ",3)
      call cmult    ("s1  ",0.5d0,"s1  ")
      call fold6    ("t3nw","s1  ",idum,2)
c   diag f6
      call cmult    ("tu  ",3.d0,"s1  ")
      call t1sqr    ("t1v ","t1v ","s1  ","s2  ")
      call t1t2     ("t1v ","s2  ","t33 ","s1  ",123456,3)
      call reord6   ("s1  ","s2  ",124356)
      call reord7   ("s2  ","s1  ",312456)
      call reord6   ("s1  ","s2  ",125346)
      call reord7   ("s2  ","s1  ",132456)
      call reord6   ("s1  ","s2  ",126345)
      call mult6    ("twf ","s2  ","s1  ",4)
      call add      ("s1  ",-1.d0/6.d0,"t3nw",1.d0,"t3nw")
c
      call zero     ("t3nw")
      call sqr3     ("t3nw")
      call dveps6   ("t3nw")
      call prtime   ("caltw3  ",0,1)
      return
      end
C  *********************************************************************
      subroutine endiag_csw
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
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

#if defined HSFS_DYN_ALLOC
      call check_modules
#endif

c  GROUP I
c   diag 1
      call reorder    ("pvvv","box1",2431,4)
c v1 Feb 95
      call fold4      ("tv  ","s1  ",3,4)
      call fold4      ("s1  ","s2  ",4,4)
      call reorder    ("s2  ","s1  ",1342,4)
      call mult       ("s1  ","ph  ","s2  ",311)
      call add        ("box1",1.d0,"s2  ",-.5d0,"box1")

#if defined HSFS_DYN_ALLOC
      call check_modules
#endif

c   diag 2
      call fold4      ("tv1 ","s1  ",2,4)
      call reorder    ("ppvh","s2  ",3142,4)
      call reorder    ("ppvh","s3  ",3241,4)
      call add        ("s2  ",2.d0,"s3  ",-1.d0,"s2  ")
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",1.d0,"box1",1.d0,"box1")
c   diag 3
      call fold4      ("tvm2","s1  ",2,4)
      call reorder    ("ppvh","s2  ",3142,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",-1.d0,"box1",1.d0,"box1")
c   diag 4
      call reorder    ("ppvh","s2  ",3241,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","s1  ",1324,4)
      call add        ("s1  ",-1.d0,"box1",1.d0,"box1")

      !call check_modules

c   diag 5
      call reorder    ("vhhp","s1  ",1432,4)
      call fold4      ("tcm ","s2  ",3,4)
      call fold4      ("s2  ","s3  ",4,4)
      call reorder    ("s3  ","s2  ",3412,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","s1  ",1342,4)
      call add        ("s1  ",1.d0,"box1",1.d0,"box1")
c   diag 6
      call fold4      ("tv  ","s1  ",3,4)
      call clear6     ("s2  ",322333,100101)
      call t1t2       ("t1c1","s1  ","s2  ","s3  ",123456,3)
      call fold4      ("s1  ","s2  ",4,4)
      call t1t2       ("t1c ","s2  ","s3  ","s1  ",123456,2)
      call reord6     ("s1  ","s2  ",124356)
      call reord7     ("s2  ","s1  ",132456)
      call reord6     ("s1  ","s2  ",126543)
      call reorder    ("hhpp","s1  ",3412,4)
      call reorder    ("hhpp","s3  ",3421,4)
      call add        ("s1  ",-2.d0,"s3  ",1.d0,"s1  ")
      call mult6      ("s2  ","s1  ","s3  ",2)
      call add        ("s3  ",1.d0,"box1",1.d0,"box1")
c  diag 7
      call fold4      ("tc  ","s1  ",3,4)
      call fold4      ("s1  ","s2  ",4,4)
      call clear6     ("s1  ",322333,100011)
      call t1t2       ("t1v ","s2  ","s1  ","s3  ",123456,1)
      call fold4      ("tv  ","s1  ",4,4)
      call t1t2       ("t1c1","s1  ","s3  ","s2  ",123456,2)
      call fold4      ("tvm ","s1  ",4,4)
      call t1t2       ("t1c1","s1  ","s2  ","s3  ",123456,3)
      call reord6     ("s3  ","s1  ",125346)
      call reord7     ("s1  ","s2  ",132456)
      call reord6     ("s2  ","s1  ",126543)
      call reorder    ("hhpp","s2  ",3412,4)
      call mult6      ("s1  ","s2  ","s3  ",2)
      call add        ("s3  ",1.d0,"box1",1.d0,"box1")
c  diag 8
      call reorder    ("vpvh","s1  ",1324,4)
      call fold4      ("t1cr","s2  ",1,2)
      call mult       ("s1  ","s2  ","s3  ",311)
      call reorder    ("s3  ","s1  ",1243,4)
      call add        ("s1  ",-1.d0,"box1",1.d0,"box1")
c  diag 9
      call reorder    ("vphv","s1  ",1423,4)
      call mult       ("s1  ","s2  ","s3  ",311)
      call reorder    ("s3  ","s1  ",1423,4)
      call add        ("s1  ",-1.d0,"box1",1.d0,"box1")
c  diag 10
      call reorder    ("ppvv","s1  ",3421,4)
      call mult       ("t1v ","s1  ","s2  ",131)
      call add        ("s2  ",1.d0,"box1",1.d0,"box1")
c
      if (iopt(17).lt.3) then
        call fold4      ("tu  ","s1  ",3,4)
        call mult6      ("s1  ","box1","s2  ",1)
        call reord6     ("s2  ","twf ",124365)
      else if (iopt(17).eq.3) then
        call wrtfile    (23,0,"box1")
      end if
c  GROUP II
c   diag 1
      call cmult      ("vvvh",-1.d0,"box2")
C v1 Feb 95 UK hp diagram
      call fold4      ("tu  ","s1  ",3,4)
      call mult       ("s1  ","hp  ","s2  ",311)
      call add        ("box2",1.d0,"s2  ",-.5d0,"box2")
cuk
c   diag 2
      call fold4      ("tv1 ","s1  ",2,4)
      call reorder    ("vhhp","s2  ",1324,4)
      call add        ("s2  ",-2.d0,"vhhp",1.d0,"s2  ")
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","s1  ",1324,4)
      call add        ("s1  ",1.d0,"box2",1.d0,"box2")
c   diag 3
      call fold4      ("tvm2","s1  ",2,4)
      call reorder    ("vhhp","s2  ",1324,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","s2  ",1324,4)
      call add        ("s2  ",1.d0,"box2",1.d0,"box2")
c   diag 4
      call mult       ("vhhp","s1  ","s2  ",222)
      call reorder    ("s2  ","s1  ",1342,4)
      call add        ("s1  ",1.d0,"box2",1.d0,"box2")
c   diag 5
      call reorder    ("ppvh","s1  ",3412,4)
      call mult       ("tum ","s1  ","s2  ",222)
      call add        ("s2  ",-1.d0,"box2",1.d0,"box2")
c   diag 6
      call reorder    ("hhpp","s1  ",1243,4)
      call add        ("s1  ",-2.d0,"hhpp",1.d0,"s1  ")
      call mult       ("t1v ","s1  ","s2  ",131)
      call fold4      ("tv1 ","s1  ",2,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","s1  ",1324,4)
      call add        ("s1  ",1.d0,"box2",1.d0,"box2")
c   diag 7
      call fold4      ("tv  ","s1  ",4,4)
      call clear6     ("s2  ",332333,110001)
      call t1t2       ("t1v ","s1  ","s2  ","s3  ",123456,1)
      call fold4      ("tv  ","s1  ",4,4)
      call t1t2       ("t1v ","s1  ","s3  ","s2  ",123456,2)
      call t1t2       ("t1c1","tum ","s2  ","s1  ",123456,3)
      call reord6     ("s1  ","s2  ",126345)
      call mult6      ("s2  ","hhpp","s1  ",2)
      call reorder    ("s1  ","s2  ",2134,4)
      call add        ("s2  ",1.d0,"box2",1.d0,"box2")
c   diag 8
      call reorder    ("vpvh","s1  ",1342,4)
      call mult       ("s1  ","t1v ","s2  ",311)
      call reorder    ("s2  ","s1  ",1423,4)
      call add        ("s1  ",-1.d0,"box2",1.d0,"box2")
c   diag 9
      call reorder    ("vphv","s1  ",1432,4)
      call mult       ("t1v ","s1  ","s2  ",131)
      call add        ("s2  ",-1.d0,"box2",1.d0,"box2")
c   diag 10
      call fold4      ("t1cr","s1  ",1,2)
      call mult       ("vvhh","s1  ","s2  ",311)
      call reorder    ("s2  ","s1  ",2143,4)
      call add        ("s1  ",1.d0,"box2",1.d0,"box2")

      !call check_modules
c
      if (iopt(17).lt.3) then
        call fold4      ("tv  ","s1  ",3,4)
        call fold4      ("s1  ","s2  ",4,4)
        call reorder    ("s2  ","s1  ",1432,4)
        call mult6      ("box2","s1  ","s2  ",1)
        call reord6     ("s2  ","s1  ",124356)
        call add        ("s1  ",1.d0,"twf ",1.d0,"twf ")
      else if (iopt(17).eq.3) then
        call wrtfile    (23,1,"box2")
      end if
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
c  diag e1
        call fold6    ("t32 ","s1  ",4,1)
        call fold6    ("s1  ","s2  ",5,1)
        call reord6   ("s2  ","s1  ",124536)
        call reorder  ("vpvh","s2  ",1342,4)
        call reorder  ("vphv","s3  ",1432,4)
        call add      ("s2  ",2.d0,"s3  ",-1.d0,"s2  ")
        call mult6    ("s1  ","s2  ","s3  ",3)
        call reord6   ("s3  ","s1  ",125346)
        call add      ("s1  ",0.5d0,"twf ",1.d0,"twf ")
c  diag e2
        call fold6    ("t32 ","s1  ",4,1)
        call fold6    ("s1  ","s2  ",6,1)
        call reord6   ("s2  ","s1  ",124635)
        call reorder  ("vpvh","s2  ",1342,4)
        call mult6    ("s1  ","s2  ","s3  ",3)
        call reord6   ("s3  ","s2  ",125346)
        call add      ("s2  ",-1.d0,"twf ",1.d0,"twf ")
c  diag e3
        call reorder  ("vphv","s2  ",1432,4)
        call mult6    ("s1  ","s2  ","s3  ",3)
        call reord6   ("s3  ","s1  ",125364)
        call add      ("s1  ",-1.d0,"twf ",1.d0,"twf ")
c  diag e4
        call fold6    ("t31 ","s1  ",4,1)
        call fold6    ("s1  ","s2  ",5,1)
        call fold6    ("s2  ","s1  ",6,1)
        call reord6   ("s1  ","s2  ",126345)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","s2  ",126543)
        call mult6    ("vvhh","s2  ","s1  ",3)
        call add      ("s1  ",0.5d0,"twf ",1.d0,"twf ")
c  diag e5
        call fold6    ("t33 ","s1  ",4,1)
        call reorder  ("ppvv","s2  ",3412,4)
        call mult6    ("s1  ","s2  ","s3  ",3)
        call add      ("s3  ",0.5d0,"twf ",1.d0,"twf ")
      end if
      if (iopt(17).ne.3) then
        !call check_modules
        call sqr3     ("twf ")
      end if
      call prtime     ("endiag  ",0,0)
      return
      end
C  ********************************************************************
      subroutine rdbox_csw
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call rdfile  (20,6,"bx7 ",-1)
      call rdfile  (21,0,"bx12",-1)
      call rdfile  (21,3,"bx15",-1)
      call rdfile  (21,4,"bx16",-1)
      call rdfile  (21,5,"bx17",-1)
      call rdfile  (21,7,"bx19",-1)
      call rdfile  (22,0,"bx24",-1)
      call rdfile  (22,1,"bx25",-1)
      return
      end
C  *******************************************************************
