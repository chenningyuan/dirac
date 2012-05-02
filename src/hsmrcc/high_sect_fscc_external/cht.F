!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine cht 
!---------------------------------------------------------------------
!
!
! some description...
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use cint_mod
      use dat_mod
      use cij1_mod
      use gen_mod
#include "implicit.h"
#include "priunit.h"
#include "csdtime.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     parameter (nbig=60000000,ndint=3000000,ndij=960000)
#include "param_ccsd.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)
#include "cfiles.h"
!     common /cint/ nint,idum1,vint(ndint),iint(ndint)
#include "cint.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /csdtime/ csdstim,csdltim
#include "csdtime.h"
#endif

      character*256 filenm,filenm6,
     1filenm51,filenm52,filenm53,filenm55,filenm56,filenm57,filenm58,
     2filenm20,filenm21,filenm31,filenm32,filenm41,filenm42
      real*4 cputime
      external fmax,fmax6
c     iopt(1) - print.
c     iopt(2) - no. of iterations.
c     iopt(3) - convergence = 10**(-iopt(3))
c     iopt(4) .ne.0 - do rle every iopt(4) iterations.
c     iopt(5) - print only elements .gt. 10**(-iopt(5))
c     iopt(6) - start avrg.
c     iopt(8) - fock.
c     iopt(9) - incomplete model space
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
      filenm58(1:1)=' '
      call getenv('tape58',filenm58)
      if (filenm58(1:1).eq.' ') filenm58='tape58'
      open(58,file=filenm58,form='UNFORMATTED')
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')
      filenm20(1:1)=' '
      call getenv('tape20',filenm20)
      if (filenm20(1:1).eq.' ') filenm20='tape20'
      open(20,file=filenm20,form='UNFORMATTED')
      filenm21(1:1)=' '
      call getenv('tape21',filenm21)
      if (filenm21(1:1).eq.' ') filenm21='tape21'
      open(21,file=filenm21,form='UNFORMATTED')
      filenm31(1:1)=' '
      call getenv('tape31',filenm31)
      if (filenm31(1:1).eq.' ') filenm31='tape31'
      open(31,file=filenm31,form='UNFORMATTED')
      filenm32(1:1)=' '
      call getenv('tape32',filenm32)
      if (filenm32(1:1).eq.' ') filenm32='tape32'
      open(32,file=filenm32,form='UNFORMATTED')
      filenm41(1:1)=' '
      call getenv('tape41',filenm41)
      if (filenm41(1:1).eq.' ') filenm41='tape41'
      open(41,file=filenm41,form='UNFORMATTED')
      filenm42(1:1)=' '
      call getenv('tape42',filenm42)
      if (filenm42(1:1).eq.' ') filenm42='tape42'
      open(42,file=filenm42,form='UNFORMATTED')
      write (*,*) '   $$$$$$$$$   cht   $$$$$$$$$  '
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
      call reads   (51,0,"s1  ")
      call reorder ("s1  ","tc  ",3412,4)
      call reads   (51,1,"s2  ")
      call reorder ("s2  ","t1c ",1243,4)
      call rdfile  (52,0,"tv  ",-1)
      call rdfile  (52,1,"t1v ",-1)
      call rdfile  (52,2,"tvf ",-1)
      if (iopt(17).ne.0) then
        call rdfile  (52,3,"t3v ",-1)
      end if
      call reads   (56,0,"s1  ")
      call reorder ("s1  ","tp  ",4321,4)
      call reads   (56,1,"s1  ")
      call reorder ("s1  ","tpf ",4321,4)
      call reads   (57,0,"s1  ")
      call reorder ("s1  ","tq  ",3412,4)
      call reads   (57,1,"s1  ")
      call reorder ("s1  ","tqf ",3412,4)
      call rdfile  (53,0,"tg  ",-1)
      call rdfile  (53,1,"t1g ",-1)    
      call rdfile  (53,2,"tgf ",-1)
      if (iopt(17).ne.0) then
        call rdfile  (53,3,"t3g ",-1)
      end if
      call rdfile  (55,0,"thu ",-1)
      call rdfile  (55,1,"thuf",-1)
      if (iopt(17).ne.0) then
        call reads   (56,2,"s2  ")
        call reord6  ("s2  ","s1  ",125346)
        call reord7  ("s1  ","s2  ",312456)
        call reord6  ("s2  ","s1  ",125346)
        call reord7  ("s1  ","s2  ",132456)
        call reord6  ("s2  ","t3p ",126435)
        call reads   (57,2,"s2  ")
        call reord6  ("s2  ","s1  ",124356)
        call reord7  ("s1  ","s2  ",312456)
        call reord6  ("s2  ","s1  ",125346)
        call reord7  ("s1  ","s2  ",132456)
        call reord6  ("s2  ","t3q ",126345)
        call rdfile  (55,2,"t3u ",-1)
      end if
      call t1sqr ("t1c ","t1c ","tc  ","tcm ")
      call t1sqr ("t1g ","t1c ","tg  ","tgm ")
      call reorder ("t1c ","s1  ",1243,4)
      call t1sqr ("t1v ","s1  ","tv  ","tvm ")
      call t1sqr ("t1g ","t1g ","thu ","thum")
      call reorder ("t1v ","s1  ",1243,4)
      call t1sqr ("t1g ","s1  ","tp  ","tpm ")
      call rdbox_cht
      call ordv_cht
      call files ("gv  ",num,new,1,-1,-1)
      if (ije1(1).gt.0) then
        call calvv10 
        call caltv1
        call fmax   ("t11 ",tmax,indmax)
        write (6,94) tmax,idec(indmax)
94      format (1x,'max t11 = ',f12.6,1x,i12)
        write (6,*) ' '
        call prt ("t11 ","t11 ",0)      
        write (6,*) ' '
      else
        call move("gv  ","t11 ")
      end if
      call ordta0_cht
      call calvv0_cht 
      call ordtb0_cht
      if (iopt(17).ne.3) then
        if (iopt(17).ne.0) then
          call calt12p
          call calt12q
        end if
        call caltv2_cht
        call fmax ("t12 ",tmax,indmax)
        write (6,99) tmax,idec(indmax)
99      format (1x,'max t12 = ',f12.6,1x,i12)
        write (6,*) ' '
        call prt ("t12 ","t12 ",0)      
        write (6,*) ' '
        call wrtfile (58,0,"t11 ",-1)
        call wrtfile (58,1,"t12 ",-1)
      else if (iopt(17).eq.3) then
        call calv30p 
        call calv30q
        call tinit_cht
        if (iopt(4).ne.0) then
          call dbdiag  ("dum1","ggvh","t12p","t12q",0111,-1)
          call wrtrle ("t12 ",2,1)    
          call wrtrle ("t12p",3,1)    
          call wrtrle ("t12q",4,1)    
        end if
        iter=0
        iquit=0
        do while (iquit.eq.0)
          iter=iter+1
          call prtime ("start iter",iter,0)
          if (iter.eq.iopt(18)) then
            iopt(4)=0
            write (6,*) 'rle switched off at iteration ',iter
          end if
          call calt12p
          call calt12q
          call dfoldp
          call dfoldq
          call caltv2_cht
          if (iopt(6).ne.0.and.iter.ge.mod(iopt(6),100)) then
            fnew=max(dfloat(iopt(6)/100),1.d0)
            if (fnew.le.6.5d0) then
              fnew=1.d0/fnew
            else
              fnew=fnew-5.d0
            endif
            call add("tnew",fnew/(1.d0+fnew),"t12 ",
     1      1.d0/(1.d0+fnew),"tnew")
            call add("tpnw",fnew/(1.d0+fnew),"t12p",
     1      1.d0/(1.d0+fnew),"tpnw")
            call add("tqnw",fnew/(1.d0+fnew),"t12q",
     1      1.d0/(1.d0+fnew),"tqnw")
            write (6,'('' t averaged with fnew '',f7.3)') fnew
          end if
          call fmax ("tnew",tmax,itmax)
          call fmax6 ("tpnw",tpmax,idpmx,indpmx)
          call fmax6 ("tqnw",tqmax,idqmx,indqmx)
          call diffmax ("t12 ","tnew","s2  ",diff,idiff)
          call dfmax3 ("t12p","tpnw","s1  ",dfp,iddfp,inddfp)
          call dfmax3 ("t12q","tqnw","s1  ",dfq,iddfq,inddfq)
          write (6,*) ' '
          write (6,'('' iteration # '',i3,'' diffmx='',e10.3,2x,i12,
     1    e10.3,2x,i12)') iter,diff,idec(idiff)
          write (6,98) dfp,idec(iddfp),idec(inddfp),
     1    dfq,idec(iddfq),idec(inddfq)
98        format (2(1x,e10.3,2x,i8,2x,i8))
          write (6,97) tmax,idec(itmax)
97        format(' max t12, t12p, t12q =',f15.9,2x,i12)
          write (6,95) tpmax,idec(idpmx),idec(indpmx),
     1    tqmax,idec(idqmx),idec(indqmx)
95        format (2(1x,f15.9,2x,i8,2x,i8))
          if (abs(tmax).gt.1.d10.or.abs(tmax).lt.1.d-10) then
            write (*,*) ' abort. tmax out of bounds'
            call abext('tmax    ')    
          endif
          call move ("tnew","t12 ")
          call t1sqr ("t11 ","t1g ","t12 ","tm  ")
          call move ("tpnw","t12p")
          call move ("tqnw","t12q")
          if (max(abs(diff),abs(dfp),abs(dfq)).lt.convrg) then      
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
              call wrtrle ("t12 ",2,itr1+1)
              call wrtrle ("t12p",3,itr1+1)
              call wrtrle ("t12q",4,itr1+1)
              if (itr1.eq.iopt(4)) then      
                call rle (itr1,0111,"dum1","t12 ","t12p","t12q",ising)
                if (ising.eq.1) then
                  iopt(4)=0
                  write (6,*) ' ***** rle singular, turned off'
                endif
                call wrtrle ("t12 ",2,1)
                call wrtrle ("t12p",3,1)
                call wrtrle ("t12q",4,1)
              end if
            end if
          end if
        end do
        if (iparm.eq.1) then
          write (6,*) ' t12 converged after iter ',iter
        else if (iparm.eq.2) then
          write (6,*) ' iteration limit'
          call abext('limit   ')     
        end if
        write (6,*) 'last files'
        call prt ("t12 ","t12 ",0)
        call prt6 ("t12p","t12p",0)
        call prt6 ("t12q","t12q",0)
        call prt6 ("ttpf","ttpf",0)
        call prt6 ("ttqf","ttqf",0)
        call move ("tgf ","thvf")
        egr0=0.d0
        call eneru2 (egr0,1,2,max(1,iopt(1)))
        call eneru2 (egr0,12,2,max(1,iopt(1)))
      end if
      write (6,*) nstrt(1,nflcur),' words used in igen,',
     1nstrt(2,nflcur),' words used in gen'
      call prtime ("end ",0,0)
!     stop
      end

C ********************************************************************
      subroutine ordv_cht
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
!     common /cint/nint,idum1,vint(ndint),iint(ndint)
#include "cint.h"
!     common /dat/ norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)
#include "fock.h"
#endif

      call readv(nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call prtime("readv   ",nint,0)
      call selt1    (32,0,"ph  ",4)
      call selt1    (22,0,"hh  ",4)
      call selt1    (32,10,"vh  ",4)
      call selt1    (23,11,"gv  ",4)
      call selt1    (23,10,"gp  ",4)
      call seltyp   (3322,0,"pphh",4)
      call seltyp   (3222,0,"phhh",4)
      call seltyp   (3322,1000,"vphh",4)
      call seltyp   (3222,1000,"vhhh",4)
      call seltyp   (2322,1100,"gvhh",4)
      call seltyp   (2332,1010,"gpvh",4)
      call seltyp   (2332,1000,"gpph",4)
      call seltyp   (2322,1000,"gphh",4)
      call seltyp   (2232,1010,"ghvh",4)
      call seltyp   (2233,1001,"ghpv",4)
      call seltyp   (2232,1000,"ghph",4)
      call seltyp   (2233,1110,"ggvp",4)
      call seltyp   (2232,1110,"ggvh",4)
      call seltyp   (2233,1100,"ggpp",4)
      call seltyp   (2232,1100,"ggph",4)
      call reorder  ("pphh","va1 ",1324,4)
      call reorder  ("pphh","va2 ",1423,4)
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
        call seltyp   (2223,0,"hhhp",4)
        call seltyp   (2222,1111,"gggg",4)
        call seltyp   (2223,1100,"gghp",4)
        call seltyp   (2222,1000,"ghhh",4)
        call seltyp   (2323,1000,"gphp",4)
        call seltyp   (2333,1010,"gpvp",4)
        call seltyp   (2323,1111,"gvgv",4)
        call seltyp   (2332,1111,"gvvg",4)
        call seltyp   (3233,1100,"vgpp",4)
        call seltyp   (3223,1000,"vhhp",4)
        call seltyp   (3232,1000,"vhph",4)
        call seltyp   (3233,1000,"vhpp",4)
      end if
      if (iopt(17).gt.0) then
        call seltyp    (2233,1000,"ghpp",4)
      end if
      if (iopt(17).eq.3) then
        call seltyp    (2222,0,"hhhh",4)
        call seltyp    (2233,1101,"ggpv",4)
        call seltyp    (2333,1000,"gppp",4)
        call reorder   ("pphh","va6 ",1243,4)
      end if
      return
      end 
C  ******************************************************************
      subroutine calvv10
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
!     common /dat/ norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

c  diag 1
      call move      ("gv  ","vv10")
c  diag 2
      call reorder   ("vphh","s1  ",1243,4)
      call add       ("vphh",-2.d0,"s1  ",1.d0,"s1  ")
      call mult      ("tgm ","s1  ","s2  ",113)
      call add       ("s2  ",1.d0,"vv10",1.d0,"vv10")
c  diag 3
      call reorder   ("tq  ","s1  ",1324,4)
      call mult      ("s1  ","ph  ","s2  ",202)
      call add       ("s2  ",2.d0,"vv10",1.d0,"vv10")
c  diag 4
      call reorder   ("tpm ","s1  ",1423,4)
      call mult      ("s1  ","ph  ","s2  ",202)
      call add       ("s2  ",-1.d0,"vv10",1.d0,"vv10")
c  diag 5
      call mult      ("t1g ","vh  ","s1  ",111)
      call add       ("s1  ",-1.d0,"vv10",1.d0,"vv10")
c  diag 6
      call reorder   ("gpph","s1  ",1432,4)
      call reorder   ("gpph","s2  ",1423,4)
      call add       ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult      ("s1  ","tvm ","s2  ",113)
      call add       ("s2  ",1.d0,"vv10",1.d0,"vv10")
c  diag 7
      call reorder   ("gpvh","s1  ",1324,4)
      call reorder   ("ghpv","s2  ",1432,4)
      call add       ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult      ("s1  ","t1c ","s2  ",202)
      call add       ("s2  ",1.d0,"vv10",1.d0,"vv10")
c  diag 8
      call mult      ("gp  ","t1v ","s1  ",111)
      call add       ("s1  ",1.d0,"vv10",1.d0,"vv10")
c  diag 9
      call reorder   ("pphh","s1  ",3412,4)
      call reorder   ("pphh","s2  ",3421,4)
      call add       ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call mult      ("t1v ","s1  ","s2  ",131)
      call reorder   ("tgm ","s1  ",1432,4)
      call mult      ("s1  ","s2  ","s3  ",113)
      call add       ("s3  ",1.d0,"vv10",1.d0,"vv10")
c  diag 10
      call reorder   ("pphh","s1  ",1243,4)
      call add       ("pphh",-2.d0,"s1  ",1.d0,"s1  ")
      call mult      ("t1g ","s1  ","s2  ",131)
      call reorder   ("tv  ","s1  ",1432,4)
      call mult      ("s2  ","s1  ","s3  ",113)
      call add       ("s3  ",1.d0,"vv10",1.d0,"vv10")
c  diags 11,12
      call reorder   ("pphh","s1  ",1324,4)
      call reorder   ("pphh","s2  ",1423,4)
      call add       ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult      ("s1  ","t1c ","s3  ",202)
      call reorder   ("tq  ","s1  ",1324,4)
      call reorder   ("tp  ","s2  ",1423,4)
      call add       ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult      ("s1  ","s3  ","s2  ",202)
      call add       ("s2  ",1.d0,"vv10",1.d0,"vv10")
      if (iopt(17).ne.0) then
c  diags 13,14
        call reord6    ("t3q ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("t3p ","s1  ",125346)
        call reord7    ("s1  ","s3  ",132456)
        call reord6    ("s3  ","s1  ",124365)
        call add       ("s2  ",1.d0,"s1  ",-1.d0,"s1  ")
        call reorder   ("pphh","s2  ",1243,4)
        call add       ("pphh",2.d0,"s2  ",-1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",0)
        call add       ("s3  ",1.d0,"vv10",1.d0,"vv10")
      end if
      call prtime    ("calvv10 ",0,0)
      return
      end
C  **********************************************************************
      subroutine caltv1
!---------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use cij1_mod
      use gen_mod
#include "implicit.h"
#include "priunit.h"
#include "strt.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     parameter (ndij=960000,nbig=60000000)
#include "param_ccsd.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1iiia(200),iiva(200)
#include "cfiles.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!    1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!    2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
#include "strt.h"
#endif

      external find2,find4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! mi- TODO: dynami allocs !!!!!!!!!!!!!!!!!!!!!!!!!
      dimension a(500,500),b(500),indlst(500),ipvt(500),z(500)

      do 10 i=1,500
        b(i)=0.d0
        do 20 j=1,500
          a(i,j)=0.d0
20      continue
10    continue
      call files ("gv  ",ntape,new,1,-1,-1)
      istrt1=nstrt(1,ntape-1)+2*nij1+4
      nequn=ije1(nij1)
      do 30 i=1,nequn
        indlst(i)=igen(istrt1+i)
30    continue
      call files ("vv10",ntape,new,1,-1,-1)
      do 60 i=1,nequn
        b(i)=-1.d0*gen(nstrt(2,ntape-1)+i)
60    continue
      do 40 i=1,nequn 
        idi=indlst(i)
        do 50 j=1,nequn
          idj=indlst(j)
c   diag f1
          call files ("tgf ",ntape,new,1,-1,-1)
          nsize1=ije1(nij1)
          istrt11=nstrt(1,ntape-1)+2*nij1+4
          istrt12=nstrt(2,ntape-1)
          if (iand(idi,z'00ff').eq.iand(idj,z'00ff')) then
            ii1=iand(ishft(idi,-8),z'ff')
            ii2=iand(ishft(idj,-8),z'ff')
            a(i,j)=a(i,j)+find2(ii1,ii2)
          end if
c   diag f2
          call files ("tvf ",ntape,new,1,-1,-1)
          nsize1=ije1(nij1)
          istrt11=nstrt(1,ntape-1)+2*nij1+4
          istrt12=nstrt(2,ntape-1)
          if (iand(idi,z'ff00').eq.iand(idj,z'ff00')) then
            ii1=iand(idi,z'ff')
            ii2=iand(idj,z'ff')
            a(i,j)=a(i,j)-find2(ii1,ii2)
          end if
c   diag f3
          call files ("tpf ",ntape,new,1,-1,-1)
          nsize2=ije1(nij1)
          istrt21=nstrt(1,ntape-1)+2*nij1+4
          istrt22=nstrt(2,ntape-1)
          ii1=iand(ishft(idi,-8),z'ff')
          ii2=iand(idj,z'ff')
          ii3=iand(ishft(idj,-8),z'ff')
          ii4=iand(idi,z'ff')
          a(i,j)=a(i,j)+find4(ii1,ii2,ii3,ii4)
c   diag f4
          call files ("tqf ",ntape,new,1,-1,-1)
          nsize2=ije1(nij1)
          istrt21=nstrt(1,ntape-1)+2*nij1+4
          istrt22=nstrt(2,ntape-1)
          ii1=iand(ishft(idi,-8),z'ff')
          ii2=iand(idj,z'ff')
          ii3=iand(idi,z'ff')
          ii4=iand(ishft(idj,-8),z'ff')
          a(i,j)=a(i,j)-2.d0*find4(ii1,ii2,ii3,ii4)
50      continue
40    continue
      call clear ("s1  ",23,11)      
      call files ("s1  ",ntape,new,1,-1,-1)
      nsize=ije1(nij1)
      istrt=nstrt(2,ntape-1)
      call dcopy (nsize,-1.d0,0,gen(istrt+1),1)
      call diveps ("s1  ")
      call cmult ("s1  ",-1.d0,"s1  ")
      do 80 i=1,nequn
        a(i,i)=a(i,i)+1.d0/gen(istrt+i)
80    continue
      call dgeco (a,500,nequn,ipvt,rcond,z)
      if (abs(rcond).lt.1.d-16) then
        write (6,*) 'Linear equations for t11 are singular.'
        stop
      end if
      call dgesl (a,500,nequn,ipvt,b,0)
      call clear ("t11 ",23,11)       
      call files ("t11 ",ntape,new,1,-1,-1)
      nsize=ije1(nij1)
      istrt=nstrt(2,ntape-1)
      do 70 i=1,nsize
        gen(istrt+i)=b(i)
70    continue
      call prtime ("caltv1  ",0,0)
      write (6,*) ' '
      return
      end
C  *********************************************************************
      subroutine calv30p 
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  GROUP I
c   diag 1
      call reads     (42,0,"s1  ")
      call reorder   ("tg  ","s2  ",1432,4)
      call mult6     ("s2  ","s1  ","s3  ",1)
      call reord6    ("s3  ","s1  ",126345)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","v30p",126435)
c   diag 2
      call reads     (31,2,"s1  ")
      call reorder   ("tp  ","s2  ",1432,4)
      call mult6     ("s2  ","s1  ","s3  ",1)
      call reord6    ("s3  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",125463)
      call add       ("s1  ",1.d0,"v30p",1.d0,"v30p")
c   diag 3
      call reads     (32,0,"s1  ")
      call mult6     ("s1  ","tv  ","s2  ",1)
      call reord6    ("s2  ","s1  ",126354)
      call add       ("s1  ",1.d0,"v30p",1.d0,"v30p")
c  GROUP II
c   diag 1
      call reads     (21,2,"s1  ")
      call mult6     ("thu ","s1  ","s2  ",1)
      call reord6    ("s2  ","s1  ",126354)
      call add       ("s1  ",1.d0,"v30p",1.d0,"v30p")
c   diag 2
      call reads     (42,1,"s1  ")
      call mult6     ("tg  ","s1  ","s2  ",1)
      call reord6    ("s2  ","s1  ",126354)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30p",1.d0,"v30p")
c   diag 3
      call reads     (31,0,"s1  ")
      call reorder   ("tp  ","s2  ",1243,4)
      call mult6     ("s2  ","s1  ","s3  ",1)
      call reord6    ("s3  ","s1  ",124653)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30p",1.d0,"v30p")
c  GROUP III
c   diag 1
      call reads     (21,3,"s1  ")
      call reord6    ("t3u ","s2  ",124563)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",126345)
      call add       ("s1  ",0.5d0,"v30p",1.d0,"v30p")
c   diag 2
      call reads     (31,3,"s1  ")
      call reord6    ("t3p ","s2  ",124536)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125364)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30p",1.d0,"v30p")
c  GROUP IV
c   diag 1
      call reads     (21,4,"s1  ")
      call reord6    ("t3u ","s2  ",124653)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",126345)
      call add       ("s1  ",1.d0,"v30p",1.d0,"v30p")
c   diag 2
      call reads     (31,4,"s1  ")
      call reord7    ("t3p ","s2  ",132456)
      call reord6    ("s2  ","s3  ",124536)
      call mult6     ("s3  ","s1  ","s2  ",3)
      call reord6    ("s2  ","s3  ",125364)
      call reord7    ("s3  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30p",1.d0,"v30p")
c   diag 3
      call reord6    ("t3p ","s2  ",125634)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125463)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30p",1.d0,"v30p")
c  GROUP V
c   diag 2
      call reads     (42,2,"s1  ")
      call reord6    ("t3g ","s2  ",124653)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",126345)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30p",1.d0,"v30p")
c   diag 3
      call reads     (31,5,"s1  ")
      call reord6    ("t3p ","s2  ",125634)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125643)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30p",1.d0,"v30p")
c  GROUP VI
c   diag 1
      call reads     (21,6,"s1  ")
      call reorder   ("s1  ","s2  ",1243,4)
      call mult6     ("t3u ","s2  ","s1  ",3)
      call reord6    ("s1  ","s2  ",123465)
      call add       ("s2  ",1.d0,"v30p",1.d0,"v30p")
c  diag 3
      call mult      ("t1g ","hhhh","s1  ",131)
      call reorder   ("tp  ","s2  ",1243,4)
      call mult6     ("s1  ","s2  ","s3  ",1)
      call reord6    ("s3  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",125436)
      call add       ("s1  ",1.d0,"v30p",1.d0,"v30p")
c  GROUP VII
c   diag 1
      call reads     (31,6,"s1  ")
      call reord6    ("t3p ","s2  ",124356)
      call reord7    ("s2  ","s3  ",132456)
      call reord6    ("s3  ","s2  ",125643)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125346)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",126354)
      call add       ("s1  ",1.d0,"v30p",1.d0,"v30p")
c   diag 2
      call reads     (32,1,"s1  ")
      call mult6     ("s1  ","t3v ","s2  ",3)
      call reord6    ("s2  ","s1  ",126453)
      call add       ("s1  ",1.d0,"v30p",1.d0,"v30p")
c   diag 3
      call mult      ("gppp","t1v ","s1  ",311)
      call reorder   ("s1  ","s2  ",1243,4)
      call reorder   ("tg  ","s1  ",1432,4)
      call mult6     ("s2  ","s1  ","s3  ",1)
      call reord6    ("s3  ","s1  ",124563)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30p",1.d0,"v30p")
c  GROUP VIII
c   diag 2
      call reads     (31,7,"s1  ")
      call reord6    ("t3p ","s2  ",124563)
      call mult6     ("s2  ","s1  ","s3  ",5)
      call reord6    ("s3  ","s1  ",126354)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30p",1.d0,"v30p")
c  GROUP IX  
c   diag 2
      call reads     (21,7,"s1  ")
      call mult6     ("t3u ","s1  ","s2  ",5)
      call add       ("s2  ",0.5d0,"v30p",1.d0,"v30p")
      call prtime    ("calv30p ",0,0)
      return
      end
C  *******************************************************************
      subroutine calv30q 
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  GROUP I
c   diag 1
      call reads     (31,1,"s1  ")
      call reorder   ("tq  ","s2  ",1432,4)
      call mult6     ("s1  ","s2  ","s3  ",1)
      call reord6    ("s3  ","s2  ",124563)
      call reord7    ("s2  ","v30q",132456)
c   diag 3
      call reorder   ("tp  ","s2  ",1432,4)
      call mult6     ("s2  ","s1  ","s3  ",1)
      call reord6    ("s3  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",125436)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c   diag 2
      call reads     (31,2,"s1  ")
      call reorder   ("tq  ","s2  ",1432,4)
      call mult6     ("s1  ","s2  ","s3  ",1)
      call reord6    ("s3  ","s1  ",124365)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30q",1.d0,"v30q")
c   diag 4
      call reads     (32,0,"s1  ")
      call reorder   ("tv  ","s2  ",1243,4)
      call mult6     ("s1  ","s2  ","s3  ",1)
      call reord6    ("s3  ","s1  ",126345)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c   diag 5
      call reads     (42,3,"s1  ")
      call reorder   ("tg  ","s2  ",1432,4)
      call mult6     ("s2  ","s1  ","s3  ",1)
      call reord6    ("s3  ","s1  ",125346)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",126453)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c  diag 6
      call box38     (0)
      call reorder   ("tc  ","s1  ",1432,4)
      call mult6     ("bx38","s1  ","s2  ",1)
      call reord6    ("s2  ","s1  ",124536)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c  GROUP II
c   diag 2
      call reads     (31,0,"s1  ")
      call mult6     ("s1  ","tq  ","s2  ",1)
      call reord6    ("s2  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",125364)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c   diag 3
      call reads     (21,1,"s1  ")
      call mult6     ("thu ","s1  ","s2  ",1)
      call reord6    ("s2  ","s1  ",126345)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c   diag 4
      call reads     (42,1,"s1  ")
      call reorder   ("s1  ","s2  ",3214,4)
      call reorder   ("tg  ","s1  ",1243,4)
      call mult6     ("s2  ","s1  ","s3  ",1)
      call reord6    ("s3  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",125346)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c   diag 5
      call reads     (42,4,"s1  ")
      call mult6     ("tg  ","s1  ","s2  ",1)
      call reord6    ("s2  ","s3  ",126345)
      call reord7    ("s3  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30q",1.d0,"v30q")
c  diag 6
      call reorder   ("tg  ","s2  ",1243,4)
      call mult6     ("s2  ","s1  ","s3  ",1)
      call reord6    ("s3  ","s1  ",126543)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30q",1.d0,"v30q")
c  GROUP III
c   diag 2
      call reads     (31,3,"s1  ")
      call reord6    ("t3q ","s2  ",124536)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord7    ("s3  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",124356)
      call add       ("s2  ",1.d0,"v30q",1.d0,"v30q")
c   diag 3
      call reads     (42,5,"s1  ")
      call reord6    ("t3g ","s2  ",124563)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",126354)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30q",1.d0,"v30q")
c  GROUP IV 
c   diag 2
      call reads     (31,4,"s1  ")
      call reord6    ("t3q ","s2  ",124635)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord7    ("s3  ","s2  ",132456)
      call reord6    ("s2  ","s3  ",124356)
      call add       ("s3  ",1.d0,"v30q",1.d0,"v30q")
c   diag 4
      call reord7    ("t3p ","s2  ",132456)
      call reord6    ("s2  ","s3  ",125634)
      call mult6     ("s1  ","s3  ","s2  ",3)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",124356)
      call add       ("s2  ",1.d0,"v30q",1.d0,"v30q")
c   diag 5
      call reads     (42,6,"s1  ")
      call reord6    ("t3g ","s2  ",124653)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s2  ",126354)
      call reord7    ("s2  ","s3  ",132456)
      call add       ("s3  ",1.d0,"v30q",1.d0,"v30q")
c   diag 6
      call reord6    ("t3g ","s2  ",125643)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",126453)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30q",1.d0,"v30q")
c  GROUP V
c   diag 2
      call reads     (31,5,"s1  ")
      call reord6    ("t3q ","s2  ",124635)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord7    ("s3  ","s2  ",132456)
      call reord6    ("s2  ","s3  ",124653)
      call add       ("s3  ",1.d0,"v30q",1.d0,"v30q")
c   diag 5
      call reord7    ("t3p ","s2  ",132456)
      call reord6    ("s2  ","s3  ",124536)
      call mult6     ("s3  ","s1  ","s2  ",3)
      call reord6    ("s2  ","s3  ",125346)
      call reord7    ("s3  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30q",1.d0,"v30q")
c   diag 6
      call reord7    ("t3p ","s2  ",132456)
      call reord6    ("s2  ","s3  ",125634)
      call mult6     ("s3  ","s1  ","s2  ",3)
      call reord6    ("s2  ","s3  ",125634)
      call reord7    ("s3  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30q",1.d0,"v30q")
c   diag 3
      call reads     (21,5,"s1  ")
      call reord6    ("t3u ","s2  ",124653)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",126354)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c   diag 4
      call reads     (42,2,"s1  ")
      call reorder   ("s1  ","s2  ",2134,4)
      call reord6    ("t3g ","s1  ",125643)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord7    ("s3  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",124635)
      call add       ("s2  ",1.d0,"v30q",1.d0,"v30q")
c  GROUP VI
c   diag 2
      call reads     (21,6,"s1  ")
      call mult6     ("t3u ","s1  ","s2  ",3)
      call add       ("s2  ",1.d0,"v30q",1.d0,"v30q")
c   diag 3
      call reord6    ("t3u ","s2  ",123654)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",123654)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c   diag 4
      call mult      ("t1g ","hhhh","s1  ",131)
      call mult6     ("s1  ","tq  ","s2  ",1)
      call reord6    ("s2  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",125463)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c  GROUP VII
c   diag 1
      call reads     (31,6,"s1  ")
      call reord6    ("t3q ","s2  ",124356)
      call reord7    ("s2  ","s3  ",132456)
      call reord6    ("s3  ","s2  ",125643)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord7    ("s3  ","s2  ",132456)
      call reord6    ("s2  ","s3  ",123645)
      call add       ("s3  ",1.d0,"v30q",1.d0,"v30q")
c   diag 2
      call reord6    ("t3p ","s2  ",124356)
      call reord7    ("s2  ","s3  ",132456)
      call reord6    ("s3  ","s2  ",125634)
      call mult6     ("s2  ","s1  ","s3  ",3)
      call reord6    ("s3  ","s1  ",125346)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",126345)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c   diag 3
      call box25a
      call reord6    ("t3v ","s1  ",123654)
      call mult6     ("b25a","s1  ","s2  ",3)
      call reord6    ("s2  ","s1  ",126435)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c   diag 4
      call reorder   ("gppp","s1  ",1243,4)
      call mult      ("s1  ","t1v ","s2  ",311)
      call reorder   ("s2  ","s1  ",1243,4)
      call reorder   ("tg  ","s2  ",1432,4)
      call mult6     ("s2  ","s1  ","s3  ",1)
      call reord6    ("s3  ","s1  ",124356)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",125463)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c   diag 5
      call mult      ("ggpp","t1v ","s1  ",311)
      call reorder   ("s1  ","s2  ",1243,4)
      call reorder   ("tc  ","s1  ",1432,4)
      call mult6     ("s2  ","s1  ","s3  ",1)
      call reord6    ("s3  ","s1  ",124536)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
c  GROUP VIII
c   diag 2
      call reads     (31,7,"s1  ")
      call reord6    ("t3q ","s2  ",124563)
      call mult6     ("s2  ","s1  ","s3  ",5)
      call reord6    ("s3  ","s2  ",126534)
      call reord7    ("s2  ","s3  ",312456)
      call add       ("s3  ",1.d0,"v30q",1.d0,"v30q")
c   diag 3
      call reord7    ("t3p ","s2  ",132456)
      call reord6    ("s2  ","s3  ",124563)
      call mult6     ("s3  ","s1  ","s2  ",5)
      call reord6    ("s2  ","s1  ",126345)
      call reord7    ("s1  ","s2  ",132456)
      call add       ("s2  ",1.d0,"v30q",1.d0,"v30q")
c  GROUP IX
c   diag 3
      call reads     (21,7,"s1  ")
      call reord6    ("t3u ","s2  ",123465)
      call mult6     ("s2  ","s1  ","s3  ",5)
      call reord6    ("s3  ","s1  ",123465)
      call add       ("s1  ",1.d0,"v30q",1.d0,"v30q")
      call prtime    ("calv30q ",0,0)
      return
      end
C  *******************************************************************
      subroutine ordta0_cht
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
!     common /dat/ norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      call clear      ("tb1 ",2322,1110)
      call clear      ("tb2 ",2223,1011)
      call clear      ("tb3 ",2223,1101)
      call clear      ("tb4 ",2322,1101)
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
        call reorder  ("tg  ","tg6 ",1243,4)
        call reorder  ("tp  ","tp3 ",1432,4)
        call reorder  ("tq  ","tq3 ",1432,4)
      end if
      return
      end 
C  *********************************************************************
      subroutine calvv0_cht
!---------------------------------------------------------------------
!
! add description
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
!     common /dat/ norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

c  diag 1
      call move      ("ggvh","vv0 ")
c  diag 2b
      call reorder   ("thu ","s1  ",1243,4)
      call mult      ("s1  ","vh  ","s2  ",311)
      call add       ("s2  ",-1.d0,"tb3 ",1.d0,"tb3 ")
c  diag 3a
      call reorder   ("tq  ","s1  ",1342,4)
      call mult      ("s1  ","gp  ","s2  ",311)
      call add       ("s2  ",1.d0,"tb4 ",1.d0,"tb4 ")
c  diag 3b
      call reorder   ("tp  ","s1  ",1432,4)
      call mult      ("gp  ","s1  ","s2  ",131)
      call add       ("s2  ",1.d0,"vv0 ",1.d0,"vv0 ")
c  diag 4a
      call reorder   ("tq  ","s1  ",1324,4)
      call reorder   ("gphh","s2  ",1324,4)
      call add       ("s2  ",2.d0,"ghph",-1.d0,"s2  ")
      call mult      ("s1  ","s2  ","s3  ",222)
      call add       ("s3  ",1.d0,"tb1 ",1.d0,"tb1 ")
c  diag 4b
      call reorder   ("tg  ","s1  ",1324,4)
      call reorder   ("gpvh","s2  ",1324,4)
      call reorder   ("ghpv","s3  ",1432,4)
      call add       ("s2  ",2.d0,"s3  ",-1.d0,"s2  ")
      call mult      ("s2  ","s1  ","s3  ",222)
      call add       ("s3  ",1.d0,"tb1 ",1.d0,"tb1 ")
c  diag 5a
      call reorder   ("tpm ","s1  ",1423,4)
      call reorder   ("gphh","s2  ",1324,4)
      call mult      ("s1  ","s2  ","s3  ",222)
      call add       ("s3  ",-1.d0,"tb1 ",1.d0,"tb1 ")
c  diag 5b
      call reorder   ("tgm ","s1  ",1423,4)
      call reorder   ("gpvh","s2  ",1324,4)
      call mult      ("s2  ","s1  ","s3  ",222)
      call add       ("s3  ",-1.d0,"tb1 ",1.d0,"tb1 ")
c  diag 6a
      call reorder   ("tgm ","s1  ",1423,4)
      call reorder   ("ghpv","s2  ",1432,4)
      call mult      ("s1  ","s2  ","s3  ",222)
      call add       ("s3  ",-1.d0,"tb2 ",1.d0,"tb2 ")
c  diag 6b
      call reorder   ("tpm ","s1  ",1423,4)
      call mult      ("ghph","s1  ","s2  ",222)
      call add       ("s2  ",-1.d0,"tb2 ",1.d0,"tb2 ")
c  diag 7 
      call mult      ("thum","vhhh","s1  ",222)
      call add       ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
c  diag 8 
      call mult      ("ggpp","tvm ","s1  ",222)
      call add       ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
c  diag 12
      call reorder   ("tpm ","s1  ",1423,4)
      call mult      ("s1  ","va2 ","s2  ",222)
      call reorder   ("tgm ","s1  ",1423,4)
      call mult      ("s1  ","s2  ","s3  ",222)
      call add       ("s3  ",1.d0,"tb2 ",1.d0,"tb2 ")
c  diag 15
      call mult      ("thum","pphh","s1  ",222)
      call mult      ("s1  ","tvm ","s2  ",222)
      call add       ("s2  ",1.d0,"vv0 ",1.d0,"vv0 ")
c  
      call boxZ1
      call reorder   ("tg  ","s1  ",1324,4)
      call mult      ("bxZ1","s1  ","s2  ",222)
      call add       ("s2  ",1.d0,"tb1 ",1.d0,"tb1 ")
c  
      call boxZ2
      call reorder   ("tgm ","s1  ",1423,4)
      call mult      ("bxZ2","s1  ","s2  ",222)
      call add       ("s2  ",1.d0,"tb1 ",1.d0,"tb1 ")
c  
      call boxZ3
      call reorder   ("tq  ","s1  ",1342,4)
      call mult      ("s1  ","bxZ3","s2  ",311)
      call add       ("s2  ",1.d0,"tb4 ",1.d0,"tb4 ")
      call reorder   ("tp  ","s1  ",1432,4)
      call mult      ("bxZ3","s1  ","s2  ",131)
      call add       ("s2  ",1.d0,"vv0 ",1.d0,"vv0 ")
c  
      call boxZ4
      call reorder   ("thu ","s1  ",1243,4)
      call mult      ("s1  ","bxZ4","s2  ",311)
      call add       ("s2  ",1.d0,"tb3 ",1.d0,"tb3 ")
c  
      call boxZ5
      call mult      ("t1g ","bxZ5","s1  ",131)
      call add       ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
c  
      call boxZ6
      call mult      ("bxZ6","t1g ","s2  ",311)
      call add       ("s2  ",1.d0,"tb4 ",1.d0,"tb4 ")
c  
      call boxZ7
      call mult      ("bxZ7","t1v ","s1  ",311)
      call add       ("s1  ",1.d0,"tb3 ",1.d0,"tb3 ")
c
      call boxZ8
      call reorder   ("t1c ","s1  ",1243,4)
      call mult      ("bxZ8","s1  ","s2  ",311)
      call add       ("s2  ",1.d0,"tb1 ",1.d0,"tb1 ")
c  
      call boxZ9
      call mult      ("bxZ9","t1v ","s1  ",311)
      call add       ("s1  ",1.d0,"tb3 ",1.d0,"tb3 ")
c  
      call boxZ0
      call reorder   ("t1c ","s1  ",1243,4)
      call mult      ("bxZ0","s1  ","s2  ",311)
      call add       ("s2  ",1.d0,"vv0 ",1.d0,"vv0 ")
c  diag f4a
      call reorder   ("t1g ","s1  ",1243,4)
      call mult      ("thuf","s1  ","s2  ",311)
      call reorder   ("s2  ","s1  ",1243,4)
      call reorder   ("t11 ","s2  ",1243,4)
      call mult      ("s1  ","s2  ","s3  ",311)
      call add       ("s3  ",-1.d0,"tb3 ",1.d0,"tb3 ")
c  diag f5 
      call reorder   ("thuf","s1  ",1243,4)
      call reorder   ("t11 ","s2  ",1243,4)
      call mult      ("s1  ","s2  ","s3  ",311)
      call unfold_ccsd("s3  ","s1  ",3)
      call add       ("s1  ",1.d0,"tb3 ",1.d0,"tb3 ")
c  diag f6a
      call reorder   ("tpf ","s1  ",1432,4)
      call mult      ("t11 ","s1  ","s2  ",131)
      call reorder   ("t1g ","s1  ",1243,4)
      call mult      ("s2  ","s1  ","s3  ",311)
      call add       ("s3  ",1.d0,"vv0 ",1.d0,"vv0 ")
c  diag f8 
      call reorder   ("tpf ","s1  ",1432,4)
      call mult      ("t11 ","s1  ","s2  ",131)
      call unfold_ccsd("s2  ","s1  ",4)
      call add       ("s1  ",-1.d0,"vv0 ",1.d0,"vv0 ")
c  diag f10a
      call reorder   ("tqf ","s1  ",1324,4)
      call reorder   ("t1g ","s2  ",1243,4)
      call mult      ("s1  ","s2  ","s3  ",311)
      call reorder   ("s3  ","s1  ",1243,4)
      call mult      ("s1  ","t11 ","s2  ",311)
      call add       ("s2  ",1.d0,"tb4 ",1.d0,"tb4 ")
c  diag f11 
      call reorder   ("tqf ","s1  ",1342,4)
      call mult      ("s1  ","t11 ","s2  ",311)
      call unfold_ccsd("s2  ","s1  ",3)
      call add       ("s1  ",-1.d0,"tb4 ",1.d0,"tb4 ")
c  dgt1p4
      call mult      ("t1g ","pphh","s1  ",131)
      call mult      ("t1g ","s1  ","s2  ",131)
      call reorder   ("t1c ","s1  ",1243,4)
      call mult      ("s1  ","s2  ","s3  ",131)
      call mult      ("t1v ","s3  ","s1  ",131)
      call reorder   ("s1  ","s2  ",3412,4)
      call add       ("s2  ",-2.d0,"vv0 ",1.d0,"vv0 ")
      if (iopt(17).ne.0) then
c  diags 30a,37a
        call reord6    ("t3p ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125634)
        call mult6     ("boxB","s1  ","s2  ",2)
        call add       ("s2  ",1.d0,"tb3 ",1.d0,"tb3 ")
c  diags 30b,37b
        call reord6    ("t3q ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125634)
        call mult6     ("s1  ","boxB","s2  ",2)
        call add       ("s2  ",1.d0,"tb4 ",1.d0,"tb4 ")
c  diags 31a,36a 
        call reord6    ("t3u ","s1  ",125643)
        call mult6     ("s1  ","boxC","s2  ",2)
        call add       ("s2  ",1.d0,"tb3 ",1.d0,"tb3 ")
c  diags 32a,39a
        call reord6    ("t3p ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126534)
        call mult6     ("boxF","s1  ","s2  ",2)
        call add       ("s2  ",1.d0,"vv0 ",1.d0,"vv0 ")
c  diags 32b,39b
        call reord6    ("t3p ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126543)
        call mult6     ("s1  ","boxF","s2  ",2)
        call add       ("s2  ",1.d0,"tb4 ",1.d0,"tb4 ")
c  diags 33a,38a
        call reord6    ("t3u ","s1  ",126543)
        call mult6     ("s1  ","boxG","s2  ",2)
        call add       ("s2  ",1.d0,"tb3 ",1.d0,"tb3 ")
      end if
      call prtime    ("calvv0  ",0,0)
      write (6,*)    ' '
      return
      end
C  *********************************************************************
      subroutine ordtb0_cht
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder    ("tb1 ","s1  ",1324,4)
      call add        ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
      call reorder    ("tb2 ","s1  ",1342,4)
      call add        ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
      call reorder    ("tb3 ","s1  ",1243,4)
      call add        ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
      call reorder    ("tb4 ","s1  ",1423,4)
      call add        ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
      call prtime   ("ordtb   ",0,1)
      return
      end
C *********************************************************************
      subroutine caltv2_cht
!---------------------------------------------------------------------
!
!
! add description
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use cij1_mod
      use gen_mod
      use dat_mod
      use tv2_mod
#include "implicit.h"
#include "priunit.h"
#include "strt.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     parameter (ndij=960000,nbig=60000000)
#include "param_ccsd.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1iiia(200),iiva(200)
#include "cfiles.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!    1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!    2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
#include "strt.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
#include "dat.h"
!     common /tv2/ a(500,500)
#include "tv2.h"
#endif
      external find2,find4

!!!!!!!!!!!!! miro: dynamic allocations needed !!!!!!!!!!!
      dimension b(500),indlst(500),ipvt(500),z(500)

      call move ("vv0 ","trhs")
      if (iopt(17).ne.0) then
        call tvsub2
      end if
      do 10 i=1,500
        b(i)=0.d0
        do 20 j=1,500
          a(i,j)=0.d0
20      continue
10    continue
      call files ("ggvh",ntape,new,1,-1,-1)
      istrt1=nstrt(1,ntape-1)+2*nij1+4
      nequn=ije1(nij1)
      do 30 i=1,nequn
        indlst(i)=igen(istrt1+i)
30    continue
      call files ("trhs",ntape,new,1,-1,-1)
      do 60 i=1,nequn
        b(i)=-1.d0*gen(nstrt(2,ntape-1)+i)
60    continue
c  diag 2a
      call move       ("hh  ","boxX")
c  diag 14a
      call reorder    ("pphh","s1  ",3412,4)
      call reorder    ("pphh","s2  ",3421,4)
      call add        ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call reorder    ("tcm ","s2  ",3412,4)
      call mult       ("s1  ","s2  ","s3  ",113)
      call add        ("s3  ",1.d0,"boxX",1.d0,"boxX")
c  diag 19b
      call reorder    ("phhh","s1  ",2413,4)
      call reorder    ("phhh","s2  ",3412,4)
      call add        ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call mult       ("s1  ","t1c ","s2  ",202)
      call add        ("s2  ",1.d0,"boxX",1.d0,"boxX")
c  diag 29b
      call reorder    ("ph  ","s1  ",1243,4)
      call reorder    ("t1c ","s2  ",1243,4)
      call mult       ("s1  ","s2  ","s3  ",111)
      call add        ("s3  ",-1.d0,"boxX",1.d0,"boxX")
      do 40 i=1,nequn 
        idi=indlst(i)
        do 50 j=1,nequn
          idj=indlst(j)
c   diags boxX
          call files ("boxX",ntape,new,1,-1,-1)
          nsize1=ije1(nij1)
          istrt11=nstrt(1,ntape-1)+2*nij1+4
          istrt12=nstrt(2,ntape-1)
          if (iand(idi,z'ffffff00').eq.iand(idj,z'ffffff00')) then
            ii1=iand(idj,z'ff')
            ii2=iand(idi,z'ff')
            a(i,j)=a(i,j)+find2(ii1,ii2)
          end if
c   diag f1
          call files ("tgf ",ntape,new,1,-1,-1)
          nsize1=ije1(nij1)
          istrt11=nstrt(1,ntape-1)+2*nij1+4
          istrt12=nstrt(2,ntape-1)
          if (iand(idi,z'ff00ffff').eq.iand(idj,z'ff00ffff')) then
            ii1=iand(ishft(idi,-16),z'ff')
            ii2=iand(ishft(idj,-16),z'ff')
            a(i,j)=a(i,j)+find2(ii1,ii2)
          end if
c   diag f2
          if (iand(idi,z'00ffffff').eq.iand(idj,z'00ffffff')) then
            ii1=iand(ishft(idi,-24),z'ff')
            ii2=iand(ishft(idj,-24),z'ff')
            a(i,j)=a(i,j)+find2(ii1,ii2)
          end if
c   diag f3
          call files ("tvf ",ntape,new,1,-1,-1)
          nsize1=ije1(nij1)
          istrt11=nstrt(1,ntape-1)+2*nij1+4
          istrt12=nstrt(2,ntape-1)
          if (iand(idi,z'ffff00ff').eq.iand(idj,z'ffff00ff')) then
            ii1=iand(ishft(idi,-8),z'ff')
            ii2=iand(ishft(idj,-8),z'ff')
            a(i,j)=a(i,j)-find2(ii1,ii2)
          end if
c   diag f4
          call files ("thuf",ntape,new,1,-1,-1)
          nsize2=ije1(nij1)
          istrt21=nstrt(1,ntape-1)+2*nij1+4
          istrt22=nstrt(2,ntape-1)
          if (iand(idi,z'0000ffff').eq.iand(idj,z'0000ffff')) then
            ii1=iand(ishft(idi,-24),z'ff')
            ii2=iand(ishft(idi,-16),z'ff')
            ii3=iand(ishft(idj,-24),z'ff')
            ii4=iand(ishft(idj,-16),z'ff')
            a(i,j)=a(i,j)-find4(ii1,ii2,ii3,ii4)
          end if
c   diag f6
          call files ("tpf ",ntape,new,1,-1,-1)
          nsize2=ije1(nij1)
          istrt21=nstrt(1,ntape-1)+2*nij1+4
          istrt22=nstrt(2,ntape-1)
          if (iand(idi,z'ff0000ff').eq.iand(idj,z'ff0000ff')) then
            ii1=iand(ishft(idi,-16),z'ff')
            ii2=iand(ishft(idj,-8),z'ff')
            ii3=iand(ishft(idj,-16),z'ff')
            ii4=iand(ishft(idi,-8),z'ff')
            a(i,j)=a(i,j)+find4(ii1,ii2,ii3,ii4)
          end if
c   diag f7
          if (iand(idi,z'00ff00ff').eq.iand(idj,z'00ff00ff')) then
            ii1=iand(ishft(idi,-24),  z'ff')
            ii2=iand(ishft(idj,-8),   z'ff')
            ii3=iand(ishft(idj,-24),  z'ff')
            ii4=iand(ishft(idi,-8),   z'ff')
            a(i,j)=a(i,j)+find4(ii1,ii2,ii3,ii4)
          end if
c   diag f9
          call files ("tqf ",ntape,new,1,-1,-1)
          nsize2=ije1(nij1)
          istrt21=nstrt(1,ntape-1)+2*nij1+4
          istrt22=nstrt(2,ntape-1)
          if (iand(idi,z'00ff00ff').eq.iand(idj,z'00ff00ff')) then
            ii1=iand(ishft(idi,-24),  z'ff')
            ii2=iand(ishft(idj,-8),   z'ff')
            ii3=iand(ishft(idi,-8),   z'ff')
            ii4=iand(ishft(idj,-24),  z'ff')
            a(i,j)=a(i,j)-2.d0*find4(ii1,ii2,ii3,ii4)
          end if
c   diag f10
          if ((iand(idi,z'ff').eq.iand(idj,z'ff')).and.  
     1    (iand(ishft(idi,-16), z'ff').eq.
     2    iand(ishft(idj,-24), z'ff'))) then
            ii1=iand(ishft(idi,-24),  z'ff')
            ii2=iand(ishft(idj,-8),   z'ff')
            ii3=iand(ishft(idi,-8),   z'ff')
            ii4=iand(ishft(idj,-16),  z'ff')
            a(i,j)=a(i,j)+find4(ii1,ii2,ii3,ii4)
          end if
          if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
            call tvsub1 (idi,idj,i,j)
          end if
50      continue
40    continue
      call clear ("s1  ",2232,1110)
      call files ("s1  ",ntape,new,1,-1,-1)
      nsize=ije1(nij1)
      istrt=nstrt(2,ntape-1)
      call dcopy (nsize,-1.d0,0,gen(istrt+1),1)
      call diveps ("s1  ")
      call cmult  ("s1  ",-1.d0,"s1  ")
      do 80 i=1,nequn
        a(i,i)=a(i,i)+1.d0/gen(istrt+i)
80    continue
      call dgeco (a,500,nequn,ipvt,rcond,z)
      if (abs(rcond).lt.1.d-16) then
        write (6,*) 'Linear equations for t12 are singular.'
        stop
      end if
      call dgesl (a,500,nequn,ipvt,b,0)
      if (iopt(17).ne.3) then
        call clear ("t12 ",2232,1110)
        call files ("t12 ",ntape,new,1,-1,-1)
      else if (iopt(17).eq.3) then
        call clear ("tnew",2232,1110)
        call files ("tnew",ntape,new,1,-1,-1)
      end if
      nsize=ije1(nij1)
      istrt=nstrt(2,ntape-1)
      do 70 i=1,nsize
        gen(istrt+i)=b(i)
70    continue
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
c  diag 2a
        call reorder  ("t12 ","s1  ",2134,4)
        call mult6    ("s1  ","phhh","s2  ",1)
        call reord6   ("s2  ","s1  ",124536)
        call zero     ("s1  ")
        call dveps6   ("s1  ")
        call cmult    ("s1  ",-1.d0,"s1  ")
        call add      ("s1  ",-1.d0,"t12q",1.d0,"t12q")
      end if
      call prtime ("caltv2  ",0,0)
      return
      end
C  *********************************************************************
      subroutine tinit_cht
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call move    ("vv0 ","t12 ")
      call diveps  ("t12 ")
      call cmult   ("t12 ",-1.d0,"t12 ")
      call t1sqr   ("t11 ","t1g ","t12 ","tm  ")
      call move    ("v30p","t12p")
      call reord6  ("t12p","s1  ",123546)
      call reord7  ("s1  ","s2  ",213456)
      call add     ("t12p",1.d0,"s2  ",1.d0,"t12p")
      call zero    ("t12p")
      call dveps6  ("t12p")
      call cmult   ("t12p",-1.d0,"t12p")
      call move    ("v30q","t12q")
      call zero    ("t12q")
      call dveps6  ("t12q")
      call cmult   ("t12q",-1.d0,"t12q")
      return
      end 
C  ********************************************************************
      subroutine calt12p
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
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
c  dg3n1g
        call reorder   ("ggph","s1  ",2143,4)
        call mult6     ("s1  ","tv  ","s3  ",1)
c  dg3n2c
        call mult6     ("thu ","vhph","s1  ",1)
        call add       ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
c  dg3nf1
        call reorder   ("tp  ","s1  ",4321,4)
        call mult6     ("gggg","s1  ","s2  ",1)
        call cmult     ("s2  ",-1.d0,"s1  ")
        call fold6     ("s3  ","s1  ",idum,2)
c
        call reord6    ("s3  ","t12p",126354)
c
c  dg3nf2
        call reorder   ("thu ","s1  ",1432,4)
        call reorder   ("gvgv","s2  ",4123,4)
        call mult6     ("s1  ","s2  ","s3  ",1)
        call cmult     ("s3  ",-1.d0,"s3  ") 
c  dg3n1c
        call reorder   ("tg  ","s1  ",1432,4)
        call mult6     ("s1  ","vgpp","s2  ",1)
        call fold6     ("s2  ","s3  ",idum,2)
c  dg3n2g
        call reorder   ("tp  ","s1  ",4123,4)
        call mult6     ("ghhh","s1  ","s3  ",1)
        call add       ("s3  ",-1.d0,"s2  ",1.d0,"s3  ")
c
        call reord6    ("s3  ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126435)
        call add       ("s1  ",1.d0,"t12p",1.d0,"t12p")
c
c  dg3n1e
        call mult6     ("tp3 ","gphp","s3  ",1)
c  dg3n2e
        call mult6     ("gvhh","tg  ","s1  ",1)
        call add       ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
c  dg3nf3
        call reorder   ("gvgv","s1  ",1432,4)
        call mult6     ("s1  ","tp  ","s2  ",1)
        call fold6     ("s3  ","s2  ",idum,2)
c
        call reord6    ("s3  ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125463)
        call add       ("s1  ",1.d0,"t12p",1.d0,"t12p")
c
        call reord6    ("t12p","s1  ",123546)
        call reord7    ("s1  ","s2  ",213456)
        call add       ("t12p",1.d0,"s2  ",1.d0,"t12p")
        call zero      ("t12p")
        call dveps6    ("t12p")
        call cmult     ("t12p",-1.d0,"t12p")
      else if (iopt(17).eq.3) then
        call move      ("v30p","tpnw")
c  GROUP V
c   diag 1
        call reord6    ("t12p","s1  ",124653)
        call mult6     ("s1  ","bx5 ","s2  ",3)
        call reord6    ("s2  ","s1  ",126354)
        call add       ("s1  ",1.d0,"tpnw",1.d0,"tpnw")
c  GROUP VI
c   diag 2
        call reord6    ("t12p","s1  ",123645)
        call mult6     ("s1  ","bx6 ","s2  ",3)
        call reord6    ("s2  ","s1  ",123564)
        call add       ("s1  ",1.d0,"tpnw",1.d0,"tpnw")
c  GROUP VIII
c   diag 1
        call reord6    ("t12p","s1  ",124563)
        call mult6     ("s1  ","bx7 ","s2  ",5)
        call reord6    ("s2  ","s1  ",126345)
        call add       ("s1  ",0.5d0,"tpnw",1.d0,"tpnw")
c  GROUP IX
c   diag 2
        call reord6    ("t12p","s1  ",123465)
        call mult6     ("s1  ","bx8 ","s2  ",5)
        call reord6    ("s2  ","s1  ",123465)
        call add       ("s1  ",1.d0,"tpnw",1.d0,"tpnw")
c
        call fold6     ("tpnw","s1  ",3,1)
        call fold6     ("s1  ","s2  ",4,1)
        call fold6     ("s2  ","ttpf",5,1)
        call reord6    ("ttpf","s1  ",123546)
        call reord7    ("s1  ","s2  ",213456)
        call add       ("ttpf",1.d0,"s2  ",1.d0,"ttpf")
        call zero      ("tpnw")
      end if
      call prtime    ("calt12p ",0,0)
      return
      end
C  *******************************************************************
      subroutine calt12q
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
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
c  dg3nf4
        call reorder   ("tq  ","s1  ",1423,4)
        call reorder   ("gvvg","s2  ",1432,4)
        call mult6     ("s1  ","s2  ","s3  ",1)
c  dg3n1a
        call reorder   ("ghpp","s1  ",1243,4)
        call mult6     ("s1  ","tq3 ","s2  ",1)
        call fold6     ("s2  ","s3  ",idum,2)
c  dg3n2i
        call reorder   ("tg  ","s1  ",1423,4)
        call mult6     ("s1  ","ghvh","s3  ",1)
        call add       ("s3  ",-1.d0,"s2  ",1.d0,"s3  ")
c
        call reord6    ("s3  ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","t12q",124563)
c
c  dg3nf7
        call reorder   ("gvgv","s1  ",1432,4)
        call reorder   ("tq  ","s2  ",1243,4)
        call mult6     ("s1  ","s2  ","s3  ",1)
c  dg3n1b
        call reorder   ("ghpp","s1  ",1423,4)
        call mult6     ("tp3 ","s1  ","s2  ",1)
        call fold6     ("s2  ","s3  ",idum,2)
c  dg3n2h
        call mult6     ("gvhh","tg6 ","s1  ",1)
        call add       ("s1  ",-1.d0,"s2  ",1.d0,"s2  ")
c
        call reord6    ("s2  ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125436)
        call add       ("s1  ",1.d0,"t12q",1.d0,"t12q")
c
c  dg3n1d
        call reorder   ("tq  ","s1  ",1432,4)
        call mult6     ("gphp","s1  ","s3  ",1)
c  dg3n2f
        call mult6     ("tg  ","ghvh","s1  ",1)
        call add       ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
c  dg3nf8
        call reorder   ("gvvg","s1  ",1432,4)
        call mult6     ("tp  ","s1  ","s2  ",1)
        call fold6     ("s3  ","s2  ",idum,2)
c
        call reord6    ("s3  ","s1  ",124365)
        call reord7    ("s1  ","s2  ",132456)
        call add       ("s2  ",1.d0,"t12q",1.d0,"t12q")
c
c  dg3n1f
        call reorder   ("tg  ","s1  ",1432,4)
        call mult6     ("s1  ","gpvp","s3  ",1)
c  dg3n2d
        call mult6     ("ghhh","tq  ","s1  ",1)
        call add       ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
c  dg3nf6
        call reorder   ("thu ","s1  ",1432,4)
        call mult6     ("s1  ","gvvg","s2  ",1)
        call cmult     ("s2  ",-1.d0,"s1  ")
        call fold6     ("s3  ","s1  ",idum,2)
c
        call reord6    ("s3  ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125463)
        call add       ("s1  ",1.d0,"t12q",1.d0,"t12q")
c
c  dg3n1h
        call reorder   ("tv  ","s1  ",1243,4)
        call mult6     ("gghp","s1  ","s3  ",1)
c  dg3n2b
        call reorder   ("vhhp","s1  ",1243,4)
        call mult6     ("thu ","s1  ","s2  ",1)
        call add       ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
c  dg3nf5
        call reorder   ("tq  ","s1  ",3421,4)
        call mult6     ("gggg","s1  ","s2  ",1)
        call cmult     ("s2  ",-1.d0,"s1  ")
        call fold6     ("s3  ","s1  ",idum,2)
c
        call reord6    ("s3  ","s1  ",126345)
        call add       ("s1  ",1.d0,"t12q",1.d0,"t12q")
c
c  dg3n1i
        call reorder   ("ggvp","s1  ",2134,4)
        call reorder   ("tc  ","s2  ",1432,4)
        call mult6     ("s1  ","s2  ","s3  ",1)
        call reord6    ("s3  ","s1  ",124536)
        call add       ("s1  ",1.d0,"t12q",1.d0,"t12q")
c
        call zero      ("t12q")
        call dveps6    ("t12q")
        call cmult     ("t12q",-1.d0,"t12q")
      else if (iopt(17).eq.3) then
        call move       ("v30q","tqnw")
c  GROUP I
c   diag 6
        call box38     (1)
        call reorder   ("tc  ","s1  ",1432,4)
        call mult6     ("bx38","s1  ","s2  ",1)
        call reord6    ("s2  ","s1  ",124536)
        call add       ("s1  ",1.d0,"tqnw",1.d0,"tqnw")
c  GROUP II
c   diag 1
        call reorder   ("t12 ","s1  ",2134,4)
        call mult6     ("s1  ","bx2 ","s2  ",1)
        call reord6    ("s2  ","s1  ",126534)
        call add       ("s1  ",1.d0,"tqnw",1.d0,"tqnw")
c  GROUP III
c   diag 1
        call reord6    ("t12q","s1  ",124563)
        call mult6     ("s1  ","bx3 ","s2  ",3)
        call reord6    ("s2  ","s1  ",126345)
        call add       ("s1  ",1.d0,"tqnw",1.d0,"tqnw")
c  GROUP IV 
c   diag 1
        call reord6    ("t12q","s1  ",125643)
        call mult6     ("s1  ","bx4 ","s2  ",3)
        call reord6    ("s2  ","s1  ",126435)
        call add       ("s1  ",1.d0,"tqnw",1.d0,"tqnw")
c   diag 3
        call reord6    ("t12p","s1  ",124653)
        call mult6     ("s1  ","bx4 ","s2  ",3)
        call reord6    ("s2  ","s1  ",126345)
        call add       ("s1  ",1.d0,"tqnw",1.d0,"tqnw")
c  GROUP V 
c   diag 1
        call reord6    ("t12q","s1  ",125643)
        call mult6     ("s1  ","bx5 ","s2  ",3)
        call reord6    ("s2  ","s1  ",126534)
        call add       ("s1  ",1.d0,"tqnw",1.d0,"tqnw")
c  GROUP VI
c   diag 1
        call reord6    ("t12q","s1  ",123546)
        call mult6     ("s1  ","bx6a","s2  ",3)
        call reord6    ("s2  ","s1  ",123546)
        call add       ("s1  ",1.d0,"tqnw",1.d0,"tqnw")
c  GROUP VIII
c   diag 1
        call reord6    ("t12q","s1  ",124563)
        call mult6     ("s1  ","bx7 ","s2  ",5)
        call reord6    ("s2  ","s1  ",126345)
        call add       ("s1  ",1.d0,"tqnw",1.d0,"tqnw")
c  GROUP IX  
c   diag 1
        call mult6     ("t12q","bx8 ","s1  ",5)
        call add       ("s1  ",1.d0,"tqnw",1.d0,"tqnw")
c   diag 2
        call reord6    ("t12q","s1  ",123564)
        call mult6     ("s1  ","bx8 ","s2  ",5)
        call reord6    ("s2  ","s1  ",123645)
        call add       ("s1  ",1.d0,"tqnw",1.d0,"tqnw")
c
        call fold6     ("tqnw","s1  ",3,1)
        call fold6     ("s1  ","s2  ",4,1)
        call fold6     ("s2  ","ttqf",6,1)
        call zero      ("tqnw")
      end if
      call prtime    ("calt12q ",0,0)
      return
      end
C  ********************************************************************
      subroutine dfoldp
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

c  diag 1
      call mult6      ("t12p","tvf ","s1  ",5)
      call add        ("s1  ",-0.5d0,"tpnw",1.d0,"tpnw")
c  diag 2
      call reord7     ("t12p","s1  ",132456)
      call reord6     ("s1  ","s2  ",124563)
      call mult6      ("s2  ","tgf ","s1  ",5)
      call reord6     ("s1  ","s2  ",126345)
      call reord7     ("s2  ","s1  ",132456)
      call add        ("s1  ",1.d0,"tpnw",1.d0,"tpnw")
c  diag 3
      call reorder    ("tp  ","s1  ",2341,4)
      call mult6      ("thuf","s1  ","s2  ",1)
      call reord6     ("s2  ","s1  ",124356)
      call cmult      ("s1  ",-1.d0,"s1  ")   
      call fold6      ("tpnw","s1  ",idum,2) 
c  diag 4
      call t1t2       ("t1g ","tp  ","t12p","s1  ",123456,1)
      call t1t2       ("t1g ","tp  ","s1  ","s2  ",123456,2)
      call reord7     ("s2  ","s1  ",312456)
      call reord6     ("s1  ","s2  ",124356)
      call reord7     ("s2  ","s1  ",132456)
      call reord6     ("s1  ","s2  ",125634)
      call mult6      ("thuf","s2  ","s1  ",3)
      call add        ("s1  ",-0.5d0,"tpnw",1.d0,"tpnw")
c  diag 5
      call reorder    ("thu ","s1  ",1432,4)
      call reorder    ("tpf ","s2  ",1243,4)
      call mult6      ("s1  ","s2  ","s3  ",1)
      call reord6     ("s3  ","s1  ",124356)
      call reord7     ("s1  ","s2  ",132456)
      call reord6     ("s2  ","s1  ",125436)
      call cmult      ("s1  ",-1.d0,"s1  ")   
      call fold6      ("tpnw","s1  ",idum,2) 
c  diag 6
      call reorder    ("tpf ","s1  ",1432,4)
      call mult6      ("tp  ","s1  ","s2  ",1)
      call reord6     ("s2  ","s1  ",124365)
      call reord7     ("s1  ","s2  ",132456)
      call fold6      ("tpnw","s2  ",idum,2) 
c  diag 7
      call t1t2       ("t1g ","tp  ","t12p","s1  ",123456,2)
      call reorder    ("t1v ","s2  ",1243,4)
      call t1t2       ("s2  ","thu ","s1  ","s3  ",123456,3)
      call reord7     ("s3  ","s1  ",132456)
      call reord6     ("s1  ","s2  ",124536)
      call reorder    ("tpf ","s1  ",1432,4)
      call mult6      ("s2  ","s1  ","s3  ",3)
      call reord6     ("s3  ","s1  ",125346)
      call reord7     ("s1  ","s2  ",132456)
      call add        ("s2  ",1.d0,"tpnw",1.d0,"tpnw")
c  diag 8
      call reord6     ("ttpf","s1  ",124563)
      call reorder    ("t1v ","s2  ",1243,4)
      call mult6      ("s1  ","s2  ","s3  ",5)
      call reord6     ("s3  ","s1  ",126345)
      call cmult      ("s1  ",0.5d0,"s1  ")   
      call fold6      ("tpnw","s1  ",idum,2) 
c  diag 9
      call reord6     ("ttpf","s1  ",123465)
      call reorder    ("t1g ","s2  ",1243,4)
      call mult6      ("s1  ","s2  ","s3  ",5)
      call reord6     ("s3  ","s1  ",123465)
      call cmult      ("s1  ",-1.d0,"s1  ")   
      call fold6      ("tpnw","s1  ",idum,2) 
c  diag 10
      call reord6     ("ttpf","s1  ",123645)
      call reorder    ("thum","s2  ",3412,4)
      call mult6      ("s1  ","s2  ","s3  ",3)
      call reord6     ("s3  ","s1  ",123564)
      call cmult      ("s1  ",0.5d0,"s1  ")   
      call fold6      ("tpnw","s1  ",idum,2) 
c  diag 11
      call reord6     ("ttpf","s1  ",124653)
      call reorder    ("tpm ","s2  ",2314,4)
      call mult6      ("s1  ","s2  ","s3  ",3)
      call reord6     ("s3  ","s1  ",125364)
      call cmult      ("s1  ",-1.d0,"s1  ")   
      call fold6      ("tpnw","s1  ",idum,2) 
c  diag 12
      call cmult      ("t1g ",2.d0,"s1  ")
      call t1t2       ("s1  ","tp  ","t12p","s2  ",123456,1)
      call reorder    ("t1v ","s1  ",1243,4)
      call t1t2       ("s1  ","thum","s2  ","s3  ",123456,3)
      call reord7     ("s3  ","s1  ",312456)
      call reord6     ("s1  ","s2  ",124356)
      call reord7     ("s2  ","s1  ",132456)
      call reord6     ("s1  ","s2  ",125346)
      call reord6     ("ttpf","s1  ",126453)
      call mult6      ("s1  ","s2  ","s3  ",4)
      call reord6     ("s3  ","s1  ",124563)
      call add        ("s1  ",-0.5d0,"tpnw",1.d0,"tpnw")
c
      call reord6     ("tpnw","s1  ",123546)
      call reord7     ("s1  ","s2  ",213456)
      call add        ("tpnw",1.d0,"s2  ",1.d0,"tpnw")
      call dveps6     ("tpnw")
      call cmult      ("tpnw",-1.d0,"tpnw")
      call prtime     ("dfoldp  ",0,0)
      return
      end
C  ********************************************************************
      subroutine dfoldq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  diag 1
      call reord6     ("t12q","s1  ",123465)
      call mult6      ("s1  ","tvf ","s2  ",5)
      call reord6     ("s2  ","s1  ",123465)
      call add        ("s1  ",-1.d0,"tqnw",1.d0,"tqnw")
c  diag 2
      call reord7     ("t12q","s1  ",132456)
      call reord6     ("s1  ","s2  ",124563)
      call mult6      ("s2  ","tgf ","s1  ",5)
      call reord6     ("s1  ","s2  ",126345)
      call reord7     ("s2  ","s1  ",132456)
      call add        ("s1  ",1.d0,"tqnw",1.d0,"tqnw")
c  diag 3
      call reord7     ("t12q","s1  ",231456)
      call reord6     ("s1  ","s2  ",124563)
      call mult6      ("s2  ","tgf ","s1  ",5)
      call reord6     ("s1  ","s2  ",126345)
      call reord7     ("s2  ","s1  ",312456)
      call add        ("s1  ",1.d0,"tqnw",1.d0,"tqnw")
c  diag 4
      call reorder    ("tq  ","s1  ",2341,4)
      call mult6      ("thuf","s1  ","s2  ",1)
      call reord6     ("s2  ","s1  ",124356)
      call cmult      ("s1  ",-1.d0,"s1  ")
      call fold6      ("tqnw","s1  ",idum,2)
c  diag 5
      call t1t2       ("t1g ","tq  ","t12q","s1  ",123456,1)
      call t1t2       ("t11 ","tg  ","s1  ","s2  ",123456,2)
      call reord7     ("s2  ","s1  ",312456)
      call reord6     ("s1  ","s2  ",124356)
      call reord7     ("s2  ","s1  ",132456)
      call reord6     ("s1  ","s2  ",125634)
      call mult6      ("thuf","s2  ","s1  ",3)
      call add        ("s1  ",-1.d0,"tqnw",1.d0,"tqnw")
c  diag 6
      call reorder    ("tq  ","s1  ",1243,4)
      call reorder    ("tqf ","s2  ",1432,4)
      call mult6      ("s1  ","s2  ","s3  ",1)
      call reord6     ("s3  ","s1  ",124563)
      call reord7     ("s1  ","s2  ",132456)
      call fold6      ("tqnw","s2  ",idum,2)
c  diag 7
      call reorder    ("tqf ","s1  ",1432,4)
      call mult6      ("tp  ","s1  ","s2  ",1)
      call reord6     ("s2  ","s1  ",124365)
      call reord7     ("s1  ","s2  ",132456)
      call fold6      ("tqnw","s2  ",idum,2)
c  diag 8
      call reorder    ("thu ","s1  ",1432,4)
      call mult6      ("s1  ","tqf ","s2  ",1)
      call reord6     ("s2  ","s1  ",124356)
      call reord7     ("s1  ","s2  ",132456)
      call reord6     ("s2  ","s1  ",125463)
      call cmult      ("s1  ",-1.d0,"s1  ")
      call fold6      ("tqnw","s1  ",idum,2)
c  diag 9
      call reord7     ("t12q","s1  ",132456)
      call reord6     ("s1  ","s2  ",124635)
      call reorder    ("tqf ","s1  ",1342,4)
      call mult6      ("s2  ","s1  ","s3  ",3)
      call reord6     ("s3  ","s1  ",125364)
      call reord7     ("s1  ","s2  ",132456)
      call add        ("s2  ",-2.d0,"tqnw",1.d0,"tqnw")
c  diag 10
      call t1t2       ("t1g ","tq  ","t12q","s1  ",123456,1)
      call t1t2       ("t11 ","tg  ","s1  ","s2  ",123456,2)
      call reord7     ("s2  ","s1  ",231456)
      call reord6     ("s1  ","s2  ",124653)
      call reorder    ("tqf ","s1  ",1324,4)
      call mult6      ("s2  ","s1  ","s3  ",3)
      call reord6     ("s3  ","s1  ",125364)
      call reord7     ("s1  ","s2  ",132456)
      call add        ("s2  ",1.d0,"tqnw",1.d0,"tqnw")
c  diag 11
      call reorder    ("t1v ","s1  ",1243,4)
      call t1t2       ("s1  ","thu ","t12p","s2  ",123456,3)
      call t1t2       ("t1g ","tp  ","s2  ","s1  ",123456,2)
      call reord7     ("s1  ","s2  ",132456)
      call reord6     ("s2  ","s1  ",124563)
      call reorder    ("tqf ","s2  ",1324,4)
      call mult6      ("s1  ","s2  ","s3  ",3)
      call reord6     ("s3  ","s1  ",125364)
      call reord7     ("s1  ","s2  ",132456)
      call add        ("s2  ",1.d0,"tqnw",1.d0,"tqnw")
c  diag 12
      call reorder    ("tpf ","s1  ",1342,4)
      call reorder    ("tq  ","s2  ",1243,4)
      call mult6      ("s1  ","s2  ","s3  ",1)
      call reord6     ("s3  ","s1  ",124356)
      call reord7     ("s1  ","s2  ",132456)
      call reord6     ("s2  ","s1  ",125346)
      call fold6      ("tqnw","s1  ",idum,2)
c  diag 13
      call reord7     ("t12q","s1  ",132456)
      call reord6     ("s1  ","s2  ",124653)
      call reorder    ("tpf ","s1  ",1423,4)
      call mult6      ("s2  ","s1  ","s3  ",3)
      call reord6     ("s3  ","s1  ",125364)
      call reord7     ("s1  ","s2  ",132456)
      call add        ("s2  ",1.d0,"tqnw",1.d0,"tqnw")
c  diag 14
      call t1t2       ("t1g ","tq  ","t12q","s1  ",123456,1)
      call t1t2       ("t11 ","tg  ","s1  ","s2  ",123456,2)
      call reord7     ("s2  ","s1  ",231456)
      call reord6     ("s1  ","s2  ",124653)
      call reorder    ("tpf ","s1  ",1423,4)
      call mult6      ("s1  ","s2  ","s3  ",3)
      call reord7     ("s3  ","s1  ",132456)
      call reord6     ("s1  ","s2  ",124536)
      call add        ("s2  ",1.d0,"tqnw",1.d0,"tqnw")
c  diag 15
      call reord6     ("ttqf","s1  ",123564)
      call reorder    ("t1g ","s2  ",1243,4)
      call mult6      ("s1  ","s2  ","s3  ",5)
      call reord6     ("s3  ","s1  ",123645)
      call cmult      ("s1  ",-1.d0,"s1  ")
      call fold6      ("tqnw","s1  ",idum,2)
c  diag 16
      call reorder    ("t1g ","s1  ",1243,4)
      call mult6      ("ttqf","s1  ","s2  ",5)
      call cmult      ("s2  ",-1.d0,"s1  ")
      call fold6      ("tqnw","s1  ",idum,2)
c  diag 17
      call reord6     ("ttqf","s1  ",124563)
      call reorder    ("t1v ","s2  ",1243,4)
      call mult6      ("s1  ","s2  ","s3  ",5)
      call reord6     ("s3  ","s1  ",126345)
      call fold6      ("tqnw","s1  ",idum,2)
c  diag 18
      call reord6     ("ttqf","s1  ",123546)
      call reorder    ("thum","s2  ",3412,4)
      call mult6      ("s1  ","s2  ","s3  ",3)
      call reord6     ("s3  ","s1  ",123546)
      call fold6      ("tqnw","s1  ",idum,2)
c  diag 19
      call reord6     ("ttqf","s1  ",125634)
      call reorder    ("tpm ","s2  ",2341,4)
      call mult6      ("s1  ","s2  ","s3  ",3)
      call reord6     ("s3  ","s1  ",125634)
      call cmult      ("s1  ",-1.d0,"s1  ")
      call fold6      ("tqnw","s1  ",idum,2)
c  diag 20
      call reord6     ("ttqf","s1  ",125643)
      call reorder    ("tq  ","s2  ",2413,4)
      call mult6      ("s1  ","s2  ","s3  ",3)
      call reord6     ("s3  ","s1  ",125436)
      call cmult      ("s1  ",-1.d0,"s1  ")
      call fold6      ("tqnw","s1  ",idum,2)
c  diag 21
      call reord6     ("ttqf","s1  ",124536)
      call reorder    ("tq  ","s2  ",2431,4)
      call mult6      ("s1  ","s2  ","s3  ",3)
      call reord6     ("s3  ","s1  ",125346)
      call cmult      ("s1  ",2.d0,"s1  ")
      call fold6      ("tqnw","s1  ",idum,2)
c  diag 22
      call reord6     ("ttqf","s1  ",124536)
      call reorder    ("tpm ","s2  ",2341,4)
      call mult6      ("s1  ","s2  ","s3  ",3)
      call reord6     ("s3  ","s1  ",125346)
      call cmult      ("s1  ",-1.d0,"s1  ")
      call fold6      ("tqnw","s1  ",idum,2)
c  diag 23
      call reorder    ("t1v ","s1  ",1243,4)
      call t1t2       ("s1  ","thu ","t12p","s2  ",123456,3)
      call t1t2       ("t1g ","tp  ","s2  ","s1  ",123456,2)
      call t1t2       ("t1g ","tpm ","s1  ","s2  ",123456,1)
      call reord7     ("s2  ","s1  ",312456)
      call reord6     ("s1  ","s2  ",124356)
      call reord7     ("s2  ","s1  ",132456)
      call reord6     ("s1  ","s2  ",125346)
      call reord6     ("ttqf","s1  ",125463)
      call mult6      ("s1  ","s2  ","s3  ",4)
      call reord6     ("s3  ","s1  ",124536)
      call add        ("s1  ",-1.d0,"tqnw",1.d0,"tqnw")
c  diag 24
      call t1t2       ("t1g ","tq  ","t12q","s1  ",123456,1)
      call t1t2       ("t11 ","tg  ","s1  ","s2  ",123456,2)
      call reord7     ("s2  ","s1  ",312456)
      call reord6     ("s1  ","s2  ",124356)
      call reord7     ("s2  ","s1  ",132456)
      call reord6     ("s1  ","s2  ",126345)
      call reord6     ("ttqf","s1  ",125463)
      call mult6      ("s1  ","s2  ","s3  ",4)
      call reord6     ("s3  ","s1  ",124536)
      call add        ("s1  ",2.d0,"tqnw",1.d0,"tqnw")
c  diag 25
      call reord6     ("s2  ","s1  ",123654)
      call reord6     ("ttqf","s2  ",125346)
      call mult6      ("s2  ","s1  ","s3  ",4)
      call reord6     ("s3  ","s1  ",124536)
      call add        ("s1  ",-1.d0,"tqnw",1.d0,"tqnw")
c  diag 26
      call reord6     ("ttpf","s1  ",124635)
      call reorder    ("tq  ","s2  ",2431,4)
      call mult6      ("s1  ","s2  ","s3  ",3)
      call reord6     ("s3  ","s1  ",125346)
      call cmult      ("s1  ",-1.d0,"s1  ")
      call fold6      ("tqnw","s1  ",idum,2)
c  diag 27
      call t1t2       ("t1g ","tq  ","t12q","s1  ",123456,1)
      call t1t2       ("t11 ","tg  ","s1  ","s2  ",123456,2)
      call reord7     ("s2  ","s1  ",312456)
      call reord6     ("s1  ","s2  ",124356)
      call reord7     ("s2  ","s1  ",132456)
      call reord6     ("s1  ","s2  ",126534)
      call reord6     ("ttpf","s1  ",126345)
      call mult6      ("s1  ","s2  ","s3  ",4)
      call reord6     ("s3  ","s1  ",124536)
      call add        ("s1  ",-1.d0,"tqnw",1.d0,"tqnw")
c
      call dveps6     ("tqnw")
      call cmult      ("tqnw",-1.d0,"tqnw")
      call prtime     ("dfoldq  ",0,0)
      return
      end 
C  *******************************************************************
      subroutine boxZ1
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  diag 9
      call reorder    ("tq  ","s1  ",1324,4)
      call add        ("va1 ",2.d0,"va2 ",-1.d0,"s2  ")
      call mult       ("s1  ","s2  ","bxZ1",222)
c  diag 10b
      call reorder    ("tpm ","s1  ",1423,4)
      call add        ("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",1.d0,"bxZ1",2.d0,"bxZ1")
c  diag 20b
      call reorder    ("gpph","s1  ",1243,4)
      call reorder    ("gpph","s2  ",1342,4)
      call add        ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult       ("s1  ","t1v ","s2  ",311)
      call reorder    ("s2  ","s1  ",1423,4)
      call add        ("s1  ",1.d0,"bxZ1",1.d0,"bxZ1")
c  diag 23a
      call reorder    ("vphh","s1  ",1243,4)
      call add        ("s1  ",-2.d0,"vphh",1.d0,"s1  ")
      call mult       ("t1g ","s1  ","s2  ",131)
      call add        ("s2  ",1.d0,"bxZ1",1.d0,"bxZ1")
      return
      end
C  *******************************************************************
      subroutine boxZ2
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  diag 10a
      call reorder    ("tq  ","s1  ",1324,4)
      call add        ("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult       ("s1  ","s2  ","bxZ2",222)
c  diag 11
      call reorder    ("tpm ","s1  ",1423,4)
      call mult       ("s1  ","va1 ","s2  ",222)
      call add        ("s2  ",1.d0,"bxZ2",1.d0,"bxZ2")
c  diag 24b
      call reorder    ("vphh","s1  ",1243,4)
      call mult       ("t1g ","s1  ","s2  ",131)
      call add        ("s2  ",1.d0,"bxZ2",1.d0,"bxZ2")
      return
      end 
C  *******************************************************************
      subroutine boxZ3
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
c  diags 13a,13b
      call reorder    ("pphh","s1  ",1243,4)
      call add        ("pphh",-2.d0,"s1  ",1.d0,"s1  ")
      call mult       ("tgm ","s1  ","bxZ3",113)
c  diags 18a,18b
      call reorder    ("gpph","s1  ",1324,4)
      call add        ("s1  ",2.d0,"gpph",-1.d0,"s1  ")
      call mult       ("s1  ","t1c ","s2  ",202)
      call add        ("s2  ",1.d0,"bxZ3",1.d0,"bxZ3")
      return
      end
C  *******************************************************************
      subroutine boxZ4
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
c  diag 14b
      call reorder    ("pphh","s1  ",3412,4)
      call reorder    ("pphh","s2  ",3421,4)
      call add        ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call mult       ("tvm ","s1  ","bxZ4",113)
c  diag 19a
      call reorder    ("vphh","s1  ",1324,4)
      call reorder    ("vphh","s2  ",1423,4)
      call add        ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call mult       ("s1  ","t1c ","s2  ",202)
      call add        ("s2  ",1.d0,"bxZ4",1.d0,"bxZ4")
      return
      end
C  *******************************************************************
      subroutine boxZ5
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
c  diag 16a
      call move       ("gvhh","bxZ5")
c  diag 25b
      call reorder    ("tp  ","s1  ",1423,4)
      call reorder    ("phhh","s2  ",2314,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",1.d0,"bxZ5",-1.d0,"bxZ5")
c  diag 26b
      call reorder    ("gpph","s1  ",1423,4)
      call mult       ("s1  ","tv  ","s2  ",222)
      call reorder    ("s2  ","s1  ",1342,4)
      call add        ("s1  ",-1.d0,"bxZ5",1.d0,"bxZ5")
c  diag 28b
      call reorder    ("tp  ","s1  ",1432,4)
      call reorder    ("ph  ","s2  ",1243,4)
      call mult       ("s1  ","s2  ","s3  ",311)
      call add        ("s3  ",-1.d0,"bxZ5",1.d0,"bxZ5")
      return
      end
C  *******************************************************************
      subroutine boxZ6
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
c  diag 16b
      call reorder    ("ghvh","bxZ6",1324,4)
c  diag 23b
      call reorder    ("tq  ","s1  ",1324,4)
      call reorder    ("phhh","s2  ",2413,4)
      call reorder    ("phhh","s3  ",2314,4)
      call add        ("s2  ",-2.d0,"s3  ",1.d0,"s2  ")
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",1.d0,"bxZ6",-1.d0,"bxZ6")
c  diag 24a
      call reorder    ("tpm ","s1  ",1423,4)
      call reorder    ("phhh","s2  ",2413,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",1.d0,"bxZ6",1.d0,"bxZ6")
c  diag 25a
      call reorder    ("tg  ","s1  ",1423,4)
      call reorder    ("vphh","s2  ",1423,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","s1  ",1324,4)
      call add        ("s1  ",1.d0,"bxZ6",1.d0,"bxZ6")
c  diag 26a
      call reorder    ("gpph","s1  ",1432,4)
      call mult       ("s1  ","tv  ","s2  ",222)
      call reorder    ("s2  ","s1  ",1342,4)
      call add        ("s1  ",-1.d0,"bxZ6",1.d0,"bxZ6")
c  diag 28a
      call reorder    ("tq  ","s1  ",1342,4)
      call reorder    ("ph  ","s2  ",1243,4)
      call mult       ("s1  ","s2  ","s3  ",311)
      call add        ("s3  ",-1.d0,"bxZ6",1.d0,"bxZ6")
      return
      end
C  *******************************************************************
      subroutine boxZ7
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
c  diag 17a
      call reorder    ("ggph","bxZ7",1243,4)
c  diag 22b
      call reorder    ("tg  ","s1  ",1423,4)
      call mult       ("s1  ","gpph","s2  ",222)
      call reorder    ("s2  ","s1  ",1324,4)
      call add        ("s1  ",-1.d0,"bxZ7",1.d0,"bxZ7")
      return
      end
C  *******************************************************************
      subroutine boxZ8
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
c  diag 17b
      call reorder    ("ggvp","bxZ8",1324,4)
c  diag 20a
      call reorder    ("tq  ","s1  ",1324,4)
      call reorder    ("gpph","s2  ",1324,4)
      call add        ("s2  ",2.d0,"gpph",-1.d0,"s2  ")
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",1.d0,"bxZ8",1.d0,"bxZ8")
c  diag 21a
      call reorder    ("tpm ","s1  ",1423,4)
      call reorder    ("gpph","s2  ",1324,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",-1.d0,"bxZ8",1.d0,"bxZ8")
      return
      end
C  *******************************************************************
      subroutine boxZ9
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
c  diag 21b
      call reorder    ("gpph","s1  ",1324,4)
      call reorder    ("tgm ","s2  ",1423,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","bxZ9",1342,4)   
c  diag 27b
      call reorder    ("phhh","s1  ",2134,4)
      call mult       ("thu ","s1  ","s2  ",222)
      call add        ("s2  ",1.d0,"bxZ9",-1.d0,"bxZ9")
c  diag 29a
      call reorder    ("thu ","s1  ",1243,4)
      call mult       ("s1  ","ph  ","s2  ",311)
      call add        ("s2  ",-1.d0,"bxZ9",1.d0,"bxZ9")
      return
      end
C  *******************************************************************
      subroutine boxZ0
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
c  diag 22a
      call reorder    ("tp  ","s1  ",1423,4)
      call mult       ("gpph","s1  ","s2  ",222)
      call reorder    ("s2  ","bxZ0",1342,4)
c  diag 27a
      call mult       ("thu ","vphh","s1  ",222)
      call add        ("s1  ",1.d0,"bxZ0",-1.d0,"bxZ0")
      return
      end
C  *******************************************************************
      subroutine rdbox_cht
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
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      if (iopt(17).ne.0) then
        call rdfile (41,1,"boxB",-1)
        call rdfile (41,2,"boxC",-1)
        call rdfile (41,3,"boxD",-1)
        call rdfile (41,5,"boxF",-1)
        call rdfile (41,6,"boxG",-1)
        call rdfile (41,7,"boxH",-1)
        call rdfile (41,8,"boxI",-1)
      end if
      if (iopt(17).eq.3) then
        call rdfile (20,1,"bx2 ",-1)
        call rdfile (20,2,"bx3 ",-1)
        call rdfile (20,3,"bx4 ",-1)
        call rdfile (20,4,"bx5 ",-1)
        call rdfile (20,5,"bx6 ",-1)
        call rdfile (31,8,"bx6a",-1)
        call rdfile (20,6,"bx7 ",-1)
        call rdfile (20,7,"bx8 ",-1)
      end if
      return
      end
C  *********************************************************************
      subroutine box25a
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
c  diag 1
      call move      ("ggpp","b25a")
c  diags 2,3
      call reorder   ("ghpp","s1  ",1432,4)
      call mult      ("t1g ","s1  ","s2  ",131)
      call add       ("s2  ",-1.d0,"b25a",1.d0,"b25a")
      call reorder   ("s2  ","s1  ",2143,4)
      call add       ("s1  ",-1.d0,"b25a",1.d0,"b25a")
c   diag 4
      call mult      ("thum","pphh","s1  ",222)
      call add       ("s1  ",1.d0,"b25a",1.d0,"b25a")
      return
      end
C  ******************************************************************
      subroutine box38 (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
      if (iparm.eq.0) then
c  diag 1
        call reorder    ("ggpv","bx38",1243,4)
c  diag 2
        call reorder    ("ghpp","s1  ",1324,4)
        call reorder    ("ghpp","s2  ",1423,4)
        call add        ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
        call reorder    ("tq  ","s2  ",1342,4)
        call mult       ("s1  ","s2  ","s3  ",222)
        call reorder    ("s3  ","s1  ",1342,4)
        call add        ("s1  ",1.d0,"bx38",1.d0,"bx38")
c  diag 3
        call reorder    ("ghpp","s1  ",1324,4)
        call reorder    ("tpm ","s2  ",1432,4)
        call mult       ("s1  ","s2  ","s3  ",222)
        call reorder    ("s3  ","s1  ",1342,4)
        call add        ("s1  ",-1.d0,"bx38",1.d0,"bx38")
c  diag 4
        call reorder    ("tpm ","s1  ",1432,4)
        call reorder    ("ghpp","s2  ",1423,4)
        call mult       ("s1  ","s2  ","s3  ",222)
        call reorder    ("s3  ","s1  ",1324,4)
        call add        ("s1  ",-1.d0,"bx38",1.d0,"bx38")
c  diag 5
        call reorder    ("vphh","s1  ",1243,4)
        call mult       ("thum","s1  ","s2  ",222)
        call add        ("s2  ",1.d0,"bx38",1.d0,"bx38")
c  diag 8
        call reorder    ("gpvh","s1  ",1324,4)
        call mult       ("t1g ","s1  ","s2  ",131)
        call add        ("s2  ",-1.d0,"bx38",1.d0,"bx38")
c  diag 9
        call reorder    ("ghpv","s1  ",1432,4)
        call mult       ("s1  ","t1g ","s2  ",311)
        call reorder    ("s2  ","s1  ",1423,4)
        call add        ("s1  ",-1.d0,"bx38",1.d0,"bx38")
      else if (iparm.eq.1) then
c  diag 6
        call clear      ("bx38",2233,1110)
        call cmult      ("t1c ",0.5d0,"s1  ")
        call reorder    ("t12 ","s2  ",2143,4)
        call t1t2       ("s1  ","s2  ","t12q","s3  ",123456,3)
        call t1t2       ("t1g ","tq  ","s3  ","s1  ",123456,1)
        call reord6     ("s1  ","s2  ",125346)
        call add        ("pphh",-2.d0,"va6 ",1.d0,"s1  ")
        call mult6      ("s2  ","s1  ","s3  ",2)
        call add        ("s3  ",1.d0,"bx38",1.d0,"bx38")
c  diag 7
        call t1t2       ("t1g ","tp  ","t12p","s1  ",123456,1)
        call t1t2       ("t1g ","tp  ","s1  ","s2  ",123456,2)
        call reorder    ("t1v ","s1  ",1243,4)
        call t1t2       ("s1  ","thum","s2  ","s3  ",123456,3)
        call reord6     ("s3  ","s1  ",126345)
        call mult6      ("s1  ","pphh","s2  ",2)
        call add        ("s2  ",1.d0,"bx38",1.d0,"bx38")
      end if
      return
      end
C  ******************************************************************
      subroutine tvsub1 (idi,idj,i,j)
!---------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC

      use dat_mod
      use cfiles_mod
      use cij1_mod
      use tv2_mod
#include "implicit.h"
#include "priunit.h"
#include "strt.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     parameter (ndij=960000)
#include "param_ccsd.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)
#include "cfiles.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!    1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!    2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
#include "strt.h"
!     common /tv2/ a(500,500)
#include "tv2.h"
#endif

!!!!!!!! needs dynamic alloc !!!!!!!!!!!!!!!
      dimension iha(120),ipa(120)

      nh=0
      do 90 ii=1,norb
        if (ic(ii).eq.2) then
          nh=nh+1
          iha(nh)=ii
        end if
90    continue
      np=0
      do 100 ii=1,120
        if (ic(ii).eq.3) then
          np=np+1
          ipa(np)=ii
        end if
100   continue
      inda=iand(ishft(idi,-24),  z'ff')
      indb=iand(ishft(idi,-16),  z'ff')
      indc=iand(ishft(idi,-8),   z'ff')
      indd=iand(idi,             z'ff')
      if (iand(idi,z'ffffff00').eq.iand(idj,z'ffffff00')) then
        ind1=iand(idj,z'ff')
        call files ("hhhp",ntape,new,1,-1,-1)
        nsizea=ije1(nij1)
        istrta1=nstrt(1,ntape-1)+2*nij1+4
        istrta2=nstrt(2,ntape-1)
c  diags 31b,36b
        call files   ("boxD",ntape,new,1,-1,-1)
        nsizeb=ije1(nij1)
        istrtb1=nstrt(1,ntape-1)+2*nij1+4
        istrtb2=nstrt(2,ntape-1)
        sum=0.d0
        do 110 k1=1,nh
          ind2=iha(k1)
          do 120 k2=1,nh
            ind3=iha(k2)
            do 130 k3=1,np
              ind4=ipa(k3)
              nsize2=nsizea
              istrt21=istrta1
              istrt22=istrta2
              term=find4(ind1,ind2,ind3,ind4)
              nsize2=nsizeb
              istrt21=istrtb1
              istrt22=istrtb2
              term=term*find4(indd,ind2,ind3,ind4)
              term=term/(eps(inda)+eps(indb)+eps(ind4)-
     1        eps(indc)-eps(ind3)-eps(ind2))
              if ((icv(ind2).eq.1).and.(icv(ind3).eq.1).
     1        and.(icv(ind4).eq.1)) then
                term=0.d0
              end if
              sum=sum+term
130         continue
120       continue
110     continue
        a(i,j)=a(i,j)+sum
c  diags 34,35a
        nsize2=nsizea
        istrt21=istrta1
        istrt22=istrta2
        call files   ("boxI",ntape,new,1,-1,-1)
        nsize1=ije1(nij1)
        istrt11=nstrt(1,ntape-1)+2*nij1+4
        istrt12=nstrt(2,ntape-1)
        sum=0.d0
        do 140 k1=1,nh
          ind2=iha(k1)
          do 150 k2=1,np
            ind3=ipa(k2)
            term1=find2(ind2,ind3)
            term2=term1
            term1=term1*find4(indd,ind2,ind1,ind3)
            term1=term1/(eps(inda)+eps(indb)+eps(ind3)-
     1      eps(indc)-eps(indd)-eps(ind2))
            term2=term2*find4(ind2,indd,ind1,ind3)
            term2=term2/(eps(inda)+eps(indb)+eps(ind3)-
     1      eps(indc)-eps(ind2)-eps(indd))
            if ((icv(indd).eq.1).and.(icv(ind2).eq.1).
     1      and.(icv(ind3).eq.1)) then
              term1=0.d0
              term2=0.d0
            end if
            sum=sum+2.d0*term1-1.d0*term2
150       continue
140     continue 
        a(i,j)=a(i,j)+sum
      end if
      return
      end
C  *********************************************************************
      subroutine tvsub2
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
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

c  diags 31b,36b
      call reord6    ("t12q","s1  ",125643)
      call mult6     ("s1  ","boxD","s2  ",2)
      call reorder   ("s2  ","s1  ",2134,4)
      call add       ("s1  ",1.d0,"trhs",1.d0,"trhs")
c  diags 33b,38b
      call reord6    ("t12p","s1  ",126543)
      call mult6     ("s1  ","boxH","s2  ",2)
      call reorder   ("s2  ","s1  ",2134,4)
      call add       ("s1  ",1.d0,"trhs",1.d0,"trhs")
c  diags 34,35a 
      call reord6    ("t12q","s1  ",124563)
      call reord6    ("t12q","s2  ",126543)
      call add       ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult6     ("s1  ","boxI","s2  ",0)
      call reorder   ("s2  ","s1  ",2143,4)
      call add       ("s1  ",1.d0,"trhs",1.d0,"trhs")
c  diag 35b 
      call reord6    ("t12p","s1  ",124653)
      call mult6     ("s1  ","boxI","s2  ",0)
      call reorder   ("s2  ","s1  ",2143,4)
      call add       ("s1  ",-1.d0,"trhs",1.d0,"trhs")
      if (iopt(17).eq.3) then
c   diag f12
        call reord6     ("ttpf","s1  ",124653)
        call mult6      ("s1  ","t11 ","s2  ",0)
        call reorder    ("s2  ","s1  ",2143,4)
        call unfold_ccsd("s1  ","s2  ",4)
        call add        ("s2  ",1.d0,"trhs",1.d0,"trhs")
c   diag f13
        call reord6     ("ttqf","s1  ",124563)
        call mult6      ("s1  ","t11 ","s2  ",0)
        call reorder    ("s2  ","s1  ",2143,4)
        call unfold_ccsd("s1  ","s2  ",4)
        call add        ("s2  ",-2.d0,"trhs",1.d0,"trhs")
c   diag f14
        call reord6     ("ttqf","s1  ",125643)
        call mult6      ("s1  ","t11 ","s2  ",0)
        call reorder    ("s2  ","s1  ",2134,4)
        call unfold_ccsd("s1  ","s2  ",4)
        call add        ("s2  ",1.d0,"trhs",1.d0,"trhs")
c   diag f15
        call reord6     ("ttpf","s1  ",126345)
        call reorder    ("tm  ","s2  ",4321,4)
        call mult6      ("s1  ","s2  ","s3  ",2)
        call reorder    ("s3  ","s1  ",2134,4)
        call add        ("s1  ",-1.d0,"trhs",1.d0,"trhs")
c   diag f16
        call reord6     ("ttqf","s1  ",125346)
        call reorder    ("tm  ","s2  ",4321,4)
        call mult6      ("s1  ","s2  ","s3  ",2)
        call reorder    ("s3  ","s1  ",2134,4)
        call add        ("s1  ",2.d0,"trhs",1.d0,"trhs")
c   diag f17
        call reord6     ("ttqf","s1  ",125346)
        call reorder    ("tm  ","s2  ",4312,4)
        call mult6      ("s1  ","s2  ","s3  ",2)
        call reorder    ("s3  ","s1  ",2134,4)
        call add        ("s1  ",-1.d0,"trhs",1.d0,"trhs")
      end if
      return
      end
C  *******************************************************************
