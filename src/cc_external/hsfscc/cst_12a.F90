!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine cst 
!---------------------------------------------------------------------
!
! (1,2)a
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use cint_mod
      use dat_mod
      use cij1_mod
#include "implicit.h"
#include "priunit.h"
#include "csdtime.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     parameter (ndint=3000000,ndij=960000)
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
!     common /csdtime/ csdstim,csdltim
#include "csdtime.h"
#endif

      character*256 filenm,filenm6,                                     &
 filenm51,filenm52,filenm53,filenm54,filenm56,filenm57,filenm58,        &
 filenm20,filenm22,filenm24
      real*4 cputime
      external idec
!     iopt(1) - print.
!     iopt(2) - no. of iterations.
!     iopt(3) - convergence = 10**(-iopt(3))
!     iopt(4) .ne.0 - do rle every iopt(4) iterations.
!     iopt(5) - print only elements .gt. 10**(-iopt(5))
!     iopt(6) - start avrg.
!     iopt(8) - fock.
!     iopt(9) - incomplete model space
!     iopt(10) .ne.0 - do not write tape9
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
      filenm20(1:1)=' '
      call getenv('tape20',filenm20)
      if (filenm20(1:1).eq.' ') filenm20='tape20'
      open(20,file=filenm20,form='UNFORMATTED')
      filenm22(1:1)=' '
      call getenv('tape22',filenm22)
      if (filenm22(1:1).eq.' ') filenm22='tape22'
      open(22,file=filenm22,form='UNFORMATTED')
      filenm24(1:1)=' '
      call getenv('tape24',filenm24)
      if (filenm24(1:1).eq.' ') filenm24='tape24'
      open(24,file=filenm24,form='UNFORMATTED')
      write (*,*) '   $$$$$$$$$   cst   $$$$$$$$$  '
      csdstim = cputime(real(0.))
      csdltim = csdstim
      call rdin
      call rdsym
      convrg = 10.**(-iopt(3))
      write (6,'('' max # of iter='',i3,''. convrg='',e10.2)')          &
  iopt(2),convrg
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
      call rdfile (53,0,"tg  ",-1)
      call rdfile (53,1,"t1g ",-1)    
      call rdfile (53,2,"tgf ",-1)
      call rdfile (54,0,"tu  ",-1)
      call rdfile (54,1,"tuf ",-1)
      call t1sqr   ("t1c ","t1c ","tc  ","tcm ")
      call t1sqr   ("t1v ","t1c ","tv  ","tvm ")
      call reorder ("t1c ","s1  ",1243,4)
      call t1sqr   ("t1g ","s1  ","tg  ","tgm ")
      call t1sqr   ("t1v ","t1v ","tu  ","tum ")
      call reorder ("t1g ","s1  ",1243,4)
      call t1sqr   ("t1v ","s1  ","tp  ","tpm ")
      call ordv_cst
      call files ("vg  ",num,new,1,-1,-1)
      if (ije1(1).gt.0) then
        call calvv10_cst 
        call caltv1_cst
        call fmax   ("t11 ",tmax,indmax)
        write (6,98) tmax,idec(indmax)
98      format (1x,'max t11 = ',f12.6,1x,i12)
        write (6,*) ' '
        call prt ("t11 ","t11 ",0)      
        write (6,*) ' '
      else
        call move ("vg  ","t11 ")
      end if
      call ordta0_cst
      call calvv0_cst 
      call ordtb0_cst
      call caltv2_cst
      call fmax ("t21 ",tmax,indmax)
      write (6,99) tmax,idec(indmax)
99    format (1x,'max t21 = ',f12.6,1x,i12)
      write (6,*) ' '
      call prt ("t21 ","t21 ",0)      
      write (6,*) ' '
      call wrtfile (58,0,"t11 ",-1)
      call wrtfile (58,1,"t21 ",-1)
      write (6,*) nstrt(1,nflcur),' words used in igen,',               &
 nstrt(2,nflcur),' words used in gen'
      call prtime ("end ",0,0)
!     stop
      end

! ********************************************************************
      subroutine ordv_cst
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
#include "dat.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)
#include "fock.h"
#endif

      call readv    (nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call prtime   ("readv   ",nint,0)
      call selt1    (23,0,"hp  ",4)
      call selt1    (33,0,"pp  ",4)
      call selt1    (23,10,"gp  ",4)
      call selt1    (32,11,"vg  ",4)
      call selt1    (32,10,"vh  ",4)
      call seltyp   (2233,0,"hhpp",4)
      call seltyp   (2333,0,"hppp",4)
      call seltyp   (2233,1000,"ghpp",4)
      call seltyp   (2333,1000,"gppp",4)
      call seltyp   (3233,1100,"vgpp",4)
      call seltyp   (3223,1010,"vhgp",4)
      call seltyp   (3223,1000,"vhhp",4)
      call seltyp   (3233,1000,"vhpp",4)
      call seltyp   (3323,1010,"vpgp",4)
      call seltyp   (3322,1001,"vphg",4)
      call seltyp   (3323,1000,"vphp",4)
      call seltyp   (3322,1110,"vvgh",4)
      call seltyp   (3323,1110,"vvgp",4)
      call seltyp   (3322,1100,"vvhh",4)
      call seltyp   (3323,1100,"vvhp",4)
      call reorder  ("hhpp","va1 ",1324,4)
      call reorder  ("hhpp","va2 ",1423,4)
      call prtime   ("ordv",0,0)
      write (6,*)   ' '
      return
      end
! *********************************************************************
      subroutine calvv10_cst
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

!  diag 1
      call move      ("vg  ","vv10")
!  diag 2
      call reorder   ("ghpp","s1  ",1243,4)
      call add       ("ghpp",2.d0,"s1  ",-1.d0,"s1  ")
      call mult      ("tvm ","s1  ","s2  ",113)
      call add       ("s2  ",1.d0,"vv10",1.d0,"vv10")
!  diag 3
      call reorder   ("tq  ","s1  ",1324,4)
      call mult      ("s1  ","hp  ","s2  ",202)
      call add       ("s2  ",2.d0,"vv10",1.d0,"vv10")
!  diag 4
      call reorder   ("tpm ","s1  ",1423,4)
      call mult      ("s1  ","hp  ","s2  ",202)
      call add       ("s2  ",-1.d0,"vv10",1.d0,"vv10")
!  diag 5
      call mult      ("t1v ","gp  ","s1  ",111)
      call add       ("s1  ",1.d0,"vv10",1.d0,"vv10")
!  diag 6
      call reorder   ("vhhp","s1  ",1432,4)
      call reorder   ("vhhp","s2  ",1423,4)
      call add       ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call mult      ("s1  ","tgm ","s2  ",113)
      call add       ("s2  ",1.d0,"vv10",1.d0,"vv10")
!  diag 7
      call reorder   ("vhgp","s1  ",1324,4)
      call reorder   ("vphg","s2  ",1432,4)
      call add       ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult      ("s1  ","t1c ","s2  ",202)
      call add       ("s2  ",1.d0,"vv10",1.d0,"vv10")
!  diag 8
      call mult      ("vh  ","t1g ","s1  ",111)
      call add       ("s1  ",-1.d0,"vv10",1.d0,"vv10")
!  diag 9
      call reorder   ("hhpp","s1  ",3412,4)
      call reorder   ("hhpp","s2  ",3421,4)
      call add       ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call mult      ("t1g ","s1  ","s2  ",131)
      call reorder   ("tvm ","s1  ",1432,4)
      call mult      ("s1  ","s2  ","s3  ",113)
      call add       ("s3  ",1.d0,"vv10",1.d0,"vv10")
!  diag 10
      call reorder   ("hhpp","s1  ",1243,4)
      call add       ("hhpp",-2.d0,"s1  ",1.d0,"s1  ")
      call mult      ("t1v ","s1  ","s2  ",131)
      call reorder   ("tg  ","s1  ",1432,4)
      call mult      ("s2  ","s1  ","s3  ",113)
      call add       ("s3  ",1.d0,"vv10",1.d0,"vv10")
!  diags 11,12
      call reorder   ("hhpp","s1  ",1324,4)
      call reorder   ("hhpp","s2  ",1423,4)
      call add       ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult      ("s1  ","t1c ","s3  ",202)
      call reorder   ("tq  ","s1  ",1324,4)
      call reorder   ("tp  ","s2  ",1423,4)
      call add       ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult      ("s1  ","s3  ","s2  ",202)
      call add       ("s2  ",1.d0,"vv10",1.d0,"vv10")
      call prtime    ("calvv10 ",0,0)
      return
      end
!  **********************************************************************
      subroutine caltv1_cst
!---------------------------------------------------------------------
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

!!!!!!!!!!!!!!!miro: dynamic alloc needed !!!!!!!!!!!!!!!!!!!!!!!!!1
      dimension a(500,500),b(500),indlst(500),ipvt(500),z(500)

      external find2,find4

      do 10 i=1,500
        b(i)=0.d0
        do 20 j=1,500
          a(i,j)=0.d0
20      continue
10    continue
      call files ("vg  ",ntape,new,1,-1,-1)
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
!   diag f1
          call files ("tvf ",ntape,new,1,-1,-1)
          nsize1=ije1(nij1)
          istrt11=nstrt(1,ntape-1)+2*nij1+4
          istrt12=nstrt(2,ntape-1)
         ! if (iand(idi,'00ff'x).eq.iand(idj,'00ff'x)) then
          if (iand(idi,z'00ff').eq.iand(idj,z'00ff')) then
           ! ii1=iand(ishft(idi,-8),'ff'x)
            ii1=iand(ishft(idi,-8),z'ff')
           ! ii2=iand(ishft(idj,-8),'ff'x)
            ii2=iand(ishft(idj,-8),z'ff')
            a(i,j)=a(i,j)-find2(ii1,ii2)
          end if
!   diag f2
          call files ("tgf ",ntape,new,1,-1,-1)
          nsize1=ije1(nij1)
          istrt11=nstrt(1,ntape-1)+2*nij1+4
          istrt12=nstrt(2,ntape-1)
         ! if (iand(idi,'ff00'x).eq.iand(idj,'ff00'x)) then
          if (iand(idi,z'ff00').eq.iand(idj,z'ff00')) then
           ! ii1=iand(idi,'ff'x)
            ii1=iand(idi,z'ff')
            ii2=iand(idj,z'ff')
            a(i,j)=a(i,j)+find2(ii1,ii2)
          end if
!   diag f3
          call files ("tpf ",ntape,new,1,-1,-1)
          nsize2=ije1(nij1)
          istrt21=nstrt(1,ntape-1)+2*nij1+4
          istrt22=nstrt(2,ntape-1)
          !ii1=iand(ishft(idi,-8),'ff'x)
          ii1=iand(ishft(idi,-8),z'ff')
          !ii2=iand(idj,'ff'x)
          ii2=iand(idj,z'ff')
          !ii3=iand(ishft(idj,-8),'ff'x)
          ii3=iand(ishft(idj,-8),z'ff')
          !ii4=iand(idi,'ff'x)
          ii4=iand(idi,z'ff')
          a(i,j)=a(i,j)+find4(ii1,ii2,ii3,ii4)
!   diag f4
          call files ("tqf ",ntape,new,1,-1,-1)
          nsize2=ije1(nij1)
          istrt21=nstrt(1,ntape-1)+2*nij1+4
          istrt22=nstrt(2,ntape-1)
          !ii1=iand(ishft(idi,-8),'ff'x)
          ii1=iand(ishft(idi,-8),z'ff')
          !ii2=iand(idj,'ff'x)
          ii2=iand(idj,z'ff')
          !ii3=iand(idi,'ff'x)
          ii3=iand(idi,z'ff')
          !ii4=iand(ishft(idj,-8),'ff'x)
          ii4=iand(ishft(idj,-8),z'ff')
          a(i,j)=a(i,j)-2.d0*find4(ii1,ii2,ii3,ii4)
50      continue
40    continue
      call clear ("s1  ",32,11)      
      call files ("s1  ",ntape,new,1,-1,-1)
      nsize=ije1(nij1)
      istrt=nstrt(2,ntape-1)
      call dcopy (nsize,-1.d0,0,gen(istrt+1),1)
      call diveps ("s1  ")
      do 80 i=1,nequn
        a(i,i)=a(i,i)+1.d0/gen(istrt+i)
80    continue
      call dgeco (a,500,nequn,ipvt,rcond,z)
      if (abs(rcond).lt.1.d-16) then
        write (6,*) 'Linear equations for t11 are singular.'
        stop
      end if
      call dgesl (a,500,nequn,ipvt,b,0)
      call clear ("t11 ",32,11)       
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
!  *********************************************************************
      subroutine ordta0_cst
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call clear      ("tb1 ",3233,1110)
      call clear      ("tb2 ",3332,1011)
      call clear      ("tb3 ",3332,1101)
      call clear      ("tb4 ",3233,1101)
      return
      end 
!  *********************************************************************
      subroutine calvv0_cst 
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

!  diag 1
      call move      ("vvgp","vv0 ")
!  diag 2b
      call reorder   ("tu  ","s1  ",1243,4)
      call mult      ("s1  ","gp  ","s2  ",311)
      call add       ("s2  ",1.d0,"tb3 ",1.d0,"tb3 ")
!  diag 3a
      call reorder   ("tq  ","s1  ",1342,4)
      call mult      ("s1  ","vh  ","s2  ",311)
      call add       ("s2  ",-1.d0,"tb4 ",1.d0,"tb4 ")
!  diag 3b
      call reorder   ("tp  ","s1  ",1432,4)
      call mult      ("vh  ","s1  ","s2  ",131)
      call add       ("s2  ",-1.d0,"vv0 ",1.d0,"vv0 ")
!  diag 4a
      call reorder   ("tq  ","s1  ",1324,4)
      call reorder   ("vhpp","s2  ",1324,4)
      call add       ("s2  ",2.d0,"vphp",-1.d0,"s2  ")
      call mult      ("s1  ","s2  ","s3  ",222)
      call add       ("s3  ",1.d0,"tb1 ",1.d0,"tb1 ")
!  diag 4b
      call reorder   ("tv  ","s1  ",1324,4)
      call reorder   ("vhgp","s2  ",1324,4)
      call reorder   ("vphg","s3  ",1432,4)
      call add       ("s2  ",2.d0,"s3  ",-1.d0,"s2  ")
      call mult      ("s2  ","s1  ","s3  ",222)
      call add       ("s3  ",1.d0,"tb1 ",1.d0,"tb1 ")
!  diag 5a
      call reorder   ("tpm ","s1  ",1423,4)
      call reorder   ("vhpp","s2  ",1324,4)
      call mult      ("s1  ","s2  ","s3  ",222)
      call add       ("s3  ",-1.d0,"tb1 ",1.d0,"tb1 ")
!  diag 5b
      call reorder   ("tvm ","s1  ",1423,4)
      call reorder   ("vhgp","s2  ",1324,4)
      call mult      ("s2  ","s1  ","s3  ",222)
      call add       ("s3  ",-1.d0,"tb1 ",1.d0,"tb1 ")
!  diag 6a
      call reorder   ("tvm ","s1  ",1423,4)
      call reorder   ("vphg","s2  ",1432,4)
      call mult      ("s1  ","s2  ","s3  ",222)
      call add       ("s3  ",-1.d0,"tb2 ",1.d0,"tb2 ")
!  diag 6b
      call reorder   ("tpm ","s1  ",1423,4)
      call mult      ("vphp","s1  ","s2  ",222)
      call add       ("s2  ",-1.d0,"tb2 ",1.d0,"tb2 ")
!  diag 7 
      call mult      ("tum ","gppp","s1  ",222)
      call add       ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
!  diag 8 
      call mult      ("vvhh","tgm ","s1  ",222)
      call add       ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
!  diag 12
      call reorder   ("tpm ","s1  ",1423,4)
      call mult      ("s1  ","va2 ","s2  ",222)
      call reorder   ("tvm ","s1  ",1423,4)
      call mult      ("s1  ","s2  ","s3  ",222)
      call add       ("s3  ",1.d0,"tb2 ",1.d0,"tb2 ")
!  diag 15
      call mult      ("tum ","hhpp","s1  ",222)
      call mult      ("s1  ","tgm ","s2  ",222)
      call add       ("s2  ",1.d0,"vv0 ",1.d0,"vv0 ")
!  
      call boxZ1_cst
      call reorder   ("tv  ","s1  ",1324,4)
      call mult      ("bxZ1","s1  ","s2  ",222)
      call add       ("s2  ",1.d0,"tb1 ",1.d0,"tb1 ")
!  
      call boxZ2_cst
      call reorder   ("tvm ","s1  ",1423,4)
      call mult      ("bxZ2","s1  ","s2  ",222)
      call add       ("s2  ",1.d0,"tb1 ",1.d0,"tb1 ")
!  
      call boxZ3_cst
      call reorder   ("tq  ","s1  ",1342,4)
      call mult      ("s1  ","bxZ3","s2  ",311)
      call add       ("s2  ",1.d0,"tb4 ",1.d0,"tb4 ")
      call reorder   ("tp  ","s1  ",1432,4)
      call mult      ("bxZ3","s1  ","s2  ",131)
      call add       ("s2  ",1.d0,"vv0 ",1.d0,"vv0 ")
!  
      call boxZ4_cst
      call reorder   ("tu  ","s1  ",1243,4)
      call mult      ("s1  ","bxZ4","s2  ",311)
      call add       ("s2  ",1.d0,"tb3 ",1.d0,"tb3 ")
!  
      call boxZ5_cst
      call mult      ("t1v ","bxZ5","s1  ",131)
      call add       ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
!  
      call boxZ6_cst
      call mult      ("bxZ6","t1v ","s2  ",311)
      call add       ("s2  ",1.d0,"tb4 ",1.d0,"tb4 ")
!  
      call boxZ7_cst
      call mult      ("bxZ7","t1g ","s1  ",311)
      call add       ("s1  ",1.d0,"tb3 ",1.d0,"tb3 ")
!
      call boxZ8_cst
      call reorder   ("t1c ","s1  ",1243,4)
      call mult      ("bxZ8","s1  ","s2  ",311)
      call add       ("s2  ",1.d0,"tb1 ",1.d0,"tb1 ")
!  
      call boxZ9_cst
      call mult      ("bxZ9","t1g ","s1  ",311)
      call add       ("s1  ",1.d0,"tb3 ",1.d0,"tb3 ")
!  
      call boxZ0_cst
      call reorder   ("t1c ","s1  ",1243,4)
      call mult      ("bxZ0","s1  ","s2  ",311)
      call add       ("s2  ",1.d0,"vv0 ",1.d0,"vv0 ")
!  diag f4a
      call reorder   ("t1v ","s1  ",1243,4)
      call mult      ("tuf ","s1  ","s2  ",311)
      call reorder   ("s2  ","s1  ",1243,4)
      call reorder   ("t11 ","s2  ",1243,4)
      call mult      ("s1  ","s2  ","s3  ",311)
      call add       ("s3  ",-1.d0,"tb3 ",1.d0,"tb3 ")
!  diag f5 
      call reorder   ("tuf ","s1  ",1243,4)
      call reorder   ("t11 ","s2  ",1243,4)
      call mult      ("s1  ","s2  ","s3  ",311)
      call unfold_ccsd("s3  ","s1  ",3)
      call add       ("s1  ",-1.d0,"tb3 ",1.d0,"tb3 ")
!  diag f6a
      call reorder   ("tpf ","s1  ",1432,4)
      call mult      ("t11 ","s1  ","s2  ",131)
      call reorder   ("t1v ","s1  ",1243,4)
      call mult      ("s2  ","s1  ","s3  ",311)
      call add       ("s3  ",1.d0,"vv0 ",1.d0,"vv0 ")
!  diag f8 
      call reorder   ("tpf ","s1  ",1432,4)
      call mult      ("t11 ","s1  ","s2  ",131)
      call unfold_ccsd("s2  ","s1  ",4)
      call add       ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
!  diag f10a
      call reorder   ("tqf ","s1  ",1324,4)
      call reorder   ("t1v ","s2  ",1243,4)
      call mult      ("s1  ","s2  ","s3  ",311)
      call reorder   ("s3  ","s1  ",1243,4)
      call mult      ("s1  ","t11 ","s2  ",311)
      call add       ("s2  ",1.d0,"tb4 ",1.d0,"tb4 ")
!  diag f11 
      call reorder   ("tqf ","s1  ",1342,4)
      call mult      ("s1  ","t11 ","s2  ",311)
      call unfold_ccsd("s2  ","s1  ",3)
      call add       ("s1  ",1.d0,"tb4 ",1.d0,"tb4 ")
!  dgt1p4
      call mult      ("t1v ","hhpp","s1  ",131)
      call mult      ("t1v ","s1  ","s2  ",131)
      call reorder   ("t1c ","s1  ",1243,4)
      call mult      ("s1  ","s2  ","s3  ",131)
      call mult      ("t1g ","s3  ","s1  ",131)
      call reorder   ("s1  ","s2  ",3412,4)
      call add       ("s2  ",-2.d0,"vv0 ",1.d0,"vv0 ")
      call prtime    ("calvv0  ",0,0)
      return
      end
!  *********************************************************************
      subroutine ordtb0_cst
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
      return
      end
! *********************************************************************
      subroutine caltv2_cst
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

!!!!!!!!!!!! needed dynamic alloc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension a(500,500),b(500),indlst(500),ipvt(500),z(500)

      external find2,find4

      do 10 i=1,500
        b(i)=0.d0
        do 20 j=1,500
          a(i,j)=0.d0
20      continue
10    continue
      call files ("vvgp",ntape,new,1,-1,-1)
      istrt1=nstrt(1,ntape-1)+2*nij1+4
      nequn=ije1(nij1)
      do 30 i=1,nequn
        indlst(i)=igen(istrt1+i)
30    continue
      call files ("vv0 ",ntape,new,1,-1,-1)
      do 60 i=1,nequn
        b(i)=-1.d0*gen(nstrt(2,ntape-1)+i)
60    continue
!  diag 2a
      call move       ("pp  ","boxX")
!  diag 14a
      call reorder    ("hhpp","s1  ",3412,4)
      call reorder    ("hhpp","s2  ",3421,4)
      call add        ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call reorder    ("tcm ","s2  ",3412,4)
      call mult       ("s1  ","s2  ","s3  ",113)
      call add        ("s3  ",1.d0,"boxX",1.d0,"boxX")
!  diag 19b
      call reorder    ("hppp","s1  ",2413,4)
      call reorder    ("hppp","s2  ",3412,4)
      call add        ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult       ("s1  ","t1c ","s2  ",202)
      call add        ("s2  ",1.d0,"boxX",1.d0,"boxX")
!  diag 29b
      call reorder    ("hp  ","s1  ",1243,4)
      call reorder    ("t1c ","s2  ",1243,4)
      call mult       ("s1  ","s2  ","s3  ",111)
      call add        ("s3  ",-1.d0,"boxX",1.d0,"boxX")
      do 40 i=1,nequn 
        idi=indlst(i)
        do 50 j=1,nequn
          idj=indlst(j)
!   diags boxX
          call files ("boxX",ntape,new,1,-1,-1)
          nsize1=ije1(nij1)
          istrt11=nstrt(1,ntape-1)+2*nij1+4
          istrt12=nstrt(2,ntape-1)
          if (iand(idi,z'ffffff00').eq.iand(idj,z'ffffff00')) then
            ii1=iand(idj,z'ff')
            ii2=iand(idi,z'ff')
            a(i,j)=a(i,j)+find2(ii1,ii2)
          end if
!   diag f1
          call files ("tvf ",ntape,new,1,-1,-1)
          nsize1=ije1(nij1)
          istrt11=nstrt(1,ntape-1)+2*nij1+4
          istrt12=nstrt(2,ntape-1)
          if (iand(idi,z'ff00ffff').eq.iand(idj,z'ff00ffff')) then
            ii1=iand(ishft(idi,-16),z'ff')
            ii2=iand(ishft(idj,-16),z'ff')
            a(i,j)=a(i,j)-find2(ii1,ii2)
          end if
!   diag f2
          if (iand(idi,z'00ffffff').eq.iand(idj,z'00ffffff')) then
            ii1=iand(ishft(idi,-24),z'ff')
            ii2=iand(ishft(idj,-24),z'ff')
            a(i,j)=a(i,j)-find2(ii1,ii2)
          end if
!   diag f3
          call files ("tgf ",ntape,new,1,-1,-1)
          nsize1=ije1(nij1)
          istrt11=nstrt(1,ntape-1)+2*nij1+4
          istrt12=nstrt(2,ntape-1)
          if (iand(idi,z'ffff00ff').eq.iand(idj,z'ffff00ff')) then
            ii1=iand(ishft(idi,-8),z'ff')
            ii2=iand(ishft(idj,-8),z'ff')
            a(i,j)=a(i,j)+find2(ii1,ii2)
          end if
!   diag f4b
          call files ("tuf ",ntape,new,1,-1,-1)
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
!   diag f6b
          call files ("tpf ",ntape,new,1,-1,-1)
          nsize2=ije1(nij1)
          istrt21=nstrt(1,ntape-1)+2*nij1+4
          istrt22=nstrt(2,ntape-1)
          if (iand(idi,z'ff0000ff').eq.iand(idj,z'ff0000ff')) then
            ii1=iand(ishft(idi,-16), z'ff')
            ii2=iand(ishft(idj,-8),  z'ff')
            ii3=iand(ishft(idj,-16), z'ff')
            ii4=iand(ishft(idi,-8),  z'ff')
            a(i,j)=a(i,j)+find4(ii1,ii2,ii3,ii4)
          end if
!   diag f7
          if (iand(idi,z'00ff00ff').eq.iand(idj,z'00ff00ff')) then
            ii1=iand(ishft(idi,-24),  z'ff')
            ii2=iand(ishft(idj,-8),   z'ff')
            ii3=iand(ishft(idj,-24),  z'ff')
            ii4=iand(ishft(idi,-8),   z'ff')
            a(i,j)=a(i,j)+find4(ii1,ii2,ii3,ii4)
          end if
!   diag f9
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
!   diag f10b
          if ((iand(idi,z'ff').eq.iand(idj,z'ff')).and.                 &
     (iand(ishft(idi,-16),z'ff').eq.                                    &
     iand(ishft(idj,-24),z'ff'))) then
            ii1=iand(ishft(idi,-24),  z'ff')
            ii2=iand(ishft(idj,-8),   z'ff')
            ii3=iand(ishft(idi,-8),   z'ff')
            ii4=iand(ishft(idj,-16),  z'ff')
            a(i,j)=a(i,j)+find4(ii1,ii2,ii3,ii4)
          end if
50      continue
40    continue
      call clear ("s1  ",3323,1110)
      call files ("s1  ",ntape,new,1,-1,-1)
      nsize=ije1(nij1)
      istrt=nstrt(2,ntape-1)
      call dcopy (nsize,-1.d0,0,gen(istrt+1),1)
      call diveps ("s1  ")
      do 80 i=1,nequn
        a(i,i)=a(i,i)+1.d0/gen(istrt+i)
80    continue
      call dgeco (a,500,nequn,ipvt,rcond,z)
      if (abs(rcond).lt.1.d-16) then
        write (6,*) 'Linear equations for t21 are singular.'
        stop
      end if
      call dgesl (a,500,nequn,ipvt,b,0)
      call clear ("t21 ",3323,1110)
      call files ("t21 ",ntape,new,1,-1,-1)
      nsize=ije1(nij1)
      istrt=nstrt(2,ntape-1)
      do 70 i=1,nsize
        gen(istrt+i)=b(i)
70    continue
      call prtime ("caltv2  ",0,0)
      write (6,*) ' '
      return
      end
!  *********************************************************************
      subroutine boxZ1_cst
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  diag 9
      call reorder    ("tq  ","s1  ",1324,4)
      call add        ("va1 ",2.d0,"va2 ",-1.d0,"s2  ")
      call mult       ("s1  ","s2  ","bxZ1",222)
!  diag 10b
      call reorder    ("tpm ","s1  ",1423,4)
      call add        ("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",1.d0,"bxZ1",2.d0,"bxZ1")
!  diag 20b
      call reorder    ("vhhp","s1  ",1243,4)
      call reorder    ("vhhp","s2  ",1342,4)
      call add        ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call mult       ("s1  ","t1g ","s2  ",311)
      call reorder    ("s2  ","s1  ",1423,4)
      call add        ("s1  ",1.d0,"bxZ1",1.d0,"bxZ1")
!  diag 23a
      call reorder    ("ghpp","s1  ",1243,4)
      call add        ("s1  ",2.d0,"ghpp",-1.d0,"s1  ")
      call mult       ("t1v ","s1  ","s2  ",131)
      call add        ("s2  ",1.d0,"bxZ1",1.d0,"bxZ1")
      return
      end
!  *******************************************************************
      subroutine boxZ2_cst
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

!  diag 10a
      call reorder    ("tq  ","s1  ",1324,4)
      call add        ("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult       ("s1  ","s2  ","bxZ2",222)
!  diag 11
      call reorder    ("tpm ","s1  ",1423,4)
      call mult       ("s1  ","va1 ","s2  ",222)
      call add        ("s2  ",1.d0,"bxZ2",1.d0,"bxZ2")
!  diag 24b
      call reorder    ("ghpp","s1  ",1243,4)
      call mult       ("t1v ","s1  ","s2  ",131)
      call add        ("s2  ",-1.d0,"bxZ2",1.d0,"bxZ2")
      return
      end 
!  *******************************************************************
      subroutine boxZ3_cst
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  diags 13a,13b
      call reorder    ("hhpp","s1  ",1243,4)
      call add        ("hhpp",-2.d0,"s1  ",1.d0,"s1  ")
      call mult       ("tvm ","s1  ","bxZ3",113)
!  diags 18a,18b
      call reorder    ("vhhp","s1  ",1324,4)
      call add        ("s1  ",-2.d0,"vhhp",1.d0,"s1  ")
      call mult       ("s1  ","t1c ","s2  ",202)
      call add        ("s2  ",1.d0,"bxZ3",1.d0,"bxZ3")
      return
      end
!  *******************************************************************
      subroutine boxZ4_cst
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  diag 14b
      call reorder    ("hhpp","s1  ",3412,4)
      call reorder    ("hhpp","s2  ",3421,4)
      call add        ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
      call mult       ("tgm ","s1  ","bxZ4",113)
!  diag 19a
      call reorder    ("ghpp","s1  ",1324,4)
      call reorder    ("ghpp","s2  ",1423,4)
      call add        ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult       ("s1  ","t1c ","s2  ",202)
      call add        ("s2  ",1.d0,"bxZ4",1.d0,"bxZ4")
      return
      end
!  *******************************************************************
      subroutine boxZ5_cst
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  diag 16a
      call move       ("vgpp","bxZ5")
!  diag 25b
      call reorder    ("tp  ","s1  ",1423,4)
      call reorder    ("hppp","s2  ",2314,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",-1.d0,"bxZ5",1.d0,"bxZ5")
!  diag 26b
      call reorder    ("vhhp","s1  ",1423,4)
      call mult       ("s1  ","tg  ","s2  ",222)
      call reorder    ("s2  ","s1  ",1342,4)
      call add        ("s1  ",1.d0,"bxZ5",1.d0,"bxZ5")
!  diag 28b
      call reorder    ("tp  ","s1  ",1432,4)
      call reorder    ("hp  ","s2  ",1243,4)
      call mult       ("s1  ","s2  ","s3  ",311)
      call add        ("s3  ",-1.d0,"bxZ5",1.d0,"bxZ5")
      return
      end
!  *******************************************************************
      subroutine boxZ6_cst
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

!  diag 16b
      call reorder    ("vpgp","bxZ6",1324,4)
!  diag 23b
      call reorder    ("tq  ","s1  ",1324,4)
      call reorder    ("hppp","s2  ",2413,4)
      call reorder    ("hppp","s3  ",2314,4)
      call add        ("s2  ",2.d0,"s3  ",-1.d0,"s2  ")
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",1.d0,"bxZ6",1.d0,"bxZ6")
!  diag 24a
      call reorder    ("tpm ","s1  ",1423,4)
      call reorder    ("hppp","s2  ",2413,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",-1.d0,"bxZ6",1.d0,"bxZ6")
!  diag 25a
      call reorder    ("tv  ","s1  ",1423,4)
      call reorder    ("ghpp","s2  ",1423,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","s1  ",1324,4)
      call add        ("s1  ",-1.d0,"bxZ6",1.d0,"bxZ6")
!  diag 26a
      call reorder    ("vhhp","s1  ",1432,4)
      call mult       ("s1  ","tg  ","s2  ",222)
      call reorder    ("s2  ","s1  ",1342,4)
      call add        ("s1  ",1.d0,"bxZ6",1.d0,"bxZ6")
!  diag 28a
      call reorder    ("tq  ","s1  ",1342,4)
      call reorder    ("hp  ","s2  ",1243,4)
      call mult       ("s1  ","s2  ","s3  ",311)
      call add        ("s3  ",-1.d0,"bxZ6",1.d0,"bxZ6")
      return
      end
!  *******************************************************************
      subroutine boxZ7_cst
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

!  diag 17a
      call reorder    ("vvhp","bxZ7",1243,4)
!  diag 22b
      call reorder    ("tv  ","s1  ",1423,4)
      call mult       ("s1  ","vhhp","s2  ",222)
      call reorder    ("s2  ","s1  ",1324,4)
      call add        ("s1  ",1.d0,"bxZ7",-1.d0,"bxZ7")
      return
      end
!  *******************************************************************
      subroutine boxZ8_cst
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

!  diag 17b
      call reorder    ("vvgh","bxZ8",1324,4)
!  diag 20a
      call reorder    ("tq  ","s1  ",1324,4)
      call reorder    ("vhhp","s2  ",1324,4)
      call add        ("s2  ",-2.d0,"vhhp",1.d0,"s2  ")
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",1.d0,"bxZ8",-1.d0,"bxZ8")
!  diag 21a
      call reorder    ("tpm ","s1  ",1423,4)
      call reorder    ("vhhp","s2  ",1324,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call add        ("s3  ",1.d0,"bxZ8",1.d0,"bxZ8")
      return
      end
!  *******************************************************************
      subroutine boxZ9_cst
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  diag 21b
      call reorder    ("vhhp","s1  ",1324,4)
      call reorder    ("tvm ","s2  ",1423,4)
      call mult       ("s1  ","s2  ","s3  ",222)
      call reorder    ("s3  ","bxZ9",1342,4)   
!  diag 27b
      call reorder    ("hppp","s1  ",2134,4)
      call mult       ("tu  ","s1  ","s2  ",222)
      call add        ("s2  ",-1.d0,"bxZ9",1.d0,"bxZ9")
!  diag 29a
      call reorder    ("tu  ","s1  ",1243,4)
      call mult       ("s1  ","hp  ","s2  ",311)
      call add        ("s2  ",-1.d0,"bxZ9",1.d0,"bxZ9")
      return
      end
!  *******************************************************************
      subroutine boxZ0_cst
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  diag 22a
      call reorder    ("tp  ","s1  ",1423,4)
      call mult       ("vhhp","s1  ","s2  ",222)
      call reorder    ("s2  ","bxZ0",1342,4)
!  diag 27a
      call mult       ("tu  ","ghpp","s1  ",222)
      call add        ("s1  ",-1.d0,"bxZ0",1.d0,"bxZ0")
      return
      end
!  *******************************************************************
