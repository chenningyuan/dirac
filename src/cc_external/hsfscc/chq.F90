!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine chq 
!---------------------------------------------------------------------
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
!    1 iiia(200),iiva(200)a
#include "cfiles.h"
!     common /cint/ nint,idum1,vint(ndint),iint(ndint)
#include "cint.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /csdtime/ csdstim,csdltim
#include "csdtime.h"
#endif

      character*256 filenm,filenm1,filenm6,                             &
 filenm51,filenm52,filenm53,filenm56,filenm57,                          &
 filenm20,filenm21,filenm31,filenm41,filenm42
      real*4 cputime

!     iopt(1) - print.
!     iopt(2) - no. of iterations.
!     iopt(3) - convergence = 10**(-iopt(3))
!     iopt(4) .ne.0 - do rle every iopt(4) iterations.
!     iopt(5) - print only elements .gt. 10**(-iopt(5))
!     iopt(6) - start avrg.
!     iopt(8) - fock.
!     iopt(9) - incomplete model space
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
      filenm56(1:1)=' '
      call getenv('tape56',filenm56)
      if (filenm56(1:1).eq.' ') filenm56='tape56'
      open(56,file=filenm56,form='UNFORMATTED')
      filenm57(1:1)=' '
      call getenv('tape57',filenm57)
      if (filenm57(1:1).eq.' ') filenm57='tape57'
      open(57,file=filenm57,form='UNFORMATTED')
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
      filenm41(1:1)=' '
      call getenv('tape41',filenm41)
      if (filenm41(1:1).eq.' ') filenm41='tape41'
      open(41,file=filenm41,form='UNFORMATTED')
      filenm42(1:1)=' '
      call getenv('tape42',filenm42)
      if (filenm42(1:1).eq.' ') filenm42='tape42'
      open(42,file=filenm42,form='UNFORMATTED')
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')
      write (*,*) '  $$$$$$$$$$$   chq    $$$$$$$$$$  '
      csdstim = cputime(real(0.))
      csdltim = csdstim
      call rdin
      call rdsym
      if (iopt(10).eq.0) then
        filenm(1:1)=' '
        call getenv('tape9',filenm)
        if (filenm(1:1).eq.' ') filenm='tape9'
        open(9,file=filenm,form='UNFORMATTED')
      endif
      convrg = 10.**(-iopt(3))
      write (6,'('' max # of iter='',i3,''. convrg='',e10.2)')          &
  iopt(2),convrg
      call prtime("rdin",0,0)
      if (iopt(9).ne.0) call rdrsp2
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      call rdfile(51,0,"tc  ",-1)
      call rdfile(51,1,"t1c ",-1)
      if (iopt(17).ne.0) then
        call rdfile(51,2,"t3c ",-1)
      end if
      call rdfile(52,0,"tv  ",-1)
      call rdfile(52,1,"t1v ",-1)
      call rdfile(52,2,"tvf ",-1)
      call rdfile(56,0,"thp ",-1)
      call rdfile(56,1,"thpf",-1)
      if (iopt(17).ne.0) then
        call rdfile(52,3,"t3v ",-1)
        call rdfile(56,2,"t3p ",-1)
      end if
      call rdfile(53,0,"tg  ",-1)
      call rdfile(53,1,"t1g ",-1)
      call rdfile(53,2,"tgf ",-1)
      if (iopt(17).ne.0) then
        call rdfile(53,3,"t3g ",-1)
      end if
      if (iopt(17).ne.0) then
!       call rdbox
        call rdbox_chq
      end if
      call t1sqr("t1c ","t1c ","tc  ","tcm ")
      call t1sqr("t1v ","t1c ","tv  ","tvm ")
      call reorder("t1c ","s2  ",1243,5)
      call t1sqr("t1g ","s2  ","tg  ","tgm ")
      call reorder("t1g ","s1  ",1243,4)
      call t1sqr("t1v ","s1  ","thp ","tpm ")
!     call ordvv
      call ordvv_chq
!     call calvv0
      call calvv0_chq
!     call pdiag
      call pdiag_chq
      iter=0
!     call tvinit
      call tvinit_chq
      if (iopt(4).ne.0) then
        if (iopt(17).eq.0) then
          call dbdiag ("dum1","vv0 ","dum3","dum4",0100,1)
        else
          call dbdiag ("dum1","vv0 ","v30 ","dum4",0110,1)
        end if
        call wrtrle ("thq ",2,1)    
        if (iopt(17).ne.0) then
          call wrtrle ("t3q ",3,1)   
        end if
      endif
      call prt("thq ","thq ",2)
      call prt("thqf","thqf",2)
      if (iopt(17).ne.0) then
        call prt("t3q ","t3q ",2)
      end if
      call fmax("thq ",tmax,itmax)
      write(6,2400) tmax,idec(itmax)
220   iter=iter+1
      call prtime("start iter",iter,0)
!     call ordtva
      call ordtva_chq
      if (iopt(17).ne.0) then
!       call caltv3
        call caltv3_chq
      end if
!     call caltv2
      call caltv2_chq
!     call ordtvb
      call ordtvb_chq
      if (iopt(6).ne.0.and.iter.ge.mod(iopt(6),100)) then
        fnew=max(dfloat(iopt(6)/100),1.d0)
        if (fnew.le.6.5d0) then
          fnew=1.d0/fnew
        else
          fnew=fnew-5.d0
        endif
        call add("tnew",fnew/(1.d0+fnew),"thq ",1.d0/(1.d0+fnew),"tnew")
        write (6,'('' t averaged with fnew '',f7.3)') fnew
      endif
      call fmax("tnew",tmax,itmax)
      t3max = 0.0d0 ! mi init
      if (iopt(17).ne.0) then
        call fmax6 ("t3nw",t3max,id3mx,ind3mx)
      end if
      call diffmax("thq ","tnew","s2  ",diff,idiff)
      if (iopt(17).ne.0) then
        call dfmax3 ("t3q ","t3nw","s2  ",df3,iddf3,inddf3)
      end if
      write (6,*) ' '
      write (6,'('' iteration # '',i3,'' diffmx='',                     &
  e10.3,2x,i12)') iter,diff
      df3=0.0d0 ! mi init
      if (iopt(17).ne.0) then
        write (6,2401) df3,idec(iddf3),idec(inddf3)
2401    format (24x,e10.3,2x,i12,2x,i12)
      end if
      write(6,2400) tmax,idec(itmax)
2400  format(' max t=',f15.9,2x,i12)
      if (iopt(17).ne.0) then
        write (6,2402) t3max,idec(id3mx),idec(ind3mx)
2402    format (' max t3=',f15.9,2x,i12,2x,i12)
      end if
!     if (abs(tmax).gt.1.d10.or.abs(tmax).lt.1.d-10) then
      if (abs(tmax).gt.1.d10) then
        write (*,*) ' abort. tmax out of bounds'
        call abext('tmax    ')    
      endif
      call move("tnew","thq ")
      if (iopt(17).ne.0) then
        call move("t3nw","t3q ")
      end if
      call prt("tnew","new t   ",2)
      if (iopt(17).ne.0) then
        call prt("t3nw","new t3  ",2)
      end if
      call prt("thqf","thqf",2)
      if (max(abs(diff),abs(df3)).lt.convrg) go to 910
      if (iopt(4).ne.0) go to 260
!   no rle
      if (iter.ge.iopt(2)) go to 920
      go to 220
!   with rle
260   continue
      itr1=mod(iter-1,iopt(4))+1
      call wrtrle ("thq ",2,itr1+1)      
      if (iopt(17).ne.0) then
        call wrtrle ("t3q ",3,itr1+1)      
      end if
      if (itr1.lt.iopt(4)) go to 220
      if (iopt(17).eq.0) then
        call rle (itr1,0100,"dum1","thq ","t3q ","dum4",ising)
      else
        call rle (itr1,0110,"dum1","thq ","t3q ","dum4",ising)
      end if
!  rle singular
      if (ising.eq.1) then
        if (iopt(17).eq.0) then
          call rle (itr1-1,0100,"dum1","thq ","t3q ","dum4",ising)
        else
          call rle (itr1-1,0110,"dum1","thq ","t3q ","dum4",ising)
        end if
        if (ising.eq.1) then
          iopt(4)=0
          write(*,*) ' ***** rle singular, turned off'
          go to 220
        endif
      endif
      call wrtrle ("thq ",2,1)      
      if (iopt(17).ne.0) then
        call wrtrle ("t3q ",3,1)      
      end if
      if (iter.ge.iopt(2)) go to 920
      go to 220
!
!  different endings
910   write (6,*) ' chq  converged after iter ',iter
      go to 970
920   write (6,*) ' iteration limit'
      call abext('limit   ')      
970   call prt("thq ","last t  ",0)
      call prt("thqf","last tf ",0)
      call ehp(iopt(1)+1)
      call wrtfile(57,0,"thq ")
      call wrtfile(57,1,"thqf")
      if (iopt(17).ne.0) then
        call wrtfile(57,2,"t3q ")
        if (iopt(17).eq.3) then
!         call wrtbox
          call wrtbox_chq
          call prt6 ("t3p ","last t3p",0)
          call prt6 ("t3q ","last t3q",0)
        end if
      end if
      write (*,*) nstrt(1,nflcur),' words used in igen,',               &
   nstrt(2,nflcur),' words used in gen'
      call prtime("end     ",iter,0)
!      stop
      end

! *********************************************
      subroutine ordvv_chq
!---------------------------------------------------------------------
! va: hhpp
! vb: hphp
! vc: hhhp
! vd: hppp
! vva: vhpp
! vvc: vhhp
!              hhpp 1324 = va1
!              hhpp 1423 = va2
!              hhpp 3412 = va4
!              hhpp 1243 = va6
!              hhpp 3421 = va8
!              hhhp 3412 = vc4
!              hhhp 3421 = vc8
!              hhhp 3124 = vc9
!              hhhp 3214 = vc10
!              hppp 2314 = vd7
!              hppp 2413 = vd11
!              hppp 2134 = vd12
!              hppp 2143 = vd14
!              vhpp 1324 = vva1
!              vhhp 1324 = vvc1
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cint_mod
      use fock_mod
      use dat_mod
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
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)
#include "fock.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

!
      call readv(nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call prtime("readv   ",nint,0)
      call seltyp(2233,0,"hhpp",4)
      call seltyp(2323,0,"hphp",4)
      call seltyp(2223,0,"hhhp",4)
      call seltyp(2322,0,"hphh",4)
      call seltyp(2333,0,"hppp",4)
      call seltyp(3323,1000,"vphp",4)
      call seltyp(3223,1000,"vhhp",4)
      call seltyp(3232,1000,"vhph",4)
      call seltyp(2322,100,"hvhh",4)
      call seltyp(3222,1000,"vhhh",4)
      call seltyp(2232,1000,"ghph",4)
      call seltyp(2233,100,"hgpp",4)
      call seltyp(2233,1000,"ghpp",4)
      call seltyp(3223,1100,"vghp",4)
      call seltyp(3223,1010,"vhgp",4)
      call seltyp(3323,1010,"vpgp",4)
      call seltyp(3222,1010,"vhgh",4)
      call reorder("hhpp","va1 ",1324,4)
      call reorder("hhpp","va2 ",1423,4)
      call reorder("hhpp","va3 ",1243,4)
      call reorder("hhpp","va4 ",3412,4)
      call reorder("hhpp","va8 ",3421,4)
      call reorder("hhhp","vc2 ",1243,4)
      call reorder("hppp","vd7 ",1324,4)
      call reorder("hppp","vd8 ",2314,4)
      call reorder("hppp","vd11",2413,4)
      call reorder("hppp","vd12",2134,4)
      call reorder("hppp","vd14",2143,4)
      call reorder("vhhp","vvc2",1243,4)
      call reorder("vhhp","vvc3",1342,4)
      call reorder("ghpp","vga1",1243,4)
      call reorder("vhgp","vgv ",1324,4)
      call selt1(33,0,"pp  ",4)
      call selt1(23,0,"hp  ",4)
      call selt1(22,0,"hh  ",4)
      call selt1(32,10,"vh  ",4)
      call selt1(23,10,"gp  ",4)
      call reorder("hp  ","ph  ",1243,4)
      if (iopt(17).ne.0) then
        call seltyp   (3322,0010,"ppgh",4)
      end if
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
        call seltyp   (3223,1111,"vggv",4)
      end if
      if (iopt(17).eq.3) then
        call seltyp4  (3333,0,"pppp",4)
        call seltyp   (2322,0001,"hphg",4)
        call seltyp   (3322,1001,"vphg",4)
        call seltyp   (2333,1000,"gppp",4)
        call reorder  ("hhpp","va6 ",1243,4)
      end if
      call prtime("ordv",0,0)
      return
      end
! ****************************************
      subroutine pdiag_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call pdig5c
      call pdig14c
      call pdig4c
      call pdig9c
      call pdig13c
      return
      end
! ******************************************
      subroutine tvinit_chq
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

      call move("vv0 ","thq ")
      call fold_ccsd("thq ","thqf")
      call diveps("thq ")
      if (iopt(17).ne.0) then
        call move     ("v30 ","t3q ")
        call reord6   ("t3q ","s1  ",123465)
        call reord7   ("s1  ","s2  ",132456)
        call add      ("t3q ",1.d0,"s2  ",1.d0,"t3q ")
        call dveps6   ("t3q ")
      end if
      return
      end
! **************************************************
      subroutine ordta0_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder("tc  ","tc1 ",1423,4)
      call reorder("tc  ","tc2 ",1324,4)
      call reorder("tcm ","tcm1",1423,4)
      call reorder("tv  ","tv1 ",1243,4)
      call reorder("tvm ","tvm2",1423,4)
      call reorder("tg  ","tg3 ",1243,4)
      call reorder("tgm ","tgm2",1432,4)
      call reorder("t1c ","t1cr",1243,4)
      call reorder("t1v ","t1vr",1243,4)
      call reorder("t1g ","t1gr",1243,4)
      call reorder("tg  ","tgr ",3412,4)
      call clear("tb1 ",3322,1010)
      call clear("tb2 ",3223,1100)
      call clear("tb4 ",3232,1001)
      call clear("tb6 ",3232,1100)
      call clear("tb15",3322,1001)
      call clear("thqf",3223,1111)
      call prtime("ordta0  ",0,1)
      return
      end
! ******************************************************
      subroutine calvv0_chq
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

      call ordta0_chq
      call move("vhgp","vv0 ")
      call diag7c
      call diag16_chq
      call diag17_chq
      call diag20
      call dgt1p4_chq
      call diag11_chq
      call diag12c
      call diag24_chq
      call diag15_chq
      call diag13
      call diag26_chq
      call diag27c
      call diag2c
      call diag3c
      call diag8_chq
      call diag4d
      call diag5c1
      call diag5d
      call diag6c
      call diag6d
      call diag14d
      call diag22c
      call diag27d
      if (iopt(17).ne.0) then
!  diags 30a,37a
        call reord6     ("t3v ","s1  ",125346)
        call mult6      ("s1  ","boxB","s2  ",2)   
        call add        ("s2  ",1.d0,"tb4 ",1.d0,"tb4 ")
!  diags 31a,36a
        call reord6     ("t3g ","s1  ",125643)
        call mult6      ("boxC","s1  ","s2  ",2)   
        call add        ("s2  ",1.d0,"tb6 ",1.d0,"tb6 ")
!  diags 32a,39a
        call reord6     ("t3v ","s1  ",126345)
        call mult6      ("s1  ","boxF","s2  ",2)   
        call add        ("s2  ",1.d0,"tb4 ",1.d0,"tb4 ")
!  diags 32d,39d
        call reord7     ("t3p ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125346)
        call mult6      ("s2  ","boxE","s1  ",2)   
        call add        ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
!  diags 33a,38a
        call reord6     ("t3g ","s1  ",126543)
        call mult6      ("boxG","s1  ","s2  ",2)   
        call add        ("s2  ",1.d0,"tb6 ",1.d0,"tb6 ")
!  diags 33d,38d
        call reord6     ("t3p ","s1  ",125346)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",126345)
        call mult6      ("s1  ","boxH","s2  ",2)   
        call add        ("s2  ",1.d0,"tb6 ",1.d0,"tb6 ")
!  diag 35a
        call reord7     ("t3p ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125634)
        call mult6      ("s2  ","boxI","s1  ",0)
        call add        ("s1  ",-1.d0,"vv0 ",1.d0,"vv0 ")
      end if
!     call ordtb0
      call ordtb0_chq
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
!  dg3n1e,2e
        call reorder    ("ppgh","s1  ",3421,4)
        call mult6      ("tv1 ","s1  ","s2  ",1)
        call reorder    ("tg  ","s1  ",1423,4)
        call mult6      ("vvc2","s1  ","s3  ",1)
        call add        ("s3  ",-1.d0,"s2  ",1.d0,"s1  ")
        call reord6     ("s1  ","v30 ",125436)
!  dg3n1i
        call mult6      ("vpgp","tc  ","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125463)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!  dg3n2i,f3
        call reorder    ("tc  ","s1  ",1432,4)
        call mult6      ("vhgh","s1  ","s2  ",1)
        call cmult      ("s2  ",-1.d0,"s2  ")
        call reorder    ("tv  ","s1  ",2341,4)
        call mult6      ("vggv","s1  ","s3  ",1)
        call fold6      ("s2  ","s3  ",idum,2)
        call reord6     ("s2  ","s1  ",124356)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
!  dg3nf4
        call reorder    ("vggv","s1  ",1342,4)
        call reorder    ("tg  ","s2  ",3421,4)
        call mult6      ("s1  ","s2  ","s3  ",1)
        call reord6     ("s3  ","s1  ",124356)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",125346)
        call cmult      ("s1  ",-1.d0,"s1  ")
        call fold6      ("v30 ","s1  ",idum,2)                
      else if (iopt(17).eq.3) then
!  GROUP I
!   diag 2
        call reads      (31,1,"s1  ")
        call mult6      ("tv1 ","s1  ","s2  ",1)
        call reord6     ("s2  ","v30 ",126435)
!   diag 3
        call box34      (0)
        call mult6      ("bx34","tc  ","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125364)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP II
!   diag 2
        call reads      (21,1,"s1  ")
        call mult6      ("s1  ","tg3 ","s2  ",1)
        call reord6     ("s2  ","s1  ",126435)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
!   diag 3
        call box35      (0)
        call reorder    ("tc  ","s1  ",1432,4)
        call mult6      ("bx35","s1  ","s2  ",1)
        call reord6     ("s2  ","s1  ",124356)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP III
!   diag 2
        call box36      (0)
        call reord6     ("t3c ","s1  ",124536)
        call mult6      ("bx36","s1  ","s2  ",3)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",124356)
        call add        ("s2  ",0.5d0,"v30 ",1.d0,"v30 ")
!  GROUP IV
!   diag 2
        call reord7     ("t3p ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125634)
        call mult6      ("s2  ","bx4 ","s1  ",3)
        call reord6     ("s1  ","s2  ",125346)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!   diag 3
        call box37      (0)
        call reord6     ("t3c ","s1  ",124635)
        call mult6      ("bx37","s1  ","s2  ",3)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",124356)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP V
!   diag 2
        call reads      (31,5,"s1  ")
        call reord6     ("t3v ","s2  ",125643)
        call mult6      ("s2  ","s1  ","s3  ",3)
        call reord6     ("s3  ","s1  ",126534)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
!   diag 3
        call reads      (21,5,"s1  ")
        call reord6     ("t3g ","s2  ",125643)
        call mult6      ("s1  ","s2  ","s3  ",3)
        call reord6     ("s3  ","s1  ",125346)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",126453)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP VI
!   diag 2
        call reads      (21,6,"s1  ")
        call reord6     ("t3g ","s2  ",123645)
        call mult6      ("s1  ","s2  ","s3  ",3)
        call reord6     ("s3  ","s1  ",126345)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP VII
!   diag 4
        call reads      (31,6,"s1  ")
        call reord6     ("t3v ","s2  ",123645)
        call mult6      ("s2  ","s1  ","s3  ",3)
        call reord6     ("s3  ","s1  ",123564)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
!   diag 5
        call reorder    ("gppp","s1  ",1243,4)
        call mult       ("t1v ","s1  ","s2  ",131)
        call mult6      ("s2  ","tc  ","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125364)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP VIII
!   diag 2
        call reads      (31,7,"s1  ")
        call reord6     ("t3v ","s2  ",123564)
        call mult6      ("s2  ","s1  ","s3  ",5)
        call reord6     ("s3  ","s1  ",123645)
        call add        ("s1  ",0.5d0,"v30 ",1.d0,"v30 ")
!  GROUP IX   
!   diag 2
        call reads      (21,7,"s1  ")
        call reord6     ("t3g ","s2  ",123564)
        call mult6      ("s2  ","s1  ","s3  ",5)
        call reord6     ("s3  ","s1  ",126345)
        call reord7     ("s1  ","s2  ",312456)
        call reord6     ("s2  ","s1  ",125346)
        call reord7     ("s1  ","s2  ",132454)
        call reord6     ("s2  ","s1  ",126345)
        call add        ("s1  ",0.5d0,"v30 ",1.d0,"v30 ")
      end if
      return
      end
!  *********************************************************************
      subroutine ordtb0_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder("tb1 ","s2  ",1432,4)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call reorder("tb2 ","s2  ",1324,4)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call reorder("tb4 ","s2  ",1243,4)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call reorder("tb6 ","s2  ",1423,4)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call reorder("tb15","s2  ",1342,4)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call prtime("ordtb0  ",0,1)
      return
      end
! *****************************************************
      subroutine ordtva_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder("thq ","ta1 ",1342,4)
      call reorder("thq ","ta3 ",1324,4)
      call clear("tb2 ",3223,1100)
      call clear("tb6 ",3232,1100)
      call clear("tb5 ",2332,1010)
      call clear("tb4 ",3232,1001)
      call clear("thqf",3223,1111)
      call prtime("ordta   ",0,1)
      return
      end
!  *******************************************************************
      subroutine caltv3_chq
!---------------------------------------------------------------------
!
!
!---------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8(a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      call move         ("v30 ","t3nw")
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
!  dg3n1a
        call mult6      ("thq ","hppp","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call add        ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  dg3n2a
        call reorder    ("hhhp","s1  ",1243,4)
        call mult6      ("ta1 ","s1  ","s2  ",1)
        call reord6     ("s2  ","s1  ",124356)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",125346)
        call add        ("s1  ",-1.d0,"t3nw",1.d0,"t3nw")
      else if (iopt(17).eq.3) then
!  GROUP I
!   diag 1
        call mult6      ("thq ","bx1 ","s1  ",1)
        call reord6     ("s1  ","s2  ",124365)
        call add        ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 3
        call box34      (1)
        call mult6      ("bx34","tc  ","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125364)
        call add        ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP II
!   diag 1
        call mult6      ("ta1 ","bx2 ","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125364)
        call add        ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 3
        call box35      (1)
        call reorder    ("tc  ","s1  ",1432,4)
        call mult6      ("bx35","s1  ","s2  ",1)
        call reord6     ("s2  ","s1  ",124356)
        call add        ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP III
!   diag 1
        call reord6     ("t3q ","s1  ",124536)
        call mult6      ("s1  ","bx3 ","s2  ",3)
        call reord6     ("s2  ","s1  ",125346)
        call add        ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 2
        call box36      (1)
        call reord6     ("t3c ","s1  ",124536)
        call mult6      ("bx36","s1  ","s2  ",3)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",124356)
        call add        ("s2  ",0.5d0,"t3nw",1.d0,"t3nw")
!  GROUP IV
!   diag 1
        call reord6     ("t3q ","s1  ",124635)
        call mult6      ("s1  ","bx4 ","s2  ",3)
        call reord6     ("s2  ","s1  ",125346)
        call add        ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 3
        call box37      (1)
        call reord6     ("t3c ","s1  ",124635)
        call mult6      ("bx37","s1  ","s2  ",3)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",124356)
        call add        ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP V
!   diag 1
        call reord6     ("t3q ","s1  ",124635)
        call mult6      ("s1  ","bx5 ","s2  ",3)
        call reord6     ("s2  ","s1  ",125364)
        call add        ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP VI
!   diag 1
        call reord6     ("t3q ","s1  ",124356)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",125634)
        call mult6      ("s1  ","bx6 ","s2  ",3)
        call reord6     ("s2  ","s1  ",125346)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",126345)
        call add        ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP VII
!   diag 1
        call t1t2       ("t1c ","thq ","t3q ","s1  ",123456,2)
        call t1t2       ("t1c ","thq ","s1  ","s2  ",123456,3)
        call mult6      ("s2  ","pppp","s1  ",3)
        call reord6     ("s1  ","s2  ",123465)
        call reord7     ("s2  ","s3  ",132456)
        call add        ("s1  ",1.d0,"s3  ",1.d0,"s1  ")
        call add        ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
!   diag 2
        call mult6      ("t3q ","vd14","s1  ",3)
        call mult6      ("s1  ","t1cr","s2  ",5)
        call add        ("s2  ",-1.d0,"t3nw",1.d0,"t3nw")
!   diag 3
        call mult6      ("t3q ","hhpp","s1  ",3)
        call reorder    ("tcm ","s2  ",3412,4)
        call mult6      ("s1  ","s2  ","s3  ",3)
        call add        ("s3  ",0.5d0,"t3nw",1.d0,"t3nw")
!  GROUP VIII
!   diag 1
        call mult6      ("t3q ","bx7 ","s1  ",5)
        call add        ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP IX
!   diag 1
        call reord6     ("t3q ","s1  ",124563)
        call mult6      ("s1  ","bx8 ","s2  ",5)
        call reord6     ("s2  ","s1  ",126345)
        call add        ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
      end if
      call prtime       ("caltv3  ",0,0)
      return
      end
!  *******************************************************************
      subroutine caltv2_chq
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

      call move("vv0 ","tnew")
      call diag2d
      call diag3d
      call diag4c
      call diag5c
      call diag9c
      call diag14c
      call diag13c
      if (iopt(17).ne.0) then
!  diags 30d,37d
        call reord6     ("t3q ","s1  ",124365)
        call mult6      ("s1  ","boxA","s2  ",2)
        call add        ("s2  ",1.d0,"tnew",1.d0,"tnew")
!  diags 31d,36d
        call reord6     ("t3q ","s1  ",124356)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",125436)
        call mult6      ("s1  ","boxD","s2  ",2)
        call add        ("s2  ",1.d0,"tb6 ",1.d0,"tb6 ")
!  diags 34a,35d
        call reord6     ("t3q ","s1  ",124536)
        call reord6     ("t3q ","s2  ",124635)
        call add        ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
        call mult6      ("s1  ","boxI","s2  ",0)
        call add        ("s2  ",1.d0,"tnew",1.d0,"tnew")    
      end if
      call prtime       ("caltv2  ",0,1)
      return
      end
! ******************************************************
      subroutine ordtvb_chq
!---------------------------------------------------------------------
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

      i=1
      call reorder("tb5 ","s2  ",3412,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
10    call reorder("tb2 ","s2  ",1324,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call reorder("tb4 ","s2  ",1243,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call reorder("tb6 ","s2  ",1423,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      if (i.eq.77) go to 11
      call fold_ccsd("tnew","thqf")
!     call heff2
      i=77
      call clear("tb2 ",3223,1100)
      call clear("tb6 ",3232,1100)
      call dfoldhq
      if (iopt(17).eq.3) then
!       call dfold3
        call dfold3_chq
      end if
      go to 10
11    call diveps("tnew")
      if (iopt(17).ne.0) then
        call reord6  ("t3nw","s1  ",123465)
        call reord7  ("s1  ","s2  ",132456)
        call add     ("t3nw",1.d0,"s2  ",1.d0,"t3nw")
        call dveps6  ("t3nw")
      end if
      call prtime("ordtvb  ",0,1)
      return
      end
! *************************************************
      subroutine diag7c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call seltyp(2333,1000,"s2  ",5)
      call mult("tvm ","s2  ","s3  ",222)
      call add("vv0 ",1.d0,"s3  ",1.d0,"vv0 ")
      call prt("s3  ","diag7c  ",3)
      call prtime("diag7c  ",0,1)
      return
      end
! ************************************************
      subroutine diag2d
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("thq ","pp  ","s2  ",311)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","diag2d  ",3)
      call prtime("diag2d  ",0,1)
      return
      end
! ************************************************
      subroutine diag3d
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("ta1 ","hh  ","s2  ",311)
      call add("tb6 ",1.d0,"s2  ",-1.d0,"tb6 ")
      call prt("s2  ","-diag3d ",3)
      call prtime("diag3d  ",0,1)
      return
      end
! **********************************************
      subroutine pdig14c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  incl 14c, 19c, 29c
!  14c
      call add("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
      call reorder("tcm ","s1  ",3412,4)
      call mult("s1  ","s2  ","p14c",113)
!  19c
      call add("vd11",2.d0,"vd8 ",-1.d0,"s2  ")
      call mult("t1c ","s2  ","s1  ",022)
      call add("p14c",1.d0,"s1  ",1.d0,"p14c")
!  29c
      call mult("t1cr","ph  ","s2  ",111)
      call add("p14c",1.d0,"s2  ",1.d0,"p14c")
      call prt("p14c","pdig14c ",3)
      call prtime("pdig14c ",0,1)
      return
      end
! ********************************************
      subroutine diag14c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("thq ","p14c","s2  ",311)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","diag14c ",3)
      call prtime("diag14c ",0,1)
      return
      end
! **************************************************
      subroutine pdig9c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 9c 10c
      call add("tc2 ",-2.d0,"tc1 ",1.d0,"s2  ")
      call add("va1 ",-2.d0,"va2 ",1.d0,"s3  ")
      call mult("s2  ","s3  ","p9c ",222)
      call prt("p9c ","pdig9c  ",3)
      call prtime("pdig9c  ",0,1)
      return
      end
! **************************************************
      subroutine diag9c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("ta3 ","p9c ","s2  ",222)
      call add("tb2 ",1.d0,"s2  ",1.d0,"tb2 ")
      call prt("s2  ","diag9c  ",3)
      call prtime("diag9c  ",0,1)
      return
      end
! ************************************************************
      subroutine pdig13c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  incl 13c 18c 28c
!  13c
      call add("hhpp",-2.d0,"va3 ",1.d0,"s2  ")
      call mult("tcm ","s2  ","p13c",113)
!  18c
      call reorder("hhhp","s2  ",1324,5)
      call add("s2  ",-2.d0,"hhhp",1.d0,"s2  ")
      call mult("t1c ","s2  ","s1  ",022)
      call add("p13c",1.d0,"s1  ",1.d0,"p13c")
!  28c
      call mult("t1c ","hp  ","s1  ",111)
      call add("p13c",1.d0,"s1  ",1.d0,"p13c")
      call prt("p13c","pdig13c ",3)
      call prtime("pdig13c ",0,1)
      return
      end
! ************************************************************
      subroutine diag13c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("ta1 ","p13c","s2  ",311)
      call add("tb6 ",1.d0,"s2  ",1.d0,"tb6 ")
      call prt("s2  ","diag13c ",3)
      call prtime("diag13c ",0,1)
      return
      end
! **************************************************
      subroutine diag12c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  incl 12c, 22d, 25c
!  12c
      call mult("tgm2","va2 ","s3  ",222)
! 22d
      call reorder("hhhp","s1  ",1342,5)
      call mult("t1g ","s1  ","s2  ",131)
      call add("s3  ",1.d0,"s2  ",1.d0,"s3  ")
!  25c
      call mult("t1c ","ghpp","s2  ",131)
      call reorder("s2  ","s1  ",2134,5)
      call add("s3  ",1.d0,"s1  ",-1.d0,"s3  ")
! --------------
      call mult("tvm2","s3  ","s2  ",222)
      call add("tb1 ",1.d0,"s2  ",1.d0,"tb1 ")
      call prt("s2  ","diag12c ",3)
      call prtime("diag12c ",0,1)
      return
      end
! *********************************************
      subroutine diag4d
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add("vgv ",2.d0,"vghp",-1.d0,"s2  ")
      call mult("s2  ","tc2 ","s1  ",222)
      call add("tb2 ",1.d0,"s1  ",1.d0,"tb2 ")
      call prt("s1  ","diag4d  ",3)
      call prtime("diag4d  ",0,1)
      return
      end
! *****************************************************
      subroutine pdig4c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  4c
      call add("va1 ",2.d0,"hphp",-1.d0,"p4c ")
!  20c
      call reorder("hphh","s1  ",1243,5)
      call add("s1  ",-2.d0,"hphh",1.d0,"s1  ")
      call mult("t1cr","s1  ","s2  ",131)
      call reorder("s2  ","s1  ",2143,5)
      call add("p4c ",1.d0,"s1  ",1.d0,"p4c ")
!  23c
      call reorder("hppp","s1  ",1243,5)
      call add("s1  ",-1.d0,"hppp",2.d0,"s1  ")
      call mult("t1c ","s1  ","s2  ",131)
      call reorder("s2  ","s1  ",1324,5)
      call add("p4c ",1.d0,"s1  ",1.d0,"p4c ")
      call prt("p4c ","pdig4c  ",3)
      call prtime("pdig4c  ",0,1)
      return
      end
! ********************************************
      subroutine diag4c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("ta3 ","p4c ","s1  ",222)
      call add("tb2 ",1.d0,"s1  ",1.d0,"tb2 ")
      call prt("s1  ","diag4c  ",3)
      call prtime("diag4c  ",0,1)
      return
      end
! ********************************************
      subroutine diag5c1
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("t1v ","va1 ","s1  ",131)
      call mult("s1  ","t1g ","s2  ",311)
      call add("tb4 ",1.d0,"s2  ",-1.d0,"tb4 ")
      call prt("s2  ","-diag5c1",3)
      call prtime("diag5c1 ",0,1)
      return
      end
! *******************************************
      subroutine pdig5c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

!  incl 5c, 10d, 11c, 21c, 24c
!  10d
      call add("va1 ",-2.d0,"va2 ",1.d0,"s1  ")
      call mult("tc2 ","s1  ","p5c ",222)
!  5c
      call add("p5c ",1.d0,"va1 ",-1.d0,"p5c ")
!  11c
      call mult("tc1 ","va1 ","s1  ",222)
      call add("p5c ",1.d0,"s1  ",1.d0,"p5c ")
!  21c
      call mult("t1cr","vc2 ","s1  ",131)
      call reorder("s1  ","s3  ",2134,5)
      call add("p5c ",1.d0,"s3  ",1.d0,"p5c ")
!  24c
      call mult("t1c ","vd12","s1  ",131)
      call add("p5c ",1.d0,"s1  ",-1.d0,"p5c ")
      call prt("p5c ","pdig5c  ",3)
      call prtime("pdig5c  ",0,1)
      return
      end
! ********************************************
      subroutine diag5c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder("thp ","s3  ",1423,4)
      call mult("s3  ","p5c ","s1  ",222)
      call add("tb2 ",1.d0,"s1  ",1.d0,"tb2 ")
      call prt("s1  ","diag5c  ",3)
      call prtime("diag5c  ",0,1)
      return
      end
! *********************************************************
      subroutine diag5d
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call mult("vgv ","tcm1","s1  ",222)
      call add("tb2 ",1.d0,"s1  ",-1.d0,"tb2 ")
      call prt("s1  ","-diag5d ",3)
      call prtime("diag5d  ",0,1)
      return
      end
! *******************************************************
      subroutine diag6c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder("ghph","s1  ",1243,5)
      call mult("tvm2","s1  ","s2  ",222)
      call add("tb1 ",1.d0,"s2  ",-1.d0,"tb1 ")
      call prt("s2  ","-diag6c ",3)
      call prtime("diag6c  ",0,1)
      return
      end
! *********************************************
      subroutine diag6d
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vphp","tgm2","s1  ",222)
      call add("tb1 ",1.d0,"s1  ",-1.d0,"tb1 ")
      call prt("s1  ","-diag6d ",3)
      call prtime("diag6d  ",0,1)
      return
      end
! ********************************************************
      subroutine diag20
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl 20d, 23d
!  20d
      call add("vvc2",-2.d0,"vvc3",1.d0,"s2  ")
      call mult("t1g ","s2  ","s1  ",131)
      call reorder("s1  ","s3  ",2134,5)
!   diag 23d
      call add("vga1",2.d0,"ghpp",-1.d0,"s2  ")
      call mult("t1v ","s2  ","s1  ",131)
      call add("s3  ",1.d0,"s1  ",1.d0,"s3  ")
! ---------------------------------
      call mult("s3  ","tc2 ","s2  ",222)
      call add("tb2 ",1.d0,"s2  ",1.d0,"tb2 ")
      call prt("s2  ","diag20d ",3)
      call prtime("diag20d ",0,1)
      return
      end
! ***************************************************
      subroutine diag26_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  incl 26c 26d
!  26c
      call mult("t1c ","vhhp","s1  ",131)
      call reorder("s1  ","s3  ",2143,5)
!  26d
      call mult("t1v ","hhhp","s1  ",131)
! ---------------------
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
      call mult("s1  ","tgm ","s2  ",222)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call prt("s2  ","diag26  ",3)
      call prtime("diag26  ",0,1)
      return
      end
! *************************************************
      subroutine diag27c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tv  ","ghpp","s1  ",222)
      call mult("s1  ","t1cr","s2  ",311)
      call add("vv0 ",1.d0,"s2  ",-1.d0,"vv0 ")
      call prt("s2  ","-diag27c",3)
      call prtime("diag27c ",0,1)
      return
      end
! *************************************************
      subroutine diag2c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tv1 ","gp  ","s1  ",311)
      call add("tb4 ",1.d0,"s1  ",1.d0,"tb4 ")
      call prt("s1  ","diag2c  ",3)
      call prtime("diag2c  ",0,1)
      return
      end
! *************************************************
      subroutine diag3c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vh  ","tg3 ","s1  ",131)
      call add("tb6 ",1.d0,"s1  ",-1.d0,"tb6 ")
      call prt("s1  ","-diag3c ",3)
      call prtime("diag3c  ",0,1)
      return
      end
! **************************************************
      subroutine diag8_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call seltyp(3222,1000,"s2  ",5)
      call mult("s2  ","tgm ","s1  ",222)
      call add("vv0 ",1.d0,"s1  ",1.d0,"vv0 ")
      call prt("s1  ","diag8   ",3)
      call prtime("diag8   ",0,1)
      return
      end
! **************************************************
      subroutine diag11_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 10d.1, 11c.1, 21c.1, 24c.1
!  10d.1
      call add("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult("tc2 ","s2  ","s3  ",222)
!  11c.1
      call mult("tcm1","va1 ","s1  ",222)
      call add("s3  ",1.d0,"s1  ",1.d0,"s3  ")
!  21c.1
      call mult("t1cr","vc2 ","s1  ",131)
      call reorder("s1  ","s2  ",2134,5)
      call add("s3  ",1.d0,"s2  ",1.d0,"s3  ")
!  24c.1
      call mult("t1c ","vd7 ","s1  ",131)
      call reorder("s1  ","s2  ",1423,5)
      call add("s3  ",1.d0,"s2  ",-1.d0,"s3  ")
! ---------------------
      call mult("t1v ","s3  ","s2  ",131)
      call mult("s2  ","t1g ","s3  ",311)
      call add("tb4 ",1.d0,"s3  ",1.d0,"tb4 ")
      call prt("s3  ","diag11  ",3)
      call prtime("diag11  ",0,1)
      return
      end
! ***************************************************
      subroutine diag24_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 21d, 24d
!  24
      call mult("t1v ","vga1","s1  ",131)
!   21
      call mult("t1g ","vvc2","s2  ",131)
      !call reorder("s2  ","s3  ",2134,5,0)
      call reorder("s2  ","s3  ",2134,5)
      call add("s3  ",1.d0,"s1  ",-1.d0,"s3  ")
! ----------------------
      call mult("s3  ","tc1 ","s2  ",222)
      call add("tb2 ",1.d0,"s2  ",1.d0,"tb2 ")
      call prt("s2  ","diag24  ",3)
      call prtime("diag24  ",0,1)
      return
      end
! ************************************************************
      subroutine diag13
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl 13d 18d 28d
!  13d
      call add("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      !call reorder("tvm ","s3  ",1324,4,0)
      call reorder("tvm ","s3  ",1324,4)
      call mult("s3  ","s2  ","s1  ",113)
!   18
      !call reorder("vhhp","s3  ",1324,4,0)
      call reorder("vhhp","s3  ",1324,4)
      call add("s3  ",-2.d0,"vhhp",1.d0,"s2  ")
      call mult("s2  ","t1c ","s3  ",202)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
!  28
      call mult("t1v ","hp  ","s2  ",111)
      call add("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
! ----------------------
      call mult("s1  ","tg3 ","s2  ",131)
      call add("tb6 ",1.d0,"s2  ",1.d0,"tb6 ")
      call prt("s2  ","diag13  ",3)
      call prtime("diag13  ",0,1)
      return
      end
! ******************************************************
      subroutine diag14d
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl 14d, 19d
      call add("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
      call mult("tgm ","s2  ","s3  ",113)
!   19
      !call reorder("ghpp","s1  ",1324,5,0)
      call reorder("ghpp","s1  ",1324,5)
      !call reorder("ghpp","s2  ",1423,5,0)
      call reorder("ghpp","s2  ",1423,5)
      call add("s1  ",2.d0,"s2  ",-1.d0,"s1  ")
      call mult("t1c ","s1  ","s2  ",022)
      call add("s2  ",1.d0,"s3  ",1.d0,"s2  ")
! ------------------
      call mult("tv1 ","s2  ","s1  ",311)
      call add("tb4 ",1.d0,"s1  ",1.d0,"tb4 ")
      call prt("s1  ","diag14d ",3)
      call prtime("diag14d ",0,1)
      return
      end
! ******************************************************
      subroutine diag22c
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  incl 22c, 25d
      call mult("vvc3","t1cr","s1  ",311)
      !call reorder("s1  ","s3  ",1423,5,0)
      call reorder("s1  ","s3  ",1423,5)
!   diag 25d
      !call reorder("hppp","s1  ",4123,5,0)
      call reorder("hppp","s1  ",4123,5)
      call mult("t1v ","s1  ","s2  ",131)
      call add("s3  ",1.d0,"s2  ",-1.d0,"s3  ")
! ----------------------------------
      !call reorder("tg  ","s2  ",1432,4,0)
      call reorder("tg  ","s2  ",1432,4)
      call mult("s3  ","s2  ","s1  ",222)
      call add("tb1 ",1.d0,"s1  ",1.d0,"tb1 ")
      call prt("s1  ","diag22c ",3)
      call prtime("diag22c ",0,1)
      return
      end
! ************************************************
      subroutine diag15_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
! 15c
      call mult("tvm ","hhpp","s1  ",222)
      call mult("s1  ","tgm ","s2  ",222)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call prt("s2  ","diag15  ",3)
      call prtime("diag15  ",0,1)
      return
      end
! ************************************************************
      subroutine diag17_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  incl 17c 17d
!  17d
      call mult("vhgh","t1cr","s2  ",311)
      call add("vv0 ",1.d0,"s2  ",-1.d0,"vv0 ")
      call prt("s2  ","-diag17d",3)
      call prtime("diag17d ",0,1)
!  17c
      !call reorder("vhhp","s1  ",1243,5,0)
      call reorder("vhhp","s1  ",1243,5)
      call mult("s1  ","t1g ","s2  ",311)
      call add("tb4 ",1.d0,"s2  ",-1.d0,"tb4 ")
      call prt("s2  ","-diag17c",3)
      call prtime("diag17c ",0,1)
      return
      end
! **********************************************************
      subroutine diag27d
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl  27d, 29d
!   27d
      call mult("tv  ","hppp","s1  ",222)
      !call reorder("s1  ","s2  ",1243,5,0)
      call reorder("s1  ","s2  ",1243,5)
!   29a
      call mult("tv1 ","hp  ","s1  ",311)
      call add("s1  ",-1.d0,"s2  ",-1.d0,"s1  ")
! ------------
      call mult("s1  ","t1g ","s2  ",311)
      call add("tb4 ",1.d0,"s2  ",1.d0,"tb4 ")
      call prt("s2  ","diag27d ",3)
      call prtime("diag27d ",0,1)
      return
      end
! ********************************************
      subroutine diag16_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  incl 16c 16d
!  16c
      call mult("t1v ","hgpp","s2  ",131)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call prt("s2  ","diag16c ",3)
      call prtime("diag16c ",0,1)
!  16d
      call mult("vpgp","t1c ","s2  ",311)
      call add("tb1 ",1.d0,"s2  ",1.d0,"tb1 ")
      call prt("s2  ","diag16d ",3)
      call prtime("diag16d ",0,1)
      return
      end
! ************************************************
      subroutine dgt1p4_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  -2*t1**4 compensates for overcounting.
!  11c.1
      call mult("t1c ","va1 ","s1  ",131)
      call mult("t1cr","s1  ","s2  ",131)
      !call reorder("s2  ","s1  ",2134,5,0)
      call reorder("s2  ","s1  ",2134,5)
!  21c.1
      call mult("t1cr","vc2 ","s2  ",131)
      !call reorder("s2  ","s3  ",2134,5,0)
      call reorder("s2  ","s3  ",2134,5)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
!  -2*t1**4 ?
      call mult("t1v ","s1  ","s3  ",131)
      call mult("s3  ","t1g ","s1  ",311)
      call add("tb4 ",1.d0,"s1  ",-2.d0,"tb4 ")
      call prt("s3  ","-.5t1**4",3)
      call prtime("dgt1p4  ",0,1)
      return
      end
! *************************************************
      subroutine dfoldhq
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

!  f2c
!     call reorder("thq ","s1  ",1243,5,1)
      call reorder("thq ","s1  ",1243,5)
      call mult("s1  ","tgf ","s2  ",311)
      call add("tb4 ",1.d0,"s2  ",1.d0,"tb4 ")
      call prt("s2  ","f2c ",3)
!  f3c
!     call reorder("thqf","s3  ",1342,5,1)
      call reorder("thqf","s3  ",1342,5)
      call mult("t1gr","s3  ","s2  ",131)
!     call reorder("s2  ","s1  ",2134,5,1)
      call reorder("s2  ","s1  ",2134,5)
      call unfold_ccsd("s1  ","s2  ",4)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s1  ","f3a ",3)
!  f1d.1
      call mult("t1gr","s3  ","s2  ",131)
      call mult("s2  ","t1vr","s1  ",311)
!     call reorder("s1  ","s2  ",2134,5,1)
      call reorder("s1  ","s2  ",2134,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","f1a ",3)
!  f2d
!     call reorder("thq ","s1  ",2341,5,1)
      call reorder("thq ","s1  ",2341,5)
      call mult("tvf ","s1  ","s2  ",131)
      call add("tnew",1.d0,"s2  ",-1.d0,"tnew")
      call prt("s2  ","f2d ",3)
!  f3d
      call mult("thqf","t1vr","s2  ",311)
      call unfold_ccsd("s2  ","s1  ",2)
      call add("tnew",1.d0,"s1  ",-1.d0,"tnew")
      call prt("s1  ","-f3d",3)
!  f1e
!     call reorder("thq ","s1  ",2413,5,1)
      call reorder("thq ","s1  ",2413,5)
      call mult("s3  ","s1  ","s2  ",222)
      call add("tb2 ",1.d0,"s2  ",-2.d0,"tb2 ")
      call prt("s2  ","f1e ",3)
!  f1c
!     call reorder("thpf","s2  ",1432,5,1)
      call reorder("thpf","s2  ",1432,5)
      call mult("s2  ","s1  ","s3  ",222)
      call add("tb2 ",1.d0,"s3  ",1.d0,"tb2 ")
      call prt("s3  ","f1c ",3)
!  f1d
!     call reorder("thp ","s1  ",2314,4,1)
      call reorder("thp ","s1  ",2314,4)
!     call reorder("thqf","s3  ",1342,5,1)
      call reorder("thqf","s3  ",1342,5)
      call mult("s3  ","s1  ","s2  ",222)
      call add("tb2 ",1.d0,"s2  ",1.d0,"tb2 ")
      call prt("s2  ","f1d ",3)
      call prtime("dfoldhq ",0,1)
      return
      end
! ********************************************************
      subroutine dfold3_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!  diag 1
      call reord7    ("t3q ","s1  ",231456)
      call reord6    ("s1  ","s2  ",124563)
      call mult6     ("s2  ","tvf ","s1  ",5)
      call reord6    ("s1  ","s2  ",126345)
      call reord7    ("s2  ","s1  ",312456)
      call add       ("s1  ",-0.5d0,"t3nw",1.d0,"t3nw")
!  diag 2
      call reord6    ("t3q ","s1  ",123564)
      call mult6     ("s1  ","tgf ","s2  ",5)
      call reord6    ("s2  ","s1  ",123645)
      call add       ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
!  diags 3,7
      call reord7    ("t3q ","s1  ",231456)
      call reord6    ("s1  ","s3  ",125634)
      call reorder   ("thpf","s1  ",1432,4)
      call reorder   ("thqf","s2  ",1342,4)
      call add       ("s1  ",0.5d0,"s2  ",-1.d0,"s1  ")
      call mult6     ("s1  ","s3  ","s2  ",3)
      call reord7    ("s2  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",124356)
      call add       ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  diag 4
      call reorder   ("tv  ","s1  ",2341,4)
      call mult6     ("thqf","s1  ","s2  ",1)
      call reord6    ("s2  ","s1  ",124356)
      call fold6     ("t3nw","s1  ",idum,2)
!  diag 5
      call reorder   ("thqf","s1  ",1342,4)
      call reorder   ("tg  ","s2  ",2341,4)
      call mult6     ("s1  ","s2  ","s3  ",1)
      call reord6    ("s3  ","s1  ",125346)
      call reord7    ("s1  ","s2  ",132456)
      call reord6    ("s2  ","s1  ",126345)
      call cmult     ("s1  ",-1.d0,"s2  ")
      call fold6     ("t3nw","s2  ",idum,2)
!  diag 6
      call t1t2      ("t1v ","tgr ","t3p ","s1  ",123456,1)
      call t1t2      ("t1gr","tv  ","s1  ","s2  ",123456,2)
      call reord7    ("s2  ","s1  ",231456)
      call reord6    ("s1  ","s2  ",124635)
      call reorder   ("thqf","s1  ",1342,4)
      call mult6     ("s1  ","s2  ","s3  ",3)
      call reord7    ("s3  ","s1  ",132456)
      call reord6    ("s1  ","s2  ",124356)
      call add       ("s2  ",-1.d0,"t3nw",1.d0,"t3nw")
      return
      end
!  *******************************************************************
      subroutine rdbox_chq
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

      call rdfile    (41,0,"boxA",-1)
      call rdfile    (41,1,"boxB",-1)
      call rdfile    (41,2,"boxC",-1)
      call rdfile    (41,3,"boxD",-1)
      call rdfile    (41,4,"boxE",-1)
      call rdfile    (41,5,"boxF",-1)
      call rdfile    (41,6,"boxG",-1)
      call rdfile    (41,7,"boxH",-1)
      call rdfile    (41,8,"boxI",-1)
      if (iopt(17).eq.3) then
        call rdfile    (20,0,"bx1 ",-1)
        call rdfile    (20,1,"bx2 ",-1)
        call rdfile    (20,2,"bx3 ",-1)
        call rdfile    (20,3,"bx4 ",-1)
        call rdfile    (20,4,"bx5 ",-1)
        call rdfile    (20,5,"bx6 ",-1)
        call rdfile    (20,6,"bx7 ",-1)
        call rdfile    (20,7,"bx8 ",-1)
        call rdfile    (42,0,"bx31",-1)
        call rdfile    (42,1,"bx32",-1)
        call rdfile    (42,2,"bx33",-1)
      end if
      return
      end 
!  ******************************************************************
      subroutine wrtbox_chq
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call wrtfile  (42,3,"bx34")
      call wrtfile  (42,4,"bx35")
      call wrtfile  (42,5,"bx36")
      call wrtfile  (42,6,"bx37")
      return
      end 
!  ******************************************************************
      subroutine box34 (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      if (iparm.eq.0) then
!  diag 1
        call reorder   ("vpgp","bx34",1324,4)
!  diag 3
        call reorder   ("tpm ","s1  ",1423,4)
        call mult      ("s1  ","vd11","s2  ",222)
        call add       ("s2  ",-1.d0,"bx34",1.d0,"bx34")
!  diag 4
        call reorder   ("ppgh","s1  ",3241,4)
        call mult      ("tvm2","s1  ","s2  ",222)
        call reorder   ("s2  ","s1  ",1324,4)
        call add       ("s1  ",-1.d0,"bx34",1.d0,"bx34")
!  diag 5
        call reorder   ("vhhp","s1  ",1432,4)
        call mult      ("s1  ","tgm ","s2  ",222)
        call reorder   ("s2  ","s1  ",1342,4)
        call add       ("s1  ",1.d0,"bx34",1.d0,"bx34")
!  diag 7
        call t1t2      ("t1v ","tgr ","t3p ","s1  ",123456,1)
        call t1t2      ("t1gr","tv  ","s1  ","s2  ",123456,2)
        call t1t2      ("t1c ","tpm ","s2  ","s1  ",123456,3)
        call reord6    ("s1  ","s2  ",125346)
        call reord7    ("s2  ","s1  ",132456)
        call reord6    ("s1  ","s2  ",126543)
        call mult6     ("s2  ","va4 ","s1  ",2)
        call add       ("s1  ",1.d0,"bx34",1.d0,"bx34")
!  diag 8
        call reorder   ("vhgp","s1  ",1342,4)
        call mult      ("s1  ","t1cr","s2  ",311)
        call reorder   ("s2  ","s1  ",1243,4)
        call add       ("s1  ",-1.d0,"bx34",1.d0,"bx34")
!  diag 9
        call reorder   ("vphp","s1  ",1243,4)
        call mult      ("s1  ","t1g ","s2  ",311)
        call reorder   ("s2  ","s1  ",1423,4)
        call add       ("s1  ",-1.d0,"bx34",1.d0,"bx34")
      else if (iparm.eq.1) then
!  diag 2
        call add       ("vd11",2.d0,"vd8 ",-1.d0,"s1  ")
        call mult      ("ta3 ","s1  ","s2  ",222)
        call move      ("s2  ","bx34")
!  diag 6
        call cmult     ("t1c ",0.5d0,"s1  ")
        call t1t2      ("s1  ","thq ","t3q ","s2  ",123456,2)
        call t1t2      ("t1c ","thq ","s2  ","s1  ",123456,3)
        call reord6    ("s1  ","s2  ",124356)
        call reord7    ("s2  ","s1  ",132456)
        call reord6    ("s1  ","s2  ",126543)
        call add       ("va4 ",-2.d0,"va8 ",1.d0,"s1  ")
        call mult6     ("s2  ","s1  ","s3  ",2)
        call add       ("s3  ",1.d0,"bx34",1.d0,"bx34")
      end if
      return
      end
!  *******************************************************************
      subroutine box35 (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      if (iparm.eq.0) then
!  diag 1
        call cmult     ("vhgh",-1.d0,"bx35") 
!  diag 3
        call reorder   ("tpm ","s1  ",1423,4)
        call reorder   ("hhhp","s2  ",1324,4)
        call mult      ("s1  ","s2  ","s3  ",222)
        call reorder   ("s3  ","s1  ",1324,4)
        call add       ("s1  ",1.d0,"bx35",1.d0,"bx35")
!  diag 4
        call mult      ("vhhp","tgm2","s1  ",222)
        call reorder   ("s1  ","s2  ",1432,4)
        call add       ("s2  ",1.d0,"bx35",1.d0,"bx35")
!  diag 5
        call reorder   ("ppgh","s1  ",3412,4)
        call mult      ("tvm ","s1  ","s2  ",222)
        call add       ("s2  ",-1.d0,"bx35",1.d0,"bx35")
!  diag 7
        call t1t2      ("t1v ","tgr ","t3p ","s1  ",123456,1)
        call t1t2      ("t1gr","tv  ","s1  ","s2  ",123456,2)
        call t1t2      ("t1c ","tpm ","s2  ","s1  ",123456,3)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125364)
        call mult6     ("s1  ","hhpp","s2  ",2)
        call add       ("s2  ",1.d0,"bx35",1.d0,"bx35")
!  diag 8
        call mult      ("vhgp","t1c ","s1  ",311)
        call reorder   ("s1  ","s2  ",1432,4)
        call add       ("s2  ",-1.d0,"bx35",1.d0,"bx35")
!  diag 9
        call reorder   ("hphg","s1  ",1432,4)
        call mult      ("t1v ","s1  ","s2  ",131)
        call add       ("s2  ",-1.d0,"bx35",1.d0,"bx35")
!  diag 10
        call reorder   ("vhhh","s1  ",1243,4)
        call mult      ("s1  ","t1g ","s2  ",311)
        call reorder   ("s2  ","s1  ",1243,4)
        call add       ("s1  ",1.d0,"bx35",1.d0,"bx35")
      else if (iparm.eq.1) then
!  diag 2
        call reorder   ("hhhp","s1  ",1324,4)
        call add       ("s1  ",-2.d0,"hhhp",1.d0,"s1  ")
        call mult      ("ta3 ","s1  ","s2  ",222)
        call reorder   ("s2  ","s1  ",1324,4)
        call move      ("s1  ","bx35")
!  diag 6
        call cmult     ("t1c ",0.5d0,"s1  ")
        call t1t2      ("s1  ","thq ","t3q ","s2  ",123456,3)
        call t1t2      ("t1c ","thq ","s2  ","s1  ",123456,2)
        call reord6    ("s1  ","s2  ",124356)
        call add       ("hhpp",-2.d0,"va6 ",1.d0,"s1  ")
        call mult6     ("s2  ","s1  ","s3  ",2)
        call add       ("s3  ",1.d0,"bx35",1.d0,"bx35")
      end if
      return
      end
!  *******************************************************************
      subroutine box36 (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      if (iparm.eq.0) then
!  diag 1
        call reorder   ("vphg","s1  ",1432,4)
        call add       ("vgv ",2.d0,"s1  ",-1.d0,"bx36")  
!  diag 2
        call reorder   ("ppgh","s1  ",3421,4)
        call reorder   ("ppgh","s2  ",3412,4)
        call add       ("s1  ",2.d0,"s2  ",-1.d0,"s1  ")  
        call mult      ("t1v ","s1  ","s2  ",131)
        call add       ("s2  ",1.d0,"bx36",1.d0,"bx36")
!  diag 3
        call add       ("vvc2",-2.d0,"vvc3",1.d0,"s1  ")  
        call mult      ("s1  ","t1g ","s2  ",311)
        call reorder   ("s2  ","s1  ",1423,4)
        call add       ("s1  ",1.d0,"bx36",1.d0,"bx36")
!  diag 5
        call reorder   ("tpm ","s1  ",1423,4)
        call add       ("va1 ",-2.d0,"va2 ",1.d0,"s2  ")  
        call mult      ("s1  ","s2  ","s3  ",222)
        call add       ("s3  ",1.d0,"bx36",1.d0,"bx36")
      else if (iparm.eq.1) then
!  diag 4
        call add       ("va1 ",2.d0,"va2 ",-1.d0,"s1  ")  
        call mult      ("ta3 ","s1  ","s2  ",222)
        call cmult     ("s2  ",2.d0,"s1  ")
        call move      ("s1  ","bx36")
      end if
      return
      end 
!  *******************************************************************
      subroutine box37 (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      if (iparm.eq.0) then
!  diag 1
        call cmult     ("vgv ",-1.d0,"bx37")              
!  diag 2
        call reorder   ("ppgh","s1  ",3421,4)
        call mult      ("t1v ","s1  ","s2  ",131)
        call add       ("s2  ",-1.d0,"bx37",1.d0,"bx37")
!  diag 3
        call mult      ("vvc2","t1g ","s2  ",311)
        call reorder   ("s2  ","s1  ",1423,4)
        call add       ("s1  ",1.d0,"bx37",1.d0,"bx37")
!  diag 5
        call reorder   ("tpm ","s1  ",1423,4)
        call mult      ("s1  ","va1 ","s2  ",222)
        call add       ("s2  ",1.d0,"bx37",1.d0,"bx37")
      else if (iparm.eq.1) then
!  diag 4
        call add       ("va1 ",-2.d0,"va2 ",1.d0,"s1  ")  
        call mult      ("ta3 ","s1  ","s2  ",222)
        call move      ("s2  ","bx37")
      end if
      return
      end 
!  *******************************************************************
