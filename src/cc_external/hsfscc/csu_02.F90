!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine csu
!-----------------------------------------------------------------------------
!
!  (0,2) sector
!
!-----------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use cint_mod
      use dat_mod
      use hsfs_parameters
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

      character*256 filenm,filenm51,filenm52,filenm54,filenm6,          &
               filenm20,filenm21,filenm22

!     iopt(1) - print.
!     iopt(2) - no. of iterations.
!     iopt(3) - convergence = 10**(-iopt(3))
!     iopt(4) .gt.0 - do rle every iopt(4) iterations.
!             .lt.0 - diis
!     iopt(5) - print only elements gt 10**(-iopt(5))
!     iopt(6) - start avrg.
!     iopt(8) - fock.
!     iopt(10) .ne.0 - do not save pppp.
      real*4 cputime

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
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')

      write (LUPRI,'(/,2X,A)') '   $$$$$$$$$   csu / (0,2)   $$$$$$$$$  '

      !call check_modules

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
      call prtime("rdin    ",0,0)
      if (iopt(9).ne.0) call rdrsp2
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      irecu=5
      call rdfile(51,0,"tc  ",-1)
      call rdfile(51,1,"t1c ",-1)
      call rdfile(52,0,"tv  ",-1)
      call rdfile(52,1,"t1v ",-1)
      call rdfile(52,2,"tvf ",-1)
      if (iopt(17).ne.0) then
        call rdfile  (52,3,"t31 ",-1)
        if (iopt(17).eq.3) then
          call rdfile (51,2,"t30 ",-1)
!         call rdbox
          call rdbox_csu
        end if
      end if
      call t1sqr("t1c ","t1c ","tc  ","tcm ")
      call t1sqr("t1v ","t1c ","tv  ","tvm ")
      call ordvu
      call calvu0_csu
!      call pdiag
!  pdiag eliminated for lack of io file numbers. to be reinstated
!  if possible, as significant time saving is achieved.
      iter=0
      ieof=0
      if (iopt(12).eq.0) then
        call tuinit_csu
      else 
        call rdfile (54,0,"tu  ",-1)
        call rdfile (54,1,"tuf ",-1)
        if (iopt(17).ne.0) then
          call rdfile (54,2,"t32 ",-1)
        end if
      end if
      call prt("tu  ","tu      ",2)
      call prt("tuf ","tuf     ",2)
      call fmax("tu  ",tmax,itmax)
      write(6,2400) tmax,idec(itmax)
      if (iopt(4).gt.0) then
        if (iopt(17).eq.0) then
          call dbdiag ("dum1","vu0 ","dum3","dum4",0100,1)
        else
          call reord6 ("v30 ","s1  ",123546)
          call reord7 ("s1  ","s2  ",213456)
          call add    ("v30 ",1.d0,"s2  ",1.d0,"s1  ")
          call dbdiag ("dum1","vu0 ","s1  ","dum4",0110,1)
        end if
        call wrtrle ("tu  ",2,1) 
        if (iopt(17).ne.0) then
          call wrtrle ("t32 ",3,1)   
        end if
      endif
      call t1sqr("t1v ","t1v ","tu  ","tm  ")
220   iter=iter+1
      call prtime("strt itr",iter,0)
      call flush (6)
      call wrtfile (54,0,"tu  ")
      call wrtfile (54,1,"tuf ")
      if (iopt(17).ne.0) then
        call wrtfile  (54,2,"t32 ")
      end if
      call ordtua_csu
      if (iopt(17).ne.0) then
        call caltu3_csu
      end if
      call caltu2_csu
      call ordtub_csu
      if (iopt(6).gt.0.and.iter.ge.mod(iopt(6),100)) then
        fnew=max(dfloat(iopt(6)/100),1.d0)
        if (fnew.le.6.5d0) then
          fnew=1.d0/fnew
        else
          fnew=fnew-5.d0
        endif
        call add("tnew",fnew/(1.d0+fnew),"tu  ",1.d0/(1.d0+fnew),"tnew")
        if (iopt(17).ne.0) then
          call add("t3nw",fnew/(1.d0+fnew),"t32 ",1.d0/(1.d0+fnew),"t3nw")
        end if
        write (6,'('' t averaged with fnew '',f7.3)') fnew
      endif
      call fmax("tnew",tmax,itmax)
      t3max = 0.0d0 ! mi initialization
      if (iopt(17).ne.0) then
        call fmax6 ("t32 ",t3max,id3mx,ind3mx)
      end if
      call diffmax("tu  ","tnew","s2  ",diff,idiff)
      df3 = 0.0d0 ! mi initialization
      if (iopt(17).ne.0) then
        call dfmax3 ("t32 ","t3nw","s2  ",df3,iddf3,inddf3)
      end if
      write (6,*) ' '
      write (6,4321) iter,diff,idec(idiff)
4321  format ("iteration # ",i3 , "diffmx= ",e10.3,2x,i12)

      if (iopt(17).ne.0) then
        write (6,2401) df3,idec(iddf3),idec(inddf3)
2401    format (24x,e10.3,2x,i12,2x,i12)
      end if
      write(6,2400) tmax,idec(itmax),t3max,idec(id3mx),idec(ind3mx)
2400  format(' max t,t3=',f15.9,2x,i12,2x,f15.9,2x,i6,2x,i10)
!      abort iteration process
      if (abs(tmax).gt.1.d10.or.abs(tmax).lt.1.d-10) then
        write (6,*) 'tmax=',tmax
        if (tmax.ne.0.0d0) then
          write (6,*) 'csu: abort - tmax out of bounds !'
          call abext('csu: abort - tmax out of bounds !')      
        endif
      endif
! diis
      if (iopt(4).lt.0) then
        call diis(iter,icnvg,01,"dumm","dumm","tnew","tu  ")
        if(iter.gt.1) then
          call diffmax("tu  ","tnew","s2  ",diff,idiff)
          write (6,*) ' '
          write (6,'('' iteration # '',i3,'' diis diffmx='',e10.3,2x,   &
      i12)')  iter,diff,idec(idiff)
        endif
        if(icnvg.eq.1) go to 910
      endif
      call move("tnew","tu  ")
      if (iopt(17).ne.0) then
        call move ("t3nw","t32 ")      
      end if
      call t1sqr("t1v ","t1v ","tu  ","tm  ")
      call prt("tnew","new t   ",2)
      if (max(abs(diff),abs(df3)).lt.convrg) go to 910
      if (iopt(4).gt.0) go to 260
!   no rle
      if (iter.ge.iopt(2)) go to 920
      go to 220
!   with rle
260   continue
      itr1=mod(iter-1,iopt(4))+1
      call wrtrle ("tu  ",2,itr1+1)
      if (iopt(17).ne.0) then
        call wrtrle ("t32 ",3,itr1+1)
      end if
      if (itr1.lt.iopt(4)) go to 220
      if (iopt(17).eq.0) then
        call rle  (itr1,0100,"dum1","tu  ","dum3","dum4",ising)
      else
        call rle  (itr1,0110,"dum1","tu  ","t32 ","dum4",ising)
      end if
!  rle singular
      if (ising.eq.1) then
        if (iopt(17).eq.0) then
          call rle  (itr1-1,0100,"dum1","tu  ","dum3","dum4",ising)
        else
          call rle  (itr1-1,0110,"dum1","tu  ","t32 ","dum4",ising)
        end if
        if (ising.eq.1) then
          iopt(4)=0
          write(*,*) ' ***** rle singular, turned off'
          go to 220
        endif
      endif
      call t1sqr("t1v ","t1v ","tu  ","tm  ")
      call wrtrle ("tu  ",2,1)   
      if (iopt(17).ne.0) then
        call wrtrle ("t32 ",3,1)      
      end if
      if (iter.ge.iopt(2)) go to 920
      go to 220
!
!  different endings
910   write (6,*) ' csu  converged after iter ',iter
      egr0=0.d0
      call eneru2 (egr0,1,3,max(1,iopt(1)))
      call eneru2 (egr0,2,3,max(1,iopt(1)))
      go to 970
920   write (6,*) ' iteration limit'
      call abext('limit   ')     
970   call prt("tu  ","last t  ",0)
      call prt("tuf ","tuf ",0)
      call wrtfile(54,0,"tu  ")
      call wrtfile(54,1,"tuf ")
      if (iopt(17).ne.0) then
        call wrtfile  (54,2,"t32 ")
        if (iopt(17).eq.3) then
          call wrtbox_csu
        end if
      end if

      write (*,*) nstrt(1,nflcur),' words used in igen,',nstrt(2,nflcur),' words used in gen'
      call prtime("end     ",iter,0)

      end
! *********************************************
      subroutine ordvu
!-----------------------------------------------------------------------------
!
! va: hhpp
! vd: hppp
! vva: vhpp
! vvc: vhhp
!              hhpp 1324 = va1
!              hhpp 1423 = va2
!              hhpp 3412 = va4
!              hhpp 1243 = va6
!              hhpp 3421 = va8
!              hppp 2314 = vd7
!              hppp 2413 = vd11
!              hppp 2134 = vd12
!              hppp 2143 = vd14
!              vhpp 1324 = vva1
!              vhhp 1324 = vvc1
!
!-----------------------------------------------------------------------------
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
#include "param_ccsd.h"
!     common /cint/nint,idum1,vint(ndint),iint(ndint)
#include "cint.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)
#include "fock.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif
!
      call readv    (nint,vint,iint,ndint,nfock,fock,ifock,ndt1)
      call prtime   ("readv   ",nint,0)
      call seltyp   (2233,0,"hhpp",4)
      call seltyp   (2333,0,"hppp",4)
      call seltyp4  (3333,0,"pppp",4)
      call seltyp   (3233,1000,"vhpp",4)
      call seltyp   (3323,1000,"vphp",4)
      call seltyp   (3223,1000,"vhhp",4)
      call seltyp   (3333,1100,"vvpp",4)
      call seltyp   (3322,1100,"vvhh",4)
!uk temporary
      call seltyp   (3333,1111,"vvvv",0)
!uk
      call reorder  ("hhpp","va1 ",1324,4)
      call reorder  ("hhpp","va2 ",1423,4)
      call reorder  ("hhpp","va4 ",3412,4)
      call reorder  ("hhpp","va6 ",1243,4)
      call reorder  ("hhpp","va8 ",3421,4)
      call reorder  ("hppp","vd7 ",2314,4)
      call reorder  ("hppp","vd11",2413,4)
      call reorder  ("hppp","vd12",2134,4)
      call reorder  ("hppp","vd14",2143,4)
      call reorder  ("vhpp","vva1",1324,4)
      call reorder  ("vhhp","vvc1",1324,4)
      call selt1    (33,0,"pp  ",4)
      call selt1    (23,0,"hp  ",4)
      call selt1    (32,10,"vh  ",4)
      if (iopt(17).ne.0) then
        call seltyp   (3333,1000,"vppp",4)
        call seltyp   (3232,1000,"vhph",4)
        call seltyp   (3332,1100,"vvph",4)
        call seltyp   (3333,1111,"vvvv",4)
      end if
      call prtime("ordu",0,0)
      return
      end
!**************************************************************
!     subroutine tuinit
      subroutine tuinit_csu
!-----------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------
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

      if (iopt(20).ge.1) then
          call rdfile (54,0,"tu  ",-1)
          call rdfile (54,1,"tuf ",-1)
          if ((iopt(17).ne.0).and.(iopt(20).eq.1))                      &
        call rdfile  (54,2,"t32 ",-1)
      else
       call move("vu0 ","tu  ")
!      call fold("tu  ","tuf ")
       call fold_ccsd("tu  ","tuf ")
       call diveps("tu  ")
      end if
      if ((iopt(17).ne.0).and.(iopt(20).ne.1)) then
        call move    ("v30 ","t32 ")
        call reord6  ("t32 ","s1  ",123546)
        call reord7  ("s1  ","s2  ",213456)
        call add     ("t32 ",1.d0,"s2  ",1.d0,"t32 ")
        call dveps6  ("t32 ")
      end if
      return
      end
! **************************************************
      subroutine ordta0_csu
!-----------------------------------------------------------------------------
!  tv1 = tv 1324
!  tvm1 = tvm 1324
!  tv2 = tv 1423
!  tvm2 = tvm 1423
!  tv3 = tv 1432
!  tcm4 = tcm 3412
!  t1cor = t1c 21
!  tb1 = tnew 1324
!  tb4 = tnew 3412
!  tb15 = tnew 1342
!-----------------------------------------------------------------------------
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call reorder("tc  ","tc1 ",1324,4)
      call reorder("tc  ","tc3 ",1423,4)
      call reorder("tc  ","tc4 ",3412,4)
      call reorder("tv  ","tv1 ",1324,4)
      call reorder("tv  ","tv2 ",1423,4)
      call reorder("tv  ","tv3 ",1432,4)
      call reorder("tv  ","tv6 ",1243,4)
      call reorder("tvm ","tvm1",1324,4)
      call reorder("tvm ","tvm2",1423,4)
      call reorder("tcm ","tcm2",1423,4)
      call reorder("tcm ","tcm4",3412,4)
      call reorder("t1c ","t1cr",1243,4)
      call reorder("t1v ","t1vr",1243,4)
      call clear("tb1 ",3333,1010)
      call clear("tb4 ",3333,11)
      call clear("tb15",3333,1010)
      call clear("tuf ",3333,1111)
      call prtime("ordta0  ",0,1)
      return
      end
! ******************************************************
      subroutine calvu0_csu
!-----------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------
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

      call ordta0_csu
      call move("vvpp","vu0 ")
      call cmult("vu0 ",0.5d0,"vu0 ")
      call diag3_csu
      call diag4_csu
      call diag5_csu
      call diag6_csu
      call diag8_csu
      call diag9_csu
      call diag11_csu
      call diag12_csu
      call diag13_csu
      call dgt1p3_csu
      call dgt1p4_csu
      call diag21_csu
      call diag26_csu
!     call diag28
      call diag28_csu
      if (iopt(17).ne.0) then
!  diag 31
        call reord6    ("t31 ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125436)
        call add       ("vhhp",-2.d0,"vvc1",1.d0,"s2  ")
        call mult6     ("s2  ","s1  ","s3  ",2)
        call add       ("s3  ",1.d0,"vu0 ",1.d0,"vu0 ")
!  diag 33
        call reord6    ("t31 ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126345)
        call mult6     ("s1  ","vhhp","s2  ",2)
        call reorder   ("s2  ","s1  ",1423,4)
        call add       ("s1  ",1.d0,"vu0 ",1.d0,"vu0 ")
      end if
      if (iopt(17).ge.2) then
!    diag 36,38
        call reord6    ("t31 ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125634)
        call add       ("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",2)
        call reord6    ("t31 ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126543)
        call mult6     ("s1  ","va8 ","s2  ",2)
        call add       ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
        call mult      ("s3  ","t1v ","s1  ",311)
        call reorder   ("s1  ","s2  ",1423,4)
        call add       ("s2  ",1.d0,"vu0 ",1.d0,"vu0 ")
      end if
      call ordtb0_csu
      call sqr("vu0 ","s2  ")
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
!  dg3n1b
        call mult6      ("tv  ","vppp","s1  ",1)
        call reord6     ("s1  ","s2  ",124365)
        call reord7     ("s2  ","v30 ",132456)
!  dg3n1c
        call mult6      ("tv6 ","vppp","s1  ",1)
        call reord6     ("s1  ","s2  ",124563)
        call reord7     ("s2  ","s1  ",132456)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
!  dg3n2a
        call reorder    ("vhhp","s1  ",1243,4)
        call mult6      ("tv3 ","s1  ","s2  ",1)
        call reord6     ("s2  ","s1  ",124356)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",125436)
        call add        ("s1  ",-1.d0,"v30 ",1.d0,"v30 ")
!  dg3n2b
        call mult6      ("tv3 ","vhph","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125463)
        call add        ("s2  ",-1.d0,"v30 ",1.d0,"v30 ")
!  dg3n2c,f1
        call reorder    ("tc  ","s1  ",1432,4)
        call mult6      ("vvph","s1  ","s3  ",1)
        call reorder    ("tv  ","s1  ",2341,4)
        call mult6      ("vvvv","s1  ","s2  ",1)
        call cmult      ("s2  ",-1.d0,"s1  ")
        call fold6      ("s3  ","s1  ",idum,2)
        call reord6     ("s3  ","s1  ",124356)
        call add        ("s1  ",-1.d0,"v30 ",1.d0,"v30 ")
      else if (iopt(17).eq.3) then
!  GROUP I
!   diag 2
        call mult6      ("tv  ","bx12","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call reord7     ("s2  ","v30 ",132456)
!   diag 3
        call mult6      ("tv6 ","bx12","s1  ",1)
        call reord6     ("s1  ","s2  ",124653)
        call reord7     ("s2  ","s1  ",132456)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP II
!   diag 1
        call mult6      ("tv3 ","bx13","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125436)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!   diag 2
        call mult6      ("tv3 ","bx14","s1  ",1)
        call reord6     ("s1  ","s2  ",124356)
        call reord7     ("s2  ","s1  ",132456)
        call reord6     ("s1  ","s2  ",125463)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!   diag 3
        call box24_csu  (0)
        call mult6      ("bx24","tc4 ","s1  ",1)
        call reord6     ("s1  ","s2  ",126354)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP III
!   diag 2
        call reord6     ("t31 ","s1  ",124536)
        call mult6      ("s1  ","bx15","s2  ",3)
        call reord6     ("s2  ","s1  ",125364)
        call reord7     ("s1  ","s2  ",132456)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP IV
!   diag 2
        call reord6     ("t31 ","s1  ",124635)
        call mult6      ("s1  ","bx16","s2  ",3)
        call reord6     ("s2  ","s1  ",125364)
        call reord7     ("s1  ","s2  ",132456)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!   diag 3
        call reord6     ("t31 ","s1  ",125634)
        call mult6      ("s1  ","bx16","s2  ",3)
        call reord6     ("s2  ","s1  ",125463)
        call reord7     ("s1  ","s2  ",132456)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP V
!   diag 2
        call reord6     ("t31 ","s1  ",124635)
        call mult6      ("s1  ","bx17","s2  ",3)
        call reord6     ("s2  ","s1  ",125346)
        call reord7     ("s1  ","s2  ",132456)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!   diag 3
        call reord6     ("t31 ","s1  ",125634)
        call mult6      ("s1  ","bx17","s2  ",3)
        call reord6     ("s2  ","s1  ",125643)
        call reord7     ("s1  ","s2  ",132456)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP VI
!   diag 1
        call reord6     ("t31 ","s1  ",124356)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",125634)
        call mult6      ("s1  ","bx18","s2  ",3)
        call reord6     ("s2  ","s1  ",125346)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",126345)
        call add        ("s1  ",1.d0,"v30 ",1.d0,"v30 ")
!   diag 2
        call box25_csu  (0)
        call reord6     ("t30 ","s1  ",126345)
        call reord7     ("s1  ","s2  ",132456)
        call reord6     ("s2  ","s1  ",126543)
        call mult6      ("bx25","s1  ","s2  ",3)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP IX
!   diag 3
        call reord6     ("t31 ","s1  ",124563)
        call mult6      ("s1  ","bx19","s2  ",5)
        call reord6     ("s2  ","s1  ",126354)
        call reord7     ("s1  ","s2  ",132456)
        call add        ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
      end if
      return
      end
!  ********************************************************************
      subroutine ordtb0_csu

!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call reorder("tb1 ","s2  ",1324,4)
      call add("vu0 ",1.d0,"s2  ",1.d0,"vu0 ")
      call reorder("tb4 ","s2  ",3412,4)
      call add("vu0 ",1.d0,"s2  ",1.d0,"vu0 ")
      call reorder("tb15","s2  ",1342,4)
      call add("vu0 ",1.d0,"s2  ",1.d0,"vu0 ")
      call prtime("ordtb0  ",0,1)
      return
      end
! *****************************************************
      subroutine ordtua_csu

!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call reorder("tu  ","ta3 ",1432,4)
      call reorder("tm  ","tm4 ",3412,4)
      call reorder("tm  ","tm6 ",1243,4)
      call clear("tb1 ",3333,1010)
      call clear("tb2 ",3333,1001)
      call clear("tb4 ",3333,11)
      call clear("tb6 ",3333,1100)
      call clear("tb15",3333,1010)
      call prtime("ordta   ",0,1)
      return
      end
! ****************************************************
      subroutine caltu3_csu
!-----------------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------------
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

      call move          ("v30 ","t3nw")
      if (iopt(17).lt.3) then
!  dg3n1a
        call mult6       ("tu  ","hppp","s1  ",1)
        call reord6      ("s1  ","s2  ",124356)
        call add         ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
      else if (iopt(17).eq.3) then
!  GROUP I
!   diag 1
        call mult6       ("tu  ","bx1 ","s1  ",1)
        call reord6      ("s1  ","s2  ",124365)
        call add         ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP II
!   diag 3
        call box24_csu   (1)
        call mult6       ("bx24","tc4 ","s1  ",1)
        call reord6      ("s1  ","s2  ",126354)
        call add         ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP III
!   diag 1
        call reord6      ("t32 ","s1  ",124536)
        call mult6       ("s1  ","bx3 ","s2  ",3)
        call reord6      ("s2  ","s1  ",125346)
        call add         ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
!  GROUP IV
!   diag 1
        call reord6      ("t32 ","s1  ",124635)
        call mult6       ("s1  ","bx4 ","s2  ",3)
        call reord6      ("s2  ","s1  ",125346)
        call add         ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP V
!   diag 1
        call reord6      ("t32 ","s1  ",124635)
        call mult6       ("s1  ","bx5 ","s2  ",3)
        call reord6      ("s2  ","s1  ",125364)
        call add         ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP VI
!   diag 2
        call box25_csu   (1)
        call reord6      ("t30 ","s1  ",126345)
        call reord7      ("s1  ","s2  ",132565)
        call reord6      ("s2  ","s1  ",126543)
        call mult6       ("bx25","s1  ","s2  ",3)
        call add         ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP VII
!   diag 1a
        call t1t2        ("t1v ","tv  ","t32 ","s1  ",123456,2)
        call t1t2        ("t1c ","tu  ","s1  ","s3  ",123456,3)
        call mult6       ("s3  ","pppp","s1  ",3)
        call add         ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
        call reord7      ("s3  ","s1  ",132456)
        call reord6      ("s1  ","s2  ",123465)
        call mult6       ("s2  ","pppp","s1  ",3)
        call reord7      ("s1  ","s2  ",132456)
        call reord6      ("s2  ","s1  ",123465)
        call add         ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 1b
        call t1t2        ("t1v ","tv  ","t32 ","s1  ",123456,1)
        call t1t2        ("t1v ","tv  ","s1  ","s2  ",123456,2)
        call reord6      ("s2  ","s1  ",123645)
        call mult6       ("s1  ","pppp","s2  ",3)
        call reord6      ("s2  ","s1  ",123564)
        call add         ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
        call reord6      ("s2  ","s1  ",123654)
        call reord7      ("s1  ","s2  ",213654)
        call add         ("s2  ",0.5d0,"t3nw",1.d0,"t3nw")
!   diag 2a
        call mult6       ("t32 ","vd12","s1  ",3)
        call mult6       ("s1  ","t1cr","s2  ",5)
        call reord6      ("s2  ","s1  ",123465)
        call add         ("s1  ",-1.d0,"t3nw",1.d0,"t3nw")
!   diag 2b
        call mult6       ("t32 ","vd14","s1  ",3)
        call mult6       ("s1  ","t1cr","s2  ",5)
        call add         ("s2  ",-1.d0,"t3nw",1.d0,"t3nw")
!   diag 2c
        call reord6      ("t32 ","s1  ",123654)
        call mult6       ("s1  ","vd12","s2  ",3)
        call mult6       ("s2  ","t1cr","s1  ",5)
        call reord6      ("s1  ","s2  ",123564)
        call add         ("s2  ",-1.d0,"t3nw",1.d0,"t3nw")
!   diag 3a
        call mult6       ("t32 ","hhpp","s1  ",3)
        call mult6       ("s1  ","tcm4","s2  ",3)
        call add         ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 3b
        call reord6      ("t32 ","s1  ",123645)
        call mult6       ("s1  ","hhpp","s2  ",3)
        call mult6       ("s2  ","tcm4","s1  ",3)
        call reord6      ("s1  ","s2  ",123564)
        call add         ("s2  ",0.5d0,"t3nw",1.d0,"t3nw")
!  GROUP VIII
!   diag 1
        call mult6       ("t32 ","bx7 ","s1  ",5)
        call add         ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
!   diag 2
        call reord6      ("t32 ","s1  ",123564)
        call mult6       ("s1  ","bx7 ","s2  ",5)
        call reord6      ("s2  ","s1  ",123645)
        call add         ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP IX
!   diag 1
        call reord6      ("t32 ","s1  ",124563)
        call mult6       ("s1  ","bx8 ","s2  ",5)
        call reord6      ("s2  ","s1  ",126345)
        call add         ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
      end if
      call prtime      ("caltu3  ",0,1)
      return
      end
!  *********************************************************************
      subroutine caltu2_csu
!-----------------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------------
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

      call move("vu0 ","tnew")
      call cmult("tnew",0.5d0,"tnew")
      call diag2_csu
      call diag7_csu
      call diag14_csu
      call diag15_csu
!     call diag27
      call diag27_csu
      if (iopt(17).ne.0) then
!  diag 30
        call reord6  ("t32 ","s1  ",124356)
        call add     ("vd14",2.d0,"vd12",-1.d0,"s2  ")
        call mult6   ("s1  ","s2  ","s3  ",2)
        call add     ("s3  ",1.d0,"tnew",1.d0,"tnew")
!  diag 32
        call reord6  ("t32 ","s1  ",126345)
        call mult6   ("s1  ","vd14","s2  ",2)
        call add     ("s2  ",-1.d0,"tb6 ",1.d0,"tb6 ")
      end if
      if (iopt(17).ge.2) then
!    diag 34,35      
        call add       ("va1 ",2.d0,"va2 ",-1.d0,"s1  ")
        call mult      ("s1  ","t1c ","s3  ",202)
! (Gal) Feb 95
        call add       ("s3  ",1.d0,"hp  ",1.d0,"s3  ")
        call reord6    ("t32 ","s1  ",124536)
        call reord6    ("t32 ","s2  ",124635)
        call add       ("s1  ",1.d0,"s2  ",-1.d0,"s1  ")
        call mult6     ("s1  ","s3  ","s2  ",0)
        call add       ("s2  ",1.d0,"tnew",1.d0,"tnew")
!    diag 37,39
        call reord6    ("t32 ","s1  ",124356)
        call add       ("hhpp",-2.d0,"va6 ",1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",2)
        call reord6    ("t32 ","s1  ",126354)
        call mult6     ("s1  ","hhpp","s2  ",2)
        call add       ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
        call mult      ("s3  ","t1cr","s1  ",311)
        call add       ("s1  ",1.d0,"tnew",1.d0,"tnew")
      end if
      call prtime  ("caltu2  ",0,1)
      return
      end
! ******************************************************
      subroutine ordtub_csu
!-----------------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------------
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

      call reorder("tb1 ","s2  ",1324,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call reorder("tb2 ","s2  ",1423,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call reorder("tb4 ","s2  ",3412,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call reorder("tb6 ","s2  ",1243,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call reorder("tb15","s2  ",1342,5)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call sqr("tnew","s2  ")
      call fold_ccsd("tnew","tuf ")
!      call heff2
      call dfoldu_csu
      call diveps("tnew")
      if (iopt(17).ne.0) then
        call reord6   ("t3nw","s1  ",123546)
        call reord7   ("s1  ","s2  ",213456)
        call add      ("t3nw",1.d0,"s2  ",1.d0,"t3nw")
        call dveps6   ("t3nw")
      end if
      call prtime("ordtub  ",0,1)
      return
      end
! *************************************************
      subroutine diag7_csu

!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call mult("pppp","tm  ","s2  ",222)
      call sqr("s2  ","s3  ")
      call add("tb4 ",1.d0,"s2  ",0.5d0,"tb4 ")
      call prt("s2  ","diag7   ",3)
      call prtime("diag7   ",0,1)
      return
      end
! **********************************************
      subroutine diag8_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vvhh","tcm4","s2  ",222)
      call add("vu0 ",1.d0,"s2  ",0.5d0,"vu0 ")
      call prt("s2  ","2diag8  ",3)
      call prtime("diag8   ",0,1)
      return
      end
! *************************************************
      subroutine diag3_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vh  ","tv3 ","s1  ",131)
      call add("vu0 ",1.d0,"s1  ",-1.d0,"vu0 ")
      call prt("s1  ","-diag3  ",3)
      call prtime("diga3   ",0,1)
      return
      end
! **************************************************
      subroutine diag2_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tu  ","pp  ","s2  ",311)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","diag2   ",3)
      call prtime("diga2   ",0,1)
      return
      end
! *******************************************************
      subroutine diag4_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add("vva1",2.d0,"vphp",-1.d0,"s2  ")
      call mult("s2  ","tv1 ","s1  ",222)
      call add("tb1 ",1.d0,"s1  ",1.d0,"tb1 ")
      call prt("s1  ","diag4   ",3)
      call prtime("diag4   ",0,1)
      return
      end
! ************************************************************
      subroutine diag5_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vva1","tvm2","s1  ",222)
      call add("tb1 ",1.d0,"s1  ",-1.d0,"tb1 ")
      call prt("s1  ","-diag5  ",3)
      call prtime("diag5   ",0,1)
      return
      end
! ************************************************************
      subroutine diag6_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vphp","tvm2","s1  ",222)
      call add("tb15",1.d0,"s1  ",-1.d0,"tb15")
      call prt("s1  ","-diag6  ",3)
      call prtime("diag6   ",0,1)
      return
      end
! *******************************************************************
      subroutine diag9_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 9 10 20 23
      call add("tv1 ",-1.d0,"tvm2",1.d0,"s2  ")
      call add("va1 ",-2.d0,"va2 ",1.d0,"s3  ")
      call mult("s2  ","s3  ","s1  ",222)
!   diag 20
      call reorder("vhhp","s2  ",1243,5)
      call reorder("vhhp","s3  ",1342,5)
      call add("s2  ",-2.d0,"s3  ",1.d0,"s2  ")
      call mult("s2  ","t1cr","s3  ",311)
      call reorder("s3  ","s2  ",1423,5)
      call add("s1  ",1.d0,"s2  ",1.d0,"s1  ")
!   diag 23
      call add("vd12",2.d0,"vd14",-1.d0,"s2  ")
      call mult("t1v ","s2  ","s3  ",131)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
! -------------
      call mult("tv1 ","s1  ","s2  ",222)
      call add("tb1 ",1.d0,"s2  ",1.d0,"tb1 ")
      call prt("s2  ","diag9   ",3)
      call prtime("diag9   ",0,1)
      return
      end
! **************************************************
      subroutine diag11_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 11 24
      call mult("tvm2","va1 ","s1  ",222)
!   diag 24
      call mult("t1v ","vd12","s2  ",131)
      call add("s1  ",0.5d0,"s2  ",-1.d0,"s1  ")
! -----------------
      call mult("s1  ","tvm2","s2  ",222)
      call add("tb1 ",1.d0,"s2  ",1.d0,"tb1 ")
      call prt("s2  ","diag11  ",3)
      call prtime("diag11  ",0,1)
      return
      end
! **************************************************
      subroutine diag12_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 12 25
      call mult("tvm2","va2 ","s1  ",222)
!   diag 25
      call mult("t1v ","vd14","s2  ",131)
      call add("s1  ",0.5d0,"s2  ",-1.d0,"s1  ")
! ------------------
      call mult("s1  ","tvm2","s2  ",222)
      call add("tb15",1.d0,"s2  ",1.d0,"tb15")
      call prt("s2  ","diag12  ",3)
      call prtime("diag12  ",0,1)
      return
      end
! ************************************************
      subroutine diag13_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 13 18
      call add("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult("tvm1","s2  ","s1  ",113)
!   diag 18
      call add("vvc1",-2.d0,"vhhp",1.d0,"s2  ")
      call mult("s2  ","t1c ","s3  ",202)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
! -----------
      call mult("s1  ","tv3 ","s3  ",131)
      call add("vu0 ",1.d0,"s3  ",1.d0,"vu0 ")
      call prt("s3  ","diag13  ",3)
      call prtime("diag13  ",0,1)
      return
      end
! ****************************************************
      subroutine diag14_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 14 19
      call add("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
      call mult("tcm4","s2  ","s1  ",113)
!   diag 19
      call add("vd11",2.d0,"vd7 ",-1.d0,"s3  ")
      call mult("s3  ","t1c ","s2  ",202)
      call add("s1  ",1.d0,"s2  ",1.d0,"s1  ")
! --- pdig14 to here
      call mult("tu  ","s1  ","s2  ",311)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","diag14  ",3)
      call prtime("diag14  ",0,1)
      return
      end
! ****************************************************
      subroutine diag15_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tm  ","hhpp","s1  ",222)
      call mult("s1  ","tcm4","s2  ",222)
      call add("tnew",1.d0,"s2  ",0.5d0,"tnew")
      call prt("s2  ","2diag15 ",3)
      call prtime("diag15  ",0,1)
      return
      end
! ****************************************************
      subroutine diag21_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 21 17 22
      call mult("tv2 ","vvc1","s1  ",222)
!   diag 17
      call seltyp(3323,1100,"s2  ",5)
      call reorder("s2  ","s3  ",2413,5)
      call add("s1  ",1.d0,"s3  ",-1.d0,"s1  ")
!   diag 22
      call mult("tv2 ","vhhp","s2  ",222)
      call reorder("s2  ","s3  ",3214,5)
      call add("s1  ",1.d0,"s3  ",1.d0,"s1  ")
! --------------
      call mult("s1  ","t1cr","s2  ",311)
      call add("tb1 ",1.d0,"s2  ",1.d0,"tb1 ")
      call prt("s2  ","diag21  ",3)
      call prtime("diag21  ",0,1)
      return
      end
! ************************************************
!     subroutine diag27
      subroutine diag27_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 27 29
!   diag 29
      call mult("tu  ","hp  ","s2  ",311)
!   diag 27
      call mult("tm6 ","hppp","s1  ",222)
      call reorder("s1  ","s3  ",1243,5)
      call add("s2  ",1.d0,"s3  ",1.d0,"s2  ")
! ---------------
      call mult("s2  ","t1cr","s1  ",311)
      call add("tnew",1.d0,"s1  ",-1.d0,"tnew")
      call prt("s1  ","-diag27 ",3)
      call prtime("diag27  ",0,1)
      return
      end
! **********************************************************
!     subroutine diag28
      subroutine diag28_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("t1v ","hp  ","s2  ",111)
      call reorder("tv  ","s3  ",1432,5)
      call mult("s2  ","s3  ","s1  ",131)
      call add("vu0 ",1.d0,"s1  ",-1.d0,"vu0 ")
      call prt("s1  ","-diag28 ",3)
      call prtime("diag28  ",0,1)
      return
      end
! *********************************************************
      subroutine diag26_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 26 16
      !call reorder("tcm ","s1  ",3412,5,0)
      call reorder("tcm ","s1  ",3412,5)
      !call reorder("vhhp","s2  ",1432,5,0)
      call reorder("vhhp","s2  ",1432,5)
      call mult("s1  ","s2  ","s3  ",222)
!   diag 16
      call seltyp(3333,10,"s2  ",5)
      call add("s3  ",1.d0,"s2  ",1.d0,"s3  ")
! -------------
      call mult("s3  ","t1v ","s2  ",311)
      call add("tb4 ",1.d0,"s2  ",1.d0,"tb4 ")
      call prt("s2  ","diag26  ",3)
      call prtime("diag26  ",0,1)
      return
      end
! *********************************************************
      subroutine dgt1p3_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   compensates for overcounting
      call mult("t1v ","vd14","s2  ",131)
      call mult("t1v ","s2  ","s3  ",131)
      call mult("s3  ","t1cr","s2  ",311)
      call add("vu0 ",1.d0,"s2  ",2.d0,"vu0 ")
      call prt("s2  ",".5t1**3 ",3)
      call prtime("dgt1p3  ",0,1)
      return
      end
! ****************************************************
      subroutine dgt1p4_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   -t1**4 compensates for overcounting
      call mult("t1v ","hhpp","s3  ",131)
      call mult("t1v ","s3  ","s2  ",131)
      call mult("t1cr","s2  ","s3  ",131)
      call mult("t1cr","s3  ","s2  ",131)
      call add("tb4 ",1.d0,"s2  ",-1.d0,"tb4 ")
      call prt("s2  ","-dgt1p4 ",3)
      call prtime("dgt1p4  ",0,1)
      return
      end
! *************************************************
      subroutine dfoldu_csu
!-----------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------
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

      call mult("tvf ","ta3 ","s2  ",131)
! --------------
      call mult("tuf ","t1vr","s3  ",311)
!     call unfold("s3  ","s1  ",3)
      call unfold_ccsd("s3  ","s1  ",3)
      call add("s2  ",1.d0,"s1  ",1.d0,"s2  ")
! --------------
      call mult("tuf ","tm4 ","s3  ",222)
      call add("s2  ",1.d0,"s3  ",0.5d0,"s2  ")
! =============
      call sqr("s2  ","s3  ")
!     if (iopt(9).ne.0) call fold("s2  ","s3  ")
      if (iopt(9).ne.0) call fold_ccsd("s2  ","s3  ")
      call add("tnew",1.d0,"s2  ",-1.d0,"tnew")
      call prt("tnew","tnew",3)
      if (iopt(17).eq.3) then
!  diag f1
        call t1t2      ("t1v ","tv  ","t32 ","s1  ",123456,1)
        call t1t2      ("t1v ","tv  ","s1  ","s2  ",123456,2)
        call reord7    ("s2  ","s1  ",312456)
        call reord6    ("s1  ","s2  ",124356)
        call reord7    ("s2  ","s1  ",132456)
        call reord6    ("s1  ","s2  ",125634)
        call mult6     ("tuf ","s2  ","s1  ",3)
        call add       ("s1  ",-0.5d0,"t3nw",1.d0,"t3nw")
!  diag f2
        call reorder   ("tv  ","s1  ",2341,4)
        call mult6     ("tuf ","s1  ","s2  ",1)
        call reord6    ("s2  ","s1  ",124356)
        call fold6     ("t3nw","s1  ",idum,2)
!  diag f3
        call reord7    ("t32 ","s1  ",231456)
        call reord6    ("s1  ","s2  ",124563)
        call mult6     ("s2  ","tvf ","s1  ",5)
        call reord6    ("s1  ","s2  ",126345)
        call reord7    ("s2  ","s1  ",312456)
        call add       ("s1  ",-1.d0,"t3nw",1.d0,"t3nw")
      end if
      call prtime("-dfoldu ",0,1)
      return
      end
! ***********************************************
!      subroutine heff2
!      implicit real*8 (a-h,o-z)
!      include 'parameters.f'
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
!      common // gen(nbig)
!      common /i/ igen(nbig)
!      dimension a(600,600),ipvt(600),z(600),iord(600)
!      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
!      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
!      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
!      common /rsp2/ nrsp2,irsp2(500)
!      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
!      if (iopt(9).eq.0) return
!      call prt("tuf ","tufbef  ",2)
!      call mult("tvf ","ta3 ","s2  ",131)
!      call fold("s2  ","s3  ")
!      call add("s3  ",-1.d0,"tuf ",1.d0,"s3  ")
!      call files("tuf ",ntapef,new,1,-1,-1)
!      call files("s3  ",ntape3,new,3,-1,-1)
!      if (new.ne.0) then
!        write (6,*) ' tuf not found in heff2'
!        stop                   
!      endif
!      nsize=ije1(nij1)
!      istrt11=nstrt(1,ntapef-1)+2*nij1+4
!      istrt12=nstrt(2,ntapef-1)
!      istrt31=nstrt(1,ntape3-1)+2*nij3+4
!      istrt32=nstrt(2,ntape3-1)
!      call dimx(nsize,600,"heff2   ")
!      call dcopy(600*600,0.d0,0,a,1)
!      do 20 i=1,nsize
!20    a(i,i)=1.d0
!      call files("t1v ",ntapev,new,2,-1,-1)
!      if (new.ne.0) then
!        write (6,*) ' t1v not found in heff2'
!        stop
!      endif
!      nsizev=ije2(nij2)
!      istrtv1=nstrt(1,ntapev-1)+2*nij2+4
!      istrtv2=nstrt(2,ntapev-1)
!      do 30 i=1,nsize
!      do 30 j=1,nsize
!      if (ishft(igen(istrt11+i),-8).eq.
!     1 ishft(igen(istrt11+j),-8)) then
!        itt=ior(ishft(iand(igen(istrt11+j),'ff'X),8),
!     1  iand(igen(istrt11+i),'ff'X))
!        a(i,j)=a(i,j)+find(itt,gen(istrtv2+1),igen(istrtv1+1),
!     1  nsizev)
!      endif
!30    continue
!      do 40 i=1,nsize
!      do 40 j=1,nsize
!      if (iand(igen(istrt11+i),'ffff00ff'X).eq.
!     1 iand(igen(istrt11+j),'ffff00ff'X)) then
!        itt=ior(ishft(iand(ishft(igen(istrt11+j),-8),'ff'X),8),
!     1  (iand(ishft(igen(istrt11+i),-8),'ff'X)))
!        a(i,j)=a(i,j)+find(itt,gen(istrtv2+1),igen(istrtv1+1),
!     1  nsizev)
!      endif
!40    continue
!      call files("tm  ",ntapem,new,2,-1,-1)
!      if (new.ne.0) then
!        write (6,*) ' tm file not found in heff2'
!        stop                     
!      endif
!      nsizev=ije2(nij2)
!      istrtv1=nstrt(1,ntapem-1)+2*nij2+4
!      istrtv2=nstrt(2,ntapem-1)
!      do 50 i=1,nsize
!      do 50 j=1,nsize
!      if (ishft(igen(istrt11+i),-16).eq.
!     1 ishft(igen(istrt11+j),-16)) then
!        itt=ior(ishft(iand(igen(istrt11+j),'ffff'X),16),
!     1  iand(igen(istrt11+i),'ffff'X))
!        a(i,j)=a(i,j)+find(itt,gen(istrtv2+1),igen(istrtv1+1),
!     1  nsizev)
!      endif
!50    continue
!c   delete transition to or from r space.
!      io=0
!      do 60 i=1,nsize
!      do 55 ir=1,nrsp2
!      if (iand(igen(istrt11+i),'ffff'X).eq.irsp2(ir)) go to 60
!      if (iand(ishft(igen(istrt11+i),-16),'ffff'X).eq.
!     1 irsp2(ir)) go to 60
!55    continue
!      io=io+1
!      iord(io)=i
!60    continue
!      nsizeo=io
!      do 70 io=1,nsizeo
!      gen(istrt32+io)=gen(istrt32+iord(io))
!      do 70 jo=1,nsizeo
!70    a(io,jo)=a(iord(io),iord(jo))
!      call dgeco(a,600,nsizeo,ipvt,rcond,z)
!      if (rcond+1.0d0.eq.1.0d0) then
!        write (6,*) ' singular lin eqs in heff2'
!        stop                     
!      endif
!      call dgesl(a,600,nsizeo,ipvt,gen(istrt32+1),0)
!c   expand back to original tf
!      do 75 i=1,nsizeo
!75    z(i)=gen(istrt32+i)
!      call dcopy(nsize,0.d0,0,gen(istrt12+1),1)
!      do 80 io=1,nsizeo
!80    gen(istrt12+iord(io))=z(io)
!      call prt("tuf ","newtuf  ",2)
!      return
!      end
!  ********************************************************************
      subroutine box24_csu (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      if (iparm.eq.0) then
!  diag 1
        call cmult    ("vvph",-1.d0,"bx24")
!  diag 2
        call add      ("vvc1",-2.d0,"vhhp",1.d0,"s1  ")
        call mult     ("tv1 ","s1  ","s3  ",222)
!  diag 3
        call mult     ("tvm2","vvc1","s1  ",222)
        call add      ("s1  ",1.d0,"s3  ",1.d0,"s1  ")
        call reorder  ("s1  ","s2  ",1324,4)
        call add      ("s2  ",1.d0,"bx24",1.d0,"bx24")
!  diag 4
        call mult     ("vhhp","tvm2","s1  ",222)
        call reorder  ("s1  ","s2  ",1342,4)
        call add      ("s2  ",1.d0,"bx24",1.d0,"bx24")
!  diag 8
        call mult     ("vhpp","t1v ","s1  ",311)
        call reorder  ("s1  ","s2  ",1432,4)
        call add      ("s2  ",-1.d0,"bx24",1.d0,"bx24")
!  diag 9
        call mult     ("t1v ","vphp","s1  ",131)
        call add      ("s1  ",-1.d0,"bx24",1.d0,"bx24")
!  diag 10
        call reorder  ("vvhh","s1  ",1243,4)
        call mult     ("s1  ","t1cr","s2  ",311)
        call reorder  ("s2  ","s1  ",1243,4)
        call add      ("s1  ",1.d0,"bx24",1.d0,"bx24")
      else if (iparm.eq.1) then
        call clear    ("bx24",3332,1100)
!  diag 5
        call mult     ("tm  ","vd14","s1  ",222)
        call add      ("s1  ",-1.d0,"bx24",1.d0,"bx24")
! v1 Feb 95
        call mult     ("tm  ","hp  ","s1  ",311)
        call add      ("bx24",1.d0,"s1  ",-.5d0,"bx24")
!  diag 6
        call t1t2     ("t1v ","tv  ","t32 ","s1  ",123456,2)
        call cmult    ("t1c ",0.5d0,"s2  ")
        call t1t2     ("s2  ","tu  ","s1  ","s3  ",123456,3)
        call reord6   ("s3  ","s1  ",124356)
        call add      ("hhpp",-2.d0,"va6 ",1.d0,"s2  ")
        call mult6    ("s1  ","s2  ","s3  ",2)
        call add      ("s3  ",1.d0,"bx24",1.d0,"bx24")
!  diag 7
        call t1t2     ("t1v ","tv  ","t32 ","s1  ",123456,1)
        call t1t2     ("t1v ","tv  ","s1  ","s2  ",123456,2)
        call t1t2     ("t1c ","tm  ","s2  ","s1  ",123456,3)
        call reord6   ("s1  ","s2  ",126354)
        call mult6    ("s2  ","hhpp","s1  ",2)
        call add      ("s1  ",1.d0,"bx24",1.d0,"bx24")
      end if
      return
      end 
!  *******************************************************************
      subroutine box25_csu (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      if (iparm.eq.0) then
!  diag 1
        call cmult    ("vvhh",0.5d0,"bx25")
!  diag 2
        call mult     ("t1v ","vhhp","s1  ",131)
        call add      ("s1  ",1.d0,"bx25",1.d0,"bx25")
      else if (iparm.eq.1) then
        call clear    ("bx25",3322,1100)
!  diag 3   
        call mult     ("tm  ","hhpp","s1  ",222)
        call add      ("s1  ",0.5d0,"bx25",1.d0,"bx25")
      end if
      return
      end
!  *******************************************************************
!     subroutine rdbox
      subroutine rdbox_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call rdfile  (20,0,"bx1 ",-1)
      call rdfile  (20,2,"bx3 ",-1)
      call rdfile  (20,3,"bx4 ",-1)
      call rdfile  (20,4,"bx5 ",-1)
      call rdfile  (20,6,"bx7 ",-1)
      call rdfile  (20,7,"bx8 ",-1)
      call rdfile  (21,0,"bx12",-1)
      call rdfile  (21,1,"bx13",-1)
      call rdfile  (21,2,"bx14",-1)
      call rdfile  (21,3,"bx15",-1)
      call rdfile  (21,4,"bx16",-1)
      call rdfile  (21,5,"bx17",-1)
      call rdfile  (21,6,"bx18",-1)
      call rdfile  (21,7,"bx19",-1)
      return
      end
!  *******************************************************************
!     subroutine wrtbox
      subroutine wrtbox_csu
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!
      call box24_csu(0)
      call wrtfile  (22,0,"bx24")
      call box24_csu(1)
      call reads    (22,0,"s1  ")
      call add      ("s1  ",1.d0,"bx24",1.d0,"bx24")
      call wrtfile  (22,0,"bx24")
!
      call box25_csu(0)
      call wrtfile  (22,1,"bx25")
      call box25_csu(1)
      call reads    (22,1,"s1  ")
      call add      ("s1  ",1.d0,"bx25",1.d0,"bx25")
      call wrtfile  (22,1,"bx25")
!
      return
      end
!  *******************************************************************
