!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine csv
!-------------------------------------------------------------------------------------
!
!  (0,1) sector
!
!-------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use gen_mod
      use cij1_mod
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
#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "cij1.h"
!      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cfiles.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)
#include "cint.h"
!     common /cint/ nint,idum1,vint(ndint),iint(ndint)
#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "csdtime.h"
!     common /csdtime/ csdstim,csdltim
#endif

      character*256 filenm,filenm51,filenm52,filenm6,filenm20,filenm21

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
      filenm20(1:1)=' '
      call getenv('tape20',filenm20)
      if (filenm20(1:1).eq.' ') filenm20='tape20'
      open(20,file=filenm20,form='UNFORMATTED')
      filenm21(1:1)=' '
      call getenv('tape21',filenm21)
      if (filenm21(1:1).eq.' ') filenm21='tape21'
      open(21,file=filenm21,form='UNFORMATTED')
      filenm(1:1)=' '
      call getenv('tape12',filenm)
      if (filenm(1:1).eq.' ') filenm='tape12'
      open(12,file=filenm,form='UNFORMATTED')
      write (*,*) '   $$$$$$$$$   csv (0,1)  $$$$$$$$$  '
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
      convrg = 10.d0**(-iopt(3))
      write (6,'('' max # of iter='',i3,''. convrg='',e10.2)')          &
  iopt(2),convrg
      call prtime("rdin    ",0,0)
      ! incomplete active space stuff
      if (iopt(9).ne.0) call rdrsp1
      nflcur=0
      nstrt(1,0)=0
      nstrt(2,0)=0
      call rdfile (51,0,"tc  ",-1)     
      call rdfile (51,1,"t1c ",-1)        
      if (iopt(17).ne.0) then
        call rdfile (51,2,"t30 ",-1)        
        if (iopt(17).eq.3) then
          call rdbox_csv
        end if
      end if
      call t1sqr("t1c ","t1c ","tc  ","tcm ")
      call ordvv_csv
      call calvv0_csv
      call pdiag_csv
      iter=0
      ieof=0
      call tvinit_csv
      call prt ("tv  ","tv      ",2)
      call fmax ("tv  ",tmax,itmax)
      call fmax ("t1v ",t1max,it1max)
      if (iopt(4).gt.0) then
        if (iopt(17).eq.0) then
          call dbdiag  ("vv10","vv0 ","dum3","dum4",1100,1)
        else if (iopt(17).ne.0) then
          call dbdiag  ("vv10","vv0 ","v30 ","dum4",1110,1)
        end if
        call wrtrle ("t1v ",1,1) 
        call wrtrle ("tv  ",2,1) 
        if (iopt(17).ne.0) then
          call wrtrle ("t31 ",3,1) 
        end if
      endif
      write(6,2400) tmax,idec(itmax),t1max,idec(it1max)
      call t1sqr("t1v ","t1c ","tv  ","tm  ")

220   iter=iter+1
      call prtime("strt itr",iter,0)
      call flush (6)
      call ordtva_csv
      if (iopt(17).ne.0) then
        call caltv3_csv
      end if
      call caltv2_csv
      call caltv1_csv
      call ordtvb_csv
      if (iopt(6).ne.0.and.iter.ge.mod(iopt(6),100)) then
        fnew=max(dfloat(iopt(6)/100),1.d0)
        if (fnew.le.6.5d0) then
          fnew=1.d0/fnew
        else
          fnew=fnew-5.d0
        endif
        call add("tnew",fnew/(1.d0+fnew),"tv  ",                        &
   1.d0/(1.d0+fnew),"tnew")
        call add("t1nw",fnew/(1.d0+fnew),"t1v ",                        &
   1.d0/(1.d0+fnew),"t1nw")
        if (iopt(17).ne.0) then
          call add("t3nw",fnew/(1.d0+fnew),"t31 ",                      &
     1.d0/(1.d0+fnew),"t3nw")
        end if
        write (6,'('' t averaged with fnew '',f7.3)') fnew
      endif
      call fmax("tnew",tmax,itmax)
      call fmax("t1nw",t1max,it1max)
      if (iopt(17).ne.0) then
        call fmax6 ("t3nw",t3max,id3mx,ind3mx)
      end if
      call diffmax("tv  ","tnew","s2  ",diff,idiff)
      call diffmax("t1v ","t1nw","s2  ",diff1,idiff1)
      df3 = 0.0d0 ! mi initialization
      if (iopt(17).ne.0) then
        call dfmax3 ("t31 ","t3nw","s2  ",df3,iddf3,inddf3)
      end if
      write (6,*) ' '
      write (6,4321) iter,diff,idec(idiff),diff1,idec(idiff1)
 4321 format (" iteration # ",i3,"  diffmx= ",e10.3,2x,i12,e10.3,2x,i12)
      if (iopt(17).ne.0) then
        write (6,2401) df3,idec(iddf3),idec(inddf3)
      end if
2401  format (24x,e10.3,2x,i12,2x,i12)
      write (6,2400) tmax,idec(itmax),t1max,idec(it1max)
2400  format(' max t, t1=',f15.9,2x,i12,f17.9,2x,i6)
      if (iopt(17).ne.0) then
        write (6,2402) t3max,idec(id3mx),idec(ind3mx)
      end if
2402  format(' maxt3 =',f15.9,2x,i12,2x,i12)
      if (abs(tmax).gt.1.d10.or.abs(tmax).lt.1.d-10) then
        write (*,*) ' abort. tmax out of bounds'
        call abext('tmax    ')    
      endif
! diis
      if (iopt(4).lt.0) then
        call diis(iter,icnvg,11,"t1nw","t1v ","tnew","tv  ")
        if(iter.gt.1) then
          call diffmax("tv  ","tnew","s2  ",diff,idiff)
          call diffmax("t1v ","t1nw","s2  ",diff1,idiff1)
          write (6,'('' iteration # '',i3,'' diis diffmx='',e10.3,2x,   &
      i12,e10.3,2x,i12)')                                               &
      iter,diff,idec(idiff),diff1,idec(idiff1)
        endif
        if(icnvg.eq.1) go to 910
      endif
      call move("tnew","tv  ")
      call move("t1nw","t1v ")
      if (iopt(17).ne.0) then
        call move ("t3nw","t31 ")    
      end if
      call t1sqr("t1v ","t1c ","tv  ","tm  ")
      call prt("tnew","new t   ",2)
      call prt("t1nw","new t1  ",2)
      call prt("tvf ","new tf  ",2)
! runtime check: The variable 'csv_$DF3' is being used without being defined
      if (max(abs(diff),abs(diff1),abs(df3)).lt.convrg) go to 910
      call wrtfile(52,0,"tv  ")
      call wrtfile(52,1,"t1v ")
      call wrtfile(52,2,"tvf ")
      if (iopt(17).ne.0) then
        call wrtfile(52,3,"t31 ")
        if (iopt(17).eq.3) then
          call wrtbox_csv
        end if
      end if
      if (iopt(4).gt.0) go to 260
!   no rle
      if (iter.ge.iopt(2)) go to 920
      go to 220
!   with rle
260   continue
      itr1=mod(iter-1,iopt(4))+1
      call wrtrle  ("t1v ",1,itr1+1)
      call wrtrle  ("tv  ",2,itr1+1)
      if (iopt(17).ne.0) then
        call wrtrle  ("t31 ",3,itr1+1)
      end if
      if (iter.ge.iopt(2)) go to 920
      if (itr1.lt.iopt(4)) go to 220
      if (iopt(17).eq.0) then
        call rle (itr1,1100,"t1v ","tv  ","dum3","dum4",ising)
      else
        call rle (itr1,1110,"t1v ","tv  ","t31 ","dum4",ising)
      end if
!  rle singular
      if (ising.eq.1) then
        if (iopt(17).eq.0) then
          call rle (itr1-1,1100,"t1v ","tv  ","dum3","dum4",ising)
        else
          call rle (itr1-1,1110,"t1v ","tv  ","t31 ","dum4",ising)
        end if
        if (ising.eq.1) then
          iopt(4)=0
          write (*,*) ' ***** rle singular, turned off'
          go to 220
        endif
      endif
      call t1sqr("t1v ","t1c ","tv  ","tm  ")
      call wrtrle ("t1v ",1,1)
      call wrtrle ("tv  ",2,1)
      if (iopt(17).ne.0) then
        call wrtrle ("t31 ",3,1)
      end if
      go to 220
!
!  different endings
910   write (6,'('' csv converged after iter '',i3)') iter
      go to 970
920   write (6,*) ' iteration limit'
      call abext('limit   ')       
970   call prt("tv  ","last t  ",0)
      call prt("t1v ","last t1 ",0)
      call prt("tvf ","last tf ",0)
      if (iopt(17).ne.0) then
        call prt6 ("t31 ","t31 ",0)
      end if
      call wrtfile(52,0,"tv  ")
      call wrtfile(52,1,"t1v ")
      call wrtfile(52,2,"tvf ")
      if (iopt(17).ne.0) then
        call wrtfile(52,3,"t31 ")
        if (iopt(17).eq.3) then
          call wrtbox_csv
        end if
      end if

      egr0=0.d0
      call eneru2 (egr0,1,3,max(1,iopt(1)))

      write (*,*) nstrt(1,nflcur),' words used in igen,',nstrt(2,nflcur),' words used in gen'

      call prtime("end",iter,0)
      call flush (6)

      end

! *********************************************
!     subroutine ordvv
      subroutine ordvv_csv
!---------------------------------------------------------------------------------------
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
!              hppp 2314 = vd8 
!              hppp 2413 = vd7 
!              hppp 2134 = vd12
!              hppp 2143 = vd14
!              vhpp 1324 = vva1
!              vhhp 1324 = vvc1
!
!---------------------------------------------------------------------------------------
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
      call seltyp   (2323,0,"hphp",4)
      call seltyp   (2223,0,"hhhp",4)
      call seltyp   (2333,0,"hppp",4)
      call seltyp4  (3333,0,"pppp",4)
      call seltyp   (3233,1000,"vhpp",4)
      call seltyp   (3323,1000,"vphp",4)
      call seltyp   (3223,1000,"vhhp",4)
      call seltyp   (3222,1000,"vhhh",4)
      call reorder  ("hhpp","va1 ",1324,4)
      call reorder  ("hhpp","va2 ",1423,4)
      call reorder  ("hhpp","va4 ",3412,4)
      call reorder  ("hhpp","va6 ",1243,4)
      call reorder  ("hhpp","va8 ",3421,4)
      call reorder  ("hhhp","vc4 ",3412,4)
      call reorder  ("hhhp","vc8 ",3421,4)
      call reorder  ("hhhp","vc9 ",3124,4)
      call reorder  ("hhhp","vc10",3214,4)
      call reorder  ("hhhp","vc11",2413,4)
      call reorder  ("hppp","vd8 ",2314,4)
      call reorder  ("hppp","vd7 ",2413,4)
      call reorder  ("hppp","vd12",2134,4)
      call reorder  ("hppp","vd14",2143,4)
      call reorder  ("vhpp","vva1",1324,4)
      call reorder  ("vhhp","vvc1",1324,4)
      call selt1    (33,0,"pp  ",4)
      call selt1    (23,0,"hp  ",4)
      call selt1    (22,0,"hh  ",4)
      call selt1    (33,10,"vp  ",4)
      call selt1    (32,10,"vh  ",4)
      call selt1    (32,0,"ph  ",4)
      if (iopt(17).ne.0) then
        call seltyp (3232,1000,"vhph",4)
        call seltyp (3333,1000,"vppp",4)
      end if
      call prtime("ordv    ",0,0)
      return
      end
!**************************************************************
      subroutine calvv0_csv
!-------------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------------
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

      call ordta0_csv
      call move("vhpp","vv0 ")
      call diag3a_csv
      call diag4b_csv
      call diag5b_csv
      call diag6b_csv
      call diag8_csv
      call diag20b_csv
      call diag18b_csv
      call diag22a_csv
      call diag26a_csv
      if (iopt(17).ge.1) then
!  diag 31a
        call reord6  ("t30 ","s1  ",124356)
        call reord7  ("s1  ","s2  ",132456)
        call reord6  ("s2  ","s1  ",126345)
        call add     ("vhhp",-2.d0,"vvc1",1.d0,"s2  ")
        call mult6   ("s2  ","s1  ","s3  ",2)
!  diag 33a
        call reord6  ("t30 ","s1  ",125346)
        call reord7  ("s1  ","s2  ",132456)
        call reord6  ("s2  ","s1  ",126345)
        call mult6   ("vhhp","s1  ","s2  ",2)
        call add     ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
        call reorder ("s3  ","s1  ",1243,4)
        call add     ("s1  ",1.d0,"vv0 ",1.d0,"vv0 ")
      end if
      if (iopt(17).ge.2) then
!  diag 36a,38a
        call reord6    ("t30 ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125634)
        call add       ("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","boxA",2)
        call reord6    ("t30 ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126543)
        call mult6     ("s1  ","va4 ","s2  ",2)
        call add       ("s2  ",1.d0,"boxA",1.d0,"boxA")
      end if
      call ordtb0_csv
      call move("vp  ","vv10")
      call dg1n6_csv
      call dg1n7_csv
      call dg1n8_csv
      if ((iopt(17).eq.1).or.(iopt(17).eq.2)) then
!  dg3n1c 
        call mult6    ("vppp","tc  ","s1  ",1)
        call reord6   ("s1  ","s2  ",124356)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","v30 ",125463)
!  dg3n2b
        call reorder  ("vhhp","s1  ",1243,4)
        call mult6    ("s1  ","tc3 ","s2  ",1)
        call reord6   ("s2  ","s1  ",124536)
        call add      ("s1  ",-1.d0,"v30 ",1.d0,"v30 ")
!  dg3n2c
        call mult6    ("vhph","tc3 ","s1  ",1)
        call reord6   ("s1  ","s2  ",124356)
        call add      ("s2  ",-1.d0,"v30 ",1.d0,"v30 ")
! 
        call reord6   ("v30 ","s1  ",123465)
        call reord7   ("s1  ","s2  ",132456)
        call add      ("v30 ",1.d0,"s2  ",1.d0,"v30 ")
      else if (iopt(17).eq.3) then
!  GROUP I
!   diag 3
!       call box12    (0)
        call box12_csv(0)
        call mult6    ("bx12","tc  ","s1  ",1)
        call reord6   ("s1  ","s2  ",124356)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","v30 ",125364)
!  GROUP II
!   diag 2
        call box13_csv(0)
        call mult6    ("bx13","tc3 ","s1  ",1)
        call reord6   ("s1  ","s2  ",124536)
        call add      ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!   diag 3
        call box14_csv(0)
        call mult6    ("bx14","tc3 ","s1  ",1)
        call reord6   ("s1  ","s2  ",124356)
        call add      ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP III
!   diag 2
        call box15_csv(0)
        call reord6   ("t30 ","s1  ",124536)
        call mult6    ("bx15","s1  ","s2  ",3)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","s2  ",124356)
        call add      ("s2  ",0.5d0,"v30 ",1.d0,"v30 ")
!  GROUP IV
!   diag 3
        call box16_csv(0)
        call reord6   ("t30 ","s1  ",124635)
        call mult6    ("bx16","s1  ","s2  ",3)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","s2  ",124356)
        call add      ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP V
!   diag 3
        call box17_csv(0)
        call reord6   ("t30 ","s1  ",124635)
        call mult6    ("bx17","s1  ","s2  ",3)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","s2  ",124653)
        call add      ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP VI
!   diag 2
        call box18_csv(0)
        call reord6   ("t30 ","s1  ",126345)
        call reord7   ("s1  ","s2  ",132456)
        call reord6   ("s2  ","s1  ",126543)
        call mult6    ("bx18","s1  ","s2  ",3)
        call add      ("s2  ",1.d0,"v30 ",1.d0,"v30 ")
!  GROUP IX
!   diag 2
        call box19_csv(0)
        call reord6   ("t30 ","s1  ",124563)
        call mult6    ("s1  ","bx19","s2  ",5)
        call reord6   ("s2  ","s1  ",126534)
        call reord7   ("s1  ","s2  ",312456)
        call add      ("s2  ",0.5d0,"v30 ",1.d0,"v30 ")
!
        call reord6   ("v30 ","s1  ",123465)
        call reord7   ("s1  ","s2  ",132456)
        call add      ("v30 ",1.d0,"s2  ",1.d0,"v30 ")
      end if
      call prtime     ("calvv0  ",0,0)
      return
      end
! ****************************************
!     subroutine pdiag
      subroutine pdiag_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call pdig9_csv
      call pdig11_csv
      call pdig12_csv
      call pdig13a_csv
      call pdig14_csv
      call pdig23b_csv
      call pdg1n10
      call pdg1n11
      call pdgt1p4
      return
      end
! ******************************************
!     subroutine tvinit
      subroutine tvinit_csv
!-------------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------------
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
          call rdfile(52,0,"tv  ",-1)
          call rdfile(52,1,"t1v ",-1)
          call rdfile(52,2,"tvf ",-1)
          if ((iopt(17).ne.0).and.(iopt(20).eq.1))                      &
        call rdfile(52,3,"t31 ",-1)
       else
        call move    ("vv0 ","tv  ")
        call diveps  ("tv  ")
        call move    ("vv10","t1v ")
        call fold_ccsd ("t1v ","tvf ")
        call diveps  ("t1v ")
       end if
          if ((iopt(17).ne.0).and.(iopt(20).ne.1)) then
           call move    ("v30 ","t31 ")
           call dveps6  ("t31 ")
          end if
      return
      end
! **************************************************
      subroutine ordta0_csv
!-------------------------------------------------------------------------------------
!  tc1 = tc 1324
!  tcm1 = tcm 1324
!  tc2 = tc 1423
!  tcm2 = tcm 1423
!  tc3 = tc 1432
!  tc4 = tc 3412
!  tcm4 = tcm 3412
!  t1cor = t1c 21
!  tb1 = tnew 1324
!  tb4 = tnew 3412
!  tb15 = tnew 1342
!-------------------------------------------------------------------------------------
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call reorder("tc  ","tc1 ",1324,4)
      call reorder("tc  ","tc2 ",1423,4)
      call reorder("tc  ","tc3 ",1432,4)
      call reorder("tc  ","tc4 ",3412,4)
      call reorder("tcm ","tcm1",1324,4)
      call reorder("tcm ","tcm2",1423,4)
      call reorder("tcm ","tcm4",3412,4)
      call reorder("t1c ","t1cr",1243,4)
      call clear("tb1 ",3323,1000)
      call clear("tb4 ",3332,10)
      call clear("tb15",3323,1000)
      call clear("tvf ",33,11)
      call prtime("ordta0  ",0,1)
      return
      end
! ******************************************************
      subroutine ordtb0_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call reorder("tb1 ","s2  ",1324,4)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call reorder("tb4 ","s2  ",3412,4)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call reorder("tb15","s2  ",1342,4)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call prtime("ordtb0  ",0,1)
      return
      end
! *****************************************************
!     subroutine ordtva
      subroutine ordtva_csv

!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call reorder("tv  ","ta1 ",1324,4)
      call reorder("tv  ","ta2 ",1423,4)
      call reorder("tm  ","tm2 ",1423,4)
      call reorder("tv  ","ta3 ",1432,4)
      call reorder("tv  ","ta6 ",1243,4)
      call reorder("tv  ","ta15",1342,4)
      call reorder("t1v ","t1or",1243,4)
      call clear("tb1 ",3323,1000)
      call clear("tb2 ",3332,1000)
      call clear("tb4 ",3332,10)
      call clear("tb6 ",3233,1000)
      call clear("tb15",3323,1000)
      call clear("tvf ",33,11)
      call prtime("ordta   ",0,1)
      return
      end
! ****************************************************
!     subroutine caltv3
      subroutine caltv3_csv
!-------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------
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
!  dg3n1a
        call mult6     ("tv  ","hppp","s1  ",1)
        call reord6    ("s1  ","t3nw",124356)
!  dg3n1b
        call mult6     ("ta6 ","hppp","s1  ",1)
        call reord6    ("s1  ","s2  ",124536)
        call add       ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  dg3n2a
        call mult6     ("ta3 ","vc11","s1  ",1)
        call reord6    ("s1  ","s2  ",124356)
        call reord7    ("s2  ","s1  ",132456)
        call reord6    ("s1  ","s2  ",126453)
        call add       ("s2  ",-1.d0,"t3nw",1.d0,"t3nw")
      else if (iopt(17).eq.3) then
!  GROUP I
!   diag 1
        call mult6    ("tv  ","bx1 ","s1  ",1)
        call reord6   ("s1  ","t3nw",124365)
!   diag 2
        call mult6    ("ta6 ","bx1 ","s1  ",1)
        call reord6   ("s1  ","s2  ",124635)
        call add      ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 3
!       call box12    (1)
        call box12_csv(1)
        call mult6    ("bx12","tc  ","s1  ",1)
        call reord6   ("s1  ","s2  ",124356)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","s2  ",125364)
        call add      ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP II
!   diag 1
        call mult6    ("ta3 ","bx2 ","s1  ",1)
        call reord6   ("s1  ","s2  ",124356)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","s2  ",125463)
        call add      ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 2
        call box13_csv(1)
        call mult6    ("bx13","tc3 ","s1  ",1)
        call reord6   ("s1  ","s2  ",124536)
        call add      ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 3
        call box14_csv(1)
        call mult6    ("bx14","tc3 ","s1  ",1)
        call reord6   ("s1  ","s2  ",124356)
        call add      ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP III
!   diag 1
        call reord6   ("t31 ","s1  ",124536)
        call mult6    ("s1  ","bx3 ","s2  ",3)
        call reord6   ("s2  ","s1  ",125346)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 2
        call box15_csv(1)
        call reord6   ("t30 ","s1  ",124536)
        call mult6    ("bx15","s1  ","s2  ",3)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","s2  ",124356)
        call add      ("s2  ",0.5d0,"t3nw",1.d0,"t3nw")
!  GROUP IV
!   diag 1
        call reord6   ("t31 ","s1  ",124635)
        call mult6    ("s1  ","bx4 ","s2  ",3)
        call reord6   ("s2  ","s1  ",125346)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 2
        call reord6   ("t31 ","s1  ",125634)
        call mult6    ("s1  ","bx4 ","s2  ",3)
        call reord6   ("s2  ","s1  ",125436)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 3
        call box16_csv(1)
        call reord6   ("t30 ","s1  ",124635)
        call mult6    ("bx16","s1  ","s2  ",3)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","s2  ",124356)
        call add      ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP V
!   diag 1
        call reord6   ("t31 ","s1  ",124635)
        call mult6    ("s1  ","bx5 ","s2  ",3)
        call reord6   ("s2  ","s1  ",125364)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 2
        call reord6   ("t31 ","s1  ",125634)
        call mult6    ("s1  ","bx5 ","s2  ",3)
        call reord6   ("s2  ","s1  ",125634)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 3
        call box17_csv(1)
        call reord6   ("t30 ","s1  ",124635)
        call mult6    ("bx17","s1  ","s2  ",3)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","s2  ",124653)
        call add      ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP VI
!   diag 1
        call reord6   ("t31 ","s1  ",124356)
        call reord7   ("s1  ","s2  ",132456)
        call reord6   ("s2  ","s1  ",125634)
        call mult6    ("s1  ","bx6 ","s2  ",3)
        call reord6   ("s2  ","s1  ",125346)
        call reord7   ("s1  ","s2  ",132456)
        call reord6   ("s2  ","s1  ",126345)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 2
        call box18_csv(1)
        call reord6   ("t30 ","s1  ",126345)
        call reord7   ("s1  ","s2  ",132456)
        call reord6   ("s2  ","s1  ",126543)
        call mult6    ("bx18","s1  ","s2  ",3)
        call add      ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP VII
!   diag 1a
        call t1t2     ("t1c ","tv  ","t31 ","s1  ",123456,2)
        call t1t2     ("t1c ","tv  ","s1  ","s2  ",123456,3)
        call mult6    ("s2  ","pppp","s1  ",3)
        call reord6   ("s1  ","s2  ",123465)
        call reord7   ("s2  ","s3  ",132456)
        call add      ("s3  ",1.d0,"s1  ",1.d0,"s1  ")
        call add      ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
!   diag 1b
        call t1t2     ("t1v ","tc  ","t31 ","s1  ",123456,1)
        call t1t2     ("t1c ","tv  ","s1  ","s3  ",123456,2)
        call reord6   ("s3  ","s1  ",123645)
        call mult6    ("s1  ","pppp","s2  ",3)
        call reord6   ("s2  ","s1  ",123564)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
        call reord6   ("s3  ","s1  ",123654)
        call reord7   ("s1  ","s2  ",213456)
        call mult6    ("s2  ","pppp","s1  ",3)
        call reord6   ("s1  ","s2  ",123654)
        call reord7   ("s2  ","s1  ",213456)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 2a
        call mult6    ("t31 ","vd14","s1  ",3)
        call mult6    ("s1  ","t1cr","s2  ",5)
        call add      ("s2  ",-1.d0,"t3nw",1.d0,"t3nw")
!   diag 2b
        call reord6   ("t31 ","s1  ",123645)
        call mult6    ("s1  ","vd12","s2  ",3)
        call mult6    ("s2  ","t1cr","s3  ",5)
        call reord6   ("s3  ","s2  ",123654)
        call add      ("s2  ",-1.d0,"t3nw",1.d0,"t3nw")
!   diag 2c
        call mult6    ("s1  ","vd14","s2  ",3)
        call mult6    ("s2  ","t1cr","s1  ",5)
        call reord6   ("s1  ","s2  ",123564)
        call add      ("s2  ",-1.d0,"t3nw",1.d0,"t3nw")
!   diag 3a
        call mult6    ("t31 ","hhpp","s1  ",3)
        call mult6    ("s1  ","tcm4","s2  ",3)
        call add      ("s2  ",0.5d0,"t3nw",1.d0,"t3nw")
!   diag 3b
        call reord6   ("t31 ","s1  ",123645)
        call mult6    ("s1  ","hhpp","s2  ",3)
        call mult6    ("s2  ","tcm4","s1  ",3)
        call reord6   ("s1  ","s2  ",123564)
        call add      ("s2  ",1.d0,"t3nw",1.d0,"t3nw")
!  GROUP VIII
!   diag 1
        call mult6    ("t31 ","bx7 ","s1  ",5)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 2
        call reord6   ("t31 ","s1  ",123564)
        call mult6    ("s1  ","bx7 ","s2  ",5)
        call reord6   ("s2  ","s1  ",123645)
        call add      ("s1  ",0.5d0,"t3nw",1.d0,"t3nw")
!  GROUP IX
!   diag 1
        call reord6   ("t31 ","s1  ",124563)
        call mult6    ("s1  ","bx8 ","s2  ",5)
        call reord6   ("s2  ","s1  ",126345)
        call add      ("s1  ",1.d0,"t3nw",1.d0,"t3nw")
!   diag 2
        call box19_csv(1)
        call reord6   ("t30 ","s1  ",124563)
        call mult6    ("s1  ","bx19","s2  ",5)
        call reord6   ("s2  ","s1  ",126534)
        call reord7   ("s1  ","s2  ",312456)
        call add      ("s2  ",0.5d0,"t3nw",1.d0,"t3nw")
      end if
      call reord6     ("t3nw","s1  ",123465)
      call reord7     ("s1  ","s2  ",132456)
      call add        ("t3nw",1.d0,"s2  ",1.d0,"t3nw")
      call add        ("v30 ",1.d0,"t3nw",1.d0,"t3nw")
      call dveps6     ("t3nw")    
      call prtime     ("caltv3  ",0,1)
      return
      end
!  *******************************************************************
      subroutine caltv2_csv
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

      call move("vv0 ","tnew")
      call diag2a_csv
      call diag4a_csv
      call diag5a_csv
      call diag6a_csv
      call diag7_csv
      call diag9_csv
      call diag11_csv
      call diag12_csv
      call diag13a_csv
      call diag13b_csv
      call diag14_csv
      call diag15_csv
      call diag21a_csv
      call diag22b_csv
      call diag23b_csv
      call diag28a_csv
      call dgt1p4_csv
      if (iopt(17).ge.1) then
!  diag 30a
        call reord6    ("t31 ","s1  ",125364)
        call add       ("vd12",2.d0,"vd14",-1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",2)   
        call add       ("s3  ",1.d0,"tb6 ",1.d0,"tb6 ")
!  diag 30b
        call reord6    ("t31 ","s1  ",124365)
        call mult6     ("s1  ","s2  ","s3  ",2)
        call add       ("s3  ",1.d0,"tnew",1.d0,"tnew")
!  diag 31b
        call reord6    ("t31 ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",125436)
        call reorder   ("hhhp","s2  ",1324,4)
        call add       ("hhhp",-2.d0,"s2  ",1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",2)
        call add       ("s3  ",1.d0,"tb2 ",1.d0,"tb2 ")
!  diag 32a
        call reord6    ("t31 ","s1  ",126354)
        call mult6     ("s1  ","vd12","s2  ",2)
        call add       ("s2  ",-1.d0,"tb6 ",1.d0,"tb6 ")
!  diag 32b
        call reord6    ("t31 ","s1  ",126345)
        call mult6     ("s1  ","vd12","s2  ",2)
        call add       ("s2  ",-1.d0,"tnew",1.d0,"tnew")
!  diag 33b
        call reord6    ("t31 ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126345)
        call mult6     ("s1  ","hhhp","s2  ",2)
        call add       ("s2  ",1.d0,"tb2 ",1.d0,"tb2 ")
      end if
      if (iopt(17).ge.2) then
!    diag 34,35a,35b
        call reord6    ("t31 ","s3  ",124536)
        call reord6    ("t31 ","s1  ",124635)
        call add       ("s1  ",-1.d0,"s3  ",2.d0,"s3  ")
        call reord6    ("t31 ","s1  ",126534)
        call add       ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
        call add       ("va1 ",2.d0,"va2 ",-1.d0,"s1  ")
        call mult      ("s1  ","t1c ","s2  ",202)
! v1 Feb 95
        call add       ("s2  ",1.d0,"hp  ",1.d0,"s2  ")
        call mult6     ("s3  ","s2  ","s1  ",0)
        call add       ("s1  ",1.d0,"tnew",1.d0,"tnew")
!    diag 36a,38a
        call mult      ("t1v ","boxA","s1  ",131)
        call add       ("s1  ",1.d0,"tb6 ",1.d0,"tb6 ")
!    diag 36b,38b
        call reord6    ("t31 ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126543)
        call add       ("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",2)
        call reord6    ("t31 ","s1  ",125346)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("s2  ","s1  ",126543)
        call mult6     ("s1  ","va4 ","s2  ",2)
        call add       ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
        call mult      ("s3  ","t1c ","s1  ",311)
        call add       ("s1  ",1.d0,"tb2 ",1.d0,"tb2 ")
!    diag 37a,39a
        call reord6    ("t31 ","s1  ",125346)
        call add       ("hhpp",-2.d0,"va6 ",1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",2)
        call reord6    ("t31 ","s1  ",126345)
        call mult6     ("s1  ","hhpp","s2  ",2)
        call add       ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
        call mult      ("s3  ","t1cr","s1  ",311)
        call add       ("s1  ",1.d0,"tb6 ",1.d0,"tb6 ")
!    diag 37b,39b
        call reord6    ("t31 ","s1  ",124356)
        call add       ("hhpp",-2.d0,"va6 ",1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",2)
        call reord6    ("t31 ","s1  ",126354)
        call mult6     ("s1  ","hhpp","s2  ",2)
        call add       ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
        call mult      ("s3  ","t1cr","s1  ",311)
        call add       ("s1  ",1.d0,"tnew",1.d0,"tnew")
      end if
      call prtime     ("caltv2  ",0,1)
      return
      end
! **********************************************************
      subroutine caltv1_csv
!-------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------
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

      call move("vv10","t1nw")
      call dg1n2_csv
      call dg1n3_csv
      call dg1n5_csv
      call dg1n9_csv
      call dg1n10_csv
      call dg1n11_csv
      if (iopt(17).ne.0) then
!  dg1n13,14
        call reord6    ("t31 ","s1  ",124356)
        call reord7    ("s1  ","s2  ",132456)
        call reord6    ("t31 ","s1  ",125346)
        call reord7    ("s1  ","s3  ",132456)
        call add       ("s2  ",1.d0,"s3  ",-1.d0,"s1  ")
        call add       ("hhpp",2.d0,"va6 ",-1.d0,"s2  ")
        call mult6     ("s1  ","s2  ","s3  ",0)
        call add       ("s3  ",1.d0,"t1nw",1.d0,"t1nw")
      end if
!
      call fold_ccsd("t1nw","tvf ")
      call heff1_csv
      call dfold1_csv
      call diveps("t1nw")
      call prtime     ("caltv1  ",0,1)
      return
      end
!  ********************************************************************
      subroutine ordtvb_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
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
      call diveps("tnew")
      call prtime("ordtvb  ",0,1)
      return
      end
! *************************************************
      subroutine diag7_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   diag7 in two parts, as only half of pppp collected
      call mult("pppp","tm  ","s2  ",222)
      call add("tb4 ",1.d0,"s2  ",1.d0,"tb4 ")
      call reorder("tm  ","s2  ",2143,4)
      call mult("pppp","s2  ","s3  ",222)
      call reorder("s3  ","s2  ",2143,4)
      call add("tb4 ",1.d0,"s2  ",1.d0,"tb4 ")
      call prt("s2  ","diag7   ",3)
      call prtime("diag7   ",0,1)
      return
      end
! **********************************************
      subroutine diag3a_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vh  ","tc3 ","s2  ",131)
      call add("vv0 ",1.d0,"s2  ",-1.d0,"vv0 ")
      call prt("s2  ","-diag3a ",3)
      call prtime("diag3a  ",0,1)
      return
      end
! ***************************************************
      subroutine diag2a_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

!   incl 2a, 2b, 3b
!   2a
      call mult("ta6 ","pp  ","s2  ",311)
      call add("tb6 ",1.d0,"s2  ",1.d0,"tb6 ")
      call prt("s2  ","diag2a  ",3)
!   2b
      call mult("tv  ","pp  ","s2  ",311)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","diag2b  ",3)
!   3b
      call mult("ta15","hh  ","s2  ",311)
      call add("tb2 ",1.d0,"s2  ",-1.d0,"tb2 ")
      call prt("s2  ","diag3b  ",3)
      call prtime("diag2a  ",0,1)
      return
      end
! ******************************************************
      subroutine diag4a_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call add("va1 ",2.d0,"hphp",-1.d0,"s2  ")
      call mult("ta1 ","s2  ","s3  ",222)
      call add("tb1 ",1.d0,"s3  ",1.d0,"tb1 ")
      call prt("s3  ","diag4a  ",3)
      call prtime("diag4a  ",0,1)
      return
      end
! *****************************************************
      subroutine diag4b_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call add("vva1",2.d0,"vphp",-1.d0,"s2  ")
      call mult("s2  ","tc1 ","s3  ",222)
      call add("tb1 ",1.d0,"s3  ",1.d0,"tb1 ")
      call prt("s3  ","diag4b  ",3)
      call prtime("diag4b  ",0,1)
      return
      end
! ********************************************
      subroutine diag5a_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call mult("tm2 ","va1 ","s3  ",222)
      call add("tb1 ",1.d0,"s3  ",-1.d0,"tb1 ")
      call prt("s3  ","-diag5a ",3)
      call prtime("diag5a  ",0,1)
      return
      end
! *********************************************************
      subroutine diag5b_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vva1","tcm2","s3  ",222)
      call add("tb1 ",1.d0,"s3  ",-1.d0,"tb1 ")
      call prt("s3  ","-diag5b ",3)
      call prtime("diag5b  ",0,1)
      return
      end
! *******************************************************
      subroutine diag6a_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tm2 ","hphp","s3  ",222)
      call add("tb15",1.d0,"s3  ",-1.d0,"tb15")
      call prt("s3  ","-diag6a ",3)
      call prtime("diag6a  ",0,1)
      return
      end
! **********************************************************
      subroutine diag6b_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vphp","tcm2","s3  ",222)
      call add("tb15",1.d0,"s3  ",-1.d0,"tb15")
      call prt("s3  ","-diag6b ",3)
      call prtime("diag6b  ",0,1)
      return
      end
! *************************************************************
      subroutine diag8_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vhhh","tcm4","s2  ",222)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call prt("s2  ","diag8   ",3)
      call prtime("diag8   ",0,1)
      return
      end
! ******************************************************
      subroutine pdig9_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 9 10a 20a 23a
      call add("tc1 ",-2.d0,"tcm2",1.d0,"s2  ")
      call add("va1 ",-2.d0,"va2 ",1.d0,"s3  ")
      call mult("s2  ","s3  ","pd9 ",222)
!   20a
      call add("vc8 ",-2.d0,"vc4 ",1.d0,"s2  ")
      call mult("t1cr","s2  ","s3  ",131)
      call reorder("s3  ","s2  ",2143,4)
      call add("pd9 ",1.d0,"s2  ",1.d0,"pd9 ")
!   23a
      call add("vd12",2.d0,"vd14",-1.d0,"s2  ")
      call mult("t1c ","s2  ","s3  ",131)
      call add("pd9 ",1.d0,"s3  ",1.d0,"pd9 ")
      call prt("pd9 ","pd9     ",3)
      call prtime("pdig9   ",0,1)
      return
      end
! **********************************************************
      subroutine diag9_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("ta1 ","pd9 ","s2  ",222)
      call add("tb1 ",1.d0,"s2  ",1.d0,"tb1 ")
      call prt("s2  ","diag9   ",3)
      call prtime("diag9   ",0,1)
      return
      end
! ************************************************************
      subroutine pdig11_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

!   incl diags 10b, 11,24a
!   10b
      call add("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult("tc1 ","s2  ","pd11",222)
!   11
      call mult("tcm2","va1 ","s2  ",222)
      call add("pd11",1.d0,"s2  ",1.d0,"pd11")
!   24a
      call mult("t1c ","vd12","s2  ",131)
      call add("pd11",1.d0,"s2  ",-1.d0,"pd11")
      call prt("pd11","pdig11  ",3)
      call prtime("pdig11  ",0,1)
      return
      end
! ************************************************
      subroutine diag11_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tm2 ","pd11","s2  ",222)
      call add("tb1 ",1.d0,"s2  ",1.d0,"tb1 ")
      call prt("s2  ","diag11  ",3)
      call prtime("diag11  ",0,1)
      return
      end
! ***************************************************
      subroutine pdig12_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl diags 12, 25a
      call mult("tcm2","va2 ","pd12",222)
!   25a
      call mult("t1c ","vd14","s2  ",131)
      call add("pd12",1.d0,"s2  ",-1.d0,"pd12")
      call prt("pd12","pdig12  ",3)
      call prtime("pdig12  ",0,1)
      return
      end
! **************************************************
      subroutine diag12_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tm2 ","pd12","s2  ",222)
      call add("tb15",1.d0,"s2  ",1.d0,"tb15")
      call prt("s2  ","diag12  ",3)
      call prtime("diag12  ",0,1)
      return
      end
! ************************************************************
      subroutine pdig13a_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl 13a 18a
      call add("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult("tcm1","s2  ","p13a",113)
!   18a
      call add("vc9 ",-2.d0,"vc10",1.d0,"s2  ")
      call mult("s2  ","t1c ","s3  ",202)
      call add("p13a",1.d0,"s3  ",1.d0,"p13a")
      call prt("p13a","pdig13a ",3)
      call prtime("pdig13a ",0,1)
      return
      end
! **********************************************************
      subroutine diag13a_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("ta15","p13a","s2  ",311)
      call add("tb2 ",1.d0,"s2  ",1.d0,"tb2 ")
      call prt("s2  ","diag13a ",3)
      call prtime("diag13a ",0,1)
      return
      end
! ******************************************************
      subroutine diag13b_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add("hhpp",-2.d0,"va6 ",1.d0,"s2  ")
      call mult("tm  ","s2  ","s3  ",113)
      call mult("s3  ","tc3 ","s2  ",131)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","diag13b ",3)
      call prtime("diag13b ",0,1)
      return
      end
! *********************************************************
      subroutine pdig14_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl 14ab, 19ab
      call add("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
      call mult("tcm4","s2  ","pd14",113)
!   19
      call add("vd7 ",2.d0,"vd8 ",-1.d0,"s2  ")
      call mult("s2  ","t1c ","s3  ",202)
      call add("pd14",1.d0,"s3  ",1.d0,"pd14")
      call prt("pd14","pdig14b ",3)
      call prtime("pdig14  ",0,1)
      return
      end
! *********************************************
      subroutine diag14_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   14a, 19a
      call mult("tv  ","pd14","s2  ",311)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","diag14a ",3)
!   14b, 19b
      call mult("ta6 ","pd14","s2  ",311)
      call add("tb6 ",1.d0,"s2  ",1.d0,"tb6 ")
      call prt("s2  ","diag14  ",3)
      call prtime("diag14  ",0,1)
      return
      end
! ******************************************************
      subroutine diag15_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("tm  ","hhpp","s3  ",222)
      call mult("s3  ","tcm4","s2  ",222)
      call add("tnew",1.d0,"s2  ",1.d0,"tnew")
      call prt("s2  ","diag15  ",3)
      call prtime("diag15  ",0,1)
      return
      end
! ************************************************************
      subroutine diag18b_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add("vvc1",-2.d0,"vhhp",1.d0,"s2  ")
      call mult("s2  ","t1c ","s3  ",202)
      call mult("s3  ","tc3 ","s2  ",131)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call prt("s2  ","diag18b ",3)
      call prtime("diag18b ",0,1)
      return
      end
! **********************************************************
      subroutine diag20b_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl 20b, 21b, 17a
      call add("vvc1",-2.d0,"vhhp",1.d0,"s2  ")
      call mult("tc1 ","s2  ","s3  ",222)
!   diag 21b
      call mult("tc2 ","vvc1","s2  ",222)
      call add("s3  ",1.d0,"s2  ",1.d0,"s3  ")
!   diag 17a
      call reorder("vhhp","s2  ",2413,5)
      call add("s3  ",1.d0,"s2  ",-1.d0,"s3  ")
! ---------------------------------
      call mult("s3  ","t1cr","s2  ",311)
      call reorder("s2  ","s3  ",3142,5)
      call add("vv0 ",1.d0,"s3  ",1.d0,"vv0 ")
      call prt("s3  ","diag20b ",3)
      call prtime("diag20b ",0,1)
      return
      end
! ***************************************************
      subroutine diag21a_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl 21a, 27a, 29a
!   27a
      call mult("tv  ","vd14","s3  ",222)
!   29a
      call mult("tv  ","hp  ","s2  ",311)
      call add("s3  ",-1.d0,"s2  ",-1.d0,"s3  ")
      call reorder("s3  ","s2  ",1324,5)
!   21a
      call mult("ta2 ","vc9 ","s3  ",222)
      call add("s2  ",1.d0,"s3  ",1.d0,"s2  ")
! ------------
      call mult("s2  ","t1cr","s3  ",311)
      call add("tb1 ",1.d0,"s3  ",1.d0,"tb1 ")
      call prt("s2  ","diag21a ",3)
      call prtime("diag21a ",0,1)
      return
      end
! ********************************************
      subroutine diag22a_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   diag 17b
      call seltyp(3232,1000,"s2  ",5)
      call reorder("s2  ","s3  ",2314,5)
!  incl 22a, 17b
      call mult("tc2 ","vhhp","s2  ",222)
      call add("s2  ",1.d0,"s3  ",-1.d0,"s2  ")
! ----------------------------------
      call mult("s2  ","t1cr","s3  ",311)
      call reorder("s3  ","s2  ",3124,5)
      call add("vv0 ",1.d0,"s2  ",1.d0,"vv0 ")
      call prt("s3  ","diag22a ",3)
      call prtime("diag22a ",0,1)
      return
      end
! **************************************************
      subroutine diag22b_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl 22b, 27b, 29b
      call mult("ta2 ","vc10","s2  ",222)
!   27b
      call mult("tv  ","hppp","s1  ",222)
      call reorder("s1  ","s3  ",1423,5)
      call add("s2  ",1.d0,"s3  ",-1.d0,"s2  ")
!   29b
      call mult("ta6 ","hp  ","s1  ",311)
      call reorder("s1  ","s3  ",1324,5)
      call add("s2  ",1.d0,"s3  ",-1.d0,"s2  ")
! -------------
      call mult("s2  ","t1cr","s1  ",311)
      call add("tb15",1.d0,"s1  ",1.d0,"tb15")
      call prt("s1  ","diag22b ",3)
      call prtime("diag22b ",0,1)
      return
      end
! ************************************************
      subroutine pdig23b_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl 23b, 24b, 25b, 26b, 16a, 18b
!   23b
      call add("vd7 ",2.d0,"vd8 ",-1.d0,"s3  ")
      call mult("tc1 ","s3  ","p23b",222)
!   24b
      call mult("tc2 ","vd7 ","s3  ",222)
      call add("p23b",1.d0,"s3  ",-1.d0,"p23b")
!   25b
      call mult("tc2 ","vd8 ","s3  ",222)
!   16a
        call add("s3  ",-1.d0,"hppp",1.d0,"s3  ")
!   25b & 16a
      call reorder("s3  ","s2  ",1324,5)
      call add("p23b",1.d0,"s2  ",1.d0,"p23b")
!   26b
      call mult("vc4 ","tcm4","s3  ",222)
      call reorder("s3  ","s2  ",1342,5)
      call add("p23b",1.d0,"s2  ",1.d0,"p23b")
!   28b
      call reorder("tc  ","s2  ",1342,5)
      call reorder("hp  ","s3  ",1243,5)
      call mult("s2  ","s3  ","s1  ",311)
      call add("p23b",1.d0,"s1  ",-1.d0,"p23b")
      call prt("p23b","p23b    ",3)
      call prtime("p23b    ",0,1)
      return
      end
! *************************************************
      subroutine diag23b_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("t1v ","p23b","s2  ",131)
      call add("tb6 ",1.d0,"s2  ",1.d0,"tb6 ")
      call prt("s2  ","diag23b ",3)
      call prtime("diag23b ",0,1)
      return
      end
! *************************************************
      subroutine diag28a_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("t1c ","hp  ","s2  ",111)
      call mult("ta15","s2  ","s1  ",311)
      call add("tb2 ",1.d0,"s1  ",-1.d0,"tb2 ")
      call prt("s3  ","-diag28a",3)
      call prtime("diag28a ",0,1)
      return
      end
! ***************************************************
      subroutine diag26a_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl 26a, 16b
      call reorder("tcm ","s1  ",3412,5)
      call reorder("vhhp","s2  ",1432,5)
      call mult("s1  ","s2  ","s3  ",222)
!   diag 16b
      call seltyp(3333,10,"s2  ",5)
      call add("s3  ",1.d0,"s2  ",1.d0,"s3  ")
! --------------------------------------
      call mult("s3  ","t1c ","s2  ",311)
      call reorder("s2  ","s3  ",3412,5)
      call add("vv0 ",1.d0,"s3  ",1.d0,"vv0 ")
      call prt("s3  ","-diag26a",3)
      call prtime("diag26a ",0,1)
      return
      end
! ************************************************
      subroutine pdgt1p4
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     sub-2*t1**4 compensates for overcounting.
      call reorder("hhpp","s2  ",3124,5)
      call mult("t1c ","s2  ","s1  ",131)
      call mult("t1cr","s1  ","s2  ",131)
      call mult("t1cr","s2  ","p1p4",131)
      call prt("p1p4","p1p4    ",3)
      call prtime("p1p4    ",0,1)
      return
      end
! **********************************************************
      subroutine dgt1p4_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("t1v ","p1p4","s2  ",131)
      call add("tb2 ",1.d0,"s2  ",-2.d0,"tb2 ")
      call prt("s2  ","-.5t1**4",3)
      call prtime("dgt1p4  ",0,1)
      return
      end
! *************************************************
      subroutine dg1n2_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add("vd14",2.d0,"vd12",-1.d0,"s2  ")
      call mult("tm  ","s2  ","s3  ",113)
      call add("t1nw",1.d0,"s3  ",1.d0,"t1nw")
      call prt("s3  ","dg1n2   ",3)
      call prtime("dg1n2   ",0,1)
      return
      end
! ***************************************************
      subroutine dg1n3_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!   incl 3 + 4
      call add("ta1 ",2.d0,"tm2 ",-1.d0,"s2  ")
      call mult("s2  ","hp  ","s3  ",202)
      call add("t1nw",1.d0,"s3  ",1.d0,"t1nw")
      call prt("s3  ","dg1n3   ",3)
      call prtime("dg1n3   ",0,1)
      return
      end
! ****************************************************
      subroutine dg1n5_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("t1v ","pp  ","s2  ",111)
      call add("t1nw",1.d0,"s2  ",1.d0,"t1nw")
      call prt("s2  ","dg1n5   ",3)
      call prtime("dg1n5   ",0,1)
      return
      end
! *****************************************************
      subroutine dg1n6_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add("vhhp",-2.d0,"vvc1",1.d0,"s2  ")
      call reorder("tcm ","s1  ",3214,5)
      call mult("s2  ","s1  ","s3  ",113)
      call add("vv10",1.d0,"s3  ",1.d0,"vv10")
      call prt("s3  ","dg1n6   ",3)
      call prtime("dg1n6   ",0,1)
      return
      end
! ********************************************************
      subroutine dg1n7_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call add("vva1",2.d0,"vphp",-1.d0,"s2  ")
      call mult("s2  ","t1c ","s3  ",202)
      call add("vv10",1.d0,"s3  ",1.d0,"vv10")
      call prt("s3  ","dg1n7   ",3)
      call prtime("dg1n7   ",0,1)
      return
      end
! *********************************************
      subroutine dg1n8_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("vh  ","t1cr","s2  ",111)
      call add("vv10",1.d0,"s2  ",-1.d0,"vv10")
      call prt("s2  ","dg1n8   ",3)
      call prtime("dg1n8   ",0,1)
      return
      end
! **********************************************
      subroutine dg1n9_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add("hhpp",-2.d0,"va6 ",1.d0,"s2  ")
      call mult("tv  ","s2  ","s3  ",113)
      call mult("s3  ","t1cr","s2  ",111)
      call add("t1nw",1.d0,"s2  ",1.d0,"t1nw")
      call prt("s2  ","dg1n9   ",3)
      call prtime("dg1n9   ",0,1)
      return
      end
! ***********************************************
      subroutine pdg1n10
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add("va4 ",-2.d0,"va8 ",1.d0,"s2  ")
      call mult("tcm4","s2  ","pdg1n10 ",113)
      call prt("pdg1n10 ","pdg1n10 ",3)
      call prtime("pdg1n10 ",0,1)
      return
      end
! ***************************************************
      subroutine dg1n10_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call mult("t1v ","pdg1n10 ","s2  ",111)
      call add("t1nw",1.d0,"s2  ",1.d0,"t1nw")
      call prt("s2  ","dg1n10  ",3)
      call prtime("dg1n10  ",0,1)
      return
      end
! ****************************************************
      subroutine pdg1n11
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

!   11 & 12
      call add("va1 ",-2.d0,"va2 ",1.d0,"s2  ")
      call mult("s2  ","t1c ","p111",202)
      call prt("p111","p111    ",3)
      call prtime("p111    ",0,1)
      return
      end
! **********************************************************
      subroutine dg1n11_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call add("ta1 ",-2.d0,"ta2 ",1.d0,"s2  ")
      call mult("s2  ","p111","s3  ",202)
      call add("t1nw",1.d0,"s3  ",1.d0,"t1nw")
      call prt("s3  ","dg1n11  ",3)
      call prtime("dg1n11  ",0,1)
      return
      end
! **************************************************
      subroutine dfold1_csv
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

      call reorder("tv  ","s2  ",2341,5)
      call mult("tvf ","s2  ","s3  ",131)
      call add("tnew",1.d0,"s3  ",-1.d0,"tnew")
      call prt("s3  ","dfold1  ",3)
! --------------------------
      call mult("tvf ","t1or","s2  ",111)
      if (iopt(9).ne.0) call fold_ccsd("s2  ","s3  ")
      call add("t1nw",1.d0,"s2  ",-1.d0,"t1nw")
      call prt("s2  ","dfold1   ",3)
      if (iopt(17).eq.3) then
!  diag f1
        call reord7   ("t31 ","s1  ",231456)
        call reord6   ("s1  ","s2  ",124563)
        call mult6    ("s2  ","tvf ","s1  ",5)
        call reord6   ("s1  ","s2  ",126345)
        call reord7   ("s2  ","s1  ",312456)
        call dveps6   ("s1  ")
        call add      ("s1  ",-1.d0,"t3nw",1.d0,"t3nw")
      end if
      call prt("s2  ","dfold1  ",3)
      call prtime("dfold1  ",0,1)
      return
      end
!  *********************************************************************
      subroutine heff1_csv
!-------------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use rsp1_mod
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)
#include "cfiles.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij2.h"
!     common /rsp1/ nrsp1,irsp1(100)
#include "rsp1.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      !integer*4 ipvt
      integer :: ipvt

!!!!!!!!!!!!!!!!!!! needed dynamic allocation !!!!!!!!!!!!!!!!!!!!!!!
      dimension a(100,100),ipvt(100),z(100),iord(100)

      if (iopt(9).eq.0) return
      call prt("tvf ","oldtvf  ",1)
      call files("tvf ",ntapef,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) ' tvf not found in heff1'
        call abext('heff1')
      endif
      nsize=ije1(nij1)
      istrt11=nstrt(1,ntapef-1)+2*nij1+4
      istrt12=nstrt(2,ntapef-1)
      call dimx(nsize,100,"heff1   ")
      call dcopy(100*100,0.d0,0,a,1)
      do 20 i=1,nsize
20    a(i,i)=1.
      call files("t1v ",ntapev,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t1v not found in heff1'
        call abext('heff1')
      endif
      nsizev=ije2(nij2)
      istrtv1=nstrt(1,ntapev-1)+2*nij2+4
      istrtv2=nstrt(2,ntapev-1)
      do 30 i=1,nsize
      do 30 j=1,nsize
      if (ishft(igen(istrt11+i),-8).eq.                                 &
  ishft(igen(istrt11+j),-8)) then
        itt=ior(ishft(iand(igen(istrt11+j),z'ff'),8),                   &
    iand(igen(istrt11+i),z'ff'))
        a(i,j)=a(i,j)+find(itt,gen(istrtv2+1),igen(istrtv1+1),          &
     nsizev)
      endif
30    continue
!   delete p to r transitions
      io=0
      do 60 i=1,nsize
      do 55 ir=1,nrsp1
      if (irsp1(ir).eq.igen(istrt11+i)) go to 60
55    continue
      io=io+1
      iord(io)=i
60    continue
      nsizeo=io
      do 70 io=1,nsizeo
      gen(istrt12+io)=gen(istrt12+iord(io))
      do 70 jo=1,nsizeo
70    a(io,jo)=a(iord(io),iord(jo))
      call dgeco(a,100,nsizeo,ipvt,rcond,z)
      if (rcond+1.0.eq.1.0) then
        write (6,*) ' singular lin eqs in heff1'
        call abext('heff1')
      endif
      call dgesl(a,100,nsizeo,ipvt,gen(istrt12+1),0)
!   expand back to original tf
      do 75 i=1,nsizeo
75    z(i)=gen(istrt12+i)
      call dcopy(nsize,0.,0,gen(istrt12+1),1)
      do 80 io=1,nsizeo
80    gen(istrt12+iord(io))=z(io)
      call prt("tvf ","newtvf  ",1)
      return
      end

!  *********************************************************************
      subroutine rdbox_csv
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

      call rdfile (20,0,"bx1 ",-1)
      call rdfile (20,1,"bx2 ",-1)
      call rdfile (20,2,"bx3 ",-1)
      call rdfile (20,3,"bx4 ",-1)
      call rdfile (20,4,"bx5 ",-1)
      call rdfile (20,5,"bx6 ",-1)
      call rdfile (20,6,"bx7 ",-1)
      call rdfile (20,7,"bx8 ",-1)
      if (iopt(17).eq.5) then
        call rdfile (20,8,"bx9 ",-1)
        call rdfile (20,9,"bx10",-1)
        call rdfile (20,10,"bx11",-1)
      end if
      return
      end 
!  *******************************************************************
      subroutine box12_csv (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      if (iparm.eq.0) then
!  diag 1
        call move     ("vppp","s3  ")
!  diag 5
        call reorder  ("vhhp","s1  ",1432,4)
        call mult     ("s1  ","tcm4","s2  ",222)
        call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!  diag 8
        call reorder  ("vhpp","s1  ",1432,4)
        call mult     ("s1  ","t1cr","s2  ",311)
        call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
!  diag 9
        call reorder  ("vphp","s1  ",1243,4)
        call mult     ("s1  ","t1cr","s2  ",311)
        call reorder  ("s2  ","s1  ",1243,4)
        call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
        call reorder  ("s3  ","bx12",1342,4)
      else if (iparm.eq.1) then
!  diag 6
        call cmult    ("t1c ",0.5d0,"s1  ")
        call t1t2     ("s1  ","tv  ","t31 ","s2  ",123456,2)
        call t1t2     ("t1c ","tv  ","s2  ","s1  ",123456,3)
        call reord6   ("s1  ","s2  ",124356)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","s2  ",126543)
        call add      ("va4 ",-2.d0,"va8 ",1.d0,"s1  ")
        call mult6    ("s2  ","s1  ","s3  ",2)
! v1 Feb 95
        call mult     ("ph  ","ta15","s1  ",131)
        call reorder  ("s1  ","s2  ",2134,4)
        call add      ("s3  ",1.d0,"s2  ",-.5d0,"s3  ")
!  diag 2
        call add      ("vd7 ",2.d0,"vd8 ",-1.d0,"s1  ")
        call mult     ("ta1 ","s1  ","s2  ",222)
        call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!  diag 3
        call mult     ("tm2 ","vd7 ","s1  ",222)
        call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!  diag 4
        call mult     ("tm2 ","vd8 ","s1  ",222)
        call reorder  ("s1  ","s2  ",1324,4)
        call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
!  diag 7
        call t1t2     ("t1v ","tc  ","t31 ","s1  ",123456,1)
        call t1t2     ("t1c ","tv  ","s1  ","s2  ",123456,2)
        call t1t2     ("t1c ","tm  ","s2  ","s1  ",123456,3)
        call reord6   ("s1  ","s2  ",125346)
        call reord7   ("s2  ","s1  ",132456)
        call reord6   ("s1  ","s2  ",126543)
        call mult6    ("s2  ","va4 ","s1  ",2)
        call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!
        call move     ("s3  ","bx12")
      end if
      return
      end
!  *******************************************************************
      subroutine box13_csv (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      if (iparm.eq.0) then
!  diag 1
        call cmult    ("vvc1",-1.d0,"s3  ")
!  diag 2
        call add      ("vvc1",-2.d0,"vhhp",1.d0,"s1  ")
        call mult     ("s1  ","tc1 ","s2  ",222)
        call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!  diag 3
        call mult     ("vvc1","tcm2","s1  ",222)
        call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!  diag 9
        call mult     ("vphp","t1c ","s1  ",311)
        call reorder  ("s1  ","s2  ",1342,4)
        call add      ("s2  ",-1.d0,"s3  ",1.d0,"s3  ")
!  diag 10
        call mult     ("vhhh","t1cr","s1  ",311)
        call reorder  ("s1  ","s2  ",1324,4)
        call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!
        call reorder  ("s3  ","bx13",1342,4)
      else if (iparm.eq.1) then
!  diag 6
        call cmult    ("t1c ",0.5d0,"s1  ")
        call t1t2     ("s1  ","tv  ","t31 ","s2  ",123456,3)
        call t1t2     ("t1v ","tc  ","s2  ","s1  ",123456,1)
        call reord6   ("s1  ","s2  ",125346)
        call add      ("hhpp",-2.d0,"va6 ",1.d0,"s1  ")
        call mult6    ("s2  ","s1  ","s3  ",2)
! v1 Feb 95
        call mult     ("ta6 ","hp  ","s1  ",311)
        call add      ("s3  ",1.d0,"s1  ",-.5d0,"s3  ")
!  diag 4
        call mult     ("tm2 ","hhhp","s1  ",222)
        call reorder  ("s1  ","s2  ",1324,4)
        call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!  diag 5
        call mult     ("tm  ","vd12","s1  ",222)
        call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!  diag 7
        call t1t2     ("t1v ","tc  ","t31 ","s1  ",123456,1)
        call t1t2     ("t1c ","tv  ","s1  ","s2  ",123456,2)
        call t1t2     ("t1c ","tm  ","s2  ","s1  ",123456,3)
        call reord6   ("s1  ","s2  ",126345)
        call mult6    ("s2  ","hhpp","s1  ",2)
        call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!  diag 8
        call mult     ("t1v ","va1 ","s1  ",131)
        call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!
        call move     ("s3  ","bx13") 
      end if
      return
      end 
!  *******************************************************************
      subroutine box14_csv (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      if (iparm.eq.0) then
!  diag 1
        call cmult    ("vhph",-1.d0,"bx14")
!  diag 4
        call mult     ("vhhp","tcm2","s1  ",222)
        call reorder  ("s1  ","s2  ",1342,4)
        call add      ("s2  ",1.d0,"bx14",1.d0,"bx14")
!  diag 8
        call mult     ("vhpp","t1c ","s1  ",311)
        call reorder  ("s1  ","s2  ",1432,4)
        call add      ("s2  ",-1.d0,"bx14",1.d0,"bx14")
!  diag 10
        call reorder  ("vhhh","s1  ",1243,4)
        call mult     ("s1  ","t1cr","s2  ",311)
        call reorder  ("s2  ","s1  ",1243,4)
        call add      ("s1  ",1.d0,"bx14",1.d0,"bx14")
      else if (iparm.eq.1) then
!  diag 6
        call cmult    ("t1c ",0.5d0,"s1  ")
        call t1t2     ("s1  ","tv  ","t31 ","s2  ",123456,3)
        call t1t2     ("t1c ","tv  ","s2  ","s1  ",123456,2)
        call reord6   ("s1  ","s2  ",124356)
        call add      ("hhpp",-2.d0,"va6 ",1.d0,"s1  ")
        call mult6    ("s2  ","s1  ","s3  ",2)
! v1 Feb 95
        call mult     ("tv  ","hp  ","s1  ",311)
        call add      ("s3  ",1.d0,"s1  ",-.5d0,"s3  ")
!  diag 2
        call add      ("vc9 ",-2.d0,"vc10",1.d0,"s1  ")
        call mult     ("ta1 ","s1  ","s2  ",222)
        call reorder  ("s2  ","s1  ",1324,4)
        call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!  diag 3
        call mult     ("tm2 ","vc9 ","s1  ",222)
        call reorder  ("s1  ","s2  ",1324,4)
        call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!  diag 5
        call mult     ("tm  ","vd14","s1  ",222)
        call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!  diag 7
        call t1t2     ("t1v ","tc  ","t31 ","s1  ",123456,1)
        call t1t2     ("t1c ","tv  ","s1  ","s2  ",123456,2)
        call t1t2     ("t1c ","tm  ","s2  ","s1  ",123456,3)
        call reord6   ("s1  ","s2  ",126354)
        call mult6    ("s2  ","hhpp","s1  ",2)
        call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!  diag 9
        call mult     ("t1v ","hphp","s1  ",131)
        call add      ("s1  ",-1.d0,"s3  ",1.d0,"s3  ")
!
        call move     ("s3  ","bx14")  
      end if
      return
      end
!  *******************************************************************
      subroutine box15_csv (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      if (iparm.eq.0) then
!  diag 1
        call add      ("vva1",2.d0,"vphp",-1.d0,"bx15")
!  diag 3
        call reorder  ("vhhp","s1  ",1243,4)
        call reorder  ("vhhp","s2  ",1342,4)
        call add      ("s1  ",-2.d0,"s2  ",1.d0,"s1  ")
        call mult     ("s1  ","t1cr","s2  ",311)
        call reorder  ("s2  ","s1  ",1423,4)
        call add      ("s1  ",1.d0,"bx15",1.d0,"bx15")
      else if (iparm.eq.1) then
!  diag 2
        call add      ("vd12",2.d0,"vd14",-1.d0,"s1  ")
        call mult     ("t1v ","s1  ","s3  ",131)
!  diag 4
        call add      ("va1 ",2.d0,"va2 ",-1.d0,"s1  ")
        call mult     ("ta1 ","s1  ","s2  ",222)
        call add      ("s2  ",2.d0,"s3  ",1.d0,"s3  ")
!  diag 5
        call add      ("va1 ",-2.d0,"va2 ",1.d0,"s1  ")
        call mult     ("tm2 ","s1  ","s2  ",222)
        call add      ("s2  ",1.d0,"s3  ",1.d0,"s3  ")
!
        call move     ("s3  ","bx15") 
      end if
      return
      end
!  *******************************************************************
      subroutine box16_csv (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      if (iparm.eq.0) then
!  diag 1
        call cmult    ("vva1",-1.d0,"bx16")
!  diag 3
        call reorder  ("vhhp","s1  ",1243,4)
        call mult     ("s1  ","t1cr","s2  ",311)
        call reorder  ("s2  ","s1  ",1423,4)
        call add      ("s1  ",1.d0,"bx16",1.d0,"bx16")
      else if (iparm.eq.1) then
!  diag 2
        call mult     ("t1v ","vd12","s3  ",131)
!  diag 4
        call add      ("va1 ",-2.d0,"va2 ",1.d0,"s1  ")
        call mult     ("ta1 ","s1  ","s2  ",222)
        call add      ("s2  ",1.d0,"s3  ",-1.d0,"s3  ")
!  diag 5
        call mult     ("tm2 ","va1 ","s1  ",222)
        call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!
        call move     ("s3  ","bx16")  
      end if
      return
      end
!  *******************************************************************
      subroutine box17_csv (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      if (iparm.eq.0) then
!  diag 1
        call cmult    ("vphp",-1.d0,"bx17")
!  diag 3
        call reorder  ("vhhp","s1  ",1432,4)
        call mult     ("s1  ","t1cr","s2  ",311)
        call reorder  ("s2  ","s1  ",1432,4)
        call add      ("s1  ",1.d0,"bx17",1.d0,"bx17")
      else if (iparm.eq.1) then
!  diag 2
        call mult      ("t1v ","vd14","s3  ",131)
!  diag 4
        call mult      ("tm2 ","va2 ","s1  ",222)
        call add       ("s1  ",1.d0,"s3  ",-1.d0,"s3  ")
!
        call move      ("s3  ","bx17")  
      end if
      return
      end 
!  ******************************************************************
      subroutine box18_csv (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      if (iparm.eq.0) then
!  diag 1
        call move     ("vhhh","bx18")
!  diag 2a
        call mult     ("vhhp","t1c ","s1  ",311)
        call reorder  ("s1  ","s2  ",1432,4)
        call add      ("s2  ",1.d0,"bx18",1.d0,"bx18")
      else if (iparm.eq.1) then
!  diag 2b
        call mult     ("t1v ","hhhp","s3  ",131)
!  diag 3
        call mult     ("tm  ","hhpp","s1  ",222)
        call add      ("s1  ",1.d0,"s3  ",1.d0,"s3  ")
!
        call move     ("s3  ","bx18")  
      end if
      return
      end
!  ******************************************************************
      subroutine box19_csv (iparm)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      if (iparm.eq.0) then
!  diag 1
        call add      ("vvc1",-2.d0,"vhhp",1.d0,"s1  ")
        call mult     ("s1  ","t1c ","bx19",202)
! v1 Feb 95
        call add      ("bx19",1.d0,"vh  ",-1.d0,"bx19")

      else if (iparm.eq.1) then
!  diag 2
        call add      ("hhpp",-2.d0,"va6 ",1.d0,"s1  ")
        call mult     ("tm  ","s1  ","s3  ",113)
! v1 Feb 95
        call mult     ("t1v ","hp  ","s1  ",111)
        call add      ("s3  ",1.d0,"s1  ",-1.d0,"s3  ")

        call move     ("s3  ","bx19")  
      end if
      return
      end
!  *******************************************************************
      subroutine wrtbox_csv
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      call box12_csv(0)
      call wrtfile  (21,0,"bx12")
      call box12_csv(1)
      call reads    (21,0,"s1  ")
      call add      ("s1  ",1.d0,"bx12",1.d0,"bx12")
      call wrtfile  (21,0,"bx12")
!
      call box13_csv(0)
      call wrtfile  (21,1,"bx13")
      call box13_csv(1)
      call reads    (21,1,"s1  ")
      call add      ("s1  ",1.d0,"bx13",1.d0,"bx13")
      call wrtfile  (21,1,"bx13")
!
      call box14_csv(0)
      call wrtfile  (21,2,"bx14")
      call box14_csv(1)
      call reads    (21,2,"s1  ")
      call add      ("s1  ",1.d0,"bx14",1.d0,"bx14")
      call wrtfile  (21,2,"bx14")
!
      call box15_csv(0)
      call wrtfile  (21,3,"bx15")
      call box15_csv(1)
      call reads    (21,3,"s1  ")
      call add      ("s1  ",1.d0,"bx15",1.d0,"bx15")
      call wrtfile  (21,3,"bx15")
!
      call box16_csv(0)
      call wrtfile  (21,4,"bx16")
      call box16_csv(1)
      call reads    (21,4,"s1  ")
      call add      ("s1  ",1.d0,"bx16",1.d0,"bx16")
      call wrtfile  (21,4,"bx16")
!
      call box17_csv(0)
      call wrtfile  (21,5,"bx17")
      call box17_csv(1)
      call reads    (21,5,"s1  ")
      call add      ("s1  ",1.d0,"bx17",1.d0,"bx17")
      call wrtfile  (21,5,"bx17")
!
      call box18_csv(0)
      call wrtfile  (21,6,"bx18")
      call box18_csv(1)
      call reads    (21,6,"s1  ")
      call add      ("s1  ",1.d0,"bx18",1.d0,"bx18")
      call wrtfile  (21,6,"bx18")
!
      call box19_csv(0)
      call wrtfile  (21,7,"bx19")
      call box19_csv(1)
      call reads    (21,7,"s1  ")
      call add      ("s1  ",1.d0,"bx19",1.d0,"bx19")
      call wrtfile  (21,7,"bx19")
!
      return
      end
!  **********************************************************************
