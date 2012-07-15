!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine clear(name,iii,iiv)
!--------------------------------------------------------------------------------------
!
!
!
!
!--------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
#include "implicit.h"
#include "priunit.h"
!     implicit real*8 (a-h,o-z)
!  zero out file (create if does not exist)
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)
#include "cij1.h"
!      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#endif

      call tmplt(iii,iiv,name,ntape)
      nsize=ije1(nij1)

!      call dcopy(nsize,0.d0,0,
!     1 gen(nstrt(2,ntape-1)+1),1)
      call dcopy(nsize,0.d0,0,gen(nstrt(2,ntape-1)+1),1)

      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine clear6 (name,iii,iiv)
!--------------------------------------------------------------------------------------
!
!
!
!
!--------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"

#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#endif

      call tmplt6 (iii,iiv,name,ntape)
      nsize=ije1(nij1)

      call dcopy (nsize,0.d0,0,gen(nstrt(2,ntape-1)+1),1)
      call prt6 (name,name,2)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine clear8 (name,iii,iiv)
!--------------------------------------------------------------------------------------
!
!
!
!
!--------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#endif

      call tmplt8 (iii,iiv,name,ntape)
      nsize=ije1(nij1)
      call dcopy (nsize,0.d0,0,gen(nstrt(2,ntape-1)+1),1)
      call prt8 (name,name,2)
      return
      end

      function find(itt,t,it,nt)
      implicit real*8 (a-h,o-z)
      dimension t(nt),it(nt)
      find=0.d0
      do 10 i=1,nt
      if (it(i).eq.itt) then
        find=t(i)
        go to 90
      endif
10    continue
90    return
      end


!  *******************************************************************
      function find2 (j1,j2) 
!----------------------------------------------------------------------
! find and returns 2-index quantity element
!----------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use gen_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#include "strt.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!      common  // gen(nbig)
!      common /i/ igen(nbig)
#include "gen.h"
!      common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!     1 istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!     2 nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
#include "strt.h"
#endif

      find2=0.d0
      id=ior(ishft(j1,8),j2)
      do 10 i=1,nsize1
        if (id.eq.igen(istrt11+i)) then
          find2=gen(istrt12+i)
        end if
10    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function find4(j1,j2,j3,j4)
!--------------------------------------------------------------------------------------
! find and returns 4-index quantity element
!--------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use gen_mod
#include "implicit.h"
#include "priunit.h"
#include "strt.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
#include "gen.h"
!      common  // gen(nbig)
!      common /i/ igen(nbig)
#include "strt.h"
!      common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!     1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!     2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
#endif

      find4=0.d0
      id=ior(ior(ishft(j1,24),ishft(j2,16)),ior(ishft(j3,8),j4))
      do 10 i=1,nsize2
        if (id.eq.igen(istrt21+i)) then
          find4=gen(istrt22+i)
        end if
10    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function find6 (j1,j2,j3,j4,j5,j6)
!--------------------------------------------------------------------------------------
! find and returns 6-index quantity element
!--------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use gen_mod
      use cij3_mod
#include "implicit.h"
#include "priunit.h"
#include "strt.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!      common  // gen(nbig)
!      common /i/ igen(nbig)
#include "gen.h"
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
!      common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!     1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!     2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
#include "strt.h"
#endif

      find6=0.d0
      id=ior(ishft(j1,16),ior(ishft(j2,8),j3))
      ind=ior(ishft(j4,16),ior(ishft(j5,8),j6))
      do 10 i=1,nij3
        if (id.eq.ij3(i)) then
          isave=i
        end if
10    continue
      lngth=ije3(isave)-ije3(isave-1)
      do 20 i=1,lngth
        ii=igen(istrt31+ije3(isave-1)+i)
        if (ind.eq.iand(ii,z'ffffff')) then
          find6=gen(istrt32+ije3(isave-1)+i)
        end if
20    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function find8 (j1,j2,j3,j4,j5,j6,j7,j8)
!--------------------------------------------------------------------------------------
! return particular array element (as real) according to his (entering) 8 indexes j1-j8
!
! 8 indexes mean quadruple excitations, like in sectors (40) or (04)
!--------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use gen_mod
      use cij4_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#include "strt.h"
#else
!      implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!     common  // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /cij4/ nij4,ij4(ndij),ije4(0:ndij),iii4,iiv4
#include "cij4.h"
!     common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!    1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!    2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
#include "strt.h"
#endif
      logical :: isave_set=.false.

      find8=0.d0
      isave = 0 ! mi-necessary initialization, otherwise might indicate ifort runtime error
      id=ior(ior(ishft(j1,24),ishft(j2,16)),ior(ishft(j3,8),j4))
      ind=ior(ior(ishft(j5,24),ishft(j6,16)),ior(ishft(j7,8),j8))
      do 10 i=1,nij4
        if (id.eq.ij4(i)) then
          isave=i
          isave_set=.true.
        end if
10    continue

!mi: problem - if isave not set, initial isave=0 causes negative index in ije4 !
      if (.not.isave_set) then
        write(LUPRI,*) 'find8: j1,j2,j3,j4:',j1,j2,j3,j4
        write(LUPRI,*) 'find8: j5,j6,j7,j8:',j5,j6,j7,j8
        write(LUPRI,*) 'nij4=',nij4,' id=',id,' ind=',ind
        call quit('find8: isave not set, must stop')
      endif
      ! if isave=0, we get index -1 what is out-of-bounds
      lngth=ije4(isave)-ije4(isave-1)
      do 20 i=1,lngth
        if (ind.eq.igen(istrt41+ije4(isave-1)+i)) then
          find8=gen(istrt42+ije4(isave-1)+i)
        endif
20    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function find10 (ja,jb,jc,jd,je,jf,jg,jh,ji,jj)
!--------------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------------
      use ipack_mod, only: ipack5
#if defined HSFS_DYN_ALLOC
      use gen_mod
      use cij5_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#include "strt.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "param_ccsd.h"
!     common  // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /cij5/ nij5,ij5(ndij),ije5(0:ndij),iii5,iiv5
#include "cij5.h"
!     common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!    1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!    2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
#include "strt.h"
#endif
      dimension jj1(10),kk(10)

      find10=0.d0
      jj1(1)=ja
      jj1(2)=jb
      jj1(3)=jc
      jj1(4)=jd
      jj1(5)=je
      jj1(6)=jf
      jj1(7)=jg
      jj1(8)=jh
      jj1(9)=ji
      jj1(10)=jj
      do 10 i=1,10
        kk(i)=ifndval(jj1(i))
10    continue
      idb=ipack5(kk(1),kk(2),kk(3),kk(4),kk(5))
      ind=ipack5(kk(6),kk(7),kk(8),kk(9),kk(10))
      do 20 i=1,nij5
        if (idb.eq.ij5(i)) then
          isave=i
        end if
20    continue
      lngth=ije5(isave)-ije5(isave-1)
      do 30 i=1,lngth
        if (ind.eq.igen(istrt51+ije5(isave-1)+i)) then
          find10=gen(istrt52+ije5(isave-1)+i)
        end if
30    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function find12 (ja,jb,jc,jd,je,jf,jg,jh,ji,jj,jk,jl)
!--------------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------------
      use ipack_mod, only: ipack6
#if defined HSFS_DYN_ALLOC
      use gen_mod
      use cij6_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#include "strt.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!     common  // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /cij6/ nij6,ij6(ndij),ije6(0:ndij),iii6,iiv6
#include "cij6.h"
!     common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!    1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!    2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
#include "strt.h"
#endif
      dimension jj1(12),kk(12)

      find12=0.d0
      jj1(1)=ja
      jj1(2)=jb
      jj1(3)=jc
      jj1(4)=jd
      jj1(5)=je
      jj1(6)=jf
      jj1(7)=jg
      jj1(8)=jh
      jj1(9)=ji
      jj1(10)=jj
      jj1(11)=jk
      jj1(12)=jl
      do 10 i=1,12
        kk(i)=ifndval(jj1(i))
10    continue
      idb=ipack6(kk(1),kk(2),kk(3),kk(4),kk(5),kk(6))
      ind=ipack6(kk(7),kk(8),kk(9),kk(10),kk(11),kk(12))
      do 20 i=1,nij6
        if (idb.eq.ij6(i)) then
          isave=i
        end if
20    continue
      lngth=ije6(isave)-ije6(isave-1)
      do 30 i=1,lngth
        if (ind.eq.igen(istrt61+ije6(isave-1)+i)) then
          find12=gen(istrt62+ije6(isave-1)+i)
        end if
30    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine fold_ccsd(name2,namef)
!-------------------------------------------------------------------------------
!
!  two-electron part modified to include chp, chq.
!  880922.
!  all four indices now checked. incomplete space implemented
!  for chp, chq only for exclusion of whole symmetries.
!
!  Renamed from "subroutine fold(name2,namef)" due to name overlaping with relccsd routine
!
!-------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use rsp1_mod
      use rsp2_mod
      use dat_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
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
!     common /rsp1/  nrsp1,irsp1(100)
#include "rsp1.h"
!     common /rsp2/  nrsp2,irsp2(500)
#include "rsp2.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      call files_ccsd(name2,ntape2,new,2,-1,-1)
      iiv1=11
      if (iii2.gt.100) iiv1=1111
      iii1=iii2
      call files_ccsd(namef,ntapef,new,1,iii1,iiv1)
      nsize1=ije1(nij1)
      istrtf1=nstrt(1,ntapef-1)+2*nij1+4
      istrtf2=nstrt(2,ntapef-1)

      call dcopy(nsize1,0.d0,0,gen(istrtf2+1),1)

      nsize2=ije2(nij2)
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      if (iii1.lt.100) then
!   1-electron t.
        do 120 i=1,nsize2
        !i1=iand(igen(istrt21+i),'ff'X)
        i1=iand(igen(istrt21+i),z'ff')
        if (icv(i1).eq.0) go to 120
        if (iopt(9).ne.0) then
          do 105 ir=1,nrsp1
105       if (igen(istrt21+i).eq.irsp1(ir)) go to 120
      endif
!   move p to p excitations to folded
        do 110 j=1,nsize1
        if (igen(istrt21+i).eq.igen(istrtf1+j)) then
          gen(istrtf2+j)=gen(istrt22+i)
          gen(istrt22+i)=0.d0
          go to 120
        endif
110     continue
        write (6,*) ' folded not found for ',igen(istrt21+i)
        call abext('fold    ')
120     continue
      else
!   2-electron t
        do 180 i=1,nsize2
!   transition from r space eliminated
        if (iopt(9).ne.0) then
         !i_4hthpf=4hthpf
         i_4hthpf=transfer("thpf",i_4hthpf)
         !i_4hthqf=4hthqf
         i_4hthqf=transfer("thqf",i_4hthqf)

!        if (namef.eq.4hthpf) then
         if (namef.eq.i_4hthpf) then
          i12=ior(iand(ishft(igen(istrt21+i),-16),z'ff00'),             &
               iand(igen(istrt21+i),z'ff'))
         else if (namef.eq.i_4hthqf) then
          i12=ior(iand(ishft(igen(istrt21+i),-16),z'ff00'),             &
               iand(ishft(igen(istrt21+i),-8),z'ff'))
         else
           i12=iand(ishft(igen(istrt21+i),-16),z'ffff')
          endif
          do 130 ir=1,nrsp2
          if (irsp2(ir).eq.i12) then
            gen(istrt22+i)=0.d0
            go to 180
          endif
130       continue
        endif
!   move p to p transitions to folded
        i1=iand(ishft(igen(istrt21+i),-24),z'ff')
        i2=iand(ishft(igen(istrt21+i),-16),z'ff')
        i3=iand(ishft(igen(istrt21+i),-8) ,z'ff')
        i4=iand(      igen(istrt21+i)     ,z'ff')
        if ((icv(i1).eq.0).or.(icv(i2).eq.0)) go to 180
        if ((icv(i3).eq.0).or.(icv(i4).eq.0)) go to 180

        if (iopt(9).eq.0) go to 170
        !i_4hthpf=4hthpf
        i_4hthpf=transfer("thpf",i_4hthpf)
        !i_4hthqf=4hthqf
        i_4hthqf=transfer("thqf",i_4hthqf)
        if (namef.eq.i_4hthpf) then
          i34=iand(ishft(igen(istrt21+i),-8),z'ffff')
        else if (namef.eq.i_4hthqf) then
          i34=ior(iand(ishft(igen(istrt21+i),-8),z'ff00'),              &
             iand(igen(istrt21+i),z'ff'))
        else
          i34=iand(igen(istrt21+i),z'ffff')
        endif
        do 150 ir=1,nrsp2
150     if (i34.eq.irsp2(ir)) go to 180
170     do 175 j=1,nsize1
        if (igen(istrt21+i).eq.igen(istrtf1+j)) then
          gen(istrtf2+j)=gen(istrt22+i)
          gen(istrt22+i)=0.d0
          go to 180
        endif
175     continue
        write (6,*) ' folded not found for ',igen(istrt21+i)
        call abext('fold    ')
180     continue
      endif
      call prt(namef,"folded  ",2)
      call prtime("folded  ",nsize1,1)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine fold4 (nameu,namef,ind,ninds)
!--------------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
!     implicit real*8 (a-h,o-z)
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use dat_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
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
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif
      dimension jj1(4)
      call files_ccsd (nameu,ntapeu,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'fold4 input tape does not exist, ',nameu
        stop                    
      end if
      do 10 i=1,4     
        jj1(i)=0
10    continue
      if (ninds.eq.2) then
        if (iiv2.ne.0) then
!        mi: inconsitent number of routine parameters !!!!!
          !call parse (iiv2,jj1(1),jj1(2))
          call parse (iiv2,jj1(1),jj1(2),0,0)
        end if
        jj1(ind)=1
        iiv1=10*jj1(1)+jj1(2)
      else if (ninds.eq.4) then
        if (iiv2.ne.0) then
          call parse (iiv2,jj1(1),jj1(2),jj1(3),jj1(4))
        end if
        jj1(ind)=1
        iiv1=1000*jj1(1)+100*jj1(2)+10*jj1(3)+jj1(4)
      end if
      iii1=iii2
      call tmplt (iii1,iiv1,namef,ntapef)
      nsizef=ije1(nij1)
      nsizeu=ije2(nij2)
      istrtf1=nstrt(1,ntapef-1)+2*nij1+4
      istrtf2=nstrt(2,ntapef-1)
      istrtu1=nstrt(1,ntapeu-1)+2*nij2+4
      istrtu2=nstrt(2,ntapeu-1)
      call dcopy (nsizef,0.d0,0,gen(istrtf2+1),1)
      j=1
      ind2=igen(istrtu1+j)
      do 20 i=1,nsizef
        ind1=igen(istrtf1+i)
        ididit=0
        do while (ididit.eq.0)
          if (ind1.eq.ind2) then
            gen(istrtf2+i)=gen(istrtu2+j)
            ididit=1
          end if
          if (ididit.eq.0) then
            j=j+1
            ind2=igen(istrtu1+j)
          end if
          if (j.gt.nsizeu) then
            write (6,*) 'folded not found - fold4'
            stop                     
          end if
        end do
20    continue
      call prtime ("fold4   ",0,2)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine fold6 (nameu,namef,ind,iparm)
!--------------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use dat_mod
      use hsfs_parameters
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
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif
      dimension jj1(6)
      call files_ccsd (nameu,ntapeu,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'fold6 input tape does not exist, ',nameu
        stop                    
      end if
      if (iparm.eq.1) then
        call parse6 (iiv2,jj1(1),jj1(2),jj1(3),jj1(4),jj1(5),jj1(6))
        jj1(ind)=1
        iiv1=100000*jj1(1)+10000*jj1(2)+1000*jj1(3)                     &
        +100*jj1(4)+10*jj1(5)+jj1(6)
        iii1=iii2
        call tmplt6 (iii1,iiv1,namef,ntapef)
      else if (iparm.eq.2) then
        call files_ccsd (namef,ntapef,new,1,-1,-1)
        if (new.ne.0) then
          write (6,*) 'fold6 input tape does not exist, ',namef
          !stop 
          call quit('fold6 input tape does not exist!')
        end if
      end if
      nsizef=ije1(nij1)
      nsizeu=ije2(nij2)
      istrtf1=nstrt(1,ntapef-1)+2*nij1+4
      istrtf2=nstrt(2,ntapef-1)
      istrtu1=nstrt(1,ntapeu-1)+2*nij2+4
      istrtu2=nstrt(2,ntapeu-1)
      if (iparm.eq.1) then
        call dcopy (nsizef,0.d0,0,gen(istrtf2+1),1)
      end if
      i2=1
      j2=1
      id2=ij2(i2)
      ind2=igen(istrtu1+j2)
      do 10 i1=1,nij1
        id1=ij1(i1)
        do 20 j1=ije1(i1-1)+1,ije1(i1)
          ind1=igen(istrtf1+j1)
          ididit=0
          do while (ididit.eq.0)
            if (id1.eq.id2) then
              if (ind1.eq.ind2) then
                if (iparm.eq.1) then
                  gen(istrtf2+j1)=gen(istrtu2+j2)
                else if (iparm.eq.2) then
                  gen(istrtu2+j2)=gen(istrtu2+j2)-gen(istrtf2+j1)
                end if
                ididit=1
              end if
            end if
            if (ididit.eq.0) then
              if (id1.ne.id2) then
                i2=i2+1
                id2=ij2(i2)
                j2=ije2(i2-1)+1
                ind2=igen(istrtu1+j2)
              else
                j2=j2+1
                ind2=igen(istrtu1+j2)
                if (j2.gt.ije2(i2)) then
                  i2=i2+1
                  id2=ij2(i2)
                end if
              end if
              if (j2.gt.nsizeu) then
                write (6,*) 'folded not found - fold6'
                call quit('folded not found - fold6!')                
              end if
            end if
          end do
20      continue
10    continue
      call prtime ("fold6   ",iparm,2)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine fold8 (nameu,namef,ind,iparm)
!--------------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use dat_mod
      use hsfs_parameters
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
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      dimension ia(8)
      call files_ccsd (nameu,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'fold8 input tape does not exist, nameu=',nameu
        call quit('fold8 input tape nameu does not exist!')
      end if
      if (iparm.eq.1) then
        call parse8(iiv2,ia(1),ia(2),ia(3),ia(4),ia(5),ia(6),ia(7),ia(8))
        ia(ind)=1
        iii1=iii2
        iiv1=10000000*ia(1)+1000000*ia(2)+100000*ia(3)+10000*ia(4)+   &
                 1000*ia(5)+100*ia(6)+10*ia(7)+ia(8)
        call tmplt8 (iii1,iiv1,namef,ntape1)
      else if (iparm.eq.2) then
        call files_ccsd (namef,ntape1,new,1,-1,-1)
        if (new.ne.0) then
          write (6,*) 'fold8 input tape does not exist, namef=',namef
          call quit('fold8 input tape namef does not exist!') 
        end if
      end if
      istrt11=nstrt(1,ntape1-1)+2*nij1+4
      istrt12=nstrt(2,ntape1-1)
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      if (iparm.eq.1) then
        call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1)
      end if
      id2=0
      idc=0
      do 10 i=1,nij1
        id1=ij1(i)
        ind2=0
        indc=0
        do 20 j=1,ije1(i)-ije1(i-1)
          ind1=igen(istrt11+ije1(i-1)+j)
          ididit=0
          do while (ididit.eq.0)
            if (id2.eq.id1) then
              if (ind2.eq.ind1) then
                if (iparm.eq.1) then
                  gen(istrt12+ije1(i-1)+j)=gen(istrt22+ije2(idc-1)+indc)
                else
                  gen(istrt22+ije2(idc-1)+indc)=           &
                       gen(istrt22+ije2(idc-1)+indc)-      &
                       gen(istrt12+ije1(i-1)+j)
                end if
                ididit=1
              else if (ind2.lt.ind1) then
                indc=indc+1
                ind2=igen(istrt21+ije2(idc-1)+indc)
              else
                write (6,*) 'folded not found - fold8 '
                call quit('folded not found - fold8')
              end if
            else if (id2.lt.id1) then
              idc=idc+1
              id2=ij2(idc)
              indc=0
            else
              write (6,*) 'folded not found - fold8 '
              call quit('folded not found - fold8')
            end if
          end do
20      continue
10    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine msub1 (name2,name3,namet)
!--------------------------------------------------------------------------------------
!
!        special type of diagrams
!
!--------------------------------------------------------------------------------------
      use ipack_mod, only: ipack3,ipack4
      use print_level
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use ctmp_mod
      use dpd_mod
      use hsfs_parameters
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
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#include "ctmp.h"
!     common /dpd/ nir,igmt(8,8),ici(200)
#include "dpd.h"
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ---- to be replaced with dynamic allocation in module ! -----
      dimension ia(8),ib(10000),ibe(0:10000)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!mi TODO: include file with generalized function
!     ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)
!     ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

#if defined HSFS_DYN_ALLOC
      call check_modules
#endif

      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub1 input file does not exist, name2=',name2
        call quit('msub1 input file name2 does not exist!')
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub1 input file does not exist, name3=',name3
        call quit('msub1 input file name3 does not exist!')
      end if
      if (((iii2-10*(iii2/10)).ne.(iii3-10*(iii3/10))).or.       &
         ((iiv2-10*(iiv2/10)).ne.(iiv3-10*(iiv3/10)))) then
        write (6,*) 'contracted indices do not match - msub1'
        call quit('contracted indices do not match - msub1!')
      end if
      if ((iii2.lt.222222).or.(iii2.gt.333333)) then
        write (6,*) 'msub1 called with inconsistent iii2=',iii2
        call quit('msub1 called with inconsistent iii2!')
      end if
      if ((iii3.lt.2222).or.(iii3.gt.3333)) then
        write (6,*) 'msub1 called with inconsistent iii3=',iii3
        call quit('msub1 called with inconsistent iii3!')
      end if
      iii1=1000*(iii2/10)+iii3/10
      iiv1=1000*(iiv2/10)+iiv3/10
      call tmplt8 (iii1,iiv1,namet,ntapet)
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)
      call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1)
      if ((ije2(nij2).eq.0).or.(ije3(nij3).eq.0)) then
        return
      end if
      call tchng6 (nij2,ij2,ije2,ntape2,4)
      call tmpcng (nij3,ij3,ije3,ntape3,2,3)
      idc=0
      idold=0 
      ibe(0)=0
      do 10 i=1,nij2
        do 20 j=0,3
          ia(4-j)=iand(ishft(ij2(i),-8*j),z'ff')
20      continue
        ii=0
        n=0
        do 30 j=1,ije2(i)-ije2(i-1)
          jj=iand(ishft(igen(istrt21+ije2(i-1)+j),-8),z'ff')
          if (ii.ne.jj) then
            n=n+1
            ib(n)=jj
            ii=jj
            ibe(n)=ibe(n-1)+1
          else
            ibe(n)=ibe(n)+1
          end if
30      continue
        call dimx (n,10000,'msub1')
        do 50 j=1,n
          ia(5)=ib(j)
          do 60 k=1,nij3 
            do 70 m=0,2
            ia(8-m)=iand(ishft(ij3(k),-8*m),z'ff')
70          continue
            ii=1
            do 80 m=1,8
              ii=igmt(ii,ici(ia(m)))
80          continue
            if (ii.eq.1) then
              lngth=ije3(k)-ije3(k-1)
              if (lngth.ne.(ibe(j)-ibe(j-1))) then
                write (LUPRI,*) 'mismatch in block lengths - msub1'
                call quit('mismatch in block lengths - msub1')
              end if
              elnew=ddot(lngth,gen(istrt22+ije2(i-1)+ibe(j-1)+1),1,     &
                         gen(istrt32+ije3(k-1)+1),1)
              idnew=ipack4(ia(1),ia(2),ia(3),ia(4))
              indnew=ipack4(ia(5),ia(6),ia(7),ia(8))
              iquit=0
              do while (iquit.eq.0)
                if (idold.eq.idnew) then
                  iquit=1   
                else if (idold.lt.idnew) then
                  idc=idc+1
                  idold=ij1(idc)
                  indc=0
                  indold=0
                else
                  write (6,*) 'msub1 - label not in template'
                  call quit('msub1 - label not in template')
                end if
              end do
              iquit=0
              do while (iquit.eq.0)
                if (indold.eq.indnew) then
                  iquit=1   
                else if (indold.lt.indnew) then
                  indc=indc+1
                  indold=igen(istrt11+ije1(idc-1)+indc)
                else
                  write (6,*) 'msub1 - label not in template 2'
                  call quit('msub1 - label not in template 2')
                end if
              end do
              gen(istrt12+ije1(idc-1)+indc)=elnew
            end if
60        continue
50      continue
10    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine msub2 (name2,name3,namet)
!--------------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use ctmp_mod
      use dpd_mod
      use hsfs_parameters
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
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#include "ctmp.h"
!     common /dpd/ nir,igmt(8,8),ici(200)
#include "dpd.h"
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension ia(8),ib(10000),ibe(0:10000)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub2 input file does not exist'
        call quit('msub2 input file does not exist')
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub2 input file does not exist'
        call quit('msub2 input file does not exist')
      end if
      if ((iii2.lt.22222222).or.(iii2.gt.33333333)) then
        write (6,*) 'msub2 called with inconsistent iii2'
        call quit('msub2 called with inconsistent iii2')
      end if
      if (iii3.gt.3333) then
        write (6,*) 'msub2 called with inconsistent iii3'
        call quit('msub2 called with inconsistent iii3')
      end if
      if (iii3.gt.33) then
        if (((iii2-10000*(iii2/10000)).ne.iii3).or.           &
            ((iiv2-10000*(iiv2/10000)).ne.iiv3)) then
          write (6,*) 'contracted indices do not match - msub2'
          call quit('contracted indices do not match - msub2')
        end if
        iii1=iii2/10000 
        iiv1=iiv2/10000 
        call tmplt (iii1,iiv1,namet,ntapet)
        istrt11=nstrt(1,ntapet-1)+2*nij1+4
        istrt12=nstrt(2,ntapet-1)
        istrt22=nstrt(2,ntape2-1)
        istrt32=nstrt(2,ntape3-1)
        call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1)
        id=0
        idc=0
        do 90 i=1,ije1(nij1)
          ind2=igen(istrt11+i)
          iquit=0
          do while (iquit.eq.0)
            if (id.eq.ind2) then
              iquit=1
            else if (id.lt.ind2) then
              idc=idc+1
              id=ij2(idc)
            else
              write (6,*) 'msub2 - id not in template'
              call quit('msub2 - id not in template')
            end if
          end do
          if (ije3(nij3).ne.(ije2(idc)-ije2(idc-1))) then
            write (LUPRI,*) 'mismatch in block lengths - msub2 '
            call quit('mismatch in block lengths - msub2')
          end if
          gen(istrt12+i)=ddot(ije3(nij3),gen(istrt22+ije2(idc-1)+1),1,  &
                                         gen(istrt32+1),1)
90      continue
      else if (iii3.le.33) then
        if (((iii2-100*(iii2/100)).ne.iii3).or.                         &
            ((iiv2-100*(iiv2/100)).ne.iiv3)) then
          write (6,*) 'contracted indices do not match - msub2'
          call quit('contracted indices do not match - msub2')
        end if
        iii1=iii2/100 
        iiv1=iiv2/100 
        call tmplt6 (iii1,iiv1,namet,ntapet)
        istrt11=nstrt(1,ntapet-1)+2*nij1+4
        istrt12=nstrt(2,ntapet-1)
        istrt22=nstrt(2,ntape2-1)
        istrt32=nstrt(2,ntape3-1)
        call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1)
        call tchng8 (nij2,ij2,ije2,ntape2)
        id4=itmp(1)
        ind4=itmp(ndp3+1)
        imark=1
        do 200 i=1,nij1
          id3=ij1(i)
          do 210 j=1,ije1(i)-ije1(i-1)
            ind3=iand(igen(istrt11+ije1(i-1)+j),z'ffffff')
            iquit=0
            do while (iquit.eq.0)
              ii=imax(id3,ind3,id4,ind4,2)
              if (ii.eq.0) then
                iquit=1
              else if (ii.eq.1) then
                imark=imark+1   
                id4=itmp(imark)
                ind4=itmp(ndp3+imark)
              else if (ii.eq.2) then
                write (6,*) 'msub2 - label not in template'
                call quit('msub2 - label not in template')
              end if
            end do
            lngth=ije2(imark)-ije2(imark-1)
            if (lngth.ne.ije3(nij3)) then
              write (LUPRI,*) 'mismatch in block lengths - msub2'
              call quit('mismatch in block lengths - msub2')
            end if
            gen(istrt12+ije1(i-1)+j)=ddot(lngth,   &
                    gen(istrt22+ije2(imark-1)+1),1,gen(istrt32+1),1)
210       continue
200     continue
      end if
     
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine msub3(name2,name3,namet)
!--------------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------------
      use ipack_mod, only: ipack3,ipack4
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use ctmp_mod
      use dpd_mod
      use hsfs_parameters
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
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#include "ctmp.h"
!     common /dpd/ nir,igmt(8,8),ici(200)
#include "dpd.h"
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension ia(8),ib(10000),ibe(0:10000)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)
!     ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub3 input file does not exist, name2=',name2
        call quit('msub3 input file name2 does not exist!')
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub3 input file does not exist, name3=',name3
        call quit('msub3 input file name3 does not exist !')
      end if
      if (((iii2-1000*(iii2/1000)).ne.(iii3-1000*(iii3/1000))).or.      &
         ((iiv2-1000*(iiv2/1000)).ne.(iiv3-1000*(iiv3/1000)))) then
        write (6,*) 'contracted indices do not match - msub3'
        call quit('contracted indices do not match - msub3')
      end if
      if ((iii2.lt.22222222).or.(iii2.gt.33333333)) then
        write (6,*) 'msub3 called with inconsistent iii2,iii2=',iii2
        call quit('msub3 called with inconsistent iii2!')
      end if
      if ((iii3.lt.2222).or.(iii3.gt.3333)) then
        write (6,*) 'msub3 called with inconsistent iii3, iii3=',iii3
        call quit('msub3 called with inconsistent iii3!')
      end if
      iii1=10*(iii2/1000)+iii3/1000
      iiv1=10*(iiv2/1000)+iiv3/1000
      call tmplt6 (iii1,iiv1,namet,ntapet)
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)
      call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1)
      if ((ije2(nij2).eq.0).or.(ije3(nij3).eq.0)) then
        return
      end if
      call tmpcng (nij3,ij3,ije3,ntape3,2,1)
      idc=0
      idold=0 
      ibe(0)=0
      do 110 i=1,nij2
        do 120 j=0,3
          ia(4-j)=iand(ishft(ij2(i),-8*j),z'ff')
120     continue
        ii=0
        n=0
        do 130 j=1,ije2(i)-ije2(i-1)
          jj=iand(ishft(igen(istrt21+ije2(i-1)+j),-24),z'ff')
          if (ii.ne.jj) then
            n=n+1
            ib(n)=jj
            ii=jj
            ibe(n)=ibe(n-1)+1
          else
            ibe(n)=ibe(n)+1
         end if
130     continue
        call dimx (n,10000,'msub3')
        do 150 j=1,n
          ia(5)=ib(j)
          do 160 k=1,nij3 
            ia(6)=ij3(k)
            ii=1
            do 180 m=1,6
              ii=igmt(ii,ici(ia(m)))
180         continue
            if (ii.eq.1) then
              lngth=ije3(k)-ije3(k-1)
              if (lngth.ne.(ibe(j)-ibe(j-1))) then
                write (LUPRI,*) 'mismatch in block lengths - msub3'
                call quit('mismatch in block lengths - msub3')
              end if
              elnew=ddot(lngth,gen(istrt22+ije2(i-1)+ibe(j-1)+1),1,     &
                         gen(istrt32+ije3(k-1)+1),1)
              idnew=ipack3(ia(1),ia(2),ia(3))
              indnew=ipack3(ia(4),ia(5),ia(6))
              iquit=0
              do while (iquit.eq.0)
                if (idold.eq.idnew) then
                  iquit=1   
                else if (idold.lt.idnew) then
                  idc=idc+1
                  idold=ij1(idc)
                  indc=0
                  indold=0
                else
                  write (6,*) 'msub3 - label not in template'
                  stop
                end if
              end do
              iquit=0
              do while (iquit.eq.0)
                if (indold.eq.indnew) then
                  iquit=1   
                else if (indold.lt.indnew) then
                  indc=indc+1
                  ii=igen(istrt11+ije1(idc-1)+indc)
                  indold=iand(ii,z'ffffff')  
                else
                  write (6,*) 'msub3 - label not in template'
                  stop
                end if
              end do
              gen(istrt12+ije1(idc-1)+indc)=elnew
            end if
160       continue
150     continue
110   continue

      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine msub4 (name2,name3,namet)
!--------------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------------
      use ipack_mod, only: ipack4
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use dpd_mod
      use hsfs_parameters
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
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
!     common /dpd/ nir,igmt(8,8),ici(200)
#include "dpd.h"
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension ia(8),ib(10000),ibe(0:10000)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub4 input file does not exist, name2=',name2
        call quit('msub4 input file name2 does not exist!')
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub4 input file does not exist,name3=',name3
        call quit('msub4 input file name3 does not exist!')
      end if
      if (((iii2-100*(iii2/100)).ne.(iii3-100*(iii3/100))).or.      &
         ((iiv2-100*(iiv2/100)).ne.(iiv3-100*(iiv3/100)))) then
         write (6,*) 'contracted indices do not match - msub4'
         call quit('contracted indices do not match - msub4!')
      end if
      if ((iii2.lt.22222222).or.(iii2.gt.33333333)) then
        write (6,*) 'msub4 called with inconsistent iii2=',iii2
        call quit('msub4 called with inconsistent iii2')
      end if
      if ((iii3.lt.2222).or.(iii3.gt.3333)) then
        write (6,*) 'msub4 called with inconsistent iii3=',iii3
        call quit('msub4 called with inconsistent iii3!')
      end if
      iii1=100*(iii2/100)+(iii3/100)
      iiv1=100*(iiv2/100)+(iiv3/100)
      call tmplt8 (iii1,iiv1,namet,ntapet)
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)
      call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1)
      if ((ije2(nij2).eq.0).or.(ije3(nij3).eq.0)) then
        return
      end if
      idc=0
      indc=0
      idold=0 
      indold=0
      ibe(0)=0
      do 300 i=1,nij2
        do 310 j=0,3
          ia(4-j)=iand(ishft(ij2(i),-8*j),z'ff')
310     continue
        ii=0
        n=0
        do 320 j=1,ije2(i)-ije2(i-1)
          jj=iand(ishft(igen(istrt21+ije2(i-1)+j),-16),z'ffff')
          if (ii.ne.jj) then
            n=n+1
            ib(n)=jj
            ii=jj
            ibe(n)=ibe(n-1)+1
          else
            ibe(n)=ibe(n)+1
          end if
320     continue
        call dimx (n,10000,'msub4')
        do 330 j=1,n
          ia(5)=iand(ishft(ib(j),-8),z'ff')
          ia(6)=iand(ib(j),z'ff')
          do 340 k=1,nij3 
            ia(7)=iand(ishft(ij3(k),-8),z'ff')
            ia(8)=iand(ij3(k),z'ff')
            ii=1
            do 350 m=1,8
              ii=igmt(ii,ici(ia(m)))
350         continue
            if (ii.eq.1) then
              lngth=ije3(k)-ije3(k-1)
              if (lngth.ne.(ibe(j)-ibe(j-1))) then
                write (LUPRI,*) 'mismatch in block lengths - msub4'
                call quit('mismatch in block lengths - msub4')
              end if
              elnew=ddot(lngth,gen(istrt22+ije2(i-1)+ibe(j-1)+1),1,     &
                         gen(istrt32+ije3(k-1)+1),1)
              idnew=ipack4(ia(1),ia(2),ia(3),ia(4))
              indnew=ipack4(ia(5),ia(6),ia(7),ia(8))
              iquit=0
              do while (iquit.eq.0)
                ii=imax(idnew,indnew,idold,indold,2)
                if (ii.eq.0) then
                  iquit=1
                else if (ii.eq.1) then
                  indc=indc+1
                  indold=igen(istrt11+indc)
                  if (indc.gt.ije1(idc)) then 
                    idc=idc+1
                    idold=ij1(idc)
                  end if
                else if (ii.eq.2) then
                  write (6,*) 'label not in template - msub4'
                  stop
                end if
              end do
              gen(istrt12+indc)=elnew
            end if
340       continue
330     continue
300   continue

      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine msub5 (name2,name3,namet)
!--------------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------------
      use ipack_mod, only: ipack4
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use ctmp_mod
      use dpd_mod
      use hsfs_parameters
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
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#include "ctmp.h"
!     common /dpd/ nir,igmt(8,8),ici(200)
#include "dpd.h"
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension ia(8),ib(10000),ibe(0:10000)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

      !call check_modules

      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub5 input file does not exist, name2=',name2
        call quit('msub5 input file name2 does not exist!')
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub5 input file does not exist,name3=',name3
        call quit('msub5 input file name3 does not exist!')
      end if
      if (((iii2-100*(iii2/100)).ne.(iii3-100*(iii3/100))).or.      &
         ((iiv2-100*(iiv2/100)).ne.(iiv3-100*(iiv3/100)))) then
        write (6,*) 'contracted indices do not match - msub5'
        call quit('contracted indices do not match - msub5!')
      end if
      if ((iii2.lt.222222).or.(iii2.gt.333333)) then
        write (6,*) 'msub5 called with inconsistent iii2=',iii2
        call quit('msub5 called with inconsistent iii2!')
      end if
      if ((iii3.lt.222222).or.(iii3.gt.333333)) then
        write (6,*) 'msub5 called with inconsistent iii3=',iii3
        call quit('msub5 called with inconsistent iii3!')
      end if
      iii1=10000*(iii2/100)+iii3/100
      iiv1=10000*(iiv2/100)+iiv3/100
      call tmplt8 (iii1,iiv1,namet,ntapet)
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)
      call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1)
      if ((ije2(nij2).eq.0).or.(ije3(nij3).eq.0)) then
        return
      end if
      call tchng6 (nij2,ij2,ije2,ntape2,4)
      call tchng6 (nij3,ij3,ije3,ntape3,4)
      idc=0
      indc=0
      idold=0 
      indold=0
      do 410 i=1,nij2
        do 420 k=0,3
          ia(4-k)=iand(ishft(ij2(i),-8*k),z'ff')
420     continue
        do 430 j=1,nij3
          do 440 k=0,3
            ia(8-k)=iand(ishft(ij3(j),-8*k),z'ff')
440       continue
          ii=1
          do 450 k=1,8
            ii=igmt(ii,ici(ia(k)))
450       continue
          if (ii.eq.1) then
            lngth=ije2(i)-ije2(i-1)
            if (lngth.ne.(ije3(j)-ije3(j-1))) then
              write (LUPRI,*) 'mismatch in block lengths - msub5'
              call quit('mismatch in block lengths - msub5')
            end if
            elnew=ddot(lngth,gen(istrt22+ije2(i-1)+1),1,       &
                       gen(istrt32+ije3(j-1)+1),1)
            idnew=ipack4(ia(1),ia(2),ia(3),ia(4))
            indnew=ipack4(ia(5),ia(6),ia(7),ia(8))
            iquit=0
            do while (iquit.eq.0)
              ii=imax(idnew,indnew,idold,indold,2)
              if (ii.eq.0) then
                iquit=1
              else if (ii.eq.1) then
                indc=indc+1
                indold=igen(istrt11+indc)
                if (indc.gt.ije1(idc)) then 
                  idc=idc+1
                  idold=ij1(idc)
                end if
              else if (ii.eq.2) then
                write (6,*) 'label not in template - msub5'
                stop
              end if
            end do
            gen(istrt12+indc)=elnew
          end if
430     continue
410   continue

      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine msub6 (name2,name3,namet)
!--------------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------------
      use ipack_mod, only: ipack4
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use ctmp_mod
      use dpd_mod
      use hsfs_parameters
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
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#include "ctmp.h"
!     common /dpd/ nir,igmt(8,8),ici(200)
#include "dpd.h"
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension ia(8),ib(10000),ibe(0:10000)
!     ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

      !call check_modules

      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub6 input file does not exist, name2=',name2
        call quit('msub6 input file name2 does not exist!')
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub6 input file does not exist, name3=',name3
        call quit('msub6 input file name3 does not exist!')
      end if
      if (((iii2-10*(iii2/10)).ne.(iii3-10*(iii3/10))).or.       &
          ((iiv2-10*(iiv2/10)).ne.(iiv3-10*(iiv3/10)))) then
        write (6,*) 'contracted indices do not match - msub6'
        call quit('contracted indices do not match - msub6!')
      end if
      if ((iii2.lt.22222222).or.(iii2.gt.33333333)) then
        write (6,*) 'msub6 called with inconsistent iii2=',iii2
        call quit('msub6 called with inconsistent iii2!')
      end if
      if (iii3.gt.33) then
        write (6,*) 'msub6 called with inconsistent iii3=',iii3
        call quit('msub6 called with inconsistent iii3!')
      end if
      iii1=10*(iii2/10)+iii3/10
      iiv1=10*(iiv2/10)+iiv3/10
      call tmplt8 (iii1,iiv1,namet,ntapet)
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)
      call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1)
      if ((ije2(nij2).eq.0).or.(ije3(nij3).eq.0)) then
        return
      end if
      call tmpcng (nij3,ij3,ije3,ntape3,2,3)
      idc=0
      indc=0
      idold=0 
      indold=0
      ibe(0)=0
      do 110 i=1,nij2
        do 120 j=0,3
          ia(4-j)=iand(ishft(ij2(i),-8*j),z'ff')
120     continue
        ii=0
        n=0
        do 130 j=ije2(i-1)+1,ije2(i)
          jj=iand(ishft(igen(istrt21+j),-8),z'ffffff')
          if (ii.ne.jj) then
            n=n+1
            ib(n)=jj
            ii=jj
            ibe(n)=ibe(n-1)+1
          else
            ibe(n)=ibe(n)+1
          end if
130     continue
        call dimx (n,10000,'msub6')
        do 150 j=1,n
          ia(5)=iand(ishft(ib(j),-16),z'ff')
          ia(6)=iand(ishft(ib(j),-8),z'ff')
          ia(7)=iand(ib(j),z'ff')
          do 160 k=1,nij3 
            ia(8)=ij3(k)
            ii=1
            do 180 m=1,8
              ii=igmt(ii,ici(ia(m)))
180         continue
            if (ii.eq.1) then
              lngth=ije3(k)-ije3(k-1)
              if (lngth.ne.(ibe(j)-ibe(j-1))) then
                write (LUPRI,*) 'mismatch in block lengths - msub6'
                call quit('mismatch in block lengths - msub6')
              end if
              elnew=ddot(lngth,gen(istrt22+ije2(i-1)+ibe(j-1)+1),1,     &
                         gen(istrt32+ije3(k-1)+1),1)
              idnew=ipack4(ia(1),ia(2),ia(3),ia(4))
              indnew=ipack4(ia(5),ia(6),ia(7),ia(8))
              iquit=0
              do while (iquit.eq.0)
                ii=imax(idnew,indnew,idold,indold,2)
                if (ii.eq.0) then
                  iquit=1
                else if (ii.eq.1) then
                  indc=indc+1
                  indold=igen(istrt11+indc)
                  if (indc.gt.ije1(idc)) then 
                    idc=idc+1
                    idold=ij1(idc)
                  end if
                else if (ii.eq.2) then
                  write (6,*) 'label not in template - msub6'
                  stop
                end if
              end do
              gen(istrt12+indc)=elnew
            end if
160       continue
150     continue
110   continue

      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine mult(name2,name3,namet,mmm)
!-----------------------------------------------------------------------
!
!   multiplication of 2 diagrams: name2 x name3 -> namet
!
! where "mmm" is type of the multiplication, mmm cosists of 3 digits m1m2m3
!
! m1- how many free lines in name2 (are not participating in multiplication)
! m2 -                       name3
! m3 - how many lines will be contracted (multiplied) to produce final diagram "namet"
!       
!-----------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use hsfs_parameters
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
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
#endif

      if (getiprint().ge.5) then
       call header ('Output from mult',-1)
       write(LUPRI,'(2x,a,a4)') 'name2=',name2
       write(LUPRI,'(2x,a,a4)') 'name3=',name3
       write(LUPRI,'(2x,a,a4)') 'namet=',namet
       write(LUPRI,'(2x,a,i4)') 'mmm=',mmm
       call print_cijx_stuff
      endif

      ! split mmm into 3 digits m1,m2,m3
      call parse(mmm,idum,m1,m2,m3)

      if (getiprint().ge.5) then
        write(LUPRI,*) 'mult: obtained m1,m2,m3:',m1,m2,m3
      endif

      ! find diagram "name2", fill cij2 common block,returns number ntape2
      call files_ccsd(name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (LUPRI,'(2x,a,a4)') ' mult input tape does not exist, ',name2
        call abext('mult    ')
      endif

      ! find diagram "name3", fill cij3 common block
      call files_ccsd(name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) ' mult input tape does not exist, ',name3
        call abext('mult    ')
      endif

      ! check if 2 multiplying files fit each other
      n=10**m3
      if (((iii2-n*(iii2/n)).ne.(iii3-n*(iii3/n))).or.((iiv2-n*(iiv2/n)).ne.(iiv3-n*(iiv3/n)))) then
        write (6,*) 'Contracted indices do not match in mult'
        write (6,*) 'iii2, iii3, iiv2, iiv3 = '
        write (6,*) iii2,iii3,iiv2,iiv3
        call abext('mult    ')
      end if

!    make characteristic of the final diagram "namet"
      iii=10**m2*(iii2/10**m3)+iii3/10**m3
      iiv=10**m2*(iiv2/10**m3)+iiv3/10**m3
!    create structure of the final diagram
      call tmplt(iii,iiv,namet,ntapet)
      istrt11=nstrt(1,ntapet-1)+2*nij1+4 
      istrt12=nstrt(2,ntapet-1)
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)

      nsize1=ije1(nij1) ! how many values in resulting diagram
      nsize2=ije2(nij2) !                in 1 diagram
      nsize3=ije3(nij3) !                in 2 diagram

      if (getiprint().ge.10) then
        write(LUPRI,*) 'mult: istrt11,istrt12:',istrt11,istrt12
        write(LUPRI,*) 'mult: istrt21,istrt22:',istrt21,istrt22
        write(LUPRI,*) 'mult: istrt31,istrt32:',istrt31,istrt32
        write(LUPRI,*) 'mult: nij1,nij2,nij3:',nij1,nij2,nij3
        write(LUPRI,*) 'mult: nsize1,nsize2,nsize3:',nsize1,nsize2,nsize3
      endif
!
      if (m3.ne.2) then
        call tmpcng(nij1,ij1,ije1,ntapet,2,4-m2)
        call tmpcng(nij2,ij2,ije2,ntape2,2,4-m3)
        call tmpcng(nij3,ij3,ije3,ntape3,2,4-m3)
      endif

!       write(LUPRI,*) 'mult: nsize1,nsize2,nsize3:',
!    &                  nsize1,nsize2,nsize3
!
      ! zeroe array for final diagram values
      call dcopy(nsize1,0.d0,0,gen(istrt12+1),1)

!       write(LUPRI,*) 'mult: nsize1,nsize2,nsize3:',
!    &                  nsize1,nsize2,nsize3

      kbb=0
      kcc=0
      ib0=z'ffff'
      msk1=256**m2-1

!===============================================
!       cycle of i1 for the final diagram
!===============================================
      do 250 i1=1,nsize1

      ib=ishft(igen(istrt11+i1),-8*m2) ! get 2 indexes at once, ib
      ic=iand(igen(istrt11+i1),msk1)   ! get 2 indexes at once, ic

      if (ib.eq.ib0) go to 225
      do 210 kb=kbb+1,nij2
      if (ib.eq.ij2(kb)) then
        ib0=ib
        kb0=kb
        go to 220
      else if (ib.lt.ij2(kb)) then
        go to 250
      endif
210   continue
      go to 250
220   if (kb0.ne.kbb) then
        kcc=0
        kbb=kb0
      endif
225   do 230 kc=kcc+1,nij3
      if (ic.eq.ij3(kc)) then
        kc0=kc
        go to 240
      else if (ic.lt.ij3(kc)) then
        go to 250
      endif
230   continue
      go to 250
240   kcc=kc0
      n2=ije2(kbb)-ije2(kbb-1)
      if (ije3(kcc)-ije3(kcc-1).ne.n2) then
        write (*,2100) ij2(kbb),ij3(kcc),name2,name3
2100    format('mismatch in length of records',2(1x,z6),' of files ',2a8,'. check isym.')
        write (*,2200) n2,(igen(istrt21+i),i=ije2(kbb-1)+1,ije2(kbb))
        n3=ije3(kcc)-ije3(kcc-1)
        write (*,2200) n3,(igen(istrt31+i),i=ije3(kcc-1)+1,ije3(kcc))
2200    format(' length and indices',i5/(10(1x,z8)))
        call abext('mult    ')
      endif

      n20=istrt22+ije2(kbb-1)+1
      n30=istrt32+ije3(kcc-1)+1

      gen(istrt12+i1)=ddot(n2,gen(n20),1,gen(n30),1)

      if (getiprint().ge.13) then
       write(LUPRI,*) 'mult: istrt12+i1,n20,n30,n2:',                   &
                 istrt12+i1,n20,n30,n2
       write(LUPRI,*) 'mult: gen(n20),gen(n30),gen(istrt12+i1):',       &
             gen(n20),gen(n30),gen(istrt12+i1)
      endif

! insert printing & checking
      if (dabs(gen(istrt12+i1)).ge.1.0d+5) then
        write(LUPRI,*) 'mult: too big value !, gen=',gen(istrt12+i1) 
         do jj=1,n2
          write(LUPRI,'(2(i10,1x,d15.7,1x))')                           &
      n20+jj-1,gen(n20+jj-1),n30+jj-1,gen(n30+jj-1)
         enddo
        !call quit('mult: error exit!') 
      endif

250   continue

      return
      end subroutine mult

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine mult6 (name2,name3,namet,id)
!--------------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------------
      use ipack_mod, only: ipack2,ipack3,ipack4
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use dpd_mod
      use hsfs_parameters
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
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
!     common /dpd/ nir,igmt(8,8),ici(200)
#include "dpd.h"
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension ia(8),itc(10000),itce(0:10000)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ipack2(i,j)=ior(ishft(i,8),j)
!   ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)
!    ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'mult6 input tape does not exist, name2=',name2
        call quit('mult6 input tape name2 does not exist!')
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'mult6 input tape does not exist, name3=',name3
        call quit('mult6 input tape name3 does not exist!')
      end if
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)
      nsize2=ije2(nij2)
      nsize3=ije3(nij3)
      if ((id.ne.0).and.(id.ne.5)) then
        if (id.eq.1) then
          if (((iii2-10*(iii2/10)).ne.(iii3-10*(iii3/10))).or.          &
             ((iiv2-10*(iiv2/10)).ne.(iiv3-10*(iiv3/10)))) then
            write (6,*) 'contracted indices do not match - mult6'
            write (6,*) 'iii2=',iii2,' iii3=',iii3
            write (6,*) 'iiv2=',iiv2,' iiv3=',iiv3
            call quit('contracted indices do not match - mult6')
          end if
          call tmpcng (nij2,ij2,ije2,ntape2,2,3)
          call tmpcng (nij3,ij3,ije3,ntape3,2,3)
          iii1=1000*(iii2/10)+iii3/10
          iiv1=1000*(iiv2/10)+iiv3/10
          call tmplt6 (iii1,iiv1,namet,ntapet)
          istrt11=nstrt(1,ntapet-1)+2*nij1+4
          istrt12=nstrt(2,ntapet-1)
          msk1=z'ffffff'
          msk2=z'ffffff'
        else if (id.eq.2) then
          if (((iii2-1000*(iii2/1000)).ne.(iii3-1000*(iii3/1000))).or.  &
              ((iiv2-1000*(iiv2/1000)).ne.(iiv3-1000*(iiv3/1000)))) then
            write (6,*) 'contracted indices do not match - mult6'
            write (6,*) 'iii2=',iii2,' iii3=',iii3
            write (6,*) 'iiv2=',iiv2,' iiv3=',iiv3
            call quit('contracted indices do not match - mult6')
          end if
          if (iii2.gt.3333) then
            call tmpcng (nij3,ij3,ije3,ntape3,2,1)
            iii1=10*(iii2/1000)+iii3/1000
            iiv1=10*(iiv2/1000)+iiv3/1000
            msk1=z'ffffff'
            msk2=z'ff'
            call tmplt (iii1,iiv1,namet,ntapet)
            istrt11=nstrt(1,ntapet-1)+2*nij1+4
            istrt12=nstrt(2,ntapet-1)
            call tmpcng (nij1,ij1,ije1,ntapet,2,3)
          else if (iii2.le.3333) then
            call tmpcng (nij2,ij2,ije2,ntape2,2,1)
            iii1=1000*(iii2/1000)+iii3/1000
            iiv1=1000*(iiv2/1000)+iiv3/1000
            msk1=z'ff'
            msk2=z'ffffff'
            call tmplt (iii1,iiv1,namet,ntapet)
            istrt11=nstrt(1,ntapet-1)+2*nij1+4
            istrt12=nstrt(2,ntapet-1)
            call tmpcng (nij1,ij1,ije1,ntapet,2,1)
          end if
        else if (id.eq.3) then
          if (((iii2-100*(iii2/100)).ne.(iii3-100*(iii3/100))).or.      &
             ((iiv2-100*(iiv2/100)).ne.(iiv3-100*(iiv3/100)))) then
            write (6,*) 'contracted indices do not match - mult6'
            write (6,*) 'iii2=',iii2,' iii3=',iii3
            write (6,*) 'iiv2=',iiv2,' iiv3=',iiv3
            call quit('contracted indices do not match - mult6, id=3')  
          end if
          if (iii2.gt.3333) then
            call tchng6 (nij2,ij2,ije2,ntape2,4)
            iii1=100*(iii2/100)+iii3/100
            iiv1=100*(iiv2/100)+iiv3/100
            call tmplt6 (iii1,iiv1,namet,ntapet)
            istrt11=nstrt(1,ntapet-1)+2*nij1+4
            istrt12=nstrt(2,ntapet-1)
            call tchng6 (nij1,ij1,ije1,ntapet,4)
            msk1=z'ffffffff'
            msk2=z'ffff'
          else if (iii2.le.3333) then
            call tchng6 (nij3,ij3,ije3,ntape3,4)
            iii1=10000*(iii2/100)+iii3/100
            iiv1=10000*(iiv2/100)+iiv3/100
            call tmplt6 (iii1,iiv1,namet,ntapet)
            istrt11=nstrt(1,ntapet-1)+2*nij1+4
            istrt12=nstrt(2,ntapet-1)
            call tchng6 (nij1,ij1,ije1,ntapet,2)
            msk1=z'ffff'
            msk2=z'ffffffff'
          end if
        else if (id.eq.4) then
          if (((iii2-1000*(iii2/1000)).ne.(iii3-1000*(iii3/1000))).or.  &
         ((iiv2-1000*(iiv2/1000)).ne.(iiv3-1000*(iiv3/1000)))) then
            write (6,*) 'contracted indices do not match - mult6'
            write (6,*) 'iii2=',iii2,' iii3=',iii3
            write (6,*) 'iiv2=',iiv2,' iiv3=',iiv3
            call quit('contracted indices do not match - mult6,id=4')
          end if
          iii1=1000*(iii2/1000)+iii3/1000
          iiv1=1000*(iiv2/1000)+iiv3/1000
          call tmplt6 (iii1,iiv1,namet,ntapet)
          istrt11=nstrt(1,ntapet-1)+2*nij1+4
          istrt12=nstrt(2,ntapet-1)
          msk1=z'ffffff'
          msk2=z'ffffff'
        end if
        ibsave=msk1      
        icsave=msk2
        kbsave=0
        irst=0
        nsize1=ije1(nij1)
        call dcopy (nsize1,0.d0,0,gen(istrt12+1),1)
        do 10 i=1,nij1
         ib=ij1(i)
         do 20 j=ije1(i-1)+1,ije1(i)
          ic=iand(igen(istrt11+j),msk2)
          if ((ib.ge.ij2(1)).and.(ic.ge.ij3(1))) then
            if ((ib.le.ij2(nij2)).and.(ic.le.ij3(nij3))) then
              if (ib.ne.ibsave) then
                irst=1
                kb=kbsave+1
                idoitb=1
                do while (idoitb.eq.1) 
                  if (ib.eq.ij2(kb)) then
                    kbsave=kb
                    ibsave=ib
                    idoitb=0
                    igotb=1
                  else if (ib.lt.ij2(kb)) then
                    idoitb=0
                    igotb=0
                  end if
                  kb=kb+1
                end do   
              end if
              if (igotb.eq.1) then
                if (irst.eq.1) then
                  kcsave=0
                  irst=0
                  icsave=msk2
                end if
                if (ic.ne.icsave) then
                  kc=kcsave+1
                  idoitc=1
                  do while (idoitc.eq.1) 
                    if (ic.eq.ij3(kc)) then
                      kcsave=kc
                      icsave=ic
                      idoitc=0
                      igotc=1
                    else if (ic.lt.ij3(kc)) then
                      idoitc=0
                      igotc=0
                    end if
                    kc=kc+1
                  end do   
                end if
              end if
              if ((igotb.eq.1).and.(igotc.eq.1)) then
                n2=ije2(kbsave)-ije2(kbsave-1)
                n3=ije3(kcsave)-ije3(kcsave-1)
                if (n2.ne.n3) then
                  write (6,96) ij2(kbsave),ij3(kcsave),name2,name3
96                format ('mismatch in length of records',2(1x,z6),     &
                          ' of files ',2a8,'.check isym.')
                  write (6,95) n2,(igen(istrt21+im),                    &
             im=ije2(kbsave-1)+1,ije2(kbsave))
95                format (' length and indices',i5/(10(1x,z8)))
                  call quit('mismatch in length of records')    
                end if
                n20=istrt22+ije2(kbsave-1)+1
                n30=istrt32+ije3(kcsave-1)+1
                gen(istrt12+j)=ddot(n2,gen(n20),1,gen(n30),1)
              end if
            end if
          end if
20       continue
10      continue
      else if (id.eq.0) then
        if (iii3.gt.33) then
          if (((iii2-10000*(iii2/10000)).ne.iii3).or.                   &
              ((iiv2-10000*(iiv2/10000)).ne.iiv3)) then
            write (6,*) 'contracted indices do not match - mult6'
            write (6,*) 'iii2=',iii2,' iii3=',iii3
            write (6,*) 'iiv2=',iiv2,' iiv3=',iiv3
            call quit('contracted indices do not match - mult6,id=0')
          end if
          call tchng6 (nij2,ij2,ije2,ntape2,2)       
          iii1=iii2/10000
          iiv1=iiv2/10000
          call tmplt (iii1,iiv1,namet,ntapet)
          istrt11=nstrt(1,ntapet-1)+2*nij1+4
          istrt12=nstrt(2,ntapet-1)
          nsize1=ije1(nij1)
          call dcopy (nsize1,0.d0,0,gen(istrt12+1),1)
          i=1
          j=1
          do while (i.le.nsize1) 
            id1=igen(istrt11+i)
            id2=ij2(j)
            if (id1.eq.id2) then
              n2=ije2(j)-ije2(j-1)
              if (n2.ne.nsize3) then
                write (LUPRI,*) 'mismatch check isym'
                call quit('mismatch check isym')
              end if
              n20=istrt22+ije2(j-1)+1
              gen(istrt12+i)=ddot(n2,gen(n20),1,gen(istrt32+1),1)
              i=i+1
              j=j+1
            else if (id1.gt.id2) then
              j=j+1
            else if (id1.lt.id2) then
              write (LUPRI,*) 'label not in template'
              call quit('label not in template')                 
            end if
          end do   
        else if (iii3.le.33) then
          if (((iii2-100*(iii2/100)).ne.iii3).or.                       &
              ((iii2-100*(iii2/100)).ne.iii3)) then
            write (6,*) 'contracted indices do not match - mult6'
            write (6,*) 'iii2=',iii2,' iii3=',iii3
            write (6,*) 'iiv2=',iiv2,' iiv3=',iiv3
            call quit('contracted indices do not match - mult6,iii3')
          end if
          call tchng6 (nij2,ij2,ije2,ntape2,4)       
          iii1=iii2/100
          iiv1=iiv2/100
          call tmplt (iii1,iiv1,namet,ntapet)
          istrt11=nstrt(1,ntapet-1)+2*nij1+4
          istrt12=nstrt(2,ntapet-1)
          nsize1=ije1(nij1)
          call dcopy (nsize1,0.d0,0,gen(istrt12+1),1)
          i=1
          j=1
          do while (i.le.nsize1) 
            id1=igen(istrt11+i)
            id2=ij2(j)
            if (id1.eq.id2) then
              n2=ije2(j)-ije2(j-1)
              if (n2.ne.nsize3) then
                write (LUPRI,*) 'mismatch check isym - mult6'
                call quit('mismatch check isym - mult6')              
              end if
              n20=istrt22+ije2(j-1)+1
              gen(istrt12+i)=ddot(n2,gen(n20),1,gen(istrt32+1),1)
              i=i+1
              j=j+1
            else if (id1.gt.id2) then
              j=j+1
            else if (id1.lt.id2) then
              write (6,*) 'label not in template - mult6'
              stop                    
            end if
          end do   
        end if
      else if (id.eq.5) then
        if (((iii2-10*(iii2/10)).ne.(iii3-10*(iii3/10))).or.            &
           ((iiv2-10*(iiv2/10)).ne.(iiv3-10*(iiv3/10)))) then
          write (6,*) 'contracted indices do not match - mult6,id=5'
          call quit('contracted indices do not match - mult6')
        end if
        iii1=10*(iii2/10)+iii3/10
        iiv1=10*(iiv2/10)+iiv3/10
        call tmplt6 (iii1,iiv1,namet,ntapet)
        istrt11=nstrt(1,ntapet-1)+2*nij1+4
        istrt12=nstrt(2,ntapet-1)
        call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1)
        if ((ije2(nij2).ne.0).and.(ije3(nij3).ne.0)) then
          call tchng6 (nij2,ij2,ije2,ntape2,4)
          call tmpcng (nij3,ij3,ije3,ntape3,2,3)
          idc=0
          idold=0 
          itce(0)=0
          do 30 i=1,nij2
            do 40 j=0,3
              ia(4-j)=iand(ishft(ij2(i),-8*j),z'ff')
40          continue
            ii=0
            n=0
            do 50 j=1,ije2(i)-ije2(i-1)
              jj=iand(ishft(igen(istrt21+ije2(i-1)+j),-8),z'ff')
              if (ii.ne.jj) then
                n=n+1
                itc(n)=jj
                ii=jj
                itce(n)=itce(n-1)+1
              else
                itce(n)=itce(n)+1
              end if
50          continue
            call dimx (n,10000,'mult6')
            do 60 j=1,n
              ia(5)=itc(j)
              do 70 k=1,nij3 
                ia(6)=ij3(k) 
                ii=1
                do 80 m=1,6
                  ii=igmt(ii,ici(ia(m)))
80              continue
                if (ii.eq.1) then
                  lngth=ije3(k)-ije3(k-1)
                  if (lngth.ne.(itce(j)-itce(j-1))) then
                    write (LUPRI,*) 'mismatch in block lengths - mult6'
                    call quit('mismatch in block lengths - mult6')
                  end if
                  elnew=ddot(lngth,gen(istrt22+ije2(i-1)+itce(j-1)+1),1,&
                             gen(istrt32+ije3(k-1)+1),1)
                  idnew=ipack3(ia(1),ia(2),ia(3))
                  indnew=ipack4(ia(3),ia(4),ia(5),ia(6))
                  iquit=0
                  do while (iquit.eq.0)
                    if (idold.eq.idnew) then
                      iquit=1   
                    else if (idold.lt.idnew) then
                      idc=idc+1
                      idold=ij1(idc)
                      indc=0
                      indold=0
                    else
                      write (6,*) 'mult6 - label not in template'
                      call quit('mult6 - label not in template')
                    end if
                  end do
                  iquit=0
                  do while (iquit.eq.0)
                    if (indold.eq.indnew) then
                      iquit=1   
                    else if (indold.lt.indnew) then
                      indc=indc+1
                      indold=igen(istrt11+ije1(idc-1)+indc)
                    else
                      write (6,*) 'mult6 - label not in template'
                      stop
                    end if
                  end do
                  gen(istrt12+ije1(idc-1)+indc)=elnew
                end if
70            continue
60          continue
30        continue
        end if
      end if
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine mult8 (name2,name3,namet,iparm)
!--------------------------------------------------------------------------------------
!  iparm=1 : 6x4=8
!  iparm=2 : 8x4=4,8x2=6
!  iparm=3 : 8x4=6
!  iparm=4 : 8x4=8
!  iparm=5 : 6x6=8
!  iparm=6 : 8x2=8
!--------------------------------------------------------------------------------------
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      if (iparm.eq.1) then
        call msub1 (name2,name3,namet)
      else if (iparm.eq.2) then
        call msub2 (name2,name3,namet)
      else if (iparm.eq.3) then
        call msub3 (name2,name3,namet)
      else if (iparm.eq.4) then
        call msub4 (name2,name3,namet)
      else if (iparm.eq.5) then
        call msub5 (name2,name3,namet)
      else if (iparm.eq.6) then
        call msub6 (name2,name3,namet)
      else
        write (6,*) 'mult8 called with bad value of iparm '
        stop
      end if
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine abext(aa)
!--------------------------------------------------------
!  ... exit routine, extended by DIRAC quit
!--------------------------------------------------------
#include "priunit.h"
      character(*) aa
      write(LUPRI,*) 'abext: aborted from ',aa
      call flush(LUPRI)
      call quit(aa)
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine add(name1,f1,name2,f2,namet)
!-------------------------------------------------------------
!
!         what does it do ?
!
!-------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use hsfs_parameters

implicit none
#include "priunit.h"
integer :: name1,name2,namet
real*8 :: f1,f2
integer :: new, nsize,ntape1,ntape2,ntapet,i
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "cfiles.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij2.h"
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
#endif

      if (getiprint()>=5) then
       call header ('Output from add',-1)
       write(LUPRI,'(2x,a,a8)')    'name1=',name1
       write(LUPRI,'(2x,a,d12.6)') 'f1=',f1
       write(LUPRI,'(2x,a,a8)')    'name2=',name2
       write(LUPRI,'(2x,a,d12.6)') 'f2=',f2
       write(LUPRI,'(2x,a,a8)')    'namet=',namet
       call flush(LUPRI)
      endif

!   namet=name1*f1+name2*f2
      if ((name1.eq.name2).and.(name1.eq.namet)) then
        write (LUPRI,*) 'all three names in add are identical'
        call abext('add     ')
      end if
      call files_ccsd(name1,ntape1,new,2,-1,-1)
      if (new.ne.0) then
        write(LUPRI,'('' add input file does not exist '',a8)') name1
        call abext('add     ')
      endif
      call files_ccsd(name2,ntape2,new,3,-1,-1)
      if (new.ne.0) then
        write(LUPRI,'('' add input file does not exist '',a8)') name2
        call abext('add     ')
      endif

      ! ......   error trap
if ( (iii2.ne.iii3).or.(iiv2.ne.iiv3).or.(ije2(nij2).ne.ije3(nij3)) ) then
  write (LUPRI,'(/,2x,a)') 'add: added files do not match ! variables:'  
  write (LUPRI,*) 'iii2,iii3:',iii2,iii3
  write (LUPRI,*) 'iiv2,iiv3:',iiv2,iiv3
  write (LUPRI,*) 'nij2,nij3:',nij2,nij3 
  write (LUPRI,*) 'ije2(nij2),ije3(nij3):',ije2(nij2),ije3(nij3)
  write(LUPRI,*) 'name1,name2:',name1,name2
  write (LUPRI,123) name1,name2,iii2,iii3,iiv2,iiv3,ije2(nij2),ije3(nij3)
123     format(2a8,6i10)
        call flush(LUPRI)
        call abext('add     ')
endif

      call files_ccsd(namet,ntapet,new,1,iii2,iiv2)

      nsize=ije2(nij2) ! picks up nsize

      if (getiprint() >= 15) then
        write(LUPRI,*) 'add: nsize=ije2(nij2)=',nsize
       write(LUPRI,'(2x,3(a,a8))') 'name1=',name1,' name2=',name2,      &
                              ' namet=',namet
      endif

      if (name1.eq.namet) then
        if (f1.ne.1.d0) call rvmuls(nsize,gen(nstrt(2,ntape1-1)+1),     &
                                    f1,gen(nstrt(2,ntape1-1)+1))
        call daxpy(nsize,f2,gen(nstrt(2,ntape2-1)+1),                   &
                   1,gen(nstrt(2,ntape1-1)+1),1)
        if (getiprint()>=5) then
          write(LUPRI,*) 'add: name1.eq.namet, f1,f2=',f1,f2
          write(LUPRI,*) 'add: nstrt(2,ntape1-1)+1:',nstrt(2,ntape1-1)+1
          if (getiprint()>=18) then
            do i=1,nsize
               write(LUPRI,*) 'i=',i,'. gen(nstrt(2,ntape1-1)+1+i)=',   &
          gen(nstrt(2,ntape1-1)+1+i)
            enddo
          endif
        endif
      elseif (name2.eq.namet) then
        if (f2.ne.1.d0) call rvmuls(nsize,gen(nstrt(2,ntape2-1)+1),     &
                                    f2,gen(nstrt(2,ntape2-1)+1))
        call daxpy(nsize,f1,gen(nstrt(2,ntape1-1)+1),                   &
                   1,gen(nstrt(2,ntape2-1)+1),1)
        if (getiprint() >= 15) then
          write(LUPRI,*) 'add: name2.eq.namet, f1,f2=',f1,f2
          write(LUPRI,*) 'add: nstrt(2,ntape2-1)+1:',nstrt(2,ntape2-1)+1
          if (getiprint()>=18) then
            do i=1,nsize
               write(LUPRI,*) 'i=',i,' gen(nstrt(2,ntape2-1)+1+i)=',    &
          gen(nstrt(2,ntape1-1)+1+i)
            enddo
          endif
        endif
      else
        call rvmuls(nsize,gen(nstrt(2,ntape1-1)+1),                     &
                     f1,gen(nstrt(2,ntapet-1)+1))
        call daxpy(nsize,f2,gen(nstrt(2,ntape2-1)+1),                   &
                   1,gen(nstrt(2,ntapet-1)+1),1)
        if (getiprint() >= 15) then
          write(LUPRI,*) 'add: else branch, f1,f2=',f1,f2
          write(LUPRI,*) 'add: nstrt(2,ntapet-1)+1:',nstrt(2,ntapet-1)+1
          if (getiprint()>=18) then
            do i=1,nsize
              write(LUPRI,*) 'i=',i,'. gen(nstrt(2,ntapet-1)+1+i)=',    &
         gen(nstrt(2,ntapet-1)+1+i)
            enddo
          endif
        endif
      endif

      return
      end subroutine add

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine cmult(namef,f,namet)
!--------------------------------------------------------------------------------------
!
!   namet=f*namef
!
!--------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "cfiles.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij2.h"
#endif

      call files_ccsd(namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) ' cmult input file does not exist'
        call abext('cmult   ')
      endif
      call files_ccsd(namet,ntapet,new,0,iii2,iiv2)
      nsize=ije2(nij2)
      istrtf=nstrt(2,ntapef-1)
      istrtt=nstrt(2,ntapet-1)
      if (namef.ne.namet) then
        call imove(ntapef,ntapet)
      endif
      call rvmuls(nsize,gen(istrtf+1),f,gen(istrtt+1))
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dbdiag (name1,name2,name3,name4,irleop,isign)
!--------------------------------------------------------------------------------------
!
!
! what does it do ?
!
!
!
!--------------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use dat_mod
      use cij1_mod
      use crle_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
!#include "param_ccsd.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "cfiles.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!UK   common /crle/ trle(nrle,4),isize(4)
!     common /crle/ trle1(nrle1),trle2(nrle2),trle3(nrle3),
!    1 trle4(nrle4),isize(4)
#include "crle.h"
#endif

      dimension iopa(4)

      call parse (irleop,iopa(1),iopa(2),iopa(3),iopa(4))

      if (iopa(1).eq.1) then
        call files_ccsd (name1,ntape,new,1,-1,-1)
        isize(1)=ije1(nij1)
        call dimx(isize(1),nrle1,"nrle1   ")
        istrt1=nstrt(1,ntape-1)+2*nij1+4
        istrt2=nstrt(2,ntape-1)
        do 20 i=1,isize(1)
          it=igen(istrt1+i)
          ka=iand(ishft(it,-8),z'ff') ! unpack 2 indexes - k,a
          kb=iand(it,z'ff')
          trle1(i)=-1.d0*isign*(eps(ka)-eps(kb))
          trle1(isize(1)+i)=gen(istrt2+i)
20      continue
      end if

      if (iopa(2).eq.1) then
        call files_ccsd (name2,ntape,new,1,-1,-1)
        isize(2)=ije1(nij1)
        call dimx(isize(2),nrle2,"nrle2   ")
        istrt1=nstrt(1,ntape-1)+2*nij1+4
        istrt2=nstrt(2,ntape-1)
        do 10 i=1,isize(2)
          it=igen(istrt1+i)
          ka=iand(ishft(it,-24), z'ff')
          kb=iand(ishft(it,-16), z'ff')
          kc=iand(ishft(it,-8),  z'ff')
          kd=iand(it,            z'ff')
          trle2(i)=-1.d0*isign*(eps(ka)+eps(kb)-eps(kc)-eps(kd))
          trle2(isize(2)+i)=gen(istrt2+i)
10      continue
      end if

      if (iopa(3).eq.1) then
        call files_ccsd (name3,ntape,new,1,-1,-1)
        isize(3)=ije1(nij1)
        call dimx(isize(3),nrle3,"nrle3   ")
        istrt1=nstrt(1,ntape-1)+2*nij1+4
        istrt2=nstrt(2,ntape-1)
        do 30 i=1,nij1 
          id=ij1(i)             
          ka=iand(ishft(id,-16), z'ff')
          kb=iand(ishft(id,-8),  z'ff')
          kc=iand(id,z'ff')
          do 40 j=ije1(i-1)+1,ije1(i)
            ind=igen(istrt1+j)
            kd=iand(ishft(ind,-16), z'ff')
            ke=iand(ishft(ind,-8),  z'ff')
            kf=iand(ind,            z'ff')
            trle3(j)=-1.d0*isign*                                       &
       (eps(ka)+eps(kb)+eps(kc)-eps(kd)-eps(ke)-eps(kf))
            trle3(isize(3)+j)=gen(istrt2+j)
40        continue
30      continue
      end if

      if (iopa(4).eq.1) then
        call files_ccsd (name4,ntape,new,1,-1,-1)
        isize(4)=ije1(nij1)
        call dimx(isize(4),nrle4,"nrle4   ")
        istrt1=nstrt(1,ntape-1)+2*nij1+4
        istrt2=nstrt(2,ntape-1)
        if (iii1.le.333333) then
          do 50 i=1,nij1 
            id=ij1(i)             
            ka=iand(ishft(id,-16),  z'ff')
            kb=iand(ishft(id,-8),   z'ff')
            kc=iand(id,  z'ff')
            do 60 j=ije1(i-1)+1,ije1(i)
              ind=igen(istrt1+j)
              kd=iand(ishft(ind,-16),  z'ff')
              ke=iand(ishft(ind,-8),   z'ff')
              kf=iand(ind,             z'ff')
              trle4(j)=-1.d0*isign*                                     &
         (eps(ka)+eps(kb)+eps(kc)-eps(kd)-eps(ke)-eps(kf))
              trle4(isize(4)+j)=gen(istrt2+j)
60          continue
50        continue
        else if (iii1.gt.333333) then
          do 70 i=1,nij1 
            id=ij1(i)             
            ka=iand(ishft(id,-24),  z'ff')
            kb=iand(ishft(id,-16),  z'ff')
            kc=iand(ishft(id,-8),   z'ff')
            kd=iand(id,             z'ff')
            do 80 j=ije1(i-1)+1,ije1(i)
              ind=igen(istrt1+j)
              ke=iand(ishft(ind,-24),  z'ff')
              kf=iand(ishft(ind,-16),  z'ff')
              kg=iand(ishft(ind,-8),   z'ff')
              kh=iand(ind,             z'ff')
              trle4(j)=-1.d0*isign*                                     &
         (eps(ka)+eps(kb)+eps(kc)+eps(kd)                               &
         -eps(ke)-eps(kf)-eps(kg)-eps(kh))
              trle4(isize(4)+j)=gen(istrt2+j)
80          continue
70        continue
        end if
      end if
      return
      end subroutine dbdiag

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine detlst(nval)
!-----------------------------------------------------------------------------
!
!
! what does it do ?
!
!-----------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use dat_mod
      use det_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "dat.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "det.h"
!     common /det/ nvlspo,lstvso(100),lstdet(1000),ndets,isplit
#endif

      if (nval.eq.1) then
        ndets=0
        do 220 i=1,nvlspo
          idum=lstvso(i)                        
          if (ispin(idum,1).gt.0) then
            ndets=ndets+1
            lstdet(ndets)=idum
          end if
220     continue
      else if (nval.eq.2) then
        ndets=0
        do 230 i=1,nvlspo-1
          do 240 j=i+1,nvlspo
            idum=ior(ishft(lstvso(i),8),lstvso(j))
            if (ispin(idum,2).ge.0) then
              ndets=ndets+1
              lstdet(ndets)=idum
            end if
240       continue
230     continue
      else if (nval.eq.3) then
        ndets=0
        do 10 i=1,nvlspo-2
          do 20 j=i+1,nvlspo-1
            do 30 k=j+1,nvlspo
              idum=(ior(ior(ishft(lstvso(i),16),                        &
         ishft(lstvso(j),8)),lstvso(k)))
              if (ispin(idum,3).gt.0) then
                ndets=ndets+1
                lstdet(ndets)=idum
              end if
30          continue
20        continue
10      continue
      else if (nval.eq.4) then
        ndets=0
        do 40 i=1,nvlspo-3
          do 50 j=i+1,nvlspo-2
            do 60 k=j+1,nvlspo-1
              do 70 l=k+1,nvlspo
                idum=ior(                                               &
           ior(ishft(lstvso(i),24),ishft(lstvso(j),16)),                &
           ior(ishft(lstvso(k),8),lstvso(l)))
                if (ispin(idum,4).ge.0) then
                  ndets=ndets+1
                  lstdet(ndets)=idum
                end if
70            continue
60          continue
50        continue
40      continue
      else if (nval.eq.5) then
        ndets=0
        do 80 i=1,nvlspo-4
          do 90 j=i+1,nvlspo-3
            do 100 k=j+1,nvlspo-2
              do 110 l=k+1,nvlspo-1
                do 120 m=l+1,nvlspo
                  idum=ior(                                             &
             ior(ior(ishft(i,16),ishft(j,12)),                          &
             ior(ishft(k,8),ishft(l,4))),m)
                  if (ispin(idum,5).gt.0) then
                    ndets=ndets+1
                    lstdet(ndets)=idum
                  end if
120             continue
110           continue
100         continue
90        continue
80      continue
      else if (nval.eq.6) then
        ndets=0
        do 130 i=1,nvlspo-5
          do 140 j=i+1,nvlspo-4
            do 150 k=j+1,nvlspo-3
              do 160 l=k+1,nvlspo-2
                do 170 m=l+1,nvlspo-1
                  do 180 n=m+1,nvlspo
                    idum=ior(ior(ishft(i,20),ishft(j,16)),              &
               ior(ior(ishft(k,12),ishft(l,8)),ior(ishft(m,4),n)))
                    if (ispin(idum,6).ge.0) then
                      ndets=ndets+1
                      lstdet(ndets)=idum
                    end if
180               continue
170             continue
160           continue
150         continue
140       continue
130     continue
      end if
      if (nval.gt.6) then
        ndets=0
        do 190 i=1,isplit
          do 200 j=i+1,isplit
            do 210 k=isplit+1,nvlspo
              idum=ior(ishft(lstvso(i),16),                             &
         ior(ishft(lstvso(j),8),lstvso(k)))
              if (ispin(idum,nval).ge.0) then
                ndets=ndets+1
                lstdet(ndets)=idum
              end if
210         continue
200       continue
190     continue
      end if       
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dfmax3 (name1,name2,names,tmax,idmax,indmax)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add   (name1,-1.d0,name2,1.d0,names)
      call fmax6 (names,tmax,idmax,indmax)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dfmax4 (name1,name2,names,tmax,idmax,indmax)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add   (name1,-1.d0,name2,1.d0,names)
      call fmax8 (names,tmax,idmax,indmax)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dge5n2 (name2,name3,namet,iord)
!--------------------------------------------------------------------------------------
!c         6x6=10
!--------------------------------------------------------------------------------------
      use ipack_mod, only: ipack4
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use idv_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "cfiles.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij2.h"
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
!     common /idv/ iv(200) ,nvalo
#include "idv.h"
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension ia(8),save1(1000),save2(1000),iord(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'dge5n2 input tape does not exist, ',name2
        stop                      
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'dge5n2 input tape does not exist, ',name3
        stop                       
      end if
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)
      call tchng6 (nij2,ij2,ije2,ntape2,4)
      call tchng6 (nij3,ij3,ije3,ntape3,4)
      call tmplt10 (namet,ntapet)                    
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1)
      do 10 i=1,nij1
        do 60 k=4,0,-1
          ii=iand(ishft(ij1(i),-4*k),z'f')
          ia(5-k)=iv(ii)
60      continue
        do 20 j=1,ije1(i)-ije1(i-1)
          ii=igen(istrt11+ije1(i-1)+j)
          do 70 k=4,0,-1
            jj=iand(ishft(ii,-4*k),z'f')
            ia(10-k)=iv(jj)
70        continue
          igot2=0
          igot3=0
          ii=ipack4(ia(iord(1)),ia(iord(2)),ia(iord(3)),ia(iord(4)))
          do 30 k=1,nij2
            if (ij2(k).eq.ii) then
              k2save=k
              igot2=1
            end if
30        continue
          ii=ipack4(ia(iord(6)),ia(iord(7)),ia(iord(8)),ia(iord(9)))
          do 40 k=1,nij3
            if (ij3(k).eq.ii) then
              k3save=k
              igot3=1
            end if
40        continue
          if ((igot2.eq.1).and.(igot3.eq.1)) then
            nel1=0
            do 50 k=1,ije2(k2save)-ije2(k2save-1)
              ii=igen(istrt21+ije2(k2save-1)+k)
              if (ia(iord(5)).eq.iand(ishft(ii,-8),z'ff')) then
                nel1=nel1+1
                call dimx (nel1,1000,'dge5n2')
                save1(nel1)=gen(istrt22+ije2(k2save-1)+k)
              end if
50          continue 
            nel2=0
            do 80 k=1,ije3(k3save)-ije3(k3save-1)
              ii=igen(istrt31+ije3(k3save-1)+k)
              if (ia(iord(10)).eq.iand(ishft(ii,-8),z'ff')) then
                nel2=nel2+1
                call dimx (nel2,1000,'dge5n2')
                save2(nel2)=gen(istrt32+ije3(k3save-1)+k)
              end if
80          continue 
            if (nel1.ne.nel2) then  
              write (LUPRI,*) 'mismatch in lengths of records - dge5n2'
              call quit('mismatch in lengths of records - dge5n2')
            else  
              gen(istrt12+ije1(i-1)+j)=ddot(nel1,save1,1,save2,1)
            end if
          end if
20      continue
10    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dge6n2 (name2,name3,namet,iord)
!-----------------------------------------------------------------------------
!
!c      6x8=12
!
!-----------------------------------------------------------------------------
      use ipack_mod, only: ipack3,ipack4
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use idv_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "cfiles.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij2.h"
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
!     common /idv/ iv(200),nvalo
#include "idv.h"
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension save1(1000),save2(1000),ia(12),iord(12)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  todo: own module !
      !ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)
      !ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'dge6n2 input tape does not exist, name2=',name2
        call quit('dge6n2 input tape name2 does not exist')
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'dge6n2 input tape does not exist, name3=',name3
        call quit('dge6n2 input tape name3 does not exist')
      end if
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)
      call tmplt12 (namet,ntapet)
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1)
      call tchng6 (nij3,ij3,ije3,ntape3,4)
      do 10 i=1,nij1
        do 60 k=5,0,-1
          ii=iand(ishft(ij1(i),-4*k),z'f')
          ia(6-k)=iv(ii)
60      continue
        do 20 j=1,ije1(i)-ije1(i-1)
          ii=igen(istrt11+ije1(i-1)+j)
          do 70 k=5,0,-1
            jj=iand(ishft(ii,-4*k),z'f')
            ia(12-k)=iv(jj)
70        continue
          igot2=0
          igot3=0
          ii=ipack4(ia(iord(1)),ia(iord(2)),ia(iord(3)),ia(iord(4)))
          do 30 k=1,nij2
            if (ij2(k).eq.ii) then
              k2save=k
              igot2=1
            end if
30        continue
          ii=ipack4(ia(iord(8)),ia(iord(9)),ia(iord(10)),ia(iord(11)))
          do 50 k=1,nij3
            if (ij3(k).eq.ii) then
              k3save=k
              igot3=1
            end if
50        continue
          if ((igot2.eq.1).and.(igot3.eq.1)) then
            nel1=0
            ii=ipack3(ia(iord(5)),ia(iord(6)),ia(iord(7)))
            do 40 k=1,ije2(k2save)-ije2(k2save-1)
              jj=igen(istrt21+ije2(k2save-1)+k)
              if (ii.eq.iand(ishft(jj,-8),z'ffffff')) then
                nel1=nel1+1
                save1(nel1)=gen(istrt22+ije2(k2save-1)+k)
                call dimx (nel1,1000,'dge6n2')
              end if
40          continue
            nel2=0
            ii=ia(iord(12))  
            do 80 k=1,ije3(k3save)-ije3(k3save-1)
              jj=igen(istrt31+ije3(k3save-1)+k)
              if (ii.eq.iand(ishft(jj,-8),z'ff')) then
                nel2=nel2+1
                save2(nel2)=gen(istrt32+ije3(k3save-1)+k)
                call dimx (nel2,1000,'dge6n2')
              end if
80          continue
            if (nel1.ne.nel2) then
              write (LUPRI,*) 'length mismatch in dge6n2'
              call quit('length mismatch in dge6n2')
            else
              gen(istrt12+ije1(i-1)+j)=ddot(nel1,save1,1,save2,1)
            end if
          end if
20      continue
10    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine diage5 (name2,name3,name4,namet,iord)
!-------------------------------------------------------------------------------------
!
!                  6x(4+4)=10
!
!-------------------------------------------------------------------------------------
      use ipack_mod, only: ipack3,ipack4
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use cij4_mod
      use idv_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "cfiles.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij2.h"
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
!     common /cij4/ nij4,ij4(ndij),ije4(0:ndij),iii4,iiv4
#include "cij4.h"
!     common /idv/ iv(200),nvalo
#include "idv.h"
#endif

      dimension iord(10),ia(10)

!todo: into module !
!     ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)
!     ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage5 input tape does not exist, name2=',name2
        call quit("diage5: input tape name2 does not exist !")
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage5 input tape does not exist, ',name3
        call quit('diage5: input tape name3 does not exist !')
      end if
      call files_ccsd (name4,ntape4,new,4,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage5 input tape does not exist, ',name4
        call quit('diage5: input tape name4 does not exist')
      end if
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)
      istrt41=nstrt(1,ntape4-1)+2*nij4+4
      istrt42=nstrt(2,ntape4-1)

      call tmplt10 (namet,ntapet)

      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1)
      call tchng6 (nij2,ij2,ije2,ntape2,4)
      call tmpcng (nij3,ij3,ije3,ntape3,2,3)
      call tmpcng (nij4,ij4,ije4,ntape4,2,3)
      do 10 i=1,nij1
        do 60 k=4,0,-1
          ii=iand(ishft(ij1(i),-4*k),z'f')
          ia(5-k)=iv(ii)
60      continue
        do 20 j=1,ije1(i)-ije1(i-1)
          ii=igen(istrt11+ije1(i-1)+j)
          do 70 k=4,0,-1
            jj=iand(ishft(ii,-4*k),z'f')
            ia(10-k)=iv(jj)
70        continue
          igot2=0
          igot3=0
          igot4=0
          ii=ipack4(ia(iord(1)),ia(iord(2)),ia(iord(3)),ia(iord(4)))
          do 30 k=1,nij2
            if (ij2(k).eq.ii) then
              k2save=k
              igot2=1
            end if
30        continue
          ii=ipack3(ia(iord(5)),ia(iord(6)),ia(iord(7)))
          do 50 k=1,nij3
            if (ij3(k).eq.ii) then
              k3save=k
              igot3=1
            end if
50        continue
          ii=ipack3(ia(iord(8)),ia(iord(9)),ia(iord(10)))
          do 80 k=1,nij4
            if (ij4(k).eq.ii) then
              k4save=k
              igot4=1
            end if
80        continue
          if ((igot2.eq.1).and.(igot3.eq.1).and.(igot4.eq.1)) then
            sum=0.d0
            do 90 k=1,ije3(k3save)-ije3(k3save-1)
              ii=iand(igen(istrt31+ije3(k3save-1)+k),z'ff')
              do 100 l=1,ije4(k4save)-ije4(k4save-1)
                jj=iand(igen(istrt41+ije4(k4save-1)+l),z'ff')
                kk=ior(ishft(ii,8),jj)
                do 110 m=1,ije2(k2save)-ije2(k2save-1)
                  if (kk.eq.                                            &
             iand(igen(istrt21+ije2(k2save-1)+m),z'ffff')) then
                    sum=sum+                                            &
               gen(istrt22+ije2(k2save-1)+m)*                           &
               gen(istrt32+ije3(k3save-1)+k)*                           &
               gen(istrt42+ije4(k4save-1)+l)
                  end if
110             continue
100           continue
90          continue
            gen(istrt12+ije1(i-1)+j)=sum
          end if
20      continue
10    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine diage6 (name2,name3,name4,name5,namet,iord)
!--------------------------------------------------------------------------------------
!
!                  6x(4+4+4)=12
!
!--------------------------------------------------------------------------------------
      use ipack_mod, only: ipack3
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use cij4_mod
      use cij5_mod
      use idv_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "cfiles.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij2.h"
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
!     common /cij4/ nij4,ij4(ndij),ije4(0:ndij),iii4,iiv4
#include "cij4.h"
!     common /cij5/ nij5,ij5(ndij),ije5(0:ndij),iii5,iiv5
#include "cij5.h"
!     common /idv/ iv(200),nvalo
#include "idv.h"
#endif
      dimension iord(12),ia(12)

      ! store into module....
      !ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)

      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage6: input tape does not exist, name2=',name2
        call quit('diage6: input tape name2 does not exist !')
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage6: input tape does not exist, name3=',name3
        call quit('diage6: input tape name3 does not exist !')
      end if
      call files_ccsd (name4,ntape4,new,4,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage6 input tape does not exist, name4=',name4
        call quit('diage6: input tape name4 does not exist !')
      end if
      call files_ccsd (name5,ntape5,new,5,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage6: input tape does not exist, name5=',name5
        call quit('diage6: input tape name5 does not exist!')
      end if
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)
      istrt41=nstrt(1,ntape4-1)+2*nij4+4
      istrt42=nstrt(2,ntape4-1)
      istrt51=nstrt(1,ntape5-1)+2*nij5+4
      istrt52=nstrt(2,ntape5-1)
      call tmplt12 (namet,ntapet)
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1)
      call tmpcng (nij3,ij3,ije3,ntape3,2,3)
      call tmpcng (nij4,ij4,ije4,ntape4,2,3)
      call tmpcng (nij5,ij5,ije5,ntape5,2,3)
      do 10 i=1,nij1
        do 60 k=5,0,-1
          ii=iand(ishft(ij1(i),-4*k),z'f')
          ia(6-k)=iv(ii)
60      continue
        do 20 j=1,ije1(i)-ije1(i-1)
          ii=igen(istrt11+ije1(i-1)+j)
          do 70 k=5,0,-1
            jj=iand(ishft(ii,-4*k),z'f')
            ia(12-k)=iv(jj)
70        continue
          igot2=0
          igot3=0
          igot4=0
          igot5=0
          ii=ipack3(ia(iord(1)),ia(iord(2)),ia(iord(3)))
          do 30 k=1,nij2
            if (ij2(k).eq.ii) then
              k2save=k
              igot2=1
            end if
30        continue
          ii=ipack3(ia(iord(4)),ia(iord(5)),ia(iord(6)))
          do 50 k=1,nij3
            if (ij3(k).eq.ii) then
              k3save=k
              igot3=1
            end if
50        continue
          ii=ipack3(ia(iord(7)),ia(iord(8)),ia(iord(9)))
          do 80 k=1,nij4
            if (ij4(k).eq.ii) then
              k4save=k
              igot4=1
            end if
80        continue
          ii=ipack3(ia(iord(10)),ia(iord(11)),ia(iord(12)))
          do 120 k=1,nij5
            if (ij5(k).eq.ii) then
              k5save=k
              igot5=1
            end if
120       continue
          if ((igot2.eq.1).and.(igot3.eq.1).and.(igot4.eq.1).and.(igot5.eq.1)) then
            sum=0.d0
            do 90 k=1,ije3(k3save)-ije3(k3save-1)
              ii=iand(igen(istrt31+ije3(k3save-1)+k),z'ff')
              do 100 l=1,ije4(k4save)-ije4(k4save-1)
                jj=iand(igen(istrt41+ije4(k4save-1)+l),z'ff')
                do 130 m=1,ije5(k5save)-ije5(k5save-1)
                  kk=iand(igen(istrt51+ije5(k5save-1)+l),z'ff')
                  ll=ipack3(ii,jj,kk)
                  do 110 n=1,ije2(k2save)-ije2(k2save-1)
                    if (ll.eq.                                          &
               iand(igen(istrt21+ije2(k2save-1)+n),z'ffffff')) then
                      sum=sum+                                          &
                 gen(istrt22+ije2(k2save-1)+n)*                         &
                 gen(istrt32+ije3(k3save-1)+k)*                         &
                 gen(istrt42+ije4(k4save-1)+l)*                         &
                 gen(istrt52+ije5(k5save-1)+m)
                    end if
110               continue
130             continue
100           continue
90          continue
            gen(istrt12+ije1(i-1)+j)=sum
          end if
20      continue
10    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine diffmax(name1,name2,names,tmax,itmax)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      call add(name2,1.d0,name1,-1.d0,names)
      call fmax(names,tmax,itmax)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine diis(iter,icnvg,i12,name1n,name1o,name2n,name2o)
!-----------------------------------------------------------------------------
!
! DIIS convergence accelerator (there is also RLE in the code).
!
! Entering parameters:  
!           iter - actual iteration
!           icnvg -
!           i12 -
!           ...
!
! called from chu,chv,csd,csu,csv
!
!-----------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use dat_mod
      use cij1_mod
      use cij2_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "cfiles.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij2.h"
#endif

! ... miro TODO: introduce dynamic alloc 
      dimension r(6,6),rr(8,8),vec(100)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      character*4 name1(6),name2(6),name1di(6),name2di(6)
      data name1 /"t101","t102","t103","t104","t105","t106"/
      data name2 /"t201","t202","t203","t204","t205","t206"/
      data name1di /"t1d1","t1d2","t1d3","t1d4","t1d5","t1d6"/
      data name2di /"t2d1","t2d2","t2d3","t2d4","t2d5","t2d6"/
      maxdim=6
      if (iter.eq.1) then
        idim=0
        irdim=0
        do i=1,maxdim
        do j=1,maxdim
          r(j,i)=0.
        enddo
        enddo
      end if
      idim=idim+1
      irdim=min0(irdim+1,maxdim)
      if(idim.gt.maxdim) idim=idim-maxdim
      do id=1,irdim
      do i=1,id
        r(i,id)=0.
      enddo
      enddo
!  store current t and delta(t) vector.
!  calc r matrix using all delta(t) vectors.
      if (i12/10.gt.0) then
        call add (name1n,1.d0,name1o,-1.d0,"s2  ")
        call move(name1o,name1(idim))
        call move("s2  ",name1di(idim))
        do id=1,irdim
        call files_ccsd(name1di(id),ntape1,new,1,-1,-1)
        nsize=ije1(nij1)
        istrt1=nstrt(2,ntape1-1)
        do i=1,id
          call files_ccsd(name1di(i),ntape2,new,2,-1,-1)
          istrt2=nstrt(2,ntape2-1)
          r(i,id)=ddot(nsize,gen(istrt1+1),1,gen(istrt2+1),1)
!         write(*,1000) i,(gen(istrt1+j),j=1,nsize)
1000      format('diis vector1',i4/(5f12.6))
!         write(*,*) 'i,id,r=',i,id,r(i,id)
        enddo
        enddo
      endif
      if (mod(i12,10).gt.0) then
        call add (name2n,1.d0,name2o,-1.d0,"s2  ")
        call move(name2o,name2(idim))
        call move("s2  ",name2di(idim))
        do id=1,irdim
        call files_ccsd(name2di(id),ntape1,new,1,-1,-1)
        nsize=ije1(nij1)
        istrt1=nstrt(2,ntape1-1)
        do i=1,id
          call files_ccsd(name2di(i),ntape2,new,2,-1,-1)
          istrt2=nstrt(2,ntape2-1)
          r(i,id)=r(i,id)+ddot(nsize,gen(istrt1+1),1,   &
                               gen(istrt2+1),1)
!         write(*,1100) i,(gen(istrt1+j),j=1,nsize)
1100      format('diis vector2',i4/(5f12.6))
          r(id,i)=r(i,id)
!         write(*,*) 'i,id,r=',i,id,r(i,id)
        enddo
        enddo
      endif
!     
      do i=1,irdim+2
      do j=1,irdim+1
        rr(j,i)=0.
      enddo
      enddo
      do i=1,irdim
      do j=1,irdim
        rr(j+1,i+1)=r(j,i)
      enddo
      enddo
      do i=2,irdim+1
        rr(i,1)=-1.
        rr(1,i)=-1.
      enddo
      rr(1,irdim+2)=-1.
!     write(*,100) irdim,((rr(j,i),j=1,irdim+1),i=1,irdim+2)
100   format(' irdim, rr=',i4,/(5f12.6))
      call minv(rr,irdim+1,maxdim+2,vec,det,0.d0,1,0)
      do i=1,irdim
      vec(i)=rr(i+1,irdim+2)
      enddo
!     write (*,*) 'diis vec=',(vec(i),i=1,irdim)
!
!    calc diis correction. 
!
      if (i12/10.gt.0) then
        do i=1,irdim
          if (i.eq.1) then
            call cmult(name1di(1),vec(1),"t1d0")
          else
            call add("t1d0",1.d0,name1di(i),vec(i),"t1d0")
          endif
!       write(*,*) 'add dt1(',i,'), vec(i)=',vec(i)
!         call prt(name1di(i),name1di(i),0)
!         call prt("t1d0","t1d0 dt1",0)
        enddo
        call fmax("t1d0",t1dimx,it1dimx)
!       call add(name1n,1.d0,name1o,-1.d0,name1di(idim))
!    calc and stor new t1.
        do i=1,irdim
          call add("t1d0",1.d0,name1(i),vec(i),"t1d0")
!       write(*,*) 'add t1(',i,'), vec(i)=',vec(i)
!         call prt(name1(i),name1(i),0)
!         call prt("t1d0","t1d0  t1",0)
        enddo
        call move("t1d0",name1n)
        call fmax(name1n,t1mx,it1mx)
      endif
      if (mod(i12,10).gt.0) then
        do i=1,irdim
          if(i.eq.1) then
            call cmult(name2di(1),vec(1),"t2d0")
          else
            call add("t2d0",1.d0,name2di(i),vec(i),"t2d0")
          endif
        enddo
        call fmax("t2d0",t2dimx,it2dimx)
      write(6,2400) t2dimx,idec(it2dimx),t1dimx,idec(it1dimx)
2400  format(' diis: diffmax t, t1=',f15.9,2x,i12,f17.9,2x,i6)
        icnvg=0
        if (max(abs(t1dimx),abs(t2dimx)).lt.10.**(-iopt(3))) icnvg=1
!       call add(name2n,1.d0,name2o,-1.d0,name2di(idim))
!    calc and stor new t2.
        do i=1,irdim
          call add("t2d0",1.d0,name2(i),vec(i),"t2d0")
        enddo
        call move("t2d0",name2n)
        call fmax(name2n,t2mx,it2mx)
        write(6,2500) t2mx,idec(it2mx),t1mx,idec(it1mx)
2500    format(' diis: max t, t1=',f15.9,2x,i12,f17.9,2x,i6)
!       call prt(name2n,"t2diis  ",0)
      endif
      return
      end

#if defined HSFSCC_USE_OLD_DIIS
!miro: I had to cut it off otherwise it's hurting dynamical allocation mode !
! ***********************************************
!      subroutine diis(iter,i12,name1n,name1o,name2n,name2o)
! old diis_o.f file
      subroutine diis_notused(iter,i12,name1n,name1o,name2n,name2o)
!
! miro : when not used, should it be deleted ?
!
      implicit real*8 (a-h,o-z)
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /dat/   norb,ic(120),icv(120),idum2,eps(120),iopt(40)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      dimension r(20,20),rr(22,22),vec(400)

      character*4 name1di(20),name2di(20)
      data name1di /"t101","t102","t103","t104","t105","t106",       &
      "t107","t108","t109","t110","t111","t112","t113","t114",       &
      "t115","t116","t117","t118","t119","t120"/
      data name2di /"td01","td02","td03","td04","td05","td06",       &
      "td07","td08","td09","td10","td11","td12","td13","td14",       &
      "td15","td16","td17","td18","td19","td20"/
      maxdim=20
      if (iter.eq.1) then
        idim=0
        irdim=0
        do i=1,20
        do j=1,20
          r(j,i)=0.
        enddo
        enddo
      end if
      idim=idim+1
      irdim=min0(irdim+1,maxdim)
      if(idim.gt.maxdim) idim=idim-maxdim
      do i=1,irdim
        r(i,idim)=0.
      enddo
      if (i12/10.gt.0) then
        call add (name1n,1.d0,name1o,-1.d0,"s2  ")
        call move("s2  ",name1di(idim))
        call files_ccsd(name1di(idim),ntape1,new,1,-1,-1)
        nsize=ije1(nij1)
        istrt1=nstrt(2,ntape1-1)
        do i=1,irdim
          call files_ccsd(name1di(i),ntape2,new,2,-1,-1)
          istrt2=nstrt(2,ntape2-1)
          r(i,idim)=r(i,idim)+                                          &
           ddot(nsize,gen(istrt1+1),1,gen(istrt2+1),1)
          write(*,1000) i,(gen(istrt1+j),j=1,nsize)
1000      format('diis vector1',i4/(5f12.6))
          write(*,*) 'i,idim,r=',i,idim,r(i,idim)
        enddo
      endif
      if (mod(i12,10).gt.0) then
        call add (name2n,1.d0,name2o,-1.d0,"s2  ")
        call move("s2  ",name2di(idim))
        call files_ccsd(name2di(idim),ntape1,new,1,-1,-1)
        nsize=ije1(nij1)
        istrt1=nstrt(2,ntape1-1)
        do i=1,irdim
          call files_ccsd(name2di(i),ntape2,new,2,-1,-1)
          istrt2=nstrt(2,ntape2-1)
          r(i,idim)=r(i,idim)+ddot(nsize,gen(istrt1+1),1,               &
      gen(istrt2+1),1)
          write(*,1100) i,(gen(istrt1+j),j=1,nsize)
1100      format('diis vector2',i4/(5f12.6))
          r(idim,i)=r(i,idim)
          write(*,*) 'i,idim,r=',i,idim,r(i,idim)
        enddo
      endif
!     
      do i=1,irdim+2
      do j=1,irdim+1
        rr(j,i)=0.
      enddo
      enddo
      do i=1,irdim
      do j=1,irdim
        rr(j+1,i+1)=r(j,i)
      enddo
      enddo
      do i=2,irdim+1
        rr(i,1)=-1.
        rr(1,i)=-1.
      enddo
      rr(1,irdim+2)=-1.
!     write(*,100) irdim,((rr(j,i),j=1,irdim+1),i=1,irdim+2)
100   format(' irdim, rr=',i4,/(5f12.6))
      call minv(rr,irdim+1,22,vec,det,0.d0,1,0)
      do i=1,irdim
      vec(i)=rr(i+1,irdim+2)
      enddo
      write (*,*) 'vec ',(vec(i),i=1,irdim)
!
!    calc diis correction. update err vector idim.
!
      if (i12/10.gt.0) then
        do i=1,irdim
          call add(name1n,1.d0,name1di(i),vec(i),name1n)
        enddo
        call add(name1n,1.d0,name1o,-1.d0,name1di(idim))
      endif
      if (mod(i12,10).gt.0) then
        do i=1,irdim
          call add(name2n,1.d0,name2di(i),vec(i),name2n)
        enddo
        call add(name2n,1.d0,name2o,-1.d0,name2di(idim))
      endif
      return
      end
#endif

