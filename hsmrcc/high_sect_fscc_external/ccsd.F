      subroutine clear(name,iii,iiv)
      implicit real*8 (a-h,o-z)
c  zero out file (create if does not exist)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)

#include "cij1.h"
!      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

      integer indx_nstrt ! mi add

      call tmplt(iii,iiv,name,ntape)
      nsize=ije1(nij1)

! problem: relocation truncated to fit: R_X86_64_32S against symbol `cfiles_' defined in COMMON...
!      call dcopy(nsize,0.d0,0,
!     1 gen(nstrt(2,ntape-1)+1),1)

!      call dcopy(nsize,0.d0,0,gen(nstrt(2,ntape-1)+1),1)
      indx_nstrt=nstrt(2,ntape-1)+1
      call dcopy(nsize,0.d0,0,gen(indx_nstrt),1)

      return
      end

      subroutine clear6 (name,iii,iiv)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)

      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

      call tmplt6 (iii,iiv,name,ntape)
      nsize=ije1(nij1)

      call dcopy (nsize,0.d0,0,gen(nstrt(2,ntape-1)+1),1)
      call prt6 (name,name,2)
      return
      end

      subroutine clear8 (name,iii,iiv)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)

#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

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


C  *******************************************************************
      function find2 (j1,j2)                
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "gen.h"

!      common  // gen(nbig)
!      common /i/ igen(nbig)

#include "strt.h"
!      common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!     1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!     2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6

      find2=0.d0
      id=ior(ishft(j1,8),j2)
      do 10 i=1,nsize1
        if (id.eq.igen(istrt11+i)) then
          find2=gen(istrt12+i)
        end if
10    continue
      return
      end


      function find4 (j1,j2,j3,j4)          
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "gen.h"
!      common  // gen(nbig)
!      common /i/ igen(nbig)

#include "strt.h"
!      common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!     1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!     2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6

      find4=0.d0
      id=ior(ior(ishft(j1,24),ishft(j2,16)),ior(ishft(j3,8),j4))
      do 10 i=1,nsize2
        if (id.eq.igen(istrt21+i)) then
          find4=gen(istrt22+i)
        end if
10    continue
      return
      end


      function find6 (j1,j2,j3,j4,j5,j6)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "gen.h"
!      common  // gen(nbig)
!      common /i/ igen(nbig)

      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3

#include "strt.h"
!      common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!     1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!     2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6


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
        if (ind.eq.iand(ii,'ffffff'x)) then
          find6=gen(istrt32+ije3(isave-1)+i)
        end if
20    continue
      return
      end


      function find8 (j1,j2,j3,j4,j5,j6,j7,j8)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common  // gen(nbig)
      common /i/ igen(nbig)
      common /cij4/ nij4,ij4(ndij),ije4(0:ndij),iii4,iiv4
      common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
     1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
     2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
      find8=0.d0
      id=ior(ior(ishft(j1,24),ishft(j2,16)),ior(ishft(j3,8),j4))
      ind=ior(ior(ishft(j5,24),ishft(j6,16)),ior(ishft(j7,8),j8))
      do 10 i=1,nij4
        if (id.eq.ij4(i)) then
          isave=i
        end if
10    continue
      lngth=ije4(isave)-ije4(isave-1)
      do 20 i=1,lngth
        if (ind.eq.igen(istrt41+ije4(isave-1)+i)) then
          find8=gen(istrt42+ije4(isave-1)+i)
        end if
20    continue
      return
      end


      function find10 (ja,jb,jc,jd,je,jf,jg,jh,ji,jj)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common  // gen(nbig)
      common /i/ igen(nbig)
      common /cij5/ nij5,ij5(ndij),ije5(0:ndij),iii5,iiv5
      common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
     1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
     2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
      dimension jj1(10),kk(10)
      ipack (i,j,k,l,m)=ior(ishft(i,16),
     1ior(ior(ishft(j,12),ishft(k,8)),ior(ishft(l,4),m)))
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
      idb=ipack(kk(1),kk(2),kk(3),kk(4),kk(5))
      ind=ipack(kk(6),kk(7),kk(8),kk(9),kk(10))
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

      function find12 (ja,jb,jc,jd,je,jf,jg,jh,ji,jj,jk,jl)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common  // gen(nbig)
      common /i/ igen(nbig)
      common /cij6/ nij6,ij6(ndij),ije6(0:ndij),iii6,iiv6
      common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
     1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
     2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
      dimension jj1(12),kk(12)
      ipack (i,j,k,l,m,n)=ior(
     1ior(ishft(i,20),ior(ishft(j,16),ishft(k,12))),
     2ior(ishft(l,8),ior(ishft(m,4),n)))
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
      idb=ipack(kk(1),kk(2),kk(3),kk(4),kk(5),kk(6))
      ind=ipack(kk(7),kk(8),kk(9),kk(10),kk(11),kk(12))
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

! renamed due to name overlaping with relccsd routine
!      subroutine fold(name2,namef)
      subroutine fold_ccsd(name2,namef)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
c  two-electron part modified to include chp, chq.
c  880922.
c  all four indices now checked. incomplete space implemented
c  for chp, chq only for exclusion of whole symmetries.
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /rsp1/  nrsp1,irsp1(100)
      common /rsp2/  nrsp2,irsp2(500)
      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
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
c   1-electron t.
        do 120 i=1,nsize2
	i1=iand(igen(istrt21+i),'ff'X)
        if (icv(i1).eq.0) go to 120
        if (iopt(9).ne.0) then
          do 105 ir=1,nrsp1
105       if (igen(istrt21+i).eq.irsp1(ir)) go to 120
      endif
c   move p to p excitations to folded
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
c   2-electron t
        do 180 i=1,nsize2
c   transition from r space eliminated
        if (iopt(9).ne.0) then
         i_4hthpf=4hthpf
         i_4hthqf=4hthqf
!        if (namef.eq.4hthpf) then
         if (namef.eq.i_4hthpf) then
          i12=ior(iand(ishft(igen(istrt21+i),-16),'ff00'X),
     1              iand(igen(istrt21+i),'ff'X))
	  else if (namef.eq.i_4hthqf) then
	    i12=ior(iand(ishft(igen(istrt21+i),-16),'ff00'X),
     1              iand(ishft(igen(istrt21+i),-8),'ff'X))
          else
	    i12=iand(ishft(igen(istrt21+i),-16),'ffff'X)
          endif
          do 130 ir=1,nrsp2
          if (irsp2(ir).eq.i12) then
            gen(istrt22+i)=0.d0
            go to 180
          endif
130       continue
        endif
c   move p to p transitions to folded
	i1=iand(ishft(igen(istrt21+i),-24),'ff'X)
	i2=iand(ishft(igen(istrt21+i),-16),'ff'X)
	i3=iand(ishft(igen(istrt21+i),-8) ,'ff'X)
	i4=iand(      igen(istrt21+i)     ,'ff'X)
	if ((icv(i1).eq.0).or.(icv(i2).eq.0)) go to 180
	if ((icv(i3).eq.0).or.(icv(i4).eq.0)) go to 180
        if (iopt(9).eq.0) go to 170
        i_4hthpf=4hthpf
        i_4hthqf=4hthqf
	if (namef.eq.i_4hthpf) then
	  i34=iand(ishft(igen(istrt21+i),-8),'ffff'X)
	else if (namef.eq.i_4hthqf) then
	  i34=ior(iand(ishft(igen(istrt21+i),-8),'ff00'X),
     1            iand(igen(istrt21+i),'ff'X))
        else
	  i34=iand(igen(istrt21+i),'ffff'X)
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
      subroutine fold4 (nameu,namef,ind,ninds)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
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
          call parse (iiv2,jj1(1),jj1(2))
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
      subroutine fold6 (nameu,namef,ind,iparm)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      dimension jj1(6)
      call files_ccsd (nameu,ntapeu,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'fold6 input tape does not exist, ',nameu
        stop                    
      end if
      if (iparm.eq.1) then
        call parse6 (iiv2,jj1(1),jj1(2),jj1(3),jj1(4),jj1(5),jj1(6))
        jj1(ind)=1
        iiv1=100000*jj1(1)+10000*jj1(2)+1000*jj1(3)
     1       +100*jj1(4)+10*jj1(5)+jj1(6)
        iii1=iii2
        call tmplt6 (iii1,iiv1,namef,ntapef)
      else if (iparm.eq.2) then
        call files_ccsd (namef,ntapef,new,1,-1,-1)
        if (new.ne.0) then
          write (6,*) 'fold6 input tape does not exist, ',namef
          stop                     
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
                stop                       
              end if
            end if
          end do
20      continue
10    continue
      call prtime ("fold6   ",iparm,2)
      return
      end
      subroutine fold8 (nameu,namef,ind,iparm)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      dimension ia(8)
      call files_ccsd (nameu,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'fold8 input tape does not exist, '
        stop                    
      end if
      if (iparm.eq.1) then
        call parse8 (iiv2,ia(1),ia(2),ia(3),ia(4),
     1  ia(5),ia(6),ia(7),ia(8))
        ia(ind)=1
        iii1=iii2
        iiv1=10000000*ia(1)+1000000*ia(2)+100000*ia(3)+10000*ia(4)+
     1  1000*ia(5)+100*ia(6)+10*ia(7)+ia(8)
        call tmplt8 (iii1,iiv1,namef,ntape1)
      else if (iparm.eq.2) then
        call files_ccsd (namef,ntape1,new,1,-1,-1)
        if (new.ne.0) then
          write (6,*) 'fold8 input tape does not exist, '
          stop                    
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
                  gen(istrt12+ije1(i-1)+j)=
     1            gen(istrt22+ije2(idc-1)+indc)
                else
                  gen(istrt22+ije2(idc-1)+indc)=
     1            gen(istrt22+ije2(idc-1)+indc)-
     2            gen(istrt12+ije1(i-1)+j)
                end if
                ididit=1
              else if (ind2.lt.ind1) then
                indc=indc+1
                ind2=igen(istrt21+ije2(idc-1)+indc)
              else
                write (6,*) 'folded not found - fold8 '
                stop
              end if
            else if (id2.lt.id1) then
              idc=idc+1
              id2=ij2(idc)
              indc=0
            else
              write (6,*) 'folded not found - fold8 '
              stop
            end if
          end do
20      continue
10    continue
      return
      end

      subroutine msub1 (name2,name3,namet)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /ctmp/ tmp(ndp4),itmp(ndp4)
      common /dpd/ nir,igmt(8,8),ici(200)
      dimension ia(8),ib(10000),ibe(0:10000)
      ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)
      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub1 input file does not exist'
        stop                       
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub1 input file does not exist'
        stop                       
      end if
      if (((iii2-10*(iii2/10)).ne.(iii3-10*(iii3/10))).or.
     1    ((iiv2-10*(iiv2/10)).ne.(iiv3-10*(iiv3/10)))) then
        write (6,*) 'contracted indices do not match - msub1'
        stop
      end if
      if ((iii2.lt.222222).or.(iii2.gt.333333)) then
        write (6,*) 'msub1 called with inconsistent iii2'
        stop
      end if
      if ((iii3.lt.2222).or.(iii3.gt.3333)) then
        write (6,*) 'msub1 called with inconsistent iii3'
        stop
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
          ia(4-j)=iand(ishft(ij2(i),-8*j),'ff'x)
20      continue
        ii=0
        n=0
        do 30 j=1,ije2(i)-ije2(i-1)
          jj=iand(ishft(igen(istrt21+ije2(i-1)+j),-8),'ff'x)
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
            ia(8-m)=iand(ishft(ij3(k),-8*m),'ff'x)
70          continue
            ii=1
            do 80 m=1,8
              ii=igmt(ii,ici(ia(m)))
80          continue
            if (ii.eq.1) then
              lngth=ije3(k)-ije3(k-1)
              if (lngth.ne.(ibe(j)-ibe(j-1))) then
                write (6,*) 'mismatch in block lengths - msub1'
                  stop
              end if
              elnew=ddot(lngth,gen(istrt22+ije2(i-1)+ibe(j-1)+1),1,
     1        gen(istrt32+ije3(k-1)+1),1)
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
                  stop
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
                  write (6,*) 'msub1 - label not in template'
                  stop
                end if
              end do
              gen(istrt12+ije1(idc-1)+indc)=elnew
            end if
60        continue
50      continue
10    continue
      return
      end
      subroutine msub2 (name2,name3,namet)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /ctmp/ tmp(ndp4),itmp(ndp4)
      common /dpd/ nir,igmt(8,8),ici(200)
      dimension ia(8),ib(10000),ibe(0:10000)
      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub2 input file does not exist'
        stop                       
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub2 input file does not exist'
        stop                       
      end if
      if ((iii2.lt.22222222).or.(iii2.gt.33333333)) then
        write (6,*) 'msub2 called with inconsistent iii2'
        stop
      end if
      if (iii3.gt.3333) then
        write (6,*) 'msub2 called with inconsistent iii3'
        stop
      end if
      if (iii3.gt.33) then
        if (((iii2-10000*(iii2/10000)).ne.iii3).or.
     1      ((iiv2-10000*(iiv2/10000)).ne.iiv3)) then
          write (6,*) 'contracted indices do not match - msub2'
          stop
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
              stop
            end if
          end do
          if (ije3(nij3).ne.(ije2(idc)-ije2(idc-1))) then
            write (6,*) 'mismatch in block lengths - msub2 '
            stop
          end if
          gen(istrt12+i)=ddot(ije3(nij3),gen(istrt22+ije2(idc-1)+1),1,
     1    gen(istrt32+1),1)
90      continue
      else if (iii3.le.33) then
        if (((iii2-100*(iii2/100)).ne.iii3).or.
     1      ((iiv2-100*(iiv2/100)).ne.iiv3)) then
          write (6,*) 'contracted indices do not match - msub2'
          stop
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
            ind3=iand(igen(istrt11+ije1(i-1)+j),'ffffff'x)
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
                stop
              end if
            end do
            lngth=ije2(imark)-ije2(imark-1)
            if (lngth.ne.ije3(nij3)) then
              write (6,*) 'mismatch in block lengths - msub2'
              stop
            end if
            gen(istrt12+ije1(i-1)+j)=ddot(lngth,
     1      gen(istrt22+ije2(imark-1)+1),1,gen(istrt32+1),1)
210       continue
200     continue
      end if
      return
      end
      subroutine msub3 (name2,name3,namet)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /ctmp/ tmp(ndp4),itmp(ndp4)
      common /dpd/ nir,igmt(8,8),ici(200)
      dimension ia(8),ib(10000),ibe(0:10000)
      ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)
      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub3 input file does not exist'
        stop                       
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub3 input file does not exist'
        stop                       
      end if
      if (((iii2-1000*(iii2/1000)).ne.(iii3-1000*(iii3/1000))).or.
     1    ((iiv2-1000*(iiv2/1000)).ne.(iiv3-1000*(iiv3/1000)))) then
        write (6,*) 'contracted indices do not match - msub3'
        stop
      end if
      if ((iii2.lt.22222222).or.(iii2.gt.33333333)) then
        write (6,*) 'msub3 called with inconsistent iii2'
        stop
      end if
      if ((iii3.lt.2222).or.(iii3.gt.3333)) then
        write (6,*) 'msub3 called with inconsistent iii3'
        stop
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
          ia(4-j)=iand(ishft(ij2(i),-8*j),'ff'x)
120     continue
        ii=0
        n=0
        do 130 j=1,ije2(i)-ije2(i-1)
          jj=iand(ishft(igen(istrt21+ije2(i-1)+j),-24),'ff'x)
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
                write (6,*) 'mismatch in block lengths - msub3'
                stop
              end if
              elnew=ddot(lngth,gen(istrt22+ije2(i-1)+ibe(j-1)+1),1,
     1        gen(istrt32+ije3(k-1)+1),1)
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
                  indold=iand(ii,'ffffff'x)  
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
      subroutine msub4 (name2,name3,namet)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /dpd/ nir,igmt(8,8),ici(200)
      dimension ia(8),ib(10000),ibe(0:10000)
      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub4 input file does not exist'
        stop                       
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub4 input file does not exist'
        stop                       
      end if
      if (((iii2-100*(iii2/100)).ne.(iii3-100*(iii3/100))).or.
     1    ((iiv2-100*(iiv2/100)).ne.(iiv3-100*(iiv3/100)))) then
        write (6,*) 'contracted indices do not match - msub4'
        stop
      end if
      if ((iii2.lt.22222222).or.(iii2.gt.33333333)) then
        write (6,*) 'msub4 called with inconsistent iii2'
        stop
      end if
      if ((iii3.lt.2222).or.(iii3.gt.3333)) then
        write (6,*) 'msub4 called with inconsistent iii3'
        stop
      end if
      iii1=100*(iii2/100)+iii3/100
      iiv1=100*(iiv2/100)+iiv3/100
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
          ia(4-j)=iand(ishft(ij2(i),-8*j),'ff'x)
310     continue
        ii=0
        n=0
        do 320 j=1,ije2(i)-ije2(i-1)
          jj=iand(ishft(igen(istrt21+ije2(i-1)+j),-16),'ffff'x)
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
          ia(5)=iand(ishft(ib(j),-8),'ff'x)
          ia(6)=iand(ib(j),'ff'x)
          do 340 k=1,nij3 
            ia(7)=iand(ishft(ij3(k),-8),'ff'x)
            ia(8)=iand(ij3(k),'ff'x)
            ii=1
            do 350 m=1,8
              ii=igmt(ii,ici(ia(m)))
350         continue
            if (ii.eq.1) then
              lngth=ije3(k)-ije3(k-1)
              if (lngth.ne.(ibe(j)-ibe(j-1))) then
                write (6,*) 'mismatch in block lengths - msub4'
                stop
              end if
              elnew=ddot(lngth,gen(istrt22+ije2(i-1)+ibe(j-1)+1),1,
     1        gen(istrt32+ije3(k-1)+1),1)
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
      subroutine msub5 (name2,name3,namet)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /ctmp/ tmp(ndp4),itmp(ndp4)
      common /dpd/ nir,igmt(8,8),ici(200)
      dimension ia(8),ib(10000),ibe(0:10000)
      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub5 input file does not exist'
        stop                       
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub5 input file does not exist'
        stop                       
      end if
      if (((iii2-100*(iii2/100)).ne.(iii3-100*(iii3/100))).or.
     1    ((iiv2-100*(iiv2/100)).ne.(iiv3-100*(iiv3/100)))) then
        write (6,*) 'contracted indices do not match - msub5'
        stop
      end if
      if ((iii2.lt.222222).or.(iii2.gt.333333)) then
        write (6,*) 'msub5 called with inconsistent iii2'
        stop
      end if
      if ((iii3.lt.222222).or.(iii3.gt.333333)) then
        write (6,*) 'msub5 called with inconsistent iii3'
        stop
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
          ia(4-k)=iand(ishft(ij2(i),-8*k),'ff'x)
420     continue
        do 430 j=1,nij3
          do 440 k=0,3
            ia(8-k)=iand(ishft(ij3(j),-8*k),'ff'x)
440       continue
          ii=1
          do 450 k=1,8
            ii=igmt(ii,ici(ia(k)))
450       continue
          if (ii.eq.1) then
            lngth=ije2(i)-ije2(i-1)
            if (lngth.ne.(ije3(j)-ije3(j-1))) then
              write (6,*) 'mismatch in block lengths - msub5'
              stop
            end if
            elnew=ddot(lngth,gen(istrt22+ije2(i-1)+1),1,
     1      gen(istrt32+ije3(j-1)+1),1)
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
      subroutine msub6 (name2,name3,namet)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /ctmp/ tmp(ndp4),itmp(ndp4)
      common /dpd/ nir,igmt(8,8),ici(200)
      dimension ia(8),ib(10000),ibe(0:10000)
      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub6 input file does not exist'
        stop                       
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'msub6 input file does not exist'
        stop                       
      end if
      if (((iii2-10*(iii2/10)).ne.(iii3-10*(iii3/10))).or.
     1    ((iiv2-10*(iiv2/10)).ne.(iiv3-10*(iiv3/10)))) then
        write (6,*) 'contracted indices do not match - msub6'
        stop
      end if
      if ((iii2.lt.22222222).or.(iii2.gt.33333333)) then
        write (6,*) 'msub6 called with inconsistent iii2'
        stop
      end if
      if (iii3.gt.33) then
        write (6,*) 'msub6 called with inconsistent iii3'
        stop
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
          ia(4-j)=iand(ishft(ij2(i),-8*j),'ff'x)
120     continue
        ii=0
        n=0
        do 130 j=ije2(i-1)+1,ije2(i)
          jj=iand(ishft(igen(istrt21+j),-8),'ffffff'x)
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
          ia(5)=iand(ishft(ib(j),-16),'ff'x)
          ia(6)=iand(ishft(ib(j),-8),'ff'x)
          ia(7)=iand(ib(j),'ff'x)
          do 160 k=1,nij3 
            ia(8)=ij3(k)
            ii=1
            do 180 m=1,8
              ii=igmt(ii,ici(ia(m)))
180         continue
            if (ii.eq.1) then
              lngth=ije3(k)-ije3(k-1)
              if (lngth.ne.(ibe(j)-ibe(j-1))) then
                write (6,*) 'mismatch in block lengths - msub6'
                stop
              end if
              elnew=ddot(lngth,gen(istrt22+ije2(i-1)+ibe(j-1)+1),1,
     1        gen(istrt32+ije3(k-1)+1),1)
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
      subroutine mult(name2,name3,namet,mmm)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      call parse(mmm,idum,m1,m2,m3)
      call files_ccsd(name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) ' mult input tape does not exist, ',name2
        call abext('mult    ')
      endif
      call files_ccsd(name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) ' mult input tape does not exist, ',name3
        call abext('mult    ')
      endif
      n=10**m3
      if (((iii2-n*(iii2/n)).ne.(iii3-n*(iii3/n))).or.
     1    ((iiv2-n*(iiv2/n)).ne.(iiv3-n*(iiv3/n)))) then
        write (6,*) 'Contracted indices do not match in mult'
        write (6,*) 'iii2, iii3, iiv2, iiv3 = '
        write (6,*) iii2,iii3,iiv2,iiv3
        call abext('mult    ')
      end if
      iii=10**m2*(iii2/10**m3)+iii3/10**m3
      iiv=10**m2*(iiv2/10**m3)+iiv3/10**m3
      call tmplt(iii,iiv,namet,ntapet)
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)
      nsize1=ije1(nij1)
      nsize2=ije2(nij2)
      nsize3=ije3(nij3)
c
      if (m3.ne.2) then
        call tmpcng(nij1,ij1,ije1,ntapet,2,4-m2)
        call tmpcng(nij2,ij2,ije2,ntape2,2,4-m3)
        call tmpcng(nij3,ij3,ije3,ntape3,2,4-m3)
      endif
c
      call dcopy(nsize1,0.d0,0,gen(istrt12+1),1)
      kbb=0
      kcc=0
      ib0='ffff'X
      msk1=256**m2-1
      do 250 i1=1,nsize1
      ib=ishft(igen(istrt11+i1),-8*m2)
      ic=iand(igen(istrt11+i1),msk1)
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
2100    format('mismatch in length of records',2(1x,z6),
     1  ' of files ',2a8,'. check isym.')
        write (*,2200) n2,(igen(istrt21+i),i=ije2(kbb-1)+1,ije2(kbb))
        n3=ije3(kcc)-ije3(kcc-1)
        write (*,2200) n3,(igen(istrt31+i),i=ije3(kcc-1)+1,ije3(kcc))
2200    format(' length and indices',i5/(10(1x,z8)))
        call abext('mult    ')
      endif
      n20=istrt22+ije2(kbb-1)+1
      n30=istrt32+ije3(kcc-1)+1
      gen(istrt12+i1)=ddot(n2,gen(n20),1,
     1 gen(n30),1)
250   continue
      return
      end
      subroutine mult6 (name2,name3,namet,id)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /dpd/ nir,igmt(8,8),ici(200)
      dimension ia(8),itc(10000),itce(0:10000)
      ipack2(i,j)=ior(ishft(i,8),j)
      ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)
      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'mult6 input tape does not exist, ',name2
        stop                      
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'mult6 input tape does not exist, ',name3
        stop                       
      end if
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)
      nsize2=ije2(nij2)
      nsize3=ije3(nij3)
      if ((id.ne.0).and.(id.ne.5)) then
        if (id.eq.1) then
          if (((iii2-10*(iii2/10)).ne.(iii3-10*(iii3/10))).or.   
     1        ((iiv2-10*(iiv2/10)).ne.(iiv3-10*(iiv3/10)))) then
            write (6,*) 'contracted indices do not match - mult6'
            write (6,*) 'iii2=',iii2,' iii3=',iii3
           write (6,*) 'iiv2=',iiv2,' iiv3=',iiv3
            stop
          end if
          call tmpcng (nij2,ij2,ije2,ntape2,2,3)
          call tmpcng (nij3,ij3,ije3,ntape3,2,3)
          iii1=1000*(iii2/10)+iii3/10
          iiv1=1000*(iiv2/10)+iiv3/10
          call tmplt6 (iii1,iiv1,namet,ntapet)
          istrt11=nstrt(1,ntapet-1)+2*nij1+4
          istrt12=nstrt(2,ntapet-1)
          msk1='ffffff'x
          msk2='ffffff'x
        else if (id.eq.2) then
          if (((iii2-1000*(iii2/1000)).ne.(iii3-1000*(iii3/1000))).or.
     1        ((iiv2-1000*(iiv2/1000)).ne.(iiv3-1000*(iiv3/1000)))) then
            write (6,*) 'contracted indices do not match - mult6'
            write (6,*) 'iii2=',iii2,' iii3=',iii3
            write (6,*) 'iiv2=',iiv2,' iiv3=',iiv3
            stop
          end if
          if (iii2.gt.3333) then
            call tmpcng (nij3,ij3,ije3,ntape3,2,1)
            iii1=10*(iii2/1000)+iii3/1000
            iiv1=10*(iiv2/1000)+iiv3/1000
            msk1='ffffff'x
            msk2='ff'x
            call tmplt (iii1,iiv1,namet,ntapet)
            istrt11=nstrt(1,ntapet-1)+2*nij1+4
            istrt12=nstrt(2,ntapet-1)
            call tmpcng (nij1,ij1,ije1,ntapet,2,3)
          else if (iii2.le.3333) then
            call tmpcng (nij2,ij2,ije2,ntape2,2,1)
            iii1=1000*(iii2/1000)+iii3/1000
            iiv1=1000*(iiv2/1000)+iiv3/1000
            msk1='ff'x
            msk2='ffffff'x
            call tmplt (iii1,iiv1,namet,ntapet)
            istrt11=nstrt(1,ntapet-1)+2*nij1+4
            istrt12=nstrt(2,ntapet-1)
            call tmpcng (nij1,ij1,ije1,ntapet,2,1)
          end if
        else if (id.eq.3) then
          if (((iii2-100*(iii2/100)).ne.(iii3-100*(iii3/100))).or.    
     1       ((iiv2-100*(iiv2/100)).ne.(iiv3-100*(iiv3/100)))) then
            write (6,*) 'contracted indices do not match - mult6'
            write (6,*) 'iii2=',iii2,' iii3=',iii3
            write (6,*) 'iiv2=',iiv2,' iiv3=',iiv3
            stop
          end if
          if (iii2.gt.3333) then
            call tchng6 (nij2,ij2,ije2,ntape2,4)
            iii1=100*(iii2/100)+iii3/100
            iiv1=100*(iiv2/100)+iiv3/100
            call tmplt6 (iii1,iiv1,namet,ntapet)
            istrt11=nstrt(1,ntapet-1)+2*nij1+4
            istrt12=nstrt(2,ntapet-1)
            call tchng6 (nij1,ij1,ije1,ntapet,4)
            msk1='ffffffff'x
            msk2='ffff'x
          else if (iii2.le.3333) then
            call tchng6 (nij3,ij3,ije3,ntape3,4)
            iii1=10000*(iii2/100)+iii3/100
            iiv1=10000*(iiv2/100)+iiv3/100
            call tmplt6 (iii1,iiv1,namet,ntapet)
            istrt11=nstrt(1,ntapet-1)+2*nij1+4
            istrt12=nstrt(2,ntapet-1)
            call tchng6 (nij1,ij1,ije1,ntapet,2)
            msk1='ffff'x
            msk2='ffffffff'x
          end if
        else if (id.eq.4) then
          if (((iii2-1000*(iii2/1000)).ne.(iii3-1000*(iii3/1000))).or.
     1        ((iiv2-1000*(iiv2/1000)).ne.(iiv3-1000*(iiv3/1000)))) then
            write (6,*) 'contracted indices do not match - mult6'
            write (6,*) 'iii2=',iii2,' iii3=',iii3
            write (6,*) 'iiv2=',iiv2,' iiv3=',iiv3
            stop
          end if
          iii1=1000*(iii2/1000)+iii3/1000
          iiv1=1000*(iiv2/1000)+iiv3/1000
          call tmplt6 (iii1,iiv1,namet,ntapet)
          istrt11=nstrt(1,ntapet-1)+2*nij1+4
          istrt12=nstrt(2,ntapet-1)
          msk1='ffffff'x
          msk2='ffffff'x
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
96                format ('mismatch in length of records',2(1x,z6),
     1            ' of files ',2a8,'.check isym.')
                  write (6,95) n2,(igen(istrt21+im),
     1            im=ije2(kbsave-1)+1,ije2(kbsave))
95                format (' length and indices',i5/(10(1x,z8)))
                  stop                     
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
          if (((iii2-10000*(iii2/10000)).ne.iii3).or.   
     1        ((iiv2-10000*(iiv2/10000)).ne.iiv3)) then
            write (6,*) 'contracted indices do not match - mult6'
            write (6,*) 'iii2=',iii2,' iii3=',iii3
            write (6,*) 'iiv2=',iiv2,' iiv3=',iiv3
            stop
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
                write (6,*) 'mismatch check isym'
                stop                      
              end if
              n20=istrt22+ije2(j-1)+1
              gen(istrt12+i)=ddot(n2,gen(n20),1,gen(istrt32+1),1)
              i=i+1
              j=j+1
            else if (id1.gt.id2) then
              j=j+1
            else if (id1.lt.id2) then
              write (6,*) 'label not in template'
              stop                    
            end if
          end do   
        else if (iii3.le.33) then
          if (((iii2-100*(iii2/100)).ne.iii3).or.    
     1        ((iii2-100*(iii2/100)).ne.iii3)) then
            write (6,*) 'contracted indices do not match - mult6'
            write (6,*) 'iii2=',iii2,' iii3=',iii3
            write (6,*) 'iiv2=',iiv2,' iiv3=',iiv3
            stop
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
                write (6,*) 'mismatch check isym - mult6'
                stop                      
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
        if (((iii2-10*(iii2/10)).ne.(iii3-10*(iii3/10))).or.
     1      ((iiv2-10*(iiv2/10)).ne.(iiv3-10*(iiv3/10)))) then
          write (6,*) 'contracted indices do not match - mult6'
          stop
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
              ia(4-j)=iand(ishft(ij2(i),-8*j),'ff'x)
40          continue
            ii=0
            n=0
            do 50 j=1,ije2(i)-ije2(i-1)
              jj=iand(ishft(igen(istrt21+ije2(i-1)+j),-8),'ff'x)
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
                    write (6,*) 'mismatch in block lengths - mult6'
                    stop
                  end if
                  elnew=ddot(lngth,gen(istrt22+ije2(i-1)+itce(j-1)+1),1,
     1            gen(istrt32+ije3(k-1)+1),1)
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
                      stop
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
      subroutine mult8 (name2,name3,namet,iparm)
      implicit real*8 (a-h,o-z)
c  iparm=1 : 6x4=8
c  iparm=2 : 8x4=4,8x2=6
c  iparm=3 : 8x4=6
c  iparm=4 : 8x4=8
c  iparm=5 : 6x6=8
c  iparm=6 : 8x2=8
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
      subroutine abext(aa)
      implicit real*8 (a-h,o-z)
      character*8 aa
      integer*4 ii
      ii=99
      write(6,*) ' aborted from ',aa
      call abort()
      end
      subroutine add(name1,f1,name2,f2,namet)
      implicit real*8 (a-h,o-z)
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
c   namet=name1*f1+name2*f2
      if ((name1.eq.name2).and.(name1.eq.namet)) then
        write (6,*) 'all three names in add are identical'
        call abext('add     ')
      end if
      call files_ccsd(name1,ntape1,new,2,-1,-1)
      if (new.ne.0) then
        write(6,'('' add input file does not exist '',a8)') name1
        call abext('add     ')
      endif
      call files_ccsd(name2,ntape2,new,3,-1,-1)
      if (new.ne.0) then
        write(6,'('' add input file does not exist '',a8)') name2
        call abext('add     ')
      endif
      if((iii2.ne.iii3).or.(iiv2.ne.iiv3).or.(ije2(nij2).ne.ije3(nij3)))
     1 then
        write (6,*) ' added files do not match. names iii iiv nsize='
        write (6,123) name1,name2,iii2,iii3,iiv2,iiv3,ije2(nij2),
     1  ije3(nij3)
123     format(2a8,6i10)
        call abext('add     ')
      endif
      call files_ccsd(namet,ntapet,new,1,iii2,iiv2)
      nsize=ije2(nij2)
      if (name1.eq.namet) then
        if (f1.ne.1.d0) call rvmuls(nsize,gen(nstrt(2,ntape1-1)+1),
     1  f1,gen(nstrt(2,ntape1-1)+1))
        call daxpy(nsize,f2,gen(nstrt(2,ntape2-1)+1),
     1  1,gen(nstrt(2,ntape1-1)+1),1)
      elseif (name2.eq.namet) then
        if (f2.ne.1.d0) call rvmuls(nsize,gen(nstrt(2,ntape2-1)+1),
     1  f2,gen(nstrt(2,ntape2-1)+1))
        call daxpy(nsize,f1,gen(nstrt(2,ntape1-1)+1),
     1  1,gen(nstrt(2,ntape2-1)+1),1)
      else
        call rvmuls(nsize,gen(nstrt(2,ntape1-1)+1),
     1  f1,gen(nstrt(2,ntapet-1)+1))
        call daxpy(nsize,f2,gen(nstrt(2,ntape2-1)+1),
     1  1,gen(nstrt(2,ntapet-1)+1),1)
      endif
      return
      end
      subroutine cmult(namef,f,namet)
      implicit real*8 (a-h,o-z)
c   namet=f*namef
      
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
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
      function cputime(x)
      cputime = mclock()/1000000.
      return
      end
      subroutine dbdiag (name1,name2,name3,name4,irleop,isign)
      implicit real*8 (a-h,o-z)
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
CUK   common /crle/ trle(nrle,4),isize(4)
      common /crle/ trle1(nrle1),trle2(nrle2),trle3(nrle3),
     1 trle4(nrle4),isize(4)
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
          ka=iand(ishft(it,-8),'ff'x)
          kb=iand(it,'ff'x)
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
          ka=iand(ishft(it,-24),'ff'x)
          kb=iand(ishft(it,-16),'ff'x)
          kc=iand(ishft(it,-8),'ff'x)
          kd=iand(it,'ff'x)
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
          ka=iand(ishft(id,-16),'ff'x)
          kb=iand(ishft(id,-8),'ff'x)
          kc=iand(id,'ff'x)
          do 40 j=ije1(i-1)+1,ije1(i)
            ind=igen(istrt1+j)
            kd=iand(ishft(ind,-16),'ff'x)
            ke=iand(ishft(ind,-8),'ff'x)
            kf=iand(ind,'ff'x)
            trle3(j)=-1.d0*isign*
     1      (eps(ka)+eps(kb)+eps(kc)-eps(kd)-eps(ke)-eps(kf))
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
            ka=iand(ishft(id,-16),'ff'x)
            kb=iand(ishft(id,-8),'ff'x)
            kc=iand(id,'ff'x)
            do 60 j=ije1(i-1)+1,ije1(i)
              ind=igen(istrt1+j)
              kd=iand(ishft(ind,-16),'ff'x)
              ke=iand(ishft(ind,-8),'ff'x)
              kf=iand(ind,'ff'x)
              trle4(j)=-1.d0*isign*
     1        (eps(ka)+eps(kb)+eps(kc)-eps(kd)-eps(ke)-eps(kf))
              trle4(isize(4)+j)=gen(istrt2+j)
60          continue
50        continue
        else if (iii1.gt.333333) then
          do 70 i=1,nij1 
            id=ij1(i)             
            ka=iand(ishft(id,-24),'ff'x)
            kb=iand(ishft(id,-16),'ff'x)
            kc=iand(ishft(id,-8),'ff'x)
            kd=iand(id,'ff'x)
            do 80 j=ije1(i-1)+1,ije1(i)
              ind=igen(istrt1+j)
              ke=iand(ishft(ind,-24),'ff'x)
              kf=iand(ishft(ind,-16),'ff'x)
              kg=iand(ishft(ind,-8),'ff'x)
              kh=iand(ind,'ff'x)
              trle4(j)=-1.d0*isign*
     1        (eps(ka)+eps(kb)+eps(kc)+eps(kd)
     2        -eps(ke)-eps(kf)-eps(kg)-eps(kh))
              trle4(isize(4)+j)=gen(istrt2+j)
80          continue
70        continue
        end if
      end if
      return
      end
      subroutine detlst(nval)
      implicit real*8 (a-h,o-z)
      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      common /det/ nvlspo,lstvso(100),lstdet(1000),ndets,isplit
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
              idum=(ior(ior(ishft(lstvso(i),16),
     1        ishft(lstvso(j),8)),lstvso(k)))
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
                idum=ior(
     1          ior(ishft(lstvso(i),24),ishft(lstvso(j),16)),
     2          ior(ishft(lstvso(k),8),lstvso(l)))
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
                  idum=ior(
     1            ior(ior(ishft(i,16),ishft(j,12)),
     2            ior(ishft(k,8),ishft(l,4))),m)
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
                    idum=ior(ior(ishft(i,20),ishft(j,16)),
     1              ior(ior(ishft(k,12),ishft(l,8)),ior(ishft(m,4),n)))
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
              idum=ior(ishft(lstvso(i),16),
     1        ior(ishft(lstvso(j),8),lstvso(k)))
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
      subroutine dfmax3 (name1,name2,names,tmax,idmax,indmax)
      implicit real*8 (a-h,o-z)
      call add   (name1,-1.d0,name2,1.d0,names)
      call fmax6 (names,tmax,idmax,indmax)
      return
      end
      subroutine dfmax4 (name1,name2,names,tmax,idmax,indmax)
      implicit real*8 (a-h,o-z)
      call add   (name1,-1.d0,name2,1.d0,names)
      call fmax8 (names,tmax,idmax,indmax)
      return
      end
      subroutine dge5n2 (name2,name3,namet,iord)
      implicit real*8 (a-h,o-z)
c  6x6=10
      
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /idv/ iv(200) ,nvalo
      dimension ia(8),save1(1000),save2(1000),iord(10)
      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
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
          ii=iand(ishft(ij1(i),-4*k),'f'x)
          ia(5-k)=iv(ii)
60      continue
        do 20 j=1,ije1(i)-ije1(i-1)
          ii=igen(istrt11+ije1(i-1)+j)
          do 70 k=4,0,-1
            jj=iand(ishft(ii,-4*k),'f'x)
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
              if (ia(iord(5)).eq.iand(ishft(ii,-8),'ff'x)) then
                nel1=nel1+1
                call dimx (nel1,1000,'dge5n2')
                save1(nel1)=gen(istrt22+ije2(k2save-1)+k)
              end if
50          continue 
            nel2=0
            do 80 k=1,ije3(k3save)-ije3(k3save-1)
              ii=igen(istrt31+ije3(k3save-1)+k)
              if (ia(iord(10)).eq.iand(ishft(ii,-8),'ff'x)) then
                nel2=nel2+1
                call dimx (nel2,1000,'dge5n2')
                save2(nel2)=gen(istrt32+ije3(k3save-1)+k)
              end if
80          continue 
            if (nel1.ne.nel2) then  
              write (6,*) 'mismatch in lengths of records - dge5n2'
              stop
            else  
              gen(istrt12+ije1(i-1)+j)=
     1        ddot(nel1,save1,1,save2,1)  
            end if
          end if
20      continue
10    continue
      return
      end
      subroutine dge6n2 (name2,name3,namet,iord)
      implicit real*8 (a-h,o-z)
c  6x8=12
      
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /idv/ iv(200),nvalo
      dimension save1(1000),save2(1000),ia(12),iord(12)
      ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)
      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'dge6n2 input tape does not exist, ',name2
        stop                      
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'dge6n2 input tape does not exist, ',name3
        stop                       
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
          ii=iand(ishft(ij1(i),-4*k),'f'x)
          ia(6-k)=iv(ii)
60      continue
        do 20 j=1,ije1(i)-ije1(i-1)
          ii=igen(istrt11+ije1(i-1)+j)
          do 70 k=5,0,-1
            jj=iand(ishft(ii,-4*k),'f'x)
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
              if (ii.eq.iand(ishft(jj,-8),'ffffff'x)) then
                nel1=nel1+1
                save1(nel1)=gen(istrt22+ije2(k2save-1)+k)
                call dimx (nel1,1000,'dge6n2')
              end if
40          continue
            nel2=0
            ii=ia(iord(12))  
            do 80 k=1,ije3(k3save)-ije3(k3save-1)
              jj=igen(istrt31+ije3(k3save-1)+k)
              if (ii.eq.iand(ishft(jj,-8),'ff'x)) then
                nel2=nel2+1
                save2(nel2)=gen(istrt32+ije3(k3save-1)+k)
                call dimx (nel2,1000,'dge6n2')
              end if
80          continue
            if (nel1.ne.nel2) then
              write (6,*) 'length mismatch in dge6n2'
              stop
            else
              gen(istrt12+ije1(i-1)+j)=ddot(nel1,save1,1,save2,1)
            end if
          end if
20      continue
10    continue
      return
      end
      subroutine diage5 (name2,name3,name4,namet,iord)
      implicit real*8 (a-h,o-z)
c   6x(4+4)=10
      
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /cij4/ nij4,ij4(ndij),ije4(0:ndij),iii4,iiv4
      common /idv/ iv(200),nvalo
      dimension iord(10),ia(10)
      ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)
      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage5 input tape does not exist, ',name2
        stop                       
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage5 input tape does not exist, ',name3
        stop                  
      end if
      call files_ccsd (name4,ntape4,new,4,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage5 input tape does not exist, ',name4
        stop                      
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
          ii=iand(ishft(ij1(i),-4*k),'f'x)
          ia(5-k)=iv(ii)
60      continue
        do 20 j=1,ije1(i)-ije1(i-1)
          ii=igen(istrt11+ije1(i-1)+j)
          do 70 k=4,0,-1
            jj=iand(ishft(ii,-4*k),'f'x)
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
              ii=iand(igen(istrt31+ije3(k3save-1)+k),'ff'x)
              do 100 l=1,ije4(k4save)-ije4(k4save-1)
                jj=iand(igen(istrt41+ije4(k4save-1)+l),'ff'x)
                kk=ior(ishft(ii,8),jj)
                do 110 m=1,ije2(k2save)-ije2(k2save-1)
                  if (kk.eq.   
     1            iand(igen(istrt21+ije2(k2save-1)+m),'ffff'x)) then
                    sum=sum+
     1              gen(istrt22+ije2(k2save-1)+m)*
     2              gen(istrt32+ije3(k3save-1)+k)*
     3              gen(istrt42+ije4(k4save-1)+l)
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
      subroutine diage6 (name2,name3,name4,name5,namet,iord)
      implicit real*8 (a-h,o-z)
c   6x(4+4+4)=12
      
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /cij4/ nij4,ij4(ndij),ije4(0:ndij),iii4,iiv4
      common /cij5/ nij5,ij5(ndij),ije5(0:ndij),iii5,iiv5
      common /idv/ iv(200),nvalo
      dimension iord(12),ia(12)
      ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)
      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage6 input tape does not exist, ',name2
        stop                       
      end if
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage6 input tape does not exist, ',name3
        stop                  
      end if
      call files_ccsd (name4,ntape4,new,4,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage6 input tape does not exist, ',name4
        stop                      
      end if
      call files_ccsd (name5,ntape5,new,5,-1,-1)
      if (new.ne.0) then
        write (6,*) 'diage6 input tape does not exist, ',name5
        stop                      
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
          ii=iand(ishft(ij1(i),-4*k),'f'x)
          ia(6-k)=iv(ii)
60      continue
        do 20 j=1,ije1(i)-ije1(i-1)
          ii=igen(istrt11+ije1(i-1)+j)
          do 70 k=5,0,-1
            jj=iand(ishft(ii,-4*k),'f'x)
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
          if ((igot2.eq.1).and.(igot3.eq.1).and.
     1    (igot4.eq.1).and.(igot5.eq.1)) then
            sum=0.d0
            do 90 k=1,ije3(k3save)-ije3(k3save-1)
              ii=iand(igen(istrt31+ije3(k3save-1)+k),'ff'x)
              do 100 l=1,ije4(k4save)-ije4(k4save-1)
                jj=iand(igen(istrt41+ije4(k4save-1)+l),'ff'x)
                do 130 m=1,ije5(k5save)-ije5(k5save-1)
                  kk=iand(igen(istrt51+ije5(k5save-1)+l),'ff'x)
                  ll=ipack3(ii,jj,kk)
                  do 110 n=1,ije2(k2save)-ije2(k2save-1)
                    if (ll.eq.
     1              iand(igen(istrt21+ije2(k2save-1)+n),'ffffff'x)) then
                      sum=sum+
     1                gen(istrt22+ije2(k2save-1)+n)*
     2                gen(istrt32+ije3(k3save-1)+k)*
     3                gen(istrt42+ije4(k4save-1)+l)*
     4                gen(istrt52+ije5(k5save-1)+m)
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
      subroutine diffmax(name1,name2,names,tmax,itmax)
      implicit real*8 (a-h,o-z)
      call add(name2,1.d0,name1,-1.d0,names)
      call fmax(names,tmax,itmax)
      return
      end
      subroutine diis(iter,icnvg,i12,name1n,name1o,name2n,name2o)
      implicit real*8 (a-h,o-z)
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      dimension r(6,6),rr(8,8),vec(100)
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
c  store current t and delta(t) vector.
c  calc r matrix using all delta(t) vectors.
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
c         write(*,1000) i,(gen(istrt1+j),j=1,nsize)
1000      format('diis vector1',i4/(5f12.6))
c         write(*,*) 'i,id,r=',i,id,r(i,id)
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
          r(i,id)=r(i,id)+ddot(nsize,gen(istrt1+1),1,
     1     gen(istrt2+1),1)
c         write(*,1100) i,(gen(istrt1+j),j=1,nsize)
1100      format('diis vector2',i4/(5f12.6))
          r(id,i)=r(i,id)
c         write(*,*) 'i,id,r=',i,id,r(i,id)
        enddo
        enddo
      endif
c     
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
c     write(*,100) irdim,((rr(j,i),j=1,irdim+1),i=1,irdim+2)
100   format(' irdim, rr=',i4,/(5f12.6))
      call minv(rr,irdim+1,maxdim+2,vec,det,0.d0,1,0)
      do i=1,irdim
      vec(i)=rr(i+1,irdim+2)
      enddo
c     write (*,*) 'diis vec=',(vec(i),i=1,irdim)
c
c    calc diis correction. 
c
      if (i12/10.gt.0) then
        do i=1,irdim
          if (i.eq.1) then
            call cmult(name1di(1),vec(1),"t1d0")
          else
            call add("t1d0",1.d0,name1di(i),vec(i),"t1d0")
          endif
c       write(*,*) 'add dt1(',i,'), vec(i)=',vec(i)
c         call prt(name1di(i),name1di(i),0)
c         call prt("t1d0","t1d0 dt1",0)
        enddo
        call fmax("t1d0",t1dimx,it1dimx)
c       call add(name1n,1.d0,name1o,-1.d0,name1di(idim))
c    calc and stor new t1.
        do i=1,irdim
          call add("t1d0",1.d0,name1(i),vec(i),"t1d0")
c       write(*,*) 'add t1(',i,'), vec(i)=',vec(i)
c         call prt(name1(i),name1(i),0)
c         call prt("t1d0","t1d0  t1",0)
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
c       call add(name2n,1.d0,name2o,-1.d0,name2di(idim))
c    calc and stor new t2.
        do i=1,irdim
          call add("t2d0",1.d0,name2(i),vec(i),"t2d0")
        enddo
        call move("t2d0",name2n)
        call fmax(name2n,t2mx,it2mx)
        write(6,2500) t2mx,idec(it2mx),t1mx,idec(it1mx)
2500    format(' diis: max t, t1=',f15.9,2x,i12,f17.9,2x,i6)
c       call prt(name2n,"t2diis  ",0)
      endif
      return
      end
c ***********************************************
!      subroutine diis(iter,i12,name1n,name1o,name2n,name2o)
! old diis_o.f file
      subroutine diis_notused(iter,i12,name1n,name1o,name2n,name2o)
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
      data name1di /"t101","t102","t103","t104","t105","t106",
     1 "t107","t108","t109","t110","t111","t112","t113","t114",
     2 "t115","t116","t117","t118","t119","t120"/
      data name2di /"td01","td02","td03","td04","td05","td06",
     1 "td07","td08","td09","td10","td11","td12","td13","td14",
     2 "td15","td16","td17","td18","td19","td20"/
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
          r(i,idim)=r(i,idim)+
     1          ddot(nsize,gen(istrt1+1),1,gen(istrt2+1),1)
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
          r(i,idim)=r(i,idim)+ddot(nsize,gen(istrt1+1),1,
     1     gen(istrt2+1),1)
          write(*,1100) i,(gen(istrt1+j),j=1,nsize)
1100      format('diis vector2',i4/(5f12.6))
          r(idim,i)=r(i,idim)
          write(*,*) 'i,idim,r=',i,idim,r(i,idim)
        enddo
      endif
c     
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
c     write(*,100) irdim,((rr(j,i),j=1,irdim+1),i=1,irdim+2)
100   format(' irdim, rr=',i4,/(5f12.6))
      call minv(rr,irdim+1,22,vec,det,0.d0,1,0)
      do i=1,irdim
      vec(i)=rr(i+1,irdim+2)
      enddo
      write (*,*) 'vec ',(vec(i),i=1,irdim)
c
c    calc diis correction. update err vector idim.
c
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
c ***********************************************
      SUBROUTINE MINV(AB,N,LDA,SCRATCH,DET,EPS,M,MODE)
C
C     A subroutine that calculates the determinant and inverse of
C          a matrix, as well as solving systems of linear equations.
C     Martin J. McBride.  11/25/85.
C     General Electric CRD, Information System Operation.
C
C  THE COMMENTS ABOVE REFER TO AN OLD "SCIPORT" VERSION OF THIS ROUTINE,
C  PRESUMABLY WRITTEN BY MR. MC BRIDE.  UNFORTUNATELY, IT SEEMS THAT 
C  HE IS A BAD PROGRAMMER AND THIS ROUTINE DID NOT WORK.  MINV IS NOW
C  A SIMPLE INTERFACE TO LINPACK ROUTINES WHICH *WORK*
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AB(LDA,1),SCRATCH(1),X(2)
C
      IOFF=1+N
      CALL DGECO(AB,LDA,N,SCRATCH,X,SCRATCH(IOFF))
      IF(M.EQ.1)CALL DGESL(AB,LDA,N,SCRATCH,AB(1,N+1),0)
      IF(M.GT.1)THEN
       WRITE(6,1000)
1000   FORMAT(T3,'@MINV-F, Only one system can be solved at a time.')
       CALL abext('minv')
      ENDIF
      IF(M.EQ.0)CALL DGEDI(AB,LDA,N,SCRATCH,X,SCRATCH(IOFF),11)
C
C THIS IS WRONG, BUT WHO CARES?
C
      DET=X(1)      
C
      RETURN
      END
      subroutine dimx(nn,n0,na)
      implicit real*8 (a-h,o-z)

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      character*8 na
      if (nn.le.n0) then   
        return
      else
        write (6,10) na,nn,n0
10      format('0 **** dimension overflow at ',a8,1h:,i8,
     1   ', n allowed ',i8)
        do 20 i=1,nflcur
          write (6,99) i,namefl(i),nstrt(2,i)-nstrt(2,i-1)
99        format (1x,i3,1x,a8,1x,i20)
20      continue
        call abext('dimx    ')
      end if
      return
      end
      subroutine diveps(namef)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      call files_ccsd(namef,ntapef,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) ' file for diveps not found'
        call abext('diveps  ')
      endif
      nsize=ije1(nij1)
      istrt1=nstrt(1,ntapef-1)+2*nij1+4
      istrt2=nstrt(2,ntapef-1)
      if (iii1.gt.100) then
        do 20 i=1,nsize
        if (gen(istrt2+i).ne.0.) then
	  igg=igen(istrt1+i)
	  ka=iand(ishft(igg,-24),'ff'X)
	  kb=iand(ishft(igg,-16),'ff'X)
	  kc=iand(ishft(igg, -8),'ff'X)
	  kd=iand(      igg     ,'ff'X)
          gen(istrt2+i)=gen(istrt2+i)/
     2    (eps(ka)+eps(kb)-eps(kc)-eps(kd))
        endif
20      continue
      else
        do 25 i=1,nsize
        if (gen(istrt2+i).ne.0.d0) then
	  kc=iand(ishft(igen(istrt1+i),-8),'ff'X)
	  kd=iand(      igen(istrt1+i)    ,'ff'X)
          gen(istrt2+i)=gen(istrt2+i)/(eps(kc)-eps(kd))
        endif
25      continue
      endif
      return
      end
      subroutine dveps6 (name)           
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      call files_ccsd (name,ntape,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) 'file for dveps6 not found'
        stop                       
      end if
      istrt1=nstrt(1,ntape-1)+2*nij1+4
      istrt2=nstrt(2,ntape-1)
      ijump=0
      do 10 i=1,nij1
        ig=ij1(i)
        ka=iand(ishft(ig,-16),'ff'x)
        kb=iand(ishft(ig,-8),'ff'x)
        kc=iand(ig,'ff'x)
        do 20 j=ije1(i-1)+1,ije1(i)
          ig=igen(istrt1+j)
          kd=iand(ishft(ig,-16),'ff'x)
          ke=iand(ishft(ig,-8),'ff'x)
          kf=iand(ig,'ff'x)
          if ((icv(ka).eq.1).and.(icv(kb).eq.1).and.
     1    (icv(kc).eq.1).and.(icv(kd).eq.1).and.
     2    (icv(ke).eq.1).and.(icv(kf).eq.1)) then
            if ((ic(ka).eq.3).and.(ic(kb).eq.3).and.
     1      (ic(kc).eq.3).and.(ic(kd).eq.3).and.
     2      (ic(ke).eq.3).and.(ic(kf).eq.3)) then
              ijump=1
            else if ((ic(ka).eq.2).and.(ic(kb).eq.2).and.
     1      (ic(kc).eq.2).and.(ic(kd).eq.2).and.
     2      (ic(ke).eq.2).and.(ic(kf).eq.2)) then
              ijump=1
            else if ((ic(ka).eq.2).and.(ic(kb).eq.2).and.
     1      (ic(kc).eq.3).and.(ic(kd).eq.2).and.
     2      (ic(ke).eq.2).and.(ic(kf).eq.3)) then
              ijump=1
            else if ((ic(ka).eq.2).and.(ic(kb).eq.2).and.
     1      (ic(kc).eq.3).and.(ic(kd).eq.2).and.
     2      (ic(ke).eq.3).and.(ic(kf).eq.2)) then
              ijump=1
            end if
          end if
          if (ijump.eq.1) then
            ijump=0
          else
            gen(istrt2+j)=gen(istrt2+j)/
     1      (eps(ka)+eps(kb)+eps(kc)-eps(kd)-eps(ke)-eps(kf))
          end if
20      continue
10    continue
      return
      end
      subroutine dveps8 (name)           
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      call files_ccsd (name,ntape,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) 'file for dveps8 not found'
        stop                       
      end if
      istrt1=nstrt(1,ntape-1)+2*nij1+4
      istrt2=nstrt(2,ntape-1)
      do 10 i=1,nij1
        id=ij1(i)
        ka=iand(ishft(id,-24),'ff'x)
        kb=iand(ishft(id,-16),'ff'x)
        kc=iand(ishft(id,-8),'ff'x)
        kd=iand(id,'ff'x)
        do 20 j=ije1(i-1)+1,ije1(i)
          ind=igen(istrt1+j)
          ke=iand(ishft(ind,-24),'ff'x)
          kf=iand(ishft(ind,-16),'ff'x)
          kg=iand(ishft(ind,-8),'ff'x)
          kh=iand(ind,'ff'x)
            gen(istrt2+j)=gen(istrt2+j)/
     1      (eps(ka)+eps(kb)+eps(kc)+eps(kd)
     2      -eps(ke)-eps(kf)-eps(kg)-eps(kh))
20      continue
10    continue
      return
      end
      subroutine ehp(ipr)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /csize/ nszv,nszhv,nszhp,nszhq,istrv1,istrv2,
     1 istrhv1,istrhv2,istrhp1,istrhp2,istrhq1,istrhq2
      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      common /rsp2/  nrsp2,irsp2(500)
      dimension ev(60),evi(60)
      dimension scr(200),vec(60,60),ham(60,60)
      common /cdet/ ldet(2,400),lsdet(400)
      dimension lldet(2,400)
      dimension ecol(400),ncol(400),iscol(400),idcol(2,400),ntr(12)
      dimension igr(400)
      data ntr /1hx,1hs,1hp,1hd,1hf,1hg,1hh,1hi,1hk,1hl,1hm,1hn/
c
      icol=0
      nd=0
      ngr=0
c
      do 270 ib=1,norb
      if (ic(ib).ne.3.or.icv(ib).eq.0) go to 270
c
      do 260 iq=1,norb
      if (ic(iq).ne.2.or.icv(iq).eq.0) go to 260
      if (iopt(9).ne.0) then
        i12=ior(ishft(ib,8),iq)
        do 250 i=1,nrsp2
250     if (irsp2(i).eq.i12) go to 260
      endif
      ldet(1,nd+1)=ib
      ldet(2,nd+1)=iq
      ldet(1,nd+2)=ib
      ldet(2,nd+2)=-iq
      ldet(1,nd+3)=-ib
      ldet(2,nd+3)=-iq
      nd=nd+3
260   continue
270   continue
      ndet=nd
      if (ipr.ne.0) print 2900, (i,(ldet(j,i),j=1,2),i=1,ndet)
2900  format (4('   det',i4,': spinorbs',2i3,1h.))
      call dimx(ndet,400,"ndet")
c
c  read folded files.
      call files_ccsd("tvf ",ntape1,new,1,-1,-1)
      nszv=ije1(nij1)
      istrv1=nstrt(1,ntape1-1)+2*nij1+4
      istrv2=nstrt(2,ntape1-1)
      call files_ccsd("tgf ",ntape2,new,1,-1,-1)
      nszhv=ije1(nij1)
      istrhv1=nstrt(1,ntape2-1)+2*nij1+4
      istrhv2=nstrt(2,ntape2-1)
      call files_ccsd("thpf",ntape3,new,1,-1,-1)
      nszhp=ije1(nij1)
      istrhp1=nstrt(1,ntape3-1)+2*nij1+4
      istrhp2=nstrt(2,ntape3-1)
      call files_ccsd("thqf",ntape4,new,1,-1,-1)
      nszhq=ije1(nij1)
      istrhq1=nstrt(1,ntape4-1)+2*nij1+4
      istrhq2=nstrt(2,ntape4-1)
c
c     lsdet(i) is -1 if det(i) has not been treated yet;
c                  0 if det(i) has already been treated;
c                  n if det(i) is now being treated, with ordinal no n.
      do 310 i=1,ndet
310   lsdet(i)=-1
320   ndt=0
      ngr=ngr+1
      do 330 i=1,60
      do 330 j=1,60
330   ham(i,j)=0.
      do 370 id1=1,ndet
      if (lsdet(id1).eq.0) go to 370
      if (ndt.ne.0.and.lsdet(id1).eq.-1) go to 370
      if (ndt.eq.0) isp=iabs(isign(1,ldet(1,id1))-isign(1,ldet(2,id1)))
     1 +1
      do 360 id2=1,ndet
      if (lsdet(id2).eq.0) go to 360
      hh=hdet11(ldet,id1,id2)
      if (abs(hh).lt.1.d-13) go to 360
      if (lsdet(id2).eq.-1) then
        ndt=ndt+1
        lsdet(id2)=ndt
        lldet(1,ndt)=ldet(1,id2)
        lldet(2,ndt)=ldet(2,id2)
      endif
      ham(lsdet(id1),lsdet(id2))=hh
360   continue
370   continue
      if (ndt.eq.0) go to 870
c
c  some matrix elements may have been overlooked. correct.
      do 390 id1=1,ndet
      if (lsdet(id1).le.0) go to 390
      do 380 id2=1,ndet
      if (lsdet(id2).le.0) go to 380
      if (ham(lsdet(id1),lsdet(id2)).eq.0.)
     1  ham(lsdet(id1),lsdet(id2))=hdet11(ldet,id1,id2)
380   continue
390   continue
      if (ipr.gt.0) write (*,3900) ngr,ndt,(i,lldet(1,i),lldet(2,i),
     1 i=1,ndt)
3900  format(' group',i3,'.',i4,' dets included:',8(i4,':',2i3)/
     1 (29x,6(i4,':',2i3)))
      call dimx(ndt,60,"ndt ")
c
      if (ipr.le.1) go to 830
      print 8100
8100  format ('  effective hamiltonian:')
      do 820 i=1,ndt
820   print 8200, i,(ham(i,j),j=1,ndt)
8200  format (' row',i3,10f12.7/(1x,10f12.7))
830   call rg(60,ndt,ham,ev,evi,1,vec,scr,scr(61),ierr)
      if (ierr.ne.0) then
	write (*,'('' rg ierr='',i4)') ierr
        call abext('ehp     ')
      endif
      if (ipr.eq.0) print 8300,(ev(i),i=1,ndt)
8300  format (' energies:',8e15.8/(1x,8e15.8))
      if (ipr.eq.0) go to 855
c   order and print eigenvalues, vectors.
      do 850 ii=1,ndt
      e=10000.
      do 845 i1=1,ndt
      if (ev(i1).ge.e) go to 845
      e=ev(i1)
      i=i1
845   continue
      print 8400, ii,ev(i),(vec(j,i),j=1,ndt)
8400  format (' energy',i3,1h:,e20.8,'. vector:',7f12.7/
     1 (40x,7f12.7))
c
c     find main determinant contributing
      cc=0
      jj=0
      do 846 j=1,ndt
      if (abs(vec(j,i)).le.cc) go to 846
      cc=abs(vec(j,i))
      jj=j
846   continue
      do 8462 id=1,ndet
      if (lsdet(id).ne.jj) go to 8462
      idd=id
      go to 8464
8462  continue
8464  continue
c
c     collect equal energy terms.
      if (icol.eq.0) go to 848
      do 847 i1=1,icol
      if (abs(ev(i)-ecol(i1)).gt.1.d-6) go to 847
      ncol(i1)=ncol(i1)+1
      if (isp.gt.1) ncol(i1)=ncol(i1)+1
      if (isp.le.iscol(i1)) go to 849
      iscol(i1)=isp
      idcol(1,i1)=ldet(1,idd)
      idcol(2,i1)=ldet(2,idd)
      go to 849
847   continue
848   icol=icol+1
      iscol(icol)=isp
      idcol(1,icol)=ldet(1,idd)
      idcol(2,icol)=ldet(2,idd)
      ecol(icol)=ev(i)
      ncol(icol)=1
      igr(icol)=ngr
849   ev(i)=10000.
850   if (evi(i).ne.0.d0) print 8500, evi(i),(vec(j,i+1),j=1,ndt)
8500  format (1x,10(1h*),' im.energy:',e20.8,'. im. vec:',5f12.7/
     1 (1x,10f12.7))
855   do 860 id=1,ndet
860   if (lsdet(id).gt.0) lsdet(id)=0
      go to 320
c     order and print collected energies.
870   if (ipr.eq.0) go to 890
      do 880 ii=1,icol
      e=10000.
      do 875 i1=1,icol
      if (ecol(i1).ge.e) go to 875
      e=ecol(i1)
      i=i1
875   continue
      itr=(ncol(i)/iscol(i)-1)/2
      if ((2*itr+1)*iscol(i).ne.ncol(i)) itr=-1
      exc=ecol(i)*27.21
      excm=exc*8066.
      print 8700, ii,ecol(i),ncol(i),exc,excm,iscol(i),ntr(itr+2),
     1 igr(i),idcol(1,i),idcol(2,i)
8700  format (' energy',i3,e16.8,'. degen',i3,'. exc en=',f7.3,
     1 'ev,',f8.0,'cm-1. term ',i1,a1,'. gr',i3,', major det=',2i3)
880   ecol(i)=10000.
890   continue
      time=cputime(0.)
      return
      end
      subroutine eneru2 (egr0,nval,ihp,ipr)
      implicit real*8 (a-h,o-z)
      external ispin,hmatel
      
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
     1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
     2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /cij4/ nij4,ij4(ndij),ije4(0:ndij),iii4,iiv4
      common /cij5/ nij5,ij5(ndij),ije5(0:ndij),iii5,iiv5
      common /cij6/ nij6,ij6(ndij),ije6(0:ndij),iii6,iiv6
      common /dat/  norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      common /det/  nvlspo,lstvso(100),lstdet(1000),ndets,isplit 
      common /idv/  iv(200),nvalo
      dimension ecol(200),eicol(200),ncol(200),iscol(200),idcol(2,600),
     1 ntr(12)
      dimension ev(60),evi(60)
      dimension scr(200),vec(60,60),veci(60,60),ham(60,60)
      dimension lldet(6,700),lsdet(700),list(6)
      equivalence (veci,ham)
      data ntr /1hx,1hs,1hp,1hd,1hf,1hg,1hh,1hi,1hk,1hl,1hm,1hn/
      if (nval.gt.6) then
        nvalp=3
        egr=egr0
      else
        nvalp=nval
      end if
      if (nval.le.4) then
        nvlspo=0
        do 10 i=1,norb
          if (icv(i).eq.1) then
            if (ic(i).eq.ihp) then
              nvlspo=nvlspo+2
              lstvso(nvlspo-1)=2*i-1
              lstvso(nvlspo)=2*i
            end if
          end if
10      continue
      else if ((nval.eq.5).or.(nval.eq.6)) then     
        nvlspo=2*nvalo
      else if (nval.eq.12) then
        nvlspo=0
        do 160 i=1,norb
          if ((icv(i).eq.1).and.(ic(i).eq.2)) then
            nvlspo=nvlspo+2
            lstvso(nvlspo-1)=2*i-1
            lstvso(nvlspo)=2*i
          end if
160     continue
        isplit=nvlspo
        do 170 i=1,norb
          if ((icv(i).eq.1).and.(ic(i).eq.3)) then
            nvlspo=nvlspo+2
            lstvso(nvlspo-1)=2*i-1
            lstvso(nvlspo)=2*i
          end if
170     continue
      else if (nval.eq.21) then
        nvlspo=0
        do 180 i=1,norb
          if ((icv(i).eq.1).and.(ic(i).eq.3)) then
            nvlspo=nvlspo+2
            lstvso(nvlspo-1)=2*i-1
            lstvso(nvlspo)=2*i
          end if
180     continue
        isplit=nvlspo
        do 190 i=1,norb
          if ((icv(i).eq.1).and.(ic(i).eq.2)) then
            nvlspo=nvlspo+2
            lstvso(nvlspo-1)=2*i-1
            lstvso(nvlspo)=2*i
          end if
190     continue
      end if
      call dimx (nvlspo,100,"nvlspo  ","eneru   ",0)
      if (ihp.eq.3) then
        if (nval.ne.21) then
          call files_ccsd ("tvf ",ntape1,new,1,-1,-1)
          if (new.ne.0) then
            write (6,*) 'eneru input file not found - tvf'
            call abext('eneru   ')
          end if
          istrt11=nstrt(1,ntape1-1)+2*nij1+4
          istrt12=nstrt(2,ntape1-1)
          nsize1=ije1(nij1)
          if (nval.ge.2) then
            call files_ccsd ("tuf ",ntape2,new,2,-1,-1)
            if (new.ne.0) then
              write (6,*) 'eneru input file not found - tuf'
              call abext('eneru   ')
            end if
            istrt21=nstrt(1,ntape2-1)+2*nij2+4
            istrt22=nstrt(2,ntape2-1)
            nsize2=ije2(nij2)
          end if
          if (nval.ge.3) then
            call files_ccsd ("twf ",ntape3,new,3,-1,-1)
            if (new.ne.0) then
              write (6,*) 'eneru input file not found - twf'
              call abext('eneru   ')
            end if
            istrt31=nstrt(1,ntape3-1)+2*nij3+4
            istrt32=nstrt(2,ntape3-1)
            nsize3=ije3(nij3)
          end if
          if (nval.ge.4) then
            call files_ccsd ("txf ",ntape4,new,4,-1,-1)
            if (new.ne.0) then
              write (6,*) 'eneru input file not found - txf'
              call abext('eneru   ')
            end if
            istrt41=nstrt(1,ntape4-1)+2*nij4+4
            istrt42=nstrt(2,ntape4-1)
            nsize4=ije4(nij4)
          end if
          if (nval.ge.5) then
            call files_ccsd ("tyf ",ntape5,new,5,-1,-1)
            if (new.ne.0) then
              write (6,*) 'eneru input file not found - tyf'
              call abext('eneru   ')
            end if
            istrt51=nstrt(1,ntape5-1)+2*nij5+4
            istrt52=nstrt(2,ntape5-1)
            nsize5=ije5(nij5)
          end if
          if (iopt(22).eq.1) then
            if (nval.eq.6) then
              call files_ccsd ("tzf ",ntape6,new,6,-1,-1)
              if (new.ne.0) then
                write (6,*) 'eneru input file not found - tzf'
                call abext('eneru   ')
              end if
              istrt61=nstrt(1,ntape6-1)+2*nij6+4
              istrt62=nstrt(2,ntape6-1)
              nsize6=ije6(nij6)
            end if
          end if
        end if
      else if (ihp.eq.2) then
        if (nval.ne.12) then
          call files_ccsd ("thvf",ntape1,new,1,-1,-1)
          if (new.ne.0) then
            write (6,*) 'eneru input file not found - thvf'
            call abext('eneru   ')
          end if
          istrt11=nstrt(1,ntape1-1)+2*nij1+4
          istrt12=nstrt(2,ntape1-1)
          nsize1=ije1(nij1)
          if (nval.ge.2) then
            call files_ccsd ("thuf",ntape2,new,2,-1,-1)
            if (new.ne.0) then
              write (6,*) 'eneru input file not found - thuf'
              call abext('eneru   ')
            end if
            istrt21=nstrt(1,ntape2-1)+2*nij2+4
            istrt22=nstrt(2,ntape2-1)
            nsize2=ije2(nij2)
          end if
          if (nval.ge.3) then
            call files_ccsd ("thwf",ntape3,new,3,-1,-1)
            if (new.ne.0) then
              write (6,*) 'eneru input file not found - thwf'
              call abext('eneru   ')
            end if
            istrt31=nstrt(1,ntape3-1)+2*nij3+4
            istrt32=nstrt(2,ntape3-1)
            nsize3=ije3(nij3)
          end if
          if (nval.ge.4) then
            call files_ccsd ("thxf",ntape4,new,4,-1,-1)
            if (new.ne.0) then
              write (6,*) 'eneru input file not found - thxf'
              call abext('eneru   ')
            end if
            istrt41=nstrt(1,ntape4-1)+2*nij4+4
            istrt42=nstrt(2,ntape4-1)
            nsize4=ije4(nij4)
          end if
          if (nval.ge.5) then
            call files_ccsd ("thyf",ntape5,new,5,-1,-1)
            if (new.ne.0) then
              write (6,*) 'eneru input file not found - thyf'
              call abext('eneru   ')
            end if
            istrt51=nstrt(1,ntape5-1)+2*nij5+4
            istrt52=nstrt(2,ntape5-1)
            nsize5=ije5(nij5)
          end if
          if (iopt(22).eq.1) then
            if (nval.eq.6) then
              call files_ccsd ("thzf",ntape6,new,6,-1,-1)
              if (new.ne.0) then
                write (6,*) 'eneru input file not found - thzf'
                call abext('eneru   ')
              end if
              istrt61=nstrt(1,ntape6-1)+2*nij6+4
              istrt62=nstrt(2,ntape6-1)
              nsize6=ije6(nij6)
            end if
          end if
        end if
      end if
c
c ... start energy calculations.
c
      if (ipr.gt.0) then
        if (nval.lt.12) then
          write (6,2300) nval
2300  format(1h0,80(1h-)/30x,i1,' valence electrons'/30x,19(1h-))
        else
          write (6,2301) nval
2301  format(1h0,80(1h-)/22x,i2,
     1' sector excitation energies'/30x,22(1h-))
        end if
      end if
      icol=0
      call detlst (nval)
      if (ndets.eq.0) then
        write (6,*) 'no determinants found'
        return
      end if
      do 240 i=1,ndets
240   lsdet(i)=-1
      if (ipr.gt.0) then  
        write (6,*) '0 determinant list:'
        if ((nval.eq.5).or.(nval.eq.6)) then
          do 1011 i=1,ndets
            jj=lstdet(i)
            do 1001 j=0,nval-1
              list(nval-j)=ifndso(mod(jj,16))
              jj=jj/16
1001        continue
            write (6,2400) i,(list(j),j=1,nvalp)
1011        continue
          else
            do 1010 i=1,ndets
            jj=lstdet(i)
            do 1000 j=0,nvalp-1
              list(nvalp-j)=mod(jj,256)
              jj=jj/256
1000        continue
            write (6,2400) i,(list(j),j=1,nvalp)
2400        format (1x,i3,':',1x,6i4)
1010      continue
        end if
      end if
      call dimx(ndets,700,"ndets   ","ener",240)
c
c     lsdet is -1 for a determinant not yet treated.
c               n for the nth determinant currently treated.
c               0 for a deteminant already treated.
320   ndt=0
      do 330 i=1,60
      do 330 j=1,60
330   ham(i,j)=0.d0
      do 370 id1=1,ndets
      if (lsdet(id1).eq.0) go to 370
      if (ndt.ne.0.and.lsdet(id1).eq.-1) go to 370
      if (ndt.eq.0) ksp=ispin(lstdet(id1),nval)+1
      do 360 id2=id1,ndets
      if (lsdet(id2).eq.0) go to 360
      hh1=hmatel(lstdet(id1),lstdet(id2),nval,ihp)
      hh2=hmatel(lstdet(id2),lstdet(id1),nval,ihp)
      if (abs(hh1)+abs(hh2).lt.1.d-13) go to 360
      if (lsdet(id2).eq.-1) then
        ndt=ndt+1
        lsdet(id2)=ndt
        if ((nval.eq.5).or.(nval.eq.6)) then
          jj=lstdet(id2)
          do 30 ii=0,nval-1
            lldet(nval-ii,ndt)=ifndso(mod(jj,16))
            jj=jj/16
30        continue
        else
          jj=lstdet(id2)
          do 20 ii=0,nvalp-1
            lldet(nvalp-ii,ndt)=mod(jj,256)
            jj=jj/256
20        continue
        end if
      end if
      ham(lsdet(id1),lsdet(id2))=hh1
      ham(lsdet(id2),lsdet(id1))=hh2
360   continue
370   continue
c  some matrix elements may have been overlooked. correct.
      do 390 id1=1,ndets
      if (lsdet(id1).le.0) go to 390
      do 380 id2=1,ndets
      if (lsdet(id2).le.0) go to 380
      if (ham(lsdet(id1),lsdet(id2)).eq.0.d0) then
        ham(lsdet(id1),lsdet(id2))=
     1  hmatel(lstdet(id1),lstdet(id2),nval,ihp)
      end if
380   continue
390   continue
      if (ipr.gt.0) then 
        write (6,3900) ndt
3900    format (1x,i4,' dets included:')
        do 3901 iiii=1,ndt
          write (6,3902) iiii,(lldet(j,iiii),j=1,nvalp)
3902      format (1x,i3,':',1x,6(1x,i3))
3901    continue
      end if 
      call dimx(ndt,60,"ndt ","ener",350)
c     if (ipr.lt.2) go to 830
      write (*,8100)
8100  format ('  effective hamiltonian:')
      do 820 i=1,ndt
820   write (*,8200) i,(ham(i,j),j=1,ndt)
8200  format (' row',i3,10f12.7/(1x,10f12.7))
830   call rg(60,ndt,ham,ev,evi,1,vec,scr,scr(61),ierr)
      if (ierr.ne.0) then                
        write (*,'('' rg ierr='',i4)') ierr
	call abext('eneru2  ')  
      endif
c      if (ipr.eq.0) write (*,8300) tit,(ev(i),i=1,ndt)
c8300  format (1x,a8,' energies:',7e15.8/(1x,8e15.8))
      if (ipr.eq.0) go to 855
c   order and write  eigenvalues, vectors.
      do 850 ii=1,ndt
      e=10000.d0
      do 845 i1=1,ndt
      if (ev(i1).ge.e) go to 845
      e=ev(i1)
      i=i1
845   continue
c   normalize
      fnor=0.
      do 8452 j=1,ndt
8452  fnor=fnor+vec(j,i)**2
      fnor=1.d0/sqrt(fnor)
      do 8455 j=1,ndt
8455  vec(j,i)=vec(j,i)*fnor
      if (ipr.ne.0) write (*,8400) ii,ev(i),(vec(j,i),j=1,ndt)
8400  format (' energy',i3,1h:,e20.8,'. vector:',7f12.7/
     1 (1x,10f12.7))
c
c     find main determinant contributing
      cc=0
      jj=0
      do 846 j=1,ndt
      if (abs(vec(j,i)).le.cc) go to 846
      cc=abs(vec(j,i))
      jj=j
846   continue
      do 8462 id=1,ndets
      if (lsdet(id).ne.jj) go to 8462
      idd=id
      go to 8464
8462  continue
8464  continue
c
c     collect equal energy terms.
      if (icol.eq.0) go to 848
      do 847 i1=1,icol
      if (abs(ev(i)-ecol(i1)).gt.1.d-7) go to 847
      ncol(i1)=ncol(i1)+1
      if (ksp.gt.1) ncol(i1)=ncol(i1)+1
      iscol(i1)=max0(iscol(i1),ksp)
      iddd=lstdet(idd)
      do 3903 iiii=0,nvalp-1
        if ((nval.eq.5).or.(nval.eq.6)) then
          idcol(nval-iiii,i1)=ifndso(mod(iddd,16))
          iddd=iddd/16
        else
          idcol(nvalp-iiii,i1)=mod(iddd,256)
          iddd=iddd/256
        end if
3903  continue
      go to 849
847   continue
848   icol=icol+1
      if (icol.eq.200) write (*,*) ' ******** 200 different energies.',
     1 ' others not collected. ***********'
      if (icol.gt.200) go to 849
      iscol(icol)=ksp
      iddd=lstdet(idd)
      do 3904 iiii=0,nvalp-1
        if ((nval.eq.5).or.(nval.eq.6)) then
          idcol(nval-iiii,icol)=ifndso(mod(iddd,16))
          iddd=iddd/16
        else
          idcol(nvalp-iiii,icol)=mod(iddd,256)
          iddd=iddd/256
        end if
3904  continue
      ecol(icol)=ev(i)
      eicol(icol)=evi(i)
      ncol(icol)=1
      if (ksp.gt.1) ncol(icol)=2
849   ev(i)=10000.d0
850   if (ipr.ne.0.and.evi(i).ne.0.) write (*,8500) evi(i),
     1 (veci(j,i),j=1,ndt)
8500  format (1x,10(1h*),' im.energy:',e20.8,'. im. vec:',5f12.7/
     1 (1x,10f12.7))
855   iend=1
      do 860 id=1,ndets
      if (lsdet(id).eq.-1) iend=0
860   if (lsdet(id).gt.0) lsdet(id)=0
      if (iend.eq.0) go to 320
c     order and write collected energies.
      if (icol.gt.200) icol=200
      do 880 ii=1,icol
      e=10000.d0
      do 875 i1=1,icol
      if (ecol(i1).ge.e) go to 875
      e=ecol(i1)
      i=i1
875   continue
      itr=(ncol(i)/iscol(i)-1)/2
      if ((2*itr+1)*iscol(i).ne.ncol(i)) itr=-1
      if (ii.ne.1) go to 877
      if (nval.gt.6) goto 877
      egr=ecol(i)
      exc=(egr0-egr)*27.21d0
      egr0=egr
      excm=exc*8066.d0
      write (*,8600) ii,ecol(i),ncol(i),exc,excm,iscol(i),ntr(itr+2),
     1 idcol(1,i),idcol(2,i)
8600  format(' energy',i3,e16.8,'. degen',i3,'. ion pot=',f8.3,
     1 'ev,',f9.0,'cm-1. term ',i1,a1,'. major det',6i3)
      if (eicol(i).ne.0.d0) write (*,8650) eicol(i)
8650  format (1x,10(1h*),' im. energy:',e20.8)
      go to 880
877   exc=(ecol(i)-egr)*27.21d0
      excm=exc*8066.d0
      write (*,8700) ii,ecol(i),ncol(i),exc,excm,iscol(i),ntr(itr+2),
     1 idcol(1,i),idcol(2,i)
8700  format (' energy',i3,e16.8,'. degen',i3,'. exc en=',f7.3,
     1 'ev,',f8.0,'cm-1. term ',i1,a1,'. major det',6i3)
      if (eicol(i).ne.0.d0) write (*,8650) eicol(i)
880   ecol(i)=10000.d0
890   continue
      return
      end
      subroutine expnd(iwork,iwork1,work,work1,nw,nw1)
      implicit real*8 (a-h,o-z)
c   move from work1 (indx iwork1) to work, filling out according to
c   indx iwork
      dimension iwork(nw),iwork1(nw1),work(nw),work1(nw1)
      if (nw1.eq.nw) then
        call dcopy(nw,work1,1,work,1)
        return
      endif
      iw1=nw1
      iw=nw
10    if (iwork1(iw1).eq.iwork(iw)) then
        work(iw)=work1(iw1)
        iw1=iw1-1
        iw=iw-1
      elseif (iwork1(iw1).lt.iwork(iw)) then
        work(iw)=0.d0
        iw=iw-1
      else
        write (*,100) iwork1(iw1)
        write (6,*) 'aborted from expnd'
        call abext('expnd   ')
      endif
      if (iw.ne.0) go to 10
100   format (' unmatched label ',z8)
      return
      end
      subroutine fmax6 (namef,tmax,idmax,indmax)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"
#include "cfiles.h"
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      call files_ccsd(namef,ntapef,new,1,-1,-1)
      if (new.ne.0) then
        write(*,*) ' file for fmax6 not found'
        stop
      endif
      tmax=0.d0
      idmax=0
      indmax=0
      istrt1=nstrt(1,ntapef-1)+2*nij1+4
      istrt2=nstrt(2,ntapef-1)
      do 10 i=1,nij1   
        do 20 j=1,ije1(i)-ije1(i-1)
          if (abs(gen(istrt2+ije1(i-1)+j)).gt.abs(tmax)) then
            tmax=gen(istrt2+ije1(i-1)+j)
            idmax=ij1(i)
            indmax=iand(igen(istrt1+ije1(i-1)+j),'ffffff'x)
          endif
20      continue
10    continue
      return
      end
      subroutine fmax8 (namef,tmax,idmax,indmax)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common // gen(nbig)
      common /i/ igen(nbig)
      call files_ccsd(namef,ntapef,new,1,-1,-1)
      if (new.ne.0) then
        write(*,*) ' file for fmax8 not found'
        stop
      endif
      tmax=0.d0
      idmax=0
      indmax=0
      istrt1=nstrt(1,ntapef-1)+2*nij1+4
      istrt2=nstrt(2,ntapef-1)
      do 10 i=1,nij1   
        do 20 j=1,ije1(i)-ije1(i-1)
          if (abs(gen(istrt2+ije1(i-1)+j)).gt.abs(tmax)) then
            tmax=gen(istrt2+ije1(i-1)+j)
            idmax=ij1(i)
            indmax=igen(istrt1+ije1(i-1)+j)
          endif
20      continue
10    continue
      return
      end
      subroutine fmax(namef,tmax,itmax)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      call files_ccsd(namef,ntapef,new,1,-1,-1)
      if (new.ne.0) then
        write(*,*) ' file for fmax not found'
        call abext('fmax    ')
      endif
      tmax=0.d0
      itmax=0
      nsize=ije1(nij1)
      istrt1=nstrt(1,ntapef-1)+2*nij1+4
      istrt2=nstrt(2,ntapef-1)
      do 20 i=1,nsize
      if (abs(gen(istrt2+i)).gt.abs(tmax)) then
        tmax=gen(istrt2+i)
        itmax=igen(istrt1+i)
      endif
20    continue
      return
      end
      function hdet11(ldet,ia,ib)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /csize/ nszv,nszhv,nszhp,nszhq,istrv1,istrv2,
     1 istrhv1,istrhv2,istrhp1,istrhp2,istrhq1,istrhq2
      dimension ldet(2,*)
      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
c  calculate heff matrix elements
      hdet11=0.
      if (isign(1,ldet(1,ia))-isign(1,ldet(2,ia)).ne.
     1 isign(1,ldet(1,ib))-isign(1,ldet(2,ib))) return
      if (isign(1,ldet(1,ia)).eq.isign(1,ldet(2,ia))) then
        ii=ior(iabs(ldet(1,ia))*256**3,ior(iabs(ldet(2,ib))*256**2,
     1   ior(iabs(ldet(2,ia))*256,iabs(ldet(1,ib)))))
        hdet11=hdet11+find(ii,gen(istrhq2+1),igen(istrhq1+1),nszhq)
      endif
      if (isign(1,ldet(1,ia)).eq.isign(1,ldet(1,ib))) then
        ii=ior(iabs(ldet(1,ia))*256**3,ior(iabs(ldet(2,ib))*256**2,
     1   ior(iabs(ldet(1,ib))*256,iabs(ldet(2,ia)))))
        hdet11=hdet11-find(ii,gen(istrhp2+1),igen(istrhp1+1),nszhp)
      endif
      if (ldet(1,ia).eq.ldet(1,ib)) then
        ii=ior(iabs(ldet(2,ia))*256,iabs(ldet(2,ib)))
        hdet11=hdet11-find(ii,gen(istrhv2+1),igen(istrhv1+1),nszhv)
        if (ldet(2,ia).eq.ldet(2,ib)) then
          hdet11=hdet11-eps(iabs(ldet(2,ia)))
        end if
      endif
      if (ldet(2,ia).eq.ldet(2,ib)) then
        ii=ior(iabs(ldet(1,ia))*256,iabs(ldet(1,ib)))
        hdet11=hdet11+find(ii,gen(istrv2+1),igen(istrv1+1),nszv)
        if (ldet(1,ia).eq.ldet(1,ib)) then
          hdet11=hdet11+eps(iabs(ldet(1,ia)))
         end if
      endif
      return
      end
      function hmatel (id1,id2,nval,ihp)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
     1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
     2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      dimension k1(6),k2(6),k1s(6),k2s(6),
     1k1save(6),k2save(6),j1(6),j2(6)
      if (nval.lt.12) then
        hmatel=0.d0
        if (nval.eq.1) then
          isign=1
          if (ispin(id1,1).ne.ispin(id2,1)) return
          hmatel=v1det(id1,id2)       
        else if (nval.eq.2) then
          if (ispin(id1,2).ne.ispin(id2,2)) return
          do 20 i=0,1 
            k1(2-i)=iand(ishft(id1,-8*i),'ff'x)
            k2(2-i)=iand(ishft(id2,-8*i),'ff'x)
20        continue
          call match6 (k1,k2,2,isign,nm)
          hmatel=v2det(k1(1),k1(2),k2(1),k2(2))                  
          if (nm.ge.1) then
            hmatel=hmatel+v1det(k1(1),k2(1))             
          end if
          if (nm.eq.2) then
            hmatel=hmatel+v1det(k1(2),k2(2))
          end if
        else if (nval.eq.3) then
          if (ispin(id1,3).ne.ispin(id2,3)) return
          do 10 i=0,2 
            k1(3-i)=iand(ishft(id1,-8*i),'ff'x)
            k2(3-i)=iand(ishft(id2,-8*i),'ff'x)
10        continue
          call match6 (k1,k2,3,isign,nm)
          hmatel=v3det(k1(1),k1(2),k1(3),k2(1),k2(2),k2(3))
          if (nm.ge.1) then
            hmatel=hmatel+v2det(k1(1),k1(2),k2(1),k2(2))
          end if
          if (nm.ge.2) then
            hmatel=hmatel+v1det(k1(1),k2(1))
            hmatel=hmatel+v2det(k1(1),k1(3),k2(1),k2(3))
          end if
          if (nm.eq.3) then
            hmatel=hmatel+v1det(k1(2),k2(2))
            hmatel=hmatel+v1det(k1(3),k2(3))
            hmatel=hmatel+v2det(k1(2),k1(3),k2(2),k2(3))
          end if
        else if (nval.eq.4) then
         if (ispin(id1,4).ne.ispin(id2,4)) return
          do 40 i=0,3 
            k1(4-i)=iand(ishft(id1,-8*i),'ff'x)
            k2(4-i)=iand(ishft(id2,-8*i),'ff'x)
40        continue
          call match6 (k1,k2,4,isign,nm)
          hmatel=v4det(k1(1),k1(2),k1(3),k1(4),k2(1),k2(2),k2(3),k2(4))
          if (nm.ge.1) then
            hmatel=hmatel+v3det(k1(1),k1(2),k1(3),k2(1),k2(2),k2(3))
          end if
          if (nm.ge.2) then
            hmatel=hmatel+v3det(k1(1),k1(2),k1(4),k2(1),k2(2),k2(4))
            hmatel=hmatel+v2det(k1(1),k1(2),k2(1),k2(2))
          end if
          if (nm.ge.3) then
            hmatel=hmatel+v3det(k1(1),k1(3),k1(4),k2(1),k2(3),k2(4))
            hmatel=hmatel+v2det(k1(1),k1(3),k2(1),k2(3))
            hmatel=hmatel+v2det(k1(1),k1(4),k2(1),k2(4))
            hmatel=hmatel+v1det(k1(1),k2(1))
          end if
          if (nm.eq.4) then
            hmatel=hmatel+v3det(k1(2),k1(3),k1(4),k2(2),k2(3),k2(4))
            hmatel=hmatel+v2det(k1(2),k1(3),k2(2),k2(3))
            hmatel=hmatel+v2det(k1(2),k1(4),k2(2),k2(4))
            hmatel=hmatel+v2det(k1(3),k1(4),k2(3),k2(4))
            do 50 i=2,4
              hmatel=hmatel+v1det(k1(i),k2(i))
50          continue
          end if
        else if (nval.eq.5) then
          if (ispin(id1,5).ne.ispin(id2,5)) return
          do 70 i=0,4 
            ii=iand(ishft(id1,-4*i),'f'x)
            k1(5-i)=ifndso(ii)
            ii=iand(ishft(id2,-4*i),'f'x)
            k2(5-i)=ifndso(ii)
70        continue
          call match6 (k1,k2,5,isign,nm)
          hmatel=v5det(k1(1),k1(2),k1(3),k1(4),k1(5),
     1    k2(1),k2(2),k2(3),k2(4),k2(5))
          if (nm.ge.1) then
            hmatel=hmatel+v4det(k1(1),k1(2),k1(3),k1(4),
     1                          k2(1),k2(2),k2(3),k2(4))
          end if
          if (nm.ge.2) then
            hmatel=hmatel+v4det(k1(1),k1(2),k1(3),k1(5),
     1                          k2(1),k2(2),k2(3),k2(5))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(3),k2(1),k2(2),k2(3))
          end if
          if (nm.ge.3) then
            hmatel=hmatel+v4det(k1(1),k1(2),k1(4),k1(5),
     1                          k2(1),k2(2),k2(4),k2(5))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(4),k2(1),k2(2),k2(4))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(5),k2(1),k2(2),k2(5))
            hmatel=hmatel+v2det(k1(1),k1(2),k2(1),k2(2))
          end if
          if (nm.ge.4) then
            hmatel=hmatel+v4det(k1(1),k1(3),k1(4),k1(5),
     1                          k2(1),k2(3),k2(4),k2(5))
            hmatel=hmatel+v3det(k1(1),k1(3),k1(4),k2(1),k2(3),k2(4))
            hmatel=hmatel+v3det(k1(1),k1(3),k1(5),k2(1),k2(3),k2(5))
            hmatel=hmatel+v3det(k1(1),k1(4),k1(5),k2(1),k2(4),k2(5))
            do 90 i=2,4
              hmatel=hmatel+v2det(k1(1),k1(1+i),k2(1),k2(1+i))
90          continue
            hmatel=hmatel+v1det(k1(1),k2(1))
          end if
          if (nm.eq.5) then
            hmatel=hmatel+v4det(k1(2),k1(3),k1(4),k1(5),
     1                          k2(2),k2(3),k2(4),k2(5))
            hmatel=hmatel+v3det(k1(2),k1(3),k1(4),k2(2),k2(3),k2(4))
            hmatel=hmatel+v3det(k1(2),k1(3),k1(5),k2(2),k2(3),k2(5))
            hmatel=hmatel+v3det(k1(2),k1(4),k1(5),k2(2),k2(4),k2(5))
            hmatel=hmatel+v3det(k1(3),k1(4),k1(5),k2(3),k2(4),k2(5))
            hmatel=hmatel+v2det(k1(2),k1(3),k2(2),k2(3))
            hmatel=hmatel+v2det(k1(2),k1(4),k2(2),k2(4))
            hmatel=hmatel+v2det(k1(2),k1(5),k2(2),k2(5))
            hmatel=hmatel+v2det(k1(3),k1(4),k2(3),k2(4))
            hmatel=hmatel+v2det(k1(3),k1(5),k2(3),k2(5))
            hmatel=hmatel+v2det(k1(4),k1(5),k2(4),k2(5))
            do 80 i=2,5
              hmatel=hmatel+v1det(k1(i),k2(i))
80          continue
          end if
        else if (nval.eq.6) then
          if (ispin(id1,6).ne.ispin(id2,6)) return
          do 100 i=0,5 
            ii=iand(ishft(id1,-4*i),'f'x)
            k1(6-i)=ifndso(ii)
            ii=iand(ishft(id2,-4*i),'f'x)
            k2(6-i)=ifndso(ii)
100       continue
          call match6 (k1,k2,6,isign,nm)
            hmatel=v6det(k1(1),k1(2),k1(3),k1(4),k1(5),k1(6),
     1                   k2(1),k2(2),k2(3),k2(4),k2(5),k2(6))
          if (nm.ge.1) then
            hmatel=hmatel+v5det(k1(1),k1(2),k1(3),k1(4),k1(5),
     1                          k2(1),k2(2),k2(3),k2(4),k2(5))
          end if
          if (nm.ge.2) then
            hmatel=hmatel+v5det(k1(1),k1(2),k1(3),k1(4),k1(6),
     1                          k2(1),k2(2),k2(3),k2(4),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(2),k1(3),k1(4),
     1                          k2(1),k2(2),k2(3),k2(4))
          end if
          if (nm.ge.3) then
            hmatel=hmatel+v5det(k1(1),k1(2),k1(3),k1(5),k1(6),
     1                          k2(1),k2(2),k2(3),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(2),k1(3),k1(5),
     1                          k2(1),k2(2),k2(3),k2(5))
            hmatel=hmatel+v4det(k1(1),k1(2),k1(3),k1(6),
     1                          k2(1),k2(2),k2(3),k2(6))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(3),k2(1),k2(2),k2(3))
          end if
          if (nm.ge.4) then
            hmatel=hmatel+v5det(k1(1),k1(2),k1(4),k1(5),k1(6),
     1                          k2(1),k2(2),k2(4),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(2),k1(4),k1(5),
     1                          k2(1),k2(2),k2(4),k2(5))
            hmatel=hmatel+v4det(k1(1),k1(2),k1(4),k1(6),
     1                          k2(1),k2(2),k2(4),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(2),k1(5),k1(6),
     1                          k2(1),k2(2),k2(5),k2(6))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(4),k2(1),k2(2),k2(4))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(5),k2(1),k2(2),k2(5))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(6),k2(1),k2(2),k2(6))
            hmatel=hmatel+v2det(k1(1),k1(2),k2(1),k2(2))
          end if
          if (nm.ge.5) then
            hmatel=hmatel+v5det(k1(1),k1(3),k1(4),k1(5),k1(6),
     1                          k2(1),k2(3),k2(4),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(3),k1(4),k1(5),
     1                          k2(1),k2(3),k2(4),k2(5))
            hmatel=hmatel+v4det(k1(1),k1(3),k1(4),k1(6),
     1                          k2(1),k2(3),k2(4),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(3),k1(5),k1(6),
     1                          k2(1),k2(3),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(4),k1(5),k1(6),
     1                          k2(1),k2(4),k2(5),k2(6))
            hmatel=hmatel+v3det(k1(1),k1(3),k1(4),k2(1),k2(3),k2(4))
            hmatel=hmatel+v3det(k1(1),k1(3),k1(5),k2(1),k2(3),k2(5))
            hmatel=hmatel+v3det(k1(1),k1(3),k1(6),k2(1),k2(3),k2(6))
            hmatel=hmatel+v3det(k1(1),k1(4),k1(5),k2(1),k2(4),k2(5))
            hmatel=hmatel+v3det(k1(1),k1(4),k1(6),k2(1),k2(4),k2(6))
            hmatel=hmatel+v3det(k1(1),k1(5),k1(6),k2(1),k2(5),k2(6))
            do 120 i=2,5
              hmatel=hmatel+v2det(k1(1),k1(1+i),k2(1),k2(1+i))
120         continue
            hmatel=hmatel+v1det(k1(1),k2(1))
          end if
          if (nm.eq.6) then
            hmatel=hmatel+v5det(k1(2),k1(3),k1(4),k1(5),k1(6),
     1                          k2(2),k2(3),k2(4),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(2),k1(3),k1(4),k1(5),
     1                          k2(2),k2(3),k2(4),k2(5))
            hmatel=hmatel+v4det(k1(2),k1(3),k1(4),k1(6),
     1                          k2(2),k2(3),k2(4),k2(6))
            hmatel=hmatel+v4det(k1(2),k1(3),k1(5),k1(6),
     1                          k2(2),k2(3),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(2),k1(4),k1(5),k1(6),
     1                          k2(2),k2(4),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(3),k1(4),k1(5),k1(6),
     1                          k2(3),k2(4),k2(5),k2(6))
            hmatel=hmatel+v3det(k1(2),k1(3),k1(4),k2(2),k2(3),k2(4))
            hmatel=hmatel+v3det(k1(2),k1(3),k1(5),k2(2),k2(3),k2(5))
            hmatel=hmatel+v3det(k1(2),k1(3),k1(6),k2(2),k2(3),k2(6))
            hmatel=hmatel+v3det(k1(2),k1(4),k1(5),k2(2),k2(4),k2(5))
            hmatel=hmatel+v3det(k1(2),k1(4),k1(6),k2(2),k2(4),k2(6))
            hmatel=hmatel+v3det(k1(2),k1(5),k1(6),k2(2),k2(5),k2(6))
            hmatel=hmatel+v3det(k1(3),k1(4),k1(5),k2(3),k2(4),k2(5))
            hmatel=hmatel+v3det(k1(3),k1(4),k1(6),k2(3),k2(4),k2(6))
            hmatel=hmatel+v3det(k1(3),k1(5),k1(6),k2(3),k2(5),k2(6))
            hmatel=hmatel+v3det(k1(4),k1(5),k1(6),k2(4),k2(5),k2(6))
            hmatel=hmatel+v2det(k1(2),k1(3),k2(2),k2(3))
            hmatel=hmatel+v2det(k1(2),k1(4),k2(2),k2(4))
            hmatel=hmatel+v2det(k1(2),k1(5),k2(2),k2(5))
            hmatel=hmatel+v2det(k1(2),k1(6),k2(2),k2(6))
            hmatel=hmatel+v2det(k1(3),k1(4),k2(3),k2(4))
            hmatel=hmatel+v2det(k1(3),k1(5),k2(3),k2(5))
            hmatel=hmatel+v2det(k1(3),k1(6),k2(3),k2(6))
            hmatel=hmatel+v2det(k1(4),k1(5),k2(4),k2(5))
            hmatel=hmatel+v2det(k1(4),k1(6),k2(4),k2(6))
            hmatel=hmatel+v2det(k1(5),k1(6),k2(5),k2(6))
            do 110 i=2,6
              hmatel=hmatel+v1det(k1(i),k2(i))
110         continue
          end if
        end if
        hmatel=hmatel*isign
      else if (nval.ge.12) then
        hmatel=0.d0
        do 200 i=0,2
          k1(3-i)=iand(ishft(id1,-8*i),'ff'x)
          k1save(3-i)=k1(3-i)
          if (mod(k1(3-i),2).eq.0) then
            k1s(3-i)=-1
          else
            k1s(3-i)=1
          end if
          k2(3-i)=iand(ishft(id2,-8*i),'ff'x)
          k2save(3-i)=k2(3-i)
          if (mod(k2(3-i),2).eq.0) then
            k2s(3-i)=-1
          else
            k2s(3-i)=1
          end if
200     continue
        if (ihp.eq.2) then
          ihpsgn=-1
        else
          ihpsgn=1
        end if
        ididq=0
        if ((k1s(1)+k1s(2)+k2s(3)).eq.(k2s(1)+k2s(2)+k1s(3))) then  
          if (k1(3).eq.k2(3)) then
            nmb=1
          else
            nmb=0
          end if
          if ((k1s(1)+k1s(2)).eq.(k2s(1)+k2s(2))) then
            call match6 (k1,k2,2,isign,nma)
            do 210 i=1,3
              if (mod(k1(i),2).eq.0) then
                k1s(i)=-1
              else
                k1s(i)=1
              end if
              j1(i)=(k1(i)+1)/2
              if (mod(k2(i),2).eq.0) then
                k2s(i)=-1
              else
                k2s(i)=1
              end if
              j2(i)=(k2(i)+1)/2
210         continue
            call files_ccsd ("ttpf",ntape3,new,3,-1,-1)
            istrt31=nstrt(1,ntape3-1)+2*nij3+4
            istrt32=nstrt(2,ntape3-1)
            nsize3=ije3(nij3)
            hmatel=hmatel-ihpsgn*isign*
     1      find6(j1(1),j1(2),j2(3),j2(1),j2(2),j1(3))
            if (k2s(1).eq.k2s(2)) then
              hmatel=hmatel+ihpsgn*isign*
     1        find6(j1(1),j1(2),j2(3),j2(2),j2(1),j1(3))
            end if
            if (nmb.eq.1) then
              if (ihp.eq.2) then
                call files_ccsd ("thuf",ntape2,new,2,-1,-1)
              else if (ihp.eq.3) then
                call files_ccsd ("tuf ",ntape2,new,2,-1,-1)
              end if
              istrt21=nstrt(1,ntape2-1)+2*nij2+4
              istrt22=nstrt(2,ntape2-1)
              nsize2=ije2(nij2)
              hmatel=hmatel+isign*find4(j1(1),j1(2),j2(1),j2(2))
              if (k2s(1).eq.k2s(2)) then
                hmatel=hmatel-isign*find4(j1(1),j1(2),j2(2),j2(1))
              end if
              if (nma.ge.1) then
                if (ihp.eq.2) then
                  call files_ccsd ("tgf ",ntape1,new,1,-1,-1)
                else if (ihp.eq.3) then
                  call files_ccsd ("tvf ",ntape1,new,1,-1,-1)
                end if
                istrt11=nstrt(1,ntape1-1)+2*nij1+4
                istrt12=nstrt(2,ntape1-1)
                nsize1=ije1(nij1)
                hmatel=hmatel+isign*v1det(k1(1),k2(1))
                if (nma.eq.2) then   
                  hmatel=hmatel+isign*v1det(k1(2),k2(2))
                end if
              end if
            end if
            if (nma.eq.2) then
              if (ihp.eq.2) then
                call files_ccsd ("tvf ",ntape1,new,1,-1,-1)
              else if (ihp.eq.3) then
                call files_ccsd ("tgf ",ntape1,new,1,-1,-1)
              end if
              istrt11=nstrt(1,ntape1-1)+2*nij1+4
              istrt12=nstrt(2,ntape1-1)
              nsize1=ije1(nij1)
              hmatel=hmatel+isign*v1det(k1(3),k2(3))
            end if
            if (nma.ge.1) then
              call files_ccsd ("tpf ",ntape2,new,2,-1,-1)
              istrt21=nstrt(1,ntape2-1)+2*nij2+4
              istrt22=nstrt(2,ntape2-1)
              nsize2=ije2(nij2)
              hmatel=hmatel-isign*find4(j1(1),j2(3),j2(1),j1(3))
              if (nma.eq.2) then
                hmatel=hmatel-isign*find4(j1(2),j2(3),j2(2),j1(3))
              end if
            end if
            if (nma.ge.1) then
              call files_ccsd ("tqf ",ntape2,new,2,-1,-1)
              istrt21=nstrt(1,ntape2-1)+2*nij2+4
              istrt22=nstrt(2,ntape2-1)
              nsize2=ije2(nij2)
              if ((k1s(1).eq.k1s(3)).and.(k2s(1).eq.k2s(3))) then
                hmatel=hmatel+isign*find4(j1(1),j2(3),j1(3),j2(1))
              end if
              if (nma.eq.2) then
                if ((k1s(2).eq.k1s(3)).and.(k2s(2).eq.k2s(3))) then
                  hmatel=hmatel+isign*find4(j1(2),j2(3),j1(3),j2(2))
                end if
              end if
              ididq=1
            end if
          end if
          do 220 i=1,3
            k1(i)=k1save(i)
            if (mod(k1(i),2).eq.0) then
              k1s(i)=-1 
            else
              k1s(i)=1
            end if
            j1(i)=(k1(i)+1)/2
            k2(i)=k2save(i)
            if (mod(k2(i),2).eq.0) then
              k2s(i)=-1 
            else
              k2s(i)=1
            end if
            j2(i)=(k2(i)+1)/2
220       continue
          if (ididq.eq.0) then
            call files_ccsd ("tqf ",ntape2,new,2,-1,-1)
            istrt21=nstrt(1,ntape2-1)+2*nij2+4
            istrt22=nstrt(2,ntape2-1)
            nsize2=ije2(nij2)
            do 230 i=1,2
              do 240 j=1,2
                isign=1
                if (k1(i).eq.k2(j)) then
                  if ((k1s(3-i).eq.k1s(3)).
     1            and.(k2s(3-j).eq.k2s(3))) then
                    if (i.ne.j) then
                      isign=-1
                    end if
                    hmatel=hmatel+isign*
     1              find4(j1(3-i),j2(3),j1(3),j2(3-j))
                  end if
                end if
240           continue
230         continue
          end if
          call files_ccsd ("ttqf",ntape3,new,3,-1,-1)
          istrt31=nstrt(1,ntape3-1)+2*nij3+4
          istrt32=nstrt(2,ntape3-1)
          nsize3=ije3(nij3)
          if ((k1s(2).eq.k1s(3)).and.(k2s(2).eq.k2s(3))) then
            hmatel=hmatel+ihpsgn*
     1      find6(j1(1),j1(2),j2(3),j2(1),j1(3),j2(2))
          end if
          if ((k1s(2).eq.k1s(3)).and.(k2s(1).eq.k2s(3))) then
            hmatel=hmatel-ihpsgn*
     1      find6(j1(1),j1(2),j2(3),j2(2),j1(3),j2(1))
          end if
          if ((k1s(1).eq.k1s(3)).and.(k2s(2).eq.k2s(3))) then
            hmatel=hmatel-ihpsgn*
     1      find6(j1(2),j1(1),j2(3),j2(1),j1(3),j2(2))
          end if
          if ((k1s(1).eq.k1s(3)).and.(k2s(1).eq.k2s(3))) then
            hmatel=hmatel+ihpsgn*
     1      find6(j1(2),j1(1),j2(3),j2(2),j1(3),j2(1))
          end if
        end if
      end if
      return
      end 
      function idec(ihex)
      implicit real*8 (a-h,o-z)
      i1=iand(ishft(ihex,-24),'ff'X)
      i2=iand(ishft(ihex,-16),'ff'X)
      i3=iand(ishft(ihex,-8) ,'ff'X)
      i4=iand(          ihex ,'ff'X)
      if (i3.eq.0) then
        idec=i1*100+i2
      else
        idec=i1*1000000+i2*10000+i3*100+i4
      endif
      return
      end
      function ifndso(ii)
      implicit real*8 (a-h,o-z)
      common /idv/ iv(200),nvalo
      if ((ii/2.d0-int(ii/2.d0)).gt.0.1d0) then
        iflag=1
      else
        iflag=0
      end if
      ifndso=2*iv((ii+1)/2)-iflag
      return
      end
      function ifndval(ii)
      implicit real*8 (a-h,o-z)
      common /idv/ iv(200),nvalo
      do 10 j=1,nvalo
        if (ii.eq.iv(j)) then
          ifndval=j             
        end if
10    continue
      return
      end 
      function imax (id1,ind1,id2,ind2,iparm)
      implicit real*8 (a-h,o-z)
	integer*4 ch1,ch2
      if (id1.gt.id2) then
        imax=1
      else if (id2.gt.id1) then
        imax=2
      else
        if (ind1.gt.ind2) then
          imax=1
        else if (ind2.gt.ind1) then
          imax=2
        else
          if (iparm.eq.1) then
            write (6,*) 'identical input values in imax'
            stop
          else if (iparm.eq.2) then
            imax=0
          end if
        end if
      end if
      return
      end
      subroutine imove(ntapef,ntapet)
      implicit real*8 (a-h,o-z)
!     integer*4 ch1,ch2
      integer ch1,ch2
      
#include "param_ccsd.h"

      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      nija=igen(nstrt(1,ntapef-1)+1)
      call kcopy(nija+1,igen(nstrt(1,ntapef-1)+1),nij1)
      call kcopy(nija+1,igen(nstrt(1,ntapef-1)+nija+2),ije1(0))
      iii1=igen(nstrt(1,ntapef-1)+2*nija+3)
      iiv1=igen(nstrt(1,ntapef-1)+2*nija+4)
      nstrt(1,ntapet)=nstrt(1,ntapet-1)+2*nij1+4+ije1(nij1)
      nstrt(2,ntapet)=nstrt(2,ntapet-1)+ije1(nij1)
      call dimx(nstrt(1,ntapet),nbig,"imove1")
      call dimx(nstrt(2,ntapet),nbig,"imove2")
c   special for scratch files, which start with "s".
           itmp=4H       
           jtmp='ffffff00'X
           ch1 = iand(itmp,jtmp)
           ch2 = iand (namefl(ntapet), '000000ff'X)
           isaa = ior (ch1, ch2)

      i_4hs=4hs
!     if (isaa.eq.4hs   ) then
      if (isaa.eq.i_4hs) then
        nstrt(1,ntapet)=nstrt(1,ntapet-1)+ndp3+2*ndij+4
        nstrt(2,ntapet)=nstrt(2,ntapet-1)+ndp3
        call dimx(nij1,ndp3,"imove   ")
        call dimx(nstrt(1,ntapet),nbig,"imove1  ")
        call dimx(nstrt(2,ntapet),nbig,"imove2  ")
      endif
      call kcopy(ije1(nij1)+2*nij1+4,igen(nstrt(1,ntapef-1)+1),
     1 igen(nstrt(1,ntapet-1)+1))
      return
      end

      function ispin(idet,nval)
      implicit real*8 (a-h,o-z)
      dimension jj(3)
      if (nval.le.4) then
        ispin=0
        do 10 i=1,nval   
          j=iand(ishft(idet,-8*(i-1)),'ff'x)
          if (mod(j,2).eq.1) then
            ispin=ispin+1
          else if (mod(j,2).eq.0) then
            ispin=ispin-1
          end if
10      continue
      else if ((nval.eq.5).or.(nval.eq.6)) then
        ispin=0
        do 20 i=0,nval-1
          j=iand(ishft(idet,-4*i),'f'x)
          if (mod(j,2).eq.1) then
            ispin=ispin+1
          else if (mod(j,2).eq.0) then
            ispin=ispin-1
          end if
20      continue
      else if (nval.gt.6) then    
        if (nval.eq.12) then
          isign=-1
        else if (nval.eq.21) then
          isign=1
        end if
        ispin=0
        jj(1)=iand(ishft(idet,-16),'ff'x)
        jj(2)=iand(ishft(idet,-8),'ff'x)
        jj(3)=iand(idet,'ff'x)
        do 30 i=1,2
          if (mod(jj(i),2).eq.1) then
            ispin=ispin+isign
          else if (mod(jj(i),2).eq.0) then
            ispin=ispin-isign
          end if
30      continue
        if (mod(jj(3),2).eq.1) then
          ispin=ispin-isign
        else if (mod(jj(3),2).eq.0) then
          ispin=ispin+isign
        end if
      end if
      return
      end
      subroutine ivsort(ivec,n8)
C
      IMPLICIT INTEGER (A-Z)
!     INTEGER*4  ivec(*),TEMP,N8
      INTEGER  ivec(*),TEMP,N8
      integer lftbit,nx4,n
C
      n=n8
      IF (N .LE. 1) RETURN
C     P0 = AND(N-1,MASK(LEADZ(N-1)+1))
      nx4=n-1
      P0 = lftbit(nx4)
      P = P0
C
   20 Q = P0
      R = 0
      D = P
      I = 0
C
   30 IF (IAND(I,P) .EQ. R) THEN
C$DIR NO_RECURRENCE
         DO 32 J=1,MIN (P,N-I-D)
c           TEMP       = MAX (VEC(J+I),VEC(J+I+D))
c           VEC(J+I)   = MIN (VEC(J+I),VEC(J+I+D))
c           VEC(J+I+D) = TEMP
        if (ishft( ivec(j+i),-1).gt.ishft( ivec(j+i+d),-1)) then
	    temp       = ivec(j+i)
	    ivec(j+i)   = ivec(j+i+d)
            ivec(J+I+D) = TEMP
        Else If ((ishft(ivec(j+i),-1).eq.ishft(ivec(j+i+d),-1))
     +  .AND. (iabs(ivec(j+i)).gt.iabs(ivec(j+i+d)))) then
           temp       = ivec(j+i)
            ivec(j+i)   = ivec(j+i+d)
            ivec(J+I+D) = TEMP
        endif
   32    CONTINUE
      ENDIF
      I = I + P
      IF ((I+D) .LT. N) GO TO 30
C
      IF (Q .NE. P) THEN
        D = Q-P
        Q = Q/2
        R = P
        I = 0
        GO TO 30
      ENDIF
C
      P = P/2
      IF (P .GT. 4) GO TO 20
C
   70 Q = P0
      R = 0
      D = P
C
   80 DO 82 I=1,P
C$DIR NO_RECURRENCE
      DO 82 J=I,N-D-R,2*P
c        TEMP       = MAX (VEC(J+R),VEC(J+R+D))
c        VEC(J+R)   = MIN (VEC(J+R),VEC(J+R+D))
c        VEC(J+R+D) = TEMP
        if (ishft( ivec(j+r),-1).gt.ishft( ivec(j+r+d),-1)) then
	    temp       = ivec(j+r)
	    ivec(j+r)   = ivec(j+r+d)
            ivec(J+R+D) = TEMP
        Else If ((ishft(ivec(j+r),-1).eq.ishft(ivec(j+r+d),-1))
     +  .AND. (iabs(ivec(j+r)).gt.iabs(ivec(j+r+d)))) then
            temp       = ivec(j+r)
            ivec(j+r)   = ivec(j+r+d)
            ivec(J+R+D) = TEMP
        endif
   82 CONTINUE
C
      IF (Q .NE. P) THEN
         D = Q - P
         Q = Q / 2
         R = P
         GO TO 80
      ENDIF
C
      P = P / 2
      IF (P .NE. 0) GO TO 70
C
      RETURN
      END
      subroutine kcopy(n,ia,ib)
      implicit real*8 (a-h,o-z)
      dimension ia(*),ib(*)
      do 10 i=1,n
10    ib(i)=ia(i)
      return
      end
	function lftbit(n)
	integer*4 n,i,ipos
	integer*4 nlim(32)
	data nlim/'0'o,'1'o,'2'o,'4'o,
     x             '10'o,'20'o,'40'o,
     x             '100'o,'200'o,'400'o,
     x             '1000'o,'2000'o,'4000'o,
     x             '10000'o,'20000'o,'40000'o,
     x             '100000'o,'200000'o,'400000'o,
     x             '1000000'o,'2000000'o,'4000000'o,
     x             '10000000'o,'20000000'o,'40000000'o,
     x             '100000000'o,'200000000'o,'400000000'o,
     x             '1000000000'o,'2000000000'o,'4000000000'o,
     x             '10000000000'o/
	ipos=0
	do i=1,32
	if (n.ge.nlim(i)) ipos=i
	enddo
	lftbit=nlim(ipos)
	return
	end
      subroutine match6 (k1,k2,nval,isign,nm)
      implicit real *8 (a-h,o-z)
      dimension k1(6),k2(6),m1(6),m2(6),mm1(6),mm2(6),
     1k1save(6),k2save(6)
      do 10 i=1,nval
        k1save(i)=k1(i)  
        k2save(i)=k2(i)  
10    continue
      im=0
      imm=0
      do 30 i=1,nval
        do 80 j=1,nval
          if ((k1(i).eq.k2(j)).and.(k1(i).ne.0)) then
            im=im+1
            m1(im)=k1(i)
            m2(im)=k2(j)
            k1(i)=0
            k2(j)=0
          end if
80      continue
30    continue
      do 40 i=1,nval
        if (k1(i).ne.0) then
          imm=imm+1
          mm1(imm)=k1(i)
          ididit=0
          do 70 j=1,nval
            if (ididit.eq.0) then
              if (k2(j).ne.0) then
                if (mod(k1(i),2).eq.mod(k2(j),2)) then
                  mm2(imm)=k2(j)
                  k2(j)=0
                  ididit=1
                end if
              end if
            end if
70        continue
        end if
40    continue
      nm=im
      nmm=imm
      do 50 i=1,nmm
        k1(i)=mm1(i)
        k2(i)=mm2(i)
50    continue
      do 60 i=1,nm
        k1(nmm+i)=m1(i)
        k2(nmm+i)=m2(i)
60    continue
      call perm (k1save,k1,isign1,nval)
      call perm (k2save,k2,isign2,nval)
      isign=isign1*isign2
      return
      end
      subroutine move(namef,namet)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      nbias=2*ndij+4
      call files_ccsd(namef,ntapef,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) ' move input file not found'
        call abext('move    ')
      endif
      call files_ccsd(namet,ntapet,new,0,iii1,iiv1)
      if (new.eq.1) then
        write (6,*) ' move output file of unknown type'
        call abext('move    ')
      endif
      nsize=ije1(nij1)+2*nij1+4
      call kcopy(nsize,igen(nstrt(1,ntapef-1)+1),
     1 igen(nstrt(1,ntapet-1)+1))
      call dcopy(ije1(nij1),gen(nstrt(2,ntapef-1)+1),1,
     1 gen(nstrt(2,ntapet-1)+1),1)
      return
      end

      subroutine rle (niter,irleop,name1,name2,name3,name4,ising)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      double precision r(31,31),z(31),a(31),rcond,rmult

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
CUK   common /crle/ trle(nrle,4),isize(4)
      common /crle/ trle1(nrle1),trle2(nrle2),trle3(nrle3),
     1 trle4(nrle4),isize(4)
      integer*4 ipvt(31)
      dimension iopa(4),name(4)
      call prtime ("rle in  ",0,0)
      name(1)=name1
      name(2)=name2
      name(3)=name3
      name(4)=name4
      do 10 i=1,niter
        a(i)=0.d0
        do 20 j=1,niter
          r(i,j)=0.d0
20      continue
10    continue
      call parse (irleop,iopa(1),iopa(2),iopa(3),iopa(4))
      if (iopa(1).eq.1) then
        do 41 i=1,niter
          a(i)=a(i)-ddot(isize(1),trle1(isize(1)+1),1,
     1    trle1((i+1)*isize(1)+1),1)
41      continue
      end if
      if (iopa(2).eq.1) then
        do 42 i=1,niter
          a(i)=a(i)-ddot(isize(2),trle2(isize(2)+1),1,
     1    trle2((i+1)*isize(2)+1),1)
42      continue
      end if
      if (iopa(3).eq.1) then
        do 43 i=1,niter
          a(i)=a(i)-ddot(isize(3),trle3(isize(3)+1),1,
     1    trle3((i+1)*isize(3)+1),1)
43      continue
      end if
      if (iopa(4).eq.1) then
        do 44 i=1,niter
          a(i)=a(i)-ddot(isize(4),trle4(isize(4)+1),1,
     1    trle4((i+1)*isize(4)+1),1)
44      continue
      end if
      do 150 i=1,niter
        do 160 j=1,niter
          r(i,j)=r(i,j)+a(i)
160     continue
150   continue
      if (iopa(1).eq.1) then
        lngth=isize(1)
        do 61 i=1,niter
          do 81 k=1,lngth
            trle1((niter+3)*lngth+k)=
     1      trle1((i+1)*lngth+k)*trle1(k)
81        continue
          do 71 j=1,niter
            if (j.eq.1) then
              rold=ddot(lngth,trle1((niter+3)*lngth+1),1,
     1        trle1(2*lngth+1),1)
            end if
            rnew=ddot(lngth,trle1((niter+3)*lngth+1),1,
     1      trle1((j+2)*lngth+1),1)
            r(i,j)=r(i,j)+rold-rnew
            rold=rnew
71        continue
61      continue
      end if
      if (iopa(2).eq.1) then
        lngth=isize(2)
        do 62 i=1,niter
          do 82 k=1,lngth
            trle2((niter+3)*lngth+k)=
     1      trle2((i+1)*lngth+k)*trle2(k)
82        continue
          do 72 j=1,niter
            if (j.eq.1) then
              rold=ddot(lngth,trle2((niter+3)*lngth+1),1,
     1        trle2(2*lngth+1),1)
            end if
            rnew=ddot(lngth,trle2((niter+3)*lngth+1),1,
     1      trle2((j+2)*lngth+1),1)
            r(i,j)=r(i,j)+rold-rnew
            rold=rnew
72        continue
62      continue
      end if
      if (iopa(3).eq.1) then
        lngth=isize(3)
        do 63 i=1,niter
          do 83 k=1,lngth
            trle3((niter+3)*lngth+k)=
     1      trle3((i+1)*lngth+k)*trle3(k)
83        continue
          do 73 j=1,niter
            if (j.eq.1) then
              rold=ddot(lngth,trle3((niter+3)*lngth+1),1,
     1        trle3(2*lngth+1),1)
            end if
            rnew=ddot(lngth,trle3((niter+3)*lngth+1),1,
     1      trle3((j+2)*lngth+1),1)
            r(i,j)=r(i,j)+rold-rnew
            rold=rnew
73        continue
63      continue
      end if
      if (iopa(4).eq.1) then
        lngth=isize(4)
        do 64 i=1,niter
          do 84 k=1,lngth
            trle4((niter+3)*lngth+k)=
     1      trle4((i+1)*lngth+k)*trle4(k)
84        continue
          do 74 j=1,niter
            if (j.eq.1) then
              rold=ddot(lngth,trle4((niter+3)*lngth+1),1,
     1        trle4(2*lngth+1),1)
            end if
            rnew=ddot(lngth,trle4((niter+3)*lngth+1),1,
     1      trle4((j+2)*lngth+1),1)
            r(i,j)=r(i,j)+rold-rnew
            rold=rnew
74        continue
64      continue
      end if
      call dgeco (r,31,niter,ipvt,rcond,z)
      write (6,*) 'rcond=',real(rcond)
      ising=0
      if (rcond.lt.1.d-16) then     
        ising=1
      end if
      if (ising.eq.1) then    
        return
      end if
      call dgesl (r,31,niter,ipvt,a,0)
      do 100 itape=1,4
        if (iopa(itape).eq.1) then
          call files_ccsd (name(itape),ntape,new,1,-1,-1)
          istrt=nstrt(2,ntape-1)
          lngth=isize(itape)
          call dcopy (lngth,0.d0,0,gen(istrt+1),1)
          if (itape.eq.1) then
            do 111 i=1,niter
              call daxpy (lngth,a(i),trle1((i+1)*lngth+1),1,
     1        gen(istrt+1),1)
111         continue
          else if (itape.eq.2) then
            do 112 i=1,niter
              call daxpy (lngth,a(i),trle2((i+1)*lngth+1),1,
     1        gen(istrt+1),1)
112         continue
          else if (itape.eq.3) then
            do 113 i=1,niter
              call daxpy (lngth,a(i),trle3((i+1)*lngth+1),1,
     1        gen(istrt+1),1)
113         continue
          else if (itape.eq.4) then
            do 114 i=1,niter
              call daxpy (lngth,a(i),trle4((i+1)*lngth+1),1,
     1        gen(istrt+1),1)
114         continue
          end if
        end if
100   continue
      call prtime ("rle out ",0,0)
      return
      end
      subroutine rvmuls(n,a,s,c)
      implicit real*8 (a-h,o-z)
      dimension a(*),c(*)
      do 10 i=1,n
10    c(i)=a(i)*s
      return
      end
      subroutine rvmulv(n,a,b,c)
      implicit real*8 (a-h,o-z)
      dimension a(*),b(*),c(*)
      do 10 i=1,n
10    c(i)=a(i)*b(i)
      return
      end

!      subroutine unfold(namef,namet,iline)
      subroutine unfold_ccsd(namef,namet,iline)
      implicit real*8 (a-h,o-z)
c   expand file namef to namet, turning line iline from valence
c   to general. (iline=4 is last line).
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      call files_ccsd(namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) ' unfold input file not found'
        call abext('unfold  ')
      endif
      iii1=iii2
      iiv1=iiv2-10**(4-iline)
      call files_ccsd(namet,ntapet,new,1,iii1,iiv1)
      nsize1=ije1(nij1)
      nsize2=ije2(nij2)
      nstrt11=nstrt(1,ntapet-1)+2*nij1+4
      nstrt12=nstrt(2,ntapet-1)
      nstrt21=nstrt(1,ntapef-1)+2*nij2+4
      nstrt22=nstrt(2,ntapef-1)
      call expnd(igen(nstrt11+1),igen(nstrt21+1),
     1  gen(nstrt12+1),gen(nstrt22+1),nsize1,nsize2)
      return
      end

      subroutine vands(n,ia,is,ic)
      implicit real*8 (a-h,o-z)
      dimension ia(*),ic(*)
      do 10 i=1,n
10    ic(i)=iand(ia(i),is)
      return
      end
      subroutine viorv(n,ia,ib,ic)
      implicit real*8 (a-h,o-z)
      dimension ia(*),ib(*),ic(*)
      do 10 i=1,n
10    ic(i)=ior(ia(i),ib(i))
      return
      end
      subroutine vshfcs(n,ia,j,ic)
      implicit real*8 (a-h,o-z)
      dimension ia(*),ic(*)
      do 10 i=1,n
10    ic(i)=ishft(ia(i),j)
      return
      end
      subroutine vsort (ivec,N8,rlvec)
C
!     IMPLICIT INTEGER*4 (A-Z)
      IMPLICIT INTEGER (A-Z)
!     INTEGER*4  ivec(*),TEMP
      INTEGER  ivec(*),TEMP
!     integer*4 N8
      integer N8
      real*8 rlvec(*),rltemp
C      integer lftbit,nx4,n
C
      n=n8
      IF (N .LE. 1) RETURN
C     P0 = AND(N-1,MASK(LEADZ(N-1)+1))
      nx4=n-1
      P0 = lftbit(nx4)
      P = P0
C
   20 Q = P0
      R = 0
      D = P
      I = 0
C
   30 IF (IAND(I,P) .EQ. R) THEN
C$DIR NO_RECURRENCE
         DO 32 J=1,MIN (P,N-I-D)
c           TEMP       = MAX (VEC(J+I),VEC(J+I+D))
c           VEC(J+I)   = MIN (VEC(J+I),VEC(J+I+D))
c           VEC(J+I+D) = TEMP
CGJ
C Now works for 32 bits
	if (ishft( ivec(j+i),-1).gt.ishft( ivec(j+i+d),-1)) then
	    temp       = ivec(j+i)
	    ivec(j+i)   = ivec(j+i+d)
            ivec(J+I+D) = TEMP
	    rltemp     = rlvec(j+i)
	    rlvec(j+i) = rlvec(j+i+d)
	    rlvec(j+i+d) = rltemp
        Else If ((ishft(ivec(j+i),-1).eq.ishft(ivec(j+i+d),-1))
     +  .AND. (iabs(ivec(j+i)).gt.iabs(ivec(j+i+d)))) then
            temp       = ivec(j+i)
            ivec(j+i)   = ivec(j+i+d)
            ivec(J+I+D) = TEMP
            rltemp     = rlvec(j+i)
            rlvec(j+i) = rlvec(j+i+d)
            rlvec(j+i+d) = rltemp

        endif
   32    CONTINUE
      ENDIF
      I = I + P
      IF ((I+D) .LT. N) GO TO 30
C
      IF (Q .NE. P) THEN
        D = Q-P
        Q = Q/2
        R = P
        I = 0
        GO TO 30
      ENDIF
C
      P = P/2
      IF (P .GT. 4) GO TO 20
C
   70 Q = P0
      R = 0
      D = P
C
   80 DO 82 I=1,P
C$DIR NO_RECURRENCE
      DO 82 J=I,N-D-R,2*P
c        TEMP       = MAX (VEC(J+R),VEC(J+R+D))
c        VEC(J+R)   = MIN (VEC(J+R),VEC(J+R+D))
c        VEC(J+R+D) = TEMP
        if (ishft( ivec(j+r),-1).gt.ishft( ivec(j+r+d),-1)) then
	    temp       = ivec(j+r)
	    ivec(j+r)   = ivec(j+r+d)
            ivec(J+R+D) = TEMP
	    rltemp     = rlvec(j+r)
	    rlvec(j+r) = rlvec(j+r+d)
	    rlvec(j+r+d) = rltemp
        Else If ((ishft(ivec(j+r),-1).eq.ishft(ivec(j+r+d),-1))
     +  .AND. (iabs(ivec(j+r)).gt.iabs(ivec(j+r+d)))) then
            temp       = ivec(j+r)
            ivec(j+r)   = ivec(j+r+d)
            ivec(J+R+D) = TEMP
            rltemp     = rlvec(j+r)
            rlvec(j+r) = rlvec(j+r+d)
            rlvec(j+r+d) = rltemp

        endif
   82 CONTINUE
C
      IF (Q .NE. P) THEN
         D = Q - P
         Q = Q / 2
         R = P
         GO TO 80
      ENDIF
C
      P = P / 2
      IF (P .NE. 0) GO TO 70
C
      RETURN
      END
      subroutine vxsort(itmp,ntmp,nt,t,it)
      implicit real*8 (a-h,o-z)
c  sort itmp, expnd to it, gather t accordingly.
      dimension itmp(ntmp),t(nt),it(nt)
      if (ntmp.gt.1) call vsort(itmp,ntmp,t)
C      print '(2hGG,Z8)',(itmp(i),i=1,ntmp)
      if (ntmp.gt.nt) then
        write (*,*) ' vxsort. ntmp gt nt! ',ntmp,nt
        call zeroup(ntmp,itmp,t,1.d-8)
        write (*,*) ntmp,' integrals left'
        if (ntmp.gt.nt) then
          call zeroup(ntmp,itmp,t,5.d-6)
          write (*,*) ntmp,' integrals left'
          if (ntmp.gt.nt) then
            write (*,*) ' vxsort. ntmp still gt nt ',ntmp,nt
            call abext('vxsort  ')
          endif
        endif
      endif
      if (ntmp.lt.nt) then
        lt=nt
        ltmp=ntmp
20      if (itmp(ltmp).eq.it(lt)) then
          t(lt)=t(ltmp)
          lt=lt-1
          ltmp=ltmp-1
        elseif ((ltmp.eq.0).or.
     +    (ishft(itmp(ltmp),-1).lt.ishft(it(lt),-1)).or.
     +    ((ishft(itmp(ltmp),-1).eq.ishft(it(lt),-1)).and.
     +    (iabs(itmp(ltmp)).lt.iabs(it(lt)))))  then
          t(lt)=0.d0
          lt=lt-1
        elseif (abs(t(ltmp)).lt.2.d-5) then
          ltmp=ltmp-1
        else
	Do i=-20,20
          write (*,100) lt,nt,ltmp+i,ntmp+i,idec(itmp(ltmp+i)),
     +                   idec(it(ntmp+i)),t(ltmp+i)
100       format(' label not in tmplt. lt nt ltmp ntmp itmp it t=',
     1    4i7,1x,2i9,e12.4)
	End Do
          call abext('vxsort  ')   
        endif
        if (lt.gt.0) go to 20
      endif
      return
      end
      subroutine wrtfile(nt,nrec,name)
      implicit real*8 (a-h,o-z)
c  write internal file "name" on tape nt, skipping nrec records.
c  disk file structure for tc etc:
c  1. nij (# of subblocks)
c  2. ij1(1:nij), subblock id's; ije1(0:nij), subblock ends; iii, iiv.
c  3. ia(1:ije1(nij)), indices.
c  4. a(1:ije1(nij)), actual values.
c    nt.eq.9 = file witout indices (pppp).
c  // contains all files. file namefl(i) starts at nstrt(i-1)+1 (see /cfiles/).
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      rewind nt
      if (nrec.gt.0) then
        do 10 i=1,4*nrec
10      read (nt)
      endif
      call files_ccsd(name,ntape,new,0,-1,-1)
      if (new.ne.0) then
        write (*,*) ' file to be written out not found'
        call abext('wrtfile ')
      endif
      nstr1=nstrt(1,ntape-1)
      nstr2=nstrt(2,ntape-1)
      nij=igen(nstr1+1)
      write (nt) nij
      write (nt) (igen(nstr1+i),i=2,2*nij+4)
      nsize=igen(nstr1+2*nij+2)
      if (nt.ne.9) then
        write (nt) (igen(nstr1+2*nij+4+i),i=1,nsize)
        write (nt) (gen(nstr2+i),i=1,nsize)
      else
        write (nt) (gen(nstr2+i),i=1,nsize)
      endif
      rewind nt
      return
      end
      subroutine wrtrle (name,itape,loc) 
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
CUK   common /crle/ trle(nrle,4),isize(4)
      common /crle/ trle1(nrle1),trle2(nrle2),trle3(nrle3),
     1 trle4(nrle4),isize(4)
      call files_ccsd (name,ntape,new,1,-1,-1)
      istrt=nstrt(2,ntape-1)
      lngth=isize(itape)
      if (itape.eq.1) then
        call dimx((loc+2)*lngth,nrle1,"nrle1wrt")
        do 11 i=1,lngth
          trle1((loc+1)*lngth+i)=gen(istrt+i)
11      continue
      else if (itape.eq.2) then
        call dimx((loc+2)*lngth,nrle2,"nrle2wrt")
        do 12 i=1,lngth
          trle2((loc+1)*lngth+i)=gen(istrt+i)
12      continue
      else if (itape.eq.3) then
        call dimx((loc+2)*lngth,nrle3,"nrle3wrt")
        do 13 i=1,lngth
          trle3((loc+1)*lngth+i)=gen(istrt+i)
13      continue
      else if (itape.eq.4) then
        call dimx((loc+2)*lngth,nrle4,"nrle4wrt")
        do 14 i=1,lngth
          trle4((loc+1)*lngth+i)=gen(istrt+i)
14      continue
      end if
      return     
      end 
      subroutine zero (name)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij,ij(ndij),ije(0:ndij),iii,iiv
      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      call files_ccsd (name,ntape,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) 'zero input tape does not exist, ',name
        stop                      
      end if
      istrt1=nstrt(1,ntape-1)+2*nij+4
      istrt2=nstrt(2,ntape-1)
      if (iii.le.333333) then
        do 10 i=1,nij   
          i1=iand(ishft(ij(i),-16),'ff'x)
          i2=iand(ishft(ij(i),-8),'ff'x)
          i3=iand(ij(i),'ff'x)
          do 20 j=ije(i-1)+1,ije(i)      
            i4=iand(ishft(igen(istrt1+j),-16),'ff'x)
            i5=iand(ishft(igen(istrt1+j),-8),'ff'x)
            i6=iand(igen(istrt1+j),'ff'x)
            if ((icv(i1).eq.1).and.(icv(i2).eq.1).and.
     1      (icv(i3).eq.1).and.(icv(i4).eq.1).and.
     2      (icv(i5).eq.1).and.(icv(i6).eq.1)) then
              gen(istrt2+j)=0.d0
            end if
20        continue
10      continue
      else if (iii.gt.333333) then
        do 30 i=1,nij   
          i1=iand(ishft(ij(i),-24),'ff'x)
          i2=iand(ishft(ij(i),-16),'ff'x)
          i3=iand(ishft(ij(i),-8),'ff'x)
          i4=iand(ij(i),'ff'x)
          do 40 j=ije(i-1)+1,ije(i)      
            i5=iand(ishft(igen(istrt1+j),-24),'ff'x)
            i6=iand(ishft(igen(istrt1+j),-16),'ff'x)
            i7=iand(ishft(igen(istrt1+j),-8),'ff'x)
            i8=iand(igen(istrt1+j),'ff'x)
            if ((icv(i1).eq.1).and.(icv(i2).eq.1).and.
     1      (icv(i3).eq.1).and.(icv(i4).eq.1).and.
     2      (icv(i5).eq.1).and.(icv(i6).eq.1).and.   
     3      (icv(i7).eq.1).and.(icv(i8).eq.1)) then
              gen(istrt2+j)=0.d0
            end if
40        continue
30      continue
      end if
      return
      end
      subroutine zeroup(ntmp,itmp,t,zero)
      implicit real*8 (a-h,o-z)
      dimension itmp(ntmp),t(ntmp)
      nn=0
      do 120 i=1,ntmp
      if (abs(t(i)).ge.zero) then
        nn=nn+1
        itmp(nn)=itmp(i)
        t(nn)=t(i)
      endif
120   continue
      ntmp=nn
      return
      end
      subroutine parse(ii,i1,i2,i3,i4)
      implicit real*8 (a-h,o-z)
      ia=ii
      i4=mod(ia,10)
      ia=ia/10
      i3=mod(ia,10)
      ia=ia/10
      i2=mod(ia,10)
      i1=ia/10
      return
      end
      subroutine parse6 (ii,i1,i2,i3,i4,i5,i6)
      implicit real*8 (a-h,o-z)
      ia=ii 
      i6=mod(ia,10)
      ia=ia/10
      i5=mod(ia,10)
      ia=ia/10
      i4=mod(ia,10)
      ia=ia/10
      i3=mod(ia,10)
      ia=ia/10
      i2=mod(ia,10)
      i1=ia/10
      return    
      end
      subroutine parse8 (ii,i1,i2,i3,i4,i5,i6,i7,i8)
      implicit real*8 (a-h,o-z)
      ia=ii 
      i8=mod(ia,10)
      ia=ia/10
      i7=mod(ia,10)
      ia=ia/10
      i6=mod(ia,10)
      ia=ia/10
      i5=mod(ia,10)
      ia=ia/10
      i4=mod(ia,10)
      ia=ia/10
      i3=mod(ia,10)
      ia=ia/10
      i2=mod(ia,10)
      i1=ia/10
      return    
      end
      subroutine perm (i1,i2,isign,nval1)
      implicit real*8 (a-h,o-z)
      dimension i1(6),i2(6),j1(6),j2(6)
      do 10 i=1,nval1
        j1(i)=i1(i)
        j2(i)=i2(i)
10    continue
      isign=1
      do 20 i=1,nval1
        do 30 j=1,nval1
          if (j1(i).eq.j2(j)) then
            jsave=j
          end if
30      continue
        nflips=jsave-i
        do 40 k=nflips,1,-1
          ksave=j2(i+k)
          j2(i+k)=j2(i+k-1)
          j2(i+k-1)=ksave
          isign=-1*isign
40      continue
20    continue
      return
      end
      subroutine perms(n)
      implicit real*8 (a-h,o-z)
      common /prm/ iperm(0:6,720)
      dimension k1(6),k2(6)
      if (n.eq.1) then
        iperm(0,1)=1
        iperm(1,1)=1
      else if (n.eq.2) then
        iperm(0,1)=1
        iperm(1,1)=1
        iperm(2,1)=2
        iperm(0,2)=-1
        iperm(1,2)=2
        iperm(2,2)=1
      else if (n.eq.3) then
        do 170 i=1,3
          k1(i)=i
170     continue
        ii=0
        do 10 i=1,3
          do 20 j=1,3
            if (i.ne.j) then 
              ii=ii+1
              k=6-i-j
              iperm(1,ii)=i
              iperm(2,ii)=j
              iperm(3,ii)=k
              k2(1)=i
              k2(2)=j
              k2(3)=k
              call perm(k1,k2,isign,3)
              iperm(0,ii)=isign
            end if
20        continue
10      continue
      else if (n.eq.4) then
        do 180 i=1,4
          k1(i)=i
180     continue
        ii=0
        do 30 i=1,4
          do 40 j=1,4
            if (i.ne.j) then 
              do 50 k=1,4
                if ((i.ne.k).and.(j.ne.k)) then
                  ii=ii+1
                  l=10-i-j-k
                  iperm(1,ii)=i
                  iperm(2,ii)=j
                  iperm(3,ii)=k
                  iperm(4,ii)=l
                  k2(1)=i
                  k2(2)=j
                  k2(3)=k
                  k2(4)=l
                  call perm (k1,k2,isign,4)
                  iperm(0,ii)=isign
                end if
50            continue
            end if
40        continue
30      continue
      else if (n.eq.5) then
        do 150 i=1,5
          k1(i)=i
150     continue
        ii=0
        do 60 i=1,5
          do 70 j=1,5
            if (i.ne.j) then 
              do 80 k=1,5
                if ((i.ne.k).and.(j.ne.k)) then
                  do 90 l=1,5
                    if ((i.ne.l).and.(j.ne.l).and.(k.ne.l)) then
                      ii=ii+1
                      m=15-i-j-k-l
                      iperm(1,ii)=i
                      iperm(2,ii)=j
                      iperm(3,ii)=k
                      iperm(4,ii)=l
                      iperm(5,ii)=m
                      k2(1)=i
                      k2(2)=j
                      k2(3)=k
                      k2(4)=l
                      k2(5)=m
                      call perm(k1,k2,isign,5)
                      iperm(0,ii)=isign
                    end if
90                continue
                end if
80            continue
            end if
70        continue
60      continue
      else if (n.eq.6) then
        do 160 i=1,6
          k1(i)=i
160     continue
        ii=0
        do 100 i=1,6
          do 110 j=1,6
            if (i.ne.j) then 
              do 120 k=1,6
                if ((i.ne.k).and.(j.ne.k)) then
                  do 130 l=1,6
                    if ((i.ne.l).and.(j.ne.l).and.(k.ne.l)) then
                      do 140 m=1,6
                        if ((i.ne.m).and.(j.ne.m).and.(k.ne.m).and.
     1                  (l.ne.m)) then
                          ii=ii+1
                          n=21-i-j-k-l-m
                          iperm(1,ii)=i
                          iperm(2,ii)=j
                          iperm(3,ii)=k
                          iperm(4,ii)=l
                          iperm(5,ii)=m
                          iperm(6,ii)=n
                          k2(1)=i
                          k2(2)=j
                          k2(3)=k
                          k2(4)=l
                          k2(5)=m
                          k2(6)=n
                          call perm(k1,k2,isign,6)    
                          iperm(0,ii)=isign
                        end if
140                   continue
                    end if
130               continue
                end if
120           continue
            end if
110       continue
100     continue
      end if
      return
      end
      subroutine prt(name,tit,iparm)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij,ij(ndij),ije(0:ndij),iii,iiv
      dimension a(6),ind(4,6)
      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      if (iparm.gt.iopt(1)) return
      flim=10.d0**(-iopt(5))
      if (iopt(5).le.0) flim=0.d0
      call files_ccsd(name,num,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) ' prt tried to read nonexistent file'
        call abext('prt     ')
      endif
      nsize=ije(nij)
      write (6,110) tit,nsize,iii,iiv,nij,(ij(i),ije(i),i=1,nij)
110   format(' file ',a8,', size iii iiv=',4i7,' ij records:'/
     1 (8(2x,z4,i8,';')))
      istrt1=nstrt(1,num-1)+2*nij+4
      istrt2=nstrt(2,num-1)
      io=0
      do 30 i=1,nsize
      if (abs(gen(istrt2+i)).ge.flim) then
        io=io+1
        a(io)=gen(istrt2+i)
	ind(1,io)=iand(ishft(igen(istrt1+i),-24),'ff'X)
	ind(2,io)=iand(ishft(igen(istrt1+i),-16),'ff'X)
	ind(3,io)=iand(ishft(igen(istrt1+i),-8) ,'ff'X)
	ind(4,io)=iand(           igen(istrt1+i),'ff'X)
        if (io.eq.6) then
          write(6,140) ((ind(jj,ii),jj=1,4),a(ii),ii=1,6)
          io=0
        endif
      endif
30    continue
      if (io.gt.0) write(6,140) ((ind(jj,ii),jj=1,4),a(ii),ii=1,io)
140   format(1x,6(1x,4i2,'=',f11.7))
      return
      end
      subroutine prt6 (name,tit,iparm)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij,ij(ndij),ije(0:ndij),iii,iiv
      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      dimension a(4),id(3),ind(3,4)
      if (iparm.gt.iopt(1)) return
      flim=10.d0**(-iopt(5))
      if (iopt(5).le.0) flim=0.d0
      call files_ccsd(name,num,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) ' prt6 tried to read nonexistent file'
        stop                  
      endif
      nsize=ije(nij)  
      write (6,99) tit,nsize,iii,iiv,nij,(ij(i),ije(i),i=1,nij)
99    format(' file ',a4,', size iii iiv=',4i9,' ij records:'/
     1 (3(2x,z8,i8,';')))
      istrt1=nstrt(1,num-1)+2*nij+4
      istrt2=nstrt(2,num-1)
      do 40 i=1,nij
        iflag=1
        io=0
        do 30 j=ije(i-1)+1,ije(i)
          if (abs(gen(istrt2+j)).ge.flim) then
            if (iflag.eq.1) then
              id(1)=iand(ishft(ij(i),-16),'ff'x)
              id(2)=iand(ishft(ij(i),-8),'ff'x)
              id(3)=iand(ij(i),'ff'x)
c             write (6,*) ' ' 
              write (6,98) i,(id(k),k=1,3)
98            format ('record no. ',i6,' is ',3i2)
c             write (6,*) ' '
              iflag=0
            end if
            io=io+1
            a(io)=gen(istrt2+j)
            ii=igen(istrt1+j)
            ind(1,io)=iand(ishft(ii,-16),'ff'x)
            ind(2,io)=iand(ishft(ii,-8),'ff'x)
            ind(3,io)=iand(ii,'ff'x)
            if (io.eq.4) then
              write (6,97) ((ind(m,n),m=1,3),a(n),n=1,4)
              io=0
            endif
          endif
30      continue
        if (io.gt.0) then
          write (6,97) ((ind(m,n),m=1,3),a(n),n=1,io)
        end if
97      format (1x,4(1x,3i2,'=',f11.7))
40    continue
      write (6,*) ' ' 
      return
      end 
      subroutine prt8 (name,tit,iparm)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij,ij(ndij),ije(0:ndij),iii,iiv
      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      dimension a(4),id(4),ind(4,4)
      if (iparm.gt.iopt(1)) return
      flim=10.d0**(-iopt(5))
      if (iopt(5).le.0) flim=0.d0
      call files_ccsd (name,num,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) ' prt8 tried to read nonexistent file'
        stop                  
      endif
      nsize=ije(nij)  
      write (6,99) tit,nsize,iii,iiv,nij,(ij(i),ije(i),i=1,nij)
99    format(' file ',a4,', size iii iiv=',4i9,' ij records:'/
     1 (3(2x,z8,i8,';')))
      istrt1=nstrt(1,num-1)+2*nij+4
      istrt2=nstrt(2,num-1)
      do 40 i=1,nij
        iflag=1
        io=0
        do 30 j=ije(i-1)+1,ije(i)
          if (abs(gen(istrt2+j)).ge.flim) then
            if (iflag.eq.1) then
              id(1)=iand(ishft(ij(i),-24),'ff'x)
              id(2)=iand(ishft(ij(i),-16),'ff'x)
              id(3)=iand(ishft(ij(i),-8),'ff'x)
              id(4)=iand(ij(i),'ff'x)
              write (6,*) ' ' 
              write (6,98) i,(id(k),k=1,4)
98            format ('record no. ',i6,' is ',4i2)
              write (6,*) ' '
              iflag=0
            end if
            io=io+1
            a(io)=gen(istrt2+j)
            ii=igen(istrt1+j)
            ind(1,io)=iand(ishft(ii,-24),'ff'x)
            ind(2,io)=iand(ishft(ii,-16),'ff'x)
            ind(3,io)=iand(ishft(ii,-8),'ff'x)
            ind(4,io)=iand(ii,'ff'x)
            if (io.eq.4) then
              write (6,97) ((ind(m,n),m=1,4),a(n),n=1,4)
              io=0
            endif
          endif
30      continue
        if (io.gt.0) then
          write (6,97) ((ind(m,n),m=1,4),a(n),n=1,io)
        end if
97      format (4(1x,4i2,'=',f10.7))
40    continue
      write (6,*) ' ' 
      return
      end 

      subroutine prtime(tit,n,iparm)
      implicit real*8 (a-h,o-z)
      common /csdtime/ csdstim,csdltim
      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      character*8 clockx
      real*4 cputime
      if (iparm.gt.iopt(1)) return
      t = mclock()/1000000.
      diff = t - csdltim
      acum = t - csdstim

!mi: gfortran does not know time...
!     call time(clockx)

      print 900,tit,acum,diff,n,clockx
900   format (' time aft ',a8,f10.3,';added ',
     1 f8.3,' n=',i8,'. ',a24)
      csdltim = t
      call flush (6)
      return
      end
! reading routine 
      subroutine rdin
      implicit real*8 (a-h,o-z)
      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      character typ(2,5),cc(200)
      character*8 title(10)
      data typ /'c','C','h','H','p','P','g','G','v','V'/
      read (5,'(10a8)') title
      write(6,'(1h1,10a8)') title
      read (5,'(40i3)') (iopt(i),i=1,40)
      write(6,'(1x,40i3)') (iopt(i),i=1,40)
      read (5,'(i3/(200a1))') norb,(cc(i),i=1,norb)
      do 110 i=1,norb
      do 105 j=1,5
      if ((cc(i).ne.typ(1,j)).and.(cc(i).ne.typ(2,j))) go to 105
      ic(i)=j
      icv(i)=0
      if (j.le.3) go to 110
      ic(i)=j-2
      icv(i)=1
      go to 110
105        continue
      write (6,1050) i,cc(i)
1050  format ('0 **** error in type of orb',i3,3x,a1)
      call abext('rdin    ')
110   continue
      write (6,'(i4,'' norbs. types''/(1x,80a1))') norb,(cc(i),i=1,norb)
      return
      end


      subroutine rdfile(nt,nrec,name,ieof)
      implicit real*8 (a-h,o-z)
c  read file nt, skipping nrec records, to internal file number
c  nflcur+1, calling it "name".
c  disk file structure for tc etc:
c  1. nij (# of subblocks)
c  2. ij1(1:nij), subblock id's; ije1(0:nij), subblock ends; iii, iiv.
c  3. ia(1:ije1(nij)), indices.
c  4. a(1:ije1(nij)), actual values.
c     ***   april 1990: 1-3 are now in igen, 4 in gen.  ****
c    nt.eq.9 = file without indices (pppp).
c  // contains all files. file namefl(i) starts at nstrt(i-1)+1 (see /cfiles/).
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)


      common // gen(nbig)
      common /i/ igen(nbig)
      rewind nt
      if (nrec.gt.0) then
        do 10 i=1,4*nrec
10      read (nt)
      endif
      read (nt,end=190) nij
      nflcur=nflcur+1
      nstr1=nstrt(1,nflcur-1)
      nstr2=nstrt(2,nflcur-1)
      igen(nstr1+1)=nij
      read (nt) (igen(nstr1+i),i=2,2*nij+4)
      iiia(nflcur)=igen(nstr1+2*nij+3)
      iiva(nflcur)=igen(nstr1+2*nij+4)
      namefl(nflcur)=name
      nsize=igen(nstr1+2*nij+2)
      if (nt.ne.9) then
        read (nt) (igen(nstr1+2*nij+4+i),i=1,nsize)
        read (nt) (gen(nstr2+i),i=1,nsize)
        nstrt(1,nflcur)=nstr1+2*nij+4+nsize
        nstrt(2,nflcur)=nstr2+nsize
      else
        read (nt) (gen(nstr2+i),i=1,nsize)
        nstrt(1,nflcur)=nstr1+2*nij+4
        nstrt(2,nflcur)=nstr2+nsize
      endif
      call dimx(nstrt(1,nflcur),nbig,"rdfile1")
      call dimx(nstrt(2,nflcur),nbig,"rdfile2")
      return
190   if (ieof.eq.0) then
        ieof=1
        return
      else
        write (*,*) 'rdfile encountered eof'
        write (*,110) nt,nrec,name
110     format(' nt nrec name=',i5,i4,2x,a8)
        call abext('rdfile  ')
      endif
      end


      subroutine rdsym
!
! reading of the input file ...
!
      implicit real*8 (a-h,o-z)
      common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
      common /dpd/ nir,igmt(8,8),ici(200)
      dimension itemp(64)
      data itemp / 1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,
     1           4,3,2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,
     2           7,8,5,6,3,4,1,2,8,7,6,5,4,3,2,1 /

      write(*,*) 'Enter nclass variable:'
      read (*,1000,iostat=ios) nclass
      if (ios.ne.0) then
        write(*,*) 'error in reading nclass!'
        stop 'error in reading nclass!'
      endif
      read (*,1000,iostat=ios) (iclass(i),i=1,nclass)
      if (ios.ne.0) then
        write(*,*) 'error in reading iclass!'
        stop 'error in reading iclass!'
      endif

1000  format(20i3)
      write (*,1100) nclass,(iclass(i),i=1,nclass)
1100  format(i4,' orb classes. # of orbs ',30i3)
c-------------------------------------------------
c  isym is  0 for orbs remaining invariant
c  isym is  1 for orbs changing sign
c  isym is -1 for orbs doing neither
c-------------------------------------------------
      nsym=1
110   read (*,1000) (isym(i,nsym),i=1,nclass)
      if (isym(1,nsym).ge.9) go to 120
      write (*,1200) nsym,(isym(i,nsym),i=1,nclass)
1200  format(' isym(',i3,')=',30i3)
      nsym=nsym+1
      go to 110
120   nsym=nsym-1
c
      ibias(1)=0
      do 140 i=1,nclass
140   ibias(i+1)=ibias(i)+iclass(i)
      do 20 i=1,nclass
        do 30 j=ibias(i)+1,ibias(i+1)
          ici(j)=i
30      continue
20    continue
      ii=0
      do 40 i=1,8
        do 50 j=1,8
          ii=ii+1
          igmt(i,j)=itemp(ii)
50      continue
40    continue
!  nir should be the number of irreducible representations ...
      if (nclass.eq.1) then
        nir=1
      else
        nir=2**nsym
      end if
  
      write(*,*) 'returning from rdsym'

      return
      end


      subroutine rdrsp1
      implicit real*8 (a-h,o-z)
      common /rsp1/  nrsp1,irsp1(100)
      dimension it(2,100)
      read (5,*) nrsp1,((it(j,i),j=1,2),i=1,nrsp1)
      call dimx(nrsp1,100,"rdrsp1")
      do 10 i=1,nrsp1
10    irsp1(i)=ior(ishft(it(1,i),8),it(2,i))
      write (6,100) nrsp1,(it(1,i),it(2,i),i=1,nrsp1)
100   format(i4,' p to r transitions: ',30(2i3,1x))
      return
      end


      subroutine rdrsp2
      implicit real*8 (a-h,o-z)
      common /rsp2/  nrsp2,irsp2(500)
      dimension it(2,500)
      read (5,*) nrsp2,((it(j,i),j=1,2),i=1,nrsp2)
      call dimx(nrsp2,500,"rdrsp2")
      do 10 i=1,nrsp2
10    irsp2(i)=ior(ishft(it(1,i),8),it(2,i))
      write (6,100) nrsp2,(it(1,i),it(2,i),i=1,nrsp2)
100   format(i4,' r space dets: ',30(2i3,1x))
      return
      end



      subroutine reads (nt,nrec,name)       
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      rewind nt
      igotit=0
      do 20 i=1,nflcur
        if (namefl(i).eq.name) then
          nfl=i
          igotit=1
        end if
20    continue
      if (igotit.eq.0) then
        nflcur=nflcur+1
        nfl=nflcur
      end if
      if (nrec.gt.0) then
        do 10 i=1,4*nrec
10      read (nt)
      endif
      read (nt,iostat=ieof) nij
      if (ieof.lt.0) then
        write (6,*) 'reads encountered end of file'
        stop
      end if
      istrt1=nstrt(1,nfl-1)
      istrt2=nstrt(2,nfl-1)
      igen(istrt1+1)=nij
      read (nt) (igen(istrt1+i),i=2,2*nij+4)
      iiia(nfl)=igen(istrt1+2*nij+3)
      iiva(nfl)=igen(istrt1+2*nij+4)
      namefl(nfl)=name
      nsize=igen(istrt1+2*nij+2)
      read (nt) (igen(istrt1+2*nij+4+i),i=1,nsize)
      read (nt) (gen(istrt2+i),i=1,nsize)
      nstrt(1,nfl)=istrt1+2*ndij+4+ndp3 
      nstrt(2,nfl)=istrt2+ndp3 
      call dimx(nstrt(1,nfl),nbig,"reads   ")
      call dimx(nstrt(2,nfl),nbig,"reads   ")
      return
      end
!
      subroutine readv(nint,v,ind,nintmx,nfock,fock,ifock,nfockmx)
!
!         basic reading subroutine called from many places
!
      implicit real*8 (a-h,o-z)
      dimension v(*),ind(*),irerdr(200),irerdrv(200),etmp(200)
      dimension fock(*),ifock(*)

      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
      common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
      dimension ilet(6)
      data ilet /1hc,1hc,1hh,1hg,1hp,1hv/
      data zero /1.d-8/

      write(*,*) 'entered readv'

c  iopt(19) is 0 for GAMESS, 1 for ACES2(4.92), 2 for ACES2(99)
      if (iopt(19).eq.0) then
        intrc=2725
        write(*,*) 'this is GAMESS connected stuff !'
      endif

      if (iopt(19).ne.0) intrc=600

      rewind(12,iostat=ios)
      if (ios.ne.0) then
       write(*,*) 'failed rewing of unit 12 !'
       call  abext('readv : failed rewing of unit 12 !')
      endif

      read (12,iostat=ios) nort,(eps(i),i=1,nort)
      if (ios.ne.0) then
        write(*,*) 'failed read (12) nort,(eps(i),i=1,nort)!'
        call  abext('readv : failed read (12) nort,(eps(i),i=1,nort)')
      endif

      if (nort.ne.norb) write (*,*) 'discrepancy in norb: ',norb,nort
      if (nort.ne.norb) then                      
        call abext('readv: nort.ne.norb ! ')
      end if

c  changes for aces99
      if (iopt(19).eq.2) then
        read (12) (irerdr(i),i=1,norb)
        do i=1,norb
          irerdrv(irerdr(i))=i
          etmp(i)=eps(irerdr(i))
        enddo
        do i=1,norb
          eps(i)=etmp(i)
        enddo
      endif

      do 340 icl=1,nclass
      write (*,3300) icl
3300  format(' irrep ',i2)
      write (*,3400) (io,ilet(2*ic(io)-1+icv(io)),eps(io),
     1 io=ibias(icl)+1,ibias(icl+1))
3400  format(3('orb ',i3,2x,a1,', e=',1pe12.5,:,';  '))
340   continue

      if (iopt(19).eq.0) read (12)
      nint=0
120   if(iopt(19).eq.0) then
        read (12) nn,(ind(nint+i),i=1,2725),
     1  (v(nint+i),i=1,2725)
      else
c       write(*,*) 'bef read'
        read (12,end=135,err=199) (v(nint+i),i=1,600),(ind(nint+i),
     1  i=1,600),nn
c       write(*,*) 'aft read. nn=',nn
        if (nn.le.0) go to 135
      endif
      call dimx(nint+intrc,nintmx,"nint")

C--------------------------------------------------------
C         UNPACK, ELIMINATE INTS WITH CORE ORBS
C--------------------------------------------------------
      nnint=nint
      do 130 i=nnint+1,nnint+abs(nn)
      if (abs(v(i)).lt.zero) go to 130
      J1=iand(ishft(ind(i),-24),'ff'X)
      J2=iand(ishft(ind(i),-16),'ff'X)
      J3=iand(ishft(ind(i),-8) ,'ff'X)
      J4=iand(ind(i)           ,'ff'X)
      if(iopt(19).eq.2) then
        j1=irerdrv(j1)
        j2=irerdrv(j2)
        j3=irerdrv(j3)
        j4=irerdrv(j4)
      endif
      if((ic(j1).eq.1).or.(ic(j2).eq.1).or.(ic(j3).eq.1).or.
     1                                      (ic(j4).eq.1)) go to 130
      nint=nint+1
      v(nint)=v(i)
      ind(nint)=ior(ishft(j1,24),ior(ishft(j3,16),ior(ishft(j2,8),j4)))
130   continue
      if (nn.gt.0) go to 120
199   write(*,*) 'Record in tape12 too short'
      call abext('12short')
135   continue
      write (6,*) nint,' integrals and indices read in '
      if (iopt(8).eq.0) then
        nfock=0
        return
      endif
      if (iopt(8).lt.0) go to 410
! ----------- Fock matrix reading --------------
      READ (12) NFOCK,(FOCK(I),I=1,NFOCK),(IFOCK(I),I=1,NFOCK)
      WRITE (6,*) NFOCK,' FOCK MATRIX ELEMENTS READ IN'
      CALL dimx(NFOCK,NFOCKMX,"NFOCK")
      lf=0
      do 220 i=1,nfock
      if (abs(fock(i)).lt.1.d-9) go to 220
      lf=lf+1
      fock(lf)=fock(i)
      ifock(lf)=ifock(i)
220   continue
      nfock=lf
      write (*,*) nfock,' nonzero fock matrix elements left'
      do 250 i=1,nfock
      i1=iand(ifock(i),'ff'X)
      i2=ishft(ifock(i),-8)
      if (i1.ne.i2) go to 250
      if (iopt(8).eq.2) eps(i1)=fock(i)
      fock(i)=fock(i)-eps(i1)
250   continue
      do 260 i=1,nfock
      i1=iand(ifock(i),'ff'X)
      i2=ishft(ifock(i),-8)
C      write (*,2600) i,i1,i2,fock(i)
C 2600  format(' Fock element no',i3,2x,2i3,f15.6)
260   continue
      go to 440
410   adde=-iopt(8)*0.1
      do 430 io=1,norb
        if (ic(io).eq.2.and.icv(io).eq.0) then
          eps(io)=eps(io)-2.*adde
          nfock=nfock+1
          ifock(nfock)=ior(ishft(io,8),io)
          fock(nfock)=2.*adde
        endif
        if (ic(io).eq.3.and.icv(io).eq.0) then
          eps(io)=eps(io)+2*adde
          nfock=nfock+1
          ifock(nfock)=ior(ishft(io,8),io)
          fock(nfock)=-2*adde
        endif
        if (ic(io).eq.2.and.icv(io).eq.1) then
          eps(io)=eps(io)-adde
          nfock=nfock+1
          ifock(nfock)=ior(ishft(io,8),io)
          fock(nfock)=adde
        endif
        if (ic(io).eq.3.and.icv(io).eq.1) then
          eps(io)=eps(io)+adde
          nfock=nfock+1
          ifock(nfock)=ior(ishft(io,8),io)
          fock(nfock)=-adde
        endif
430   continue
440   write (*,*) ' modified orbital energies'
      do 450 icl=1,nclass
        write (*,3300) icl
        write (*,3400) (io,ilet(2*ic(io)-1+icv(io)),eps(io),
     1  io=ibias(icl)+1,ibias(icl+1))
450   continue
      return
      end

!      subroutine files(name,num,new,ijnum,iii,iiv)
      subroutine files_ccsd(name,num,new,ijnum,iii,iiv)
      implicit real*8 (a-h,o-z)
!     integer*4 ch1,ch2
      integer ch1,ch2
      
#include "param_ccsd.h"

      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
      common /cij4/ nij4,ij4(ndij),ije4(0:ndij),iii4,iiv4
      common /cij5/ nij5,ij5(ndij),ije5(0:ndij),iii5,iiv5
      common /cij6/ nij6,ij6(ndij),ije6(0:ndij),iii6,iiv6

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
c   subroutine finds # of file named NAME. if not found,
c   increases nflcur by 1, starts new file.
c   nflcur set to 0 at start of program.
c   new is 1 for new file, -1 for new file with indices copied from
c   another file, 0 for existing file. ijnum=0 used for existing
c   scratch file, i.e. file to be completely rewritten.
      if (nflcur.eq.0) go to 20
      num=1000
      do 10 i=1,nflcur
10    if (name.eq.namefl(i)) num=i
      if (num.le.nflcur) then
        new=0
c  next lines for case where temp file reopened with same name,
c  different type (must be sn files). iii=-1 means ignore iii.

	if (iii.gt.-1.and.((iii.ne.iiia(num)).or.(iiv.ne.iiva(num)))) then
C Gal 95
           itmp=4H       
           jtmp='ffffff00'X
           ch1 = iand(itmp,jtmp)
           ch2 = iand (namefl(num), '000000ff'X)
           isaa = ior (ch1, ch2)

          i_4hs=4hs 
!         if (isaa.ne.4hs   ) then
          if (isaa.ne.i_4hs) then
            write (*,120) name,iii,iiv,iiia(num),iiva(num)
120         format(' file ',a8,' with iii iiv',2i5,' previously',2i5)
            call abext('files   ')
          endif
        endif
      if (iii.gt.-1.and.((iii.ne.iiia(num)).or.(iiv.ne.iiva(num)))) then
          new=1
          iiia(num)=iii
          iiva(num)=iiv
          do 15 i=1,nflcur
         if ((i.ne.num).and.(iiia(i).eq.iii).and.(iiva(i).eq.iiv)) then
            new=-1
            call imove(i,num)
            return
          endif
15        continue
          return
        endif
        if (ijnum.eq.0) return
        nija=igen(nstrt(1,num-1)+1)
        if (ijnum.eq.1) then
          call kcopy(nija+1,igen(nstrt(1,num-1)+1),nij1)
          call kcopy(nija+1,igen(nstrt(1,num-1)+nija+2),ije1(0))
          iii1=igen(nstrt(1,num-1)+2*nija+3)
          iiv1=igen(nstrt(1,num-1)+2*nija+4)
        elseif (ijnum.eq.2) then
          call kcopy(nija+1,igen(nstrt(1,num-1)+1),nij2)
          call kcopy(nija+1,igen(nstrt(1,num-1)+nija+2),ije2(0))
          iii2=igen(nstrt(1,num-1)+2*nija+3)
          iiv2=igen(nstrt(1,num-1)+2*nija+4)
        elseif (ijnum.eq.3) then
          call kcopy(nija+1,igen(nstrt(1,num-1)+1),nij3)
          call kcopy(nija+1,igen(nstrt(1,num-1)+nija+2),ije3(0))
          iii3=igen(nstrt(1,num-1)+2*nija+3)
          iiv3=igen(nstrt(1,num-1)+2*nija+4)
        elseif (ijnum.eq.4) then
          call kcopy(nija+1,igen(nstrt(1,num-1)+1),nij4)
          call kcopy(nija+1,igen(nstrt(1,num-1)+nija+2),ije4(0))
          iii4=igen(nstrt(1,num-1)+2*nija+3)
          iiv4=igen(nstrt(1,num-1)+2*nija+4)
        elseif (ijnum.eq.5) then
          call kcopy(nija+1,igen(nstrt(1,num-1)+1),nij5)
          call kcopy(nija+1,igen(nstrt(1,num-1)+nija+2),ije5(0))
          iii5=igen(nstrt(1,num-1)+2*nija+3)
          iiv5=igen(nstrt(1,num-1)+2*nija+4)
        elseif (ijnum.eq.6) then
          call kcopy(nija+1,igen(nstrt(1,num-1)+1),nij6)
          call kcopy(nija+1,igen(nstrt(1,num-1)+nija+2),ije6(0))
          iii6=igen(nstrt(1,num-1)+2*nija+3)
          iiv6=igen(nstrt(1,num-1)+2*nija+4)
        else
          write (6,*) ' files called with incorrect ijnum, ',ijnum
          call abext('files   ')  
        endif
        return
      endif
c  new file
20    nflcur=nflcur+1
      namefl(nflcur)=name
      num=nflcur
      iiia(nflcur)=iii
      iiva(nflcur)=iiv
      new=1
      if (iii.eq.-2) return
      if (nflcur.eq.1) go to 40
      do 30 i=1,nflcur-1
      if ((iiia(i).eq.iii).and.(iiva(i).eq.iiv)) then
        new=-1
        call imove(i,num)
c     write (6,*) ' transf indx from ',i,' to ',num
        go to 40
      endif
30    continue
40    continue
      return
      end

      subroutine reorder(namef,namet,iord,iparm)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /ctmp/ tmp(ndp4),itmp(ndp4)

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      dimension j(4),msk(4)
      data msk /'ff000000'X,'ff0000'X,'ff00'X,'ff'X/
c   reorder from file namef (info in /cij2/) to namet (/cij1/)
      ipack(i1,i2,i3,i4)=ior(ior(ishft(i1,24),ishft(i2,16)),
     1                       ior(ishft(i3,8) ,i4)         )
      call files_ccsd(namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (6,10) namef
10      format(' reorder input tape ',a6,' does not exist.')
        call abext('reorder ')
      endif
      call parse(iord,i1,i2,i3,i4)
      call parse(iii2,j(1),j(2),j(3),j(4))
      iii1=1000*j(i1)+100*j(i2)+10*j(i3)+j(i4)
      call parse(iiv2,j(1),j(2),j(3),j(4))
      iiv1=1000*j(i1)+100*j(i2)+10*j(i3)+j(i4)
      call tmplt(iii1,iiv1,namet,ntapet)
      if (ije2(nij2).eq.0) return
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      istrt21=nstrt(1,ntapef-1)+2*nij2+4
      istrt22=nstrt(2,ntapef-1)
      lngth1=ije1(nij1)
      lngth2=ije2(nij2)
      call vands(lngth2,igen(istrt21+1),msk(i1),itmp)
      call vshfcs(lngth2,itmp,8*i1-8,itmp)
      call vands(lngth2,igen(istrt21+1),msk(i2),itmp(lngth2+600))
      call vshfcs(lngth2,itmp(lngth2+600),8*i2-16,itmp(lngth2+600))
      call viorv(lngth2,itmp,itmp(lngth2+600),itmp)
      call vands(lngth2,igen(istrt21+1),msk(i3),itmp(lngth2+600))
      call vshfcs(lngth2,itmp(lngth2+600),8*i3-24,itmp(lngth2+600))
      call viorv(lngth2,itmp,itmp(lngth2+600),itmp)
      call vands(lngth2,igen(istrt21+1),msk(i4),itmp(lngth2+600))
      call vshfcs(lngth2,itmp(lngth2+600),8*i4-32,itmp(lngth2+600))
      call viorv(lngth2,itmp,itmp(lngth2+600),itmp)
      call dcopy(lngth2,gen(istrt22+1),1,gen(istrt12+1),1)
      call vxsort(itmp,lngth2,lngth1,gen(istrt12+1),
     1 igen(istrt11+1))
      call prt(namet,namet,iparm)
      call prtime(namet,ije1(nij1),1)
      return
      end


      subroutine reordx (namef,namet,iord)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
      common /ctmp/ tmp(ndp4),itmp(ndp4)
      common /dpd/ nir,igmt(8,8),ici(200)
      dimension ia(8),ja(8)
      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
      do 500 i=1,8
        ja(i)=0
500   continue
      call parse8 (iord,ia(1),ia(2),ia(3),ia(4),ia(5),ia(6),ia(7),ia(8))
      do 510 i=1,8
        do 520 j=1,8
          if (ia(i).eq.j) then
            ja(j)=ja(j)+1
          end if
520     continue
510   continue
      do 530 i=1,8
        if (ja(i).ne.1) then
          write (6,*) 'reordx called with bad value of iord'
          stop
        end if
530   continue
      call files_ccsd (namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'reordx input tape does not exist, ',name2
        stop                      
      end if
      call parse8 (iord,i1,i2,i3,i4,i5,i6,i7,i8)
      call parse8 (iii2,ia(1),ia(2),ia(3),ia(4),ia(5),ia(6),ia(7),ia(8))
      call parse8 (iiv2,ja(1),ja(2),ja(3),ja(4),ja(5),ja(6),ja(7),ja(8))
      iii1=10000000*ia(i1)+1000000*ia(i2)+100000*ia(i3)+10000*ia(i4)+
     1     1000*ia(i5)+100*ia(i6)+10*ia(i7)+ia(i8)
      iiv1=10000000*ja(i1)+1000000*ja(i2)+100000*ja(i3)+10000*ja(i4)+
     1     1000*ja(i5)+100*ja(i6)+10*ja(i7)+ja(i8)
      call tmplt8 (iii1,iiv1,namet,ntapet)
      istrt12=nstrt(2,ntapet-1)
      istrt21=nstrt(1,ntapef-1)+2*nij2+4
      istrt22=nstrt(2,ntapef-1)
      do 10 i=1,nij2
        do 20 j=ije2(i-1)+1,ije2(i)
          itmp(j)=ij2(i)
          itmp(ndp3+j)=igen(istrt21+j)
          tmp(j)=gen(istrt22+j)
20      continue
10    continue
      nsize=ije2(nij2)
      do 30 i=1,nsize     
        do 40 j=0,3
          ia(4-j)=iand(ishft(itmp(i),-8*j),'ff'x)
          ia(8-j)=iand(ishft(itmp(ndp3+i),-8*j),'ff'x)
40      continue
        itmp(i)=ipack4(ia(i1),ia(i2),ia(i3),ia(i4))
        itmp(ndp3+i)=ipack4(ia(i5),ia(i6),ia(i7),ia(i8))
30    continue
      k=nsize/2+1
      ir=nsize
      iquit=0
      do while (iquit.eq.0)
        if (k.gt.1) then
          k=k-1
          idtmp=itmp(k)
          indtmp=itmp(ndp3+k)
          vtmp=tmp(k)
        else
          idtmp=itmp(ir)
          indtmp=itmp(ndp3+ir)
          vtmp=tmp(ir)
          itmp(ir)=itmp(1)
          itmp(ndp3+ir)=itmp(ndp3+1)
          tmp(ir)=tmp(1)
          ir=ir-1
          if (ir.eq.1) then
            itmp(1)=idtmp
            itmp(ndp3+1)=indtmp
            tmp(1)=vtmp
            iquit=1
          end if
        end if
        if (iquit.eq.0) then
          i=k
          j=2*k
          do while (j.le.ir)
            if (j.lt.ir) then
              ii=imax(itmp(j),itmp(ndp3+j),
     1        itmp(j+1),itmp(ndp3+j+1),1)
              if (ii.eq.2) then
                j=j+1
              end if
            end if
            ii=imax(idtmp,indtmp,itmp(j),itmp(ndp3+j),1)
            if (ii.eq.2) then         
              itmp(i)=itmp(j)
              itmp(ndp3+i)=itmp(ndp3+j)
              tmp(i)=tmp(j)
              i=j
              j=2*j
            else
              j=ir+1
            end if
          end do  
          itmp(i)=idtmp
          itmp(ndp3+i)=indtmp
          tmp(i)=vtmp   
        end if
      end do  
      call dcopy (nsize,tmp(1),1,gen(istrt12+1),1)
      return
      end


      subroutine reord6 (namef,namet,iord)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
      common /ctmp/ tmp(ndp4),itmp(ndp4)
      common /dpd/ nir,igmt(8,8),ici(200)
      dimension ia(6),ja(6),msk(4),igot(8)
      data msk /'ff000000'x,'ff0000'x,'ff00'x,'ff'x/
      ndp8=ndp3/nir
      call files_ccsd (namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'reord6 input tape does not exist, ',name2
        stop                      
      end if
      call parse6 (iord,idum1,idum2,i1,i2,i3,i4)
      if ((idum1.ne.1).or.(idum2.ne.2)) then
        write (6,*) 'reord6 called with incorrect iord'
        stop
      end if
      call parse6 (iii2,ia(1),ia(2),ia(3),ia(4),ia(5),ia(6))
      call parse6 (iiv2,ja(1),ja(2),ja(3),ja(4),ja(5),ja(6))
      iii1=100000*ia(1)+10000*ia(2)+1000*ia(i1)+
     1     100*ia(i2)+10*ia(i3)+ia(i4)
      iiv1=100000*ja(1)+10000*ja(2)+1000*ja(i1)+
     1     100*ja(i2)+10*ja(i3)+ja(i4)
      i1=i1-2
      i2=i2-2
      i3=i3-2
      i4=i4-2
      call tmplt6 (iii1,iiv1,namet,ntapet)
      istrtt2=nstrt(2,ntapet-1)
      istrtf1=nstrt(1,ntapef-1)+2*nij2+4
      istrtf2=nstrt(2,ntapef-1)
      call tchng6 (nij1,ij1,ije1,ntapet,2)
      call tchng6 (nij2,ij2,ije2,ntapef,2)
      do 100 i=1,nir
        igot(i)=0
100   continue
      do 10 i=1,nij2
        lngth=ije1(i)-ije1(i-1)
        call dimx (lngth,ndp8,"reord6  ")
        if (ije2(i)-ije2(i-1).ne.lngth) then
          write (6,*) 'mismatch in block lengths - reord6'
          stop                    
        end if
        ii=iand(ishft(ij2(i),-8),'ff'x)
        jj=iand(ij2(i),'ff'x)
        ii=igmt(ici(ii),ici(jj))
        if (igot(ii).eq.0) then
          igot(ii)=1
          istrt=istrtf1+ije2(i-1)+1
          call vands  (lngth,igen(istrt),msk(i1),itmp(1))
          call vshfcs (lngth,itmp(1),8*i1-8,itmp(1))
          call vands  (lngth,igen(istrt),msk(i2),itmp(ndp8))
          call vshfcs (lngth,itmp(ndp8),8*i2-16,itmp(ndp8))
          call viorv  (lngth,itmp(1),itmp(ndp8),itmp(1))
          call vands  (lngth,igen(istrt),msk(i3),itmp(ndp8))
          call vshfcs (lngth,itmp(ndp8),8*i3-24,itmp(ndp8))
          call viorv  (lngth,itmp(1),itmp(ndp8),itmp(1))
          call vands  (lngth,igen(istrt),msk(i4),itmp(ndp8))
          call vshfcs (lngth,itmp(ndp8),8*i4-32,itmp(ndp8))
          call viorv  (lngth,itmp(1),itmp(ndp8),itmp(1))
          do 60 j=1,lngth       
            tmp(j)=j+0.01d0
60        continue
          call vsort (itmp(1),lngth,tmp)
          do 70 j=1,lngth
           itmp(ndp3+(ii-1)*ndp8+j)=int(tmp(j))
70        continue
        end if
        do 80 j=1,lngth
          gen(istrtt2+ije1(i-1)+j)=
     1    gen(istrtf2+ije2(i-1)+itmp(ndp3+(ii-1)*ndp8+j))
80      continue
10    continue
      call prtime ("reord6   ",0,1)
      return
      end


      subroutine reord7 (namef,namet,iord)         
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /ctmp/ tmp(ndp4),itmp(ndp4)
      dimension ia(6),ja(6),msk(3)
      data msk / 'ff0000'x,'ff00'x,'ff'x /
      call files_ccsd (namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'reord7 input tape does not exist, ',name2
        stop                     
      end if
      call parse6 (iord,i1,i2,i3,idum1,idum2,idum3)
      call parse6 (iii2,ia(1),ia(2),ia(3),ia(4),ia(5),ia(6))
      call parse6 (iiv2,ja(1),ja(2),ja(3),ja(4),ja(5),ja(6))
      iii1=100000*ia(i1)+10000*ia(i2)+1000*ia(i3)+
     1     100*ia(4)+10*ia(5)+ia(6)
      iiv1=100000*ja(i1)+10000*ja(i2)+1000*ja(i3)+
     1     100*ja(4)+10*ja(5)+ja(6)
      call tmplt6 (iii1,iiv1,namet,ntapet)
      if (nij2.ne.nij1) then
        write (6,*) 'mismatch in numbers of blocks - reord7'
        stop
      end if
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      istrt21=nstrt(1,ntapef-1)+2*nij2+4
      istrt22=nstrt(2,ntapef-1)
      call vands  (nij2,ij2,msk(i1),itmp(1))
      call vshfcs (nij2,itmp(1),8*i1-8,itmp(1))
      call vands  (nij2,ij2,msk(i2),itmp(ndp3+1))
      call vshfcs (nij2,itmp(ndp3+1),8*i2-16,itmp(ndp3+1))
      call viorv  (nij2,itmp(1),itmp(ndp3+1),itmp(1))
      call vands  (nij2,ij2,msk(i3),itmp(ndp3+1))
      call vshfcs (nij2,itmp(ndp3+1),8*i3-24,itmp(ndp3+1))
      call viorv  (nij2,itmp(1),itmp(ndp3+1),itmp(1))
      do 10 i=1,nij2        
        tmp(i)=i+0.01d0
10    continue
      call vsort (itmp(1),nij2,tmp(1))
      do 20 i=1,nij2
        itmp(i)=tmp(i)     
20    continue
      do 30 i=1,nij2
        lngth=ije2(itmp(i))-ije2(itmp(i)-1)
        call dcopy (lngth,gen(istrt22+ije2(itmp(i)-1)+1),1,
     1  gen(istrt12+ije1(i-1)+1),1)
30    continue
      call prtime ("reord7  ",0,1)
      return
      end

      subroutine reord10 (namef,namet,iord)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

      common // gen(nbig)
      common /i/ igen(nbig)
      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
      common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
      common /ijtemp/ ijtmp1(ndp6)
      common /ctmp/ tmp(ndp4),itmp(ndp4)
      dimension ia(10),ja(10),msk(5),iord(10)
      data msk / 'f0000'x,'f000'x,'f00'x,'f0'x,'f'x /
      call files_ccsd (namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'reord10 input tape does not exist, ',name2
        stop                       
      end if
      i1=iord(1)
      i2=iord(2)
      i3=iord(3)
      i4=iord(4)
      i5=iord(5)
      call tmplt10 (namet,ntapet)
      istrtt1=nstrt(1,ntapet-1)+2*nij1+4
      istrtt2=nstrt(2,ntapet-1)
      istrtf1=nstrt(1,ntapef-1)+2*nij2+4
      istrtf2=nstrt(2,ntapef-1)
      do 10 i=1,nij2
        lngth=ije1(i)-ije1(i-1)
        call dimx (lngth,ndp5,"reord10 ")
        if (ije2(i)-ije2(i-1).ne.lngth) then
          write (6,*) 'mismatch in block lengths - reord10'
          stop                        
        end if
        do 50 j=1,lngth
          ijtmp1(j)=igen(istrtf1+ije2(i-1)+j)
50      continue
        call vands  (lngth,ijtmp1(1),msk(i1),ijtmp1(ndp5))
        call vshfcs (lngth,ijtmp1(ndp5),4*i1-4,ijtmp1(ndp5))
        call vands  (lngth,ijtmp1(1),msk(i2),ijtmp1(2*ndp5))
        call vshfcs (lngth,ijtmp1(2*ndp5),4*i2-8,ijtmp1(2*ndp5))
        call viorv  (lngth,ijtmp1(ndp5),ijtmp1(2*ndp5),ijtmp1(ndp5))
        call vands  (lngth,ijtmp1(1),msk(i3),ijtmp1(2*ndp5))
        call vshfcs (lngth,ijtmp1(2*ndp5),4*i3-12,ijtmp1(2*ndp5))
        call viorv  (lngth,ijtmp1(ndp5),ijtmp1(2*ndp5),ijtmp1(ndp5))
        call vands  (lngth,ijtmp1(1),msk(i4),ijtmp1(2*ndp5))
        call vshfcs (lngth,ijtmp1(2*ndp5),4*i4-16,ijtmp1(2*ndp5))
        call viorv  (lngth,ijtmp1(ndp5),ijtmp1(2*ndp5),ijtmp1(ndp5))
        call vands  (lngth,ijtmp1(1),msk(i5),ijtmp1(2*ndp5))
        call vshfcs (lngth,ijtmp1(2*ndp5),4*i5-20,ijtmp1(2*ndp5))
        call viorv  (lngth,ijtmp1(ndp5),ijtmp1(2*ndp5),ijtmp1(ndp5))
        call dcopy  (lngth,gen(istrtf2+ije2(i-1)+1),1,tmp,1)
        call vsort  (ijtmp1(ndp5),lngth,tmp)
        call dcopy  (lngth,tmp,1,gen(istrtt2+ije1(i-1)+1),1)
10    continue
      call prtime ("reord10  ",0,1)
      return
      end


      subroutine reord11 (namef,namet,iord)         
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)

      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2

#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)

      dimension ia(10),ja(10),msk(5),indx(ndij),dummy(ndij),iord(10)
      data msk / 'f0000'x,'f000'x,'f00'x,'f0'x,'f'x /

      i1=iord(1)
      i2=iord(2)
      i3=iord(3)
      i4=iord(4)
      i5=iord(5)
      call tmplt10 (namet,ntapet)
      call files_ccsd (namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'reorder input tape does not exist, ',name2
        stop                      
      end if
      if (nij2.ne.nij1) then
        write (6,*) 'mismatch in numbers of blocks - reord11'
        stop                      
      end if
      istrtt1=nstrt(1,ntapet-1)+2*nij1+4
      istrtt2=nstrt(2,ntapet-1)
      istrtf1=nstrt(1,ntapef-1)+2*nij2+4
      istrtf2=nstrt(2,ntapef-1)
      call vands  (nij2,ij2,msk(i1),itmp)
      call vshfcs (nij2,itmp,4*i1-4,itmp)
      call vands  (nij2,ij2,msk(i2),itmp(ndp5))
      call vshfcs (nij2,itmp(ndp5),4*i2-8,itmp(ndp5))
      call viorv  (nij2,itmp,itmp(ndp5),itmp)
      call vands  (nij2,ij2,msk(i3),itmp(ndp5))
      call vshfcs (nij2,itmp(ndp5),4*i3-12,itmp(ndp5))
      call viorv  (nij2,itmp,itmp(ndp5),itmp)
      call vands  (nij2,ij2,msk(i4),itmp(ndp5))
      call vshfcs (nij2,itmp(ndp5),4*i4-16,itmp(ndp5))
      call viorv  (nij2,itmp,itmp(ndp5),itmp)
      call vands  (nij2,ij2,msk(i5),itmp(ndp5))
      call vshfcs (nij2,itmp(ndp5),4*i5-20,itmp(ndp5))
      call viorv  (nij2,itmp,itmp(ndp5),itmp)
      do 30 i=1,nij2        
        dummy(i)=i*1.d0
30    continue
      call vsort (itmp,nij2,dummy)
      do 40 i=1,nij2
        d=dummy(i)-int(dummy(i))
        if (d.gt.0.5d0) then
          id=int(dummy(i))+1
        else   
          id=int(dummy(i))
        end if
        indx(i)=id
40    continue
      k=1
      do 10 i=1,nij2
        n1=ije1(i)-ije1(i-1)
        n2=ije2(indx(i))-ije2(indx(i)-1)
        if (n1.ne.n2) then
          write (6,*) 'mismatch in block lengths - reord11'
        end if
        do 20 j=ije2(indx(i)-1)+1,ije2(indx(i))
          tmp(k)=gen(istrtf2+j)
          k=k+1
20      continue
10    continue
      call dcopy (k,tmp,1,gen(istrtt2+1),1)   
      call prtime ("reord11 ",0,2)
      return
      end


      subroutine reord12 (namef,namet,iord)
      implicit real*8 (a-h,o-z)
       
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)

#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2

#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)

      common /ijtemp/ ijtmp1(ndp6)

#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)

      dimension ia(12),ja(12),msk(6),iord(12)
      data msk / 'f00000'x,'f0000'x,'f000'x,'f00'x,'f0'x,'f'x /
      call files_ccsd (namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'reord12 input tape does not exist, ',name2
        stop                     
      end if
      i1=iord(1)
      i2=iord(2)
      i3=iord(3)
      i4=iord(4)
      i5=iord(5)
      i6=iord(6)
      call tmplt12 (namet,ntapet)
      istrtt1=nstrt(1,ntapet-1)+2*nij1+4
      istrtt2=nstrt(2,ntapet-1)
      istrtf1=nstrt(1,ntapef-1)+2*nij2+4
      istrtf2=nstrt(2,ntapef-1)
      do 10 i=1,nij2
        lngth=ije1(i)-ije1(i-1)
        call dimx (lngth,ndp5,"reord12 ")
        if (ije2(i)-ije2(i-1).ne.lngth) then
          write (6,*) 'mismatch in block lengths - reord12'
          stop                         
        end if
        do 50 j=1,lngth
          ijtmp1(j)=igen(istrtf1+ije2(i-1)+j)
50      continue
        call vands  (lngth,ijtmp1(1),msk(i1),ijtmp1(ndp5))
        call vshfcs (lngth,ijtmp1(ndp5),4*i1-4,ijtmp1(ndp5))
        call vands  (lngth,ijtmp1(1),msk(i2),ijtmp1(2*ndp5))
        call vshfcs (lngth,ijtmp1(2*ndp5),4*i2-8,ijtmp1(2*ndp5))
        call viorv  (lngth,ijtmp1(ndp5),ijtmp1(2*ndp5),ijtmp1(ndp5))
        call vands  (lngth,ijtmp1(1),msk(i3),ijtmp1(2*ndp5))
        call vshfcs (lngth,ijtmp1(2*ndp5),4*i3-12,ijtmp1(2*ndp5))
        call viorv  (lngth,ijtmp1(ndp5),ijtmp1(2*ndp5),ijtmp1(ndp5))
        call vands  (lngth,ijtmp1(1),msk(i4),ijtmp1(2*ndp5))
        call vshfcs (lngth,ijtmp1(2*ndp5),4*i4-16,ijtmp1(2*ndp5))
        call viorv  (lngth,ijtmp1(ndp5),ijtmp1(2*ndp5),ijtmp1(ndp5))
        call vands  (lngth,ijtmp1(1),msk(i5),ijtmp1(2*ndp5))
        call vshfcs (lngth,ijtmp1(2*ndp5),4*i5-20,ijtmp1(2*ndp5))
        call viorv  (lngth,ijtmp1(ndp5),ijtmp1(2*ndp5),ijtmp1(ndp5))
        call vands  (lngth,ijtmp1(1),msk(i6),ijtmp1(2*ndp5))
        call vshfcs (lngth,ijtmp1(2*ndp5),4*i6-24,ijtmp1(2*ndp5))
        call viorv  (lngth,ijtmp1(ndp5),ijtmp1(2*ndp5),ijtmp1(ndp5))
        call dcopy  (lngth,gen(istrtf2+ije2(i-1)+1),1,tmp,1)
        call vsort  (ijtmp1(ndp5),lngth,tmp)
        call dcopy  (lngth,tmp,1,gen(istrtt2+ije1(i-1)+1),1)
10    continue
      call prtime ("reord12  ",0,1)
      return
      end


      subroutine reord13 (namef,namet,iord)         
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)

#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2

#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)

      dimension ia(12),ja(12),msk(6),indx(ndij),dummy(ndij),iord(12)
      data msk / 'f00000'x,'f0000'x,'f000'x,'f00'x,'f0'x,'f'x /
      i1=iord(1)
      i2=iord(2)
      i3=iord(3)
      i4=iord(4)
      i5=iord(5)
      i6=iord(6)
      call tmplt12 (namet,ntapet)
      call files_ccsd (namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'reorder input tape does not exist, ',name2
        stop                    
      end if
      if (nij2.ne.nij1) then
        write (6,*) 'mismatch in numbers of blocks - reord13'
        stop                        
      end if
      istrtt1=nstrt(1,ntapet-1)+2*nij1+4
      istrtt2=nstrt(2,ntapet-1)
      istrtf1=nstrt(1,ntapef-1)+2*nij2+4
      istrtf2=nstrt(2,ntapef-1)
      call vands  (nij2,ij2,msk(i1),itmp)
      call vshfcs (nij2,itmp,4*i1-4,itmp)
      call vands  (nij2,ij2,msk(i2),itmp(ndp5))
      call vshfcs (nij2,itmp(ndp5),4*i2-8,itmp(ndp5))
      call viorv  (nij2,itmp,itmp(ndp5),itmp)
      call vands  (nij2,ij2,msk(i3),itmp(ndp5))
      call vshfcs (nij2,itmp(ndp5),4*i3-12,itmp(ndp5))
      call viorv  (nij2,itmp,itmp(ndp5),itmp)
      call vands  (nij2,ij2,msk(i4),itmp(ndp5))
      call vshfcs (nij2,itmp(ndp5),4*i4-16,itmp(ndp5))
      call viorv  (nij2,itmp,itmp(ndp5),itmp)
      call vands  (nij2,ij2,msk(i5),itmp(ndp5))
      call vshfcs (nij2,itmp(ndp5),4*i5-20,itmp(ndp5))
      call viorv  (nij2,itmp,itmp(ndp5),itmp)
      call vands  (nij2,ij2,msk(i6),itmp(ndp5))
      call vshfcs (nij2,itmp(ndp5),4*i6-24,itmp(ndp5))
      call viorv  (nij2,itmp,itmp(ndp5),itmp)
      do 30 i=1,nij2        
        dummy(i)=i*1.d0
30    continue
      call vsort (itmp,nij2,dummy)
      do 40 i=1,nij2
        d=dummy(i)-int(dummy(i))
        if (d.gt.0.5d0) then
          id=int(dummy(i))+1
        else   
          id=int(dummy(i))
        end if
        indx(i)=id
40    continue
      k=1
      do 10 i=1,nij2
        n1=ije1(i)-ije1(i-1)
        n2=ije2(indx(i))-ije2(indx(i)-1)
        if (n1.ne.n2) then
          write (6,*) 'mismatch in block lengths - reord13'
          stop                      
        end if
        do 20 j=ije2(indx(i)-1)+1,ije2(indx(i))
          tmp(k)=gen(istrtf2+j)
          k=k+1
20      continue
10    continue
      call dcopy (k,tmp,1,gen(istrtt2+1),1)   
      call prtime ("reord13 ",0,2)
      return
      end


      subroutine seltyp(iii,iiv,name,iparm)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)

c  ints collected in gen, indices in itmp. at the end, arrays are
c  sorted and expanded according to indices (in igen) generated
c  by tmplt.
      common /cij1/ nij,ij(ndij),ije(0:ndij),iiii,iivi

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)

#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)

#include "rsp2.h"
!     common /rsp2/  nrsp2,irsp2(500)

#include "cint.h"
!     common /cint/ nint,idum1,vint(ndint),iint(ndint)

      logical lrsp
      data zero / 1.d-13/
      ipack(i1,i2,i3,i4)=ior(ior(ishft(i1,24),ishft(i2,16)),
     1                       ior(ishft(i3,8) ,i4)           )
      call tmplt(iii,iiv,name,num)
      lngth=ije(nij)
      istrt1=nstrt(1,num-1)+2*nij+4
      istrt2=nstrt(2,num-1)
      call parse(iii,i1,i2,i3,i4)
      call parse(iiv,iv1,iv2,iv3,iv4)
      lrsp=(iopt(9).ne.0).and.(iv1.eq.1).and.(iv2.eq.1)
      nh=0
      nv=0
      if (i1.eq.2) nh=nh+1
      if (iv1.eq.1) nv=nv+1
      if (i2.eq.2) nh=nh+1
      if (iv2.eq.1) nv=nv+1
      if (i3.eq.2) nh=nh+1
      if (iv3.eq.1) nv=nv+1
      if (i4.eq.2) nh=nh+1
      if (iv4.eq.1) nv=nv+1
      i1234=ipack(i1,i2,i3,i4)
      iv1234=ipack(iv1,iv2,iv3,iv4)
      nout=0
      do 190 i=1,nint
      if (abs(vint(i)).lt.zero) go to 190
      j1=iand(ishft(iint(i),-24),'ff'X)
      j2=iand(ishft(iint(i),-16),'ff'X)
      j3=iand(ishft(iint(i),-8) ,'ff'X)
      j4=iand(          iint(i) ,'ff'X)
      njh=0
      njv=0
      if (ic(j1).eq.2) njh=njh+1
      if (icv(j1).eq.1) njv=njv+1
      if (ic(j2).eq.2) njh=njh+1
      if (icv(j2).eq.1) njv=njv+1
      if (ic(j3).eq.2) njh=njh+1
      if (icv(j3).eq.1) njv=njv+1
      if (ic(j4).eq.2) njh=njh+1
      if (icv(j4).eq.1) njv=njv+1
      if (njh.ne.nh) go to 190
      if (njv.lt.nv) go to 190
      if (ipack(ic(j1),ic(j2),ic(j3),ic(j4)).ne.i1234) go to 110
      if (iand(ipack(icv(j1),icv(j2),icv(j3),icv(j4)),iv1234).ne.iv1234)
     1 go to 110
      if (lrsp) then
	j12=ior(ishft(j1,8),j2)
        do 105 ii=1,nrsp2
105     if (j12.eq.irsp2(ii)) go to 110
      endif
      nout=nout+1
      gen(istrt2+nout)=vint(i)
      itmp(nout)=ipack(j1,j2,j3,j4)
110   if (j1.eq.j3) go to 120
      if (ipack(ic(j3),ic(j2),ic(j1),ic(j4)).ne.i1234) go to 120
      if (iand(ipack(icv(j3),icv(j2),icv(j1),icv(j4)),iv1234).ne.iv1234)
     1 go to 120
      if (lrsp) then
	j12=ior(ishft(j3,8),j2)
        do 115 ii=1,nrsp2
115     if (j12.eq.irsp2(ii)) go to 120
      endif
      nout=nout+1
      gen(istrt2+nout)=vint(i)
      itmp(nout)=ipack(j3,j2,j1,j4)
120   if (j2.eq.j4) go to 130
      if (ipack(ic(j1),ic(j4),ic(j3),ic(j2)).ne.i1234) go to 130
      if (iand(ipack(icv(j1),icv(j4),icv(j3),icv(j2)),iv1234).ne.iv1234)
     1 go to 130
      if (lrsp) then
	j12=ior(ishft(j1,8),j4)
        do 125 ii=1,nrsp2
125     if (j12.eq.irsp2(ii)) go to 130
      endif
      nout=nout+1
      gen(istrt2+nout)=vint(i)
      itmp(nout)=ipack(j1,j4,j3,j2)
130   if ((j1.eq.j3).or.(j2.eq.j4)) go to 140
      if (ipack(ic(j3),ic(j4),ic(j1),ic(j2)).ne.i1234) go to 140
      if (iand(ipack(icv(j3),icv(j4),icv(j1),icv(j2)),iv1234).ne.iv1234)
     1 go to 140
      if (lrsp) then
	j12=ior(ishft(j3,8),j4)
        do 135 ii=1,nrsp2
135     if (j12.eq.irsp2(ii)) go to 140
      endif
      nout=nout+1
      gen(istrt2+nout)=vint(i)
      itmp(nout)=ipack(j3,j4,j1,j2)
140   if ((j1.eq.j2).and.(j3.eq.j4)) go to 190
      if (ipack(ic(j2),ic(j1),ic(j4),ic(j3)).ne.i1234) go to 150
      if (iand(ipack(icv(j2),icv(j1),icv(j4),icv(j3)),iv1234).ne.iv1234)
     1 go to 150
      if (lrsp) then
	j12=ior(ishft(j2,8),j1)
        do 145 ii=1,nrsp2
145     if (j12.eq.irsp2(ii)) go to 150
      endif
      nout=nout+1
      gen(istrt2+nout)=vint(i)
      itmp(nout)=ipack(j2,j1,j4,j3)
150   if (j1.eq.j3) go to 160
      if (ipack(ic(j2),ic(j3),ic(j4),ic(j1)).ne.i1234) go to 160
      if (iand(ipack(icv(j2),icv(j3),icv(j4),icv(j1)),iv1234).ne.iv1234)
     1 go to 160
      if (lrsp) then
	j12=ior(ishft(j2,8),j3)
        do 155 ii=1,nrsp2
155     if (j12.eq.irsp2(ii)) go to 160
      endif
      nout=nout+1
      gen(istrt2+nout)=vint(i)
      itmp(nout)=ipack(j2,j3,j4,j1)
160   if (j2.eq.j4) go to 170
      if (ipack(ic(j4),ic(j1),ic(j2),ic(j3)).ne.i1234) go to 170
      if (iand(ipack(icv(j4),icv(j1),icv(j2),icv(j3)),iv1234).ne.iv1234)
     1 go to 170
      if (lrsp) then
	j12=ior(ishft(j4,8),j1)
        do 165 ii=1,nrsp2
165     if (j12.eq.irsp2(ii)) go to 170
      endif
      nout=nout+1
      gen(istrt2+nout)=vint(i)
      itmp(nout)=ipack(j4,j1,j2,j3)
170   if ((j1.eq.j3).or.(j2.eq.j4)) go to 190
      if (ipack(ic(j4),ic(j3),ic(j2),ic(j1)).ne.i1234) go to 190
      if (iand(ipack(icv(j4),icv(j3),icv(j2),icv(j1)),iv1234).ne.iv1234)
     1 go to 190
      if (lrsp) then
	j12=ior(ishft(j4,8),j3)
        do 175 ii=1,nrsp2
175     if (j12.eq.irsp2(ii)) go to 190
      endif
      nout=nout+1
      gen(istrt2+nout)=vint(i)
      itmp(nout)=ipack(j4,j3,j2,j1)
190   continue
      call vxsort(itmp,nout,ije(nij),gen(istrt2+1),
     1 igen(istrt1+1))
      write (6,1900) ije(nij),iii,iiv
1900  format('collected ',i9,' integrals of type ',2i6)
      call prt(name,name,iparm)
      call prtime(name,ije(nij),0)
      return
      end

      subroutine selt1(iii,iiv,name,iparm)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)

      common /cij1/ nij,ij(ndij),ije(0:ndij),iiii,iivi

#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)

#include "fock.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)

#include "rsp1.h"
!     common /rsp1/  nrsp1,irsp1(100)

      logical lrsp
      data zero / 1.d-13/
      ipack(i1,i2)=ior(ishft(i1,8),i2)
      if (iopt(8).eq.0) then
        call clear(name,iii,iiv)
        return
      endif
      call tmplt(iii,iiv,name,num)
      lngth=ije(nij)
      istrt1=nstrt(1,num-1)+2*nij+4
      istrt2=nstrt(2,num-1)
      call parse(iii,id,id,i1,i2)
      call parse(iiv,id,id,iv1,iv2)
      lrsp=(iopt(9).ne.0).and.(iv1.eq.1)
      i12=ipack(i1,i2)
      iv12=ipack(iv1,iv2)
      nout=0
      do 190 i=1,nfock
      if (abs(fock(i)).lt.zero) go to 190
      j1=iand(ishft(ifock(i),-8),'ff'X)
      j2=iand(ifock(i),'ff'X)
      if(ipack(ic(j1),ic(j2)).ne.i12) go to 110
      if(iand(ipack(icv(j1),icv(j2)),iv12).ne.iv12) go to 110
      if (lrsp) then
	j12=ior(ishft(j1,8),j2)
        do 105 ii=1,nrsp1
105     if (j12.eq.irsp1(ii)) go to 110
      endif
      nout=nout+1
      gen(istrt2+nout)=fock(i)
      itmp(nout)=ipack(j1,j2)
110   if (j1.eq.j2) go to 190
      if(ipack(ic(j2),ic(j1)).ne.i12) go to 190
      if(iand(ipack(icv(j2),icv(j1)),iv12).ne.iv12) go to 190
      if (lrsp) then
	j12=ior(ishft(j2,8),j1)
        do 115 ii=1,nrsp1
115     if (j12.eq.irsp1(ii)) go to 190
      endif
      nout=nout+1
      gen(istrt2+nout)=fock(i)
      itmp(nout)=ipack(j2,j1)
190   continue
      call vxsort(itmp,nout,ije(nij),gen(istrt2+1),
     1 igen(istrt1+1))
      write (6,1900) ije(nij),iii,iiv
1900  format('collected ',i9,' integrals of type ',2i6)
      call prt(name,name,iparm)
      call prtime(name,ije(nij),0)
      return
      end

      subroutine seltyp4(iii,iiv,name,iparm)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

c  this program for pppp only. integrals <ij|kl> collected, i.ge.j .

#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)

c  ints collected in gen, indices in itmp. at the end, arrays are
c  sorted and expanded according to indices (in igen) generated
c  by tmplt.
      common /cij1/ nij,ij(ndij),ije(0:ndij),iiii,iivi

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)

#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)

#include "cint.h"
!      common /cint/ nint,idum1,vint(ndint),iint(ndint)

      data zero / 1.e-13/

      ipack(i1,i2,i3,i4)=ior(ior(ishft(i1,24),ishft(i2,16)),
     1                       ior(ishft(i3,8),i4)           )
      if ((iii.ne.3333).or.(iiv.ne.0)) then
        write (*,*) ' seltyp4 called with incorrect iii, iiv:',iii,iiv
        call abext('seltyp4 ')
      endif
      ieof=0
      call rdfile(9,0,"pppp",ieof)
      if (ieof.eq.0) return
      call tmplt(-iii,iiv,name,num)
      lngth=ije(nij)
      istrt1=nstrt(1,num-1)+2*nij+4
      istrt2=nstrt(2,num-1)
      nout=0
      do 190 i=1,nint
      if (abs(vint(i)).lt.zero) go to 190
c     call unpack(iint(i),j1,j2,j3,j4)
      j1=iand(ishft(iint(i),-24),'ff'X)
      j2=iand(ishft(iint(i),-16),'ff'X)
      j3=iand(ishft(iint(i),-8) ,'ff'X)
      j4=iand(           iint(i),'ff'X)
      if (ic(j1).ne.3) go to 190
      if (ic(j2).ne.3) go to 190
      if (ic(j3).ne.3) go to 190
      if (ic(j4).ne.3) go to 190
      if (j1.lt.j2) go to 110
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j1.eq.j2) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack(j1,j2,j3,j4)
110   if (j1.eq.j3) go to 120
      if (j3.lt.j2) go to 120
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j2.eq.j3) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack(j3,j2,j1,j4)
120   if (j2.eq.j4) go to 130
      if (j1.lt.j4) go to 130
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j1.eq.j4) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack(j1,j4,j3,j2)
130   if ((j1.eq.j3).or.(j2.eq.j4)) go to 140
      if (j3.lt.j4) go to 140
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j3.eq.j4) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack(j3,j4,j1,j2)
140   if ((j1.eq.j2).and.(j3.eq.j4)) go to 190
      if (j2.lt.j1) go to 150
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j1.eq.j2) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack(j2,j1,j4,j3)
150   if (j1.eq.j3) go to 160
      if (j2.lt.j3) go to 160
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j3.eq.j2) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack(j2,j3,j4,j1)
160   if (j2.eq.j4) go to 170
      if (j4.lt.j1) go to 170
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j1.eq.j4) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack(j4,j1,j2,j3)
170   if ((j1.eq.j3).or.(j2.eq.j4)) go to 190
      if (j4.lt.j3) go to 190
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j3.eq.j4) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack(j4,j3,j2,j1)
190   continue
      call vxsort(itmp,nout,lngth,gen(istrt2+1),
     1 igen(istrt1+1))
c  indices not retained for pppp.
      nstrt(1,nflcur)=nstrt(1,nflcur)-lngth
      write (6,1900) lngth,iii,iiv
1900  format('collected ',i9,' integrals of type ',2i6)
      if (iopt(10).eq.0) call wrtfile(9,0,"pppp")
      call prt(name,name,iparm)
      call prtime(name,lngth,0)
      return
      end

      subroutine sqr3 (name)       
      implicit real*8 (a-h,o-z)
#include "prm.h"
!     common /prm/ iperm(0:6,720)

      call perms (3)
      call move    (name,"s1  ")
      do 10 i=2,6  
        iord1=123000+100*(iperm(1,i)+3)+10*(iperm(2,i)+3)+iperm(3,i)+3
        iord2=100000*iperm(1,i)+10000*iperm(2,i)+1000*iperm(3,i)+456
        call reord6  ("s1  ","s2  ",iord1)
        call reord7  ("s2  ","s3  ",iord2)
        call add     ("s3  ",1.d0,name,1.d0,name)
10    continue
      call prtime ("sqr3",0,2)
      return
      end

      subroutine sqr4 (name)         
      implicit real*8 (a-h,o-z)
#include "prm.h"
!     common /prm/ iperm(0:6,720)

      call perms (4)
      call move  (name,"s1  ")
      do 10 i=2,24
        iord=10000000*iperm(1,i)+1000000*iperm(2,i)+
     1  100000*iperm(3,i)+10000*iperm(4,i)+
     2  1000*(iperm(1,i)+4)+100*(iperm(2,i)+4)+
     3  10*(iperm(3,i)+4)+iperm(4,i)+4
        call reordx  ("s1  ","s2  ",iord)
        call add     ("s2  ",1.d0,name,1.d0,name)
10    continue
      call prtime ("sqr4",0,2)
      return
      end

      subroutine sqr4l3 (name)
      implicit real*8 (a-h,o-z)
#include "prm.h"
!      common /prm/ iperm(0:6,720)
      call perms (3)
      call move  (name,"s1  ")
      do 10 i=2,6  
        iord=10000000+1000000*(iperm(1,i)+1)+
     1  100000*(iperm(2,i)+1)+10000*(iperm(3,i)+1)+5000+
     2  100*(iperm(1,i)+5)+10*(iperm(2,i)+5)+iperm(3,i)+5
        call reordx  ("s1  ","s2  ",iord)
        call add     ("s2  ",1.d0,name,1.d0,name)
10    continue
      call prtime ("sqr4l3  ",0,2)
      return
      end
      subroutine sqr(namef,names)
      implicit real*8 (a-h,o-z)
      call reorder(namef,names,2143,5)
      call add(namef,1.d0,names,1.d0,namef)
      return
      end

      subroutine tmplt(ii1234,j1234,name,num)
      implicit real*8 (a-h,o-z)
!     integer*4 ch1,ch2
      integer ch1,ch2
c   ii1234<0 - collect 1/2 file (i>j)
      
#include "param_ccsd.h"

      common /cij1/ nij,ij(ndij),ije(0:ndij),iii,iiv

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)

#include "dat.h"
!      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)

#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)

      ipack(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1               ior(ishft(k,8),l))
      i1234=iabs(ii1234)
      call parse(i1234,i1,i2,i3,i4)
      call parse(j1234,j1,j2,j3,j4)
      call files_ccsd(name,num,new,1,i1234,j1234)
      if (new.le.0) return
      nij=0
      ije(0)=0
      nt=nstrt(1,num-1)+2*ndij+4
      nt0=nt
      if ((i1.eq.0).and.(i2.eq.0)) go to 310
      do 290 ica=1,nclass
      ja1=ibias(ica)+1
      ja2=ibias(ica+1)
      do 290 ja=ja1,ja2
      if ((ic(ja).ne.i1).or.(icv(ja).lt.j1)) go to 290
      do 280 icb=1,nclass
      jb1=ibias(icb)+1
      jb2=ibias(icb+1)
      do 280 jb=jb1,jb2
      if ((ic(jb).ne.i2).or.(icv(jb).lt.j2)) go to 280
      if ((ii1234.lt.0).and.(ja.lt.jb)) go to 280
      nt00=nt
      do 240 icc=1,nclass
      do 240 icd=1,nclass
c  symmetry operations
      do 170 i=1,nsym
      if (mod(isym(ica,i)+isym(icb,i)+isym(icc,i)+isym(icd,i),2).
     1 ne.0) go to 240
170   continue
      jc1=ibias(icc)+1
      jc2=ibias(icc+1)
      jd1=ibias(icd)+1
      jd2=ibias(icd+1)
      do 230 jc=jc1,jc2
      if ((ic(jc).ne.i3).or.(icv(jc).lt.j3)) go to 230
      do 220 jd=jd1,jd2
      if ((ic(jd).ne.i4).or.(icv(jd).lt.j4)) go to 220
      nt=nt+1
      igen(nt)=ipack(ja,jb,jc,jd)
220   continue
230   continue
240   continue
      if (nt.gt.nt00) then
        ije1=ije(nij)+1
        nij=nij+1
        ije(nij)=nt-nt0
	ij(nij)=ior(ishft(ja,8),jb)
        call ivsort(igen(nt00+1),nt-nt00)
C	print '(Z8)',(igen(i),i=nt00+1,nt)
        if (nij.gt.ndij) write (*,*) ' nij reached ',nij,' at jab=',ja,jb
        call dimx(nij,ndij,"tmpndij")
      endif
280   continue
290   continue
c
      iii=i1234
      iiv=j1234
      call dimx(nij,ndij,"tmpndij")
      go to 410
c  one electron
310   do 360 icc=1,nclass
      do 360 icd=1,nclass
      do 320 i=1,nsym
      if (mod(isym(icc,i)+isym(icd,i),2).ne.0) go to 360
320   continue
      jc1=ibias(icc)+1
      jc2=ibias(icc+1)
      jd1=ibias(icd)+1
      jd2=ibias(icd+1)
      do 350 jc=jc1,jc2
      if ((ic(jc).ne.i3).or.(icv(jc).lt.j3)) go to 350
      do 340 jd=jd1,jd2
      if ((ic(jd).ne.i4).or.(icv(jd).lt.j4)) go to 340
      nt=nt+1
      igen(nt)=ipack(0,0,jc,jd)
340   continue
350   continue
360   continue
      nij=1
      ije(1)=nt-nt0
      ij(1)=0
      call ivsort(igen(nt0+1),nt-nt0)
      iii=i1234
      iiv=j1234
410   call kcopy(nij+1,nij,igen(nstrt(1,num-1)+1))
      call kcopy(nij+1,ije(0),igen(nstrt(1,num-1)+nij+2))
      igen(nstrt(1,num-1)+2*nij+3)=iii
      igen(nstrt(1,num-1)+2*nij+4)=iiv
      do 420 i=1,nt-nt0
420   igen(nt0-2*(ndij-nij)+i)=igen(nt0+i)
      nstrt(1,num)=nstrt(1,num-1)+(nt-nt0)+2*nij+4
      nstrt(2,num)=nstrt(2,num-1)+(nt-nt0)
      call dimx(nstrt(1,num),nbig,"tmplt1")
      call dimx(nstrt(2,num),nbig,"tmplt2")
c   special for scratch files, which start with "s".
           itmp=4H       
           jtmp='ffffff00'X
           ch1 = iand(itmp,jtmp)
           ch2 = iand (name, '000000ff'X)
           isaa = ior (ch1, ch2)
      i_4hs=4hs
!     if (isaa.eq.4hs   ) then
      if (isaa.eq.i_4hs) then
        call dimx(nt-nt0,ndp3,"tmplt s")
        nstrt(1,num)=nstrt(1,num-1)+ndp3+2*ndij+4
        nstrt(2,num)=nstrt(2,num-1)+ndp3
        call dimx(nstrt(1,num),nbig,"tmplt1")
        call dimx(nstrt(2,num),nbig,"tmplt2")
      endif
      return
      end

      subroutine tmplt6 (iii,iiv,name,num)
      implicit real*8 (a-h,o-z)
!     integer*4 ch1,ch2
      integer ch1,ch2
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)

      common /cij1/ nij,ij(ndij),ije(0:ndij),iii1,iiv1

#include "dat.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)

#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)

      ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)
      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
      call parse6 (iii,i1,i2,i3,i4,i5,i6)
      call parse6 (iiv,j1,j2,j3,j4,j5,j6)
      call files_ccsd (name,num,new,1,iii,iiv)
      if (new.le.0) return
      nij=0
      ije(0)=0
      nt=nstrt(1,num-1)+2*ndij+4
      nt0=nt
      do 10 ica=1,nclass
      ja1=ibias(ica)+1
      ja2=ibias(ica+1)
      do 10 ja=ja1,ja2
       if ((ic(ja).eq.i1).and.(icv(ja).ge.j1)) then
        do 20 icb=1,nclass
        jb1=ibias(icb)+1
        jb2=ibias(icb+1)
        do 20 jb=jb1,jb2
         if ((ic(jb).eq.i2).and.(icv(jb).ge.j2)) then
          do 30 icc=1,nclass
          jc1=ibias(icc)+1
          jc2=ibias(icc+1)
          do 30 jc=jc1,jc2
           if ((ic(jc).eq.i3).and.(icv(jc).ge.j3)) then
            nt00=nt
            do 40 icd=1,nclass
            do 40 ice=1,nclass
            do 40 icf=1,nclass
             iflag=1
             do 50 ii=1,nsym
              if (mod(isym(ica,ii)+isym(icb,ii)+
     1        isym(icc,ii)+isym(icd,ii)+
     2        isym(ice,ii)+isym(icf,ii),2).ne.0) then
               iflag=0
              end if
50           continue
             if (iflag.eq.1) then
              jd1=ibias(icd)+1
              jd2=ibias(icd+1)
              je1=ibias(ice)+1
              je2=ibias(ice+1)
              jf1=ibias(icf)+1
              jf2=ibias(icf+1)
              do 60 jd=jd1,jd2
               if ((ic(jd).eq.i4).and.(icv(jd).ge.j4)) then
                do 70 je=je1,je2
                 if ((ic(je).eq.i5).and.(icv(je).ge.j5)) then
                  do 80 jf=jf1,jf2
                   if ((ic(jf).eq.i6).and.(icv(jf).ge.j6)) then
                    nt=nt+1
                    igen(nt)=ipack4(jc,jd,je,jf)
                   end if
80                continue
                 end if
70              continue
               end if
60            continue
             end if
40          continue
            if (nt.gt.nt00) then
             call dimx (nt,nbig,"tmp6    ")
             nij=nij+1
             ije(nij)=nt-nt0
             ij(nij)=ipack3(ja,jb,jc)
             call ivsort (igen(nt00+1),nt-nt00)
             if (nij.gt.ndij) write (6,*) 'nij reached ',
     1       nij,' at jabc=',ja,jb,jc
             call dimx (nij,ndij,"tmpndij ")
            end if
           end if
30        continue
         end if
20      continue
       end if
10    continue
      call dimx (ije(nij),ndp3,"tmp6    ")
      igen(nstrt(1,num-1)+1)=nij
      call kcopy (nij,ij(1),igen(nstrt(1,num-1)+2))
      call kcopy (nij+1,ije(0),igen(nstrt(1,num-1)+nij+2))
      igen(nstrt(1,num-1)+2*nij+3)=iii
      igen(nstrt(1,num-1)+2*nij+4)=iiv
      do 90 i=1,nt-nt0
        igen(nt0-2*(ndij-nij)+i)=igen(nt0+i)
90    continue
      nstrt(1,num)=nstrt(1,num-1)+2*nij+4+(nt-nt0)
      nstrt(2,num)=nstrt(2,num-1)+(nt-nt0)
      call dimx (nstrt(1,num),nbig,"tmplt61 ")
      call dimx (nstrt(2,num),nbig,"tmplt62 ")
c   special for scratch files, which start with s
           itmp=4H       
           jtmp='ffffff00'X
           ch1 = iand(itmp,jtmp)
           ch2 = iand (name, '000000ff'X)
           isaa = ior (ch1, ch2)

       i_4hs=4hs
!     if (isaa.eq.4hs   ) then
      if (isaa.eq.i_4hs) then
        call dimx (nt-nt0,ndp3,"tmplt s ")
        nstrt(1,num)=nstrt(1,num-1)+ndp3+2*ndij+4
        nstrt(2,num)=nstrt(2,num-1)+ndp3
        call dimx (nstrt(1,num),nbig,"tmplt1   ")
        call dimx (nstrt(2,num),nbig,"tmplt2   ")
      end if
      return
      end

      subroutine tmplt8 (iii,iiv,name,num)
      implicit real*8 (a-h,o-z)
!     integer*4 ch1,ch2
      integer ch1,ch2
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)

      common /cij1/ nij,ij(ndij),ije(0:ndij),iii1,iiv1

#include "dat.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)

#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)

#include "dpd.h"
!      common /dpd/ nir,igmt(8,8),ici(200)

      ipack(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
      call parse8 (iii,i1,i2,i3,i4,i5,i6,i7,i8)
      call parse8 (iiv,j1,j2,j3,j4,j5,j6,j7,j8)
      call files_ccsd (name,num,new,1,iii,iiv)
      if (new.le.0) return
      nij=0
      ije(0)=0
      nt=nstrt(1,num-1)+2*ndij+4
      nt0=nt
      do 10 ica=1,nclass
      ja1=ibias(ica)+1
      ja2=ibias(ica+1)
      do 10 ja=ja1,ja2
       if ((ic(ja).eq.i1).and.(icv(ja).ge.j1)) then
        do 20 icb=1,nclass
        jb1=ibias(icb)+1
        jb2=ibias(icb+1)
        do 20 jb=jb1,jb2
         if ((ic(jb).eq.i2).and.(icv(jb).ge.j2)) then
          do 30 icc=1,nclass
          jc1=ibias(icc)+1
          jc2=ibias(icc+1)
          do 30 jc=jc1,jc2
           if ((ic(jc).eq.i3).and.(icv(jc).ge.j3)) then
            do 40 icd=1,nclass
            jd1=ibias(icd)+1
            jd2=ibias(icd+1)
            do 40 jd=jd1,jd2
             if ((ic(jd).eq.i4).and.(icv(jd).ge.j4)) then
              nt00=nt
              irrout=igmt(ica,icb)        
              irrout=igmt(irrout,icc)    
              irrout=igmt(irrout,icd)     
              do 50 ice=1,nclass
              do 50 icf=1,nclass
              do 50 icg=1,nclass
              do 50 ich=1,nclass
               irrin=igmt(irrout,ice)
               irrin=igmt(irrin,icf)
               irrin=igmt(irrin,icg)
               irrin=igmt(irrin,ich)
               if (irrin.eq.1) then
                je1=ibias(ice)+1
                je2=ibias(ice+1)
                jf1=ibias(icf)+1
                jf2=ibias(icf+1)
                jg1=ibias(icg)+1
                jg2=ibias(icg+1)
                jh1=ibias(ich)+1
                jh2=ibias(ich+1)
                do 70 je=je1,je2
                 if ((ic(je).eq.i5).and.
     1           (icv(je).ge.j5)) then
                  do 80 jf=jf1,jf2
                   if ((ic(jf).eq.i6).and.
     1             (icv(jf).ge.j6)) then
                    do 90 jg=jg1,jg2
                     if ((ic(jg).eq.i7).and.
     1               (icv(jg).ge.j7)) then
                      do 100 jh=jh1,jh2
                       if ((ic(jh).eq.i8).and.
     1                 (icv(jh).ge.j8)) then
                        nt=nt+1
                        call dimx (nt,nbig,"tmp8")
                        igen(nt)=ipack(je,jf,jg,jh)
                       end if
100                   continue
                     end if
90                  continue
                   end if
80                continue
                 end if
70              continue
               end if
50            continue
              if (nt.gt.nt00) then
               nij=nij+1
               ije(nij)=nt-nt0
               ij(nij)=ipack(ja,jb,jc,jd)
               call ivsort (igen(nt00+1),nt-nt00)
               if (nij.gt.ndij) write (6,*) 'nij reached ',
     1         nij,' at jabcd=',ja,jb,jc,jd
               call dimx (nij,ndij,"tmpndij ")
              end if
             end if
40          continue
           end if
30        continue
         end if
20      continue
       end if
10    continue
      igen(nstrt(1,num-1)+1)=nij
      call kcopy (nij,ij(1),igen(nstrt(1,num-1)+2))
      call kcopy (nij+1,ije(0),igen(nstrt(1,num-1)+nij+2))
      igen(nstrt(1,num-1)+2*nij+3)=iii
      igen(nstrt(1,num-1)+2*nij+4)=iiv
      do 110 i=1,nt-nt0
        igen(nt0-2*(ndij-nij)+i)=igen(nt0+i)
110   continue
      nstrt(1,num)=nstrt(1,num-1)+2*nij+4+(nt-nt0)
      nstrt(2,num)=nstrt(2,num-1)+(nt-nt0)
      call dimx (nstrt(1,num),nbig,"tmplt61 ")
      call dimx (nstrt(2,num),nbig,"tmplt62 ")
c   special for scratch files, which start with s
           itmp=4H      
           jtmp='ffffff00'X
           ch1 = iand(itmp,jtmp)
           ch2 = iand (name, '000000ff'X)
           isaa = ior (ch1, ch2)

      i_4hs=4hs
      if (isaa.eq.i_4hs) then
        call dimx (nt-nt0,ndp3,"tmplt s ")
        nstrt(1,num)=nstrt(1,num-1)+ndp3+2*ndij+4
        nstrt(2,num)=nstrt(2,num-1)+ndp3
        call dimx (nstrt(1,num),nbig,"tmplt1   ")
        call dimx (nstrt(2,num),nbig,"tmplt2   ")
      end if
      return
      end

      subroutine tmplt10 (name,num)
      implicit real*8 (a-h,o-z)
!     integer*4 ch1,ch2
      integer ch1,ch2
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)

      common /cij1/ nij,ij(ndij),ije(0:ndij),iii1,iiv1

#include "dat.h"
!      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)

#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)

#include "idv.h"
!     common /idv/ iv(200),nvalo

      dimension icsym(200)
      ipack(i,j,k,l,m)=ior(ishft(i,16),
     1ior(ior(ishft(j,12),ishft(k,8)),ior(ishft(l,4),m)))
      iii=5
      iiv=1
      call files_ccsd (name,num,new,1,iii,iiv)
      if (new.le.0) return
      nij=0
      ije(0)=0
      nt=nstrt(1,num-1)+2*ndij+4
      nt0=nt
      isum=0
      do 10 i=1,nclass
        do 20 j=1,iclass(i)
          icsym(isum+j)=i
20      continue
        isum=isum+iclass(i)
10    continue
      do 30 ja=1,nvalo
      do 30 jb=1,nvalo
      do 30 jc=1,nvalo
      do 30 jd=1,nvalo
      do 30 je=1,nvalo
       nt00=nt
       do 90 jf=1,nvalo
       do 90 jg=1,nvalo
       do 90 jh=1,nvalo
       do 90 ji=1,nvalo
       do 90 jj=1,nvalo
        iflag=1
        do 100 ii=1,nsym
         if (mod(isym((icsym(iv(ja))),ii)+
     1   isym((icsym(iv(jb))),ii)+isym((icsym(iv(jc))),ii)+
     2   isym((icsym(iv(jd))),ii)+isym((icsym(iv(je))),ii)+
     3   isym((icsym(iv(jf))),ii)+isym((icsym(iv(jg))),ii)+
     4   isym((icsym(iv(jh))),ii)+isym((icsym(iv(ji))),ii)+
     5   isym((icsym(iv(jj))),ii),2).ne.0) then
          iflag=0
         end if
100     continue
        if (iflag.eq.1) then
         nt=nt+1
         call dimx (nt,nbig,"tmp10    ")
         igen(nt)=ipack(jf,jg,jh,ji,jj)
        end if
90     continue
       if (nt.gt.nt00) then
        nij=nij+1
        ije(nij)=nt-nt0
        ij(nij)=ipack(ja,jb,jc,jd,je)          
        call ivsort (igen(nt00+1),nt-nt00)
        if (nij.gt.ndij) write (6,*) 'nij reached ',
     1  nij,' at jabcde=',ja,jb,jc,jd,je
        call dimx (nij,ndij,"tmpndij ")
       end if
30    continue
      igen(nstrt(1,num-1)+1)=nij
      call kcopy (nij,ij(1),igen(nstrt(1,num-1)+2))
      call kcopy (nij+1,ije(0),igen(nstrt(1,num-1)+nij+2))
      igen(nstrt(1,num-1)+2*nij+3)=iii
      igen(nstrt(1,num-1)+2*nij+4)=iiv
      do 150 i=1,nt-nt0
        igen(nt0-2*(ndij-nij)+i)=igen(nt0+i)
150   continue
      nstrt(1,num)=nstrt(1,num-1)+2*nij+4+(nt-nt0)
      nstrt(2,num)=nstrt(2,num-1)+(nt-nt0)
      call dimx (nstrt(1,num),nbig,"tmplt61 ")
      call dimx (nstrt(2,num),nbig,"tmplt62 ")
c   special for scratch files, which start with s
           itmp=4H       
           jtmp='ffffff00'X
           ch1 = iand(itmp,jtmp)
           ch2 = iand (name, '000000ff'X)
           isaa = ior (ch1, ch2)

      i_4hs=4hs
!     if (isaa.eq.4hs   ) then
      if (isaa.eq.i_4hs ) then
        call dimx (nt-nt0,ndp3,"tmplt s ")
        nstrt(1,num)=nstrt(1,num-1)+ndp3+2*ndij+4
        nstrt(2,num)=nstrt(2,num-1)+ndp3
        call dimx (nstrt(1,num),nbig,"tmplt1   ")
        call dimx (nstrt(2,num),nbig,"tmplt2   ")
      end if
      return
      end

      subroutine tmplt12 (name,num)
      implicit real*8 (a-h,o-z)
!     integer*4 ch1,ch2
      integer ch1,ch2
      
#include "param_ccsd.h"

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)

#include "cij1.h" ! attention, variable changed names !
!     common /cij1/ nij,ij(ndij),ije(0:ndij),iii1,iiv1

#include "dat.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)

      common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)

      common /idv/ iv(200),nvalo

      dimension icsym(200)

      ipack(i,j,k,l,m,n)=ior(ior(ishft(i,20),ishft(j,16)),
     1ior(ior(ishft(k,12),ishft(l,8)),ior(ishft(m,4),n)))

      nij=nij1 ! mi: from common block

      iii=6
      iiv=1
      call files_ccsd (name,num,new,1,iii,iiv)
      if (new.le.0) return
      nij=0
!     ije(0)=0
      ije1(0)=0
      nt=nstrt(1,num-1)+2*ndij+4
      nt0=nt
      isum=0
      do 10 i=1,nclass
        do 20 j=1,iclass(i)
          icsym(isum+j)=i
20      continue
        isum=isum+iclass(i)
10    continue
      do 30 ja=1,nvalo
      do 30 jb=1,nvalo
      do 30 jc=1,nvalo
      do 30 jd=1,nvalo
      do 30 je=1,nvalo
      do 30 jf=1,nvalo
       nt00=nt
       do 90 jg=1,nvalo
       do 90 jh=1,nvalo
       do 90 ji=1,nvalo
       do 90 jj=1,nvalo
       do 90 jk=1,nvalo
       do 90 jl=1,nvalo
        iflag=1
        do 100 ii=1,nsym
         if (mod(isym((icsym(iv(ja))),ii)+
     1   isym((icsym(iv(jb))),ii)+isym((icsym(iv(jc))),ii)+
     2   isym((icsym(iv(jd))),ii)+isym((icsym(iv(je))),ii)+
     3   isym((icsym(iv(jf))),ii)+isym((icsym(iv(jg))),ii)+
     4   isym((icsym(iv(jh))),ii)+isym((icsym(iv(ji))),ii)+
     5   isym((icsym(iv(jj))),ii)+isym((icsym(iv(jk))),ii)+
     6   isym((icsym(iv(jl))),ii),2).ne.0) then
          iflag=0
         end if
100     continue
        if (iflag.eq.1) then
         nt=nt+1
         call dimx (nt,nbig,"tmp12   ")
         igen(nt)=ipack(jg,jh,ji,jj,jk,jl)
        end if
90     continue
       if (nt.gt.nt00) then
        nij=nij+1
!       ije(nij)=nt-nt0
        ije1(nij)=nt-nt0
!       ij(nij)=ipack(ja,jb,jc,jd,je,jf) 
        ij1(nij)=ipack(ja,jb,jc,jd,je,jf) 
        call ivsort (igen(nt00+1),nt-nt00)
        if (nij.gt.ndij) write (6,*) 'nij reached ',
     1  nij,' at jabcdef=',ja,jb,jc,jd,je,jf
        call dimx (nij,ndij,"tmpndij ")
       end if
30    continue
      igen(nstrt(1,num-1)+1)=nij
!     call kcopy (nij,ij(1),igen(nstrt(1,num-1)+2))
      call kcopy (nij,ij1(1),igen(nstrt(1,num-1)+2))
!     call kcopy (nij+1,ije(0),igen(nstrt(1,num-1)+nij+2))
      call kcopy (nij+1,ije1(0),igen(nstrt(1,num-1)+nij+2))
      igen(nstrt(1,num-1)+2*nij+3)=iii
      igen(nstrt(1,num-1)+2*nij+4)=iiv
      do 150 i=1,nt-nt0
        igen(nt0-2*(ndij-nij)+i)=igen(nt0+i)
150   continue
      nstrt(1,num)=nstrt(1,num-1)+2*nij+4+(nt-nt0)
      nstrt(2,num)=nstrt(2,num-1)+(nt-nt0)
      call dimx (nstrt(1,num),nbig,"tmplt61 ")
      call dimx (nstrt(2,num),nbig,"tmplt62 ")
c   special for scratch files, which start with s
           itmp=4H       
           jtmp='ffffff00'X
           ch1 = iand(itmp,jtmp)
           ch2 = iand (name, '000000ff'X)
           isaa = ior (ch1, ch2)

       i_4hs=4hs
!     if (isaa.eq.4hs   ) then
      if (isaa.eq.i_4hs) then
        call dimx (nt-nt0,ndp3,"tmplt s ")
        nstrt(1,num)=nstrt(1,num-1)+ndp3+2*ndij+4
        nstrt(2,num)=nstrt(2,num-1)+ndp3
        call dimx (nstrt(1,num),nbig,"tmplt1   ")
        call dimx (nstrt(2,num),nbig,"tmplt2   ")
      end if
      return
      end

      subroutine t1sqr(name1a,name1b,name2,namet)
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)

#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij3.h"
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3

      call files_ccsd(name1a,ntape1a,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t1a file in t1sqr does not exist'
        call abext('t1sqr   ')
      endif
      call files_ccsd(name1b,ntape1b,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t1b file in t1sqr does not exist'
        call abext('t1sqr   ')
      endif
      iii=(iii2/10)*1000+(iii3/10)*100+mod(iii2,10)*10+mod(iii3,10)
      iiv=(iiv2/10)*1000+(iiv3/10)*100+mod(iiv2,10)*10+mod(iiv3,10)
      call tmplt(iii,iiv,namet,ntapet)
c
      ntmp=0
      nsize1=ije1(nij1)
      nsize2=ije2(nij2)
      nsize3=ije3(nij3)
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      istrt21=nstrt(1,ntape1a-1)+2*nij2+4
      istrt22=nstrt(2,ntape1a-1)
      istrt31=nstrt(1,ntape1b-1)+2*nij3+4
      istrt32=nstrt(2,ntape1b-1)
      do 120 i2=1,nsize2
      do 120 i3=1,nsize3
      ntmp=ntmp+1
      gen(istrt12+ntmp)=gen(istrt22+i2)*gen(istrt32+i3)
      it2=igen(istrt21+i2)
      it3=igen(istrt31+i3)
120   itmp(ntmp)=ior(ior(ishft(iand(ishft(it2,-8),'ff'X),24),
     1                   ishft(iand(ishft(it3,-8),'ff'X),16)),
     2               ior(ishft(iand(it2,'ff'X)          , 8),
     3                   iand(it3,'ff'X)                    ))
      call vxsort(itmp,ntmp,nsize1,gen(istrt12+1),igen(istrt11+1))
      call add(namet,1.d0,name2,1.d0,namet)
      return
      end

      subroutine t1t2 (name1,name2,name3,names,iord,iparm) 
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /ctmp/ tmp(ndp4),itmp(ndp4)

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)

#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2

#include "cij3.h"
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3

      dimension ia(5,6),ind(12),ninds(ndij)
      data ia / 5,1,1,0,0,1,5,2,0,0,2,2,5,0,0,
     1          6,3,3,0,0,3,6,4,0,0,4,4,6,0,0 /
      call files_ccsd (name3,ntape3,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t3 file in t1t2 does not exist'
        stop
      endif
      call tmplt6 (iii1,iiv1,names,ntapes)
      call files_ccsd (names,ntapes,new,1,-1,-1)
      istrt11=nstrt(1,ntapes-1)+2*nij1+4
      istrt12=nstrt(2,ntapes-1)
      call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1) 
      call tchng6 (nij1,ij1,ije1,ntapes,2)
      do 60 i=1,nij1
        ninds(i)=0
60    continue
      call files_ccsd (name1,ntape1,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t1 file in t1t2 does not exist'
        stop
      endif
      call files_ccsd (name2,ntape2,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t2 file in t1t2 does not exist'
        stop
      endif
      istrt21=nstrt(1,ntape2-1)+2*nij2+4
      istrt22=nstrt(2,ntape2-1)
      istrt31=nstrt(1,ntape1-1)+2*nij3+4
      istrt32=nstrt(2,ntape1-1)
      call parse6 (iord,ia(4,1),ia(4,2),ia(4,3),ia(4,4),ia(4,5),ia(4,6))
      do 10 i=1,6
        ia(5,i)=ia(iparm,ia(4,i))
10    continue
      do 20 i=1,ije2(nij2)
        do 30 j=1,ije3(nij3)
          v=gen(istrt22+i)*gen(istrt32+j)
          if (abs(v).gt.1.d-13) then
            ii=igen(istrt21+i)
            do 40 k=1,4
              ind(k+6)=iand(ishft(ii,(k-4)*8),'ff'x)
40          continue
            ii=igen(istrt31+j)
            ind(11)=iand(ishft(ii,-8),'ff'x)
            ind(12)=iand(ii,'ff'x)
            do 50 k=1,6
              ind(k)=ind(ia(5,k)+6)
50          continue
            idnew=ior(ishft(ind(1),8),ind(2))
            indnew=ior(ior(ishft(ind(3),24),ishft(ind(4),16)),
     1      ior(ishft(ind(5),8),ind(6)))
            iquit=0
            k=1
            do while (iquit.eq.0)
              if (ij1(k).lt.idnew) then
                k=k+1
              else if (ij1(k).eq.idnew) then
                ninds(k)=ninds(k)+1
                itmp(ije1(k-1)+ninds(k))=indnew
                gen(istrt12+ije1(k-1)+ninds(k))=v
                iquit=1
              else
                write (6,*) 't1t2 - block id not found '
                stop
              end if
            end do
          end if
30      continue
20    continue
      do 70 i=1,nij1
        if (ninds(i).gt.0) then
          call vxsort (itmp(ije1(i-1)+1),ninds(i),ije1(i)-ije1(i-1),
     1    gen(istrt12+ije1(i-1)+1),igen(istrt11+ije1(i-1)+1))
        end if
70    continue
      call add (names,1.d0,name3,1.d0,names)
      return
      end

      subroutine t1t3 (name1,name3,name4,namet,iparm) 
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)

#include "cij1.h"
!      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2

#include "cij3.h"
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3

      dimension ia(4,8),ind(16),ninds(ndij)
      data ia / 1,3,3,3,3,1,4,4,4,4,1,5,5,5,5,1, 
     1          2,6,6,6,6,2,7,7,7,7,2,8,8,8,8,2 /
      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
      call files_ccsd (name4,ntape4,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t4 file in t1t3 does not exist'
        stop
      endif
      call tmplt8 (iii1,iiv1,namet,ntapet)
      call files_ccsd (namet,ntapet,new,1,-1,-1)
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1) 
      do 60 i=1,nij1
        ninds(i)=0
60    continue
      call files_ccsd (name1,ntape1,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t1 file in t1t3 does not exist'
        stop
      endif
      call files_ccsd (name3,ntape3,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t3 file in t1t3 does not exist'
        stop
      endif
      istrt21=nstrt(1,ntape1-1)+2*nij2+4
      istrt22=nstrt(2,ntape1-1)
      istrt31=nstrt(1,ntape3-1)+2*nij3+4
      istrt32=nstrt(2,ntape3-1)
      do 20 i=1,ije2(nij2)
        do 30 j=1,nij3
          do 40 k=ije3(j-1)+1,ije3(j)
            v=gen(istrt22+i)*gen(istrt32+k)
            if (abs(v).gt.1.d-13) then
              ii=igen(istrt21+i)
              ind(9)=iand(ishft(ii,-8),'ff'x)
              ind(10)=iand(ii,'ff'x)
              ii=ij3(j) 
              jj=igen(istrt31+k)
              do 50 m=1,3
                ind(m+10)=iand(ishft(ii,(m-3)*8),'ff'x)
                ind(m+13)=iand(ishft(jj,(m-3)*8),'ff'x)
50            continue
              do 70 m=1,8
                ind(m)=ind(ia(iparm,m)+8)
70            continue
              idnew=ipack4(ind(1),ind(2),ind(3),ind(4))
              indnew=ipack4(ind(5),ind(6),ind(7),ind(8))
              iquit=0
              m=1
              do while (iquit.eq.0)
                if (ij1(m).lt.idnew) then
                  m=m+1
                else if (ij1(m).eq.idnew) then
                  ninds(m)=ninds(m)+1
                  itmp(ije1(m-1)+ninds(m))=indnew
                  gen(istrt12+ije1(m-1)+ninds(m))=v
                  iquit=1
                else
                  write (6,*) 't1t3 - block id not found '
                  stop
                end if
              end do
            end if
40        continue
30      continue
20    continue
      do 80 i=1,nij1
        if (ninds(i).gt.0) then
          call vxsort (itmp(ije1(i-1)+1),ninds(i),ije1(i)-ije1(i-1),
     1    gen(istrt12+ije1(i-1)+1),igen(istrt11+ije1(i-1)+1))
        end if
80    continue
      call add (namet,1.d0,name4,1.d0,namet)
      return
      end

      subroutine t2sqr (name2a,name2b,name4,namet,iparm) 
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      common /ctmp/ tmp(ndp4),itmp(ndp4)

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)

#include "cij1.h"
!      common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

#include "cij2.h"
!      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2

#include "cij3.h"
!      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3

      dimension ia(3,8),ind(16),ninds(ndij)

      data ia / 1,1,1,2,5,5,5,2,6,6,6,2,3,3,3,4,7,7,7,4,8,8,8,4 /
      ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),
     1ior(ishft(k,8),l))
      iparm=iparm-1
      call files_ccsd (name4,ntape4,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t4 file in t2sqr does not exist'
        stop
      endif
      call tmplt8 (iii1,iiv1,namet,ntapet)
      call files_ccsd (namet,ntapet,new,1,-1,-1)
      istrt11=nstrt(1,ntapet-1)+2*nij1+4
      istrt12=nstrt(2,ntapet-1)
      call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1) 
      do 60 i=1,nij1
        ninds(i)=0
60    continue
      call files_ccsd (name2a,ntape2a,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t2a file in t2sqr does not exist'
        stop
      endif
      call files_ccsd (name2b,ntape2b,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t2b file in t2sqr does not exist'
        stop
      endif
      istrt21=nstrt(1,ntape2a-1)+2*nij2+4
      istrt22=nstrt(2,ntape2a-1)
      istrt31=nstrt(1,ntape2b-1)+2*nij3+4
      istrt32=nstrt(2,ntape2b-1)
      do 20 i=1,ije2(nij2)
        do 30 j=1,ije3(nij3)
          v=gen(istrt22+i)*gen(istrt32+j)
          if (abs(v).gt.1.d-13) then
            ii=igen(istrt21+i)
            jj=igen(istrt31+j)
            do 40 k=1,4
              ind(k+8)=iand(ishft(ii,(k-4)*8),'ff'x)
              ind(k+12)=iand(ishft(jj,(k-4)*8),'ff'x)
40          continue
            do 50 k=1,8
              ind(k)=ind(ia(iparm,k)+8)
50          continue
            idnew=ipack4(ind(1),ind(2),ind(3),ind(4))
            indnew=ipack4(ind(5),ind(6),ind(7),ind(8))
            iquit=0
            k=1
            do while (iquit.eq.0)
              if (ij1(k).lt.idnew) then
                k=k+1
              else if (ij1(k).eq.idnew) then
                ninds(k)=ninds(k)+1
                itmp(ije1(k-1)+ninds(k))=indnew
                gen(istrt12+ije1(k-1)+ninds(k))=v
                iquit=1
              else
                write (6,*) 't2sqr - block id not found '
                stop
              end if
            end do
          end if
30      continue
20    continue
      do 70 i=1,nij1
        if (ninds(i).gt.0) then
          call vxsort (itmp(ije1(i-1)+1),ninds(i),ije1(i)-ije1(i-1),
     1    gen(istrt12+ije1(i-1)+1),igen(istrt11+ije1(i-1)+1))
        end if
70    continue
      call add (namet,1.d0,name4,1.d0,namet)
      return
      end
      subroutine tchng6 (nij,ij,ije,ntapef,iparm)       
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      dimension ij(ndij),ije(0:ndij),ij1(ndij),ije1(0:ndij),itmp(ndij)

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1iiia(200),iiva(200)

#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)

      call dimx (ije(nij),ndp3,"tchng6  ")
      if (iparm.eq.2) then
        nij1=0
        iirp='ffff'x
        do 10 i=1,nij
          ii=ishft(ij(i),-8)
          if (ii.ne.iirp) then
            nij1=nij1+1
            iirp=ii
            ij1(nij1)=ii
            ije1(nij1-1)=ije(i-1)
          end if
10      continue
        ije1(nij1)=ije(nij)
        nij=nij1
        call kcopy (nij,ij1(1),ij(1))
        call kcopy (nij+1,ije1(0),ije(0))
      else if (iparm.eq.4) then
        nij1=0
        istrt=nstrt(1,ntapef-1)+2*nij+4
        do 30 i=1,nij
          iirp='ffffffff'x
          do 40 j=ije(i-1)+1,ije(i)
            ii=iand(ishft(igen(istrt+j),-16),'ff'x)
            if (ii.ne.iirp) then
              nij1=nij1+1
              iirp=ii
              ij1(nij1)=ior(ishft(ij(i),8),ii)
              ije1(nij1-1)=j-1
            end if
40        continue
30      continue
        call dimx (nij1,ndij,"tchng6  ")
        ije1(nij1)=ije(nij)
        nij=nij1
        call kcopy (nij,ij1(1),ij(1))
        call kcopy (nij+1,ije1(0),ije(0))
      end if
      return
      end
      subroutine tchng8 (nij,ij,ije,ntapef)       
      implicit real*8 (a-h,o-z)
      
#include "param_ccsd.h"

      dimension ij(ndij),ije(0:ndij),ije1(0:ndij)

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1iiia(200),iiva(200)

#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)

#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)

      call dimx (ije(nij),ndp3,"tchng8  ")
      nij1=0
      istrt=nstrt(1,ntapef-1)+2*nij+4
      do 10 i=1,nij
        iirp='ffffffff'x
        do 20 j=ije(i-1)+1,ije(i)
          ii=iand(ishft(igen(istrt+j),-16),'ffff'x)
          if (ii.ne.iirp) then
            nij1=nij1+1
            iirp=ii
            itmp(nij1)=iand(ishft(ij(i),-8),'ffffff'x)
            jj=iand(ij(i),'ff'x)
            itmp(ndp3+nij1)=ior(ishft(jj,16),ii)
            ije1(nij1-1)=j-1
          end if
20      continue
10    continue
      call dimx (nij1,ndij,"tchng8  ")
      ije1(nij1)=ije(nij)
      nij=nij1
      call kcopy (nij+1,ije1(0),ije(0))
      return
      end

      subroutine tmpcng(nij,ij,ije,ntapef,nf,nt)
      implicit real*8 (a-h,o-z)
!     integer*4 ch1,ch2
      integer ch1,ch2
c   change nf-index grouping to nt-index grouping.
      
#include "param_ccsd.h"

      dimension ij(nij),ije(0:nij),ij1(ndij),ije1(0:ndij)

#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)

#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)

      if (nt.eq.nf) return
      if (nt.gt.nf) go to 310
      nij1=0
      iirp='ffff'X
      do 160 i=1,nij
      i12=ishft(ij(i),-8*(nf-nt))
      if (i12.ne.iirp) then
        nij1=nij1+1
        ij1(nij1)=i12
        iirp=i12
        ije1(nij1-1)=ije(i-1)
      endif
160   continue
      ije1(nij1)=ije(nij)
      nij=nij1
      call kcopy(nij,ij1,ij)
      call kcopy(nij+1,ije1(0),ije(0))
      return
310   istrt=nstrt(1,ntapef-1)+2*nij+4
      nsize=igen(istrt-2)
      iirp='ffff'X
      nij=0
      do 320 i=1,nsize
      i12=ishft(igen(istrt+i),8*nt-32)
      if (i12.ne.iirp) then
        nij=nij+1
        ij(nij)=i12
        iirp=i12
        ije(nij-1)=i-1
      endif
320   continue
      ije(nij)=nsize
      call dimx(nij,ndij,"tmpcng")
      return
      end

      function v1det (ja,jb)
      implicit real*8 (a-h,o-z)

#include "dat.h"
!      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)

      v1det=0.d0
      j1=(ja+1)/2
      j2=(jb+1)/2
      v1det=v1det+find2(j1,j2)                                       
      if (ja.eq.jb) then
        v1det=v1det+eps(j1)
      end if
      if (ic(j1).eq.2) then
        v1det=-1.d0*v1det
      end if
      return
      end
      function v2det (ja,jb,jc,jd)
      implicit real*8 (a-h,o-z)
      v2det=0.d0
      j1=(ja+1)/2
      j2=(jb+1)/2
      j3=(jc+1)/2
      j4=(jd+1)/2
      v2det=v2det+find4(j1,j2,j3,j4)                                
      if ((mod(ja,2)).eq.(mod(jb,2))) then
        v2det=v2det-find4(j1,j2,j4,j3)                               
      end if
      return
      end

      function v3det (ja,jb,jc,jd,je,jf)                   
      implicit real*8 (a-h,o-z)
#include "prm.h"
!     common /prm/ iperm(0:6,720)

#include "dat.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)

      dimension jj1(6),ii(3)
      jj1(1)=(ja+1)/2
      jj1(2)=(jb+1)/2
      jj1(3)=(jc+1)/2 
      jj1(4)=(jd+1)/2 
      jj1(5)=(je+1)/2 
      jj1(6)=(jf+1)/2 
      ii(1)=mod(ja,2)
      ii(2)=mod(jb,2)
      ii(3)=mod(jc,2)
      v3det=0.d0
      call perms (3)
      do 10 i=1,6  
        iflag=1
        do 20 j=1,3
          if (ii(j).ne.ii(iperm(j,i))) then
            iflag=0
          end if
20      continue
        if (iflag.eq.1) then
          v3det=v3det+iperm(0,i)*find6(jj1(1),jj1(2),jj1(3),
     1    jj1(iperm(1,i)+3),jj1(iperm(2,i)+3),jj1(iperm(3,i)+3))
        end if
10    continue
      if ((ic(jj1(1)).eq.2).and.(ic(jj1(2)).eq.2)
     1.and.(ic(jj1(3)).eq.2)) then
        v3det=-1.d0*v3det
      end if
      return
      end
      function v4det (ja,jb,jc,jd,je,jf,jg,jh)         
      implicit real*8 (a-h,o-z)
      common /prm/ iperm(0:6,720)
      dimension jj1(8),ii(4)
      jj1(1)=(ja+1)/2
      jj1(2)=(jb+1)/2 
      jj1(3)=(jc+1)/2 
      jj1(4)=(jd+1)/2 
      jj1(5)=(je+1)/2 
      jj1(6)=(jf+1)/2 
      jj1(7)=(jg+1)/2 
      jj1(8)=(jh+1)/2 
      ii(1)=mod(ja,2)
      ii(2)=mod(jb,2)
      ii(3)=mod(jc,2)
      ii(4)=mod(jd,2)
      v4det=0.d0
      call perms (4)
      do 10 i=1,24  
        iflag=1
        do 20 j=1,4
          if (ii(j).ne.ii(iperm(j,i))) then
            iflag=0
          end if
20      continue
        if (iflag.eq.1) then
          v4det=v4det+iperm(0,i)*find8(jj1(1),jj1(2),jj1(3),jj1(4),
     1    jj1(iperm(1,i)+4),jj1(iperm(2,i)+4),
     2    jj1(iperm(3,i)+4),jj1(iperm(4,i)+4))                    
          end if
10    continue
      return
      end 

      function v5det (ja,jb,jc,jd,je,jf,jg,jh,ji,jj)
      implicit real*8 (a-h,o-z)
#include "prm.h"
!      common /prm/ iperm(0:6,720)
#include "dat.h"
!      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)

      dimension jj1(10),ii(5)
      jj1(1)=(ja+1)/2
      jj1(2)=(jb+1)/2
      jj1(3)=(jc+1)/2
      jj1(4)=(jd+1)/2 
      jj1(5)=(je+1)/2 
      jj1(6)=(jf+1)/2 
      jj1(7)=(jg+1)/2 
      jj1(8)=(jh+1)/2 
      jj1(9)=(ji+1)/2 
      jj1(10)=(jj+1)/2 
      ii(1)=mod(ja,2)
      ii(2)=mod(jb,2)
      ii(3)=mod(jc,2)
      ii(4)=mod(jd,2)
      ii(5)=mod(je,2)
      v5det=0.d0
      call perms (5)
      do 10 i=1,120
        iflag=1
        do 20 j=1,5
          if (ii(j).ne.ii(iperm(j,i))) then
            iflag=0
          end if
20      continue
        if (iflag.eq.1) then
          v5det=v5det+iperm(0,i)*find10(jj1(1),jj1(2),jj1(3),
     1    jj1(4),jj1(5),jj1(iperm(1,i)+5),jj1(iperm(2,i)+5),
     2    jj1(iperm(3,i)+5),jj1(iperm(4,i)+5),jj1(iperm(5,i)+5))
        end if
10    continue
      if ((ic(jj1(1)).eq.2).and.(ic(jj1(2)).eq.2).and.
     1(ic(jj1(3)).eq.2).and.(ic(jj1(4)).eq.2).and.
     2(ic(jj1(5)).eq.2)) then
        v5det=-1.d0*v5det
      end if
      return
      end

      function v6det (ja,jb,jc,jd,je,jf,jg,jh,ji,jj,jk,jl)
      implicit real*8 (a-h,o-z)
#include "prm.h"
!      common /prm/ iperm(0:6,720)
      dimension jj1(12),ii(6)
      jj1(1)=(ja+1)/2 
      jj1(2)=(jb+1)/2 
      jj1(3)=(jc+1)/2 
      jj1(4)=(jd+1)/2 
      jj1(5)=(je+1)/2 
      jj1(6)=(jf+1)/2 
      jj1(7)=(jg+1)/2 
      jj1(8)=(jh+1)/2 
      jj1(9)=(ji+1)/2 
      jj1(10)=(jj+1)/2 
      jj1(11)=(jk+1)/2 
      jj1(12)=(jl+1)/2 
      ii(1)=mod(ja,2)
      ii(2)=mod(jb,2)
      ii(3)=mod(jc,2)
      ii(4)=mod(jd,2)
      ii(5)=mod(je,2)
      ii(6)=mod(jf,2)
      v6det=0.d0
      n=6
      call perms (n)
      do 10 i=1,720
        iflag=1
        do 20 j=1,6
          if (ii(j).ne.ii(iperm(j,i))) then
            iflag=0
          end if
20      continue
        if (iflag.eq.1) then
          v6det=v6det+iperm(0,i)*find12(jj1(1),jj1(2),jj1(3),jj1(4),
     1    jj1(5),jj1(6),jj1(iperm(1,i)+5),jj1(iperm(2,i)+5),
     2    jj1(iperm(3,i)+5),jj1(iperm(4,i)+5),jj1(iperm(5,i)+5),
     3    jj1(iperm(6,i)+5))                                         
          end if
10    continue
      return
      end
