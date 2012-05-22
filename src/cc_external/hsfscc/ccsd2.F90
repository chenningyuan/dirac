!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dimx(nn,n0,na)
!------------------------------------------------------------------------
!
!    check whether there is enough space of arrays (nn<=n0)
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
#include "implicit.h"
#include "priunit.h"
#else
#include "implicit.h"
#include "priunit.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#endif

      character*8 na

      if (nn.le.n0) then   
        ! ... proper return, no overflow
        return
      else ! treat overflow

        write (LUPRI,10) na,nn,n0
10      format(/,'dimx: dimension overflow at ',a8,1h:,i18,             &
    ', size allowed (allocated):',i18,/)

        write(LUPRI,'(4x,a)') 'List of array consumption:'
        do 20 i=1,nflcur
          write (LUPRI,99) i,namefl(i),nstrt(2,i)-nstrt(2,i-1)
99        format (1x,i3,1x,a8,1x,i20)
20      continue
        call abext('dimx    ')
      end if
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine diveps(namef)
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use dat_mod
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
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#endif

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
          ka=iand(ishft(igg,-24),z'ff')
          kb=iand(ishft(igg,-16),z'ff')
          kc=iand(ishft(igg, -8),z'ff')
          kd=iand(      igg     ,z'ff')
          gen(istrt2+i)=gen(istrt2+i)/                                  &
     (eps(ka)+eps(kb)-eps(kc)-eps(kd))
        endif
20      continue
      else
        do 25 i=1,nsize
        if (gen(istrt2+i).ne.0.d0) then
          kc=iand(ishft(igen(istrt1+i),-8),z'ff')
          kd=iand(      igen(istrt1+i)    ,z'ff')
          gen(istrt2+i)=gen(istrt2+i)/(eps(kc)-eps(kd))
        endif
25      continue
      endif
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dveps6 (name)           
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use dat_mod
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
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "dat.h"
!      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#endif

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
        ka=iand(ishft(ig,-16),  z'ff')
        kb=iand(ishft(ig,-8),   z'ff')
        kc=iand(ig,             z'ff')
        do 20 j=ije1(i-1)+1,ije1(i)
          ig=igen(istrt1+j)
          kd=iand(ishft(ig,-16),  z'ff')
          ke=iand(ishft(ig,-8),   z'ff')
          kf=iand(ig,             z'ff')
          if ((icv(ka).eq.1).and.(icv(kb).eq.1).and.                    &
     (icv(kc).eq.1).and.(icv(kd).eq.1).and.                             &
     (icv(ke).eq.1).and.(icv(kf).eq.1)) then
            if ((ic(ka).eq.3).and.(ic(kb).eq.3).and.                    &
       (ic(kc).eq.3).and.(ic(kd).eq.3).and.                             &
       (ic(ke).eq.3).and.(ic(kf).eq.3)) then
              ijump=1
            else if ((ic(ka).eq.2).and.(ic(kb).eq.2).and.               &
       (ic(kc).eq.2).and.(ic(kd).eq.2).and.                             &
       (ic(ke).eq.2).and.(ic(kf).eq.2)) then
              ijump=1
            else if ((ic(ka).eq.2).and.(ic(kb).eq.2).and.               &
       (ic(kc).eq.3).and.(ic(kd).eq.2).and.                             &
       (ic(ke).eq.2).and.(ic(kf).eq.3)) then
              ijump=1
            else if ((ic(ka).eq.2).and.(ic(kb).eq.2).and.               &
       (ic(kc).eq.3).and.(ic(kd).eq.2).and.                             &
       (ic(ke).eq.3).and.(ic(kf).eq.2)) then
              ijump=1
            end if
          end if
          if (ijump.eq.1) then
            ijump=0
          else
            gen(istrt2+j)=gen(istrt2+j)/                                &
       (eps(ka)+eps(kb)+eps(kc)-eps(kd)-eps(ke)-eps(kf))
          end if
20      continue
10    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine dveps8 (name) 
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use dat_mod
      use cij1_mod
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
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
#endif

      call files_ccsd (name,ntape,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) 'file for dveps8 not found'
        stop                       
      end if
      istrt1=nstrt(1,ntape-1)+2*nij1+4
      istrt2=nstrt(2,ntape-1)
      do 10 i=1,nij1
        id=ij1(i)
        ka=iand(ishft(id,-24),  z'ff')
        kb=iand(ishft(id,-16),  z'ff')
        kc=iand(ishft(id,-8),   z'ff')
        kd=iand(id,             z'ff')
        do 20 j=ije1(i-1)+1,ije1(i)
          ind=igen(istrt1+j)
          ke=iand(ishft(ind,-24),  z'ff')
          kf=iand(ishft(ind,-16),  z'ff')
          kg=iand(ishft(ind,-8),   z'ff')
          kh=iand(ind,             z'ff')
            gen(istrt2+j)=gen(istrt2+j)/                                &
       (eps(ka)+eps(kb)+eps(kc)+eps(kd)                                 &
       -eps(ke)-eps(kf)-eps(kg)-eps(kh))
20      continue
10    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine ehp(ipr)
!---------------------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use cij1_mod
      use dat_mod
      use rsp2_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#include "csize.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "cfiles.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common /csize/ nszv,nszhv,nszhp,nszhq,istrv1,istrv2,
!    1 istrhv1,istrhv2,istrhp1,istrhp2,istrhq1,istrhq2
#include "csize.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /rsp2/  nrsp2,irsp2(500)
#include "rsp2.h"
#endif

!!!!!!!!!!!!!!!!!!!!!!!!TODO: replace with modules !!!!!!
      dimension ev(60),evi(60)
      dimension scr(200),vec(60,60),ham(60,60)
      common /cdet/ ldet(2,400),lsdet(400)
      dimension lldet(2,400)
      dimension ecol(400),ncol(400),iscol(400),idcol(2,400),ntr(12)
      dimension igr(400)

      ! mi: must set for g95
      real*4 time,cputime

      data ntr /1hx,1hs,1hp,1hd,1hf,1hg,1hh,1hi,1hk,1hl,1hm,1hn/
!
      icol=0
      nd=0
      ngr=0
!
      do 270 ib=1,norb
      if (ic(ib).ne.3.or.icv(ib).eq.0) go to 270
!
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
!
!  read folded files.
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
!
!     lsdet(i) is -1 if det(i) has not been treated yet;
!                  0 if det(i) has already been treated;
!                  n if det(i) is now being treated, with ordinal no n.
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
      if (ndt.eq.0) isp=iabs(isign(1,ldet(1,id1))-isign(1,ldet(2,id1))) &
  +1
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
!
!  some matrix elements may have been overlooked. correct.
      do 390 id1=1,ndet
      if (lsdet(id1).le.0) go to 390
      do 380 id2=1,ndet
      if (lsdet(id2).le.0) go to 380
      if (ham(lsdet(id1),lsdet(id2)).eq.0.)                             &
   ham(lsdet(id1),lsdet(id2))=hdet11(ldet,id1,id2)
380   continue
390   continue
      if (ipr.gt.0) write (*,3900) ngr,ndt,(i,lldet(1,i),lldet(2,i),    &
  i=1,ndt)
3900  format(' group',i3,'.',i4,' dets included:',8(i4,':',2i3)/        &
  (29x,6(i4,':',2i3)))
      call dimx(ndt,60,"ndt ")
!
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
!   order and print eigenvalues, vectors.
      do 850 ii=1,ndt
      e=10000.
      do 845 i1=1,ndt
      if (ev(i1).ge.e) go to 845
      e=ev(i1)
      i=i1
845   continue
      print 8400, ii,ev(i),(vec(j,i),j=1,ndt)
8400  format (' energy',i3,1h:,e20.8,'. vector:',7f12.7/                &
  (40x,7f12.7))
!
!     find main determinant contributing
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
!
!     collect equal energy terms.
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
8500  format (1x,10(1h*),' im.energy:',e20.8,'. im. vec:',5f12.7/       &
  (1x,10f12.7))
855   do 860 id=1,ndet
860   if (lsdet(id).gt.0) lsdet(id)=0
      go to 320
!     order and print collected energies.
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
      print 8700, ii,ecol(i),ncol(i),exc,excm,iscol(i),ntr(itr+2),      &
  igr(i),idcol(1,i),idcol(2,i)
8700  format (' energy',i3,e16.8,'. degen',i3,'. exc en=',f7.3,         &
  'ev,',f8.0,'cm-1. term ',i1,a1,'. gr',i3,', major det=',2i3)
880   ecol(i)=10000.
890   continue
      !time=cputime(0.0)
      time=cputime(0.0d0)

      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine eneru2 (egr0,nval,ihp,ipr)
!--------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use cij4_mod
      use cij5_mod
      use cij6_mod
      use dat_mod
      use det_mod
      use idv_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#include "strt.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     external ispin,hmatel
#include "param_ccsd.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "cfiles.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!    1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!    2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
#include "strt.h"
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
!     common /cij6/ nij6,ij6(ndij),ije6(0:ndij),iii6,iiv6
#include "cij6.h"
!     common /dat/  norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!     common /det/  nvlspo,lstvso(100),lstdet(1000),ndets,isplit 
#include "det.h"
!     common /idv/  iv(200),nvalo
#include "idv.h"
#endif

      external ispin,hmatel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO: replace with dyn .alloc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension ecol(200),eicol(200),ncol(200),iscol(200),idcol(2,600), &
  ntr(12)
      dimension ev(60),evi(60)
      dimension scr(200),vec(60,60),veci(60,60),ham(60,60)
      dimension lldet(6,700),lsdet(700),list(6)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   be carefull when replacing with F90 structures !!!!
      equivalence (veci,ham) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!mi          ... watch up possible index leak of lstvso !!!            
              nvlspo=nvlspo+2
              !write(LUPRI,*) 'nvlspo=',nvlspo
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
      !call dimx (nvlspo,100,"nvlspo  ","eneru   ",0)
      call dimx (nvlspo,100,"nvlspo  ")
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
!
! ... start energy calculations.
!
      if (ipr.gt.0) then
        if (nval.lt.12) then
          write (6,2300) nval
2300  format(1h0,80(1h-)/30x,i1,' valence electrons'/30x,19(1h-))
        else
          write (6,2301) nval
2301  format(1h0,80(1h-)/22x,i2,                                        &
 ' sector excitation energies'/30x,22(1h-))
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
      !call dimx(ndets,700,"ndets   ","ener",240)
      call dimx(ndets,700,"ndets   ")
!
!     lsdet is -1 for a determinant not yet treated.
!               n for the nth determinant currently treated.
!               0 for a deteminant already treated.
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
!  some matrix elements may have been overlooked. correct.
      do 390 id1=1,ndets
      if (lsdet(id1).le.0) go to 390
      do 380 id2=1,ndets
      if (lsdet(id2).le.0) go to 380
      if (ham(lsdet(id1),lsdet(id2)).eq.0.d0) then
        ham(lsdet(id1),lsdet(id2))=                                     &
   hmatel(lstdet(id1),lstdet(id2),nval,ihp)
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
      !call dimx(ndt,60,"ndt ","ener",350)
      call dimx(ndt,60,"ndt ")
!     if (ipr.lt.2) go to 830
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
!      if (ipr.eq.0) write (*,8300) tit,(ev(i),i=1,ndt)
!8300  format (1x,a8,' energies:',7e15.8/(1x,8e15.8))
      if (ipr.eq.0) go to 855
!   order and write  eigenvalues, vectors.
      do 850 ii=1,ndt
      e=10000.d0
      do 845 i1=1,ndt
      if (ev(i1).ge.e) go to 845
      e=ev(i1)
      i=i1
845   continue
!   normalize
      fnor=0.
      do 8452 j=1,ndt
8452  fnor=fnor+vec(j,i)**2
      fnor=1.d0/sqrt(fnor)
      do 8455 j=1,ndt
8455  vec(j,i)=vec(j,i)*fnor
      if (ipr.ne.0) write (*,8400) ii,ev(i),(vec(j,i),j=1,ndt)
8400  format (' energy',i3,1h:,e20.8,'. vector:',7f12.7/                &
  (1x,10f12.7))
!
!     find main determinant contributing
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
!
!     collect equal energy terms.
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
!miro: fix idcol index leak...
          if (nvalp-iiii.le.2) then
            idcol(nvalp-iiii,i1)=mod(iddd,256)
          endif
          iddd=iddd/256
        end if
3903  continue
      go to 849
847   continue
848   icol=icol+1
      if (icol.eq.200) write (*,*) ' ******** 200 different energies.', &
  ' others not collected. ***********'
      if (icol.gt.200) go to 849
      iscol(icol)=ksp
      iddd=lstdet(idd)
      do 3904 iiii=0,nvalp-1
        if ((nval.eq.5).or.(nval.eq.6)) then
          idcol(nval-iiii,icol)=ifndso(mod(iddd,16))
          iddd=iddd/16
        else
!mi: idcol first index leak, introduce "if-"fix
          if (nvalp-iiii.le.2) then
            idcol(nvalp-iiii,icol)=mod(iddd,256)
          endif
          iddd=iddd/256
        end if
3904  continue
      ecol(icol)=ev(i)
      eicol(icol)=evi(i)
      ncol(icol)=1
      if (ksp.gt.1) ncol(icol)=2
849   ev(i)=10000.d0
850   if (ipr.ne.0.and.evi(i).ne.0.) write (*,8500) evi(i),             &
  (veci(j,i),j=1,ndt)
8500  format (1x,10(1h*),' im.energy:',e20.8,'. im. vec:',5f12.7/       &
  (1x,10f12.7))
855   iend=1
      do 860 id=1,ndets
      if (lsdet(id).eq.-1) iend=0
860   if (lsdet(id).gt.0) lsdet(id)=0
      if (iend.eq.0) go to 320
!     order and write collected energies.
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
      write (*,8600) ii,ecol(i),ncol(i),exc,excm,iscol(i),ntr(itr+2),   &
  idcol(1,i),idcol(2,i)
8600  format(' energy',i3,e16.8,'. degen',i3,'. ion pot=',f8.3,         &
  'ev,',f9.0,'cm-1. term ',i1,a1,'. major det',6i3)
      if (eicol(i).ne.0.d0) write (*,8650) eicol(i)
8650  format (1x,10(1h*),' im. energy:',e20.8)
      go to 880
877   exc=(ecol(i)-egr)*27.21d0
      excm=exc*8066.d0
      write (*,8700) ii,ecol(i),ncol(i),exc,excm,iscol(i),ntr(itr+2),   &
  idcol(1,i),idcol(2,i)
8700  format (' energy',i3,e16.8,'. degen',i3,'. exc en=',f7.3,         &
  'ev,',f8.0,'cm-1. term ',i1,a1,'. major det',6i3)
      if (eicol(i).ne.0.d0) write (*,8650) eicol(i)
880   ecol(i)=10000.d0
890   continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine expnd(iwork,iwork1,work,work1,nw,nw1)
!--------------------------------------------------------------------------------
!
!   move from work1 (indx iwork1) to work, filling out according to
!   indx iwork
!
!--------------------------------------------------------------------------------

!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension iwork(nw),iwork1(nw1),work(nw),work1(nw1)

      if (nw1.eq.nw) then
        call dcopy(nw,work1,1,work,1)
        return
      endif
      iw1=nw1
      iw=nw
!mi: leak complain
!Fortran runtime error: Array reference out of bounds for array 'iwork1', lower bound of dimension 1 exceeded (0 < 1)
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
!     if (iw.ne.0) go to 10
! mi: fix the iwork1 index leak ...
      if (iw.ne.0.and.iw1.gt.0) go to 10

100   format (' unmatched label ',z8)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine fmax6 (namef,tmax,idmax,indmax)
!--------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cij1_mod
      use cfiles_mod
      use gen_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#endif

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
            indmax=iand(igen(istrt1+ije1(i-1)+j),z'ffffff')
          endif
20      continue
10    continue

      if (getiprint()>= 5) then
        write(LUPRI,*) 'fmax6: tmax,namef=',tmax,namef
        write(LUPRI,*) 'fmax6: idmax,indmax=',idmax,indmax
      endif

      return
      end subroutine fmax6

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine fmax8 (namef,tmax,idmax,indmax)
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cij1_mod
      use cfiles_mod
      use gen_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "cfiles.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
#endif

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

      if (getiprint()>= 5) then
       write(LUPRI,*) 'fmax8: tmax,namef=',tmax,namef
       write(LUPRI,*) 'fmax8: idmax,indmax=',idmax,indmax
      endif

      return
      end subroutine fmax8

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine fmax(namef,tmax,itmax)
!--------------------------------------------------------------------------------
!
! returns tmax, itmax
!
!--------------------------------------------------------------------------------
      use print_level
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
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "cfiles.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
#endif

      call files_ccsd(namef,ntapef,new,1,-1,-1)
      if (new.ne.0) then
        write(LUPRI,*) 'file for fmax not found, namef=',namef
        call abext('fmax    ')
      endif
! find out maximal element and its index
      tmax=0.d0
      itmax=0
!    ... size of array 
      nsize=ije1(nij1)
!    ... starting position of elements
      istrt1=nstrt(1,ntapef-1)+2*nij1+4
!    ... starting position of indexes
      istrt2=nstrt(2,ntapef-1)

      if (getiprint()>= 5) then
       write(LUPRI,*) 'fmax: nsize=',nsize,' istrt1=',istrt1,           &
  ' istrt2=',istrt2 
      endif

      do 20 i=1,nsize

        if (getiprint()>= 12) then
          write(LUPRI,*) 'fmax: i=',i,' gen(istrt2+i)=',gen(istrt2+i)
        endif

      if (abs(gen(istrt2+i)).gt.abs(tmax)) then
        tmax=gen(istrt2+i)
        itmax=igen(istrt1+i)
      endif

20    continue

      if (getiprint()>= 5) then
        write(LUPRI,*) 'fmax: tmax,namef=',tmax,namef
        write(LUPRI,*) 'fmax: itmax=',itmax
      endif

      return
      end subroutine fmax

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function hdet11(ldet,ia,ib)
!---------------------------------------------------------------------------------
!
!
!---------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use dat_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#include "csize.h"
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
!     common /csize/ nszv,nszhv,nszhp,nszhq,istrv1,istrv2,
!    1 istrhv1,istrhv2,istrhp1,istrhp2,istrhq1,istrhq2
#include "csize.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

      dimension ldet(2,*)

!  calculate heff matrix elements
      hdet11=0.
      if (isign(1,ldet(1,ia))-isign(1,ldet(2,ia)).ne.                   &
  isign(1,ldet(1,ib))-isign(1,ldet(2,ib))) return
      if (isign(1,ldet(1,ia)).eq.isign(1,ldet(2,ia))) then
        ii=ior(iabs(ldet(1,ia))*256**3,ior(iabs(ldet(2,ib))*256**2,     &
    ior(iabs(ldet(2,ia))*256,iabs(ldet(1,ib)))))
        hdet11=hdet11+find(ii,gen(istrhq2+1),igen(istrhq1+1),nszhq)
      endif
      if (isign(1,ldet(1,ia)).eq.isign(1,ldet(1,ib))) then
        ii=ior(iabs(ldet(1,ia))*256**3,ior(iabs(ldet(2,ib))*256**2,     &
    ior(iabs(ldet(1,ib))*256,iabs(ldet(2,ia)))))
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function hmatel (id1,id2,nval,ihp)
!--------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use dat_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#include "strt.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij2.h"
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij3.h"
!     common /strt/ istrt11,istrt12,istrt21,istrt22,istrt31,istrt32,
!    1istrt41,istrt42,istrt51,istrt52,istrt61,istrt62,
!    2nsize1,nsize2,nsize3,nsize4,nsize5,nsize6
#include "strt.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
#endif

!!! miro: to replace by dynamic alloc ?
      dimension k1(6),k2(6),k1s(6),k2s(6),                              &
 k1save(6),k2save(6),j1(6),j2(6)

      if (nval.lt.12) then
        hmatel=0.d0
        if (nval.eq.1) then
          isign=1
          if (ispin(id1,1).ne.ispin(id2,1)) return
          hmatel=v1det(id1,id2)       
        else if (nval.eq.2) then
          if (ispin(id1,2).ne.ispin(id2,2)) return
          do 20 i=0,1 
            k1(2-i)=iand(ishft(id1,-8*i),z'ff')
            k2(2-i)=iand(ishft(id2,-8*i),z'ff')
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
            k1(3-i)=iand(ishft(id1,-8*i),z'ff')
            k2(3-i)=iand(ishft(id2,-8*i),z'ff')
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
            k1(4-i)=iand(ishft(id1,-8*i),z'ff')
            k2(4-i)=iand(ishft(id2,-8*i),z'ff')
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
            ii=iand(ishft(id1,-4*i),z'f')
            k1(5-i)=ifndso(ii)
            ii=iand(ishft(id2,-4*i),z'f')
            k2(5-i)=ifndso(ii)
70        continue
          call match6 (k1,k2,5,isign,nm)
          hmatel=v5det(k1(1),k1(2),k1(3),k1(4),k1(5),                   &
     k2(1),k2(2),k2(3),k2(4),k2(5))
          if (nm.ge.1) then
            hmatel=hmatel+v4det(k1(1),k1(2),k1(3),k1(4),                &
                           k2(1),k2(2),k2(3),k2(4))
          end if
          if (nm.ge.2) then
            hmatel=hmatel+v4det(k1(1),k1(2),k1(3),k1(5),                &
                           k2(1),k2(2),k2(3),k2(5))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(3),k2(1),k2(2),k2(3))
          end if
          if (nm.ge.3) then
            hmatel=hmatel+v4det(k1(1),k1(2),k1(4),k1(5),                &
                           k2(1),k2(2),k2(4),k2(5))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(4),k2(1),k2(2),k2(4))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(5),k2(1),k2(2),k2(5))
            hmatel=hmatel+v2det(k1(1),k1(2),k2(1),k2(2))
          end if
          if (nm.ge.4) then
            hmatel=hmatel+v4det(k1(1),k1(3),k1(4),k1(5),                &
                           k2(1),k2(3),k2(4),k2(5))
            hmatel=hmatel+v3det(k1(1),k1(3),k1(4),k2(1),k2(3),k2(4))
            hmatel=hmatel+v3det(k1(1),k1(3),k1(5),k2(1),k2(3),k2(5))
            hmatel=hmatel+v3det(k1(1),k1(4),k1(5),k2(1),k2(4),k2(5))
            do 90 i=2,4
              hmatel=hmatel+v2det(k1(1),k1(1+i),k2(1),k2(1+i))
90          continue
            hmatel=hmatel+v1det(k1(1),k2(1))
          end if
          if (nm.eq.5) then
            hmatel=hmatel+v4det(k1(2),k1(3),k1(4),k1(5),                &
                           k2(2),k2(3),k2(4),k2(5))
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
            ii=iand(ishft(id1,-4*i),z'f')
            k1(6-i)=ifndso(ii)
            ii=iand(ishft(id2,-4*i),z'f')
            k2(6-i)=ifndso(ii)
100       continue
          call match6 (k1,k2,6,isign,nm)
            hmatel=v6det(k1(1),k1(2),k1(3),k1(4),k1(5),k1(6),           &
                    k2(1),k2(2),k2(3),k2(4),k2(5),k2(6))
          if (nm.ge.1) then
            hmatel=hmatel+v5det(k1(1),k1(2),k1(3),k1(4),k1(5),          &
                           k2(1),k2(2),k2(3),k2(4),k2(5))
          end if
          if (nm.ge.2) then
            hmatel=hmatel+v5det(k1(1),k1(2),k1(3),k1(4),k1(6),          &
                           k2(1),k2(2),k2(3),k2(4),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(2),k1(3),k1(4),                &
                           k2(1),k2(2),k2(3),k2(4))
          end if
          if (nm.ge.3) then
            hmatel=hmatel+v5det(k1(1),k1(2),k1(3),k1(5),k1(6),          &
                           k2(1),k2(2),k2(3),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(2),k1(3),k1(5),                &
                           k2(1),k2(2),k2(3),k2(5))
            hmatel=hmatel+v4det(k1(1),k1(2),k1(3),k1(6),                &
                           k2(1),k2(2),k2(3),k2(6))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(3),k2(1),k2(2),k2(3))
          end if
          if (nm.ge.4) then
            hmatel=hmatel+v5det(k1(1),k1(2),k1(4),k1(5),k1(6),          &
                           k2(1),k2(2),k2(4),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(2),k1(4),k1(5),                &
                           k2(1),k2(2),k2(4),k2(5))
            hmatel=hmatel+v4det(k1(1),k1(2),k1(4),k1(6),                &
                           k2(1),k2(2),k2(4),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(2),k1(5),k1(6),                &
                           k2(1),k2(2),k2(5),k2(6))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(4),k2(1),k2(2),k2(4))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(5),k2(1),k2(2),k2(5))
            hmatel=hmatel+v3det(k1(1),k1(2),k1(6),k2(1),k2(2),k2(6))
            hmatel=hmatel+v2det(k1(1),k1(2),k2(1),k2(2))
          end if
          if (nm.ge.5) then
            hmatel=hmatel+v5det(k1(1),k1(3),k1(4),k1(5),k1(6),          &
                           k2(1),k2(3),k2(4),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(3),k1(4),k1(5),                &
                           k2(1),k2(3),k2(4),k2(5))
            hmatel=hmatel+v4det(k1(1),k1(3),k1(4),k1(6),                &
                           k2(1),k2(3),k2(4),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(3),k1(5),k1(6),                &
                           k2(1),k2(3),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(1),k1(4),k1(5),k1(6),                &
                           k2(1),k2(4),k2(5),k2(6))
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
            hmatel=hmatel+v5det(k1(2),k1(3),k1(4),k1(5),k1(6),          &
                           k2(2),k2(3),k2(4),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(2),k1(3),k1(4),k1(5),                &
                           k2(2),k2(3),k2(4),k2(5))
            hmatel=hmatel+v4det(k1(2),k1(3),k1(4),k1(6),                &
                           k2(2),k2(3),k2(4),k2(6))
            hmatel=hmatel+v4det(k1(2),k1(3),k1(5),k1(6),                &
                           k2(2),k2(3),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(2),k1(4),k1(5),k1(6),                &
                           k2(2),k2(4),k2(5),k2(6))
            hmatel=hmatel+v4det(k1(3),k1(4),k1(5),k1(6),                &
                           k2(3),k2(4),k2(5),k2(6))
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
          k1(3-i)=iand(ishft(id1,-8*i), z'ff')
          k1save(3-i)=k1(3-i)
          if (mod(k1(3-i),2).eq.0) then
            k1s(3-i)=-1
          else
            k1s(3-i)=1
          end if
          k2(3-i)=iand(ishft(id2,-8*i),  z'ff')
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
            hmatel=hmatel-ihpsgn*isign*                                 &
       find6(j1(1),j1(2),j2(3),j2(1),j2(2),j1(3))
            if (k2s(1).eq.k2s(2)) then
              hmatel=hmatel+ihpsgn*isign*                               &
         find6(j1(1),j1(2),j2(3),j2(2),j2(1),j1(3))
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
                  if ((k1s(3-i).eq.k1s(3)).and.(k2s(3-j).eq.k2s(3))) then
                    if (i.ne.j) then
                      isign=-1
                    end if
                    hmatel=hmatel+isign*find4(j1(3-i),j2(3),j1(3),j2(3-j))
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
            hmatel=hmatel+ihpsgn*                                       &
       find6(j1(1),j1(2),j2(3),j2(1),j1(3),j2(2))
          end if
          if ((k1s(2).eq.k1s(3)).and.(k2s(1).eq.k2s(3))) then
            hmatel=hmatel-ihpsgn*                                       &
       find6(j1(1),j1(2),j2(3),j2(2),j1(3),j2(1))
          end if
          if ((k1s(1).eq.k1s(3)).and.(k2s(2).eq.k2s(3))) then
            hmatel=hmatel-ihpsgn*                                       &
       find6(j1(2),j1(1),j2(3),j2(1),j1(3),j2(2))
          end if
          if ((k1s(1).eq.k1s(3)).and.(k2s(1).eq.k2s(3))) then
            hmatel=hmatel+ihpsgn*                                       &
       find6(j1(2),j1(1),j2(3),j2(2),j1(3),j2(1))
          end if
        end if
      end if
      return
      end
 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function idec(ihex)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

      i1=iand(ishft(ihex,-24),  z'ff')
      i2=iand(ishft(ihex,-16),  z'ff')
      i3=iand(ishft(ihex,-8) ,  z'ff')
      i4=iand(          ihex ,  z'ff')
      if (i3.eq.0) then
        idec=i1*100+i2
      else
        idec=i1*1000000+i2*10000+i3*100+i4
      endif
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function ifndso(ii)
!---------------------------------------------------------------------------------
!
!
!---------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use idv_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /idv/ iv(200),nvalo
#include "idv.h"
#endif

      if ((ii/2.d0-int(ii/2.d0)).gt.0.1d0) then
        iflag=1
      else
        iflag=0
      end if
      ifndso=2*iv((ii+1)/2)-iflag
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function ifndval(ii)
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use idv_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
!     common /idv/ iv(200),nvalo
#include "idv.h"
#endif

      do 10 j=1,nvalo
        if (ii.eq.iv(j)) then
          ifndval=j             
        end if
10    continue
      return
      end 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function imax (id1,ind1,id2,ind2,iparm)
!--------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      !integer*4 ch1,ch2
      integer :: ch1,ch2

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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine imove(ntapef,ntapet)
!------------------------------------------------------------------
!
!   copy from 1 dim array to another on "igen"
!
!  "ntapef" - file name to be copied into another file "ntapet"
!
!------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cij1_mod
      use cfiles_mod
      use gen_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
!#include "param_ccsd.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#endif

      !integer*4 ch1,ch2
      integer :: ch1,ch2

      nija=igen(nstrt(1,ntapef-1)+1)

      if (getiprint() >= 6) then
        write(LUPRI,'(/,4x,a)') 'Output from "imove"'
        write(LUPRI,*) 'ntapef=',ntapef,' ntapet=',ntapet
        write(LUPRI,*) 'nija=',nija,' nij1=',nij1,' ije1(0)=',ije1(0)
        call print_cijx_stuff
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     call kcopy(nija+1,igen(nstrt(1,ntapef-1)+1),nij1) ! why copying to single variable ?? WRONG ! MUST BE ARRAY !!!
!          using Ephraim's substitution:
      nij1=nija
      do k=1,nija
        ij1(k)=igen(nstrt(1,ntapef-1)+1+k)
      enddo

      if (getiprint() >= 6) then
        write(LUPRI,'(4x,a)') '"imove" - after kcopy (1)'
        call print_cijx_stuff
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call kcopy(nija+1,igen(nstrt(1,ntapef-1)+nija+2),ije1(0))

      if (getiprint() >= 6) then
        write(LUPRI,'(4x,a)') '"imove" - after kcopy (2)'
        call print_cijx_stuff
      endif

      iii1=igen(nstrt(1,ntapef-1)+2*nija+3)
      iiv1=igen(nstrt(1,ntapef-1)+2*nija+4)
      nstrt(1,ntapet)=nstrt(1,ntapet-1)+2*nij1+4+ije1(nij1)
      nstrt(2,ntapet)=nstrt(2,ntapet-1)+ije1(nij1)

      call dimx(nstrt(1,ntapet),nbig,"imove1")
      call dimx(nstrt(2,ntapet),nbig,"imove2")

!   special for scratch files, which start with "s".
           itmp=4H       
           jtmp=z'ffffff00'
           ch1 = iand(itmp,jtmp)
           ch2 = iand (namefl(ntapet), z'000000ff')
           isaa = ior (ch1, ch2)

#if defined VAR_IFORT
      if (isaa.eq.4hs   ) then
#else
      i_4hs=4hs
      if (isaa.eq.i_4hs) then
#endif

        nstrt(1,ntapet)=nstrt(1,ntapet-1)+ndp3+2*ndij+4
        nstrt(2,ntapet)=nstrt(2,ntapet-1)+ndp3
        call dimx(nij1,ndp3,"imove   ")
        call dimx(nstrt(1,ntapet),nbig,"imove1  ")
        call dimx(nstrt(2,ntapet),nbig,"imove2  ")
      endif

      call kcopy(ije1(nij1)+2*nij1+4,igen(nstrt(1,ntapef-1)+1),         &
  igen(nstrt(1,ntapet-1)+1))

      if (getiprint() >= 6) then
        write(LUPRI,'(4x,a)') '"imove" - on exit'
        call print_cijx_stuff
      endif

      return
      end subroutine imove

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function ispin(idet,nval)
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension jj(3)

      if (nval.le.4) then
        ispin=0
        do 10 i=1,nval   
          j=iand(ishft(idet,-8*(i-1)),z'ff')
          if (mod(j,2).eq.1) then
            ispin=ispin+1
          else if (mod(j,2).eq.0) then
            ispin=ispin-1
          end if
10      continue
      else if ((nval.eq.5).or.(nval.eq.6)) then
        ispin=0
        do 20 i=0,nval-1
          j=iand(ishft(idet,-4*i),z'f')
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
        jj(1)=iand(ishft(idet,-16),z'ff')
        jj(2)=iand(ishft(idet,-8), z'ff')
        jj(3)=iand(idet,z'ff')
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine ivsort(ivec,n8)
!--------------------------------------------------------
!
! makes some sorting
!
!--------------------------------------------------------
!     miro: replace with implicit none !!!
      IMPLICIT INTEGER (A-Z)
#include "priunit.h"
      !INTEGER*4  ivec(*),TEMP,N8

      INTEGER  ivec(*),TEMP,N8
      integer lftbit,nx4,n

!
      !if (iopt(1).ge.5) then
!      if (iopt(1).ge.4) then
!       call header("Output from ivsort",-1)
!       write(LUPRI,*) 'ivsort: entering n8=',n8
!      endif

      n=n8
      IF (N .LE. 1) RETURN
!     P0 = AND(N-1,MASK(LEADZ(N-1)+1))
      nx4=n-1
      P0 = lftbit(nx4)
      P = P0
!
   20 Q = P0
      R = 0
      D = P
      I = 0
!
   30 IF (IAND(I,P) .EQ. R) THEN
!$DIR NO_RECURRENCE
         DO 32 J=1,MIN (P,N-I-D)
!           TEMP       = MAX (VEC(J+I),VEC(J+I+D))
!           VEC(J+I)   = MIN (VEC(J+I),VEC(J+I+D))
!           VEC(J+I+D) = TEMP
        if (ishft( ivec(j+i),-1).gt.ishft( ivec(j+i+d),-1)) then
            temp       = ivec(j+i)
            ivec(j+i)   = ivec(j+i+d)
            ivec(J+I+D) = TEMP
        Else If ((ishft(ivec(j+i),-1).eq.ishft(ivec(j+i+d),-1))         &
   .AND. (iabs(ivec(j+i)).gt.iabs(ivec(j+i+d)))) then
             temp       = ivec(j+i)
            ivec(j+i)   = ivec(j+i+d)
            ivec(J+I+D) = TEMP
        endif
   32    CONTINUE
      ENDIF
      I = I + P
      IF ((I+D) .LT. N) GO TO 30
!
      IF (Q .NE. P) THEN
        D = Q-P
        Q = Q/2
        R = P
        I = 0
        GO TO 30
      ENDIF
!
      P = P/2
      IF (P .GT. 4) GO TO 20
!
   70 Q = P0
      R = 0
      D = P
!
   80 DO 82 I=1,P
!$DIR NO_RECURRENCE
      DO 82 J=I,N-D-R,2*P
!        TEMP       = MAX (VEC(J+R),VEC(J+R+D))
!        VEC(J+R)   = MIN (VEC(J+R),VEC(J+R+D))
!        VEC(J+R+D) = TEMP
        if (ishft( ivec(j+r),-1).gt.ishft( ivec(j+r+d),-1)) then
            temp       = ivec(j+r)
            ivec(j+r)   = ivec(j+r+d)
            ivec(J+R+D) = TEMP
        Else If ((ishft(ivec(j+r),-1).eq.ishft(ivec(j+r+d),-1))         &
   .AND. (iabs(ivec(j+r)).gt.iabs(ivec(j+r+d)))) then
            temp       = ivec(j+r)
            ivec(j+r)   = ivec(j+r+d)
            ivec(J+R+D) = TEMP
        endif
   82 CONTINUE
!
      IF (Q .NE. P) THEN
         D = Q - P
         Q = Q / 2
         R = P
         GO TO 80
      ENDIF
!
      P = P / 2
      IF (P .NE. 0) GO TO 70
!
      RETURN
      END SUBROUTINE ivsort

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine kcopy(n,ia,ib)
!--------------------------------------------------------------------------------
! copy n members of integer array "ia" to array "ib"
!     implicit real*8 (a-h,o-z)
!--------------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"

      dimension ia(*),ib(*)
      do 10 i=1,n
10    ib(i)=ia(i)
      return
      end subroutine kcopy

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function lftbit(n)
!
! Bit manipulating function
!
      !integer*4 n,i,ipos
      integer :: n,i,ipos
      !integer*4 nlim(32)
      integer :: nlim(32)

      data nlim/   o'0',o'1',o'2',o'4',                                 &
              o'10',o'20',o'40',                                        &
              o'100',o'200',o'400',                                     &
              o'1000',o'2000',o'4000',                                  &
              o'10000',o'20000',o'40000',                               &
              o'100000',o'200000',o'400000',                            &
              o'1000000',o'2000000',o'4000000',                         &
              o'10000000',o'20000000',o'40000000',                      &
              o'100000000',o'200000000',o'400000000',                   &
              o'1000000000',o'2000000000',o'4000000000',                &
              o'10000000000'/
        ipos=0
        do i=1,32
        if (n.ge.nlim(i)) ipos=i
        enddo
        lftbit=nlim(ipos)
        return
        end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine match6 (k1,k2,nval,isign,nm)
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
!     implicit real *8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension k1(6),k2(6),m1(6),m2(6),mm1(6),mm2(6),                  &
 k1save(6),k2save(6)

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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine move(namef,namet)
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
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
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "cfiles.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#include "gen.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
#endif

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
      call kcopy(nsize,igen(nstrt(1,ntapef-1)+1),                       &
  igen(nstrt(1,ntapet-1)+1))
      call dcopy(ije1(nij1),gen(nstrt(2,ntapef-1)+1),1,                 &
  gen(nstrt(2,ntapet-1)+1),1)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine rle (niter,irleop,name1,name2,name3,name4,ising)
!-----------------------------------------------------------------------
!
!
! miro: please add comments...
!
!-----------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use crle_mod
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
!UK   common /crle/ trle(nrle,4),isize(4)
!     common /crle/ trle1(nrle1),trle2(nrle2),trle3(nrle3),
!    1 trle4(nrle4),isize(4)
#include "crle.h"
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! dynamic alloc !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision r(31,31),z(31),a(31),rcond,rmult
      !integer*4 ipvt(31)
      integer :: ipvt(31)
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
          a(i)=a(i)-ddot(isize(1),trle1(isize(1)+1),1,                  &
     trle1((i+1)*isize(1)+1),1)
41      continue
      end if
      if (iopa(2).eq.1) then
        do 42 i=1,niter
          a(i)=a(i)-ddot(isize(2),trle2(isize(2)+1),1,                  &
     trle2((i+1)*isize(2)+1),1)
42      continue
      end if
      if (iopa(3).eq.1) then
        do 43 i=1,niter
          a(i)=a(i)-ddot(isize(3),trle3(isize(3)+1),1,                  &
     trle3((i+1)*isize(3)+1),1)
43      continue
      end if
      if (iopa(4).eq.1) then
        do 44 i=1,niter
          a(i)=a(i)-ddot(isize(4),trle4(isize(4)+1),1,                  &
     trle4((i+1)*isize(4)+1),1)
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
            trle1((niter+3)*lngth+k)=                                   &
       trle1((i+1)*lngth+k)*trle1(k)
81        continue
          do 71 j=1,niter
            if (j.eq.1) then
              rold=ddot(lngth,trle1((niter+3)*lngth+1),1,               &
         trle1(2*lngth+1),1)
            end if
            rnew=ddot(lngth,trle1((niter+3)*lngth+1),1,                 &
       trle1((j+2)*lngth+1),1)
            r(i,j)=r(i,j)+rold-rnew
            rold=rnew
71        continue
61      continue
      end if
      if (iopa(2).eq.1) then
        lngth=isize(2)
        do 62 i=1,niter
          do 82 k=1,lngth
            trle2((niter+3)*lngth+k)=                                   &
       trle2((i+1)*lngth+k)*trle2(k)
82        continue
          do 72 j=1,niter
            if (j.eq.1) then
              rold=ddot(lngth,trle2((niter+3)*lngth+1),1,               &
         trle2(2*lngth+1),1)
            end if
            rnew=ddot(lngth,trle2((niter+3)*lngth+1),1,                 &
       trle2((j+2)*lngth+1),1)
            r(i,j)=r(i,j)+rold-rnew
            rold=rnew
72        continue
62      continue
      end if
      if (iopa(3).eq.1) then
        lngth=isize(3)
        do 63 i=1,niter
          do 83 k=1,lngth
            trle3((niter+3)*lngth+k)=                                   &
       trle3((i+1)*lngth+k)*trle3(k)
83        continue
          do 73 j=1,niter
            if (j.eq.1) then
              rold=ddot(lngth,trle3((niter+3)*lngth+1),1,               &
         trle3(2*lngth+1),1)
            end if
            rnew=ddot(lngth,trle3((niter+3)*lngth+1),1,                 &
       trle3((j+2)*lngth+1),1)
            r(i,j)=r(i,j)+rold-rnew
            rold=rnew
73        continue
63      continue
      end if
      if (iopa(4).eq.1) then
        lngth=isize(4)
        do 64 i=1,niter
          do 84 k=1,lngth
            trle4((niter+3)*lngth+k)=                                   &
       trle4((i+1)*lngth+k)*trle4(k)
84        continue
          do 74 j=1,niter
            if (j.eq.1) then
              rold=ddot(lngth,trle4((niter+3)*lngth+1),1,               &
         trle4(2*lngth+1),1)
            end if
            rnew=ddot(lngth,trle4((niter+3)*lngth+1),1,                 &
       trle4((j+2)*lngth+1),1)
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
              call daxpy (lngth,a(i),trle1((i+1)*lngth+1),1,            &
         gen(istrt+1),1)
111         continue
          else if (itape.eq.2) then
            do 112 i=1,niter
              call daxpy (lngth,a(i),trle2((i+1)*lngth+1),1,            &
         gen(istrt+1),1)
112         continue
          else if (itape.eq.3) then
            do 113 i=1,niter
              call daxpy (lngth,a(i),trle3((i+1)*lngth+1),1,            &
         gen(istrt+1),1)
113         continue
          else if (itape.eq.4) then
            do 114 i=1,niter
              call daxpy (lngth,a(i),trle4((i+1)*lngth+1),1,            &
         gen(istrt+1),1)
114         continue
          end if
        end if
100   continue
      call prtime ("rle out ",0,0)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine rvmuls(n,a,s,c)
!  ... multiplication of a vector by a scalar
!  ...  c(i)=a(i)*s
! miro: why not using "daxpy" ?
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension a(*),c(*)
      do 10 i=1,n
10    c(i)=a(i)*s
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine rvmulv(n,a,b,c)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension a(*),b(*),c(*)
      do 10 i=1,n
10    c(i)=a(i)*b(i)
      return
      end


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!      subroutine unfold(namef,namet,iline)
      subroutine unfold_ccsd(namef,namet,iline)
!--------------------------------------------------------------------------------
!
!c   expand file namef to namet, turning line iline from valence
!c   to general. (iline=4 is last line).
!
!
!--------------------------------------------------------------------------------
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
      call expnd(igen(nstrt11+1),igen(nstrt21+1),                       &
   gen(nstrt12+1),gen(nstrt22+1),nsize1,nsize2)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine vands(n,ia,is,ic)
!------------------------------------------------------------ 
!
! do and() of array "ia(*)" with scalar "is"
!
!------------------------------------------------------------ 
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension ia(*),ic(*)
      do 10 i=1,n
10    ic(i)=iand(ia(i),is)
      return
      end

      subroutine viorv(n,ia,ib,ic)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension ia(*),ib(*),ic(*)
      do 10 i=1,n
10    ic(i)=ior(ia(i),ib(i))
      return
      end

      subroutine vshfcs(n,ia,j,ic)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
      dimension ia(*),ic(*)
      do 10 i=1,n
10    ic(i)=ishft(ia(i),j)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine vsort (ivec,N8,rlvec)
!--------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------
      !IMPLICIT INTEGER*4 (A-Z)
      IMPLICIT INTEGER (A-Z)
! miro: changed to integer
      !INTEGER*4  ivec(*),TEMP
      INTEGER ::  ivec(*),TEMP
      !integer*4 N8
      integer :: N8

      real*8 rlvec(*),rltemp
!C      integer lftbit,nx4,n
!C
      n=n8
      IF (N .LE. 1) RETURN

!C     P0 = AND(N-1,MASK(LEADZ(N-1)+1))
      nx4=n-1
      P0 = lftbit(nx4)
      P = P0
!C
   20 Q = P0
      R = 0
      D = P
      I = 0
!C
   30 IF (IAND(I,P) .EQ. R) THEN
!C$DIR NO_RECURRENCE
         DO 32 J=1,MIN (P,N-I-D)
!c           TEMP       = MAX (VEC(J+I),VEC(J+I+D))
!c           VEC(J+I)   = MIN (VEC(J+I),VEC(J+I+D))
!c           VEC(J+I+D) = TEMP
!CGJ
!C Now works for 32 bits
        if (ishft( ivec(j+i),-1).gt.ishft( ivec(j+i+d),-1)) then
            temp       = ivec(j+i)
            ivec(j+i)   = ivec(j+i+d)
            ivec(J+I+D) = TEMP
            rltemp     = rlvec(j+i)
            rlvec(j+i) = rlvec(j+i+d)
            rlvec(j+i+d) = rltemp
        Else If ((ishft(ivec(j+i),-1).eq.ishft(ivec(j+i+d),-1))         &
   .AND. (iabs(ivec(j+i)).gt.iabs(ivec(j+i+d)))) then
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
!C
      IF (Q .NE. P) THEN
        D = Q-P
        Q = Q/2
        R = P
        I = 0
        GO TO 30
      ENDIF
!C
      P = P/2
      IF (P .GT. 4) GO TO 20
!C
   70 Q = P0
      R = 0
      D = P
!C
   80 DO 82 I=1,P
!C$DIR NO_RECURRENCE
      DO 82 J=I,N-D-R,2*P
!c        TEMP       = MAX (VEC(J+R),VEC(J+R+D))
!c        VEC(J+R)   = MIN (VEC(J+R),VEC(J+R+D))
!c        VEC(J+R+D) = TEMP
        if (ishft( ivec(j+r),-1).gt.ishft( ivec(j+r+d),-1)) then
            temp       = ivec(j+r)
            ivec(j+r)   = ivec(j+r+d)
            ivec(J+R+D) = TEMP
            rltemp     = rlvec(j+r)
            rlvec(j+r) = rlvec(j+r+d)
            rlvec(j+r+d) = rltemp
        Else If ((ishft(ivec(j+r),-1).eq.ishft(ivec(j+r+d),-1))         &
   .AND. (iabs(ivec(j+r)).gt.iabs(ivec(j+r+d)))) then
            temp       = ivec(j+r)
            ivec(j+r)   = ivec(j+r+d)
            ivec(J+R+D) = TEMP
            rltemp     = rlvec(j+r)
            rlvec(j+r) = rlvec(j+r+d)
            rlvec(j+r+d) = rltemp
        endif
   82 CONTINUE
!C
      IF (Q .NE. P) THEN
         D = Q - P
         Q = Q / 2
         R = P
         GO TO 80
      ENDIF
!C
      P = P / 2
      IF (P .NE. 0) GO TO 70
!C
      RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine vxsort(itmp,ntmp,nt,t,it)
!---------------------------------------------------------------------------------------
!  ...  what does it do ?
!
!  sort according indexes in "itmp", expand to "it", gather "t" accordingly.
!
!  t, it is usually gen, igen
!
!  checks nt,ntmp
!
!  called many times from ccsd3.F
!
!---------------------------------------------------------------------------------------
      use print_level
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
! miro : watch up for array leaks !
      dimension itmp(ntmp),t(nt),it(nt)
!     dimension itmp(*),t(nt),it(nt)

      if (getiprint() >= 5) then
       call header('Output from vxsort',-1)
       write(LUPRI,'(2x,a,2i6)') 'ntmp,nt=',ntmp,nt
      endif
    
      if (ntmp.gt.1) then
         call vsort(itmp,ntmp,t)
         !write(LUPRI,*) 'after vsort called'
      endif
!C      print '(2hGG,Z8)',(itmp(i),i=1,ntmp)

      if (ntmp.gt.nt) then
        write (LUPRI,*) ' vxsort. ntmp gt nt! ',ntmp,nt
        call zeroup(ntmp,itmp,t,1.d-8)
        write (LUPRI,*) ntmp,' integrals left'
        if (ntmp.gt.nt) then
          call zeroup(ntmp,itmp,t,5.d-6)
          write (LUPRI,*) ntmp,' integrals left'
          if (ntmp.gt.nt) then
            write (LUPRI,'(2x,a,2i5)')                                  &
       ' vxsort. ntmp still gt nt ',ntmp,nt
            call abext('vxsort  ')
          endif
        endif
      endif

!     if (ntmp.lt.nt) then
      if (ntmp.lt.nt.and.ntmp.gt.0) then ! mi fix

        !write (LUPRI,'(1x,a,2i5)') 'ntmp.lt.nt !, ntmp,nt:',ntmp,nt
        lt=nt
        ltmp=ntmp

!       ... cycle
! mi: Fortran runtime error: 
!  Array reference out of bounds for array 'itmp', lower bound of dimension 1 exceeded (0 < 1)
!
20      if (itmp(ltmp).eq.it(lt)) then
          t(lt)=t(ltmp)
          lt=lt-1
          ltmp=ltmp-1
        elseif ((ltmp.eq.0).or.                                         &
     (ishft(itmp(ltmp),-1).lt.ishft(it(lt),-1)).or.                     &
     ((ishft(itmp(ltmp),-1).eq.ishft(it(lt),-1)).and.                   &
     (iabs(itmp(ltmp)).lt.iabs(it(lt)))))  then
          t(lt)=0.d0
          lt=lt-1
        elseif (abs(t(ltmp)).lt.2.d-5) then
          ltmp=ltmp-1
        else
         ! ... error branch
         Do i=-20,20
          write (LUPRI,100) lt,nt,ltmp+i,ntmp+i,idec(itmp(ltmp+i)),     &
                       idec(it(ntmp+i)),t(ltmp+i)
100       format(' label not in tmplt!; lt nt ltmp ntmp itmp it t=',    &
     4i7,1x,2i9,e12.4)
         EndDo
         call abext('vxsort  ')   
        endif
!    mi: fix for ltmp
!       if (lt.gt.0) go to 20
        if (lt.gt.0.and.ltmp.gt.0) go to 20
      endif ! of ntmp < nt

      return
      end subroutine vxsort

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine wrtfile(nt,nrec,name)
!----------------------------------------------------------------------------------
!
!c  write internal file "name" on tape nt, skipping nrec records.
!c  disk file structure for tc etc:
!c  1. nij (# of subblocks)
!c  2. ij1(1:nij), subblock id's; ije1(0:nij), subblock ends; iii, iiv.
!c  3. ia(1:ije1(nij)), indices.
!c  4. a(1:ije1(nij)), actual values.
!c    nt.eq.9 = file witout indices (pppp).
!c  // contains all files. file namefl(i) starts at nstrt(i-1)+1 (see /cfiles/).
!
!----------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
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
!     common // gen(nbig)
!     common /i/ igen(nbig)
#endif

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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine wrtrle (name,itape,loc) 
!--------------------------------------------------------------------------------
!
!
!  what does it do ?
!
!--------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use dat_mod
      use crle_mod
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
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "dat.h"
!UK   common /crle/ trle(nrle,4),isize(4)
!     common /crle/ trle1(nrle1),trle2(nrle2),trle3(nrle3),
!    1 trle4(nrle4),isize(4)
#include "crle.h"
#endif

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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine zero (name)
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use dat_mod
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
! mi: original common block
!     common /cij1/ nij, ij(ndij),ije(0:ndij),iii,iiva
! mi : replaces with true cij1 block:
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
#include "dat.h"
#endif

      call files_ccsd (name,ntape,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) 'zero input tape does not exist, ',name
        stop                      
      end if
!     istrt1=nstrt(1,ntape-1)+2*nij+4
      istrt1=nstrt(1,ntape-1)+2*nij1+4
      istrt2=nstrt(2,ntape-1)
!     if (iii.le.333333) then
      if (iii1.le.333333) then
!       do 10 i=1,nij   
        do 10 i=1,nij1
!         i1=iand(ishft(ij(i),-16),'ff'x)
          i1=iand(ishft(ij1(i),-16),z'ff')
!         i2=iand(ishft(ij(i),-8),'ff'x)
          i2=iand(ishft(ij1(i),-8),z'ff')
!         i3=iand(ij(i),'ff'x)
          i3=iand(ij1(i),z'ff')
!         do 20 j=ije(i-1)+1,ije(i)      
          do 20 j=ije1(i-1)+1,ije1(i)      
            i4=iand(ishft(igen(istrt1+j),-16),z'ff')
            i5=iand(ishft(igen(istrt1+j),-8),z'ff')
            i6=iand(igen(istrt1+j),z'ff')
            if ((icv(i1).eq.1).and.(icv(i2).eq.1).and.                  &
       (icv(i3).eq.1).and.(icv(i4).eq.1).and.                           &
       (icv(i5).eq.1).and.(icv(i6).eq.1)) then
              gen(istrt2+j)=0.d0
            end if
20        continue
10      continue
!     else if (iii.gt.333333) then
      else if (iii1.gt.333333) then
!       do 30 i=1,nij   
        do 30 i=1,nij1
!         i1=iand(ishft(ij(i),-24),'ff'x)
          i1=iand(ishft(ij1(i),-24),z'ff')
!         i2=iand(ishft(ij(i),-16),'ff'x)
          i2=iand(ishft(ij1(i),-16),z'ff')
!         i3=iand(ishft(ij(i),-8),'ff'x)
          i3=iand(ishft(ij1(i),-8),z'ff')
!         i4=iand(ij(i),'ff'x)
          i4=iand(ij1(i),z'ff')
!         do 40 j=ije(i-1)+1,ije(i)      
          do 40 j=ije1(i-1)+1,ije1(i)      
            i5=iand(ishft(igen(istrt1+j),-24),z'ff')
            i6=iand(ishft(igen(istrt1+j),-16),z'ff')
            i7=iand(ishft(igen(istrt1+j),-8),z'ff')
            i8=iand(igen(istrt1+j),z'ff')
            if ((icv(i1).eq.1).and.(icv(i2).eq.1).and.                  &
       (icv(i3).eq.1).and.(icv(i4).eq.1).and.                           &
       (icv(i5).eq.1).and.(icv(i6).eq.1).and.                           &
       (icv(i7).eq.1).and.(icv(i8).eq.1)) then
              gen(istrt2+j)=0.d0
            end if
40        continue
30      continue
      end if
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine zeroup(ntmp,itmp,t,zero)
!--------------------------------------------------------------------------------
!
!       moves 'zero' elements & their indexes at the end of the array
!
!--------------------------------------------------------------------------------
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
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
      end subroutine zeroup

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine parse(ii,i1,i2,i3,i4)
!--------------------------------------------------------------------------------
!
! extracts 4 digits from the 4-digit integer number
!     implicit real*8 (a-h,o-z)
!
!--------------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"

      ia=ii
      i4=mod(ia,10)
      ia=ia/10
      i3=mod(ia,10)
      ia=ia/10
      i2=mod(ia,10)
      i1=ia/10
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine parse6 (ii,i1,i2,i3,i4,i5,i6)
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
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
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine perm (i1,i2,isign,nval1)
!--------------------------------------------------------------------------------
!
!                         returns "isign"
!
!--------------------------------------------------------------------------------
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine perms(n)
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
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
              !call perm(k1,k2,isign,3)
              call perm(k1,k2,isign_x,3)
              !iperm(0,ii)=isign
              iperm(0,ii)=isign_x
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
                  call perm (k1,k2,isign_x,4)
                  iperm(0,ii)=isign_x
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
                      call perm(k1,k2,isign_x,5)
                      iperm(0,ii)=isign_x
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
                        if ((i.ne.m).and.(j.ne.m).and.(k.ne.m).and.     &
                   (l.ne.m)) then
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
                          call perm(k1,k2,isign_x,6)    
                          iperm(0,ii)=isign_x
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine prt(name,tit,iparm)
!---------------------------------------------------------------------------------
!
!            prints integrals & indexes ?
!
!  name : 4-byte (integer) character 
!  tit  :  8-byte  (real)  character
!
!  iparm : printout parameter 
!
!---------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use dat_mod
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
!     common // gen(nbig)
!     common /i/ igen(nbig)
!     ... needs own include file ...
!     common /cij1/ nij, ij(ndij), ije(0:ndij), iii, iiv
! mi: replaced
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
#include "dat.h"
!      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif

      dimension a(6),ind(4,6)
!miro: try eliminate non-ascii characters from routine name by declaration
      character(*) tit

      if (iopt(1).ge.5) then
        call header('Output from "prt"',-1)
        write(LUPRI,'(2x,a,a8)') 'tit=',tit
        write(LUPRI,'(2x,a,i5)') 'iparm=',iparm
        write(LUPRI,'(2x,a,a4)') 'name=',name
      endif

      if (iparm.gt.iopt(1)) return

      flim=10.d0**(-iopt(5))
      if (iopt(5).le.0) flim=0.d0

      call files_ccsd(name,num,new,1,-1,-1)

      if (iopt(1).ge.5) then
        write(LUPRI,*) 'prt: flim=',flim
        write(LUPRI,*) 'prt: num=',num
      endif

      if (new.ne.0) then
        write (LUPRI,*) ' prt tried to read nonexistent file'
        call abext('prt     ')
      endif

!     nsize=ije(nij)
      nsize=ije1(nij1)

!     write (LUPRI,110) tit,nsize,iii,iiv,nij,(ij(i),ije(i),i=1,nij)
      write (LUPRI,110)                                                 &
    tit,nsize,iii1,iiv1,nij1,(ij1(i),ije1(i),i=1,nij1)
110   format(' file ',a4,', size iii1 iiv1=',4i7,' ij records:'/        &
  (8(2x,z4,i8,';')))

!     istrt1=nstrt(1,num-1)+2*nij+4
      istrt1=nstrt(1,num-1)+2*nij1+4
      istrt2=nstrt(2,num-1)
      io=0
      do 30 i=1,nsize
      if (abs(gen(istrt2+i)).ge.flim) then
        io=io+1
        a(io)=gen(istrt2+i)
        ind(1,io)=iand(ishft(igen(istrt1+i),-24),z'ff')
        ind(2,io)=iand(ishft(igen(istrt1+i),-16),z'ff')
        ind(3,io)=iand(ishft(igen(istrt1+i),-8) ,z'ff')
        ind(4,io)=iand(           igen(istrt1+i),z'ff')
        if (io.eq.6) then
          write(LUPRI,140) ((ind(jj,ii),jj=1,4),a(ii),ii=1,6)
          io=0
        endif

!       ... mi:  printout/error trap !
!       if (dabs(a(io)).gt.1.0d+7) then
!        write(LUPRI,*) 'error ! big value ! index istrt2+i=',istrt2+i
!        write(LUPRI,*) 'istrt1=nstrt(1,num-1)+2*nij+4,where num=',num
!        write(LUPRI,*) 'a(io)=gen(istrt2+i)=',gen(istrt2+i)
!        call quit('prt: big values of integral!')
!       endif

      endif
30    continue

      if (io.gt.0) then
        write(LUPRI,140) ((ind(jj,ii),jj=1,4),a(ii),ii=1,io)
140     format(1x,6(1x,4i2,'=',d11.5))
      endif

      return
      end subroutine prt

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine prt6 (name,tit,iparm)
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use dat_mod
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
!     common // gen(nbig)
!     common /i/ igen(nbig)
!miro : variables replacement
!     common /cij1/ nij ,ij(ndij) ,ije(0:ndij), iii ,iiv
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
#include "dat.h"
!      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif
      dimension a(4),id(3),ind(3,4)

      if (iparm.gt.iopt(1)) return
      flim=10.d0**(-iopt(5))
      if (iopt(5).le.0) flim=0.d0

      call files_ccsd(name,num,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) ' prt6 tried to read nonexistent file'
        stop                  
      endif
!     nsize=ije(nij)  
      nsize=ije1(nij1)
!     write (6,99) tit,nsize,iii,iiv,nij,(ij(i),ije(i),i=1,nij)
      write (6,99) tit,nsize,iii1,iiv1,nij1,(ij1(i),ije1(i),i=1,nij1)
99    format(' file ',a4,', size iii1 iiv1=',4i9,' ij records:'/        &
  (3(2x,z8,i8,';')))
!     istrt1=nstrt(1,num-1)+2*nij+4
      istrt1=nstrt(1,num-1)+2*nij1+4
      istrt2=nstrt(2,num-1)
!     do 40 i=1,nij
      do 40 i=1,nij1
        iflag=1
        io=0
!       do 30 j=ije(i-1)+1,ije(i)
        do 30 j=ije1(i-1)+1,ije1(i)
          if (abs(gen(istrt2+j)).ge.flim) then
            if (iflag.eq.1) then
!             id(1)=iand(ishft(ij(i),-16),'ff'x)
              id(1)=iand(ishft(ij1(i),-16),z'ff')
!             id(2)=iand(ishft(ij(i),-8),'ff'x)
              id(2)=iand(ishft(ij1(i),-8),z'ff')
!             id(3)=iand(ij(i),'ff'x)
              id(3)=iand(ij1(i),z'ff')
!             write (6,*) ' ' 
              write (6,98) i,(id(k),k=1,3)
98            format ('record no. ',i6,' is ',3i2)
!             write (6,*) ' '
              iflag=0
            end if
            io=io+1
            a(io)=gen(istrt2+j)
            ii=igen(istrt1+j)
            ind(1,io)=iand(ishft(ii,-16),z'ff')
            ind(2,io)=iand(ishft(ii,-8),z'ff')
            ind(3,io)=iand(ii,z'ff')
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine prt8 (name,tit,iparm)
!--------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use dat_mod
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
!     common // gen(nbig)
!     common /i/ igen(nbig)
! miro: for replacement
!     common /cij1/ nij ,ij(ndij), ije(0:ndij),iii,iiv
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
#include "dat.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif

      dimension a(4),id(4),ind(4,4)

      if (iparm.gt.iopt(1)) return
      flim=10.d0**(-iopt(5))
      if (iopt(5).le.0) flim=0.d0
      call files_ccsd (name,num,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) ' prt8 tried to read nonexistent file'
        stop                  
      endif
!     nsize=ije(nij)  
      nsize=ije1(nij1)  
!     write (6,99) tit,nsize,iii,iiv,nij,(ij(i),ije(i),i=1,nij)
      write (6,99) tit,nsize,iii1,iiv1,nij1,(ij1(i),ije1(i),i=1,nij1)
99    format(' file ',a4,', size iii iiv=',4i9,' ij records:'/          &
  (3(2x,z8,i8,';')))
!     istrt1=nstrt(1,num-1)+2*nij+4
      istrt1=nstrt(1,num-1)+2*nij1+4
      istrt2=nstrt(2,num-1)
!     do 40 i=1,nij
      do 40 i=1,nij1
        iflag=1
        io=0
!       do 30 j=ije(i-1)+1,ije(i)
        do 30 j=ije1(i-1)+1,ije1(i)
          if (abs(gen(istrt2+j)).ge.flim) then
            if (iflag.eq.1) then
!             id(1)=iand(ishft(ij(i),-24),'ff'x)
              id(1)=iand(ishft(ij1(i),-24),z'ff')
!             id(2)=iand(ishft(ij(i),-16),'ff'x)
              id(2)=iand(ishft(ij1(i),-16),z'ff')
!             id(3)=iand(ishft(ij(i),-8),'ff'x)
              id(3)=iand(ishft(ij1(i),-8),z'ff')
!             id(4)=iand(ij(i),'ff'x)
              id(4)=iand(ij1(i),z'ff')
              write (6,*) ' ' 
              write (6,98) i,(id(k),k=1,4)
98            format ('record no. ',i6,' is ',4i2)
              write (6,*) ' '
              iflag=0
            end if
            io=io+1
            a(io)=gen(istrt2+j)
            ii=igen(istrt1+j)
            ind(1,io)=iand(ishft(ii,-24), z'ff')
            ind(2,io)=iand(ishft(ii,-16), z'ff')
            ind(3,io)=iand(ishft(ii,-8),  z'ff')
            ind(4,io)=iand(ii,            z'ff')
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine prtime(tit,n,iparm)
!--------------------------------------------------------------------------
!
!  print time using system specific routine
!
!  tit - "title" of the routine, entering as the string among parameters list
!
!  n - hm, like additional specifikation number  (?)
!  
!  iparm - parameter, which must be lower than printing value
!
!--------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#include "csdtime.h"
#else
#include "implicit.h"
!     implicit real*8 (a-h,o-z)
#include "priunit.h"
#include "csdtime.h"
!     common /csdtime/ csdstim,csdltim
#include "dat.h"
!      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif

!mi: warning - sometimes I could have to deactivate it otherwise getting crash in printout
      character(*) tit ! to cut off trash characters

      real*4 cputime
      CHARACTER*8  TIMEX,DATEX*10
      character*12 temp

      ! generalized function
!      real*4 cputime

!       ... return if "param > print level"
      if (iparm.gt.iopt(1)) return

      t=cputime(real(0.0))
!      t=cputime(double(0.0))
!     t = mclock()/1000000.

      diff = t - csdltim
      acum = t - csdstim

!mi: take DIRAC general routine ( gp/gpluuk.F )
      CALL DAYTIME(DATEX,TIMEX)

!      write(LUPRI,*)        ' * tit=',tit
!      write(LUPRI,*)        ' * adjustl(tit)=',adjustl(tit)
!      write(LUPRI,'(a,a8)') 'a8 tit=',tit
!      write(LUPRI,'(a,a8)') 'a8 temp_tit=',temp_tit

!     print 900,tit,acum,diff,n,clockx
!900   format (' time aft ',a8,f10.3,';added ',
!     1 f8.3,' n=',i8,'. ',a24)

!     temp=tit
!     call clean_string(temp)
      !call clean_string(tit)

!      write(LUPRI,*) 'prtime....'
!      write(LUPRI,*) 'len(tit)=',len(tit)
!      write(LUPRI,*) 'trim(tit)=',trim(tit)
!     temp=tit
!     call flush(LUPRI)
!     write(LUPRI,'(a,a9)') 'temp=',temp
      call flush(LUPRI)
!     write(LUPRI,'(a,a9)') 'tit=',tit
!     write(LUPRI,*) 'diff=',diff
!     write(LUPRI,*) 'n=',n
!     write(LUPRI,*) 'datex=',datex
!     write(LUPRI,*) 'timex=',timex

!!!! might be prioblem...
!     print 900,temp,acum,diff,n,datex,timex
      print 900,tit,acum,diff,n,datex,timex
900   format (' time after :',a9,1x,f10.3,'; added ',f10.3,'. n=',i8,'. ',2A12)

      csdltim = t
      call flush (LUPRI)

      return
      end subroutine prtime

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
      subroutine clean_string(temp)
! ... to cut off trash characters      
#include "implicit.h"
#include "priunit.h"
      character(*) temp

      !write(LUPRI,*) 'string length:',sizeof(temp)
!      write(LUPRI,'(a,i4)') 'string length:',len_trim(temo)
      !write(LUPRI,'(a,i4)') 'string length:',len(temp)
      !write(LUPRI,*) 'string length:',len(temp)

      do i=1,len(temp)
        ia=iachar(temp(i:i))
        ! could get problem with printout !
        !write(LUPRI,*) i,"ascii=",ia,' char=',char(ia)

        !write(LUPRI,*) i,"ascii=",ia

       if (.not.( (ia.ge.48.and.ia.le.57).or.                           &
             (ia.ge.65.and.ia.le.90).or.                                &
             (ia.ge.97.and.ia.le.122).or.(ia.eq.32)  )) then 

          !tit(i:i)=char(32,kind=1) ! place 'space'
          !tit(i:i)=" " ! place 'space'
          !s=tit(i:i); write(LUPRI,*) s
          !s=char(32); write(LUPRI,*) s
          !write(LUPRI,*) 'trim:',trim(tit)
          !write(tit(i:i),'(a1)') char(32)
          !write(tit(i:i),'(a)') s
          !tit=trim(tit)
           !temp(i:i)=char(32,kind=1)
           temp(i:i)='\0' ! NULL terminating character
           !write(LUPRI,*) 'added null terminating character !'

       endif
       !write(LUPRI,'(a1,4x,i4)') tit(i:i),iachar(tit(i:i))
      enddo
      !write(LUPRI,*) 'tit=',temp

      end subroutine clean_string

! *************************************************
!     function cputime(x)
!     real*8 function cputime(x)
      real*4 function cputime(x)
!DEC  Real*8 Runsec
      Real*8 Runsec
      real*8 x
!
!NVX  cputime = CPUTIME(0.0)
!RAY  cputime = Second()
!HP9  cputime = RUNSEC()
!DEC  cputime = RUNSEC()
      cputime = secnds(0.0)
!IBM  cputime = MCLOCK()/100.D0
!SGI  cputime = MCLOCK()/100.D0
!
!mi   - old version 
!     t = mclock()/1000000.
      return
      end
! *************************************************

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine rdin
!--------------------------------------------------------------------------------
!
!   input reading routine called from each method routine (formerly executable)
!
!--------------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
#include "implicit.h"
#include "priunit.h"
!     implicit real*8 (a-h,o-z)
#include "dat.h"
!      common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif

      character typ(2,5),cc(200)
      character*8 title(10)
      data typ /'c','C','h','H','p','P','g','G','v','V'/
      logical :: is_input_file,is_open
 
      logical :: dirac_input = .false.
       
      logical :: rdin_called = .false.

      if (rdin_called) then
       write(LUPRI,*) 'rdin was already called, leaving'
       return
      else
       write(LUPRI,*) '...continuing with the rdin procedure...'
      endif

      dirac_input = .true.

      if (.not.dirac_input) then
!!! first reading of the input file - open it !!!
!!!! anyway, input file presence must be checked !!!
      inquire(file='input',exist=is_input_file,opened=is_open)
      ! needed for repeating reading by subsequent executables
      if (.not.is_input_file) then
         call quit('input file not present !')
      else
         write(LUPRI,*) 'input file present, is_open=',is_open
      endif
      rewind(LUCMD)
      endif

!     read (5,'(10a8)')     title
      read (LUCMD,'(10a8)')     title
!     write(6,'(1h1,10a8)') title
      write(LUPRI,'(1h1,10a8)') title

!     read (5,'(40i3)')    (iopt(i),i=1,40)
      read (LUCMD,'(40i3)')    (iopt(i),i=1,40)
      write(LUPRI,'(1x,40i3)') (iopt(i),i=1,40)

      write(LUPRI,'(4x,a)') "iopt(1..40) read in"

      call extract_print_level
      write(LUPRI,'(2x,a,i3,/)') 'Print level (in module):',getiprint()

!     read (5,'(i3/(200a1))') norb,(cc(i),i=1,norb)
      read (LUCMD,'(i3/(200a1))') norb,(cc(i),i=1,norb)
      write(LUPRI,'(4x,a,i3)') "norb read in, norb=",norb

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

!     ... error branch
      write (LUPRI,1050) i,cc(i)
1050  format ('0 **** error in type of orb',i3,3x,a1)
      call abext('rdin    ')
   
!     ... standard output branch 
110   continue
      write (LUPRI,'(i4,'' norbs. types''/(1x,80a1))')                  &
         norb,(cc(i),i=1,norb)

      write(LUPRI,'(/,2x,a)') 'end of rdin, returning'
      call flush (LUPRI)

      rdin_called = .true.

      return
      end subroutine rdin

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine rdfile(nt,nrec,name,ieof)
!-------------------------------------------------------------------------------------
!
!
!  read file "nt", skipping "nrec" records, to internal file number
!  "nflcur"+1, calling it "name".
!
!  disk file structure for tc etc:
!
!  1. nij (# of subblocks)
!  2. ij1(1:nij), subblock id's; ije1(0:nij), subblock ends; iii, iiv.
!  3. ia(1:ije1(nij)), indices.
!  4. a(1:ije1(nij)), actual values.
!
!     ***   april 1990: 1-3 are now in igen, 4 in gen.  ****
!    nt.eq.9 = file without indices (pppp).
!  // contains all files. file namefl(i) starts at nstrt(i-1)+1 (see /cfiles/).
!
!-------------------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
#include "implicit.h"
#include "priunit.h"
!     implicit real*8 (a-h,o-z)
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
! mi: added for control printouts
!#include "dat.h"
!     logical :: is_open,is_present
#endif

! mi: it's useless
!     inquire(unit=nt,exist=is_present,opened=is_open,recl=irecl )

      if (getiprint().ge.5) then
        call header('Output from rdfile',-1)
        write(LUPRI,'(2x,a,i5)') 'file nt=',nt
        write(LUPRI,'(2x,a,i5)') 'nrec=',nrec
        write(LUPRI,'(2x,a,a4)') 'name=',name
      endif

      ! write(LUPRI,*) 'rdile:',is_present,is_open,irecl 

!     if (.not.(is_present.or.is_open.or.irecl==0)) then
!       call quit('rdfile: file tapeXX is NOT present ! ')
!     endif

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

      if (getiprint() >= 8) then
       write(LUPRI,'(2x,a,i7)') 'nflcur=',nflcur
       write(LUPRI,'(2x,a,i7)') 'nstr2=',nstr2
       ! print out extra big values 
       do i=1,nsize
         if (dabs(gen(nstr2+i)).ge.1.0d+5) then
           write(LUPRI,*) nstr2+i,gen(nstr2+i)
         endif
       enddo
      endif

      call dimx(nstrt(1,nflcur),nbig,"rdfile1")
      call dimx(nstrt(2,nflcur),nbig,"rdfile2")
      return

190   if (ieof.eq.0) then
        ieof=1
        return
      else
        write (LUPRI,*) 'rdfile encountered eof'
        write (LUPRI,110) nt,nrec,name
110     format(' nt nrec name=',i5,i4,2x,a8)
        call abext('rdfile  ')
      endif
      end subroutine rdfile

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine rdsym
!----------------------------------------------------------------------------
!
!  continue the reading from  the input file "input"
!
!  gets in symmetry information
!
!  Called from each method (executable)
!
!----------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use sym_mod
      use dpd_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
#include "dpd.h"
!     common /dpd/ nir,igmt(8,8),ici(200)
#endif

      dimension itemp(64) ! 8x8
      data itemp / 1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,     &
            4,3,2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,            &
            7,8,5,6,3,4,1,2,8,7,6,5,4,3,2,1 /

      ! mi: save time - no typing of input, data ready
      !logical :: no_inter_rd=.true.
      logical :: no_inter_rd=.false.
      logical :: rdsym_called=.false.  

      if (rdsym_called) then    
        write(LUPRI,'(/,2x,a,/)') 'rdsym was already read, leaving'
        return
      else
        write(LUPRI,*) '...continuing rdsym procedure...'
      endif

      if (no_inter_rd) then
       nclass=4
       write(LUPRI,'(/,2x,a,i3)')                                       &
  'rdsym: nclass variable set to',nclass
      else
        write(LUPRI,'(/,2x,a,$)') 'Reading "nclass" variable:'
       ! read (*,1000,iostat=ios) nclass
        read (LUCMD,*,iostat=ios) nclass
        if (ios.ne.0) then
          write(LUPRI,*) 'error in reading nclass, read nclass=',nclass
          call quit('rdsym: error in reading nclass!')
        endif
      endif
      write(LUPRI,'(2x,a,i2)') 'read nclass=',nclass

      write(LUPRI,'(/,2x,a)')                                           &
  'Number of of orbitals (20i3) per symmetry(1-nclass):'

      if (no_inter_rd) then
        do i=1,nclass
          iclass(i)=0
        enddo
        iclass(1)=7; iclass(2)=0; iclass(3)=1; iclass(4)=1
        write(LUPRI,'(/,4x,a,/)') 'Hard set orbitals  7 0 1 1 now in !!'
      else

!     read (*,1000,iostat=ios) (iclass(i),i=1,nclass)
!     read (LUCMD,1000,iostat=ios) (iclass(i),i=1,nclass)
      read (LUCMD,*,iostat=ios) (iclass(i),i=1,nclass)
1000  format(20i3)
      if (ios.ne.0) then
        write(LUPRI,*) 'error in reading iclass, so far read iclass(i):'
        write(LUPRI,*) (iclass(i),i=1,nclass)
        call quit('error in reading iclass!')
      endif

      endif

      write (LUPRI,1100) nclass,(iclass(i),i=1,nclass)
1100  format(/,i4,' orb classes; # of orbs ',30i3)

!c-------------------------------------------------
!c  isym is  0 for orbs remaining invariant
!c  isym is  1 for orbs changing sign
!c  isym is -1 for orbs doing neither
!c-------------------------------------------------
      nsym=1

110   continue
      if (no_inter_rd) then
         do i=1,nclass
           isym(i,nsym)=0
         enddo
        ! isym(1,nsym)=9
      else
!       read (*,1000) (isym(i,nsym),i=1,nclass)
        read (LUCMD,1000) (isym(i,nsym),i=1,nclass)
      endif

      !if (isym(1,nsym).ge.9) go to 120
      if (isym(1,nsym).ge.9.or.no_inter_rd) go to 120

      write (LUPRI,1200) nsym,(isym(i,nsym),i=1,nclass)
1200  format(' isym(',i3,')=',30i3)

      nsym=nsym+1
      go to 110
120   nsym=nsym-1


!  ----   fill ibias, ici   ----- 
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
          igmt(i,j)=itemp(ii) ! MI: igmt multiplication table ?
50      continue
40    continue

! MI: nir should be the number of irreducible representations ?
      if (nclass.eq.1) then
        nir=1
      else
        nir=2**nsym
      end if
 
      rdsym_called = .true. 
      write(LUPRI,'(/,2x,a,/)') 'returning from rdsym'
      call flush (LUPRI)

       ! write(LUPRI,'(/,2x,a,/)') 'rdsym was already read, leaving'

      return
      end subroutine rdsym

