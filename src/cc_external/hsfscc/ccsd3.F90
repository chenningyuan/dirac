!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine rdrsp1
!---------------------------------------------------------------------------
!
!
!---------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use rsp1_mod 
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "rsp1.h"
!     common /rsp1/  nrsp1,irsp1(100)
#endif

!!!!!!!!!!!!!!!!!!!!!!!! miro: use dynam.alloc !!!!!!!!!!!!!!!!!!!!!!!!!
      dimension it(2,100)

      read (5,*) nrsp1,((it(j,i),j=1,2),i=1,nrsp1)
      call dimx(nrsp1,100,"rdrsp1")
      do 10 i=1,nrsp1
10    irsp1(i)=ior(ishft(it(1,i),8),it(2,i))
      write (6,100) nrsp1,(it(1,i),it(2,i),i=1,nrsp1)
100   format(i4,' p to r transitions: ',30(2i3,1x))
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine rdrsp2
!---------------------------------------------------------------------------
!
!
!---------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use rsp2_mod 
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "rsp2.h"
!      common /rsp2/  nrsp2,irsp2(500)
#endif

!!!!!!!!!!!!!!!!!!!!!!!! miro: use dynam.alloc !!!!!!!!!!!!!!!!!!!!!!!!!
      dimension it(2,500)

      read (5,*) nrsp2,((it(j,i),j=1,2),i=1,nrsp2)
      call dimx(nrsp2,500,"rdrsp2")
      do 10 i=1,nrsp2
10    irsp2(i)=ior(ishft(it(1,i),8),it(2,i))
      write (6,100) nrsp2,(it(1,i),it(2,i),i=1,nrsp2)
100   format(i4,' r space dets: ',30(2i3,1x))
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine reads (nt,nrec,namein) 
!---------------------------------------------------------------------------
!
! Note:  entering parameter "name" renamed to "namein"
!
!---------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod 
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!      implicit real*8 (a-h,o-z)
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
      igotit=0
      do 20 i=1,nflcur
        if (namefl(i).eq.namein) then
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
        write (LUPRI,*) 'reads encountered end of file'
        call quit('reads encountered end of file!')
      endif
      istrt1=nstrt(1,nfl-1)
      istrt2=nstrt(2,nfl-1)
      igen(istrt1+1)=nij
      read (nt) (igen(istrt1+i),i=2,2*nij+4)
      iiia(nfl)=igen(istrt1+2*nij+3)
      iiva(nfl)=igen(istrt1+2*nij+4)
      namefl(nfl)=namein
      nsize=igen(istrt1+2*nij+2)
      read (nt) (igen(istrt1+2*nij+4+i),i=1,nsize)
      read (nt) (gen(istrt2+i),i=1,nsize)
      nstrt(1,nfl)=istrt1+2*ndij+4+ndp3 
      nstrt(2,nfl)=istrt2+ndp3 
      call dimx(nstrt(1,nfl),nbig,"reads   ")
      call dimx(nstrt(2,nfl),nbig,"reads   ")
      return
      end subroutine reads
!

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine readv(nint,v,ind,nintmx,nfock,fock,ifock,nfockmx)
!---------------------------------------------------------------------------
!
! Basic reading subroutine for integrals & indexes, fock matrix
!
!         nint - total number of 2-el integrals
!         v    - array for 2-el. integrals
!         ind  - array with packed indexes of 2-el. integrals
!
! Called from all unique hsfscc methods
!
!---------------------------------------------------------------------------
     ! miro : deactivate
     ! use get_gamess_mo_ints ! my own module for getting GAMESS new integrals
     use ipack_mod, only: ipack4
#if defined HSFS_DYN_ALLOC
      use dat_mod 
      use sym_mod 
#include "implicit.h"
#include "priunit.h"
#else
!      implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
#endif

      dimension v(*),ind(*),irerdr(200),irerdrv(200),etmp(200)
      dimension fock(*),ifock(*)

      dimension ilet(6)

      logical :: skip_old_read=.true.
      logical :: use_new_gamess=.true.

      data ilet /1hc,1hc,1hh,1hg,1hp,1hv/
      data zero /1.d-8/

! function definition - packing of indexes
    !  ipack(i1,i2,i3,i4)=ior(ior(ishft(i1,24),ishft(i2,16)),ior(ishft(i3,8) ,i4))

      if (iopt(1).ge.5) then
         call header('Output from "readv"',-1)
      endif

!  iopt(19) is 0 for GAMESS, 1 for ACES2(4.92), 2 for ACES2(99)
      if (iopt(19).eq.0) then
        intrc=2725
        write(LUPRI,'(/,2x,a)') 'GAMESS is providing integrals.'
      else if (iopt(19).eq.7) then
        write(LUPRI,'(/,2x,a)') 'DIRAC is providing integrals.'
      else ! (iopt(19).ne.0) then
         call quit('miro: aces ii not YET supported !')
         intrc=600
      endif

!   =================================
!     ==== Old reading - skip  ====
!   =================================
      nort = - 1 ! mi-to satisfy runtime check
      if (.not.skip_old_read) then

      rewind(12,iostat=ios)
      if (ios.ne.0) then
       write(*,*) 'failed rewing of unit 12 !'
       call  abext('readv : failed rewing of unit 12 !')
      endif

!     ... read number of orbitals and orbital energies
      read (12,iostat=ios) nort,(eps(i),i=1,nort)
      if (ios.ne.0) then
        write(LUPRI,'(/,4x,a)')                                         &
   'failed read (12) nort,(eps(i),i=1,nort)!'
        call  abext('readv : failed read (12) nort,(eps(i),i=1,nort)')
      endif

        write(LUPRI,*) 'readv - nort=',nort
        CALL FLSHFO(LUPRI)

      endif

       write(LUPRI,*) 'readv - nort,norb,skip_old_read=',               &
                   nort,norb,skip_old_read

      if (nort.ne.norb.and..not.skip_old_read) then
        write (*,*) 'discrepancy in norb: ',norb,nort 
        call abext('readv: nort.ne.norb ! ')
      end if

!  changes for aces99
      if (iopt(19).eq.2) then
        ! read # of orbs pre irreducible reps or indexes ?
        read (12) (irerdr(i),i=1,norb)
        do i=1,norb
          irerdrv(irerdr(i))=i
          etmp(i)=eps(irerdr(i))
        enddo
        do i=1,norb
          eps(i)=etmp(i)
        enddo
      endif

!     ... printout - orb.energies per symmetries
!      do 340 icl=1,nclass
!      write (*,3300) icl
!3300  format(' irrep ',i2)
!      write (*,3400) (io,ilet(2*ic(io)-1+icv(io)),eps(io),
!     1 io=ibias(icl)+1,ibias(icl+1))
!3400  format(3('orb ',i3,2x,a1,', e=',1pe12.5,:,';  '))
!340   continue

      if (.not.skip_old_read) then
        if (iopt(19).eq.0) read (12)
      endif

      nint=0 ! initialization for the loop

      if (.not.use_new_gamess) then

120   continue
!
      if(iopt(19).eq.0) then ! GAMESS
       if (.not.skip_old_read) then
        read (12) nn,(ind(nint+i),i=1,2725),                            &
   (v(nint+i),i=1,2725)
       endif
      else ! ACESS II
!       write(*,*) 'bef read'
        read (12,end=135,err=199) (v(nint+i),i=1,600),(ind(nint+i),     &
   i=1,600),nn
!       write(*,*) 'aft read. nn=',nn
        if (nn.le.0) go to 135
      endif
!
      call dimx(nint+intrc,nintmx,"nint")
!--------------------------------------------------------
!         UNPACK, ELIMINATE INTS WITH CORE ORBS
!--------------------------------------------------------
      nnint=nint
      do 130 i=nnint+1,nnint+abs(nn)

      if (abs(v(i)).lt.zero) go to 130

      J1=iand(ishft(ind(i),-24),z'ff')
      J2=iand(ishft(ind(i),-16),z'ff')
      J3=iand(ishft(ind(i),-8) ,z'ff')
      J4=iand(ind(i)           ,z'ff')

      if(iopt(19).eq.2) then ! ACESII
        j1=irerdrv(j1)
        j2=irerdrv(j2)
        j3=irerdrv(j3)
        j4=irerdrv(j4)
      endif

!     miro: what is ic ?
      if((ic(j1).eq.1).or.(ic(j2).eq.1).or.(ic(j3).eq.1).or.            &
                                       (ic(j4).eq.1)) go to 130

      nint=nint+1 ! increase index

      v(nint)=v(i)
      ind(nint)=ior(ishft(j1,24),ior(ishft(j3,16),ior(ishft(j2,8),j4)))

130   continue

      if (nn.gt.0) go to 120 ! loop

199   write(*,*) 'Record in tape12 too short'
      call abext('12short')

135   continue

      else if (iopt(19).ne.7) then

        call quit('GAMESS connection deactivated for now ...')

        ! ... get in new GAMESS integrals ... 
        !call get_gamess_ints(nort,eps,nn,ind,v,nint)
        if (iopt(1).ge.5) then ! control prinout
          write(LUPRI,'(/,2x,a,i5)')                 &
     'after call get_gamess_ints, nint=',nint
          write(LUPRI,*) 'after, nort=',nort
          write(LUPRI,*) '         nn=',nn
          write(LUPRI,*) 'eps: ',(eps(i),i=1,nort)
        endif

          ! very important - use switch to DIRAC notation of indexes
          do i=1, nint
            J1=iand(ishft(ind(i),-24),z'ff')
            J2=iand(ishft(ind(i),-16),z'ff')
            J3=iand(ishft(ind(i),-8) ,z'ff')
            J4=iand(ind(i)           ,z'ff')

            !ipckindx=ipack(j1,j2,j3,j4)
            ipckindx=ipack4(j1,j3,j2,j4) ! dirac notation !
            ind(i)=ipckindx

           if (iopt(1).ge.4) then ! control printout
            write(LUPRI,'(i4,a,f18.12,a,2i11,a,4i3)')                   &
       i,'. int:',v(i),' pck in:',ind(i),ipckindx,                      &
       ' indx:',J1,J2,J3,J4
           endif

          enddo
      else 
        ! ... get DIRAC integrals
        call get_dirac_ints(nort,eps,nn,ind,v,nint) 
      endif

!-----------------------------------------------------------------------------------------------
!
!                  place here the printout - orbital energies per symmetries
!
!-----------------------------------------------------------------------------------------------
      WRITE(LUPRI,'(/,2X,A)')                                           &
   'Orbitals per symmetry, their irreps, energies...'
      call flush (LUPRI)

      do 340 icl=1,nclass

      write (*,3300) icl
3300  format(4x,'* irrep ',i2)

      do io=ibias(icl)+1,ibias(icl+1)

!     write (*,3400) (io,ilet(2*ic(io)-1+icv(io)),eps(io),
!    1 io=ibias(icl)+1,ibias(icl+1))

       ilet_indx = 2*ic(io)-1+icv(io)

       if (ilet_indx.le.0.or.ilet_indx.gt.6) then
         write(LUPRI,*)
         write(LUPRI,*) 'error, io,ic(io),icv(io): ',                   &
                    io,ic(io),icv(io)
         write(LUPRI,*) 'ibias(icl)+1,ibias(icl+1): ',                  &
                    ibias(icl)+1,ibias(icl+1)
         write(LUPRI,*) ' ilet_indx=',ilet_indx
         call flush(LUPRI)
          write(LUPRI,*) 'Check the order of MO !!!  '
         call quit('ilet_indx error!')
       endif
       write (*,3400) io,ilet(ilet_indx),eps(io)
       call flush (LUPRI)
      enddo

3400  format(1x,'orb.: ',i3,2x,a1,', orb.energy=',f18.10)

340   continue

!------------------------------------------------------------
      write (LUPRI,'(/,4x,i6,a,/)')                                     &
       nint,' integrals and indices read in !'

      if (iopt(8).eq.0) then
        nfock=0
        write (LUPRI,'(/,a,/)')                                         &
  'Cannonical Hartree-Fock orbitals'//                                  &
  ' (-orbital energies needed only) are used, leaving.'

      !  call quit('end in ccsd3.f')
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
      i1=iand(ifock(i),z'ff')
      i2=ishft(ifock(i),-8)
      if (i1.ne.i2) go to 250
      if (iopt(8).eq.2) eps(i1)=fock(i)
      fock(i)=fock(i)-eps(i1)
250   continue
      do 260 i=1,nfock
      i1=iand(ifock(i),z'ff')
      i2=ishft(ifock(i),-8)
!      write (*,2600) i,i1,i2,fock(i)
! 2600  format(' Fock element no',i3,2x,2i3,f15.6)
260   continue
      go to 440
! ... modification of fock, orbital energies
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
        write (*,3400) (io,ilet(2*ic(io)-1+icv(io)),eps(io),            &
   io=ibias(icl)+1,ibias(icl+1))
450   continue

      return
      end subroutine readv

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine files_ccsd(name,num,new,ijnum,iii,iiv)
!-----------------------------------------------------------------------------------
! 
! FYI: previously known as  "subroutine files(name,num,new,ijnum,iii,iiv)", name changed
! due to overlaping with other subroutine
!
! PURPOSE:
!  subroutine finds # of file named NAME; if not found,
!  increases "nflcur" by 1 and starts new file.
!
!  nflcur set to 0 at start of program.
!
!  new is 1 for new file, -1 for new file with indices copied from
!  another file, 0 for existing file. 
!
!  ijnum=0 used for existing
!  scratch file, i.e. file to be completely rewritten.
!
!  Not dealing with true files, here file=diagram;
!  organizes room for particular diagram on array gen/igen
!
!  Subroutine fills common blocks cij[1-6]
!
!  returns num
!
! miro TODO: make only one return with control printout before
!
!-----------------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cij1_mod 
      use cij2_mod 
      use cij3_mod 
      use cij4_mod 
      use cij5_mod 
      use cij6_mod 
      use cfiles_mod 
      use gen_mod 
      use hsfs_parameters
 implicit none
#include "priunit.h"
 integer :: name, num, new, ijnum, iii, iiv, i
 integer :: i_4hs, isaa, itmp, jtmp, k, nija
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij3.h"
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij4.h"
!     common /cij4/ nij4,ij4(ndij),ije4(0:ndij),iii4,iiv4
#include "cij5.h"
!     common /cij5/ nij5,ij5(ndij),ije5(0:ndij),iii5,iiv5
#include "cij6.h"
!     common /cij6/ nij6,ij6(ndij),ije6(0:ndij),iii6,iiv6
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#endif

      !integer*4 ch1,ch2
      integer :: ch1,ch2

      if (getiprint() >= 14) then
        call header("Output from files_ccsd",-1)
        write(LUPRI,'(4x,a)')       '>>>>> entering data:'
        write(LUPRI,'(2x,a,i7)')    'nflcur=',nflcur
        write(LUPRI,'(2x,a,a8)')    'name=',name
        write(LUPRI,'(2x,2(a,i12))') 'num=',num,' ijnum=',ijnum
        write(LUPRI,'(2x,a,2i5)')   'iii,iiv=',iii,iiv
        call print_cijx_stuff
      endif

!      ... skips filling common blocks variables (no diagram ready now)
      if (nflcur.eq.0) go to 20 

! namefl(:) - array with all diagram names
      num=1000
      do 10 i=1, nflcur
10    if (name.eq.namefl(i)) num=i

      if (num.le.nflcur) then
           new=0 ! diagram already exists

!  next lines for case where temp file reopened with same name,
!  different type (must be sn files). iii=-1 means ignore iii.
      if (iii.gt.-1.and.((iii.ne.iiia(num)).or.(iiv.ne.iiva(num)))) then
           if (getiprint() >= 5) then
             write(LUPRI,*) '1) num,iiia(num)=',num,iiia(num),'...iiva(num)=',iiva(num)
           endif
! Gal 95
           ! take care of "    " as Hollerith constant, "4H    "
           itmp=transfer("    ",itmp) 
           jtmp=z'ffffff00'
           ch1 = iand(itmp,jtmp)
           ch2 = iand(namefl(num),z'000000ff')
           isaa = ior(ch1,ch2)
           if (getiprint() >= 5) then
             ! take care of "s   " as Hollerith constant
             i_4hs=transfer("s   ",i_4hs)
             write(LUPRI,*) "isaa.ne.i_4hs ?",(isaa.ne.i_4hs)
             write(LUPRI,*) "ch1=",ch1," ch2=",ch2
             write(LUPRI,*) "itmp=",itmp," jtmp=",jtmp
           endif

          ! take care of "s   " as Hollerith constant "4hs   "
          i_4hs=transfer("s   ",i_4hs)
          if (isaa.ne.i_4hs) then
            write(LUPRI,*) 
            write(LUPRI,*) "error exit ! isaa, i_4hs=",isaa,i_4hs
            write(LUPRI,*) 'ch1=',ch1,' ch2=',ch2
            write(LUPRI,*) 'itmp=',itmp,' jtmp=',jtmp
            write (LUPRI,120) name,iii,iiv,iiia(num),iiva(num)
120         format(' file ',a8,' with iii iiv',2i5,' previously',2i5)
            call abext('files   ') 
          endif
      endif

       if (iii.gt.-1.and.((iii.ne.iiia(num)).or.(iiv.ne.iiva(num))))then
          if (getiprint() >= 5) then
             write(LUPRI,*) '2) num,iiia(num)=',num,iiia(num),'...iiva(num)=',iiva(num)
          endif

          new=1 ! new file
          iiia(num)=iii
          iiva(num)=iiv
          do 15 i=1,nflcur
          if ((i.ne.num).and.(iiia(i).eq.iii).and.(iiva(i).eq.iiv)) then
            new=-1
            call imove(i,num)
            if (getiprint() >= 14) then
              write(LUPRI,*) 'files_ccsd exit 1; nflcur,num=',nflcur,num
              call print_cijx_stuff
            endif
            return
          endif
15        continue

          if (getiprint() >= 14) then
             write(LUPRI,*) 'files_ccsd exit 2; nflcur,num=',nflcur,num
             call print_cijx_stuff
          endif
          return

        endif

        if (ijnum.eq.0) then 
          if (getiprint() >= 14) then
             write(LUPRI,*) 'files_ccsd exit 3; nflcur,num=',nflcur,num
             call print_cijx_stuff
          endif
          return
        endif

        nija=igen(nstrt(1,num-1)+1)

!         ...  filling  cij* common blocks with diagrams created earlier
        if (ijnum.eq.1) then
         ! this is wrong 
!         call kcopy(nija+1,igen(nstrt(1,num-1)+1),nij1) ! fill nij1,ij1(:)
          nij1=nija ! Ephraim's substitution
          do k=1,nija
            ij1(k)=igen(nstrt(1,num-1)+1+k)
          enddo
          call kcopy(nija+1,igen(nstrt(1,num-1)+nija+2),ije1(0)) ! fill ije1(:)
          iii1=igen(nstrt(1,num-1)+2*nija+3)
          iiv1=igen(nstrt(1,num-1)+2*nija+4)

          if (getiprint().ge.14) then
            write(LUPRI,*) 'ijnum=1, nij1=',nij1
            write(LUPRI,*) 'num=',num
            write(LUPRI,*) 'nija,iii1,iiv1:',nija,iii1,iiv1
          endif

        elseif (ijnum.eq.2) then
!         call kcopy(nija+1,igen(nstrt(1,num-1)+1),nij2) ! wrong !
          nij2=nija ! Ephraim's substitution
          do k=1,nija
            ij2(k)=igen(nstrt(1,num-1)+1+k)
            if (getiprint().ge.15) then
              write(LUPRI,*) nstrt(1,num-1)+1+k,                        &
         '. igen(nstrt(1,num-1)+1+k)=',                                 &
         igen(nstrt(1,num-1)+1+k)
            endif
          enddo
          call kcopy(nija+1,igen(nstrt(1,num-1)+nija+2),ije2(0))
          iii2=igen(nstrt(1,num-1)+2*nija+3)
          iiv2=igen(nstrt(1,num-1)+2*nija+4)

          if (getiprint().ge.14) then
            write(LUPRI,*) 'ijnum=2, nij2=',nij2
            write(LUPRI,*) 'num=',num
            write(LUPRI,*) 'nija,iii2,iiv2:',nija,iii2,iiv2
          endif

        elseif (ijnum.eq.3) then
!         call kcopy(nija+1,igen(nstrt(1,num-1)+1),nij3)
          nij3=nija ! Ephraim's substitution
          do k=1,nija
            ij3(k)=igen(nstrt(1,num-1)+1+k)
          enddo
          call kcopy(nija+1,igen(nstrt(1,num-1)+nija+2),ije3(0))
          iii3=igen(nstrt(1,num-1)+2*nija+3)
          iiv3=igen(nstrt(1,num-1)+2*nija+4)

          if (getiprint().ge.14) then
            write(LUPRI,*) 'ijnum=3, nij3=',nij3
            write(LUPRI,*) 'num=',num
            write(LUPRI,*) 'nija,iii3,iiv3:',nija,iii3,iiv3
          endif

        elseif (ijnum.eq.4) then
         !call kcopy(nija+1,igen(nstrt(1,num-1)+1),nij4) ! wrong !
          nij4=nija ! Ephraim's substitution
          do k=1,nija
            ij4(k)=igen(nstrt(1,num-1)+1+k)
          enddo
          call kcopy(nija+1,igen(nstrt(1,num-1)+nija+2),ije4(0))
          iii4=igen(nstrt(1,num-1)+2*nija+3)
          iiv4=igen(nstrt(1,num-1)+2*nija+4)

          if (getiprint().ge.14) then
            write(LUPRI,*) 'ijnum=4, nij4=',nij4
            write(LUPRI,*) 'num=',num
            write(LUPRI,*) 'nija,iii4,iiv4:',nija,iii4,iiv4
          endif

        elseif (ijnum.eq.5) then
!         call kcopy(nija+1,igen(nstrt(1,num-1)+1),nij5)
          nij5=nija ! Ephraim's substitution
          do k=1,nija
            ij5(k)=igen(nstrt(1,num-1)+1+k)
          enddo
          call kcopy(nija+1,igen(nstrt(1,num-1)+nija+2),ije5(0))
          iii5=igen(nstrt(1,num-1)+2*nija+3)
          iiv5=igen(nstrt(1,num-1)+2*nija+4)

          if (getiprint().ge.14) then
            write(LUPRI,*) 'ijnum=5, nij5=',nij5
            write(LUPRI,*) 'num=',num
            write(LUPRI,*) 'nija,iii5,iiv5:',nija,iii5,iiv5
          endif

        elseif (ijnum.eq.6) then
!         call kcopy(nija+1,igen(nstrt(1,num-1)+1),nij6)
          nij6=nija ! Ephraim's substitution
          do k=1,nija
            ij6(k)=igen(nstrt(1,num-1)+1+k)
          enddo
          call kcopy(nija+1,igen(nstrt(1,num-1)+nija+2),ije6(0))
          iii6=igen(nstrt(1,num-1)+2*nija+3)
          iiv6=igen(nstrt(1,num-1)+2*nija+4)

          if (getiprint().ge.14) then
            write(LUPRI,*) 'ijnum=6, nij6=',nij6
            write(LUPRI,*) 'num=',num
            write(LUPRI,*) 'nija,iii6,iiv6:',nija,iii6,iiv6
          endif
        else
          write (LUPRI,*) 'ccsd_files called with incorrect ijnum, ',ijnum
          call abext('ccsd_files   ')
        endif

        if (getiprint() >= 14) then
          write(LUPRI,*) 'files_ccsd exit 4; nflcur,num=',nflcur,num
          call print_cijx_stuff
        endif
        return

      endif

!  new file when nflcur is 0
20    nflcur=nflcur+1
      namefl(nflcur)=name
      num=nflcur

      iiia(nflcur)=iii
      iiva(nflcur)=iiv

      new=1

      if (iii.eq.-2) then
        if (getiprint() >= 14) then
           write(LUPRI,*) 'files_ccsd exit 5; nflcur,num=',nflcur,num
           call print_cijx_stuff
        endif
        return ! return
      endif

      if (nflcur.eq.1) go to 40 ! skip the cycle

!     ... cycle ...
      do 30 i=1, nflcur-1
      if ((iiia(i).eq.iii).and.(iiva(i).eq.iiv)) then
        new=-1
        call imove(i,num)

!       write (6,*) ' transf indx from ',i,' to ',num
        if (getiprint() >= 14) then
         write (LUPRI,*) ' transf indx from ',i,' to ',num
        endif

        go to 40
      endif
30    continue

40    continue ! bye step

      if (getiprint() >= 14 ) then
        write(LUPRI,*) 'files_ccsd exit 6; nflcur,num=',nflcur,num
        write(LUPRI,*) 'files_ccsd exit: iii,iiv=',iii,iiv
        call print_cijx_stuff
      endif

      return

      end subroutine files_ccsd

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine reorder(namef,namet,iord,iparm)
!---------------------------------------------------------------------------------------------------
!
!  iparm - print level
!
!  namet, namef
!
! change characteristics of the diagram, for example 2143 as iord means changing order of indexes
!
!---------------------------------------------------------------------------------------------------
      use print_level
      !use ipack_mod, only: ipack4
#if defined HSFS_DYN_ALLOC
      use ctmp_mod 
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
#include "ctmp.h"
!      common /ctmp/ tmp(ndp4),itmp(ndp4)
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
#endif

      dimension j(4),msk(4)
      ! ff - place where indexes are placed
      data msk /z'ff000000',z'ff0000',z'ff00',z'ff'/

!     character(*) namet, namef ! mi fix - data specification for prtime
      character*4 namet, namef ! mi fix - data specification for prtime

      if (getiprint() >= 5) then
        call header('Output from "reorder"',-1)
        write(LUPRI,*) 'namef=',namef
        write(LUPRI,*) 'namet=',namet
        write(LUPRI,*) 'iord=',iord
        write(LUPRI,*) 'iparm=',iparm
        ! call check_ctmp_mod
         call print_cijx_stuff
      endif

      call files_ccsd(namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (6,10) namef
10      format(' reorder input tape ',a6,' does not exist.')
        call abext('reorder ')
      endif

      if (getiprint() >= 5) then
         write(LUPRI,'(4x,a)') 'reorder: cijn printout:'
         call print_cijx_stuff
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

!miro: "lngth2+600" could cause index leak for "itmp" array ???
! todo: replace 600 with  some internal parameter ...
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

!     ... acting on gen !!!
      call dcopy(lngth2,gen(istrt22+1),1,gen(istrt12+1),1)
      call vxsort(itmp,lngth2,lngth1,gen(istrt12+1),igen(istrt11+1))

      ! ... print out ...
      !call prt(namet,namet,iparm)
      call prt(namet,"reorder",iparm)
      call prtime(namet,ije1(nij1),1)

      return
      end subroutine reorder

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine reordx (namef,namet,iord)
!---------------------------------------------------------------------------
!
! reordering of indexes with different kind of order
!
!---------------------------------------------------------------------------
      use print_level
      use ipack_mod, only: ipack4
#if defined HSFS_DYN_ALLOC
      use cfiles_mod 
      use gen_mod 
      use cij1_mod 
      use cij2_mod 
      use sym_mod 
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
#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#include "dpd.h"
!     common /dpd/ nir,igmt(8,8),ici(200)
#endif

      dimension ia(8),ja(8)

! ... to be replaced by include file ...
      !ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

      if (getiprint() >= 5) then
        call header('Output from "reordx"',-1)
        write(LUPRI,*) 'namef=',namef
        write(LUPRI,*) 'namet=',namet
        write(LUPRI,*) 'iord=',iord
        !call check_modules
      endif

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
        call quit('reordx input tape name2 does not exist')
      end if
      call parse8 (iord,i1,i2,i3,i4,i5,i6,i7,i8)
      call parse8 (iii2,ia(1),ia(2),ia(3),ia(4),ia(5),ia(6),ia(7),ia(8))
      call parse8 (iiv2,ja(1),ja(2),ja(3),ja(4),ja(5),ja(6),ja(7),ja(8))
      iii1=10000000*ia(i1)+1000000*ia(i2)+100000*ia(i3)+10000*ia(i4)+   &
               1000*ia(i5)+100*ia(i6)+10*ia(i7)+ia(i8)
      iiv1=10000000*ja(i1)+1000000*ja(i2)+100000*ja(i3)+10000*ja(i4)+   &
               1000*ja(i5)+100*ja(i6)+10*ja(i7)+ja(i8)
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
          ia(4-j)=iand(ishft(itmp(i),-8*j),z'ff')  ! unpack 
          ia(8-j)=iand(ishft(itmp(ndp3+i),-8*j),z'ff')
40      continue
        ! another packing with different order
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
              ii=imax(itmp(j),itmp(ndp3+j),           &
                      itmp(j+1),itmp(ndp3+j+1),1)
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine reord6 (namef,namet,iord)
!---------------------------------------------------------------------------
!
! used for triples
!
! Entering "namef" - initial diagram 
!          namet - final, reordered diagram
!          iord - type of reordering, like example above
!
! Example: call reord6("s1  ","t3p ",126345)
!
!---------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cfiles_mod 
      use gen_mod 
      use cij1_mod 
      use cij2_mod 
      use sym_mod 
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
#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#include "dpd.h"
!     common /dpd/ nir,igmt(8,8),ici(200)
#endif

      dimension ia(6),ja(6),msk(4),igot(8)
      data msk /z'ff000000',z'ff0000',z'ff00',z'ff'/

      ndp8=ndp3/nir
      if (getiprint() >= 5) then
        call header('Output from "reord6"',-1)
        write(LUPRI,*) 'namef=',namef
        write(LUPRI,*) 'namet=',namet
        write(LUPRI,*) 'iord=',iord
        write(LUPRI,*) 'ndp8=ndp3/nir=',ndp8
        write(LUPRI,*) 'ndp4=',ndp4
      endif

      if (ndp8.gt.ndp4) then 
       call quit('ndp8 > ndp4 !')
      endif

      call files_ccsd (namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'reord6 input tape does not exist, ',name2
!       stop                      
        call quit('reord6 input tape does not exist')
      end if
      call parse6 (iord,idum1,idum2,i1,i2,i3,i4)
      if ((idum1.ne.1).or.(idum2.ne.2)) then
        write (6,*) 'reord6 called with incorrect iord'
!       stop
        call quit('reord6 called with incorrect iord')
      end if
      call parse6 (iii2,ia(1),ia(2),ia(3),ia(4),ia(5),ia(6))
      call parse6 (iiv2,ja(1),ja(2),ja(3),ja(4),ja(5),ja(6))
      iii1=100000*ia(1)+10000*ia(2)+1000*ia(i1)+                        &
      100*ia(i2)+10*ia(i3)+ia(i4)
      iiv1=100000*ja(1)+10000*ja(2)+1000*ja(i1)+                        &
      100*ja(i2)+10*ja(i3)+ja(i4)
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
          !stop
          call quit('mismatch in block lengths - reord6')
        end if
        ii=iand(ishft(ij2(i),-8),z'ff')
        jj=iand(ij2(i),z'ff')
        ii=igmt(ici(ii),ici(jj))
        if (igot(ii).eq.0) then
          igot(ii)=1
          istrt=istrtf1+ije2(i-1)+1
          call vands  (lngth,igen(istrt),msk(i1),itmp(1))
          call vshfcs (lngth,itmp(1),8*i1-8,itmp(1))
!         miro: ndp8 must be lower than ndp4, checked above
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
!miro: watch possible index leak ! (make check)
           itmp(ndp3+(ii-1)*ndp8+j)=int(tmp(j))
70        continue
        end if
        do 80 j=1,lngth
! miro: might get index leak for the "itmp" array !
          itmp_indx=ndp3+((ii-1)*ndp8)+j
          if (itmp_indx.gt.ndp4) then
            write(6,*) 'reord6:ndp4=',ndp4
            write(6,*) 'reord6:ndp3,ii,ndp8,j',ndp3,ii,ndp8,j
            write(6,*) 'reord6:itmp_indx=ndp3+((ii-1)*ndp8)+j=',        &
       itmp_indx
            call quit('reord6: itmp ugly index leak! Increase ndp4')
          endif
          gen(istrtt2+ije1(i-1)+j)=gen(istrtf2+ije2(i-1)+itmp(itmp_indx))
!    1    gen(istrtf2+ije2(i-1)+itmp(ndp3+(ii-1)*ndp8+j))               &
80      continue
10    continue

      !call check_modules

      call prtime ("reord6   ",0,1)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine reord7 (namef,namet,iord)   
!---------------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cfiles_mod 
      use gen_mod 
      use cij1_mod 
      use cij2_mod 
      use ctmp_mod 
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
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#endif

      dimension ia(6),ja(6),msk(3)
      data msk / z'ff0000',z'ff00',z'ff' /

      !call check_modules

      call files_ccsd (namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (6,*) 'reord7 input tape does not exist, ',name2
        call quit('reord7: input tape does not exist,')   
      end if
      call parse6 (iord,i1,i2,i3,idum1,idum2,idum3)
      call parse6 (iii2,ia(1),ia(2),ia(3),ia(4),ia(5),ia(6))
      call parse6 (iiv2,ja(1),ja(2),ja(3),ja(4),ja(5),ja(6))
      iii1=100000*ia(i1)+10000*ia(i2)+1000*ia(i3)+                      &
      100*ia(4)+10*ia(5)+ia(6)
      iiv1=100000*ja(i1)+10000*ja(i2)+1000*ja(i3)+                      &
      100*ja(4)+10*ja(5)+ja(6)
      call tmplt6 (iii1,iiv1,namet,ntapet)
      if (nij2.ne.nij1) then
        write (LUPRI,*)                                                 &
   'mismatch in numbers of blocks - reord7: nij2,nij1',nij2,nij1
        call quit('mismatch in numbers of blocks - reord7')
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
        call dcopy (lngth,gen(istrt22+ije2(itmp(i)-1)+1),1,             &
   gen(istrt12+ije1(i-1)+1),1)
30    continue

      !call check_modules

      call prtime ("reord7  ",0,1)

      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine reord10 (namef,namet,iord)
!---------------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cfiles_mod 
      use gen_mod 
      use cij1_mod 
      use cij2_mod 
      use sym_mod 
      use ijtemp_mod 
      use ctmp_mod 
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
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
#include "ijtemp.h"
!      common /ijtemp/ ijtmp1(ndp6)
#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#endif

      dimension ia(10),ja(10),msk(5),iord(10)
      data msk / z'f0000',z'f000',z'f00',z'f0',z'f' /

      !call check_modules

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

      !call check_modules

      do 10 i=1,nij2
        lngth=ije1(i)-ije1(i-1)
        call dimx (lngth,ndp5,"reord10 ")
        if (ije2(i)-ije2(i-1).ne.lngth) then
          write (LUPRI,*) 'mismatch in block lengths - reord10'
          call quit('mismatch in block lengths - reord10')
        end if
        do 50 j=1,lngth
          !miro: index leak may occur here
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

      !call check_modules

      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine reord11 (namef,namet,iord)         
!-----------------------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cfiles_mod 
      use gen_mod 
      use cij1_mod 
      use cij2_mod 
      use ctmp_mod 
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
!#include "param_ccsd.h"
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
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#endif

!!!!!! miro: ndij very big, needed dynamic allocs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     dimension ia(10),ja(10),msk(5),indx(ndij),dummy(ndij),iord(10)
      dimension indx(ndij),dummy(ndij)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension ia(10),ja(10),msk(5),iord(10)

      data msk / z'f0000',z'f000',z'f00',z'f0',z'f' /

      !call check_modules

      i1=iord(1)
      i2=iord(2)
      i3=iord(3)
      i4=iord(4)
      i5=iord(5)
      call tmplt10 (namet,ntapet)
      call files_ccsd (namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (LUPRI,*) 'reorder input tape does not exist, ',name2
        call quit('reorder input tape does not exist')             
      end if
      if (nij2.ne.nij1) then
        write (LUPRI,*) 'mismatch in numbers of blocks - reord11'
        call quit('mismatch in numbers of blocks - reord11')
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
          write (LUPRI,*) 'mismatch in block lengths - reord11'
          write (LUPRI,*) 'n1,n2:',n1,n2
          write (LUPRI,*) 'ije1(i),ije1(i-1):',ije1(i),ije1(i-1)
          write (LUPRI,*) 'ije2(indx(i)),ije2(indx(i)-1):',             &
                   ije2(indx(i)),ije2(indx(i)-1)
          call quit('mismatch in block lengths - reord11')
        end if
        do 20 j=ije2(indx(i)-1)+1,ije2(indx(i))
          tmp(k)=gen(istrtf2+j)
          k=k+1
20      continue
10    continue
      call dcopy (k,tmp,1,gen(istrtt2+1),1)   

      !call check_modules

      call prtime ("reord11 ",0,2)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine reord12 (namef,namet,iord)
!-----------------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cfiles_mod 
      use gen_mod 
      use cij1_mod 
      use cij2_mod 
      use sym_mod 
      use ijtemp_mod 
      use ctmp_mod 
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
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
#include "ijtemp.h"
!     common /ijtemp/ ijtmp1(ndp6)
#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#endif

      dimension ia(12),ja(12),msk(6),iord(12)
      data msk / z'f00000',z'f0000',z'f000',z'f00',z'f0',z'f' /

      !call check_modules

      call files_ccsd (namef,ntapef,new,2,-1,-1)

      if (new.ne.0) then
        write (LUPRI,*) 'reord12 input tape does not exist, ',name2
        call quit('reord12 input tape does not exist')
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
          call quit('mismatch in block lengths - reord12')
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

      !call check_modules

      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine reord13 (namef,namet,iord)         
!---------------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use cfiles_mod 
      use gen_mod 
      use cij1_mod 
      use cij2_mod 
      use ctmp_mod 
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
!#include "param_ccsd.h"
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
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#endif

!!!!!!!!! to be replaced by dynamic allocs !!!!!!!!!!!!!!!!!!!!!!!!
!     dimension ia(12),ja(12),msk(6),indx(ndij),dummy(ndij),iord(12)
      dimension indx(ndij),dummy(ndij)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension ia(12),ja(12),msk(6),iord(12)

      data msk / z'f00000',z'f0000',z'f000',z'f00',z'f0',z'f' /

      !call check_modules

      i1=iord(1)
      i2=iord(2)
      i3=iord(3)
      i4=iord(4)
      i5=iord(5)
      i6=iord(6)
      call tmplt12 (namet,ntapet)
      call files_ccsd (namef,ntapef,new,2,-1,-1)
      if (new.ne.0) then
        write (LUPRI,*) 'reorder input tape does not exist, ',name2
        call quit('reorder input tape does not exist')  
      end if
      if (nij2.ne.nij1) then
        write (LUPRI,*) 'mismatch in numbers of blocks - reord13'
        call quit('mismatch in numbers of blocks - reord13')       
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
          write (LUPRI,*) 'mismatch in block lengths - reord13'
          call quit('mismatch in block lengths - reord13') 
        end if
        do 20 j=ije2(indx(i)-1)+1,ije2(indx(i))
          tmp(k)=gen(istrtf2+j)
          k=k+1
20      continue
10    continue
      call dcopy (k,tmp,1,gen(istrtt2+1),1)   

      !call check_modules

      call prtime ("reord13 ",0,2)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine seltyp(iii,iiv,name,iparm)
!------------------------------------------------------------------------
!
!  integrals are collected in "gen", indices are in "itmp"
!
!  at the end, arrays are sorted and expanded according to indices (in "igen") generated by tmplt.
!
!  routine  is called many times from high.sector executables
!
!  it assumes vint,iint arrays filled with integrals/indexes !
!
!  parameters: iparm - printing level for "prt"
!
!  Important: integrals must be unique in Dirac notation ! No repeating
!
!------------------------------------------------------------------------
      use print_level, only:getiprint
      use ipack_mod, only: ipack4
#if defined HSFS_DYN_ALLOC
      use ctmp_mod 
      use cij1_mod 
      use cfiles_mod 
      use gen_mod 
      use dat_mod 
      use rsp2_mod 
      use cint_mod 
      use hsfs_parameters
implicit none
#include "priunit.h"
#else
!      implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
!     common /cij1/ nij, ij(ndij), ije(0:ndij), iiii,iivi
!miro: replaced by proper /cij1/ block
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
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
#endif

!    ... miro: introduced to have strings clean, without garbage
      character*4, intent(in) :: name
      integer, intent(in):: iparm
      integer, intent(in) :: iii, iiv

      logical :: lrsp
      !data zero / 1.d-13/ ! miro: to generalize in a module
      real*8, parameter :: zero=1.0d-13

      integer :: num,lngth,istrt1,istrt2,i1,i2,i3,i4,iv1,iv2,iv3,iv4
      integer :: nh,nv,i1234,nout,i,j1,j2,j3,j4,njh,njv,j12,ii,iv1234

      ! call template function
      call tmplt(iii,iiv,name,num)

!        ... assumes ije1(:) is filled
!     lngth=ije(nij)
      lngth=ije1(nij1)
!     istrt1=nstrt(1,num-1)+2*nij+4
      istrt1=nstrt(1,num-1)+2*nij1+4
      istrt2=nstrt(2,num-1)

      if (getiprint() >= 5) then ! control printout
        write(LUPRI,'(/,2x,a,i5)') 'seltyp, lngth=',lngth
        write(LUPRI,'(2x,a,i5)')   '       istrt1=',istrt1
        write(LUPRI,'(2x,a,i5)')   '       istrt2=',istrt2
      endif

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
      i1234=ipack4(i1,i2,i3,i4)
      iv1234=ipack4(iv1,iv2,iv3,iv4)

      nout=0
      do 190 i=1,nint
      if (abs(vint(i)).lt.zero) go to 190

!      TODO: replace with module function !
      j1=iand(ishft(iint(i),-24),z'ff')
      j2=iand(ishft(iint(i),-16),z'ff')
      j3=iand(ishft(iint(i),-8) ,z'ff')
      j4=iand(      iint(i),     z'ff')

      if (getiprint() >= 5) then ! control printout
         write(LUPRI,'(1x,i5,a,i10,a,4i4,a,f18.12)')   &
         i,'-',iint(i),' -',j1,j2,j3,j4,' :',vint(i)
      endif

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
      if (ipack4(ic(j1),ic(j2),ic(j3),ic(j4)).ne.i1234) go to 110
      if (iand(ipack4(icv(j1),icv(j2),icv(j3),icv(j4)),iv1234).ne.iv1234) go to 110

      if (lrsp) then
        j12=ior(ishft(j1,8),j2)
        do 105 ii=1,nrsp2
105     if (j12.eq.irsp2(ii)) go to 110
      endif

      nout=nout+1
      gen(istrt2+nout)=vint(i)
      ! <j1 j2|j3 j4> = <j1 j2|j3 j4>
      itmp(nout)=ipack4(j1,j2,j3,j4)

      if (getiprint() >= 5) then ! control printout
        write(LUPRI,*) nout,j1,j2,j3,j4,gen(istrt2+nout)
      endif

110   if (j1.eq.j3) go to 120

      if (ipack4(ic(j3),ic(j2),ic(j1),ic(j4)).ne.i1234) go to 120
      if (iand(ipack4(icv(j3),icv(j2),icv(j1),icv(j4)),iv1234).ne.iv1234) go to 120

      if (lrsp) then
        j12=ior(ishft(j3,8),j2)
        do 115 ii=1,nrsp2
115     if (j12.eq.irsp2(ii)) go to 120
      endif

      nout=nout+1
      gen(istrt2+nout)=vint(i)
      !<j1 j2|j3 j4>=<j3 j2|j1 j4>
      itmp(nout)=ipack4(j3,j2,j1,j4)

      if (getiprint() >= 5) then ! control printout
        write(LUPRI,*) nout,j3,j2,j1,j4,gen(istrt2+nout)
      endif


120   if (j2.eq.j4) go to 130
      if (ipack4(ic(j1),ic(j4),ic(j3),ic(j2)).ne.i1234) go to 130

      if (iand(ipack4(icv(j1),icv(j4),icv(j3),icv(j2)),iv1234).ne.iv1234) go to 130

      if (lrsp) then
        j12=ior(ishft(j1,8),j4)
        do 125 ii=1,nrsp2
125     if (j12.eq.irsp2(ii)) go to 130
      endif

      nout=nout+1
      gen(istrt2+nout)=vint(i)
      ! <j1 j2 |j3 j4> = <j1 j4|j3 j2>
      itmp(nout)=ipack4(j1,j4,j3,j2)

      if (getiprint() >= 5) then ! control printout
        write(LUPRI,*) nout,'-',j1,j4,j3,j2,gen(istrt2+nout)
      endif

130   if ((j1.eq.j3).or.(j2.eq.j4)) go to 140
      if (ipack4(ic(j3),ic(j4),ic(j1),ic(j2)).ne.i1234) go to 140

      if (iand(ipack4(icv(j3),icv(j4),icv(j1),icv(j2)),iv1234).ne.iv1234) go to 140

      if (lrsp) then
        j12=ior(ishft(j3,8),j4)
        do 135 ii=1,nrsp2
135     if (j12.eq.irsp2(ii)) go to 140
      endif

      nout=nout+1
      gen(istrt2+nout)=vint(i)
      ! <j1 j2 |j3 j4> = <j3 j4|j1 j2>
      itmp(nout)=ipack4(j3,j4,j1,j2)

      if (getiprint() >= 5) then ! control printout
        write(LUPRI,*) nout,j3,j4,j1,j2,gen(istrt2+nout)
      endif

140   if ((j1.eq.j2).and.(j3.eq.j4)) go to 190

      if (ipack4(ic(j2),ic(j1),ic(j4),ic(j3)).ne.i1234) go to 150

      if (iand(ipack4(icv(j2),icv(j1),icv(j4),icv(j3)),iv1234).ne.iv1234) go to 150

      if (lrsp) then
       j12=ior(ishft(j2,8),j1)
        do 145 ii=1,nrsp2
145     if (j12.eq.irsp2(ii)) go to 150
      endif

      nout=nout+1
      gen(istrt2+nout)=vint(i)
      ! <j1 j2 |j3 j4> = <j2 j1|j4 j3>
      itmp(nout)=ipack4(j2,j1,j4,j3)

      if (getiprint() >= 5) then ! control printout
        write(LUPRI,*) nout,'-',j2,j1,j4,j3,gen(istrt2+nout)
      endif

150   if (j1.eq.j3) go to 160
      if (ipack4(ic(j2),ic(j3),ic(j4),ic(j1)).ne.i1234) go to 160
      if (iand(ipack4(icv(j2),icv(j3),icv(j4),icv(j1)),iv1234).ne.iv1234) go to 160

      if (lrsp) then
        j12=ior(ishft(j2,8),j3)
        do 155 ii=1,nrsp2
155     if (j12.eq.irsp2(ii)) go to 160
      endif

      nout=nout+1
      gen(istrt2+nout)=vint(i)
      ! <j1 j2 |j3 j4> = <j2 j3|j4 j1>
      itmp(nout)=ipack4(j2,j3,j4,j1)

      if (getiprint() >= 5) then ! control printout
        write(LUPRI,*) nout,'-',j2,j3,j4,j1,gen(istrt2+nout)
      endif

160   if (j2.eq.j4) go to 170
      if (ipack4(ic(j4),ic(j1),ic(j2),ic(j3)).ne.i1234) go to 170
      if (iand(ipack4(icv(j4),icv(j1),icv(j2),icv(j3)),iv1234).ne.iv1234)  go to 170
      if (lrsp) then
        j12=ior(ishft(j4,8),j1)
        do 165 ii=1,nrsp2
165     if (j12.eq.irsp2(ii)) go to 170
      endif

      nout=nout+1
      gen(istrt2+nout)=vint(i)
      ! <j1 j2 |j3 j4> = <j4 j1|j2 j3>
      itmp(nout)=ipack4(j4,j1,j2,j3)

      if (getiprint() >= 5) then ! control printout
        write(LUPRI,*) nout,'-',j4,j1,j2,j3,gen(istrt2+nout)
        call check_modules
      endif

170   if ((j1.eq.j3).or.(j2.eq.j4)) go to 190
      if (ipack4(ic(j4),ic(j3),ic(j2),ic(j1)).ne.i1234) go to 190
      if (iand(ipack4(icv(j4),icv(j3),icv(j2),icv(j1)),iv1234).ne.iv1234) go to 190
      if (lrsp) then
        j12=ior(ishft(j4,8),j3)
        do 175 ii=1,nrsp2
175     if (j12.eq.irsp2(ii)) go to 190
      endif

      nout=nout+1
      gen(istrt2+nout)=vint(i)
      ! <j1 j2 |j3 j4> = <j4 j3|j2 j1>
      itmp(nout)=ipack4(j4,j3,j2,j1)

      if (getiprint() >= 5) then ! control printout
        write(LUPRI,*) nout,'-',j4,j3,j2,j1,gen(istrt2+nout)
      endif

190   continue ! end of cycle

!  control printout before calling vxort
      if (getiprint() >= 7) then ! control printout
       write(LUPRI,'(/,2x,a,i6)')   &
            'seltyp-before calling vxsort, nout=',nout
       write(LUPRI,'(2x,a,i6)') 'ije1(nij1)=',ije1(nij1)
       write(LUPRI,'(2x,a,i6)') 'istrt2+1=',istrt2+1
       write(LUPRI,'(2x,a,i6)') 'istrt1+1=',istrt1+1
      endif

!     call vxsort(itmp,nout,ije(nij),gen(istrt2+1),
      call vxsort(itmp,nout,ije1(nij1),gen(istrt2+1),igen(istrt1+1))

!     write (LUPRI,1900) ije(nij),iii,iiv
      write (LUPRI,1900) ije1(nij1),iii,iiv
1900  format(1x,'seltyp collected ',i9,' integrals of type ',2i6)
      call prt(name,'seltyp',iparm) ! ... print integrals ...

!     call prtime(name,ije(nij),0) ! ... print time ...
      call prtime(name,ije1(nij1),0) ! ... print time ...

      return
      end subroutine seltyp

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine selt1(iii,iiv,name,iparm)
!------------------------------------------------------------------------------
!
! Extracts and expands Fock matrix elements into gen,igen.
!
!------------------------------------------------------------------------------
      use ipack_mod, only: ipack2
#if defined HSFS_DYN_ALLOC
      use ctmp_mod 
      use cfiles_mod 
      use gen_mod 
      use cij1_mod 
      use dat_mod 
      use fock_mod 
      use rsp1_mod 
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)

!     common /cij1/ nij,ij(ndij),ije(0:ndij),iiii,iivi
!to be replaced with proper /cij1/ common block
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
#include "dat.h"
!     common /dat/   norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "fock.h"
!     common /fock/ nfock,idum3,fock(ndt1),ifock(ndt1)
#include "rsp1.h"
!     common /rsp1/  nrsp1,irsp1(100)
#endif

      character(*) :: name
      logical lrsp
      data zero / 1.d-13/

!      32-bit function (to be moved into module)
      !ipack(i1,i2)=ior(ishft(i1,8),i2)

        call check_modules

      if (iopt(8).eq.0) then
        call clear(name,iii,iiv)
        return
      endif

      call tmplt(iii,iiv,name,num)

!     lngth=ije(nij)
      lngth=ije1(nij1)

!     istrt1=nstrt(1,num-1)+2*nij+4
      istrt1=nstrt(1,num-1)+2*nij1+4

      istrt2=nstrt(2,num-1)
      call parse(iii,id,id,i1,i2)
      call parse(iiv,id,id,iv1,iv2)
      lrsp=(iopt(9).ne.0).and.(iv1.eq.1)
      i12=ipack2(i1,i2)
      iv12=ipack2(iv1,iv2)
      nout=0
      do 190 i=1,nfock
      if (abs(fock(i)).lt.zero) go to 190
      j1=iand(ishft(ifock(i),-8),z'ff')
      j2=iand(ifock(i),z'ff')
      if(ipack2(ic(j1),ic(j2)).ne.i12) go to 110
      if(iand(ipack2(icv(j1),icv(j2)),iv12).ne.iv12) go to 110
      if (lrsp) then
        j12=ior(ishft(j1,8),j2)
        do 105 ii=1,nrsp1
105     if (j12.eq.irsp1(ii)) go to 110
      endif
      nout=nout+1
      gen(istrt2+nout)=fock(i)
      itmp(nout)=ipack2(j1,j2)
110   if (j1.eq.j2) go to 190
      if(ipack2(ic(j2),ic(j1)).ne.i12) go to 190
      if(iand(ipack2(icv(j2),icv(j1)),iv12).ne.iv12) go to 190
      if (lrsp) then
        j12=ior(ishft(j2,8),j1)
        do 115 ii=1,nrsp1
115     if (j12.eq.irsp1(ii)) go to 190
      endif
      nout=nout+1
      gen(istrt2+nout)=fock(i)
      itmp(nout)=ipack2(j2,j1)
190   continue ! end of cycle i=1,nfock

!     call vxsort(itmp,nout,ije(nij),gen(istrt2+1),
      call vxsort(itmp,nout,ije1(nij1),gen(istrt2+1),igen(istrt1+1))

!     write (6,1900) ije(nij),iii,iiv
      write (6,1900) ije1(nij1),iii,iiv
1900  format('selt1 collected ',i9,' integrals of type ',2i6)

      call check_modules

      call prt(name,'selt1 ',iparm)
!     call prtime(name,ije(nij),0)
      call prtime(name,ije1(nij1),0)

      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine seltyp4(iii,iiv,name,iparm)
!---------------------------------------------------------------------------------
!
!  this routine for pppp only. integrals <ij|kl> collected, i.ge.j .
!
!  integrals are collected in "gen", indices in "itmp" 
!   at the end, arrays are
!  sorted and expanded according to indices (in igen) generated
!  in routine tmplt.
!
!----------------------------------------------------------------------------------
      use ipack_mod, only: ipack4
#if defined HSFS_DYN_ALLOC
      use ctmp_mod 
      use cij1_mod 
      use cfiles_mod 
      use gen_mod 
      use dat_mod 
      use cint_mod 
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
!     common /cij1/ nij,ij(ndij),ije(0:ndij),iiii,iivi
!miro: replaced by proper /cij1/ block
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
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
#endif

!    miro: to prevent ugly crash in prtime routine due to rubbish in
!    "name"
      character(*) name

      data zero / 1.e-13/

      !ipack(i1,i2,i3,i4)=ior(ior(ishft(i1,24),ishft(i2,16)),ior(ishft(i3,8),i4)           )

      call check_modules

      if ((iii.ne.3333).or.(iiv.ne.0)) then
        write (*,*) ' seltyp4 called with incorrect iii, iiv:',iii,iiv
        call abext('seltyp4 ')
      endif
      ieof=0
      call rdfile(9,0,"pppp",ieof)
      if (ieof.eq.0) return
      call tmplt(-iii,iiv,name,num)
!     lngth=ije(nij)
      lngth=ije1(nij1)
!     istrt1=nstrt(1,num-1)+2*nij+4
      istrt1=nstrt(1,num-1)+2*nij1+4
      istrt2=nstrt(2,num-1)
      nout=0
      do 190 i=1,nint
      if (abs(vint(i)).lt.zero) go to 190
!     call unpack(iint(i),j1,j2,j3,j4)
      j1=iand(ishft(iint(i),-24),z'ff')
      j2=iand(ishft(iint(i),-16),z'ff')
      j3=iand(ishft(iint(i),-8) ,z'ff')
      j4=iand(           iint(i),z'ff')
      if (ic(j1).ne.3) go to 190
      if (ic(j2).ne.3) go to 190
      if (ic(j3).ne.3) go to 190
      if (ic(j4).ne.3) go to 190
      if (j1.lt.j2) go to 110
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j1.eq.j2) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack4(j1,j2,j3,j4)
110   if (j1.eq.j3) go to 120
      if (j3.lt.j2) go to 120
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j2.eq.j3) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack4(j3,j2,j1,j4)
120   if (j2.eq.j4) go to 130
      if (j1.lt.j4) go to 130
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j1.eq.j4) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack4(j1,j4,j3,j2)
130   if ((j1.eq.j3).or.(j2.eq.j4)) go to 140
      if (j3.lt.j4) go to 140
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j3.eq.j4) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack4(j3,j4,j1,j2)
140   if ((j1.eq.j2).and.(j3.eq.j4)) go to 190
      if (j2.lt.j1) go to 150
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j1.eq.j2) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack4(j2,j1,j4,j3)
150   if (j1.eq.j3) go to 160
      if (j2.lt.j3) go to 160
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j3.eq.j2) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack4(j2,j3,j4,j1)
160   if (j2.eq.j4) go to 170
      if (j4.lt.j1) go to 170
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j1.eq.j4) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack4(j4,j1,j2,j3)
170   if ((j1.eq.j3).or.(j2.eq.j4)) go to 190
      if (j4.lt.j3) go to 190
         nout=nout+1
         gen(istrt2+nout)=vint(i)
         if (j3.eq.j4) gen(istrt2+nout)=vint(i)*0.5d0
         itmp(nout)=ipack4(j4,j3,j2,j1)
190   continue

      call vxsort(itmp,nout,lngth,gen(istrt2+1),igen(istrt1+1))

!  indices not retained for pppp.
      nstrt(1,nflcur)=nstrt(1,nflcur)-lngth
      write (6,1900) lngth,iii,iiv
1900  format('seltyp4 collected ',i9,' integrals of type ',2i6)
      if (iopt(10).eq.0) call wrtfile(9,0,"pppp")

      call prt(name,'seltyp4 ',iparm)
      call prtime(name,lngth,0)

 return
 end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine sqr3 (name) 
!------------------------------------------------------------------------
!
!  makes final symmetrization of the "twf " intermediate, what can be
!  wave-function or effective hamiltonian
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use prm_mod 
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "prm.h"
!     common /prm/ iperm(0:6,720)
#endif

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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine sqr4 (name)       
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use prm_mod 
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "prm.h"
!     common /prm/ iperm(0:6,720)
#endif
     
      call perms (4)
      call move  (name,"s1  ")
      do 10 i=2,24
        iord=10000000*iperm(1,i)+1000000*iperm(2,i)+                    &
               100000*iperm(3,i)+10000*iperm(4,i)+                      &
                 1000*(iperm(1,i)+4)+100*(iperm(2,i)+4)+                &
                    10*(iperm(3,i)+4)+iperm(4,i)+4
        call reordx  ("s1  ","s2  ",iord)
        call add     ("s2  ",1.d0,name,1.d0,name)
10    continue
      call prtime ("sqr4",0,2)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine sqr4l3 (name)
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use prm_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "prm.h"
!      common /prm/ iperm(0:6,720)
#endif

      call perms (3)
      call move  (name,"s1  ")
      do 10 i=2,6  
        iord=10000000+1000000*(iperm(1,i)+1)+                     &
        100000*(iperm(2,i)+1)+10000*(iperm(3,i)+1)+5000+          &
           100*(iperm(1,i)+5)+10*(iperm(2,i)+5)+iperm(3,i)+5
        call reordx  ("s1  ","s2  ",iord)
        call add     ("s2  ",1.d0,name,1.d0,name)
10    continue
      call prtime ("sqr4l3  ",0,2)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  subroutine sqr(namef,names)
!--------------------------------------------------------------------------------
!
! symmetrization; "namef" name of the nonsymmetrized diagram to be symmetrized 
! into "names"
!
!--------------------------------------------------------------------------------
implicit none
#include "priunit.h"
      character*4, intent(in) :: namef,names
      call reorder(namef,names,2143,5) ! switch indexes
      call add(namef,1.d0,names,1.d0,namef)
      return
 end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine tmplt(ii1234,j1234,name,num)
!--------------------------------------------------------------------------------
!
! What does it do ?
!
! Creates the structure of diagrams with iii,iiv numbers; make space for diagrams
! on gen,igen
!
! Calculates number of integrals, pairs, sorts integrals according to two first (ij) indexes 
!
! ii1234<0 - collect 1/2 file (i>j)
!
!--------------------------------------------------------------------------------
      use ipack_mod, only: ipack4
      use print_level
#if defined HSFS_DYN_ALLOC
      use cij1_mod
      use cfiles_mod
      use gen_mod
      use dat_mod
      use sym_mod
      use hsfs_parameters
implicit none
#include "priunit.h"
integer :: i,j,k,l,ii1234,j1234,name,num,i1,i1234,i2,i3,i4
integer :: i_4hs,ica,icb,icc,icd,ije_1,isaa,itmpx,j1,j2,j3,j4
integer :: ja,ja1,ja2,jb,jb1,jb2,jc,jc1,jc2,jd,jd1,jd2,jtmp
integer :: new,nt,nt0,nt00
!integer :: ipack
#else
!      implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
!     common /cij1/ nij,ij(ndij),ije(0:ndij),iii,iiv
! miro - replaced with proper common block
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
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
#endif

      !integer*4 ch1,ch2
      integer :: ch1,ch2

!     4-indexes (up 256) packed into one number
      !ipack(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

      if (getiprint() >= 5) then ! control printout
        call header('Output from "tmplt"',-1)
        !write(LUPRI,'(2x,a,i3)') 'nint=',nint
      endif

      i1234=iabs(ii1234)
      call parse(i1234,i1,i2,i3,i4)
      call parse(j1234,j1,j2,j3,j4)

      call files_ccsd(name,num,new,1,i1234,j1234)

      if (getiprint() >= 5) then ! control printout
        write(LUPRI,*) "tmplt: after files_ccsd, new=",new
        write(LUPRI,*) "tmplt: after files_ccsd, num=",num
      endif

      if (new.le.0) then
       ! write(LUPRI,*) "tmplt: new.le.0, return"
        if (getiprint() >= 8) then
          write(LUPRI,'(2x,a)') "tmplt: new <= 0, leaving"
        endif
        return
      endif

!     nij=0
      nij1=0
!     ije(0)=0
      ije1(0)=0
      nt=nstrt(1,num-1)+2*ndij+4
      nt0=nt

      if (getiprint() >= 5) then ! ... control printout
        write(LUPRI,'(/,2x,a,i12)') 'tmplt, nt=',nt
      endif

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
!  symmetry operations
      do 170 i=1,nsym
      if (mod(isym(ica,i)+isym(icb,i)+isym(icc,i)+isym(icd,i),2).ne.0) go to 240
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
! mi: complain due to leak :
! Array reference out of bounds for array 'igen', upper bound of dimension 1 exceeded (10098089 > 9000000)
! by changing parameters also ...upper bound of dimension 1 exceeded (1920005 > -1094967296)
      igen(nt)=ipack4(ja,jb,jc,jd)

      if (getiprint() >= 15) then
        write(LUPRI,*) 'tmplt: nt,igen(nt):',nt,igen(nt)
      endif

220   continue
230   continue
240   continue

      if (nt.gt.nt00) then
!       ije1=ije(nij)+1
        ije_1=ije1(nij1)+1
!       nij=nij+1
        nij1=nij1+1
!       ije(nij)=nt-nt0
        ije1(nij1)=nt-nt0
!       ij(nij)=ior(ishft(ja,8),jb)
        ij1(nij1)=ior(ishft(ja,8),jb)

        if (getiprint() >= 5) then
!         write(LUPRI,*) 'nij,ije(nij),ij(nij)=',nij,ije(nij),ij(nij)
          write(LUPRI,*) 'tmplt: nij1,ije1(nij1),ij1(nij1)=',           &
      nij1,ije1(nij1),ij1(nij1)
        endif

        call ivsort(igen(nt00+1),nt-nt00)
!	print '(Z8)',(igen(i),i=nt00+1,nt)
!       if (nij.gt.ndij) then 
        if (nij1.gt.ndij) then 
          write (LUPRI,*) 'tmplt: nij1 reached ',nij1,' at jab=',ja,jb
        endif
!       call dimx(nij,ndij,"tmpndij")
        call dimx(nij1,ndij,"tmpndij")
      endif

280   continue
290   continue
!
!     iii=i1234
      iii1=i1234
!     iiv=j1234
      iiv1=j1234
!     call dimx(nij,ndij,"tmpndij")
      call dimx(nij1,ndij,"tmpndij")
      go to 410
!  one electron
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
      igen(nt)=ipack4(0,0,jc,jd)
      if (getiprint() >= 15) then
        write(LUPRI,*) 'tmplt (2): nt,igen(nt): ',nt,igen(nt)
      endif
340   continue
350   continue
360   continue

!     nij=1
      nij1=1
!     ije(1)=nt-nt0
      ije1(1)=nt-nt0
!     ij(1)=0
      ij1(1)=0
      call ivsort(igen(nt0+1),nt-nt0)
!     iii=i1234
      iii1=i1234
!     iiv=j1234
      iiv1=j1234


!410   call kcopy(nij+1,nij,igen(nstrt(1,num-1)+1))
410    continue
       ! call kcopy(nij1+1,nij1,igen(nstrt(1,num-1)+1))
       ! miro: proper substitution due to replacement of common blocks by F90 allocatable arrays
       igen(nstrt(1,num-1)+1)=nij1
       do k=1,nij1
         igen(nstrt(1,num-1)+1+k)=ij1(k)
       enddo

!     call kcopy(nij+1,ije(0),igen(nstrt(1,num-1)+nij+2))
      call kcopy(nij1+1,ije1(0),igen(nstrt(1,num-1)+nij1+2))

!     igen(nstrt(1,num-1)+2*nij+3)=iii
      igen(nstrt(1,num-1)+2*nij1+3)=iii1
!     igen(nstrt(1,num-1)+2*nij+4)=iiv
      igen(nstrt(1,num-1)+2*nij1+4)=iiv1

      do 420 i=1,nt-nt0
!420   igen(nt0-2*(ndij-nij)+i)=igen(nt0+i)
420   igen(nt0-2*(ndij-nij1)+i)=igen(nt0+i)
!     nstrt(1,num)=nstrt(1,num-1)+(nt-nt0)+2*nij+4
      nstrt(1,num)=nstrt(1,num-1)+(nt-nt0)+2*nij1+4
      nstrt(2,num)=nstrt(2,num-1)+(nt-nt0)
      call dimx(nstrt(1,num),nbig,"tmplt1")
      call dimx(nstrt(2,num),nbig,"tmplt2")

!   special for scratch files, which start with "s".
!          itmp=4H       
           !itmpx=4H       
           itmpx=transfer("    ",itmpx)
           jtmp=z'ffffff00'
!          ch1 = iand(itmp,jtmp)
           ch1 = iand(itmpx,jtmp)
           ch2 = iand (name, z'000000ff')
           isaa = ior (ch1, ch2)

      ! fix of i_4hs=4hs
      i_4hs=transfer("s   ",i_4hs)
      if (isaa.eq.i_4hs) then
        call dimx(nt-nt0,ndp3,"tmplt s")
        nstrt(1,num)=nstrt(1,num-1)+ndp3+2*ndij+4
        nstrt(2,num)=nstrt(2,num-1)+ndp3
        call dimx(nstrt(1,num),nbig,"tmplt1")
        call dimx(nstrt(2,num),nbig,"tmplt2")
      endif

      if (getiprint() >= 10) then ! control printout
        write(LUPRI,'(4x,a)') 'tmplt - cijn control printout:'
        call print_cijx_stuff
      endif

      return
      end subroutine tmplt

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine tmplt6 (iii,iiv,name,num)
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
      use ipack_mod, only: ipack3,ipack4
#if defined HSFS_DYN_ALLOC
      use cij1_mod
      use cfiles_mod
      use gen_mod
      use dat_mod
      use sym_mod
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
!     common /cij1/ nij,ij(ndij),ije(0:ndij),iii1,iiv1
! miro - replaced with proper common block
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
#include "dat.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
#endif

      !integer*4 ch1,ch2
      integer :: ch1,ch2

      !ipack3(i,j,k)=ior(ior(ishft(i,16),ishft(j,8)),k)
      !ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

      call parse6 (iii,i1,i2,i3,i4,i5,i6)
      call parse6 (iiv,j1,j2,j3,j4,j5,j6)
      call files_ccsd (name,num,new,1,iii,iiv)
      if (new.le.0) return
!     nij=0
      nij1=0
!     ije(0)=0
      ije1(0)=0
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
              if (mod(isym(ica,ii)+isym(icb,ii)+                        &
         isym(icc,ii)+isym(icd,ii)+                                     &
         isym(ice,ii)+isym(icf,ii),2).ne.0) then
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
!            nij=nij+1
             nij1=nij1+1
!            ije(nij)=nt-nt0
             ije1(nij1)=nt-nt0
!            ij(nij)=ipack3(ja,jb,jc)
             ij1(nij1)=ipack3(ja,jb,jc)
             call ivsort (igen(nt00+1),nt-nt00)

!            if (nij.gt.ndij) then 
             if (nij1.gt.ndij) then 
!              write (6,*) 'nij reached ',
               write (6,*) 'nij1 reached ',nij1,' at jabc=',ja,jb,jc
!    1       nij,' at jabc=',ja,jb,jc                                   &
             endif

!            call dimx (nij,ndij,"tmpndij ")
             call dimx (nij1,ndij,"tmpndij ")
            end if
           end if
30        continue
         end if
20      continue
       end if
10    continue
!     call dimx (ije(nij),ndp3,"tmp6    ")
      call dimx (ije1(nij1),ndp3,"tmp6    ")
!     igen(nstrt(1,num-1)+1)=nij
      igen(nstrt(1,num-1)+1)=nij1
!     call kcopy (nij,ij(1),igen(nstrt(1,num-1)+2))
      call kcopy (nij1,ij1(1),igen(nstrt(1,num-1)+2))
!     call kcopy (nij+1,ije(0),igen(nstrt(1,num-1)+nij+2))
      call kcopy (nij1+1,ije1(0),igen(nstrt(1,num-1)+nij1+2))
!          iii, iiv entering
!     igen(nstrt(1,num-1)+2*nij+3)=iii
      igen(nstrt(1,num-1)+2*nij1+3)=iii
!     igen(nstrt(1,num-1)+2*nij+4)=iiv
      igen(nstrt(1,num-1)+2*nij1+4)=iiv
      do 90 i=1,nt-nt0
!       igen(nt0-2*(ndij-nij)+i)=igen(nt0+i)
        igen(nt0-2*(ndij-nij1)+i)=igen(nt0+i)
90    continue
!     nstrt(1,num)=nstrt(1,num-1)+2*nij+4+(nt-nt0)
      nstrt(1,num)=nstrt(1,num-1)+2*nij1+4+(nt-nt0)
      nstrt(2,num)=nstrt(2,num-1)+(nt-nt0)
      call dimx (nstrt(1,num),nbig,"tmplt61 ")
      call dimx (nstrt(2,num),nbig,"tmplt62 ")
!   special for scratch files, which start with s
           !itmp=4H       
           itmp=transfer("    ",itmp)
           jtmp=z'ffffff00'
           ch1 = iand(itmp,jtmp)
           ch2 = iand (name, z'000000ff')
           isaa = ior (ch1, ch2)

      !fix of i_4hs=4hs
      i_4hs=transfer("s   ",i_4hs)
      if (isaa.eq.i_4hs) then
        call dimx (nt-nt0,ndp3,"tmplt s ")
        nstrt(1,num)=nstrt(1,num-1)+ndp3+2*ndij+4
        nstrt(2,num)=nstrt(2,num-1)+ndp3
        call dimx (nstrt(1,num),nbig,"tmplt1   ")
        call dimx (nstrt(2,num),nbig,"tmplt2   ")
      end if
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine tmplt8 (iii,iiv,name,num)
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
      use ipack_mod, only: ipack4
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use dat_mod
      use sym_mod
      use dpd_mod
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
!#include "cij1_0.h"
!     common /cij1/ nij,ij(ndij),ije(0:ndij),iii1,iiv1
! miro - replaced with proper common block
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
#include "dat.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
#include "dpd.h"
!      common /dpd/ nir,igmt(8,8),ici(200)
#endif

      !integer*4 ch1,ch2
      integer :: ch1,ch2

!!! miro - into module !!!
      !ipack(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

      call parse8 (iii,i1,i2,i3,i4,i5,i6,i7,i8)
      call parse8 (iiv,j1,j2,j3,j4,j5,j6,j7,j8)
      call files_ccsd (name,num,new,1,iii,iiv)
      if (new.le.0) return
!     nij=0
      nij1=0
!     ije(0)=0
      ije1(0)=0
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
                 if ((ic(je).eq.i5).and.                                &
            (icv(je).ge.j5)) then
                  do 80 jf=jf1,jf2
                   if ((ic(jf).eq.i6).and.                              &
              (icv(jf).ge.j6)) then
                    do 90 jg=jg1,jg2
                     if ((ic(jg).eq.i7).and.                            &
                (icv(jg).ge.j7)) then
                      do 100 jh=jh1,jh2
                       if ((ic(jh).eq.i8).and.                          &
                  (icv(jh).ge.j8)) then
                        nt=nt+1
                        call dimx (nt,nbig,"tmp8")
                        igen(nt)=ipack4(je,jf,jg,jh)
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
!              nij=nij+1
               nij1=nij1+1
!              ije(nij)=nt-nt0
               ije1(nij1)=nt-nt0
!              ij(nij)=ipack(ja,jb,jc,jd)
               ij1(nij1)=ipack4(ja,jb,jc,jd)
               call ivsort (igen(nt00+1),nt-nt00)
!              if (nij.gt.ndij) write (6,*) 'nij reached ',
!    1         nij,' at jabcd=',ja,jb,jc,jd                             &
               if (nij1.gt.ndij) write (6,*) 'nij1 reached ',nij1,' at jabcd=',ja,jb,jc,jd
!              call dimx (nij,ndij,"tmpndij ")
               call dimx (nij1,ndij,"tmpndij ")
              end if
             end if
40          continue
           end if
30        continue
         end if
20      continue
       end if
10    continue

!     igen(nstrt(1,num-1)+1)=nij
      igen(nstrt(1,num-1)+1)=nij1

!     call kcopy (nij,ij(1),igen(nstrt(1,num-1)+2))
      call kcopy (nij1,ij1(1),igen(nstrt(1,num-1)+2))
!     call kcopy (nij+1,ije(0),igen(nstrt(1,num-1)+nij+2))
      call kcopy (nij1+1,ije1(0),igen(nstrt(1,num-1)+nij1+2))
!     igen(nstrt(1,num-1)+2*nij+3)=iii
      igen(nstrt(1,num-1)+2*nij1+3)=iii
!     igen(nstrt(1,num-1)+2*nij+4)=iiv
      igen(nstrt(1,num-1)+2*nij1+4)=iiv
      do 110 i=1,nt-nt0
!       igen(nt0-2*(ndij-nij)+i)=igen(nt0+i)
        igen(nt0-2*(ndij-nij1)+i)=igen(nt0+i)
110   continue
!     nstrt(1,num)=nstrt(1,num-1)+2*nij+4+(nt-nt0)
      nstrt(1,num)=nstrt(1,num-1)+2*nij1+4+(nt-nt0)
      nstrt(2,num)=nstrt(2,num-1)+(nt-nt0)
      call dimx (nstrt(1,num),nbig,"tmplt61 ")
      call dimx (nstrt(2,num),nbig,"tmplt62 ")
!   special for scratch files, which start with s
           itmp=transfer("    ",itmp) ! fix of itmp=4H
           jtmp=z'ffffff00'
           ch1 = iand(itmp,jtmp)
           ch2 = iand (name, z'000000ff')
           isaa = ior (ch1, ch2)

      !fix of i_4hs=4hs
      i_4hs=transfer("s   ",i_4hs)
      if (isaa.eq.i_4hs) then

        call dimx (nt-nt0,ndp3,"tmplt s ")
        nstrt(1,num)=nstrt(1,num-1)+ndp3+2*ndij+4
        nstrt(2,num)=nstrt(2,num-1)+ndp3
        call dimx (nstrt(1,num),nbig,"tmplt1   ")
        call dimx (nstrt(2,num),nbig,"tmplt2   ")
      end if
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine tmplt10 (name,num)
!------------------------------------------------------------------------
!
! Called from: diage5
!
!------------------------------------------------------------------------
      use ipack_mod, only: ipack5
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use dat_mod
      use sym_mod
      use idv_mod
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
!     common /cij1/ nij,ij(ndij),ije(0:ndij),iii1,iiv1
! miro - replaced with proper common block
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
#include "dat.h"
!      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
#include "idv.h"
!     common /idv/ iv(200),nvalo
#endif

      !integer*4 ch1,ch2
      integer :: ch1,ch2

!!!!!!!!!!!!!!!! assign own dynamic allocation !!!!!!!!!!!!!!!!!!
      dimension icsym(200)

     ! ipack(i,j,k,l,m)=ior(ishft(i,16),ior(ior(ishft(j,12),ishft(k,8)),ior(ishft(l,4),m)))

      iii=5
      iiv=1
      call files_ccsd (name,num,new,1,iii,iiv)
      if (new.le.0) return
!     nij=0
      nij1=0
!     ije(0)=0
      ije1(0)=0
      nt=nstrt(1,num-1)+2*ndij+4
      nt0=nt
      isum=0
      ! make symmetry table
      do 10 i=1,nclass
        do 20 j=1,iclass(i)
          icsym(isum+j)=i
20      continue
        isum=isum+iclass(i)
10    continue

      call check_modules

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

!miro: fixing first index leaks (shows when range-check compiler flag included)
        if (icsym(iv(jj)).ge.1.and.icsym(iv(jj)).le.16.and.             &
       icsym(iv(jh)).ge.1.and.icsym(iv(jh)).le.16.and.                  &
       icsym(iv(jf)).ge.1.and.icsym(iv(jf)).le.16.and.                  &
       icsym(iv(jd)).ge.1.and.icsym(iv(jd)).le.16.and.                  &
       icsym(iv(ji)).ge.1.and.icsym(iv(ji)).le.16.and.                  &
       icsym(iv(jg)).ge.1.and.icsym(iv(jg)).le.16.and.                  &
       icsym(iv(je)).ge.1.and.icsym(iv(je)).le.16.and.                  &
       icsym(iv(jc)).ge.1.and.icsym(iv(jc)).le.16.and.                  &
       icsym(iv(ja)).ge.1.and.icsym(iv(ja)).le.16.and.                  &
       icsym(iv(jb)).ge.1.and.icsym(iv(jb)).le.16) then

         if (mod(isym((icsym(iv(ja))),ii)+                              &
    isym((icsym(iv(jb))),ii)+isym((icsym(iv(jc))),ii)+                  &
    isym((icsym(iv(jd))),ii)+isym((icsym(iv(je))),ii)+                  &
    isym((icsym(iv(jf))),ii)+isym((icsym(iv(jg))),ii)+                  &
    isym((icsym(iv(jh))),ii)+isym((icsym(iv(ji))),ii)+                  &
    isym((icsym(iv(jj))),ii),2).ne.0) then

          iflag=0

         end if
         end if

100     continue

        if (iflag.eq.1) then
         nt=nt+1
         call dimx (nt,nbig,"tmp10    ")
         igen(nt)=ipack5(jf,jg,jh,ji,jj)
        end if

90     continue
       if (nt.gt.nt00) then
!       nij=nij+1
        nij1=nij1+1
!       ije(nij)=nt-nt0
        ije1(nij1)=nt-nt0
!       ij(nij)=ipack(ja,jb,jc,jd,je)          
        ij1(nij1)=ipack5(ja,jb,jc,jd,je)          
        call ivsort (igen(nt00+1),nt-nt00)
!       if (nij.gt.ndij) write (6,*) 'nij reached ',
        if (nij1.gt.ndij) write (6,*) 'nij1 reached ',nij1,' at jabcde=',ja,jb,jc,jd,je
!    1  nij,' at jabcde=',ja,jb,jc,jd,je                                &
!       call dimx (nij,ndij,"tmpndij ")
        call dimx (nij1,ndij,"tmpndij ")
       end if
30    continue

!     igen(nstrt(1,num-1)+1)=nij
      igen(nstrt(1,num-1)+1)=nij1
!     call kcopy (nij,ij(1),igen(nstrt(1,num-1)+2))
      call kcopy (nij1,ij1(1),igen(nstrt(1,num-1)+2))
!     call kcopy (nij+1,ije(0),igen(nstrt(1,num-1)+nij+2))
      call kcopy (nij1+1,ije1(0),igen(nstrt(1,num-1)+nij1+2))
!     igen(nstrt(1,num-1)+2*nij+3)=iii
      igen(nstrt(1,num-1)+2*nij1+3)=iii
!     igen(nstrt(1,num-1)+2*nij+4)=iiv
      igen(nstrt(1,num-1)+2*nij1+4)=iiv
      do 150 i=1,nt-nt0
!       igen(nt0-2*(ndij-nij)+i)=igen(nt0+i)
        igen(nt0-2*(ndij-nij1)+i)=igen(nt0+i)
150   continue
!     nstrt(1,num)=nstrt(1,num-1)+2*nij+4+(nt-nt0)
      nstrt(1,num)=nstrt(1,num-1)+2*nij1+4+(nt-nt0)
      nstrt(2,num)=nstrt(2,num-1)+(nt-nt0)
      call dimx (nstrt(1,num),nbig,"tmplt61 ")
      call dimx (nstrt(2,num),nbig,"tmplt62 ")
!   special for scratch files, which start with s
           ! fix of itmp=4H       
           itmp=transfer("    ",itmp)
           jtmp=z'ffffff00'
           ch1 = iand(itmp,jtmp)
           ch2 = iand (name, z'000000ff')
           isaa = ior (ch1, ch2)

      ! fix of "4hs   "
      i_4hs=transfer("s   ",i_4hs) 
      if (isaa.eq.i_4hs ) then
        call dimx (nt-nt0,ndp3,"tmplt s ")
        nstrt(1,num)=nstrt(1,num-1)+ndp3+2*ndij+4
        nstrt(2,num)=nstrt(2,num-1)+ndp3
        call dimx (nstrt(1,num),nbig,"tmplt1   ")
        call dimx (nstrt(2,num),nbig,"tmplt2   ")
      end if
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine tmplt12 (name,num)
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
      use ipack_mod, only: ipack6
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use dat_mod
      use sym_mod
      use idv_mod
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
!     common /cij1/ nij,ij(ndij),ije(0:ndij),iii1,iiv1
! miro - replaced with proper common block
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij1.h"
#include "dat.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#include "sym.h"
!     common /sym/ nsym,nclass,iclass(16),isym(16,8),ibias(17)
#include "idv.h"
!     common /idv/ iv(200),nvalo
#endif

      !integer*4 ch1,ch2
      integer :: ch1,ch2

!!!!!!!!!!!!!!!! miro: need replace !!!!!!!!!!!!!!!!!!!
      dimension icsym(200)

    ! ipack(i,j,k,l,m,n)=ior(ior(ishft(i,20),ishft(j,16)),  &
    !                    ior(ior(ishft(k,12),ishft(l,8)),ior(ishft(m,4),n)))

      iii=6
      iiv=1
      call files_ccsd (name,num,new,1,iii,iiv)
      if (new.le.0) return
!     nij=0
      nij1=0
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
         if (mod(isym((icsym(iv(ja))),ii)+                              &
    isym((icsym(iv(jb))),ii)+isym((icsym(iv(jc))),ii)+                  &
    isym((icsym(iv(jd))),ii)+isym((icsym(iv(je))),ii)+                  &
    isym((icsym(iv(jf))),ii)+isym((icsym(iv(jg))),ii)+                  &
    isym((icsym(iv(jh))),ii)+isym((icsym(iv(ji))),ii)+                  &
    isym((icsym(iv(jj))),ii)+isym((icsym(iv(jk))),ii)+                  &
    isym((icsym(iv(jl))),ii),2).ne.0) then
          iflag=0
         end if
100     continue
        if (iflag.eq.1) then
         nt=nt+1
         call dimx (nt,nbig,"tmp12   ")
         igen(nt)=ipack6(jg,jh,ji,jj,jk,jl)
        end if
90     continue

       if (nt.gt.nt00) then
!       nij=nij+1
        nij1=nij1+1
!       ije(nij)=nt-nt0
        ije1(nij1)=nt-nt0
!       ij(nij)=ipack(ja,jb,jc,jd,je,jf) 
        ij1(nij1)=ipack6(ja,jb,jc,jd,je,jf) 
        call ivsort (igen(nt00+1),nt-nt00)
!       if (nij.gt.ndij) write (6,*) 'nij reached ',
        if (nij1.gt.ndij) write (6,*) 'nij1 reached ',nij1,' at jabcdef=',ja,jb,jc,jd,je,jf
!    1  nij,' at jabcdef=',ja,jb,jc,jd,je,jf                            &
!       call dimx (nij,ndij,"tmpndij ")
        call dimx (nij1,ndij,"tmpndij ")
       end if

30    continue
!     igen(nstrt(1,num-1)+1)=nij
      igen(nstrt(1,num-1)+1)=nij1
!     call kcopy (nij,ij(1),igen(nstrt(1,num-1)+2))
      call kcopy (nij1,ij1(1),igen(nstrt(1,num-1)+2))
!     call kcopy (nij+1,ije(0),igen(nstrt(1,num-1)+nij+2))
      call kcopy (nij1+1,ije1(0),igen(nstrt(1,num-1)+nij1+2))
!     igen(nstrt(1,num-1)+2*nij+3)=iii
      igen(nstrt(1,num-1)+2*nij1+3)=iii
!     igen(nstrt(1,num-1)+2*nij+4)=iiv
      igen(nstrt(1,num-1)+2*nij1+4)=iiv
      do 150 i=1,nt-nt0
!       igen(nt0-2*(ndij-nij)+i)=igen(nt0+i)
        igen(nt0-2*(ndij-nij1)+i)=igen(nt0+i)
150   continue
!     nstrt(1,num)=nstrt(1,num-1)+2*nij+4+(nt-nt0)
      nstrt(1,num)=nstrt(1,num-1)+2*nij1+4+(nt-nt0)
      nstrt(2,num)=nstrt(2,num-1)+(nt-nt0)
      call dimx (nstrt(1,num),nbig,"tmplt61 ")
      call dimx (nstrt(2,num),nbig,"tmplt62 ")
!   special for scratch files, which start with s
           !itmp=4H       
           itmp=transfer("    ",itmp)
           jtmp=z'ffffff00'
           ch1 = iand(itmp,jtmp)
           ch2 = iand (name, z'000000ff')
           isaa = ior (ch1, ch2)
      !fix of i_4hs=4hs
      i_4hs=transfer("s   ",i_4hs)
      if (isaa.eq.i_4hs) then
        call dimx (nt-nt0,ndp3,"tmplt s ")
        nstrt(1,num)=nstrt(1,num-1)+ndp3+2*ndij+4
        nstrt(2,num)=nstrt(2,num-1)+ndp3
        call dimx (nstrt(1,num),nbig,"tmplt1   ")
        call dimx (nstrt(2,num),nbig,"tmplt2   ")
      end if
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine t1sqr(name1a,name1b,name2,namet)
!------------------------------------------------------------------------
!
! what does it do ?
!
!------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use ctmp_mod
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
#endif

      if (getiprint() >= 5) then
        call header('Output from t1sqr',-1)
        write(LUPRI,*) 'name1a=',name1a
        write(LUPRI,*) 'name1b=',name1b
        write(LUPRI,*) 'name2=',name2
        write(LUPRI,*) 'namet=',namet
        !call check_ctmp_mod
      endif

      call files_ccsd(name1a,ntape1a,new,2,-1,-1)
      if (new.ne.0) then
        write (LUPRI,*) ' t1a file in t1sqr does not exist'
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

      if (getiprint() >= 9) then
        write(LUPRI,'(/,4x,a)')  't1sqr data:'
        write(LUPRI,*) 'nij1,ije1(nij1),iii1,iiv1=',                    &
                   nij1,ije1(nij1),iii1,iiv1
        write(LUPRI,*) 'nij2,ije2(nij2),iii2,iiv2=',                    &
                   nij2,ije2(nij2),iii2,iiv2
        write(LUPRI,*) 'nij3,ije3(nij3),iii3,iiv3=',                    &
                   nij3,ije3(nij3),iii3,iiv3
        call check_modules
      endif
!
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
!   .......................................................
      do 120 i2=1,nsize2
      do 120 i3=1,nsize3
      ntmp=ntmp+1
      gen(istrt12+ntmp)=gen(istrt22+i2)*gen(istrt32+i3)
      it2=igen(istrt21+i2)
      it3=igen(istrt31+i3)
120   itmp(ntmp)=ior(ior(ishft(iand(ishft(it2,-8),z'ff'),24),           &
                    ishft(iand(ishft(it3,-8),z'ff'),16)),               &
                ior(ishft(iand(it2,z'ff')          , 8),                &
                    iand(it3,z'ff')                    ))

      if (getiprint() >= 8) then
        write(LUPRI,'(/,2x,a)') 't1sqr: before call vxsort'
        write(LUPRI,*) 'ntmp=',ntmp
        write(LUPRI,*) 'nsize1=',nsize1
        write(LUPRI,*) 'nsize2=',nsize2
        write(LUPRI,*) 'nsize3=',nsize3
        write(LUPRI,*) 'istrt11=',istrt11
        write(LUPRI,*) 'istrt12=',istrt12
        call check_modules
      endif

      call vxsort(itmp,ntmp,nsize1,gen(istrt12+1),igen(istrt11+1))

      if (getiprint() >= 8) then
        write(LUPRI,*) 't1sqr: after vxsort, before call add'
        write(LUPRI,*) 'namet=',namet 
        write(LUPRI,*) 'name2=',name2
        call check_modules
      endif

      call add(namet,1.d0,name2,1.d0,namet)

      if (getiprint() >= 8) then
        write(LUPRI,*) 't1sqr: after call add'
        write(LUPRI,*) 'namet=',namet 
        write(LUPRI,*) 'name2=',name2
        call check_modules
      endif

      return
      end subroutine t1sqr

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine t1t2 (name1,name2,name3,names,iord,iparm) 
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use ctmp_mod
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
!#include "param_ccsd.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
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
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij3.h"
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#endif

!!!!!!!!!!!!!!!!!! ninds needs dynamic alloc !!!!!!!!!!!!!!!!!!!!
      dimension ia(5,6),ind(12),ninds(ndij)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      data ia / 5,1,1,0,0,1,5,2,0,0,2,2,5,0,0,                          &
           6,3,3,0,0,3,6,4,0,0,4,4,6,0,0 /

      call files_ccsd (name3,ntape3,new,1,-1,-1)
      if (new.ne.0) then
        write (LUPRI,*) ' t3 file in t1t2 does not exist'
        call quit('t3 file in t1t2 does not exist')
      endif
      call tmplt6 (iii1,iiv1,names,ntapes)
      call files_ccsd (names,ntapes,new,1,-1,-1)
      istrt11=nstrt(1,ntapes-1)+2*nij1+4
      istrt12=nstrt(2,ntapes-1)

      call dcopy (ije1(nij1),0.d0,0,gen(istrt12+1),1) 

      ! miro - problematic call here ... SIGSEGV
      continue
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
              ind(k+6)=iand(ishft(ii,(k-4)*8),z'ff')
40          continue
            ii=igen(istrt31+j)
            ind(11)=iand(ishft(ii,-8),z'ff')
            ind(12)=iand(ii,z'ff')
            do 50 k=1,6
              ind(k)=ind(ia(5,k)+6)
50          continue
            idnew=ior(ishft(ind(1),8),ind(2))
            indnew=ior(ior(ishft(ind(3),24),ishft(ind(4),16)),          &
       ior(ishft(ind(5),8),ind(6)))
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
          call vxsort (itmp(ije1(i-1)+1),ninds(i),ije1(i)-ije1(i-1),    &
     gen(istrt12+ije1(i-1)+1),igen(istrt11+ije1(i-1)+1))
        end if
70    continue

      !  call check_modules

      call add (names,1.d0,name3,1.d0,names)
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine t1t3 (name1,name3,name4,namet,iparm) 
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
      use ipack_mod, only: ipack4
#if defined HSFS_DYN_ALLOC
      use ctmp_mod
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
!#include "param_ccsd.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
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
#endif

!!!!!!!!!!!!! needs dynamic alloc for ninds !!!!!!!!!!!!!!!
!     dimension ia(4,8),ind(16),ninds(ndij)
      dimension ninds(ndij)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension ia(4,8),ind(16)

      data ia / 1,3,3,3,3,1,4,4,4,4,1,5,5,5,5,1,   &
                2,6,6,6,6,2,7,7,7,7,2,8,8,8,8,2 /

!!!! TODO: into module !
      !ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

        !call check_ctmp_mod

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
        write (LUPRI,*) ' t3 file in t1t3 does not exist'
        call quit(' t3 file in t1t3 does not exist')
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
              ind(9)=iand(ishft(ii,-8),z'ff')
              ind(10)=iand(ii,z'ff')
              ii=ij3(j) 
              jj=igen(istrt31+k)
              do 50 m=1,3
                ind(m+10)=iand(ishft(ii,(m-3)*8),z'ff')
                ind(m+13)=iand(ishft(jj,(m-3)*8),z'ff')
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
                  write (LUPRI,*) 't1t3 - block id not found '
                  call quit('t1t3 - block id not found ')
                end if
              end do
            end if
40        continue
30      continue
20    continue
      do 80 i=1,nij1
        if (ninds(i).gt.0) then
          call vxsort (itmp(ije1(i-1)+1),ninds(i),ije1(i)-ije1(i-1),    &
                     gen(istrt12+ije1(i-1)+1),igen(istrt11+ije1(i-1)+1))
        end if
80    continue

      call add (namet,1.d0,name4,1.d0,namet)

      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine t2sqr (name2a,name2b,name4,namet,iparm) 
!------------------------------------------------------------------------
!
!
!
! add some description...
!
!------------------------------------------------------------------------
      use ipack_mod, only: ipack4
      use print_level
#if defined HSFS_DYN_ALLOC
      use ctmp_mod
      use cfiles_mod
      use gen_mod
      use cij1_mod
      use cij2_mod
      use cij3_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
!#include "param_ccsd.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
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
!      common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij3.h"
!      common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#endif

!!!!!!!!!!!!!!!!!!!!!!!!! need dynamic alloc !!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! dimension ia(3,8),ind(16),ninds(ndij)
      dimension ninds(ndij)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension ia(3,8),ind(16)

      data ia / 1,1,1,2,5,5,5,2,6,6,6,2,3,3,3,4,7,7,7,4,8,8,8,4 /

     ! ipack4(i,j,k,l)=ior(ior(ishft(i,24),ishft(j,16)),ior(ishft(k,8),l))

      if (getiprint() >= 5) then ! control printout
        call header('Output from "t2sqr"',-1)
        write(LUPRI,*) 'iparm=',iparm     
      endif
        !call check_ctmp_mod

      iparm=iparm-1
      call files_ccsd (name4,ntape4,new,1,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t4 file in t2sqr does not exist'
        call quit('t4 file in t2sqr does not exist')
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
        call quit(' t2a file in t2sqr does not exist')
      endif
      call files_ccsd (name2b,ntape2b,new,3,-1,-1)
      if (new.ne.0) then
        write (6,*) ' t2b file in t2sqr does not exist'
        call quit(' t2b file in t2sqr does not exist')
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
              ind(k+8)=iand(ishft(ii,(k-4)*8),z'ff')
              ind(k+12)=iand(ishft(jj,(k-4)*8),z'ff')
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
                call quit('t2sqr - block id not found')
              end if
            end do
          end if
30      continue
20    continue

      do 70 i=1,nij1
        if (ninds(i).gt.0) then
          call vxsort (itmp(ije1(i-1)+1),ninds(i),ije1(i)-ije1(i-1),    &
                  gen(istrt12+ije1(i-1)+1),igen(istrt11+ije1(i-1)+1))
        end if
70    continue

      call add (namet,1.d0,name4,1.d0,namet)

      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine tchng6 (nij,ij,ije,ntapef,iparm) 
!------------------------------------------------------------------------
! 
!   what does it do ? 
! 
!   entering arrays ij(*),ije(0:*)
!
!------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
! ... this is for ij(ndij),ije...
!#include "param_ccsd.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1iiia(200),iiva(200)
#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)
#endif
      !dimension ij(ndij),ije(0:ndij) ! entering !!!
      dimension ij(*),ije(0:*) ! entering !!!

!!!!!!!!! needs own dynamic allocation for ij1,ije1, itmp !!!!!!!!!!!!!!!!!
     !dimension ij1(ndij),ije1(0:ndij),itmp(ndij)
      dimension ij1(ndij),ije1(0:ndij) ! note: itmp not used here ! 

      if (getiprint() >= 5) then ! control printout
        call header('Output from "tchng6"',-1)
        write(LUPRI,*) 'ntapef=',ntapef
        write(LUPRI,*) ' iparm=',iparm
      endif

      call dimx (ije(nij),ndp3,"tchng6  ")
      if (iparm.eq.2) then
        nij1=0
        iirp=z'ffff'
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
          iirp=z'ffffffff'
          do 40 j=ije(i-1)+1,ije(i)
            ii=iand(ishft(igen(istrt+j),-16),z'ff')
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
      if (getiprint() >= 5) then ! control printout
        write(LUPRI,*) 'nij,ij(nij),ije(nij)=',nij,ij(nij),ije(nij)
      endif
      return
      end subroutine tchng6

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine tchng8 (nij,ij,ije,ntapef) 
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use ctmp_mod
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
!  ... needed for static alloc of ij(ndij),ije...
!#include "param_ccsd.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1iiia(200),iiva(200)
#include "gen.h"
!      common // gen(nbig)
!      common /i/ igen(nbig)
#include "ctmp.h"
!     common /ctmp/ tmp(ndp4),itmp(ndp4)
#endif

!!!!!!!!!!! need dynamic alloc of ije1 !!!!!!!!!!
      dimension ij(ndij),ije(0:ndij),ije1(0:ndij)

      if (getiprint() >= 5) then ! control printout
        call header('Output from "tchng8"',-1)
        write(LUPRI,*)  'nij=',nij
        !call check_ctmp_mod
      endif

      call dimx (ije(nij),ndp3,"tchng8  ")
      nij1=0
      istrt=nstrt(1,ntapef-1)+2*nij+4
      do 10 i=1,nij
        iirp=z'ffffffff'
        do 20 j=ije(i-1)+1,ije(i)
          ii=iand(ishft(igen(istrt+j),-16),z'ffff')
          if (ii.ne.iirp) then
            nij1=nij1+1
            iirp=ii
            itmp(nij1)=iand(ishft(ij(i),-8),z'ffffff')
            jj=iand(ij(i),z'ff')
            itmp(ndp3+nij1)=ior(ishft(jj,16),ii)
            ije1(nij1-1)=j-1
          end if
20      continue
10    continue

      !call check_ctmp_mod

      call dimx (nij1,ndij,"tchng8  ")
      ije1(nij1)=ije(nij)
      nij=nij1
      call kcopy (nij+1,ije1(0),ije(0))

      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine tmpcng(nij,ij,ije,ntapef,nf,nt)
!------------------------------------------------------------------------
!
!           change nf-index grouping to nt-index grouping
!
!           entering integer arrays: ij, ije
!
!------------------------------------------------------------------------
      use print_level
#if defined HSFS_DYN_ALLOC
      use cfiles_mod
      use gen_mod
      use hsfs_parameters
      implicit none
#include "priunit.h"
! ... needed for static alloc ij(nij),ije...
!#include "param_ccsd.h"
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
#endif

   integer :: nij,nt,nf
   integer :: ntapef
   integer :: nij1,iirp,i,i12,istrt,nsize

!     dimension ij(nij),ije(0:nij)
! mi - to fix index leak 
      !dimension ij(*),ije(0:*)
       integer :: ij(*),ije(0:*)
!!!!!!!!!!!!!! needs dynamic alloc of "ij1, ije1" !!!!!!!!!!!!!!!!!!!
      !dimension ij1(ndij),ije1(0:ndij)
      integer :: ij1(ndij),ije1(0:ndij)

!!!! miro: changed to general integer !!!
      !integer*4 ch1,ch2
      integer :: ch1,ch2

      if (getiprint() >= 5) then ! control printout
        !call header('Output from "tmpcng"',-1)
        write(LUPRI,*) 'Output from "tmpcng" - nij,nt,nf:',nij,nt,nf
        write(LUPRI,*) 'ntapef=',ntapef
      endif
      
      if (nt.eq.nf) then
         if (getiprint() >= 5) then 
           write(LUPRI,*) 'tmpcng exit (nt=nf): nt,nf:',nt,nf
         endif
         return
      endif

      if (nt.gt.nf) go to 310

      nij1=0
      iirp=z'ffff'
      do 160 i=1,nij
      i12=ishft(ij(i),-8*(nf-nt))
      if (i12.ne.iirp) then
        nij1=nij1+1
        ij1(nij1)=i12
        iirp=i12
        ije1(nij1-1)=ije(i-1)
      endif

! miro: important printout to check arrays !
      if (getiprint() >= 15) then 
        write(LUPRI,*) 'i,ij(i),nij1,ije(i-1)',i,ij(i),nij1,ije(i-1)
      endif

160   continue

      ije1(nij1)=ije(nij)

!    refresh nij
      nij=nij1  

!     write(LUPRI,*) 'tmpcng: nij=',nij

      call kcopy(nij,ij1,ij)
      call kcopy(nij+1,ije1(0),ije(0))

      if (getiprint() >= 5) then 
        write(LUPRI,*) 'tmpcng exit 1: nij,ije(nij)=',nij,ije(nij)
      endif
      return

310   istrt=nstrt(1,ntapef-1)+2*nij+4
      nsize=igen(istrt-2)
      iirp=z'ffff'
      nij=0
      do 320 i=1,nsize
      i12=ishft(igen(istrt+i),8*nt-32)
      if (i12.ne.iirp) then
        nij=nij+1
!mi:  runtime-check complains due to "ij" index leak !
        ij(nij)=i12
        iirp=i12
!mi:  runtime-check complains due to "ije" index leak !
        ije(nij-1)=i-1
      endif
320   continue

      ije(nij)=nsize
      call dimx(nij,ndij,"tmpcng")

      if (getiprint() >= 5) then 
        write(LUPRI,*) 'tmpcng exit 2: nij,ije(nij)=',nij,ije(nij)
      endif

      return
      end subroutine tmpcng

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function v1det (ja,jb)
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "dat.h"
!      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif

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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function v2det (ja,jb,jc,jd)
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"

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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function v3det (ja,jb,jc,jd,je,jf) 
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use prm_mod
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "prm.h"
!     common /prm/ iperm(0:6,720)
#include "dat.h"
!     common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif
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
          v3det=v3det+iperm(0,i)*find6(jj1(1),jj1(2),jj1(3),            &
     jj1(iperm(1,i)+3),jj1(iperm(2,i)+3),jj1(iperm(3,i)+3))
        end if
10    continue
      if ((ic(jj1(1)).eq.2).and.(ic(jj1(2)).eq.2)                       &
 .and.(ic(jj1(3)).eq.2)) then
        v3det=-1.d0*v3det
      end if
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function v4det (ja,jb,jc,jd,je,jf,jg,jh) 
!------------------------------------------------------------------------
!
!  Calculates matrix element between two Slater determinants ???
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use prm_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "prm.h"
!      common /prm/ iperm(0:6,720)
#endif

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
          v4det=v4det+iperm(0,i)*find8(jj1(1),jj1(2),jj1(3),jj1(4),     &
                     jj1(iperm(1,i)+4),jj1(iperm(2,i)+4),               &
                     jj1(iperm(3,i)+4),jj1(iperm(4,i)+4))                    
          end if
10    continue
      return
      end 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function v5det (ja,jb,jc,jd,je,jf,jg,jh,ji,jj)
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use prm_mod
      use dat_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "prm.h"
!      common /prm/ iperm(0:6,720)
#include "dat.h"
!      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
#endif

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
          v5det=v5det+iperm(0,i)*find10(jj1(1),jj1(2),jj1(3),           &
     jj1(4),jj1(5),jj1(iperm(1,i)+5),jj1(iperm(2,i)+5),                 &
     jj1(iperm(3,i)+5),jj1(iperm(4,i)+5),jj1(iperm(5,i)+5))
        end if
10    continue
      if ((ic(jj1(1)).eq.2).and.(ic(jj1(2)).eq.2).and.                  &
 (ic(jj1(3)).eq.2).and.(ic(jj1(4)).eq.2).and.                           &
 (ic(jj1(5)).eq.2)) then
        v5det=-1.d0*v5det
      end if
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function v6det (ja,jb,jc,jd,je,jf,jg,jh,ji,jj,jk,jl)
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
#if defined HSFS_DYN_ALLOC
      use prm_mod
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "prm.h"
!      common /prm/ iperm(0:6,720)
#endif
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
          v6det=v6det+iperm(0,i)*find12(jj1(1),jj1(2),jj1(3),jj1(4),    &
     jj1(5),jj1(6),jj1(iperm(1,i)+5),jj1(iperm(2,i)+5),                 &
     jj1(iperm(3,i)+5),jj1(iperm(4,i)+5),jj1(iperm(5,i)+5),             &
     jj1(iperm(6,i)+5))                                         
          end if
10    continue
      return
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine print_cijx_stuff
!---------------------------------------------------------------------------------------
!                       control printout
!---------------------------------------------------------------------------------------
!     use print_level
#if defined HSFS_DYN_ALLOC
      use cij1_mod 
      use cij2_mod 
      use cij3_mod 
      use cij4_mod 
      use cij5_mod 
      use cij6_mod 
!      use cfiles_mod 
!      use gen_mod 
      use hsfs_parameters
#include "implicit.h"
#include "priunit.h"
#else
!     implicit real*8 (a-h,o-z)
#include "implicit.h"
#include "priunit.h"
#include "param_ccsd.h"
#include "cij1.h"
!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
#include "cij2.h"
!     common /cij2/ nij2,ij2(ndij),ije2(0:ndij),iii2,iiv2
#include "cij3.h"
!     common /cij3/ nij3,ij3(ndij),ije3(0:ndij),iii3,iiv3
#include "cij4.h"
!     common /cij4/ nij4,ij4(ndij),ije4(0:ndij),iii4,iiv4
#include "cij5.h"
!     common /cij5/ nij5,ij5(ndij),ije5(0:ndij),iii5,iiv5
#include "cij6.h"
!     common /cij6/ nij6,ij6(ndij),ije6(0:ndij),iii6,iiv6
!#include "cfiles.h"
!      common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!     1 iiia(200),iiva(200)
!#include "gen.h"
!     common // gen(nbig)
!     common /i/ igen(nbig)
#endif
        write(LUPRI,'(4x,a)') 'Content of cijn blocks/modules:'

        if (nij1.ne.0) then
        write(LUPRI,*) 'nij1,ij1(nij1),ije1(nij1),iii1,iiv1=',          &
                   nij1,ij1(nij1),ije1(nij1),iii1,iiv1
        endif

        if (nij2.ne.0) then
        write(LUPRI,*) 'nij2,ij2(nij2),ije2(nij2),iii2,iiv2=',          &
                   nij2,ij2(nij2),ije2(nij2),iii2,iiv2
        endif

        if (nij3.ne.0) then
        write(LUPRI,*) 'nij3,ij3(nij3),ije3(nij3),iii3,iiv3=',          &
                   nij3,ij3(nij3),ije3(nij3),iii3,iiv3
        endif

   !    write(LUPRI,*) 'nij4,ij4(nij4),ije4(nij4),iii4,iiv4=',          &
   !               nij4,ij4(nij4),ije4(nij4),iii4,iiv4
   !    write(LUPRI,*) 'nij5,ij5(nij5),ije5(nij5),iii5,iiv5=',          &
   !               nij5,ij5(nij5),ije5(nij5),iii5,iiv5
   !    write(LUPRI,*) 'nij6,ij6(nij6),ije6(nij6),iii6,iiv6=',          &
   !               nij6,ij6(nij6),ije6(nij6),iii6,iiv6

      end subroutine print_cijx_stuff

