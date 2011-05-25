*deck readcoor
c     **************************************************************
      subroutine readcoor (i, nont, mline, nmline, namnar) 
c     **************************************************************
      implicit none
#include "argosparm.h"
      integer   nont(mxatom)
      character namnar(mxatom,mxatom)*4
      do j = 1, nont(i)
         read (mline(nmline),'(a4)') namnar(i,j)
         IF (MLINE(NMLINE)(1:4) .EQ. '    ')
     &      NAMNAR(i,j) = NAMNAR(i,j-1)
c        # Coordinate will be taken directly from DIRAC
         nmline = nmline + 1  
      enddo
 2699 continue
      end
*deck readbskey 
c     ******************************************************************
      subroutine readbskey (mline, nmline, bset, bskeyword, i, iqm, jco,
     &                      basfil,bscode)
c     ******************************************************************
#include "argosparm.h"
#include "argosbas.h"
      bscode = 0
c
c check large basis set
c
 2710 continue
      read (mline(nmline)(1:5),'(a5)') bset
      if (bset .eq. 'LARGE') then
         goto 2720
      else
         write (*,*) bset
         write (*,*) 'error BSET is not LARGE'
         goto 2791
      endif
c
c check basis set keyword
c
 2720 continue
c     explicit
      if ((mline(nmline)(7:12).eq.'INTGRL')) then
             read (mline(nmline)(1:80),'(6x,a6,12i5)',err=2791)
     &            bskeyword,iqm(i),(jco(i,j),j=1,iqm(i))
             goto 2799
      elseif ((mline(nmline)(7:14).eq.'EXPLICIT')) then
             read (mline(nmline)(1:80),'(6x,a8,12i5)',err=2791)
     &            bskeyword,iqm(i),(jco(i,j),j=1,iqm(i))
             goto 2799
      elseif ((mline(nmline)(7:8).eq.'  ')) then
c            assume it is 'LARGE   2    1....'
             read (mline(nmline)(1:65),'(5x,12i5)',err=2791)
     &            iqm(i),(jco(i,j),j=1,iqm(i))
             bskeyword='NONE'
             goto 2799
c     library 
      elseif ((mline(nmline)(7:11).eq.'BASIS')) then
             read (mline(nmline)(1:80),'(6x,a5,1x,a20)',err=2791)
     &            bskeyword, basfil(i)
             bscode = 5
             goto 2799
c     else
      else
             write (*,*) 'Stopped in readbskey(argosin) routine'
             write (*,*) 'Use INTGRL, EXPLICIT, BASIS keyword',
     &                   'or explicit description of iqm, jco',
     &                   'are also allowed'
             stop     
      endif
      goto 2799
 2791 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readbskey...'
      write(*,*) 'Cannot read basis keyword'
      stop
 2799 continue
      nmline = nmline + 1
      return
      end
*deck readexplibas 
c     ******************************************************************  
      subroutine readexbas(i, mline, nmline, iqm, jco, nucij, nrcij,
     &                     expn, cont, bfree)
c     ******************************************************************
c     read explicit basis set
      implicit none
#include "argosparm.h"
#include "argosbas.h"
c     *** basis lines ***
      integer   nbline, mxbline
      parameter (mxbline=150)
      character bline(mxbline)*90
      integer   p
      do j = 1, iqm(i)                            !j for
         do k = 1, jco(i,j)                       !k for
c           pass (i,j,k)th block
            nbline = 1
            read(mline(nmline),'(a80)') bline(nbline)
            read(bline(nbline),'(a1,i4)') frmt, nucij(i,j,k)
            nmline = nmline + 1
            nbline = nbline + 1
            do p = 1, nucij(i,j,k)
               read(mline(nmline),'(a80)') bline(nbline)
               nmline = nmline + 1
               nbline = nbline + 1
            enddo
c           get expn & cont from (i,j,k)th block
            basexpl = .true.
            nbline  = 1
            call acporbag(bline,nbline,i,j,k,iqm,jco,nucij,nrcij,
     &                    expn,cont,basexpl,bfree)
         enddo
      enddo
c     write(*,*) 'end of basis'
      return
      end
*deck readlibbas 
c     **************************************************************
      subroutine readlibbas(basfil, i, iqm, jco, q,
     &                      nucij, nrcij, expn, cont, bfree)
c     **************************************************************
#include "argosparm.h"
#include "argosbas.h"
c     *** basis lines ***
      integer   nbline, mxbline
      parameter (mxbline=150)
      character bline(mxbline)*90

      real*8    q(mxatom)
      integer   intq, nucel, bdstart, bdend, bdcolon
      character string*300, sign1*1, basdir*600
      logical   exst
      integer   t
c
      basexpl = .false.
c     check basis directory
      call getenv ('BASDIR', basdir)
      bdstart = 1
 2803 continue 
c     (find next colon in basdir)
      bdcolon = bdstart - 1 + index (basdir(bdstart:), ':' )
c     if there are no colons, we are the last directory in basdir
      if ( bdcolon .eq. bdstart - 1 ) then
c        find end of string
         bdcolon = bdstart - 1 + index (basdir(bdstart:), ' ' )
c        unless find directory
         if (bdcolon .eq. bdstart) goto 2894
      endif
      bdend = bdcolon - 1
      string = basdir(bdstart:bdend)
      bdstart = bdcolon + 1
c     (add trailing slash if not present)
      bdend = index (string(1:), ' ') - 1
      if (string(bdend:bdend).ne.'/') then
         bdend = bdend + 1 
         string(bdend:bdend) = '/'
      end if
      string = string(1:bdend)//basfil(i)
      bdend = index (string(1:),' ') - 1
c     (inquire if "/path/basis-set" exists.)
      inquire (file=string, exist=exst)
c     if (.not.exst) goto 2803
      if (exst) then
         goto 2810
      else
         goto 2803
      endif   
c     check element
 2810 continue
      lufile6=26
      open (lufile6, file=string, status='old', err=2891)
      rewind (lufile6)
 2811 continue
      nbline = 1
      read (lufile6,'(a90)') bline(nbline)
      if ((bline(nbline)(1:1).eq.'a').or.
     &    (bline(nbline)(1:1).eq.'A')) then
         intq = int(q(i))
         read (bline(nbline),'(bn, a1, i4)') sign1, nucel
         if (nucel.eq.intq) then
            goto 2820
         endif
         goto 2811
      else
         goto 2811
      endif
c
c     read basis file
c
 2820 continue
      read (lufile6,'(a90)') bline(nbline)
c     check end of reading
      if ((bline(nbline)(1:1).eq.'a').or.
     &    (bline(nbline)(1:1).eq.'A')) goto 2821
c     check commentary
      if ((bline(nbline)(1:1).eq.'!').or.
     &    (bline(nbline)(1:1).eq.'$').or.
     &    (bline(nbline)(1:1).eq.'#')) goto 2820
      nbline = nbline + 1
      goto 2820
 2821 continue
      close(lufile6)
c         
c     get exponent contration
c
      bline(nbline)(1:3)='eof'
      do t = 1, nbline
      enddo
      nbline = 1
      do j = 1, mxiqm
         k = 1  !in case of library
         if (bline(nbline)(1:3).eq.'eof') goto 2822
         call acporbag(bline,nbline,i,j,k,iqm,jco,nucij,nrcij,
     &                 expn,cont,basexpl,bfree)
      enddo
c
c     set iqm, jco for library
c
 2822 continue
      iqm(i) = j - 1
      do j = 1, iqm(i)
         jco(i,j) = 1
      enddo
      goto 2899
c
c     error message
c
 2891 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readlibbas...'
      write(*,*) 'Wrong path or wrong name'
      stop
 2892 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readlibbas...'
      write(*,*) 'Cannot read from basis file'
 2893 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readlibbas...'
      write(*,*) 'Unsupported element for this basis'
      stop
 2894 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readlibbas...'
      write(*,*) 'Cannot find basis directory'
      stop
 2899 continue
      end
*deck acporbag
c     **********************************************************
      subroutine acporbag(bline,nbline,i,j,k,iqm,jco,nucij,nrcij,
     &                    expn,cont,basexpl,bfree)
c     **********************************************************
#include "argosparm.h"
#include "argosbas.h"
c     *** basis lines ***
      integer   nbline, mxbline
      parameter (mxbline=150)
      character bline(mxbline)*90

      integer   lenexp, temp1, p 
      read (bline(nbline)(1:15),'(a1,i4,i5)')
     &     frmt, nucij(i,j,k), nrcij(i,j,k)
      nbline = nbline + 1
c
c     check format
c
      if (frmt.eq.'F' .or. frmt.eq.'f') then
         goto 100
      elseif (frmt.eq.'H' .or. frmt.eq.'h') then
c        free(i,j,k) = .false.
         if (nrcij(i,j,k).ge.4) then
            write(*,*) 'Basis Set Error :'
            write(*,*) 'For ARGOS routine, orbital exponent and '
            write(*,*) 'contraction coefficient should be in same line '
            write(*,*) 'For high precision, fix number of contraction  '
            write(*,*) 'less than 4 for each block or use free format!!'
            stop ! *** fortran stopped ***
         endif
      else
         if (nrcij(i,j,k).ge.7) then
            write(*,*) 'Basis Set Error :'
            write(*,*) 'For ARGOS routine, orbital exponent and '
            write(*,*) 'contraction coefficient should be in same line'
            write(*,*) 'Fix number of contraction less than 7 for each'
            write(*,*) 'block or use free format !!!'
            stop ! *** fortran stopped ***
         endif
      endif
c
c     check contraction
c
 100  continue
      if(nrcij(i,j,k).ne.0) then
        do l = 1, nucij(i,j,k)
           read(bline(nbline),'(a80)') bfree(i,j,k,l)
           nbline = nbline + 1
        enddo 
c
c     no contraction
c
      else
c
        nrcij(i,j,k) = nucij(i,j,k)
        do 300 l = 1, nucij(i,j,k)
c          read primitive
           temp1=0
           do p=1,80
            if (bline(nbline)(p:p).ne.' ') then
               temp1 = temp1 + 1
               expnc(temp1:temp1)=bline(nbline)(p:p)
            endif
           enddo
c          make contraction
           do m = 1, nrcij(i,j,k)
c             cont(i,j,k,l,m) = D0
              contc(l,m) = '0'
c             if (l.eq.m) cont(i,j,k,l,m) = D1
              if (l.eq.m) contc(l,m) = '1'
           enddo
c          write contraction
           do m = 1, nrcij(i,j,k)
              sumcontc(2*m-1:2*m)=' '//contc(l,m)
           enddo
           btemp(1:1)=' '
           btemp(2:temp1+1)=expnc(1:temp1)
           btemp(temp1+2:temp1+1+2*nrcij(i,j,k))
     &          =sumcontc(1:2*nrcij(i,j,k))
           bfree(i,j,k,l)
     &            = btemp(1:temp1+1+2*nrcij(i,j,k))
           nbline = nbline + 1
  300   continue
      endif
      end
c
*deck readecpexp  
c     *******************************************************
      subroutine readecpexp(mline, nmline, i, namatm, namnar,
     &                      lcru, llsu, nbfcr, nbfls, soindex,
     &                      ncr, zcr, ccr, nls, zls, cls,
     &                      mcrs, nmcrs, nontyp)
c     *******************************************************
#include "argosparm.h"
#include "argosecp.h"
#include "argosatm.h"
c
c     for SO-index
c
      csoindex = 0
      read(mline(nmline)(5:8),'(a4)') namatm(i) ! name of atom
      write(6,*) 'yc-xxxx', namatm(i)
      nmline = nmline + 1
      do j = 1, nontyp
         if (namatm(i)(1:3).eq.namnar(j,1)(1:3)) then
            if (mcrs(j).ne.0) goto 3227
            mcrs(j) = soindex
            nmcrs = nmcrs + 1
            csoindex = 1
            goto 3221
         endif
      enddo
      if (csoindex.ne.1) goto 3226
c
c     read ARPP / SOPP
c
 3221 continue
      read(mline(nmline),*,err=3231) lcru(i), llsu(i) ! read number of ARPP/SOPP
      nmline = nmline + 1
      do j = 1, lcru(i) + 1                           ! read ARPPs
         read (mline(nmline),*,err=3232) nbfcr(i,j)
         nmline = nmline + 1
         do k = 1, nbfcr(i,j)
            read (mline(nmline),*,err=3233) 
     &           ncr(i,j,k),zcr(i,j,k),ccr(i,j,k)
            nmline = nmline + 1
         enddo
      enddo
      do j = 1, llsu(i)                      ! read SOPPs
         read (mline(nmline),*,err=3234) nbfls(i,j)
         nmline = nmline + 1
         do k = 1, nbfls(i,j)
            read (mline(nmline),*,err=3235) 
     &           nls(i,j,k),zls(i,j,k),cls(i,j,k)
            nmline = nmline + 1
         enddo
      enddo
      goto 3999
c
c     Error
c
 3226 continue ! namatm
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecpexp...'
      write(*,'(a9,a4,a15)') 'ECP Name ',namatm(i),' dose not match'
      write(*,*) 'It is case sensitive'
      stop
 3227 continue ! namatm
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecpexp...'
      write(*,'(a8,a4,a20)') 'ECP for ',namatm(i),' is already declared'
      stop
 3231 continue ! lcru, llsu
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecpexp...'
      write(*,*) 'Cannot read number of ARPP or SOPP'
      write(*,'(a14,a4)') 'Check ECP for ',namatm(i)
      stop
 3232 continue ! nbfcr
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecpexp...'
      write(*,*) 'Cannot read number of function'
      write(*,*) 'in potential expansion'
      write(*,'(a14,a4)') 'Check ECP for ',namatm(i)
      stop
 3233 continue ! ncr, zcr, ccr
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecpexp...'
      write(*,*) 'Cannot read exponent and coefficient'
      write(*,'(a14,a4)') 'Check ECP for ',namatm(i)
      stop
 3234 continue ! bnfls
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecpexp...'
      write(*,*) 'Cannot read number of function'
      write(*,*) 'in Spin-Orbit potential expansion'
      write(*,'(a14,a4)') 'Check ECP for ',namatm(i)
      stop
 3235 continue ! nls, zls, cls
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecpexp...'
      write(*,*) 'Cannot read exponent and coefficient'
      write(*,'(a14,a4)') 'Check ECP for ',namatm(i)
      stop
 3999 continue
      end
c
*deck readecplib
c     ********************************************************
      subroutine readecplib(mline, nmline, i, namatm, namnar,
     &                      lcru, llsu, nbfcr, nbfls, soindex,
     &                      ncr, zcr, ccr, nls, zls, cls,
     &                      mcrs, nmcrs, nontyp)
c     ********************************************************
#include "argosparm.h"
#include "argosecp.h"
#include "argosatm.h"
      parameter (mxeline=100)
      integer   ecpidtmp, 
     &          ecpidn, ncore, ncoretmp, bdstart, bdend, bdcolon,
     &          epend
      character eline(mxeline)*90, basdir*300,
     &          ecpid*1, ecpso*2, ecpsotmp*2, 
     &          ecpmthd*3, ecpmthdtmp*3, string*300
      logical   ecpidl,exst
c
c     read keyword
c
      read (mline(nmline)(8:10),*) namatm(i)
      read (mline(nmline),'(11x,a1)') ecpid
          if (ecpid.eq.'I') then
             read(mline(nmline),'(bn,12x,i3)') ecpidn
             ecpidl=.true.
          elseif (ecpid.eq.'C') then
             read(mline(nmline),'(bn,12x,i3,1x,a2,1x,a3)')
     &           ncore,ecpso,ecpmthd
             ecpidl=.false.
          else
             goto 2897
          endif
      nmline = nmline + 1
c
c     find element 
c
c     check basis directory
      call getenv ('BASDIR', basdir)
c     print *, 'BASDIR=',basdir
      bdstart = 1
 2803 continue 
c     find next colon in basdir
      bdcolon = bdstart - 1 + index (basdir(bdstart:), ':' )
c     if there are no colons, we are the last directory in basdir
      if ( bdcolon .eq. bdstart - 1 ) then
c        find end of string
         bdcolon = bdstart - 1 + index (basdir(bdstart:), ' ' )
c        unless find directory
         if (bdcolon .eq. bdstart) goto 2894
      endif
      bdend = bdcolon - 1
      string = basdir(bdstart:bdend)
      bdstart = bdcolon + 1
c     (add trailing slash if not present)
      bdend = index (string(1:), ' ') - 1
      if (string(bdend:bdend).ne.'/') then
         bdend = bdend + 1 
         string(bdend:bdend) = '/'
      end if
c     ecpname check
      epend = 1
 2804 continue
      if (namatm(i)(epend:epend).ne.' ') then
         epend = epend + 1
         goto 2804
      endif
      epend = epend - 1
c
      string = string(1:bdend)//namatm(i)(1:epend)//'.ecp' ! file
      bdend = index (string(1:),' ') - 1
c     (inquire if "/path/basis-set" exists)
      inquire (file=string, exist=exst)
      if (.not.exst) goto 2803
      if (exst) then
c        write (*,*) 'we got string'
c        write (*,*) srting 
         goto 2810
      endif
c   
c     open file
c
 2810 continue
      lufile6=26
      open (lufile6, file=string, status='old', err=2891)
      rewind (lufile6)
 2811 continue
      neline = 1
      read(lufile6,'(a80)',iostat=ioerr,err=2892,end=2893)eline(neline)
      if (eline(neline)(1:3).eq.'eof') goto 2892 ! No such ecp found
      if ((eline(neline)(1:3).eq.'ecp').or.
     &    (eline(neline)(1:3).eq.'ECP')) then
c          
         if (ecpidl .eqv. .true.) then
            read(eline(neline)(6:8),'(i3)') ecpidtmp
            if (ecpidtmp.eq.ecpidn) goto 2819
            goto 2811  
         else
            read(eline(neline)(11:20),'(i3,1x,a2,1x,a3)')
     &           ncoretmp, ecpsotmp, ecpmthdtmp  
            if ((ncoretmp.eq.ncore).and.
     &          (ecpsotmp.eq.ecpso).and.
     &          (ecpmthdtmp.eq.ecpmthd)) goto 2819
            goto 2811
         endif
      else
         goto 2811
      endif
c
c     read ecp lines
c
 2819 continue
      neline = 2
 2820 continue
      read(lufile6,'(a90)') eline(neline)
c     check commentary
      if ((eline(neline)(1:1).eq.'!').or.
     &    (eline(neline)(1:1).eq.'$').or.
     &    (eline(neline)(1:1).eq.'#')) then
         goto 2820
      endif
      if ((eline(neline)(1:3).eq.'ECP').or.
     &    (eline(neline)(1:3).eq.'eof')) goto 2821
      neline = neline + 1
      goto 2820
c
c     number of ARPP,SOPP
c
 2821 continue
      lcru(i) = 0
      llsu(i) = 0
      neline  = 2
      read(eline(neline),*) lcru(i), llsu(i)
      neline = neline + 1
c     for SO potential index
      csoindex = 0
      do j = 1, nontyp
         if (namatm(i)(1:3).eq.namnar(j,1)(1:3)) then
            if (mcrs(j).ne.0) goto 2895 
            mcrs(j) = soindex
            nmcrs = nmcrs + 1
            csoindex = 1
            goto 2822
         endif
      enddo
      if (csoindex.ne.1) goto 2896
c     read ARPP/SOPP
 2822 continue
      if (lcru(i).gt.0) then
         do j = 1, lcru(i)+1         ! ARPP
            read(eline(neline),*) nbfcr(i,j)
            neline = neline + 1
            do k = 1, nbfcr(i,j) 
               read(eline(neline),*) ncr(i,j,k),zcr(i,j,k),ccr(i,j,k)
               neline = neline + 1
            enddo
         enddo
      endif
      if (llsu(i).gt.0) then
         do j = 1, llsu(i)           ! SOPP
            read (eline(neline),*) nbfls(i,j)
            neline = neline + 1
            do k = 1, nbfls(i,j)
               read (eline(neline),*) nls(i,j,k),zls(i,j,k),cls(i,j,k)
               neline = neline + 1
            enddo
         enddo
      endif
      goto 2899 
c
c     error message
 2891 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecplib...'
      write(*,*) 'Cannot open file'
      stop
 2892 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecplib...'
      write(*,*) 'Cannot read ecp from basis file'
      stop 
 2893 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecplib...'
      write(*,*) 'Unsupported type for this ECP'
      stop
 2894 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecplib...'
      write(*,*) 'ECP directory :', string
      write(*,*) 'Cannot find ecp directory'
      stop
 2895 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecplib...'
      write(*,*) 'ECP for this atom is already declared'
      stop
 2896 continue       
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecplib...'
      write(*,*) 'Atom names are not match' 
      stop
 2897 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'in subroutine readecplib...'
      write(*,*) 'Specify "C" for Core index '
      write(*,*) '     or "I" for Index number'
      stop
 2899 continue
c     write (*,*) 'end of reading ecp library'
      close(lufile6)
      end
c
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     function
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*deck sqdif
      function sqdif( a, b, n )
c     # returns the sum of the squares of the differences between
c     # the elements of the vectors a(:) and b(:)
      implicit none
      integer n
      real*8 a(n), b(n), sqdif
      integer it
      real*8 zero
      parameter(zero=0d0)
      sqdif = zero
      do it = 1, n
         sqdif = sqdif + (a(it) - b(it))**2
      enddo
      return
      end
c
*deck chkopt 
c     *******************************************************
      subroutine chkopt ( ns, naords, ncons, ngcs,
     &                    icut, inrm,
     &                    ncrs, 
     &                    symtxt, nontyp, nont, iqm, jco, nmcrs)
c     *******************************************************
c     card 2 : options
c     (dirac symmetry is only low level symmetry for argos)
#include "argosparm.h"
#include "argosecp.h"
#include "argosatm.h"
#include "argosopti.h"
#include "argoscom.h"
#include "argoscomch.h"
      integer iqm(mxatom), jco(mxatom,mxiqm)
c     ****ngen****
c     'ngen' is substituted to 'ngenp'
c     check subroutine permgen 
c     *****ns*****
      ns     = nsag
c     ***naords***
      naords = 0
      do i = 1, ns
         do j = 1, nont(i)
            do k = 1, iqm(i)
               naords = naords + 1
            enddo
         enddo
      enddo
c     ***ncons***
      ncons  = 0
c     do i = 1, ns
      do i = 1, nontyp
         do j = 1, iqm(i)
            ncons = ncons + jco(i,j)
         enddo
      enddo
c     ***ngcs****
      ngcs   = naords
      icut   = 15
      inrm   = 0       
c     ***ncrs****
c     number of total ECP center
      ncrs   = nmcrs   
c     ***********
      end

*deck permatm
c     *****************************************************
      subroutine permatm (isend,nontyp,nont,ngenl,
     &                    xyz3,mj2,symtxt,kasym)
c     *****************************************************
#include "argosparm.h"
#include "argosatm.h"
      integer   lgen(3,3), nirep, ltrans(3,8),
     &          iclat, natm, iopl, icentar, j2, j3, lth
      real*8    distsq, nclose, sqdif 
c        
c     read rotation 
c           
      do i = 1, symtxt
         do j = 1, 3
            lgen(j,i) = 1
         enddo 
         do j = 1, 3
            if(kasym(i,j).eq.'X') lgen(1,i)=-1
            if(kasym(i,j).eq.'Y') lgen(2,i)=-1
            if(kasym(i,j).eq.'Z') lgen(3,i)=-1
         enddo
      enddo
c
c     generate symetry point group operations
c     
 20   continue
      ngenl = symtxt 
      do j=1,3
         ltrans(j,1)=1
         if (ngenl.ge.1) then
            ltrans(j,2)=lgen(j,1)
            if (ngenl.ge.2) then
               ltrans(j,3)=lgen(j,2)
               ltrans(j,4)=lgen(j,2)*lgen(j,1)
               if (ngenl.ge.3) then
                  ltrans(j,5)=lgen(j,3)
                  ltrans(j,6)=lgen(j,3)*lgen(j,1)
                  ltrans(j,7)=lgen(j,3)*lgen(j,2)
                  ltrans(j,8)=lgen(j,3)*lgen(j,2)*lgen(j,1)
               endif
            endif
         endif
      enddo
c
c     set coordinate
c
 30   continue
      icentar=1
      do i=1,nontyp
         natm=nont(i)
         do j=1,natm
            do k=1,3
               xyz1(k,j,i)=get_cord(k,icentar) ! get CORD from DIRAC
               if(symtxt.eq.0)then
                 xyz3(k,1,j,i)=xyz1(k,j,i)
                 mj2(i,j)=2
               endif
            enddo
            icentar=icentar+1
         enddo
      enddo
      if(symtxt.eq.0) goto 999
c
c     operation on to origin coordinate
c
      nirep = 2 ** symtxt
      do j2= 1, nirep
      enddo
c
      do i = 1, nontyp       ! number of atom type
         natm = nont(i)
         do j = 1, natm      ! number of symmetry-distinct atom of this type
            mj2(i,j)=1       !  
            do j2=1,nirep
               do k=1,3
                  xyz2(k,j2,j,i)=ltrans(k,j2)*xyz1(k,j,i)
                  xyz4(k)=xyz2(k,j2,j,i)
               enddo
c              check 
               do lth=1,j2-1
                  distsq=sqdif(xyz4(1),xyz2(1,lth,j,i),3)
                  if (distsq.lt.tolsq) then
                     goto 50
                  endif
               enddo
               do k=1,3
                  j3=mj2(i,j)
                  xyz3(k,j3,j,i)=xyz2(k,j2,j,i)
               enddo
               mj2(i,j)=mj2(i,j)+1
 50            continue
            enddo
         enddo
      enddo
c
c  operation  
c
      do i=1,nontyp
         do iopl=1,nirep
            do j=1,nont(i)
               do j2=1,mj2(i,j)-1  ! natom(nc)
                  do k=1,3
                     xyz4(k)=ltrans(k,iopl)*xyz3(k,j2,j,i)
                  enddo
c                 # check closest atom 
                  nclose=10000
                  do j3=1,mj2(i,j)-1
                     distsq=sqdif(xyz4(1),xyz3(1,j3,j,i),3)
                     if (distsq.lt.nclose) then
                        nclose=distsq
                        iclat=j3    ! closest atom so far
                     endif
                  enddo
                  isend(j2,j,iopl,i)=iclat
               enddo
            enddo
         enddo
      enddo
 999  continue
      end
c
*deck permgen
c     *******************************************************
      subroutine permgen (isend,ngenl,ngenp,nont,nsag,mj2)
c     *******************************************************
#include  "argosparm.h"
#include  "argosatm.h"
      integer   nsag, is, iatm
c
      ngenp=0
      igenq(1)=2
      igenq(2)=3
      igenq(3)=5
      do jgen=1,ngenl              !ngenl=symtxt
         igen=igenq(jgen)
         do is=1,nsag              !nsag=ns
            do j=1,nont(is)
               do iatm=1,mj2(is,j)-1          !mj2(i,j)-1 = nc(is)
                  if (isend(iatm,j,igen,is).ne.iatm) goto 550
               enddo
            enddo
         enddo
         goto 559
 550     continue               ! this is not the identity
         do kgen=1,ngenp
            kgen1=igenp(kgen)
            do is=1,nsag
               do j=1,nont(is)
                  do iatm=1,mj2(is,j)-1
                     if(isend(iatm,j,igen,is).ne.isend(iatm,j,kgen1,is))
     &               goto 551
                  enddo
               enddo
            enddo
            goto 559
 551        continue            ! this is not generator kgen
         enddo
         if (ngenp.ge.2) then
            kgen1=igenp(1)
            kgen2=igenp(2)
            do is=1,nsag
               do j=1,nont(is)
                  do iatm=1,mj2(is,j)-1
                     if (isend(iatm,j,igen,is).ne.
     &                isend(isend(iatm,j,kgen2,is),j,kgen1,is)) goto 552
                  enddo
               enddo
            enddo
            goto 559
         endif
 552     continue     ! this is not a product of previous generators
         ngenp=ngenp+1
         igenp(ngenp)=igenq(jgen)
 559     continue
      enddo
      end

!     *******************************************************
      SUBROUTINE RECP_INPUT_READ(FILE_NUM, INPUT_LINE)
      IMPLICIT NONE
      INTEGER   FILE_NUM
      CHARACTER INPUT_LINE*90

 1000 CONTINUE    
      READ(FILE_NUM,FMT='(a90)',end=9000) INPUT_LINE
      if ((INPUT_LINE(1:1).eq.'!') .or. 
     &    (INPUT_LINE(1:1).eq.'$') .or.
     &    (INPUT_LINE(1:1).eq.'#')) goto 1000
      if ((INPUT_LINE(1:6).eq.'FINISH').or.
     &    (INPUT_LINE(1:6).eq.'finish')) goto 9000
 9000 CONTINUE
      CLOSE(FILE_NUM)
 9999 CONTINUE
      END


