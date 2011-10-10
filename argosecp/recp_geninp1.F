      subroutine argosin
!
!     read DIRAC input files and convert to ARGOS input file
!

      USE RECP_NTR

#include "mxcent.h"          
#include "nuclei.h"          

#include "argosparm.h" 
#include "argosbas.h"  
#include "argosecp.h"    
#include "argosatm.h"    
#include "argosopti.h"   
 
#include "argoscom.h"    
#include "argoscomch.h"  
#include "pgroup.h"          

      integer dimmr, dimmc 
      integer blkcnt
      LOGICAL FILE_EXIST
      INTEGER ECPLINE, MAX_ECPLINE
      MAX_ECPLINE = 1000
 
! files 
      lufile2 = 9  !molecule.inp
      lufile3 = 51 !argosin

      IF (PRINT_LEVEL.GE.5)
     &   PRINT *,'* ECPINT : entering... ARGOS_GENINPUT'

!     # Delete RECP_INPUT file if exist
      INQUIRE(FILE='RECP_INPUT',EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
         OPEN (UNIT=51,FILE='RECP_INPUT',FORM='FORMATTED',STATUS='OLD')
         CLOSE(UNIT=51,STATUS='DELETE')
      ENDIF

 
!****************************************
!**                                    **
!**      read 'molecule.mol' file      **
!**                                    **
!****************************************
 
 2000 continue
      OPEN (lufile2,file='MOLECULE.MOL',status='old',err=5200)
      rewind (lufile2)
      nmline = 1
 2010 continue 
      read(lufile2,fmt='(a90)',end=2011) mline(nmline)
      if ((mline(nmline)(1:1).eq.'!').or.
     &    (mline(nmline)(1:1).eq.'$').or.
     &    (mline(nmline)(1:1).eq.'#')) goto 2010 
      if ((mline(nmline)(1:6).eq.'FINISH').or.
     &    (mline(nmline)(1:6).eq.'finish')) then
 2011    continue
         close(lufile2)
         goto 2200
      endif
      nmline = nmline + 1 
      goto 2010
 
!****************************************
!**                                    **
!**                                    **
!**      processing data .....1/2      **
!**                                    **
!**                                    **
!****************************************
 2200 continue
!
! Section 4. Point group generator
!
      nmline = 4
      read (mline(nmline)(1:80),'(bn,a1, i4, i3, i2, 9a1, a1, d10.2)')
     &              crt, RECP_NONTYP, kcharg, symtxt, 
     &              ((kasym(i,j),j=1,3),i=1,3), id3, thrs
      nmline = nmline + 1

      nmcrs = 0  ! index for the core and spin-orbit expansion set
      soindex = 0
!
! Section 5. charge, atoms, exponent   
!
      do i = 1, RECP_NONTYP
         read (mline(nmline)(1:80),'(1x,f9.0,i5,f20.5)') 
     &        q(i),RECP_NONT(i),qexp1
         nmline = nmline + 1
!
! Section 6. Atom type, coordinates
!
!        Coordinate will be taken directly from DIRAC
!        Only atom name will be taken here
         call readcoor(i, RECP_NONT, mline, nmline, namnar)
c
c Section 7. Type of Basis, # of AO blocks
c
         call readbskey(mline,nmline,bset,bskeyword,i,iqm,jco,
     &                  basfil,bscode)
c
c Section 8. Exponents and coeeficients blocks
c
              if (bscode.eq.0)           ! explicit basis (18-1)
     &           call readexbas(i, mline, nmline, iqm, jco,
     &                          nucij, nrcij, expn, cont, bfree)
              if (bscode.eq.5)           ! basis library
     &           call readlibbas(basfil, i, iqm, jco, q,
     &                           nucij, nrcij, expn, cont, bfree)
!        --------
!        Read ECP    
!        --------
         mcrs(i) = 0
         IF (MLINE(NMLINE)(1:3).EQ.'ECP') THEN
            ecpchk = .true.
            if ((mline(nmline)(4:6).ne.'LIB')) THEN
!              (3-2. explicit ECP)
               call readecpexp(mline, nmline, i, namatm, namnar,
     &                      lcru, llsu, nbfcr, nbfls,
     &                      ncr, zcr, ccr, nls, zls, cls,
     &                      mcrs, nmcrs, RECP_NONTYP, q)
            else
!              (3-3. ECP library)
               call readecplib(mline, nmline, i, namatm, namnar,
     &                      lcru, llsu, nbfcr, nbfls, soindex,
     &                      ncr, zcr, ccr, nls, zls, cls,
     &                      mcrs, nmcrs, RECP_NONTYP )
            endif
            soindex = soindex + 1
            mcrs(i) = soindex
            nmcrs   = nmcrs + 1
         ENDIF
!
      enddo

      CALL RECP_NTRD_READCP


      goto 3000
c error
 2991 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'In subroutine argosin...'
      write(*,'(a6,a10,a1)') 'line (',mline(nmline)(1:10),')'
      write(*,*) 'Please check your input'
      write(*,*) '"START ECP", "FINISH" or "ECP" keyword is needed'
      stop
 
 3000 continue
      call chkopt (ns, naords, ncons, ngcs,
     &             icut, inrm, ncrs, 
     &             symtxt, RECP_NONTYP, RECP_NONT, iqm, jco, nmcrs)
c
c     coordinate copy
c     check symmetry operation
c     (for card 2 & 9) 
      call permatm (isend,RECP_NONTYP,RECP_NONT,ngenl,
     &              xyz3,mj2,symtxt,kasym)
      call permgen (isend,ngenl,ngenp,RECP_NONT,nsag,mj2)
 
!****************************************
!     Generate RECP_INPUT
!****************************************
 
 4000 continue
      open (lufile3,file='RECP_INPUT',status='new',err=5300)
 
!     Card 2: Options 
      write (lufile3,'(6i3)') ngenp, nsag, naordsag, ncons, 
     &                         ngcsag, nmcrs
 
!     Card 3: Point group symmetry
      write (lufile3,'(i3,8(i2,1x,a3))')
     &       nstag, (ndag(i), itypag(i),i=0,nstag-1)
 
!     Card 4: NDPT, # of distinct products of irrep
      write(lufile3,'(i3)')ndptag
 
!     Card 5: P1 X P2 = P3
!     ----------------------------------------------------
!     *  THIS NEEDS ATTENTION FOR D2h POINT GROUP !!!!!  *
!     ----------------------------------------------------
      if (ndptag .eq. 0) then
!        write(lufile3,'(a)'),'  0'
      elseif (ndptag .eq. 1) then
         write(lufile3,'(a)'),'  2  3  4'
      elseif (ndptag .eq. 7) then
         write(lufile3,'(a)'),'  2  3  4'
         write(lufile3,'(a)'),'  2  5  6'
         write(lufile3,'(a)'),'  2  7  8'
         write(lufile3,'(a)'),'  3  5  7'
         write(lufile3,'(a)'),'  3  6  8'
         write(lufile3,'(a)'),'  4  5  8'
         write(lufile3,'(a)'),'  4  6  7'
      endif

!     Card 6: AO to SO transformation matrices
      do icent = 1, nsag
!        do k = 1, RECP_NONT(icent) 
            do lang = 1, maxang(icent)
               write(lufile3,'(20i3)'),irowmax(icent,lang),
     &           (ireparg(icent,lang,lmang),lmang=1,irowmax(icent,lang))
            enddo   
!        enddo
      enddo

      blkcnt = 0
      do 51 icent = 1, NUCIND
         do 52 lang = 1, maxang(icent)
            blkcnt = blkcnt + 1
            write(lufile3,'(20i3)'),
     &         irowmax(icent,lang),irowmax(icent,lang),blkcnt
            do 53 dimmr = 1, irowmax(icent,lang)
               write(lufile3,'(20i3)') (nint(CTRANarg(icent,lang,dimmr,
     &               dimmc)),dimmc=1,irowmax(icent,lang))
  53        continue
  52     continue
  51  continue
      deallocate(ctranarg)

c----------------------------------------------
c Card 7 : exp, contracted coeff
c----------------------------------------------
      do i = 1, nsag                  !atom
         do j = 1, iqm(i)             !s,p,d...orbital
            lmnp1(i) = j
            do k = 1, jco(i,j)
               write(lufile3,'(3i5)') 
     &              nucij(i,j,k),lmnp1(i),nrcij(i,j,k)
               do l = 1, nucij(i,j,k)
                  write(lufile3,*) bfree(i,j,k,l) 
               enddo
            enddo
         enddo
      enddo
c----------------------------------------------
c Card 8 : ecp
c----------------------------------------------
      if (ecpchk.EQV..false.) goto 2900
!     do i = 1, soindex-1                   ! (nmcrs)
      do i = 1, soindex                     ! (nmcrs)
         write(lufile3,'(2i5)') lcru(i),llsu(i)
         do j = 1, lcru(i) + 1                 ! ARPP
            write (lufile3,'(i5)') nbfcr(i,j)
            do k = 1, nbfcr(i,j)
               write (lufile3,'(i4,2f15.7)')
     &               ncr(i,j,k),zcr(i,j,k),ccr(i,j,k)
            enddo
         enddo
         do j = 1, llsu(i)                     ! SOPP
            write (lufile3,'(i5)') nbfls(i,j)
            do k = 1, nbfls(i,j)
               write (lufile3,'(i4,2f16.8)')
     &               nls(i,j,k),zls(i,j,k),cls(i,j,k)
            enddo
         enddo
      enddo
c----------------------------------------------
c Card 9 : geo
c----------------------------------------------
 2900 continue
!     number sym-related atoms
      k = 1
      do i = 1, nsag
         do j = 1, RECP_NONT(i)
!           nc(i,j)= get_nucdeg(k)
            nc(i,j)= nucdeg(k)
            k = k + 1
         enddo
      enddo

      igcs = 0
      do i = 1, nsag
         do j = 1, RECP_NONT(i)
            mtype(i) = namnar(i,1)(1:3)
            nf(i) = 0
            do k = 1, iqm(i)
               nf(i) = nf(i) + jco(i,k)
            enddo 
            chg(i) = q(i) 
c           write(lufile3,'(a3,2i3,f3.0)') mtype(i),nf(i),nc(i,j),chg(i)
            write(lufile3,*) mtype(i),nf(i),nc(i,j),chg(i)

!           Geometry
            do j2=1,mj2(i,j)-1
               write(lufile3,'(3f14.10)') (xyz3(k,j2,j,i),k=1,3)
            enddo

c           (igen)
            igenq(1)=2
            igenq(2)=3
            igenq(3)=5
c           if (nc(i,j).gt.1) then
            if (mj2(i,j).gt.2) then
               do jgenp=1,ngenp           
                  write (lufile3,'(26i3)')
     &            (isend(k,j,igenq(jgenp),i),k=1,mj2(i,j)-1)
               enddo
            endif
c           (mcons, igcs)
            mcons = 0
            do k = 1, i-1
               do l = 1, iqm(k)
                  mcons = mcons + jco(k,l)
               enddo
            enddo
            do k = 1, iqm(i)
               igcs = igcs + 1
               do l = 1, jco(i,k) 
                  mcons = mcons + 1
                  write (lufile3,'(2i3)') mcons,igcs
               enddo
           enddo
c           (mcrs)
            if (mcrs(i).eq.0) then 
               do l = 1, i-1
                  if (mtype(i).eq.namnar(l,1)(1:3)) mcrs(i)=mcrs(l)
               enddo
            endif
            if (ncrs.gt.0) then
               write (lufile3,'(i3)') mcrs(i)
            endif
         enddo
      enddo
 
 3999 continue 
      close (lufile3)
      goto 6000 
c
c error message 
c
 5200 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'In subroutine argosin...'
      write(*,*) 'Cannot open MOLECULE.INP file'
      stop     
 5300 continue
      write(*,*) '*** ERROR ***'
      write(*,*) 'In subroutine argosin...'
      write(*,*) 'Cannot make argosin file'
      stop     
c
c end program   
c
 6000 continue

      IF (PRINT_LEVEL.GE.5)
     &   PRINT *,' * ECPINT : exiting...  ARGOS_GENINPUT'
      return
      end

      function get_nucdeg(k)
!radovan: this function was added because of overlapping header files
!         in calling subroutine
         integer :: get_nucdeg
         integer :: k
#include "mxcent.h"      
#include "nuclei.h"      

         get_nucdeg = nucdeg(k)
         return
      end function

      function get_cord(k, i)
!radovan: this function was added because of overlapping header files
!         in calling subroutine
         integer :: get_cord
         integer :: k
         integer :: i
#include "mxcent.h"      
#include "nuclei.h"      

         get_cord = cord(k, i)
         return
      end function



