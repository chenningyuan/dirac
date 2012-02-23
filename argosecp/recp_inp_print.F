      subroutine wthead( nsym,   ns,     nbft,
     &   ityp,   nbpsy,  mtype,  ms,     mnl,  nlist )
!  write the header info to the output file.
!  input:
!   nsym = number of symmetry blocks.
!   ns   = number of symmetry inequivalent atom labels.
!   nbft = total number of basis functions.
!   repnuc = nuclear repulsion energy.
!   ityp(1:nsym) = input symmetry labels.
!   nbpsy(1:nsym) = number of basis functions in each symmetry block.
!   mtype(1:ns) = atom labels.
!   ms(1:nbft) = basis_function-to-center mapping vector.
!   mnl(1:nbft) = integer codes for basis function types.
!                 1:s, 2:p, 3:3s, 4:d, 5:4p, 6:f, 7:5s, 8:5d, 9:g
!   nlist = listing file output unit.
!   map(1:2,1:nbft) = temporary scratch array for writing the output
!                     map(*) vectors.
      USE RECP_IPT
      implicit logical(a-z)
#include "inc_print.h"
      integer    nbfmx,     nmap,   nenrgy,   nrtype
      parameter( nbfmx=1000, nmap=2, nenrgy=1, nrtype=-1 )
!     parameter( nbfmx=500, nmap=2, nenrgy=1, nrtype=-1 )
      integer nsym, ns, nbft, nlist
!     real*8        repnuc
      character*3   ityp(1:nsym)
      integer       nbpsy(1:nsym)
      character*3   mtype(1:ns)
      integer       ms(1:nbft)
      integer       mnl(1:nbft)
!     # local:
      integer    bmin,   bmax
      parameter( bmin=0, bmax=17 )
      integer      clip, i, j
      character*8  bfnlab(nbfmx)
      character*4  slabel(8)
!     real*8       energy(nenrgy)
      integer      ietype(nenrgy)
      integer      imtype(nmap)
      character*2  btype(bmin:bmax)
      real*8   sifsce
      external sifsce
!     # map vector types.
      data imtype/ 3, 4 /
!     # basis function types.
      data btype/ '??', '_s', '_p', '3s', '_d', '4p', '_f', '5s', '5d',
     &      '_g', '6p', '6f', '_h', '7s', '7d', '7g', '_i', '??' /
 
      ALLOCATE(IPT_WTHEAD_MAP(nbft,2))

      if ( nbft .gt. nbfmx ) then
        WRITE (NLIST,1000) nbfmx,(nbft-nbfmx)
 1000   FORMAT(/
     &  1x,'In this version of ECP integral routine,'/
     &  1x,'the number of basis function are limited to',i4,/
     &  1x,i4,' basis functions are exceed')
        CALL QUIT('Number of basis function are exceed.')
      endif
 
!     # set map(*,*) and create bfnlab(*) from the input arrays.
 
      do 20 i = 1, nbft
        IPT_WTHEAD_MAP(i,1) = ms(i)
        IPT_WTHEAD_MAP(i,2) = mnl(i)
        clip = min( max( bmin, mnl(i) ), bmax )
        write( bfnlab(i), fmt='(i3,a3,a2)' )
     &   i, mtype(ms(i)), btype(clip)
        do 10 j = 4, 8
!         # replace embedded spaces with a printing character to
!         # avoid confusion when printing labels later.
          if ( bfnlab(i)(j:j) .eq. ' ' ) then
            bfnlab(i)(j:j) = '_'
          endif
10      continue
20    continue
 
!     # create the output symmetry labels.
      do i = 1, nsym
        slabel(i)(1:1) = ' '
        slabel(i)(2:4) = ityp(i)
      enddo

      IF (PRINT_LEVEL .GE. 2) THEN
 
!     # write the headers.

!     write(nlist,6100)'RECP integral information:'
!     write(nlist,6100)'output energy(*) values:'
!6100 format(/1x,a,(1x,10i10))
!     call sifpre( 6, nenrgy, energy, ietype )
!     write(nlist,6200)'total core energy :',
!    & sifsce( nenrgy, energy, ietype )
!6200 format(/1x,a,(1x,1p4e20.12))

      write(nlist,'(/,a,2x,i2)') 'number of symmetry operator = ',nsym
      write(nlist,6042) nbft
6042  format( 1x,'number of basis functions   = ',i4)
      write(nlist,6051) (i,i=1,nsym)
6051  format(/1x,'symmetry number  = ',8i5)
      write(nlist,6052) (slabel(i),i=1,nsym)
6052  format( 1x,'symmetry label   = ',1x,8(1x,a4))
      write(nlist,6053) (nbpsy(i),i=1,nsym)
6053  format( 1x,'basis function   = ',8i5)

!    # orbital labels
      write(nlist,'(/,a)') ' output orbital labels (i:bfnlab(i))'
      write(nlist,'(5(i4,a1,a8))') (i,':',bfnlab(i),i=1,nbft)

!    # basis function to atomic center map
      write(nlist,'(/,a)') ' basis function to center (i:map(i))'
      write(nlist,'(5(i4,a1,i4,4x))') 
     &            (i,':',IPT_WTHEAD_MAP(i,1),i=1,nbft)

!    # basis function to orbital type map
      write(nlist,'(/,a)') ' basis function to orbital type (i:map(i))'
      write(nlist,'(5(i4,a1,i4,4x))') 
     &            (i,':',IPT_WTHEAD_MAP(i,2),i=1,nbft)

      ENDIF

 
      DEALLOCATE(IPT_WTHEAD_MAP)
      return
      end


      SUBROUTINE RECP_INP_INFO1(
     &           nlist,is,nc,mcu,msu,ica,ng,ngenp1,chg,x,y,z,mtype)
      IMPLICIT NONE
      INTEGER     nlist,is,nc(*),mcu,msu,ica(mcu,msu,*),ng,ngenp1
      REAL*8      chg(*),x(mcu,*),y(mcu,*),z(mcu,*)
      CHARACTER*3 mtype(*)
!     local
      INTEGER I,J

      WRITE (nlist,'(/36X,A3,A)') mtype(is),' atoms'
      WRITE (nlist,'(/30X,A,F7.2)') 'nuclear charge',chg(is)
      WRITE (nlist,'(/11x,A,12X,A,15X,A,15X,A)') 'center','x','y','z'
      WRITE (nlist,'(I14,7X,3F16.8)')
     &        (I, x(I,is), y(I,is), z(I,is), I=1,nc(is))

      IF ( nc(is) .GT. 1 ) THEN
         WRITE (nlist,'(A,6X,A)') ' operator','center interchanges'
         WRITE (nlist,'(A,6X,24I3)') 
     &         '  1(id.)',(ica(I,is,1), I=1,nc(is))

         DO J = 2,ngenp1
            WRITE (nlist,'(I3,A,5X,24I3)')
     &            J,'(gen.)',(ica(I,is,J), I=1,nc(is))
         ENDDO

         DO J = ngenp1+1,ng
            WRITE(nlist,'(I3,11X,24I3)') J,(ica(I,is,J), I=1,nc(is))
         ENDDO
      ENDIF

      END SUBROUTINE RECP_INP_INFO1


      SUBROUTINE RECP_INP_INFO2(
     &           LMNP1, ICONS, IRCRU, ICONU, ZET, ETA, ICSU, IGCS,
     &           ISFR, IRU, LA, ND, LB, ITYP, ICU, ITL, ITU,
     &           LMNV, IAORDS)
      USE RECP_IPT
      IMPLICIT NONE

      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst

      integer
     & mconsu, mru,   mcru, msfu, 
     & msfru,  mnru,   ngcs,  nu,   
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & mconsu, mru,   mcru, msfu,
     & msfru,  mnru,   ngcs,  nu,   
     & lxyzir(3),      inam(5),     nnam, mdum(32)

      INTEGER     NLIST

      INTEGER     LMNP1(*), I
      REAL*8      ZET(MCONU,*), ETA(MRCRU,MCONU,*)
      INTEGER     ICONS, ICONU, ICON
      INTEGER     IRCRU, IRCRH, IRCR, IRCRL
      INTEGER     JSFR, ISFR, IR, IRU
      INTEGER     LAI, LA(MRU,*), IAORDS, IA, IAU, ND(*), LB(*)
      INTEGER     ICSR, ICOL, NCOL, LBLA(16), ICT, IC, IT
      INTEGER     ICU, ITL, ITU, LMNV(3,*), ICSU, IGCS
      CHARACTER*3 ITYP(*), LTYP(16)
      CHARACTER*1 LBLSH(21)
      DATA        LBLSH
     &             /'s','p','d','f','g','h','i','k','l','m','n',
     &              'o','q','r','t','u','v','w','x','y','z'/

      NLIST = 6

      write (nlist,'(/,10X,60A)') ('-', I=1,60)
!     # orbital exponents  contraction coefficients for each orbital
      write (nlist,708) LMNP1(ICONS), lblsh(lmnp1(icons))
  708 format(/i21,a1,' orbitals'//
     &       ' orbital exponents  contraction coefficients')
!     use larger value than 4 for output lines longer than 80
      ircrh = min(4,ircru)
      do 712 icon = 1,iconu
        write (nlist,'(8(1pg16.7))') zet(icon,icons),
     &    (eta(ircr,icon,icons), ircr=1,ircrh)
  712 continue

!     # symmetry orbital labels
      write (nlist,"(/21x,'symmetry orbital labels')")
      jsfr = isfr
      do 720 ir = 1,iru
        jsfr = jsfr + 1
        lai = la(ir,iaords)
        do 716 ia = 1,nd(lai)
          write (nlist,'(10x,7(i12,a3,i1))') (lb(jsfr+(ircr-1)*
     &      icsu),ityp(lai),ia, ircr=1,ircrh)
  716   continue
  720 continue

  724 if(ircrh.lt.ircru) then
        ircrl = ircrh + 1
        ircrh = min(ircrh+4,ircru)
        write (6,"(/20x,'contraction coefficients')")
        do 728 icon = 1,iconu
          write (nlist,'(16x,7(1pg16.7))')
     &      (eta(ircr,icon,icons), ircr=ircrl,ircrh)
  728   continue
        write (nlist,"(/21x,'symmetry orbital labels')")
        jsfr = isfr
        do 736 ir = 1,iru
          jsfr = jsfr + 1
          lai = la(ir,iaords)
          do 732 ia = 1,nd(lai)
            write (nlist,'(10x,7(i12,a3,i1))') (lb(jsfr+(ircr-1)*
     &        icsu),ityp(lai),ia, ircr=ircrl,ircrh)
  732     continue
  736   continue
        go to 724
      endif

      write (nlist,"(/11x,'symmetry orbitals')")
      icsr = 0
      ncol = 0
      do 760 ir = 1,iru
        lai = la(ir,iaords)
        iau = nd(lai)
        do 756 ia = 1,iau
          ncol = ncol+1
          ltyp(ncol) = ityp(lai)
          lbla(ncol) = ia
          if ( ncol.ne.10 .and.
     &      (ir.ne.iru .or. ia.ne.iau) ) go to 756
          write (nlist,740) (ltyp(icol),lbla(icol), icol=1,ncol)
  740     format(' ctr, ao',16(3x,a3,i1))
          ict = 0
          do 752 ic = 1,icu
            do 748 it = itl,itu
              ict = ict+1
              write (nlist,744) ic, (lmnv(i,it), i=1,3),
     &          (IPT_AO2SO(ict,icsr+icol,igcs), icol=1,ncol)
  744         format(i3,', ',3i1,16f7.3)
  748       continue
  752     continue
          icsr = icsr+ncol
          ncol = 0
  756   continue
  760 continue

      WRITE(NLIST,*)
      RETURN
      END

