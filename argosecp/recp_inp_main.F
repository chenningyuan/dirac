      SUBROUTINE RECP_INP_MAIN(  a,     ccr,    chg,    eta,
     &   ica,    idp,    ipt,    ityp,  la,     lb,     lcr,
     &   lls,    lmnp1,   lmn1u, lproju, maords, mcons,
     &   mcrs,   mgcs,   mnl,    ms,    mtype,  nblpr,  nc,
     &   ncon,   ncr,    ncru,   nct,   nd,     nf,     nir,
     &   nkcrl,  nkcru,  nklsl,  nklsu, nrcr,   nso,    nsopr,
     &   nt,     ntl,    ntu,    x,     y,      z,      zcr,
     &   zet,    lcore,  N1MAX,  L1REC,CX_TEMP )
 
!  read and verify the user input: symmetry info, basis sets,
!  geometry, etc.
      USE RECP_INP_READ
      USE RECP_INP_ORBITAL 
      USE RECP_IPT
      USE RECP_OUTPUT
      
      implicit logical(a-z)

#include "inc_mxvalue.h"
#include "inc_print.h" 

!     # /bufout/ holds some output integral file parameters.
      integer
     & ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
      common /bufout/
     & ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
      real*8         tol
      common /parmr/ tol
      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      integer
     & mconsu, mru,  mcru, msfu, 
     & msfru,  mnru, ngcs, nu,   
     & lxyzir, inam, nnam, mdum
      common /ntgr/
     & mconsu, mru,  mcru, msfu, 
     & msfru,  mnru, ngcs, nu,   
     & lxyzir(3),inam(5),nnam,mdum(32)
!     # dummy:
      integer lcore,ierr
      integer
     & ica(mcu,msu,*),idp(mstu,mstu,*),ipt(*),la(mru,*),
     & lb(*),lcr(*),lls(*),lmnp1(*),lmn1u,lproju,maords(*),
     & mcons(*),mcrs(*),mgcs(*),mnl(*),ms(*),nblpr(*),nc(*),ncon(*),
     & ncr(*),ncru,nct(*),nd(*),nf(*),nir(*),nkcrl(6,*),nkcru(6,*),
     & nklsl(4,*),nklsu(4,*),nrcr(*),nso(*),
     & nsopr(*),nt(*),ntl(*),ntu(*)
      real*8
     & a(*),ccr(*),chg(*),eta(mrcru,mconu,*),x(mcu,*),y(mcu,*),z(mcu,*),
     & zcr(*),zet(mconu,*)
      character*3 ityp(*), mtype(*)
 
!     # local:
      integer I,J,K, itol, inrm, ncrs, l1rec,
     & ngen, ncons, naords, ist, 
     & iaords, iru, ir, igcs, icsu, ictu, ics, ict, ic, icons, iconu,
     & ircru, icon, ircr, 
     & l, ngenp1, isf, isfr, is, icu, ig, if, ilmnp1, lai,
     & itl, itu, jcts, it, jct, jt, ixyz, i2, i3,
     & js, ngp1, iso, n1max, n2max, iau, ia, leig, lmn, 
     & lmnvmx,lmnvmx2,lmnvmx3 
      integer ixyzir(3), lblso(16),CX_TEMP
      real*8 c_sum, t, eval
      character*1 ibl, ipc, isk, lblsh(21)
      INTEGER  ECP_INPUT,RECPIN1,nlist
      integer  iodfac
      external iodfac
      LOGICAL  OLDINP
      data ibl, isk /' ','0'/, lblsh /'s','p','d','f','g','h','i','k',
     &  'l','m','n','o','q','r','t','u','v','w','x','y','z'/

      CALL OUTPUT_LOC('RECP_INP_MAIN','E')

!     assign the default unit numbers and filenames.
      ECP_INPUT = 5  ! RECP_INPUT :argos input file
      RECPIN1   = 51 ! RECP_INPUT :argos input file
      nlist     = 6  ! DIRAC output file 
!     intitalize some input variables.
      itol   = 20
      tol    = 2.30258d0 * itol   ! ln10 *itol

!     ----------------------------------------------------------------
      OPEN(UNIT=RECPIN1,STATUS='OLD',FILE='RECP_INPUT2')

!     # read general information 
      CALL RECP_INP2_READGEN(RECPIN1,ngen,ns,naords,ncons,ngcs,ncrs)

!     # Number of distinct products of irrep 
      CALL RECP_INP2_READIRREP(RECPIN1,idp,nd,ityp)

!     # Read AO to SO transformation matrices 
      CALL RECP_INP2_READAO2SO
     &    (RECPIN1,naords,ngcs,maords,la,nir,CX_TEMP)

!     # read in exponents and contraction coefficients.
      CALL RECP_INP2_READORB(RECPIN1,zet,eta,lmnp1,nrcr,ncon)

!     read in expansions for core and spin-orbit potentials.
      CALL RECP_INP2_READCSO(RECPIN1,lcr,lls,nkcrl,nkcru,nklsl,
     &     nklsu,ncr,ncrs,ncru,lproju,ns,mcrs,zcr,ccr)

      ngenp1 = ngen + 1
      ng     = ngenp1

!     read in atomic labels, charges, coordinates
      CALL RECP_INP2_READGEO(RECPIN1,ns,ng,nf,nc,ica,X,Y,Z,chg,mtype,
     &     nst,nso,nbft,mcons,mgcs,nt,ntl,ntu,la,lb,nrcr,nir,lmn1u,
     &     lmnp1,maords,ncrs,mcrs)

      CLOSE(UNIT=RECPIN1)

!     ----------------------------------------------------------------

!     FIXME : lmnvmx value
      lmnvmx = (lmn1u*(lmn1u + 1)*(lmn1u + 2))/6
      lmnvmx2=  84
      lmnvmx3=  max(lmnvmx*3,84)
      ipt(2) = ipt(1) 
      ALLOCATE(IPT_LMNV(3,lmnvmx3))

      DO I = 1, 3
         DO J = 1, lmnvmx3
!        DO J = 1, 84
            IPT_LMNV(I,J) = 0
         ENDDO
      ENDDO
!     # ipt(*)-->1:lmnv(1:3,1:lmnvmx)
      call mklmnv(lmn1u)

!     ----------------------------------------------------------------
 
!     # determine the number of types of one-electron integrals.
 
!     # initialize nu for s(*),t(*),v(*).
      nu = 3
 
      if ( ncrs .gt. 0 ) then
        do 196 is = 1,ns
          if ( mcrs(is) .ne. 0 ) then
            if ( lcr(mcrs(is)) .gt. 0 ) then
              nu = 4
              inam(4) = 4
              go to 200
            endif
          endif
  196   continue
c
  200   continue
        do 204 is = 1, ns
          if ( mcrs(is) .ne. 0 ) then
            if ( lls(mcrs(is)) .gt. 0 ) then
              nu = nu + 1
              inam(nu) = 5
              go to 208
            endif
          endif
  204   continue
c
  208   continue
      endif
 
!     # normalization of symmetry orbitals. (deleted)
!     if ( inrm .ne. 0 ) CALL RECP_INP_ORB_SONORM
 
      DO I = 1, 3
         ixyzir(I) = 0
      ENDDO
 
!     contracted orbital and lower-triangle-packed matrix symmetry
!     offsets.
      CALL RECP_INP_ORBITAL_OFFSET(nsopr,nso,nblpr,nst)
 
!     # nnbft = total number of unique so function pairs.
!     # mpru = the number of elements in a diagonal,
!     #        symmetry-blocked, lower-triangle-packed matrix.  this
!     #        is also an upper bound to the size of any
!     #        symmetry-blocked array in which only unique
!     #        off-diagonal elements are stored and which is
!     #        associated with an operator that transforms as an
!     #        irreducible representation.
 
      nnbft = (nbft * (nbft + 1)) / 2
      mpru  = nblpr(nst) + (nso(nst)+1)*((nso(nst)+1)-1)/2
!     print *,'bbnft',nnbft 
!     print *,'mpru ',mpru  

!     # normalize the contraction coefficients.
      CALL RECP_INP_ORBITAL_CCNORM
     &   ( mrcru,mconu,ncons,nrcr,ncon,lmnp1,eta,zet)
!     DO I = 1, ncons
!        print *,'From RECP_INP(CCNORM)'
!        DO J = 1, ncon(I)
!           print *, zet(J,I),(eta(K,J,I),K=1,nrcr(I))
!        ENDDO
!     ENDDO

!     # calculate the nuclear repulsion energy.
!     CALL RECP_INP_NUCREP
!    &   ( ic,nc,chg,NLIST,PRINT_LEVEL,x,y,z )

      call gcentr( ica, nc )
 
      ngp1 = ng+1
      if ( ngp1 .gt. mgup ) then
        call bummer('change mgup (two places) to ',ngp1,2)
      endif

      IF (PRINT_LEVEL .GE. 2) THEN
         write (nlist,'(a,a,i2,a)')
     &         ' * primitive ao integrals neglected',
     &         'if exponential factor below 10**(-',itol,')'
      ENDIF 

      isf = 0
      isfr = 0

      do 840 is = 1,ns
        icu = nc(is)

        IF (PRINT_LEVEL .GE. 2) THEN        
           write(nlist,656) mtype(is), chg(is),
     &                     (ic, x(ic,is), y(ic,is), z(ic,is), ic=1,icu)
  656      format(//36x,a3,' atoms'
     &            //30x,'nuclear charge',f7.2
     &            //11x,'center',12x,'x',15x,'y',15x,'z'
     &             /(i14,7x,3f16.8))
           if ( icu .gt. 1 ) then
              write (nlist,660) (ica(ic,is,1), ic=1,icu)
  660         format(' operator',6x,'center interchanges'/
     &               '  1(id.)',6x,24i3)

              do ig = 2,ngenp1
                 write (nlist,664) ig, (ica(ic,is,ig), ic=1,icu)
  664            format(i3,'(gen.)',5x,24i3)
              enddo   

              do ig = ngenp1+1,ng
                 write(nlist,'(i3,11x,24i3)') 
     &                ig,(ica(ic,is,ig), ic=1,icu)
              enddo   
           endif
        ENDIF

        do 796 if = 1,nf(is)
          isf = isf+1
          itl = ntl(isf)
          itu = ntu(isf)
          igcs = mgcs(isf)
          iaords = maords(igcs)
          iru = nir(iaords)
!         ictu = nct(isf)
          icons = mcons(isf)
          ircru = nrcr(icons)
          iconu = ncon(icons)
          icsu = 0
          do ir = 1,iru
            icsu = icsu + nd(la(ir,iaords))
          enddo

          IF (PRINT_LEVEL .GE. 2) THEN        
              CALL ARGOS_INPUT_INFO2(
     &             LMNP1, ICONS, IRCRU, ICONU, ZET, ETA, ICSU, IGCS,
     &             ISFR, IRU, LA, ND, LB, ITYP, ICU, ITL, ITU, 
     &             IPT_LMNV, IAORDS)
          ENDIF

          do 792 ircr = 1,ircru
            ics = 0
            do 788 ir = 1,iru
              isfr = isfr+1
              lai = la(ir,iaords)
              iso = nsopr(lai)+lb(isfr)
              ms(iso) = is
 
!             # assign the mnl(*) values to the so functions.
!             # 1:s, 2:p, 3:3s, 4:d, 5:4p, 6:f, 7:5s, 8:5d, 9:g,...
!             #
!             # compute the L^2-operator expectation value and
!             # determine if it is an eigenvalue.
 
              lmn = lmnp1(icons) - 1
!             numxyz = ((lmn + 1) * (lmn + 2)) / 2
!             # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:scr(1:xyzdim)
 
!             call l2opxv(lmn,IPT_AO2SO(1,ics+1,igcs), itl, itu, 
              call l2opxv(lmn, itl, itu, 
     &                  IPT_LMNV,  eval,leig,ics,igcs )

              if ( leig .lt. 0 ) then
                IF (PRINT_LEVEL .GE. 2) THEN 
!                  # not an eigenvalue. print a warning.
                   write(nlist,'(a,i3,a3,a)') 
     &               '   s.o. basis function',iso,ityp(lai),
     &               ' is not an L^2 eigenvector.'
                   write(nlist,'(a,f12.8)') '   <v|L^2|v>/<v|v> =', eval
                ENDIF 
!               # assign the maximum value.
                leig = lmn
              endif

              mnl(iso) = ( lmnp1(icons)**2 + 2 * leig + 5 ) / 4
              iau = nd(lai)
              do 784 ia = 1,iau
                ics = ics+1
 
                if ( inam(nu) .ge. 5 .and. leig .eq. 1 ) then
 
!                 # determine angular momentum symmetries.
                  do 780 ixyz = 1, 3
                    i2 = mod(ixyz,  3) + 1
                    i3 = mod(ixyz+1,3) + 1
                    c_sum = 0.0d0
                    ict = 0
                    do ic = 1, icu
                      do it = itl,itu
                        ict = ict + 1
                        if( mod(IPT_LMNV(ixyz,it),2) .eq. 1 .and.
     &                      mod(IPT_LMNV(i2,it)  ,2) .eq. 0 .and.
     &                      mod(IPT_LMNV(i3,it)  ,2) .eq. 0 ) then
                          c_sum = c_sum + IPT_AO2SO(ict,ics,igcs)
                        endif
                      enddo   
                    enddo   

                    if (  (c_sum .eq. 0.0d0)
     &               .or. (lai .eq. ixyzir(ixyz)) ) go to 780
                    if ( ixyzir(ixyz) .ne. 0 ) then
                      write (nlist,*) 'axis choice is unsuitable for'
     &                              //' spin-orbit integrals'
                      call bummer('syminp: ixyzir(ixyz)=',
     &                 ixyzir(ixyz),2)
                    endif
                    ixyzir(ixyz) = lai
  780             continue
                endif
  784         continue
  788       continue
  792     continue
  796   continue

c
        IF (PRINT_LEVEL .GE. 2) THEN        
           CALL ARGOS_INPUT_PRINT_RECP
     &        ( is, mcrs, lcr,nkcrl, nkcru,
     &          NLIST, ncr,nklsl, nklsu, lls,
     &          zcr, ccr, mtype, lblsh )
        ENDIF

  840 continue

!     -------------------------------

      if ( inam(nu) .eq. 5 ) then
        do 860 ixyz = 1, 3
          i2 = mod(ixyz,  3) + 1
          i3 = mod(ixyz+1,3) + 1
          do ist = 1, nst
            if ( idp(ist,ixyzir(i2),ixyzir(i3)) .ne. 0 ) then
              lxyzir(ixyz) = ist
              goto 860
            endif
          enddo
  860   continue

        IF (PRINT_LEVEL.GE.5) THEN
           write (nlist,864) (ityp(lxyzir(ixyz)), ixyz=1,3)
  864      format(/'lx: ',a3,10x,'ly: ',a3,10x,'lz: ',a3)
        ENDIF
      endif
 
c     # set up the 1-e record parameters.
      l1rec  = -1
      call sifcfg( nbft, ifmt1, l1rec, n1max, ierr )
      IF (PRINT_LEVEL.GE.5) THEN
         PRINT *, '   L1RECX :',L1RECX
         PRINT *, '   N1MAX  :',N1MAX
      ENDIF
 
      if ( ierr .ne. 0 ) then
        call bummer('syminp: from 1-e sifcfg(), ierr=',ierr,2)
      endif
 
!     # write out the header information
!     # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:map(1:2,1:nbft)
!     # push map(*) allocation onto the stack.
 
      call wthead( nst,    ns,     isfr,   ityp,
     &   nso,      mtype,  ms,     mnl,    nlist)
 
!     # multiply normalization constants into contraction coefficients.
      do icons = 1,ncons
!        print *,'From RECP_INP'
         do icon = 1,ncon(icons)
            do ircr = 1,nrcr(icons)
               t = 4.0d0 * zet(icon,icons)
               eta(ircr,icon,icons) =
     &          sqrt( t**lmnp1(icons) / 2.0d0 * sqrt(t / 2.0d0) )
     &          * eta(ircr,icon,icons)
            enddo   
!           print *,zet(icon,icons),
!    &              (eta(ircr,icon,icons),ircr=1,nrcr(icons))
         enddo   
      enddo   

      IF (PRINT_LEVEL.GE.2) THEN
         write (nlist,880) isfr, (ityp(ist), nso(ist), ist=1,nst)
  880    format(//i9,' symmetry orbitals,',3x,4(3x,a3,':',i4),
     &                                   (/31x,4(3x,a3,':',i4)))
      ENDIF

 1000 FORMAT(/,10X,
     & '-----------------------------------------------------------')

      CALL OUTPUT_LOC('RECP_INP_MAIN','X')
      END


