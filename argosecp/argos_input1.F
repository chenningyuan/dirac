*deck syminp_argos
      subroutine argos_input( a,      c,          ccr,    chg,    eta,
     &   ica,    idp,    ipt,    ityp,       la,     lb,     lcr,
     &   lls,    lmnp1,  lmnv,   lmn1u,      lproju, maords, mcons,
     &   mcrs,   mgcs,   mnl,    ms,         mtype,  nblpr,  nc,
     &   ncon,   ncr,    ncru,   nct,        nd,     nf,     nir,
     &   nkcrl,  nkcru,  nklsl,  nklsu,      nrcr,   nso,    nsopr,
     &   nt,     ntl,    ntu,    x,          y,      z,      zcr,
     &   zet,    lcore,  N1MAX,  L1REC )
 
!  read and verify the user input: symmetry info, basis sets,
!  geometry, etc.

      USE RECP_INPUT_READ
      USE RECP_ORBITAL 

      implicit logical(a-z)

#include "argos_mxvalue.h"
#include "argos_files.h" 
#include "argos_error.h" 
#include "argos_print.h" 

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
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
 
      integer
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir(3),      inam(5),     nnam, mdum(32)
 
      integer        ipq
      common /dim21/ ipq(1000) !256->1000
 
!     # dummy:
      integer lcore
      integer
     & ica(mcu,msu,*), idp(mstu,mstu,*), ipt(*),     la(mru,*),
     & lb(*),          lcr(*),           lls(*),     lmnp1(*),
     & lmnv(3,*),      lmn1u,            lproju,     maords(*),
     & mcons(*),       mcrs(*),          mgcs(*),    mnl(*),
     & ms(*),          nblpr(*),         nc(*),      ncon(*),
     & ncr(*),         ncru,             nct(*),     nd(*),
     & nf(*),          nir(*),           nkcrl(6,*), nkcru(6,*),
     & nklsl(4,*),     nklsu(4,*),       nrcr(*),    nso(*),
     & nsopr(*),       nt(*),            ntl(*),     ntu(*)
      real*8
     & a(*),               c(mctu,mcsu,*), ccr(*),   chg(*),
     & eta(mrcru,mconu,*), x(mcu,*),       y(mcu,*), z(mcu,*),
     & zcr(*),             zet(mconu,*)
      character*3 ityp(*), mtype(*)
c
c     # local:
      integer i, itol, inrm, ncrs, l1rec,
     & ngen, ncons, naords, ist, 
     & iaords, iru, ir, igcs, icsu, ictu, ics, ict, ic, icons, iconu,
     & ircru, icon, ircr, 
     & l, ngenp1, isf, isfr, is, icu, ig, if, ilmnp1, lai,
     & itl, itu, jcts, it, jct, jt, ixyz, i2, i3,
     & js, ngp1, iso, n1max, n2max, iau, ia, leig, lmn, 
     & N_LMNV !lmnvmx
      integer ixyzir(3), lblso(16)
      real*8 c_sum, t, eval

      character*1 ibl, ipc, isk, lblsh(21)

      INTEGER  ECP_INPUT
c     # for workspace allocation.
      integer  atebyt, forbyt
      external atebyt, forbyt

      integer  iodfac
      external iodfac

      data ibl, isk /' ','0'/, lblsh /'s','p','d','f','g','h','i','k',
     &  'l','m','n','o','q','r','t','u','v','w','x','y','z'/


!     # Start RECP_INPUT
      IF (PRINT_LEVEL.GE.5) 
     &   PRINT *,'* ECPINT : entering... ARGOS_INPUT'

!     assign the default unit numbers and filenames.
      ECP_INPUT = 5  ! RECP_INPUT :argos input file
      nlist  = 6  ! DIRAC output file 
!     intitalize some input variables.
      itol   = 20
      tol    = 2.30258d0 * itol   ! ln10 *itol
 
!     # Open the input file.
      OPEN(UNIT=ECP_INPUT,STATUS='OLD',FILE='RECP_INPUT')
      call ibummr(nlist)
      rewind ECP_INPUT
 
!     # read general information 
      CALL RECP_INPUT_READ_GENERAL
     &    (ECP_INPUT,ngen,ns,naords,ncons,ngcs,ncrs)

 
!     # Number of distinct products of irrep 
      call izero_argos( nst, nso, 1 )
      call izero_argos( mstu*mstu*nst, idp, 1 )
      CALL ECPINT_INPUT_READ_IRREP (ECP_INPUT,idp,nst,nd,ityp)
 
!     # Read AO to SO transformation matrices 
      CALL ECPINT_INPUT_READ_AO2SO
     &    (naords, ngcs, maords, la, nir, ECP_INPUT, c)
 

!     # read in exponents and contraction coefficients.
      CALL ECPINT_INPUT_READ_BASISSET
     &    (ECP_INPUT, lmnp1, nrcr, ncon, zet, eta, ncons)
 
!     # read in expansions for core and spin-orbit potentials.
      CALL ECPINT_INPUT_READ_RECP
     &   ( ncru, lproju, ncrs, ECP_INPUT, lcr,lls,nkcrl, nkcru,
     &     nklsl,nklsu, ncr, zcr, ccr )
 
      ngenp1 = ngen + 1
      ng     = ngenp1
 
!     # read in atomic labels, charges, coordinates, basis sets,
!       and core and spin-orbit potentials
      call izero_argos(ns,mcrs,1)
      CALL ECPINT_INPUT_READ_GEO
     &           (mcu, msu, ng, nbft, ECP_INPUT, ns, nf,
     &            nc, ica, mcons, mgcs, maords, lmnp1, lmn1u,
     &            ipq, nrcr, nir, nso, nt, ntl, ntu, la, lb, mcrs,
     &            chg, x, y, z, mtype, mru, ncrs)

      N_LMNV = (lmn1u*(lmn1u + 1)*(lmn1u + 2))/6
      ipt(2) = ipt(1) + forbyt( 3*N_LMNV )

!     # ipt(*)-->1:lmnv(1:3,1:N_LMNV)
      call mklmnv(lmn1u,lmnv)

!     # end of reading input file.
      CLOSE(UNIT=ECP_INPUT)

!     ---------------------------

 
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
!     if ( inrm .ne. 0 ) CALL ARGOS_ORBITAL_NORM
 
      call izero_argos(3,ixyzir,1)
 
!     # compute contracted orbital and lower-triangle-packed
!     # matrix symmetry offsets.
      nsopr(1) = 0
      nblpr(1) = 0
      do ist = 2,nst
         nsopr(ist) = nsopr(ist-1) + nso(ist-1)
         nblpr(ist) = nblpr(ist-1) + ipq( nso(ist-1) + 1 )
      enddo
 
c     # nnbft = total number of unique so function pairs.
c     # mpru = the number of elements in a diagonal,
c     #        symmetry-blocked, lower-triangle-packed matrix.  this
c     #        is also an upper bound to the size of any
c     #        symmetry-blocked array in which only unique
c     #        off-diagonal elements are stored and which is
c     #        associated with an operator that transforms as an
c     #        irreducible representation.
 
      nnbft = (nbft * (nbft + 1)) / 2
      mpru  = nblpr(nst) + ipq( nso(nst) + 1 )
 
!     # normalize the contraction coefficients.
      CALL RECP_ORBITAL_CONT_COEFF_NORM
     &   ( mrcru,mconu,ncons,nrcr,ncon,lmnp1,eta,zet)
!     # calculate the nuclear repulsion energy.
      CALL ARGOS_INPUT_CALC_NUCREP
     &   ( ic, nc, chg, NLIST, PRINT_LEVEL, x, y, z )


      call gcentr( ica, nc )
c
      ngp1 = ng+1
      if ( ngp1 .gt. mgu ) then
        call bummer('change mgup (two places) to ',ngp1,faterr)
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
     &             ISFR, IRU, LA, ND, LB, ITYP, C, ICU, ITL, ITU, 
     &             LMNV, IAORDS)
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
!             # ipt(*)-->1:lmnv(1:3,1:N_LMNV), 2:scr(1:xyzdim)
              call l2opxv( lmn, c(1,ics+1,igcs), itl, itu, lmnv,
     &         ipt(55)-ipt(2), a(ipt(2)), eval, leig )

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
                        if( mod(lmnv(ixyz,it),2) .eq. 1 .and.
     &                      mod(lmnv(i2,  it),2) .eq. 0 .and.
     &                      mod(lmnv(i3,  it),2) .eq. 0 ) then
                          c_sum = c_sum + c(ict,ics,igcs)
                        endif
                      enddo   
                    enddo   

                    if (  (c_sum .eq. 0.0d0)
     &               .or. (lai .eq. ixyzir(ixyz)) ) go to 780
                    if ( ixyzir(ixyz) .ne. 0 ) then
                      write (nlist,*) 'axis choice is unsuitable for'
     &                 //' spin-orbit integrals'
                      call bummer('syminp: ixyzir(ixyz)=',
     &                 ixyzir(ixyz),faterr)
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
      WRITE(NLIST,1000)

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
        call bummer('syminp: from 1-e sifcfg(), ierr=',ierr,faterr)
      endif
 
!     # write out the header information
!     # ipt(*)-->1:lmnv(1:3,1:N_LMNV), 2:map(1:2,1:nbft)
!     # push map(*) allocation onto the stack.
 
      call wthead( nst,    ns,     isfr,   ityp,
     &   nso,      mtype,  ms,     mnl,    nlist,
     &   a(ipt(2)) )
 
!     # multiply normalization constants into contraction coefficients.
      do icons = 1,ncons
         iconu = ncon(icons)
         ircru = nrcr(icons)
         do icon = 1,iconu
            do ircr = 1,ircru
               t = 4.0d0 * zet(icon,icons)
               eta(ircr,icon,icons) =
     &          sqrt( t**lmnp1(icons) / 2.0d0 * sqrt(t / 2.0d0) )
     &          * eta(ircr,icon,icons)
            enddo   
         enddo   
      enddo   

      IF (PRINT_LEVEL .GE. 2) THEN
         write (nlist,880) isfr, (ityp(ist), nso(ist), ist=1,nst)
  880    format(//i9,' symmetry orbitals,',3x,4(3x,a3,':',i4),
     &                                   (/31x,4(3x,a3,':',i4)))
      ENDIF
 
      IF (PRINT_LEVEL.GE.5) 
     &   PRINT *,' * ECPINT : exiting...  ARGOS_INPUT'

      return
 1000 FORMAT(/,10X,
     & '-----------------------------------------------------------')
      end


