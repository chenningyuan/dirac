      SUBROUTINE RECP_INP_MAIN(  a,     ccr,    chg,    eta,
     &   ica,    idp,    ityp,  la,     lb,     lcr,
     &   lls,    lmnp1,  lmn1u, lproju, maords, mcons,
     &   mcrs,   mgcs,   mnl,    ms,    mtype,  nblpr,  nc,
     &   ncon,   ncr,    ncru,   nct,   nd,     nf,     nir,
     &   nkcrl,  nkcru,  nklsl,  nklsu, nrcr,   nso,    nsopr,
     &   nt,     ntl,    ntu,    x,     y,      z,      zcr,
     &   zet,    N1MAX,  L1REC,CX_TEMP )
!  read and verify the user input: symmetry info, basis sets,
!  geometry, etc.
      USE RECP_INP_READ
      USE RECP_INP_ORBITAL 
      USE RECP_IPT
      USE RECP_NTR
      USE RECP_OUTPUT
      USE RECP_FUNCTION
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
      integer ierr
      integer
     & ica(mcu,msu,*),idp(mstu,mstu,*),la(mru,*),
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
     & iaords, iru, ir, igcs, icsu, ics, ic, icons, iconu,
     & ircru, icon, ircr, 
     & l, ngenp1, isf, isfr, is, icu, ig, if, ilmnp1, lai,
     & itl, itu, jcts, it, jct, jt, ixyz, i2, i3,
     & js, ngp1, iso, n1max, n2max, iau, ia, leig, lmn, 
     & lmnvmx,lmnvmx2,lmnvmx3 
      integer ixyzir(3), lblso(16),CX_TEMP
      real*8  t, eval
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
      CALL OUTPUT_DBG_TITLE('INP_MAIN : start process input data')

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
     &     nklsu,ncr,ncrs,ncru,lproju,ns,zcr,ccr)

      ngenp1 = ngen + 1
      ng     = ngenp1

!     read in atomic labels, charges, coordinates
      CALL RECP_INP2_READGEO(RECPIN1,ns,ng,mcu,msu,
     &     nf,nc,ica,X,Y,Z,chg,mtype,
     &     nst,nso,nbft,mcons,mgcs,nt,ntl,ntu,la,lb,nrcr,nir,lmn1u,
     &     lmnp1,maords,ncrs,mcrs)

      DEALLOCATE ( RECP_SET ) 
      CALL OUTPUT_DBG_TITLE('INP_MAIN : end process input data')
!     ----------------------------------------------------------------

!     # generate cartesian gaussian exponent array. (mklmnv)
      CALL RECP_INP_CARTEXP
 
!     # determine the number of types of one-electron integrals.
      CALL RECP_INP_INTCHECK(nu,ns,ncrs,mcrs,lcr,lls,inam) 
 
!     # normalization of symmetry orbitals. (deleted)
!     if (inrm.ne.0) CALL RECP_INP_ORB_SONORM
 
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

!     # normalize the contraction coefficients.
      CALL RECP_INP_ORBITAL_CCNORM
     &   ( mrcru,mconu,ncons,nrcr,ncon,lmnp1,eta,zet)

!     # calculate the nuclear repulsion energy.
!     CALL RECP_INP_NUCREP
!    &   (ic,nc,chg,NLIST,PRINT_LEVEL,x,y,z)

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

      CALL OUTPUT_DBG_TITLE('INP_MAIN : start print input info')
      CALL RECP_SETZERO_I1(ixyzir,3)

      do 840 is = 1,ns
        icu = nc(is)

        IF (PRINT_LEVEL.GE.2) THEN        
           CALL RECP_INP_INFO1(
     &          nlist,is,nc,mcu,msu,ica,ng,ngenp1,chg,x,y,z,mtype)
        ENDIF

        do 796 if = 1,nf(is)
          isf = isf+1
          itl = ntl(isf)
          itu = ntu(isf)
          igcs = mgcs(isf)
          iaords = maords(igcs)
          iru = nir(iaords)
          icons = mcons(isf)
          ircru = nrcr(icons)
          iconu = ncon(icons)
          icsu = 0
          do ir = 1,iru
            icsu = icsu + nd(la(ir,iaords))
          enddo

          IF (PRINT_LEVEL.GE.2) THEN        
             CALL RECP_INP_INFO2(
     &            lmnp1,icons,ircru,iconu,zet,eta,icsu,igcs,
     &            isfr,iru,la,nd,lb,ityp,icu,itl,itu,IPT_LMNV,iaords)
          ENDIF

          do 792 ircr = 1,ircru
            ics = 0
            do 788 ir = 1,iru
              isfr = isfr+1
              lai = la(ir,iaords)
              iso = nsopr(lai)+lb(isfr)
              ms(iso) = is
 
!             # compute the L^2-operator expectation value and
!             # determine if it is an eigenvalue.
              lmn = lmnp1(icons) - 1
              call l2opxv(lmn,itl,itu,IPT_LMNV,eval,leig,ics,igcs) 

!             # not an eigenvalue. print a warning.
              if ( leig .lt. 0 ) then
                IF (PRINT_LEVEL .GE. 2) THEN 
                   write(nlist,'(3x,a,i3,a3,a,a,f12.8)') 
     &               's.o. basis function',iso,ityp(lai),
     &               ' is not an L^2 eigenvector. : ',
     &               '<v|L^2|v>/<v|v> =', eval
                ENDIF 
                leig = lmn  !assign the maximum value
              endif

!             # assign the mnl(*) values to the so functions.
!             # 1:s, 2:p, 3:3s, 4:d, 5:4p, 6:f, 7:5s, 8:5d, 9:g,...
              mnl(iso) = ( lmnp1(icons)**2 + 2 * leig + 5 ) / 4

              iau = nd(lai)
              do ia = 1,iau
                 ics = ics+1
                 if ( (inam(nu).ge.5).and.(leig.eq.1) ) then
!                   # determine angular momentum symmetries.
                    CALL RECP_INP_ORBANGSYM
     &                   (icu,itl,itu,ics,igcs,lai,ixyzir)
                 endif
              enddo

  788       continue
  792     continue
  796   continue
 
        IF (PRINT_LEVEL.GE.2) THEN        
           CALL RECP_INP_PRINTCSO
     &        ( is,mcrs,lcr,nkcrl,nkcru,NLIST,ncr,nklsl,nklsu,lls,
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
 
!     # set up the 1-e record parameters.
      l1rec  = -1
      call sifcfg( nbft, ifmt1, l1rec, n1max, ierr )
      if ( ierr .ne. 0 ) then
        call bummer('syminp: from 1-e sifcfg(), ierr=',ierr,2)
      endif
 
!     # write out the header information
      call wthead(nst,ns,isfr,ityp,nso,mtype,ms,mnl,nlist)
 
!     # multiply normalization constants into contraction coefficients.
      do icons = 1,ncons
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

      CALL OUTPUT_LOC('RECP_INP_MAIN','X')
      END


