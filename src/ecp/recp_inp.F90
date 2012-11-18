!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

MODULE RECP_INP
CONTAINS
SUBROUTINE RECP_INP_MAIN(a,ccr,chg,eta,ica,idp,ityp,la,lb,lcr,  &
           lls,lmnp1,lmn1u,lproju,maords,mcons,mcrs,mgcs,mnl,ms,mtype,nblpr,nc, &
           ncon,ncr,ncru,nct,nd,nf,nir,nkcrl,nkcru,nklsl,nklsu,nrcr,nso,nsopr,  &
           nt,ntl,ntu,x,y,z,zcr,zet,CX_TEMP)
! read and verify the user input: symmetry info, basis sets,
! geometry, etc.
  USE RECP_INP_READ
  USE RECP_INP_ORBITAL 
  USE RECP_INP_PRINT
  USE RECP_IPT
  USE RECP_NTR
  USE RECP_OUTPUT
  USE RECP_FUNCTION1
  USE RECP_FUNCTION2
  USE RECP_WRITE
  implicit logical(a-z)
#include "inc_mxvalue.h"
#include "inc_print.h" 
! /bufout/ holds some output integral file parameters.
  integer         itypea,itypeb,ibuf,numout,nrec,ntape
  common /bufout/ itypea,itypeb,ibuf,numout,nrec,ntape
  real*8         tol
  common /parmr/ tol
  integer        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /parmi/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  integer        mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir,   inam,   nnam,mdum
  common /ntgr/  mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir(3),inam(5),nnam,mdum(32)
  integer ierr
  integer ica(mcu,msu,*),idp(mstu,mstu,*),la(mru,*),lb(*),lcr(*),lls(*),  &
          lmnp1(*),lmn1u,lproju,maords(*),mcons(*),mcrs(*),mgcs(*),mnl(*),ms(*),   &
          nblpr(*),nc(*),ncon(*),ncr(*),ncru,nct(*),nd(*),nf(*),nir(*),nkcrl(6,*), &
          nkcru(6,*),nklsl(4,*),nklsu(4,*),nrcr(*),nso(*),nsopr(*),nt(*),ntl(*),ntu(*)
  real*8  a(*),ccr(*),chg(*),eta(mrcru,mconu,*),x(mcu,*),y(mcu,*),z(mcu,*),zcr(*),zet(mconu,*)
  character*3 ityp(*), mtype(*)
! local variables
  integer I,J,K,itol,inrm,ncrs,ngen,ncons,naords,ist,iaords,iru,ir,igcs,icsu, &
   ics,ic,icons,iconu,ircru,icon,ircr,l,isf,isfr,is,icu,ig,if,ilmnp1,lai,itl,itu,   &
   jcts,it,jct,jt,ixyz,i2,i3,js,iso,iau,ia,leig,lmn,lmnvmx,lmnvmx2,lmnvmx3 
  integer ixyzir(3), lblso(16),CX_TEMP
  real*8  t, eval
  character*1 ibl, ipc, isk, lblsh(21)
  INTEGER  FILE_NUM(2)
  CHARACTER(10) FILE_NAME(2)

  DATA FILE_NUM  / 64, 65 /
  DATA FILE_NAME /'RECP_INP_0','RECP_INP_P'/
  DATA ibl, isk /' ','0'/ 
  DATA lblsh /'s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u','v','w','x','y','z'/

  CALL OUTPUT_LOC('RECP_INP_MAIN','E')
! write file in debug mode
  IF (RECP_DBG.GE.1) THEN
     CALL RECP_WRITE_FILEOPEN(FILE_NUM(1),FILE_NAME(1),'REPLACE')
  ENDIF

! intitalize some input variables.
  itol   = 20
  tol    = 2.30258d0 * itol   ! ln10 *itol

! ----------------------------------------------------------------
  CALL OUTPUT_DBG_TITLE('INP_MAIN : start process input data')

! # read general information 
  CALL RECP_INP2_READGEN(ngen,ns,naords,ncons,ngcs,ncrs)

! # Number of distinct products of irrep 
  CALL RECP_INP2_READIRREP(idp,nd,ityp)

! # Read AO to SO transformation matrices 
  CALL RECP_INP2_READAO2SO(naords,ngcs,maords,la,nir,CX_TEMP)

! # read in exponents and contraction coefficients.
  CALL RECP_INP2_READORB(zet,eta,lmnp1,nrcr,ncon)

! read in expansions for core and spin-orbit potentials.
  CALL RECP_INP2_READCSO(lcr,lls,nkcrl,nkcru,nklsl,nklsu,ncr,ncrs,ncru,lproju,ns,zcr,ccr)

  ng = ngen + 1

! read in atomic labels, charges, coordinates
  CALL RECP_INP2_READGEO(ns,ng,mcu,msu,nf,nc,ica,X,Y,Z,chg,mtype,nst,nso,nbft, &
       mcons,mgcs,nt,ntl,ntu,la,lb,nrcr,nir,lmn1u,lmnp1,maords,ncrs,mcrs)

  DEALLOCATE ( RECP_SET ) 
  CALL OUTPUT_DBG_TITLE('INP_MAIN : end process input data')
! ----------------------------------------------------------------

! # generate cartesian gaussian exponent array. (mklmnv)
  CALL RECP_INP_CARTEXP
 
! # determine the number of types of one-electron integrals.
  CALL RECP_INP_INTCHECK(nu,ns,ncrs,mcrs,lcr,lls,inam) 
 
! # normalization of symmetry orbitals. (deleted)
! if (inrm.ne.0) CALL RECP_INP_ORB_SONORM
 
! contracted orbital and lower-triangle-packed matrix symmetry
! offsets.
  CALL RECP_INP_ORBITAL_OFFSET(nsopr,nso,nblpr,nst,FILE_NUM)

! nbft  = total number of basis functions. 
! nnbft = total number of unique so function pairs.
! mpru  = the number of elements in a diagonal,
!         symmetry-blocked, lower-triangle-packed matrix.  this
!         is also an upper bound to the size of any
!         symmetry-blocked array in which only unique
!         off-diagonal elements are stored and which is
!         associated with an operator that transforms as an
!         irreducible representation.
  nnbft = (nbft * (nbft + 1)) / 2
  mpru  = nblpr(nst) + (nso(nst)+1)*((nso(nst)+1)-1)/2

! # normalize the contraction coefficients.
  CALL RECP_INP_ORBITAL_CCNORM(mrcru,mconu,ncons,nrcr,ncon,lmnp1,eta,zet)

! # calculate the nuclear repulsion energy. (deleted)
! CALL RECP_INP_NUCREP(ic,nc,chg,x,y,z)

  call gcentr(ica,nc)
 
  IF ((ng+1).GT.mgup) call bummer('change mgup (two places) to ',(ng+1),2)
  IF (PRINT_LEVEL.GE.2) THEN
     WRITE (RECP_OUT,'(a,a,i2,a)')' * primitive ao integrals neglected', &
           'if exponential factor below 10**(-',itol,')'
  ENDIF 

  CALL OUTPUT_DBG_TITLE('INP_MAIN : start print input info')
  CALL RECP_SETZERO_I1(ixyzir,3)

  isf = 0
  isfr = 0
  do 840 is = 1,ns
     icu = nc(is)

     CALL RECP_INP_INFO1(is,nc,mcu,msu,ica,ng,ngen,chg,x,y,z,mtype)

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

        CALL RECP_INP_INFO2(lmnp1,icons,ircru,iconu,zet,eta,icsu,igcs, &
             isfr,iru,la,nd,lb,ityp,icu,itl,itu,IPT_LMNV,iaords)

        DO ircr = 1,ircru
           ics = 0
           DO ir = 1,iru
              isfr = isfr+1
              lai = la(ir,iaords)
              iso = nsopr(lai)+lb(isfr)
              ms(iso) = is
 
!             # compute the L^2-operator expectation value and
!             # determine if it is an eigenvalue.
              lmn = lmnp1(icons) - 1
              CALL RECP_INP_ORBITAL_L2OPXV(lmn,itl,itu,IPT_LMNV,eval,leig,ics,igcs) 

!             # not an eigenvalue. print a warning.
              IF (leig.LT.0) THEN
                 IF (PRINT_LEVEL.GE.2) THEN 
                    WRITE(RECP_OUT,'(3x,a,i3,a3,a,f12.8)')'s.o. basis function',iso,ityp(lai), &
                    ' is not an L^2 eigenvector. : <v|L^2|v>/<v|v> =', eval
                 ENDIF 
                 leig = lmn  !assign the maximum value
              ENDIF

!             # assign the mnl(*) values to the so functions.
!             # 1:s, 2:p, 3:3s, 4:d, 5:4p, 6:f, 7:5s, 8:5d, 9:g,...
              mnl(iso) = ( lmnp1(icons)**2 + 2 * leig + 5 ) / 4

              iau = nd(lai)
              DO ia = 1,iau
                 ics = ics+1
                 IF ( (inam(nu).ge.5).and.(leig.eq.1) ) then
!                   # determine angular momentum symmetries.
                    CALL RECP_INP_ORBITAL_ANGSYM(icu,itl,itu,ics,igcs,lai,ixyzir)
                 ENDIF
              ENDDO
           ENDDO    
        ENDDO
  796   continue
        CALL RECP_INP_PRINTCSO(is,mcrs,lcr,nkcrl,nkcru,ncr,nklsl,nklsu,lls,zcr,ccr,mtype,lblsh)
  840 continue

! # multiply normalization constants into contraction coefficients.
  CALL RECP_INP_ORBITAL_MULTNORMCOEF(ncons,ncon,nrcr,mrcru,mconu,zet,eta,lmnp1)

!     -------------------------------

  if (inam(nu).eq.5) then
     do ixyz = 1, 3
        i2 = mod(ixyz,  3) + 1
        i3 = mod(ixyz+1,3) + 1
        do ist = 1, nst
          if ( idp(ist,ixyzir(i2),ixyzir(i3)) .ne. 0 ) then
            lxyzir(ixyz) = ist
            goto 860
          endif
        enddo
        860 continue
     enddo

     IF (PRINT_LEVEL.GE.5) THEN
        write (RECP_OUT,864) (ityp(lxyzir(ixyz)), ixyz=1,3)
        864 format(/'lx: ',a3,10x,'ly: ',a3,10x,'lz: ',a3)
     ENDIF
  endif
 
! # write out the header information
  CALL WTHEAD(nst,ns,isfr,ityp,nso,mtype,ms,mnl)
  IF (PRINT_LEVEL.GE.2) THEN
     WRITE(RECP_OUT,880) isfr, (ityp(I), nso(I), I=1,nst)
     880 FORMAT(//i9,' symmetry orbitals,',3x,4(3x,a3,':',i4),(/31x,4(3x,a3,':',i4)))
  ENDIF
 
  IF (RECP_DBG.GE.1) CLOSE (74) !close file (RECP_DBG_INP)
  CALL OUTPUT_LOC('RECP_INP_MAIN','X')
END SUBROUTINE RECP_INP_MAIN
END MODULE RECP_INP
