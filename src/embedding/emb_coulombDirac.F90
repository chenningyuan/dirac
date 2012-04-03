!dirac_copyright_start
!      copyright (c) 2011 by the authors of dirac.
!      all rights reserved.
!
!      this source code is part of the dirac program package.
!      it is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      in particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      this means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      for information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      dirac program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

! file    : emb1int.f
!
!
!  /* deck oneprm */
      module vc_grid

      implicit none

      private

      public emb_vc
      public get_esp
      public esp_nuc
      public esp_el
      
!     module for calculating the coulomb potential on a grid.
!     contains:
!     emb_vc    the interface for users. 
!     loop_prim loops over primitive orbitals.
!     vc_1prim  the routine that calls for one primitive shell vc and vc_cart
!               to calculate the shell's contribution to the coulomb potential.
!     vc        calculates the coulomb integrals in hermite basis
!     vc_cart   transforms the coulomb integrals to a cartesian basis.
!     
!
!     calls the dirac f77 routines: getodc in abacus/her1odc.f
!     based on the routines: hernai, cint0 in abacus/her1car.f
!     the loop over the primitives is based on abacus/her1drv.f and abacus/her1sym.f
!
!     the include files that need to be set:
!     onecom.h
!     sdepre.h
!

      contains

      subroutine get_esp(esp,              &
                         irep,             &
                         mat_dim,          &
                         dmat,             &
                         px,               &
                         py,               &
                         pz,               &
                         include_nuc_part, &
                         include_el_part)

!     --------------------------------------------------------------------------
      real(8), intent(out)          :: esp
      integer, intent(in)           :: irep
      integer, intent(in)           :: mat_dim
      real(8), intent(in)           :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)           :: px
      real(8), intent(in)           :: py
      real(8), intent(in)           :: pz
      logical, intent(in), optional :: include_nuc_part
      logical, intent(in), optional :: include_el_part
!     --------------------------------------------------------------------------
      logical                       :: nuc_part
      logical                       :: el_part
      integer, parameter            :: lwork = 10000
      real(8)                       :: r_p(3,1), work(lwork), vc(1)
!     --------------------------------------------------------------------------

      r_p(1, 1) = px
      r_p(2, 1) = py
      r_p(3, 1) = pz

      nuc_part = .true.
      if (present(include_nuc_part)) nuc_part = include_nuc_part
      el_part = .true.
      if (present(include_el_part)) el_part = include_el_part

      vc = 0.0d0
      if (nuc_part) then
         call esp_nuc(1, r_p, vc, 0)
      end if
      if (el_part) then
         call esp_el(1, dmat, r_p, vc, work, lwork, 0, .false.)
      end if

      esp = vc(1)

      end subroutine

      subroutine emb_vc(nr_points, dmat, points, vc, iprint, work, lwork)
!     interface for the user. calles loop_prim
!     calls vc_prim.
!     needed to make contact with f77.
!     jhs dmat is in mo basis. transformations needed.

#include "priunit.h"
#include "maxaqn.h"
#include "ccom.h"

!     passed variables
      integer, intent(in) :: nr_points
      real(8) :: points(3,*),vc(*)
      integer :: iprint
      real(8) :: dmat(*)
      real(8) :: work(*)
      integer :: lwork
      
      if (iprint .ge. 5) call titler('output from emb_vc','*',103)

!     electrostatic potential: electronic contribution
      call esp_el(nr_points, dmat, points, vc, work, & 
                  lwork, iprint, .false.)

!     electrostatic potential: nuclear contribution
      call esp_nuc(nr_points, points,vc,iprint)

      end subroutine

      subroutine esp_nuc(nr_points, points, vc, iprint)
!      
! Routine for the calculation of the nuclear part of the molecular 
! electrostatic potential on a certain grid of points {r_i}:
!     V_nuc(r_i) = sum_K Z_K/|R_K - r_i| 
! 
!  array vc(r_i) contains the nuclear MEP at point r_i.
! 
! We also take care of symmetry, since it is straightforward 
! for the nuclei.      
! This routine is essential for the Polarizable Continuum Model.
!
! Originally written by JHS. Cleaned by RDR 0312.
!

#include "mxcent.h"
#include "nuclei.h"

!     Passed variables
      integer, intent(in) :: nr_points
      real(8) :: points(3, *), vc(nr_points)
      integer :: iprint

!     Local variables
      real(8) :: coora(3, nucdep), charges(nucdep), dist
      integer :: ipoint, i, j, k

!     Get coordinates and charges of all centers, not only the ones
!     that are symmetry independent      
      call getacord(coora)
      i = 0
      do j=1,nucind
         do k = 1,nucdeg(j)
            i = i + 1
            charges(i) = charge(j)
         enddo
      enddo

      vc = 0.0d0
      do i = 1, nucdep 
        do ipoint = 1, nr_points
            dist = (coora(1, i) - points(1, ipoint))**2 +               &
                   (coora(2, i) - points(2, ipoint))**2 +               &
                   (coora(3, i) - points(3, ipoint))**2
            vc(ipoint) = vc(ipoint) + charges(i) / (sqrt(dist))       
        end do
      end do

      if (iprint > 5) then
          do ipoint = 1, nr_points
            print *, "Nuclear ESP @point", ipoint, vc(ipoint)
          end do
      end if
                                                                 
      end subroutine

      subroutine esp_el(nr_points, fock_or_dmat, points, pot_or_chg,    &
                        work, lwork, iprint, get_fock)
!      
! Driver routine for the calculation of the electronic part of the molecular 
! electrostatic potential on a certain grid of points {r_i}:
!     V_el(r_i) = tr(DV_i)
! tr is the trace operator, D is the density matrix, V^i is the matrix of the 
! "nuclear attraction" integrals calculated at point r_i of the grid:
!     V_mu,nu,i =  - <mu|1/|r-r_i||nu>
! 
!  array ADER(mu, nu, r_i) contains these integrals.
!  array vc(r_i) contains the electronic MEP at point r_i.
!
! This routine is essential for the Polarizable Continuum Model.
!
! Originally written by JHS. Debugged by RB and RDR 0312.
! 
! RDR 050312 CANNOT yet handle symmetry. This is unimportant at the moment
!            since also PCMSolver cannot handle symmetry.
! RDR 220312 This routine will be used to form both potentials and Fock
!            matrix contribution for PCM.
!            fock_or_dmat is Fock or density matrix, pot_or_chg is potentials 
!            or charges vector.
!            get_fock logical is present and TRUE:
!            charges vector as input, Fock matrix contribution as output.
!            get_fock logical is absent or is present and FALSE:
!            density matrix as input, potentials vector as output.                        
!                        

#include "priunit.h"
#include "pi.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#ifdef prg_dirac
#include "dcbgrd.h"
#endif
#include "ccom.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbham.h"
#include "cbisol.h"
#include "shells.h"
#include "onecom.h"
#include "primit.h"
#include "symind.h"
#include "symmet.h"
#include "lmns.h"
#include "nuclei.h"

! Passed variables
      integer, intent(in) :: nr_points
      real(8)             :: fock_or_dmat(ntbas(0), ntbas(0), nz)
! Triggers whether we want PCM Fock matrix contribution or just the
! potential.      
      logical, optional   :: get_fock 
      real(8)             :: pot_or_chg(nr_points)
      real(8)             :: points(3, nr_points), work(*)
      integer             :: iprint, lwork

! Local variables
      real(8), allocatable :: ader(:, :, :)
      real(8) :: tolog, tols, factor = 1.0d0
      integer :: ishela, ishelb, ica, icb, ia, ib
      integer :: multa, multb, nhktab, kab
      integer :: ipoint, idxmax
      logical :: do_fock

      if (nbsym > 1) then
         call quit('Electrostatic potential calculation cannot          &
      handle symmetry!')
      end if
      if (BSS.or.x2c) then
         call quit('Electrostatic potential calculation not working     &
      with X2C Hamiltonian!')
      end if

! Fock matrix contribution or potential calculation?      
      do_fock = .false. 
      if (present(get_fock)) then
        if (get_fock == .true.) then
            do_fock = .true.
        end if
      end if
                
      tols = thrs**2
      tolog = 2 * log(thrs)
! mc is the number of components, as defined in include/dcbham.h
! Decide whether we loop on LL block only or not
      if (mc == 1) then 
          idxmax = nlrgsh
      else 
          idxmax = kmax
      end if
! Loop over bras <mu|      
      idena = 0
      do ishela = 1, idxmax
        ica   = lclass(ishela)
        nhkta = nhkt(ishela)
        khkta = khkt(ishela)
        kckta = kckt(ishela)
        call lmnval(nhkta, kckta, lvalua, mvalua, nvalua)
        ncenta = ncent(ishela)
        icenta = nucnum(ncenta, 1)
        mula   = istbao(ishela)
        multa  = mult(mula)
        nuca   = nuco(ishela)
        numcfa = numcf(ishela)
        jsta   = jstrt(ishela)
        corax  = cent(ishela, 1, 1)
        coray  = cent(ishela, 2, 1)
        coraz  = cent(ishela, 3, 1)
! Loop over kets |nu>
        idenb = 0        
        do ishelb = 1, ishela
          icb   = lclass(ishelb)
          ldiag = ishela .eq. ishelb
          nhktb = nhkt(ishelb)
          khktb = khkt(ishelb)
          kcktb = kckt(ishelb)
          call lmnval(nhktb, kcktb, lvalub, mvalub, nvalub)
          ncentb = ncent(ishelb)
          nhktab = nhkta + nhktb
          mulb   = istbao(ishelb)
          multb  = mult(mulb)
          nucb   = nuco(ishelb)
          numcfb = numcf(ishelb)
          jstb   = jstrt(ishelb)
          corbx  = cent(ishelb, 1, 1)
          corby  = cent(ishelb, 2, 1)
          corbz  = cent(ishelb, 3, 1)
          khktab = khkta * khktb
          kcktab = kckta * kcktb
          mab    = ior(mula, mulb)
          kab    = iand(mula ,mulb)
          hkab   = fmult(kab)
! Calculate -<mu|1/|r-r_i||nu> integrals for each shell pair mu, nu and point r_i,
! then contract with the density matrix. The minus sign accounts for the
! charge of the electron.
          if (ica == icb) then 
              allocate(ader(kckta, kcktb, nr_points))
              ader = 0.0d0
              call vc_shell(nr_points, ader, tolog, tols,               &
                    points, iprint, work, lwork)
              if (ishela /= ishelb) then 
                  ader = 2.0d0*ader
              end if
! We multiply by +2.0d0. The D matrix refers to the alpha part only.
              do ipoint = 1, nr_points
              if (do_fock) then
                 do ib = 1, kcktb
                   do ia = 1, kckta 
               fock_or_dmat(idena+ia,idenb+ib, 1) = fock_or_dmat(idena+ia,idenb+ib, 1) &
                    - ader(ia, ib, ipoint) * pot_or_chg(ipoint)
                   end do
                 end do
               else
               pot_or_chg(ipoint) = pot_or_chg(ipoint) + 2.0d0                    * &
                    sum(fock_or_dmat(idena+1:idena+kckta, idenb+1:idenb+kcktb, 1) * &
                            ader(1:kckta, 1:kcktb, ipoint))
                end if
              end do 
              deallocate(ader)
          end if  
        idenb = idenb + khktb * multb
       end do
       idena = idena + khkta * multa
      end do

      if (iprint > 5) then
          do ipoint = 1, nr_points
            print *, "Electronic ESP @point", ipoint, pot_or_chg(ipoint)
          end do
      end if

      end subroutine

      subroutine vc_shell(nr_points, ader, tolog, tols,                 &
     &                    points, iprint, work, lwork)
!      
! Calculates the contribution for one primitive orbital set.
!
! RDR 090312 Clean-up.      
!
#include "priunit.h"
#include "pi.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#ifdef prg_dirac
#include "dcbgrd.h"
#endif
#include "cbisol.h"
#include "onecom.h"
#include "ader.h"
#include "primit.h"

! Passed variables
      integer, intent(in) :: nr_points, iprint, lwork
      real(8) :: points(3, nr_points), ader(kckta, kcktb, nr_points)
      real(8) :: work(*)
      real(8) :: tolog, tols, factor=1.0

! Local variables
      real(8), allocatable :: ahgtf(:, :, :), odc(:, :, :, :, :, :)
      real(8) :: difab(3), difcp(3), corp(3), origin(3),                &
     &           cora(3), corb(3)
      real(8) :: distab, conta, expa, contb, expb, expp, exppi,         &
     &           expabq, saab, asaab, saab13, expapi, expbpi
      integer :: jmaxa, jmaxb, jmaxd, jmaxt, jmaxm,               &
     &           ipoint, iprima, iprimb, jprima, jprimb
      integer :: idummy
      real(8) :: pval
      

! Allocation
      jmaxd = 2
      jmaxa = nhkta - 1
      jmaxb = nhktb - 1
      jmaxt = jmaxa + jmaxb + jmaxd
      jmaxm = 0
      jmax = jmaxa + jmaxb
      allocate(odc(0:jmaxa,0:jmaxb,0:jmaxt,0:jmaxd,0:jmaxm,3))
      ader(:, :, :) = 0.0d0

      cora(1) = corax
      cora(2) = coray
      cora(3) = coraz
      corb(1) = corbx
      corb(2) = corby
      corb(3) = corbz
      difab(:) = cora(:) - corb(:)
      distab = difab(1)*difab(1) + difab(2)*difab(2) + difab(3)*difab(3)

! Loop over primitive orbitals
! Shell a       
      do iprima = 1, nuca
        jprima = jsta + iprima
        conta = priccf(jprima, numcfa)
        expa = priexp(jprima)
! Shell b        
        do iprimb = 1, nucb
          jprimb = jstb + iprimb
          contb = priccf(jprimb, numcfb)
          expb = priexp(jprimb)
          expp = expa + expb
          exppi = 1.0d0 / expp
          expabq = expa * expb * exppi 
          saab = conta * contb * exp(-expabq * distab) 
          asaab = abs(saab)
!          if (expabq * distab < tolog) then 
!                cycle
!          end if
!          if (asaab < tols) then 
!                cycle
!          end if
          saab13 = sign(asaab**(1.0d0/3.0d0), saab)
          expapi  = expa * exppi
          expbpi  = expb * exppi
          corp(:) = expapi * cora(:) + expbpi * corb(:)
! Calculate the Overlap Distribution Coefficients
          idummy = 0
          call getodc(odc, jmaxa, jmaxb, jmaxt, jmaxd, jmaxm, .false.,  &
     &           .false., onecen, expa, expb, iprint, saab13, exppi,    &
     &           work, lwork, corp(1), corp(2), corp(3),                &
     &           .true.,.false.,origin,idummy)
          do ipoint = 1, nr_points
            difcp(:) = points(1:3,ipoint) - corp(:)
            allocate(ahgtf(0:jmax, 0:jmax, 0:jmax))
            ahgtf = 0.0d0
            call vnuc(ahgtf, jmax, expp, difcp, iprint)
            call cart_vc(odc, jmaxa, jmaxb, jmaxt, jmaxd, jmaxm,        &
     &            ader(:, :, ipoint), ahgtf, iprint)
            deallocate(ahgtf)
          end do ! Close loop over grid points
        end do ! Close loop over second shell
      end do ! Close loop over first shell
      deallocate(odc)

      return

      end subroutine

      subroutine vnuc(ahgtf, jmax, pval, cp, iprint)
!
! This subroutine calculates the R integrals as defined by
! McMurchie and Davidson in J. Comp. Phys. 26 (1978) 218.
! The recursion formulas (4.6) - (4.8) are used.
!
! JHS 260308 Only slightly slower than the implementation hernai 
!            in abacus/her1car.f of TUH
! RDR 090312 Clean-up.      
!     
#include "priunit.h"
#include "pi.h"
#include "maxaqn.h"
#include "gamcom.h"
! Passed variables
      integer :: jmax, iprint
      real(8) :: ahgtf(0:jmax, 0:jmax, 0:jmax), cp(3)
      real(8) :: pval

! Local variables
      real(8), allocatable :: r(:, :, :, :)
      real(8) :: factor, prod
      integer :: jval, t, u, v

! Initialization
      allocate(r(0:jmax, 0:jmax, 0:jmax, 0:jmax))
      r = 0.0d0

! Incomplete gamma function
      wval = pval * (cp(1) * cp(1) + cp(2) * cp(2) + cp(3) * cp(3))
      jmax0 = jmax
      call gamfun
      
! Calculate r(0, 0, 0, jval)
      factor = (2.0d0 * pi) / pval
      do jval = 0, jmax
         fjw(jval)         =   factor * fjw(jval)
         factor            = - 2.0d0 * pval * factor
         r(0, 0, 0, jval) =   fjw(jval)
      end do

! Calculate r(t, u, v, jval) 
      do jval = jmax, 1, -1
         do v = 0, jmax - jval
            prod = -cp(3) * r(0, 0, v, jval)
            if (v > 0) then 
                prod = prod + v * r(0, 0, v - 1, jval)
            end if    
            r(0, 0, v + 1, jval - 1) = prod
            do u = 0, jmax - jval - v
               prod = -cp(2) * r(0, u, v, jval)
               if (u > 0) then 
                   prod = prod + u * r(0, u - 1, v, jval)
               end if    
               r(0, u + 1, v, jval - 1) = prod
               do t = 0, jmax - jval - u - v
                  prod = -cp(1) * r(t, u, v, jval)
                  if (t > 0) then 
                      prod = prod + t * r(t - 1, u, v, jval)
                  end if
                  r(t + 1, u, v, jval - 1) = prod
               end do
            end do
         end do
      end do
! The nuclear attraction integrals are given as r(t, u, v, 0)
      do t = 0, jmax
         do u = 0, jmax
            do v = 0, jmax
               ahgtf(t, u, v) = - r(t, u, v, 0)
            end do
         end do
      end do
      deallocate(r)

      end subroutine

      subroutine cart_vc(odc, jmaxa, jmaxb, jmaxt, jmaxd, jmaxm,        &
     &                 ader, ahgtf, iprint)
!      
! Originally written by JHS.
!      
! RDR 060312 Clean-up.
!

#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "onecom.h"
#include "lmns.h"

! Passed variables
      real(8) :: ahgtf(0:jmax, 0:jmax, 0:jmax)
      real(8) :: ader(kckta, kcktb, 1),                                 &
     &           odc(0:jmaxa, 0:jmaxb, 0:jmaxt, 0:jmaxd, 0:jmaxm, 3)
      integer :: jmaxa, jmaxb, jmaxt, jmaxd, jmaxm, iprint

! Local variables
      real(8) :: naint, ev, ee, eee
      integer :: icompa, lvala, mvala, nvala,                           &
     &           icompb, lvalb, mvalb, nvalb,                           &
     &           t, u, v


      do icompa = 1,kckta
        lvala = lvalua(icompa)
        mvala = mvalua(icompa)
        nvala = nvalua(icompa)
        do icompb = 1,kcktb
          lvalb = lvalub(icompb)
          mvalb = mvalub(icompb)
          nvalb = nvalub(icompb)
          naint = 0.0d0
          do v = 0, nvala + nvalb
            ev = odc(nvala,nvalb,v,0,0,3)
            do u = 0, mvala + mvalb
              ee = odc(mvala,mvalb,u,0,0,2)*ev
              do t = 0, lvala + lvalb
                eee = odc(lvala,lvalb,t,0,0,1)*ee
                naint = naint + eee*ahgtf(t,u,v)
              end do
            end do
          end do
          ader(icompa, icompb, 1) = ader(icompa, icompb, 1) + naint
        end do
      end do

      return
      end subroutine

      end module
