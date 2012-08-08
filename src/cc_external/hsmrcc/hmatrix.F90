!-------------------------------------------------------------------------------------------------
!
!     This routine calculates matrix elements of one and two body operator
!
!     H = H1 + V2
!
!     H1 = sum(i,j) h1(i <- j) i+ j
!
!     V2 = sum(i>j, k>l) <i,j|v2|k,l> i+ j+ l k
!
!     between any two determinantal states |left> and |right>
!
!        hmat(left,right) = <left| H |right>
!
!                         = delta(left,right)*core_energy  +
!
!                           <left| F1 |right> +
!
!                           <left| F2 |right>
!
!        where F1 and F2 are second-quantized operators of the
!        same form as H1 and V2, but defined over only non-core
!        part of orbital space. F2 matrix elements are same as
!        that of V2, and F1 is the modified fock-operator
!        with respect to core.
!
!        h1(i <- j)
!
!            is stored as H1(i,j)
!
!        <i,j|v2|k,l> for i>j, and k>l
!
!            is stored as V2(ij,kl)
!
!            ij = MXMN(i,j)
!            kl = MXMN(k,l)
!
!     Special Notes :
!
!     1. Core used in this subroutine only means core part of all
!        occupation vectors passed in occstate array. It has nothing
!        to do with core orbitals used else where in the program.
!        Infact, core orbitals may be anywhere in the occupation vectors.
!
!     2. This routine has no dependency on "gtci.com". Therefore, this
!        routine is not a part of gtci and can be used elsewhere without
!        any problems.
!
!
!  mi: Called from modules/vacmod.F90 - subr. GenVac
!
!
!-------------------------------------------------------------------------------------------------
      subroutine hmatrix           &
                 (                 &
                    nstate_in,     &
                    nspinorb,      &
                    nspinorb2,     &
                    ldhmat,        &
                    hmat,          &
                    occstate,      &
                    H1, V2         &
                 )

      use IOModule
      implicit none
      integer nstate_in
      integer nspinorb
      integer nspinorb2
      integer ldhmat
      double precision  hmat(ldhmat,nstate_in)
      integer occstate(nspinorb,nstate_in)
      double precision  H1(nspinorb,nspinorb)
      double precision  V2(nspinorb2,nspinorb2)

!     Important Local Variables

      integer ncore, nopen
      double precision  energy_core
      integer occcore(nspinorb)
      integer coreorb(nspinorb)
      integer openorb(nspinorb)
      integer occright(nspinorb)
      integer occleft(nspinorb)
      double precision  F1(nspinorb,nspinorb)

!     Local variables

      double precision  dtmp
      integer i, j, k, l
      integer kcore
      integer ik, jk
      integer istate, ileft, iright
      integer norbdiff
      integer idum
      integer iij, ikl
      integer iopen, jopen, lopen
      integer idiff
      integer ij, kl, il
      integer ist, ien, ist1, ien1, ist2, ien2
      integer s, is
      integer sgn, sgn1, sgn2

!      logical dbg

!     statement function

      integer ix, iy
      integer MXMN
      MXMN(ix,iy) = ((MAX(ix,iy)-2)*(MAX(ix,iy)-1))/2 + MIN(ix,iy)

!     set dbg flag

!      dbg = .false.

!     print banner

!     if (dbg) then
      if ( get_iprint() >= 6 ) then
         write(UOUT,9000)
      endif
 9000 format(/,&
             "********************************************",/, &
             "Program to Calculate Full Hamiltonian Matrix",/, &
             "    over a small number of determinants     ",/, &
             "********************************************",/,/)

!     first find occupation pattern of core

      do i = 1, nspinorb
         occcore(i) = 1
      enddo

      do istate = 1, nstate_in
         do i = 1, nspinorb
            if ((occstate(i,istate).ne.0).and.(occstate(i,istate).ne.1)) then
                write(UOUT,*) "hmatrix: invalid entry in occstate array"
                !call flush(UOUT)
                call flushout
               ! call abort
                call quit('hmatrix: invalid entry in occstate array')
               ! return
             endif
            occcore(i) = min(occcore(i),occstate(i,istate))
         enddo
      enddo

!     define core and open orbitals

      do i = 1, nspinorb
         coreorb(i) = 0
         openorb(i) = 0
      enddo

      ncore = 0
      nopen = 0
      do i = 1, nspinorb
         if (occcore(i).eq.1) then
            ncore = ncore + 1
            coreorb(ncore) = i
         else
            nopen = nopen + 1
            openorb(nopen) = i
         endif
      enddo

!     calculate modified fock-matrix
!     also calculate core energy

      energy_core = zero
      do i = 1, nspinorb
      do j = 1, nspinorb
         F1(j,i) = NaN
      enddo
      enddo

      do i = 1, nspinorb
      do j = 1, nspinorb
         dtmp = zero
         do 100 kcore = 1, ncore
            k = coreorb(kcore)
            if ((k.eq.i).or.(k.eq.j)) goto 100
            ik = MXMN(i,k)
            jk = MXMN(j,k)
            if ((k.gt.max(i,j)).or.(k.lt.min(i,j))) then
               dtmp = dtmp + V2(jk,ik)
            else
               dtmp = dtmp - V2(jk,ik)
            endif
 100     continue
         F1(j,i) = H1(j,i) + dtmp
         if (i.eq.j) then
            if (occcore(i).eq.1) energy_core = energy_core + half*(F1(j,i)+H1(j,i))
         endif
      enddo
      enddo

!     now we are all set to build hamiltonian matrix
!     using slater rules.
!
!     norbdiff is the number of orbitals by which
!     right and left state differ.

      do iright = 1, nstate_in
      do ileft  = 1, nstate_in
         hmat(ileft,iright) = zero
      enddo
      enddo

      do iright = 1, nstate_in
         do i = 1, nspinorb
            occright(i) = occstate(i,iright)
         enddo
         do ileft = 1, nstate_in
            do i = 1, nspinorb
               occleft(i) = occstate(i,ileft)
            enddo
            i = 0
            j = 0
            k = 0
            l = 0
            iij = 0
            ikl = 0
            norbdiff = 0
            do 200 iopen =1 , nopen
               idum = openorb(iopen)
               idiff = occleft(idum) - occright(idum)
               if (idiff.eq.0) goto 200
               norbdiff = norbdiff + abs(idiff)
               if (norbdiff.gt.4) goto 200
               if (idiff.eq.1) then
                  ikl = ikl + 1
                  if (ikl.eq.1) k = idum
                  if (ikl.eq.2) then
                     l = k
                     k = idum
                  endif
               else
                  iij = iij + 1
                  if (iij.eq.1) i = idum
                  if (iij.eq.2) then
                     j = i
                     i = idum
                  endif
               endif
 200        continue

            if ((mod(norbdiff,2).eq.0).and.((norbdiff/2).eq.0)) goto 300
            if ((mod(norbdiff,2).eq.0).and.((norbdiff/2).eq.1)) goto 400
            if ((mod(norbdiff,2).eq.0).and.((norbdiff/2).eq.2)) goto 500
            if ((mod(norbdiff,2).eq.0).and.((norbdiff/2).gt.2)) goto 600

!           it should not get here

            write(UOUT,*) "hmatrix: norbdiff has to be even"
            !call flush(UOUT)
            call flushout
            !call abort
            call quit('hmatrix: norbdiff has to be even!')
            !return

!           left and right are same
!
!           <left|H|right> = E_core
!
!                               +
!
!                            sum(iopen, occupied in |right>) * [ <i|F|i> ]
!
!                               +
!
!                            sum(iopen > jopen, occupied in |right>) * [ <ij||ij> ]

 300        continue

            if (ileft.ne.iright) then
               write(UOUT,*) "same state appears in "
               write(UOUT,*) "state list more than once"
               write(UOUT,*) "ileft  =",ileft
               write(UOUT,*) "iright =",iright
               do i = 1, nspinorb
                  write(UOUT,fmt='(3(i5,1x))') i, occstate(i,ileft), occstate(i,iright)
               enddo
               !call flush(UOUT)
               call flushout
               !call abort
               call quit('hmatrix: state list more than once !')
               !return
            endif

!           calculate matrix element

            dtmp = energy_core
            do 320 iopen = 1, nopen
               i = openorb(iopen)
               if (occright(i).ne.1) goto 320
               dtmp = dtmp + F1(i,i)
               do 310 jopen = 1, (iopen-1)
                  j = openorb(jopen)
                  if (occright(j).ne.1) goto 310
                  ij = MXMN(i,j)
                  dtmp = dtmp + V2(ij,ij)
 310           continue
 320        continue

!           add to corresponding matrix element

            hmat(ileft,iright) = dtmp

            goto 700


!           left and right differ by one orbital
!
!           <left|H|right> = sgn * [ <k|F|i> + sum(lopen) {<kl||il>} ]
!
!           where,
!
!                sgn = <left|(k+ i)|right>

 400        continue

          ! if (dbg) then
               if (k.eq.i) then
                  write(UOUT,*) "hmatrix: occupancy error"
                  !call flush(UOUT)
                  call flushout
                  !call abort
                  call quit('hmatrix: occupancy error!')
                  !return
               endif
         !   endif

!           now we have k /= i . relaitive order of
!           k and i determines how sign is calculated.

            if (k.lt.i) goto 410
            if (k.gt.i) goto 420

!           should not get here

            write(UOUT,*) "hmatrix: occupancy error"
           ! call flush(UOUT)
            call flushout
            !call abort
            call quit('hmatrix: occupancy error')
            !return

!           k < i
!
!           k+ i => k+ i

 410        continue
            ist = k+1
            ien = i-1
            sgn = (+1)* ( (-1)**0 )
            goto 450

!           k > i
!
!           k+ i => i  k+

 420        continue
            ist = i+1
            ien = k-1
            sgn = (-1)* ( (-1)**1 )
            goto 450

!           calculate overall sign

 450        continue

            s = 0
            do is = ist, ien
               s = s + occright(is)
            enddo
            sgn1 = (-1)**s

            sgn = sgn*sgn1

!           sign is now determined
!           calculate matrix element

            dtmp = F1(k,i)
            do 460 lopen = 1, nopen
               l = openorb(lopen)
               if (occright(l).eq.0) goto 460
               if (occleft(l).eq.0 ) goto 460
               kl = MXMN(k,l)
               il = MXMN(i,l)
               if ((l.gt.max(k,i)).or.(l.lt.min(k,i))) then
                  dtmp = dtmp + V2(kl,il)
               else
                  dtmp = dtmp - V2(kl,il)
               endif
 460        continue

!           add to corresponding matrix element

            if (sgn.eq.1) then
               hmat(ileft,iright) = +dtmp
            else
               hmat(ileft,iright) = -dtmp
            endif

            goto 700

!           left and right differ by two orbitals
!
!           <left|H|right> = sgn * <kl||ij>
!
!           where,
!
!                sgn = <left|(k+ l+ j  i)|right>

 500        continue

           ! if (dbg) then
             if (get_iprint() > 5) then
               if ((i.le.j).or.(k.le.l)) then
                  write(UOUT,*) "hmatrix: (i,j) or (k,l) not ordered"
                  !call flush(UOUT)
                  call flushout
                  !call abort
                  call quit('hmatrix: (i,j) or (k,l) not ordered')
                  !return
               endif
               if ((i.eq.k).or.(i.eq.l).or.(j.eq.k).or.(j.eq.l)) then
                  write(UOUT,*) "hmatrix: occupancy error"
                 ! call flush(UOUT)
                  call flushout
                  !call abort
                  call quit('hmatrix: occupancy error !')
                  !return
               endif
             endif

            ij = MXMN(i,j)
            kl = MXMN(k,l)
            if (V2(kl,ij).eq.zero) goto 600


!           now we have i > j and k > l, and are all distinct.
!           the order of k and l relative to i and j now
!           determines how sign is calculated.

            if (l.gt.i) goto 510
            if (l.gt.j) goto 520
            if (l.lt.j) goto 540
            
!           should not get here

            write(UOUT,*) "hmatrix: Something wrong when norbdiff = 4"
           ! call flush(UOUT)
            call flushout
           ! call abort
            call quit('hmatrix: Something wrong when norbdiff = 4 !')
            return

!           j < i < l < k
!
!           k+ l+ j  i => j  i  l+ k+

 510        continue
            ist1 = j+1
            ien1 = i-1
            ist2 = l+1
            ien2 = k-1
            sgn = (-1)*( ((-1)**1) * ((-1)**2) * ((-1)**2) )
            goto 580

!           j < l < i < k
!
!           k+ l+ j  i => j  l+ i  k+

 520        continue
            if (k.lt.i) goto 530
            ist1 = j+1
            ien1 = l-1
            ist2 = i+1
            ien2 = k-1
            sgn = (+1)*( ((-1)**1) * ((-1)**1) * ((-1)**2) )
            goto 580

!           j < l < k < i
!
!           k+ l+ j  i => j  l+ k+ i 

 530        continue
            ist1 = j+1
            ien1 = l-1
            ist2 = k+1
            ien2 = i-1
            sgn = (-1)*( ((-1)**1) * ((-1)**1) * ((-1)**1) )
            goto 580

!           l < j < i < k
!
!           k+ l+ j  i => l+ j  i  k+

 540        continue
            if (k.lt.i) goto 550
            ist1 = l+1
            ien1 = j-1
            ist2 = i+1
            ien2 = k-1
            sgn = (-1)*( ((-1)**0) * ((-1)**1) * ((-1)**2) )
            goto 580

!           l < j < k < i
!
!           k+ l+ j  i => l+ j  k+ i 

 550        continue
            if (k.lt.j) goto 560
            ist1 = l+1
            ien1 = j-1
            ist2 = k+1
            ien2 = i-1
            sgn = (+1)*( ((-1)**0) * ((-1)**1) * ((-1)**1) )
            goto 580

!           l < k < j < i
!
!           k+ l+ j  i => l+ k+ j  i 

 560        continue
            ist1 = l+1
            ien1 = k-1
            ist2 = j+1
            ien2 = i-1
            sgn = (-1)*( ((-1)**0) * ((-1)**0) * ((-1)**1) )
            goto 580

!           calculate overall sign

 580        continue

            s = 0
            do is = ist1, ien1
               s = s + occright(is)
            enddo
            sgn1 = (-1)**s

            s = 0
            do is = ist2, ien2
               s = s + occright(is)
            enddo
            sgn2 = (-1)**s

            sgn = sgn*sgn1*sgn2

!           sign is now determined, add to corresponding matrix element

            if (sgn.eq.1) then
               hmat(ileft,iright) = +V2(kl,ij)
            else
               hmat(ileft,iright) = -V2(kl,ij)
            endif

            goto 700


!           left and right differ by more than two orbitals
!
!           <left|H|right> = 0

 600        continue

            hmat(ileft,iright) = zero

            goto 700

!           matrix element is now computed, print it

 700        continue

            !if (dbg) then
            if ( get_iprint() >= 3 ) then
              write(UOUT,8990) norbdiff,ileft, iright, hmat(ileft,iright)
            endif

         enddo
      enddo

 8990 format(I3,"< ",I8," | H | ",I8," > = ",F25.18)

      return
      end subroutine hmatrix


