!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PH_TESTER()
      IMPLICIT INTEGER (A-Z)

#include "../relccsd/symm.inc"   ! provide NREP, MULTB and S-arrays

      integer, dimension(:), allocatable     :: inverse

      call pst('Performing ph/ph symmetry test (main block)+')
!
! determine the inverse f/b ireps for each f/b irep
! and scrutinize mutab laws
!
      allocate(inverse(2*nrep))
      inverse = 0
      do i=1,nrep
        do j=1,nrep
          if(multb(i,j,1).eq.1) then ! f x f inverse
            inverse(i)=j
            write(*,*) 'F inverse of',i,' is ',j
          endif
        enddo
        do j=1,nrep
          if(multb(i+nrep,j+nrep,1).eq.1) then  ! b x b inverse
            inverse(i+nrep)=j
            write(*,*) 'B inverse of',i,' is ',j
          endif
        enddo
      enddo
      do i=1,nrep
        if(multb(inverse(i),i,1).ne.1) stop 'PH error #01.'
        if(multb(i,i,2).ne.1) stop 'PH error #02.'
        if(multb(inverse(i+nrep)+nrep,i+nrep,1).ne.1) 
     &                             stop 'PH error #03.'
        if(multb(i+nrep,i+nrep,2).ne.1) stop 'PH error #04.'
      enddo
      write(*,*) '  ... inverse test passed.'
      write(*,*)
      write(*,*)
!
! determine the irep order on the bra and ket side for a given
! final state symmetry irep 'desrep' (destination irep)
!
! bra side: <a||i> = desrep*
! ket side: |b><j| = desrep
! We strictly keep the convention that bra sides are to be complex
! conjugated.
! At the end for each matrix block the following holds:
!
!          desrep*  x  desrep = Gamma_0 (totally symm. irep)
!
!  we first loop through the ket side |b><i|
!
!
      Do desrep = 1,nrep

        write(*,*) 'Testing <P|H>|P><H| block for desrep:',desrep

        Do jrep = 1,nrep    !jrep is to be interpreted as the gamma*_j number!
          brep = MULTB(desrep+nrep,inverse(jrep),1)
          if(multb(brep,jrep,1).ne.desrep) 
     &              Stop 'PH_TESTER error #10.'


          airep = inverse(desrep+nrep)     !airep must now be desrep* (bosonic!)
          locnum=0
          DO IREP = 1, NREP
            DO AREP = 1, NREP
              xrep = MULTB(AREP,IREP,2) ! form <gamma_a| x |gamma_i>
              If(xrep.eq.airep) then
                 write(*,*) ' G*_a:',arep,' G_i: ',irep,
     &                      ' G_b:',brep,' G*_j: ',jrep
!                do i=1,no(IREP)
!                do a=1,nv(AREP)
!                  write(*,*) '    a(absindx),i(absindx):',
!    &               io(nrep)+iv(arep)+a,io(irep)+i
!                  locnum=locnum+1
!                enddo
!                enddo
                 if(MULTB(xrep+nrep,desrep+nrep,1).ne.1) then
                    stop 'PH_TESTER error #11.'
                 endif
              Endif  !xrep.eq.airep
            ENDDO
          ENDDO
!         write(*,*) 'with',locnum,' valid p/h combinations:'
!         write(*,*)
        Enddo !jrep
      Enddo !desrep
      write(*,*) '... passed'


      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PPHH_TESTER()
      IMPLICIT INTEGER (A-Z)

#include "../relccsd/symm.inc"   ! provide NREP, MULTB and S-arrays

      integer, dimension(:), allocatable     :: inverse

      call pst('Performing 2p2h/ph symmetry test (coupling block)+')
!
! determine the inverse f/b ireps for each f/b irep
!
      allocate(inverse(2*nrep))
      inverse = 0
      do i=1,nrep
        do j=1,nrep
          if(multb(i,j,1).eq.1) then ! f x f inverse
            inverse(i)=j
            write(*,*) 'F inverse of',i,' is ',j
          endif
        enddo
        do j=1,nrep
          if(multb(i+nrep,j+nrep,1).eq.1) then  ! b x b inverse
            inverse(i+nrep)=j
            write(*,*) 'B inverse of',i,' is ',j
          endif
        enddo
      enddo
      do i=1,nrep
        if(multb(inverse(i),i,1).ne.1) stop 'PH error #01.'
        if(multb(i,i,2).ne.1) stop 'PH error #02.'
        if(multb(inverse(i+nrep)+nrep,i+nrep,1).ne.1) 
     &                             stop 'PH error #03.'
        if(multb(i+nrep,i+nrep,2).ne.1) stop 'PH error #04.'
      enddo
      write(*,*) '  ... inverse test passed.'
      write(*,*)
!
! We loop over the ket side: |b><j| = desrep
! and scrutinize the <c|<d| |k>|l> loop structure. Since cdkl is on the
! bra side (row index) we must have
! G*c x G*d x Gk x Gl = desrep*
! As always, for the coupling block the following holds:
!
!    cdkl = desrep*  x  desrep = Gamma_0 (totally symm. irep)
!
!
      Do desrep = 1,nrep

        write(*,*) 'Testing <P|<P|H>|H> block for bra(desrep):',desrep

        desrepi = inverse(desrep+nrep)
        do klrep = 1,nrep
          cdrep = MULTB(desrepi+nrep,inverse(klrep+nrep)+nrep,1)  !cdrep is a bra number!
          if(MULTB(cdrep+nrep,klrep+nrep,1).ne.desrepi) then
            stop 'PPHH error #001.'
          endif

          DO LREP = 1, NREP
            KREP = MULTB(LREP,KLREP+NREP,2)
            if(multb(krep,lrep,1).ne.klrep) then
              stop 'PPHH error #002.'
            endif
            DO DREP = 1, NREP  !drep is understood as a bra loop
              CREP = MULTB(cdrep+nrep,inverse(drep),1)   !crep is a bra number!
              if(multb(crep,drep,1).ne.cdrep) then
                stop 'PPHH error #003.'
              endif
            ENDDO  ! DREP
          ENDDO  ! LREP
        enddo  ! klrep

      Enddo !desrep
      write(*,*) '... passed'
      write(*,*)

      write(*,*) 'Now testing full coupling block symmetry with deltas.'
!
! all previously tested irep multiplications can now be taken over
! and need not be tested again. This shortens the loops considerably and
! we can concentrate on the deltas.
!
      Do desrep = 1,nrep

        write(*,*) 'Testing <P|H>|P><H| |P><H| block for desrep:',desrep

!************* main program loop
        Do jrep = 1,nrep    !jrep is a bra index
          brep = MULTB(desrep+nrep,inverse(jrep),1) !brep is fermionic ket

!************* xcolmake loop
          desrepi = inverse(desrep+nrep) ! desrepi is a bosonic bra index
          do klrep = 1,nrep ! klrep is a bosonic ket index
            cdrep = MULTB(desrepi+nrep,inverse(klrep+nrep)+nrep,1)  !cdrep is a bra index!
            DO LREP = 1, NREP ! lrep is a ket index!
              KREP = MULTB(LREP,KLREP+NREP,2) ! krep is a ket index!
              DO DREP = 1, NREP  !drep is a bra index!
                CREP = MULTB(cdrep+nrep,inverse(drep),1)   !crep is a bra index!

!************ from here we have all ireps for the delta treatment

!|  Case A:  < dj||kl > delta_bc, we calc a VOOO-offset.

                if(brep.eq.inverse(crep)) then
!                 write(*,*) 'delta A condition met.'
                  djrep=multb(drep,jrep,1)  ! djrep is a bra index!
                  if(multb(djrep+nrep,klrep+nrep,1).ne.1) then
                    stop 'PPHH error #044.'
                  endif
                  if(NV(brep).ne.NV(crep)) stop 'PPHH error #744.'
                endif

!|  Case B:  < cj||kl > delta_bd, we calc a VOOO-offset.

                if(brep.eq.inverse(drep)) then
!                 write(*,*) 'delta B condition met.'
                  cjrep=multb(crep,jrep,1)  ! cjrep is a bra index!
                  if(multb(cjrep+nrep,klrep+nrep,1).ne.1) then
                    stop 'PPHH error #045.'
                  endif
                  if(NV(brep).ne.NV(drep)) stop 'PPHH error #745.'
                endif

!|  Case C:  < cd||bl > delta_jk, we calc a VOVV offset

                if(krep.eq.inverse(jrep)) then
!                 write(*,*) 'delta C condition met.'
                  cdrep=multb(crep,drep,1)  ! cdrep is a bra index!
                  blrep=multb(brep,lrep,1)  ! blrep is a ket index!
                  if(multb(cdrep+nrep,blrep+nrep,1).ne.1) then
                    stop 'PPHH error #046.'
                  endif
                  if(NO(krep).ne.NO(jrep)) stop 'PPHH error #746.'
                endif

!|  Case D:  < cd||bk > delta_jl, we calc a VOVV offset

                if(lrep.eq.inverse(jrep)) then
!                 write(*,*) 'delta D condition met.'
                  cdrep=multb(crep,drep,1)  ! cdrep is a bra index!
                  bkrep=multb(brep,krep,1)  ! bkrep is a ket index!
                  if(multb(cdrep+nrep,bkrep+nrep,1).ne.1) then
                    stop 'PPHH error #047.'
                  endif
                  if(NO(lrep).ne.NO(jrep)) stop 'PPHH error #746.'
                endif

              ENDDO  ! DREP
            ENDDO  ! LREP
          enddo  ! klrep

        Enddo  ! jrep

      Enddo  ! desrep
      write(*,*) '... passed'
      write(*,*)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
