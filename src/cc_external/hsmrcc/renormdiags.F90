module RenormDiags

   use Parameters
   use Options
   use IOModule
   use VacMod
   use TAmps
   use TIndexMod
   use HEFFModule
   implicit none
   private
   public :: DrvRenorm

contains

subroutine DrvRenorm(vac)
!
! called from source/diagrams.F90/diagrams
!
   use TimingMod
   implicit none
   integer, intent(in) :: vac

   call QEnter('DrvRenorm')

   call Renorm_HEFF_OneBody(TimeTaken(:,RENORM),vac)
   call Renorm_HEFF_TwoBody(TimeTaken(:,RENORM),vac)

   call QExit('DrvRenorm')
   return

end subroutine DrvRenorm
   
subroutine Renorm_HEFF_OneBody(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: lamda, gaama
   integer :: tvac, exc
   double precision :: WINT
   real :: dtim, dtime, tarray(2)

   dtim = dtime(tarray)

   do lamda = 1, NVMO
   do gaama = 1, NVMO
      if (.not.IsHole(lamda+NCMO,vac)) cycle
      if (.not.IsPtcl(gaama,vac)) cycle
      exc  = Exct(lamda,gaama,vac)
      WINT = HEFF(exc,vac)
      if (DABS(WINT) <= 1.0D-06) cycle
      tvac = ExcVaccum(lamda,gaama,vac)
      read (UTOLD,rec=tvac) TOLD
      call Rnr_Sing_HEFF_OneBody(lamda,gaama,WINT)
      call Rnr_Dble_HEFF_OneBody(lamda,gaama,WINT)
   enddo
   enddo

   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine Renorm_HEFF_OneBody

subroutine Renorm_HEFF_TwoBody(TimeTaken,vac)

   implicit none
   real, intent(inout) :: TimeTaken(2)
   integer, intent(in) :: vac
   integer :: lamda, gaama, sigma, delta
   integer :: tvac, exc
   double precision :: WINT
   real :: dtim, dtime, tarray(2)

   call QEnter('Renorm_HEFF_TwoBody')

   if ( get_iprint() >= 8) then
    write(UOUT,'(2x,a)') 'Renorm_HEFF_TwoBody: beginning of 4-fold cycle'
   endif

   dtim = dtime(tarray)

   do lamda = 1, NVMO
   do sigma = 1, NVMO
   do gaama = 1, NVMO
   do delta = 1, NVMO
      if (.not.IsHole(lamda+NCMO,vac)) cycle
      if (.not.IsHole(sigma+NCMO,vac)) cycle
      if (.not.IsPtcl(gaama,vac)) cycle
      if (.not.IsPtcl(delta,vac)) cycle
      if ( lamda <= sigma ) cycle
      if ( gaama <= delta ) cycle
      tvac = ExcVaccum(lamda,sigma,gaama,delta,vac)
      if (Ms2Vac(tvac).ne.Ms2Vac(vac)) cycle
      if (IrpVac(tvac).ne.IrpVac(vac)) cycle
      exc  = Exct(lamda,sigma,gaama,delta,vac)
      WINT = HEFF(exc,vac)
      read (UTOLD,rec=tvac) TOLD
      call Rnr_Sing_HEFF_TwoBody(lamda,sigma,gaama,delta,+WINT)
      call Rnr_Sing_HEFF_TwoBody(lamda,sigma,delta,gaama,-WINT)
      call Rnr_Sing_HEFF_TwoBody(sigma,lamda,gaama,delta,-WINT)
      call Rnr_Sing_HEFF_TwoBody(sigma,lamda,delta,gaama,+WINT)
      call Rnr_Dble_HEFF_TwoBody(lamda,sigma,gaama,delta,+WINT)
      call Rnr_Dble_HEFF_TwoBody(lamda,sigma,delta,gaama,-WINT)
      call Rnr_Dble_HEFF_TwoBody(sigma,lamda,gaama,delta,-WINT)
      call Rnr_Dble_HEFF_TwoBody(sigma,lamda,delta,gaama,+WINT)
   enddo
   enddo
   enddo
   enddo

   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   if ( get_iprint() >= 8) then
    write(UOUT,'(2x,a)') 'Renorm_HEFF_TwoBody: end of the 4-fold cycle'
   endif

   call QExit('Renorm_HEFF_TwoBody')
   return

end subroutine Renorm_HEFF_TwoBody

subroutine Rnr_Sing_HEFF_OneBody(lamda,gaama,WINT)

   implicit none
   integer, intent(in) :: lamda, gaama
   double precision, intent(in) :: WINT
   double precision :: Rnr_Sum, TmpR, TmpT
   integer :: eqn, ti1, ti2
   integer :: i, a, j, b

   do i = 1, NHMO
   do a = 1, NPMO

      Rnr_Sum = dble(0)
      eqn = TINDEX(i,a)

      ! Diagram 01
      ! TERM : WT1T1
      !
      ! - sum(lamda,gaama)[w(lamda->gaama)*t(gaama,a)*t(i,lamda)]

      ti1 = TINDEX(gaama+NCMO,a)
      ti2 = TINDEX(i,lamda)
      Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)

      ! Diagram 02
      ! TERM : W1T2
      ! 
      ! + sum(lamda,gaama)[w(lamda->gaama)*t(i,gaama->a,lamda)]

      ti1 = TINDEX(i,gaama+NCMO,a,lamda)
      if ( ti1 > 0 ) then
         Rnr_Sum = Rnr_Sum + TOLD(ti1)
      else
         ti1 = IABS(ti1)
         Rnr_Sum = Rnr_Sum - TOLD(ti1)
      endif

      !
      ! Diagram 03
      ! TERM : W1T1
      !
      ! + sum(gaama)[w(lamda->gaama)*t(gaama->a)]

      if ( i == (lamda+NCMO) ) then
         ti1 = TINDEX(gaama+NCMO,a)
         Rnr_Sum = Rnr_Sum + TOLD(ti1)
      endif

      ! Diagram 04
      ! TERM : W1T1
      !
      ! - sum(lamda)[w(lamda->gaama)*t(i->lamda)]

      if ( a == gaama ) then
         ti1 = TINDEX(i,lamda)
         Rnr_Sum = Rnr_Sum - TOLD(ti1)
      endif

      ! All diagram done
      ! Now we do the multiplication
      ! Renormalization terms come with
      ! a negative sign, so make WINT negative.

      TNEW(eqn) = TNEW(eqn) + Rnr_Sum*(-WINT)

      ! Now Disconnected diagrams are summed.
      ! The procedure is standard. TmpR holds
      ! the contribution by Singles part.
      ! This is to be multiplied by
      ! [TOLD(j->b,vac)-TOLD(j->b,tvac)] = TmpT
      ! TOLD(j->b,vac) is stored in TSING.

      TmpR = -(Rnr_Sum*(-WINT))
      do j = 1, NHMO
      do b = 1, NPMO
         eqn = TINDEX(i,j,a,b)
         ti1 = TINDEX(j,b)
         TmpT = TSING(ti1) - TOLD(ti1)
         if ( eqn > 0 ) then
            TNEW(eqn) = TNEW(eqn) + TmpT*TmpR
         else
            eqn = IABS(eqn)
            TNEW(eqn) = TNEW(eqn) - TmpT*TmpR
         endif
      enddo
      enddo
   enddo
   enddo

   return

end subroutine Rnr_Sing_HEFF_OneBody

subroutine Rnr_Sing_HEFF_TwoBody(lamda,sigma,gaama,delta,WINT)

   implicit none
   integer, intent(in) :: lamda, sigma, gaama, delta
   double precision, intent(in) :: WINT
   double precision :: Rnr_Sum, TmpR, TmpT
   integer :: eqn, ti1, ti2
   integer :: i, a, j, b

   do i = 1, NHMO
   do a = 1, NPMO

      Rnr_Sum = dble(0)
      eqn = TINDEX(i,a)

      ! Diagram 01
      ! TERM : W2T1T2
      !
      ! - 1/2 sum(lamda,sigma,gaama,delta)[
      !          w(lamda,sigma->gaama,delta)*t(delta->a)*
      !          t(i,gaama->sigma,lamda)]
      !  1/2 factor removed by lamda > sigma

      if ( lamda > sigma ) then
         ti1 = TINDEX(delta+NCMO,a)
         ti2 = TINDEX(i,gaama+NCMO,sigma,lamda)
         if ( ti2 > 0 ) then
            Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
         else
            ti2 = IABS(ti2)
            Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
         endif
      endif

      ! Diagram 02
      ! TERM : W2T1T2
      !
      ! - 1/2 sum(lamda,sigma,gaama,delta)[
      !          w(lamda,sigma->gaama,delta)*t(i->sigma)*
      !          t(gaama,delta->lamda,a)]
      !  1/2 factor removed by gaama > delta
      !

      if ( gaama > delta ) then
         ti1 = TINDEX(i,sigma)
         ti2 = TINDEX(gaama+NCMO,delta+NCMO,lamda,a)
         if ( ti2 > 0 ) then
            Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
         else
            ti2 = IABS(ti2)
            Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
         endif
      endif

      ! Diagram 03
      ! TERM : W2T2
      !
      ! + 1/2 sum(gaama,delta,sigma)[w(lamda,sigma->gaama,delta)*
      !                          t(gaama,delta->a,sigma)]
      !
      !  1/2 factor removed by gaama > delta

      if ( i == lamda+NCMO ) then
      if ( gaama > delta ) then
         ti1 = TINDEX(gaama+NCMO,delta+NCMO,a,sigma)
         if ( ti1 > 0 ) then
            Rnr_Sum = Rnr_Sum + TOLD(ti1)
         else
            ti1 = IABS(ti1)
            Rnr_Sum = Rnr_Sum - TOLD(ti1)
         endif
      endif
      endif

      ! Diagram 04
      ! TERM : W2T2
      !
      ! - 1/2 sum(lamda,sigma,delta)[w(lamda,sigma->gaama,delta)*
      !                          t(i,delta->lamda,sigma)]
      !
      !  1/2 factor removed by lamda > sigma

      if ( a == gaama ) then
      if ( lamda > sigma ) then
         ti1 = TINDEX(i,delta+NCMO,lamda,sigma)
         if ( ti1 > 0 ) then
            Rnr_Sum = Rnr_Sum - TOLD(ti1)
         else
            ti1 = IABS(ti1)
            Rnr_Sum = Rnr_Sum + TOLD(ti1)
         endif
      endif
      endif

      ! All diagram done
      ! Now we do the multiplication
      ! Renormalization terms come with
      ! a negative sign, so make WINT negative.

      TNEW(eqn) = TNEW(eqn) + Rnr_Sum*(-WINT)

      ! Now Disconnected diagrams are summed.
      ! The procedure is standard. TmpR holds
      ! the contribution by Singles part.
      ! This is to be multiplied by
      ! [TOLD(j->b,vac)-TOLD(j->b,tvac)] = TmpT
      ! TOLD(j->b,vac) is stored in TSING.

      TmpR = -(Rnr_Sum*(-WINT))
      do j = 1, NHMO
      do b = 1, NPMO
         eqn = TINDEX(i,j,a,b)
         ti1 = TINDEX(j,b)
         TmpT = TSING(ti1) - TOLD(ti1)
         if ( eqn > 0 ) then
            TNEW(eqn) = TNEW(eqn) + TmpT*TmpR
         else
            eqn = IABS(eqn)
            TNEW(eqn) = TNEW(eqn) - TmpT*TmpR
         endif
      enddo
      enddo
   enddo
   enddo

   return

end subroutine Rnr_Sing_HEFF_TwoBody

subroutine Rnr_Dble_HEFF_OneBody(lamda,gaama,WINT)

   implicit none
   integer, intent(in) :: lamda, gaama
   double precision, intent(in) :: WINT
   double precision :: Rnr_Sum, TmpR, TmpT
   integer :: eqn, ti1, ti2
   integer :: i, a, j, b

   do i = 1, NHMO
   do j = 1, NHMO
   do a = 1, NPMO
   do b = 1, NPMO

      Rnr_Sum = dble(0)
      eqn = TINDEX(i,j,a,b)

      ! Diagram 01
      ! TERM : W1T2
      !
      ! + sum(gaama)[w(lamda->gaama)*t(gaama,j->a,b)
      !   (a,b) are canonical

      if ( i == lamda+NCMO ) then
      if ( a > b ) then
         ti1 = TINDEX(gaama+NCMO,j,a,b)
         if ( ti1 > 0 ) then
            Rnr_Sum = Rnr_Sum + TOLD(ti1)
         else
            ti1 = IABS(ti1)
            Rnr_Sum = Rnr_Sum - TOLD(ti1)
         endif
      endif
      endif

      ! Diagram 02
      ! TERM : W1T2
      !
      ! - sum(lamda)[w(lamda->gaama)*t(i,j->a,lamda)
      !   (i,j) are canonical

      if ( b == gaama ) then
      if ( i > j ) then
         ti1 = TINDEX(i,j,a,lamda)
         if ( ti1 > 0 ) then
            Rnr_Sum = Rnr_Sum - TOLD(ti1)
         else
            ti1 = IABS(ti1)
            Rnr_Sum = Rnr_Sum + TOLD(ti1)
         endif
      endif
      endif

      ! Diagram 03
      ! TERM : W1T1T2
      !
      ! -sum(lamda,gaama)[w(lamda->gaama)*t(i,lamda)*t(gaama,j->a,b)]
      !  (a,b) is canonical

      if ( a > b ) then
         ti1 = TINDEX(i,lamda)
         ti2 = TINDEX(gaama+NCMO,j,a,b)
         if ( ti2 > 0 ) then
            Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
         else
            ti2 = IABS(ti2)
            Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
         endif
      endif

      ! Diagram 04
      ! TERM : W1T1T2
      !
      ! -sum(lamda,gaama)[w(lamda->gaama)*t(gaama,a)*t(i,j->lamda,b)]
      !  (i,j) are canonical

      if ( i > j ) then
         ti1 = TINDEX(gaama+NCMO,a)
         ti2 = TINDEX(i,j,lamda,b)
         if ( ti2 > 0 ) then
            Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
         else
            ti2 = IABS(ti2)
            Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
         endif
      endif

      ! All diagram done
      ! Now we do the multiplication
      ! Renormalization terms come with
      ! a negative sign, so make WINT negative.

      if ( eqn > 0 ) then
         TNEW(eqn) = TNEW(eqn) + Rnr_Sum*(-WINT)
      else
         eqn = IABS(eqn)
         TNEW(eqn) = TNEW(eqn) - Rnr_Sum*(-WINT)
      endif

   enddo
   enddo
   enddo
   enddo

   return

end subroutine Rnr_Dble_HEFF_OneBody

subroutine Rnr_Dble_HEFF_TwoBody(lamda,sigma,gaama,delta,WINT)

   implicit none
   integer, intent(in) :: lamda, sigma, gaama, delta
   double precision, intent(in) :: WINT
   double precision :: Rnr_Sum, TmpR, TmpT
   integer :: eqn, ti1, ti2, ti3, ti4
   integer :: i, a, j, b

   do i = 1, NHMO
   do j = 1, NHMO
   do a = 1, NPMO
   do b = 1, NPMO

      Rnr_Sum = dble(0)
      eqn = TINDEX(i,j,a,b)

      ! TERM : W2T1
      !    Diagram 01
      !
      ! -sum(sigma)[w(lamda,sigma->gaama,delta)*t(j,sigma)]
      !  (gaama,delta) are canonical

      if ( gaama > delta ) then
      if ( a == gaama ) then
      if ( b == delta ) then
      if ( i == lamda+NCMO ) then
         ti1 = TINDEX(j,sigma)
         Rnr_Sum = Rnr_Sum - TOLD(ti1)
      endif
      endif
      endif
      endif

      ! TERM : W2T1
      !    Diagram 02
      !
      ! + sum(gaama)[w(lamda,sigma->gaama,delta)*t(gaama,a)
      !   (lamda,sigma) are canonical

      if ( lamda > sigma ) then
      if ( i == lamda+NCMO ) then
      if ( j == sigma+NCMO ) then
      if ( b == delta ) then
         ti1 = TINDEX(gaama+NCMO,a)
         Rnr_Sum = Rnr_Sum + TOLD(ti1)
      endif
      endif
      endif
      endif

      ! TERM : W2T1T1
      !    Diagram 01
      !
      ! -sum(sigma,gaama)[w(lamda,sigma->gaama,delta)*
      !                   t(j->sigma)*t(gaama->a)]

      if ( i == lamda+NCMO ) then
      if ( b == delta ) then
         ti1 = TINDEX(j,sigma)
         ti2 = TINDEX(gaama+NCMO,a)
         Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
      endif
      endif

      ! TERM : W2T1T1
      !    Diagram 02
      !
      ! +sum(gaama,delta)*[w(lamda,sigma->gaama,delta)*
      !                    t(gaama->a)*t(delta->b)]
      !  (a,b) are canonical
      !  (lamda,sigma) are canonical

      if ( lamda > sigma ) then
      if ( i == lamda+NCMO ) then
      if ( j == sigma+NCMO ) then
      if ( a > b ) then
         ti1 = TINDEX(gaama+NCMO,a)
         ti2 = TINDEX(delta+NCMO,b)
         Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
      endif
      endif
      endif
      endif

      ! TERM : W2T1T1
      !    Diagram 03
      ! +sum(lamda,sigma)[w(lamda,sigma->gaama,delta)*
      !                   t(i->lamda)*t(j->sigma)]
      !
      ! (i,j) are canonical
      ! (gaama,delta) are canonical

      if ( gaama > delta ) then
      if ( a == gaama ) then
      if ( b == delta ) then
      if ( i > j ) then
         ti1 = TINDEX(i,lamda)
         ti2 = TINDEX(j,sigma)
         Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
      endif
      endif
      endif
      endif

      ! TERM : W2T1T1T1
      !    Diagram 01
      !
      ! -sum(sigma,gaama,delta)[w(lamda,sigma->gaama,delta)*
      !                         t(gaama->a)*t(delta->b)*
      !                         t(j->sigma)]
      !  (a,b) are canonical

      !if (.not.OPTION(QUAD)) then
      if ( a > b ) then
      if ( i == lamda+NCMO ) then
         ti1 = TINDEX(gaama+NCMO,a)
         ti2 = TINDEX(delta+NCMO,b)
         ti3 = TINDEX(j,sigma)
         Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)*TOLD(ti3)
      endif
      endif
      !endif

      ! TERM : W2T1T1T1
      !    Diagram 02
      !
      ! +sum(lamda,sigma,gaama)[w(lamda,sigma->gaama,delta)*
      !                         t(i->lamda)*t(j->sigma)*
      !                         t(gaama->a)]
      !  (i,j) are canonical

      !if (.not.OPTION(QUAD)) then
      if ( i > j ) then
      if ( b == delta ) then
         ti1 = TINDEX(i,lamda)
         ti2 = TINDEX(j,sigma)
         ti3 = TINDEX(gaama+NCMO,a)
         Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)*TOLD(ti3)
      endif
      endif
      !endif

      ! TERM : W2T1T1T1T1
      !    Diagram 01
      !
      ! +sum(lamda,sigma,gaama,delta)[w(lamda,sigma->gaama,delta)*
      !                               t(i->lamda)*t(j->sigma)*
      !                               t(gaama,a)*t(delta->b)]
      !  (i,j) are canonical
      !  (a,b) are canonical

      if ( i > j ) then
      if ( a > b ) then
         ti1 = TINDEX(i,lamda)
         ti2 = TINDEX(j,sigma)
         ti3 = TINDEX(gaama+NCMO,a)
         ti4 = TINDEX(delta+NCMO,b)
         Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)*TOLD(ti3)*TOLD(ti4)
      endif
      endif

      ! TERM : W2T2
      !    Diagram 01
      !
      ! -sum(sigma,gaama)[w(lamda,sigma->gaama,delta)*
      !                   t(gaama,j->a,sigma)]

      if ( i == lamda+NCMO ) then
      if ( b == delta ) then
         ti1 = TINDEX(gaama+NCMO,j,a,sigma)
         if ( ti1 > 0 ) then
            Rnr_Sum = Rnr_Sum - TOLD(ti1)
         else
            ti1 = IABS(ti1)
            Rnr_Sum = Rnr_Sum + TOLD(ti1)
         endif
      endif
      endif

      ! TERM : W2T2
      !    Diagram 02
      !
      ! + 1/2 sum(gaama,delta)[w(lamda,sigma->gaama,delta)*
      !                   t(gaama,delta->a,b)]
      !  (a,b) are canonical
      !  (lamda,sigma) are canonical
      !  (gamma,delta) canonical to remove 1/2 factor.

      if ( lamda > sigma ) then
      if ( gaama > delta ) then
      if ( i == lamda+NCMO ) then
      if ( j == sigma+NCMO ) then
      if ( a > b ) then
         ti1 = TINDEX(gaama+NCMO,delta+NCMO,a,b)
         if ( ti1 > 0 ) then
            Rnr_Sum = Rnr_Sum + TOLD(ti1)
         else
            ti1 = IABS(ti1)
            Rnr_Sum = Rnr_Sum - TOLD(ti1)
         endif
      endif
      endif
      endif
      endif
      endif

      ! TERM : W2T2
      !    Diagram 03
      !
      ! + 1/2 sum(lamda,sigma)[w(lamda,sigma->gaama,delta)*
      !                   t(i,j->lamda,sigma)]
      !  (i,j) are canonical
      !  (gamma,delta) are canonical
      !  (lamda,sigma) canonical to remove 1/2 factor.

      if ( gaama > delta ) then
      if ( lamda > sigma ) then
      if ( a == gaama ) then
      if ( b == delta ) then
      if ( i > j ) then
         ti1 = TINDEX(i,j,lamda,sigma)
         if ( ti1 > 0 ) then
            Rnr_Sum = Rnr_Sum + TOLD(ti1)
         else
            ti1 = IABS(ti1)
            Rnr_Sum = Rnr_Sum - TOLD(ti1)
         endif
      endif
      endif
      endif
      endif
      endif

      ! TERM : W2T2T2
      !    Diagram 01
      !
      ! +sum(lamda,sigma,gaama,delta)[w(lamda,sigma->gaama,delta)*
      !                               t(delta,i->b,lamda)*
      !                               t(gaama,j->a,sigma)]
      ! Either (i,j) or (a,b) are canonical, but not both.

      if ( i > j ) then
         ti1 = TINDEX(delta+NCMO,i,b,lamda)
         ti2 = TINDEX(gaama+NCMO,j,a,sigma)
         if ( ti1 > 0 ) then
            if ( ti2 > 0 ) then
               Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
            else
               ti2 = IABS(ti2)
               Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
            endif
         else
            ti1 = IABS(ti1)
            if ( ti2 > 0 ) then
               Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
            else
               ti2 = IABS(ti2)
               Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
            endif
         endif
      endif

      ! TERM : W2T2T2
      !    Diagram 02
      !
      ! + 1/4 sum(lamda,sigma,gaama,delta)[w(lamda,sigma->gaama,delta)*
      !                                    t(gaama,delta->a,b)*
      !                                    t(i,j->lamda,sigma)]
      !   (i,j) and (a,b) are canonical
      !   (gaama,delta) and (lamda,sigma) are canonical to remove 1/4.

      if ( lamda > sigma ) then
      if ( gaama > delta ) then
      if ( i > j ) then
      if ( a > b ) then
         ti1 = TINDEX(gaama+NCMO,delta+NCMO,a,b)
         ti2 = TINDEX(i,j,lamda,sigma)
         if ( ti1 > 0 ) then
            if ( ti2 > 0 ) then
               Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
            else
               ti2 = IABS(ti2)
               Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
            endif
         else
            ti1 = IABS(ti1)
            if ( ti2 > 0 ) then
               Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
            else
               ti2 = IABS(ti2)
               Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
            endif
         endif
      endif
      endif
      endif
      endif

      ! TERM : W2T2T2
      !    Diagram 03
      !
      ! - 1/2 sum(lamda,sigma,gaama,delta)[w(lamda,sigma->gaama,delta)*
      !                                    t(gaama,i->lamda,sigma)*
      !                                    t(delta,j->a,b)]
      !   (a,b) are canonical
      !   (lamda,sigma) are canonical to remove 1/2.

      if ( lamda > sigma ) then
      if ( a > b ) then
         ti1 = TINDEX(gaama+NCMO,i,lamda,sigma)
         ti2 = TINDEX(delta+NCMO,j,a,b)
         if ( ti1 > 0 ) then
            if ( ti2 > 0 ) then
               Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
            else
               ti2 = IABS(ti2)
               Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
            endif
         else
            ti1 = IABS(ti1)
            if ( ti2 > 0 ) then
               Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
            else
               ti2 = IABS(ti2)
               Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
            endif
         endif
      endif
      endif

      ! TERM : W2T2T2
      !    Diagram 04
      !
      ! - 1/2 sum(lamda,sigma,gaama,delta)[w(lamda,sigma->gaama,delta)*
      !                                    t(gaama,delta->lamda,a)*
      !                                    t(i,j->sigma,b)]
      !   (i,j) are canonical
      !   (gaama,delta) are canonical to remove 1/2.

      if ( gaama > delta ) then
      if ( i > j ) then
         ti1 = TINDEX(gaama+NCMO,delta+NCMO,lamda,a)
         ti2 = TINDEX(i,j,sigma,b)
         if ( ti1 > 0 ) then
            if ( ti2 > 0 ) then
               Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
            else
               ti2 = IABS(ti2)
               Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
            endif
         else
            ti1 = IABS(ti1)
            if ( ti2 > 0 ) then
               Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
            else
               ti2 = IABS(ti2)
               Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
            endif
         endif
      endif
      endif

      ! TERM : W2T1T2
      !    Diagram 01
      !
      ! - 1/2 sum(sigma,gaama,delta)*[w(lamda,sigma->gaama,delta)*
      !                               t(gaama,delta->a,b)*
      !                               t(j->sigma)]
      !   (a,b) are canonical
      !   (gaama,delta) are canonical to remove 1/2.

      if ( gaama > delta ) then
      if ( i == lamda+NCMO ) then
      if ( a > b ) then
         ti1 = TINDEX(gaama+NCMO,delta+NCMO,a,b)
         ti2 = TINDEX(j,sigma)
         if ( ti1 > 0 ) then
            Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
         else
            ti1 = IABS(ti1)
            Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
         endif
      endif
      endif
      endif

      ! TERM : W2T1T2
      !    Diagram 02
      !
      ! + 1/2 sum(lamda,sigma,gaama)*[w(lamda,sigma->gaama,delta)*
      !                               t(i,j->lamda,sigma)*
      !                               t(gaama->a)]
      !   (i,j) are canonical
      !   (lamda,sigma) are canonical to remove 1/2.

      if ( lamda > sigma ) then
      if ( b == delta ) then
      if ( i > j ) then
         ti1 = TINDEX(i,j,lamda,sigma)
         ti2 = TINDEX(gaama+NCMO,a)
         if ( ti1 > 0 ) then
            Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
         else
            ti1 = IABS(ti1)
            Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
         endif
      endif
      endif
      endif

      ! TERM : W2T1T2
      !    Diagram 03
      !
      ! - sum(sigma,gaama,delta)[w(lamda,sigma->gaama,delta)*
      !                          t(delta->b)*t(gaama,j->a,sigma)]

      if ( i == lamda+NCMO ) then
         ti1 = TINDEX(delta+NCMO,b)
         ti2 = TINDEX(gaama+NCMO,j,a,sigma)
         if ( ti2 > 0 ) then
            Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
         else
            ti2 = IABS(ti2)
            Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
         endif
      endif

      ! TERM : W2T1T2
      !    Diagram 04
      !
      ! + sum(lamda,sigma,gaama)[w(lamda,sigma->gaama,delta)*
      !                          t(i->lamda)*t(gaama,j->a,sigma)]

      if ( b == delta ) then
         ti1 = TINDEX(i,lamda)
         ti2 = TINDEX(gaama+NCMO,j,a,sigma)
         if ( ti2 > 0 ) then
            Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)
         else
            ti2 = IABS(ti2)
            Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)
         endif
      endif

      ! TERM : W2T1T1T2
      !    Diagram 01
      !
      ! + 1/2 sum(lamda,sigma,gaama,delta)[w(lamda,sigma->gaama,delta)*
      !                                    t(i->lamda)*t(j->sigma)*
      !                                    t(gaama,delta->a,b)]
      !   (a,b) and (i,j) are canonical
      !   (gaama,delta) are canonical to remove 1/2.

      !if (.not.OPTION(QUAD)) then
      if ( gaama > delta ) then
      if ( a > b ) then
      if ( i > j ) then
         ti1 = TINDEX(i,lamda)
         ti2 = TINDEX(j,sigma)
         ti3 = TINDEX(gaama+NCMO,delta+NCMO,a,b)
         if ( ti3 > 0 ) then
            Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)*TOLD(ti3)
         else
            ti3 = IABS(ti3)
            Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)*TOLD(ti3)
         endif
      endif
      endif
      endif
      !endif

      ! TERM : W2T1T1T2
      !    Diagram 02
      !
      ! + 1/2 sum(lamda,sigma,gaama,delta)[w(lamda,sigma->gaama,delta)*
      !                                    t(gaama->a)*t(delta->b)*
      !                                    t(i,j->lamda,sigma)]
      !   (a,b) and (i,j) are canonical
      !   (lamda,sigma) are canonical to remove 1/2.

      !if (.not.OPTION(QUAD)) then
      if ( lamda > sigma ) then
      if ( i > j ) then
      if ( a > b ) then
         ti1 = TINDEX(gaama+NCMO,a)
         ti2 = TINDEX(delta+NCMO,b)
         ti3 = TINDEX(i,j,lamda,sigma)
         if ( ti3 > 0 ) then
            Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)*TOLD(ti3)
         else
            ti3 = IABS(ti3)
            Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)*TOLD(ti3)
         endif
      endif
      endif
      endif
      !endif

      ! TERM : W2T1T1T2
      !    Diagram 03
      !
      ! + sum(lamda,sigma,gaama,delta)[w(lamda,sigma->gaama,delta)*
      !                                  t(i->lamda)*t(delta->b)*
      !                                  t(gaama,j->a,sigma)]
      !if (.not.OPTION(QUAD)) then
         ti1 = TINDEX(i,lamda)
         ti2 = TINDEX(delta+NCMO,b)
         ti3 = TINDEX(gaama+NCMO,j,a,sigma)
         if ( ti3 > 0 ) then
            Rnr_Sum = Rnr_Sum + TOLD(ti1)*TOLD(ti2)*TOLD(ti3)
         else
            ti3 = IABS(ti3)
            Rnr_Sum = Rnr_Sum - TOLD(ti1)*TOLD(ti2)*TOLD(ti3)
         endif
      !endif

      ! All diagram done
      ! Now we do the multiplication
      ! Renormalization terms come with
      ! a negative sign, so make WINT negative.

      if ( eqn > 0 ) then
         TNEW(eqn) = TNEW(eqn) + Rnr_Sum*(-WINT)
      else
         eqn = IABS(eqn)
         TNEW(eqn) = TNEW(eqn) - Rnr_Sum*(-WINT)
      endif

   enddo
   enddo
   enddo
   enddo

   return

end subroutine Rnr_Dble_HEFF_TwoBody

end module RenormDiags
