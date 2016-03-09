module timing
!
  double precision, private :: ti,tx,tim
!
  contains
      subroutine timit(indix)
!
!***********************************************************************
!          ----- compute and write the elapsed CPU time -----
!***********************************************************************
!
      implicit none
!
!     ------ arguments -----
!
      integer, intent(in) :: indix
!
!     ------ local variables -----
!
      double precision tu,sec
!
!     ----- formats -----
!
 1000 format(/,10x,'elapsed time = ',f10.3,5x,'total time = ',f10.3,//)
!
!     +++++
!     start
!     +++++
!
      if(indix.lt.0) then
         call seconde(sec)
         tx = sec
         ti = tx
      else
         call seconde(sec)
         tim = sec-tx
         tu = tim-ti
         ti = tim
      endif
      if(indix.gt.0) write(6,1000) tu,tim
      return
      end subroutine timit
      subroutine seconde(sec)
!
!***********************************************************************
!                   ----- get the elapsed CPU time -----
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
      double precision, intent(out) :: sec
!
!     ----- local variables -----
!
      real :: ut
!
!     +++++
!     start
!     +++++
!
   call cpu_time(ut)
   sec = ut
   return
   end subroutine seconde
   subroutine dateur(jour,mois,an,heure,minute,secnde)
!
!*******************************************************************************
!               ----- find date and time of the call -----
!-------------------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
#ifdef _MYDEBUG_
     use mdebug
#endif
     implicit none
!
!     ----- arguments -----
!
     integer, intent(out) :: jour,mois,an,heure,minute,secnde
!
!     ----- local variables -----
!
     character (len=9) :: dada
     character (len=10) :: titi
     character (len=3), dimension (12) :: months=(/ 'Jan','Feb','Mar','Apr',  &
                                          'May','Jun','Jul','Aug','Sep','Oct', &
                                          'Nov','Dec' /)
!
 9999 format(/,5x,'Job started ',i2.2,' ',a3,' ',i2.2,'  at ',i2.2,':',  &
                i2.2,':',i2.2,//)
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
     call subent('DATEUR')
#endif
     call date_and_time(dada,titi)
     read(dada(3:4),'(i2)') an
     read(dada(5:6),'(i2)') mois
     read(dada(7:8),'(i2)') jour
     read(titi(1:2),'(i2)') heure
     read(titi(3:4),'(i2)') minute
     read(titi(5:6),'(i2)') secnde
     write(6,9999) jour,months(mois),an,heure,minute,secnde
#ifdef _MYDEBUG_
     call subret('DATEUR',1)
#endif
     return
   end subroutine dateur
end module timing
