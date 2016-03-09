module jpsort
!
!***********************************************************************
!     ----- generic module for sorting data in ascending order -----
!=======================================================================
!     contains:
!       subroutine gensort(n,crit,indix)    generic subroutine
!       subroutine sort_i4(n,crit,indix)    to sort integer data
!       subroutine sort_r8(n,crit,indix)    to sort real*8 data
!       subroutine sort_ch(n,crit,indix,a)  to sort chracter data
!-----------------------------------------------------------------------
!   usage:
!   ------
!   call gensort(n,crit,indix)  generic subroutine
!
!      N(I4) = number of data to be sorted (INPUT)
!   CRIT(??) = array of dimension N at least, containing the
!              data to be sorted (INPUT/OUTPUT))
!              *NOTE* character data are not sorted.
!  INDIX(I4) = array of dimension N at least N (INPUT/OUPUT)
!              on input contains indices
!              on exit contains the permutted indices
!
!   call gensort(n,crit,indix,a)  generic subroutine
!
!      N(I4) = number of data to be sorted (INPUT)
!   CRIT(??) = array of dimension N at least, containing the
!              data to be sorted (INPUT/OUTPUT))
!              *NOTE* character data are not sorted.
!  INDIX(I4) = array of dimension N at least N (INPUT/OUPUT)
!              on input contains indices
!              on exit contains the permutted indices
!      A(??) = a character variable of the same type (charcter*(X))
!              as crit (SCRATCH) for sorting of character strings
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
interface gensort
   module procedure sort_i4,sort_r8,sort_ch
end interface
contains
! subroutine sort_ch(n,crit,indix,a)
  subroutine sort_ch(n,crit,indix,a)
!
!**********************************************************************
!      quick-sort pour 1 tableau crit. on conserve un tableau d'index
!**********************************************************************
!
   implicit none
!
!     ----- parameters -----
!
   integer, parameter :: M=7,NSTACK=50
   double precision, parameter :: FM=7.875D+03,FA=2.11D+02,FC=1.663D+03,  &
                                  FMI=1.2698413D-04,ZERO=0.0D+00
!
!     ----- arguments ----
!
   character (len=*), intent(in), dimension (*) :: crit
   character (len=*), intent(inout) :: a
   integer, intent(in) :: n
   integer, intent(inout), dimension (*) :: indix
!
!     ----- local variables -----
!
   integer :: jstack,l,ir,ind,i,j,iq
   integer, dimension(NSTACK) :: istack
   double precision :: fx
!
!     ----- formats -----
!
 9999 format(//,5x,'******* QUICK-SORT erreur *******',  &
                   ' le tableau istack est trop petit',/,  &
                5x,' c''est fini pour aujourd''hui')
      save
!
!     +++++
!     start
!     +++++
!
   jstack = 0
   l = 1
   ir = n
   fx = ZERO
10 if(ir-l.lt.M) then
     do j=l+1,ir
       ind = indix(j)
       a = crit(ind)
       do i=j-1,1,-1
         if(crit(indix(i)).le.a) goto 12
!           crit(i+1) = crit(i)
         indix(i+1) = indix(i)
       enddo
       i = 0
!  12     crit(i+1) = a
12     indix(i+1) = ind
     enddo
     if(jstack.eq.0) return
     ir = istack(jstack)
     l = istack(jstack-1)
     jstack = jstack-2
   else
     i = l
     j = ir
     fx = mod(fx*FA+FC,FM)
     iq = l+(ir-l+1)*(fx*FMI)
     a = crit(indix(iq))
!       crit(iq) = crit(l)
     ind = indix(iq)
     indix(iq) = indix(l)
20   continue
21     if(j.gt.0) then
          if(a.lt.crit(indix(j))) then
             j = j-1
             goto 21
          endif
       endif
       if(j.le.i) then
!            crit(i) = a
          indix(i) = ind
          goto 30
       endif
!         crit(i) = crit(j)
       indix(i) = indix(j)
       i = i+1
22     if(i.le.n) then
          if(a.gt.crit(indix(i))) then
             i = i+1
             goto 22
          endif
       endif
       if(j.le.i) then
!            crit(j) = a
          indix(j) = ind
          i = j
          goto 30
      endif
!        crit(j) = crit(i)
      indix(j) = indix(i)
      j = j-1
      goto 20
30    jstack = jstack+2
      if(jstack.gt.NSTACK) then
         write(6,9999)
         stop
      endif
      if(ir-i.ge.i-l) then
         istack(jstack) = ir
         istack(jstack-1) = i+1
         ir = i-1
      else
         istack(jstack) = i-1
         istack(jstack-1) = l
         l = i+1
      endif
   endif
   goto 10
   end subroutine sort_ch
   subroutine sort_i4(n,crit,indix)
!
!**********************************************************************
!     ----- Quick-sort for the array crit. We keep the index table ----
!**********************************************************************
!
   implicit none
!
!     ----- parameters -----
!
   integer, parameter :: M=7,NSTACK=50
   double precision, parameter :: FM=7.875D+03,FA=2.11D+02,FC=1.663D+03,  &
                                  FMI=1.2698413D-04,ZERO=0.0D+00
!
!     ----- arguments ----
!
   integer, intent(in) :: n
   integer, intent(inout), dimension (*) :: crit
   integer, intent(inout), dimension (*) :: indix
!
!     ----- local variables -----
!
   integer :: jstack,l,ir,ind,i,j,iq,a
   integer, dimension(NSTACK) :: istack
   double precision :: fx
!
!     ----- formats -----
!
 9999 format(//,5x,'******* QUICK-SORT erreur *******',  &
                   ' le tableau istack est trop petit',/,  &
                5x,' c''est fini pour aujourd''hui')
      save
!
!     +++++
!     start
!     +++++
!
      jstack = 0
      l = 1
      ir = n
      fx = ZERO
   10 if(ir-l.lt.M) then
        do 13 j=l+1,ir
          a = crit(j)
          ind = indix(j)
          do 11 i=j-1,1,-1
            if(crit(i).le.a) goto 12
            crit(i+1) = crit(i)
            indix(i+1) = indix(i)
   11     continue
          i = 0
   12     crit(i+1) = a
          indix(i+1) = ind
   13   continue
        if(jstack.eq.0) return
        ir = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack-2
      else
        i = l
        j = ir
        fx = mod(fx*FA+FC,FM)
        iq = l+(ir-l+1)*(fx*FMI)
        a = crit(iq)
        crit(iq) = crit(l)
        ind = indix(iq)
        indix(iq) = indix(l)
   20   continue
   21     if(j.gt.0) then
             if(a.lt.crit(j)) then
                j = j-1
                goto 21
             endif
          endif
          if(j.le.i) then
             crit(i) = a
             indix(i) = ind
             goto 30
          endif
          crit(i) = crit(j)
          indix(i) = indix(j)
          i = i+1
   22     if(i.le.n) then
             if(a.gt.crit(i)) then
                i = i+1
                goto 22
             endif
          endif
          if(j.le.i) then
             crit(j) = a
             indix(j) = ind
             i = j
             goto 30
         endif
         crit(j) = crit(i)
         indix(j) = indix(i)
         j = j-1
         goto 20
   30    jstack = jstack+2
         if(jstack.gt.NSTACK) then
            write(6,9999)
            stop
         endif
         if(ir-i.ge.i-l) then
            istack(jstack) = ir
            istack(jstack-1) = i+1
            ir = i-1
         else
            istack(jstack) = i-1
            istack(jstack-1) = l
            l = i+1
         endif
      endif
      goto 10
      end subroutine sort_i4
      subroutine sort_r8(n,tab,indix)
!
!**********************************************************************
!      quick-sort pour 1 tableau tab. on conserve un tableau d'index
!**********************************************************************
!
   implicit none
!
!     ----- parameters -----
!
   integer, parameter :: M=7,NSTACK=50
   double precision, parameter :: FM=7.875D+03,FA=2.11D+02,FC=1.663D+03,  &
                                  FMI=1.2698413D-04,ZERO=0.0D+00
!
!     ----- arguments ----
!
   integer, intent(in) :: n
   integer, intent(inout), dimension (*) :: indix
   double precision, intent(inout), dimension (*) :: tab
!
!     ----- local variables -----
!
   integer :: jstack,l,ir,ind,i,j,iq
   integer, dimension(NSTACK) :: istack
   double precision :: fx,a
!
!     ----- formats -----
!
 9999 format(//,5x,'******* QUICK-SORT erreur *******',  &
                   ' le tableau istack est trop petit',/,  &
                5x,' c''est fini pour aujourd''hui')
      save
!
!     +++++
!     start
!     +++++
!
      jstack = 0
      l = 1
      ir = n
      fx = ZERO
   10 if(ir-l.lt.M) then
        do j=l+1,ir
          a = tab(j)
          ind = indix(j)
          do i=j-1,1,-1
            if(tab(i).le.a) goto 12
            tab(i+1) = tab(i)
            indix(i+1) = indix(i)
          enddo
          i = 0
   12     tab(i+1) = a
          indix(i+1) = ind
        enddo
        if(jstack.eq.0) return
        ir = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack-2
      else
        i = l
        j = ir
        fx = mod(fx*FA+FC,FM)
        iq = l+(ir-l+1)*(fx*FMI)
        a = tab(iq)
        tab(iq) = tab(l)
        ind = indix(iq)
        indix(iq) = indix(l)
   20   continue
   21     if(j.gt.0) then
             if(a.lt.tab(j)) then
                j = j-1
                goto 21
             endif
          endif
          if(j.le.i) then
             tab(i) = a
             indix(i) = ind
             goto 30
          endif
          tab(i) = tab(j)
          indix(i) = indix(j)
          i = i+1
   22     if(i.le.n) then
             if(a.gt.tab(i)) then
                i = i+1
                goto 22
             endif
          endif
          if(j.le.i) then
             tab(j) = a
             indix(j) = ind
             i = j
             goto 30
         endif
         tab(j) = tab(i)
         indix(j) = indix(i)
         j = j-1
         goto 20
   30    jstack = jstack+2
         if(jstack.gt.NSTACK) then
            write(6,9999)
            stop
         endif
         if(ir-i.ge.i-l) then
            istack(jstack) = ir
            istack(jstack-1) = i+1
            ir = i-1
         else
            istack(jstack) = i-1
            istack(jstack-1) = l
            l = i+1
         endif
      endif
      goto 10
      end subroutine sort_r8
end module jpsort
