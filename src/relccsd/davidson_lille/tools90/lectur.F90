module lectur
  use caps
!
!***********************************************************************
!          ----- Module for reading data in free format -----
!=======================================================================
!     contains:
!         subroutine inpa   read a character string
!                    inpi   read an integer
!                    inpe   read a real*4
!                    inpf   read a real*8
!                    inpl   read a logical*4
!                    inpnum return the number of data on a the line
!                    inpprt print the current line (for debug)
!                    fetchl get a line of data
!                    lecerr (private, flush routine in case of error)
!         function   rval   (private, used by inpe, inpf)
!-----------------------------------------------------------------------
!   usage:
!          first a call to fetchl(unit,ier) to get a line of data
!                       return IER=0 if succesful.
!                       UNIT (I4) is the logical unit where to read
!          then calls to inpx(var) where x = a i e f l
!               return a value of the given type in var
!
!       the data on the line are sperated by blanks, missing data
!       at the end of the line are set to the "NULL" value of the
!       type:
!             0     for integers
!            0.0    for reals
!           blank   for strings
!          .FALSE.  for logicals
!   note * "NULL" values in the middle of the list must be explicitly
!                    typed. A blank string may not be read in the middle
!                    use any thing, like the string ...
!        * strings may not contain imbedded blanks
!        * real values may be given in F E or D format, or as rational
!          numbers as in 3/4 (ratio of 2 INTEGERS)
!-----------------------------------------------------------------------
!    exemple:
!            integer i,ier
!            real*8 val
!            character*8 string
!
!            call fetchl(5,ier)
!            call inpa(string)
!            call inpi(i)
!            call inpf(val)
!
!     the line of data looks like
!
!              data  6 3.876
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
  implicit none
!
!   ---- parameters -----
!
  character (len=1) , parameter :: BLANK1=' '
  character (len=13), parameter :: jchar = '0123456789+-.'
  integer , parameter :: ALEN=80, MAXDAT=40
  double precision , parameter :: ZERO=0.0D+00,ONE=1.0D+00,  &
                                  TEN=10.0D+00,PT1=0.1D+00
  private:: ZERO,ONE,TEN,PT1,BLANK1
!
!   ----- global variables -----
!
  character(len=ALEN) :: ia
  integer ::  jrec,jump
  integer , dimension (MAXDAT) :: istrt, inumb
  private  :: ia, jrec,jump, istrt, inumb
!
  contains
    function rval(cuf) result(value)
!
!***********************************************************************
!          ----- read an integer*4 for rational numbers -----
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character(len=*) , intent(in) :: cuf
      double precision  :: value
!
!     ----- local variables -----
!
      character (len=1) :: itemp
      integer :: n,ifact,nstrt,j,jbuf
!
!     ----- format -----
!
 9999 format(3x,'error while reading ',a12)
!
!     +++++
!     start
!     +++++
!
      n = len(cuf)
      ifact = 1
      jbuf = 0
      do nstrt=n,1,-1
         itemp = cuf(nstrt:nstrt)
         j = index(jchar,itemp)
         if(j.eq.0) then
            write(6,9999) cuf
            stop
         endif
         if(j.ge.11) then
            if(nstrt.ne.1) then
               write(6,9999) cuf
               stop
            endif
            if(j.eq.12) then
               jbuf = -jbuf
               value = dble(jbuf)
               return
            endif
         else
            jbuf = jbuf+(j-1)*ifact
            ifact = ifact*10
         endif
      enddo
      value = dble(jbuf)
      return
    end function rval
    subroutine inpa(guf)
!
!***********************************************************************
!                   ----- read a string -----
!***********************************************************************
!
       implicit none
!
!   ----- arguments -----
!
       character(len=*) , intent(out) :: guf
!
!   ----- local variables -----
!
       integer :: ll,i,n,nstrt,ifin
!
!   +++++
!   start
!   +++++
!
       ll = len(guf)
       do i=1,ll
          guf(i:i) = BLANK1
       enddo
       jrec = jrec+1
       if(jrec.gt.jump) then
          return
       else
          n = min(inumb(jrec),ll)
          nstrt = istrt(jrec)
          ifin = nstrt+n-1
          guf(1:n) = ia(nstrt:ifin)
       endif
       return
    end subroutine inpa

    subroutine inpe(buf)
!
!***********************************************************************
!                   ----- read a real*8 -----
!***********************************************************************
!
  implicit none
!
!     ----- arguments -----
!
  real, intent(out) :: buf
!
!     ----- local variables -----
!
  character (len=1) :: itemp
  integer :: iexp1,ifact,jexp,nstrt,n,ist,iend,i,if2,ifact2, &
             islash,is2,j,limit,if1
  double precision :: expo,fact,dum
!
!   +++++
!   start
!   +++++
!
  jrec = jrec+1
  if(jrec.gt.jump) then
     buf = ZERO
     return
  endif
  limit = 13
  ist = istrt(jrec)
  n = inumb(jrec)
  iend = ist+n-1
!
!     ----- look for an exponent -----
!
  iexp1 = index(ia(ist:iend),'d')
  if(iexp1.eq.0) iexp1 = index(ia(ist:iend),'D')
  if(iexp1.eq.0) iexp1 = index(ia(ist:iend),'e')
  if(iexp1.eq.0) iexp1 = index(ia(ist:iend),'E')
  if(iexp1.ne.0) then
     jexp = 0
     ifact = 1
     do nstrt=iend,ist+iexp1,-1
        itemp = ia(nstrt:nstrt)
        j = index(jchar(1:12),itemp)
        if(j.eq.0) call lecerr(nstrt,1)
        if(j.ge.11) then
           if(nstrt.ne.ist+iexp1) call lecerr(nstrt,1)
           if(j.eq.12) jexp = -jexp
        else
           jexp = jexp+(j-1)*ifact
           ifact = ifact*10
        endif
     enddo
     expo = TEN**jexp
     iend = iexp1-2+ist
  else
     expo = ONE
  endif
!
  ifact2 = 0
  fact = ONE
  i = 0
  dum = ZERO
  if2 = iend
  islash = index(ia(ist:if2),'/')
  if(islash.ne.0) then
     if1 = islash+ist-2
     is2 = ist+islash
     buf = real(rval(ia(ist:if1))/rval(ia(is2:if2)))
  else
     do nstrt=iend,ist,-1
        i = i+1
        itemp = ia(nstrt:nstrt)
        j = index(jchar(1:limit),itemp)
        if(j.eq.0) call lecerr(nstrt,1)
        if(j.ge.11) then
           if(j.eq.13) then
              ifact2 = i-1
              limit = 12
           else
              if(nstrt.ne.ist) call lecerr(nstrt,1)
              if(j.eq.12) dum = -dum
              dum = (PT1**ifact2)*dum*expo
              buf = real(dum)
              return
           endif
        else
           dum = dum+dble(j-1)*fact
           fact = fact*TEN
        endif
     enddo
     dum = (PT1**ifact2)*dum*expo
     buf = real(dum)
  endif
  return
end subroutine inpe
    subroutine inpf(buf)
!
!***********************************************************************
!                   ----- read a real*8 -----
!***********************************************************************
!
  implicit none
!
!     ----- arguments -----
!
  double precision , intent(out) :: buf
!
!     ----- local variables -----
!
  character (len=1) :: itemp
  integer :: iexp1,ifact,jexp,nstrt,n,ist,iend,i,if2,ifact2, &
             islash,is2,j,limit,if1
  double precision :: expo,fact
!
!   +++++
!   start
!   +++++
!
  jrec = jrec+1
  if(jrec.gt.jump) then
     buf = ZERO
     return
  endif
  limit = 13
  ist = istrt(jrec)
  n = inumb(jrec)
  iend = ist+n-1
!
!     ----- look for an exponent -----
!
  iexp1 = index(ia(ist:iend),'d')
  if(iexp1.eq.0) iexp1 = index(ia(ist:iend),'D')
  if(iexp1.eq.0) iexp1 = index(ia(ist:iend),'e')
  if(iexp1.eq.0) iexp1 = index(ia(ist:iend),'E')
  if(iexp1.ne.0) then
     jexp = 0
     ifact = 1
     do nstrt=iend,ist+iexp1,-1
        itemp = ia(nstrt:nstrt)
        j = index(jchar(1:12),itemp)
        if(j.eq.0) call lecerr(nstrt,1)
        if(j.ge.11) then
           if(nstrt.ne.ist+iexp1) call lecerr(nstrt,1)
           if(j.ge.12) jexp = -jexp
        else
           jexp = jexp+(j-1)*ifact
           ifact = ifact*10
        endif
     enddo
     expo = TEN**jexp
     iend = iexp1-2+ist
  else
     expo = ONE
  endif
!
  ifact2 = 0
  fact = ONE
  i = 0
  buf = ZERO
  if2 = iend
  islash = index(ia(ist:if2),'/')
  if(islash.ne.0) then
     if1 = islash+ist-2
     is2 = ist+islash
     buf = rval(ia(ist:if1))/rval(ia(is2:if2))
  else
     do nstrt=iend,ist,-1
        i = i+1
        itemp = ia(nstrt:nstrt)
        j = index(jchar(1:limit),itemp)
        if(j.eq.0) call lecerr(nstrt,1)
        if(j.ge.11) then
           if(j.eq.13) then
              ifact2 = i-1
              limit = 12
           else
              if(nstrt.ne.ist) call lecerr(nstrt,1)
              if(j.eq.12) buf = -buf
              buf = (PT1**ifact2)*buf*expo
              return
           endif
        else
           buf = buf+dble(j-1)*fact
           fact = fact*TEN
        endif
     enddo
     buf = (PT1**ifact2)*buf*expo
  endif
  return
end subroutine inpf
subroutine inpi(jbuf)
!
!***********************************************************************
!                 ----- read an integer*4 -----
!***********************************************************************
!
   implicit none
!
!    ----- arguments -----
!
    integer , intent(out) :: jbuf
!
!    ----- local variables -----
!
    character (len=1) :: itemp
    integer n,ifact,ist,j,nstrt
!
!    +++++
!    start
!    +++++
!
      jrec = jrec+1
      if(jrec.gt.jump) then
         jbuf = 0
      else
         n = inumb(jrec)
         ifact = 1
         ist = istrt(jrec)
         jbuf = 0
         do nstrt=ist+n-1,ist,-1
            itemp = ia(nstrt:nstrt)
            j = index(jchar,itemp)
            if(j.eq.0) call lecerr(nstrt,2)
            if(j.ge.11) then
               if(nstrt.ne.ist) call lecerr(nstrt,2)
               if(j.eq.12) then
                  jbuf = -jbuf
                  return
               endif
            else
               jbuf = jbuf+(j-1)*ifact
               ifact = ifact*10
            endif
         enddo
      endif
      return
      end subroutine inpi
subroutine inpl(jbuf)
!
!***********************************************************************
!                   ----- read a logical*4 -----
!***********************************************************************
!
    implicit none
!
!     ----- arguments -----
!
      logical , intent(out) :: jbuf
!
!     ----- local variables -----
!
      character (len=1) itemp
      integer :: ist
!
!     +++++
!     start
!     +++++
!
      jrec = jrec+1
      if(jrec.gt.jump) then
         jbuf = .FALSE.
      else
         ist = istrt(jrec)
         itemp = ia(ist:ist+1)
         call capson(itemp)
         if(itemp.eq.'T') then
            jbuf = .TRUE.
         else if(itemp.eq.'F') then
            jbuf = .FALSE.
         else
            call lecerr(ist,3)
         endif
      endif
      return
      end subroutine inpl
subroutine inpnum(num)
!
!***********************************************************************
!         ----- returns the  number of data of the line -----
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer , intent(out) :: num
!
!     +++++
!     start
!     +++++
!
   num=jump
   return
end subroutine inpnum
subroutine fetchl(iunit,ier)
!
!***********************************************************************
!     ----- read a line of data as a character string -----
!***********************************************************************
!
   implicit none
!
!   ----- arguments -----
!
   integer , intent(in) :: iunit
   integer , intent(out) :: ier
!
!  ----- local variables -----
!
   integer :: ios,i
   logical :: isw
!
   jump = 0
   jrec = 0
   isw = .TRUE.
   read(iunit,fmt='(a)',iostat=ios) ia
   if(ios.lt.0) then
      ier = 1
      return
   endif
   ier = 0
   do i=1,ALEN
      if(ia(i:i).ne.BLANK1) then
         if(isw) then
            jump = jump+1
            istrt(jump) = i
            inumb(jump) = 1
            isw = .FALSE.
         else
            inumb(jump) = inumb(jump)+1
         endif
      else
         isw = .TRUE.
      endif
   enddo
   return
  end subroutine fetchl
  subroutine lecerr(m,ier)
!
!***********************************************************************
!               ----- error routine. Kick you off -----
!***********************************************************************
!
    implicit none
!
!     ----- parameters -----
!
    character(len=1), parameter :: HYPHEN='-', STAR='*'
!
!     ----- arguments -----
!
    integer , intent (in) :: m,ier
!
!     ----- local variables -----
!
    character(len=ALEN) :: ib
    character(len=7), dimension(3) :: typev=(/ 'real   ','integer','logical' /)
    integer :: i
!
!     ----- formats -----
!
 9999 format(/,3x,' === error encountered while reading a ',a7, &
             //,3x,a80,/,3x,a80)
!
!     +++++
!     start
!     +++++
!
       do i=1,m-1
          ib(i:i) = HYPHEN
       enddo
       ib(m:m) = STAR
       do i=m+1,ALEN
          ib(i:i) = BLANK1
       enddo
       write(6,9999) typev(ier),ia,ib
       stop
 end subroutine lecerr
 subroutine inpprt(iw)
!
!***********************************************************************
!     lecture d'une carte de donnees comme un chaine de caracteres
!***********************************************************************
!
!     ----- arguments -----
!
       integer, intent(in) :: iw     ! unite d'ecriture
!
  100 format(2x,' LINE=',a80)
  110 format(2x,' NREC=',i4)
!
       write(iw,100) ia
       write(iw,110) jump
       return
  end subroutine inpprt
 subroutine inpdup(iw)
!
!***********************************************************************
!     lecture d'une carte de donnees comme un chaine de caracteres
!***********************************************************************
!
!     ----- arguments -----
!
       integer, intent(in) :: iw     ! unite d'ecriture
!
  100 format(a80)
!
       write(iw,100) ia
       return
  end subroutine inpdup
  subroutine getcom(icom)
!
!***********************************************************************
!     lecture d'une carte de donnees comme un chaine de caracteres
!***********************************************************************
!
!     ----- arguments -----
!
   character (len=*), intent(out) :: icom
!
   icom = ia
   return
  end subroutine getcom
end module lectur
