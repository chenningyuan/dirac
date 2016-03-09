module rdmfil
!
!***********************************************************************
!    ----- a series of generic routines to open, close, read and -----
!    ----- a direct access file                                  -----
!=======================================================================
!  usage:
!       call openms(iunit,irndm,nd,iflag) to open file iunit.
!                                         only iunit is used, and
!                                         iflag must be eqal to 0
!
!       call writms(iunit,a,n,irec,iflag) to write N data of array A
!                                         with flag IREC.
!            If iflag = 0  the data is written to the end of the file
!                     = -1 the data is written at the same place if there
!                          is room enough, otherwise at the end of the file
!                     = 1 the data is written at the same place if there
!                          is room enough, otherwise stop
!
!       call readms(iunit,a,n,irec) to read N data of array A with flag
!                                   IREC.
!
!         array A may be be of type integer, real, double precision or
!         complex*16
!
!    the record length is in bytes, to be remembered when compiling
!-----------------------------------------------------------------------
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
!
!     ----- parameters -----
!
   use kindvar
   use mercod
   implicit none
   integer (kind=K_I4), parameter :: I_ONE=1
   integer (kind=K_I4), parameter :: NRECL=10,MASK=z'ffff'
   integer (kind=K_I4), parameter :: LRECL=2**NRECL
   integer (kind=K_I4), parameter :: NBI4=LRECL*2
   integer (kind=K_I4), parameter :: NBR4=LRECL*2
   integer (kind=K_I4), parameter :: NBR8=LRECL
   integer (kind=K_I4), parameter :: NBC8=LRECL/2
   integer (kind=K_I4), parameter :: NBS1=LRECL*8
!
!     ----- LRECL = 1024                                    -----
!     ----- record length = LRECL*8 = 8096 bytes            -----
!     ----- for a variable occupying NBYTE bytes, there are -----
!     ----- LRECL*8/NBYTE of them in a record               -----
!
!     ----- global variables -----
!
  integer (kind=K_I4), private, save :: last,mxrec
  integer (kind=K_I4), save, allocatable, dimension (:) :: irndm
!
 interface readms
    module procedure rdms_i4,rdms_r4,rdms_r8,rdms_c8,rdms_ch1
 end interface
 interface writms
    module procedure wrms_i4, wrms_r4,wrms_r8,wrms_c8,wrms_ch1
 end interface
   contains
   subroutine openms(iunit,jrmdm,nd,iflag)
!
!***********************************************************************
!    ---- open direct acces unit iunit -----
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer (kind=K_I4), intent(in) :: iunit,nd,iflag
   integer (kind=K_I4), intent(inout), dimension (:) :: jrmdm
!
!     ----- local variables -----
!
   integer (kind=K_I4) :: ios
!  save last
!
!     +++++
!     start
!     +++++
!
      ios = 0
      if(.not.allocated(irndm)) allocate(irndm(nd),stat=ios)
      if(ios.ne.0) then
         write(6,'(5x,a)') 'OPENMS failed to allocate table of records'
         write(6,"(5x,'IOS=',i4)") ios
         call pricod(ios)
         stop
      endif
      open(unit=iunit,access='direct',recl=LRECL*8,status='scratch',  &
           form = 'unformatted')
      mxrec = nd
      if(iflag.ne.0) then
         read(iunit,rec=1) last,irndm
      else
         irndm = -1
         last = 1
      endif
      return
   end subroutine openms
   subroutine  rdms_i4(iunit,a,n,irec)
!
!***********************************************************************
!***********************************************************************
!
   implicit none
!
!     ----- parameters -----
!
      integer, parameter :: NBYTE=4
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: iunit,n,irec
      integer (kind=K_I4), intent(out), dimension (:) :: a
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: nrec,nrest,krec,i,i1,i2
!
!     save irndm,last
!
!     +++++
!     start
!     +++++
!
      nrec = n/NBI4
      nrest = n-nrec*NBI4
!
      krec = ishft(irndm(irec),-16)
      i2 = 0
      do i=1,nrec
         i1 = i2+I_ONE
         i2 = i2+NBI4
         read(unit=iunit,rec=krec+i) a(i1:i2)
      enddo
      if(nrest.ne.0) then
         i1 = i2+I_ONE
         i2 = i2+nrest
         read(unit=iunit,rec=krec+nrec+I_ONE) a(i1:i2)
      endif
      return
   end subroutine rdms_i4
   subroutine  rdms_r4(iunit,a,n,irec)
!
!***********************************************************************
!***********************************************************************
!
   implicit none
!
!     ----- parameters -----
!
      integer, parameter :: NBYTE=4
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: iunit,n,irec
      real, intent(out), dimension (:) :: a
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: nrec,nrest,krec,i,i1,i2
!
!     save irndm,last
!
!     +++++
!     start
!     +++++
!
      nrec = n/NBR4
      nrest = n-nrec*NBR4
!
      krec = ishft(irndm(irec),-16)
      i2 = 0
      do i=1,nrec
         i1 = i2+I_ONE
         i2 = i2+NBR4
         read(unit=iunit,rec=krec+i) a(i1:i2)
      enddo
      if(nrest.ne.0) then
         i1 = i2+I_ONE
         i2 = i2+nrest
         read(unit=iunit,rec=krec+nrec+I_ONE) a(i1:i2)
      endif
      return
   end subroutine rdms_r4
   subroutine  rdms_r8(iunit,a,n,irec)
!
!***********************************************************************
!***********************************************************************
!
   implicit none
!
!     ----- parameters -----
!
      integer, parameter :: NBYTE=8
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: iunit,n,irec
      double precision, intent(out), dimension (:) :: a
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: nrec,nrest,krec,i,i1,i2
!
!     save irndm,last
!
!     +++++
!     start
!     +++++
!
      nrec = n/NBR8
      nrest = n-nrec*NBR8
!
      krec = ishft(irndm(irec),-16)
      i2 = 0
      do i=1,nrec
         i1 = i2+I_ONE
         i2 = i2+NBR8
         read(unit=iunit,rec=krec+i) a(i1:i2)
      enddo
      if(nrest.ne.0) then
         i1 = i2+I_ONE
         i2 = i2+nrest
         read(unit=iunit,rec=krec+nrec+I_ONE) a(i1:i2)
      endif
      return
   end subroutine rdms_r8
   subroutine  rdms_c8(iunit,a,n,irec)
!
!***********************************************************************
!***********************************************************************
!
   implicit none
!
!     ----- parameters -----
!
      integer, parameter :: NBYTE=16
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: iunit,n,irec
      complex (kind=K_C8), intent(out), dimension (:) :: a
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: nrec,nrest,krec,i,i1,i2
!
!     save irndm,last
!
!     +++++
!     start
!     +++++
!
      nrec = n/NBC8
      nrest = n-nrec*NBC8
!
      krec = ishft(irndm(irec),-16)
      i2 = 0
      do i=1,nrec
         i1 = i2+I_ONE
         i2 = i2+NBC8
         read(unit=iunit,rec=krec+i) a(i1:i2)
      enddo
      if(nrest.ne.0) then
         i1 = i2+I_ONE
         i2 = i2+nrest
         read(unit=iunit,rec=krec+nrec+I_ONE) a(i1:i2)
      endif
      return
   end subroutine rdms_c8
   subroutine wrms_i4(iunit,a,n,irec,iflag)
!
!***********************************************************************
!     iflag  =  0   on ecrit l'enregistrement en fin de fichier
!            = -1   on ecrit l'enregistrement au meme endroit s'il y a
!                 la place,sinon on l'ecrit a la fin.
!            =  1   on reecrit en place. si la place manque... stop.
!***********************************************************************
!
      implicit none
!
!     ----- parameters -----
!
      integer, parameter :: NBYTE=4
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: iunit,n,irec,iflag
      integer (kind=K_I4), intent(in), dimension (:) :: a
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: nrec,nrest,mrec,krec,i,i1,i2
!     save last
!
!     +++++
!     start
!     +++++
!
      if(irec.gt.size(irndm)) then
         write(6,"(5x,2(a,i4))") 'WRITMS error: IREC=',irec,  &
                    ' exceeds the maximum number of possible data:',size(irndm)
         stop
      endif
      nrec = n/NBI4
      nrest = n-nrec*NBI4
      if(nrest.ne.0) then
         mrec = nrec+I_ONE
      else
         mrec = nrec
      endif
      if((iflag.eq.1).and.(iand(irndm(irec),MASK).lt.mrec)) then
         stop '***  stop  ***  error in WRITMS... record too long ***'
      else if((iflag.lt.0).and.(iand(irndm(irec),MASK).lt.mrec)  &
             .or.(iflag.eq.0).or.(irndm(irec).lt.0)) then
         irndm(irec) = ior(ishft(last,16),mrec)
         last = last+mrec
      endif
!
      krec = ishft(irndm(irec),-16)
      i2 = 0
      do i=1,nrec
         i1 = i2+I_ONE
         i2 = i2+NBI4
         write(unit=iunit,rec=krec+i) a(i1:i2)
      enddo
      if(nrest.ne.0) then
         i1 = i2+I_ONE
         i2 = i2+nrest
         write(unit=iunit,rec=krec+mrec) a(i1:i2)
      endif
      return
   end subroutine wrms_i4
   subroutine wrms_r4(iunit,a,n,irec,iflag)
!
!***********************************************************************
!     iflag  =  0   on ecrit l'enregistrement en fin de fichier
!            = -1   on ecrit l'enregistrement au meme endroit s'il y a
!                 la place,sinon on l'ecrit a la fin.
!            =  1   on reecrit en place. si la place manque... stop.
!***********************************************************************
!
      implicit none
!
!     ----- parameters -----
!
      integer, parameter :: NBYTE=4
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: iunit,n,irec,iflag
      real, intent(in), dimension (:) :: a
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: nrec,nrest,mrec,krec,i,i1,i2
!     save last
!
!     +++++
!     start
!     +++++
!
      if(irec.gt.size(irndm)) then
         write(6,"(5x,2(a,i4))") 'WRITMS error: IREC=',irec,  &
                    ' exceeds the maximum number of possible data:',size(irndm)
         stop
      endif
      nrec = n/NBR4
      nrest = n-nrec*NBR4
      if(nrest.ne.0) then
         mrec = nrec+I_ONE
      else
         mrec = nrec
      endif
      if((iflag.eq.1).and.(iand(irndm(irec),MASK).lt.mrec)) then
         stop '***  stop  ***  error in WRITMS... record too long ***'
      else if((iflag.lt.0).and.(iand(irndm(irec),MASK).lt.mrec)  &
               .or.(iflag.eq.0).or.(irndm(irec).lt.0)) then
         irndm(irec) = ior(ishft(last,16),mrec)
         last = last+mrec
      endif
!
      krec = ishft(irndm(irec),-16)
      i2 = 0
      do i=1,nrec
         i1 = i2+I_ONE
         i2 = i2+NBR4
         write(unit=iunit,rec=krec+i) a(i1:i2)
      enddo
      if(nrest.ne.0) then
         i1 = i2+I_ONE
         i2 = i2+nrest
         write(unit=iunit,rec=krec+mrec) a(i1:i2)
      endif
      return
   end subroutine wrms_r4
   subroutine wrms_r8(iunit,a,n,irec,iflag)
!
!***********************************************************************
!     iflag  =  0   on ecrit l'enregistrement en fin de fichier
!            = -1   on ecrit l'enregistrement au meme endroit s'il y a
!                 la place,sinon on l'ecrit a la fin.
!            =  1   on reecrit en place. si la place manque... stop.
!***********************************************************************
!
      implicit none
!
!     ----- parameters -----
!
      integer, parameter :: NBYTE=8
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: iunit,n,irec,iflag
      double precision, intent(in), dimension (:) :: a
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: nrec,nrest,mrec,krec,i,i1,i2
!     save last
!
!     +++++
!     start
!     +++++
!
      if(irec.gt.size(irndm)) then
         write(6,"(5x,2(a,i4))") 'WRITMS error: IREC=',irec,  &
                    ' exceeds the maximum number of possible data:',size(irndm)
         stop
      endif
      nrec = n/NBR8
      nrest = n-nrec*NBR8
      if(nrest.ne.0) then
         mrec = nrec+I_ONE
      else
         mrec = nrec
      endif
      if((iflag.eq.1).and.(iand(irndm(irec),MASK).lt.mrec)) then
         stop '***  stop  ***  error in WRITMS... record too long ***'
      else if((iflag.lt.0).and.(iand(irndm(irec),MASK).lt.mrec)  &
               .or.(iflag.eq.0).or.(irndm(irec).lt.0)) then
         irndm(irec) = ior(ishft(last,16),mrec)
         last = last+mrec
      endif
!
      krec = ishft(irndm(irec),-16)
      i2 = 0
      do i=1,nrec
         i1 = i2+I_ONE
         i2 = i2+NBR8
         write(unit=iunit,rec=krec+i) a(i1:i2)
      enddo
      if(nrest.ne.0) then
         i1 = i2+I_ONE
         i2 = i2+nrest
         write(unit=iunit,rec=krec+mrec) a(i1:i2)
      endif
      return
   end subroutine wrms_r8
   subroutine wrms_c8(iunit,a,n,irec,iflag)
!
!***********************************************************************
!     iflag  =  0   on ecrit l'enregistrement en fin de fichier
!            = -1   on ecrit l'enregistrement au meme endroit s'il y a
!                 la place,sinon on l'ecrit a la fin.
!            =  1   on reecrit en place. si la place manque... stop.
!***********************************************************************
!
      implicit none
!
!     ----- parameters -----
!
      integer, parameter :: NBYTE=16
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: iunit,n,irec,iflag
      complex (kind=K_C8), intent(in), dimension (:) :: a
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: nrec,nrest,mrec,krec,i,i1,i2
!     save last
!
!     +++++
!     start
!     +++++
!
      if(irec.gt.size(irndm)) then
         write(6,"(5x,2(a,i4))") 'WRITMS error: IREC=',irec,  &
                    ' exceeds the maximum number of possible data:',size(irndm)
         stop
      endif
      nrec = n/NBC8
      nrest = n-nrec*NBC8
      if(nrest.ne.0) then
         mrec = nrec+I_ONE
      else
         mrec = nrec
      endif
      if((iflag.eq.1).and.(iand(irndm(irec),MASK).lt.mrec)) then
         stop '***  stop  ***  error in WRITMS... record too long ***'
      else if((iflag.lt.0).and.(iand(irndm(irec),MASK).lt.mrec)  &
               .or.(iflag.eq.0).or.(irndm(irec).lt.0)) then
         irndm(irec) = ior(ishft(last,16),mrec)
         last = last+mrec
      endif
!
      krec = ishft(irndm(irec),-16)
      i2 = 0
      do i=1,nrec
         i1 = i2+I_ONE
         i2 = i2+NBC8
         write(unit=iunit,rec=krec+i) a(i1:i2)
      enddo
      if(nrest.ne.0) then
         i1 = i2+I_ONE
         i2 = i2+nrest
         write(unit=iunit,rec=krec+mrec) a(i1:i2)
      endif
      return
   end subroutine wrms_c8
   subroutine  rdms_ch1(iunit,a,n,irec)
!
!***********************************************************************
!     ----- read data of type character*1 -----
!***********************************************************************
!
   implicit none
!
!     ----- parameters -----
!
      integer, parameter :: NBYTE=1
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: iunit,n,irec
      character (len=1), intent(out), dimension (:) :: a
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: nrec,nrest,krec,i,i1,i2
!
!     save irndm,last
!
!     +++++
!     start
!     +++++
!
      nrec = n/NBS1
      nrest = n-nrec*NBS1
!
      krec = ishft(irndm(irec),-16)
      i2 = 0
      do i=1,nrec
         i1 = i2+I_ONE
         i2 = i2+NBS1
         read(unit=iunit,rec=krec+i) a(i1:i2)
      enddo
      if(nrest.ne.0) then
         i1 = i2+I_ONE
         i2 = i2+nrest
         read(unit=iunit,rec=krec+nrec+1) a(i1:i2)
      endif
      return
   end subroutine rdms_ch1
   subroutine wrms_ch1(iunit,a,n,irec,iflag)
!
!***********************************************************************
!     iflag  =  0   on ecrit l'enregistrement en fin de fichier
!            = -1   on ecrit l'enregistrement au meme endroit s'il y a
!                 la place,sinon on l'ecrit a la fin.
!            =  1   on reecrit en place. si la place manque... stop.
!***********************************************************************
!
      implicit none
!
!     ----- parameters -----
!
      integer, parameter :: NBYTE=1
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: iunit,n,irec,iflag
      character (len=1), intent(in), dimension (:) :: a
!
!     ----- local variables -----
!
      integer (kind=K_I4) :: nrec,nrest,mrec,krec,i,i1,i2
!     save last
!
!     +++++
!     start
!     +++++
!
      if(irec.gt.size(irndm)) then
         write(6,"(5x,2(a,i4))") 'WRITMS error: IREC=',irec,  &
                    ' exceeds the maximum number of possible data:',size(irndm)
         stop
      endif
      nrec = n/NBS1
      nrest = n-nrec*NBS1
      if(nrest.ne.0) then
         mrec = nrec+I_ONE
      else
         mrec = nrec
      endif
      if((iflag.eq.1).and.(iand(irndm(irec),MASK).lt.mrec)) then
         stop '***  stop  ***  error in WRITMS... record too long ***'
      else if((iflag.lt.0).and.(iand(irndm(irec),MASK).lt.mrec)  &
               .or.(iflag.eq.0).or.(irndm(irec).lt.0)) then
         irndm(irec) = ior(ishft(last,16),mrec)
         last = last+mrec
      endif
!
      krec = ishft(irndm(irec),-16)
      i2 = 0
      do i=1,nrec
         i1 = i2+I_ONE
         i2 = i2+NBS1
         write(unit=iunit,rec=krec+i) a(i1:i2)
      enddo
      if(nrest.ne.0) then
         i1 = i2+I_ONE
         i2 = i2+nrest
         write(unit=iunit,rec=krec+mrec) a(i1:i2)
      endif
      return
   end subroutine wrms_ch1
   subroutine closms(iunit)
!
!***********************************************************************
!    fermeture du fichier iunit. on ecrit irndm sur le premier bloc
!***********************************************************************
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: iunit
!
!     +++++
!     start
!     +++++
!
      write(unit=iunit,rec=1) last,irndm
      close(unit=iunit)
      deallocate(irndm)
      return
   end subroutine closms
end module rdmfil
