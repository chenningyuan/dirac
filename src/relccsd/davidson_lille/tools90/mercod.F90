module mercod
!
!*******************************************************************************
!     ----- this module sets the run-time error code of the -----
!     ----- compiler                                        -----
!===============================================================================
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
!     ----- parameters -----
!
   integer, parameter :: MICODE=-2,MACODE=600,iw=6
   private :: iw
!
!     ----- local variables -----
!
   character (len=120), save, dimension (MICODE:MACODE) :: txtcod
   contains
   subroutine setcod
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!===============================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
      implicit none
!
!     ----- arguments -----
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
      call subent('SETCOD')
#endif
      txtcod = ' '
#ifdef _INTEL_
      txtcod(11) = 'Unit not connected'
      txtcod(29) = 'File not found'
      txtcod(34) = 'Unit already open'
      txtcod(41) = 'Insufficient virtual memory'
      txtcod(151) = 'Allocatable array is already allocated'
      txtcod(153) = 'Allocatable array or pointer is not allocated'
      txtcod(173) = 'A pointer passed to DEALLOCATE points to an array that&
                   & cannot be deallocated'
      txtcod(179) = 'Cannot allocate array - overflow on array size calculation'
      txtcod(585) = 'Array not allocated'
#endif
#ifdef _G95_64_
      txtcod(-2) = "End of record"
      txtcod(-1) = "End of file"
      txtcod(0) = "Successful return"
      txtcod(1:199) = ' '
      txtcod(151) = 'Allocatable array is already allocated'
      txtcod(200) = "Conflicting statement options"
      txtcod(201) = "Bad statement option"
      txtcod(202) = "Missing statement option"
      txtcod(203) = "File already opened in another unit"
      txtcod(204) = "Unattached unit"
      txtcod(205) = "FORMAT error"
      txtcod(206) = "Incorrect ACTION specified"
      txtcod(207) = "Read past ENDFILE record"
      txtcod(208) = "Bad value during read"
      txtcod(209) = "Numeric overflow on read"
      txtcod(210) = "Out of memory"
      txtcod(211) = "Array already allocated"
      txtcod(212) = "Deallocated a bad pointer"
      txtcod(214) = "Corrupt record in unformatted sequential-access file"
      txtcod(215) = "Reading more data than the record size (RECL)"
      txtcod(216) = "Writing more data than the record size (RECL)"
#endif
#ifdef _G95_32_
      txtcod(-2) = "End of record"
      txtcod(-1) = "End of file"
      txtcod(0) = "Successful return"
      txtcod(1:199) = ' '
      txtcod(151) = 'Allocatable array is already allocated'
      txtcod(200) = "Conflicting statement options"
      txtcod(201) = "Bad statement option"
      txtcod(202) = "Missing statement option"
      txtcod(203) = "File already opened in another unit"
      txtcod(204) = "Unattached unit"
      txtcod(205) = "FORMAT error"
      txtcod(206) = "Incorrect ACTION specified"
      txtcod(207) = "Read past ENDFILE record"
      txtcod(208) = "Bad value during read"
      txtcod(209) = "Numeric overflow on read"
      txtcod(210) = "Out of memory"
      txtcod(211) = "Array already allocated"
      txtcod(212) = "Deallocated a bad pointer"
      txtcod(214) = "Corrupt record in unformatted sequential-access file"
      txtcod(215) = "Reading more data than the record size (RECL)"
      txtcod(216) = "Writing more data than the record size (RECL)"
#endif
#ifdef _GNU_32_
      txtcod(0) = 'this error is not known by the program.&
                   & Consult the list of run-time errors of the compiler'
#endif
#ifdef _GNU_64_
      txtcod(0) = 'this error is not known by the program.&
                   & Consult the list of run-time errors of the compiler'
#endif
#ifdef _PGF_
      txtcod(0) = 'this error is not known by the program.&
                   & Consult the list of run-time errors of the compiler'
#endif
#ifdef _XLF_
      txtcod(0) = 'this error is not known by the program.&
                   & Consult the list of run-time errors of the compiler'
#endif
#ifdef _MYDEBUG_
      call subret('SETCOD',1)
#endif
      return
   end subroutine setcod
   subroutine pricod(icode)
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!===============================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      integer, intent(in) :: icode
!
!     ----- formats -----
!
    9999 format(/,5x,a)
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
      call subent('PRICOD')
#endif
      if(icode.lt.MICODE.or.icode.gt.MACODE) then
         write(iw,9999) 'error code unknown to the program'
      else
         write(iw,9999) trim(txtcod(icode))
      endif
#ifdef _MYDEBUG_
      call subret('PRICOD',1)
#endif
      return
   end subroutine pricod
end module mercod
