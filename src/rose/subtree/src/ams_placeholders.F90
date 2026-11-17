module CoreKF

contains

!
!    Placeholders for routines of CoreKF in case we can not link to ams.
!
 
    logical function kfFileExists(filename)
        character(*) :: filename

        print*, "kfFileExists called but AMS is not available"
        kfFileExists = .false.
    end function kfFileExists

    integer function kfOpenFile(filename)
        character(*) :: filename

        print*, "kfOpenFile called but AMS is not available"
        kfOpenfile = -1
    end function kfOpenFile

    integer function kfCloseFile(iunit)
        integer      :: iunit

        print*, "kfCloseFile called but AMS is not available"
        kfClosefile = -1
    end function kfCloseFile

    integer function kfReadInt(iunit,text, array, n, i)
        integer      :: iunit,n,i
        character(*) :: text
        integer      :: array(*)

        print*, "kfReadInt called but AMS is not available"
        kfReadInt = -1
    endfunction  kfReadInt

    integer function kfReadReal(iunit,text, array, n, i)
        integer      :: iunit,n,i
        character(*) :: text
        real(8)      :: array(*)

        print*, "kfReadReal called but AMS is not available"
        kfReadReal = -1
    end function kfReadReal

    integer function kfWriteInt(iunit,text, array, n, i)
        integer      :: iunit,n,i
        character(*) :: text
        integer      :: array(*)

        print*, "kfWriteInt called but AMS is not available"
        kfWriteInt = -1
    endfunction  kfWriteInt

    integer function kfWriteReal(iunit,text, array, n, i)
        integer      :: iunit,n,i
        character(*) :: text
        real(8)      :: array(*)

        print*, "kfWriteReal called but AMS is not available"
        kfWriteReal = -1
    end function kfWriteReal

    integer function kfQuerySection(iunit, section)
        integer,       intent(in) :: iunit
        character(*),  intent(in) :: section

        print*, "kfQuerySection called but AMS is not available"
        kfQuerySection = -1
    end function kfQuerySection


    integer function kfDeleteSection(iunit, section)
         integer,       intent(in) :: iunit
         character(*),  intent(in) :: section
 
         print*, "kfDeleteSection called but AMS is not available"
         kfDeleteSection = -1
    end function kfDeleteSection
    
     integer function kfCreateSection (iunit, section)
        integer,       intent(in) :: iunit
        character(*),  intent(in) :: section

        print*, "kfCreateSection called but AMS is not available"
        kfCreateSection = -1
      end function kfCreateSection
      
end module CoreKF
