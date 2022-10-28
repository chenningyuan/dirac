!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end
!
!
! file-I/O routines used in the X2C routines.
!
! written by sknecht april 2010
!
module x2c_fio

  implicit none

  public x2c_write
  public x2c_read

  private

contains

!**********************************************************************
  subroutine x2c_write(flabel,fmat,nlen,funit)
!**********************************************************************
     use checkpoint
     real(8),            intent(in)     :: fmat(*)
     integer,            intent(in)     :: nlen
     integer,            intent(in)     :: funit
     character (len=12), intent(in)     :: flabel
!----------------------------------------------------------------------
     logical                            :: fndlab12
     character (len=800)                :: funit_name
     integer                            :: lfunit_name
     character*(:), allocatable         :: group
!**********************************************************************

     ! write essential x2c information to CHECKPOINT
     if (flabel(1:8)=='prpint2c'.or.flabel(1:7)=='pctmtAO'.or.flabel(1:4)=='h12c') then
       group = '/result/hamiltonian/x2c/ao_matrices/'
       if (nlen>0) then
          call checkpoint_write(group//flabel,fmat(1:nlen))
       endif
       return
     end if

       if(nlen >= 0)then
         rewind funit                                                    
         if(fndlab12('EOFLABEL-x2c',funit))then
           backspace(funit)
           call newlab12(flabel,funit,6)
           call writt(funit,nlen,fmat)
           call newlab12('EOFLABEL-x2c',funit,6)
         else
           inquire(unit=funit, name=funit_name)
           lfunit_name = len_trim(funit_name)
           print '(/2x,a,a,a)', '*** error in x2c_write: end-of-file label is missing in the file ',            &
                              funit_name(1:lfunit_name),'. The program therefore stops. ***'                            
           call quit('*** error in x2c_write: end-of-file label is missing. ***')
         end if
       else 
         call newlab12('EOFLABEL-x2c',funit,6)
       end if

  end subroutine x2c_write

!**********************************************************************
  subroutine x2c_read(flabel,fmat,nlen,funit)
!**********************************************************************
     use checkpoint
     real(8),            intent(inout)  :: fmat(*)
     integer,            intent(in)     :: nlen
     integer,            intent(in)     :: funit
     character (len=12), intent(in)     :: flabel
!----------------------------------------------------------------------
     logical                            :: fndlab12
     character (len=800)                :: funit_name
     character (len=1000)               :: error_stream
     integer                            :: lfunit_name
     integer                            :: lerror_stream
     character*(:), allocatable         :: group
!**********************************************************************
 
     ! read essential x2c information from CHECKPOINT
     if (flabel(1:8)=='prpint2c'.or.flabel(1:7)=='pctmtAO'.or.flabel(1:4)=='h12c') then
       group = '/result/hamiltonian/x2c/ao_matrices/'
       if (nlen>0) then
          call checkpoint_read(group//flabel,fmat(1:nlen))
       endif
       return
     end if

       rewind(funit)
       if(nlen > 0)then
         if(fndlab12(flabel,funit))then
           call readt(funit,nlen,fmat)
         else
           inquire(unit=funit, name=funit_name)
           lfunit_name = len_trim(funit_name)
           print '(a,a,a,a)', '   *** error in x2c_read: data with label ',flabel,' is not present in the file ',  &
                                  funit_name(1:lfunit_name)
           write(error_stream,'(a,a,a,a)') & 
           '   *** error in x2c_read: data with label ',flabel,' is not present in the file ',funit_name(1:lfunit_name)
           lerror_stream = len_trim(error_stream)
           call quit(error_stream(1:lerror_stream))
         end if
       else
         print '(2x,a)', ' *** warning in x2c_read: attempt to read data array which is supposed to be zero. ***'
!        call quit('bla bla')
       end if
       
  end subroutine x2c_read
!**********************************************************************

end module x2c_fio
