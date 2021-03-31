!  Collection of simple utility functions used in exacorr

        real function relative_error(quantity_1,quantity_2)
        ! Gives percentage error for values > 1, otherwise absolute error

         implicit none
         real(8),intent(in):: quantity_1,quantity_2
         real(8)           :: error

         error = quantity_1 - quantity_2

         if (abs(quantity_1).gt.1.0) then
            relative_error = abs (error / quantity_1)
         else
            relative_error = abs (error)
         end if

        end function relative_error
           
        interface
            procedure print_date
        contains
        subroutine print_date(text,output)
        ! Prints time and date with a string
         implicit none
         character(len=*), intent(in)  :: text
         logical,intent(in), optional  :: output

         character(len=10) :: datex
         character(len=8)  :: timex
         integer           :: lentext, unit
         logical           :: output_local

         lentext = len(text)
         call daytime(datex,timex)
         if (present(output)) then
            output_local = output
         else
            output_local = .false.
         endif
         if (output_local) then
             write (*,'(1x,a10,1x,a8,1x,a)') datex,timex,text(1:lentext)
         else
             call get_free_fileunit(unit)
             open (unit,file='EXACORR.log',position='APPEND',form='FORMATTED')
             write (unit,'(1x,a10,1x,a8,1x,a)') datex,timex,text(1:lentext)
             close (unit,status='KEEP')
         end if
         
        end subroutine print_date
        end interface 

        subroutine get_free_fileunit(unit)
        ! idea from http://fortranwiki.org, will become obsolete once fortran2008 is supported everywhere
        ! gets free file unit
          integer, intent(out) :: unit
          integer, parameter   :: LUN_MIN=100, LUN_MAX=1000
          logical :: opened
          integer :: lun
          unit=-1
          do lun=LUN_MIN,LUN_MAX
            inquire(unit=lun,opened=opened)
            if (.not. opened) then
              unit=lun
              exit
            end if
          end do
        end subroutine get_free_fileunit

