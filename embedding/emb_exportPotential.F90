!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

module FDE_exportPotentials

   use memory_allocator
   use fde_types

#ifdef MOD_UNRELEASED 

   use fde_cfg
   use fde_data
   use fde_io
   use xml_parser
   use xml_file
   use fde_eval_dirac
   
   private
   
   public fde_export_write_data_to_file
   public fde_export_write_data_to_memory

   contains


!   ----------------------------------------------------------------------------
      subroutine fde_export_write_data_to_memory(level)
!   ----------------------------------------------------------------------------
         character(len=3) :: level
         real(kind=8), allocatable :: dmat(:)

         call fde_test_active_density
         call fde_old_density_matrix_generator(dmat)
         call fde_test_active_density
         call fde_calculate_active_elpot(dmat,gf_active%elpot)
         call dealloc(dmat)
      end subroutine fde_export_write_data_to_memory


!   ----------------------------------------------------------------------------
      subroutine fde_export_write_data_to_file(level)
!   ----------------------------------------------------------------------------
         character(len=3) :: level
         character(len=4) :: export_format
         type(fde_export) :: etmp

         call fde_get_export_info(etmp)
         
         if (etmp%do_grid_out .and. (level.eq.etmp%ex_level)) then
! aspg: the one below has to be reactivated, now with the old
!       setup for calculating the density matrix its not neded
!           call fde_prepare_export(etmp)
            call fde_prepare_export(etmp%ex_level)
         
            call fde_get_export_info(etmp) 
            if (etmp%ex_format.eq.'XML') call fde_export_data_as_xml
         endif
      end subroutine fde_export_write_data_to_file


!   ----------------------------------------------------------------------------
!     subroutine fde_prepare_export(export_control)
      subroutine fde_prepare_export(level)
!   ----------------------------------------------------------------------------
!       type(fde_export) :: export_control
        character(len=3) :: level
        real(kind=8), allocatable :: dmat(:)
        integer :: mat_dim

         call fde_test_active_density

!! calculate vc, density and gradient
!         call fde_get_density_matrix(dmat,mat_dim,export_control%ex_level)
!         call fde_get_active_density(dmat, mat_dim, .true.)
!
! not so nice solution: while i don't get around to writing the ones above,
! just reuse the old code in the routine below
         call fde_old_density_matrix_generator(dmat)
         
         call fde_test_active_density
         call fde_calculate_active_elpot(dmat,fde_active_elpot)

         call dealloc(dmat)

      end subroutine fde_prepare_export


!   ----------------------------------------------------------------------------
      subroutine fde_export_data_as_xml
!   ----------------------------------------------------------------------------
         integer :: export_unit, log_unit
         character(len=60)  :: export_name, log_name
         character(len=132) :: csizept1
         character(len=132) :: csizept2
         character(len=132) :: csizegrdrho1
!not handling hessian now
!        character(len=132) :: csizehesrho1
         type(xml_tag),pointer :: tag
         real(kind=8), pointer :: gp(:,:), np(:), gnp(:,:), vcp(:)
         type(fde_files)  :: ftmp
         type(fde_export) :: etmp


         if (allocated(fde_grid)) then                
            gp  => fde_grid
         else
            call quit('grid not available in fde export')
         endif

         if (allocated(fde_active_n)) then
            np  => fde_active_n
         else
            call quit('density not available in fde export')
         endif
         
         if (allocated(fde_active_gn)) then
            gnp => fde_active_gn
         else
            call quit('density gradient not available in fde export')
         endif   
         
         if (allocated(fde_active_elpot)) then
            vcp => fde_active_elpot
         else
            call quit('potential not available in fde export')
         endif
         
         call fde_get_files_info(ftmp)
         call fde_get_export_info(etmp)
      
         export_unit = ftmp%export%unit
         export_name = ftmp%export%name

         log_unit    = ftmp%logfile%unit
         log_name    = ftmp%logfile%name
         

         CALL TITLER('Output of density to XML file','*',118)
         WRITE(log_unit,*) 'Name of the outputfile:', trim(export_name)
         WRITE(log_unit,*) 'Writing data to xml_file:', export_name

         write (csizept1,*)     size(gp,1)        ! ; trim (csizept1) 
         write (csizept2,*)     size(gp,2)        ! ; trim (csizept2) 
! not handling hessian now
!        write (csizehesrho1,*) size(fde_active_hn,1)   ! ; trim (csizehesrho1) 
         write (csizegrdrho1,*) size(gnp,1)  ! ; trim (csizegrdrho1) 

         tag=>xml_open(export_name)
           write(log_unit,*) 'tst'
           tag=>xml_tag_open(tag,'grid')
             call set_attribute(tag,'size',csizept2)
             tag=>xml_tag_open(tag,'dataset','gridpoints')
               call set_attribute(tag,'size',  csizept2)
               call set_attribute(tag,'width', csizept1)
               call xml_add_data(tag,gp)
             tag=>xml_tag_close(tag)

             if (etmp%ex_n) then
                tag=>xml_tag_open(tag,'dataset','density')
                  call set_attribute(tag,'size',  csizept2)
                  call set_attribute(tag,'width', '1')
                  call xml_add_data(tag,np)
                tag=>xml_tag_close(tag)
             endif

             if (etmp%ex_gn) then
                tag=>xml_tag_open(tag,'dataset','gradient')
                  call set_attribute(tag,'size', csizept2)
                  call set_attribute(tag,'width',csizegrdrho1)
                  call xml_add_data(tag,gnp)
                tag=>xml_tag_close(tag)
             endif

!
! while it's possible to do it, we don't export the hessian at the
! moment since there's no use for it...
!
!            tag=>xml_tag_open(tag,'dataset','hessian')
!              call set_attribute(tag,'size',  csizept2)
!              call set_attribute(tag,'width', csizehesrho1)    
!              call xml_add_data(tag,hes_rho)
!            tag=>xml_tag_close(tag)

             if (etmp%ex_coulomb) then
                tag=>xml_tag_open(tag,'dataset','vc')
                  call set_attribute(tag,'size',  csizept2)
                  call set_attribute(tag,'width', '1')
                  call xml_add_data(tag,vcp)
                tag=>xml_tag_close(tag)
             endif

           tag=>xml_tag_close(tag)

           write(log_unit,*) 'WRITING:'       
           call xml_write(export_name)
           write(log_unit,*) 'CLOSING:'
         tag=>xml_close(export_name)
         write(*,*) 'CLOSED!'

!aspg come back to this later         CALL TIME()

      end subroutine fde_export_data_as_xml

#endif /* MOD_UNRELEASED */

end module FDE_exportPotentials
