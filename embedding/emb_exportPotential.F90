!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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
   use fde_types
   use fde_cfg
   use fde_data
   use fde_io
   use fde_evaluators
   use xml_parser
   use xml_file

   private
   
   public fde_export_write_data_to_file

   contains

      subroutine fde_export_write_data_to_file(level)
         character(len=3) :: level
         character(len=4) :: export_format
         type(fde_export) :: etmp

         call fde_get_export_info(etmp)
         
         if (etmp%do_grid_out .and. (level.eq.etmp%ex_level)) then
            call fde_prepare_export
         
            call fde_get_export_info(etmp) 
            if (etmp%ex_format.eq.'XML') call fde_export_data_as_xml
         endif
      end subroutine fde_export_write_data_to_file


      subroutine fde_prepare_export
         real(kind=8), allocatable :: dmat(:,:,:)

! calculate vc, density and gradient
         call fde_get_density_matrix(dmat)

! calculate vc, density and gradient, will be stored at the fde_data module
         call fde_calculate_active_n_gn(dmat)
         call fde_calculate_active_elpot(dmat)

         call fde_cleanup_density_matrix(dmat)
      end subroutine fde_prepare_export

      subroutine fde_export_data_as_xml
         integer :: file_unit
         character(len=60)  :: file_name
         character(len=132) :: csizept1
         character(len=132) :: csizept2
         character(len=132) :: csizegrdrho1
!not handling hessian now
!        character(len=132) :: csizehesrho1
         type(xml_tag),pointer :: tag
         real(kind=8), pointer :: gp(:,:)
         type(fde_files)  :: ftmp
         type(fde_export) :: etmp

         call set_pointer_fde_grid(gp)

         call fde_get_files_info(ftmp)
         call fde_get_export_info(etmp)
      
         file_unit = ftmp%export%unit
         file_name = ftmp%export%name

!         call fde_get_log_unit(lu)
!         call fde_get_name(export_file)
!         call fde_get_what_to_export(export_n, export_gn, export_coulomb)

         CALL TITLER('Output of density to XML file','*',118)
         WRITE(lu,*) 'Name of the outputfile:', trim(file_name)
         WRITE(lu,*) 'Writing data to xml_file:', file_name

         write (csizept1,*)     size(gp,1)        ! ; trim (csizept1) 
         write (csizept2,*)     size(gp,2)        ! ; trim (csizept2) 
! not handling hessian now
!        write (csizehesrho1,*) size(fde_active_hn,1)   ! ; trim (csizehesrho1) 
         write (csizegrdrho1,*) size(fde_active_gn,1)  ! ; trim (csizegrdrho1) 

         tag=>xml_open(file_name)
           write(lu,*) 'tst'
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
                  call xml_add_data(tag,fde_active_n)
                tag=>xml_tag_close(tag)
             endif

             if (etmp%ex_gn) then
                tag=>xml_tag_open(tag,'dataset','gradient')
                  call set_attribute(tag,'size', csizept2)
                  call set_attribute(tag,'width',csizegrdrho1)
                  call xml_add_data(tag,fde_active_gn)
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
                  call xml_add_data(tag,fde_active_elpot)
                tag=>xml_tag_close(tag)
             endif

           tag=>xml_tag_close(tag)

           write(*,*) 'WRITING:'       
           call xml_write(file_name)
           write(*,*) 'CLOSING:'
         tag=>xml_close(file_name)
         write(*,*) 'CLOSED!'

!aspg come back to this later         CALL TIME()

      end subroutine fde_export_data_as_xml

end module FDE_exportPotentials
