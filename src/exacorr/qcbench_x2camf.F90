module qcbench_x2camf
! Canonical X2CAMF-CCSD/CCSD(T)/EOM-CCSD paper binaries and operand recording for qcbench.

        use tensor_algebra
        use exacorr_tensor_interface
        use exacorr_datatypes
        use exacorr_scaling_parameters
        use exacorr_mo
        use, intrinsic :: ISO_C_BINDING

        implicit none
        private

        public qcbench_pack_chol_ao
        public qcbench_pack_chol_mo
        public qcbench_pack_coefficients
        public qcbench_x2camf_ccsd_residual
        public qcbench_x2camf_ccsd_energy
        public qcbench_x2camf_ccsd_t_energy
        public qcbench_record_residual
        public qcbench_record_triples
        public qcbench_x2camf_eom_probe

        character(len=512) :: record_dir = ''
        integer :: record_calls = 0
        integer :: record_nocc = 0, record_nvir = 0, record_naux = 0, record_nao = 0
        character(len=64) :: record_case_ids(10000)
        character(len=64) :: record_call_sites(10000)
        character(len=16) :: record_routes(10000)
        character(len=8) :: record_coefficients(10000)
        character(len=16) :: record_components(10000)
        character(len=128) :: record_equations(10000)
        character(len=64) :: record_shapes(3,10000)

        complex(8), parameter :: Q_ZERO=(0.D0,0.D0), Q_ONE=(1.D0,0.D0)

       contains

        function coefficient(case_id) result(value)
          character(len=*), intent(in) :: case_id
          character(len=8) :: value
          value = '+1'
          select case (trim(case_id))
          case ('fmi.tau_oovv','wmnij.tau','t2.oooo','t2.ao_occ', &
                't2.ao_project_a','t2.foo_dress','response.t3','response.t4')
            value = '+1/2'
          case ('fmi.fme','fae.fme','fae.tau_oovv','t1.ooov','t2.fvv_dress', &
                'sigma1.i6','sigma1.i7','ao.t1_correction')
            value = '-1/2'
          case ('energy.tau_oovv','quartic.tau')
            value = '+1/4'
          case ('T.energy.w_v')
            value = '+1/36'
          case ('wmbej.t1_oo','wmbej.tau','wmbej.cd_outer_m','fae.cd_outer_o', &
                't1.foo','t1.cd_outer','t1.ovov','t2.foo','t2.last_bare','t2.ao_last', &
                'T.connected.cd_exchange_outer','T.connected.occupied', &
                'T.integral.bcjk_exchange','T.integral.jkma_exchange', &
                'i6.t1_oovv','i7.t1_oovv','sigma2.wovoo','sigma2.i3','sigma2.foo')
            value = '-1'
          case default
            value = '+1'
          end select
        end function coefficient

        subroutine recorder_start(outdir, nocc, nvir, naux, nao)
          character(len=*), intent(in) :: outdir
          integer, intent(in) :: nocc, nvir, naux, nao
          record_dir = trim(outdir)
          record_nocc = nocc
          record_nvir = nvir
          record_naux = naux
          record_nao = nao
          record_calls = 0
          call execute_command_line('mkdir -p "'//trim(record_dir)//'"')
          call write_manifest()
        end subroutine recorder_start

        subroutine write_json_array_int(values, unit)
          integer, intent(in) :: values(:), unit
          integer :: i
          write(unit,'(a)',advance='no') '['
          do i=1,size(values)
            if (i>1) write(unit,'(a)',advance='no') ','
            write(unit,'(i0)',advance='no') values(i)
          end do
          write(unit,'(a)',advance='no') ']'
        end subroutine write_json_array_int

        function shape_string(tensor) result(text)
          type(tensor_t), intent(inout) :: tensor
          character(len=64) :: text
          integer(C_INT) :: dims(6), ierr, nmodes, i
          text = '['
          ierr = get_nmodes(tensor,nmodes)
          if (nmodes>0) then
            ierr = get_extents_static(tensor,dims)
            do i=1,nmodes
              if (i>1) text = trim(text)//','
              write(text,'(a,i0)') trim(text)//' ',dims(i)
            end do
          end if
          text = trim(text)//' ]'
        end function shape_string

        subroutine write_manifest()
          integer :: unit, i, role
          if (len_trim(record_dir)==0) return
          open(newunit=unit,file=trim(record_dir)//'/manifest.json',status='replace',action='write')
          write(unit,'(a)') '{'
          write(unit,'(a)') '  "route": "DIRAC X2CAMF-CCSD/CCSD(T)/EOM-CCSD paper routes (Z2023 Eqs. 5-17 and Table 2)",'
          write(unit,'(a)') '  "corpus": "cd_x2camf_ccsd_corpus",'
          write(unit,'(a,i0,a)') '  "n_contractions": ',record_calls,','
          write(unit,'(a,i0,a,i0,a,i0,a,i0,a)') '  "dimensions": {"o": ',record_nocc,', "v": ',record_nvir, &
               ', "ao": ',record_nao,', "naux": ',record_naux,'},'
          write(unit,'(a)') '  "calls": ['
          do i=1,record_calls
            write(unit,'(a)') '    {'
            write(unit,'(a)') '      "case": "'//trim(record_case_ids(i))//'",'
            write(unit,'(a)') '      "corpus_case": "'//trim(record_case_ids(i))//'",'
            write(unit,'(a)') '      "call_site": "'//trim(record_call_sites(i))//'",'
            write(unit,'(a)') '      "route": "'//trim(record_routes(i))//'",'
            write(unit,'(a)') '      "equation": "'//trim(record_equations(i))//'",'
            write(unit,'(a)') '      "coefficient": "'//trim(record_coefficients(i))//'",'
            write(unit,'(a)') '      "component": "'//trim(record_components(i))//'",'
            write(unit,'(a,i0,a)') '      "n_call": ',i,','
            do role=1,3
              write(unit,'(a)') '      "'//achar(ichar('A')+role-1)//'": {'
              write(unit,'(a)') '        "file": "'//trim(record_case_ids(i))//write_suffix(i,role)//'",'
              write(unit,'(a)') '        "dtype": "complex128",'
              write(unit,'(a)') '        "shape": '//trim(record_shapes(role,i))
              if (role<3) then
                write(unit,'(a)') '      },'
              else
                write(unit,'(a)') '      }'
              end if
            end do
            if (i<record_calls) then
              write(unit,'(a)') '    },'
            else
              write(unit,'(a)') '    }'
            end if
          end do
          write(unit,'(a)') '  ]'
          write(unit,'(a)') '}'
          close(unit)
        end subroutine write_manifest

        function write_suffix(call_index, role) result(suffix)
          integer, intent(in) :: call_index, role
          character(len=16) :: suffix
          write(suffix,'(a,i6.6,a)') '__',call_index,'__'//achar(ichar('A')+role-1)//'.c8'
        end function write_suffix

        subroutine record_array(tensor, case_id, role, call_index)
          type(tensor_t), intent(inout) :: tensor
          character(len=*), intent(in) :: case_id, role
          integer, intent(in) :: call_index
          integer :: unit, n, ierr
          integer(8) :: numel
          type(C_PTR) :: body_p
          complex(8), pointer, contiguous :: data(:)
          character(len=32) :: suffix
          suffix = '__'//trim(role)//'.c8'
          n = get_numel(tensor,numel)
          ierr = direct_data_pointer(tensor,body_p)
          call c_f_pointer(body_p,data,(/numel/))
          suffix = write_suffix(call_index,ichar(role)-ichar('A')+1)
          open(newunit=unit,file=trim(record_dir)//'/'//trim(case_id)//trim(suffix), &
               form='unformatted',access='stream',status='replace')
          write(unit) data
          close(unit)
        end subroutine record_array

        subroutine record_complex_scalar(value, name)
          complex(8), intent(in) :: value
          character(len=*), intent(in) :: name
          integer :: unit
          if (len_trim(record_dir)==0) return
          open(newunit=unit,file=trim(record_dir)//'/'//trim(name)//'__000000__A.c8', &
               form='unformatted',access='stream',status='replace')
          write(unit) value
          close(unit)
        end subroutine record_complex_scalar

        function clone_tensor(source) result(dest)
          type(tensor_t), intent(inout) :: source
          type(tensor_t) :: dest
          integer :: ierr
          ierr = copy_tensor_t(dest,source)
        end function clone_tensor

        function clone_zero(source) result(dest)
          type(tensor_t), intent(inout) :: source
          type(tensor_t) :: dest
          integer :: ierr
          ierr = copy_tensor_t(dest,source)
          ierr = initialize_tensor_t(dest,Q_ZERO)
        end function clone_zero

        function comma_labels(labels) result(with_commas)
          character(len=*), intent(in) :: labels
          character(len=64) :: with_commas
          integer :: i
          with_commas = ''
          do i=1,len_trim(labels)
            if (i>1) with_commas = trim(with_commas)//','
            if (labels(i:i)>='A' .and. labels(i:i)<='Z') then
              with_commas = trim(with_commas)//achar(iachar(labels(i:i))+32)
            else
              with_commas = trim(with_commas)//labels(i:i)
            end if
          end do
        end function comma_labels

        function native_expression(equation) result(native)
          character(len=*), intent(in) :: equation
          character(len=256) :: native
          character(len=64) :: a_labels, b_labels, c_labels
          integer :: arrow, comma
          arrow = index(equation,'->')
          comma = index(equation,',')
          if (arrow==0 .or. comma==0) then
            native = trim(equation)
            return
          end if
          a_labels = equation(1:comma-1)
          b_labels = equation(comma+1:arrow-1)
          c_labels = equation(arrow+2:len_trim(equation))
          if (len_trim(c_labels)==0) then
            native = 'E()+=A('//trim(comma_labels(a_labels))//')*B('//trim(comma_labels(b_labels))//')'
          else
            native = 'C('//trim(comma_labels(c_labels))//')+=A('//trim(comma_labels(a_labels))// &
                     ')*B('//trim(comma_labels(b_labels))//')'
          end if
        end function native_expression

        function construct_output(A, B, equation) result(C)
          type(tensor_t), intent(inout) :: A, B
          character(len=*), intent(in) :: equation
          type(tensor_t) :: C
          character(len=64) :: a_labels, b_labels, c_labels
          integer(C_INT) :: dims_a(6), dims_b(6), dims_c(6)
          integer :: arrow, comma, i, pos, ierr
          ierr = get_extents_static(A,dims_a)
          ierr = get_extents_static(B,dims_b)
          arrow = index(equation,'->')
          comma = index(equation,',')
          a_labels = equation(1:comma-1)
          b_labels = equation(comma+1:arrow-1)
          c_labels = equation(arrow+2:len_trim(equation))
          do i=1,len_trim(c_labels)
            pos = index(a_labels,c_labels(i:i))
            if (pos>0) then
              dims_c(i) = dims_a(pos)
            else
              pos = index(b_labels,c_labels(i:i))
              dims_c(i) = dims_b(pos)
            end if
          end do
          ierr = construct_tensor_t(C,C8,dims_c(1:len_trim(c_labels)),init_val=Q_ZERO)
        end function construct_output

        function paper_binary(case_id, expression, A, B, record, component, call_site, route_name, &
                              coefficient_value) result(C)
          character(len=*), intent(in) :: case_id, expression
          type(tensor_t), intent(inout) :: A, B
          logical, intent(in), optional :: record
          character(len=*), intent(in), optional :: component, call_site, route_name, coefficient_value
          type(tensor_t) :: C, at, bt, ct
          integer :: ierr, contraction_ierr
          logical :: record_this
          character(len=16) :: this_component, this_route
          character(len=64) :: this_call_site
          character(len=8) :: this_coefficient
          record_this = .true.
          if (present(record)) record_this = record
          this_component = 'full'
          if (present(component)) this_component = component
          this_call_site = case_id
          if (present(call_site)) this_call_site = call_site
          this_route = 'ground'
          if (present(route_name)) this_route = route_name
          this_coefficient = coefficient(case_id)
          if (present(coefficient_value)) this_coefficient = coefficient_value
          C = construct_output(A,B,expression)
          ierr = convert(A,at,out_executor=1,destroy=0)
          ierr = convert(B,bt,out_executor=1,destroy=0)
          ierr = convert(C,ct,out_executor=1,destroy=0)
          contraction_ierr = t_einsum(native_expression(expression),ct,at,bt)
          ierr = destroy_tensor_t(C)
          ierr = convert(ct,C,out_executor=A%which_executor,destroy=1)
          ierr = destroy_tensor_t(at)
          ierr = destroy_tensor_t(bt)
          ierr = contraction_ierr
          if (ierr/=0) then
            write(*,*) ' QCBENCH contraction failed: ',trim(native_expression(expression)),' ierr=',ierr
            stop 'QCBENCH binary contraction failed'
          end if
          if (len_trim(record_dir)>0 .and. record_this) then
            record_calls = record_calls + 1
            record_case_ids(record_calls) = trim(case_id)
            record_call_sites(record_calls) = trim(this_call_site)
            record_routes(record_calls) = trim(this_route)
            record_coefficients(record_calls) = trim(this_coefficient)
            record_components(record_calls) = trim(this_component)
            record_equations(record_calls) = trim(expression)
            record_shapes(1,record_calls) = shape_string(A)
            record_shapes(2,record_calls) = shape_string(B)
            record_shapes(3,record_calls) = shape_string(C)
            call record_array(A,case_id,'A',record_calls)
            call record_array(B,case_id,'B',record_calls)
            call record_array(C,case_id,'C',record_calls)
            call write_manifest()
          end if
        end function paper_binary

        subroutine qcbench_pack_chol_ao(native, packed)
          type(tensor_t), intent(inout) :: native, packed
          integer(C_INT) :: dims(3), ierr, u, v, p
          integer :: nao
          type(C_PTR) :: body_p
          complex(8), pointer :: src(:,:,:), dst(:,:,:)
          ierr = get_extents_static(native,dims)
          nao = dims(1)
          ierr = construct_tensor_t(packed,C8,(/dims(3),2*nao,2*nao/),init_val=Q_ZERO)
          ierr = direct_data_pointer(native,body_p)
          call c_f_pointer(body_p,src,dims)
          ierr = direct_data_pointer(packed,body_p)
          call c_f_pointer(body_p,dst,(/dims(3),2*nao,2*nao/))
          do v=1,dims(2)
            do u=1,dims(1)
              do p=1,dims(3)
                dst(p,u,v) = src(u,v,p)
                dst(p,nao+u,nao+v) = src(u,v,p)
              end do
            end do
          end do
        end subroutine qcbench_pack_chol_ao

        subroutine qcbench_pack_chol_mo(uoo, uov, uvo, uvv, packed)
          type(tensor_t), intent(inout) :: uoo, uov, uvo, uvv, packed
          integer(C_INT) :: dims_oo(3), dims_ov(3), dims_vo(3), dims_vv(3), dims(3), ierr
          integer :: i, j, a, b, p, nocc, nvir
          type(C_PTR) :: body_p
          complex(8), pointer :: oo(:,:,:), ov(:,:,:), vo(:,:,:), vv(:,:,:), dst(:,:,:)
          ierr = get_extents_static(uoo,dims_oo)
          ierr = get_extents_static(uov,dims_ov)
          ierr = get_extents_static(uvo,dims_vo)
          ierr = get_extents_static(uvv,dims_vv)
          nocc = dims_oo(1)
          nvir = dims_vv(1)
          dims = (/dims_oo(3),nocc+nvir,nocc+nvir/)
          ierr = construct_tensor_t(packed,C8,dims,init_val=Q_ZERO)
          ierr = direct_data_pointer(uoo,body_p)
          call c_f_pointer(body_p,oo,dims_oo)
          ierr = direct_data_pointer(uov,body_p)
          call c_f_pointer(body_p,ov,dims_ov)
          ierr = direct_data_pointer(uvo,body_p)
          call c_f_pointer(body_p,vo,dims_vo)
          ierr = direct_data_pointer(uvv,body_p)
          call c_f_pointer(body_p,vv,dims_vv)
          ierr = direct_data_pointer(packed,body_p)
          call c_f_pointer(body_p,dst,dims)
          do j=1,nocc
            do i=1,nocc
              do p=1,dims(1)
                dst(p,i,j) = oo(i,j,p)
              end do
            end do
          end do
          do a=1,nvir
            do i=1,nocc
              do p=1,dims(1)
                dst(p,i,nocc+a) = ov(i,a,p)
                dst(p,nocc+a,i) = vo(a,i,p)
              end do
            end do
          end do
          do b=1,nvir
            do a=1,nvir
              do p=1,dims(1)
                dst(p,nocc+a,nocc+b) = vv(a,b,p)
              end do
            end do
          end do
        end subroutine qcbench_pack_chol_mo

        function cholesky_block(chol, first_1, count_1, first_2, count_2) result(block)
          type(tensor_t), intent(inout) :: chol
          integer, intent(in) :: first_1, count_1, first_2, count_2
          type(tensor_t) :: block
          integer(C_INT) :: dims(3), out_dims(3), ierr
          type(C_PTR) :: body_p
          complex(8), pointer :: src(:,:,:), dst(:,:,:)
          ierr = get_extents_static(chol,dims)
          out_dims = (/dims(1),count_1,count_2/)
          ierr = construct_tensor_t(block,C8,out_dims,init_val=Q_ZERO)
          ierr = direct_data_pointer(chol,body_p)
          call c_f_pointer(body_p,src,dims)
          ierr = direct_data_pointer(block,body_p)
          call c_f_pointer(body_p,dst,out_dims)
          dst = src(:,first_1:first_1+count_1-1,first_2:first_2+count_2-1)
        end function cholesky_block

        subroutine make_permutation_ij(input, output)
          type(tensor_t), intent(inout) :: input, output
          type(tensor_t) :: unit
          integer :: ierr
          integer(C_INT) :: scalar_dims(1)
          scalar_dims = 1
          ierr = construct_tensor_t(unit,C8,scalar_dims(1:0),init_val=-Q_ONE)
          output = clone_tensor(input)
          ierr = t_einsum('P(j,i,a,b)+=X(i,j,a,b)*K()',output,input,unit)
          if (ierr/=0) stop 'QCBENCH P_ij failed'
          ierr = destroy_tensor_t(unit)
        end subroutine make_permutation_ij

        subroutine make_permutation_ab(input, output)
          type(tensor_t), intent(inout) :: input, output
          type(tensor_t) :: unit
          integer :: ierr
          integer(C_INT) :: scalar_dims(1)
          scalar_dims = 1
          ierr = construct_tensor_t(unit,C8,scalar_dims(1:0),init_val=-Q_ONE)
          output = clone_tensor(input)
          ierr = t_einsum('P(i,j,b,a)+=X(i,j,a,b)*K()',output,input,unit)
          if (ierr/=0) then
            write(*,*) 'QCBENCH P_ab ierr=',ierr,' input=',trim(shape_string(input)), &
                       ' output=',trim(shape_string(output))
            stop 'QCBENCH P_ab failed'
          end if
          ierr = destroy_tensor_t(unit)
        end subroutine make_permutation_ab

        subroutine permute_ij_ab(input, output)
          type(tensor_t), intent(inout) :: input, output
          type(tensor_t) :: tmp
          integer :: ierr
          call make_permutation_ab(input,tmp)
          call make_permutation_ij(tmp,output)
          ierr = destroy_tensor_t(tmp)
        end subroutine permute_ij_ab

        subroutine add_scaled(target, source, factor)
          type(tensor_t), intent(inout) :: target, source
          complex(8), intent(in) :: factor
          integer :: ierr
          ierr = t_einsum('T(a,b,c,d)+=S(a,b,c,d)',target,source,scale=factor)
          if (ierr/=0) stop 'QCBENCH rank-4 scaled addition failed'
        end subroutine add_scaled

        subroutine add_scaled_2d(target, source, factor)
          type(tensor_t), intent(inout) :: target, source
          complex(8), intent(in) :: factor
          integer :: ierr
          ierr = t_einsum('T(a,b)+=S(a,b)',target,source,scale=factor)
          if (ierr/=0) stop 'QCBENCH rank-2 scaled addition failed'
        end subroutine add_scaled_2d

        subroutine add_binary_2d(case_id, expression, A, B, target, factor, call_site, route_name, &
                                 coefficient_value)
          character(len=*), intent(in) :: case_id, expression
          type(tensor_t), intent(inout) :: A, B, target
          complex(8), intent(in) :: factor
          character(len=*), intent(in), optional :: call_site, route_name, coefficient_value
          type(tensor_t) :: C
          integer :: ierr
          C = paper_binary(case_id,expression,A,B,call_site=call_site,route_name=route_name, &
                           coefficient_value=coefficient_value)
          call add_scaled_2d(target,C,factor)
          ierr = destroy_tensor_t(C)
        end subroutine add_binary_2d

        subroutine add_binary_4d(case_id, expression, A, B, target, factor, permutation, &
                                 call_site, route_name, coefficient_value)
          character(len=*), intent(in) :: case_id, expression
          type(tensor_t), intent(inout) :: A, B, target
          complex(8), intent(in) :: factor
          integer, intent(in) :: permutation
          character(len=*), intent(in), optional :: call_site, route_name, coefficient_value
          type(tensor_t) :: C, P
          integer :: ierr
          C = paper_binary(case_id,expression,A,B,call_site=call_site,route_name=route_name, &
                           coefficient_value=coefficient_value)
          if (permutation==0) then
            call add_scaled(target,C,factor)
          else
            if (permutation==1) call make_permutation_ab(C,P)
            if (permutation==2) call make_permutation_ij(C,P)
            if (permutation==3) call permute_ij_ab(C,P)
            call add_scaled(target,P,factor)
            ierr = destroy_tensor_t(P)
          end if
          ierr = destroy_tensor_t(C)
        end subroutine add_binary_4d

        subroutine destroy_many(tensors)
          type(tensor_t), intent(inout) :: tensors(:)
          integer :: i, ierr
          do i=1,size(tensors)
            ierr = destroy_tensor_t(tensors(i))
          end do
        end subroutine destroy_many

        subroutine qcbench_pack_coefficients(cspinor, packed)
          type(cmo), intent(in) :: cspinor
          type(tensor_t), intent(inout) :: packed
          integer(C_INT) :: dims(2), ierr
          type(C_PTR) :: body_p
          complex(8), pointer :: dst(:,:)
          dims(1) = 2*cspinor%nao
          dims(2) = cspinor%nmo
          ierr = construct_tensor_t(packed,C8,dims,init_val=Q_ZERO)
          ierr = direct_data_pointer(packed,body_p)
          call c_f_pointer(body_p,dst,dims)
          dst(1:cspinor%nao,:) = cspinor%coeff(:,:,1)
          dst(cspinor%nao+1:2*cspinor%nao,:) = cspinor%coeff(:,:,2)
        end subroutine qcbench_pack_coefficients

        function coefficient_block(coefficients, first, count, conjugate) result(block)
          type(tensor_t), intent(inout) :: coefficients
          integer, intent(in) :: first, count
          logical, intent(in), optional :: conjugate
          type(tensor_t) :: block
          integer(C_INT) :: dims(2), out_dims(2), ierr
          type(C_PTR) :: body_p
          complex(8), pointer :: src(:,:), dst(:,:)
          logical :: do_conjugate
          ierr = get_extents_static(coefficients,dims)
          out_dims = (/dims(1),count/)
          ierr = construct_tensor_t(block,C8,out_dims,init_val=Q_ZERO)
          ierr = direct_data_pointer(coefficients,body_p)
          call c_f_pointer(body_p,src,dims)
          ierr = direct_data_pointer(block,body_p)
          call c_f_pointer(body_p,dst,out_dims)
          do_conjugate = .false.
          if (present(conjugate)) do_conjugate = conjugate
          if (do_conjugate) then
            dst = conjg(src(:,first:first+count-1))
          else
            dst = src(:,first:first+count-1)
          end if
        end function coefficient_block

        subroutine qcbench_x2camf_ccsd_residual(comm_t, int_t, chol_ao, chol_mo, coefficients, &
                                               nocc, nvir, r1, r2, outdir)
          type(eti_comm_tens_t), intent(inout) :: comm_t
          type(eti_intg_tens_t), intent(inout) :: int_t
          type(tensor_t), intent(inout) :: chol_ao, chol_mo, coefficients
          integer, intent(in) :: nocc, nvir
          type(tensor_t), intent(inout) :: r1, r2
          character(len=*), intent(in), optional :: outdir
          type(tensor_t) :: t1, t2, outer, tau_tilde, tau, fmi, fme, fae, wmnij, wmbej
          type(tensor_t) :: inner_v, inner_o, inner_b, inner_m, wtmp, tmp, tmp2
          type(tensor_t) :: fvv_dressed, foo_dressed, ring_outer
          type(tensor_t) :: cvn, cv, co, ladder, projected, ao_occ, ao_mb
          type(tensor_t) :: direct, exchange, ao_g
          type(tensor_t) :: unit, lov, lvo, lvv, goovo, govov, gooov_conj
          integer :: ierr, o, nv
          integer(C_INT) :: chol_dims(3), ao_dims(3), unit_dims(1), t1_dims(2)
          o = nocc
          nv = nvir
          unit_dims = 1
          ierr = construct_tensor_t(unit,C8,unit_dims(1:0),init_val=Q_ONE)
          if (present(outdir)) then
            ierr = get_extents_static(chol_mo,chol_dims)
            ierr = get_extents_static(chol_ao,ao_dims)
            call recorder_start(outdir,nocc,nvir,chol_dims(1),ao_dims(2))
            call record_array(int_t%oooo,'input_oooo','A',0)
            call record_array(int_t%foo,'input_foo','A',0)
            call record_array(int_t%fov,'input_fov','A',0)
            call record_array(int_t%fvv,'input_fvv','A',0)
            call record_array(chol_mo,'input_chol_mo','A',0)
            call record_array(chol_ao,'input_chol_ao','A',0)
            call record_array(coefficients,'input_coefficients','A',0)
          end if

          lov = cholesky_block(chol_mo,1,nocc,nocc+1,nvir)
          lvo = cholesky_block(chol_mo,nocc+1,nvir,1,nocc)
          lvv = cholesky_block(chol_mo,nocc+1,nvir,nocc+1,nvir)

          t1_dims = (/nocc,nvir/)
          ierr = construct_tensor_t(t1,C8,t1_dims,init_val=Q_ZERO)
          ierr = t_einsum('T(i,a)+=X(a,i)*K()',t1,comm_t%t1,unit)
          ierr = construct_tensor_t(t2,C8,(/nocc,nocc,nvir,nvir/),init_val=Q_ZERO)
          ierr = t_einsum('T(i,j,a,b)+=X(a,b,i,j)*K()',t2,comm_t%t2,unit)

          outer = paper_binary('tau.outer','ia,jb->ijab',t1,t1)
          tau_tilde = clone_tensor(t2)
          tmp = clone_zero(outer)
          ierr = t_einsum('T(i,j,b,a)+=X(i,j,a,b)*K()',tmp,outer,unit)
          call add_scaled(tau_tilde,tmp,(-0.5D0)*Q_ONE)
          call add_scaled(tau_tilde,outer,0.5D0*Q_ONE)
          tau = clone_tensor(t2)
          call add_scaled(tau,outer,Q_ONE)
          call add_scaled(tau,tmp,(-1.D0)*Q_ONE)

          fmi = clone_tensor(int_t%foo)
          tmp = paper_binary('fmi.fme','me,ie->mi',int_t%fov,t1)
          call add_scaled_2d(fmi,tmp,(-0.5D0)*Q_ONE)
          tmp = paper_binary('fmi.t1_ooov','ne,mnie->mi',t1,int_t%ooov)
          call add_scaled_2d(fmi,tmp,Q_ONE)
          tmp = paper_binary('fmi.tau_oovv','inef,mnef->mi',tau_tilde,int_t%oovv)
          call add_scaled_2d(fmi,tmp,0.5D0*Q_ONE)

          inner_v = paper_binary('fae.cd_inner_v','mf,Pmf->fP',t1,lov)
          inner_o = paper_binary('fae.cd_inner_o','mf,Paf->amP',t1,lvv)
          fae = clone_tensor(int_t%fvv)
          tmp = paper_binary('fae.fme','me,ma->ae',int_t%fov,t1)
          call add_scaled_2d(fae,tmp,(-0.5D0)*Q_ONE)
          tmp = paper_binary('fae.tau_oovv','mnaf,mnef->ae',tau_tilde,int_t%oovv)
          call add_scaled_2d(fae,tmp,(-0.5D0)*Q_ONE)
          tmp = paper_binary('fae.cd_outer_v','fP,Pae->ae',inner_v,lvv)
          call add_scaled_2d(fae,tmp,Q_ONE)
          tmp = paper_binary('fae.cd_outer_o','amP,Pme->ae',inner_o,lov)
          call add_scaled_2d(fae,tmp,(-1.D0)*Q_ONE)
          fme = clone_tensor(int_t%fov)
          tmp = paper_binary('fme.t1_oovv','nf,mnef->me',t1,int_t%oovv)
          call add_scaled_2d(fme,tmp,Q_ONE)

          ierr = construct_tensor_t(goovo,C8,(/o,o,nv,o/),init_val=Q_ZERO)
          ierr = t_einsum('G(m,n,e,j)+=V(m,n,j,e)*K()',goovo,int_t%ooov,unit,scale=-Q_ONE)
          ierr = construct_tensor_t(govov,C8,(/o,nv,o,nv/),init_val=Q_ZERO)
          ierr = t_einsum('G(m,a,i,f)+=V(a,m,f,i)*K()',govov,int_t%vovo,unit)
          gooov_conj = clone_zero(int_t%ooov)
          ierr = t_einsum('G(m,n,i,e)+=V+(m,n,i,e)*K()',gooov_conj,int_t%ooov,unit)

          wtmp = paper_binary('wmnij.t1','je,mnie->mnij',t1,int_t%ooov)
          wmnij = clone_tensor(int_t%oooo)
          call add_scaled(wmnij,wtmp,Q_ONE)
          tmp2 = clone_zero(wtmp)
          ierr = t_einsum('T(m,n,j,i)+=X(m,n,i,j)*K()',tmp2,wtmp,unit)
          call add_scaled(wmnij,tmp2,(-1.D0)*Q_ONE)
          tmp = paper_binary('wmnij.tau','ijef,mnef->mnij',tau,int_t%oovv)
          call add_scaled(wmnij,tmp,0.5D0*Q_ONE)

          ierr = construct_tensor_t(wmbej,C8,(/o,nv,nv,o/),init_val=Q_ZERO)
          ierr = t_einsum('W(m,b,e,j)+=V(b,m,e,j)*K()',wmbej,int_t%vovo,unit,scale=-Q_ONE)
          tmp = paper_binary('wmbej.t1_oo','nb,mnej->mbej',t1,goovo)
          call add_scaled(wmbej,tmp,(-1.D0)*Q_ONE)
          tmp = clone_zero(t2)
          ierr = t_einsum('T(i,j,a,b)+=X(i,j,a,b)',tmp,t2,scale=0.5D0*Q_ONE)
          call add_scaled(tmp,outer,Q_ONE)
          tmp = paper_binary('wmbej.tau','jnfb,mnef->mbej',tmp,int_t%oovv)
          call add_scaled(wmbej,tmp,(-1.D0)*Q_ONE)
          inner_b = paper_binary('wmbej.cd_inner_b','jf,Pbf->jbP',t1,lvv)
          inner_m = paper_binary('wmbej.cd_inner_m','jf,Pmf->jmP',t1,lov)
          tmp = paper_binary('wmbej.cd_outer_b','jbP,Pme->mbej',inner_b,lov)
          call add_scaled(wmbej,tmp,Q_ONE)
          tmp = paper_binary('wmbej.cd_outer_m','jmP,Pbe->mbej',inner_m,lvv)
          call add_scaled(wmbej,tmp,(-1.D0)*Q_ONE)

          r1 = clone_tensor(int_t%fov)
          tmp = paper_binary('t1.fvv','ie,ae->ia',t1,fae)
          call add_scaled_2d(r1,tmp,Q_ONE)
          tmp = paper_binary('t1.foo','ma,mi->ia',t1,fmi)
          call add_scaled_2d(r1,tmp,(-1.D0)*Q_ONE)
          tmp = paper_binary('t1.fov','imae,me->ia',t2,fme)
          call add_scaled_2d(r1,tmp,Q_ONE)
          tmp = paper_binary('t1.ovov','nf,naif->ia',t1,govov)
          call add_scaled_2d(r1,tmp,(-1.D0)*Q_ONE)
          tmp = paper_binary('t1.cd_inner','imef,Pme->ifP',t2,lov)
          tmp2 = paper_binary('t1.cd_outer','ifP,Paf->ia',tmp,lvv)
          call add_scaled_2d(r1,tmp2,(-1.D0)*Q_ONE)
          tmp = paper_binary('t1.ooov','mnae,nmei->ia',t2,goovo)
          call add_scaled_2d(r1,tmp,(-0.5D0)*Q_ONE)

          fvv_dressed = clone_tensor(fae)
          tmp = paper_binary('t2.fvv_dress','mb,me->be',t1,fme)
          call add_scaled_2d(fvv_dressed,tmp,(-0.5D0)*Q_ONE)
          foo_dressed = clone_tensor(fmi)
          tmp = paper_binary('t2.foo_dress','je,me->mj',t1,fme)
          call add_scaled_2d(foo_dressed,tmp,0.5D0*Q_ONE)

          r2 = clone_zero(int_t%oovv)
          ierr = t_einsum('R(i,j,a,b)+=V+(i,j,a,b)*K()',r2,int_t%oovv,unit)
          tmp = paper_binary('t2.fvv','be,ijae->ijab',fvv_dressed,t2)
          call make_permutation_ab(tmp,tmp2)
          call add_scaled(r2,tmp2,Q_ONE)
          tmp = paper_binary('t2.foo','mj,imab->ijab',foo_dressed,t2)
          call make_permutation_ij(tmp,tmp2)
          call add_scaled(r2,tmp2,(-1.D0)*Q_ONE)
          tmp = paper_binary('t2.oooo','mnab,mnij->ijab',tau,wmnij)
          call add_scaled(r2,tmp,0.5D0*Q_ONE)

          direct = paper_binary('ao.eri_direct','Pus,Pvr->uvsr',chol_ao,chol_ao)
          exchange = paper_binary('ao.eri_exchange','Pur,Pvs->uvsr',chol_ao,chol_ao)
          ao_g = clone_tensor(direct)
          call add_scaled(ao_g,exchange,(-1.D0)*Q_ONE)
          cvn = coefficient_block(coefficients,nocc+1,nvir)
          cv = coefficient_block(coefficients,nocc+1,nvir,conjugate=.true.)
          co = coefficient_block(coefficients,1,nocc,conjugate=.true.)
          call ao_ladder(tau,cvn,chol_ao,ladder,shared_g=ao_g)
          tmp = paper_binary('t2.ao_project_a','ijuv,ua->ijva',ladder,cv)
          projected = paper_binary('t2.ao_project_b','ijva,vb->ijab',tmp,cv)
          call add_scaled(r2,projected,0.5D0*Q_ONE)
          ao_occ = paper_binary('t2.ao_occ','ijuv,um->ijvm',ladder,co)
          tmp = clone_tensor(ao_occ)
          call add_scaled(tmp,ao_occ,(-0.5D0)*Q_ONE)
          ao_mb = paper_binary('t2.ao_mb','ijvm,vb->ijmb',tmp,cv)
          tmp = paper_binary('t2.ao_last','ma,ijmb->ijba',t1,ao_mb)
          tmp2 = clone_zero(tmp)
          ierr = t_einsum('T(i,j,a,b)+=X(i,j,b,a)*K()',tmp2,tmp,unit)
          call make_permutation_ab(tmp2,tmp)
          call add_scaled(r2,tmp,(-1.D0)*Q_ONE)
          tmp = paper_binary('t2.ring_t2','imae,mbej->ijab',t2,wmbej)
          call permute_ij_ab(tmp,tmp2)
          call add_scaled(r2,tmp2,Q_ONE)
          ring_outer = paper_binary('t2.ring_t1t1_outer','ie,ma->imae',t1,t1)
          tmp = paper_binary('t2.ring_t1t1','imae,mbje->iajb',ring_outer,govov)
          ierr = destroy_tensor_t(tmp2)
          ierr = construct_tensor_t(tmp2,C8,(/o,o,nv,nv/),init_val=Q_ZERO)
          ierr = t_einsum('T(i,a,j,b)+=X(i,j,a,b)*K()',tmp2,tmp,unit)
          call permute_ij_ab(tmp2,tmp)
          call add_scaled(r2,tmp,Q_ONE)
          tmp = paper_binary('t2.tprime_inner','ie,Pae->iaP',t1,lvv)
          tmp2 = paper_binary('t2.tprime','iaP,Pbj->iajb',tmp,lvo)
          ierr = destroy_tensor_t(tmp)
          ierr = construct_tensor_t(tmp,C8,(/o,o,nv,nv/),init_val=Q_ZERO)
          ierr = t_einsum('T(i,a,j,b)+=X(i,j,a,b)*K()',tmp,tmp2,unit)
          call permute_ij_ab(tmp,tmp2)
          call add_scaled(r2,tmp2,Q_ONE)
          tmp = paper_binary('t2.last_bare','ma,ijmb->ijab',t1,gooov_conj)
          call make_permutation_ab(tmp,tmp2)
          call add_scaled(r2,tmp2,(-1.D0)*Q_ONE)

          if (present(outdir)) then
            call record_array(r1,'paper_r1','A',0)
            call record_array(r2,'paper_r2','A',0)
          end if

          ierr = destroy_tensor_t(t1)
          ierr = destroy_tensor_t(t2)
          ierr = destroy_tensor_t(outer)
          ierr = destroy_tensor_t(tau_tilde)
          ierr = destroy_tensor_t(tau)
          ierr = destroy_tensor_t(fmi)
          ierr = destroy_tensor_t(fme)
          ierr = destroy_tensor_t(fae)
          ierr = destroy_tensor_t(wmnij)
          ierr = destroy_tensor_t(wmbej)
          ierr = destroy_tensor_t(inner_v)
          ierr = destroy_tensor_t(inner_o)
          ierr = destroy_tensor_t(inner_b)
          ierr = destroy_tensor_t(inner_m)
          ierr = destroy_tensor_t(wtmp)
          ierr = destroy_tensor_t(tmp)
          ierr = destroy_tensor_t(tmp2)
          ierr = destroy_tensor_t(cvn)
          ierr = destroy_tensor_t(cv)
          ierr = destroy_tensor_t(co)
          ierr = destroy_tensor_t(ladder)
          ierr = destroy_tensor_t(projected)
          ierr = destroy_tensor_t(ao_occ)
          ierr = destroy_tensor_t(ao_mb)
          ierr = destroy_tensor_t(direct)
          ierr = destroy_tensor_t(exchange)
          ierr = destroy_tensor_t(ao_g)
          ierr = destroy_tensor_t(fvv_dressed)
          ierr = destroy_tensor_t(foo_dressed)
          ierr = destroy_tensor_t(ring_outer)
          ierr = destroy_tensor_t(lov)
          ierr = destroy_tensor_t(lvo)
          ierr = destroy_tensor_t(lvv)
          ierr = destroy_tensor_t(goovo)
          ierr = destroy_tensor_t(govov)
          ierr = destroy_tensor_t(gooov_conj)
          ierr = destroy_tensor_t(unit)
        end subroutine qcbench_x2camf_ccsd_residual

        subroutine ao_ladder(amplitude, coefficients, chol_ao, result, shared_g, eom_route)
          type(tensor_t), intent(inout) :: amplitude, coefficients, chol_ao, result
          type(tensor_t), intent(inout), optional :: shared_g
          logical, intent(in), optional :: eom_route
          type(tensor_t) :: ao_amp, direct, exchange, ao_g, tmp
          integer :: ierr
          logical :: use_eom_route
          use_eom_route = .false.
          if (present(eom_route)) use_eom_route = eom_route
          if (use_eom_route) then
            ao_amp = paper_binary('ao.virtual_transform_e','ijef,se->ijsf',amplitude,coefficients, &
                                  route_name='eom')
          else
            ao_amp = paper_binary('ao.virtual_transform_e','ijef,se->ijsf',amplitude,coefficients)
          end if
          tmp = ao_amp
          if (use_eom_route) then
            ao_amp = paper_binary('ao.virtual_transform_f','ijsf,rf->ijsr',tmp,coefficients, &
                                  route_name='eom')
          else
            ao_amp = paper_binary('ao.virtual_transform_f','ijsf,rf->ijsr',tmp,coefficients)
          end if
          ierr = destroy_tensor_t(tmp)
          if (present(shared_g)) then
            if (use_eom_route) then
              result = paper_binary('ao.ladder','uvsr,ijsr->ijuv',shared_g,ao_amp,route_name='eom')
            else
              result = paper_binary('ao.ladder','uvsr,ijsr->ijuv',shared_g,ao_amp)
            end if
          else
            direct = paper_binary('ao.eri_direct','Pus,Pvr->uvsr',chol_ao,chol_ao)
            exchange = paper_binary('ao.eri_exchange','Pur,Pvs->uvsr',chol_ao,chol_ao)
            ao_g = clone_tensor(direct)
            call add_scaled(ao_g,exchange,(-1.D0)*Q_ONE)
            result = paper_binary('ao.ladder','uvsr,ijsr->ijuv',ao_g,ao_amp)
            ierr = destroy_tensor_t(direct)
            ierr = destroy_tensor_t(exchange)
            ierr = destroy_tensor_t(ao_g)
          end if
          ierr = destroy_tensor_t(ao_amp)
        end subroutine ao_ladder

        subroutine qcbench_record_residual(native_r1, native_r2)
          type(tensor_t), intent(inout) :: native_r1, native_r2
          if (len_trim(record_dir)==0) return
          call record_array(native_r1,'native_r1','A',0)
          call record_array(native_r2,'native_r2','A',0)
        end subroutine qcbench_record_residual

        subroutine qcbench_x2camf_ccsd_energy(comm_t, int_t, energy)
          type(eti_comm_tens_t), intent(inout) :: comm_t
          type(eti_intg_tens_t), intent(inout) :: int_t
          real(8), intent(out) :: energy
          type(tensor_t) :: t1, outer, tau, tmp, result, unit
          complex(8), pointer :: scalar(:)
          type(C_PTR) :: body_p
          integer :: ierr
          integer(C_INT) :: unit_dims(1), native_dims(2), native_t2_dims(4)
          unit_dims = 1
          ierr = construct_tensor_t(unit,C8,unit_dims(1:0),init_val=Q_ONE)
          ierr = get_extents_static(comm_t%t1,native_dims)
          ierr = construct_tensor_t(t1,C8,(/native_dims(2),native_dims(1)/),init_val=Q_ZERO)
          ierr = t_einsum('T(i,a)+=X(a,i)*K()',t1,comm_t%t1,unit)
          outer = paper_binary('tau.outer','ia,jb->ijab',t1,t1,record=.false.)
          ierr = get_extents_static(comm_t%t2,native_t2_dims)
          ierr = construct_tensor_t(tau,C8,(/native_t2_dims(3),native_t2_dims(4), &
                                               native_t2_dims(1),native_t2_dims(2)/),init_val=Q_ZERO)
          ierr = t_einsum('T(i,j,a,b)+=X(a,b,i,j)*K()',tau,comm_t%t2,unit)
          call add_scaled(tau,outer,Q_ONE)
          tmp = clone_zero(outer)
          ierr = t_einsum('T(i,j,b,a)+=X(i,j,a,b)*K()',tmp,outer,unit)
          call add_scaled(tau,tmp,(-1.D0)*Q_ONE)
          result = paper_binary('energy.t1_fov','ia,ia->',t1,int_t%fov)
          tmp = paper_binary('energy.tau_oovv','ijab,ijab->',tau,int_t%oovv)
          ierr = direct_data_pointer(result,body_p)
          call c_f_pointer(body_p,scalar,(/1/))
          energy = real(scalar(1),8)
          ierr = direct_data_pointer(tmp,body_p)
          call c_f_pointer(body_p,scalar,(/1/))
          energy = energy + real(scalar(1),8) * 0.25D0
          ierr = destroy_tensor_t(t1)
          ierr = destroy_tensor_t(outer)
          ierr = destroy_tensor_t(tau)
          ierr = destroy_tensor_t(tmp)
          ierr = destroy_tensor_t(result)
          ierr = destroy_tensor_t(unit)
        end subroutine qcbench_x2camf_ccsd_energy

        subroutine qcbench_x2camf_ccsd_t_energy(comm_t, chol_mo, eps_occ, eps_vir, nocc, nvir, energy)
          type(eti_comm_tens_t), intent(inout) :: comm_t
          type(tensor_t), intent(inout) :: chol_mo
          real(8), intent(in) :: eps_occ(:), eps_vir(:)
          integer, intent(in) :: nocc, nvir
          real(8), intent(out) :: energy
          type(tensor_t) :: t1, t2, unit, loo, lov, lvo, lvv
          type(tensor_t) :: direct_inner, direct_outer, exchange_inner, exchange_outer
          type(tensor_t) :: e_jkma_direct, e_jkma_exchange, e_jkma, e_jkma_conj, occupied
          type(tensor_t) :: e_bcjk_direct, e_bcjk_exchange, e_bcjk, disconnected
          type(tensor_t) :: connected, cyclic, w, v, result
          complex(8), pointer :: connected_data(:,:,:,:,:,:), disconnected_data(:,:,:,:,:,:)
          complex(8), pointer :: direct_data(:,:,:,:,:,:), exchange_data(:,:,:,:,:,:)
          complex(8), pointer :: occupied_data(:,:,:,:,:,:), cyclic_data(:,:,:,:,:,:)
          complex(8), pointer :: w_data(:,:,:,:,:,:), v_data(:,:,:,:,:,:), scalar(:)
          type(C_PTR) :: body_p
          integer :: ierr, i, j, k, a, b, c, o, nv
          integer(C_INT) :: unit_dims(1), triple_dims(6)
          real(8) :: den
          o = nocc
          nv = nvir
          triple_dims = (/o,o,o,nv,nv,nv/)
          unit_dims = 1
          ierr = construct_tensor_t(unit,C8,unit_dims(1:0),init_val=Q_ONE)
          ierr = construct_tensor_t(t1,C8,(/o,nv/),init_val=Q_ZERO)
          ierr = t_einsum('T(i,a)+=X(a,i)*K()',t1,comm_t%t1,unit)
          ierr = construct_tensor_t(t2,C8,(/o,o,nv,nv/),init_val=Q_ZERO)
          ierr = t_einsum('T(i,j,a,b)+=X(a,b,i,j)*K()',t2,comm_t%t2,unit)
          if (len_trim(record_dir)>0) then
            call record_array(t1,'input_triples_t1','A',0)
            call record_array(t2,'input_triples_t2','A',0)
          end if

          loo = cholesky_block(chol_mo,1,o,1,o)
          lov = cholesky_block(chol_mo,1,o,o+1,nv)
          lvo = cholesky_block(chol_mo,o+1,nv,1,o)
          lvv = cholesky_block(chol_mo,o+1,nv,o+1,nv)
          direct_inner = paper_binary('T.connected.cd_direct_inner', &
                                      'jkae,Pbe->jkabP',t2,lvv)
          direct_outer = paper_binary('T.connected.cd_direct_outer', &
                                      'jkabP,Pci->ijkabc',direct_inner,lvo)
          exchange_inner = paper_binary('T.connected.cd_exchange_inner', &
                                        'jkae,Pce->jkacP',t2,lvv)
          exchange_outer = paper_binary('T.connected.cd_exchange_outer', &
                                        'jkacP,Pbi->ijkabc',exchange_inner,lvo)
          e_jkma_direct = paper_binary('T.integral.jkma_direct', &
                                       'Pjm,Pka->jkma',loo,lov)
          e_jkma_exchange = paper_binary('T.integral.jkma_exchange', &
                                         'Pja,Pkm->jkma',lov,loo)
          e_jkma = clone_tensor(e_jkma_direct)
          call add_scaled(e_jkma,e_jkma_exchange,(-1.D0)*Q_ONE)
          e_jkma_conj = clone_zero(e_jkma)
          ierr = t_einsum('G(j,k,m,a)+=V+(j,k,m,a)*K()',e_jkma_conj,e_jkma,unit)
          occupied = paper_binary('T.connected.occupied', &
                                  'imbc,jkma->ijkabc',t2,e_jkma_conj)
          e_bcjk_direct = paper_binary('T.integral.bcjk_direct', &
                                       'Pbj,Pck->bcjk',lvo,lvo)
          e_bcjk_exchange = paper_binary('T.integral.bcjk_exchange', &
                                         'Pbk,Pcj->bcjk',lvo,lvo)
          e_bcjk = clone_tensor(e_bcjk_direct)
          call add_scaled(e_bcjk,e_bcjk_exchange,(-1.D0)*Q_ONE)
          disconnected = paper_binary('T.disconnected.t1_ovov', &
                                      'ia,bcjk->ijkabc',t1,e_bcjk)

          ierr = construct_tensor_t(connected,C8,triple_dims,init_val=Q_ZERO)
          ierr = direct_data_pointer(connected,body_p)
          call c_f_pointer(body_p,connected_data,triple_dims)
          ierr = direct_data_pointer(direct_outer,body_p)
          call c_f_pointer(body_p,direct_data,triple_dims)
          ierr = direct_data_pointer(exchange_outer,body_p)
          call c_f_pointer(body_p,exchange_data,triple_dims)
          ierr = direct_data_pointer(occupied,body_p)
          call c_f_pointer(body_p,occupied_data,triple_dims)
          connected_data = direct_data - exchange_data - occupied_data
          ierr = construct_tensor_t(cyclic,C8,triple_dims,init_val=Q_ZERO)
          ierr = construct_tensor_t(w,C8,triple_dims,init_val=Q_ZERO)
          ierr = construct_tensor_t(v,C8,triple_dims,init_val=Q_ZERO)
          ierr = direct_data_pointer(cyclic,body_p)
          call c_f_pointer(body_p,cyclic_data,triple_dims)
          ierr = direct_data_pointer(disconnected,body_p)
          call c_f_pointer(body_p,disconnected_data,triple_dims)
          ierr = direct_data_pointer(w,body_p)
          call c_f_pointer(body_p,w_data,triple_dims)
          ierr = direct_data_pointer(v,body_p)
          call c_f_pointer(body_p,v_data,triple_dims)
          do c=1,nv
            do b=1,nv
              do a=1,nv
                do k=1,o
                  do j=1,o
                    do i=1,o
                      cyclic_data(i,j,k,a,b,c) = connected_data(i,j,k,a,b,c) + &
                           connected_data(k,i,j,a,b,c) + connected_data(j,k,i,a,b,c)
                    end do
                  end do
                end do
              end do
            end do
          end do
          do c=1,nv
            do b=1,nv
              do a=1,nv
                do k=1,o
                  do j=1,o
                    do i=1,o
                      den = eps_occ(i)+eps_occ(j)+eps_occ(k)-eps_vir(a)-eps_vir(b)-eps_vir(c)
                      if (a>b .and. b>c) then
                        w_data(i,j,k,a,b,c) = (cyclic_data(i,j,k,a,b,c) + &
                             cyclic_data(i,j,k,c,a,b) - cyclic_data(i,j,k,b,a,c)) / den
                      end if
                      v_data(i,j,k,a,b,c) = conjg(connected_data(i,j,k,a,b,c) + &
                           disconnected_data(i,j,k,a,b,c) + connected_data(i,j,k,c,a,b) + &
                           disconnected_data(i,j,k,c,a,b) - connected_data(i,j,k,b,a,c) - &
                           disconnected_data(i,j,k,b,a,c))
                    end do
                  end do
                end do
              end do
            end do
          end do
          result = paper_binary('T.energy.w_v','ijkabc,ijkabc->',w,v)
          ierr = direct_data_pointer(result,body_p)
          call c_f_pointer(body_p,scalar,(/1/))
          energy = 0.5D0*real(scalar(1),8)

          ierr = destroy_tensor_t(t1)
          ierr = destroy_tensor_t(t2)
          ierr = destroy_tensor_t(unit)
          ierr = destroy_tensor_t(loo)
          ierr = destroy_tensor_t(lov)
          ierr = destroy_tensor_t(lvo)
          ierr = destroy_tensor_t(lvv)
          ierr = destroy_tensor_t(direct_inner)
          ierr = destroy_tensor_t(direct_outer)
          ierr = destroy_tensor_t(exchange_inner)
          ierr = destroy_tensor_t(exchange_outer)
          ierr = destroy_tensor_t(e_jkma_direct)
          ierr = destroy_tensor_t(e_jkma_exchange)
          ierr = destroy_tensor_t(e_jkma)
          ierr = destroy_tensor_t(e_jkma_conj)
          ierr = destroy_tensor_t(occupied)
          ierr = destroy_tensor_t(e_bcjk_direct)
          ierr = destroy_tensor_t(e_bcjk_exchange)
          ierr = destroy_tensor_t(e_bcjk)
          ierr = destroy_tensor_t(disconnected)
          ierr = destroy_tensor_t(connected)
          ierr = destroy_tensor_t(cyclic)
          ierr = destroy_tensor_t(w)
          ierr = destroy_tensor_t(v)
          ierr = destroy_tensor_t(result)
        end subroutine qcbench_x2camf_ccsd_t_energy

        subroutine qcbench_record_triples(native_energy)
          real(8), intent(in) :: native_energy
          integer :: unit
          complex(8) :: value
          if (len_trim(record_dir)==0) return
          value = dcmplx(native_energy,0.D0)
          open(newunit=unit,file=trim(record_dir)//'/native_T__000000__A.c8', &
               form='unformatted',access='stream',status='replace')
          write(unit) value
          close(unit)
        end subroutine qcbench_record_triples

        subroutine qcbench_x2camf_eom_probe(int_t, interm_t, chol_ao, coefficients)
          use eti_sigma_eomccsd, only: eti_sigma_eomccsd_rhs_ee
          type(eti_intg_tens_t), intent(inout) :: int_t
          type(eti_intermediates_tens_t), intent(inout) :: interm_t
          type(tensor_t), intent(inout) :: chol_ao, coefficients
          type(eti_sa_tens_t) :: r_tensor, native_sigma
          type(tensor_t) :: raw_r2
          complex(8), pointer :: r1_data(:,:), r2_data(:,:,:,:), raw_data(:,:,:,:)
          type(C_PTR) :: body_p
          integer :: ierr, i, j, a, b, o, nv
          o = interm_t%nocc
          nv = interm_t%nvir
          ierr = construct_tensor_t(r_tensor%s1,C8,(/nv,o/),init_val=Q_ZERO)
          ierr = construct_tensor_t(r_tensor%s2,C8,(/nv,nv,o,o/),init_val=Q_ZERO)
          ierr = construct_tensor_t(raw_r2,C8,(/nv,nv,o,o/),init_val=Q_ZERO)
          ierr = construct_tensor_t(native_sigma%s1,C8,(/nv,o/),init_val=Q_ZERO)
          ierr = construct_tensor_t(native_sigma%s2,C8,(/nv,nv,o,o/),init_val=Q_ZERO)
          ierr = direct_data_pointer(r_tensor%s1,body_p)
          call c_f_pointer(body_p,r1_data,(/nv,o/))
          ierr = direct_data_pointer(r_tensor%s2,body_p)
          call c_f_pointer(body_p,r2_data,(/nv,nv,o,o/))
          ierr = direct_data_pointer(raw_r2,body_p)
          call c_f_pointer(body_p,raw_data,(/nv,nv,o,o/))
          do i=1,o
            do a=1,nv
              r1_data(a,i) = dcmplx(0.003D0*(a+2*i),0.002D0*(2*a-i))
            end do
          end do
          do j=1,o
            do i=1,o
              do b=1,nv
                do a=1,nv
                  raw_data(a,b,i,j) = dcmplx(0.001D0*(a+3*b+5*i+7*j), &
                                             0.0007D0*(2*a-b+3*i-j))
                end do
              end do
            end do
          end do
          do j=1,o
            do i=1,o
              do b=1,nv
                do a=1,nv
                  r2_data(a,b,i,j) = 0.25D0*(raw_data(a,b,i,j)-raw_data(b,a,i,j) &
                       -raw_data(a,b,j,i)+raw_data(b,a,j,i))
                end do
              end do
            end do
          end do
          call eti_sigma_eomccsd_rhs_ee(int_t,interm_t,r_tensor,native_sigma)
          call paper_eom_sigma(int_t,interm_t,chol_ao,coefficients,r_tensor,native_sigma)
          ierr = destroy_tensor_t(raw_r2)
          ierr = destroy_tensor_t(r_tensor%s1)
          ierr = destroy_tensor_t(r_tensor%s2)
          ierr = destroy_tensor_t(native_sigma%s1)
          ierr = destroy_tensor_t(native_sigma%s2)
        end subroutine qcbench_x2camf_eom_probe

        subroutine paper_eom_sigma(int_t, interm_t, chol_ao, coefficients, r_tensor, native_sigma)
          type(eti_intg_tens_t), intent(inout) :: int_t
          type(eti_intermediates_tens_t), intent(inout) :: interm_t
          type(tensor_t), intent(inout) :: chol_ao, coefficients
          type(eti_sa_tens_t), intent(inout) :: r_tensor, native_sigma
          type(tensor_t) :: unit, t1, t2, r1, r2, native_s1, native_s2, s1, s2
          type(tensor_t) :: outer, outer_swap, tau, fvv, foo, fov, w_mnij, w_mbej, w_mbij, i3, ovvv
          type(tensor_t) :: tmp, i6, i7, resp1, resp2, resp3, resp3_raw, resp4, resp4_raw
          type(tensor_t) :: scalar_a, scalar_b, cvn, cv, co, direct, exchange, ao_g, response_ao
          type(tensor_t) :: project_a, ao_virtual, ao_am, quartic_response
          complex(8), pointer :: scalar(:)
          complex(8) :: s0
          type(C_PTR) :: body_p
          integer :: ierr, o, nv
          integer(C_INT) :: unit_dims(1)
          o = interm_t%nocc
          nv = interm_t%nvir
          unit_dims = 1
          ierr = construct_tensor_t(unit,C8,unit_dims(1:0),init_val=Q_ONE)
          ierr = construct_tensor_t(t1,C8,(/o,nv/),init_val=Q_ZERO)
          ierr = t_einsum('T(i,a)+=X(a,i)*K()',t1,interm_t%gs%r1,unit)
          ierr = construct_tensor_t(t2,C8,(/o,o,nv,nv/),init_val=Q_ZERO)
          ierr = t_einsum('T(i,j,a,b)+=X(a,b,i,j)*K()',t2,interm_t%gs%r2,unit)
          ierr = construct_tensor_t(r1,C8,(/o,nv/),init_val=Q_ZERO)
          ierr = t_einsum('T(i,a)+=X(a,i)*K()',r1,r_tensor%s1,unit)
          ierr = construct_tensor_t(r2,C8,(/o,o,nv,nv/),init_val=Q_ZERO)
          ierr = t_einsum('T(i,j,a,b)+=X(a,b,i,j)*K()',r2,r_tensor%s2,unit)
          ierr = construct_tensor_t(native_s1,C8,(/o,nv/),init_val=Q_ZERO)
          ierr = t_einsum('T(i,a)+=X(a,i)*K()',native_s1,native_sigma%s1,unit)
          ierr = construct_tensor_t(native_s2,C8,(/o,o,nv,nv/),init_val=Q_ZERO)
          ierr = t_einsum('T(i,j,a,b)+=X(a,b,i,j)*K()',native_s2,native_sigma%s2,unit)
          s1 = clone_zero(r1)
          s2 = clone_zero(r2)

          outer = paper_binary('tau.outer','ia,jb->ijab',t1,t1,record=.false.)
          outer_swap = clone_zero(outer)
          ierr = t_einsum('T(i,j,b,a)+=X(i,j,a,b)*K()',outer_swap,outer,unit)
          tau = clone_tensor(t2)
          call add_scaled(tau,outer,Q_ONE)
          call add_scaled(tau,outer_swap,(-1.D0)*Q_ONE)
          fvv = clone_tensor(interm_t%fbar_vv)
          foo = clone_tensor(interm_t%fbar_oo)
          fov = clone_tensor(interm_t%fbar_ov)
          w_mnij = clone_tensor(interm_t%w_oooo)
          w_mbij = clone_tensor(interm_t%w_ovoo)
          ierr = construct_tensor_t(w_mbej,C8,(/o,nv,nv,o/),init_val=Q_ZERO)
          ierr = t_einsum('W(m,b,e,j)+=X(b,m,e,j)*K()',w_mbej,interm_t%w_vovo,unit,scale=-Q_ONE)
          ierr = construct_tensor_t(i3,C8,(/o,nv,nv,nv/),init_val=Q_ZERO)
          ierr = t_einsum('I(j,c,a,b)+=X(a,b,c,j)*K()',i3,interm_t%w_vvvo,unit,scale=-Q_ONE)
          ierr = construct_tensor_t(ovvv,C8,(/o,nv,nv,nv/),init_val=Q_ZERO)
          ierr = t_einsum('O(m,a,e,f)+=X(a,m,e,f)*K()',ovvv,int_t%vovv,unit)
          if (len_trim(record_dir)>0) then
            call record_array(r1,'input_eom_r1','A',0)
            call record_array(r2,'input_eom_r2','A',0)
            call record_array(t1,'input_eom_t1','A',0)
            call record_array(t2,'input_eom_t2','A',0)
            call record_array(tau,'input_eom_tau','A',0)
            call record_array(int_t%fov,'input_eom_bare_fov','A',0)
            call record_array(fov,'input_eom_fov','A',0)
            call record_array(fvv,'input_eom_fvv','A',0)
            call record_array(foo,'input_eom_foo','A',0)
            call record_array(w_mnij,'input_eom_w_mnij','A',0)
            call record_array(w_mbej,'input_eom_w_mbej','A',0)
            call record_array(w_mbij,'input_eom_w_mbij','A',0)
            call record_array(i3,'input_eom_i3','A',0)
            call record_array(int_t%ooov,'input_eom_ooov','A',0)
            call record_array(int_t%oovv,'input_eom_oovv','A',0)
            call record_array(ovvv,'input_eom_ovvv','A',0)
          end if

          i6 = clone_tensor(int_t%ooov)
          tmp = paper_binary('i6.t1_oovv','id,klcd->klic',t1,int_t%oovv,route_name='eom')
          call add_scaled(i6,tmp,(-1.D0)*Q_ONE)
          ierr = destroy_tensor_t(tmp)
          i7 = clone_tensor(ovvv)
          tmp = paper_binary('i7.t1_oovv','la,klcd->kacd',t1,int_t%oovv,route_name='eom')
          call add_scaled(i7,tmp,(-1.D0)*Q_ONE)
          ierr = destroy_tensor_t(tmp)
          resp1 = paper_binary('response.t1','kc,jkic->ij',r1,i6,route_name='eom')
          resp2 = paper_binary('response.t2','kc,kabc->ab',r1,i7,route_name='eom')
          resp3_raw = paper_binary('response.t3','jkab,ikab->ij',int_t%oovv,r2,route_name='eom')
          resp3 = clone_tensor(resp3_raw)
          call add_scaled_2d(resp3,resp3_raw,(-0.5D0)*Q_ONE)
          ierr = destroy_tensor_t(resp3_raw)
          resp4_raw = paper_binary('response.t4','ijbc,ijac->ab',int_t%oovv,r2,route_name='eom')
          resp4 = clone_tensor(resp4_raw)
          call add_scaled_2d(resp4,resp4_raw,(-0.5D0)*Q_ONE)
          ierr = destroy_tensor_t(resp4_raw)

          scalar_a = paper_binary('energy.t1_fov','ia,ia->',r1,int_t%fov, &
                                  call_site='sigma0.singles',route_name='eom')
          scalar_b = paper_binary('energy.tau_oovv','ijab,ijab->',r2,int_t%oovv, &
                                  call_site='sigma0.doubles',route_name='eom')
          ierr = direct_data_pointer(scalar_a,body_p)
          call c_f_pointer(body_p,scalar,(/1/))
          s0 = scalar(1)
          ierr = direct_data_pointer(scalar_b,body_p)
          call c_f_pointer(body_p,scalar,(/1/))
          s0 = s0 + 0.25D0*scalar(1)
          call record_complex_scalar(s0,'paper_eom_sigma0')
          ierr = destroy_tensor_t(scalar_a)
          ierr = destroy_tensor_t(scalar_b)

          call add_binary_2d('t1.fvv','ie,ae->ia',r1,fvv,s1,Q_ONE, &
                             call_site='sigma1.fvv',route_name='eom')
          call add_binary_2d('t1.foo','ma,mi->ia',r1,foo,s1,-Q_ONE, &
                             call_site='sigma1.foo',route_name='eom')
          call add_binary_2d('sigma1.wovvo','maei,me->ia',w_mbej,r1,s1,Q_ONE,route_name='eom')
          call add_binary_2d('t1.fov','imae,me->ia',r2,fov,s1,Q_ONE, &
                             call_site='sigma1.fov',route_name='eom')
          call add_binary_2d('sigma1.i6','jkib,jkab->ia',i6,r2,s1,(-0.5D0)*Q_ONE,route_name='eom')
          call add_binary_2d('sigma1.i7','jabc,ijbc->ia',i7,r2,s1,(-0.5D0)*Q_ONE,route_name='eom')

          call add_binary_4d('sigma2.wovoo','mbij,ma->ijab',w_mbij,r1,s2,-Q_ONE,1,route_name='eom')
          call add_binary_4d('sigma2.i3','jcab,ic->ijab',i3,r1,s2,-Q_ONE,2,route_name='eom')
          call add_binary_4d('sigma2.t1','il,jlab->ijab',resp1,t2,s2,Q_ONE,2,route_name='eom')
          call add_binary_4d('sigma2.t2','ad,ijbd->ijab',resp2,t2,s2,Q_ONE,1,route_name='eom')
          call add_binary_4d('sigma2.foo','mjab,mi->ijab',r2,foo,s2,-Q_ONE,2,route_name='eom')
          call add_binary_4d('sigma2.fvv','ijac,bc->ijab',r2,fvv,s2,Q_ONE,1,route_name='eom')
          call add_binary_4d('sigma2.wovvo','mbej,imae->ijab',w_mbej,r2,s2,Q_ONE,3,route_name='eom')
          call add_binary_4d('t2.oooo','mnab,mnij->ijab',r2,w_mnij,s2,0.5D0*Q_ONE,0, &
                             call_site='sigma2.woooo',route_name='eom')
          call add_binary_4d('sigma2.t3','il,jlab->ijab',resp3,t2,s2,Q_ONE,2,route_name='eom')
          call add_binary_4d('sigma2.t4','ad,ijbd->ijab',resp4,t2,s2,Q_ONE,1,route_name='eom')

          cvn = coefficient_block(coefficients,o+1,nv)
          cv = coefficient_block(coefficients,o+1,nv,conjugate=.true.)
          co = coefficient_block(coefficients,1,o,conjugate=.true.)
          direct = paper_binary('ao.eri_direct','Pus,Pvr->uvsr',chol_ao,chol_ao,route_name='eom')
          exchange = paper_binary('ao.eri_exchange','Pur,Pvs->uvsr',chol_ao,chol_ao, &
                                  coefficient_value='-1',route_name='eom')
          ao_g = clone_tensor(direct)
          call add_scaled(ao_g,exchange,(-1.D0)*Q_ONE)
          call ao_ladder(r2,cvn,chol_ao,response_ao,shared_g=ao_g,eom_route=.true.)
          project_a = paper_binary('t2.ao_project_a','ijuv,ua->ijva',response_ao,cv, &
                                   call_site='ao.project_a',route_name='eom')
          call add_binary_4d('t2.ao_project_b','ijva,vb->ijab',project_a,cv,s2,0.5D0*Q_ONE,0, &
                             call_site='ao.project_b',route_name='eom')
          ao_virtual = paper_binary('t2.ao_project_a','ijuv,ua->ijva',response_ao,cv, &
                                    call_site='ao.virtual_for_correction',route_name='eom', &
                                    coefficient_value='+1')
          ao_am = paper_binary('ao.occupied_for_correction','ijva,vm->ijam',ao_virtual,co,route_name='eom')
          call add_binary_4d('ao.t1_correction','mb,ijam->ijab',t1,ao_am,s2,(-0.5D0)*Q_ONE,1, &
                             route_name='eom')
          quartic_response = paper_binary('quartic.response_oovv','ijef,mnef->ijmn',r2,int_t%oovv, &
                                          route_name='eom')
          call add_binary_4d('quartic.tau','mnab,ijmn->ijab',tau,quartic_response,s2,0.25D0*Q_ONE,0, &
                             route_name='eom')

          if (len_trim(record_dir)>0) then
            call record_array(native_s1,'native_eom_sigma1','A',0)
            call record_array(native_s2,'native_eom_sigma2','A',0)
            call record_array(s1,'paper_eom_sigma1','A',0)
            call record_array(s2,'paper_eom_sigma2','A',0)
          end if
          ierr = destroy_tensor_t(unit)
          ierr = destroy_tensor_t(t1)
          ierr = destroy_tensor_t(t2)
          ierr = destroy_tensor_t(r1)
          ierr = destroy_tensor_t(r2)
          ierr = destroy_tensor_t(native_s1)
          ierr = destroy_tensor_t(native_s2)
          ierr = destroy_tensor_t(s1)
          ierr = destroy_tensor_t(s2)
          ierr = destroy_tensor_t(outer)
          ierr = destroy_tensor_t(outer_swap)
          ierr = destroy_tensor_t(tau)
          ierr = destroy_tensor_t(fvv)
          ierr = destroy_tensor_t(foo)
          ierr = destroy_tensor_t(fov)
          ierr = destroy_tensor_t(w_mnij)
          ierr = destroy_tensor_t(w_mbej)
          ierr = destroy_tensor_t(w_mbij)
          ierr = destroy_tensor_t(i3)
          ierr = destroy_tensor_t(ovvv)
          ierr = destroy_tensor_t(i6)
          ierr = destroy_tensor_t(i7)
          ierr = destroy_tensor_t(resp1)
          ierr = destroy_tensor_t(resp2)
          ierr = destroy_tensor_t(resp3)
          ierr = destroy_tensor_t(resp4)
          ierr = destroy_tensor_t(cvn)
          ierr = destroy_tensor_t(cv)
          ierr = destroy_tensor_t(co)
          ierr = destroy_tensor_t(direct)
          ierr = destroy_tensor_t(exchange)
          ierr = destroy_tensor_t(ao_g)
          ierr = destroy_tensor_t(response_ao)
          ierr = destroy_tensor_t(project_a)
          ierr = destroy_tensor_t(ao_virtual)
          ierr = destroy_tensor_t(ao_am)
          ierr = destroy_tensor_t(quartic_response)
        end subroutine paper_eom_sigma

end module qcbench_x2camf
