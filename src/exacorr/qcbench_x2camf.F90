module qcbench_x2camf
! Canonical X2CAMF-CCSD paper binaries and operand recording for qcbench.

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
        public qcbench_record_residual

        character(len=512) :: record_dir = ''
        integer :: record_calls = 0
        integer :: record_nocc = 0, record_nvir = 0, record_naux = 0, record_nao = 0
        character(len=64) :: record_case_ids(10000)
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
                't2.ao_project_a','t2.foo_dress')
            value = '+1/2'
          case ('fmi.fme','fae.fme','fae.tau_oovv','t1.ooov','t2.fvv_dress')
            value = '-1/2'
          case ('energy.tau_oovv')
            value = '+1/4'
          case ('wmbej.t1_oo','wmbej.tau','wmbej.cd_outer_m','fae.cd_outer_o', &
                't1.foo','t1.cd_outer','t1.ovov','t2.foo','t2.last_bare','t2.ao_last')
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
          write(unit,'(a)') '  "route": "DIRAC X2CAMF-CCSD paper route (Z2023 Eqs. 5-12 and Table 2)",'
          write(unit,'(a)') '  "corpus": "cd_x2camf_ccsd_corpus",'
          write(unit,'(a,i0,a)') '  "n_contractions": ',record_calls,','
          write(unit,'(a,i0,a,i0,a,i0,a,i0,a)') '  "dimensions": {"o": ',record_nocc,', "v": ',record_nvir, &
               ', "ao": ',record_nao,', "naux": ',record_naux,'},'
          write(unit,'(a)') '  "calls": ['
          do i=1,record_calls
            write(unit,'(a)') '    {'
            write(unit,'(a)') '      "case": "'//trim(record_case_ids(i))//'",'
            write(unit,'(a)') '      "corpus_case": "'//trim(record_case_ids(i))//'",'
            write(unit,'(a)') '      "equation": "'//trim(record_equations(i))//'",'
            write(unit,'(a)') '      "coefficient": "'//trim(coefficient(record_case_ids(i)))//'",'
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

        function paper_binary(case_id, expression, A, B, record, component) result(C)
          character(len=*), intent(in) :: case_id, expression
          type(tensor_t), intent(inout) :: A, B
          logical, intent(in), optional :: record
          character(len=*), intent(in), optional :: component
          type(tensor_t) :: C, at, bt, ct
          integer :: ierr, contraction_ierr
          logical :: record_this
          character(len=16) :: this_component
          record_this = .true.
          if (present(record)) record_this = record
          this_component = 'full'
          if (present(component)) this_component = component
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

        subroutine ao_ladder(amplitude, coefficients, chol_ao, result, shared_g)
          type(tensor_t), intent(inout) :: amplitude, coefficients, chol_ao, result
          type(tensor_t), intent(inout), optional :: shared_g
          type(tensor_t) :: ao_amp, direct, exchange, ao_g, tmp
          integer :: ierr
          ao_amp = paper_binary('ao.virtual_transform_e','ijef,se->ijsf',amplitude,coefficients)
          tmp = ao_amp
          ao_amp = paper_binary('ao.virtual_transform_f','ijsf,rf->ijsr',tmp,coefficients)
          ierr = destroy_tensor_t(tmp)
          if (present(shared_g)) then
            result = paper_binary('ao.ladder','uvsr,ijsr->ijuv',shared_g,ao_amp)
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

end module qcbench_x2camf
