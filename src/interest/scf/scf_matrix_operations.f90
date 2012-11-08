! ------------------------------------------------------------------------------------
!
! Program:      Dirac 
!
! Module:       Kramers unrestricted SCF module 
!
! Description:  
!
! Contains:     
!
! Licensing:
!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end
!
! Author:       Michal Repisky (michal.repisky@uit.no)
!
! Logs:         
!
! Comments:     
!                            
! ------------------------------------------------------------------------------------

MODULE SCF_MATRIX_OPERATIONS

  use scf_data_types

  implicit none


  !> public

  public index
  public matrix_allocate
  public matrix_deallocate
  public matrix_dotproduct
  public matrix_multiplication 

  public test_association 


  !> private

  private

  interface matrix_allocate 
    module procedure matrix_allocate_1comp
    module procedure matrix_allocate_4comp
  end interface

  interface matrix_deallocate 
    module procedure matrix_deallocate_1comp
    module procedure matrix_deallocate_4comp
  end interface

  interface matrix_dotproduct
    module procedure matrix_dotproduct_1comp
    module procedure matrix_dotproduct_4comp
  end interface

  interface matrix_multiplication
    module procedure matrix_multiplication_1comp
    module procedure matrix_multiplication_4comp
  end interface

CONTAINS
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE matrix_allocate_1comp( A )

    !> input/output
       type(real_matrix), intent(inout) :: A
    !> end
 
    !> local


    !> if associated then nullify the matrix otherwise perform the allocation
    !> note:  real_matrix has pointer type elements => association check is relevant
       if( associated(A%elms) )then
         A%nrow = scf%ndim 
         A%ncol = scf%ndim
         A%elms = 0.0d0
       else
         A%nrow = scf%ndim
         A%ncol = scf%ndim
         allocate( A%elms( A%nrow * A%ncol, scf%nshl ) ); A%elms = 0.0d0
       endif
    !> end
    
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE matrix_allocate_4comp( A )

    !> input
    type(complex_matrix), intent(inout) :: A
 
    !> local



    !> if associated then nullify the matrix otherwise perform the allocation
    !> note: complex_matrix has allocatable type elements => allocated check is relevant
       if( allocated(A%elms) )then
         A%nrow = scf%ndim 
         A%ncol = scf%ndim
         A%elms = ( 0.0d0, 0.0d0 )
       else
         A%nrow = scf%ndim 
         A%ncol = scf%ndim
         allocate( A%elms( A%nrow * A%ncol, 4, 4 ) ); A%elms = ( 0.0d0, 0.0d0 )
       endif
    !> end
    
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE matrix_deallocate_1comp( A )

    !> input
    type(real_matrix), intent(inout) :: A
 
    !> local


    A%nrow = 0
    A%ncol = 0

    deallocate( A%elms )
    
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE matrix_deallocate_4comp( A )

    !> input
    type(complex_matrix), intent(inout) :: A
 
    !> local


    A%nrow = 0
    A%ncol = 0

    call deassociate( A = A, full = .true. )
    deallocate( A%elms )
    
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE matrix_multiplication_1comp( alpha, opa, A, opb, B, beta, C)

    !> input
    type(real_matrix),          intent(in) :: A
    type(real_matrix),          intent(in) :: B
    character(len=1), optional, intent(in) :: opa
    character(len=1), optional, intent(in) :: opb
    real(8),          optional, intent(in) :: beta
    real(8),          optional, intent(in) :: alpha

    !> output 
    type(real_matrix), intent(inout) :: C

    !> local 
    integer          :: i
    character(len=1) :: opa_final
    character(len=1) :: opb_final
    real(8)          :: beta_final 
    real(8)          :: alpha_final 


    opa_final   = 'n'
    opb_final   = 'n'
    beta_final  = 0.0d0
    alpha_final = 1.0d0

    if( present(opa)   ) opa_final   = opa 
    if( present(opb)   ) opb_final   = opb
    if( present(beta)  ) beta_final  = beta
    if( present(alpha) ) alpha_final = alpha

    !> 'n','n' case: A(i,k)*B(k,j)
    if( opa_final=='n' .and. opb_final=='n' )then
      do i=1,ubound(A%elms,2) 
        call dgemm( opa_final,                      & 
                    opb_final,                      &
                    A%nrow,B%ncol,A%ncol,           &
                    alpha_final,A%elms(1,i),A%nrow, &
                                B%elms(1,i),B%nrow, &
                    beta_final, C%elms(1,i),C%nrow  )
      enddo
      return

    !> 't/c','n' case: A(k,i)*B(k,j)
    elseif( (opa_final=='t'.or.opa_final=='c') .and. opb_final=='n' )then
      do i=1,ubound(A%elms,2) 
        call dgemm( opa_final,                      & 
                    opb_final,                      &
                    A%ncol,B%ncol,A%nrow,           &
                    alpha_final,A%elms(1,i),A%nrow, &
                                B%elms(1,i),B%nrow, &
                    beta_final, C%elms(1,i),C%nrow  )
      enddo
      return

    !> 'n','t/c' case: A(i,k)*B(j,k)
    elseif( opa_final=='n' .and. (opb_final=='t'.or.opb_final=='c') )then
      do i=1,ubound(A%elms,2) 
        call dgemm( opa_final,                      & 
                    opb_final,                      &
                    A%nrow,B%nrow,A%ncol,           &
                    alpha_final,A%elms(1,i),A%nrow, &
                                B%elms(1,i),B%nrow, &
                    beta_final, C%elms(1,i),C%nrow  )
      enddo
      return

    !> 't/c','t/c' case: A(k,i)*B(j,k)
    elseif( (opa_final=='t'.or.opa_final=='c') .and. (opb_final=='t'.or.opb_final=='c') )then
      do i=1,ubound(A%elms,2) 
        call dgemm( opa_final,                      & 
                    opb_final,                      &
                    A%ncol,B%nrow,A%nrow,           &
                    alpha_final,A%elms(1,i),A%nrow, &
                                B%elms(1,i),B%nrow, &
                    beta_final, C%elms(1,i),C%nrow  )
      enddo
      return
    endif

    write(6,'(2x,a)') 'A matrix operation on A/B is badly defined => stop'
    stop
     
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE matrix_multiplication_4comp( alpha, opa, A, opb, B, beta, C)

    !> input
    type(complex_matrix),       intent(in) :: A
    type(complex_matrix),       intent(in) :: B
    character(len=1), optional, intent(in) :: opa
    character(len=1), optional, intent(in) :: opb
    real(8),          optional, intent(in) :: beta
    real(8),          optional, intent(in) :: alpha

    !> output 
    type(complex_matrix), target :: C

    !> local 
    integer          :: i
    integer          :: j
    integer          :: k
    character(len=1) :: opa_final
    character(len=1) :: opb_final
    complex(8)       :: beta_final 
    complex(8)       :: alpha_final 


    opa_final   = 'n'
    opb_final   = 'n'
    beta_final  = ( 0.0d0, 0.0d0 )
    alpha_final = ( 1.0d0, 0.0d0 )

    if( present(opa)   ) opa_final   = opa 
    if( present(opb)   ) opb_final   = opb
    if( present(beta)  ) beta_final  = dcmplx( beta,  0.0d0 )
    if( present(alpha) ) alpha_final = dcmplx( alpha, 0.0d0 )


    !> protection agains uncorrect zeroing in the block-matrix multiplication 
    if( beta_final == (0.0d0,0.0d0) )then
      call deassociate( A=C, full=.true. )
      C%elms = ( 0.0d0, 0.0d0 )  
    else
      C%elms = beta_final*C%elms 
    endif
    beta_final = ( 1.0d0, 0.0d0 )


    !> 'n','n' case: A(i,k)*B(k,j)
    if( opa_final=='n' .and. opb_final=='n' )then
      do j=1,4
        do i=1,4
          do k=1,4
            if( .not.associated( A%info(i,k)%p ) .or. .not.associated( B%info(k,j)%p ) )cycle
            call zgemm( opa_final,                        & 
                        opb_final,                        &
                        A%nrow,B%ncol,A%ncol,             &
                        alpha_final,A%elms(1,i,k),A%nrow, &
                                    B%elms(1,k,j),B%nrow, &
                        beta_final, C%elms(1,i,j),C%nrow  )
            C%info(i,j)%p => C%elms(1,i,j)
          enddo
        enddo
      enddo
      return

    !> 't/c','n' case: A(k,i)*B(k,j)
    elseif( (opa_final=='t'.or.opa_final=='c') .and. opb_final=='n' )then
      do j=1,4
        do i=1,4
          do k=1,4
            if( .not.associated( A%info(k,i)%p ) .or. .not.associated( B%info(k,j)%p ) )cycle
            call zgemm( opa_final,                        & 
                        opb_final,                        &
                        A%ncol,B%ncol,A%nrow,             &
                        alpha_final,A%elms(1,k,i),A%nrow, &
                                    B%elms(1,k,j),B%nrow, &
                        beta_final, C%elms(1,i,j),C%nrow  )
            C%info(i,j)%p => C%elms(1,i,j)
          enddo
        enddo
      enddo
      return

    !> 'n','t/c' case: A(i,k)*B(j,k)
    elseif( opa_final=='n' .and. (opb_final=='t'.or.opb_final=='c') )then
      do j=1,4
        do i=1,4
          do k=1,4
            if( .not.associated( A%info(i,k)%p ) .or. .not.associated( B%info(j,k)%p ) )cycle
            call zgemm( opa_final,                        & 
                        opb_final,                        &
                        A%nrow,B%nrow,A%ncol,             &
                        alpha_final,A%elms(1,i,k),A%nrow, &
                                    B%elms(1,j,k),B%nrow, &
                        beta_final, C%elms(1,i,j),C%nrow  )
            C%info(i,j)%p => C%elms(1,i,j)
          enddo
        enddo
      enddo
      return

    !> 't/c','t/c' case: A(k,i)*B(j,k)
    elseif( (opa_final=='t'.or.opa_final=='c') .and. (opb_final=='t'.or.opb_final=='c') )then
      do j=1,4
        do i=1,4
          do k=1,4
            if( .not.associated( A%info(k,i)%p ) .or. .not.associated( B%info(j,k)%p ) )cycle
            call zgemm( opa_final,                        & 
                        opb_final,                        &
                        A%ncol,B%nrow,A%nrow,             &
                        alpha_final,A%elms(1,k,i),A%nrow, &
                                    B%elms(1,j,k),B%nrow, &
                        beta_final, C%elms(1,i,j),C%nrow  )
            C%info(i,j)%p => C%elms(1,i,j)
          enddo
        enddo
      enddo
      return
    endif

    write(6,'(2x,a)') 'A matrix operation on A/B is badly defined => stop'
    stop
     
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  FUNCTION matrix_dotproduct_1comp( opa, A, opb, B ) result(trace)

    !> input
    type(real_matrix),          intent(in) :: A
    type(real_matrix),          intent(in) :: B
    character(len=1), optional, intent(in) :: opa
    character(len=1), optional, intent(in) :: opb

    !> output 
    real(8) :: trace 

    !> local 
    integer          :: i
    integer          :: m
    integer          :: n
    character(len=1) :: opa_final
    character(len=1) :: opb_final


    opa_final = 'n'
    opb_final = 'n'

    if( present(opa) ) opa_final = opa 
    if( present(opb) ) opb_final = opb


    !> 'n','n' case: A(i,j)*B(j,i)
    if( opa_final=='n' .and. opb_final=='n' )then
      trace = 0.0d0
      do i=1,ubound(A%elms,2)
        do n=1,A%ncol 
          do m=1,A%nrow 
            trace = trace + A%elms(index(n,m,A%ncol),i)*B%elms(index(m,n,B%nrow),i)
          enddo
        enddo
      enddo
      return
    endif

    write(6,'(2x,a)') 'A trace operation on A/B is badly defined => stop'
    stop
     
  END FUNCTION
! ------------------------------------------------------------------------------------
!> 
  FUNCTION matrix_dotproduct_4comp( opa, A, opb, B ) result(trace)

    !> input
    type(complex_matrix),       intent(in) :: A
    type(complex_matrix),       intent(in) :: B
    character(len=1), optional, intent(in) :: opa
    character(len=1), optional, intent(in) :: opb

    !> output 
    complex(8) :: trace 

    !> local 
    integer          :: i
    integer          :: j
    integer          :: m
    integer          :: n
    character(len=1) :: opa_final
    character(len=1) :: opb_final


    opa_final = 'n'
    opb_final = 'n'

    if( present(opa) ) opa_final = opa 
    if( present(opb) ) opb_final = opb


    !> 'n','n' case: A(i,j)*B(j,i)
    if( opa_final=='n' .and. opb_final=='n' )then
      trace = ( 0.0d0, 0.0d0 )
      do j=1,4
        do i=1,4
          if( .not.associated( A%info(i,j)%p ) .or. .not.associated( B%info(j,i)%p ) )cycle
          do n=1,A%nrow 
            do m=1,A%nrow 
              trace = trace + A%elms(index(n,m,A%nrow),i,j)*B%elms(index(m,n,B%nrow),i,j)
            enddo
          enddo
        enddo
      enddo
      return
    endif

    write(6,'(2x,a)') 'A trace operation on A/B is badly defined => stop'
    stop
     
  END FUNCTION
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE associate( A, aall, aals, aasl, aass, &
                           abll, abls, absl, abss, &
                           ball, bals, basl, bass, & 
                           bbll, bbls, bbsl, bbss, &
                           full                    )
    !> input
    type(complex_matrix), target, intent(inout) :: A
    logical, optional,            intent(in)    :: full 
    logical, optional,            intent(in)    :: aall, aals, aasl, aass 
    logical, optional,            intent(in)    :: abll, abls, absl, abss
    logical, optional,            intent(in)    :: ball, bals, basl, bass
    logical, optional,            intent(in)    :: bbll, bbls, bbsl, bbss

    !> local
    integer :: i
    integer :: j


    if( present(full) )then
      do j=1,4
        do i=1,4
          A%info(i,j)%p => A%elms(1,i,j)
        enddo
      enddo
      return
    endif

    if( present(aall) ) A%info(1,1)%p => A%elms(1,1,1)
    if( present(aals) ) A%info(1,2)%p => A%elms(1,1,2)
    if( present(aasl) ) A%info(2,1)%p => A%elms(1,2,1)
    if( present(aass) ) A%info(2,2)%p => A%elms(1,2,2)

    if( present(abll) ) A%info(1,3)%p => A%elms(1,1,3)
    if( present(abls) ) A%info(1,4)%p => A%elms(1,1,4)
    if( present(absl) ) A%info(2,3)%p => A%elms(1,2,3)
    if( present(abss) ) A%info(2,4)%p => A%elms(1,2,4)

    if( present(ball) ) A%info(3,1)%p => A%elms(1,3,1)
    if( present(bals) ) A%info(3,2)%p => A%elms(1,3,2)
    if( present(basl) ) A%info(4,1)%p => A%elms(1,4,1)
    if( present(bass) ) A%info(4,2)%p => A%elms(1,4,2)

    if( present(bbll) ) A%info(3,3)%p => A%elms(1,3,3)
    if( present(bbls) ) A%info(3,4)%p => A%elms(1,3,4)
    if( present(bbsl) ) A%info(4,3)%p => A%elms(1,4,3)
    if( present(bbss) ) A%info(4,4)%p => A%elms(1,4,4)

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE deassociate( A, aall, aals, aasl, aass, &
                             abll, abls, absl, abss, &
                             ball, bals, basl, bass, & 
                             bbll, bbls, bbsl, bbss, &
                             full                    )
    !> input
    type(complex_matrix), target, intent(inout) :: A
    logical, optional,            intent(in)    :: full 
    logical, optional,            intent(in)    :: aall, aals, aasl, aass 
    logical, optional,            intent(in)    :: abll, abls, absl, abss
    logical, optional,            intent(in)    :: ball, bals, basl, bass
    logical, optional,            intent(in)    :: bbll, bbls, bbsl, bbss

    !> local
    integer :: i
    integer :: j


    if( present(full) )then
      do j=1,4
        do i=1,4
          A%info(i,j)%p => null() 
        enddo
      enddo
      return
    endif

    if( present(aall) ) A%info(1,1)%p => null() 
    if( present(aals) ) A%info(1,2)%p => null() 
    if( present(aasl) ) A%info(2,1)%p => null() 
    if( present(aass) ) A%info(2,2)%p => null() 

    if( present(abll) ) A%info(1,3)%p => null() 
    if( present(abls) ) A%info(1,4)%p => null() 
    if( present(absl) ) A%info(2,3)%p => null() 
    if( present(abss) ) A%info(2,4)%p => null() 

    if( present(ball) ) A%info(3,1)%p => null() 
    if( present(bals) ) A%info(3,2)%p => null() 
    if( present(basl) ) A%info(4,1)%p => null() 
    if( present(bass) ) A%info(4,2)%p => null() 

    if( present(bbll) ) A%info(3,3)%p => null() 
    if( present(bbls) ) A%info(3,4)%p => null() 
    if( present(bbsl) ) A%info(4,3)%p => null() 
    if( present(bbss) ) A%info(4,4)%p => null() 

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  FUNCTION index( i, j, n )

    !> input
    integer :: i
    integer :: j
    integer :: n

    !> output 
    integer :: index

    index = (j-1)*n + i
     
  END FUNCTION
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE test_association( A, name )

    !> input
    type(complex_matrix), intent(in) :: A
    character(len=*),     intent(in) :: name

    !> local
    integer :: i
    integer :: j


    write(6,*)
    do j=1,4
      do i=1,4
        if(       associated( A%info(i,j)%p ) ) write(6,'(2x,2(2a,2(i1,a),5x),2d20.10)')    &
                                                      name ,'(',i,',',j,')     associated', &
                                                      name ,'(',i,',',j,') = ',A%info(i,j)%p   
        
        if( .not. associated( A%info(i,j)%p ) ) write(6,'(2x,1(2a,2(i1,a))       )')        &
                                                      name ,'(',i,',',j,') not associated'
      enddo
    enddo
    write(6,*)
     
  END SUBROUTINE
! ------------------------------------------------------------------------------------
END MODULE
