module exatensor

use, intrinsic:: ISO_C_BINDING

private

!BASIC TYPE KINDS of TALSH/EXATENSOR:
       integer(C_INT), parameter:: INTB=1   !byte integer size
       integer(C_INT), parameter:: INTS=2   !short integer size
       integer(C_INT), parameter, public:: INTD=4   !default integer size
       integer(C_INT), parameter, public:: INTL=8   !long integer size
       integer(C_INT), parameter:: REALH=2  !half real size
       integer(C_INT), parameter:: REALS=4  !short real size
       integer(C_INT), parameter:: REALD=8  !default real size
       integer(C_INT), parameter:: REALL=16 !long real size
!BASIC NUMERIC DATA KINDS (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter:: NO_TYPE=0 !no type/kind
        integer(C_INT), parameter:: R2=2      !half-precision float tensor data kind
        integer(C_INT), parameter:: R4=4      !single-precision float tensor data kind
        integer(C_INT), parameter:: R8=8      !double-precision float tensor data kind
!       integer(C_INT), parameter:: R16=10    !quadruple-precision float tensor data kind
        integer(C_INT), parameter:: C2=12     !half-precision complex float tensor data kind
        integer(C_INT), parameter:: C4=14     !single-precision complex float tensor data kind
        integer(C_INT), parameter:: C8=18     !double-precision complex float tensor data kind
!       integer(C_INT), parameter:: C16=20    !quadruple-precision complex float tensor data kind
!Tensor data kinds of EXATENSOR which all should be made public::
       integer(INTD), parameter, public:: EXA_DATA_KIND_NN=NO_TYPE !none
       integer(INTD), parameter, public:: EXA_DATA_KIND_R4=R4      !single precision real
       integer(INTD), parameter, public:: EXA_DATA_KIND_R8=R8      !double precision real
       integer(INTD), parameter, public:: EXA_DATA_KIND_C4=C4      !single precision complex
       integer(INTD), parameter, public:: EXA_DATA_KIND_C8=C8      !double precision complex


!    Error code thrown when calling one of the placeholder routines
       integer(INTD), parameter :: EXA_NOT_AVAILABLE = 1
!
!    Placeholders for ExaTensor data types
!
       type, public:: tens_rcrsv_t
       end type
 
!
!    Placeholders for ExaTensor interfaces
!
       public exatns_tensor_contract, exatns_tensor_init, exatns_tensor_create, exatns_tensor_destroy

  !Create a tensor:
       interface exatns_tensor_create
        module procedure exatns_tensor_create_scalar
        module procedure exatns_tensor_create_tensor
       end interface exatns_tensor_create
  !Initialize a tensor:
       interface exatns_tensor_init
        module procedure exatns_tensor_init_scalar
        module procedure exatns_tensor_init_method_static
       !module procedure exatns_tensor_init_method_dynamic
       end interface exatns_tensor_init
  !Copy a tensor:
       interface exatns_tensor_copy
       !module procedure exatns_tensor_copy_us
       !module procedure exatns_tensor_copy_eu
       end interface exatns_tensor_copy
  !Contract tensors:
       interface exatns_tensor_contract
        module procedure exatns_tensor_contract_eu
       end interface exatns_tensor_contract
  !Destroy a tensor:
       interface exatns_tensor_destroy
        module procedure exatns_tensor_destroy
       end interface exatns_tensor_destroy

contains
!
!    Placeholders for ExaTensor procedures
!
      function exatns_tensor_contract_eu(pattern,tensor0,tensor1,tensor2,prefactor,restrictions,sync) result(ierr)
      !Operation: Tensor contraction: tensor0 += tensor1 * tensor2 * prefactor
        implicit none
        integer(INTD):: ierr
        type(tens_rcrsv_t), intent(inout):: tensor0       !inout: tensor
        type(tens_rcrsv_t), intent(inout):: tensor1       !in: tensor
        type(tens_rcrsv_t), intent(inout):: tensor2       !in: tensor
        character(*), intent(in):: pattern                !in: symbolic tensor contraction pattern
        complex(8), intent(in), optional:: prefactor      !in: scalar prefactor
        character(*), intent(in), optional:: restrictions !in: symbolic specification of operational index restrictions
        logical, intent(in), optional:: sync              !in: if FALSE, the operation will run asynchronously, requiring a separate exatns_sync() call later (defaults to TRUE)
      
        ierr=EXA_NOT_AVAILABLE
        return
      end function exatns_tensor_contract_eu

      function exatns_tensor_init_scalar(tensor,scalar,sync) result(ierr)
        !Initializes all tensor elements to a given scalar value. If <scalar> is
        !absent, tensor elements will be initialized to random values.
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !inout: tensor
        complex(8), intent(in), optional:: scalar  !in: scalar
        logical, intent(in), optional:: sync       !in: if FALSE, the operation will run asynchronously, requiring a separate exatns_sync() call later (defaults to TRUE)

        ierr=EXA_NOT_AVAILABLE
        return
      end function exatns_tensor_init_scalar

      function exatns_tensor_init_method_static(tensor,method,sync) result(ierr)
        !Initializes the tensor via a static user-defined (registered) method.
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !inout: tensor
        character(*), intent(in):: method          !in: registered (static) method name
        logical, intent(in), optional:: sync       !in: if FALSE, the operation will run asynchronously, requiring a separate exatns_sync() call later (defaults to TRUE)

        ierr=EXA_NOT_AVAILABLE
        return
      end function exatns_tensor_init_method_static

      function exatns_tensor_create_scalar(tensor,tens_name,data_kind,sync) result(ierr)
        !Creates a scalar (order-0 tensor).
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !out: tensor (must be empty on entrance)
        character(*), intent(in):: tens_name       !in: symbolic tensor name
        integer(INTD), intent(in):: data_kind      !in: data kind of tensor elements: {EXA_DATA_KIND_XX: XX = R4,R8,C4,C8}
        logical, intent(in), optional:: sync       !in: if FALSE, the operation will run asynchronously, requiring a separate exatns_sync() call later (defaults to TRUE)

        ierr=EXA_NOT_AVAILABLE
        return
      end function exatns_tensor_create_scalar

      function exatns_tensor_create_tensor(tensor,tens_name,hspaces,subspaces,data_kind,dim_extent,dim_group,group_spec,sync)&
                                           &result(ierr)
        !Creates a tensor.
        implicit none
        integer(INTD):: ierr                                 !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor           !out: tensor (must be empty on entrance)
        character(*), intent(in):: tens_name                 !in: symbolic tensor name
        integer(INTD), intent(in):: hspaces(1:)              !in: registered id of the defining vector space for each tensor dimension
        integer(INTL), intent(in):: subspaces(1:)            !in: id of the defining subspace for each tensor dimension
        integer(INTD), intent(in):: data_kind                !in: data kind of tensor elements: {EXA_DATA_KIND_XX: XX = R4,R8,C4,C8}
        integer(INTL), intent(in), optional:: dim_extent(1:) !in: dimension extent for each tensor dimension (0 means deferred, to be set later)
        integer(INTD), intent(in), optional:: dim_group(1:)  !in: symmetric group (>=0) for each tensor dimension (0 means default)
        integer(INTD), intent(in), optional:: group_spec(1:) !in: symmetric group specification for non-trivial symmetric groups (see tensor_recursive.F90)
        logical, intent(in), optional:: sync                 !in: if FALSE, the operation will run asynchronously, requiring a separate exatns_sync() call later (defaults to TRUE)

        ierr=EXA_NOT_AVAILABLE
        return
      end function exatns_tensor_create_tensor

      function exatns_tensor_destroy(tensor,sync) result(ierr)
        !Destroys a tensor.
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !inout: tensor
        logical, intent(in), optional:: sync       !in: if FALSE, the operation will run asynchronously, requiring a separate exatns_sync() call later (defaults to TRUE)

        ierr=EXA_NOT_AVAILABLE
        return
      end function exatns_tensor_destroy


end module exatensor

