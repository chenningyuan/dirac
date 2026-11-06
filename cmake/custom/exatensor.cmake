#.rst:
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--exatensor=<ENABLE_EXATENSOR> Toggle use of ExaTensor <ON/OFF> [default: ON]."
#   define:
#     - "'-DENABLE_EXATENSOR={0}'.format(arguments['--exatensor'])"

option(ENABLE_EXATENSOR "Enable ExaTENSOR library" ON)
option(TALSH_ONLY "Enable only TALSH component of ExaTENSOR library" OFF)

#Check whether ExaTensor can be used (if enabled)
if(ENABLE_EXATENSOR)
   if(ENABLE_64BIT_INTEGERS)
      message(STATUS "ExaTensor library cannot work with 64 bit integers, switched off")
      set(ENABLE_EXATENSOR OFF)
   endif()

endif()

if(ENABLE_EXATENSOR)

   # Provide information about the location of the repository and the specific version (HASH) that should be used. 
   # This is configurable to make it possible to configure on machines
   # on which the repository resides locally but where network cannot be accessed
   # in which case EXATENSOR_GIT_REPO_LOCATION can be set to point to a path on the hard disk
   # ./setup --cmake-options="-DEXATENSOR_GIT_REPO_LOCATION='/path/to/exatensor/'"
   set (EXATENSOR_GIT_REPO_LOCATION "https://github.com/RelMBdev/ExaTENSOR.git" CACHE STRING "ExaTENSOR Git repository location")
   set (EXATENSOR_GIT_HASH a1efa9efdcae13caacdc4776af839c3cd1957e85 CACHE STRING "ExaTENSOR Git repository hash in use")
   set (EXATENSOR_INSTALL_DIR ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/lib)

   # user chooses whether to install specific tensor libraries used through TAPP
   if(NOT DEFINED ENABLE_TBLIS)
       set(_val "$ENV{ENABLE_TBLIS}")
       if(_val STREQUAL "")
           set(_val "OFF")
       endif()
       set(ENABLE_TBLIS "${_val}" CACHE BOOL "Build and link TBLIS and TBLIS bindings")
   endif()

   if(DEFINED ENV{ENABLE_TBLIS} AND NOT "${ENABLE_TBLIS}" STREQUAL "$ENV{ENABLE_TBLIS}")
      set(ENABLE_TBLIS "$ENV{ENABLE_TBLIS}" CACHE BOOL "" FORCE)
   endif()

   # set default executor (1=tapp with tblis, 2=talsh, 3=exatensor, 4=tapp reference implementation) 
   # of tensor operations in exacorr_tensor_interface
   if(NOT DEFINED TENSOR_EXECUTOR)
       set(_val "$ENV{TENSOR_EXECUTOR}")
       if(_val STREQUAL "")
           set(_val "2")
       endif()
       set(TENSOR_EXECUTOR "${_val}" CACHE STRING "Tensor lib.: 1=tblis, 2=talsh, 3=exatensor, 4=ref. impl.")
   endif()

   if(DEFINED ENV{TENSOR_EXECUTOR} AND NOT "${TENSOR_EXECUTOR}" STREQUAL "$ENV{TENSOR_EXECUTOR}")
      set(TENSOR_EXECUTOR "$ENV{TENSOR_EXECUTOR}" CACHE STRING "" FORCE)
   endif()

   if (TENSOR_EXECUTOR STREQUAL "1" AND ENABLE_TBLIS)
     add_definitions(-DTENSOR_EXECUTOR_=1)
     add_definitions(-DDEFAULT_EXECUTOR_WITHIN_TAPP=2)
   elseif(TENSOR_EXECUTOR STREQUAL "4")
     add_definitions(-DTENSOR_EXECUTOR_=1)
     add_definitions(-DDEFAULT_EXECUTOR_WITHIN_TAPP=1) # reference implementation if tblis not there
   else()
     add_definitions(-DTENSOR_EXECUTOR_=${TENSOR_EXECUTOR})
     add_definitions(-DDEFAULT_EXECUTOR_WITHIN_TAPP=1)
   endif()


   include(FetchContent)

   FetchContent_Declare(
     tapp
     GIT_REPOSITORY https://github.com/TAPPorg/reference-implementation.git
     GIT_TAG        55e5b1e
   )

   FetchContent_MakeAvailable(tapp)




   # check build type
   if (DEFINED ENV{BUILD_TYPE})
	   set(EXATENSOR_BUILD BUILD_TYPE=$ENV{BUILD_TYPE})
   else()
	   set(EXATENSOR_BUILD BUILD_TYPE=OPT)
   endif()

   # check wrapping
   if (DEFINED ENV{WRAP})
	   set(EXATENSOR_WRAP WRAP=$ENV{WRAP})
   else()
	   set(EXATENSOR_WRAP WRAP=NOWRAP)
   endif()

   if (DEFINED ENV{WRAP})
      set(EXATENSOR_WRAP WRAP=$ENV{WRAP})
   else()
      set(EXATENSOR_WRAP WRAP=NOWRAP)
   endif()

   # The supported operating systems for ExaTensor are Linux(default) and MacOSX
   if (DEFINED ENV{EXA_OS})
      set (EXATENSOR_OS EXA_OS=$ENV{EXA_OS})
   elseif (CMAKE_HOST_APPLE) 
      set (EXATENSOR_OS EXA_OS=NO_LINUX)
   else()
      set (EXATENSOR_OS EXA_OS=LINUX)
   endif()

   # Check whether we can use GPUs (may need to be more precise as to minimum requirements)

   set (EXA_GPU_ENABLED FALSE)
   if (DEFINED ENV{EXA_GPU})
      if ("$ENV{EXA_GPU}" MATCHES "AMD")
         set (EXATENSOR_HIP USE_HIP=YES)
         set (EXATENSOR_GPUENV GPU_CUDA=CUDA)
	 set (EXATENSOR_ROCM PATH_ROCM=$ENV{ROCM_PATH})
         # cuda_arch will not matter for AMD gpus, TAL-SH will handle that
         # but we set it anyway now 
	 set (EXATENSOR_CUDA_ARCH GPU_SM_ARCH=80)
         set (EXA_GPU_ENABLED TRUE)
      elseif ("$ENV{EXA_GPU}" MATCHES "NVIDIA")
         find_package(CUDAToolkit QUIET)
         if (DEFINED ENV{GPU_CUDA})
            set (EXATENSOR_GPUENV GPU_CUDA=$ENV{GPU_CUDA})
            set (EXA_GPU_ENABLED TRUE)
            set (EXATENSOR_CUDA_ARCH GPU_SM_ARCH=$ENV{CUDA_SM_ARCH})
	    set (EXATENSOR_CUDA PATH_CUDA=$ENV{CUDA_ROOT})
         elseif (CUDAToolkit_FOUND)
            set (EXATENSOR_GPUENV GPU_CUDA=CUDA)
            set (EXA_GPU_ENABLED TRUE)
            set (EXATENSOR_CUDA_ARCH GPU_SM_ARCH=$ENV{CUDA_SM_ARCH})
	    set (EXATENSOR_CUDA PATH_CUDA=$ENV{CUDA_ROOT})
         else()
            set (EXATENSOR_GPUENV GPU_CUDA=NOCUDA)
         endif()
      else()
         message(STATUS "WARNING: unsupported GPU architecture. Disabling GPUs in ExaCorr.") 
         set (EXATENSOR_GPUENV GPU_CUDA=NOCUDA)
      endif()
   else()
      message(STATUS "WARNING: Assuming no GPUs present for ExaCorr.") 
      message(STATUS "  If that is incorrect define shell variable EXA_GPU to AMD or NVIDIA and reconfigure.")
      set (EXATENSOR_GPUENV GPU_CUDA=NOCUDA)
   endif()

   # check if TALSH is requested
   if ($ENV{TALSH_ONLY} MATCHES "YES") 
      set(EXATENSOR_TALSH EXA_TALSH_ONLY=YES)
      set(EXATENSOR_MPIENV MPILIB=NONE)
   # check if MPI is available
   elseif(NOT ENABLE_MPI)
      set(EXATENSOR_TALSH EXA_TALSH_ONLY=YES)
      set(EXATENSOR_MPIENV MPILIB=NONE)
   # get MPI information
   else()
       if (DEFINED ENV{MPILIB})
         set(EXATENSOR_MPIENV MPILIB=$ENV{MPILIB})
       else()
         # get MPI TYPE from Python script
         execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cmake/custom/IdentifyMPI.py
                                 ${MPIEXEC} OUTPUT_VARIABLE PYTHON_OUTPUT)
         get_filename_component(MPIDIR ${MPIEXEC} DIRECTORY)
         get_filename_component(PATH_MPI ${MPIDIR} DIRECTORY)
         string(STRIP "${PYTHON_OUTPUT}" PYTHON_OUTPUT)
         if (${PYTHON_OUTPUT} MATCHES "MPICH")
           set(EXATENSOR_MPIENV "MPILIB=MPICH PATH_MPICH=${PATH_MPI}")
         elseif (${PYTHON_OUTPUT} MATCHES "OPENMPI")
           set(EXATENSOR_MPIENV "MPILIB=OPENMPI PATH_OPENMPI=${PATH_MPI}")
         elseif (${PYTHON_OUTPUT} MATCHES "INTEL")
           set(EXATENSOR_MPIENV "MPILIB=MPICH PATH_MPICH=${PATH_MPI}")
         else()
           message(STATUS "WARNING: The MPI compiler ${PYTHON_OUTPUT} is not supported by ExaTENSOR")
           message(STATUS "Only MPICH (3.2.1+) and OpenMPI (4.1.1+) are currently supported by ExaTENSOR")
           set(EXATENSOR_TALSH EXA_TALSH_ONLY=YES)
           set(EXATENSOR_MPIENV MPILIB=NONE)
         endif()
       endif()
       if (DEFINED ENV{INTEL_MPI_WRAPPER_FC})
          set(EXATENSOR_MPIENV "${EXATENSOR_MPIENV} INTEL_MPI_WRAPPER_FC=$ENV{INTEL_MPI_WRAPPER_FC}") 
       endif()
       if (DEFINED ENV{INTEL_MPI_WRAPPER_CC})
          set(EXATENSOR_MPIENV "${EXATENSOR_MPIENV} INTEL_MPI_WRAPPER_CC=$ENV{INTEL_MPI_WRAPPER_CC}") 
       endif()
       if (DEFINED ENV{INTEL_MPI_WRAPPER_CPP})
          set(EXATENSOR_MPIENV "${EXATENSOR_MPIENV} INTEL_MPI_WRAPPER_CPP=$ENV{INTEL_MPI_WRAPPER_CPP}") 
       endif()
   endif()

   # Find out which compiler family we are using (this can be tricky, tweak after the configure step if necessary)
   if (DEFINED ENV{TOOLKIT})
      set(EXATENSOR_TOOLKIT TOOLKIT=$ENV{TOOLKIT})
   elseif (CMAKE_Fortran_COMPILER_ID MATCHES GNU) 
      set (EXATENSOR_TOOLKIT TOOLKIT=GNU)
      if (   "${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS "8"
         OR "${CMAKE_Fortran_COMPILER_VERSION}" VERSION_EQUAL "9"
         OR "${CMAKE_Fortran_COMPILER_VERSION}" VERSION_EQUAL "10" )
         message(STATUS "WARNING: GNU compiler ${CMAKE_Fortran_COMPILER_VERSION} is not supported by ExaTENSOR") 
         set(EXATENSOR_TALSH EXA_TALSH_ONLY=YES)
      endif()
   elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
      set (EXATENSOR_TOOLKIT TOOLKIT=INTEL)
      if ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS "18")
         message(STATUS "WARNING: INTEL compiler ${CMAKE_Fortran_COMPILER_VERSION} is not supported by ExaTENSOR")
         set(EXATENSOR_TALSH EXA_TALSH_ONLY=YES)
      endif()
   elseif(CMAKE_Fortran_COMPILER_ID MATCHES Cray)
      set (EXATENSOR_TOOLKIT TOOLKIT=CRAY)
      set(EXATENSOR_TALSH EXA_TALSH_ONLY=YES)
   elseif(CMAKE_Fortran_COMPILER_ID MATCHES Flang)
      set (EXATENSOR_TOOLKIT TOOLKIT=CRAY)
   elseif(CMAKE_Fortran_COMPILER_ID MATCHES XL)
      set (EXATENSOR_TOOLKIT TOOLKIT=IBM)
   elseif(CMAKE_Fortran_COMPILER_ID MATCHES NVHPC)
      set (EXATENSOR_TOOLKIT TOOLKIT=NVHPC)
   else()
      message(STATUS "WARNING: ExaTensor TOOLKIT not found, setting to GNU")
      set (EXATENSOR_TOOLKIT TOOLKIT=GNU)
      message(STATUS "WARNING: ExaTensor will probably not work")
   endif()

   #set environment needed for TALSH
   if (EXATENSOR_TALSH) 
      add_definitions(-DEXA_TALSH_ONLY)
      message(STATUS "Enabling only the TALSH (serial) component of the ExaTENSOR library")
      # The following weird construction is needed because these cmake variables are otherwise not active inside the environment and make that is called below
      set(TALSH_COMPILERS "CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} CMAKE_C_COMPILER=${CMAKE_C_COMPILER} CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
   endif()
  
   # We need linking to the math library for cholesky
   if (DEFINED ENV{BLASLIB})
	   set(EXATENSOR_BLAS BLASLIB=$ENV{BLASLIB})
   else()
      if ( BLAS_TYPE MATCHES ATLAS )
	      set(EXATENSOR_BLAS BLASLIB=ATLAS)
	      string(REPLACE "." ";" EXATENSOR_BLAS_PATH ${BLAS_LIBRARIES})
	      list (GET EXATENSOR_BLAS_PATH 0 EXATENSOR_TEMP)
	      get_filename_component(EXATENSOR_BLAS_PATH ${EXATENSOR_TEMP} PATH)
	      set(EXATENSOR_BLAS "${EXATENSOR_BLAS} PATH_BLAS_ATLAS=${EXATENSOR_BLAS_PATH}")
      elseif( BLAS_TYPE MATCHES OPENBLAS )
         set(EXATENSOR_BLAS BLASLIB=OPENBLAS)
	      string(REPLACE "." ";" EXATENSOR_BLAS_PATH ${BLAS_LIBRARIES})
	      list (GET EXATENSOR_BLAS_PATH 0 EXATENSOR_TEMP)
	      get_filename_component(EXATENSOR_BLAS_PATH ${EXATENSOR_TEMP} PATH)
	      set(EXATENSOR_BLAS "${EXATENSOR_BLAS} PATH_BLAS_OPENBLAS=${EXATENSOR_BLAS_PATH}")
      elseif( BLAS_TYPE MATCHES LIBSCI )
         set(EXATENSOR_BLAS BLASLIB=LIBSCI)
         string(REPLACE "." ";" EXATENSOR_BLAS_PATH ${BLAS_LIBRARIES})
         list (GET EXATENSOR_BLAS_PATH 0 EXATENSOR_TEMP)	 
	 get_filename_component(EXATENSOR_BLAS_PATH ${EXATENSOR_TEMP} PATH)
	 set(EXATENSOR_BLAS "${EXATENSOR_BLAS} PATH_BLAS_LIBSCI=${EXATENSOR_BLAS_PATH}")
      elseif( BLAS_TYPE MATCHES MKL )
	      set(EXATENSOR_BLAS BLASLIB=MKL)
      elseif( BLAS_TYPE MATCHES ESSL )
	      set(EXATENSOR_BLAS BLASLIB=ESSL)
      elseif( BLAS_TYPE MATCHES ACML )
         set(EXATENSOR_BLAS BLASLIB=ACML)
      else()
         set(EXATENSOR_BLAS BLASLIB=NONE)
	      message(STATUS "WARNING: No BLAS library for TALSH / EXATENSOR - certain functionalties will not work (cholesky) ")
      endif()
   endif()

   # Collect everything in one string (also including the hardwires ones that need not be changed) and store this in a file (direct passing appears to be impossible within cmake)
   set(EXATENSOR_ENV "${EXATENSOR_WRAP} ${EXATENSOR_BUILD} ${EXATENSOR_TALSH} ${TALSH_COMPILERS} ${EXATENSOR_TOOLKIT} ${EXATENSOR_OS} ${EXATENSOR_GPUENV} ${EXATENSOR_MPIENV} ${EXATENSOR_BLAS} ${EXATENSOR_CUDA} ${EXATENSOR_ROCM} ${EXATENSOR_CUDA_ARCH} ${EXATENSOR_AMD_ARCH} ${EXATENSOR_HIP}")
   message(STATUS "The environment variables used to build ExaTensor are collected in the file ExaTensor_ENV (can be inspected/changed if necessary)")
   file(WRITE ${PROJECT_BINARY_DIR}/ExaTensor_ENV ${EXATENSOR_ENV})
   file(WRITE ${PROJECT_BINARY_DIR}/ExaTensor_ENV_UP "${EXATENSOR_ENV}  EXA_NO_BUILD=YES")

   ExternalProject_Add(exatensor
        PREFIX "${PROJECT_BINARY_DIR}/exatensor"
        GIT_REPOSITORY ${EXATENSOR_GIT_REPO_LOCATION}
        GIT_TAG ${EXATENSOR_GIT_HASH}
        GIT_CONFIG advice.detachedHead=false
        CONFIGURE_COMMAND true  # currently no configure command, but this is needed for cmake to function
        BUILD_COMMAND cd ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/ && set -a && . ${PROJECT_BINARY_DIR}/ExaTensor_ENV && set +a && make && cp ${PROJECT_BINARY_DIR}/ExaTensor_ENV_UP ${PROJECT_BINARY_DIR}/ExaTensor_ENV
        INSTALL_DIR ${EXATENSOR_INSTALL_DIR}
        INSTALL_COMMAND true
        )

   set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${EXATENSOR_INSTALL_DIR}/libtalsh.a)
   if (ENABLE_MPI AND NOT EXATENSOR_TALSH)
       set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${EXATENSOR_INSTALL_DIR}/libexatensor.a)
   endif()
   
   #In case of problems with the build, one may manually add ExaTENSOR dependencies here (BLAS, OpenMP, CUDA, C++):
   #set(EXTERNAL_LIBS ${EXTERNAL_LIBS} -L/sw/summit/essl/6.1.0-2/essl/6.1/lib64 -L/sw/summit/xl/16.1.1-5/xlC/16.1.1/lib -L/sw/summit/xl/16.1.1-5/xlf/16.1.1/lib -lessl -lxlf90_r -lxlfmath -L/sw/summit/cuda/10.1.243/lib64 -lcublas -lcudart -lnvToolsExt -lstdc++ -lgomp)
   # note: the same line could be added to CMAKE_EXE_LINKER_FLAGS

   include_directories(${PROJECT_BINARY_DIR}/exatensor/src/exatensor/include)

   #Add also the tests (weird to do this here, but this whole test set up needs an overhaul).
   #disabled some tests waiting for code verification/efficiency improvement
    dirac_test(exacorr_ccsdt_mbauto "cc;exacorr;mbauto;long" "")
    if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC") # nvfortran does not accept arrays with 8 or more indices
       dirac_test(exacorr_ccsdtq_tenpi "cc;exacorr;tenpi;ccsdtq;long" "")
       dirac_test(exacorr_ccsdtq_mbauto "cc;exacorr;mbauto;long" "")
    if (EXA_GPU_ENABLED)
# this test is currently very long so disable it unless running on GPU nodes
       dirac_test(exacorr_ccsdtqp_mbauto "cc;exacorr;mbauto;long" "")
    endif()
    endif()
    dirac_test(exacorr_ccrsp_eomtpa "ccrsp;exacorr;long" "")
    dirac_test(exacorr_lambda "cc;exacorr;long" "")
    dirac_test(exacorr_ccrsp_cpp "ccrsp;exacorr;long" "")
    dirac_test(exacorr_ccrsp_qr_static "ccrsp;exacorr;long" "")
    dirac_test(exacorr_cholesky "cc;exacorr;long" "")
    dirac_test(exacorr_mp2no "mp2no;exacorr;long" "")
    dirac_test(exacorr_mp2no_restart "mp2no_re;exacorr;long" "")
# disabling mp2neo for now, until further testing
#    dirac_test(exacorr_mp2neo "mp2neo;exacorr;long" "")
    dirac_test(exacorr_mp2fno++ "mp2fno++;exacorr;long" "")
    #dirac_test(exacorr_lambda_cc2 "cc;exacorr" "")
    dirac_test(exacorr_mp2lap "mp2lap;exacorr;long" "")
    dirac_test(exacorr_tripl_lap "cc;exacorr;long" "")
# reactivate this when lambda is handled by mbauto; autokern for lambda is tested with lambda_tenpi
#   dirac_test(exacorr_lambda_mbauto "cc;exacorr;medium" "")
    dirac_test(exacorr_lambda_tapp "cc;exacorr;tapp;medium" "")
    dirac_test(exacorr_cc2 "cc;exacorr;medium" "")
    dirac_test(exacorr_spinfree "cc;exacorr;medium" "")
    dirac_test(exacorr_respect "cc;exacorr;respect;medium" "")
    dirac_test(exacorr_ccrsp_spsp "ccrsp;exacorr;medium" "")
    dirac_test(exacorr_open_lambda "cc;exacorr;medium" "")
    dirac_test(exacorr_open_ccrsp_lr_cc-cc "cc;exacorr;medium" "")
    dirac_test(exacorr_ccrsp_lr_cc-cc "ccrsp;exacorr;medium" "")
    dirac_test(exacorr_ccrsp_qr_spsp_pol "ccrsp;exacorr;medium" "")
    dirac_test(exacorr_symmstart "cc;exacorr;short" "")
    dirac_test(exacorr_debug "cc;exacorr;short" "")
    dirac_test(exacorr_nr_contracted "cc;exacorr;medium" "")
    dirac_test(exacorr_standalone "cc;exacorr;medium" "")
    dirac_test(exacorr_lambda_2rdm "cc;exacorr;medium" "")
    dirac_test(exacorr_fock "cc;exacorr;medium" "")
    dirac_test(exacorr_open "cc;exacorr;short" "")
    dirac_test(exacorr_eomip "cc;exacorr;short" "")
    dirac_test(exacorr_eomip_full "cc;exacorr;medium" "")
    dirac_test(exacorr_eomee "cc;exacorr;short" "")
    dirac_test(exacorr_eomee_full "cc;exacorr;short" "")
    dirac_test(exacorr_eomea "cc;exacorr;short" "")
    dirac_test(exacorr_eomea_full "cc;exacorr;short" "")
    dirac_test(exacorr_cis "cis;exacorr;short" "")
    dirac_test(exacorr_cis_d "cis_d;exacorr;short" "")
# disabling for now, until further testing
#   dirac_test(exacorr_rt_cis "rt_cis;exacorr;short" "")
    dirac_test(exacorr_ccrsp_lr_static "ccrsp;exacorr;short" "")
    dirac_test(exacorr_ccrsp_lr_cc-ci "ccrsp;exacorr;short" "")
    dirac_test(exacorr_ccrsp_qr_cc-ci "ccrsp;exacorr;short" "")
    dirac_test(exacorr_ccrsp_qr_cc-cc "ccrsp;exacorr;short" "")
    dirac_test(exacorr_finite_field "cc;exacorr;medium" "")
    dirac_test(exacorr_ccsdt_tenpi "cc;exacorr;tenpi;ccsdt;short" "")
    dirac_test(exacorr_ccsdt_tenpi_tapp "cc;exacorr;tenpi:tapp;ccsdt;short" "")
    dirac_test(exacorr_ccsdt_lambda_tenpi_tapp "cc;exacorr;tenpi:tapp;ccsdt;lambda;short" "")
   if (ENABLE_MPI AND NOT EXATENSOR_TALSH)
        if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
        if (EXA_GPU_ENABLED)
# this test is currently very long so disable it unless running on GPU nodes
           dirac_test(exacorr_exatensor_ccsdtqp_mbauto "cc;exatensor;mbauto;long" "")
        endif()
           dirac_test(exacorr_exatensor_ccsdtq_mbauto "cc;exatensor;mbauto" "")
           dirac_test(exacorr_exatensor_ccsdtq_tenpi "cc;exatensor;tenpi;ccsdtq" "")
        endif()
        dirac_test(exacorr_exatensor_debug "cc;exatensor" "")
        dirac_test(exacorr_exatensor_nr_contracted "cc;exatensor" "")
        dirac_test(exacorr_exatensor_fock "cc;exatensor" "")
        dirac_test(exacorr_exatensor_open "cc;exatensor" "")
        dirac_test(exacorr_exatensor_lambda "cc;exatensor" "")
        dirac_test(exacorr_exatensor_cc2 "cc;exatensor" "")
        dirac_test(exacorr_exatensor_spinfree "cc;exatensor" "")
        dirac_test(exacorr_exatensor_ao2mo "cc;exatensor;long" "")
        dirac_test(exacorr_exatensor_finite_field "cc;exatensor" "")
        #works 05.2021, but maximum search is very slow
        #dirac_test(exacorr_exatensor_cholesky "cc;exatensor;long" "")
        dirac_test(exacorr_exatensor_mp2no "mp2no;exatensor" "")
        dirac_test(exacorr_exatensor_mp2no_restart "mp2no_re;exatensor" "")
# unlike the tal-sh test, this exatensor_lambda_mbauto test does not yet have the
# automatic code generation for lambda integrated. once that's done, this test should
# also be deactivated unless mbauto lambda is also implemented
        dirac_test(exacorr_exatensor_lambda_mbauto "cc;exatensor;mbauto" "")
        dirac_test(exacorr_exatensor_ccsdt_mbauto "cc;exatensor;mbauto" "")
        dirac_test(exacorr_exatensor_ccsdt_tenpi "cc;exatensor;tenpi;ccsdt" "")
    endif()

endif()

# redefine ENABLE_EXATENSOR with current value in CmakeCache.txt, used by e.g. ccmake and cmake-gui
# otherwise it would still apear as ON to ccmake .., even if it has been set OFF above
set(ENABLE_EXATENSOR ${ENABLE_EXATENSOR} CACHE BOOL "Enable ExaTENSOR library" FORCE)
set(TALSH_ONLY ${TALSH_ONLY} CACHE BOOL "Enable only TALSH component of ExaTENSOR library" FORCE)

message(STATUS "Enable ExaTENSOR library: ${ENABLE_EXATENSOR}")
