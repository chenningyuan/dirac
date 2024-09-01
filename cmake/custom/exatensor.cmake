#.rst:
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--exatensor=<ENABLE_EXATENSOR> Toggle use of ExaTensor <ON/OFF> [default: ON]."
#   define:
#     - "'-DENABLE_EXATENSOR={0}'.format(arguments['--exatensor'])"

#Check whether ExaTensor can be used
option(ENABLE_EXATENSOR "Enable ExaTENSOR library" ON)
option(TALSH_ONLY "Enable only TALSH component of ExaTENSOR library" OFF)

if(ENABLE_64BIT_INTEGERS)
   message(STATUS "ExaTensor library can not work with 64 bit integers, switched off")
   set(ENABLE_EXATENSOR OFF)
endif()

if(ENABLE_EXATENSOR)

   # Provide information about the location of the repository and the specific version (HASH) that should be used. 
   # This is configurable to make it possible to configure on machines
   # on which the repository resides locally but where network cannot be accessed
   # in which case EXATENSOR_GIT_REPO_LOCATION can be set to point to a path on the hard disk
   # ./setup --cmake-options="-DEXATENSOR_GIT_REPO_LOCATION='/path/to/exatensor/'"
   set (EXATENSOR_GIT_REPO_LOCATION "https://github.com/RelMBdev/ExaTENSOR.git" CACHE STRING "ExaTENSOR Git repository location")
   set (EXATENSOR_GIT_HASH 1a4ae3e4db8ae6843e27899e72cbc2de046498cc CACHE STRING "ExaTENSOR Git repository hash in use")
   set (EXATENSOR_INSTALL_DIR ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/lib)

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
	 set (EXATENSOR_CUDA_ARCH CUDA_ARCH=80)
         set (EXA_GPU_ENABLED TRUE)
      elseif ("$ENV{EXA_GPU}" MATCHES "NVIDIA")
         find_package(CUDAToolkit QUIET)
         if (DEFINED ENV{GPU_CUDA})
            set (EXATENSOR_GPUENV GPU_CUDA=$ENV{GPU_CUDA})
            set (EXA_GPU_ENABLED TRUE)
         elseif (CUDAToolkit_FOUND)
            set (EXATENSOR_GPUENV GPU_CUDA=CUDA)
            set (EXA_GPU_ENABLED TRUE)
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
   elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
      set (EXATENSOR_TOOLKIT TOOLKIT=PGI)
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
   set(EXATENSOR_ENV "${EXATENSOR_WRAP} ${EXATENSOR_BUILD} ${EXATENSOR_TALSH} ${TALSH_COMPILERS} ${EXATENSOR_TOOLKIT} ${EXATENSOR_OS} ${EXATENSOR_GPUENV} ${EXATENSOR_MPIENV} ${EXATENSOR_BLAS} ${EXATENSOR_ROCM} ${EXATENSOR_CUDA_ARCH} ${EXATENSOR_AMD_ARCH} ${EXATENSOR_HIP}")
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
    dirac_test(exacorr_symmstart "cc;talsh;short" "")
    dirac_test(exacorr_talsh_memory "cc;talsh;short" "")
    dirac_test(exacorr_talsh_debug "cc;talsh;short" "")
    dirac_test(exacorr_talsh_nr_contracted "cc;talsh;short" "")
    dirac_test(exacorr_talsh_standalone "cc;talsh;short" "")
    dirac_test(exacorr_talsh_lambda "cc;talsh;short" "")
    dirac_test(exacorr_talsh_lambda_2rdm "cc;talsh;short" "")
    dirac_test(exacorr_talsh_open_lambda "cc;talsh;medium" "")
    dirac_test(exacorr_talsh_fock "cc;talsh;short" "")
    dirac_test(exacorr_talsh_open "cc;talsh;short" "")
    dirac_test(exacorr_talsh_eomip "cc;exacorr" "")
    dirac_test(exacorr_talsh_eomip_full "cc;exacorr" "")
    dirac_test(exacorr_talsh_eomee "cc;exacorr" "")
    dirac_test(exacorr_talsh_eomee_full "cc;exacorr" "")
    dirac_test(exacorr_talsh_eomea "cc;exacorr;short" "")
    dirac_test(exacorr_talsh_eomea_full "cc;exacorr" "")
    dirac_test(exacorr_talsh_cis "cis;talsh;short" "")
    dirac_test(exacorr_talsh_ccrsp_lr_static "ccrsp;talsh" "")
    dirac_test(exacorr_talsh_ccrsp_lr_cc-cc "ccrsp;talsh" "")
    dirac_test(exacorr_talsh_ccrsp_lr_cc-ci "ccrsp;talsh" "")
    dirac_test(exacorr_talsh_ccrsp_cpp "ccrsp;talsh" "")
    dirac_test(exacorr_talsh_ccrsp_spsp "ccrsp;talsh" "")
    dirac_test(exacorr_talsh_ccrsp_qr_static "ccrsp;talsh" "")
    dirac_test(exacorr_talsh_ccrsp_qr_cc-ci "ccrsp;talsh;short" "")
    dirac_test(exacorr_talsh_ccrsp_eomtpa "ccrsp;talsh" "")
    dirac_test(exacorr_talsh_cc2 "cc;talsh;medium" "")
    dirac_test(exacorr_talsh_spinfree "cc;talsh;medium" "")
    dirac_test(exacorr_talsh_respect "cc;talsh;respect;medium" "")
    dirac_test(exacorr_talsh_finite_field "cc;talsh;short" "")
    dirac_test(exacorr_talsh_cholesky "cc;talsh;long" "")
    dirac_test(exacorr_talsh_mp2no "mp2no;talsh;long" "")
    dirac_test(exacorr_talsh_mp2no_restart "mp2no_re;talsh;long" "")
    #dirac_test(exacorr_talsh_lambda_cc2 "cc;talsh" "")
    if(ENABLE_LAPLACE)
      dirac_test(exacorr_talsh_mp2lap "mp2lap;talsh" "")
      dirac_test(exacorr_talsh_tripl_lap "cc;talsh" "")
    endif()
    dirac_test(exacorr_talsh_lambda_mbauto "cc;talsh" "")
    dirac_test(exacorr_talsh_ccsdt_mbauto "cc;exacorr;talsh;mbauto" "")
    dirac_test(exacorr_talsh_ccsdtq_mbauto "cc;exacorr;talsh;mbauto" "")
# this test is currently very long so disable it unless running on GPU nodes
    if (EXA_GPU_ENABLED)
       dirac_test(exacorr_talsh_ccsdtqp_mbauto "cc;exacorr;talsh;mbauto;long" "")
    endif()
    dirac_test(exacorr_talsh_ccsdt_tenpi "cc;exacorr;talsh;tenpi;ccsdt" "")
    dirac_test(exacorr_talsh_ccsdtq_tenpi "cc;exacorr;talsh;tenpi;ccsdtq" "")

   if (ENABLE_MPI AND NOT EXATENSOR_TALSH)
        dirac_test(exacorr_exatensor_memory "cc;exatensor;short" "")
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
        dirac_test(exacorr_exatensor_lambda_mbauto "cc;exatensor;mbauto" "")
        dirac_test(exacorr_exatensor_ccsdt_mbauto "cc;exatensor;mbauto" "")
        dirac_test(exacorr_exatensor_ccsdtq_mbauto "cc;exatensor;mbauto" "")
# this test is currently very long so disable it unless running on GPU nodes
        if (EXA_GPU_ENABLED)
           dirac_test(exacorr_exatensor_ccsdtqp_mbauto "cc;exatensor;mbauto;long" "")
        endif()
        dirac_test(exacorr_exatensor_ccsdt_tenpi "cc;exatensor;tenpi;ccsdt" "")
        dirac_test(exacorr_exatensor_ccsdtq_tenpi "cc;exatensor;tenpi;ccsdtq" "")
    endif()

endif()

message(STATUS "Enable ExaTENSOR library: ${ENABLE_EXATENSOR}")
