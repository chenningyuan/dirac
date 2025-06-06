option(ENABLE_PELIB "Enable polarizabile embedding library" ON)

if(ENABLE_PELIB)

    add_definitions(-DHAS_PELIB)
    set(PE_HOST_PROGRAM "DIRAC")
    if(ENABLE_GEN1INT)
        set(PE_INTEGRAL_LIBRARY "GEN1INT")
    else()
        message(FATAL_ERROR "-- PElib requires Gen1Int, use -DENABLE_GEN1INT=ON to enable Gen1Int or -DENABLE_PELIB=OFF to disable PElib")
    endif()

    set(PE_INCLUDE_DIRS)

    # here we control whether PElib will use mpi.mod or mpif.h
    # MPI_COMPILER_MATCHES is defined in cmake/custom/mpi.cmake
    # TODO: we should use generator expressions instead but this probably requires
    # moving away from add_external() which we should also do at some point
    if(MPI_COMPILER_MATCHES)
      set(_enable_mpif "0")
    else()
      set(_enable_mpif "1")
    endif()

    if(ENABLE_64BIT_INTEGERS)
      set(INTEGER_PRECISION "INT64")
    else()
      set(INTEGER_PRECISION "INT32")
    endif()

    # Include PDE by default, if HDF5 with Fortran component available
    find_package(HDF5 COMPONENTS Fortran)
    if(HDF5_FOUND)
       if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
          set(ENABLE_PDE on)
          message(STATUS "PDE feature included in PElib.")
       else()
       # libhdf5-fortran*.lib is (always?) compiled with gnu when installed on linux systems
       # if compiler is e.g. Intel, then the .mod files are not OK
          set(ENABLE_PDE off)
          message(STATUS "PDE feature not included in PElib, it is only available for GNU compilers.")
       endif()
    elseif(ENABLE_PDE)
       message(FATAL_ERROR "-- PDE feature in PElib requires HDF5 for Fortran. Either fix that or -DENABLE_PDE=OFF")
    else()
       message(STATUS "No PDE in PElib, a HDF5 library callable from Fortran is not available.")
    endif()
    if(ENABLE_PDE)
        add_definitions(-DVAR_PDE)
        find_package(HDF5 REQUIRED COMPONENTS Fortran)
        set(PE_INCLUDE_DIRS ${PE_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
        add_definitions(${HDF5_Fortran_DEFINITIONS})
    endif()

    set(ExternalProjectCMakeArgs
        -DCMAKE_WARN_DEPRECATED=FALSE
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_INSTALL_INCLUDEDIR=${PROJECT_BINARY_DIR}/modules
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
	      -DPARENT_INCLUDE_DIRS=${PE_INCLUDE_DIRS}
	      -DPARENT_MODULE_DIRS=${PROJECT_BINARY_DIR}/modules
        -DINTEGRAL_LIBRARY=${PE_INTEGRAL_LIBRARY}
	      -DENABLE_MPI=${ENABLE_MPI}
        -DINTEGER_PRECISION=${INTEGER_PRECISION}
        -DUSE_MPIF=${_enable_mpif}
        -DENABLE_PDE=${ENABLE_PDE}
        -DHOST_PROGRAM=${PE_HOST_PROGRAM}
    )

    unset(_enable_mpif)

    add_external(pelib)
    set(EXTERNAL_LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libPElib.a
        ${EXTERNAL_LIBS})
      if(ENABLE_PDE)
        set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${HDF5_Fortran_LIBRARIES})
    endif()
endif()

message(STATUS "PElib module: ${ENABLE_PELIB}")
