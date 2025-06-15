# definitions - from the CMake perspective this is not a good way to do this
# but this is kept here during iterative improvement of the code

if(MPI_COMPILER_MATCHES)
  set(USE_MPIF "0")
else()
  set(USE_MPIF "1")
endif()

if(ENABLE_64BIT_INTEGERS)
  set(INTEGER_PRECISION "INT64")
else()
  set(INTEGER_PRECISION "INT32")
endif()

set(HOST_PROGRAM "DIRAC")
if(ENABLE_GEN1INT)
  set(INTEGRAL_LIBRARY "GEN1INT")
else()
  message(FATAL_ERROR "-- PElib requires Gen1Int, use -DENABLE_GEN1INT=ON to enable Gen1Int or -DENABLE_PELIB=OFF to disable PElib")
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

#-------------------------------------------------------------------------------
# the interface

add_library(pelib_interface)

target_sources(
  pelib_interface
  PRIVATE
    interface/pelib_dirac_interface.F90
  )

target_include_directories(
  pelib_interface
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_compile_definitions(
  pelib_interface
  PUBLIC
    HAS_PELIB
  )

if(ENABLE_PDE)
  target_compile_definitions(
    pelib_interface
    PRIVATE
      VAR_PDE
    )
endif()

#-------------------------------------------------------------------------------
# the library

add_subdirectory(subtree)

target_link_libraries(
  pelib_interface
  PRIVATE
    PElib
  )

if(ENABLE_PDE)
  target_link_libraries(
    PElib
    PRIVATE
      hdf5::hdf5_fortran
    )

  target_compile_definitions(
    PElib
    PRIVATE
      VAR_PDE
    )
endif()

#-------------------------------------------------------------------------------
# definition of tests

dirac_test(pe_energy "medium;scf;pelib" "")
dirac_test(pe_response "medium;polarizabilities;pelib" "")
dirac_test(pe_exc "medium;pelib" "")
dirac_test(pe_cpp "short;pelib" "")
dirac_test(pe_solvat "short;pelib" "")
if(ENABLE_PDE)
  dirac_test(pde_exc "short;pelib" "")
endif()
