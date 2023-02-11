# generates git information
add_custom_command(
  OUTPUT
    ${CMAKE_CURRENT_BINARY_DIR}/git_info.h
    ALL
  COMMAND
    ${CMAKE_COMMAND} -D INPUT_DIR=${PROJECT_SOURCE_DIR}/cmake/downloaded -D OUTPUT_DIR=${PROJECT_BINARY_DIR} -P ${PROJECT_SOURCE_DIR}/cmake/downloaded/git-hash.cmake
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_custom_target(git_info
  ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/git_info.h
  )

#-------------------------------------------------------------------------------

add_library(quaternion_algebra)

target_sources(quaternion_algebra
  PRIVATE
    quaternion_algebra.F90
  )

#-------------------------------------------------------------------------------

add_library(hdf5_util)

option(ENABLE_HDF5 "Enable HDF5" OFF)

if(ENABLE_HDF5)
  find_package(HDF5 COMPONENTS C REQUIRED)

  target_sources(hdf5_util
    PRIVATE
     hdf5_util/mh5.F90
     hdf5_util/mh5_capi.c
   )

  target_link_libraries(hdf5_util
    PRIVATE
      HDF5::HDF5
    )
else()
  target_sources(hdf5_util
    PRIVATE
     hdf5_util/mh5_nohdf5.F90
   )
endif()

#-------------------------------------------------------------------------------

add_library(gp)

target_sources(gp
  PRIVATE
    execinfo_bsd.c
    gpc.c
    memory_sys_if.c
    print_info.c
    selftest_fortran_c_interoperability.c
    stacktraverse_bsd.c
    trace.c
    finesort_scheme6.F
    global_file.F
    gphjj.F
    gpjeppe.F
    gpjth.F
    gpkrmc.F
    gplupar.F
    gpluuk.F
    gpsaue.F
    gptrygve.F
    gpulfek.F
    graphic.F
    io.F
    matexport.F
    math.F
    mempkg.F
    mpi2macro.F
    mxmpack.F
    order2.F
    qalgebra.F
    qpack.F
    spline.F
    time.F
    gplucia_io.F
    gplucia.F
    basis_datatypes.F90
    checkpoint.F90
    codata.F90
    debug_mem_alloc.F90
    dirac_cfg.F90
    dirac_mocoefficients.F90
    dynamic_string.F90
    electrostatic_potential.F90
    # FIXME
    quaternion_algebra.F90
    exception_handling.F90
    file_handler.F90
    f77_work_interface.F90
    gp_q_diag.F90
    gpsk.F90
    heapsort.F90
    include_h_modules.F90
    kinds.F90
    luci_wrkspc.F90
    machine_parameters.F90
    memory_allocator.F90
    memory_control.F90
    memory_errorhandler.F90
    memory_internal_special_i8.F90
    memory_internal_ada_1d.F90
    memory_internal_ada_2d.F90
    memory_internal_ada_3d.F90
    memory_internal_ada_4d.F90
    memory_internal_init.F90
    memory_legacy_dirac.F90
    memory_parameters.F90
    memory_track_if.F90
    memory_tracker.F90
    moltra_labeling.F90
    num_der.F90
    orbital_rotation_indices.F90
    osutils.F90
    picture_change_operations.F90
    selftests.F90
    string_manipulations.f90
    type_laplace.F90
    unit_test_generator.f90
    unit_testing.F90
    vector_functions.F90
    xmlout.F90
    xml/xml_file.F90
    xml/xml_structure.F90
    xml/xml_parser.F90
    labeled_storage.F90
    dircmo.F90
  )

target_include_directories(gp
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_BINARY_DIR}
  )

target_link_libraries(gp
  PRIVATE
    hdf5_util
#   quaternion_algebra
    interface_mpi
    pdpack
    hergp
  )

#-------------------------------------------------------------------------------

# FIXME move this somewhere else and adapt for internal math
find_package(BLAS)

if(BLAS_FOUND)
  target_link_libraries(gp
    PRIVATE
      BLAS::BLAS
    )
endif()

add_dependencies(gp git_info)
