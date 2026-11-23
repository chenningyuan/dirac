find_package(HDF5 COMPONENTS C) 

add_library(
  rose_interest
    interest/module_interest_one.f90
  )

# the auxiliaries library contains modules that are independent of ROSE and are (or can be) used elsewhere as well
add_library(
  rose_auxiliaries
    # module to set some globals like output file units
    rose_global.F90
    # module to support various types of (quaternion/complex/real) linear algebra operations
    linear_algebra.F90
    # module to work with arrays with variably sized rows
    jagged_array.F90
    # module to generate coefficients from cartesian to spherical transformation (taken from DIRAC)
    cartesian_to_spherical.F90
  )

add_library(
  rose_ams_placeholders
    # Placeholders for routines of CoreKF in case we can not link to AMS at all
    ams_placeholders.F90
  )

if(CMAKE_PROJECT_NAME STREQUAL "rose")
  target_sources(
    rose_auxiliaries
    PRIVATE
      # module to read xyz input and store this as a linked list
      read_xyzfile.F90
      # labeled storage support, provided by DIRAC's gp library in the full build
      labeled_storage.F90
  )
  if(ENABLE_HDF5)
    target_sources(
      rose_auxiliaries
      PRIVATE
        hdf5_util/mh5.F90
        hdf5_util/mh5_capi.c
    )
  else()
    target_sources(
      rose_auxiliaries
      PRIVATE
        hdf5_util/mh5_nohdf5.F90
    )
  endif()
  find_library (scm_base libscm_base.a)
  find_library (scm_stubs libscm_stubs.a)
  find_library (scm_ftl  libftl.a)
  find_library (scm_z    libz.a)

  if(scm_base MATCHES scm_base-NOTFOUND OR scm_stubs MATCHES scm_stubs-NOTFOUND OR scm_ftl MATCHES scm_ftl-NOTFOUND OR scm_z MATCHES scm_z-NOTFOUND)
    set(ROSE_AMS_LIBRARIES rose_ams_placeholders)
  else()
    set(ROSE_AMS_LIBRARIES ${scm_base} ${scm_stubs} ${scm_ftl} ${scm_z})
  endif()
  target_link_libraries(
    rose_auxiliaries
    PUBLIC
      ${ROSE_LINALG_LIBRARIES}
  )
  set(ROSE_HDF5_BACKEND hdf5_interface_standalone.F90)
  set(ROSE_CORRELATING_BASIS_BACKEND correlating_basisfunctions_standalone.F90)
else()
  add_dependencies(rose_auxiliaries abacus)
  set(ROSE_AMS_LIBRARIES rose_ams_placeholders)
  set(ROSE_HDF5_BACKEND hdf5_interface_dirac.F90)
  set(ROSE_CORRELATING_BASIS_BACKEND correlating_basisfunctions_dirac.F90)
endif()

# ROSE specific subroutines and modules
add_library(
  rose
    # data types and methods to deal with AO and MO bases
    rose_ao.F90
    rose_mo.F90
    rose_structure.F90
    # files containing the implementation of intrinsic fragment orbitals and localization
    generate_ibo.F90
    rose_input_data.F90
    import_fragments.F90
    get_overlap.F90
    change_aobasis.F90
    make_iao.F90
    make_ibo.F90
    avas.F90
    properties.F90
    utils.F90
    check.F90
    # implementation of the addition of correlating functions option
    correlating_basisfunctions.F90
    ${ROSE_CORRELATING_BASIS_BACKEND}
    project_corrfuns.F90
    # interfaces to MO-containing files
    hdf5_interface.F90
    ${ROSE_HDF5_BACKEND}
    ADF_interface.F90
    fchk_interface.F90
    read_corrfuns.F90
    # Placeholders for AMS routines in case we only can link to CoreKF
    rose_noAMS.f90
    rose_main.F90
  )

target_link_libraries(
  rose
    rose_auxiliaries
    ${ROSE_AMS_LIBRARIES}
  )

add_dependencies(rose rose_interest)
add_dependencies(rose rose_auxiliaries)

if(CMAKE_PROJECT_NAME STREQUAL "rose")
  target_sources(
    rose
    PRIVATE
      input_parser.F90
      input_reader.F90
      input_reader_sections.F90
  )

  add_executable(
    rose.x
      rose_standalone.F90
  )

  target_link_libraries(
    rose.x
      rose
      rose_interest
      rose_auxiliaries
      ${ROSE_AMS_LIBRARIES}
  )
  if(ENABLE_HDF5)
    target_link_libraries(
      rose.x
        HDF5::HDF5
    )
  endif()

  install(TARGETS rose.x DESTINATION ${PROJECT_SOURCE_DIR}/bin)
endif()
