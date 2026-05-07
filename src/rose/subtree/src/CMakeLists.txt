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
    # module to read xyz input and store this as a linked list
    read_xyz.F90
    # modules to deal with hdf5(-like) files
    labeled_storage.F90
    hdf5_util/mh5.F90
    hdf5_util/mh5_capi.c
  )

add_library(
  ams_placeholders
    # Placeholders for routines of CoreKF in case we can not link to AMS at all
    ams_placeholders.F90
  )

target_link_libraries(
  rose_auxiliaries
    ${LAPACK_LIBRARIES}
  )

# ROSE specific subroutines and modules
add_library(
  rose
    # data types and methods to deal with AO and MO bases
    rose_ao.F90
    rose_mo.F90
    # modules to read DIRAC keyword-style input (need to be here rather than in auxiliaries, unfortunately)
    input_parser.F90
    input_reader.F90
    # files containing the implementation of intrinsic fragment orbitals and localization
    generate_ibo.F90
    rose_input_data.F90
    input_reader_sections.F90
    import_fragments.F90
    get_overlap.F90
    make_iao.F90
    make_ibo.F90
    properties.F90
    utils.F90
    check.F90
    # interfaces to MO-containing files
    hdf5_interface.F90
    ADF_interface.F90
    fchk_interface.F90
    # Placeholders for AMS routines in case we only can link to CoreKF
    rose_noAMS.f90
  )

find_library (scm_base libscm_base.a)
find_library (scm_stubs libscm_stubs.a)
find_library (scm_ftl  libftl.a)
find_library (scm_z    libz.a)

if(scm_base MATCHES scm_base-NOTFOUND OR scm_stubs MATCHES scm_stubs-NOTFOUND OR scm_ftl MATCHES scm_ftl-NOTFOUND OR scm_z MATCHES scm_z-NOTFOUND)
  set(AMS ams_placeholders)
else()
  set(AMS ${scm_base} ${scm_stubs} ${scm_ftl} ${scm_z})
endif()

target_link_libraries(
  rose
    rose_auxiliaries
    ${AMS}
  )

add_executable(
  project_corrfuns.x
    project_corrfuns.F90
    read_corrfuns.F90
  )

target_link_libraries(
  project_corrfuns.x
    rose
    rose_interest
    rose_auxiliaries
    HDF5::HDF5
    ${AMS}
  )


add_executable(
  rose.x
    rose_main.F90
  )

target_link_libraries(
  rose.x
    rose
    rose_interest
    rose_auxiliaries
    HDF5::HDF5
    ${AMS}
  )

add_dependencies(rose rose_interest)
add_dependencies(rose rose_auxiliaries)

install(TARGETS rose.x project_corrfuns.x DESTINATION ${PROJECT_SOURCE_DIR}/bin)
