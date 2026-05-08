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
    # module to generate coefficients from cartesian to spherical transformation (taken from DIRAC)
    cartesian_to_spherical.F90
    # Placeholders for routines of CoreKF in case we can not link to AMS at all
    ams_placeholders.F90
  )

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
    project_corrfuns.F90
    # interfaces to MO-containing files
    hdf5_interface.F90
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
    ${AMS}
  )

add_dependencies(rose rose_interest)
add_dependencies(rose rose_auxiliaries)
