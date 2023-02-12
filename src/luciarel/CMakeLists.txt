add_library(luciarel)

target_sources(luciarel
  PRIVATE
    io_r.F
    par_lushmem_r.F
    par_solver_r.F
    utils_r.F
    vecmat_r.F
    par_lucia_r.F
    property_r.F
    diagonalize_r.F
    diagonal_r.F
    gasci_main_driver_r.F
    memory_r.F
    strings_r.F
    ciinfo_r.F
    sigden_r.F
    interface_r.F
    resolution_utils_r.F
    io_ttss_r.F90
    krci_cfg.F90
    krci_prop.F90
    mospinor_info.F90
    string_info.F90
    symmetry_setup.F90
  )

target_include_directories(luciarel
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(luciarel
  PRIVATE
    gp
  )
