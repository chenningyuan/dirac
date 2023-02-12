include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/cfun
    )

include(autocmake_CMakeDependentOption)

# FIXME re-introduce
# set_source_files_properties(${PROJECT_BINARY_DIR}/config_info.F90 PROPERTIES GENERATED 1)


#Miro: remove problematic linking library
list(REMOVE_ITEM CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES  "pthread")

# Ivan Hrasko's determination of compiler versions - very good tool !
include(FindCompilerVersion)

# Miro Ilias restored CMake printout, associated with build_info.h generation
include(build_info)

add_executable(dirac.x)

target_sources(dirac.x PRIVATE src/main/main.F90)

target_link_libraries(dirac.x
  PRIVATE
    abacus
    amfi
    cfun
    dirac
    dirrci
    ecp
    embedding
    esr
    gp
    grid
    input
    interest
    interest
    interface_mpi
    krmc
    laplace
    localization
    luciarel
    lucita
    main
    matrix_operations
    moltra
    prp
    qcorr
    reladc
    relccsd
    visual
    x2c
    xcfun_fortran_bindings
    xcint
  )

if(ENABLE_PCMSOLVER)
  target_link_libraries(dirac.x
    PRIVATE
      pcm_dirac
    )
endif()
if(ENABLE_SRDFT)
  target_link_libraries(dirac.x
    PRIVATE
      srdft
    )
endif()
if(ENABLE_KRCC)
  target_link_libraries(dirac.x
    PRIVATE
      krcc
    )
endif()
target_link_libraries(dirac.x
  PRIVATE
    ${EXTERNAL_LIBS}
  )

#synchronize basis set directories of trunk/basis* into the installation directory
add_custom_command(TARGET dirac.x PRE_LINK COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/basis ${PROJECT_BINARY_DIR}/basis COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/basis_ecp ${PROJECT_BINARY_DIR}/basis_ecp  COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/basis_dalton ${PROJECT_BINARY_DIR}/basis_dalton COMMAND ${CMAKE_COMMAND} -E echo "INFO:basis set directories, basis*, synchronized into current installation directory")

find_package(HDF5 COMPONENTS C HL)

if(ENABLE_PELIB)
    add_dependencies(pelib gen1int_interface)
    add_dependencies(dirac pelib)
endif()

#printout of all linked libraries to dirac.x executable
get_target_property(OUT dirac.x  LINK_LIBRARIES)
message(STATUS "For checking, linked libraries to dirac.x: ${OUT} ")
