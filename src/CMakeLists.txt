option(ENABLE_GEN1INT "Enable Gen1Int library" ON)

if(ENABLE_GEN1INT)
  add_definitions(-DBUILD_GEN1INT)
  add_subdirectory(gen1int)
endif()

#-------------------------------------------------------------------------------

option(ENABLE_OPENRSP "Enable OpenRSP" OFF)

if(ENABLE_OPENRSP)
  add_definitions(-DMOD_OPENRSP)
  add_subdirectory(openrsp)
endif()

#-------------------------------------------------------------------------------

option(ENABLE_KRCC "KRCC module" ON)

if(ENABLE_KRCC)
  add_definitions(-DMOD_KRCC)
  add_subdirectory(krcc)
endif()

#-------------------------------------------------------------------------------

option(ENABLE_XCFUN "Enable XCFun library" ON)

if(ENABLE_XCFUN)
  add_subdirectory(xcfun)
else()
  add_definitions(-DEXCLUDE_XCFUN)
endif()

#-------------------------------------------------------------------------------

option(ENABLE_PELIB "Enable polarizabile embedding library" ON)

if(ENABLE_PELIB)
  add_definitions(-DHAS_PELIB)
  add_subdirectory(${PROJECT_SOURCE_DIR}/src/pelib)

  if(NOT ENABLE_GEN1INT)
    message(FATAL_ERROR "-- PElib requires Gen1Int, use -DENABLE_GEN1INT=ON to enable Gen1Int or -DENABLE_PELIB=OFF to disable PElib")
  endif()
endif()

#-------------------------------------------------------------------------------

option(ENABLE_MCSCF_spinfree "Spinfree MCSCF module" ON)

if(ENABLE_MCSCF_spinfree)
  add_definitions(-DMOD_MCSCF_spinfree)
endif()

#-------------------------------------------------------------------------------

# FIXME
option(ENABLE_XAMFI "(e)amf 2-electron scalar and spin-orbit picture-change correction module" ON)

if(ENABLE_XAMFI AND HDF5_FOUND)
  add_definitions(-DMOD_XAMFI)
else()
  if(ENABLE_XAMFI AND NOT HDF5_FOUND)
      message(STATUS "X-AMFI, the (e)amf 2-electron scalar and spin-orbit picture-change correction module requires a working HDF5 library")
  endif()
  set(ENABLE_XAMFI OFF)
endif()

#-------------------------------------------------------------------------------

add_subdirectory(pdpack)
add_subdirectory(interface_mpi)
add_subdirectory(gp)
add_subdirectory(abacus)
add_subdirectory(matrix_operations)

#-------------------------------------------------------------------------------

option(ENABLE_EXATENSOR "Enable ExaTENSOR library" ON)

if(ENABLE_64BIT_INTEGERS)
  message(STATUS "ExaTensor library can not work with 64 bit integers, switched off")
  set(ENABLE_EXATENSOR OFF)
endif()

if(ENABLE_EXATENSOR)
  add_subdirectory(exacorr)
  add_definitions(-DMOD_EXACORR)
endif()

#-------------------------------------------------------------------------------

add_subdirectory(moltra)
add_subdirectory(embedding)
add_subdirectory(x2c)
add_subdirectory(localization)
add_subdirectory(cfun)
add_subdirectory(xcint)
add_subdirectory(qcorr)
add_subdirectory(luciarel)
add_subdirectory(lucita)
add_subdirectory(krmc)
add_subdirectory(ecp)
add_subdirectory(amfi)
add_subdirectory(dirrci)
add_subdirectory(reladc)
add_subdirectory(visual)
add_subdirectory(grid)

#-------------------------------------------------------------------------------

# pcmsolver as separate external library
if(ENABLE_PCMSOLVER)
  include_directories(${PROJECT_BINARY_DIR}/external/pcmsolver/src/pcmsolver-build/modules)
  #    set related pcm_dirac as separate library
  # ... needs HAS_PCMSOLVER definition !
  add_subdirectory(${PROJECT_SOURCE_DIR}/src/pcm)
  set(EXTERNAL_LIBS
    ${PCMSOLVER_LIBS}
    ${EXTERNAL_LIBS}
    )
  if (ENABLE_STATIC_LINKING)
    list(REMOVE_ITEM EXTERNAL_LIBS "c")
    message(STATUS "PCM accompanying libs cleared of 'libc' due to static linking")
  endif()
  dirac_test(pcm_energy "short;scf;pcm" "")
  dirac_test(pcm_response "medium;polarizabilities;pcm" "")

  set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:pcm_dirac> ${EXTERNAL_OBJECTS})
endif()

#-------------------------------------------------------------------------------

add_subdirectory(prp)
add_subdirectory(interest)
add_subdirectory(relccsd)

#-------------------------------------------------------------------------------

option(ENABLE_SRDFT "srDFT module" ON)

if(ENABLE_SRDFT)
  add_subdirectory(srdft)
endif()

#-------------------------------------------------------------------------------

add_subdirectory(laplace-minimax)
add_subdirectory(laplace)
add_subdirectory(dirac)
add_subdirectory(input)
add_subdirectory(main)

#-------------------------------------------------------------------------------

add_executable(dirac.x)

target_sources(dirac.x
  PRIVATE
    main/main.F90
  )

target_include_directories(dirac.x
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

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
    pelib_interface
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

#-------------------------------------------------------------------------------

# copy basis set libraries to build directory
add_custom_command(
  TARGET
    dirac.x
  COMMAND
    ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/basis ${PROJECT_BINARY_DIR}/basis
  COMMAND
    ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/basis_ecp ${PROJECT_BINARY_DIR}/basis_ecp
  COMMAND
    ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/basis_dalton ${PROJECT_BINARY_DIR}/basis_dalton
  )
