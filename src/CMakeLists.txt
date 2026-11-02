#################################################################################
#
# Big substitution of many autocmake.cfg's custom/FUNCTIONALITY.cmake files.
#
# Written by Miro Ilias, February 2016, GSI
#
#################################################################################

#----------------------------------------------------------------------------------
# Set all module options and definitions first
# so they can be DEFINED for all source subdirectories below
#----------------------------------------------------------------------------------



option(ENABLE_GEN1INT "Enable Gen1Int library" ON)
if(ENABLE_GEN1INT)
  # fixme: later this needs to become a target_compile_definitions
  add_definitions(-DBUILD_GEN1INT)

  add_subdirectory(${PROJECT_SOURCE_DIR}/src/gen1int)

  set(EXTERNAL_LIBS gen1int_interface ${EXTERNAL_LIBS})
endif()



option(ENABLE_PELIB "Enable polarizabile embedding library" ON)
option(ENABLE_PELIB_PDE "Enable polarizabile density embedding library" ON)
if(ENABLE_PELIB)
  # fixme: later this needs to become a target_compile_definitions
  add_definitions(-DHAS_PELIB)

  add_subdirectory(${PROJECT_SOURCE_DIR}/src/pelib)

  # fixme: this will later be done using targets
  add_dependencies(PElib gen1int_interface)
  add_dependencies(pelib_interface interface_mpi)

  set(EXTERNAL_LIBS pelib_interface ${EXTERNAL_LIBS})
endif()



option(ENABLE_RELCCSD_STANDALONE "Enable standalone RELCCSD" OFF)
option(ENABLE_XCFUN "Enable XCFun library" ON)
if(NOT ENABLE_XCFUN)
    add_definitions(-DEXCLUDE_XCFUN)
endif()

option(ENABLE_STIELTJES "Enable stieltjes library" OFF)
if(ENABLE_STIELTJES)
    add_definitions(-DHAS_STIELTJES)
endif()
option(ENABLE_OPENRSP "Enable OpenRSP" OFF)
if(ENABLE_OPENRSP)
    add_definitions(-DMOD_OPENRSP)
endif()
option(ENABLE_LAO_REARRANGED "LAO properties without connection matrices" ON)
if(ENABLE_LAO_REARRANGED)
    add_definitions(-DMOD_LAO_REARRANGED)
endif()
option(ENABLE_MCSCF_spinfree "Spinfree MCSCF module" ON)
if(ENABLE_MCSCF_spinfree)
    add_definitions(-DMOD_MCSCF_spinfree)
endif()
option(ENABLE_XAMFI "(e)amf 2-electron scalar and spin-orbit picture-change correction module" ON)
if(ENABLE_XAMFI AND HDF5_FOUND)
    add_definitions(-DMOD_XAMFI)
else()
    if(ENABLE_XAMFI AND NOT HDF5_FOUND)
        message(STATUS "X-AMFI, the (e)amf 2-electron scalar and spin-orbit picture-change correction module requires a working HDF5 library")
    endif()
    set(ENABLE_XAMFI OFF)
endif()
option(ENABLE_ESR "ESR property module" ON)
if(ENABLE_ESR)
    add_definitions(-DMOD_ESR)
endif()
option(ENABLE_KRCC "KRCC module" ON)
if(ENABLE_KRCC)
    add_definitions(-DMOD_KRCC)
endif()
option(ENABLE_SRDFT "srDFT module" ON)
if(ENABLE_SRDFT)
    add_definitions(-DMOD_SRDFT)
endif()

#--------------------------------------------------------------------------------------------------------------
#
# WARNING: The order of the functionality modules is crucial for the proper linking of all DIRAC executables !
#
# The first module here (pdpack) is the last in the linking command due way of filling the ${EXTERNAL_LIBS}
# variable.
#
# The resulting ${EXTERNAL_LIBS} variable contains localized source libraries of all functionality modules
# of DIRAC (dirac_main,dirac,relccsd,prp,reladc, interest ...)
#
# Due to entanglements between libraries (like dirac depends localization and vise versa) it is necessary
# to keep source files (modules) as localized as possible and with minimal interdependencie on other modules.
#
#----------------------------------------------------------------------------------------------------------

# make the EXTERNAL_LIBS variable visible to the scope up !
set(EXTERNAL_LIBS PARENT_SCOPE)
set(EXTERNAL_OBJECTS PARENT_SCOPE)

# pdpack as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/pdpack)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:pdpack> ${EXTERNAL_OBJECTS})

# Making interface_mpi as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/interface_mpi)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:interface_mpi> ${EXTERNAL_OBJECTS})




# Making mfsolver as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/mfsolver)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:mfsolver> ${EXTERNAL_OBJECTS})

# gp as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/gp)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:gp> ${EXTERNAL_OBJECTS})

# Abacus as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/abacus)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:abacus> ${EXTERNAL_OBJECTS})

# Making matrix_operations  as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/matrix_operations)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:matrix_operations> ${EXTERNAL_OBJECTS})

# laplace-minimax
add_subdirectory(${PROJECT_SOURCE_DIR}/src/laplace-minimax)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/laplace)

# exacorr as separate library
if(ENABLE_EXATENSOR)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/exacorr)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:exacorr> ${EXTERNAL_OBJECTS})
add_definitions(-DMOD_EXACORR)
endif()

#    Moltra as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/moltra)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:moltra> ${EXTERNAL_OBJECTS})

#    frozen density embedding as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/fde)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:fde> ${EXTERNAL_OBJECTS})

# XAMFI as separate library; before x2c because x2c depends on XAMFI
if(ENABLE_XAMFI)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/x-amfi)
    set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:x-amfi> ${EXTERNAL_OBJECTS})
    dirac_test(atomic_x2c_2esoc_models_energy "medium;xamfi" "")
    dirac_test(atomic_x2c_spinfree_energy "medium;xamfi" "")
    dirac_test(x_amfi "medium;xamfi" "")
    message(STATUS "NOTE: tests for the (e)amfX2C 2e-scalar and spin-orbit picture-change correction module require h5py in your python env!")
endif()
message(STATUS "(e)amfX2C 2e-scalar and spin-orbit picture-change correction module: ${ENABLE_XAMFI}")

# Making x2c as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/x2c)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:x2c> ${EXTERNAL_OBJECTS})

# Making localization  as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/localization)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:localization> ${EXTERNAL_OBJECTS})

# cfun as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/cfun)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:cfun> ${EXTERNAL_OBJECTS})

# xcfun as external separate library
if(ENABLE_XCFUN)
    message(STATUS "The XCFun submodule ENABLED")
    include(ExternalProject)
    # miro: pass carefully all necessary arguments !
    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
        -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS} )
    ExternalProject_Add(xcfun
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/xcfun
        BINARY_DIR ${PROJECT_BINARY_DIR}/external/xcfun-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/xcfun-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/xcfun-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${ExternalProjectCMakeArgs}  )
    if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
        SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
        message(STATUS "For libxcfun_fortran_bindings.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
    endif()
    add_library(
        xcfun_fortran_bindings
        ${PROJECT_SOURCE_DIR}/src/xcfun/fortran/xcfun_autogen.F90
        ${PROJECT_SOURCE_DIR}/src/xcfun/fortran/xcfun_module.F90  )
    # Miro: how will this cope with OBJECT libraries ?
    set(EXTERNAL_LIBS
        xcfun_fortran_bindings
        ${PROJECT_BINARY_DIR}/external/xcfun-build/src/libxcfun.a
        ${EXTERNAL_LIBS} )
    add_dependencies(xcfun_fortran_bindings xcfun)
    include_directories(${PROJECT_BINARY_DIR}/external/xcfun-build/fortran)
else()
    message(STATUS "The XCFun submodule DISABLED")
endif()

# xcint as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/xcint)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:xcint> ${EXTERNAL_OBJECTS})

add_subdirectory(${PROJECT_SOURCE_DIR}/src/qcorr)

# luciarel as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/luciarel)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:luciarel> ${EXTERNAL_OBJECTS})

# lucita as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/lucita)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:lucita> ${EXTERNAL_OBJECTS})

# krmc as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/krmc)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:krmc> ${EXTERNAL_OBJECTS})

#     esr as separate library
if(ENABLE_ESR)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/prp/esr)
    set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:esr> ${EXTERNAL_OBJECTS})
    dirac_test(esr "medium;gascip;ci" "")
endif()
message(STATUS "ESR property module: ${ENABLE_ESR}")

# ecp as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/ecp)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:ecp> ${EXTERNAL_OBJECTS})

# amfi as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/amfi)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:amfi> ${EXTERNAL_OBJECTS})

# dirrci as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/dirrci)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:dirrci> ${EXTERNAL_OBJECTS})

# Making reladc as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/reladc)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:reladc> ${EXTERNAL_OBJECTS})

# Making visual as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/visual)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:visual> ${EXTERNAL_OBJECTS})

# grid as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/grid)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:grid> ${EXTERNAL_OBJECTS})

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

# prp as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/prp)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:prp> ${EXTERNAL_OBJECTS})

# rose as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/rose/subtree/src)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:rose> ${EXTERNAL_OBJECTS})
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:rose_auxiliaries> ${EXTERNAL_OBJECTS})
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:rose_ams_placeholders> ${EXTERNAL_OBJECTS})
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:rose_interest> ${EXTERNAL_OBJECTS})

#  stieltjes as separate external library
if(ENABLE_STIELTJES)
  message(STATUS "Stieltjes external module ENABLED")
  include_directories (${PROJECT_BINARY_DIR}/external/stieltjes-build)
  set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}  )
  add_external(stieltjes)
  # Miro: how will this cope with OBJECT libraries ?
  set(EXTERNAL_LIBS
    ${PROJECT_BINARY_DIR}/external/lib/libstieltjes.a
    ${EXTERNAL_LIBS})
else()
  message(STATUS "Stieltjes external module DISABLED")
endif()

# interest as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/interest)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:interest> ${EXTERNAL_OBJECTS})

#   krcc as separate library
if(ENABLE_KRCC)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/krcc)
    set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:krcc> ${EXTERNAL_OBJECTS})
endif()
message(STATUS "KRCC module: ${ENABLE_KRCC}")

# Making caspt2 as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/caspt2)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:caspt2> ${EXTERNAL_OBJECTS})

# Making relccsd as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/relccsd)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:relccsd> ${EXTERNAL_OBJECTS})
#  relccsd.x as standalone executable by default
message(STATUS "Enable compilation of standalone relccsd.x: ${ENABLE_RELCCSD_STANDALONE}")

#  openrsp as separate library
if(ENABLE_OPENRSP)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/openrsp)
    set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:openrsp> ${EXTERNAL_OBJECTS})
    dirac_test(legacy_openrsp_cars "openrsp" "")
    dirac_test(legacy_openrsp_cme "openrsp" "")
    dirac_test(legacy_openrsp_delta "openrsp" "")
    dirac_test(legacy_openrsp_efgb "openrsp" "")
    dirac_test(legacy_openrsp_efgb_dft "openrsp" "")
    dirac_test(legacy_openrsp_electric "openrsp" "")
    dirac_test(legacy_openrsp_jones "openrsp" "")
    dirac_test(legacy_openrsp_magnetizability "openrsp" "")
    dirac_test(legacy_openrsp_molecular_gradient "openrsp" "")
    dirac_test(legacy_openrsp_pv "openrsp" "")
    dirac_test(legacy_openrsp_quadrupole "openrsp" "")
    dirac_test(legacy_openrsp_roa "openrsp" "")
    dirac_test(legacy_openrsp_pnc_gradient "openrsp" "")
endif()
message(STATUS "OpenRSP library: ${ENABLE_OPENRSP}")

#    lao_rearranged as separate library
if(ENABLE_LAO_REARRANGED)
    dirac_test(lao_shield_rearranged "long" "")
endif()
message(STATUS "LAO properties without connection matrices: ${ENABLE_LAO_REARRANGED}")

#  mcscf_sf as separate library
message(STATUS "Spinfree MCSCF module: ${ENABLE_MCSCF_spinfree}")

#  srdft module as separate library
if(ENABLE_SRDFT)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/srdft)
    set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:srdft> ${EXTERNAL_OBJECTS})
endif()
message(STATUS "srDFT module: ${ENABLE_SRDFT}")

#  dirac sources as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/dirac)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:dirac> ${EXTERNAL_OBJECTS})

#  input reader directory as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/input)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:input> ${EXTERNAL_OBJECTS})

# Modules with variables
add_subdirectory(${PROJECT_SOURCE_DIR}/src/modules)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:modules> ${EXTERNAL_OBJECTS})
### Compile before other targets
function(link_modules target_name)
    if (TARGET ${target_name})
        target_link_libraries(${target_name} PRIVATE modules)
        add_dependencies(${target_name} modules)
    endif()
endfunction()
### TARGETS that depend on modules
link_modules(interface_mpi)
link_modules(interest)
link_modules(mfsolver)
link_modules(abacus)
link_modules(gp)
link_modules(matrix_operations)
link_modules(exacorr)
link_modules(moltra)
link_modules(esr)
link_modules(ecp)
link_modules(amfi)
link_modules(dirac)
link_modules(dirrci)
link_modules(fde)
link_modules(grid)
link_modules(input)
link_modules(krcc)
link_modules(krmc)
link_modules(localization)
link_modules(luciarel)
link_modules(lucita)
link_modules(openrsp)
link_modules(pdpack)
link_modules(pelib_interface)
link_modules(prp)
link_modules(qcorr)
link_modules(reladc)
link_modules(relccsd)

# main directory as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/main)
link_modules(dirac_main)
set(EXTERNAL_OBJECTS  $<TARGET_OBJECTS:dirac_main> ${EXTERNAL_OBJECTS})

# make these variables visible to the cmake one level up
set(EXTERNAL_LIBS ${EXTERNAL_LIBS} PARENT_SCOPE)
set(EXTERNAL_OBJECTS ${EXTERNAL_OBJECTS} PARENT_SCOPE)
