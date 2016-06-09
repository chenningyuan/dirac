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
option(ENABLE_RELCCSD_STANDALONE "Enable standalone RELCCSD" ON)
option(ENABLE_XCFUN "Enable XCFun library" ON)
if(NOT ENABLE_XCFUN)
    add_definitions(-DEXCLUDE_XCFUN)
endif()
option(ENABLE_QCORR "Davidson-type +Q corrections for (MR)CISD" ON)
if(ENABLE_QCORR)
    add_definitions(-DMOD_QCORR)
endif()
option(ENABLE_STIELTJES "Enable stieltjes library" ON)
if(ENABLE_STIELTJES)
    add_definitions(-DHAS_STIELTJES)
endif()
# interest as separate library
cmake_dependent_option(ENABLE_INTEREST "Enable Interest" ON  "ENABLE_UNRELEASED" OFF)
if(ENABLE_INTEREST)
    add_definitions(-DMOD_INTEREST)
endif()
option(ENABLE_OPENRSP "Enable OpenRSP" OFF)
if(ENABLE_OPENRSP)
    add_definitions(-DMOD_OPENRSP)
endif()
cmake_dependent_option(ENABLE_LAO_REARRANGED "LAO properties without connection matrices" ON "ENABLE_UNRELEASED" OFF)
if(ENABLE_LAO_REARRANGED)
    add_definitions(-DMOD_LAO_REARRANGED)
endif()
cmake_dependent_option(ENABLE_MCSCF_spinfree "Spinfree MCSCF module" ON "ENABLE_UNRELEASED" OFF)
if(ENABLE_MCSCF_spinfree)
    add_definitions(-DMOD_MCSCF_spinfree)
endif()
cmake_dependent_option(ENABLE_XAMFI "Atomic 2-electron scalar and spin-orbit correction module" ON "ENABLE_UNRELEASED" OFF)
if(ENABLE_XAMFI)
    add_definitions(-DMOD_XAMFI)
endif()
cmake_dependent_option(ENABLE_ERI "ERI integral routines" ON "ENABLE_UNRELEASED" OFF)
if(ENABLE_ERI)
    add_definitions(-DMOD_ERI)
endif()
cmake_dependent_option(ENABLE_DENSFIT "Density fitting module" ON "ENABLE_UNRELEASED" OFF)
if(ENABLE_DENSFIT)
    add_definitions(-DMOD_DNF)
endif()
cmake_dependent_option(ENABLE_ESR "ESR property module" ON "ENABLE_UNRELEASED" OFF)
if(ENABLE_ESR)
    add_definitions(-DMOD_ESR)
endif()
cmake_dependent_option(ENABLE_KRCC "KRCC module" ON "ENABLE_UNRELEASED" OFF)
if(ENABLE_KRCC)
    add_definitions(-DMOD_KRCC)
endif()
cmake_dependent_option(ENABLE_SRDFT "srDFT module" ON  "ENABLE_UNRELEASED" OFF)
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
# of DIRAC (main,dirac,relccsd,prp,reladc, interest ...)
#
# Due to entanglements between libraries (like dirac depends localization and vise versa) it is necessary
# to keep source files (modules) as localized as possible and with minimal interdependencie on other modules.
#
#----------------------------------------------------------------------------------------------------------

# make the EXTERNAL_LIBS variable visible to the scope up !
set(EXTERNAL_LIBS PARENT_SCOPE)

# pdpack as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/pdpack)
set(EXTERNAL_LIBS pdpack ${EXTERNAL_LIBS})

# Making interface_mpi as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/interface_mpi)
set(EXTERNAL_LIBS interface_mpi ${EXTERNAL_LIBS})

# gp as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/gp)
set(EXTERNAL_LIBS gp ${EXTERNAL_LIBS})

# Abacus as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/abacus)
set(EXTERNAL_LIBS abacus ${EXTERNAL_LIBS})

# Making matrix_operations  as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/matrix_operations)
set(EXTERNAL_LIBS matrix_operations ${EXTERNAL_LIBS})

#    Moltra as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/moltra)
set(EXTERNAL_LIBS moltra ${EXTERNAL_LIBS})

#    embedding as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/embedding)
set(EXTERNAL_LIBS embedding ${EXTERNAL_LIBS})

# Making x2c as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/x2c)
set(EXTERNAL_LIBS x2c ${EXTERNAL_LIBS})

# Making localization  as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/localization)
set(EXTERNAL_LIBS localization ${EXTERNAL_LIBS})

# cfun as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/cfun)
set(EXTERNAL_LIBS cfun ${EXTERNAL_LIBS})

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
set(EXTERNAL_LIBS xcint ${EXTERNAL_LIBS})

if(ENABLE_QCORR)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/qcorr)
    set(EXTERNAL_LIBS qcorr ${EXTERNAL_LIBS})
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_QCORR)
endif()
message(STATUS "Davidson-type +Q corrections for (MR)CISD: ${ENABLE_QCORR}")

# luciarel as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/luciarel)
set(EXTERNAL_LIBS luciarel ${EXTERNAL_LIBS})

# lucita as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/lucita)
set(EXTERNAL_LIBS lucita ${EXTERNAL_LIBS})

# krmc as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/krmc)
set(EXTERNAL_LIBS krmc ${EXTERNAL_LIBS})

#     esr as separate library
if(ENABLE_ESR)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/prp/esr)
    set(EXTERNAL_LIBS esr ${EXTERNAL_LIBS})
    dirac_test(esr "short;gascip;ci")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_ESR)
endif()
message(STATUS "ESR property module: ${ENABLE_ESR}")

# ecp as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/ecp)
set(EXTERNAL_LIBS ecp ${EXTERNAL_LIBS})

# amfi as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/amfi)
set(EXTERNAL_LIBS amfi ${EXTERNAL_LIBS})

# dirrci as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/dirrci)
set(EXTERNAL_LIBS dirrci ${EXTERNAL_LIBS})

# Making reladc as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/reladc)
set(EXTERNAL_LIBS reladc ${EXTERNAL_LIBS})

# Making visual as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/visual)
set(EXTERNAL_LIBS visual ${EXTERNAL_LIBS})

# grid as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/grid)
set(EXTERNAL_LIBS grid ${EXTERNAL_LIBS})

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
  dirac_test(pcm_energy "short;scf;pcm")
  dirac_test(pcm_response "short;polarizabilities;pcm")

  set(EXTERNAL_LIBS pcm_dirac ${EXTERNAL_LIBS})
endif()

# prp as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/prp)
set(EXTERNAL_LIBS prp ${EXTERNAL_LIBS})

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
  set(EXTERNAL_LIBS
    ${PROJECT_BINARY_DIR}/external/lib/libstieltjes.a
    ${EXTERNAL_LIBS})
else()
  message(STATUS "Stieltjes external module DISABLED")
endif()

# interest as separate library
if(ENABLE_INTEREST)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/interest)
    set(EXTERNAL_LIBS interest ${EXTERNAL_LIBS})
    dirac_test(mdirac_closed-shell_dft "interest")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_INTEREST)
endif()
message(STATUS "Interest library: ${ENABLE_INTEREST}")

#   krcc as separate library
if(ENABLE_KRCC)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/krcc)
    set(EXTERNAL_LIBS krcc ${EXTERNAL_LIBS})
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_KRCC)
endif()
message(STATUS "KRCC module: ${ENABLE_KRCC}")

# Making relccsd as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/relccsd)
set(EXTERNAL_LIBS relccsd ${EXTERNAL_LIBS})
#  relccsd.x as standalone executable by default
message(STATUS "Enable compilation of standalone relccsd.x: ${ENABLE_RELCCSD_STANDALONE}")

#  openrsp as separate library
if(ENABLE_OPENRSP)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/openrsp)
    set(EXTERNAL_LIBS openrsp ${EXTERNAL_LIBS})
    dirac_test(legacy_openrsp_cars "openrsp")
    dirac_test(legacy_openrsp_cme "openrsp")
    dirac_test(legacy_openrsp_delta "openrsp")
    dirac_test(legacy_openrsp_efgb "openrsp")
    dirac_test(legacy_openrsp_efgb_dft "openrsp")
    dirac_test(legacy_openrsp_electric "openrsp")
    dirac_test(legacy_openrsp_jones "openrsp")
    dirac_test(legacy_openrsp_magnetizability "openrsp")
    dirac_test(legacy_openrsp_molecular_gradient "openrsp")
    dirac_test(legacy_openrsp_pv "openrsp")
    dirac_test(legacy_openrsp_quadrupole "openrsp")
    dirac_test(legacy_openrsp_roa "openrsp")
    dirac_test(legacy_openrsp_pnc_gradient "openrsp")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_OPENRSP)
endif()
message(STATUS "OpenRSP library: ${ENABLE_OPENRSP}")

#    lao_rearranged as separate library
if(ENABLE_LAO_REARRANGED)
    dirac_test(lao_shield_rearranged "")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_LAO_REARRANGED)
endif()
message(STATUS "LAO properties without connection matrices: ${ENABLE_LAO_REARRANGED}")

#  mcscf_sf as separate library
if(NOT ENABLE_MCSCF_spinfree)
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_MCSCF_spinfree)
endif()
message(STATUS "Spinfree MCSCF module: ${ENABLE_MCSCF_spinfree}")

# XAMFI as separate library
if(ENABLE_XAMFI)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/x-amfi)
    set(EXTERNAL_LIBS x-amfi ${EXTERNAL_LIBS})
    dirac_test(atomic_x2c_2esoc_models_energy "medium,xamfi")
    dirac_test(atomic_x2c_spinfree_energy "medium")
    dirac_test(density_at_nuclei_2c_4c "medium,xamfi")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_XAMFI)
endif()
message(STATUS "Atomic oo-order spin-orbit correction module: ${ENABLE_XAMFI}")

#  eri as separate library
if(ENABLE_ERI)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/eri)
    set(EXTERNAL_LIBS eri ${EXTERNAL_LIBS})
    dirac_test(integrals_eri "")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_ERI)
endif()
message(STATUS "ERI integral routines: ${ENABLE_ERI}")

#  densfit as separate library
if(ENABLE_DENSFIT)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/densfit)
    set(EXTERNAL_LIBS densfit ${EXTERNAL_LIBS})
    dirac_test(integrals_densfit "")
    dirac_test(integrals_jengine "")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_DNF)
endif()
message(STATUS "Density fitting module: ${ENABLE_DENSFIT}")

#  srdft module as separate library
if(ENABLE_SRDFT)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/srdft)
    set(EXTERNAL_LIBS srdft ${EXTERNAL_LIBS})
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_SRDFT)
endif()
message(STATUS "srDFT module: ${ENABLE_SRDFT}")

#  dirac sources as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/dirac)
set(EXTERNAL_LIBS dirac ${EXTERNAL_LIBS})

#  input reader directory as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/input)
set(EXTERNAL_LIBS input ${EXTERNAL_LIBS})

# main directory as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/main)
set(EXTERNAL_LIBS main ${EXTERNAL_LIBS})

# make the EXTERNAL_LIBS variable visible to the scope up !
set(EXTERNAL_LIBS ${EXTERNAL_LIBS} PARENT_SCOPE)
