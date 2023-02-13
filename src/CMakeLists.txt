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
# of DIRAC (main,dirac,relccsd,prp,reladc, interest ...)
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

# gp as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/gp)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:gp> ${EXTERNAL_OBJECTS})

# Abacus as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/abacus)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:abacus> ${EXTERNAL_OBJECTS})

# Making matrix_operations  as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/matrix_operations)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:matrix_operations> ${EXTERNAL_OBJECTS})

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

#    Moltra as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/moltra)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:moltra> ${EXTERNAL_OBJECTS})

#    embedding as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/embedding)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:embedding> ${EXTERNAL_OBJECTS})

# Making x2c as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/x2c)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:x2c> ${EXTERNAL_OBJECTS})

# Making localization  as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/localization)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:localization> ${EXTERNAL_OBJECTS})

# cfun as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/cfun)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:cfun> ${EXTERNAL_OBJECTS})

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

# interest as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/interest)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:interest> ${EXTERNAL_OBJECTS})


# Making relccsd as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/relccsd)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:relccsd> ${EXTERNAL_OBJECTS})

#  mcscf_sf as separate library
message(STATUS "Spinfree MCSCF module: ${ENABLE_MCSCF_spinfree}")

# XAMFI as separate library
if(ENABLE_XAMFI)
    dirac_test(atomic_x2c_2esoc_models_energy "medium;xamfi" "")
    dirac_test(atomic_x2c_spinfree_energy "medium;xamfi" "")
    dirac_test(x_amfi "short;xamfi" "")
endif()

#  srdft module as separate library
if(ENABLE_SRDFT)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/srdft)
    set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:srdft> ${EXTERNAL_OBJECTS})
endif()
message(STATUS "srDFT module: ${ENABLE_SRDFT}")

# laplace-minimax
add_subdirectory(${PROJECT_SOURCE_DIR}/src/laplace-minimax)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/laplace)

#  dirac sources as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/dirac)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:dirac> ${EXTERNAL_OBJECTS})

#  input reader directory as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/input)
set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:input> ${EXTERNAL_OBJECTS})

# main directory as separate library
add_subdirectory(${PROJECT_SOURCE_DIR}/src/main)
set(EXTERNAL_OBJECTS  $<TARGET_OBJECTS:main> ${EXTERNAL_OBJECTS})

# make these variables visible to the cmake one level up
set(EXTERNAL_LIBS ${EXTERNAL_LIBS} PARENT_SCOPE)
set(EXTERNAL_OBJECTS ${EXTERNAL_OBJECTS} PARENT_SCOPE)
