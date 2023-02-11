add_executable(vibcal.x)

target_sources(vibcal.x
  PRIVATE
    vibcal.F90
  )

target_include_directories(vibcal.x
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(vibcal.x
  PRIVATE
    gp
  )

#-------------------------------------------------------------------------------

add_executable(dirac_mointegral_export.x)

target_sources(dirac_mointegral_export.x
  PRIVATE
    dirac_mointegral_export.F90
  )

#-------------------------------------------------------------------------------

add_executable(rspread.x)

target_sources(rspread.x
  PRIVATE
    rspread.F90
  )

target_include_directories(rspread.x
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

#-------------------------------------------------------------------------------

add_executable(labread.x)

target_sources(labread.x
  PRIVATE
    labread.F90
  )

target_include_directories(labread.x
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

#-------------------------------------------------------------------------------

add_executable(twofit.x)

target_sources(twofit.x
  PRIVATE
    twofit.F90
  )

target_include_directories(twofit.x
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(twofit.x
  PRIVATE
    gp
  )

#-------------------------------------------------------------------------------

add_executable(cfread.x)

target_sources(cfread.x
  PRIVATE
    cfread.F90
  )

target_include_directories(cfread.x
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(cfread.x
  PRIVATE
    gp
  )

#-------------------------------------------------------------------------------

add_executable(cf_addlabels.x)

target_sources(cf_addlabels.x
  PRIVATE
    cf_addlabels.F90
  )

target_link_libraries(cf_addlabels.x
  PRIVATE
    gp
  )

#-------------------------------------------------------------------------------

add_executable(pcmo_addlabels.x)

target_sources(pcmo_addlabels.x
  PRIVATE
    pcmo_addlabels.F90
  )

target_link_libraries(pcmo_addlabels.x
  PRIVATE
    gp
  )

#-------------------------------------------------------------------------------

add_executable(polfit.x)

target_sources(polfit.x
  PRIVATE
    polfit.F90
  )

target_include_directories(polfit.x
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(polfit.x
  PRIVATE
    gp
  )

#-------------------------------------------------------------------------------

add_executable(mx2fit.x)

target_sources(mx2fit.x
  PRIVATE
    mx2fit.F90
  )

target_include_directories(mx2fit.x
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(mx2fit.x
  PRIVATE
    gp
  )

#-------------------------------------------------------------------------------

# FIXME
# this binary has so many dependencies
# I will have to re-activate it later
# after disentangling the rest

# add_executable(diag.x)

# target_sources(diag.x
#   PRIVATE
#     diag.F90
#   )

# target_include_directories(diag.x
#   PRIVATE
#     ${PROJECT_SOURCE_DIR}/src/include
#   )

# target_link_libraries(diag.x
#   PRIVATE
#     gp
#     input
#   )

#-------------------------------------------------------------------------------

add_executable(test_allocator.x)

target_sources(test_allocator.x
  PRIVATE
    test_allocator.F90
  )

target_link_libraries(test_allocator.x
  PRIVATE
    gp
  )

#-------------------------------------------------------------------------------

# FIXME move this somewhere else
find_package(HDF5 COMPONENTS C HL)

if(HDF5_FOUND)
  add_executable(h5towfx.x)

  target_sources(h5towfx.x
    PRIVATE
      h5towfx.cpp
    )

  target_link_libraries(h5towfx.x
    PRIVATE
      HDF5::HDF5
    )
endif()

#-------------------------------------------------------------------------------

install(
  TARGETS
    vibcal.x
    dirac_mointegral_export.x
    rspread.x
    labread.x
    twofit.x
    cfread.x
    cf_addlabels.x
    pcmo_addlabels.x
    polfit.x
    mx2fit.x
#   FIXME
#   diag.x
    test_allocator.x
    h5towfx.x
  DESTINATION
    share/dirac
  PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
  )
