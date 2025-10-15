add_library(localization)

target_sources(
  localization
  PRIVATE
    func_PipekMezey.F90
    localization_interface.F90
    second_order_minimization.f90
  )

target_include_directories(
  localization
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(
  localization
  PRIVATE
    gp
  )
