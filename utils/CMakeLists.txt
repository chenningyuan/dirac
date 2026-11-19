add_executable(labread.x)

target_sources(
  labread.x
  PRIVATE
    labread.F90
  )

target_include_directories(
  labread.x
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

set_target_properties(
  labread.x
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
  )
