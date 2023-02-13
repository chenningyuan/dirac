add_library(pelib_interface)

target_sources(pelib_interface
  PRIVATE
    pe_dirac_interfaces.F
  )

target_include_directories(pelib_interface
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

#-------------------------------------------------------------------------------

add_subdirectory(pelib)

target_link_libraries(pelib_interface
  PRIVATE
    pelib
  )

#-------------------------------------------------------------------------------

dirac_test(pe_energy "medium;scf;pelib" "")
dirac_test(pe_response "medium;polarizabilities;pelib" "")
dirac_test(pe_exc "medium;pelib" "")
dirac_test(pe_cpp "short;pelib" "")
