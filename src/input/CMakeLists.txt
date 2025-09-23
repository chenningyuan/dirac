add_library(input STATIC)

target_sources(input
  PRIVATE
    input_reader.F90
    input_reader_sections.F90
    parse_input.F90
  )

#target_link_libraries(input
#  PRIVATE
#    reladc
#    gp
#    dirac
#    dirac_main
#    visual
#    relccsd
#    grid
#    qcorr
#    caspt2
#  )

# dependencies, because input needs fortran modules compiled in the listed libraries
add_dependencies(input reladc gp dirac dirac_main visual relccsd grid qcorr caspt2)

if(ENABLE_EXATENSOR)
#  target_link_libraries(input
#    PRIVATE
#      exacorr
#    )
  add_dependencies(input exacorr)
endif()

target_include_directories(input
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )


