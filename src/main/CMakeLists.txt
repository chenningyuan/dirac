add_library(main)

target_sources(main
  PRIVATE
    dirac.F
    logo.F90
    mpi_framework.F90
  )

target_include_directories(main
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(main
  PRIVATE
    gp
    xcint
    visual
    embedding
#   interest
#   laplace-minimax
#   interface_mpi
  )

if(ENABLE_PCMSOLVER)
  target_link_libraries(main
    PRIVATE
      pcm_dirac
    )
endif()

# FIXME
# if(ENABLE_STIELTJES)
#   target_link_libraries(main
#     PRIVATE
#       stieltjes
#     )
# endif()
