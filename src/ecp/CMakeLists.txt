add_library(ecp)

target_sources(ecp
  PRIVATE
    recp1.F
    recp_lnk.F
    recp_chk.F90
    recp_cso.F90
    recp_cso_ang.F90
    recp_cso_ps.F90
    recp_cso_rad.F90
    recp_cso_sub.F90
    recp_cso_table.F90
    recp_function1.F90
    recp_function2.F90
    recp_inp.F90
    recp_inp_orb.F90
    recp_inp_print.F90
    recp_inp_read.F90
    recp_ipt.F90
    recp_ntr.F90
    recp_one.F90
    recp_one_stvcz.F90
    recp_one_sub.F90
    recp_output.F90
    recp_socfpd.F90
    recp_write.F90
  )

target_include_directories(ecp
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

# target_link_libraries(interface_mpi
#   PRIVATE
#     gp
#     xcint
#     visual
#     embedding
#   )


