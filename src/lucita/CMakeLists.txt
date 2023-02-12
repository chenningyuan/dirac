add_library(lucita)

target_sources(lucita
  PRIVATE
    doubledir.F
    fock.F
    io.F
    lucia_cc.F
    lucia_icci.F
    property.F
    transform.F
    vecmat.F
    lucita_diagonalize.F
    density.F
    sigma.F
    diagonal.F
    strings.F
    ciinfo.F
    interface.F
    memory.F
    program.F
    utils.F
    par_lucia.F #Miro: Subscript #1 of the array WORK has value 5905696052053 which is greater than the upper bound of 1
    par_solver.F
    flush_stdout.c
  )

target_include_directories(lucita
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(lucita
  PRIVATE
    gp
    qcorr
  )
