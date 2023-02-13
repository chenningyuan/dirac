add_library(dirac)

target_sources(dirac
  PRIVATE
    atom2molecule_utils.F
    dftaux.F
    dirana.F
    diraos.F
    dirbss.F
    dircnv.F
    dirden.F
    dirgeo.F
    dirgp.F
    dirgrd.F
    dirmp2no.F
    dirone.F
    diropen.F
    dirout.F
    dirrdn.F
    dirscf.F
    dirset.F
    dirsol.F
    dirspf.F
    dirtra.F
    dirtwo.F
    dirv1c.F
    huckel.F
    dirimpmos.F90
 )

target_link_libraries(dirac
  PRIVATE
    amfi
    localization
    lucita
    relccsd
    x2c
    srdft
    abacus
    # FIXME
    ${PROJECT_BINARY_DIR}/src/exacorr/libexacorr.a
    ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/lib/libtalsh.a
  )

target_include_directories(dirac
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

if(ENABLE_SRDFT)
  target_compile_definitions(dirac PRIVATE MOD_SRDFT)
endif()
