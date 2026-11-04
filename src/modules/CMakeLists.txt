# Fortran modules
set(MODULE_SOURCES
    precision.F90
    aovec.F90
    primit.F90
    #shells.F90
    maxaqn.F90
    maxorb.F90
    mxcent.F90
    symmet.F90
    #fitshells.F90
    #fitprimit.F90
    #nqcc.F90 # From here, importing only used variables
    #pincom.F90
    #parint.F90
    #dorps.F90
    #maxmom.F90
    #xyzpow.F90
    #rotg.F90
    #infdim.F90
    #nucdata.F90
    #symind.F90
    #dcbesr.F90
    #dcbnmr.F90
    #ccom.F90
    #dipole.F90
    #consts.F90
    priunit.F90
    #twocom.F90
    #twosta.F90
    #dirprt.F90
    #doxyz.F90
)

# Create TARGET
add_library(modules OBJECT ${MODULE_SOURCES})

# Export directory that contains .mod files
target_include_directories(modules INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
