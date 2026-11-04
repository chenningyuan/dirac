# Fortran modules
set(MODULE_SOURCES
    aovec.F90
    cbione.F90
    chrnos.F90
    chrxyz.F90
    consts.F90
    cotabs.F90
    crossd.F90
    crsdir.F90
    csym1.F90
    derzer.F90
    dirprt.F90
    disbuf.F90
    doxyz.F90
    dummy.F90
    efield.F90
    expcom.F90
    gencon.F90
    hertop.F90
    ibtpar.F90
    iratdef.F90
    maxaqn.F90
    maxmom.F90
    maxorb.F90
    mxcent.F90
    mxsymm.F90
    nodint.F90
    nucdata.F90
    oneadr.F90
    pi.F90
    precision.F90
    primit.F90
    priunit.F90
    shells.F90
    sotabs.F90
    subdir.F90
    symmet.F90
    twoao.F90
    twocom.F90
    twosta.F90
    # drafts, not migrated yet
    #ccom.F90
    #dcbesr.F90
    #dcbnmr.F90
    #dipole.F90
    #dorps.F90
    #fitprimit.F90
    #fitshells.F90
    #infdim.F90
    #nqcc.F90
    #parint.F90
    #pincom.F90
    #rotg.F90
    #symind.F90
    #xyzpow.F90
)

# Create TARGET
add_library(modules OBJECT ${MODULE_SOURCES})

# Export directory that contains .mod files
target_include_directories(modules INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
