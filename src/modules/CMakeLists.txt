# Fortran modules
set(MODULE_SOURCES
    primitives.F90
    shells.F90
    maxaqn.F90
)

# Create TARGET
add_library(modules OBJECT ${MODULE_SOURCES})

# Export directory that contains .mod files
target_include_directories(modules INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
