#
# Making relccsd as separate library
#
# (agreed in issue #143, Miro, Rado, Luuk)
#

option(ENABLE_RELCCSD_STANDALONE "Enable standalone RELCCSD" ON)

add_subdirectory(${PROJECT_SOURCE_DIR}/src/relccsd)
# miro: "relccsd" library must be before math containing "external_libs"
set(EXTERNAL_LIBS relccsd ${EXTERNAL_LIBS})
#
message(STATUS "Enable compilation of standalone relccsd.x: ${ENABLE_RELCCSD_STANDALONE}")
