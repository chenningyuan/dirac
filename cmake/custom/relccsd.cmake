#
# Making relccsd as separate library
#
# (agreed in issue #143, Miro, Rado, Luuk)
#

option(ENABLE_RELCCSD_STANDALONE "Enable standalone RELCCSD" ON)

add_subdirectory(${PROJECT_SOURCE_DIR}/src/relccsd)
set(EXTERNAL_LIBS ${EXTERNAL_LIBS} relccsd)
#
message(STATUS "Enable compilation of standalone relccsd.x: ${ENABLE_RELCCSD_STANDALONE}")
