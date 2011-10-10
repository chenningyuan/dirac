include(InstallRequiredSystemLibraries)

set(CPACK_SOURCE_GENERATOR "TGZ")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "DIRAC")
set(CPACK_PACKAGE_VENDOR "T. Saue, L. Visscher, H. J. Aa. Jensen, and contributors")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "R. Bast")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/INSTALL.rst")
set(CPACK_RESOURCE_FILE_LICENSE    "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR    "${VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR    "${VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH    "${VERSION_PATCH}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY
    "DIRAC ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}"
    )

SET(EXPORT_DIR ${CMAKE_BINARY_DIR}/export)
SET(CPACK_SOURCE_INSTALLED_DIRECTORIES
    "${EXPORT_DIR};/"
    )

include(CPack)

file(GLOB_RECURSE H_FILES   RELATIVE ${CMAKE_SOURCE_DIR} "*.h")
file(GLOB_RECURSE INC_FILES RELATIVE ${CMAKE_SOURCE_DIR} "*.inc")

set(INCLUDE_FILES
    ${H_FILES}
    ${INC_FILES}
    )

if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_custom_target(release
        COMMAND mkdir -p ${EXPORT_DIR}
        COMMAND mkdir -p ${EXPORT_DIR}/doc
        COMMAND tar cf ${EXPORT_DIR}/sources.tar -C ${CMAKE_SOURCE_DIR} ${DIRAC_SOURCES} ${INCLUDE_FILES}
        COMMAND tar xf ${EXPORT_DIR}/sources.tar -C ${EXPORT_DIR}
        COMMAND rm     ${EXPORT_DIR}/sources.tar
        COMMAND rm -rf ${EXPORT_DIR}/build
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/cmake              ${EXPORT_DIR}
        COMMAND ${CMAKE_SOURCE_DIR}/maintenance/remove_unreleased ${EXPORT_DIR}
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/basis             ${EXPORT_DIR}
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/basis_dalton      ${EXPORT_DIR}
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/test              ${EXPORT_DIR}
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/utils             ${EXPORT_DIR}
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/xcint/xcfun       ${EXPORT_DIR}/xcint
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/cbuild            ${EXPORT_DIR}
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/modules           ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/pam               ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/runtest           ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/setup             ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/CMakeLists.txt    ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/COPYING           ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/CTestConfig.cmake ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/INSTALL.rst       ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/config.guess      ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/config.status     ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/config.sub        ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/configure         ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/Makefile          ${EXPORT_DIR}
        COMMAND make package_source
        COMMENT "Packaging source files"
        VERBATIM
        )
endif()
