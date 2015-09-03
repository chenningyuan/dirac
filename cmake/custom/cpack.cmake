include(InstallRequiredSystemLibraries)

set(CPACK_SOURCE_GENERATOR "TGZ")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "DIRAC")
set(CPACK_PACKAGE_VENDOR "R. Bast, H. J. Aa. Jensen, T. Saue, and L. Visscher, and contributors")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "R. Bast")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/doc/installation/general.rst")
set(CPACK_RESOURCE_FILE_LICENSE    "${PROJECT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_FILE_NAME        "DIRAC")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "DIRAC-${PROGRAM_VERSION}-Source")
set(CPACK_PACKAGE_INSTALL_DIRECTORY
    "DIRAC ${PROGRAM_VERSION}"
    )

SET(EXPORT_DIR ${PROJECT_BINARY_DIR}/export)
SET(CPACK_SOURCE_INSTALLED_DIRECTORIES
    "${EXPORT_DIR};/"
    )

include(CPack)
include(IncludeFiles)

if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_custom_target(release
        COMMAND ${PROJECT_SOURCE_DIR}/maintenance/fetch_external_sources ${PROJECT_SOURCE_DIR}
        COMMAND mkdir -p ${EXPORT_DIR}
        COMMAND tar cf ${EXPORT_DIR}/sources.tar -C ${PROJECT_SOURCE_DIR} ${DIRAC_SOURCES} ${INCLUDE_FILES}
        COMMAND tar xf ${EXPORT_DIR}/sources.tar -C ${EXPORT_DIR}
        COMMAND rm     ${EXPORT_DIR}/sources.tar
        COMMAND rm -rf ${EXPORT_DIR}/build
        COMMAND ${PROJECT_SOURCE_DIR}/maintenance/remove_unreleased_code  ${EXPORT_DIR} ${PROJECT_BINARY_DIR}
        COMMAND cp -r  ${PROJECT_SOURCE_DIR}/basis                        ${EXPORT_DIR}
        COMMAND cp -r  ${PROJECT_SOURCE_DIR}/basis_dalton                 ${EXPORT_DIR}
        COMMAND cp -r  ${PROJECT_SOURCE_DIR}/basis_ecp                    ${EXPORT_DIR}
        COMMAND ${PROJECT_SOURCE_DIR}/maintenance/remove_basis_sets       ${EXPORT_DIR}
        COMMAND cp -r  ${PROJECT_SOURCE_DIR}/test                         ${EXPORT_DIR}
        COMMAND ${PROJECT_SOURCE_DIR}/maintenance/remove_unreleased_tests ${EXPORT_DIR} ${PROJECT_BINARY_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/CMakeLists.txt               ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/COPYING                      ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/CHANGELOG.rst                ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/CTestConfig.cmake            ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/atomic-start-x2c-pam         ${EXPORT_DIR}
        COMMAND cp -r  ${PROJECT_SOURCE_DIR}/cmake                        ${EXPORT_DIR}
        COMMAND cp -r  ${PROJECT_SOURCE_DIR}/external                     ${EXPORT_DIR}
        COMMAND rm -f  ${EXPORT_DIR}/external/pcmsolver/.git  # prevens "git init; git add ." inside tarball
        COMMAND rm -f  ${EXPORT_DIR}/external/unit-tests/.git # prevens "git init; git add ." inside tarball
        COMMAND cp -r  ${PROJECT_BINARY_DIR}/pamadm-generated           ${EXPORT_DIR}/cmake
        COMMAND cp -rL ${PROJECT_SOURCE_DIR}/doc                          ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/pam.in                       ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/setup.py                     ${EXPORT_DIR}
        COMMAND cp -r  ${PROJECT_SOURCE_DIR}/utils                        ${EXPORT_DIR}
        COMMAND cp -r  ${PROJECT_SOURCE_DIR}/src/xcfun                    ${EXPORT_DIR}/src
        COMMAND cp     ${PROJECT_SOURCE_DIR}/src/main/main.F90            ${EXPORT_DIR}/src/main
        COMMAND cp     ${PROJECT_SOURCE_DIR}/VERSION                      ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/.gitignore                   ${EXPORT_DIR}
        COMMAND ${PROJECT_SOURCE_DIR}/maintenance/remove_copyright_tags   ${EXPORT_DIR}
        COMMAND echo "${GIT_COMMIT_HASH}"                                > ${EXPORT_DIR}/cmake/GIT_HASH
        COMMAND make package_source
        COMMENT "Packaging source files"
        VERBATIM
        )
endif()
