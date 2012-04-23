foreach(
    EXECUTABLE
    ${LIST_OF_EXECUTABLES}
    )
    install(
        TARGETS ${EXECUTABLE}
        DESTINATION share/dirac
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ             GROUP_EXECUTE
        WORLD_READ             WORLD_EXECUTE
        )
endforeach()

install(
    FILES ${PROJECT_BINARY_DIR}/pam-dirac
    DESTINATION bin
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )

if(ENABLE_PCMMOD)
    install( # Install PCMSolver input parsing script
        FILES ${PROJECT_SOURCE_DIR}/external/pcmsolver-build/src/solverlib/pcmsolver
        DESTINATION bin
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ             GROUP_EXECUTE
        WORLD_READ             WORLD_EXECUTE
        )
    install( # Install libraries 
        DIRECTORY ${PROJECT_SOURCE_DIR}/external/share/libgetkw
        DESTINATION share/dirac
        )
endif()

install(
    FILES ${PROJECT_SOURCE_DIR}/utils/wrapper.py
    DESTINATION share/dirac
    PERMISSIONS
    OWNER_READ OWNER_WRITE
    GROUP_READ
    WORLD_READ
    )

install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/basis ${PROJECT_SOURCE_DIR}/basis_dalton
    DESTINATION share/dirac
    PATTERN .git EXCLUDE
    )

include(InstallRequiredSystemLibraries)

set(CPACK_SOURCE_GENERATOR "TGZ")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "DIRAC")
set(CPACK_PACKAGE_VENDOR "R. Bast, H. J. Aa. Jensen, T. Saue, and L. Visscher, and contributors")
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
include(IncludeFilesList)

if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_custom_target(release
        COMMAND mkdir -p ${EXPORT_DIR}
        COMMAND mkdir -p ${EXPORT_DIR}/doc
        COMMAND tar cf ${EXPORT_DIR}/sources.tar -C ${CMAKE_SOURCE_DIR} ${DIRAC_SOURCES} ${INCLUDE_FILES}
        COMMAND tar xf ${EXPORT_DIR}/sources.tar -C ${EXPORT_DIR}
        COMMAND rm     ${EXPORT_DIR}/sources.tar
        COMMAND rm -rf ${EXPORT_DIR}/build
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/cmake                        ${EXPORT_DIR}
        COMMAND ${CMAKE_SOURCE_DIR}/maintenance/remove_unreleased_code  ${EXPORT_DIR}
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/basis                        ${EXPORT_DIR}
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/basis_dalton                 ${EXPORT_DIR}
        COMMAND ${CMAKE_SOURCE_DIR}/maintenance/remove_basis_sets       ${EXPORT_DIR}
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/test                         ${EXPORT_DIR}
        COMMAND ${CMAKE_SOURCE_DIR}/maintenance/remove_unreleased_tests ${EXPORT_DIR}
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/utils                        ${EXPORT_DIR}
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/xcfun                        ${EXPORT_DIR}
        COMMAND cp -rL ${CMAKE_SOURCE_DIR}/doc                          ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/pam                          ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/runtest                      ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/setup                        ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/CMakeLists.txt               ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/COPYING                      ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/CTestConfig.cmake            ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/INSTALL.rst                  ${EXPORT_DIR}
        COMMAND ${CMAKE_SOURCE_DIR}/maintenance/remove_copyright_tags   ${EXPORT_DIR}
        COMMAND make package_source
        COMMENT "Packaging source files"
        VERBATIM
        )
endif()
