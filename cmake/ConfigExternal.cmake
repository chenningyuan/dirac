add_custom_target(
    git_update
    COMMAND git submodule init
    COMMAND git submodule update
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

include(ExternalProject)

macro(add_external _project)
    ExternalProject_Add(${_project}
        DOWNLOAD_COMMAND git submodule update
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
        PREFIX ${PROJECT_SOURCE_DIR}/external
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/${_project}
        BINARY_DIR ${PROJECT_SOURCE_DIR}/external/${_project}-build
        INSTALL_DIR ${PROJECT_SOURCE_DIR}/external
        CMAKE_ARGS ${ExternalProjectCMakeArgs}
        )
    include_directories(${PROJECT_SOURCE_DIR}/external/${_project}-build)
    include_directories(${PROJECT_SOURCE_DIR}/external/${_project}-build/modules)
    link_directories(${PROJECT_SOURCE_DIR}/external/lib)
    link_directories(${CMAKE_SOURCE_DIR}/external/${_project}/external/lib)
    add_dependencies(${_project} git_update)

  # # remove stamps for external builds so that they are rebuilt every time
  # add_custom_command(
  #     TARGET ${_project}
  #     PRE_BUILD
  #     COMMAND rm -rf ${PROJECT_SOURCE_DIR}/external/src/${_project}-stamp
  #     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  #     )
endmacro()
