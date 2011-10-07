if (ENABLE_MPI) 
if (ENABLE_SLOWDOWN_SLAVE)
# mi: attempt to impose "-DSLOWDOWN_SLAVES" command through COMPILE_FLAGS, but it's not working with cmake 2.8.0...
# see http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:set_target_properties
# "The COMPILE_FLAGS property sets additional compiler flags used to build sources within the target. 
# It may also be used to pass additional preprocessor definitions."
    set_target_properties(dirac.x PROPERTIES COMPILE_FLAGS "-DSLOWDOWN_SLAVES")
# I use rather this substitution 
    add_definitions(-DSLOWDOWN_SLAVES)
endif()
# touch this file in both Linux & Windows in connection with "-DSLOWDOWN_SLAVES" directive
execute_process( 
    #COMMAND python -c "import os; os.utime('${CMAKE_SOURCE_DIR}/main/dirac.F', None)"
    COMMAND python -c "import os; os.utime('main/dirac.F', None)"
    TIMEOUT 1
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_QUIET
)
endif()
