if (ENABLE_MPI) 
if (ENABLE_SLOWDOWN_SLAVE)
    set_target_properties(dirac.x PROPERTIES COMPILE_FLAGS "-DSLOWDOWN_SLAVES")
    add_definitions(-DSLOWDOWN_SLAVES)
endif()
# touch this file in both Linux & Windows in connection with "-DSLOWDOWN_SLAVES" directive
execute_process( 
    COMMAND python -c "import os; os.utime('main/dirac.F', None)"
    TIMEOUT 1
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_QUIET
)
endif()
