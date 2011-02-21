# this file is useful when calibrating cdash runs
# do not use it in production cdash
add_test(acmoin ${CMAKE_SOURCE_DIR}/runtest ${RUNTEST_FLAGS} --tests=acmoin)
add_test(basis_automatic_augmentation ${CMAKE_SOURCE_DIR}/runtest ${RUNTEST_FLAGS} --tests=basis_automatic_augmentation)
add_test(basis_contraction ${CMAKE_SOURCE_DIR}/runtest ${RUNTEST_FLAGS} --tests=basis_contraction)
add_test(basis_input ${CMAKE_SOURCE_DIR}/runtest ${RUNTEST_FLAGS} --tests=basis_input)
