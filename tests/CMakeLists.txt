file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/generated)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/generated/testSuites.inc "")

set(PFUNIT_DIR $ENV{PFUNIT})

include_directories(
    ${PROJECT_BINARY_DIR}/src
    ${CMAKE_CURRENT_BINARY_DIR}/generated
    ${PFUNIT_DIR}/mod
    )

file(GLOB test_files test_*.pf)

foreach(test ${test_files})

    get_filename_component(test_we ${test} NAME_WE)    

    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/generated/${test_we}.F90
        COMMAND ${PFUNIT_DIR}/bin/pFUnitParser.py ${test} ${CMAKE_CURRENT_BINARY_DIR}/generated/${test_we}.F90
    )

    list(APPEND test_sources ${CMAKE_CURRENT_BINARY_DIR}/generated/${test_we}.F90)
    message("test_source: ${test_sources}")
    
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/generated/testSuites.inc "ADD_TEST_SUITE(${test_we}_suite)\n")

endforeach()

add_executable(
    pftest_alltests
    ${PFUNIT_DIR}/include/driver.F90
    ${test_sources}
    )

target_link_libraries(
    pftest_alltests
    libtreecode
    ${PFUNIT_DIR}/lib/libpfunit.a
    )

file(COPY test_files DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()
add_test(pftest_alltests ${PROJECT_BINARY_DIR}/pftest_alltests)
