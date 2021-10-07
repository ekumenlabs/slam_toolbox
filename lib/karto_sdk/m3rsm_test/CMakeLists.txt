include_directories (
  ${PROJECT_SOURCE_DIR}/m3rsm_test/gtest/include
  ${PROJECT_SOURCE_DIR}/m3rsm_test/gtest
  ${PROJECT_SOURCE_DIR}/m3rsm_test
  ${PROJECT_BINARY_DIR}
)

# configure_file (test_config.h.in ${PROJECT_BINARY_DIR}/include/rndf_gazebo_plugin/test_config.h)

# Build gtest
add_library(gtest STATIC gtest/src/gtest-all.cc)
add_library(gtest_main STATIC gtest/src/gtest_main.cc)
target_link_libraries(gtest_main gtest)

set_target_properties(gtest PROPERTIES CXX_CPPLINT "")
set_target_properties(gtest_main PROPERTIES CXX_CPPLINT "")

set(GTEST_LIBRARY "${PROJECT_BINARY_DIR}/test/libgtest.a")
set(GTEST_MAIN_LIBRARY "${PROJECT_BINARY_DIR}/test/libgtest_main.a")

execute_process(COMMAND cmake -E remove_directory ${CMAKE_BINARY_DIR}/test_results)
execute_process(COMMAND cmake -E make_directory ${CMAKE_BINARY_DIR}/test_results)
include_directories(${GTEST_INCLUDE_DIRS})

# This should be migrated to more fine control solution based on set_property APPEND
# directories. It's present on cmake 2.8.8 while precise version is 2.8.7
link_directories(${PROJECT_BINARY_DIR}/test)
include_directories("${PROJECT_SOURCE_DIR}/test/gtest/include")

macro (m3rsm_tests_build)
  # Build all the tests
  foreach(GTEST_SOURCE_file ${ARGN})
    string(REGEX REPLACE ".cc" "" BINARY_NAME ${GTEST_SOURCE_file})
    message(${BINARY_NAME})
    set(BINARY_NAME ${TEST_TYPE}_${BINARY_NAME})
    if(USE_LOW_MEMORY_TESTS)
      add_definitions(-DUSE_LOW_MEMORY_TESTS=1)
    endif(USE_LOW_MEMORY_TESTS)

	include_directories(
		test/gtest/include
	)
    add_executable(${BINARY_NAME} ${GTEST_SOURCE_file})

    add_dependencies(${BINARY_NAME}
      gtest
    	gtest_main
    )

    target_link_libraries(${BINARY_NAME}
    	kartoSlamToolbox
    )

    target_link_libraries(${BINARY_NAME}
        libgtest.a
        libgtest_main.a
        pthread
    )

    add_test(${BINARY_NAME} ${CMAKE_CURRENT_BINARY_DIR}/${BINARY_NAME}
        --gtest_output=xml:${CMAKE_BINARY_DIR}/test_results/${BINARY_NAME}.xml)
    set_tests_properties(${BINARY_NAME} PROPERTIES TIMEOUT 240)
  endforeach()
endmacro()

add_subdirectory(src)
