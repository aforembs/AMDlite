cmake_minimum_required(VERSION 3.10)
project(amd_lite VERSION 1.0.0 LANGUAGES C)

option(AMD_LITE_BUILD_TESTS "Build amd_lite tests" ON)
option(AMD_LITE_BUILD_DEMO  "Build amd_lite demo"  ON)

add_library(amd_lite src/amd_lite.c)
target_include_directories(amd_lite PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(amd_lite PUBLIC c_std_99)
if (NOT MSVC)
    target_link_libraries(amd_lite PRIVATE m)
    target_compile_options(amd_lite PRIVATE -Wall -Wextra -Wpedantic)
endif()

if (AMD_LITE_BUILD_TESTS)
    enable_testing()
    add_executable(test_amd_lite tests/test_amd_lite.c)
    target_link_libraries(test_amd_lite PRIVATE amd_lite)
    add_test(NAME amd_lite_tests COMMAND test_amd_lite)
endif()

if (AMD_LITE_BUILD_DEMO)
    add_executable(amd_lite_demo demo/amd_lite_demo.c)
    target_link_libraries(amd_lite_demo PRIVATE amd_lite)
endif()

include(GNUInstallDirs)
install(TARGETS amd_lite
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES include/amd_lite.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
