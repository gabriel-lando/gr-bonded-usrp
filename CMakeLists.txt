cmake_minimum_required(VERSION 3.16)
project(gr-bonded_usrp CXX)

enable_testing()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake/Modules)

find_package(Gnuradio "3.10" REQUIRED COMPONENTS runtime)
find_package(UHD REQUIRED)

# Boost.Test for the hardware-free C++ unit tests (optional).
option(ENABLE_TESTING "Build the C++ unit tests" ON)
find_package(Boost COMPONENTS unit_test_framework)

include(GrVersion)
include(GrPlatform)
include(GrMinReq)
include(GrCompilerSettings)
include(GrPython)
include(GrPybind)

set(GR_PKG_DATA_DIR ${GR_DATA_DIR}/${CMAKE_PROJECT_NAME})

if(NOT CMAKE_MODULES_DIR)
    set(CMAKE_MODULES_DIR lib${LIB_SUFFIX}/cmake)
endif()

set(GR_INCLUDE_DIR include/gnuradio/bonded_usrp)
set(GR_CMAKE_DIR ${CMAKE_MODULES_DIR}/gnuradio-bonded_usrp)

add_subdirectory(include/gnuradio/bonded_usrp)
add_subdirectory(lib)
if(ENABLE_PYTHON)
    add_subdirectory(python/bonded_usrp)
endif()
add_subdirectory(grc)

if(ENABLE_TESTING AND Boost_FOUND)
    add_subdirectory(tests)
elseif(ENABLE_TESTING)
    message(STATUS "Boost.Test not found - skipping C++ unit tests")
endif()

install(FILES cmake/Modules/gnuradio-bonded_usrpConfig.cmake
    DESTINATION ${GR_CMAKE_DIR})
