# Hardware-free C++ unit tests for the bonded engine (Boost.Test, ctest).
# Mirrors the former in-UHD setup: each test is its own executable with a
# Boost-provided main.
add_definitions(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

list(APPEND bonded_tests
    overlap_reconstructor_test   # alignment recovery, band-flatness, seam, streaming
    bonded_receiver_test         # burst-alignment / ring-buffer extraction (header-only)
    bonded_transmitter_test)     # config / freq-plan geometry

foreach(test_name ${bonded_tests})
    add_executable(${test_name} ${test_name}.cpp)
    target_include_directories(${test_name} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../lib
        ${UHD_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS})
    target_link_libraries(${test_name}
        gnuradio-bonded_usrp
        Boost::unit_test_framework)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
