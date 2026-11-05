# Hardware-free C++ unit tests (Boost.Test, ctest); each test is its own
# executable with a Boost-provided main.
add_definitions(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

list(APPEND bonded_tests
    overlap_reconstructor_test   # alignment, flatness, seam, frozen path, threads, direct output, drift
    bonded_resample_test         # filter design, dot kernels, polyphase resamplers, nco
    bonded_receiver_test         # timestamp alignment of the receive rings
    bonded_transmitter_test      # config checks, phase API, frequency/phase correction
    band_splitter_test           # TX split: flatness, half assignment, seam width, streaming, drift
    ofdm_splitter_test           # OFDM TX split: any relative phase, subcarrier assignment, timeline
    bonded_csv_test)             # block CSV parameters

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
