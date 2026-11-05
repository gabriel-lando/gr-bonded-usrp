list(APPEND bonded_usrp_sources
    # GNU Radio block wrappers
    bonded_source_impl.cc
    overlap_reconstructor_impl.cc
    bonded_sink_impl.cc
    band_splitter_impl.cc
    # DSP and device engine (public UHD API only)
    overlap_reconstructor.cpp
    bonded_receiver.cpp
    bonded_transmitter.cpp
    band_splitter.cpp)

set(bonded_usrp_sources
    "${bonded_usrp_sources}"
    PARENT_SCOPE)

add_library(gnuradio-bonded_usrp SHARED ${bonded_usrp_sources})

target_link_libraries(gnuradio-bonded_usrp
    PUBLIC
    gnuradio::gnuradio-runtime
    ${UHD_LIBRARIES})

target_include_directories(gnuradio-bonded_usrp
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    PUBLIC $<INSTALL_INTERFACE:include>
    PRIVATE ${UHD_INCLUDE_DIRS}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_definitions(gnuradio-bonded_usrp
    PRIVATE gnuradio_bonded_usrp_EXPORTS)

# VOLK is optional: it accelerates the one-shot dot products and the NCO
# rotators; scalar fallbacks are used without it.
find_library(VOLK_LIB NAMES volk)
find_path(VOLK_INC NAMES volk/volk.h)
if(VOLK_LIB AND VOLK_INC)
    message(STATUS "Found VOLK: ${VOLK_LIB}")
    target_link_libraries(gnuradio-bonded_usrp PRIVATE ${VOLK_LIB})
    target_include_directories(gnuradio-bonded_usrp PRIVATE ${VOLK_INC})
    set_property(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/overlap_reconstructor.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/band_splitter.cpp
        APPEND PROPERTY COMPILE_DEFINITIONS BONDED_HAVE_VOLK=1)
endif()

# Aggressive vectorization, scoped to the reconstructor hot loops only.
# -march=native: the build is tuned for the machine it is compiled on.
if(NOT MSVC)
    set_property(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/overlap_reconstructor.cpp
        APPEND PROPERTY COMPILE_OPTIONS
        -O3 -ffast-math -funroll-loops -march=native)
endif()

set_target_properties(gnuradio-bonded_usrp PROPERTIES
    OUTPUT_NAME gnuradio-bonded_usrp)

install(TARGETS gnuradio-bonded_usrp
    EXPORT gnuradio-bonded_usrpTargets
    DESTINATION ${GR_LIBRARY_DIR})

install(EXPORT gnuradio-bonded_usrpTargets
    FILE gnuradio-bonded_usrpTargets.cmake
    NAMESPACE gnuradio::
    DESTINATION ${GR_CMAKE_DIR})
