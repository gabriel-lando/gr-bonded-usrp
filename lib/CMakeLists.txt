list(APPEND bonded_usrp_sources
    # GNU Radio block wrappers
    bonded_source_impl.cc
    overlap_reconstructor_impl.cc
    bonded_sink_impl.cc
    # Bonded DSP / acquisition engine (formerly compiled into libuhd;
    # now self-contained in this module, links only public UHD).
    overlap_reconstructor.cpp
    bonded_receiver.cpp
    bonded_transmitter.cpp)

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

# ---------------------------------------------------------------------------
# Optional accelerators for the overlap_reconstructor DSP hot path.
# Ported from the former in-libuhd build. All optional: the reconstructor
# falls back to a built-in FFT / scalar loops when absent. Macro NAMES
# (UHD_HAVE_FFTW3F / UHD_HAVE_VOLK) are kept because the engine sources
# #ifdef on them.
# ---------------------------------------------------------------------------
find_library(FFTW3F_LIB NAMES fftw3f libfftw3f)
find_path(FFTW3F_INC NAMES fftw3.h)
if(FFTW3F_LIB AND FFTW3F_INC)
    message(STATUS "Found FFTW3f: ${FFTW3F_LIB} (overlap reconstructor accelerated)")
    target_link_libraries(gnuradio-bonded_usrp PRIVATE ${FFTW3F_LIB})
    target_include_directories(gnuradio-bonded_usrp PRIVATE ${FFTW3F_INC})
    set_property(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/overlap_reconstructor.cpp
        APPEND PROPERTY COMPILE_DEFINITIONS UHD_HAVE_FFTW3F=1)
endif()

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP (overlap reconstructor frame-parallel)")
    target_link_libraries(gnuradio-bonded_usrp PRIVATE OpenMP::OpenMP_CXX)
endif()

find_library(VOLK_LIB NAMES volk)
find_path(VOLK_INC NAMES volk/volk.h)
if(VOLK_LIB AND VOLK_INC)
    message(STATUS "Found VOLK: ${VOLK_LIB} (overlap reconstructor time-domain SIMD)")
    target_link_libraries(gnuradio-bonded_usrp PRIVATE ${VOLK_LIB})
    target_include_directories(gnuradio-bonded_usrp PRIVATE ${VOLK_INC})
    set_property(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/overlap_reconstructor.cpp
        APPEND PROPERTY COMPILE_DEFINITIONS UHD_HAVE_VOLK=1)
endif()

# Aggressive vectorization, scoped to the reconstructor hot loops only.
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
