include(CMakeFindDependencyMacro)
find_dependency(Gnuradio COMPONENTS runtime)
find_dependency(UHD)

if(NOT TARGET gnuradio-bonded_usrp)
    include("${CMAKE_CURRENT_LIST_DIR}/gnuradio-bonded_usrpTargets.cmake")
endif()
