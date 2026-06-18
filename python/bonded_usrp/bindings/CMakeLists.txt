include(GrPybind)

list(APPEND bonded_usrp_python_files
    bonded_source_python.cc
    overlap_reconstructor_python.cc
    bonded_sink_python.cc
    python_bindings.cc)

pybind11_add_module(bonded_usrp_python ${bonded_usrp_python_files})

target_include_directories(bonded_usrp_python PRIVATE
    ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(bonded_usrp_python PRIVATE gnuradio-bonded_usrp)

install(TARGETS bonded_usrp_python
    DESTINATION ${GR_PYTHON_DIR}/gnuradio/bonded_usrp
    COMPONENT pythonapi)
