include(GrPython)
if(NOT PYTHONINTERP_FOUND)
    return()
endif()

add_subdirectory(bindings)

gr_python_install(FILES __init__.py
    DESTINATION ${GR_PYTHON_DIR}/gnuradio/bonded_usrp)

# Flowgraph tests of the blocks, run against the module in the build tree.
if(ENABLE_TESTING)
    foreach(qa qa_band_splitter qa_overlap_reconstructor)
        add_test(NAME ${qa}
            COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${qa}.py)
        set_tests_properties(${qa} PROPERTIES ENVIRONMENT
            "BONDED_USRP_PYTHON_DIR=${CMAKE_CURRENT_BINARY_DIR}/bindings")
    endforeach()
endif()
