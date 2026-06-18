include(GrPython)
if(NOT PYTHONINTERP_FOUND)
    return()
endif()

add_subdirectory(bindings)

gr_python_install(FILES __init__.py
    DESTINATION ${GR_PYTHON_DIR}/gnuradio/bonded_usrp)
