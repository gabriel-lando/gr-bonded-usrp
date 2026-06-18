#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_bonded_source(py::module& m);
void bind_overlap_reconstructor(py::module& m);
void bind_bonded_sink(py::module& m);

PYBIND11_MODULE(bonded_usrp_python, m)
{
    py::module::import("gnuradio.gr");
    bind_bonded_source(m);
    bind_overlap_reconstructor(m);
    bind_bonded_sink(m);
}
