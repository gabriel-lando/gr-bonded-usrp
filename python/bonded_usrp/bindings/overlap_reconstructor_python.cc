#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <gnuradio/bonded_usrp/overlap_reconstructor.h>

namespace py = pybind11;

void bind_overlap_reconstructor(py::module& m)
{
    using overlap_reconstructor = ::gr::bonded_usrp::overlap_reconstructor;

    py::class_<overlap_reconstructor,
               gr::block,
               gr::basic_block,
               std::shared_ptr<overlap_reconstructor>>(m, "overlap_reconstructor")
        .def(py::init(&overlap_reconstructor::make),
             py::arg("rate_in"),
             py::arg("rate_out"),
             py::arg("target_bw"),
             py::arg("overlap_bw"),
             py::arg("radio1_offset"),
             py::arg("radio2_offset"),
             py::arg("fft_size") = 4096);
}
