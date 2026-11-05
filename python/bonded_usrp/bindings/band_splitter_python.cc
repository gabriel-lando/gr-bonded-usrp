#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <gnuradio/bonded_usrp/band_splitter.h>

namespace py = pybind11;

void bind_band_splitter(py::module& m)
{
    using band_splitter = ::gr::bonded_usrp::band_splitter;

    py::class_<band_splitter,
               gr::block,
               gr::basic_block,
               std::shared_ptr<band_splitter>>(m, "band_splitter")
        .def(py::init(&band_splitter::make),
             py::arg("rate_in"),
             py::arg("rate_out"),
             py::arg("radio1_offset"),
             py::arg("radio2_offset"),
             py::arg("taps_per_phase") = 64);
}
