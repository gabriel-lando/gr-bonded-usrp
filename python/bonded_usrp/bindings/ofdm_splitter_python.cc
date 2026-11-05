#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <gnuradio/bonded_usrp/ofdm_splitter.h>

namespace py = pybind11;

void bind_ofdm_splitter(py::module& m)
{
    using ofdm_splitter = ::gr::bonded_usrp::ofdm_splitter;

    py::class_<ofdm_splitter,
               gr::block,
               gr::basic_block,
               std::shared_ptr<ofdm_splitter>>(m, "ofdm_splitter")
        .def(py::init(&ofdm_splitter::make),
             py::arg("rate_in"),
             py::arg("rate_out"),
             py::arg("radio1_offset"),
             py::arg("radio2_offset"),
             py::arg("fft_len")     = 64,
             py::arg("cp_len")      = 16,
             py::arg("window_len")  = 8,
             py::arg("len_tag_key") = "packet_len");
}
