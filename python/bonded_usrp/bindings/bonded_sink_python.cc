#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <gnuradio/bonded_usrp/bonded_sink.h>

namespace py = pybind11;

void bind_bonded_sink(py::module& m)
{
    using bonded_sink = ::gr::bonded_usrp::bonded_sink;

    py::class_<bonded_sink,
               gr::sync_block,
               gr::block,
               gr::basic_block,
               std::shared_ptr<bonded_sink>>(m, "bonded_sink")
        .def(py::init(&bonded_sink::make),
             py::arg("serials_csv"),
             py::arg("rate"),
             py::arg("freq"),
             py::arg("gain"),
             py::arg("clock_source")   = "external",
             py::arg("time_source")    = "external",
             py::arg("strict_lock")    = false,
             py::arg("lock_timeout")   = 5.0,
             py::arg("subdev")         = "",
             py::arg("freq_plan_csv")  = "",
             py::arg("stream_args")    = "",
             py::arg("delay_trim_a")   = 0,
             py::arg("delay_trim_b")   = 0,
             py::arg("gain_plan_csv")  = "",
             py::arg("freq_offset_csv")  = "",
             py::arg("phase_offset_csv") = "");
}
