#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <gnuradio/bonded_usrp/bonded_source.h>

namespace py = pybind11;

void bind_bonded_source(py::module& m)
{
    using bonded_source = ::gr::bonded_usrp::bonded_source;

    py::class_<bonded_source,
               gr::sync_block,
               gr::block,
               gr::basic_block,
               std::shared_ptr<bonded_source>>(m, "bonded_source")
        .def(py::init(&bonded_source::make),
             py::arg("serials_csv"),
             py::arg("rate"),
             py::arg("freq"),
             py::arg("gain"),
             py::arg("clock_source") = "external",
             py::arg("time_source") = "external",
             py::arg("strict_lock") = false,
             py::arg("lock_timeout") = 5.0,
             py::arg("fetch_timeout") = 1.0,
             py::arg("subdev") = "",
             py::arg("freq_plan_csv") = "",
             py::arg("stream_args") = "");
}
