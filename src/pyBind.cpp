#include <pybind11/pybind11.h>
#include "Solver.h"

namespace py = pybind11;

PYBIND11_MODULE(bkevol, m) {

    py::class_<Settings>(m, "Settings")
        .def(py::init<>())
        .def_readwrite("asMZ", &Settings::asMZ)
        .def_readwrite("eps", &Settings::eps)
        .def_readwrite("Nint", &Settings::Nint)
        .def_readwrite("N", &Settings::N)
        .def_readwrite("Nrap", &Settings::Nrap)
        .def_readwrite("inputDir", &Settings::inputDir)
        .def_readwrite("funStr", &Settings::funStr);


    py::class_<Solver>(m, "Solver")
        .def(py::init<int>())
        .def(py::init<Settings>())

        .def("InitMat", &Solver::InitMat)
        .def("EvolveNew", &Solver::EvolveNew)
        .def("PrintGrid", &Solver::PrintGrid)
        .def("PrintBaseGrid", &Solver::PrintBaseGrid)
        .def("PrintReduce", &Solver::PrintReduce)
        .def("CalcF2L", &Solver::CalcF2L);




}
