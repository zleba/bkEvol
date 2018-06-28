#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Solver.h"

namespace py = pybind11;

PYBIND11_MODULE(bkevol, m) {

    py::class_<Settings>(m, "Settings")
        .def(py::init<>())
        .def("printInfo", &Settings::printInfo)
        .def_readwrite("alphaS", &Settings::asMZ)
        .def_readwrite("freezingScale", &Settings::freezingScale)
        .def_readwrite("eps", &Settings::eps)
        .def_readwrite("mu2", &Settings::mu2)

        .def_readwrite("NkT2", &Settings::N)
        .def_readwrite("NkT2int", &Settings::Nint)
        .def_readwrite("kT2Min", &Settings::kT2Min)
        .def_readwrite("kT2Max", &Settings::kT2Max)

        .def_readwrite("Nrap", &Settings::Nrap)
        .def_readwrite("xMin", &Settings::xMin)
        .def_readwrite("xMax", &Settings::xMax)

        .def_readwrite("kernelType", &Settings::kernelType)

        .def_readwrite("pars", &Settings::pars)
        .def_readwrite("funStr", &Settings::funStr);


    py::class_<Solver>(m, "Solver")
        .def(py::init<Settings>())

        .def("CalcEvolKernel", &Solver::CalcEvolKernel)
        .def("EvolveAll", &Solver::EvolveAll)
        .def("PrintGrid", &Solver::PrintGrid)
        .def("PrintBaseGrid", &Solver::PrintBaseGrid)
        .def("PrintReduce", &Solver::PrintReduce)
        .def("CalcF2L", &Solver::CalcF2L)
        .def("SaveEvolKernels", &Solver::SaveEvolKernels)
        .def("LoadEvolKernels", &Solver::LoadEvolKernels);


}
