#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "particles.h"


namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(cParticle, m) {

    py::class_<Box, std::shared_ptr<Box>>(m, "Box")
        .def(py::init<double, double, int>())
        .def_readwrite("width", &Box::width)
        .def_readwrite("height", &Box::height)
        .def_readwrite("dt", &Box::dt)
        .def_readwrite("particleSteps", &Box::particleSteps)
        .def_readwrite("pressureSteps", &Box::pressureSteps)
        // .def_readwrite("collisionSteps", &Box::collisionSteps)
        .def("addRandomParticle", &Box::addRandomParticle)
        .def("runSim", &Box::runSim);
}

