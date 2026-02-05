#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "robot_dynamics.h"

namespace py = pybind11;

PYBIND11_MODULE(robot_dynamics_module, m) {
    py::class_<RobotDynamics>(m, "RobotDynamics")
        .def(py::init<>())
        .def("pmin_calc", &RobotDynamics::Pmin_calc, 
             "Calculate Pmin matrix",
             py::arg("q"), py::arg("dq"), py::arg("ddq"))
        .def("Yc_calc", &RobotDynamics::Yc_calc, 
             "Calculate Yc matrix",
             py::arg("dq"));
}