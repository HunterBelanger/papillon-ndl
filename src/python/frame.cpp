/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2023, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Papillon Nuclear Data Library (PapillonNDL).
 *
 * PapillonNDL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PapillonNDL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PapillonNDL. If not, see <https://www.gnu.org/licenses/>.
 *
 * */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/frame.hpp>
#include <tuple>

namespace py = pybind11;

using namespace pndl;

void init_Frame(py::module& m) {
  py::enum_<Frame>(m, "Frame").value("Lab", Frame::Lab).value("CM", Frame::CM);

  // Because fundamental types (i.e. float) are imutable in Python, we need
  // to bind a lambda function which returns a tuple for the first transform
  // overload.
  py::class_<CMToLab>(m, "CMToLab")
      .def("transform",
           [](double Ein, double A, double mu, double Eout) {
             CMToLab::transform(Ein, A, mu, Eout);
             return std::make_tuple(mu, Eout);
           })
      .def("transform", py::overload_cast<double, double, AngleEnergyPacket&>(
                            &CMToLab::transform))
      .def("angle_jacobian", py::overload_cast<double, double, double, double>(
                                 &CMToLab::angle_jacobian))
      .def("angle_jacobian",
           py::overload_cast<double, double, double, double, double>(
               &CMToLab::angle_jacobian))
      .def("angle_jacobian",
           py::overload_cast<double, double, AngleEnergyPacket>(
               &CMToLab::angle_jacobian))
      .def("jacobian", &CMToLab::jacobian);

  py::class_<LabToCM>(m, "LabToCM")
      .def("transform",
           [](double Ein, double A, double mu, double Eout) {
             LabToCM::transform(Ein, A, mu, Eout);
             return std::make_tuple(mu, Eout);
           })
      .def("transform", py::overload_cast<double, double, AngleEnergyPacket&>(
                            &LabToCM::transform))
      .def("angle", &LabToCM::angle);
}
