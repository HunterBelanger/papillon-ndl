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
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/angle_distribution.hpp>

namespace py = pybind11;

using namespace pndl;

void init_AngleDistribution(py::module& m) {
  py::class_<AngleDistribution, std::shared_ptr<AngleDistribution>>(
      m, "AngleDistribution")
      .def(py::init<>())
      .def(py::init<const ACE&, int>())
      .def(py::init<const std::vector<double>&,
                    const std::vector<std::shared_ptr<AngleLaw>>&>())
      .def("sample_angle", &AngleDistribution::sample_angle)
      .def("pdf", &AngleDistribution::pdf)
      .def("size", &AngleDistribution::size)
      .def("energy",
           py::overload_cast<>(&AngleDistribution::energy, py::const_))
      .def("energy",
           py::overload_cast<size_t>(&AngleDistribution::energy, py::const_))
      .def("law", &AngleDistribution::law,
           py::return_value_policy::reference_internal);
}
