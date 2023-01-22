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

#include <PapillonNDL/pctable.hpp>

namespace py = pybind11;

using namespace pndl;

void init_PCTable(py::module& m) {
  py::class_<PCTable>(m, "PCTable")
      .def(py::init<const ACE&, size_t, double>())
      .def(py::init<const std::vector<double>&, const std::vector<double>&,
                    const std::vector<double>&, Interpolation>())
      .def("sample_value", &PCTable::sample_value)
      .def("min_value", &PCTable::min_value)
      .def("max_value", &PCTable::max_value)
      .def("interpolation", &PCTable::interpolation)
      .def("values", &PCTable::values)
      .def("pdf", py::overload_cast<double>(&PCTable::pdf, py::const_))
      .def("pdf", py::overload_cast<>(&PCTable::pdf, py::const_))
      .def("cdf", &PCTable::cdf)
      .def("size", &PCTable::size);
}
