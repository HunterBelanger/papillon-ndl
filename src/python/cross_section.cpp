/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2022, Hunter Belanger
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

#include <PapillonNDL/cross_section.hpp>

namespace py = pybind11;

using namespace pndl;

void init_CrossSection(py::module& m) {
  py::class_<CrossSection, std::shared_ptr<CrossSection>>(m, "CrossSection")
      .def(py::init<const ACE&, size_t, std::shared_ptr<EnergyGrid>, bool,
                    bool>(),
           py::arg("ace"), py::arg("i"), py::arg("E_grid"),
           py::arg("get_index") = true, py::arg("is_heating") = false)
      .def(py::init<const std::vector<double>&, std::shared_ptr<EnergyGrid>,
                    size_t>())
      .def(py::init<double, std::shared_ptr<EnergyGrid>>())
      .def("__getitem__", &CrossSection::operator[])
      .def("__call__",
           py::overload_cast<double>(&CrossSection::operator(), py::const_))
      .def("__call__", py::overload_cast<double, size_t>(
                           &CrossSection::operator(), py::const_))
      .def("__call__", py::overload_cast<double, size_t, double, double>(
                           &CrossSection::operator(), py::const_))
      .def("evaluate",
           py::overload_cast<double>(&CrossSection::evaluate, py::const_))
      .def("evaluate", py::overload_cast<double, size_t>(
                           &CrossSection::evaluate, py::const_))
      .def("evaluate", py::overload_cast<double, size_t, double, double>(
                           &CrossSection::evaluate, py::const_))
      .def("size", &CrossSection::size)
      .def("index", &CrossSection::index)
      .def("xs", py::overload_cast<size_t>(&CrossSection::xs, py::const_))
      .def("xs", py::overload_cast<>(&CrossSection::xs, py::const_))
      .def("energy",
           py::overload_cast<size_t>(&CrossSection::energy, py::const_))
      .def("energy", py::overload_cast<>(&CrossSection::energy, py::const_))
      .def("energy_grid", &CrossSection::energy_grid);
}
